"""Background task manager: queue, progress, cancel, restart recovery."""
from __future__ import annotations

import queue
import threading

from .designer import CancelledError as DesignCancelled
from .designer import design_primers
from .models import DesignParams, TaskRecord, TaskStatus
from .specificity import CancelledError as EvalCancelled
from .specificity import ReferenceLibrary, SpecificityEvaluator
from .storage import Storage


class TaskManager:
    """Single-consumer (default) task queue with durable status transitions.

    Persistence policy: status transitions are written to storage, but
    intermediate progress/results of a QUEUED or RUNNING task are not. On
    restart, unfinished tasks are re-queued and start from progress 0;
    only COMPLETED tasks carry a result, so no partial output can ever be
    served as a final answer.
    """

    def __init__(
        self,
        library: ReferenceLibrary,
        storage: Storage | None = None,
        max_workers: int = 1,
    ) -> None:
        self.library = library
        self.storage = storage
        self.max_workers = max(1, max_workers)
        self._tasks: dict[str, TaskRecord] = {}
        self._cancel_events: dict[str, threading.Event] = {}
        self._queue: queue.Queue[str] = queue.Queue()
        self._lock = threading.RLock()
        self._workers: list[threading.Thread] = []
        self._stop = threading.Event()

    # -- lifecycle ------------------------------------------------------

    def start(self) -> None:
        if self._workers:
            return
        self._stop.clear()
        # NOTE: do not clear the queue here -- tasks may already have been
        # submitted (or loaded from disk) before workers existed.
        for i in range(self.max_workers):
            t = threading.Thread(target=self._worker_loop, name=f"primer-worker-{i}",
                                 daemon=True)
            t.start()
            self._workers.append(t)

    def shutdown(self, wait: bool = True) -> None:
        self._stop.set()
        # drop queued-but-unstarted work; one sentinel per worker wakes them
        drained = 0
        while True:
            try:
                item = self._queue.get_nowait()
            except queue.Empty:
                break
            if item is not None:
                drained += 1
            self._queue.task_done()
        for _ in self._workers:
            self._queue.put(None)
        for t in self._workers:
            if wait:
                t.join(timeout=5.0)
        self._workers = []
        return drained

    def recover(self) -> int:
        """Re-queue persisted QUEUED/RUNNING tasks. Returns re-queued count.

        Task records are loaded into memory regardless of status (completed
        and cancelled records become queryable immediately); only unfinished
        ones are re-queued.
        """
        if self.storage is None:
            return 0
        self.storage.cleanup_temp()
        persisted = self.storage.load_tasks()
        requeued = 0
        with self._lock:
            for tid, task in persisted.items():
                self._tasks[tid] = task
                if task.status in (TaskStatus.QUEUED, TaskStatus.RUNNING):
                    task.status = TaskStatus.QUEUED
                    task.progress = 0.0
                    task.result = None
                    self._cancel_events[tid] = threading.Event()
                    self._queue.put(tid)
                    requeued += 1
            self._persist_unlocked()
        return requeued

    def load_persisted(self) -> int:
        """Load persisted task records into memory without re-queuing."""
        if self.storage is None:
            return 0
        with self._lock:
            for tid, task in self.storage.load_tasks().items():
                self._tasks.setdefault(tid, task)
            return len(self._tasks)

    # -- submission / queries ------------------------------------------

    def submit(self, params: DesignParams) -> TaskRecord:
        task = TaskRecord.new(params)
        with self._lock:
            self._tasks[task.task_id] = task
            self._cancel_events[task.task_id] = threading.Event()
            self._persist_unlocked()
        self._queue.put(task.task_id)
        return task

    def get(self, task_id: str) -> TaskRecord | None:
        with self._lock:
            return self._tasks.get(task_id)

    def list_tasks(self) -> list[TaskRecord]:
        with self._lock:
            return sorted(self._tasks.values(), key=lambda t: t.created_at)

    def cancel(self, task_id: str) -> bool:
        """Request cancellation. Returns True when a task was cancelled."""
        with self._lock:
            task = self._tasks.get(task_id)
            if task is None or task.status in (
                TaskStatus.COMPLETED, TaskStatus.CANCELLED, TaskStatus.FAILED
            ):
                return False
            self._cancel_events[task_id].set()
            if task.status is TaskStatus.QUEUED:
                task.status = TaskStatus.CANCELLED
                task.progress = 0.0
                task.message = "cancelled before execution"
                self._touch(task)
                self._persist_unlocked()
            else:
                task.message = "cancellation requested"
                self._touch(task)
            return True

    # -- worker ---------------------------------------------------------

    def _worker_loop(self) -> None:
        while not self._stop.is_set():
            tid = self._queue.get()  # blocking; shutdown enqueues sentinel
            if tid is None:
                return
            try:
                self._run_one(tid)
            finally:
                self._queue.task_done()

    def _run_one(self, tid: str) -> None:
        with self._lock:
            task = self._tasks.get(tid)
            event = self._cancel_events.get(tid)
            if task is None:
                return
            if task.status is TaskStatus.CANCELLED:
                return
            task.status = TaskStatus.RUNNING
            task.progress = 0.0
            task.message = "running"
            self._touch(task)
            self._persist_unlocked()
            params = task.params

        evaluator = SpecificityEvaluator(
            self.library,
            max_mismatches=params.max_mismatches,
            min_three_prime_match=params.min_three_prime_match,
            chunk_size=params.chunk_size,
        )

        def progress_cb(value: float) -> None:
            with self._lock:
                t = self._tasks[tid]
                t.progress = value
                self._touch(t)

        try:
            result = design_primers(
                params,
                evaluator=evaluator,
                progress_cb=progress_cb,
                cancel_cb=event.is_set if event else None,
            )
        except (DesignCancelled, EvalCancelled):
            with self._lock:
                t = self._tasks[tid]
                t.status = TaskStatus.CANCELLED
                t.progress = 0.0
                t.result = None
                t.message = "cancelled"
                self._touch(t)
                self._persist_unlocked()
            return
        except Exception as exc:  # noqa: BLE001 - surfaced as task failure
            with self._lock:
                t = self._tasks[tid]
                t.status = TaskStatus.FAILED
                t.error = f"{type(exc).__name__}: {exc}"
                t.message = "failed"
                self._touch(t)
                self._persist_unlocked()
            return

        with self._lock:
            t = self._tasks[tid]
            if event is not None and event.is_set():
                t.status = TaskStatus.CANCELLED
                t.progress = 0.0
                t.message = "cancelled"
            else:
                t.status = TaskStatus.COMPLETED
                t.progress = 1.0
                t.result = result
                t.message = f"{len(result.pairs)} pair(s)"
            self._touch(t)
            self._persist_unlocked()

    # -- helpers --------------------------------------------------------

    @staticmethod
    def _touch(task: TaskRecord) -> None:
        import time

        task.updated_at = time.time()

    def _persist_unlocked(self) -> None:
        if self.storage is not None:
            self.storage.save_tasks(self._tasks)

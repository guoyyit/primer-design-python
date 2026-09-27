"""Task lifecycle management: progress, cancellation, crash recovery.

Each task lives in one JSON record written via temp-file + atomic
replace, so a cancelled or crashed task can never expose a partially
written final result.  On restart, tasks that were PENDING/RUNNING are
recovered as FAILED("interrupted") with no result; only records that
reached COMPLETED before the crash are returned with their result.
"""
from __future__ import annotations

import json
import os
import tempfile
import threading
import time
import uuid
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace

from .models import (
    DesignParams,
    OffTargetHit,
    Primer,
    PrimerPair,
    StructureReport,
    TaskRecord,
    TaskStatus,
    TemplateTopology,
)

ProgressCb = Callable[[int, int], None]
TaskFn = Callable[[TaskRecord, threading.Event, ProgressCb], list[PrimerPair]]


# --------------------------------------------------------------------------
# JSON serialization (models stay behaviour-free dataclasses)
def params_to_dict(p: DesignParams) -> dict:
    d = dict(p.__dict__)
    d["topology"] = p.topology.value
    return d


def params_from_dict(d: dict) -> DesignParams:
    d = dict(d)
    d["topology"] = TemplateTopology(d["topology"])
    return DesignParams(**d)


def pair_to_dict(pair: PrimerPair) -> dict:
    def primer_d(p: Primer) -> dict:
        return {
            "sequence": p.sequence,
            "start": p.start,
            "length": p.length,
            "is_forward": p.is_forward,
            "report": p.report.__dict__,
        }

    def primer_o(d: dict) -> Primer:
        return Primer(
            d["sequence"], d["start"], d["length"], d["is_forward"],
            StructureReport(**d["report"]),
        )

    return {
        "forward": primer_d(pair.forward),
        "reverse": primer_d(pair.reverse),
        "product_size": pair.product_size,
        "heterodimer_dg": pair.heterodimer_dg,
        "off_targets": [h.__dict__ for h in pair.off_targets],
        "penalty": pair.penalty,
        "rank_score": pair.rank_score,
    }


def pair_from_dict(d: dict) -> PrimerPair:
    def primer_o(q: dict) -> Primer:
        return Primer(
            q["sequence"], q["start"], q["length"], q["is_forward"],
            StructureReport(**q["report"]),
        )

    return PrimerPair(
        primer_o(d["forward"]),
        primer_o(d["reverse"]),
        d["product_size"],
        d["heterodimer_dg"],
        tuple(OffTargetHit(**h) for h in d["off_targets"]),
        d["penalty"],
        d["rank_score"],
    )


def record_to_dict(r: TaskRecord) -> dict:
    return {
        "task_id": r.task_id,
        "status": r.status.value,
        "progress_done": r.progress_done,
        "progress_total": r.progress_total,
        "error": r.error,
        "params": params_to_dict(r.params) if r.params else None,
        "result": [pair_to_dict(p) for p in r.result] if r.result is not None else None,
    }


def record_from_dict(d: dict) -> TaskRecord:
    result = None
    if d.get("result") is not None:
        result = [pair_from_dict(p) for p in d["result"]]
    return TaskRecord(
        task_id=d["task_id"],
        status=TaskStatus(d["status"]),
        params=params_from_dict(d["params"]) if d.get("params") else None,
        progress_done=d.get("progress_done", 0),
        progress_total=d.get("progress_total", 0),
        result=result,
        error=d.get("error"),
    )


# --------------------------------------------------------------------------
class TaskManager:
    def __init__(self, workdir: str, max_workers: int = 2) -> None:
        self.workdir = workdir
        self.tasks_dir = os.path.join(workdir, "tasks")
        os.makedirs(self.tasks_dir, exist_ok=True)
        self._lock = threading.RLock()
        self._records: dict[str, TaskRecord] = {}
        self._cancel: dict[str, threading.Event] = {}
        self._last_persist: dict[str, float] = {}
        self._executor = ThreadPoolExecutor(max_workers=max_workers)
        self.recover()

    # -- persistence ------------------------------------------------------
    def _path(self, task_id: str) -> str:
        return os.path.join(self.tasks_dir, f"{task_id}.json")

    def _persist(self, record: TaskRecord) -> None:
        path = self._path(record.task_id)
        directory = os.path.dirname(path)
        fd, tmp = tempfile.mkstemp(prefix=".task-", suffix=".tmp", dir=directory)
        try:
            with os.fdopen(fd, "w") as fh:
                json.dump(record_to_dict(record), fh, ensure_ascii=False)
                fh.flush()
                os.fsync(fh.fileno())
            os.replace(tmp, path)
        except BaseException:
            try:
                os.unlink(tmp)
            except FileNotFoundError:
                pass
            raise

    def recover(self) -> None:
        """Load task records from disk; abort interrupted tasks."""
        with self._lock:
            self._records = {}
            for name in sorted(os.listdir(self.tasks_dir)):
                if not name.endswith(".json"):
                    continue
                with open(os.path.join(self.tasks_dir, name), encoding="utf-8") as fh:
                    record = record_from_dict(json.load(fh))
                if record.status in (TaskStatus.PENDING, TaskStatus.RUNNING):
                    record = replace(
                        record,
                        status=TaskStatus.FAILED,
                        error="interrupted by service restart; please resubmit",
                        result=None,
                    )
                    self._persist(record)
                self._records[record.task_id] = record

    # -- API --------------------------------------------------------------
    def submit(self, params: DesignParams, fn: TaskFn) -> str:
        task_id = uuid.uuid4().hex
        record = TaskRecord(task_id=task_id, status=TaskStatus.PENDING, params=params)
        event = threading.Event()
        with self._lock:
            self._records[task_id] = record
            self._cancel[task_id] = event
            self._persist(record)
        self._executor.submit(self._run, record, event, fn)
        return task_id

    def _run(self, record: TaskRecord, event: threading.Event, fn: TaskFn) -> None:
        task_id = record.task_id
        with self._lock:
            record = replace(record, status=TaskStatus.RUNNING)
            self._records[task_id] = record
            self._persist(record)

        def progress(done: int, total: int) -> None:
            with self._lock:
                cur = self._records[task_id]
                updated = replace(cur, progress_done=done, progress_total=total)
                self._records[task_id] = updated
                # in-memory state is always fresh; throttle fsync-heavy
                # persistence so fine-grained progress stays cheap
                now = time.monotonic()
                if now - self._last_persist.get(task_id, 0.0) >= 0.25 or done >= total:
                    self._persist(updated)
                    self._last_persist[task_id] = now

        try:
            if event.is_set():
                raise TaskCancelled()
            result = fn(record, event, progress)
            if event.is_set():
                raise TaskCancelled()
            with self._lock:
                done = replace(
                    self._records[task_id],
                    status=TaskStatus.COMPLETED,
                    result=result,
                    progress_done=self._records[task_id].progress_total or 1,
                    progress_total=max(self._records[task_id].progress_total, 1),
                )
                self._records[task_id] = done
                self._persist(done)
        except TaskCancelled:
            with self._lock:
                cancelled = replace(
                    self._records[task_id],
                    status=TaskStatus.CANCELLED,
                    result=None,
                    error="cancelled",
                )
                self._records[task_id] = cancelled
                self._persist(cancelled)
        except BaseException as exc:  # noqa: BLE001 - record every failure
            with self._lock:
                if event.is_set():
                    cancelled = replace(
                        self._records[task_id],
                        status=TaskStatus.CANCELLED,
                        result=None,
                        error="cancelled",
                    )
                    self._records[task_id] = cancelled
                    self._persist(cancelled)
                else:
                    failed = replace(
                        self._records[task_id],
                        status=TaskStatus.FAILED,
                        result=None,
                        error=f"{type(exc).__name__}: {exc}",
                    )
                    self._records[task_id] = failed
                    self._persist(failed)
        finally:
            with self._lock:
                self._cancel.pop(task_id, None)

    def cancel(self, task_id: str) -> bool:
        with self._lock:
            record = self._records.get(task_id)
            event = self._cancel.get(task_id)
            if record is None:
                return False
            if record.status in (TaskStatus.COMPLETED, TaskStatus.CANCELLED, TaskStatus.FAILED):
                return False
            if event is not None:
                event.set()
            else:
                rec = replace(record, status=TaskStatus.CANCELLED, result=None, error="cancelled")
                self._records[task_id] = rec
                self._persist(rec)
            return True

    def get(self, task_id: str) -> TaskRecord:
        with self._lock:
            if task_id not in self._records:
                raise KeyError(task_id)
            return self._records[task_id]

    def list_tasks(self) -> list[TaskRecord]:
        with self._lock:
            return list(self._records.values())

    def shutdown(self) -> None:
        self._executor.shutdown(wait=False, cancel_futures=True)


class TaskCancelled(Exception):
    """Raised by task functions to signal cooperative cancellation."""

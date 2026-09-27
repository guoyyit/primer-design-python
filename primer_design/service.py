"""Service facade wiring reference library, task manager and storage."""
from __future__ import annotations

import threading

from .models import DesignParams, ReferenceSequence, TaskRecord
from .specificity import ReferenceLibrary
from .storage import Storage
from .tasks import TaskManager


class PrimerDesignService:
    def __init__(self, storage_dir: str | None = None, max_workers: int = 1) -> None:
        self.storage = Storage(storage_dir) if storage_dir else None
        if self.storage is not None:
            self.storage.cleanup_temp()
            self.library = self.storage.load_library() or ReferenceLibrary()
        else:
            self.library = ReferenceLibrary()
        self.manager = TaskManager(self.library, self.storage, max_workers=max_workers)
        # Make all persisted task records queryable immediately; unfinished
        # ones are re-queued by start()/recover().
        self.manager.load_persisted()
        self._lib_lock = threading.RLock()

    def start(self) -> None:
        # Workers must exist before re-queued tasks are published, so no
        # recovered task can be missed by a not-yet-polling worker.
        self.manager.start()
        self.manager.recover()

    def shutdown(self) -> None:
        self.manager.shutdown()

    def submit_design(self, params: DesignParams) -> TaskRecord:
        return self.manager.submit(params)

    def get_task(self, task_id: str) -> TaskRecord | None:
        return self.manager.get(task_id)

    def list_tasks(self) -> list[TaskRecord]:
        return self.manager.list_tasks()

    def cancel_task(self, task_id: str) -> bool:
        return self.manager.cancel(task_id)

    # -- reference library ---------------------------------------------

    def add_references(self, refs: list[ReferenceSequence]) -> int:
        with self._lib_lock:
            version = self.library.add(refs)
            self._persist_library()
        return version.version

    def replace_references(self, refs: list[ReferenceSequence]) -> int:
        with self._lib_lock:
            version = self.library.replace(refs)
            self._persist_library()
        return version.version

    def clear_references(self) -> None:
        with self._lib_lock:
            self.library.clear()
            self._persist_library()

    def library_info(self) -> dict:
        v = self.library.version
        return {"version": v.version, "size": v.size, "created_at": v.created_at}

    def _persist_library(self) -> None:
        if self.storage is not None:
            self.storage.save_library(self.library)

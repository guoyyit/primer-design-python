"""Crash-safe persistence for tasks and the reference library.

Writes are atomic: serialize to a unique ``.tmp`` file, fsync it, then
``os.replace`` onto the final path. A crash mid-write therefore leaves only
a stale ``.tmp`` file which is discarded at startup -- the previous final
file (if any) remains intact. Completed task results are only ever persisted
inside a task whose status is ``completed``; partial progress of a running
task is not serialized, so intermediate results can never be read as final.
"""
from __future__ import annotations

import json
import os
import tempfile
import time

from .models import (
    DesignParams,
    DesignResult,
    OffTargetHit,
    PrimerMetrics,
    PrimerPair,
    Strand,
    TaskRecord,
    TaskStatus,
)
from .specificity import ReferenceLibrary

TASKS_FILE = "tasks.json"
LIBRARY_FILE = "library.json"


class Storage:
    TMP_SUFFIX = ".partial"

    def __init__(self, directory: str) -> None:
        self.directory = directory
        os.makedirs(directory, exist_ok=True)

    # -- atomic primitives ---------------------------------------------

    def _path(self, name: str) -> str:
        return os.path.join(self.directory, name)

    def _atomic_write_json(self, name: str, payload) -> None:
        # Retry: a concurrent startup cleanup in another process could in
        # rare races unlink our temp file before os.replace.
        for attempt in range(3):
            final = self._path(name)
            fd, tmp = tempfile.mkstemp(prefix="." + name + ".",
                                       suffix=self.TMP_SUFFIX,
                                       dir=self.directory)
            try:
                with os.fdopen(fd, "w", encoding="utf-8") as fh:
                    json.dump(payload, fh, ensure_ascii=False,
                              separators=(",", ":"))
                    fh.flush()
                    os.fsync(fh.fileno())
                os.replace(tmp, final)
                return
            except FileNotFoundError:
                if attempt == 2:
                    raise
                time.sleep(0.05 * (attempt + 1))
            except BaseException:
                try:
                    os.unlink(tmp)
                except OSError:
                    pass
                raise

    def _read_json(self, name: str):
        path = self._path(name)
        if not os.path.exists(path):
            return None
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)

    def cleanup_temp(self, max_age_seconds: float = 10.0) -> list[str]:
        """Remove stale partial writes older than ``max_age_seconds``.

        Recent files are left untouched: another (live) process may be
        writing them right now, and deleting an in-flight temp file would
        make its ``os.replace`` fail. Truly stale files come from crashes.
        """
        removed = []
        cutoff = time.time() - max_age_seconds
        for entry in os.listdir(self.directory):
            if not entry.endswith(self.TMP_SUFFIX):
                continue
            path = self._path(entry)
            try:
                if os.path.getmtime(path) > cutoff:
                    continue
                os.unlink(path)
                removed.append(entry)
            except OSError:
                pass
        return removed

    # -- tasks ----------------------------------------------------------

    def save_task(self, task: TaskRecord) -> None:
        all_tasks = self.load_tasks_raw()
        all_tasks[task.task_id] = self._dump_task(task)
        self._atomic_write_json(TASKS_FILE, all_tasks)

    def save_tasks(self, tasks: dict[str, TaskRecord]) -> None:
        self._atomic_write_json(
            TASKS_FILE, {tid: self._dump_task(t) for tid, t in tasks.items()}
        )

    def load_tasks_raw(self) -> dict:
        return self._read_json(TASKS_FILE) or {}

    def load_tasks(self) -> dict[str, TaskRecord]:
        raw = self._read_json(TASKS_FILE)
        if not raw:
            return {}
        out: dict[str, TaskRecord] = {}
        for tid, blob in raw.items():
            try:
                out[tid] = self._load_task(blob)
            except (KeyError, TypeError, ValueError):
                # A corrupt record must not prevent the rest from recovering.
                continue
        return out

    def _dump_task(self, task: TaskRecord) -> dict:
        blob = {
            "task_id": task.task_id,
            "params": vars(task.params).copy(),
            "status": task.status.value,
            "created_at": task.created_at,
            "updated_at": task.updated_at,
            "error": task.error,
            "message": task.message,
        }
        if task.status is TaskStatus.COMPLETED and task.result is not None:
            blob["progress"] = 1.0
            blob["result"] = self._dump_result(task.result)
        else:
            # Never persist intermediate results of unfinished tasks.
            blob["progress"] = task.progress if task.status is TaskStatus.CANCELLED else 0.0
        return blob

    def _load_task(self, blob: dict) -> TaskRecord:
        params = DesignParams(**blob["params"])
        status = TaskStatus(blob["status"])
        task = TaskRecord(
            task_id=blob["task_id"],
            params=params,
            status=status,
            progress=blob.get("progress", 0.0),
            message=blob.get("message", ""),
            created_at=blob["created_at"],
            updated_at=blob.get("updated_at", blob["created_at"]),
            error=blob.get("error"),
        )
        if status is TaskStatus.COMPLETED and "result" in blob:
            task.result = self._load_result(blob["result"])
        return task

    def _dump_result(self, result: DesignResult) -> dict:
        return {
            "pairs": [self._dump_pair(p) for p in result.pairs],
            "evaluated_pairs": result.evaluated_pairs,
            "generated_at": result.generated_at,
        }

    def _load_result(self, blob: dict) -> DesignResult:
        return DesignResult(
            pairs=tuple(self._load_pair(p) for p in blob["pairs"]),
            evaluated_pairs=blob["evaluated_pairs"],
            generated_at=blob["generated_at"],
        )

    def _dump_pair(self, pair: PrimerPair) -> dict:
        return {
            "pair_id": pair.pair_id,
            "forward": self._dump_metrics(pair.forward),
            "reverse": self._dump_metrics(pair.reverse),
            "product_len": pair.product_len,
            "heterodimer_dg": pair.heterodimer_dg,
            "score": pair.score,
            "rank": pair.rank,
            "off_targets": [vars(h) | {"strand": h.strand.value} for h in pair.off_targets],
            "penalty_detail": pair.penalty_detail,
        }

    def _load_pair(self, blob: dict) -> PrimerPair:
        hits = tuple(
            OffTargetHit(
                reference_id=h["reference_id"],
                strand=Strand(h["strand"]),
                position=h["position"],
                mismatches=h["mismatches"],
                three_prime_match=h["three_prime_match"],
                aligned_primer=h["aligned_primer"],
                aligned_reference=h["aligned_reference"],
            )
            for h in blob.get("off_targets", ())
        )
        return PrimerPair(
            pair_id=blob["pair_id"],
            forward=self._load_metrics(blob["forward"]),
            reverse=self._load_metrics(blob["reverse"]),
            product_len=blob["product_len"],
            heterodimer_dg=blob["heterodimer_dg"],
            score=blob["score"],
            rank=blob["rank"],
            off_targets=hits,
            penalty_detail=blob.get("penalty_detail", {}),
        )

    def _dump_metrics(self, m: PrimerMetrics) -> dict:
        return vars(m) | {"strand": m.strand.value}

    def _load_metrics(self, blob: dict) -> PrimerMetrics:
        return PrimerMetrics(
            sequence=blob["sequence"],
            strand=Strand(blob["strand"]),
            position=blob["position"],
            length=blob["length"],
            gc_content=blob["gc_content"],
            tm=blob["tm"],
            self_dimer_dg=blob["self_dimer_dg"],
            hairpin_dg=blob["hairpin_dg"],
            homopolymer_run=blob["homopolymer_run"],
            dinucleotide_run=blob["dinucleotide_run"],
            terminal_palindrome=blob["terminal_palindrome"],
            penalties=blob.get("penalties", {}),
        )

    # -- library --------------------------------------------------------

    def save_library(self, library: ReferenceLibrary) -> None:
        self._atomic_write_json(LIBRARY_FILE, library.dump())

    def load_library(self) -> ReferenceLibrary | None:
        payload = self._read_json(LIBRARY_FILE)
        if payload is None:
            return None
        return ReferenceLibrary.load(payload)

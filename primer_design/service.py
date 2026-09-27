"""End-to-end orchestration: search -> filter -> specificity -> rank.

Specificity is evaluated on an immutable library snapshot taken when the
task starts, so concurrent library updates never produce mixed state.

Candidate pairs are streamed through a bounded top-K heap (memory stays
O(K) regardless of how many pairs the template generates).  The
heterodimer term and the off-target scan are computed only for the
window; the window is then *widened until provably exact*: because both
the heterodimer penalty and the off-target penalty are >= 0, any pair
outside the window has final penalty >= its pre-penalty, so once the
K-th smallest pre-penalty is >= the current max_pairs-th best final
penalty, no outside pair can enter the final top-max_pairs ranking.
"""
from __future__ import annotations

import heapq
import os
import threading
from dataclasses import replace

from .reflib import ReferenceLibrary
from .scoring import (
    W_OFFTARGET,
    heterodimer_penalty,
    pair_passes,
    pair_pre_penalty,
    passes_filters,
    pre_penalty,
    primer_penalty,
    rank_pairs,
)
from .search import generate_primers, iter_pairs
from .specificity import evaluate_primer
from .thermo import heterodimer_dg
from .tasks import TaskCancelled, TaskManager
from .models import DesignParams, OffTargetHit, PrimerPair

SPECIFICITY_PREFACTOR = 10
MIN_WINDOW = 100
MAX_WINDOW = 10_000  # beyond this the top-n is exact within the window only
CHUNK_SIZE = 100_000


class PrimerDesignService:
    def __init__(self, workdir: str) -> None:
        self.workdir = workdir
        self.lib_dir = os.path.join(workdir, "libraries")
        os.makedirs(self.lib_dir, exist_ok=True)
        self._lib_lock = threading.RLock()
        self._libraries: dict[str, ReferenceLibrary] = {}
        for name in os.listdir(self.lib_dir):
            if name.endswith(".json"):
                self._libraries[name[:-5]] = ReferenceLibrary.load(
                    os.path.join(self.lib_dir, name)
                )
        self.tasks = TaskManager(workdir)

    # -- reference libraries ---------------------------------------------
    def _lib_path(self, name: str) -> str:
        if not name.replace("_", "").replace("-", "").isalnum():
            raise ValueError("invalid library name")
        return os.path.join(self.lib_dir, f"{name}.json")

    def library(self, name: str) -> ReferenceLibrary:
        with self._lib_lock:
            lib = self._libraries.get(name)
            if lib is None:
                raise KeyError(f"unknown reference library: {name}")
            return lib

    def add_references(self, name: str, refs) -> None:
        with self._lib_lock:
            path = self._lib_path(name)
            lib = self._libraries.get(name) or ReferenceLibrary()
            lib.add_many(refs)
            lib.save(path)
            self._libraries[name] = lib

    def remove_references(self, name: str, ref_ids) -> None:
        with self._lib_lock:
            lib = self.library(name)
            lib.remove_many(ref_ids)
            lib.save(self._lib_path(name))

    # -- task submission --------------------------------------------------
    def submit_design(self, params: DesignParams) -> str:
        if params.reference_library is not None:
            self.library(params.reference_library)  # fail fast if unknown
        if params.target_end is not None:
            n_tpl = len(params.template)
            ts, te = params.target_start, params.target_end
            width = te - ts if te > ts else (n_tpl - ts + te if n_tpl else 0)
            if width > params.product_size_max:
                raise ValueError(
                    "target interval is wider than product_size_max; "
                    "no amplicon can cover it"
                )
        return self.tasks.submit(params, self._run_design)

    def get_task(self, task_id: str):
        return self.tasks.get(task_id)

    def cancel_task(self, task_id: str) -> bool:
        return self.tasks.cancel(task_id)

    # -- pipeline ---------------------------------------------------------
    def _run_design(self, record, cancel_event, progress) -> list[PrimerPair]:
        params: DesignParams = record.params
        # Phase 1: generation + structure filtering
        progress(0, 1)
        all_primers = generate_primers(params)
        good = [p for p in all_primers if passes_filters(p, params)]
        progress(1, 1)
        if not good:
            return []
        if cancel_event.is_set():
            raise TaskCancelled()

        snapshot = (
            self.library(params.reference_library).snapshot()
            if params.reference_library is not None
            else None
        )
        hetero_cache: dict[tuple[str, str], float] = {}
        hits_by_seq: dict[str, tuple[OffTargetHit, ...]] = {}
        pen_cache: dict[str, float] = {}

        def ppen(seq: str, report) -> float:
            pen = pen_cache.get(seq)
            if pen is None:
                pen = pen_cache[seq] = primer_penalty(report, params)
            return pen

        n = params.max_pairs
        k = max(MIN_WINDOW, n * SPECIFICITY_PREFACTOR)

        while True:
            # Stream pairs through a bounded heap of the k smallest pre-penalties
            heap: list = []
            total = 0
            counter = 0
            for pair in iter_pairs(good, params, heterodimer=False):
                f, r = pair.forward, pair.reverse
                pre = pair_pre_penalty(
                    ppen(f.sequence, f.report),
                    ppen(r.sequence, r.report),
                    f.report.tm,
                    r.report.tm,
                    pair.product_size,
                    params,
                )
                counter += 1
                item = (-pre, counter, pair)
                if len(heap) < k:
                    heapq.heappush(heap, item)
                elif pre < -heap[0][0]:
                    heapq.heapreplace(heap, item)
                total += 1
                if total % 20_000 == 0 and cancel_event.is_set():
                    raise TaskCancelled()
            if total == 0:
                return []
            window = [it[2] for it in heap]
            pre_k = -heap[0][0] if total > k else float("inf")

            # Heterodimer (memoized across expansions), then constraint filter
            scored: list[tuple[PrimerPair, float, float]] = []
            for pair in window:
                key = (pair.forward.sequence, pair.reverse.sequence)
                dg = hetero_cache.get(key)
                if dg is None:
                    dg = hetero_cache[key] = heterodimer_dg(*key)
                if dg < params.max_dimer_dg:
                    continue
                full = pre_penalty(pair, params) + heterodimer_penalty(dg)
                scored.append((pair, dg, full))

            # Off-targets on the committed snapshot (cached per sequence)
            if snapshot is not None:
                distinct = sorted({
                    seq
                    for pair, _, _ in scored
                    for seq in (pair.forward.sequence, pair.reverse.sequence)
                    if seq not in hits_by_seq
                })
                for i, seq in enumerate(distinct, start=1):
                    hits: dict[tuple[str, int], OffTargetHit] = {}
                    for ref_id, ref_seq in snapshot.references():
                        for hit in evaluate_primer(
                            seq,
                            [(ref_id, ref_seq)],
                            params.max_mismatches,
                            params.max_3p_consecutive,
                            CHUNK_SIZE,
                        ):
                            hits.setdefault((hit.reference_id, hit.position), hit)
                        progress(i, max(1, len(distinct)))
                        if cancel_event.is_set():
                            raise TaskCancelled()
                    hits_by_seq[seq] = tuple(hits.values())

            # Final penalties + exactness check
            finals: list[tuple[float, PrimerPair, float]] = []
            for pair, dg, full in scored:
                merged: dict[tuple[str, int], OffTargetHit] = {}
                for seq in (pair.forward.sequence, pair.reverse.sequence):
                    for hit in hits_by_seq.get(seq, ()):
                        merged.setdefault((hit.reference_id, hit.position), hit)
                finals.append((full + W_OFFTARGET * len(merged), pair, dg))
            finals.sort(key=lambda x: x[0])
            if total <= k or len(finals) < n:
                break
            cutoff = finals[n - 1][0]  # n-th best final penalty
            if pre_k >= cutoff:
                break  # exact: any outside pair costs at least pre_k
            if k >= MAX_WINDOW:
                break  # degenerate input: exact within the window only
            k = min(k * 4, MAX_WINDOW)

        enriched = []
        for _final, pair, dg in finals:
            merged: dict[tuple[str, int], OffTargetHit] = {}
            for seq in (pair.forward.sequence, pair.reverse.sequence):
                for hit in hits_by_seq.get(seq, ()):
                    merged.setdefault((hit.reference_id, hit.position), hit)
            enriched.append(replace(
                pair,
                heterodimer_dg=dg,
                off_targets=tuple(merged.values()),
            ))
        progress(1, 1)
        return rank_pairs(enriched, params)

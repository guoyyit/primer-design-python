"""Off-target binding evaluation against a reference sequence library.

Concurrency model
-----------------
:class:`ReferenceLibrary` keeps its contents in *immutable snapshots*. An
update (add/replace/clear) builds a brand-new tuple and swaps it in under a
lock. An evaluation that has taken a snapshot keeps reading the same tuple
until it finishes, so an update arriving mid-evaluation can never produce a
mixed new/old view.

Search algorithm
----------------
Direct window scan: every length-L window of a reference strand is compared
to the primer; a hit is recorded when mismatches <= ``max_mismatches`` AND
the consecutive matched run from the primer 3' end is >=
``min_three_prime_match``. The scan is O(L*R) per primer; long references
are processed in ``chunk_size`` blocks (constant extra memory) and progress
is reported per block. This deliberately favors correctness and bounded
memory over building a large k-mer index.
"""
from __future__ import annotations

import time
from .models import OffTargetHit, ReferenceSequence, SnapshotVersion, Strand
from .sequence import normalize_sequence, reverse_complement

_MAX_HITS_PER_PRIMER = 100


class ReferenceLibrary:
    """Versioned, append/replace-able store of reference sequences."""

    def __init__(self) -> None:
        import threading
        self._lock = threading.RLock()
        self._refs: tuple[ReferenceSequence, ...] = ()
        self._version = 0

    def add(self, refs: list[ReferenceSequence]) -> SnapshotVersion:
        norm = [self._normalize(r) for r in refs]
        with self._lock:
            existing = {r.ref_id for r in self._refs}
            merged = list(self._refs)
            for r in norm:
                if r.ref_id in existing:
                    raise ValueError(f"duplicate reference id: {r.ref_id}")
                merged.append(r)
                existing.add(r.ref_id)
            return self._publish(merged)

    def replace(self, refs: list[ReferenceSequence]) -> SnapshotVersion:
        norm = [self._normalize(r) for r in refs]
        ids = [r.ref_id for r in norm]
        if len(set(ids)) != len(ids):
            raise ValueError("duplicate reference id in replacement batch")
        with self._lock:
            return self._publish(norm)

    def clear(self) -> SnapshotVersion:
        with self._lock:
            return self._publish([])

    def _publish(self, refs: list[ReferenceSequence]) -> SnapshotVersion:
        self._version += 1
        self._refs = tuple(refs)
        return SnapshotVersion(self._version, time.time(), len(refs))

    @staticmethod
    def _normalize(ref: ReferenceSequence) -> ReferenceSequence:
        return ReferenceSequence(
            ref_id=ref.ref_id,
            sequence=normalize_sequence(ref.sequence),
            description=ref.description,
        )

    @property
    def version(self) -> SnapshotVersion:
        with self._lock:
            return SnapshotVersion(self._version, time.time(), len(self._refs))

    def snapshot(self) -> tuple[tuple[ReferenceSequence, ...], int]:
        with self._lock:
            return self._refs, self._version

    def dump(self) -> dict:
        with self._lock:
            return {
                "version": self._version,
                "references": [
                    {"ref_id": r.ref_id, "sequence": r.sequence,
                     "description": r.description}
                    for r in self._refs
                ],
            }

    @classmethod
    def load(cls, payload: dict) -> "ReferenceLibrary":
        lib = cls()
        refs = [
            ReferenceSequence(r["ref_id"], r["sequence"], r.get("description", ""))
            for r in payload.get("references", [])
        ]
        with lib._lock:
            lib._refs = tuple(refs)
            lib._version = int(payload.get("version", 0))
        return lib


class SpecificityEvaluator:
    def __init__(
        self,
        library: ReferenceLibrary,
        max_mismatches: int = 2,
        min_three_prime_match: int = 8,
        chunk_size: int = 1000,
    ) -> None:
        self.library = library
        self.max_mismatches = max(0, max_mismatches)
        self.min_three_prime_match = max(1, min_three_prime_match)
        self.chunk_size = max(64, chunk_size)

    # -- single-reference scan ------------------------------------------

    @staticmethod
    def count_hits(
        seq: str,
        ref: str,
        max_mismatches: int,
        min_three_prime_match: int,
    ) -> list[tuple[int, int, int]]:
        """Scan ``seq`` (5'->3') over ``ref`` (+ strand, 5'->3').

        Returns sorted ``(position, mismatches, 3' match)`` hits.
        """
        hits: list[tuple[int, int, int]] = []
        L, R = len(seq), len(ref)
        if L == 0 or R < L:
            return hits
        for start in range(R - L + 1):
            mm = 0
            for i in range(L):
                if seq[i] != ref[start + i]:
                    mm += 1
                    if mm > max_mismatches:
                        break
            if mm > max_mismatches:
                continue
            tpm = 0
            for i in range(L - 1, -1, -1):
                if seq[i] == ref[start + i]:
                    tpm += 1
                else:
                    break
            if tpm >= min_three_prime_match:
                hits.append((start, mm, tpm))
        return hits

    # -- whole-library evaluation ---------------------------------------

    def evaluate_pair(
        self,
        fwd_seq: str,
        rev_seq: str,
        progress_cb=None,
        cancel_cb=None,
    ) -> tuple[OffTargetHit, ...]:
        # One snapshot for the whole evaluation: updates during the run are
        # invisible and the two primers always see the same library version.
        refs, _version = self.library.snapshot()
        jobs = (
            (fwd_seq, Strand.FORWARD, False),
            (rev_seq, Strand.REVERSE, True),
        )
        out: list[OffTargetHit] = []
        total_bp = max(1, sum(len(r.sequence) for r in refs) * 2)
        done_bp = 0

        def scan_ref(seq: str, ref_seq: str):
            """Yield (global_pos, mm, tpm, aligned_segment) hits, scanning
            the reference in bounded-memory chunks; neighboring chunks
            overlap by L-1 so windows spanning a block boundary are kept."""
            L, R = len(seq), len(ref_seq)
            if R < L:
                return
            step = max(self.chunk_size, L)
            for block_start in range(0, R - L + 1, step):
                if cancel_cb is not None and cancel_cb():
                    raise CancelledError("specificity evaluation cancelled")
                block_end = min(R, block_start + step + L - 1)
                segment = ref_seq[block_start:block_end]
                for local, mm, tpm in self.count_hits(
                    seq, segment, self.max_mismatches,
                    self.min_three_prime_match,
                ):
                    yield (block_start + local, mm, tpm,
                           segment[local:local + L])
                # only windows whose *start* is in [block_start, block_end-L]
                # are produced by this block; account progress accordingly.
                consumed = max(1, block_end - block_start - (L - 1))
                yield (consumed,)  # progress marker

        for primer, strand, is_reverse in jobs:
            seq = normalize_sequence(primer)
            for ref in refs:
                target = (
                    reverse_complement(ref.sequence) if is_reverse
                    else ref.sequence
                )
                hits_for_ref = 0
                for item in scan_ref(seq, target):
                    if len(item) == 1:
                        done_bp += item[0]
                        if progress_cb is not None:
                            progress_cb(min(1.0, done_bp / total_bp))
                        continue
                    if hits_for_ref >= _MAX_HITS_PER_PRIMER:
                        continue
                    pos, mm, tpm, aligned = item
                    hits_for_ref += 1
                    if is_reverse:
                        # window at pos on target(=revcomp ref) maps to a
                        # primer 5' end at R-pos-L on the reference + strand.
                        coordinate = len(ref.sequence) - pos - len(seq)
                    else:
                        coordinate = pos
                    out.append(
                        OffTargetHit(
                            reference_id=ref.ref_id,
                            strand=strand,
                            position=coordinate,
                            mismatches=mm,
                            three_prime_match=tpm,
                            aligned_primer=seq,
                            aligned_reference=aligned,
                        )
                    )

        out.sort(key=lambda h: (h.reference_id, h.strand.value, h.position))
        return tuple(out)


class CancelledError(Exception):
    """Raised when an evaluation observes a cancellation request."""

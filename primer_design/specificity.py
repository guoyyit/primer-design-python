"""Off-target evaluation of a primer against reference sequences.

An alignment is reported as a hit when *both* hold:

* total mismatches <= ``max_mismatches``
* 3'-end consecutive matches > ``max_3p_consecutive`` (strictly greater;
  a 3' run exactly equal to the limit is considered safe)

Both strands of every reference are scanned.  References are processed
in overlapping windows of ``chunk_size`` bases so memory stays bounded
regardless of reference length.
"""
from __future__ import annotations

from collections.abc import Iterable

from .models import OffTargetHit
from .sequtil import revcomp


def _count(primer: str, window: str) -> tuple[int, int]:
    """(total mismatches, 3'-end consecutive matches) of one alignment."""
    mismatches = 0
    for a, b in zip(primer, window):
        if a != b:
            mismatches += 1
    run = 0
    for a, b in zip(reversed(primer), reversed(window)):
        if a != b:
            break
        run += 1
    return mismatches, run


def _scan_strand(
    primer: str,
    ref_id: str,
    strand_seq: str,
    strand: str,
    ref_len: int,
    max_mismatches: int,
    max_3p_consecutive: int,
    chunk_size: int,
) -> list[OffTargetHit]:
    k = len(primer)
    hits: list[OffTargetHit] = []
    step = max(1, chunk_size - k + 1)
    for wstart in range(0, max(1, len(strand_seq) - k + 1), step):
        window = strand_seq[wstart : wstart + chunk_size]
        for i in range(len(window) - k + 1):
            mism, run3 = _count(primer, window[i : i + k])
            if mism <= max_mismatches and run3 > max_3p_consecutive:
                pos = wstart + i
                if strand == "-":
                    pos = ref_len - pos - k
                hits.append(OffTargetHit(ref_id, pos, strand, mism, run3))
    return hits


def evaluate_primer(
    seq: str,
    references: Iterable[tuple[str, str]],
    max_mismatches: int,
    max_3p_consecutive: int,
    chunk_size: int = 100_000,
) -> list[OffTargetHit]:
    primer = seq.upper()
    best: dict[tuple[str, int], OffTargetHit] = {}
    for ref_id, ref_seq in references:
        ref = ref_seq.upper()
        raw = _scan_strand(primer, ref_id, ref, "+", len(ref),
                           max_mismatches, max_3p_consecutive, chunk_size)
        raw += _scan_strand(primer, ref_id, revcomp(ref), "-", len(ref),
                            max_mismatches, max_3p_consecutive, chunk_size)
        for hit in raw:
            # A palindromic site is seen on both strands at the same
            # forward coordinate; keep the more dangerous report only.
            key = (hit.reference_id, hit.position)
            old = best.get(key)
            if old is None or (hit.mismatches, -hit.three_prime_consecutive) < (
                old.mismatches, -old.three_prime_consecutive
            ):
                best[key] = hit
    return list(best.values())

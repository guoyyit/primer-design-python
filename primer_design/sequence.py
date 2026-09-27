"""Sequence normalization and coordinate helpers.

Coordinates are 0-based indices along the forward strand of the template.
For circular templates indices wrap around modulo the template length.
"""
from __future__ import annotations

# IUPAC ambiguous base -> complement (kept as-is if it appears; primer search
# only uses unambiguous ACGT windows).
IUPAC_COMPLEMENT = str.maketrans("ACGTNRYSWKMBDHVacgtnryswkmbdhv",
                                 "TGCANYRSWMKVHDBtgcanyrswmkvhdb")
BASES = frozenset("ACGT")


def normalize_sequence(seq: str) -> str:
    """Uppercase, strip whitespace/newlines; validate IUPAC characters."""
    cleaned = "".join(seq.split()).upper()
    if not cleaned:
        raise ValueError("sequence must not be empty")
    bad = {ch for ch in cleaned if ch not in "ACGTNRYSWKMBDHV"}
    if bad:
        raise ValueError(f"unsupported sequence characters: {''.join(sorted(bad))}")
    return cleaned


def reverse_complement(seq: str) -> str:
    return seq.translate(IUPAC_COMPLEMENT)[::-1]


def gc_content(seq: str) -> float:
    if not seq:
        return 0.0
    gc = sum(1 for ch in seq if ch in "GC")
    return gc / len(seq)


def slice_template(template: str, start: int, length: int, circular: bool) -> str:
    """Return ``length`` bases starting at ``start``.

    For circular templates coordinates wrap; for linear templates bases past
    the end are truncated (returns a shorter string or "").
    """
    n = len(template)
    if n == 0 or length <= 0:
        return ""
    if circular:
        return "".join(template[(start + i) % n] for i in range(length))
    if start < 0 or start >= n:
        return ""
    return template[start:start + length]

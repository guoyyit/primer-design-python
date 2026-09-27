"""Small sequence helpers shared across modules."""
from __future__ import annotations

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def normalize_seq(seq: str) -> str:
    """Upper-case and validate a DNA sequence (ACGT only)."""
    s = seq.strip().upper()
    bad = set(s) - set("ACGT")
    if bad:
        raise ValueError(f"invalid bases in sequence: {sorted(bad)}")
    if not s:
        raise ValueError("empty sequence")
    return s


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1].upper()


def is_palindrome(seq: str) -> bool:
    """True if the (double-stranded DNA) sequence equals its reverse complement."""
    return seq.upper() == revcomp(seq)

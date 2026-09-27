"""Self-dimer, heterodimer and hairpin assessment.

These are *approximate* structural scores rather than a full mfold/UNAFold
calculation: every end-overlap alignment of the two strands is enumerated
and scored with nearest-neighbor stacking energies; the most stable
alignment (lowest dG, kcal/mol at 37 C) is reported. A more negative value
means a more stable (more undesirable) structure.

Alignment convention: the top string is written 5'->3' left-to-right; the
bottom string is written 3'->5' left-to-right. Column ``i`` therefore
represents an anti-parallel base pair.
"""
from __future__ import annotations

# Nearest-neighbor deltaG at 37 C, kcal/mol, indexed by the 5'->3'
# dinucleotide of the top (5'->3') strand (SantaLucia 1998).
_PAIR_DG = {
    "AA": -1.00,
    "AT": -0.88,
    "TA": -0.58,
    "CA": -1.45,
    "GT": -1.44,
    "CT": -1.28,
    "GA": -1.30,
    "CG": -2.17,
    "GC": -2.24,
    "GG": -1.84,
}

_COMP = str.maketrans("ACGT", "TGCA")


def _complement(ch: str) -> str:
    return ch.translate(_COMP)


def align_energy(top: str, bottom: str) -> float:
    """Sum stacking energy of two equal-length anti-parallel alignment rows.

    A column only contributes when both bases are paired complements;
    mismatches and gaps (internal bulges are not modeled) contribute 0.
    """
    if len(top) != len(bottom) or len(top) < 2:
        return 0.0
    total = 0.0
    for i in range(len(top) - 1):
        t1, t2 = top[i], top[i + 1]
        b1, b2 = bottom[i], bottom[i + 1]
        if "-" in (t1, t2, b1, b2):
            continue
        if b1 != _complement(t1) or b2 != _complement(t2):
            continue
        total += _PAIR_DG.get(t1 + t2, 0.0)
    return total


def _best_end_alignment(seq_a: str, seq_b: str) -> tuple[float, str, str]:
    """Best dG over every anti-parallel overlap that includes an end.

    Dimer nucleation is driven by 3'-end overlaps, and enumerating all
    relative shifts covers them as a subset (plus internal seed alignments).
    ``k`` is the offset of B's start relative to A's start; the B segment
    overlapping A is shown 3'->5' on the bottom row.
    """
    best = (0.0, "", "")
    la, lb = len(seq_a), len(seq_b)
    for k in range(-(lb - 1), la):
        left = max(0, k)
        right = min(la, lb + k)
        if right - left < 2:
            continue
        top = seq_a[left:right]
        bot_5to3 = seq_b[left - k:right - k]
        bot_3to5 = bot_5to3[::-1]
        dg = align_energy(top, bot_3to5)
        if dg < best[0] - 1e-12:
            best = (dg, top, bot_3to5)
    return best


def best_self_dimer(seq: str) -> tuple[float, str, str]:
    return _best_end_alignment(seq, seq)


def best_heterodimer(seq_a: str, seq_b: str) -> tuple[float, str, str]:
    return _best_end_alignment(seq_a, seq_b)


def best_hairpin(seq: str, min_stem: int = 3, min_loop: int = 3) -> tuple[float, str]:
    """Best (most negative) hairpin stem dG and a textual structure.

    Enumerates outer stem pairs (a, e): stem closes from the outside while
    seq[a+m] complements seq[e-m]; the unpaired loop has length
    e - a - 2k + 1 >= min_loop.
    """
    best_dg = 0.0
    best_struct = ""
    n = len(seq)
    for a in range(n):
        for e in range(a + 1, n):
            k = 0
            while a + k < e - k and _complement(seq[a + k]) == seq[e - k]:
                k += 1
            if k < min_stem:
                continue
            loop_len = e - a - 2 * k + 1
            if loop_len < min_loop:
                continue
            stem5 = seq[a:a + k]
            loop = seq[a + k:e - k + 1]
            stem3 = seq[e - k + 1:e + 1]
            dg = align_energy(stem5, stem3[::-1]) + 0.1 * loop_len
            if dg < best_dg - 1e-12:
                best_dg = dg
                best_struct = f"{stem5}({loop}){stem3}"
    return best_dg, best_struct


def longest_homopolymer(seq: str) -> int:
    if not seq:
        return 0
    best = run = 1
    for i in range(1, len(seq)):
        run = run + 1 if seq[i] == seq[i - 1] else 1
        best = max(best, run)
    return best


def longest_dinucleotide_repeat(seq: str) -> int:
    """Largest count of a consecutively repeated 2-base unit (ATATAT -> 3)."""
    n = len(seq)
    best = 0
    for phase in (0, 1):
        i = phase
        while i + 1 < n:
            unit = seq[i:i + 2]
            count = 0
            j = i
            while j + 1 < n and seq[j:j + 2] == unit:
                count += 1
                j += 2
            best = max(best, count)
            i = j if j > i else i + 2
    return best


def terminal_palindrome(seq: str, min_len: int = 4) -> bool:
    """True if ``min_len`` bases at each end are reverse-complementary."""
    if len(seq) < min_len * 2:
        return False
    five = seq[:min_len]
    three = seq[-min_len:]
    return five == three[::-1].translate(_COMP)

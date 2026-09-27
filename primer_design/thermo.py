"""Thermodynamic and structural calculations for primers.

Model: SantaLucia (1998) nearest-neighbour parameters.

* Tm uses the NN enthalpy/entropy sum with initiation terms, a salt
  correction of 16.6*log10([Na+]) and primer concentration Ct.
* Dimer / hairpin stability is reported as the most negative (most
  stable) contiguous complementary duplex stretch, evaluated as
  dG = dH - T*dS at 37 C.  This is a deliberately simple, fully
  deterministic approximation; it is monotonic in duplex length and
  GC content which is what the ranking and filtering rely on.
"""
from __future__ import annotations

import math

from .models import StructureReport
from .sequtil import is_palindrome, revcomp

# --- SantaLucia 1998 NN parameters (dH kcal/mol, dS cal/mol*K) -------------
_NN: dict[str, tuple[float, float]] = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}
_INIT_AT = (2.3, 4.1)    # terminal A/T pair
_INIT_GC = (0.1, -2.8)   # terminal G/C pair
_SYMMETRY_DS = -1.4      # self-complementary correction

R_CAL = 1.987            # gas constant cal/(mol*K)
CT_M = 250e-9            # primer concentration
NA_M = 0.05              # monovalent salt concentration
DG_TEMP_K = 310.15       # 37 C for structure dG evaluation
MIN_DUPLEX = 2           # shortest complementary stretch scored


def gc_percent(seq: str) -> float:
    seq = seq.upper()
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq)


def _nn_sums(seq: str) -> tuple[float, float]:
    dh = ds = 0.0
    for i in range(len(seq) - 1):
        h, s = _NN[seq[i : i + 2]]
        dh += h
        ds += s
    for end in (seq[0], seq[-1]):
        h, s = _INIT_AT if end in "AT" else _INIT_GC
        dh += h
        ds += s
    if is_palindrome(seq):
        ds += _SYMMETRY_DS
    return dh, ds


def tm_celsius(seq: str) -> float:
    seq = seq.upper()
    if len(seq) < 2:
        return 0.0
    dh, ds = _nn_sums(seq)
    ct = CT_M if is_palindrome(seq) else CT_M / 4.0
    tm_k = (1000.0 * dh) / (ds + R_CAL * math.log(ct))
    return tm_k - 273.15 + 16.6 * math.log10(NA_M)


def annealing_temp(seq: str) -> float:
    """Rule of thumb: anneal 5 C below the primer Tm."""
    return tm_celsius(seq) - 5.0


def _step_dg(seq: str) -> float:
    """dG (kcal/mol) of a fully complementary duplex stretch at 37 C.

    No self-complementarity symmetry correction: dimer/hairpin stretches
    are evaluated as plain duplex segments."""
    dh = ds = 0.0
    seq = seq.upper()
    for i in range(len(seq) - 1):
        h, s = _NN[seq[i : i + 2]]
        dh += h
        ds += s
    for end in (seq[0], seq[-1]):
        h, s = _INIT_AT if end in "AT" else _INIT_GC
        dh += h
        ds += s
    return dh - DG_TEMP_K * ds / 1000.0


def _best_stretch_dg(a: str, b_rc: str) -> float:
    """Most negative dG over contiguous complementary stretches when
    sequence ``a`` is aligned against ``b_rc`` (= revcomp of the partner)
    at every offset.  In this orientation a base pair forms exactly when
    ``a[i] == b_rc[j]`` (Watson-Crick antiparallel pairing).

    The scan accumulates NN (dH, dS) incrementally along each diagonal so
    per-cell work is O(1)."""
    best = 0.0
    na, nb = len(a), len(b_rc)
    for off in range(-(nb - 1), na):
        i = max(0, off)
        j = i - off
        run_len = 0
        run_dh = run_ds = 0.0
        first = prev = ""
        while i < na and j < nb:
            if a[i] == b_rc[j]:
                base = a[i]
                if run_len:
                    h, s = _NN[prev + base]
                    run_dh += h
                    run_ds += s
                else:
                    first = base
                run_len += 1
                prev = base
                if run_len >= MIN_DUPLEX:
                    dh = run_dh
                    ds = run_ds
                    for end in (first, base):
                        h, s = _INIT_AT if end in "AT" else _INIT_GC
                        dh += h
                        ds += s
                    dg = dh - DG_TEMP_K * ds / 1000.0
                    if dg < best:
                        best = dg
            else:
                run_len = 0
                run_dh = run_ds = 0.0
            i += 1
            j += 1
    return best


def self_dimer_dg(seq: str) -> float:
    seq = seq.upper()
    return _best_stretch_dg(seq, revcomp(seq))


def heterodimer_dg(seq_a: str, seq_b: str) -> float:
    return _best_stretch_dg(seq_a.upper(), revcomp(seq_b.upper()))


def hairpin_dg(seq: str) -> float:
    """Most stable intramolecular stem with a loop of >= 3 nt."""
    seq = seq.upper()
    n = len(seq)
    best = 0.0
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    for i in range(n):
        for j in range(i + 2 * MIN_DUPLEX + 3 - 1, n):
            run: list[str] = []
            k = 0
            while i + k < j - k and (j - i + 1 - 2 * k) >= 3 + 2:
                if comp[seq[i + k]] == seq[j - k]:
                    run.append(seq[i + k])
                    k += 1
                    if len(run) >= MIN_DUPLEX and (j - i + 1 - 2 * len(run)) >= 3:
                        best = min(best, _step_dg("".join(run)))
                else:
                    break
    return best


def longest_homopolymer(seq: str) -> int:
    best = run = 1
    for i in range(1, len(seq)):
        run = run + 1 if seq[i] == seq[i - 1] else 1
        best = max(best, run)
    return best


def has_palindromic_3p(seq: str, min_len: int = 4, max_len: int = 8) -> bool:
    seq = seq.upper()
    for k in range(min_len, min(max_len, len(seq)) + 1):
        if is_palindrome(seq[-k:]):
            return True
    return False


def analyze(seq: str) -> StructureReport:
    seq = seq.upper()
    return StructureReport(
        gc_percent=gc_percent(seq),
        tm=tm_celsius(seq),
        annealing_temp=annealing_temp(seq),
        self_dimer_dg=self_dimer_dg(seq),
        hairpin_dg=hairpin_dg(seq),
        homopolymer_run=longest_homopolymer(seq),
        has_palindromic_3p=has_palindromic_3p(seq),
    )

"""Deterministic ranking weights and penalty aggregation.

Every valid primer pair receives a non-negative penalty (lower is better).
The penalty is a weighted sum of:

* Tm distance of each primer from the middle of the allowed Tm window,
* GC distance from 0.50 (ideal 50%),
* pair Tm difference (uneven amplification),
* product length distance from the middle of the requested range,
* self-dimer / hairpin stability beyond a zero-penalty floor,
* heterodimer stability,
* homopolymer and dinucleotide-satellite runs,
* terminal palindromes and strong 3' GC clamps,
* off-target hits weighted by 3'-end match and mismatch count.

All weights are plain constants below; the sort key additionally breaks
ties on product length, forward sequence and reverse sequence so that equal
inputs (and equal seeds) always yield identical ordering.
"""
from __future__ import annotations

# ---------------------------------------------------------------- weights
WEIGHTS: dict[str, float] = {
    # per-primer terms
    "tm_distance": 1.0,          # per deg C from Tm window center
    "gc_distance": 8.0,          # per unit fraction from 0.50
    "self_dimer": 1.0,           # per kcal/mol below floor
    "hairpin": 1.0,              # per kcal/mol below floor
    "homopolymer": 1.5,          # per base beyond 4
    "dinucleotide_run": 1.0,     # per repeat beyond 3
    "terminal_palindrome": 4.0,  # fixed flag
    "gc_clamp": 1.0,             # per extra GC beyond 3 in last 5 bases
    # pair terms
    "tm_delta": 1.5,
    "product_len": 0.02,
    "heterodimer": 1.5,
    # specificity terms (per off-target hit)
    "off_target_base": 3.0,
    "off_target_mm": -2.0,       # fewer mismatches -> heavier (added: mm-2)
    "off_target_tpm": 0.5,       # per consecutive 3' matched base beyond 8
}

# Stability floors: dG >= floor costs zero; below floor penalized.
SELF_DIMER_FLOOR = -3.0
HAIRPIN_FLOOR = -1.0
HETERODIMER_FLOOR = -3.0
HOMOPOLYMER_LIMIT = 4
DINUCLEOTIDE_LIMIT = 3
GC_CLAMP_LIMIT = 3


def _primer_terms(m, params) -> dict[str, float]:
    tm_center = (params.min_tm + params.max_tm) / 2.0
    terms: dict[str, float] = {
        "tm_distance": abs(m.tm - tm_center),
        "gc_distance": abs(m.gc_content - 0.50),
        "self_dimer": max(0.0, SELF_DIMER_FLOOR - m.self_dimer_dg),
        "hairpin": max(0.0, HAIRPIN_FLOOR - m.hairpin_dg),
        "homopolymer": max(0, m.homopolymer_run - HOMOPOLYMER_LIMIT),
        "dinucleotide_run": max(0, m.dinucleotide_run - DINUCLEOTIDE_LIMIT),
        "terminal_palindrome": 1.0 if m.terminal_palindrome else 0.0,
        "gc_clamp": max(0, _gc_clamp(m.sequence) - GC_CLAMP_LIMIT),
    }
    return terms


def _gc_clamp(seq: str, n: int = 5) -> int:
    return sum(1 for ch in seq[-n:] if ch in "GC")


def penalty_for_pair(
    fwd,
    rev,
    product_len: int,
    heterodimer_dg: float,
    off_targets,
    params,
) -> tuple[float, dict[str, float]]:
    """Return (total penalty, per-term breakdown)."""
    detail: dict[str, float] = {}
    total = 0.0

    for label, m in (("fwd", fwd), ("rev", rev)):
        for name, value in _primer_terms(m, params).items():
            w = WEIGHTS[name]
            key = f"{label}_{name}"
            detail[key] = w * value
            total += detail[key]

    tm_delta = abs(fwd.tm - rev.tm)
    product_center = (params.min_product_len + params.max_product_len) / 2.0
    pair_terms = {
        "tm_delta": tm_delta,
        "product_len": abs(product_len - product_center),
        "heterodimer": max(0.0, HETERODIMER_FLOOR - heterodimer_dg),
    }
    for name, value in pair_terms.items():
        detail[name] = WEIGHTS[name] * value
        total += detail[name]

    ot_penalty = 0.0
    for hit in off_targets:
        hit_penalty = (
            WEIGHTS["off_target_base"]
            + WEIGHTS["off_target_mm"] * (hit.mismatches - params.max_mismatches)
            + WEIGHTS["off_target_tpm"]
            * max(0, hit.three_prime_match - params.min_three_prime_match)
        )
        ot_penalty += max(WEIGHTS["off_target_base"], hit_penalty)
    detail["off_target"] = ot_penalty
    total += ot_penalty
    detail["total"] = total
    return total, detail


def sort_key(total: float, product_len: int, fwd_seq: str, rev_seq: str):
    """Deterministic tie-breaker order."""
    return (round(total, 6), product_len, fwd_seq, rev_seq)

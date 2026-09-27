"""Constraint filtering and deterministic weighted ranking.

Penalty model (lower is better; rank_score = -penalty):

===================  ======  ===========================================
term                 weight  definition
===================  ======  ===========================================
Tm deviation         1.0     |Tm - Tm_opt| per primer, Tm_opt = midpoint
                             of the allowed Tm window
Tm difference        2.0     |Tm_fwd - Tm_rev|
GC deviation         0.05    |GC% - 50| per primer
self-dimer           1.0     -min(dG, 0) per primer
heterodimer          1.0     -min(dG, 0)
hairpin              1.0     -min(dG, 0) per primer
homopolymer          2.0     (run - 3) per primer when run >= 4
palindromic 3' end   3.0     flat per primer
off-target hit       5.0     per hit (both primers of the pair)
product size         0.01    |size - midpoint of allowed range|
===================  ======  ===========================================

Ties are broken by a SHA-256 hash of (seed, sequences, positions), so
output is byte-for-byte reproducible for the same input and seed.
"""
from __future__ import annotations

import hashlib
from dataclasses import replace

from .models import DesignParams, Primer, PrimerPair

W_TM = 1.0
W_TM_DIFF = 2.0
W_GC = 0.05
W_DIMER = 1.0
W_HETERODIMER = 1.0
W_HAIRPIN = 1.0
W_HOMOPOLYMER = 2.0
W_PALINDROME_3P = 3.0
W_OFFTARGET = 5.0
W_PRODUCT = 0.01
HOMOPOLYMER_SOFT_LIMIT = 3


def passes_filters(primer: Primer, params: DesignParams) -> bool:
    r = primer.report
    return (
        params.gc_min <= r.gc_percent <= params.gc_max
        and params.tm_min <= r.tm <= params.tm_max
        and r.self_dimer_dg >= params.max_dimer_dg
        and r.hairpin_dg >= params.max_hairpin_dg
    )


def pair_passes(pair: PrimerPair, params: DesignParams) -> bool:
    return pair.heterodimer_dg >= params.max_dimer_dg


def heterodimer_penalty(dg: float) -> float:
    return W_HETERODIMER * (-min(dg, 0.0))


def primer_penalty(report, params: DesignParams) -> float:
    """Per-primer penalty part (cacheable by sequence)."""
    tm_opt = 0.5 * (params.tm_min + params.tm_max)
    pen = W_TM * abs(report.tm - tm_opt)
    pen += W_GC * abs(report.gc_percent - 50.0)
    pen += W_DIMER * (-min(report.self_dimer_dg, 0.0))
    pen += W_HAIRPIN * (-min(report.hairpin_dg, 0.0))
    if report.homopolymer_run > HOMOPOLYMER_SOFT_LIMIT:
        pen += W_HOMOPOLYMER * (report.homopolymer_run - HOMOPOLYMER_SOFT_LIMIT)
    if report.has_palindromic_3p:
        pen += W_PALINDROME_3P
    return pen


def pair_pre_penalty(
    f_pen: float,
    r_pen: float,
    tm_f: float,
    tm_r: float,
    product_size: int,
    params: DesignParams,
) -> float:
    size_opt = 0.5 * (params.product_size_min + params.product_size_max)
    return (
        f_pen
        + r_pen
        + W_TM_DIFF * abs(tm_f - tm_r)
        + W_PRODUCT * abs(product_size - size_opt)
    )


def pre_penalty(pair: PrimerPair, params: DesignParams) -> float:
    """Structure penalty without the heterodimer term (cheap lower bound)."""
    return pair_pre_penalty(
        primer_penalty(pair.forward.report, params),
        primer_penalty(pair.reverse.report, params),
        pair.forward.report.tm,
        pair.reverse.report.tm,
        pair.product_size,
        params,
    )


def structure_penalty(pair: PrimerPair, params: DesignParams) -> float:
    return pre_penalty(pair, params) + heterodimer_penalty(pair.heterodimer_dg)


def offtarget_penalty(pair: PrimerPair) -> float:
    return W_OFFTARGET * len(pair.off_targets)


def _tie_key(pair: PrimerPair, seed: int) -> int:
    h = hashlib.sha256(
        f"{seed}|{pair.forward.sequence}|{pair.forward.start}|"
        f"{pair.reverse.sequence}|{pair.reverse.start}".encode()
    )
    return int.from_bytes(h.digest(), "big")


def rank_pairs(pairs: list[PrimerPair], params: DesignParams) -> list[PrimerPair]:
    """Score, sort (deterministic under ``params.seed``) and truncate."""
    scored = []
    for pair in pairs:
        pen = structure_penalty(pair, params) + offtarget_penalty(pair)
        scored.append(replace(pair, penalty=pen, rank_score=-pen))
    scored.sort(key=lambda p: (-p.rank_score, _tie_key(p, params.seed)))
    return scored[: params.max_pairs]

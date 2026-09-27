"""Candidate primer generation and pair assembly.

Circular templates are handled with modular arithmetic: a primer or an
amplicon may wrap the origin.  Positions are always reported on the
original 0-based top-strand coordinate system.
"""
from __future__ import annotations

from .models import DesignParams, Primer, PrimerPair, TemplateTopology
from .sequtil import normalize_seq, revcomp
from .thermo import analyze, heterodimer_dg


def _template(params: DesignParams) -> str:
    return normalize_seq(params.template)


def _slice_wrap(seq: str, start: int, length: int) -> str:
    """Substring of ``seq`` allowing wrap-around (circular templates)."""
    n = len(seq)
    start %= n
    end = start + length
    if end <= n:
        return seq[start:end]
    return (seq + seq)[start:end]


def generate_primers(params: DesignParams) -> list[Primer]:
    """All position/length primer candidates (both orientations)."""
    tpl = _template(params)
    n_tpl = len(tpl)
    circular = params.topology is TemplateTopology.CIRCULAR
    primers: list[Primer] = []
    report_cache: dict[str, object] = {}

    def cached(seq: str):
        rep = report_cache.get(seq)
        if rep is None:
            rep = report_cache[seq] = analyze(seq)
        return rep

    starts = range(n_tpl) if circular else range(n_tpl - params.primer_len_min + 1)
    for s in starts:
        for plen in range(params.primer_len_min, params.primer_len_max + 1):
            if not circular and s + plen > n_tpl:
                break
            fwd_seq = _slice_wrap(tpl, s, plen) if circular else tpl[s : s + plen]
            primers.append(Primer(fwd_seq, s, plen, True, cached(fwd_seq)))
            rev_seg = fwd_seq
            rev_seq = revcomp(rev_seg)
            primers.append(Primer(rev_seq, s, plen, False, cached(rev_seq)))
    return primers


def covers_target(start: int, product: int, params: DesignParams, n_tpl: int) -> bool:
    """True if the amplicon arc [start, start+product) covers the target.

    On circular templates a target interval with ``target_end <
    target_start`` is interpreted as wrapping the origin:
    [target_start, L) + [0, target_end).
    """
    if params.target_end is None:
        return True
    ts, te = params.target_start, params.target_end
    if not (0 <= ts < n_tpl and 0 < te <= n_tpl) or ts == te:
        raise ValueError("invalid target interval")
    if params.topology is TemplateTopology.CIRCULAR:
        length = te - ts if te > ts else n_tpl - ts + te
        return (ts - start) % n_tpl + length <= product
    if te <= ts:
        raise ValueError("target_end must exceed target_start on a linear template")
    return start <= ts and start + product >= te


def iter_pairs(
    primers: list[Primer],
    params: DesignParams,
    heterodimer: bool = True,
):
    """Yield pairs satisfying product-size and target-coverage constraints.

    With ``heterodimer=False`` the heterodimer term is left as 0.0 so the
    caller can compute it lazily for a promising subset only."""
    tpl = _template(params)
    n_tpl = len(tpl)
    circular = params.topology is TemplateTopology.CIRCULAR
    p_max = min(params.product_size_max, n_tpl)
    p_min = min(params.product_size_min, p_max)

    fwd_by_start: dict[int, list[Primer]] = {}
    rev_by_end: dict[int, list[Primer]] = {}
    for p in primers:
        if p.is_forward:
            fwd_by_start.setdefault(p.start, []).append(p)
        else:
            rev_by_end.setdefault((p.start + p.length) % n_tpl if circular else p.start + p.length, []).append(p)

    hetero_cache: dict[tuple[str, str], float] = {}

    def hetero(a: str, b: str) -> float:
        key = (a, b)
        dg = hetero_cache.get(key)
        if dg is None:
            dg = hetero_cache[key] = heterodimer_dg(a, b)
        return dg

    start_iter = range(n_tpl) if circular else range(n_tpl - p_min + 1)
    for s in start_iter:
        fwds = fwd_by_start.get(s)
        if not fwds:
            continue
        for size in range(p_min, p_max + 1):
            if not circular and s + size > n_tpl:
                break
            if not covers_target(s, size, params, n_tpl):
                continue
            end = (s + size) % n_tpl if circular else s + size
            for f in fwds:
                for r in rev_by_end.get(end, ()):
                    # primers must not overlap within the amplicon
                    gap = (r.start - s) % n_tpl if circular else r.start - s
                    if gap < f.length:
                        continue
                    yield PrimerPair(
                        forward=f,
                        reverse=r,
                        product_size=size,
                        heterodimer_dg=hetero(f.sequence, r.sequence) if heterodimer else 0.0,
                        off_targets=(),
                        penalty=0.0,
                        rank_score=0.0,
                    )


def assemble_pairs(primers: list[Primer], params: DesignParams) -> list[PrimerPair]:
    """Eager variant of :func:`iter_pairs` (heterodimer included)."""
    return list(iter_pairs(primers, params, heterodimer=True))

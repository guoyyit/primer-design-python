"""Primer candidate generation, filtering and deterministic pair ranking."""
from __future__ import annotations

from .models import DesignParams, DesignResult, PrimerMetrics, PrimerPair, Strand
from .scoring import penalty_for_pair, sort_key
from .sequence import gc_content, reverse_complement, slice_template
from .structures import (
    best_hairpin,
    best_heterodimer,
    best_self_dimer,
    longest_dinucleotide_repeat,
    longest_homopolymer,
    terminal_palindrome,
)
from .thermo import nearest_neighbor_tm

_MAX_PAIRS_EVALUATED = 200000


def _metrics(seq: str, strand: Strand, position: int, params: DesignParams) -> PrimerMetrics:
    dg, _, _ = best_self_dimer(seq)
    hg, _ = best_hairpin(seq)
    return PrimerMetrics(
        sequence=seq,
        strand=strand,
        position=position,
        length=len(seq),
        gc_content=gc_content(seq),
        tm=nearest_neighbor_tm(seq, params.tm_concentration_nm),
        self_dimer_dg=dg,
        hairpin_dg=hg,
        homopolymer_run=longest_homopolymer(seq),
        dinucleotide_run=longest_dinucleotide_repeat(seq),
        terminal_palindrome=terminal_palindrome(seq),
    )


def _passes_single(m: PrimerMetrics, params: DesignParams) -> bool:
    return (
        params.min_gc - 1e-12 <= m.gc_content <= params.max_gc + 1e-12
        and params.min_tm - 1e-12 <= m.tm <= params.max_tm + 1e-12
        and m.self_dimer_dg >= params.max_dimer_dg - 1e-12
        and m.hairpin_dg >= params.max_hairpin_dg - 1e-12
    )


def generate_candidates(params: DesignParams) -> tuple[list[PrimerMetrics], list[PrimerMetrics]]:
    """Return filtered (forward, reverse) candidates in deterministic order.

    Coordinates are forward-strand 0-based start positions. A reverse primer
    at ``position`` is the reverse complement of the template window starting
    there (i.e. it extends leftward in forward-strand coordinates).
    """
    template = params.normalized_template()
    n = len(template)
    target_start = params.target_start
    target_end = n if params.target_end is None else params.target_end
    if not (0 <= target_start <= target_end <= n):
        raise ValueError("target interval outside template bounds")
    if params.circular:
        span = n
    else:
        span = target_end - target_start
    if span < params.min_primer_len:
        return [], []

    forwards: list[PrimerMetrics] = []
    reverses: list[PrimerMetrics] = []

    def usable(seq: str) -> bool:
        return set(seq) <= set("ACGT")

    # Forward primers: start inside the target region (circular: starts in
    # the forward-coordinate interval, bases may wrap past n).
    f_positions = range(target_start, target_end)
    for pos in f_positions:
        if not params.circular and pos + params.min_primer_len > n:
            break
        for length in range(params.min_primer_len, params.max_primer_len + 1):
            if not params.circular and pos + length > n:
                break
            seq = slice_template(template, pos, length, params.circular)
            if not usable(seq):
                continue
            m = _metrics(seq, Strand.FORWARD, pos, params)
            if _passes_single(m, params):
                forwards.append(m)

    # Reverse primers: window [pos, pos+length) starts inside the target
    # region; circular windows may wrap across the origin.
    r_positions = range(target_start, target_end)
    for end in r_positions:
        for length in range(params.min_primer_len, params.max_primer_len + 1):
            pos = end - length + 1
            if not params.circular:
                if pos < target_start:
                    continue
            else:
                pos %= n
            win = slice_template(template, pos, length, params.circular)
            if not usable(win):
                continue
            seq = reverse_complement(win)
            m = _metrics(seq, Strand.REVERSE, pos, params)
            if _passes_single(m, params):
                reverses.append(m)

    forwards.sort(key=lambda m: (m.position, m.length))
    reverses.sort(key=lambda m: (m.position, m.length))
    return forwards, reverses


def evaluate_pair_metrics(fwd: PrimerMetrics, rev: PrimerMetrics, product_len: int) -> dict:
    dg, _, _ = best_heterodimer(fwd.sequence, rev.sequence)
    return {"heterodimer_dg": dg}


def design_primers(
    params: DesignParams,
    evaluator=None,
    progress_cb=None,
    cancel_cb=None,
) -> DesignResult:
    """Run the full deterministic design pipeline.

    Progress (progress_cb receives a float in [0, 1]):
        0.00-0.20  candidate generation/filtering
        0.20-0.60  pair enumeration & structural scoring
        0.60-1.00  specificity evaluation of the surviving shortlist
    """

    def report(value: float) -> None:
        if progress_cb is not None:
            progress_cb(min(1.0, max(0.0, value)))

    def cancelled() -> bool:
        return bool(cancel_cb and cancel_cb())

    import time as _time

    report(0.0)
    forwards, reverses = generate_candidates(params)
    n = len(params.normalized_template())
    report(0.2)

    scored: list[tuple] = []
    evaluated = 0
    total_pairs = max(1, len(forwards) * len(reverses))
    pair_step = max(1, total_pairs // 40)

    for fi, fwd in enumerate(forwards):
        if cancelled():
            raise CancelledError("design cancelled during pair enumeration")
        for ri, rev in enumerate(reverses):
            evaluated += 1
            if evaluated > _MAX_PAIRS_EVALUATED:
                break
            # Amplicon coordinates: fwd points right, rev points left.
            if params.circular:
                # product wraps: distance from fwd start to rev window end.
                rev_end = (rev.position + rev.length) % n
                product_len = (rev_end - fwd.position) % n
                if product_len == 0:
                    product_len = n
            else:
                if rev.position < fwd.position:
                    continue
                product_len = rev.position + rev.length - fwd.position
            if not (params.min_product_len <= product_len <= params.max_product_len):
                continue
            if product_len < fwd.length + rev.length:
                continue
            extra = evaluate_pair_metrics(fwd, rev, product_len)
            dg = extra["heterodimer_dg"]
            if dg < params.max_dimer_dg - 1e-12:
                continue
            total, detail = penalty_for_pair(
                fwd, rev, product_len, dg, (), params
            )
            scored.append((total, product_len, fwd, rev, dg, detail))
        if evaluated > _MAX_PAIRS_EVALUATED:
            break
        if fi % pair_step == 0:
            report(0.2 + 0.4 * (fi / max(1, len(forwards))))

    # structural-only ordering selects the specificity shortlist
    scored.sort(key=lambda row: sort_key(row[0], row[1], row[2].sequence, row[3].sequence))
    shortlist = scored[: max(1, params.max_candidates * 3)]
    report(0.6)

    pairs: list[PrimerPair] = []
    shortlist_total = max(1, len(shortlist))
    for idx, (pre_total, product_len, fwd, rev, dg, _detail) in enumerate(shortlist):
        if cancelled():
            raise CancelledError("design cancelled during specificity evaluation")
        hits = ()
        if evaluator is not None:
            def _sub_progress(frac, _base=0.6 + 0.4 * idx / shortlist_total,
                              _span=0.4 / shortlist_total):
                report(_base + _span * frac)

            hits = evaluator.evaluate_pair(
                fwd.sequence, rev.sequence, progress_cb=_sub_progress,
                cancel_cb=cancel_cb,
            )
        total, detail = penalty_for_pair(
            fwd, rev, product_len, dg, hits, params
        )
        pairs.append(
            PrimerPair(
                pair_id=f"P{idx + 1:04d}",
                forward=fwd,
                reverse=rev,
                product_len=product_len,
                heterodimer_dg=dg,
                score=total,
                rank=0,
                off_targets=hits,
                penalty_detail=detail,
            )
        )

    pairs.sort(
        key=lambda p: sort_key(p.score, p.product_len, p.forward.sequence,
                               p.reverse.sequence)
    )
    final = pairs[: params.max_candidates]
    ranked = [
        PrimerPair(
            pair_id=f"P{i + 1:04d}",
            forward=p.forward,
            reverse=p.reverse,
            product_len=p.product_len,
            heterodimer_dg=p.heterodimer_dg,
            score=p.score,
            rank=i + 1,
            off_targets=p.off_targets,
            penalty_detail=p.penalty_detail,
        )
        for i, p in enumerate(final)
    ]
    report(1.0)
    return DesignResult(
        pairs=tuple(ranked),
        evaluated_pairs=evaluated,
        generated_at=_time.time(),
    )


class CancelledError(Exception):
    """Raised when a design run observes a cancellation request."""

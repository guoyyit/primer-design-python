"""Candidate generation: circular wrap-around, target coverage, extremes."""
import pytest

from primer_design.models import DesignParams, TemplateTopology
from primer_design import search, scoring
from primer_design.sequtil import revcomp

TPL = "ACGT" * 100  # 400 bp circular molecule


def test_linear_template_primer_in_bounds():
    params = DesignParams(template="A" * 100 + "CG" * 50, primer_len_min=18)
    primers = search.generate_primers(params)
    for p in primers:
        if p.is_forward:
            assert 0 <= p.start
            assert p.start + p.length <= 200


def test_circular_primer_wraps_origin():
    params = DesignParams(
        template=TPL, topology=TemplateTopology.CIRCULAR,
        primer_len_min=18, primer_len_max=22,
    )
    primers = search.generate_primers(params)
    # every origin position appears as a forward start
    fwd_starts = {p.start for p in primers if p.is_forward}
    assert fwd_starts == set(range(400))
    # a start near the origin yields a wrap-around sequence
    near = next(p for p in primers if p.is_forward and p.start == 395 and p.length == 20)
    assert near.sequence == (TPL * 2)[395:415]


def test_circular_amplicon_covers_target_near_origin():
    # target wraps the origin
    params = DesignParams(
        template=TPL, topology=TemplateTopology.CIRCULAR,
        target_start=390, target_end=20,
        product_size_min=60, product_size_max=120,
    )
    assert search.covers_target(380, 60, params, 400) is True
    assert search.covers_target(0, 50, params, 400) is False  # misses 390-399
    assert search.covers_target(100, 100, params, 400) is False
    primers = search.generate_primers(params)
    pairs = search.assemble_pairs(primers, params)
    assert pairs
    for pair in pairs:
        # target [390,400) U [0,20) must be inside the arc
        s = pair.forward.start
        assert search.covers_target(s, pair.product_size, params, 400)


def test_invalid_target_interval_rejected():
    params = DesignParams(template=TPL, target_start=50, target_end=10)
    with pytest.raises(ValueError):
        search.covers_target(0, 80, params, len(TPL))


def test_extreme_gc_all_filtered_or_accepted():
    all_at = "AT" * 200
    params = DesignParams(
        template=all_at, gc_min=30.0, gc_max=70.0,
        product_size_min=40, product_size_max=80,
    )
    primers = search.generate_primers(params)
    good = [p for p in primers if scoring.passes_filters(p, params)]
    assert good == []  # 0% GC cannot pass a 30% floor
    all_gc = "GC" * 200
    params_gc = DesignParams(
        template=all_gc, gc_min=30.0, gc_max=70.0,
        product_size_min=40, product_size_max=80,
    )
    good_gc = [p for p in search.generate_primers(params_gc) if scoring.passes_filters(p, params_gc)]
    assert good_gc == []  # 100% GC above ceiling
    # relaxed window keeps candidates
    relaxed = DesignParams(
        template=all_gc, gc_min=0.0, gc_max=100.0, tm_min=0.0, tm_max=100.0,
        max_dimer_dg=-100.0, max_hairpin_dg=-100.0,
        product_size_min=40, product_size_max=80,
    )
    assert any(scoring.passes_filters(p, relaxed) for p in search.generate_primers(relaxed))

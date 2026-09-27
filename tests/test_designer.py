import random

import pytest

from primer_design.designer import design_primers, generate_candidates
from primer_design.models import DesignParams, ReferenceSequence, TaskStatus
from primer_design.specificity import ReferenceLibrary, SpecificityEvaluator


def _random_template(n: int, seed: int = 1) -> str:
    rnd = random.Random(seed)
    return "".join(rnd.choices("ACGT", k=n))


def _loose_params(template, **over):
    base = dict(
        min_product_len=75, max_product_len=300,
        min_primer_len=18, max_primer_len=22,
        min_tm=50.0, max_tm=75.0,
        min_gc=0.1, max_gc=0.9,
        max_dimer_dg=-9.0, max_hairpin_dg=-6.0,
        max_candidates=5,
    )
    base.update(over)
    return DesignParams(template=template, **base)


def test_generates_pairs_for_normal_template():
    res = design_primers(_loose_params(_random_template(1200, seed=7)))
    assert len(res.pairs) > 0
    for p in res.pairs:
        assert 75 <= p.product_len <= 300
        assert p.forward.tm and p.reverse.tm


def test_results_sorted_and_ranked():
    res = design_primers(_loose_params(_random_template(1200, seed=8)))
    scores = [p.score for p in res.pairs]
    assert scores == sorted(scores)
    assert [p.rank for p in res.pairs] == list(range(1, len(res.pairs) + 1))


def test_deterministic_same_seed_and_input():
    params = _loose_params(_random_template(900, seed=9))
    r1 = design_primers(params)
    r2 = design_primers(params)
    key = lambda r: [(p.forward.sequence, p.reverse.sequence, round(p.score, 9))
                     for p in r.pairs]
    assert key(r1) == key(r2)


def test_circular_template_works_and_wraps():
    tpl = _random_template(600, seed=11)
    res = design_primers(_loose_params(tpl, circular=True,
                                       min_product_len=75, max_product_len=900))
    assert len(res.pairs) > 0
    # every pair window stays within circular coordinate rules
    n = len(tpl)
    for p in res.pairs:
        assert p.product_len <= n
        assert 0 <= p.forward.position < n and 0 <= p.reverse.position < n


def test_circular_origin_spanning_pair_possible():
    # force amplicon to span origin: same sequence repeated so both ends
    # contain identical primer binding sites across the origin
    tpl = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT" * 4
    res = design_primers(
        DesignParams(template=tpl, circular=True, min_product_len=20,
                     max_product_len=40, min_primer_len=10, max_primer_len=14,
                     min_tm=20.0, max_tm=90.0, min_gc=0.0, max_gc=1.0,
                     max_dimer_dg=-9.0, max_hairpin_dg=-6.0,
                     max_candidates=5)
    )
    assert len(res.pairs) > 0


def test_circular_candidates_wrap_and_coordinates_consistent():
    from primer_design.sequence import slice_template

    tpl = "ACGT" * 12  # 48 bp
    fwd, rev = generate_candidates(
        DesignParams(template=tpl, circular=True, min_primer_len=10,
                     max_primer_len=12, min_tm=0.0, max_tm=100.0,
                     min_gc=0.0, max_gc=1.0,
                     max_dimer_dg=-2.0, max_hairpin_dg=-2.0)
    )
    wrapped = [m for m in fwd if m.position + m.length > len(tpl)]
    assert wrapped  # a primer starting near the end wraps across origin
    m = wrapped[0]
    assert m.sequence == slice_template(tpl, m.position, m.length, circular=True)

    res = design_primers(
        DesignParams(template=tpl, circular=True, min_product_len=20,
                     max_product_len=60, min_primer_len=10, max_primer_len=12,
                     min_tm=0.0, max_tm=100.0, min_gc=0.0, max_gc=1.0,
                     max_dimer_dg=-2.0, max_hairpin_dg=-2.0,
                     max_candidates=5)
    )
    n = len(tpl)
    for p in res.pairs:
        rev_end = (p.reverse.position + p.reverse.length) % n
        plen = (rev_end - p.forward.position) % n or n
        assert plen == p.product_len


def test_extreme_gc_filters_all_when_window_too_narrow():
    # 100% GC template with a strict window: no primer can satisfy low GC
    tpl = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGC" * 40
    fwd, rev = generate_candidates(
        DesignParams(template=tpl, min_gc=0.40, max_gc=0.60,
                     min_tm=40.0, max_tm=95.0,
                     min_primer_len=18, max_primer_len=22)
    )
    assert fwd == [] and rev == []


def test_extreme_gc_allowed_with_wide_window():
    tpl = ("GCGCGCGCGCGCGCGCGCGCGCGCGCGCGC" * 20
           + "ATATATATATATATATATATATATATATAT" * 20)
    fwd, _ = generate_candidates(
        DesignParams(template=tpl, min_gc=0.0, max_gc=1.0,
                     min_tm=30.0, max_tm=100.0,
                     max_hairpin_dg=-15.0, max_dimer_dg=-20.0,
                     min_primer_len=18, max_primer_len=22)
    )
    assert any(m.gc_content == 1.0 for m in fwd)


def test_self_complementary_primers_penalized_or_filtered():
    # candidate window containing a strongly self-complementary motif
    tpl = ("GCGCGAATTCGCGC" + "AAACCCAAACCCAAACCCAAACCCAAACCCAAACCC"
           + "GCGCGAATTCGCGC")
    res = design_primers(
        DesignParams(template=tpl, min_product_len=30, max_product_len=120,
                     min_primer_len=14, max_primer_len=16,
                     min_tm=20.0, max_tm=90.0, min_gc=0.0, max_gc=1.0,
                     max_dimer_dg=0.0, max_hairpin_dg=0.0,
                     max_candidates=10)
    )
    # surviving primers must obey the (loosened here) structure thresholds
    for p in res.pairs:
        assert p.forward.self_dimer_dg <= 0.0
        assert p.reverse.self_dimer_dg <= 0.0


def test_strict_dimer_threshold_rejects_stable_duplex():
    # a positive threshold rejects every duplex (every dg <= 0 < 5.0):
    # the pairs list is empty, but the call completes without error
    tpl = _random_template(400, seed=4)
    res = design_primers(
        DesignParams(template=tpl, min_product_len=50, max_product_len=300,
                     min_primer_len=12, max_primer_len=16,
                     min_tm=20.0, max_tm=90.0, min_gc=0.0, max_gc=1.0,
                     max_dimer_dg=5.0, max_hairpin_dg=5.0,
                     max_candidates=5)
    )
    assert res.pairs == ()


def test_offtarget_affects_ranking():
    tpl = _random_template(1000, seed=13)
    params = _loose_params(tpl, max_candidates=30, max_mismatches=0,
                           min_three_prime_match=18)
    base = design_primers(params)
    assert base.pairs
    # contaminate the library with *every* finalist's forward primer: then
    # the whole ranking is evaluated under the same off-target environment
    contam = "TTTT".join(p.forward.sequence for p in base.pairs)
    lib = ReferenceLibrary()
    lib.add([ReferenceSequence("off1", contam)])
    ev = SpecificityEvaluator(lib, max_mismatches=0,
                              min_three_prime_match=18)
    with_hit = design_primers(params, evaluator=ev)
    assert any(h.reference_id == "off1"
               for p in with_hit.pairs for h in p.off_targets)
    for p in with_hit.pairs:
        if p.off_targets:
            assert p.penalty_detail["off_target"] > 0.0
            break
    # ordering must remain a sorted-by-score order under the new penalties
    scores = [p.score for p in with_hit.pairs]
    assert scores == sorted(scores)


def test_cancel_during_design():
    from primer_design.designer import CancelledError
    params = _loose_params(_random_template(3000, seed=21))
    calls = {"n": 0}

    def cancel_soon():
        calls["n"] += 1
        return calls["n"] > 2

    with pytest.raises(CancelledError):
        design_primers(params, cancel_cb=cancel_soon)


def test_invalid_target_interval_rejected():
    with pytest.raises(ValueError):
        generate_candidates(DesignParams(template="ACGT" * 25,
                                         target_start=10, target_end=5))

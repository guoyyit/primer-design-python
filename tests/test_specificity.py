"""Off-target scanning: mismatch limits, 3' rule, both strands, chunking."""
from primer_design.specificity import evaluate_primer
from primer_design.sequtil import revcomp


def test_mismatch_boundary():
    primer = "ACGTACGT"
    # exactly max_mismatches differences, strong 3' end -> hit
    ref = "AAA" + "ACGT" + "ACGT" + "AAA"
    hits = evaluate_primer(primer, [("r", ref)], max_mismatches=0, max_3p_consecutive=3)
    assert len(hits) == 1 and hits[0].mismatches == 0
    # one extra mismatch beyond the limit -> no hit
    ref2 = "AAA" + "TCGTACGT" + "AAA"  # 1 mismatch at 5' end
    assert evaluate_primer(primer, [("r", ref2)], 0, 3) == []
    assert len(evaluate_primer(primer, [("r", ref2)], 1, 3)) == 1


def test_three_prime_consecutive_boundary():
    primer = "ACGTACGT"
    # 3' run of exactly 4 with 4 mismatches at the 5' half
    # (GGCC mismatches ACGT at every position, T vs C caps the run at 4)
    ref = "AAA" + "GGCC" + "ACGT" + "AAA"
    # limit 4 -> safe (strictly greater required)
    assert evaluate_primer(primer, [("r", ref)], 4, 4) == []
    # limit 3 -> reported
    hits = evaluate_primer(primer, [("r", ref)], 4, 3)
    assert len(hits) == 1
    assert hits[0].three_prime_consecutive == 4
    assert hits[0].mismatches == 4


def test_minus_strand_scanned_and_position_mapped():
    primer = "AACCGGTA"  # non-palindromic
    ref = "GGGG" + revcomp(primer) + "CCCC"
    hits = evaluate_primer(primer, [("r", ref)], 0, 3)
    assert len(hits) == 1
    hit = hits[0]
    assert hit.strand == "-"
    # binding site on forward coordinates starts at 4
    assert hit.position == 4


def test_palindromic_site_deduplicated():
    primer = "GAATTC"  # EcoRI palindrome: matches both strands at same locus
    ref = "TTT" + "GAATTC" + "TTT"
    hits = evaluate_primer(primer, [("r", ref)], 0, 2)
    assert len(hits) == 1


def test_chunked_scan_matches_unchunked():
    import random
    rng = random.Random(5)
    primer = "".join(rng.choice("ACGT") for _ in range(15))
    ref = "".join(rng.choice("ACGT") for _ in range(5000))
    big = evaluate_primer(primer, [("r", ref)], 3, 4, chunk_size=1_000_000)
    small = evaluate_primer(primer, [("r", ref)], 3, 4, chunk_size=64)
    key = lambda h: (h.reference_id, h.position, h.strand, h.mismatches,
                     h.three_prime_consecutive)
    assert sorted(map(key, big)) == sorted(map(key, small))


def test_multiple_references():
    primer = "ACGTACGT"
    refs = [("a", "ACGTACGT"), ("b", "TTTTTTTT"), ("c", "GGACGTACGTGG")]
    hits = evaluate_primer(primer, refs, 0, 3)
    assert {h.reference_id for h in hits} == {"a", "c"}

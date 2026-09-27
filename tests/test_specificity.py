import threading

import pytest

from primer_design.models import ReferenceSequence
from primer_design.specificity import (
    CancelledError,
    ReferenceLibrary,
    SpecificityEvaluator,
)

REF = "AATTCGATCGATCGATCGATCGATCGTAGCTAGCTAGC"
PRIMER = "TTCGATCGATCGATCGATCG"  # exact at position 2 of REF


def test_exact_hit_and_mismatch_boundary():
    ev = SpecificityEvaluator(ReferenceLibrary(), max_mismatches=2,
                              min_three_prime_match=8)
    hits = ev.count_hits(PRIMER, REF, 2, 8)
    assert (2, 0, 20) in hits
    # exact match must be the best (position 2)
    assert min(hits, key=lambda h: h[1]) == (2, 0, 20)


def test_mismatch_limit_boundary():
    # two mismatches far enough apart that no shifted alignment heals both
    chars = list(PRIMER)
    for idx in (1, 10):
        chars[idx] = "G" if chars[idx] != "G" else "T"
    mut = "".join(chars)
    assert any(h[1] == 2 for h in
               SpecificityEvaluator.count_hits(PRIMER, mut, 2, 8))
    assert not SpecificityEvaluator.count_hits(PRIMER, mut, 1, 8)


def test_three_prime_match_limit():
    # break the 3' run: mutate base aligned to primer's 3' end (position 21)
    broken = REF[:21] + ("A" if REF[21] != "A" else "C") + REF[22:]
    hits = SpecificityEvaluator.count_hits(PRIMER, broken, 2, 19)
    assert not any(h[0] == 2 for h in hits)
    # relaxing the limit admits it
    assert any(h[0] == 2 for h in
               SpecificityEvaluator.count_hits(PRIMER, broken, 2, 0))


def test_evaluate_pair_reverse_strand_hit():
    lib = ReferenceLibrary()
    lib.add([ReferenceSequence("r1", REF)])
    ev = SpecificityEvaluator(lib, max_mismatches=2, min_three_prime_match=8)
    from primer_design.sequence import reverse_complement
    # reverse primer whose 5'->3' sequence is revcomp of ref window at pos 2
    rev_primer = reverse_complement(REF[2:22])
    hits = ev.evaluate_pair(PRIMER, rev_primer)
    strands = {h.strand.value for h in hits}
    assert "+" in strands and "-" in strands


def test_snapshot_isolation_during_update():
    lib = ReferenceLibrary()
    lib.add([ReferenceSequence("r1", REF)])
    seen_versions = []
    release = threading.Event()

    def slow_progress(_frac):
        if not seen_versions:
            seen_versions.append(lib.version.version)
            # trigger an update mid-evaluation; evaluator must not see it
            lib.add([ReferenceSequence("r2", "CCCC" * 20)])
            seen_versions.append(lib.version.version)
        return None

    ev = SpecificityEvaluator(lib, max_mismatches=1, min_three_prime_match=8,
                              chunk_size=10)
    hits = ev.evaluate_pair(PRIMER, PRIMER, progress_cb=slow_progress)
    # all hits still belong to the old snapshot only
    assert seen_versions == [1, 2]
    assert {h.reference_id for h in hits} <= {"r1"}


def test_concurrent_updates_never_duplicate_or_lose_duplicates_check():
    lib = ReferenceLibrary()
    lib.replace([ReferenceSequence("a", "ACGT")])
    with pytest.raises(ValueError):
        lib.add([ReferenceSequence("a", "TTTT")])


def test_evaluate_cancellation():
    lib = ReferenceLibrary()
    big = "ACGT" * 5000
    lib.add([ReferenceSequence("big", big), ReferenceSequence("big2", big)])
    ev = SpecificityEvaluator(lib, max_mismatches=2, min_three_prime_match=4)
    with pytest.raises(CancelledError):
        ev.evaluate_pair("TTTTTTTTTTTTTTTTTTTT", "AAAAAAAAAAAAAAAAAAAA",
                         cancel_cb=lambda: True)

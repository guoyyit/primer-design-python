from primer_design.thermo import gc_clamped, nearest_neighbor_tm


def test_tm_ordering_gc_vs_at():
    at = nearest_neighbor_tm("AAAATTTTAAAATTTTAAAA")
    mixed = nearest_neighbor_tm("GCGCAGTCGATCGATCGATC")
    gc = nearest_neighbor_tm("GCGCGCGCGCGCGCGCGCGC")
    assert at < mixed < gc


def test_tm_concentration_effect():
    seq = "GCGCAGTCGATCGATCGATC"
    assert nearest_neighbor_tm(seq, 100) < nearest_neighbor_tm(seq, 500)


def test_short_sequence_tm_zero():
    assert nearest_neighbor_tm("A") == 0.0


def test_gc_clamp():
    assert gc_clamped("AAAAAAGG") == 2
    assert gc_clamped("AAAAAGCGC") == 4

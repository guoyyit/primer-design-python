import pytest

from primer_design.sequence import (
    gc_content,
    normalize_sequence,
    reverse_complement,
    slice_template,
)


def test_normalize_strips_and_uppercases():
    assert normalize_sequence(" a c\ngt\n") == "ACGT"


def test_normalize_rejects_garbage():
    with pytest.raises(ValueError):
        normalize_sequence("ZZZ")
    with pytest.raises(ValueError):
        normalize_sequence("")


def test_reverse_complement():
    assert reverse_complement("AACGT") == "ACGTT"
    assert reverse_complement("N") == "N"


def test_gc_content():
    assert gc_content("") == 0.0
    assert gc_content("ACGT") == 0.5
    assert gc_content("GGCC") == 1.0


def test_slice_linear_truncates():
    assert slice_template("ABCDEF", 4, 4, circular=False) == "EF"
    assert slice_template("ABCDEF", 7, 2, circular=False) == ""


def test_slice_circular_wraps():
    assert slice_template("ABCDEF", 4, 4, circular=True) == "EFAB"
    assert slice_template("ABCDEF", 0, 8, circular=True) == "ABCDEFAB"

from primer_design.structures import (
    align_energy,
    best_hairpin,
    best_heterodimer,
    best_self_dimer,
    longest_dinucleotide_repeat,
    longest_homopolymer,
    terminal_palindrome,
)


def test_align_energy_complement_only():
    assert align_energy("AAAA", "TTTT") < 0
    assert align_energy("AAAA", "AAAA") == 0.0
    # isolated single pairs next to a gap cannot stack
    assert align_energy("A-A", "T-T") == 0.0


def test_self_dimer_complementary_ends():
    dg, top, bot = best_self_dimer("GCGCGAATTCGCGC")
    assert dg < 0
    assert top and bot


def test_heterodimer_stable_for_complement():
    dg, _, _ = best_heterodimer("AAAAAAAAAAAA", "TTTTTTTTTTTT")
    assert dg < -10


def test_heterodimer_zero_for_nonpairing():
    dg, _, _ = best_heterodimer("AAAAAAAAAAAA", "AAAAAAAAAAAA")
    # AA/TT stacking impossible -> no stabilizing alignment
    assert dg == 0.0


def test_hairpin_detected_and_scored():
    dg, struct = best_hairpin("GGGCCCTTTGGGCCC")
    assert dg < 0
    assert "(" in struct


def test_no_hairpin_for_polya():
    dg, struct = best_hairpin("AAAAAAAAAAAAAAAA")
    assert dg == 0.0
    assert struct == ""


def test_runs():
    assert longest_homopolymer("AAAGGGGGGTT") == 6
    assert longest_dinucleotide_repeat("ATATATATGG") == 4
    assert longest_dinucleotide_repeat("GGCATATATTT") == 3
    assert longest_dinucleotide_repeat("GGCATCATTT") == 1


def test_terminal_palindrome():
    assert terminal_palindrome("GCGCAAAAGCGC") is True
    assert terminal_palindrome("AAACAAAAGCGC") is False
    assert terminal_palindrome("ACGT") is False  # too short

"""Thermodynamics: Tm/GC scale, dimers, hairpins, repeats/palindromes."""
import math

from primer_design import thermo


def test_gc_percent():
    assert thermo.gc_percent("ACGT") == 50.0
    assert thermo.gc_percent("AAAA") == 0.0
    assert thermo.gc_percent("GGCC") == 100.0


def test_tm_monotonic_in_gc_and_length():
    at = "AT" * 10
    gc = "GC" * 10
    assert thermo.tm_celsius(gc) > thermo.tm_celsius(at)
    assert thermo.tm_celsius("ACGT" * 6) > thermo.tm_celsius("ACGT" * 5)
    assert thermo.annealing_temp("ACGT" * 5) == thermo.tm_celsius("ACGT" * 5) - 5.0


def test_perfectly_complementary_heterodimer_is_very_stable():
    a = "GCGATCGATCGATCGATCG"  # 20 nt, GC-bearing
    from primer_design.sequtil import revcomp
    dg = thermo.heterodimer_dg(a, revcomp(a))
    assert dg < -20.0
    # 17 nt AT/GT-rich perfect duplex is clearly negative too
    b = "AGTCAGTCAGTCAGTC"
    assert thermo.heterodimer_dg(b, revcomp(b)) < -10.0
    # non-complementary control is ~0
    assert thermo.heterodimer_dg("AAAAA", "CCCCC") == 0.0


def test_self_dimer_palindrome():
    assert thermo.self_dimer_dg("GAATTC" * 3) < thermo.self_dimer_dg(
        "AAACCCGGGTTTAAA"
    )


def test_hairpin_stem_loop_detected():
    # GGGG + AAA loop + CCCC complementary -> stable hairpin
    stable = thermo.hairpin_dg("GGGGAAAACCCC")
    assert stable < 0.0
    # shuffled control with no long inverted repeat is much weaker
    control = thermo.hairpin_dg("AAAGGGAAAACCC")
    assert stable < control


def test_homopolymer_and_palindromic_end():
    assert thermo.longest_homopolymer("AAATTTCCCC") == 4
    assert thermo.longest_homopolymer("ACGTACGT") == 1
    # GAATTC is a restriction-site palindrome at the 3' end
    assert thermo.has_palindromic_3p("TTTTTGAATTC") is True
    assert thermo.has_palindromic_3p("AGTAGTAGTAGT") is False


def test_analyze_report_complete():
    report = thermo.analyze("ACGTACGTACGTACGT")
    assert math.isfinite(report.tm)
    assert report.annealing_temp < report.tm
    assert math.isfinite(report.self_dimer_dg)
    assert math.isfinite(report.hairpin_dg)

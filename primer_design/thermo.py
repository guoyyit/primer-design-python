"""Melting-temperature calculations.

Uses the unified nearest-neighbor parameters of SantaLucia (1998) with
initiation and terminal AT penalties. No salt correction is applied (the
service treats the reaction as salt-free); primer concentration is the
total strand concentration Ct (non-self-complementary approximation).
"""
from __future__ import annotations
import math

# deltaH kcal/mol, deltaS cal/(mol*K)
_NN: dict[str, tuple[float, float]] = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AC": (-8.4, -22.4), "GT": (-8.4, -22.4),
    "AG": (-7.8, -21.0), "CT": (-7.8, -21.0),
    "AT": (-7.2, -20.4),
    "CA": (-8.5, -22.8), "TG": (-7.8, -21.0),
    "CC": (-8.0, -19.9), "GG": (-8.0, -19.9),
    "CG": (-10.6, -27.2),
    "GA": (-8.2, -22.4), "TC": (-8.2, -22.4),
    "GC": (-9.8, -24.4),
    "TA": (-7.2, -21.3),
}
_INIT_DH = 0.2
_INIT_DS = -5.7
_TERM_AT_DH = 2.2
_TERM_AT_DS = 6.9
_R = 1.987  # cal/(mol*K)


def nearest_neighbor_tm(seq: str, conc_nm: float = 250.0) -> float:
    """Return duplex melting temperature in degrees Celsius."""
    seq = seq.upper()
    if len(seq) < 2:
        return 0.0
    dh = _INIT_DH
    ds = _INIT_DS
    for i in range(len(seq) - 1):
        step = seq[i:i + 2]
        if step not in _NN:
            # Ambiguous base: neutral contribution.
            continue
        h, s = _NN[step]
        dh += h
        ds += s
    if seq[0] in "AT":
        dh += _TERM_AT_DH
        ds += _TERM_AT_DS
    if seq[-1] in "AT":
        dh += _TERM_AT_DH
        ds += _TERM_AT_DS
    ct = max(conc_nm, 1e-9) * 1e-9  # nM -> M
    tm_k = (1000.0 * dh) / (ds + _R * math.log(ct / 4.0))
    return tm_k - 273.15


def gc_clamped(seq: str, n: int = 5) -> int:
    """Number of G/C bases among the last ``n`` bases at the 3' end."""
    tail = seq[-n:] if n > 0 else ""
    return sum(1 for ch in tail if ch in "GC")

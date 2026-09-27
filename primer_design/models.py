"""Core data types for the primer design service.

Only type definitions live here; no behaviour.
"""
from __future__ import annotations

import enum
from dataclasses import dataclass, field
from typing import Optional


class TemplateTopology(str, enum.Enum):
    LINEAR = "linear"
    CIRCULAR = "circular"


class TaskStatus(str, enum.Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    CANCELLED = "cancelled"
    FAILED = "failed"


@dataclass(frozen=True)
class DesignParams:
    """Parameters controlling primer candidate generation and filtering."""

    template: str
    topology: TemplateTopology = TemplateTopology.LINEAR
    # Target interval [start, end) on the template (0-based). The amplicon
    # must overlap / cover this interval.
    target_start: int = 0
    target_end: Optional[int] = None
    product_size_min: int = 80
    product_size_max: int = 300
    primer_len_min: int = 18
    primer_len_max: int = 25
    gc_min: float = 30.0          # percent
    gc_max: float = 70.0          # percent
    tm_min: float = 45.0          # Celsius (model-relative, see docs)
    tm_max: float = 65.0          # Celsius
    max_dimer_dg: float = -5.0    # kcal/mol; reject if dG < this (more stable)
    max_hairpin_dg: float = -3.0  # kcal/mol; reject if dG < this
    max_pairs: int = 10
    seed: int = 0                 # for deterministic tie-breaking
    # Specificity evaluation
    reference_library: Optional[str] = None  # library name; None => skip
    max_mismatches: int = 3
    max_3p_consecutive: int = 4   # max allowed 3'-end consecutive matches


@dataclass(frozen=True)
class StructureReport:
    """Thermodynamic / structural annotation of one primer."""

    gc_percent: float
    tm: float
    annealing_temp: float
    self_dimer_dg: float
    hairpin_dg: float
    homopolymer_run: int          # longest mono-nucleotide run
    has_palindromic_3p: bool      # 3' end palindrome detected


@dataclass(frozen=True)
class Primer:
    sequence: str
    start: int                    # 0-based position on template (forward coord)
    length: int
    is_forward: bool
    report: StructureReport


@dataclass(frozen=True)
class OffTargetHit:
    reference_id: str
    position: int                 # 0-based on the forward reference strand
    strand: str                   # "+" or "-" (strand the primer matches)
    mismatches: int
    three_prime_consecutive: int


@dataclass(frozen=True)
class PrimerPair:
    forward: Primer
    reverse: Primer
    product_size: int
    heterodimer_dg: float
    off_targets: tuple[OffTargetHit, ...]
    penalty: float                # composite penalty (lower is better)
    rank_score: float             # final score (higher is better)


@dataclass
class TaskRecord:
    """Persisted state of a design task."""

    task_id: str
    status: TaskStatus
    params: Optional[DesignParams] = None
    progress_done: int = 0
    progress_total: int = 0
    result: Optional[list[PrimerPair]] = None
    error: Optional[str] = None

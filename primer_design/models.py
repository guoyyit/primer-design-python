"""Core data types for the primer design service."""
from __future__ import annotations

import time
import uuid
from dataclasses import dataclass, field
from enum import Enum


class Strand(str, Enum):
    FORWARD = "+"
    REVERSE = "-"


class TaskStatus(str, Enum):
    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


@dataclass(frozen=True)
class DesignParams:
    """Parameters controlling one primer-design task."""

    template: str
    target_start: int = 0
    target_end: int | None = None  # exclusive; None -> whole template
    circular: bool = False
    min_product_len: int = 75
    max_product_len: int = 300
    min_primer_len: int = 18
    max_primer_len: int = 27
    min_gc: float = 0.40
    max_gc: float = 0.60
    min_tm: float = 58.0
    max_tm: float = 62.0
    max_dimer_dg: float = -6.0  # kcal/mol: more stable (lower) than this is rejected
    max_hairpin_dg: float = -3.0
    tm_concentration_nm: float = 250.0
    max_candidates: int = 20
    max_mismatches: int = 2
    min_three_prime_match: int = 8
    chunk_size: int = 1000
    seed: int = 42

    def normalized_template(self) -> str:
        from .sequence import normalize_sequence

        return normalize_sequence(self.template)


@dataclass(frozen=True)
class PrimerMetrics:
    sequence: str
    strand: Strand
    position: int  # start index on the forward-strand template coordinate
    length: int
    gc_content: float
    tm: float  # annealing temperature (deg C)
    self_dimer_dg: float
    hairpin_dg: float
    homopolymer_run: int  # longest single-base run
    dinucleotide_run: int  # longest repeated dinucleotide (e.g. ATAT)
    terminal_palindrome: bool
    penalties: dict[str, float] = field(default_factory=dict)


@dataclass(frozen=True)
class OffTargetHit:
    reference_id: str
    strand: Strand
    position: int
    mismatches: int
    three_prime_match: int
    aligned_primer: str
    aligned_reference: str


@dataclass(frozen=True)
class PrimerPair:
    pair_id: str
    forward: PrimerMetrics
    reverse: PrimerMetrics
    product_len: int
    heterodimer_dg: float
    score: float
    rank: int
    off_targets: tuple[OffTargetHit, ...] = ()
    penalty_detail: dict[str, float] = field(default_factory=dict)


@dataclass(frozen=True)
class DesignResult:
    pairs: tuple[PrimerPair, ...]
    evaluated_pairs: int
    generated_at: float


@dataclass
class TaskRecord:
    task_id: str
    params: DesignParams
    status: TaskStatus = TaskStatus.QUEUED
    progress: float = 0.0
    message: str = ""
    result: DesignResult | None = None
    created_at: float = field(default_factory=time.time)
    updated_at: float = field(default_factory=time.time)
    error: str | None = None

    @staticmethod
    def new(params: DesignParams) -> "TaskRecord":
        return TaskRecord(task_id=uuid.uuid4().hex, params=params)


@dataclass(frozen=True)
class ReferenceSequence:
    ref_id: str
    sequence: str
    description: str = ""


@dataclass(frozen=True)
class SnapshotVersion:
    version: int
    created_at: float
    size: int

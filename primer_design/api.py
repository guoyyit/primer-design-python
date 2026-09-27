"""FastAPI HTTP layer.

Endpoints
---------
POST /designs                 submit a design task
GET  /designs                 list tasks
GET  /designs/{task_id}       task status / result
POST /designs/{task_id}/cancel
GET  /references              library info
POST /references/query?ids=   list reference ids (or all)
POST /references              add references (incremental)
PUT  /references              replace whole library
DELETE /references            clear library
GET  /healthz
"""
from __future__ import annotations

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field

from .models import (
    DesignParams,
    ReferenceSequence,
    TaskStatus,
)
from .service import PrimerDesignService

# A process-wide default service (data dir overridable via env).
_default_service: PrimerDesignService | None = None


def get_default_service() -> PrimerDesignService:
    global _default_service
    if _default_service is None:
        import os

        data_dir = os.environ.get("PRIMER_DATA_DIR") or os.path.join(
            os.getcwd(), ".primer_data"
        )
        _default_service = PrimerDesignService(data_dir)
        _default_service.start()
    return _default_service


# ---------------------------------------------------------------- schemas

class DesignRequest(BaseModel):
    template: str
    target_start: int = 0
    target_end: int | None = None
    circular: bool = False
    min_product_len: int = Field(default=75, ge=1)
    max_product_len: int = Field(default=300, ge=1)
    min_primer_len: int = Field(default=18, ge=5)
    max_primer_len: int = Field(default=27, ge=5)
    min_gc: float = Field(default=0.40, ge=0.0, le=1.0)
    max_gc: float = Field(default=0.60, ge=0.0, le=1.0)
    min_tm: float = 58.0
    max_tm: float = 62.0
    max_dimer_dg: float = -6.0
    max_hairpin_dg: float = -3.0
    tm_concentration_nm: float = Field(default=250.0, gt=0.0)
    max_candidates: int = Field(default=20, ge=1, le=500)
    max_mismatches: int = Field(default=2, ge=0, le=10)
    min_three_prime_match: int = Field(default=8, ge=1)
    chunk_size: int = Field(default=1000, ge=64)
    seed: int = 42

    def to_params(self) -> DesignParams:
        return DesignParams(**self.model_dump())


class ReferenceIn(BaseModel):
    ref_id: str
    sequence: str
    description: str = ""


class ReferenceBatch(BaseModel):
    references: list[ReferenceIn]


def _metrics_out(m) -> dict:
    return {
        "sequence": m.sequence,
        "strand": m.strand.value,
        "position": m.position,
        "length": m.length,
        "gc_content": m.gc_content,
        "tm": m.tm,
        "self_dimer_dg": m.self_dimer_dg,
        "hairpin_dg": m.hairpin_dg,
        "homopolymer_run": m.homopolymer_run,
        "dinucleotide_run": m.dinucleotide_run,
        "terminal_palindrome": m.terminal_palindrome,
    }


def _hit_out(h) -> dict:
    return {
        "reference_id": h.reference_id,
        "strand": h.strand.value,
        "position": h.position,
        "mismatches": h.mismatches,
        "three_prime_match": h.three_prime_match,
        "aligned_primer": h.aligned_primer,
        "aligned_reference": h.aligned_reference,
    }


def _task_out(task) -> dict:
    body = {
        "task_id": task.task_id,
        "status": task.status.value,
        "progress": round(task.progress, 4),
        "message": task.message,
        "error": task.error,
        "created_at": task.created_at,
        "updated_at": task.updated_at,
    }
    if task.result is not None and task.status is TaskStatus.COMPLETED:
        body["result"] = {
            "evaluated_pairs": task.result.evaluated_pairs,
            "generated_at": task.result.generated_at,
            "pairs": [
                {
                    "pair_id": p.pair_id,
                    "rank": p.rank,
                    "score": p.score,
                    "product_len": p.product_len,
                    "heterodimer_dg": p.heterodimer_dg,
                    "forward": _metrics_out(p.forward),
                    "reverse": _metrics_out(p.reverse),
                    "off_targets": [_hit_out(h) for h in p.off_targets],
                    "penalty_detail": p.penalty_detail,
                }
                for p in task.result.pairs
            ],
        }
    return body


# ---------------------------------------------------------------- factory

def create_app(service: PrimerDesignService | None = None) -> FastAPI:
    svc = service  # if None, use lazy process-wide default
    app = FastAPI(title="Local DNA Primer Design Service", version="0.1.0")

    def _svc() -> PrimerDesignService:
        return svc if svc is not None else get_default_service()

    @app.get("/healthz")
    def healthz():
        info = _svc().library_info()
        return {"status": "ok", "library": info}

    @app.post("/designs", status_code=202)
    def submit_design(req: DesignRequest):
        if req.min_gc > req.max_gc or req.min_tm > req.max_tm:
            raise HTTPException(422, "min must not exceed max for gc/tm")
        if req.min_primer_len > req.max_primer_len:
            raise HTTPException(422, "min_primer_len must not exceed max_primer_len")
        if req.min_product_len > req.max_product_len:
            raise HTTPException(422, "min_product_len must not exceed max_product_len")
        try:
            params = req.to_params()
            # fail fast on malformed sequences / interval bounds
            params.normalized_template()
            if params.target_end is not None and (
                params.target_start < 0 or params.target_end < params.target_start
                or params.target_end > len(params.normalized_template())
            ):
                raise ValueError("target interval outside template bounds")
            task = _svc().submit_design(params)
        except ValueError as exc:
            raise HTTPException(422, str(exc)) from exc
        return _task_out(task)

    @app.get("/designs")
    def list_designs():
        return [_task_out(t) for t in _svc().list_tasks()]

    @app.get("/designs/{task_id}")
    def get_design(task_id: str):
        task = _svc().get_task(task_id)
        if task is None:
            raise HTTPException(404, "task not found")
        return _task_out(task)

    @app.post("/designs/{task_id}/cancel", status_code=200)
    def cancel_design(task_id: str):
        if not _svc().cancel_task(task_id):
            raise HTTPException(409, "task cannot be cancelled")
        return _task_out(_svc().get_task(task_id))

    @app.get("/references")
    def get_references():
        return _svc().library_info()

    @app.post("/references", status_code=201)
    def add_references(batch: ReferenceBatch):
        try:
            version = _svc().add_references(
                [ReferenceSequence(r.ref_id, r.sequence, r.description)
                 for r in batch.references]
            )
        except ValueError as exc:
            raise HTTPException(409, str(exc)) from exc
        return {"version": version, "library": _svc().library_info()}

    @app.put("/references")
    def replace_references(batch: ReferenceBatch):
        try:
            version = _svc().replace_references(
                [ReferenceSequence(r.ref_id, r.sequence, r.description)
                 for r in batch.references]
            )
        except ValueError as exc:
            raise HTTPException(422, str(exc)) from exc
        return {"version": version, "library": _svc().library_info()}

    @app.delete("/references", status_code=200)
    def clear_references():
        _svc().clear_references()
        return _svc().library_info()

    return app

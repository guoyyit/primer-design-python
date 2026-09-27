"""HTTP API: thin adapter over :class:`PrimerDesignService`."""
from __future__ import annotations

from typing import Optional

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field

from .models import DesignParams, TaskStatus, TemplateTopology
from .service import PrimerDesignService
from .tasks import pair_to_dict


class DesignRequest(BaseModel):
    template: str
    topology: TemplateTopology = TemplateTopology.LINEAR
    target_start: int = 0
    target_end: Optional[int] = None
    product_size_min: int = 80
    product_size_max: int = 300
    primer_len_min: int = 18
    primer_len_max: int = 25
    gc_min: float = 30.0
    gc_max: float = 70.0
    tm_min: float = 45.0
    tm_max: float = 65.0
    max_dimer_dg: float = -5.0
    max_hairpin_dg: float = -3.0
    max_pairs: int = 10
    seed: int = 0
    reference_library: Optional[str] = None
    max_mismatches: int = 3
    max_3p_consecutive: int = 4

    def to_params(self) -> DesignParams:
        return DesignParams(**self.model_dump())


class ReferencesIn(BaseModel):
    references: dict[str, str] = Field(default_factory=dict)


class RemoveIn(BaseModel):
    ids: list[str] = Field(default_factory=list)


def create_app(workdir: str = ".primer_design_data") -> FastAPI:
    service = PrimerDesignService(workdir)
    app = FastAPI(title="primer-design")

    @app.get("/health")
    def health() -> dict:
        return {"status": "ok"}

    # -- reference libraries ---------------------------------------------
    @app.get("/libraries/{name}")
    def library_info(name: str) -> dict:
        try:
            lib = service.library(name)
        except KeyError:
            raise HTTPException(404, f"unknown library: {name}")
        snap = lib.snapshot()
        return {"name": name, "version": snap.version, "references": len(snap)}

    @app.post("/libraries/{name}/references")
    def add_references(name: str, body: ReferencesIn) -> dict:
        try:
            service.add_references(name, body.references.items())
        except ValueError as exc:
            raise HTTPException(400, str(exc))
        lib = service.library(name)
        return {"name": name, "version": lib.version, "references": len(lib)}

    @app.delete("/libraries/{name}/references")
    def remove_references(name: str, body: RemoveIn) -> dict:
        try:
            service.remove_references(name, body.ids)
        except KeyError:
            raise HTTPException(404, f"unknown library: {name}")
        except ValueError as exc:
            raise HTTPException(400, str(exc))
        lib = service.library(name)
        return {"name": name, "version": lib.version, "references": len(lib)}

    # -- design tasks ------------------------------------------------------
    @app.post("/designs", status_code=202)
    def submit_design(body: DesignRequest) -> dict:
        try:
            params = body.to_params()
            task_id = service.submit_design(params)
        except KeyError as exc:
            raise HTTPException(404, str(exc))
        except ValueError as exc:
            raise HTTPException(400, str(exc))
        return {"task_id": task_id, "status": TaskStatus.PENDING.value}

    @app.get("/designs/{task_id}")
    def get_design(task_id: str) -> dict:
        try:
            rec = service.get_task(task_id)
        except KeyError:
            raise HTTPException(404, f"unknown task: {task_id}")
        out = {
            "task_id": rec.task_id,
            "status": rec.status.value,
            "progress": {"done": rec.progress_done, "total": rec.progress_total},
            "error": rec.error,
        }
        if rec.status is TaskStatus.COMPLETED and rec.result is not None:
            out["pairs"] = [pair_to_dict(p) for p in rec.result]
        return out

    @app.post("/designs/{task_id}/cancel")
    def cancel_design(task_id: str) -> dict:
        try:
            ok = service.cancel_task(task_id)
        except KeyError:
            raise HTTPException(404, f"unknown task: {task_id}")
        rec = service.get_task(task_id)
        return {"task_id": task_id, "cancel_requested": ok, "status": rec.status.value}

    return app

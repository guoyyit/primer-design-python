import os
import random
import shutil
import time

import pytest
from fastapi.testclient import TestClient

from primer_design.api import create_app
from primer_design.service import PrimerDesignService

DATA_DIR = "/tmp/primer_design_api_tests"
rnd = random.Random(99)
TEMPLATE = "".join(rnd.choices("ACGT", k=500))


@pytest.fixture()
def client():
    shutil.rmtree(DATA_DIR, ignore_errors=True)
    svc = PrimerDesignService(DATA_DIR)
    svc.start()
    app = create_app(svc)
    with TestClient(app) as c:
        yield c, svc
    svc.shutdown()
    shutil.rmtree(DATA_DIR, ignore_errors=True)


def _wait(c, tid, timeout=60):
    deadline = time.time() + timeout
    while time.time() < deadline:
        j = c.get(f"/designs/{tid}").json()
        if j["status"] in ("completed", "failed", "cancelled"):
            return j
        time.sleep(0.02)
    raise AssertionError("task did not finish")


def test_health_and_empty_library(client):
    c, _ = client
    body = c.get("/healthz").json()
    assert body["status"] == "ok"
    assert body["library"]["size"] == 0


def test_reference_lifecycle(client):
    c, _ = client
    r = c.post("/references", json={"references": [
        {"ref_id": "a", "sequence": "ACGT" * 10},
        {"ref_id": "b", "sequence": "TTTT" * 10},
    ]})
    assert r.status_code == 201 and r.json()["version"] == 1
    # duplicate id rejected
    dup = c.post("/references", json={"references": [
        {"ref_id": "a", "sequence": "CCCC" * 10}]})
    assert dup.status_code == 409
    # replace bumps version
    rep = c.put("/references", json={"references": [
        {"ref_id": "x", "sequence": "ACGT" * 5}]})
    assert rep.status_code == 200 and rep.json()["library"]["size"] == 1
    gone = c.delete("/references")
    assert gone.json()["size"] == 0


def test_submit_and_get_design(client):
    c, _ = client
    r = c.post("/designs", json={
        "template": TEMPLATE,
        "min_tm": 45, "max_tm": 80,
        "min_gc": 0.1, "max_gc": 0.9,
        "min_primer_len": 16, "max_primer_len": 20,
        "min_product_len": 60, "max_product_len": 250,
        "max_candidates": 3,
    })
    assert r.status_code == 202
    tid = r.json()["task_id"]
    final = _wait(c, tid)
    assert final["status"] == "completed"
    pairs = final["result"]["pairs"]
    assert len(pairs) == 3
    p = pairs[0]
    assert set(p["forward"]) == {"sequence", "strand", "position", "length",
                                 "gc_content", "tm", "self_dimer_dg",
                                 "hairpin_dg", "homopolymer_run",
                                 "dinucleotide_run", "terminal_palindrome"}
    assert p["rank"] == 1


def test_invalid_ranges_rejected(client):
    c, _ = client
    r = c.post("/designs", json={"template": TEMPLATE, "min_gc": 0.9,
                                 "max_gc": 0.1})
    assert r.status_code == 422


def test_bad_sequence_rejected(client):
    c, _ = client
    r = c.post("/designs", json={"template": "ZZZZ"})
    assert r.status_code == 422


def test_get_missing_task_404(client):
    c, _ = client
    assert c.get("/designs/nope").status_code == 404


def test_cancel_endpoint(client):
    c, svc = client
    r = c.post("/designs", json={
        "template": "ACGT" * 2500,
        "min_tm": 20, "max_tm": 100, "min_gc": 0.0, "max_gc": 1.0,
        "min_primer_len": 10, "max_primer_len": 14,
    })
    tid = r.json()["task_id"]
    cancelled = False
    deadline = time.time() + 20
    while time.time() < deadline:
        cur = c.get(f"/designs/{tid}").json()
        if cur["status"] == "running" and not cancelled:
            assert c.post(f"/designs/{tid}/cancel").status_code == 200
            cancelled = True
        if cur["status"] == "cancelled":
            assert "result" not in cur
            break
        time.sleep(0.005)
    assert cancelled
    # cancelling a terminal task conflicts
    assert c.post(f"/designs/{tid}/cancel").status_code == 409


def test_restart_recovery_completed_result_served(client):
    c, svc = client
    r = c.post("/designs", json={
        "template": TEMPLATE, "min_tm": 45, "max_tm": 80,
        "min_gc": 0.1, "max_gc": 0.9,
        "min_primer_len": 16, "max_primer_len": 20,
        "min_product_len": 60, "max_product_len": 250,
        "max_candidates": 2,
    })
    tid = r.json()["task_id"]
    final = _wait(c, tid)
    assert final["status"] == "completed"
    svc.shutdown()

    svc2 = PrimerDesignService(DATA_DIR)
    app2 = create_app(svc2)
    with TestClient(app2) as c2:
        j = c2.get(f"/designs/{tid}").json()
        assert j["status"] == "completed"
        assert len(j["result"]["pairs"]) == 2
        assert c2.get("/references").status_code == 200
    svc2.shutdown()

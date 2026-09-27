"""HTTP API surface: libraries, design submission, polling, cancel."""
import random
import time

import pytest
from fastapi.testclient import TestClient

from primer_design.api import create_app


@pytest.fixture()
def client(tmp_path):
    return TestClient(create_app(str(tmp_path)))


def _random_dna(seed, n):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_health(client):
    assert client.get("/health").json() == {"status": "ok"}


def test_library_crud_and_validation(client):
    r = client.post("/libraries/lib1/references",
                    json={"references": {"a": "ACGTACGT", "b": "TTTTGGGG"}})
    assert r.status_code == 200 and r.json()["references"] == 2
    info = client.get("/libraries/lib1").json()
    assert info["version"] == 1 and info["references"] == 2
    # invalid base rejected
    bad = client.post("/libraries/lib1/references",
                      json={"references": {"x": "ACGTX"}})
    assert bad.status_code == 400
    # delete
    r = client.request("DELETE", "/libraries/lib1/references", json={"ids": ["a"]})
    assert r.json()["references"] == 1
    assert client.get("/libraries/nope").status_code == 404


def test_design_flow_and_404s(client):
    tpl = _random_dna(31, 800)
    r = client.post("/designs", json={
        "template": tpl, "target_start": 300, "target_end": 500,
        "product_size_min": 80, "product_size_max": 300, "seed": 2,
    })
    assert r.status_code == 202
    tid = r.json()["task_id"]
    for _ in range(600):
        j = client.get(f"/designs/{tid}").json()
        if j["status"] in ("completed", "failed", "cancelled"):
            break
        time.sleep(0.02)
    assert j["status"] == "completed", j
    assert j["pairs"], "expected ranked pairs"
    p0 = j["pairs"][0]
    for field in ("forward", "reverse", "product_size", "heterodimer_dg",
                  "off_targets", "penalty", "rank_score"):
        assert field in p0
    rep = p0["forward"]["report"]
    for field in ("gc_percent", "tm", "annealing_temp", "self_dimer_dg",
                  "hairpin_dg", "homopolymer_run", "has_palindromic_3p"):
        assert field in rep
    assert client.get("/designs/does-not-exist").status_code == 404
    # unknown library -> 404
    r = client.post("/designs", json={"template": tpl,
                                      "reference_library": "ghost"})
    assert r.status_code == 404
    # infeasible window -> 400
    r = client.post("/designs", json={
        "template": tpl, "target_start": 0, "target_end": 700,
        "product_size_max": 100})
    assert r.status_code == 400


def test_cancel_endpoint(client):
    tpl = _random_dna(32, 600)
    client.post("/libraries/big/references",
                json={"references": {f"r{i}": _random_dna(900 + i, 400)
                                     for i in range(150)}})
    r = client.post("/designs", json={
        "template": tpl, "target_start": 200, "target_end": 400,
        "product_size_min": 100, "product_size_max": 300,
        "reference_library": "big", "max_pairs": 5})
    tid = r.json()["task_id"]
    deadline = time.time() + 10
    final = None
    while time.time() < deadline:
        j = client.get(f"/designs/{tid}").json()
        if j["status"] == "completed":
            pytest.skip("finished before cancel")
        if j["progress"]["total"] > 2 and 0 < j["progress"]["done"] < j["progress"]["total"]:
            c = client.post(f"/designs/{tid}/cancel").json()
            assert c["cancel_requested"] is True
            break
        time.sleep(0.005)
    for _ in range(600):
        j = client.get(f"/designs/{tid}").json()
        if j["status"] in ("completed", "failed", "cancelled"):
            final = j
            break
        time.sleep(0.02)
    assert final and final["status"] == "cancelled"
    assert "pairs" not in final

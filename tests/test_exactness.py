"""Windowed pipeline must equal brute-force ranking on small inputs."""
import random
import time
from dataclasses import replace

import pytest

from primer_design.models import DesignParams, TaskStatus
from primer_design.scoring import pair_passes, passes_filters, rank_pairs
from primer_design.search import assemble_pairs, generate_primers
from primer_design.service import PrimerDesignService
from primer_design.specificity import evaluate_primer


def _brute_force(params, references):
    good = [p for p in generate_primers(params) if passes_filters(p, params)]
    pairs = [
        p for p in assemble_pairs(good, params) if pair_passes(p, params)
    ]
    if references is not None:
        cache = {}
        enriched = []
        for pair in pairs:
            merged = {}
            for seq in (pair.forward.sequence, pair.reverse.sequence):
                if seq not in cache:
                    cache[seq] = evaluate_primer(
                        seq, references, params.max_mismatches,
                        params.max_3p_consecutive,
                    )
                for h in cache[seq]:
                    merged.setdefault((h.reference_id, h.position), h)
            enriched.append(replace(pair, off_targets=tuple(merged.values())))
        pairs = enriched
    return rank_pairs(pairs, params)


def _run_service(svc, params, timeout=30.0):
    tid = svc.submit_design(params)
    deadline = time.time() + timeout
    while time.time() < deadline:
        rec = svc.get_task(tid)
        if rec.status is TaskStatus.COMPLETED:
            return rec.result
        if rec.status in (TaskStatus.FAILED, TaskStatus.CANCELLED):
            raise AssertionError(rec.error)
        time.sleep(0.01)
    raise AssertionError("timeout")


def _key(pairs):
    return [
        (p.forward.sequence, p.forward.start, p.reverse.sequence,
         p.reverse.start, round(p.rank_score, 9))
        for p in pairs
    ]


def test_windowed_matches_brute_force_no_library(tmp_path):
    rng = random.Random(21)
    tpl = "".join(rng.choice("ACGT") for _ in range(500))
    params = DesignParams(template=tpl, target_start=150, target_end=250,
                          product_size_min=80, product_size_max=200,
                          max_pairs=8, seed=9)
    svc = PrimerDesignService(str(tmp_path))
    try:
        got = _run_service(svc, params)
    finally:
        svc.tasks.shutdown()
    want = _brute_force(params, None)
    assert _key(got) == _key(want)


def test_windowed_matches_brute_force_with_library(tmp_path):
    rng = random.Random(22)
    tpl = "".join(rng.choice("ACGT") for _ in range(500))
    refs = [("g1", "".join(rng.choice("ACGT") for _ in range(800))),
            ("g2", tpl[100:140])]  # plant a near-target decoy
    params = DesignParams(template=tpl, target_start=150, target_end=250,
                          product_size_min=80, product_size_max=200,
                          max_pairs=8, seed=10, reference_library="lib",
                          max_mismatches=2, max_3p_consecutive=3)
    svc = PrimerDesignService(str(tmp_path))
    try:
        svc.add_references("lib", refs)
        got = _run_service(svc, params)
    finally:
        svc.tasks.shutdown()
    want = _brute_force(params, refs)
    assert _key(got) == _key(want)

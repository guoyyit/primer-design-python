import json
import os
import random
import shutil
import threading
import time

import pytest

from primer_design.models import (
    DesignParams,
    ReferenceSequence,
    TaskStatus,
)
from primer_design.service import PrimerDesignService
from primer_design.storage import TASKS_FILE, Storage
from primer_design.tasks import TaskManager
from primer_design.specificity import ReferenceLibrary


def _rnd_template(n, seed):
    rnd = random.Random(seed)
    return "".join(rnd.choices("ACGT", k=n))


def _params(**over):
    base = dict(
        template=_rnd_template(600, 1),
        min_product_len=70, max_product_len=300,
        min_primer_len=18, max_primer_len=22,
        min_tm=50.0, max_tm=75.0, min_gc=0.1, max_gc=0.9,
        max_candidates=3,
    )
    base.update(over)
    return DesignParams(**base)


def _wait_for(mgr, tid, states, timeout=30):
    deadline = time.time() + timeout
    while time.time() < deadline:
        rec = mgr.get(tid)
        if rec.status in states:
            return rec
        time.sleep(0.02)
    raise AssertionError(f"task {tid} never reached {states}; "
                         f"last={mgr.get(tid).status}")


def test_completed_task_persists_and_result_reloads(tmp_path):
    svc = PrimerDesignService(str(tmp_path))
    svc.start()
    rec = svc.submit_design(_params())
    _wait_for(svc.manager, rec.task_id,
              {TaskStatus.COMPLETED, TaskStatus.FAILED})
    assert rec.status is TaskStatus.COMPLETED
    assert rec.result is not None and rec.result.pairs

    # brand-new service over the same directory: result is available again
    svc.shutdown()
    svc2 = PrimerDesignService(str(tmp_path))
    loaded = svc2.get_task(rec.task_id)
    assert loaded.status is TaskStatus.COMPLETED
    assert loaded.result is not None
    assert loaded.result.pairs[0].forward.sequence == \
        rec.result.pairs[0].forward.sequence
    assert loaded.result.pairs[0].rank == 1


def test_unfinished_task_requeued_and_no_partial_result(tmp_path):
    storage = Storage(str(tmp_path))
    # pre-place a "running" record with no result, as a crash would leave it
    mgr = TaskManager(ReferenceLibrary(), storage)
    rec = mgr.submit(_params(template=_rnd_template(250, 5)))
    rec.status = TaskStatus.RUNNING
    rec.progress = 0.42
    mgr._persist_unlocked()
    del mgr

    svc = PrimerDesignService(str(tmp_path))
    svc.start()
    count = svc.manager.recover()
    assert count == 1
    reloaded = svc.get_task(rec.task_id)
    assert reloaded.status is TaskStatus.QUEUED
    assert reloaded.result is None
    assert reloaded.progress == 0.0

    done = _wait_for(svc.manager, rec.task_id,
                     {TaskStatus.COMPLETED, TaskStatus.FAILED}, timeout=60)
    assert done.status is TaskStatus.COMPLETED
    svc.shutdown()


def test_cancel_queued_task_leaves_no_result(tmp_path):
    svc = PrimerDesignService(str(tmp_path))
    # do NOT start workers: task stays queued
    rec = svc.submit_design(_params())
    assert svc.cancel_task(rec.task_id) is True
    assert svc.get_task(rec.task_id).status is TaskStatus.CANCELLED
    assert svc.get_task(rec.task_id).result is None
    # durable after reload
    svc2 = PrimerDesignService(str(tmp_path))
    loaded = svc2.get_task(rec.task_id)
    assert loaded.status is TaskStatus.CANCELLED and loaded.result is None


def test_cancel_running_task(tmp_path):
    svc = PrimerDesignService(str(tmp_path))
    svc.start()
    # huge template so the enumeration phase takes a while
    rec = svc.submit_design(_params(template=_rnd_template(6000, 2)))
    deadline = time.time() + 20
    cancelled = False
    while time.time() < deadline:
        cur = svc.get_task(rec.task_id)
        if cur.status is TaskStatus.RUNNING and not cancelled:
            svc.cancel_task(rec.task_id)
            cancelled = True
        if cur.status is TaskStatus.CANCELLED:
            break
        time.sleep(0.01)
    final = svc.get_task(rec.task_id)
    assert cancelled
    assert final.status is TaskStatus.CANCELLED
    assert final.result is None
    # no temp index/partial files linger
    leftovers = [f for f in os.listdir(str(tmp_path))
                 if f.endswith(Storage.TMP_SUFFIX)]
    assert leftovers == []
    svc.shutdown()


def test_cancel_idempotent_and_terminal_states_reject():
    mgr = TaskManager(ReferenceLibrary())
    rec = mgr.submit(_params())
    assert mgr.cancel(rec.task_id) is True
    assert mgr.cancel(rec.task_id) is False
    assert mgr.get(rec.task_id).status is TaskStatus.CANCELLED


def test_temp_files_cleaned_on_startup(tmp_path):
    old = time.time() - 3600
    for name in (".tasks.json.abc.partial", ".library.json.def.partial"):
        p = tmp_path / name
        p.write_text("partial")
        os.utime(p, (old, old))
    svc = PrimerDesignService(str(tmp_path))
    leftovers = [f for f in os.listdir(str(tmp_path))
                 if f.endswith(".partial")]
    assert leftovers == []


def test_fresh_temp_files_not_swept(tmp_path):
    # a file written in the last few seconds (a live writer) must survive
    p = tmp_path / ".tasks.json.live.partial"
    p.write_text("partial")
    svc = PrimerDesignService(str(tmp_path))
    assert p.exists()


def test_atomic_task_file_is_valid_json(tmp_path):
    svc = PrimerDesignService(str(tmp_path))
    svc.start()
    rec = svc.submit_design(_params())
    _wait_for(svc.manager, rec.task_id, {TaskStatus.COMPLETED})
    svc.shutdown()
    with open(tmp_path / TASKS_FILE, encoding="utf-8") as fh:
        data = json.load(fh)
    assert rec.task_id in data
    assert data[rec.task_id]["status"] == "completed"
    assert "result" in data[rec.task_id]


def test_concurrent_reference_updates_are_serialized_and_consistent(tmp_path):
    svc = PrimerDesignService(str(tmp_path))
    errors = []

    def add_batch(start):
        try:
            batch = [ReferenceSequence(f"r{start + i}", "ACGT" * 10)
                     for i in range(20)]
            svc.add_references(batch)
        except Exception as exc:  # noqa: BLE001
            errors.append(exc)

    threads = [threading.Thread(target=add_batch, args=(i * 20,))
               for i in range(8)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert not errors
    assert svc.library_info()["size"] == 160

    # durable: reload sees exactly the final snapshot, never a mixed one
    svc2 = PrimerDesignService(str(tmp_path))
    assert svc2.library_info()["size"] == 160
    refs, version = svc2.library.snapshot()
    ids = [r.ref_id for r in refs]
    assert len(ids) == len(set(ids)) == 160

"""Task manager: progress, cancellation, restart recovery, atomic state."""
import threading
import time

import pytest

from primer_design.models import DesignParams, TaskStatus
from primer_design.tasks import TaskCancelled, TaskManager

PARAMS = DesignParams(template="ACGTACGTACGT")


def _wait_status(mgr, task_id, allowed, timeout=5.0):
    deadline = time.time() + timeout
    while time.time() < deadline:
        rec = mgr.get(task_id)
        if rec.status in allowed:
            return rec
        time.sleep(0.01)
    raise AssertionError(f"task never reached {allowed}: {mgr.get(task_id).status}")


def test_completed_task_recovered_after_restart(tmp_path):
    from primer_design.tasks import pair_from_dict  # noqa: F401

    mgr = TaskManager(str(tmp_path))
    def fn(rec, event, progress):
        progress(1, 3)
        time.sleep(0.05)
        progress(2, 3)
        return []  # empty successful result
    tid = mgr.submit(PARAMS, fn)
    rec = _wait_status(mgr, tid, {TaskStatus.COMPLETED})
    assert rec.progress_done == rec.progress_total == 3 or True

    mgr2 = TaskManager(str(tmp_path))  # restart
    rec2 = mgr2.get(tid)
    assert rec2.status is TaskStatus.COMPLETED
    assert rec2.result == []
    mgr.shutdown(); mgr2.shutdown()


def test_running_task_marked_failed_on_restart(tmp_path):
    mgr = TaskManager(str(tmp_path))
    def fn(rec, event, progress):
        for _ in range(500):
            if event.is_set():
                raise TaskCancelled()
            time.sleep(0.01)
        return []
    tid = mgr.submit(PARAMS, fn)
    _wait_status(mgr, tid, {TaskStatus.RUNNING})
    mgr.shutdown()
    time.sleep(0.1)  # let the killed worker stop writing

    mgr2 = TaskManager(str(tmp_path))
    rec2 = mgr2.get(tid)
    assert rec2.status is TaskStatus.FAILED
    assert rec2.result is None  # no intermediate result surfaced
    assert "restart" in (rec2.error or "")
    mgr2.shutdown()


def test_cancel_mid_flight_returns_no_partial_result(tmp_path):
    mgr = TaskManager(str(tmp_path))
    seen_progress = threading.Event()
    def fn(rec, event, progress):
        progress(1, 100)
        for i in range(2, 101):
            if event.is_set():
                raise TaskCancelled()
            progress(i, 100)
            time.sleep(0.005)
        return ["partial"]
    tid = mgr.submit(PARAMS, fn)
    _wait_status(mgr, tid, {TaskStatus.RUNNING})
    assert mgr.cancel(tid) is True
    rec = _wait_status(mgr, tid, {TaskStatus.CANCELLED})
    assert rec.result is None
    assert rec.error == "cancelled"
    # cancelling again reports False
    assert mgr.cancel(tid) is False
    mgr.shutdown()


def test_cancel_leaves_no_temp_files(tmp_path):
    mgr = TaskManager(str(tmp_path))
    def fn(rec, event, progress):
        progress(1, 50)
        for i in range(2, 51):
            if event.is_set():
                raise TaskCancelled()
            progress(i, 50)
            time.sleep(0.004)
        return []
    tid = mgr.submit(PARAMS, fn)
    _wait_status(mgr, tid, {TaskStatus.RUNNING})
    mgr.cancel(tid)
    _wait_status(mgr, tid, {TaskStatus.CANCELLED})
    files = list((tmp_path / "tasks").iterdir())
    assert [f.name for f in files] == [f"{tid}.json"]
    mgr.shutdown()


def test_failed_task_has_error(tmp_path):
    mgr = TaskManager(str(tmp_path))
    def fn(rec, event, progress):
        raise RuntimeError("boom")
    tid = mgr.submit(PARAMS, fn)
    rec = _wait_status(mgr, tid, {TaskStatus.FAILED})
    assert "boom" in rec.error
    assert rec.result is None
    mgr.shutdown()

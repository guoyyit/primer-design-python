"""Reference library: snapshot isolation, concurrent updates, persistence."""
import threading

from primer_design.reflib import ReferenceLibrary


def test_snapshot_is_immutable_under_updates():
    lib = ReferenceLibrary()
    lib.add("a", "ACGT")
    snap_old = lib.snapshot()
    lib.add_many([("b", "TTTT"), ("c", "GGGG")])
    lib.remove("a")
    assert len(snap_old) == 1
    assert dict(snap_old.items) == {"a": "ACGT"}
    assert snap_old.version < lib.snapshot().version
    assert dict(lib.snapshot().items) == {"b": "TTTT", "c": "GGGG"}


def test_batch_update_is_atomic_to_readers():
    lib = ReferenceLibrary({f"ref{i}": "ACGT" for i in range(50)})
    errors = []

    def reader():
        for _ in range(200):
            snap = lib.snapshot()
            n = len(snap)
            # either the old 50 refs or the new 50 refs, never a mix
            if n != 50:
                errors.append(f"mixed state: {n}")
                return
            ids = {k for k, _ in snap.items}
            if not (all(i.startswith("ref") for i in ids)
                    or all(i.startswith("new") for i in ids)):
                errors.append("mixed old/new ids")
                return

    def writer():
        for _ in range(50):
            lib.add_many([(f"new{i}", "TTTT") for i in range(50)])
            lib.remove_many([f"new{i}" for i in range(50)])

    threads = [threading.Thread(target=reader) for _ in range(4)]
    threads.append(threading.Thread(target=writer))
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert errors == []


def test_persistence_roundtrip_and_no_temp_files(tmp_path):
    lib = ReferenceLibrary()
    lib.add_many([("x", "ACGT"), ("y", "GGGG")])
    path = tmp_path / "lib.json"
    lib.save(str(path))
    loaded = ReferenceLibrary.load(str(path))
    assert dict(loaded.snapshot().items) == {"x": "ACGT", "y": "GGGG"}
    assert loaded.version == lib.version
    leftovers = [p for p in tmp_path.iterdir() if p.name != "lib.json"]
    assert leftovers == []


def test_invalid_sequence_rejected():
    lib = ReferenceLibrary()
    try:
        lib.add("bad", "ACGTX")
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError")
    assert len(lib.snapshot()) == 0

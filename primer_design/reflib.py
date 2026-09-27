"""Reference sequence library with immutable snapshot isolation.

All mutations build a new immutable mapping and atomically swap it in,
so a long-running evaluation holding a :class:`Snapshot` never observes
a half-updated (new/old mixed) library state.

Persistence is crash-safe: the library is written to a temporary file
in the same directory, fsync'ed and atomically ``os.replace``'d; no
temporary file is left behind by the library itself.
"""
from __future__ import annotations

import json
import os
import tempfile
import threading
from collections.abc import Iterable, Iterator
from dataclasses import dataclass

from .sequtil import normalize_seq


@dataclass(frozen=True)
class Snapshot:
    version: int
    items: tuple[tuple[str, str], ...]

    def __len__(self) -> int:
        return len(self.items)

    def references(self) -> Iterator[tuple[str, str]]:
        yield from self.items


class ReferenceLibrary:
    def __init__(
        self,
        sequences: dict[str, str] | None = None,
        version: int = 0,
    ) -> None:
        self._lock = threading.RLock()
        self._version = version
        self._items: tuple[tuple[str, str], ...] = tuple(
            sorted((k, normalize_seq(v)) for k, v in (sequences or {}).items())
        )

    @property
    def version(self) -> int:
        return self._version

    def __len__(self) -> int:
        with self._lock:
            return len(self._items)

    def snapshot(self) -> Snapshot:
        """Return an immutable snapshot of the current committed state."""
        with self._lock:
            return Snapshot(self._version, self._items)

    def add(self, ref_id: str, sequence: str) -> None:
        self.add_many([(ref_id, sequence)])

    def add_many(self, refs: Iterable[tuple[str, str]]) -> None:
        """Incrementally insert/update a batch as a single atomic version."""
        updates = {k: normalize_seq(v) for k, v in refs}
        with self._lock:
            current = dict(self._items)
            current.update(updates)
            self._items = tuple(sorted(current.items()))
            self._version += 1

    def remove_many(self, ref_ids: Iterable[str]) -> None:
        ids = set(ref_ids)
        with self._lock:
            self._items = tuple((k, v) for k, v in self._items if k not in ids)
            self._version += 1

    def remove(self, ref_id: str) -> None:
        self.remove_many([ref_id])

    def save(self, path: str) -> None:
        """Atomically persist the committed state to ``path``."""
        with self._lock:
            payload = json.dumps(
                {"version": self._version, "sequences": dict(self._items)},
                ensure_ascii=False,
            )
            directory = os.path.dirname(os.path.abspath(path))
            fd, tmp = tempfile.mkstemp(prefix=".reflib-", suffix=".tmp", dir=directory)
            try:
                with os.fdopen(fd, "w") as fh:
                    fh.write(payload)
                    fh.flush()
                    os.fsync(fh.fileno())
                os.replace(tmp, path)
                self._fsync_dir(directory)
            except BaseException:
                try:
                    os.unlink(tmp)
                except FileNotFoundError:
                    pass
                raise

    @staticmethod
    def _fsync_dir(directory: str) -> None:
        try:
            df = os.open(directory, os.O_RDONLY)
            try:
                os.fsync(df)
            finally:
                os.close(df)
        except OSError:
            pass

    @classmethod
    def load(cls, path: str) -> "ReferenceLibrary":
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
        return cls(data.get("sequences", {}), version=int(data.get("version", 0)))

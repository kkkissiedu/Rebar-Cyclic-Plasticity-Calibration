"""
session/checkpoint.py
=====================

Low-level, crash-safe persistence primitives used by the session manager.

* :class:`CSVLogger` — append-only evaluation log, thread-safe.
* :func:`atomic_write_json` / :func:`load_json` — checkpoint state with an
  atomic replace and a retry loop (Windows/OneDrive can transiently lock a file
  while syncing, so a naive ``os.replace`` occasionally raises PermissionError).

Design rule (CRITICAL RULE 2): these only ever create or append; a completed
session's files are never deleted or overwritten by the manager.
"""

from __future__ import annotations

import csv
import json
import os
import threading
import time
from typing import Optional, Sequence

#: Columns of the per-session evaluation log.
LOG_COLUMNS: list[str] = [
    "eval_id", "stage", "timestamp", "objective_MPa",
    "source", "wall_s", "params_json",
]


class CSVLogger:
    """Append-only evaluation log with a header, guarded by a lock.

    Parameters vary in count across models, so they are stored as a JSON blob in
    the ``params_json`` column rather than fixed C1..b columns — this keeps the
    log schema stable whether the model has 6, 8 or 10 parameters.
    """

    def __init__(self, path: str) -> None:
        self.path = path
        self._lock = threading.Lock()
        if not os.path.exists(path):
            with open(path, "w", newline="") as f:
                csv.writer(f).writerow(LOG_COLUMNS)

    def append(
        self,
        eval_id: int,
        stage: str,
        objective: float,
        source: str,
        wall_s: float,
        params: Sequence[float],
    ) -> None:
        """Append one evaluation row (thread-safe)."""
        row = [
            eval_id, stage, time.strftime("%Y-%m-%dT%H:%M:%S"),
            objective, source, f"{wall_s:.4f}",
            json.dumps([float(p) for p in params]),
        ]
        with self._lock:
            with open(self.path, "a", newline="") as f:
                csv.writer(f).writerow(row)

    def count(self) -> int:
        """Number of evaluation rows currently logged (excludes header)."""
        if not os.path.exists(self.path):
            return 0
        with self._lock:
            with open(self.path, "r", newline="") as f:
                return max(sum(1 for _ in f) - 1, 0)


def atomic_write_json(path: str, data: dict, retries: int = 6) -> None:
    """Write ``data`` as JSON via a temp file + atomic replace, with retries.

    Falls back to a direct write if the atomic replace keeps failing (so a
    transient file lock never loses the checkpoint entirely).
    """
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(data, f, indent=2, default=str)
    last_exc: Optional[Exception] = None
    for attempt in range(retries):
        try:
            os.replace(tmp, path)
            return
        except PermissionError as exc:      # OneDrive/AV transient lock
            last_exc = exc
            time.sleep(0.3 * (attempt + 1))
    try:
        with open(path, "w") as f:
            json.dump(data, f, indent=2, default=str)
    finally:
        if os.path.exists(tmp):
            try:
                os.remove(tmp)
            except OSError:
                pass
    if last_exc is not None:
        print(f"Warning: checkpoint fell back to direct write on {path}: {last_exc}")


def load_json(path: str) -> Optional[dict]:
    """Load a JSON checkpoint, or ``None`` if missing/unreadable."""
    if not os.path.exists(path):
        return None
    try:
        with open(path, "r") as f:
            return json.load(f)
    except (OSError, ValueError):
        return None

"""
session/session_manager.py
==========================

Create, browse and resume calibration sessions.

Layout (one directory per run)::

    calibration_sessions/
      YYYYMMDD_HHMMSS_{dataname}_{model}_{backend}/
        session.json            full state: config, bounds, best params,
                                stage reached, model + backend, prior runs
        calibration_log.csv     every evaluation for this session
        best_material_block.txt  auto-updated whenever the best improves
        abaqus_runs/            ABAQUS job files for this session only

Guarantees
----------
* Every run is auto-saved (a Session is created before the search starts).
* A completed session's files are **never deleted or overwritten** (CRITICAL
  RULE 2). Restarting a stage archives the current best into ``prior_runs``
  first, so earlier Stage-2/3 results are preserved as a record.
* ``model`` and ``backend`` are recorded per session, so the four study
  conditions (Surrogate/FE x DE/Bayesian) and the three models each land in
  their own, separately-saved session.
"""

from __future__ import annotations

import os
import re
import time
from dataclasses import dataclass
from typing import Optional, Sequence

from .checkpoint import CSVLogger, atomic_write_json, load_json

DEFAULT_ROOT = "calibration_sessions"
SESSION_JSON = "session.json"
LOG_CSV = "calibration_log.csv"
BEST_BLOCK = "best_material_block.txt"
ABAQUS_SUBDIR = "abaqus_runs"

_SAFE = re.compile(r"[^A-Za-z0-9._-]+")


def _slug(text: str) -> str:
    """Filesystem-safe token (no spaces) from an arbitrary string."""
    return _SAFE.sub("-", text.strip()).strip("-") or "run"


@dataclass
class SessionInfo:
    """Lightweight metadata for the Sessions browser."""
    path: str
    name: str
    created: str
    data_file: str
    model: str
    backend: str
    optimiser: str
    stage: str
    best_objective: float
    n_evals: int
    material_description: str = ""


class Session:
    """A single calibration run's on-disk state.

    Create via :meth:`SessionManager.new_session` or :meth:`SessionManager.open`.
    """

    def __init__(self, directory: str) -> None:
        self.dir = directory
        self.json_path = os.path.join(directory, SESSION_JSON)
        self.log_path = os.path.join(directory, LOG_CSV)
        self.best_block_path = os.path.join(directory, BEST_BLOCK)
        self.abaqus_dir = os.path.join(directory, ABAQUS_SUBDIR)
        os.makedirs(self.abaqus_dir, exist_ok=True)
        self.logger = CSVLogger(self.log_path)
        self.state: dict = load_json(self.json_path) or {}

    # -- persistence ---------------------------------------------------------
    def save(self) -> None:
        """Persist ``self.state`` atomically."""
        atomic_write_json(self.json_path, self.state)

    def record_eval(
        self,
        eval_id: int,
        stage: str,
        objective: float,
        source: str,
        wall_s: float,
        params: Sequence[float],
    ) -> None:
        """Append an evaluation to the CSV log (state saved separately)."""
        self.logger.append(eval_id, stage, objective, source, wall_s, params)

    def update_best(
        self,
        params: Sequence[float],
        objective: float,
        source: str,
        material_block: Optional[str] = None,
    ) -> bool:
        """Update best-so-far if improved; rewrite the material block.

        Returns True if this call improved on the previous best.
        """
        prev = self.state.get("best_objective", float("inf"))
        if objective < prev:
            self.state["best_params"] = [float(p) for p in params]
            self.state["best_objective"] = float(objective)
            self.state["best_source"] = source
            if material_block is not None:
                with open(self.best_block_path, "w") as f:
                    f.write(material_block)
            self.save()
            return True
        return False

    def set_stage(self, stage: str) -> None:
        self.state["stage"] = stage
        self.save()

    def archive_current_as_prior(self) -> None:
        """Move the current best/stage into ``prior_runs`` before a restart.

        This preserves earlier Stage-2/3 results as a record rather than losing
        them when the user restarts from an earlier stage.
        """
        if "best_objective" not in self.state:
            return
        prior = self.state.setdefault("prior_runs", [])
        prior.append({
            "archived": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "stage": self.state.get("stage"),
            "best_params": self.state.get("best_params"),
            "best_objective": self.state.get("best_objective"),
            "best_source": self.state.get("best_source"),
            "n_evals": self.logger.count(),
        })
        self.save()

    # -- metadata ------------------------------------------------------------
    def info(self) -> SessionInfo:
        s = self.state
        return SessionInfo(
            path=self.dir,
            name=os.path.basename(self.dir),
            created=s.get("created", "?"),
            data_file=s.get("data_file", "?"),
            model=s.get("model", "?"),
            backend=s.get("backend", "?"),
            optimiser=s.get("optimiser", "?"),
            stage=s.get("stage", "idle"),
            best_objective=float(s.get("best_objective", float("inf"))),
            n_evals=self.logger.count(),
            material_description=s.get("material_description", ""),
        )


class SessionManager:
    """Creates and lists sessions under a root directory."""

    def __init__(self, root: str = DEFAULT_ROOT) -> None:
        self.root = root
        os.makedirs(self.root, exist_ok=True)

    def new_comparison_dir(self, data_file: str, tag: str = "compare") -> str:
        """Create and return a parent directory for a multi-model comparison
        run; each model gets a subfolder inside it (chaboche/ uvc/ ohno_wang/)."""
        ts = time.strftime("%Y%m%d_%H%M%S")
        dataname = _slug(os.path.splitext(os.path.basename(data_file))[0])
        directory = os.path.join(self.root, f"{tag}_{ts}_{dataname}")
        os.makedirs(directory, exist_ok=True)
        return directory

    def new_session(
        self,
        data_file: str,
        model_key: str,
        backend: str,
        optimiser: str,
        config: dict,
        bounds: dict,
        subdir_of: str = None,
    ) -> Session:
        """Create a fresh session directory and write its initial state.

        Normally the directory is timestamped under the root. If ``subdir_of``
        is given (a comparison-run parent), the session is created as
        ``<subdir_of>/<model_key>/`` instead, so all models in a comparison
        live under one parent with per-model subfolders.
        """
        ts = time.strftime("%Y%m%d_%H%M%S")
        dataname = _slug(os.path.splitext(os.path.basename(data_file))[0])
        if subdir_of:
            name = _slug(model_key)
            directory = os.path.join(subdir_of, name)
        else:
            name = f"{ts}_{dataname}_{_slug(model_key)}_{_slug(backend)}"
            directory = os.path.join(self.root, name)
        os.makedirs(directory, exist_ok=True)
        sess = Session(directory)
        sess.state = {
            "name": name,
            "created": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "data_file": data_file,
            "material_description": config.get("material_description", ""),
            "model": model_key,
            "backend": backend,
            "optimiser": optimiser,
            "config": config,
            "bounds": bounds,
            "stage": "idle",
            "prior_runs": [],
        }
        sess.save()
        return sess

    def open(self, directory: str) -> Session:
        """Open an existing session directory."""
        return Session(directory)

    def list_sessions(self) -> list[SessionInfo]:
        """All sessions under the root, newest first."""
        infos: list[SessionInfo] = []
        for entry in sorted(os.listdir(self.root), reverse=True):
            d = os.path.join(self.root, entry)
            if os.path.isdir(d) and os.path.exists(os.path.join(d, SESSION_JSON)):
                try:
                    infos.append(Session(d).info())
                except Exception:
                    continue
        return infos

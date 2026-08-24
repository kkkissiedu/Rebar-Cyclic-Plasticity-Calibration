"""
ui/progress_panel.py
====================

Independent progress display for the three pipeline stages:

* Stage 1 - Sobol seeding      (X of Y evaluations, %, ETA)
* Stage 2 - DE / Bayesian      (iteration, current best objective, ETA)
* Stage 3 - ABAQUS verify      (candidate X of M, elapsed)

Each stage shows an animated bar and a status badge:
PENDING / RUNNING / DONE / FAILED.
"""

from __future__ import annotations

from tkinter import ttk
from typing import Optional

PENDING, RUNNING, DONE, FAILED = "PENDING", "RUNNING", "DONE", "FAILED"
_BADGE_STYLE = {PENDING: "Muted.TLabel", RUNNING: "Accent.TLabel",
                DONE: "Ok.TLabel", FAILED: "Fail.TLabel"}


class _Stage:
    """One stage row: title, badge, bar, detail line."""

    def __init__(self, parent: ttk.Frame, row: int, title: str) -> None:
        ttk.Label(parent, text=title).grid(row=row, column=0, sticky="w",
                                           padx=4, pady=(4, 0))
        self.badge = ttk.Label(parent, text=PENDING, style=_BADGE_STYLE[PENDING])
        self.badge.grid(row=row, column=1, sticky="e", padx=4, pady=(4, 0))
        self.bar = ttk.Progressbar(parent, mode="determinate", maximum=100)
        self.bar.grid(row=row + 1, column=0, columnspan=2, sticky="we", padx=4)
        self.detail = ttk.Label(parent, text="-", style="Muted.TLabel")
        self.detail.grid(row=row + 2, column=0, columnspan=2, sticky="w",
                         padx=4, pady=(0, 6))
        self._animating = False

    def set(self, status: str, done: int = 0, total: int = 0,
            detail: str = "") -> None:
        """Update this stage's badge, bar and detail line."""
        self.badge.configure(text=status, style=_BADGE_STYLE.get(status, "Muted.TLabel"))
        if status == RUNNING and total <= 0:
            # unknown length -> indeterminate marquee
            if not self._animating:
                self.bar.configure(mode="indeterminate")
                self.bar.start(60)
                self._animating = True
        else:
            if self._animating:
                self.bar.stop()
                self.bar.configure(mode="determinate")
                self._animating = False
            pct = (100.0 * done / total) if total > 0 else (100.0 if status == DONE else 0.0)
            self.bar.configure(value=max(0.0, min(100.0, pct)))
        if detail:
            self.detail.configure(text=detail)


def _eta(done: int, total: int, elapsed: float) -> str:
    """Human ETA string from progress + elapsed seconds."""
    if done <= 0 or total <= 0 or elapsed <= 0:
        return "ETA --"
    rate = done / elapsed
    remain = (total - done) / rate if rate > 0 else 0
    return f"ETA {remain:5.0f}s"


class ProgressPanel:
    """Container holding the three stage rows."""

    def __init__(self, parent) -> None:
        self.frame = ttk.LabelFrame(parent, text="Progress")
        self.frame.columnconfigure(0, weight=1)
        self.stage1 = _Stage(self.frame, 0, "Stage 1 - Sobol seeding")
        self.stage2 = _Stage(self.frame, 3, "Stage 2 - Refinement (DE / Bayesian)")
        self.stage3 = _Stage(self.frame, 6, "Stage 3 - ABAQUS verification")
        self._stages = {"sobol": self.stage1, "de": self.stage2,
                        "bayesian": self.stage2, "abaqus": self.stage3}

    def update_stage(self, stage_key: str, *, status: str, done: int = 0,
                     total: int = 0, elapsed: float = 0.0,
                     best: Optional[float] = None) -> None:
        """Update one stage from a pipeline progress message."""
        stage = self._stages.get(stage_key)
        if stage is None:
            return
        parts = []
        if total > 0:
            parts.append(f"{done}/{total} ({100.0*done/max(total,1):.0f}%)")
        elif done > 0:
            parts.append(f"{done} evals")
        if best is not None and best < 1e8:
            parts.append(f"best {best:.2f} MPa")
        parts.append(_eta(done, total, elapsed))
        stage.set(status, done, total, "  ".join(parts))

    def reset(self) -> None:
        for s in (self.stage1, self.stage2, self.stage3):
            s.set(PENDING, 0, 0, "-")

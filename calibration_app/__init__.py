"""
calibration_app
================

Chaboche / UVC / Ohno-Wang cyclic-plasticity calibration application.

Rebuilt package (Phase 4). Layout:

    core/     model surrogates, objective, search, ABAQUS runner, data loader
    utils/    grade-agnostic physics gate
    session/  checkpoint + session management
    ui/       tkinter GUI

The B500C / Kashani 2019 dataset is used for pipeline validation only; the
intended target material is locally manufactured (Ghanaian scrap-metal) rebar
with batch-variable properties, so the whole stack is grade-agnostic: yield
stress and bounds are per-run inputs, and the physics gate checks consistency
against the loaded data rather than any fixed steel grade.
"""

__version__ = "0.2.0"

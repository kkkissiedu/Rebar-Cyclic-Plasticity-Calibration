"""
core/cae_transfer.py
====================

"Transfer to CAE": push a completed session's calibrated material, loading
amplitude, step controls and output requests into a SESSION-SPECIFIC COPY of
FATIGUE_rebuilt.cae (``FATIGUE_{session}_{D}mm_LD{ratio}.cae``, same folder).
Neither FATIGUE.cae nor FATIGUE_rebuilt.cae is ever modified. When the
session geometry differs from the base 12mm/LD5, the copy's cylinder is
rebuilt (new circle sketch + extrusion, re-seeded, remeshed C3D8R, quality-
checked but never aborted, section reassigned, BC regions rebound, SS_AUTO /
TOP_AUTO recreated).

Flow (all driven from the GUI, but usable headless):

1. :func:`build_transfer_config` -- collect everything the CAE script needs
   from a session's ``session.json`` state (best params, sigma_y0, geometry,
   N_fit) plus the user-configured .cae path. All three models supported
   (2026-07-14): Chaboche writes the native combined-hardening tables;
   UVC / Ohno-Wang write *USER MATERIAL + *DEPVAR, switch elements to C3D8
   (user materials have no hourglass stiffness for C3D8R) and prebuild the
   multiaxial C++ UMAT into standardU.dll (:func:`build_umat_library`,
   loaded via ``usub_lib_dir`` -- submit needs no compiler). Refuses any
   attempt to target FATIGUE.cae (hard rule: never modified).
2. :func:`write_transfer_script` -- write ``transfer_to_cae.py`` into the
   session folder with the config baked in as JSON. The script is idempotent:
   every operation *sets* state (tables, amplitude data, step controls,
   output requests), so running it twice yields the same model.
3. :func:`run_transfer_script` -- ``abaqus cae noGUI=transfer_to_cae.py``
   (shell=True on Windows, confirmed requirement). Success is detected by the
   ``TRANSFER_OK`` marker on stdout, not just the exit code (abaqus cae can
   exit 0 after a scripting error).
4. :func:`write_transfer_log` -- ``transfer_log.json`` in the session folder
   recording timestamp, session, model, every value transferred and the .cae
   path. ``plot_fe_comparison.py`` reads this log after the job runs.

What the generated script updates in the copy (spec 2026-07-08 + 07-11):
geometry rebuild when the session is not 12mm/LD5 (see above), then material
(*PLASTIC HARDENING=COMBINED + *CYCLIC HARDENING tables), first
TabularAmplitude (triangular wave, period 20 s, N_fit cycles), the
Cyclic_Loading step (timePeriod / increments / nlgeom=ON), the displacement
BC's amplitude reference, and fresh history (SS_AUTO: S33/NE33/LE33,
TOP_AUTO: RF3/U3, freq=1) + field (S,E,PE,PEEQ,U,RF, freq=10) output requests.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import time

#: Kashani test frequency -- fixed loading period for the full-model run (s).
PERIOD_S: float = 20.0

#: Name of the cyclic step expected in FATIGUE_rebuilt.cae.
STEP_NAME: str = "Cyclic_Loading"

TRANSFER_SCRIPT = "transfer_to_cae.py"
TRANSFER_LOG = "transfer_log.json"

#: Geometry FATIGUE_rebuilt.cae was meshed for (warning only, never changed).
CAE_DIAMETER_MM: float = 12.0
CAE_LD_RATIO: float = 5.0


class TransferError(RuntimeError):
    """A transfer cannot proceed (wrong model, bad target, missing data)."""


def build_transfer_config(state: dict, cae_path: str,
                          session_dir: str = None) -> dict:
    """Assemble the transfer config from a session's saved state.

    Only native-Chaboche sessions are fully automated. UVC / Ohno-Wang use a
    UMAT ``*USER MATERIAL`` card that the CAE material editor cannot express
    the same way -- rather than silently writing wrong values, refuse.

    ``session_dir``: when given, the working copy goes into
    ``<session_dir>/cae/`` so the session folder holds EVERYTHING for the
    run (calibration + cae + job outputs + comparison) -- submit the ABAQUS
    job with the work directory set to that ``cae/`` folder. Without it the
    copy lands next to the base template (legacy behaviour).
    """
    model_key = state.get("model")
    if model_key not in ("chaboche", "uvc", "ohno_wang"):
        raise TransferError(f"Unknown model '{model_key}'.")

    base = os.path.basename(cae_path).lower()
    if base == "fatigue.cae":
        raise TransferError("FATIGUE.cae is the protected original and is "
                            "never modified. Point the setting at "
                            "FATIGUE_rebuilt.cae.")
    if not os.path.exists(cae_path):
        raise TransferError(f"CAE file not found: {cae_path}")

    best = state.get("best_params")
    if not best:
        raise TransferError("Session has no best parameters yet.")
    cfg = state.get("config") or {}
    n_bs = int(cfg.get("n_backstresses", 3))
    expected = {"chaboche": 2 * n_bs + 2,     # C1,g1,...,Q,b
                "uvc": 4 + 2 * n_bs,          # Q,b,D,a,C1,g1,...
                "ohno_wang": 2 + 3 * n_bs}    # Q,b,C1,r1,m1,...
    if len(best) != expected[model_key]:
        raise TransferError(
            f"best_params length {len(best)} does not match "
            f"{n_bs} backstresses for {model_key} "
            f"(expected {expected[model_key]}).")

    sy0 = state.get("best_sigma_y0") or cfg.get("sigma_y0")
    if not sy0:
        raise TransferError("No calibrated sigma_y0 in the session.")

    E_MOD, NU = 200000.0, 0.3
    best_f = [float(v) for v in best]
    if model_key == "chaboche":
        # Param vector: (C1, g1, ..., CN, gN, Q, b) -- surrogate.py order.
        material_mode = "native"
        cg = best_f[: 2 * n_bs]
        Q, b = best_f[2 * n_bs], best_f[2 * n_bs + 1]
        plastic_table = [float(sy0)] + cg                # sy0, C1, g1, ...
        cyclic_table = [float(sy0), Q, b]
        user_props, nstatv, umat_src = None, None, None
    else:
        # UMAT models: the CAE copy gets *USER MATERIAL + *DEPVAR and the
        # multiaxial C++ UMAT (umats/*multiaxial*.cpp) prebuilt into a
        # standardU.dll loaded via usub_lib_dir (see build_umat_library).
        # PROPS layouts match the UMAT headers (E, nu leading — the 3D
        # subroutines need Poisson's ratio, unlike the uniaxial FE cards).
        material_mode = "user"
        plastic_table, cyclic_table = None, None
        if model_key == "uvc":
            user_props = [E_MOD, NU, float(sy0)] + best_f       # Q,b,D,a,C,g
        else:
            user_props = [E_MOD, NU, float(sy0)] + best_f       # Q,b,C,r,m
        nstatv = 7 + 6 * n_bs
        proj_root = os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))))
        umat_rel = ("calibration_app/umats/UVCmultiaxial.cpp"
                    if model_key == "uvc"
                    else "calibration_app/umats/OhnoWangMultiaxial.cpp")
        umat_src = os.path.join(proj_root, umat_rel)
        if not os.path.exists(umat_src):
            raise TransferError(f"Multiaxial UMAT missing: {umat_src}")

    strain_pct = float(cfg.get("strain_pct"))
    gauge = float(cfg.get("gauge_length_mm"))
    diameter = float(cfg.get("diameter_mm", 0.0))
    n_fit = int(cfg.get("n_fit"))
    amp_mm = strain_pct / 100.0 * gauge

    # Session-specific working copy: the base template is NEVER modified.
    # Geometry differing from the base (12mm/LD5) triggers a rebuild of the
    # cylinder + mesh inside the copy.
    if diameter <= 0:
        raise TransferError("Session has no bar diameter - cannot size the "
                            "working copy geometry.")
    ld_ratio = gauge / diameter
    rebuild = (abs(diameter - CAE_DIAMETER_MM) > 1e-6
               or abs(ld_ratio - CAE_LD_RATIO) > 1e-6)
    safe_session = re.sub(r"[^A-Za-z0-9._-]+", "-",
                          str(state.get("name", "run"))).strip("-") or "run"
    copy_name = f"FATIGUE_{safe_session}_{diameter:g}mm_LD{ld_ratio:g}.cae"
    if session_dir:
        cae_dir = os.path.join(os.path.abspath(session_dir), "cae")
        os.makedirs(cae_dir, exist_ok=True)
    else:
        cae_dir = os.path.dirname(os.path.abspath(cae_path))
    copy_path = os.path.join(cae_dir, copy_name)

    # Cycle-1 burn-in path (machine ramp) prepended to the loading amplitude
    # so the full model reaches cycle 2 with conditioned hardening state --
    # the same treatment Stage 3 gives the single element. Without it the
    # virgin first FE cycle is scored against the already-cycled experimental
    # cycle 2 (observed: cycle-2 RMSE 114 vs ~65 MPa steady-state, 3data).
    burnin_table: list[list[float]] = []
    data_file = state.get("data_file")
    if data_file and os.path.exists(data_file):
        try:
            from . import data_loader
            d = data_loader.load_experimental(
                data_file, gauge_length_mm=gauge, bar_diameter_mm=diameter)
            bs = d.get("burnin_strain")
            if bs is not None and len(bs) >= 4:
                stride = max(1, len(bs) // 48)
                samples = [float(v) for v in bs[::stride]]
                if samples[-1] != float(bs[-1]):
                    samples.append(float(bs[-1]))
                n_pts = len(samples)
                for i, sv in enumerate(samples):
                    burnin_table.append(
                        [PERIOD_S * i / float(n_pts), sv * gauge])
        except Exception:
            burnin_table = []          # burn-in is an enhancement, not a gate

    return {
        "session_name": state.get("name", "?"),
        "model": model_key,
        "burnin_table": burnin_table,
        "burnin_cycles": 1 if burnin_table else 0,
        # cae_path = the working COPY (what the script edits and the job
        # runs from); base_cae_path = the untouched template.
        "cae_path": copy_path,
        "base_cae_path": os.path.abspath(cae_path),
        "rebuild_geometry": rebuild,
        "ld_ratio": ld_ratio,
        # Paper-grade density (2026-07-11, replaces the original coarse
        # max(3,D/4)/max(10,gauge/3) formulas): ~24+ elements around the
        # circumference, ~2 elements/mm axially, global size D/10 driving
        # the radial layers. 16mm/LD5 -> a few thousand C3D8R.
        "circum_seeds": max(24, int(round(1.5 * diameter))),
        "axial_seeds": max(40, int(round(gauge / 2.0))),
        "global_seed_mm": diameter / 10.0,
        "n_backstresses": n_bs,
        "sigma_y0_cyclic": float(sy0),
        "material_mode": material_mode,              # "native" | "user"
        "plastic_table": plastic_table,              # native only
        "cyclic_hardening_table": cyclic_table,      # native only
        "user_props": user_props,                    # user only
        "nstatv": nstatv,                            # user only (*DEPVAR)
        "umat_source": umat_src,                     # user only (abs .cpp)
        "E": E_MOD,
        "nu": NU,
        "strain_pct": strain_pct,
        "gauge_length_mm": gauge,
        "diameter_mm": diameter,
        "amp_mm": amp_mm,
        "period_s": PERIOD_S,
        "n_cycles": n_fit,
        "step_name": STEP_NAME,
        "weight_reversals": bool(cfg.get("weight_reversals", True)),
        "surrogate_objective_MPa": state.get("best_objective"),
        "surrogate_peak_tension_MPa": state.get("best_peak_tension"),
        "surrogate_peak_compression_MPa": state.get("best_peak_compression"),
        "data_file": state.get("data_file"),
    }


def parse_mesh_stats(output: str) -> dict:
    """Pull the rebuild markers the generated script prints from its stdout."""
    stats = {"n_elements": None, "n_failed": None}
    m = re.search(r"N_ELEMENTS=(\d+)", output or "")
    if m:
        stats["n_elements"] = int(m.group(1))
    m = re.search(r"MESH_FAILED=(\d+)", output or "")
    if m:
        stats["n_failed"] = int(m.group(1))
    return stats


def transfer_summary(tcfg: dict, output: str = "") -> str:
    """Dialog text after a successful transfer (spec 2026-07-11 wording)."""
    stats = parse_mesh_stats(output)
    lines = [
        f"Working on copy: {os.path.basename(tcfg['cae_path'])}",
        f"Base template unchanged: {os.path.basename(tcfg['base_cae_path'])}",
        f"Geometry: {tcfg['diameter_mm']:g}mm dia x "
        f"{tcfg['gauge_length_mm']:g}mm gauge (LD{tcfg['ld_ratio']:g})",
    ]
    if tcfg.get("rebuild_geometry"):
        n = stats["n_elements"]
        etype = "C3D8" if tcfg.get("material_mode") == "user" else "C3D8R"
        lines.append(f"Mesh rebuilt: {n if n is not None else '?'} "
                     f"{etype} elements")
        if stats["n_failed"]:
            lines.append(f"Mesh quality: {stats['n_failed']} element(s) "
                         f"failed analysis checks - review in CAE before "
                         f"submitting (transfer NOT aborted)")
    else:
        lines.append("Geometry: matches base (no rebuild needed)")
    if tcfg.get("burnin_cycles"):
        lines.append("Burn-in: experimental cycle-1 ramp prepended "
                     "(1 extra cycle before the scored cycles)")
    if tcfg.get("material_mode") == "user":
        lines.append("Material: *USER MATERIAL (%d constants, DEPVAR=%d) "
                     "from %s" % (len(tcfg["user_props"]), tcfg["nstatv"],
                                  os.path.basename(tcfg["umat_source"])))
        lines.append("Elements switched to C3D8 (full integration - user "
                     "materials have no hourglass stiffness for C3D8R)")
        lines.append("UMAT precompiled to standardU.dll in cae/ "
                     "(usub_lib_dir in abaqus_v6.env) - submit the job "
                     "NORMALLY: no compiler and no user-subroutine job "
                     "setting needed")
    lines.append("Submit the job with ABAQUS work directory set to: "
                 + os.path.dirname(tcfg["cae_path"]))
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# The CAE-side script. Runs under `abaqus cae noGUI=` (ABAQUS 2024 kernel
# Python). Idempotent: every mutation is a set-to-value. Never touches
# geometry or mesh. Emits TRANSFER_OK / TRANSFER_FAIL markers for the caller.
# ---------------------------------------------------------------------------
_SCRIPT_TEMPLATE = r'''# -*- coding: utf-8 -*-
# transfer_to_cae.py -- generated by the calibration app ("Transfer to CAE").
# Run with:  abaqus cae noGUI=transfer_to_cae.py
# Idempotent: safe to run repeatedly; every operation sets absolute values.
# NEVER points at FATIGUE.cae (guarded below and at generation time).
import json
import os
import shutil
import sys
import traceback

CONFIG = json.loads(r"""__CONFIG_JSON__""")

# On Windows the `abaqus cae noGUI` launcher runs the kernel in a child
# process whose stdout is NOT piped back to the caller (only the license
# banner is). So every message is also accumulated and flushed to
# transfer_to_cae_output.log in the working dir -- the app reads THAT file
# for the TRANSFER_OK/TRANSFER_FAIL marker and the mesh stats.
LOG_LINES = []
OUTPUT_LOG = "transfer_to_cae_output.log"


def log(msg):
    msg = str(msg)
    print(msg)
    LOG_LINES.append(msg)


def _write_result():
    try:
        f = open(os.path.join(os.getcwd(), OUTPUT_LOG), "w")
        try:
            f.write("\n".join(LOG_LINES) + "\n")
        finally:
            f.close()
    except Exception:
        pass


def fail(msg):
    log("TRANSFER_FAIL: %s" % msg)
    sys.exit(1)


def main():
    from abaqus import openMdb, mdb  # noqa: F401 (openMdb rebinds mdb)
    from abaqusConstants import COMBINED, PARAMETERS, ON, UNSET
    # Preload the CAE kernel feature modules: some model repositories (e.g.
    # historyOutputRequests) only materialise once their module is imported.
    for _mod in ("part", "material", "section", "assembly", "step",
                 "interaction", "load", "mesh", "job", "sketch",
                 "regionToolset"):
        try:
            __import__(_mod)
        except Exception:
            pass

    base_path = CONFIG["base_cae_path"]
    copy_path = CONFIG["cae_path"]
    for p in (base_path, copy_path):
        if os.path.basename(p).lower() == "fatigue.cae":
            fail("refusing to touch protected FATIGUE.cae")
    if not os.path.exists(base_path):
        fail("base CAE file not found: %s" % base_path)
    # ALWAYS work on a fresh session-specific copy: the base template is
    # never opened for writing, and re-running re-copies (idempotent).
    shutil.copyfile(base_path, copy_path)
    log("COPY: %s" % os.path.basename(copy_path))
    log("BASE_UNCHANGED: %s" % os.path.basename(base_path))

    my_mdb = openMdb(pathName=copy_path)

    # --- locate the model + combined-hardening material -------------------
    target_model = None
    target_mat = None
    for mname in my_mdb.models.keys():
        model = my_mdb.models[mname]
        for matname in model.materials.keys():
            mat = model.materials[matname]
            plastic = getattr(mat, "plastic", None)
            if plastic is not None and plastic.hardening == COMBINED:
                target_model, target_mat = model, mat
                break
        if target_mat is not None:
            break
    if target_mat is None:
        fail("no material with *PLASTIC HARDENING=COMBINED found in %s"
             % cae_path)
    log("Material: %s (model %s)" % (target_mat.name, target_model.name))

    # --- 0. geometry: rebuild the copy's cylinder if the session differs ---
    if CONFIG["rebuild_geometry"]:
        from abaqusConstants import (C3D8R, C3D8, STRUCTURED, HEX,
                                     ANALYSIS_CHECKS)
        import mesh as mesh_mod
        # user materials cannot use C3D8R (no hourglass stiffness)
        elem_code = C3D8 if CONFIG.get("material_mode") == "user" else C3D8R
        import regionToolset
        D = float(CONFIG["diameter_mm"])
        gauge = float(CONFIG["gauge_length_mm"])
        pnames = list(target_model.parts.keys())
        if not pnames:
            fail("model has no parts to rebuild")
        pname = "Rebar_Coupon" if "Rebar_Coupon" in pnames else pnames[0]
        part = target_model.parts[pname]
        # Remember the section so the new mesh can be re-covered afterwards.
        sect_name = None
        if len(part.sectionAssignments) > 0:
            sect_name = part.sectionAssignments[0].sectionName
        elif len(target_model.sections.keys()) > 0:
            sect_name = list(target_model.sections.keys())[0]
        # Delete the base extrusion (with its sketch and every dependent
        # feature), then recreate the circle + extrusion at session geometry.
        part.deleteFeatures(tuple(part.features.keys()))
        sk = target_model.ConstrainedSketch(name="rebuild_profile",
                                            sheetSize=max(4.0 * D, 20.0))
        sk.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(D / 2.0, 0.0))
        part.BaseSolidExtrude(sketch=sk, depth=gauge)
        # Seeding: curved (circumferential) edges by the circum count,
        # straight (axial) edges by the axial count; the global seed size
        # (D/10) drives the through-radius layer count.
        circum = int(CONFIG["circum_seeds"])
        axial = int(CONFIG["axial_seeds"])
        part.seedPart(size=float(CONFIG["global_seed_mm"]),
                      deviationFactor=0.1, minSizeFactor=0.1)
        for i in range(len(part.edges)):
            try:
                part.edges[i].getRadius()
                n = circum
            except Exception:
                n = axial
            part.seedEdgeByNumber(edges=part.edges[i:i + 1], number=n)
        try:
            part.setMeshControls(regions=part.cells, technique=STRUCTURED,
                                 elemShape=HEX)
        except Exception as exc:
            log("WARNING: structured hex not applicable (%s); using the "
                  "default meshing technique" % exc)
        part.setElementType(regions=(part.cells,),
                            elemTypes=(mesh_mod.ElemType(elemCode=elem_code),))
        part.generateMesh()
        n_elements = len(part.elements)
        if n_elements == 0:
            fail("remesh produced 0 elements")
        log("N_ELEMENTS=%d" % n_elements)
        # Quality check: report failures but NEVER abort (user decides).
        try:
            mq = part.verifyMeshQuality(criterion=ANALYSIS_CHECKS)
            failed = mq.get("failedElements", ())
            log("MESH_FAILED=%d C3D8R" % len(failed))
        except Exception as exc:
            log("WARNING: verifyMeshQuality failed: %s" % exc)
        # Section assignment must cover the new mesh.
        if sect_name is not None:
            try:
                for _ in range(len(part.sectionAssignments)):
                    del part.sectionAssignments[0]
            except Exception:
                pass
            part.SectionAssignment(
                region=regionToolset.Region(cells=part.cells),
                sectionName=sect_name)
            log("Section '%s' reassigned to the new mesh" % sect_name)
        else:
            log("WARNING: no section found to reassign - assign manually")
        target_model.rootAssembly.regenerate()
        log("Geometry rebuilt: %gmm dia x %gmm gauge (LD%g)"
              % (D, gauge, CONFIG["ld_ratio"]))
    else:
        log("GEOMETRY_MATCHES_BASE")

    # --- 1. material parameters (calibrated, incl. cyclic sigma_y0) -------
    n_bs = int(CONFIG["n_backstresses"])
    if CONFIG.get("material_mode") == "user":
        # UMAT models (UVC / Ohno-Wang): replace the native behaviours with
        # *USER MATERIAL + *DEPVAR. The UMAT supplies elasticity too, so the
        # Elastic/Plastic definitions must go (leaving them in makes the
        # input processor reject the material as doubly defined).
        for beh in ("plastic", "elastic", "cyclicHardening"):
            try:
                exec("del target_mat.%s" % beh)
                log("Removed native behaviour: %s" % beh)
            except Exception:
                pass
        target_mat.Depvar(n=int(CONFIG["nstatv"]))
        target_mat.UserMaterial(
            mechanicalConstants=tuple(CONFIG["user_props"]))
        log("USER MATERIAL: %d constants, DEPVAR=%d (UMAT %s)"
            % (len(CONFIG["user_props"]), int(CONFIG["nstatv"]),
               os.path.basename(CONFIG.get("umat_source") or "?")))
        # C3D8R + user material has no hourglass stiffness (fatal input
        # error) -> switch every meshed part to fully-integrated C3D8.
        from abaqusConstants import C3D8 as _C3D8
        import mesh as _mesh_mod
        for _pn in target_model.parts.keys():
            _pt = target_model.parts[_pn]
            try:
                if len(_pt.elements) > 0:
                    _pt.setElementType(
                        regions=(_pt.cells,),
                        elemTypes=(_mesh_mod.ElemType(elemCode=_C3D8),))
                    log("ELEM_TYPE=C3D8 on part '%s' (user material has no "
                        "hourglass stiffness for C3D8R)" % _pn)
            except Exception as exc:
                log("WARNING: could not switch part '%s' to C3D8: %s - set "
                    "the element type manually before submitting" % (_pn, exc))
        target_model.rootAssembly.regenerate()
    else:
        target_mat.Plastic(hardening=COMBINED, dataType=PARAMETERS,
                           numBackstresses=n_bs,
                           table=(tuple(CONFIG["plastic_table"]),))
        target_mat.plastic.CyclicHardening(
            parameters=ON, table=(tuple(CONFIG["cyclic_hardening_table"]),))
        log("Plastic table: %s" % (CONFIG["plastic_table"],))
        log("Cyclic hardening table: %s" % (CONFIG["cyclic_hardening_table"],))

    # --- 2. loading amplitude (triangular, matches the test file) ---------
    # If a burn-in table is present, the experimental cycle-1 machine ramp
    # is prepended over [0, T) so the scored cycles start with conditioned
    # hardening state (mirrors Stage 3's burn-in handling).
    amp = float(CONFIG["amp_mm"])
    T = float(CONFIG["period_s"])
    n_cyc = int(CONFIG["n_cycles"])
    burnin = CONFIG.get("burnin_table") or []
    pts = []
    t_off = 0.0
    if burnin:
        for tb, vb in burnin:
            pts.append((float(tb), float(vb)))
        t_off = T
        log("Burn-in ramp prepended: %d points over %.0f s" % (len(burnin), T))
    for c in range(n_cyc):
        t0 = t_off + c * T
        pts.append((t0, 0.0))
        pts.append((t0 + T / 4.0, amp))
        pts.append((t0 + 3.0 * T / 4.0, -amp))
    pts.append((t_off + n_cyc * T, 0.0))
    total_time = t_off + n_cyc * T

    amp_name = None
    for aname in target_model.amplitudes.keys():
        a = target_model.amplitudes[aname]
        if type(a).__name__ == "TabularAmplitude":
            a.setValues(data=tuple(pts))
            amp_name = aname
            break
    if amp_name is None:
        # No tabular amplitude in the model yet: create one.
        amp_name = "AMP_TRI"
        target_model.TabularAmplitude(name=amp_name, data=tuple(pts))
    log("Amplitude '%s': %d points, peak %.4f mm, %d cycles"
          % (amp_name, len(pts), amp, n_cyc))

    # --- 3. step time period + increments ---------------------------------
    step_name = CONFIG["step_name"]
    if step_name not in target_model.steps.keys():
        # Fall back to the first non-Initial step.
        others = [s for s in target_model.steps.keys() if s != "Initial"]
        if not others:
            fail("no analysis step found (expected '%s')" % step_name)
        step_name = others[0]
        log("Step '%s' not found; using '%s'" % (CONFIG["step_name"],
                                                   step_name))
    step = target_model.steps[step_name]
    step.setValues(timePeriod=total_time,
                   initialInc=T / 200.0,
                   minInc=T / 1.0e6,
                   maxInc=T / 50.0,
                   maxNumInc=100000,
                   nlgeom=ON)
    log("Step '%s': timePeriod=%.1f s (%d scored + %d burn-in cycles), "
        "nlgeom=ON" % (step_name, total_time, n_cyc,
                       1 if burnin else 0))

    # --- 4. displacement BC references the amplitude ----------------------
    # A BC's amplitude (and prescribed value) is STEP-DEPENDENT state, not an
    # attribute of the BC object itself -- read it from the step's
    # boundaryConditionStates. "Driven" = has an amplitude in this step OR a
    # nonzero prescribed u3 (the fixed bottom BC has u3 = 0).
    def bc_is_driven(bcname):
        try:
            st = target_model.steps[step_name].boundaryConditionStates[bcname]
        except Exception:
            return False
        ampv = getattr(st, "amplitude", UNSET)
        if ampv not in (UNSET, "UNSET", "", None):
            return True
        try:
            return abs(float(getattr(st, "u3", 0.0))) > 0.0
        except (TypeError, ValueError):
            return False

    bc_driven = []
    for bcname in target_model.boundaryConditions.keys():
        bc = target_model.boundaryConditions[bcname]
        if type(bc).__name__ != "DisplacementBC":
            continue
        if bc_is_driven(bcname):
            bc.setValuesInStep(stepName=step_name, amplitude=amp_name)
            bc_driven.append(bcname)
    if not bc_driven:
        log("WARNING: no driven DisplacementBC found in step '%s'; "
            "check the loading BC manually." % step_name)
    else:
        log("Displacement BC(s) using '%s': %s" % (amp_name, bc_driven))

    # After a geometry rebuild the BC regions reference deleted geometry:
    # rebind driven BCs (amplitude set) to the new top face and fixed BCs to
    # the bottom face. Best effort -- failures are reported, never fatal.
    if CONFIG["rebuild_geometry"]:
        import regionToolset
        a0 = target_model.rootAssembly
        inames = list(a0.instances.keys())
        if inames:
            inst0 = a0.instances[inames[0]]
            gauge = float(CONFIG["gauge_length_mm"])
            tol = 1e-4 * max(1.0, gauge)
            top_faces = inst0.faces.getByBoundingBox(
                xMin=-1e9, yMin=-1e9, zMin=gauge - tol,
                xMax=1e9, yMax=1e9, zMax=gauge + tol)
            bot_faces = inst0.faces.getByBoundingBox(
                xMin=-1e9, yMin=-1e9, zMin=-tol,
                xMax=1e9, yMax=1e9, zMax=tol)
            # Rebind EVERY BC type -- a fixed/encastre BC left pointing at a
            # wiped part set makes the input writer emit an undefined set
            # (seen: "Unknown part instance node set SPECIMEN-1.SET-BOTTOM"
            # -> Analysis Input File Processor abort).
            for bcname in target_model.boundaryConditions.keys():
                bc = target_model.boundaryConditions[bcname]
                driven = bc_is_driven(bcname)
                faces = top_faces if driven else bot_faces
                if len(faces) == 0:
                    continue
                try:
                    bc.setValues(region=regionToolset.Region(faces=faces))
                    log("BC '%s' (%s) region rebound to %s face"
                        % (bcname, type(bc).__name__,
                           "top" if driven else "bottom"))
                except Exception as exc:
                    log("WARNING: could not rebind BC '%s' (%s) region: %s - "
                        "fix it in CAE before submitting, the old region no "
                        "longer exists" % (bcname, type(bc).__name__, exc))

    # --- 5. SS_AUTO / TOP_AUTO sets (created only if absent) --------------
    a = target_model.rootAssembly
    inst_names = list(a.instances.keys())
    if not inst_names:
        fail("assembly has no instances")
    inst = a.instances[inst_names[0]]
    gauge_mid = float(CONFIG["gauge_length_mm"]) / 2.0

    # A rebuilt mesh invalidates any pre-existing sets: recreate them fresh.
    if CONFIG["rebuild_geometry"]:
        for sname in ("SS_AUTO", "TOP_AUTO"):
            if sname in a.sets.keys():
                del a.sets[sname]

    if "TOP_AUTO" not in a.sets.keys():
        zmax = max(n.coordinates[2] for n in inst.nodes)
        tol = 1e-6 * max(1.0, abs(zmax))
        top_nodes = inst.nodes.getByBoundingBox(
            xMin=-1e9, yMin=-1e9, zMin=zmax - tol,
            xMax=1e9, yMax=1e9, zMax=zmax + tol)
        if len(top_nodes) == 0:
            fail("could not find top-face nodes (z = %.4f)" % zmax)
        a.Set(name="TOP_AUTO", nodes=top_nodes)
        log("Created TOP_AUTO: %d nodes at z=%.4f" % (len(top_nodes), zmax))

    if "SS_AUTO" not in a.sets.keys():
        node_z = {}
        for n in inst.nodes:
            node_z[n.label] = n.coordinates[2]
        best_lab, best_d = None, None
        for el in inst.elements:
            zs = [node_z[lab] for lab in el.connectivity if lab in node_z]
            if not zs:
                continue
            cz = sum(zs) / float(len(zs))
            d = abs(cz - gauge_mid)
            if best_d is None or d < best_d:
                best_lab, best_d = el.label, d
        if best_lab is None:
            fail("could not locate a mid-specimen element for SS_AUTO")
        a.Set(name="SS_AUTO",
              elements=inst.elements.sequenceFromLabels((best_lab,)))
        log("Created SS_AUTO: element %d (centroid %.3f mm from z=%.3f)"
              % (best_lab, best_d, gauge_mid))

    # --- 6. output requests (fresh, comparison-optimised) ------------------
    # Guarded: some kernel sessions / imported models do not expose the
    # output-request repositories -- report clearly instead of dying after
    # the material/amplitude/step work is already done.
    try:
        for name in list(target_model.historyOutputRequests.keys()):
            del target_model.historyOutputRequests[name]
        for name in list(target_model.fieldOutputRequests.keys()):
            del target_model.fieldOutputRequests[name]
        target_model.HistoryOutputRequest(
            name="H_SS_AUTO", createStepName=step_name,
            region=a.sets["SS_AUTO"],
            variables=("S33", "NE33", "LE33"), frequency=1)
        target_model.HistoryOutputRequest(
            name="H_TOP_AUTO", createStepName=step_name,
            region=a.sets["TOP_AUTO"],
            variables=("RF3", "U3"), frequency=1)
        target_model.FieldOutputRequest(
            name="F_COMPARISON", createStepName=step_name,
            variables=("S", "E", "PE", "PEEQ", "U", "RF"), frequency=10)
        log("Output requests rebuilt: H_SS_AUTO (S33/NE33/LE33, freq 1), "
            "H_TOP_AUTO (RF3/U3, freq 1), F_COMPARISON (freq 10)")
    except AttributeError as exc:
        log("WARNING: output-request repositories unavailable in this "
            "kernel session (%s). Add them manually in CAE: history freq=1 "
            "S33/NE33/LE33 on set SS_AUTO + RF3/U3 on TOP_AUTO; field "
            "freq=10 S,E,PE,PEEQ,U,RF. Sets already exist." % exc)

    my_mdb.save()
    log("Saved %s" % copy_path)
    log("TRANSFER_OK")


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        _write_result()
        raise
    except Exception:
        log("TRANSFER_FAIL: unhandled exception")
        log(traceback.format_exc())
        _write_result()
        sys.exit(1)
    else:
        _write_result()
'''


def write_transfer_script(session_dir: str, tcfg: dict) -> str:
    """Write ``transfer_to_cae.py`` (config baked in) into the session dir."""
    script = _SCRIPT_TEMPLATE.replace(
        "__CONFIG_JSON__", json.dumps(tcfg, indent=1))
    path = os.path.join(session_dir, TRANSFER_SCRIPT)
    with open(path, "w") as f:
        f.write(script)
    return path


#: Kernel-side log written by the generated script (see template note: on
#: Windows the CAE kernel's stdout is not piped back to the caller).
OUTPUT_LOG = "transfer_to_cae_output.log"


def run_transfer_script(script_path: str, abaqus_cmd: str = "abaqus",
                        timeout_s: int = 900) -> tuple[bool, str]:
    """Run the transfer script via ``abaqus cae noGUI=``.

    Returns ``(ok, combined_output)``; ok requires the TRANSFER_OK marker.
    The marker is read from ``transfer_to_cae_output.log`` written by the
    script itself -- on Windows the CAE kernel's prints never reach the
    captured stdout (only the license banner does), so stdout alone is
    NOT trusted. shell=True on Windows (abaqus is a .bat launcher).
    """
    workdir = os.path.dirname(script_path)
    log_path = os.path.join(workdir, OUTPUT_LOG)
    if os.path.exists(log_path):                 # stale marker from a prior run
        try:
            os.remove(log_path)
        except OSError:
            pass
    cmd = [abaqus_cmd, "cae", f"noGUI={os.path.basename(script_path)}"]
    try:
        proc = subprocess.run(cmd, cwd=workdir, timeout=timeout_s,
                              capture_output=True, text=True,
                              shell=(os.name == "nt"))
    except Exception as exc:
        return False, f"Could not launch abaqus cae: {exc!r}"
    out = (proc.stdout or "") + "\n" + (proc.stderr or "")
    if os.path.exists(log_path):
        try:
            with open(log_path, "r") as f:
                out += "\n--- kernel log (%s) ---\n%s" % (OUTPUT_LOG, f.read())
        except OSError:
            pass
    else:
        out += ("\n--- kernel log missing: the CAE kernel never ran the "
                "script (launcher/licence problem above) ---")
    ok = "TRANSFER_OK" in out
    return ok, out


def build_umat_library(tcfg: dict, abaqus_cmd: str = "abaqus",
                       timeout_s: int = 900) -> tuple[bool, str]:
    """Prebuild the user-subroutine DLL for a UMAT-model transfer.

    Runs ``abaqus make library=<multiaxial .cpp>`` (inside the MSVC vcvars
    environment) in the session's ``cae/`` folder and writes an
    ``abaqus_v6.env`` there containing the Intel-free ``link_sl`` override
    plus ``usub_lib_dir`` pointing at that folder. The user then submits the
    job from CAE normally: the solver loads the precompiled standardU.dll
    from ``usub_lib_dir`` — no compiler in CAE's environment and no
    user-subroutine job setting needed. No-op for native (Chaboche) mode.
    """
    if tcfg.get("material_mode") != "user":
        return True, "native material - no user-subroutine library needed"
    from .abaqus_runner import _UMAT_ENV_OVERRIDE, _find_vcvars
    cae_dir = os.path.dirname(os.path.abspath(tcfg["cae_path"]))
    umat = tcfg["umat_source"]
    vcvars = _find_vcvars()
    if vcvars is None:
        return False, ("vcvars64.bat not found (no MSVC C++ toolchain): "
                       "cannot compile the UMAT. Install VS Build Tools "
                       "with the C++ workload.")
    env_path = os.path.join(cae_dir, "abaqus_v6.env")
    with open(env_path, "w") as f:
        f.write(_UMAT_ENV_OVERRIDE)
        f.write("\n# Precompiled user-subroutine library (abaqus make): the\n"
                "# job loads standardU.dll from here - do NOT also set a\n"
                "# user subroutine on the job.\n")
        f.write("usub_lib_dir=r'%s'\n" % cae_dir)
    cmd = ('cd /d "%s" && call "%s" >nul 2>&1 && %s make library="%s" '
           % (cae_dir, vcvars, abaqus_cmd, umat))
    try:
        proc = subprocess.run('cmd /s /c "%s"' % cmd, timeout=timeout_s,
                              capture_output=True, text=True)
    except Exception as exc:
        return False, f"Could not run abaqus make: {exc!r}"
    out = (proc.stdout or "") + "\n" + (proc.stderr or "")
    # abaqus make drops the shared libraries in the working directory.
    dlls = [f for f in os.listdir(cae_dir) if f.lower().endswith(".dll")]
    if not dlls:
        return False, ("abaqus make produced no DLL in %s.\n%s"
                       % (cae_dir, out[-3000:]))
    return True, out + "\nUMAT library built: %s (usub_lib_dir=%s)" % (
        ", ".join(sorted(dlls)), cae_dir)


def write_transfer_log(session_dir: str, tcfg: dict, output: str = "") -> str:
    """Persist transfer_log.json after a successful transfer.

    ``cae_path`` is the session working COPY (the file the user submits the
    job from -- what post-run comparison should look next to); the untouched
    template is recorded as ``base_cae_path``.
    """
    stats = parse_mesh_stats(output)
    log = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "session_name": tcfg["session_name"],
        "model": tcfg["model"],
        "cae_path": tcfg["cae_path"],
        "base_cae_path": tcfg["base_cae_path"],
        "geometry": {
            "diameter_mm": tcfg["diameter_mm"],
            "gauge_length_mm": tcfg["gauge_length_mm"],
            "ld_ratio": tcfg["ld_ratio"],
            "rebuilt": bool(tcfg.get("rebuild_geometry")),
            "n_elements": stats["n_elements"],
            "n_failed_quality": stats["n_failed"],
        },
        "values_transferred": {
            "material_mode": tcfg.get("material_mode", "native"),
            "plastic_table": tcfg["plastic_table"],
            "cyclic_hardening_table": tcfg["cyclic_hardening_table"],
            "user_material_constants": tcfg.get("user_props"),
            "depvar": tcfg.get("nstatv"),
            "umat_source": tcfg.get("umat_source"),
            "usub_lib_dir": (os.path.dirname(os.path.abspath(
                tcfg["cae_path"]))
                if tcfg.get("material_mode") == "user" else None),
            "sigma_y0_cyclic": tcfg["sigma_y0_cyclic"],
            "amp_mm": tcfg["amp_mm"],
            "period_s": tcfg["period_s"],
            "n_cycles": tcfg["n_cycles"],
            "burnin_cycles": tcfg.get("burnin_cycles", 0),
            "step_name": tcfg["step_name"],
            "timePeriod": (tcfg["n_cycles"] + tcfg.get("burnin_cycles", 0))
                          * tcfg["period_s"],
            "initialInc": tcfg["period_s"] / 200.0,
            "minInc": tcfg["period_s"] / 1.0e6,
            "maxInc": tcfg["period_s"] / 50.0,
            "maxNumInc": 100000,
            "nlgeom": "ON",
        },
        # Extras plot_fe_comparison.py reads after the job completes.
        "strain_pct": tcfg["strain_pct"],
        "gauge_length_mm": tcfg["gauge_length_mm"],
        "diameter_mm": tcfg["diameter_mm"],
        "weight_reversals": tcfg["weight_reversals"],
        "surrogate_objective_MPa": tcfg["surrogate_objective_MPa"],
        "surrogate_peak_tension_MPa": tcfg["surrogate_peak_tension_MPa"],
        "surrogate_peak_compression_MPa": tcfg["surrogate_peak_compression_MPa"],
        "data_file": tcfg["data_file"],
    }
    path = os.path.join(session_dir, TRANSFER_LOG)
    with open(path, "w") as f:
        json.dump(log, f, indent=1)
    # Duplicate next to the working copy: the job's ODB lands there (work
    # dir = cae/), and plot_fe_comparison.py reads the log beside the ODB.
    cae_dir = os.path.dirname(tcfg["cae_path"])
    if os.path.abspath(cae_dir) != os.path.abspath(session_dir):
        try:
            with open(os.path.join(cae_dir, TRANSFER_LOG), "w") as f:
                json.dump(log, f, indent=1)
        except OSError:
            pass
    return path

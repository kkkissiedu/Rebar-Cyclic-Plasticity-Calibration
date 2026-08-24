"""
core/abaqus_runner.py
=====================

FE-in-the-loop evaluation backend.

For each candidate parameter set this writes a small single-element uniaxial
ABAQUS deck (axial along Z to match the real specimen, so it records **S33/E33**
- NOT the S22/E22 the legacy extractor wrongly used), submits it, extracts the
mid-point stress-strain history from the ODB, and scores it with the *same*
objective the surrogate backend uses (`objective.score_from_abaqus`).

The single-element mesh is used during the search for speed; the full cylinder
(`FATIGUE_rebuilt.cae` / rebuild script) is reserved for a final verification
run triggered separately from the GUI.

Nothing here is model-specific beyond asking the model for its
``material_block`` - so UVC / Ohno-Wang (UMAT) decks drop in unchanged once
those models exist, because ``material_block`` emits the right ``*USER MATERIAL``
block for them.

ABAQUS is invoked through a configurable command (default ``abaqus``; on this
machine that resolves via the SIMULIA Commands dir). Concurrency is configurable
(default 1) and applied to batch operations (Sobol seeding, top-M verification)
via a thread pool; the optimisers themselves call :meth:`evaluate` serially.
"""

from __future__ import annotations

import os
import subprocess
import threading
from concurrent.futures import ThreadPoolExecutor
from typing import Optional, Sequence

import numpy as np

from . import objective as objective_mod

# ---------------------------------------------------------------------------
# C++ UMAT support (2026-07-14). This machine has no Intel Fortran, so the
# UMAT models ship C++ ports (umats/*.cpp) compiled by MSVC. Two things are
# needed per job directory:
#   1. MSVC on PATH  -> the solve command is wrapped in `call vcvars64.bat`.
#   2. A local abaqus_v6.env overriding link_sl: the site env pulls Intel
#      runtime import libs via /DEFAULTLIB (LIBIFCOREMD etc.) that do not
#      exist here (LNK1104). The override /NODEFAULTLIBs the full Intel set;
#      the three libifcoremd symbols the ABAQUS stub objects still reference
#      (c_f_pointer_set_scalar / for_trim / for_concat) are shimmed inside
#      each UMAT .cpp. Toolchain proven end-to-end 2026-07-14 (elastic C++
#      UMAT on T3D2: S11 = E*NE11 exactly).
# ---------------------------------------------------------------------------
_UMAT_ENV_OVERRIDE = """\
# Auto-written by calibration_app (core/abaqus_runner.py) for C++ UMAT jobs.
# Local override: link user-subroutine DLL without Intel Fortran runtime
# import libs (none installed; user subroutines are MSVC-compiled C++).
link_sl=['LINK', '/nologo', '/NOENTRY', '/INCREMENTAL:NO',
         '/subsystem:console', '/machine:AMD64',
         '/NODEFAULTLIB:LIBC.LIB', '/NODEFAULTLIB:LIBCMT.LIB',
         '/DEFAULTLIB:OLDNAMES.LIB',
         '/NODEFAULTLIB:LIBIFCOREMD.LIB', '/NODEFAULTLIB:LIBIFPORTMD.LIB',
         '/NODEFAULTLIB:LIBMMD.LIB', '/NODEFAULTLIB:ifmodintr.lib',
         '/NODEFAULTLIB:ifconsol.lib', '/NODEFAULTLIB:libirc.lib',
         '/NODEFAULTLIB:svml_dispmd.lib', '/NODEFAULTLIB:IFWIN.LIB',
         '/DEFAULTLIB:kernel32.lib', '/DEFAULTLIB:user32.lib',
         '/DEFAULTLIB:advapi32.lib',
         '/FIXED:NO', '/dll', '/def:%E', '/out:%U', '%F', '%A', '%L', '%B',
         'oldnames.lib', 'user32.lib', 'ws2_32.lib', 'netapi32.lib',
         'advapi32.lib', 'msvcrt.lib', 'vcruntime.lib', 'ucrt.lib']
"""

_VCVARS_CANDIDATES = (
    r"C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat",
    r"C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvars64.bat",
    r"C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvars64.bat",
    r"C:\Program Files\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat",
    r"C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvars64.bat",
    r"C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Auxiliary\Build\vcvars64.bat",
)


def _find_vcvars() -> Optional[str]:
    """Locate vcvars64.bat (MSVC env) for C++ UMAT compilation."""
    for p in _VCVARS_CANDIDATES:
        if os.path.exists(p):
            return p
    return None

# ---------------------------------------------------------------------------
# Runtime-written ABAQUS-Python extractor. The app writes this into the work
# dir at run time (so it is NOT a tracked source file) and calls it with
# `abaqus python`. It pulls a configurable axial component (default S33/E33).
# Python 2.7 compatible (ABAQUS kernel).
# ---------------------------------------------------------------------------
_EXTRACTOR_SRC = r'''
from __future__ import print_function
import sys, os, math
from odbAccess import openOdb

def main():
    args = sys.argv[sys.argv.index("--")+1:] if "--" in sys.argv else sys.argv[1:]
    odb_path, csv_path = args[0], args[1]
    s_comp = args[2] if len(args) > 2 else "S33"
    e_comp = args[3] if len(args) > 3 else "E33"
    digs = "".join(c for c in e_comp if c.isdigit())
    odb = openOdb(path=odb_path, readOnly=True)
    try:
        step = odb.steps[list(odb.steps.keys())[-1]]
        s_data = {}
        strain = {"NE": {}, "E": {}, "LE": {}}   # engineering / small / logarithmic
        for rname, region in step.historyRegions.items():
            for key, ho in region.historyOutputs.items():
                if ho.data is None:
                    continue
                ku = key.upper()
                if ku == s_comp or ku.startswith(s_comp + " "):
                    for t, v in ho.data: s_data.setdefault(t, []).append(v)
                    continue
                for pre in ("NE", "LE", "E"):
                    tag = pre + digs
                    if ku == tag or ku.startswith(tag + " "):
                        for t, v in ho.data: strain[pre].setdefault(t, []).append(v)
                        break
        # Report ENGINEERING strain so the FE strain axis matches the
        # experimental engineering strain (displacement / gauge) used by the
        # surrogate objective -- otherwise the comparison drifts at 5-6% amp.
        # NE (nominal) and E (small) are already engineering; LE (logarithmic /
        # true) is converted: eng = exp(LE) - 1.
        if strain["NE"]:
            e_src, convert, kind = strain["NE"], False, "NE(engineering)"
        elif strain["E"]:
            e_src, convert, kind = strain["E"], False, "E(small)"
        else:
            e_src, convert, kind = strain["LE"], True, "LE->engineering"
        times = sorted(set(s_data.keys()) & set(e_src.keys()))
        f = open(csv_path, "w")
        try:
            f.write("time,strain,stress_MPa\n")
            for t in times:
                sv = sum(s_data[t])/float(len(s_data[t]))
                ev = sum(e_src[t])/float(len(e_src[t]))
                if convert:
                    ev = math.exp(ev) - 1.0
                f.write("%.10g,%.10g,%.10g\n" % (t, ev, sv))
        finally:
            f.close()
        sys.stdout.write("wrote %d rows (strain=%s)\n" % (len(times), kind))
    finally:
        odb.close()

if __name__ == "__main__":
    main()
'''


class FEBackend:
    """Real-ABAQUS evaluation backend, API-compatible with the surrogate one.

    Parameters
    ----------
    model : CyclicPlasticityModel
        Supplies the ABAQUS material block for each candidate.
    window : dict
        Calibration window (from ``data_loader.slice_window``) used to score and
        to size the amplitude / number of cycles.
    work_dir : str
        Directory for per-candidate job files (a session's ``abaqus_runs/``).
    abaqus_cmd : str
        Command to invoke ABAQUS (default ``abaqus``).
    E, nu, sy0 : float
        Elastic constants and per-run yield stress.
    gauge_length_mm, amp_mm, period_s, n_cycles : float / int
        Loading definition; ``amp_mm`` is the peak displacement.
    axial_stress, axial_strain : str
        History components to extract (default ``S33``/``E33`` for the Z-axis
        specimen).
    concurrent_jobs : int
        Max simultaneous ABAQUS jobs for batch operations (default 1).
    weight_reversals : bool
        Passed through to the shared objective.
    """

    def __init__(
        self,
        model,
        window: dict,
        *,
        work_dir: str,
        abaqus_cmd: str = "abaqus",
        E: float,
        nu: float,
        sy0: float,
        gauge_length_mm: float,
        amp_mm: float,
        period_s: float = 20.0,
        n_cycles: int = 5,
        axial_stress: str = "S33",
        axial_strain: str = "E33",
        concurrent_jobs: int = 1,
        weight_reversals: bool = True,
        timeout_s: int = 1800,
    ) -> None:
        self.model = model
        self.window = window
        self.work_dir = work_dir
        self.abaqus_cmd = abaqus_cmd
        self.E, self.nu, self.sy0 = E, nu, sy0
        self.gauge_length_mm = gauge_length_mm
        self.amp_mm = amp_mm
        self.period_s = period_s
        self.n_cycles = int(n_cycles)
        self.axial_stress = axial_stress
        self.axial_strain = axial_strain
        self.concurrent_jobs = max(1, int(concurrent_jobs))
        self.weight_reversals = weight_reversals
        self.timeout_s = timeout_s
        self._job_counter = 0
        self._lock = threading.Lock()
        # UMAT models (UVC, Ohno-Wang) expose a `umat_file` (relative to the
        # project root). For those the FE deck is a uniaxial truss (T3D2, one
        # stress component -> matches the uniaxial UMAT) submitted with `user=`.
        umat_rel = getattr(model, "umat_file", None)
        if umat_rel:
            proj_root = os.path.dirname(os.path.dirname(os.path.dirname(
                os.path.abspath(__file__))))
            self.umat_abs = os.path.join(proj_root, umat_rel)
            self.axial_stress, self.axial_strain = "S11", "E11"   # truss axial
        else:
            self.umat_abs = None

        os.makedirs(self.work_dir, exist_ok=True)
        # C++ UMAT jobs need MSVC (vcvars) + the local link_sl env override;
        # both are no-ops for the native-Chaboche deck.
        self.vcvars = _find_vcvars() if self.umat_abs else None
        if self.umat_abs and str(self.umat_abs).lower().endswith(
                (".cpp", ".c", ".cxx", ".cc")):
            env_path = os.path.join(self.work_dir, "abaqus_v6.env")
            if not os.path.exists(env_path):
                with open(env_path, "w") as f:
                    f.write(_UMAT_ENV_OVERRIDE)
        # Write the extractor once per backend instance.
        self._extractor = os.path.join(self.work_dir, "_extract_component.py")
        with open(self._extractor, "w") as f:
            f.write(_EXTRACTOR_SRC)

    # -- deck generation -----------------------------------------------------
    #: One idealised triangular cycle run and discarded before the scored
    #: cycles, mirroring data_loader's cycle-1 burn-in convention: the
    #: experimental scoring window starts mid-hysteresis (already cyclically
    #: loaded), not from a virgin, zero-stress state. Without this the FE
    #: deck's first scored cycle starts at 0 MPa while the matching
    #: experimental cycle starts at its already-built-up stress, which the
    #: arc-length resampler cannot correct -- inflating the FE objective by
    #: ~4x (206 MPa vs ~46 MPa surrogate on the same params).
    BURNIN_CYCLES: int = 1

    def _amplitude_lines(self) -> list[str]:
        """Tabular triangle amplitude, peak = amp_mm; BURNIN_CYCLES + n_cycles
        total cycles are written, `run_candidate` strips the burn-in portion
        from the extracted trace before it reaches the objective."""
        strain_amp = self.amp_mm / self.gauge_length_mm      # unit-cube disp
        total_cycles = self.BURNIN_CYCLES + self.n_cycles
        pts = []
        for c in range(total_cycles):
            t0 = c * self.period_s
            pts += [(t0, 0.0),
                    (t0 + self.period_s / 4, strain_amp),
                    (t0 + self.period_s / 2, 0.0),
                    (t0 + 3 * self.period_s / 4, -strain_amp)]
        pts.append((total_cycles * self.period_s, 0.0))
        return [f"{t:.6f}, {a:.8f}" for t, a in pts]

    def _write_single_element_inp(self, path: str, params: Sequence[float],
                                  sy0: float = None) -> None:
        """Single C3D8R unit cube, axial along Z (dof 3) to match the specimen.

        Edge = 1 mm so prescribed top-face Z displacement equals the strain
        amplitude; history records S33/E33 on the single element. ``sy0`` may be
        a per-candidate value when the yield stress is being calibrated.
        """
        mat = self.model.material_block(params, E=self.E, nu=self.nu,
                                        sy0=self.sy0 if sy0 is None else sy0,
                                        name="STEEL")
        L = [
            "*HEADING",
            "single-element uniaxial (axial = Z)",
            "*PREPRINT, ECHO=NO, MODEL=NO, HISTORY=NO",
            "*NODE",
            "1, 0.0, 0.0, 0.0", "2, 1.0, 0.0, 0.0", "3, 1.0, 1.0, 0.0",
            "4, 0.0, 1.0, 0.0", "5, 0.0, 0.0, 1.0", "6, 1.0, 0.0, 1.0",
            "7, 1.0, 1.0, 1.0", "8, 0.0, 1.0, 1.0",
            "*ELEMENT, TYPE=C3D8R, ELSET=E_ALL",
            "1, 1, 2, 3, 4, 5, 6, 7, 8",
            "*NSET, NSET=N_BOT",       # z = 0 face
            "1, 2, 3, 4",
            "*NSET, NSET=N_TOP",       # z = 1 face
            "5, 6, 7, 8",
            "*SOLID SECTION, ELSET=E_ALL, MATERIAL=STEEL",
            "1.0",
            mat.rstrip("\n"),
            "*AMPLITUDE, NAME=AMP_TRI, DEFINITION=TABULAR",
        ]
        L += self._amplitude_lines()
        L += [
            "*BOUNDARY",
            "N_BOT, 3, 3, 0.0",         # fix axial (z) at bottom
            "1, 1, 2",                  # pin one node in-plane (no rigid body)
            "2, 2, 2",
            "*STEP, NLGEOM=YES, INC=100000",
            "*STATIC",
            f"{self.period_s/200.0}, {(self.BURNIN_CYCLES + self.n_cycles)*self.period_s}, "
            f"{self.period_s/1.0e6}, {self.period_s/50.0}",
            "*BOUNDARY, AMPLITUDE=AMP_TRI",
            "N_TOP, 3, 3, 1.0",         # displacement-controlled top face (z)
            "*OUTPUT, HISTORY, FREQ=1",
            "*ELEMENT OUTPUT, ELSET=E_ALL",
            # NLGEOM=ON => request NE (nominal = engineering strain) as the
            # primary axial strain so it matches the experimental engineering
            # strain used by the surrogate objective; LE (logarithmic) is also
            # requested as a fallback and the extractor converts it if needed.
            "%s, NE%s, LE%s" % (self.axial_stress,
                                "".join(c for c in self.axial_strain if c.isdigit()),
                                "".join(c for c in self.axial_strain if c.isdigit())),
            "*END STEP",
        ]
        with open(path, "w") as f:
            f.write("\n".join(L) + "\n")

    def _write_umat_truss_inp(self, path: str, params: Sequence[float],
                              sy0: float = None) -> None:
        """Uniaxial truss (T3D2) deck driven by the model's UMAT (submitted with
        `user=`). A truss has a single axial stress component (NTENS=1), matching
        the uniaxial UMAT; records S11 + NE11 (engineering)."""
        mat = self.model.material_block(params, E=self.E, nu=self.nu,
                                        sy0=self.sy0 if sy0 is None else sy0,
                                        name="STEEL")
        L = [
            "*HEADING",
            "single truss uniaxial (UMAT)",
            "*PREPRINT, ECHO=NO, MODEL=NO, HISTORY=NO",
            "*NODE",
            "1, 0.0, 0.0, 0.0",
            "2, 0.0, 0.0, 1.0",
            "*ELEMENT, TYPE=T3D2, ELSET=E_ALL",
            "1, 1, 2",
            "*SOLID SECTION, ELSET=E_ALL, MATERIAL=STEEL",
            "1.0",
            mat.rstrip("\n"),
            "*AMPLITUDE, NAME=AMP_TRI, DEFINITION=TABULAR",
        ]
        L += self._amplitude_lines()
        L += [
            "*BOUNDARY",
            "1, 1, 3, 0.0",             # fix node 1 fully
            "2, 1, 2, 0.0",             # node 2: fix x,y; axial (z) is driven
            "*STEP, NLGEOM=YES, INC=100000",
            "*STATIC",
            f"{self.period_s/200.0}, {(self.BURNIN_CYCLES + self.n_cycles)*self.period_s}, "
            f"{self.period_s/1.0e6}, {self.period_s/50.0}",
            "*BOUNDARY, AMPLITUDE=AMP_TRI",
            "2, 3, 3, 1.0",
            "*OUTPUT, HISTORY, FREQ=1",
            "*ELEMENT OUTPUT, ELSET=E_ALL",
            "S11, NE11, LE11",
            "*END STEP",
        ]
        with open(path, "w") as f:
            f.write("\n".join(L) + "\n")

    # -- job submission ------------------------------------------------------
    def _next_jobname(self) -> str:
        with self._lock:
            self._job_counter += 1
            return f"fe_cand_{self._job_counter:05d}"

    def _write_stderr_log(self, job: str, stage: str, stdout: str,
                          stderr: str) -> None:
        """Persist captured stdout+stderr to ``{job}_abaqus_stderr.log`` on the
        failure path so a future UMAT compile failure (or any FE crash) leaves a
        readable trace on disk instead of vanishing. Best-effort: never raises.
        """
        log_path = os.path.join(self.work_dir, job + "_abaqus_stderr.log")
        try:
            with open(log_path, "w") as f:
                f.write("=== FE failure trace for job %s (stage: %s) ===\n" % (job, stage))
                f.write("\n--- STDOUT ---\n")
                f.write(stdout or "")
                f.write("\n--- STDERR ---\n")
                f.write(stderr or "")
                f.write("\n")
        except Exception:
            pass

    def run_candidate(self, params: Sequence[float], sy0: float = None) -> Optional[tuple[np.ndarray, np.ndarray]]:
        """Run one single-element job; return (strain, stress) or None on failure.

        On any failure path the captured ABAQUS stdout+stderr is written to
        ``{jobname}_abaqus_stderr.log`` in ``work_dir`` before returning None.
        """
        job = self._next_jobname()
        inp = os.path.join(self.work_dir, job + ".inp")
        if self.umat_abs:
            if not os.path.exists(self.umat_abs):
                self._write_stderr_log(job, "umat_missing", "",
                                       "UMAT file not found: %s" % self.umat_abs)
                return None                        # UMAT missing -> FE unavailable
            if self.vcvars is None:
                self._write_stderr_log(
                    job, "msvc_missing", "",
                    "vcvars64.bat not found (no MSVC): cannot compile the "
                    "C++ UMAT %s. Install VS Build Tools with C++ workload."
                    % self.umat_abs)
                return None
            self._write_umat_truss_inp(inp, params, sy0)
            # single cmd string: load the MSVC env, then run ABAQUS in it
            cmd = ('call "%s" >nul 2>&1 && %s job=%s input=%s.inp '
                   'user="%s" interactive ask_delete=OFF'
                   % (self.vcvars, self.abaqus_cmd, job, job, self.umat_abs))
        else:
            self._write_single_element_inp(inp, params, sy0)
            cmd = [self.abaqus_cmd, f"job={job}", f"input={job}.inp",
                   "interactive", "ask_delete=OFF"]
        try:
            proc = subprocess.run(cmd, cwd=self.work_dir, timeout=self.timeout_s,
                                  capture_output=True, text=True,
                                  shell=(os.name == "nt"))
        except Exception as exc:
            self._write_stderr_log(job, "solve_exception", "", repr(exc))
            return None
        odb = os.path.join(self.work_dir, job + ".odb")
        if not os.path.exists(odb):
            self._write_stderr_log(job, "solve_no_odb", proc.stdout, proc.stderr)
            return None
        csv_out = os.path.join(self.work_dir, job + ".csv")
        ex = [self.abaqus_cmd, "python", os.path.basename(self._extractor), "--",
              job + ".odb", job + ".csv", self.axial_stress, self.axial_strain]
        try:
            eproc = subprocess.run(ex, cwd=self.work_dir, timeout=300,
                                   capture_output=True, text=True,
                                   shell=(os.name == "nt"))
        except Exception as exc:
            self._write_stderr_log(job, "extract_exception", "", repr(exc))
            return None
        if not os.path.exists(csv_out):
            self._write_stderr_log(job, "extract_no_csv", eproc.stdout, eproc.stderr)
            return None
        arr = np.loadtxt(csv_out, delimiter=",", skiprows=1)
        if arr.ndim != 2 or arr.shape[0] < 4:
            self._write_stderr_log(job, "bad_csv",
                                   "csv rows=%s" % (getattr(arr, "shape", None),), "")
            return None
        # Drop the burn-in cycle(s) run for pre-conditioning (see
        # BURNIN_CYCLES) so the caller only ever sees the n_cycles that match
        # the calibration window -- same contract as before this was added.
        cutoff = self.BURNIN_CYCLES * self.period_s
        scored = arr[arr[:, 0] >= cutoff - 1e-9]
        if scored.shape[0] < 4:
            self._write_stderr_log(job, "bad_csv_after_burnin_strip",
                                   "rows=%d cutoff=%s" % (scored.shape[0], cutoff), "")
            return None
        return scored[:, 1], scored[:, 2]    # strain, stress (scored cycles only)

    def evaluate(self, params: Sequence[float], sy0: float = None) -> float:
        """Objective for one candidate - the callable handed to the optimiser.

        Uses the identical scorer as the surrogate backend so the four study
        conditions are directly comparable. ``sy0`` overrides the yield stress
        per candidate when the yield is being calibrated.
        """
        res = self.run_candidate(params, sy0)
        if res is None:
            return objective_mod.PENALTY
        strain, stress = res
        return objective_mod.score_from_abaqus(strain, stress, self.window,
                                               self.weight_reversals)

    def evaluate_batch(self, param_sets: Sequence[Sequence[float]]) -> list[float]:
        """Evaluate several candidates concurrently (up to ``concurrent_jobs``).

        Used for Sobol seeding / top-M verification where independent jobs can
        run in parallel; the optimiser loops still call :meth:`evaluate` serially.
        """
        if self.concurrent_jobs == 1:
            return [self.evaluate(p) for p in param_sets]
        with ThreadPoolExecutor(max_workers=self.concurrent_jobs) as ex:
            return list(ex.map(self.evaluate, param_sets))

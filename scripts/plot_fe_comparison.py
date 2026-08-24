# -*- coding: utf-8 -*-
"""
plot_fe_comparison.py
=====================

Standalone post-run FE-vs-experiment comparison for FATIGUE_rebuilt.cae jobs.

Usage (after the ABAQUS job completes):

    abaqus python plot_fe_comparison.py -- <odb_path> <exp_csv>

Runs under ABAQUS Python (odbAccess), NOT the cuda_pt env: pure stdlib only
(no numpy/scipy). matplotlib is attempted for the PNG and skipped gracefully
if the ABAQUS Python lacks it -- the calibration app can render the PNG from
fe_comparison.csv instead.

Outputs (saved alongside the ODB):
    fe_comparison.csv   time, strain_FE, stress_FE_MPa,
                        strain_exp_interp, stress_exp_MPa
    fe_comparison.json  machine-readable summary (objective, peaks) for the
                        app's "Run post-run comparison" button
    fe_comparison.png   publication plot (DPI=150), if matplotlib available

The objective is the SAME per-cycle weighted RMSE as the calibration app
(core/objective.py): FE trace split into N_fit equal chunks, resampled onto
the experimental strain points per cycle by normalised arc-length, then
mean over cycles of sqrt(weighted mean squared error); samples in the
first/last 15% of each cycle's strain arc-length weighted 1.5x.

Metadata (model, N_fit, geometry, surrogate objective) is read from
transfer_log.json in the ODB's folder or, missing that, from the
experimental file path conventions ({N}data.csv, "{D}mm dia", "LD{ratio}").
"""

from __future__ import print_function

import csv
import json
import math
import os
import re
import sys

# -- objective constants (MUST match calibration_app/core/objective.py) ------
REVERSAL_FRACTION = 0.15
REVERSAL_WEIGHT = 1.5
PENALTY = 1.0e9

_STRAIN_RE = re.compile(r"(\d+)\s*data\.csv$", re.IGNORECASE)
_DIA_RE = re.compile(r"(\d+)\s*mm\s*dia", re.IGNORECASE)
_LD_RE = re.compile(r"LD\s*(\d+)", re.IGNORECASE)


# ---------------------------------------------------------------------------
# small pure-python numerics (no numpy under ABAQUS python by requirement)
# ---------------------------------------------------------------------------
def interp(x, xp, fp):
    """Linear interpolation of fp(xp) at points x. xp must be ascending."""
    out = []
    n = len(xp)
    j = 0
    for xi in x:
        if xi <= xp[0]:
            out.append(fp[0]); continue
        if xi >= xp[-1]:
            out.append(fp[-1]); continue
        while j < n - 2 and xp[j + 1] < xi:
            j += 1
        while j > 0 and xp[j] > xi:
            j -= 1
        x0, x1 = xp[j], xp[j + 1]
        f0, f1 = fp[j], fp[j + 1]
        out.append(f0 if x1 == x0 else f0 + (f1 - f0) * (xi - x0) / (x1 - x0))
    return out


def arc_norm(strain):
    """Normalised cumulative |d(strain)| arc-length, or None if degenerate."""
    arc = [0.0]
    for i in range(1, len(strain)):
        arc.append(arc[-1] + abs(strain[i] - strain[i - 1]))
    total = arc[-1]
    if total <= 0:
        return None
    return [a / total for a in arc]


def cycle_rmse(strain_c, sim_c, exp_c, weight_reversals):
    """Weighted RMSE (MPa) over one cycle -- mirrors objective._cycle_rmse."""
    if not weight_reversals:
        s = sum((a - b) ** 2 for a, b in zip(sim_c, exp_c))
        return math.sqrt(s / float(len(sim_c)))
    frac = arc_norm(strain_c)
    if frac is None:
        frac = [0.0] * len(strain_c)
    num = den = 0.0
    for f, a, b in zip(frac, sim_c, exp_c):
        w = (REVERSAL_WEIGHT
             if (f <= REVERSAL_FRACTION or f >= 1.0 - REVERSAL_FRACTION)
             else 1.0)
        num += w * (a - b) ** 2
        den += w
    return math.sqrt(num / den) if den > 0 else float("nan")


# ---------------------------------------------------------------------------
# experimental CSV (Kashani layout -- mirrors core/data_loader.py, stdlib only)
# ---------------------------------------------------------------------------
def load_experimental(path, gauge_mm, diameter_mm, ref_cycle=2):
    """Return dict of lists: time, cycle, strain, stress for cycles >= 2."""
    area = math.pi * (diameter_mm / 2.0) ** 2
    times, cycles, positions, loads = [], [], [], []
    f = open(path, "r")
    try:
        reader = csv.reader(f)
        header = None
        for raw in reader:
            if not raw:
                continue
            joined = ",".join(s.strip() for s in raw).lower()
            if joined.startswith("notes") or joined.startswith('"notes"'):
                header = raw
                break
        if header is None:
            raise RuntimeError("Header row (containing 'Notes') not found "
                               "in %s" % path)
        clean = [c.strip().strip('"').lower() for c in header]
        i_t = clean.index("time")
        i_c = clean.index("cycle")
        i_p = clean.index("position mm")
        i_l = clean.index("load kn")
        for raw in reader:
            if not raw or len(raw) <= max(i_t, i_c, i_p, i_l):
                continue
            try:
                t = float(raw[i_t]); c = int(float(raw[i_c]))
                pos = float(raw[i_p]); ld = float(raw[i_l])
            except ValueError:
                continue
            times.append(t); cycles.append(c)
            positions.append(pos); loads.append(ld)
    finally:
        f.close()
    if not cycles:
        raise RuntimeError("No numeric rows parsed from %s" % path)
    ref_position = None
    for c, p in zip(cycles, positions):
        if c == ref_cycle:
            ref_position = p
            break
    if ref_position is None:
        raise RuntimeError("Reference cycle %d not present in %s"
                           % (ref_cycle, path))
    out = {"time": [], "cycle": [], "strain": [], "stress": []}
    for t, c, p, ld in zip(times, cycles, positions, loads):
        if c < ref_cycle:
            continue
        out["time"].append(t)
        out["cycle"].append(c)
        out["strain"].append((p - ref_position) / gauge_mm)
        out["stress"].append(ld * 1000.0 / area)
    return out


# ---------------------------------------------------------------------------
# ODB history extraction (S33 + NE33/E33/LE33, SS_AUTO single element)
# ---------------------------------------------------------------------------
def extract_odb(odb_path, gauge_mm, area_mm2):
    """Return (nominal, local) traces from the last step's history.

    Each is a (time, strain_eng, stress_MPa) tuple or None.

    nominal -- stress = sum(RF3 over TOP_AUTO)/A0, strain = mean(U3)/gauge.
      This is the direct analogue of the experiment (load cell / crosshead
      over gauge) and is IMMUNE to strain localisation: with nlgeom + a
      slender specimen the mid-gauge element's local strain drifts away
      from the nominal cycle (barrelling/incipient buckling), which is a
      structural effect, not a calibration error.
    local -- the SS_AUTO element's S33 vs NE33 (LE converted to
      engineering), kept as a localisation diagnostic.
    """
    from odbAccess import openOdb
    odb = openOdb(path=odb_path, readOnly=True)
    try:
        step = odb.steps[list(odb.steps.keys())[-1]]
        s_data = {}
        strain = {"NE": {}, "E": {}, "LE": {}}
        rf3 = {}
        u3 = {}
        for rname in step.historyRegions.keys():
            region = step.historyRegions[rname]
            for key in region.historyOutputs.keys():
                ho = region.historyOutputs[key]
                if ho.data is None:
                    continue
                ku = key.upper()
                if ku == "S33" or ku.startswith("S33 "):
                    for t, v in ho.data:
                        s_data.setdefault(t, []).append(v)
                    continue
                if ku == "RF3" or ku.startswith("RF3 "):
                    for t, v in ho.data:
                        rf3.setdefault(t, []).append(v)
                    continue
                if ku == "U3" or ku.startswith("U3 "):
                    for t, v in ho.data:
                        u3.setdefault(t, []).append(v)
                    continue
                for pre in ("NE", "LE", "E"):
                    tag = pre + "33"
                    if ku == tag or ku.startswith(tag + " "):
                        for t, v in ho.data:
                            strain[pre].setdefault(t, []).append(v)
                        break

        nominal = None
        times_n = sorted(set(rf3.keys()) & set(u3.keys()))
        if times_n:
            t_out, e_out, s_out = [], [], []
            for t in times_n:
                force_n = sum(rf3[t])                       # N (sum over face)
                disp = sum(u3[t]) / float(len(u3[t]))       # mm (mean)
                t_out.append(t)
                e_out.append(disp / gauge_mm)
                s_out.append(force_n / area_mm2)            # N/mm^2 = MPa
            nominal = (t_out, e_out, s_out)

        local = None
        if strain["NE"]:
            e_src, convert = strain["NE"], False
        elif strain["E"]:
            e_src, convert = strain["E"], False
        else:
            e_src, convert = strain["LE"], True
        times_l = sorted(set(s_data.keys()) & set(e_src.keys()))
        if times_l:
            t_out, e_out, s_out = [], [], []
            for t in times_l:
                sv = sum(s_data[t]) / float(len(s_data[t]))
                ev = sum(e_src[t]) / float(len(e_src[t]))
                if convert:
                    ev = math.exp(ev) - 1.0
                t_out.append(t); e_out.append(ev); s_out.append(sv)
            local = (t_out, e_out, s_out)

        if nominal is None and local is None:
            raise RuntimeError(
                "No usable history in the ODB (need RF3/U3 on TOP_AUTO or "
                "S33/NE33 on SS_AUTO) - was the job run after 'Transfer to "
                "CAE' set up the output requests?")
        return nominal, local
    finally:
        odb.close()


# ---------------------------------------------------------------------------
# metadata
# ---------------------------------------------------------------------------
def load_metadata(odb_path, exp_csv):
    """Merge transfer_log.json (if present next to the ODB) with path-derived
    geometry. Returns a dict with model, n_fit, gauge, dia, strain_pct, etc."""
    meta = {"model": "chaboche", "n_fit": None, "gauge_length_mm": None,
            "diameter_mm": None, "strain_pct": None,
            "weight_reversals": True, "surrogate_objective_MPa": None,
            "sigma_y0": None, "period_s": 20.0, "burnin_cycles": 0}
    log_path = os.path.join(os.path.dirname(os.path.abspath(odb_path)),
                            "transfer_log.json")
    if os.path.exists(log_path):
        f = open(log_path, "r")
        try:
            log = json.load(f)
        finally:
            f.close()
        meta["model"] = log.get("model", meta["model"])
        vt = log.get("values_transferred", {})
        meta["n_fit"] = vt.get("n_cycles")
        meta["period_s"] = vt.get("period_s", 20.0)
        meta["burnin_cycles"] = int(vt.get("burnin_cycles", 0))
        meta["sigma_y0"] = vt.get("sigma_y0_cyclic")
        meta["gauge_length_mm"] = log.get("gauge_length_mm")
        meta["diameter_mm"] = log.get("diameter_mm")
        meta["strain_pct"] = log.get("strain_pct")
        meta["weight_reversals"] = log.get("weight_reversals", True)
        meta["surrogate_objective_MPa"] = log.get("surrogate_objective_MPa")
        print("Read transfer_log.json: model=%s, N_fit=%s"
              % (meta["model"], meta["n_fit"]))
    # Path-derived fallbacks (Kashani conventions).
    norm = exp_csv.replace("\\", "/")
    if meta["strain_pct"] is None:
        m = _STRAIN_RE.search(os.path.basename(exp_csv))
        meta["strain_pct"] = float(m.group(1)) if m else None
    if meta["diameter_mm"] is None:
        m = _DIA_RE.search(norm)
        meta["diameter_mm"] = float(m.group(1)) if m else None
    if meta["gauge_length_mm"] is None:
        m = _LD_RE.search(norm)
        if m and meta["diameter_mm"]:
            meta["gauge_length_mm"] = float(m.group(1)) * meta["diameter_mm"]
    if not meta["gauge_length_mm"] or not meta["diameter_mm"]:
        raise RuntimeError(
            "Could not determine geometry: no transfer_log.json next to the "
            "ODB and the experimental path does not match "
            "'{D}mm dia/LD{ratio}/{N}data.csv'.")
    return meta


# ---------------------------------------------------------------------------
# objective + resampling
# ---------------------------------------------------------------------------
def group_cycles(exp):
    """Ordered list of (cycle_no, [indices]) with >= 4 samples."""
    order, groups = [], {}
    for i, c in enumerate(exp["cycle"]):
        if c not in groups:
            groups[c] = []
            order.append(c)
        groups[c].append(i)
    return [(c, groups[c]) for c in order if len(groups[c]) >= 4]


def time_chunks(fe_time, n_cyc, period_s):
    """Index ranges [(lo, hi), ...] splitting the FE trace into cycles BY
    TIME (cycle k = t in [k*T, (k+1)*T]). Equal-count chunking is WRONG
    here: automatic incrementation packs points unevenly (dense at
    reversals), so count-based boundaries drift across cycles and the
    per-cycle pairing degrades into garbage RMSE despite a perfect fit."""
    out = []
    lo = 0
    n = len(fe_time)
    for k in range(n_cyc):
        t_end = (k + 1) * period_s
        hi = lo
        while hi < n and (fe_time[hi] <= t_end + 1e-9 or k == n_cyc - 1):
            hi += 1
        out.append((lo, hi))
        lo = hi
    return out


def fe_objective(fe_time, fe_strain, fe_stress, exp, target,
                 weight_reversals, period_s):
    """Per-cycle weighted RMSE, FE resampled onto exp points by arc-length --
    the same function of the response as core/objective.score_from_abaqus."""
    n_cyc = len(target)
    if n_cyc == 0 or len(fe_strain) < 4:
        return PENALTY
    rmses = []
    for k, ((c, idx), (lo, hi)) in enumerate(
            zip(target, time_chunks(fe_time, n_cyc, period_s))):
        if hi - lo < 4:
            continue
        e_exp = [exp["strain"][i] for i in idx]
        s_exp = [exp["stress"][i] for i in idx]
        e_sim = fe_strain[lo:hi]
        s_sim = fe_stress[lo:hi]
        a_exp = arc_norm(e_exp)
        a_sim = arc_norm(e_sim)
        if a_exp is None or a_sim is None:
            continue
        resampled = interp(a_exp, a_sim, s_sim)
        r = cycle_rmse(e_exp, resampled, s_exp, weight_reversals)
        if not math.isnan(r):
            rmses.append(r)
    if not rmses:
        return PENALTY
    return sum(rmses) / float(len(rmses))


def exp_onto_fe(fe_time, fe_strain, exp, target, period_s):
    """Experimental (strain, stress) resampled onto the FE points per cycle
    by arc-length, for the side-by-side fe_comparison.csv columns."""
    n_cyc = len(target)
    n = len(fe_strain)
    e_out = [""] * n
    s_out = [""] * n
    if not n_cyc:
        return e_out, s_out
    for (c, idx), (lo, hi) in zip(target,
                                  time_chunks(fe_time, n_cyc, period_s)):
        if hi - lo < 2:
            continue
        e_sim = fe_strain[lo:hi]
        a_sim = arc_norm(e_sim)
        e_exp = [exp["strain"][i] for i in idx]
        s_exp = [exp["stress"][i] for i in idx]
        a_exp = arc_norm(e_exp)
        if a_sim is None or a_exp is None:
            continue
        ei = interp(a_sim, a_exp, e_exp)
        si = interp(a_sim, a_exp, s_exp)
        for j in range(lo, hi):
            e_out[j] = ei[j - lo]
            s_out[j] = si[j - lo]
    return e_out, s_out


# ---------------------------------------------------------------------------
# outputs
# ---------------------------------------------------------------------------
def write_csv(path, fe_time, fe_strain, fe_stress, e_exp_i, s_exp_i):
    f = open(path, "w")
    try:
        f.write("time,strain_FE,stress_FE_MPa,strain_exp_interp,"
                "stress_exp_MPa\n")
        for row in zip(fe_time, fe_strain, fe_stress, e_exp_i, s_exp_i):
            f.write(",".join("" if v == "" else "%.10g" % v for v in row)
                    + "\n")
    finally:
        f.close()


def try_png(path, fe_strain, fe_stress, exp, target, meta, annotation):
    """Publication plot; returns True if written (matplotlib available)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        print("PNG skipped (matplotlib unavailable in this Python): %s" % exc)
        return False
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    first = True
    for c, idx in target:
        ax.plot([exp["strain"][i] * 100.0 for i in idx],
                [exp["stress"][i] for i in idx],
                color="0.6", lw=0.9,
                label=("Experiment (cycles %s-%s)"
                       % (target[0][0], target[-1][0])) if first else None)
        first = False
    ax.plot([e * 100.0 for e in fe_strain], fe_stress,
            color="tab:red", lw=1.2,
            label="FE simulation (%s)" % meta.get("channel", "nominal"))
    ld = (meta["gauge_length_mm"] / meta["diameter_mm"]
          if meta["diameter_mm"] else 0)
    ax.set_title("FE verification - %s - %g%% strain - %gmm LD%g"
                 % (meta["model"], meta["strain_pct"] or 0,
                    meta["diameter_mm"] or 0, ld))
    ax.text(0.02, 0.02, annotation, transform=ax.transAxes, fontsize=8,
            va="bottom",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
    ax.set_xlabel("Strain (%)")
    ax.set_ylabel("Stress (MPa)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return True


def main():
    args = (sys.argv[sys.argv.index("--") + 1:]
            if "--" in sys.argv else sys.argv[1:])
    if len(args) < 2:
        print("Usage: abaqus python plot_fe_comparison.py -- "
              "<odb_path> <exp_csv>")
        sys.exit(2)
    odb_path, exp_csv = args[0], args[1]
    out_dir = os.path.dirname(os.path.abspath(odb_path))

    meta = load_metadata(odb_path, exp_csv)
    exp = load_experimental(exp_csv, meta["gauge_length_mm"],
                            meta["diameter_mm"])
    area = math.pi * (meta["diameter_mm"] / 2.0) ** 2
    nominal, local = extract_odb(odb_path, meta["gauge_length_mm"], area)

    # Drop the burn-in cycle (state conditioning only, never scored) and
    # rebase time so the scored cycles start at t=0 for the cycle windows.
    off = meta["burnin_cycles"] * meta["period_s"]
    if off > 0:
        def _trim(trace):
            if trace is None:
                return None
            tt, ee, ss = trace
            keep = [i for i, tv in enumerate(tt) if tv >= off]
            return ([tt[i] - off for i in keep], [ee[i] for i in keep],
                    [ss[i] for i in keep])
        nominal = _trim(nominal)
        local = _trim(local)
        print("Burn-in trimmed: first %.0f s (%d cycle(s)) excluded from "
              "scoring" % (off, meta["burnin_cycles"]))

    # Primary channel: nominal (load-cell analogue). The local mid-gauge
    # element is a localisation diagnostic only -- under nlgeom a slender
    # specimen's local strain drifts from the nominal cycle.
    if nominal is not None:
        channel = "nominal (RF3/U3)"
        fe_time, fe_strain, fe_stress = nominal
    else:
        channel = "local element (S33/NE33) - RF3/U3 missing in ODB"
        fe_time, fe_strain, fe_stress = local
    print("Channel: %s" % channel)
    print("FE trace: %d points; experiment: %d points"
          % (len(fe_time), len(exp["time"])))

    cycles = group_cycles(exp)
    n_fit = meta["n_fit"] or len(cycles)
    target = cycles[:int(n_fit)]

    obj = fe_objective(fe_time, fe_strain, fe_stress, exp, target,
                       meta["weight_reversals"], meta["period_s"])
    obj_local = None
    if nominal is not None and local is not None:
        obj_local = fe_objective(local[0], local[1], local[2], exp, target,
                                 meta["weight_reversals"], meta["period_s"])
        drift = max(abs(v) for v in local[1]) - max(abs(v) for v in fe_strain)
        if drift > 0.005:
            print("NOTE: local element strain exceeds nominal by %.2f%% - "
                  "strain localisation in the specimen (structural effect); "
                  "nominal channel used for the comparison." % (drift * 100))

    # Area-corrected objective: the surrogate is calibrated so its TRUE
    # (Cauchy) stress output matches the ENGINEERING experimental data, but
    # the FE nominal channel is F/A0 = sigma_true/(1+eps) under nlgeom. So
    # recover the FE's true stress (sigma_eng*(1+eps)) and score it against
    # the raw engineering experiment -- that reproduces the calibration's
    # own convention and isolates the finite-strain area effect. (Converting
    # BOTH sides is a no-op: same factor at matched strain.)
    obj_true = None
    if nominal is not None:
        fe_true = [s * (1.0 + e) for s, e in zip(fe_stress, fe_strain)]
        obj_true = fe_objective(fe_time, fe_strain, fe_true, exp,
                                target, meta["weight_reversals"],
                                meta["period_s"])
    peak_t = max(fe_stress)
    peak_c = min(fe_stress)

    csv_path = os.path.join(out_dir, "fe_comparison.csv")
    e_exp_i, s_exp_i = exp_onto_fe(fe_time, fe_strain, exp, target,
                                   meta["period_s"])
    write_csv(csv_path, fe_time, fe_strain, fe_stress, e_exp_i, s_exp_i)
    print("Wrote %s" % csv_path)

    sy0 = meta["sigma_y0"]
    annotation = ("Objective: %.1f MPa | Peak T: %+.0f MPa | "
                  "Peak C: %+.0f MPa | sigma_y0: %s MPa"
                  % (obj, peak_t, peak_c,
                     "%.0f" % sy0 if sy0 else "?"))
    png_path = os.path.join(out_dir, "fe_comparison.png")
    meta["channel"] = "nominal" if nominal is not None else "local element"
    png_ok = try_png(png_path, fe_strain, fe_stress, exp, target, meta,
                     annotation)
    if png_ok:
        print("Wrote %s" % png_path)

    surr = meta["surrogate_objective_MPa"]
    summary = {
        "channel": channel,
        "burnin_cycles": meta["burnin_cycles"],
        "fe_objective_local_MPa": obj_local,
        "fe_objective_truestress_MPa": obj_true,
        "fe_objective_MPa": obj,
        "fe_peak_tension_MPa": peak_t,
        "fe_peak_compression_MPa": peak_c,
        "surrogate_objective_MPa": surr,
        "difference_MPa": (obj - surr) if surr is not None else None,
        "sigma_y0": sy0,
        "model": meta["model"],
        "strain_pct": meta["strain_pct"],
        "diameter_mm": meta["diameter_mm"],
        "gauge_length_mm": meta["gauge_length_mm"],
        "n_fit": int(n_fit),
        "odb": os.path.abspath(odb_path),
        "exp_csv": os.path.abspath(exp_csv),
        "csv": csv_path,
        "png": png_path if png_ok else None,
    }
    json_path = os.path.join(out_dir, "fe_comparison.json")
    f = open(json_path, "w")
    try:
        json.dump(summary, f, indent=1)
    finally:
        f.close()
    print("Wrote %s" % json_path)

    print("FE peak tension:    %+.1f MPa" % peak_t)
    print("FE peak compression: %+.1f MPa" % peak_c)
    if surr is not None:
        print("Surrogate objective: %.1f MPa (from transfer_log.json)" % surr)
    else:
        print("Surrogate objective: n/a (no transfer_log.json found)")
    print("FE objective:        %.1f MPa  [%s]" % (obj, channel))
    if obj_true is not None:
        print("FE objective (area-corrected, true stress): %.1f MPa "
              "(geometric share of the gap: %.1f MPa)"
              % (obj_true, obj - obj_true))
    if obj_local is not None:
        print("FE objective (local element, diagnostic): %.1f MPa"
              % obj_local)
    if surr is not None:
        print("Difference:          %.1f MPa" % (obj - surr))
    print("COMPARISON_OK")


if __name__ == "__main__":
    main()

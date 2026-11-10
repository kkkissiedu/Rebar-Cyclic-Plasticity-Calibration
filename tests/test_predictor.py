"""Regression and sanity tests for the V7 predictor (calibration_app.predictor)."""

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from calibration_app.predictor import apparent_modulus, load_params, predict_bar  # noqa: E402


def test_reference_values():
    # values from the validated build, 16 mm, L/D 8, 3 % strain (checked equal to the research code to 3 decimals)
    r = predict_bar(16, 8, 3)
    assert math.isclose(r["life"]["cycles_to_failure"], 18.1489, rel_tol=1e-3)
    assert math.isclose(r["cycle2_peaks"]["tension_mpa"], 558.08, rel_tol=1e-3)
    assert math.isclose(r["cycle2_peaks"]["compression_mpa"], 418.21, rel_tol=1e-3)
    assert math.isclose(r["apparent_modulus_gpa"], 80.0, abs_tol=0.5)
    assert not r["warnings"]


def test_loops_match_peaks_and_decay():
    r = predict_bar(16, 8, 3)
    l2 = r["loops"][2]
    assert math.isclose(l2["tension_peak"], r["cycle2_peaks"]["tension_mpa"], rel_tol=0.01)
    t = [r["loops"][c]["tension_peak"] for c in sorted(r["loops"])]
    assert all(a >= b for a, b in zip(t, t[1:]))              # peaks fall with cycling
    assert np.all(np.isfinite(l2["stress"]))


def test_life_falls_with_amplitude_and_slenderness():
    lives = [predict_bar(12, 8, a)["life"]["cycles_to_failure"] for a in (1, 2, 3, 4, 6)]
    assert all(a > b for a, b in zip(lives, lives[1:]))
    assert predict_bar(12, 5, 3)["life"]["cycles_to_failure"] > predict_bar(12, 15, 3)["life"]["cycles_to_failure"]


def test_factors_scale():
    base = predict_bar(12, 8, 3)
    f = predict_bar(12, 8, 3, factors=(1.1, 2.0))
    assert f["life"]["cycles_to_failure"] > base["life"]["cycles_to_failure"]
    assert math.isclose(f["cycle2_peaks"]["tension_mpa"], 1.1 * base["cycle2_peaks"]["tension_mpa"], rel_tol=1e-6)


def test_warnings_outside_range():
    assert any("diameter" in s for s in predict_bar(25, 8, 3)["warnings"])
    assert any("fy" in s for s in predict_bar(12, 8, 3, fy_mpa=600)["warnings"])
    assert any("Extensometer" in s for s in predict_bar(12, 8, 3, strain_basis="extensometer")["warnings"])


def test_compliance_formula():
    assert math.isclose(apparent_modulus(16, 5, 1e30), 200000.0, rel_tol=1e-9)
    assert apparent_modulus(16, 5, 1e5) < apparent_modulus(16, 15, 1e5) < 200000.0


def test_export_csv_and_png(tmp_path=None):
    import csv
    import tempfile
    from calibration_app.predictor.export import export_csv, export_png
    d = str(tmp_path) if tmp_path else tempfile.mkdtemp()
    r = predict_bar(12, 8, 3, cycles=[2, 10])
    p = export_csv(r, os.path.join(d, "x.csv"))
    rows = [x for x in csv.reader(open(p, encoding="utf-8")) if x and not x[0].startswith("#")]
    assert rows[0] == ["Time", "Cycle", "Position mm", "Load kN", "Strain %", "Stress MPa"]
    assert len(rows) == 1 + 2 * len(r["loops"][2]["strain"])
    gauge, area = 12 * 8.0, math.pi * 12 ** 2 / 4
    row = rows[50]
    assert math.isclose(float(row[2]), float(row[4]) / 100 * gauge, rel_tol=1e-3, abs_tol=1e-4)
    assert math.isclose(float(row[3]), float(row[5]) * area / 1000, rel_tol=1e-3, abs_tol=1e-4)
    for f in export_png(r, os.path.join(d, "x.png")):
        assert os.path.getsize(f) > 5000


def test_scrap_presets_and_factor_warning():
    import json
    p = json.load(open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "calibration_app", "predictor", "scrap_presets.json"),
                       encoding="utf-8"))
    assert "HYPOTHESIS" in p["status"] and len(p["presets"]) == 4
    for q in p["presets"]:
        assert math.isclose(q["f_s_rm"], q["rm_mean_mpa"] / p["anchor"]["fu_mpa"], abs_tol=2e-3)
        r = predict_bar(12, 5, 3, factors=(q["f_s_rm"], q["f_e"]))
        assert any("untested" in w for w in r["warnings"])
    assert not any("untested" in w for w in predict_bar(12, 5, 3)["warnings"])


def test_params_file_complete():
    p = load_params()
    for k in ("life", "peak2", "decay", "core", "buckling", "constants"):
        assert k in p


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)

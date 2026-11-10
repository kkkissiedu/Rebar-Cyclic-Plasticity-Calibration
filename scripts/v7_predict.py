"""
v7_predict.py
=============

The V7 framework as one object: fit every layer on a list of training files, then predict any file's life, cycle-2 peaks and loops, optionally with
the two material factors of the transfer step (V7_PROTOCOL.md sections 3 and 6, Amendments 1 and 2).

    models = fit_models(train_files)
    f_e, f_s = calib_factors(models, calibration_files)        # from a few tests of a new bar
    out = predict_file(models, file, f_s=f_s, f_e=f_e)

Material factors: f_e multiplies the fatigue ductility coefficient eps_f of the Koh-Stephens relation (ln eps = c + a ln 2Nf, so ln c -> ln c + ln f_e);
f_s multiplies the predicted cycle-2 peak stresses (and with them every later peak, because the degradation curve is normalised by cycle 2).
Loops are the expected loop over the log-normal life distribution (Amendment 2), with the scatter taken from the life model's leave-one-out residuals.
"""

from __future__ import annotations

import math

import numpy as np

import v7_buckle as B
import v7_decay as DC
import v7_life as L1
import v7_loop as V
from v7_loops_cv import fit_epp, loo_log_rms

ZQ, WQ = np.polynomial.hermite.hermgauss(9)


def fit_models(train: list) -> dict:
    rtr = np.array([(f["d"], f["ld"], f["n"], f["nf"]) for f in train], float)
    lv = L1.select_variant(rtr)
    dv = DC.select_decay(train)
    pv = DC.select_p2(train)
    short_ld = min(f["ld"] for f in train)
    core = V.fit_joint([f for f in train if f["ld"] == short_ld], max_nfev=150)
    return {"lv": lv, "lc": L1.fit(rtr, lv), "dv": dv, "th_t": DC.fit_decay(train, "yt", dv), "th_k": DC.fit_decay(train, "yk", dv),
            "pv": pv, "p2c": DC.fit_p2(train, pv), "core": core, "km": math.exp(core[0]),
            "params": B.fit(train, core, mode="joint"), "s_log": loo_log_rms(rtr, lv)}


def _ca(m: dict, f: dict):
    rt = np.array([[f["d"], f["ld"], f["n"], f["nf"]]], float)
    x = L1._x(rt[:, 0], rt[:, 1], m["lv"])
    k = x.shape[1]
    return float((x @ m["lc"][:k])[0]), float(np.minimum(x @ m["lc"][k:], L1.ALPHA_MAX)[0])


def life_pred(m: dict, f: dict, f_e: float = 1.0) -> float:
    c, a = _ca(m, f)
    return float(math.exp((math.log(f["n"] / 100.0) - c - math.log(f_e)) / a) / 2.0)


def calib_factors(m: dict, calib: list) -> tuple:
    """(f_e, f_s) from a few tests of the new bar: mean log offsets of the life relation and of the cycle-2 peaks."""
    ln_e, ln_s = [], []
    for f in calib:
        c, a = _ca(m, f)
        ln_e.append(math.log(f["n"] / 100.0) - c - a * math.log(2.0 * f["nf"]))
        t2p, k2p = DC.pred_p2(m["p2c"], [f], m["pv"])
        ln_s += [math.log(f["t2"] / t2p[0]), math.log(f["k2"] / k2p[0])]
    return math.exp(float(np.mean(ln_e))), math.exp(float(np.mean(ln_s)))


def predict_file(m: dict, f: dict, f_s: float = 1.0, f_e: float = 1.0) -> dict:
    nfp = life_pred(m, f, f_e)
    t2p, k2p = (float(a[0]) * f_s for a in DC.pred_p2(m["p2c"], [f], m["pv"]))
    rows = []
    for c, lp in f["loops"].items():
        one = dict(f, c=np.array([float(c)]))
        spec = {"d": f["d"], "ld": f["ld"], "n": f["n"], "nf": nfp, "cycle": c, "loop": lp}
        avg = np.zeros(lp["s"].shape, float)
        for z, w in zip(ZQ, WQ):
            nfi = nfp * 10.0 ** (m["s_log"] * math.sqrt(2.0) * z)
            gt = float(DC.pred_g(one, m["th_t"], m["dv"], nfi)[0]); gk = float(DC.pred_g(one, m["th_k"], m["dv"], nfi)[0])
            ti, ki = (t2p, k2p) if c == 2 else (t2p * gt, k2p * gk)
            avg += w / math.sqrt(math.pi) * B.predict(dict(spec, nf=nfi), m["core"], m["params"], {"t": ti, "k": ki})
        b1 = B.predict(dict(spec, cycle=2), m["core"], m["params"], {"t": t2p, "k": k2p})
        eo = lp["w"]
        rows.append({"c": int(c), "rmse_fw": float(np.sqrt(np.mean((avg - lp["s"]) ** 2))),
                     "rmse_b1": float(np.sqrt(np.mean((b1 - lp["s"]) ** 2))),
                     "en_fw": float(abs(abs(np.trapz(avg, lp["e"])) - eo) / eo * 100.0), "en_b1": float(abs(abs(np.trapz(b1, lp["e"])) - eo) / eo * 100.0)})
    return {"nf_pred": nfp, "log_life_err": math.log10(nfp / f["nf"]), "t2_err": abs(t2p - f["t2"]) / f["t2"],
            "k2_err": abs(k2p - f["k2"]) / f["k2"], "loops": rows}

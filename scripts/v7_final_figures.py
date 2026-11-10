"""
v7_final_figures.py
===================

Descriptive figures of the single sealed evaluation (post-evaluation; changes nothing in the framework). Needs V7_UNSEAL=1 because it reads the 16 mm bars.

    set V7_UNSEAL=1 & python v7_final_figures.py

Figure 1: predicted against observed life for the 30 sealed bars, and the median loop error per bar (framework against the no-degradation baseline B1).
Figure 2: observed and predicted loops of six sealed bars at three points of life (cycle 2, about half of life, about 90 % of life).
"""

from __future__ import annotations

import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import v7_buckle as B
import v7_data as D
import v7_decay as DC
import v7_final as F
import v7_loops_cv as LC
import v7_predict as P

OUT = os.path.join(D.PROJ, "figures")


def predicted_loop(m: dict, f: dict, c: int, lp: dict) -> np.ndarray:
    nfp = P.life_pred(m, f)
    t2p, k2p = (float(a[0]) for a in DC.pred_p2(m["p2c"], [f], m["pv"]))
    one = dict(f, c=np.array([float(c)]))
    spec = {"d": f["d"], "ld": f["ld"], "n": f["n"], "nf": nfp, "cycle": c, "loop": lp}
    avg = np.zeros(lp["s"].shape, float)
    for z, w in zip(P.ZQ, P.WQ):
        nfi = nfp * 10.0 ** (m["s_log"] * math.sqrt(2.0) * z)
        gt = float(DC.pred_g(one, m["th_t"], m["dv"], nfi)[0]); gk = float(DC.pred_g(one, m["th_k"], m["dv"], nfi)[0])
        ti, ki = (t2p, k2p) if c == 2 else (t2p * gt, k2p * gk)
        avg += w / math.sqrt(math.pi) * B.predict(dict(spec, nf=nfi), m["core"], m["params"], {"t": ti, "k": ki})
    return avg


def main() -> int:
    dev = LC.files()
    test = F.build_files(D.conditions(sealed=True))
    m = P.fit_models(dev)
    outs = [P.predict_file(m, f) for f in test]
    os.makedirs(OUT, exist_ok=True)

    fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))
    obs = np.array([f["nf"] for f in test]); pred = np.array([o["nf_pred"] for o in outs])
    sc = ax[0].scatter(obs, pred, c=[f["ld"] for f in test], cmap="viridis", s=45)
    lim = [min(obs.min(), pred.min()) * 0.7, max(obs.max(), pred.max()) * 1.4]
    ax[0].plot(lim, lim, "k-", lw=1); ax[0].plot(lim, [2 * x for x in lim], "k--", lw=0.8); ax[0].plot(lim, [x / 2 for x in lim], "k--", lw=0.8)
    ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlim(lim); ax[0].set_ylim(lim)
    ax[0].set_xlabel("observed cycles to failure"); ax[0].set_ylabel("predicted cycles to failure"); ax[0].set_title("Life of the sealed 16 mm bars (dashed: factor of 2)")
    fig.colorbar(sc, ax=ax[0], label="L/D"); ax[0].grid(alpha=0.25, which="both")
    order = sorted(range(len(test)), key=lambda i: (test[i]["ld"], test[i]["n"]))
    fw = [np.median([r["rmse_fw"] for r in outs[i]["loops"]]) for i in order]
    b1 = [np.median([r["rmse_b1"] for r in outs[i]["loops"]]) for i in order]
    x = np.arange(len(order))
    ax[1].bar(x - 0.2, fw, 0.4, label="framework", color="#0072B2"); ax[1].bar(x + 0.2, b1, 0.4, label="no degradation (B1)", color="#D55E00")
    ax[1].set_xticks(x); ax[1].set_xticklabels([f"{test[i]['ld']}/{test[i]['n']}%" for i in order], rotation=90, fontsize=7)
    ax[1].set_ylabel("median loop RMSE (MPa)"); ax[1].set_title("Loop error per sealed bar (L/D / strain amplitude)"); ax[1].legend(); ax[1].grid(axis="y", alpha=0.25)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_v7_sealed_summary.png"), dpi=130); plt.close(fig)

    want = [(5, 2), (5, 4), (8, 1), (10, 3), (12, 5), (15, 3)]
    fig, axs = plt.subplots(2, 3, figsize=(16, 9.5))
    cols = {"early": "#009E73", "mid": "#0072B2", "late": "#D55E00"}
    for a, (ld, n) in zip(axs.ravel(), want):
        f = next(x for x in test if x["ld"] == ld and x["n"] == n)
        cs = sorted(f["loops"])
        picks = {"early": 2, "mid": min(cs, key=lambda c: abs(c / f["nf"] - 0.5)), "late": min(cs, key=lambda c: abs(c / f["nf"] - 0.9))}
        for lab, c in picks.items():
            lp = f["loops"][c]
            a.plot(lp["e"] * 100, lp["s"], color="k", lw=1.1, alpha=0.75)
            a.plot(lp["e"] * 100, predicted_loop(m, f, c, lp), color=cols[lab], lw=1.1, label=f"cycle {c} (predicted)")
        a.set_title(f"16 mm, L/D {ld}, {n} % strain, observed life {f['nf']} cycles", fontsize=10)
        a.set_xlabel("engineering strain (%)"); a.set_ylabel("stress (MPa)"); a.grid(alpha=0.25); a.legend(fontsize=7)
    fig.suptitle("Sealed 16 mm bars: observed loops (black) and framework predictions (colour); no data of these bars was used", fontsize=12)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_v7_sealed_loops.png"), dpi=120)
    print("figures saved")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

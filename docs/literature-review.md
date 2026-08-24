# Literature Review - Model & Algorithm Selection for Rebar LCF Calibration

**Prepared:** 2026-07-02, as the evidence base for the model/algorithm choices
implemented in this repository (see `README.md`). Companion dataset summary:
`kashani-2019-dataset.md`.

**Method:** independent literature queries, cross-checked against each other
before comparison to the Kashani (2019) dataset. Every load-bearing claim
carries a confidence tag (HIGH / MEDIUM / LOW); nothing is cited that a search
did not return, and gaps are marked "not found" rather than filled in.

---

## 1. Executive summary

1. **Keep the Chaboche model, but move from 2 → 3 backstresses.** ABAQUS
   natively supports *only* classic combined-hardening Chaboche; every
   competing model (Ohno-Wang, Abdel-Karim-Ohno, Delobelle, McDowell,
   Bari-Hassan) requires a UMAT. The broad steel-Chaboche literature converges
   on **3 backstresses as the minimum-sufficient default** for wide-range cyclic
   strain. **Confidence: HIGH for general steel; MEDIUM-HIGH transferred to
   rebar** (no rebar-specific 2-vs-3-vs-4 study exists).
2. **The Stage-1 grid search must be replaced.** At 6 parameters the grid is
   already 15,625 points; at 3 backstresses (8–9 free parameters) a grid is
   combinatorially infeasible (N⁹). Sobol/LHS space-filling + DE (or CMA-ES) is
   the evidence-backed replacement. **Confidence: HIGH.**
3. **Parameter bounds should be reset to a literature envelope**, replacing the
   app's trial-and-error bounds and its duplicated/contradictory
   `physically_valid()` gate. **Confidence: HIGH for the envelope shape; LOW for
   B500C-specific point values** (no directly-calibrated B500C Chaboche set was
   found in full text).
4. **Kashani 2019 remains an adequate dataset.** Two newer candidates are
   flagged (Moodley et al. 2026 - tests B500C directly; Egger et al. 2021 - 1–5%
   amplitude, cleaner strain data) but **no switch is recommended**.
   **Confidence: MEDIUM-HIGH.**
5. **A ~67 MPa tension/compression asymmetry in the data is unmodelable by any
   symmetric Chaboche model** (2, 3, or 4 backstresses). This is a genuine model
   limitation, already documented in the project README. **Confidence: HIGH.**

---

## 2. Cyclic plasticity model landscape (Q1A)

| Model | Free params | ABAQUS native? | Used for rebar? | Note | Conf. |
|---|---|---|---|---|---|
| Chaboche combined (2/3/4 backstress) | 2/backstress + Q∞,b,σ₀ | **YES** (`*PLASTIC, HARDENING=COMBINED`) | TMT/Fe500 rebar via UMAT comparison only | Baseline; symmetric; over-predicts ratcheting | HIGH |
| Ohno-Wang (1993) | ~Chaboche + 1 exponent/term | NO (UMAT) | Yes (TMT rebar, 1 study) | Better loop shape/area for TMT rebar (single study) | MED |
| Abdel-Karim-Ohno (2000) | OW + μᵢ/term | NO (UMAT) | Not found | Good multiaxial, worse uniaxial | HIGH |
| Delobelle (1995) | ~2-backstress + coupling | NO (UMAT) | Not found | **Uniaxially transparent** → no benefit here | HIGH |
| McDowell (1995) | ~9–12 | NO (UMAT) | Not found | Overestimates multiaxial ratcheting | MED |
| Bari-Hassan (2000/02) | ~Chaboche + 1 | NO (UMAT) | Not found | Multiaxial only | HIGH |
| **Updated Voce-Chaboche (UVC)** (Hartloper/Lignos 2021) | 3 backstress + updated Voce | NO (open-source UMAT) | Structural steel (not rebar) | Fixes yield-plateau + 10–30% yield underestimate | HIGH |
| Memory-surface Chaboche (2025) | Chaboche + surface | NO (UMAT) | Structural steel | Fixes non-saturating isotropic term | HIGH |
| Rebar uniaxial laws (Steel02, Hysteretic, ReinforcingSteel + fatigue) | varies | N/A (OpenSees) | **Yes - dominant rebar practice** | Not 3-D; not ABAQUS; different tool | HIGH |

**Key facts.**
- *Only classic Chaboche is native to ABAQUS.* Switching model means writing/
  validating a UMAT - significant cost and a new verification burden. (HIGH)
- The rebar-specific literature (Kashani, Tripathi & Dhakal, Moodley et al.)
  overwhelmingly uses **OpenSees uniaxial laws, not 3-D Chaboche.** So there is
  no published B500C Chaboche calibration to benchmark against - this project
  sits at the intersection of two mature-but-separate literatures. (HIGH)
- **UVC is the most credible alternative** (open-source UMAT, addresses the
  yield-plateau typical of hot-rolled rebar, reported convergence ≥ ABAQUS
  built-in). But it is still a UMAT, and the existing pure-Python surrogate +
  the existing `FATIGUE.cae` are built on native Chaboche. (HIGH)

**Citations (selected):** Chaboche 1986 (doi:10.1016/0749-6419(86)90010-0);
Ohno & Wang 1993 (doi:10.1016/0749-6419(93)90042-o); Abdel-Karim & Ohno 2000
(doi:10.1016/s0749-6419(99)00052-2); Bari & Hassan 2000/02
(S0749641901000122); Hartloper, de Castro e Sousa & Lignos 2021, UVC
(doi:10.1061/(ASCE)ST.1943-541X.0002964); de Castro e Sousa, Suzuki & Lignos
2020, inverse-problem consistency (doi:10.1061/(ASCE)EM.1943-7889.0001839);
Bakkar et al. 2020 TMT rebar Chaboche/OW (doi:10.1007/s11668-020-00911-z).

---

## 3. Number of backstresses (Q1B)

| Count | Evidence for wide-range (1–6%) rebar LCF | Conf. |
|---|---|---|
| 2 | Works for narrow-range / bearing-steel ratcheting (Koo et al. 2019); likely under-fits the full 1–6% strain span | LOW-MED |
| **3** | **Consensus minimum-sufficient**: Chaboche 1986, ANSYS docs, ASTRJ 2024 review ("above three … does not improve the fitting"), Du et al. 2022, Krolo et al. 2016, Song et al. 2025 | **HIGH** (general steel) / MED-HIGH (rebar) |
| 4 | Only when near-elastic-limit accuracy across many stabilized cycles matters (Santus et al. 2024) | MED (situational) |

**Resolution of the key contradiction (2 vs 3).** The current app and
`FATIGUE.cae` use **2** backstresses. The literature favours **3** because each
backstress maps to a physical region of the loop (fast-saturating near yield;
transient nonlinear mid-loop; near-linear at high strain) - and the project's
new scope (**all diameters × all L/D × 1–6% amplitude**) is precisely the
wide-strain-range regime where the 3rd backstress earns its place. The 2-vs-3
advantage is *not* mainly about the cycle-2-to-6 hysteresis of a single
amplitude (where 2 can suffice), so this is a judgement call, not a proven
necessity for rebar. **Recommendation: default to 3, keep 2 as a fast option,
expose 4 as an optional refinement.**

---

## 4. Calibration algorithms (Q1C)

| Method | Cost @ ~9–11 params | Python | Rebar use | Conf. |
|---|---|---|---|---|
| Grid (brute) | ~2,500+ combos; **N⁹ infeasible at 3 backstresses** | `scipy.optimize.brute` | none | HIGH |
| Differential Evolution | pop 60–150 × 1000+ iters; hours | `scipy…differential_evolution` | none | HIGH |
| Sobol/LHS + Nelder-Mead/Powell | space-filling DOE + local refine (used for 11-var problems) | `scipy.stats.qmc.Sobol` + `minimize` | none | MED-HIGH |
| Bayesian opt | ~120–400 evals | GPyOpt, Optuna+BoTorch | none (closest: Do & Ohsaki 2022, structural steel) | HIGH |
| Genetic algorithm | 11-param Chaboche: 37 iters / ~18.5 h (Dvoršek 2023) | PyGAD | none | HIGH |
| Gradient/adjoint | fewer evals, init-sensitive | `least_squares`, OpenTURNS | none; adjoint-Chaboche literature thin | MED |
| PINN / NN surrogate | ~100× throughput over FE-in-loop DE (Morand 2024) | bespoke PyTorch/TF | none | HIGH (exists) / MED (applicability) |

**Findings.** DE "reliably finds the global optimum" where Nelder-Mead alone
gets trapped (Dorward et al. 2024, doi:10.1016/j.matdes.2024.113409) - so the
app's DE stage is well-chosen. The **grid stage is the weak link**: it does not
scale to 3 backstresses. Sobol/LHS sampling gives the same "seed the global
search" role at any dimensionality. Bayesian optimisation and CMA-ES are strong
alternatives but add dependencies; the surrogate is fast (pure-Python
integrator, ms-scale), so sample-hungry DE remains practical here. No
rebar-specific Chaboche calibration used *any* of these algorithms - a genuine
gap this project would help close.

---

## 5. Parameter bounds envelope (Q1D)

No directly-calibrated **B500C** classical-Chaboche parameter set was found in
full text (LOW confidence on B500C-specific values). The envelope below is drawn
from structural-steel-class calibrations (σ_y 270–600 MPa) and is **HIGH
confidence for order-of-magnitude bounds**, LOW for material specificity.

| Param | Min observed | Max observed | Representative source |
|---|---|---|---|
| C1 (MPa) | 2,303 (SA508) | 28,528 (S355J2+N, UVC) | Krolo 2016; Hartloper 2021 |
| γ1 | 16 (SA508) | 765 (S275) | Krolo 2016; Song 2025 |
| C2 (MPa) | 800 (Q235) | 4,240 (S275) | Wang 2021; Krolo 2016 |
| γ2 | 2.1 (SA508) | 315 (UVC) | Song 2025; Hartloper 2021 |
| C3 (MPa) | 1,120 (S355) | 1,999 (SA508) | Krolo 2016; Song 2025 |
| γ3 | 10 (S355) | 72 (SA508) | Krolo 2016 |
| Q∞ (MPa) | 20.8 (S355) | 228 (A500 HSS) | Krolo 2016; Hartloper 2021 |
| b | 3.2 (S355) | 40 (Q235) | Krolo 2016; Wang 2021 |
| E (GPa) | 185 | 207 | consistent with fixed E=200 GPa |
| σ_y0 (MPa) | - | - | B500C ≈ 500–575 (EN 10080 / EN 1992-1-1) |

**Two full-text papers likely contain true rebar values** (flagged, not
fabricated): Bakkar et al. 2020 (TMT/Fe500, Chaboche UMAT) and Zhu et al. 2024
(HRB400/HTRB600) - obtain in full text for exact numbers.

**Note vs the inspected model.** `FATIGUE.cae` currently holds
C1=2122, γ1=1.097, C2=62358, γ2=138, Q∞=−502, b=18. **γ1≈1.1 and C2≈62 GPa sit
outside the envelope above**, and the implied kinematic saturation
(C1/γ1 + C2/γ2 ≈ 2386 MPa) is physically very large - evidence these are
placeholder/loosely-fit values, reinforcing the need for a literature-anchored
bounds/realism gate.

---

## 6. Dataset assessment (Q1E)

Kashani, Cai, Davis & Vardanega (2019), *ASCE J. Mater. Civ. Eng.* 31(4),
doi:10.1061/(ASCE)MT.1943-5533.0002637 - 120 LCF tests, 4 diameters × 5 L/D,
open Bristol dataset (doi:10.5523/bris.1kz5015zjoel92ueb97kwxd4ps). Adequate and
already in hand.

Flagged (not recommended switches): **Moodley, De Risi & Afshan (2026)**,
J. Building Eng., doi:10.1016/j.jobe.2026.115378 - tests **B500C directly**, same
diameter/L/D matrix, higher IF; **Egger, Rojas & Massone (2021)**,
doi:10.1186/s40069-021-00474-9 - 1–5% amplitude, RGB-photogrammetry strain
(cleaner time-series), open access.

**Open gap:** the "ordered time-series vs sorted" property could not be
confirmed for *any* dataset from metadata alone - including Kashani's. The app's
`data_loader` already reconstructs strain from ordered `Position mm` and excludes
cycle 1 as machine ramp, which is the correct handling; this should be preserved.

---

## 7. Critical appraisal

Applying the peer-review discipline (prioritise the few decisive issues; mark
what is confirmed vs weak; resolve contradictions):

1. **Source quality is uneven and honestly disclosed.** The strongest
   backstress-count and algorithm evidence is from reputable venues (Int. J.
   Plasticity, J. Struct. Eng., Materials, Int. J. Fatigue). The rebar-specific
   evidence is weakest: 3 of 4 attempted full-text fetches were paywalled, and
   the closest rebar-Chaboche source (Bakkar 2020) is abstract-level only.
   *This does not invalidate the recommendation* - it means the 3-backstress
   default is inherited from general steel practice, which is appropriate given
   the absence of rebar-specific data, but it should be stated as such.
2. **The main contradiction (2 vs 3) is resolvable.** Koo et al. (2 backstress)
   is a hardened bearing steel under stress-controlled ratcheting - a different
   regime. It does not undermine the 3-backstress default for strain-controlled
   wide-amplitude LCF. No source argues *against* 3 for this regime.
3. **Is 2-backstress Chaboche actually insufficient for B500C at these
   amplitudes?** *No direct evidence that it is.* For a single amplitude's
   cycles 2–6, a 2-backstress fit can be adequate (the existing app achieves
   usable RMSE). The case for 3 is strongest **because the scope is now all
   amplitudes 1–6% at once**, where one backstress must cover the near-linear
   high-strain branch. This is a well-motivated upgrade, not a correction of a
   broken model. Confidence: MEDIUM-HIGH.
4. **The cited bounds are physically consistent with B500C** at the
   envelope level (E≈200 GPa, σ_y≈500 MPa, Q∞ tens–hundreds of MPa, b single-to
   -tens). They are *not* B500C-validated point values. The realism gate must be
   an envelope + a data-driven check (e.g., σ_y0 not exceeding measured peak
   stress), not hard-coded guesses.
5. **Unresolved external dependencies** (flagged, not blocking): full-text of
   Bakkar 2020 and Zhu 2024 for real rebar Chaboche numbers; license terms of
   the Bristol and Soton datasets; direct inspection of raw-file time ordering.

**Overall:** the evidence supports an *incremental, native-ABAQUS* upgrade
(2→3 backstresses + a scalable search + literature bounds), not a model switch.
A model switch (UVC/Ohno-Wang) is defensible only if yield-plateau fidelity or
loop-shape accuracy proves inadequate after the 3-backstress fit - and it
carries a UMAT development + verification cost.

---

## 8. Recommendation

- **Model:** stay with native Chaboche; **upgrade 2 → 3 backstresses** (default),
  keep 2 (fast) and 4 (near-yield refinement) selectable. Reconsider UVC (UMAT)
  only if the 3-backstress fit is inadequate. *Confidence: MEDIUM-HIGH.*
- **Algorithm:** replace the Stage-1 grid with **Sobol/LHS sampling**; keep
  **Differential Evolution** for Stage-2 refinement (optionally CMA-ES). Grid does
  not scale past 2 backstresses. *Confidence: HIGH.*
- **Bounds:** adopt the §5 literature envelope as named, cited constants; add a
  startup data-consistency check. *Confidence: HIGH (envelope) / LOW (B500C).* 
- **Dataset:** keep Kashani; note Moodley 2026 & Egger 2021 as future options.
  *Confidence: MEDIUM-HIGH.*
- **Known limitation to carry forward:** ~67 MPa T/C asymmetry is unmodelable by
  symmetric Chaboche; document, do not chase. *Confidence: HIGH.*

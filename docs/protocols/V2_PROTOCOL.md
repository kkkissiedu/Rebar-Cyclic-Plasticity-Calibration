# v2 study protocol (written before any v2 result was computed)

Date written: 2026-10-19. Status: frozen. Any change after results are seen must be
listed in a "Deviations" section at the end of this file, with the reason.

## Why a v2
The first ablation (see Method_Document_Monotonic_Calibration.docx) tuned the prior r,
the modulus treatment and the gamma floor on the same 12 bars that were then scored, used a
non-standard reversal weighting, and compared against an ABAQUS coupon whose structural
softening confounds the calibration. v2 fixes these in the surrogate study. ABAQUS is out of
scope for v2 headline numbers.

## Data
12 bars: D in {10, 12, 16} mm, LD5, strain amplitude N in {3, 4, 5, 6} %. Files
`input_data/{D}mm dia/LD5/{N}data.csv`. Never modified.

## Strain treatment (replaces the "apparent modulus" workaround)
Series-spring machine compliance: eps = (Position - C*Load)/L0, C in mm/kN, L0 = D*LD,
Load in kN. Model modulus E = 200000 MPa (Kashani et al. 2015 report monotonic E of 191.7 to
204.5 GPa for B500-class bars).
C is derived from the bar's own monotonic ramp: slope E_ramp = LSQ slope over stress 10 to 50 %
of ramp peak (same function as before), then C_own = (1/E_ramp - 1/200000) * 1000 * L0 / A.
Variants (sensitivity only, not primary):
- `global`: mean C_own over all 12 bars.
- `lodo`: mean C_own over the 8 bars of the other two diameters.
Note (stated in advance): for `own` C, the elastic slope of the corrected ramp is 200 GPa by
construction, and for a linear compliance the result is mathematically identical to the old
apparent-modulus approach. Experiment E0 checks this numerically.

## Objective
Primary: unweighted RMSE, mean over scored cycles (objective.score_hysteresis with
weight_reversals=False). Cycle 1 is run as conditioning history (it is part of the real load
history), never scored. The model is simulated over cycles 1..max(scored cycle) contiguously.
Search: Sobol 4096 points + Bayesian TPE 200 trials, seed 42, physics gate as in the app,
default Chaboche bounds. No tuning of search settings after seeing results.

## Experiments
E0 Equivalence check. For 3 bars (10mm/4%, 12mm/5%, 16mm/3%) and N=1: same parameters scored
   (a) E = E_ramp on raw strain, (b) E = 200 GPa on strain corrected with C_own. Pass criterion:
   absolute difference < 0.01 MPa. Reported whatever the outcome.
E1 Standard cyclic baseline. For each bar and N in {1, 2, 3}: calibrate on cycles 2 and 3 with
   sigma_y0 free in [100, 600] MPa and Q_inf, b free (default bounds); validate on cycles 4 to 6
   (cycles that exist). C = C_own. Report calibration and validation error.
E2 Monotonic route. For each bar, N in {1,2,3}, r in {0.3,...,1.0 step 0.1}: calibrate on the
   corrected ramp only, sigma_y0 in [0.995*r*500, r*500], Q_inf in (0, 1e-6). Score cycles 4 to 6
   (primary) and 2 to 6 (secondary).
   Held-out choice of r (primary result): for each held-out diameter and each N, r* is the r
   minimising the mean cycles-4-to-6 error over the bars of the OTHER two diameters; the reported
   error is on the held-out diameter's four bars at r*. Result per N = mean over the 12 bars.
   Also report: oracle r per bar (lower bound on achievable error, not a result), and the
   fixed-r curve.
E3 Reference models (what error does a no-information model get?). Elastic-perfectly-plastic
   (Chaboche N=1 with C1=1e-3, gamma1=1, Q=0), E = 200 GPa, sigma_y0 on a grid 200 to 600 step
   20, best grid value per bar (generous to the reference). Scored cycles 4 to 6.
E4 Sensitivity (secondary): C = global and C = lodo for E2 at the r* from E2; weighted
   objective (weight_reversals=True) for E1/E2; E2 with Q_inf and b free instead of fixed.

## Decision rules (stated in advance)
- "One backstress is enough" is supported only if, in E1 (validation) and in E2 (held-out r),
  the mean N=1 error is within 10 MPa of the mean N=3 error AND N=1 is not worse than N=3 by
  more than 10 MPa on more than 3 of the 12 bars.
- "The monotonic route is useful" requires E2 held-out error to be clearly below the E3 reference
  and to be reported next to the E1 baseline. If E2 is no better than E3 the route is reported as
  not useful in this form.
- Results are reported whether or not they support the claims.

## Deviations
1. E0 result (2026-10-19, run before any E1/E2/E3 result): FAILED the < 0.01 MPa criterion.
   Same parameters scored under the two treatments differ by 6.4 (10mm/4%), 11.4 (12mm/5%) and
   44.3 (16mm/3%) MPa unweighted. The statement in "Strain treatment" that the two are
   "mathematically identical" for a linear compliance was WRONG. Reason: the corrected strain
   subtracts C times the MEASURED load, whereas in the apparent-modulus formulation the series
   spring acts on the MODEL's stress; they coincide only if model stress equals measured stress.
   Consequences: (a) results of the first study (apparent modulus, raw strain) are not directly
   comparable to v2 numbers; (b) v2 uses the corrected-strain (standard data-preprocessing)
   treatment as specified; no other protocol element is changed.

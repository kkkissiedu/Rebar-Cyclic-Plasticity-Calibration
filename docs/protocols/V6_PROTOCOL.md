# V6: rebar re-calibration with the lessons of V3-V5 (frozen 2026-10-20, before any V6 result exists)

## What changes relative to V3, and why (each change is standard practice, none is tuned to this dataset)
1. Cycle 1 (the virgin loading) is part of the objective. The UVC objective integrates over whole protocols
   including the initial monotonic loading; V3 scored only cycles 2-3 and its fits overshoot the virgin ramp
   (at 5 % strain the model reached 643-722 MPa against 561-580 MPa measured).
2. Non-softening isotropic term (Q_inf >= 1e-4), as in the published UVC/JEM fits (positivity constraints).
   A secondary run lets Q_inf be negative (Q = q*sigma0, q in [-1, 1]) so that its effect can be measured.
3. E is a fitted parameter (as in UVC/JEM), not measured from the calibration cycles. At the point level it is the
   lumped stiffness of bar plus machine; the series-spring form follows: bar E = 200 GPa and machine spring
   k_m from 1/k_m = (L0/A0)(1/E_fit - 1/200000). (Point-level model and series spring are the same model
   for a uniform bar; the spring matters only in the coupon.)
4. Data in true stress and true strain (UVC practice): eps_true = ln(1 + eps_eng), sigma_true = sigma_eng*(1 + eps_eng),
   with eps_eng = (Position - Position[0])/L0 and sigma_eng = (Load - Load[0])*1000/A0 (strain origin = start of test,
   not start of cycle 2). L0 = D*5 mm.
5. Strain path refined to increments <= 1e-4 (linear interpolation) for the simulation, with the stress read back at the
   original samples (inverse problems need fine increments, Suchocki 2025; the coarse path changed the stress
   by up to 5 MPa).
6. Objective and metric are the UVC ones (accumulated-strain weighted, normalised): MSE = integral of
   (sigma_model - sigma_test)^2 d eps* / integral d eps*, S likewise with sigma_test^2, phi_bar = 100*sqrt(MSE/S),
   evaluated over the window as one test.
The model is Chaboche N = 1, 2, 3 backstresses + Voce. No gate, no weighting, no prior. Backstresses ordered by gamma.

## Data and windows
12 files input_data/{D}mm dia/LD5/{N}data.csv (D = 10, 12, 16 mm; N = 3, 4, 5, 6 %), columns Cycle, Position mm, Load kN.
Calibration window = all samples of cycles 1-3. Validation window = cycles 4-6 of the same test (the simulation runs
from the first sample through cycle 6); a cycle is scored only if it has >= 90 % of the median sample count of the
file's cycles 2-6 (this drops the truncated cycle 4 of 16 mm/6 %).

## Parameter box
E_fit in [20, 200] GPa (linear); sigma0 in [1, smax] MPa (smax = largest |true stress| of the calibration window);
C in [1, 1e6] MPa, gamma in [0.1, 1e4], b in [0.01, 100] (log10); Q_inf in [1e-4, smax] (variant "nonsoft") or
Q_inf = q*sigma0 with q in [-1, 1] (variant "soft"). Optimiser: Sobol 8192, TPE 500, Nelder-Mead 2000; seeds 1, 2, 3.

## Experiments
V6-A  Calibrate each of the 12 bars on cycles 1-3 for N = 1, 2, 3 and both variants (3 seeds). Report:
      calibration phi_bar (all of cycles 1-3, and the cycle-1 part alone); the error on the virgin ramp (model peak
      stress at the first tension peak against measured); validation on cycles 4-6 as phi_bar and as mean per-cycle
      RMSE in engineering stress (model stress divided by 1 + eps_eng) so that it is comparable to V3 (43.4 / 38.9 /
      39.4 MPa for N = 1, 2, 3 on 11 bars); and an elastic-perfectly-plastic reference with sigma0 and E fitted on
      the same window (validated on cycles 4-6).
V6-B  Cross-amplitude: parameters of each diameter's 3 % bar (the least buckled record) applied to its 4, 5 and
      6 % bars, with E taken from the test bar's own unloading slopes in its cycles 2-3 (as in V3-B); scored on
      cycles 1-3 and 4-6 separately.
V6-C  Coupon check: for 10 and 16 mm bars at 4 and 5 % and N = 1, 2, 3 (12 jobs), the seed-1 "nonsoft" parameters of
      V6-A in the ABAQUS coupon with bar E = 200 GPa and the series spring (k_m from E_fit), same set-up as V4-T,
      scored as V3-F (cycles 4-6, per-cycle tension peaks). Baselines: V3-F coupon 89.9 MPa, V4-T 105.6 MPa
      (8 N = 1/3 jobs), tension-peak decay gap to experiment 12.9 points (V3-F) and 16.2 (V4-T).

## Decision rules (fixed now)
- V6-A is a success at the point level if the "nonsoft" N = 2 validation error (engineering RMSE, mean over the bars
  with a valid validation cycle) is within 25 % of V3's 38.9 MPa (i.e. <= 48.6 MPa) AND the mean virgin-ramp
  overshoot is below 5 % of the measured peak stress. Otherwise it is reported as a point-level cost of the standard
  constraints.
- One backstress is sufficient only if its mean validation error is <= 1.10 times that of N = 3 and it is more than
  20 % worse than N = 3 on at most 3 of the bars (same rule as V3/V5, applied to the variant stated).
- V6-C supports the new calibration if the mean cycles-4-6 coupon error is <= 67.4 MPa (25 % below V3-F) AND the mean
  absolute difference in tension-peak decay to the experiment is below 12.9 points. Partial if only one holds.
- Failed jobs are reported, not re-run with other parameters. All outcomes are reported.

## Not in V6
Buckling in the geometry, imperfections, a series-spring load dependence, other bar lengths.

## Deviations
(none yet)

## Results (2026-10-20)
V6-A/R (216 + 36 fits), cycles 1-3 calibration, cycles 4-6 validation (engineering RMSE, MPa; N = 1/2/3):
- non-softening: 78.2 / 78.8 / 79.1 (validation phi_bar 19.4 / 18.9 / 18.8 %); calibration phi_bar 9.26 / 8.86 / 8.71 %;
  virgin-peak error -6.8 / -4.5 / -3.9 %; box-limit hits 36/36, 35/36, 36/36; gamma collapse 0/0/2.
- softening allowed: 56.3 / 50.5 / 49.8 (phi_bar 12.6 / 11.2 / 10.9 %); virgin-peak error -0.2 / +2.4 / +2.2 %;
  limit hits 33 / 29 / 32 of 36.
- elastic-perfectly-plastic reference: 87.3 MPa.
- Point-level success rule (non-softening N=2 <= 48.6 MPa and virgin overshoot < 5 %): NOT MET (78.8 MPa); reported
  as the point-level cost of the non-softening constraint. Virgin overshoot criterion alone is met (-4.5 %).
- N rule: non-softening N1/N3 = 0.988, N2/N3 = 0.997 (sufficient, but all N equal because the constraint dominates);
  softening N1/N3 = 1.129 with 5 bars > 20 % worse (N=1 NOT sufficient), N2/N3 = 1.013 (N=2 sufficient).
V6-B (3 % bar parameters to 4/5/6 % bars, cycles 4-6 RMSE): non-softening 115.6 / 126.3 / 128.6; softening
92.6 / 106.2 / 120.8 (N = 1 / 2 / 3).
V6-C (12 coupon jobs, all solved; seed-1 non-softening parameters, bar E = 200 GPa + series spring), mean cycles 4-6
coupon error over 10/16 mm x 4/5 %: N=1 63.3, N=2 50.5, N=3 41.8 MPa (V3-F 89.9, V4-T 105.6). Mean tension-peak change
cycle 2 to 6: FE -17.5 / -14.3 / -8.1 % for N=1/2/3 against -13.8 % measured; mean absolute decay gap 9.0 / 8.1 / 9.1 points
(V3-F 12.9, V4-T 16.2). Both criteria (<= 67.4 MPa, gap < 12.9) are met for every N: V6-C SUPPORTS the new calibration.
Caveats: 4 bars; the coupon error (nominal engineering stress) and the point-level error use different metrics and are
not directly comparable; Q_inf and b sit on their box limits in nearly all fits (isotropic term effectively inactive);
validation is on later cycles of the same test; per-bar coupon errors range 33-88 MPa and the 5 % bars are worst.
In the coupon N=1 is 1.52 times N=3 (not sufficient by the 1.10 rule).

## Addendum V6-C2 (2026-10-20, written after V6-C results and before any V6-C2 job was run)
Reason: V6-C used 4 bars and only the non-softening variant; the user asked whether N=1 is really worse.
Jobs (same set-up, seed 1, bar E = 200 GPa + series spring, templates from V3-F input files):
- V6-C2a: non-softening, 12 mm bar at 4 % and 5 %, N = 1, 2, 3 (6 jobs).
- V6-C2b: softening variant (Q_inf = q*sigma0, q in [-1, 1]), 10 and 16 mm at 4 % and 5 %, N = 1, 2, 3 (12 jobs).
The 3 % and 6 % bars have no V3-F template input files and are not run. Decision rules as in V6-C (coupon error <= 67.4 MPa,
decay gap < 12.9 points); N=1 sufficiency judged by the same 1.10 ratio and 3-bar rule. Failed jobs are reported, not re-run.
Result files: calibration_sessions/v6/c.jsonl (variant field) and c_soft.jsonl.
Deviation (2026-10-20): the first 6 V6-C2a jobs (12 mm) failed in the input processor (8 s) because v4_spring_test.edit_inp
hard-coded the lateral top constraint set "_PickedSet5", which the 12 mm deck names "Specimen-1.Set-Top". edit_inp now reads the set
name from the deck (output for the 10 and 16 mm decks is unchanged). The failed records are kept in
calibration_sessions/v6/c_failed_first_attempt_12mm_bc_setname.jsonl; the same 6 jobs are re-run with the corrected script.
Deviation 2 (2026-10-20): the second attempt of the 12 mm jobs solved but is INVALID. The 12 mm deck lists the top-face nodes
explicitly (not "generate"), and edit_inp read the node list as a range, tying only 2 of 515 top nodes, so the bar was barely loaded
(errors 351-366 MPa for every N, stress about a tenth of measured). edit_inp now expands both formats; the 10 and 16 mm output is
byte-identical to before (checked). The invalid records are kept in
calibration_sessions/v6/c_invalid_second_attempt_12mm_top_nodes.jsonl; the 6 jobs are re-run a third time.

## Results V6-C2 (2026-10-20; 12 mm jobs from the third attempt, after both deviations)
V6-C2a, non-softening, 12 mm: coupon error 35.6 / 32.0 / 35.3 MPa (4 %) and 34.3 / 33.6 / 32.1 MPa (5 %) for N = 1 / 2 / 3.
Over all 6 non-softening bars (10, 12, 16 mm x 4, 5 %): mean coupon error N=1 53.9, N=2 44.6, N=3 39.1 MPa; N1/N3 = 1.378
(4 bars > 20 % worse), N2/N3 = 1.141 (1 bar > 20 % worse); mean tension-peak decay gap 6.3 / 6.5 / 7.2 points. Both V6-C criteria are met
for every N. One backstress is NOT sufficient by the frozen rule (ratio 1.378 > 1.10, 4 bars > 3), N=2 narrowly misses it
(1.141 > 1.10). Per bar, N=1 minus N=3: +0.3 and +2 MPa on the 12 mm bars, +12 and +13 on 16 mm, +21 and +40 on 10 mm (4 % and 5 %).
V6-C2b, softening allowed (10, 16 mm x 4, 5 %): mean coupon error N=1 233.1, N=2 113.8, N=3 139.4 MPa; decay gap 46.3 / 17.9 / 28.5 points.
Criteria are not met for any N. The fits with Q_inf near -sigma0 let the yield stress fall toward zero in the coupon (large cyclic softening
in the loops, see figures/v6-coupon-loops-soft.png); the softening variant is not usable in the coupon without a constraint.
Not run: 3 % and 6 % bars (no V3-F template input files).

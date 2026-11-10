# V7 (FROZEN 2026-11-10): a predictive framework for rebar hysteresis across diameter, slenderness, strain amplitude and cycle

Status: frozen on 2026-11-10 after the user confirmed the choices in section 7. Nothing below has been fitted. The only look at the data so far is one exploratory plot on the development set (section 2): the decay of peak stress against life
fraction n/Nf.

## 1. Goal
Predict, for a bar of given diameter D, slenderness L/D, strain amplitude and cycle number: (a) the cycles to failure, (b) the peak tension
and compression stress of every cycle (the degradation), and (c) the loop shape, so that the loops can be reconstructed. The framework is
fitted on the Kashani B500C data and transferred to a new bar by a small number of material factors calibrated from a few tests.

## 2. Data and the sealed final test
120 files: D = 10, 12, 16, 20 mm; L/D = 5, 8, 10, 12, 15; amplitude 1 to 6 %; one specimen each; strain from the actuator LVDT (so it includes
machine compliance; the same machine for every bar); cycles to failure recorded.
- Development set: D = 10, 12, 20 mm (90 files). All model choices are made here, by cross-validation.
- Sealed final test: all D = 16 mm files (30). Not read by any script until the framework and its parameters are frozen. It is then evaluated once.
- Exploratory look already taken (development set only): peak stress, normalised by the cycle-2 value, falls along a common curve against n/Nf that is
  ordered by L/D, with compression below tension. This motivated the architecture below; it is not a result.

## 3. Architecture (hypothesis; variants allowed only as listed)
- L1 Life: eps_a = eps_f * (2 Nf)^alpha (Koh-Stephens form, as in Kashani 2019 and Moodley 2024), with ln eps_f and alpha linear in lambda_p = (L/D) *
  sqrt(fy/100) and in D. Variants: with and without the D term; lambda_p replaced by L/D alone.
- L2 Degradation: tension and compression peak stress of cycle n, normalised by their cycle-2 values, as a master curve g(x; lambda_p, eps_a) of the life
  fraction x = n / Nf. Variants (compared only by cross-validation on the development set): a power term plus a terminal drop; a Weibull-type
  survival form; a sum of two exponentials. The cycle-2 peaks follow from a cyclic stress-strain curve of the amplitude.
- L3 Loop shape: the Chaboche (N = 2) core, with E fitted as in V6 (apparent modulus), gives the unit loop; its stress range is scaled to the predicted
  peaks of L2. Loops are reconstructed by this mapping.
- Transfer: material factors on the cycle-2 peak stress level and on eps_f, calibrated from a few tests of the new bar (section 6). The scale factor
  uses fy in lambda_p.

## 4. Validation on the development set (all run, all reported)
- Leave one condition out (90 folds).
- Leave one diameter out (3), one L/D out (5), one amplitude out (6).
- Leave one pair of factors out: D x L/D (15 cells), D x amplitude (18), L/D x amplitude (30).
- Extrapolation blocks: amplitudes >= 5 % held out; L/D >= 12 held out; D = 20 mm held out (fit on 10 and 12 mm).
- Repeated random splits: 20 % of conditions held out, 200 repeats, fixed seed.
Every fold refits everything (including the choice among variants inside the training data only).

## 5. Metrics and baselines (fixed now)
Metrics, per held-out condition: (i) loop RMSE per cycle (MPa) against the test, (ii) relative error of the tension and compression peak per cycle (%),
(iii) log10(Nf predicted / Nf observed) and the share within a factor of 2, (iv) dissipated energy per cycle error (%).
Baselines: B0 per-condition oracle (the framework's loop model fitted to the condition itself; the floor); B1 no degradation (cycle-2 loop repeated);
B2 life from the Kashani-form regression on lambda_p fitted the same way; B3 elastic-perfectly-plastic with fitted E and yield.
Decision rules (relative, fixed now): the framework is accepted for interpolation if, over the leave-one-condition-out and leave-one-amplitude-out folds,
the median peak-stress error is below that of B1 on at least 80 % of folds and the median |log10 life ratio| is at most 0.30 (a factor of 2, the usual
scatter band for fatigue-life prediction). For the extrapolation blocks the same quantities are reported without a pass/fail rule.

## 6. Transfer test (new-bar claim)
On the sealed 16 mm set after freezing: fit the two material factors from the L/D = 5 bars at 2 % and 4 % only, predict the other 28 files; also
predict with no factors. Reported with the baselines. Additionally, on the development set: treat each diameter in turn as a new bar with the same
two calibration tests.

## 7. Choices confirmed by the user on 2026-11-10
1. The sealed final test is D = 16 mm (all 30 files).
2. The factor-of-2 life criterion and the 80 % of folds criterion in section 5 are accepted.

## 8. Integrity and limits
No result from the sealed set before freezing; deviations logged and failed runs archived; the choice among variants only inside training folds.
Limits stated in advance: strain is LVDT-based with machine compliance (same machine for all bars; transfer to extensometer-based tests is untested);
one specimen per condition, so no scatter estimate; the framework describes constant-amplitude tests to fracture, not variable amplitude.


## Amendment 1 (2026-11-10), written before any further evaluation and before the sealed set is touched
Reason: a loop layer made of the Chaboche core scaled to the peaks (section 3, L3) was tested on development data only and gave loop-shape errors of
90 to 152 MPa on slender bars (L/D 8 to 15) even when the true peaks were supplied, because those loops are strongly distorted by buckling. L3 is therefore
replaced by two parts, both fitted on development data only:
- L3a core: the Chaboche (N = 2) core with the machine-stiffness compliance model, fitted on the training files of the SHORTEST available L/D (L/D 5 unless it
  is held out) inside every fold.
- L3b buckling distortion (v7_buckle.py): a compression branch limited by a power-law capacity envelope (the closed form of a three-hinge mechanism, as in
  Gomes and Appleton 1997), and a tension branch that straightens back to the tension peak; five shape parameters, linear in (1, ln lambda_p, ln amplitude, n/Nf),
  20 coefficients, fitted on the training loops of each fold. The form was chosen on development data by comparing alternatives under leave-one-L/D-out and
  leave-one-condition-out cross-validation; that selection used development data only.
Loop metrics (sections 4 and 5 unchanged otherwise): evaluated on cycle 2 and on the cycles nearest to life fractions 0.1, 0.3, 0.5, 0.7 and 0.9 of each held-out
file; loop RMSE (MPa) against the observed loop, and energy per cycle error (%). Full chain = predicted Nf (L1), predicted cycle-2 peaks and decay (L2).
Baselines: B1 (the cycle-2 loop repeated, no degradation), B3 (elastic-perfectly-plastic with the same compliance model, yield fitted on the training files),
and the core alone scaled to the same predicted peaks. Pass rule fixed now: over the leave-one-condition-out and leave-one-amplitude-out folds, the full-chain
median loop RMSE of the framework is below that of B1 and below that of the core alone on at least 80 % of folds, and its pooled median energy error is below
that of B1. Groups run for L3: leave one condition, diameter, L/D, amplitude out and the three extrapolation blocks; the cell groups and the random splits are
added if time allows and reported either way.


## Amendment 2 (2026-11-10), written after seeing development-set results of Amendment 1 and before the sealed set is touched
Observation (development set only): for L/D 5 bars at 4 to 6 % strain, which live only 4 to 12 cycles, the full-chain loop error on later cycles was 93 to 266 MPa
against 42 to 59 MPa for the no-degradation baseline, while the shape-only error stayed at 22 to 26 MPa. Cause: a small error in predicted life places a late cycle
on the fracture cliff of the degradation curve.
Change: the loop of every cycle is the EXPECTED loop under the predictive distribution of life, not the loop at the point estimate. Life is log-normal (ASTM E739:
log life normally distributed with constant scatter), centred on the L1 prediction, with the standard deviation of log10 life taken as the RMS of the L1 model's own
leave-one-out residuals on the training files of the fold (never a constant chosen by hand). The expectation is computed with 9-point Gauss-Hermite quadrature. The
expected value is the minimiser of squared error, which is what the loop RMSE measures.
Rule fixed now: the averaged version replaces the point estimate only if (i) its median full-chain loop RMSE over the leave-one-condition-out and leave-one-amplitude-out
folds is not higher than the point estimate's, and (ii) for L/D 5 bars at 4 to 6 % strain, cycles after cycle 2, its median loop RMSE is lower. Otherwise the point estimate
is kept and the result is reported. Both versions are reported in every case.


## Amendment 3 (2026-11-10): the framework is frozen; rules for the single evaluation on the sealed set
Frozen framework: L1 life (Koh-Stephens form, variants chosen inside the training data), L2 degradation (three forms, chosen inside the training data), cycle-2 peaks (log-linear
in amplitude), L3a core (Chaboche N = 2 fitted on the shortest available L/D with the machine-stiffness compliance model), L3b buckling distortion (Amendment 1), and loops as the
expected loop over the log-normal life distribution (Amendment 2). The code that defines it is fixed by the SHA-256 hashes recorded in V7_FREEZE.txt.
Evaluation: fit on all 90 development files; predict all 30 sealed files (D = 16 mm) using none of their data. Pass rules, fixed now: (1) median |log10 (predicted / observed life)| is at
most log10 2; (2) the median loop RMSE of the framework is below that of the cycle-2-repeated baseline B1 for at least 80 % of the sealed files; (3) the pooled median energy-per-cycle
error is below that of B1. Reported whatever the outcome, with breakdowns by L/D and amplitude. Separately, and without a pass/fail rule: the material factors calibrated from the
sealed L/D 5 tests at 2 % and 4 % applied to the other 28 files, compared with no factors.
The sealed set is used once. v7_final.py records its use in V7_FREEZE.txt and refuses to run again; any change to the framework after the result is seen, and any rerun, must be reported as
post-final and cannot be presented as an out-of-sample validation. A dry run (the 20 mm development bars standing in) has already proved the code path.


## Results of the single evaluation on the sealed set (2026-11-10, after Amendment 3; frozen framework, fitted on the 90 development files only)
All three pre-set rules are met on the 30 sealed 16 mm bars.
- Life: median |log10 (predicted / observed)| = 0.101 (limit 0.301), 100 % of bars within a factor of 2.
- Loops (cycle 2 and four later cycles per bar): median loop RMSE 39.1 MPa for the framework against 93.7 MPa for the cycle-2-repeated baseline B1; the framework
  is better on 87 % of the bars (needs 80 %).
- Energy per cycle: median error 11.6 % against 34.1 % for B1.
- Cycle-2 peak stress error: 4.5 % in tension, 8.1 % in compression.
By L/D (loop RMSE, framework against B1, MPa): L/D 5: 57 against 68; L/D 8: 35 against 93; L/D 10: 42 against 90; L/D 12: 37 against 115; L/D 15: 30 against 105. By amplitude: 1 %: 50 against 110; 2 %: 41 against 103; 3 %: 42 against 83; 4 %: 34 against 90; 5 %: 36 against 78; 6 %: 32 against 93.
Known limit confirmed on sealed data: for the short L/D 5 bars at 4, 5 and 6 % strain (lives of 4 to 10 cycles) the framework is worse than B1 (about 62 against 41, 150 against 54 and 139 against 49 MPa);
the loop shape itself is adequate there, the degradation curve near fracture is not.
Material factors from the sealed L/D 5 tests at 2 % and 4 % (f_e = 0.948, f_s = 1.004) applied to the other 28 bars: loop RMSE 38.0 without and
39.2 with the factors; life error 0.101 and 0.097. Within one steel the factors are close to 1 and change little; this does not
test transfer to a different steel. Development-set evidence for the factors: a synthetic 1.2 stress change was recovered accurately.
Limits (unchanged): one steel (B500C), one machine, actuator-LVDT strain with machine compliance, one specimen per condition, constant amplitude to fracture, diameters 10 to 20 mm, L/D 5 to 15, amplitudes 1 to 6 %.

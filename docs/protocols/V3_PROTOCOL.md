# V3 protocol: the textbook method, nothing tuned to this dataset

Date written: 2026-10-19. Status: frozen before any V3 result exists. Changes after results are
seen must be listed under "Deviations", with the reason.

## Principle
V3 uses only what a standard calibration of a Chaboche combined-hardening model uses: the
recorded test (strain, stress), the model equations, a generic parameter box and a generic
optimiser. Everything that was specific to this dataset or to our route in V1/V2 is removed:

| Removed from V1/V2 | Reason |
|---|---|
| Monotonic-ramp-only calibration and the Bauschinger prior r | needs an assumed value tuned on the data |
| Hard-coded 500 MPa "monotonic proof stress" | dataset-specific number |
| Assumed E = 200 GPa and the compliance correction | an assumption about the test machine |
| Reversal weighting (1.5x) | project-specific, no literature precedent found |
| Physics gate (gamma-separation floor, saturation cap, yield floor) | project-specific constraints tuned on this data |
| Default bounds taken from other steels (C1 <= 30000 etc.), which were hit as limits | dataset-dependent prior |
| Special treatment of cycle 1 as "burn-in" | replaced by simulating the whole recorded test |
| Selecting scope (LD5) after seeing failures | scope is fixed here in advance |

"No assumptions" cannot be literally true. The irreducible choices are listed here so that nothing
is hidden: (1) E is measured from the test itself (see below); (2) the parameter box is very wide
and dimensionless in shape (log scale), not taken from any steel; (3) which cycles calibrate and
which validate; (4) optimiser settings. None was adjusted after seeing results.

## Data and scope
Kashani et al. 2019 B500C files `input_data/{D}mm dia/LD5/{N}data.csv`, D in {10,12,16} mm,
N in {3,4,5,6} % (12 files). Strain is used exactly as recorded (actuator position / gauge length,
gauge = D*LD, strain zero at the first sample of cycle 2; only increments matter to the model).
Stress = load / original area. No file is modified. Longer bars (LD8-15) are out of scope for V3
(compressive buckling is a structural effect, not a constitutive one); this scope is fixed now.

## Model
Chaboche, N = 1, 2, 3 backstresses (Armstrong-Frederick), isotropic Voce hardening (Q_inf, b),
initial yield sigma0 as a calibrated parameter, exactly the equations in the existing surrogate
(uniaxial radial return). The whole recorded test is simulated from its first sample (cycle 1
included), so no special conditioning step exists.

## Elastic modulus (measured, not assumed)
E is the median of least-squares slopes of stress against strain on the unloading branches
following each strain extremum of the calibration cycles (cycles 2 and 3; 2 extrema each),
using the points where stress has dropped by 5 % to 40 % of that extremum's stress. It is taken
from the calibration cycles only (never from validation cycles) and is held fixed during the
search. This is the usual "E from the experiment" step. Because strain here is actuator-based,
E is an effective stiffness of the test system, and it is reported as such.

## Parameter box (generic, same for every file; no steel-specific values)
Let S = max |stress| over the calibration cycles.
- sigma0 in [1 MPa, S], linear.
- C_k in [1, 1e6] MPa, log10 scale. gamma_k in [0.1, 1e4], log10 scale. b in [0.01, 100], log10.
- Q_inf in [-0.5 S, +S] MPa, linear.
Backstresses are ordered by gamma (descending) before simulation (removes label switching; no
physical constraint is imposed). No other constraint, projection or gate.

## Objective
Unweighted RMSE (MPa) per scored cycle, mean over scored cycles (objective.score_hysteresis with
weight_reversals=False). Data-quality rule, fixed in advance: a cycle is scored only if it contains
at least 90 % of the median number of samples of the file's cycles 2..6.

## Optimiser
Sobol 8192 points (scrambled, seed s) -> Bayesian TPE 500 trials (seed s, best Sobol point
enqueued) -> Nelder-Mead polish 2000 evaluations from the best point (inside the box).
Three seeds s in {1, 2, 3}. The reported value of a bar is the mean over the three seeds;
seed-to-seed SD is reported.

## Experiments
V3-A  Within-test validation. Calibrate on cycles 2-3; validate on cycles 4-6 (the same test,
      later cycles, never seen by the optimiser nor by the E measurement). N = 1, 2, 3; 12 bars.
V3-B  Cross-amplitude validation. For each diameter and every ordered pair (calibration bar a,
      test bar b, a != b) of its four amplitudes: apply the parameters calibrated on bar a
      (V3-A record) to bar b using bar b's own measured E, and score bar b's cycles 2-6.
V3-R  Reference model. Elastic-perfectly-plastic (C_k -> 0, no isotropic hardening) with sigma0
      fitted by the same optimiser on cycles 2-3 only (no use of validation data); validated on
      cycles 4-6. Gives the error of a model with no hardening, to show what the Chaboche terms buy.

## Reporting rules (fixed now)
- Report mean, SD over seeds, and per-bar values; report by diameter (the 12 bars are clustered
  in 3 diameters); never hide a bar. Report with all 12 bars (cycles failing the quality rule are
  dropped by rule, not by hand).
- Parameter-identifiability report: fraction of fits with any parameter within 1 % (log-range) of
  its box limit; gamma collapse (gamma_k/gamma_{k+1} < 1.1) counts.
- "One backstress is enough" is claimed only if, in V3-A, mean N=1 error <= 1.10 x mean N=3 error
  AND N=1 exceeds N=3 by more than 20 % on at most 3 of 12 bars; the same comparison is reported
  for V3-B. Otherwise the result is stated as "N=1 costs X %".
- V3-R is the yardstick: a model that is not clearly better than V3-R is reported as such.
- All failures, boundary hits and non-converged fits are reported.

## V3-F  Full-coupon ABAQUS check (frozen 2026-10-19, before any V3-F job was submitted)
Purpose: does the surrogate-calibrated material, put into the real 3-D coupon, reproduce the
measured loops, and how large is the surrogate-versus-coupon gap?
- Jobs: the frozen V3-A parameter sets of seed 1 for N = 1 and N = 3, bars 10, 12, 16 mm at 4 % and
  5 % strain (12 jobs). No other seed or parameter set; no parameter is changed after seeing a result.
  Failed or non-converged jobs are reported as failures, not re-run with other parameters.
- Material: native ABAQUS Chaboche (combined hardening table sigma0, C_k, gamma_k; cyclic hardening
  sigma0, Q_inf, b), E = the bar's measured V3 E, nu = 0.3.
- Model: existing project coupon (FATIGUE_rebuilt.cae copy; mesh rebuilt for other diameters, C3D8R,
  nlgeom ON, bottom fixed, top displacement-driven), loading = the measured cycle-1 ramp followed by a
  triangular strain wave of the file's nominal amplitude, period 20 s, cycles 2-6. Only change from
  the earlier coupon runs: maximum time increment = period/200 (was period/50), so that the
  amplitude peaks are reached closely (the earlier setting clipped them, e.g. 3.86 % for a 4 % bar).
- Response compared: nominal channel (total reaction force / original area against mean top
  displacement / gauge), as in the project's comparison script.
- Scoring: unweighted per-cycle RMSE (same function as the surrogate). Headline: cycles 4-6 (as
  V3-A). Secondary: cycles 2-6. Reported next to the surrogate value for the same parameters, per
  bar and per cycle, with FE and measured peak stresses per cycle and the FE/measured peak strains.
- Interpretation rule (conventional thresholds, not derived): surrogate and coupon are called
  consistent if the mean absolute difference over the bars (cycles 4-6) is <= 15 MPa and no bar
  differs by more than 30 MPa. Otherwise the gap is reported with its diagnostics (per-cycle peak
  decay, local versus nominal strain) and no accuracy claim is made for the coupon.
- Waves: at most 6 jobs run concurrently (memory).

## Deviations
1. (2026-10-19, after the frozen V3-A/V3-R/V3-B runs) The protocol run found Q_inf at its lower box
   limit (-0.5 S) in 30 of 36 fits for every N, i.e. the data want more cyclic softening than the
   box allowed. The frozen results are reported unchanged and are conditional on that limit.
   A sensitivity run V3-S (added AFTER seeing this, so exploratory, not part of the frozen protocol)
   repeats V3-A with Q_inf = q * sigma0, q in [-1, 1] (the generic physical limit: yield radius
   stays non-negative). Same data, seeds, optimiser and objective; env var V3_QREL=1; output
   calibration_sessions/v3/a_qrel.jsonl. It does not replace the frozen result.

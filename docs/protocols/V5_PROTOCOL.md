# V5: independent benchmark on extensometer coupon data (frozen 2026-10-20, before any V5 result exists)

## Why
All rebar-data results so far are confounded: strain is actuator LVDT (apparent modulus), the bars buckle
and degrade (they never stabilise), and the data come from one dataset. V5 removes those confounds by
using the UVC authors' open true stress-true strain coupon database (extensometer strain, several
protocols per steel; folder data/Clean_Data (set UVC_DATA), from the Zenodo "Database of Uniaxial Cyclic and Tensile
Coupon Tests for Structural Metallic Materials"), for which published fit errors exist. Aims:
(1) check that our calibration implementation reproduces published Voce-Chaboche fit errors;
(2) answer, on proper data, how many backstresses are needed (N = 1, 2, 3), in-sample and out-of-sample.

## Data
Six steel sets, each with the same five protocols (files *_05inc, *_2const, *_2inc, *_random, *_tensile;
columns Sigma_true [MPa], e_true):
A500 (UVC paper ID 9), A992_Gr50/A992_W14X82_web (ID 7), A992_Gr50/A992_W14X82_flange (ID 8),
BCP325 (ID 10), BCR295 (ID 11), HYP400 (ID 12). Files are never modified. All five protocols of a set are
the five used in the UVC paper (LP 1, 2, 5, 6, 9). Other sets in the folder are out of scope for V5.

## Model (standard Voce-Chaboche, no additions)
Chaboche with N backstresses (Armstrong-Frederick) plus Voce isotropic hardening; E and sigma_y0 are
fitted parameters (as in UVC/JEM); uniaxial radial return of the existing surrogate, driven by the recorded
true strain, output true stress. N = 1, 2, 3. No physics gate, no weighting, no yield plateau term.

## Objective and metric (UVC Eqs. 30-33)
Per protocol j, accumulated strain eps*_j = integral |d eps|. MSE_j = integral (sigma_model - sigma_test)^2
d eps* / integral d eps*; S_j = integral sigma_test^2 d eps* / integral d eps* (trapezoidal on samples).
Fitted objective f = sum_j MSE_j. Reported metric phi_bar (%) = 100 * sqrt( sum_j MSE_j / sum_j S_j ).
For a single held-out protocol: phi_bar_j = 100 * sqrt(MSE_j / S_j).

## Parameter box (generic, positivity as in the UVC/JEM papers)
E in [100, 300] GPa (linear); sigma_y0 in [1 MPa, max |sigma| of the fitted protocols]; C_k in [1, 1e6] MPa
and gamma_k in [0.1, 1e4] on a log scale; Q_inf in [1e-4, max |sigma|] (non-softening, as the published
fits); b in [0.01, 100] on a log scale. Backstresses ordered by gamma (label switching only). Nothing else.

## Optimiser
Sobol 8192 points, then Bayesian TPE 500 trials (best Sobol point enqueued), then Nelder-Mead polish 2000
evaluations; seeds 1, 2, 3; the reported value is the mean over seeds, with the seed spread.

## Experiments
V5-A  Joint calibration, all five protocols (in-sample), per steel and N. Compare with the published
      Voce-Chaboche phi_bar of the same sets (UVC paper Table 3, N = 2) and with the published UVC values.
V5-B  Leave-one-protocol-out: fit on four protocols, predict the fifth (each protocol held out in turn),
      per steel and N. Report phi_bar_j on the held-out protocol and the in-sample phi_bar of the four.
V5-C  Identifiability report for V5-A: fraction of fits touching a box limit, gamma collapse
      (gamma_k / gamma_{k+1} < 1.1), and spread of parameters across seeds and across N.

## Published reference values (UVC paper Table 3, phi_bar %, plain Voce-Chaboche with 2 backstresses / UVC)
ID 7 A992 web 7.01 / 5.21; ID 8 A992 flange 7.31 / 4.76; ID 9 A500 6.86 / 4.99; ID 10 BCP325 5.05 / 3.89;
ID 11 BCR295 8.38 / 5.21; ID 12 HYP400 5.31 / 4.67.

## Decision rules (fixed now)
- Implementation check: our V5-A N = 2 phi_bar is "consistent with the published Voce-Chaboche value" if it is
  within +-1.0 percentage point on at least 4 of the 6 steels. If it is worse by more than 1.5 points on any
  steel, that steel is reported as a discrepancy and investigated before any N conclusion is drawn.
- Number of backstresses (for V5-A and, separately, V5-B): N = 1 is called sufficient only if its mean phi_bar
  over the six steels is <= 1.10 times that of N = 3 AND it is more than 20 % worse than N = 3 on at most
  1 of the 6 steels. The same test is applied to N = 2. Otherwise the relative cost is reported.
- Results are reported whether or not they support the claims.

## Not in V5
UVC yield-plateau term, joint calibration across steels, the rebar data, ABAQUS.

## Results (2026-10-20; 54 joint fits, 270 leave-one-protocol-out fits, no deviations from the protocol)
V5-A joint, phi_bar % (mean of 3 seeds): N=1 7.83, N=2 6.67, N=3 6.26 (mean over the six steels).
Implementation check: our N=2 minus published Voce-Chaboche = -0.12, +0.21, -0.34, +0.51, -0.51, +0.36 points
(A992 web, A992 flange, A500, BCP325, BCR295, HYP400); mean 6.67 vs published 6.65. Within +-1.0 point on 6/6
steels: CONSISTENT. Published UVC mean 4.79 (UVC adds a yield-plateau term that we do not have).
Number of backstresses, V5-A: N=1 vs N=3 mean ratio 1.251, N=1 more than 20 % worse on 4 of 6 steels: NOT
sufficient. N=2 vs N=3 ratio 1.066, 0 steels more than 20 % worse: SUFFICIENT.
V5-B held-out phi_bar %: N=1 10.40, N=2 8.55, N=3 8.30 (in-sample on the four fitted protocols 7.54 / 6.27 /
6.10). N=1 vs N=3 ratio 1.253, 4 of 6 steels more than 20 % worse: NOT sufficient. N=2 vs N=3 ratio 1.030, 0
steels: SUFFICIENT. Hardest held-out protocols: tensile (11.8-15.1 %) and random (8.4-12.1 %).
V5-C identifiability: fits touching a box limit 5/18 (N=1), 11/18 (N=2), 14/18 (N=3); gamma collapse 0, 0, 1 of
18; fitted E 167-201 GPa, sigma_y0 212-402 MPa across steels and N.

## Deviations
(none)

## Note added 2026-11-09 (the text above is unchanged)
The six 'published' Voce-Chaboche values listed above could not be traced to a source we can open. The open companion paper
(de Castro e Sousa, Hartloper and Lignos, J. Eng. Mech. 2021, Table 1 for the steel IDs and protocols) gives Voce-Chaboche N=2 values
of 6.23, 5.75, 5.43, 4.88, 6.24 and 5.22 for IDs 7 to 12. Against these, our N=2 fits average 6.67 (mean +1.05 points, range +0.45 to
+1.77), so the implementation-check rule above ('within 1.0 point on at least 4 of 6 steels') is met on 3 of 6 steels only. The
statement that the Voce-Chaboche fits of the original paper were reproduced to within 0.5 points is withdrawn. The statement that
bars buckle in the Context paragraph is also unresolved in the literature (Kashani et al. 2019 describe their L/D 5 bars as not
buckled).

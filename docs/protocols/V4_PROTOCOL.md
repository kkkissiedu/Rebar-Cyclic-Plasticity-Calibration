# V4-T: series-spring test of the coupon model (frozen 2026-10-19, before any job was run)

## Question
In V3-F the full coupon used the lumped effective modulus (54-74 GPa, measured from the unloading
slopes of strain = actuator position / gauge) as the BAR's Young's modulus. The bars' own tension tests
(Kashani 2017, same B500C bars, extensometer) give E of about 194-206 GPa, and the lab's LCF strain is the
grip LVDT displacement (Kashani 2015), i.e. it includes load-train compliance. Physically the compliance
belongs in a spring in series with the bar, not in the bar. Does representing it that way (bar E = 200 GPa
plus a series spring) reduce the gap between the full-coupon response and the measured loops, and the
coupon's excess cyclic peak-stress decay?

This is a hypothesis test, not a calibration: the material parameters are NOT changed.

## Design
- Jobs: the frozen V3-F parameter sets (seed 1) for N = 1 and N = 3, 10 mm and 16 mm bars at 4 % and
  5 % strain (8 jobs). Material hardening cards, mesh, geometry, clamped bottom, lateral top constraints,
  loading history and the step definition (maximum increment period/200) are copied unchanged from the
  corresponding V3-F input files.
- Changes, and only these:
  1. Bar elastic modulus E = 200000 MPa (nu = 0.3), instead of the measured effective modulus.
  2. Machine spring: an axial spring (SPRING2, DOF 3) in series between the actuator node and the top face
     of the bar. The top-face nodes are tied in DOF 3 so the grip stays plane (all top nodes keep U3 equal to
     that of the first top node); lateral DOFs of the top face stay fixed as before.
  3. The prescribed displacement history (the same amplitude table) is applied to the ACTUATOR node,
     not to the top face.
  4. Spring stiffness from the measured system stiffness of that bar by series addition:
     1/K_m = (L0/A0) * (1/E_eff - 1/200000), with E_eff the bar's measured V3 E, L0 = D*5 mm,
     A0 = pi*D^2/4 (constant linear compliance, as in the V2/V3 treatment).
- Response compared: force on the actuator node over A0 (stress) against actuator displacement over L0
  (strain). This is what the experiment recorded (load cell and actuator position). The existing comparison
  script is used unchanged (it sums RF3 and averages U3 over the history set, which is now the actuator node).
- Scoring: unweighted per-cycle RMSE as in V3-F; headline cycles 4-6; also cycles 2-6; per-cycle tension
  and compression peaks; peak strains reached.

## Reference values (already known, from V3-F, same 8 parameter sets)
Coupon with E_eff in the bar (cycles 4-6): N=1: 83.9, 98.8, 105.4, 107.5; N=3: 123.5, 59.0, 77.0, 64.3
(mean 89.9 MPa). Surrogate (same parameters): mean 44.5 MPa. Coupon tension-peak change from cycle 2 to
cycle 6 averages -26 % in V3-F against -13 % in the experiment.

## Decision rule (fixed now)
The hypothesis "lumping the compliance into the bar's E is a main cause of the surrogate-coupon gap" is
SUPPORTED if (a) the mean cycles-4-6 coupon error over the 8 jobs falls to 67 MPa or below (a reduction of
at least 25 % from 89.9) AND (b) the mean absolute difference between the coupon's and the experiment's
tension-peak change from cycle 2 to cycle 6 is smaller than in V3-F. Otherwise it is NOT SUPPORTED. If it is
partly supported (only one of the two) it is reported as partial. Failed jobs are reported, not re-run.

## Not part of this test
Recalibration of the material parameters, the compliance being load dependent, imperfections, other meshes.
Those belong to a later V4 calibration protocol, written only if this test supports the approach.

## Result (2026-10-20): hypothesis NOT SUPPORTED
All 8 jobs solved. Mean cycles-4-6 coupon error: 105.6 MPa with E = 200 GPa + series spring, against 89.9 in
V3-F (E_eff in the bar) and 44.5 for the surrogate on the same parameters (N=1: 123.3 vs 98.9; N=3: 88.0 vs
81.0). Mean tension-peak change cycle 2 to 6: -29.1 % (V4-T), -25.8 % (V3-F), -13.8 % (experiment). Neither
criterion (<= 67 MPa; decay gap smaller than V3-F) was met. Conclusion: lumping the compliance into the bar's
E is not the main cause of the surrogate-coupon gap; the series spring is nevertheless the physically correct
representation (Moodley 2024 and Kashani 2015: LCF strain is the Instron integral LVDT; Kashani 2017 and
Moodley 2024: bar E 193-207 GPa from clip-on extensometer tension tests) and is kept as the baseline.
Side finding (not part of the frozen test): the frozen V3 fits overshoot the virgin cycle-1 ramp (at 5 % strain
the model reaches 643-722 MPa against 561-580 MPa measured), because V3 never scored cycle 1.

## Deviations
(none yet)

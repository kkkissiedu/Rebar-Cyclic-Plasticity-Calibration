# Ohno-Wang (1993) Cyclic-Plasticity Model - Reference for Python Surrogate + ABAQUS UMAT

Prepared: 2026-07-02

Primary source: Ohno, N. and Wang, J.D. (1993). "Kinematic hardening rules with critical
state of dynamic recovery, part I: formulation and basic features for ratchetting
behavior." *International Journal of Plasticity*, 9(3), 375–390.
DOI: [10.1016/0749-6419(93)90042-O](https://doi.org/10.1016/0749-6419(93)90042-O)

Companion paper: Ohno, N. and Wang, J.D. (1993). "Kinematic hardening rules with
critical state of dynamic recovery, part II: application to experiments of ratchetting
behavior." *International Journal of Plasticity*, 9(3), 391–403.
DOI: [10.1016/0749-6419(93)90043-P](https://doi.org/10.1016/0749-6419(93)90043-P)

> **Note on access.** The 1993 Part I/II papers are paywalled on ScienceDirect
> (`S0749641993900042O` / `074964199390043P`) and could not be retrieved as
> machine-readable text in this session (WebFetch returned HTTP 403 on ScienceDirect,
> Academia.edu, and ResearchGate mirrors; the only full PDF that downloaded successfully
> - via an academia.edu link - turned out to be a different, unrelated paper).
> Consequently the equations below are given with the **numbering convention used
> throughout the secondary/derivative literature** (Chen, Jiao & Kim 2005; Kobayashi &
> Ohno 2002; Bandyopadhyay et al. 2021; Chaboche's own review papers), which is
> internally consistent and matches the structure of the original paper's Eqs. (1)-(10),
> but the exact page-numbering of the original 1993 typesetting is **not independently
> confirmed** from primary-source text extraction. Where this matters, it is flagged
> explicitly. The physics/structure is corroborated from **three independent
> re-derivations**: (a) the LaTeX documentation and AceGen-generated Fortran of the
> KnutAM/MaterialModels open-source UMAT (Section 2 below), (b) Bandyopadhyay et al.,
> *Int. J. Plasticity* 136 (2021) 102887, which explicitly restates the Ohno-Wang
> evolution law (their Eq. 7) citing Ohno & Wang (1993a), and (c) the Abdel-Karim &
> Ohno (2000) generalization reproduced in a piecewise-linearization paper (CMES 2016),
> which recovers the Ohno-Wang rule as a special case of its more general evolution law.

---

## 1. Governing equations

### 1.1 Backstress decomposition

The total (kinematic hardening) backstress **α** is decomposed into *M* components,
each with its own evolution law:

```
α = Σ(i=1 to M) α_i                                          ... (1)
```

This is the same decomposition strategy Chaboche (1986) uses, but Ohno-Wang replaces
Chaboche's single dynamic-recovery term (linear in α_i) with a **critical-state**
dynamic-recovery term that only activates once the norm of α_i approaches a
component-specific saturation ("critical") value r_i.

### 1.2 Multiaxial evolution rule for each backstress component (Model II - the one universally implemented in FE codes)

Ohno & Wang presented two versions:

- **Model I (OW-I):** uses a Heaviside step function so the dynamic-recovery term
  switches on discontinuously at the critical surface ‖α_i‖ = r_i. This reproduces
  the "multilinear" (piecewise-linear) hardening rule exactly, but it predicts
  perfect uniaxial ratchetting shakedown (closed loops) and has a discontinuous
  tangent - it is rarely used in FE practice.
- **Model II (OW-II):** replaces the Heaviside function with a **power-law
  (Macaulay-bracket) nonlinearity of exponent m_i**, giving a smooth, differentiable
  transition into the critical-recovery regime. This is *the* "Ohno-Wang model" as
  normally cited and coded (Chaboche-compatible, differentiable, used in essentially
  every FE implementation, including the UMAT audited in Section 2).

The multiaxial rate form of Model II is:

```
dα_i = C_i * [ (2/3) r_i * dε_p  -  α_i * <dε_p : (α_i/‖α_i‖)> / r_i * (‖α_i‖/r_i)^{m_i} ]
                                                                ... (2)  [equiv. to OW Eq. (9)]
```

with the norm ‖α_i‖ = sqrt( (3/2) α_i:α_i ), dε_p the plastic strain-rate tensor, and
⟨·⟩ the Macaulay bracket, ⟨x⟩ = (x + |x|)/2 (i.e. ⟨x⟩ = x if x > 0, else 0).

An equivalent and very common way to write it (used e.g. by Chen, Jiao & Kim 2005 and
by Bandyopadhyay et al. 2021, their Eq. 7, generalized to M components) is:

```
dα_i = (2/3) C_i * dε_p  -  C_i * <n : α_i> / r_i * (‖α_i‖/r_i)^{m_i} * α_i * dp     ... (2')
```

where **n = (3/2) s'/σ_eq** is the flow-direction (unit normal to the yield surface,
see §1.4) and dp = sqrt((2/3) dε_p:dε_p) is the equivalent plastic strain increment.
Equations (2) and (2') are algebraically identical statements of the same rule; both
appear in the literature depending on whether the author normalizes by dp or works
directly with the tensor dε_p.

### 1.3 Exact uniaxial (1-D) form - the form to hard-code and unit-test against

For uniaxial loading, α_i, dε_p, and n collapse to scalars (α_i, dε_p, sign(dε_p)),
and the rule reduces to the textbook uniaxial Ohno-Wang equation:

```
dα_i = C_i * dε_p  -  C_i * <sign(dε_p) * α_i> * |α_i|^{m_i-1} / r_i^{m_i} * dε_p        ... (3)
```

or, written with the Macaulay bracket applied to the whole recovery multiplier
(the most commonly quoted uniaxial form):

```
dα_i = C_i [ dε_p  -  <(α_i/r_i) * sign(dε_p)> * (|α_i|/r_i)^{m_i-1} * dε_p ]            ... (3')
```

Interpretation: the direct-hardening term C_i·dε_p always acts; the dynamic-recovery
term only activates (via the Macaulay bracket) when α_i and the plastic-strain
increment have the same sign (loading is pushing α_i further from zero, i.e.
"unloading-then-reloading in the same sense"), and even then it is scaled down by
(|α_i|/r_i)^{m_i-1}, so it is negligible until |α_i| approaches the critical value
r_i, then grows sharply (a smooth approximation of a hard limit ‖α_i‖ ≤ r_i).

Sanity checks used by every re-implementation (including the KnutAM UMAT, confirmed
by reading its AceGen-derived residual, see §2.3):
- At α_i = 0: recovery term is 0, so dα_i = C_i·dε_p (pure linear/Prager hardening
  at the start of the first loading branch).
- As ‖α_i‖ → r_i with same-sign loading: the bracket saturates the recovery term,
  capping the hardening - the classic "sum of piecewise-linear segments" limit.
- m_i → ∞ collapses Model II onto Model I (the multilinear/piecewise-linear
  rule with hard corners at ‖α_i‖ = r_i); this is proven asymptotically in the
  1993 paper and repeated in every derivative work, e.g. Chen–Jiao–Kim (2005) and
  Kobayashi & Ohno (2002).
- m_i small (formally m_i → 1, or more precisely dropping the power-law
  sharpening) recovers the **Armstrong–Frederick (1966)** linear dynamic-recovery
  rule dα_i = C_i dε_p − γ_i α_i |dε_p| (equivalently the single-term
  **Chaboche (1986)** rule); i.e. Armstrong–Frederick/Chaboche is the m_i = 1
  special case of the Ohno-Wang recovery term with r_i = C_i/γ_i. This is the
  standard justification cited across the literature (e.g. Bari & Hassan 2002;
  Kang 2008 review) for why Ohno-Wang "contains" Armstrong-Frederick/Chaboche as
  a nested/limiting model - useful for surrogate-model bounds (start optimizer
  searches for m_i in [1, 20] and let the fit tell you whether the material wants
  AF-like or multilinear-like behavior).

### 1.4 Yield function and flow rule

Von Mises yield surface with combined isotropic/kinematic hardening:

```
Φ = f(σ - α) - (σ_y0 + R) ≤ 0                                  ... (4)
f(σ - α) = sqrt( (3/2) (s - α_dev):(s - α_dev) )                ... (5)
```

where s = dev(σ) is the deviatoric stress, σ_y0 the initial (virgin) yield stress,
and R the isotropic hardening stress (often modeled with a saturating exponential,
R = Q(1 - exp(-b·p)), independent of the choice of kinematic rule).

Associated (normality) flow rule:

```
dε_p = dλ * n,           n = ∂Φ/∂σ = (3/2) (s - α_dev) / f(σ - α)        ... (6)
dp   = dλ = sqrt( (2/3) dε_p:dε_p )                                       ... (7)
```

with the standard rate-independent Kuhn-Tucker (KKT) loading/unloading conditions:

```
Φ ≤ 0,   dλ ≥ 0,   dλ·Φ = 0                                               ... (8)
```

Consistency (Φ̇ = 0 during plastic flow) closes the system and yields the plastic
multiplier dλ in terms of the hardening moduli - standard return-mapping mechanics,
unaffected by which kinematic-hardening sub-rule (AF, Chaboche, Ohno-Wang) is used;
only Eq. (2)/(2') changes between models.

### 1.5 Parameters required per backstress component

| Symbol | Meaning | Units | Notes |
|---|---|---|---|
| σ_y0 | initial (virgin) yield stress | MPa | one value, not per-component |
| C_i | initial kinematic-hardening modulus of component i | MPa | same role as Chaboche's C_i |
| r_i (≡ ζ_i in some papers) | critical ("saturation") value of ‖α_i‖ for component i | MPa | sometimes written β_∞,i (e.g. KnutAM docs) |
| m_i | dynamic-recovery sharpness exponent for component i | dimensionless | m_i ≥ 1 typically; m_i → ∞ gives multilinear OW-I |
| M | number of backstress components | integer | typically 2–5; 4–5 needed to capture both LCF hysteresis shape and ratchet rate |
| Q, b (optional) | isotropic hardening saturation stress / rate | MPa, - | independent of kinematic model, needed for cyclic hardening/softening |

Total kinematic-hardening modulus at α_i = 0 equals Σ C_i, which must reproduce the
measured initial (or cyclic) hardening slope; Σ r_i (loosely) bounds the ultimate
kinematic contribution to flow stress, so **C_i and r_i are not independent of the
monotonic/cyclic stress-strain curve** - they are usually fit simultaneously with a
Voce-type isotropic term via nonlinear least squares/genetic algorithms against a
stabilized hysteresis loop plus at least one ratcheting test (uniaxial mean-stress
cycling), because a single symmetric loop alone under-constrains m_i.

**Practical guidance for a Python calibration surrogate:**
- Use M = 3–5 backstresses. With M ≤ 2 the model cannot simultaneously match the
  hysteresis-loop curvature and the ratchet rate; this is one of the central
  findings of Chen, Jiao & Kim (2005) and is reiterated in essentially every
  ratcheting-calibration paper since.
- Set bounds m_i ∈ [1, 30] (m_i = 1 ~ Armstrong-Frederick/Chaboche limit; m_i > 15
  is numerically close to the multilinear/Model-I limit and tends to be
  ill-conditioned in gradient-based fits - genetic/evolutionary optimizers are
  commonly used for this reason, e.g. the multi-objective GA calibration paper in
  *Int. J. Comp. Mat. Sci. Surf. Eng.*, 2014).
- Constrain Σ C_i·r_i to be consistent with the (σ_ult − σ_y0) monotonic hardening
  range as a warm-start / regularization, then let the optimizer refine.

---

## 2. Open-source ABAQUS UMAT implementation (verified)

### 2.1 Repository, file, and license

- **Repository:** `KnutAM/MaterialModels` - https://github.com/KnutAM/MaterialModels
  (author: Knut Andreas Meyer, Chalmers University of Technology; archived on Zenodo,
  DOI badge via `https://zenodo.org/badge/latestdoi/191778601`)
- **Cloned locally** (shallow clone, `git clone --depth 1`) to verify contents directly
  - confirmed present and readable.
- **License: MIT** (file `LICENSE` at repo root, copyright 2019 Knut Andreas Meyer).
  This is a **permissive license**, compatible with reuse/adaptation in this project,
  **with one caveat**: the LICENSE file itself states that portions of the code
  generated by the AceGen symbolic-code-generation tool are subject to AceGen's own
  redistribution restriction ("you may not include the codes generated by AceGen into
  other code if [the] other code is later use[d] for resale, rent or lease" -
  http://symech.fgg.uni-lj.si/Download.htm). The Ohno-Wang model's residual/Jacobian
  files (`ohnowang_acegen_mod.f90`, and the `AGfiles` under `GenFiniteStrain`) are
  AceGen-generated; the hand-written driver/dispatch files
  (`GeneralSmallStrain.f90`, `gss_module.f90`, `ohnowang.f90` - the parameter-count
  checker) are plain MIT. For a non-commercial academic/research surrogate + UMAT
  calibration tool this is not a practical obstacle, but it should not be
  redistributed as part of a for-resale/lease commercial product without checking
  the AceGen terms.
- **Specific files implementing Ohno-Wang:**
  - `models/GenSmallStrain/src/GeneralSmallStrain.f90` - the actual `SUBROUTINE UMAT(...)`
    entry point (standard ABAQUS UMAT argument list: `stress, statev, ddsdde, sse, spd,
    scd, rpl, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, dtemp, predef,
    dpred, cmname, ndi, nshr, ntens, nstatv, props, nprops, coords, drot, pnewdt,
    celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep, kinc`), which uses
    `use model_module` to select the active hardening rule.
  - `models/GenSmallStrain/src/ohnowang.f90` - the Ohno-Wang-specific `model_module`:
    defines `checkinput()`, which validates `nprops = 6 + 3*nback` and
    `nstatv = 2 + 6*nback`, and maps `E, ν → G, K`.
  - `models/GenSmallStrain/src/ohnowang_acegen_mod.f90` - AceGen-generated residual
    (`RF1`-`RF4`) and Jacobian (`dRdX1`-`dRdX4`) subroutines implementing the
    backward-Euler-integrated Ohno-Wang evolution law for 1–4 backstresses (license
    caveat above applies to this file specifically).
  - `models/GenSmallStrain/doc/ohnowang.md` and
    `models/GenSmallStrain/doc/latex/description.tex` - human-readable
    documentation stating the exact model equations (reproduced/verified in §2.2).
  - Rate-dependent (viscoplastic) variants also exist in the same folder:
    `ohnowang_rdep.f90`, `ohnowang_norton_acegen_mod.f90`,
    `ohnowang_cowsym_acegen_mod.f90`, `ohnowang_delobelle_acegen_mod.f90` (different
    overstress functions η(Φ) layered on top of the same Ohno-Wang kinematic rule).

### 2.2 Equation match against Ohno & Wang (1993) - verified from `description.tex`

The repo's own LaTeX documentation (`models/GenSmallStrain/doc/latex/description.tex`,
read in full) defines the kinematic-hardening-strain-conjugated evolution law as:

```
db_i/dλ = g(σ, β, β_i)
g_ohnowang(σ, β, β_i) = ν − ⟨ν:β_i⟩/β_{∞,i} · (f(β_i)/β_{∞,i})^{m_i} · ( δ · 3β_i^dev/(2f(β_i)) + (1-δ)·ν )
```

with ν = ∂Φ/∂σ (the flow-direction tensor, i.e. the "n" of §1.4 above), β_i the
backstress-conjugate variable, β_{∞,i} the saturation value (≡ r_i in §1.2-1.3
notation), and δ a blending parameter: **δ = 1 recovers the original Ohno-Wang
(1993) rule exactly** (the file states this explicitly: "The original Ohno Wang
model is obtained by setting δ=1"); δ = 0 gives the Burlet–Cailletaud (1986)
multiaxial modification instead, and intermediate δ gives the Delobelle blend. This
matches Eq. (2)/(2′) above term-for-term once δ is set to 1 and the Macaulay-bracket
power-law structure (⟨·⟩ and exponent m_i acting on the normalized backstress
magnitude) is compared side-by-side - confirmed by direct inspection, not taken on
faith from the repo's own claim.

### 2.3 PROPS (material-constant) ordering - confirmed from `ohnowang.md`

For a model with `nback` backstress components, `nprops = 6 + 3*nback`
(checked programmatically in `checkinput()` in `ohnowang.f90`, which aborts the
analysis with a diagnostic message if `nprops` or `nstatv` don't match this formula):

| PROPS index | Symbol | Meaning |
|---|---|---|
| 1 | E | Young's modulus |
| 2 | ν | Poisson's ratio |
| 3 | σ_y0 | Initial yield limit |
| 4 | H_iso | Isotropic hardening modulus |
| 5 | 1/Y_iso | Inverse of isotropic saturation stress |
| 6 | δ | Blend factor: δ=1 → pure Ohno-Wang (1993); δ=0 → Burlet-Cailletaud |
| 7 | H_k1 | Kinematic hardening modulus, backstress 1 (≡ C_1) |
| 8 | 1/Y_k1 | Inverse kinematic saturation stress, backstress 1 (≡ 1/r_1) |
| 9 | m_k1 | Exponent, backstress 1 (≡ m_1) |
| 10 | H_k2 | Kinematic hardening modulus, backstress 2 (≡ C_2) |
| 11 | 1/Y_k2 | Inverse kinematic saturation stress, backstress 2 (≡ 1/r_2) |
| 12 | m_k2 | Exponent, backstress 2 (≡ m_2) |
| ... | ... | pattern repeats in groups of 3 for each further backstress, up to `nback = 4` |

Note the repo parameterizes r_i as its **reciprocal** (`invYk_i = 1/r_i`) rather than
r_i directly - a common numerical-conditioning choice; this must be inverted when
mapping to/from the C_i, r_i, m_i notation used in §1.

### 2.4 State-variable (SDV/STATEV) layout - confirmed from `ohnowang.md`

`nstatv = 2 + 6*nback`:

| STATEV index | Meaning |
|---|---|
| 1 | κ (isotropic hardening stress) |
| 2 | λ (accumulated plastic multiplier, i.e. equivalent plastic strain) |
| 3–8 | backstress tensor component 1, α_1 (full 6-component tensor, ABAQUS ordering: 11,22,33,12,13,23) |
| 9–14 | backstress tensor component 2, α_2 |
| ... | pattern repeats in blocks of 6 for each further backstress |

`nback = (nstatv − 2)/6` is back-computed inside `checkinput()`, and only
`nback ∈ {1, 2, 3, 4}` is accepted (`nback_allowed` array in `ohnowang.f90`) - so if
a surrogate/calibration workflow wants M > 4 backstresses, this particular UMAT
would need to be extended (the AceGen generation script would need to be re-run for
`RF5`/`dRdX5`, etc.), or a different UMAT selected.

### 2.5 Other candidates checked (for completeness / due diligence)

- **KnutAM/umat** - an older, now-superseded standalone repo referenced in the
  `ref.bib` of `MaterialModels` (`url = https://github.com/KnutAM/umat`); superseded
  by `KnutAM/MaterialModels`, not used as the primary citation here.
- **dithoap/RESSForLab** - Voce-Chaboche UMATs only (no Ohno-Wang critical-state
  recovery term); not applicable to this task.
- **theysy/mml_subroutine_public**, **jpsferreira/UMAT-ABAQUS_library**,
  **sreepatiballa/UMAT_Lectures** - general UMAT teaching/collection repos found in
  search; spot-checked descriptions and none advertise an Ohno-Wang critical-state
  dynamic-recovery implementation specifically (mostly J2/Chaboche/Armstrong-Frederick
  teaching examples). Not cloned/verified in depth because KnutAM/MaterialModels
  already gave a directly-confirmed, permissively-licensed, line-by-line-readable
  match.
- No Zenodo/ResearchGate standalone "Ohno-Wang UMAT" code deposit (as opposed to a
  paper) was found independent of the KnutAM repository during this search.

**Conclusion: KnutAM/MaterialModels is the recommended reference implementation** -
permissively licensed (MIT, with the narrow AceGen-file caveat noted above), openly
cloneable, textually verified against its own equation documentation, and the only
candidate found with an explicit, auditable Ohno-Wang δ=1 special case plus a
documented PROPS/STATEV contract.

---

## 3. Typical calibrated parameter values (structural / rebar / carbon steel)

Numeric values for C_i, r_i (ζ_i), m_i are scattered across paywalled papers; the
exact tables in several of the papers below (Chen–Jiao–Kim 2005; Bari & Hassan 2002)
could **not** be extracted as machine-readable text in this session (ScienceDirect
and ResearchGate blocked programmatic fetch with HTTP 403). The ranges/orders of
magnitude below are triangulated from what **was** retrievable (abstracts, citing
papers, and the general consensus reported across the ratcheting-calibration
literature); treat these as **starting search bounds for an optimizer, not
ground truth for any single heat/grade of steel** - always recalibrate against your
own test data when available.

| Material | C_i range (MPa) | r_i / ζ_i range (MPa) | m_i range | M (# backstresses) | σ_y0 (MPa) | Source |
|---|---|---|---|---|---|---|
| AISI 316L / 316FR austenitic stainless steel | O(10³–10⁵) decreasing with i | O(10–200) | **m₁=1, m₂=2, m₃≈4-7, m₄→ large** (increasing with i; m→∞ for the last/"anchor" component to enforce a hard ratchet limit) | 4–5 | ≈120–200 | Kobayashi, M. & Ohno, N. (2002). "Implementation of cyclic plasticity models based on a general form of kinematic hardening." *Int. J. Numer. Methods Eng.*, 53(9), 2217–2238. DOI: 10.1002/nme.377 |
| Medium carbon steel S45C (JIS, ≈0.45%C) | multi-term, decreasing C_i with increasing i (typical O-W calibration pattern) | increasing r_i with i, last component largest (10²–10³ MPa range) | increasing m_i with i; last (highest-r) component given the largest m_i to sharply cap ratchet strain | 3–5 | ≈250–380 | Chen, X., Jiao, R., Kim, K.S. (2005). "On the Ohno–Wang kinematic hardening rules for multiaxial ratcheting modeling of medium carbon steel." *Int. J. Plasticity*, 21(1), 161–184. DOI: 10.1016/j.ijplas.2004.05.017 |
| Carbon/structural steel (general calibration methodology, mild steel σ_F≈250 N/mm²) | fit via 2–3 backstresses to a stabilized loop | fit alongside C_i | O(1–5) typically sufficient for structural-steel-grade LCF/seismic work (values much beyond ~10 rarely improve fit without ratchet data) | 2–4 (fatigue-only); 4–5 (if ratcheting must also be captured) | E≈2.0–2.1×10⁵ MPa, σ_y≈235–420 (mild/structural grades) | Halama, R., Sedlák, J., Šofer, M. "Choice and Calibration of Cyclic Plasticity Model with Regard to Subsequent Fatigue Analysis." *Engineering Mechanics*, 19(2), 87–97 (2012); and general practice summarized in Bari, S. & Hassan, T. (2002). "An advancement in cyclic plasticity modeling for multiaxial ratcheting simulation." *Int. J. Plasticity*, 18(7), 873–894. DOI: 10.1016/S0749-6419(01)00012-2 |
| Reinforcing steel / rebar (TMT, mild-steel-grade, seismic/LCF) | not tabulated in retrievable sources this session | not tabulated | not tabulated | typically 3–5, per general O-W LCF practice | ≈400–500 (TMT Fe500-class) | Application confirmed (Ohno-Wang used via UMAT in ABAQUS for TMT rebar low-cycle-fatigue and seismic performance studies), but the specific numeric C_i/r_i/m_i table was not accessible in this session. See search hits for "Low Cycle Fatigue Performance and Failure Analysis of Reinforcing Bar" and "Seismic Performance Assessment of a TMT Rebar" (ResearchGate/Academia.edu records located, full text not retrieved - **flag for manual follow-up** if rebar-specific numbers are required). |

### 3.1 Kobayashi & Ohno (2002) 316L m-values - what could be confirmed

The task specifically asks about the AISI 316L m-value sequence
(m₁=1, m₂=2, m₃=..., from Kobayashi & Ohno 2002). This pattern - **assigning
successively larger integer/near-integer exponents to successively higher-index
(larger-r) backstress components** - is a well-documented, widely-cited convention
in the Ohno-Wang literature (it lets the low-r, low-m components govern the smooth
hysteresis-loop shape while the high-r, high-m components act as a nearly-rigid
"backstop" that caps ratcheting), and Kobayashi & Ohno (2002), *Int. J. Numer.
Methods Eng.* 53(9):2217–2238 (DOI 10.1002/nme.377) is the standard citation for it.
However, this session's tools could not retrieve the paper's actual Table (Wiley
Online Library was not directly fetchable, and no open PDF mirror was found), so
**the exact m₁=1, m₂=2, ... sequence and the associated C_i/r_i values in MPa should
be treated as reported-in-secondary-literature / not independently verified
character-for-character in this session.** If exact reproduction of that table is
required for a paper or a defensible default-bounds file, obtain the PDF via
institutional access (Wiley) and re-extract; flagging this explicitly per the
instruction not to fabricate values.

### 3.2 Recommended default search bounds for the Python calibration surrogate

Given the above (and consistent with the general Ohno-Wang calibration methodology
described in Chen–Jiao–Kim 2005, Bari & Hassan 2002, and the multi-objective-GA
calibration paper in *Int. J. Comp. Mat. Sci. Surf. Eng.*, 2014), reasonable
**generic structural/carbon-steel starting bounds** for an M=4 backstress Ohno-Wang
fit (to be tightened once real test data is available) are:

- σ_y0: 250–500 MPa (structural/rebar-grade carbon steel)
- C_i: 1,000–200,000 MPa, one to two orders of magnitude apart per component,
  largest C on the lowest-index (fastest-saturating) backstress
- r_i (≡ ζ_i, ≡ β_∞,i): 10–400 MPa, increasing with component index i
- m_i: 1–20, increasing with component index i (m₁ ≈ 1, last component m_M as high
  as the optimizer/numerics will tolerate, often 10–50, to emulate the
  multilinear/Model-I ratchet cap)
- Isotropic hardening (if included): Q ∈ [-150, +150] MPa (softening negative,
  hardening positive), b ∈ [1, 50]

These are intentionally wide **search bounds**, not point estimates - the actual
optimizer (genetic algorithm / least-squares) should be left to explore within them
against the project's own stabilized-hysteresis-loop and/or ratchet-test data.

---

## 4. Reference list (full)

1. Ohno, N. & Wang, J.D. (1993). Kinematic hardening rules with critical state of
   dynamic recovery, part I: formulation and basic features for ratchetting
   behavior. *International Journal of Plasticity*, 9(3), 375–390.
   DOI: 10.1016/0749-6419(93)90042-O
2. Ohno, N. & Wang, J.D. (1993). Kinematic hardening rules with critical state of
   dynamic recovery, part II: application to experiments of ratchetting behavior.
   *International Journal of Plasticity*, 9(3), 391–403.
   DOI: 10.1016/0749-6419(93)90043-P
3. Kobayashi, M. & Ohno, N. (2002). Implementation of cyclic plasticity models
   based on a general form of kinematic hardening. *International Journal for
   Numerical Methods in Engineering*, 53(9), 2217–2238. DOI: 10.1002/nme.377
4. Chen, X., Jiao, R., & Kim, K.S. (2005). On the Ohno–Wang kinematic hardening
   rules for multiaxial ratcheting modeling of medium carbon steel. *International
   Journal of Plasticity*, 21(1), 161–184. DOI: 10.1016/j.ijplas.2004.05.017
5. Bari, S. & Hassan, T. (2002). An advancement in cyclic plasticity modeling for
   multiaxial ratcheting simulation. *International Journal of Plasticity*, 18(7),
   873–894. DOI: 10.1016/S0749-6419(01)00012-2
6. Halama, R., Sedlák, J., & Šofer, M. (2012). Choice and calibration of cyclic
   plasticity model with regard to subsequent fatigue analysis. *Engineering
   Mechanics*, 19(2), 87–97.
7. Bandyopadhyay, R., Gustafson, S.E., Kapoor, K., Naragani, D., Pagan, D.C., &
   Sangid, M.D. (2021). Comparative assessment of backstress models using
   high-energy X-ray diffraction microscopy experiments and crystal plasticity
   finite element simulations. *International Journal of Plasticity*, 136, 102887.
   DOI: 10.1016/j.ijplas.2020.102887 (reproduces the Ohno-Wang evolution law as
   their Eq. 7, citing Ohno & Wang 1993a/1993b; verified by direct text extraction
   in this session).
8. Meyer, K.A. (2019). `KnutAM/MaterialModels` - User material models for Abaqus
   (UMAT). GitHub repository, MIT License, Zenodo-archived (DOI badge:
   `zenodo.org/badge/latestdoi/191778601`). https://github.com/KnutAM/MaterialModels
   - cloned and verified directly in this session (LICENSE file, `ohnowang.f90`,
   `ohnowang_acegen_mod.f90`, `ohnowang.md`, `description.tex`).
9. Armstrong, P.J. & Frederick, C.O. (1966). A mathematical representation of the
   multiaxial Bauschinger effect. CEGB Report RD/B/N731, Berkeley Nuclear
   Laboratories (the m_i=1 / linear-recovery limiting case referenced in §1.3).
10. Chaboche, J.L. (1986/1989/1991). Various papers establishing the multi-component
    nonlinear-kinematic-hardening decomposition that both Chaboche's own model and
    Ohno-Wang build on (Σ α_i decomposition, §1.1).

---

## 5. Open items / follow-up if higher precision is needed

- The primary 1993 Part I/II papers were not obtainable as extractable text in this
  session (institutional/paywall access needed) - the equation numbering above
  follows the near-universal secondary-literature convention, not confirmed
  character-for-character against the original typeset pages.
- The exact Kobayashi & Ohno (2002) 316L parameter table (C_i, r_i, and the
  m_i = 1, 2, ... sequence) needs Wiley institutional access to confirm numerically;
  flagged in §3.1.
- No numeric Ohno-Wang table specific to reinforcing bar/rebar steel was retrievable
  this session; flagged in the rebar row of the §3 table - recommend a targeted
  follow-up search once institutional journal access is available (candidates:
  the "Low Cycle Fatigue Performance and Failure Analysis of Reinforcing Bar" and
  "Seismic Performance Assessment of a TMT Rebar" papers found on
  ResearchGate/Academia.edu but not full-text-accessible here).

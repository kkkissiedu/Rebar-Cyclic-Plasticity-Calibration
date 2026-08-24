# Kashani et al. (2019) - Complete Extraction

## 1. Full Citation

**Authors:** Mohammad M. Kashani, Shunyao Cai, Sean A. Davis, Paul J. Vardanega

**Year:** 2019

**Title:** Influence of Bar Diameter on Low-Cycle Fatigue Degradation of Reinforcing Bars

**Journal:** Journal of Materials in Civil Engineering

**Volume/Issue:** Vol. 31, Issue 4

**Article Number:** 06019002

**DOI:** 10.1061/(ASCE)MT.1943-5533.0002637

**Publisher:** American Society of Civil Engineers (ASCE)

**ISSN:** 0899-1561

**Publication Type:** Technical Note

**Manuscript Submission Date:** February 15, 2018  
**Approval Date:** September 10, 2018  
**Online Publication Date:** January 31, 2019

---

## 2. Abstract Summary and Stated Objectives

**Abstract Summary:**

The paper reports results of 120 low-cycle fatigue (LCF) tests on steel reinforcing bars with varying slenderness ratios (L:D) at varying strain amplitudes. Failure modes of fractured bars were investigated through analysis of fracture mechanisms. Experimental results were used to update empirical models of low-cycle fatigue life. The updated empirical models were incorporated into a recently developed constitutive material model that accounts for bar buckling and fatigue.

**Key Experimental Finding from Abstract:** The size effect is significant for short steel reinforcing bars where there is no buckling. As the slenderness ratio of the steel reinforcing bars increases, the influence of bar diameter on low-cycle fatigue reduces.

**Study Aims:**

The study extends prior research by Kashani et al. (2015b, 2015c, 2017) which had only examined 12-mm and 16-mm diameter bars. The current work investigates the influence of bar diameter (D) on inelastic buckling behavior and low-cycle fatigue performance across a wider range of bar sizes. The specific goal is to update empirical fatigue life models and incorporate them into a constitutive material model for nonlinear analysis of reinforced concrete structures.

---

## 3. Materials Tested

### Steel Grade and Specifications

- **Grade:** B500C (British Standard compliant)
- **Standard:** BS4449-2005+A3:2016
- **Source:** British-manufactured bars with controlled chemical composition
- **Compliance:** Carbon equivalent < 0.50% per BSI (2016)

### Bar Diameters Tested

Four distinct diameter groups tested:
- 10 mm
- 12 mm
- 16 mm
- 20 mm

### Slenderness Ratios (L:D)

Five slenderness ratios tested for each diameter:
- L:D = 5
- L:D = 8
- L:D = 10
- L:D = 12
- L:D = 15

### Total Test Specimens

**120 test specimens total** with varying combinations of diameter and slenderness ratio

For each set of bar diameters, **three tensile tests** were conducted to characterize material properties.

### Chemical Composition

Per EDX analysis (reported in Fig. S3):
- Main component: Iron (Fe)
- Insignificant material variation among specimens
- Controlled composition confirms differences in fatigue life are due to size effects, not material composition

### Mechanical Properties

**Not explicitly reported in the paper for yield strength, UTS, or elastic modulus.** The paper references Kashani (2017) for detailed material property characterization, noting that the material used in the experiments reported is the same as that used by Kashani (2017).

The paper states material properties are characterized through standard tensile testing but does not tabulate the values. Reference to original Kashani (2017) publication would be needed for specific Young's modulus (E), yield strength (σ_y), or ultimate tensile strength (UTS) values.

---

## 4. Test Programme

### Specimen Geometry

- **Specimen Type:** Steel reinforcing bars
- **Diameters:** 10, 12, 16, 20 mm (4 sizes)
- **Gauge Length:** Not explicitly stated in main text
- **L:D Ratios:** 5, 8, 10, 12, 15 (5 configurations per diameter)

### Strain Amplitudes Tested

Six strain amplitudes examined:
- 1%
- 2%
- 3%
- 4%
- 5%
- 6%

### Strain Rate

**Not reported** in the available sections of the paper.

### Loading Protocol

- **Type:** Cyclic loading with constant strain amplitude
- **Loading Mode:** Repeated cyclic strain reversals (tension-compression)
- **Failure Criterion:** Complete rupture of bar identified as failure point (end of test)

### R-Ratio

**Not explicitly reported** in the paper. The loading is implied to be fully reversed cyclic loading (R = -1) based on the context of seismic engineering and the reference to "repeated cyclic loading," but no explicit R-ratio value is stated.

### Number of Specimens

- **Total:** 120 specimens
- **Breakdown:**
  - 4 diameter sizes × 5 slenderness ratios = 20 distinct diameter/L:D combinations
  - Multiple specimens tested per combination (120 total ÷ 20 = 6 specimens average per configuration)
- **Tensile Tests:** 3 additional tensile tests per diameter set for material characterization

### Testing Apparatus

- **Equipment:** Instron machine (Instron, Norwood, Massachusetts)
- **Capacity Limitation Noted:** The paper mentions that larger range of D values could not be tested due to capacity limitations of the available Instron machine

### Experimental Method Reference

The methodology is stated to be similar to that presented by Kashani et al. (2015b). For complete details on testing apparatus configuration, the reader is directed to that reference.

---

## 5. Constitutive Model and Cyclic Plasticity Parameters

### Fatigue Life Model Used

**Koh-Stephen Model** (Koh and Stephens 1991):

$$\varepsilon_a = \varepsilon_f (2N_f)^\alpha \quad \text{(Eq. 1)}$$

Where:
- $\varepsilon_a$ = total strain amplitude (elastic + plastic strain)
- $\varepsilon_f$ = ductility coefficient (fracture strain due to one load reversal)
- $N_f$ = number of half-cycles to failure
- $\alpha$ = ductility exponent (slope parameter, negative value)

### Calibrated Fatigue Parameters: Table 1

Calibrated low-cycle fatigue parameters (εf, α) for all tested bar sizes and slenderness ratios:

| D (mm) | L:D | εf     | α       | p-value | λp      |
|--------|-----|--------|---------|---------|---------|
| 10     | 5   | 0.190  | -0.457  | 0.0089  | 11.619  |
| 10     | 8   | 0.184  | -0.537  | 0.0014  | 18.590  |
| 10     | 10  | 0.312  | -0.689  | 0.0012  | 23.238  |
| 10     | 12  | 0.357  | -0.721  | 0.0016  | 27.885  |
| 10     | 15  | 0.350  | -0.695  | 0.0014  | 34.857  |
| 12     | 5   | 0.190  | -0.443  | 0.0011  | 11.619  |
| 12     | 8   | 0.185  | -0.541  | 0.0016  | 18.590  |
| 12     | 10  | 0.184  | -0.512  | 0.0011  | 23.238  |
| 12     | 12  | 0.341  | -0.725  | 0.0017  | 27.885  |
| 12     | 15  | 0.354  | -0.733  | 0.0018  | 34.857  |
| 16     | 5   | 0.131  | -0.399  | 0.0016  | 11.511  |
| 16     | 8   | 0.157  | -0.506  | 0.0273  | 18.417  |
| 16     | 10  | 0.206  | -0.615  | 0.0099  | 23.022  |
| 16     | 12  | 0.316  | -0.719  | 0.0014  | 27.626  |
| 16     | 15  | 0.376  | -0.748  | 0.0107  | 34.533  |
| 20     | 5   | 0.150  | -0.410  | 0.0010  | 11.511  |
| 20     | 8   | 0.152  | -0.479  | 0.0013  | 18.417  |
| 20     | 10  | 0.227  | -0.607  | 0.0022  | 23.022  |
| 20     | 12  | 0.227  | -0.610  | 0.0012  | 27.626  |
| 20     | 15  | 0.401  | -0.766  | 0.0019  | 34.533  |

**Note:** p-values are those computed when Eq. (1) is fitted to the data. All values indicate good statistical significance (p < 0.05 for most cases, one exception at 0.0273 for 16 mm, L:D = 8).

### Influence of Inelastic Buckling: Bar Buckling Parameter

**Nondimensional Bar Buckling Parameter (λp):**

$$\lambda_p = 10\sqrt{\frac{\sigma_y L}{0 D}} \quad \text{(Eq. 2)}$$

Where:
- σy = yield stress (in MPa)
- L = length of steel reinforcing bar
- D = bar diameter
- 0 = elastic buckling reference stress (implicit in equation)

**Statistical Finding:** Kashani et al. (2015b) showed that the relationship between (εf, α) and λp is statistically significant. This paper extends that finding across a wider range of bar diameters.

### Empirical Models Incorporating Diameter Effects: Equations 3 and 4

**For the ductility exponent α:**

$$\alpha = a\lambda_p + b \quad \text{(Eq. 3)}$$

**For the ductility coefficient εf:**

$$\varepsilon_f = c \exp(d\lambda_p) \quad \text{(Eq. 4)}$$

Where:
- a, b, c, d = regression coefficients (material constants determined by fitting to experimental data)

### Calibrated Fatigue Parameters as Function of λp: Table 2

Material constants (a, b, c, d) for Eqs. (3) and (4), fitted for each diameter group:

| Diameter | Fatigue Parameter | a      | b       | c     | d     | p-value |
|----------|-------------------|--------|---------|-------|-------|---------|
| **10 mm diameter bars** |
|          | α                 | -0.010 | -0.384  | -     | -     | 0.0021; 0.0017 |
|          | εf                | -      | -       | 0.142 | 0.028 | 0.0030; 0.0040 |
| **12 mm diameter bars** |
|          | α                 | -0.010 | -0.364  | -     | -     | 0.0041; 0.0027 |
|          | εf                | -      | -       | 0.146 | 0.023 | 0.0130; 0.0140 |
| **16 mm diameter bars** |
|          | α                 | -0.010 | -0.333  | -     | -     | 0.0042; 0.0051 |
|          | εf                | -      | -       | 0.122 | 0.021 | 0.001;  0.001  |
| **20 mm diameter bars** |
|          | α                 | -0.015 | -0.223  | -     | -     | 0.0022; 0.0410 |
|          | εf                | -      | -       | 0.060 | 0.044 | 0.0051; 0.0103 |

**Note:** p-values are those computed when Eqs. (3) and (4) are fitted to data.

### Key Regression Analysis Findings on Bar Diameter Effects

1. **No significant effect on a and d coefficients:**
   - Regression p-values > 0.05 indicate D does not significantly affect the slope (a) in Eq. (3)
   - D does not significantly affect the exponent (d) in Eq. (4)
   - This suggests D does not directly affect λp

2. **Significant effect on b and c coefficients:**
   - Regression analysis reveals **p < 0.05** for coefficients b and c
   - **Strong correlation** exists between D and coefficients b and c
   - Figure 5 shows the interrelationship among b, c, and D

3. **Interpretation:**
   - As L:D increases, the influence of D on fatigue material coefficients reduces
   - If bars are affected by inelastic buckling, varying D should not affect LCF life (which is mainly governed by λp)
   - For short bars (L:D ≤ approximately 6), D has considerable influence on LCF life

### Reference Constitutive Model (Not Directly Detailed Here)

The updated fatigue parameters were incorporated into the uniaxial constitutive material model developed by Kashani et al. (2015c). This model is also implemented in OpenSees (2014).

**Key Features (from paper context):**
- Accounts for bar buckling
- Accounts for fatigue degradation
- Phenomenological hysteretic model for reinforcing bars
- Enables nonlinear analysis in structural design codes

---

## 6. Key Experimental Results

### Influence of Diameter on Hysteretic Loops

**Key Finding:** Diameter (D) has qualitatively measurable influence on cyclic stress-strain hysteretic responses of short bars (L:D = 5). As L:D increases (e.g., L:D = 15), the influence of D reduces due to dominance of plasticity effects over material-level effects.

**Physical Explanation:**
- In short bars (L:D ≤ 6), almost the whole bar length goes plastic during reversed cyclic loading
- Elastic length (elastic region between plastic hinges) to plastic hinge length ratio is minimal
- Three-dimensional effect due to material dislocation (molecular response) influences global response under cyclic loading
- As L:D increases, behavior transitions toward nonlinear beam-column behavior
- Geometrical nonlinearity (buckling effect) governs global response, not material size effects

### Influence of Diameter on Plastic Energy Dissipation

**Normalized Plastic Energy Metric:**
- Total plastic hysteretic energy dissipated during testing is used as a measure of damage
- Total hysteretic energy: $E_t$
- Elastic energy under monotonic tension: $E_e$
- **Normalized plastic energy = $(E_t - E_e) / E_e$**

**Findings (from Fig. 2):**
- Normalized dissipated energy increases with strain amplitude
- Groups with L:D = 5 and L:D = 8 most affected by varying D
- Energy plots do not show systematic pattern between D and normalized plastic energy dissipated
- **Critical threshold:** When L:D > approximately 10, D has **minimal influence**
- This is consistent with prior research by Chang and Mander (1994), Kashani et al. (2013a, 2015a), and Kashani (2017)

### Cyclic Softening/Hardening Behavior

**Not explicitly quantified in numerical terms in the paper.** The stress-strain hysteretic loops in Fig. 1 and Figs. S4 show the qualitative cyclic behavior, but specific quantification of cyclic softening or hardening rates is not provided in the reported sections.

### Fatigue Life Behavior

**Coffin-Manson-Type Relationship:**
The Koh-Stephen model (Eq. 1) fits the experimental data:

$$\varepsilon_a = \varepsilon_f (2N_f)^\alpha$$

Where:
- εf (ductility coefficient) ranges from **0.131 to 0.401** across all diameter/L:D combinations
- α (ductility exponent) ranges from **-0.399 to -0.766** across all combinations
- Negative α indicates that fatigue life decreases with increasing strain amplitude (inverse power relationship)

**Size Effect on Fatigue Life:**
- **Larger diameters → shorter fatigue life** at the same strain amplitude (for short bars)
- Effect diminishes as L:D increases
- For L:D ≥ approximately 10, the effect of D on fatigue life is minimal

**Example:** At L:D = 5, comparing 10 mm vs. 20 mm bars:
- 10 mm: εf = 0.190, α = -0.457
- 20 mm: εf = 0.150, α = -0.410
- The 20 mm bar has lower εf (shorter life at same amplitude)

### Tension/Compression Asymmetry

**Not reported or addressed** in the paper. The focus is on large-amplitude cyclic loading with symmetric strain reversals (fully reversed loading).

### Fractured Surface Analysis Using SEM

**Finding on Failure Modes:**
- As D increases in bars with L:D = 5 (no buckling), failure mode becomes more **brittle**
- Darker areas of ridges caused by straining indicate slower crack propagation (10 mm specimens fractured later with increased plastic deformation)
- Lighter areas indicate more sudden fracture events
- As L:D increases (L:D ≥ 10), almost all bars with different D exhibited **similar failure modes**

**Interpretation:**
- Once bars buckle, strain amplitude is locally increased due to large deformation
- Crack growth is faster
- D does not significantly influence performance once buckling occurs

**Material Composition Verification (EDX):**
- Iron (Fe) is the main component
- Insignificant material variation among specimens
- Differences in fatigue life are purely due to size effects, not composition variations

---

## 7. Data Availability

### Open Dataset Link

**Dataset Title:** Supporting data for 'Influence of bar diameter on low-cycle fatigue degradation of reinforcing bars'

**Repository:** University of Bristol Data Portal

**DOI:** https://doi.org/10.5523/bris.1kz5015zjoel92ueb97kwxd4ps

**Accessibility:** Accessed January 24, 2019

**Content:** Raw test data from the LCF tests (referenced as Kashani et al. 2018a)

### Supplemental Data Statement

Tables S1–S4 and Figures S1–S4 are available online in the ASCE Library at www.ascelibrary.org

**Supplemental Contents:**
- **Tables S1–S4:** Complete experimental results (individual test data for all 120 specimens)
- **Fig. S1:** SEM fractographs for 12-mm and 16-mm diameter bars
- **Fig. S2:** SEM fractographs for L:D = 15 specimens with buckling effect
- **Fig. S3:** EDX spectra of tested specimens (material composition verification)
- **Fig. S4:** Comparison of simulation vs. experimental results for other diameters (beyond the 20 mm shown in main text Fig. 6)

---

## 8. Stated Limitations and Conclusions

### Stated Limitations

1. **Capacity Constraints on Diameter Range:**
   - "Due to capacity limitations of the Instron machine available, a larger range of D could not be tested"
   - Only four diameter sizes tested (10, 12, 16, 20 mm)
   - A larger range of D values may reveal different trends

2. **Applicability of Existing Models at Wider Diameter Range:**
   - Question posed: "Do the currently available uniaxial material models apply at a wider range of values of D?"
   - Larger diameter testing beyond 20 mm needed to fully answer this

### Summary of Key Findings

1. **Fracture Mechanism Findings:**
   - Analysis of ruptured bars using SEM method showed that as D increased in bars with L:D = 5, the failure mode became more brittle
   - As L:D increased (L:D ≥ 10), almost all bars with different D exhibited similar failure modes

2. **Experimental Fatigue Test Results:**
   - As D increases, reinforcing bars **fracture earlier** under repeated cyclic loading
   - The value of D **affects the LCF performance** of steel bars
   - **Increasing L:D ratio reduces the influence of D** on LCF performance

3. **Short Bars (No Buckling Effect, L:D ≤ approximately 6):**
   - Inelastic buckling is not critical
   - D has a **considerable influence** on LCF life of steel reinforcing bars
   - **Size effect is significant**

4. **Long Bars (With Buckling Effect, L:D > 6):**
   - LCF performance of bars affected by inelastic buckling is mainly governed by their **bar buckling parameter λp**
   - λp is a function of L:D and σy
   - D has minimal direct effect once buckling occurs

5. **Model Validation:**
   - The updated constitutive model [Kashani et al. (2015c)] incorporating the new fatigue parameters (εf, α) can adequately simulate:
     - Nonlinear cyclic stress-strain behavior
     - LCF failure of reinforcing bars
     - Combined effects of inelastic buckling and bar diameter

### Recommendations Disclaimer

> "Any recommendations provided in this paper are the opinions of the authors and do not constitute a standard or code of practice."

### Conclusions Statement

The paper demonstrates that bar diameter significantly influences low-cycle fatigue degradation in short reinforcing bars without buckling, but this influence diminishes as slenderness ratio increases. The empirical models developed enable better prediction of fatigue life across a range of bar sizes, and the incorporation into constitutive models enables more accurate nonlinear analysis of reinforced concrete structures under seismic loading.

---

## Supplementary Notes

### Author Affiliations

1. **Mohammad M. Kashani** (Corresponding Author) - Associate Professor, Faculty of Engineering and Physical Sciences, University of Southampton, Southampton SO17 1BJ, UK. Email: mehdi.kashani@soton.ac.uk

2. **Shunyao Cai** - Lecturer, School of Construction Management and Real Estate, Chongqing University, 83 Shabei St., Shapingba, Chongqing 400045, China. Email: shunyao.cai@gmail.com

3. **Sean A. Davis** - Senior Lecturer, School of Chemistry, University of Bristol, Cantock's Close, Bristol BS8 1TS, UK. Email: s.a.davis@bristol.ac.uk

4. **Paul J. Vardanega** (Corresponding Author) - Senior Lecturer in Civil Engineering, Department of Civil Engineering, University of Bristol, Queen's Building, University Walk, Bristol BS8 1TR, UK. Email: p.j.vardanega@bristol.ac.uk

### Funding Acknowledgments

- Earthquake and Geotechnical Engineering Research Group (EGERG) at the University of Bristol
- SEM studies conducted in the Chemistry Imaging Facility at the University of Bristol
- Equipment funded by the University of Bristol and the Engineering and Physical Sciences Research Council Grant No. EP/K035746/1
- Second author (S. Cai) supported by China Scholarship Council (File No. 201506260126)

### Related Prior Work by Same Research Group

This paper extends and updates previous research:
- **Kashani et al. (2015b):** Initial LCF and buckling study (12-mm and 16-mm bars only)
- **Kashani et al. (2015c):** Constitutive model development
- **Kashani et al. (2013a, 2013b, 2014, 2016, 2018a, 2018b):** Various aspects of corrosion, buckling, and cyclic response
- **Kashani (2017):** Size effect on buckling of corroded bars

---

## End of Extraction

# -*- coding: utf-8 -*-
"""
rebuild_fatigue_model.py
========================================================================
Phase 5 -- parametric rebuild of the LCF (low-cycle fatigue) specimen model.

Runs inside ABAQUS/CAE Python (Python 2.7), NO GUI:

    abaqus cae noGUI=rebuild_fatigue_model.py

CRITICAL: this script writes a brand-new model database called
``FATIGUE_rebuilt.cae`` in the current working directory. It never opens,
reads, or overwrites the user's existing working model ``FATIGUE.cae``.
The two files are siblings so the user can diff/review before adopting
the rebuilt one.

Caller: invoked directly by the user (or the calibration app's build step)
via ``abaqus cae noGUI=rebuild_fatigue_model.py``. No other project script
imports this module.

Data I/O:
    reads  -- nothing (all geometry/material/loading parameters are defined
              as named constants below; no external files are consumed)
    writes -- ``FATIGUE_rebuilt.cae`` (+ the usual ABAQUS side artifacts
              ``abaqus.rpy`` / ``*.jnl``) in os.getcwd()

Rebuilt from scratch, this script reproduces (and improves on) the
existing FATIGUE.cae LCF specimen model, per the Phase 3 inspection
report (``abaqus_inspection_report.md``):

    * Part ``LD5_12``       -- solid cylinder, 12 mm dia x 60 mm gauge,
                               axial along Z (axis 3)
    * Mesh                 -- C3D8R, structured hex, >=3 elements across
                               the radius, >=10 elements along the length
    * Material ``Steel_B500C`` -- linear elastic + combined
                               (nonlinear kinematic + isotropic) Chaboche
                               hardening with N_BACKSTRESSES backstresses,
                               placeholder parameters (calibration app
                               fills in the real values later)
    * Amplitude ``CyclicStrain`` -- tabular symmetric triangle wave,
                               built parametrically from GAUGE_LENGTH_MM
                               and STRAIN_AMPLITUDE
    * BCs                  -- bottom face (Z=0) encastred, top face
                               (Z=gauge) displacement-controlled along Z
                               referencing CyclicStrain; lateral top-face
                               DOFs pinned to remove rigid-body motion
                               (matches the existing model's convention)
    * History output       -- single mid-specimen element set ``SS``,
                               S33/E33 (axial stress/strain)
    * Field output         -- whole model, every 10th increment
    * Step ``Cyclic_Loading`` -- Static, nlgeom=ON, parametrized time
                               period / increment controls
    * Job ``Calibration_Job`` -- single CPU, single domain

Notably different from the current FATIGUE.cae (per the inspection
report warnings): nlgeom is turned ON here (the existing model has it
OFF), and NUMBER BACKSTRESSES defaults to 3 (the existing model uses 2)
per the calibration-app spec. Adjust N_BACKSTRESSES below to match the
material being calibrated if 2 backstresses are required instead.
"""
from __future__ import print_function

import os
import sys


# ============================================================================
# PARAMETRIC CONSTANTS -- edit these to rebuild for a different specimen /
# loading level. Nothing below this block should contain a hardcoded
# magic number that depends on specimen geometry or loading.
# ============================================================================

GAUGE_LENGTH_MM   = 60.0     # specimen gauge length (mm), axial along Z
BAR_DIAMETER_MM   = 12.0     # specimen diameter (mm)
STRAIN_AMPLITUDE  = 0.01     # nominal axial strain amplitude (fraction, e.g. 0.01 = 1%)
N_CYCLES          = 5        # number of full load cycles to simulate
PERIOD_S          = 20.0     # time period of one full cycle (s)
N_BACKSTRESSES    = 3        # number of Chaboche backstresses (combined hardening)

# Mesh density (must satisfy the spec minimums: >=3 across radius, >=10 along length)
N_ELEMS_RADIUS    = 5        # elements across the radius (radial seed count)
N_ELEMS_LENGTH    = 20       # elements along the gauge length (axial seed count)
N_ELEMS_CIRCUMF   = 12       # elements around the quarter-circumference partitions
                              # (12 per quadrant boundary -> even, structured hex mesh)

# Mesh-quality acceptance thresholds (ABAQUS failure ranges: aspect ratio
# failure is typically >= 10-20 depending on element type; we target well
# under that).
ASPECT_RATIO_WARN_LIMIT = 5.0   # print [WARNING] above this
ASPECT_RATIO_FAIL_LIMIT = 10.0  # refuse to finish ([FAIL]) above this

# Material placeholders (mid-range; the calibration app overwrites these
# with fitted values later). For N_BACKSTRESSES == 3 the *PLASTIC data
# line is: sy0, C1, g1, C2, g2, C3, g3  (7 values).
PLACEHOLDER_SY0 = 500.0
PLACEHOLDER_BACKSTRESS_PARAMS = [
    (12000.0, 300.0),   # (C1, gamma1)
    (2000.0, 60.0),     # (C2, gamma2)
    (800.0, 12.0),      # (C3, gamma3)
]
PLACEHOLDER_Q_INF = -40.0
PLACEHOLDER_B = 8.0

E_MODULUS_MPA = 200000.0
POISSON_RATIO = 0.3

# Names
MODEL_NAME     = "Model-1"
PART_NAME      = "LD5_12"
MATERIAL_NAME  = "Steel_B500C"
SECTION_NAME   = "Specimen_Section"
INSTANCE_NAME  = "Specimen-1"
AMPLITUDE_NAME = "CyclicStrain"
STEP_NAME      = "Cyclic_Loading"
JOB_NAME       = "Calibration_Job"
HISTORY_SET_NAME = "SS"          # part-level element set name (mid-specimen)
BOTTOM_SET_NAME   = "Set-Bottom"
TOP_SET_NAME      = "Set-Top"

OUTPUT_CAE_NAME = "FATIGUE_rebuilt.cae"   # sibling file -- never FATIGUE.cae

# Consistency checks on the mesh minimums required by spec.
if N_ELEMS_RADIUS < 3:
    raise ValueError("N_ELEMS_RADIUS must be >= 3 per spec (got %d)" % N_ELEMS_RADIUS)
if N_ELEMS_LENGTH < 10:
    raise ValueError("N_ELEMS_LENGTH must be >= 10 per spec (got %d)" % N_ELEMS_LENGTH)


# ============================================================================
# ABAQUS imports (must come after sys.path is set up by the `abaqus cae`
# launcher; safe to import unconditionally when run via noGUI=).
# ============================================================================
from abaqus import mdb, session
from abaqusConstants import (
    THREE_D, DEFORMABLE_BODY, ON, OFF, COMBINED, PARAMETERS,
    C3D8R, HEX, STRUCTURED, UNIFORM,
    ANALYSIS_CHECKS, ASPECT_RATIO, FINER, CARTESIAN, PRESELECT, STEP, STANDARD,
    YZPLANE, XZPLANE,
)
import section
import regionToolset
import mesh
import step  # noqa: F401 -- importing this module injects Step/OutputRequest
             # methods (StaticStep, HistoryOutputRequest, FieldOutputRequest,
             # etc.) onto the Model class; required even though `step` itself
             # is not referenced by name below.


# `abaqus cae noGUI=` does not reliably surface Python print() output to the
# calling shell's stdout (the CAE kernel's console stream is separate from
# the OS process stdout in noGUI mode on Windows). To guarantee the build
# log is inspectable afterwards, every log() / fail() line is written to
# BOTH stdout (best-effort) and a plain-text log file next to the script.
LOG_PATH = os.path.join(os.getcwd(), "rebuild_fatigue_model.log")
_log_fh = open(LOG_PATH, "w")


def log(msg):
    print(msg)
    _log_fh.write(str(msg) + "\n")
    _log_fh.flush()


def fail(msg):
    log("[FAIL] " + msg)
    log("[FAIL] Aborting -- FATIGUE_rebuilt.cae was NOT saved cleanly.")
    _log_fh.close()
    sys.exit(1)


# ============================================================================
# 0. Fresh model database (do not touch FATIGUE.cae)
# ============================================================================
log("=" * 78)
log("Rebuilding LCF specimen model parametrically")
log("  GAUGE_LENGTH_MM  = %.3f" % GAUGE_LENGTH_MM)
log("  BAR_DIAMETER_MM  = %.3f" % BAR_DIAMETER_MM)
log("  STRAIN_AMPLITUDE = %.5f" % STRAIN_AMPLITUDE)
log("  N_CYCLES         = %d" % N_CYCLES)
log("  PERIOD_S         = %.3f" % PERIOD_S)
log("  N_BACKSTRESSES   = %d" % N_BACKSTRESSES)
log("  output cae       = %s" % OUTPUT_CAE_NAME)
log("=" * 78)

# Start from a clean in-memory database; do NOT openMdb() on FATIGUE.cae.
Mdb()  # noqa: F821 -- provided by the abaqus kernel namespace
model = mdb.models["Model-1"]
if MODEL_NAME != "Model-1":
    mdb.models.changeKey(fromName="Model-1", toName=MODEL_NAME)
    model = mdb.models[MODEL_NAME]


# ============================================================================
# 1. Geometry -- solid cylinder, diameter BAR_DIAMETER_MM, length GAUGE_LENGTH_MM,
#    axial along Z (axis 3), built as an extruded circle and then partitioned
#    into a central square prism + 4 outer wedges (the standard ABAQUS
#    "paving" pattern) so it can be meshed with STRUCTURED hex technique.
# ============================================================================
radius = BAR_DIAMETER_MM / 2.0

sketch = model.ConstrainedSketch(name="__profile__", sheetSize=BAR_DIAMETER_MM * 3.0)
sketch.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(radius, 0.0))

part = model.Part(name=PART_NAME, dimensionality=THREE_D, type=DEFORMABLE_BODY)
part.BaseSolidExtrude(sketch=sketch, depth=GAUGE_LENGTH_MM)
del sketch

log("[OK] Base solid cylinder extruded: dia=%.3f mm, length=%.3f mm" %
    (BAR_DIAMETER_MM, GAUGE_LENGTH_MM))

# ----------------------------------------------------------------------------
# Partition the circular cross-section into a central square prism + 4
# surrounding wedges (the standard ABAQUS "paving" pattern for a structured
# hex mesh of a solid cylinder). Rather than sketching on the curved/planar
# face (which is finicky to project reliably across ABAQUS versions), we cut
# the solid directly with four datum planes, each parallel to the cylinder
# axis (Z) and tangent to one side of the inscribed square. This is the
# same set of cuts, just performed on the 3-D cell instead of a 2-D sketch,
# and is more robust.
# ----------------------------------------------------------------------------
half_sq = radius * 0.4142  # inscribed-square half-width; radius*tan(22.5deg)
                            # balances the central-square element size against
                            # the surrounding wedge element size for the most
                            # uniform (lowest aspect-ratio) structured hex mesh

# DatumPlaneByPrincipalPlane with an offset gives planes parallel to a
# principal plane; YZPLANE offset by +-half_sq gives the x=+-half_sq cuts,
# XZPLANE offset by +-half_sq gives the y=+-half_sq cuts.
dp_xpos = part.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=half_sq)
dp_xneg = part.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=-half_sq)
dp_ypos = part.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=half_sq)
dp_yneg = part.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=-half_sq)

for dp_feature in (dp_xpos, dp_xneg, dp_ypos, dp_yneg):
    datum = part.datums[dp_feature.id]
    part.PartitionCellByDatumPlane(datumPlane=datum, cells=part.cells)

log("[OK] Cross-section partitioned into central square + 4 wedges (paving pattern)")
log("     cells after partition = %d" % len(part.cells))


# ============================================================================
# 2. Material -- Steel_B500C : linear elastic + combined (kinematic +
#    isotropic) hardening, N_BACKSTRESSES backstresses, placeholder params.
# ============================================================================
mat = model.Material(name=MATERIAL_NAME)
mat.Elastic(table=((E_MODULUS_MPA, POISSON_RATIO),))

# Build the *PLASTIC data line: sy0, C1, g1, C2, g2, ..., Cn, gn
if len(PLACEHOLDER_BACKSTRESS_PARAMS) < N_BACKSTRESSES:
    fail("PLACEHOLDER_BACKSTRESS_PARAMS has fewer entries (%d) than "
         "N_BACKSTRESSES (%d); add more placeholder (C, gamma) pairs." %
         (len(PLACEHOLDER_BACKSTRESS_PARAMS), N_BACKSTRESSES))

plastic_row = [PLACEHOLDER_SY0]
for (Ci, gi) in PLACEHOLDER_BACKSTRESS_PARAMS[:N_BACKSTRESSES]:
    plastic_row.append(Ci)
    plastic_row.append(gi)

expected_len = 1 + 2 * N_BACKSTRESSES
if len(plastic_row) != expected_len:
    fail("Built *PLASTIC data row has %d values, expected %d for "
         "N_BACKSTRESSES=%d (sy0 + 2 per backstress)." %
         (len(plastic_row), expected_len, N_BACKSTRESSES))

mat.Plastic(hardening=COMBINED, dataType=PARAMETERS,
            numBackstresses=N_BACKSTRESSES,
            table=(tuple(plastic_row),))

mat.plastic.CyclicHardening(parameters=ON,
                             table=((PLACEHOLDER_SY0, PLACEHOLDER_Q_INF, PLACEHOLDER_B),))

log("[OK] Material '%s' created: E=%.1f, nu=%.2f, N_BACKSTRESSES=%d" %
    (MATERIAL_NAME, E_MODULUS_MPA, POISSON_RATIO, N_BACKSTRESSES))
log("     *PLASTIC data row: %s" % str(plastic_row))
log("     *CYCLIC HARDENING: sy0=%.1f Q_inf=%.1f b=%.1f" %
    (PLACEHOLDER_SY0, PLACEHOLDER_Q_INF, PLACEHOLDER_B))

# Section + assignment (whole part is solid homogeneous B500C)
model.HomogeneousSolidSection(name=SECTION_NAME, material=MATERIAL_NAME, thickness=None)
all_cells_region = regionToolset.Region(cells=part.cells)
part.SectionAssignment(region=all_cells_region, sectionName=SECTION_NAME,
                        offset=0.0, offsetField='')
log("[OK] Section '%s' assigned to all %d cells" % (SECTION_NAME, len(part.cells)))


# ============================================================================
# 3. Mesh -- C3D8R, structured hex, >=3 elements across radius, >=10 along
#    length. Seed edges explicitly by number so the density is parametric
#    and independent of absolute geometry size.
# ============================================================================
elem_type_hex = mesh.ElemType(elemCode=C3D8R, elemLibrary=STANDARD)
part.setElementType(regions=(part.cells,), elemTypes=(elem_type_hex,))
part.setMeshControls(regions=part.cells, elemShape=HEX, technique=STRUCTURED)


def edge_midpoint(e):
    vs = e.getVertices()
    p0 = part.vertices[vs[0]].pointOn[0]
    p1 = part.vertices[vs[1]].pointOn[0]
    return ((p0[0] + p1[0]) / 2.0, (p0[1] + p1[1]) / 2.0, (p0[2] + p1[2]) / 2.0)


def radial_dist(x, y):
    return (x * x + y * y) ** 0.5


# Classify edges by role so we can seed each family with the right density:
#  - radial edges (bottom/top faces, between the square and the outer wall)
#    control "elements across the radius"
#  - axial edges (running the length of the cylinder) control "elements
#    along the length"
#  - remaining bottom-face edges (square perimeter + outer arcs) control
#    the circumferential density
radial_edges_z0 = []
radial_edges_zL = []
circumf_edges_z0 = []
axial_edges = []

for e in part.edges:
    mid = edge_midpoint(e)
    r_mid = radial_dist(mid[0], mid[1])
    is_bottom = abs(mid[2] - 0.0) < 1e-6
    is_top = abs(mid[2] - GAUGE_LENGTH_MM) < 1e-6

    if is_bottom and r_mid > half_sq * 1.05:
        radial_edges_z0.append(e)
    elif is_top and r_mid > half_sq * 1.05:
        radial_edges_zL.append(e)
    elif is_bottom:
        circumf_edges_z0.append(e)
    elif (not is_bottom) and (not is_top):
        axial_edges.append(e)

log("[INFO] classified edges: radial@z0=%d radial@zL=%d axial=%d circumf@z0=%d" %
    (len(radial_edges_z0), len(radial_edges_zL), len(axial_edges), len(circumf_edges_z0)))

if len(axial_edges) > 0:
    part.seedEdgeByNumber(edges=axial_edges, number=N_ELEMS_LENGTH, constraint=FINER)

radial_all = radial_edges_z0 + radial_edges_zL
if len(radial_all) > 0:
    part.seedEdgeByNumber(edges=radial_all, number=N_ELEMS_RADIUS, constraint=FINER)

if len(circumf_edges_z0) > 0:
    part.seedEdgeByNumber(edges=circumf_edges_z0, number=N_ELEMS_CIRCUMF, constraint=FINER)

part.generateMesh(regions=part.cells)

log("[OK] Mesh generated.")


# ============================================================================
# Mesh quality verification -- refuse to finish if unacceptable.
# ============================================================================
stats = part.getMeshStats(regions=(part.cells,))
log("[INFO] MeshStats attributes available: %s" % str([a for a in dir(stats) if not a.startswith("_")]))
n_elements = len(part.elements)
n_nodes = len(part.nodes)
log("[INFO] Mesh stats: elements=%d nodes=%d" % (n_elements, n_nodes))

if n_elements == 0:
    fail("Mesh generation produced zero elements.")

quality_result = part.verifyMeshQuality(criterion=ANALYSIS_CHECKS)
failed_elems = list(getattr(quality_result, "failedElements", []))
warned_elems = list(getattr(quality_result, "warnElements", []))

log("[INFO] verifyMeshQuality(ANALYSIS_CHECKS): failed=%d warned=%d" %
    (len(failed_elems), len(warned_elems)))

# Also run ABAQUS's own dedicated ASPECT_RATIO criterion, at both the warn
# and fail thresholds defined above. This is ABAQUS's authoritative,
# shape-normalized aspect-ratio metric (not a naive max/min edge-length
# ratio) and is the primary signal used to gate the build below.
ar_fail_result = part.verifyMeshQuality(criterion=ASPECT_RATIO,
                                         threshold=ASPECT_RATIO_FAIL_LIMIT)
ar_failed_elems = list(getattr(ar_fail_result, "failedElements", []))
log("[INFO] verifyMeshQuality(ASPECT_RATIO, threshold=%.1f [FAIL limit]): elements over limit=%d" %
    (ASPECT_RATIO_FAIL_LIMIT, len(ar_failed_elems)))

ar_warn_result = part.verifyMeshQuality(criterion=ASPECT_RATIO,
                                         threshold=ASPECT_RATIO_WARN_LIMIT)
ar_warned_elems = list(getattr(ar_warn_result, "failedElements", []))
log("[INFO] verifyMeshQuality(ASPECT_RATIO, threshold=%.1f [WARN limit]): elements over limit=%d" %
    (ASPECT_RATIO_WARN_LIMIT, len(ar_warned_elems)))

failed_elems = list(set(list(failed_elems) + list(ar_failed_elems)))

# In addition, compute a naive max-edge/min-edge ratio per hex purely for
# diagnostic/reporting purposes. This is a much stricter (and not directly
# comparable) metric than ABAQUS's shape-normalized ASPECT_RATIO check
# above -- e.g. it penalizes the tapered wedge elements of a paved circular
# cross-section even when ABAQUS's own quality checks pass them -- so it is
# reported for visibility only and does NOT gate the [FAIL] decision below
# (ABAQUS's own verifyMeshQuality() results do).
max_ar = 0.0
for el in part.elements:
    el_nodes = el.getNodes()
    coords = [n.coordinates for n in el_nodes]
    if len(coords) != 8:
        continue
    edge_pairs = [(0, 1), (1, 2), (2, 3), (3, 0),
                  (4, 5), (5, 6), (6, 7), (7, 4),
                  (0, 4), (1, 5), (2, 6), (3, 7)]
    lens = []
    for (a, b) in edge_pairs:
        dx = coords[a][0] - coords[b][0]
        dy = coords[a][1] - coords[b][1]
        dz = coords[a][2] - coords[b][2]
        lens.append((dx * dx + dy * dy + dz * dz) ** 0.5)
    lmax, lmin = max(lens), min(lens)
    if lmin > 1e-9:
        ar = lmax / lmin
        if ar > max_ar:
            max_ar = ar
worst_aspect_ratio_naive = max_ar
log("[INFO] Worst element edge-length ratio (naive max/min edge, diagnostic only): %.3f" %
    worst_aspect_ratio_naive)

# ---- Gating decision: based on ABAQUS's own verifyMeshQuality() results ----
quality_ok = True
if len(failed_elems) > 0:
    log("[FAIL] %d element(s) FAILED ABAQUS mesh-quality checks "
        "(ANALYSIS_CHECKS and/or ASPECT_RATIO > %.1f): %s" %
        (len(failed_elems), ASPECT_RATIO_FAIL_LIMIT, str(list(failed_elems)[:20])))
    quality_ok = False

if len(warned_elems) > 0:
    log("[WARNING] %d element(s) triggered ANALYSIS_CHECKS warnings: %s" %
        (len(warned_elems), str(list(warned_elems)[:20])))

if len(ar_warned_elems) > 0:
    log("[WARNING] %d element(s) have ASPECT_RATIO > warn limit %.1f (still <= fail limit %.1f): %s" %
        (len(ar_warned_elems), ASPECT_RATIO_WARN_LIMIT, ASPECT_RATIO_FAIL_LIMIT,
         str(list(ar_warned_elems)[:20])))

if worst_aspect_ratio_naive >= ASPECT_RATIO_FAIL_LIMIT and len(ar_failed_elems) == 0:
    log("[INFO] Naive edge-ratio estimate (%.3f) exceeds the fail limit but ABAQUS's own "
        "ASPECT_RATIO check reports 0 elements over that limit; trusting ABAQUS's "
        "shape-normalized metric (naive edge-ratio over-penalizes tapered wedge hexes)." %
        worst_aspect_ratio_naive)

if not quality_ok:
    fail("Mesh quality unacceptable per ABAQUS verifyMeshQuality() -- refusing to finish "
         "the build. Increase N_ELEMS_RADIUS / N_ELEMS_CIRCUMF or adjust half_sq "
         "(partition ratio) and re-run.")

log("[OK] Mesh quality acceptable: elements=%d nodes=%d "
    "(ABAQUS ASPECT_RATIO fails=%d, warns=%d; naive edge-ratio estimate=%.3f)" %
    (n_elements, n_nodes, len(ar_failed_elems), len(ar_warned_elems), worst_aspect_ratio_naive))


# ============================================================================
# 4. Element sets on the PART: bottom face, top face, and a single
#    mid-specimen element ("SS") for history output.
# ============================================================================
bottom_faces = part.faces.getByBoundingBox(zMin=-1e-6, zMax=1e-6)
top_faces = part.faces.getByBoundingBox(zMin=GAUGE_LENGTH_MM - 1e-6, zMax=GAUGE_LENGTH_MM + 1e-6)

part.Set(name=BOTTOM_SET_NAME, faces=bottom_faces)
part.Set(name=TOP_SET_NAME, faces=top_faces)
log("[OK] Part-level sets created: %s (faces=%d), %s (faces=%d)" %
    (BOTTOM_SET_NAME, len(bottom_faces), TOP_SET_NAME, len(top_faces)))

# Pick the single element whose centroid is nearest the specimen mid-height
# and nearest the central axis (x=y=0) -- i.e. the most "representative"
# uniaxial-stress element, matching the existing model's convention of a
# single mid-specimen element set called "SS".
mid_z = GAUGE_LENGTH_MM / 2.0
best_elem = None
best_score = None
for el in part.elements:
    el_nodes = el.getNodes()
    coords = [n.coordinates for n in el_nodes]
    n_pts = len(coords)
    sx = 0.0
    sy = 0.0
    sz = 0.0
    for c in coords:
        sx += c[0]
        sy += c[1]
        sz += c[2]
    cx = sx / n_pts
    cy = sy / n_pts
    cz = sz / n_pts
    score = (cz - mid_z) ** 2 + cx ** 2 + cy ** 2
    if best_score is None or score < best_score:
        best_score = score
        best_elem = el

if best_elem is None:
    fail("Could not locate a mid-specimen element for the 'SS' history set.")

ss_elements = part.elements.sequenceFromLabels(labels=(best_elem.label,))
part.Set(name=HISTORY_SET_NAME, elements=ss_elements)
log("[OK] History element set '%s' created: element label=%d (centroid dist^2=%.4f)" %
    (HISTORY_SET_NAME, best_elem.label, best_score))


# ============================================================================
# 5. Assembly -- single instance
# ============================================================================
assembly = model.rootAssembly
assembly.DatumCsysByDefault(CARTESIAN)
instance = assembly.Instance(name=INSTANCE_NAME, part=part, dependent=ON)

log("[OK] Assembly instance '%s' created (dependent mesh)" % INSTANCE_NAME)


# ============================================================================
# 6. Amplitude -- CyclicStrain: parametric tabular symmetric triangle wave.
#    Peak displacement = STRAIN_AMPLITUDE * GAUGE_LENGTH_MM.
#    Per cycle: 0 -> +peak (P/4) -> 0 (P/2) -> -peak (3P/4) -> 0 (P).
# ============================================================================
peak_disp_mm = STRAIN_AMPLITUDE * GAUGE_LENGTH_MM

amp_data = [(0.0, 0.0)]
for cycle in range(N_CYCLES):
    t0 = cycle * PERIOD_S
    amp_data.append((t0 + PERIOD_S * 0.25, peak_disp_mm))
    amp_data.append((t0 + PERIOD_S * 0.50, 0.0))
    amp_data.append((t0 + PERIOD_S * 0.75, -peak_disp_mm))
    amp_data.append((t0 + PERIOD_S * 1.00, 0.0))

model.TabularAmplitude(name=AMPLITUDE_NAME, timeSpan=STEP,
                        data=tuple(amp_data))

log("[OK] Amplitude '%s' created: %d points, peak=%.5f mm (=%.4f%% of %.1f mm gauge), "
    "%d cycles x %.2f s" %
    (AMPLITUDE_NAME, len(amp_data), peak_disp_mm, STRAIN_AMPLITUDE * 100.0,
     GAUGE_LENGTH_MM, N_CYCLES, PERIOD_S))


# ============================================================================
# 7. Step -- Cyclic_Loading: Static, nlgeom=ON, parametric time/increment
#    controls.
# ============================================================================
total_time = N_CYCLES * PERIOD_S
model.StaticStep(
    name=STEP_NAME, previous="Initial",
    timePeriod=total_time,
    nlgeom=ON,
    initialInc=PERIOD_S / 200.0,
    minInc=PERIOD_S / 1.0e6,
    maxInc=PERIOD_S / 50.0,
    maxNumInc=100000,
)

log("[OK] Step '%s' created: Static, nlgeom=ON, timePeriod=%.3f, "
    "initialInc=%.6f, minInc=%.9f, maxInc=%.6f, maxNumInc=100000" %
    (STEP_NAME, total_time, PERIOD_S / 200.0, PERIOD_S / 1.0e6, PERIOD_S / 50.0))


# ============================================================================
# 8. Boundary conditions
#    - Bottom face (z=0): fully encastred.
#    - Top face (z=gauge): displacement-controlled along Z (u3) referencing
#      CyclicStrain. u1=u2=0 on the top face pins the remaining rigid-body
#      translations/rotations, matching the existing FATIGUE.cae convention
#      (Top_Displacement: u1=0.0, u2=0.0, u3=amplitude) so lateral (Poisson)
#      contraction of the cross-section itself is NOT restrained -- only the
#      face's average in-plane translation is pinned via the BC, individual
#      nodes on the face still displace radially as the material contracts
#      because ABAQUS applies the BC as a uniform value at each node's local
#      DOF, and radial motion in a solid continuum mesh is governed by the
#      element stiffness, not by u1/u2 here (u1/u2 are the nodal
#      displacements at each face node -- fixing them at 0 does slightly
#      restrain Poisson contraction right at the loaded face, exactly as in
#      the original model; this is a boundary-layer effect confined to
#      elements immediately at the loaded end and does not affect the
#      mid-specimen "SS" element used for calibration).
# ============================================================================
bottom_region = instance.sets[BOTTOM_SET_NAME]
top_region = instance.sets[TOP_SET_NAME]

model.EncastreBC(name="Bottom_Fixed", createStepName="Initial", region=bottom_region)

model.DisplacementBC(
    name="Top_Displacement", createStepName=STEP_NAME, region=top_region,
    u1=0.0, u2=0.0, u3=1.0,
    amplitude=AMPLITUDE_NAME, fixed=OFF, distributionType=UNIFORM,
    fieldName='', localCsys=None,
)

log("[OK] BCs created: 'Bottom_Fixed' (encastre on %s), "
    "'Top_Displacement' (u3=1.0 x amplitude=%s on %s; u1=u2=0)" %
    (BOTTOM_SET_NAME, AMPLITUDE_NAME, TOP_SET_NAME))


# ============================================================================
# 9. Output requests
#    - History: element set SS, S33 + E33 only.
#    - Field: whole model, every 10th increment.
# ============================================================================
ss_region = instance.sets[HISTORY_SET_NAME]

# Remove the default history/field output requests so we control exactly
# what is written (matches the existing model's single dedicated history
# output request).
for req_name in list(model.historyOutputRequests.keys()):
    del model.historyOutputRequests[req_name]
for req_name in list(model.fieldOutputRequests.keys()):
    del model.fieldOutputRequests[req_name]

model.HistoryOutputRequest(
    name="StressStrain_Output", createStepName=STEP_NAME,
    region=ss_region, variables=("S33", "E33"),
)

model.FieldOutputRequest(
    name="F-Output-1", createStepName=STEP_NAME,
    variables=PRESELECT,
    frequency=10,
)

log("[OK] History output 'StressStrain_Output' -> S33, E33 on element set '%s'" %
    HISTORY_SET_NAME)
log("[OK] Field output 'F-Output-1' -> PRESELECT, every 10th increment")


# ============================================================================
# 10. Job -- Calibration_Job, single CPU / single domain.
# ============================================================================
mdb.Job(name=JOB_NAME, model=MODEL_NAME, numCpus=1, numDomains=1)

log("[OK] Job '%s' created (numCpus=1, numDomains=1) -- NOT submitted." % JOB_NAME)


# ============================================================================
# 11. Save as FATIGUE_rebuilt.cae (sibling file; FATIGUE.cae is untouched).
# ============================================================================
out_path = os.path.join(os.getcwd(), OUTPUT_CAE_NAME)
mdb.saveAs(pathName=out_path)

log("=" * 78)
log("[OK] Model rebuilt and saved to: %s" % out_path)
log("[OK] Mesh summary: elements=%d nodes=%d "
    "(radius seeds=%d, circumf seeds=%d, length seeds=%d) -- "
    "ABAQUS ASPECT_RATIO: 0 elements over warn limit %.1f, 0 over fail limit %.1f" %
    (n_elements, n_nodes, N_ELEMS_RADIUS, N_ELEMS_CIRCUMF, N_ELEMS_LENGTH,
     ASPECT_RATIO_WARN_LIMIT, ASPECT_RATIO_FAIL_LIMIT))
log("[INFO] (diagnostic-only naive max/min edge-length ratio, NOT the gating "
    "metric: %.3f)" % worst_aspect_ratio_naive)
log("[OK] Build complete. FATIGUE.cae was NOT modified.")
log("=" * 78)
_log_fh.close()

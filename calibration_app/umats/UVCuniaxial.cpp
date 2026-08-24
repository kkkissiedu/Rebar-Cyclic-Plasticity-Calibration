/* ==========================================================================
 * UVCuniaxial.cpp - Updated Voce-Chaboche uniaxial UMAT (C++ / MSVC build)
 * ==========================================================================
 *
 * C++ port of UVCuniaxial.for (Hartloper, de Castro e Sousa & Lignos 2021,
 * "Constitutive Modeling of Structural Steels...", J. Struct. Eng.,
 * doi:10.1061/(ASCE)ST.1943-541X.0002964; original UMAT MIT-licensed,
 * github.com/ahartloper/UVC_MatMod - see UVC_LICENSE_MIT.txt).
 *
 * WHY C++: this machine has no Intel Fortran; ABAQUS 2024 compiles C++ user
 * subroutines with MSVC (cl) via the documented aba_for_c.h interface. The
 * job work directory must contain the abaqus_v6.env link_sl override written
 * by core/abaqus_runner.py (no Intel runtime import libs on this machine).
 *
 * Constitutive model (uniaxial, NTENS = 1, e.g. T3D2 truss):
 *   Yield:      phi = (sigma - alpha)^2 - sigma_y^2        (paper Eq. 1-2)
 *   Isotropic:  sigma_y(p) = sy0 + Q(1-e^{-b p}) - D(1-e^{-a p})   (Eq. 3;
 *               second Voce term = yield-plateau extension, D=0 disables)
 *   Kinematic:  Chaboche decomposition, closed-form update over the plastic
 *               increment (Eq. 4-5):
 *               alpha_k = s*C_k/g_k - (s*C_k/g_k - alpha_k^0) e^{-g_k dp}
 *
 * Return mapping: fixed-branch Newton on the LINEAR residual f = |yr| - sy
 * with the branch sign s = sign(yr_trial) frozen at the elastic trial -
 * the same scheme as the app's validated surrogate
 * (core/uvc_model.py::_uvc_integrate_impl, correctness fix of 2026-07-03).
 * It converges for coarse increments where the original .for phi-squared
 * Newton needed ABAQUS substepping, and reaches the identical converged
 * state (surrogate verified < 1 MPa vs substepped UMAT over 500 param sets).
 * On non-convergence PNEWDT = 0.25 requests a smaller increment, matching
 * the original UMAT's behaviour.
 *
 * PROPS (nprops = 6 + 2N, matches uvc_model.material_block):
 *   1: E   2: sy0   3: Q   4: b   5: D   6: a   7..: C_k, gamma_k pairs
 * STATEV (nstatv = 1 + N):
 *   1: equivalent plastic strain p, 2..: alpha_k
 * ========================================================================== */

#include <aba_for_c.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* ---- Intel Fortran runtime shims -------------------------------------
 * standardU_static.lib's stub objects (ifort-compiled) reference three
 * libifcoremd symbols on their error/string-formatting paths. No Intel
 * runtime exists on this machine; the paths never execute when the model
 * only uses the UMAT supplied here - abort loudly if one is ever called. */
extern "C" void c_f_pointer_set_scalar(void)
{
    fprintf(stderr, "FATAL: Intel RTL shim c_f_pointer_set_scalar called\n");
    abort();
}
extern "C" void for_trim(void)
{
    fprintf(stderr, "FATAL: Intel RTL shim for_trim called\n");
    abort();
}
extern "C" void for_concat(void)
{
    fprintf(stderr, "FATAL: Intel RTL shim for_concat called\n");
    abort();
}

#define UVC_MAX_BACK 8
#define UVC_MAXIT    200

/* NOTE: ABAQUS's job driver detects C++ user subroutines with the regex
 *   ^\s*extern.*"C"\s*void\s*FOR_NAME\((\w+)
 * so the declaration below must keep exactly this one-line shape. */
extern "C" void FOR_NAME(umat, UMAT)(
    double *stress, double *statev, double *ddsdde, double *sse,
    double *spd, double *scd, double *rpl, double *ddsddt,
    double *drplde, double *drpldt, double *stran, double *dstran,
    double *time, double *dtime, double *temp, double *dtemp,
    double *predef, double *dpred, char *cmname, int *ndi, int *nshr,
    int *ntens, int *nstatv, double *props, int *nprops,
    double *coords, double *drot, double *pnewdt, double *celent,
    double *dfgrd0, double *dfgrd1, int *noel, int *npt, int *layer,
    int *kspt, int *jstep, int *kinc, const Length cmname_len)
{
    const double E   = props[0];
    const double sy0 = props[1];
    const double Q   = props[2];
    const double b   = props[3];
    const double D   = props[4];
    const double a   = props[5];
    const int    nb  = (*nprops - 6) / 2;

    double C[UVC_MAX_BACK], g[UVC_MAX_BACK];
    double ak[UVC_MAX_BACK], ak0[UVC_MAX_BACK];
    int k;

    if (nb < 1 || nb > UVC_MAX_BACK || *nstatv < 1 + nb) {
        fprintf(stderr, "UVC UMAT: bad nprops=%d/nstatv=%d\n",
                *nprops, *nstatv);
        abort();
    }
    for (k = 0; k < nb; ++k) {
        C[k]   = props[6 + 2 * k];
        g[k]   = props[7 + 2 * k];
        ak[k]  = statev[1 + k];
        ak0[k] = ak[k];
    }

    const double p0 = statev[0];

    /* elastic trial (UVCuniaxial.for L62) */
    const double sig0 = stress[0] + E * dstran[0];
    double sy = sy0 + Q * (1.0 - exp(-b * p0)) - D * (1.0 - exp(-a * p0));
    double alpha = 0.0;
    for (k = 0; k < nb; ++k) alpha += ak[k];
    double yr = sig0 - alpha;
    const double phi = yr * yr - sy * sy;

    if (phi <= 1.0e-10) {                       /* elastic step */
        stress[0] = sig0;
        ddsdde[0] = E;
        return;
    }

    /* plastic: fixed-branch Newton on f = |yr| - sy (branch s frozen) */
    const double s = (yr >= 0.0) ? 1.0 : -1.0;
    double dp = 0.0, pp = p0, al = 0.0, aux = E;
    int it, converged = 0;
    for (it = 0; it < UVC_MAXIT; ++it) {
        const double sig_c = sig0 - E * s * dp;
        al = 0.0;
        aux = E;                                /* d|yr|/ddp accumulator */
        for (k = 0; k < nb; ++k) {
            double akk;
            if (g[k] > 1.0e-12) {
                const double sat = s * C[k] / g[k];
                akk = sat - (sat - ak0[k]) * exp(-g[k] * dp);
                aux += C[k] - s * g[k] * akk;
            } else {                            /* g -> 0 linear limit */
                akk = ak0[k] + s * C[k] * dp;
                aux += C[k];
            }
            ak[k] = akk;
            al += akk;
        }
        pp = p0 + dp;
        sy = sy0 + Q * (1.0 - exp(-b * pp)) - D * (1.0 - exp(-a * pp));
        yr = sig_c - al;
        const double f = s * yr - sy;           /* residual = |yr| - sy */
        if (fabs(f) < 1.0e-9 * (sy0 + fabs(yr))) { converged = 1; break; }
        const double dfa = aux + Q * b * exp(-b * pp) - D * a * exp(-a * pp);
        if (dfa <= 0.0) break;                  /* softening pathologies */
        double dp_new = dp + f / dfa;
        if (dp_new < 0.0) dp_new = 0.5 * dp;
        dp = dp_new;
    }

    if (!converged) {                           /* ask ABAQUS to cut back */
        *pnewdt = 0.25;
        stress[0] = sig0 - E * s * dp;
        ddsdde[0] = E;
        return;
    }

    /* commit state */
    stress[0] = sig0 - E * s * dp;
    statev[0] = pp;
    for (k = 0; k < nb; ++k) statev[1 + k] = ak[k];

    /* consistent tangent, same closed form as UVCuniaxial.for L152-160:
     * A = b(Q - isoQ) - a(D - isoD) + sum_k (C_k - s g_k alpha_k);
     * ddsdde = E A / (E + A). aux already holds E + sum(...) so
     * A = (aux - E) + hardening-rate terms. */
    {
        const double A = (aux - E)
                       + Q * b * exp(-b * pp) - D * a * exp(-a * pp);
        const double den = E + A;
        ddsdde[0] = (fabs(den) > 1.0e-8 * E) ? (E * A / den) : E;
    }
}

/* ==========================================================================
 * OhnoWang.cpp — Ohno-Wang model I uniaxial UMAT (C++ / MSVC build)
 * ==========================================================================
 *
 * Ohno & Wang 1993, "Kinematic hardening rules with critical state of
 * dynamic recovery, Part I", Int. J. Plasticity 9(3):375-390,
 * doi:10.1016/0749-6419(93)90042-O. Optional Voce isotropic extension
 * (Abdel-Karim 2010, doi:10.1016/j.ijpvp.2010.02.003); Q = 0 recovers the
 * pure OW-I 1993 model. Written for this project (no permissive OW UMAT
 * was available); mirrors the app's validated surrogate integrator
 * core/ohno_wang_model.py::_ow_integrate_impl exactly, plus internal
 * substepping so ABAQUS's coarser auto-increments match the surrogate's
 * fine experimental sampling.
 *
 * Constitutive model (uniaxial, NTENS = 1, e.g. T3D2 truss):
 *   Yield/flow: f = |sigma - alpha| - sigma_y <= 0,
 *               d eps_p = dp * n, n = sign(sigma - alpha)   (OW93 Eq. 1-2)
 *   Backstress: alpha = sum alpha_i                         (OW93 Eq. 5)
 *               dalpha_i = C_i n dp
 *                 - C_i (|alpha_i|/r_i)^{m_i} <n sgn(alpha_i)> sgn(alpha_i) dp
 *                                                           (OW93 Eq. 10)
 *   Isotropic:  sigma_y(p) = sy0 + Q(1 - e^{-b p})
 *
 * Integration: per substep, radial return with flow direction fixed at the
 * elastic trial; dp solved by Newton on isotropic consistency with the
 * backstress frozen; then one explicit Eq. 10 backstress step — identical
 * to the surrogate. Substep size limit 1e-4 strain (OW_SUBSTEP_EPS, same
 * constant as the surrogate's _OW_SUBSTEP_EPS) keeps the explicit
 * backstress integration within ~2 MPa of the converged ODE solution.
 *
 * WHY C++: no Intel Fortran on this machine; ABAQUS 2024 compiles C++ user
 * subroutines with MSVC (cl). Job dir needs the abaqus_v6.env link_sl
 * override written by core/abaqus_runner.py.
 *
 * PROPS (nprops = 4 + 3N, matches ohno_wang_model.material_block):
 *   1: E   2: sy0   3: Q   4: b   5..: C_k, r_k, m_k triplets
 * STATEV (nstatv = 1 + N):
 *   1: equivalent plastic strain p, 2..: alpha_k
 * ========================================================================== */

#include <aba_for_c.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* ---- Intel Fortran runtime shims (see UVCuniaxial.cpp for rationale) --- */
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

#define OW_MAX_BACK    8
#define OW_NEWTON_ITER 30
#define OW_TOL         1.0e-9
#define OW_SUBSTEP_EPS 1.0e-4
#define OW_MAX_SUBSTEP 1000

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
    const int    nb  = (*nprops - 4) / 3;

    double C[OW_MAX_BACK], r[OW_MAX_BACK], m[OW_MAX_BACK];
    double ak[OW_MAX_BACK];
    int k;

    if (nb < 1 || nb > OW_MAX_BACK || *nstatv < 1 + nb) {
        fprintf(stderr, "OW UMAT: bad nprops=%d/nstatv=%d\n",
                *nprops, *nstatv);
        abort();
    }
    for (k = 0; k < nb; ++k) {
        C[k]  = props[4 + 3 * k];
        r[k]  = props[5 + 3 * k];
        m[k]  = props[6 + 3 * k];
        ak[k] = statev[1 + k];
    }

    double p     = statev[0];
    double sigma = stress[0];

    /* substep the strain increment so the explicit Eq. 10 backstress
     * update sees increments no coarser than the surrogate's data path */
    int nsub = (int)ceil(fabs(dstran[0]) / OW_SUBSTEP_EPS);
    if (nsub < 1) nsub = 1;
    if (nsub > OW_MAX_SUBSTEP) nsub = OW_MAX_SUBSTEP;
    const double de = dstran[0] / (double)nsub;

    int any_plastic = 0;
    double A_last = E;                          /* hardening modulus */

    for (int isub = 0; isub < nsub; ++isub) {
        const double sig_tr = sigma + E * de;
        double alpha = 0.0;
        for (k = 0; k < nb; ++k) alpha += ak[k];
        const double eta = sig_tr - alpha;
        double sy = sy0 + Q * (1.0 - exp(-b * p));

        if (fabs(eta) - sy <= 0.0) {
            sigma = sig_tr;                     /* elastic substep */
            continue;
        }

        const double ndir    = (eta >= 0.0) ? 1.0 : -1.0;
        const double eta_abs = fabs(eta);

        /* Newton for dp, backstress frozen (isotropic consistency):
         * F(dp) = eta_abs - E dp - sy0 - Q(1 - e^{-b(p+dp)}) */
        double dp = (eta_abs - sy) / E;
        if (dp < 1.0e-12) dp = 1.0e-12;
        for (int it = 0; it < OW_NEWTON_ITER; ++it) {
            const double ex = exp(-b * (p + dp));
            const double F  = eta_abs - E * dp - sy0 - Q * (1.0 - ex);
            if (fabs(F) < OW_TOL * (sy0 + eta_abs)) break;
            const double dF = -E - Q * b * ex;
            if (dF == 0.0) break;
            double dp_new = dp - F / dF;
            if (dp_new <= 0.0) dp_new = dp * 0.5;
            dp = dp_new;
        }

        p += dp;
        sigma = sig_tr - E * ndir * dp;
        any_plastic = 1;

        /* explicit Ohno-Wang I backstress update (Eq. 10) + tangent accum */
        A_last = Q * b * exp(-b * p);
        for (k = 0; k < nb; ++k) {
            const double a0 = ak[k];
            const double sk = (a0 > 0.0) ? 1.0 : ((a0 < 0.0) ? -1.0 : 0.0);
            const double recov =
                (r[k] > 0.0) ? pow(fabs(a0) / r[k], m[k]) : 0.0;
            const double active = ((ndir * sk) > 0.0) ? 1.0 : 0.0;
            ak[k] = a0 + C[k] * ndir * dp - C[k] * recov * active * sk * dp;
            A_last += C[k] * (1.0 - recov * active);
        }
    }

    stress[0] = sigma;
    statev[0] = p;
    for (k = 0; k < nb; ++k) statev[1 + k] = ak[k];

    if (!any_plastic) {
        ddsdde[0] = E;
    } else {
        /* elastoplastic tangent E*A/(E+A); A may be negative (cyclic
         * softening, Q < 0) — guard the denominator, floor the magnitude
         * so ABAQUS's Newton stays bounded */
        const double den = E + A_last;
        double t = (fabs(den) > 1.0e-3 * E) ? (E * A_last / den)
                                            : 1.0e-3 * E;
        if (t > E) t = E;
        if (t < -E) t = -E;
        if (fabs(t) < 1.0e-4 * E) t = (t >= 0.0 ? 1.0 : -1.0) * 1.0e-4 * E;
        ddsdde[0] = t;
    }
}

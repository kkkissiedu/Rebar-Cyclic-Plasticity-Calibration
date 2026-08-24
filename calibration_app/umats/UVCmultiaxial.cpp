/* ==========================================================================
 * UVCmultiaxial.cpp — Updated Voce-Chaboche 3D UMAT (C++ / MSVC build)
 * ==========================================================================
 *
 * C++ port of UVCmultiaxial.for (Hartloper, de Castro e Sousa & Lignos 2021,
 * J. Struct. Eng., doi:10.1061/(ASCE)ST.1943-541X.0002964; original UMAT
 * MIT-licensed, github.com/ahartloper/UVC_MatMod — UVC_LICENSE_MIT.txt).
 * Used by "Transfer to CAE" for the full C3D8R coupon (the FE search backend
 * uses the uniaxial UVCuniaxial.cpp on a T3D2 truss instead — faster and
 * exactly equal to the surrogate).
 *
 * Constitutive model (J2 plasticity, radial return; paper Eq. 1-5):
 *   Yield:      f = ||s - a|| - sqrt(2/3) sigma_y(p)
 *   Isotropic:  sigma_y(p) = sy0 + Q(1-e^{-b p}) - D(1-e^{-a p})
 *   Kinematic:  Chaboche backstresses, closed-form update over the plastic
 *               increment with fixed flow normal n (as in the .for):
 *               a_k <- e_k a_k + sqrt(2/3) (C_k/g_k)(1-e_k) n,
 *               e_k = exp(-g_k (p - p0))
 *   Newton on the plastic multiplier (UVCmultiaxial.for L164-204).
 *
 * Differences vs the .for (deliberate, robustness-neutral):
 *   - state tensors are rotated by DROT with an explicit R T R^T (the .for
 *     calls the ABAQUS utility ROTSIG; avoiding it keeps this file free of
 *     Fortran-runtime link dependencies);
 *   - gamma_k ~ 0 handled by the linear limit instead of dividing by g_k.
 *
 * PROPS (nprops = 7 + 2N — note nu, unlike the uniaxial card):
 *   1: E  2: nu  3: sy0  4: Q  5: b  6: D  7: a  8..: C_k, gamma_k pairs
 * STATEV (nstatv = 7 + 6N):
 *   1: p, 2-7: plastic strain (engineering shears), 8..: a_k (6 each)
 *
 * Voigt order (ABAQUS): 11, 22, 33, 12, 13, 23; strains use engineering
 * shear, stresses tensor shear. dotprod6 doubles the shear terms (tensor
 * contraction of symmetric tensors stored as 6-vectors).
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

#define UVC_MAX_BACK 8
#define UVC_MAXIT    1000
#define UVC_TOL      1.0e-10

/* tensor contraction of two symmetric 6-vectors (shears doubled) */
static double dotprod6(const double *A, const double *B)
{
    return A[0] * B[0] + A[1] * B[1] + A[2] * B[2]
         + 2.0 * (A[3] * B[3] + A[4] * B[4] + A[5] * B[5]);
}

/* Rotate a symmetric tensor 6-vector by the (Fortran column-major) 3x3
 * incremental rotation DROT: T' = R T R^T. is_strain: shears stored as
 * engineering (2x tensor) components. Order: 11,22,33,12,13,23. */
static void rot6(const double *drot, double *v, int is_strain)
{
    const double half = is_strain ? 0.5 : 1.0;
    double T[3][3], RT[3][3], out[3][3];
    int i, j, k;
    T[0][0] = v[0]; T[1][1] = v[1]; T[2][2] = v[2];
    T[0][1] = T[1][0] = v[3] * half;
    T[0][2] = T[2][0] = v[4] * half;
    T[1][2] = T[2][1] = v[5] * half;
    /* R(i,j) = drot[j*3 + i] (Fortran column-major) */
    for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j) {
            RT[i][j] = 0.0;
            for (k = 0; k < 3; ++k)
                RT[i][j] += drot[k * 3 + i] * T[k][j];
        }
    for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j) {
            out[i][j] = 0.0;
            for (k = 0; k < 3; ++k)
                out[i][j] += RT[i][k] * drot[k * 3 + j];
        }
    v[0] = out[0][0]; v[1] = out[1][1]; v[2] = out[2][2];
    v[3] = out[0][1] / half; v[4] = out[0][2] / half;
    v[5] = out[1][2] / half;
}

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
    const double SQRT23 = sqrt(2.0 / 3.0);
    const double E    = props[0];
    const double nu   = props[1];
    const double sy0  = props[2];
    const double Q    = props[3];
    const double b    = props[4];
    const double D    = props[5];
    const double a    = props[6];
    const int    nb   = (*nprops - 7) / 2;
    const int    nt   = *ntens;        /* 6 for C3D8R */

    double C[UVC_MAX_BACK], g[UVC_MAX_BACK];
    double ak[UVC_MAX_BACK][6];
    int i, j, k;

    if (nt != 6 || nb < 1 || nb > UVC_MAX_BACK || *nstatv < 7 + 6 * nb) {
        fprintf(stderr, "UVC 3D UMAT: bad ntens=%d/nprops=%d/nstatv=%d\n",
                nt, *nprops, *nstatv);
        abort();
    }
    for (k = 0; k < nb; ++k) {
        C[k] = props[7 + 2 * k];
        g[k] = props[8 + 2 * k];
    }

    /* state: rotate plastic strain + backstresses by the increment rotation */
    const double p0 = statev[0];
    double eps_p[6];
    for (i = 0; i < 6; ++i) eps_p[i] = statev[1 + i];
    rot6(drot, eps_p, 1);
    double alpha[6] = {0, 0, 0, 0, 0, 0};
    for (k = 0; k < nb; ++k) {
        for (i = 0; i < 6; ++i) ak[k][i] = statev[7 + 6 * k + i];
        rot6(drot, ak[k], 0);
        for (i = 0; i < 6; ++i) alpha[i] += ak[k][i];
    }

    /* elastic moduli tensor (UVCmultiaxial.for L118-125) */
    const double mu   = E / (2.0 * (1.0 + nu));
    const double bulk = E / (3.0 * (1.0 - 2.0 * nu));
    const double mu2  = 2.0 * mu;
    double cmat[6][6];
    for (i = 0; i < 6; ++i)
        for (j = 0; j < 6; ++j) {
            const double id2 = (i < 3 && j < 3) ? 1.0 : 0.0;
            const double id4 = (i == j) ? (i < 3 ? 1.0 : 0.5) : 0.0;
            cmat[i][j] = id2 * bulk + mu2 * (id4 - (1.0 / 3.0) * id2);
        }

    /* elastic trial */
    double sig[6];
    for (i = 0; i < 6; ++i) {
        sig[i] = stress[i];
        for (j = 0; j < 6; ++j) sig[i] += cmat[i][j] * dstran[j];
    }
    const double hyd = (sig[0] + sig[1] + sig[2]) / 3.0;
    double srel[6];
    for (i = 0; i < 6; ++i)
        srel[i] = sig[i] - (i < 3 ? hyd : 0.0) - alpha[i];
    const double srn = sqrt(dotprod6(srel, srel));

    double iso_Q = Q * (1.0 - exp(-b * p0));
    double iso_D = D * (1.0 - exp(-a * p0));
    double sy = sy0 + iso_Q - iso_D;
    const double fyield = srn - SQRT23 * sy;

    double n[6];
    for (i = 0; i < 6; ++i) n[i] = srel[i] / (UVC_TOL + srn);

    if (fyield <= UVC_TOL) {                     /* elastic step */
        for (i = 0; i < 6; ++i) stress[i] = sig[i];
        for (i = 0; i < 6; ++i)
            for (j = 0; j < 6; ++j) ddsdde[j * 6 + i] = cmat[i][j];
        statev[0] = p0;
        for (i = 0; i < 6; ++i) statev[1 + i] = eps_p[i];
        for (k = 0; k < nb; ++k)
            for (i = 0; i < 6; ++i) statev[7 + 6 * k + i] = ak[k][i];
        return;
    }

    /* Newton on the plastic multiplier (UVCmultiaxial.for L164-204) */
    double dl = 0.0, p = p0;
    double iso_mod = 0.0, kin_mod = 0.0;
    double alpha_upd[6];
    int it, converged = 0;
    for (it = 0; it < UVC_MAXIT && !converged; ++it) {
        iso_Q = Q * (1.0 - exp(-b * p));
        iso_D = D * (1.0 - exp(-a * p));
        sy = sy0 + iso_Q - iso_D;
        iso_mod = b * (Q - iso_Q) - a * (D - iso_D);
        kin_mod = 0.0;
        for (i = 0; i < 6; ++i) alpha_upd[i] = 0.0;
        for (k = 0; k < nb; ++k) {
            if (g[k] > 1.0e-12) {
                const double ek = exp(-g[k] * (p - p0));
                kin_mod += C[k] * ek - sqrt(1.5) * g[k] * ek
                           * dotprod6(n, ak[k]);
                for (i = 0; i < 6; ++i)
                    alpha_upd[i] += ek * ak[k][i]
                        + SQRT23 * (C[k] / g[k]) * (1.0 - ek) * n[i];
            } else {                             /* g -> 0 linear limit */
                kin_mod += C[k];
                for (i = 0; i < 6; ++i)
                    alpha_upd[i] += ak[k][i]
                        + SQRT23 * C[k] * (p - p0) * n[i];
            }
        }
        double da_n;
        {
            double diff[6];
            for (i = 0; i < 6; ++i) diff[i] = alpha_upd[i] - alpha[i];
            da_n = dotprod6(diff, n);
        }
        const double numer = srn - (da_n + SQRT23 * sy + mu2 * dl);
        const double denom = -mu2 * (1.0 + (kin_mod + iso_mod) / (3.0 * mu));
        dl = dl - numer / denom;
        p = p0 + SQRT23 * dl;
        if (fabs(numer) < UVC_TOL) converged = 1;
    }

    if (!converged) {
        *pnewdt = 0.25;                          /* cut the increment */
    }

    /* commit: plastic strain (engineering shears doubled), stress */
    double dpe[6];
    for (i = 0; i < 6; ++i) dpe[i] = dl * n[i] * (i < 3 ? 1.0 : 2.0);
    for (i = 0; i < 6; ++i) eps_p[i] += dpe[i];
    for (i = 0; i < 6; ++i) {
        double corr = 0.0;
        for (j = 0; j < 6; ++j) corr += cmat[i][j] * dpe[j];
        stress[i] = sig[i] - corr;
    }

    /* backstress closed-form update (same expressions as the Newton) */
    double alpha_new_total[6] = {0, 0, 0, 0, 0, 0};
    for (k = 0; k < nb; ++k) {
        if (g[k] > 1.0e-12) {
            const double ek = exp(-g[k] * (p - p0));
            for (i = 0; i < 6; ++i)
                ak[k][i] = ek * ak[k][i]
                    + SQRT23 * (C[k] / g[k]) * (1.0 - ek) * n[i];
        } else {
            for (i = 0; i < 6; ++i)
                ak[k][i] += SQRT23 * C[k] * (p - p0) * n[i];
        }
        for (i = 0; i < 6; ++i) alpha_new_total[i] += ak[k][i];
    }
    double alpha_diff[6];
    for (i = 0; i < 6; ++i) alpha_diff[i] = alpha_new_total[i] - alpha[i];

    /* consistent tangent (UVCmultiaxial.for L243-263) */
    {
        const double beta = 1.0 + (kin_mod + iso_mod) / (3.0 * mu);
        const double th1 = 1.0 - mu2 * dl / srn;
        const double th3 = 1.0 / (beta * srn);
        const double th2 = 1.0 / beta + dotprod6(n, alpha_diff) * th3
                         - (1.0 - th1);
        double dd[6][6];
        for (i = 0; i < 6; ++i)
            for (j = 0; j < 6; ++j) {
                const double id2 = (i < 3 && j < 3) ? 1.0 : 0.0;
                const double id4 = (i == j) ? (i < 3 ? 1.0 : 0.5) : 0.0;
                dd[i][j] = bulk * id2
                         + mu2 * th1 * (id4 - (1.0 / 3.0) * id2)
                         - mu2 * th2 * n[i] * n[j]
                         + mu2 * th3 * alpha_diff[i] * n[j];
            }
        for (i = 0; i < 6; ++i)
            for (j = 0; j < 6; ++j)
                ddsdde[j * 6 + i] = 0.5 * (dd[i][j] + dd[j][i]);
    }

    /* state variables */
    statev[0] = p;
    for (i = 0; i < 6; ++i) statev[1 + i] = eps_p[i];
    for (k = 0; k < nb; ++k)
        for (i = 0; i < 6; ++i) statev[7 + 6 * k + i] = ak[k][i];
}

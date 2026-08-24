/* ==========================================================================
 * OhnoWangMultiaxial.cpp — Ohno-Wang model I 3D UMAT (C++ / MSVC build)
 * ==========================================================================
 *
 * Ohno & Wang 1993, Int. J. Plasticity 9(3):375-390,
 * doi:10.1016/0749-6419(93)90042-O (Eq. 1-2 flow, Eq. 5 decomposition,
 * Eq. 10 kinematic rule) + optional Voce isotropic extension (Abdel-Karim
 * 2010, doi:10.1016/j.ijpvp.2010.02.003; Q = 0 recovers pure OW-I).
 * Written for this project; used by "Transfer to CAE" for the full C3D8R
 * coupon (the FE search backend uses the uniaxial OhnoWang.cpp on a T3D2
 * truss — exactly equal to the surrogate).
 *
 * Multiaxial kinematic rule, written so that its uniaxial reduction is
 * IDENTICAL to the calibrated uniaxial law in core/ohno_wang_model.py
 * (da = C dp - C (|a|/r)^m <..> dp on the axial stress-space component):
 *   da_i = (2/3) C_i de_p
 *          - C_i (abar_i/r_i)^{m_i} <de_p : nhat_i> nhat_i
 *   abar_i = sqrt(3/2 a_i:a_i),  nhat_i = a_i/abar_i,  <.> = Macaulay.
 *
 * Integration mirrors the uniaxial UMAT/surrogate scheme, tensorised:
 * per substep (equivalent-strain cap OW_SUBSTEP_EPS, same constant as
 * the uniaxial pair): elastic predictor -> radial-return direction n
 * fixed at the trial -> scalar Newton for the plastic multiplier with
 * backstress frozen (isotropic consistency) -> one explicit Eq. 10
 * backstress step. Approximate elastoplastic tangent (radial-return
 * structure) — ABAQUS just takes a few more equilibrium iterations.
 *
 * PROPS (nprops = 5 + 3N — note nu, unlike the uniaxial card):
 *   1: E  2: nu  3: sy0  4: Q  5: b  6..: C_k, r_k, m_k triplets
 * STATEV (nstatv = 7 + 6N):
 *   1: p, 2-7: plastic strain (engineering shears), 8..: a_k (6 each)
 *
 * Voigt order (ABAQUS): 11, 22, 33, 12, 13, 23; strains engineering shear,
 * stresses tensor shear; dotprod6 doubles shear terms.
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
#define OW_NEWTON_ITER 50
#define OW_TOL         1.0e-9
#define OW_SUBSTEP_EPS 1.0e-4
#define OW_MAX_SUBSTEP 1000

static double dotprod6(const double *A, const double *B)
{
    return A[0] * B[0] + A[1] * B[1] + A[2] * B[2]
         + 2.0 * (A[3] * B[3] + A[4] * B[4] + A[5] * B[5]);
}

/* Rotate a symmetric tensor 6-vector by the (Fortran column-major) DROT.
 * is_strain: shears stored engineering. Order: 11,22,33,12,13,23. */
static void rot6(const double *drot, double *v, int is_strain)
{
    const double half = is_strain ? 0.5 : 1.0;
    double T[3][3], RT[3][3], out[3][3];
    int i, j, k;
    T[0][0] = v[0]; T[1][1] = v[1]; T[2][2] = v[2];
    T[0][1] = T[1][0] = v[3] * half;
    T[0][2] = T[2][0] = v[4] * half;
    T[1][2] = T[2][1] = v[5] * half;
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
    const double E   = props[0];
    const double nu  = props[1];
    const double sy0 = props[2];
    const double Q   = props[3];
    const double b   = props[4];
    const int    nb  = (*nprops - 5) / 3;
    const int    nt  = *ntens;

    double C[OW_MAX_BACK], r[OW_MAX_BACK], m[OW_MAX_BACK];
    double ak[OW_MAX_BACK][6];
    int i, j, k;

    if (nt != 6 || nb < 1 || nb > OW_MAX_BACK || *nstatv < 7 + 6 * nb) {
        fprintf(stderr, "OW 3D UMAT: bad ntens=%d/nprops=%d/nstatv=%d\n",
                nt, *nprops, *nstatv);
        abort();
    }
    for (k = 0; k < nb; ++k) {
        C[k] = props[5 + 3 * k];
        r[k] = props[6 + 3 * k];
        m[k] = props[7 + 3 * k];
    }

    /* rotate state tensors by the increment rotation */
    double p = statev[0];
    double eps_p[6];
    for (i = 0; i < 6; ++i) eps_p[i] = statev[1 + i];
    rot6(drot, eps_p, 1);
    for (k = 0; k < nb; ++k) {
        for (i = 0; i < 6; ++i) ak[k][i] = statev[7 + 6 * k + i];
        rot6(drot, ak[k], 0);
    }

    /* elastic moduli */
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

    /* substep on the equivalent strain increment */
    double de_eq = SQRT23 * sqrt(dotprod6(dstran, dstran));
    int nsub = (int)ceil(de_eq / OW_SUBSTEP_EPS);
    if (nsub < 1) nsub = 1;
    if (nsub > OW_MAX_SUBSTEP) nsub = OW_MAX_SUBSTEP;
    double de[6];
    for (i = 0; i < 6; ++i) de[i] = dstran[i] / (double)nsub;

    double sig[6];
    for (i = 0; i < 6; ++i) sig[i] = stress[i];

    int any_plastic = 0;
    double H_last = 0.0, srn_last = 1.0, dl_last = 0.0;
    double n[6] = {0, 0, 0, 0, 0, 0};
    double adiff_last[6] = {0, 0, 0, 0, 0, 0};

    for (int isub = 0; isub < nsub; ++isub) {
        double sig_tr[6];
        for (i = 0; i < 6; ++i) {
            sig_tr[i] = sig[i];
            for (j = 0; j < 6; ++j) sig_tr[i] += cmat[i][j] * de[j];
        }
        double alpha[6] = {0, 0, 0, 0, 0, 0};
        for (k = 0; k < nb; ++k)
            for (i = 0; i < 6; ++i) alpha[i] += ak[k][i];
        const double hyd = (sig_tr[0] + sig_tr[1] + sig_tr[2]) / 3.0;
        double srel[6];
        for (i = 0; i < 6; ++i)
            srel[i] = sig_tr[i] - (i < 3 ? hyd : 0.0) - alpha[i];
        const double srn = sqrt(dotprod6(srel, srel));
        double sy = sy0 + Q * (1.0 - exp(-b * p));

        if (srn - SQRT23 * sy <= 0.0) {           /* elastic substep */
            for (i = 0; i < 6; ++i) sig[i] = sig_tr[i];
            continue;
        }

        for (i = 0; i < 6; ++i) n[i] = srel[i] / srn;

        /* scalar Newton for dl, backstress frozen:
         * F(dl) = srn - 2 mu dl - sqrt(2/3) sy(p + sqrt(2/3) dl) */
        double dl = (srn - SQRT23 * sy) / mu2;
        if (dl < 1.0e-14) dl = 1.0e-14;
        for (int it = 0; it < OW_NEWTON_ITER; ++it) {
            const double pp = p + SQRT23 * dl;
            const double ex = exp(-b * pp);
            const double F  = srn - mu2 * dl
                            - SQRT23 * (sy0 + Q * (1.0 - ex));
            if (fabs(F) < OW_TOL * (sy0 + srn)) break;
            const double dF = -mu2 - (2.0 / 3.0) * Q * b * ex;
            if (dF == 0.0) break;
            double dl_new = dl - F / dF;
            if (dl_new <= 0.0) dl_new = dl * 0.5;
            dl = dl_new;
        }

        p += SQRT23 * dl;
        any_plastic = 1;

        /* stress update */
        double dpe[6];
        for (i = 0; i < 6; ++i) dpe[i] = dl * n[i] * (i < 3 ? 1.0 : 2.0);
        for (i = 0; i < 6; ++i) {
            eps_p[i] += dpe[i];
            double corr = 0.0;
            for (j = 0; j < 6; ++j) corr += cmat[i][j] * dpe[j];
            sig[i] = sig_tr[i] - corr;
        }

        /* explicit OW-I Eq. 10 backstress update (de_p = dl * n, tensor
         * components) + hardening-modulus accumulation for the tangent */
        H_last = Q * b * exp(-b * p);
        double aold[6];
        for (i = 0; i < 6; ++i) aold[i] = 0.0;
        for (k = 0; k < nb; ++k)
            for (i = 0; i < 6; ++i) aold[i] += ak[k][i];
        for (k = 0; k < nb; ++k) {
            const double abar = sqrt(1.5 * dotprod6(ak[k], ak[k]));
            double recov = 0.0, proj = 0.0;
            double nhat[6] = {0, 0, 0, 0, 0, 0};
            if (abar > 1.0e-14 && r[k] > 0.0) {
                for (i = 0; i < 6; ++i) nhat[i] = ak[k][i] / abar;
                recov = pow(abar / r[k], m[k]);
                proj = dl * dotprod6(n, nhat);    /* de_p : nhat */
            }
            const double active = (proj > 0.0) ? 1.0 : 0.0;
            for (i = 0; i < 6; ++i)
                ak[k][i] += (2.0 / 3.0) * C[k] * dl * n[i]
                          - C[k] * recov * active * proj * nhat[i];
            /* uniaxial-equivalent modulus contribution */
            H_last += C[k] * (1.0 - recov * active);
        }
        for (i = 0; i < 6; ++i) {
            adiff_last[i] = -aold[i];
            for (k = 0; k < nb; ++k) adiff_last[i] += ak[k][i];
        }
        srn_last = srn;
        dl_last = dl;
    }

    for (i = 0; i < 6; ++i) stress[i] = sig[i];
    statev[0] = p;
    for (i = 0; i < 6; ++i) statev[1 + i] = eps_p[i];
    for (k = 0; k < nb; ++k)
        for (i = 0; i < 6; ++i) statev[7 + 6 * k + i] = ak[k][i];

    /* tangent: elastic, or radial-return structure with the last substep's
     * hardening modulus (approximate — stress is exact, ABAQUS's global
     * Newton just converges a little slower than with a consistent one) */
    if (!any_plastic) {
        for (i = 0; i < 6; ++i)
            for (j = 0; j < 6; ++j) ddsdde[j * 6 + i] = cmat[i][j];
    } else {
        double H = H_last;
        if (H < -2.0 * mu) H = -2.0 * mu;         /* keep denom sane */
        const double beta = 1.0 + H / (3.0 * mu);
        const double th1 = 1.0 - mu2 * dl_last / srn_last;
        const double th3 = 1.0 / (beta * srn_last);
        const double th2 = 1.0 / beta + dotprod6(n, adiff_last) * th3
                         - (1.0 - th1);
        double dd[6][6];
        for (i = 0; i < 6; ++i)
            for (j = 0; j < 6; ++j) {
                const double id2 = (i < 3 && j < 3) ? 1.0 : 0.0;
                const double id4 = (i == j) ? (i < 3 ? 1.0 : 0.5) : 0.0;
                dd[i][j] = bulk * id2
                         + mu2 * th1 * (id4 - (1.0 / 3.0) * id2)
                         - mu2 * th2 * n[i] * n[j]
                         + mu2 * th3 * adiff_last[i] * n[j];
            }
        for (i = 0; i < 6; ++i)
            for (j = 0; j < 6; ++j)
                ddsdde[j * 6 + i] = 0.5 * (dd[i][j] + dd[j][i]);
    }
}

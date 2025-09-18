#ifndef OME_H
#define OME_H
#include "constants.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>
#define FSQUARE(F) (F) * (F)
#define EPSILON (1e-5)
#define DIMIM (20)
#define DIMRE (40)
#define ZI (0.2)
#define FACPI (g_b * g_b / f_pi / f_pi / 24)
#include <complex.h>

// #define DEBUG
#define RECOIL
// #define PIIIII
// #define ETA
//
double complex ln1(double complex E, double complex p1, double complex p2,
                   double m0, double m1, double m2);
double complex ln2(double complex E, double complex p1, double complex p2,
                   double m0, double m1, double m2);

static inline double complex csquare(double complex x) { return x * x; }
static inline double fsquare(double x) { return x * x; }

struct OME {
    int set;
    double complex xxpiup[2 * DIMIM + DIMRE];
    double complex wwpiup[2 * DIMIM + DIMRE];
    double complex xxpidn[2 * DIMIM + DIMRE];
    double complex wwpidn[2 * DIMIM + DIMRE];
    double xxz[DIMRE];
    double wwz[DIMRE];
    double xx0ii[DIMRE];
    double ww0ii[DIMRE];
};

void ome_build(struct OME *self);
void ome_free(struct OME *self);
void printp(struct OME *self);

// Omega functions
static inline double complex omega_00(double complex p, double complex pprime) {
#ifdef RECOIL

    return 2 * m_B + (p * p + pprime * pprime) / (2 * m_B);
#else
    return 2 * m_B;
#endif // RECOIL
}

static inline double complex omega_01(double complex p, double complex pprime) {
#ifdef RECOIL
    return m_B + pprime * pprime / (2 * m_B) + m_B_s + p * p / (2 * m_B_s);
#else
    return m_B + m_B_s;
#endif
}

static inline double complex omega_10(double complex p, double complex pprime) {
#ifdef RECOIL
    return m_B_s + pprime * pprime / (2 * m_B_s) + m_B + p * p / (2 * m_B);
#else
    return m_B_s + m_B;
#endif
}

static inline double complex omega_11(double complex p, double complex pprime) {
#ifdef RECOIL
    return 2 * m_B_s + (p * p + pprime * pprime) / (2 * m_B_s);
#else
    return 2 * m_B_s;
#endif
}

static inline double complex omegaprime_00(double complex p,
                                           double complex pprime) {
#ifdef RECOIL
    return 2 * m_B_star + (p * p + pprime * pprime) / (2 * m_B_star);
#else
    return 2 * m_B_star;
#endif
}

static inline double complex omegaprime_01(double complex p,
                                           double complex pprime) {
#ifdef RECOIL
    return m_B_star + pprime * pprime / (2 * m_B_star) + m_B_star_s +
           p * p / (2 * m_B_star_s);
#else
    return m_B_star + m_B_star_s;
#endif
}

static inline double complex omegaprime_10(double complex p,
                                           double complex pprime) {
#ifdef RECOIL
    return m_B_star_s + pprime * pprime / (2 * m_B_star_s) + m_B_star +
           p * p / (2 * m_B_star);
#else
    return m_B_star_s + m_B_star;
#endif
}

static inline double complex omegaprime_11(double complex p,
                                           double complex pprime) {
#ifdef RECOIL
    return 2 * m_B_star_s + (p * p + pprime * pprime) / (2 * m_B_star_s);
#else
    return 2 * m_B_star_s;
#endif
}

static inline double complex Epi(double complex z, double complex p1,
                                 double complex p2, double m0) {
    auto A = p1 * p1 + p2 * p2 + m0 * m0;
    auto B = -2 * p1 * p2;
    return xsqrtleft(A + B * z);
    if (fabs(cimag(B)) < 1e-8) {
        return xsqrtright(cabs(A + B * z) + 0I);
    } else {
        return creal(A) >= 0 ? xsqrtleft(A + B * z) : xsqrtright(A + B * z);
    }
}

static inline double complex Dij(double complex E, double complex z,
                                 double complex p1, double complex p2,
                                 double complex mi, double complex mj,
                                 double m0) {
#ifdef RECOIL
    return E - (mi + p1 * p1 / (2 * mi)) - (mj + p2 * p2 / (2 * mj)) -
           Epi(z, p1, p2, m0) + I * EPSILON;
#else
    return E - mi - mj - Epi(z, p1, p2, m0) + I * EPSILON;
#endif
}

static inline double complex p1p2(double complex p1, double complex p2,
                                  double complex z) {
#ifdef DEBUG
    // return 1;
    return 1;
#else
    return p1 * p1 + p2 * p2 - 2 * p1 * p2 * z;
#endif
}

static inline double complex TOPTintegrand(double complex E, double complex z,
                                           double complex p1, double complex p2,
                                           double complex m1, double complex m2,
                                           double m0) {
    auto D1 = Dij(E, z, p1, p2, m1, m2, m0);
#ifdef DEBUG
    return (1 / D1) / (2 * Epi(z, p1, p2, m0)) * p1p2(p1, p2, z);
#else
    return FACPI * 1 / D1 / (2 * Epi(z, p1, p2, m0)) * p1p2(p1, p2, z);
#endif
}

static inline double complex z0(double complex E, double complex m1,
                                double complex m2, double complex p1,
                                double complex p2, double m0) {
#ifdef RECOIL
    auto res = (csquare(E - m1 - m2 - p1 * p1 / 2 / m1 - p2 * p2 / 2 / m2 +
                        I * EPSILON) -
                m0 * m0 - p1 * p1 - p2 * p2) /
               (-2 * p2 * p1);
#else
    auto res =
        (csquare(E - m1 - m2 + I * EPSILON) - m0 * m0 - p1 * p1 - p2 * p2) /
        (-2 * p2 * p1);
#endif
#ifdef DEBUG
    printf("z0: %s\n", formatC(res));
    printf("Dij: %e\n",
           cabs(xsqrt(p1 * p1 + p2 * p2 - 2 * p1 * p2 * res + m0 * m0) + m1 +
                m2 + p1 * p1 / 2 / m1 + p2 * p2 / 2 / m2 - E - I * EPSILON));
#endif
    if (cabs(Dij(E, res, p1, p2, m1, m2, m0)) < 1e-6) {
        return res;
    } else {
        return 5;
    }
}

static inline double complex z0E(double complex p1, double complex p2,
                                 double m0) {
    return (p1 * p1 + p2 * p2 + m0 * m0) / (2 * p2 * p1);
}

double complex Vpiu(struct OME ome, double complex E, double complex p1,
                    double complex p2, double complex m1, double complex m2,
                    double complex m3, double complex m4, double m0,
                    double fac);

double complex quad(double complex E, double complex p1, double complex p2);

static inline double complex quadreal(struct OME ome, double complex E,
                                      double complex p1, double complex p2,
                                      double complex m1, double complex m2,
                                      double m0, double fac) {
    double complex res = 0;
    for (size_t i = 0; i < DIMRE; i += 1) {
        auto z = ome.xxz[i];
        auto w = ome.wwz[i];
        // auto D1 = Dij(E, z, p1, p2, m1 - I * gam1 / 2, m3 - I * gam3 / 2,
        // m0); auto D2 = Dij(E, z, p1, p2, m2 - I * gam2 / 2, m4 - I * gam4 /
        // 2, m0); auto Dint = FACPI * fac * (1 / D1 + 1 / D2) / (2 * Epi(z, p1,
        // p2, m0)) * p1p2(p1, p2, z);
        auto Dint = fac * TOPTintegrand(E, z, p1, p2, m1, m2, m0);
        res += Dint * w;
    }
    return res;
}
static inline double complex quadup(struct OME ome, double complex E,
                                    double complex p1, double complex p2,
                                    double complex m1, double complex m2,
                                    double m0, double fac) {
    double complex res = 0;
    for (size_t i = 0; i < 2 * DIMIM + DIMRE; i += 1) {
        auto z = ome.xxpiup[i];
        auto w = ome.wwpiup[i];
        // auto D1 = Dij(E, z, p1, p2, m1 - I * gam1 / 2, m3 - I * gam3 / 2,
        // m0); auto D2 = Dij(E, z, p1, p2, m2 - I * gam2 / 2, m4 - I * gam4 /
        // 2, m0); auto Dint = FACPI * fac * (1 / D1 + 1 / D2) / (2 * Epi(z, p1,
        // p2, m0)) * p1p2(p1, p2, z);
        auto Dint = fac * TOPTintegrand(E, z, p1, p2, m1, m2, m0);
        res += Dint * w;
    }
    return res;
}
static inline double complex quaddn(struct OME ome, double complex E,
                                    double complex p1, double complex p2,
                                    double complex m1, double complex m2,
                                    double m0, double fac) {
    double complex res = 0;
    for (size_t i = 0; i < 2 * DIMIM + DIMRE; i += 1) {
        auto z = ome.xxpidn[i];
        auto w = ome.wwpidn[i];
        // auto D1 = Dij(E, z, p1, p2, m1 - I * gam1 / 2, m3 - I * gam3 / 2,
        // m0); auto D2 = Dij(E, z, p1, p2, m2 - I * gam2 / 2, m4 - I * gam4 /
        // 2, m0); auto Dint = FACPI * fac * (1 / D1 + 1 / D2) / (2 * Epi(z, p1,
        // p2, m0)) * p1p2(p1, p2, z);
        auto Dint = fac * TOPTintegrand(E, z, p1, p2, m1, m2, m0);
        res += Dint * w;
    }
    return res;
}

static inline double complex quadii(struct OME ome, double complex E,
                                    double complex p1, double complex p2,
                                    double complex m1, double complex m2,
                                    double m0, double fac, double complex _z,
                                    double complex offset) {
    double complex res = 0;
    for (size_t i = 0; i < DIMRE; i += 1) {
        auto z = ome.xx0ii[i] * _z + offset;
        auto w = ome.ww0ii[i] * _z;
        // auto D1 = Dij(E, z, p1, p2, m1 - I * gam1 / 2, m3 - I * gam3 / 2,
        // m0); auto D2 = Dij(E, z, p1, p2, m2 - I * gam2 / 2, m4 - I * gam4 /
        // 2, m0); auto Dint = FACPI * fac * (1 / D1 + 1 / D2) / (2 * Epi(z, p1,
        // p2, m0)) * p1p2(p1, p2, z);
        auto Dint = fac * TOPTintegrand(E, z, p1, p2, m1, m2, m0);
        res += Dint * w;
    }
    return res;
}

static inline double complex OME_00(struct OME ome, double complex E,
                                    double complex p, double complex pprime) {
#ifdef PIIIIII
    return Vpiu(ome, E, p, pprime, m_B_star, gamma_B_star, m_B, 0, m_B_star,
                gamma_B_star, m_B, 0, m_pi, 1);
#else
    return Vpiu(ome, E, p, pprime, m_B_star, m_B, m_B_star, m_B, m_pi, 3) +
           Vpiu(ome, E, p, pprime, m_B_star, m_B, m_B_star, m_B, m_eta, 1. / 3);
#endif
}

static inline double complex OME_01(struct OME ome, double complex E,
                                    double complex p, double complex pprime) {
#ifdef PIIIIII
    return Vpiu(ome, E, p, pprime, m_B_star_s, gamma_B_star_s, m_B_s, 0,
                m_B_star, gamma_B_star, m_B, 0, m_pi, 1);
#else
    return Vpiu(ome, E, p, pprime, m_B_star_s, m_B_s, m_B_star, m_B, m_K,
                pow(2, 3. / 2));
#endif
}

static inline double complex OME_10(struct OME ome, double complex E,
                                    double complex p, double complex pprime) {
#ifdef PIIIIII

    return Vpiu(ome, E, p, pprime, m_B_star, gamma_B_star, m_B, 0, m_B_star_s,
                gamma_B_star_s, m_B_s, 0, m_pi, 1);
#else
    return Vpiu(ome, E, p, pprime, m_B_star, m_B, m_B_star_s, m_B_s, m_K,
                pow(2, 3. / 2));
#endif
}

static inline double complex OME_11(struct OME ome, double complex E,
                                    double complex p, double complex pprime) {
#ifdef PIIIIII

    return Vpiu(ome, E, p, pprime, m_B_star_s, gamma_B_star_s, m_B_s, 0,
                m_B_star_s, gamma_B_star_s, m_B_s, 0, m_pi, 1);
#else
    return Vpiu(ome, E, p, pprime, m_B_star_s, m_B_s, m_B_star_s, m_B_s, m_eta,
                4. / 3);
#endif
}

#define DEFINE_DELTA0(suffix)                                                  \
    static inline double complex Delta0_##suffix(                              \
        double complex e, double complex p, double complex pprime,             \
        double m0) {                                                           \
        auto A = p * p + pprime * pprime + m0 * m0;                            \
        auto B = 2 * p * pprime;                                               \
        auto C = omega_##suffix(pprime, p) - e - I * EPSILON;                  \
        auto D = omegaprime_##suffix(p, pprime) - e - I * EPSILON;             \
        auto a = Epi(1, p, pprime, m0);                                        \
        auto b = Epi(-1, p, pprime, m0);                                       \
        auto log1 = (xlog((a + C)) - xlog((b + C))) / B;                       \
        auto log2 = (xlog((a + D)) - xlog((b + D))) / B;                       \
        return -(log1 + log2);                                                 \
    }

#define DEFINE_DELTA1(suffix)                                                  \
    static inline double complex Delta1_##suffix(                              \
        double complex E, double complex p, double complex pprime,             \
        double m0) {                                                           \
        auto A = p * p + pprime * pprime + m0 * m0;                            \
        auto B = 2 * p * pprime;                                               \
        auto C = omega_##suffix(pprime, p) - E - I * EPSILON;                  \
        auto D = omegaprime_##suffix(p, pprime) - E - I * EPSILON;             \
        auto a = Epi(1, p, pprime, m0);                                        \
        auto b = Epi(-1, p, pprime, m0);                                       \
        auto log1 = (xlog((a + C)) - xlog((b + C))) / B;                       \
        auto log2 = (xlog((a + D)) - xlog((b + D))) / B;                       \
        auto B0 = B;                                                           \
        auto ret = -((C + D) * (a - b) / B0 / B + 2 / B +                      \
                     (A - C * C) * (log1) / B + (A - D * D) * (log2) / B);     \
        return ret;                                                            \
    }
#define DEFINE_ANA(suffix)                                                     \
    static inline double complex ANA_##suffix(                                 \
        double complex E, double complex p, double complex pprime,             \
        double m0) {                                                           \
        return FACPI * (2 * p * pprime * Delta1_##suffix(E, p, pprime, m0) -   \
                        (p * p + pprime * pprime) *                            \
                            Delta0_##suffix(E, p, pprime, m0));                \
    }

double complex juliana(double complex E, double complex p,
                       double complex pprime);

DEFINE_DELTA0(00);
DEFINE_DELTA0(01);
DEFINE_DELTA0(10);
DEFINE_DELTA0(11);

DEFINE_DELTA1(00);
DEFINE_DELTA1(01);
DEFINE_DELTA1(10);
DEFINE_DELTA1(11);

DEFINE_ANA(00);
DEFINE_ANA(01);
DEFINE_ANA(10);
DEFINE_ANA(11);

static inline double complex OMEANA_00(double complex E, double complex p,
                                       double complex pprime) {
#ifdef PIIIII
    return ANA_00(E, p, pprime, m_pi);
#elif defined(ETA)
    return ANA_00(E, p, pprime, m_eta);
#else
    return 3 * ANA_00(E, p, pprime, m_pi) +
           1. / 3 * ANA_00(E, p, pprime, m_eta);
#endif
}

static inline double complex OMEANA_01(double complex E, double complex p,
                                       double complex pprime) {
#ifdef PIIIII
    return ANA_01(E, p, pprime, m_pi);
#elif defined(ETA)
    return ANA_00(E, p, pprime, m_eta);
#else
    return 2 * sqrt(2) * ANA_01(E, p, pprime, m_K);
#endif
}

static inline double complex OMEANA_10(double complex E, double complex p,
                                       double complex pprime) {
#ifdef PIIIII
    return ANA_10(E, p, pprime, m_pi);
#elif defined(ETA)
    return ANA_00(E, p, pprime, m_eta);
#else
    return 2 * sqrt(2) * ANA_10(E, p, pprime, m_K);
#endif
}

static inline double complex OMEANA_11(double complex E, double complex p,
                                       double complex pprime) {
#ifdef PIIIII
    return ANA_11(E, p, pprime, m_pi);
#elif defined(ETA)
    return ANA_00(E, p, pprime, m_eta);
#else
    return 4. / 3 * ANA_11(E, p, pprime, m_eta);
#endif
}

#endif // OME_H

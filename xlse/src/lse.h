#ifndef LSE_H
#define LSE_H
#include "constants.h"
#include "ome.h"
#include "utils.h"
#include "wavefunction.h"
#include <complex.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>

#define NCHANNELS 3

typedef gsl_matrix_complex matrix;
#define inverse gsl_complex_inverse
#define add gsl_complex_add
#define mul gsl_complex_mul
#define mul_real gsl_complex_mul_real
#define sub gsl_complex_sub
#define matrix_alloc gsl_matrix_complex_alloc
#define matrix_free gsl_matrix_complex_free
#define matrix_set_zero gsl_matrix_complex_set_zero
#define matrix_set gsl_matrix_complex_set
#define matrix_memcpy gsl_matrix_complex_memcpy
#define matrix_scale gsl_matrix_complex_scale
#define matrix_get gsl_matrix_complex_get

// #define TPO

#define vlarray(ptr, ...) ((double complex(*) __VA_ARGS__)(ptr))

typedef struct {
    size_t pNgauss;
    double Lambda;
    double epsilon;
    double complex E;
    double complex det;
    matrix *TOME;
    matrix *VOME;
    matrix *G;
    matrix *reg;
    double complex x0[NCHANNELS];
    gsl_integration_glfixed_table *table;
    WaveFunction *wf;
    double *xi;
    double *wi;
    struct OME ome;
    void *psi_n_mat;
    double *E_vec;
    void *Xin;
    void *Xout;
    double complex onshellT[NCHANNELS][NCHANNELS];
    void *v;
    void *sigmat;
    double V0;
    double C00;
    double C01;
    double C10;
    double C11;
} LSE;

typedef enum {
    NN = 0,
    PN = 1,
    NP = 2,
    PP = 3,
} RS;

LSE *lse_malloc(size_t pNgauss, double Lambda, double epsilon);
int lse_compute(LSE *self, double complex E, double C[NCHANNELS * NCHANNELS],
                RS rs);
int lse_compute_single(LSE *self, double complex E,
                       double C[NCHANNELS * NCHANNELS], RS rs);
void lse_free(LSE *self);
double complex *lse_get_g_data(LSE *self);
void lse_get_g_size(LSE *self, unsigned int *rows, unsigned int *cols);
double complex *lse_get_v_data(LSE *self);
void lse_get_v_size(LSE *self, unsigned int *rows, unsigned int *cols);
double complex *lse_get_t_data(LSE *self);
void lse_get_t_size(LSE *self, unsigned int *rows, unsigned int *cols);
double complex *lse_get_ivg_data(LSE *self);
void lse_get_ivg_size(LSE *self, unsigned int *rows, unsigned int *cols);
double complex *lse_get_iivg_data(LSE *self);
void lse_get_iivg_size(LSE *self, unsigned int *rows, unsigned int *cols);
double complex *lse_get_psi(LSE *self);
void lse_get_psi_size(LSE *self, unsigned int *rows, unsigned int *cols);
double *lse_get_E(LSE *self);
void lse_get_E_size(unsigned int *levels);
void lse_get_M_size(LSE *self, unsigned int *rows, unsigned int *cols);
double complex *lse_get_onshellT(LSE *self);
double complex pole(LSE *lse, double complex E, double C[NCHANNELS * NCHANNELS],
                    RS rs);
double cost(const gsl_vector *x, void *params);
double *minimize(LSE *lse, double C[NCHANNELS * NCHANNELS]);

// LSE methods
int lse_gmat(LSE *self);
int lse_vmat(LSE *self);
int lse_tmat(LSE *self);
int lse_tmat_single(LSE *self);
double complex lse_invT(LSE *self, double complex E,
                        double C[NCHANNELS * NCHANNELS], RS rs);
double complex lse_detImVG(LSE *self, double complex E,
                           double C[NCHANNELS * NCHANNELS], RS rs);
double complex lse_detVG(LSE *self, double complex E,
                         double C[NCHANNELS * NCHANNELS], RS rs);
double lse_cost(LSE *self, double C[NCHANNELS * NCHANNELS], RS rs);
void lse_refresh(LSE *self, double complex E, double C[NCHANNELS * NCHANNELS],
                 RS rs);
void lse_X(LSE *self);
void lse_XtX(LSE *self);

static inline double min(double x, double y) { return x > y ? y : x; }

static inline double max(double x, double y) { return x > y ? x : y; }

double complex V_QM_00(LSE *self, size_t p, size_t pprime);
double complex V_QM_01(LSE *self, size_t p, size_t pprime);
double complex V_QM_10(LSE *self, size_t p, size_t pprime);
double complex V_QM_11(LSE *self, size_t p, size_t pprime);

#define DEFINE_VQMTEST(alpha, beta)                                            \
    static inline double complex V_QM_TEST_##alpha##beta(LSE *self, size_t p,  \
                                                         size_t pprime) {      \
        double E[6] = {-0.2, -0.8, -1.2, -0.1, -0.5, 0.4};                     \
        double complex res = 0;                                                \
        for (size_t i = 0; i < 6; i += 1) {                                    \
            res -= 1 / (self->E - E[i]);                                       \
        }                                                                      \
        return res * g##alpha * g##beta;                                       \
    }

DEFINE_VQMTEST(0, 0)
DEFINE_VQMTEST(0, 1)
DEFINE_VQMTEST(1, 0)
DEFINE_VQMTEST(1, 1)

#define DEFINE_V_TEST(alpha, beta)                                             \
    static inline double complex V_TEST_##alpha##beta(                         \
        double complex E, double complex p, double complex pprime) {           \
        return csin(E) * (p + pprime * pprime - 1 / p) * g##alpha * g##beta;   \
    }

DEFINE_V_TEST(0, 0)
DEFINE_V_TEST(0, 1)
DEFINE_V_TEST(1, 0)
DEFINE_V_TEST(1, 1)

#define DEFINE_VQM(alpha, beta)                                                \
    double complex V_QM_##alpha##beta(LSE *self, size_t p, size_t pprime) {    \
        double complex res = 0 + 0I;                                           \
        auto E = self->E;                                                      \
        size_t pNgauss = self->pNgauss;                                        \
        auto psi =                                                             \
            (double complex(*)[N_MAX + 1][pNgauss + 1]) self->psi_n_mat;       \
        size_t chan0 = p / (pNgauss + 1);                                      \
        size_t chan1 = pprime / (pNgauss + 1);                                 \
        for (size_t i = 0; i < N_TOWER; i++) {                                 \
            res += psi[chan0][i][p % (pNgauss + 1)] *                          \
                   conj(psi[chan1][i][pprime % (pNgauss + 1)]) /               \
                   (E - self->E_vec[i] - self->V0);                            \
        }                                                                      \
        return res * g##alpha * g##beta;                                       \
    }

#define DEFINE_V_FUNCTION(suffix)                                              \
    static inline gsl_complex V##suffix(LSE *self, double complex p,           \
                                        double complex pprime, size_t pi,      \
                                        size_t ppi) {                          \
        auto E = self->E;                                                      \
        E += m11 + m12;                                                        \
        auto res = OMEANA_##suffix(E, p, pprime);                              \
        return res;                                                            \
    }

DEFINE_V_FUNCTION(00);
DEFINE_V_FUNCTION(01);
DEFINE_V_FUNCTION(02);
DEFINE_V_FUNCTION(10);
DEFINE_V_FUNCTION(11);
DEFINE_V_FUNCTION(12);
DEFINE_V_FUNCTION(20);
DEFINE_V_FUNCTION(21);
DEFINE_V_FUNCTION(22);

#define DEFINE_EFFMOM(alpha, beta)                                             \
    static inline double complex effmom##alpha##beta(double complex E,         \
                                                     double complex p) {       \
        double m1 = m_B##alpha;                                                \
        double m2 = m_B##beta;                                                 \
        double complex m2star = m_B##beta##_star;                              \
        double m0 = m_##alpha##beta;                                           \
        auto muBphi = m2 * m0 / (m2 + m0);                                     \
        auto muBB = m1 * m2star / (m1 + m2star);                               \
        return xsqrtright(2 * muBphi * (E - m1 - m2 - m0 - p * p / 2 / muBB)); \
    }

DEFINE_EFFMOM(0, 1);
DEFINE_EFFMOM(0, 2);
DEFINE_EFFMOM(1, 0);
DEFINE_EFFMOM(1, 2);
DEFINE_EFFMOM(2, 0);
DEFINE_EFFMOM(2, 1);
DEFINE_EFFMOM(2, 2);

static inline double complex effmom(double complex E, double complex p,
                                    double m1, double m2, double complex m2star,
                                    double m0) {
    auto muBphi = m2 * m0 / (m2 + m0);
    auto muBB = m1 * m2star / (m1 + m2star);
    return xsqrtright(2 * muBphi * (E - m1 - m2 - m0 - p * p / 2 / muBB));
}

static inline double complex effmom00pi(double complex E, double complex p) {
    return effmom(E, p, m_B0, m_B0, m_B0_star, m_pi);
}

static inline double complex effmom00eta(double complex E, double complex p) {
    return effmom(E, p, m_B0, m_B0, m_B0_star, m_eta);
}

static inline double complex effmom11pi(double complex E, double complex p) {
    return effmom(E, p, m_B1, m_B1, m_B1_star, m_pi);
}

static inline double complex effmom11eta(double complex E, double complex p) {
    return effmom(E, p, m_B1, m_B1, m_B1_star, m_eta);
}

#define DEFINE_WIDTH(alpha, beta)                                              \
    static inline double complex Gamma##alpha##beta(double complex E,          \
                                                    double complex p) {        \
        auto fac = g_b * g_b * m_B##beta * 2 / 24 / M_PI / f_pi / f_pi /       \
                   m_B##alpha##_star;                                          \
        return fac * (xcube(effmom##alpha##beta(E, p)) -                       \
                      creal(xcube(effmom##alpha##beta(                         \
                          m_B##alpha + m_B##alpha##_star, 0))));               \
    }

DEFINE_WIDTH(0, 1);
DEFINE_WIDTH(0, 2);
DEFINE_WIDTH(1, 0);
DEFINE_WIDTH(1, 2);
DEFINE_WIDTH(2, 0);
DEFINE_WIDTH(2, 1);

static inline double complex Gamma00pi(double complex E, double complex p) {
    auto fac = g_b * g_b * m_B0 / 24 / M_PI / f_pi / f_pi / m_B0_star;
    return fac * (xcube(effmom00pi(E, p)) -
                  creal(xcube(effmom00pi(m_B0 + m_B0_star, 0))));
}

static inline double complex Gamma00eta(double complex E, double complex p) {
    auto fac = g_b * g_b * m_B0 / 3 / 24 / M_PI / f_pi / f_pi / m_B0_star;
    return fac * (xcube(effmom00eta(E, p)) -
                  creal(xcube(effmom00eta(m_B0 + m_B0_star, 0))));
}

static inline double complex Gamma11pi(double complex E, double complex p) {
    auto fac = g_b * g_b * m_B1 / 24 / M_PI / f_pi / f_pi / m_B1_star;
    return fac * (xcube(effmom11pi(E, p)) -
                  creal(xcube(effmom11pi(m_B1 + m_B1_star, 0))));
}

static inline double complex Gamma11eta(double complex E, double complex p) {
    auto fac = g_b * g_b * m_B1 / 24 / M_PI / f_pi / f_pi / m_B1_star;
    return fac * (xcube(effmom11eta(E, p)) -
                  creal(xcube(effmom11eta(m_B1 + m_B1_star, 0))));
}

static inline double complex Gamma22(double complex E, double complex p) {
    auto fac = g_b * g_b * m_B2 * 4 / 3 / 24 / M_PI / f_pi / f_pi / m_B2_star;
    return fac * (xcube(effmom22(E, p)) -
                  creal(xcube(effmom22(m_B2 + m_B2_star, 0))));
}

static inline double complex Gamma0(double complex E, double complex p) {
    return Gamma00pi(E, p) + Gamma00eta(E, p) + Gamma01(E, p) + Gamma02(E, p);
}

static inline double complex Gamma1(double complex E, double complex p) {
    return Gamma10(E, p) + Gamma11pi(E, p) + Gamma11eta(E, p) + Gamma12(E, p);
}

static inline double complex Gamma2(double complex E, double complex p) {
    return Gamma20(E, p) + Gamma21(E, p) + Gamma22(E, p);
}

#endif // !LSE_H

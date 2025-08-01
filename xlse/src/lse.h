#ifndef LSE_H
#define LSE_H
#include "constants.h"
#include "ome.h"
#include "wavefunction.h"
#include <complex.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>

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
      double complex x0[2];
      gsl_integration_glfixed_table *table;
      WaveFunction *wf;
      double *xi;
      double *wi;
      struct OME ome;
      void *psi_n_mat;
      double *E_vec;
      void *Xin;
      void *Xout;
      double complex onshellT[2][2];
      void *v;
      void *sigmat;
      double V0;
      double C00;
      double C01;
      double C10;
      double C11;
      double deviation;
} LSE;

typedef enum {
      NN = 0,
      PN = 1,
      NP = 2,
      PP = 3,
} RS;

LSE *lse_malloc(size_t pNgauss, double Lambda, double epsilon, double devi);
int lse_compute(LSE *self, double complex E, double C[4], RS rs, double devi);
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
double complex pole(LSE *lse, double complex E, double C[4], RS rs);
double cost(const gsl_vector *x, void *params);
double *minimize(LSE *lse, double C[4]);

// LSE methods
int lse_gmat(LSE *self);
int lse_vmat(LSE *self);
int lse_tmat(LSE *self);
double complex lse_detImVG(LSE *self, double complex E, double C[4], RS rs, double devi);
double complex lse_detVG(LSE *self, double complex E, double C[4], RS rs, double devi);
double lse_cost(LSE *self, double C[4], RS rs);
void lse_refresh(LSE *self, double complex E, double C[4], RS rs, double devi);
void lse_X(LSE *self);
void lse_XtX(LSE *self);

static inline double min(double x, double y) { return x > y ? y : x; }

static inline double max(double x, double y) { return x > y ? x : y; }

// Complex square root function
static inline double complex xsqrt(double complex x)
{
      if (cimag(x) >= 0)
	    return csqrt(x);
      else
	    return -csqrt(x - 0 * I);
}

// Contact term functions
static inline double Ctct_00(double g_C) { return -3 * 2 * g_C; }

static inline double Ctct_01(double g_C) { return pow(1.0 / 2.0, 3.0 / 2.0) * 4 * g_C; }

static inline double Ctct_10(double g_C) { return pow(1.0 / 2.0, 3.0 / 2.0) * 4 * g_C; }

static inline double Ctct_11(double g_C) { return g_C; }

#define DEFINE_O_FUNCTION(suffix)                                                                                              \
      static inline double complex O_##suffix(double complex E, double complex p, double complex pprime, double m)             \
      {                                                                                                                        \
	    if (cabs(p) <= 1e-8) {                                                                                             \
		  p += 1e-6;                                                                                                   \
	    }                                                                                                                  \
	    if (cabs(pprime) <= 1e-8) {                                                                                        \
		  pprime += 1e-6;                                                                                              \
	    }                                                                                                                  \
	    return 4 * fsquare(g_b) / fsquare(f_pi) * -1. / 4. / p / pprime *                                                  \
		   (clog((E - (m + csquare(p - pprime) / 2 / m) - omega_##suffix(p, pprime)) /                                 \
			 (E - (m + csquare(p + pprime) / 2 / m) - omega_##suffix(p, pprime))) +                                \
		    clog((E - (m + csquare(p - pprime) / 2 / m) - omegaprime_##suffix(p, pprime)) /                            \
			 (E - (m + csquare(p + pprime) / 2 / m) - omegaprime_##suffix(p, pprime))));                           \
      }

// DEFINE_O_FUNCTION(00);
// DEFINE_O_FUNCTION(01);
// DEFINE_O_FUNCTION(10);
// DEFINE_O_FUNCTION(11);

// static inline double complex V_OME_00(double complex E, double complex p, double complex pprime)
// {
//       return -3 * (3 * O_00(E, p, pprime, m_pi) + O_00(E, p, pprime, m_eta) / 3);
// }
//
// static inline double complex V_OME_01(double complex E, double complex p, double complex pprime)
// {
//       return pow(2, 3. / 2) * O_01(E, p, pprime, m_K);
// }
//
// static inline double complex V_OME_10(double complex E, double complex p, double complex pprime)
// {
//       return pow(2, 3. / 2) * O_10(E, p, pprime, m_K);
// }
//
// static inline double complex V_OME_11(double complex E, double complex p, double complex pprime)
// {
//       return 2. / 3 * O_11(E, p, pprime, m_eta);
// }
//
double complex V_QM_00(LSE *self, size_t p, size_t pprime);
double complex V_QM_01(LSE *self, size_t p, size_t pprime);
double complex V_QM_10(LSE *self, size_t p, size_t pprime);
double complex V_QM_11(LSE *self, size_t p, size_t pprime);

#define DEFINE_VQMTEST(alpha, beta)                                                                                            \
      static inline double complex V_QM_TEST_##alpha##beta(LSE *self, size_t p, size_t pprime)                                 \
      {                                                                                                                        \
	    double E[6] = { -0.2, -0.8, -1.2, -0.1, -0.5, 0.4 };                                                               \
	    double complex res = 0;                                                                                            \
	    for (size_t i = 0; i < 6; i += 1) {                                                                                \
		  res -= 1 / (self->E - E[i]);                                                                                 \
	    }                                                                                                                  \
	    return res * g##alpha * g##beta;                                                                                   \
      }

DEFINE_VQMTEST(0, 0)
DEFINE_VQMTEST(0, 1)
DEFINE_VQMTEST(1, 0)
DEFINE_VQMTEST(1, 1)

#define DEFINE_V_TEST(alpha, beta)                                                                                             \
      static inline double complex V_TEST_##alpha##beta(double complex E, double complex p, double complex pprime)             \
      {                                                                                                                        \
	    return csin(E) * (p + pprime * pprime - 1 / p) * g##alpha * g##beta;                                               \
      }

DEFINE_V_TEST(0, 0)
DEFINE_V_TEST(0, 1)
DEFINE_V_TEST(1, 0)
DEFINE_V_TEST(1, 1)

#define DEFINE_VQM(alpha, beta)                                                                                                \
      double complex V_QM_##alpha##beta(LSE *self, size_t p, size_t pprime)                                                    \
      {                                                                                                                        \
	    double complex res = 0 + 0I;                                                                                       \
	    auto E = self->E;                                                                                                  \
	    size_t pNgauss = self->pNgauss;                                                                                    \
	    auto psi = (double complex(*)[N_MAX + 1][pNgauss + 1]) self->psi_n_mat;                                            \
	    size_t chan0 = p / (pNgauss + 1);                                                                                  \
	    size_t chan1 = pprime / (pNgauss + 1);                                                                             \
	    for (size_t i = 0; i < N_TOWER; i++) {                                                                             \
		  res += psi[chan0][i][p % (pNgauss + 1)] * conj(psi[chan1][i][pprime % (pNgauss + 1)]) /                      \
			 (E - self->E_vec[i] - self->V0);                                                                      \
	    }                                                                                                                  \
	    return res * g##alpha * g##beta;                                                                                   \
      }

static inline double complex curlO(double complex p, double complex pprime, double m)
{
      if (cabs(p) < 1e-8)
	    p += 1e-6;
      if (cabs(pprime) < 1e-8)
	    pprime += 1e-6;
      return 4 * fsquare(g_b) / fsquare(f_pi) *
	     (1 - fsquare(m) / (4 * p * pprime) *
		      (clog((csquare(p) + csquare(pprime) + 2 * p * pprime + fsquare(m))) -
		       clog((csquare(p) + csquare(pprime) - 2 * p * pprime + fsquare(m)))));
}

static inline double complex V_curlOME_00(double complex E, double complex p, double complex pprime)
{
      return -3 * (3 * curlO(p, pprime, m_pi) + curlO(p, pprime, m_eta) / 3);
}

static inline double complex V_curlOME_01(double complex E, double complex p, double complex pprime)
{
      return pow(2, 3. / 2) * curlO(p, pprime, m_K);
}

static inline double complex V_curlOME_10(double complex E, double complex p, double complex pprime)
{
      return pow(2, 3. / 2) * curlO(p, pprime, m_K);
}

static inline double complex V_curlOME_11(double complex E, double complex p, double complex pprime)
{
      return 2. / 3 * curlO(p, pprime, m_eta);
}

#define DEFINE_V_FUNCTION(suffix)                                                                                              \
      static inline gsl_complex V##suffix(LSE *self, double complex p, double complex pprime, size_t pi, size_t ppi)           \
      {                                                                                                                        \
	    auto E = self->E;                                                                                                  \
	    E += m11 + m12;                                                                                                    \
	    return OMEANA_##suffix(E, p, pprime, self->deviation);                                                             \
      }

DEFINE_V_FUNCTION(00);
DEFINE_V_FUNCTION(01);
DEFINE_V_FUNCTION(10);
DEFINE_V_FUNCTION(11);

#endif // !LSE_H

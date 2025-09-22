#ifndef CONSTANTS_H
#define CONSTANTS_H
#include <complex.h>
#define PNGAUSS 64
#define RNGAUSS 256
#define RLAMBDA 20

#define f_pi (0.092)
#define g_b (0.5704)
#define g_c (-1. / 4.)

#define m_pi (0.138039407)
#define m_K (0.498)
#define m_eta (0.548)
#define m_eta_s (0.69)

#define WIDTH (-1e-6 * I)
// #define DEVI (-0.2)
#define m_B (5.27934)
#define m_B_star (5.32471)
#define m_B_s (5.36691)
#define m_B_star_s (5.4154 + WIDTH)
#define gamma_B_star (1e-6)
#define gamma_B_star_s (1e-6)
#define m_B0 (5.27941)
#define m_B1 (5.27963)
#define m_B2 (5.36691)
#define m_B0_star (5.32475 + WIDTH)
#define m_B1_star (5.32475 + WIDTH)
#define m_B2_star (5.4154 + WIDTH)

#define m_01 (m_pi)
#define m_02 (m_K)
#define m_10 (m_pi)
#define m_12 (m_K)
#define m_20 (m_K)
#define m_21 (m_K)
#define m_22 (m_eta)

#define m11 (5.27934)
#define m12 (5.32471)

#define m21 (5.36692)
#define m22 (5.4154)

#define m_Xb13P (10.5134 - (m11 + m12))
#define m_Xb12P (10.25546 - (m11 + m12))
#define m_Xb11P (9.89278 - (m11 + m12))

#define a_l (0.06426 * 5.068)
#define gcoupling (0.898794378386677 / a_l * 10)
#define g0 (gcoupling * 0.014926616931653945)
#define g1 (gcoupling * 0.006467550544943349)
// #define g1 (36. / 1000 * 1.4142135623730951)
#define delta0 0
#define delta1 (m_B1 + m_B1_star - m_B0 - m_B0_star)
#define delta2 (m_B2 + m_B2_star - m_B0 - m_B0_star)
#define mu0 (m_B0 * m_B0_star / (m_B0 + m_B0_star))
#define mu1 (m_B1 * m_B1_star / (m_B1 + m_B1_star))
#define mu2 (m_B2 * m_B2_star / (m_B2 + m_B2_star))
#define m_c 1.85
#define a_Cornell 1.95
#define g_qm (1)
#define partialwave 1

// wavefunction constants
#define N_MAX 64
#define N_TOWER 3
#define R_1 (0.02 * 5.068)      // GeV^-1
#define R_N_MAX (N_MAX * 5.068) // GeV^-1
#define C_T ((-1) * 0.19732 * 0.19732 / (0.5 * 4.18) * 5.068)
#define PI 3.14159265358979323846
#define V0FIT 0.4
#define SIGMA_L 0.0199179973550142
#define SIGMA (SIGMA_L / a_l / a_l)
#define ALPHA 0.476814032326273

constexpr double complex mB[3] = {m_B0, m_B1, m_B2};
constexpr double complex mBstar[3] = {m_B0_star, m_B1_star, m_B2_star};
constexpr double complex delta[3] = {delta0, delta1, delta2};
constexpr double complex mu[3] = {mu0, mu1, mu2};
constexpr double Delta[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
constexpr double gab[3] = {g0, g1};
constexpr double mdiff[2] = {m_Xb12P - m_Xb11P, m_Xb13P - m_Xb11P};

#endif // !DEBUG

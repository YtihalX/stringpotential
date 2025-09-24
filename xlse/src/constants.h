#ifndef CONSTANTS_H
#define CONSTANTS_H
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
#define m_B_star (5.32471 + WIDTH)
#define m_B_s (5.36692)
#define m_B_star_s (5.4154 + WIDTH)
#define gamma_B_star (1e-6)
#define gamma_B_star_s (1e-6)

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
#define delta1 (m21 + m22 - m11 - m12)
#define mu0 (m11 * m12 / (m11 + m12))
#define mu1 (m21 * m22 / (m21 + m22))
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

extern const double delta[2];
extern const double mu[2];
extern const double Delta[2][2];
extern const double gab[2];
extern const double mdiff[2];

#endif // !DEBUG

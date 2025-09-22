#include "autofree.h"
#include "constants.h"
#include "lse.h"
#include "ome.h"
#include "script.h"
#include "utils.h"
#include "wavefunction.h"
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix_double.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

double complex delta0lim(double complex E, double complex p, double m0,
                         double m1, double m2) {
    auto sq = csqrt(p * p + m0 * m0);
    // auto omega = m1 + m2 + p * p / 2 / m2 - E;
    auto omega = 2 * m_B + p * p / 2 / m_B - E - m11 - m12;
    auto ret = 2 / sq / (sq + omega);
    return ret;
}

double complex delta01(double complex E, double complex p1, double complex p2,
                       double m0, double complex m1, double complex m2,
                       double complex m3, double complex m4) {
    printf("E: %s\n", formatC(E));
    printf("p1: %s\n", formatC(p1));
    printf("p2: %s\n", formatC(p2));
    auto B0 = p1 * p1 + p2 * p2 + m0 * m0;
    printf("B0: %s\n", formatC(B0));
    auto C0 = 2 * p1 * p2;
    printf("C0: %s\n", formatC(C0));
    auto D0 =
        -(m1 + m3 + p1 * p1 / 2 / m1 + p2 * p2 / 2 / m3 - E - I * EPSILON);
    printf("D0: %s\n", formatC(D0));
    auto E0 =
        -(m2 + m4 + p1 * p1 / 2 / m2 + p2 * p2 / 2 / m4 - E - I * EPSILON);
    printf("E0: %s\n", formatC(E0));
    auto a0 = xsqrtright(B0 + C0);
    printf("a0: %s\n", formatC(a0));
    auto b0 = xsqrtright(B0 - C0);
    printf("b0: %s\n", formatC(b0));
    auto xlog1 = xlog((a0 - D0) / (b0 - D0));
    auto xlog2 = xlog((a0 - E0) / (b0 - E0));
    printf("log1: %s\n", formatC(xlog1));
    printf("log2: %s\n", formatC(xlog2));
    auto ret = (xlog1 + xlog2) / C0;
    printf("ret: %s\n", formatC(ret));
    // printf("FACPI*-1*(p1*p1+p2*p2): %s\n", formatC(FACPI * -1 * (p1 * p1 + p2
    // * p2)));
    puts("");
    puts("");
    return ret;
}

double complex delta11(double complex E, double complex p1, double complex p2,
                       double m0, double complex m1, double complex m2,
                       double complex m3, double complex m4) {
    auto A1 = p1 * p1 + p2 * p2 + m0 * m0;
    printf("A1: %s\n\n", formatC(A1));
    auto B1 = 2 * p1 * p2;
    printf("B1: %s\n\n", formatC(B1));
    auto C1 = m1 + m3 + p1 * p1 / 2 / m1 + p2 * p2 / 2 / m3 - E - I * EPSILON;
    printf("C1: %s\n\n", formatC(C1));
    auto D1 = m2 + m4 + p1 * p1 / 2 / m2 + p2 * p2 / 2 / m4 - E - I * EPSILON;
    printf("D1: %s\n\n", formatC(D1));
    auto a1 = xsqrtright(A1 - B1);
    printf("a1: %s\n\n", formatC(a1));
    auto b1 = xsqrtright(A1 + B1);
    printf("b1: %s\n\n", formatC(b1));
    auto CplusDtimesaminusb = (C1 + D1) * (a1 - b1);
    printf("(C+D)*(a-b): %s\n\n", formatC(CplusDtimesaminusb));
    auto log1 = (A1 - C1 * C1) * xlog((a1 + C1) / (b1 + C1));
    printf("log1: %s\n\n", formatC(log1));
    auto log2 = (A1 - D1 * D1) * xlog((a1 + D1) / (b1 + D1));
    printf("log2: %s\n\n", formatC(log2));
    auto ret = -(CplusDtimesaminusb + 2 * B1 + log1 + log2) / B1 / B1;
    printf("ret: %s\n\n", formatC(ret));
    return ret;
}

void testconvergence() {
    size_t Ngauss = 512;
    double Lambda;
    WaveFunction *wf [[gnu::cleanup(wffree)]] = WFnew(partialwave, 4, Ngauss);
    auto self = wf;
    double complex p = 1.3 * I;
    size_t n = 1;
    do {
        scanf("%lf", &Lambda);
        wf_refresh(wf, Lambda);
        const uint64_t Ngauss = self->rNgauss;
        const uint64_t l = self->l;
        double *xi = self->xi;
        double *wi = self->wi;
        double complex psi = 0.0;
        auto nidx = 40;
        double complex quad = 0 + 0 * I;
        printf("nu: %10.6e\n", nu_n(nidx + 1));
        for (size_t i = 0; i < Ngauss; i += 1) {
            quad += integrand_complex(xi[i], p, nidx + 1, l) * wi[i];
        }
        printf("%10.6e%+10.6eim  ", creal(quad), cimag(quad));
        if (((nidx + 1) % 8) == 0) {
            puts("");
        }
        psi += gsl_matrix_get(self->c_solution, nidx, n - 1) *
               N_nl(nidx + 1, l) * quad;
        printf("%10.4e%+10.4eim\n", creal(psi), cimag(psi));
    } while (true);
}

double *linspace(double start, double end, size_t len) {
    double *res = (double *)malloc(sizeof(double) * len);
    double step = (end - start) / len;
    for (size_t i = 0; i < len; i += 1) {
        res[i] = start + step * i;
    }
    return res;
}

void testwf() {
    // WaveFunction *wf __attribute__((cleanup(auto_wffree))) = WFnew(0, 20,
    // 64); for (size_t i = 0; i < N_MAX; i += 1) {
    //   for (size_t j = 0; j < N_MAX; j += 1) {
    //     const size_t idx = i * wf->c_solution->tda + j;
    //     printf("%.4e\t", wf->c_solution->data[idx]);
    //   }
    //   puts("");
    // }
    size_t pNgauss = 64;
    LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(pNgauss, 4, 1e-6);
    double complex(*psi)[N_MAX + 1][pNgauss + 1] = lse->psi_n_mat;
    lse_refresh(lse, -0.3, (double[]){1, 1, 1, 1}, PP);
    // printf("%25s %28s\n", "Column 1", "Column 2");
    // char delim[55];
    // for (size_t i = 0; i < 54; i += 1) {
    //   delim[i] = '-';
    // }
    // delim[54] = '\0';
    // printf("%s\n", delim);
    const int width = 11;    // Field width for each number
    const int precision = 4; // Decimal places

    for (size_t i = 0; i < N_MAX + 1; i += 1) {
        auto x = psi[0][i][pNgauss];
        auto y = psi[1][i][pNgauss];
        printf("(%*.*e %*.*ei)   ", width, precision, creal(x), width,
               precision, cimag(x));

        // Format second complex number
        printf("(%*.*e %*.*ei)\n", width, precision, creal(y), width, precision,
               cimag(y));
    }
    // puts("energy");
    // for (size_t i = 0; i < N_MAX; i += 1) {
    //   auto e = lse->E_vec[i];
    //   printf("%*.*e\n", width, precision, e);
    // }
}

void printmat(matrix *m) {
    auto row = m->size1;
    auto col = m->size2;
    for (size_t i = 0; i < row; i += 1) {
        for (size_t j = 0; j < col; j += 1) {
            auto val = matrix_get(m, i, j);
            printf("%5.2e%+5.2eim  ", creal(val), cimag(val));
        }
        puts("");
    }
}

void printabsmat(matrix *m) {
    auto row = m->size1;
    auto col = m->size2;
    for (size_t i = 0; i < row; i += 1) {
        for (size_t j = 0; j < col; j += 1) {
            auto x = matrix_get(m, i, j);
            printf("%9.2e ", cabs(x));
        }
        puts("");
    }
}

void testwopotential() {
    LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(5, 4, 1e-6);
    lse_compute(lse, 0.1, (double[]){1, 1, 1, 1}, 1);
    // lse_X(lse);
    // lse_XtX(lse);
    auto pNgauss = lse->pNgauss;
    auto psi = (double complex(*)[N_MAX + 1][pNgauss + 1]) lse->psi_n_mat;
    auto v = V_QM_00(lse, pNgauss, pNgauss);
    printf("%f %f(im)\n", creal(v), cimag(v));
    // for (size_t i = 0; i < N_MAX + 1; i += 1) {
    //   for (size_t pi = 0; pi < pNgauss + 1; pi += 1) {
    //     auto x = psi[0][i][pi];
    //     printf("%f %f(im) ", creal(x), cimag(x));
    //   }
    //   puts("");
    // }
    // for (size_t i = 0; i < 2; i += 1) {
    //   for (size_t j = 0; j < 2; j += 1) {
    //     size_t idx = (pNgauss + i * (pNgauss + 1)) * 2 * (pNgauss + 1) +
    //     pNgauss +
    //                  j * (pNgauss + 1);
    //     auto ptr = lse->TOME->data;
    //     double complex x = ptr[2 * idx] + ptr[2 * idx + 1] * I;
    //     printf("%f %f(im) ", creal(x), cimag(x));
    //   }
    // }
}

double complex qmhelp(LSE *lse, size_t pi, size_t ppi) {
    double complex res = 0;
    auto E = lse->E;
    size_t pNgauss = lse->pNgauss;
    auto psi = (double complex(*)[N_MAX + 1][pNgauss + 1]) lse->psi_n_mat;
    size_t chan0 = pi / (pNgauss + 1);
    size_t chan1 = ppi / (pNgauss + 1);
    printf("chan0: %lu, chan1: %lu\n", chan0, chan1);
    for (size_t i = 0; i < N_MAX; i += 1) {
        auto x = psi[chan0][i][pi % (pNgauss + 1)] *
                 conj(psi[chan1][i][ppi % (pNgauss + 1)]) / (E - lse->E_vec[i]);
        res += x;
        if (i < 5) {
            printf("(i = %lu) x: %f %f(im)\n", i, creal(x), cimag(x));
            printf("\tE - lse->E_vec[i]: %f %f(im)\n", creal(E - lse->E_vec[i]),
                   cimag(E - lse->E_vec[i]));
            auto si = psi[chan0][i][pi % (pNgauss + 1)];
            auto sj = conj(psi[chan1][i][ppi % (pNgauss + 1)]);
            printf("\tpsi[chan0]: %f %f(im)\n", creal(si), cimag(si));
            printf("\tpsi[chan1]: %f %f(im)\n", creal(sj), cimag(sj));
            printf("\tchan1 = %lu, i = %lu, ppi %% (pNgauss + 1) = "
                   "%lu\n",
                   chan1, i, ppi % (pNgauss + 1));
            puts("");
        }
    }
    return res;
}

void testome() {
    struct OME *ome = malloc(sizeof(*ome));
    ome_build(ome);
    auto x = OME_01(*ome, 0.1, 0.3, 0.8);
    printf("%f%+f\n", creal(x), cimag(x));
}

void testlse() {
    const size_t Ngauss = 4;
    LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(Ngauss, 4, 1e-9);
    if (!lse) {
        fprintf(stderr, "Failed to create LSE solver\n");
        exit(1);
    }
    lse_compute(lse, 0.1, (double[]){0, 0, 0, 0}, 3);
    auto row = Ngauss + 1;
    auto col = Ngauss + 1;
    auto m = lse->VOME;
    for (size_t i = 0; i < row; i += 1) {
        for (size_t j = 0; j < col; j += 1) {
            auto val = matrix_get(m, i + Ngauss + 1, j);
            printf("%5.2e%+5.2eim  ", creal(val), cimag(val));
        }
        puts("");
    }
    // auto Espace = linspace(-2, -0.3, 8);
    // for (size_t Ei = 0; Ei < 1; Ei += 1) {
    //   auto E = Espace[Ei];
    //   lse_compute(lse, E, (double[]){1, 1, 1, 1}, 1);
    //   // lse_gmat(lse);
    //   // lse_vmat(lse);
    //   // puts("G matrix");
    //   // gsl_matrix_complex_fprintf(stdout, lse->G, "%.2e");
    //   puts("V matrix");
    //   printabsmat(lse->VOME);
    //   puts("T matrix");
    //   printabsmat(lse->TOME);
    //   puts("onshell T");
    //   for (size_t alpha = 0; alpha < 2; alpha += 1) {
    //     for (size_t beta = 0; beta < 2; beta += 1) {
    //       auto x = matrix_get(lse->TOME, Ngauss + alpha * (Ngauss + 1),
    //                           Ngauss + beta * (Ngauss + 1));
    //       printf("%lu%lu: %.3e %.3e(im)\n", alpha, beta, creal(x),
    //       cimag(x));
    //     }
    //   }
    // }
    // puts("det");
    // auto x = lse_detImVG(lse, E);
    // printf("%f %f(im)\n", creal(x), cimag(x));
    // for (size_t pi = 0; pi < Ngauss + 1; pi += 1) {
    //   auto x = V_QM_00(lse, pi, Ngauss);
    //   printf("%f %f(im)\n", creal(x), cimag(x));
    // }
    // auto s = qmhelp(lse, 0, Ngauss);
    // gsl_matrix_complex_fprintf(stdout, lse->T, "%.2e");
    // puts("I - VG matrix");
    // gsl_matrix_complex_fprintf(stdout, lse->iIVG, "%.2e");
    // int result = lse_run(lse);
    // if (result != 0) {
    //   fprintf(stderr, "LSE solver failed with error code %d\n", result);
    //   lse_deinit(lse);
    //   return 1;
    // }
    //
    // printf("LSE solver ran successfully\n");
}

void testonshellpsi() {
    const size_t Ngauss = 2;
    LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(Ngauss, 4, 1e-9);
    if (!lse) {
        fprintf(stderr, "Failed to create LSE solver\n");
        exit(1);
    }
    double complex E = -0.8;
    puts("");
    puts("");
    puts("");
    puts("");
    lse_refresh(lse, E, (double[]){1, 1, 1, 1}, 1);
}

void testonshell() {
    const size_t Ngauss = 2;
    LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(Ngauss, 4, 1e-9);
    if (!lse) {
        fprintf(stderr, "Failed to create LSE solver\n");
        exit(1);
    }
    lse_refresh(lse, -0.578802970, (double[]){1, 1, 1, 1}, 1);
    __auto_type q = lse->x0[0];
    printf("(%.12e) + Im(%.12e)\n", creal(q), cimag(q));
    q = lse->x0[1];
    printf("(%.12e) + Im(%.12e)\n", creal(q), cimag(q));
    lse_vmat(lse);
    __auto_type onshellV = matrix_get(lse->VOME, Ngauss, Ngauss);
    printf("(%.12e) + Im(%.12e)\n", creal(onshellV), cimag(onshellV));
    lse_refresh(lse, -0.578802965, (double[]){1, 1, 1, 1}, 1);
    lse_vmat(lse);
    q = lse->x0[0];
    printf("(%.12e) + Im(%.12e)\n", creal(q), cimag(q));
    q = lse->x0[1];
    printf("(%.12e) + Im(%.12e)\n", creal(q), cimag(q));
    onshellV = matrix_get(lse->VOME, Ngauss, Ngauss);
    printf("(%.12e) + Im(%.12e)\n", creal(onshellV), cimag(onshellV));
}

double complex Ofunc(double complex E, double complex p,
                     double complex pprime) {
    auto m = m_pi;
    return 4 * fsquare(g_b) / fsquare(f_pi) * -1. / 4. / p / pprime *
           (clog(
                (E - (m + csquare(p - pprime) / 2 / m) - omega_00(p, pprime)) /
                (E - (m + csquare(p + pprime) / 2 / m) - omega_00(p, pprime))) +
            clog((E - (m + csquare(p - pprime) / 2 / m) -
                  omegastar_00(p, pprime)) /
                 (E - (m + csquare(p + pprime) / 2 / m) -
                  omegastar_00(p, pprime))));
}

void unitest() {
    size_t pNgauss = 64;
    double Lambda = 4;
    double epsilon = 1e-6;
    double E = 1.2;
    LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(pNgauss, Lambda, epsilon);
    lse_compute(lse, E, (double[]){0, 0, 0, 0}, PP);
    printf("m_B: %.8f  det: %.8f\n", m_B, cabs(lse->det));
    // E = 0.2;
    E += m11 + m12;
    // 4 * square(g_pi) / square(f_pi) * -1 / 4 / p / pprime *
    // (clog((E - (m + csquare(p - pprime) / 2 / m) - omega_00(p, pprime)) /
    //       (E - (m + csquare(p + pprime) / 2 / m) - omega_00(p, pprime)))
    //       +
    //  clog(
    //      (E - (m + csquare(p - pprime) / 2 / m) - omegaprime_00(p,
    //      pprime)) / (E - (m + csquare(p + pprime) / 2 / m) -
    //      omegaprime_00(p, pprime))));
}
void ptrfree(double **ptr) { free(*ptr); }
void testscript() {
    size_t len = 32;
    double *E [[gnu::cleanup(ptrfree)]] = malloc(sizeof(double[len]));
    double start = -2;
    double end = -0.31;
    double step = (end - start) / len;
    for (size_t i = 0; i < len; i += 1) {
        E[i] = start + step * i;
    }
    // free(onshellT(E, len, 32, 4, 1e-6));
}

void testpole() {
    double Er[3] = {-1, 0, 1};
    double Ei[3] = {-1, 0, 1};
    free(Poles(Er, 3, Ei, 3, (double[]){1, 1, 1, 1}, 64, 4, 1e-6));
}

void test() {
    double Lambda = 2;
    size_t Ngauss = 64;
    LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(Ngauss, Lambda, 0);
    double E = 0.1;
    // lse_vmat(lse);
    printf("Energy: %s\n", formatC(E + m11 + m12));
    lse_refresh(lse, E, (double[4]){0, 0, 0, 0}, PP);
    lse_gmat(lse);
    E += 2 * m_B + m_pi;
    // E += -0.3+m11 + m12;
    double complex p1;
    double complex p2;
    p1 = lse->xi[12];
    p2 = lse->xi[1];
    printf("EE: %s\n", formatC(E));
    printf("p1: %s\n", formatC(p1));
    printf("p2: %s\n", formatC(p2));
    double complex v1;
    double complex temp;
    double complex v0;
    double complex quadresult, anaresult;
    puts("----------------------\n");
    quadresult =
        Vpiu(lse->ome, E, p1, p2, m_B_star, m_B, m_B_star, m_B, m_pi, 1);
    anaresult = ANA_00(E, p1, p2, m_pi);
    printf("V: %s\n", formatC(anaresult));
    puts("---------------------\n");
    // printf("quadrature: %s\n", formatC(quadresult));
    // printf("analytic: %s\n", formatC(anaresult));
    // puts("----------------------\n");
    auto Delta0 = delta01(E + 0I, p1, p2, m_pi, m_B, m_B_star, m_B, m_B_star);
    Delta0 = Delta0_00(E, p1, p2, m_pi);
    puts("\n\n----------Delta1---------------\n");
    auto Delta1 = delta11(E + 0I, p1, p2, m_pi, m_B, m_B_star, m_B, m_B_star);
    Delta1 = Delta1_00(E, p1, p2, m_pi);
    auto v = g_b * g_b / 24 / f_pi / f_pi *
             (2 * p1 * p2 * Delta1 - (p1 * p1 + p2 * p2) * Delta0);
    printf("------------V-------------\n%s\n", formatC(v));

    printf("Delta0: %s\n", formatC(Delta0));
    printf("Delta1: %s\n", formatC(Delta1));
    //   for (size_t i = 0; i < 2 * DIMIM + DIMRE; i += 1) {
    //       temp = TOPTintegrand(E, lse->ome.xxpiup[i], p1, p2, m_B, m_B,
    //       m_pi);
    // printf("%s\n", formatC(temp));
    //   }
    size_t ngauss = 64;
    gsl_integration_glfixed_table *t =
        gsl_integration_glfixed_table_alloc(ngauss);
    double complex res = 0;
    double x, w;
    printf("explict quadrature: %s\n", formatC(res));
    gsl_integration_glfixed_table_free(t);
}

void testv() {
    double Lambda = 4;
    size_t Ngauss = 64;
    LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(Ngauss, Lambda, 1e-6);
    lse_refresh(lse, -0.1, (double[4]){0, 0, 0, 0}, PP);
    lse_vmat(lse);
    getV(-0.1, Ngauss, 4, 1e-6);
}

int main(int argc, char *argv[]) {
    if (argc == 1) {
        puts("what do you want?");
        return EXIT_SUCCESS;
    }
    if (strcmp(argv[1], "lse") == 0) {
        testlse();
    } else if (strcmp(argv[1], "wf") == 0) {
        testwf();
    } else if (strcmp(argv[1], "unit") == 0) {
        unitest();
    } else if (strcmp(argv[1], "onshell") == 0) {
        testonshell();
    } else if (strcmp(argv[1], "script") == 0) {
        testscript();
    } else if (strcmp(argv[1], "twop") == 0) {
        testwopotential();
    } else if (strcmp(argv[1], "onshellpsi") == 0) {
        testonshellpsi();
    } else if (strcmp(argv[1], "testpole") == 0) {
        testpole();
    } else if (strcmp(argv[1], "convergence") == 0) {
        testconvergence();
    } else if (strcmp(argv[1], "testome") == 0) {
        testome();
    } else if (strcmp(argv[1], "repl") == 0) {
        test();
    } else if (strcmp(argv[1], "v") == 0) {
        testv();
    }
    return EXIT_SUCCESS;
}

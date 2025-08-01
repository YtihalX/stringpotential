#include "autofree.h"
#include "lse.h"
#include "ome.h"
#include "script.h"
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

void testconvergence()
{
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
	    psi += gsl_matrix_get(self->c_solution, nidx, n - 1) * N_nl(nidx + 1, l) * quad;
	    printf("%10.4e%+10.4eim\n", creal(psi), cimag(psi));
      } while (true);
}

double *linspace(double start, double end, size_t len)
{
      double *res = (double *)malloc(sizeof(double) * len);
      double step = (end - start) / len;
      for (size_t i = 0; i < len; i += 1) {
	    res[i] = start + step * i;
      }
      return res;
}

void testwf()
{
      // WaveFunction *wf __attribute__((cleanup(auto_wffree))) = WFnew(0, 20,
      // 64); for (size_t i = 0; i < N_MAX; i += 1) {
      //   for (size_t j = 0; j < N_MAX; j += 1) {
      //     const size_t idx = i * wf->c_solution->tda + j;
      //     printf("%.4e\t", wf->c_solution->data[idx]);
      //   }
      //   puts("");
      // }
      size_t pNgauss = 64;
      LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(pNgauss, 4, 1e-6, 0);
      double complex(*psi)[N_MAX + 1][pNgauss + 1] = lse->psi_n_mat;
      lse_refresh(lse, -0.3, (double[]) { 1, 1, 1, 1 }, PP, 0);
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
	    printf("(%*.*e %*.*ei)   ", width, precision, creal(x), width, precision, cimag(x));

	    // Format second complex number
	    printf("(%*.*e %*.*ei)\n", width, precision, creal(y), width, precision, cimag(y));
      }
      // puts("energy");
      // for (size_t i = 0; i < N_MAX; i += 1) {
      //   auto e = lse->E_vec[i];
      //   printf("%*.*e\n", width, precision, e);
      // }
}

void printmat(matrix *m)
{
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

void printabsmat(matrix *m)
{
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

void testwopotential()
{
      LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(5, 4, 1e-6, 0);
      lse_compute(lse, 0.1, (double[]) { 1, 1, 1, 1 }, 1, 0);
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

double complex qmhelp(LSE *lse, size_t pi, size_t ppi)
{
      double complex res = 0;
      auto E = lse->E;
      size_t pNgauss = lse->pNgauss;
      auto psi = (double complex(*)[N_MAX + 1][pNgauss + 1]) lse->psi_n_mat;
      size_t chan0 = pi / (pNgauss + 1);
      size_t chan1 = ppi / (pNgauss + 1);
      printf("chan0: %lu, chan1: %lu\n", chan0, chan1);
      for (size_t i = 0; i < N_MAX; i += 1) {
	    auto x = psi[chan0][i][pi % (pNgauss + 1)] * conj(psi[chan1][i][ppi % (pNgauss + 1)]) / (E - lse->E_vec[i]);
	    res += x;
	    if (i < 5) {
		  printf("(i = %lu) x: %f %f(im)\n", i, creal(x), cimag(x));
		  printf("\tE - lse->E_vec[i]: %f %f(im)\n", creal(E - lse->E_vec[i]), cimag(E - lse->E_vec[i]));
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

void testome()
{
      struct OME *ome = malloc(sizeof(*ome));
      ome_build(ome);
      auto x = OME_01(*ome, 0.1, 0.3, 0.8);
      printf("%f%+f\n", creal(x), cimag(x));
}

void testlse()
{
      const size_t Ngauss = 4;
      LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(Ngauss, 4, 1e-9, 0);
      if (!lse) {
	    fprintf(stderr, "Failed to create LSE solver\n");
	    exit(1);
      }
      lse_compute(lse, 0.1, (double[]) { 0, 0, 0, 0 }, 3, 0);
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

void testonshellpsi()
{
      const size_t Ngauss = 2;
      LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(Ngauss, 4, 1e-9, 0);
      if (!lse) {
	    fprintf(stderr, "Failed to create LSE solver\n");
	    exit(1);
      }
      double complex E = -0.8;
      puts("");
      puts("");
      puts("");
      puts("");
      lse_refresh(lse, E, (double[]) { 1, 1, 1, 1 }, 1, 0);
}

void testonshell()
{
      const size_t Ngauss = 2;
      LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(Ngauss, 4, 1e-9, 0);
      if (!lse) {
	    fprintf(stderr, "Failed to create LSE solver\n");
	    exit(1);
      }
      lse_refresh(lse, -0.578802970, (double[]) { 1, 1, 1, 1 }, 1, 0);
      __auto_type q = lse->x0[0];
      printf("(%.12e) + Im(%.12e)\n", creal(q), cimag(q));
      q = lse->x0[1];
      printf("(%.12e) + Im(%.12e)\n", creal(q), cimag(q));
      lse_vmat(lse);
      __auto_type onshellV = matrix_get(lse->VOME, Ngauss, Ngauss);
      printf("(%.12e) + Im(%.12e)\n", creal(onshellV), cimag(onshellV));
      lse_refresh(lse, -0.578802965, (double[]) { 1, 1, 1, 1 }, 1, 0);
      lse_vmat(lse);
      q = lse->x0[0];
      printf("(%.12e) + Im(%.12e)\n", creal(q), cimag(q));
      q = lse->x0[1];
      printf("(%.12e) + Im(%.12e)\n", creal(q), cimag(q));
      onshellV = matrix_get(lse->VOME, Ngauss, Ngauss);
      printf("(%.12e) + Im(%.12e)\n", creal(onshellV), cimag(onshellV));
}

double complex Ofunc(double complex E, double complex p, double complex pprime)
{
      auto m = m_pi;
      return 4 * fsquare(g_b) / fsquare(f_pi) * -1. / 4. / p / pprime *
	     (clog((E - (m + csquare(p - pprime) / 2 / m) - omega_00(p, pprime, 0)) /
		   (E - (m + csquare(p + pprime) / 2 / m) - omega_00(p, pprime, 0))) +
	      clog((E - (m + csquare(p - pprime) / 2 / m) - omegaprime_00(p, pprime, 0)) /
		   (E - (m + csquare(p + pprime) / 2 / m) - omegaprime_00(p, pprime, 0))));
}

void unitest()
{
      size_t pNgauss = 64;
      double Lambda = 4;
      double epsilon = 1e-6;
      double E = 1.2;
      LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(pNgauss, Lambda, epsilon, 0);
      lse_compute(lse, E, (double[]) { 0, 0, 0, 0 }, PP, 0);
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
void testscript()
{
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

void testpole()
{
      double Er[3] = { -1, 0, 1 };
      double Ei[3] = { -1, 0, 1 };
      free(Poles(Er, 3, Ei, 3, (double[]) { 1, 1, 1, 1 }, 64, 4, 1e-6));
}

int main(int argc, char *argv[])
{
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
      }
      return EXIT_SUCCESS;
}

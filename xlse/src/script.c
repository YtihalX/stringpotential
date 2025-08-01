#include "script.h"
#include "autofree.h"
#include "lse.h"
#include "ome.h"
#include "wavefunction.h"
#include <gsl/gsl_matrix_complex_double.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <threads.h>

double complex *onshellT(double *E, size_t len, double C[4], size_t pNgauss, double Lambda, double epsilon)
{
      // puts("in C");
      // for (size_t i = 0; i < len; i += 1) {
      //   printf("%f,  ", E[i]);
      // }
      // puts("");
      double complex *res = malloc(sizeof(double complex) * len * 4);
      onshellElements onshellBuf = {
	    .ose00 = res,
	    .ose01 = res + len,
	    .ose10 = res + 2 * len,
	    .ose11 = res + 3 * len,
      };
      thrd_t tid[NTHREADS];
      argstruct args[NTHREADS];
      size_t ntasks = len / NTHREADS;
      size_t residue = len % NTHREADS;
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    args[i].pNgauss = pNgauss;
	    args[i].Lambda = Lambda;
	    args[i].epsilon = epsilon;
	    args[i].res = &onshellBuf;
	    args[i].rs = PP, args[i].E = E;
	    args[i].id = i;
	    for (size_t cc = 0; cc < 4; cc += 1) {
		  args[i].C[cc] = C[cc];
	    }
	    if (i < residue) {
		  args[i].len = ntasks + 1;
		  args[i].start = (ntasks + 1) * i;
	    } else {
		  args[i].len = ntasks;
		  args[i].start = ntasks * i + residue;
	    }
      }
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    thrd_create(&tid[i], oT, &args[i]);
      }
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    thrd_join(tid[i], NULL);
      }
      return res;
}

double complex *onshellG(double *E, size_t len, double C[4], size_t pNgauss, double Lambda, double epsilon)
{
      double complex *res = malloc(sizeof(double complex) * len * 4);
      onshellElements onshellBuf = {
	    .ose00 = res,
	    .ose01 = res + len,
	    .ose10 = res + 2 * len,
	    .ose11 = res + 3 * len,
      };
      thrd_t tid[NTHREADS];
      argstruct args[NTHREADS];
      size_t ntasks = len / NTHREADS;
      size_t residue = len % NTHREADS;
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    args[i].pNgauss = pNgauss;
	    args[i].Lambda = Lambda;
	    args[i].epsilon = epsilon;
	    args[i].res = &onshellBuf;
	    args[i].rs = PP, args[i].E = E;
	    args[i].id = i;
	    for (size_t cc = 0; cc < 4; cc += 1) {
		  args[i].C[cc] = C[cc];
	    }
	    if (i < residue) {
		  args[i].len = ntasks + 1;
		  args[i].start = (ntasks + 1) * i;
	    } else {
		  args[i].len = ntasks;
		  args[i].start = ntasks * i + residue;
	    }
      }
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    thrd_create(&tid[i], oG, &args[i]);
      }
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    thrd_join(tid[i], NULL);
      }
      return res;
}

double complex *onshellV(double *E, size_t len, double C[4], size_t pNgauss, double Lambda, double epsilon)
{
      // puts("in C");
      // for (size_t i = 0; i < len; i += 1) {
      //   printf("%f,  ", E[i]);
      // }
      // puts("");
      double complex *res = malloc(sizeof(double complex) * len * 4);
      onshellElements onshellBuf = {
	    .ose00 = res,
	    .ose01 = res + len,
	    .ose10 = res + 2 * len,
	    .ose11 = res + 3 * len,
      };
      thrd_t tid[NTHREADS];
      argstruct args[NTHREADS];
      size_t ntasks = len / NTHREADS;
      size_t residue = len % NTHREADS;
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    args[i].pNgauss = pNgauss;
	    args[i].Lambda = Lambda;
	    args[i].epsilon = epsilon;
	    args[i].res = &onshellBuf;
	    args[i].rs = PP, args[i].E = E;
	    args[i].id = i;
	    for (size_t cc = 0; cc < 4; cc += 1) {
		  args[i].C[cc] = C[cc];
	    }
	    if (i < residue) {
		  args[i].len = ntasks + 1;
		  args[i].start = (ntasks + 1) * i;
	    } else {
		  args[i].len = ntasks;
		  args[i].start = ntasks * i + residue;
	    }
      }
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    thrd_create(&tid[i], oV, &args[i]);
      }
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    thrd_join(tid[i], NULL);
      }
      return res;
}

double complex *onshellTV(double *E, size_t len, double C[4], size_t pNgauss, double Lambda, double epsilon)
{
      double complex *res = malloc(sizeof(double complex) * len * 8);
      thrd_t tid[NTHREADS];
      argstruct args[NTHREADS];
      size_t ntasks = len / NTHREADS;
      size_t residue = len % NTHREADS;
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    args[i].pNgauss = pNgauss;
	    args[i].Lambda = Lambda;
	    args[i].epsilon = epsilon;
	    args[i].res = res;
	    args[i].size = len;
	    args[i].rs = PP;
	    args[i].E = E;
	    args[i].id = i;
	    for (size_t cc = 0; cc < 4; cc += 1) {
		  args[i].C[cc] = C[cc];
	    }
	    if (i < residue) {
		  args[i].len = ntasks + 1;
		  args[i].start = (ntasks + 1) * i;
	    } else {
		  args[i].len = ntasks;
		  args[i].start = ntasks * i + residue;
	    }
      }
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    thrd_create(&tid[i], oTV, &args[i]);
      }
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    thrd_join(tid[i], NULL);
      }
      return res;
}

struct OME *ome_malloc()
{
      struct OME *ome = malloc(sizeof(*ome));
      ome_build(ome);
      return ome;
}

double complex V(struct OME *ome, double E, double p, double pprime) { return OMEANA_00(E, p, pprime, 0); }

double complex *traceG(double *E, size_t len, double C[4], size_t pNgauss, double Lambda, double epsilon)
{
      double complex *res = malloc(sizeof(double complex) * len * 2);
      thrd_t tid[NTHREADS];
      argstruct args[NTHREADS];
      size_t ntasks = len / NTHREADS;
      size_t residue = len % NTHREADS;
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    args[i].pNgauss = pNgauss;
	    args[i].Lambda = Lambda;
	    args[i].epsilon = epsilon;
	    args[i].res = res;
	    args[i].rs = PP;
	    args[i].E = E;
	    args[i].id = i;
	    args[i].size = len;
	    for (size_t cc = 0; cc < 4; cc += 1) {
		  args[i].C[cc] = C[cc];
	    }
	    if (i < residue) {
		  args[i].len = ntasks + 1;
		  args[i].start = (ntasks + 1) * i;
	    } else {
		  args[i].len = ntasks;
		  args[i].start = ntasks * i + residue;
	    }
      }
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    thrd_create(&tid[i], trG, &args[i]);
      }
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    thrd_join(tid[i], NULL);
      }
      return res;
}

double complex *Det(double *E, size_t len, double C[4], size_t pNgauss, double Lambda, double epsilon)
{
      double complex *res = malloc(sizeof(double complex) * len);
      thrd_t tid[NTHREADS];
      argstruct args[NTHREADS];
      size_t ntasks = len / NTHREADS;
      size_t residue = len % NTHREADS;
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    args[i].pNgauss = pNgauss;
	    args[i].Lambda = Lambda;
	    args[i].epsilon = epsilon;
	    args[i].res = res;
	    args[i].rs = PP;
	    args[i].E = E;
	    args[i].id = i;
	    for (size_t cc = 0; cc < 4; cc += 1) {
		  args[i].C[cc] = C[cc];
	    }
	    if (i < residue) {
		  args[i].len = ntasks + 1;
		  args[i].start = (ntasks + 1) * i;
	    } else {
		  args[i].len = ntasks;
		  args[i].start = ntasks * i + residue;
	    }
      }
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    thrd_create(&tid[i], det, &args[i]);
      }
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    thrd_join(tid[i], NULL);
      }
      return res;
}

double complex *mass(double *M, size_t len, double E, double C[4], size_t pNgauss, double Lambda, double epsilon)
{
      double complex *res = malloc(sizeof(double complex) * len);
      thrd_t tid[NTHREADS];
      argstruct args[NTHREADS];
      size_t ntasks = len / NTHREADS;
      size_t residue = len % NTHREADS;
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    args[i].pNgauss = pNgauss;
	    args[i].Lambda = Lambda;
	    args[i].epsilon = epsilon;
	    args[i].res = res;
	    args[i].rs = PP;
	    args[i].E = M;
	    args[i].id = i;
	    for (size_t cc = 0; cc < 4; cc += 1) {
		  args[i].C[cc] = C[cc];
	    }
	    args[i].C[0] = E;
	    if (i < residue) {
		  args[i].len = ntasks + 1;
		  args[i].start = (ntasks + 1) * i;
	    } else {
		  args[i].len = ntasks;
		  args[i].start = ntasks * i + residue;
	    }
      }
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    thrd_create(&tid[i], ma, &args[i]);
      }
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    thrd_join(tid[i], NULL);
      }
      return res;
}

double complex *Both(double *E, size_t len, double C[4], size_t pNgauss, double Lambda, double epsilon)
{
      double complex *res = malloc(sizeof(double complex) * len * 5);
      thrd_t tid[NTHREADS];
      argstruct args[NTHREADS];
      size_t ntasks = len / NTHREADS;
      size_t residue = len % NTHREADS;
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    args[i].pNgauss = pNgauss;
	    args[i].Lambda = Lambda;
	    args[i].epsilon = epsilon;
	    args[i].res = res;
	    args[i].rs = PP;
	    args[i].E = E;
	    args[i].id = i;
	    for (size_t cc = 0; cc < 4; cc += 1) {
		  args[i].C[cc] = C[cc];
	    }
	    if (i < residue) {
		  args[i].len = ntasks + 1;
		  args[i].start = (ntasks + 1) * i;
	    } else {
		  args[i].len = ntasks;
		  args[i].start = ntasks * i + residue;
	    }
      }
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    thrd_create(&tid[i], both, &args[i]);
      }
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    thrd_join(tid[i], NULL);
      }
      return res;
}

void Evec(double *E, double V0)
{
      WaveFunction *wf [[gnu::cleanup(wffree)]] = WFnew(partialwave, RLAMBDA, RNGAUSS);
      // LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(10, 4, 1e-6);
      for (size_t i = 0; i < N_TOWER; i += 1) {
	    E[i] = wf->E_solution->data[i] + V0;
      }
}

void Ntower(size_t *len) { *len = N_TOWER; }

double complex *Poles(double *Er, size_t rlen, double *Ei, size_t ilen, double C[4], size_t pNgauss, double Lambda,
		      double epsilon)
{
      size_t len = rlen * ilen;
      double complex *res = malloc(sizeof(double complex) * len * 4);
      thrd_t tid[NTHREADS];
      struct polestruct args[NTHREADS];
      size_t ntasks = len / NTHREADS;
      size_t residue = len % NTHREADS;
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    args[i].pNgauss = pNgauss;
	    args[i].Lambda = Lambda;
	    args[i].epsilon = epsilon;
	    args[i].res = res;
	    args[i].rs = PP;
	    args[i].Er = Er;
	    args[i].rlen = rlen;
	    args[i].Ei = Ei;
	    args[i].ilen = ilen;
	    args[i].id = i;
	    for (size_t cc = 0; cc < 4; cc += 1) {
		  args[i].C[cc] = C[cc];
	    }
	    if (i < residue) {
		  args[i].len = ntasks + 1;
		  args[i].start = (ntasks + 1) * i;
	    } else {
		  args[i].len = ntasks;
		  args[i].start = ntasks * i + residue;
	    }
      }
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    thrd_create(&tid[i], poles, &args[i]);
      }
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    thrd_join(tid[i], NULL);
      }
      return res;
}

double *Fit(double C[4], size_t pNgauss, double Lambda, double epsilon)
{
      LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(pNgauss, Lambda, epsilon, 0);
      return minimize(lse, C);
}

double *Cost(double *c, size_t len, double complex resonance, size_t pNgauss, double Lambda, double epsilon)
{
      double *res = malloc(sizeof(double) * len);
      thrd_t tid[NTHREADS];
      argstruct args[NTHREADS];
      size_t ntasks = len / NTHREADS;
      size_t residue = len % NTHREADS;
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    args[i].pNgauss = pNgauss;
	    args[i].Lambda = Lambda;
	    args[i].epsilon = epsilon;
	    args[i].res = res;
	    args[i].rs = PP;
	    args[i].E = c;
	    args[i].id = i;
	    if (i < residue) {
		  args[i].len = ntasks + 1;
		  args[i].start = (ntasks + 1) * i;
	    } else {
		  args[i].len = ntasks;
		  args[i].start = ntasks * i + residue;
	    }
      }
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    thrd_create(&tid[i], cst, &args[i]);
      }
      for (size_t i = 0; i < NTHREADS; i += 1) {
	    thrd_join(tid[i], NULL);
      }
      return res;
}

void ose_free(onshellElements ose) { free(ose.ose00); }

int oT(void *arg)
{
      argstruct foo = *(argstruct *)arg;
      LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon, 0);
      size_t ngauss = foo.pNgauss;
      onshellElements *res = (onshellElements *)foo.res;
      double complex *ose00 = (double complex *)res->ose00;
      double complex *ose01 = (double complex *)res->ose01;
      double complex *ose10 = (double complex *)res->ose10;
      double complex *ose11 = (double complex *)res->ose11;
      // printf("start: %lu, len: %lu\n", foo.start, foo.len);
      for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
	    lse_compute(lse, foo.E[i], foo.C, foo.rs, 0);
// printf("%f\n", foo.E[i]);
#ifdef TPO
	    lse_X(lse);
	    lse_XtX(lse);
	    ose00[i] = lse->onshellT[0][0];
	    ose01[i] = lse->onshellT[0][1];
	    ose10[i] = lse->onshellT[1][0];
	    ose11[i] = lse->onshellT[1][1];
#else
	    auto T = (double complex(*)[2 * ngauss + 2]) lse->TOME->data;
	    ose00[i] = T[ngauss][ngauss];
	    ose01[i] = T[ngauss][2 * ngauss + 1];
	    ose10[i] = T[2 * ngauss + 1][ngauss];
	    ose11[i] = T[2 * ngauss + 1][2 * ngauss + 1];
	    // size_t idx = ngauss * 2 * (ngauss + 1) + ngauss;
	    // ose00[i] = lse->TOME->data[2 * idx] + lse->TOME->data[2 * idx
	    // + 1] * I; idx += ngauss + 1; ose01[i] = lse->TOME->data[2 *
	    // idx] + lse->TOME->data[2 * idx + 1] * I; idx = (2 * ngauss +
	    // 1) * 2 * (ngauss + 1) + ngauss; ose10[i] = lse->TOME->data[2
	    // * idx] + lse->TOME->data[2 * idx + 1] * I; idx += ngauss + 1;
	    // ose11[i] = lse->TOME->data[2 * idx] + lse->TOME->data[2 * idx
	    // + 1] * I;
#endif
      }
      return 0;
}

int oG(void *arg)
{
      argstruct foo = *(argstruct *)arg;
      LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon, 0);
      size_t ngauss = foo.pNgauss;
      onshellElements *res = (onshellElements *)foo.res;
      double complex *ose00 = (double complex *)res->ose00;
      double complex *ose01 = (double complex *)res->ose01;
      double complex *ose10 = (double complex *)res->ose10;
      double complex *ose11 = (double complex *)res->ose11;
      // printf("start: %lu, len: %lu\n", foo.start, foo.len);
      for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
	    lse_refresh(lse, foo.E[i], foo.C, foo.rs, 0);
	    lse_gmat(lse);
	    auto G = (double complex(*)[2 * ngauss + 2]) lse->G->data;
	    ose00[i] = G[ngauss][ngauss];
	    ose01[i] = G[ngauss][2 * ngauss + 1];
	    ose10[i] = G[2 * ngauss + 1][ngauss];
	    ose11[i] = G[2 * ngauss + 1][2 * ngauss + 1];
	    // size_t idx = ngauss * 2 * (ngauss + 1) + ngauss;
	    // ose00[i] = lse->TOME->data[2 * idx] + lse->TOME->data[2 * idx
	    // + 1] * I; idx += ngauss + 1; ose01[i] = lse->TOME->data[2 *
	    // idx] + lse->TOME->data[2 * idx + 1] * I; idx = (2 * ngauss +
	    // 1) * 2 * (ngauss + 1) + ngauss; ose10[i] = lse->TOME->data[2
	    // * idx] + lse->TOME->data[2 * idx + 1] * I; idx += ngauss + 1;
	    // ose11[i] = lse->TOME->data[2 * idx] + lse->TOME->data[2 * idx
	    // + 1] * I;
      }
      return 0;
}

int oV(void *arg)
{
      argstruct foo = *(argstruct *)arg;
      LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon, 0);
      size_t ngauss = foo.pNgauss;
      if (foo.id == 0) {
	    printf("%f\n", lse->xi[0]);
      }
      onshellElements *res = (onshellElements *)foo.res;
      double complex *ose00 = (double complex *)res->ose00;
      double complex *ose01 = (double complex *)res->ose01;
      double complex *ose10 = (double complex *)res->ose10;
      double complex *ose11 = (double complex *)res->ose11;
      int64_t xoffset = -ngauss;
      int64_t yoffset = -ngauss;
      // printf("start: %lu, len: %lu\n", foo.start, foo.len);
      for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
	    // lse_refresh(lse, foo.E[i], foo.C, foo.rs);
	    // lse_vmat(lse);
	    // auto V = (double complex(*)[2 * ngauss + 2]) lse->VOME->data;
	    // ose00[i] = V[ngauss + xoffset][ngauss + yoffset];
	    // ose01[i] = V[ngauss + xoffset][2 * ngauss + 1 + yoffset];
	    // ose10[i] = V[2 * ngauss + 1 + xoffset][ngauss + yoffset];
	    // ose11[i] = V[2 * ngauss + 1 + xoffset][2 * ngauss + 1 + yoffset];
	    ose00[i] = OME_00(lse->ome, foo.E[i] + m11 + m12, 0.003525, 0.003525);
	    ose01[i] = OME_00(lse->ome, foo.E[i] + m11 + m12, 0.003525, 0.003525);
	    ose10[i] = OME_00(lse->ome, foo.E[i] + m11 + m12, 0.003525, 0.003525);
	    ose11[i] = OME_00(lse->ome, foo.E[i] + m11 + m12, 0.003525, 0.003525);
      }
      return 0;
}

int oTV(void *arg)
{
      argstruct foo = *(argstruct *)arg;
      LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon, 0);
      // auto buf = (double complex *)foo.res;
      auto buf = (double complex(*)[4][foo.size])foo.res;
      size_t pNgauss = foo.pNgauss;
      // printf("thread %lu, foo.res %p\n", foo.id, foo.res);
      // for (size_t i = 0; i < 11; i += 1) {
      //   printf("%.2f ", creal(buf[i]));
      // }
      // puts("");
      for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
	    lse_compute(lse, foo.E[i], foo.C, foo.rs, 0);
	    auto T = (double complex(*)[2 * pNgauss + 2]) lse->TOME->data;
	    auto V = (double complex(*)[2 * pNgauss + 2]) lse->VOME->data;
	    // printf("thread %lu, i %lu, det[i] %.1f\n", foo.id, i,
	    // creal(det[i])); puts("");
	    buf[0][0][i] = T[pNgauss][pNgauss];
	    buf[0][1][i] = T[pNgauss][2 * pNgauss + 1];
	    buf[0][2][i] = T[2 * pNgauss + 1][pNgauss];
	    buf[0][3][i] = T[2 * pNgauss + 1][2 * pNgauss + 1];
	    buf[1][0][i] = V[pNgauss][pNgauss];
	    buf[1][1][i] = V[pNgauss][2 * pNgauss + 1];
	    buf[1][2][i] = V[2 * pNgauss + 1][pNgauss];
	    buf[1][3][i] = V[2 * pNgauss + 1][2 * pNgauss + 1];
      }
      return 0;
}

int det(void *arg)
{
      argstruct foo = *(argstruct *)arg;
      LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon, 0);
      double complex *res = (double complex *)foo.res;
      for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
	    res[i] = lse_detImVG(lse, foo.E[i], foo.C, PP, 0);
	    // lse_compute(lse, foo.E[i], foo.rs);
	    // res[i] = lse->det;
      }
      return 0;
}

int ma(void *arg)
{
      auto foo = *(argstruct *)arg;
      LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon, 0);
      double complex *res = (double complex *)foo.res;
      for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
	    res[i] = lse_detImVG(lse, foo.C[0], foo.C, PP, foo.E[i]);
      }
      return 0;
}

int both(void *arg)
{
      argstruct foo = *(argstruct *)arg;
      LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon, 0);
      auto buf = (double complex *)foo.res;
      size_t pNgauss = foo.pNgauss;
      // printf("thread %lu, foo.res %p\n", foo.id, foo.res);
      // for (size_t i = 0; i < 11; i += 1) {
      //   printf("%.2f ", creal(buf[i]));
      // }
      // puts("");
      for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
	    lse_compute(lse, foo.E[i], foo.C, foo.rs, 0);
	    auto res = buf + 5 * i;
	    auto T = (double complex(*)[2 * pNgauss + 2]) lse->TOME->data;
	    // printf("thread %lu, i %lu, det[i] %.1f\n", foo.id, i,
	    // creal(det[i])); puts("");
	    res[0] = lse->det;
	    res[1] = T[pNgauss][pNgauss];
	    res[2] = T[pNgauss][2 * pNgauss + 1];
	    res[3] = T[2 * pNgauss + 1][pNgauss];
	    res[4] = T[2 * pNgauss + 1][2 * pNgauss + 1];
      }
      return 0;
}

int poles(void *arg)
{
      auto foo = *(struct polestruct *)arg;
      size_t len = foo.ilen * foo.rlen;
      LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon, 0);
      double complex *res = foo.res;
      for (size_t idx = foo.start; idx < foo.start + foo.len; idx += 1) {
	    size_t re = idx / foo.ilen;
	    size_t im = idx % foo.ilen;
	    double complex E = foo.Er[re] + foo.Ei[im] * I;
	    // res[idx] = pole(lse, E, foo.C, NN);
	    // res[idx + len] = pole(lse, E, foo.C, PN);
	    // res[idx + len * 2] = pole(lse, E, foo.C, NP);
	    res[idx + len * 3] = pole(lse, E, foo.C, PP);
      }
      return 0;
}

int cst(void *arg)
{
      auto foo = *(argstruct *)arg;
      LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon, 0);
      double *res = foo.res;
      double C[4] = { 0 };
      C[0] = -1.010589943548671;
      C[2] = -1.220749787118462;
      for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
	    C[3] = foo.E[i];
	    res[i] = lse_cost(lse, C, PP);
      }
      return 0;
}

int trG(void *arg)
{
      auto foo = *(argstruct *)arg;
      LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon, 0);
      auto res = (double complex(*)[foo.size])foo.res;
      auto pNgauss = foo.pNgauss;
      for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
	    lse_refresh(lse, foo.E[i], foo.C, foo.rs, 0);
	    lse_gmat(lse);
	    auto G = (double complex(*)[2 * pNgauss + 2]) lse->G->data;
	    for (size_t ch = 0; ch < 2; ch += 1) {
		  res[ch][i] = 0;
		  for (size_t p = 0; p < pNgauss + 1; p += 1) {
			res[ch][i] += G[p + ch * (pNgauss + 1)][p + ch * (pNgauss + 1)];
		  }
	    }
      }
      return 0;
}

void Free(void *ptr) { free(ptr); }

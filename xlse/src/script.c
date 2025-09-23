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

double complex *onshell(double *E, size_t len, size_t pNgauss, double Lambda,
                        double epsilon) {
    double complex(*res)[len] = malloc(sizeof(double complex) * len * 2);
    LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(pNgauss, Lambda, epsilon);
    for (size_t i = 0; i < len; i += 1) {
        lse_refresh(lse, E[i], (double[4]){0, 0, 0, 0}, PP);
        res[0][i] = lse->x0[0];
        res[1][i] = lse->x0[1];
    }
    return (double complex *)res;
}
double complex *onshellT(double *E, size_t len, double C[4], size_t pNgauss,
                         double Lambda, double epsilon) {
    // puts("in C");
    // for (size_t i = 0; i < len; i += 1) {
    //   printf("%f,  ", E[i]);
    // }
    // puts("");
    double complex *res =
        malloc(sizeof(double complex) * len * NCHANNELS * NCHANNELS);
    thrd_t tid[NTHREADS];
    argstruct args[NTHREADS];
    size_t ntasks = len / NTHREADS;
    size_t residue = len % NTHREADS;
    for (size_t i = 0; i < NTHREADS; i += 1) {
        args[i].size = len;
        args[i].pNgauss = pNgauss;
        args[i].Lambda = Lambda;
        args[i].epsilon = epsilon;
        args[i].res = res;
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

double complex *onshellT_single(double *E, size_t len, double C[4],
                                size_t pNgauss, double Lambda, double epsilon) {
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
        thrd_create(&tid[i], oTsing, &args[i]);
    }
    for (size_t i = 0; i < NTHREADS; i += 1) {
        thrd_join(tid[i], NULL);
    }
    return res;
}

double complex *onshellG(double *E, size_t len, double C[4], size_t pNgauss,
                         double Lambda, double epsilon) {
    double complex *res =
        malloc(sizeof(double complex) * len * NCHANNELS * NCHANNELS);
    thrd_t tid[NTHREADS];
    argstruct args[NTHREADS];
    size_t ntasks = len / NTHREADS;
    size_t residue = len % NTHREADS;
    for (size_t i = 0; i < NTHREADS; i += 1) {
        args[i].pNgauss = pNgauss;
        args[i].Lambda = Lambda;
        args[i].epsilon = epsilon;
        args[i].size = len;
        args[i].res = res;
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

double complex *onshellV(double *E, size_t len, double C[4], size_t pNgauss,
                         double Lambda, double epsilon) {
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
        thrd_create(&tid[i], oV, &args[i]);
    }
    for (size_t i = 0; i < NTHREADS; i += 1) {
        thrd_join(tid[i], NULL);
    }
    return res;
}

double complex *onshellTV(double *E, size_t len, double C[4], size_t pNgauss,
                          double Lambda, double epsilon) {
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

struct OME *ome_malloc() {
    struct OME *ome = malloc(sizeof(*ome));
    ome_build(ome);
    return ome;
}

double complex V(double E, double complex p, double complex pprime) {
    return ANA_00(E, p, pprime, m_pi);
}

double complex Vquad(double E, double complex p, double complex pprime) {
    static struct OME ome = {0};
    ome_build(&ome);
    return Vpiu(ome, E, p, pprime, m_B_star, m_B, m_B_star, m_B, m_pi, 1);
}

double complex *V3d(double E, size_t pNgauss, double Lambda, double epsilon) {
    LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(pNgauss, Lambda, epsilon);
    lse_refresh(lse, E, (double[4]){0, 0, 0, 0}, PP);
    lse_vmat(lse);
    double complex(*vmat)[2 * pNgauss + 2] =
        malloc(sizeof(double complex) * 2 * (pNgauss + 1) * 2 * (pNgauss + 1));
    for (size_t i = 0; i < pNgauss + 1; i += 1) {
        for (size_t j = 0; j < pNgauss + 1; j += 1) {
            vmat[i][j] = matrix_get(lse->VOME, i, j);
        }
    }
    return vmat;
    // return (double complex*) lse->VOME->data;
}

double complex *traceG(double *E, size_t len, double C[4], size_t pNgauss,
                       double Lambda, double epsilon) {
    double complex *res = malloc(sizeof(double complex) * len * NCHANNELS);
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

double complex *Det(double *E, size_t len, double C[4], size_t pNgauss,
                    double Lambda, double epsilon) {
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

double complex *Both(double *E, size_t len, double C[4], size_t pNgauss,
                     double Lambda, double epsilon) {
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

void Evec(double *E, double V0) {
    WaveFunction *wf [[gnu::cleanup(wffree)]] =
        WFnew(partialwave, RLAMBDA, RNGAUSS);
    // LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(10, 4, 1e-6);
    for (size_t i = 0; i < N_TOWER; i += 1) {
        E[i] = wf->E_solution->data[i] + V0;
    }
}

void Ntower(size_t *len) { *len = N_TOWER; }

double complex *Poles(double *Er, size_t rlen, double *Ei, size_t ilen,
                      double C[4], size_t pNgauss, double Lambda,
                      double epsilon) {
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

double *Fit(double C[4], size_t pNgauss, double Lambda, double epsilon) {
    LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(pNgauss, Lambda, epsilon);
    return minimize(lse, C);
}

double *Cost(double *c, size_t len, double complex resonance, size_t pNgauss,
             double Lambda, double epsilon) {
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

int oT(void *arg) {
    argstruct foo = *(argstruct *)arg;
    LSE *lse [[gnu::cleanup(lsefree)]] =
        lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon);
    size_t ngauss = foo.pNgauss;
    double complex(*ose)[NCHANNELS][foo.size] = foo.res;
    // printf("start: %lu, len: %lu\n", foo.start, foo.len);
    for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
        lse_compute(lse, foo.E[i], foo.C, foo.rs);
// printf("%f\n", foo.E[i]);
#ifdef TPO
        lse_X(lse);
        lse_XtX(lse);
        ose00[i] = lse->onshellT[0][0];
        ose01[i] = lse->onshellT[0][1];
        ose10[i] = lse->onshellT[1][0];
        ose11[i] = lse->onshellT[1][1];
#else
        // auto T = (double complex(*)[2 * ngauss + 2]) lse->TOME->data;
        for (size_t alpha = 0; alpha < NCHANNELS; alpha += 1) {
            for (size_t beta = 0; beta < NCHANNELS; beta += 1) {
                ose[alpha][beta][i] =
                    matrix_get(lse->TOME, ngauss + alpha * (ngauss + 1),
                               ngauss + beta * (ngauss + 1));
            }
        }
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

int oTsing(void *arg) {
    argstruct foo = *(argstruct *)arg;
    LSE *lse [[gnu::cleanup(lsefree)]] =
        lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon);
    size_t ngauss = foo.pNgauss;
    double complex *res = foo.res;
    // printf("start: %lu, len: %lu\n", foo.start, foo.len);
    for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
        lse_compute_single(lse, foo.E[i], foo.C, foo.rs);
// printf("%f\n", foo.E[i]);
#ifdef TPO
        lse_X(lse);
        lse_XtX(lse);
        ose00[i] = lse->onshellT[0][0];
        ose01[i] = lse->onshellT[0][1];
        ose10[i] = lse->onshellT[1][0];
        ose11[i] = lse->onshellT[1][1];
#else
        res[i] = matrix_get(lse->TOME, ngauss, ngauss);
        // res[i] = lse->onshellT[0][0];
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

int oG(void *arg) {
    argstruct foo = *(argstruct *)arg;
    LSE *lse [[gnu::cleanup(lsefree)]] =
        lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon);
    size_t ngauss = foo.pNgauss;
    double complex(*ose)[NCHANNELS][foo.size] = foo.res;
    // printf("start: %lu, len: %lu\n", foo.start, foo.len);
    for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
        lse_refresh(lse, foo.E[i], foo.C, foo.rs);
        lse_gmat(lse);
        auto G = (double complex(*)[NCHANNELS * (ngauss + 1)]) lse->G->data;
        for (size_t alpha = 0; alpha < NCHANNELS; alpha += 1) {
            for (size_t beta = 0; beta < NCHANNELS; beta += 1) {
                ose[alpha][beta][i] = G[ngauss + alpha * (ngauss + 1)]
                                       [ngauss + beta * (ngauss + 1)];
            }
        }
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

int oV(void *arg) {
    argstruct foo = *(argstruct *)arg;
    LSE *lse [[gnu::cleanup(lsefree)]] =
        lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon);
    size_t ngauss = foo.pNgauss;
    onshellElements *res = (onshellElements *)foo.res;
    double complex *ose00 = (double complex *)res->ose00;
    double complex *ose01 = (double complex *)res->ose01;
    double complex *ose10 = (double complex *)res->ose10;
    double complex *ose11 = (double complex *)res->ose11;
    int64_t xoffset = 0;
    int64_t yoffset = 0;
    // printf("start: %lu, len: %lu\n", foo.start, foo.len);
    for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
        lse_refresh(lse, foo.E[i], foo.C, foo.rs);
        lse_vmat(lse);
        ose00[i] = matrix_get(lse->VOME, ngauss + xoffset, ngauss + yoffset);
        ose01[i] =
            matrix_get(lse->VOME, ngauss + xoffset, 2 * ngauss + 1 + yoffset);
        ose10[i] =
            matrix_get(lse->VOME, 2 * ngauss + 1 + xoffset, ngauss + yoffset);
        ose11[i] = matrix_get(lse->VOME, 2 * ngauss + 1 + xoffset,
                              2 * ngauss + 1 + yoffset);
        // ose00[i] = OME_00(lse->ome, foo.E[i] + m11 + m12, 0.003525,
        // 0.003525); ose01[i] = OME_00(lse->ome, foo.E[i] + m11 + m12,
        // 0.003525, 0.003525); ose10[i] = OME_00(lse->ome, foo.E[i] + m11 +
        // m12, 0.003525, 0.003525); ose11[i] = OME_00(lse->ome, foo.E[i] + m11
        // + m12, 0.003525, 0.003525); ose00[i] = OMEANA_00(foo.E[i] + m11 + m12
        // , lse->x0[0], lse->x0[0]);
    }
    return 0;
}

int oTV(void *arg) {
    argstruct foo = *(argstruct *)arg;
    LSE *lse [[gnu::cleanup(lsefree)]] =
        lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon);
    // auto buf = (double complex *)foo.res;
    auto buf = (double complex(*)[4][foo.size])foo.res;
    size_t pNgauss = foo.pNgauss;
    // printf("thread %lu, foo.res %p\n", foo.id, foo.res);
    // for (size_t i = 0; i < 11; i += 1) {
    //   printf("%.2f ", creal(buf[i]));
    // }
    // puts("");
    for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
        lse_compute(lse, foo.E[i], foo.C, foo.rs);
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

int det(void *arg) {
    argstruct foo = *(argstruct *)arg;
    LSE *lse [[gnu::cleanup(lsefree)]] =
        lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon);
    double complex *res = (double complex *)foo.res;
    for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
        // res[i] = lse_detImVG(lse, foo.E[i], foo.C, PP);
        lse_compute(lse, foo.E[i], foo.C, foo.rs);
        res[i] = lse->det;
    }
    return 0;
}

int both(void *arg) {
    argstruct foo = *(argstruct *)arg;
    LSE *lse [[gnu::cleanup(lsefree)]] =
        lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon);
    auto buf = (double complex *)foo.res;
    size_t pNgauss = foo.pNgauss;
    // printf("thread %lu, foo.res %p\n", foo.id, foo.res);
    // for (size_t i = 0; i < 11; i += 1) {
    //   printf("%.2f ", creal(buf[i]));
    // }
    // puts("");
    for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
        lse_compute(lse, foo.E[i], foo.C, foo.rs);
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

int poles(void *arg) {
    auto foo = *(struct polestruct *)arg;
    size_t len = foo.ilen * foo.rlen;
    LSE *lse [[gnu::cleanup(lsefree)]] =
        lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon);
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

int cst(void *arg) {
    auto foo = *(argstruct *)arg;
    LSE *lse [[gnu::cleanup(lsefree)]] =
        lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon);
    double *res = foo.res;
    double C[4] = {0};
    C[0] = -1.010589943548671;
    C[2] = -1.220749787118462;
    for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
        C[3] = foo.E[i];
        res[i] = lse_cost(lse, C, PP);
    }
    return 0;
}

int trG(void *arg) {
    auto foo = *(argstruct *)arg;
    LSE *lse [[gnu::cleanup(lsefree)]] =
        lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon);
    auto res = (double complex(*)[foo.size])foo.res;
    auto pNgauss = foo.pNgauss;
    for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
        lse_refresh(lse, foo.E[i], foo.C, foo.rs);
        lse_gmat(lse);
        auto G = (double complex(*)[NCHANNELS * (pNgauss + 1)]) lse->G->data;
        for (size_t ch = 0; ch < NCHANNELS; ch += 1) {
            res[ch][i] = 0;
            for (size_t p = 0; p < pNgauss + 1; p += 1) {
                res[ch][i] += G[p + ch * (pNgauss + 1)][p + ch * (pNgauss + 1)];
            }
        }
    }
    return 0;
}

void Free(void *ptr) { free(ptr); }
double complex *getV(double E, size_t pNgauss, double Lambda, double epsilon) {
    // puts("getV get called");
    // printf("E: %f, pNgauss: %zu, Lambda: %f, epsilon: %f\n", E, pNgauss,
    // Lambda, epsilon);
    LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(pNgauss, Lambda, epsilon);
    lse_refresh(lse, E, (double[4]){0, 0, 0, 0}, PP);
    lse_vmat(lse);
    double complex(*V)[2 * pNgauss + 2] =
        malloc(sizeof(*V) * (2 * pNgauss + 2));
    for (size_t i = 0; i < 2 * pNgauss + 2; i += 1) {
        for (size_t j = 0; j < 2 * pNgauss + 2; j += 1) {
            V[i][j] = matrix_get(lse->VOME, i, j);
        }
    }
    return (double complex *)V;
}

double complex *getIntegrand() {
    double Lambda = 2;
    size_t Ngauss = 64;
    LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(Ngauss, Lambda, 1e-6);
    double E = 0.1;
    lse_refresh(lse, E - m11 - m12 + 2 * m_B + m_pi, (double[4]){0, 0, 0, 0},
                PP);
    double complex p1;
    double complex p2;
    p1 = lse->xi[12];
    p2 = lse->xi[1];
    double complex *res = malloc(sizeof(double complex) * (2 * DIMIM + DIMRE));
    for (size_t i = 0; i < 2 * DIMIM + DIMRE; i += 1) {
        res[i] = TOPTintegrand(E, lse->ome.xxpiup[i], p1, p2, m_B, m_B, m_pi);
    }
    return res;
}

double complex Gamma(double complex E, double complex p) {
    return Gamma0(E, p);
}

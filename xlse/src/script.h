#ifndef SCRIPT_H
#define SCRIPT_H


#include "complex.h"
#include "lse.h"
#include <stddef.h>
#include <threads.h>

typedef struct {
      void *ose00; // Using void* to avoid complex.h in C++ code
      void *ose01;
      void *ose10;
      void *ose11;
} onshellElements;

// Declare the function with void* to avoid complex.h dependency in C++
double complex *onshell(double *E, size_t len, size_t pNgauss, double Lambda, double epsilon);
double complex *onshellT(double *E, size_t len, double C[4], size_t pNgauss, double Lambda, double epsilon);
double complex *onshellT_single(double *E, size_t len, double C[4], size_t pNgauss, double Lambda, double epsilon);
double complex *onshellG(double *E, size_t len, double C[4], size_t pNgauss, double Lambda, double epsilon);
double complex *onshellV(double *E, size_t len, double C[4], size_t pNgauss, double Lambda, double epsilon);
double complex *onshellTV(double *E, size_t len, double C[4], size_t pNgauss, double Lambda, double epsilon);
double complex *traceG(double *E, size_t len, double C[4], size_t pNgauss, double Lambda, double epsilon);
double complex *Det(double *E, size_t len, double C[4], size_t pNgauss, double Lambda, double epsilon);
double complex *Both(double *E, size_t len, double C[4], size_t pNgauss, double Lambda, double epsilon);
double complex V(double E, double complex p, double complex pprime);
double complex *V3d(double E, size_t pNgauss, double Lambda, double epsilon);
double complex *getV(double E, size_t pNgauss, double Lambda, double epsilon);
double complex *getIntegrand();
double complex Gamma(double complex E, double complex p);
struct OME* ome_malloc();
void Evec(double *, double);
void Ntower(size_t *);
double complex *Poles(double *Er, size_t rlen, double *Ei, size_t ilen, double C[4], size_t pNgauss, double Lambda,
		      double epsilon);
double *Fit(double C[4], size_t pNgauss, double Lambda, double epsilon);
double *Cost(double *C, size_t len, double complex resonance, size_t pNgauss, double Lambda, double epsilon);
void Free(void *ptr);

void ose_free(onshellElements ose);
typedef struct {
      size_t pNgauss;
      double Lambda;
      double epsilon;
      size_t start;
      size_t len;
      double *E;
      double C[4];
      void *res;
      RS rs;
      size_t id;
      size_t size;
} argstruct;

struct polestruct {
      size_t pNgauss;
      double Lambda;
      double epsilon;
      size_t start;
      size_t len;
      double *Er;
      size_t rlen;
      double *Ei;
      size_t ilen;
      double C[4];
      void *res;
      RS rs;
      size_t id;
};

int oT(void *arg);
int oTsing(void *arg);
int oG(void *arg);
int oV(void *arg);
int oTV(void *arg);
int trG(void *arg);
int det(void *arg);
int both(void *arg);
int poles(void *arg);
int cst(void *arg);

#endif // !SCRIPT_H

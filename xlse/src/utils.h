#ifndef UTILS_H
#define UTILS_H

#include <complex.h>
#include <math.h>

char *formatC(double complex);
static inline double complex xsqrtleft(double complex x) { return csqrt(x); }

static inline double complex xsqrtright(double complex x) {
    return cimag(x) >= 0 ? csqrt(x + 0I) : -csqrt(x - 0I);
}

static inline double complex xsqrtup(double complex x) {
    return creal(x) < 0 && cimag >= 0 ? -csqrt(x + 0I) : csqrt(x);
}

static inline double complex xsqrtdown(double complex x) {
    return cimag(x) < 0 && creal(x) < 0 ? -csqrt(x) : csqrt(x + 0I);
}
static inline double complex sqrtdn(double complex x) {
    return csqrt(x * cexp(-I * M_PI / 2)) / csqrt(cexp(-I * M_PI / 2));
}

static inline double complex xlog(double complex x) {
    // return cimag(x) >= 0 ? clog(x + 0I) : clog(x + 0I) + 2 * M_PI * I;
    return fabs(cimag(x)) <= 1e-7 ? clog(creal(x) + 0I) : clog(x);
    return clog(x);
    auto res = clog(x);
    // return res;
    // return clog(cabs(x));
    if (cimag(res) < 0 && fabs(cimag(x)) < 1e-7 && creal(x) < 0) {
        res += 2 * M_PI * I;
    }
    return res;
}

static inline double complex xcube(double complex x) { return x * x * x; }
#endif // UTILS_H

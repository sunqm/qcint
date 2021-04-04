// This function is implemented based on libslater library
// https://github.com/nubakery/libslater

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "simd.h"
#include "rys_roots.h"
#include "erfc_roots_xw.dat"

void _CINT_clenshaw_d1(double *rr, const double *x, double u, int nroots);
void _CINT_clenshaw_dc(double *rr, const double *x, double u, int nroots);
void _CINT_matmul_14_14(double *imc, double *im, int nroots);

void CINTsr_rys_polyfits(int nroots, double x, double lower, double* u, double* w)
{
        const double* dx = DATA_X + (nroots-1)*nroots/2 * 1960;
        const double* dw = DATA_W + (nroots-1)*nroots/2 * 1960;
        double ll, t, tt;
        int k, it;
        int offset;
        double im[14*nroots+10];

        t = x;
        if (t > 19682.99) t = 19682.99;
        if (t > 1.0) {
                tt = log(t) * 0.9102392266268373 + 1.0; // log(3)+1 
        } else {
                tt = sqrt(t);
        }

        it = (int)tt;
        tt = tt - it;
        tt = 2.0 * tt - 1.0;
        ll = 2.0 * lower - 1.0;

        offset = nroots * 196 * it;
        _CINT_clenshaw_dc(im, dx+offset, ll, nroots);
        _CINT_clenshaw_d1(u, im, tt, nroots);
        _CINT_clenshaw_dc(im, dw+offset, ll, nroots);
        _CINT_clenshaw_d1(w, im, tt, nroots);

        double fac = exp(-x * lower * lower);
        for (k = 0; k < nroots; k++) {
                u[k*SIMDD] = u[k*SIMDD] / (1 - u[k*SIMDD]);
                w[k*SIMDD] *= fac;
        }
}

/*
 * Qcint is a general GTO integral library for computational chemistry
 * Copyright (C) 2014- Qiming Sun <osirpt.sun@gmail.com>
 *
 * This file is part of Qcint.
 *
 * Qcint is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <math.h>
#include <complex.h>
#include "simd.h"

void CINTdcmplx_re(const int n, double complex *RESTRICT z,
                   const double *RESTRICT re)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = re[i] + 0 * _Complex_I;
        }
}

void CINTdcmplx_im(const int n, double complex *RESTRICT z,
                   const double *RESTRICT im)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = 0 + im[i] * _Complex_I;
        }
}

void CINTdcmplx_pp(const int n, double complex *RESTRICT z,
                   const double *RESTRICT re, const double *RESTRICT im)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = re[i] + im[i] * _Complex_I;
        }
}
void CINTdcmplx_pn(const int n, double complex *RESTRICT z,
                   const double *RESTRICT re, const double *RESTRICT im)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = re[i] - im[i] * _Complex_I;
        }
}
void CINTdcmplx_np(const int n, double complex *RESTRICT z,
                   const double *RESTRICT re, const double *RESTRICT im)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = -re[i] + im[i] * _Complex_I;
        }
}
void CINTdcmplx_nn(const int n, double complex *RESTRICT z,
                   const double *RESTRICT re, const double *RESTRICT im)
{
        int i;
        for (i = 0; i < n; i++) {
                z[i] = -re[i] - im[i] * _Complex_I;
        }
}


double CINTsquare_dist(const double *r1, const double *r2)
{
        double r12[3];

        r12[0] = r1[0] - r2[0];
        r12[1] = r1[1] - r2[1];
        r12[2] = r1[2] - r2[2];

        return r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2];
}

static int factorial(int n)
{
        int i, fact = 1;
        for (i = 1; i <= n; i++) {
                fact *= i;
        }
        return fact;
}

/*
 * Normalized factor for GTO radial part g=r^l e^{-\alpha r^2}
 *
 * \frac{1}{\sqrt{\int g^2 r^2 dr}}
 *   = \sqrt{\frac{2^{2l+3} (l+1)! (2a)^{l+1.5}}{(2l+2)!\sqrt{\pi}}}
 *
 * Ref: H. B. Schlegel and M. J. Frisch, Int. J. Quant.  Chem., 54(1995), 83-87.
 */
double CINTgto_norm(int n, double a)
{
        double nn = pow(2, (2*n+3)) * factorial(n+1) * pow((2*a), (n+1.5)) \
                / (factorial(2*n+2) * sqrt(M_PI));
        return sqrt(nn);
}

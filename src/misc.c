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
#include "config.h"
#include "cint_const.h"

#if defined(__GNUC__)
#define ALIGN16 __attribute__((aligned(16)))
#define RESTRICT __restrict__
#else
#define ALIGN16
#define RESTRICT
#endif

void CINTdcmplx_re(const FINT n,
                   double complex *RESTRICT z,
                   const double *RESTRICT re)
{
        FINT i;
        for (i = 0; i < n; i++) {
                z[i] = re[i] + 0 * _Complex_I;
        }
}

void CINTdcmplx_im(const FINT n, double complex *z, const double *im)
{
        FINT i;
        for (i = 0; i < n; i++) {
                z[i] = 0 + im[i] * _Complex_I;
        }
}

void CINTdcmplx_pp(const FINT n, double complex *z,
                   const double *re, const double *im)
{
        FINT i;
        for (i = 0; i < n; i++) {
                z[i] = re[i] + im[i] * _Complex_I;
        }
}
void CINTdcmplx_pn(const FINT n, double complex *z,
                   const double *re, const double *im)
{
        FINT i;
        for (i = 0; i < n; i++) {
                z[i] = re[i] - im[i] * _Complex_I;
        }
}
void CINTdcmplx_np(const FINT n, double complex *z,
                   const double *re, const double *im)
{
        FINT i;
        for (i = 0; i < n; i++) {
                z[i] = -re[i] + im[i] * _Complex_I;
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

static FINT factorial(FINT n)
{
        FINT i, fact = 1;
        for (i = 1; i <= n; i++) {
                fact *= i;
        }
        return fact;
}
double CINTgto_norm(FINT n, double a)
{
        double nn = pow(2, (2*n+3)) * factorial(n+1) * pow((2*a), (n+1.5)) \
                / (factorial(2*n+2) * sqrt(M_PI));
        return sqrt(nn);
}
double CINTgto_norm_(FINT *n, double *a)
{
        return CINTgto_norm(*n, *a);
}

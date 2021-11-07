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

/*
 * blas-like functions
 */

#include <string.h>
#include <complex.h>
#include "fblas.h"

#define OF_CMPLX        2

void CINTdset0(int n, double *x)
{
        int i;
        for (i = 0; i < n; i++) {
                x[i] = 0;
        }
}


/*
 * v = a * x + y
 */
void CINTdaxpy2v(int n, double a, double *x, double *y, double *v)
{
        //cblas_dcopy(n, y, 1, v, 1);
        //cblas_daxpy(n, a, x, 1, v, 1);
        int i;
        for (i = 0; i < n; i++) {
                v[i] = a * x[i] + y[i];
        }
}

/*
 * a[m,n] -> a_t[n,m]
 */
void CINTdmat_transpose(double *a_t, double *a, int m, int n)
{
        int i, j, k;
        double *pa1, *pa2, *pa3;

        for (j = 0; j < n-3; j+=4) {
                pa1 = a_t + m;
                pa2 = pa1 + m;
                pa3 = pa2 + m;
                for (i = 0, k = j; i < m; i++, k+=n) {
                        a_t[i] = a[k+0];
                        pa1[i] = a[k+1];
                        pa2[i] = a[k+2];
                        pa3[i] = a[k+3];
                }
                a_t += m * 4;
        }

        switch (n-j) {
        case 1:
                for (i = 0, k = j; i < m; i++, k+=n) {
                        a_t[i] = a[k];
                }
                break;
        case 2:
                pa1 = a_t + m;
                for (i = 0, k = j; i < m; i++, k+=n) {
                        a_t[i] = a[k+0];
                        pa1[i] = a[k+1];
                }
                break;
        case 3:
                pa1 = a_t + m;
                pa2 = pa1 + m;
                for (i = 0, k = j; i < m; i++, k+=n) {
                        a_t[i] = a[k+0];
                        pa1[i] = a[k+1];
                        pa2[i] = a[k+2];
                }
                break;
        }
}

/*
 * a_t[n,m] += a[m,n]
 */
void CINTdplus_transpose(double *a_t, double *a, int m, int n)
{
        int i, j, k;

        for (j = 0; j < n-3; j+=4) {
#pragma GCC ivdep
                for (i = 0; i < m; i++) {
                        a_t[(j+0)*m+i] += a[i*n+j+0];
                        a_t[(j+1)*m+i] += a[i*n+j+1];
                        a_t[(j+2)*m+i] += a[i*n+j+2];
                        a_t[(j+3)*m+i] += a[i*n+j+3];
                }
        }

        switch (n-j) {
        case 1:
#pragma GCC ivdep
                for (i = 0; i < m; i++) {
                        a_t[j*m+i] += a[i*n+j];
                }
                break;
        case 2:
#pragma GCC ivdep
                for (i = 0; i < m; i++) {
                        a_t[(j+0)*m+i] += a[i*n+j+0];
                        a_t[(j+1)*m+i] += a[i*n+j+1];
                }
                break;
        case 3:
#pragma GCC ivdep
                for (i = 0; i < m; i++) {
                        a_t[(j+0)*m+i] += a[i*n+j+0];
                        a_t[(j+1)*m+i] += a[i*n+j+1];
                        a_t[(j+2)*m+i] += a[i*n+j+2];
                }
                break;
        }
}

/*
 * a[m,n] -> a_t[n,m]
 */
void CINTzmat_transpose(double complex *a_t, double complex *a, int m, int n)
{
        int i, j;

        switch (n) {
        case 2:
                for (i = 0; i < m; i++) {
                        a_t[i  ] = a[2*i+0];
                        a_t[i+m] = a[2*i+1];
                }
                break;
        default:
                switch (m) {
                case 2: for (i = 0; i < n; i++) {
                                a_t[2*i+0] = a[i  ];
                                a_t[2*i+1] = a[i+n];
                        }
                        break;
                default:
                        for (i = 0; i < n; i++) {
                                for (j = 0; j < m; j++) {
                                        a_t[i*m+j] = a[j*n+i];
                                }
                        }
                }
        }
}

void CINTzmat_dagger(double complex *a_t, double complex *a, int m, int n)
{
        int i, j;

        for (i = 0; i < n; i++) {
                for (j = 0; j < m; j++) {
                        a_t[i*m+j] = conj(a[j*n+i]);
                }
        }
}

void CINTdgemm_NN1(int m, int n, int k,
                   double *a, double *b, double *c, int ldc)
{
        int i, j, kp;
        double bi;
        for (j = 0; j < n; j++) {
                for (i = 0; i < m; i++) {
                        c[i+ldc*j] = 0;
                }
                for (kp = 0; kp < k; kp++) {
                        bi = b[kp+k*j];
#pragma GCC ivdep
                        for (i = 0; i < m; i++) {
                                c[i+ldc*j] += a[i+m*kp] * bi;
                        }
                }
        }
}

void CINTdgemm_NN(int m, int n, int k,
                  double *a, double *b, double *c)
{
        CINTdgemm_NN1(m, n, k, a, b, c, m);
}

void CINTdgemm_TN(int m, int n, int k,
                  double *a, double *b, double *c)
{
        int i, j, kp;
        double ci;
        for (j = 0; j < n; j++) {
                for (i = 0; i < m; i++) {
                        ci = 0;
#pragma GCC ivdep
                        for (kp = 0; kp < k; kp++) {
                                ci += a[kp+k*i] * b[kp+k*j];
                        }
                        c[i+m*j] = ci;
                }
        }
}

void CINTdgemm_NT(int m, int n, int k,
                  double *a, double *b, double *c)
{
        int i, j, kp;
        double bi;
        for (j = 0; j < n; j++) {
                for (i = 0; i < m; i++) {
                        c[i+m*j] = 0;
                }
                for (kp = 0; kp < k; kp++) {
                        bi = b[j+n*kp];
#pragma GCC ivdep
                        for (i = 0; i < m; i++) {
                                c[i+m*j] += a[i+m*kp] * bi;
                        }
                }
        }
}

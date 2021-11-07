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

#if defined __cplusplus
extern "C" {
#endif
#include <complex.h>
void CINTdset0(int n, double *x);
void CINTdaxpy2v(const int n, double a, double *x, double *y, double *v);
void CINTdmat_transpose(double *a_t, double *a, int m, int n);
void CINTdplus_transpose(double *a_t, double *a, int m, int n);
void CINTzmat_transpose(double complex *a_t, double complex *a, int m, int n);
void CINTzmat_dagger(double complex *a_c, double complex *a, int m, int n);

void CINTdgemm_NN(int m, int n, int k,
                  double *a, double *b, double *c);
void CINTdgemm_NN1(int m, int n, int k,
                   double *a, double *b, double *c, int ldc);
void CINTdgemm_TN(int m, int n, int k,
                  double *a, double *b, double *c);
void CINTdgemm_NT(int m, int n, int k,
                  double *a, double *b, double *c);
#if defined __cplusplus
} // end extern "C"
#endif

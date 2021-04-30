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

#include "config.h"
#include "g1e.h"

void CINTrys_roots(int nroots, double x, double *u, double *w);
void _CINTrys_roots_batch(int nroots, double *x, double *u, double *w, int count);
void CINTsr_rys_roots(int nroots, double x, double lower, double *u, double *w);
void CINTsr_rys_polyfits(int nroots, double x, double lower, double *u, double *w);
int _CINTsr_rys_roots_batch(CINTEnvVars *envs, double *x, double *theta, double *u, double *w, int count);

void CINTstg_roots(int nroots, double ta, double ua, double* rr, double* ww);
void _CINTstg_roots_batch(int nroots, double *ta, double *ua, double* rr, double* ww, int count);

int CINTrys_schmidt(int nroots, double x, double lower, double *roots, double *weights);
int CINTlrys_schmidt(int nroots, double x, double lower, double *roots, double *weights);
int CINTrys_laguerre(int n, double x, double lower, double *roots, double *weights);
int CINTlrys_laguerre(int n, double x, double lower, double *roots, double *weights);
int CINTrys_jacobi(int n, double x, double lower, double *roots, double *weights);
int CINTlrys_jacobi(int n, double x, double lower, double *roots, double *weights);
#ifdef HAVE_QUADMATH_H
int CINTqrys_schmidt(int nroots, double x, double lower, double *roots, double *weights);
int CINTqrys_laguerre(int n, double x, double lower, double *roots, double *weights);
int CINTqrys_jacobi(int n, double x, double lower, double *roots, double *weights);
#endif

void gamma_inc_like(double *f, double t, int m);
void lgamma_inc_like(long double *f, long double t, int m);
//void fmt1_gamma_inc_like(double *f, double t, int m);
//void fmt1_lgamma_inc_like(long double *f, long double t, int m);
void fmt_erfc_like(double *f, double t, double lower, int m);
void fmt1_erfc_like(double *f, double t, double lower, int m);
void fmt1_lerfc_like(long double *f, long double t, long double lower, int m);
#ifdef HAVE_QUADMATH_H
void qgamma_inc_like(__float128 *f, __float128 t, int m);
void fmt1_qerfc_like(__float128 *f, __float128 t, __float128 lower, int m);
#endif

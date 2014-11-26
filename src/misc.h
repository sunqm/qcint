/*
 * Qcint is a general GTO integral library for computational chemistry
 * Copyright (C) 2014 Qiming Sun
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

#include "fblas.h"

void CINTdcmplx_re(const int n, double complex *z, const double *re);
void CINTdcmplx_im(const int n, double complex *z, const double *im);
void CINTdcmplx_pp(const int n, double complex *z, const double *re, const double *im);
void CINTdcmplx_pn(const int n, double complex *z, const double *re, const double *im);
void CINTdcmplx_np(const int n, double complex *z, const double *re, const double *im);

double CINTsquare_dist(const double *r1, const double *r2);

void CINTrys_roots(const int nroots, double x, double *u, double *w);

double CINTgto_norm(int n, double a);


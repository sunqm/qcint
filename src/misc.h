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
#include "fblas.h"

void CINTdcmplx_re(const FINT n, double complex *z, const double *re);
void CINTdcmplx_im(const FINT n, double complex *z, const double *im);
void CINTdcmplx_pp(const FINT n, double complex *z, const double *re, const double *im);
void CINTdcmplx_pn(const FINT n, double complex *z, const double *re, const double *im);
void CINTdcmplx_np(const FINT n, double complex *z, const double *re, const double *im);

double CINTsquare_dist(const double *r1, const double *r2);

void CINTrys_roots(const FINT nroots, double x, double *u, double *w);

double CINTgto_norm(FINT n, double a);


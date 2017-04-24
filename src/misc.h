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

#include "fblas.h"

#define MIN(X,Y) ((X)<(Y)?(X):(Y))
#define MAX(X,Y) ((X)>(Y)?(X):(Y))
#define SQUARE(r)       ((r)[0]*(r)[0] + (r)[1]*(r)[1] + (r)[2]*(r)[2])

void CINTdcmplx_re(const int n, double complex *z, const double *re);
void CINTdcmplx_im(const int n, double complex *z, const double *im);
void CINTdcmplx_pp(const int n, double complex *z, const double *re, const double *im);
void CINTdcmplx_pn(const int n, double complex *z, const double *re, const double *im);
void CINTdcmplx_np(const int n, double complex *z, const double *re, const double *im);
void CINTdcmplx_nn(const int n, double complex *z, const double *re, const double *im);

double CINTsquare_dist(const double *r1, const double *r2);

double CINTgto_norm(int n, double a);


#define INT2CINT(NAME, SUFFIX) \
int c##NAME##SUFFIX(double *out, int *shls, int *atm, int natm, \
                    int *bas, int nbas, double *env, CINTOpt *opt) { \
        return NAME##SUFFIX(out, NULL, shls, atm, natm, bas, nbas, env, opt, NULL); \
} \
void c##NAME##SUFFIX##_optimizer(CINTOpt **opt, int *atm, int natm, \
                                 int *bas, int nbas, double *env) { \
        NAME##_optimizer(opt, atm, natm, bas, nbas, env); \
}

#define ALL_CINT(NAME) \
    INT2CINT(NAME, _cart) \
    INT2CINT(NAME, _sph) \
int c##NAME(double *out, int *shls, int *atm, int natm, \
            int *bas, int nbas, double *env, CINTOpt *opt) { \
        return NAME##_spinor((double complex *)out, NULL, shls, \
                             atm, natm, bas, nbas, env, opt, NULL); \
} \
void c##NAME##_optimizer(CINTOpt **opt, int *atm, int natm, \
                         int *bas, int nbas, double *env) { \
        NAME##_optimizer(opt, atm, natm, bas, nbas, env); \
}

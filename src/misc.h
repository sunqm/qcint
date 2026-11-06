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

#include "cint_config.h"
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

double CINTgto_norm(FINT n, double a);


#ifdef WITH_CINT2_INTERFACE
#define ALL_CINT(NAME) \
FINT c##NAME##_cart(double *out, FINT *shls, FINT *atm, FINT natm, \
            FINT *bas, FINT nbas, double *env, CINTOpt *opt) { \
        return NAME##_cart(out, NULL, shls, atm, natm, bas, nbas, env, opt, NULL); \
} \
void c##NAME##_cart_optimizer(CINTOpt **opt, FINT *atm, FINT natm, \
                         FINT *bas, FINT nbas, double *env) { \
        NAME##_optimizer(opt, atm, natm, bas, nbas, env); \
} \
FINT c##NAME##_sph(double *out, FINT *shls, FINT *atm, FINT natm, \
            FINT *bas, FINT nbas, double *env, CINTOpt *opt) { \
        return NAME##_sph(out, NULL, shls, atm, natm, bas, nbas, env, opt, NULL); \
} \
void c##NAME##_sph_optimizer(CINTOpt **opt, FINT *atm, FINT natm, \
                         FINT *bas, FINT nbas, double *env) { \
        NAME##_optimizer(opt, atm, natm, bas, nbas, env); \
} \
FINT c##NAME(double *out, FINT *shls, FINT *atm, FINT natm, \
            FINT *bas, FINT nbas, double *env, CINTOpt *opt) { \
        return NAME##_spinor((double complex *)out, NULL, shls, \
                             atm, natm, bas, nbas, env, opt, NULL); \
} \
void c##NAME##_optimizer(CINTOpt **opt, FINT *atm, FINT natm, \
                         FINT *bas, FINT nbas, double *env) { \
        NAME##_optimizer(opt, atm, natm, bas, nbas, env); \
}


#define ALL_CINT1E(NAME) \
FINT c##NAME##_cart(double *out, FINT *shls, FINT *atm, FINT natm, \
            FINT *bas, FINT nbas, double *env) { \
        return NAME##_cart(out, NULL, shls, atm, natm, bas, nbas, env, NULL, NULL); \
} \
FINT c##NAME##_sph(double *out, FINT *shls, FINT *atm, FINT natm, \
            FINT *bas, FINT nbas, double *env) { \
        return NAME##_sph(out, NULL, shls, atm, natm, bas, nbas, env, NULL, NULL); \
} \
FINT c##NAME(double *out, FINT *shls, FINT *atm, FINT natm, \
            FINT *bas, FINT nbas, double *env) { \
        return NAME##_spinor((double complex *)out, NULL, shls, \
                             atm, natm, bas, nbas, env, NULL, NULL); \
}

#else

#define ALL_CINT(NAME)
#define ALL_CINT1E(NAME)

#endif  // WITH_CINT2_INTERFACE

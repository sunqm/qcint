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

#define ALL_CINT_FORTRAN_(NAME) \
int c##NAME##_sph_(double *out, int *shls, int *atm, int *natm, \
                    int *bas, int *nbas, double *env, int64_t optptr_as_integer8) { \
        CINTOpt **opt = (CINTOpt **)optptr_as_integer8; \
        return NAME##_sph(out, NULL, shls, \
                          atm, *natm, bas, *nbas, env, *opt, NULL); \
} \
void c##NAME##_sph_optimizer_(CINTOpt **opt, int *atm, int *natm, \
                              int *bas, int *nbas, double *env) { \
        NAME##_optimizer(opt, atm, *natm, bas, *nbas, env); \
} \
int c##NAME##_cart_(double *out, int *shls, int *atm, int *natm, \
                    int *bas, int *nbas, double *env, int64_t optptr_as_integer8) { \
        CINTOpt **opt = (CINTOpt **)optptr_as_integer8; \
        return NAME##_cart(out, NULL, shls, \
                           atm, *natm, bas, *nbas, env, *opt, NULL); \
} \
void c##NAME##_cart_optimizer_(CINTOpt **opt, int *atm, int *natm, \
                               int *bas, int *nbas, double *env) { \
        NAME##_optimizer(opt, atm, *natm, bas, *nbas, env); \
} \
int c##NAME##_(double *out, int *shls, int *atm, int *natm, \
               int *bas, int *nbas, double *env, int64_t optptr_as_integer8) { \
        CINTOpt **opt = (CINTOpt **)optptr_as_integer8; \
        return NAME##_spinor((double complex *)out, NULL, shls, \
                             atm, *natm, bas, *nbas, env, *opt, NULL); \
} \
void c##NAME##_optimizer_(int64_t optptr_as_integer8, int *atm, int *natm, \
                         int *bas, int *nbas, double *env) { \
        CINTOpt **opt = (CINTOpt **)optptr_as_integer8; \
        NAME##_optimizer(opt, atm, *natm, bas, *nbas, env); \
}

#define ALL_CINT1E_FORTRAN_(NAME) \
int c##NAME##_sph_(double *out, int *shls, int *atm, int *natm, \
                    int *bas, int *nbas, double *env) { \
        return NAME##_sph(out, NULL, shls, atm, *natm, bas, *nbas, env, NULL, NULL); \
} \
int c##NAME##_cart_(double *out, int *shls, int *atm, int *natm, \
                    int *bas, int *nbas, double *env) { \
        return NAME##_cart(out, NULL, shls, \
                           atm, *natm, bas, *nbas, env, NULL, NULL); \
} \
int c##NAME##_(double *out, int *shls, int *atm, int *natm, \
               int *bas, int *nbas, double *env) { \
        return NAME##_spinor((double complex *)out, NULL, shls, \
                             atm, *natm, bas, *nbas, env, NULL, NULL); \
}

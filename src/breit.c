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
 * Breit = Gaunt + gauge
 * Gaunt ~ - \sigma1\dot\sigma2/r12
 * gauge ~  1/2 \sigma1\dot\sigma2/r12 - 1/2 (\sigma1\dot r12) (\sigma2\dot r12)/r12^3
 * Breit ~ -1/2 \sigma1\dot\sigma2/r12 - 1/2 (\sigma1\dot r12) (\sigma2\dot r12)/r12^3
 */

#include <stdlib.h>
#include <complex.h>
#include "cint_bas.h"
#include "optimizer.h"

#define DECLARE(X)      int X(double complex *opijkl, int *shls, \
                              int *atm, int natm, \
                              int *bas, int nbas, double *env, CINTOpt *opt)

#define BREIT0(X) \
DECLARE(cint2e_##X); \
DECLARE(cint2e_gauge_r1_##X); \
DECLARE(cint2e_gauge_r2_##X); \
void cint2e_breit_##X##_optimizer(CINTOpt **opt, int *atm, int natm, \
                                  int *bas, int nbas, double *env) \
{ \
        *opt = NULL; \
} \
int cint2e_breit_##X(double complex *opijkl, int *shls, \
                          int *atm, int natm, \
                          int *bas, int nbas, double *env, CINTOpt *opt) \
{ \
        int has_value = cint2e_##X(opijkl, shls, atm, natm, bas, nbas, env, NULL); \
 \
        const int ip = CINTcgto_spinor(shls[0], bas); \
        const int jp = CINTcgto_spinor(shls[1], bas); \
        const int kp = CINTcgto_spinor(shls[2], bas); \
        const int lp = CINTcgto_spinor(shls[3], bas); \
        const int nop = ip * jp * kp * lp; \
        double complex *buf = malloc(sizeof(double complex) * nop); \
        int i; \
        has_value = (cint2e_gauge_r1_##X(buf, shls, atm, natm, bas, nbas, env, NULL) || \
                     has_value); \
        /* [1/2 gaunt] - [1/2 xxx*\sigma\dot r1] */ \
        if (has_value) { \
                for (i = 0; i < nop; i++) { \
                        opijkl[i] = -opijkl[i] - buf[i]; \
                } \
        } \
        /* ... [- 1/2 xxx*\sigma\dot(-r2)] */ \
        has_value = (cint2e_gauge_r2_##X(buf, shls, atm, natm, bas, nbas, env, NULL) || \
                     has_value); \
        if (has_value) { \
                for (i = 0; i < nop; i++) { \
                        opijkl[i] = (opijkl[i] + buf[i]) * .5; \
                } \
        } \
        free(buf); \
        return has_value; \
}


BREIT0(ssp1ssp2);
BREIT0(ssp1sps2);
BREIT0(sps1ssp2);
BREIT0(sps1sps2);

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
#include "simd.h"
#include "cart2sph.h"
#include "g1e.h"
#include "g2e.h"
#include "optimizer.h"
#include "cint1e.h"
#include "cint2e.h"
#include "misc.h"
#include "c2f.h"

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

/* modfied from
 * '("int2e_breit_r1p2"  ( nabla \, r0 \| dot nabla-r12 \| \, nabla ))
 */
static void CINTgout2e_int2e_breit_r1p2(double *RESTRICT gout, double *RESTRICT g,
                                        int *RESTRICT idx, CINTEnvVars *envs) {
        int nf = envs->nf;
        int nfc = nf * 9;
        int nrys_roots = envs->nrys_roots;
        int ix, iy, iz, i, n;
        DECLARE_GOUT;
        double *RESTRICT g0 = g;
        double *RESTRICT g1 = g0 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g2 = g1 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g3 = g2 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g4 = g3 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g5 = g4 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g6 = g5 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g7 = g6 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g8 = g7 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g9 = g8 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g10 = g9 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g11 = g10 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g12 = g11 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g13 = g12 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g14 = g13 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g15 = g14 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g16 = g15 + envs->g_size * 3 * SIMDD;
        G2E_D_L(g1, g0, envs->i_l+2, envs->j_l+2, envs->k_l+0, envs->l_l+0);
        G2E_R0J(g3, g1, envs->i_l+1, envs->j_l+0, envs->k_l, envs->l_l);
        G2E_D_J(g4, g0, envs->i_l+1, envs->j_l+1, envs->k_l, envs->l_l);
        G2E_D_I(g5, g0, envs->i_l+1, envs->j_l+1, envs->k_l, envs->l_l);
        for (ix = 0; ix < envs->g_size * 3 * SIMDD; ix++) {g4[ix] += g5[ix];}
        G2E_D_J(g5, g1, envs->i_l+1, envs->j_l+1, envs->k_l, envs->l_l);
        G2E_D_I(g6, g1, envs->i_l+1, envs->j_l+1, envs->k_l, envs->l_l);
        for (ix = 0; ix < envs->g_size * 3 * SIMDD; ix++) {g5[ix] += g6[ix];}
        G2E_R0J(g7, g5, envs->i_l+1, envs->j_l+0, envs->k_l, envs->l_l);
        G2E_D_I(g12, g4, envs->i_l+0, envs->j_l, envs->k_l, envs->l_l);
        G2E_D_I(g15, g7, envs->i_l+0, envs->j_l, envs->k_l, envs->l_l);
        __MD r1;
        for (n = 0; n < nf; n++) {
                ix = idx[0+n*3];
                iy = idx[1+n*3];
                iz = idx[2+n*3];
                r1 = MM_SET1(0.);
                for (i = 0; i < nrys_roots; i++) {
                        r1 += MM_LOAD(g15+(ix+i)*SIMDD) * MM_LOAD(g0+(iy+i)*SIMDD) * MM_LOAD(g0+(iz+i)*SIMDD);
                        r1 += MM_LOAD(g12+(ix+i)*SIMDD) * MM_LOAD(g3+(iy+i)*SIMDD) * MM_LOAD(g0+(iz+i)*SIMDD);
                        r1 += MM_LOAD(g12+(ix+i)*SIMDD) * MM_LOAD(g0+(iy+i)*SIMDD) * MM_LOAD(g3+(iz+i)*SIMDD);
                        r1 += MM_LOAD(g3+(ix+i)*SIMDD) * MM_LOAD(g12+(iy+i)*SIMDD) * MM_LOAD(g0+(iz+i)*SIMDD);
                        r1 += MM_LOAD(g0+(ix+i)*SIMDD) * MM_LOAD(g15+(iy+i)*SIMDD) * MM_LOAD(g0+(iz+i)*SIMDD);
                        r1 += MM_LOAD(g0+(ix+i)*SIMDD) * MM_LOAD(g12+(iy+i)*SIMDD) * MM_LOAD(g3+(iz+i)*SIMDD);
                        r1 += MM_LOAD(g3+(ix+i)*SIMDD) * MM_LOAD(g0+(iy+i)*SIMDD) * MM_LOAD(g12+(iz+i)*SIMDD);
                        r1 += MM_LOAD(g0+(ix+i)*SIMDD) * MM_LOAD(g3+(iy+i)*SIMDD) * MM_LOAD(g12+(iz+i)*SIMDD);
                        r1 += MM_LOAD(g0+(ix+i)*SIMDD) * MM_LOAD(g0+(iy+i)*SIMDD) * MM_LOAD(g15+(iz+i)*SIMDD);
                }
                GOUT_SCATTER(gout, n, r1);
        }
}
static void CINTgout2e_int2e_breit_r1p2_simd1(double *RESTRICT gout, double *RESTRICT g,
                                              int *RESTRICT idx, CINTEnvVars *envs) {
        int nf = envs->nf;
        int nrys_roots = envs->nrys_roots;
        int ix, iy, iz, i, n;
        double *RESTRICT g0 = g;
        double *RESTRICT g1 = g0 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g2 = g1 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g3 = g2 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g4 = g3 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g5 = g4 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g6 = g5 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g7 = g6 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g8 = g7 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g9 = g8 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g10 = g9 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g11 = g10 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g12 = g11 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g13 = g12 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g14 = g13 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g15 = g14 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g16 = g15 + envs->g_size * 3 * SIMDD;
        G2E_D_L_SIMD1(g1, g0, envs->i_l+2, envs->j_l+2, envs->k_l+0, envs->l_l+0);
        G2E_R0J_SIMD1(g3, g1, envs->i_l+1, envs->j_l+0, envs->k_l, envs->l_l);
        G2E_D_J_SIMD1(g4, g0, envs->i_l+1, envs->j_l+1, envs->k_l, envs->l_l);
        G2E_D_I_SIMD1(g5, g0, envs->i_l+1, envs->j_l+1, envs->k_l, envs->l_l);
        for (ix = 0; ix < envs->g_size * 3 * SIMDD; ix++) {g4[ix] += g5[ix];}
        G2E_D_J_SIMD1(g5, g1, envs->i_l+1, envs->j_l+1, envs->k_l, envs->l_l);
        G2E_D_I_SIMD1(g6, g1, envs->i_l+1, envs->j_l+1, envs->k_l, envs->l_l);
        for (ix = 0; ix < envs->g_size * 3 * SIMDD; ix++) {g5[ix] += g6[ix];}
        G2E_R0J_SIMD1(g7, g5, envs->i_l+1, envs->j_l+0, envs->k_l, envs->l_l);
        G2E_D_I_SIMD1(g12, g4, envs->i_l+0, envs->j_l, envs->k_l, envs->l_l);
        G2E_D_I_SIMD1(g15, g7, envs->i_l+0, envs->j_l, envs->k_l, envs->l_l);
        double s;
        for (n = 0; n < nf; n++) {
                ix = idx[0+n*3];
                iy = idx[1+n*3];
                iz = idx[2+n*3];
                s = 0;
                for (i = 0; i < nrys_roots; i++) {
                        s += g15[ix+i] * g0[iy+i] * g0[iz+i];
                        s += g12[ix+i] * g3[iy+i] * g0[iz+i];
                        s += g12[ix+i] * g0[iy+i] * g3[iz+i];
                        s += g3[ix+i] * g12[iy+i] * g0[iz+i];
                        s += g0[ix+i] * g15[iy+i] * g0[iz+i];
                        s += g0[ix+i] * g12[iy+i] * g3[iz+i];
                        s += g3[ix+i] * g0[iy+i] * g12[iz+i];
                        s += g0[ix+i] * g3[iy+i] * g12[iz+i];
                        s += g0[ix+i] * g0[iy+i] * g15[iz+i];
                }
                gout[n] = s;
        }
}
void int2e_breit_r1p2_optimizer(CINTOpt **opt, int *atm, int natm, int *bas, int nbas, double *env) {
        int ng[] = {2, 2, 0, 1, 4, 1, 1, 9};
        CINTall_2e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
int int2e_breit_r1p2_cart(double *out, int *dims, int *shls,
                int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {2, 2, 0, 1, 4, 1, 1, 9};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int2e_breit_r1p2;
        envs.f_gout_simd1 = &CINTgout2e_int2e_breit_r1p2_simd1;
        return CINT2e_cart_drv(out, dims, &envs, opt, cache);
} // int2e_breit_r1p2_cart
int int2e_breit_r1p2_sph(double *out, int *dims, int *shls,
                int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {2, 2, 0, 1, 4, 1, 1, 9};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int2e_breit_r1p2;
        envs.f_gout_simd1 = &CINTgout2e_int2e_breit_r1p2_simd1;
        return CINT2e_spheric_drv(out, dims, &envs, opt, cache);
} // int2e_breit_r1p2_sph
int int2e_breit_r1p2_spinor(double complex *out, int *dims, int *shls,
                int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {2, 2, 0, 1, 4, 1, 1, 9};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int2e_breit_r1p2;
        envs.f_gout_simd1 = &CINTgout2e_int2e_breit_r1p2_simd1;
        return CINT2e_spinor_drv(out, dims, &envs, opt, cache, &c2s_sf_2e1i, &c2s_sf_2e2i);
} // int2e_breit_r1p2_spinor
ALL_CINT(int2e_breit_r1p2)
//ALL_CINT_FORTRAN_(cint2e_breit_r1p2)

/* modfied from
 * '("int2e_breit_r2p2"  ( nabla \, r0 \| dot nabla-r12 \| \, nabla ))
 */
static void CINTgout2e_int2e_breit_r2p2(double *RESTRICT gout, double *RESTRICT g,
                                        int *RESTRICT idx, CINTEnvVars *envs) {
        int nf = envs->nf;
        int nfc = nf * 1;
        int nrys_roots = envs->nrys_roots;
        int ix, iy, iz, i, n;
        DECLARE_GOUT;
        double *RESTRICT g0 = g;
        double *RESTRICT g1 = g0 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g2 = g1 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g3 = g2 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g4 = g3 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g5 = g4 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g6 = g5 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g7 = g6 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g8 = g7 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g9 = g8 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g10 = g9 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g11 = g10 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g12 = g11 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g13 = g12 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g14 = g13 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g15 = g14 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g16 = g15 + envs->g_size * 3 * SIMDD;
        G2E_R0L(g2, g0, envs->i_l+2, envs->j_l+1, envs->k_l+0, envs->l_l+1);
        G2E_D_L(g3, g2, envs->i_l+2, envs->j_l+1, envs->k_l+0, envs->l_l+0);
        G2E_D_J(g4, g0, envs->i_l+1, envs->j_l+0, envs->k_l, envs->l_l);
        G2E_D_I(g5, g0, envs->i_l+1, envs->j_l+0, envs->k_l, envs->l_l);
        for (ix = 0; ix < envs->g_size * 3 * SIMDD; ix++) {g4[ix] += g5[ix];}
        G2E_D_J(g7, g3, envs->i_l+1, envs->j_l+0, envs->k_l, envs->l_l);
        G2E_D_I(g8, g3, envs->i_l+1, envs->j_l+0, envs->k_l, envs->l_l);
        for (ix = 0; ix < envs->g_size * 3 * SIMDD; ix++) {g7[ix] += g8[ix];}
        G2E_D_I(g12, g4, envs->i_l+0, envs->j_l, envs->k_l, envs->l_l);
        G2E_D_I(g15, g7, envs->i_l+0, envs->j_l, envs->k_l, envs->l_l);
        __MD r1;
        for (n = 0; n < nf; n++) {
                ix = idx[0+n*3];
                iy = idx[1+n*3];
                iz = idx[2+n*3];
                r1 = MM_SET1(0.);
                for (i = 0; i < nrys_roots; i++) {
                        r1 += MM_LOAD(g15+(ix+i)*SIMDD) * MM_LOAD(g0+(iy+i)*SIMDD) * MM_LOAD(g0+(iz+i)*SIMDD);
                        r1 += MM_LOAD(g12+(ix+i)*SIMDD) * MM_LOAD(g3+(iy+i)*SIMDD) * MM_LOAD(g0+(iz+i)*SIMDD);
                        r1 += MM_LOAD(g12+(ix+i)*SIMDD) * MM_LOAD(g0+(iy+i)*SIMDD) * MM_LOAD(g3+(iz+i)*SIMDD);
                        r1 += MM_LOAD(g3+(ix+i)*SIMDD) * MM_LOAD(g12+(iy+i)*SIMDD) * MM_LOAD(g0+(iz+i)*SIMDD);
                        r1 += MM_LOAD(g0+(ix+i)*SIMDD) * MM_LOAD(g15+(iy+i)*SIMDD) * MM_LOAD(g0+(iz+i)*SIMDD);
                        r1 += MM_LOAD(g0+(ix+i)*SIMDD) * MM_LOAD(g12+(iy+i)*SIMDD) * MM_LOAD(g3+(iz+i)*SIMDD);
                        r1 += MM_LOAD(g3+(ix+i)*SIMDD) * MM_LOAD(g0+(iy+i)*SIMDD) * MM_LOAD(g12+(iz+i)*SIMDD);
                        r1 += MM_LOAD(g0+(ix+i)*SIMDD) * MM_LOAD(g3+(iy+i)*SIMDD) * MM_LOAD(g12+(iz+i)*SIMDD);
                        r1 += MM_LOAD(g0+(ix+i)*SIMDD) * MM_LOAD(g0+(iy+i)*SIMDD) * MM_LOAD(g15+(iz+i)*SIMDD);
                }
                GOUT_SCATTER(gout, n, r1);
        }
}
static void CINTgout2e_int2e_breit_r2p2_simd1(double *RESTRICT gout, double *RESTRICT g,
                                              int *RESTRICT idx, CINTEnvVars *envs) {
        int nf = envs->nf;
        int nrys_roots = envs->nrys_roots;
        int ix, iy, iz, i, n;
        double *RESTRICT g0 = g;
        double *RESTRICT g1 = g0 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g2 = g1 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g3 = g2 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g4 = g3 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g5 = g4 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g6 = g5 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g7 = g6 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g8 = g7 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g9 = g8 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g10 = g9 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g11 = g10 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g12 = g11 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g13 = g12 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g14 = g13 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g15 = g14 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g16 = g15 + envs->g_size * 3 * SIMDD;
        G2E_R0L_SIMD1(g2, g0, envs->i_l+2, envs->j_l+1, envs->k_l+0, envs->l_l+1);
        G2E_D_L_SIMD1(g3, g2, envs->i_l+2, envs->j_l+1, envs->k_l+0, envs->l_l+0);
        G2E_D_J_SIMD1(g4, g0, envs->i_l+1, envs->j_l+0, envs->k_l, envs->l_l);
        G2E_D_I_SIMD1(g5, g0, envs->i_l+1, envs->j_l+0, envs->k_l, envs->l_l);
        for (ix = 0; ix < envs->g_size * 3 * SIMDD; ix++) {g4[ix] += g5[ix];}
        G2E_D_J_SIMD1(g7, g3, envs->i_l+1, envs->j_l+0, envs->k_l, envs->l_l);
        G2E_D_I_SIMD1(g8, g3, envs->i_l+1, envs->j_l+0, envs->k_l, envs->l_l);
        for (ix = 0; ix < envs->g_size * 3 * SIMDD; ix++) {g7[ix] += g8[ix];}
        G2E_D_I_SIMD1(g12, g4, envs->i_l+0, envs->j_l, envs->k_l, envs->l_l);
        G2E_D_I_SIMD1(g15, g7, envs->i_l+0, envs->j_l, envs->k_l, envs->l_l);
        double s;
        for (n = 0; n < nf; n++) {
                ix = idx[0+n*3];
                iy = idx[1+n*3];
                iz = idx[2+n*3];
                s = 0;
                for (i = 0; i < nrys_roots; i++) {
                        s += g15[ix+i] * g0[iy+i] * g0[iz+i];
                        s += g12[ix+i] * g3[iy+i] * g0[iz+i];
                        s += g12[ix+i] * g0[iy+i] * g3[iz+i];
                        s += g3[ix+i] * g12[iy+i] * g0[iz+i];
                        s += g0[ix+i] * g15[iy+i] * g0[iz+i];
                        s += g0[ix+i] * g12[iy+i] * g3[iz+i];
                        s += g3[ix+i] * g0[iy+i] * g12[iz+i];
                        s += g0[ix+i] * g3[iy+i] * g12[iz+i];
                        s += g0[ix+i] * g0[iy+i] * g15[iz+i];
                }
                gout[n] = s;
        }
}
void int2e_breit_r2p2_optimizer(CINTOpt **opt, int *atm, int natm, int *bas, int nbas, double *env) {
        int ng[] = {2, 1, 0, 2, 4, 1, 1, 1};
        CINTall_2e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
int int2e_breit_r2p2_cart(double *out, int *dims, int *shls,
                          int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {2, 1, 0, 2, 4, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int2e_breit_r2p2;
        envs.f_gout_simd1 = &CINTgout2e_int2e_breit_r2p2_simd1;
        return CINT2e_cart_drv(out, dims, &envs, opt, cache);
} // int2e_breit_r2p2_cart
int int2e_breit_r2p2_sph(double *out, int *dims, int *shls,
                         int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {2, 1, 0, 2, 4, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int2e_breit_r2p2;
        envs.f_gout_simd1 = &CINTgout2e_int2e_breit_r2p2_simd1;
        return CINT2e_spheric_drv(out, dims, &envs, opt, cache);
} // int2e_breit_r2p2_sph
int int2e_breit_r2p2_spinor(double complex *out, int *dims, int *shls,
                            int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {2, 1, 0, 2, 4, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int2e_breit_r2p2;
        envs.f_gout_simd1 = &CINTgout2e_int2e_breit_r2p2_simd1;
        return CINT2e_spinor_drv(out, dims, &envs, opt, cache, &c2s_sf_2e1i, &c2s_sf_2e2i);
} // int2e_breit_r2p2_spinor
ALL_CINT(int2e_breit_r2p2)
//ALL_CINT_FORTRAN_(cint2e_breit_r2p2)

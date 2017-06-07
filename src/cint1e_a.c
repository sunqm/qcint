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

#include <stdlib.h>
#include "cint_bas.h"
#include "simd.h"
#include "cart2sph.h"
#include "g1e.h"
#include "optimizer.h"
#include "cint1e.h"
#include "misc.h"
#include "c2f.h"

/* based on
 * '("int1e_r2_origi" (r dot r \| )) 
 */
static void CINTgout1e_int1e_r2_origi(double *gout, double *g, int *idx, CINTEnvVars *envs, int count) {
        CINTg1e_ovlp(g, envs, count);
        int nf = envs->nf;
        int nfc = nf * 1;
        int ix, iy, iz, n;
        DECLARE_GOUT;
        double *RESTRICT g0 = g;
        double *RESTRICT g1 = g0  + envs->g_size * 3 * SIMDD;
        double *RESTRICT g2 = g1  + envs->g_size * 3 * SIMDD;
        double *RESTRICT g3 = g2  + envs->g_size * 3 * SIMDD;
        __MD r1;
        G1E_R_I(g1, g0, envs->i_l+1, envs->j_l, 0);
        G1E_R_I(g3, g1, envs->i_l+0, envs->j_l, 0);
        for (n = 0; n < nf; n++) {
                ix = idx[0+n*3];
                iy = idx[1+n*3];
                iz = idx[2+n*3];
                r1 = MM_LOAD(g3+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
                r1+= MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g3+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
                r1+= MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g3+iz*SIMDD);
                GOUT_SCATTER(gout, n, r1);
        }
}
void int1e_r2_origi_optimizer(CINTOpt **opt, int *atm, int natm, int *bas, int nbas, double *env) {
        int ng[] = {2, 0, 0, 0, 2, 1, 1, 1};
        CINTall_1e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
int int1e_r2_origi_cart(double *out, int *dims, int *shls,
                int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {2, 0, 0, 0, 2, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_int1e_r2_origi;
        return CINT1e_drv(out, dims, &envs, opt, cache, &c2s_cart_1e);
} // int1e_r2_origi_cart
int int1e_r2_origi_sph(double *out, int *dims, int *shls,
                int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {2, 0, 0, 0, 2, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_int1e_r2_origi;
        return CINT1e_drv(out, dims, &envs, opt, cache, &c2s_sph_1e);
} // int1e_r2_origi_sph
int int1e_r2_origi_spinor(double complex *out, int *dims, int *shls,
                int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {2, 0, 0, 0, 2, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_int1e_r2_origi;
        return CINT1e_spinor_drv(out, dims, &envs, opt, cache, &c2s_sf_1e);
} // int1e_r2_origi_spinor
ALL_CINT1E(int1e_r2_origi)
//ALL_CINT1E_FORTRAN_(cint1e_r2_origi)

/* based on
 * '("int1e_r4_origi" (r dot r r dot r \| )) 
 */
static void CINTgout1e_int1e_r4_origi(double *gout, double *g, int *idx, CINTEnvVars *envs, int count) {
        CINTg1e_ovlp(g, envs, count);
        int nf = envs->nf;
        int nfc = nf * 1;
        int ix, iy, iz, n;
        DECLARE_GOUT;
        double *RESTRICT g0 = g;
        double *RESTRICT g1 = g0  + envs->g_size * 3 * SIMDD;
        double *RESTRICT g2 = g1  + envs->g_size * 3 * SIMDD;
        double *RESTRICT g3 = g2  + envs->g_size * 3 * SIMDD;
        double *RESTRICT g4 = g3  + envs->g_size * 3 * SIMDD;
        double *RESTRICT g5 = g4  + envs->g_size * 3 * SIMDD;
        double *RESTRICT g6 = g5  + envs->g_size * 3 * SIMDD;
        double *RESTRICT g7 = g6  + envs->g_size * 3 * SIMDD;
        double *RESTRICT g8 = g7  + envs->g_size * 3 * SIMDD;
        double *RESTRICT g9 = g8  + envs->g_size * 3 * SIMDD;
        double *RESTRICT g10 = g9  + envs->g_size * 3 * SIMDD;
        double *RESTRICT g11 = g10  + envs->g_size * 3 * SIMDD;
        double *RESTRICT g12 = g11  + envs->g_size * 3 * SIMDD;
        double *RESTRICT g13 = g12  + envs->g_size * 3 * SIMDD;
        double *RESTRICT g14 = g13  + envs->g_size * 3 * SIMDD;
        double *RESTRICT g15 = g14  + envs->g_size * 3 * SIMDD;
        __MD r1;
        __MD r2 = MM_SET1(2.);
        G1E_R_I(g1, g0, envs->i_l+3, envs->j_l, 0);
        G1E_R_I(g3, g1, envs->i_l+2, envs->j_l, 0);
        G1E_R_I(g4, g0, envs->i_l+1, envs->j_l, 0);
        G1E_R_I(g7, g3, envs->i_l+1, envs->j_l, 0);
        G1E_R_I(g12, g4, envs->i_l+0, envs->j_l, 0);
        G1E_R_I(g15, g7, envs->i_l+0, envs->j_l, 0);
        for (n = 0; n < nf; n++) {
                ix = idx[0+n*3];
                iy = idx[1+n*3];
                iz = idx[2+n*3];
                r1 = MM_LOAD(g15+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
                r1+= MM_LOAD(g12+ix*SIMDD) * MM_LOAD(g3+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD) * r2;
                r1+= MM_LOAD(g12+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g3+iz*SIMDD) * r2;
                r1+= MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g15+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
                r1+= MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g12+iy*SIMDD) * MM_LOAD(g3+iz*SIMDD) * r2;
                r1+= MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g15+iz*SIMDD);
                GOUT_SCATTER(gout, n, r1);
        }
}
void int1e_r4_origi_optimizer(CINTOpt **opt, int *atm, int natm, int *bas, int nbas, double *env) {
        int ng[] = {4, 0, 0, 0, 4, 1, 1, 1};
        CINTall_1e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
int int1e_r4_origi_cart(double *out, int *dims, int *shls,
                int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {4, 0, 0, 0, 4, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_int1e_r4_origi;
        return CINT1e_drv(out, dims, &envs, opt, cache, &c2s_cart_1e);
} // int1e_r4_origi_cart
int int1e_r4_origi_sph(double *out, int *dims, int *shls,
                int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {4, 0, 0, 0, 4, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_int1e_r4_origi;
        return CINT1e_drv(out, dims, &envs, opt, cache, &c2s_sph_1e);
} // int1e_r4_origi_sph
int int1e_r4_origi_spinor(double complex *out, int *dims, int *shls,
                int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {4, 0, 0, 0, 4, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_int1e_r4_origi;
        return CINT1e_spinor_drv(out, dims, &envs, opt, cache, &c2s_sf_1e);
} // int1e_r4_origi_spinor
ALL_CINT1E(int1e_r4_origi)
//ALL_CINT1E_FORTRAN_(cint1e_r4_origi)

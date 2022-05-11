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
 * '("int3c1e_r2_origk" ( \, \, r dot r r dot r r dot r))
 */
static void CINTgout1e_int3c1e_r2_origk(double *RESTRICT gout,
                double *RESTRICT g, int *RESTRICT idx, CINTEnvVars *envs, int count) {
        CINTg3c1e_ovlp(g, envs, count);
        int nf = envs->nf;
        int nfc = nf * 1;
        int ix, iy, iz, n;
        DECLARE_GOUT;
        double *RESTRICT g0 = g;
        double *RESTRICT g1 = g0 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g2 = g1 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g3 = g2 + envs->g_size * 3 * SIMDD;
        __MD r1;
        G1E_R_K(g1, g0, envs->i_l+0, envs->j_l+0, envs->k_l+1);
        G1E_R_K(g3, g1, envs->i_l+0, envs->j_l+0, envs->k_l+0);
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
void int3c1e_r2_origk_optimizer(CINTOpt **opt, int *atm, int natm, int *bas, int nbas, double *env) {
        int ng[] = {0, 0, 2, 0, 2, 1, 1, 1};
        CINTall_3c1e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
CACHE_SIZE_T int3c1e_r2_origk_cart(double *out, int *dims, int *shls,
                int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {0, 0, 2, 0, 2, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_int3c1e_r2_origk;
        return CINT3c1e_drv(out, dims, &envs, opt, cache, &c2s_cart_3c1e);
} // int3c1e_r2_origk_cart
CACHE_SIZE_T int3c1e_r2_origk_sph(double *out, int *dims, int *shls,
                int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {0, 0, 2, 0, 2, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_int3c1e_r2_origk;
        return CINT3c1e_drv(out, dims, &envs, opt, cache, &c2s_sph_3c1e);
} // int3c1e_r2_origk_sph
CACHE_SIZE_T int3c1e_r2_origk_spinor(double complex *out, int *dims, int *shls,
                int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {0, 0, 2, 0, 2, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_int3c1e_r2_origk;
        return CINT3c1e_spinor_drv(out, dims, &envs, opt, cache, &c2s_sf_3c2e1);
} // int3c1e_r2_origk_spinor
ALL_CINT(int3c1e_r2_origk)
ALL_CINT_FORTRAN_(int3c1e_r2_origk)

/* based on
 * '("int3c1e_r4_origk" ( \, \, r dot r r dot r))
 */
static void CINTgout1e_int3c1e_r4_origk(double *RESTRICT gout,
                double *RESTRICT g, int *RESTRICT idx, CINTEnvVars *envs, int count) {
        CINTg3c1e_ovlp(g, envs, count);
        int nf = envs->nf;
        int nfc = nf * 1;
        int ix, iy, iz, n;
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
        __MD r1;
        __MD r2 = MM_SET1(2.);
        G1E_R_K(g1, g0, envs->i_l+0, envs->j_l+0, envs->k_l+3);
        G1E_R_K(g3, g1, envs->i_l+0, envs->j_l+0, envs->k_l+2);
        G1E_R_K(g4, g0, envs->i_l+0, envs->j_l+0, envs->k_l+1);
        G1E_R_K(g7, g3, envs->i_l+0, envs->j_l+0, envs->k_l+1);
        G1E_R_K(g12, g4, envs->i_l+0, envs->j_l+0, envs->k_l+0);
        G1E_R_K(g15, g7, envs->i_l+0, envs->j_l+0, envs->k_l+0);
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
void int3c1e_r4_origk_optimizer(CINTOpt **opt, int *atm, int natm, int *bas, int nbas, double *env) {
        int ng[] = {0, 0, 4, 0, 4, 1, 1, 1};
        CINTall_3c1e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
CACHE_SIZE_T int3c1e_r4_origk_cart(double *out, int *dims, int *shls,
                int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {0, 0, 4, 0, 4, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_int3c1e_r4_origk;
        return CINT3c1e_drv(out, dims, &envs, opt, cache, &c2s_cart_3c1e);
} // int3c1e_r4_origk_cart
CACHE_SIZE_T int3c1e_r4_origk_sph(double *out, int *dims, int *shls,
                int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {0, 0, 4, 0, 4, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_int3c1e_r4_origk;
        return CINT3c1e_drv(out, dims, &envs, opt, cache, &c2s_sph_3c1e);
} // int3c1e_r4_origk_sph
CACHE_SIZE_T int3c1e_r4_origk_spinor(double complex *out, int *dims, int *shls,
                int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {0, 0, 4, 0, 4, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_int3c1e_r4_origk;
return CINT3c1e_spinor_drv(out, dims, &envs, opt, cache, &c2s_sf_3c2e1);
} // int3c1e_r4_origk_spinor
ALL_CINT(int3c1e_r4_origk)
ALL_CINT_FORTRAN_(int3c1e_r4_origk)

/* based on
 * '("int3c1e_r6_origk" ( \, \, r dot r r dot r r dot r))
 */
static void CINTgout1e_int3c1e_r6_origk(double *RESTRICT gout,
                double *RESTRICT g, int *RESTRICT idx, CINTEnvVars *envs, int count) {
        CINTg3c1e_ovlp(g, envs, count);
        int nf = envs->nf;
        int nfc = nf * 1;
        int ix, iy, iz, n;
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
        double *RESTRICT g17 = g16 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g18 = g17 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g19 = g18 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g20 = g19 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g21 = g20 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g22 = g21 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g23 = g22 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g24 = g23 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g25 = g24 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g26 = g25 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g27 = g26 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g28 = g27 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g29 = g28 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g30 = g29 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g31 = g30 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g32 = g31 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g33 = g32 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g34 = g33 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g35 = g34 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g36 = g35 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g37 = g36 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g38 = g37 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g39 = g38 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g40 = g39 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g41 = g40 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g42 = g41 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g43 = g42 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g44 = g43 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g45 = g44 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g46 = g45 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g47 = g46 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g48 = g47 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g49 = g48 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g50 = g49 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g51 = g50 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g52 = g51 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g53 = g52 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g54 = g53 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g55 = g54 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g56 = g55 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g57 = g56 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g58 = g57 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g59 = g58 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g60 = g59 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g61 = g60 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g62 = g61 + envs->g_size * 3 * SIMDD;
        double *RESTRICT g63 = g62 + envs->g_size * 3 * SIMDD;
        __MD r1;
        __MD r3 = MM_SET1(3.);
        __MD r6 = MM_SET1(6.);
        G1E_R_K(g1, g0, envs->i_l+0, envs->j_l+0, envs->k_l+5);
        G1E_R_K(g3, g1, envs->i_l+0, envs->j_l+0, envs->k_l+4);
        G1E_R_K(g4, g0, envs->i_l+0, envs->j_l+0, envs->k_l+3);
        G1E_R_K(g7, g3, envs->i_l+0, envs->j_l+0, envs->k_l+3);
        G1E_R_K(g12, g4, envs->i_l+0, envs->j_l+0, envs->k_l+2);
        G1E_R_K(g15, g7, envs->i_l+0, envs->j_l+0, envs->k_l+2);
        G1E_R_K(g16, g0, envs->i_l+0, envs->j_l+0, envs->k_l+1);
        G1E_R_K(g28, g12, envs->i_l+0, envs->j_l+0, envs->k_l+1);
        G1E_R_K(g31, g15, envs->i_l+0, envs->j_l+0, envs->k_l+1);
        G1E_R_K(g48, g16, envs->i_l+0, envs->j_l+0, envs->k_l+0);
        G1E_R_K(g60, g28, envs->i_l+0, envs->j_l+0, envs->k_l+0);
        G1E_R_K(g63, g31, envs->i_l+0, envs->j_l+0, envs->k_l+0);
        for (n = 0; n < nf; n++) {
                ix = idx[0+n*3];
                iy = idx[1+n*3];
                iz = idx[2+n*3];
                r1 = MM_LOAD(g63+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
                r1+= MM_LOAD(g60+ix*SIMDD) * MM_LOAD(g3+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD) * r3;
                r1+= MM_LOAD(g60+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g3+iz*SIMDD) * r3;
                r1+= MM_LOAD(g48+ix*SIMDD) * MM_LOAD(g15+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD) * r3;
                r1+= MM_LOAD(g48+ix*SIMDD) * MM_LOAD(g12+iy*SIMDD) * MM_LOAD(g3+iz*SIMDD) * r6;
                r1+= MM_LOAD(g48+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g15+iz*SIMDD) * r3;
                r1+= MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g63+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
                r1+= MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g60+iy*SIMDD) * MM_LOAD(g3+iz*SIMDD) * r3;
                r1+= MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g48+iy*SIMDD) * MM_LOAD(g15+iz*SIMDD) * r3;
                r1+= MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g63+iz*SIMDD);
                GOUT_SCATTER(gout, n, r1);
        }
}
void int3c1e_r6_origk_optimizer(CINTOpt **opt, int *atm, int natm, int *bas, int nbas, double *env) {
        int ng[] = {0, 0, 6, 0, 6, 1, 1, 1};
        CINTall_3c1e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
CACHE_SIZE_T int3c1e_r6_origk_cart(double *out, int *dims, int *shls,
                int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {0, 0, 6, 0, 6, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_int3c1e_r6_origk;
        return CINT3c1e_drv(out, dims, &envs, opt, cache, &c2s_cart_3c1e);
} // int3c1e_r6_origk_cart
CACHE_SIZE_T int3c1e_r6_origk_sph(double *out, int *dims, int *shls,
                int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {0, 0, 6, 0, 6, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_int3c1e_r6_origk;
        return CINT3c1e_drv(out, dims, &envs, opt, cache, &c2s_sph_3c1e);
} // int3c1e_r6_origk_sph
CACHE_SIZE_T int3c1e_r6_origk_spinor(double complex *out, int *dims, int *shls,
                int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {0, 0, 6, 0, 6, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_int3c1e_r6_origk;
        return CINT3c1e_spinor_drv(out, dims, &envs, opt, cache, &c2s_sf_3c2e1);
} // int3c1e_r6_origk_spinor
ALL_CINT(int3c1e_r6_origk)
ALL_CINT_FORTRAN_(int3c1e_r6_origk)

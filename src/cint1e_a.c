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
void int1e_r2_origi_optimizer(CINTOpt **opt, FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env) {
        int ng[] = {2, 0, 0, 0, 2, 1, 1, 1};
        CINTall_1e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
CACHE_SIZE_T int1e_r2_origi_cart(double *out, FINT *dims, FINT *shls,
                FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {2, 0, 0, 0, 2, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_int1e_r2_origi;
        return CINT1e_drv(out, dims, &envs, opt, cache, &c2s_cart_1e);
} // int1e_r2_origi_cart
CACHE_SIZE_T int1e_r2_origi_sph(double *out, FINT *dims, FINT *shls,
                FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {2, 0, 0, 0, 2, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_int1e_r2_origi;
        return CINT1e_drv(out, dims, &envs, opt, cache, &c2s_sph_1e);
} // int1e_r2_origi_sph
CACHE_SIZE_T int1e_r2_origi_spinor(double complex *out, FINT *dims, FINT *shls,
                FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache) {
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
void int1e_r4_origi_optimizer(CINTOpt **opt, FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env) {
        int ng[] = {4, 0, 0, 0, 4, 1, 1, 1};
        CINTall_1e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
CACHE_SIZE_T int1e_r4_origi_cart(double *out, FINT *dims, FINT *shls,
                FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {4, 0, 0, 0, 4, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_int1e_r4_origi;
        return CINT1e_drv(out, dims, &envs, opt, cache, &c2s_cart_1e);
} // int1e_r4_origi_cart
CACHE_SIZE_T int1e_r4_origi_sph(double *out, FINT *dims, FINT *shls,
                FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {4, 0, 0, 0, 4, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_int1e_r4_origi;
        return CINT1e_drv(out, dims, &envs, opt, cache, &c2s_sph_1e);
} // int1e_r4_origi_sph
CACHE_SIZE_T int1e_r4_origi_spinor(double complex *out, FINT *dims, FINT *shls,
                FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {4, 0, 0, 0, 4, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_int1e_r4_origi;
        return CINT1e_spinor_drv(out, dims, &envs, opt, cache, &c2s_sf_1e);
} // int1e_r4_origi_spinor
ALL_CINT1E(int1e_r4_origi)
//ALL_CINT1E_FORTRAN_(cint1e_r4_origi)

/* <R DOT R i|OVLP |NABLA j> */
void CINTgout1e_int1e_r2_origi_ip2(double *gout, double *g, int *idx, CINTEnvVars *envs, int count) {
CINTg1e_ovlp(g, envs, count);
int nf = envs->nf;
int nfc = nf * 3;
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
__MD r1;
__MD rs[27];
G1E_D_J(g1, g0, envs->i_l+2, envs->j_l+0, 0);
G1E_R_I(g2, g0, envs->i_l+1, envs->j_l, 0);
G1E_R_I(g3, g1, envs->i_l+1, envs->j_l, 0);
G1E_R_I(g4, g0, envs->i_l+0, envs->j_l, 0);
G1E_R_I(g5, g1, envs->i_l+0, envs->j_l, 0);
G1E_R_I(g6, g2, envs->i_l+0, envs->j_l, 0);
G1E_R_I(g7, g3, envs->i_l+0, envs->j_l, 0);
for (n = 0; n < nf; n++) {
ix = idx[0+n*3];
iy = idx[1+n*3];
iz = idx[2+n*3];
rs[0] = MM_LOAD(g7+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[1] = MM_LOAD(g6+ix*SIMDD) * MM_LOAD(g1+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[2] = MM_LOAD(g6+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g1+iz*SIMDD);
rs[3] = MM_LOAD(g5+ix*SIMDD) * MM_LOAD(g2+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[4] = MM_LOAD(g4+ix*SIMDD) * MM_LOAD(g3+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[5] = MM_LOAD(g4+ix*SIMDD) * MM_LOAD(g2+iy*SIMDD) * MM_LOAD(g1+iz*SIMDD);
rs[6] = MM_LOAD(g5+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g2+iz*SIMDD);
rs[7] = MM_LOAD(g4+ix*SIMDD) * MM_LOAD(g1+iy*SIMDD) * MM_LOAD(g2+iz*SIMDD);
rs[8] = MM_LOAD(g4+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g3+iz*SIMDD);
rs[9] = MM_LOAD(g3+ix*SIMDD) * MM_LOAD(g4+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[10] = MM_LOAD(g2+ix*SIMDD) * MM_LOAD(g5+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[11] = MM_LOAD(g2+ix*SIMDD) * MM_LOAD(g4+iy*SIMDD) * MM_LOAD(g1+iz*SIMDD);
rs[12] = MM_LOAD(g1+ix*SIMDD) * MM_LOAD(g6+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[13] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g7+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[14] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g6+iy*SIMDD) * MM_LOAD(g1+iz*SIMDD);
rs[15] = MM_LOAD(g1+ix*SIMDD) * MM_LOAD(g4+iy*SIMDD) * MM_LOAD(g2+iz*SIMDD);
rs[16] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g5+iy*SIMDD) * MM_LOAD(g2+iz*SIMDD);
rs[17] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g4+iy*SIMDD) * MM_LOAD(g3+iz*SIMDD);
rs[18] = MM_LOAD(g3+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g4+iz*SIMDD);
rs[19] = MM_LOAD(g2+ix*SIMDD) * MM_LOAD(g1+iy*SIMDD) * MM_LOAD(g4+iz*SIMDD);
rs[20] = MM_LOAD(g2+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g5+iz*SIMDD);
rs[21] = MM_LOAD(g1+ix*SIMDD) * MM_LOAD(g2+iy*SIMDD) * MM_LOAD(g4+iz*SIMDD);
rs[22] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g3+iy*SIMDD) * MM_LOAD(g4+iz*SIMDD);
rs[23] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g2+iy*SIMDD) * MM_LOAD(g5+iz*SIMDD);
rs[24] = MM_LOAD(g1+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g6+iz*SIMDD);
rs[25] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g1+iy*SIMDD) * MM_LOAD(g6+iz*SIMDD);
rs[26] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g7+iz*SIMDD);
r1 = + rs[0] + rs[12] + rs[24]; GOUT_SCATTER(gout, n*3+0, r1);
r1 = + rs[1] + rs[13] + rs[25]; GOUT_SCATTER(gout, n*3+1, r1);
r1 = + rs[2] + rs[14] + rs[26]; GOUT_SCATTER(gout, n*3+2, r1);
}}
void int1e_r2_origi_ip2_optimizer(CINTOpt **opt, FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env) {
int ng[] = {2, 1, 0, 0, 3, 1, 1, 3};
CINTall_1e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
CACHE_SIZE_T int1e_r2_origi_ip2_cart(double *out, FINT *dims, FINT *shls,
FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache) {
int ng[] = {2, 1, 0, 0, 3, 1, 1, 3};
CINTEnvVars envs;
CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
envs.f_gout = &CINTgout1e_int1e_r2_origi_ip2;
return CINT1e_drv(out, dims, &envs, opt, cache, &c2s_cart_1e);
} // int1e_r2_origi_ip2_cart
CACHE_SIZE_T int1e_r2_origi_ip2_sph(double *out, FINT *dims, FINT *shls,
FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache) {
int ng[] = {2, 1, 0, 0, 3, 1, 1, 3};
CINTEnvVars envs;
CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
envs.f_gout = &CINTgout1e_int1e_r2_origi_ip2;
return CINT1e_drv(out, dims, &envs, opt, cache, &c2s_sph_1e);
} // int1e_r2_origi_ip2_sph
CACHE_SIZE_T int1e_r2_origi_ip2_spinor(double complex *out, FINT *dims, FINT *shls,
FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache) {
int ng[] = {2, 1, 0, 0, 3, 1, 1, 3};
CINTEnvVars envs;
CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
envs.f_gout = &CINTgout1e_int1e_r2_origi_ip2;
return CINT1e_spinor_drv(out, dims, &envs, opt, cache, &c2s_sf_1e);
} // int1e_r2_origi_ip2_spinor
ALL_CINT1E(int1e_r2_origi_ip2)
//ALL_CINT1E_FORTRAN_(cint1e_r2_origi_ip2)
/* <R DOT R R DOT R i|OVLP |NABLA j> */
void CINTgout1e_int1e_r4_origi_ip2(double *gout, double *g, int *idx, CINTEnvVars *envs, int count) {
CINTg1e_ovlp(g, envs, count);
int nf = envs->nf;
int nfc = nf * 3;
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
double *RESTRICT g16 = g15  + envs->g_size * 3 * SIMDD;
double *RESTRICT g17 = g16  + envs->g_size * 3 * SIMDD;
double *RESTRICT g18 = g17  + envs->g_size * 3 * SIMDD;
double *RESTRICT g19 = g18  + envs->g_size * 3 * SIMDD;
double *RESTRICT g20 = g19  + envs->g_size * 3 * SIMDD;
double *RESTRICT g21 = g20  + envs->g_size * 3 * SIMDD;
double *RESTRICT g22 = g21  + envs->g_size * 3 * SIMDD;
double *RESTRICT g23 = g22  + envs->g_size * 3 * SIMDD;
double *RESTRICT g24 = g23  + envs->g_size * 3 * SIMDD;
double *RESTRICT g25 = g24  + envs->g_size * 3 * SIMDD;
double *RESTRICT g26 = g25  + envs->g_size * 3 * SIMDD;
double *RESTRICT g27 = g26  + envs->g_size * 3 * SIMDD;
double *RESTRICT g28 = g27  + envs->g_size * 3 * SIMDD;
double *RESTRICT g29 = g28  + envs->g_size * 3 * SIMDD;
double *RESTRICT g30 = g29  + envs->g_size * 3 * SIMDD;
double *RESTRICT g31 = g30  + envs->g_size * 3 * SIMDD;
__MD r1;
__MD rs[243];
G1E_D_J(g1, g0, envs->i_l+4, envs->j_l+0, 0);
G1E_R_I(g2, g0, envs->i_l+3, envs->j_l, 0);
G1E_R_I(g3, g1, envs->i_l+3, envs->j_l, 0);
G1E_R_I(g4, g0, envs->i_l+2, envs->j_l, 0);
G1E_R_I(g5, g1, envs->i_l+2, envs->j_l, 0);
G1E_R_I(g6, g2, envs->i_l+2, envs->j_l, 0);
G1E_R_I(g7, g3, envs->i_l+2, envs->j_l, 0);
G1E_R_I(g8, g0, envs->i_l+1, envs->j_l, 0);
G1E_R_I(g9, g1, envs->i_l+1, envs->j_l, 0);
G1E_R_I(g10, g2, envs->i_l+1, envs->j_l, 0);
G1E_R_I(g11, g3, envs->i_l+1, envs->j_l, 0);
G1E_R_I(g12, g4, envs->i_l+1, envs->j_l, 0);
G1E_R_I(g13, g5, envs->i_l+1, envs->j_l, 0);
G1E_R_I(g14, g6, envs->i_l+1, envs->j_l, 0);
G1E_R_I(g15, g7, envs->i_l+1, envs->j_l, 0);
G1E_R_I(g16, g0, envs->i_l+0, envs->j_l, 0);
G1E_R_I(g17, g1, envs->i_l+0, envs->j_l, 0);
G1E_R_I(g18, g2, envs->i_l+0, envs->j_l, 0);
G1E_R_I(g19, g3, envs->i_l+0, envs->j_l, 0);
G1E_R_I(g20, g4, envs->i_l+0, envs->j_l, 0);
G1E_R_I(g21, g5, envs->i_l+0, envs->j_l, 0);
G1E_R_I(g22, g6, envs->i_l+0, envs->j_l, 0);
G1E_R_I(g23, g7, envs->i_l+0, envs->j_l, 0);
G1E_R_I(g24, g8, envs->i_l+0, envs->j_l, 0);
G1E_R_I(g25, g9, envs->i_l+0, envs->j_l, 0);
G1E_R_I(g26, g10, envs->i_l+0, envs->j_l, 0);
G1E_R_I(g27, g11, envs->i_l+0, envs->j_l, 0);
G1E_R_I(g28, g12, envs->i_l+0, envs->j_l, 0);
G1E_R_I(g29, g13, envs->i_l+0, envs->j_l, 0);
G1E_R_I(g30, g14, envs->i_l+0, envs->j_l, 0);
G1E_R_I(g31, g15, envs->i_l+0, envs->j_l, 0);
for (n = 0; n < nf; n++) {
ix = idx[0+n*3];
iy = idx[1+n*3];
iz = idx[2+n*3];
rs[0] = MM_LOAD(g31+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[1] = MM_LOAD(g30+ix*SIMDD) * MM_LOAD(g1+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[2] = MM_LOAD(g30+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g1+iz*SIMDD);
rs[3] = MM_LOAD(g29+ix*SIMDD) * MM_LOAD(g2+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[4] = MM_LOAD(g28+ix*SIMDD) * MM_LOAD(g3+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[5] = MM_LOAD(g28+ix*SIMDD) * MM_LOAD(g2+iy*SIMDD) * MM_LOAD(g1+iz*SIMDD);
rs[6] = MM_LOAD(g29+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g2+iz*SIMDD);
rs[7] = MM_LOAD(g28+ix*SIMDD) * MM_LOAD(g1+iy*SIMDD) * MM_LOAD(g2+iz*SIMDD);
rs[8] = MM_LOAD(g28+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g3+iz*SIMDD);
rs[9] = MM_LOAD(g27+ix*SIMDD) * MM_LOAD(g4+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[10] = MM_LOAD(g26+ix*SIMDD) * MM_LOAD(g5+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[11] = MM_LOAD(g26+ix*SIMDD) * MM_LOAD(g4+iy*SIMDD) * MM_LOAD(g1+iz*SIMDD);
rs[12] = MM_LOAD(g25+ix*SIMDD) * MM_LOAD(g6+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[13] = MM_LOAD(g24+ix*SIMDD) * MM_LOAD(g7+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[14] = MM_LOAD(g24+ix*SIMDD) * MM_LOAD(g6+iy*SIMDD) * MM_LOAD(g1+iz*SIMDD);
rs[15] = MM_LOAD(g25+ix*SIMDD) * MM_LOAD(g4+iy*SIMDD) * MM_LOAD(g2+iz*SIMDD);
rs[16] = MM_LOAD(g24+ix*SIMDD) * MM_LOAD(g5+iy*SIMDD) * MM_LOAD(g2+iz*SIMDD);
rs[17] = MM_LOAD(g24+ix*SIMDD) * MM_LOAD(g4+iy*SIMDD) * MM_LOAD(g3+iz*SIMDD);
rs[18] = MM_LOAD(g27+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g4+iz*SIMDD);
rs[19] = MM_LOAD(g26+ix*SIMDD) * MM_LOAD(g1+iy*SIMDD) * MM_LOAD(g4+iz*SIMDD);
rs[20] = MM_LOAD(g26+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g5+iz*SIMDD);
rs[21] = MM_LOAD(g25+ix*SIMDD) * MM_LOAD(g2+iy*SIMDD) * MM_LOAD(g4+iz*SIMDD);
rs[22] = MM_LOAD(g24+ix*SIMDD) * MM_LOAD(g3+iy*SIMDD) * MM_LOAD(g4+iz*SIMDD);
rs[23] = MM_LOAD(g24+ix*SIMDD) * MM_LOAD(g2+iy*SIMDD) * MM_LOAD(g5+iz*SIMDD);
rs[24] = MM_LOAD(g25+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g6+iz*SIMDD);
rs[25] = MM_LOAD(g24+ix*SIMDD) * MM_LOAD(g1+iy*SIMDD) * MM_LOAD(g6+iz*SIMDD);
rs[26] = MM_LOAD(g24+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g7+iz*SIMDD);
rs[27] = MM_LOAD(g23+ix*SIMDD) * MM_LOAD(g8+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[28] = MM_LOAD(g22+ix*SIMDD) * MM_LOAD(g9+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[29] = MM_LOAD(g22+ix*SIMDD) * MM_LOAD(g8+iy*SIMDD) * MM_LOAD(g1+iz*SIMDD);
rs[30] = MM_LOAD(g21+ix*SIMDD) * MM_LOAD(g10+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[31] = MM_LOAD(g20+ix*SIMDD) * MM_LOAD(g11+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[32] = MM_LOAD(g20+ix*SIMDD) * MM_LOAD(g10+iy*SIMDD) * MM_LOAD(g1+iz*SIMDD);
rs[33] = MM_LOAD(g21+ix*SIMDD) * MM_LOAD(g8+iy*SIMDD) * MM_LOAD(g2+iz*SIMDD);
rs[34] = MM_LOAD(g20+ix*SIMDD) * MM_LOAD(g9+iy*SIMDD) * MM_LOAD(g2+iz*SIMDD);
rs[35] = MM_LOAD(g20+ix*SIMDD) * MM_LOAD(g8+iy*SIMDD) * MM_LOAD(g3+iz*SIMDD);
rs[36] = MM_LOAD(g19+ix*SIMDD) * MM_LOAD(g12+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[37] = MM_LOAD(g18+ix*SIMDD) * MM_LOAD(g13+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[38] = MM_LOAD(g18+ix*SIMDD) * MM_LOAD(g12+iy*SIMDD) * MM_LOAD(g1+iz*SIMDD);
rs[39] = MM_LOAD(g17+ix*SIMDD) * MM_LOAD(g14+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[40] = MM_LOAD(g16+ix*SIMDD) * MM_LOAD(g15+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[41] = MM_LOAD(g16+ix*SIMDD) * MM_LOAD(g14+iy*SIMDD) * MM_LOAD(g1+iz*SIMDD);
rs[42] = MM_LOAD(g17+ix*SIMDD) * MM_LOAD(g12+iy*SIMDD) * MM_LOAD(g2+iz*SIMDD);
rs[43] = MM_LOAD(g16+ix*SIMDD) * MM_LOAD(g13+iy*SIMDD) * MM_LOAD(g2+iz*SIMDD);
rs[44] = MM_LOAD(g16+ix*SIMDD) * MM_LOAD(g12+iy*SIMDD) * MM_LOAD(g3+iz*SIMDD);
rs[45] = MM_LOAD(g19+ix*SIMDD) * MM_LOAD(g8+iy*SIMDD) * MM_LOAD(g4+iz*SIMDD);
rs[46] = MM_LOAD(g18+ix*SIMDD) * MM_LOAD(g9+iy*SIMDD) * MM_LOAD(g4+iz*SIMDD);
rs[47] = MM_LOAD(g18+ix*SIMDD) * MM_LOAD(g8+iy*SIMDD) * MM_LOAD(g5+iz*SIMDD);
rs[48] = MM_LOAD(g17+ix*SIMDD) * MM_LOAD(g10+iy*SIMDD) * MM_LOAD(g4+iz*SIMDD);
rs[49] = MM_LOAD(g16+ix*SIMDD) * MM_LOAD(g11+iy*SIMDD) * MM_LOAD(g4+iz*SIMDD);
rs[50] = MM_LOAD(g16+ix*SIMDD) * MM_LOAD(g10+iy*SIMDD) * MM_LOAD(g5+iz*SIMDD);
rs[51] = MM_LOAD(g17+ix*SIMDD) * MM_LOAD(g8+iy*SIMDD) * MM_LOAD(g6+iz*SIMDD);
rs[52] = MM_LOAD(g16+ix*SIMDD) * MM_LOAD(g9+iy*SIMDD) * MM_LOAD(g6+iz*SIMDD);
rs[53] = MM_LOAD(g16+ix*SIMDD) * MM_LOAD(g8+iy*SIMDD) * MM_LOAD(g7+iz*SIMDD);
rs[54] = MM_LOAD(g23+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g8+iz*SIMDD);
rs[55] = MM_LOAD(g22+ix*SIMDD) * MM_LOAD(g1+iy*SIMDD) * MM_LOAD(g8+iz*SIMDD);
rs[56] = MM_LOAD(g22+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g9+iz*SIMDD);
rs[57] = MM_LOAD(g21+ix*SIMDD) * MM_LOAD(g2+iy*SIMDD) * MM_LOAD(g8+iz*SIMDD);
rs[58] = MM_LOAD(g20+ix*SIMDD) * MM_LOAD(g3+iy*SIMDD) * MM_LOAD(g8+iz*SIMDD);
rs[59] = MM_LOAD(g20+ix*SIMDD) * MM_LOAD(g2+iy*SIMDD) * MM_LOAD(g9+iz*SIMDD);
rs[60] = MM_LOAD(g21+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g10+iz*SIMDD);
rs[61] = MM_LOAD(g20+ix*SIMDD) * MM_LOAD(g1+iy*SIMDD) * MM_LOAD(g10+iz*SIMDD);
rs[62] = MM_LOAD(g20+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g11+iz*SIMDD);
rs[63] = MM_LOAD(g19+ix*SIMDD) * MM_LOAD(g4+iy*SIMDD) * MM_LOAD(g8+iz*SIMDD);
rs[64] = MM_LOAD(g18+ix*SIMDD) * MM_LOAD(g5+iy*SIMDD) * MM_LOAD(g8+iz*SIMDD);
rs[65] = MM_LOAD(g18+ix*SIMDD) * MM_LOAD(g4+iy*SIMDD) * MM_LOAD(g9+iz*SIMDD);
rs[66] = MM_LOAD(g17+ix*SIMDD) * MM_LOAD(g6+iy*SIMDD) * MM_LOAD(g8+iz*SIMDD);
rs[67] = MM_LOAD(g16+ix*SIMDD) * MM_LOAD(g7+iy*SIMDD) * MM_LOAD(g8+iz*SIMDD);
rs[68] = MM_LOAD(g16+ix*SIMDD) * MM_LOAD(g6+iy*SIMDD) * MM_LOAD(g9+iz*SIMDD);
rs[69] = MM_LOAD(g17+ix*SIMDD) * MM_LOAD(g4+iy*SIMDD) * MM_LOAD(g10+iz*SIMDD);
rs[70] = MM_LOAD(g16+ix*SIMDD) * MM_LOAD(g5+iy*SIMDD) * MM_LOAD(g10+iz*SIMDD);
rs[71] = MM_LOAD(g16+ix*SIMDD) * MM_LOAD(g4+iy*SIMDD) * MM_LOAD(g11+iz*SIMDD);
rs[72] = MM_LOAD(g19+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g12+iz*SIMDD);
rs[73] = MM_LOAD(g18+ix*SIMDD) * MM_LOAD(g1+iy*SIMDD) * MM_LOAD(g12+iz*SIMDD);
rs[74] = MM_LOAD(g18+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g13+iz*SIMDD);
rs[75] = MM_LOAD(g17+ix*SIMDD) * MM_LOAD(g2+iy*SIMDD) * MM_LOAD(g12+iz*SIMDD);
rs[76] = MM_LOAD(g16+ix*SIMDD) * MM_LOAD(g3+iy*SIMDD) * MM_LOAD(g12+iz*SIMDD);
rs[77] = MM_LOAD(g16+ix*SIMDD) * MM_LOAD(g2+iy*SIMDD) * MM_LOAD(g13+iz*SIMDD);
rs[78] = MM_LOAD(g17+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g14+iz*SIMDD);
rs[79] = MM_LOAD(g16+ix*SIMDD) * MM_LOAD(g1+iy*SIMDD) * MM_LOAD(g14+iz*SIMDD);
rs[80] = MM_LOAD(g16+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g15+iz*SIMDD);
rs[81] = MM_LOAD(g15+ix*SIMDD) * MM_LOAD(g16+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[82] = MM_LOAD(g14+ix*SIMDD) * MM_LOAD(g17+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[83] = MM_LOAD(g14+ix*SIMDD) * MM_LOAD(g16+iy*SIMDD) * MM_LOAD(g1+iz*SIMDD);
rs[84] = MM_LOAD(g13+ix*SIMDD) * MM_LOAD(g18+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[85] = MM_LOAD(g12+ix*SIMDD) * MM_LOAD(g19+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[86] = MM_LOAD(g12+ix*SIMDD) * MM_LOAD(g18+iy*SIMDD) * MM_LOAD(g1+iz*SIMDD);
rs[87] = MM_LOAD(g13+ix*SIMDD) * MM_LOAD(g16+iy*SIMDD) * MM_LOAD(g2+iz*SIMDD);
rs[88] = MM_LOAD(g12+ix*SIMDD) * MM_LOAD(g17+iy*SIMDD) * MM_LOAD(g2+iz*SIMDD);
rs[89] = MM_LOAD(g12+ix*SIMDD) * MM_LOAD(g16+iy*SIMDD) * MM_LOAD(g3+iz*SIMDD);
rs[90] = MM_LOAD(g11+ix*SIMDD) * MM_LOAD(g20+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[91] = MM_LOAD(g10+ix*SIMDD) * MM_LOAD(g21+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[92] = MM_LOAD(g10+ix*SIMDD) * MM_LOAD(g20+iy*SIMDD) * MM_LOAD(g1+iz*SIMDD);
rs[93] = MM_LOAD(g9+ix*SIMDD) * MM_LOAD(g22+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[94] = MM_LOAD(g8+ix*SIMDD) * MM_LOAD(g23+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[95] = MM_LOAD(g8+ix*SIMDD) * MM_LOAD(g22+iy*SIMDD) * MM_LOAD(g1+iz*SIMDD);
rs[96] = MM_LOAD(g9+ix*SIMDD) * MM_LOAD(g20+iy*SIMDD) * MM_LOAD(g2+iz*SIMDD);
rs[97] = MM_LOAD(g8+ix*SIMDD) * MM_LOAD(g21+iy*SIMDD) * MM_LOAD(g2+iz*SIMDD);
rs[98] = MM_LOAD(g8+ix*SIMDD) * MM_LOAD(g20+iy*SIMDD) * MM_LOAD(g3+iz*SIMDD);
rs[99] = MM_LOAD(g11+ix*SIMDD) * MM_LOAD(g16+iy*SIMDD) * MM_LOAD(g4+iz*SIMDD);
rs[100] = MM_LOAD(g10+ix*SIMDD) * MM_LOAD(g17+iy*SIMDD) * MM_LOAD(g4+iz*SIMDD);
rs[101] = MM_LOAD(g10+ix*SIMDD) * MM_LOAD(g16+iy*SIMDD) * MM_LOAD(g5+iz*SIMDD);
rs[102] = MM_LOAD(g9+ix*SIMDD) * MM_LOAD(g18+iy*SIMDD) * MM_LOAD(g4+iz*SIMDD);
rs[103] = MM_LOAD(g8+ix*SIMDD) * MM_LOAD(g19+iy*SIMDD) * MM_LOAD(g4+iz*SIMDD);
rs[104] = MM_LOAD(g8+ix*SIMDD) * MM_LOAD(g18+iy*SIMDD) * MM_LOAD(g5+iz*SIMDD);
rs[105] = MM_LOAD(g9+ix*SIMDD) * MM_LOAD(g16+iy*SIMDD) * MM_LOAD(g6+iz*SIMDD);
rs[106] = MM_LOAD(g8+ix*SIMDD) * MM_LOAD(g17+iy*SIMDD) * MM_LOAD(g6+iz*SIMDD);
rs[107] = MM_LOAD(g8+ix*SIMDD) * MM_LOAD(g16+iy*SIMDD) * MM_LOAD(g7+iz*SIMDD);
rs[108] = MM_LOAD(g7+ix*SIMDD) * MM_LOAD(g24+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[109] = MM_LOAD(g6+ix*SIMDD) * MM_LOAD(g25+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[110] = MM_LOAD(g6+ix*SIMDD) * MM_LOAD(g24+iy*SIMDD) * MM_LOAD(g1+iz*SIMDD);
rs[111] = MM_LOAD(g5+ix*SIMDD) * MM_LOAD(g26+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[112] = MM_LOAD(g4+ix*SIMDD) * MM_LOAD(g27+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[113] = MM_LOAD(g4+ix*SIMDD) * MM_LOAD(g26+iy*SIMDD) * MM_LOAD(g1+iz*SIMDD);
rs[114] = MM_LOAD(g5+ix*SIMDD) * MM_LOAD(g24+iy*SIMDD) * MM_LOAD(g2+iz*SIMDD);
rs[115] = MM_LOAD(g4+ix*SIMDD) * MM_LOAD(g25+iy*SIMDD) * MM_LOAD(g2+iz*SIMDD);
rs[116] = MM_LOAD(g4+ix*SIMDD) * MM_LOAD(g24+iy*SIMDD) * MM_LOAD(g3+iz*SIMDD);
rs[117] = MM_LOAD(g3+ix*SIMDD) * MM_LOAD(g28+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[118] = MM_LOAD(g2+ix*SIMDD) * MM_LOAD(g29+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[119] = MM_LOAD(g2+ix*SIMDD) * MM_LOAD(g28+iy*SIMDD) * MM_LOAD(g1+iz*SIMDD);
rs[120] = MM_LOAD(g1+ix*SIMDD) * MM_LOAD(g30+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[121] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g31+iy*SIMDD) * MM_LOAD(g0+iz*SIMDD);
rs[122] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g30+iy*SIMDD) * MM_LOAD(g1+iz*SIMDD);
rs[123] = MM_LOAD(g1+ix*SIMDD) * MM_LOAD(g28+iy*SIMDD) * MM_LOAD(g2+iz*SIMDD);
rs[124] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g29+iy*SIMDD) * MM_LOAD(g2+iz*SIMDD);
rs[125] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g28+iy*SIMDD) * MM_LOAD(g3+iz*SIMDD);
rs[126] = MM_LOAD(g3+ix*SIMDD) * MM_LOAD(g24+iy*SIMDD) * MM_LOAD(g4+iz*SIMDD);
rs[127] = MM_LOAD(g2+ix*SIMDD) * MM_LOAD(g25+iy*SIMDD) * MM_LOAD(g4+iz*SIMDD);
rs[128] = MM_LOAD(g2+ix*SIMDD) * MM_LOAD(g24+iy*SIMDD) * MM_LOAD(g5+iz*SIMDD);
rs[129] = MM_LOAD(g1+ix*SIMDD) * MM_LOAD(g26+iy*SIMDD) * MM_LOAD(g4+iz*SIMDD);
rs[130] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g27+iy*SIMDD) * MM_LOAD(g4+iz*SIMDD);
rs[131] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g26+iy*SIMDD) * MM_LOAD(g5+iz*SIMDD);
rs[132] = MM_LOAD(g1+ix*SIMDD) * MM_LOAD(g24+iy*SIMDD) * MM_LOAD(g6+iz*SIMDD);
rs[133] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g25+iy*SIMDD) * MM_LOAD(g6+iz*SIMDD);
rs[134] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g24+iy*SIMDD) * MM_LOAD(g7+iz*SIMDD);
rs[135] = MM_LOAD(g7+ix*SIMDD) * MM_LOAD(g16+iy*SIMDD) * MM_LOAD(g8+iz*SIMDD);
rs[136] = MM_LOAD(g6+ix*SIMDD) * MM_LOAD(g17+iy*SIMDD) * MM_LOAD(g8+iz*SIMDD);
rs[137] = MM_LOAD(g6+ix*SIMDD) * MM_LOAD(g16+iy*SIMDD) * MM_LOAD(g9+iz*SIMDD);
rs[138] = MM_LOAD(g5+ix*SIMDD) * MM_LOAD(g18+iy*SIMDD) * MM_LOAD(g8+iz*SIMDD);
rs[139] = MM_LOAD(g4+ix*SIMDD) * MM_LOAD(g19+iy*SIMDD) * MM_LOAD(g8+iz*SIMDD);
rs[140] = MM_LOAD(g4+ix*SIMDD) * MM_LOAD(g18+iy*SIMDD) * MM_LOAD(g9+iz*SIMDD);
rs[141] = MM_LOAD(g5+ix*SIMDD) * MM_LOAD(g16+iy*SIMDD) * MM_LOAD(g10+iz*SIMDD);
rs[142] = MM_LOAD(g4+ix*SIMDD) * MM_LOAD(g17+iy*SIMDD) * MM_LOAD(g10+iz*SIMDD);
rs[143] = MM_LOAD(g4+ix*SIMDD) * MM_LOAD(g16+iy*SIMDD) * MM_LOAD(g11+iz*SIMDD);
rs[144] = MM_LOAD(g3+ix*SIMDD) * MM_LOAD(g20+iy*SIMDD) * MM_LOAD(g8+iz*SIMDD);
rs[145] = MM_LOAD(g2+ix*SIMDD) * MM_LOAD(g21+iy*SIMDD) * MM_LOAD(g8+iz*SIMDD);
rs[146] = MM_LOAD(g2+ix*SIMDD) * MM_LOAD(g20+iy*SIMDD) * MM_LOAD(g9+iz*SIMDD);
rs[147] = MM_LOAD(g1+ix*SIMDD) * MM_LOAD(g22+iy*SIMDD) * MM_LOAD(g8+iz*SIMDD);
rs[148] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g23+iy*SIMDD) * MM_LOAD(g8+iz*SIMDD);
rs[149] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g22+iy*SIMDD) * MM_LOAD(g9+iz*SIMDD);
rs[150] = MM_LOAD(g1+ix*SIMDD) * MM_LOAD(g20+iy*SIMDD) * MM_LOAD(g10+iz*SIMDD);
rs[151] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g21+iy*SIMDD) * MM_LOAD(g10+iz*SIMDD);
rs[152] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g20+iy*SIMDD) * MM_LOAD(g11+iz*SIMDD);
rs[153] = MM_LOAD(g3+ix*SIMDD) * MM_LOAD(g16+iy*SIMDD) * MM_LOAD(g12+iz*SIMDD);
rs[154] = MM_LOAD(g2+ix*SIMDD) * MM_LOAD(g17+iy*SIMDD) * MM_LOAD(g12+iz*SIMDD);
rs[155] = MM_LOAD(g2+ix*SIMDD) * MM_LOAD(g16+iy*SIMDD) * MM_LOAD(g13+iz*SIMDD);
rs[156] = MM_LOAD(g1+ix*SIMDD) * MM_LOAD(g18+iy*SIMDD) * MM_LOAD(g12+iz*SIMDD);
rs[157] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g19+iy*SIMDD) * MM_LOAD(g12+iz*SIMDD);
rs[158] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g18+iy*SIMDD) * MM_LOAD(g13+iz*SIMDD);
rs[159] = MM_LOAD(g1+ix*SIMDD) * MM_LOAD(g16+iy*SIMDD) * MM_LOAD(g14+iz*SIMDD);
rs[160] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g17+iy*SIMDD) * MM_LOAD(g14+iz*SIMDD);
rs[161] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g16+iy*SIMDD) * MM_LOAD(g15+iz*SIMDD);
rs[162] = MM_LOAD(g15+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g16+iz*SIMDD);
rs[163] = MM_LOAD(g14+ix*SIMDD) * MM_LOAD(g1+iy*SIMDD) * MM_LOAD(g16+iz*SIMDD);
rs[164] = MM_LOAD(g14+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g17+iz*SIMDD);
rs[165] = MM_LOAD(g13+ix*SIMDD) * MM_LOAD(g2+iy*SIMDD) * MM_LOAD(g16+iz*SIMDD);
rs[166] = MM_LOAD(g12+ix*SIMDD) * MM_LOAD(g3+iy*SIMDD) * MM_LOAD(g16+iz*SIMDD);
rs[167] = MM_LOAD(g12+ix*SIMDD) * MM_LOAD(g2+iy*SIMDD) * MM_LOAD(g17+iz*SIMDD);
rs[168] = MM_LOAD(g13+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g18+iz*SIMDD);
rs[169] = MM_LOAD(g12+ix*SIMDD) * MM_LOAD(g1+iy*SIMDD) * MM_LOAD(g18+iz*SIMDD);
rs[170] = MM_LOAD(g12+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g19+iz*SIMDD);
rs[171] = MM_LOAD(g11+ix*SIMDD) * MM_LOAD(g4+iy*SIMDD) * MM_LOAD(g16+iz*SIMDD);
rs[172] = MM_LOAD(g10+ix*SIMDD) * MM_LOAD(g5+iy*SIMDD) * MM_LOAD(g16+iz*SIMDD);
rs[173] = MM_LOAD(g10+ix*SIMDD) * MM_LOAD(g4+iy*SIMDD) * MM_LOAD(g17+iz*SIMDD);
rs[174] = MM_LOAD(g9+ix*SIMDD) * MM_LOAD(g6+iy*SIMDD) * MM_LOAD(g16+iz*SIMDD);
rs[175] = MM_LOAD(g8+ix*SIMDD) * MM_LOAD(g7+iy*SIMDD) * MM_LOAD(g16+iz*SIMDD);
rs[176] = MM_LOAD(g8+ix*SIMDD) * MM_LOAD(g6+iy*SIMDD) * MM_LOAD(g17+iz*SIMDD);
rs[177] = MM_LOAD(g9+ix*SIMDD) * MM_LOAD(g4+iy*SIMDD) * MM_LOAD(g18+iz*SIMDD);
rs[178] = MM_LOAD(g8+ix*SIMDD) * MM_LOAD(g5+iy*SIMDD) * MM_LOAD(g18+iz*SIMDD);
rs[179] = MM_LOAD(g8+ix*SIMDD) * MM_LOAD(g4+iy*SIMDD) * MM_LOAD(g19+iz*SIMDD);
rs[180] = MM_LOAD(g11+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g20+iz*SIMDD);
rs[181] = MM_LOAD(g10+ix*SIMDD) * MM_LOAD(g1+iy*SIMDD) * MM_LOAD(g20+iz*SIMDD);
rs[182] = MM_LOAD(g10+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g21+iz*SIMDD);
rs[183] = MM_LOAD(g9+ix*SIMDD) * MM_LOAD(g2+iy*SIMDD) * MM_LOAD(g20+iz*SIMDD);
rs[184] = MM_LOAD(g8+ix*SIMDD) * MM_LOAD(g3+iy*SIMDD) * MM_LOAD(g20+iz*SIMDD);
rs[185] = MM_LOAD(g8+ix*SIMDD) * MM_LOAD(g2+iy*SIMDD) * MM_LOAD(g21+iz*SIMDD);
rs[186] = MM_LOAD(g9+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g22+iz*SIMDD);
rs[187] = MM_LOAD(g8+ix*SIMDD) * MM_LOAD(g1+iy*SIMDD) * MM_LOAD(g22+iz*SIMDD);
rs[188] = MM_LOAD(g8+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g23+iz*SIMDD);
rs[189] = MM_LOAD(g7+ix*SIMDD) * MM_LOAD(g8+iy*SIMDD) * MM_LOAD(g16+iz*SIMDD);
rs[190] = MM_LOAD(g6+ix*SIMDD) * MM_LOAD(g9+iy*SIMDD) * MM_LOAD(g16+iz*SIMDD);
rs[191] = MM_LOAD(g6+ix*SIMDD) * MM_LOAD(g8+iy*SIMDD) * MM_LOAD(g17+iz*SIMDD);
rs[192] = MM_LOAD(g5+ix*SIMDD) * MM_LOAD(g10+iy*SIMDD) * MM_LOAD(g16+iz*SIMDD);
rs[193] = MM_LOAD(g4+ix*SIMDD) * MM_LOAD(g11+iy*SIMDD) * MM_LOAD(g16+iz*SIMDD);
rs[194] = MM_LOAD(g4+ix*SIMDD) * MM_LOAD(g10+iy*SIMDD) * MM_LOAD(g17+iz*SIMDD);
rs[195] = MM_LOAD(g5+ix*SIMDD) * MM_LOAD(g8+iy*SIMDD) * MM_LOAD(g18+iz*SIMDD);
rs[196] = MM_LOAD(g4+ix*SIMDD) * MM_LOAD(g9+iy*SIMDD) * MM_LOAD(g18+iz*SIMDD);
rs[197] = MM_LOAD(g4+ix*SIMDD) * MM_LOAD(g8+iy*SIMDD) * MM_LOAD(g19+iz*SIMDD);
rs[198] = MM_LOAD(g3+ix*SIMDD) * MM_LOAD(g12+iy*SIMDD) * MM_LOAD(g16+iz*SIMDD);
rs[199] = MM_LOAD(g2+ix*SIMDD) * MM_LOAD(g13+iy*SIMDD) * MM_LOAD(g16+iz*SIMDD);
rs[200] = MM_LOAD(g2+ix*SIMDD) * MM_LOAD(g12+iy*SIMDD) * MM_LOAD(g17+iz*SIMDD);
rs[201] = MM_LOAD(g1+ix*SIMDD) * MM_LOAD(g14+iy*SIMDD) * MM_LOAD(g16+iz*SIMDD);
rs[202] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g15+iy*SIMDD) * MM_LOAD(g16+iz*SIMDD);
rs[203] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g14+iy*SIMDD) * MM_LOAD(g17+iz*SIMDD);
rs[204] = MM_LOAD(g1+ix*SIMDD) * MM_LOAD(g12+iy*SIMDD) * MM_LOAD(g18+iz*SIMDD);
rs[205] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g13+iy*SIMDD) * MM_LOAD(g18+iz*SIMDD);
rs[206] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g12+iy*SIMDD) * MM_LOAD(g19+iz*SIMDD);
rs[207] = MM_LOAD(g3+ix*SIMDD) * MM_LOAD(g8+iy*SIMDD) * MM_LOAD(g20+iz*SIMDD);
rs[208] = MM_LOAD(g2+ix*SIMDD) * MM_LOAD(g9+iy*SIMDD) * MM_LOAD(g20+iz*SIMDD);
rs[209] = MM_LOAD(g2+ix*SIMDD) * MM_LOAD(g8+iy*SIMDD) * MM_LOAD(g21+iz*SIMDD);
rs[210] = MM_LOAD(g1+ix*SIMDD) * MM_LOAD(g10+iy*SIMDD) * MM_LOAD(g20+iz*SIMDD);
rs[211] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g11+iy*SIMDD) * MM_LOAD(g20+iz*SIMDD);
rs[212] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g10+iy*SIMDD) * MM_LOAD(g21+iz*SIMDD);
rs[213] = MM_LOAD(g1+ix*SIMDD) * MM_LOAD(g8+iy*SIMDD) * MM_LOAD(g22+iz*SIMDD);
rs[214] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g9+iy*SIMDD) * MM_LOAD(g22+iz*SIMDD);
rs[215] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g8+iy*SIMDD) * MM_LOAD(g23+iz*SIMDD);
rs[216] = MM_LOAD(g7+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g24+iz*SIMDD);
rs[217] = MM_LOAD(g6+ix*SIMDD) * MM_LOAD(g1+iy*SIMDD) * MM_LOAD(g24+iz*SIMDD);
rs[218] = MM_LOAD(g6+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g25+iz*SIMDD);
rs[219] = MM_LOAD(g5+ix*SIMDD) * MM_LOAD(g2+iy*SIMDD) * MM_LOAD(g24+iz*SIMDD);
rs[220] = MM_LOAD(g4+ix*SIMDD) * MM_LOAD(g3+iy*SIMDD) * MM_LOAD(g24+iz*SIMDD);
rs[221] = MM_LOAD(g4+ix*SIMDD) * MM_LOAD(g2+iy*SIMDD) * MM_LOAD(g25+iz*SIMDD);
rs[222] = MM_LOAD(g5+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g26+iz*SIMDD);
rs[223] = MM_LOAD(g4+ix*SIMDD) * MM_LOAD(g1+iy*SIMDD) * MM_LOAD(g26+iz*SIMDD);
rs[224] = MM_LOAD(g4+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g27+iz*SIMDD);
rs[225] = MM_LOAD(g3+ix*SIMDD) * MM_LOAD(g4+iy*SIMDD) * MM_LOAD(g24+iz*SIMDD);
rs[226] = MM_LOAD(g2+ix*SIMDD) * MM_LOAD(g5+iy*SIMDD) * MM_LOAD(g24+iz*SIMDD);
rs[227] = MM_LOAD(g2+ix*SIMDD) * MM_LOAD(g4+iy*SIMDD) * MM_LOAD(g25+iz*SIMDD);
rs[228] = MM_LOAD(g1+ix*SIMDD) * MM_LOAD(g6+iy*SIMDD) * MM_LOAD(g24+iz*SIMDD);
rs[229] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g7+iy*SIMDD) * MM_LOAD(g24+iz*SIMDD);
rs[230] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g6+iy*SIMDD) * MM_LOAD(g25+iz*SIMDD);
rs[231] = MM_LOAD(g1+ix*SIMDD) * MM_LOAD(g4+iy*SIMDD) * MM_LOAD(g26+iz*SIMDD);
rs[232] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g5+iy*SIMDD) * MM_LOAD(g26+iz*SIMDD);
rs[233] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g4+iy*SIMDD) * MM_LOAD(g27+iz*SIMDD);
rs[234] = MM_LOAD(g3+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g28+iz*SIMDD);
rs[235] = MM_LOAD(g2+ix*SIMDD) * MM_LOAD(g1+iy*SIMDD) * MM_LOAD(g28+iz*SIMDD);
rs[236] = MM_LOAD(g2+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g29+iz*SIMDD);
rs[237] = MM_LOAD(g1+ix*SIMDD) * MM_LOAD(g2+iy*SIMDD) * MM_LOAD(g28+iz*SIMDD);
rs[238] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g3+iy*SIMDD) * MM_LOAD(g28+iz*SIMDD);
rs[239] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g2+iy*SIMDD) * MM_LOAD(g29+iz*SIMDD);
rs[240] = MM_LOAD(g1+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g30+iz*SIMDD);
rs[241] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g1+iy*SIMDD) * MM_LOAD(g30+iz*SIMDD);
rs[242] = MM_LOAD(g0+ix*SIMDD) * MM_LOAD(g0+iy*SIMDD) * MM_LOAD(g31+iz*SIMDD);
r1 = + rs[0] + MM_SET1(2)*rs[12] + MM_SET1(2)*rs[24] + rs[120] + MM_SET1(2)*rs[132] + rs[240]; GOUT_SCATTER(gout, n*3+0, r1);
r1 = + rs[1] + MM_SET1(2)*rs[13] + MM_SET1(2)*rs[25] + rs[121] + MM_SET1(2)*rs[133] + rs[241]; GOUT_SCATTER(gout, n*3+1, r1);
r1 = + rs[2] + MM_SET1(2)*rs[14] + MM_SET1(2)*rs[26] + rs[122] + MM_SET1(2)*rs[134] + rs[242]; GOUT_SCATTER(gout, n*3+2, r1);
}}
void int1e_r4_origi_ip2_optimizer(CINTOpt **opt, FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env) {
int ng[] = {4, 1, 0, 0, 5, 1, 1, 3};
CINTall_1e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
CACHE_SIZE_T int1e_r4_origi_ip2_cart(double *out, FINT *dims, FINT *shls,
FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache) {
int ng[] = {4, 1, 0, 0, 5, 1, 1, 3};
CINTEnvVars envs;
CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
envs.f_gout = &CINTgout1e_int1e_r4_origi_ip2;
return CINT1e_drv(out, dims, &envs, opt, cache, &c2s_cart_1e);
} // int1e_r4_origi_ip2_cart
CACHE_SIZE_T int1e_r4_origi_ip2_sph(double *out, FINT *dims, FINT *shls,
FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache) {
int ng[] = {4, 1, 0, 0, 5, 1, 1, 3};
CINTEnvVars envs;
CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
envs.f_gout = &CINTgout1e_int1e_r4_origi_ip2;
return CINT1e_drv(out, dims, &envs, opt, cache, &c2s_sph_1e);
} // int1e_r4_origi_ip2_sph
CACHE_SIZE_T int1e_r4_origi_ip2_spinor(double complex *out, FINT *dims, FINT *shls,
FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache) {
int ng[] = {4, 1, 0, 0, 5, 1, 1, 3};
CINTEnvVars envs;
CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
envs.f_gout = &CINTgout1e_int1e_r4_origi_ip2;
return CINT1e_spinor_drv(out, dims, &envs, opt, cache, &c2s_sf_1e);
} // int1e_r4_origi_ip2_spinor
ALL_CINT1E(int1e_r4_origi_ip2)
//ALL_CINT1E_FORTRAN_(cint1e_r4_origi_ip2)

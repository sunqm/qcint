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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "config.h"
#include "cint_bas.h"
#include "simd.h"
#include "misc.h"
#include "g2e.h"

void CINTinit_int3c2e_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                              int *atm, int natm, int *bas, int nbas, double *env)
{
        envs->natm = natm;
        envs->nbas = nbas;
        envs->atm = atm;
        envs->bas = bas;
        envs->env = env;
        envs->shls = shls;

        int i_sh = shls[0];
        int j_sh = shls[1];
        int k_sh = shls[2];
        envs->i_l = bas(ANG_OF, i_sh);
        envs->j_l = bas(ANG_OF, j_sh);
        envs->k_l = bas(ANG_OF, k_sh);
        envs->l_l = 0;
        envs->x_ctr[0] = bas(NCTR_OF, i_sh);
        envs->x_ctr[1] = bas(NCTR_OF, j_sh);
        envs->x_ctr[2] = bas(NCTR_OF, k_sh);
        envs->x_ctr[3] = 1;
        envs->nfi = (envs->i_l+1)*(envs->i_l+2)/2;
        envs->nfj = (envs->j_l+1)*(envs->j_l+2)/2;
        envs->nfk = (envs->k_l+1)*(envs->k_l+2)/2;
        envs->nfl = 1;
        envs->nf = envs->nfi * envs->nfk * envs->nfj;
        envs->common_factor = 1;
        if (env[PTR_EXPCUTOFF] == 0) {
                envs->expcutoff = EXPCUTOFF;
        } else {
                envs->expcutoff = MAX(MIN_EXPCUTOFF, env[PTR_EXPCUTOFF]);
        }

        envs->gbits = ng[GSHIFT];
        envs->ncomp_e1 = ng[POS_E1];
        envs->ncomp_e2 = ng[POS_E2];
        envs->ncomp_tensor = ng[TENSOR];

        envs->li_ceil = envs->i_l + ng[IINC];
        envs->lj_ceil = envs->j_l + ng[JINC];
        envs->lk_ceil = 0;
        envs->ll_ceil = envs->k_l + ng[KINC];

        envs->ri = env + atm(PTR_COORD, bas(ATOM_OF, i_sh));
        envs->rj = env + atm(PTR_COORD, bas(ATOM_OF, j_sh));
        envs->rk = env + atm(PTR_COORD, bas(ATOM_OF, k_sh));

        int nroots = (envs->li_ceil + envs->lj_ceil + envs->ll_ceil)/2 + 1;
        envs->nrys_roots = nroots;
        assert(nroots < MXRYSROOTS);

        int dli, dlj, dlk;
        int ibase = envs->li_ceil > envs->lj_ceil;
        if (envs->nrys_roots <= 2) { // use the fully optimized lj_4d algorithm
                ibase = 0;
        }

        if (ibase) {
                dli = envs->li_ceil + envs->lj_ceil + 1;
                dlj = envs->lj_ceil + 1;
        } else {
                dli = envs->li_ceil + 1;
                dlj = envs->li_ceil + envs->lj_ceil + 1;
        }
        dlk = envs->ll_ceil + 1;

        envs->g_stride_i = nroots;
        envs->g_stride_k = nroots * dli;
        envs->g_stride_l = nroots * dli;
        envs->g_stride_j = nroots * dli * dlk;
        envs->g_size     = nroots * dli * dlk * dlj;

        MM_STORE(envs->al, MM_SET1(0.));
        MM_STORE(envs->rkl+0*SIMDD, MM_SET1(envs->rk[0]));
        MM_STORE(envs->rkl+1*SIMDD, MM_SET1(envs->rk[1]));
        MM_STORE(envs->rkl+2*SIMDD, MM_SET1(envs->rk[2]));
        envs->g2d_klmax = envs->g_stride_k;
        envs->rkrl[0] = envs->rk[0];
        envs->rkrl[1] = envs->rk[1];
        envs->rkrl[2] = envs->rk[2];
        // in g0_2d rklrx = rkl - rx = 0 => rkl = rx
        envs->rx_in_rklrx = envs->rk;

        if (ibase) {
                envs->g2d_ijmax = envs->g_stride_i;
                envs->rx_in_rijrx = envs->ri;
                envs->rirj[0] = envs->ri[0] - envs->rj[0];
                envs->rirj[1] = envs->ri[1] - envs->rj[1];
                envs->rirj[2] = envs->ri[2] - envs->rj[2];
        } else {
                envs->g2d_ijmax = envs->g_stride_j;
                envs->rx_in_rijrx = envs->rj;
                envs->rirj[0] = envs->rj[0] - envs->ri[0];
                envs->rirj[1] = envs->rj[1] - envs->ri[1];
                envs->rirj[2] = envs->rj[2] - envs->ri[2];
        }

        if (ibase) {
                envs->f_g0_2d4d = &CINTg0_2e_il2d4d;
                envs->f_g0_2d4d_simd1 = &CINTg0_2e_il2d4d_simd1;
        } else {
                envs->f_g0_2d4d = &CINTg0_2e_lj2d4d;
                envs->f_g0_2d4d_simd1 = &CINTg0_2e_lj2d4d_simd1;
        }
        envs->f_g0_2e = &CINTg0_2e;
        envs->f_g0_2e_simd1 = &CINTg0_2e_simd1;
}


#ifdef WITH_GTG
void CINTg0_2e_lj2d4d_regular(double *g, Rys2eT *bc, CINTEnvVars *envs);
void CINTg0_2e_lj2d4d_simd1_regular(double *g, Rys2eT *bc, CINTEnvVars *envs);
int CINTg0_2e_gtg(double *g, Rys2eT *bc, CINTEnvVars *envs, int count);
int CINTg0_2e_gtg_simd1(double *g, Rys2eT *bc, CINTEnvVars *envs, int idsimd);

void CINTinit_int3c2e_gtg_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                                  int *atm, int natm, int *bas, int nbas, double *env)
{
        envs->natm = natm;
        envs->nbas = nbas;
        envs->atm = atm;
        envs->bas = bas;
        envs->env = env;
        envs->shls = shls;

        int i_sh = shls[0];
        int j_sh = shls[1];
        int k_sh = shls[2];
        envs->i_l = bas(ANG_OF, i_sh);
        envs->j_l = bas(ANG_OF, j_sh);
        envs->k_l = bas(ANG_OF, k_sh);
        envs->l_l = 0;
        envs->x_ctr[0] = bas(NCTR_OF, i_sh);
        envs->x_ctr[1] = bas(NCTR_OF, j_sh);
        envs->x_ctr[2] = bas(NCTR_OF, k_sh);
        envs->x_ctr[3] = 1;
        envs->nfi = (envs->i_l+1)*(envs->i_l+2)/2;
        envs->nfj = (envs->j_l+1)*(envs->j_l+2)/2;
        envs->nfk = (envs->k_l+1)*(envs->k_l+2)/2;
        envs->nfl = 1;
        envs->nf = envs->nfi * envs->nfk * envs->nfj;
        envs->common_factor = SQRTPI * .5;
        if (env[PTR_EXPCUTOFF] == 0) {
                envs->expcutoff = EXPCUTOFF;
        } else {
                // +1 to ensure accuracy. See comments in libcint/cint2e.c
                envs->expcutoff = MAX(MIN_EXPCUTOFF, env[PTR_EXPCUTOFF]);
        }

        envs->gbits = ng[GSHIFT];
        envs->ncomp_e1 = ng[POS_E1];
        envs->ncomp_e2 = ng[POS_E2];
        envs->ncomp_tensor = ng[TENSOR];

        envs->li_ceil = envs->i_l + ng[IINC];
        envs->lj_ceil = envs->j_l + ng[JINC];
        envs->lk_ceil = 0;
        envs->ll_ceil = envs->k_l + ng[KINC];

        envs->ri = env + atm(PTR_COORD, bas(ATOM_OF, i_sh));
        envs->rj = env + atm(PTR_COORD, bas(ATOM_OF, j_sh));
        envs->rk = env + atm(PTR_COORD, bas(ATOM_OF, k_sh));
        envs->nrys_roots = 1;

        int dli, dlj, dlk;
        int ibase = envs->li_ceil > envs->lj_ceil;

        if (ibase) {
                dli = envs->li_ceil + envs->lj_ceil + 1;
                dlj = envs->lj_ceil + 1;
        } else {
                dli = envs->li_ceil + 1;
                dlj = envs->li_ceil + envs->lj_ceil + 1;
        }
        dlk = envs->ll_ceil + 1;

        envs->g_stride_i = 1;
        envs->g_stride_k = dli;
        envs->g_stride_l = dli;
        envs->g_stride_j = dli * dlk;
        envs->g_size     = dli * dlk * dlj;

        MM_STORE(envs->al, MM_SET1(0.));
        MM_STORE(envs->rkl+0*SIMDD, MM_SET1(envs->rk[0]));
        MM_STORE(envs->rkl+1*SIMDD, MM_SET1(envs->rk[1]));
        MM_STORE(envs->rkl+2*SIMDD, MM_SET1(envs->rk[2]));
        envs->g2d_klmax = envs->g_stride_k;
        envs->rkrl[0] = envs->rk[0];
        envs->rkrl[1] = envs->rk[1];
        envs->rkrl[2] = envs->rk[2];
        // in g0_2d rklrx = rkl - rx = 0 => rkl = rx
        envs->rx_in_rklrx = envs->rk;

        if (ibase) {
                envs->g2d_ijmax = envs->g_stride_i;
                envs->rx_in_rijrx = envs->ri;
                envs->rirj[0] = envs->ri[0] - envs->rj[0];
                envs->rirj[1] = envs->ri[1] - envs->rj[1];
                envs->rirj[2] = envs->ri[2] - envs->rj[2];
        } else {
                envs->g2d_ijmax = envs->g_stride_j;
                envs->rx_in_rijrx = envs->rj;
                envs->rirj[0] = envs->rj[0] - envs->ri[0];
                envs->rirj[1] = envs->rj[1] - envs->ri[1];
                envs->rirj[2] = envs->rj[2] - envs->ri[2];
        }

        if (ibase) {
                envs->f_g0_2d4d = &CINTg0_2e_il2d4d;
                envs->f_g0_2d4d_simd1 = &CINTg0_2e_il2d4d_simd1;
        } else {
                envs->f_g0_2d4d = &CINTg0_2e_lj2d4d_regular;
                envs->f_g0_2d4d_simd1 = &CINTg0_2e_lj2d4d_simd1_regular;
        }
        envs->f_g0_2e = &CINTg0_2e_gtg;
        envs->f_g0_2e_simd1 = &CINTg0_2e_gtg_simd1;
}
#endif

/*
 * Qcint is a general GTO integral library for computational chemistry
 * Copyright (C) 2014 Qiming Sun <osirpt.sun@gmail.com>
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
#include <math.h>
#include <assert.h>
#include "cint_const.h"
#include "cint_bas.h"
#include "g2e.h"

struct _BC {
        double c00[MXRYSROOTS*3];
        double c0p[MXRYSROOTS*3];
        double b01[MXRYSROOTS];
        double b00[MXRYSROOTS];
        double b10[MXRYSROOTS];
};

void CINTg0_2e_lj2d4d(double *g, const CINTEnvVars *envs,struct _BC *bc);
void CINTg0_2e_kj2d4d(double *g, const CINTEnvVars *envs,struct _BC *bc);
void CINTg0_2e_il2d4d(double *g, const CINTEnvVars *envs,struct _BC *bc);
void CINTg0_2e_ik2d4d(double *g, const CINTEnvVars *envs,struct _BC *bc);
static void CINTset_g2e_params(CINTEnvVars *envs);

int CINTinit_int2e_EnvVars(CINTEnvVars *envs, const int *ng, const int *shls,
                           const int *atm, const int natm,
                           const int *bas, const int nbas, const double *env)
{
        envs->natm = natm;
        envs->nbas = nbas;
        envs->atm = atm;
        envs->bas = bas;
        envs->env = env;
        envs->shls = shls;

        const int i_sh = shls[0];
        const int j_sh = shls[1];
        const int k_sh = shls[2];
        const int l_sh = shls[3];
        envs->i_l = bas(ANG_OF, i_sh);
        envs->j_l = bas(ANG_OF, j_sh);
        envs->k_l = bas(ANG_OF, k_sh);
        envs->l_l = bas(ANG_OF, l_sh);
        envs->i_prim = bas(NPRIM_OF, i_sh);
        envs->j_prim = bas(NPRIM_OF, j_sh);
        envs->k_prim = bas(NPRIM_OF, k_sh);
        envs->l_prim = bas(NPRIM_OF, l_sh);
        envs->i_ctr = bas(NCTR_OF, i_sh);
        envs->j_ctr = bas(NCTR_OF, j_sh);
        envs->k_ctr = bas(NCTR_OF, k_sh);
        envs->l_ctr = bas(NCTR_OF, l_sh);
        envs->nfi = CINTlen_cart(envs->i_l);
        envs->nfj = CINTlen_cart(envs->j_l);
        envs->nfk = CINTlen_cart(envs->k_l);
        envs->nfl = CINTlen_cart(envs->l_l);
        envs->nf = envs->nfi * envs->nfk * envs->nfl * envs->nfj;

        envs->ri = env + atm(PTR_COORD, bas(ATOM_OF, i_sh));
        envs->rj = env + atm(PTR_COORD, bas(ATOM_OF, j_sh));
        envs->rk = env + atm(PTR_COORD, bas(ATOM_OF, k_sh));
        envs->rl = env + atm(PTR_COORD, bas(ATOM_OF, l_sh));

        envs->common_factor = (M_PI*M_PI*M_PI)*2/SQRTPI
                * CINTcommon_fac_sp(envs->i_l) * CINTcommon_fac_sp(envs->j_l)
                * CINTcommon_fac_sp(envs->k_l) * CINTcommon_fac_sp(envs->l_l);

        envs->gbits = ng[GSHIFT];
        envs->ncomp_e1 = ng[POS_E1];
        envs->ncomp_e2 = ng[POS_E2];
        envs->ncomp_tensor = ng[TENSOR];

        envs->li_ceil = envs->i_l + ng[IINC]; // ng[0] has different meaning in cint1e
        envs->lj_ceil = envs->j_l + ng[JINC];
        envs->lk_ceil = envs->k_l + ng[KINC];
        envs->ll_ceil = envs->l_l + ng[LINC];
        envs->nrys_roots =(envs->li_ceil + envs->lj_ceil
                         + envs->lk_ceil + envs->ll_ceil)/2 + 1;

        assert(i_sh < SHLS_MAX);
        assert(j_sh < SHLS_MAX);
        assert(k_sh < SHLS_MAX);
        assert(l_sh < SHLS_MAX);
        assert(envs->i_l < ANG_MAX);
        assert(envs->j_l < ANG_MAX);
        assert(envs->k_l < ANG_MAX);
        assert(envs->l_l < ANG_MAX);
        assert(envs->i_ctr < NCTR_MAX);
        assert(envs->j_ctr < NCTR_MAX);
        assert(envs->k_ctr < NCTR_MAX);
        assert(envs->l_ctr < NCTR_MAX);
        assert(envs->i_prim < NPRIM_MAX);
        assert(envs->j_prim < NPRIM_MAX);
        assert(envs->k_prim < NPRIM_MAX);
        assert(envs->l_prim < NPRIM_MAX);
        assert(envs->i_prim >= envs->i_ctr);
        assert(envs->j_prim >= envs->j_ctr);
        assert(envs->k_prim >= envs->k_ctr);
        assert(envs->l_prim >= envs->l_ctr);
        assert(bas(ATOM_OF,i_sh) >= 0);
        assert(bas(ATOM_OF,j_sh) >= 0);
        assert(bas(ATOM_OF,k_sh) >= 0);
        assert(bas(ATOM_OF,l_sh) >= 0);
        assert(bas(ATOM_OF,i_sh) < natm);
        assert(bas(ATOM_OF,j_sh) < natm);
        assert(bas(ATOM_OF,k_sh) < natm);
        assert(bas(ATOM_OF,l_sh) < natm);
        assert(envs->nrys_roots < MXRYSROOTS);

        CINTset_g2e_params(envs);
        return 0;
}

/* set strides and parameters for g0_2d4d algorithm */
static void CINTset_g2e_params(CINTEnvVars *envs)
{
        int dli, dlj, dlk, dll;
        int ibase = envs->li_ceil > envs->lj_ceil;
        int kbase = envs->lk_ceil > envs->ll_ceil;
        if (envs->nrys_roots <= 2) { // use the fully optimized lj_4d algorithm
                ibase = 0;
                kbase = 0;
        }
        if (kbase) {
                dlk = envs->lk_ceil + envs->ll_ceil + 1;
                dll = envs->ll_ceil + 1;
        } else {
                dlk = envs->lk_ceil + 1;
                dll = envs->lk_ceil + envs->ll_ceil + 1;
        }

        if (ibase) {
                dli = envs->li_ceil + envs->lj_ceil + 1;
                dlj = envs->lj_ceil + 1;
        } else {
                dli = envs->li_ceil + 1;
                dlj = envs->li_ceil + envs->lj_ceil + 1;
        }
        if (envs->nrys_roots > 2 && envs->nrys_roots%2) {
                envs->g_stride_i = (envs->nrys_roots+1);
                envs->g_stride_k = (envs->nrys_roots+1) * dli;
                envs->g_stride_l = (envs->nrys_roots+1) * dli * dlk;
                envs->g_stride_j = (envs->nrys_roots+1) * dli * dlk * dll;
                envs->g_size     = (envs->nrys_roots+1) * dli * dlk * dll * dlj;
        } else {
                envs->g_stride_i = envs->nrys_roots;
                envs->g_stride_k = envs->nrys_roots * dli;
                envs->g_stride_l = envs->nrys_roots * dli * dlk;
                envs->g_stride_j = envs->nrys_roots * dli * dlk * dll;
                envs->g_size     = envs->nrys_roots * dli * dlk * dll * dlj;
        }

        if (kbase) {
                envs->g2d_klmax = envs->g_stride_k;
                envs->rx_in_rklrx = envs->rk;
                envs->rkrl[0] = envs->rk[0] - envs->rl[0];
                envs->rkrl[1] = envs->rk[1] - envs->rl[1];
                envs->rkrl[2] = envs->rk[2] - envs->rl[2];
        } else {
                envs->g2d_klmax = envs->g_stride_l;
                envs->rx_in_rklrx = envs->rl;
                envs->rkrl[0] = envs->rl[0] - envs->rk[0];
                envs->rkrl[1] = envs->rl[1] - envs->rk[1];
                envs->rkrl[2] = envs->rl[2] - envs->rk[2];
        }

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

        if (kbase) {
                if (ibase) {
                        envs->f_g0_2d4d = &CINTg0_2e_ik2d4d;
                } else {
                        envs->f_g0_2d4d = &CINTg0_2e_kj2d4d;
                }
        } else {
                if (ibase) {
                        envs->f_g0_2d4d = &CINTg0_2e_il2d4d;
                } else {
                        envs->f_g0_2d4d = &CINTg0_2e_lj2d4d;
                }
        }
}

void CINTg2e_index_xyz(int *idx, const CINTEnvVars *envs)
{
        const int i_l = envs->i_l;
        const int j_l = envs->j_l;
        const int k_l = envs->k_l;
        const int l_l = envs->l_l;
        const int nfi = envs->nfi;
        const int nfj = envs->nfj;
        const int nfk = envs->nfk;
        const int nfl = envs->nfl;
        const int di = envs->g_stride_i;
        const int dk = envs->g_stride_k;
        const int dl = envs->g_stride_l;
        const int dj = envs->g_stride_j;
        int i, j, k, l, n;
        int ofx, ofkx, oflx;
        int ofy, ofky, ofly;
        int ofz, ofkz, oflz;
        int i_nx[CART_MAX], i_ny[CART_MAX], i_nz[CART_MAX];
        int j_nx[CART_MAX], j_ny[CART_MAX], j_nz[CART_MAX];
        int k_nx[CART_MAX], k_ny[CART_MAX], k_nz[CART_MAX];
        int l_nx[CART_MAX], l_ny[CART_MAX], l_nz[CART_MAX];

        CINTcart_comp(i_nx, i_ny, i_nz, i_l);
        CINTcart_comp(j_nx, j_ny, j_nz, j_l);
        CINTcart_comp(k_nx, k_ny, k_nz, k_l);
        CINTcart_comp(l_nx, l_ny, l_nz, l_l);

        ofx = 0;
        ofy = envs->g_size;
        ofz = envs->g_size * 2;
        n = 0;
        for (j = 0; j < nfj; j++) {
                for (l = 0; l < nfl; l++) {
                        oflx = ofx + dj * j_nx[j] + dl * l_nx[l];
                        ofly = ofy + dj * j_ny[j] + dl * l_ny[l];
                        oflz = ofz + dj * j_nz[j] + dl * l_nz[l];
                        for (k = 0; k < nfk; k++) {
                                ofkx = oflx + dk * k_nx[k];
                                ofky = ofly + dk * k_ny[k];
                                ofkz = oflz + dk * k_nz[k];
                                switch (i_l) {
                                        case 0:
                                                idx[n+0] = ofkx;
                                                idx[n+1] = ofky;
                                                idx[n+2] = ofkz;
                                                n += 3;
                                                break;
                                        case 1:
                                                idx[n+0] = ofkx + di;
                                                idx[n+1] = ofky;
                                                idx[n+2] = ofkz;
                                                idx[n+3] = ofkx;
                                                idx[n+4] = ofky + di;
                                                idx[n+5] = ofkz;
                                                idx[n+6] = ofkx;
                                                idx[n+7] = ofky;
                                                idx[n+8] = ofkz + di;
                                                n += 9;
                                                break;
                                        case 2:
                                                idx[n+0 ] = ofkx + di*2;
                                                idx[n+1 ] = ofky;
                                                idx[n+2 ] = ofkz;
                                                idx[n+3 ] = ofkx + di;
                                                idx[n+4 ] = ofky + di;
                                                idx[n+5 ] = ofkz;
                                                idx[n+6 ] = ofkx + di;
                                                idx[n+7 ] = ofky;
                                                idx[n+8 ] = ofkz + di;
                                                idx[n+9 ] = ofkx;
                                                idx[n+10] = ofky + di*2;
                                                idx[n+11] = ofkz;
                                                idx[n+12] = ofkx;
                                                idx[n+13] = ofky + di;
                                                idx[n+14] = ofkz + di;
                                                idx[n+15] = ofkx;
                                                idx[n+16] = ofky;
                                                idx[n+17] = ofkz + di*2;
                                                n += 18;
                                                break;
                                        default:
                                                for (i = 0; i < nfi; i++) {
                                                        idx[n+0] = ofkx + di * i_nx[i]; //(:,ix,kx,lx,jx,1)
                                                        idx[n+1] = ofky + di * i_ny[i]; //(:,iy,ky,ly,jy,2)
                                                        idx[n+2] = ofkz + di * i_nz[i]; //(:,iz,kz,lz,jz,3)
                                                        n += 3;
                                                } // i
                                }
                        } // k
                } // l
        } // j
}


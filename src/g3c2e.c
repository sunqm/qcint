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
#include "cint_const.h"
#include "cint_bas.h"
#include "simd.h"
#include "rys_roots.h"
#include "g2e.h"

#define DEF_GXYZ(type, G, GX, GY, GZ) \
        type *GX = G; \
        type *GY = G + envs->g_size     * SIMDD; \
        type *GZ = G + envs->g_size * 2 * SIMDD

void CINTg0_3c2e_kj2d3d(double *g, Rys2eT *bc, CINTEnvVars *envs);
void CINTg0_3c2e_ik2d3d(double *g, Rys2eT *bc, CINTEnvVars *envs);
void CINTg0_3c2e_kj2d3d_simd1(double *g, Rys2eT *bc, CINTEnvVars *envs);
void CINTg0_3c2e_ik2d3d_simd1(double *g, Rys2eT *bc, CINTEnvVars *envs);


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
        envs->x_ctr[0] = bas(NCTR_OF, i_sh);
        envs->x_ctr[1] = bas(NCTR_OF, j_sh);
        envs->x_ctr[2] = bas(NCTR_OF, k_sh);
        envs->x_ctr[3] = 1;
        envs->nfi = (envs->i_l+1)*(envs->i_l+2)/2;
        envs->nfj = (envs->j_l+1)*(envs->j_l+2)/2;
        envs->nfk = (envs->k_l+1)*(envs->k_l+2)/2;
        envs->nf = envs->nfi * envs->nfk * envs->nfj;
        envs->common_factor = 1;

        envs->gbits = ng[GSHIFT];
        envs->ncomp_e1 = ng[POS_E1];
        envs->ncomp_e2 = ng[POS_E2];
        envs->ncomp_tensor = ng[TENSOR];

        envs->li_ceil = envs->i_l + ng[IINC];
        envs->lj_ceil = envs->j_l + ng[JINC];
        envs->lk_ceil = envs->k_l + ng[KINC];
        envs->ll_ceil = 0;

        envs->ri = env + atm(PTR_COORD, bas(ATOM_OF, i_sh));
        envs->rj = env + atm(PTR_COORD, bas(ATOM_OF, j_sh));
        envs->rk = env + atm(PTR_COORD, bas(ATOM_OF, k_sh));

        int nroots = (envs->li_ceil + envs->lj_ceil + envs->lk_ceil)/2 + 1;
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
        dlk = envs->lk_ceil + 1;

        envs->g_stride_i = nroots;
        envs->g_stride_k = nroots * dli;
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
                envs->f_g0_2d4d = &CINTg0_3c2e_ik2d3d;
                envs->f_g0_2d4d_simd1 = &CINTg0_3c2e_ik2d3d_simd1;
        } else {
                envs->f_g0_2d4d = &CINTg0_3c2e_kj2d3d;
                envs->f_g0_2d4d_simd1 = &CINTg0_3c2e_kj2d3d_simd1;
        }
}

void CINTg3c_index_xyz(int *idx, CINTEnvVars *envs)
{
        int i_l = envs->i_l;
        int j_l = envs->j_l;
        int k_l = envs->k_l;
        int nfi = envs->nfi;
        int nfj = envs->nfj;
        int nfk = envs->nfk;
        int di = envs->g_stride_i;
        int dk = envs->g_stride_k;
        int dj = envs->g_stride_j;
        int i, j, k, n;
        int ofx, ofkx;
        int ofy, ofky;
        int ofz, ofkz;
        int i_nx[CART_MAX], i_ny[CART_MAX], i_nz[CART_MAX];
        int j_nx[CART_MAX], j_ny[CART_MAX], j_nz[CART_MAX];
        int k_nx[CART_MAX], k_ny[CART_MAX], k_nz[CART_MAX];

        CINTcart_comp(i_nx, i_ny, i_nz, i_l);
        CINTcart_comp(j_nx, j_ny, j_nz, j_l);
        CINTcart_comp(k_nx, k_ny, k_nz, k_l);

        ofx = 0;
        ofy = envs->g_size;
        ofz = envs->g_size * 2;
        n = 0;
        for (j = 0; j < nfj; j++) {
                for (k = 0; k < nfk; k++) {
                        ofkx = ofx + dj * j_nx[j] + dk * k_nx[k];
                        ofky = ofy + dj * j_ny[j] + dk * k_ny[k];
                        ofkz = ofz + dj * j_nz[j] + dk * k_nz[k];
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
                                                idx[n+0] = ofkx + di * i_nx[i]; //(:,ix,kx,jx,1)
                                                idx[n+1] = ofky + di * i_ny[i]; //(:,iy,ky,jy,2)
                                                idx[n+2] = ofkz + di * i_nz[i]; //(:,iz,kz,jz,3)
                                                n += 3;
                                        } // i
                        }
                } // k
        } // j
}

/* 2d is based on i,k */
void CINTg0_3c2e_ik2d3d(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        CINTg0_2e_2d(g, bc, envs);
        int lj = envs->lj_ceil;
        if (lj == 0) {
                return;
        }
        //int li = envs->li_ceil;
        int lk = envs->lk_ceil;
        int j, k, ptr, n;
        int di = envs->g_stride_i;
        int dk = envs->g_stride_k;
        int dj = envs->g_stride_j;
        double *rirj = envs->rirj;
        DEF_GXYZ(double, g, gx, gy, gz);
        double *p1x, *p1y, *p1z, *p2x, *p2y, *p2z;
        __MD rx = MM_SET1(rirj[0]);
        __MD ry = MM_SET1(rirj[1]);
        __MD rz = MM_SET1(rirj[2]);
        // g(i,...,j) = rirj * g(i,...,j-1) +  g(i+1,...,j-1)
        p1x = gx  - dj * SIMDD;
        p1y = gy  - dj * SIMDD;
        p1z = gz  - dj * SIMDD;
        p2x = p1x + di * SIMDD;
        p2y = p1y + di * SIMDD;
        p2z = p1z + di * SIMDD;
        for (j = 1; j <= lj; j++) {
        for (k = 0; k <= lk; k++) {
                ptr = j*dj + k*dk;
                for (n = ptr; n < ptr+dk-di*j; n++) {
MM_STORE(gx+n*SIMDD, MM_FMA(rx, MM_LOAD(p1x+n*SIMDD), MM_LOAD(p2x+n*SIMDD)));
MM_STORE(gy+n*SIMDD, MM_FMA(ry, MM_LOAD(p1y+n*SIMDD), MM_LOAD(p2y+n*SIMDD)));
MM_STORE(gz+n*SIMDD, MM_FMA(rz, MM_LOAD(p1z+n*SIMDD), MM_LOAD(p2z+n*SIMDD)));
                }
        } }
}

void CINTg0_3c2e_ik2d3d_simd1(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        CINTg0_2e_2d_simd1(g, bc, envs);
        int lj = envs->lj_ceil;
        if (lj == 0) {
                return;
        }
        //const int li = envs->li_ceil;
        const int lk = envs->lk_ceil;
        int j, k, ptr, n;
        const int di = envs->g_stride_i;
        const int dk = envs->g_stride_k;
        const int dj = envs->g_stride_j;
        const double *rirj = envs->rirj;
        double *gx = g;
        double *gy = g + envs->g_size;
        double *gz = g + envs->g_size * 2;
        const double *p1x, *p1y, *p1z, *p2x, *p2y, *p2z;

        // g(i,...,j) = rirj * g(i,...,j-1) +  g(i+1,...,j-1)
        p1x = gx  - dj;
        p1y = gy  - dj;
        p1z = gz  - dj;
        p2x = p1x + di;
        p2y = p1y + di;
        p2z = p1z + di;
        for (j = 1; j <= lj; j++) {
        for (k = 0; k <= lk; k++) {
                ptr = j*dj + k*dk;
                for (n = ptr; n < ptr+dk-di*j; n++) {
                        gx[n] = rirj[0] * p1x[n] + p2x[n];
                        gy[n] = rirj[1] * p1y[n] + p2y[n];
                        gz[n] = rirj[2] * p1z[n] + p2z[n];
                }
        } }
}

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
#include "misc.h"
#include "simd.h"
#include "rys_roots.h"
#include "g2e.h"

#define DEF_GXYZ(type, G, GX, GY, GZ) \
        type *RESTRICT GX = G; \
        type *RESTRICT GY = G + envs->g_size     * SIMDD; \
        type *RESTRICT GZ = G + envs->g_size * 2 * SIMDD


void CINTinit_int2e_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
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
        int l_sh = shls[3];
        envs->i_l = bas(ANG_OF, i_sh);
        envs->j_l = bas(ANG_OF, j_sh);
        envs->k_l = bas(ANG_OF, k_sh);
        envs->l_l = bas(ANG_OF, l_sh);
        envs->x_ctr[0] = bas(NCTR_OF, i_sh);
        envs->x_ctr[1] = bas(NCTR_OF, j_sh);
        envs->x_ctr[2] = bas(NCTR_OF, k_sh);
        envs->x_ctr[3] = bas(NCTR_OF, l_sh);
        envs->nfi = (envs->i_l+1)*(envs->i_l+2)/2;
        envs->nfj = (envs->j_l+1)*(envs->j_l+2)/2;
        envs->nfk = (envs->k_l+1)*(envs->k_l+2)/2;
        envs->nfl = (envs->l_l+1)*(envs->l_l+2)/2;
        envs->nf = envs->nfi * envs->nfk * envs->nfl * envs->nfj;
        envs->common_factor = 1;
        if (env[PTR_EXPCUTOFF] == 0) {
                envs->expcutoff = EXPCUTOFF;
        } else {
                // +1 to ensure accuracy. See comments in libcint/cint2e.c
                envs->expcutoff = MAX(MIN_EXPCUTOFF, env[PTR_EXPCUTOFF]) + 1;
        }

        envs->gbits = ng[GSHIFT];
        envs->ncomp_e1 = ng[POS_E1];
        envs->ncomp_e2 = ng[POS_E2];
        envs->ncomp_tensor = ng[TENSOR];

        envs->li_ceil = envs->i_l + ng[IINC];
        envs->lj_ceil = envs->j_l + ng[JINC];
        envs->lk_ceil = envs->k_l + ng[KINC];
        envs->ll_ceil = envs->l_l + ng[LINC];

        envs->ri = env + atm(PTR_COORD, bas(ATOM_OF, i_sh));
        envs->rj = env + atm(PTR_COORD, bas(ATOM_OF, j_sh));
        envs->rk = env + atm(PTR_COORD, bas(ATOM_OF, k_sh));
        envs->rl = env + atm(PTR_COORD, bas(ATOM_OF, l_sh));

        int nroots = (envs->li_ceil + envs->lj_ceil +
                      envs->lk_ceil + envs->ll_ceil)/2 + 1;
        envs->nrys_roots = nroots;
        assert(nroots < MXRYSROOTS);

        int dli, dlj, dlk, dll;
        int ibase = envs->li_ceil > envs->lj_ceil;
        int kbase = envs->lk_ceil > envs->ll_ceil;
        if (nroots <= 2) {
                envs->f_g0_2d4d = &CINTg0_2e_2d4d_unrolled;
                envs->f_g0_2d4d_simd1 = &CINTg0_2e_2d4d_unrolled_simd1;
        } else if (kbase) {
                if (ibase) {
                        envs->f_g0_2d4d = &CINTg0_2e_ik2d4d;
                        envs->f_g0_2d4d_simd1 = &CINTg0_2e_ik2d4d_simd1;
                } else {
                        envs->f_g0_2d4d = &CINTg0_2e_kj2d4d;
                        envs->f_g0_2d4d_simd1 = &CINTg0_2e_kj2d4d_simd1;
                }
        } else {
                if (ibase) {
                        envs->f_g0_2d4d = &CINTg0_2e_il2d4d;
                        envs->f_g0_2d4d_simd1 = &CINTg0_2e_il2d4d_simd1;
                } else {
                        envs->f_g0_2d4d = &CINTg0_2e_lj2d4d;
                        envs->f_g0_2d4d_simd1 = &CINTg0_2e_lj2d4d_simd1;
                }
        }
        envs->f_g0_2e = &CINTg0_2e;
        envs->f_g0_2e_simd1 = &CINTg0_2e_simd1;

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
        envs->g_stride_i = nroots;
        envs->g_stride_k = nroots * dli;
        envs->g_stride_l = nroots * dli * dlk;
        envs->g_stride_j = nroots * dli * dlk * dll;
        envs->g_size     = nroots * dli * dlk * dll * dlj;

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
}

void CINTg4c_index_xyz(int *idx, CINTEnvVars *envs)
{
        int i_l = envs->i_l;
        int j_l = envs->j_l;
        int k_l = envs->k_l;
        int l_l = envs->l_l;
        int nfi = envs->nfi;
        int nfj = envs->nfj;
        int nfk = envs->nfk;
        int nfl = envs->nfl;
        int di = envs->g_stride_i;
        int dk = envs->g_stride_k;
        int dl = envs->g_stride_l;
        int dj = envs->g_stride_j;
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
                                                idx[n+0] = ofkx + di * i_nx[i];
                                                idx[n+1] = ofky + di * i_ny[i];
                                                idx[n+2] = ofkz + di * i_nz[i];
                                                n += 3;
                                        } // i
                                }
                        } // k
                } // l
        } // j
}


void CINTg0_2e_2d(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        int nroots = envs->nrys_roots;
        int i, m, n;
        DEF_GXYZ(double, g, gx, gy, gz);
        int nmax = envs->li_ceil + envs->lj_ceil;
        int mmax = envs->lk_ceil + envs->ll_ceil;
        int dm = envs->g2d_klmax;
        int dn = envs->g2d_ijmax;
        double *RESTRICT c00x = bc->c00x;
        double *RESTRICT c00y = bc->c00y;
        double *RESTRICT c00z = bc->c00z;
        double *RESTRICT c0px = bc->c0px;
        double *RESTRICT c0py = bc->c0py;
        double *RESTRICT c0pz = bc->c0pz;
        double *RESTRICT b01 = bc->b01;
        double *RESTRICT b00 = bc->b00;
        double *RESTRICT b10 = bc->b10;
        double *RESTRICT p0x, *RESTRICT p0y, *RESTRICT p0z;
        double *RESTRICT p1x, *RESTRICT p1y, *RESTRICT p1z;
        __MD r0, r1;
        __MD r0x, r0y, r0z;
        __MD r1x, r1y, r1z;
        __MD r2x, r2y, r2z;
        __MD rcx, rcy, rcz;

        if (nmax > 0) {
                p0x = gx + dn * SIMDD;
                p0y = gy + dn * SIMDD;
                p0z = gz + dn * SIMDD;
                p1x = gx - dn * SIMDD;
                p1y = gy - dn * SIMDD;
                p1z = gz - dn * SIMDD;
                for (i = 0; i < nroots; i++) {
                        r2x = MM_LOAD(gx+i*SIMDD);
                        r2y = MM_LOAD(gy+i*SIMDD);
                        r2z = MM_LOAD(gz+i*SIMDD);
                        rcx = MM_LOAD(c00x+i*SIMDD);
                        rcy = MM_LOAD(c00y+i*SIMDD);
                        rcz = MM_LOAD(c00z+i*SIMDD);
                        r0x = rcx * r2x;
                        r0y = rcy * r2y;
                        r0z = rcz * r2z;
                        MM_STORE(p0x+i*SIMDD, r0x);
                        MM_STORE(p0y+i*SIMDD, r0y);
                        MM_STORE(p0z+i*SIMDD, r0z);
                        for (n = 1; n < nmax; n++) {
                                r0 = MM_MUL(MM_SET1(n), MM_LOAD(b10+i*SIMDD));
                                r1x = rcx * r0x + r0 * r2x;
                                r1y = rcy * r0y + r0 * r2y;
                                r1z = rcz * r0z + r0 * r2z;
                                MM_STORE(p0x+(i+n*dn)*SIMDD, r1x);
                                MM_STORE(p0y+(i+n*dn)*SIMDD, r1y);
                                MM_STORE(p0z+(i+n*dn)*SIMDD, r1z);
                                r2x = r0x;
                                r2y = r0y;
                                r2z = r0z;
                                r0x = r1x;
                                r0y = r1y;
                                r0z = r1z;
                        }
                }
        }

        if (mmax > 0) {
                p0x = gx + dm * SIMDD;
                p0y = gy + dm * SIMDD;
                p0z = gz + dm * SIMDD;
                for (i = 0; i < nroots; i++) {
                        r2x = MM_LOAD(gx+i*SIMDD);
                        r2y = MM_LOAD(gy+i*SIMDD);
                        r2z = MM_LOAD(gz+i*SIMDD);
                        rcx = MM_LOAD(c0px+i*SIMDD);
                        rcy = MM_LOAD(c0py+i*SIMDD);
                        rcz = MM_LOAD(c0pz+i*SIMDD);
                        r0x = rcx * r2x;
                        r0y = rcy * r2y;
                        r0z = rcz * r2z;
                        MM_STORE(p0x+i*SIMDD, r0x);
                        MM_STORE(p0y+i*SIMDD, r0y);
                        MM_STORE(p0z+i*SIMDD, r0z);
                        for (m = 1; m < mmax; m++) {
                                r0 = MM_MUL(MM_SET1(m), MM_LOAD(b01+i*SIMDD));
                                r1x = rcx * r0x + r0 * r2x;
                                r1y = rcy * r0y + r0 * r2y;
                                r1z = rcz * r0z + r0 * r2z;
                                MM_STORE(p0x+(i+m*dm)*SIMDD, r1x);
                                MM_STORE(p0y+(i+m*dm)*SIMDD, r1y);
                                MM_STORE(p0z+(i+m*dm)*SIMDD, r1z);
                                r2x = r0x;
                                r2y = r0y;
                                r2z = r0z;
                                r0x = r1x;
                                r0y = r1y;
                                r0z = r1z;
                        }
                }
        }

        if (nmax > 0 && mmax > 0) {
                p0x = gx + (dm + dn) * SIMDD;
                p0y = gy + (dm + dn) * SIMDD;
                p0z = gz + (dm + dn) * SIMDD;
                p1x = gx + dn * SIMDD;
                p1y = gy + dn * SIMDD;
                p1z = gz + dn * SIMDD;
                for (i = 0; i < nroots; i++) {
                        r0 = MM_LOAD(b00+i*SIMDD);
                        r1x = MM_LOAD(c0px+i*SIMDD) * MM_LOAD(gx+(i+dn)*SIMDD) + r0 * MM_LOAD(gx+i*SIMDD);
                        r1y = MM_LOAD(c0py+i*SIMDD) * MM_LOAD(gy+(i+dn)*SIMDD) + r0 * MM_LOAD(gy+i*SIMDD);
                        r1z = MM_LOAD(c0pz+i*SIMDD) * MM_LOAD(gz+(i+dn)*SIMDD) + r0 * MM_LOAD(gz+i*SIMDD);
                        MM_STORE(p0x+i*SIMDD, r1x);
                        MM_STORE(p0y+i*SIMDD, r1y);
                        MM_STORE(p0z+i*SIMDD, r1z);
                        r2x = MM_LOAD(gx+(i+dm)*SIMDD);
                        r2y = MM_LOAD(gy+(i+dm)*SIMDD);
                        r2z = MM_LOAD(gz+(i+dm)*SIMDD);
                        rcx = MM_LOAD(c00x+i*SIMDD);
                        rcy = MM_LOAD(c00y+i*SIMDD);
                        rcz = MM_LOAD(c00z+i*SIMDD);
                        for (n = 1; n < nmax; n++) {
                                r1 = MM_SET1(n) * MM_LOAD(b10+i*SIMDD);
                                r0x = rcx * r1x + r1 * r2x + r0 * MM_LOAD(gx+(i+n*dn)*SIMDD);
                                r0y = rcy * r1y + r1 * r2y + r0 * MM_LOAD(gy+(i+n*dn)*SIMDD);
                                r0z = rcz * r1z + r1 * r2z + r0 * MM_LOAD(gz+(i+n*dn)*SIMDD);
                                MM_STORE(p0x+(i+n*dn)*SIMDD, r0x);
                                MM_STORE(p0y+(i+n*dn)*SIMDD, r0y);
                                MM_STORE(p0z+(i+n*dn)*SIMDD, r0z);
                                r2x = r1x;
                                r2y = r1y;
                                r2z = r1z;
                                r1x = r0x;
                                r1y = r0y;
                                r1z = r0z;
                        }

                }

                for (m = 1; m < mmax; m++) {
                for (i = 0; i < nroots; i++) {
                        r0 = MM_LOAD(b00+i*SIMDD);
                        r1 = MM_SET1(m) * MM_LOAD(b01+i*SIMDD);
                        r1x = MM_LOAD(c0px+i*SIMDD) * MM_LOAD(p1x+(i+m*dm)*SIMDD) + r1 * MM_LOAD(p1x+(i+(m-1)*dm)*SIMDD) + r0 * MM_LOAD(gx+(i+m*dm)*SIMDD);
                        r1y = MM_LOAD(c0py+i*SIMDD) * MM_LOAD(p1y+(i+m*dm)*SIMDD) + r1 * MM_LOAD(p1y+(i+(m-1)*dm)*SIMDD) + r0 * MM_LOAD(gy+(i+m*dm)*SIMDD);
                        r1z = MM_LOAD(c0pz+i*SIMDD) * MM_LOAD(p1z+(i+m*dm)*SIMDD) + r1 * MM_LOAD(p1z+(i+(m-1)*dm)*SIMDD) + r0 * MM_LOAD(gz+(i+m*dm)*SIMDD);
                        MM_STORE(p0x+(i+m*dm)*SIMDD, r1x);
                        MM_STORE(p0y+(i+m*dm)*SIMDD, r1y);
                        MM_STORE(p0z+(i+m*dm)*SIMDD, r1z);

                        r2x = MM_LOAD(gx+(i+(m+1)*dm)*SIMDD);
                        r2y = MM_LOAD(gy+(i+(m+1)*dm)*SIMDD);
                        r2z = MM_LOAD(gz+(i+(m+1)*dm)*SIMDD);
                        rcx = MM_LOAD(c00x+i*SIMDD);
                        rcy = MM_LOAD(c00y+i*SIMDD);
                        rcz = MM_LOAD(c00z+i*SIMDD);
                        r0 = MM_SET1(m + 1) * r0;
                        for (n = 1; n < nmax; n++) {
                                r1 = MM_SET1(n) * MM_LOAD(b10+i*SIMDD);
                                r0x = rcx * r1x + r1 * r2x + r0 * MM_LOAD(gx+(i+m*dm+n*dn)*SIMDD);
                                r0y = rcy * r1y + r1 * r2y + r0 * MM_LOAD(gy+(i+m*dm+n*dn)*SIMDD);
                                r0z = rcz * r1z + r1 * r2z + r0 * MM_LOAD(gz+(i+m*dm+n*dn)*SIMDD);
                                MM_STORE(p0x+(i+m*dm+n*dn)*SIMDD, r0x);
                                MM_STORE(p0y+(i+m*dm+n*dn)*SIMDD, r0y);
                                MM_STORE(p0z+(i+m*dm+n*dn)*SIMDD, r0z);
                                r2x = r1x;
                                r2y = r1y;
                                r2z = r1z;
                                r1x = r0x;
                                r1y = r0y;
                                r1z = r0z;
                        }
                } }
        }
}


/*
 * g0[i,k,l,j] = < ik | lj > = ( i j | k l )
 */
/* 2d is based on l,j */
void CINTg0_lj_4d(double *g, CINTEnvVars *envs)
{
        int li = envs->li_ceil;
        int lk = envs->lk_ceil;
        if (li == 0 && lk == 0) {
                return;
        }
        int nmax = envs->li_ceil + envs->lj_ceil;
        int mmax = envs->lk_ceil + envs->ll_ceil;
        //int ll = envs->ll_ceil;
        int lj = envs->lj_ceil;
        int nroots = envs->nrys_roots;
        int i, j, k, l, ptr, n;
        int di = envs->g_stride_i;
        int dk = envs->g_stride_k;
        int dl = envs->g_stride_l;
        int dj = envs->g_stride_j;
        double *rirj = envs->rirj;
        double *rkrl = envs->rkrl;
        DEF_GXYZ(double, g, gx, gy, gz);
        double *RESTRICT p1x;
        double *RESTRICT p1y;
        double *RESTRICT p1z;
        double *RESTRICT p2x;
        double *RESTRICT p2y;
        double *RESTRICT p2z;
        __MD rx = MM_SET1(rirj[0]);
        __MD ry = MM_SET1(rirj[1]);
        __MD rz = MM_SET1(rirj[2]);

        // g(i,...,j) = rirj * g(i-1,...,j) +  g(i-1,...,j+1)
        p1x = gx  - di * SIMDD;
        p1y = gy  - di * SIMDD;
        p1z = gz  - di * SIMDD;
        p2x = p1x + dj * SIMDD;
        p2y = p1y + dj * SIMDD;
        p2z = p1z + dj * SIMDD;
        for (i = 1; i <= li; i++) {
        for (j = 0; j <= nmax-i; j++) {
        for (l = 0; l <= mmax; l++) {
                ptr = j*dj + l*dl + i*di;
                for (n = ptr; n < ptr+nroots; n++) {
//for (s = 0; s < SIMDD; s++) {
//        gx[n*SIMDD+s] = rirj[0] * p1x[n*SIMDD+s] + p2x[n*SIMDD+s];
//        gy[n*SIMDD+s] = rirj[1] * p1y[n*SIMDD+s] + p2y[n*SIMDD+s];
//        gz[n*SIMDD+s] = rirj[2] * p1z[n*SIMDD+s] + p2z[n*SIMDD+s];
//}
MM_STORE(gx+n*SIMDD, MM_FMA(rx, MM_LOAD(p1x+n*SIMDD), MM_LOAD(p2x+n*SIMDD)));
MM_STORE(gy+n*SIMDD, MM_FMA(ry, MM_LOAD(p1y+n*SIMDD), MM_LOAD(p2y+n*SIMDD)));
MM_STORE(gz+n*SIMDD, MM_FMA(rz, MM_LOAD(p1z+n*SIMDD), MM_LOAD(p2z+n*SIMDD)));
                }
        } } }

        rx = MM_SET1(rkrl[0]);
        ry = MM_SET1(rkrl[1]);
        rz = MM_SET1(rkrl[2]);
        // g(...,k,l,..) = rkrl * g(...,k-1,l,..) + g(...,k-1,l+1,..)
        p1x = gx  - dk * SIMDD;
        p1y = gy  - dk * SIMDD;
        p1z = gz  - dk * SIMDD;
        p2x = p1x + dl * SIMDD;
        p2y = p1y + dl * SIMDD;
        p2z = p1z + dl * SIMDD;
        for (j = 0; j <= lj; j++) {
        for (k = 1; k <= lk; k++) {
        for (l = 0; l <= mmax-k; l++) {
                ptr = j*dj + l*dl + k*dk;
                for (n = ptr; n < ptr+dk; n++) {
//for (s = 0; s < SIMDD; s++) {
//        gx[n*SIMDD+s] = rkrl[0] * p1x[n*SIMDD+s] + p2x[n*SIMDD+s];
//        gy[n*SIMDD+s] = rkrl[1] * p1y[n*SIMDD+s] + p2y[n*SIMDD+s];
//        gz[n*SIMDD+s] = rkrl[2] * p1z[n*SIMDD+s] + p2z[n*SIMDD+s];
//}
MM_STORE(gx+n*SIMDD, MM_FMA(rx, MM_LOAD(p1x+n*SIMDD), MM_LOAD(p2x+n*SIMDD)));
MM_STORE(gy+n*SIMDD, MM_FMA(ry, MM_LOAD(p1y+n*SIMDD), MM_LOAD(p2y+n*SIMDD)));
MM_STORE(gz+n*SIMDD, MM_FMA(rz, MM_LOAD(p1z+n*SIMDD), MM_LOAD(p2z+n*SIMDD)));
                }
        } } }
}
/* 2d is based on k,j */
void CINTg0_kj_4d(double *g, CINTEnvVars *envs)
{
        int li = envs->li_ceil;
        int ll = envs->ll_ceil;
        if (li == 0 && ll == 0) {
                return;
        }
        int nmax = envs->li_ceil + envs->lj_ceil;
        int mmax = envs->lk_ceil + envs->ll_ceil;
        //int lk = envs->lk_ceil;
        int lj = envs->lj_ceil;
        int nroots = envs->nrys_roots;
        int i, j, k, l, ptr, n;
        int di = envs->g_stride_i;
        int dk = envs->g_stride_k;
        int dl = envs->g_stride_l;
        int dj = envs->g_stride_j;
        double *rirj = envs->rirj;
        double *rkrl = envs->rkrl;
        DEF_GXYZ(double, g, gx, gy, gz);
        double *RESTRICT p1x;
        double *RESTRICT p1y;
        double *RESTRICT p1z;
        double *RESTRICT p2x;
        double *RESTRICT p2y;
        double *RESTRICT p2z;
        __MD rx = MM_SET1(rirj[0]);
        __MD ry = MM_SET1(rirj[1]);
        __MD rz = MM_SET1(rirj[2]);

        // g(i,...,j) = rirj * g(i-1,...,j) +  g(i-1,...,j+1)
        p1x = gx  - di * SIMDD;
        p1y = gy  - di * SIMDD;
        p1z = gz  - di * SIMDD;
        p2x = p1x + dj * SIMDD;
        p2y = p1y + dj * SIMDD;
        p2z = p1z + dj * SIMDD;
        for (i = 1; i <= li; i++) {
        for (j = 0; j <= nmax-i; j++) {
        for (k = 0; k <= mmax; k++) {
                ptr = j*dj + k*dk + i*di;
                for (n = ptr; n < ptr+nroots; n++) {
//for (s = 0; s < SIMDD; s++) {
//        gx[n*SIMDD+s] = rirj[0] * p1x[n*SIMDD+s] + p2x[n*SIMDD+s];
//        gy[n*SIMDD+s] = rirj[1] * p1y[n*SIMDD+s] + p2y[n*SIMDD+s];
//        gz[n*SIMDD+s] = rirj[2] * p1z[n*SIMDD+s] + p2z[n*SIMDD+s];
//}
MM_STORE(gx+n*SIMDD, MM_FMA(rx, MM_LOAD(p1x+n*SIMDD), MM_LOAD(p2x+n*SIMDD)));
MM_STORE(gy+n*SIMDD, MM_FMA(ry, MM_LOAD(p1y+n*SIMDD), MM_LOAD(p2y+n*SIMDD)));
MM_STORE(gz+n*SIMDD, MM_FMA(rz, MM_LOAD(p1z+n*SIMDD), MM_LOAD(p2z+n*SIMDD)));
                }
        } } }

        rx = MM_SET1(rkrl[0]);
        ry = MM_SET1(rkrl[1]);
        rz = MM_SET1(rkrl[2]);
        // g(...,k,l,..) = rkrl * g(...,k,l-1,..) + g(...,k+1,l-1,..)
        p1x = gx  - dl * SIMDD;
        p1y = gy  - dl * SIMDD;
        p1z = gz  - dl * SIMDD;
        p2x = p1x + dk * SIMDD;
        p2y = p1y + dk * SIMDD;
        p2z = p1z + dk * SIMDD;
        for (j = 0; j <= lj; j++) {
        for (l = 1; l <= ll; l++) {
        for (k = 0; k <= mmax-l; k++) {
                ptr = j*dj + l*dl + k*dk;
                for (n = ptr; n < ptr+dk; n++) {
//for (s = 0; s < SIMDD; s++) {
//        gx[n*SIMDD+s] = rkrl[0] * p1x[n*SIMDD+s] + p2x[n*SIMDD+s];
//        gy[n*SIMDD+s] = rkrl[1] * p1y[n*SIMDD+s] + p2y[n*SIMDD+s];
//        gz[n*SIMDD+s] = rkrl[2] * p1z[n*SIMDD+s] + p2z[n*SIMDD+s];
//}
MM_STORE(gx+n*SIMDD, MM_FMA(rx, MM_LOAD(p1x+n*SIMDD), MM_LOAD(p2x+n*SIMDD)));
MM_STORE(gy+n*SIMDD, MM_FMA(ry, MM_LOAD(p1y+n*SIMDD), MM_LOAD(p2y+n*SIMDD)));
MM_STORE(gz+n*SIMDD, MM_FMA(rz, MM_LOAD(p1z+n*SIMDD), MM_LOAD(p2z+n*SIMDD)));
                }
        } } }
}
/* 2d is based on i,l */
void CINTg0_il_4d(double *g, CINTEnvVars *envs)
{
        int lj = envs->lj_ceil;
        int lk = envs->lk_ceil;
        if (lj == 0 && lk == 0) {
                return;
        }
        int nmax = envs->li_ceil + envs->lj_ceil;
        int mmax = envs->lk_ceil + envs->ll_ceil;
        //int li = envs->li_ceil;
        int ll = envs->ll_ceil;
        int nroots = envs->nrys_roots;
        int i, j, k, l, ptr, n;
        int di = envs->g_stride_i;
        int dk = envs->g_stride_k;
        int dl = envs->g_stride_l;
        int dj = envs->g_stride_j;
        double *rirj = envs->rirj;
        double *rkrl = envs->rkrl;
        DEF_GXYZ(double, g, gx, gy, gz);
        double *RESTRICT p1x;
        double *RESTRICT p1y;
        double *RESTRICT p1z;
        double *RESTRICT p2x;
        double *RESTRICT p2y;
        double *RESTRICT p2z;
        __MD rx = MM_SET1(rkrl[0]);
        __MD ry = MM_SET1(rkrl[1]);
        __MD rz = MM_SET1(rkrl[2]);

        // g(...,k,l,..) = rkrl * g(...,k-1,l,..) + g(...,k-1,l+1,..)
        p1x = gx  - dk * SIMDD;
        p1y = gy  - dk * SIMDD;
        p1z = gz  - dk * SIMDD;
        p2x = p1x + dl * SIMDD;
        p2y = p1y + dl * SIMDD;
        p2z = p1z + dl * SIMDD;
        for (k = 1; k <= lk; k++) {
        for (l = 0; l <= mmax-k; l++) {
        for (i = 0; i <= nmax; i++) {
                ptr = l*dl + k*dk + i*di;
                for (n = ptr; n < ptr+nroots; n++) {
//for (s = 0; s < SIMDD; s++) {
//        gx[n*SIMDD+s] = rkrl[0] * p1x[n*SIMDD+s] + p2x[n*SIMDD+s];
//        gy[n*SIMDD+s] = rkrl[1] * p1y[n*SIMDD+s] + p2y[n*SIMDD+s];
//        gz[n*SIMDD+s] = rkrl[2] * p1z[n*SIMDD+s] + p2z[n*SIMDD+s];
//}
MM_STORE(gx+n*SIMDD, MM_FMA(rx, MM_LOAD(p1x+n*SIMDD), MM_LOAD(p2x+n*SIMDD)));
MM_STORE(gy+n*SIMDD, MM_FMA(ry, MM_LOAD(p1y+n*SIMDD), MM_LOAD(p2y+n*SIMDD)));
MM_STORE(gz+n*SIMDD, MM_FMA(rz, MM_LOAD(p1z+n*SIMDD), MM_LOAD(p2z+n*SIMDD)));
                }
        } } }

        rx = MM_SET1(rirj[0]);
        ry = MM_SET1(rirj[1]);
        rz = MM_SET1(rirj[2]);
        // g(i,...,j) = rirj * g(i,...,j-1) +  g(i+1,...,j-1)
        p1x = gx  - dj * SIMDD;
        p1y = gy  - dj * SIMDD;
        p1z = gz  - dj * SIMDD;
        p2x = p1x + di * SIMDD;
        p2y = p1y + di * SIMDD;
        p2z = p1z + di * SIMDD;
        for (j = 1; j <= lj; j++) {
        for (l = 0; l <= ll; l++) {
        for (k = 0; k <= lk; k++) {
                ptr = j*dj + l*dl + k*dk;
                for (n = ptr; n < ptr+dk-di*j; n++) {
//for (s = 0; s < SIMDD; s++) {
//        gx[n*SIMDD+s] = rirj[0] * p1x[n*SIMDD+s] + p2x[n*SIMDD+s];
//        gy[n*SIMDD+s] = rirj[1] * p1y[n*SIMDD+s] + p2y[n*SIMDD+s];
//        gz[n*SIMDD+s] = rirj[2] * p1z[n*SIMDD+s] + p2z[n*SIMDD+s];
//}
MM_STORE(gx+n*SIMDD, MM_FMA(rx, MM_LOAD(p1x+n*SIMDD), MM_LOAD(p2x+n*SIMDD)));
MM_STORE(gy+n*SIMDD, MM_FMA(ry, MM_LOAD(p1y+n*SIMDD), MM_LOAD(p2y+n*SIMDD)));
MM_STORE(gz+n*SIMDD, MM_FMA(rz, MM_LOAD(p1z+n*SIMDD), MM_LOAD(p2z+n*SIMDD)));
                }
        } } }
}
/* 2d is based on i,k */
void CINTg0_ik_4d(double *g, CINTEnvVars *envs)
{
        int lj = envs->lj_ceil;
        int ll = envs->ll_ceil;
        if (lj == 0 && ll == 0) {
                return;
        }
        int nmax = envs->li_ceil + envs->lj_ceil;
        int mmax = envs->lk_ceil + envs->ll_ceil;
        //int li = envs->li_ceil;
        int lk = envs->lk_ceil;
        int nroots = envs->nrys_roots;
        int i, j, k, l, ptr, n;
        int di = envs->g_stride_i;
        int dk = envs->g_stride_k;
        int dl = envs->g_stride_l;
        int dj = envs->g_stride_j;
        double *rirj = envs->rirj;
        double *rkrl = envs->rkrl;
        DEF_GXYZ(double, g, gx, gy, gz);
        double *RESTRICT p1x;
        double *RESTRICT p1y;
        double *RESTRICT p1z;
        double *RESTRICT p2x;
        double *RESTRICT p2y;
        double *RESTRICT p2z;
        __MD rx = MM_SET1(rkrl[0]);
        __MD ry = MM_SET1(rkrl[1]);
        __MD rz = MM_SET1(rkrl[2]);

        // g(...,k,l,..) = rkrl * g(...,k,l-1,..) + g(...,k+1,l-1,..)
        p1x = gx  - dl * SIMDD;
        p1y = gy  - dl * SIMDD;
        p1z = gz  - dl * SIMDD;
        p2x = p1x + dk * SIMDD;
        p2y = p1y + dk * SIMDD;
        p2z = p1z + dk * SIMDD;
        for (l = 1; l <= ll; l++) {
                // (:,i) is full, so loop:k and loop:n can be merged to
                // for(n = l*dl; n < ptr+dl-dk*l; n++)
                for (k = 0; k <= mmax-l; k++) {
                for (i = 0; i <= nmax; i++) {
                        ptr = l*dl + k*dk + i*di;
                        for (n = ptr; n < ptr+nroots; n++) {
//for (s = 0; s < SIMDD; s++) {
//        gx[n*SIMDD+s] = rkrl[0] * p1x[n*SIMDD+s] + p2x[n*SIMDD+s];
//        gy[n*SIMDD+s] = rkrl[1] * p1y[n*SIMDD+s] + p2y[n*SIMDD+s];
//        gz[n*SIMDD+s] = rkrl[2] * p1z[n*SIMDD+s] + p2z[n*SIMDD+s];
//}
MM_STORE(gx+n*SIMDD, MM_FMA(rx, MM_LOAD(p1x+n*SIMDD), MM_LOAD(p2x+n*SIMDD)));
MM_STORE(gy+n*SIMDD, MM_FMA(ry, MM_LOAD(p1y+n*SIMDD), MM_LOAD(p2y+n*SIMDD)));
MM_STORE(gz+n*SIMDD, MM_FMA(rz, MM_LOAD(p1z+n*SIMDD), MM_LOAD(p2z+n*SIMDD)));
                        }
                } }
        }

        rx = MM_SET1(rirj[0]);
        ry = MM_SET1(rirj[1]);
        rz = MM_SET1(rirj[2]);
        // g(i,...,j) = rirj * g(i,...,j-1) +  g(i+1,...,j-1)
        p1x = gx  - dj * SIMDD;
        p1y = gy  - dj * SIMDD;
        p1z = gz  - dj * SIMDD;
        p2x = p1x + di * SIMDD;
        p2y = p1y + di * SIMDD;
        p2z = p1z + di * SIMDD;
        for (j = 1; j <= lj; j++) {
        for (l = 0; l <= ll; l++) {
        for (k = 0; k <= lk; k++) {
                ptr = j*dj + l*dl + k*dk;
                for (n = ptr; n < ptr+dk-di*j; n++) {
//for (s = 0; s < SIMDD; s++) {
//        gx[n*SIMDD+s] = rirj[0] * p1x[n*SIMDD+s] + p2x[n*SIMDD+s];
//        gy[n*SIMDD+s] = rirj[1] * p1y[n*SIMDD+s] + p2y[n*SIMDD+s];
//        gz[n*SIMDD+s] = rirj[2] * p1z[n*SIMDD+s] + p2z[n*SIMDD+s];
//}
MM_STORE(gx+n*SIMDD, MM_FMA(rx, MM_LOAD(p1x+n*SIMDD), MM_LOAD(p2x+n*SIMDD)));
MM_STORE(gy+n*SIMDD, MM_FMA(ry, MM_LOAD(p1y+n*SIMDD), MM_LOAD(p2y+n*SIMDD)));
MM_STORE(gz+n*SIMDD, MM_FMA(rz, MM_LOAD(p1z+n*SIMDD), MM_LOAD(p2z+n*SIMDD)));
                }
        } } }
}
/************* some special g0_4d results *************/
/* 4 digits stand for i_ceil, k_ceil, l_ceil, j_ceil */
static inline void _make_rc(double *rc, double *cx, double *cy, double *cz, double *r)
{
        __MD r0 = MM_SET1(r[0]);
        __MD r1 = MM_SET1(r[1]);
        __MD r2 = MM_SET1(r[2]);
        MM_STORE(rc+0*SIMDD, MM_ADD(r0, MM_LOAD(cx+0*SIMDD)));
        MM_STORE(rc+1*SIMDD, MM_ADD(r0, MM_LOAD(cx+1*SIMDD)));
        MM_STORE(rc+2*SIMDD, MM_ADD(r1, MM_LOAD(cy+0*SIMDD)));
        MM_STORE(rc+3*SIMDD, MM_ADD(r1, MM_LOAD(cy+1*SIMDD)));
        MM_STORE(rc+4*SIMDD, MM_ADD(r2, MM_LOAD(cz+0*SIMDD)));
        MM_STORE(rc+5*SIMDD, MM_ADD(r2, MM_LOAD(cz+1*SIMDD)));
}

static inline void _g0_2d4d_0000(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        //MM_STORE(g+0*SIMDD, MM_SET1(1.));
        //MM_STORE(g+1*SIMDD, MM_SET1(1.));
        //g[2] = w[0];
}

static inline void _g0_2d4d_0001(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *cx = bc->c0px;
        double *cy = bc->c0py;
        double *cz = bc->c0pz;
        MM_STORE(g+1*SIMDD, MM_LOAD(cx+0*SIMDD));
        MM_STORE(g+3*SIMDD, MM_LOAD(cy+0*SIMDD));
        MM_STORE(g+5*SIMDD, MM_MUL(MM_LOAD(cz+0*SIMDD), MM_LOAD(g+4*SIMDD)));
}

static inline void _g0_2d4d_0002(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *cx = bc->c0px;
        double *cy = bc->c0py;
        double *cz = bc->c0pz;
        double *b =  bc->b01;
        __MD cx0 = MM_LOAD(cx+0*SIMDD);
        __MD cx1 = MM_LOAD(cx+1*SIMDD);
        __MD cy0 = MM_LOAD(cy+0*SIMDD);
        __MD cy1 = MM_LOAD(cy+1*SIMDD);
        __MD cz0 = MM_LOAD(cz+0*SIMDD);
        __MD cz1 = MM_LOAD(cz+1*SIMDD);
        __MD b0 = MM_LOAD(b+0*SIMDD);
        __MD b1 = MM_LOAD(b+1*SIMDD);
        __MD g12 = MM_LOAD(g+12*SIMDD);
        __MD g13 = MM_LOAD(g+13*SIMDD);
        MM_STORE(g+2 *SIMDD, cx0);
        MM_STORE(g+3 *SIMDD, cx1);
        MM_STORE(g+8 *SIMDD, cy0);
        MM_STORE(g+9 *SIMDD, cy1);
        MM_STORE(g+14*SIMDD, MM_MUL(cz0, g12));
        MM_STORE(g+15*SIMDD, MM_MUL(cz1, g13));
        MM_STORE(g+4 *SIMDD, MM_FMA(cx0, cx0, b0));
        MM_STORE(g+5 *SIMDD, MM_FMA(cx1, cx1, b1));
        MM_STORE(g+10*SIMDD, MM_FMA(cy0, cy0, b0));
        MM_STORE(g+11*SIMDD, MM_FMA(cy1, cy1, b1));
        MM_STORE(g+16*SIMDD, MM_MUL(MM_FMA(cz0, cz0, b0), g12));
        MM_STORE(g+17*SIMDD, MM_MUL(MM_FMA(cz1, cz1, b1), g13));
}

static inline void _g0_2d4d_0003(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *cx = bc->c0px;
        double *cy = bc->c0py;
        double *cz = bc->c0pz;
        double *b =  bc->b01;
        __MD b0 = MM_LOAD(b+0*SIMDD);
        __MD b1 = MM_LOAD(b+1*SIMDD);
        __MD g16 = MM_LOAD(g+16*SIMDD);
        __MD g17 = MM_LOAD(g+17*SIMDD);
        __MD i3 = MM_SET1(3.);
        MM_STORE(g+2 *SIMDD, MM_LOAD(cx+0*SIMDD));
        MM_STORE(g+3 *SIMDD, MM_LOAD(cx+1*SIMDD));
        MM_STORE(g+10*SIMDD, MM_LOAD(cy+0*SIMDD));
        MM_STORE(g+11*SIMDD, MM_LOAD(cy+1*SIMDD));
        MM_STORE(g+18*SIMDD, MM_LOAD(cz+0*SIMDD) * g16);
        MM_STORE(g+19*SIMDD, MM_LOAD(cz+1*SIMDD) * g17);
        MM_STORE(g+4 *SIMDD, MM_LOAD(cx+0*SIMDD) * MM_LOAD(cx+0*SIMDD) + b0);
        MM_STORE(g+5 *SIMDD, MM_LOAD(cx+1*SIMDD) * MM_LOAD(cx+1*SIMDD) + b1);
        MM_STORE(g+12*SIMDD, MM_LOAD(cy+0*SIMDD) * MM_LOAD(cy+0*SIMDD) + b0);
        MM_STORE(g+13*SIMDD, MM_LOAD(cy+1*SIMDD) * MM_LOAD(cy+1*SIMDD) + b1);
        MM_STORE(g+20*SIMDD,(MM_LOAD(cz+0*SIMDD) * MM_LOAD(cz+0*SIMDD) + b0)* g16);
        MM_STORE(g+21*SIMDD,(MM_LOAD(cz+1*SIMDD) * MM_LOAD(cz+1*SIMDD) + b1)* g17);
        MM_STORE(g+6 *SIMDD, MM_LOAD(cx+0*SIMDD) *(MM_LOAD(cx+0*SIMDD) * MM_LOAD(cx+0*SIMDD) + i3 * b0));
        MM_STORE(g+7 *SIMDD, MM_LOAD(cx+1*SIMDD) *(MM_LOAD(cx+1*SIMDD) * MM_LOAD(cx+1*SIMDD) + i3 * b1));
        MM_STORE(g+14*SIMDD, MM_LOAD(cy+0*SIMDD) *(MM_LOAD(cy+0*SIMDD) * MM_LOAD(cy+0*SIMDD) + i3 * b0));
        MM_STORE(g+15*SIMDD, MM_LOAD(cy+1*SIMDD) *(MM_LOAD(cy+1*SIMDD) * MM_LOAD(cy+1*SIMDD) + i3 * b1));
        MM_STORE(g+22*SIMDD,(MM_LOAD(cz+0*SIMDD) * MM_LOAD(cz+0*SIMDD) + i3 * b0)* MM_LOAD(g+18*SIMDD));
        MM_STORE(g+23*SIMDD,(MM_LOAD(cz+1*SIMDD) * MM_LOAD(cz+1*SIMDD) + i3 * b1)* MM_LOAD(g+19*SIMDD));
}

static inline void _g0_2d4d_0010(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *cx = bc->c0px;
        double *cy = bc->c0py;
        double *cz = bc->c0pz;
        MM_STORE(g+1*SIMDD, MM_LOAD(cx+0*SIMDD));
        MM_STORE(g+3*SIMDD, MM_LOAD(cy+0*SIMDD));
        MM_STORE(g+5*SIMDD, MM_MUL(MM_LOAD(cz+0*SIMDD), MM_LOAD(g+4*SIMDD)));
}

static inline void _g0_2d4d_0011(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *cx = bc->c0px;
        double *cy = bc->c0py;
        double *cz = bc->c0pz;
        double *b  = bc->b01;
        double *r  = envs->rkrl;
        ALIGNMM double rc[6*SIMDD];
        _make_rc(rc, cx, cy, cz, r);
        __MD r0 = MM_LOAD(rc+0*SIMDD);
        __MD r1 = MM_LOAD(rc+1*SIMDD);
        __MD r2 = MM_LOAD(rc+2*SIMDD);
        __MD r3 = MM_LOAD(rc+3*SIMDD);
        __MD r4 = MM_LOAD(rc+4*SIMDD);
        __MD r5 = MM_LOAD(rc+5*SIMDD);
        __MD cx0 = MM_LOAD(cx+0*SIMDD);
        __MD cx1 = MM_LOAD(cx+1*SIMDD);
        __MD cy0 = MM_LOAD(cy+0*SIMDD);
        __MD cy1 = MM_LOAD(cy+1*SIMDD);
        __MD cz0 = MM_LOAD(cz+0*SIMDD);
        __MD cz1 = MM_LOAD(cz+1*SIMDD);
        __MD b0 = MM_LOAD(b+0*SIMDD);
        __MD b1 = MM_LOAD(b+1*SIMDD);
        __MD g24 = MM_LOAD(g+24*SIMDD);
        __MD g25 = MM_LOAD(g+25*SIMDD);
        MM_STORE(g+2 *SIMDD, r0);
        MM_STORE(g+3 *SIMDD, r1);
        MM_STORE(g+4 *SIMDD, cx0);
        MM_STORE(g+5 *SIMDD, cx1);
        MM_STORE(g+14*SIMDD, r2);
        MM_STORE(g+15*SIMDD, r3);
        MM_STORE(g+16*SIMDD, cy0);
        MM_STORE(g+17*SIMDD, cy1);
        MM_STORE(g+26*SIMDD, MM_MUL(r4, g24));
        MM_STORE(g+27*SIMDD, MM_MUL(r5, g25));
        MM_STORE(g+28*SIMDD, MM_MUL(cz0,g24));
        MM_STORE(g+29*SIMDD, MM_MUL(cz1,g25));
        MM_STORE(g+6 *SIMDD, MM_FMA(r0, cx0, b0));
        MM_STORE(g+7 *SIMDD, MM_FMA(r1, cx1, b1));
        MM_STORE(g+18*SIMDD, MM_FMA(r2, cy0, b0));
        MM_STORE(g+19*SIMDD, MM_FMA(r3, cy1, b1));
        MM_STORE(g+30*SIMDD, MM_MUL(MM_FMA(r4, cz0, b0), g24));
        MM_STORE(g+31*SIMDD, MM_MUL(MM_FMA(r5, cz1, b1), g25));
}

static inline void _g0_2d4d_0012(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *cx = bc->c0px;
        double *cy = bc->c0py;
        double *cz = bc->c0pz;
        double *b =  bc->b01;
        double *r  = envs->rkrl;
        ALIGNMM double rc[6*SIMDD];
        _make_rc(rc, cx, cy, cz, r);
        __MD r0 = MM_LOAD(rc+0*SIMDD);
        __MD r1 = MM_LOAD(rc+1*SIMDD);
        __MD r2 = MM_LOAD(rc+2*SIMDD);
        __MD r3 = MM_LOAD(rc+3*SIMDD);
        __MD r4 = MM_LOAD(rc+4*SIMDD);
        __MD r5 = MM_LOAD(rc+5*SIMDD);
        __MD cx0 = MM_LOAD(cx+0*SIMDD);
        __MD cx1 = MM_LOAD(cx+1*SIMDD);
        __MD cy0 = MM_LOAD(cy+0*SIMDD);
        __MD cy1 = MM_LOAD(cy+1*SIMDD);
        __MD cz0 = MM_LOAD(cz+0*SIMDD);
        __MD cz1 = MM_LOAD(cz+1*SIMDD);
        __MD b0 = MM_LOAD(b+0*SIMDD);
        __MD b1 = MM_LOAD(b+1*SIMDD);
        __MD i2 = MM_SET1(2.);
        __MD g32 = MM_LOAD(g+32*SIMDD);
        __MD g33 = MM_LOAD(g+33*SIMDD);
        MM_STORE(g+2 *SIMDD, r0 );
        MM_STORE(g+3 *SIMDD, r1 );
        MM_STORE(g+4 *SIMDD, cx0);
        MM_STORE(g+5 *SIMDD, cx1);
        MM_STORE(g+18*SIMDD, r2 );
        MM_STORE(g+19*SIMDD, r3 );
        MM_STORE(g+20*SIMDD, cy0);
        MM_STORE(g+21*SIMDD, cy1);
        MM_STORE(g+34*SIMDD, MM_MUL(r4 , g32));
        MM_STORE(g+35*SIMDD, MM_MUL(r5 , g33));
        MM_STORE(g+36*SIMDD, MM_MUL(cz0, g32));
        MM_STORE(g+37*SIMDD, MM_MUL(cz1, g33));
        MM_STORE(g+6 *SIMDD, MM_FMA(r0 , cx0, b0));
        MM_STORE(g+7 *SIMDD, MM_FMA(r1 , cx1, b1));
        MM_STORE(g+8 *SIMDD, MM_FMA(cx0, cx0, b0));
        MM_STORE(g+9 *SIMDD, MM_FMA(cx1, cx1, b1));
        MM_STORE(g+22*SIMDD, MM_FMA(r2 , cy0, b0));
        MM_STORE(g+23*SIMDD, MM_FMA(r3 , cy1, b1));
        MM_STORE(g+24*SIMDD, MM_FMA(cy0, cy0, b0));
        MM_STORE(g+25*SIMDD, MM_FMA(cy1, cy1, b1));
        MM_STORE(g+38*SIMDD, MM_MUL(MM_FMA(r4 , cz0, b0), g32));
        MM_STORE(g+39*SIMDD, MM_MUL(MM_FMA(r5 , cz1, b1), g33));
        MM_STORE(g+40*SIMDD, MM_MUL(MM_FMA(cz0, cz0, b0), g32));
        MM_STORE(g+41*SIMDD, MM_MUL(MM_FMA(cz1, cz1, b1), g33));
        MM_STORE(g+10*SIMDD, r0 * MM_LOAD(g+8 *SIMDD) + i2 * b0 * cx0                );
        MM_STORE(g+11*SIMDD, r1 * MM_LOAD(g+9 *SIMDD) + i2 * b1 * cx1                );
        MM_STORE(g+26*SIMDD, r2 * MM_LOAD(g+24*SIMDD) + i2 * b0 * cy0                );
        MM_STORE(g+27*SIMDD, r3 * MM_LOAD(g+25*SIMDD) + i2 * b1 * cy1                );
        MM_STORE(g+42*SIMDD, r4 * MM_LOAD(g+40*SIMDD) + i2 * b0 * MM_LOAD(g+36*SIMDD));
        MM_STORE(g+43*SIMDD, r5 * MM_LOAD(g+41*SIMDD) + i2 * b1 * MM_LOAD(g+37*SIMDD));
}

static inline void _g0_2d4d_0020(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *cx = bc->c0px;
        double *cy = bc->c0py;
        double *cz = bc->c0pz;
        double *b  = bc->b01;
        __MD r0 = MM_LOAD(cx+0*SIMDD);
        __MD r1 = MM_LOAD(cx+1*SIMDD);
        __MD r2 = MM_LOAD(cy+0*SIMDD);
        __MD r3 = MM_LOAD(cy+1*SIMDD);
        __MD r4 = MM_LOAD(cz+0*SIMDD);
        __MD r5 = MM_LOAD(cz+1*SIMDD);
        __MD b0 = MM_LOAD(b+0*SIMDD);
        __MD b1 = MM_LOAD(b+1*SIMDD);
        __MD g12 = MM_LOAD(g+12*SIMDD);
        __MD g13 = MM_LOAD(g+13*SIMDD);
        MM_STORE(g+2 *SIMDD, r0);
        MM_STORE(g+3 *SIMDD, r1);
        MM_STORE(g+8 *SIMDD, r2);
        MM_STORE(g+9 *SIMDD, r3);
        MM_STORE(g+14*SIMDD, MM_MUL(r4, g12));
        MM_STORE(g+15*SIMDD, MM_MUL(r5, g13));
        MM_STORE(g+4 *SIMDD, MM_FMA(r0, r0, b0));
        MM_STORE(g+5 *SIMDD, MM_FMA(r1, r1, b1));
        MM_STORE(g+10*SIMDD, MM_FMA(r2, r2, b0));
        MM_STORE(g+11*SIMDD, MM_FMA(r3, r3, b1));
        MM_STORE(g+16*SIMDD, MM_MUL(MM_FMA(r4, r4, b0), g12));
        MM_STORE(g+17*SIMDD, MM_MUL(MM_FMA(r5, r5, b1), g13));
}

static inline void _g0_2d4d_0021(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *cx = bc->c0px;
        double *cy = bc->c0py;
        double *cz = bc->c0pz;
        double *b1 = bc->b01;
        double *rkl = envs->rkrl;
        ALIGNMM double rc[6*SIMDD];
        _make_rc(rc, cx, cy, cz, rkl);
        __MD b10 = MM_LOAD(b1+0*SIMDD);
        __MD b11 = MM_LOAD(b1+1*SIMDD);
        __MD g32 = MM_LOAD(g+32*SIMDD);
        __MD g33 = MM_LOAD(g+33*SIMDD);
        __MD i2 = MM_SET1(2.);
        __MD r0 = MM_LOAD(rc+0*SIMDD);
        __MD r1 = MM_LOAD(rc+1*SIMDD);
        __MD r2 = MM_LOAD(rc+2*SIMDD);
        __MD r3 = MM_LOAD(rc+3*SIMDD);
        __MD r4 = MM_LOAD(rc+4*SIMDD);
        __MD r5 = MM_LOAD(rc+5*SIMDD);
        __MD s0 = MM_LOAD(cx+0*SIMDD);
        __MD s1 = MM_LOAD(cx+1*SIMDD);
        __MD s2 = MM_LOAD(cy+0*SIMDD);
        __MD s3 = MM_LOAD(cy+1*SIMDD);
        __MD s4 = MM_LOAD(cz+0*SIMDD);
        __MD s5 = MM_LOAD(cz+1*SIMDD);
        MM_STORE(g+2 *SIMDD, s0);
        MM_STORE(g+3 *SIMDD, s1);
        MM_STORE(g+8 *SIMDD, r0);
        MM_STORE(g+9 *SIMDD, r1);
        MM_STORE(g+18*SIMDD, s2);
        MM_STORE(g+19*SIMDD, s3);
        MM_STORE(g+24*SIMDD, r2);
        MM_STORE(g+25*SIMDD, r3);
        MM_STORE(g+34*SIMDD, s4 * g32);
        MM_STORE(g+35*SIMDD, s5 * g33);
        MM_STORE(g+40*SIMDD, r4 * g32);
        MM_STORE(g+41*SIMDD, r5 * g33);
        MM_STORE(g+4 *SIMDD, s0 * s0 + b10);
        MM_STORE(g+5 *SIMDD, s1 * s1 + b11);
        MM_STORE(g+10*SIMDD, s0 * r0 + b10);
        MM_STORE(g+11*SIMDD, s1 * r1 + b11);
        MM_STORE(g+20*SIMDD, s2 * s2 + b10);
        MM_STORE(g+21*SIMDD, s3 * s3 + b11);
        MM_STORE(g+26*SIMDD, s2 * r2 + b10);
        MM_STORE(g+27*SIMDD, s3 * r3 + b11);
        MM_STORE(g+36*SIMDD,(s4 * s4 + b10) * g32);
        MM_STORE(g+37*SIMDD,(s5 * s5 + b11) * g33);
        MM_STORE(g+42*SIMDD,(s4 * r4 + b10) * g32);
        MM_STORE(g+43*SIMDD,(s5 * r5 + b11) * g33);
        MM_STORE(g+12*SIMDD, r0 * MM_LOAD(g+4 *SIMDD ) + i2 * b10 * s0);
        MM_STORE(g+13*SIMDD, r1 * MM_LOAD(g+5 *SIMDD ) + i2 * b11 * s1);
        MM_STORE(g+28*SIMDD, r2 * MM_LOAD(g+20*SIMDD ) + i2 * b10 * s2);
        MM_STORE(g+29*SIMDD, r3 * MM_LOAD(g+21*SIMDD ) + i2 * b11 * s3);
        MM_STORE(g+44*SIMDD, r4 * MM_LOAD(g+36*SIMDD ) + i2 * b10 * MM_LOAD(g+34*SIMDD ));
        MM_STORE(g+45*SIMDD, r5 * MM_LOAD(g+37*SIMDD ) + i2 * b11 * MM_LOAD(g+35*SIMDD ));
}

static inline void _g0_2d4d_0030(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *cx = bc->c0px;
        double *cy = bc->c0py;
        double *cz = bc->c0pz;
        double *b  = bc->b01;
        __MD b0 = MM_LOAD(b+0*SIMDD);
        __MD b1 = MM_LOAD(b+1*SIMDD);
        __MD i3 = MM_SET1(3.);
        __MD r0 = MM_LOAD(cx+0*SIMDD);
        __MD r1 = MM_LOAD(cx+1*SIMDD);
        __MD r2 = MM_LOAD(cy+0*SIMDD);
        __MD r3 = MM_LOAD(cy+1*SIMDD);
        __MD r4 = MM_LOAD(cz+0*SIMDD);
        __MD r5 = MM_LOAD(cz+1*SIMDD);
        __MD g16 = MM_LOAD(g+16*SIMDD);
        __MD g17 = MM_LOAD(g+17*SIMDD);
        MM_STORE(g+2 *SIMDD, r0);
        MM_STORE(g+3 *SIMDD, r1);
        MM_STORE(g+10*SIMDD, r2);
        MM_STORE(g+11*SIMDD, r3);
        MM_STORE(g+18*SIMDD, r4 * g16);
        MM_STORE(g+19*SIMDD, r5 * g17);
        MM_STORE(g+4 *SIMDD, r0 * r0 + b0);
        MM_STORE(g+5 *SIMDD, r1 * r1 + b1);
        MM_STORE(g+12*SIMDD, r2 * r2 + b0);
        MM_STORE(g+13*SIMDD, r3 * r3 + b1);
        MM_STORE(g+20*SIMDD,(r4 * r4 + b0)* g16);
        MM_STORE(g+21*SIMDD,(r5 * r5 + b1)* g17);
        MM_STORE(g+6 *SIMDD, r0 *(r0 * r0 + i3 * b0));
        MM_STORE(g+7 *SIMDD, r1 *(r1 * r1 + i3 * b1));
        MM_STORE(g+14*SIMDD, r2 *(r2 * r2 + i3 * b0));
        MM_STORE(g+15*SIMDD, r3 *(r3 * r3 + i3 * b1));
        MM_STORE(g+22*SIMDD,(r4 * r4 + i3 * b0) * MM_LOAD(g+18*SIMDD));
        MM_STORE(g+23*SIMDD,(r5 * r5 + i3 * b1) * MM_LOAD(g+19*SIMDD));
}

static inline void _g0_2d4d_0100(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *cx = bc->c00x;
        double *cy = bc->c00y;
        double *cz = bc->c00z;
        MM_STORE(g+1*SIMDD, MM_LOAD(cx+0*SIMDD));
        MM_STORE(g+3*SIMDD, MM_LOAD(cy+0*SIMDD));
        MM_STORE(g+5*SIMDD, MM_MUL(MM_LOAD(cz+0*SIMDD), MM_LOAD(g+4*SIMDD)));
}

static inline void _g0_2d4d_0101(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *c0x = bc->c00x;
        double *c0y = bc->c00y;
        double *c0z = bc->c00z;
        double *cpx = bc->c0px;
        double *cpy = bc->c0py;
        double *cpz = bc->c0pz;
        double *b   = bc->b00;
        __MD cx0 = MM_LOAD(c0x+0*SIMDD);
        __MD cx1 = MM_LOAD(c0x+1*SIMDD);
        __MD cy0 = MM_LOAD(c0y+0*SIMDD);
        __MD cy1 = MM_LOAD(c0y+1*SIMDD);
        __MD cz0 = MM_LOAD(c0z+0*SIMDD);
        __MD cz1 = MM_LOAD(c0z+1*SIMDD);
        __MD px0 = MM_LOAD(cpx+0*SIMDD);
        __MD px1 = MM_LOAD(cpx+1*SIMDD);
        __MD py0 = MM_LOAD(cpy+0*SIMDD);
        __MD py1 = MM_LOAD(cpy+1*SIMDD);
        __MD pz0 = MM_LOAD(cpz+0*SIMDD);
        __MD pz1 = MM_LOAD(cpz+1*SIMDD);
        __MD b0 = MM_LOAD(b+0*SIMDD);
        __MD b1 = MM_LOAD(b+1*SIMDD);
        __MD g16 = MM_LOAD(g+16*SIMDD);
        __MD g17 = MM_LOAD(g+17*SIMDD);
        MM_STORE(g+2 *SIMDD, px0);
        MM_STORE(g+3 *SIMDD, px1);
        MM_STORE(g+4 *SIMDD, cx0);
        MM_STORE(g+5 *SIMDD, cx1);
        MM_STORE(g+10*SIMDD, py0);
        MM_STORE(g+11*SIMDD, py1);
        MM_STORE(g+12*SIMDD, cy0);
        MM_STORE(g+13*SIMDD, cy1);
        MM_STORE(g+18*SIMDD, MM_MUL(pz0, g16));
        MM_STORE(g+19*SIMDD, MM_MUL(pz1, g17));
        MM_STORE(g+20*SIMDD, MM_MUL(cz0, g16));
        MM_STORE(g+21*SIMDD, MM_MUL(cz1, g17));
        MM_STORE(g+6 *SIMDD, MM_FMA(px0, cx0, b0));
        MM_STORE(g+7 *SIMDD, MM_FMA(px1, cx1, b1));
        MM_STORE(g+14*SIMDD, MM_FMA(py0, cy0, b0));
        MM_STORE(g+15*SIMDD, MM_FMA(py1, cy1, b1));
        MM_STORE(g+22*SIMDD, MM_MUL(MM_FMA(pz0, cz0, b0), g16));
        MM_STORE(g+23*SIMDD, MM_MUL(MM_FMA(pz1, cz1, b1), g17));
}

static inline void _g0_2d4d_0102(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *c0x = bc->c00x;
        double *c0y = bc->c00y;
        double *c0z = bc->c00z;
        double *cpx = bc->c0px;
        double *cpy = bc->c0py;
        double *cpz = bc->c0pz;
        double *b0 = bc->b00;
        double *b1 = bc->b01;
        __MD r0 = MM_LOAD(c0x+0*SIMDD);
        __MD r1 = MM_LOAD(c0x+1*SIMDD);
        __MD r2 = MM_LOAD(c0y+0*SIMDD);
        __MD r3 = MM_LOAD(c0y+1*SIMDD);
        __MD r4 = MM_LOAD(c0z+0*SIMDD);
        __MD r5 = MM_LOAD(c0z+1*SIMDD);
        __MD s0 = MM_LOAD(cpx+0*SIMDD);
        __MD s1 = MM_LOAD(cpx+1*SIMDD);
        __MD s2 = MM_LOAD(cpy+0*SIMDD);
        __MD s3 = MM_LOAD(cpy+1*SIMDD);
        __MD s4 = MM_LOAD(cpz+0*SIMDD);
        __MD s5 = MM_LOAD(cpz+1*SIMDD);
        __MD b00 = MM_LOAD(b0+0*SIMDD);
        __MD b01 = MM_LOAD(b0+1*SIMDD);
        __MD b10 = MM_LOAD(b1+0*SIMDD);
        __MD b11 = MM_LOAD(b1+1*SIMDD);
        __MD g24 = MM_LOAD(g+24*SIMDD);
        __MD g25 = MM_LOAD(g+25*SIMDD);
        MM_STORE(g+2 *SIMDD, s0);
        MM_STORE(g+3 *SIMDD, s1);
        MM_STORE(g+6 *SIMDD, r0);
        MM_STORE(g+7 *SIMDD, r1);
        MM_STORE(g+14*SIMDD, s2);
        MM_STORE(g+15*SIMDD, s3);
        MM_STORE(g+18*SIMDD, r2);
        MM_STORE(g+19*SIMDD, r3);
        MM_STORE(g+26*SIMDD, MM_MUL(s4, g24));
        MM_STORE(g+27*SIMDD, MM_MUL(s5, g25));
        MM_STORE(g+30*SIMDD, MM_MUL(r4, g24));
        MM_STORE(g+31*SIMDD, MM_MUL(r5, g25));
        MM_STORE(g+4 *SIMDD, MM_FMA(s0, s0, b10));
        MM_STORE(g+5 *SIMDD, MM_FMA(s1, s1, b11));
        MM_STORE(g+8 *SIMDD, MM_FMA(r0, s0, b00));
        MM_STORE(g+9 *SIMDD, MM_FMA(r1, s1, b01));
        MM_STORE(g+16*SIMDD, MM_FMA(s2, s2, b10));
        MM_STORE(g+17*SIMDD, MM_FMA(s3, s3, b11));
        MM_STORE(g+20*SIMDD, MM_FMA(r2, s2, b00));
        MM_STORE(g+21*SIMDD, MM_FMA(r3, s3, b01));
        MM_STORE(g+28*SIMDD, MM_MUL(MM_FMA(s4, s4, b10), g24));
        MM_STORE(g+29*SIMDD, MM_MUL(MM_FMA(s5, s5, b11), g25));
        MM_STORE(g+32*SIMDD, MM_MUL(MM_FMA(r4, s4, b00), g24));
        MM_STORE(g+33*SIMDD, MM_MUL(MM_FMA(r5, s5, b01), g25));
        MM_STORE(g+10*SIMDD, s0 * (MM_LOAD(g+8 *SIMDD) + b00) + b10 * r0);
        MM_STORE(g+11*SIMDD, s1 * (MM_LOAD(g+9 *SIMDD) + b01) + b11 * r1);
        MM_STORE(g+22*SIMDD, s2 * (MM_LOAD(g+20*SIMDD) + b00) + b10 * r2);
        MM_STORE(g+23*SIMDD, s3 * (MM_LOAD(g+21*SIMDD) + b01) + b11 * r3);
        MM_STORE(g+34*SIMDD, s4 * MM_LOAD(g+32*SIMDD) + b00 * MM_LOAD(g+26*SIMDD) + b10 * MM_LOAD(g+30*SIMDD));
        MM_STORE(g+35*SIMDD, s5 * MM_LOAD(g+33*SIMDD) + b01 * MM_LOAD(g+27*SIMDD) + b11 * MM_LOAD(g+31*SIMDD));
}

static inline void _g0_2d4d_0110(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *c0x = bc->c00x;
        double *c0y = bc->c00y;
        double *c0z = bc->c00z;
        double *cpx = bc->c0px;
        double *cpy = bc->c0py;
        double *cpz = bc->c0pz;
        double *b   = bc->b00;
        __MD b0  = MM_LOAD(b+0*SIMDD);
        __MD b1  = MM_LOAD(b+1*SIMDD);
        __MD r0 = MM_LOAD(c0x+0*SIMDD);
        __MD r1 = MM_LOAD(c0x+1*SIMDD);
        __MD r2 = MM_LOAD(c0y+0*SIMDD);
        __MD r3 = MM_LOAD(c0y+1*SIMDD);
        __MD r4 = MM_LOAD(c0z+0*SIMDD);
        __MD r5 = MM_LOAD(c0z+1*SIMDD);
        __MD s0 = MM_LOAD(cpx+0*SIMDD);
        __MD s1 = MM_LOAD(cpx+1*SIMDD);
        __MD s2 = MM_LOAD(cpy+0*SIMDD);
        __MD s3 = MM_LOAD(cpy+1*SIMDD);
        __MD s4 = MM_LOAD(cpz+0*SIMDD);
        __MD s5 = MM_LOAD(cpz+1*SIMDD);
        __MD g16 = MM_LOAD(g+16*SIMDD);
        __MD g17 = MM_LOAD(g+17*SIMDD);
        MM_STORE(g+2 *SIMDD, s0);
        MM_STORE(g+3 *SIMDD, s1);
        MM_STORE(g+10*SIMDD, s2);
        MM_STORE(g+11*SIMDD, s3);
        MM_STORE(g+18*SIMDD, s4 * g16);
        MM_STORE(g+19*SIMDD, s5 * g17);
        MM_STORE(g+4 *SIMDD, r0);
        MM_STORE(g+5 *SIMDD, r1);
        MM_STORE(g+12*SIMDD, r2);
        MM_STORE(g+13*SIMDD, r3);
        MM_STORE(g+20*SIMDD, r4 * g16);
        MM_STORE(g+21*SIMDD, r5 * g17);
        MM_STORE(g+6 *SIMDD, s0 * r0 + b0);
        MM_STORE(g+7 *SIMDD, s1 * r1 + b1);
        MM_STORE(g+14*SIMDD, s2 * r2 + b0);
        MM_STORE(g+15*SIMDD, s3 * r3 + b1);
        MM_STORE(g+22*SIMDD,(s4 * r4 + b0) * g16);
        MM_STORE(g+23*SIMDD,(s5 * r5 + b1) * g17);
}

static inline void _g0_2d4d_0111(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *c0x = bc->c00x;
        double *c0y = bc->c00y;
        double *c0z = bc->c00z;
        double *cpx = bc->c0px;
        double *cpy = bc->c0py;
        double *cpz = bc->c0pz;
        double *b0  = bc->b00;
        double *b1  = bc->b01;
        double *rkl = envs->rkrl;
        ALIGNMM double rcp[6*SIMDD];
        _make_rc(rcp, cpx, cpy, cpz, rkl);
        __MD b00 = MM_LOAD(b0+0*SIMDD);
        __MD b01 = MM_LOAD(b0+1*SIMDD);
        __MD b10 = MM_LOAD(b1+0*SIMDD);
        __MD b11 = MM_LOAD(b1+1*SIMDD);
        __MD g48 = MM_LOAD(g+48*SIMDD);
        __MD g49 = MM_LOAD(g+49*SIMDD);
        MM_STORE(g+2 *SIMDD, MM_LOAD(rcp+0*SIMDD));
        MM_STORE(g+3 *SIMDD, MM_LOAD(rcp+1*SIMDD));
        MM_STORE(g+4 *SIMDD, MM_LOAD(cpx+0*SIMDD));
        MM_STORE(g+5 *SIMDD, MM_LOAD(cpx+1*SIMDD));
        MM_STORE(g+12*SIMDD, MM_LOAD(c0x+0*SIMDD));
        MM_STORE(g+13*SIMDD, MM_LOAD(c0x+1*SIMDD));
        MM_STORE(g+26*SIMDD, MM_LOAD(rcp+2*SIMDD));
        MM_STORE(g+27*SIMDD, MM_LOAD(rcp+3*SIMDD));
        MM_STORE(g+28*SIMDD, MM_LOAD(cpy+0*SIMDD));
        MM_STORE(g+29*SIMDD, MM_LOAD(cpy+1*SIMDD));
        MM_STORE(g+36*SIMDD, MM_LOAD(c0y+0*SIMDD));
        MM_STORE(g+37*SIMDD, MM_LOAD(c0y+1*SIMDD));
        MM_STORE(g+50*SIMDD, MM_LOAD(rcp+4*SIMDD) * g48);
        MM_STORE(g+51*SIMDD, MM_LOAD(rcp+5*SIMDD) * g49);
        MM_STORE(g+52*SIMDD, MM_LOAD(cpz+0*SIMDD) * g48);
        MM_STORE(g+53*SIMDD, MM_LOAD(cpz+1*SIMDD) * g49);
        MM_STORE(g+60*SIMDD, MM_LOAD(c0z+0*SIMDD) * g48);
        MM_STORE(g+61*SIMDD, MM_LOAD(c0z+1*SIMDD) * g49);
        MM_STORE(g+14*SIMDD, MM_LOAD(c0x+0*SIMDD) * MM_LOAD(rcp+0*SIMDD) + b00);
        MM_STORE(g+15*SIMDD, MM_LOAD(c0x+1*SIMDD) * MM_LOAD(rcp+1*SIMDD) + b01);
        MM_STORE(g+16*SIMDD, MM_LOAD(c0x+0*SIMDD) * MM_LOAD(cpx+0*SIMDD) + b00);
        MM_STORE(g+17*SIMDD, MM_LOAD(c0x+1*SIMDD) * MM_LOAD(cpx+1*SIMDD) + b01);
        MM_STORE(g+6 *SIMDD, MM_LOAD(cpx+0*SIMDD) * MM_LOAD(rcp+0*SIMDD) + b10);
        MM_STORE(g+7 *SIMDD, MM_LOAD(cpx+1*SIMDD) * MM_LOAD(rcp+1*SIMDD) + b11);
        MM_STORE(g+30*SIMDD, MM_LOAD(cpy+0*SIMDD) * MM_LOAD(rcp+2*SIMDD) + b10);
        MM_STORE(g+31*SIMDD, MM_LOAD(cpy+1*SIMDD) * MM_LOAD(rcp+3*SIMDD) + b11);
        MM_STORE(g+38*SIMDD, MM_LOAD(c0y+0*SIMDD) * MM_LOAD(rcp+2*SIMDD) + b00);
        MM_STORE(g+39*SIMDD, MM_LOAD(c0y+1*SIMDD) * MM_LOAD(rcp+3*SIMDD) + b01);
        MM_STORE(g+40*SIMDD, MM_LOAD(c0y+0*SIMDD) * MM_LOAD(cpy+0*SIMDD) + b00);
        MM_STORE(g+41*SIMDD, MM_LOAD(c0y+1*SIMDD) * MM_LOAD(cpy+1*SIMDD) + b01);
        MM_STORE(g+54*SIMDD,(MM_LOAD(cpz+0*SIMDD) * MM_LOAD(rcp+4*SIMDD) + b10) * g48);
        MM_STORE(g+55*SIMDD,(MM_LOAD(cpz+1*SIMDD) * MM_LOAD(rcp+5*SIMDD) + b11) * g49);
        MM_STORE(g+62*SIMDD,(MM_LOAD(c0z+0*SIMDD) * MM_LOAD(rcp+4*SIMDD) + b00) * g48);
        MM_STORE(g+63*SIMDD,(MM_LOAD(c0z+1*SIMDD) * MM_LOAD(rcp+5*SIMDD) + b01) * g49);
        MM_STORE(g+64*SIMDD,(MM_LOAD(c0z+0*SIMDD) * MM_LOAD(cpz+0*SIMDD) + b00) * g48);
        MM_STORE(g+65*SIMDD,(MM_LOAD(c0z+1*SIMDD) * MM_LOAD(cpz+1*SIMDD) + b01) * g49);
        MM_STORE(g+18*SIMDD, MM_LOAD(c0x+0*SIMDD) * MM_LOAD(g+6 *SIMDD) + b00 * (MM_LOAD(rcp+0*SIMDD) + MM_LOAD(cpx+0*SIMDD)));
        MM_STORE(g+19*SIMDD, MM_LOAD(c0x+1*SIMDD) * MM_LOAD(g+7 *SIMDD) + b01 * (MM_LOAD(rcp+1*SIMDD) + MM_LOAD(cpx+1*SIMDD)));
        MM_STORE(g+42*SIMDD, MM_LOAD(c0y+0*SIMDD) * MM_LOAD(g+30*SIMDD) + b00 * (MM_LOAD(rcp+2*SIMDD) + MM_LOAD(cpy+0*SIMDD)));
        MM_STORE(g+43*SIMDD, MM_LOAD(c0y+1*SIMDD) * MM_LOAD(g+31*SIMDD) + b01 * (MM_LOAD(rcp+3*SIMDD) + MM_LOAD(cpy+1*SIMDD)));
        MM_STORE(g+66*SIMDD, MM_LOAD(c0z+0*SIMDD) * MM_LOAD(g+54*SIMDD) + b00 * (MM_LOAD(g+50*SIMDD) + MM_LOAD(g+52*SIMDD)));
        MM_STORE(g+67*SIMDD, MM_LOAD(c0z+1*SIMDD) * MM_LOAD(g+55*SIMDD) + b01 * (MM_LOAD(g+51*SIMDD) + MM_LOAD(g+53*SIMDD)));
}

static inline void _g0_2d4d_0120(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *c0x = bc->c00x;
        double *c0y = bc->c00y;
        double *c0z = bc->c00z;
        double *cpx = bc->c0px;
        double *cpy = bc->c0py;
        double *cpz = bc->c0pz;
        double *b0 = bc->b00;
        double *b1 = bc->b01;
        __MD r0 = MM_LOAD(c0x+0*SIMDD);
        __MD r1 = MM_LOAD(c0x+1*SIMDD);
        __MD r2 = MM_LOAD(c0y+0*SIMDD);
        __MD r3 = MM_LOAD(c0y+1*SIMDD);
        __MD r4 = MM_LOAD(c0z+0*SIMDD);
        __MD r5 = MM_LOAD(c0z+1*SIMDD);
        __MD s0 = MM_LOAD(cpx+0*SIMDD);
        __MD s1 = MM_LOAD(cpx+1*SIMDD);
        __MD s2 = MM_LOAD(cpy+0*SIMDD);
        __MD s3 = MM_LOAD(cpy+1*SIMDD);
        __MD s4 = MM_LOAD(cpz+0*SIMDD);
        __MD s5 = MM_LOAD(cpz+1*SIMDD);
        __MD b00 = MM_LOAD(b0+0*SIMDD);
        __MD b01 = MM_LOAD(b0+1*SIMDD);
        __MD b10 = MM_LOAD(b1+0*SIMDD);
        __MD b11 = MM_LOAD(b1+1*SIMDD);
        __MD g24 = MM_LOAD(g+24*SIMDD);
        __MD g25 = MM_LOAD(g+25*SIMDD);
        MM_STORE(g+2 *SIMDD, s0);
        MM_STORE(g+3 *SIMDD, s1);
        MM_STORE(g+6 *SIMDD, r0);
        MM_STORE(g+7 *SIMDD, r1);
        MM_STORE(g+14*SIMDD, s2);
        MM_STORE(g+15*SIMDD, s3);
        MM_STORE(g+18*SIMDD, r2);
        MM_STORE(g+19*SIMDD, r3);
        MM_STORE(g+26*SIMDD, MM_MUL(s4, g24));
        MM_STORE(g+27*SIMDD, MM_MUL(s5, g25));
        MM_STORE(g+30*SIMDD, MM_MUL(r4, g24));
        MM_STORE(g+31*SIMDD, MM_MUL(r5, g25));
        MM_STORE(g+4 *SIMDD, MM_FMA(s0, s0, b10));
        MM_STORE(g+5 *SIMDD, MM_FMA(s1, s1, b11));
        MM_STORE(g+8 *SIMDD, MM_FMA(r0, s0, b00));
        MM_STORE(g+9 *SIMDD, MM_FMA(r1, s1, b01));
        MM_STORE(g+16*SIMDD, MM_FMA(s2, s2, b10));
        MM_STORE(g+17*SIMDD, MM_FMA(s3, s3, b11));
        MM_STORE(g+20*SIMDD, MM_FMA(r2, s2, b00));
        MM_STORE(g+21*SIMDD, MM_FMA(r3, s3, b01));
        MM_STORE(g+28*SIMDD, MM_MUL(MM_FMA(s4, s4, b10), g24));
        MM_STORE(g+29*SIMDD, MM_MUL(MM_FMA(s5, s5, b11), g25));
        MM_STORE(g+32*SIMDD, MM_MUL(MM_FMA(r4, s4, b00), g24));
        MM_STORE(g+33*SIMDD, MM_MUL(MM_FMA(r5, s5, b01), g25));
        MM_STORE(g+10*SIMDD, s0 * (MM_LOAD(g+8 *SIMDD) + b00) + b10 * r0);
        MM_STORE(g+11*SIMDD, s1 * (MM_LOAD(g+9 *SIMDD) + b01) + b11 * r1);
        MM_STORE(g+22*SIMDD, s2 * (MM_LOAD(g+20*SIMDD) + b00) + b10 * r2);
        MM_STORE(g+23*SIMDD, s3 * (MM_LOAD(g+21*SIMDD) + b01) + b11 * r3);
        MM_STORE(g+34*SIMDD, s4 * MM_LOAD(g+32*SIMDD) + b00 * MM_LOAD(g+26*SIMDD) + b10 * MM_LOAD(g+30*SIMDD));
        MM_STORE(g+35*SIMDD, s5 * MM_LOAD(g+33*SIMDD) + b01 * MM_LOAD(g+27*SIMDD) + b11 * MM_LOAD(g+31*SIMDD));
}

static inline void _g0_2d4d_0200(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *cx = bc->c00x;
        double *cy = bc->c00y;
        double *cz = bc->c00z;
        double *b =  bc->b10;
        __MD cx0 = MM_LOAD(cx+0*SIMDD);
        __MD cx1 = MM_LOAD(cx+1*SIMDD);
        __MD cy0 = MM_LOAD(cy+0*SIMDD);
        __MD cy1 = MM_LOAD(cy+1*SIMDD);
        __MD cz0 = MM_LOAD(cz+0*SIMDD);
        __MD cz1 = MM_LOAD(cz+1*SIMDD);
        __MD b0 = MM_LOAD(b+0*SIMDD);
        __MD b1 = MM_LOAD(b+1*SIMDD);
        __MD g12 = MM_LOAD(g+12*SIMDD);
        __MD g13 = MM_LOAD(g+13*SIMDD);
        MM_STORE(g+2 *SIMDD, cx0);
        MM_STORE(g+3 *SIMDD, cx1);
        MM_STORE(g+8 *SIMDD, cy0);
        MM_STORE(g+9 *SIMDD, cy1);
        MM_STORE(g+14*SIMDD, MM_MUL(cz0, g12));
        MM_STORE(g+15*SIMDD, MM_MUL(cz1, g13));
        MM_STORE(g+4 *SIMDD, MM_FMA(cx0, cx0, b0));
        MM_STORE(g+5 *SIMDD, MM_FMA(cx1, cx1, b1));
        MM_STORE(g+10*SIMDD, MM_FMA(cy0, cy0, b0));
        MM_STORE(g+11*SIMDD, MM_FMA(cy1, cy1, b1));
        MM_STORE(g+16*SIMDD, MM_MUL(MM_FMA(cz0, cz0, b0), g12));
        MM_STORE(g+17*SIMDD, MM_MUL(MM_FMA(cz1, cz1, b1), g13));
}

static inline void _g0_2d4d_0201(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *c0x = bc->c00x;
        double *c0y = bc->c00y;
        double *c0z = bc->c00z;
        double *cpx = bc->c0px;
        double *cpy = bc->c0py;
        double *cpz = bc->c0pz;
        double *b0  = bc->b00;
        double *b1  = bc->b10;
        __MD b00 = MM_LOAD(b0+0*SIMDD);
        __MD b01 = MM_LOAD(b0+1*SIMDD);
        __MD b10 = MM_LOAD(b1+0*SIMDD);
        __MD b11 = MM_LOAD(b1+1*SIMDD);
        __MD g24 = MM_LOAD(g+24*SIMDD);
        __MD g25 = MM_LOAD(g+25*SIMDD);
        __MD i2 = MM_SET1(2.);
        MM_STORE(g+2 *SIMDD, MM_LOAD(cpx+0*SIMDD));
        MM_STORE(g+3 *SIMDD, MM_LOAD(cpx+1*SIMDD));
        MM_STORE(g+4 *SIMDD, MM_LOAD(c0x+0*SIMDD));
        MM_STORE(g+5 *SIMDD, MM_LOAD(c0x+1*SIMDD));
        MM_STORE(g+14*SIMDD, MM_LOAD(cpy+0*SIMDD));
        MM_STORE(g+15*SIMDD, MM_LOAD(cpy+1*SIMDD));
        MM_STORE(g+16*SIMDD, MM_LOAD(c0y+0*SIMDD));
        MM_STORE(g+17*SIMDD, MM_LOAD(c0y+1*SIMDD));
        MM_STORE(g+26*SIMDD, MM_LOAD(cpz+0*SIMDD) * g24);
        MM_STORE(g+27*SIMDD, MM_LOAD(cpz+1*SIMDD) * g25);
        MM_STORE(g+28*SIMDD, MM_LOAD(c0z+0*SIMDD) * g24);
        MM_STORE(g+29*SIMDD, MM_LOAD(c0z+1*SIMDD) * g25);
        MM_STORE(g+6 *SIMDD, MM_LOAD(cpx+0*SIMDD) * MM_LOAD(c0x+0*SIMDD) + b00);
        MM_STORE(g+7 *SIMDD, MM_LOAD(cpx+1*SIMDD) * MM_LOAD(c0x+1*SIMDD) + b01);
        MM_STORE(g+8 *SIMDD, MM_LOAD(c0x+0*SIMDD) * MM_LOAD(c0x+0*SIMDD) + b10);
        MM_STORE(g+9 *SIMDD, MM_LOAD(c0x+1*SIMDD) * MM_LOAD(c0x+1*SIMDD) + b11);
        MM_STORE(g+18*SIMDD, MM_LOAD(cpy+0*SIMDD) * MM_LOAD(c0y+0*SIMDD) + b00);
        MM_STORE(g+19*SIMDD, MM_LOAD(cpy+1*SIMDD) * MM_LOAD(c0y+1*SIMDD) + b01);
        MM_STORE(g+20*SIMDD, MM_LOAD(c0y+0*SIMDD) * MM_LOAD(c0y+0*SIMDD) + b10);
        MM_STORE(g+21*SIMDD, MM_LOAD(c0y+1*SIMDD) * MM_LOAD(c0y+1*SIMDD) + b11);
        MM_STORE(g+30*SIMDD,(MM_LOAD(cpz+0*SIMDD) * MM_LOAD(c0z+0*SIMDD) + b00) * g24);
        MM_STORE(g+31*SIMDD,(MM_LOAD(cpz+1*SIMDD) * MM_LOAD(c0z+1*SIMDD) + b01) * g25);
        MM_STORE(g+32*SIMDD,(MM_LOAD(c0z+0*SIMDD) * MM_LOAD(c0z+0*SIMDD) + b10) * g24);
        MM_STORE(g+33*SIMDD,(MM_LOAD(c0z+1*SIMDD) * MM_LOAD(c0z+1*SIMDD) + b11) * g25);
        MM_STORE(g+10*SIMDD, MM_LOAD(cpx+0*SIMDD) * MM_LOAD(g+8 *SIMDD ) + i2 * b00 * MM_LOAD(c0x+0*SIMDD));
        MM_STORE(g+11*SIMDD, MM_LOAD(cpx+1*SIMDD) * MM_LOAD(g+9 *SIMDD ) + i2 * b01 * MM_LOAD(c0x+1*SIMDD));
        MM_STORE(g+22*SIMDD, MM_LOAD(cpy+0*SIMDD) * MM_LOAD(g+20*SIMDD ) + i2 * b00 * MM_LOAD(c0y+0*SIMDD));
        MM_STORE(g+23*SIMDD, MM_LOAD(cpy+1*SIMDD) * MM_LOAD(g+21*SIMDD ) + i2 * b01 * MM_LOAD(c0y+1*SIMDD));
        MM_STORE(g+34*SIMDD, MM_LOAD(cpz+0*SIMDD) * MM_LOAD(g+32*SIMDD ) + i2 * b00 * MM_LOAD(g+28*SIMDD));
        MM_STORE(g+35*SIMDD, MM_LOAD(cpz+1*SIMDD) * MM_LOAD(g+33*SIMDD ) + i2 * b01 * MM_LOAD(g+29*SIMDD));
}

static inline void _g0_2d4d_0210(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *c0x = bc->c00x;
        double *c0y = bc->c00y;
        double *c0z = bc->c00z;
        double *cpx = bc->c0px;
        double *cpy = bc->c0py;
        double *cpz = bc->c0pz;
        double *b0 = bc->b00;
        double *b1 = bc->b01;
        __MD r0 = MM_LOAD(c0x+0*SIMDD);
        __MD r1 = MM_LOAD(c0x+1*SIMDD);
        __MD r2 = MM_LOAD(c0y+0*SIMDD);
        __MD r3 = MM_LOAD(c0y+1*SIMDD);
        __MD r4 = MM_LOAD(c0z+0*SIMDD);
        __MD r5 = MM_LOAD(c0z+1*SIMDD);
        __MD s0 = MM_LOAD(cpx+0*SIMDD);
        __MD s1 = MM_LOAD(cpx+1*SIMDD);
        __MD s2 = MM_LOAD(cpy+0*SIMDD);
        __MD s3 = MM_LOAD(cpy+1*SIMDD);
        __MD s4 = MM_LOAD(cpz+0*SIMDD);
        __MD s5 = MM_LOAD(cpz+1*SIMDD);
        __MD b00 = MM_LOAD(b0+0*SIMDD);
        __MD b01 = MM_LOAD(b0+1*SIMDD);
        __MD b10 = MM_LOAD(b1+0*SIMDD);
        __MD b11 = MM_LOAD(b1+1*SIMDD);
        __MD g24 = MM_LOAD(g+24*SIMDD);
        __MD g25 = MM_LOAD(g+25*SIMDD);
        MM_STORE(g+2 *SIMDD, s0);
        MM_STORE(g+3 *SIMDD, s1);
        MM_STORE(g+4 *SIMDD, r0);
        MM_STORE(g+5 *SIMDD, r1);
        MM_STORE(g+14*SIMDD, s2);
        MM_STORE(g+15*SIMDD, s3);
        MM_STORE(g+16*SIMDD, r2);
        MM_STORE(g+17*SIMDD, r3);
        MM_STORE(g+26*SIMDD, MM_MUL(s4, g24));
        MM_STORE(g+27*SIMDD, MM_MUL(s5, g25));
        MM_STORE(g+28*SIMDD, MM_MUL(r4, g24));
        MM_STORE(g+29*SIMDD, MM_MUL(r5, g25));
        MM_STORE(g+6 *SIMDD, MM_FMA(s0, r0, b00));
        MM_STORE(g+7 *SIMDD, MM_FMA(s1, r1, b01));
        MM_STORE(g+8 *SIMDD, MM_FMA(r0, r0, b10));
        MM_STORE(g+9 *SIMDD, MM_FMA(r1, r1, b11));
        MM_STORE(g+18*SIMDD, MM_FMA(s2, r2, b00));
        MM_STORE(g+19*SIMDD, MM_FMA(s3, r3, b01));
        MM_STORE(g+20*SIMDD, MM_FMA(r2, r2, b10));
        MM_STORE(g+21*SIMDD, MM_FMA(r3, r3, b11));
        MM_STORE(g+30*SIMDD, MM_MUL(MM_FMA(s4, r4, b00), g24));
        MM_STORE(g+31*SIMDD, MM_MUL(MM_FMA(s5, r5, b01), g25));
        MM_STORE(g+32*SIMDD, MM_MUL(MM_FMA(r4, r4, b10), g24));
        MM_STORE(g+33*SIMDD, MM_MUL(MM_FMA(r5, r5, b11), g25));
        MM_STORE(g+10*SIMDD, r0 * (MM_LOAD(g+6 *SIMDD) + b00) + b10 * s0                 );
        MM_STORE(g+11*SIMDD, r1 * (MM_LOAD(g+7 *SIMDD) + b01) + b11 * s1                 );
        MM_STORE(g+22*SIMDD, r2 * (MM_LOAD(g+18*SIMDD) + b00) + b10 * s2                 );
        MM_STORE(g+23*SIMDD, r3 * (MM_LOAD(g+19*SIMDD) + b01) + b11 * s3                 );
        MM_STORE(g+34*SIMDD, r4 * MM_LOAD(g+30*SIMDD) + b10 * MM_LOAD(g+26*SIMDD) + b00 * MM_LOAD(g+28*SIMDD));
        MM_STORE(g+35*SIMDD, r5 * MM_LOAD(g+31*SIMDD) + b11 * MM_LOAD(g+27*SIMDD) + b01 * MM_LOAD(g+29*SIMDD));
}

static inline void _g0_2d4d_0300(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *cx = bc->c00x;
        double *cy = bc->c00y;
        double *cz = bc->c00z;
        double *b =  bc->b10;
        __MD b0 = MM_LOAD(b+0*SIMDD);
        __MD b1 = MM_LOAD(b+1*SIMDD);
        __MD g16 = MM_LOAD(g+16*SIMDD);
        __MD g17 = MM_LOAD(g+17*SIMDD);
        __MD i3 = MM_SET1(3.);
        MM_STORE(g+2 *SIMDD, MM_LOAD(cx+0*SIMDD));
        MM_STORE(g+3 *SIMDD, MM_LOAD(cx+1*SIMDD));
        MM_STORE(g+10*SIMDD, MM_LOAD(cy+0*SIMDD));
        MM_STORE(g+11*SIMDD, MM_LOAD(cy+1*SIMDD));
        MM_STORE(g+18*SIMDD, MM_LOAD(cz+0*SIMDD) * g16);
        MM_STORE(g+19*SIMDD, MM_LOAD(cz+1*SIMDD) * g17);
        MM_STORE(g+4 *SIMDD, MM_LOAD(cx+0*SIMDD) * MM_LOAD(cx+0*SIMDD) + b0);
        MM_STORE(g+5 *SIMDD, MM_LOAD(cx+1*SIMDD) * MM_LOAD(cx+1*SIMDD) + b1);
        MM_STORE(g+12*SIMDD, MM_LOAD(cy+0*SIMDD) * MM_LOAD(cy+0*SIMDD) + b0);
        MM_STORE(g+13*SIMDD, MM_LOAD(cy+1*SIMDD) * MM_LOAD(cy+1*SIMDD) + b1);
        MM_STORE(g+20*SIMDD,(MM_LOAD(cz+0*SIMDD) * MM_LOAD(cz+0*SIMDD) + b0)* g16);
        MM_STORE(g+21*SIMDD,(MM_LOAD(cz+1*SIMDD) * MM_LOAD(cz+1*SIMDD) + b1)* g17);
        MM_STORE(g+6 *SIMDD, MM_LOAD(cx+0*SIMDD) *(MM_LOAD(cx+0*SIMDD) * MM_LOAD(cx+0*SIMDD) + i3 * b0));
        MM_STORE(g+7 *SIMDD, MM_LOAD(cx+1*SIMDD) *(MM_LOAD(cx+1*SIMDD) * MM_LOAD(cx+1*SIMDD) + i3 * b1));
        MM_STORE(g+14*SIMDD, MM_LOAD(cy+0*SIMDD) *(MM_LOAD(cy+0*SIMDD) * MM_LOAD(cy+0*SIMDD) + i3 * b0));
        MM_STORE(g+15*SIMDD, MM_LOAD(cy+1*SIMDD) *(MM_LOAD(cy+1*SIMDD) * MM_LOAD(cy+1*SIMDD) + i3 * b1));
        MM_STORE(g+22*SIMDD,(MM_LOAD(cz+0*SIMDD) * MM_LOAD(cz+0*SIMDD) + i3 * b0)* MM_LOAD(g+18*SIMDD));
        MM_STORE(g+23*SIMDD,(MM_LOAD(cz+1*SIMDD) * MM_LOAD(cz+1*SIMDD) + i3 * b1)* MM_LOAD(g+19*SIMDD));
}

static inline void _g0_2d4d_1000(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *cx = bc->c00x;
        double *cy = bc->c00y;
        double *cz = bc->c00z;
        MM_STORE(g+1*SIMDD, MM_LOAD(cx+0*SIMDD));
        MM_STORE(g+3*SIMDD, MM_LOAD(cy+0*SIMDD));
        MM_STORE(g+5*SIMDD, MM_MUL(MM_LOAD(cz+0*SIMDD), MM_LOAD(g+4*SIMDD)));
}

static inline void _g0_2d4d_1001(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *c0x = bc->c00x;
        double *c0y = bc->c00y;
        double *c0z = bc->c00z;
        double *cpx = bc->c0px;
        double *cpy = bc->c0py;
        double *cpz = bc->c0pz;
        double *b   = bc->b00;
        __MD cx0 = MM_LOAD(c0x+0*SIMDD);
        __MD cx1 = MM_LOAD(c0x+1*SIMDD);
        __MD cy0 = MM_LOAD(c0y+0*SIMDD);
        __MD cy1 = MM_LOAD(c0y+1*SIMDD);
        __MD cz0 = MM_LOAD(c0z+0*SIMDD);
        __MD cz1 = MM_LOAD(c0z+1*SIMDD);
        __MD px0 = MM_LOAD(cpx+0*SIMDD);
        __MD px1 = MM_LOAD(cpx+1*SIMDD);
        __MD py0 = MM_LOAD(cpy+0*SIMDD);
        __MD py1 = MM_LOAD(cpy+1*SIMDD);
        __MD pz0 = MM_LOAD(cpz+0*SIMDD);
        __MD pz1 = MM_LOAD(cpz+1*SIMDD);
        __MD b0 = MM_LOAD(b+0*SIMDD);
        __MD b1 = MM_LOAD(b+1*SIMDD);
        __MD g16 = MM_LOAD(g+16*SIMDD);
        __MD g17 = MM_LOAD(g+17*SIMDD);
        MM_STORE(g+2 *SIMDD, px0);
        MM_STORE(g+3 *SIMDD, px1);
        MM_STORE(g+4 *SIMDD, cx0);
        MM_STORE(g+5 *SIMDD, cx1);
        MM_STORE(g+10*SIMDD, py0);
        MM_STORE(g+11*SIMDD, py1);
        MM_STORE(g+12*SIMDD, cy0);
        MM_STORE(g+13*SIMDD, cy1);
        MM_STORE(g+18*SIMDD, MM_MUL(pz0, g16));
        MM_STORE(g+19*SIMDD, MM_MUL(pz1, g17));
        MM_STORE(g+20*SIMDD, MM_MUL(cz0, g16));
        MM_STORE(g+21*SIMDD, MM_MUL(cz1, g17));
        MM_STORE(g+6 *SIMDD, MM_FMA(px0, cx0, b0));
        MM_STORE(g+7 *SIMDD, MM_FMA(px1, cx1, b1));
        MM_STORE(g+14*SIMDD, MM_FMA(py0, cy0, b0));
        MM_STORE(g+15*SIMDD, MM_FMA(py1, cy1, b1));
        MM_STORE(g+22*SIMDD, MM_MUL(MM_FMA(pz0, cz0, b0), g16));
        MM_STORE(g+23*SIMDD, MM_MUL(MM_FMA(pz1, cz1, b1), g17));
}

static inline void _g0_2d4d_1002(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *c0x = bc->c00x;
        double *c0y = bc->c00y;
        double *c0z = bc->c00z;
        double *cpx = bc->c0px;
        double *cpy = bc->c0py;
        double *cpz = bc->c0pz;
        double *b0 = bc->b00;
        double *b1 = bc->b01;
        __MD r0 = MM_LOAD(c0x+0*SIMDD);
        __MD r1 = MM_LOAD(c0x+1*SIMDD);
        __MD r2 = MM_LOAD(c0y+0*SIMDD);
        __MD r3 = MM_LOAD(c0y+1*SIMDD);
        __MD r4 = MM_LOAD(c0z+0*SIMDD);
        __MD r5 = MM_LOAD(c0z+1*SIMDD);
        __MD s0 = MM_LOAD(cpx+0*SIMDD);
        __MD s1 = MM_LOAD(cpx+1*SIMDD);
        __MD s2 = MM_LOAD(cpy+0*SIMDD);
        __MD s3 = MM_LOAD(cpy+1*SIMDD);
        __MD s4 = MM_LOAD(cpz+0*SIMDD);
        __MD s5 = MM_LOAD(cpz+1*SIMDD);
        __MD b00 = MM_LOAD(b0+0*SIMDD);
        __MD b01 = MM_LOAD(b0+1*SIMDD);
        __MD b10 = MM_LOAD(b1+0*SIMDD);
        __MD b11 = MM_LOAD(b1+1*SIMDD);
        __MD g24 = MM_LOAD(g+24*SIMDD);
        __MD g25 = MM_LOAD(g+25*SIMDD);
        MM_STORE(g+2 *SIMDD, r0);
        MM_STORE(g+3 *SIMDD, r1);
        MM_STORE(g+4 *SIMDD, s0);
        MM_STORE(g+5 *SIMDD, s1);
        MM_STORE(g+14*SIMDD, r2);
        MM_STORE(g+15*SIMDD, r3);
        MM_STORE(g+16*SIMDD, s2);
        MM_STORE(g+17*SIMDD, s3);
        MM_STORE(g+26*SIMDD, MM_MUL(r4, g24));
        MM_STORE(g+27*SIMDD, MM_MUL(r5, g25));
        MM_STORE(g+28*SIMDD, MM_MUL(s4, g24));
        MM_STORE(g+29*SIMDD, MM_MUL(s5, g25));
        MM_STORE(g+6 *SIMDD, MM_FMA(r0, s0, b00));
        MM_STORE(g+7 *SIMDD, MM_FMA(r1, s1, b01));
        MM_STORE(g+8 *SIMDD, MM_FMA(s0, s0, b10));
        MM_STORE(g+9 *SIMDD, MM_FMA(s1, s1, b11));
        MM_STORE(g+18*SIMDD, MM_FMA(r2, s2, b00));
        MM_STORE(g+19*SIMDD, MM_FMA(r3, s3, b01));
        MM_STORE(g+20*SIMDD, MM_FMA(s2, s2, b10));
        MM_STORE(g+21*SIMDD, MM_FMA(s3, s3, b11));
        MM_STORE(g+30*SIMDD, MM_MUL(MM_FMA(r4, s4, b00), g24));
        MM_STORE(g+31*SIMDD, MM_MUL(MM_FMA(r5, s5, b01), g25));
        MM_STORE(g+32*SIMDD, MM_MUL(MM_FMA(s4, s4, b10), g24));
        MM_STORE(g+33*SIMDD, MM_MUL(MM_FMA(s5, s5, b11), g25));
        MM_STORE(g+10*SIMDD, s0 * (MM_LOAD(g+6 *SIMDD) + b00) + b10 * r0);
        MM_STORE(g+11*SIMDD, s1 * (MM_LOAD(g+7 *SIMDD) + b01) + b11 * r1);
        MM_STORE(g+22*SIMDD, s2 * (MM_LOAD(g+18*SIMDD) + b00) + b10 * r2);
        MM_STORE(g+23*SIMDD, s3 * (MM_LOAD(g+19*SIMDD) + b01) + b11 * r3);
        MM_STORE(g+34*SIMDD, s4 * MM_LOAD(g+30*SIMDD) + b10 * MM_LOAD(g+26*SIMDD) + b00 * MM_LOAD(g+28*SIMDD));
        MM_STORE(g+35*SIMDD, s5 * MM_LOAD(g+31*SIMDD) + b11 * MM_LOAD(g+27*SIMDD) + b01 * MM_LOAD(g+29*SIMDD));
}

static inline void _g0_2d4d_1010(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *c0x = bc->c00x;
        double *c0y = bc->c00y;
        double *c0z = bc->c00z;
        double *cpx = bc->c0px;
        double *cpy = bc->c0py;
        double *cpz = bc->c0pz;
        double *b   = bc->b00;
        __MD cx0 = MM_LOAD(c0x+0*SIMDD);
        __MD cx1 = MM_LOAD(c0x+1*SIMDD);
        __MD cy0 = MM_LOAD(c0y+0*SIMDD);
        __MD cy1 = MM_LOAD(c0y+1*SIMDD);
        __MD cz0 = MM_LOAD(c0z+0*SIMDD);
        __MD cz1 = MM_LOAD(c0z+1*SIMDD);
        __MD px0 = MM_LOAD(cpx+0*SIMDD);
        __MD px1 = MM_LOAD(cpx+1*SIMDD);
        __MD py0 = MM_LOAD(cpy+0*SIMDD);
        __MD py1 = MM_LOAD(cpy+1*SIMDD);
        __MD pz0 = MM_LOAD(cpz+0*SIMDD);
        __MD pz1 = MM_LOAD(cpz+1*SIMDD);
        __MD b0 = MM_LOAD(b+0*SIMDD);
        __MD b1 = MM_LOAD(b+1*SIMDD);
        __MD g16 = MM_LOAD(g+16*SIMDD);
        __MD g17 = MM_LOAD(g+17*SIMDD);
        MM_STORE(g+2 *SIMDD, px0);
        MM_STORE(g+3 *SIMDD, px1);
        MM_STORE(g+4 *SIMDD, cx0);
        MM_STORE(g+5 *SIMDD, cx1);
        MM_STORE(g+10*SIMDD, py0);
        MM_STORE(g+11*SIMDD, py1);
        MM_STORE(g+12*SIMDD, cy0);
        MM_STORE(g+13*SIMDD, cy1);
        MM_STORE(g+18*SIMDD, MM_MUL(pz0, g16));
        MM_STORE(g+19*SIMDD, MM_MUL(pz1, g17));
        MM_STORE(g+20*SIMDD, MM_MUL(cz0, g16));
        MM_STORE(g+21*SIMDD, MM_MUL(cz1, g17));
        MM_STORE(g+6 *SIMDD, MM_FMA(px0, cx0, b0));
        MM_STORE(g+7 *SIMDD, MM_FMA(px1, cx1, b1));
        MM_STORE(g+14*SIMDD, MM_FMA(py0, cy0, b0));
        MM_STORE(g+15*SIMDD, MM_FMA(py1, cy1, b1));
        MM_STORE(g+22*SIMDD, MM_MUL(MM_FMA(pz0, cz0, b0), g16));
        MM_STORE(g+23*SIMDD, MM_MUL(MM_FMA(pz1, cz1, b1), g17));
}

static inline void _g0_2d4d_1011(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *c0x = bc->c00x;
        double *c0y = bc->c00y;
        double *c0z = bc->c00z;
        double *cpx = bc->c0px;
        double *cpy = bc->c0py;
        double *cpz = bc->c0pz;
        double *b0  = bc->b00;
        double *b1  = bc->b01;
        double *rkl = envs->rkrl;
        ALIGNMM double rcp[6*SIMDD];
        _make_rc(rcp, cpx, cpy, cpz, rkl);
        __MD b00 = MM_LOAD(b0+0*SIMDD);
        __MD b01 = MM_LOAD(b0+1*SIMDD);
        __MD b10 = MM_LOAD(b1+0*SIMDD);
        __MD b11 = MM_LOAD(b1+1*SIMDD);
        __MD g48 = MM_LOAD(g+48*SIMDD);
        __MD g49 = MM_LOAD(g+49*SIMDD);
        MM_STORE(g+2 *SIMDD, MM_LOAD(c0x+0*SIMDD));
        MM_STORE(g+3 *SIMDD, MM_LOAD(c0x+1*SIMDD));
        MM_STORE(g+4 *SIMDD, MM_LOAD(rcp+0*SIMDD));
        MM_STORE(g+5 *SIMDD, MM_LOAD(rcp+1*SIMDD));
        MM_STORE(g+8 *SIMDD, MM_LOAD(cpx+0*SIMDD));
        MM_STORE(g+9 *SIMDD, MM_LOAD(cpx+1*SIMDD));
        MM_STORE(g+26*SIMDD, MM_LOAD(c0y+0*SIMDD));
        MM_STORE(g+27*SIMDD, MM_LOAD(c0y+1*SIMDD));
        MM_STORE(g+28*SIMDD, MM_LOAD(rcp+2*SIMDD));
        MM_STORE(g+29*SIMDD, MM_LOAD(rcp+3*SIMDD));
        MM_STORE(g+32*SIMDD, MM_LOAD(cpy+0*SIMDD));
        MM_STORE(g+33*SIMDD, MM_LOAD(cpy+1*SIMDD));
        MM_STORE(g+50*SIMDD, MM_LOAD(c0z+0*SIMDD) * g48);
        MM_STORE(g+51*SIMDD, MM_LOAD(c0z+1*SIMDD) * g49);
        MM_STORE(g+52*SIMDD, MM_LOAD(rcp+4*SIMDD) * g48);
        MM_STORE(g+53*SIMDD, MM_LOAD(rcp+5*SIMDD) * g49);
        MM_STORE(g+56*SIMDD, MM_LOAD(cpz+0*SIMDD) * g48);
        MM_STORE(g+57*SIMDD, MM_LOAD(cpz+1*SIMDD) * g49);
        MM_STORE(g+6 *SIMDD , MM_LOAD(rcp+0*SIMDD) * MM_LOAD(c0x+0*SIMDD) + b00);
        MM_STORE(g+7 *SIMDD , MM_LOAD(rcp+1*SIMDD) * MM_LOAD(c0x+1*SIMDD) + b01);
        MM_STORE(g+10*SIMDD , MM_LOAD(cpx+0*SIMDD) * MM_LOAD(c0x+0*SIMDD) + b00);
        MM_STORE(g+11*SIMDD , MM_LOAD(cpx+1*SIMDD) * MM_LOAD(c0x+1*SIMDD) + b01);
        MM_STORE(g+12*SIMDD , MM_LOAD(cpx+0*SIMDD) * MM_LOAD(rcp+0*SIMDD) + b10);
        MM_STORE(g+13*SIMDD , MM_LOAD(cpx+1*SIMDD) * MM_LOAD(rcp+1*SIMDD) + b11);
        MM_STORE(g+30*SIMDD , MM_LOAD(rcp+2*SIMDD) * MM_LOAD(c0y+0*SIMDD) + b00);
        MM_STORE(g+31*SIMDD , MM_LOAD(rcp+3*SIMDD) * MM_LOAD(c0y+1*SIMDD) + b01);
        MM_STORE(g+34*SIMDD , MM_LOAD(cpy+0*SIMDD) * MM_LOAD(c0y+0*SIMDD) + b00);
        MM_STORE(g+35*SIMDD , MM_LOAD(cpy+1*SIMDD) * MM_LOAD(c0y+1*SIMDD) + b01);
        MM_STORE(g+36*SIMDD , MM_LOAD(cpy+0*SIMDD) * MM_LOAD(rcp+2*SIMDD) + b10);
        MM_STORE(g+37*SIMDD , MM_LOAD(cpy+1*SIMDD) * MM_LOAD(rcp+3*SIMDD) + b11);
        MM_STORE(g+54*SIMDD,(MM_LOAD(rcp+4*SIMDD) * MM_LOAD(c0z+0*SIMDD) + b00)* g48);
        MM_STORE(g+55*SIMDD,(MM_LOAD(rcp+5*SIMDD) * MM_LOAD(c0z+1*SIMDD) + b01)* g49);
        MM_STORE(g+58*SIMDD,(MM_LOAD(cpz+0*SIMDD) * MM_LOAD(c0z+0*SIMDD) + b00)* g48);
        MM_STORE(g+59*SIMDD,(MM_LOAD(cpz+1*SIMDD) * MM_LOAD(c0z+1*SIMDD) + b01)* g49);
        MM_STORE(g+60*SIMDD,(MM_LOAD(cpz+0*SIMDD) * MM_LOAD(rcp+4*SIMDD) + b10)* g48);
        MM_STORE(g+61*SIMDD,(MM_LOAD(cpz+1*SIMDD) * MM_LOAD(rcp+5*SIMDD) + b11)* g49);
        MM_STORE(g+14*SIMDD, MM_LOAD(rcp+0*SIMDD) * MM_LOAD(g+10*SIMDD) + b00 * MM_LOAD(cpx+0*SIMDD) + b10 * MM_LOAD(c0x+0*SIMDD));
        MM_STORE(g+15*SIMDD, MM_LOAD(rcp+1*SIMDD) * MM_LOAD(g+11*SIMDD) + b01 * MM_LOAD(cpx+1*SIMDD) + b11 * MM_LOAD(c0x+1*SIMDD));
        MM_STORE(g+38*SIMDD, MM_LOAD(rcp+2*SIMDD) * MM_LOAD(g+34*SIMDD) + b00 * MM_LOAD(cpy+0*SIMDD) + b10 * MM_LOAD(c0y+0*SIMDD));
        MM_STORE(g+39*SIMDD, MM_LOAD(rcp+3*SIMDD) * MM_LOAD(g+35*SIMDD) + b01 * MM_LOAD(cpy+1*SIMDD) + b11 * MM_LOAD(c0y+1*SIMDD));
        MM_STORE(g+62*SIMDD, MM_LOAD(rcp+4*SIMDD) * MM_LOAD(g+58*SIMDD) + b00 * MM_LOAD(g+56*SIMDD) + b10 * MM_LOAD(g+50*SIMDD));
        MM_STORE(g+63*SIMDD, MM_LOAD(rcp+5*SIMDD) * MM_LOAD(g+59*SIMDD) + b01 * MM_LOAD(g+57*SIMDD) + b11 * MM_LOAD(g+51*SIMDD));
}

static inline void _g0_2d4d_1020(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *c0x = bc->c00x;
        double *c0y = bc->c00y;
        double *c0z = bc->c00z;
        double *cpx = bc->c0px;
        double *cpy = bc->c0py;
        double *cpz = bc->c0pz;
        double *b0 = bc->b00;
        double *b1 = bc->b01;
        __MD r0 = MM_LOAD(c0x+0*SIMDD);
        __MD r1 = MM_LOAD(c0x+1*SIMDD);
        __MD r2 = MM_LOAD(c0y+0*SIMDD);
        __MD r3 = MM_LOAD(c0y+1*SIMDD);
        __MD r4 = MM_LOAD(c0z+0*SIMDD);
        __MD r5 = MM_LOAD(c0z+1*SIMDD);
        __MD s0 = MM_LOAD(cpx+0*SIMDD);
        __MD s1 = MM_LOAD(cpx+1*SIMDD);
        __MD s2 = MM_LOAD(cpy+0*SIMDD);
        __MD s3 = MM_LOAD(cpy+1*SIMDD);
        __MD s4 = MM_LOAD(cpz+0*SIMDD);
        __MD s5 = MM_LOAD(cpz+1*SIMDD);
        __MD b00 = MM_LOAD(b0+0*SIMDD);
        __MD b01 = MM_LOAD(b0+1*SIMDD);
        __MD b10 = MM_LOAD(b1+0*SIMDD);
        __MD b11 = MM_LOAD(b1+1*SIMDD);
        __MD g24 = MM_LOAD(g+24*SIMDD);
        __MD g25 = MM_LOAD(g+25*SIMDD);
        MM_STORE(g+2 *SIMDD, r0);
        MM_STORE(g+3 *SIMDD, r1);
        MM_STORE(g+4 *SIMDD, s0);
        MM_STORE(g+5 *SIMDD, s1);
        MM_STORE(g+14*SIMDD, r2);
        MM_STORE(g+15*SIMDD, r3);
        MM_STORE(g+16*SIMDD, s2);
        MM_STORE(g+17*SIMDD, s3);
        MM_STORE(g+26*SIMDD, MM_MUL(r4, g24));
        MM_STORE(g+27*SIMDD, MM_MUL(r5, g25));
        MM_STORE(g+28*SIMDD, MM_MUL(s4, g24));
        MM_STORE(g+29*SIMDD, MM_MUL(s5, g25));
        MM_STORE(g+6 *SIMDD, MM_FMA(r0, s0, b00));
        MM_STORE(g+7 *SIMDD, MM_FMA(r1, s1, b01));
        MM_STORE(g+8 *SIMDD, MM_FMA(s0, s0, b10));
        MM_STORE(g+9 *SIMDD, MM_FMA(s1, s1, b11));
        MM_STORE(g+18*SIMDD, MM_FMA(r2, s2, b00));
        MM_STORE(g+19*SIMDD, MM_FMA(r3, s3, b01));
        MM_STORE(g+20*SIMDD, MM_FMA(s2, s2, b10));
        MM_STORE(g+21*SIMDD, MM_FMA(s3, s3, b11));
        MM_STORE(g+30*SIMDD, MM_MUL(MM_FMA(r4, s4, b00), g24));
        MM_STORE(g+31*SIMDD, MM_MUL(MM_FMA(r5, s5, b01), g25));
        MM_STORE(g+32*SIMDD, MM_MUL(MM_FMA(s4, s4, b10), g24));
        MM_STORE(g+33*SIMDD, MM_MUL(MM_FMA(s5, s5, b11), g25));
        MM_STORE(g+10*SIMDD, s0 * (MM_LOAD(g+6 *SIMDD) + b00) + b10 * r0);
        MM_STORE(g+11*SIMDD, s1 * (MM_LOAD(g+7 *SIMDD) + b01) + b11 * r1);
        MM_STORE(g+22*SIMDD, s2 * (MM_LOAD(g+18*SIMDD) + b00) + b10 * r2);
        MM_STORE(g+23*SIMDD, s3 * (MM_LOAD(g+19*SIMDD) + b01) + b11 * r3);
        MM_STORE(g+34*SIMDD, s4 * MM_LOAD(g+30*SIMDD) + b10 * MM_LOAD(g+26*SIMDD) + b00 * MM_LOAD(g+28*SIMDD));
        MM_STORE(g+35*SIMDD, s5 * MM_LOAD(g+31*SIMDD) + b11 * MM_LOAD(g+27*SIMDD) + b01 * MM_LOAD(g+29*SIMDD));
}

static inline void _g0_2d4d_1100(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *cx = bc->c00x;
        double *cy = bc->c00y;
        double *cz = bc->c00z;
        double *b  = bc->b10;
        double *r  = envs->rirj;
        ALIGNMM double rc[6*SIMDD];
        _make_rc(rc, cx, cy, cz, r);
        __MD r0 = MM_LOAD(rc+0*SIMDD);
        __MD r1 = MM_LOAD(rc+1*SIMDD);
        __MD r2 = MM_LOAD(rc+2*SIMDD);
        __MD r3 = MM_LOAD(rc+3*SIMDD);
        __MD r4 = MM_LOAD(rc+4*SIMDD);
        __MD r5 = MM_LOAD(rc+5*SIMDD);
        __MD cx0 = MM_LOAD(cx+0*SIMDD);
        __MD cx1 = MM_LOAD(cx+1*SIMDD);
        __MD cy0 = MM_LOAD(cy+0*SIMDD);
        __MD cy1 = MM_LOAD(cy+1*SIMDD);
        __MD cz0 = MM_LOAD(cz+0*SIMDD);
        __MD cz1 = MM_LOAD(cz+1*SIMDD);
        __MD b0 = MM_LOAD(b+0*SIMDD);
        __MD b1 = MM_LOAD(b+1*SIMDD);
        __MD g24 = MM_LOAD(g+24*SIMDD);
        __MD g25 = MM_LOAD(g+25*SIMDD);
        MM_STORE(g+2 *SIMDD, r0);
        MM_STORE(g+3 *SIMDD, r1);
        MM_STORE(g+4 *SIMDD, cx0);
        MM_STORE(g+5 *SIMDD, cx1);
        MM_STORE(g+14*SIMDD, r2);
        MM_STORE(g+15*SIMDD, r3);
        MM_STORE(g+16*SIMDD, cy0);
        MM_STORE(g+17*SIMDD, cy1);
        MM_STORE(g+26*SIMDD, MM_MUL(r4, g24));
        MM_STORE(g+27*SIMDD, MM_MUL(r5, g25));
        MM_STORE(g+28*SIMDD, MM_MUL(cz0,g24));
        MM_STORE(g+29*SIMDD, MM_MUL(cz1,g25));
        MM_STORE(g+6 *SIMDD, MM_FMA(r0, cx0, b0));
        MM_STORE(g+7 *SIMDD, MM_FMA(r1, cx1, b1));
        MM_STORE(g+18*SIMDD, MM_FMA(r2, cy0, b0));
        MM_STORE(g+19*SIMDD, MM_FMA(r3, cy1, b1));
        MM_STORE(g+30*SIMDD, MM_MUL(MM_FMA(r4, cz0, b0), g24));
        MM_STORE(g+31*SIMDD, MM_MUL(MM_FMA(r5, cz1, b1), g25));
}

static inline void _g0_2d4d_1101(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *c0x = bc->c00x;
        double *c0y = bc->c00y;
        double *c0z = bc->c00z;
        double *cpx = bc->c0px;
        double *cpy = bc->c0py;
        double *cpz = bc->c0pz;
        double *b0  = bc->b00;
        double *b1  = bc->b10;
        double *rij = envs->rirj;
        ALIGNMM double rc0[6*SIMDD];
        _make_rc(rc0, c0x, c0y, c0z, rij);
        __MD b00 = MM_LOAD(b0+0*SIMDD);
        __MD b01 = MM_LOAD(b0+1*SIMDD);
        __MD b10 = MM_LOAD(b1+0*SIMDD);
        __MD b11 = MM_LOAD(b1+1*SIMDD);
        __MD g48 = MM_LOAD(g+48*SIMDD);
        __MD g49 = MM_LOAD(g+49*SIMDD);
        MM_STORE(g+2 *SIMDD, MM_LOAD(rc0+0*SIMDD));
        MM_STORE(g+3 *SIMDD, MM_LOAD(rc0+1*SIMDD));
        MM_STORE(g+4 *SIMDD, MM_LOAD(cpx+0*SIMDD));
        MM_STORE(g+5 *SIMDD, MM_LOAD(cpx+1*SIMDD));
        MM_STORE(g+8 *SIMDD, MM_LOAD(c0x+0*SIMDD));
        MM_STORE(g+9 *SIMDD, MM_LOAD(c0x+1*SIMDD));
        MM_STORE(g+26*SIMDD, MM_LOAD(rc0+2*SIMDD));
        MM_STORE(g+27*SIMDD, MM_LOAD(rc0+3*SIMDD));
        MM_STORE(g+28*SIMDD, MM_LOAD(cpy+0*SIMDD));
        MM_STORE(g+29*SIMDD, MM_LOAD(cpy+1*SIMDD));
        MM_STORE(g+32*SIMDD, MM_LOAD(c0y+0*SIMDD));
        MM_STORE(g+33*SIMDD, MM_LOAD(c0y+1*SIMDD));
        MM_STORE(g+50*SIMDD, MM_LOAD(rc0+4*SIMDD) * g48);
        MM_STORE(g+51*SIMDD, MM_LOAD(rc0+5*SIMDD) * g49);
        MM_STORE(g+52*SIMDD, MM_LOAD(cpz+0*SIMDD) * g48);
        MM_STORE(g+53*SIMDD, MM_LOAD(cpz+1*SIMDD) * g49);
        MM_STORE(g+56*SIMDD, MM_LOAD(c0z+0*SIMDD) * g48);
        MM_STORE(g+57*SIMDD, MM_LOAD(c0z+1*SIMDD) * g49);
        MM_STORE(g+6 *SIMDD, MM_LOAD(cpx+0*SIMDD) * MM_LOAD(rc0+0*SIMDD) + b00);
        MM_STORE(g+7 *SIMDD, MM_LOAD(cpx+1*SIMDD) * MM_LOAD(rc0+1*SIMDD) + b01);
        MM_STORE(g+10*SIMDD, MM_LOAD(c0x+0*SIMDD) * MM_LOAD(rc0+0*SIMDD) + b10);
        MM_STORE(g+11*SIMDD, MM_LOAD(c0x+1*SIMDD) * MM_LOAD(rc0+1*SIMDD) + b11);
        MM_STORE(g+12*SIMDD, MM_LOAD(c0x+0*SIMDD) * MM_LOAD(cpx+0*SIMDD) + b00);
        MM_STORE(g+13*SIMDD, MM_LOAD(c0x+1*SIMDD) * MM_LOAD(cpx+1*SIMDD) + b01);
        MM_STORE(g+30*SIMDD, MM_LOAD(cpy+0*SIMDD) * MM_LOAD(rc0+2*SIMDD) + b00);
        MM_STORE(g+31*SIMDD, MM_LOAD(cpy+1*SIMDD) * MM_LOAD(rc0+3*SIMDD) + b01);
        MM_STORE(g+34*SIMDD, MM_LOAD(c0y+0*SIMDD) * MM_LOAD(rc0+2*SIMDD) + b10);
        MM_STORE(g+35*SIMDD, MM_LOAD(c0y+1*SIMDD) * MM_LOAD(rc0+3*SIMDD) + b11);
        MM_STORE(g+36*SIMDD, MM_LOAD(c0y+0*SIMDD) * MM_LOAD(cpy+0*SIMDD) + b00);
        MM_STORE(g+37*SIMDD, MM_LOAD(c0y+1*SIMDD) * MM_LOAD(cpy+1*SIMDD) + b01);
        MM_STORE(g+54*SIMDD,(MM_LOAD(cpz+0*SIMDD) * MM_LOAD(rc0+4*SIMDD) + b00)* g48);
        MM_STORE(g+55*SIMDD,(MM_LOAD(cpz+1*SIMDD) * MM_LOAD(rc0+5*SIMDD) + b01)* g49);
        MM_STORE(g+58*SIMDD,(MM_LOAD(c0z+0*SIMDD) * MM_LOAD(rc0+4*SIMDD) + b10)* g48);
        MM_STORE(g+59*SIMDD,(MM_LOAD(c0z+1*SIMDD) * MM_LOAD(rc0+5*SIMDD) + b11)* g49);
        MM_STORE(g+60*SIMDD,(MM_LOAD(c0z+0*SIMDD) * MM_LOAD(cpz+0*SIMDD) + b00)* g48);
        MM_STORE(g+61*SIMDD,(MM_LOAD(c0z+1*SIMDD) * MM_LOAD(cpz+1*SIMDD) + b01)* g49);
        MM_STORE(g+14*SIMDD, MM_LOAD(cpx+0*SIMDD) * MM_LOAD(g+10*SIMDD ) + b00 *(MM_LOAD(rc0+0*SIMDD) + MM_LOAD(c0x+0*SIMDD)));
        MM_STORE(g+15*SIMDD, MM_LOAD(cpx+1*SIMDD) * MM_LOAD(g+11*SIMDD ) + b01 *(MM_LOAD(rc0+1*SIMDD) + MM_LOAD(c0x+1*SIMDD)));
        MM_STORE(g+38*SIMDD, MM_LOAD(cpy+0*SIMDD) * MM_LOAD(g+34*SIMDD ) + b00 *(MM_LOAD(rc0+2*SIMDD) + MM_LOAD(c0y+0*SIMDD)));
        MM_STORE(g+39*SIMDD, MM_LOAD(cpy+1*SIMDD) * MM_LOAD(g+35*SIMDD ) + b01 *(MM_LOAD(rc0+3*SIMDD) + MM_LOAD(c0y+1*SIMDD)));
        MM_STORE(g+62*SIMDD, MM_LOAD(cpz+0*SIMDD) * MM_LOAD(g+58*SIMDD ) + b00 *(MM_LOAD(g+50*SIMDD ) + MM_LOAD(g+56*SIMDD )));
        MM_STORE(g+63*SIMDD, MM_LOAD(cpz+1*SIMDD) * MM_LOAD(g+59*SIMDD ) + b01 *(MM_LOAD(g+51*SIMDD ) + MM_LOAD(g+57*SIMDD )));
}

static inline void _g0_2d4d_1110(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *c0x = bc->c00x;
        double *c0y = bc->c00y;
        double *c0z = bc->c00z;
        double *cpx = bc->c0px;
        double *cpy = bc->c0py;
        double *cpz = bc->c0pz;
        double *b0  = bc->b00;
        double *b1  = bc->b10;
        double *r0  = envs->rirj;
        ALIGNMM double rc0[6*SIMDD];
        _make_rc(rc0, c0x, c0y, c0z, r0);
        __MD b00 = MM_LOAD(b0+0*SIMDD);
        __MD b01 = MM_LOAD(b0+1*SIMDD);
        __MD b10 = MM_LOAD(b1+0*SIMDD);
        __MD b11 = MM_LOAD(b1+1*SIMDD);
        __MD g48 = MM_LOAD(g+48*SIMDD);
        __MD g49 = MM_LOAD(g+49*SIMDD);
        MM_STORE(g+2 *SIMDD, MM_LOAD(rc0+0*SIMDD));
        MM_STORE(g+3 *SIMDD, MM_LOAD(rc0+1*SIMDD));
        MM_STORE(g+4 *SIMDD, MM_LOAD(cpx+0*SIMDD));
        MM_STORE(g+5 *SIMDD, MM_LOAD(cpx+1*SIMDD));
        MM_STORE(g+8 *SIMDD, MM_LOAD(c0x+0*SIMDD));
        MM_STORE(g+9 *SIMDD, MM_LOAD(c0x+1*SIMDD));
        MM_STORE(g+26*SIMDD, MM_LOAD(rc0+2*SIMDD));
        MM_STORE(g+27*SIMDD, MM_LOAD(rc0+3*SIMDD));
        MM_STORE(g+28*SIMDD, MM_LOAD(cpy+0*SIMDD));
        MM_STORE(g+29*SIMDD, MM_LOAD(cpy+1*SIMDD));
        MM_STORE(g+32*SIMDD, MM_LOAD(c0y+0*SIMDD));
        MM_STORE(g+33*SIMDD, MM_LOAD(c0y+1*SIMDD));
        MM_STORE(g+50*SIMDD, MM_LOAD(rc0+4*SIMDD) * g48);
        MM_STORE(g+51*SIMDD, MM_LOAD(rc0+5*SIMDD) * g49);
        MM_STORE(g+52*SIMDD, MM_LOAD(cpz+0*SIMDD) * g48);
        MM_STORE(g+53*SIMDD, MM_LOAD(cpz+1*SIMDD) * g49);
        MM_STORE(g+56*SIMDD, MM_LOAD(c0z+0*SIMDD) * g48);
        MM_STORE(g+57*SIMDD, MM_LOAD(c0z+1*SIMDD) * g49);
        MM_STORE(g+6 *SIMDD, MM_LOAD(cpx+0*SIMDD) * MM_LOAD(rc0+0*SIMDD) + b00);
        MM_STORE(g+7 *SIMDD, MM_LOAD(cpx+1*SIMDD) * MM_LOAD(rc0+1*SIMDD) + b01);
        MM_STORE(g+10*SIMDD, MM_LOAD(c0x+0*SIMDD) * MM_LOAD(rc0+0*SIMDD) + b10);
        MM_STORE(g+11*SIMDD, MM_LOAD(c0x+1*SIMDD) * MM_LOAD(rc0+1*SIMDD) + b11);
        MM_STORE(g+12*SIMDD, MM_LOAD(c0x+0*SIMDD) * MM_LOAD(cpx+0*SIMDD) + b00);
        MM_STORE(g+13*SIMDD, MM_LOAD(c0x+1*SIMDD) * MM_LOAD(cpx+1*SIMDD) + b01);
        MM_STORE(g+30*SIMDD, MM_LOAD(cpy+0*SIMDD) * MM_LOAD(rc0+2*SIMDD) + b00);
        MM_STORE(g+31*SIMDD, MM_LOAD(cpy+1*SIMDD) * MM_LOAD(rc0+3*SIMDD) + b01);
        MM_STORE(g+34*SIMDD, MM_LOAD(c0y+0*SIMDD) * MM_LOAD(rc0+2*SIMDD) + b10);
        MM_STORE(g+35*SIMDD, MM_LOAD(c0y+1*SIMDD) * MM_LOAD(rc0+3*SIMDD) + b11);
        MM_STORE(g+36*SIMDD, MM_LOAD(c0y+0*SIMDD) * MM_LOAD(cpy+0*SIMDD) + b00);
        MM_STORE(g+37*SIMDD, MM_LOAD(c0y+1*SIMDD) * MM_LOAD(cpy+1*SIMDD) + b01);
        MM_STORE(g+54*SIMDD,(MM_LOAD(cpz+0*SIMDD) * MM_LOAD(rc0+4*SIMDD) + b00)* g48);
        MM_STORE(g+55*SIMDD,(MM_LOAD(cpz+1*SIMDD) * MM_LOAD(rc0+5*SIMDD) + b01)* g49);
        MM_STORE(g+58*SIMDD,(MM_LOAD(c0z+0*SIMDD) * MM_LOAD(rc0+4*SIMDD) + b10)* g48);
        MM_STORE(g+59*SIMDD,(MM_LOAD(c0z+1*SIMDD) * MM_LOAD(rc0+5*SIMDD) + b11)* g49);
        MM_STORE(g+60*SIMDD,(MM_LOAD(c0z+0*SIMDD) * MM_LOAD(cpz+0*SIMDD) + b00)* g48);
        MM_STORE(g+61*SIMDD,(MM_LOAD(c0z+1*SIMDD) * MM_LOAD(cpz+1*SIMDD) + b01)* g49);
        MM_STORE(g+14*SIMDD, MM_LOAD(cpx+0*SIMDD) * MM_LOAD(g+10*SIMDD ) + b00 *(MM_LOAD(rc0+0*SIMDD) + MM_LOAD(c0x+0*SIMDD)));
        MM_STORE(g+15*SIMDD, MM_LOAD(cpx+1*SIMDD) * MM_LOAD(g+11*SIMDD ) + b01 *(MM_LOAD(rc0+1*SIMDD) + MM_LOAD(c0x+1*SIMDD)));
        MM_STORE(g+38*SIMDD, MM_LOAD(cpy+0*SIMDD) * MM_LOAD(g+34*SIMDD ) + b00 *(MM_LOAD(rc0+2*SIMDD) + MM_LOAD(c0y+0*SIMDD)));
        MM_STORE(g+39*SIMDD, MM_LOAD(cpy+1*SIMDD) * MM_LOAD(g+35*SIMDD ) + b01 *(MM_LOAD(rc0+3*SIMDD) + MM_LOAD(c0y+1*SIMDD)));
        MM_STORE(g+62*SIMDD, MM_LOAD(cpz+0*SIMDD) * MM_LOAD(g+58*SIMDD ) + b00 *(MM_LOAD(g+50*SIMDD ) + MM_LOAD(g+56*SIMDD )));
        MM_STORE(g+63*SIMDD, MM_LOAD(cpz+1*SIMDD) * MM_LOAD(g+59*SIMDD ) + b01 *(MM_LOAD(g+51*SIMDD ) + MM_LOAD(g+57*SIMDD )));
}

static inline void _g0_2d4d_1200(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *cx = bc->c00x;
        double *cy = bc->c00y;
        double *cz = bc->c00z;
        double *b  = bc->b10;
        double *r  = envs->rirj;
        ALIGNMM double rc[6*SIMDD];
        _make_rc(rc, cx, cy, cz, r);
        __MD r0 = MM_LOAD(rc+0*SIMDD);
        __MD r1 = MM_LOAD(rc+1*SIMDD);
        __MD r2 = MM_LOAD(rc+2*SIMDD);
        __MD r3 = MM_LOAD(rc+3*SIMDD);
        __MD r4 = MM_LOAD(rc+4*SIMDD);
        __MD r5 = MM_LOAD(rc+5*SIMDD);
        __MD cx0 = MM_LOAD(cx+0*SIMDD);
        __MD cx1 = MM_LOAD(cx+1*SIMDD);
        __MD cy0 = MM_LOAD(cy+0*SIMDD);
        __MD cy1 = MM_LOAD(cy+1*SIMDD);
        __MD cz0 = MM_LOAD(cz+0*SIMDD);
        __MD cz1 = MM_LOAD(cz+1*SIMDD);
        __MD b0 = MM_LOAD(b+0*SIMDD);
        __MD b1 = MM_LOAD(b+1*SIMDD);
        __MD i2 = MM_SET1(2.);
        __MD g32 = MM_LOAD(g+32*SIMDD);
        __MD g33 = MM_LOAD(g+33*SIMDD);
        MM_STORE(g+2 *SIMDD, r0 );
        MM_STORE(g+3 *SIMDD, r1 );
        MM_STORE(g+4 *SIMDD, cx0);
        MM_STORE(g+5 *SIMDD, cx1);
        MM_STORE(g+18*SIMDD, r2 );
        MM_STORE(g+19*SIMDD, r3 );
        MM_STORE(g+20*SIMDD, cy0);
        MM_STORE(g+21*SIMDD, cy1);
        MM_STORE(g+34*SIMDD, MM_MUL(r4 , g32));
        MM_STORE(g+35*SIMDD, MM_MUL(r5 , g33));
        MM_STORE(g+36*SIMDD, MM_MUL(cz0, g32));
        MM_STORE(g+37*SIMDD, MM_MUL(cz1, g33));
        MM_STORE(g+6 *SIMDD, MM_FMA(r0 , cx0, b0));
        MM_STORE(g+7 *SIMDD, MM_FMA(r1 , cx1, b1));
        MM_STORE(g+8 *SIMDD, MM_FMA(cx0, cx0, b0));
        MM_STORE(g+9 *SIMDD, MM_FMA(cx1, cx1, b1));
        MM_STORE(g+22*SIMDD, MM_FMA(r2 , cy0, b0));
        MM_STORE(g+23*SIMDD, MM_FMA(r3 , cy1, b1));
        MM_STORE(g+24*SIMDD, MM_FMA(cy0, cy0, b0));
        MM_STORE(g+25*SIMDD, MM_FMA(cy1, cy1, b1));
        MM_STORE(g+38*SIMDD, MM_MUL(MM_FMA(r4 , cz0, b0), g32));
        MM_STORE(g+39*SIMDD, MM_MUL(MM_FMA(r5 , cz1, b1), g33));
        MM_STORE(g+40*SIMDD, MM_MUL(MM_FMA(cz0, cz0, b0), g32));
        MM_STORE(g+41*SIMDD, MM_MUL(MM_FMA(cz1, cz1, b1), g33));
        MM_STORE(g+10*SIMDD, r0 * MM_LOAD(g+8 *SIMDD) + i2 * b0 * cx0                );
        MM_STORE(g+11*SIMDD, r1 * MM_LOAD(g+9 *SIMDD) + i2 * b1 * cx1                );
        MM_STORE(g+26*SIMDD, r2 * MM_LOAD(g+24*SIMDD) + i2 * b0 * cy0                );
        MM_STORE(g+27*SIMDD, r3 * MM_LOAD(g+25*SIMDD) + i2 * b1 * cy1                );
        MM_STORE(g+42*SIMDD, r4 * MM_LOAD(g+40*SIMDD) + i2 * b0 * MM_LOAD(g+36*SIMDD));
        MM_STORE(g+43*SIMDD, r5 * MM_LOAD(g+41*SIMDD) + i2 * b1 * MM_LOAD(g+37*SIMDD));
}

static inline void _g0_2d4d_2000(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *cx = bc->c00x;
        double *cy = bc->c00y;
        double *cz = bc->c00z;
        double *b  = bc->b10;
        __MD r0 = MM_LOAD(cx+0*SIMDD);
        __MD r1 = MM_LOAD(cx+1*SIMDD);
        __MD r2 = MM_LOAD(cy+0*SIMDD);
        __MD r3 = MM_LOAD(cy+1*SIMDD);
        __MD r4 = MM_LOAD(cz+0*SIMDD);
        __MD r5 = MM_LOAD(cz+1*SIMDD);
        __MD b0 = MM_LOAD(b+0*SIMDD);
        __MD b1 = MM_LOAD(b+1*SIMDD);
        __MD g12 = MM_LOAD(g+12*SIMDD);
        __MD g13 = MM_LOAD(g+13*SIMDD);
        MM_STORE(g+2 *SIMDD, r0);
        MM_STORE(g+3 *SIMDD, r1);
        MM_STORE(g+8 *SIMDD, r2);
        MM_STORE(g+9 *SIMDD, r3);
        MM_STORE(g+14*SIMDD, MM_MUL(r4, g12));
        MM_STORE(g+15*SIMDD, MM_MUL(r5, g13));
        MM_STORE(g+4 *SIMDD, MM_FMA(r0, r0, b0));
        MM_STORE(g+5 *SIMDD, MM_FMA(r1, r1, b1));
        MM_STORE(g+10*SIMDD, MM_FMA(r2, r2, b0));
        MM_STORE(g+11*SIMDD, MM_FMA(r3, r3, b1));
        MM_STORE(g+16*SIMDD, MM_MUL(MM_FMA(r4, r4, b0), g12));
        MM_STORE(g+17*SIMDD, MM_MUL(MM_FMA(r5, r5, b1), g13));
}

static inline void _g0_2d4d_2001(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *c0x = bc->c00x;
        double *c0y = bc->c00y;
        double *c0z = bc->c00z;
        double *cpx = bc->c0px;
        double *cpy = bc->c0py;
        double *cpz = bc->c0pz;
        double *b0 = bc->b00;
        double *b1 = bc->b01;
        __MD r0 = MM_LOAD(c0x+0*SIMDD);
        __MD r1 = MM_LOAD(c0x+1*SIMDD);
        __MD r2 = MM_LOAD(c0y+0*SIMDD);
        __MD r3 = MM_LOAD(c0y+1*SIMDD);
        __MD r4 = MM_LOAD(c0z+0*SIMDD);
        __MD r5 = MM_LOAD(c0z+1*SIMDD);
        __MD s0 = MM_LOAD(cpx+0*SIMDD);
        __MD s1 = MM_LOAD(cpx+1*SIMDD);
        __MD s2 = MM_LOAD(cpy+0*SIMDD);
        __MD s3 = MM_LOAD(cpy+1*SIMDD);
        __MD s4 = MM_LOAD(cpz+0*SIMDD);
        __MD s5 = MM_LOAD(cpz+1*SIMDD);
        __MD b00 = MM_LOAD(b0+0*SIMDD);
        __MD b01 = MM_LOAD(b0+1*SIMDD);
        __MD b10 = MM_LOAD(b1+0*SIMDD);
        __MD b11 = MM_LOAD(b1+1*SIMDD);
        __MD g24 = MM_LOAD(g+24*SIMDD);
        __MD g25 = MM_LOAD(g+25*SIMDD);
        MM_STORE(g+2 *SIMDD, r0);
        MM_STORE(g+3 *SIMDD, r1);
        MM_STORE(g+6 *SIMDD, s0);
        MM_STORE(g+7 *SIMDD, s1);
        MM_STORE(g+14*SIMDD, r2);
        MM_STORE(g+15*SIMDD, r3);
        MM_STORE(g+18*SIMDD, s2);
        MM_STORE(g+19*SIMDD, s3);
        MM_STORE(g+26*SIMDD, MM_MUL(r4, g24));
        MM_STORE(g+27*SIMDD, MM_MUL(r5, g25));
        MM_STORE(g+30*SIMDD, MM_MUL(s4, g24));
        MM_STORE(g+31*SIMDD, MM_MUL(s5, g25));
        MM_STORE(g+4 *SIMDD, MM_FMA(r0, r0, b10));
        MM_STORE(g+5 *SIMDD, MM_FMA(r1, r1, b11));
        MM_STORE(g+8 *SIMDD, MM_FMA(s0, r0, b00));
        MM_STORE(g+9 *SIMDD, MM_FMA(s1, r1, b01));
        MM_STORE(g+16*SIMDD, MM_FMA(r2, r2, b10));
        MM_STORE(g+17*SIMDD, MM_FMA(r3, r3, b11));
        MM_STORE(g+20*SIMDD, MM_FMA(s2, r2, b00));
        MM_STORE(g+21*SIMDD, MM_FMA(s3, r3, b01));
        MM_STORE(g+28*SIMDD, MM_MUL(MM_FMA(r4, r4, b10), g24));
        MM_STORE(g+29*SIMDD, MM_MUL(MM_FMA(r5, r5, b11), g25));
        MM_STORE(g+32*SIMDD, MM_MUL(MM_FMA(s4, r4, b00), g24));
        MM_STORE(g+33*SIMDD, MM_MUL(MM_FMA(s5, r5, b01), g25));
        MM_STORE(g+10*SIMDD, r0 * (MM_LOAD(g+8 *SIMDD) + b00) + b10 * s0);
        MM_STORE(g+11*SIMDD, r1 * (MM_LOAD(g+9 *SIMDD) + b01) + b11 * s1);
        MM_STORE(g+22*SIMDD, r2 * (MM_LOAD(g+20*SIMDD) + b00) + b10 * s2);
        MM_STORE(g+23*SIMDD, r3 * (MM_LOAD(g+21*SIMDD) + b01) + b11 * s3);
        MM_STORE(g+34*SIMDD, r4 * MM_LOAD(g+32*SIMDD) + b00 * MM_LOAD(g+26*SIMDD) + b10 * MM_LOAD(g+30*SIMDD));
        MM_STORE(g+35*SIMDD, r5 * MM_LOAD(g+33*SIMDD) + b01 * MM_LOAD(g+27*SIMDD) + b11 * MM_LOAD(g+31*SIMDD));
}

static inline void _g0_2d4d_2010(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *c0x = bc->c00x;
        double *c0y = bc->c00y;
        double *c0z = bc->c00z;
        double *cpx = bc->c0px;
        double *cpy = bc->c0py;
        double *cpz = bc->c0pz;
        double *b0 = bc->b00;
        double *b1 = bc->b01;
        __MD r0 = MM_LOAD(c0x+0*SIMDD);
        __MD r1 = MM_LOAD(c0x+1*SIMDD);
        __MD r2 = MM_LOAD(c0y+0*SIMDD);
        __MD r3 = MM_LOAD(c0y+1*SIMDD);
        __MD r4 = MM_LOAD(c0z+0*SIMDD);
        __MD r5 = MM_LOAD(c0z+1*SIMDD);
        __MD s0 = MM_LOAD(cpx+0*SIMDD);
        __MD s1 = MM_LOAD(cpx+1*SIMDD);
        __MD s2 = MM_LOAD(cpy+0*SIMDD);
        __MD s3 = MM_LOAD(cpy+1*SIMDD);
        __MD s4 = MM_LOAD(cpz+0*SIMDD);
        __MD s5 = MM_LOAD(cpz+1*SIMDD);
        __MD b00 = MM_LOAD(b0+0*SIMDD);
        __MD b01 = MM_LOAD(b0+1*SIMDD);
        __MD b10 = MM_LOAD(b1+0*SIMDD);
        __MD b11 = MM_LOAD(b1+1*SIMDD);
        __MD g24 = MM_LOAD(g+24*SIMDD);
        __MD g25 = MM_LOAD(g+25*SIMDD);
        MM_STORE(g+2 *SIMDD, r0);
        MM_STORE(g+3 *SIMDD, r1);
        MM_STORE(g+6 *SIMDD, s0);
        MM_STORE(g+7 *SIMDD, s1);
        MM_STORE(g+14*SIMDD, r2);
        MM_STORE(g+15*SIMDD, r3);
        MM_STORE(g+18*SIMDD, s2);
        MM_STORE(g+19*SIMDD, s3);
        MM_STORE(g+26*SIMDD, MM_MUL(r4, g24));
        MM_STORE(g+27*SIMDD, MM_MUL(r5, g25));
        MM_STORE(g+30*SIMDD, MM_MUL(s4, g24));
        MM_STORE(g+31*SIMDD, MM_MUL(s5, g25));
        MM_STORE(g+4 *SIMDD, MM_FMA(r0, r0, b10));
        MM_STORE(g+5 *SIMDD, MM_FMA(r1, r1, b11));
        MM_STORE(g+8 *SIMDD, MM_FMA(s0, r0, b00));
        MM_STORE(g+9 *SIMDD, MM_FMA(s1, r1, b01));
        MM_STORE(g+16*SIMDD, MM_FMA(r2, r2, b10));
        MM_STORE(g+17*SIMDD, MM_FMA(r3, r3, b11));
        MM_STORE(g+20*SIMDD, MM_FMA(s2, r2, b00));
        MM_STORE(g+21*SIMDD, MM_FMA(s3, r3, b01));
        MM_STORE(g+28*SIMDD, MM_MUL(MM_FMA(r4, r4, b10), g24));
        MM_STORE(g+29*SIMDD, MM_MUL(MM_FMA(r5, r5, b11), g25));
        MM_STORE(g+32*SIMDD, MM_MUL(MM_FMA(s4, r4, b00), g24));
        MM_STORE(g+33*SIMDD, MM_MUL(MM_FMA(s5, r5, b01), g25));
        MM_STORE(g+10*SIMDD, r0 * (MM_LOAD(g+8 *SIMDD) + b00) + b10 * s0);
        MM_STORE(g+11*SIMDD, r1 * (MM_LOAD(g+9 *SIMDD) + b01) + b11 * s1);
        MM_STORE(g+22*SIMDD, r2 * (MM_LOAD(g+20*SIMDD) + b00) + b10 * s2);
        MM_STORE(g+23*SIMDD, r3 * (MM_LOAD(g+21*SIMDD) + b01) + b11 * s3);
        MM_STORE(g+34*SIMDD, r4 * MM_LOAD(g+32*SIMDD) + b00 * MM_LOAD(g+26*SIMDD) + b10 * MM_LOAD(g+30*SIMDD));
        MM_STORE(g+35*SIMDD, r5 * MM_LOAD(g+33*SIMDD) + b01 * MM_LOAD(g+27*SIMDD) + b11 * MM_LOAD(g+31*SIMDD));
}

static inline void _g0_2d4d_2100(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *cx = bc->c00x;
        double *cy = bc->c00y;
        double *cz = bc->c00z;
        double *b1 = bc->b10;
        double *rij = envs->rirj;
        ALIGNMM double rc[6*SIMDD];
        _make_rc(rc, cx, cy, cz, rij);
        __MD b10 = MM_LOAD(b1+0*SIMDD);
        __MD b11 = MM_LOAD(b1+1*SIMDD);
        __MD g32 = MM_LOAD(g+32*SIMDD);
        __MD g33 = MM_LOAD(g+33*SIMDD);
        __MD i2 = MM_SET1(2.);
        __MD r0 = MM_LOAD(rc+0*SIMDD);
        __MD r1 = MM_LOAD(rc+1*SIMDD);
        __MD r2 = MM_LOAD(rc+2*SIMDD);
        __MD r3 = MM_LOAD(rc+3*SIMDD);
        __MD r4 = MM_LOAD(rc+4*SIMDD);
        __MD r5 = MM_LOAD(rc+5*SIMDD);
        __MD s0 = MM_LOAD(cx+0*SIMDD);
        __MD s1 = MM_LOAD(cx+1*SIMDD);
        __MD s2 = MM_LOAD(cy+0*SIMDD);
        __MD s3 = MM_LOAD(cy+1*SIMDD);
        __MD s4 = MM_LOAD(cz+0*SIMDD);
        __MD s5 = MM_LOAD(cz+1*SIMDD);
        MM_STORE(g+2 *SIMDD, s0);
        MM_STORE(g+3 *SIMDD, s1);
        MM_STORE(g+8 *SIMDD, r0);
        MM_STORE(g+9 *SIMDD, r1);
        MM_STORE(g+18*SIMDD, s2);
        MM_STORE(g+19*SIMDD, s3);
        MM_STORE(g+24*SIMDD, r2);
        MM_STORE(g+25*SIMDD, r3);
        MM_STORE(g+34*SIMDD, s4 * g32);
        MM_STORE(g+35*SIMDD, s5 * g33);
        MM_STORE(g+40*SIMDD, r4 * g32);
        MM_STORE(g+41*SIMDD, r5 * g33);
        MM_STORE(g+4 *SIMDD, s0 * s0 + b10);
        MM_STORE(g+5 *SIMDD, s1 * s1 + b11);
        MM_STORE(g+10*SIMDD, s0 * r0 + b10);
        MM_STORE(g+11*SIMDD, s1 * r1 + b11);
        MM_STORE(g+20*SIMDD, s2 * s2 + b10);
        MM_STORE(g+21*SIMDD, s3 * s3 + b11);
        MM_STORE(g+26*SIMDD, s2 * r2 + b10);
        MM_STORE(g+27*SIMDD, s3 * r3 + b11);
        MM_STORE(g+36*SIMDD,(s4 * s4 + b10) * g32);
        MM_STORE(g+37*SIMDD,(s5 * s5 + b11) * g33);
        MM_STORE(g+42*SIMDD,(s4 * r4 + b10) * g32);
        MM_STORE(g+43*SIMDD,(s5 * r5 + b11) * g33);
        MM_STORE(g+12*SIMDD, r0 * MM_LOAD(g+4 *SIMDD ) + i2 * b10 * s0);
        MM_STORE(g+13*SIMDD, r1 * MM_LOAD(g+5 *SIMDD ) + i2 * b11 * s1);
        MM_STORE(g+28*SIMDD, r2 * MM_LOAD(g+20*SIMDD ) + i2 * b10 * s2);
        MM_STORE(g+29*SIMDD, r3 * MM_LOAD(g+21*SIMDD ) + i2 * b11 * s3);
        MM_STORE(g+44*SIMDD, r4 * MM_LOAD(g+36*SIMDD ) + i2 * b10 * MM_LOAD(g+34*SIMDD ));
        MM_STORE(g+45*SIMDD, r5 * MM_LOAD(g+37*SIMDD ) + i2 * b11 * MM_LOAD(g+35*SIMDD ));
}

static inline void _g0_2d4d_3000(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        double *cx = bc->c00x;
        double *cy = bc->c00y;
        double *cz = bc->c00z;
        double *b  = bc->b10;
        __MD b0 = MM_LOAD(b+0*SIMDD);
        __MD b1 = MM_LOAD(b+1*SIMDD);
        __MD i3 = MM_SET1(3.);
        __MD r0 = MM_LOAD(cx+0*SIMDD);
        __MD r1 = MM_LOAD(cx+1*SIMDD);
        __MD r2 = MM_LOAD(cy+0*SIMDD);
        __MD r3 = MM_LOAD(cy+1*SIMDD);
        __MD r4 = MM_LOAD(cz+0*SIMDD);
        __MD r5 = MM_LOAD(cz+1*SIMDD);
        __MD g16 = MM_LOAD(g+16*SIMDD);
        __MD g17 = MM_LOAD(g+17*SIMDD);
        MM_STORE(g+2 *SIMDD, r0);
        MM_STORE(g+3 *SIMDD, r1);
        MM_STORE(g+10*SIMDD, r2);
        MM_STORE(g+11*SIMDD, r3);
        MM_STORE(g+18*SIMDD, r4 * g16);
        MM_STORE(g+19*SIMDD, r5 * g17);
        MM_STORE(g+4 *SIMDD, r0 * r0 + b0);
        MM_STORE(g+5 *SIMDD, r1 * r1 + b1);
        MM_STORE(g+12*SIMDD, r2 * r2 + b0);
        MM_STORE(g+13*SIMDD, r3 * r3 + b1);
        MM_STORE(g+20*SIMDD,(r4 * r4 + b0)* g16);
        MM_STORE(g+21*SIMDD,(r5 * r5 + b1)* g17);
        MM_STORE(g+6 *SIMDD, r0 *(r0 * r0 + i3 * b0));
        MM_STORE(g+7 *SIMDD, r1 *(r1 * r1 + i3 * b1));
        MM_STORE(g+14*SIMDD, r2 *(r2 * r2 + i3 * b0));
        MM_STORE(g+15*SIMDD, r3 *(r3 * r3 + i3 * b1));
        MM_STORE(g+22*SIMDD,(r4 * r4 + i3 * b0) * MM_LOAD(g+18*SIMDD));
        MM_STORE(g+23*SIMDD,(r5 * r5 + i3 * b1) * MM_LOAD(g+19*SIMDD));
}
/************** end special g0_4d results *************/



void CINTg0_2e_2d4d_unrolled(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        int type_ijkl = ((envs->li_ceil << 6) | (envs->lj_ceil << 4) |
                         (envs->lk_ceil << 2) | (envs->ll_ceil));
        switch (type_ijkl) {
                case 0b00000000: _g0_2d4d_0000(g, bc, envs); return;
                case 0b00000001: _g0_2d4d_0001(g, bc, envs); return;
                case 0b00000010: _g0_2d4d_0002(g, bc, envs); return;
                case 0b00000011: _g0_2d4d_0003(g, bc, envs); return;
                case 0b00000100: _g0_2d4d_0010(g, bc, envs); return;
                case 0b00000101: _g0_2d4d_0011(g, bc, envs); return;
                case 0b00000110: _g0_2d4d_0012(g, bc, envs); return;
                case 0b00001000: _g0_2d4d_0020(g, bc, envs); return;
                case 0b00001001: _g0_2d4d_0021(g, bc, envs); return;
                case 0b00001100: _g0_2d4d_0030(g, bc, envs); return;
                case 0b00010000: _g0_2d4d_0100(g, bc, envs); return;
                case 0b00010001: _g0_2d4d_0101(g, bc, envs); return;
                case 0b00010010: _g0_2d4d_0102(g, bc, envs); return;
                case 0b00010100: _g0_2d4d_0110(g, bc, envs); return;
                case 0b00010101: _g0_2d4d_0111(g, bc, envs); return;
                case 0b00011000: _g0_2d4d_0120(g, bc, envs); return;
                case 0b00100000: _g0_2d4d_0200(g, bc, envs); return;
                case 0b00100001: _g0_2d4d_0201(g, bc, envs); return;
                case 0b00100100: _g0_2d4d_0210(g, bc, envs); return;
                case 0b00110000: _g0_2d4d_0300(g, bc, envs); return;
                case 0b01000000: _g0_2d4d_1000(g, bc, envs); return;
                case 0b01000001: _g0_2d4d_1001(g, bc, envs); return;
                case 0b01000010: _g0_2d4d_1002(g, bc, envs); return;
                case 0b01000100: _g0_2d4d_1010(g, bc, envs); return;
                case 0b01000101: _g0_2d4d_1011(g, bc, envs); return;
                case 0b01001000: _g0_2d4d_1020(g, bc, envs); return;
                case 0b01010000: _g0_2d4d_1100(g, bc, envs); return;
                case 0b01010001: _g0_2d4d_1101(g, bc, envs); return;
                case 0b01010100: _g0_2d4d_1110(g, bc, envs); return;
                case 0b01100000: _g0_2d4d_1200(g, bc, envs); return;
                case 0b10000000: _g0_2d4d_2000(g, bc, envs); return;
                case 0b10000001: _g0_2d4d_2001(g, bc, envs); return;
                case 0b10000100: _g0_2d4d_2010(g, bc, envs); return;
                case 0b10010000: _g0_2d4d_2100(g, bc, envs); return;
                case 0b11000000: _g0_2d4d_3000(g, bc, envs); return;
        }
        fprintf(stderr, "Dimension error for CINTg0_2e_lj2d4d: iklj = %d %d %d %d",
               (int)envs->li_ceil, (int)envs->lk_ceil,
               (int)envs->ll_ceil, (int)envs->lj_ceil);
}

void CINTg0_2e_lj2d4d(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        CINTg0_2e_2d(g, bc, envs);
        CINTg0_lj_4d(g, envs);
}
void CINTg0_2e_kj2d4d(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        CINTg0_2e_2d(g, bc, envs);
        CINTg0_kj_4d(g, envs);
}
void CINTg0_2e_ik2d4d(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        CINTg0_2e_2d(g, bc, envs);
        CINTg0_ik_4d(g, envs);
}
void CINTg0_2e_il2d4d(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        CINTg0_2e_2d(g, bc, envs);
        CINTg0_il_4d(g, envs);
}

/*
 * g[i,k,l,j] = < ik | lj > = ( i j | k l )
 */
int CINTg0_2e(double *g, double *cutoff, Rys2eT *bc, CINTEnvVars *envs, int count)
{
        ALIGNMM double aij[SIMDD];
        ALIGNMM double akl[SIMDD];
        ALIGNMM double a0[SIMDD];
        ALIGNMM double a1[SIMDD];
        ALIGNMM double aijkl[SIMDD];
        ALIGNMM double fac1[SIMDD];
        ALIGNMM double x[SIMDD];
        ALIGNMM double rijrkl[SIMDD*3];
        ALIGNMM double rijrx[SIMDD*3];
        ALIGNMM double rklrx[SIMDD*3];
        double *rij = envs->rij;
        double *rkl = envs->rkl;
        double *u = bc->u;
        double *w = bc->w;
        __MD ra, r0, r1, r2, r3, r4, r5, r6, r7, r8;
        int nroots = envs->nrys_roots;
        int i;

        //:for (int k = 0; k < count; k++) {
        //:        aij[k] = envs->ai[k] + envs->aj[k];
        //:        akl[k] = envs->ak[k] + envs->al[k];
        //:        aijkl[k] = aij[k] + akl[k];
        //:        a1[k] = aij[k] * akl[k];
        //:        a0[k] = a1[k] / aijkl[k];
        //:        //fac1[k] = sqrt(a0[k] / (a1[k] * a1[k] * a1[k])) * envs->fac[k];
        //:        fac1[k] = envs->fac[k] / (sqrt(aijakl[k]) * a1[k]);
        //:}
        r2 = MM_ADD(MM_LOAD(envs->ai),  MM_LOAD(envs->aj));
        r3 = MM_ADD(MM_LOAD(envs->ak),  MM_LOAD(envs->al));
        MM_STORE(aij, r2);
        MM_STORE(akl, r3);
        r1 = MM_MUL(r2, r3);
        MM_STORE(a1, r1);
        ra = MM_ADD(r2, r3);
        MM_STORE(aijkl, ra);
        r0 = MM_DIV(r1, ra);
        MM_STORE(a0, r0);

#ifdef WITH_RANGE_COULOMB
// Not recommended to mix range-separated Coulomb with regular Coulomb operator.
// Keep this for backward compatibility to cint2
        const double omega = envs->env[PTR_RANGE_OMEGA];
        ALIGNMM double theta[SIMDD];
        if (omega != 0) {
                //:theta = omega * omega / (omega * omega + a0);
                r0 = MM_SET1(omega);
                r0 = MM_MUL(r0, r0);
                r1 = MM_LOAD(a0);
                r0 = MM_DIV(r0, MM_ADD(r0, r1));
                MM_STORE(theta, r0);
                if (omega > 0) { // long-range part of range-separated Coulomb operator
                        //:a0 *= theta;
                        MM_STORE(a0, MM_MUL(r0, r1));
                }
        }

        r1 = MM_LOAD(a1);
        r0 = MM_DIV(MM_LOAD(a0), MM_MUL(r1, MM_MUL(r1, r1)));
        MM_STORE(fac1, MM_MUL(MM_SQRT(r0), MM_LOAD(envs->fac)));
#else
        MM_STORE(fac1, MM_DIV(MM_LOAD(envs->fac), MM_MUL(MM_SQRT(ra), r1)));
#endif

        r0 = MM_LOAD(rij+0*SIMDD);
        r1 = MM_LOAD(rij+1*SIMDD);
        r2 = MM_LOAD(rij+2*SIMDD);
        r3 = MM_LOAD(rkl+0*SIMDD);
        r4 = MM_LOAD(rkl+1*SIMDD);
        r5 = MM_LOAD(rkl+2*SIMDD);

        //:for (k = 0; k < count; k++) {
        //:        rijrkl[0*SIMDD+k] = rij[0*SIMDD+k] - rkl[0*SIMDD+k];
        //:        rijrkl[1*SIMDD+k] = rij[1*SIMDD+k] - rkl[1*SIMDD+k];
        //:        rijrkl[2*SIMDD+k] = rij[2*SIMDD+k] - rkl[2*SIMDD+k];
        //:        rijrx[0*SIMDD+k] = rij[0*SIMDD+k] - envs->rx_in_rijrx[0];
        //:        rijrx[1*SIMDD+k] = rij[1*SIMDD+k] - envs->rx_in_rijrx[1];
        //:        rijrx[2*SIMDD+k] = rij[2*SIMDD+k] - envs->rx_in_rijrx[2];
        //:        rklrx[0*SIMDD+k] = rkl[0*SIMDD+k] - envs->rx_in_rklrx[0];
        //:        rklrx[1*SIMDD+k] = rkl[1*SIMDD+k] - envs->rx_in_rklrx[1];
        //:        rklrx[2*SIMDD+k] = rkl[2*SIMDD+k] - envs->rx_in_rklrx[2];
        //:        x[k] = a0[k] *(rijrkl[0*SIMDD+k] * rijrkl[0*SIMDD+k]
        //:                     + rijrkl[1*SIMDD+k] * rijrkl[1*SIMDD+k]
        //:                     + rijrkl[2*SIMDD+k] * rijrkl[2*SIMDD+k]);
        //:}
        r6 = MM_SUB(r0, r3); MM_STORE(rijrkl+0*SIMDD, r6);
        r7 = MM_SUB(r1, r4); MM_STORE(rijrkl+1*SIMDD, r7);
        r8 = MM_SUB(r2, r5); MM_STORE(rijrkl+2*SIMDD, r8);
        ra = MM_FMA(r6, r6, MM_FMA(r7, r7, MM_MUL(r8, r8)));
        MM_STORE(x, MM_MUL(MM_LOAD(a0), ra));
        MM_STORE(rijrx+0*SIMDD, MM_SUB(r0, MM_SET1(envs->rx_in_rijrx[0])));
        MM_STORE(rijrx+1*SIMDD, MM_SUB(r1, MM_SET1(envs->rx_in_rijrx[1])));
        MM_STORE(rijrx+2*SIMDD, MM_SUB(r2, MM_SET1(envs->rx_in_rijrx[2])));
        MM_STORE(rklrx+0*SIMDD, MM_SUB(r3, MM_SET1(envs->rx_in_rklrx[0])));
        MM_STORE(rklrx+1*SIMDD, MM_SUB(r4, MM_SET1(envs->rx_in_rklrx[1])));
        MM_STORE(rklrx+2*SIMDD, MM_SUB(r5, MM_SET1(envs->rx_in_rklrx[2])));

//ABORT        ALIGNMM double erfx[SIMDD];
//ABORT        if (envs->g_size == 1) {
//ABORT                // gz = erf(sqrt(rr*aij*akl/(aij+akl)))/sqrt(rr)
//ABORT                MM_STORE(erfx, MM_SQRT(MM_LOAD(x)));
//ABORT                for (k = 0; k < count; k++) {
//ABORT                        erfx[k] = erf(erfx[k]);
//ABORT                }
//ABORT                MM_STORE(gz, MM_DIV(MM_LOAD(erfx), MM_SQRT(MM_LOAD(x))));
//ABORT                return;
//ABORT        }
#ifdef WITH_RANGE_COULOMB
        if (omega < 0) {
                int all_negligible = _CINTsr_rys_roots_batch(
                        envs, x, theta, u, w, cutoff, count);
                if (all_negligible) {
                        // g still has to be evaluated since iempty (which
                        // indicates whether g is zero) in cint2e is determined
                        // before calling the g0_2e function.
                        for (i = 0; i < envs->g_size * SIMDD * 3; i++) {
                                g[i] = 0;
                        }
                        return 0;
                }
        } else {
                _CINTrys_roots_batch(nroots, x, u, w, count);
        }
#else
        _CINTrys_roots_batch(nroots, x, u, w, count);
#endif

        double *gx = g;
        double *gy = gx + envs->g_size * SIMDD;
        double *gz = gy + envs->g_size * SIMDD;
        //:for (i = 0; i < nroots; i++) {
        //:for (k = 0; k < count; k++) {
        //:        gx[i*SIMDD+k] = 1;
        //:        gy[i*SIMDD+k] = 1;
        //:        gz[i*SIMDD+k] = w[k+i*SIMDD] * fac1[k];
        //:} }
        r0 = MM_LOAD(fac1);
        r1 = MM_SET1(1.);
        for (i = 0; i < nroots; i++) {
                //MM_STORE(gx+i*SIMDD, r1);
                //MM_STORE(gy+i*SIMDD, r1);
                MM_STORE(gz+i*SIMDD, MM_MUL(MM_LOAD(w+i*SIMDD), r0));
        }
        if (envs->g_size == 1) { // ssss
                return 1;
        }

#ifdef WITH_RANGE_COULOMB
        if (omega > 0) {
                /* u[:] = tau^2 / (1 - tau^2)
                 * transform u[:] to theta^-1 tau^2 / (theta^-1 - tau^2)
                 * so the rest code can be reused.
                 */
                //:for (irys = 0; irys < nroots; irys++) {
                //:        u[irys] /= u[irys] + 1 - u[irys] * theta;
                //:}
                r0 = MM_LOAD(theta);
                r1 = MM_SET1(1.);
                for (i = 0; i < nroots; i++) {
                        r2 = MM_LOAD(u+i*SIMDD);
                        r3 = r2 + r1 - r2 * r0;
                        MM_STORE(u+i*SIMDD, MM_DIV(r2, r3));
                }
        }
#endif

        double *b00 = bc->b00;
        double *b10 = bc->b10;
        double *b01 = bc->b01;
        double *c00x = bc->c00x;
        double *c00y = bc->c00y;
        double *c00z = bc->c00z;
        double *c0px = bc->c0px;
        double *c0py = bc->c0py;
        double *c0pz = bc->c0pz;

        ALIGNMM double tmp1[MXRYSROOTS*SIMDD];
        ALIGNMM double tmp4[MXRYSROOTS*SIMDD];
        //:double u2[SIMDD];
        //:double tmp2, tmp3;
        //:double div[SIMDD];
        //:for (i = 0; i < nroots; i++) {
        //:for (k = 0; k < count; k++) {
        //:        u2[k] = a0[k] * u[k+i*SIMDD];
        //:        div[k] = 1 / (u2[k] * aijkl[k] + a1[k]);
        //:        tmp1[i*SIMDD+k] = u2[k] * div[k];
        //:        tmp4[i*SIMDD+k] = .5 * div[k];
        //:
        //:        b00[i*SIMDD+k] = 0.5 * tmp1[i*SIMDD+k];
        //:        tmp2 = tmp1[i*SIMDD+k] * akl[k];
        //:        tmp3 = tmp1[i*SIMDD+k] * aij[k];
        //:        b10[i*SIMDD+k] = b00[i*SIMDD+k] + tmp4[i*SIMDD+k] * akl[k];
        //:        b01[i*SIMDD+k] = b00[i*SIMDD+k] + tmp4[i*SIMDD+k] * aij[k];
        //:        c00x[i*SIMDD+k] = rijrx[0*SIMDD+k] - tmp2 * rijrkl[0*SIMDD+k];
        //:        c00y[i*SIMDD+k] = rijrx[1*SIMDD+k] - tmp2 * rijrkl[1*SIMDD+k];
        //:        c00z[i*SIMDD+k] = rijrx[2*SIMDD+k] - tmp2 * rijrkl[2*SIMDD+k];
        //:        c0px[i*SIMDD+k] = rklrx[0*SIMDD+k] + tmp3 * rijrkl[0*SIMDD+k];
        //:        c0py[i*SIMDD+k] = rklrx[1*SIMDD+k] + tmp3 * rijrkl[1*SIMDD+k];
        //:        c0pz[i*SIMDD+k] = rklrx[2*SIMDD+k] + tmp3 * rijrkl[2*SIMDD+k];
        //:} }

        ra = MM_LOAD(aijkl);
        r0 = MM_LOAD(a0);
        r1 = MM_LOAD(a1);
        r2 = MM_SET1(.5);
        r3 = MM_SET1(1.);
        for (i = 0; i < nroots; i++) {
                r4 = MM_MUL(r0, MM_LOAD(u+i*SIMDD));
                r5 = MM_DIV(r3, MM_FMA(r4, ra, r1));
                MM_STORE(tmp4+i*SIMDD, MM_MUL(r2, r5));
                MM_STORE(tmp1+i*SIMDD, MM_MUL(r4, r5));
        }
        ra = MM_SET1(.5);
        r2 = MM_LOAD(akl);
        r3 = MM_LOAD(aij);
        for (i = 0; i < nroots; i++) {
                r0 = MM_MUL(ra, MM_LOAD(tmp1+i*SIMDD));
                MM_STORE(b00+i*SIMDD, r0);
                r1 = MM_LOAD(tmp4+i*SIMDD);
                MM_STORE(b10+i*SIMDD, MM_FMA(r1, r2, r0));
                MM_STORE(b01+i*SIMDD, MM_FMA(r1, r3, r0));
        }
        r4 = MM_LOAD(rijrkl+0*SIMDD);
        r5 = MM_LOAD(rijrkl+1*SIMDD);
        r6 = MM_LOAD(rijrkl+2*SIMDD);
        ra = MM_LOAD(akl);
        r1 = MM_LOAD(rijrx+0*SIMDD);
        r2 = MM_LOAD(rijrx+1*SIMDD);
        r3 = MM_LOAD(rijrx+2*SIMDD);
        for (i = 0; i < nroots; i++) {
                r0 = MM_MUL(MM_LOAD(tmp1+i*SIMDD), ra);
                MM_STORE(c00x+i*SIMDD, MM_FNMA(r0, r4, r1));
                MM_STORE(c00y+i*SIMDD, MM_FNMA(r0, r5, r2));
                MM_STORE(c00z+i*SIMDD, MM_FNMA(r0, r6, r3));
        }
        ra = MM_LOAD(aij);
        r1 = MM_LOAD(rklrx+0*SIMDD);
        r2 = MM_LOAD(rklrx+1*SIMDD);
        r3 = MM_LOAD(rklrx+2*SIMDD);
        for (i = 0; i < nroots; i++) {
                r0 = MM_MUL(MM_LOAD(tmp1+i*SIMDD), ra);
                MM_STORE(c0px+i*SIMDD, MM_FMA (r0, r4, r1));
                MM_STORE(c0py+i*SIMDD, MM_FMA (r0, r5, r2));
                MM_STORE(c0pz+i*SIMDD, MM_FMA (r0, r6, r3));
        }

        (*envs->f_g0_2d4d)(g, bc, envs);
        return 1;
}


/*
 * ( \nabla i j | kl )
 */
void CINTnabla1i_2e(double *f, double *g,
                    int li, int lj, int lk, int ll, CINTEnvVars *envs)
{
        int nroots = envs->nrys_roots;
        int di = envs->g_stride_i;
        DEF_GXYZ(double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);
        __MD ai2 = MM_MUL(MM_SET1(-2.), MM_LOAD(envs->ai));

        if (nroots == 1) { // nabla_i ssss
                assert(li == 0);
                assert(lj == 0);
                assert(lk == 0);
                assert(ll == 0);
                MM_STORE(fx, MM_MUL(ai2, MM_LOAD(gx+di*SIMDD)));
                MM_STORE(fy, MM_MUL(ai2, MM_LOAD(gy+di*SIMDD)));
                MM_STORE(fz, MM_MUL(ai2, MM_LOAD(gz+di*SIMDD)));
        } else {
                int i, j, k, l, n, ptr;
                int dk = envs->g_stride_k;
                int dl = envs->g_stride_l;
                int dj = envs->g_stride_j;
                double *RESTRICT p1x = gx - di * SIMDD;
                double *RESTRICT p1y = gy - di * SIMDD;
                double *RESTRICT p1z = gz - di * SIMDD;
                double *RESTRICT p2x = gx + di * SIMDD;
                double *RESTRICT p2y = gy + di * SIMDD;
                double *RESTRICT p2z = gz + di * SIMDD;
                __MD ri;

                for (j = 0; j <= lj; j++)
                for (l = 0; l <= ll; l++)
                for (k = 0; k <= lk; k++) {
                        ptr = dj * j + dl * l + dk * k;
                        //f(...,0,...) = -2*ai*g(...,1,...)
                        for (n = ptr; n < ptr+nroots; n++) {
//fx[n*SIMDD] = ai2 * p2x[n*SIMDD];
//fy[n*SIMDD] = ai2 * p2y[n*SIMDD];
//fz[n*SIMDD] = ai2 * p2z[n*SIMDD];
MM_STORE(fx+n*SIMDD, MM_MUL(ai2, MM_LOAD(p2x+n*SIMDD)));
MM_STORE(fy+n*SIMDD, MM_MUL(ai2, MM_LOAD(p2y+n*SIMDD)));
MM_STORE(fz+n*SIMDD, MM_MUL(ai2, MM_LOAD(p2z+n*SIMDD)));
                        }
                        ptr += di;
                        //f(...,i,...) = i*g(...,i-1,...)-2*ai*g(...,i+1,...)
                        ri = MM_SET1(1.);
                        for (i = 1; i <= li; i++) {
                                for (n = ptr; n < ptr+nroots; n++) {
//fx[n*SIMDD] = i*p1x[n*SIMDD] + ai2*p2x[n*SIMDD];
//fy[n*SIMDD] = i*p1y[n*SIMDD] + ai2*p2y[n*SIMDD];
//fz[n*SIMDD] = i*p1z[n*SIMDD] + ai2*p2z[n*SIMDD];
MM_STORE(fx+n*SIMDD, ri * MM_LOAD(p1x+n*SIMDD) + ai2 * MM_LOAD(p2x+n*SIMDD));
MM_STORE(fy+n*SIMDD, ri * MM_LOAD(p1y+n*SIMDD) + ai2 * MM_LOAD(p2y+n*SIMDD));
MM_STORE(fz+n*SIMDD, ri * MM_LOAD(p1z+n*SIMDD) + ai2 * MM_LOAD(p2z+n*SIMDD));
                                }
                                ptr += di;
                                ri = MM_ADD(ri, MM_SET1(1.));
                        }
                }
        }
}


/*
 * ( i \nabla j | kl )
 */
void CINTnabla1j_2e(double *f, double *g,
                    int li, int lj, int lk, int ll, CINTEnvVars *envs)
{
        int nroots = envs->nrys_roots;
        int dj = envs->g_stride_j;
        DEF_GXYZ(double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);
        __MD aj2 = MM_MUL(MM_SET1(-2.), MM_LOAD(envs->aj));

        if (nroots == 1) { // nabla_j ssss
                assert(li == 0);
                assert(lj == 0);
                assert(lk == 0);
                assert(ll == 0);
                MM_STORE(fx, MM_MUL(aj2, MM_LOAD(gx+dj*SIMDD)));
                MM_STORE(fy, MM_MUL(aj2, MM_LOAD(gy+dj*SIMDD)));
                MM_STORE(fz, MM_MUL(aj2, MM_LOAD(gz+dj*SIMDD)));
        } else {
                int i, j, k, l, n, ptr;
                int di = envs->g_stride_i;
                int dk = envs->g_stride_k;
                int dl = envs->g_stride_l;
                double *RESTRICT p1x = gx - dj * SIMDD;
                double *RESTRICT p1y = gy - dj * SIMDD;
                double *RESTRICT p1z = gz - dj * SIMDD;
                double *RESTRICT p2x = gx + dj * SIMDD;
                double *RESTRICT p2y = gy + dj * SIMDD;
                double *RESTRICT p2z = gz + dj * SIMDD;
                __MD rj;

                //f(...,0,...) = -2*aj*g(...,1,...)
                for (l = 0; l <= ll; l++) {
                for (k = 0; k <= lk; k++) {
                        ptr = dl * l + dk * k;
                        for (i = 0; i <= li; i++) {
                                for (n = ptr; n < ptr+nroots; n++) {
//fx[n*SIMDD] = aj2 * p2x[n*SIMDD];
//fy[n*SIMDD] = aj2 * p2y[n*SIMDD];
//fz[n*SIMDD] = aj2 * p2z[n*SIMDD];
MM_STORE(fx+n*SIMDD, MM_MUL(aj2, MM_LOAD(p2x+n*SIMDD)));
MM_STORE(fy+n*SIMDD, MM_MUL(aj2, MM_LOAD(p2y+n*SIMDD)));
MM_STORE(fz+n*SIMDD, MM_MUL(aj2, MM_LOAD(p2z+n*SIMDD)));
                                }
                                ptr += di;
                        }
                } }
                //f(...,j,...) = j*g(...,j-1,...)-2*aj*g(...,j+1,...)
                for (j = 1; j <= lj; j++) {
                        rj = MM_SET1(j);
                        for (l = 0; l <= ll; l++) {
                        for (k = 0; k <= lk; k++) {
                                ptr = dj * j + dl * l + dk * k;
                                for (i = 0; i <= li; i++) {
                                        for (n = ptr; n < ptr+nroots; n++) {
//fx[n*SIMDD] = j*p1x[n*SIMDD] + aj2*p2x[n*SIMDD];
//fy[n*SIMDD] = j*p1y[n*SIMDD] + aj2*p2y[n*SIMDD];
//fz[n*SIMDD] = j*p1z[n*SIMDD] + aj2*p2z[n*SIMDD];
MM_STORE(fx+n*SIMDD, rj * MM_LOAD(p1x+n*SIMDD) + aj2 * MM_LOAD(p2x+n*SIMDD));
MM_STORE(fy+n*SIMDD, rj * MM_LOAD(p1y+n*SIMDD) + aj2 * MM_LOAD(p2y+n*SIMDD));
MM_STORE(fz+n*SIMDD, rj * MM_LOAD(p1z+n*SIMDD) + aj2 * MM_LOAD(p2z+n*SIMDD));
                                        }
                                        ptr += di;
                                }
                        } }
                }
        }
}


/*
 * ( ij | \nabla k l )
 */
void CINTnabla1k_2e(double *f, double *g,
                    int li, int lj, int lk, int ll, CINTEnvVars *envs)
{
        int nroots = envs->nrys_roots;
        int dk = envs->g_stride_k;
        DEF_GXYZ(double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);
        __MD ak2 = MM_MUL(MM_SET1(-2.), MM_LOAD(envs->ak));

        if (nroots == 1) { // nabla_k ssss
                assert(li == 0);
                assert(lj == 0);
                assert(lk == 0);
                assert(ll == 0);
                MM_STORE(fx, MM_MUL(ak2, MM_LOAD(gx+dk*SIMDD)));
                MM_STORE(fy, MM_MUL(ak2, MM_LOAD(gy+dk*SIMDD)));
                MM_STORE(fz, MM_MUL(ak2, MM_LOAD(gz+dk*SIMDD)));
        } else {
                int i, j, k, l, n, ptr;
                int di = envs->g_stride_i;
                int dl = envs->g_stride_l;
                int dj = envs->g_stride_j;
                double *RESTRICT p1x = gx - dk * SIMDD;
                double *RESTRICT p1y = gy - dk * SIMDD;
                double *RESTRICT p1z = gz - dk * SIMDD;
                double *RESTRICT p2x = gx + dk * SIMDD;
                double *RESTRICT p2y = gy + dk * SIMDD;
                double *RESTRICT p2z = gz + dk * SIMDD;
                __MD rk;

                for (j = 0; j <= lj; j++)
                for (l = 0; l <= ll; l++) {
                        ptr = dj * j + dl * l;
                        //f(...,0,...) = -2*ak*g(...,1,...)
                        for (i = 0; i <= li; i++) {
                                for (n = ptr; n < ptr+nroots; n++) {
//fx[n*SIMDD] = ak2 * p2x[n*SIMDD];
//fy[n*SIMDD] = ak2 * p2y[n*SIMDD];
//fz[n*SIMDD] = ak2 * p2z[n*SIMDD];
MM_STORE(fx+n*SIMDD, MM_MUL(ak2, MM_LOAD(p2x+n*SIMDD)));
MM_STORE(fy+n*SIMDD, MM_MUL(ak2, MM_LOAD(p2y+n*SIMDD)));
MM_STORE(fz+n*SIMDD, MM_MUL(ak2, MM_LOAD(p2z+n*SIMDD)));
                                }
                                ptr += di;
                        }
                        //f(...,k,...) = k*g(...,k-1,...)-2*ak*g(...,k+1,...)
                        rk = MM_SET1(1.);
                        for (k = 1; k <= lk; k++) {
                                ptr = dj * j + dl * l + dk * k;
                                for (i = 0; i <= li; i++) {
                                        for (n = ptr; n < ptr+nroots; n++) {
//fx[n*SIMDD] = k*p1x[n*SIMDD] + ak2*p2x[n*SIMDD];
//fy[n*SIMDD] = k*p1y[n*SIMDD] + ak2*p2y[n*SIMDD];
//fz[n*SIMDD] = k*p1z[n*SIMDD] + ak2*p2z[n*SIMDD];
MM_STORE(fx+n*SIMDD, rk * MM_LOAD(p1x+n*SIMDD) + ak2 * MM_LOAD(p2x+n*SIMDD));
MM_STORE(fy+n*SIMDD, rk * MM_LOAD(p1y+n*SIMDD) + ak2 * MM_LOAD(p2y+n*SIMDD));
MM_STORE(fz+n*SIMDD, rk * MM_LOAD(p1z+n*SIMDD) + ak2 * MM_LOAD(p2z+n*SIMDD));
                                        }
                                        ptr += di;
                                }
                                rk = MM_ADD(rk, MM_SET1(1.));
                        }
                }
        }
}


/*
 * ( ij | k \nabla l )
 */
void CINTnabla1l_2e(double *f, double *g,
                    int li, int lj, int lk, int ll, CINTEnvVars *envs)
{
        int nroots = envs->nrys_roots;
        int dl = envs->g_stride_l;
        DEF_GXYZ(double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);
        __MD al2 = MM_MUL(MM_SET1(-2.), MM_LOAD(envs->al));

        if (nroots == 1) { // nabla_l ssss
                assert(li == 0);
                assert(lj == 0);
                assert(lk == 0);
                assert(ll == 0);
                MM_STORE(fx, MM_MUL(al2, MM_LOAD(gx+dl*SIMDD)));
                MM_STORE(fy, MM_MUL(al2, MM_LOAD(gy+dl*SIMDD)));
                MM_STORE(fz, MM_MUL(al2, MM_LOAD(gz+dl*SIMDD)));
        } else {
                int i, j, k, l, n, ptr;
                int di = envs->g_stride_i;
                int dk = envs->g_stride_k;
                int dj = envs->g_stride_j;
                double *RESTRICT p1x = gx - dl * SIMDD;
                double *RESTRICT p1y = gy - dl * SIMDD;
                double *RESTRICT p1z = gz - dl * SIMDD;
                double *RESTRICT p2x = gx + dl * SIMDD;
                double *RESTRICT p2y = gy + dl * SIMDD;
                double *RESTRICT p2z = gz + dl * SIMDD;
                __MD rl;

                for (j = 0; j <= lj; j++) {
                        //f(...,0,...) = -2*al*g(...,1,...)
                        for (k = 0; k <= lk; k++) {
                                ptr = dj * j + dk * k;
                                for (i = 0; i <= li; i++) {
                                        for (n = ptr; n < ptr+nroots; n++) {
//fx[n*SIMDD] = al2 * p2x[n*SIMDD];
//fy[n*SIMDD] = al2 * p2y[n*SIMDD];
//fz[n*SIMDD] = al2 * p2z[n*SIMDD];
MM_STORE(fx+n*SIMDD, MM_MUL(al2, MM_LOAD(p2x+n*SIMDD)));
MM_STORE(fy+n*SIMDD, MM_MUL(al2, MM_LOAD(p2y+n*SIMDD)));
MM_STORE(fz+n*SIMDD, MM_MUL(al2, MM_LOAD(p2z+n*SIMDD)));
                                        }
                                        ptr += di;
                                }
                        }
                        //f(...,l,...) = l*g(...,l-1,...)-2*al*g(...,l+1,...)
                        rl = MM_SET1(1.);
                        for (l = 1; l <= ll; l++) {
                                for (k = 0; k <= lk; k++) {
                                        ptr = dj * j + dl * l + dk * k;
                                        for (i = 0; i <= li; i++, ptr += di) {
                                        for (n = ptr; n < ptr+nroots; n++) {
//fx[n*SIMDD] = l*p1x[n*SIMDD] + al2*p2x[n*SIMDD];
//fy[n*SIMDD] = l*p1y[n*SIMDD] + al2*p2y[n*SIMDD];
//fz[n*SIMDD] = l*p1z[n*SIMDD] + al2*p2z[n*SIMDD];
MM_STORE(fx+n*SIMDD, rl * MM_LOAD(p1x+n*SIMDD) + al2 * MM_LOAD(p2x+n*SIMDD));
MM_STORE(fy+n*SIMDD, rl * MM_LOAD(p1y+n*SIMDD) + al2 * MM_LOAD(p2y+n*SIMDD));
MM_STORE(fz+n*SIMDD, rl * MM_LOAD(p1z+n*SIMDD) + al2 * MM_LOAD(p2z+n*SIMDD));
                                        } }
                                }
                                rl = MM_ADD(rl, MM_SET1(1.));
                        }
                }
        }
}

/*
 * ( x^1 i j | kl )
 * ri is the shift from the center R_O to the center of |i>
 * r - R_O = (r-R_i) + ri, ri = R_i - R_O
 */
void CINTx1i_2e(double *f, double *g, double *ri,
                int li, int lj, int lk, int ll, CINTEnvVars *envs)
{
        int nroots = envs->nrys_roots;
        int di = envs->g_stride_i;
        DEF_GXYZ(double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);
        __MD r0 = MM_SET1(ri[0]);
        __MD r1 = MM_SET1(ri[1]);
        __MD r2 = MM_SET1(ri[2]);

        if (nroots == 1) { // x_i ssss
                assert(li == 0);
                assert(lj == 0);
                assert(lk == 0);
                assert(ll == 0);
                MM_STORE(fx, MM_FMA(r0, MM_LOAD(gx), MM_LOAD(gx+di*SIMDD)));
                MM_STORE(fy, MM_FMA(r1, MM_LOAD(gy), MM_LOAD(gy+di*SIMDD)));
                MM_STORE(fz, MM_FMA(r2, MM_LOAD(gz), MM_LOAD(gz+di*SIMDD)));
        } else {
                int i, j, k, l, n, ptr;
                int dk = envs->g_stride_k;
                int dl = envs->g_stride_l;
                int dj = envs->g_stride_j;
                double *RESTRICT p1x = gx + di * SIMDD;
                double *RESTRICT p1y = gy + di * SIMDD;
                double *RESTRICT p1z = gz + di * SIMDD;

                for (j = 0; j <= lj; j++)
                for (l = 0; l <= ll; l++) {
                for (k = 0; k <= lk; k++) {
                        //f(...,0:li,...) = g(...,1:li+1,...) + ri(1)*g(...,0:li,...)
                        ptr = dj * j + dl * l + dk * k;
                        for (i = 0; i <= li; i++) {
                                for (n = ptr; n < ptr+nroots; n++) {
//fx[n*SIMDD] = p1x[n*SIMDD] + ri[0] * gx[n*SIMDD];
//fy[n*SIMDD] = p1y[n*SIMDD] + ri[1] * gy[n*SIMDD];
//fz[n*SIMDD] = p1z[n*SIMDD] + ri[2] * gz[n*SIMDD];
MM_STORE(fx+n*SIMDD, MM_FMA(r0, MM_LOAD(gx+n*SIMDD), MM_LOAD(p1x+n*SIMDD)));
MM_STORE(fy+n*SIMDD, MM_FMA(r1, MM_LOAD(gy+n*SIMDD), MM_LOAD(p1y+n*SIMDD)));
MM_STORE(fz+n*SIMDD, MM_FMA(r2, MM_LOAD(gz+n*SIMDD), MM_LOAD(p1z+n*SIMDD)));
                                }
                                ptr += di;
                        }
                } }
        }
}


/*
 * ( i x^1 j | kl )
 */
void CINTx1j_2e(double *f, double *g, double *rj,
                int li, int lj, int lk, int ll, CINTEnvVars *envs)
{
        int nroots = envs->nrys_roots;
        int dj = envs->g_stride_j;
        DEF_GXYZ(double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);
        __MD r0 = MM_SET1(rj[0]);
        __MD r1 = MM_SET1(rj[1]);
        __MD r2 = MM_SET1(rj[2]);

        if (nroots == 1) { // x_j ssss
                assert(li == 0);
                assert(lj == 0);
                assert(lk == 0);
                assert(ll == 0);
                MM_STORE(fx, MM_FMA(r0, MM_LOAD(gx), MM_LOAD(gx+dj*SIMDD)));
                MM_STORE(fy, MM_FMA(r1, MM_LOAD(gy), MM_LOAD(gy+dj*SIMDD)));
                MM_STORE(fz, MM_FMA(r2, MM_LOAD(gz), MM_LOAD(gz+dj*SIMDD)));
        } else {
                int i, j, k, l, n, ptr;
                int di = envs->g_stride_i;
                int dk = envs->g_stride_k;
                int dl = envs->g_stride_l;
                double *RESTRICT p1x = gx + dj * SIMDD;
                double *RESTRICT p1y = gy + dj * SIMDD;
                double *RESTRICT p1z = gz + dj * SIMDD;

                for (j = 0; j <= lj; j++)
                for (l = 0; l <= ll; l++) {
                for (k = 0; k <= lk; k++) {
                        // f(...,0:lj,...) = g(...,1:lj+1,...) + rj(1)*g(...,0:lj,...)
                        ptr = dj * j + dl * l + dk * k;
                        for (i = 0; i <= li; i++) {
                                for (n = ptr; n < ptr+nroots; n++) {
//fx[n*SIMDD] = p1x[n*SIMDD] + rj[0] * gx[n*SIMDD];
//fy[n*SIMDD] = p1y[n*SIMDD] + rj[1] * gy[n*SIMDD];
//fz[n*SIMDD] = p1z[n*SIMDD] + rj[2] * gz[n*SIMDD];
MM_STORE(fx+n*SIMDD, MM_FMA(r0, MM_LOAD(gx+n*SIMDD), MM_LOAD(p1x+n*SIMDD)));
MM_STORE(fy+n*SIMDD, MM_FMA(r1, MM_LOAD(gy+n*SIMDD), MM_LOAD(p1y+n*SIMDD)));
MM_STORE(fz+n*SIMDD, MM_FMA(r2, MM_LOAD(gz+n*SIMDD), MM_LOAD(p1z+n*SIMDD)));
                                }
                                ptr += di;
                        }
                } }
        }
}


/*
 * ( ij | x^1 k l )
 */
void CINTx1k_2e(double *f, double *g, double *rk,
                int li, int lj, int lk, int ll, CINTEnvVars *envs)
{
        int nroots = envs->nrys_roots;
        int dk = envs->g_stride_k;
        DEF_GXYZ(double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);
        __MD r0 = MM_SET1(rk[0]);
        __MD r1 = MM_SET1(rk[1]);
        __MD r2 = MM_SET1(rk[2]);

        if (nroots == 1) { // x_k ssss
                assert(li == 0);
                assert(lj == 0);
                assert(lk == 0);
                assert(ll == 0);
                MM_STORE(fx, MM_FMA(r0, MM_LOAD(gx), MM_LOAD(gx+dk*SIMDD)));
                MM_STORE(fy, MM_FMA(r1, MM_LOAD(gy), MM_LOAD(gy+dk*SIMDD)));
                MM_STORE(fz, MM_FMA(r2, MM_LOAD(gz), MM_LOAD(gz+dk*SIMDD)));
        } else {
                int i, j, k, l, n, ptr;
                int di = envs->g_stride_i;
                int dl = envs->g_stride_l;
                int dj = envs->g_stride_j;
                double *RESTRICT p1x = gx + dk * SIMDD;
                double *RESTRICT p1y = gy + dk * SIMDD;
                double *RESTRICT p1z = gz + dk * SIMDD;

                for (j = 0; j <= lj; j++)
                for (l = 0; l <= ll; l++) {
                for (k = 0; k <= lk; k++) {
                        // f(...,0:lk,...) = g(...,1:lk+1,...) + rk(1)*g(...,0:lk,...)
                        ptr = dj * j + dl * l + dk * k;
                        for (i = 0; i <= li; i++) {
                                for (n = ptr; n < ptr+nroots; n++) {
//fx[n*SIMDD] = p1x[n*SIMDD] + rk[0] * gx[n*SIMDD];
//fy[n*SIMDD] = p1y[n*SIMDD] + rk[1] * gy[n*SIMDD];
//fz[n*SIMDD] = p1z[n*SIMDD] + rk[2] * gz[n*SIMDD];
MM_STORE(fx+n*SIMDD, MM_FMA(r0, MM_LOAD(gx+n*SIMDD), MM_LOAD(p1x+n*SIMDD)));
MM_STORE(fy+n*SIMDD, MM_FMA(r1, MM_LOAD(gy+n*SIMDD), MM_LOAD(p1y+n*SIMDD)));
MM_STORE(fz+n*SIMDD, MM_FMA(r2, MM_LOAD(gz+n*SIMDD), MM_LOAD(p1z+n*SIMDD)));
                                }
                                ptr += di;
                        }
                } }
        }
}


/*
 * ( i j | x^1 kl )
 */
void CINTx1l_2e(double *f, double *g, double *rl,
                int li, int lj, int lk, int ll, CINTEnvVars *envs)
{
        int nroots = envs->nrys_roots;
        int dl = envs->g_stride_l;
        DEF_GXYZ(double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);
        __MD r0 = MM_SET1(rl[0]);
        __MD r1 = MM_SET1(rl[1]);
        __MD r2 = MM_SET1(rl[2]);

        if (nroots == 1) { // x_l ssss
                assert(li == 0);
                assert(lj == 0);
                assert(lk == 0);
                assert(ll == 0);
                MM_STORE(fx, MM_FMA(r0, MM_LOAD(gx), MM_LOAD(gx+dl*SIMDD)));
                MM_STORE(fy, MM_FMA(r1, MM_LOAD(gy), MM_LOAD(gy+dl*SIMDD)));
                MM_STORE(fz, MM_FMA(r2, MM_LOAD(gz), MM_LOAD(gz+dl*SIMDD)));
        } else {
                int i, j, k, l, n, ptr;
                int di = envs->g_stride_i;
                int dk = envs->g_stride_k;
                int dj = envs->g_stride_j;
                double *RESTRICT p1x = gx + dl * SIMDD;
                double *RESTRICT p1y = gy + dl * SIMDD;
                double *RESTRICT p1z = gz + dl * SIMDD;

                for (j = 0; j <= lj; j++)
                for (l = 0; l <= ll; l++) {
                for (k = 0; k <= lk; k++) {
                        // f(...,0:ll,...) = g(...,1:ll+1,...) + rl(1)*g(...,0:ll,...)
                        ptr = dj * j + dl * l + dk * k;
                        for (i = 0; i <= li; i++) {
                                for (n = ptr; n < ptr+nroots; n++) {
//fx[n*SIMDD] = p1x[n*SIMDD] + rl[0] * gx[n*SIMDD];
//fy[n*SIMDD] = p1y[n*SIMDD] + rl[1] * gy[n*SIMDD];
//fz[n*SIMDD] = p1z[n*SIMDD] + rl[2] * gz[n*SIMDD];
MM_STORE(fx+n*SIMDD, MM_FMA(r0, MM_LOAD(gx+n*SIMDD), MM_LOAD(p1x+n*SIMDD)));
MM_STORE(fy+n*SIMDD, MM_FMA(r1, MM_LOAD(gy+n*SIMDD), MM_LOAD(p1y+n*SIMDD)));
MM_STORE(fz+n*SIMDD, MM_FMA(r2, MM_LOAD(gz+n*SIMDD), MM_LOAD(p1z+n*SIMDD)));
                                }
                                ptr += di;
                        }
                } }
        }
}


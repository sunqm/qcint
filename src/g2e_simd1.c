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
#include "simd.h"
#include "rys_roots.h"
#include "misc.h"
#include "g2e.h"

#define DEF_GXYZ(type, G, GX, GY, GZ) \
        type *GX = G; \
        type *GY = G + envs->g_size; \
        type *GZ = G + envs->g_size * 2

/*
 * g(nroots,0:nmax,0:mmax)
 */
void CINTg0_2e_2d_simd1(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        int nroots = envs->nrys_roots;
        int i, j, m, n;
        DEF_GXYZ(double, g, gx, gy, gz);
        int nmax = envs->li_ceil + envs->lj_ceil;
        int mmax = envs->lk_ceil + envs->ll_ceil;
        int dm = envs->g2d_klmax;
        int dn = envs->g2d_ijmax;
        double *c00x = bc->c00x;
        double *c00y = bc->c00y;
        double *c00z = bc->c00z;
        double *c0px = bc->c0px;
        double *c0py = bc->c0py;
        double *c0pz = bc->c0pz;
        double *b01 = bc->b01;
        double *b00 = bc->b00;
        double *b10 = bc->b10;
        double *p0x, *p0y, *p0z;
        double *p1x, *p1y, *p1z;
        double *p2x, *p2y, *p2z;
        double xm, xn;
        if (nmax > 0) {
                p0x = gx + dn;
                p0y = gy + dn;
                p0z = gz + dn;
                p1x = gx - dn;
                p1y = gy - dn;
                p1z = gz - dn;
                // gx(irys,0,1) = c00(irys) * gx(irys,0,0)
                for (i = 0; i < nroots; i++) {
                        p0x[i] = c00x[i] * gx[i];
                        p0y[i] = c00y[i] * gy[i];
                        p0z[i] = c00z[i] * gz[i];
                }
                // gx(irys,0,n+1) = c00(irys)*gx(irys,0,n)
                // + n*b10(irys)*gx(irys,0,n-1)
                for (n = 1, xn = 1; n < nmax; n++, xn+=1) {
                for (i = 0, j = n*dn; i < nroots; i++, j++) {
                        p0x[j] = c00x[i] * gx[j] + xn * b10[i] * p1x[j];
                        p0y[j] = c00y[i] * gy[j] + xn * b10[i] * p1y[j];
                        p0z[j] = c00z[i] * gz[j] + xn * b10[i] * p1z[j];
                } }
        }

        if (mmax > 0) {
                p0x = gx + dm;
                p0y = gy + dm;
                p0z = gz + dm;
                p1x = gx - dm;
                p1y = gy - dm;
                p1z = gz - dm;
                // gx(irys,1,0) = c0p(irys) * gx(irys,0,0)
                for (i = 0; i < nroots; i++) {
                        p0x[i] = c0px[i] * gx[i];
                        p0y[i] = c0py[i] * gy[i];
                        p0z[i] = c0pz[i] * gz[i];
                }
                // gx(irys,m+1,0) = c0p(irys)*gx(irys,m,0)
                // + m*b01(irys)*gx(irys,m-1,0)
                for (m = 1, xm = 1; m < mmax; m++, xm+=1) {
                for (i = 0, j = m*dm; i < nroots; i++, j++) {
                        p0x[j] = c0px[i] * gx[j] + xm * b01[i] * p1x[j];
                        p0y[j] = c0py[i] * gy[j] + xm * b01[i] * p1y[j];
                        p0z[j] = c0pz[i] * gz[j] + xm * b01[i] * p1z[j];
                } }
        }

        if (nmax > 0 && mmax > 0) {
                p0x = gx + dn;
                p0y = gy + dn;
                p0z = gz + dn;
                p1x = p0x + dm;
                p1y = p0y + dm;
                p1z = p0z + dm;
                // gx(irys,1,1) = c0p(irys)*gx(irys,0,1)
                // + b00(irys)*gx(irys,0,0)
                for (i = 0; i < nroots; i++) {
                        p1x[i] = c0px[i] * p0x[i] +  b00[i] * gx[i];
                        p1y[i] = c0py[i] * p0y[i] +  b00[i] * gy[i];
                        p1z[i] = c0pz[i] * p0z[i] +  b00[i] * gz[i];
                }

                p0x = gx + dm;
                p0y = gy + dm;
                p0z = gz + dm;
                p1x = gx - dn;
                p1y = gy - dn;
                p1z = gz - dn;
                p2x = gx - dm;
                p2y = gy - dm;
                p2z = gz - dm;
                // gx(irys,m+1,1) = c0p(irys)*gx(irys,m,1)
                // + m*b01(irys)*gx(irys,m-1,1)
                // + b00(irys)*gx(irys,m,0)
                for (m = 1, xm = 1; m < mmax; m++, xm+=1) {
                for (i = 0, j = m*dm+dn; i < nroots; i++, j++) {
                         p0x[j] = c0px[i] * gx[j] + xm * b01[i] * p2x[j] + b00[i] * p1x[j];
                         p0y[j] = c0py[i] * gy[j] + xm * b01[i] * p2y[j] + b00[i] * p1y[j];
                         p0z[j] = c0pz[i] * gz[j] + xm * b01[i] * p2z[j] + b00[i] * p1z[j];
                } }

                p0x = gx + dn;
                p0y = gy + dn;
                p0z = gz + dn;
                // gx(irys,m,n+1) = c00(irys)*gx(irys,m,n)
                // + n*b10(irys)*gx(irys,m,n-1)
                // + m*b00(irys)*gx(irys,m-1,n)
                for (m = 1, xm = 1; m <= mmax; m++, xm+=1) {
                for (n = 1, xn = 1; n < nmax; n++, xn+=1) {
                for (i = 0, j = m*dm+n*dn; i < nroots; i++, j++) {
                         p0x[j] = c00x[i] * gx[j] + xn * b10[i] * p1x[j] + xm * b00[i] * p2x[j];
                         p0y[j] = c00y[i] * gy[j] + xn * b10[i] * p1y[j] + xm * b00[i] * p2y[j];
                         p0z[j] = c00z[i] * gz[j] + xn * b10[i] * p1z[j] + xm * b00[i] * p2z[j];
                } } }
        }
}


/*
 * g0[i,k,l,j] = < ik | lj > = ( i j | k l )
 */
/* 2d is based on l,j */
void CINTg0_lj_4d_simd1(double *g, CINTEnvVars *envs)
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
        double *p1x, *p1y, *p1z, *p2x, *p2y, *p2z;

        // g(i,...,j) = rirj * g(i-1,...,j) +  g(i-1,...,j+1)
        p1x = gx - di;
        p1y = gy - di;
        p1z = gz - di;
        p2x = gx - di + dj;
        p2y = gy - di + dj;
        p2z = gz - di + dj;
        for (i = 1; i <= li; i++) {
        for (j = 0; j <= nmax-i; j++) {
        for (l = 0; l <= mmax; l++) {
                ptr = j*dj + l*dl + i*di;
                for (n = ptr; n < ptr+nroots; n++) {
                        gx[n] = rirj[0] * p1x[n] + p2x[n];
                        gy[n] = rirj[1] * p1y[n] + p2y[n];
                        gz[n] = rirj[2] * p1z[n] + p2z[n];
                }
        } } }

        // g(...,k,l,..) = rkrl * g(...,k-1,l,..) + g(...,k-1,l+1,..)
        p1x = gx - dk;
        p1y = gy - dk;
        p1z = gz - dk;
        p2x = gx - dk + dl;
        p2y = gy - dk + dl;
        p2z = gz - dk + dl;
        for (j = 0; j <= lj; j++) {
        for (k = 1; k <= lk; k++) {
        for (l = 0; l <= mmax-k; l++) {
                ptr = j*dj + l*dl + k*dk;
                for (n = ptr; n < ptr+dk; n++) {
                        gx[n] = rkrl[0] * p1x[n] + p2x[n];
                        gy[n] = rkrl[1] * p1y[n] + p2y[n];
                        gz[n] = rkrl[2] * p1z[n] + p2z[n];
                }
        } } }
}
/* 2d is based on k,j */
void CINTg0_kj_4d_simd1(double *g, CINTEnvVars *envs)
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
        double *p1x, *p1y, *p1z, *p2x, *p2y, *p2z;

        // g(i,...,j) = rirj * g(i-1,...,j) +  g(i-1,...,j+1)
        p1x = gx - di;
        p1y = gy - di;
        p1z = gz - di;
        p2x = gx - di + dj;
        p2y = gy - di + dj;
        p2z = gz - di + dj;
        for (i = 1; i <= li; i++) {
        for (j = 0; j <= nmax-i; j++) {
        for (k = 0; k <= mmax; k++) {
                ptr = j*dj + k*dk + i*di;
                for (n = ptr; n < ptr+nroots; n++) {
                        gx[n] = rirj[0] * p1x[n] + p2x[n];
                        gy[n] = rirj[1] * p1y[n] + p2y[n];
                        gz[n] = rirj[2] * p1z[n] + p2z[n];
                }
        } } }

        // g(...,k,l,..) = rkrl * g(...,k,l-1,..) + g(...,k+1,l-1,..)
        p1x = gx - dl;
        p1y = gy - dl;
        p1z = gz - dl;
        p2x = gx - dl + dk;
        p2y = gy - dl + dk;
        p2z = gz - dl + dk;
        for (j = 0; j <= lj; j++) {
        for (l = 1; l <= ll; l++) {
        for (k = 0; k <= mmax-l; k++) {
                ptr = j*dj + l*dl + k*dk;
                for (n = ptr; n < ptr+dk; n++) {
                        gx[n] = rkrl[0] * p1x[n] + p2x[n];
                        gy[n] = rkrl[1] * p1y[n] + p2y[n];
                        gz[n] = rkrl[2] * p1z[n] + p2z[n];
                }
        } } }
}
/* 2d is based on i,l */
void CINTg0_il_4d_simd1(double *g, CINTEnvVars *envs)
{
        int lk = envs->lk_ceil;
        int lj = envs->lj_ceil;
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
        double *p1x, *p1y, *p1z, *p2x, *p2y, *p2z;

        // g(...,k,l,..) = rkrl * g(...,k-1,l,..) + g(...,k-1,l+1,..)
        p1x = gx - dk;
        p1y = gy - dk;
        p1z = gz - dk;
        p2x = gx - dk + dl;
        p2y = gy - dk + dl;
        p2z = gz - dk + dl;
        for (k = 1; k <= lk; k++) {
        for (l = 0; l <= mmax-k; l++) {
        for (i = 0; i <= nmax; i++) {
                ptr = l*dl + k*dk + i*di;
                for (n = ptr; n < ptr+nroots; n++) {
                        gx[n] = rkrl[0] * p1x[n] + p2x[n];
                        gy[n] = rkrl[1] * p1y[n] + p2y[n];
                        gz[n] = rkrl[2] * p1z[n] + p2z[n];
                }
        } } }

        // g(i,...,j) = rirj * g(i,...,j-1) +  g(i+1,...,j-1)
        p1x = gx - dj;
        p1y = gy - dj;
        p1z = gz - dj;
        p2x = gx - dj + di;
        p2y = gy - dj + di;
        p2z = gz - dj + di;
        for (j = 1; j <= lj; j++) {
        for (l = 0; l <= ll; l++) {
        for (k = 0; k <= lk; k++) {
                ptr = j*dj + l*dl + k*dk;
                for (n = ptr; n < ptr+dk-di*j; n++) {
                        gx[n] = rirj[0] * p1x[n] + p2x[n];
                        gy[n] = rirj[1] * p1y[n] + p2y[n];
                        gz[n] = rirj[2] * p1z[n] + p2z[n];
                }
        } } }
}
/* 2d is based on i,k */
void CINTg0_ik_4d_simd1(double *g, CINTEnvVars *envs)
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
        double *p1x, *p1y, *p1z, *p2x, *p2y, *p2z;

        // g(...,k,l,..) = rkrl * g(...,k,l-1,..) + g(...,k+1,l-1,..)
        p1x = gx - dl;
        p1y = gy - dl;
        p1z = gz - dl;
        p2x = gx - dl + dk;
        p2y = gy - dl + dk;
        p2z = gz - dl + dk;
        for (l = 1; l <= ll; l++) {
                // (:,i) is full, so loop:k and loop:n can be merged to
                // for(n = l*dl; n < ptr+dl-dk*l; n++)
                for (k = 0; k <= mmax-l; k++) {
                for (i = 0; i <= nmax; i++) {
                        ptr = l*dl + k*dk + i*di;
                        for (n = ptr; n < ptr+nroots; n++) {
                                gx[n] = rkrl[0] * p1x[n] + p2x[n];
                                gy[n] = rkrl[1] * p1y[n] + p2y[n];
                                gz[n] = rkrl[2] * p1z[n] + p2z[n];
                        }
                } }
        }

        // g(i,...,j) = rirj * g(i,...,j-1) +  g(i+1,...,j-1)
        p1x = gx - dj;
        p1y = gy - dj;
        p1z = gz - dj;
        p2x = gx - dj + di;
        p2y = gy - dj + di;
        p2z = gz - dj + di;
        for (j = 1; j <= lj; j++) {
        for (l = 0; l <= ll; l++) {
        for (k = 0; k <= lk; k++) {
                ptr = j*dj + l*dl + k*dk;
                for (n = ptr; n < ptr+dk-di*j; n++) {
                        gx[n] = rirj[0] * p1x[n] + p2x[n];
                        gy[n] = rirj[1] * p1y[n] + p2y[n];
                        gz[n] = rirj[2] * p1z[n] + p2z[n];
                }
        } } }
}
/************* some special g0_4d results *************/
/* 4 digits stand for i_ceil, k_ceil, l_ceil, j_ceil */
static inline void _make_rc(double *rc, double *cx, double *cy, double *cz, double *r)
{
        rc[0] = r[0] + cx[0];
        rc[1] = r[0] + cx[1];
        rc[2] = r[1] + cy[0];
        rc[3] = r[1] + cy[1];
        rc[4] = r[2] + cz[0];
        rc[5] = r[2] + cz[1];
}
static inline void _g0_lj_4d_0001(double *g, double *cx, double *cy, double *cz,
                                  double *r)
{
        g[1] = cx[0];
        g[3] = cy[0];
        g[5] = cz[0] * g[4];
}
static inline void _g0_lj_4d_1000(double *g, double *cx, double *cy, double *cz,
                                  double *r)
{
        g[1] = r[0] + cx[0];
        g[5] = r[1] + cy[0];
        g[9] =(r[2] + cz[0]) * g[8];
}
static inline void _g0_lj_4d_0002(double *g, double *cx, double *cy, double *cz,
                                  double *b, double *r)
{
        g[2 ] = cx[0];
        g[3 ] = cx[1];
        g[8 ] = cy[0];
        g[9 ] = cy[1];
        g[14] = cz[0] * g[12];
        g[15] = cz[1] * g[13];
        g[4 ] = cx[0] * cx[0] + b[0];
        g[5 ] = cx[1] * cx[1] + b[1];
        g[10] = cy[0] * cy[0] + b[0];
        g[11] = cy[1] * cy[1] + b[1];
        g[16] =(cz[0] * cz[0] + b[0])* g[12];
        g[17] =(cz[1] * cz[1] + b[1])* g[13];
}
static inline void _g0_lj_4d_1001(double *g, double *cx, double *cy, double *cz,
                                  double *b, double *r)
{
        double rc[6];
        _make_rc(rc, cx, cy, cz, r);
        g[4 ] = cx[0];
        g[5 ] = cx[1];
        g[16] = cy[0];
        g[17] = cy[1];
        g[28] = cz[0] * g[24];
        g[29] = cz[1] * g[25];
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[14] = rc[2];
        g[15] = rc[3];
        g[26] = rc[4] * g[24];
        g[27] = rc[5] * g[25];
        g[6 ] = rc[0] * cx[0] + b[0];
        g[7 ] = rc[1] * cx[1] + b[1];
        g[18] = rc[2] * cy[0] + b[0];
        g[19] = rc[3] * cy[1] + b[1];
        g[30] =(rc[4] * cz[0] + b[0])* g[24];
        g[31] =(rc[5] * cz[1] + b[1])* g[25];
}
static inline void _g0_lj_4d_2000(double *g, double *cx, double *cy, double *cz,
                                  double *b, double *r)
{
        double rc[6];
        _make_rc(rc, cx, cy, cz, r);
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[20] = rc[2];
        g[21] = rc[3];
        g[38] = rc[4] * g[36];
        g[39] = rc[5] * g[37];
        g[4 ] = rc[0] * rc[0] + b[0];
        g[5 ] = rc[1] * rc[1] + b[1];
        g[22] = rc[2] * rc[2] + b[0];
        g[23] = rc[3] * rc[3] + b[1];
        g[40] =(rc[4] * rc[4] + b[0])* g[36];
        g[41] =(rc[5] * rc[5] + b[1])* g[37];
}
static inline void _g0_lj_4d_0003(double *g, double *cx, double *cy, double *cz,
                                  double *b, double *r)
{
        g[2 ] = cx[0];
        g[3 ] = cx[1];
        g[10] = cy[0];
        g[11] = cy[1];
        g[18] = cz[0] * g[16];
        g[19] = cz[1] * g[17];
        g[4 ] = cx[0] * cx[0] + b[0];
        g[5 ] = cx[1] * cx[1] + b[1];
        g[12] = cy[0] * cy[0] + b[0];
        g[13] = cy[1] * cy[1] + b[1];
        g[20] =(cz[0] * cz[0] + b[0])* g[16];
        g[21] =(cz[1] * cz[1] + b[1])* g[17];
        g[6 ] = cx[0] *(cx[0] * cx[0] + 3 * b[0]);
        g[7 ] = cx[1] *(cx[1] * cx[1] + 3 * b[1]);
        g[14] = cy[0] *(cy[0] * cy[0] + 3 * b[0]);
        g[15] = cy[1] *(cy[1] * cy[1] + 3 * b[1]);
        g[22] =(cz[0] * cz[0] + 3 *b[0])* g[18];
        g[23] =(cz[1] * cz[1] + 3 *b[1])* g[19];
}
static inline void _g0_lj_4d_1002(double *g, double *cx, double *cy, double *cz,
                                  double *b, double *r)
{
        double rc[6];
        _make_rc(rc, cx, cy, cz, r);
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[4 ] = cx[0];
        g[5 ] = cx[1];
        g[18] = rc[2];
        g[19] = rc[3];
        g[20] = cy[0];
        g[21] = cy[1];
        g[34] = rc[4] * g[32];
        g[35] = rc[5] * g[33];
        g[36] = cz[0] * g[32];
        g[37] = cz[1] * g[33];
        g[6 ] = rc[0] * cx[0] + b[0];
        g[7 ] = rc[1] * cx[1] + b[1];
        g[8 ] = cx[0] * cx[0] + b[0];
        g[9 ] = cx[1] * cx[1] + b[1];
        g[22] = rc[2] * cy[0] + b[0];
        g[23] = rc[3] * cy[1] + b[1];
        g[24] = cy[0] * cy[0] + b[0];
        g[25] = cy[1] * cy[1] + b[1];
        g[38] =(rc[4] * cz[0] + b[0])* g[32];
        g[39] =(rc[5] * cz[1] + b[1])* g[33];
        g[40] =(cz[0] * cz[0] + b[0])* g[32];
        g[41] =(cz[1] * cz[1] + b[1])* g[33];
        g[10] = rc[0] * g[8 ] + 2* b[0] * cx[0];
        g[11] = rc[1] * g[9 ] + 2* b[1] * cx[1];
        g[26] = rc[2] * g[24] + 2* b[0] * cy[0];
        g[27] = rc[3] * g[25] + 2* b[1] * cy[1];
        g[42] = rc[4] * g[40] + 2* b[0] * g[36];
        g[43] = rc[5] * g[41] + 2* b[1] * g[37];
}
static inline void _g0_lj_4d_2001(double *g, double *cx, double *cy, double *cz,
                                  double *b, double *r)
{
        double rc[6];
        _make_rc(rc, cx, cy, cz, r);
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[6 ] = cx[0];
        g[7 ] = cx[1];
        g[26] = rc[2];
        g[27] = rc[3];
        g[30] = cy[0];
        g[31] = cy[1];
        g[50] = rc[4] * g[48];
        g[51] = rc[5] * g[49];
        g[54] = cz[0] * g[48];
        g[55] = cz[1] * g[49];
        g[4 ] = rc[0] * rc[0] + b[0];
        g[5 ] = rc[1] * rc[1] + b[1];
        g[8 ] = cx[0] * rc[0] + b[0];
        g[9 ] = cx[1] * rc[1] + b[1];
        g[28] = rc[2] * rc[2] + b[0];
        g[29] = rc[3] * rc[3] + b[1];
        g[32] = cy[0] * rc[2] + b[0];
        g[33] = cy[1] * rc[3] + b[1];
        g[52] =(rc[4] * rc[4] + b[0])* g[48];
        g[53] =(rc[5] * rc[5] + b[1])* g[49];
        g[56] =(cz[0] * rc[4] + b[0])* g[48];
        g[57] =(cz[1] * rc[5] + b[1])* g[49];
        g[10] = cx[0] * g[4 ] + 2 *b[0] * rc[0];
        g[11] = cx[1] * g[5 ] + 2 *b[1] * rc[1];
        g[34] = cy[0] * g[28] + 2 *b[0] * rc[2];
        g[35] = cy[1] * g[29] + 2 *b[1] * rc[3];
        g[58] = cz[0] * g[52] + 2 *b[0] * g[50];
        g[59] = cz[1] * g[53] + 2 *b[1] * g[51];
}
static inline void _g0_lj_4d_3000(double *g, double *cx, double *cy, double *cz,
                                  double *b, double *r)
{
        double rc[6];
        _make_rc(rc, cx, cy, cz, r);
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[34] = rc[2];
        g[35] = rc[3];
        g[66] = rc[4] * g[64];
        g[67] = rc[5] * g[65];
        g[4 ] = rc[0] * rc[0] + b[0];
        g[5 ] = rc[1] * rc[1] + b[1];
        g[36] = rc[2] * rc[2] + b[0];
        g[37] = rc[3] * rc[3] + b[1];
        g[68] =(rc[4] * rc[4] + b[0])* g[64];
        g[69] =(rc[5] * rc[5] + b[1])* g[65];
        g[6 ] = rc[0] *(rc[0] * rc[0] + 3 * b[0]);
        g[7 ] = rc[1] *(rc[1] * rc[1] + 3 * b[1]);
        g[38] = rc[2] *(rc[2] * rc[2] + 3 * b[0]);
        g[39] = rc[3] *(rc[3] * rc[3] + 3 * b[1]);
        g[70] =(rc[4] * rc[4] + 3 *b[0]) * g[66];
        g[71] =(rc[5] * rc[5] + 3 *b[1]) * g[67];
}
static inline void _g0_lj_4d_0011(double *g, double *c0x, double *c0y, double *c0z,
                                  double *cpx, double *cpy, double *cpz, double *b,
                                  double *r0, double *rp)
{
        g[2 ] = cpx[0];
        g[3 ] = cpx[1];
        g[4 ] = c0x[0];
        g[5 ] = c0x[1];
        g[10] = cpy[0];
        g[11] = cpy[1];
        g[12] = c0y[0];
        g[13] = c0y[1];
        g[18] = cpz[0] * g[16];
        g[19] = cpz[1] * g[17];
        g[20] = c0z[0] * g[16];
        g[21] = c0z[1] * g[17];
        g[6 ] = cpx[0] * c0x[0] + b[0];
        g[7 ] = cpx[1] * c0x[1] + b[1];
        g[14] = cpy[0] * c0y[0] + b[0];
        g[15] = cpy[1] * c0y[1] + b[1];
        g[22] =(cpz[0] * c0z[0] + b[0]) * g[16];
        g[23] =(cpz[1] * c0z[1] + b[1]) * g[17];
}
static inline void _g0_lj_4d_1010(double *g, double *c0x, double *c0y, double *c0z,
                                  double *cpx, double *cpy, double *cpz, double *b,
                                  double *r0, double *rp)
{
        double rc[6];
        _make_rc(rc, c0x, c0y, c0z, r0);
        g[4 ] = cpx[0];
        g[5 ] = cpx[1];
        g[20] = cpy[0];
        g[21] = cpy[1];
        g[36] = cpz[0] * g[32];
        g[37] = cpz[1] * g[33];
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[18] = rc[2];
        g[19] = rc[3];
        g[34] = rc[4] * g[32];
        g[35] = rc[5] * g[33];
        g[6 ] = rc[0] * cpx[0] + b[0];
        g[7 ] = rc[1] * cpx[1] + b[1];
        g[22] = rc[2] * cpy[0] + b[0];
        g[23] = rc[3] * cpy[1] + b[1];
        g[38] =(rc[4] * cpz[0] + b[0]) * g[32];
        g[39] =(rc[5] * cpz[1] + b[1]) * g[33];
}
static inline void _g0_lj_4d_0101(double *g, double *c0x, double *c0y, double *c0z,
                                  double *cpx, double *cpy, double *cpz, double *b,
                                  double *r0, double *rp)
{
        double rc[6];
        _make_rc(rc, cpx, cpy, cpz, rp);
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[18] = rc[2];
        g[19] = rc[3];
        g[34] = rc[4] * g[32];
        g[35] = rc[5] * g[33];
        g[8 ] = c0x[0];
        g[9 ] = c0x[1];
        g[24] = c0y[0];
        g[25] = c0y[1];
        g[40] = c0z[0] * g[32];
        g[41] = c0z[1] * g[33];
        g[10] = rc[0] * c0x[0] + b[0];
        g[11] = rc[1] * c0x[1] + b[1];
        g[26] = rc[2] * c0y[0] + b[0];
        g[27] = rc[3] * c0y[1] + b[1];
        g[42] =(rc[4] * c0z[0] + b[0]) * g[32];
        g[43] =(rc[5] * c0z[1] + b[1]) * g[33];
}
static inline void _g0_lj_4d_1100(double *g, double *c0x, double *c0y, double *c0z,
                                  double *cpx, double *cpy, double *cpz, double *b,
                                  double *r0, double *rp)
{
        double rc0[6];
        double rcp[6];
        _make_rc(rc0, c0x, c0y, c0z, r0);
        _make_rc(rcp, cpx, cpy, cpz, rp);
        g[2 ] = rc0[0];
        g[3 ] = rc0[1];
        g[4 ] = rcp[0];
        g[5 ] = rcp[1];
        g[34] = rc0[2];
        g[35] = rc0[3];
        g[36] = rcp[2];
        g[37] = rcp[3];
        g[66] = rc0[4] * g[64];
        g[67] = rc0[5] * g[65];
        g[68] = rcp[4] * g[64];
        g[69] = rcp[5] * g[65];
        g[6 ] = rc0[0] * rcp[0] + b[0];
        g[7 ] = rc0[1] * rcp[1] + b[1];
        g[38] = rc0[2] * rcp[2] + b[0];
        g[39] = rc0[3] * rcp[3] + b[1];
        g[70] =(rc0[4] * rcp[4] + b[0]) * g[64];
        g[71] =(rc0[5] * rcp[5] + b[1]) * g[65];
}
static inline void _g0_lj_4d_0021(double *g, double *c0x, double *c0y, double *c0z,
                                  double *cpx, double *cpy, double *cpz, 
                                  double *b0, double *b1)
{
        g[2 ] = cpx[0];
        g[3 ] = cpx[1];
        g[6 ] = c0x[0];
        g[7 ] = c0x[1];
        g[14] = cpy[0];
        g[15] = cpy[1];
        g[18] = c0y[0];
        g[19] = c0y[1];
        g[26] = cpz[0] * g[24];
        g[27] = cpz[1] * g[25];
        g[30] = c0z[0] * g[24];
        g[31] = c0z[1] * g[25];
        g[4 ] = cpx[0] * cpx[0] + b1[0];
        g[5 ] = cpx[1] * cpx[1] + b1[1];
        g[8 ] = cpx[0] * c0x[0] + b0[0];
        g[9 ] = cpx[1] * c0x[1] + b0[1];
        g[16] = cpy[0] * cpy[0] + b1[0];
        g[17] = cpy[1] * cpy[1] + b1[1];
        g[20] = cpy[0] * c0y[0] + b0[0];
        g[21] = cpy[1] * c0y[1] + b0[1];
        g[28] =(cpz[0] * cpz[0] + b1[0]) * g[24];
        g[29] =(cpz[1] * cpz[1] + b1[1]) * g[25];
        g[32] =(cpz[0] * c0z[0] + b0[0]) * g[24];
        g[33] =(cpz[1] * c0z[1] + b0[1]) * g[25];
        g[10] = c0x[0] * g[4 ] + 2 * b0[0] * cpx[0];
        g[11] = c0x[1] * g[5 ] + 2 * b0[1] * cpx[1];
        g[22] = c0y[0] * g[16] + 2 * b0[0] * cpy[0];
        g[23] = c0y[1] * g[17] + 2 * b0[1] * cpy[1];
        g[34] = c0z[0] * g[28] + 2 * b0[0] * g[26];
        g[35] = c0z[1] * g[29] + 2 * b0[1] * g[27];
}
static inline void _g0_lj_4d_1020(double *g, double *c0x, double *c0y, double *c0z,
                                  double *cpx, double *cpy, double *cpz,
                                  double *b0, double *b1, double *r0, double *rp)
{
        double rc0[6];
        _make_rc(rc0, c0x, c0y, c0z, r0);
        g[2 ] = rc0[0];
        g[3 ] = rc0[1];
        g[4 ] = cpx[0];
        g[5 ] = cpx[1];
        g[26] = rc0[2];
        g[27] = rc0[3];
        g[28] = cpy[0];
        g[29] = cpy[1];
        g[50] = rc0[4] * g[48];
        g[51] = rc0[5] * g[49];
        g[52] = cpz[0] * g[48];
        g[53] = cpz[1] * g[49];
        g[6 ] = cpx[0] * rc0[0] + b0[0];
        g[7 ] = cpx[1] * rc0[1] + b0[1];
        g[8 ] = cpx[0] * cpx[0] + b1[0];
        g[9 ] = cpx[1] * cpx[1] + b1[1];
        g[30] = cpy[0] * rc0[2] + b0[0];
        g[31] = cpy[1] * rc0[3] + b0[1];
        g[32] = cpy[0] * cpy[0] + b1[0];
        g[33] = cpy[1] * cpy[1] + b1[1];
        g[54] =(cpz[0] * rc0[4] + b0[0]) * g[48];
        g[55] =(cpz[1] * rc0[5] + b0[1]) * g[49];
        g[56] =(cpz[0] * cpz[0] + b1[0]) * g[48];
        g[57] =(cpz[1] * cpz[1] + b1[1]) * g[49];
        g[10] = rc0[0] * g[8 ] + 2 * b0[0] * cpx[0];
        g[11] = rc0[1] * g[9 ] + 2 * b0[1] * cpx[1];
        g[34] = rc0[2] * g[32] + 2 * b0[0] * cpy[0];
        g[35] = rc0[3] * g[33] + 2 * b0[1] * cpy[1];
        g[58] = rc0[4] * g[56] + 2 * b0[0] * g[52];
        g[59] = rc0[5] * g[57] + 2 * b0[1] * g[53];
}
static inline void _g0_lj_4d_0111(double *g, double *c0x, double *c0y, double *c0z,
                                  double *cpx, double *cpy, double *cpz,
                                  double *b0, double *b1, double *r0, double *rp)
{
        double rcp[6];
        _make_rc(rcp, cpx, cpy, cpz, rp);
        g[2 ] = rcp[0];
        g[3 ] = rcp[1];
        g[4 ] = cpx[0];
        g[5 ] = cpx[1];
        g[12] = c0x[0];
        g[13] = c0x[1];
        g[26] = rcp[2];
        g[27] = rcp[3];
        g[28] = cpy[0];
        g[29] = cpy[1];
        g[36] = c0y[0];
        g[37] = c0y[1];
        g[50] = rcp[4] * g[48];
        g[51] = rcp[5] * g[49];
        g[52] = cpz[0] * g[48];
        g[53] = cpz[1] * g[49];
        g[60] = c0z[0] * g[48];
        g[61] = c0z[1] * g[49];
        g[14] = c0x[0] * rcp[0] + b0[0];
        g[15] = c0x[1] * rcp[1] + b0[1];
        g[16] = c0x[0] * cpx[0] + b0[0];
        g[17] = c0x[1] * cpx[1] + b0[1];
        g[6 ] = cpx[0] * rcp[0] + b1[0];
        g[7 ] = cpx[1] * rcp[1] + b1[1];
        g[30] = cpy[0] * rcp[2] + b1[0];
        g[31] = cpy[1] * rcp[3] + b1[1];
        g[38] = c0y[0] * rcp[2] + b0[0];
        g[39] = c0y[1] * rcp[3] + b0[1];
        g[40] = c0y[0] * cpy[0] + b0[0];
        g[41] = c0y[1] * cpy[1] + b0[1];
        g[54] =(cpz[0] * rcp[4] + b1[0]) * g[48];
        g[55] =(cpz[1] * rcp[5] + b1[1]) * g[49];
        g[62] =(c0z[0] * rcp[4] + b0[0]) * g[48];
        g[63] =(c0z[1] * rcp[5] + b0[1]) * g[49];
        g[64] =(c0z[0] * cpz[0] + b0[0]) * g[48];
        g[65] =(c0z[1] * cpz[1] + b0[1]) * g[49];
        g[18] = c0x[0] * g[6 ] + b0[0] * (rcp[0] + cpx[0]);
        g[19] = c0x[1] * g[7 ] + b0[1] * (rcp[1] + cpx[1]);
        g[42] = c0y[0] * g[30] + b0[0] * (rcp[2] + cpy[0]);
        g[43] = c0y[1] * g[31] + b0[1] * (rcp[3] + cpy[1]);
        g[66] = c0z[0] * g[54] + b0[0] * (g[50] + g[52]);
        g[67] = c0z[1] * g[55] + b0[1] * (g[51] + g[53]);
}
static inline void _g0_lj_4d_1110(double *g, double *c0x, double *c0y, double *c0z,
                                  double *cpx, double *cpy, double *cpz,
                                  double *b0, double *b1, double *r0, double *rp)
{
        double rc0[6];
        double rcp[6];
        _make_rc(rc0, c0x, c0y, c0z, r0);
        _make_rc(rcp, cpx, cpy, cpz, rp);
        g[2 ] = rc0[0];
        g[3 ] = rc0[1];
        g[4 ] = rcp[0];
        g[5 ] = rcp[1];
        g[8 ] = cpx[0];
        g[9 ] = cpx[1];
        g[50] = rc0[2];
        g[51] = rc0[3];
        g[52] = rcp[2];
        g[53] = rcp[3];
        g[56] = cpy[0];
        g[57] = cpy[1];
        g[98 ] = rc0[4] * g[96 ];
        g[99 ] = rc0[5] * g[97 ];
        g[100] = rcp[4] * g[96 ];
        g[101] = rcp[5] * g[97 ];
        g[104] = cpz[0] * g[96 ];
        g[105] = cpz[1] * g[97 ];
        g[6 ] = rcp[0] * rc0[0] + b0[0];
        g[7 ] = rcp[1] * rc0[1] + b0[1];
        g[10] = cpx[0] * rc0[0] + b0[0];
        g[11] = cpx[1] * rc0[1] + b0[1];
        g[12] = cpx[0] * rcp[0] + b1[0];
        g[13] = cpx[1] * rcp[1] + b1[1];
        g[54] = rcp[2] * rc0[2] + b0[0];
        g[55] = rcp[3] * rc0[3] + b0[1];
        g[58] = cpy[0] * rc0[2] + b0[0];
        g[59] = cpy[1] * rc0[3] + b0[1];
        g[60] = cpy[0] * rcp[2] + b1[0];
        g[61] = cpy[1] * rcp[3] + b1[1];
        g[102] =(rcp[4] *rc0[4] + b0[0])* g[96 ];
        g[103] =(rcp[5] *rc0[5] + b0[1])* g[97 ];
        g[106] =(cpz[0] *rc0[4] + b0[0])* g[96 ];
        g[107] =(cpz[1] *rc0[5] + b0[1])* g[97 ];
        g[108] =(cpz[0] *rcp[4] + b1[0])* g[96 ];
        g[109] =(cpz[1] *rcp[5] + b1[1])* g[97 ];
        g[14] = rc0[0] * g[12 ] + b0[0] *(rcp[0] + cpx[0]);
        g[15] = rc0[1] * g[13 ] + b0[1] *(rcp[1] + cpx[1]);
        g[62] = rc0[2] * g[60 ] + b0[0] *(rcp[2] + cpy[0]);
        g[63] = rc0[3] * g[61 ] + b0[1] *(rcp[3] + cpy[1]);
        g[110] = rc0[4] *g[108] + b0[0] *(g[100] + g[104]);
        g[111] = rc0[5] *g[109] + b0[1] *(g[101] + g[105]);
}
static inline void _g0_lj_4d_0201(double *g, double *c0x, double *c0y, double *c0z,
                                  double *cpx, double *cpy, double *cpz,
                                  double *b0, double *b1, double *r0, double *rp)
{
        double rcp[6];
        _make_rc(rcp, cpx, cpy, cpz, rp);
        g[2 ] = rcp[0];
        g[3 ] = rcp[1];
        g[18] = c0x[0];
        g[19] = c0x[1];
        g[38] = rcp[2];
        g[39] = rcp[3];
        g[54] = c0y[0];
        g[55] = c0y[1];
        g[74] = rcp[4] * g[72];
        g[75] = rcp[5] * g[73];
        g[90] = c0z[0] * g[72];
        g[91] = c0z[1] * g[73];
        g[4 ] = rcp[0] * rcp[0] + b1[0];
        g[5 ] = rcp[1] * rcp[1] + b1[1];
        g[20] = rcp[0] * c0x[0] + b0[0];
        g[21] = rcp[1] * c0x[1] + b0[1];
        g[40] = rcp[2] * rcp[2] + b1[0];
        g[41] = rcp[3] * rcp[3] + b1[1];
        g[56] = rcp[2] * c0y[0] + b0[0];
        g[57] = rcp[3] * c0y[1] + b0[1];
        g[76] =(rcp[4] * rcp[4] + b1[0])* g[72];
        g[77] =(rcp[5] * rcp[5] + b1[1])* g[73];
        g[92] =(rcp[4] * c0z[0] + b0[0])* g[72];
        g[93] =(rcp[5] * c0z[1] + b0[1])* g[73];
        g[22] = c0x[0] * g[4 ] + 2 * b0[0] * rcp[0];
        g[23] = c0x[1] * g[5 ] + 2 * b0[1] * rcp[1];
        g[58] = c0y[0] * g[40] + 2 * b0[0] * rcp[2];
        g[59] = c0y[1] * g[41] + 2 * b0[1] * rcp[3];
        g[94] = c0z[0] * g[76] + 2 * b0[0] * g[74];
        g[95] = c0z[1] * g[77] + 2 * b0[1] * g[75];
}
static inline void _g0_lj_4d_1200(double *g, double *c0x, double *c0y, double *c0z,
                                  double *cpx, double *cpy, double *cpz,
                                  double *b0, double *b1, double *r0, double *rp)
{
        double rc0[6];
        double rcp[6];
        _make_rc(rc0, c0x, c0y, c0z, r0);
        _make_rc(rcp, cpx, cpy, cpz, rp);
        g[2 ] = rc0[0];
        g[3 ] = rc0[1];
        g[4 ] = rcp[0];
        g[5 ] = rcp[1];
        g[6 ] = rcp[0] * rc0[0] + b0[0];
        g[7 ] = rcp[1] * rc0[1] + b0[1];
        g[8 ] = rcp[0] * rcp[0] + b1[0];
        g[9 ] = rcp[1] * rcp[1] + b1[1];
        g[10] = rc0[0] * g[8] + 2 * b0[0] * rcp[0];
        g[11] = rc0[1] * g[9] + 2 * b0[1] * rcp[1];
        g[74] = rc0[2];
        g[75] = rc0[3];
        g[76] = rcp[2];
        g[77] = rcp[3];
        g[78] = rcp[2] * rc0[2] + b0[0];
        g[79] = rcp[3] * rc0[3] + b0[1];
        g[80] = rcp[2] * rcp[2] + b1[0];
        g[81] = rcp[3] * rcp[3] + b1[1];
        g[82] = rc0[2] * g[80] + 2 * b0[0] * rcp[2];
        g[83] = rc0[3] * g[81] + 2 * b0[1] * rcp[3];
        g[146] = rc0[4] * g[144];
        g[147] = rc0[5] * g[145];
        g[148] = rcp[4] * g[144];
        g[149] = rcp[5] * g[145];
        g[150] =(rcp[4] * rc0[4] +b0[0])* g[144];
        g[151] =(rcp[5] * rc0[5] +b0[1])* g[145];
        g[152] =(rcp[4] * rcp[4] +b1[0])* g[144];
        g[153] =(rcp[5] * rcp[5] +b1[1])* g[145];
        g[154] = rc0[4] * g[152] + 2 * b0[0] * g[148];
        g[155] = rc0[5] * g[153] + 2 * b0[1] * g[149];
}
static inline void _g0_lj_4d_0012(double *g, double *c0x, double *c0y, double *c0z,
                                  double *cpx, double *cpy, double *cpz,
                                  double *b0, double *b1)
{
        g[2 ] = cpx[0];
        g[3 ] = cpx[1];
        g[4 ] = c0x[0];
        g[5 ] = c0x[1];
        g[14] = cpy[0];
        g[15] = cpy[1];
        g[16] = c0y[0];
        g[17] = c0y[1];
        g[26] = cpz[0] * g[24];
        g[27] = cpz[1] * g[25];
        g[28] = c0z[0] * g[24];
        g[29] = c0z[1] * g[25];
        g[6 ] = cpx[0] * c0x[0] + b0[0];
        g[7 ] = cpx[1] * c0x[1] + b0[1];
        g[8 ] = c0x[0] * c0x[0] + b1[0];
        g[9 ] = c0x[1] * c0x[1] + b1[1];
        g[18] = cpy[0] * c0y[0] + b0[0];
        g[19] = cpy[1] * c0y[1] + b0[1];
        g[20] = c0y[0] * c0y[0] + b1[0];
        g[21] = c0y[1] * c0y[1] + b1[1];
        g[30] =(cpz[0] * c0z[0] + b0[0]) * g[24];
        g[31] =(cpz[1] * c0z[1] + b0[1]) * g[25];
        g[32] =(c0z[0] * c0z[0] + b1[0]) * g[24];
        g[33] =(c0z[1] * c0z[1] + b1[1]) * g[25];
        g[10] = cpx[0] * g[8 ] + 2 * b0[0] * c0x[0];
        g[11] = cpx[1] * g[9 ] + 2 * b0[1] * c0x[1];
        g[22] = cpy[0] * g[20] + 2 * b0[0] * c0y[0];
        g[23] = cpy[1] * g[21] + 2 * b0[1] * c0y[1];
        g[34] = cpz[0] * g[32] + 2 * b0[0] * g[28];
        g[35] = cpz[1] * g[33] + 2 * b0[1] * g[29];
}
static inline void _g0_lj_4d_1011(double *g, double *c0x, double *c0y, double *c0z,
                                  double *cpx, double *cpy, double *cpz,
                                  double *b0, double *b1, double *r0, double *rp)
{
        double rc0[6];
        _make_rc(rc0, c0x, c0y, c0z, r0);
        g[2 ] = rc0[0];
        g[3 ] = rc0[1];
        g[4 ] = cpx[0];
        g[5 ] = cpx[1];
        g[8 ] = c0x[0];
        g[9 ] = c0x[1];
        g[26] = rc0[2];
        g[27] = rc0[3];
        g[28] = cpy[0];
        g[29] = cpy[1];
        g[32] = c0y[0];
        g[33] = c0y[1];
        g[50] = rc0[4] * g[48];
        g[51] = rc0[5] * g[49];
        g[52] = cpz[0] * g[48];
        g[53] = cpz[1] * g[49];
        g[56] = c0z[0] * g[48];
        g[57] = c0z[1] * g[49];
        g[6 ] = cpx[0] * rc0[0] + b0[0];
        g[7 ] = cpx[1] * rc0[1] + b0[1];
        g[10] = c0x[0] * rc0[0] + b1[0];
        g[11] = c0x[1] * rc0[1] + b1[1];
        g[12] = c0x[0] * cpx[0] + b0[0];
        g[13] = c0x[1] * cpx[1] + b0[1];
        g[30] = cpy[0] * rc0[2] + b0[0];
        g[31] = cpy[1] * rc0[3] + b0[1];
        g[34] = c0y[0] * rc0[2] + b1[0];
        g[35] = c0y[1] * rc0[3] + b1[1];
        g[36] = c0y[0] * cpy[0] + b0[0];
        g[37] = c0y[1] * cpy[1] + b0[1];
        g[54] =(cpz[0] * rc0[4] + b0[0])* g[48];
        g[55] =(cpz[1] * rc0[5] + b0[1])* g[49];
        g[58] =(c0z[0] * rc0[4] + b1[0])* g[48];
        g[59] =(c0z[1] * rc0[5] + b1[1])* g[49];
        g[60] =(c0z[0] * cpz[0] + b0[0])* g[48];
        g[61] =(c0z[1] * cpz[1] + b0[1])* g[49];
        g[14] = cpx[0] * g[10] + b0[0] *(rc0[0] + c0x[0]);
        g[15] = cpx[1] * g[11] + b0[1] *(rc0[1] + c0x[1]);
        g[38] = cpy[0] * g[34] + b0[0] *(rc0[2] + c0y[0]);
        g[39] = cpy[1] * g[35] + b0[1] *(rc0[3] + c0y[1]);
        g[62] = cpz[0] * g[58] + b0[0] *(g[50] + g[56]);
        g[63] = cpz[1] * g[59] + b0[1] *(g[51] + g[57]);
}
static inline void _g0_lj_4d_2010(double *g, double *c0x, double *c0y, double *c0z,
                                  double *cpx, double *cpy, double *cpz,
                                  double *b0, double *b1, double *r0, double *rp)
{
        double rc0[6];
        _make_rc(rc0, c0x, c0y, c0z, r0);
        g[2 ] = rc0[0];
        g[3 ] = rc0[1];
        g[6 ] = cpx[0];
        g[7 ] = cpx[1];
        g[38] = rc0[2];
        g[39] = rc0[3];
        g[42] = cpy[0];
        g[43] = cpy[1];
        g[74] = rc0[4] * g[72];
        g[75] = rc0[5] * g[73];
        g[78] = cpz[0] * g[72];
        g[79] = cpz[1] * g[73];
        g[4 ] = rc0[0] * rc0[0] + b1[0];
        g[5 ] = rc0[1] * rc0[1] + b1[1];
        g[8 ] = cpx[0] * rc0[0] + b0[0];
        g[9 ] = cpx[1] * rc0[1] + b0[1];
        g[40] = rc0[2] * rc0[2] + b1[0];
        g[41] = rc0[3] * rc0[3] + b1[1];
        g[44] = cpy[0] * rc0[2] + b0[0];
        g[45] = cpy[1] * rc0[3] + b0[1];
        g[76] =(rc0[4] * rc0[4] + b1[0]) * g[72];
        g[77] =(rc0[5] * rc0[5] + b1[1]) * g[73];
        g[80] =(cpz[0] * rc0[4] + b0[0]) * g[72];
        g[81] =(cpz[1] * rc0[5] + b0[1]) * g[73];
        g[10] = cpx[0] * g[4 ] + 2 * b0[0] * rc0[0];
        g[11] = cpx[1] * g[5 ] + 2 * b0[1] * rc0[1];
        g[46] = cpy[0] * g[40] + 2 * b0[0] * rc0[2];
        g[47] = cpy[1] * g[41] + 2 * b0[1] * rc0[3];
        g[82] = cpz[0] * g[76] + 2 * b0[0] * g[74];
        g[83] = cpz[1] * g[77] + 2 * b0[1] * g[75];
}
static inline void _g0_lj_4d_0102(double *g, double *c0x, double *c0y, double *c0z,
                                  double *cpx, double *cpy, double *cpz,
                                  double *b0, double *b1, double *r0, double *rp)
{
        double rcp[6];
        _make_rc(rcp, cpx, cpy, cpz, rp);
        g[2 ] = rcp[0];
        g[3 ] = rcp[1];
        g[8 ] = c0x[0];
        g[9 ] = c0x[1];
        g[26] = rcp[2];
        g[27] = rcp[3];
        g[32] = c0y[0];
        g[33] = c0y[1];
        g[50] = rcp[4] * g[48];
        g[51] = rcp[5] * g[49];
        g[56] = c0z[0] * g[48];
        g[57] = c0z[1] * g[49];
        g[10] = rcp[0] * c0x[0] + b0[0];
        g[11] = rcp[1] * c0x[1] + b0[1];
        g[16] = c0x[0] * c0x[0] + b1[0];
        g[17] = c0x[1] * c0x[1] + b1[1];
        g[34] = rcp[2] * c0y[0] + b0[0];
        g[35] = rcp[3] * c0y[1] + b0[1];
        g[40] = c0y[0] * c0y[0] + b1[0];
        g[41] = c0y[1] * c0y[1] + b1[1];
        g[58] =(rcp[4] * c0z[0] + b0[0]) * g[48];
        g[59] =(rcp[5] * c0z[1] + b0[1]) * g[49];
        g[64] =(c0z[0] * c0z[0] + b1[0]) * g[48];
        g[65] =(c0z[1] * c0z[1] + b1[1]) * g[49];
        g[18] = rcp[0] * g[16] + 2 * b0[0] * c0x[0];
        g[19] = rcp[1] * g[17] + 2 * b0[1] * c0x[1];
        g[42] = rcp[2] * g[40] + 2 * b0[0] * c0y[0];
        g[43] = rcp[3] * g[41] + 2 * b0[1] * c0y[1];
        g[66] = rcp[4] * g[64] + 2 * b0[0] * g[56];
        g[67] = rcp[5] * g[65] + 2 * b0[1] * g[57];
}
static inline void _g0_lj_4d_1101(double *g, double *c0x, double *c0y, double *c0z,
                                  double *cpx, double *cpy, double *cpz,
                                  double *b0, double *b1, double *r0, double *rp)
{
        double rc0[6];
        double rcp[6];
        _make_rc(rc0, c0x, c0y, c0z, r0);
        _make_rc(rcp, cpx, cpy, cpz, rp);
        g[2 ] = rc0[0];
        g[3 ] = rc0[1];
        g[4 ] = rcp[0];
        g[5 ] = rcp[1];
        g[16] = c0x[0];
        g[17] = c0x[1];
        g[50] = rc0[2];
        g[51] = rc0[3];
        g[52] = rcp[2];
        g[53] = rcp[3];
        g[64] = c0y[0];
        g[65] = c0y[1];
        g[98 ] = rc0[4] * g[96];
        g[99 ] = rc0[5] * g[97];
        g[100] = rcp[4] * g[96];
        g[101] = rcp[5] * g[97];
        g[112] = c0z[0] * g[96];
        g[113] = c0z[1] * g[97];
        g[6 ] = rcp[0] * rc0[0] + b0[0];
        g[7 ] = rcp[1] * rc0[1] + b0[1];
        g[18] = c0x[0] * rc0[0] + b1[0];
        g[19] = c0x[1] * rc0[1] + b1[1];
        g[20] = c0x[0] * rcp[0] + b0[0];
        g[21] = c0x[1] * rcp[1] + b0[1];
        g[54] = rcp[2] * rc0[2] + b0[0];
        g[55] = rcp[3] * rc0[3] + b0[1];
        g[66] = c0y[0] * rc0[2] + b1[0];
        g[67] = c0y[1] * rc0[3] + b1[1];
        g[68] = c0y[0] * rcp[2] + b0[0];
        g[69] = c0y[1] * rcp[3] + b0[1];
        g[102] = (rcp[4] * rc0[4] + b0[0])* g[96];
        g[103] = (rcp[5] * rc0[5] + b0[1])* g[97];
        g[114] = (c0z[0] * rc0[4] + b1[0])* g[96];
        g[115] = (c0z[1] * rc0[5] + b1[1])* g[97];
        g[116] = (c0z[0] * rcp[4] + b0[0])* g[96];
        g[117] = (c0z[1] * rcp[5] + b0[1])* g[97];
        g[22] = rcp[0] * g[18] + b0[0] *(rc0[0] + c0x[0]);
        g[23] = rcp[1] * g[19] + b0[1] *(rc0[1] + c0x[1]);
        g[70] = rcp[2] * g[66] + b0[0] *(rc0[2] + c0y[0]);
        g[71] = rcp[3] * g[67] + b0[1] *(rc0[3] + c0y[1]);
        g[118] = rcp[4] *g[114] + b0[0] *(g[98] + g[112]);
        g[119] = rcp[5] *g[115] + b0[1] *(g[99] + g[113]);
}
static inline void _g0_lj_4d_2100(double *g, double *c0x, double *c0y, double *c0z,
                                  double *cpx, double *cpy, double *cpz,
                                  double *b0, double *b1, double *r0, double *rp)
{
        double rc0[6];
        double rcp[6];
        _make_rc(rc0, c0x, c0y, c0z, r0);
        _make_rc(rcp, cpx, cpy, cpz, rp);
        g[2 ] = rc0[0];
        g[3 ] = rc0[1];
        g[6 ] = rcp[0];
        g[7 ] = rcp[1];
        g[74] = rc0[2];
        g[75] = rc0[3];
        g[78] = rcp[2];
        g[79] = rcp[3];
        g[146] = rc0[4] * g[144];
        g[147] = rc0[5] * g[145];
        g[150] = rcp[4] * g[144];
        g[151] = rcp[5] * g[145];
        g[4 ] = rc0[0] * rc0[0] + b1[0];
        g[5 ] = rc0[1] * rc0[1] + b1[1];
        g[8 ] = rcp[0] * rc0[0] + b0[0];
        g[9 ] = rcp[1] * rc0[1] + b0[1];
        g[76] = rc0[2] * rc0[2] + b1[0];
        g[77] = rc0[3] * rc0[3] + b1[1];
        g[80] = rcp[2] * rc0[2] + b0[0];
        g[81] = rcp[3] * rc0[3] + b0[1];
        g[148] = (rc0[4] * rc0[4] + b1[0])* g[144];
        g[149] = (rc0[5] * rc0[5] + b1[1])* g[145];
        g[152] = (rcp[4] * rc0[4] + b0[0])* g[144];
        g[153] = (rcp[5] * rc0[5] + b0[1])* g[145];
        g[10 ] = rcp[0] * g[4  ] + 2 * b0[0] * rc0[0];
        g[11 ] = rcp[1] * g[5  ] + 2 * b0[1] * rc0[1];
        g[82 ] = rcp[2] * g[76 ] + 2 * b0[0] * rc0[2];
        g[83 ] = rcp[3] * g[77 ] + 2 * b0[1] * rc0[3];
        g[154] = rcp[4] * g[148] + 2 * b0[0] * g[146];
        g[155] = rcp[5] * g[149] + 2 * b0[1] * g[147];
}
/************** end special g0_4d results *************/



void CINTg0_2e_lj2d4d_simd1(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        int nmax = envs->li_ceil + envs->lj_ceil;
        int mmax = envs->lk_ceil + envs->ll_ceil;
        switch (nmax) {
        case 0: switch(mmax) {
                case 0: goto _g0_4d_default; // ssss
                case 1: switch (envs->lk_ceil) {
                        case 0: _g0_lj_4d_0001(g, bc->c0px, bc->c0py, bc->c0pz, envs->rkrl); return;
                        case 1: _g0_lj_4d_1000(g, bc->c0px, bc->c0py, bc->c0pz, envs->rkrl); return;
                        default: goto error; }
                case 2: switch (envs->lk_ceil) {
                        case 0: _g0_lj_4d_0002(g, bc->c0px, bc->c0py, bc->c0pz, bc->b01, envs->rkrl); return;
                        case 1: _g0_lj_4d_1001(g, bc->c0px, bc->c0py, bc->c0pz, bc->b01, envs->rkrl); return;
                        case 2: _g0_lj_4d_2000(g, bc->c0px, bc->c0py, bc->c0pz, bc->b01, envs->rkrl); return;
                        default: goto error; }
                case 3: switch (envs->lk_ceil) {
                        case 0: _g0_lj_4d_0003(g, bc->c0px, bc->c0py, bc->c0pz, bc->b01, envs->rkrl); return;
                        case 1: _g0_lj_4d_1002(g, bc->c0px, bc->c0py, bc->c0pz, bc->b01, envs->rkrl); return;
                        case 2: _g0_lj_4d_2001(g, bc->c0px, bc->c0py, bc->c0pz, bc->b01, envs->rkrl); return;
                        case 3: _g0_lj_4d_3000(g, bc->c0px, bc->c0py, bc->c0pz, bc->b01, envs->rkrl); return;
                        default: goto error; }
                default: goto _g0_4d_default; }
        case 1: switch(mmax) {
                case 0: switch (envs->li_ceil) {
                        case 0: _g0_lj_4d_0001(g, bc->c00x, bc->c00y, bc->c00z, envs->rirj); return;
                        case 1: _g0_lj_4d_1000(g, bc->c00x, bc->c00y, bc->c00z, envs->rirj); return;
                        default: goto error; }
                case 1: switch (envs->lk_ceil) {
                        case 0: switch (envs->li_ceil) {
                                case 0: _g0_lj_4d_0011(g, bc->c00x, bc->c00y, bc->c00z, bc->c0px, bc->c0py, bc->c0pz, bc->b00, envs->rirj, envs->rkrl); return;
                                case 1: _g0_lj_4d_1010(g, bc->c00x, bc->c00y, bc->c00z, bc->c0px, bc->c0py, bc->c0pz, bc->b00, envs->rirj, envs->rkrl); return;
                                default: goto error; }
                        case 1: switch (envs->li_ceil) {
                                case 0: _g0_lj_4d_0101(g, bc->c00x, bc->c00y, bc->c00z, bc->c0px, bc->c0py, bc->c0pz, bc->b00, envs->rirj, envs->rkrl); return;
                                case 1: _g0_lj_4d_1100(g, bc->c00x, bc->c00y, bc->c00z, bc->c0px, bc->c0py, bc->c0pz, bc->b00, envs->rirj, envs->rkrl); return;
                                default: goto error; }
                        default: goto error; }
                case 2: switch (envs->lk_ceil) {
                        case 0: switch (envs->li_ceil) {
                                case 0: _g0_lj_4d_0021(g, bc->c00x, bc->c00y, bc->c00z, bc->c0px, bc->c0py, bc->c0pz, bc->b00, bc->b01); return;
                                case 1: _g0_lj_4d_1020(g, bc->c00x, bc->c00y, bc->c00z, bc->c0px, bc->c0py, bc->c0pz, bc->b00, bc->b01, envs->rirj, envs->rkrl); return;
                                default: goto error; }
                        case 1: switch (envs->li_ceil) {
                                case 0: _g0_lj_4d_0111(g, bc->c00x, bc->c00y, bc->c00z, bc->c0px, bc->c0py, bc->c0pz, bc->b00, bc->b01, envs->rirj, envs->rkrl); return;
                                case 1: _g0_lj_4d_1110(g, bc->c00x, bc->c00y, bc->c00z, bc->c0px, bc->c0py, bc->c0pz, bc->b00, bc->b01, envs->rirj, envs->rkrl); return;
                                default: goto error; }
                        case 2: switch (envs->li_ceil) {
                                case 0: _g0_lj_4d_0201(g, bc->c00x, bc->c00y, bc->c00z, bc->c0px, bc->c0py, bc->c0pz, bc->b00, bc->b01, envs->rirj, envs->rkrl); return;
                                case 1: _g0_lj_4d_1200(g, bc->c00x, bc->c00y, bc->c00z, bc->c0px, bc->c0py, bc->c0pz, bc->b00, bc->b01, envs->rirj, envs->rkrl); return;
                                default: goto error; }
                        default: goto error; }
                default: goto _g0_4d_default; }
        case 2: switch(mmax) {
                case 0: switch (envs->li_ceil) {
                        case 0: _g0_lj_4d_0002(g, bc->c00x, bc->c00y, bc->c00z, bc->b10, envs->rirj); return;
                        case 1: _g0_lj_4d_1001(g, bc->c00x, bc->c00y, bc->c00z, bc->b10, envs->rirj); return;
                        case 2: _g0_lj_4d_2000(g, bc->c00x, bc->c00y, bc->c00z, bc->b10, envs->rirj); return;
                        default: goto error; }
                case 1: switch (envs->lk_ceil) {
                        case 0: switch (envs->li_ceil) {
                                case 0: _g0_lj_4d_0012(g, bc->c00x, bc->c00y, bc->c00z, bc->c0px, bc->c0py, bc->c0pz, bc->b00, bc->b10); return;
                                case 1: _g0_lj_4d_1011(g, bc->c00x, bc->c00y, bc->c00z, bc->c0px, bc->c0py, bc->c0pz, bc->b00, bc->b10, envs->rirj, envs->rkrl); return;
                                case 2: _g0_lj_4d_2010(g, bc->c00x, bc->c00y, bc->c00z, bc->c0px, bc->c0py, bc->c0pz, bc->b00, bc->b10, envs->rirj, envs->rkrl); return;
                                default: goto error; }
                        case 1: switch (envs->li_ceil) {
                                case 0: _g0_lj_4d_0102(g, bc->c00x, bc->c00y, bc->c00z, bc->c0px, bc->c0py, bc->c0pz, bc->b00, bc->b10, envs->rirj, envs->rkrl); return;
                                case 1: _g0_lj_4d_1101(g, bc->c00x, bc->c00y, bc->c00z, bc->c0px, bc->c0py, bc->c0pz, bc->b00, bc->b10, envs->rirj, envs->rkrl); return;
                                case 2: _g0_lj_4d_2100(g, bc->c00x, bc->c00y, bc->c00z, bc->c0px, bc->c0py, bc->c0pz, bc->b00, bc->b10, envs->rirj, envs->rkrl); return;
                                default: goto error; }
                        default: goto error; }
                default: goto _g0_4d_default; }
        case 3: switch(mmax) {
                case 0: switch (envs->li_ceil) {
                        case 0: _g0_lj_4d_0003(g, bc->c00x, bc->c00y, bc->c00z, bc->b10, envs->rirj); return;
                        case 1: _g0_lj_4d_1002(g, bc->c00x, bc->c00y, bc->c00z, bc->b10, envs->rirj); return;
                        case 2: _g0_lj_4d_2001(g, bc->c00x, bc->c00y, bc->c00z, bc->b10, envs->rirj); return;
                        case 3: _g0_lj_4d_3000(g, bc->c00x, bc->c00y, bc->c00z, bc->b10, envs->rirj); return;
                        default: goto error; }
                default: goto _g0_4d_default; }
        default:
_g0_4d_default:
                        CINTg0_2e_2d_simd1(g, bc, envs);
                        CINTg0_lj_4d_simd1(g, envs);
                        return;
        }
error:
        fprintf(stderr, "Dimension error for CINTg0_2e_lj2d4d_simd1: iklj = %d %d %d %d\n",
                envs->li_ceil, envs->lk_ceil, envs->ll_ceil, envs->lj_ceil);
        exit(1);
}

void CINTg0_2e_kj2d4d_simd1(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        CINTg0_2e_2d_simd1(g, bc, envs);
        CINTg0_kj_4d_simd1(g, envs);
}
void CINTg0_2e_ik2d4d_simd1(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        CINTg0_2e_2d_simd1(g, bc, envs);
        CINTg0_ik_4d_simd1(g, envs);
}
void CINTg0_2e_il2d4d_simd1(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        CINTg0_2e_2d_simd1(g, bc, envs);
        CINTg0_il_4d_simd1(g, envs);
}

#ifdef WITH_F12
void CINTg0_2e_stg_lj2d4d_simd1(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        CINTg0_2e_2d_simd1(g, bc, envs);
        CINTg0_lj_4d_simd1(g, envs);
}
#endif


/*
 * g[i,k,l,j] = < ik | lj > = ( i j | k l )
 */
void CINTg0_2e_simd1(double *g, Rys2eT *bc, CINTEnvVars *envs, int idsimd)
{
        double aij, akl, a0, a1, fac1;
        ALIGNMM double x[SIMDD];
        double *rij = envs->rij;
        double *rkl = envs->rkl;
        double rijrkl[3];
        double rijrx[3];
        double rklrx[3];
        double *u = bc->u;
        double *w = bc->w;

        int i;
        aij = envs->ai[idsimd] + envs->aj[idsimd];
        akl = envs->ak[idsimd] + envs->al[idsimd];
        a1 = aij * akl;
        a0 = a1 / (aij + akl);
#ifdef WITH_RANGE_COULOMB
        const double omega = envs->env[PTR_RANGE_OMEGA];
        double theta = 1;
        if (omega > 0) {
// For long-range part of range-separated Coulomb operator
                theta = omega * omega / (omega * omega + a0);
                a0 *= theta;
        }
#endif
        fac1 = sqrt(a0 / (a1 * a1 * a1)) * envs->fac[idsimd];

        rijrkl[0] = rij[0*SIMDD+idsimd] - rkl[0*SIMDD+idsimd];
        rijrkl[1] = rij[1*SIMDD+idsimd] - rkl[1*SIMDD+idsimd];
        rijrkl[2] = rij[2*SIMDD+idsimd] - rkl[2*SIMDD+idsimd];
        x[0] = a0 * SQUARE(rijrkl);
        CINTrys_roots(envs->nrys_roots, x, u, w, 1);

        double *gx = g;
        double *gy = g + envs->g_size;
        double *gz = g + envs->g_size * 2;
        for (i = 0; i < envs->nrys_roots; i++) {
                gx[i] = 1;
                gy[i] = 1;
                gz[i] = w[i*SIMDD] * fac1;
        }
#ifdef WITH_RANGE_COULOMB
        if (omega > 0) {
                /* u[:] = tau^2 / (1 - tau^2)
                 * transform u[:] to theta^-1 tau^2 / (theta^-1 - tau^2)
                 * so the rest code can be reused.
                 */
                for (i = 0; i < envs->nrys_roots; i++) {
                        u[i*SIMDD] /= u[i*SIMDD] + 1 - u[i*SIMDD] * theta;
                }
        }
#endif
        if (envs->g_size == 1) {
                return;
        }

        double u2, div, tmp1, tmp2, tmp3, tmp4;
        double *b00 = bc->b00;
        double *b10 = bc->b10;
        double *b01 = bc->b01;
        double *c00x = bc->c00x;
        double *c00y = bc->c00y;
        double *c00z = bc->c00z;
        double *c0px = bc->c0px;
        double *c0py = bc->c0py;
        double *c0pz = bc->c0pz;

        rijrx[0] = rij[0*SIMDD+idsimd] - envs->rx_in_rijrx[0];
        rijrx[1] = rij[1*SIMDD+idsimd] - envs->rx_in_rijrx[1];
        rijrx[2] = rij[2*SIMDD+idsimd] - envs->rx_in_rijrx[2];
        rklrx[0] = rkl[0*SIMDD+idsimd] - envs->rx_in_rklrx[0];
        rklrx[1] = rkl[1*SIMDD+idsimd] - envs->rx_in_rklrx[1];
        rklrx[2] = rkl[2*SIMDD+idsimd] - envs->rx_in_rklrx[2];
        for (i = 0; i < envs->nrys_roots; i++) {
                /*
                 *t2 = u(i)/(1+u(i))
                 *u2 = aij*akl/(aij+akl)*t2/(1-t2)
                 */
                u2 = a0 * u[i*SIMDD];
                div = 1 / (u2 * (aij + akl) + a1);
                tmp1 = u2 * div;
                tmp4 = .5 * div;
                b00[i] = 0.5 * tmp1;
                tmp2 = tmp1 * akl;
                tmp3 = tmp1 * aij;
                b10[i] = b00[i] + tmp4 * akl;
                b01[i] = b00[i] + tmp4 * aij;
                c00x[i] = rijrx[0] - tmp2 * rijrkl[0];
                c00y[i] = rijrx[1] - tmp2 * rijrkl[1];
                c00z[i] = rijrx[2] - tmp2 * rijrkl[2];
                c0px[i] = rklrx[0] + tmp3 * rijrkl[0];
                c0py[i] = rklrx[1] + tmp3 * rijrkl[1];
                c0pz[i] = rklrx[2] + tmp3 * rijrkl[2];
        }

        (*envs->f_g0_2d4d_simd1)(g, bc, envs);
}


/*
 * ( \nabla i j | kl )
 */
void CINTnabla1i_2e_simd1(double *f, double *g,
                          int li, int lj, int lk, int ll, CINTEnvVars *envs)
{
        int nroots = envs->nrys_roots;
        int di = envs->g_stride_i;
        double ai2 = -2 * envs->ai[0];
        DEF_GXYZ(double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        if (nroots == 1) { // nabla_i ssss
                assert(li == 0);
                assert(lj == 0);
                assert(lk == 0);
                assert(ll == 0);
                fx[0] = ai2 * gx[di];
                fy[0] = ai2 * gy[di];
                fz[0] = ai2 * gz[di];
        } else {
                int i, j, k, l, n, ptr;
                int dk = envs->g_stride_k;
                int dl = envs->g_stride_l;
                int dj = envs->g_stride_j;
                double *p1x = gx - di;
                double *p1y = gy - di;
                double *p1z = gz - di;
                double *p2x = gx + di;
                double *p2y = gy + di;
                double *p2z = gz + di;
                __m128d r0, r1, r2, r3;
                r0 = _mm_load1_pd(&ai2);
                for (j = 0; j <= lj; j++)
                for (l = 0; l <= ll; l++)
                for (k = 0; k <= lk; k++) {
                        ptr = dj * j + dl * l + dk * k;
                        //f(...,0,...) = -2*ai*g(...,1,...)
                        for (n = ptr; n < ptr+nroots; n += 2) {
                                //fx[n] = ai2 * p2x[n];
                                //fy[n] = ai2 * p2y[n];
                                //fz[n] = ai2 * p2z[n];
                                r3 = _mm_loadu_pd(&p2x[n]);
                                r3 = _mm_mul_pd(r3, r0);
                                _mm_storeu_pd(&fx[n], r3);
                                r3 = _mm_loadu_pd(&p2y[n]);
                                r3 = _mm_mul_pd(r3, r0);
                                _mm_storeu_pd(&fy[n], r3);
                                r3 = _mm_loadu_pd(&p2z[n]);
                                r3 = _mm_mul_pd(r3, r0);
                                _mm_storeu_pd(&fz[n], r3);
                        }
                        ptr += di;
                        //f(...,i,...) = i*g(...,i-1,...)-2*ai*g(...,i+1,...)
                        for (i = 1; i <= li; i++) {
                                for (n = ptr; n < ptr+nroots; n += 2) {
                                        //They are always even iterations
                                        //fx[n] = i*p1x[n] + ai2*p2x[n];
                                        //fy[n] = i*p1y[n] + ai2*p2y[n];
                                        //fz[n] = i*p1z[n] + ai2*p2z[n];
                                        r1 = _mm_set1_pd(i);
                                        r2 = _mm_loadu_pd(&p1x[n]);
                                        r3 = _mm_loadu_pd(&p2x[n]);
                                        r2 = _mm_mul_pd(r2, r1);
                                        r3 = _mm_mul_pd(r3, r0);
                                        r2 = _mm_add_pd(r2, r3);
                                        _mm_storeu_pd(&fx[n], r2);
                                        r2 = _mm_loadu_pd(&p1y[n]);
                                        r3 = _mm_loadu_pd(&p2y[n]);
                                        r2 = _mm_mul_pd(r2, r1);
                                        r3 = _mm_mul_pd(r3, r0);
                                        r2 = _mm_add_pd(r2, r3);
                                        _mm_storeu_pd(&fy[n], r2);
                                        r2 = _mm_loadu_pd(&p1z[n]);
                                        r3 = _mm_loadu_pd(&p2z[n]);
                                        r2 = _mm_mul_pd(r2, r1);
                                        r3 = _mm_mul_pd(r3, r0);
                                        r2 = _mm_add_pd(r2, r3);
                                        _mm_storeu_pd(&fz[n], r2);
                                }
                                ptr += di;
                        }
                }
        }
}


/*
 * ( i \nabla j | kl )
 */
void CINTnabla1j_2e_simd1(double *f, double *g,
                          int li, int lj, int lk, int ll, CINTEnvVars *envs)
{
        int nroots = envs->nrys_roots;
        int dj = envs->g_stride_j;
        double aj2 = -2 * envs->aj[0];
        DEF_GXYZ(double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        if (nroots == 1) { // nabla_j ssss
                assert(li == 0);
                assert(lj == 0);
                assert(lk == 0);
                assert(ll == 0);
                fx[0] = aj2 * gx[dj];
                fy[0] = aj2 * gy[dj];
                fz[0] = aj2 * gz[dj];
        } else {
                int i, j, k, l, n, ptr;
                int di = envs->g_stride_i;
                int dk = envs->g_stride_k;
                int dl = envs->g_stride_l;
                double *p1x = gx - dj;
                double *p1y = gy - dj;
                double *p1z = gz - dj;
                double *p2x = gx + dj;
                double *p2y = gy + dj;
                double *p2z = gz + dj;
                __m128d r0, r1, r2, r3;
                r0 = _mm_load1_pd(&aj2);
                //f(...,0,...) = -2*aj*g(...,1,...)
                for (l = 0; l <= ll; l++) {
                for (k = 0; k <= lk; k++) {
                        ptr = dl * l + dk * k;
                        for (i = 0; i <= li; i++) {
                                for (n = ptr; n < ptr+nroots; n += 2) {
                                        //fx[n] = aj2 * p2x[n];
                                        //fy[n] = aj2 * p2y[n];
                                        //fz[n] = aj2 * p2z[n];
                                        r3 = _mm_loadu_pd(&p2x[n]);
                                        r3 = _mm_mul_pd(r3, r0);
                                        _mm_storeu_pd(&fx[n], r3);
                                        r3 = _mm_loadu_pd(&p2y[n]);
                                        r3 = _mm_mul_pd(r3, r0);
                                        _mm_storeu_pd(&fy[n], r3);
                                        r3 = _mm_loadu_pd(&p2z[n]);
                                        r3 = _mm_mul_pd(r3, r0);
                                        _mm_storeu_pd(&fz[n], r3);
                                }
                                ptr += di;
                        }
                } }
                //f(...,j,...) = j*g(...,j-1,...)-2*aj*g(...,j+1,...)
                for (j = 1; j <= lj; j++) {
                        r1 = _mm_set1_pd(j);
                        for (l = 0; l <= ll; l++) {
                        for (k = 0; k <= lk; k++) {
                                ptr = dj * j + dl * l + dk * k;
                                for (i = 0; i <= li; i++) {
                                        for (n = ptr; n < ptr+nroots; n += 2) {
                                                //fx[n] = j*p1x[n] + aj2*p2x[n];
                                                //fy[n] = j*p1y[n] + aj2*p2y[n];
                                                //fz[n] = j*p1z[n] + aj2*p2z[n];
                                                r2 = _mm_loadu_pd(&p1x[n]);
                                                r3 = _mm_loadu_pd(&p2x[n]);
                                                r2 = _mm_mul_pd(r2, r1);
                                                r3 = _mm_mul_pd(r3, r0);
                                                r2 = _mm_add_pd(r2, r3);
                                                _mm_storeu_pd(&fx[n], r2);
                                                r2 = _mm_loadu_pd(&p1y[n]);
                                                r3 = _mm_loadu_pd(&p2y[n]);
                                                r2 = _mm_mul_pd(r2, r1);
                                                r3 = _mm_mul_pd(r3, r0);
                                                r2 = _mm_add_pd(r2, r3);
                                                _mm_storeu_pd(&fy[n], r2);
                                                r2 = _mm_loadu_pd(&p1z[n]);
                                                r3 = _mm_loadu_pd(&p2z[n]);
                                                r2 = _mm_mul_pd(r2, r1);
                                                r3 = _mm_mul_pd(r3, r0);
                                                r2 = _mm_add_pd(r2, r3);
                                                _mm_storeu_pd(&fz[n], r2);
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
void CINTnabla1k_2e_simd1(double *f, double *g,
                          int li, int lj, int lk, int ll, CINTEnvVars *envs)
{
        int nroots = envs->nrys_roots;
        double ak2 = -2 * envs->ak[0];
        int dk = envs->g_stride_k;
        DEF_GXYZ(double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        if (nroots == 1) { // nabla_k ssss
                assert(li == 0);
                assert(lj == 0);
                assert(lk == 0);
                assert(ll == 0);
                fx[0] = ak2 * gx[dk];
                fy[0] = ak2 * gy[dk];
                fz[0] = ak2 * gz[dk];
        } else {
                int i, j, k, l, n, ptr;
                int di = envs->g_stride_i;
                int dl = envs->g_stride_l;
                int dj = envs->g_stride_j;
                double *p1x = gx - dk;
                double *p1y = gy - dk;
                double *p1z = gz - dk;
                double *p2x = gx + dk;
                double *p2y = gy + dk;
                double *p2z = gz + dk;
                __m128d r0, r1, r2, r3;
                r0 = _mm_load1_pd(&ak2);
                for (j = 0; j <= lj; j++)
                for (l = 0; l <= ll; l++) {
                        ptr = dj * j + dl * l;
                        //f(...,0,...) = -2*ak*g(...,1,...)
                        for (i = 0; i <= li; i++) {
                                for (n = ptr; n < ptr+nroots; n += 2) {
                                        //fx[n] = ak2 * p2x[n];
                                        //fy[n] = ak2 * p2y[n];
                                        //fz[n] = ak2 * p2z[n];
                                        r3 = _mm_loadu_pd(&p2x[n]);
                                        r3 = _mm_mul_pd(r3, r0);
                                        _mm_storeu_pd(&fx[n], r3);
                                        r3 = _mm_loadu_pd(&p2y[n]);
                                        r3 = _mm_mul_pd(r3, r0);
                                        _mm_storeu_pd(&fy[n], r3);
                                        r3 = _mm_loadu_pd(&p2z[n]);
                                        r3 = _mm_mul_pd(r3, r0);
                                        _mm_storeu_pd(&fz[n], r3);
                                }
                                ptr += di;
                        }
                        //f(...,k,...) = k*g(...,k-1,...)-2*ak*g(...,k+1,...)
                        for (k = 1; k <= lk; k++) {
                                r1 = _mm_set1_pd(k);
                                ptr = dj * j + dl * l + dk * k;
                                for (i = 0; i <= li; i++) {
                                        for (n = ptr; n < ptr+nroots; n += 2) {
                                                //fx[n] = k*p1x[n] + ak2*p2x[n];
                                                //fy[n] = k*p1y[n] + ak2*p2y[n];
                                                //fz[n] = k*p1z[n] + ak2*p2z[n];
                                                r2 = _mm_loadu_pd(&p1x[n]);
                                                r3 = _mm_loadu_pd(&p2x[n]);
                                                r2 = _mm_mul_pd(r2, r1);
                                                r3 = _mm_mul_pd(r3, r0);
                                                r2 = _mm_add_pd(r2, r3);
                                                _mm_storeu_pd(&fx[n], r2);
                                                r2 = _mm_loadu_pd(&p1y[n]);
                                                r3 = _mm_loadu_pd(&p2y[n]);
                                                r2 = _mm_mul_pd(r2, r1);
                                                r3 = _mm_mul_pd(r3, r0);
                                                r2 = _mm_add_pd(r2, r3);
                                                _mm_storeu_pd(&fy[n], r2);
                                                r2 = _mm_loadu_pd(&p1z[n]);
                                                r3 = _mm_loadu_pd(&p2z[n]);
                                                r2 = _mm_mul_pd(r2, r1);
                                                r3 = _mm_mul_pd(r3, r0);
                                                r2 = _mm_add_pd(r2, r3);
                                                _mm_storeu_pd(&fz[n], r2);
                                        }
                                        ptr += di;
                                }
                        }
                }
        }
}


/*
 * ( ij | k \nabla l )
 */
void CINTnabla1l_2e_simd1(double *f, double *g,
                          int li, int lj, int lk, int ll, CINTEnvVars *envs)
{
        int nroots = envs->nrys_roots;
        int dl = envs->g_stride_l;
        double al2 = -2 * envs->al[0];
        DEF_GXYZ(double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        if (nroots == 1) { // nabla_l ssss
                assert(li == 0);
                assert(lj == 0);
                assert(lk == 0);
                assert(ll == 0);
                fx[0] = al2 * gx[dl];
                fy[0] = al2 * gy[dl];
                fz[0] = al2 * gz[dl];
        } else {
                int i, j, k, l, n, ptr;
                int di = envs->g_stride_i;
                int dk = envs->g_stride_k;
                int dj = envs->g_stride_j;
                double *p1x = gx - dl;
                double *p1y = gy - dl;
                double *p1z = gz - dl;
                double *p2x = gx + dl;
                double *p2y = gy + dl;
                double *p2z = gz + dl;
                __m128d r0, r1, r2, r3;
                r0 = _mm_load1_pd(&al2);
                for (j = 0; j <= lj; j++) {
                        //f(...,0,...) = -2*al*g(...,1,...)
                        for (k = 0; k <= lk; k++) {
                                ptr = dj * j + dk * k;
                                for (i = 0; i <= li; i++) {
                                        for (n = ptr; n < ptr+nroots; n += 2) {
                                                //fx[n] = al2 * p2x[n];
                                                //fy[n] = al2 * p2y[n];
                                                //fz[n] = al2 * p2z[n];
                                                r3 = _mm_loadu_pd(&p2x[n]);
                                                r3 = _mm_mul_pd(r3, r0);
                                                _mm_storeu_pd(&fx[n], r3);
                                                r3 = _mm_loadu_pd(&p2y[n]);
                                                r3 = _mm_mul_pd(r3, r0);
                                                _mm_storeu_pd(&fy[n], r3);
                                                r3 = _mm_loadu_pd(&p2z[n]);
                                                r3 = _mm_mul_pd(r3, r0);
                                                _mm_storeu_pd(&fz[n], r3);
                                        }
                                        ptr += di;
                                }
                        }
                        //f(...,l,...) = l*g(...,l-1,...)-2*al*g(...,l+1,...)
                        for (l = 1; l <= ll; l++) {
                                r1 = _mm_set1_pd(l);
                                for (k = 0; k <= lk; k++) {
                                        ptr = dj * j + dl * l + dk * k;
                                        for (i = 0; i <= li; i++, ptr += di) {
                                        for (n = ptr; n < ptr+nroots; n += 2) {
                                                //fx[n] = l*p1x[n] + al2*p2x[n];
                                                //fy[n] = l*p1y[n] + al2*p2y[n];
                                                //fz[n] = l*p1z[n] + al2*p2z[n];
                                                r2 = _mm_loadu_pd(&p1x[n]);
                                                r3 = _mm_loadu_pd(&p2x[n]);
                                                r2 = _mm_mul_pd(r2, r1);
                                                r3 = _mm_mul_pd(r3, r0);
                                                r2 = _mm_add_pd(r2, r3);
                                                _mm_storeu_pd(&fx[n], r2);
                                                r2 = _mm_loadu_pd(&p1y[n]);
                                                r3 = _mm_loadu_pd(&p2y[n]);
                                                r2 = _mm_mul_pd(r2, r1);
                                                r3 = _mm_mul_pd(r3, r0);
                                                r2 = _mm_add_pd(r2, r3);
                                                _mm_storeu_pd(&fy[n], r2);
                                                r2 = _mm_loadu_pd(&p1z[n]);
                                                r3 = _mm_loadu_pd(&p2z[n]);
                                                r2 = _mm_mul_pd(r2, r1);
                                                r3 = _mm_mul_pd(r3, r0);
                                                r2 = _mm_add_pd(r2, r3);
                                                _mm_storeu_pd(&fz[n], r2);
                                        } }
                                }
                        }
                }
        }
}

/*
 * ( x^1 i j | kl )
 * ri is the shift from the center R_O to the center of |i>
 * r - R_O = (r-R_i) + ri, ri = R_i - R_O
 */
void CINTx1i_2e_simd1(double *f, double *g, double *ri,
                      int li, int lj, int lk, int ll, CINTEnvVars *envs)
{
        int nroots = envs->nrys_roots;
        int di = envs->g_stride_i;
        DEF_GXYZ(double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        if (nroots == 1) { // x_i ssss
                assert(li == 0);
                assert(lj == 0);
                assert(lk == 0);
                assert(ll == 0);
                fx[0] = gx[di] + ri[0] * gx[0];
                fy[0] = gy[di] + ri[1] * gy[0];
                fz[0] = gz[di] + ri[2] * gz[0];
        } else {
                int i, j, k, l, n, ptr;
                int dk = envs->g_stride_k;
                int dl = envs->g_stride_l;
                int dj = envs->g_stride_j;
                double *p1x = gx + di;
                double *p1y = gy + di;
                double *p1z = gz + di;
                __m128d r0, r1, r2, r3, r4;
                r0 = _mm_load1_pd(&ri[0]);
                r1 = _mm_load1_pd(&ri[1]);
                r2 = _mm_load1_pd(&ri[2]);
                for (j = 0; j <= lj; j++)
                for (l = 0; l <= ll; l++) {
                for (k = 0; k <= lk; k++) {
                        //f(...,0:li,...) = g(...,1:li+1,...) + ri(1)*g(...,0:li,...)
                        ptr = dj * j + dl * l + dk * k;
                        for (i = 0; i <= li; i++) {
                                for (n = ptr; n < ptr+nroots; n += 2) {
                                        //fx[n] = p1x[n] + ri[0] * gx[n];
                                        //fy[n] = p1y[n] + ri[1] * gy[n];
                                        //fz[n] = p1z[n] + ri[2] * gz[n];
                                        r3 = _mm_loadu_pd(&p1x[n]);
                                        r4 = _mm_loadu_pd(&gx [n]);
                                        r4 = _mm_mul_pd(r4, r0);
                                        r3 = _mm_add_pd(r3, r4);
                                        _mm_storeu_pd(&fx[n], r3);
                                        r3 = _mm_loadu_pd(&p1y[n]);
                                        r4 = _mm_loadu_pd(&gy [n]);
                                        r4 = _mm_mul_pd(r4, r1);
                                        r3 = _mm_add_pd(r3, r4);
                                        _mm_storeu_pd(&fy[n], r3);
                                        r3 = _mm_loadu_pd(&p1z[n]);
                                        r4 = _mm_loadu_pd(&gz [n]);
                                        r4 = _mm_mul_pd(r4, r2);
                                        r3 = _mm_add_pd(r3, r4);
                                        _mm_storeu_pd(&fz[n], r3);
                                }
                                ptr += di;
                        }
                } }
        }
}


/*
 * ( i x^1 j | kl )
 */
void CINTx1j_2e_simd1(double *f, double *g, double *rj,
                      int li, int lj, int lk, int ll, CINTEnvVars *envs)
{
        int nroots = envs->nrys_roots;
        int dj = envs->g_stride_j;
        DEF_GXYZ(double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        if (nroots == 1) { // x_j ssss
                assert(li == 0);
                assert(lj == 0);
                assert(lk == 0);
                assert(ll == 0);
                fx[0] = gx[dj] + rj[0] * gx[0];
                fy[0] = gy[dj] + rj[1] * gy[0];
                fz[0] = gz[dj] + rj[2] * gz[0];
        } else {
                int i, j, k, l, n, ptr;
                int di = envs->g_stride_i;
                int dk = envs->g_stride_k;
                int dl = envs->g_stride_l;
                double *p1x = gx + dj;
                double *p1y = gy + dj;
                double *p1z = gz + dj;
                __m128d r0, r1, r2, r3, r4;
                r0 = _mm_load1_pd(&rj[0]);
                r1 = _mm_load1_pd(&rj[1]);
                r2 = _mm_load1_pd(&rj[2]);
                for (j = 0; j <= lj; j++)
                for (l = 0; l <= ll; l++) {
                for (k = 0; k <= lk; k++) {
                        // f(...,0:lj,...) = g(...,1:lj+1,...) + rj(1)*g(...,0:lj,...)
                        ptr = dj * j + dl * l + dk * k;
                        for (i = 0; i <= li; i++) {
                                for (n = ptr; n < ptr+nroots; n += 2) {
                                        //fx[n] = p1x[n] + rj[0] * gx[n];
                                        //fy[n] = p1y[n] + rj[1] * gy[n];
                                        //fz[n] = p1z[n] + rj[2] * gz[n];
                                        r3 = _mm_loadu_pd(&p1x[n]);
                                        r4 = _mm_loadu_pd(&gx [n]);
                                        r4 = _mm_mul_pd(r4, r0);
                                        r3 = _mm_add_pd(r3, r4);
                                        _mm_storeu_pd(&fx[n], r3);
                                        r3 = _mm_loadu_pd(&p1y[n]);
                                        r4 = _mm_loadu_pd(&gy [n]);
                                        r4 = _mm_mul_pd(r4, r1);
                                        r3 = _mm_add_pd(r3, r4);
                                        _mm_storeu_pd(&fy[n], r3);
                                        r3 = _mm_loadu_pd(&p1z[n]);
                                        r4 = _mm_loadu_pd(&gz [n]);
                                        r4 = _mm_mul_pd(r4, r2);
                                        r3 = _mm_add_pd(r3, r4);
                                        _mm_storeu_pd(&fz[n], r3);
                                }
                                ptr += di;
                        }
                } }
        }
}


/*
 * ( ij | x^1 k l )
 */
void CINTx1k_2e_simd1(double *f, double *g, double *rk,
                      int li, int lj, int lk, int ll, CINTEnvVars *envs)
{
        int nroots = envs->nrys_roots;
        int dk = envs->g_stride_k;
        DEF_GXYZ(double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        if (nroots == 1) { // x_k ssss
                assert(li == 0);
                assert(lj == 0);
                assert(lk == 0);
                assert(ll == 0);
                fx[0] = gx[dk] + rk[0] * gx[0];
                fy[0] = gy[dk] + rk[1] * gy[0];
                fz[0] = gz[dk] + rk[2] * gz[0];
        } else {
                int i, j, k, l, n, ptr;
                int di = envs->g_stride_i;
                int dl = envs->g_stride_l;
                int dj = envs->g_stride_j;
                double *p1x = gx + dk;
                double *p1y = gy + dk;
                double *p1z = gz + dk;
                __m128d r0, r1, r2, r3, r4;
                r0 = _mm_load1_pd(&rk[0]);
                r1 = _mm_load1_pd(&rk[1]);
                r2 = _mm_load1_pd(&rk[2]);
                for (j = 0; j <= lj; j++)
                for (l = 0; l <= ll; l++) {
                for (k = 0; k <= lk; k++) {
                        // f(...,0:lk,...) = g(...,1:lk+1,...) + rk(1)*g(...,0:lk,...)
                        ptr = dj * j + dl * l + dk * k;
                        for (i = 0; i <= li; i++) {
                                for (n = ptr; n < ptr+nroots; n += 2) {
                                        //fx[n] = p1x[n] + rk[0] * gx[n];
                                        //fy[n] = p1y[n] + rk[1] * gy[n];
                                        //fz[n] = p1z[n] + rk[2] * gz[n];
                                        r3 = _mm_loadu_pd(&p1x[n]);
                                        r4 = _mm_loadu_pd(&gx [n]);
                                        r4 = _mm_mul_pd(r4, r0);
                                        r3 = _mm_add_pd(r3, r4);
                                        _mm_storeu_pd(&fx[n], r3);
                                        r3 = _mm_loadu_pd(&p1y[n]);
                                        r4 = _mm_loadu_pd(&gy [n]);
                                        r4 = _mm_mul_pd(r4, r1);
                                        r3 = _mm_add_pd(r3, r4);
                                        _mm_storeu_pd(&fy[n], r3);
                                        r3 = _mm_loadu_pd(&p1z[n]);
                                        r4 = _mm_loadu_pd(&gz [n]);
                                        r4 = _mm_mul_pd(r4, r2);
                                        r3 = _mm_add_pd(r3, r4);
                                        _mm_storeu_pd(&fz[n], r3);
                                }
                                ptr += di;
                        }
                } }
        }
}


/*
 * ( i j | x^1 kl )
 */
void CINTx1l_2e_simd1(double *f, double *g, double *rl,
                      int li, int lj, int lk, int ll, CINTEnvVars *envs)
{
        int nroots = envs->nrys_roots;
        int dl = envs->g_stride_l;
        DEF_GXYZ(double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        if (nroots == 1) { // x_l ssss
                assert(li == 0);
                assert(lj == 0);
                assert(lk == 0);
                assert(ll == 0);
                fx[0] = gx[dl] + rl[0] * gx[0];
                fy[0] = gy[dl] + rl[1] * gy[0];
                fz[0] = gz[dl] + rl[2] * gz[0];
        } else {
                int i, j, k, l, n, ptr;
                int di = envs->g_stride_i;
                int dk = envs->g_stride_k;
                int dj = envs->g_stride_j;
                double *p1x = gx + dl;
                double *p1y = gy + dl;
                double *p1z = gz + dl;
                __m128d r0, r1, r2, r3, r4;
                r0 = _mm_load1_pd(&rl[0]);
                r1 = _mm_load1_pd(&rl[1]);
                r2 = _mm_load1_pd(&rl[2]);
                for (j = 0; j <= lj; j++)
                for (l = 0; l <= ll; l++) {
                for (k = 0; k <= lk; k++) {
                        // f(...,0:ll,...) = g(...,1:ll+1,...) + rl(1)*g(...,0:ll,...)
                        ptr = dj * j + dl * l + dk * k;
                        for (i = 0; i <= li; i++) {
                                for (n = ptr; n < ptr+nroots; n += 2) {
                                        //fx[n] = p1x[n] + rl[0] * gx[n];
                                        //fy[n] = p1y[n] + rl[1] * gy[n];
                                        //fz[n] = p1z[n] + rl[2] * gz[n];
                                        r3 = _mm_loadu_pd(&p1x[n]);
                                        r4 = _mm_loadu_pd(&gx [n]);
                                        r4 = _mm_mul_pd(r4, r0);
                                        r3 = _mm_add_pd(r3, r4);
                                        _mm_storeu_pd(&fx[n], r3);
                                        r3 = _mm_loadu_pd(&p1y[n]);
                                        r4 = _mm_loadu_pd(&gy [n]);
                                        r4 = _mm_mul_pd(r4, r1);
                                        r3 = _mm_add_pd(r3, r4);
                                        _mm_storeu_pd(&fy[n], r3);
                                        r3 = _mm_loadu_pd(&p1z[n]);
                                        r4 = _mm_loadu_pd(&gz [n]);
                                        r4 = _mm_mul_pd(r4, r2);
                                        r3 = _mm_add_pd(r3, r4);
                                        _mm_storeu_pd(&fz[n], r3);
                                }
                                ptr += di;
                        }
                } }
        }
}


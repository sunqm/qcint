/*
 * Qcint is a general GTO integral library for computational chemistry
 * Copyright (C) 2014 Qiming Sun
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


#include <pmmintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "cint_const.h"
#include "cint_bas.h"
#include "misc.h"
#include "g2e.h"

#if defined(__GNUC__)
#define ALIGN16 __attribute__((aligned(16)))
#define RESTRICT __restrict__
#else
#define ALIGN16
#define RESTRICT
#endif

#define DEF_GXYZ(type, G, GX, GY, GZ) \
        type *GX = G; \
        type *GY = G + envs->g_size; \
        type *GZ = G + envs->g_size * 2

struct _BC {
        double c00[MXRYSROOTS*3];
        double c0p[MXRYSROOTS*3];
        double b01[MXRYSROOTS];
        double b00[MXRYSROOTS];
        double b10[MXRYSROOTS];
};

void CINTg0_2e_2d(double *g, struct _BC *bc, const CINTEnvVars *envs);

/*
 * g0[i,k,l,j] = < ik | lj > = ( i j | k l )
 */
/* 2d is based on l,j */
void CINTg0_lj2d_4d(double *g, const CINTEnvVars *envs)
{
        const int nmax = envs->li_ceil + envs->lj_ceil;
        const int mmax = envs->lk_ceil + envs->ll_ceil;
        const int li = envs->li_ceil;
        const int lk = envs->lk_ceil;
        //const int ll = envs->ll_ceil;
        const int lj = envs->lj_ceil;
        const int nroots = envs->nrys_roots;
        int i, j, k, l, ptr, n;
        const int di = envs->g_stride_i;
        const int dk = envs->g_stride_k;
        const int dl = envs->g_stride_l;
        const int dj = envs->g_stride_j;
        const double *rirj = envs->rirj;
        const double *rkrl = envs->rkrl;
        DEF_GXYZ(double, g, gx, gy, gz);
        const double *p1x, *p1y, *p1z, *p2x, *p2y, *p2z;

        // g(i,...,j) = rirj * g(i-1,...,j) +  g(i-1,...,j+1)
        p1x = gx - di;
        p1y = gy - di;
        p1z = gz - di;
        p2x = gx - di + dj;
        p2y = gy - di + dj;
        p2z = gz - di + dj;
        __m128d r0, r1, r2, r3, r4;
        r0 = _mm_load1_pd(&rirj[0]);
        r1 = _mm_load1_pd(&rirj[1]);
        r2 = _mm_load1_pd(&rirj[2]);
        for (i = 1; i <= li; i++) {
        for (j = 0; j <= nmax-i; j++) {
        for (l = 0; l <= mmax; l++) {
                ptr = j*dj + l*dl + i*di;
                for (n = ptr; n < ptr+nroots; n += 2) {
                        //gx[n] = rirj[0] * p1x[n] + p2x[n];
                        //gy[n] = rirj[1] * p1y[n] + p2y[n];
                        //gz[n] = rirj[2] * p1z[n] + p2z[n];
                        r3 = _mm_load_pd(&p1x[n]);
                        r4 = _mm_load_pd(&p2x[n]);
                        r3 = _mm_mul_pd(r3, r0);
                        r3 = _mm_add_pd(r3, r4);
                        _mm_store_pd(&gx[n], r3);
                        r3 = _mm_load_pd(&p1y[n]);
                        r4 = _mm_load_pd(&p2y[n]);
                        r3 = _mm_mul_pd(r3, r1);
                        r3 = _mm_add_pd(r3, r4);
                        _mm_store_pd(&gy[n], r3);
                        r3 = _mm_load_pd(&p1z[n]);
                        r4 = _mm_load_pd(&p2z[n]);
                        r3 = _mm_mul_pd(r3, r2);
                        r3 = _mm_add_pd(r3, r4);
                        _mm_store_pd(&gz[n], r3);
                }
        } } }

        // g(...,k,l,..) = rkrl * g(...,k-1,l,..) + g(...,k-1,l+1,..)
        p1x = gx - dk;
        p1y = gy - dk;
        p1z = gz - dk;
        p2x = gx - dk + dl;
        p2y = gy - dk + dl;
        p2z = gz - dk + dl;
        r0 = _mm_load1_pd(&rkrl[0]);
        r1 = _mm_load1_pd(&rkrl[1]);
        r2 = _mm_load1_pd(&rkrl[2]);
        for (j = 0; j <= lj; j++) {
        for (k = 1; k <= lk; k++) {
        for (l = 0; l <= mmax-k; l++) {
                ptr = j*dj + l*dl + k*dk;
                for (n = ptr; n < ptr+dk; n += 2) {
                        //gx[n] = rkrl[0] * p1x[n] + p2x[n];
                        //gy[n] = rkrl[1] * p1y[n] + p2y[n];
                        //gz[n] = rkrl[2] * p1z[n] + p2z[n];
                        r3 = _mm_load_pd(&p1x[n]);
                        r4 = _mm_load_pd(&p2x[n]);
                        r3 = _mm_mul_pd(r3, r0);
                        r3 = _mm_add_pd(r3, r4);
                        _mm_store_pd(&gx[n], r3);
                        r3 = _mm_load_pd(&p1y[n]);
                        r4 = _mm_load_pd(&p2y[n]);
                        r3 = _mm_mul_pd(r3, r1);
                        r3 = _mm_add_pd(r3, r4);
                        _mm_store_pd(&gy[n], r3);
                        r3 = _mm_load_pd(&p1z[n]);
                        r4 = _mm_load_pd(&p2z[n]);
                        r3 = _mm_mul_pd(r3, r2);
                        r3 = _mm_add_pd(r3, r4);
                        _mm_store_pd(&gz[n], r3);
                }
        } } }
}
/* 2d is based on k,j */
static void CINTg0_kj2d_4d(double *g, const CINTEnvVars *envs)
{
        const int nmax = envs->li_ceil + envs->lj_ceil;
        const int mmax = envs->lk_ceil + envs->ll_ceil;
        const int li = envs->li_ceil;
        //const int lk = envs->lk_ceil;
        const int ll = envs->ll_ceil;
        const int lj = envs->lj_ceil;
        const int nroots = envs->nrys_roots;
        int i, j, k, l, ptr, n;
        const int di = envs->g_stride_i;
        const int dk = envs->g_stride_k;
        const int dl = envs->g_stride_l;
        const int dj = envs->g_stride_j;
        const double *rirj = envs->rirj;
        const double *rkrl = envs->rkrl;
        DEF_GXYZ(double, g, gx, gy, gz);
        const double *p1x, *p1y, *p1z, *p2x, *p2y, *p2z;

        // g(i,...,j) = rirj * g(i-1,...,j) +  g(i-1,...,j+1)
        p1x = gx - di;
        p1y = gy - di;
        p1z = gz - di;
        p2x = gx - di + dj;
        p2y = gy - di + dj;
        p2z = gz - di + dj;
        __m128d r0, r1, r2, r3, r4;
        r0 = _mm_load1_pd(&rirj[0]);
        r1 = _mm_load1_pd(&rirj[1]);
        r2 = _mm_load1_pd(&rirj[2]);
        for (i = 1; i <= li; i++) {
        for (j = 0; j <= nmax-i; j++) {
        for (k = 0; k <= mmax; k++) {
                ptr = j*dj + k*dk + i*di;
                for (n = ptr; n < ptr+nroots; n += 2) {
                        //gx[n] = rirj[0] * p1x[n] + p2x[n];
                        //gy[n] = rirj[1] * p1y[n] + p2y[n];
                        //gz[n] = rirj[2] * p1z[n] + p2z[n];
                        r3 = _mm_load_pd(&p1x[n]);
                        r4 = _mm_load_pd(&p2x[n]);
                        r3 = _mm_mul_pd(r3, r0);
                        r3 = _mm_add_pd(r3, r4);
                        _mm_store_pd(&gx[n], r3);
                        r3 = _mm_load_pd(&p1y[n]);
                        r4 = _mm_load_pd(&p2y[n]);
                        r3 = _mm_mul_pd(r3, r1);
                        r3 = _mm_add_pd(r3, r4);
                        _mm_store_pd(&gy[n], r3);
                        r3 = _mm_load_pd(&p1z[n]);
                        r4 = _mm_load_pd(&p2z[n]);
                        r3 = _mm_mul_pd(r3, r2);
                        r3 = _mm_add_pd(r3, r4);
                        _mm_store_pd(&gz[n], r3);
                }
        } } }

        // g(...,k,l,..) = rkrl * g(...,k,l-1,..) + g(...,k+1,l-1,..)
        p1x = gx - dl;
        p1y = gy - dl;
        p1z = gz - dl;
        p2x = gx - dl + dk;
        p2y = gy - dl + dk;
        p2z = gz - dl + dk;
        r0 = _mm_load1_pd(&rkrl[0]);
        r1 = _mm_load1_pd(&rkrl[1]);
        r2 = _mm_load1_pd(&rkrl[2]);
        for (j = 0; j <= lj; j++) {
        for (l = 1; l <= ll; l++) {
        for (k = 0; k <= mmax-l; k++) {
                ptr = j*dj + l*dl + k*dk;
                for (n = ptr; n < ptr+dk; n += 2) {
                        //gx[n] = rkrl[0] * p1x[n] + p2x[n];
                        //gy[n] = rkrl[1] * p1y[n] + p2y[n];
                        //gz[n] = rkrl[2] * p1z[n] + p2z[n];
                        r3 = _mm_load_pd(&p1x[n]);
                        r4 = _mm_load_pd(&p2x[n]);
                        r3 = _mm_mul_pd(r3, r0);
                        r3 = _mm_add_pd(r3, r4);
                        _mm_store_pd(&gx[n], r3);
                        r3 = _mm_load_pd(&p1y[n]);
                        r4 = _mm_load_pd(&p2y[n]);
                        r3 = _mm_mul_pd(r3, r1);
                        r3 = _mm_add_pd(r3, r4);
                        _mm_store_pd(&gy[n], r3);
                        r3 = _mm_load_pd(&p1z[n]);
                        r4 = _mm_load_pd(&p2z[n]);
                        r3 = _mm_mul_pd(r3, r2);
                        r3 = _mm_add_pd(r3, r4);
                        _mm_store_pd(&gz[n], r3);
                }
        } } }
}
/* 2d is based on i,l */
static void CINTg0_il2d_4d(double *g, const CINTEnvVars *envs)
{
        const int nmax = envs->li_ceil + envs->lj_ceil;
        const int mmax = envs->lk_ceil + envs->ll_ceil;
        //const int li = envs->li_ceil;
        const int lk = envs->lk_ceil;
        const int ll = envs->ll_ceil;
        const int lj = envs->lj_ceil;
        const int nroots = envs->nrys_roots;
        int i, j, k, l, ptr, n;
        const int di = envs->g_stride_i;
        const int dk = envs->g_stride_k;
        const int dl = envs->g_stride_l;
        const int dj = envs->g_stride_j;
        const double *rirj = envs->rirj;
        const double *rkrl = envs->rkrl;
        DEF_GXYZ(double, g, gx, gy, gz);
        const double *p1x, *p1y, *p1z, *p2x, *p2y, *p2z;

        // g(...,k,l,..) = rkrl * g(...,k-1,l,..) + g(...,k-1,l+1,..)
        p1x = gx - dk;
        p1y = gy - dk;
        p1z = gz - dk;
        p2x = gx - dk + dl;
        p2y = gy - dk + dl;
        p2z = gz - dk + dl;
        __m128d r0, r1, r2, r3, r4;
        r0 = _mm_load1_pd(&rkrl[0]);
        r1 = _mm_load1_pd(&rkrl[1]);
        r2 = _mm_load1_pd(&rkrl[2]);
        for (k = 1; k <= lk; k++) {
        for (l = 0; l <= mmax-k; l++) {
        for (i = 0; i <= nmax; i++) {
                ptr = l*dl + k*dk + i*di;
                for (n = ptr; n < ptr+nroots; n += 2) {
                        //gx[n] = rkrl[0] * p1x[n] + p2x[n];
                        //gy[n] = rkrl[1] * p1y[n] + p2y[n];
                        //gz[n] = rkrl[2] * p1z[n] + p2z[n];
                        r3 = _mm_load_pd(&p1x[n]);
                        r4 = _mm_load_pd(&p2x[n]);
                        r3 = _mm_mul_pd(r3, r0);
                        r3 = _mm_add_pd(r3, r4);
                        _mm_store_pd(&gx[n], r3);
                        r3 = _mm_load_pd(&p1y[n]);
                        r4 = _mm_load_pd(&p2y[n]);
                        r3 = _mm_mul_pd(r3, r1);
                        r3 = _mm_add_pd(r3, r4);
                        _mm_store_pd(&gy[n], r3);
                        r3 = _mm_load_pd(&p1z[n]);
                        r4 = _mm_load_pd(&p2z[n]);
                        r3 = _mm_mul_pd(r3, r2);
                        r3 = _mm_add_pd(r3, r4);
                        _mm_store_pd(&gz[n], r3);
                }
        } } }

        // g(i,...,j) = rirj * g(i,...,j-1) +  g(i+1,...,j-1)
        p1x = gx - dj;
        p1y = gy - dj;
        p1z = gz - dj;
        p2x = gx - dj + di;
        p2y = gy - dj + di;
        p2z = gz - dj + di;
        r0 = _mm_load1_pd(&rirj[0]);
        r1 = _mm_load1_pd(&rirj[1]);
        r2 = _mm_load1_pd(&rirj[2]);
        for (j = 1; j <= lj; j++) {
        for (l = 0; l <= ll; l++) {
        for (k = 0; k <= lk; k++) {
                ptr = j*dj + l*dl + k*dk;
                for (n = ptr; n < ptr+dk-di*j; n += 2) {
                        //gx[n] = rirj[0] * p1x[n] + p2x[n];
                        //gy[n] = rirj[1] * p1y[n] + p2y[n];
                        //gz[n] = rirj[2] * p1z[n] + p2z[n];
                        r3 = _mm_load_pd(&p1x[n]);
                        r4 = _mm_load_pd(&p2x[n]);
                        r3 = _mm_mul_pd(r3, r0);
                        r3 = _mm_add_pd(r3, r4);
                        _mm_store_pd(&gx[n], r3);
                        r3 = _mm_load_pd(&p1y[n]);
                        r4 = _mm_load_pd(&p2y[n]);
                        r3 = _mm_mul_pd(r3, r1);
                        r3 = _mm_add_pd(r3, r4);
                        _mm_store_pd(&gy[n], r3);
                        r3 = _mm_load_pd(&p1z[n]);
                        r4 = _mm_load_pd(&p2z[n]);
                        r3 = _mm_mul_pd(r3, r2);
                        r3 = _mm_add_pd(r3, r4);
                        _mm_store_pd(&gz[n], r3);
                }
        } } }
}
/* 2d is based on i,k */
static void CINTg0_ik2d_4d(double *g, const CINTEnvVars *envs)
{
        const int nmax = envs->li_ceil + envs->lj_ceil;
        const int mmax = envs->lk_ceil + envs->ll_ceil;
        //const int li = envs->li_ceil;
        const int lk = envs->lk_ceil;
        const int ll = envs->ll_ceil;
        const int lj = envs->lj_ceil;
        const int nroots = envs->nrys_roots;
        int i, j, k, l, ptr, n;
        const int di = envs->g_stride_i;
        const int dk = envs->g_stride_k;
        const int dl = envs->g_stride_l;
        const int dj = envs->g_stride_j;
        const double *rirj = envs->rirj;
        const double *rkrl = envs->rkrl;
        DEF_GXYZ(double, g, gx, gy, gz);
        const double *p1x, *p1y, *p1z, *p2x, *p2y, *p2z;

        // g(...,k,l,..) = rkrl * g(...,k,l-1,..) + g(...,k+1,l-1,..)
        p1x = gx - dl;
        p1y = gy - dl;
        p1z = gz - dl;
        p2x = gx - dl + dk;
        p2y = gy - dl + dk;
        p2z = gz - dl + dk;
        __m128d r0, r1, r2, r3, r4;
        r0 = _mm_load1_pd(&rkrl[0]);
        r1 = _mm_load1_pd(&rkrl[1]);
        r2 = _mm_load1_pd(&rkrl[2]);
        for (l = 1; l <= ll; l++) {
                // (:,i) is full, so loop:k and loop:n can be merged to
                // for(n = l*dl; n < ptr+dl-dk*l; n++)
                for (k = 0; k <= mmax-l; k++) {
                for (i = 0; i <= nmax; i++) {
                        ptr = l*dl + k*dk + i*di;
                        for (n = ptr; n < ptr+nroots; n += 2) {
                                //gx[n] = rkrl[0] * p1x[n] + p2x[n];
                                //gy[n] = rkrl[1] * p1y[n] + p2y[n];
                                //gz[n] = rkrl[2] * p1z[n] + p2z[n];
                                r3 = _mm_load_pd(&p1x[n]);
                                r4 = _mm_load_pd(&p2x[n]);
                                r3 = _mm_mul_pd(r3, r0);
                                r3 = _mm_add_pd(r3, r4);
                                _mm_store_pd(&gx[n], r3);
                                r3 = _mm_load_pd(&p1y[n]);
                                r4 = _mm_load_pd(&p2y[n]);
                                r3 = _mm_mul_pd(r3, r1);
                                r3 = _mm_add_pd(r3, r4);
                                _mm_store_pd(&gy[n], r3);
                                r3 = _mm_load_pd(&p1z[n]);
                                r4 = _mm_load_pd(&p2z[n]);
                                r3 = _mm_mul_pd(r3, r2);
                                r3 = _mm_add_pd(r3, r4);
                                _mm_store_pd(&gz[n], r3);
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
        r0 = _mm_load1_pd(&rirj[0]);
        r1 = _mm_load1_pd(&rirj[1]);
        r2 = _mm_load1_pd(&rirj[2]);
        for (j = 1; j <= lj; j++) {
        for (l = 0; l <= ll; l++) {
        for (k = 0; k <= lk; k++) {
                ptr = j*dj + l*dl + k*dk;
                for (n = ptr; n < ptr+dk-di*j; n += 2) {
                        //gx[n] = rirj[0] * p1x[n] + p2x[n];
                        //gy[n] = rirj[1] * p1y[n] + p2y[n];
                        //gz[n] = rirj[2] * p1z[n] + p2z[n];
                        r3 = _mm_load_pd(&p1x[n]);
                        r4 = _mm_load_pd(&p2x[n]);
                        r3 = _mm_mul_pd(r3, r0);
                        r3 = _mm_add_pd(r3, r4);
                        _mm_store_pd(&gx[n], r3);
                        r3 = _mm_load_pd(&p1y[n]);
                        r4 = _mm_load_pd(&p2y[n]);
                        r3 = _mm_mul_pd(r3, r1);
                        r3 = _mm_add_pd(r3, r4);
                        _mm_store_pd(&gy[n], r3);
                        r3 = _mm_load_pd(&p1z[n]);
                        r4 = _mm_load_pd(&p2z[n]);
                        r3 = _mm_mul_pd(r3, r2);
                        r3 = _mm_add_pd(r3, r4);
                        _mm_store_pd(&gz[n], r3);
                }
        } } }
}

void CINTg0_2e_kj2d4d(double *g, const CINTEnvVars *envs,struct _BC *bc)
{
        CINTg0_2e_2d(g, bc, envs);
        CINTg0_kj2d_4d(g, envs);
}
void CINTg0_2e_ik2d4d(double *g, const CINTEnvVars *envs,struct _BC *bc)
{
        CINTg0_2e_2d(g, bc, envs);
        CINTg0_ik2d_4d(g, envs);
}
void CINTg0_2e_il2d4d(double *g, const CINTEnvVars *envs,struct _BC *bc)
{
        CINTg0_2e_2d(g, bc, envs);
        CINTg0_il2d_4d(g, envs);
}

/*
 * g[i,k,l,j] = < ik | lj > = ( i j | k l )
 */
void CINTg0_2e(double *g, const double fac, const CINTEnvVars *envs)
{
        const double aij = envs->aij;
        const double akl = envs->akl;
        double a0, a1, fac1, x;
        double u[MXRYSROOTS];
        double *w = g + envs->g_size * 2; // ~ gz
        double rijrkl[3];
        rijrkl[0] = envs->rij[0] - envs->rkl[0];
        rijrkl[1] = envs->rij[1] - envs->rkl[1];
        rijrkl[2] = envs->rij[2] - envs->rkl[2];

        a1 = aij * akl;
        a0 = a1 / (aij + akl);
        fac1 = sqrt(a0 / (a1 * a1 * a1)) * fac;
        x = a0 *(rijrkl[0] * rijrkl[0]
               + rijrkl[1] * rijrkl[1]
               + rijrkl[2] * rijrkl[2]);
        CINTrys_roots(envs->nrys_roots, x, u, w);

        int irys;
        double u2, div, tmp1, tmp2, tmp3, tmp4;
        const double *RESTRICT rijrx = envs->rijrx;
        const double *RESTRICT rklrx = envs->rklrx;
//FIXME: I assumed the 16-byte alignment of stack by
// gcc -mpreferred-stack-boundary=4.
// Maybe not, maybe depend on compiler/platform (old gcc does not support?).
// if not, change _mm_load_pd(c00/c0p) to _mm_loadu_pd(c00/c0p)
        struct _BC bc ALIGN16;
        double *RESTRICT b00 = bc.b00;
        double *RESTRICT b10 = bc.b10;
        double *RESTRICT b01 = bc.b01;
        double *RESTRICT c00x = bc.c00;
        double *RESTRICT c00y = bc.c00+MXRYSROOTS;
        double *RESTRICT c00z = bc.c00+MXRYSROOTS*2;
        double *RESTRICT c0px = bc.c0p;
        double *RESTRICT c0py = bc.c0p+MXRYSROOTS;
        double *RESTRICT c0pz = bc.c0p+MXRYSROOTS*2;

        for (irys = 0; irys < envs->nrys_roots; irys++)
        {
                /*
                 *t2 = u(irys)/(1+u(irys))
                 *u2 = aij*akl/(aij+akl)*t2/(1-t2)
                 */
                u2 = a0 * u[irys];
                div = 1 / (u2 * (aij + akl) + a1);
                tmp1 = u2 * div;
                tmp2 = tmp1 * akl;
                tmp3 = tmp1 * aij;
                tmp4 = .5 * div;
                b00[irys] = 0.5 * tmp1;
                b10[irys] = b00[irys] + tmp4 * akl;
                b01[irys] = b00[irys] + tmp4 * aij;
                c00x[irys] = rijrx[0] - tmp2 * rijrkl[0];
                c00y[irys] = rijrx[1] - tmp2 * rijrkl[1];
                c00z[irys] = rijrx[2] - tmp2 * rijrkl[2];
                c0px[irys] = rklrx[0] + tmp3 * rijrkl[0];
                c0py[irys] = rklrx[1] + tmp3 * rijrkl[1];
                c0pz[irys] = rklrx[2] + tmp3 * rijrkl[2];
                w[irys] *= fac1;
        }

        (*envs->f_g0_2d4d)(g, envs, &bc);
}


/*
 * ( \nabla i j | kl )
 */
void CINTnabla1i_2e(double *f, const double *g,
                    const int li, const int lj, const int lk, const int ll,
                    const CINTEnvVars *envs)
{
        const int nroots = envs->nrys_roots;
        const int di = envs->g_stride_i;
        const double ai2 = -2 * envs->ai;
        DEF_GXYZ(const double, g, gx, gy, gz);
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
                const int dk = envs->g_stride_k;
                const int dl = envs->g_stride_l;
                const int dj = envs->g_stride_j;
                const double *p1x = gx - di;
                const double *p1y = gy - di;
                const double *p1z = gz - di;
                const double *p2x = gx + di;
                const double *p2y = gy + di;
                const double *p2z = gz + di;
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
                                r3 = _mm_load_pd(&p2x[n]);
                                r3 = _mm_mul_pd(r3, r0);
                                _mm_store_pd(&fx[n], r3);
                                r3 = _mm_load_pd(&p2y[n]);
                                r3 = _mm_mul_pd(r3, r0);
                                _mm_store_pd(&fy[n], r3);
                                r3 = _mm_load_pd(&p2z[n]);
                                r3 = _mm_mul_pd(r3, r0);
                                _mm_store_pd(&fz[n], r3);
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
                                        r2 = _mm_load_pd(&p1x[n]);
                                        r3 = _mm_load_pd(&p2x[n]);
                                        r2 = _mm_mul_pd(r2, r1);
                                        r3 = _mm_mul_pd(r3, r0);
                                        r2 = _mm_add_pd(r2, r3);
                                        _mm_store_pd(&fx[n], r2);
                                        r2 = _mm_load_pd(&p1y[n]);
                                        r3 = _mm_load_pd(&p2y[n]);
                                        r2 = _mm_mul_pd(r2, r1);
                                        r3 = _mm_mul_pd(r3, r0);
                                        r2 = _mm_add_pd(r2, r3);
                                        _mm_store_pd(&fy[n], r2);
                                        r2 = _mm_load_pd(&p1z[n]);
                                        r3 = _mm_load_pd(&p2z[n]);
                                        r2 = _mm_mul_pd(r2, r1);
                                        r3 = _mm_mul_pd(r3, r0);
                                        r2 = _mm_add_pd(r2, r3);
                                        _mm_store_pd(&fz[n], r2);
                                }
                                ptr += di;
                        }
                }
        }
}


/*
 * ( i \nabla j | kl )
 */
void CINTnabla1j_2e(double *f, const double *g,
                    const int li, const int lj, const int lk, const int ll,
                    const CINTEnvVars *envs)
{
        const int nroots = envs->nrys_roots;
        const int dj = envs->g_stride_j;
        const double aj2 = -2 * envs->aj;
        DEF_GXYZ(const double, g, gx, gy, gz);
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
                const int di = envs->g_stride_i;
                const int dk = envs->g_stride_k;
                const int dl = envs->g_stride_l;
                const double *p1x = gx - dj;
                const double *p1y = gy - dj;
                const double *p1z = gz - dj;
                const double *p2x = gx + dj;
                const double *p2y = gy + dj;
                const double *p2z = gz + dj;
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
                                        r3 = _mm_load_pd(&p2x[n]);
                                        r3 = _mm_mul_pd(r3, r0);
                                        _mm_store_pd(&fx[n], r3);
                                        r3 = _mm_load_pd(&p2y[n]);
                                        r3 = _mm_mul_pd(r3, r0);
                                        _mm_store_pd(&fy[n], r3);
                                        r3 = _mm_load_pd(&p2z[n]);
                                        r3 = _mm_mul_pd(r3, r0);
                                        _mm_store_pd(&fz[n], r3);
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
                                                r2 = _mm_load_pd(&p1x[n]);
                                                r3 = _mm_load_pd(&p2x[n]);
                                                r2 = _mm_mul_pd(r2, r1);
                                                r3 = _mm_mul_pd(r3, r0);
                                                r2 = _mm_add_pd(r2, r3);
                                                _mm_store_pd(&fx[n], r2);
                                                r2 = _mm_load_pd(&p1y[n]);
                                                r3 = _mm_load_pd(&p2y[n]);
                                                r2 = _mm_mul_pd(r2, r1);
                                                r3 = _mm_mul_pd(r3, r0);
                                                r2 = _mm_add_pd(r2, r3);
                                                _mm_store_pd(&fy[n], r2);
                                                r2 = _mm_load_pd(&p1z[n]);
                                                r3 = _mm_load_pd(&p2z[n]);
                                                r2 = _mm_mul_pd(r2, r1);
                                                r3 = _mm_mul_pd(r3, r0);
                                                r2 = _mm_add_pd(r2, r3);
                                                _mm_store_pd(&fz[n], r2);
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
void CINTnabla1k_2e(double *f, const double *g,
                    const int li, const int lj, const int lk, const int ll,
                    const CINTEnvVars *envs)
{
        const int nroots = envs->nrys_roots;
        const double ak2 = -2 * envs->ak;
        const int dk = envs->g_stride_k;
        DEF_GXYZ(const double, g, gx, gy, gz);
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
                const int di = envs->g_stride_i;
                const int dl = envs->g_stride_l;
                const int dj = envs->g_stride_j;
                const double *p1x = gx - dk;
                const double *p1y = gy - dk;
                const double *p1z = gz - dk;
                const double *p2x = gx + dk;
                const double *p2y = gy + dk;
                const double *p2z = gz + dk;
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
                                        r3 = _mm_load_pd(&p2x[n]);
                                        r3 = _mm_mul_pd(r3, r0);
                                        _mm_store_pd(&fx[n], r3);
                                        r3 = _mm_load_pd(&p2y[n]);
                                        r3 = _mm_mul_pd(r3, r0);
                                        _mm_store_pd(&fy[n], r3);
                                        r3 = _mm_load_pd(&p2z[n]);
                                        r3 = _mm_mul_pd(r3, r0);
                                        _mm_store_pd(&fz[n], r3);
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
                                                r2 = _mm_load_pd(&p1x[n]);
                                                r3 = _mm_load_pd(&p2x[n]);
                                                r2 = _mm_mul_pd(r2, r1);
                                                r3 = _mm_mul_pd(r3, r0);
                                                r2 = _mm_add_pd(r2, r3);
                                                _mm_store_pd(&fx[n], r2);
                                                r2 = _mm_load_pd(&p1y[n]);
                                                r3 = _mm_load_pd(&p2y[n]);
                                                r2 = _mm_mul_pd(r2, r1);
                                                r3 = _mm_mul_pd(r3, r0);
                                                r2 = _mm_add_pd(r2, r3);
                                                _mm_store_pd(&fy[n], r2);
                                                r2 = _mm_load_pd(&p1z[n]);
                                                r3 = _mm_load_pd(&p2z[n]);
                                                r2 = _mm_mul_pd(r2, r1);
                                                r3 = _mm_mul_pd(r3, r0);
                                                r2 = _mm_add_pd(r2, r3);
                                                _mm_store_pd(&fz[n], r2);
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
void CINTnabla1l_2e(double *f, const double *g,
                    const int li, const int lj, const int lk, const int ll,
                    const CINTEnvVars *envs)
{
        const int nroots = envs->nrys_roots;
        const int dl = envs->g_stride_l;
        const double al2 = -2 * envs->al;
        DEF_GXYZ(const double, g, gx, gy, gz);
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
                const int di = envs->g_stride_i;
                const int dk = envs->g_stride_k;
                const int dj = envs->g_stride_j;
                const double *p1x = gx - dl;
                const double *p1y = gy - dl;
                const double *p1z = gz - dl;
                const double *p2x = gx + dl;
                const double *p2y = gy + dl;
                const double *p2z = gz + dl;
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
                                                r3 = _mm_load_pd(&p2x[n]);
                                                r3 = _mm_mul_pd(r3, r0);
                                                _mm_store_pd(&fx[n], r3);
                                                r3 = _mm_load_pd(&p2y[n]);
                                                r3 = _mm_mul_pd(r3, r0);
                                                _mm_store_pd(&fy[n], r3);
                                                r3 = _mm_load_pd(&p2z[n]);
                                                r3 = _mm_mul_pd(r3, r0);
                                                _mm_store_pd(&fz[n], r3);
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
                                                r2 = _mm_load_pd(&p1x[n]);
                                                r3 = _mm_load_pd(&p2x[n]);
                                                r2 = _mm_mul_pd(r2, r1);
                                                r3 = _mm_mul_pd(r3, r0);
                                                r2 = _mm_add_pd(r2, r3);
                                                _mm_store_pd(&fx[n], r2);
                                                r2 = _mm_load_pd(&p1y[n]);
                                                r3 = _mm_load_pd(&p2y[n]);
                                                r2 = _mm_mul_pd(r2, r1);
                                                r3 = _mm_mul_pd(r3, r0);
                                                r2 = _mm_add_pd(r2, r3);
                                                _mm_store_pd(&fy[n], r2);
                                                r2 = _mm_load_pd(&p1z[n]);
                                                r3 = _mm_load_pd(&p2z[n]);
                                                r2 = _mm_mul_pd(r2, r1);
                                                r3 = _mm_mul_pd(r3, r0);
                                                r2 = _mm_add_pd(r2, r3);
                                                _mm_store_pd(&fz[n], r2);
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
void CINTx1i_2e(double *f, const double *g,
                const int li, const int lj, const int lk, const int ll,
                const double *ri, const CINTEnvVars *envs)
{
        const int nroots = envs->nrys_roots;
        const int di = envs->g_stride_i;
        DEF_GXYZ(const double, g, gx, gy, gz);
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
                const int dk = envs->g_stride_k;
                const int dl = envs->g_stride_l;
                const int dj = envs->g_stride_j;
                const double *p1x = gx + di;
                const double *p1y = gy + di;
                const double *p1z = gz + di;
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
                                        r3 = _mm_load_pd(&p1x[n]);
                                        r4 = _mm_load_pd(&gx [n]);
                                        r4 = _mm_mul_pd(r4, r0);
                                        r3 = _mm_add_pd(r3, r4);
                                        _mm_store_pd(&fx[n], r3);
                                        r3 = _mm_load_pd(&p1y[n]);
                                        r4 = _mm_load_pd(&gy [n]);
                                        r4 = _mm_mul_pd(r4, r1);
                                        r3 = _mm_add_pd(r3, r4);
                                        _mm_store_pd(&fy[n], r3);
                                        r3 = _mm_load_pd(&p1z[n]);
                                        r4 = _mm_load_pd(&gz [n]);
                                        r4 = _mm_mul_pd(r4, r2);
                                        r3 = _mm_add_pd(r3, r4);
                                        _mm_store_pd(&fz[n], r3);
                                }
                                ptr += di;
                        }
                } }
        }
}


/*
 * ( i x^1 j | kl )
 */
void CINTx1j_2e(double *f, const double *g,
                const int li, const int lj, const int lk, const int ll,
                const double *rj, const CINTEnvVars *envs)
{
        const int nroots = envs->nrys_roots;
        const int dj = envs->g_stride_j;
        DEF_GXYZ(const double, g, gx, gy, gz);
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
                const int di = envs->g_stride_i;
                const int dk = envs->g_stride_k;
                const int dl = envs->g_stride_l;
                const double *p1x = gx + dj;
                const double *p1y = gy + dj;
                const double *p1z = gz + dj;
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
                                        r3 = _mm_load_pd(&p1x[n]);
                                        r4 = _mm_load_pd(&gx [n]);
                                        r4 = _mm_mul_pd(r4, r0);
                                        r3 = _mm_add_pd(r3, r4);
                                        _mm_store_pd(&fx[n], r3);
                                        r3 = _mm_load_pd(&p1y[n]);
                                        r4 = _mm_load_pd(&gy [n]);
                                        r4 = _mm_mul_pd(r4, r1);
                                        r3 = _mm_add_pd(r3, r4);
                                        _mm_store_pd(&fy[n], r3);
                                        r3 = _mm_load_pd(&p1z[n]);
                                        r4 = _mm_load_pd(&gz [n]);
                                        r4 = _mm_mul_pd(r4, r2);
                                        r3 = _mm_add_pd(r3, r4);
                                        _mm_store_pd(&fz[n], r3);
                                }
                                ptr += di;
                        }
                } }
        }
}


/*
 * ( ij | x^1 k l )
 */
void CINTx1k_2e(double *f, const double *g,
                const int li, const int lj, const int lk, const int ll,
                const double *rk, const CINTEnvVars *envs)
{
        const int nroots = envs->nrys_roots;
        const int dk = envs->g_stride_k;
        DEF_GXYZ(const double, g, gx, gy, gz);
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
                const int di = envs->g_stride_i;
                const int dl = envs->g_stride_l;
                const int dj = envs->g_stride_j;
                const double *p1x = gx + dk;
                const double *p1y = gy + dk;
                const double *p1z = gz + dk;
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
                                        r3 = _mm_load_pd(&p1x[n]);
                                        r4 = _mm_load_pd(&gx [n]);
                                        r4 = _mm_mul_pd(r4, r0);
                                        r3 = _mm_add_pd(r3, r4);
                                        _mm_store_pd(&fx[n], r3);
                                        r3 = _mm_load_pd(&p1y[n]);
                                        r4 = _mm_load_pd(&gy [n]);
                                        r4 = _mm_mul_pd(r4, r1);
                                        r3 = _mm_add_pd(r3, r4);
                                        _mm_store_pd(&fy[n], r3);
                                        r3 = _mm_load_pd(&p1z[n]);
                                        r4 = _mm_load_pd(&gz [n]);
                                        r4 = _mm_mul_pd(r4, r2);
                                        r3 = _mm_add_pd(r3, r4);
                                        _mm_store_pd(&fz[n], r3);
                                }
                                ptr += di;
                        }
                } }
        }
}


/*
 * ( i j | x^1 kl )
 */
void CINTx1l_2e(double *f, const double *g,
                const int li, const int lj, const int lk, const int ll,
                const double *rl, const CINTEnvVars *envs)
{
        const int nroots = envs->nrys_roots;
        const int dl = envs->g_stride_l;
        DEF_GXYZ(const double, g, gx, gy, gz);
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
                const int di = envs->g_stride_i;
                const int dk = envs->g_stride_k;
                const int dj = envs->g_stride_j;
                const double *p1x = gx + dl;
                const double *p1y = gy + dl;
                const double *p1z = gz + dl;
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
                                        r3 = _mm_load_pd(&p1x[n]);
                                        r4 = _mm_load_pd(&gx [n]);
                                        r4 = _mm_mul_pd(r4, r0);
                                        r3 = _mm_add_pd(r3, r4);
                                        _mm_store_pd(&fx[n], r3);
                                        r3 = _mm_load_pd(&p1y[n]);
                                        r4 = _mm_load_pd(&gy [n]);
                                        r4 = _mm_mul_pd(r4, r1);
                                        r3 = _mm_add_pd(r3, r4);
                                        _mm_store_pd(&fy[n], r3);
                                        r3 = _mm_load_pd(&p1z[n]);
                                        r4 = _mm_load_pd(&gz [n]);
                                        r4 = _mm_mul_pd(r4, r2);
                                        r3 = _mm_add_pd(r3, r4);
                                        _mm_store_pd(&fz[n], r3);
                                }
                                ptr += di;
                        }
                } }
        }
}


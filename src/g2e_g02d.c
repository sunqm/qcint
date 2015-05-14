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


#include <pmmintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "cint_const.h"
#include "cint_bas.h"
#include "g2e.h"

#if defined(__GNUC__)
#define RESTRICT __restrict__
#else
#define RESTRICT
#endif

/*
 * g(nroots,0:nmax,0:mmax)
 */
void CINTg0_2e_2d(double *g, struct _BC *bc, const CINTEnvVars *envs)
{
        const FINT nroots = envs->nrys_roots;
        FINT i, j, m, n, off;
        double *gx = g;
        double *gy = g + envs->g_size;
        double *gz = g + envs->g_size * 2;

        for (i = 0; i < nroots; i++) {
                gx[i] = 1;
                gy[i] = 1;
                //gz[i] = w[i];
        }

        const FINT nmax = envs->li_ceil + envs->lj_ceil;
        const FINT mmax = envs->lk_ceil + envs->ll_ceil;
        const FINT dm = envs->g2d_klmax;
        const FINT dn = envs->g2d_ijmax;
        const double *RESTRICT c00x = bc->c00;
        const double *RESTRICT c00y = bc->c00+MXRYSROOTS;
        const double *RESTRICT c00z = bc->c00+MXRYSROOTS*2;
        const double *RESTRICT c0px = bc->c0p;
        const double *RESTRICT c0py = bc->c0p+MXRYSROOTS;
        const double *RESTRICT c0pz = bc->c0p+MXRYSROOTS*2;
        const double *RESTRICT b01 = bc->b01;
        const double *RESTRICT b00 = bc->b00;
        const double *RESTRICT b10 = bc->b10;
        double *p0x, *p0y, *p0z;
        const double *p1x, *p1y, *p1z, *p2x, *p2y, *p2z;
        __m128d r0, r1, r2, r3, r4;
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
                for (n = 1; n < nmax; n++) {
                        off = n * dn;
                        for (i = 0, j = off; i < nroots-1; i+=2, j+=2) {
                                r1 = _mm_load_pd(&b10[i]);
                                r2 = _mm_set1_pd(n);
                                r1 = _mm_mul_pd(r1, r2);
//FIXME: I assumed the alignment of stack.
// if not, change _mm_load_pd to _mm_loadu_pd
                                r4 = _mm_load_pd(&c00x[i]);
                                r2 = _mm_load_pd(&gx[j]);
                                r2 = _mm_mul_pd(r2, r4);
                                r3 = _mm_load_pd(&p1x[j]);
                                r3 = _mm_mul_pd(r3, r1);
                                r2 = _mm_add_pd(r2, r3);
                                _mm_store_pd(&p0x[j], r2);
                                r4 = _mm_load_pd(&c00y[i]);
                                r2 = _mm_load_pd(&gy[j]);
                                r2 = _mm_mul_pd(r2, r4);
                                r3 = _mm_load_pd(&p1y[j]);
                                r3 = _mm_mul_pd(r3, r1);
                                r2 = _mm_add_pd(r2, r3);
                                _mm_store_pd(&p0y[j], r2);
                                r4 = _mm_load_pd(&c00z[i]);
                                r2 = _mm_load_pd(&gz[j]);
                                r2 = _mm_mul_pd(r2, r4);
                                r3 = _mm_load_pd(&p1z[j]);
                                r3 = _mm_mul_pd(r3, r1);
                                r2 = _mm_add_pd(r2, r3);
                                _mm_store_pd(&p0z[j], r2);
                        }
                        if (i < nroots) {
                                p0x[j] = c00x[i] * gx[j] + n * b10[i] * p1x[j];
                                p0y[j] = c00y[i] * gy[j] + n * b10[i] * p1y[j];
                                p0z[j] = c00z[i] * gz[j] + n * b10[i] * p1z[j];
                        }
                }
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
                for (m = 1; m < mmax; m++) {
                        off = m * dm;
                        for (i = 0, j = off; i < nroots-1; i+=2, j+=2) {
                                r2 = _mm_set1_pd(m);
                                r1 = _mm_load_pd(&b01[i]);
                                r1 = _mm_mul_pd(r1, r2);
                                r4 = _mm_load_pd(&c0px[i]);
                                r2 = _mm_load_pd(&gx[j]);
                                r2 = _mm_mul_pd(r2, r4);
                                r3 = _mm_load_pd(&p1x[j]);
                                r3 = _mm_mul_pd(r3, r1);
                                r2 = _mm_add_pd(r2, r3);
                                _mm_store_pd(&p0x[j], r2);
                                r4 = _mm_load_pd(&c0py[i]);
                                r2 = _mm_load_pd(&gy[j]);
                                r2 = _mm_mul_pd(r2, r4);
                                r3 = _mm_load_pd(&p1y[j]);
                                r3 = _mm_mul_pd(r3, r1);
                                r2 = _mm_add_pd(r2, r3);
                                _mm_store_pd(&p0y[j], r2);
                                r4 = _mm_load_pd(&c0pz[i]);
                                r2 = _mm_load_pd(&gz[j]);
                                r2 = _mm_mul_pd(r2, r4);
                                r3 = _mm_load_pd(&p1z[j]);
                                r3 = _mm_mul_pd(r3, r1);
                                r2 = _mm_add_pd(r2, r3);
                                _mm_store_pd(&p0z[j], r2);
                        }
                        if (i < nroots) {
                                p0x[j] = c0px[i] * gx[j] + m * b01[i] * p1x[j];
                                p0y[j] = c0py[i] * gy[j] + m * b01[i] * p1y[j];
                                p0z[j] = c0pz[i] * gz[j] + m * b01[i] * p1z[j];
                        }
                }
        }

        if (nmax > 0 && mmax > 0) {
                p0x = gx + dn;
                p0y = gy + dn;
                p0z = gz + dn;
                p1x = gx - dn;
                p1y = gy - dn;
                p1z = gz - dn;
                p2x = gx - dm;
                p2y = gy - dm;
                p2z = gz - dm;
                // gx(irys,1,1) = c0p(irys)*gx(irys,0,1)
                // + b00(irys)*gx(irys,0,0)
                for (i = 0; i < nroots; i++) {
                        p0x[i+dm] = c0px[i] * p0x[i] + b00[i] * gx[i];
                        p0y[i+dm] = c0py[i] * p0y[i] + b00[i] * gy[i];
                        p0z[i+dm] = c0pz[i] * p0z[i] + b00[i] * gz[i];
                }

                // gx(irys,m+1,1) = c0p(irys)*gx(irys,m,1)
                // + m*b01(irys)*gx(irys,m-1,1)
                // + b00(irys)*gx(irys,m,0)
                for (m = 1; m < mmax; m++) {
                        off = m * dm + dn;
                        for (i = 0, j = off; i < nroots-1; i+=2, j+=2) {
                                r2 = _mm_set1_pd(m);
                                r1 = _mm_load_pd(&b01[i]);
                                r1 = _mm_mul_pd(r1, r2);
                                r0 = _mm_load_pd(&b00[i]);
                                r4 = _mm_load_pd(&c0px[i]);
                                r2 = _mm_load_pd(&gx[j]);
                                r2 = _mm_mul_pd(r2, r4);
                                r3 = _mm_load_pd(&p1x[j]);
                                r3 = _mm_mul_pd(r3, r0);
                                r2 = _mm_add_pd(r2, r3);
                                r3 = _mm_load_pd(&p2x[j]);
                                r3 = _mm_mul_pd(r3, r1);
                                r2 = _mm_add_pd(r2, r3);
                                _mm_store_pd(&gx[j+dm], r2);
                                r4 = _mm_load_pd(&c0py[i]);
                                r2 = _mm_load_pd(&gy[j]);
                                r2 = _mm_mul_pd(r2, r4);
                                r3 = _mm_load_pd(&p1y[j]);
                                r3 = _mm_mul_pd(r3, r0);
                                r2 = _mm_add_pd(r2, r3);
                                r3 = _mm_load_pd(&p2y[j]);
                                r3 = _mm_mul_pd(r3, r1);
                                r2 = _mm_add_pd(r2, r3);
                                _mm_store_pd(&gy[j+dm], r2);
                                r4 = _mm_load_pd(&c0pz[i]);
                                r2 = _mm_load_pd(&gz[j]);
                                r2 = _mm_mul_pd(r2, r4);
                                r3 = _mm_load_pd(&p1z[j]);
                                r3 = _mm_mul_pd(r3, r0);
                                r2 = _mm_add_pd(r2, r3);
                                r3 = _mm_load_pd(&p2z[j]);
                                r3 = _mm_mul_pd(r3, r1);
                                r2 = _mm_add_pd(r2, r3);
                                _mm_store_pd(&gz[j+dm], r2);
                        }
                        if (i < nroots) {
                                gx[j+dm] = c0px[i]*gx[j] + m*b01[i]*p2x[j] +b00[i]*p1x[j];
                                gy[j+dm] = c0py[i]*gy[j] + m*b01[i]*p2y[j] +b00[i]*p1y[j];
                                gz[j+dm] = c0pz[i]*gz[j] + m*b01[i]*p2z[j] +b00[i]*p1z[j];
                        }
                }

                // gx(irys,m,n+1) = c00(irys)*gx(irys,m,n)
                // + n*b10(irys)*gx(irys,m,n-1)
                // + m*b00(irys)*gx(irys,m-1,n)
                for (m = 1; m <= mmax; m++) {
                        for (n = 1; n < nmax; n++) {
                                off = m * dm + n * dn;
                                for (i = 0, j = off; i < nroots-1; i+=2, j+=2) {
                                        r2 = _mm_set1_pd(n);
                                        r1 = _mm_load_pd(&b10[i]);
                                        r1 = _mm_mul_pd(r1, r2);
                                        r2 = _mm_set1_pd(m);
                                        r0 = _mm_load_pd(&b00[i]);
                                        r0 = _mm_mul_pd(r0, r2);
                                        r4 = _mm_load_pd(&c00x[i]);
                                        r2 = _mm_load_pd(&gx[j]);
                                        r2 = _mm_mul_pd(r2, r4);
                                        r3 = _mm_load_pd(&p1x[j]);
                                        r3 = _mm_mul_pd(r3, r1);
                                        r2 = _mm_add_pd(r2, r3);
                                        r3 = _mm_load_pd(&p2x[j]);
                                        r3 = _mm_mul_pd(r3, r0);
                                        r2 = _mm_add_pd(r2, r3);
                                        _mm_store_pd(&p0x[j], r2);
                                        r4 = _mm_load_pd(&c00y[i]);
                                        r2 = _mm_load_pd(&gy[j]);
                                        r2 = _mm_mul_pd(r2, r4);
                                        r3 = _mm_load_pd(&p1y[j]);
                                        r3 = _mm_mul_pd(r3, r1);
                                        r2 = _mm_add_pd(r2, r3);
                                        r3 = _mm_load_pd(&p2y[j]);
                                        r3 = _mm_mul_pd(r3, r0);
                                        r2 = _mm_add_pd(r2, r3);
                                        _mm_store_pd(&p0y[j], r2);
                                        r4 = _mm_load_pd(&c00z[i]);
                                        r2 = _mm_load_pd(&gz[j]);
                                        r2 = _mm_mul_pd(r2, r4);
                                        r3 = _mm_load_pd(&p1z[j]);
                                        r3 = _mm_mul_pd(r3, r1);
                                        r2 = _mm_add_pd(r2, r3);
                                        r3 = _mm_load_pd(&p2z[j]);
                                        r3 = _mm_mul_pd(r3, r0);
                                        r2 = _mm_add_pd(r2, r3);
                                        _mm_store_pd(&p0z[j], r2);
                                }
                                if (i < nroots) {
                                        p0x[j] = c00x[i]*gx[j] +n*b10[i]*p1x[j] + m*b00[i]*p2x[j];
                                        p0y[j] = c00y[i]*gy[j] +n*b10[i]*p1y[j] + m*b00[i]*p2y[j];
                                        p0z[j] = c00z[i]*gz[j] +n*b10[i]*p1z[j] + m*b00[i]*p2z[j];
                                }
                        }
                }
        }
}


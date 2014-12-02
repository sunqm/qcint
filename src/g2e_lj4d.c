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
#include <stdlib.h>
#include "cint_const.h"
#include "g2e.h"

struct _BC {
        double c00[MXRYSROOTS*3];
        double c0p[MXRYSROOTS*3];
        double b01[MXRYSROOTS];
        double b00[MXRYSROOTS];
        double b10[MXRYSROOTS];
};

void CINTg0_2e_2d(double *g, struct _BC *bc, const CINTEnvVars *envs);
void CINTg0_lj2d_4d(double *g, const CINTEnvVars *envs);

/************* some special g0_4d results *************/
/* 4 digits stand for i_ceil, k_ceil, l_ceil, j_ceil */
static inline void _g0_lj_4d_0001(double *g, double *c,
                                  const double *r)
{
        g[0] = 1;
        g[1] = c[0];
        g[2] = 1;
        g[3] = c[MXRYSROOTS];
        //g[4] = w[0];
        g[5] = c[MXRYSROOTS*2] * g[4];
        //g[0] = 1;
        //g[2] = c[0];
        //g[4] = 1;
        //g[6] = c[MXRYSROOTS];
        ////g[8] = w[0];
        //g[10] = c[MXRYSROOTS*2] * g[8];
}
static inline void _g0_lj_4d_1000(double *g, double *c,
                                  const double *r)
{
        g[0] = 1;
        g[1] = r[0] + c[0];
        g[4] = 1;
        g[5] = r[1] + c[MXRYSROOTS];
        //g[8] = w[0];
        g[9] =(r[2] + c[MXRYSROOTS*2]) * g[8];
        //g[0] = 1;
        //g[2] = r[0] + c[0];
        //g[8] = 1;
        //g[10] = r[1] + c[MXRYSROOTS];
        ////g[16] = w[0];
        //g[18] =(r[2] + c[MXRYSROOTS*2]) * g[16];
}
static inline void _g0_lj_4d_0002(double *g, double *c, double *b,
                                  const double *r)
{
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = c[0];
        g[3 ] = c[1];
        g[4 ] = c[0] * c[0] + b[0];
        g[5 ] = c[1] * c[1] + b[1];
        g[6 ] = 1;
        g[7 ] = 1;
        g[8 ] = c[MXRYSROOTS];
        g[9 ] = c[MXRYSROOTS+1];
        g[10] = c[MXRYSROOTS] * c[MXRYSROOTS] + b[0];
        g[11] = c[MXRYSROOTS+1] * c[MXRYSROOTS+1] + b[1];
        //g[12] = w[0];
        //g[13] = w[1];
        g[14] = c[MXRYSROOTS*2] * g[12];
        g[15] = c[MXRYSROOTS*2+1] * g[13];
        g[16] =(c[MXRYSROOTS*2] * c[MXRYSROOTS*2] + b[0])* g[12];
        g[17] =(c[MXRYSROOTS*2+1] * c[MXRYSROOTS*2+1] + b[1])* g[13];
}
static inline void _g0_lj_4d_1001(double *g, double *c, double *b,
                                  const double *r)
{
        double rc[] = {r[0]+c[0], r[0]+c[1],
                       r[1]+c[MXRYSROOTS], r[1]+c[MXRYSROOTS+1],
                       r[2]+c[MXRYSROOTS*2], r[2]+c[MXRYSROOTS*2+1]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[4 ] = c[0];
        g[5 ] = c[1];
        g[6 ] = rc[0] * c[0] + b[0];
        g[7 ] = rc[1] * c[1] + b[1];
        g[12] = 1;
        g[13] = 1;
        g[14] = rc[2];
        g[15] = rc[3];
        g[16] = c[MXRYSROOTS];
        g[17] = c[MXRYSROOTS+1];
        g[18] = rc[2] * c[MXRYSROOTS] + b[0];
        g[19] = rc[3] * c[MXRYSROOTS+1] + b[1];
        //g[24] = w[0];
        //g[25] = w[1];
        g[26] = rc[4] * g[24];
        g[27] = rc[5] * g[25];
        g[28] = c[MXRYSROOTS*2] * g[24];
        g[29] = c[MXRYSROOTS*2+1] * g[25];
        g[30] =(rc[4] * c[MXRYSROOTS*2] + b[0])* g[24];
        g[31] =(rc[5] * c[MXRYSROOTS*2+1] + b[1])* g[25];
}
static inline void _g0_lj_4d_2000(double *g, double *c, double *b,
                                  const double *r)
{
        double rc[] = {r[0]+c[0], r[0]+c[1],
                       r[1]+c[MXRYSROOTS], r[1]+c[MXRYSROOTS+1],
                       r[2]+c[MXRYSROOTS*2], r[2]+c[MXRYSROOTS*2+1]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[4 ] = rc[0] * rc[0] + b[0];
        g[5 ] = rc[1] * rc[1] + b[1];
        g[18] = 1;
        g[19] = 1;
        g[20] = rc[2];
        g[21] = rc[3];
        g[22] = rc[2] * rc[2] + b[0];
        g[23] = rc[3] * rc[3] + b[1];
        //g[36] = w[0];
        //g[37] = w[1];
        g[38] = rc[4] * g[36];
        g[39] = rc[5] * g[37];
        g[40] =(rc[4] * rc[4] + b[0])* g[36];
        g[41] =(rc[5] * rc[5] + b[1])* g[37];
}
static inline void _g0_lj_4d_0003(double *g, double *c, double *b,
                                  const double *r)
{
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = c[0];
        g[3 ] = c[1];
        g[4 ] = c[0] * c[0] + b[0];
        g[5 ] = c[1] * c[1] + b[1];
        g[6 ] = c[0] *(c[0] * c[0] + 3 * b[0]);
        g[7 ] = c[1] *(c[1] * c[1] + 3 * b[1]);
        g[8 ] = 1;
        g[9 ] = 1;
        g[10] = c[MXRYSROOTS];
        g[11] = c[MXRYSROOTS+1];
        g[12] = c[MXRYSROOTS] * c[MXRYSROOTS] + b[0];
        g[13] = c[MXRYSROOTS+1] * c[MXRYSROOTS+1] + b[1];
        g[14] = c[MXRYSROOTS] *(c[MXRYSROOTS] * c[MXRYSROOTS] + 3 * b[0]);
        g[15] = c[MXRYSROOTS+1] *(c[MXRYSROOTS+1] * c[MXRYSROOTS+1] + 3 * b[1]);
        //g[16] = w[0];
        //g[17] = w[1];
        g[18] = c[MXRYSROOTS*2] * g[16];
        g[19] = c[MXRYSROOTS*2+1] * g[17];
        g[20] =(c[MXRYSROOTS*2] * c[MXRYSROOTS*2] + b[0])* g[16];
        g[21] =(c[MXRYSROOTS*2+1] * c[MXRYSROOTS*2+1] + b[1])* g[17];
        g[22] =(c[MXRYSROOTS*2] * c[MXRYSROOTS*2] + 3 * b[0])* c[MXRYSROOTS*2] * g[16];
        g[23] =(c[MXRYSROOTS*2+1] * c[MXRYSROOTS*2+1] + 3 * b[1])* c[MXRYSROOTS*2+1] * g[17];
}
static inline void _g0_lj_4d_1002(double *g, double *c, double *b,
                                  const double *r)
{
/*
        __m128d rcx = _mm_set_pd(r[0]+c[1], r[0]+c[0]);
        __m128d cx = _mm_set_pd(c[1], c[0]);
        __m128d rb = _mm_loadu_pd(&b[0]);
        __m128d m0, m1, m2;
        g[0 ] = 1;
        g[1 ] = 1;
        //g[2 ] = rc[0];
        //g[3 ] = rc[1];
        _mm_store_pd(&g[2], rcx);
        g[4 ] = c[0];
        g[5 ] = c[1];
        //g[6 ] = rc[0] * c[0] + b[0];
        //g[7 ] = rc[1] * c[1] + b[1];
        m1 = _mm_mul_pd(rcx, cx);
        m1 = _mm_add_pd(m1, rb);
        _mm_store_pd(&g[6], m1);
        //g[8 ] = c[0] * c[0] + b[0];
        //g[9 ] = c[1] * c[1] + b[1];
        m0 = _mm_mul_pd(cx, cx);
        m1 = _mm_add_pd(m0, rb);
        _mm_store_pd(&g[8], m1);
        //g[10] = rc[0]*c[0]*c[0] + b[0]*(rc[0]+2*c[0]);
        //g[11] = rc[1]*c[1]*c[1] + b[1]*(rc[1]+2*c[1]);
        m0 = _mm_mul_pd(rcx, m0);
        m1 = _mm_add_pd(cx, cx);
        m1 = _mm_add_pd(m1, rcx);
        m1 = _mm_mul_pd(m1, rb);
        m1 = _mm_add_pd(m1, m0);
        _mm_store_pd(&g[10], m1);
        __m128d rcy = _mm_set_pd(r[1]+c[MXRYSROOTS+1], r[1]+c[MXRYSROOTS]);
        __m128d cy = _mm_set_pd(c[MXRYSROOTS+1], c[MXRYSROOTS]);
        g[16] = 1;
        g[17] = 1;
        //g[18] = rc[2];
        //g[19] = rc[3];
        _mm_store_pd(&g[18], rcy);
        g[20] = c[MXRYSROOTS];
        g[21] = c[MXRYSROOTS+1];
        //g[22] = rc[2] * c[MXRYSROOTS] + b[0];
        //g[23] = rc[3] * c[MXRYSROOTS+1] + b[1];
        m1 = _mm_mul_pd(rcy, cy);
        m1 = _mm_add_pd(m1, rb);
        _mm_store_pd(&g[22], m1);
        //g[24] = c[MXRYSROOTS] * c[MXRYSROOTS] + b[0];
        //g[25] = c[MXRYSROOTS+1] * c[MXRYSROOTS+1] + b[1];
        m0 = _mm_mul_pd(cy, cy);
        m1 = _mm_add_pd(m0, rb);
        _mm_store_pd(&g[24], m1);
        //g[26] = rc[2]*c[MXRYSROOTS]*c[MXRYSROOTS] + b[0]*(rc[2]+2*c[MXRYSROOTS]);
        //g[27] = rc[3]*c[MXRYSROOTS+1]*c[MXRYSROOTS+1] + b[1]*(rc[3]+2*c[MXRYSROOTS+1]);
        m0 = _mm_mul_pd(rcy, m0);
        m1 = _mm_add_pd(cy, cy);
        m1 = _mm_add_pd(m1, rcy);
        m1 = _mm_mul_pd(m1, rb);
        m1 = _mm_add_pd(m1, m0);
        _mm_store_pd(&g[26], m1);
        __m128d rcz = _mm_set_pd(r[2]+c[MXRYSROOTS*2+1], r[2]+c[MXRYSROOTS*2]);
        __m128d cz = _mm_set_pd(c[MXRYSROOTS*2+1], c[MXRYSROOTS*2]);
        //g[32] = w[0];
        //g[33] = w[1];
        //g[34] = rc[4] * g[32];
        //g[35] = rc[5] * g[33];
        m2 = _mm_load_pd(&g[32]);
        m1 = _mm_mul_pd(m2, rcz);
        _mm_store_pd(&g[34], m1);
        //g[36] = c[MXRYSROOTS*2] * g[32];
        //g[37] = c[MXRYSROOTS*2+1] * g[33];
        m1 = _mm_mul_pd(m2, cz);
        _mm_store_pd(&g[36], m1);
        //g[38] =(rc[4] * c[MXRYSROOTS*2] + b[0])* g[32];
        //g[39] =(rc[5] * c[MXRYSROOTS*2+1] + b[1])* g[33];
        m1 = _mm_mul_pd(rcz, cz);
        m1 = _mm_add_pd(m1, rb);
        m1 = _mm_mul_pd(m1, m2);
        _mm_store_pd(&g[38], m1);
        //g[40] =(c[MXRYSROOTS*2] * c[MXRYSROOTS*2] + b[0])* g[32];
        //g[41] =(c[MXRYSROOTS*2+1] * c[MXRYSROOTS*2+1] + b[1])* g[33];
        m0 = _mm_mul_pd(cz, cz);
        m1 = _mm_add_pd(m0, rb);
        m1 = _mm_mul_pd(m1, m2);
        _mm_store_pd(&g[40], m1);
        //g[42] =(rc[4]*c[MXRYSROOTS*2]*c[MXRYSROOTS*2]+b[0]*(rc[4]+2*c[MXRYSROOTS*2]))*g[32];
        //g[43] =(rc[5]*c[MXRYSROOTS*2+1]*c[MXRYSROOTS*2+1]+b[1]*(rc[5]+2*c[MXRYSROOTS*2+1]))*g[33];
        m0 = _mm_mul_pd(rcz, m0);
        m1 = _mm_add_pd(cz, cz);
        m1 = _mm_add_pd(m1, rcz);
        m1 = _mm_mul_pd(m1, rb);
        m1 = _mm_add_pd(m1, m0);
        m1 = _mm_mul_pd(m1, m2);
        _mm_store_pd(&g[42], m1);
*/
        double rc[] = {r[0]+c[0], r[0]+c[1],
                       r[1]+c[MXRYSROOTS], r[1]+c[MXRYSROOTS+1],
                       r[2]+c[MXRYSROOTS*2], r[2]+c[MXRYSROOTS*2+1]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[4 ] = c[0];
        g[5 ] = c[1];
        g[6 ] = rc[0] * c[0] + b[0];
        g[7 ] = rc[1] * c[1] + b[1];
        g[8 ] = c[0] * c[0] + b[0];
        g[9 ] = c[1] * c[1] + b[1];
        g[10] = rc[0]*c[0]*c[0] + b[0]*(rc[0]+2*c[0]);
        g[11] = rc[1]*c[1]*c[1] + b[1]*(rc[1]+2*c[1]);
        g[16] = 1;
        g[17] = 1;
        g[18] = rc[2];
        g[19] = rc[3];
        g[20] = c[MXRYSROOTS];
        g[21] = c[MXRYSROOTS+1];
        g[22] = rc[2] * c[MXRYSROOTS] + b[0];
        g[23] = rc[3] * c[MXRYSROOTS+1] + b[1];
        g[24] = c[MXRYSROOTS] * c[MXRYSROOTS] + b[0];
        g[25] = c[MXRYSROOTS+1] * c[MXRYSROOTS+1] + b[1];
        g[26] = rc[2]*c[MXRYSROOTS]*c[MXRYSROOTS] + b[0]*(rc[2]+2*c[MXRYSROOTS]);
        g[27] = rc[3]*c[MXRYSROOTS+1]*c[MXRYSROOTS+1] + b[1]*(rc[3]+2*c[MXRYSROOTS+1]);
        //g[32] = w[0];
        //g[33] = w[1];
        g[34] = rc[4] * g[32];
        g[35] = rc[5] * g[33];
        g[36] = c[MXRYSROOTS*2] * g[32];
        g[37] = c[MXRYSROOTS*2+1] * g[33];
        g[38] =(rc[4] * c[MXRYSROOTS*2] + b[0])* g[32];
        g[39] =(rc[5] * c[MXRYSROOTS*2+1] + b[1])* g[33];
        g[40] =(c[MXRYSROOTS*2] * c[MXRYSROOTS*2] + b[0])* g[32];
        g[41] =(c[MXRYSROOTS*2+1] * c[MXRYSROOTS*2+1] + b[1])* g[33];
        g[42] =(rc[4]*c[MXRYSROOTS*2]*c[MXRYSROOTS*2]+b[0]*(rc[4]+2*c[MXRYSROOTS*2]))*g[32];
        g[43] =(rc[5]*c[MXRYSROOTS*2+1]*c[MXRYSROOTS*2+1]+b[1]*(rc[5]+2*c[MXRYSROOTS*2+1]))*g[33];
}
static inline void _g0_lj_4d_2001(double *g, double *c, double *b,
                                  const double *r)
{
        double rc[] = {r[0]+c[0], r[0]+c[1],
                       r[1]+c[MXRYSROOTS], r[1]+c[MXRYSROOTS+1],
                       r[2]+c[MXRYSROOTS*2], r[2]+c[MXRYSROOTS*2+1]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[4 ] = rc[0] * rc[0] + b[0];
        g[5 ] = rc[1] * rc[1] + b[1];
        g[6 ] = c[0];
        g[7 ] = c[1];
        g[8 ] = c[0] * rc[0] + b[0];
        g[9 ] = c[1] * rc[1] + b[1];
        g[10] = c[0]*rc[0]*rc[0] + b[0]*(2*rc[0]+c[0]);
        g[11] = c[1]*rc[1]*rc[1] + b[1]*(2*rc[1]+c[1]);
        g[24] = 1;
        g[25] = 1;
        g[26] = rc[2];
        g[27] = rc[3];
        g[28] = rc[2] * rc[2] + b[0];
        g[29] = rc[3] * rc[3] + b[1];
        g[30] = c[MXRYSROOTS];
        g[31] = c[MXRYSROOTS+1];
        g[32] = c[MXRYSROOTS] * rc[2] + b[0];
        g[33] = c[MXRYSROOTS+1] * rc[3] + b[1];
        g[34] = c[MXRYSROOTS]*rc[2]*rc[2] + b[0]*(2*rc[2]+c[MXRYSROOTS]);
        g[35] = c[MXRYSROOTS+1]*rc[3]*rc[3] + b[1]*(2*rc[3]+c[MXRYSROOTS+1]);
        //g[48] = w[0];
        //g[49] = w[1];
        g[50] = rc[4] * g[48];
        g[51] = rc[5] * g[49];
        g[52] =(rc[4] * rc[4] + b[0])* g[48];
        g[53] =(rc[5] * rc[5] + b[1])* g[49];
        g[54] = c[MXRYSROOTS*2] * g[48];
        g[55] = c[MXRYSROOTS*2+1] * g[49];
        g[56] =(c[MXRYSROOTS*2] * rc[4] + b[0])* g[48];
        g[57] =(c[MXRYSROOTS*2+1] * rc[5] + b[1])* g[49];
        g[58] =(c[MXRYSROOTS*2]*rc[4]*rc[4] + b[0]*(2*rc[4]+c[MXRYSROOTS*2]))* g[48];
        g[59] =(c[MXRYSROOTS*2+1]*rc[5]*rc[5] + b[1]*(2*rc[5]+c[MXRYSROOTS*2+1]))* g[49];
}
static inline void _g0_lj_4d_3000(double *g, double *c, double *b,
                                  const double *r)
{
        double rc[] = {r[0]+c[0], r[0]+c[1],
                       r[1]+c[MXRYSROOTS], r[1]+c[MXRYSROOTS+1],
                       r[2]+c[MXRYSROOTS*2], r[2]+c[MXRYSROOTS*2+1]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[4 ] = rc[0] * rc[0] + b[0];
        g[5 ] = rc[1] * rc[1] + b[1];
        g[6 ] = rc[0] *(rc[0] * rc[0] + 3 * b[0]);
        g[7 ] = rc[1] *(rc[1] * rc[1] + 3 * b[1]);
        g[32] = 1;
        g[33] = 1;
        g[34] = rc[2];
        g[35] = rc[3];
        g[36] = rc[2] * rc[2] + b[0];
        g[37] = rc[3] * rc[3] + b[1];
        g[38] = rc[2] *(rc[2] * rc[2] + 3 * b[0]);
        g[39] = rc[3] *(rc[3] * rc[3] + 3 * b[1]);
        //g[64] = w[0];
        //g[65] = w[1];
        g[66] = rc[4] * g[64];
        g[67] = rc[5] * g[65];
        g[68] =(rc[4] * rc[4] + b[0])* g[64];
        g[69] =(rc[5] * rc[5] + b[1])* g[65];
        g[70] =(rc[4] * rc[4] + 3 * b[0])* rc[4] * g[64];
        g[71] =(rc[5] * rc[5] + 3 * b[1])* rc[5] * g[65];
}
static inline void _g0_lj_4d_0011(double *g, double *c0, double *cp, double *b,
                                  const double *r0, const double *rp)
{
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = cp[0];
        g[3 ] = cp[1];
        g[4 ] = c0[0];
        g[5 ] = c0[1];
        g[6 ] = cp[0] * c0[0] + b[0];
        g[7 ] = cp[1] * c0[1] + b[1];
        g[8 ] = 1;
        g[9 ] = 1;
        g[10] = cp[MXRYSROOTS];
        g[11] = cp[MXRYSROOTS+1];
        g[12] = c0[MXRYSROOTS];
        g[13] = c0[MXRYSROOTS+1];
        g[14] = cp[MXRYSROOTS] * c0[MXRYSROOTS] + b[0];
        g[15] = cp[MXRYSROOTS+1] * c0[MXRYSROOTS+1] + b[1];
        //g[16] = w[0];
        //g[17] = w[1];
        g[18] = cp[MXRYSROOTS*2] * g[16];
        g[19] = cp[MXRYSROOTS*2+1] * g[17];
        g[20] = c0[MXRYSROOTS*2] * g[16];
        g[21] = c0[MXRYSROOTS*2+1] * g[17];
        g[22] =(cp[MXRYSROOTS*2] * c0[MXRYSROOTS*2] + b[0]) * g[16];
        g[23] =(cp[MXRYSROOTS*2+1] * c0[MXRYSROOTS*2+1] + b[1]) * g[17];
}
static inline void _g0_lj_4d_1010(double *g, double *c0, double *cp, double *b,
                                  const double *r0, const double *rp)
{
        double rc[] = {r0[0]+c0[0], r0[0]+c0[1],
                       r0[1]+c0[MXRYSROOTS], r0[1]+c0[MXRYSROOTS+1],
                       r0[2]+c0[MXRYSROOTS*2], r0[2]+c0[MXRYSROOTS*2+1]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[4 ] = cp[0];
        g[5 ] = cp[1];
        g[6 ] = rc[0] * cp[0] + b[0];
        g[7 ] = rc[1] * cp[1] + b[1];
        g[16] = 1;
        g[17] = 1;
        g[18] = rc[2];
        g[19] = rc[3];
        g[20] = cp[MXRYSROOTS];
        g[21] = cp[MXRYSROOTS+1];
        g[22] = rc[2] * cp[MXRYSROOTS] + b[0];
        g[23] = rc[3] * cp[MXRYSROOTS+1] + b[1];
        //g[32] = w[0];
        //g[33] = w[1];
        g[34] = rc[4] * g[32];
        g[35] = rc[5] * g[33];
        g[36] = cp[MXRYSROOTS*2] * g[32];
        g[37] = cp[MXRYSROOTS*2+1] * g[33];
        g[38] =(rc[4]*cp[MXRYSROOTS*2] + b[0]) * g[32];
        g[39] =(rc[5]*cp[MXRYSROOTS*2+1] + b[1]) * g[33];
}
static inline void _g0_lj_4d_0101(double *g, double *c0, double *cp, double *b,
                                  const double *r0, const double *rp)
{
        double rc[] = {rp[0]+cp[0], rp[0]+cp[1],
                       rp[1]+cp[MXRYSROOTS], rp[1]+cp[MXRYSROOTS+1],
                       rp[2]+cp[MXRYSROOTS*2], rp[2]+cp[MXRYSROOTS*2+1]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[8 ] = c0[0];
        g[9 ] = c0[1];
        g[10] = rc[0] * c0[0] + b[0];
        g[11] = rc[1] * c0[1] + b[1];
        g[16] = 1;
        g[17] = 1;
        g[18] = rc[2];
        g[19] = rc[3];
        g[24] = c0[MXRYSROOTS];
        g[25] = c0[MXRYSROOTS+1];
        g[26] = rc[2] * c0[MXRYSROOTS] + b[0];
        g[27] = rc[3] * c0[MXRYSROOTS+1] + b[1];
        //g[32] = w[0];
        //g[33] = w[1];
        g[34] = rc[4] * g[32];
        g[35] = rc[5] * g[33];
        g[40] = c0[MXRYSROOTS*2] * g[32];
        g[41] = c0[MXRYSROOTS*2+1] * g[33];
        g[42] =(rc[4]*c0[MXRYSROOTS*2] + b[0]) * g[32];
        g[43] =(rc[5]*c0[MXRYSROOTS*2+1] + b[1]) * g[33];
}
static inline void _g0_lj_4d_1100(double *g, double *c0, double *cp, double *b,
                                  const double *r0, const double *rp)
{
        double rc0[] = {r0[0]+c0[0], r0[0]+c0[1],
                        r0[1]+c0[MXRYSROOTS], r0[1]+c0[MXRYSROOTS+1],
                        r0[2]+c0[MXRYSROOTS*2], r0[2]+c0[MXRYSROOTS*2+1]};
        double rcp[] = {rp[0]+cp[0], rp[0]+cp[1],
                        rp[1]+cp[MXRYSROOTS], rp[1]+cp[MXRYSROOTS+1],
                        rp[2]+cp[MXRYSROOTS*2], rp[2]+cp[MXRYSROOTS*2+1]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc0[0];
        g[3 ] = rc0[1];
        g[4 ] = rcp[0];
        g[5 ] = rcp[1];
        g[6 ] = rc0[0] * rcp[0] + b[0];
        g[7 ] = rc0[1] * rcp[1] + b[1];
        g[32] = 1;
        g[33] = 1;
        g[34] = rc0[2];
        g[35] = rc0[3];
        g[36] = rcp[2];
        g[37] = rcp[3];
        g[38] = rc0[2] * rcp[2] + b[0];
        g[39] = rc0[3] * rcp[3] + b[1];
        //g[64] = w[0];
        //g[65] = w[1];
        g[66] = rc0[4] * g[64];
        g[67] = rc0[5] * g[65];
        g[68] = rcp[4] * g[64];
        g[69] = rcp[5] * g[65];
        g[70] =(rc0[4]*rcp[4] + b[0]) * g[64];
        g[71] =(rc0[5]*rcp[5] + b[1]) * g[65];
}
static inline void _g0_lj_4d_0021(double *g, double *c0, double *cp,
                                  double *b0, double *b1)
{
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = cp[0];
        g[3 ] = cp[1];
        g[4 ] = cp[0] * cp[0] + b1[0];
        g[5 ] = cp[1] * cp[1] + b1[1];
        g[6 ] = c0[0];
        g[7 ] = c0[1];
        g[8 ] = cp[0] * c0[0] + b0[0];
        g[9 ] = cp[1] * c0[1] + b0[1];
        g[10] = c0[0] * g[4] + 2 * b0[0] * cp[0];
        g[11] = c0[1] * g[5] + 2 * b0[1] * cp[1];
        g[12] = 1;
        g[13] = 1;
        g[14] = cp[MXRYSROOTS];
        g[15] = cp[MXRYSROOTS+1];
        g[16] = cp[MXRYSROOTS] * cp[MXRYSROOTS] + b1[0];
        g[17] = cp[MXRYSROOTS+1] * cp[MXRYSROOTS+1] + b1[1];
        g[18] = c0[MXRYSROOTS];
        g[19] = c0[MXRYSROOTS+1];
        g[20] = cp[MXRYSROOTS] * c0[MXRYSROOTS] + b0[0];
        g[21] = cp[MXRYSROOTS+1] * c0[MXRYSROOTS+1] + b0[1];
        g[22] = c0[MXRYSROOTS] * g[16] + 2 * b0[0] * cp[MXRYSROOTS];
        g[23] = c0[MXRYSROOTS+1] * g[17] + 2 * b0[1] * cp[MXRYSROOTS+1];
        //g[24] = w[0];
        //g[25] = w[1];
        g[26] = cp[MXRYSROOTS*2] * g[24];
        g[27] = cp[MXRYSROOTS*2+1] * g[25];
        g[28] =(cp[MXRYSROOTS*2] * cp[MXRYSROOTS*2] + b1[0]) * g[24];
        g[29] =(cp[MXRYSROOTS*2+1] * cp[MXRYSROOTS*2+1] + b1[1]) * g[25];
        g[30] = c0[MXRYSROOTS*2] * g[24];
        g[31] = c0[MXRYSROOTS*2+1] * g[25];
        g[32] =(cp[MXRYSROOTS*2] * c0[MXRYSROOTS*2] + b0[0]) * g[24];
        g[33] =(cp[MXRYSROOTS*2+1] * c0[MXRYSROOTS*2+1] + b0[1]) * g[25];
        g[34] = c0[MXRYSROOTS*2] * g[28] + 2 * b0[0] * g[26];
        g[35] = c0[MXRYSROOTS*2+1] * g[29] + 2 * b0[1] * g[27];
}
static inline void _g0_lj_4d_1020(double *g, double *c0, double *cp,
                                  double *b0, double *b1,
                                  const double *r0, const double *rp)
{
        double rc[] = {r0[0]+c0[0], r0[0]+c0[1],
                       r0[1]+c0[MXRYSROOTS], r0[1]+c0[MXRYSROOTS+1],
                       r0[2]+c0[MXRYSROOTS*2], r0[2]+c0[MXRYSROOTS*2+1]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[4 ] = cp[0];
        g[5 ] = cp[1];
        g[6 ] = cp[0] * rc[0] + b0[0];
        g[7 ] = cp[1] * rc[1] + b0[1];
        g[8 ] = cp[0] * cp[0] + b1[0];
        g[9 ] = cp[1] * cp[1] + b1[1];
        g[10] = rc[0] * g[8] + 2 * b0[0] * cp[0];
        g[11] = rc[1] * g[9] + 2 * b0[1] * cp[1];
        g[24] = 1;
        g[25] = 1;
        g[26] = rc[2];
        g[27] = rc[3];
        g[28] = cp[MXRYSROOTS];
        g[29] = cp[MXRYSROOTS+1];
        g[30] = cp[MXRYSROOTS] * rc[2] + b0[0];
        g[31] = cp[MXRYSROOTS+1] * rc[3] + b0[1];
        g[32] = cp[MXRYSROOTS] * cp[MXRYSROOTS] + b1[0];
        g[33] = cp[MXRYSROOTS+1] * cp[MXRYSROOTS+1] + b1[1];
        g[34] = rc[2] * g[32] + 2 * b0[0] * cp[MXRYSROOTS];
        g[35] = rc[3] * g[33] + 2 * b0[1] * cp[MXRYSROOTS+1];
        //g[48] = w[0];
        //g[49] = w[1];
        g[50] = rc[4] * g[48];
        g[51] = rc[5] * g[49];
        g[52] = cp[MXRYSROOTS*2] * g[48];
        g[53] = cp[MXRYSROOTS*2+1] * g[49];
        g[54] =(cp[MXRYSROOTS*2] * rc[4] + b0[0]) * g[48];
        g[55] =(cp[MXRYSROOTS*2+1] * rc[5] + b0[1]) * g[49];
        g[56] =(cp[MXRYSROOTS*2] * cp[MXRYSROOTS*2] + b1[0]) * g[48];
        g[57] =(cp[MXRYSROOTS*2+1] * cp[MXRYSROOTS*2+1] + b1[1]) * g[49];
        g[58] = rc[4] * g[56] + 2 * b0[0] * g[52];
        g[59] = rc[5] * g[57] + 2 * b0[1] * g[53];
}
static inline void _g0_lj_4d_0111(double *g, double *c0, double *cp,
                                  double *b0, double *b1,
                                  const double *r0, const double *rp)
{
        double rc[] = {rp[0]+cp[0], rp[0]+cp[1],
                       rp[1]+cp[MXRYSROOTS], rp[1]+cp[MXRYSROOTS+1],
                       rp[2]+cp[MXRYSROOTS*2], rp[2]+cp[MXRYSROOTS*2+1]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[4 ] = cp[0];
        g[5 ] = cp[1];
        g[6 ] = cp[0] * rc[0] + b1[0];
        g[7 ] = cp[1] * rc[1] + b1[1];
        g[12] = c0[0];
        g[13] = c0[1];
        g[14] = c0[0] * rc[0] + b0[0];
        g[15] = c0[1] * rc[1] + b0[1];
        g[16] = c0[0] * cp[0] + b0[0];
        g[17] = c0[1] * cp[1] + b0[1];
        g[18] = c0[0] * g[6] + b0[0] *(rc[0] + cp[0]);
        g[19] = c0[1] * g[7] + b0[1] *(rc[1] + cp[1]);
        g[24] = 1;
        g[25] = 1;
        g[26] = rc[2];
        g[27] = rc[3];
        g[28] = cp[MXRYSROOTS];
        g[29] = cp[MXRYSROOTS+1];
        g[30] = cp[MXRYSROOTS] * rc[2] + b1[0];
        g[31] = cp[MXRYSROOTS+1] * rc[3] + b1[1];
        g[36] = c0[MXRYSROOTS];
        g[37] = c0[MXRYSROOTS+1];
        g[38] = c0[MXRYSROOTS] * rc[2] + b0[0];
        g[39] = c0[MXRYSROOTS+1] * rc[3] + b0[1];
        g[40] = c0[MXRYSROOTS] * cp[MXRYSROOTS] + b0[0];
        g[41] = c0[MXRYSROOTS+1] * cp[MXRYSROOTS+1] + b0[1];
        g[42] = c0[MXRYSROOTS] * g[30] + b0[0] *(rc[2] + cp[MXRYSROOTS]);
        g[43] = c0[MXRYSROOTS+1] * g[31] + b0[1] *(rc[3] + cp[MXRYSROOTS+1]);
        //g[48] = w[0];
        //g[49] = w[1];
        g[50] = rc[4] * g[48];
        g[51] = rc[5] * g[49];
        g[52] = cp[MXRYSROOTS*2] * g[48];
        g[53] = cp[MXRYSROOTS*2+1] * g[49];
        g[54] =(cp[MXRYSROOTS*2] * rc[4] + b1[0]) * g[48];
        g[55] =(cp[MXRYSROOTS*2+1] * rc[5] + b1[1]) * g[49];
        g[60] = c0[MXRYSROOTS*2] * g[48];
        g[61] = c0[MXRYSROOTS*2+1] * g[49];
        g[62] =(c0[MXRYSROOTS*2] * rc[4] + b0[0]) * g[48];
        g[63] =(c0[MXRYSROOTS*2+1] * rc[5] + b0[1]) * g[49];
        g[64] =(c0[MXRYSROOTS*2] * cp[MXRYSROOTS*2] + b0[0]) * g[48];
        g[65] =(c0[MXRYSROOTS*2+1] * cp[MXRYSROOTS*2+1] + b0[1]) * g[49];
        g[66] = c0[MXRYSROOTS*2] * g[54] + b0[0] *(g[50] + g[52]);
        g[67] = c0[MXRYSROOTS*2+1] * g[55] + b0[1] *(g[51] + g[53]);
}
static inline void _g0_lj_4d_1110(double *g, double *c0, double *cp, 
                                  double *b0, double *b1,
                                  const double *r0, const double *rp)
{
        double rc0[] = {r0[0]+c0[0], r0[0]+c0[1],
                        r0[1]+c0[MXRYSROOTS], r0[1]+c0[MXRYSROOTS+1],
                        r0[2]+c0[MXRYSROOTS*2], r0[2]+c0[MXRYSROOTS*2+1]};
        double rcp[] = {rp[0]+cp[0], rp[0]+cp[1],
                        rp[1]+cp[MXRYSROOTS], rp[1]+cp[MXRYSROOTS+1],
                        rp[2]+cp[MXRYSROOTS*2], rp[2]+cp[MXRYSROOTS*2+1]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc0[0];
        g[3 ] = rc0[1];
        g[4 ] = rcp[0];
        g[5 ] = rcp[1];
        g[6 ] = rcp[0] * rc0[0] + b0[0];
        g[7 ] = rcp[1] * rc0[1] + b0[1];
        g[8 ] = cp[0];
        g[9 ] = cp[1];
        g[10] = cp[0] * rc0[0] + b0[0];
        g[11] = cp[1] * rc0[1] + b0[1];
        g[12] = cp[0] * rcp[0] + b1[0];
        g[13] = cp[1] * rcp[1] + b1[1];
        g[14] = rc0[0] * g[12] + b0[0] *(rcp[0] + cp[0]);
        g[15] = rc0[1] * g[13] + b0[1] *(rcp[1] + cp[1]);
        g[48] = 1;
        g[49] = 1;
        g[50] = rc0[2];
        g[51] = rc0[3];
        g[52] = rcp[2];
        g[53] = rcp[3];
        g[54] = rcp[2] * rc0[2] + b0[0];
        g[55] = rcp[3] * rc0[3] + b0[1];
        g[56] = cp[MXRYSROOTS];
        g[57] = cp[MXRYSROOTS+1];
        g[58] = cp[MXRYSROOTS] * rc0[2] + b0[0];
        g[59] = cp[MXRYSROOTS+1] * rc0[3] + b0[1];
        g[60] = cp[MXRYSROOTS] * rcp[2] + b1[0];
        g[61] = cp[MXRYSROOTS+1] * rcp[3] + b1[1];
        g[62] = rc0[2] * g[60] + b0[0] *(rcp[2] + cp[MXRYSROOTS]);
        g[63] = rc0[3] * g[61] + b0[1] *(rcp[3] + cp[MXRYSROOTS+1]);
        //g[96 ] = w[0];
        //g[97 ] = w[1];
        g[98 ] = rc0[4] * g[96 ];
        g[99 ] = rc0[5] * g[97 ];
        g[100] = rcp[4] * g[96 ];
        g[101] = rcp[5] * g[97 ];
        g[102] =(rcp[4] * rc0[4] + b0[0])* g[96 ];
        g[103] =(rcp[5] * rc0[5] + b0[1])* g[97 ];
        g[104] = cp[MXRYSROOTS*2] * g[96 ];
        g[105] = cp[MXRYSROOTS*2+1] * g[97 ];
        g[106] =(cp[MXRYSROOTS*2] * rc0[4] + b0[0])* g[96 ];
        g[107] =(cp[MXRYSROOTS*2+1] * rc0[5] + b0[1])* g[97 ];
        g[108] =(cp[MXRYSROOTS*2] * rcp[4] + b1[0])* g[96 ];
        g[109] =(cp[MXRYSROOTS*2+1] * rcp[5] + b1[1])* g[97 ];
        g[110] = rc0[4] * g[108] + b0[0] *(g[100] + g[104]);
        g[111] = rc0[5] * g[109] + b0[1] *(g[101] + g[105]);
}
static inline void _g0_lj_4d_0201(double *g, double *c0, double *cp,
                                  double *b0, double *b1,
                                  const double *r0, const double *rp)
{
        double rc[] = {rp[0]+cp[0], rp[0]+cp[1],
                       rp[1]+cp[MXRYSROOTS], rp[1]+cp[MXRYSROOTS+1],
                       rp[2]+cp[MXRYSROOTS*2], rp[2]+cp[MXRYSROOTS*2+1]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[4 ] = rc[0] * rc[0] + b1[0];
        g[5 ] = rc[1] * rc[1] + b1[1];
        g[18] = c0[0];
        g[19] = c0[1];
        g[20] = rc[0] * c0[0] + b0[0];
        g[21] = rc[1] * c0[1] + b0[1];
        g[22] = c0[0] * g[4] + 2 * b0[0] * rc[0];
        g[23] = c0[1] * g[5] + 2 * b0[1] * rc[1];
        g[36] = 1;
        g[37] = 1;
        g[38] = rc[2];
        g[39] = rc[3];
        g[40] = rc[2] * rc[2] + b1[0];
        g[41] = rc[3] * rc[3] + b1[1];
        g[54] = c0[MXRYSROOTS];
        g[55] = c0[MXRYSROOTS+1];
        g[56] = rc[2] * c0[MXRYSROOTS] + b0[0];
        g[57] = rc[3] * c0[MXRYSROOTS+1] + b0[1];
        g[58] = c0[MXRYSROOTS] * g[40] + 2 * b0[0] * rc[2];
        g[59] = c0[MXRYSROOTS+1] * g[41] + 2 * b0[1] * rc[3];
        //g[72] = w[0];
        //g[73] = w[1];
        g[74] = rc[4] * g[72];
        g[75] = rc[5] * g[73];
        g[76] =(rc[4] * rc[4] + b1[0])* g[72];
        g[77] =(rc[5] * rc[5] + b1[1])* g[73];
        g[90] = c0[MXRYSROOTS*2] * g[72];
        g[91] = c0[MXRYSROOTS*2+1] * g[73];
        g[92] =(rc[4] * c0[MXRYSROOTS*2] + b0[0])* g[72];
        g[93] =(rc[5] * c0[MXRYSROOTS*2+1] + b0[1])* g[73];
        g[94] = c0[MXRYSROOTS*2] * g[76] + 2 * b0[0] * g[74];
        g[95] = c0[MXRYSROOTS*2+1] * g[77] + 2 * b0[1] * g[75];
}
static inline void _g0_lj_4d_1200(double *g, double *c0, double *cp, 
                                  double *b0, double *b1,
                                  const double *r0, const double *rp)
{
        double rc0[] = {r0[0]+c0[0], r0[0]+c0[1],
                        r0[1]+c0[MXRYSROOTS], r0[1]+c0[MXRYSROOTS+1],
                        r0[2]+c0[MXRYSROOTS*2], r0[2]+c0[MXRYSROOTS*2+1]};
        double rcp[] = {rp[0]+cp[0], rp[0]+cp[1],
                        rp[1]+cp[MXRYSROOTS], rp[1]+cp[MXRYSROOTS+1],
                        rp[2]+cp[MXRYSROOTS*2], rp[2]+cp[MXRYSROOTS*2+1]};
        g[0 ] = 1;
        g[1 ] = 1;
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
        g[72] = 1;
        g[73] = 1;
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
        //g[144] = w[0];
        //g[145] = w[1];
        g[146] = rc0[4] * g[144];
        g[147] = rc0[5] * g[145];
        g[148] = rcp[4] * g[144];
        g[149] = rcp[5] * g[145];
        g[150] =(rcp[4] * rc0[4] + b0[0])* g[144];
        g[151] =(rcp[5] * rc0[5] + b0[1])* g[145];
        g[152] =(rcp[4] * rcp[4] + b1[0])* g[144];
        g[153] =(rcp[5] * rcp[5] + b1[1])* g[145];
        g[154] = rc0[4] * g[152] + 2 * b0[0] * g[148];
        g[155] = rc0[5] * g[153] + 2 * b0[1] * g[149];
}
static inline void _g0_lj_4d_0012(double *g, double *c0, double *cp,
                                  double *b0, double *b1)
{
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = cp[0];
        g[3 ] = cp[1];
        g[4 ] = c0[0];
        g[5 ] = c0[1];
        g[6 ] = cp[0] * c0[0] + b0[0];
        g[7 ] = cp[1] * c0[1] + b0[1];
        g[8 ] = c0[0] * c0[0] + b1[0];
        g[9 ] = c0[1] * c0[1] + b1[1];
        g[10] = cp[0] * g[8] + 2 * b0[0] * c0[0];
        g[11] = cp[1] * g[9] + 2 * b0[1] * c0[1];
        g[12] = 1;
        g[13] = 1;
        g[14] = cp[MXRYSROOTS];
        g[15] = cp[MXRYSROOTS+1];
        g[16] = c0[MXRYSROOTS];
        g[17] = c0[MXRYSROOTS+1];
        g[18] = cp[MXRYSROOTS] * c0[MXRYSROOTS] + b0[0];
        g[19] = cp[MXRYSROOTS+1] * c0[MXRYSROOTS+1] + b0[1];
        g[20] = c0[MXRYSROOTS] * c0[MXRYSROOTS] + b1[0];
        g[21] = c0[MXRYSROOTS+1] * c0[MXRYSROOTS+1] + b1[1];
        g[22] = cp[MXRYSROOTS] * g[20] + 2 * b0[0] * c0[MXRYSROOTS];
        g[23] = cp[MXRYSROOTS+1] * g[21] + 2 * b0[1] * c0[MXRYSROOTS+1];
        //g[24] = w[0];
        //g[25] = w[1];
        g[26] = cp[MXRYSROOTS*2] * g[24];
        g[27] = cp[MXRYSROOTS*2+1] * g[25];
        g[28] = c0[MXRYSROOTS*2] * g[24];
        g[29] = c0[MXRYSROOTS*2+1] * g[25];
        g[30] =(cp[MXRYSROOTS*2] * c0[MXRYSROOTS*2] + b0[0]) * g[24];
        g[31] =(cp[MXRYSROOTS*2+1] * c0[MXRYSROOTS*2+1] + b0[1]) * g[25];
        g[32] =(c0[MXRYSROOTS*2] * c0[MXRYSROOTS*2] + b1[0]) * g[24];
        g[33] =(c0[MXRYSROOTS*2+1] * c0[MXRYSROOTS*2+1] + b1[1]) * g[25];
        g[34] = cp[MXRYSROOTS*2] * g[32] + 2 * b0[0] * g[28];
        g[35] = cp[MXRYSROOTS*2+1] * g[33] + 2 * b0[1] * g[29];
}
static inline void _g0_lj_4d_1011(double *g, double *c0, double *cp,
                                  double *b0, double *b1,
                                  const double *r0, const double *rp)
{
        double rc[] = {r0[0]+c0[0], r0[0]+c0[1],
                       r0[1]+c0[MXRYSROOTS], r0[1]+c0[MXRYSROOTS+1],
                       r0[2]+c0[MXRYSROOTS*2], r0[2]+c0[MXRYSROOTS*2+1]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[4 ] = cp[0];
        g[5 ] = cp[1];
        g[6 ] = cp[0] * rc[0] + b0[0];
        g[7 ] = cp[1] * rc[1] + b0[1];
        g[8 ] = c0[0];
        g[9 ] = c0[1];
        g[10] = c0[0] * rc[0] + b1[0];
        g[11] = c0[1] * rc[1] + b1[1];
        g[12] = c0[0] * cp[0] + b0[0];
        g[13] = c0[1] * cp[1] + b0[1];
        g[14] = cp[0] * g[10] + b0[0] *(rc[0] + c0[0]);
        g[15] = cp[1] * g[11] + b0[1] *(rc[1] + c0[1]);
        g[24] = 1;
        g[25] = 1;
        g[26] = rc[2];
        g[27] = rc[3];
        g[28] = cp[MXRYSROOTS];
        g[29] = cp[MXRYSROOTS+1];
        g[30] = cp[MXRYSROOTS] * rc[2] + b0[0];
        g[31] = cp[MXRYSROOTS+1] * rc[3] + b0[1];
        g[32] = c0[MXRYSROOTS];
        g[33] = c0[MXRYSROOTS+1];
        g[34] = c0[MXRYSROOTS] * rc[2] + b1[0];
        g[35] = c0[MXRYSROOTS+1] * rc[3] + b1[1];
        g[36] = c0[MXRYSROOTS] * cp[MXRYSROOTS] + b0[0];
        g[37] = c0[MXRYSROOTS+1] * cp[MXRYSROOTS+1] + b0[1];
        g[38] = cp[MXRYSROOTS] * g[34] + b0[0] *(rc[2] + c0[MXRYSROOTS]);
        g[39] = cp[MXRYSROOTS+1] * g[35] + b0[1] *(rc[3] + c0[MXRYSROOTS+1]);
        //g[48] = w[0];
        //g[49] = w[1];
        g[50] = rc[4] * g[48];
        g[51] = rc[5] * g[49];
        g[52] = cp[MXRYSROOTS*2] * g[48];
        g[53] = cp[MXRYSROOTS*2+1] * g[49];
        g[54] =(cp[MXRYSROOTS*2] * rc[4] + b0[0])* g[48];
        g[55] =(cp[MXRYSROOTS*2+1] * rc[5] + b0[1])* g[49];
        g[56] = c0[MXRYSROOTS*2] * g[48];
        g[57] = c0[MXRYSROOTS*2+1] * g[49];
        g[58] =(c0[MXRYSROOTS*2] * rc[4] + b1[0])* g[48];
        g[59] =(c0[MXRYSROOTS*2+1] * rc[5] + b1[1])* g[49];
        g[60] =(c0[MXRYSROOTS*2] * cp[MXRYSROOTS*2] + b0[0])* g[48];
        g[61] =(c0[MXRYSROOTS*2+1] * cp[MXRYSROOTS*2+1] + b0[1])* g[49];
        g[62] = cp[MXRYSROOTS*2] * g[58] + b0[0] *(g[50] + g[56]);
        g[63] = cp[MXRYSROOTS*2+1] * g[59] + b0[1] *(g[51] + g[57]);
}
static inline void _g0_lj_4d_2010(double *g, double *c0, double *cp,
                                  double *b0, double *b1,
                                  const double *r0, const double *rp)
{
        double rc[] = {r0[0]+c0[0], r0[0]+c0[1],
                       r0[1]+c0[MXRYSROOTS], r0[1]+c0[MXRYSROOTS+1],
                       r0[2]+c0[MXRYSROOTS*2], r0[2]+c0[MXRYSROOTS*2+1]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[4 ] = rc[0] * rc[0] + b1[0];
        g[5 ] = rc[1] * rc[1] + b1[1];
        g[6 ] = cp[0];
        g[7 ] = cp[1];
        g[8 ] = cp[0] * rc[0] + b0[0];
        g[9 ] = cp[1] * rc[1] + b0[1];
        g[10] = cp[0] * g[4] + 2 * b0[0] * rc[0];
        g[11] = cp[1] * g[5] + 2 * b0[1] * rc[1];
        g[36] = 1;
        g[37] = 1;
        g[38] = rc[2];
        g[39] = rc[3];
        g[40] = rc[2] * rc[2] + b1[0];
        g[41] = rc[3] * rc[3] + b1[1];
        g[42] = cp[MXRYSROOTS];
        g[43] = cp[MXRYSROOTS+1];
        g[44] = cp[MXRYSROOTS] * rc[2] + b0[0];
        g[45] = cp[MXRYSROOTS+1] * rc[3] + b0[1];
        g[46] = cp[MXRYSROOTS] * g[40] + 2 * b0[0] * rc[2];
        g[47] = cp[MXRYSROOTS+1] * g[41] + 2 * b0[1] * rc[3];
        //g[72] = w[0];
        //g[73] = w[1];
        g[74] = rc[4] * g[72];
        g[75] = rc[5] * g[73];
        g[76] =(rc[4] * rc[4] + b1[0]) * g[72];
        g[77] =(rc[5] * rc[5] + b1[1]) * g[73];
        g[78] = cp[MXRYSROOTS*2] * g[72];
        g[79] = cp[MXRYSROOTS*2+1] * g[73];
        g[80] =(cp[MXRYSROOTS*2] * rc[4] + b0[0]) * g[72];
        g[81] =(cp[MXRYSROOTS*2+1] * rc[5] + b0[1]) * g[73];
        g[82] = cp[MXRYSROOTS*2] * g[76] + 2 * b0[0] * g[74];
        g[83] = cp[MXRYSROOTS*2+1] * g[77] + 2 * b0[1] * g[75];
}
static inline void _g0_lj_4d_0102(double *g, double *c0, double *cp,
                                  double *b0, double *b1,
                                  const double *r0, const double *rp)
{
        double rc[] = {rp[0]+cp[0], rp[0]+cp[1],
                       rp[1]+cp[MXRYSROOTS], rp[1]+cp[MXRYSROOTS+1],
                       rp[2]+cp[MXRYSROOTS*2], rp[2]+cp[MXRYSROOTS*2+1]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[8 ] = c0[0];
        g[9 ] = c0[1];
        g[10] = rc[0] * c0[0] + b0[0];
        g[11] = rc[1] * c0[1] + b0[1];
        g[16] = c0[0] * c0[0] + b1[0];
        g[17] = c0[1] * c0[1] + b1[1];
        g[18] = rc[0] * g[16] + 2 * b0[0] * c0[0];
        g[19] = rc[1] * g[17] + 2 * b0[1] * c0[1];
        g[24] = 1;
        g[25] = 1;
        g[26] = rc[2];
        g[27] = rc[3];
        g[32] = c0[MXRYSROOTS];
        g[33] = c0[MXRYSROOTS+1];
        g[34] = rc[2] * c0[MXRYSROOTS] + b0[0];
        g[35] = rc[3] * c0[MXRYSROOTS+1] + b0[1];
        g[40] = c0[MXRYSROOTS] * c0[MXRYSROOTS] + b1[0];
        g[41] = c0[MXRYSROOTS+1] * c0[MXRYSROOTS+1] + b1[1];
        g[42] = rc[2] * g[40] + 2 * b0[0] * c0[MXRYSROOTS];
        g[43] = rc[3] * g[41] + 2 * b0[1] * c0[MXRYSROOTS+1];
        //g[48] = w[0];
        //g[49] = w[1];
        g[50] = rc[4] * g[48];
        g[51] = rc[5] * g[49];
        g[56] = c0[MXRYSROOTS*2] * g[48];
        g[57] = c0[MXRYSROOTS*2+1] * g[49];
        g[58] =(rc[4] * c0[MXRYSROOTS*2] + b0[0]) * g[48];
        g[59] =(rc[5] * c0[MXRYSROOTS*2+1] + b0[1]) * g[49];
        g[64] =(c0[MXRYSROOTS*2] * c0[MXRYSROOTS*2] + b1[0]) * g[48];
        g[65] =(c0[MXRYSROOTS*2+1] * c0[MXRYSROOTS*2+1] + b1[1]) * g[49];
        g[66] = rc[4] * g[64] + 2 * b0[0] * g[56];
        g[67] = rc[5] * g[65] + 2 * b0[1] * g[57];
}
static inline void _g0_lj_4d_1101(double *g, double *c0, double *cp, 
                                  double *b0, double *b1,
                                  const double *r0, const double *rp)
{
        double rc0[] = {r0[0]+c0[0], r0[0]+c0[1],
                        r0[1]+c0[MXRYSROOTS], r0[1]+c0[MXRYSROOTS+1],
                        r0[2]+c0[MXRYSROOTS*2], r0[2]+c0[MXRYSROOTS*2+1]};
        double rcp[] = {rp[0]+cp[0], rp[0]+cp[1],
                        rp[1]+cp[MXRYSROOTS], rp[1]+cp[MXRYSROOTS+1],
                        rp[2]+cp[MXRYSROOTS*2], rp[2]+cp[MXRYSROOTS*2+1]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc0[0];
        g[3 ] = rc0[1];
        g[4 ] = rcp[0];
        g[5 ] = rcp[1];
        g[6 ] = rcp[0] * rc0[0] + b0[0];
        g[7 ] = rcp[1] * rc0[1] + b0[1];
        g[16] = c0[0];
        g[17] = c0[1];
        g[18] = c0[0] * rc0[0] + b1[0];
        g[19] = c0[1] * rc0[1] + b1[1];
        g[20] = c0[0] * rcp[0] + b0[0];
        g[21] = c0[1] * rcp[1] + b0[1];
        g[22] = rcp[0] * g[18] + b0[0] *(rc0[0] + c0[0]);
        g[23] = rcp[1] * g[19] + b0[1] *(rc0[1] + c0[1]);
        g[48] = 1;
        g[49] = 1;
        g[50] = rc0[2];
        g[51] = rc0[3];
        g[52] = rcp[2];
        g[53] = rcp[3];
        g[54] = rcp[2] * rc0[2] + b0[0];
        g[55] = rcp[3] * rc0[3] + b0[1];
        g[64] = c0[MXRYSROOTS];
        g[65] = c0[MXRYSROOTS+1];
        g[66] = c0[MXRYSROOTS] * rc0[2] + b1[0];
        g[67] = c0[MXRYSROOTS+1] * rc0[3] + b1[1];
        g[68] = c0[MXRYSROOTS] * rcp[2] + b0[0];
        g[69] = c0[MXRYSROOTS+1] * rcp[3] + b0[1];
        g[70] = rcp[2] * g[66] + b0[0] *(rc0[2] + c0[MXRYSROOTS]);
        g[71] = rcp[3] * g[67] + b0[1] *(rc0[3] + c0[MXRYSROOTS+1]);
        //g[96 ] = w[0];
        //g[97 ] = w[1];
        g[98 ] = rc0[4] * g[96];
        g[99 ] = rc0[5] * g[97];
        g[100] = rcp[4] * g[96];
        g[101] = rcp[5] * g[97];
        g[102] =(rcp[4] * rc0[4] + b0[0])* g[96];
        g[103] =(rcp[5] * rc0[5] + b0[1])* g[97];
        g[112] = c0[MXRYSROOTS*2] * g[96];
        g[113] = c0[MXRYSROOTS*2+1] * g[97];
        g[114] =(c0[MXRYSROOTS*2] * rc0[4] + b1[0])* g[96];
        g[115] =(c0[MXRYSROOTS*2+1] * rc0[5] + b1[1])* g[97];
        g[116] =(c0[MXRYSROOTS*2] * rcp[4] + b0[0])* g[96];
        g[117] =(c0[MXRYSROOTS*2+1] * rcp[5] + b0[1])* g[97];
        g[118] = rcp[4] * g[114] + b0[0] *(g[98] + g[112]);
        g[119] = rcp[5] * g[115] + b0[1] *(g[99] + g[113]);
}
static inline void _g0_lj_4d_2100(double *g, double *c0, double *cp, 
                                  double *b0, double *b1,
                                  const double *r0, const double *rp)
{
/*
        //double rc0[] = {r0[0]+c0[0], r0[0]+c0[1],
        //                r0[1]+c0[MXRYSROOTS], r0[1]+c0[MXRYSROOTS+1],
        //                r0[2]+c0[MXRYSROOTS*2], r0[2]+c0[MXRYSROOTS*2+1]};
        //double rcp[] = {rp[0]+cp[0], rp[0]+cp[1],
        //                rp[1]+cp[MXRYSROOTS], rp[1]+cp[MXRYSROOTS+1],
        //                rp[2]+cp[MXRYSROOTS*2], rp[2]+cp[MXRYSROOTS*2+1]};
        __m128d rcx = _mm_set_pd(r0[0]+c0[1], r0[0]+c0[0]);
        __m128d rpx = _mm_set_pd(rp[0]+cp[1], rp[0]+cp[0]);
        __m128d rb0 = _mm_loadu_pd(&b0[0]);
        __m128d rb2 = _mm_add_pd(rb0, rb0);
        __m128d rb1 = _mm_loadu_pd(&b1[0]);
        __m128d m0, m1, m2, m3;
        g[0 ] = 1;
        g[1 ] = 1;
        //g[2 ] = rc0[0];
        //g[3 ] = rc0[1];
        _mm_store_pd(&g[2 ], rcx);
        //g[4 ] = rc0[0] * rc0[0] + b1[0];
        //g[5 ] = rc0[1] * rc0[1] + b1[1];
        m0 = _mm_mul_pd(rcx, rcx);
        m0 = _mm_add_pd(rb1, m0);
        _mm_store_pd(&g[4 ], m0);
        //g[6 ] = rcp[0];
        //g[7 ] = rcp[1];
        _mm_store_pd(&g[6 ], rpx);
        //g[8 ] = rcp[0] * rc0[0] + b0[0];
        //g[9 ] = rcp[1] * rc0[1] + b0[1];
        m1 = _mm_mul_pd(rpx, rcx);
        m1 = _mm_add_pd(rb0, m1);
        _mm_store_pd(&g[8 ], m1);
        //g[10] = rcp[0] * g[4] + 2 * b0[0] * rc0[0];
        //g[11] = rcp[1] * g[5] + 2 * b0[1] * rc0[1];
        m0 = _mm_mul_pd(m0, rpx);
        m1 = _mm_mul_pd(rb2, rcx);
        m1 = _mm_add_pd(m0, m1);
        _mm_store_pd(&g[10], m1);
        __m128d rcy = _mm_set_pd(r0[1]+c0[MXRYSROOTS+1], r0[1]+c0[MXRYSROOTS]);
        __m128d rpy = _mm_set_pd(rp[1]+cp[MXRYSROOTS+1], rp[1]+cp[MXRYSROOTS]);
        g[72] = 1;
        g[73] = 1;
        //g[74] = rc0[2];
        //g[75] = rc0[3];
        _mm_store_pd(&g[74], rcy);
        //g[76] = rc0[2] * rc0[2] + b1[0];
        //g[77] = rc0[3] * rc0[3] + b1[1];
        m0 = _mm_mul_pd(rcy, rcy);
        m0 = _mm_add_pd(rb1, m0);
        _mm_store_pd(&g[76], m0);
        //g[78] = rcp[2];
        //g[79] = rcp[3];
        _mm_store_pd(&g[78], rpy);
        //g[80] = rcp[2] * rc0[2] + b0[0];
        //g[81] = rcp[3] * rc0[3] + b0[1];
        m1 = _mm_mul_pd(rpy, rcy);
        m1 = _mm_add_pd(rb0, m1);
        _mm_store_pd(&g[80], m1);
        //g[82] = rcp[2] * g[76] + 2 * b0[0] * rc0[2];
        //g[83] = rcp[3] * g[77] + 2 * b0[1] * rc0[3];
        m0 = _mm_mul_pd(m0, rpy);
        m1 = _mm_mul_pd(rb2, rcy);
        m1 = _mm_add_pd(m0, m1);
        _mm_store_pd(&g[82], m1);
        __m128d rcz = _mm_set_pd(r0[2]+c0[MXRYSROOTS*2+1], r0[2]+c0[MXRYSROOTS*2]);
        __m128d rpz = _mm_set_pd(rp[2]+cp[MXRYSROOTS*2+1], rp[2]+cp[MXRYSROOTS*2]);
        //g[144] = w[0];
        //g[145] = w[1];
        //g[146] = rc0[4] * g[144];
        //g[147] = rc0[5] * g[145];
        m0 = _mm_load_pd(&g[144]);
        m2 = _mm_mul_pd(rcz, m0);
        _mm_store_pd(&g[146], m2);
        //g[148] =(rc0[4] * rc0[4] + b1[0])* g[144];
        //g[149] =(rc0[5] * rc0[5] + b1[1])* g[145];
        m3 = _mm_mul_pd(rcz, rcz);
        m3 = _mm_add_pd(rb1, m3);
        m3 = _mm_mul_pd(m3, m0);
        _mm_store_pd(&g[148], m3);
        //g[150] = rcp[4] * g[144];
        //g[151] = rcp[5] * g[145];
        m1 = _mm_mul_pd(rpz, m0);
        _mm_store_pd(&g[150], m1);
        //g[152] =(rcp[4] * rc0[4] + b0[0])* g[144];
        //g[153] =(rcp[5] * rc0[5] + b0[1])* g[145];
        m1 = _mm_mul_pd(rcz, rpz);
        m1 = _mm_add_pd(rb0, m1);
        m1 = _mm_mul_pd(m1, m0);
        _mm_store_pd(&g[152], m1);
        //g[154] = rcp[4] * g[148] + 2 * b0[0] * g[146];
        //g[155] = rcp[5] * g[149] + 2 * b0[1] * g[147];
        m3 = _mm_mul_pd(m3, rpz);
        m2 = _mm_mul_pd(m2, rb2);
        m2 = _mm_add_pd(m2, m3);
        _mm_store_pd(&g[154], m2);
*/
        double rc0[] = {r0[0]+c0[0], r0[0]+c0[1],
                        r0[1]+c0[MXRYSROOTS], r0[1]+c0[MXRYSROOTS+1],
                        r0[2]+c0[MXRYSROOTS*2], r0[2]+c0[MXRYSROOTS*2+1]};
        double rcp[] = {rp[0]+cp[0], rp[0]+cp[1],
                        rp[1]+cp[MXRYSROOTS], rp[1]+cp[MXRYSROOTS+1],
                        rp[2]+cp[MXRYSROOTS*2], rp[2]+cp[MXRYSROOTS*2+1]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc0[0];
        g[3 ] = rc0[1];
        g[4 ] = rc0[0] * rc0[0] + b1[0];
        g[5 ] = rc0[1] * rc0[1] + b1[1];
        g[6 ] = rcp[0];
        g[7 ] = rcp[1];
        g[8 ] = rcp[0] * rc0[0] + b0[0];
        g[9 ] = rcp[1] * rc0[1] + b0[1];
        g[10] = rcp[0] * g[4] + 2 * b0[0] * rc0[0];
        g[11] = rcp[1] * g[5] + 2 * b0[1] * rc0[1];
        g[72] = 1;
        g[73] = 1;
        g[74] = rc0[2];
        g[75] = rc0[3];
        g[76] = rc0[2] * rc0[2] + b1[0];
        g[77] = rc0[3] * rc0[3] + b1[1];
        g[78] = rcp[2];
        g[79] = rcp[3];
        g[80] = rcp[2] * rc0[2] + b0[0];
        g[81] = rcp[3] * rc0[3] + b0[1];
        g[82] = rcp[2] * g[76] + 2 * b0[0] * rc0[2];
        g[83] = rcp[3] * g[77] + 2 * b0[1] * rc0[3];
        //g[144] = w[0];
        //g[145] = w[1];
        g[146] = rc0[4] * g[144];
        g[147] = rc0[5] * g[145];
        g[148] =(rc0[4] * rc0[4] + b1[0])* g[144];
        g[149] =(rc0[5] * rc0[5] + b1[1])* g[145];
        g[150] = rcp[4] * g[144];
        g[151] = rcp[5] * g[145];
        g[152] =(rcp[4] * rc0[4] + b0[0])* g[144];
        g[153] =(rcp[5] * rc0[5] + b0[1])* g[145];
        g[154] = rcp[4] * g[148] + 2 * b0[0] * g[146];
        g[155] = rcp[5] * g[149] + 2 * b0[1] * g[147];
}
/************** end special g0_4d results *************/



void CINTg0_2e_lj2d4d(double *g, const CINTEnvVars *envs,struct _BC *bc)
{
        const FINT nmax = envs->li_ceil + envs->lj_ceil;
        const FINT mmax = envs->lk_ceil + envs->ll_ceil;
        switch (nmax) {
                case 0: switch(mmax) {
                        case 0: goto _g0_4d_default; // ssss
                        case 1: switch (envs->lk_ceil) {
                                case 0: _g0_lj_4d_0001(g, bc->c0p, envs->rkrl); goto normal_end;
                                case 1: _g0_lj_4d_1000(g, bc->c0p, envs->rkrl); goto normal_end;
                                default: goto error; }
                        case 2: switch (envs->lk_ceil) {
                                case 0: _g0_lj_4d_0002(g, bc->c0p, bc->b01, envs->rkrl); goto normal_end;
                                case 1: _g0_lj_4d_1001(g, bc->c0p, bc->b01, envs->rkrl); goto normal_end;
                                case 2: _g0_lj_4d_2000(g, bc->c0p, bc->b01, envs->rkrl); goto normal_end;
                                default: goto error; }
                        case 3: switch (envs->lk_ceil) {
                                case 0: _g0_lj_4d_0003(g, bc->c0p, bc->b01, envs->rkrl); goto normal_end;
                                case 1: _g0_lj_4d_1002(g, bc->c0p, bc->b01, envs->rkrl); goto normal_end;
                                case 2: _g0_lj_4d_2001(g, bc->c0p, bc->b01, envs->rkrl); goto normal_end;
                                case 3: _g0_lj_4d_3000(g, bc->c0p, bc->b01, envs->rkrl); goto normal_end;
                                default: goto error; }
                        default: goto _g0_4d_default; }
                case 1: switch(mmax) {
                        case 0: switch (envs->li_ceil) {
                                case 0: _g0_lj_4d_0001(g, bc->c00, envs->rirj); goto normal_end;
                                case 1: _g0_lj_4d_1000(g, bc->c00, envs->rirj); goto normal_end;
                                default: goto error; }
                        case 1: switch (envs->lk_ceil) {
                                case 0: switch (envs->li_ceil) {
                                        case 0: _g0_lj_4d_0011(g, bc->c00, bc->c0p, bc->b00, envs->rirj, envs->rkrl); goto normal_end;
                                        case 1: _g0_lj_4d_1010(g, bc->c00, bc->c0p, bc->b00, envs->rirj, envs->rkrl); goto normal_end;
                                        default: goto error; }
                                case 1: switch (envs->li_ceil) {
                                        case 0: _g0_lj_4d_0101(g, bc->c00, bc->c0p, bc->b00, envs->rirj, envs->rkrl); goto normal_end;
                                        case 1: _g0_lj_4d_1100(g, bc->c00, bc->c0p, bc->b00, envs->rirj, envs->rkrl); goto normal_end;
                                        default: goto error; }
                                default: goto error; }
                        case 2: switch (envs->lk_ceil) {
                                case 0: switch (envs->li_ceil) {
                                        case 0: _g0_lj_4d_0021(g, bc->c00, bc->c0p, bc->b00, bc->b01); goto normal_end;
                                        case 1: _g0_lj_4d_1020(g, bc->c00, bc->c0p, bc->b00, bc->b01, envs->rirj, envs->rkrl); goto normal_end;
                                        default: goto error; }
                                case 1: switch (envs->li_ceil) {
                                        case 0: _g0_lj_4d_0111(g, bc->c00, bc->c0p, bc->b00, bc->b01, envs->rirj, envs->rkrl); goto normal_end;
                                        case 1: _g0_lj_4d_1110(g, bc->c00, bc->c0p, bc->b00, bc->b01, envs->rirj, envs->rkrl); goto normal_end;
                                        default: goto error; }
                                case 2: switch (envs->li_ceil) {
                                        case 0: _g0_lj_4d_0201(g, bc->c00, bc->c0p, bc->b00, bc->b01, envs->rirj, envs->rkrl); goto normal_end;
                                        case 1: _g0_lj_4d_1200(g, bc->c00, bc->c0p, bc->b00, bc->b01, envs->rirj, envs->rkrl); goto normal_end;
                                        default: goto error; }
                                default: goto error; }
                        default: goto _g0_4d_default; }
                case 2: switch(mmax) {
                        case 0: switch (envs->li_ceil) {
                                case 0: _g0_lj_4d_0002(g, bc->c00, bc->b10, envs->rirj); goto normal_end;
                                case 1: _g0_lj_4d_1001(g, bc->c00, bc->b10, envs->rirj); goto normal_end;
                                case 2: _g0_lj_4d_2000(g, bc->c00, bc->b10, envs->rirj); goto normal_end;
                                default: goto error; }
                        case 1: switch (envs->lk_ceil) {
                                case 0: switch (envs->li_ceil) {
                                        case 0: _g0_lj_4d_0012(g, bc->c00, bc->c0p, bc->b00, bc->b10); goto normal_end;
                                        case 1: _g0_lj_4d_1011(g, bc->c00, bc->c0p, bc->b00, bc->b10, envs->rirj, envs->rkrl); goto normal_end;
                                        case 2: _g0_lj_4d_2010(g, bc->c00, bc->c0p, bc->b00, bc->b10, envs->rirj, envs->rkrl); goto normal_end;
                                        default: goto error; }
                                case 1: switch (envs->li_ceil) {
                                        case 0: _g0_lj_4d_0102(g, bc->c00, bc->c0p, bc->b00, bc->b10, envs->rirj, envs->rkrl); goto normal_end;
                                        case 1: _g0_lj_4d_1101(g, bc->c00, bc->c0p, bc->b00, bc->b10, envs->rirj, envs->rkrl); goto normal_end;
                                        case 2: _g0_lj_4d_2100(g, bc->c00, bc->c0p, bc->b00, bc->b10, envs->rirj, envs->rkrl); goto normal_end;
                                        default: goto error; }
                                default: goto error; }
                        default: goto _g0_4d_default; }
                case 3: switch(mmax) {
                        case 0: switch (envs->li_ceil) {
                                case 0: _g0_lj_4d_0003(g, bc->c00, bc->b10, envs->rirj); goto normal_end;
                                case 1: _g0_lj_4d_1002(g, bc->c00, bc->b10, envs->rirj); goto normal_end;
                                case 2: _g0_lj_4d_2001(g, bc->c00, bc->b10, envs->rirj); goto normal_end;
                                case 3: _g0_lj_4d_3000(g, bc->c00, bc->b10, envs->rirj); goto normal_end;
                                default: goto error; }
                        default: goto _g0_4d_default; }
                default:
_g0_4d_default:
                        CINTg0_2e_2d(g, bc, envs);
                        CINTg0_lj2d_4d(g, envs);
        }
normal_end:
        return;
error:
        printf("Dimension error for CINTg0_2e_lj2d4d: iklj = %d %d %d %d",
               (int)envs->li_ceil, (int)envs->lk_ceil,
               (int)envs->ll_ceil, (int)envs->lj_ceil);
        exit(1);
}


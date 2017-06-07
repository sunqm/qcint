/*
 *      Fast exp(x) computation (with SSE2 optimizations).
 *
 * Copyright (c) 2010, Naoaki Okazaki
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the names of the authors nor the names of its contributors
 *       may be used to endorse or promote products derived from this
 *       software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * Modified from   http://www.chokkan.org/blog/archives/352
 *                 http://www.chokkan.org/software/dist/fastexp.c
 */

#include <math.h>
#include "simd.h"

typedef union {
    double d;
    unsigned short s[4];
} ieee754;

static const double MAXLOG =  7.08396418532264106224E2;     /* log 2**1022 */
static const double MINLOG = -7.08396418532264106224E2;     /* log 2**-1022 */
static const double LOG2E  =  1.4426950408889634073599;     /* 1/log(2) */
//static const double INFINITY = 1.79769313486231570815E308;
static const double C1 = 6.93145751953125E-1;
static const double C2 = 1.42860682030941723212E-6;

/*
 * Note this function cannot handle double precision overflow.
 * It has huge error when |x| > 709
 */
void CINTexp_cephes(double *expx, double *x)
{
        int i, n;
        ALIGNMM double a[SIMDD];
        __MD xx, px, qx;
        __MD rx = MM_LOAD(x);
        ALIGNMM ieee754 u[SIMDD] = {0,};

        /* n = round(x / log 2) */
        MM_STORE(a, MM_SET1(LOG2E) * rx + MM_SET1(0.5));
        for (i = 0; i < SIMDD; i++) {
                n  = (int)a[i];
                n -= (a[i] < 0);
                a[i] = n;
                /* Build 2^n in double. */
                (u+i)->s[3] = (unsigned short)(((n+1023) << 4) & 0x7FF0);
        }

        /* x -= n * log2 */
        //px = (double)n;
        px = MM_LOAD(a);
        rx -= px * MM_SET1(C1);
        rx -= px * MM_SET1(C2);
        xx = rx * rx;

        /* px = x * P(x**2). */
        px = MM_SET1(1.26177193074810590878E-4) * xx + MM_SET1(3.02994407707441961300E-2);
        px = px * xx + MM_SET1(9.99999999999999999910E-1);
        px *= rx;

        /* Evaluate Q(x**2). */
        qx = MM_SET1(3.00198505138664455042E-6) * xx + MM_SET1(2.52448340349684104192E-3);
        qx = qx * xx + MM_SET1(2.27265548208155028766E-1);
        qx = qx * xx + MM_SET1(2.00000000000000000009E0);

        /* e**x = 1 + 2x P(x**2)/( Q(x**2) - P(x**2) ) */
        rx = px / (qx - px);
        rx = MM_SET1(1.0) + MM_SET1(2.0) * rx;
        MM_STORE(expx, rx * MM_LOAD((double *)u));

        for (i = 0; i < SIMDD; i++) {
                if (x[i] < -740) {
                        expx[i] = 0;
                }
        }
}


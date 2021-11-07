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

#if !defined HAVE_DEFINED_ALIGN
#define HAVE_DEFINED_ALIGN

#include <stdint.h>
#include <immintrin.h>
#include <mm_malloc.h>
#include "cint.h"

#ifdef __AVX512F__
#define __MD            __m512d
#define MM_LOAD         _mm512_load_pd
#define MM_LOADU        _mm512_loadu_pd
#define MM_MUL          _mm512_mul_pd
#define MM_ADD          _mm512_add_pd
#define MM_SUB          _mm512_sub_pd
#define MM_DIV          _mm512_div_pd
#define MM_SQRT         _mm512_sqrt_pd
#define MM_SET1         _mm512_set1_pd
#define MM_STORE        _mm512_store_pd
#define MM_STOREU       _mm512_storeu_pd
#define MM_GATHER(base_addr, vindex, scale)     _mm512_i32gather_pd(vindex, base_addr, scale)
#define MM_FMA          _mm512_fmadd_pd
#define MM_FNMA         _mm512_fnmadd_pd
#define MM_CMP          _mm512_cmp_pd_mask
//#define MM_EXPN(y,x,rx) _mm512_store_pd(y, _mm512_exp_pd(rx))
#define MM_EXPN(y,x,rx) y[0] = exp(-x[0]); y[1] = exp(-x[1]); y[2] = exp(-x[2]); y[3] = exp(-x[3]); \
                        y[4] = exp(-x[4]); y[5] = exp(-x[5]); y[6] = exp(-x[6]); y[7] = exp(-x[7])

#elif __AVX__
#define __MD            __m256d
#define MM_LOAD         _mm256_load_pd
#define MM_LOADU        _mm256_loadu_pd
#define MM_MUL          _mm256_mul_pd
#define MM_ADD          _mm256_add_pd
#define MM_SUB          _mm256_sub_pd
#define MM_DIV          _mm256_div_pd
#define MM_SQRT         _mm256_sqrt_pd
#define MM_SET1         _mm256_set1_pd
#define MM_STORE        _mm256_store_pd
#define MM_STOREU       _mm256_storeu_pd
#define MM_GATHER       _mm256_i32gather_pd
#ifdef __FMA__
#define MM_FMA          _mm256_fmadd_pd
#define MM_FNMA         _mm256_fnmadd_pd
#else
#define MM_FMA(a,b,c)   _mm256_add_pd(_mm256_mul_pd(a, b), c)
#define MM_FNMA(a,b,c)  _mm256_sub_pd(c, _mm256_mul_pd(a, b))
#endif
#define MM_CMP(a,b,c)   _mm256_movemask_pd(_mm256_cmp_pd(a,b,c))
#define MM_EXPN(y,x,rx) y[0] = exp(-x[0]); y[1] = exp(-x[1]); y[2] = exp(-x[2]); y[3] = exp(-x[3])

#elif __SSE3__
#define __MD            __m128d
#define MM_LOAD         _mm_load_pd
#define MM_LOADU        _mm_loadu_pd
#define MM_MUL          _mm_mul_pd
#define MM_ADD          _mm_add_pd
#define MM_SUB          _mm_sub_pd
#define MM_DIV          _mm_div_pd
#define MM_SQRT         _mm_sqrt_pd
#define MM_SET1         _mm_set1_pd
#define MM_STORE        _mm_store_pd
#define MM_STOREU       _mm_storeu_pd
#define MM_GATHER       _mm_i32gather_pd
#ifdef __FMA__
#define MM_FMA          _mm_fmadd_pd
#define MM_FNMA         _mm_fnmadd_pd
#else
#define MM_FMA(a,b,c)   _mm_add_pd(_mm_mul_pd(a, b), c)
#define MM_FNMA(a,b,c)  _mm_sub_pd(c, _mm_mul_pd(a, b))
#endif
#ifdef __SSE4_2__
#define MM_CMP(a,b,c)   _mm_movemask_pd(_mm_cmp_pd(a,b,c))
#endif
#define MM_EXPN(y,x,rx) y[0] = exp(-x[0]); y[1] = exp(-x[1])

#endif

#define ALIGN_UP(x, align)  (((uintptr_t)(x) + (uintptr_t)(align) - 1) & (-(uintptr_t)(align)))

#define MALLOC_DOUBLE_WITHCACHE(var, n, cache) \
        double *var; \
        int __##var##len = 0; \
        if (cache != NULL) { \
                __##var##len = n; \
                var = (void *)ALIGN_UP(cache, sizeof(double)*SIMDD); \
                cache = var + __##var##len; \
        } else { \
                var = (double *)_mm_malloc(sizeof(double)*(n), sizeof(double)*SIMDD); \
        }

// Simple malloc(sizeof(type(data)))
#define MALLOC_DATA_INSTACK(var, n) \
        var = (void *)ALIGN_UP(cache, sizeof(double)); \
        cache = (double *)(var + (n));

// malloc(sizeof(double)) with alignment
#define MALLOC_ALIGNED_DOUBLE_INSTACK(var, n) \
        var = (void *)ALIGN_UP(cache, sizeof(double)*SIMDD); \
        cache = (double *)(var + (n));
#define MALLOC_INSTACK(var, n)  MALLOC_ALIGNED_DOUBLE_INSTACK(var, n)

#endif

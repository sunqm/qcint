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

#include "cint.h"
#include "simd.h"

void CINTinit_int1e_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                            int *atm, int natm, int *bas, int nbas, double *env);
void CINTinit_int3c1e_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                              int *atm, int natm, int *bas, int nbas, double *env);

void CINTg2c_index_xyz(int *idx, CINTEnvVars *envs);

void CINTg1e_ovlp(double *g, CINTEnvVars *envs, int count);
void CINTg1e_nuc(double *g, CINTEnvVars *envs, int count, int nuc_id);
void CINTg3c1e_ovlp(double *g, CINTEnvVars *envs, int count);
void CINTg3c1e_nuc(double *g, CINTEnvVars *envs, int count, int nuc_id);
double CINTnuc_mod(double aij, int nuc_id, int *atm, double *env);

void CINTnabla1i_1e(double *f, double *g,
                    int li, int lj, int lk, CINTEnvVars *envs);
void CINTnabla1j_1e(double *f, double *g,
                    int li, int lj, int lk, CINTEnvVars *envs);
void CINTnabla1k_1e(double *f, double *g,
                    int li, int lj, int lk, CINTEnvVars *envs);
void CINTx1i_1e(double *f, double *g, double *ri,
                int li, int lj, int lk, CINTEnvVars *envs);
void CINTx1j_1e(double *f, double *g, double *rj,
                int li, int lj, int lk, CINTEnvVars *envs);
void CINTx1k_1e(double *f, double *g, double *rk,
                int li, int lj, int lk, CINTEnvVars *envs);

void CINTprim_to_ctr(double *gc, int nf, double *gp,
                     int inc, int nprim, int nctr, double *pcoeff);
void CINTprim_to_ctr_0(double *gc, double *gp, double *coeff, size_t nf,
                       int nprim, int nctr, int non0ctr, int *sortedidx);
void CINTprim_to_ctr_1(double *gc, double *gp, double *coeff, size_t nf,
                       int nprim, int nctr, int non0ctr, int *sortedidx);
void CINTiprim_to_ctr_0(double *gc, double *gp, double *coeff, size_t nf,
                        int nprim, int nctr, int non0ctr, int *sortedidx);
void CINTiprim_to_ctr_1(double *gc, double *gp, double *coeff, size_t nf,
                        int nprim, int nctr, int non0ctr, int *sortedidx);
void CINTsort_gout(double *sout, double *gout, int nf, int count);

double CINTcommon_fac_sp(int l);

#define G1E_D_I(f, g, li, lj, lk)   CINTnabla1i_1e(f, g, li, lj, lk, envs)
#define G1E_D_J(f, g, li, lj, lk)   CINTnabla1j_1e(f, g, li, lj, lk, envs)
#define G1E_D_K(f, g, li, lj, lk)   CINTnabla1k_1e(f, g, li, lj, lk, envs)
/* r-R_0, R_0 is (0,0,0) */
#define G1E_R0I(f, g, li, lj, lk)   CINTx1i_1e(f, g, envs->ri, li, lj, lk, envs)
#define G1E_R0J(f, g, li, lj, lk)   CINTx1j_1e(f, g, envs->rj, li, lj, lk, envs)
#define G1E_R0K(f, g, li, lj, lk)   CINTx1k_1e(f, g, envs->rk, li, lj, lk, envs)
/* r-R_C, R_C is common origin */
#define G1E_RCI(f, g, li, lj, lk)   CINTx1i_1e(f, g, dri, li, lj, lk, envs)
#define G1E_RCJ(f, g, li, lj, lk)   CINTx1j_1e(f, g, drj, li, lj, lk, envs)
#define G1E_RCK(f, g, li, lj, lk)   CINTx1k_1e(f, g, drk, li, lj, lk, envs)
/* origin from center of each basis
 * x1[ij]_1e(f, g, ng, li, lj, 0d0) */
#define G1E_R_I(f, g, li, lj, lk)   f = g + envs->g_stride_i * SIMDD
#define G1E_R_J(f, g, li, lj, lk)   f = g + envs->g_stride_j * SIMDD
#define G1E_R_K(f, g, li, lj, lk)   f = g + envs->g_stride_k * SIMDD


#if (SIMDD == 8)

#define DECLARE_GOUT \
        __m256i vindex = _mm256_set_epi32( \
                nfc*7, nfc*6, nfc*5, nfc*4, nfc*3, nfc*2, nfc*1,    0);

#define GOUT_SCATTER(gout, n, r0) _mm512_i32scatter_pd(gout+n, vindex, r0, 8)

#elif (SIMDD == 4)

#define DECLARE_GOUT \
        ALIGNMM double gc[SIMDD*SIMDD]; \
        double *RESTRICT gout1 = gout + nfc*1; \
        double *RESTRICT gout2 = gout + nfc*2; \
        double *RESTRICT gout3 = gout + nfc*3;

#define GOUT_SCATTER(addr, n, r0) \
        _mm256_store_pd(gc, r0); \
        gout [(n)] = gc[0]; \
        gout1[(n)] = gc[1]; \
        gout2[(n)] = gc[2]; \
        gout3[(n)] = gc[3];

#else // SSE3

#define DECLARE_GOUT \
        ALIGNMM double gc[SIMDD*SIMDD]; \
        double *RESTRICT gout1 = gout + nfc*1; \

#define GOUT_SCATTER(addr, n, r0) \
        _mm_store_pd(gc, r0); \
        gout [(n)] = gc[0]; \
        gout1[(n)] = gc[1];
#endif

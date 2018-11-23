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

#include "simd.h"
#if !defined HAVE_DEFINED_CINTENVVARS_H
#define HAVE_DEFINED_CINTENVVARS_H
// ref to CINTinit_int1e_EnvVars, CINTinit_int2e_EnvVars
typedef struct {
        int *atm;
        int *bas;
        double *env;
        int *shls;
        int natm;
        int nbas;

        int i_l;
        int j_l;
        int k_l;
        int l_l;
        int nfi;  // number of cartesion components
        int nfj;
        int nfk;
        int nfl;
        int nf;  // = nfi*nfj*nfk*nfl;
        int _padding;
        int x_ctr[4];

        int gbits;
        int ncomp_e1; // = 1 if spin free, = 4 when spin included, it
        int ncomp_e2; // corresponds to POSX,POSY,POSZ,POS1, see cint_const.h
        int ncomp_tensor; // e.g. = 3 for gradients

        /* values may diff based on the g0_2d4d algorithm */
        int li_ceil; // power of x, == i_l if nabla is involved, otherwise == i_l
        int lj_ceil;
        int lk_ceil;
        int ll_ceil;
        int g_stride_i; // nrys_roots * shift of (i++,k,l,j)
        int g_stride_k; // nrys_roots * shift of (i,k++,l,j)
        int g_stride_l; // nrys_roots * shift of (i,k,l++,j)
        int g_stride_j; // nrys_roots * shift of (i,k,l,j++)
        int nrys_roots;
        int g_size;  // ref to cint2e.c g = malloc(sizeof(double)*g_size)

        int g2d_ijmax;
        int g2d_klmax;
        double common_factor;
        double _padding1;
        double rirj[3]; // diff by sign in different g0_2d4d algorithm
        double rkrl[3];
        double *rx_in_rijrx;
        double *rx_in_rklrx;

        double *ri;
        double *rj;
        double *rk;
        double *rl;

        void (*f_g0_2e)();
        void (*f_g0_2e_simd1)();
        void (*f_g0_2d4d)();
        void (*f_g0_2d4d_simd1)();
        void (*f_gout)();
        void (*f_gout_simd1)();

        /* values are assigned during calculation */
        ALIGNMM double ai[SIMDD];
        ALIGNMM double aj[SIMDD];
        ALIGNMM double ak[SIMDD];
        ALIGNMM double al[SIMDD];
        ALIGNMM double fac[SIMDD];
        ALIGNMM double rij[SIMDD*3];
        ALIGNMM double rkl[SIMDD*3];
} CINTEnvVars;
#endif

#define RYS_ROOTS       6
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
void CINTprim_to_ctr_0(double *gc, double *gp, double *coeff, int nf,
                       int nprim, int nctr, int non0ctr, int *sortedidx);
void CINTprim_to_ctr_1(double *gc, double *gp, double *coeff, int nf,
                       int nprim, int nctr, int non0ctr, int *sortedidx);
void CINTiprim_to_ctr_0(double *gc, double *gp, double *coeff, int nf,
                        int nprim, int nctr, int non0ctr, int *sortedidx);
void CINTiprim_to_ctr_1(double *gc, double *gp, double *coeff, int nf,
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

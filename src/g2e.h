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

#include "config.h"
#include "g1e.h"

#ifndef HAVE_RYS2E
#define HAVE_RYS2E
typedef struct {
        ALIGNMM double c00x[SIMDD*MXRYSROOTS];
        ALIGNMM double c00y[SIMDD*MXRYSROOTS];
        ALIGNMM double c00z[SIMDD*MXRYSROOTS];
        ALIGNMM double c0px[SIMDD*MXRYSROOTS];
        ALIGNMM double c0py[SIMDD*MXRYSROOTS];
        ALIGNMM double c0pz[SIMDD*MXRYSROOTS];
        ALIGNMM double b01[SIMDD*MXRYSROOTS];
        ALIGNMM double b00[SIMDD*MXRYSROOTS];
        ALIGNMM double b10[SIMDD*MXRYSROOTS];
        ALIGNMM double u[MXRYSROOTS*SIMDD];
        ALIGNMM double w[MXRYSROOTS*SIMDD];
} Rys2eT;
#endif


void CINTg4c_index_xyz(int *idx, CINTEnvVars *envs);
void CINTinit_int2e_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                            int *atm, int natm, int *bas, int nbas, double *env);
void CINTinit_int3c2e_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                              int *atm, int natm, int *bas, int nbas, double *env);
void CINTinit_int2c2e_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                              int *atm, int natm, int *bas, int nbas, double *env);
void CINTinit_int2e_coulerf_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                                    int *atm, int natm, int *bas, int nbas, double *env);
void CINTinit_int2e_stg_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                                int *atm, int natm, int *bas, int nbas, double *env);
void CINTinit_int2e_yp_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                               int *atm, int natm, int *bas, int nbas, double *env);

int CINTg0_2e(double *g, double *cutoff,
              Rys2eT *bc, CINTEnvVars *envs, int count);
int CINTg0_2e_simd1(double *g, double *cutoff,
                    Rys2eT *bc, CINTEnvVars *envs, int idsimd);
void CINTg0_2e_2d(double *g, Rys2eT *bc, CINTEnvVars *envs);
void CINTg0_2e_2d_simd1(double *g, Rys2eT *bc, CINTEnvVars *envs);
void CINTg0_2e_2d4d_unrolled(double *g, Rys2eT *bc, CINTEnvVars *envs);
void CINTg0_2e_2d4d_unrolled_simd1(double *g, Rys2eT *bc, CINTEnvVars *envs);
void CINTg0_2e_lj2d4d(double *g, Rys2eT *bc, CINTEnvVars *envs);
void CINTg0_2e_kj2d4d(double *g, Rys2eT *bc, CINTEnvVars *envs);
void CINTg0_2e_il2d4d(double *g, Rys2eT *bc, CINTEnvVars *envs);
void CINTg0_2e_ik2d4d(double *g, Rys2eT *bc, CINTEnvVars *envs);
void CINTg0_2e_lj2d4d_simd1(double *g, Rys2eT *bc, CINTEnvVars *envs);
void CINTg0_2e_kj2d4d_simd1(double *g, Rys2eT *bc, CINTEnvVars *envs);
void CINTg0_2e_il2d4d_simd1(double *g, Rys2eT *bc, CINTEnvVars *envs);
void CINTg0_2e_ik2d4d_simd1(double *g, Rys2eT *bc, CINTEnvVars *envs);

void CINTinit_int2e_stg_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                                int *atm, int natm, int *bas, int nbas, double *env);
void CINTinit_int2e_yp_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                               int *atm, int natm, int *bas, int nbas, double *env);
void CINTinit_int2e_coulerf_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                                    int *atm, int natm, int *bas, int nbas, double *env);

void CINTnabla1i_2e(double *f, double *g,
                    int li, int lj, int lk, int ll, CINTEnvVars *envs);
void CINTnabla1j_2e(double *f, double *g,
                    int li, int lj, int lk, int ll, CINTEnvVars *envs);
void CINTnabla1k_2e(double *f, double *g,
                    int li, int lj, int lk, int ll, CINTEnvVars *envs);
void CINTnabla1l_2e(double *f, double *g,
                    int li, int lj, int lk, int ll, CINTEnvVars *envs);

void CINTx1i_2e(double *f, double *g, double *ri,
                int li, int lj, int lk, int ll, CINTEnvVars *envs);
void CINTx1j_2e(double *f, double *g, double *rj,
                int li, int lj, int lk, int ll, CINTEnvVars *envs);
void CINTx1k_2e(double *f, double *g, double *rk,
                int li, int lj, int lk, int ll, CINTEnvVars *envs);
void CINTx1l_2e(double *f, double *g, double *rl,
                int li, int lj, int lk, int ll, CINTEnvVars *envs);

void CINTnabla1i_2e_simd1(double *f, double *g,
                          int li, int lj, int lk, int ll, CINTEnvVars *envs);
void CINTnabla1j_2e_simd1(double *f, double *g,
                          int li, int lj, int lk, int ll, CINTEnvVars *envs);
void CINTnabla1k_2e_simd1(double *f, double *g,
                          int li, int lj, int lk, int ll, CINTEnvVars *envs);
void CINTnabla1l_2e_simd1(double *f, double *g,
                          int li, int lj, int lk, int ll, CINTEnvVars *envs);

void CINTx1i_2e_simd1(double *f, double *g, double *ri,
                      int li, int lj, int lk, int ll, CINTEnvVars *envs);
void CINTx1j_2e_simd1(double *f, double *g, double *rj,
                      int li, int lj, int lk, int ll, CINTEnvVars *envs);
void CINTx1k_2e_simd1(double *f, double *g, double *rk,
                      int li, int lj, int lk, int ll, CINTEnvVars *envs);
void CINTx1l_2e_simd1(double *f, double *g, double *rl,
                      int li, int lj, int lk, int ll, CINTEnvVars *envs);

#ifdef WITH_F12
void CINTinit_int2e_yp_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                               int *atm, int natm, int *bas, int nbas, double *env);
void CINTinit_int2e_stg_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                                int *atm, int natm, int *bas, int nbas, double *env);
#endif
 
#ifdef WITH_GTG
void CINTinit_int2e_gtg_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                                int *atm, int natm, int *bas, int nbas, double *env);
void CINTinit_int3c2e_gtg_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                                  int *atm, int natm, int *bas, int nbas, double *env);
void CINTinit_int2c2e_gtg_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                                  int *atm, int natm, int *bas, int nbas, double *env);
#endif


#define G2E_D_I(f, g, li, lj, lk, ll)   CINTnabla1i_2e(f, g, li, lj, lk, ll, envs)
#define G2E_D_J(f, g, li, lj, lk, ll)   CINTnabla1j_2e(f, g, li, lj, lk, ll, envs)
#define G2E_D_K(f, g, li, lj, lk, ll)   CINTnabla1k_2e(f, g, li, lj, lk, ll, envs)
#define G2E_D_L(f, g, li, lj, lk, ll)   CINTnabla1l_2e(f, g, li, lj, lk, ll, envs)
/* r-R_0, R_0 is (0,0,0) */
#define G2E_R0I(f, g, li, lj, lk, ll)   CINTx1i_2e(f, g, envs->ri, li, lj, lk, ll, envs)
#define G2E_R0J(f, g, li, lj, lk, ll)   CINTx1j_2e(f, g, envs->rj, li, lj, lk, ll, envs)
#define G2E_R0K(f, g, li, lj, lk, ll)   CINTx1k_2e(f, g, envs->rk, li, lj, lk, ll, envs)
#define G2E_R0L(f, g, li, lj, lk, ll)   CINTx1l_2e(f, g, envs->rl, li, lj, lk, ll, envs)
/* r-R_C, R_C is common origin */
#define G2E_RCI(f, g, li, lj, lk, ll)   CINTx1i_2e(f, g, dri, li, lj, lk, ll, envs)
#define G2E_RCJ(f, g, li, lj, lk, ll)   CINTx1j_2e(f, g, drj, li, lj, lk, ll, envs)
#define G2E_RCK(f, g, li, lj, lk, ll)   CINTx1k_2e(f, g, drk, li, lj, lk, ll, envs)
#define G2E_RCL(f, g, li, lj, lk, ll)   CINTx1l_2e(f, g, drl, li, lj, lk, ll, envs)
/* origin from center of each basis
 * x1[ijkl]_2e(f, g, ng, li, lj, lk, ll, 0d0) */
#define G2E_R_I(f, g, li, lj, lk, ll)   f = g + envs->g_stride_i * SIMDD
#define G2E_R_K(f, g, li, lj, lk, ll)   f = g + envs->g_stride_k * SIMDD
#define G2E_R_L(f, g, li, lj, lk, ll)   f = g + envs->g_stride_l * SIMDD
#define G2E_R_J(f, g, li, lj, lk, ll)   f = g + envs->g_stride_j * SIMDD


#define G2E_D_I_SIMD1(f, g, li, lj, lk, ll)   CINTnabla1i_2e_simd1(f, g, li, lj, lk, ll, envs)
#define G2E_D_J_SIMD1(f, g, li, lj, lk, ll)   CINTnabla1j_2e_simd1(f, g, li, lj, lk, ll, envs)
#define G2E_D_K_SIMD1(f, g, li, lj, lk, ll)   CINTnabla1k_2e_simd1(f, g, li, lj, lk, ll, envs)
#define G2E_D_L_SIMD1(f, g, li, lj, lk, ll)   CINTnabla1l_2e_simd1(f, g, li, lj, lk, ll, envs)
/* r-R_0, R_0 is (0,0,0) */
#define G2E_R0I_SIMD1(f, g, li, lj, lk, ll)   CINTx1i_2e_simd1(f, g, envs->ri, li, lj, lk, ll, envs)
#define G2E_R0J_SIMD1(f, g, li, lj, lk, ll)   CINTx1j_2e_simd1(f, g, envs->rj, li, lj, lk, ll, envs)
#define G2E_R0K_SIMD1(f, g, li, lj, lk, ll)   CINTx1k_2e_simd1(f, g, envs->rk, li, lj, lk, ll, envs)
#define G2E_R0L_SIMD1(f, g, li, lj, lk, ll)   CINTx1l_2e_simd1(f, g, envs->rl, li, lj, lk, ll, envs)
/* r-R_C, R_C is common origin */
#define G2E_RCI_SIMD1(f, g, li, lj, lk, ll)   CINTx1i_2e_simd1(f, g, dri, li, lj, lk, ll, envs)
#define G2E_RCJ_SIMD1(f, g, li, lj, lk, ll)   CINTx1j_2e_simd1(f, g, drj, li, lj, lk, ll, envs)
#define G2E_RCK_SIMD1(f, g, li, lj, lk, ll)   CINTx1k_2e_simd1(f, g, drk, li, lj, lk, ll, envs)
#define G2E_RCL_SIMD1(f, g, li, lj, lk, ll)   CINTx1l_2e_simd1(f, g, drl, li, lj, lk, ll, envs)
/* origin from center of each basis
 * x1[ijkl]_2e(f, g, ng, li, lj, lk, ll, 0d0) */
#define G2E_R_I_SIMD1(f, g, li, lj, lk, ll)   f = g + envs->g_stride_i
#define G2E_R_K_SIMD1(f, g, li, lj, lk, ll)   f = g + envs->g_stride_k
#define G2E_R_L_SIMD1(f, g, li, lj, lk, ll)   f = g + envs->g_stride_l
#define G2E_R_J_SIMD1(f, g, li, lj, lk, ll)   f = g + envs->g_stride_j


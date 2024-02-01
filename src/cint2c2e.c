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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "cint_bas.h"
#include "misc.h"
#include "g2e.h"
#include "optimizer.h"
#include "cint1e.h"
#include "cint2e.h"
#include "cart2sph.h"
#include "c2f.h"

#define SHLTYPi       0
#define SHLTYPk       1

#define ALIAS_ADDR_IF_EQUAL(x, y) \
        if (y##_ctr == 1) { \
                gctr[SHLTYP##x] = gctr[SHLTYP##y]; \
                x##empty = y##empty; \
        } else { \
                gctr[SHLTYP##x] = g1; \
                g1 += len##x; \
        }

#define PRIM2CTR(psymb, csymb) \
        if (csymb##_ctr > 1) {\
                if (*csymb##empty) { \
                        fp2c[np2c] = CINTprim_to_ctr_0; \
                } else { \
                        fp2c[np2c] = CINTprim_to_ctr_1; \
                } \
                shltyp[np2c] = SHLTYP##csymb; \
                gprim[np2c] = gctr[SHLTYP##psymb]; \
                iprim[np2c] = csymb##p; \
                np2c++; \
        } \
        *csymb##empty = 0; \

#define POP_PRIM2CTR \
        for (i = 0; i < np2c; i++) { \
                it = shltyp[i]; \
                im = iprim[i]; \
                (*(fp2c[i]))(gctr[it], gprim[i], coeff[it]+im, \
                             ngp[it], x_prim[it], x_ctr[it], \
                             non0ctr[it][im], non0idx[it]+im*x_ctr[it]); \
                empty_overall = 0; \
        } \
        cum = 0; \
        np2c = 0;

#define POP_PRIM2CTR_AND_SET0 \
        for (i = 0; i < np2c; i++) { \
                it = shltyp[i]; \
                if (it != SHLTYPi) { \
                        im = iprim[i]; \
                        (*(fp2c[i]))(gctr[it], gprim[i], coeff[it]+im, \
                                     ngp[it], x_prim[it], x_ctr[it], \
                                     non0ctr[it][im], non0idx[it]+im*x_ctr[it]); \
                        empty_overall = 0; \
                } else if (fp2c[i] == CINTiprim_to_ctr_0) { \
                        double *pout = gctr[SHLTYPi]; \
                        int k; \
                        for (k = 0; k < ngp[SHLTYPk]; k++) { \
                                pout[k] = 0.; \
                        } \
                } \
        } \
        cum = 0; \
        np2c = 0;

#define PUSH \
        if (cum == SIMDD) { \
                if ((*envs->f_g0_2e)(g, cutoff, &bc, envs, cum)) { \
                        (*envs->f_gout)(gout, g, idx, envs); \
                        POP_PRIM2CTR; \
                } else { \
                        POP_PRIM2CTR_AND_SET0; \
                } \
        } \
        envs->ai[cum] = ai[ip]; \
        envs->ak[cum] = ak[kp]; \
        envs->fac[cum] = fac1k; \
        cutoff[cum] = expcutoff; \
        if (*iempty) { \
                fp2c[np2c] = CINTiprim_to_ctr_0; \
                *iempty = 0; \
        } else { \
                fp2c[np2c] = CINTiprim_to_ctr_1; \
        } \
        gprim[np2c] = gout + cum * ngp[0]; \
        iprim[np2c] = ip; \
        shltyp[np2c] = 0; \
        cum++; \
        np2c++;

#define INITSIMD \
        double *gx = g; \
        double *gy = g + envs->g_size * SIMDD; \
        __MD r1 = MM_SET1(1.); \
        for (i = 0; i < envs->nrys_roots; i++) { \
                MM_STORE(gx+i*SIMDD, r1); \
                MM_STORE(gy+i*SIMDD, r1); \
        } \
        int cum = 0; \
        int np2c = 0; \
        double *gprim[SIMDD*2]; \
        int shltyp[SIMDD*2]; \
        int iprim[SIMDD*2]; \
        void (*fp2c[SIMDD*2])(); \
        MM_STORE(envs->ai, MM_SET1(1.)); \
        MM_STORE(envs->ak, MM_SET1(1.)); \
        MM_STORE(envs->fac, MM_SET1(0.));

#define RUN_REST \
        if (cum == 1) { \
                if ((*envs->f_g0_2e_simd1)(g, cutoff, &bc, envs, 0)) { \
                        (*envs->f_gout_simd1)(gout, g, idx, envs); \
                        POP_PRIM2CTR; \
                } else { \
                        POP_PRIM2CTR_AND_SET0; \
                } \
        } else if (cum > 1) { \
                if ((*envs->f_g0_2e)(g, cutoff, &bc, envs, cum)) { \
                        (*envs->f_gout)(gout, g, idx, envs); \
                        POP_PRIM2CTR; \
                } else { \
                        POP_PRIM2CTR_AND_SET0; \
                } \
        } else { \
                assert(np2c == 0); \
        }

#define TRANSPOSE(a) \
        if (*empty) { \
                CINTdmat_transpose(out, a, nf*nc, n_comp); \
                *empty = 0; \
        } else { \
                CINTdplus_transpose(out, a, nf*nc, n_comp); \
        }

int CINT2c2e_loop_nopt(double *out, CINTEnvVars *envs, double *cache, int *empty)
{
        int *shls = envs->shls;
        int *bas = envs->bas;
        double *env = envs->env;
        int i_sh = shls[0];
        int k_sh = shls[1];
        int i_ctr = envs->x_ctr[0];
        int k_ctr = envs->x_ctr[1];
        int i_prim = bas(NPRIM_OF, i_sh);
        int k_prim = bas(NPRIM_OF, k_sh);
        int *x_ctr = envs->x_ctr;
        int x_prim[2] = {i_prim, k_prim};
        double *ai = env + bas(PTR_EXP, i_sh);
        double *ak = env + bas(PTR_EXP, k_sh);
        double *ci = env + bas(PTR_COEFF, i_sh);
        double *ck = env + bas(PTR_COEFF, k_sh);
        double *coeff[2] = {ci, ck};
        double expcutoff = envs->expcutoff;
        int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        int nf = envs->nf;
        double fac1k;
        int ip, kp, i, it, im;
        int _empty[2] = {1, 1};
        int *iempty = _empty + 0;
        int *kempty = _empty + 1;
        int empty_overall = 1;
        int ngp[2];
        ngp[0] = nf * n_comp;
        ngp[1] = ngp[0] * i_ctr;
        int leng = envs->g_size * 3 * ((1<<envs->gbits)+1) * SIMDD;
        int len0 = ngp[0] * SIMDD;
        int leni = ALIGN_UP(ngp[1], SIMDD);
        int lenk = ALIGN_UP(ngp[1] * k_ctr, SIMDD);

        double *gout, *g, *g1;
        double *gctr[2];
        if (n_comp == 1) {
                MALLOC_INSTACK(g1, lenk);
                gctr[SHLTYPk] = out;
                *kempty = *empty;
        } else {
                MALLOC_INSTACK(gctr[SHLTYPk], lenk);
                g1 = out;
        }
        ALIAS_ADDR_IF_EQUAL(i, k);
        MALLOC_INSTACK(gout, leng+len0);
        g = gout + len0;  // for gx, gy, gz

        ALIGNMM Rys2eT bc;
        ALIGNMM double cutoff[SIMDD];
        int *idx;
        MALLOC_INSTACK(idx, envs->nf * 3);
        int *non0ctri, *non0ctrk;
        int *non0idxi, *non0idxk;
        MALLOC_INSTACK(non0ctri, i_prim+k_prim+i_prim*i_ctr+k_prim*k_ctr);
        non0ctrk = non0ctri + i_prim;
        non0idxi = non0ctrk + k_prim;
        non0idxk = non0idxi + i_prim*i_ctr;
        CINTg2c_index_xyz(idx, envs);
        CINTOpt_non0coeff_byshell(non0idxi, non0ctri, coeff[0], i_prim, i_ctr);
        CINTOpt_non0coeff_byshell(non0idxk, non0ctrk, coeff[1], k_prim, k_ctr);
        int *non0ctr[2] = {non0ctri, non0ctrk};
        int *non0idx[2] = {non0idxi, non0idxk};
        double common_factor = envs->common_factor;
        INITSIMD;

        for (kp = 0; kp < k_prim; kp++) {
                if (k_ctr == 1) {
                        fac1k = common_factor * ck[kp];
                } else {
                        fac1k = common_factor;
                        *iempty = 1;
                }
                for (ip = 0; ip < i_prim; ip++) {
                        PUSH;
                } // end loop i_prim
                if (!*iempty) {
                        PRIM2CTR(i, k);
                }
        } // end loop k_prim
        RUN_REST;

        if (n_comp > 1 && !empty_overall) {
                int nc = i_ctr * k_ctr;
                TRANSPOSE(gctr[SHLTYPk]);
        }
        *empty &= empty_overall;
        return !empty_overall;
}


int CINT2c2e_loop(double *out, CINTEnvVars *envs, double *cache, int *empty)
{
        int *shls  = envs->shls;
        int *bas = envs->bas;
        double *env = envs->env;
        int i_sh = shls[0];
        int k_sh = shls[1];
        int i_ctr  = envs->x_ctr[0];
        int k_ctr  = envs->x_ctr[1];
        int i_prim = bas(NPRIM_OF, i_sh);
        int k_prim = bas(NPRIM_OF, k_sh);
        int *x_ctr = envs->x_ctr;
        int x_prim[2] = {i_prim, k_prim};
        double *ai = env + bas(PTR_EXP, i_sh);
        double *ak = env + bas(PTR_EXP, k_sh);
        double *ci = env + bas(PTR_COEFF, i_sh);
        double *ck = env + bas(PTR_COEFF, k_sh);
        double *coeff[2] = {ci, ck};
        double expcutoff = envs->expcutoff;
        int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        int nf = envs->nf;
        double fac1k;
        int ip, kp, i, it, im;
        int _empty[2] = {1, 1};
        int *iempty = _empty + 0;
        int *kempty = _empty + 1;
        int empty_overall = 1;
        int ngp[2];
        ngp[0] = nf * n_comp;
        ngp[1] = ngp[0] * i_ctr;
        // (irys,i,j,k,l,coord,0:1); +1 for nabla-r12
        int leng = envs->g_size * 3 * ((1<<envs->gbits)+1) * SIMDD;
        int len0 = ngp[0] * SIMDD;
        int leni = ALIGN_UP(ngp[1], SIMDD);
        int lenk = ALIGN_UP(ngp[1] * k_ctr, SIMDD);
        double *gout, *g, *g1;
        double *gctr[2];
        if (n_comp == 1) {
                MALLOC_INSTACK(g1, lenk);
                gctr[SHLTYPk] = out;
                *kempty = *empty;
        } else {
                MALLOC_INSTACK(gctr[SHLTYPk], lenk);
                g1 = out;
        }
        ALIAS_ADDR_IF_EQUAL(i, k);
        MALLOC_INSTACK(gout, leng+len0);
        g = gout + len0;  // for gx, gy, gz

        ALIGNMM Rys2eT bc;
        ALIGNMM double cutoff[SIMDD];
        double common_factor = envs->common_factor;

        CINTOpt *opt = envs->opt;
        int *idx = opt->index_xyz_array[envs->i_l*LMAX1+envs->k_l];
        int *non0ctr[2] = {opt->non0ctr[i_sh], opt->non0ctr[k_sh]};
        int *non0idx[2] = {opt->sortedidx[i_sh], opt->sortedidx[k_sh]};

        INITSIMD;

        for (kp = 0; kp < k_prim; kp++) {
                if (k_ctr == 1) {
                        fac1k = common_factor * ck[kp];
                } else {
                        fac1k = common_factor;
                        *iempty = 1;
                }
                for (ip = 0; ip < i_prim; ip++) {
                        PUSH;
                } // end loop i_prim
                if (!*iempty) {
                        PRIM2CTR(i, k);
                }
        } // end loop k_prim
        RUN_REST;

        if (n_comp > 1 && !empty_overall) {
                int nc = i_ctr * k_ctr;
                TRANSPOSE(gctr[SHLTYPk]);
        }
        *empty &= empty_overall;
        return !empty_overall;
}

void CINTinit_int2c2e_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                              int *atm, int natm, int *bas, int nbas, double *env)
{
        envs->natm = natm;
        envs->nbas = nbas;
        envs->atm = atm;
        envs->bas = bas;
        envs->env = env;
        envs->shls = shls;

        int i_sh = shls[0];
        int k_sh = shls[1];
        envs->i_l = bas(ANG_OF, i_sh);
        envs->j_l = 0;
        envs->k_l = bas(ANG_OF, k_sh);
        envs->l_l = 0;
        envs->x_ctr[0] = bas(NCTR_OF, i_sh);
        envs->x_ctr[1] = bas(NCTR_OF, k_sh);
        envs->x_ctr[2] = 1;
        envs->x_ctr[3] = 1;
        envs->nfi = (envs->i_l+1)*(envs->i_l+2)/2;
        envs->nfj = 1;
        envs->nfk = (envs->k_l+1)*(envs->k_l+2)/2;
        envs->nfl = 1;
        envs->nf = envs->nfi * envs->nfk;

        envs->ri = env + atm(PTR_COORD, bas(ATOM_OF, i_sh));
        envs->rk = env + atm(PTR_COORD, bas(ATOM_OF, k_sh));

        envs->common_factor = (M_PI*M_PI*M_PI)*2/SQRTPI
                * CINTcommon_fac_sp(envs->i_l) * CINTcommon_fac_sp(envs->k_l);
        if (env[PTR_EXPCUTOFF] == 0) {
                envs->expcutoff = EXPCUTOFF;
        } else {
                envs->expcutoff = MAX(MIN_EXPCUTOFF, env[PTR_EXPCUTOFF]);
        }

        envs->gbits = ng[GSHIFT];
        envs->ncomp_e1 = ng[POS_E1];
        envs->ncomp_e2 = ng[POS_E2];
        envs->ncomp_tensor = ng[TENSOR];

        envs->li_ceil = envs->i_l + ng[IINC];
        envs->lj_ceil = 0;
        envs->lk_ceil = envs->k_l + ng[KINC];
        envs->ll_ceil = 0;
        int rys_order =(envs->li_ceil + envs->lk_ceil)/2 + 1;
        int nrys_roots = rys_order;
        double omega = env[PTR_RANGE_OMEGA];
        if (omega < 0 && rys_order <= 3) {
                nrys_roots *= 2;
        }
        envs->rys_order = rys_order;
        envs->nrys_roots = nrys_roots;

        int dli = envs->li_ceil + 1;
        int dlk = envs->lk_ceil + 1;
        envs->g_stride_i = nrys_roots;
        envs->g_stride_k = nrys_roots * dli;
        envs->g_stride_l = envs->g_stride_k;
        envs->g_size     = nrys_roots * dli * dlk;

        MM_STORE(envs->aj, MM_SET1(0.));
        MM_STORE(envs->al, MM_SET1(0.));
        MM_STORE(envs->rij+0*SIMDD, MM_SET1(envs->ri[0]));
        MM_STORE(envs->rij+1*SIMDD, MM_SET1(envs->ri[1]));
        MM_STORE(envs->rij+2*SIMDD, MM_SET1(envs->ri[2]));
        MM_STORE(envs->rkl+0*SIMDD, MM_SET1(envs->rk[0]));
        MM_STORE(envs->rkl+1*SIMDD, MM_SET1(envs->rk[1]));
        MM_STORE(envs->rkl+2*SIMDD, MM_SET1(envs->rk[2]));
        envs->g2d_ijmax = envs->g_stride_i;
        envs->g2d_klmax = envs->g_stride_k;
        envs->rkrl[0] = envs->rk[0];
        envs->rkrl[1] = envs->rk[1];
        envs->rkrl[2] = envs->rk[2];
        envs->rirj[0] = envs->ri[0];
        envs->rirj[1] = envs->ri[1];
        envs->rirj[2] = envs->ri[2];
        envs->rx_in_rklrx = envs->rk;
        envs->rx_in_rijrx = envs->ri;

        if (rys_order <= 2) {
                envs->f_g0_2d4d = &CINTg0_2e_2d4d_unrolled;
                envs->f_g0_2d4d_simd1 = &CINTg0_2e_2d4d_unrolled_simd1;
                if (rys_order != nrys_roots) {
                        envs->f_g0_2d4d = &CINTsrg0_2e_2d4d_unrolled;
                        envs->f_g0_2d4d_simd1 = &CINTsrg0_2e_2d4d_unrolled_simd1;
                }
        } else {
                envs->f_g0_2d4d = &CINTg0_2e_2d;
                envs->f_g0_2d4d_simd1 = &CINTg0_2e_2d_simd1;
        }
        envs->f_g0_2e = &CINTg0_2e;
        envs->f_g0_2e_simd1 = &CINTg0_2e_simd1;

        // for CINTg2c_index_xyz and c2s_sph_1e function
        envs->j_l = envs->k_l;
        envs->nfj = envs->nfk;
        envs->g_stride_j = envs->g_stride_k;
}


CACHE_SIZE_T CINT2c2e_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
                          double *cache, void (*f_c2s)())
{
        if (out == NULL) {
                return int1e_cache_size(envs);
        }
        int *x_ctr = envs->x_ctr;
        int nc = envs->nf * x_ctr[0] * x_ctr[1];
        int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        double *stack = NULL;
        if (cache == NULL) {
                int cache_size = int1e_cache_size(envs);
                stack = _mm_malloc(sizeof(double)*cache_size, sizeof(double)*SIMDD);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, nc*n_comp);

        int n;
        int empty = 1;
        if (opt != NULL) {
                envs->opt = opt;
                CINT2c2e_loop(gctr, envs, cache, &empty);
        } else {
                CINT2c2e_loop_nopt(gctr, envs, cache, &empty);
        }

        int counts[4];
        if (f_c2s == c2s_sph_1e) {
                counts[0] = (envs->i_l*2+1) * x_ctr[0];
                counts[1] = (envs->k_l*2+1) * x_ctr[1];
        } else {
                counts[0] = envs->nfi * x_ctr[0];
                counts[1] = envs->nfk * x_ctr[1];
        }
        counts[2] = 1;
        counts[3] = 1;
        if (dims == NULL) {
                dims = counts;
        }
        int nout = dims[0] * dims[1];
        if (!empty) {
                for (n = 0; n < n_comp; n++) {
                        (*f_c2s)(out+nout*n, gctr+nc*n, dims, envs, cache);
                }
        } else {
                for (n = 0; n < n_comp; n++) {
                        c2s_dset0(out+nout*n, dims, counts);
                }
        }
        if (stack != NULL) {
                free(stack);
        }
        return !empty;
}
// (spinor|spinor)
CACHE_SIZE_T CINT2c2e_spinor_drv(double complex *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
                        double *cache, void (*f_e1_c2s)())
{
        if (envs->ncomp_e1 > 1 || envs->ncomp_e2 > 1) {
                fprintf(stderr, "CINT2c2e_spinor_drv not implemented\n");
                exit(1);
        }
        if (out == NULL) {
                return int1e_cache_size(envs);
        }
        int *x_ctr = envs->x_ctr;
        int counts[4];
        counts[0] = CINTcgto_spinor(envs->shls[0], envs->bas);
        counts[1] = CINTcgto_spinor(envs->shls[1], envs->bas);
        counts[2] = 1;
        counts[3] = 1;
        int nc = envs->nf * x_ctr[0] * x_ctr[1];
        int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        double *stack = NULL;
        if (cache == NULL) {
                int cache_size = int1e_cache_size(envs);
                stack = _mm_malloc(sizeof(double)*cache_size, sizeof(double)*SIMDD);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, nc*n_comp);

        int n;
        int empty = 1;
        if (opt != NULL) {
                envs->opt = opt;
                CINT2c2e_loop(gctr, envs, cache, &empty);
        } else {
                CINT2c2e_loop_nopt(gctr, envs, cache, &empty);
        }

        if (dims == NULL) {
                dims = counts;
        }
        int nout = dims[0] * dims[1];
        if (!empty) {
                for (n = 0; n < n_comp; n++) {
                        (*f_e1_c2s)(out+nout*n, gctr, dims, envs, cache);
                        gctr += nc;
                }
        } else {
                for (n = 0; n < n_comp; n++) {
                        c2s_zset0(out+nout*n, dims, counts);
                }
        }
        if (stack != NULL) {
                free(stack);
        }
        return !empty;
}

CACHE_SIZE_T int2c2e_sph(double *out, int *dims, int *shls, int *atm, int natm,
                int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        envs.f_gout_simd1 = &CINTgout2e_simd1;
        return CINT2c2e_drv(out, dims, &envs, opt, cache, &c2s_sph_1e);
}
void int2c2e_optimizer(CINTOpt **opt, int *atm, int natm,
                       int *bas, int nbas, double *env)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTall_2c2e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}

CACHE_SIZE_T int2c2e_cart(double *out, int *dims, int *shls, int *atm, int natm,
                 int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        envs.f_gout_simd1 = &CINTgout2e_simd1;
        return CINT2c2e_drv(out, dims, &envs, opt, cache, &c2s_cart_1e);
}

CACHE_SIZE_T int2c2e_spinor(double complex *out, int *dims, int *shls, int *atm, int natm,
                   int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        envs.f_gout_simd1 = &CINTgout2e_simd1;
        return CINT2c2e_spinor_drv(out, dims, &envs, opt, cache, &c2s_sf_1e);
}


ALL_CINT(int2c2e)
//ALL_CINT_FORTRAN_(cint2c2e)


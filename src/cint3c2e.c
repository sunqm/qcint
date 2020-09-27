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
#include <string.h>
#include <math.h>
#include "cint_bas.h"
#include "misc.h"
#include "g2e.h"
#include "optimizer.h"
#include "cint2e.h"
#include "cart2sph.h"
#include "c2f.h"

#define SHLTYPi       0
#define SHLTYPj       1
#define SHLTYPk       2

#define ALIAS_ADDR_IF_EQUAL(x, y) \
        if (y##_ctr == 1) { \
                gctr[SHLTYP##x] = gctr[SHLTYP##y]; \
                x##empty = y##empty; \
        } else { \
                gctr[SHLTYP##x] = g1; \
                g1 += len##x; \
        }

#define PRIM2CTR(ctrsymb, gp) \
        if (ctrsymb##_ctr > 1) {\
                if (*ctrsymb##empty) { \
                        fp2c[np2c] = CINTprim_to_ctr_0; \
                } else { \
                        fp2c[np2c] = CINTprim_to_ctr_1; \
                } \
                shltyp[np2c] = SHLTYP##ctrsymb; \
                gprim[np2c] = gp; \
                iprim[np2c] = ctrsymb##p; \
                np2c++; \
        } \
        *ctrsymb##empty = 0; \

#define POP_PRIM2CTR \
        for (i = 0; i < np2c; i++) { \
                it = shltyp[i]; \
                im = iprim[i]; \
                (*(fp2c[i]))(gctr[it], gprim[i], coeff[it]+im, \
                             ngp[it], x_prim[it], x_ctr[it], \
                             non0ctr[it][im], non0idx[it]+im*x_ctr[it]); \
        } \
        cum = 0; \
        np2c = 0;

#define PUSH(RIJ, EXPIJ) \
        if (cum == SIMDD) { \
                if ((*envs->f_g0_2e)(g, &bc, envs, cum)) { \
                        (*envs->f_gout)(gout, g, idx, envs); \
                        POP_PRIM2CTR; \
                } else { \
                        cum = 0; \
                        np2c = 0; \
                } \
        } \
        envs->ai[cum] = ai[ip]; \
        envs->aj[cum] = aj[jp]; \
        envs->ak[cum] = ak[kp]; \
        envs->rij[0*SIMDD+cum] = *(RIJ+0); \
        envs->rij[1*SIMDD+cum] = *(RIJ+1); \
        envs->rij[2*SIMDD+cum] = *(RIJ+2); \
        fac1i = fac1j * EXPIJ; \
        envs->fac[cum] = fac1i; \
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
        double *gprim[SIMDD*3]; \
        int shltyp[SIMDD*3]; \
        int iprim[SIMDD*3]; \
        void (*fp2c[SIMDD*3])(); \
        MM_STORE(envs->ai, MM_SET1(1.)); \
        MM_STORE(envs->aj, MM_SET1(1.)); \
        MM_STORE(envs->ak, MM_SET1(1.)); \
        MM_STORE(envs->fac, MM_SET1(0.));

#define RUN_REST \
        if (cum == 1) { \
                if ((*envs->f_g0_2e_simd1)(g, &bc, envs, 0)) { \
                        (*envs->f_gout_simd1)(gout, g, idx, envs); \
                        POP_PRIM2CTR; \
                } \
        } else if (cum > 1) { \
                r1 = MM_SET1(1.); \
                for (i = 0; i < envs->nrys_roots; i++) { \
                        MM_STORE(bc.u+i*SIMDD, r1); \
                        MM_STORE(bc.w+i*SIMDD, r1); \
                } \
                if ((*envs->f_g0_2e)(g, &bc, envs, cum)) { \
                        (*envs->f_gout)(gout, g, idx, envs); \
                        POP_PRIM2CTR; \
                } \
        }

int CINT3c2e_loop_nopt(double *out, CINTEnvVars *envs, double *cache)
{
        int *shls = envs->shls;
        int *bas = envs->bas;
        double *env = envs->env;
        int i_sh = shls[0];
        int j_sh = shls[1];
        int k_sh = shls[2];
        int i_ctr = envs->x_ctr[0];
        int j_ctr = envs->x_ctr[1];
        int k_ctr = envs->x_ctr[2];
        int i_prim = bas(NPRIM_OF, i_sh);
        int j_prim = bas(NPRIM_OF, j_sh);
        int k_prim = bas(NPRIM_OF, k_sh);
        int *x_ctr = envs->x_ctr;
        int x_prim[3] = {i_prim, j_prim, k_prim};
        double *ai = env + bas(PTR_EXP, i_sh);
        double *aj = env + bas(PTR_EXP, j_sh);
        double *ak = env + bas(PTR_EXP, k_sh);
        double *ci = env + bas(PTR_COEFF, i_sh);
        double *cj = env + bas(PTR_COEFF, j_sh);
        double *ck = env + bas(PTR_COEFF, k_sh);
        double *coeff[3] = {ci, cj, ck};
        int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        int nf = envs->nf;
        double fac1i, fac1j, fac1k;
        int ip, jp, kp, i, it, im;
        int empty[3] = {1, 1, 1};
        int *iempty = empty + 0;
        int *jempty = empty + 1;
        int *kempty = empty + 2;
        int ngp[4];
        ngp[0] = nf * n_comp;
        ngp[1] = ngp[0] * i_ctr;
        ngp[2] = ngp[1] * j_ctr;
        ngp[3] = ngp[2] * k_ctr;
        int leng = envs->g_size * 3 * ((1<<envs->gbits)+1) * SIMDD;
        int len0 = ngp[0] * SIMDD;
        int leni = ALIGN_UP(ngp[1], SIMDD);
        int lenj = ALIGN_UP(ngp[2], SIMDD);
        int lenk = ALIGN_UP(ngp[3], SIMDD);

        double *gout, *g, *g1;
        double *gctr[3];
        if (n_comp == 1) {
                // patch SIMDD for leni, lenj with s functions
                MALLOC_INSTACK(g1, lenk+SIMDD);
                gctr[SHLTYPk] = out;
        } else {
                // enlarge out by SIMDD, for leni, lenj with s functions
                cache += SIMDD;
                g1 = out;  // Use out as cache for gctrk, gctrj, gctri
                MALLOC_INSTACK(gctr[SHLTYPk], lenk);
        }
        ALIAS_ADDR_IF_EQUAL(j, k);
        ALIAS_ADDR_IF_EQUAL(i, j);
        MALLOC_INSTACK(gout, leng+len0);
        g = gout + len0;  // for gx, gy, gz

        ALIGNMM Rys2eT bc;
        PairData pdata[i_prim*j_prim];
        int idx[nf*3];
        int non0ctri[i_prim];
        int non0ctrj[j_prim];
        int non0ctrk[k_prim];
        int non0idxi[i_prim*i_ctr];
        int non0idxj[j_prim*j_ctr];
        int non0idxk[k_prim*k_ctr];
        double expcutoff = envs->expcutoff;

        *kempty = CINTset_pairdata(pdata, ai, aj, envs->ri, envs->rj,
                                   envs->li_ceil, envs->lj_ceil,
                                   i_prim, j_prim, SQUARE(envs->rirj), expcutoff);
        if (*kempty) {
                goto normal_end;
        }

        CINTg4c_index_xyz(idx, envs);
        CINTOpt_non0coeff_byshell(non0idxi, non0ctri, coeff[0], i_prim, i_ctr);
        CINTOpt_non0coeff_byshell(non0idxj, non0ctrj, coeff[1], j_prim, j_ctr);
        CINTOpt_non0coeff_byshell(non0idxk, non0ctrk, coeff[2], k_prim, k_ctr);
        int *non0ctr[3] = {non0ctri, non0ctrj, non0ctrk};
        int *non0idx[3] = {non0idxi, non0idxj, non0idxk};
        double common_factor = envs->common_factor * (M_PI*M_PI*M_PI)*2/SQRTPI
                * CINTcommon_fac_sp(envs->i_l) * CINTcommon_fac_sp(envs->j_l)
                * CINTcommon_fac_sp(envs->k_l);
        INITSIMD;

        PairData *pdata_ij;
        *kempty = 1;
        for (kp = 0; kp < k_prim; kp++) {
                if (k_ctr == 1) {
                        fac1k = common_factor * ck[kp];
                } else {
                        fac1k = common_factor;
                        *jempty = 1;
                }

                pdata_ij = pdata;
                for (jp = 0; jp < j_prim; jp++) {
                        if (j_ctr == 1) {
                                fac1j = fac1k * cj[jp];
                        } else {
                                fac1j = fac1k;
                                *iempty = 1;
                        }
                        for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                if (pdata_ij->cceij < expcutoff) {
                                        PUSH(pdata_ij->rij, pdata_ij->eij);
                                }
                        } // end loop i_prim
                        if (!*iempty) {
                                PRIM2CTR(j, gctr[SHLTYPi]);
                        }
                } // end loop j_prim
                if (!*jempty) {
                        PRIM2CTR(k, gctr[SHLTYPj]);
                }
        } // end loop k_prim
        RUN_REST;

        if (n_comp > 1 && !*kempty) {
                int nc = i_ctr * j_ctr * k_ctr;
                CINTdmat_transpose(out, gctr[SHLTYPk], nf*nc, n_comp);
        }
normal_end:
        return !*kempty;
}


int CINT3c2e_loop(double *out, CINTEnvVars *envs, CINTOpt *opt, double *cache)
{
        int *shls  = envs->shls;
        int i_sh = shls[0];
        int j_sh = shls[1];
        int k_sh = shls[2];
        if (opt->data_ptr != NULL &&
            opt->data_ptr[i_sh*opt->nbas+j_sh] == NOVALUE) {
                return 0;
        }
        int *bas = envs->bas;
        double *env = envs->env;
        int i_ctr  = envs->x_ctr[0];
        int j_ctr  = envs->x_ctr[1];
        int k_ctr  = envs->x_ctr[2];
        int i_prim = bas(NPRIM_OF, i_sh);
        int j_prim = bas(NPRIM_OF, j_sh);
        int k_prim = bas(NPRIM_OF, k_sh);
        int *x_ctr = envs->x_ctr;
        int x_prim[3] = {i_prim, j_prim, k_prim};
        double *ai = env + bas(PTR_EXP, i_sh);
        double *aj = env + bas(PTR_EXP, j_sh);
        double *ak = env + bas(PTR_EXP, k_sh);
        double *ci = env + bas(PTR_COEFF, i_sh);
        double *cj = env + bas(PTR_COEFF, j_sh);
        double *ck = env + bas(PTR_COEFF, k_sh);
        double *coeff[3] = {ci, cj, ck};
        int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        int nf = envs->nf;
        double fac1i, fac1j, fac1k;
        int ip, jp, kp, i, it, im;
        int empty[3] = {1, 1, 1};
        int *iempty = empty + 0;
        int *jempty = empty + 1;
        int *kempty = empty + 2;
        int ngp[4];
        ngp[0] = nf * n_comp;
        ngp[1] = ngp[0] * i_ctr;
        ngp[2] = ngp[1] * j_ctr;
        ngp[3] = ngp[2] * k_ctr;
        // (irys,i,j,k,l,coord,0:1); +1 for nabla-r12
        int leng = envs->g_size * 3 * ((1<<envs->gbits)+1) * SIMDD;
        int len0 = ngp[0] * SIMDD;
        int leni = ALIGN_UP(ngp[1], SIMDD);
        int lenj = ALIGN_UP(ngp[2], SIMDD);
        int lenk = ALIGN_UP(ngp[3], SIMDD);

        double *gout, *g, *g1;
        double *gctr[3];
        if (n_comp == 1) {
                // patch SIMDD for leni, lenj with s functions
                MALLOC_INSTACK(g1, lenk+SIMDD);
                gctr[SHLTYPk] = out;
        } else {
                // enlarge out by SIMDD, for leni, lenj with s functions
                cache += SIMDD;
                g1 = out;  // Use out as cache for gctrk, gctrj, gctri
                MALLOC_INSTACK(gctr[SHLTYPk], lenk);
        }
        ALIAS_ADDR_IF_EQUAL(j, k);
        ALIAS_ADDR_IF_EQUAL(i, j);
        MALLOC_INSTACK(gout, leng+len0);
        g = gout + len0;  // for gx, gy, gz

        ALIGNMM Rys2eT bc;
        double common_factor = envs->common_factor * (M_PI*M_PI*M_PI)*2/SQRTPI
                * CINTcommon_fac_sp(envs->i_l) * CINTcommon_fac_sp(envs->j_l)
                * CINTcommon_fac_sp(envs->k_l);

        int *idx = opt->index_xyz_array[envs->i_l*LMAX1*LMAX1
                                       +envs->j_l*LMAX1
                                       +envs->k_l];
        int allocated_idx = 0;
        if (idx == NULL) {
                idx = malloc(sizeof(int) * nf * 3);
                CINTg4c_index_xyz(idx, envs);
                allocated_idx = 1;
        }

        int non0ctrk[k_prim];
        int non0idxk[k_prim*k_ctr];
        CINTOpt_non0coeff_byshell(non0idxk, non0ctrk, ck, k_prim, k_ctr);
        int *non0ctr[3] = {opt->non0ctr[i_sh], opt->non0ctr[j_sh], non0ctrk};
        int *non0idx[3] = {opt->sortedidx[i_sh], opt->sortedidx[j_sh], non0idxk};
        double expcutoff = envs->expcutoff;

        INITSIMD;

        PairData pdata[i_prim*j_prim];
        PairData *pdata_base, *pdata_ij;
        if (opt->data_ptr == NULL) {
                *kempty = CINTset_pairdata(pdata, ai, aj, envs->ri, envs->rj,
                                           envs->li_ceil, envs->lj_ceil,
                                           i_prim, j_prim, SQUARE(envs->rirj), expcutoff);
                if (*kempty) {
                        goto normal_end;
                }
                pdata_base = pdata;
        } else {
                pdata_base = opt->data + opt->data_ptr[i_sh*opt->nbas+j_sh];
        }
        *kempty = 1;
        for (kp = 0; kp < k_prim; kp++) {
                if (k_ctr == 1) {
                        fac1k = common_factor * ck[kp];
                } else {
                        fac1k = common_factor;
                        *jempty = 1;
                }

                pdata_ij = pdata_base;
                for (jp = 0; jp < j_prim; jp++) {
                        if (j_ctr == 1) {
                                fac1j = fac1k * cj[jp];
                        } else {
                                fac1j = fac1k;
                                *iempty = 1;
                        }
                        for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                if (pdata_ij->cceij < expcutoff) {
                                        PUSH(pdata_ij->rij, pdata_ij->eij);
                                }
                        } // end loop i_prim
                        if (!*iempty) {
                                PRIM2CTR(j, gctr[SHLTYPi]);
                        }
                } // end loop j_prim
                if (!*jempty) {
                        PRIM2CTR(k, gctr[SHLTYPj]);
                }
        } // end loop k_prim
        RUN_REST;

        if (n_comp > 1 && !*kempty) {
                int nc = i_ctr * j_ctr * k_ctr;
                CINTdmat_transpose(out, gctr[SHLTYPk], nf*nc, n_comp);
        }
normal_end:
        if (allocated_idx) {
                free(idx);
        }
        return !*kempty;
}


int CINT3c2e_cart_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
                      double *cache)
{
        int *x_ctr = envs->x_ctr;
        int nc = envs->nf * x_ctr[0] * x_ctr[1] * x_ctr[2];
        int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        if (out == NULL) {
                int leng = envs->g_size*3*((1<<envs->gbits)+1)*SIMDD;
                int len0 = envs->nf*n_comp * SIMDD;
                int cache_size = leng + len0 + nc*n_comp*2 + SIMDD*4;
                return cache_size;
        }
        double *stack = NULL;
        if (cache == NULL) {
                int leng = envs->g_size*3*((1<<envs->gbits)+1)*SIMDD;
                int len0 = envs->nf*n_comp * SIMDD;
                int cache_size = leng + len0 + nc*n_comp*2 + SIMDD*4;
                stack = _mm_malloc(sizeof(double)*cache_size, sizeof(double)*SIMDD);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, nc*n_comp);

        int n, has_value;
        if (opt != NULL) {
                has_value = CINT3c2e_loop(gctr, envs, opt, cache);
        } else {
                has_value = CINT3c2e_loop_nopt(gctr, envs, cache);
        }

        int counts[4];
        counts[0] = envs->nfi * x_ctr[0];
        counts[1] = envs->nfj * x_ctr[1];
        counts[2] = envs->nfk * x_ctr[2];
        counts[3] = 1;
        if (dims == NULL) {
                dims = counts;
        }
        int nout = dims[0] * dims[1] * dims[2];
        if (has_value) {
                for (n = 0; n < n_comp; n++) {
                        c2s_cart_3c2e1(out+nout*n, gctr+nc*n, dims, envs, cache);
                }
        } else {
                for (n = 0; n < n_comp; n++) {
                        c2s_dset0(out+nout*n, dims, counts);
                }
        }
        if (stack != NULL) {
                free(stack);
        }
        return has_value;
}
int CINT3c2e_spheric_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
                         double *cache)
{
        int *x_ctr = envs->x_ctr;
        int nc = envs->nf * x_ctr[0] * x_ctr[1] * x_ctr[2];
        int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        if (out == NULL) {
                int leng = envs->g_size*3*((1<<envs->gbits)+1)*SIMDD;
                int len0 = envs->nf*n_comp * SIMDD;
                int cache_size = MAX(leng+len0+nc*n_comp*2, nc*n_comp+envs->nf*3) + SIMDD*4;
                return cache_size;
        }
        double *stack = NULL;
        if (cache == NULL) {
                int leng = envs->g_size*3*((1<<envs->gbits)+1)*SIMDD;
                int len0 = envs->nf*n_comp * SIMDD;
                int cache_size = MAX(leng+len0+nc*n_comp*2, nc*n_comp+envs->nf*3) + SIMDD*4;
                stack = _mm_malloc(sizeof(double)*cache_size, sizeof(double)*SIMDD);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, nc*n_comp);

        int n, has_value;
        if (opt != NULL) {
                has_value = CINT3c2e_loop(gctr, envs, opt, cache);
        } else {
                has_value = CINT3c2e_loop_nopt(gctr, envs, cache);
        }

        int counts[4];
        counts[0] = (envs->i_l*2+1) * x_ctr[0];
        counts[1] = (envs->j_l*2+1) * x_ctr[1];
        counts[2] = (envs->k_l*2+1) * x_ctr[2];
        counts[3] = 1;
        if (dims == NULL) {
                dims = counts;
        }
        int nout = dims[0] * dims[1] * dims[2];
        if (has_value) {
                for (n = 0; n < n_comp; n++) {
                        c2s_sph_3c2e1(out+nout*n, gctr+nc*n, dims, envs, cache);
                }
        } else {
                for (n = 0; n < n_comp; n++) {
                        c2s_dset0(out+nout*n, dims, counts);
                }
        }
        if (stack != NULL) {
                free(stack);
        }
        return has_value;
}
// (spinor,spinor|spherical)
int CINT3c2e_spinor_drv(double complex *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
                        double *cache, void (*f_e1_c2s)())
{
        int *x_ctr = envs->x_ctr;
        int counts[4];
        counts[0] = CINTcgto_spinor(envs->shls[0], envs->bas);
        counts[1] = CINTcgto_spinor(envs->shls[1], envs->bas);
        counts[2] = (envs->k_l*2+1) * x_ctr[2];
        counts[3] = 1;
        int nc = envs->nf * x_ctr[0] * x_ctr[1] * x_ctr[2];
        int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        if (out == NULL) {
                int leng = envs->g_size*3*((1<<envs->gbits)+1)*SIMDD;
                int len0 = envs->nf*n_comp * SIMDD;
                int cache_size = MAX(leng+len0+nc*n_comp*2,
                                     nc*n_comp + envs->nf*14*OF_CMPLX) + SIMDD*4;
                return cache_size;
        }
        double *stack = NULL;
        if (cache == NULL) {
                int leng = envs->g_size*3*((1<<envs->gbits)+1)*SIMDD;
                int len0 = envs->nf*n_comp * SIMDD;
                int cache_size = MAX(leng+len0+nc*n_comp*2,
                                     nc*n_comp + envs->nf*14*OF_CMPLX) + SIMDD*4;
                stack = _mm_malloc(sizeof(double)*cache_size, sizeof(double)*SIMDD);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, nc*n_comp);

        int n, has_value;

        if (opt != NULL) {
                has_value = CINT3c2e_loop(gctr, envs, opt, cache);
        } else {
                has_value = CINT3c2e_loop_nopt(gctr, envs, cache);
        }

        if (dims == NULL) {
                dims = counts;
        }
        int nout = dims[0] * dims[1] * dims[2];
        if (has_value) {
                for (n = 0; n < envs->ncomp_e2 * envs->ncomp_tensor; n++) {
                        (*f_e1_c2s)(out+nout*n, gctr, dims, envs, cache);
                        gctr += nc * envs->ncomp_e1;
                }
        } else {
                for (n = 0; n < envs->ncomp_e2 * envs->ncomp_tensor; n++) {
                        c2s_zset0(out+nout*n, dims, counts);
                }
        }
        if (stack != NULL) {
                free(stack);
        }
        return has_value;
}

int int3c2e_sph(double *out, int *dims, int *shls, int *atm, int natm,
                int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        envs.f_gout_simd1 = &CINTgout2e_simd1;
        return CINT3c2e_spheric_drv(out, dims, &envs, opt, cache);
}
void int3c2e_optimizer(CINTOpt **opt, int *atm, int natm,
                       int *bas, int nbas, double *env)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTall_3c2e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}

int int3c2e_cart(double *out, int *dims, int *shls, int *atm, int natm,
                 int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        envs.f_gout_simd1 = &CINTgout2e_simd1;
        return CINT3c2e_cart_drv(out, dims, &envs, opt, cache);
}

int int3c2e_spinor(double complex *out, int *dims, int *shls, int *atm, int natm,
                   int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        envs.f_gout_simd1 = &CINTgout2e_simd1;
        return CINT3c2e_spinor_drv(out, dims, &envs, opt, cache, &c2s_sf_3c2e1);
}


ALL_CINT(int3c2e)
//ALL_CINT_FORTRAN_(cint3c2e)


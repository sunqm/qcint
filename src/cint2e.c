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
#define SHLTYPl       3

#define ALIAS_ADDR_IF_EQUAL(x, y) \
        if (y##_ctr == 1) { \
                bufctr[SHLTYP##x] = bufctr[SHLTYP##y]; \
                x##empty = y##empty; \
        } else { \
                bufctr[SHLTYP##x] = g1; \
                g1 += len##x; \
        }

#define INIT_GCTR_ADDR(x, y, fac) \
        if (y##_ctr == 1) { \
                fac1##y = (fac) * c##y[y##p]; \
                gctr[SHLTYP##x] = gctr[SHLTYP##y]; \
        } else if (*y##empty || non0ctr##y[y##p] > 1) { \
                fac1##y = (fac); \
                gctr[SHLTYP##x] = bufctr[SHLTYP##x]; \
                *x##empty = 1; \
        } else { \
                fac1##y = (fac) * c##y[y##p+y##_prim*non0idx##y[y##p*y##_ctr]]; \
                gctr[SHLTYP##x] = gctr[SHLTYP##y] + ngp[SHLTYP##y] * non0idx##y[y##p*y##_ctr]; \
                *x##empty = *y##empty; \
        }

#define PRIM2CTR(psymb, csymb) \
        if (csymb##_ctr > 1) {\
                if (*csymb##empty) { \
                        fp2c[np2c] = CINTprim_to_ctr_0; \
                        shltyp[np2c] = SHLTYP##csymb; \
                        gprim[np2c] = gctr[SHLTYP##psymb]; \
                        gp2c [np2c] = gctr[SHLTYP##csymb]; \
                        iprim[np2c] = csymb##p; \
                        np2c++; \
                } else if (non0ctr##csymb[csymb##p] > 1) { \
                        fp2c[np2c] = CINTprim_to_ctr_1; \
                        shltyp[np2c] = SHLTYP##csymb; \
                        gprim[np2c] = gctr[SHLTYP##psymb]; \
                        gp2c [np2c] = gctr[SHLTYP##csymb]; \
                        iprim[np2c] = csymb##p; \
                        np2c++; \
                } \
        } \
        *csymb##empty = 0;

#define POP_PRIM2CTR \
        for (i = 0; i < np2c; i++) { \
                it = shltyp[i]; \
                im = iprim[i]; \
                (*(fp2c[i]))(gp2c[i], gprim[i], coeff[it]+im, \
                             ngp[it], x_prim[it], x_ctr[it], \
                             non0ctr[it][im], non0idx[it]+im*x_ctr[it]); \
        } \
        cum = 0; \
        np2c = 0;

#define PUSH(RIJ, RKL) \
        if (cum == SIMDD) { \
                (*envs->f_g0_2e)(g, &bc, envs, cum); \
                (*envs->f_gout)(gout, g, idx, envs); \
                POP_PRIM2CTR; \
        } \
        envs->ai[cum] = ai[ip]; \
        envs->aj[cum] = aj[jp]; \
        envs->ak[cum] = ak[kp]; \
        envs->al[cum] = al[lp]; \
        envs->rij[0*SIMDD+cum] = *(RIJ+0); \
        envs->rij[1*SIMDD+cum] = *(RIJ+1); \
        envs->rij[2*SIMDD+cum] = *(RIJ+2); \
        envs->rkl[0*SIMDD+cum] = *(RKL+0); \
        envs->rkl[1*SIMDD+cum] = *(RKL+1); \
        envs->rkl[2*SIMDD+cum] = *(RKL+2); \
        fac1i = fac1j * expijkl; \
        envs->fac[cum] = fac1i; \
        if (*iempty) { \
                fp2c[np2c] = CINTiprim_to_ctr_0; \
                *iempty = 0; \
        } else { \
                fp2c[np2c] = CINTiprim_to_ctr_1; \
        } \
        gprim[np2c] = gout + cum * ngp[0]; \
        gp2c [np2c] = gctr[SHLTYPi]; \
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
        double *gp2c [SIMDD*4]; \
        double *gprim[SIMDD*4]; \
        int shltyp[SIMDD*4]; \
        int iprim[SIMDD*4]; \
        void (*fp2c[SIMDD*4])(); \
        MM_STORE(envs->ai, MM_SET1(1.)); \
        MM_STORE(envs->aj, MM_SET1(1.)); \
        MM_STORE(envs->ak, MM_SET1(1.)); \
        MM_STORE(envs->al, MM_SET1(1.)); \
        MM_STORE(envs->fac, MM_SET1(0.));

#define RUN_REST \
        if (cum == 1) { \
                (*envs->f_g0_2e_simd1)(g, &bc, envs, 0); \
                (*envs->f_gout_simd1)(gout, g, idx, envs); \
        } else if (cum > 1) { \
                r1 = MM_SET1(1.); \
                for (i = 0; i < envs->nrys_roots; i++) { \
                        MM_STORE(bc.u+i*SIMDD, r1); \
                        MM_STORE(bc.w+i*SIMDD, r1); \
                } \
                (*envs->f_g0_2e)(g, &bc, envs, cum); \
                (*envs->f_gout)(gout, g, idx, envs); \
        } \
        POP_PRIM2CTR

// little endian on x86
//typedef union {
//    double d;
//    unsigned short s[4];
//} type_IEEE754;
static double approx_log(double x)
{
        //type_IEEE754 y;
        //y.d = x;
        //return ((double)(y.s[3] >> 4) - 1023) * 0.7;
        //return log(x);
        return 2.5;
}

int CINT2e_loop_nopt(double *out, CINTEnvVars *envs, double *cache)
{
        int *shls = envs->shls;
        int *bas = envs->bas;
        double *env = envs->env;
        int i_sh = shls[0];
        int j_sh = shls[1];
        int k_sh = shls[2];
        int l_sh = shls[3];
        int i_ctr = envs->x_ctr[0];
        int j_ctr = envs->x_ctr[1];
        int k_ctr = envs->x_ctr[2];
        int l_ctr = envs->x_ctr[3];
        int i_prim = bas(NPRIM_OF, i_sh);
        int j_prim = bas(NPRIM_OF, j_sh);
        int k_prim = bas(NPRIM_OF, k_sh);
        int l_prim = bas(NPRIM_OF, l_sh);
        int *x_ctr = envs->x_ctr;
        int x_prim[4] = {i_prim, j_prim, k_prim, l_prim};
        double *ai = env + bas(PTR_EXP, i_sh);
        double *aj = env + bas(PTR_EXP, j_sh);
        double *ak = env + bas(PTR_EXP, k_sh);
        double *al = env + bas(PTR_EXP, l_sh);
        double *ci = env + bas(PTR_COEFF, i_sh);
        double *cj = env + bas(PTR_COEFF, j_sh);
        double *ck = env + bas(PTR_COEFF, k_sh);
        double *cl = env + bas(PTR_COEFF, l_sh);
        double *coeff[4] = {ci, cj, ck, cl};
        double *ri = envs->ri;
        double *rj = envs->rj;
        double *rk = envs->rk;
        double *rl = envs->rl;
        int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        int nf = envs->nf;
        double fac1i, fac1j, fac1k, fac1l;
        int ip, jp, kp, lp, ij, i, it, im;
        int empty[4] = {1, 1, 1, 1};
        int *iempty = empty + 0;
        int *jempty = empty + 1;
        int *kempty = empty + 2;
        int *lempty = empty + 3;
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
        int lenl = ALIGN_UP(ngp[3] * l_ctr, SIMDD);

        double *gout, *g, *g1;
        double *gctr[4];
        double *bufctr[4];
        if (n_comp == 1) {
                // patch SIMDD*2 for leni, lenj, lenk with s functions
                MALLOC_INSTACK(g1, lenl+SIMDD*2);
                bufctr[SHLTYPl] = out;
        } else {
                // enlarge out by SIMDD*2, for leni, lenj, lenk with s functions
                cache += SIMDD*2;
                g1 = out;  // Use out as cache for gctrk, gctrj, gctri
                MALLOC_INSTACK(bufctr[SHLTYPl], lenl);
        }
        gctr[SHLTYPl] = bufctr[SHLTYPl];
        ALIAS_ADDR_IF_EQUAL(k, l);
        ALIAS_ADDR_IF_EQUAL(j, k);
        ALIAS_ADDR_IF_EQUAL(i, j);
        MALLOC_INSTACK(gout, len0+MAX(len0,leng));
        g = gout + len0;  // for gx, gy, gz

        ALIGNMM Rys2eT bc;
        int len_ijprim = i_prim * j_prim;
        double rr_ij = SQUARE(envs->rirj);
        double rr_kl = SQUARE(envs->rkrl);
        double log_rr_ij = (envs->li_ceil+envs->lj_ceil+1)*approx_log(rr_ij+1)/2;
        double log_rr_kl = (envs->lk_ceil+envs->ll_ceil+1)*approx_log(rr_kl+1)/2;
        double aij, akl, eij, ekl, ccekl, expijkl;
        ALIGNMM double rij[len_ijprim*4];  // rij[::4] is exp(-eij)
        ALIGNMM double rkl[4];
        double cceij[len_ijprim];
        int non0ctri[i_prim];
        int non0ctrj[j_prim];
        int non0ctrk[k_prim];
        int non0ctrl[l_prim];
        int non0idxi[i_prim*i_ctr];
        int non0idxj[j_prim*j_ctr];
        int non0idxk[k_prim*k_ctr];
        int non0idxl[l_prim*l_ctr];

        *lempty = 1;
        for (ij = 0, jp = 0; jp < j_prim; jp++) {
                for (ip = 0; ip < i_prim; ip++, ij++) {
                        aij = 1/(ai[ip] + aj[jp]);
                        eij = rr_ij * ai[ip] * aj[jp] * aij;
                        for (i = 0; i < 3; i++) {
                                rij[ij*4+i] = (ai[ip]*ri[i] +
                                               aj[jp]*rj[i]) * aij;
                        }
                        cceij[ij] = eij - log_rr_ij;
                        if (cceij[ij] <= CUTOFF15) {
                                *lempty = 0;
                                rij[ij*4+3] = exp(-eij);
                        }
                }
        }
        if (*lempty) {
                goto normal_end;
        }

        int *idx = malloc(sizeof(int) * nf * 3);
        CINTg4c_index_xyz(idx, envs);
        CINTOpt_non0coeff_byshell(non0idxi, non0ctri, coeff[0], i_prim, i_ctr);
        CINTOpt_non0coeff_byshell(non0idxj, non0ctrj, coeff[1], j_prim, j_ctr);
        CINTOpt_non0coeff_byshell(non0idxk, non0ctrk, coeff[2], k_prim, k_ctr);
        CINTOpt_non0coeff_byshell(non0idxl, non0ctrl, coeff[3], l_prim, l_ctr);
        int *non0ctr[4] = {non0ctri, non0ctrj, non0ctrk, non0ctrl};
        int *non0idx[4] = {non0idxi, non0idxj, non0idxk, non0idxl};
        double common_factor = envs->common_factor * (M_PI*M_PI*M_PI)*2/SQRTPI
                * CINTcommon_fac_sp(envs->i_l) * CINTcommon_fac_sp(envs->j_l)
                * CINTcommon_fac_sp(envs->k_l) * CINTcommon_fac_sp(envs->l_l);
        INITSIMD;

        *lempty = 1;
        for (lp = 0; lp < l_prim; lp++) {
                INIT_GCTR_ADDR(k, l, common_factor);
                for (kp = 0; kp < k_prim; kp++) {
                        akl = 1 / (ak[kp] + al[lp]);
                        ekl = rr_kl * ak[kp] * al[lp] * akl;
                        ccekl = ekl - log_rr_kl;
                        if (ccekl > CUTOFF15) {
                                goto k_contracted;
                        }
                        for (i = 0; i < 3; i++) {
                                rkl[i] = (ak[kp]*rk[i] + al[lp]*rl[i]) * akl;
                        }
                        INIT_GCTR_ADDR(j, k, fac1l);
                        ekl = exp(-ekl);

                        for (ij = 0, jp = 0; jp < j_prim; jp++) {
                                INIT_GCTR_ADDR(i, j, fac1k);
                                for (ip = 0; ip < i_prim; ip++, ij++) {
                                        if (cceij[ij] > CUTOFF15 ||
                                            cceij[ij]+ccekl > CUTOFF15) {
                                                goto i_contracted;
                                        }
                                        eij = rij[ij*4+3];
                                        expijkl = eij * ekl;
                                        PUSH(rij+ij*4, rkl);
i_contracted: ;
                                } // end loop i_prim
                                if (!*iempty) {
                                        PRIM2CTR(i, j);
                                }
                        } // end loop j_prim
                        if (!*jempty) {
                                PRIM2CTR(j, k);
                        }
k_contracted: ;
                } // end loop k_prim
                if (!*kempty) {
                        PRIM2CTR(k, l);
                }
        } // end loop l_prim
        RUN_REST;

        if (n_comp > 1 && !*lempty) {
                int nc = i_ctr * j_ctr * k_ctr * l_ctr;
                CINTdmat_transpose(out, gctr[SHLTYPl], nf*nc, n_comp);
        }
        free(idx);
normal_end:
        return !*lempty;
}


int CINT2e_loop(double *out, CINTEnvVars *envs, CINTOpt *opt, double *cache)
{
        int *shls  = envs->shls;
        int i_sh = shls[0];
        int j_sh = shls[1];
        int k_sh = shls[2];
        int l_sh = shls[3];
        if ((opt->data_ptr[i_sh*opt->nbas+j_sh] == NOVALUE) ||
            (opt->data_ptr[k_sh*opt->nbas+l_sh] == NOVALUE)) {
                return 0;
        }
        int *bas = envs->bas;
        double *env = envs->env;
        int i_ctr  = envs->x_ctr[0];
        int j_ctr  = envs->x_ctr[1];
        int k_ctr  = envs->x_ctr[2];
        int l_ctr  = envs->x_ctr[3];
        int i_prim = bas(NPRIM_OF, i_sh);
        int j_prim = bas(NPRIM_OF, j_sh);
        int k_prim = bas(NPRIM_OF, k_sh);
        int l_prim = bas(NPRIM_OF, l_sh);
        int *x_ctr = envs->x_ctr;
        int x_prim[4] = {i_prim, j_prim, k_prim, l_prim};
        double *ai = env + bas(PTR_EXP, i_sh);
        double *aj = env + bas(PTR_EXP, j_sh);
        double *ak = env + bas(PTR_EXP, k_sh);
        double *al = env + bas(PTR_EXP, l_sh);
        double *ci = env + bas(PTR_COEFF, i_sh);
        double *cj = env + bas(PTR_COEFF, j_sh);
        double *ck = env + bas(PTR_COEFF, k_sh);
        double *cl = env + bas(PTR_COEFF, l_sh);
        double *coeff[4] = {ci, cj, ck, cl};
        int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        int nf = envs->nf;
        double fac1i, fac1j, fac1k, fac1l;
        int ip, jp, kp, lp, i, it, im;
        int empty[4] = {1, 1, 1, 1};
        int *iempty = empty + 0;
        int *jempty = empty + 1;
        int *kempty = empty + 2;
        int *lempty = empty + 3;
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
        int lenl = ALIGN_UP(ngp[3] * l_ctr, SIMDD);

        double *gout, *g, *g1;
        double *gctr[4];
        double *bufctr[4];
        if (n_comp == 1) {
                // patch SIMDD*2 for leni, lenj, lenk with s functions
                MALLOC_INSTACK(g1, lenl+SIMDD*2);
                bufctr[SHLTYPl] = out;
        } else {
                // enlarge out by SIMDD*2, for leni, lenj, lenk with s functions
                cache += SIMDD*2;
                g1 = out;  // Use out as cache for gctrk, gctrj, gctri
                MALLOC_INSTACK(bufctr[SHLTYPl], lenl);
        }
        gctr[SHLTYPl] = bufctr[SHLTYPl];
        ALIAS_ADDR_IF_EQUAL(k, l);
        ALIAS_ADDR_IF_EQUAL(j, k);
        ALIAS_ADDR_IF_EQUAL(i, j);
        MALLOC_INSTACK(gout, len0+MAX(len0,leng));
        g = gout + len0;  // for gx, gy, gz

        ALIGNMM Rys2eT bc;
        double common_factor = envs->common_factor * (M_PI*M_PI*M_PI)*2/SQRTPI
                * CINTcommon_fac_sp(envs->i_l) * CINTcommon_fac_sp(envs->j_l)
                * CINTcommon_fac_sp(envs->k_l) * CINTcommon_fac_sp(envs->l_l);
        double expijkl;

        int *idx = opt->index_xyz_array[envs->i_l*LMAX1*LMAX1*LMAX1
                                       +envs->j_l*LMAX1*LMAX1
                                       +envs->k_l*LMAX1
                                       +envs->l_l];
        int allocated_idx = 0;
        if (idx == NULL) {
                idx = malloc(sizeof(int) * nf * 3);
                CINTg4c_index_xyz(idx, envs);
                allocated_idx = 1;
        }
        int *non0ctri = opt->non0ctr[i_sh];
        int *non0ctrj = opt->non0ctr[j_sh];
        int *non0ctrk = opt->non0ctr[k_sh];
        int *non0ctrl = opt->non0ctr[l_sh];
        int *non0idxi = opt->sortedidx[i_sh];
        int *non0idxj = opt->sortedidx[j_sh];
        int *non0idxk = opt->sortedidx[k_sh];
        int *non0idxl = opt->sortedidx[l_sh];
        int *non0ctr[4] = {non0ctri, non0ctrj, non0ctrk, non0ctrl};
        int *non0idx[4] = {non0idxi, non0idxj, non0idxk, non0idxl};

        INITSIMD;

        PairData *pdata_kl, *pdata_ij;
        pdata_kl = opt->data + opt->data_ptr[k_sh*opt->nbas+l_sh];
        for (lp = 0; lp < l_prim; lp++) {
                INIT_GCTR_ADDR(k, l, common_factor);
                for (kp = 0; kp < k_prim; kp++, pdata_kl++) {
                        if (pdata_kl->cceij > CUTOFF15) {
                                goto k_contracted;
                        }
                        INIT_GCTR_ADDR(j, k, fac1l);

                        pdata_ij = opt->data + opt->data_ptr[i_sh*opt->nbas+j_sh];
                        for (jp = 0; jp < j_prim; jp++) {
                                INIT_GCTR_ADDR(i, j, fac1k);
                                for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                        if (pdata_ij->cceij > CUTOFF15 ||
                                            pdata_ij->cceij+pdata_kl->cceij > CUTOFF15) {
                                                goto i_contracted;
                                        }
                                        expijkl = pdata_ij->eij * pdata_kl->eij;
                                        PUSH(pdata_ij->rij, pdata_kl->rij);
i_contracted: ;
                                } // end loop i_prim
                                if (!*iempty) {
                                        PRIM2CTR(i, j);
                                }
                        } // end loop j_prim
                        if (!*jempty) {
                                PRIM2CTR(j, k);
                        }
k_contracted: ;
                } // end loop k_prim
                if (!*kempty) {
                        PRIM2CTR(k, l);
                }
        } // end loop l_prim
        RUN_REST;

        if (n_comp > 1 && !*lempty) {
                int nc = i_ctr * j_ctr * k_ctr * l_ctr;
                CINTdmat_transpose(out, gctr[SHLTYPl], nf*nc, n_comp);
        }
        if (allocated_idx) {
                free(idx);
        }
        return !*lempty;
}


int CINT2e_cart_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
                    double *cache)
{
        int *x_ctr = envs->x_ctr;
        int nc = envs->nf * x_ctr[0] * x_ctr[1] * x_ctr[2] * x_ctr[3];
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
                has_value = CINT2e_loop(gctr, envs, opt, cache);
        } else {
                has_value = CINT2e_loop_nopt(gctr, envs, cache);
        }

        int counts[4];
        counts[0] = envs->nfi * x_ctr[0];
        counts[1] = envs->nfj * x_ctr[1];
        counts[2] = envs->nfk * x_ctr[2];
        counts[3] = envs->nfl * x_ctr[3];
        if (dims == NULL) {
                dims = counts;
        }
        int nout = dims[0] * dims[1] * dims[2] * dims[3];
        if (has_value) {
                for (n = 0; n < n_comp; n++) {
                        c2s_cart_2e1(out+nout*n, gctr+nc*n, dims, envs, cache);
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
int CINT2e_spheric_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
                       double *cache)
{
        int *x_ctr = envs->x_ctr;
        int nc = envs->nf * x_ctr[0] * x_ctr[1] * x_ctr[2] * x_ctr[3];
        int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        if (out == NULL) {
                int leng = envs->g_size*3*((1<<envs->gbits)+1)*SIMDD;
                int len0 = envs->nf*n_comp * SIMDD;
                int cache_size = MAX(leng+len0+nc*n_comp*2, nc*n_comp+envs->nf*4) + SIMDD*4;
                return cache_size;
        }
        double *stack = NULL;
        if (cache == NULL) {
                int leng = envs->g_size*3*((1<<envs->gbits)+1)*SIMDD;
                int len0 = envs->nf*n_comp * SIMDD;
                int cache_size = MAX(leng+len0+nc*n_comp*2, nc*n_comp+envs->nf*4) + SIMDD*4;
                stack = _mm_malloc(sizeof(double)*cache_size, sizeof(double)*SIMDD);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, nc*n_comp);

        int n, has_value;
        if (opt != NULL) {
                has_value = CINT2e_loop(gctr, envs, opt, cache);
        } else {
                has_value = CINT2e_loop_nopt(gctr, envs, cache);
        }

        int counts[4];
        counts[0] = (envs->i_l*2+1) * x_ctr[0];
        counts[1] = (envs->j_l*2+1) * x_ctr[1];
        counts[2] = (envs->k_l*2+1) * x_ctr[2];
        counts[3] = (envs->l_l*2+1) * x_ctr[3];
        if (dims == NULL) {
                dims = counts;
        }
        int nout = dims[0] * dims[1] * dims[2] * dims[3];
        if (has_value) {
                for (n = 0; n < n_comp; n++) {
                        c2s_sph_2e1(out+nout*n, gctr+nc*n, dims, envs, cache);
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
int CINT2e_spinor_drv(double complex *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
                      double *cache, void (*f_e1_c2s)(), void (*f_e2_c2s)())
{
        int *shls = envs->shls;
        int *bas = envs->bas;
        int counts[4];
        counts[0] = CINTcgto_spinor(shls[0], bas);
        counts[1] = CINTcgto_spinor(shls[1], bas);
        counts[2] = CINTcgto_spinor(shls[2], bas);
        counts[3] = CINTcgto_spinor(shls[3], bas);
        int *x_ctr = envs->x_ctr;
        int nc = envs->nf * x_ctr[0] * x_ctr[1] * x_ctr[2] * x_ctr[3];
        int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        int n1 = counts[0] * envs->nfk * x_ctr[2]
                           * envs->nfl * x_ctr[3] * counts[1];
        if (out == NULL) {
                int leng = envs->g_size*3*((1<<envs->gbits)+1)*SIMDD;
                int len0 = envs->nf*n_comp * SIMDD;
                int cache_size = MAX(leng+len0+nc*n_comp*2,
                                     nc*n_comp + n1*envs->ncomp_e2*OF_CMPLX
                                     + envs->nf*32*OF_CMPLX) + SIMDD*4;
                return cache_size;
        }
        double *stack = NULL;
        if (cache == NULL) {
                int leng = envs->g_size*3*((1<<envs->gbits)+1)*SIMDD;
                int len0 = envs->nf*n_comp * SIMDD;
                int cache_size = MAX(leng+len0+nc*n_comp*2,
                                     nc*n_comp + n1*envs->ncomp_e2*OF_CMPLX
                                     + envs->nf*32*OF_CMPLX) + SIMDD*4;
                stack = _mm_malloc(sizeof(double)*cache_size, sizeof(double)*SIMDD);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, nc*n_comp);

        int n, m, has_value;

        if (opt != NULL) {
                has_value = CINT2e_loop(gctr, envs, opt, cache);
        } else {
                has_value = CINT2e_loop_nopt(gctr, envs, cache);
        }

        if (dims == NULL) {
                dims = counts;
        }
        int nout = dims[0] * dims[1] * dims[2] * dims[3];
        if (has_value) {
                double complex *opij;
                MALLOC_INSTACK(opij, n1*envs->ncomp_e2);
                for (n = 0; n < envs->ncomp_tensor; n++) {
                        for (m = 0; m < envs->ncomp_e2; m++) {
                                (*f_e1_c2s)(opij+n1*m, gctr, dims, envs, cache);
                                gctr += nc * envs->ncomp_e1;
                        }
                        (*f_e2_c2s)(out+nout*n, opij, dims, envs, cache);
                }
        } else {
                for (n = 0; n < envs->ncomp_tensor; n++) {
                        c2s_zset0(out+nout*n, dims, counts);
                }
        }
        if (stack != NULL) {
                free(stack);
        }
        return has_value;
}

int int2e_sph(double *out, int *dims, int *shls, int *atm, int natm,
              int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        envs.f_gout_simd1 = &CINTgout2e_simd1;
        return CINT2e_spheric_drv(out, dims, &envs, opt, cache);
}
void int2e_optimizer(CINTOpt **opt, int *atm, int natm,
                     int *bas, int nbas, double *env)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTall_2e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}

int int2e_cart(double *out, int *dims, int *shls, int *atm, int natm,
               int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        envs.f_gout_simd1 = &CINTgout2e_simd1;
        return CINT2e_cart_drv(out, dims, &envs, opt, cache);
}

/*
 * spinor <ki|jl> = (ij|kl); i,j\in electron 1; k,l\in electron 2
 */
int int2e_spinor(double complex *out, int *dims, int *shls, int *atm, int natm,
                 int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        envs.f_gout_simd1 = &CINTgout2e_simd1;
        return CINT2e_spinor_drv(out, dims, &envs, opt, cache,
                                 &c2s_sf_2e1, &c2s_sf_2e2);
}


ALL_CINT(int2e)
ALL_CINT_FORTRAN_(int2e)


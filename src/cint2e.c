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
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include "cint_bas.h"
#include "g2e.h"
#include "optimizer.h"
#include "cint2e.h"
#include "misc.h"
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
                empty_overall = 0; \
        } \
        cum = 0; \
        np2c = 0;

#define POP_PRIM2CTR_AND_SET0 \
        for (i = 0; i < np2c; i++) { \
                it = shltyp[i]; \
                if (it != SHLTYPi) { \
                        im = iprim[i]; \
                        (*(fp2c[i]))(gp2c[i], gprim[i], coeff[it]+im, \
                                     ngp[it], x_prim[it], x_ctr[it], \
                                     non0ctr[it][im], non0idx[it]+im*x_ctr[it]); \
                        empty_overall = 0; \
                } else if (fp2c[i] == CINTiprim_to_ctr_0) { \
                        double *pout = gctr[SHLTYPi]; \
                        int k; \
                        for (k = 0; k < ngp[SHLTYPj]; k++) { \
                                pout[k] = 0.; \
                        } \
                } \
        } \
        cum = 0; \
        np2c = 0;

#define TRANSPOSE(a) \
        if (*empty) { \
                CINTdmat_transpose(out, a, nf*nc, n_comp); \
                *empty = 0; \
        } else { \
                CINTdplus_transpose(out, a, nf*nc, n_comp); \
        } \

#define PUSH(RIJ, RKL) \
        if (cum == SIMDD) { \
                if ((*envs->f_g0_2e)(g, cutoff, &bc, envs, cum)) { \
                        (*envs->f_gout)(gout, g, idx, envs); \
                        POP_PRIM2CTR; \
                } else { \
                        POP_PRIM2CTR_AND_SET0; \
                } \
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
        cutoff[cum] = eijcutoff - pdata_ij->cceij; \
        if (*iempty) { \
                fp2c[np2c] = CINTiprim_to_ctr_0; \
                *iempty = 0; \
        } else { \
                fp2c[np2c] = CINTiprim_to_ctr_1; \
        } \
        gprim[np2c] = gout + cum * ngp[0]; \
        gp2c [np2c] = gctr[SHLTYPi]; \
        iprim[np2c] = ip; \
        shltyp[np2c] = SHLTYPi; \
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
        double *gprim[SIMDD*4]; \
        double *gp2c [SIMDD*4]; \
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

int CINT2e_loop_nopt(double *out, CINTEnvVars *envs, double *cache, int *empty)
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
        double *rk = envs->rk;
        double *rl = envs->rl;
        double *ai = env + bas(PTR_EXP, i_sh);
        double *aj = env + bas(PTR_EXP, j_sh);
        double *ak = env + bas(PTR_EXP, k_sh);
        double *al = env + bas(PTR_EXP, l_sh);
        double *ci = env + bas(PTR_COEFF, i_sh);
        double *cj = env + bas(PTR_COEFF, j_sh);
        double *ck = env + bas(PTR_COEFF, k_sh);
        double *cl = env + bas(PTR_COEFF, l_sh);
        double expcutoff = envs->expcutoff;
        double rr_ij = SQUARE(envs->rirj);
        double rr_kl = SQUARE(envs->rkrl);
        double *log_maxci, *log_maxcj, *log_maxck, *log_maxcl;
        PairData *pdata_base, *pdata_ij;
        MALLOC_DATA_INSTACK(log_maxci, i_prim+j_prim+k_prim+l_prim);
        MALLOC_DATA_INSTACK(pdata_base, i_prim*j_prim);
        log_maxcj = log_maxci + i_prim;
        log_maxck = log_maxcj + j_prim;
        log_maxcl = log_maxck + k_prim;
        CINTOpt_log_max_pgto_coeff(log_maxci, ci, i_prim, i_ctr);
        CINTOpt_log_max_pgto_coeff(log_maxcj, cj, j_prim, j_ctr);
        if (CINTset_pairdata(pdata_base, ai, aj, envs->ri, envs->rj,
                             log_maxci, log_maxcj, envs->li_ceil, envs->lj_ceil,
                             i_prim, j_prim, rr_ij, expcutoff, env)) {
                return 0;
        }
        CINTOpt_log_max_pgto_coeff(log_maxck, ck, k_prim, k_ctr);
        CINTOpt_log_max_pgto_coeff(log_maxcl, cl, l_prim, l_ctr);

        int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        size_t nf = envs->nf;
        double fac1i, fac1j, fac1k, fac1l;
        int ip, jp, kp, lp, i, it, im;
        int _empty[4] = {1, 1, 1, 1};
        int *iempty = _empty + 0;
        int *jempty = _empty + 1;
        int *kempty = _empty + 2;
        int *lempty = _empty + 3;
        int empty_overall = 1;
        int lkl = envs->lk_ceil + envs->ll_ceil;
        double akl, ekl, expijkl, ccekl, log_rr_kl, eijcutoff;

        akl = ak[k_prim-1] + al[l_prim-1];
        log_rr_kl = 1.7 - 1.5 * approx_log(akl);
        double omega = env[PTR_RANGE_OMEGA];
        if (omega < 0) {
                // Normally the factor
                //    (aj*d/aij+theta*R)^li * (ai*d/aij+theta*R)^lj * pi^1.5/aij^{(li+lj+3)/2}
                // is a good approximation for polynomial parts in SR-ERIs.
                //    <~ (aj*d/aij+theta*R)^li * (ai*d/aij+theta*R)^lj * (pi/aij)^1.5
                //    <~ (d+theta*R)^li * (d+theta*R)^lj * (pi/aij)^1.5
                if (envs->rys_order > 1) {
                        double r_guess = 8.;
                        double omega2 = omega * omega;
                        int lij = envs->li_ceil + envs->lj_ceil;
                        if (lij > 0) {
                                double aij = ai[i_prim-1] + aj[j_prim-1];
                                double dist_ij = sqrt(rr_ij);
                                double theta = omega2 / (omega2 + aij);
                                expcutoff += lij * approx_log(
                                        (dist_ij+theta*r_guess+1.)/(dist_ij+1.));
                        }
                        if (lkl > 0) {
                                double theta = omega2 / (omega2 + akl);
                                log_rr_kl += lkl * approx_log(
                                        sqrt(rr_kl) + theta*r_guess + 1.);
                        }
                }
        } else {
                if (lkl > 0) {
                        log_rr_kl += lkl * approx_log(sqrt(rr_kl) + 1.);
                }
        }

        int *idx;
        MALLOC_DATA_INSTACK(idx, nf * 3);
        CINTg4c_index_xyz(idx, envs);

        double *coeff[4] = {ci, cj, ck, cl};
        int *non0ctri, *non0ctrj, *non0ctrk, *non0ctrl;
        int *non0idxi, *non0idxj, *non0idxk, *non0idxl;
        MALLOC_DATA_INSTACK(non0ctri, i_prim+j_prim+k_prim+l_prim+i_prim*i_ctr+j_prim*j_ctr+k_prim*k_ctr+l_prim*l_ctr);
        non0ctrj = non0ctri + i_prim;
        non0ctrk = non0ctrj + j_prim;
        non0ctrl = non0ctrk + k_prim;
        non0idxi = non0ctrl + l_prim;
        non0idxj = non0idxi + i_prim*i_ctr;
        non0idxk = non0idxj + j_prim*j_ctr;
        non0idxl = non0idxk + k_prim*k_ctr;
        CINTOpt_non0coeff_byshell(non0idxi, non0ctri, coeff[0], i_prim, i_ctr);
        CINTOpt_non0coeff_byshell(non0idxj, non0ctrj, coeff[1], j_prim, j_ctr);
        CINTOpt_non0coeff_byshell(non0idxk, non0ctrk, coeff[2], k_prim, k_ctr);
        CINTOpt_non0coeff_byshell(non0idxl, non0ctrl, coeff[3], l_prim, l_ctr);
        int *non0ctr[4] = {non0ctri, non0ctrj, non0ctrk, non0ctrl};
        int *non0idx[4] = {non0idxi, non0idxj, non0idxk, non0idxl};
        double common_factor = envs->common_factor;

        size_t ngp[4];
        ngp[0] = nf * n_comp;
        ngp[1] = ngp[0] * i_ctr;
        ngp[2] = ngp[1] * j_ctr;
        ngp[3] = ngp[2] * k_ctr;
        // (irys,i,j,k,l,coord,0:1); +1 for nabla-r12
        int leng = envs->g_size * 3 * ((1<<envs->gbits)+1) * SIMDD;
        size_t len0 = ngp[0] * SIMDD;
        size_t leni = ALIGN_UP(ngp[1], SIMDD);
        size_t lenj = ALIGN_UP(ngp[2], SIMDD);
        size_t lenk = ALIGN_UP(ngp[3], SIMDD);
        size_t lenl = ALIGN_UP(ngp[3] * l_ctr, SIMDD);

        double *gout, *g, *g1;
        double *gctr[4];
        double *bufctr[4];
        if (n_comp == 1) {
                // patch SIMDD*2 for leni, lenj, lenk with s functions
                MALLOC_INSTACK(g1, lenl+SIMDD*2);
                bufctr[SHLTYPl] = out;
                *lempty = *empty;
        } else {
                // enlarge sizeo of out by SIMDD*2, for leni, lenj, lenk with s functions
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
        ALIGNMM double cutoff[SIMDD];
        ALIGNMM double rkl[4];

        INITSIMD;

        for (lp = 0; lp < l_prim; lp++) {
                INIT_GCTR_ADDR(k, l, common_factor);
                for (kp = 0; kp < k_prim; kp++) {
                        akl = 1 / (ak[kp] + al[lp]);
                        ekl = rr_kl * ak[kp] * al[lp] * akl;
                        ccekl = ekl - log_rr_kl - log_maxck[kp] - log_maxcl[lp];
                        if (ccekl > expcutoff) {
                                goto k_contracted;
                        }
                        eijcutoff = expcutoff - MAX(ccekl, 0);
                        rkl[0] = (ak[kp]*rk[0] + al[lp]*rl[0]) * akl;
                        rkl[1] = (ak[kp]*rk[1] + al[lp]*rl[1]) * akl;
                        rkl[2] = (ak[kp]*rk[2] + al[lp]*rl[2]) * akl;
                        INIT_GCTR_ADDR(j, k, fac1l);
                        ekl = exp(-ekl);

                        pdata_ij = pdata_base;
                        for (jp = 0; jp < j_prim; jp++) {
                                INIT_GCTR_ADDR(i, j, fac1k);
                                for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                        if (pdata_ij->cceij > eijcutoff) {
                                                goto i_contracted;
                                        }
                                        expijkl = pdata_ij->eij * ekl;
                                        PUSH(pdata_ij->rij, rkl);
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

        if (n_comp > 1 && !empty_overall) {
                int nc = i_ctr * j_ctr * k_ctr * l_ctr;
                TRANSPOSE(gctr[SHLTYPl]);
        }
        *empty &= empty_overall;
        return !empty_overall;
}


int CINT2e_loop(double *out, CINTEnvVars *envs, double *cache, int *empty)
{
        int *shls  = envs->shls;
        int i_sh = shls[0];
        int j_sh = shls[1];
        int k_sh = shls[2];
        int l_sh = shls[3];
        CINTOpt *opt = envs->opt;
        if (opt->pairdata != NULL &&
            ((opt->pairdata[i_sh*opt->nbas+j_sh] == NOVALUE) ||
             (opt->pairdata[k_sh*opt->nbas+l_sh] == NOVALUE))) {
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
        double expcutoff = envs->expcutoff;
        PairData *_pdata_ij, *_pdata_kl, *pdata_kl, *pdata_ij;
        if (opt->pairdata != NULL) {
                _pdata_ij = opt->pairdata[i_sh*opt->nbas+j_sh];
                _pdata_kl = opt->pairdata[k_sh*opt->nbas+l_sh];
        } else {
                double *log_maxci = opt->log_max_coeff[i_sh];
                double *log_maxcj = opt->log_max_coeff[j_sh];
                MALLOC_DATA_INSTACK(_pdata_ij, i_prim*j_prim + k_prim*l_prim);
                if (CINTset_pairdata(_pdata_ij, ai, aj, envs->ri, envs->rj,
                                     log_maxci, log_maxcj, envs->li_ceil, envs->lj_ceil,
                                     i_prim, j_prim, SQUARE(envs->rirj), expcutoff, env)) {
                        return 0;
                }

                double *log_maxck = opt->log_max_coeff[k_sh];
                double *log_maxcl = opt->log_max_coeff[l_sh];
                _pdata_kl = _pdata_ij + i_prim*j_prim;
                if (CINTset_pairdata(_pdata_kl, ak, al, envs->rk, envs->rl,
                                     log_maxck, log_maxcl, envs->lk_ceil, envs->ll_ceil,
                                     k_prim, l_prim, SQUARE(envs->rkrl), expcutoff, env)) {
                        return 0;
                }
        }

        int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        size_t nf = envs->nf;
        double fac1i, fac1j, fac1k, fac1l;
        int ip, jp, kp, lp, i, it, im;
        int _empty[4] = {1, 1, 1, 1};
        int *iempty = _empty + 0;
        int *jempty = _empty + 1;
        int *kempty = _empty + 2;
        int *lempty = _empty + 3;
        int empty_overall = 1;

        int *idx = opt->index_xyz_array[envs->i_l*LMAX1*LMAX1*LMAX1
                                       +envs->j_l*LMAX1*LMAX1
                                       +envs->k_l*LMAX1
                                       +envs->l_l];
        if (idx == NULL) {
                MALLOC_DATA_INSTACK(idx, nf * 3);
                CINTg4c_index_xyz(idx, envs);
        }

        size_t ngp[4];
        ngp[0] = nf * n_comp;
        ngp[1] = ngp[0] * i_ctr;
        ngp[2] = ngp[1] * j_ctr;
        ngp[3] = ngp[2] * k_ctr;
        // (irys,i,j,k,l,coord,0:1); +1 for nabla-r12
        int leng = envs->g_size * 3 * ((1<<envs->gbits)+1) * SIMDD;
        size_t len0 = ngp[0] * SIMDD;
        size_t leni = ALIGN_UP(ngp[1], SIMDD);
        size_t lenj = ALIGN_UP(ngp[2], SIMDD);
        size_t lenk = ALIGN_UP(ngp[3], SIMDD);
        size_t lenl = ALIGN_UP(ngp[3] * l_ctr, SIMDD);
        double *gout, *g, *g1;
        double *gctr[4];
        double *bufctr[4];
        if (n_comp == 1) {
                // patch SIMDD*2 for leni, lenj, lenk with s functions
                MALLOC_INSTACK(g1, lenl+SIMDD*2);
                bufctr[SHLTYPl] = out;
                *lempty = *empty;
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

        double *coeff[4] = {ci, cj, ck, cl};
        ALIGNMM Rys2eT bc;
        double common_factor = envs->common_factor;
        double expijkl, eijcutoff;

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
        ALIGNMM double cutoff[SIMDD];

        INITSIMD;

        pdata_kl = _pdata_kl;
        for (lp = 0; lp < l_prim; lp++) {
                INIT_GCTR_ADDR(k, l, common_factor);
                for (kp = 0; kp < k_prim; kp++, pdata_kl++) {
                        if (pdata_kl->cceij > expcutoff) {
                                goto k_contracted;
                        }
                        INIT_GCTR_ADDR(j, k, fac1l);

                        eijcutoff = expcutoff - MAX(pdata_kl->cceij, 0);
                        pdata_ij = _pdata_ij;
                        for (jp = 0; jp < j_prim; jp++) {
                                INIT_GCTR_ADDR(i, j, fac1k);
                                for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                        if (pdata_ij->cceij > eijcutoff) {
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

        if (n_comp > 1 && !empty_overall) {
                int nc = i_ctr * j_ctr * k_ctr * l_ctr;
                TRANSPOSE(gctr[SHLTYPl]);
        }
        *empty &= empty_overall;
        return !empty_overall;
}

int (*CINT2e_1111_loop)(double *, CINTEnvVars *, double *, int *) = &CINT2e_loop;

#define PAIRDATA_NON0IDX_SIZE(ps) \
                int *bas = envs->bas; \
                int *shls  = envs->shls; \
                int i_prim = bas(NPRIM_OF, shls[0]); \
                int j_prim = bas(NPRIM_OF, shls[1]); \
                int k_prim = bas(NPRIM_OF, shls[2]); \
                int l_prim = bas(NPRIM_OF, shls[3]); \
                size_t ps = ((i_prim*j_prim + k_prim*l_prim) * 5 \
                          + i_prim * x_ctr[0] \
                          + j_prim * x_ctr[1] \
                          + k_prim * x_ctr[2] \
                          + l_prim * x_ctr[3] \
                          +(i_prim+j_prim+k_prim+l_prim)*2 + nf*3);

CACHE_SIZE_T CINT2e_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
                        double *cache, void (*f_c2s)())
{
        int *x_ctr = envs->x_ctr;
        size_t nf = envs->nf;
        size_t nc = nf * x_ctr[0] * x_ctr[1] * x_ctr[2] * x_ctr[3];
        int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        if (out == NULL) {
                PAIRDATA_NON0IDX_SIZE(pdata_size);
                size_t leng = envs->g_size*3*((1<<envs->gbits)+1)*SIMDD;
                size_t len0 = nf*n_comp * SIMDD;
                size_t cache_size = MAX(leng+len0+nc*n_comp*3 + pdata_size,
                                        nc*n_comp+nf*4) + SIMDD*4;
#ifndef CACHE_SIZE_I8
                if (cache_size >= INT32_MAX) {
                        fprintf(stderr, "CINT2e_drv cache_size overflow: "
                                "cache_size %zu > %d, nf %zu, nc %zu, n_comp %d\n",
                                cache_size, INT32_MAX, nf, nc, n_comp);
                        cache_size = 0;
                }
#endif
                return cache_size;
        }
        double *stack = NULL;
        if (cache == NULL) {
                PAIRDATA_NON0IDX_SIZE(pdata_size);
                size_t leng = envs->g_size*3*((1<<envs->gbits)+1)*SIMDD;
                size_t len0 = nf*n_comp * SIMDD;
                size_t cache_size = MAX(leng+len0+nc*n_comp*3 + pdata_size,
                                        nc*n_comp+nf*4) + SIMDD*4;
                stack = _mm_malloc(sizeof(double)*cache_size, sizeof(double)*SIMDD);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, nc*n_comp);

        int n;
        int empty = 1;
        if (opt != NULL) {
                envs->opt = opt;
                CINT2e_loop(gctr, envs, cache, &empty);
        } else {
                CINT2e_loop_nopt(gctr, envs, cache, &empty);
        }

        int counts[4];
        if (f_c2s == &c2s_sph_2e1) {
                counts[0] = (envs->i_l*2+1) * x_ctr[0];
                counts[1] = (envs->j_l*2+1) * x_ctr[1];
                counts[2] = (envs->k_l*2+1) * x_ctr[2];
                counts[3] = (envs->l_l*2+1) * x_ctr[3];
        } else {
                counts[0] = envs->nfi * x_ctr[0];
                counts[1] = envs->nfj * x_ctr[1];
                counts[2] = envs->nfk * x_ctr[2];
                counts[3] = envs->nfl * x_ctr[3];
        }
        if (dims == NULL) {
                dims = counts;
        }
        int nout = dims[0] * dims[1] * dims[2] * dims[3];
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
CACHE_SIZE_T CINT2e_spinor_drv(double complex *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
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
        size_t nf = envs->nf;
        size_t nc = nf * x_ctr[0] * x_ctr[1] * x_ctr[2] * x_ctr[3];
        int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        int n1 = counts[0] * envs->nfk * x_ctr[2]
                           * envs->nfl * x_ctr[3] * counts[1];
        if (out == NULL) {
                PAIRDATA_NON0IDX_SIZE(pdata_size);
                size_t leng = envs->g_size*3*((1<<envs->gbits)+1)*SIMDD;
                size_t len0 = nf*n_comp * SIMDD;
                size_t cache_size = MAX(leng+len0+nc*n_comp*3 + pdata_size,
                                     nc*n_comp + n1*envs->ncomp_e2*OF_CMPLX
                                     + nf*32*OF_CMPLX) + SIMDD*4;
#ifndef CACHE_SIZE_I8
                if (cache_size >= INT32_MAX) {
                        fprintf(stderr, "CINT2e_spinor_drv cache_size overflow: "
                                "cache_size %zu > %d, nf %zu, nc %zu, n_comp %d\n",
                                cache_size, INT32_MAX, nf, nc, n_comp);
                        cache_size = 0;
                }
#endif
                return cache_size;
        }
        double *stack = NULL;
        if (cache == NULL) {
                PAIRDATA_NON0IDX_SIZE(pdata_size);
                size_t leng = envs->g_size*3*((1<<envs->gbits)+1)*SIMDD;
                size_t len0 = nf*n_comp * SIMDD;
                size_t cache_size = MAX(leng+len0+nc*n_comp*3 + pdata_size,
                                     nc*n_comp + n1*envs->ncomp_e2*OF_CMPLX
                                     + nf*32*OF_CMPLX) + SIMDD*4;
                stack = _mm_malloc(sizeof(double)*cache_size, sizeof(double)*SIMDD);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, nc*n_comp);

        int n, m;
        int empty = 1;
        if (opt != NULL) {
                envs->opt = opt;
                CINT2e_loop(gctr, envs, cache, &empty);
        } else {
                CINT2e_loop_nopt(gctr, envs, cache, &empty);
        }

        if (dims == NULL) {
                dims = counts;
        }
        int nout = dims[0] * dims[1] * dims[2] * dims[3];
        if (!empty) {
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
        return !empty;
}

CACHE_SIZE_T int2e_sph(double *out, int *dims, int *shls, int *atm, int natm,
              int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        envs.f_gout_simd1 = &CINTgout2e_simd1;
        return CINT2e_drv(out, dims, &envs, opt, cache, &c2s_sph_2e1);
}
void int2e_optimizer(CINTOpt **opt, int *atm, int natm,
                     int *bas, int nbas, double *env)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTall_2e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}

CACHE_SIZE_T int2e_cart(double *out, int *dims, int *shls, int *atm, int natm,
               int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        envs.f_gout_simd1 = &CINTgout2e_simd1;
        return CINT2e_drv(out, dims, &envs, opt, cache, &c2s_cart_2e1);
}

/*
 * spinor <ki|jl> = (ij|kl); i,j\in electron 1; k,l\in electron 2
 */
CACHE_SIZE_T int2e_spinor(double complex *out, int *dims, int *shls, int *atm, int natm,
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


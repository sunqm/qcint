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
 *
 * 4-center 1-electron integrals
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "cint_bas.h"
#include "misc.h"
#include "g1e.h"
#include "g2e.h"
#include "optimizer.h"
#include "cint1e.h"
#include "cart2sph.h"
#include "c2f.h"

#define SHLTYPi       0
#define SHLTYPj       1
#define SHLTYPk       2
#define SHLTYPl       3

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

#define PUSH(RIJ, RKL, EIJKL) \
        if (cum == SIMDD) { \
                for (i = 0; i < cum; i++) { \
                        envs->fac[i] *= exp(-eijks[i]); \
                } \
                CINTg4c1e_ovlp(g, bufg, envs, cum); \
                (*envs->f_gout)(gout, g, idx, envs, cum); \
                POP_PRIM2CTR; \
        } \
        envs->ai[cum] = ai[ip]; \
        envs->aj[cum] = aj[jp]; \
        envs->ak[cum] = ak[kp]; \
        envs->al[cum] = al[lp]; \
        envs->rij[0*SIMDD+cum] = (RIJ)[0]; \
        envs->rij[1*SIMDD+cum] = (RIJ)[1]; \
        envs->rij[2*SIMDD+cum] = (RIJ)[2]; \
        envs->rkl[0*SIMDD+cum] = (RKL)[0]; \
        envs->rkl[1*SIMDD+cum] = (RKL)[1]; \
        envs->rkl[2*SIMDD+cum] = (RKL)[2]; \
        envs->fac[cum] = fac1j; \
        eijks[cum] = EIJKL; \
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
        int cum = 0; \
        int np2c = 0; \
        double eijks[SIMDD]; \
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
        if (cum > 0) { \
                for (i = 0; i < cum; i++) { \
                        envs->fac[i] *= exp(-eijks[i]); \
                } \
                CINTg4c1e_ovlp(g, bufg, envs, cum); \
                (*envs->f_gout)(gout, g, idx, envs, cum); \
        } \
        POP_PRIM2CTR

#define TRANSPOSE(a) \
        if (*empty) { \
                CINTdmat_transpose(out, a, nf*nc, n_comp); \
        } else { \
                CINTdplus_transpose(out, a, nf*nc, n_comp); \
        } \
        *empty = 0;

void CINTg4c1e_ovlp(double *g, double *buf, CINTEnvVars *envs, int count);
void CINTinit_int4c1e_EnvVars(CINTEnvVars *envs, int *ng, FINT *shls,
                              FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env);

int CINT4c1e_loop_nopt(double *out, CINTEnvVars *envs, double *cache, int *empty)
{
        FINT *shls = envs->shls;
        FINT *bas = envs->bas;
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
        size_t nf = envs->nf;
        double fac1j, fac1k, fac1l;
        int ip, jp, kp, lp, i, it, im;
        int _empty[4] = {1, 1, 1, 1};
        int *iempty = _empty + 0;
        int *jempty = _empty + 1;
        int *kempty = _empty + 2;
        int *lempty = _empty + 3;

        int *idx;
        MALLOC_DATA_INSTACK(idx, nf * 3);
        CINTg4c_index_xyz(idx, envs);

        int *non0ctri, *non0ctrj, *non0ctrk, *non0ctrl;
        int *non0idxi, *non0idxj, *non0idxk, *non0idxl;
        MALLOC_DATA_INSTACK(non0ctri, i_prim+j_prim+k_prim+l_prim
                            +i_prim*i_ctr+j_prim*j_ctr+k_prim*k_ctr+l_prim*l_ctr);
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

        double common_factor = envs->common_factor * SQRTPI * M_PI
                * CINTcommon_fac_sp(envs->i_l) * CINTcommon_fac_sp(envs->j_l)
                * CINTcommon_fac_sp(envs->k_l) * CINTcommon_fac_sp(envs->l_l);
        double expcutoff = envs->expcutoff;

        size_t ngp[4];
        ngp[0] = nf * n_comp;
        ngp[1] = ngp[0] * i_ctr;
        ngp[2] = ngp[1] * j_ctr;
        ngp[3] = ngp[2] * k_ctr;
        int leng = envs->g_size * 3 * ((1<<envs->gbits)+1) * SIMDD;
        size_t len0 = ngp[0] * SIMDD;
        size_t leni = ALIGN_UP(ngp[1], SIMDD);
        size_t lenj = ALIGN_UP(ngp[2], SIMDD);
        size_t lenk = ALIGN_UP(ngp[3], SIMDD);
        size_t lenl = ALIGN_UP(ngp[3] * l_ctr, SIMDD);

        // workspace for the base 2d recurrences in CINTg4c1e_ovlp
        int nmax = envs->li_ceil + envs->lj_ceil;
        int mmax = envs->lk_ceil + envs->ll_ceil;
        int b_size = (nmax+mmax+1) * (MAX(nmax,mmax)+1);
        double *bufg;
        MALLOC_INSTACK(bufg, b_size * 3 * SIMDD);

        double *gout, *g, *g1;
        double *gctr[4];
        if (n_comp == 1) {
                // patch SIMDD for leni, lenj, lenk with s functions
                MALLOC_INSTACK(g1, lenl+SIMDD);
                gctr[SHLTYPl] = out;
                lempty = empty;
        } else {
                // enlarge out by SIMDD, for leni, lenj, lenk with s functions
                cache += SIMDD;
                g1 = out;  // Use out as cache for gctrk, gctrj, gctri
                MALLOC_INSTACK(gctr[SHLTYPl], lenl);
        }
        ALIAS_ADDR_IF_EQUAL(k, l);
        ALIAS_ADDR_IF_EQUAL(j, k);
        ALIAS_ADDR_IF_EQUAL(i, j);
        MALLOC_INSTACK(gout, leng+len0);
        g = gout + len0;  // for gx, gy, gz

        double aij, akl, a0, eij, ekl, eijkl;
        double rijrkl[3];
        ALIGNMM double rij[4];
        ALIGNMM double rkl[4];
        const double rr_ij = SQUARE(envs->rirj);
        const double rr_kl = SQUARE(envs->rkrl);
        INITSIMD;

        for (lp = 0; lp < l_prim; lp++) {
                if (l_ctr == 1) {
                        fac1l = common_factor * cl[lp];
                } else {
                        fac1l = common_factor;
                        *kempty = 1;
                }
                for (kp = 0; kp < k_prim; kp++) {
                        akl = ak[kp] + al[lp];
                        ekl = rr_kl * ak[kp] * al[lp] / akl;
                        if (ekl > expcutoff) {
                                goto k_contracted;
                        }
                        rkl[0] = (ak[kp]*rk[0] + al[lp]*rl[0]) / akl;
                        rkl[1] = (ak[kp]*rk[1] + al[lp]*rl[1]) / akl;
                        rkl[2] = (ak[kp]*rk[2] + al[lp]*rl[2]) / akl;
                        if (k_ctr == 1) {
                                fac1k = fac1l * ck[kp];
                        } else {
                                fac1k = fac1l;
                                *jempty = 1;
                        }

                        for (jp = 0; jp < j_prim; jp++) {
                                if (j_ctr == 1) {
                                        fac1j = fac1k * cj[jp];
                                } else {
                                        fac1j = fac1k;
                                        *iempty = 1;
                                }
                                for (ip = 0; ip < i_prim; ip++) {
                                        aij = ai[ip] + aj[jp];
                                        eij = rr_ij * ai[ip] * aj[jp] / aij;
                                        if (eij+ekl > expcutoff) {
                                                goto i_contracted;
                                        }
                                        rij[0] = (ai[ip]*ri[0] + aj[jp]*rj[0]) / aij;
                                        rij[1] = (ai[ip]*ri[1] + aj[jp]*rj[1]) / aij;
                                        rij[2] = (ai[ip]*ri[2] + aj[jp]*rj[2]) / aij;
                                        rijrkl[0] = rij[0] - rkl[0];
                                        rijrkl[1] = rij[1] - rkl[1];
                                        rijrkl[2] = rij[2] - rkl[2];
                                        a0 = aij * akl / (aij + akl);
                                        eijkl = a0 * SQUARE(rijrkl) + eij + ekl;
                                        if (eijkl > expcutoff) {
                                                goto i_contracted;
                                        }
                                        PUSH(rij, rkl, eijkl);
i_contracted: ;
                                } // end loop i_prim
                                if (!*iempty) {
                                        PRIM2CTR(j, gctr[SHLTYPi]);
                                }
                        } // end loop j_prim
                        if (!*jempty) {
                                PRIM2CTR(k, gctr[SHLTYPj]);
                        }
k_contracted: ;
                } // end loop k_prim
                if (!*kempty) {
                        PRIM2CTR(l, gctr[SHLTYPk]);
                }
        } // end loop l_prim
        RUN_REST;

        if (n_comp > 1 && !*lempty) {
                int nc = i_ctr * j_ctr * k_ctr * l_ctr;
                TRANSPOSE(gctr[SHLTYPl]);
        }
        return !*empty;
}


#define PAIRDATA_NON0IDX_SIZE(ps) \
                FINT *bas = envs->bas; \
                FINT *shls  = envs->shls; \
                int i_prim = bas(NPRIM_OF, shls[0]); \
                int j_prim = bas(NPRIM_OF, shls[1]); \
                int k_prim = bas(NPRIM_OF, shls[2]); \
                int l_prim = bas(NPRIM_OF, shls[3]); \
                size_t ps = (i_prim * x_ctr[0] \
                           + j_prim * x_ctr[1] \
                           + k_prim * x_ctr[2] \
                           + l_prim * x_ctr[3] \
                           + (i_prim+j_prim+k_prim+l_prim)*2 \
                           + envs->nf*3);

CACHE_SIZE_T CINT4c1e_drv(double *out, FINT *dims, CINTEnvVars *envs, CINTOpt *opt,
                          double *cache, void (*f_c2s)())
{
        int *x_ctr = envs->x_ctr;
        size_t nf = envs->nf;
        size_t nc = nf * x_ctr[0] * x_ctr[1] * x_ctr[2] * x_ctr[3];
        int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        int nmax = envs->li_ceil + envs->lj_ceil;
        int mmax = envs->lk_ceil + envs->ll_ceil;
        size_t bufg_size = (size_t)(nmax+mmax+1) * (MAX(nmax,mmax)+1) * 3 * SIMDD;
        if (out == NULL) {
                PAIRDATA_NON0IDX_SIZE(pdata_size);
                size_t leng = envs->g_size*3*((1<<envs->gbits)+1)*SIMDD;
                size_t len0 = nf*n_comp * SIMDD;
                size_t cache_size = MAX(leng+len0+nc*n_comp*3+pdata_size+bufg_size,
                                        nc*n_comp+nf*4) + SIMDD*4;
#ifndef CACHE_SIZE_I8
                if (cache_size >= INT32_MAX) {
                        fprintf(stderr, "CINT4c1e_drv cache_size overflow: "
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
                size_t cache_size = MAX(leng+len0+nc*n_comp*3+pdata_size+bufg_size,
                                        nc*n_comp+nf*4) + SIMDD*4;
                stack = _mm_malloc(sizeof(double)*cache_size, sizeof(double)*SIMDD);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, nc*n_comp);

        int n;
        int empty = 1;
        CINT4c1e_loop_nopt(gctr, envs, cache, &empty);

        FINT counts[4];
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

void CINTgout4c1e(double *gout, double *g, int *idx, CINTEnvVars *envs, int count)
{
        // CINTg4c1e_ovlp has been called by the primitive loop
        int nf = envs->nf;
        int nfc = nf;
        int n;
        double *gx, *gy, *gz;
        DECLARE_GOUT;
        __MD r0;
        for (n = 0; n < nf; n++) {
                gx = g + idx[n*3+0] * SIMDD;
                gy = g + idx[n*3+1] * SIMDD;
                gz = g + idx[n*3+2] * SIMDD;
                r0 = MM_LOAD(gx) * MM_LOAD(gy) * MM_LOAD(gz);
                GOUT_SCATTER(gout, n, r0);
        }
}

CACHE_SIZE_T int4c1e_sph(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                         FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int4c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout4c1e;
        return CINT4c1e_drv(out, dims, &envs, opt, cache, &c2s_sph_2e1);
}
void int4c1e_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                       FINT *bas, FINT nbas, double *env)
{
        *opt = NULL;
}

CACHE_SIZE_T int4c1e_cart(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                          FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int4c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout4c1e;
        return CINT4c1e_drv(out, dims, &envs, opt, cache, &c2s_cart_2e1);
}

CACHE_SIZE_T int4c1e_spinor(double complex *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                            FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        fprintf(stderr, "int4c1e_spinor not implemented\n");
        exit(1);
}

ALL_CINT(int4c1e)
ALL_CINT_FORTRAN_(int4c1e)

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

#define PUSH(EIJK) \
        if (cum == SIMDD) { \
                for (i = 0; i < cum; i++) { \
                        envs->fac[i] *= exp(-eijks[i]); \
                } \
                (*envs->f_gout)(gout, g, idx, envs, cum); \
                POP_PRIM2CTR; \
        } \
        envs->ai[cum] = ai[ip]; \
        envs->aj[cum] = aj[jp]; \
        envs->ak[cum] = ak[kp]; \
        envs->fac[cum] = fac1j; \
        eijks[cum] = EIJK; \
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
        double *gprim[SIMDD*3]; \
        int shltyp[SIMDD*3]; \
        int iprim[SIMDD*3]; \
        void (*fp2c[SIMDD*3])(); \
        MM_STORE(envs->ai, MM_SET1(1.)); \
        MM_STORE(envs->aj, MM_SET1(1.)); \
        MM_STORE(envs->ak, MM_SET1(1.)); \
        MM_STORE(envs->fac, MM_SET1(0.));

#define RUN_REST \
        if (cum > 0) { \
                for (i = 0; i < cum; i++) { \
                        envs->fac[i] *= exp(-eijks[i]); \
                } \
                (*envs->f_gout)(gout, g, idx, envs, cum); \
        } \
        POP_PRIM2CTR

int int3c1e_cache_size(CINTEnvVars *envs);
void CINTg3c1e_ovlp(double *g, CINTEnvVars *envs, int count);

int CINT3c1e_loop_nopt(double *out, CINTEnvVars *envs, double *cache)
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
        double *ri = envs->ri;
        double *rk = envs->rk;
        int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        int nf = envs->nf;
        double fac1j, fac1k;
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

        double aijk, eijk;
        double aiajrr, aiakrr, ajakrr;
        double rirk[3];
        rirk[0] = ri[0] - rk[0];
        rirk[1] = ri[1] - rk[1];
        rirk[2] = ri[2] - rk[2];
        const double rr_ij = SQUARE(envs->rirj);
        const double rr_ik = SQUARE(      rirk);
        const double rr_jk = SQUARE(envs->rkrl);
        int idx[nf*3];
        int non0ctri[i_prim];
        int non0ctrj[j_prim];
        int non0ctrk[k_prim];
        int non0idxi[i_prim*i_ctr];
        int non0idxj[j_prim*j_ctr];
        int non0idxk[k_prim*k_ctr];

        CINTg4c_index_xyz(idx, envs);
        CINTOpt_non0coeff_byshell(non0idxi, non0ctri, coeff[0], i_prim, i_ctr);
        CINTOpt_non0coeff_byshell(non0idxj, non0ctrj, coeff[1], j_prim, j_ctr);
        CINTOpt_non0coeff_byshell(non0idxk, non0ctrk, coeff[2], k_prim, k_ctr);
        int *non0ctr[3] = {non0ctri, non0ctrj, non0ctrk};
        int *non0idx[3] = {non0idxi, non0idxj, non0idxk};
        double common_factor = envs->common_factor * M_PI*SQRTPI
                * CINTcommon_fac_sp(envs->i_l) * CINTcommon_fac_sp(envs->j_l)
                * CINTcommon_fac_sp(envs->k_l);
        INITSIMD;

        *kempty = 1;
        for (kp = 0; kp < k_prim; kp++) {
                if (k_ctr == 1) {
                        fac1k = common_factor * ck[kp];
                } else {
                        fac1k = common_factor;
                        *jempty = 1;
                }

                for (jp = 0; jp < j_prim; jp++) {
                        if (j_ctr == 1) {
                                fac1j = fac1k * cj[jp];
                        } else {
                                fac1j = fac1k;
                                *iempty = 1;
                        }
                        ajakrr = aj[jp] * ak[kp] * rr_jk;
                        for (ip = 0; ip < i_prim; ip++) {
                                aijk = ai[ip] + aj[jp] + ak[kp];
                                aiakrr = ai[ip] * ak[kp] * rr_ik;
                                aiajrr = ai[ip] * aj[jp] * rr_ij;
                                eijk = (aiajrr+aiakrr+ajakrr) / aijk;
                                if (eijk > CUTOFF15) {
                                        goto i_contracted;
                                }
                                PUSH(eijk);
i_contracted: ;
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
        return !*kempty;
}


int CINT3c1e_loop(double *out, CINTEnvVars *envs, CINTOpt *opt, double *cache)
{
        int *shls  = envs->shls;
        int i_sh = shls[0];
        int j_sh = shls[1];
        int k_sh = shls[2];
        if (opt->data_ptr[i_sh*envs->nbas+j_sh] == NOVALUE ||
            opt->data_ptr[i_sh*envs->nbas+k_sh] == NOVALUE ||
            opt->data_ptr[j_sh*envs->nbas+k_sh] == NOVALUE) {
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
        double *ri = envs->ri;
        double *rk = envs->rk;
        int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        int nf = envs->nf;
        double fac1j, fac1k;
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

        double aijk, eijk;
        double aiajrr, aiakrr, ajakrr;
        double rirk[3];
        rirk[0] = ri[0] - rk[0];
        rirk[1] = ri[1] - rk[1];
        rirk[2] = ri[2] - rk[2];
        const double rr_ij = SQUARE(envs->rirj);
        const double rr_ik = SQUARE(      rirk);
        const double rr_jk = SQUARE(envs->rkrl);
        double common_factor = envs->common_factor * M_PI*SQRTPI
                * CINTcommon_fac_sp(envs->i_l) * CINTcommon_fac_sp(envs->j_l)
                * CINTcommon_fac_sp(envs->k_l);
        int *idx = opt->index_xyz_array[envs->i_l*LMAX1*LMAX1
                                       +envs->j_l*LMAX1
                                       +envs->k_l];
        int *non0ctr[3] = {opt->non0ctr[i_sh], opt->non0ctr[j_sh], opt->non0ctr[k_sh]};
        int *non0idx[3] = {opt->sortedidx[i_sh], opt->sortedidx[j_sh], opt->sortedidx[k_sh]};

        INITSIMD;

        PairData *pdata_ij;
        PairData *pdata_jk = opt->data + opt->data_ptr[j_sh*envs->nbas+k_sh];
        *kempty = 1;
        for (kp = 0; kp < k_prim; kp++) {
                if (k_ctr == 1) {
                        fac1k = common_factor * ck[kp];
                } else {
                        fac1k = common_factor;
                        *jempty = 1;
                }

                for (jp = 0; jp < j_prim; jp++, pdata_jk++) {
                        if (j_ctr == 1) {
                                fac1j = fac1k * cj[jp];
                        } else {
                                fac1j = fac1k;
                                *iempty = 1;
                        }
                        if (pdata_jk->cceij > CUTOFF15) {
                                goto j_contracted;
                        }
                        ajakrr = aj[jp] * ak[kp] * rr_jk;
                        pdata_ij = opt->data + opt->data_ptr[i_sh*envs->nbas+j_sh] + jp*i_prim;
                        for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                if (pdata_ij->cceij > CUTOFF15) {
                                        goto i_contracted;
                                }
                                aijk = ai[ip] + aj[jp] + ak[kp];
                                aiakrr = ai[ip] * ak[kp] * rr_ik;
                                aiajrr = ai[ip] * aj[jp] * rr_ij;
                                eijk = (aiajrr+aiakrr+ajakrr) / aijk;
                                if (eijk > CUTOFF15) {
                                        goto i_contracted;
                                }
                                PUSH(eijk);
i_contracted: ;
                        } // end loop i_prim
                        if (!*iempty) {
                                PRIM2CTR(j, gctr[SHLTYPi]);
                        }
j_contracted: ;
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
        return !*kempty;
}


int CINT3c1e_cart_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
                      double *cache)
{
        int *x_ctr = envs->x_ctr;
        int nc = envs->nf * x_ctr[0] * x_ctr[1] * x_ctr[2];
        int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        if (out == NULL) {
                int leng = envs->g_size*3*((1<<envs->gbits)+1)*SIMDD;
                int len0 = envs->nf*n_comp * SIMDD;
                int cache_size = leng + len0 + nc*n_comp*2 + SIMDD*3;
                return cache_size;
        }
        double *stack = NULL;
        if (cache == NULL) {
                int leng = envs->g_size*3*((1<<envs->gbits)+1)*SIMDD;
                int len0 = envs->nf*n_comp * SIMDD;
                int cache_size = leng + len0 + nc*n_comp*2 + SIMDD*3;
                stack = _mm_malloc(sizeof(double)*cache_size, sizeof(double)*SIMDD);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, nc*n_comp);

        int n, has_value;
        if (opt != NULL) {
                has_value = CINT3c1e_loop(gctr, envs, opt, cache);
        } else {
                has_value = CINT3c1e_loop_nopt(gctr, envs, cache);
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
int CINT3c1e_spheric_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
                         double *cache)
{
        int *x_ctr = envs->x_ctr;
        int nc = envs->nf * x_ctr[0] * x_ctr[1] * x_ctr[2];
        int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        if (out == NULL) {
                int leng = envs->g_size*3*((1<<envs->gbits)+1)*SIMDD;
                int len0 = envs->nf*n_comp * SIMDD;
                int cache_size = MAX(leng+len0+nc*n_comp*2, nc*n_comp+envs->nf*3) + SIMDD*3;
                return cache_size;
        }
        double *stack = NULL;
        if (cache == NULL) {
                int leng = envs->g_size*3*((1<<envs->gbits)+1)*SIMDD;
                int len0 = envs->nf*n_comp * SIMDD;
                int cache_size = MAX(leng+len0+nc*n_comp*2, nc*n_comp+envs->nf*3) + SIMDD*3;
                stack = _mm_malloc(sizeof(double)*cache_size, sizeof(double)*SIMDD);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, nc*n_comp);

        int n, has_value;
        if (opt != NULL) {
                has_value = CINT3c1e_loop(gctr, envs, opt, cache);
        } else {
                has_value = CINT3c1e_loop_nopt(gctr, envs, cache);
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
int CINT3c1e_spinor_drv(double complex *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
                        double *cache, void (*f_e1_c2s)())
{
        fprintf(stderr, "CINT3c1e_spinor_drv not implemented");
        exit(1);
}

void CINTgout3c1e(double *gout, double *g, int *idx, CINTEnvVars *envs, int count)
{
        CINTg3c1e_ovlp(g, envs, count);
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

int int3c1e_sph(double *out, int *dims, int *shls, int *atm, int natm,
                int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout3c1e;
        return CINT3c1e_spheric_drv(out, dims, &envs, opt, cache);
}
void int3c1e_optimizer(CINTOpt **opt, int *atm, int natm,
                       int *bas, int nbas, double *env)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTall_3c1e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}

int int3c1e_cart(double *out, int *dims, int *shls, int *atm, int natm,
                 int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout3c1e;
        return CINT3c1e_cart_drv(out, dims, &envs, opt, cache);
}

int int3c1e_spinor(double complex *out, int *dims, int *shls, int *atm, int natm,
                   int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout3c1e;
        return CINT3c1e_spinor_drv(out, dims, &envs, opt, cache, &c2s_sf_3c2e1);
}


ALL_CINT(int3c1e)
//ALL_CINT_FORTRAN_(cint3c1e)


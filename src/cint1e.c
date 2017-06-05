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
#include "g1e.h"
#include "optimizer.h"
#include "cint1e.h"
#include "cart2sph.h"
#include "c2f.h"

#define SHLTYPi       0
#define SHLTYPj       1


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
                (*envs->f_gout)(gout, g, idx, envs, cum); \
                POP_PRIM2CTR; \
        } \
        envs->ai[cum] = ai[ip]; \
        envs->aj[cum] = aj[jp]; \
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
        int cum = 0; \
        int np2c = 0; \
        double *gprim[SIMDD*2]; \
        int shltyp[SIMDD*2]; \
        int iprim[SIMDD*2]; \
        void (*fp2c[SIMDD*2])(); \
        MM_STORE(envs->ai, MM_SET1(1.)); \
        MM_STORE(envs->aj, MM_SET1(1.)); \
        MM_STORE(envs->fac, MM_SET1(0.));

#define RUN_REST \
        if (cum > 0) { \
                (*envs->f_gout)(gout, g, idx, envs, cum); \
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

int CINT1e_loop_nopt(double *out, CINTEnvVars *envs, double *cache)
{
        int *shls = envs->shls;
        int *bas = envs->bas;
        double *env = envs->env;
        int i_sh = shls[0];
        int j_sh = shls[1];
        int i_ctr = envs->x_ctr[0];
        int j_ctr = envs->x_ctr[1];
        int i_prim = bas(NPRIM_OF, i_sh);
        int j_prim = bas(NPRIM_OF, j_sh);
        int *x_ctr = envs->x_ctr;
        int x_prim[2] = {i_prim, j_prim};
        double *ai = env + bas(PTR_EXP, i_sh);
        double *aj = env + bas(PTR_EXP, j_sh);
        double *ci = env + bas(PTR_COEFF, i_sh);
        double *cj = env + bas(PTR_COEFF, j_sh);
        double *coeff[4] = {ci, cj};
        double *ri = envs->ri;
        double *rj = envs->rj;
        int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        int nf = envs->nf;
        double fac1i, fac1j;
        int ip, jp, i, it, im;
        int empty[2] = {1, 1};
        int *iempty = empty + 0;
        int *jempty = empty + 1;
        int ngp[2];
        ngp[0] = nf * n_comp;
        ngp[1] = ngp[0] * i_ctr;
        int leng = envs->g_size * 3 * ((1<<envs->gbits)+1) * SIMDD;
        int len0 = ngp[0] * SIMDD;
        int leni = ALIGN_UP(ngp[1], SIMDD);
        int lenj = ALIGN_UP(ngp[1] * j_ctr, SIMDD);

        double *gout, *g, *g1;
        double *gctr[2];
        if (n_comp == 1) {
                MALLOC_INSTACK(g1, lenj);
                gctr[SHLTYPj] = out;
        } else {
                MALLOC_INSTACK(gctr[SHLTYPj], lenj);
                g1 = out;
        }
        ALIAS_ADDR_IF_EQUAL(i, j);
        MALLOC_INSTACK(gout, len0*2+leng);
        g = gout + len0*2;

        double rr_ij = SQUARE(envs->rirj);
        double log_rr_ij = (envs->li_ceil+envs->lj_ceil+1)*approx_log(rr_ij+1)/2;
        double aij, eij, cceij;
        double rij[4];
        int idx[nf*3];
        int non0ctri[i_prim];
        int non0ctrj[j_prim];
        int non0idxi[i_prim*i_ctr];
        int non0idxj[j_prim*j_ctr];

        CINTg2c_index_xyz(idx, envs);
        CINTOpt_non0coeff_byshell(non0idxi, non0ctri, coeff[0], i_prim, i_ctr);
        CINTOpt_non0coeff_byshell(non0idxj, non0ctrj, coeff[1], j_prim, j_ctr);
        int *non0ctr[2] = {non0ctri, non0ctrj};
        int *non0idx[2] = {non0idxi, non0idxj};
        double common_factor = envs->common_factor
                * CINTcommon_fac_sp(envs->i_l) * CINTcommon_fac_sp(envs->j_l);
        INITSIMD;

        for (jp = 0; jp < j_prim; jp++) {
                if (j_ctr == 1) {
                        fac1j = common_factor * cj[jp];
                } else {
                        fac1j = common_factor;
                        *iempty = 1;
                }
                for (ip = 0; ip < i_prim; ip++) {
                        aij = 1 / (ai[ip] + aj[jp]);
                        eij = rr_ij * ai[ip] * aj[jp] * aij;
                        cceij = eij - log_rr_ij;
                        if (cceij > CUTOFF15) {
                                goto i_contracted;
                        }
                        for (i = 0; i < 3; i++) {
                                rij[i] = (ai[ip]*ri[i] + aj[jp]*rj[i]) * aij;
                        }
                        PUSH(rij, exp(-eij));
i_contracted: ;
                } // end loop i_prim
                if (!*iempty) {
                        PRIM2CTR(j, gctr[SHLTYPi]);
                }
        } // end loop j_prim
        RUN_REST;

        if (n_comp > 1 && !*jempty) {
                int nc = i_ctr * j_ctr;
                CINTdmat_transpose(out, gctr[SHLTYPj], nf*nc, n_comp);
        }
        return !*jempty;
}


int CINT1e_loop(double *out, CINTEnvVars *envs, CINTOpt *opt, double *cache)
{
        int *shls  = envs->shls;
        int i_sh = shls[0];
        int j_sh = shls[1];
        if (opt->data_ptr[i_sh*envs->nbas+j_sh] == NOVALUE) {
                return 0;
        }
        int *bas = envs->bas;
        double *env = envs->env;
        int i_ctr  = envs->x_ctr[0];
        int j_ctr  = envs->x_ctr[1];
        int i_prim = bas(NPRIM_OF, i_sh);
        int j_prim = bas(NPRIM_OF, j_sh);
        int *x_ctr = envs->x_ctr;
        int x_prim[2] = {i_prim, j_prim};
        double *ai = env + bas(PTR_EXP, i_sh);
        double *aj = env + bas(PTR_EXP, j_sh);
        double *ci = env + bas(PTR_COEFF, i_sh);
        double *cj = env + bas(PTR_COEFF, j_sh);
        double *coeff[2] = {ci, cj};
        int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        int nf = envs->nf;
        double fac1i, fac1j;
        int ip, jp, i, it, im;
        int empty[2] = {1, 1};
        int *iempty = empty + 0;
        int *jempty = empty + 1;
        int ngp[2];
        ngp[0] = nf * n_comp;
        ngp[1] = ngp[0] * i_ctr;
        // (irys,i,j,k,l,coord,0:1); +1 for nabla-r12
        int leng = envs->g_size * 3 * ((1<<envs->gbits)+1) * SIMDD;
        int len0 = ngp[0] * SIMDD;
        int leni = ALIGN_UP(ngp[1], SIMDD);
        int lenj = ALIGN_UP(ngp[1] * j_ctr, SIMDD);

        double *gout, *g, *g1;
        double *gctr[2];
        if (n_comp == 1) {
                MALLOC_INSTACK(g1, lenj);
                gctr[SHLTYPj] = out;
        } else {
                MALLOC_INSTACK(gctr[SHLTYPj], lenj);
                g1 = out;  // Use out as cache for gctrk, gctrj, gctri
        }
        ALIAS_ADDR_IF_EQUAL(i, j);
        MALLOC_INSTACK(gout, len0*2+leng);
        g = gout + len0*2;

        double common_factor = envs->common_factor
                * CINTcommon_fac_sp(envs->i_l) * CINTcommon_fac_sp(envs->j_l);

        int *idx = opt->index_xyz_array[envs->i_l*ANG_MAX+envs->j_l];
        int *non0ctr[2] = {opt->non0ctr[i_sh], opt->non0ctr[j_sh]};
        int *non0idx[2] = {opt->sortedidx[i_sh], opt->sortedidx[j_sh]};

        INITSIMD;

        PairData *pdata_ij = opt->data + opt->data_ptr[i_sh*envs->nbas+j_sh];
        for (jp = 0; jp < j_prim; jp++) {
                if (j_ctr == 1) {
                        fac1j = common_factor * cj[jp];
                } else {
                        fac1j = common_factor;
                        *iempty = 1;
                }
                for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                        if (pdata_ij->cceij > CUTOFF15) {
                                goto i_contracted;
                        }
                        PUSH(pdata_ij->rij, pdata_ij->eij);
i_contracted: ;
                } // end loop i_prim
                if (!*iempty) {
                        PRIM2CTR(j, gctr[SHLTYPi]);
                }
        } // end loop j_prim
        RUN_REST;

        if (n_comp > 1 && !*jempty) {
                int nc = i_ctr * j_ctr;
                CINTdmat_transpose(out, gctr[SHLTYPj], nf*nc, n_comp);
        }
        return !*jempty;
}

int int1e_cache_size(CINTEnvVars *envs)
{
        int *x_ctr = envs->x_ctr;
        int nc = envs->nf * x_ctr[0] * x_ctr[1];
        int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        int leng = envs->g_size*3*((1<<envs->gbits)+1)*SIMDD;
        int len0 = envs->nf*n_comp * SIMDD;
        int cache_size = MAX(leng+len0*2+nc*n_comp*2,
                             nc*n_comp + envs->nf*8*OF_CMPLX) + SIMDD*2;
        return cache_size;
}

int CINT1e_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
               double *cache, void (*f_c2s)())
{
        if (out == NULL) {
                return int1e_cache_size(envs);
        }
        int *x_ctr = envs->x_ctr;
        int nc = envs->nf * x_ctr[0] * x_ctr[1];
        int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        double *stack = NULL;
        if (cache == NULL) {
                int cache_size = int1e_cache_size(envs);
                stack = _mm_malloc(sizeof(double)*cache_size, sizeof(double)*SIMDD);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, nc*n_comp);

        int n, has_value;
        if (opt != NULL) {
                has_value = CINT1e_loop(gctr, envs, opt, cache);
        } else {
                has_value = CINT1e_loop_nopt(gctr, envs, cache);
        }

        int counts[4];
        if (dims == NULL) {
                dims = counts;
        }
        if (f_c2s == c2s_sph_1e) {
                counts[0] = (envs->i_l*2+1) * x_ctr[0];
                counts[1] = (envs->j_l*2+1) * x_ctr[1];
        } else {
                counts[0] = envs->nfi * x_ctr[0];
                counts[1] = envs->nfj * x_ctr[1];
        }
        counts[2] = 1;
        counts[3] = 1;
        int nout = dims[0] * dims[1];
        if (has_value) {
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
        return has_value;
}
int CINT1e_spinor_drv(double complex *out, int *dims, CINTEnvVars *envs,
                      CINTOpt *opt, double *cache, void (*f_c2s)())
{
        if (out == NULL) {
                return int1e_cache_size(envs);
        }
        int *x_ctr = envs->x_ctr;
        int nc = envs->nf * x_ctr[0] * x_ctr[1] * envs->ncomp_e1;
        double *stack = NULL;
        if (cache == NULL) {
                int cache_size = int1e_cache_size(envs);
                stack = _mm_malloc(sizeof(double)*cache_size, sizeof(double)*SIMDD);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, nc*envs->ncomp_tensor);

        int n, has_value;
        if (opt != NULL) {
                has_value = CINT1e_loop(gctr, envs, opt, cache);
        } else {
                has_value = CINT1e_loop_nopt(gctr, envs, cache);
        }

        int counts[4];
        if (dims == NULL) {
                dims = counts;
        }
        counts[0] = CINTcgto_spinor(envs->shls[0], envs->bas);
        counts[1] = CINTcgto_spinor(envs->shls[1], envs->bas);
        counts[2] = 1;
        counts[3] = 1;
        int nout = dims[0] * dims[1];
        if (has_value) {
                for (n = 0; n < envs->ncomp_tensor; n++) {
                        (*f_c2s)(out+nout*n, gctr+nc*n, dims, envs, cache);
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

void CINTgout1e(double *gout, double *g, int *idx, CINTEnvVars *envs, int count)
{
        CINTg1e_ovlp(g, envs, count);
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

void CINTgout1e_nuc(double *gout, double *g, int *idx, CINTEnvVars *envs, int count)
{
        int nf = envs->nf;
        int nfc = nf;
        double *gtmp = gout + nf * SIMDD;
        int nrys_roots = envs->nrys_roots;
        int ia, n, i;
        double *gx, *gy, *gz;
        __MD r0;
        for (n = 0; n < nf*SIMDD; n++) {
                gtmp[n] = 0;
        }
        for (ia = 0; ia < envs->natm; ia++) {
                CINTg1e_nuc(g, envs, count, ia);
                for (n = 0; n < nf; n++) {
                        gx = g + idx[n*3+0] * SIMDD;
                        gy = g + idx[n*3+1] * SIMDD;
                        gz = g + idx[n*3+2] * SIMDD;
                        r0 = MM_LOAD(gtmp+n*SIMDD);
                        for (i = 0; i < nrys_roots; i++) {
                                r0 += MM_LOAD(gx+i*SIMDD) * MM_LOAD(gy+i*SIMDD) * MM_LOAD(gz+i*SIMDD);
                        }
                        MM_STORE(gtmp+n*SIMDD, r0);
                }
        }
        CINTsort_gout(gout, gtmp, nfc, SIMDD);
}

int int1e_ovlp_sph(double *out, int *dims, int *shls, int *atm, int natm,
                   int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e;
        return CINT1e_drv(out, dims, &envs, opt, cache, &c2s_sph_1e);
}
void int1e_ovlp_optimizer(CINTOpt **opt, int *atm, int natm,
                          int *bas, int nbas, double *env)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTall_1e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}

int int1e_ovlp_cart(double *out, int *dims, int *shls, int *atm, int natm,
                    int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e;
        return CINT1e_drv(out, dims, &envs, opt, cache, &c2s_cart_1e);
}

int int1e_ovlp_spinor(double complex *out, int *dims, int *shls, int *atm, int natm,
                      int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e;
        return CINT1e_spinor_drv(out, dims, &envs, opt, cache, &c2s_sf_1e);
}

int int1e_nuc_sph(double *out, int *dims, int *shls, int *atm, int natm,
                  int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 0, 1};
        CINTEnvVars envs;
        CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_nuc;
        return CINT1e_drv(out, dims, &envs, opt, cache, &c2s_sph_1e);
}
void int1e_nuc_optimizer(CINTOpt **opt, int *atm, int natm,
                         int *bas, int nbas, double *env)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 0, 1};
        CINTall_1e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
int int1e_nuc_cart(double *out, int *dims, int *shls, int *atm, int natm,
                   int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 0, 1};
        CINTEnvVars envs;
        CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_nuc;
        return CINT1e_drv(out, dims, &envs, opt, cache, &c2s_cart_1e);
}

int int1e_nuc_spinor(double complex *out, int *dims, int *shls, int *atm, int natm,
                     int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 0, 1};
        CINTEnvVars envs;
        CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_nuc;
        return CINT1e_spinor_drv(out, dims, &envs, opt, cache, &c2s_sf_1e);
}


ALL_CINT1E(int1e_ovlp)
ALL_CINT1E(int1e_nuc)
ALL_CINT1E_FORTRAN_(int1e_ovlp)
ALL_CINT1E_FORTRAN_(int1e_nuc)


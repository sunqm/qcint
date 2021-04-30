/*
 * Copyright (C) 2021  Qiming Sun <osirpt.sun@gmail.com>
 *
 * <i|1/r|j> integrals for multiple grids
 */

#include <stdlib.h>
#include <math.h>
#include "cint_bas.h"
#include "optimizer.h"
#include "g1e.h"
#include "g1e_grids.h"
#include "misc.h"
#include "cart2sph.h"
#include "c2f.h"
#include "simd.h"

#define PRIM2CTR(ctrsymb, gp, ngp) \
        if (ctrsymb##_ctr > 1) {\
                if (*ctrsymb##empty) { \
                        CINTprim_to_ctr_0(gctr##ctrsymb, gp, c##ctrsymb+ctrsymb##p, \
                                          ngp, ctrsymb##_prim, ctrsymb##_ctr, \
                                          non0ctr##ctrsymb[ctrsymb##p], \
                                          non0idx##ctrsymb+ctrsymb##p*ctrsymb##_ctr); \
                } else { \
                        CINTprim_to_ctr_1(gctr##ctrsymb, gp, c##ctrsymb+ctrsymb##p, \
                                          ngp, ctrsymb##_prim, ctrsymb##_ctr, \
                                          non0ctr##ctrsymb[ctrsymb##p], \
                                          non0idx##ctrsymb+ctrsymb##p*ctrsymb##_ctr); \
                } \
        } \
        *ctrsymb##empty = 0

static void _transpose_comps(double *gctr, double *gctrj,
                             int bgrids, int dij, int ngrids, int n_comp);

int CINT1e_grids_loop(double *gctr, CINTEnvVars *envs, double *cache)
{
        int *shls  = envs->shls;
        int *bas = envs->bas;
        double *env = envs->env;
        int i_sh = shls[0];
        int j_sh = shls[1];
        int i_ctr = envs->x_ctr[0];
        int j_ctr = envs->x_ctr[1];
        int i_prim = bas(NPRIM_OF, i_sh);
        int j_prim = bas(NPRIM_OF, j_sh);
        int nf = envs->nf;
        int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        int ngrids = ALIGN_UP(envs->ngrids, SIMDD);
        double *grids = envs->grids;
        double *ai = env + bas(PTR_EXP, i_sh);
        double *aj = env + bas(PTR_EXP, j_sh);
        double *ci = env + bas(PTR_COEFF, i_sh);
        double *cj = env + bas(PTR_COEFF, j_sh);

        double expcutoff = envs->expcutoff;
        double *log_maxci, *log_maxcj;
        PairData *pdata_base, *pdata_ij;
        MALLOC_INSTACK(log_maxci, i_prim+j_prim);
        MALLOC_INSTACK(pdata_base, i_prim*j_prim);
        log_maxcj = log_maxci + i_prim;
        CINTOpt_log_max_pgto_coeff(log_maxci, ci, i_prim, i_ctr);
        CINTOpt_log_max_pgto_coeff(log_maxcj, cj, j_prim, j_ctr);
        if (CINTset_pairdata(pdata_base, ai, aj, envs->ri, envs->rj,
                             log_maxci, log_maxcj, envs->li_ceil, envs->lj_ceil,
                             i_prim, j_prim, SQUARE(envs->rirj), expcutoff)) {
                return 0;
        }

        double fac1i, fac1j, expij;
        double *rij;
        int ip, jp, i, grids_offset, bgrids;
        int empty[4] = {1, 1, 1, 1};
        int *gempty = empty + 0;
        int *iempty = empty + 1;
        int *jempty = empty + 2;
        int all_empty = 1;

        int *idx;
        MALLOC_DATA_INSTACK(idx, nf * 3);
        CINTg2c_index_xyz(idx, envs);

        int *non0ctri, *non0ctrj;
        int *non0idxi, *non0idxj;
        MALLOC_INSTACK(non0ctri, i_prim+j_prim+i_prim*i_ctr+j_prim*j_ctr);
        non0ctrj = non0ctri + i_prim;
        non0idxi = non0ctrj + j_prim;
        non0idxj = non0idxi + i_prim*i_ctr;
        CINTOpt_non0coeff_byshell(non0idxi, non0ctri, ci, i_prim, i_ctr);
        CINTOpt_non0coeff_byshell(non0idxj, non0ctrj, cj, j_prim, j_ctr);

        const int nc = i_ctr * j_ctr;
        const int leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const int lenj = GRID_BLKSIZE * nf * nc * n_comp; // gctrj
        const int leni = GRID_BLKSIZE * nf * i_ctr * n_comp; // gctri
        const int len0 = GRID_BLKSIZE * nf * n_comp; // gout
        const int len = leng + lenj + leni + len0;
        double *gridsT;
        MALLOC_ALIGNED_DOUBLE_INSTACK(gridsT, len + GRID_BLKSIZE * 3);
        double *g = gridsT + GRID_BLKSIZE * 3;
        double *g1 = g + leng;
        double *gout, *gctri, *gctrj;
        if (n_comp == 1) {
                gctrj = gctr;
        } else {
                gctrj = g1;
                g1 += lenj;
        }
        if (j_ctr == 1) {
                gctri = gctrj;
                iempty = jempty;
        } else {
                gctri = g1;
                g1 += leni;
        }
        if (i_ctr == 1) {
                gout = gctri;
                gempty = iempty;
        } else {
                gout = g1;
        }

#if (SIMDD == 8)
        __m256i vindex = _mm256_set_epi32(21, 18, 15, 12, 9, 6, 3, 0);
#elif __AVX2__
        __m128i vindex = _mm_set_epi32(9, 6, 3, 0);
#endif
        for (grids_offset = 0; grids_offset < ngrids; grids_offset += GRID_BLKSIZE) {
                envs->grids_offset = grids_offset;
                bgrids = MIN(ngrids - grids_offset, GRID_BLKSIZE);
                bgrids = ALIGN_UP(bgrids, SIMDD);
#if (SIMDD == 8) || __AVX2__
                for (i = 0; i < bgrids; i += SIMDD) {
                        MM_STORE(gridsT+i+GRID_BLKSIZE*0, MM_GATHER(grids+(grids_offset+i)*3+0, vindex, 8));
                        MM_STORE(gridsT+i+GRID_BLKSIZE*1, MM_GATHER(grids+(grids_offset+i)*3+1, vindex, 8));
                        MM_STORE(gridsT+i+GRID_BLKSIZE*2, MM_GATHER(grids+(grids_offset+i)*3+2, vindex, 8));
                }
#else
                for (i = 0; i < bgrids; i++) {
                        gridsT[i+GRID_BLKSIZE*0] = grids[(grids_offset+i)*3+0];
                        gridsT[i+GRID_BLKSIZE*1] = grids[(grids_offset+i)*3+1];
                        gridsT[i+GRID_BLKSIZE*2] = grids[(grids_offset+i)*3+2];
                }
#endif

                empty[0] = 1;
                empty[1] = 1;
                empty[2] = 1;
                if (n_comp == 1) {
                        gctrj = gctr + grids_offset * nf * nc;
                        if (j_ctr == 1) {
                                gctri = gctrj;
                        }
                        if (i_ctr == 1) {
                                gout = gctri;
                        }
                }
                pdata_ij = pdata_base;
                for (jp = 0; jp < j_prim; jp++) {
                        envs->aj[0] = aj[jp];
                        if (j_ctr == 1) {
                                fac1j = envs->common_factor * cj[jp];
                        } else {
                                fac1j = envs->common_factor;
                                *iempty = 1;
                        }
                        for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                if (pdata_ij->cceij > expcutoff) {
                                        continue;
                                }
                                envs->ai[0] = ai[ip];
                                expij = pdata_ij->eij;
                                rij = pdata_ij->rij;
                                envs->rij[0] = rij[0];
                                envs->rij[1] = rij[1];
                                envs->rij[2] = rij[2];
                                if (i_ctr == 1) {
                                        fac1i = fac1j*ci[ip]*expij;
                                } else {
                                        fac1i = fac1j*expij;
                                }

                                CINTg0_1e_grids(g, fac1i, envs, cache, gridsT);
                                (*envs->f_gout)(gout, g, idx, envs, *gempty);
                                PRIM2CTR(i, gout, bgrids * nf * n_comp);
                        }
                        if (!*iempty) {
                                PRIM2CTR(j, gctri, bgrids * nf * i_ctr * n_comp);
                        }
                }
                if (n_comp > 1 && !*jempty) {
                        _transpose_comps(gctr+grids_offset*nf*nc, gctrj,
                                         bgrids, nf*nc, ngrids, n_comp);
                }
                all_empty &= *jempty;
        }
        return !all_empty;
}

static void _transpose_comps(double *RESTRICT gctr, double *RESTRICT gctrj,
                             int bgrids, int dij, int ngrids, int n_comp)
{
        int n, ic, ig;
        double *pgctr, *pgctrj;
        for (ic = 0; ic < n_comp; ic++) {
                pgctr = gctr + ic * dij * ngrids;
                for (n = 0; n < dij; n++) {
                        pgctrj = gctrj + (n * n_comp + ic) * bgrids;
                        for (ig = 0; ig < bgrids; ig++) {
                                pgctr[ig + n * bgrids] = pgctrj[ig];
                        }
                }
        } 
}

size_t int1e_grids_cache_size(CINTEnvVars *envs)
{
        int *bas = envs->bas;
        int *shls  = envs->shls;
        int *x_ctr = envs->x_ctr;
        int ngrids = ALIGN_UP(envs->ngrids, SIMDD);
        int nroots = envs->nrys_roots;
        int nf = envs->nf;
        int nc = ngrids * nf * x_ctr[0] * x_ctr[1];
        int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        int i_prim = bas(NPRIM_OF, shls[0]);
        int j_prim = bas(NPRIM_OF, shls[1]);
        int pdata_size = (i_prim*j_prim * 5
                           + i_prim * x_ctr[0]
                           + j_prim * x_ctr[1]
                           +(i_prim+j_prim)*2 + envs->nf*3);
        size_t leng = envs->g_size*3*((1<<envs->gbits)+1);
        size_t len0 = GRID_BLKSIZE * nf * n_comp;
        size_t leni = len0 * x_ctr[0];
        size_t lenj = leni * x_ctr[1];
        size_t cache_size = MAX(nc*n_comp + leng + len0 + leni + lenj + pdata_size +
                                GRID_BLKSIZE*MAX(n_comp, nroots+6),
                                nc*n_comp + GRID_BLKSIZE * nf*8*OF_CMPLX);
        return cache_size + SIMDD*8;
}

/*
 * 1e integrals <i|O|j> without 1/r
 */
CACHE_SIZE_T CINT1e_grids_drv(double *out, int *dims, CINTEnvVars *envs,
                              double *cache, void (*f_c2s)())
{
        if (out == NULL) {
                return int1e_grids_cache_size(envs);
        }
        int *x_ctr = envs->x_ctr;
        int ngrids_nf = ALIGN_UP(envs->ngrids, SIMDD) * envs->nf;
        int nc = ngrids_nf * x_ctr[0] * x_ctr[1];
        int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        double *stack = NULL;
        if (cache == NULL) {
                size_t cache_size = int1e_grids_cache_size(envs);
                stack = malloc(sizeof(double)*cache_size);
                cache = stack;
        }
        double *gctr;
        MALLOC_ALIGNED_DOUBLE_INSTACK(gctr, nc * n_comp);

        int has_value = CINT1e_grids_loop(gctr, envs, cache);

        int counts[4];
        if (dims == NULL) {
                dims = counts;
        }
        if (f_c2s == &c2s_sph_1e_grids) {
                counts[0] = (envs->i_l*2+1) * x_ctr[0];
                counts[1] = (envs->j_l*2+1) * x_ctr[1];
                counts[2] = envs->ngrids;
                counts[3] = 1;
        } else if (f_c2s == &c2s_cart_1e_grids) {
                counts[0] = envs->nfi * x_ctr[0];
                counts[1] = envs->nfj * x_ctr[1];
                counts[2] = envs->ngrids;
                counts[3] = 1;
        }
        int nout = dims[0] * dims[1] * dims[2];
        int n;
        if (has_value) {
                for (n = 0; n < n_comp; n++) {
                        (*f_c2s)(out+nout*n, gctr+nc*n, dims, envs, cache);
                }
        } else {
                for (n = 0; n < n_comp; n++) {
                        c2s_grids_dset0(out+nout*n, dims, counts);
                }
        }
        if (stack != NULL) {
                free(stack);
        }
        return has_value;
}

CACHE_SIZE_T CINT1e_grids_spinor_drv(double complex *out, int *dims, CINTEnvVars *envs,
                                     double *cache, void (*f_c2s)())
{
        if (out == NULL) {
                return int1e_grids_cache_size(envs);
        }
        int *x_ctr = envs->x_ctr;
        int ngrids_nf = ALIGN_UP(envs->ngrids, SIMDD) * envs->nf;
        int nc = ngrids_nf * x_ctr[0] * x_ctr[1] * envs->ncomp_e1;
        double *stack = NULL;
        if (cache == NULL) {
                size_t cache_size = int1e_grids_cache_size(envs);
                stack = malloc(sizeof(double)*cache_size);
                cache = stack;
        }
        double *gctr;
        MALLOC_ALIGNED_DOUBLE_INSTACK(gctr, nc * envs->ncomp_tensor);

        int has_value = CINT1e_grids_loop(gctr, envs, cache);

        int counts[4];
        if (dims == NULL) {
                dims = counts;
        }
        counts[0] = CINTcgto_spinor(envs->shls[0], envs->bas);
        counts[1] = CINTcgto_spinor(envs->shls[1], envs->bas);
        counts[2] = envs->ngrids;
        counts[3] = 1;
        int nout = dims[0] * dims[1] * dims[2];
        int n;
        if (has_value) {
                for (n = 0; n < envs->ncomp_tensor; n++) {
                        (*f_c2s)(out+nout*n, gctr+nc*n, dims, envs, cache);
                }
        } else {
                for (n = 0; n < envs->ncomp_tensor; n++) {
                        c2s_grids_zset0(out+nout*n, dims, counts);
                }
        }
        if (stack != NULL) {
                free(stack);
        }
        return has_value;
}

CACHE_SIZE_T int1e_grids_sph(double *out, int *dims, int *shls, int *atm, int natm,
                     int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 0, 1};
        CINTEnvVars envs;
        CINTinit_int1e_grids_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_grids;
        return CINT1e_grids_drv(out, dims, &envs, cache, &c2s_sph_1e_grids);
}

void int1e_grids_optimizer(CINTOpt **opt, int *atm, int natm,
                           int *bas, int nbas, double *env)
{
        *opt = NULL;
}


CACHE_SIZE_T int1e_grids_cart(double *out, int *dims, int *shls, int *atm, int natm,
                              int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 0, 1};
        CINTEnvVars envs;
        CINTinit_int1e_grids_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_grids;
        return CINT1e_grids_drv(out, dims, &envs, cache, &c2s_cart_1e_grids);
}

CACHE_SIZE_T int1e_grids_spinor(double complex *out, int *dims, int *shls, int *atm, int natm,
                                int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 0, 1};
        CINTEnvVars envs;
        CINTinit_int1e_grids_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_grids;
        return CINT1e_grids_spinor_drv(out, dims, &envs, cache, &c2s_sf_1e_grids);
}

ALL_CINT(int1e_grids)


/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 */

#include <stdlib.h>
#include <stdio.h>
#include "cint_bas.h"
#include "g2e.h"
#include "optimizer.h"
#include "cint2e.h"
#include "misc.h"
#include "cart2sph.h"


#ifdef WITH_GTG
int CINTg0_2e_gtg(double *g, Rys2eT *bc, CINTEnvVars *envs, int count);
int CINTg0_2e_gtg_simd1(double *g, Rys2eT *bc, CINTEnvVars *envs, int idsimd);

void CINTinit_int2c2e_gtg_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
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
        envs->common_factor = SQRTPI * .5;
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

        envs->ri = env + atm(PTR_COORD, bas(ATOM_OF, i_sh));
        envs->rk = env + atm(PTR_COORD, bas(ATOM_OF, k_sh));

        envs->nrys_roots = 1;

        int dli = envs->li_ceil + 1;
        int dlk = envs->lk_ceil + 1;
        envs->g_stride_i = 1;
        envs->g_stride_k = dli;
        envs->g_stride_l = envs->g_stride_k;
        envs->g_size     = dli * dlk;

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

        envs->f_g0_2d4d = &CINTg0_2e_2d;
        envs->f_g0_2d4d_simd1 = &CINTg0_2e_2d_simd1;
        envs->f_g0_2e = &CINTg0_2e_gtg;
        envs->f_g0_2e_simd1 = &CINTg0_2e_gtg_simd1;

        // for CINTg2c_index_xyz and c2s_sph_1e function
        envs->j_l = envs->k_l;
        envs->nfj = envs->nfk;
        envs->g_stride_j = envs->g_stride_k;
}
#endif

CACHE_SIZE_T int2c2e_gtg_sph(double *out, int *dims, int *shls, int *atm, int natm,
                  int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2c2e_gtg_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        envs.f_gout_simd1 = &CINTgout2e_simd1;
        return CINT2c2e_drv(out, dims, &envs, opt, cache, &c2s_sph_1e);
}
void int2c2e_gtg_optimizer(CINTOpt **opt, int *atm, int natm,
                         int *bas, int nbas, double *env)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTall_2c2e_gtg_optimizer(opt, ng, atm, natm, bas, nbas, env);
}

CACHE_SIZE_T int2c2e_gtg_cart(double *out, int *dims, int *shls, int *atm, int natm,
                     int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2c2e_gtg_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        envs.f_gout_simd1 = &CINTgout2e_simd1;
        return CINT2c2e_drv(out, dims, &envs, opt, cache, &c2s_cart_1e);
}

CACHE_SIZE_T int2c2e_gtg_spinor(double complex *out, int *dims, int *shls, int *atm, int natm,
                     int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        fprintf(stderr, "int2c2e_gtg_spinor not implemented\n");
}

ALL_CINT(int2c2e_gtg)


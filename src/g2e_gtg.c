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
#include <math.h>
#include <assert.h>
#include "cint_const.h"
#include "cint_bas.h"
#include "simd.h"
#include "g2e.h"
#include "rys_roots.h"

void CINTg0_2e_lj2d4d_regular(double *g, Rys2eT *bc, CINTEnvVars *envs);
void CINTg0_2e_lj2d4d_simd1_regular(double *g, Rys2eT *bc, CINTEnvVars *envs);
void CINTg0_2e_gtg(double *g, Rys2eT *bc, CINTEnvVars *envs, int count);
void CINTg0_2e_gtg_simd1(double *g, Rys2eT *bc, CINTEnvVars *envs, int idsimd);
void CINTg0_lj_4d(double *g, CINTEnvVars *envs);
void CINTg0_lj_4d_simd1(double *g, CINTEnvVars *envs);

void CINTinit_int2e_gtg_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                                int *atm, int natm, int *bas, int nbas, double *env)
{
        envs->natm = natm;
        envs->nbas = nbas;
        envs->atm = atm;
        envs->bas = bas;
        envs->env = env;
        envs->shls = shls;

        const int i_sh = shls[0];
        const int j_sh = shls[1];
        const int k_sh = shls[2];
        const int l_sh = shls[3];
        envs->i_l = bas(ANG_OF, i_sh);
        envs->j_l = bas(ANG_OF, j_sh);
        envs->k_l = bas(ANG_OF, k_sh);
        envs->l_l = bas(ANG_OF, l_sh);
        envs->x_ctr[0] = bas(NCTR_OF, i_sh);
        envs->x_ctr[1] = bas(NCTR_OF, j_sh);
        envs->x_ctr[2] = bas(NCTR_OF, k_sh);
        envs->x_ctr[3] = bas(NCTR_OF, l_sh);
        envs->nfi = (envs->i_l+1)*(envs->i_l+2)/2;
        envs->nfj = (envs->j_l+1)*(envs->j_l+2)/2;
        envs->nfk = (envs->k_l+1)*(envs->k_l+2)/2;
        envs->nfl = (envs->l_l+1)*(envs->l_l+2)/2;
        envs->nf = envs->nfi * envs->nfk * envs->nfl * envs->nfj;

        envs->ri = env + atm(PTR_COORD, bas(ATOM_OF, i_sh));
        envs->rj = env + atm(PTR_COORD, bas(ATOM_OF, j_sh));
        envs->rk = env + atm(PTR_COORD, bas(ATOM_OF, k_sh));
        envs->rl = env + atm(PTR_COORD, bas(ATOM_OF, l_sh));

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
        envs->lj_ceil = envs->j_l + ng[JINC];
        envs->lk_ceil = envs->k_l + ng[KINC];
        envs->ll_ceil = envs->l_l + ng[LINC];
        envs->nrys_roots = 1;

        assert(i_sh < SHLS_MAX);
        assert(j_sh < SHLS_MAX);
        assert(k_sh < SHLS_MAX);
        assert(l_sh < SHLS_MAX);
        assert(envs->i_l < ANG_MAX);
        assert(envs->j_l < ANG_MAX);
        assert(envs->k_l < ANG_MAX);
        assert(envs->l_l < ANG_MAX);
        assert(bas(ATOM_OF,i_sh) >= 0);
        assert(bas(ATOM_OF,j_sh) >= 0);
        assert(bas(ATOM_OF,k_sh) >= 0);
        assert(bas(ATOM_OF,l_sh) >= 0);
        assert(bas(ATOM_OF,i_sh) < natm);
        assert(bas(ATOM_OF,j_sh) < natm);
        assert(bas(ATOM_OF,k_sh) < natm);
        assert(bas(ATOM_OF,l_sh) < natm);
        assert(envs->nrys_roots < MXRYSROOTS);

        int dli, dlj, dlk, dll;
        int ibase = envs->li_ceil > envs->lj_ceil;
        int kbase = envs->lk_ceil > envs->ll_ceil;

        if (kbase) {
                dlk = envs->lk_ceil + envs->ll_ceil + 1;
                dll = envs->ll_ceil + 1;
        } else {
                dlk = envs->lk_ceil + 1;
                dll = envs->lk_ceil + envs->ll_ceil + 1;
        }

        if (ibase) {
                dli = envs->li_ceil + envs->lj_ceil + 1;
                dlj = envs->lj_ceil + 1;
        } else {
                dli = envs->li_ceil + 1;
                dlj = envs->li_ceil + envs->lj_ceil + 1;
        }
        envs->g_stride_i = envs->nrys_roots;
        envs->g_stride_k = envs->nrys_roots * dli;
        envs->g_stride_l = envs->nrys_roots * dli * dlk;
        envs->g_stride_j = envs->nrys_roots * dli * dlk * dll;
        envs->g_size     = envs->nrys_roots * dli * dlk * dll * dlj;

        if (kbase) {
                envs->g2d_klmax = envs->g_stride_k;
                envs->rx_in_rklrx = envs->rk;
                envs->rkrl[0] = envs->rk[0] - envs->rl[0];
                envs->rkrl[1] = envs->rk[1] - envs->rl[1];
                envs->rkrl[2] = envs->rk[2] - envs->rl[2];
        } else {
                envs->g2d_klmax = envs->g_stride_l;
                envs->rx_in_rklrx = envs->rl;
                envs->rkrl[0] = envs->rl[0] - envs->rk[0];
                envs->rkrl[1] = envs->rl[1] - envs->rk[1];
                envs->rkrl[2] = envs->rl[2] - envs->rk[2];
        }

        if (ibase) {
                envs->g2d_ijmax = envs->g_stride_i;
                envs->rx_in_rijrx = envs->ri;
                envs->rirj[0] = envs->ri[0] - envs->rj[0];
                envs->rirj[1] = envs->ri[1] - envs->rj[1];
                envs->rirj[2] = envs->ri[2] - envs->rj[2];
        } else {
                envs->g2d_ijmax = envs->g_stride_j;
                envs->rx_in_rijrx = envs->rj;
                envs->rirj[0] = envs->rj[0] - envs->ri[0];
                envs->rirj[1] = envs->rj[1] - envs->ri[1];
                envs->rirj[2] = envs->rj[2] - envs->ri[2];
        }

        if (kbase) {
                if (ibase) {
                        envs->f_g0_2d4d = &CINTg0_2e_ik2d4d;
                        envs->f_g0_2d4d_simd1 = &CINTg0_2e_ik2d4d_simd1;
                } else {
                        envs->f_g0_2d4d = &CINTg0_2e_kj2d4d;
                        envs->f_g0_2d4d_simd1 = &CINTg0_2e_kj2d4d_simd1;
                }
        } else {
                if (ibase) {
                        envs->f_g0_2d4d = &CINTg0_2e_il2d4d;
                        envs->f_g0_2d4d_simd1 = &CINTg0_2e_il2d4d_simd1;
                } else {
                        envs->f_g0_2d4d = &CINTg0_2e_lj2d4d_regular;
                        envs->f_g0_2d4d_simd1 = &CINTg0_2e_lj2d4d_simd1_regular;
                }
        }
        envs->f_g0_2e = &CINTg0_2e_gtg;
        envs->f_g0_2e_simd1 = &CINTg0_2e_gtg_simd1;
}


void CINTg0_2e_lj2d4d_regular(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        CINTg0_2e_2d(g, bc, envs);
        CINTg0_lj_4d(g, envs);
        return;
}

void CINTg0_2e_lj2d4d_simd1_regular(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        CINTg0_2e_2d_simd1(g, bc, envs);
        CINTg0_lj_4d_simd1(g, envs);
        return;
}

/*
 * g[i,k,l,j] = < ik | lj > = ( i j | k l )
 */
void CINTg0_2e_gtg(double *g, Rys2eT *bc, CINTEnvVars *envs, int count)
{
        const double zeta = envs->env[PTR_GTG_ZETA];
        double *gx = g;
        double *gy = gx + envs->g_size * SIMDD;
        double *gz = gy + envs->g_size * SIMDD;
        ALIGNMM double aij[SIMDD];
        ALIGNMM double akl[SIMDD];
        ALIGNMM double a0[SIMDD];
        ALIGNMM double a1[SIMDD];
        ALIGNMM double aijkl[SIMDD];
        ALIGNMM double fac1[SIMDD];
        ALIGNMM double x[SIMDD];
        ALIGNMM double t[SIMDD];
        ALIGNMM double rijrkl[SIMDD*3];
        ALIGNMM double rijrx[SIMDD*3];
        ALIGNMM double rklrx[SIMDD*3];
        double *rij = envs->rij;
        double *rkl = envs->rkl;
        __MD ra, r0, r1, r2, r3, r4, r5, r6, r7, r8;

        //:for (int k = 0; k < count; k++) {
        //:        aij[k] = envs->ai[k] + envs->aj[k];
        //:        akl[k] = envs->ak[k] + envs->al[k];
        //:        aijkl[k] = aij[k] + akl[k];
        //:        a1[k] = aij[k] * akl[k];
        //:        a0[k] = a1[k] / aijkl[k];
        //:}
        r2 = MM_ADD(MM_LOAD(envs->ai),  MM_LOAD(envs->aj));
        r3 = MM_ADD(MM_LOAD(envs->ak),  MM_LOAD(envs->al));
        MM_STORE(aij, r2);
        MM_STORE(akl, r3);
        r1 = MM_MUL(r2, r3);
        MM_STORE(a1, r1);
        ra = MM_ADD(r2, r3);
        MM_STORE(aijkl, ra);
        r0 = MM_DIV(r1, ra);
        MM_STORE(a0, r0);

        // t = zeta / (zeta + a0);
        // fac1 = (1-t) / a1;
        r2 = MM_SET1(zeta);
        r2 = MM_DIV(r2, MM_ADD(r2, r0));
        MM_STORE(t, r2);
        MM_STORE(fac1, MM_DIV(MM_SUB(MM_SET1(1), r2), r1));

        r0 = MM_LOAD(rij+0*SIMDD);
        r1 = MM_LOAD(rij+1*SIMDD);
        r2 = MM_LOAD(rij+2*SIMDD);
        r3 = MM_LOAD(rkl+0*SIMDD);
        r4 = MM_LOAD(rkl+1*SIMDD);
        r5 = MM_LOAD(rkl+2*SIMDD);

        //:for (k = 0; k < count; k++) {
        //:        rijrkl[0*SIMDD+k] = rij[0*SIMDD+k] - rkl[0*SIMDD+k];
        //:        rijrkl[1*SIMDD+k] = rij[1*SIMDD+k] - rkl[1*SIMDD+k];
        //:        rijrkl[2*SIMDD+k] = rij[2*SIMDD+k] - rkl[2*SIMDD+k];
        //:        rijrx[0*SIMDD+k] = rij[0*SIMDD+k] - envs->rx_in_rijrx[0];
        //:        rijrx[1*SIMDD+k] = rij[1*SIMDD+k] - envs->rx_in_rijrx[1];
        //:        rijrx[2*SIMDD+k] = rij[2*SIMDD+k] - envs->rx_in_rijrx[2];
        //:        rklrx[0*SIMDD+k] = rkl[0*SIMDD+k] - envs->rx_in_rklrx[0];
        //:        rklrx[1*SIMDD+k] = rkl[1*SIMDD+k] - envs->rx_in_rklrx[1];
        //:        rklrx[2*SIMDD+k] = rkl[2*SIMDD+k] - envs->rx_in_rklrx[2];
        //:        x[k] = a0[k] *(rijrkl[0*SIMDD+k] * rijrkl[0*SIMDD+k]
        //:                     + rijrkl[1*SIMDD+k] * rijrkl[1*SIMDD+k]
        //:                     + rijrkl[2*SIMDD+k] * rijrkl[2*SIMDD+k]);
        //:}
        MM_STORE(rijrx+0*SIMDD, MM_SUB(r0, MM_SET1(envs->rx_in_rijrx[0])));
        MM_STORE(rijrx+1*SIMDD, MM_SUB(r1, MM_SET1(envs->rx_in_rijrx[1])));
        MM_STORE(rijrx+2*SIMDD, MM_SUB(r2, MM_SET1(envs->rx_in_rijrx[2])));
        MM_STORE(rklrx+0*SIMDD, MM_SUB(r3, MM_SET1(envs->rx_in_rklrx[0])));
        MM_STORE(rklrx+1*SIMDD, MM_SUB(r4, MM_SET1(envs->rx_in_rklrx[1])));
        MM_STORE(rklrx+2*SIMDD, MM_SUB(r5, MM_SET1(envs->rx_in_rklrx[2])));

        r6 = MM_SUB(r0, r3); MM_STORE(rijrkl+0*SIMDD, r6);
        r7 = MM_SUB(r1, r4); MM_STORE(rijrkl+1*SIMDD, r7);
        r8 = MM_SUB(r2, r5); MM_STORE(rijrkl+2*SIMDD, r8);
        ra = MM_FMA(r6, r6, MM_FMA(r7, r7, MM_MUL(r8, r8)));

        // x = -x * t
        // gz[0] = fac1*sqrt(fac1) * exp(x) * fac;
        ra = MM_MUL(MM_LOAD(a0), ra);
        ra = MM_MUL(MM_LOAD(t), ra);
        MM_STORE(x, -ra);
        CINTexp_cephes(t, x);
        ra = MM_LOAD(fac1);
        r1 = MM_MUL(MM_SQRT(ra), ra);
        r2 = MM_MUL(MM_LOAD(t), MM_LOAD(envs->fac));
        MM_STORE(gz, MM_MUL(r1, r2));
        r1 = MM_SET1(1.);
        MM_STORE(gx, r1);
        MM_STORE(gy, r1);
        if (envs->g_size == 1) {
                return 1;
        }

        double *b00 = bc->b00;
        double *b10 = bc->b10;
        double *b01 = bc->b01;
        double *c00x = bc->c00x;
        double *c00y = bc->c00y;
        double *c00z = bc->c00z;
        double *c0px = bc->c0px;
        double *c0py = bc->c0py;
        double *c0pz = bc->c0pz;

        ALIGNMM double tmp1[SIMDD];
        ALIGNMM double tmp4[SIMDD];
        //:double u2[SIMDD];
        //:double tmp2, tmp3;
        //:double div[SIMDD];
        //:for (k = 0; k < count; k++) {
        //:        div[k] = 1 / (zeta * aijkl[k] + a1[k]);
        //:        tmp1[k] = zeta * div[k];
        //:        tmp4[k] = .5 * div[k];
        //:
        //:        b00[k] = 0.5 * tmp1[k];
        //:        tmp2 = tmp1[k] * akl[k];
        //:        tmp3 = tmp1[k] * aij[k];
        //:        b10[k] = b00[k] + tmp4[k] * akl[k];
        //:        b01[k] = b00[k] + tmp4[k] * aij[k];
        //:        c00x[k] = rijrx[0*SIMDD+k] - tmp2 * rijrkl[0*SIMDD+k];
        //:        c00y[k] = rijrx[1*SIMDD+k] - tmp2 * rijrkl[1*SIMDD+k];
        //:        c00z[k] = rijrx[2*SIMDD+k] - tmp2 * rijrkl[2*SIMDD+k];
        //:        c0px[k] = rklrx[0*SIMDD+k] + tmp3 * rijrkl[0*SIMDD+k];
        //:        c0py[k] = rklrx[1*SIMDD+k] + tmp3 * rijrkl[1*SIMDD+k];
        //:        c0pz[k] = rklrx[2*SIMDD+k] + tmp3 * rijrkl[2*SIMDD+k];
        //:}

        ra = MM_LOAD(aijkl);
        r0 = MM_SET1(zeta);
        r1 = MM_LOAD(a1);
        r2 = MM_SET1(.5);
        r3 = MM_SET1(1.);

        r5 = MM_DIV(r3, MM_FMA(r0, ra, r1));
        r1 = MM_MUL(r0, r5);
        r4 = MM_MUL(r2, r5);
        MM_STORE(tmp1, r1);
        MM_STORE(tmp4, r4);

        r0 = MM_MUL(r2, r1);
        MM_STORE(b00, r0);
        MM_STORE(b10, MM_FMA(r4, MM_LOAD(akl), r0));
        MM_STORE(b01, MM_FMA(r4, MM_LOAD(aij), r0));

        r4 = MM_LOAD(rijrkl+0*SIMDD);
        r5 = MM_LOAD(rijrkl+1*SIMDD);
        r6 = MM_LOAD(rijrkl+2*SIMDD);
        ra = MM_LOAD(akl);
        r1 = MM_LOAD(rijrx+0*SIMDD);
        r2 = MM_LOAD(rijrx+1*SIMDD);
        r3 = MM_LOAD(rijrx+2*SIMDD);
        r0 = MM_MUL(MM_LOAD(tmp1), ra);
        MM_STORE(c00x, MM_FNMA(r0, r4, r1));
        MM_STORE(c00y, MM_FNMA(r0, r5, r2));
        MM_STORE(c00z, MM_FNMA(r0, r6, r3));

        ra = MM_LOAD(aij);
        r1 = MM_LOAD(rklrx+0*SIMDD);
        r2 = MM_LOAD(rklrx+1*SIMDD);
        r3 = MM_LOAD(rklrx+2*SIMDD);
        r0 = MM_MUL(MM_LOAD(tmp1), ra);
        MM_STORE(c0px, MM_FMA (r0, r4, r1));
        MM_STORE(c0py, MM_FMA (r0, r5, r2));
        MM_STORE(c0pz, MM_FMA (r0, r6, r3));

        (*envs->f_g0_2d4d)(g, bc, envs);
        return 1;
}

void CINTg0_2e_gtg_simd1(double *g, Rys2eT *bc, CINTEnvVars *envs, int idsimd)
{
        const double aij = envs->ai[idsimd] + envs->aj[idsimd];
        const double akl = envs->ak[idsimd] + envs->al[idsimd];
        const double zeta = envs->env[PTR_GTG_ZETA];
        double *gx = g;
        double *gy = gx + envs->g_size;
        double *gz = gy + envs->g_size;
        double a0, a1, fac1, x, t;
        double *rij = envs->rij;
        double *rkl = envs->rkl;
        double rijrkl[3];
        double rijrx[3];
        double rklrx[3];
        rijrkl[0] = rij[0*SIMDD+idsimd] - rkl[0*SIMDD+idsimd];
        rijrkl[1] = rij[1*SIMDD+idsimd] - rkl[1*SIMDD+idsimd];
        rijrkl[2] = rij[2*SIMDD+idsimd] - rkl[2*SIMDD+idsimd];

        a1 = aij * akl;
        a0 = a1 / (aij + akl);
        t = zeta / (zeta + a0);
        x = a0 *(rijrkl[0] * rijrkl[0]
               + rijrkl[1] * rijrkl[1]
               + rijrkl[2] * rijrkl[2]);
        fac1 = (1-t) / a1;
        gx[0] = 1;
        gy[0] = 1;
        gz[0] = fac1*sqrt(fac1) * exp(-t * x) * envs->fac[idsimd];
        if (envs->g_size == 1) {
                return 1;
        }

        double div, tmp1, tmp2, tmp3, tmp4;
        double *b00 = bc->b00;
        double *b10 = bc->b10;
        double *b01 = bc->b01;
        double *c00x = bc->c00x;
        double *c00y = bc->c00y;
        double *c00z = bc->c00z;
        double *c0px = bc->c0px;
        double *c0py = bc->c0py;
        double *c0pz = bc->c0pz;

        rijrx[0] = rij[0*SIMDD+idsimd] - envs->rx_in_rijrx[0];
        rijrx[1] = rij[1*SIMDD+idsimd] - envs->rx_in_rijrx[1];
        rijrx[2] = rij[2*SIMDD+idsimd] - envs->rx_in_rijrx[2];
        rklrx[0] = rkl[0*SIMDD+idsimd] - envs->rx_in_rklrx[0];
        rklrx[1] = rkl[1*SIMDD+idsimd] - envs->rx_in_rklrx[1];
        rklrx[2] = rkl[2*SIMDD+idsimd] - envs->rx_in_rklrx[2];
        div = 1 / (zeta * (aij + akl) + a1);
        tmp1 = zeta * div;
        tmp2 = tmp1 * akl;
        tmp3 = tmp1 * aij;
        tmp4 = .5 * div;
        b00[0] = 0.5 * tmp1;
        b10[0] = b00[0] + tmp4 * akl;
        b01[0] = b00[0] + tmp4 * aij;
        c00x[0] = rijrx[0] - tmp2 * rijrkl[0];
        c00y[0] = rijrx[1] - tmp2 * rijrkl[1];
        c00z[0] = rijrx[2] - tmp2 * rijrkl[2];
        c0px[0] = rklrx[0] + tmp3 * rijrkl[0];
        c0py[0] = rklrx[1] + tmp3 * rijrkl[1];
        c0pz[0] = rklrx[2] + tmp3 * rijrkl[2];

        (*envs->f_g0_2d4d_simd1)(g, bc, envs);
        return 1;
}

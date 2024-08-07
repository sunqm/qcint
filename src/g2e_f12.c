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
#include "cint_config.h"
#include "cint_bas.h"
#include "simd.h"
#include "rys_roots.h"
#include "misc.h"
#include "g2e.h"

int CINTg0_2e_yp(double *g, double *cutoff,
                 Rys2eT *bc, CINTEnvVars *envs, int count);
int CINTg0_2e_yp_simd1(double *g, double *cutoff,
                       Rys2eT *bc, CINTEnvVars *envs, int idsimd);
int CINTg0_2e_stg(double *g, double *cutoff,
                  Rys2eT *bc, CINTEnvVars *envs, int count);
int CINTg0_2e_stg_simd1(double *g, double *cutoff,
                        Rys2eT *bc, CINTEnvVars *envs, int idsimd);
void CINTg0_2e_stg_lj2d4d(double *g, Rys2eT *bc, CINTEnvVars *envs);
void CINTg0_2e_stg_lj2d4d_simd1(double *g, Rys2eT *bc, CINTEnvVars *envs);
void CINTg0_lj_4d(double *g, CINTEnvVars *envs);
void CINTg0_lj_4d_simd1(double *g, CINTEnvVars *envs);

void CINTinit_int2e_yp_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                               int *atm, int natm, int *bas, int nbas, double *env)
{
        envs->natm = natm;
        envs->nbas = nbas;
        envs->atm = atm;
        envs->bas = bas;
        envs->env = env;
        envs->shls = shls;

        int i_sh = shls[0];
        int j_sh = shls[1];
        int k_sh = shls[2];
        int l_sh = shls[3];
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
        envs->common_factor = 1;
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

        envs->ri = env + atm(PTR_COORD, bas(ATOM_OF, i_sh));
        envs->rj = env + atm(PTR_COORD, bas(ATOM_OF, j_sh));
        envs->rk = env + atm(PTR_COORD, bas(ATOM_OF, k_sh));
        envs->rl = env + atm(PTR_COORD, bas(ATOM_OF, l_sh));

        // ceil(L_tot/2) + 1
        int nroots = (envs->li_ceil + envs->lj_ceil +
                      envs->lk_ceil + envs->ll_ceil + 3)/2;
        envs->nrys_roots = nroots;
        assert(nroots < MXRYSROOTS);

        int dli, dlj, dlk, dll;
        int ibase = envs->li_ceil > envs->lj_ceil;
        int kbase = envs->lk_ceil > envs->ll_ceil;
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
                        envs->f_g0_2d4d = &CINTg0_2e_stg_lj2d4d;
                        envs->f_g0_2d4d_simd1 = &CINTg0_2e_stg_lj2d4d_simd1;
                }
        }
        envs->f_g0_2e = &CINTg0_2e_yp;
        envs->f_g0_2e_simd1 = &CINTg0_2e_yp_simd1;

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
        envs->g_stride_i = nroots;
        envs->g_stride_k = nroots * dli;
        envs->g_stride_l = nroots * dli * dlk;
        envs->g_stride_j = nroots * dli * dlk * dll;
        envs->g_size     = nroots * dli * dlk * dll * dlj;

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
}

void CINTinit_int2e_stg_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                                int *atm, int natm, int *bas, int nbas, double *env)
{
        CINTinit_int2e_yp_EnvVars(envs, ng, shls, atm, natm, bas, nbas, env);
        envs->f_g0_2e = &CINTg0_2e_stg;
        envs->f_g0_2e_simd1 = &CINTg0_2e_stg_simd1;
}

void CINTg0_2e_stg_lj2d4d(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        CINTg0_2e_2d(g, bc, envs);
        CINTg0_lj_4d(g, envs);
}

void CINTg0_2e_stg_lj2d4d_simd1(double *g, Rys2eT *bc, CINTEnvVars *envs)
{
        CINTg0_2e_2d_simd1(g, bc, envs);
        CINTg0_lj_4d_simd1(g, envs);
}

int CINTg0_2e_yp(double *g, double *cutoff,
                 Rys2eT *bc, CINTEnvVars *envs, int count)
{
        ALIGNMM double aij[SIMDD];
        ALIGNMM double akl[SIMDD];
        ALIGNMM double a0[SIMDD];
        ALIGNMM double a1[SIMDD];
        ALIGNMM double aijkl[SIMDD];
        ALIGNMM double fac1[SIMDD];
        ALIGNMM double x[SIMDD];
        ALIGNMM double rijrkl[SIMDD*3];
        ALIGNMM double rijrx[SIMDD*3];
        ALIGNMM double rklrx[SIMDD*3];
        ALIGNMM double u[MXRYSROOTS*SIMDD];
        double *rij = envs->rij;
        double *rkl = envs->rkl;
        double *w = g + envs->g_size * 2; // ~ gz
        __MD ra, r0, r1, r2, r3, r4, r5, r6, r7, r8;
        double zeta = envs->env[PTR_F12_ZETA];
        int nroots = envs->nrys_roots;
        int i;

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
        MM_STORE(fac1, MM_DIV(MM_LOAD(envs->fac), MM_MUL(MM_SQRT(ra), r1)));

        ALIGNMM double ua[SIMDD];
        //:ua[k] = zeta*zeta / a0[k];
        if (zeta > 0) {
                r0 = MM_SET1(zeta);
                MM_STORE(ua, r0 * r0 * MM_SET1(.25) / MM_LOAD(a0));
        }

        r0 = MM_LOAD(rij+0*SIMDD);
        r1 = MM_LOAD(rij+1*SIMDD);
        r2 = MM_LOAD(rij+2*SIMDD);
        r3 = MM_LOAD(rkl+0*SIMDD);
        r4 = MM_LOAD(rkl+1*SIMDD);
        r5 = MM_LOAD(rkl+2*SIMDD);

        r6 = MM_SUB(r0, r3); MM_STORE(rijrkl+0*SIMDD, r6);
        r7 = MM_SUB(r1, r4); MM_STORE(rijrkl+1*SIMDD, r7);
        r8 = MM_SUB(r2, r5); MM_STORE(rijrkl+2*SIMDD, r8);
        ra = MM_FMA(r6, r6, MM_FMA(r7, r7, MM_MUL(r8, r8)));
        MM_STORE(x, MM_MUL(MM_LOAD(a0), ra));
        MM_STORE(rijrx+0*SIMDD, MM_SUB(r0, MM_SET1(envs->rx_in_rijrx[0])));
        MM_STORE(rijrx+1*SIMDD, MM_SUB(r1, MM_SET1(envs->rx_in_rijrx[1])));
        MM_STORE(rijrx+2*SIMDD, MM_SUB(r2, MM_SET1(envs->rx_in_rijrx[2])));
        MM_STORE(rklrx+0*SIMDD, MM_SUB(r3, MM_SET1(envs->rx_in_rklrx[0])));
        MM_STORE(rklrx+1*SIMDD, MM_SUB(r4, MM_SET1(envs->rx_in_rklrx[1])));
        MM_STORE(rklrx+2*SIMDD, MM_SUB(r5, MM_SET1(envs->rx_in_rklrx[2])));

        double *gx = g;
        double *gy = gx + envs->g_size * SIMDD;
        double *gz = gy + envs->g_size * SIMDD;
        if (zeta > 0) {
                _CINTstg_roots_batch(nroots, x, ua, u, w, count);
                //:w *= t;
                //:uu = t/(1-t);
                r1 = MM_SET1(1.);
                for (i = 0; i < nroots; i++) {
                        r0 = MM_LOAD(u+i*SIMDD);
                        r2 = MM_LOAD(w+i*SIMDD);
                        MM_STORE(w+i*SIMDD, MM_MUL(r0, r2));
                        MM_STORE(u+i*SIMDD, MM_DIV(r0, r1 - r0));
                }
        } else {
                _CINTrys_roots_batch(nroots, x, u, w, count);
        }

        r0 = MM_LOAD(fac1);
        r1 = MM_SET1(1.);
        for (i = 0; i < nroots; i++) {
                //MM_STORE(gx+i*SIMDD, r1);
                //MM_STORE(gy+i*SIMDD, r1);
                MM_STORE(gz+i*SIMDD, MM_MUL(MM_LOAD(w+i*SIMDD), r0));
        }

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

        ALIGNMM double tmp1[MXRYSROOTS*SIMDD];
        ALIGNMM double tmp4[MXRYSROOTS*SIMDD];

        ra = MM_LOAD(aijkl);
        r0 = MM_LOAD(a0);
        r1 = MM_LOAD(a1);
        r2 = MM_SET1(.5);
        r3 = MM_SET1(1.);
        for (i = 0; i < nroots; i++) {
                r4 = MM_MUL(r0, MM_LOAD(u+i*SIMDD));
                r5 = MM_DIV(r3, MM_FMA(r4, ra, r1));
                MM_STORE(tmp4+i*SIMDD, MM_MUL(r2, r5));
                MM_STORE(tmp1+i*SIMDD, MM_MUL(r4, r5));
        }
        ra = MM_SET1(.5);
        r2 = MM_LOAD(akl);
        r3 = MM_LOAD(aij);
        for (i = 0; i < nroots; i++) {
                r0 = MM_MUL(ra, MM_LOAD(tmp1+i*SIMDD));
                MM_STORE(b00+i*SIMDD, r0);
                r1 = MM_LOAD(tmp4+i*SIMDD);
                MM_STORE(b10+i*SIMDD, MM_FMA(r1, r2, r0));
                MM_STORE(b01+i*SIMDD, MM_FMA(r1, r3, r0));
        }
        r4 = MM_LOAD(rijrkl+0*SIMDD);
        r5 = MM_LOAD(rijrkl+1*SIMDD);
        r6 = MM_LOAD(rijrkl+2*SIMDD);
        ra = MM_LOAD(akl);
        r1 = MM_LOAD(rijrx+0*SIMDD);
        r2 = MM_LOAD(rijrx+1*SIMDD);
        r3 = MM_LOAD(rijrx+2*SIMDD);
        for (i = 0; i < nroots; i++) {
                r0 = MM_MUL(MM_LOAD(tmp1+i*SIMDD), ra);
                MM_STORE(c00x+i*SIMDD, MM_FNMA(r0, r4, r1));
                MM_STORE(c00y+i*SIMDD, MM_FNMA(r0, r5, r2));
                MM_STORE(c00z+i*SIMDD, MM_FNMA(r0, r6, r3));
        }
        ra = MM_LOAD(aij);
        r1 = MM_LOAD(rklrx+0*SIMDD);
        r2 = MM_LOAD(rklrx+1*SIMDD);
        r3 = MM_LOAD(rklrx+2*SIMDD);
        for (i = 0; i < nroots; i++) {
                r0 = MM_MUL(MM_LOAD(tmp1+i*SIMDD), ra);
                MM_STORE(c0px+i*SIMDD, MM_FMA (r0, r4, r1));
                MM_STORE(c0py+i*SIMDD, MM_FMA (r0, r5, r2));
                MM_STORE(c0pz+i*SIMDD, MM_FMA (r0, r6, r3));
        }

        (*envs->f_g0_2d4d)(g, bc, envs);
        return 1;
}

int CINTg0_2e_yp_simd1(double *g, double *cutoff,
                       Rys2eT *bc, CINTEnvVars *envs, int idsimd)
{
        double aij = envs->ai[idsimd] + envs->aj[idsimd];
        double akl = envs->ak[idsimd] + envs->al[idsimd];
        double zeta = envs->env[PTR_F12_ZETA];
        int nroots = envs->nrys_roots;
        double a0, a1, fac1, x, ua;
        double *rij = envs->rij;
        double *rkl = envs->rkl;
        double rijrkl[3];
        double rijrx[3];
        double rklrx[3];
        ALIGNMM double u[MXRYSROOTS];
        double *w = g + envs->g_size * 2; // ~ gz
        int i;

        a1 = aij * akl;
        a0 = a1 / (aij + akl);
        //fac1 = sqrt(a0 / (a1 * a1 * a1)) * envs->fac[idsimd];
        fac1 = envs->fac[idsimd] / (sqrt(aij+akl) * a1);
        rijrkl[0] = rij[0*SIMDD+idsimd] - rkl[0*SIMDD+idsimd];
        rijrkl[1] = rij[1*SIMDD+idsimd] - rkl[1*SIMDD+idsimd];
        rijrkl[2] = rij[2*SIMDD+idsimd] - rkl[2*SIMDD+idsimd];
        x = a0 * SQUARE(rijrkl);

        //:ua[k] = zeta*zeta / a0[k];
        if (zeta > 0) {
                ua = .25 * zeta * zeta / a0;
                CINTstg_roots(nroots, x, ua, u, w);
                //:w *= t;
                //:uu = t/(1-t);
                for (i = 0; i < nroots; i++) {
                        w[i] *= u[i];
                        u[i] = u[i] / (1 - u[i]);
                }
        } else {
                CINTrys_roots(nroots, x, u, w);
        }

        double *gx = g;
        double *gy = g + envs->g_size;
        double *gz = g + envs->g_size * 2;
        for (i = 0; i < nroots; i++) {
                gx[i] = 1;
                gy[i] = 1;
                gz[i] = w[i] * fac1;
        }

        if (envs->g_size == 1) {
                return 1;
        }

        double u2, div, tmp1, tmp2, tmp3, tmp4;
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
        for (i = 0; i < nroots; i++) {
                /*
                 *u(i) = t2/(1-t2)
                 *t2 = u(i)/(1+u(i))
                 *u2 = aij*akl/(aij+akl)*t2/(1-t2)
                 */
                u2 = a0 * u[i];
                div = 1 / (u2 * (aij + akl) + a1);
                tmp1 = u2 * div;
                tmp4 = .5 * div;
                b00[i] = 0.5 * tmp1;
                tmp2 = tmp1 * akl;
                tmp3 = tmp1 * aij;
                b10[i] = b00[i] + tmp4 * akl;
                b01[i] = b00[i] + tmp4 * aij;
                c00x[i] = rijrx[0] - tmp2 * rijrkl[0];
                c00y[i] = rijrx[1] - tmp2 * rijrkl[1];
                c00z[i] = rijrx[2] - tmp2 * rijrkl[2];
                c0px[i] = rklrx[0] + tmp3 * rijrkl[0];
                c0py[i] = rklrx[1] + tmp3 * rijrkl[1];
                c0pz[i] = rklrx[2] + tmp3 * rijrkl[2];
        }

        (*envs->f_g0_2d4d_simd1)(g, bc, envs);
        return 1;
}

int CINTg0_2e_stg(double *g, double *cutoff,
                  Rys2eT *bc, CINTEnvVars *envs, int count)
{
        ALIGNMM double aij[SIMDD];
        ALIGNMM double akl[SIMDD];
        ALIGNMM double a0[SIMDD];
        ALIGNMM double a1[SIMDD];
        ALIGNMM double aijkl[SIMDD];
        ALIGNMM double fac1[SIMDD];
        ALIGNMM double x[SIMDD];
        ALIGNMM double ua[SIMDD];
        ALIGNMM double rijrkl[SIMDD*3];
        ALIGNMM double rijrx[SIMDD*3];
        ALIGNMM double rklrx[SIMDD*3];
        ALIGNMM double u[MXRYSROOTS*SIMDD];
        double *rij = envs->rij;
        double *rkl = envs->rkl;
        double *w = g + envs->g_size * 2; // ~ gz
        __MD ra, r0, r1, r2, r3, r4, r5, r6, r7, r8;
        double zeta = envs->env[PTR_F12_ZETA];
        int nroots = envs->nrys_roots;
        int i;

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
        MM_STORE(fac1, MM_DIV(MM_LOAD(envs->fac), MM_MUL(MM_SQRT(ra), r1)));

        r0 = MM_LOAD(rij+0*SIMDD);
        r1 = MM_LOAD(rij+1*SIMDD);
        r2 = MM_LOAD(rij+2*SIMDD);
        r3 = MM_LOAD(rkl+0*SIMDD);
        r4 = MM_LOAD(rkl+1*SIMDD);
        r5 = MM_LOAD(rkl+2*SIMDD);

        r6 = MM_SUB(r0, r3); MM_STORE(rijrkl+0*SIMDD, r6);
        r7 = MM_SUB(r1, r4); MM_STORE(rijrkl+1*SIMDD, r7);
        r8 = MM_SUB(r2, r5); MM_STORE(rijrkl+2*SIMDD, r8);
        ra = MM_FMA(r6, r6, MM_FMA(r7, r7, MM_MUL(r8, r8)));
        MM_STORE(x, MM_MUL(MM_LOAD(a0), ra));
        MM_STORE(rijrx+0*SIMDD, MM_SUB(r0, MM_SET1(envs->rx_in_rijrx[0])));
        MM_STORE(rijrx+1*SIMDD, MM_SUB(r1, MM_SET1(envs->rx_in_rijrx[1])));
        MM_STORE(rijrx+2*SIMDD, MM_SUB(r2, MM_SET1(envs->rx_in_rijrx[2])));
        MM_STORE(rklrx+0*SIMDD, MM_SUB(r3, MM_SET1(envs->rx_in_rklrx[0])));
        MM_STORE(rklrx+1*SIMDD, MM_SUB(r4, MM_SET1(envs->rx_in_rklrx[1])));
        MM_STORE(rklrx+2*SIMDD, MM_SUB(r5, MM_SET1(envs->rx_in_rklrx[2])));

        //:ua[k] = zeta*zeta / a0[k];
        if (zeta > 0) {
                r0 = MM_SET1(zeta);
                MM_STORE(ua, r0 * r0 * MM_SET1(.25) / MM_LOAD(a0));
                // T = x, U = ua, rr = u
                _CINTstg_roots_batch(nroots, x, ua, u, w, count);
                //:w *= (1-t) * 2*ua/zeta;
                //:uu = t/(1-t);
                r1 = MM_SET1(1.);
                r2 = MM_DIV(MM_LOAD(ua) * MM_SET1(2.), MM_SET1(zeta));
                for (i = 0; i < nroots; i++) {
                        r0 = MM_LOAD(u+i*SIMDD);
                        r3 = r1 - r0;
                        r4 = MM_LOAD(w+i*SIMDD);
                        MM_STORE(w+i*SIMDD, r4 * r3 * r2);
                        MM_STORE(u+i*SIMDD, MM_DIV(r0, r3));
                }
        } else {
                _CINTrys_roots_batch(nroots, x, u, w, count);
        }

        double *gx = g;
        double *gy = gx + envs->g_size * SIMDD;
        double *gz = gy + envs->g_size * SIMDD;
        r0 = MM_LOAD(fac1);
        r1 = MM_SET1(1.);
        for (i = 0; i < nroots; i++) {
                //MM_STORE(gx+i*SIMDD, r1);
                //MM_STORE(gy+i*SIMDD, r1);
                MM_STORE(gz+i*SIMDD, MM_MUL(MM_LOAD(w+i*SIMDD), r0));
        }

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

        ALIGNMM double tmp1[MXRYSROOTS*SIMDD];
        ALIGNMM double tmp4[MXRYSROOTS*SIMDD];

        ra = MM_LOAD(aijkl);
        r0 = MM_LOAD(a0);
        r1 = MM_LOAD(a1);
        r2 = MM_SET1(.5);
        r3 = MM_SET1(1.);
        for (i = 0; i < nroots; i++) {
                r4 = MM_MUL(r0, MM_LOAD(u+i*SIMDD));
                r5 = MM_DIV(r3, MM_FMA(r4, ra, r1));
                MM_STORE(tmp4+i*SIMDD, MM_MUL(r2, r5));
                MM_STORE(tmp1+i*SIMDD, MM_MUL(r4, r5));
        }
        ra = MM_SET1(.5);
        r2 = MM_LOAD(akl);
        r3 = MM_LOAD(aij);
        for (i = 0; i < nroots; i++) {
                r0 = MM_MUL(ra, MM_LOAD(tmp1+i*SIMDD));
                MM_STORE(b00+i*SIMDD, r0);
                r1 = MM_LOAD(tmp4+i*SIMDD);
                MM_STORE(b10+i*SIMDD, MM_FMA(r1, r2, r0));
                MM_STORE(b01+i*SIMDD, MM_FMA(r1, r3, r0));
        }
        r4 = MM_LOAD(rijrkl+0*SIMDD);
        r5 = MM_LOAD(rijrkl+1*SIMDD);
        r6 = MM_LOAD(rijrkl+2*SIMDD);
        ra = MM_LOAD(akl);
        r1 = MM_LOAD(rijrx+0*SIMDD);
        r2 = MM_LOAD(rijrx+1*SIMDD);
        r3 = MM_LOAD(rijrx+2*SIMDD);
        for (i = 0; i < nroots; i++) {
                r0 = MM_MUL(MM_LOAD(tmp1+i*SIMDD), ra);
                MM_STORE(c00x+i*SIMDD, MM_FNMA(r0, r4, r1));
                MM_STORE(c00y+i*SIMDD, MM_FNMA(r0, r5, r2));
                MM_STORE(c00z+i*SIMDD, MM_FNMA(r0, r6, r3));
        }
        ra = MM_LOAD(aij);
        r1 = MM_LOAD(rklrx+0*SIMDD);
        r2 = MM_LOAD(rklrx+1*SIMDD);
        r3 = MM_LOAD(rklrx+2*SIMDD);
        for (i = 0; i < nroots; i++) {
                r0 = MM_MUL(MM_LOAD(tmp1+i*SIMDD), ra);
                MM_STORE(c0px+i*SIMDD, MM_FMA (r0, r4, r1));
                MM_STORE(c0py+i*SIMDD, MM_FMA (r0, r5, r2));
                MM_STORE(c0pz+i*SIMDD, MM_FMA (r0, r6, r3));
        }

        (*envs->f_g0_2d4d)(g, bc, envs);
        return 1;
}

int CINTg0_2e_stg_simd1(double *g, double *cutoff,
                        Rys2eT *bc, CINTEnvVars *envs, int idsimd)
{
        double aij = envs->ai[idsimd] + envs->aj[idsimd];
        double akl = envs->ak[idsimd] + envs->al[idsimd];
        double zeta = envs->env[PTR_F12_ZETA];
        int nroots = envs->nrys_roots;
        double a0, a1, fac1, x, ua;
        double *rij = envs->rij;
        double *rkl = envs->rkl;
        double rijrkl[3];
        double rijrx[3];
        double rklrx[3];
        ALIGNMM double u[MXRYSROOTS*SIMDD];
        double *w = g + envs->g_size * 2; // ~ gz
        int i;

        a1 = aij * akl;
        a0 = a1 / (aij + akl);
        //fac1 = sqrt(a0 / (a1 * a1 * a1)) * envs->fac[idsimd];
        fac1 = envs->fac[idsimd] / (sqrt(aij+akl) * a1);

        rijrkl[0] = rij[0*SIMDD+idsimd] - rkl[0*SIMDD+idsimd];
        rijrkl[1] = rij[1*SIMDD+idsimd] - rkl[1*SIMDD+idsimd];
        rijrkl[2] = rij[2*SIMDD+idsimd] - rkl[2*SIMDD+idsimd];
        x = a0 * SQUARE(rijrkl);

        //:ua[k] = zeta*zeta / a0[k];
        if (zeta > 0) {
                ua = .25 * zeta * zeta / a0;
                CINTstg_roots(nroots, x, ua, u, w);
                //:w *= (1-t) * 2*ua/zeta;
                //:uu = t/(1-t);
                double ua2 = 2. * ua / zeta;
                for (i = 0; i < nroots; i++) {
                        w[i] *= (1-u[i]) * ua2;
                        u[i] = u[i] / (1 - u[i]);
                }
        } else {
                CINTrys_roots(nroots, x, u, w);
        }

        double *gx = g;
        double *gy = g + envs->g_size;
        double *gz = g + envs->g_size * 2;
        for (i = 0; i < nroots; i++) {
                gx[i] = 1;
                gy[i] = 1;
                gz[i] = w[i] * fac1;
        }

        if (envs->g_size == 1) {
                return 1;
        }

        double u2, div, tmp1, tmp2, tmp3, tmp4;
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
        for (i = 0; i < nroots; i++) {
                /*
                 *t2 = u(i)/(1+u(i))
                 *u2 = aij*akl/(aij+akl)*t2/(1-t2)
                 */
                u2 = a0 * u[i];
                div = 1 / (u2 * (aij + akl) + a1);
                tmp1 = u2 * div;
                tmp4 = .5 * div;
                b00[i] = 0.5 * tmp1;
                tmp2 = tmp1 * akl;
                tmp3 = tmp1 * aij;
                b10[i] = b00[i] + tmp4 * akl;
                b01[i] = b00[i] + tmp4 * aij;
                c00x[i] = rijrx[0] - tmp2 * rijrkl[0];
                c00y[i] = rijrx[1] - tmp2 * rijrkl[1];
                c00z[i] = rijrx[2] - tmp2 * rijrkl[2];
                c0px[i] = rklrx[0] + tmp3 * rijrkl[0];
                c0py[i] = rklrx[1] + tmp3 * rijrkl[1];
                c0pz[i] = rklrx[2] + tmp3 * rijrkl[2];
        }

        (*envs->f_g0_2d4d_simd1)(g, bc, envs);
        return 1;
}



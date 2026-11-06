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

#include <string.h>
#include <math.h>
#include "cint_bas.h"
#include "simd.h"
#include "misc.h"
#include "g1e.h"
#include "g2e.h"

void CINTg0_lj_4d(double *g, CINTEnvVars *envs);
void CINTg0_kj_4d(double *g, CINTEnvVars *envs);
void CINTg0_il_4d(double *g, CINTEnvVars *envs);
void CINTg0_ik_4d(double *g, CINTEnvVars *envs);

void CINTinit_int4c1e_EnvVars(CINTEnvVars *envs, int *ng, FINT *shls,
                              FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env)
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

        envs->ri = env + atm(PTR_COORD, bas(ATOM_OF, i_sh));
        envs->rj = env + atm(PTR_COORD, bas(ATOM_OF, j_sh));
        envs->rk = env + atm(PTR_COORD, bas(ATOM_OF, k_sh));
        envs->rl = env + atm(PTR_COORD, bas(ATOM_OF, l_sh));

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
        envs->rys_order = 1;
        envs->nrys_roots = 1;

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
        envs->g_stride_i = 1;
        envs->g_stride_k = dli;
        envs->g_stride_l = dli * dlk;
        envs->g_stride_j = dli * dlk * dll;
        envs->g_size     = dli * dlk * dll * dlj;

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
                        envs->f_g0_2d4d = &CINTg0_ik_4d;
                } else {
                        envs->f_g0_2d4d = &CINTg0_kj_4d;
                }
        } else {
                if (ibase) {
                        envs->f_g0_2d4d = &CINTg0_il_4d;
                } else {
                        envs->f_g0_2d4d = &CINTg0_lj_4d;
                }
        }
}

/*
 * SIMD version of the 4-center 1-electron overlap-type g matrix.
 * Each SIMD lane holds one primitive quartet.  The per-lane exponents
 * and pair centers are read from envs->ai,aj,ak,al and envs->rij,rkl.
 * buf provides workspace of 3 * db*(MAX(nmax,mmax)+1) * SIMDD doubles,
 * with db = nmax+mmax+1, nmax = li_ceil+lj_ceil, mmax = lk_ceil+ll_ceil.
 */
void CINTg4c1e_ovlp(double *g, double *buf, CINTEnvVars *envs, int count)
{
        int li = envs->li_ceil;
        int lj = envs->lj_ceil;
        int lk = envs->lk_ceil;
        int ll = envs->ll_ceil;
        int nmax = li + lj;
        int mmax = lk + ll;
        int db = nmax + mmax + 1;
        int b_size = db * (MAX(nmax, mmax) + 1);
        double *bufx = buf;
        double *bufy = bufx + b_size * SIMDD;
        double *bufz = bufy + b_size * SIMDD;
        double *rij = envs->rij;
        double *rkl = envs->rkl;
        double *r1;
        double *r2;
        __MD aij = MM_LOAD(envs->ai) + MM_LOAD(envs->aj);
        __MD akl = MM_LOAD(envs->ak) + MM_LOAD(envs->al);
        __MD aijkl = aij + akl;
        __MD rt = MM_SET1(1.) / aijkl;
        __MD half = MM_SET1(.5) * rt;
        __MD r12r1[3];  // (aij*rij+akl*rkl)/aijkl - r1, the negative of libcint's r1r12
        __MD r1r2[3];
        __MD ren;
        int i, j, ptr;

        double *gx = g;
        double *gy = g + envs->g_size     * SIMDD;
        double *gz = g + envs->g_size * 2 * SIMDD;
        int dn = envs->g2d_ijmax;
        int dm = envs->g2d_klmax;

        MM_STORE(bufx, MM_SET1(1.));
        MM_STORE(bufy, MM_SET1(1.));
        MM_STORE(bufz, MM_LOAD(envs->fac) / (aijkl * MM_SQRT(aijkl)));

        if (nmax >= mmax) {
                r1 = envs->rx_in_rijrx;
                r2 = envs->rx_in_rklrx;
        } else {
                r1 = envs->rx_in_rklrx;
                r2 = envs->rx_in_rijrx;
        }
        for (i = 0; i < 3; i++) {
                r12r1[i] = (aij * MM_LOAD(rij+i*SIMDD)
                          + akl * MM_LOAD(rkl+i*SIMDD)) * rt - MM_SET1(r1[i]);
                r1r2[i] = MM_SET1(r1[i] - r2[i]);
        }

        if (nmax + mmax > 0) {
                MM_STORE(bufx+SIMDD, r12r1[0] * MM_LOAD(bufx));
                MM_STORE(bufy+SIMDD, r12r1[1] * MM_LOAD(bufy));
                MM_STORE(bufz+SIMDD, r12r1[2] * MM_LOAD(bufz));
        }

        for (i = 1; i < nmax+mmax; i++) {
                ren = MM_SET1(i) * half;
                MM_STORE(bufx+(i+1)*SIMDD, ren * MM_LOAD(bufx+(i-1)*SIMDD) + r12r1[0] * MM_LOAD(bufx+i*SIMDD));
                MM_STORE(bufy+(i+1)*SIMDD, ren * MM_LOAD(bufy+(i-1)*SIMDD) + r12r1[1] * MM_LOAD(bufy+i*SIMDD));
                MM_STORE(bufz+(i+1)*SIMDD, ren * MM_LOAD(bufz+(i-1)*SIMDD) + r12r1[2] * MM_LOAD(bufz+i*SIMDD));
        }

        // horizontal transfer along the axis of the second center
        for (j = 1; j <= MIN(nmax, mmax); j++) {
                ptr = db * j;
                for (i = ptr; i <= ptr + nmax+mmax - j; i++) {
                        MM_STORE(bufx+i*SIMDD, MM_LOAD(bufx+(i+1-db)*SIMDD) + r1r2[0] * MM_LOAD(bufx+(i-db)*SIMDD));
                        MM_STORE(bufy+i*SIMDD, MM_LOAD(bufy+(i+1-db)*SIMDD) + r1r2[1] * MM_LOAD(bufy+(i-db)*SIMDD));
                        MM_STORE(bufz+i*SIMDD, MM_LOAD(bufz+(i+1-db)*SIMDD) + r1r2[2] * MM_LOAD(bufz+(i-db)*SIMDD));
                }
        }

        if (nmax >= mmax) {
                for (j = 0; j <= mmax; j++) {
                for (i = 0; i <= nmax; i++) {
                        MM_STORE(gx+(i*dn+j*dm)*SIMDD, MM_LOAD(bufx+(i+j*db)*SIMDD));
                        MM_STORE(gy+(i*dn+j*dm)*SIMDD, MM_LOAD(bufy+(i+j*db)*SIMDD));
                        MM_STORE(gz+(i*dn+j*dm)*SIMDD, MM_LOAD(bufz+(i+j*db)*SIMDD));
                } }
        } else {
                for (j = 0; j <= mmax; j++) {
                for (i = 0; i <= nmax; i++) {
                        MM_STORE(gx+(i*dn+j*dm)*SIMDD, MM_LOAD(bufx+(i*db+j)*SIMDD));
                        MM_STORE(gy+(i*dn+j*dm)*SIMDD, MM_LOAD(bufy+(i*db+j)*SIMDD));
                        MM_STORE(gz+(i*dn+j*dm)*SIMDD, MM_LOAD(bufz+(i*db+j)*SIMDD));
                } }
        }

        (*envs->f_g0_2d4d)(g, envs);
}

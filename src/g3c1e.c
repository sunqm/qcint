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

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "cint_bas.h"
#include "simd.h"
#include "misc.h"
#include "g1e.h"
#include "rys_roots.h"

void CINTinit_int3c1e_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
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
        envs->i_l = bas(ANG_OF, i_sh);
        envs->j_l = bas(ANG_OF, j_sh);
        envs->k_l = bas(ANG_OF, k_sh);
        envs->l_l = 0;
        envs->x_ctr[0] = bas(NCTR_OF, i_sh);
        envs->x_ctr[1] = bas(NCTR_OF, j_sh);
        envs->x_ctr[2] = bas(NCTR_OF, k_sh);
        envs->x_ctr[3] = 1;
        envs->nfi = (envs->i_l+1)*(envs->i_l+2)/2;
        envs->nfj = (envs->j_l+1)*(envs->j_l+2)/2;
        envs->nfk = (envs->k_l+1)*(envs->k_l+2)/2;
        envs->nfl = 1;
        envs->nf = envs->nfi * envs->nfk * envs->nfj;
        envs->common_factor = 1;

        envs->gbits = ng[GSHIFT];
        envs->ncomp_e1 = ng[POS_E1];
        envs->ncomp_e2 = 0;
        envs->ncomp_tensor = ng[TENSOR];

        envs->li_ceil = envs->i_l + ng[IINC];
        envs->lj_ceil = envs->j_l + ng[JINC];
        envs->lk_ceil = envs->k_l + ng[KINC];
        envs->ll_ceil = 0;

        envs->ri = env + atm(PTR_COORD, bas(ATOM_OF, i_sh));
        envs->rj = env + atm(PTR_COORD, bas(ATOM_OF, j_sh));
        envs->rk = env + atm(PTR_COORD, bas(ATOM_OF, k_sh));

        if (ng[RYS_ROOTS] == 0) {
                envs->nrys_roots = (envs->li_ceil + envs->lj_ceil + envs->lk_ceil)/2 + 1;
        } else {
                envs->nrys_roots = ng[RYS_ROOTS];
        }
        int dli = envs->li_ceil + 1;
        int dlj = envs->lj_ceil + envs->lk_ceil + 1;
        int dlk = envs->lk_ceil + 1;
        int nmax = envs->li_ceil + dlj;
        envs->g_stride_i = envs->nrys_roots;
        envs->g_stride_j = dli * envs->nrys_roots;
        envs->g_stride_k = dli * dlj * envs->nrys_roots;
        envs->g_stride_l = envs->g_stride_k;
        envs->g_size     = MAX(dli*dlj*dlk, dli*nmax) * envs->nrys_roots;

        envs->rirj[0] = envs->ri[0] - envs->rj[0];
        envs->rirj[1] = envs->ri[1] - envs->rj[1];
        envs->rirj[2] = envs->ri[2] - envs->rj[2];
        envs->rkrl[0] = envs->rj[0] - envs->rk[0];
        envs->rkrl[1] = envs->rj[1] - envs->rk[1];
        envs->rkrl[2] = envs->rj[2] - envs->rk[2];
}

void CINTg3c1e_ovlp(double *g, CINTEnvVars *envs, int count)
{
        int li = envs->li_ceil;
        int lj = envs->lj_ceil;
        int lk = envs->lk_ceil;
        int nmax = li + lj + lk;
        int mmax = lj + lk;
        double *gx = g;
        double *gy = g + envs->g_size     * SIMDD;
        double *gz = g + envs->g_size * 2 * SIMDD;

        __MD aijk;
        aijk = MM_LOAD(envs->ai) + MM_LOAD(envs->aj) + MM_LOAD(envs->ak);
        MM_STORE(gx, MM_SET1(1.));
        MM_STORE(gy, MM_SET1(1.));
        MM_STORE(gz, MM_LOAD(envs->fac) / (aijk * MM_SQRT(aijk)));

        if (nmax == 0) {
                return;
        }

        int dj = envs->g_stride_j;
        int dk = envs->g_stride_k;
        double *ri = envs->ri;
        double *rj = envs->rj;
        double *rk = envs->rk;
        double *rirj = envs->rirj;
        double *rjrk = envs->rkrl;
        __MD r0, r1, r2, rt;

        rt = MM_SET1(1.) / aijk;
        r0 = MM_SET1(rj[0]);
        r1 = MM_SET1(rj[1]);
        r2 = MM_SET1(rj[2]);
        r0 = (MM_LOAD(envs->ai) * MM_SET1(ri[0]) + MM_LOAD(envs->aj) * r0 + MM_LOAD(envs->ak) * MM_SET1(rk[0])) * rt - r0;
        r1 = (MM_LOAD(envs->ai) * MM_SET1(ri[1]) + MM_LOAD(envs->aj) * r1 + MM_LOAD(envs->ak) * MM_SET1(rk[1])) * rt - r1;
        r2 = (MM_LOAD(envs->ai) * MM_SET1(ri[2]) + MM_LOAD(envs->aj) * r2 + MM_LOAD(envs->ak) * MM_SET1(rk[2])) * rt - r2;

        MM_STORE(gx+dj*SIMDD, r0 * MM_LOAD(gx+0));
        MM_STORE(gy+dj*SIMDD, r1 * MM_LOAD(gy+0));
        MM_STORE(gz+dj*SIMDD, r2 * MM_LOAD(gz+0));

        int i, j, k, off;
        double *p0x, *p0y, *p0z;
        double *p1x, *p1y, *p1z;
        double *p2x, *p2y, *p2z;
        p0x = gx + dj * SIMDD;
        p0y = gy + dj * SIMDD;
        p0z = gz + dj * SIMDD;
        p1x = gx - dj * SIMDD;
        p1y = gy - dj * SIMDD;
        p1z = gz - dj * SIMDD;

        aijk = MM_SET1(.5) / aijk;
        for (j = 1; j < nmax; j++) {
                MM_STORE(p0x+j*dj*SIMDD, aijk * MM_SET1(j) * MM_LOAD(p1x+j*dj*SIMDD) + r0 * MM_LOAD(gx+j*dj*SIMDD));
                MM_STORE(p0y+j*dj*SIMDD, aijk * MM_SET1(j) * MM_LOAD(p1y+j*dj*SIMDD) + r1 * MM_LOAD(gy+j*dj*SIMDD));
                MM_STORE(p0z+j*dj*SIMDD, aijk * MM_SET1(j) * MM_LOAD(p1z+j*dj*SIMDD) + r2 * MM_LOAD(gz+j*dj*SIMDD));
        }

        r0 = MM_SET1(rirj[0]);
        r1 = MM_SET1(rirj[1]);
        r2 = MM_SET1(rirj[2]);

        for (i = 1; i <= li; i++) {
        for (j = 0; j <= nmax-i; j++) { // upper limit lj+lk
                MM_STORE(gx+(i+j*dj)*SIMDD, MM_LOAD(p0x+(i-1+j*dj)*SIMDD) - r0 * MM_LOAD(gx+(i-1+j*dj)*SIMDD));
                MM_STORE(gy+(i+j*dj)*SIMDD, MM_LOAD(p0y+(i-1+j*dj)*SIMDD) - r1 * MM_LOAD(gy+(i-1+j*dj)*SIMDD));
                MM_STORE(gz+(i+j*dj)*SIMDD, MM_LOAD(p0z+(i-1+j*dj)*SIMDD) - r2 * MM_LOAD(gz+(i-1+j*dj)*SIMDD));
        } }

        r0 = MM_SET1(rjrk[0]);
        r1 = MM_SET1(rjrk[1]);
        r2 = MM_SET1(rjrk[2]);
        for (k = 1; k <= lk; k++) {
        for (j = 0; j <= mmax-k; j++) {
                off = k * dk + j * dj;
                p0x = gx + off * SIMDD;
                p0y = gy + off * SIMDD;
                p0z = gz + off * SIMDD;
                p2x = p0x - dk * SIMDD;
                p2y = p0y - dk * SIMDD;
                p2z = p0z - dk * SIMDD;
                p1x = p2x + dj * SIMDD;
                p1y = p2y + dj * SIMDD;
                p1z = p2z + dj * SIMDD;
                for (i = 0; i <= li; i++) {
                        MM_STORE(p0x+i*SIMDD, MM_LOAD(p1x+i*SIMDD) + r0 * MM_LOAD(p2x+i*SIMDD));
                        MM_STORE(p0y+i*SIMDD, MM_LOAD(p1y+i*SIMDD) + r1 * MM_LOAD(p2y+i*SIMDD));
                        MM_STORE(p0z+i*SIMDD, MM_LOAD(p1z+i*SIMDD) + r2 * MM_LOAD(p2z+i*SIMDD));
                }
        } }
}


void CINTg3c1e_nuc(double *g, CINTEnvVars *envs, int count, int nuc_id)
{
        int li = envs->li_ceil;
        int lj = envs->lj_ceil;
        int lk = envs->lk_ceil;
        int nmax = li + lj + lk;
        int mmax = lj + lk;
        int nrys_roots = envs->nrys_roots;
        int *atm = envs->atm;
        double *env = envs->env;
        int di = envs->g_stride_i;
        int dj = envs->g_stride_j;
        int dk = envs->g_stride_k;
        double *ri = envs->ri;
        double *rj = envs->rj;
        double *rk = envs->rk;
        double *rirj = envs->rirj;
        double *rjrk = envs->rkrl;
        ALIGNMM double tau[SIMDD];
        ALIGNMM double x[SIMDD];
        ALIGNMM double u[MXRYSROOTS*SIMDD];
        ALIGNMM double w[MXRYSROOTS*SIMDD];
        ALIGNMM double t2[MXRYSROOTS*SIMDD];
        ALIGNMM double r0jx[MXRYSROOTS*SIMDD];
        ALIGNMM double r0jy[MXRYSROOTS*SIMDD];
        ALIGNMM double r0jz[MXRYSROOTS*SIMDD];
        ALIGNMM double rijk[3*SIMDD];
        double *cr, *p0x, *p0y, *p0z, *p1x, *p1y, *p1z, *p2x, *p2y, *p2z;
        int i, j, k, n, ptr, off;
        __MD crijk[3];
        __MD r0, r1, r2, rt2, fac1, aijk;

        aijk = MM_LOAD(envs->ai) + MM_LOAD(envs->aj) + MM_LOAD(envs->ak);
        MM_STORE(tau, aijk);
        for (k = 0; k < count; k++) {
                tau[k] = CINTnuc_mod(tau[k], nuc_id, atm, env);
        }

        if (nuc_id < 0) {
                fac1 = 2./SQRTPI * MM_LOAD(envs->fac) * MM_LOAD(tau) / aijk;
                cr = env + PTR_RINV_ORIG;
        } else {
                fac1 = 2./SQRTPI * MM_SET1(-fabs(atm[CHARGE_OF+nuc_id*ATM_SLOTS]));
                fac1 = fac1 * MM_LOAD(envs->fac) * MM_LOAD(tau) / aijk;
                cr = env + atm(PTR_COORD, nuc_id);
        }

        r0 = MM_LOAD(envs->ai) * MM_SET1(ri[0]) + MM_LOAD(envs->aj) * MM_SET1(rj[0]) + MM_LOAD(envs->ak) * MM_SET1(rk[0]);
        r1 = MM_LOAD(envs->ai) * MM_SET1(ri[1]) + MM_LOAD(envs->aj) * MM_SET1(rj[1]) + MM_LOAD(envs->ak) * MM_SET1(rk[1]);
        r2 = MM_LOAD(envs->ai) * MM_SET1(ri[2]) + MM_LOAD(envs->aj) * MM_SET1(rj[2]) + MM_LOAD(envs->ak) * MM_SET1(rk[2]);
        MM_STORE(rijk+0*SIMDD, r0 / aijk);
        MM_STORE(rijk+1*SIMDD, r1 / aijk);
        MM_STORE(rijk+2*SIMDD, r2 / aijk);
        crijk[0] = MM_SET1(cr[0]) - MM_LOAD(rijk+0*SIMDD);
        crijk[1] = MM_SET1(cr[1]) - MM_LOAD(rijk+1*SIMDD);
        crijk[2] = MM_SET1(cr[2]) - MM_LOAD(rijk+2*SIMDD);
        MM_STORE(x, aijk * MM_LOAD(tau) * MM_LOAD(tau) * SQUARE(crijk));
        for (i = 0; i < nrys_roots*SIMDD; i++) {
                u[i] = 0;
                w[i] = 0;
        }
        CINTrys_roots(nrys_roots, x, u, w, count);

        double *gx = g;
        double *gy = g + envs->g_size     * SIMDD;
        double *gz = g + envs->g_size * 2 * SIMDD;
        r1 = MM_SET1(1.);
        for (i = 0; i < nrys_roots; i++) {
                MM_STORE(gx+i*SIMDD, r1);
                MM_STORE(gy+i*SIMDD, r1);
                MM_STORE(gz+i*SIMDD, fac1 * MM_LOAD(w+i*SIMDD));
        }
        if (envs->g_size == 1) {
                return;
        }

        r1 = MM_SET1(1.) * MM_LOAD(tau) * MM_LOAD(tau);
        r0 = MM_SET1(0.5);
        for (i = 0; i < nrys_roots; i++) {
                rt2 = MM_DIV(r1 * MM_LOAD(u+i*SIMDD), r1 + MM_LOAD(u+i*SIMDD));
                MM_STORE(r0jx+i*SIMDD, MM_LOAD(rijk+0*SIMDD) + rt2 * crijk[0] - MM_SET1(rj[0]));
                MM_STORE(r0jy+i*SIMDD, MM_LOAD(rijk+1*SIMDD) + rt2 * crijk[1] - MM_SET1(rj[1]));
                MM_STORE(r0jz+i*SIMDD, MM_LOAD(rijk+2*SIMDD) + rt2 * crijk[2] - MM_SET1(rj[2]));
                MM_STORE(t2+i*SIMDD, MM_DIV(r0 - r0 * rt2, aijk));
        }

        p0x = gx + dj * SIMDD;
        p0y = gy + dj * SIMDD;
        p0z = gz + dj * SIMDD;
        for (n = 0; n < nrys_roots; n++) {
MM_STORE(p0x+n*SIMDD, MM_LOAD(r0jx+n*SIMDD) * MM_LOAD(gx+n*SIMDD));
MM_STORE(p0y+n*SIMDD, MM_LOAD(r0jy+n*SIMDD) * MM_LOAD(gy+n*SIMDD));
MM_STORE(p0z+n*SIMDD, MM_LOAD(r0jz+n*SIMDD) * MM_LOAD(gz+n*SIMDD));
        }

        p0x = gx + dj * SIMDD;
        p0y = gy + dj * SIMDD;
        p0z = gz + dj * SIMDD;
        p1x = gx - dj * SIMDD;
        p1y = gy - dj * SIMDD;
        p1z = gz - dj * SIMDD;
        for (j = 1; j < nmax; j++) {
                ptr = j * dj;
                for (n = 0; n < nrys_roots; n++) {
                        r1 = MM_LOAD(t2+n*SIMDD) * MM_SET1(j);
MM_STORE(p0x+(ptr+n)*SIMDD, r1 * MM_LOAD(p1x+(ptr+n)*SIMDD) + MM_LOAD(r0jx+n*SIMDD) * MM_LOAD(gx+(ptr+n)*SIMDD));
MM_STORE(p0y+(ptr+n)*SIMDD, r1 * MM_LOAD(p1y+(ptr+n)*SIMDD) + MM_LOAD(r0jy+n*SIMDD) * MM_LOAD(gy+(ptr+n)*SIMDD));
MM_STORE(p0z+(ptr+n)*SIMDD, r1 * MM_LOAD(p1z+(ptr+n)*SIMDD) + MM_LOAD(r0jz+n*SIMDD) * MM_LOAD(gz+(ptr+n)*SIMDD));
                }
        }

        r0 = MM_SET1(rirj[0]);
        r1 = MM_SET1(rirj[1]);
        r2 = MM_SET1(rirj[2]);
        
        for (i = 1; i <= li; i++) {
        for (j = 0; j <= nmax-i; j++) { // upper limit lj+lk
                off = i * di + j * dj;
                p0x = gx + off * SIMDD;
                p0y = gy + off * SIMDD;
                p0z = gz + off * SIMDD;
                p2x = p0x - di * SIMDD;
                p2y = p0y - di * SIMDD;
                p2z = p0z - di * SIMDD;
                p1x = p2x + dj * SIMDD;
                p1y = p2y + dj * SIMDD;
                p1z = p2z + dj * SIMDD;
                for (n = 0; n < nrys_roots; n++) {
                        MM_STORE(p0x+n*SIMDD, MM_LOAD(p1x+n*SIMDD) - r0 * MM_LOAD(p2x+n*SIMDD));
                        MM_STORE(p0y+n*SIMDD, MM_LOAD(p1y+n*SIMDD) - r1 * MM_LOAD(p2y+n*SIMDD));
                        MM_STORE(p0z+n*SIMDD, MM_LOAD(p1z+n*SIMDD) - r2 * MM_LOAD(p2z+n*SIMDD));
                }
        } }

        r0 = MM_SET1(rjrk[0]);
        r1 = MM_SET1(rjrk[1]);
        r2 = MM_SET1(rjrk[2]);
        for (k = 1; k <= lk; k++) {
        for (j = 0; j <= mmax-k; j++) {
                off = k * dk + j * dj;
                p0x = gx + off * SIMDD;
                p0y = gy + off * SIMDD;
                p0z = gz + off * SIMDD;
                p2x = p0x - dk * SIMDD;
                p2y = p0y - dk * SIMDD;
                p2z = p0z - dk * SIMDD;
                p1x = p2x + dj * SIMDD;
                p1y = p2y + dj * SIMDD;
                p1z = p2z + dj * SIMDD;
                for (i = 0; i <= li; i++) {
                for (n = i*di; n < i*di+nrys_roots; n++) {
                        MM_STORE(p0x+n*SIMDD, MM_LOAD(p1x+n*SIMDD) + r0 * MM_LOAD(p2x+n*SIMDD));
                        MM_STORE(p0y+n*SIMDD, MM_LOAD(p1y+n*SIMDD) + r1 * MM_LOAD(p2y+n*SIMDD));
                        MM_STORE(p0z+n*SIMDD, MM_LOAD(p1z+n*SIMDD) + r2 * MM_LOAD(p2z+n*SIMDD));
                } }
        } }
}


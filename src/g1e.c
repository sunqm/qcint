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
#include <assert.h>
#include "cint_bas.h"
#include "simd.h"
#include "misc.h"
#include "g1e.h"
#include "rys_roots.h"

#define DEF_GXYZ(type, G, GX, GY, GZ) \
        type *RESTRICT GX = G; \
        type *RESTRICT GY = G + envs->g_size     * SIMDD; \
        type *RESTRICT GZ = G + envs->g_size * 2 * SIMDD

void CINTinit_int1e_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
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
        envs->i_l = bas(ANG_OF, i_sh);
        envs->j_l = bas(ANG_OF, j_sh);
        envs->x_ctr[0] = bas(NCTR_OF, i_sh);
        envs->x_ctr[1] = bas(NCTR_OF, j_sh);
        envs->nfi = (envs->i_l+1)*(envs->i_l+2)/2;
        envs->nfj = (envs->j_l+1)*(envs->j_l+2)/2;
        envs->nf = envs->nfi * envs->nfj;
        envs->common_factor = 1;
        if (env[PTR_EXPCUTOFF] == 0) {
                envs->expcutoff = EXPCUTOFF;
        } else {
                envs->expcutoff = MAX(MIN_EXPCUTOFF, env[PTR_EXPCUTOFF]);
        }

        envs->li_ceil = envs->i_l + ng[IINC];
        envs->lj_ceil = envs->j_l + ng[JINC];
        envs->ri = env + atm(PTR_COORD, bas(ATOM_OF, i_sh));
        envs->rj = env + atm(PTR_COORD, bas(ATOM_OF, j_sh));

        envs->gbits = ng[GSHIFT];
        envs->ncomp_e1 = ng[POS_E1];
        envs->ncomp_tensor = ng[TENSOR];
        if (ng[SLOT_RYS_ROOTS] > 0) {
                envs->nrys_roots = ng[SLOT_RYS_ROOTS];
        } else {
                envs->nrys_roots = (envs->li_ceil + envs->lj_ceil)/2 + 1;
        }

        int dli, dlj;
        int ibase = envs->li_ceil > envs->lj_ceil;
        if (ibase) {
                dli = envs->li_ceil + envs->lj_ceil + 1;
                dlj = envs->lj_ceil + 1;
                envs->rirj[0] = envs->ri[0] - envs->rj[0];
                envs->rirj[1] = envs->ri[1] - envs->rj[1];
                envs->rirj[2] = envs->ri[2] - envs->rj[2];
        } else {
                dli = envs->li_ceil + 1;
                dlj = envs->li_ceil + envs->lj_ceil + 1;
                envs->rirj[0] = envs->rj[0] - envs->ri[0];
                envs->rirj[1] = envs->rj[1] - envs->ri[1];
                envs->rirj[2] = envs->rj[2] - envs->ri[2];
        }
        envs->g_stride_i = envs->nrys_roots;
        envs->g_stride_j = envs->nrys_roots * dli;
        envs->g_size     = envs->nrys_roots * dli * dlj;

        envs->lk_ceil = 1;
        envs->ll_ceil = 1;
        envs->g_stride_k = 0;
        envs->g_stride_l = 0;
}

void CINTg2c_index_xyz(int *idx, CINTEnvVars *envs)
{
        int i_l = envs->i_l;
        int j_l = envs->j_l;
        int nfi = envs->nfi;
        int nfj = envs->nfj;
        int di = envs->g_stride_i;
        int dj = envs->g_stride_j;
        int i, j, n;
        int ofx, ofjx;
        int ofy, ofjy;
        int ofz, ofjz;
        int i_nx[CART_MAX], i_ny[CART_MAX], i_nz[CART_MAX];
        int j_nx[CART_MAX], j_ny[CART_MAX], j_nz[CART_MAX];

        CINTcart_comp(i_nx, i_ny, i_nz, i_l);
        CINTcart_comp(j_nx, j_ny, j_nz, j_l);

        ofx = 0;
        ofy = envs->g_size;
        ofz = envs->g_size * 2;
        n = 0;
        for (j = 0; j < nfj; j++) {
                ofjx = ofx + dj * j_nx[j];
                ofjy = ofy + dj * j_ny[j];
                ofjz = ofz + dj * j_nz[j];
                switch (i_l) {
                case 0:
                        idx[n+0] = ofjx;
                        idx[n+1] = ofjy;
                        idx[n+2] = ofjz;
                        n += 3;
                        break;
                case 1:
                        idx[n+0] = ofjx + di;
                        idx[n+1] = ofjy;
                        idx[n+2] = ofjz;
                        idx[n+3] = ofjx;
                        idx[n+4] = ofjy + di;
                        idx[n+5] = ofjz;
                        idx[n+6] = ofjx;
                        idx[n+7] = ofjy;
                        idx[n+8] = ofjz + di;
                        n += 9;
                        break;
                case 2:
                        idx[n+0 ] = ofjx + di*2;
                        idx[n+1 ] = ofjy;
                        idx[n+2 ] = ofjz;
                        idx[n+3 ] = ofjx + di;
                        idx[n+4 ] = ofjy + di;
                        idx[n+5 ] = ofjz;
                        idx[n+6 ] = ofjx + di;
                        idx[n+7 ] = ofjy;
                        idx[n+8 ] = ofjz + di;
                        idx[n+9 ] = ofjx;
                        idx[n+10] = ofjy + di*2;
                        idx[n+11] = ofjz;
                        idx[n+12] = ofjx;
                        idx[n+13] = ofjy + di;
                        idx[n+14] = ofjz + di;
                        idx[n+15] = ofjx;
                        idx[n+16] = ofjy;
                        idx[n+17] = ofjz + di*2;
                        n += 18;
                        break;
                default:
                        for (i = 0; i < nfi; i++) {
                                idx[n+0] = ofjx + di * i_nx[i];
                                idx[n+1] = ofjy + di * i_ny[i];
                                idx[n+2] = ofjz + di * i_nz[i];
                                n += 3;
                        }
                }
        }
}
void CINTg1e_index_xyz(int *idx, CINTEnvVars *envs)
{
        CINTg2c_index_xyz(idx, envs);
}


void CINTg1e_ovlp(double *g, CINTEnvVars *envs, int count)
{
        double *gx = g;
        double *gy = g + envs->g_size     * SIMDD;
        double *gz = g + envs->g_size * 2 * SIMDD;
        double *p0x, *p0y, *p0z, *p1x, *p1y, *p1z;
        __MD r1, aij;

        aij = MM_LOAD(envs->ai) + MM_LOAD(envs->aj);
        MM_STORE(gx, MM_SET1(1.));
        MM_STORE(gy, MM_SET1(1.));
        MM_STORE(gz, MM_LOAD(envs->fac) * MM_SET1(SQRTPI*M_PI) / (aij * MM_SQRT(aij)));

        int nmax = envs->li_ceil + envs->lj_ceil;
        if (nmax == 0) {
                return;
        }

        double *RESTRICT rij = envs->rij;
        double *RESTRICT rirj = envs->rirj;
        int lj, di, dj;
        int i, j, n, ptr;
        double *rx;
        if (envs->li_ceil > envs->lj_ceil) {
                // li = envs->li_ceil;
                lj = envs->lj_ceil;
                di = envs->g_stride_i;
                dj = envs->g_stride_j;
                rx = envs->ri;
        } else {
                // li = envs->lj_ceil;
                lj = envs->li_ceil;
                di = envs->g_stride_j;
                dj = envs->g_stride_i;
                rx = envs->rj;
        }
        __MD rijrx[3];
        rijrx[0] = MM_LOAD(rij+0*SIMDD) - MM_SET1(rx[0]);
        rijrx[1] = MM_LOAD(rij+1*SIMDD) - MM_SET1(rx[1]);
        rijrx[2] = MM_LOAD(rij+2*SIMDD) - MM_SET1(rx[2]);

        //gx[1] = rijrx[0] * gx[0];
        //gy[1] = rijrx[1] * gy[0];
        //gz[1] = rijrx[2] * gz[0];
        MM_STORE(gx+di*SIMDD, rijrx[0] * MM_LOAD(gx));
        MM_STORE(gy+di*SIMDD, rijrx[1] * MM_LOAD(gy));
        MM_STORE(gz+di*SIMDD, rijrx[2] * MM_LOAD(gz));

        p0x = gx + di * SIMDD;
        p0y = gy + di * SIMDD;
        p0z = gz + di * SIMDD;
        p1x = gx - di * SIMDD;
        p1y = gy - di * SIMDD;
        p1z = gz - di * SIMDD;
        r1 = MM_SET1(.5) / aij;
        for (i = 1; i < nmax; i++) {
MM_STORE(p0x+i*di*SIMDD, r1 * MM_SET1(i) * MM_LOAD(p1x+i*di*SIMDD) + rijrx[0] * MM_LOAD(gx+i*di*SIMDD));
MM_STORE(p0y+i*di*SIMDD, r1 * MM_SET1(i) * MM_LOAD(p1y+i*di*SIMDD) + rijrx[1] * MM_LOAD(gy+i*di*SIMDD));
MM_STORE(p0z+i*di*SIMDD, r1 * MM_SET1(i) * MM_LOAD(p1z+i*di*SIMDD) + rijrx[2] * MM_LOAD(gz+i*di*SIMDD));
        }

        p0x = gx  - dj * SIMDD;
        p0y = gy  - dj * SIMDD;
        p0z = gz  - dj * SIMDD;
        p1x = p0x + di * SIMDD;
        p1y = p0y + di * SIMDD;
        p1z = p0z + di * SIMDD;
        for (j = 1; j <= lj; j++) {
                ptr = dj * j;
                for (i = 0, n = ptr; i <= nmax-j; i++, n+=di) {
MM_STORE(gx+n*SIMDD, MM_LOAD(p1x+n*SIMDD) + MM_SET1(rirj[0]) * MM_LOAD(p0x+n*SIMDD));
MM_STORE(gy+n*SIMDD, MM_LOAD(p1y+n*SIMDD) + MM_SET1(rirj[1]) * MM_LOAD(p0y+n*SIMDD));
MM_STORE(gz+n*SIMDD, MM_LOAD(p1z+n*SIMDD) + MM_SET1(rirj[2]) * MM_LOAD(p0z+n*SIMDD));
                }
        }
}

/*
 * Calculate temporary parameter tau for nuclear charge distribution.
 * The charge parameter zeta is defined as    rho(r) = Norm * exp(-zeta*r^2)
 */
double CINTnuc_mod(double aij, int nuc_id, int *atm, double *env)
{
        double zeta;
        if (nuc_id < 0) {
                zeta = env[PTR_RINV_ZETA];
        } else if (atm(NUC_MOD_OF, nuc_id) == GAUSSIAN_NUC) {
                zeta = env[atm(PTR_ZETA, nuc_id)];
        } else {
                zeta = 0;
        }

        if (zeta > 0) {
                return sqrt(zeta / (aij + zeta));
        } else {
                return 1;
        }
}

void CINTg1e_nuc(double *g, CINTEnvVars *envs, int count, int nuc_id)
{
        int nrys_roots = envs->nrys_roots;
        int *atm = envs->atm;
        double *env = envs->env;
        double *RESTRICT rij = envs->rij;
        double *RESTRICT gx = g;
        double *RESTRICT gy = g + envs->g_size     * SIMDD;
        double *RESTRICT gz = g + envs->g_size * 2 * SIMDD;
        ALIGNMM double tau[SIMDD];
        ALIGNMM double x[SIMDD];
        ALIGNMM double u[MXRYSROOTS*SIMDD];
        double *RESTRICT w = gz;
        double *cr;
        int i, j, k, n;
        __MD crij[3];
        __MD r0, r1, r2, fac1, aij;

        aij = MM_LOAD(envs->ai) + MM_LOAD(envs->aj);
        MM_STORE(tau, aij);
        for (k = 0; k < count; k++) {
                tau[k] = CINTnuc_mod(tau[k], nuc_id, atm, env);
        }

        if (nuc_id < 0) {
                fac1 = MM_SET1(2*M_PI) * MM_LOAD(envs->fac) * MM_LOAD(tau) / aij;
                cr = env + PTR_RINV_ORIG;
        } else if (atm(NUC_MOD_OF,nuc_id) == FRAC_CHARGE_NUC) {
                fac1 = MM_SET1(2*M_PI) * MM_SET1(-env[atm[PTR_FRAC_CHARGE+nuc_id*ATM_SLOTS]]);
                fac1 = fac1 * MM_LOAD(envs->fac) * MM_LOAD(tau) / aij;
                cr = env + atm(PTR_COORD, nuc_id);
        } else {
                fac1 = MM_SET1(2*M_PI) * MM_SET1(-fabs(atm[CHARGE_OF+nuc_id*ATM_SLOTS]));
                fac1 = fac1 * MM_LOAD(envs->fac) * MM_LOAD(tau) / aij;
                cr = env + atm(PTR_COORD, nuc_id);
        }
        crij[0] = MM_SET1(cr[0]) - MM_LOAD(rij+0*SIMDD);
        crij[1] = MM_SET1(cr[1]) - MM_LOAD(rij+1*SIMDD);
        crij[2] = MM_SET1(cr[2]) - MM_LOAD(rij+2*SIMDD);
        MM_STORE(x, aij * MM_LOAD(tau) * MM_LOAD(tau) * SQUARE(crij));
        _CINTrys_roots_batch(nrys_roots, x, u, w, count);

        r1 = MM_SET1(1.);
        for (i = 0; i < nrys_roots; i++) {
                MM_STORE(gx+i*SIMDD, r1);
                MM_STORE(gy+i*SIMDD, r1);
                MM_STORE(gz+i*SIMDD, fac1 * MM_LOAD(gz+i*SIMDD));
        }
        int nmax = envs->li_ceil + envs->lj_ceil;
        if (nmax == 0) {
                return;
        }

        double *RESTRICT p0x, *RESTRICT p0y, *RESTRICT p0z;
        double *RESTRICT p1x, *RESTRICT p1y, *RESTRICT p1z;
        double *RESTRICT p2x, *RESTRICT p2y, *RESTRICT p2z;
        int lj, di, dj;
        double *rx;
        if (envs->li_ceil > envs->lj_ceil) {
                // li = envs->li_ceil;
                lj = envs->lj_ceil;
                di = envs->g_stride_i;
                dj = envs->g_stride_j;
                rx = envs->ri;
        } else {
                // li = envs->lj_ceil;
                lj = envs->li_ceil;
                di = envs->g_stride_j;
                dj = envs->g_stride_i;
                rx = envs->rj;
        }
        __MD rijrx = MM_LOAD(rij+0*SIMDD) - MM_SET1(rx[0]);
        __MD rijry = MM_LOAD(rij+1*SIMDD) - MM_SET1(rx[1]);
        __MD rijrz = MM_LOAD(rij+2*SIMDD) - MM_SET1(rx[2]);
        __MD aij2 = MM_SET1(0.5) / aij;
        __MD ru, rt;
        __MD r0x, r0y, r0z;
        __MD r1x, r1y, r1z;
        __MD r2x, r2y, r2z;

        p0x = gx + di * SIMDD;
        p0y = gy + di * SIMDD;
        p0z = gz + di * SIMDD;
        for (n = 0; n < nrys_roots; n++) {
                ru = MM_DIV(MM_LOAD(tau) * MM_LOAD(tau) * MM_LOAD(u+n*SIMDD), MM_SET1(1.) + MM_LOAD(u+n*SIMDD));
                rt = aij2 - aij2 * ru;
                r0 = rijrx + ru * crij[0];
                r1 = rijry + ru * crij[1];
                r2 = rijrz + ru * crij[2];

                r2x = MM_LOAD(gx+n*SIMDD);
                r2y = MM_LOAD(gy+n*SIMDD);
                r2z = MM_LOAD(gz+n*SIMDD);
                r1x = r0 * r2x;
                r1y = r1 * r2y;
                r1z = r2 * r2z;
                MM_STORE(p0x+n*SIMDD, r1x);
                MM_STORE(p0y+n*SIMDD, r1y);
                MM_STORE(p0z+n*SIMDD, r1z);
                for (i = 1; i < nmax; i++) {
                        r0x = MM_SET1(i) * rt * r2x + r0 * r1x;
                        r0y = MM_SET1(i) * rt * r2y + r1 * r1y;
                        r0z = MM_SET1(i) * rt * r2z + r2 * r1z;
                        MM_STORE(p0x+(n+i*di)*SIMDD, r0x);
                        MM_STORE(p0y+(n+i*di)*SIMDD, r0y);
                        MM_STORE(p0z+(n+i*di)*SIMDD, r0z);
                        r2x = r1x;
                        r2y = r1y;
                        r2z = r1z;
                        r1x = r0x;
                        r1y = r0y;
                        r1z = r0z;
                }
        }

        __MD rirjx = MM_SET1(envs->rirj[0]);
        __MD rirjy = MM_SET1(envs->rirj[1]);
        __MD rirjz = MM_SET1(envs->rirj[2]);
        for (j = 1; j <= lj; j++) {
        for (n = 0; n < nrys_roots; n++) {
                p0x = gx + (n + j * dj) * SIMDD;
                p0y = gy + (n + j * dj) * SIMDD;
                p0z = gz + (n + j * dj) * SIMDD;
                p1x = p0x - dj * SIMDD;
                p1y = p0y - dj * SIMDD;
                p1z = p0z - dj * SIMDD;
                p2x = p1x + di * SIMDD;
                p2y = p1y + di * SIMDD;
                p2z = p1z + di * SIMDD;
                r1x = MM_LOAD(p1x);
                r1y = MM_LOAD(p1y);
                r1z = MM_LOAD(p1z);
                for (i = 0; i <= nmax - j; i++) {
                        r2x = MM_LOAD(p2x+i*di*SIMDD);
                        r2y = MM_LOAD(p2y+i*di*SIMDD);
                        r2z = MM_LOAD(p2z+i*di*SIMDD);
                        MM_STORE(p0x+i*di*SIMDD, r2x + rirjx * r1x);
                        MM_STORE(p0y+i*di*SIMDD, r2y + rirjy * r1y);
                        MM_STORE(p0z+i*di*SIMDD, r2z + rirjz * r1z);
                        r1x = r2x;
                        r1y = r2y;
                        r1z = r2z;
                }
        } }
}

void CINTnabla1i_1e(double *f, double *g,
                    int li, int lj, int lk, CINTEnvVars *envs)
{
        int i, j, k, n, ptr;
        int di = envs->g_stride_i;
        int dj = envs->g_stride_j;
        int dk = envs->g_stride_k;
        DEF_GXYZ(double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);
        double *RESTRICT p1x = gx - di * SIMDD;
        double *RESTRICT p1y = gy - di * SIMDD;
        double *RESTRICT p1z = gz - di * SIMDD;
        double *RESTRICT p2x = gx + di * SIMDD;
        double *RESTRICT p2y = gy + di * SIMDD;
        double *RESTRICT p2z = gz + di * SIMDD;
        __MD ai2 = MM_MUL(MM_SET1(-2.), MM_LOAD(envs->ai));
        __MD ri;

        for (j = 0; j <= lj; j++)
        for (k = 0; k <= lk; k++) {
                ptr = dj * j + dk * k;
                //f(...,0,...) = -2*ai*g(...,1,...)
MM_STORE(fx+ptr*SIMDD, MM_MUL(ai2, MM_LOAD(p2x+ptr*SIMDD)));
MM_STORE(fy+ptr*SIMDD, MM_MUL(ai2, MM_LOAD(p2y+ptr*SIMDD)));
MM_STORE(fz+ptr*SIMDD, MM_MUL(ai2, MM_LOAD(p2z+ptr*SIMDD)));
                //f(...,i,...) = i*g(...,i-1,...)-2*ai*g(...,i+1,...)
                ri = MM_SET1(1.);
                for (i = 1; i <= li; i++) {
                        n = ptr + di * i;
MM_STORE(fx+n*SIMDD, ri * MM_LOAD(p1x+n*SIMDD) + ai2 * MM_LOAD(p2x+n*SIMDD));
MM_STORE(fy+n*SIMDD, ri * MM_LOAD(p1y+n*SIMDD) + ai2 * MM_LOAD(p2y+n*SIMDD));
MM_STORE(fz+n*SIMDD, ri * MM_LOAD(p1z+n*SIMDD) + ai2 * MM_LOAD(p2z+n*SIMDD));
                        ri = MM_ADD(ri, MM_SET1(1.));
                }
        }
}

void CINTnabla1j_1e(double *f, double *g,
                    int li, int lj, int lk, CINTEnvVars *envs)
{
        int i, j, k, n, ptr;
        int di = envs->g_stride_i;
        int dj = envs->g_stride_j;
        int dk = envs->g_stride_k;
        DEF_GXYZ(double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);
        double *RESTRICT p1x = gx - dj * SIMDD;
        double *RESTRICT p1y = gy - dj * SIMDD;
        double *RESTRICT p1z = gz - dj * SIMDD;
        double *RESTRICT p2x = gx + dj * SIMDD;
        double *RESTRICT p2y = gy + dj * SIMDD;
        double *RESTRICT p2z = gz + dj * SIMDD;
        __MD aj2 = MM_MUL(MM_SET1(-2.), MM_LOAD(envs->aj));
        __MD rj;

        //f(...,0,...) = -2*aj*g(...,1,...)
        for (k = 0; k <= lk; k++) {
                ptr = dk * k;
                for (i = 0; i <= li; i++) {
                        n = ptr + di * i;
MM_STORE(fx+n*SIMDD, MM_MUL(aj2, MM_LOAD(p2x+n*SIMDD)));
MM_STORE(fy+n*SIMDD, MM_MUL(aj2, MM_LOAD(p2y+n*SIMDD)));
MM_STORE(fz+n*SIMDD, MM_MUL(aj2, MM_LOAD(p2z+n*SIMDD)));
                }
        }
        //f(...,j,...) = j*g(...,j-1,...)-2*aj*g(...,j+1,...)
        for (j = 1; j <= lj; j++) {
                rj = MM_SET1(j);
                for (k = 0; k <= lk; k++) {
                        ptr = dj * j + dk * k;
                        for (i = 0; i <= li; i++) {
                                n = ptr + di * i;
MM_STORE(fx+n*SIMDD, rj * MM_LOAD(p1x+n*SIMDD) + aj2 * MM_LOAD(p2x+n*SIMDD));
MM_STORE(fy+n*SIMDD, rj * MM_LOAD(p1y+n*SIMDD) + aj2 * MM_LOAD(p2y+n*SIMDD));
MM_STORE(fz+n*SIMDD, rj * MM_LOAD(p1z+n*SIMDD) + aj2 * MM_LOAD(p2z+n*SIMDD));
                        }
                }
        }
}

void CINTnabla1k_1e(double *f, double *g,
                    int li, int lj, int lk, CINTEnvVars *envs)
{
        int i, j, k, n, ptr;
        int di = envs->g_stride_i;
        int dj = envs->g_stride_j;
        int dk = envs->g_stride_k;
        DEF_GXYZ(double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);
        double *RESTRICT p1x = gx - dk * SIMDD;
        double *RESTRICT p1y = gy - dk * SIMDD;
        double *RESTRICT p1z = gz - dk * SIMDD;
        double *RESTRICT p2x = gx + dk * SIMDD;
        double *RESTRICT p2y = gy + dk * SIMDD;
        double *RESTRICT p2z = gz + dk * SIMDD;
        __MD ak2 = MM_MUL(MM_SET1(-2.), MM_LOAD(envs->ak));
        __MD rk;

        for (j = 0; j <= lj; j++) {
                ptr = dj * j;
                //f(...,0,...) = -2*ak*g(...,1,...)
                for (i = 0; i <= li; i++) {
                        n = ptr + di * i;
MM_STORE(fx+n*SIMDD, MM_MUL(ak2, MM_LOAD(p2x+n*SIMDD)));
MM_STORE(fy+n*SIMDD, MM_MUL(ak2, MM_LOAD(p2y+n*SIMDD)));
MM_STORE(fz+n*SIMDD, MM_MUL(ak2, MM_LOAD(p2z+n*SIMDD)));
                }
                //f(...,k,...) = k*g(...,k-1,...)-2*ak*g(...,k+1,...)
                rk = MM_SET1(1.);
                for (k = 1; k <= lk; k++) {
                        ptr = dj * j + dk * k;
                        for (i = 0; i <= li; i++) {
                                n = ptr + di * i;
MM_STORE(fx+n*SIMDD, rk * MM_LOAD(p1x+n*SIMDD) + ak2 * MM_LOAD(p2x+n*SIMDD));
MM_STORE(fy+n*SIMDD, rk * MM_LOAD(p1y+n*SIMDD) + ak2 * MM_LOAD(p2y+n*SIMDD));
MM_STORE(fz+n*SIMDD, rk * MM_LOAD(p1z+n*SIMDD) + ak2 * MM_LOAD(p2z+n*SIMDD));
                        }
                        rk = MM_ADD(rk, MM_SET1(1.));
                }
        }
}


/*
 * < x^1 i | j >
 * ri is the shift from the center R_O to the center of |i>
 * r - R_O = (r-R_i) + ri, ri = R_i - R_O
 */
void CINTx1i_1e(double *f, double *g, double *ri,
                int li, int lj, int lk, CINTEnvVars *envs)
{
        int i, j, k, n, ptr;
        int di = envs->g_stride_i;
        int dj = envs->g_stride_j;
        int dk = envs->g_stride_k;
        DEF_GXYZ(double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);
        double *RESTRICT p1x = gx + di * SIMDD;
        double *RESTRICT p1y = gy + di * SIMDD;
        double *RESTRICT p1z = gz + di * SIMDD;
        __MD r0 = MM_SET1(ri[0]);
        __MD r1 = MM_SET1(ri[1]);
        __MD r2 = MM_SET1(ri[2]);

        for (j = 0; j <= lj; j++)
        for (k = 0; k <= lk; k++) {
                //f(...,0:li,...) = g(...,1:li+1,...) + ri(1)*g(...,0:li,...)
                ptr = dj * j + dk * k;
                for (i = 0; i <= li; i++) {
                        n = ptr + di * i;
MM_STORE(fx+n*SIMDD, MM_FMA(r0, MM_LOAD(gx+n*SIMDD), MM_LOAD(p1x+n*SIMDD)));
MM_STORE(fy+n*SIMDD, MM_FMA(r1, MM_LOAD(gy+n*SIMDD), MM_LOAD(p1y+n*SIMDD)));
MM_STORE(fz+n*SIMDD, MM_FMA(r2, MM_LOAD(gz+n*SIMDD), MM_LOAD(p1z+n*SIMDD)));
                }
        }
}

void CINTx1j_1e(double *f, double *g, double *rj,
                int li, int lj, int lk, CINTEnvVars *envs)
{
        int i, j, k, n, ptr;
        int di = envs->g_stride_i;
        int dj = envs->g_stride_j;
        int dk = envs->g_stride_k;
        DEF_GXYZ(double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);
        double *RESTRICT p1x = gx + dj * SIMDD;
        double *RESTRICT p1y = gy + dj * SIMDD;
        double *RESTRICT p1z = gz + dj * SIMDD;
        __MD r0 = MM_SET1(rj[0]);
        __MD r1 = MM_SET1(rj[1]);
        __MD r2 = MM_SET1(rj[2]);

        for (j = 0; j <= lj; j++)
        for (k = 0; k <= lk; k++) {
                // f(...,0:lj,...) = g(...,1:lj+1,...) + rj(1)*g(...,0:lj,...)
                ptr = dj * j + dk * k;
                for (i = 0; i <= li; i++) {
                        n = ptr + di * i;
MM_STORE(fx+n*SIMDD, MM_FMA(r0, MM_LOAD(gx+n*SIMDD), MM_LOAD(p1x+n*SIMDD)));
MM_STORE(fy+n*SIMDD, MM_FMA(r1, MM_LOAD(gy+n*SIMDD), MM_LOAD(p1y+n*SIMDD)));
MM_STORE(fz+n*SIMDD, MM_FMA(r2, MM_LOAD(gz+n*SIMDD), MM_LOAD(p1z+n*SIMDD)));
                }
        }
}

void CINTx1k_1e(double *f, double *g, double *rk,
                int li, int lj, int lk, CINTEnvVars *envs)
{
        int i, j, k, n, ptr;
        int di = envs->g_stride_i;
        int dj = envs->g_stride_j;
        int dk = envs->g_stride_k;
        DEF_GXYZ(double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);
        double *RESTRICT p1x = gx + dk * SIMDD;
        double *RESTRICT p1y = gy + dk * SIMDD;
        double *RESTRICT p1z = gz + dk * SIMDD;
        __MD r0 = MM_SET1(rk[0]);
        __MD r1 = MM_SET1(rk[1]);
        __MD r2 = MM_SET1(rk[2]);

        for (j = 0; j <= lj; j++)
        for (k = 0; k <= lk; k++) {
                // f(...,0:lk,...) = g(...,1:lk+1,...) + rk(1)*g(...,0:lk,...)
                ptr = dj * j + dk * k;
                for (i = 0; i <= li; i++) {
                        n = ptr + di * i;
MM_STORE(fx+n*SIMDD, MM_FMA(r0, MM_LOAD(gx+n*SIMDD), MM_LOAD(p1x+n*SIMDD)));
MM_STORE(fy+n*SIMDD, MM_FMA(r1, MM_LOAD(gy+n*SIMDD), MM_LOAD(p1y+n*SIMDD)));
MM_STORE(fz+n*SIMDD, MM_FMA(r2, MM_LOAD(gz+n*SIMDD), MM_LOAD(p1z+n*SIMDD)));
                }
        }
}


void CINTprim_to_ctr_0(double *RESTRICT gc, double *RESTRICT gp, double *RESTRICT coeff, size_t nf,
                       int nprim, int nctr, int non0ctr, int *sortedidx)
{
        int i;
        size_t n;
        double c0, c1, c2, c3;
        double *RESTRICT p0;
        double *RESTRICT p1;
        double *RESTRICT p2;
        double *RESTRICT p3;
        __MD r0, r1, r2, r3, rg;

        switch (non0ctr) {
        case 1:
                c0 = coeff[nprim*sortedidx[0]];
                p0 = gc + nf * sortedidx[0];
                r0 = MM_SET1(c0);
                for (n = 0; n < nf/SIMDD; n++) {
                        MM_STOREU(p0+n*SIMDD, MM_MUL(r0, MM_LOADU(gp+n*SIMDD)));
                }
                for (n = n*SIMDD; n < nf; n++) {
                        p0[n] = c0 * gp[n];
                }
                break;
        case 2:
                c0 = coeff[nprim*sortedidx[0]];
                c1 = coeff[nprim*sortedidx[1]];
                p0 = gc + nf * sortedidx[0];
                p1 = gc + nf * sortedidx[1];
                r0 = MM_SET1(c0);
                r1 = MM_SET1(c1);
                for (n = 0; n < nf/SIMDD; n++) {
                        rg = MM_LOADU(gp+n*SIMDD);
                        MM_STOREU(p0+n*SIMDD, MM_MUL(r0, rg));
                        MM_STOREU(p1+n*SIMDD, MM_MUL(r1, rg));
                }
                for (n = n*SIMDD; n < nf; n++) {
                        p0[n] = c0 * gp[n];
                        p1[n] = c1 * gp[n];
                }
                break;
        case 3:
                c0 = coeff[nprim*sortedidx[0]];
                c1 = coeff[nprim*sortedidx[1]];
                c2 = coeff[nprim*sortedidx[2]];
                p0 = gc + nf * sortedidx[0];
                p1 = gc + nf * sortedidx[1];
                p2 = gc + nf * sortedidx[2];
                r0 = MM_SET1(c0);
                r1 = MM_SET1(c1);
                r2 = MM_SET1(c2);
                for (n = 0; n < nf/SIMDD; n++) {
                        rg = MM_LOADU(gp+n*SIMDD);
                        MM_STOREU(p0+n*SIMDD, MM_MUL(r0, rg));
                        MM_STOREU(p1+n*SIMDD, MM_MUL(r1, rg));
                        MM_STOREU(p2+n*SIMDD, MM_MUL(r2, rg));
                }
                for (n = n*SIMDD; n < nf; n++) {
                        p0[n] = c0 * gp[n];
                        p1[n] = c1 * gp[n];
                        p2[n] = c2 * gp[n];
                }
                break;
        case 4:
                c0 = coeff[nprim*sortedidx[0]];
                c1 = coeff[nprim*sortedidx[1]];
                c2 = coeff[nprim*sortedidx[2]];
                c3 = coeff[nprim*sortedidx[3]];
                p0 = gc + nf * sortedidx[0];
                p1 = gc + nf * sortedidx[1];
                p2 = gc + nf * sortedidx[2];
                p3 = gc + nf * sortedidx[3];
                r0 = MM_SET1(c0);
                r1 = MM_SET1(c1);
                r2 = MM_SET1(c2);
                r3 = MM_SET1(c3);
                for (n = 0; n < nf/SIMDD; n++) {
                        rg = MM_LOADU(gp+n*SIMDD);
                        MM_STOREU(p0+n*SIMDD, MM_MUL(r0, rg));
                        MM_STOREU(p1+n*SIMDD, MM_MUL(r1, rg));
                        MM_STOREU(p2+n*SIMDD, MM_MUL(r2, rg));
                        MM_STOREU(p3+n*SIMDD, MM_MUL(r3, rg));
                }
                for (n = n*SIMDD; n < nf; n++) {
                        p0[n] = c0 * gp[n];
                        p1[n] = c1 * gp[n];
                        p2[n] = c2 * gp[n];
                        p3[n] = c3 * gp[n];
                }
                break;
        default:
                for (i = 0; i < non0ctr-1; i+=2) {
                        c0 = coeff[nprim*sortedidx[i]];
                        c1 = coeff[nprim*sortedidx[i+1]];
                        p0 = gc + nf * sortedidx[i];
                        p1 = gc + nf * sortedidx[i+1];
                        r0 = MM_SET1(c0);
                        r1 = MM_SET1(c1);
                        for (n = 0; n < nf/SIMDD; n++) {
                                rg = MM_LOADU(gp+n*SIMDD);
                                MM_STOREU(p0+n*SIMDD, MM_MUL(r0, rg));
                                MM_STOREU(p1+n*SIMDD, MM_MUL(r1, rg));
                        }
                        for (n = n*SIMDD; n < nf; n++) {
                                p0[n] = c0 * gp[n];
                                p1[n] = c1 * gp[n];
                        }
                }
                if (i < non0ctr) {
                        c0 = coeff[nprim*sortedidx[i]];
                        p0 = gc + nf * sortedidx[i];
                        r0 = MM_SET1(c0);
                        for (n = 0; n < nf/SIMDD; n++) {
                                MM_STOREU(p0+n*SIMDD, MM_MUL(r0, MM_LOADU(gp+n*SIMDD)));
                        }
                        for (n = n*SIMDD; n < nf; n++) {
                                p0[n] = c0 * gp[n];
                        }
                }
        }
        r0 = MM_SET1(0.);
        for (i = non0ctr; i < nctr; i++) {
                p0 = gc + nf * sortedidx[i];
                for (n = 0; n < nf/SIMDD; n++) {
                        MM_STOREU(p0+n*SIMDD, r0);
                }
                for (n = n*SIMDD; n < nf; n++) {
                        p0[n] = 0;
                }
        }
}

void CINTprim_to_ctr_1(double *RESTRICT gc, double *RESTRICT gp, double *RESTRICT coeff, size_t nf,
                       int nprim, int nctr, int non0ctr, int *sortedidx)
{
        int i;
        size_t n;
        double c0, c1, c2, c3;
        double *RESTRICT p0;
        double *RESTRICT p1;
        double *RESTRICT p2;
        double *RESTRICT p3;
        __MD r0, r1, r2, r3, rg;

        switch (non0ctr) {
        case 1:
                c0 = coeff[nprim*sortedidx[0]];
                p0 = gc + nf * sortedidx[0];
                r0 = MM_SET1(c0);
                for (n = 0; n < nf/SIMDD; n++) {
                        rg = MM_LOADU(gp+n*SIMDD);
                        MM_STOREU(p0+n*SIMDD, MM_FMA(r0, rg, MM_LOADU(p0+n*SIMDD)));
                }
                for (n = n*SIMDD; n < nf; n++) {
                        p0[n] += c0 * gp[n];
                }
                break;
        case 2:
                c0 = coeff[nprim*sortedidx[0]];
                c1 = coeff[nprim*sortedidx[1]];
                p0 = gc + nf * sortedidx[0];
                p1 = gc + nf * sortedidx[1];
                r0 = MM_SET1(c0);
                r1 = MM_SET1(c1);
                for (n = 0; n < nf/SIMDD; n++) {
                        rg = MM_LOADU(gp+n*SIMDD);
                        MM_STOREU(p0+n*SIMDD, MM_FMA(r0, rg, MM_LOADU(p0+n*SIMDD)));
                        MM_STOREU(p1+n*SIMDD, MM_FMA(r1, rg, MM_LOADU(p1+n*SIMDD)));
                }
                for (n = n*SIMDD; n < nf; n++) {
                        p0[n] += c0 * gp[n];
                        p1[n] += c1 * gp[n];
                }
                break;
        case 3:
                c0 = coeff[nprim*sortedidx[0]];
                c1 = coeff[nprim*sortedidx[1]];
                c2 = coeff[nprim*sortedidx[2]];
                p0 = gc + nf * sortedidx[0];
                p1 = gc + nf * sortedidx[1];
                p2 = gc + nf * sortedidx[2];
                r0 = MM_SET1(c0);
                r1 = MM_SET1(c1);
                r2 = MM_SET1(c2);
                for (n = 0; n < nf/SIMDD; n++) {
                        rg = MM_LOADU(gp+n*SIMDD);
                        MM_STOREU(p0+n*SIMDD, MM_FMA(r0, rg, MM_LOADU(p0+n*SIMDD)));
                        MM_STOREU(p1+n*SIMDD, MM_FMA(r1, rg, MM_LOADU(p1+n*SIMDD)));
                        MM_STOREU(p2+n*SIMDD, MM_FMA(r2, rg, MM_LOADU(p2+n*SIMDD)));
                }
                for (n = n*SIMDD; n < nf; n++) {
                        p0[n] += c0 * gp[n];
                        p1[n] += c1 * gp[n];
                        p2[n] += c2 * gp[n];
                }
                break;
        case 4:
                c0 = coeff[nprim*sortedidx[0]];
                c1 = coeff[nprim*sortedidx[1]];
                c2 = coeff[nprim*sortedidx[2]];
                c3 = coeff[nprim*sortedidx[3]];
                p0 = gc + nf * sortedidx[0];
                p1 = gc + nf * sortedidx[1];
                p2 = gc + nf * sortedidx[2];
                p3 = gc + nf * sortedidx[3];
                r0 = MM_SET1(c0);
                r1 = MM_SET1(c1);
                r2 = MM_SET1(c2);
                r3 = MM_SET1(c3);
                for (n = 0; n < nf/SIMDD; n++) {
                        rg = MM_LOADU(gp+n*SIMDD);
                        MM_STOREU(p0+n*SIMDD, MM_FMA(r0, rg, MM_LOADU(p0+n*SIMDD)));
                        MM_STOREU(p1+n*SIMDD, MM_FMA(r1, rg, MM_LOADU(p1+n*SIMDD)));
                        MM_STOREU(p2+n*SIMDD, MM_FMA(r2, rg, MM_LOADU(p2+n*SIMDD)));
                        MM_STOREU(p3+n*SIMDD, MM_FMA(r3, rg, MM_LOADU(p3+n*SIMDD)));
                }
                for (n = n*SIMDD; n < nf; n++) {
                        p0[n] += c0 * gp[n];
                        p1[n] += c1 * gp[n];
                        p2[n] += c2 * gp[n];
                        p3[n] += c3 * gp[n];
                }
                break;
        default:
                for (i = 0; i < non0ctr-1; i+=2) {
                        c0 = coeff[nprim*sortedidx[i]];
                        c1 = coeff[nprim*sortedidx[i+1]];
                        p0 = gc + nf * sortedidx[i];
                        p1 = gc + nf * sortedidx[i+1];
                        r0 = MM_SET1(c0);
                        r1 = MM_SET1(c1);
                        for (n = 0; n < nf/SIMDD; n++) {
                                rg = MM_LOADU(gp+n*SIMDD);
                                MM_STOREU(p0+n*SIMDD, MM_FMA(r0, rg, MM_LOADU(p0+n*SIMDD)));
                                MM_STOREU(p1+n*SIMDD, MM_FMA(r1, rg, MM_LOADU(p1+n*SIMDD)));
                        }
                        for (n = n*SIMDD; n < nf; n++) {
                                p0[n] += c0 * gp[n];
                                p1[n] += c1 * gp[n];
                        }
                }
                if (i < non0ctr) {
                        c0 = coeff[nprim*sortedidx[i]];
                        p0 = gc + nf * sortedidx[i];
                        r0 = MM_SET1(c0);
                        for (n = 0; n < nf/SIMDD; n++) {
                                rg = MM_LOADU(gp+n*SIMDD);
                                MM_STOREU(p0+n*SIMDD, MM_FMA(r0, rg, MM_LOADU(p0+n*SIMDD)));
                        }
                        for (n = n*SIMDD; n < nf; n++) {
                                p0[n] += c0 * gp[n];
                        }
                }
        }
}

/*
 * to optimize memory copy in cart2sph.c, remove the common factor for s
 * and p function in cart2sph
 */
double CINTcommon_fac_sp(int l)
{
        switch (l) {
                case 0: return 0.282094791773878143;
                case 1: return 0.488602511902919921;
                default: return 1;
        }
}


void CINTiprim_to_ctr_0(double *RESTRICT gc, double *RESTRICT gp, double *RESTRICT coeff, size_t nf,
                        int nprim, int nctr, int non0ctr, int *sortedidx)
{
        int i;
        double c0;

        switch (nf) {
        case 1:
                for (i = 0; i < nctr; i++) {
                        gc[i] = coeff[nprim*i] * gp[0];
                }
                break;
        case 3:
                for (i = 0; i < nctr; i++) {
                        c0 = coeff[nprim*i];
                        gc[i*3+0] = c0 * gp[0];
                        gc[i*3+1] = c0 * gp[1];
                        gc[i*3+2] = c0 * gp[2];
                }
                break;
        case 6:
                for (i = 0; i < nctr; i++) {
                        c0 = coeff[nprim*i];
                        gc[i*6+0] = c0 * gp[0];
                        gc[i*6+1] = c0 * gp[1];
                        gc[i*6+2] = c0 * gp[2];
                        gc[i*6+3] = c0 * gp[3];
                        gc[i*6+4] = c0 * gp[4];
                        gc[i*6+5] = c0 * gp[5];
                }
                break;
        case 9:
                for (i = 0; i < nctr; i++) {
                        c0 = coeff[nprim*i];
                        gc[i*9+0] = c0 * gp[0];
                        gc[i*9+1] = c0 * gp[1];
                        gc[i*9+2] = c0 * gp[2];
                        gc[i*9+3] = c0 * gp[3];
                        gc[i*9+4] = c0 * gp[4];
                        gc[i*9+5] = c0 * gp[5];
                        gc[i*9+6] = c0 * gp[6];
                        gc[i*9+7] = c0 * gp[7];
                        gc[i*9+8] = c0 * gp[8];
                }
                break;
        default:
                CINTprim_to_ctr_0(gc, gp, coeff, nf, nprim, nctr, non0ctr, sortedidx);
        }
}

void CINTiprim_to_ctr_1(double *RESTRICT gc, double *RESTRICT gp, double *RESTRICT coeff, size_t nf,
                        int nprim, int nctr, int non0ctr, int *sortedidx)
{
        int i;
        double c0;
        double *RESTRICT p0;
        switch (nf) {
        case 1:
                for (i = 0; i < nctr; i++) {
                        gc[i] += coeff[nprim*i] * gp[0];
                }
                break;
        case 3:
                for (i = 0; i < nctr; i++) {
                        gc[i*3+0] += coeff[nprim*i] * gp[0];
                        gc[i*3+1] += coeff[nprim*i] * gp[1];
                        gc[i*3+2] += coeff[nprim*i] * gp[2];
                }
                break;
        case 6:
                for (i = 0; i < non0ctr; i++) {
                        c0 = coeff[nprim*sortedidx[i]];
                        p0 = gc + nf * sortedidx[i];
                        p0[0] += c0 * gp[0];
                        p0[1] += c0 * gp[1];
                        p0[2] += c0 * gp[2];
                        p0[3] += c0 * gp[3];
                        p0[4] += c0 * gp[4];
                        p0[5] += c0 * gp[5];
                }
                break;
        case 9:
                for (i = 0; i < non0ctr; i++) {
                        c0 = coeff[nprim*sortedidx[i]];
                        p0 = gc + nf * sortedidx[i];
                        p0[0] += c0 * gp[0];
                        p0[1] += c0 * gp[1];
                        p0[2] += c0 * gp[2];
                        p0[3] += c0 * gp[3];
                        p0[4] += c0 * gp[4];
                        p0[5] += c0 * gp[5];
                        p0[6] += c0 * gp[6];
                        p0[7] += c0 * gp[7];
                        p0[8] += c0 * gp[8];
                }
                break;
        default:
                CINTprim_to_ctr_1(gc, gp, coeff, nf,
                                  nprim, nctr, non0ctr, sortedidx);
        }
}

void CINTsort_gout(double *sout, double *gout, int nf, int count)
{
        if (nf == 1) {
                MM_STORE(sout, MM_LOAD(gout));
                return;
        }

        int i = 0;
#if __AVX2__
#if (SIMDD == 8)
        __m256i vindex = _mm256_set_epi32(7*SIMDD, 6*SIMDD, 5*SIMDD, 4*SIMDD,
                                          3*SIMDD, 2*SIMDD, 1*SIMDD, 0);
#else
        __m128i vindex = _mm_set_epi32(3*SIMDD, 2*SIMDD, 1*SIMDD, 0);
#endif
        switch (count) {
        case 1:
                for (i = 0; i < nf+1-SIMDD; i+=SIMDD) {
                        MM_STOREU(sout     +i, MM_GATHER(gout+i*SIMDD  , vindex, 8));
                }
                break;
        case 2:
                for (i = 0; i < nf+1-SIMDD; i+=SIMDD) {
                        MM_STOREU(sout+0*nf+i, MM_GATHER(gout+i*SIMDD+0, vindex, 8));
                        MM_STOREU(sout+1*nf+i, MM_GATHER(gout+i*SIMDD+1, vindex, 8));
                }
                break;
        case 3:
                for (i = 0; i < nf+1-SIMDD; i+=SIMDD) {
                        MM_STOREU(sout+0*nf+i, MM_GATHER(gout+i*SIMDD+0, vindex, 8));
                        MM_STOREU(sout+1*nf+i, MM_GATHER(gout+i*SIMDD+1, vindex, 8));
                        MM_STOREU(sout+2*nf+i, MM_GATHER(gout+i*SIMDD+2, vindex, 8));
                }
                break;
        case 4:
                for (i = 0; i < nf+1-SIMDD; i+=SIMDD) {
                        MM_STOREU(sout+0*nf+i, MM_GATHER(gout+i*SIMDD+0, vindex, 8));
                        MM_STOREU(sout+1*nf+i, MM_GATHER(gout+i*SIMDD+1, vindex, 8));
                        MM_STOREU(sout+2*nf+i, MM_GATHER(gout+i*SIMDD+2, vindex, 8));
                        MM_STOREU(sout+3*nf+i, MM_GATHER(gout+i*SIMDD+3, vindex, 8));
                }
                break;
#if (SIMDD == 8)
        case 5:
                for (i = 0; i < nf+1-SIMDD; i+=SIMDD) {
                        MM_STOREU(sout+0*nf+i, MM_GATHER(gout+i*SIMDD+0, vindex, 8));
                        MM_STOREU(sout+1*nf+i, MM_GATHER(gout+i*SIMDD+1, vindex, 8));
                        MM_STOREU(sout+2*nf+i, MM_GATHER(gout+i*SIMDD+2, vindex, 8));
                        MM_STOREU(sout+3*nf+i, MM_GATHER(gout+i*SIMDD+3, vindex, 8));
                        MM_STOREU(sout+4*nf+i, MM_GATHER(gout+i*SIMDD+4, vindex, 8));
                }
                break;
        case 6:
                for (i = 0; i < nf+1-SIMDD; i+=SIMDD) {
                        MM_STOREU(sout+0*nf+i, MM_GATHER(gout+i*SIMDD+0, vindex, 8));
                        MM_STOREU(sout+1*nf+i, MM_GATHER(gout+i*SIMDD+1, vindex, 8));
                        MM_STOREU(sout+2*nf+i, MM_GATHER(gout+i*SIMDD+2, vindex, 8));
                        MM_STOREU(sout+3*nf+i, MM_GATHER(gout+i*SIMDD+3, vindex, 8));
                        MM_STOREU(sout+4*nf+i, MM_GATHER(gout+i*SIMDD+4, vindex, 8));
                        MM_STOREU(sout+5*nf+i, MM_GATHER(gout+i*SIMDD+5, vindex, 8));
                }
                break;
        case 7:
                for (i = 0; i < nf+1-SIMDD; i+=SIMDD) {
                        MM_STOREU(sout+0*nf+i, MM_GATHER(gout+i*SIMDD+0, vindex, 8));
                        MM_STOREU(sout+1*nf+i, MM_GATHER(gout+i*SIMDD+1, vindex, 8));
                        MM_STOREU(sout+2*nf+i, MM_GATHER(gout+i*SIMDD+2, vindex, 8));
                        MM_STOREU(sout+3*nf+i, MM_GATHER(gout+i*SIMDD+3, vindex, 8));
                        MM_STOREU(sout+4*nf+i, MM_GATHER(gout+i*SIMDD+4, vindex, 8));
                        MM_STOREU(sout+5*nf+i, MM_GATHER(gout+i*SIMDD+5, vindex, 8));
                        MM_STOREU(sout+6*nf+i, MM_GATHER(gout+i*SIMDD+6, vindex, 8));
                }
                break;
        case 8:
                for (i = 0; i < nf+1-SIMDD; i+=SIMDD) {
                        MM_STOREU(sout+0*nf+i, MM_GATHER(gout+i*SIMDD+0, vindex, 8));
                        MM_STOREU(sout+1*nf+i, MM_GATHER(gout+i*SIMDD+1, vindex, 8));
                        MM_STOREU(sout+2*nf+i, MM_GATHER(gout+i*SIMDD+2, vindex, 8));
                        MM_STOREU(sout+3*nf+i, MM_GATHER(gout+i*SIMDD+3, vindex, 8));
                        MM_STOREU(sout+4*nf+i, MM_GATHER(gout+i*SIMDD+4, vindex, 8));
                        MM_STOREU(sout+5*nf+i, MM_GATHER(gout+i*SIMDD+5, vindex, 8));
                        MM_STOREU(sout+6*nf+i, MM_GATHER(gout+i*SIMDD+6, vindex, 8));
                        MM_STOREU(sout+7*nf+i, MM_GATHER(gout+i*SIMDD+7, vindex, 8));
                }
#endif
        }

        if (i < nf) {
#if (SIMDD == 8)
                vindex = _mm256_set_epi32(nf*7, nf*6, nf*5, nf*4, nf*3, nf*2, nf*1, 0);
                for (; i < nf; i++) {
                        _mm512_i32scatter_pd(sout+i, vindex, MM_LOAD(gout+i*SIMDD), 8);
                }
#else
                switch (count) {
                case 1:
                        for (; i < nf; i++) {
                                sout[i] = gout[i*SIMDD];
                        }
                        break;
                case 2:
                        for (; i < nf; i++) {
                                sout[0*nf+i] = gout[i*SIMDD+0];
                                sout[1*nf+i] = gout[i*SIMDD+1];
                        }
                        break;
                case 3:
                        for (; i < nf; i++) {
                                sout[0*nf+i] = gout[i*SIMDD+0];
                                sout[1*nf+i] = gout[i*SIMDD+1];
                                sout[2*nf+i] = gout[i*SIMDD+2];
                        }
                        break;
                case 4:
                        for (; i < nf; i++) {
                                sout[0*nf+i] = gout[i*SIMDD+0];
                                sout[1*nf+i] = gout[i*SIMDD+1];
                                sout[2*nf+i] = gout[i*SIMDD+2];
                                sout[3*nf+i] = gout[i*SIMDD+3];
                        }
                }
#endif
        }

#else // AVX or SSE3
        switch (count) {
        case 1:
                for (i = 0; i < nf; i++) {
                        sout[i] = gout[i*SIMDD];
                }
                break;
        case 2:
                for (i = 0; i < nf; i++) {
                        sout[0*nf+i] = gout[i*SIMDD+0];
                        sout[1*nf+i] = gout[i*SIMDD+1];
                }
                break;
#if (SIMDD > 2)
        case 3:
                for (i = 0; i < nf; i++) {
                        sout[0*nf+i] = gout[i*SIMDD+0];
                        sout[1*nf+i] = gout[i*SIMDD+1];
                        sout[2*nf+i] = gout[i*SIMDD+2];
                }
                break;
        case 4:
                for (i = 0; i < nf; i++) {
                        sout[0*nf+i] = gout[i*SIMDD+0];
                        sout[1*nf+i] = gout[i*SIMDD+1];
                        sout[2*nf+i] = gout[i*SIMDD+2];
                        sout[3*nf+i] = gout[i*SIMDD+3];
                }
#endif
        }
#endif
}

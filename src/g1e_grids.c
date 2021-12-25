/*
 * Copyright (C) 2021  Qiming Sun <osirpt.sun@gmail.com>
 *
 */

#include <string.h>
#include <math.h>
#include <assert.h>
#include "cint_bas.h"
#include "misc.h"
#include "g1e.h"
#include "g1e_grids.h"
#include "rys_roots.h"
#include "simd.h"


void CINTinit_int1e_grids_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                                  int *atm, int natm,
                                  int *bas, int nbas, double *env)
{
        CINTinit_int1e_EnvVars(envs, ng, shls, atm, natm, bas, nbas, env);
        int ngrids = shls[3] - shls[2];
        double *grids = env + (size_t)env[PTR_GRIDS] + shls[2] * 3;

        envs->ngrids = ngrids;
        envs->grids = grids;
        envs->common_factor = 2 * M_PI
                * CINTcommon_fac_sp(envs->i_l) * CINTcommon_fac_sp(envs->j_l);

        int nroots = envs->nrys_roots;
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
        envs->g_stride_i = GRID_BLKSIZE * nroots;
        envs->g_stride_j = GRID_BLKSIZE * nroots * dli;
        envs->g_size     = GRID_BLKSIZE * nroots * dli * dlj;
        envs->g_stride_k = envs->g_size;
        envs->g_stride_l = envs->g_size;
}

#define RGSQUARE(r, ig)     (r[ig+GRID_BLKSIZE*0]*r[ig+GRID_BLKSIZE*0] + \
                             r[ig+GRID_BLKSIZE*1]*r[ig+GRID_BLKSIZE*1] + \
                             r[ig+GRID_BLKSIZE*2]*r[ig+GRID_BLKSIZE*2])

int CINTg0_1e_grids(double *RESTRICT g, CINTEnvVars *envs,
                    double *cache, double *RESTRICT gridsT)
{
        int ngrids = envs->ngrids;
        int bgrids = MIN(ngrids - envs->grids_offset, GRID_BLKSIZE);
        int nroots = envs->nrys_roots;
        double *RESTRICT gx = g;
        double *RESTRICT gy = g + envs->g_size;
        double *RESTRICT gz = g + envs->g_size * 2;
        double *RESTRICT w = gz;
        double *rij = envs->rij;
        double ubuf[MXRYSROOTS];
        double wbuf[MXRYSROOTS];
        double *RESTRICT u;
        MALLOC_ALIGNED_DOUBLE_INSTACK(u, GRID_BLKSIZE*nroots);
        double *RESTRICT rijrg;
        MALLOC_ALIGNED_DOUBLE_INSTACK(rijrg, GRID_BLKSIZE*3);
        double aij = envs->ai[0] + envs->aj[0];
        int n, i, j, ig;
        double x, fac1;

        __MD r0, r1, r2;
        r1 = MM_SET1(1.);
        for (i = 0; i < nroots; i++) {
                for (ig = 0; ig < bgrids; ig += SIMDD) {
                        MM_STORE(gx+ig+GRID_BLKSIZE*i, r1);
                        MM_STORE(gy+ig+GRID_BLKSIZE*i, r1);
                }
        }

        r0 = MM_SET1(rij[0]);
        r1 = MM_SET1(rij[1]);
        r2 = MM_SET1(rij[2]);
        for (ig = 0; ig < bgrids; ig += SIMDD) {
                MM_STORE(rijrg+ig+GRID_BLKSIZE*0, MM_LOAD(gridsT+ig+GRID_BLKSIZE*0) - r0);
                MM_STORE(rijrg+ig+GRID_BLKSIZE*1, MM_LOAD(gridsT+ig+GRID_BLKSIZE*1) - r1);
                MM_STORE(rijrg+ig+GRID_BLKSIZE*2, MM_LOAD(gridsT+ig+GRID_BLKSIZE*2) - r2);
        }

#ifdef WITH_RANGE_COULOMB
        const double omega = envs->env[PTR_RANGE_OMEGA];
        double theta, sqrt_theta, aij_theta;

        if (omega == 0) {
                fac1 = envs->fac[0] / aij;
                for (ig = 0; ig < bgrids; ig++) {
                        x = aij * RGSQUARE(rijrg, ig);
                        CINTrys_roots(nroots, x, ubuf, wbuf);
                        for (i = 0; i < nroots; i++) {
                                // transform to t^2
                                u[ig+GRID_BLKSIZE*i] = ubuf[i] / (ubuf[i] + 1);
                                w[ig+GRID_BLKSIZE*i] = wbuf[i] * fac1;
                        }
                }
        } else if (omega < 0) { // short-range part of range-separated Coulomb
                theta = omega * omega / (omega * omega + aij);
                sqrt_theta = sqrt(theta);
                fac1 = envs->fac[0] / aij;
                // FIXME:
                // very small erfc() leads to ~0 weights. They can cause
                // numerical issue in sr_rys_roots Use this cutoff as a
                // temporary solution to avoid the numerical issue
                double temp_cutoff = MIN(envs->expcutoff, EXPCUTOFF_SR - nroots);
                for (ig = 0; ig < bgrids; ig++) {
                        x = aij * RGSQUARE(rijrg, ig);
                        if (theta * x > temp_cutoff) {
                                // very small erfc() leads to ~0 weights
                                for (i = 0; i < nroots; i++) {
                                        u[ig+GRID_BLKSIZE*i] = 0;
                                        w[ig+GRID_BLKSIZE*i] = 0;
                                }
                        } else {
                                CINTsr_rys_roots(nroots, x, sqrt_theta, ubuf, wbuf);
                                for (i = 0; i < nroots; i++) {
                                        u[ig+GRID_BLKSIZE*i] = ubuf[i] / (ubuf[i] + 1);
                                        w[ig+GRID_BLKSIZE*i] = wbuf[i] * fac1;
                                }
                        }
                }
        } else {  // long-range part of range-separated Coulomb
                theta = omega * omega / (omega * omega + aij);
                fac1 = envs->fac[0] * sqrt(theta) / aij;
                aij_theta = aij * theta;
                for (ig = 0; ig < bgrids; ig++) {
                        x = aij_theta * RGSQUARE(rijrg, ig);
                        CINTrys_roots(nroots, x, ubuf, wbuf);
                        for (i = 0; i < nroots; i++) {
                                // u stores t^2 = tau^2 * theta
                                u[ig+GRID_BLKSIZE*i] = ubuf[i] / (ubuf[i] + 1) * theta;
                                w[ig+GRID_BLKSIZE*i] = wbuf[i] * fac1;
                        }
                }
        }
#else
        fac1 = envs->fac[0] / aij;
        for (ig = 0; ig < bgrids; ig++) {
                x = aij * RGSQUARE(rijrg, ig);
                CINTrys_roots(nroots, x, ubuf, wbuf);
                for (i = 0; i < nroots; i++) {
                        u[ig+GRID_BLKSIZE*i] = ubuf[i] / (ubuf[i] + 1);
                        w[ig+GRID_BLKSIZE*i] = wbuf[i] * fac1;
                }
        }
#endif
        int nmax = envs->li_ceil + envs->lj_ceil;
        if (nmax == 0) {
                return 1;
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
        __MD rijrx = MM_SET1(rij[0] - rx[0]);
        __MD rijry = MM_SET1(rij[1] - rx[1]);
        __MD rijrz = MM_SET1(rij[2] - rx[2]);
        __MD aij2 = MM_SET1(0.5 / aij);
        __MD ru, rt;
        __MD r0x, r0y, r0z;
        __MD r1x, r1y, r1z;
        __MD r2x, r2y, r2z;

        for (n = 0; n < nroots; n++) {
                p0x = gx + GRID_BLKSIZE*n;
                p0y = gy + GRID_BLKSIZE*n;
                p0z = gz + GRID_BLKSIZE*n;
                p1x = p0x + di;
                p1y = p0y + di;
                p1z = p0z + di;
                p2x = p0x - di;
                p2y = p0y - di;
                p2z = p0z - di;
                for (ig = 0; ig < bgrids; ig += SIMDD) {
                        ru = MM_LOAD(u+ig+GRID_BLKSIZE*n);
                        rt = aij2 - aij2 * ru;
                        r0 = rijrx + ru * MM_LOAD(rijrg+ig+GRID_BLKSIZE*0);
                        r1 = rijry + ru * MM_LOAD(rijrg+ig+GRID_BLKSIZE*1);
                        r2 = rijrz + ru * MM_LOAD(rijrg+ig+GRID_BLKSIZE*2);

                        r2x = MM_LOAD(p0x+ig);
                        r2y = MM_LOAD(p0y+ig);
                        r2z = MM_LOAD(p0z+ig);
                        r0x = r0 * r2x;
                        r0y = r1 * r2y;
                        r0z = r2 * r2z;
                        MM_STORE(p1x+ig, r0x);
                        MM_STORE(p1y+ig, r0y);
                        MM_STORE(p1z+ig, r0z);
                        for (i = 1; i < nmax; i++) {
                                //p1x[ig+i*di] = i * rt * p2x[ig+i*di] + r0 * p0x[ig+i*di];
                                //p1y[ig+i*di] = i * rt * p2y[ig+i*di] + r1 * p0y[ig+i*di];
                                //p1z[ig+i*di] = i * rt * p2z[ig+i*di] + r2 * p0z[ig+i*di];
                                r1x = MM_SET1(i) * rt * r2x + r0 * r0x;
                                r1y = MM_SET1(i) * rt * r2y + r1 * r0y;
                                r1z = MM_SET1(i) * rt * r2z + r2 * r0z;
                                MM_STORE(p1x+ig+i*di, r1x);
                                MM_STORE(p1y+ig+i*di, r1y);
                                MM_STORE(p1z+ig+i*di, r1z);
                                r2x = r0x;
                                r2y = r0y;
                                r2z = r0z;
                                r0x = r1x;
                                r0y = r1y;
                                r0z = r1z;
                        }
                }
        }

        __MD rirjx = MM_SET1(envs->rirj[0]);
        __MD rirjy = MM_SET1(envs->rirj[1]);
        __MD rirjz = MM_SET1(envs->rirj[2]);
        for (j = 1; j <= lj; j++) {
        for (n = 0; n < nroots; n++) {
                p0x = gx + j * dj + GRID_BLKSIZE*n;
                p0y = gy + j * dj + GRID_BLKSIZE*n;
                p0z = gz + j * dj + GRID_BLKSIZE*n;
                p1x = p0x - dj;
                p1y = p0y - dj;
                p1z = p0z - dj;
                p2x = p1x + di;
                p2y = p1y + di;
                p2z = p1z + di;

                for (ig = 0; ig < bgrids; ig += SIMDD) {
                        r1x = MM_LOAD(p1x+ig);
                        r1y = MM_LOAD(p1y+ig);
                        r1z = MM_LOAD(p1z+ig);
                        for (i = 0; i <= nmax - j; i++) {
                                //p0x[ig+i*di] = p2x[ig+i*di] + rirjx * p1x[ig+i*di];
                                //p0y[ig+i*di] = p2y[ig+i*di] + rirjy * p1y[ig+i*di];
                                //p0z[ig+i*di] = p2z[ig+i*di] + rirjz * p1z[ig+i*di];
                                r2x = MM_LOAD(p2x+ig+i*di);
                                r2y = MM_LOAD(p2y+ig+i*di);
                                r2z = MM_LOAD(p2z+ig+i*di);
                                MM_STORE(p0x+ig+i*di, r2x + rirjx * r1x);
                                MM_STORE(p0y+ig+i*di, r2y + rirjy * r1y);
                                MM_STORE(p0z+ig+i*di, r2z + rirjz * r1z);
                                r1x = r2x;
                                r1y = r2y;
                                r1z = r2z;
                        }
                }
        } }
        return 1;
}

void CINTgout1e_grids(double *gout, double *g, int *idx,
                      CINTEnvVars *envs, int gout_empty)
{
        int ngrids = envs->ngrids;
        int bgrids = MIN(ngrids - envs->grids_offset, GRID_BLKSIZE);
        bgrids = ALIGN_UP(bgrids, SIMDD);
        int nroots = envs->nrys_roots;
        int nf = envs->nf;
        int i, n, ig;
        double *gx, *gy, *gz;
        __MD s;

        if (gout_empty) {
                switch(nroots) {
                case 1:
                        for (n = 0; n < nf; n++) {
                                gx = g + idx[n*3+0];
                                gy = g + idx[n*3+1];
                                gz = g + idx[n*3+2];
                                for (ig = 0; ig < bgrids; ig += SIMDD) {
                                        s = MM_LOAD(gx+ig) * MM_LOAD(gy+ig) * MM_LOAD(gz+ig);
                                        MM_STORE(gout+ig+bgrids*n, s);
                                }
                        }
                        return;
                case 2:
                        for (n = 0; n < nf; n++) {
                                gx = g + idx[n*3+0];
                                gy = g + idx[n*3+1];
                                gz = g + idx[n*3+2];
                                for (ig = 0; ig < bgrids; ig += SIMDD) {
                                        s  = MM_LOAD(gx+ig             ) * MM_LOAD(gy+ig             ) * MM_LOAD(gz+ig             );
                                        s += MM_LOAD(gx+ig+GRID_BLKSIZE) * MM_LOAD(gy+ig+GRID_BLKSIZE) * MM_LOAD(gz+ig+GRID_BLKSIZE);
                                        MM_STORE(gout+ig+bgrids*n, s);
                                }
                        }
                        break;
                case 3:
                        for (n = 0; n < nf; n++) {
                                gx = g + idx[n*3+0];
                                gy = g + idx[n*3+1];
                                gz = g + idx[n*3+2];
                                for (ig = 0; ig < bgrids; ig += SIMDD) {
                                        s  = MM_LOAD(gx+ig               ) * MM_LOAD(gy+ig               ) * MM_LOAD(gz+ig               );
                                        s += MM_LOAD(gx+ig+GRID_BLKSIZE  ) * MM_LOAD(gy+ig+GRID_BLKSIZE  ) * MM_LOAD(gz+ig+GRID_BLKSIZE  );
                                        s += MM_LOAD(gx+ig+GRID_BLKSIZE*2) * MM_LOAD(gy+ig+GRID_BLKSIZE*2) * MM_LOAD(gz+ig+GRID_BLKSIZE*2);
                                        MM_STORE(gout+ig+bgrids*n, s);
                                }
                        }
                        break;
                case 4:
                        for (n = 0; n < nf; n++) {
                                gx = g + idx[n*3+0];
                                gy = g + idx[n*3+1];
                                gz = g + idx[n*3+2];
                                for (ig = 0; ig < bgrids; ig += SIMDD) {
                                        s  = MM_LOAD(gx+ig               ) * MM_LOAD(gy+ig               ) * MM_LOAD(gz+ig               );
                                        s += MM_LOAD(gx+ig+GRID_BLKSIZE  ) * MM_LOAD(gy+ig+GRID_BLKSIZE  ) * MM_LOAD(gz+ig+GRID_BLKSIZE  );
                                        s += MM_LOAD(gx+ig+GRID_BLKSIZE*2) * MM_LOAD(gy+ig+GRID_BLKSIZE*2) * MM_LOAD(gz+ig+GRID_BLKSIZE*2);
                                        s += MM_LOAD(gx+ig+GRID_BLKSIZE*3) * MM_LOAD(gy+ig+GRID_BLKSIZE*3) * MM_LOAD(gz+ig+GRID_BLKSIZE*3);
                                        MM_STORE(gout+ig+bgrids*n, s);
                                }
                        }
                        break;
                default:
                        for (n = 0; n < nf; n++) {
                                gx = g + idx[n*3+0];
                                gy = g + idx[n*3+1];
                                gz = g + idx[n*3+2];
                                for (ig = 0; ig < bgrids; ig += SIMDD) {
                                        s = MM_SET1(0);
                                        for (i = 0; i < nroots; i++) {
                                                s += MM_LOAD(gx+ig+GRID_BLKSIZE*i) * MM_LOAD(gy+ig+GRID_BLKSIZE*i) * MM_LOAD(gz+ig+GRID_BLKSIZE*i);
                                        }
                                        MM_STORE(gout+ig+bgrids*n, s);
                                }
                        }
                }
        } else {
                switch(nroots) {
                case 1:
                        for (n = 0; n < nf; n++) {
                                gx = g + idx[n*3+0];
                                gy = g + idx[n*3+1];
                                gz = g + idx[n*3+2];
                                for (ig = 0; ig < bgrids; ig += SIMDD) {
                                        s = MM_LOAD(gout+ig+bgrids*n);
                                        s += MM_LOAD(gx+ig) * MM_LOAD(gy+ig) * MM_LOAD(gz+ig);
                                        MM_STORE(gout+ig+bgrids*n, s);
                                }
                        }
                        return;
                case 2:
                        for (n = 0; n < nf; n++) {
                                gx = g + idx[n*3+0];
                                gy = g + idx[n*3+1];
                                gz = g + idx[n*3+2];
                                for (ig = 0; ig < bgrids; ig += SIMDD) {
                                        s = MM_LOAD(gout+ig+bgrids*n);
                                        s += MM_LOAD(gx+ig             ) * MM_LOAD(gy+ig             ) * MM_LOAD(gz+ig             );
                                        s += MM_LOAD(gx+ig+GRID_BLKSIZE) * MM_LOAD(gy+ig+GRID_BLKSIZE) * MM_LOAD(gz+ig+GRID_BLKSIZE);
                                        MM_STORE(gout+ig+bgrids*n, s);
                                }
                        }
                        break;
                case 3:
                        for (n = 0; n < nf; n++) {
                                gx = g + idx[n*3+0];
                                gy = g + idx[n*3+1];
                                gz = g + idx[n*3+2];
                                for (ig = 0; ig < bgrids; ig += SIMDD) {
                                        s = MM_LOAD(gout+ig+bgrids*n);
                                        s += MM_LOAD(gx+ig               ) * MM_LOAD(gy+ig               ) * MM_LOAD(gz+ig               );
                                        s += MM_LOAD(gx+ig+GRID_BLKSIZE  ) * MM_LOAD(gy+ig+GRID_BLKSIZE  ) * MM_LOAD(gz+ig+GRID_BLKSIZE  );
                                        s += MM_LOAD(gx+ig+GRID_BLKSIZE*2) * MM_LOAD(gy+ig+GRID_BLKSIZE*2) * MM_LOAD(gz+ig+GRID_BLKSIZE*2);
                                        MM_STORE(gout+ig+bgrids*n, s);
                                }
                        }
                        break;
                case 4:
                        for (n = 0; n < nf; n++) {
                                gx = g + idx[n*3+0];
                                gy = g + idx[n*3+1];
                                gz = g + idx[n*3+2];
                                for (ig = 0; ig < bgrids; ig += SIMDD) {
                                        s = MM_LOAD(gout+ig+bgrids*n);
                                        s += MM_LOAD(gx+ig               ) * MM_LOAD(gy+ig               ) * MM_LOAD(gz+ig               );
                                        s += MM_LOAD(gx+ig+GRID_BLKSIZE  ) * MM_LOAD(gy+ig+GRID_BLKSIZE  ) * MM_LOAD(gz+ig+GRID_BLKSIZE  );
                                        s += MM_LOAD(gx+ig+GRID_BLKSIZE*2) * MM_LOAD(gy+ig+GRID_BLKSIZE*2) * MM_LOAD(gz+ig+GRID_BLKSIZE*2);
                                        s += MM_LOAD(gx+ig+GRID_BLKSIZE*3) * MM_LOAD(gy+ig+GRID_BLKSIZE*3) * MM_LOAD(gz+ig+GRID_BLKSIZE*3);
                                        MM_STORE(gout+ig+bgrids*n, s);
                                }
                        }
                        break;
                default:
                        for (n = 0; n < nf; n++) {
                                gx = g + idx[n*3+0];
                                gy = g + idx[n*3+1];
                                gz = g + idx[n*3+2];
                                for (ig = 0; ig < bgrids; ig += SIMDD) {
                                        s = MM_LOAD(gout+ig+bgrids*n);
                                        for (i = 0; i < nroots; i++) {
                                                s += MM_LOAD(gx+ig+GRID_BLKSIZE*i) * MM_LOAD(gy+ig+GRID_BLKSIZE*i) * MM_LOAD(gz+ig+GRID_BLKSIZE*i);
                                        }
                                        MM_STORE(gout+ig+bgrids*n, s);
                                }
                        }
                }
        }
}

void CINTnabla1i_grids(double *f, double *g,
                       int li, int lj, CINTEnvVars *envs)
{
        int ngrids = envs->ngrids;
        int bgrids = MIN(ngrids - envs->grids_offset, GRID_BLKSIZE);
        bgrids = ALIGN_UP(bgrids, SIMDD);
        int nroots = envs->nrys_roots;
        const int di = envs->g_stride_i;
        const int dj = envs->g_stride_j;
        int i, j, n, ig, ptr;
        const double *gx = g;
        const double *gy = g + envs->g_size;
        const double *gz = g + envs->g_size * 2;
        double *fx = f;
        double *fy = f + envs->g_size;
        double *fz = f + envs->g_size * 2;
        __MD ai2 = MM_SET1(-2 * envs->ai[0]);

        for (j = 0; j <= lj; j++) {
                //f(...,0,...) = -2*ai*g(...,1,...)
                for (n = 0; n < nroots; n++) {
                        ptr = dj * j + n * GRID_BLKSIZE;
                        for (ig = ptr; ig < ptr+bgrids; ig += SIMDD) {
                                MM_STORE(fx+ig, ai2 * MM_LOAD(gx+ig+di));
                                MM_STORE(fy+ig, ai2 * MM_LOAD(gy+ig+di));
                                MM_STORE(fz+ig, ai2 * MM_LOAD(gz+ig+di));
                        }
                }
                //f(...,i,...) = i*g(...,i-1,...)-2*ai*g(...,i+1,...)
                for (i = 1; i <= li; i++) {
                for (n = 0; n < nroots; n++) {
                        ptr = dj * j + di * i + n * GRID_BLKSIZE;
                        for (ig = ptr; ig < ptr+bgrids; ig += SIMDD) {
                                MM_STORE(fx+ig, MM_SET1(i) * MM_LOAD(gx+ig-di) + ai2 * MM_LOAD(gx+ig+di));
                                MM_STORE(fy+ig, MM_SET1(i) * MM_LOAD(gy+ig-di) + ai2 * MM_LOAD(gy+ig+di));
                                MM_STORE(fz+ig, MM_SET1(i) * MM_LOAD(gz+ig-di) + ai2 * MM_LOAD(gz+ig+di));
                        }
                } }
        }
}

void CINTnabla1j_grids(double *f, double *g,
                       int li, int lj, CINTEnvVars *envs)
{
        int ngrids = envs->ngrids;
        int bgrids = MIN(ngrids - envs->grids_offset, GRID_BLKSIZE);
        bgrids = ALIGN_UP(bgrids, SIMDD);
        int nroots = envs->nrys_roots;
        const int di = envs->g_stride_i;
        const int dj = envs->g_stride_j;
        int i, j, n, ig, ptr;
        const double *gx = g;
        const double *gy = g + envs->g_size;
        const double *gz = g + envs->g_size * 2;
        double *fx = f;
        double *fy = f + envs->g_size;
        double *fz = f + envs->g_size * 2;
        __MD aj2 = MM_SET1(-2 * envs->aj[0]);

        //f(...,0,...) = -2*aj*g(...,1,...)
        for (i = 0; i <= li; i++) {
        for (n = 0; n < nroots; n++) {
                ptr = di * i + n * GRID_BLKSIZE;
                for (ig = ptr; ig < ptr+bgrids; ig += SIMDD) {
                        MM_STORE(fx+ig, aj2 * MM_LOAD(gx+ig+dj));
                        MM_STORE(fy+ig, aj2 * MM_LOAD(gy+ig+dj));
                        MM_STORE(fz+ig, aj2 * MM_LOAD(gz+ig+dj));
                }
        } }
        //f(...,j,...) = j*g(...,j-1,...)-2*aj*g(...,j+1,...)
        for (j = 1; j <= lj; j++) {
        for (i = 0; i <= li; i++) {
        for (n = 0; n < nroots; n++) {
                ptr = dj * j + di * i + n * GRID_BLKSIZE;
                for (ig = ptr; ig < ptr+bgrids; ig += SIMDD) {
                        MM_STORE(fx+ig, MM_SET1(j) * MM_LOAD(gx+ig-dj) + aj2 * MM_LOAD(gx+ig+dj));
                        MM_STORE(fy+ig, MM_SET1(j) * MM_LOAD(gy+ig-dj) + aj2 * MM_LOAD(gy+ig+dj));
                        MM_STORE(fz+ig, MM_SET1(j) * MM_LOAD(gz+ig-dj) + aj2 * MM_LOAD(gz+ig+dj));
                }
        } } }
}


void CINTx1i_grids(double *f, double *g, double *ri,
                   int li, int lj, CINTEnvVars *envs)
{
        int ngrids = envs->ngrids;
        int bgrids = MIN(ngrids - envs->grids_offset, GRID_BLKSIZE);
        bgrids = ALIGN_UP(bgrids, SIMDD);
        int nroots = envs->nrys_roots;
        int i, j, n, ig, ptr;
        const int di = envs->g_stride_i;
        const int dj = envs->g_stride_j;
        const double *gx = g;
        const double *gy = g + envs->g_size;
        const double *gz = g + envs->g_size * 2;
        double *fx = f;
        double *fy = f + envs->g_size;
        double *fz = f + envs->g_size * 2;
        __MD r0 = MM_SET1(ri[0]);
        __MD r1 = MM_SET1(ri[1]);
        __MD r2 = MM_SET1(ri[2]);

        for (j = 0; j <= lj; j++) {
        for (i = 0; i <= li; i++) {
        for (n = 0; n < nroots; n++) {
                ptr = dj * j + di * i + n * GRID_BLKSIZE;
                for (ig = ptr; ig < ptr+bgrids; ig += SIMDD) {
                        MM_STORE(fx+ig, MM_LOAD(gx+ig+di) + r0 * MM_LOAD(gx+ig));
                        MM_STORE(fy+ig, MM_LOAD(gy+ig+di) + r1 * MM_LOAD(gy+ig));
                        MM_STORE(fz+ig, MM_LOAD(gz+ig+di) + r2 * MM_LOAD(gz+ig));
                }
        } } }
}

void CINTx1j_grids(double *f, double *g, double *rj,
                   int li, int lj, CINTEnvVars *envs)
{
        int ngrids = envs->ngrids;
        int bgrids = MIN(ngrids - envs->grids_offset, GRID_BLKSIZE);
        bgrids = ALIGN_UP(bgrids, SIMDD);
        int nroots = envs->nrys_roots;
        int i, j, n, ig, ptr;
        const int di = envs->g_stride_i;
        const int dj = envs->g_stride_j;
        const double *gx = g;
        const double *gy = g + envs->g_size;
        const double *gz = g + envs->g_size * 2;
        double *fx = f;
        double *fy = f + envs->g_size;
        double *fz = f + envs->g_size * 2;
        __MD r0 = MM_SET1(rj[0]);
        __MD r1 = MM_SET1(rj[1]);
        __MD r2 = MM_SET1(rj[2]);

        for (j = 0; j <= lj; j++) {
        for (i = 0; i <= li; i++) {
        for (n = 0; n < nroots; n++) {
                ptr = dj * j + di * i + n * GRID_BLKSIZE;
                for (ig = ptr; ig < ptr+bgrids; ig += SIMDD) {
                        MM_STORE(fx+ig, MM_LOAD(gx+ig+dj) + r0 * MM_LOAD(gx+ig));
                        MM_STORE(fy+ig, MM_LOAD(gy+ig+dj) + r1 * MM_LOAD(gy+ig));
                        MM_STORE(fz+ig, MM_LOAD(gz+ig+dj) + r2 * MM_LOAD(gz+ig));
                }
        } } }
}

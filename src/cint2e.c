/*
 * Qcint is a general GTO integral library for computational chemistry
 * Copyright (C) 2014 Qiming Sun
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
#include <pmmintrin.h>
#include <mm_malloc.h>
#include "cint_bas.h"
#include "g2e.h"
#include "optimizer.h"
#include "cint2e.h"
#include "misc.h"
#include "cart2sph.h"
#include "c2f.h"

#define SQUARE(r)       (r)[0]*(r)[0] + (r)[1]*(r)[1] + (r)[2]*(r)[2]

#define PRIM2CTR0(ctrsymb, gp, ngp) \
        if (ctrsymb##_ctr > 1) {\
                if (*ctrsymb##empty) { \
                        CINTprim_to_ctr_0(gctr##ctrsymb, ngp, gp, \
                                          envs->ctrsymb##_prim, \
                                          ctrsymb##_ctr, c##ctrsymb+ctrsymb##p); \
                } else { \
                        CINTprim_to_ctr_1(gctr##ctrsymb, ngp, gp, \
                                          envs->ctrsymb##_prim, \
                                          ctrsymb##_ctr, c##ctrsymb+ctrsymb##p); \
                } \
        } \
        *ctrsymb##empty = 0

void CINTprim_to_ctr_0(double *gc, const int nf, const double *gp,
                       const int nprim, const int nctr, const double *coeff)
{
        int n, i;
        double *p0, *p1, *p2;
        double non0coeff[32];
        int non0idx[32];
        int non0ctr = 0;

        for (i = 0; i < nctr; i++) {
                if (coeff[nprim*i] != 0) {
                        non0coeff[non0ctr] = coeff[nprim*i];
                        non0idx[non0ctr] = i;
                        non0ctr++;
                } else { // need to initialize the memory, since += is used later
                        memset(gc+nf*i, 0, sizeof(double)*nf);
                }
        }

        __m128d r0, r1, r2, r3, r4;
        switch (non0ctr) {
                case 1:
                        r0 = _mm_load1_pd(non0coeff);
                        p0 = gc + nf*non0idx[0];
                        for (n = 0; n < nf-1; n+=2) {
                                r3 = _mm_loadu_pd(&gp[n]);
                                r4 = _mm_mul_pd(r0, r3);
                                _mm_storeu_pd(p0+n, r4);
                        }
                        if (n < nf) {
                                p0[n] = non0coeff[0] * gp[n];
                        }
                        break;
                case 2:
                        r0 = _mm_load1_pd(non0coeff);
                        r1 = _mm_load1_pd(non0coeff+1);
                        p0 = gc + nf*non0idx[0];
                        p1 = gc + nf*non0idx[1];
                        for (n = 0; n < nf-1; n+=2) {
                                r3 = _mm_loadu_pd(&gp[n]);
                                r4 = _mm_mul_pd(r0, r3);
                                _mm_storeu_pd(p0+n, r4);
                                r4 = _mm_mul_pd(r1, r3);
                                _mm_storeu_pd(p1+n, r4);
                        }
                        if (n < nf) {
                                p0[n] = non0coeff[0] * gp[n];
                                p1[n] = non0coeff[1] * gp[n];
                        }
                        break;
                case 3:
                        r0 = _mm_load1_pd(&non0coeff[0]);
                        r1 = _mm_load1_pd(&non0coeff[1]);
                        r2 = _mm_load1_pd(&non0coeff[2]);
                        p0 = gc + nf*non0idx[0];
                        p1 = gc + nf*non0idx[1];
                        p2 = gc + nf*non0idx[2];
                        for (n = 0; n < nf-1; n+=2) {
                                r3 = _mm_loadu_pd(&gp[n]);
                                r4 = _mm_mul_pd(r0, r3);
                                _mm_storeu_pd(p0+n, r4);
                                r4 = _mm_mul_pd(r1, r3);
                                _mm_storeu_pd(p1+n, r4);
                                r4 = _mm_mul_pd(r2, r3);
                                _mm_storeu_pd(p2+n, r4);
                        }
                        if (n < nf) {
                                p0[n] = non0coeff[0] * gp[n];
                                p1[n] = non0coeff[1] * gp[n];
                                p2[n] = non0coeff[2] * gp[n];
                        }
                        break;
                default:
                        for (i = 0; i < non0ctr-1; i+=2) {
                                r0 = _mm_load1_pd(&non0coeff[i  ]);
                                r1 = _mm_load1_pd(&non0coeff[i+1]);
                                p0 = gc + nf*non0idx[i  ];
                                p1 = gc + nf*non0idx[i+1];
                                for (n = 0; n < nf-1; n+=2) {
                                        r3 = _mm_loadu_pd(&gp[n]);
                                        r4 = _mm_mul_pd(r0, r3);
                                        _mm_storeu_pd(p0+n, r4);
                                        r4 = _mm_mul_pd(r1, r3);
                                        _mm_storeu_pd(p1+n, r4);
                                }
                                if (n < nf) {
                                        p0[n] = non0coeff[i  ] * gp[n];
                                        p1[n] = non0coeff[i+1] * gp[n];
                                }
                        }
                        if (i < non0ctr) {
                                r0 = _mm_load1_pd(&non0coeff[i]);
                                p0 = gc + nf*non0idx[i];
                                for (n = 0; n < nf-1; n+=2) {
                                        r3 = _mm_loadu_pd(&gp[n]);
                                        r4 = _mm_mul_pd(r0, r3);
                                        _mm_storeu_pd(p0+n, r4);
                                }
                                if (n < nf) {
                                        p0[n] = non0coeff[i] * gp[n];
                                }
                        }
        }
}

static void prim_to_ctr_opt(double *gc, const int nf, const double *gp,
                            double *non0coeff, int *non0idx, int non0ctr)
{
        int n, i;
        double *p0, *p1, *p2;

        __m128d r0, r1, r2, r3, r4, r5;
        switch (non0ctr) {
                case 1:
                        r0 = _mm_load1_pd(non0coeff);
                        p0 = gc + nf*non0idx[0];
                        for (n = 0; n < nf-1; n+=2) {
                                r3 = _mm_loadu_pd(&gp[n]);
                                r4 = _mm_loadu_pd(&p0[n]);
                                r5 = _mm_mul_pd(r0, r3);
                                r4 = _mm_add_pd(r5, r4);
                                _mm_storeu_pd(p0+n, r4);
                        }
                        if (n < nf) {
                                p0[n] += non0coeff[0] * gp[n];
                        }
                        break;
                case 2:
                        r0 = _mm_load1_pd(non0coeff);
                        r1 = _mm_load1_pd(non0coeff+1);
                        p0 = gc + nf*non0idx[0];
                        p1 = gc + nf*non0idx[1];
                        for (n = 0; n < nf-1; n+=2) {
                                r3 = _mm_loadu_pd(&gp[n]);
                                r4 = _mm_loadu_pd(&p0[n]);
                                r5 = _mm_mul_pd(r0, r3);
                                r4 = _mm_add_pd(r5, r4);
                                _mm_storeu_pd(p0+n, r4);
                                r4 = _mm_loadu_pd(&p1[n]);
                                r5 = _mm_mul_pd(r1, r3);
                                r4 = _mm_add_pd(r5, r4);
                                _mm_storeu_pd(p1+n, r4);
                        }
                        if (n < nf) {
                                p0[n] += non0coeff[0] * gp[n];
                                p1[n] += non0coeff[1] * gp[n];
                        }
                        break;
                case 3:
                        r0 = _mm_load1_pd(&non0coeff[0]);
                        r1 = _mm_load1_pd(&non0coeff[1]);
                        r2 = _mm_load1_pd(&non0coeff[2]);
                        p0 = gc + nf*non0idx[0];
                        p1 = gc + nf*non0idx[1];
                        p2 = gc + nf*non0idx[2];
                        for (n = 0; n < nf-1; n+=2) {
                                r3 = _mm_loadu_pd(&gp[n]);
                                r4 = _mm_loadu_pd(&p0[n]);
                                r5 = _mm_mul_pd(r0, r3);
                                r4 = _mm_add_pd(r5, r4);
                                _mm_storeu_pd(p0+n, r4);
                                r4 = _mm_loadu_pd(&p1[n]);
                                r5 = _mm_mul_pd(r1, r3);
                                r4 = _mm_add_pd(r5, r4);
                                _mm_storeu_pd(p1+n, r4);
                                r4 = _mm_loadu_pd(&p2[n]);
                                r5 = _mm_mul_pd(r2, r3);
                                r4 = _mm_add_pd(r5, r4);
                                _mm_storeu_pd(p2+n, r4);
                        }
                        if (n < nf) {
                                p0[n] += non0coeff[0] * gp[n];
                                p1[n] += non0coeff[1] * gp[n];
                                p2[n] += non0coeff[2] * gp[n];
                        }
                        break;
                default:
                        for (i = 0; i < non0ctr-1; i+=2) {
                                r0 = _mm_load1_pd(&non0coeff[i  ]);
                                r1 = _mm_load1_pd(&non0coeff[i+1]);
                                p0 = gc + nf*non0idx[i  ];
                                p1 = gc + nf*non0idx[i+1];
                                for (n = 0; n < nf-1; n+=2) {
                                        r3 = _mm_loadu_pd(&gp[n]);
                                        r4 = _mm_loadu_pd(&p0[n]);
                                        r5 = _mm_mul_pd(r0, r3);
                                        r4 = _mm_add_pd(r5, r4);
                                        _mm_storeu_pd(p0+n, r4);
                                        r4 = _mm_loadu_pd(&p1[n]);
                                        r3 = _mm_mul_pd(r1, r3);
                                        r4 = _mm_add_pd(r3, r4);
                                        _mm_storeu_pd(p1+n, r4);
                                }
                                if (n < nf) {
                                        p0[n] += non0coeff[i  ] * gp[n];
                                        p1[n] += non0coeff[i+1] * gp[n];
                                }
                        }
                        if (i < non0ctr) {
                                r0 = _mm_load1_pd(&non0coeff[i]);
                                p0 = gc + nf*non0idx[i];
                                for (n = 0; n < nf-1; n+=2) {
                                        r3 = _mm_loadu_pd(&gp[n]);
                                        r4 = _mm_loadu_pd(&p0[n]);
                                        r5 = _mm_mul_pd(r0, r3);
                                        r4 = _mm_add_pd(r5, r4);
                                        _mm_storeu_pd(p0+n, r4);
                                }
                                if (n < nf) {
                                        p0[n] += non0coeff[i] * gp[n];
                                }
                        }
        }
}

void CINTprim_to_ctr_1(double *gc, const int nf, const double *gp,
                       const int nprim, const int nctr, const double *coeff)
{
        int i;
        double non0coeff[32];
        int non0idx[32];
        int non0ctr = 0;

        for (i = 0; i < nctr; i++) {
                if (coeff[nprim*i] != 0) {
                        non0coeff[non0ctr] = coeff[nprim*i];
                        non0idx[non0ctr] = i;
                        non0ctr++;
                }
        }
        prim_to_ctr_opt(gc, nf, gp, non0coeff, non0idx, non0ctr);
}

/*
static int *_init_index_array(CINTEnvVars *envs, const CINTOpt *opt, const int nf)
{
        int *idx;
        if (opt) {
                if (opt->index_xyz_array) {
                        idx = opt->index_xyz_array[envs->i_l*ANG_MAX*ANG_MAX*ANG_MAX
                                                  +envs->j_l*ANG_MAX*ANG_MAX
                                                  +envs->k_l*ANG_MAX
                                                  +envs->l_l];
                } else {
                        idx = (int *)malloc(sizeof(int) * nf * 3);
                        CINTg2e_index_xyz(idx, envs);
                }
        } else {
                idx = (int *)malloc(sizeof(int) * nf * 3);
                CINTg2e_index_xyz(idx, envs);
        }
        return idx;
} */

inline void CINT2e_core(double *gout, double *g, double fac1i,
                        CINTEnvVars *envs, int empty)
{
        if (envs->g_size == 1) {
                if (empty) {
                        *gout = CINTg0_2e_ssss(fac1i, envs);
                } else { // if same address for gctri and gout and
                         // it has been initialized.
                        *gout+= CINTg0_2e_ssss(fac1i, envs);
                }
        } else {
                CINTg0_2e(g, fac1i, envs);
                (*envs->f_gout)(g, gout, envs->idx, envs, empty);
        }
}

int CINT2e_loop_nopt(double *gctr, CINTEnvVars *envs)
{
        /* COMMON_ENVS_AND_DECLARE */
        const int *shls  = envs->shls;
        const int *bas = envs->bas;
        const double *env = envs->env;
        const int i_ctr  = envs->i_ctr;
        const int j_ctr  = envs->j_ctr;
        const int k_ctr  = envs->k_ctr;
        const int l_ctr  = envs->l_ctr;
        const double *ri = envs->ri;
        const double *rj = envs->rj;
        const double *rk = envs->rk;
        const double *rl = envs->rl;
        const int i_sh = shls[0];
        const int j_sh = shls[1];
        const int k_sh = shls[2];
        const int l_sh = shls[3];
        const double *ai = env + bas(PTR_EXP, i_sh);
        const double *aj = env + bas(PTR_EXP, j_sh);
        const double *ak = env + bas(PTR_EXP, k_sh);
        const double *al = env + bas(PTR_EXP, l_sh);
        const double *ci = env + bas(PTR_COEFF, i_sh);
        const double *cj = env + bas(PTR_COEFF, j_sh);
        const double *ck = env + bas(PTR_COEFF, k_sh);
        const double *cl = env + bas(PTR_COEFF, l_sh);
        const int n_comp = envs->ncomp_e1 * envs->ncomp_e2
                                  * envs->ncomp_tensor;
        double fac1i, fac1j, fac1k, fac1l;
        int ip, jp, kp, lp, n;
        int empty[5] = {1, 1, 1, 1, 1};
        int *iempty = empty + 0;
        int *jempty = empty + 1;
        int *kempty = empty + 2;
        int *lempty = empty + 3;
        int *gempty = empty + 4;
        /* COMMON_ENVS_AND_DECLARE end */
        const int nc = i_ctr * j_ctr * k_ctr * l_ctr;
        // (irys,i,j,k,l,coord,0:1); +1 for nabla-r12
        const int leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const int lenl = envs->nf * nc * n_comp; // gctrl
        const int lenk = envs->nf * i_ctr * j_ctr * k_ctr * n_comp; // gctrk
        const int lenj = envs->nf * i_ctr * j_ctr * n_comp; // gctrj
        const int leni = envs->nf * i_ctr * n_comp; // gctri
        const int len0 = envs->nf * n_comp; // gout
        const int len = leng + lenl + lenk + lenj + leni + len0;
        double *const g = (double *)_mm_malloc(sizeof(double)*len, 16);
        double *g1 = g + leng;
        double *gout, *gctri, *gctrj, *gctrk, *gctrl;

        if (n_comp == 1) {
                gctrl = gctr;
        } else {
                gctrl = g1;
                g1 += lenl;
        }
        if (l_ctr == 1) {
                gctrk = gctrl;
                kempty = lempty;
        } else {
                gctrk = g1;
                g1 += lenk;
        }
        if (k_ctr == 1) {
                gctrj = gctrk;
                jempty = kempty;
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

        double eij, ekl, expijkl;
        const double dist_ij = SQUARE(envs->rirj);
        const double dist_kl = SQUARE(envs->rkrl);
        envs->idx = (int *)malloc(sizeof(int) * envs->nf * 3);
        CINTg2e_index_xyz(envs->idx, envs);

        *lempty = 1;
        for (lp = 0; lp < envs->l_prim; lp++) {
                envs->al = al[lp];
                if (l_ctr == 1) {
                        fac1l = envs->common_factor * cl[lp];
                } else {
                        fac1l = envs->common_factor;
                        *kempty = 1;
                }
                for (kp = 0; kp < envs->k_prim; kp++) {
                        envs->ak = ak[kp];
                        envs->akl = ak[kp] + al[lp];
                        ekl = dist_kl * ak[kp] * al[lp] / envs->akl;
                        if (ekl > EXPCUTOFF) {
                                goto k_contracted;
                        }
                        envs->rkl[0] = (ak[kp]*rk[0] + al[lp]*rl[0]) / envs->akl;
                        envs->rkl[1] = (ak[kp]*rk[1] + al[lp]*rl[1]) / envs->akl;
                        envs->rkl[2] = (ak[kp]*rk[2] + al[lp]*rl[2]) / envs->akl;
                        envs->rklrx[0] = envs->rkl[0] - envs->rx_in_rklrx[0];
                        envs->rklrx[1] = envs->rkl[1] - envs->rx_in_rklrx[1];
                        envs->rklrx[2] = envs->rkl[2] - envs->rx_in_rklrx[2];
                        if (k_ctr == 1) {
                                fac1k = fac1l * ck[kp];
                        } else {
                                fac1k = fac1l;
                                *jempty = 1;
                        }

                        for (jp = 0; jp < envs->j_prim; jp++) {
                                envs->aj = aj[jp];
                                if (j_ctr == 1) {
                                        fac1j = fac1k * cj[jp];
                                } else {
                                        fac1j = fac1k;
                                        *iempty = 1;
                                }
                                for (ip = 0; ip < envs->i_prim; ip++) {
                                        envs->ai = ai[ip];
                                        envs->aij = ai[ip] + aj[jp];
                                        eij = dist_ij * ai[ip] * aj[jp] / envs->aij;
                                        if (eij > EXPCUTOFF) {
                                                goto i_contracted;
                                        }
                                        envs->rij[0] = (ai[ip]*ri[0] + aj[jp]*rj[0]) / envs->aij;
                                        envs->rij[1] = (ai[ip]*ri[1] + aj[jp]*rj[1]) / envs->aij;
                                        envs->rij[2] = (ai[ip]*ri[2] + aj[jp]*rj[2]) / envs->aij;
                                        envs->rijrx[0] = envs->rij[0] - envs->rx_in_rijrx[0];
                                        envs->rijrx[1] = envs->rij[1] - envs->rx_in_rijrx[1];
                                        envs->rijrx[2] = envs->rij[2] - envs->rx_in_rijrx[2];
                                        expijkl = exp(-(eij+ekl));
                                        if (i_ctr == 1) {
                                                fac1i = fac1j*ci[ip]*expijkl;
                                        } else {
                                                fac1i = fac1j*expijkl;
                                        }
                                        CINT2e_core(gout, g, fac1i, envs, *gempty);
                                        PRIM2CTR0(i, gout, envs->nf*n_comp);
i_contracted: ;
                                } // end loop i_prim
                                if (!*iempty) {
                                        PRIM2CTR0(j, gctri, envs->nf*i_ctr*n_comp);
                                }
                        } // end loop j_prim
                        if (!*jempty) {
                                PRIM2CTR0(k, gctrj,envs->nf*i_ctr*j_ctr*n_comp);
                        }
k_contracted: ;
                } // end loop k_prim
                if (!*kempty) {
                        PRIM2CTR0(l, gctrk, envs->nf*i_ctr*j_ctr*k_ctr*n_comp);
                }
        } // end loop l_prim

        /* COPY_AND_CLOSING(gctrl, *lempty); */
        if (n_comp > 1 && !*lempty) {
                const int INC1 = 1;
                double *gctr1, *gctr2, *gctr3;
                switch (n_comp) {
                case 3:
                        gctr1 = gctr  + envs->nf*nc;
                        gctr2 = gctr1 + envs->nf*nc;
                        for (n = 0, ip = 0; n < envs->nf*nc; n++, ip+=3) {
                                gctr [n] = gctrl[ip+0];
                                gctr1[n] = gctrl[ip+1];
                                gctr2[n] = gctrl[ip+2];
                        }
                        break;
                default:
                        for (kp = 0; kp < n_comp-3; kp+=4) {
                                gctr1 = gctr  + envs->nf*nc;
                                gctr2 = gctr1 + envs->nf*nc;
                                gctr3 = gctr2 + envs->nf*nc;
                                for (n = 0, ip = kp; n < envs->nf*nc; n++,ip+=n_comp) {
                                        gctr [n] = gctrl[ip+0];
                                        gctr1[n] = gctrl[ip+1];
                                        gctr2[n] = gctrl[ip+2];
                                        gctr3[n] = gctrl[ip+3];
                                }
                                gctr += envs->nf*nc * 4;
                        }
                        for (; kp < n_comp; kp++) {
                                n = envs->nf*nc;
                                dcopy_(&n, gctrl+kp, &n_comp, gctr, &INC1);
                                gctr += envs->nf*nc;
                        }
                }
        }
        _mm_free(g);
        free(envs->idx);
        return !*lempty;
        /* COPY_AND_CLOSING(gctrl, *lempty); end */
}


#define COMMON_ENVS_AND_DECLARE \
        const int *shls = envs->shls; \
        const int *bas = envs->bas; \
        const double *env = envs->env; \
        const int i_ctr  = envs->i_ctr; \
        const int j_ctr  = envs->j_ctr; \
        const int k_ctr  = envs->k_ctr; \
        const int l_ctr  = envs->l_ctr; \
        const int i_sh = shls[0]; \
        const int j_sh = shls[1]; \
        const int k_sh = shls[2]; \
        const int l_sh = shls[3]; \
        const double *ai = env + bas(PTR_EXP, i_sh); \
        const double *aj = env + bas(PTR_EXP, j_sh); \
        const double *ak = env + bas(PTR_EXP, k_sh); \
        const double *al = env + bas(PTR_EXP, l_sh); \
        const double *ci = env + bas(PTR_COEFF, i_sh); \
        const double *cj = env + bas(PTR_COEFF, j_sh); \
        const double *ck = env + bas(PTR_COEFF, k_sh); \
        const double *cl = env + bas(PTR_COEFF, l_sh); \
        const int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor; \
        double fac1i, fac1j, fac1k, fac1l; \
        int ip, jp, kp, lp, n; \
        int empty[5] = {1, 1, 1, 1, 1}; \
        int *iempty = empty + 0; \
        int *jempty = empty + 1; \
        int *kempty = empty + 2; \
        int *lempty = empty + 3; \
        int *gempty = empty + 4;

#define USE_OPT \
        int off; \
        double expij, expkl; \
        double *prij; \
        const int io = opt->prim_offset[i_sh]; \
        const int jo = opt->prim_offset[j_sh]; \
        const int ko = opt->prim_offset[k_sh]; \
        const int lo = opt->prim_offset[l_sh]; \
        envs->idx = opt->index_xyz_array[envs->i_l*ANG_MAX*ANG_MAX*ANG_MAX \
                                        +envs->j_l*ANG_MAX*ANG_MAX \
                                        +envs->k_l*ANG_MAX \
                                        +envs->l_l]

#define SET_RIJ(I,J)    \
        envs->a##I = a##I[I##p]; \
        envs->a##I##J = a##I[I##p] + a##J[J##p]; \
        off = I##o + I##p; \
        if (opt->cceij[J##o+J##p][off] > CUTOFF15) { \
                goto I##_contracted; } \
        exp##I##J = opt->expij[J##o+J##p][off]; \
        prij = opt->rij[J##o+J##p]; \
        envs->r##I##J[0] = prij[off*3+0]; \
        envs->r##I##J[1] = prij[off*3+1]; \
        envs->r##I##J[2] = prij[off*3+2]; \
        envs->r##I##J##rx[0] = envs->r##I##J[0] - envs->rx_in_r##I##J##rx[0]; \
        envs->r##I##J##rx[1] = envs->r##I##J[1] - envs->rx_in_r##I##J##rx[1]; \
        envs->r##I##J##rx[2] = envs->r##I##J[2] - envs->rx_in_r##I##J##rx[2];

#define PRIM2CTR(ctrsymb, gp, ngp) \
        if (ctrsymb##_ctr > 1) {\
                if (*ctrsymb##empty) { \
                        CINTprim_to_ctr_0(gctr##ctrsymb, ngp, gp, \
                                          envs->ctrsymb##_prim, \
                                          ctrsymb##_ctr, c##ctrsymb+ctrsymb##p); \
                } else { \
                        off = ctrsymb##o + ctrsymb##p; \
                        prim_to_ctr_opt(gctr##ctrsymb, ngp, gp, \
                                        opt->non0coeff[off], \
                                        opt->non0idx[off], \
                                        opt->non0ctr[off]); \
                } \
        } \
        *ctrsymb##empty = 0

#define COPY_AND_CLOSING(GCTRL, EMPTY) \
        if (n_comp > 1 && !(EMPTY)) { \
                const int INC1 = 1; \
                double *gctr1, *gctr2, *gctr3; \
                switch (n_comp) { \
                case 3: \
                        gctr1 = gctr  +envs->nf*nc; \
                        gctr2 = gctr1 +envs->nf*nc; \
                        for (n = 0, ip = 0; n < envs->nf*nc; n++, ip+=3) { \
                                gctr [n] = GCTRL[ip+0]; \
                                gctr1[n] = GCTRL[ip+1]; \
                                gctr2[n] = GCTRL[ip+2]; \
                        } \
                        break; \
                default: \
                        for (kp = 0; kp < n_comp-3; kp+=4) { \
                                gctr1 = gctr  +envs->nf*nc; \
                                gctr2 = gctr1 +envs->nf*nc; \
                                gctr3 = gctr2 +envs->nf*nc; \
                                for (n = 0, ip = kp; n < envs->nf*nc; n++,ip+=n_comp) { \
                                        gctr [n] = GCTRL[ip+0]; \
                                        gctr1[n] = GCTRL[ip+1]; \
                                        gctr2[n] = GCTRL[ip+2]; \
                                        gctr3[n] = GCTRL[ip+3]; \
                                } \
                                gctr +=envs->nf*nc * 4; \
                        } \
                        for (; kp < n_comp; kp++) { \
                                n =envs->nf*nc; \
                                dcopy_(&n, GCTRL+kp, &n_comp, gctr, &INC1); \
                                gctr +=envs->nf*nc; \
                        } \
                } \
        } \
        _mm_free(g); \
        return !(EMPTY);

// i_ctr = j_ctr = k_ctr = l_ctr = 1;
int CINT2e_1111_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;
        const int nc = 1;
        const int leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const int len0 = envs->nf * n_comp;
        const int len = leng + len0;
        double *const g = (double *)_mm_malloc(sizeof(double)*len, 16);
        double *gout;
        if (n_comp == 1) {
                gout = gctr;
        } else {
                gout = g + leng;
        }

        USE_OPT;

        for (lp = 0; lp < envs->l_prim; lp++) {
                envs->al = al[lp];
                fac1l = envs->common_factor * cl[lp];
                for (kp = 0; kp < envs->k_prim; kp++) {
                        SET_RIJ(k, l);
                        fac1k = fac1l * ck[kp];

                        for (jp = 0; jp < envs->j_prim; jp++) {
                                envs->aj = aj[jp];
                                fac1j = fac1k * cj[jp];
                                for (ip = 0; ip < envs->i_prim; ip++) {
                                        if (opt->cceij[lo+lp][ko+kp]
                                            +opt->cceij[jo+jp][io+ip]
                                            > CUTOFF15) {
                                                goto i_contracted;
                                        }
                                        SET_RIJ(i, j);
                                        fac1i = fac1j*ci[ip]*expij*expkl;
                                        CINT2e_core(gout, g, fac1i, envs, *empty);
                                        *empty = 0;
i_contracted: ;
                                } // end loop i_prim
                        } // end loop j_prim
k_contracted: ;
                } // end loop k_prim
        } // end loop l_prim

        COPY_AND_CLOSING(gout, *empty);
}

// i_ctr = n; j_ctr = k_ctr = l_ctr = 1;
int CINT2e_n111_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;

        const int nc = i_ctr;
        const int leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const int leni = envs->nf * i_ctr * n_comp; // gctri
        const int len0 = envs->nf * n_comp; // gout
        const int len = leng + leni + len0;
        double *const g = (double *)_mm_malloc(sizeof(double)*len, 16);
        double *g1 = g + leng;
        double *gout, *gctri;
        if (n_comp == 1) {
                gctri = gctr;
        } else {
                gctri = g1;
                g1 += leni;
        }
        gout = g1;

        USE_OPT;

        for (lp = 0; lp < envs->l_prim; lp++) {
                envs->al = al[lp];
                fac1l = envs->common_factor * cl[lp];
                for (kp = 0; kp < envs->k_prim; kp++) {
                        SET_RIJ(k, l);
                        fac1k = fac1l * ck[kp];

                        for (jp = 0; jp < envs->j_prim; jp++) {
                                envs->aj = aj[jp];
                                fac1j = fac1k * cj[jp];
                                for (ip = 0; ip < envs->i_prim; ip++) {
                                        if (opt->cceij[lo+lp][ko+kp]
                                            +opt->cceij[jo+jp][io+ip]
                                            > CUTOFF15) {
                                                goto i_contracted;
                                        }
                                        SET_RIJ(i, j);
                                        fac1i = fac1j*expij*expkl;
                                        CINT2e_core(gout, g, fac1i, envs, 1);
                                        PRIM2CTR(i, gout,envs->nf*n_comp);
i_contracted: ;
                                } // end loop i_prim
                        } // end loop j_prim
k_contracted: ;
                } // end loop k_prim
        } // end loop l_prim

        COPY_AND_CLOSING(gctri, *iempty);
}

// j_ctr = n; i_ctr = k_ctr = l_ctr = 1;
int CINT2e_1n11_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;

        const int nc = j_ctr;
        const int leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const int lenj = envs->nf * j_ctr * n_comp; // gctrj
        const int len0 = envs->nf * n_comp; // gout
        const int len = leng + lenj + len0;
        double *const g = (double *)_mm_malloc(sizeof(double)*len, 16);
        double *g1 = g + leng;
        double *gout, *gctrj;
        if (n_comp == 1) {
                gctrj = gctr;
        } else {
                gctrj = g1;
                g1 += lenj;
        }
        gout = g1;

        USE_OPT;

        for (lp = 0; lp < envs->l_prim; lp++) {
                envs->al = al[lp];
                fac1l = envs->common_factor * cl[lp];
                for (kp = 0; kp < envs->k_prim; kp++) {
                        SET_RIJ(k, l);
                        fac1k = fac1l * ck[kp];

                        for (jp = 0; jp < envs->j_prim; jp++) {
                                envs->aj = aj[jp];
                                fac1j = fac1k;
                                *iempty = 1;
                                for (ip = 0; ip < envs->i_prim; ip++) {
                                        if (opt->cceij[lo+lp][ko+kp]
                                            +opt->cceij[jo+jp][io+ip]
                                            > CUTOFF15) {
                                                goto i_contracted;
                                        }
                                        SET_RIJ(i, j);
                                        fac1i = fac1j*ci[ip]*expij*expkl;
                                        CINT2e_core(gout, g, fac1i, envs, *iempty);
                                        *iempty = 0;
i_contracted: ;
                                } // end loop i_prim
                                if (!*iempty) {
                                        PRIM2CTR(j, gout,envs->nf*n_comp);
                                }
                        } // end loop j_prim
k_contracted: ;
                } // end loop k_prim
        } // end loop l_prim

        COPY_AND_CLOSING(gctrj, *jempty);
}

// k_ctr = n; i_ctr = j_ctr = l_ctr = 1;
int CINT2e_11n1_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;

        const int nc = k_ctr;
        const int leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const int lenk = envs->nf * k_ctr * n_comp; // gctrk
        const int len0 = envs->nf * n_comp; // gout
        const int len = leng + lenk + len0;
        double *const g = (double *)_mm_malloc(sizeof(double)*len, 16);
        double *g1 = g + leng;
        double *gout, *gctrk;
        if (n_comp == 1) {
                gctrk = gctr;
        } else {
                gctrk = g1;
                g1 += lenk;
        }
        gout = g1;

        USE_OPT;

        for (lp = 0; lp <envs-> l_prim; lp++) {
                envs->al = al[lp];
                fac1l = envs->common_factor * cl[lp];
                for (kp = 0; kp <envs-> k_prim; kp++) {
                        SET_RIJ(k, l);
                        fac1k = fac1l;
                        *jempty = 1;
                        for (jp = 0; jp <envs-> j_prim; jp++) {
                                envs->aj = aj[jp];
                                fac1j = fac1k * cj[jp];
                                for (ip = 0; ip <envs-> i_prim; ip++) {
                                        if (opt->cceij[lo+lp][ko+kp]
                                            +opt->cceij[jo+jp][io+ip]
                                            > CUTOFF15) {
                                                goto i_contracted;
                                        }
                                        SET_RIJ(i, j);
                                        fac1i = fac1j*ci[ip]*expij*expkl;
                                        CINT2e_core(gout, g, fac1i, envs, *jempty);
                                        *jempty = 0;
i_contracted: ;
                                } // end loop i_prim
                        } // end loop j_prim
                        if (!*jempty) {
                                PRIM2CTR(k, gout,envs->nf*n_comp);
                        }
k_contracted: ;
                } // end loop k_prim
        } // end loop l_prim

        COPY_AND_CLOSING(gctrk, *kempty);
}

// l_ctr = n; i_ctr = j_ctr = k_ctr = 1;
int CINT2e_111n_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;

        const int nc = l_ctr;
        const int leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const int lenl = envs->nf * l_ctr * n_comp; // gctrl
        const int len0 = envs->nf * n_comp; // gout
        const int len = leng + lenl + len0;
        double *const g = (double *)_mm_malloc(sizeof(double)*len, 16);
        double *g1 = g + leng;
        double *gout, *gctrl;
        if (n_comp == 1) {
                gctrl = gctr;
        } else {
                gctrl = g1;
                g1 += lenl;
        }
        gout = g1;

        USE_OPT;

        for (lp = 0; lp <envs-> l_prim; lp++) {
                envs->al = al[lp];
                fac1l = envs->common_factor;
                *kempty = 1;
                for (kp = 0; kp <envs-> k_prim; kp++) {
                        SET_RIJ(k, l);
                        fac1k = fac1l * ck[kp];
                        for (jp = 0; jp <envs-> j_prim; jp++) {
                                envs->aj = aj[jp];
                                fac1j = fac1k * cj[jp];
                                for (ip = 0; ip <envs-> i_prim; ip++) {
                                        if (opt->cceij[lo+lp][ko+kp]
                                            +opt->cceij[jo+jp][io+ip]
                                            > CUTOFF15) {
                                                goto i_contracted;
                                        }
                                        SET_RIJ(i, j);
                                        fac1i = fac1j*ci[ip]*expij*expkl;
                                        CINT2e_core(gout, g, fac1i, envs, *kempty);
                                        *kempty = 0;
i_contracted: ;
                                } // end loop i_prim
                        } // end loop j_prim
k_contracted: ;
                } // end loop k_prim
                if (!*kempty) {
                        PRIM2CTR(l, gout,envs->nf*n_comp);
                }
        } // end loop l_prim

        COPY_AND_CLOSING(gctrl, *lempty);
}


int CINT2e_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;
        const int nc = i_ctr * j_ctr * k_ctr * l_ctr;
        // (irys,i,j,k,l,coord,0:1); +1 for nabla-r12
        const int leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const int lenl = envs->nf * nc * n_comp; // gctrl
        const int lenk = envs->nf * i_ctr * j_ctr * k_ctr * n_comp; // gctrk
        const int lenj = envs->nf * i_ctr * j_ctr * n_comp; // gctrj
        const int leni = envs->nf * i_ctr * n_comp; // gctri
        const int len0 = envs->nf * n_comp; // gout
        const int len = leng + lenl + lenk + lenj + leni + len0;
        double *const g = (double *)_mm_malloc(sizeof(double)*len, 16);
        double *g1 = g + leng;
        double *gout, *gctri, *gctrj, *gctrk, *gctrl;

        if (n_comp == 1) {
                gctrl = gctr;
        } else {
                gctrl = g1;
                g1 += lenl;
        }
        if (l_ctr == 1) {
                gctrk = gctrl;
                kempty = lempty;
        } else {
                gctrk = g1;
                g1 += lenk;
        }
        if (k_ctr == 1) {
                gctrj = gctrk;
                jempty = kempty;
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

        /* USE_OPT */
        int off;
        double expij, expkl, expijkl;
        double *prij;
        const int io = opt->prim_offset[i_sh];
        const int jo = opt->prim_offset[j_sh];
        const int ko = opt->prim_offset[k_sh];
        const int lo = opt->prim_offset[l_sh];
        envs->idx = opt->index_xyz_array[envs->i_l*ANG_MAX*ANG_MAX*ANG_MAX
                                        +envs->j_l*ANG_MAX*ANG_MAX
                                        +envs->k_l*ANG_MAX
                                        +envs->l_l];
        /* USE_OPT end */

        *lempty = 1;
        for (lp = 0; lp < envs->l_prim; lp++) {
                envs->al = al[lp];
                if (l_ctr == 1) {
                        fac1l = envs->common_factor * cl[lp];
                } else {
                        fac1l = envs->common_factor;
                        *kempty = 1;
                }
                for (kp = 0; kp < envs->k_prim; kp++) {
                        /* SET_RIJ(k, l); */
                        envs->ak = ak[kp];
                        envs->akl = ak[kp] + al[lp];
                        off = ko + kp;
                        if (opt->cceij[lo+lp][off] > CUTOFF15) {
                                goto k_contracted;
                        }
                        expkl = opt->expij[lo+lp][off];
                        prij = opt->rij[lo+lp];
                        envs->rkl[0] = prij[off*3+0];
                        envs->rkl[1] = prij[off*3+1];
                        envs->rkl[2] = prij[off*3+2];
                        envs->rklrx[0] = envs->rkl[0] - envs->rx_in_rklrx[0];
                        envs->rklrx[1] = envs->rkl[1] - envs->rx_in_rklrx[1];
                        envs->rklrx[2] = envs->rkl[2] - envs->rx_in_rklrx[2];
                        /* SET_RIJ(k, l); end */
                        if (k_ctr == 1) {
                                fac1k = fac1l * ck[kp];
                        } else {
                                fac1k = fac1l;
                                *jempty = 1;
                        }

                        for (jp = 0; jp < envs->j_prim; jp++) {
                                envs->aj = aj[jp];
                                if (j_ctr == 1) {
                                        fac1j = fac1k * cj[jp];
                                } else {
                                        fac1j = fac1k;
                                        *iempty = 1;
                                }
                                for (ip = 0; ip < envs->i_prim; ip++) {
                                        if (opt->cceij[lo+lp][ko+kp]
                                            +opt->cceij[jo+jp][io+ip]
                                            > CUTOFF15) {
                                                goto i_contracted;
                                        }
                                        /* SET_RIJ(i, j); */
                                        envs->ai = ai[ip];
                                        envs->aij = ai[ip] + aj[jp];
                                        off = io + ip;
                                        if (opt->cceij[jo+jp][off] > CUTOFF15) {
                                                goto i_contracted;
                                        }
                                        expij = opt->expij[jo+jp][off];
                                        prij = opt->rij[jo+jp];
                                        envs->rij[0] = prij[off*3+0];
                                        envs->rij[1] = prij[off*3+1];
                                        envs->rij[2] = prij[off*3+2];
                                        envs->rijrx[0] = envs->rij[0] - envs->rx_in_rijrx[0];
                                        envs->rijrx[1] = envs->rij[1] - envs->rx_in_rijrx[1];
                                        envs->rijrx[2] = envs->rij[2] - envs->rx_in_rijrx[2];
                                        /* SET_RIJ(i, j); end */
                                        expijkl = expij*expkl;
                                        if (i_ctr == 1) {
                                                fac1i = fac1j*ci[ip]*expijkl;
                                        } else {
                                                fac1i = fac1j*expijkl;
                                        }
                                        CINT2e_core(gout, g, fac1i, envs, *gempty);
                                        PRIM2CTR(i, gout, envs->nf*n_comp);
i_contracted: ;
                                } // end loop i_prim
                                if (!*iempty) {
                                        PRIM2CTR(j, gctri, envs->nf*i_ctr*n_comp);
                                }
                        } // end loop j_prim
                        if (!*jempty) {
                                PRIM2CTR(k, gctrj, envs->nf*i_ctr*j_ctr*n_comp);
                        }
k_contracted: ;
                } // end loop k_prim
                if (!*kempty) {
//TODO: merge this contraction with COPY_AND_CLOSING for n_comp>1
                        PRIM2CTR(l, gctrk, envs->nf*i_ctr*j_ctr*k_ctr*n_comp);
                }
        } // end loop l_prim

        COPY_AND_CLOSING(gctrl, *lempty);
}


static int (*CINTf_2e_loop[16])() = {
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_n111_loop,
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_1n11_loop,
        CINT2e_loop,
        CINT2e_11n1_loop,
        CINT2e_111n_loop,
        CINT2e_1111_loop,
};

int CINT2e_cart_drv(double *opijkl, CINTEnvVars *envs, const CINTOpt *opt)
{
        const int ip = CINTcgto_cart(envs->shls[0], envs->bas);
        const int jp = CINTcgto_cart(envs->shls[1], envs->bas);
        const int kp = CINTcgto_cart(envs->shls[2], envs->bas);
        const int lp = CINTcgto_cart(envs->shls[3], envs->bas);
        const int nop = ip * jp * kp * lp;
        const int nc = envs->nf * envs->i_ctr * envs->j_ctr
                                * envs->k_ctr * envs->l_ctr * envs->ncomp_e1;
        double *const gctr = _mm_malloc(sizeof(double) * nc * envs->ncomp_e1
                                        * envs->ncomp_tensor, 16);
        double *pgctr = gctr;
        int n;
        int has_value;

        n = ((envs->i_ctr==1) << 3) + ((envs->j_ctr==1) << 2)
          + ((envs->k_ctr==1) << 1) +  (envs->l_ctr==1);
        if (opt) {
                has_value = CINTf_2e_loop[n](gctr, envs, opt);
        } else {
                has_value = CINT2e_loop_nopt(gctr, envs);
        }

        if (has_value) {
                for (n = 0; n < envs->ncomp_tensor; n++) {
                        c2s_cart_2e1(opijkl, pgctr, envs->shls, envs->bas);
                        opijkl += nop;
                        pgctr += nc;
                }
        } else {
                CINTdset0(nop * envs->ncomp_tensor, opijkl);
        }
        _mm_free(gctr);
        return has_value;
}
int CINT2e_spheric_drv(double *opijkl, CINTEnvVars *envs, const CINTOpt *opt)
{
        const int ip = CINTcgto_spheric(envs->shls[0], envs->bas);
        const int jp = CINTcgto_spheric(envs->shls[1], envs->bas);
        const int kp = CINTcgto_spheric(envs->shls[2], envs->bas);
        const int lp = CINTcgto_spheric(envs->shls[3], envs->bas);
        const int nop = ip * jp * kp * lp;
        const int nc = envs->nf * envs->i_ctr * envs->j_ctr
                                * envs->k_ctr * envs->l_ctr * envs->ncomp_e1;
        double *const gctr = _mm_malloc(sizeof(double) * nc * envs->ncomp_e2
                                        * envs->ncomp_tensor, 16);
        double *pgctr = gctr;
        int n;
        int has_value;

        n = ((envs->i_ctr==1) << 3) + ((envs->j_ctr==1) << 2)
          + ((envs->k_ctr==1) << 1) +  (envs->l_ctr==1);
        if (opt) {
                has_value = CINTf_2e_loop[n](gctr, envs, opt);
        } else {
                has_value = CINT2e_loop_nopt(gctr, envs);
        }

        if (has_value) {
                for (n = 0; n < envs->ncomp_tensor; n++) {
                        c2s_sph_2e1(opijkl, pgctr, envs->shls, envs->bas);
                        opijkl += nop;
                        pgctr += nc;
                }
        } else {
                CINTdset0(nop * envs->ncomp_tensor, opijkl);
        }
        _mm_free(gctr);
        return has_value;
}
int CINT2e_spinor_drv(double *opijkl, CINTEnvVars *envs, const CINTOpt *opt,
                      void (*const f_e1_c2s)(), void (*const f_e2_c2s)())
{
        const int ip = CINTcgto_spinor(envs->shls[0], envs->bas);
        const int jp = CINTcgto_spinor(envs->shls[1], envs->bas);
        const int kp = CINTcgto_spinor(envs->shls[2], envs->bas);
        const int lp = CINTcgto_spinor(envs->shls[3], envs->bas);
        const int nop = ip * jp * kp * lp;
        const int nc = envs->nf * envs->i_ctr * envs->j_ctr
                                * envs->k_ctr * envs->l_ctr * envs->ncomp_e1;
        const int n1 = ip * envs->nfk * envs->k_ctr
                * envs->nfl * envs->l_ctr * jp * OF_CMPLX;
        const int len1 = (nc*envs->ncomp_e2*envs->ncomp_tensor+16)&0xfffffff0;
        double *gctr = _mm_malloc(sizeof(double)*(len1+n1*envs->ncomp_e2), 16);
        int n, m;
        int has_value;

        n = ((envs->i_ctr==1) << 3) + ((envs->j_ctr==1) << 2)
          + ((envs->k_ctr==1) << 1) +  (envs->l_ctr==1);
        if (opt) {
                has_value = CINTf_2e_loop[n](gctr, envs, opt);
        } else {
                has_value = CINT2e_loop_nopt(gctr, envs);
        }

        if (has_value) {
                double *pgctr = gctr;
                double *opij = gctr + len1;
                for (n = 0; n < envs->ncomp_tensor; n++) {
                        for (m = 0; m < envs->ncomp_e2; m++) {
                                (*f_e1_c2s)(opij+n1*m, pgctr, envs->shls,
                                            envs->bas);
                                pgctr += nc;
                        }
                        (*f_e2_c2s)(opijkl, opij, envs->shls, envs->bas);
                        opijkl += nop * OF_CMPLX;
                }
        } else {
                CINTdset0(nop * OF_CMPLX * envs->ncomp_tensor, opijkl);
        }
        _mm_free(gctr);
        return has_value;
}


static void r1_gout2e(double *g, double *gout, const int *idx,
                      const CINTEnvVars *envs, int gout_empty)
{
        int ix, iy, iz, n;

        if (gout_empty) {
                for (n = 0; n < envs->nf; n++, idx+=3) {
                        ix = idx[0];
                        iy = idx[1];
                        iz = idx[2];
                        gout[n] = g[ix] * g[iy] * g[iz];
                }
        } else {
                for (n = 0; n < envs->nf; n++, idx+=3) {
                        ix = idx[0];
                        iy = idx[1];
                        iz = idx[2];
                        gout[n] += g[ix] * g[iy] * g[iz];
                }
        }
}
static void r2_gout2e(double *g, double *gout, const int *idx,
                      const CINTEnvVars *envs, int gout_empty)
{
        int ix, iy, iz, n;
        __m128d r0, r1;
        double s;

        if (gout_empty) {
                for (n = 0; n < envs->nf; n++, idx+=3) {
                        ix = idx[0];
                        iy = idx[1];
                        iz = idx[2];
                        r0 = _mm_load_pd(g+ix  );
                        r1 = _mm_load_pd(g+iy  );
                        r0 = _mm_mul_pd(r0, r1);
                        r1 = _mm_load_pd(g+iz  );
                        r0 = _mm_mul_pd(r0, r1);
                        r0 = _mm_hadd_pd(r0, r0);
                        _mm_store_sd(gout+n, r0);
                }
        } else {
                for (n = 0; n < envs->nf; n++, idx+=3) {
                        ix = idx[0];
                        iy = idx[1];
                        iz = idx[2];
                        r0 = _mm_load_pd(g+ix  );
                        r1 = _mm_load_pd(g+iy  );
                        r0 = _mm_mul_pd(r0, r1);
                        r1 = _mm_load_pd(g+iz  );
                        r0 = _mm_mul_pd(r0, r1);
                        r0 = _mm_hadd_pd(r0, r0);
                        _mm_store_sd(&s, r0);
                        gout[n] += s;
                }
        }
}
static void r3_gout2e(double *g, double *gout, const int *idx,
                      const CINTEnvVars *envs, int gout_empty)
{
        int ix, iy, iz, n;
        __m128d r0, r1;
        double s;

        if (gout_empty) {
                for (n = 0; n < envs->nf; n++, idx+=3) {
                        ix = idx[0];
                        iy = idx[1];
                        iz = idx[2];
                        r0 = _mm_load_pd(g+ix  );
                        r1 = _mm_load_pd(g+iy  );
                        r0 = _mm_mul_pd(r0, r1);
                        r1 = _mm_load_pd(g+iz  );
                        r0 = _mm_mul_pd(r0, r1);
                        r0 = _mm_hadd_pd(r0, r0);
                        _mm_store_sd(gout+n, r0);
                        gout[n] += g[ix+2] * g[iy+2] * g[iz+2];
                }
        } else {
                for (n = 0; n < envs->nf; n++, idx+=3) {
                        ix = idx[0];
                        iy = idx[1];
                        iz = idx[2];
                        r0 = _mm_load_pd(g+ix  );
                        r1 = _mm_load_pd(g+iy  );
                        r0 = _mm_mul_pd(r0, r1);
                        r1 = _mm_load_pd(g+iz  );
                        r0 = _mm_mul_pd(r0, r1);
                        r0 = _mm_hadd_pd(r0, r0);
                        _mm_store_sd(&s, r0);
                        gout[n] += s + g[ix+2] * g[iy+2] * g[iz+2];
                }
        }
}
static void r4_gout2e(double *g, double *gout, const int *idx,
                      const CINTEnvVars *envs, int gout_empty)
{
        int ix, iy, iz, n;
        __m128d r0, r1, r2;
        double s;

        if (gout_empty) {
                for (n = 0; n < envs->nf; n++, idx+=3) {
                        ix = idx[0];
                        iy = idx[1];
                        iz = idx[2];
                        r0 = _mm_load_pd(g+ix  );
                        r1 = _mm_load_pd(g+iy  );
                        r0 = _mm_mul_pd(r0, r1);
                        r1 = _mm_load_pd(g+iz  );
                        r0 = _mm_mul_pd(r0, r1);
                        r2 = _mm_load_pd(g+ix+2);
                        r1 = _mm_load_pd(g+iy+2);
                        r2 = _mm_mul_pd(r2, r1);
                        r1 = _mm_load_pd(g+iz+2);
                        r2 = _mm_mul_pd(r2, r1);
                        r0 = _mm_add_pd(r0, r2);
                        r0 = _mm_hadd_pd(r0, r0);
                        _mm_store_sd(gout+n, r0);
                }
        } else {
                for (n = 0; n < envs->nf; n++, idx+=3) {
                        ix = idx[0];
                        iy = idx[1];
                        iz = idx[2];
                        r0 = _mm_load_pd(g+ix  );
                        r1 = _mm_load_pd(g+iy  );
                        r0 = _mm_mul_pd(r0, r1);
                        r1 = _mm_load_pd(g+iz  );
                        r0 = _mm_mul_pd(r0, r1);
                        r2 = _mm_load_pd(g+ix+2);
                        r1 = _mm_load_pd(g+iy+2);
                        r2 = _mm_mul_pd(r2, r1);
                        r1 = _mm_load_pd(g+iz+2);
                        r2 = _mm_mul_pd(r2, r1);
                        r0 = _mm_add_pd(r0, r2);
                        r0 = _mm_hadd_pd(r0, r0);
                        _mm_store_sd(&s, r0);
                        gout[n] += s;
                }
        }
}
static void r5_gout2e(double *g, double *gout, const int *idx,
                      const CINTEnvVars *envs, int gout_empty)
{
        int ix, iy, iz, n;
        __m128d r0, r1, r2;
        double s;

        if (gout_empty) {
                for (n = 0; n < envs->nf; n++, idx+=3) {
                        ix = idx[0];
                        iy = idx[1];
                        iz = idx[2];
                        r0 = _mm_load_pd(g+ix  );
                        r1 = _mm_load_pd(g+iy  );
                        r0 = _mm_mul_pd(r0, r1);
                        r1 = _mm_load_pd(g+iz  );
                        r0 = _mm_mul_pd(r0, r1);
                        r2 = _mm_load_pd(g+ix+2);
                        r1 = _mm_load_pd(g+iy+2);
                        r2 = _mm_mul_pd(r2, r1);
                        r1 = _mm_load_pd(g+iz+2);
                        r2 = _mm_mul_pd(r2, r1);
                        r0 = _mm_add_pd(r0, r2);
                        r0 = _mm_hadd_pd(r0, r0);
                        _mm_store_sd(gout+n, r0);
                        gout[n] += g[ix+4] * g[iy+4] * g[iz+4];
                }
        } else {
                for (n = 0; n < envs->nf; n++, idx+=3) {
                        ix = idx[0];
                        iy = idx[1];
                        iz = idx[2];
                        r0 = _mm_load_pd(g+ix  );
                        r1 = _mm_load_pd(g+iy  );
                        r0 = _mm_mul_pd(r0, r1);
                        r1 = _mm_load_pd(g+iz  );
                        r0 = _mm_mul_pd(r0, r1);
                        r2 = _mm_load_pd(g+ix+2);
                        r1 = _mm_load_pd(g+iy+2);
                        r2 = _mm_mul_pd(r2, r1);
                        r1 = _mm_load_pd(g+iz+2);
                        r2 = _mm_mul_pd(r2, r1);
                        r0 = _mm_add_pd(r0, r2);
                        r0 = _mm_hadd_pd(r0, r0);
                        _mm_store_sd(&s, r0);
                        gout[n] += s + g[ix+4] * g[iy+4] * g[iz+4];
                }
        }
}
static void r6_gout2e(double *g, double *gout, const int *idx,
                      const CINTEnvVars *envs, int gout_empty)
{
        int ix, iy, iz, n;
        __m128d r0, r1, r2;
        double s;

        if (gout_empty) {
                for (n = 0; n < envs->nf; n++, idx+=3) {
                        ix = idx[0];
                        iy = idx[1];
                        iz = idx[2];
                        r0 = _mm_load_pd(g+ix  );
                        r1 = _mm_load_pd(g+iy  );
                        r0 = _mm_mul_pd(r0, r1);
                        r1 = _mm_load_pd(g+iz  );
                        r0 = _mm_mul_pd(r0, r1);
                        r2 = _mm_load_pd(g+ix+2);
                        r1 = _mm_load_pd(g+iy+2);
                        r2 = _mm_mul_pd(r2, r1);
                        r1 = _mm_load_pd(g+iz+2);
                        r2 = _mm_mul_pd(r2, r1);
                        r0 = _mm_add_pd(r0, r2);
                        r2 = _mm_load_pd(g+ix+4);
                        r1 = _mm_load_pd(g+iy+4);
                        r2 = _mm_mul_pd(r2, r1);
                        r1 = _mm_load_pd(g+iz+4);
                        r2 = _mm_mul_pd(r2, r1);
                        r0 = _mm_add_pd(r0, r2);
                        r0 = _mm_hadd_pd(r0, r0);
                        _mm_store_sd(gout+n, r0);
                }
        } else {
                for (n = 0; n < envs->nf; n++, idx+=3) {
                        ix = idx[0];
                        iy = idx[1];
                        iz = idx[2];
                        r0 = _mm_load_pd(g+ix  );
                        r1 = _mm_load_pd(g+iy  );
                        r0 = _mm_mul_pd(r0, r1);
                        r1 = _mm_load_pd(g+iz  );
                        r0 = _mm_mul_pd(r0, r1);
                        r2 = _mm_load_pd(g+ix+2);
                        r1 = _mm_load_pd(g+iy+2);
                        r2 = _mm_mul_pd(r2, r1);
                        r1 = _mm_load_pd(g+iz+2);
                        r2 = _mm_mul_pd(r2, r1);
                        r0 = _mm_add_pd(r0, r2);
                        r2 = _mm_load_pd(g+ix+4);
                        r1 = _mm_load_pd(g+iy+4);
                        r2 = _mm_mul_pd(r2, r1);
                        r1 = _mm_load_pd(g+iz+4);
                        r2 = _mm_mul_pd(r2, r1);
                        r0 = _mm_add_pd(r0, r2);
                        r0 = _mm_hadd_pd(r0, r0);
                        _mm_store_sd(&s, r0);
                        gout[n] += s;
                }
        }
}
static void r7_gout2e(double *g, double *gout, const int *idx,
                      const CINTEnvVars *envs, int gout_empty)
{
        int ix, iy, iz, n;
        __m128d r0, r1, r2;
        double s;

        if (gout_empty) {
                for (n = 0; n < envs->nf; n++, idx+=3) {
                        ix = idx[0];
                        iy = idx[1];
                        iz = idx[2];
                        r0 = _mm_load_pd(g+ix  );
                        r1 = _mm_load_pd(g+iy  );
                        r0 = _mm_mul_pd(r0, r1);
                        r1 = _mm_load_pd(g+iz  );
                        r0 = _mm_mul_pd(r0, r1);
                        r2 = _mm_load_pd(g+ix+2);
                        r1 = _mm_load_pd(g+iy+2);
                        r2 = _mm_mul_pd(r2, r1);
                        r1 = _mm_load_pd(g+iz+2);
                        r2 = _mm_mul_pd(r2, r1);
                        r0 = _mm_add_pd(r0, r2);
                        r2 = _mm_load_pd(g+ix+4);
                        r1 = _mm_load_pd(g+iy+4);
                        r2 = _mm_mul_pd(r2, r1);
                        r1 = _mm_load_pd(g+iz+4);
                        r2 = _mm_mul_pd(r2, r1);
                        r0 = _mm_add_pd(r0, r2);
                        r0 = _mm_hadd_pd(r0, r0);
                        _mm_store_sd(gout+n, r0);
                        gout[n] += g[ix+6] * g[iy+6] * g[iz+6];
                }
        } else {
                for (n = 0; n < envs->nf; n++, idx+=3) {
                        ix = idx[0];
                        iy = idx[1];
                        iz = idx[2];
                        r0 = _mm_load_pd(g+ix  );
                        r1 = _mm_load_pd(g+iy  );
                        r0 = _mm_mul_pd(r0, r1);
                        r1 = _mm_load_pd(g+iz  );
                        r0 = _mm_mul_pd(r0, r1);
                        r2 = _mm_load_pd(g+ix+2);
                        r1 = _mm_load_pd(g+iy+2);
                        r2 = _mm_mul_pd(r2, r1);
                        r1 = _mm_load_pd(g+iz+2);
                        r2 = _mm_mul_pd(r2, r1);
                        r0 = _mm_add_pd(r0, r2);
                        r2 = _mm_load_pd(g+ix+4);
                        r1 = _mm_load_pd(g+iy+4);
                        r2 = _mm_mul_pd(r2, r1);
                        r1 = _mm_load_pd(g+iz+4);
                        r2 = _mm_mul_pd(r2, r1);
                        r0 = _mm_add_pd(r0, r2);
                        r0 = _mm_hadd_pd(r0, r0);
                        _mm_store_sd(&s, r0);
                        gout[n] += s + g[ix+6] * g[iy+6] * g[iz+6];
                }
        }
}
static void ra_gout2e(double *g, double *gout, const int *idx,
                      const CINTEnvVars *envs, int gout_empty)
{
        int i, ix, iy, iz, n;
        __m128d r0, r1, r2, r3;
        double s;

        if (gout_empty) {
                for (n = 0; n < envs->nf; n++, idx+=3) {
                        ix = idx[0];
                        iy = idx[1];
                        iz = idx[2];
                        r0 = _mm_load_pd(g+ix  );
                        r1 = _mm_load_pd(g+iy  );
                        r0 = _mm_mul_pd(r0, r1);
                        r1 = _mm_load_pd(g+iz  );
                        r0 = _mm_mul_pd(r0, r1);
                        r2 = _mm_load_pd(g+ix+2);
                        r1 = _mm_load_pd(g+iy+2);
                        r2 = _mm_mul_pd(r2, r1);
                        r1 = _mm_load_pd(g+iz+2);
                        r2 = _mm_mul_pd(r2, r1);
                        r0 = _mm_add_pd(r0, r2);
                        r2 = _mm_load_pd(g+ix+4);
                        r1 = _mm_load_pd(g+iy+4);
                        r2 = _mm_mul_pd(r2, r1);
                        r1 = _mm_load_pd(g+iz+4);
                        r2 = _mm_mul_pd(r2, r1);
                        r0 = _mm_add_pd(r0, r2);
                        r2 = _mm_load_pd(g+ix+6);
                        r1 = _mm_load_pd(g+iy+6);
                        r2 = _mm_mul_pd(r2, r1);
                        r1 = _mm_load_pd(g+iz+6);
                        r2 = _mm_mul_pd(r2, r1);
                        r3 = _mm_add_pd(r0, r2);
                        for (i = 8; i < envs->nrys_roots-1; i+=2) {
                                r0 = _mm_load_pd(g+ix+i);
                                r1 = _mm_load_pd(g+iy+i);
                                r0 = _mm_mul_pd(r0, r1);
                                r1 = _mm_load_pd(g+iz+i);
                                r0 = _mm_mul_pd(r0, r1);
                                r3 = _mm_add_pd(r3, r0);
                        }
                        r3 = _mm_hadd_pd(r3, r3);
                        _mm_store_sd(gout+n, r3);
                        if (i < envs->nrys_roots) {
                                gout[n] += g[ix+i] * g[iy+i] * g[iz+i];
                        }
                }
        } else {
                for (n = 0; n < envs->nf; n++, idx+=3) {
                        ix = idx[0];
                        iy = idx[1];
                        iz = idx[2];
                        r0 = _mm_load_pd(g+ix  );
                        r1 = _mm_load_pd(g+iy  );
                        r0 = _mm_mul_pd(r0, r1);
                        r1 = _mm_load_pd(g+iz  );
                        r0 = _mm_mul_pd(r0, r1);
                        r2 = _mm_load_pd(g+ix+2);
                        r1 = _mm_load_pd(g+iy+2);
                        r2 = _mm_mul_pd(r2, r1);
                        r1 = _mm_load_pd(g+iz+2);
                        r2 = _mm_mul_pd(r2, r1);
                        r0 = _mm_add_pd(r0, r2);
                        r2 = _mm_load_pd(g+ix+4);
                        r1 = _mm_load_pd(g+iy+4);
                        r2 = _mm_mul_pd(r2, r1);
                        r1 = _mm_load_pd(g+iz+4);
                        r2 = _mm_mul_pd(r2, r1);
                        r0 = _mm_add_pd(r0, r2);
                        r2 = _mm_load_pd(g+ix+6);
                        r1 = _mm_load_pd(g+iy+6);
                        r2 = _mm_mul_pd(r2, r1);
                        r1 = _mm_load_pd(g+iz+6);
                        r2 = _mm_mul_pd(r2, r1);
                        r3 = _mm_add_pd(r0, r2);
                        for (i = 8; i < envs->nrys_roots-1; i+=2) {
                                r0 = _mm_load_pd(g+ix+i);
                                r1 = _mm_load_pd(g+iy+i);
                                r0 = _mm_mul_pd(r0, r1);
                                r1 = _mm_load_pd(g+iz+i);
                                r0 = _mm_mul_pd(r0, r1);
                                r3 = _mm_add_pd(r3, r0);
                        }
                        r3 = _mm_hadd_pd(r3, r3);
                        _mm_store_sd(&s, r3);
                        if (i < envs->nrys_roots) {
                                gout[n] += s + g[ix+i]*g[iy+i]*g[iz+i];
                        } else {
                                gout[n] += s;
                        }
                }
        }
}
static void (*fgout2e[MXRYSROOTS])() = {
        r1_gout2e,
        r2_gout2e,
        r3_gout2e,
        r4_gout2e,
        r5_gout2e,
        r6_gout2e,
        r7_gout2e,
        ra_gout2e,
        ra_gout2e,
        ra_gout2e,
        ra_gout2e,
        ra_gout2e,
        ra_gout2e,
        ra_gout2e,
        ra_gout2e,
        ra_gout2e,
};


int cint2e_sph(double *opijkl, const int *shls,
               const int *atm, const int natm,
               const int *bas, const int nbas, const double *env,
               const CINTOpt *opt)
{
        int ng[] = {0, 0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        //envs.f_gout = fgout2e[envs.nrys_roots-1];
        return CINT2e_spheric_drv(opijkl, &envs, opt);
}
void cint2e_sph_optimizer(CINTOpt **opt, const int *atm, const int natm,
                          const int *bas, const int nbas, const double *env)
{
        int ng[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        CINTuse_all_optimizer(opt, ng, atm, natm, bas, nbas, env);
}

int cint2e_cart(double *opijkl, const int *shls,
                const int *atm, const int natm,
                const int *bas, const int nbas, const double *env,
                const CINTOpt *opt)
{
        int ng[] = {0, 0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        //envs.f_gout = fgout2e[envs.nrys_roots-1];
        return CINT2e_cart_drv(opijkl, &envs, opt);
}
void cint2e_cart_optimizer(CINTOpt **opt, const int *atm, const int natm,
                           const int *bas, const int nbas, const double *env)
{
        cint2e_sph_optimizer(opt, atm, natm, bas, nbas, env);
}


/*
 * spinor <ki|jl> = (ij|kl); i,j\in electron 1; k,l\in electron 2
 */
int cint2e(double *opijkl, const int *shls,
           const int *atm, const int natm,
           const int *bas, const int nbas, const double *env,
           const CINTOpt *opt)
{
        int ng[] = {0, 0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        //envs.f_gout = fgout2e[envs.nrys_roots-1];
        return CINT2e_spinor_drv(opijkl, &envs, opt, &c2s_sf_2e1, &c2s_sf_2e2);
}
void cint2e_optimizer(CINTOpt **opt, const int *atm, const int natm,
                      const int *bas, const int nbas, const double *env)
{
        cint2e_sph_optimizer(opt, atm, natm, bas, nbas, env);
}


/*
 * * * * * * * * * * * * * * * * * * * * *
 * c to fortran interface
 */

C2Fo_(cint2e_cart);
C2Fo_(cint2e_sph);
C2Fo_(cint2e);
OPTIMIZER2F_(cint2e_cart_optimizer);
OPTIMIZER2F_(cint2e_sph_optimizer);
OPTIMIZER2F_(cint2e_optimizer);


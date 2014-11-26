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


#include <pmmintrin.h>
#include <mm_malloc.h>
#include "g2e.h"

/*
 * <ki|jl> = (ij|kl); i,j\in electron 1; k,l\in electron 2
void CINTgout2e(double *g, double *gout, const int *idx,
                const CINTEnvVars *envs, int gout_empty)
{
        int i, ix, iy, iz, n;

        if (gout_empty) {
                switch (envs->nrys_roots) {
                        case 1:
                                for (n = 0; n < envs->nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix] * g[iy] * g[iz];
                                }
                                break;
                        case 2:
                                for (n = 0; n < envs->nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1];
                                }
                                break;
                        case 3:
                                for (n = 0; n < envs->nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2];
                                }
                                break;
                        case 4:
                                for (n = 0; n < envs->nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3];
                                }
                                break;
                        case 5:
                                for (n = 0; n < envs->nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3]
                                                + g[ix+4] * g[iy+4] * g[iz+4];
                                }
                                break;
                        case 6:
                                for (n = 0; n < envs->nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3]
                                                + g[ix+4] * g[iy+4] * g[iz+4]
                                                + g[ix+5] * g[iy+5] * g[iz+5];
                                }
                                break;
                        case 7:
                                for (n = 0; n < envs->nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3]
                                                + g[ix+4] * g[iy+4] * g[iz+4]
                                                + g[ix+5] * g[iy+5] * g[iz+5]
                                                + g[ix+6] * g[iy+6] * g[iz+6];
                                }
                                break;
                        case 8:
                                for (n = 0; n < envs->nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3]
                                                + g[ix+4] * g[iy+4] * g[iz+4]
                                                + g[ix+5] * g[iy+5] * g[iz+5]
                                                + g[ix+6] * g[iy+6] * g[iz+6]
                                                + g[ix+7] * g[iy+7] * g[iz+7];
                                }
                                break;
                        default:
                                for (n = 0; n < envs->nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = 0;
                                        for (i = 0; i < envs->nrys_roots; i++)
                                                gout[n] += g[ix+i] * g[iy+i] * g[iz+i];
                                }
                                break;
                } // end switch nroots
        } else { // not flag_acc
                switch (envs->nrys_roots) {
                        case 1:
                                for (n = 0; n < envs->nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] += g[ix] * g[iy] * g[iz];
                                }
                                break;
                        case 2:
                                for (n = 0; n < envs->nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] +=g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1];
                                }
                                break;
                        case 3:
                                for (n = 0; n < envs->nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] +=g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2];
                                }
                                break;
                        case 4:
                                for (n = 0; n < envs->nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] +=g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3];
                                }
                                break;
                        case 5:
                                for (n = 0; n < envs->nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] +=g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3]
                                                + g[ix+4] * g[iy+4] * g[iz+4];
                                }
                                break;
                        case 6:
                                for (n = 0; n < envs->nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] +=g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3]
                                                + g[ix+4] * g[iy+4] * g[iz+4]
                                                + g[ix+5] * g[iy+5] * g[iz+5];
                                }
                                break;
                        case 7:
                                for (n = 0; n < envs->nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] +=g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3]
                                                + g[ix+4] * g[iy+4] * g[iz+4]
                                                + g[ix+5] * g[iy+5] * g[iz+5]
                                                + g[ix+6] * g[iy+6] * g[iz+6];
                                }
                                break;
                        case 8:
                                for (n = 0; n < envs->nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] +=g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3]
                                                + g[ix+4] * g[iy+4] * g[iz+4]
                                                + g[ix+5] * g[iy+5] * g[iz+5]
                                                + g[ix+6] * g[iy+6] * g[iz+6]
                                                + g[ix+7] * g[iy+7] * g[iz+7];
                                }
                                break;
                        default:
                                for (n = 0; n < envs->nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        for (i = 0; i < envs->nrys_roots; i++)
                                                gout[n] += g[ix+i] * g[iy+i] * g[iz+i];
                                }
                                break;
                } // end switch nroots
        }
}
 */

void CINTgout2e(double *g, double *gout, const int *idx,
                const CINTEnvVars *envs, int gout_empty)
{
        int i, ix, iy, iz, n;
        int jx, jy, jz;
        __m128d r0, r1, r2, r3;
        double s;

        if (gout_empty) {
                switch (envs->nrys_roots) {
                        case 1:
                                for (n = 0; n < envs->nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix] * g[iy] * g[iz];
                                }
                                break;
                        case 2:
                                for (n = 0; n < envs->nf-1; n+=2, idx+=6) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        jx = idx[3];
                                        jy = idx[4];
                                        jz = idx[5];
                                        r0 = _mm_load_pd(g+ix  );
                                        r1 = _mm_load_pd(g+iy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+iz  );
                                        r3 = _mm_mul_pd (r0, r1);
                                        r0 = _mm_load_pd(g+jx  );
                                        r1 = _mm_load_pd(g+jy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+jz  );
                                        r0 = _mm_mul_pd (r0, r1);

                                        r3 = _mm_hadd_pd(r3, r0);
                                        _mm_storeu_pd(gout+n, r3);
                                }
                                if (n < envs->nf) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        r0 = _mm_load_pd(g+ix  );
                                        r1 = _mm_load_pd(g+iy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+iz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r0 = _mm_hadd_pd(r0, r0);
                                        _mm_store_sd(gout+n, r0);
                                }
                                break;
                        case 3:
                                for (n = 0; n < envs->nf-1; n+=2, idx+=6) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        jx = idx[3];
                                        jy = idx[4];
                                        jz = idx[5];
                                        r0 = _mm_load_pd(g+ix  );
                                        r1 = _mm_load_pd(g+iy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+iz  );
                                        r3 = _mm_mul_pd (r0, r1);

                                        r0 = _mm_load_pd(g+jx  );
                                        r1 = _mm_load_pd(g+jy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+jz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r3 = _mm_hadd_pd(r3, r0);

                                        r0 = _mm_loadl_pd(r0, g+ix+2);
                                        r0 = _mm_loadh_pd(r0, g+jx+2);
                                        r1 = _mm_loadl_pd(r1, g+iy+2);
                                        r1 = _mm_loadh_pd(r1, g+jy+2);
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_loadl_pd(r1, g+iz+2);
                                        r1 = _mm_loadh_pd(r1, g+jz+2);
                                        r0 = _mm_mul_pd (r0, r1);
                                        r3 = _mm_add_pd (r0, r3);
                                        _mm_storeu_pd(gout+n, r3);
                                }
                                if (n < envs->nf) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        r0 = _mm_load_pd(g+ix  );
                                        r1 = _mm_load_pd(g+iy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+iz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r0 = _mm_hadd_pd(r0, r0);
                                        _mm_store_sd(gout+n, r0);
                                        gout[n] += g[ix+2] * g[iy+2] * g[iz+2];
                                }
                                break;
                        case 4:
                                for (n = 0; n < envs->nf-1; n+=2, idx+=6) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        jx = idx[3];
                                        jy = idx[4];
                                        jz = idx[5];
                                        r0 = _mm_load_pd(g+ix  );
                                        r1 = _mm_load_pd(g+iy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+iz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r2 = _mm_load_pd(g+ix+2);
                                        r1 = _mm_load_pd(g+iy+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+iz+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r3 = _mm_add_pd (r0, r2);

                                        r0 = _mm_load_pd(g+jx  );
                                        r1 = _mm_load_pd(g+jy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+jz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r2 = _mm_load_pd(g+jx+2);
                                        r1 = _mm_load_pd(g+jy+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+jz+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r0 = _mm_add_pd (r0, r2);

                                        r3 = _mm_hadd_pd(r3, r0);
                                        _mm_storeu_pd(gout+n, r3);
                                }
                                if (n < envs->nf) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        r0 = _mm_load_pd(g+ix  );
                                        r1 = _mm_load_pd(g+iy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+iz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r2 = _mm_load_pd(g+ix+2);
                                        r1 = _mm_load_pd(g+iy+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+iz+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r0 = _mm_add_pd (r0, r2);
                                        r0 = _mm_hadd_pd(r0, r0);
                                        _mm_store_sd(gout+n, r0);
                                }
                                break;
                        case 5:
                                for (n = 0; n < envs->nf-1; n+=2, idx+=6) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        jx = idx[3];
                                        jy = idx[4];
                                        jz = idx[5];
                                        r0 = _mm_load_pd(g+ix  );
                                        r1 = _mm_load_pd(g+iy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+iz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r2 = _mm_load_pd(g+ix+2);
                                        r1 = _mm_load_pd(g+iy+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+iz+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r3 = _mm_add_pd (r0, r2);

                                        r0 = _mm_load_pd(g+jx  );
                                        r1 = _mm_load_pd(g+jy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+jz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r2 = _mm_load_pd(g+jx+2);
                                        r1 = _mm_load_pd(g+jy+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+jz+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r0 = _mm_add_pd (r0, r2);
                                        r3 = _mm_hadd_pd(r3, r0);

                                        r0 = _mm_loadl_pd(r0, g+ix+4);
                                        r0 = _mm_loadh_pd(r0, g+jx+4);
                                        r1 = _mm_loadl_pd(r1, g+iy+4);
                                        r1 = _mm_loadh_pd(r1, g+jy+4);
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_loadl_pd(r1, g+iz+4);
                                        r1 = _mm_loadh_pd(r1, g+jz+4);
                                        r0 = _mm_mul_pd (r0, r1);
                                        r3 = _mm_add_pd (r0, r3);
                                        _mm_storeu_pd(gout+n, r3);
                                }
                                if (n < envs->nf) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        r0 = _mm_load_pd(g+ix  );
                                        r1 = _mm_load_pd(g+iy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+iz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r2 = _mm_load_pd(g+ix+2);
                                        r1 = _mm_load_pd(g+iy+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+iz+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r0 = _mm_add_pd (r0, r2);
                                        r0 = _mm_hadd_pd(r0, r0);
                                        _mm_store_sd(gout+n, r0);
                                        gout[n] += g[ix+4] * g[iy+4] * g[iz+4];
                                }
                                break;
                        case 6:
                                for (n = 0; n < envs->nf-1; n+=2, idx+=6) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        jx = idx[3];
                                        jy = idx[4];
                                        jz = idx[5];
                                        r0 = _mm_load_pd(g+ix  );
                                        r1 = _mm_load_pd(g+iy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+iz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r2 = _mm_load_pd(g+ix+2);
                                        r1 = _mm_load_pd(g+iy+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+iz+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r0 = _mm_add_pd (r0, r2);
                                        r2 = _mm_load_pd(g+ix+4);
                                        r1 = _mm_load_pd(g+iy+4);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+iz+4);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r3 = _mm_add_pd (r0, r2);

                                        r0 = _mm_load_pd(g+jx  );
                                        r1 = _mm_load_pd(g+jy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+jz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r2 = _mm_load_pd(g+jx+2);
                                        r1 = _mm_load_pd(g+jy+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+jz+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r0 = _mm_add_pd (r0, r2);
                                        r2 = _mm_load_pd(g+jx+4);
                                        r1 = _mm_load_pd(g+jy+4);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+jz+4);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r0 = _mm_add_pd (r0, r2);

                                        r3 = _mm_hadd_pd(r3, r0);
                                        _mm_storeu_pd(gout+n, r3);
                                }
                                if (n < envs->nf) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        r0 = _mm_load_pd(g+ix  );
                                        r1 = _mm_load_pd(g+iy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+iz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r2 = _mm_load_pd(g+ix+2);
                                        r1 = _mm_load_pd(g+iy+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+iz+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r0 = _mm_add_pd (r0, r2);
                                        r2 = _mm_load_pd(g+ix+4);
                                        r1 = _mm_load_pd(g+iy+4);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+iz+4);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r0 = _mm_add_pd (r0, r2);
                                        r0 = _mm_hadd_pd(r0, r0);
                                        _mm_store_sd(gout+n, r0);
                                }
                                break;
                        case 7:
                                for (n = 0; n < envs->nf-1; n+=2, idx+=6) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        jx = idx[3];
                                        jy = idx[4];
                                        jz = idx[5];
                                        r0 = _mm_load_pd(g+ix  );
                                        r1 = _mm_load_pd(g+iy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+iz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r2 = _mm_load_pd(g+ix+2);
                                        r1 = _mm_load_pd(g+iy+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+iz+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r0 = _mm_add_pd (r0, r2);
                                        r2 = _mm_load_pd(g+ix+4);
                                        r1 = _mm_load_pd(g+iy+4);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+iz+4);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r3 = _mm_add_pd (r0, r2);

                                        r0 = _mm_load_pd(g+jx  );
                                        r1 = _mm_load_pd(g+jy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+jz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r2 = _mm_load_pd(g+jx+2);
                                        r1 = _mm_load_pd(g+jy+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+jz+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r0 = _mm_add_pd (r0, r2);
                                        r2 = _mm_load_pd(g+jx+4);
                                        r1 = _mm_load_pd(g+jy+4);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+jz+4);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r0 = _mm_add_pd (r0, r2);
                                        r3 = _mm_hadd_pd(r3, r0);

                                        r0 = _mm_loadl_pd(r0, g+ix+6);
                                        r0 = _mm_loadh_pd(r0, g+jx+6);
                                        r1 = _mm_loadl_pd(r1, g+iy+6);
                                        r1 = _mm_loadh_pd(r1, g+jy+6);
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_loadl_pd(r1, g+iz+6);
                                        r1 = _mm_loadh_pd(r1, g+jz+6);
                                        r0 = _mm_mul_pd (r0, r1);
                                        r3 = _mm_add_pd (r0, r3);
                                        _mm_storeu_pd(gout+n, r3);
                                }
                                if (n < envs->nf) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        r0 = _mm_load_pd(g+ix  );
                                        r1 = _mm_load_pd(g+iy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+iz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r2 = _mm_load_pd(g+ix+2);
                                        r1 = _mm_load_pd(g+iy+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+iz+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r0 = _mm_add_pd (r0, r2);
                                        r2 = _mm_load_pd(g+ix+4);
                                        r1 = _mm_load_pd(g+iy+4);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+iz+4);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r0 = _mm_add_pd (r0, r2);
                                        r0 = _mm_hadd_pd(r0, r0);
                                        _mm_store_sd(gout+n, r0);
                                        gout[n] += g[ix+6] * g[iy+6] * g[iz+6];
                                }
                                break;
                        default:
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
                                break;
                } // end switch nroots
        } else { // not flag_acc
                switch (envs->nrys_roots) {
                        case 1:
                                for (n = 0; n < envs->nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] += g[ix] * g[iy] * g[iz];
                                }
                                break;
                        case 2:
                                for (n = 0; n < envs->nf-1; n+=2, idx+=6) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        jx = idx[3];
                                        jy = idx[4];
                                        jz = idx[5];
                                        r0 = _mm_load_pd(g+ix  );
                                        r1 = _mm_load_pd(g+iy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+iz  );
                                        r3 = _mm_mul_pd (r0, r1);

                                        r0 = _mm_load_pd(g+jx  );
                                        r1 = _mm_load_pd(g+jy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+jz  );
                                        r0 = _mm_mul_pd (r0, r1);

                                        r1 = _mm_loadu_pd(gout+n);
                                        r3 = _mm_hadd_pd(r3, r0);
                                        r3 = _mm_add_pd(r1, r3);
                                        _mm_storeu_pd(gout+n, r3);
                                }
                                if (n < envs->nf) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        r0 = _mm_load_pd(g+ix  );
                                        r1 = _mm_load_pd(g+iy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+iz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r0 = _mm_hadd_pd(r0, r0);
                                        _mm_store_sd(&s, r0);
                                        gout[n] += s;
                                }
                                break;
                        case 3:
                                for (n = 0; n < envs->nf-1; n+=2, idx+=6) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        jx = idx[3];
                                        jy = idx[4];
                                        jz = idx[5];
                                        r0 = _mm_load_pd(g+ix  );
                                        r1 = _mm_load_pd(g+iy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+iz  );
                                        r3 = _mm_mul_pd (r0, r1);

                                        r0 = _mm_load_pd(g+jx  );
                                        r1 = _mm_load_pd(g+jy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+jz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r3 = _mm_hadd_pd(r3, r0);

                                        r0 = _mm_loadl_pd(r0, g+ix+2);
                                        r0 = _mm_loadh_pd(r0, g+jx+2);
                                        r1 = _mm_loadl_pd(r1, g+iy+2);
                                        r1 = _mm_loadh_pd(r1, g+jy+2);
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_loadl_pd(r1, g+iz+2);
                                        r1 = _mm_loadh_pd(r1, g+jz+2);
                                        r0 = _mm_mul_pd (r0, r1);
                                        r3 = _mm_add_pd (r0, r3);

                                        r1 = _mm_loadu_pd(gout+n);
                                        r3 = _mm_add_pd(r1, r3);
                                        _mm_storeu_pd(gout+n, r3);
                                }
                                if (n < envs->nf) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        r0 = _mm_load_pd(g+ix  );
                                        r1 = _mm_load_pd(g+iy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+iz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r0 = _mm_hadd_pd(r0, r0);
                                        _mm_store_sd(&s, r0);
                                        gout[n] += s + g[ix+2] * g[iy+2] * g[iz+2];
                                }
                                break;
                        case 4:
                                for (n = 0; n < envs->nf-1; n+=2, idx+=6) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        jx = idx[3];
                                        jy = idx[4];
                                        jz = idx[5];
                                        r0 = _mm_load_pd(g+ix  );
                                        r1 = _mm_load_pd(g+iy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+iz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r2 = _mm_load_pd(g+ix+2);
                                        r1 = _mm_load_pd(g+iy+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+iz+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r3 = _mm_add_pd (r0, r2);

                                        r0 = _mm_load_pd(g+jx  );
                                        r1 = _mm_load_pd(g+jy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+jz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r2 = _mm_load_pd(g+jx+2);
                                        r1 = _mm_load_pd(g+jy+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+jz+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r0 = _mm_add_pd (r0, r2);

                                        r1 = _mm_loadu_pd(gout+n);
                                        r3 = _mm_hadd_pd(r3, r0);
                                        r3 = _mm_add_pd(r1, r3);
                                        _mm_storeu_pd(gout+n, r3);
                                }
                                if (n < envs->nf) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        r0 = _mm_load_pd(g+ix  );
                                        r1 = _mm_load_pd(g+iy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+iz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r2 = _mm_load_pd(g+ix+2);
                                        r1 = _mm_load_pd(g+iy+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+iz+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r0 = _mm_add_pd (r0, r2);
                                        r0 = _mm_hadd_pd(r0, r0);
                                        _mm_store_sd(&s, r0);
                                        gout[n] += s;
                                }
                                break;
                        case 5:
                                for (n = 0; n < envs->nf-1; n+=2, idx+=6) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        jx = idx[3];
                                        jy = idx[4];
                                        jz = idx[5];
                                        r0 = _mm_load_pd(g+ix  );
                                        r1 = _mm_load_pd(g+iy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+iz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r2 = _mm_load_pd(g+ix+2);
                                        r1 = _mm_load_pd(g+iy+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+iz+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r3 = _mm_add_pd (r0, r2);

                                        r0 = _mm_load_pd(g+jx  );
                                        r1 = _mm_load_pd(g+jy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+jz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r2 = _mm_load_pd(g+jx+2);
                                        r1 = _mm_load_pd(g+jy+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+jz+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r0 = _mm_add_pd (r0, r2);
                                        r3 = _mm_hadd_pd(r3, r0);

                                        r0 = _mm_loadl_pd(r0, g+ix+4);
                                        r0 = _mm_loadh_pd(r0, g+jx+4);
                                        r1 = _mm_loadl_pd(r1, g+iy+4);
                                        r1 = _mm_loadh_pd(r1, g+jy+4);
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_loadl_pd(r1, g+iz+4);
                                        r1 = _mm_loadh_pd(r1, g+jz+4);
                                        r0 = _mm_mul_pd (r0, r1);
                                        r3 = _mm_add_pd (r0, r3);

                                        r1 = _mm_loadu_pd(gout+n);
                                        r3 = _mm_add_pd(r1, r3);
                                        _mm_storeu_pd(gout+n, r3);
                                }
                                if (n < envs->nf) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        r0 = _mm_load_pd(g+ix  );
                                        r1 = _mm_load_pd(g+iy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+iz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r2 = _mm_load_pd(g+ix+2);
                                        r1 = _mm_load_pd(g+iy+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+iz+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r0 = _mm_add_pd (r0, r2);
                                        r0 = _mm_hadd_pd(r0, r0);
                                        _mm_store_sd(&s, r0);
                                        gout[n] += s + g[ix+4] * g[iy+4] * g[iz+4];
                                }
                                break;
                        case 6:
                                for (n = 0; n < envs->nf-1; n+=2, idx+=6) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        jx = idx[3];
                                        jy = idx[4];
                                        jz = idx[5];
                                        r0 = _mm_load_pd(g+ix  );
                                        r1 = _mm_load_pd(g+iy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+iz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r2 = _mm_load_pd(g+ix+2);
                                        r1 = _mm_load_pd(g+iy+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+iz+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r0 = _mm_add_pd (r0, r2);
                                        r2 = _mm_load_pd(g+ix+4);
                                        r1 = _mm_load_pd(g+iy+4);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+iz+4);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r3 = _mm_add_pd (r0, r2);

                                        r0 = _mm_load_pd(g+jx  );
                                        r1 = _mm_load_pd(g+jy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+jz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r2 = _mm_load_pd(g+jx+2);
                                        r1 = _mm_load_pd(g+jy+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+jz+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r0 = _mm_add_pd (r0, r2);
                                        r2 = _mm_load_pd(g+jx+4);
                                        r1 = _mm_load_pd(g+jy+4);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+jz+4);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r0 = _mm_add_pd (r0, r2);

                                        r1 = _mm_loadu_pd(gout+n);
                                        r3 = _mm_hadd_pd(r3, r0);
                                        r3 = _mm_add_pd(r1, r3);
                                        _mm_storeu_pd(gout+n, r3);
                                }
                                if (n < envs->nf) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        r0 = _mm_load_pd(g+ix  );
                                        r1 = _mm_load_pd(g+iy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+iz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r2 = _mm_load_pd(g+ix+2);
                                        r1 = _mm_load_pd(g+iy+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+iz+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r0 = _mm_add_pd (r0, r2);
                                        r2 = _mm_load_pd(g+ix+4);
                                        r1 = _mm_load_pd(g+iy+4);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+iz+4);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r0 = _mm_add_pd (r0, r2);
                                        r0 = _mm_hadd_pd(r0, r0);
                                        _mm_store_sd(&s, r0);
                                        gout[n] += s;
                                }
                                break;
                        case 7:
                                for (n = 0; n < envs->nf-1; n+=2, idx+=6) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        jx = idx[3];
                                        jy = idx[4];
                                        jz = idx[5];
                                        r0 = _mm_load_pd(g+ix  );
                                        r1 = _mm_load_pd(g+iy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+iz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r2 = _mm_load_pd(g+ix+2);
                                        r1 = _mm_load_pd(g+iy+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+iz+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r0 = _mm_add_pd (r0, r2);
                                        r2 = _mm_load_pd(g+ix+4);
                                        r1 = _mm_load_pd(g+iy+4);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+iz+4);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r3 = _mm_add_pd (r0, r2);

                                        r0 = _mm_load_pd(g+jx  );
                                        r1 = _mm_load_pd(g+jy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+jz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r2 = _mm_load_pd(g+jx+2);
                                        r1 = _mm_load_pd(g+jy+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+jz+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r0 = _mm_add_pd (r0, r2);
                                        r2 = _mm_load_pd(g+jx+4);
                                        r1 = _mm_load_pd(g+jy+4);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+jz+4);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r0 = _mm_add_pd (r0, r2);
                                        r3 = _mm_hadd_pd(r3, r0);

                                        r0 = _mm_loadl_pd(r0, g+ix+6);
                                        r0 = _mm_loadh_pd(r0, g+jx+6);
                                        r1 = _mm_loadl_pd(r1, g+iy+6);
                                        r1 = _mm_loadh_pd(r1, g+jy+6);
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_loadl_pd(r1, g+iz+6);
                                        r1 = _mm_loadh_pd(r1, g+jz+6);
                                        r0 = _mm_mul_pd (r0, r1);
                                        r3 = _mm_add_pd (r0, r3);

                                        r1 = _mm_loadu_pd(gout+n);
                                        r3 = _mm_add_pd(r1, r3);
                                        _mm_storeu_pd(gout+n, r3);
                                }
                                if (n < envs->nf) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        r0 = _mm_load_pd(g+ix  );
                                        r1 = _mm_load_pd(g+iy  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r1 = _mm_load_pd(g+iz  );
                                        r0 = _mm_mul_pd (r0, r1);
                                        r2 = _mm_load_pd(g+ix+2);
                                        r1 = _mm_load_pd(g+iy+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+iz+2);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r0 = _mm_add_pd (r0, r2);
                                        r2 = _mm_load_pd(g+ix+4);
                                        r1 = _mm_load_pd(g+iy+4);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r1 = _mm_load_pd(g+iz+4);
                                        r2 = _mm_mul_pd (r2, r1);
                                        r0 = _mm_add_pd (r0, r2);
                                        r0 = _mm_hadd_pd(r0, r0);
                                        _mm_store_sd(&s, r0);
                                        gout[n] += s + g[ix+6] * g[iy+6] * g[iz+6];
                                }
                                break;
                        default:
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
                                break;
                } // end switch nroots
        }
}


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

#include <complex.h>

void CINT2e_core(double *gout, double *g, double fac1i,
                 CINTEnvVars *envs, int empty);

void CINTgout2e      (double *gout, double *g, int *idx, CINTEnvVars *envs);
void CINTgout2e_simd1(double *gout, double *g, int *idx, CINTEnvVars *envs);

int CINT2e_loop(double *gctr, CINTEnvVars *envs, CINTOpt *opt, double *cache);

CACHE_SIZE_T CINT2e_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
                        double *cache, void (*f_c2s)());
CACHE_SIZE_T CINT2e_spinor_drv(double complex *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
                        double *cache, void (*f_e1_c2s)(), void (*f_e2_c2s)());

CACHE_SIZE_T CINT3c2e_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
                          double *cache, void (*f_c2s)(), int is_ssc);
CACHE_SIZE_T CINT3c2e_spinor_drv(double complex *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
                          double *cache, void (*f_e1_c2s)());
CACHE_SIZE_T CINT2c2e_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
                          double *cache, void (*f_c2s)());
CACHE_SIZE_T CINT2c2e_spinor_drv(double complex *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
                          double *cache, void (*f_e1_c2s)());

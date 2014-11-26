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


void CINTprim_to_ctr_0(double *gc, const int nf, const double *gp,
                       const int nprim, const int nctr, const double *coeff);
void CINTprim_to_ctr_1(double *gc, const int nf, const double *gp,
                       const int nprim, const int nctr, const double *coeff);

void CINTgout2e(double *g, double *gout, const int *idx,
                const CINTEnvVars *envs, int gout_empty);

int CINT2e_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt);

int CINT2e_drv(double *opijkl, CINTEnvVars *envs, const CINTOpt *opt,
               void (*const f_e1_c2s)(), void (*const f_e2_c2s)());

int CINT2e_cart_drv(double *opijkl, CINTEnvVars *envs, const CINTOpt *opt);
int CINT2e_spheric_drv(double *opijkl, CINTEnvVars *envs, const CINTOpt *opt);
int CINT2e_spinor_drv(double *opijkl, CINTEnvVars *envs, const CINTOpt *opt,
                      void (*const f_e1_c2s)(), void (*const f_e2_c2s)());

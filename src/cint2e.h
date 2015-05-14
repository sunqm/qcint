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

#include "config.h"

void CINT2e_core(double *gout, double *g, double fac1i,
                 CINTEnvVars *envs, FINT empty);

void CINTgout2e(double *g, double *gout, const FINT *idx,
                const CINTEnvVars *envs, FINT gout_empty);

FINT CINT2e_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt);

FINT CINT2e_drv(double *opijkl, CINTEnvVars *envs, const CINTOpt *opt,
               void (*const f_e1_c2s)(), void (*const f_e2_c2s)());

FINT CINT2e_cart_drv(double *opijkl, CINTEnvVars *envs, const CINTOpt *opt);
FINT CINT2e_spheric_drv(double *opijkl, CINTEnvVars *envs, const CINTOpt *opt);
FINT CINT2e_spinor_drv(double *opijkl, CINTEnvVars *envs, const CINTOpt *opt,
                      void (*const f_e1_c2s)(), void (*const f_e2_c2s)());

/*
 * Qcint is a general GTO integral library for computational chemistry
 * Copyright (C) 2014 Qiming Sun <osirpt.sun@gmail.com>
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

int CINT1e_loop(double *gctr, CINTEnvVars *envs, double fac);

int CINT1e_nuc_loop(double *gctr, CINTEnvVars *envs, double fac, int nuc_id);

int CINT1e_drv(double *opij, CINTEnvVars *envs, double fac,
               void (*const f_c2s)());

int CINT1e_rinv_drv(double *opij, CINTEnvVars *envs, double fac,
                    void (*const f_c2s)());

int CINT1e_nuc_drv(double *opij, CINTEnvVars *envs, double fac,
                    void (*const f_c2s)());


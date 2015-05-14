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

void CINTrys_roots(FINT nroots, double x, double *u, double *w);
static void Root123(FINT n, double x, double roots[], double weights[]);
static void Root4(double x, double roots[], double weights[]);
static void Root5(double x, double roots[], double weights[]);
static void R_droot(FINT nroots, double x, double roots[], double weights[]);
static void R_qroot(FINT nroots, double x, double roots[], double weights[]);

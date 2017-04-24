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

#include "cint_const.h"

#define bas(SLOT,I)     bas[BAS_SLOTS * (I) + (SLOT)]
#define atm(SLOT,I)     atm[ATM_SLOTS * (I) + (SLOT)]

int CINTlen_cart(const int l);
int CINTlen_spinor(const int bas_id, const int *bas);

int CINTcgtos_cart(const int bas_id, const int *bas);
int CINTcgtos_spheric(const int bas_id, const int *bas);
int CINTcgtos_spinor(const int bas_id, const int *bas);
int CINTcgto_cart(const int bas_id, const int *bas);
int CINTcgto_spheric(const int bas_id, const int *bas);
int CINTcgto_spinor(const int bas_id, const int *bas);

int CINTtot_pgto_spheric(const int *bas, const int nbas);
int CINTtot_pgto_spinor(const int *bas, const int nbas);

int CINTtot_cgto_cart(const int *bas, const int nbas);
int CINTtot_cgto_spheric(const int *bas, const int nbas);
int CINTtot_cgto_spinor(const int *bas, const int nbas);

void CINTshells_cart_offset(int ao_loc[], const int *bas, const int nbas);
void CINTshells_spheric_offset(int ao_loc[], const int *bas, const int nbas);
void CINTshells_spinor_offset(int ao_loc[], const int *bas, const int nbas);

void CINTcart_comp(int *nx, int *ny, int *nz, const int lmax);


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

#define __CVAR_CALL__   atm, *natm, bas, *nbas, env
#define __FVAR_FUNC__   const int *atm, const int *natm, \
                        const int *bas, const int *nbas, const double *env
#define C2F_(NAME)      int NAME##_(double *op, const int *shls, \
                                    __FVAR_FUNC__) \
                        { return NAME(op, shls, __CVAR_CALL__); }
// C2Fo for 2e integrals with optimizer
#define C2Fo_(NAME)     int NAME##_(double *op, const int *shls, \
                                    __FVAR_FUNC__, unsigned long *optptr) \
                        { CINTOpt *opt = (CINTOpt *)*optptr; \
                            return NAME(op, shls, __CVAR_CALL__, opt); }

typedef long CINTOptPtrAsInteger8;
#define OPTIMIZER2F_(NAME) \
    void NAME##_(CINTOptPtrAsInteger8 *optptr, __FVAR_FUNC__) { \
        CINTOpt **opt = (CINTOpt **)optptr; \
        NAME(opt, __CVAR_CALL__); }

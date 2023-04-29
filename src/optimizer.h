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

#include "cint.h"

#define NOVALUE                 ((void *)0xffffffffffffffffuL)
#define MAX_PGTO_FOR_PAIRDATA   2048

void CINTinit_2e_optimizer(CINTOpt **opt, int *atm, int natm,
                           int *bas, int nbas, double *env);
void CINTinit_optimizer(CINTOpt **opt, int *atm, int natm,
                        int *bas, int nbas, double *env);
void CINTdel_2e_optimizer(CINTOpt **opt);
void CINTdel_optimizer(CINTOpt **opt);
void CINTdel_pairdata_optimizer(CINTOpt *cintopt);
void CINTOpt_log_max_pgto_coeff(double *log_maxc, double *coeff, int nprim, int nctr);
void CINTOpt_set_log_maxc(CINTOpt *opt, int *atm, int natm,
                          int *bas, int nbas, double *env);
void CINTOpt_setij(CINTOpt *opt, int *ng,
                   int *atm, int natm, int *bas, int nbas, double *env);
void CINTOpt_non0coeff_byshell(int *sortedidx, int *non0ctr, double *ci,
                               int iprim, int ictr);
void CINTOpt_set_non0coeff(CINTOpt *opt, int *atm, int natm,
                           int *bas, int nbas, double *env);
int CINTset_pairdata(PairData *pairdata, double *ai, double *aj, double *ri, double *rj,
                     double *log_maxci, double *log_maxcj,
                     int li_ceil, int lj_ceil, int iprim, int jprim,
                     double rr_ij, double expcutoff, double *env);

// optimizer examples
void CINTno_optimizer(CINTOpt **opt, int *atm, int natm,
                      int *bas, int nbas, double *env);
void CINTall_1e_optimizer(CINTOpt **opt, int *ng,
                          int *atm, int natm, int *bas, int nbas, double *env);
void CINTall_2e_optimizer(CINTOpt **opt, int *ng,
                          int *atm, int natm, int *bas, int nbas, double *env);
void CINTall_3c2e_optimizer(CINTOpt **opt, int *ng,
                            int *atm, int natm, int *bas, int nbas, double *env);
void CINTall_2c2e_optimizer(CINTOpt **opt, int *ng,
                            int *atm, int natm, int *bas, int nbas, double *env);
void CINTall_3c1e_optimizer(CINTOpt **opt, int *ng,
                            int *atm, int natm, int *bas, int nbas, double *env);
void CINTall_1e_grids_optimizer(CINTOpt **opt, int *ng,
                                int *atm, int natm, int *bas, int nbas, double *env);

#ifdef WITH_F12
void CINTall_2e_stg_optimizer(CINTOpt **opt, int *ng,
                              int *atm, int natm, int *bas, int nbas, double *env);
#endif

#ifdef WITH_GTG
void CINTall_2e_gtg_optimizer(CINTOpt **opt, int *ng,
                              int *atm, int natm, int *bas, int nbas, double *env);
void CINTall_3c2e_gtg_optimizer(CINTOpt **opt, int *ng,
                                int *atm, int natm, int *bas, int nbas, double *env);
void CINTall_2c2e_gtg_optimizer(CINTOpt **opt, int *ng,
                                int *atm, int natm, int *bas, int nbas, double *env);
#endif

#ifndef HAVE_DEFINED_APPROX_LOG
#define HAVE_DEFINED_APPROX_LOG
#ifdef __X86__
//// little endian on x86
//typedef union {
//    double d;
//    unsigned short s[4];
//} type_IEEE754;
//// ~4 times faster than built-in log
//static inline double approx_log(double x)
//{
//        type_IEEE754 y;
//        y.d = x;
//        return ((y.s[3] >> 4) - 1023 + 1) * 0.693145751953125;
//}
#define approx_log      log
#else
#define approx_log      log
#endif
#endif

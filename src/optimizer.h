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

#if !defined HAVE_DEFINED_CINTOPT_H
#define HAVE_DEFINED_CINTOPT_H
typedef struct {
    double rij[3];
    double eij;
    int cceij;
    int _padding;
} PairData;
typedef struct {
    int **index_xyz_array; // LMAX1**4 pointers to index_xyz
    int **non0ctr;
    int **sortedidx;
    int nbas;

    int *data_ptr;
    PairData *data;
} CINTOpt;
#endif

#define NOVALUE    0xffffffff

void CINTinit_2e_optimizer(CINTOpt **opt, int *atm, int natm,
                           int *bas, int nbas, double *env);
void CINTinit_optimizer(CINTOpt **opt, int *atm, int natm,
                        int *bas, int nbas, double *env);
void CINTdel_2e_optimizer(CINTOpt **opt);
void CINTdel_optimizer(CINTOpt **opt);
void CINTOpt_setij(CINTOpt *opt, int *ng,
                   int *atm, int natm, int *bas, int nbas, double *env);
void CINTOpt_non0coeff_byshell(int *sortedidx, int *non0ctr, double *ci,
                               int iprim, int ictr);
void CINTOpt_set_non0coeff(CINTOpt *opt, int *atm, int natm,
                           int *bas, int nbas, double *env);

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

void CINTall_2e_stg_optimizer(CINTOpt **opt, int *ng,
                              int *atm, int natm, int *bas, int nbas, double *env);

int CINTset_pairdata(PairData *pdata, double *ai, double *aj, double *ri, double *rj,
                     int li_ceil, int lj_ceil, int iprim, int jprim, double rr_ij);

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

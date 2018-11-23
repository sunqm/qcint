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


#include <stdlib.h>
#include "cint_bas.h"
#include "g2e.h"
#include "optimizer.h"
#include "cint2e.h"

int int2e_stg_sph(double *out, int *dims, int *shls, int *atm, int natm,
                  int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_stg_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        envs.f_gout_simd1 = &CINTgout2e_simd1;
        return CINT2e_spheric_drv(out, dims, &envs, opt, cache);
}
void int2e_stg_optimizer(CINTOpt **opt, int *atm, int natm,
                         int *bas, int nbas, double *env)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTall_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env);
}

int int2e_yp_sph(double *out, int *dims, int *shls, int *atm, int natm,
                 int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_yp_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        envs.f_gout_simd1 = &CINTgout2e_simd1;
        return CINT2e_spheric_drv(out, dims, &envs, opt, cache);
}
void int2e_yp_optimizer(CINTOpt **opt, int *atm, int natm,
                        int *bas, int nbas, double *env)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTall_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env);
}

#define ALL_CINT(NAME) \
int c##NAME##_sph(double *out, int *shls, int *atm, int natm, \
            int *bas, int nbas, double *env, CINTOpt *opt) { \
        return NAME##_sph(out, NULL, shls, atm, natm, bas, nbas, env, opt, NULL); \
} \
void c##NAME##_sph_optimizer(CINTOpt **opt, int *atm, int natm, \
                         int *bas, int nbas, double *env) { \
        NAME##_optimizer(opt, atm, natm, bas, nbas, env); \
}

ALL_CINT(int2e_yp)
ALL_CINT(int2e_stg)



/*
 * ((NABLA i) j| F12 |k l)
 */
void CINTgout2e_int2e_ip1(double *RESTRICT gout,
        double *RESTRICT g, int *RESTRICT idx, CINTEnvVars *envs);
void CINTgout2e_int2e_ip1_simd1(double *RESTRICT gout,
        double *RESTRICT g, int *RESTRICT idx, CINTEnvVars *envs);

void int2e_yp_ip1_optimizer(CINTOpt **opt, int *atm, int natm,
                             int *bas, int nbas, double *env) {
        int ng[] = {1, 0, 0, 0, 1, 1, 1, 3};
        CINTall_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
int int2e_yp_ip1_sph(double *out, int *dims, int *shls, int *atm, int natm,
                      int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {1, 0, 0, 0, 1, 1, 1, 3};
        CINTEnvVars envs;
        CINTinit_int2e_yp_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int2e_ip1;
        envs.f_gout_simd1 = &CINTgout2e_int2e_ip1_simd1;
        return CINT2e_spheric_drv(out, dims, &envs, opt, cache);
}
ALL_CINT(int2e_yp_ip1)

void int2e_stg_ip1_optimizer(CINTOpt **opt, int *atm, int natm,
                             int *bas, int nbas, double *env) {
        int ng[] = {1, 0, 0, 0, 1, 1, 1, 3};
        CINTall_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
int int2e_stg_ip1_sph(double *out, int *dims, int *shls, int *atm, int natm,
                      int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {1, 0, 0, 0, 1, 1, 1, 3};
        CINTEnvVars envs;
        CINTinit_int2e_stg_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int2e_ip1;
        envs.f_gout_simd1 = &CINTgout2e_int2e_ip1_simd1;
        return CINT2e_spheric_drv(out, dims, &envs, opt, cache);
}
ALL_CINT(int2e_stg_ip1)

/*
 * ((NABLA NABLA i) j| F12 |k l)
 */
void CINTgout2e_int2e_ipip1(double *RESTRICT gout,
        double *RESTRICT g, int *RESTRICT idx, CINTEnvVars *envs);
void CINTgout2e_int2e_ipip1_simd1(double *RESTRICT gout,
        double *RESTRICT g, int *RESTRICT idx, CINTEnvVars *envs);

void int2e_yp_ipip1_optimizer(CINTOpt **opt, int *atm, int natm,
                              int *bas, int nbas, double *env) {
        int ng[] = {2, 0, 0, 0, 2, 1, 1, 9};
        CINTall_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
int int2e_yp_ipip1_sph(double *out, int *dims, int *shls, int *atm, int natm,
                       int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {2, 0, 0, 0, 2, 1, 1, 9};
        CINTEnvVars envs;
        CINTinit_int2e_yp_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int2e_ipip1;
        envs.f_gout_simd1 = &CINTgout2e_int2e_ipip1_simd1;
        return CINT2e_spheric_drv(out, dims, &envs, opt, cache);
}
ALL_CINT(int2e_yp_ipip1)

void int2e_stg_ipip1_optimizer(CINTOpt **opt, int *atm, int natm,
                               int *bas, int nbas, double *env) {
        int ng[] = {2, 0, 0, 0, 2, 1, 1, 9};
        CINTall_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
int int2e_stg_ipip1_sph(double *out, int *dims, int *shls, int *atm, int natm,
                        int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {2, 0, 0, 0, 2, 1, 1, 9};
        CINTEnvVars envs;
        CINTinit_int2e_stg_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int2e_ipip1;
        envs.f_gout_simd1 = &CINTgout2e_int2e_ipip1_simd1;
        return CINT2e_spheric_drv(out, dims, &envs, opt, cache);
}
ALL_CINT(int2e_stg_ipip1)

/*
 * ((NABLA i) (NABLA j)| F12 |k l)
 */
void CINTgout2e_int2e_ipvip1(double *RESTRICT gout,
        double *RESTRICT g, int *RESTRICT idx, CINTEnvVars *envs);
void CINTgout2e_int2e_ipvip1_simd1(double *RESTRICT gout,
        double *RESTRICT g, int *RESTRICT idx, CINTEnvVars *envs);

void int2e_yp_ipvip1_optimizer(CINTOpt **opt, int *atm, int natm,
                               int *bas, int nbas, double *env) {
        int ng[] = {1, 1, 0, 0, 2, 1, 1, 9};
        CINTall_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
int int2e_yp_ipvip1_sph(double *out, int *dims, int *shls, int *atm, int natm,
                        int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {1, 1, 0, 0, 2, 1, 1, 9};
        CINTEnvVars envs;
        CINTinit_int2e_yp_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int2e_ipvip1;
        envs.f_gout_simd1 = &CINTgout2e_int2e_ipvip1_simd1;
        return CINT2e_spheric_drv(out, dims, &envs, opt, cache);
}
ALL_CINT(int2e_yp_ipvip1)

void int2e_stg_ipvip1_optimizer(CINTOpt **opt, int *atm, int natm,
                                int *bas, int nbas, double *env) {
        int ng[] = {1, 1, 0, 0, 2, 1, 1, 9};
        CINTall_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
int int2e_stg_ipvip1_sph(double *out, int *dims, int *shls, int *atm, int natm,
                         int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {1, 1, 0, 0, 2, 1, 1, 9};
        CINTEnvVars envs;
        CINTinit_int2e_stg_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int2e_ipvip1;
        envs.f_gout_simd1 = &CINTgout2e_int2e_ipvip1_simd1;
        return CINT2e_spheric_drv(out, dims, &envs, opt, cache);
}
ALL_CINT(int2e_stg_ipvip1)

/*
 * ((NABLA i) j| F12 |(NABLA k) l)
 */
void CINTgout2e_int2e_ip1ip2(double *RESTRICT gout,
        double *RESTRICT g, int *RESTRICT idx, CINTEnvVars *envs);
void CINTgout2e_int2e_ip1ip2_simd1(double *RESTRICT gout,
        double *RESTRICT g, int *RESTRICT idx, CINTEnvVars *envs);

void int2e_yp_ip1ip2_optimizer(CINTOpt **opt, int *atm, int natm,
                               int *bas, int nbas, double *env) {
        int ng[] = {1, 0, 1, 0, 2, 1, 1, 9};
        CINTall_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
int int2e_yp_ip1ip2_sph(double *out, int *dims, int *shls, int *atm, int natm,
                        int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {1, 0, 1, 0, 2, 1, 1, 9};
        CINTEnvVars envs;
        CINTinit_int2e_yp_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int2e_ip1ip2;
        envs.f_gout_simd1 = &CINTgout2e_int2e_ip1ip2_simd1;
        return CINT2e_spheric_drv(out, dims, &envs, opt, cache);
}
ALL_CINT(int2e_yp_ip1ip2)

void int2e_stg_ip1ip2_optimizer(CINTOpt **opt, int *atm, int natm,
                                int *bas, int nbas, double *env) {
        int ng[] = {1, 0, 1, 0, 2, 1, 1, 9};
        CINTall_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
int int2e_stg_ip1ip2_sph(double *out, int *dims, int *shls, int *atm, int natm,
                         int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {
        int ng[] = {1, 0, 1, 0, 2, 1, 1, 9};
        CINTEnvVars envs;
        CINTinit_int2e_stg_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int2e_ip1ip2;
        envs.f_gout_simd1 = &CINTgout2e_int2e_ip1ip2_simd1;
        return CINT2e_spheric_drv(out, dims, &envs, opt, cache);
}
ALL_CINT(int2e_stg_ip1ip2)


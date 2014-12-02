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


#include <math.h>
#include "g2e.h"

static double rys_root1(double x)
{
        const double pie4 = 7.85398163397448e-01;
        double ww1;
        double f1,e,y;
        if (x < 3.e-7){
                ww1 = 1.0e+00 -x/3.0e+00;
        } else if (x < 1.) {
                f1 = ((((((((-8.36313918003957e-08*x+1.21222603512827e-06 )*x-
                            1.15662609053481e-05 )*x+9.25197374512647e-05 )*x-
                          6.40994113129432e-04 )*x+3.78787044215009e-03 )*x-
                        1.85185172458485e-02 )*x+7.14285713298222e-02 )*x-
                      1.99999999997023e-01 )*x+3.33333333333318e-01;
                ww1 = (x+x)*f1+exp(-x);
        } else if (x < 3.) {
                y = x-2.0e+00;
                f1 = ((((((((((-1.61702782425558e-10*y+1.96215250865776e-09 )*y-
                              2.14234468198419e-08 )*y+2.17216556336318e-07 )*y-
                            1.98850171329371e-06 )*y+1.62429321438911e-05 )*y-
                          1.16740298039895e-04 )*y+7.24888732052332e-04 )*y-
                        3.79490003707156e-03 )*y+1.61723488664661e-02 )*y-
                      5.29428148329736e-02 )*y+1.15702180856167e-01;
                ww1 = (x+x)*f1+exp(-x);
        } else if (x < 5.) {
                y = x-4.0e+00;
                f1 = ((((((((((-2.62453564772299e-11*y+3.24031041623823e-10 )*y-
                              3.614965656163e-09)*y+3.760256799971e-08)*y-
                            3.553558319675e-07)*y+3.022556449731e-06)*y-
                          2.290098979647e-05)*y+1.526537461148e-04)*y-
                        8.81947375894379e-04 )*y+4.33207949514611e-03 )*y-
                      1.75257821619926e-02 )*y+5.28406320615584e-02;
                ww1 = (x+x)*f1+exp(-x);
        } else if (x < 10) {
                e = exp(-x);
                ww1 = (((((( 4.6897511375022e-01/x-6.9955602298985e-01)/x +
                           5.3689283271887e-01)/x-3.2883030418398e-01)/x +
                         2.4645596956002e-01)/x-4.9984072848436e-01)/x -
                       3.1501078774085e-06)*e + sqrt(pie4/x);
        } else if (x < 15) {
                e = exp(-x);
                ww1 = (((-1.8784686463512e-01/x+2.2991849164985e-01)/x -
                        4.9893752514047e-01)/x-2.1916512131607e-05)*e
                        + sqrt(pie4/x);
        } else if (x < 33) {
                e = exp(-x);
                ww1 = (( 1.9623264149430e-01/x-4.9695241464490e-01)/x -
                       6.0156581186481e-05)*e + sqrt(pie4/x);
        } else {
                ww1 = sqrt(pie4/x);
        }
        return ww1;
}

double CINTg0_2e_ssss(const double fac, const CINTEnvVars *envs)
{
        const double aij = envs->aij;
        const double akl = envs->akl;
        double a0, a1, x;
        double rijrkl[3];
        rijrkl[0] = envs->rij[0] - envs->rkl[0];
        rijrkl[1] = envs->rij[1] - envs->rkl[1];
        rijrkl[2] = envs->rij[2] - envs->rkl[2];

        a1 = aij * akl;
        a0 = a1 / (aij + akl);
        x = a0 *(rijrkl[0] * rijrkl[0]
               + rijrkl[1] * rijrkl[1]
               + rijrkl[2] * rijrkl[2]);

        return sqrt(a0 / (a1 * a1 * a1)) * fac * rys_root1(x);
}


/*#    Copyright (C) 2012 Daniel Gamermann <gamermann@gmail.com>
#
#    This file is part of ExGUtils
#
#    ExGUtils is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    ExGUtils is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with ExGUtils.  If not, see <http://www.gnu.org/licenses/>.
#
#    
#    Please, cite us in your reasearch!
# */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "integral.h"
#include "routs.h"

/* ******************
      Integrals
  ****************** */


double const Xs[10] = {.993128599185094924786,.963971927277913791268,.912234428251325905868,
       .839116971822218823395,.74633190646015079614,.636053680726515025453,
       .510867001950827098004,.373706088715419560673,.227785851141645078080,
                     .076526521133497333755};

double const Ws[10] = {.017614007139152118312,.040601429800386941331,.062672048334109063570,
      .083276741576704748725,.101930119817240435037,.118194531961518417312,
      .131688638449176626898,.142096109318382051329,.149172986472603746788,
                   .152753387130725850698};


void gen_points(double ini, double fin, int N, double *XXX) {
/* Generates the points in the X axis */
    int ii, jj;
    double dxi, bet, alp, A, B;
    dxi = (fin - ini) / N;
    bet = (dxi)*.5;
    for (ii=0; ii<N; ii++) {
        A = ini + ii*dxi;
        B = A + dxi;
        alp = (A+B)*.5;
        for (jj=0; jj<10; jj++){
            *(XXX + jj + ii*20) = -bet*Xs[jj]+alp;
        }
        for (jj=0; jj<10; jj++){
            *(XXX + jj + 10 + ii*20) = bet*Xs[9-jj]+alp;
        }
    }
}


double integrate(double ini, double fin, int N, double *YYY){
/* Integrates the function */
    int ii, jj;
    double sum = 0., dxi;
    dxi = (fin - ini) / N;
    for (ii=0; ii<N; ii++) {
        for (jj=0; jj<10; jj++){
            sum += Ws[jj]* *(YYY+20*ii+jj);
        }
        for (jj=0; jj<10; jj++){
            sum += Ws[9-jj]* *(YYY+20*ii+jj+10);
        }
    }
    return(.5*dxi*sum);
}






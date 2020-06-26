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
#include "randoms.h"
#include <time.h>

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-14
#define RNMX (1.0-EPS)
#define PI 3.14159265358979323846


long *idum = NULL;
double *gauss_number = NULL;

void seed(void){
    idum = malloc(sizeof(long));
    *idum = -time(NULL);
}



double drand(void) {
    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    double temp;
    if (!idum) {
        seed();
    }
    if (*idum <= 0 || !iy) {
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        for (j=NTAB+7;j>=0;j--) { 
            k=(*idum)/IQ;
            *idum=IA*(*idum-k*IQ)-IR*k;
            if (*idum < 0) *idum += IM;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ; 
    *idum=IA*(*idum-k*IQ)-IR*k; 
    if (*idum < 0) *idum += IM; 
    j=iy/NDIV; 
    iy=iv[j]; 
    iv[j] = *idum; 
    if ((temp=AM*iy) > RNMX) return RNMX; 
    else return temp;
}



double drand_exp(double tau) {
    return -tau*log(1.-drand());
}




double drand_gauss(double mu, double sig) {
    double r1, r2;
    double z1, z2, r, theta;
    if (gauss_number) {
        r1 = *gauss_number;
        //free(gauss_number);
        gauss_number = NULL;
        return mu + r1*sig;
    } else {
        gauss_number = malloc(sizeof(double));
        r1 = drand();
        r2 = drand();
        r = sqrt(-2.*log(r1));
        theta = 2.*PI*r2;
        z1 = r*cos(theta);
        z2 = r*sin(theta);
        *gauss_number = z2;
        return mu + z1*sig;
    }
}
















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
#include "matrices.h"




void matmat(double *A, double *B, double *C, int N) {
        int i, j, k;
        double sum;
        
        for (i=0; i<N; i++) {
                for (j=0; j<N; j++) {
                        sum = 0;
                        for (k=0; k<N; k++) {
                                sum += *(A + i*N + k) * *(B + k*N + j);
                        }
                        *(C + i*N + j) = sum;
                }
        }
}

void matvec(double *A, double *B, double *C, int N) {
        int i, k;
        double sum;
        
        for (i=0; i<N; i++) {
                sum = 0;
                for (k=0; k<N; k++) {
                        sum += *(A + i*N + k) * *(B + k);
                }
                *(C + i) = sum;
        }
}


int matr_inv(double *A, int N) {
        /*
        **  Matrix inversion. Return 0 if matrix is singular and 1 otherwise.
        */
        int i, j, k, l, ll;
        int ip[N], iR=0, iC=0;
        int indR[N], indC[N];
        double big, dum, pinv;
        
        for (j=0; j<N; j++) {
                ip[j] = 0;
        }
        
        for (i=0; i<N; i++) {
                big = 0.;
                for (j=0; j<N; j++) {
                        if (ip[j] != 1) {
                                for (k=0; k<N; k++) {
                                        if (ip[k] == 0) {
                                                if (fabs( *(A + j*N + k)) > big) {
                                                        big = abs( *(A + N*j + k) );
                                                        iR = j;
                                                        iC = k;
                                                }
                                        }
                                        else if (ip[k] > 1) {
                                                return(0);
                                        }
                                }
                        }
                }
                ip[iC]++;
                if (iR != iC) {
                        for (l=0; l<N; l++) {
                                dum = *(A + iR*N + l);
                                *(A + iR*N + l) = *(A + iC*N + l);
                                *(A + iC*N + l) = dum;
                        }
                }
                indR[i] = iR;
                indC[i] = iC;
                if ( *(A + iC*N + iC) == 0.) {
                        return(0);
                }
                pinv = 1./ *(A + iC*N + iC);
                *(A + iC*N + iC) = 1.;
                for (l=0; l<N; l++) {
                        *(A + iC*N + l) *= pinv;
                }
                for (ll=0; ll<N; ll++) {
                        if (ll != iC) {
                                dum = *(A + ll*N + iC);
                                *(A + ll*N + iC) = 0.;
                                for (l=0; l<N; l++) {
                                        *(A + ll*N + l) -= *(A + iC*N + l)*dum;
                                }
                        }
                }
        }
        for (l=N-1; l>=0; l--) {
                if (indR[l] != indC[l]) {
                        for (k=0; k<N; k++) {
                                dum = *(A + k*N + indR[l]);
                                *(A + k*N + indR[l]) = *(A + k*N + indC[l]);
                                *(A + k*N + indC[l]) = dum;
                        }
                }
        }
        return(1);
}



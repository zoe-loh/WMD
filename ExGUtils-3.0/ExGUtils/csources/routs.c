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
#include "routs.h"
#include "randoms.h"


#define SQ2 1.414213562373095048801689
#define PI 3.1415926535897932384626433


double phig(double x, double mu, double sig) {
    double A;
    A = erf( (x-mu)/(sig*SQ2) );
    return .5*(1+A);
}

double exgCDF(double x, double mu, double sig, double tau) {
    double p, v, v2, A, B, C;
    p = (x-mu)/tau;
    v = sig/tau;
    v2 = v*v;
    A = phig(p, 0., v);
    B = phig(p, v2, v);
    C = exp(-p+.5*v2);
    return A-C*B;
}


double exg_q(double Q, double mu, double sig, double tau, double eps, double xi) {
    double x, F;
    if (xi==0.) {
        x = mu+tau;
    } else {
        x = xi;
    }
    F = exgCDF(x, mu, sig, tau);
    while(fabs(F-Q)>eps) {
        x = x - (F-Q)/exgauss(x, mu, sig, tau);
        F = exgCDF(x, mu, sig, tau);
    }
    return x;
}


double exgCDF_lamb(double z, double lamb) {
    double p, v, v2, A, B, C;
    p = z/lamb+1.;
    v2 = 1./(lamb*lamb)-1.;
    v = sqrt(v2);
    A = phig(p, 0., v);
    B = phig(p, v2, v);
    C = exp(-p+.5*v2);
    return A-C*B;
}


double exg_q_lamb(double Q, double lamb, double eps, double xi) {
    double x, F;
    x = xi;
    F = exgCDF_lamb(x, lamb);
    while(fabs(F-Q)>eps) {
        x = x - (F-Q)/exgauss_lamb(x, lamb);
        F = exgCDF_lamb(x, lamb);
    }
    return x;
}


void histogram(double *lista, int N, double ini, double fin, double Nint, int accu, int norm, double dell, double *xi, double *ni) {
    int i, j;
    double sum, dx, ele, inter, suma;
        
    dx = (fin-ini)/Nint;
    if (norm>0) {
        suma = 1./(dx*N);
    } else if (norm==0) {
        suma = 1.;
    } else {
        suma = 1./N;
    }
    for (i=0; i<Nint; i++) {
        *(ni+i) = 0.;
        *(xi+i) = ini + dx*(i+dell);
    }
    for (i=0; i<N; i++) {
        ele = *(lista+i);
        if ( (ele>=ini) & (ele<fin) ) {
            inter = ((ele-ini)/dx);
            j = (int) inter;
            *(ni+j) += 1.*suma;
        }
    }
    if (accu==-1) {
        for (i=0; i<Nint; i++) {
            sum = 0;
            for (j=i; j<Nint; j++) {
                sum += *(ni+j);
            }
            *(ni+i) = sum;
        }
    } else if (accu==1) {
        for (i=Nint-1; i>=0; i--) {
            sum = 0;
            for (j=0; j<=i; j++) {
                sum += *(ni+j);
            }
            *(ni+i) = sum;
        }
    }
}







void stats(double *nums, int N, double *xb, double *sig, double *t, int assy) {
    int i;
    double sum=0., sum2=0., sum3=0., xi;
    
    if (assy) {
        for (i=0; i<N; i++) {
            xi = *(nums+i);
            sum += xi;
            sum2 += xi*xi;
            sum3 += xi*xi*xi;
        }        
        *xb = sum*1./N;
        *sig = sqrt((sum2 - N* *xb * *xb)/(N-1));
        *t = (sum3/N - 3.* *xb * *sig * *sig - *xb * *xb * *xb)/(*sig * *sig * *sig);
    } else {    
        for (i=0; i<N; i++) {
            xi = *(nums+i);
            sum += xi;
            sum2 += xi*xi;
        }
        *xb = sum*1./N;
        *sig = sqrt((sum2 - N* *xb * *xb)/(N-1));
    }
}




void stats_his(double *xx, double *yi, int NN, double *xb, double *sig, double *t, int assy, int norm, int NNN) {
    int i;
    double sum=0., sum2=0., sum3=0., xi, nn, N=0., dx, bla;
    
    if (assy) {
        for (i=0; i<NN; i++) {
            xi = *(xx+i);
            nn = *(yi+i);
            N += nn;
            sum += xi*nn;
            sum2 += xi*xi*nn;
            sum3 += xi*xi*xi*nn;
        }
        if (norm==0) {
            *xb = sum*1./N;
            *sig = sqrt((sum2 - N* *xb * *xb)/(N-1));
            *t = (sum3/N - 3.* *xb * *sig * *sig - *xb * *xb * *xb)/(*sig * *sig * *sig);
        } else if (norm==-1) {
            *xb = sum;
            *sig = sqrt((sum2 - *xb * *xb));
            bla = sqrt(1.*NNN)/sqrt(1.*NNN-1.);  // N should be the total amount of data, not the number of intervals!!!!!
            *sig = bla* *sig;
            *t = (sum3 - 3.* *xb * *sig * *sig - *xb * *xb * *xb)/(*sig * *sig * *sig);
        } else if (norm==1) {
            dx = *(xx+1) - *xx;
            *xb = sum*dx;
            *sig = sqrt((sum2*dx - *xb * *xb));
            bla = sqrt(1.*NNN)/sqrt(1.*NNN-1.);
            *sig = bla* *sig;
            *t = (sum3*dx - 3.* *xb * *sig * *sig - *xb * *xb * *xb)/(*sig * *sig * *sig);
        }
    } else {    
        for (i=0; i<NN; i++) {
            xi = *(xx+i);
            nn = *(yi+i);
            N += nn;
            sum += xi*nn;
            sum2 += xi*xi*nn;
        }
        if (norm==0) {
            *xb = sum*1./N;
            *sig = sqrt((sum2 - N* *xb * *xb)/(N-1));
        } else if (norm==-1) {
            *xb = sum;
            *sig = sqrt(1.*(sum2 - *xb * *xb));
            bla = sqrt(1.*NNN)/sqrt(1.*NNN-1.);
            *sig = bla* *sig;
        } else if (norm==1) {
            dx = *(xx+1) - *xx;
            *xb = sum*dx;
            *sig = sqrt(1.*(sum2*dx - *xb * *xb));
            bla = sqrt(1.*NNN)/sqrt(1.*NNN-1.);
            *sig = bla* *sig;
        }
    }
}





double gaussian(double x, double mu, double sig) {
    return (1./(sig*SQ2*sqrt(PI)))*exp(-.5*((x-mu)/sig)*((x-mu)/sig));
}



double exgauss(double x, double mu, double sig, double tau) {
    double arg1, arg2, bla1, bla2, sig2;
    sig2 = sig*sig;
    arg1 = 2.*mu+sig2/tau-2.*x;
    arg2 = mu+sig2/tau-x;
    bla1 = (0.5/tau)*exp((0.5/tau)*arg1);
    bla2 = 1.-erf(arg2/(SQ2*sig));
    //printf("%f  %f  %f  %f\n", arg1, arg2, bla1, bla2);
    if (isinf(bla1)) {
        return 0.;
    } else {
        return bla1*bla2;
    }
}


double exgauss_lamb(double z, double lamb) {
    double arg1, arg2, bla1, bla2, lamb2;
    lamb2 = lamb*lamb;
    arg1 = -2.*z*lamb-3*lamb2+1; 
    arg1 = arg1/(2*lamb2);
    arg2 = -z + 1./lamb - 2.*lamb;
    arg2 = arg2/sqrt(1.-lamb2);
    bla1 = exp(arg1);
    bla2 = 1.-erf(arg2/SQ2);
    if (isinf(bla1)) {
        return 0.;
    } else {
        return bla1*bla2/(2.*lamb);
    }
}


void pars_to_stats(double mu, double sig, double tau, double *M, double *S, double *lamb) {
    *M = mu+tau;
    *S = sqrt(sig*sig+tau*tau);
    *lamb = tau/ *S;
}


void stats_to_pars(double M, double S, double lamb, double *mu, double *sig, double *tau) {
    *mu = M - S*lamb;
    *sig = S*sqrt(1-lamb*lamb);
    *tau = S*lamb;
}



void exgLKHD(double *xi, int N, double mu, double sig, double tau, double *lnL, double *grad) {
    double A1, A2, eA3, sig2, sqpi, dfdm, dfds, dfdt, fi, x;
    double sumL=0., sumdm=0., sumds=0., sumdt=0.;
    int i;
    sqpi = pow(PI, 0.5);
    sig2 = sig*sig;
    for (i=0; i<N; i++) {
        x = *(xi + i);
        A1 = (2.*mu+sig2/tau-2.*x)*0.5/tau;
        A2 = (mu+sig2/tau-x)/(SQ2*sig);
        eA3 = exp(A1-A2*A2);
        fi = (0.5/tau)*exp(A1)*(1.-erf(A2)); /// Careful if fi = 0 !!!!!!!!
        //printf(" %f, %f, %f, %f, %f, %d - ", x, sumL, fi, A2, erf(A2), i);
        dfdm = fi/tau - (1./(SQ2*sqpi*sig*tau))*eA3;
        dfds = sig*fi/(tau*tau) + (A2/(sqpi*sig*tau) - SQ2/(sqpi*tau*tau))*eA3;
        dfdt = (-(A1+1)/tau - sig2/(2*tau*tau*tau))*fi + (sig/(SQ2*sqpi*tau*tau*tau))*eA3;
        if (fi!=0.0) {
        sumL += log(fi);
        sumdm += dfdm/fi;
        sumds += dfds/fi;
        sumdt += dfdt/fi;}
    }
    //printf("---------%f---------\n", sumL);
    *lnL = sumL;
    *grad = sumdm;
    *(grad+1) = sumds;
    *(grad+2) = sumdt;
}


void maxLKHD(double *xi, int N, double *mu, double *sig, double *tau, double lamb, double eps) {
    double mup, sigp, taup;
    double lkhd, lnL;
    double grad[3], ngrad[3], normgrad;
    
    exgLKHD(xi, N, *mu, *sig, *tau, &lnL, &grad[0]);
    normgrad = sqrt(grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2]);
    //printf("\n____________\n%6.18f\n", normgrad);
    while (fabs(lamb*normgrad)>eps) {
        mup = *mu + lamb* grad[0];
        sigp = *sig + lamb* grad[1];
        taup = *tau + lamb* grad[2];
        exgLKHD(xi, N, mup, sigp, taup, &lkhd, &ngrad[0]);
        //printf("lkhd: %f lnl: %f mu: %f sig: %f tau: %f g1: %f  g2: %f  g3: %f lamb: %f\n", lkhd, lnL, *mu, *sig, *tau, grad[0], grad[1], grad[2], lamb);
        if (lkhd>lnL) {
            //printf("up!\n");
            *mu = mup;
            *sig = sigp;
            *tau = taup;
            lnL = lkhd;
            grad[0] = ngrad[0];
            grad[1] = ngrad[1];
            grad[2] = ngrad[2];
            normgrad = sqrt(grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2]);
            lamb = lamb * 1.1;
            //printf("%6.18f\n", normgrad);
        } else {
            //printf("down!\n");
            lamb = lamb*0.9;
        }
    }
    //printf("Saiu\n");
}





void exgSQR(double *xi, double *yi, int N, double mu, double sig, double tau, double *lnL, double *grad) {
    double A1, A2, eA3, sig2, sqpi, dfdm, dfds, dfdt, fi, x, y;
    double sumL=0., sumdm=0., sumds=0., sumdt=0.;
    int i;
    sqpi = pow(PI, 0.5);
    sig2 = sig*sig;
    for (i=0; i<N; i++) {
        x = *(xi + i);
        y = *(yi + i);
        A1 = (2.*mu+sig2/tau-2.*x)*0.5/tau;
        A2 = (mu+sig2/tau-x)/(SQ2*sig);
        eA3 = exp(A1-A2*A2);
        fi = (0.5/tau)*exp(A1)*(1.-erf(A2));
        dfdm = fi/tau - (1./(SQ2*sqpi*sig*tau))*eA3;
        dfds = sig*fi/(tau*tau) + (A2/(sqpi*sig*tau) - SQ2/(sqpi*tau*tau))*eA3;
        dfdt = (-(A1+1)/tau - sig2/(2*tau*tau*tau))*fi + (sig/(SQ2*sqpi*tau*tau*tau))*eA3;
        sumL += (fi - y) * (fi - y);
        sumdm += dfdm*2.*(fi - y);
        sumds += dfds*2.*(fi - y);
        sumdt += dfdt*2.*(fi - y);
        //printf(" %f, %f, %f, %f, %f, %d - ", x, sumL, fi, A2, erf(A2), i);
    }
    *lnL = sumL;
    *grad = sumdm;
    *(grad+1) = sumds;
    *(grad+2) = sumdt;
}



void minSQR(double *xi, double *yi, int N, double *mu, double *sig, double *tau, double lamb, double eps) {
    double mup, sigp, taup;
    double lkhd, lnL;
    double grad[3], ngrad[3], normgrad;
    
    exgSQR(xi, yi, N, *mu, *sig, *tau, &lnL, &grad[0]);
    normgrad = sqrt(grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2]);
    while (fabs(lamb*normgrad)>eps) {
        mup = *mu - lamb* grad[0];
        sigp = *sig - lamb* grad[1];
        taup = *tau - lamb* grad[2];
        exgSQR(xi, yi, N, mup, sigp, taup, &lkhd, &ngrad[0]);
        //printf("lkhd: %E lnl: %E mu: %E sig: %E tau: %f g1: %E  g2: %E  g3: %E lamb: %E\n", lkhd, lnL, *mu, *sig, *tau, grad[0], grad[1], grad[2], lamb);
        if (lkhd<lnL) {
            //printf("up!\n");
            *mu = mup;
            *sig = sigp;
            *tau = taup;
            lnL = lkhd;
            grad[0] = ngrad[0];
            grad[1] = ngrad[1];
            grad[2] = ngrad[2];
            normgrad = sqrt(grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2]);
            lamb = lamb * 1.1;
        } else {
            //printf("down!\n");
            lamb = lamb*0.9;
        }
    }
    //printf("Saiu\n");
}














double poly(double *a, double x, int N) {
    double f;
    int i;
    f = *a + *(a+1) * x;
    for (i=2; i<N; i++) {
        f += *(a+i) * pow( x, i );
    }
    return f;
}







double correlation(double *xs1, double *xs2, int N) {
    double xb1, s1, xb2, s2;
    double bla=0.;
    int i;
    stats(xs1, N, &xb1, &s1, NULL, 0);
    stats(xs2, N, &xb2, &s2, NULL, 0);
    for (i=0; i<N; i++) {
        bla += (*(xs1 + i)-xb1) * (*(xs2 + i)-xb2);
    }
    return bla/(s1*s2*(N-1.));
}





double minsquare(double *xs, double *ys, double *sigs, int N, int deg, double *as) {
    int i, j;
    double *csi, *phi, *lamb;
    double sum1, sum2, bla;
    double chi2=0;

    csi = malloc( (2*deg+1) * sizeof(double));
    phi = malloc( (deg+1) * sizeof(double));
    lamb = malloc( (deg+1)*(deg+1) * sizeof(double));

    if (sigs) {
        for (i=0; i<(deg+1); i++) {
            sum1=0;
            sum2=0;
            for (j=0; j<N; j++) {
                bla = pow( *(xs+j), i)/ *(sigs+j);
                sum1 += bla;
                sum2 += bla* *(ys+j);
            }
            *(csi+i) = sum1;
            *(phi+i) = sum2;
        }
        for (i=(deg+1); i<(2*deg+1); i++) {
            sum1=0;
            for (j=0; j<N; j++) {
                sum1 += pow( *(xs+j), i)/ *(sigs+j);
            }
            *(csi+i) = sum1;
        }
        for (i=0; i<(deg+1); i++) {
            for (j=0; j<(deg+1); j++) {
                *(lamb + i*(deg+1) + j) = *(csi + i + j);
            }
        }
        matr_inv(lamb, deg+1);
        matvec(lamb, phi, as, deg+1);
        for (i=0; i<N; i++) {
            chi2 += pow( *(ys+i) - poly(as, *(xs+i), deg+1), 2 )/ *(sigs+i) ;
        }
    } else {
        for (i=0; i<(deg+1); i++) {
            sum1=0;
            sum2=0;
            for (j=0; j<N; j++) {
                bla = pow( *(xs+j), i);
                sum1 += bla;
                sum2 += bla* *(ys+j);
            }
            *(csi+i) = sum1;
            *(phi+i) = sum2;
        }
        for (i=(deg+1); i<(2*deg+1); i++) {
            sum1=0;
            for (j=0; j<N; j++) {
                sum1 += pow( *(xs+j), i);
            }
            *(csi+i) = sum1;
        }
        for (i=0; i<(deg+1); i++) {
            for (j=0; j<(deg+1); j++) {
                *(lamb + i*(deg+1) + j) = *(csi + i + j);
            }
        }
        matr_inv(lamb, deg+1);
        matvec(lamb, phi, as, deg+1);
        for (i=0; i<N; i++) {
            chi2 += pow( *(ys+i) - poly(as, *(xs+i), deg+1), 2 ) ;
        }
    }
    free(csi);
    free(phi);
    free(lamb);
    return chi2/(N-deg-1);
}





double min_list(double *lista, int N) {
    int i;
    double min=*lista;
    for (i=1; i<N; i++) {
        if (*(lista+i)<min) {
            min = *(lista+i);
        }
    }
    return min;
}


double max_list(double *lista, int N) {
    int i;
    double max=*lista;
    for (i=1; i<N; i++) {
        if (*(lista+i)>max) {
            max = *(lista+i);
        }
    }
    return max;
}










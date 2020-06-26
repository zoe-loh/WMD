#    Copyright (C) 2012 Daniel Gamermann <gamermann@gmail.com>
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
#


from numpy import exp
from numpy import matrix
from numpy import array
from scipy.special import erf
from scipy.stats.distributions import f as fsnede
from math import log, pi
from random import random as drand
from random import gauss as gauss_rvs
from random import expovariate as exp_rvs



def exg_rvs(mu, sig, tau):
    """Generates a random number with ex-gaussian distribution\n\
    usage: exg_rvs(mu, sig, tau)"""
    a = gauss_rvs(mu, sig)
    b = exp_rvs(tau)
    return a+b



SQ2 = 2.**.5

def exgLKHD(nums, mu, sig, tau):
    """Evaluates the Log likelyhood and its gradient for the points in the list given the ex gaussian parameters.\n\
    Usage: [logLKHD, grad] = exgLKHD(lista, mu, sig, tau)"""
    sumL = 0.
    sumdm = 0.
    sumds = 0.
    sumdt = 0.
    sig2 = sig**2
    sqpi = pi**.5
    for x in nums:
        A1 = (2.*mu+sig2/tau-2.*x)*0.5/tau
        A2 = (mu+sig2/tau-x)/(SQ2*sig)
        eA3 = exp(A1-A2*A2)
        fi = (0.5/tau)*exp(A1)*(1.-erf(A2)) # Careful if fi = 0 !!!!!!!!
        dfdm = fi/tau - (1./(SQ2*sqpi*sig*tau))*eA3
        dfds = sig*fi/(tau*tau) + (A2/(sqpi*sig*tau) - SQ2/(sqpi*tau*tau))*eA3
        dfdt = (-(A1+1)/tau - sig2/(2*tau*tau*tau))*fi + (sig/(SQ2*sqpi*tau*tau*tau))*eA3
        if fi>0.0:
            sumL += log(fi)
            sumdm += dfdm/fi
            sumds += dfds/fi
            sumdt += dfdt/fi
    return sumL, (sumdm, sumds, sumdt)


### This algorithm is very unstable for some small sets of numbers
def maxLKHD(nums, mu=0., sig=0., tau=0., lamb=1., eps=1.e-10):
    """Finds the parameters that maximize the likelyhood for the values in a list.\n\
    Usage: [mu, sig, tau] = maxLKHD(lista, mu=0., sig=0., tau=0., lamb=1., eps=1.e-10)\n\
    The first three optional arguments are the initial point for searching (actually if the three \n\
    are zero, it evaluates mu, sig and tau with respect to the distribution average, standard deviation \n\
    and assymmetry), lamb is the the starting update in the parameters with respect to the gradient module\n\
    and eps is the precision (the interactions stop when the absolute value of the gradient \n\
    times lamb is smaller than eps)."""
    if not mu and not sig and not tau:
        M, S, t = stats(nums, 1)
        lambi = (.5*t)**(1./3)
        if lambi<0.:
            lambi = .2
        elif lambi>1.:
            lambi = .99
        mu, sig, tau = stats_to_pars(M, S, lambi)
    lkhd, grad = exgLKHD(nums, mu, sig, tau)
    norm = (grad[0]**2+grad[1]**2+grad[2]**2)**.5
    if not lamb:
        lamb = 5./norm
    elif lamb >0.:
        lamb = lamb/norm
    else:
        lamb = -lamb
    nmu = mu+lamb*grad[0]/norm
    nsig = sig+lamb*grad[1]/norm
    ntau = tau+lamb*grad[2]/norm
    nlkhd, ngrad = exgLKHD(nums, mu, sig, tau)
    while abs(lamb*norm) > eps:
        #print lamb, norm
        if nlkhd>lkhd: # accepted
            grad = ngrad
            mu = nmu
            sig = nsig
            tau = ntau
            lkhd = nlkhd
            lamb *= 1.1
        else:
            lamb *= .9
        norm = (grad[0]**2+grad[1]**2+grad[2]**2)**.5
        nmu = mu+lamb*grad[0]/norm
        nsig = sig+lamb*grad[1]/norm
        ntau = tau+lamb*grad[2]/norm
        nlkhd, ngrad = exgLKHD(nums, nmu, nsig, ntau)
    return mu, sig, tau


def exgSQR(XX, YY, mu, sig, tau):
    """Sum of the squares (yi-f(xi))**2 for the set of points and its gradient.\n\
    Usage: [chi2, grad] = exgSQR(listax, listay, mu, sig, tau)"""
    sumL = 0.
    sumdm = 0.
    sumds = 0.
    sumdt = 0.
    sig2 = sig**2
    sqpi = pi**.5
    for ii, x in enumerate(XX):
        y = YY[ii]
        A1 = (2.*mu+sig2/tau-2.*x)*0.5/tau
        A2 = (mu+sig2/tau-x)/(SQ2*sig)
        #print A1, A2,
        eA3 = exp(A1-A2*A2)
        fi = (0.5/tau)*exp(A1)*(1.-erf(A2))
        dfdm = fi/tau - (1./(SQ2*sqpi*sig*tau))*eA3
        dfds = sig*fi/(tau*tau) + (A2/(sqpi*sig*tau) - SQ2/(sqpi*tau*tau))*eA3
        dfdt = (-(A1+1)/tau - sig2/(2*tau*tau*tau))*fi + (sig/(SQ2*sqpi*tau*tau*tau))*eA3
        sumL += (fi - y) * (fi - y)
        sumdm += dfdm*2.*(fi - y)
        sumds += dfds*2.*(fi - y)
        sumdt += dfdt*2.*(fi - y)
    return sumL, (sumdm, sumds, sumdt)


def minSQR(XX, YY, mu=0., sig=0., tau=0., lamb=0., eps=1.e-10, NN=0):
    """Finds the parameters that minimize the sum of the squares. \n\
    Usage: [mu, sig, tau] = minSQR(listax, listay, mu=0., sig=0., tau=0., lamb=1., eps=1.e-10, NN=0)\n\
    The first three optional arguments are the initial point for searching (actually if the three \n\
    are zero, it evaluates mu, sig and tau with respect to the distribution average, standard deviation \n\
    and assymmetry), lamb is the the starting update in the parameters with respect to the gradient module\n\
    and eps is the precision (the interactions stop when the absolute value of the gradient \n\
    times lamb is smaller than eps).\n\
    Note that listax and listay should represent a normalized histogram of the data!"""
    N = len(XX)
    if not mu and not sig and not tau:
        if not NN:
            NN = int(N*N/4. + N/4.)
        M, S, t = stats_his(XX, YY, 1, 1, NN)
        lambi = (.5*t)**(1./3)
        if lambi<0.:
            lambi = .2
        elif lambi>1.:
            lambi = .99
        mu, sig, tau = stats_to_pars(M, S, lambi)
    suma = sum(YY)*(XX[1]-XX[0])
    YYY = [ele/suma for ele in YY]
    lkhd, grad = exgSQR(XX, YYY, mu, sig, tau)
    norm = (grad[0]**2+grad[1]**2+grad[2]**2)**.5
    if not lamb:
        lamb = 5./norm
    elif lamb >0.:
        lamb = lamb/norm
    else:
        lamb = -lamb
    nmu = mu-lamb*grad[0]
    nsig = sig-lamb*grad[1]
    ntau = tau-lamb*grad[2]
    nlkhd, ngrad = exgSQR(XX, YYY, mu, sig, tau)
    while abs(lamb*norm) > eps:
        #print lamb, norm
        if nlkhd<lkhd: # accepted
            grad = ngrad
            mu = nmu
            sig = nsig
            tau = ntau
            lkhd = nlkhd
            lamb *= 1.1
        else:
            lamb *= .9
        norm = (grad[0]**2+grad[1]**2+grad[2]**2)**.5
        nmu = mu-lamb*grad[0]
        nsig = sig-lamb*grad[1]
        ntau = tau-lamb*grad[2]
        nlkhd, ngrad = exgSQR(XX, YYY, nmu, nsig, ntau)
    return mu, sig, tau

def correlation(x1, x2):
    """Evaluates the linear correlation coefficient between the numbers in two lists.\n\
    Usage: corr = correlation(lista1, lista2)"""
    N = len(x1)
    [xb1, s1] = stats(x1)
    [xb2, s2] = stats(x2)
    bla = sum([(x1[ii]-xb1)*(x2[ii]-xb2) for ii in xrange(N)])
    return bla/(s1*s2*(N-1))






def gauss_pdf(x, mu, sig):
    """Evaluates the gaussian density function\n\
    usage: gauss_pdf(x, mu, sig)\n\
    returns the value of the gaussian density function at point x."""
    sqpi = pi**.5
    num = (1./(sig*SQ2*sqpi))*exp(-.5*((x-mu)/sig)**2)
    return num


def gauss_cdf(x, mu, sig):
    """Evaluates the gaussian cumulative distribution\n\
    usage: gauss_cdf(x, mu, sig)\n\
    returns the values of the gaussian cumulative distribution at point x."""
    A = erf( (x-mu)/(sig*(2.**.5)) )
    return .5*(1+A)

def exg_cdf(x, mu, sig, tau):
    """Evaluates the ex-gaussian cumulative distribution\n\
    usage: exg_cdf(x, mu, sig, tau)\n\
    returns the values of the ex-gaussian cumulative distribution at point x."""
    p = (x-mu)/tau
    v = sig/tau
    v2 = v**2
    A = gauss_cdf(p, 0., v)
    B = gauss_cdf(p, v2, v)
    C = exp(-p+.5*v2)
    return A-C*B

def exg_ppf(Q, mu, sig, tau, eps=1.e-14, xi=0.):
    """Evaluates the position of the quartile Q for the exgaussian distribution\n\
    usage: exg_ppf(Q, mu, sig, tau, eps=1.e-14, xi=0.)\n\
    returns the position where the ex-gaussian distribution accumulates an amount Q.
    xi is the point where newton's methos starts the search (if 0. xi=M)."""
    if not xi:
        x = mu+tau
    else:
        x = xi
    F = exg_cdf(x, mu, sig, tau)
    while abs(F-Q)>eps:
        x -= (F-Q)/exg_pdf(x, mu, sig, tau)
        F = exg_cdf(x, mu, sig, tau)
    return x

def exg_lamb_cdf(z, lamb):
    """Evaluates the ex-gaussian cumulative distribution\n\
    usage: exg_lamb_cdf(x, lamb)\n\
    returns the values of the ex-gaussian cumulative distribution at point x."""
    p = z/lamb+1.
    v2 = 1./(lamb**2)-1.
    v = v2**.5
    A = gauss_cdf(p, 0., v)
    B = gauss_cdf(p, v2, v)
    C = exp(-p+.5*v2)
    return A-C*B

def exg_lamb_ppf(Q, lamb, eps=1.e-14):
    """Evaluates the position of the quartile Q for the exgaussian distribution\n\
    usage: exg_lamb_ppf(Q, lamb, eps=1.e-14)\n\
    returns the position where the ex-gaussian distribution accumulates an amount Q."""
    x = 0.
    F = exg_lamb_cdf(x, lamb)
    while abs(F-Q)>eps:
        x -= (F-Q)/exg_lamb_pdf(x, lamb)
        F = exg_lamb_cdf(x, lamb)
    return x








def stats(lista, assymetry=False):
    """Evaluates the statisticals for a distribution\n\
    usage: [M, S] = stats(lista)      or\n\
    [M, S, t] = stats(lista, assymetry=True)  \n\
    returns M, S and (optionally) t, the average, standard deviation and skewness of the numbers in lista."""
    N = float(len(lista))
    xmed = sum(lista)/N
    if N == 1:
        return [xmed, -1]
    s2 = [(xmed-ele)**2 for ele in lista]
    s = (sum(s2)/(N-1.))**.5
    if assymetry:
        s3 = sum([ele**3 for ele in lista])
        t = (s3/N - 3*xmed*s*s - xmed**3)/(s**3)
        return [xmed, s, t]
    return [xmed, s]


def stats_his(xi, ni, assymetry=False, norm=0, N=0):
    """Evaluates the statisticals for an histogram\n\
    usage: [M, S] = stats_his(xi, yi, assymetry=True, norm=0, N=0)      or\n\
    [M, S, t] = stats_his(xi, yi, 1)  \n\
    returns M, S and (optionally) t, the average, standard deviation and skewness of the numbers in lista."""
    nc = ni[:]
    xn = [nc[ii]*xi[ii] for ii in xrange(len(xi))]
    xn2 = [nc[ii]*(xi[ii]**2) for ii in xrange(len(xi))]
    if assymetry: xn3 = [nc[ii]*(xi[ii]**3) for ii in xrange(len(xi))]
    Nd = float(sum(ni))
    if norm==0:
        xmed = sum(xn)/Nd
        s = ((sum(xn2) - Nd*xmed*xmed)/(Nd-1.))**.5
        if assymetry: t = (sum(xn3)/Nd - 3.*xmed*s*s -xmed**3)/(s**3)
    elif norm==-1:
        xmed = sum(xn)
        s = ((sum(xn2) - xmed*xmed))**.5
        s *= (N*1./(N-1))**.5
        if assymetry: t = (sum(xn3) - 3.*xmed*s*s -xmed**3)/(s**3)
    else:
        dx = xi[1]-xi[0]
        xmed = sum(xn)*dx
        s = ((sum(xn2)*dx - xmed*xmed))**.5
        s *= (N*1./(N-1))**.5
        if assymetry: t = (sum(xn3)*dx - 3.*xmed*s*s -xmed**3)/(s**3)
    if assymetry: return [xmed, s, t]
    return [xmed, s]    



def histogram(lista, ini=None, fin=None, Nint=None, dell=.5, accu=0, norm=0., eps=1.e-10):
    """Produces an histogram for a list of numbers\n\
    usage: [x, y] = histogram(lista, ini=ini, fin=fin, Nint=nints, dell=dell, accu=accu, norm=norm)\n\
    -lista is a list of numbers,\n\
    -ini is the initial point of the first interval (default is the smallest number in lista),\n\
    -fin is the end point of the last interval (default is the biggest number in lista),\n\
    -nints is the number of intervals (default is 2*sqrt(#number of points in lista)),\n\
    -dell is the gap from the begining of the interval where the coordinate of the midle point should be (default is 0.5)\n\
    -accu should be 0, 1 or -1 for distribution, right or left cummulative distributions, respectively (default is 0),\n\
    -norm should be 0, 1 or -1 for non-normalized, normalized or frequency distribution, respectively (default is 0).\n\
    The function returns two lists, x and y with the x and y coordinates of the points in the histogram."""
    if ini==None:
        ini=min(lista)
    if fin==None:
        fin=max(lista)
    if Nint==None:
        Nint = int(2*len(lista)**.5)
    fin += eps
    dx = (fin-ini)/Nint
    N = len(lista)
    if norm > 0:
        suma = 1./(dx*N)
    elif norm == 0:
        suma = 1.
    else:
        suma = 1./N
    fin = float(fin)
    ini = float(ini)
    anch = (fin-ini)*1.0/Nint
    hist = [0. for ele in xrange(Nint)]
    for ele in lista:
        if ele >= ini and ele < fin:
            Int = int((ele-ini)/anch)
            hist[Int] += suma
    dx=(fin-ini)/Nint
    if accu==-1:
        hist = [sum(hist[ii:]) for ii in xrange(Nint)]
    if accu==1:
        hist = [sum(hist[:ii+1]) for ii in xrange(Nint)]
    xxx = [ini+dx*(ii+dell) for ii in xrange(Nint)]
    return [xxx, hist]



def exg_pdf(x,mu,sig,tau):
    """Evaluates the ex-gaussian density function\n\
    usage: exg_pdf(x, mu, sig, tau)\n\
    returns the values of the ex-gaussian density function at point x."""
    arg1 = 2.*mu+(sig**2)/tau-2.*x
    arg2 = mu+(sig**2)/tau-x
    sq2 = 2.**.5
    bla1 = (0.5/tau)*exp((0.5/tau)*arg1)
    bla2 = 1.-erf(arg2/(sq2*sig))
    return bla1*bla2

def exg_lamb_pdf(z,lamb):
    """Evaluates the ex-gaussian density function in terms of the assymetry lambda\n\
    usage: exg_lamb_pdf(x, lamb)\n\
    returns the values of the ex-gaussian density function parametrized in terms of its assymetry at point x."""
    arg1 = -2.*z*lamb-3*lamb**2+1
    arg1 /= (2*lamb**2)
    arg2 = -z + 1./lamb - 2.*lamb
    arg2 /= (1.-lamb**2)**.5
    sq2 = 2.**.5
    bla1 = exp(arg1)
    bla2 = 1.-erf(arg2/sq2)
    return bla1*bla2/(2.*lamb)

def pars_to_stats(mu, sig, tau):
    """converts parameters to statisticals\n\
    usage: [M, S, t] = pars_to_stats(mu, sig, tau)"""
    M = mu+tau
    S = (sig**2+tau**2)**.5
    lamb = tau/S
    t = 2.*(lamb**3)
    return [M, S, t]

def stats_to_pars(M, S, t):
    """converts statisticals to parameters\n\
    usage: [mu, sig, tau] = stats_to_pars(M, S, t)"""
    if t<0. or t>2.:
        print "Warning: t must be in the interval (0; 2)!"
        return None, None, None
    lamb = (0.5*t)**(1./3)
    mu = M - S*lamb
    sig = S*(1-lamb**2)**.5
    tau = S*lamb
    return [mu, sig, tau]





def minsquare(X, Y, err=[], deg=1):
    """Adjust a set of points to fit a pollynomial of any degree by the minimum square method.\n\
    Usage: [coefs_list, chi2] = minsquare(listax, listay, err=[], deg=1)\n\
    listax and listay are lists with the x and y values of the points, respectivelly. If the \n\
    keyword argument err is given, it should be a list containing the uncertainties (error bars) in \n\
    the y points. the optional argument deg is the degree of the pollynomial to be fitted (1 \n\
    means a line). The routine returns coefs_list which is a list with the coefficients of the \n\
    fitted pollynomial:\n\
    y = coefs_list[0] + coefs_list[1]* x + coefs_list[2]* x**2 + ... \n\
    chi2 is the value of the sum (y_i - y[x_i])**2 / err_i, in other words, it's the minimized sum."""
    if not err:
        err = [1. for ii in xrange(len(X))]
    np=len(X)
    csi=[sum([x**i for x in X]) for i in range(2*deg+1)]
    phi=matrix([sum([Y[j]*X[j]**i for j in range(np)]) for i in range(deg+1)])
    lamb=matrix([[csi[i+j] for j in range(deg+1)]for i in range(deg+1)])
    alpha=lamb.I*phi.T
    coefs = alpha.T.tolist()[0]
    chi2 = sum([(ele - sum([coef*(X[ii]**jj) for jj, coef in enumerate(coefs)]))**2 / err[ii] for ii, ele in enumerate(Y)])
    return coefs, chi2




def integral(func, ini, fin, Nints=20):
    """ Calculates the integral of func between ini and fin
    with 20 points gaussian method dividing the interval [ini; fin] 
    in Nint intervals."""
    Y=[.993128599185094924786,.963971927277913791268,.912234428251325905868,
       .839116971822218823395,.74633190646015079614,.636053680726515025453,
       .510867001950827098004,.373706088715419560673,.227785851141645078080,
                     .076526521133497333755]
    W=[.017614007139152118312,.040601429800386941331,.062672048334109063570,
      .083276741576704748725,.101930119817240435037,.118194531961518417312,
      .131688638449176626898,.142096109318382051329,.149172986472603746788,
                   .152753387130725850698]
    stepint = (fin-ini)/float(Nints)
    suma = 0.
    for ii in xrange(Nints):
        A = ini + ii*stepint
        B = ini + (ii+1)*stepint
        bma = (B-A)*.5
        apb = (A+B)*.5
        tosum = [W[jj]*(func(bma*Y[jj]+apb)+func(-bma*Y[jj]+apb)) for jj in xrange(10)]
        suma += bma*sum(tosum)
    return suma



def exp_rvs(tau):
    """ Generates random number with exponential distribution."""
    nrand = drand()
    return -tau*log(1.-nrand)






def zero(func, x0, eps=1.e-9, delt=.01):
    """
       Finds the zero of function func starting at x0 with precision eps.
       delt is the dx for calculation the derivative (func(x+delt)-func(x))/delt)
       Uses Newton's method.
    """
    diff = func(x0)
    xx = x0
    while abs(diff) > eps:
        der = (func(xx+delt)-func(xx))/delt
        ddd = -diff/der
        xx += ddd
        diff = func(xx)
    return xx



def ANOVA(tab):
    """
    ANOVA test for table tab (tab should be a list of lists).
     Values returned are in order:
       Fs: Value for the variable F (F of snedecor).
       glentre: degrees of freedom in between.
       gldentro: degrees of fredom inside.
       1-fsnede: left tail (p-value for the test).
    """
    r = len(tab)
    ni = [len(ele) for ele in tab]
    xbi = [stats(ele)[0] for ele in tab]
    N = sum(ni)
    XB = sum([ni[ii]*xbi[ii] for ii in xrange(r)])/N
    ssi = [sum([(ele - xbi[ii])**2 for ele in ele2]) for ii, ele2 in enumerate(tab)]
    SSdentro = sum(ssi)
    gldentro = N-r
    MSdentro = SSdentro/gldentro
    SSentre = sum([ni[ii]*(ele-XB)**2 for ii, ele in enumerate(xbi)])
    glentre = r-1
    MSentre = SSentre/glentre
    Fs = MSentre/MSdentro
    return [Fs, glentre, gldentro, 1.-fsnede.cdf(Fs,glentre,gldentro)]




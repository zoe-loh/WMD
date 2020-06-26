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

#include <Python.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <structmember.h>
#include "routs.h"
#include "matrices.h"
#include "randoms.h"
#include "integral.h"




// *******************************************
/*           Doc strings                    */
// *******************************************
//
// Module 
static char utils_docstring[] = "Uts functions: functions programmed in C (python C API).";
// Functions
static char Function_drand_docstring[] = "Generates a random number between 0 and 1\n\
                                          usage: drand()";
static char Function_histogram_docstring[] = "Produces an histogram for a list of numbers\n\
                                              usage: [x, y] = histogram(lista, ini=ini, fin=fin, Nint=nints, dell=dell, accu=accu, norm=norm)\n\
                                               -lista is a list of numbers,\n\
                                               -ini is the initial point of the first interval (default is the smallest number in lista),\n\
                                               -fin is the end point of the last interval (default is the biggest number in lista),\n\
                                               -nints is the number of intervals (default is 2*sqrt(#number of points in lista)),\n\
                                               -dell is the gap from the begining of the interval where the coordinate of the midle point should be (default is 0.5)\n\
                                               -accu should be 0, 1 or -1 for distribution, right or left cummulative distributions, respectively (default is 0),\n\
                                               -norm should be 0, 1 or -1 for non-normalized, normalized or frequency distribution, respectively (default is 0).\n\
                                              The function returns two lists, x and y with the x and y coordinates of the points in the histogram.";
static char Function_stats_docstring[] = "Evaluates the statisticals for a distribution\n\
                                          usage: [M, S] = stats(lista)      or\n\
                                                 [M, S, t] = stats(lista, assymetry=True)  \n\
                                          returns M, S and (optionally) t, the average, standard deviation and skewness of the numbers in lista.";
static char Function_stats_his_docstring[] = "Evaluates the statisticals for an histogram\n\
                                          usage: [M, S] = stats_his(xi, yi, assymetry=True, norm=0, N=0)      or\n\
                                                 [M, S, t] = stats_his(xi, yi, 1)  \n\
                                          returns M, S and (optionally) t, the average, standard deviation and skewness of the numbers in lista.";
static char Function_drand_exp_docstring[] = "Generates a random number with exponential distribution\n\
                                              usage: exp_rvs(tau)";
static char Function_drand_gauss_docstring[] = "Generates a random number with gaussian distribution\n\
                                                usage: gauss_rvs(mu, sig)";
static char Function_drand_exg_docstring[] = "Generates a random number with ex-gaussian distribution\n\
                                              usage: exg_rvs(mu, sig, tau)";
static char Function_gaussian_docstring[] = "Evaluates the gaussian density function\n\
                                             usage: gauss_pdf(x, mu, sig)\n\
                                             returns the value of the gaussian density function at point x.";
static char Function_exgauss_docstring[] = "Evaluates the ex-gaussian density function\n\
                                            usage: exg_pdf(x, mu, sig, tau)\n\
                                            returns the values of the ex-gaussian density function at point x.";
static char Function_exgauss_lamb_docstring[] = "Evaluates the ex-gaussian density function in terms of the assymetry lambda\n\
                                                 usage: exg_lamb_pdf(x, lamb)\n\
                                                 returns the values of the ex-gaussian density function parametrized in terms of its assymetry at point x.";
static char Function_pars_to_stats_docstring[] = "converts parameters to statisticals\n\
                                                  usage: [M, S, t] = pars_to_stats(mu, sig, tau)";
static char Function_stats_to_pars_docstring[] = "converts statisticals to parameters\n\
                                                  usage: [mu, sig, tau] = stats_to_pars(M, S, t)";
static char Function_int_points_gauss_docstring[] = "Generates gauss partition for an interval\n\
                                                     usage: int_points_gauss(ini, fin, N=N)\n\
                                                     generates a list of numbers between ini and fin for gaussian integration with 20 points.\n\
                                                     -N is the number of gaussian intervals (20 by default).\n\
                                                     The returned list will have 20*N elements.";
static char Function_intsum_docstring[] = "Integral by gauss quadrature.\n\
                                           usage: intsum(ini, fin, y)\n\
                                           calculates an integral. y should be a list of numbers with the function values calculated at \n\
                                           the partition points x given by a int_point_gauss call. Ex (given a function func(x) is deffined):\n\
                                           x = int_point_gauss(-10, 10) # creates a partition\n\
                                           y = [func(ele) for ele in x] # evaluates the function at the points\n\
                                           integral = intsum(-10, 10, y) # evaluates the integral of func(x) between -10 and 10.";
static char Function_exgLKHD_docstring[] = "Evaluates the Log likelyhood and its gradient for the points in the list given the ex gaussian parameters.\n\
                                            Usage: [logLKHD, grad] = exgLKHD(lista, mu, sig, tau)";
static char Function_maxLKHD_docstring[] = "Finds the parameters that maximize the likelyhood for the values in a list.\n\
                                            Usage: [mu, sig, tau] = maxLKHD(lista, mu=0., sig=0., tau=0., lamb=1., eps=1.e-10)\n\
                                            The first three optional arguments are the initial point for searching (actually if the three \n\
                                            are zero, it evaluates mu, sig and tau with respect to the distribution average, standard deviation \n\
                                            and assymmetry), lamb is the the starting update in the parameters with respect to the gradient module\n\
                                            and eps is the precision (the interactions stop when the absolute value of the gradient \n\
                                            times lamb is smaller than eps).";
static char Function_exgSQR_docstring[] = "Sum of the squares (yi-f(xi))**2 for the set of points and its gradient.\n\
                                           Usage: [chi2, grad] = exgSQR(listax, listay, mu, sig, tau)";
static char Function_minSQR_docstring[] = "Finds the parameters that minimize the sum of the squares. \n\
                                           Usage: [mu, sig, tau] = minSQR(listax, listay, mu=0., sig=0., tau=0., lamb=1., eps=1.e-10, NN=0)\n\
                                           The first three optional arguments are the initial point for searching (actually if the three \n\
                                           are zero, it evaluates mu, sig and tau with respect to the distribution average, standard deviation \n\
                                           and assymmetry), lamb is the the starting update in the parameters with respect to the gradient module\n\
                                           and eps is the precision (the interactions stop when the absolute value of the gradient \n\
                                           times lamb is smaller than eps).\n\
                                           Note that listax and listay should represent a normalized histogram of the data!";
static char Function_correlation_docstring[] = "Evaluates the linear correlation coefficient between the numbers in two lists.\n\
                                                Usage: corr = correlation(lista1, lista2)";
static char Function_minsquare_docstring[] = "Adjust a set of points to fit a pollynomial of any degree by the minimum square method.\n\
                                              Usage: [coefs_list, chi2] = minsquare(listax, listay, err=[], deg=1)\n\
                                              listax and listay are lists with the x and y values of the points, respectivelly. If the \n\
                                              keyword argument err is given, it should be a list containing the uncertainties (error bars) in \n\
                                              the y points. the optional argument deg is the degree of the pollynomial to be fitted (1 \n\
                                              means a line). The routine returns coefs_list which is a list with the coefficients of the \n\
                                              fitted pollynomial:\n\
                                                 y = coefs_list[0] + coefs_list[1]* x + coefs_list[2]* x**2 + ... \n\
                                              chi2 is the value of the sum (y_i - y[x_i])**2 / err_i, in other words, it's the minimized sum.";
static char Function_phig_docstring[] = "Evaluates the gaussian cumulative distribution\n\
                                            usage: gauss_cdf(x, mu, sig)\n\
                                            returns the values of the gaussian cumulative distribution at point x.";
static char Function_exgCDF_docstring[] = "Evaluates the ex-gaussian cumulative distribution\n\
                                            usage: exg_cdf(x, mu, sig, tau)\n\
                                            returns the values of the ex-gaussian cumulative distribution at point x.";
static char Function_exg_q_docstring[] = "Evaluates the position of the quartile Q for the exgaussian distribution\n\
                                            usage: exg_ppf(Q, mu, sig, tau, eps=1.e-14, xi=0.)\n\
                                            returns the position where the ex-gaussian distribution accumulates an amount Q.\n\
                                            xi is the point where newton's methos starts the search (if 0. xi=M).";
static char Function_exgCDF_lamb_docstring[] = "Evaluates the ex-gaussian cumulative distribution\n\
                                            usage: exg_lamb_cdf(x, lamb)\n\
                                            returns the values of the ex-gaussian cumulative distribution at point x.";
static char Function_exg_q_lamb_docstring[] = "Evaluates the position of the quartile Q for the exgaussian distribution\n\
                                            usage: exg_lamb_ppf(Q, lamb, eps=1.e-14, xi=0.)\n\
                                            returns the position where the ex-gaussian distribution accumulates an amount Q.\n\
                                            xi is the point where newton's methos starts the search (if 0. xi=M).";


// *******************************************
/*             Module functions             */
// *******************************************
/* function header */
static PyObject *Function_drand(PyObject *self);
static PyObject *Function_histogram(PyObject *self, PyObject *args, PyObject *kwds);
static PyObject *Function_stats(PyObject *self, PyObject *args, PyObject *kwds);
static PyObject *Function_stats_his(PyObject *self, PyObject *args, PyObject *kwds);
static PyObject *Function_drand_exp(PyObject *self, PyObject *args);
static PyObject *Function_drand_gauss(PyObject *self, PyObject *args, PyObject *kwds);
static PyObject *Function_drand_exg(PyObject *self, PyObject *args);
static PyObject *Function_gaussian(PyObject *self, PyObject *args);
static PyObject *Function_exgauss(PyObject *self, PyObject *args);
static PyObject *Function_exgauss_lamb(PyObject *self, PyObject *args);
static PyObject *Function_pars_to_stats(PyObject *self, PyObject *args);
static PyObject *Function_stats_to_pars(PyObject *self, PyObject *args);
static PyObject *Function_int_points_gauss(PyObject *self, PyObject *args, PyObject *kwds);
static PyObject *Function_intsum(PyObject *self, PyObject *args);
static PyObject *Function_exgLKHD(PyObject *self, PyObject *args);
static PyObject *Function_maxLKHD(PyObject *self, PyObject *args, PyObject *kwds);
static PyObject *Function_exgSQR(PyObject *self, PyObject *args);
static PyObject *Function_minSQR(PyObject *self, PyObject *args, PyObject *kwds);
static PyObject *Function_correlation(PyObject *self, PyObject *args);
static PyObject *Function_minsquare(PyObject *self, PyObject *args, PyObject *kwds);
static PyObject *Function_phig(PyObject *self, PyObject *args);
static PyObject *Function_exgCDF(PyObject *self, PyObject *args);
static PyObject *Function_exg_q(PyObject *self, PyObject *args, PyObject *kwds);
static PyObject *Function_exgCDF_lamb(PyObject *self, PyObject *args);
static PyObject *Function_exg_q_lamb(PyObject *self, PyObject *args, PyObject *kwds);
/* Functions that will be callable from python */
static PyMethodDef module_methods[] = {
    {"drand", (PyCFunction)Function_drand, METH_NOARGS, Function_drand_docstring},
    {"histogram", (PyCFunction)Function_histogram, METH_VARARGS | METH_KEYWORDS, Function_histogram_docstring},
    {"stats", (PyCFunction)Function_stats, METH_VARARGS | METH_KEYWORDS, Function_stats_docstring},
    {"stats_his", (PyCFunction)Function_stats_his, METH_VARARGS | METH_KEYWORDS, Function_stats_his_docstring},
    {"exp_rvs", (PyCFunction)Function_drand_exp, METH_VARARGS, Function_drand_exp_docstring},
    {"gauss_rvs", (PyCFunction)Function_drand_gauss, METH_VARARGS | METH_KEYWORDS, Function_drand_gauss_docstring},
    {"exg_rvs", (PyCFunction)Function_drand_exg, METH_VARARGS, Function_drand_exg_docstring},
    {"gauss_pdf", (PyCFunction)Function_gaussian, METH_VARARGS, Function_gaussian_docstring},
    {"exg_pdf", (PyCFunction)Function_exgauss, METH_VARARGS, Function_exgauss_docstring},
    {"exg_lamb_pdf", (PyCFunction)Function_exgauss_lamb, METH_VARARGS, Function_exgauss_lamb_docstring},
    {"pars_to_stats", (PyCFunction)Function_pars_to_stats, METH_VARARGS, Function_pars_to_stats_docstring},
    {"stats_to_pars", (PyCFunction)Function_stats_to_pars, METH_VARARGS, Function_stats_to_pars_docstring},
    {"int_points_gauss", (PyCFunction)Function_int_points_gauss, METH_VARARGS | METH_KEYWORDS, Function_int_points_gauss_docstring},
    {"intsum", (PyCFunction)Function_intsum, METH_VARARGS, Function_intsum_docstring},
    {"exgLKHD", (PyCFunction)Function_exgLKHD, METH_VARARGS, Function_exgLKHD_docstring},
    {"maxLKHD", (PyCFunction)Function_maxLKHD, METH_VARARGS | METH_KEYWORDS, Function_maxLKHD_docstring},
    {"exgSQR", (PyCFunction)Function_exgSQR, METH_VARARGS, Function_exgSQR_docstring},
    {"minSQR", (PyCFunction)Function_minSQR, METH_VARARGS | METH_KEYWORDS, Function_minSQR_docstring},
    {"correlation", (PyCFunction)Function_correlation, METH_VARARGS, Function_correlation_docstring},
    {"minsquare", (PyCFunction)Function_minsquare, METH_VARARGS | METH_KEYWORDS, Function_minsquare_docstring},
    {"gauss_cdf", (PyCFunction)Function_phig, METH_VARARGS, Function_phig_docstring},
    {"exg_cdf", (PyCFunction)Function_exgCDF, METH_VARARGS, Function_exgCDF_docstring},
    {"exg_ppf", (PyCFunction)Function_exg_q, METH_VARARGS | METH_KEYWORDS, Function_exg_q_docstring},
    {"exg_lamb_cdf", (PyCFunction)Function_exgCDF_lamb, METH_VARARGS, Function_exgCDF_lamb_docstring},
    {"exg_lamb_ppf", (PyCFunction)Function_exg_q_lamb, METH_VARARGS | METH_KEYWORDS, Function_exg_q_lamb_docstring},
    {NULL, NULL, 0, NULL}
};


// *******************************************
/*           Initialization                 */
// *******************************************
PyMODINIT_FUNC inituts(void) {
    PyObject* mod;
    // Create the module
    mod = Py_InitModule3("uts", module_methods, utils_docstring);
    if (mod == NULL) {
       return;
    }
}




// *******************************************
/*              Functions                   */
// *******************************************



static PyObject *Function_drand(PyObject *self) {
    double num;    
    num = drand();        
    return Py_BuildValue("d", num);
}








static PyObject *Function_histogram(PyObject *self, PyObject *args, PyObject *kwds) {   
    int i, N, Nint=-1, accu=0, norm=0;
    double *lista;
    double ini=-999999999999999., fin=999999999999999., dell=.5, eps=1.e-10;
    double *xi, *ni;
    PyObject *Numbers, *XX, *YY;

    // Parse the input
    static char *kwlist[] = {"lista", "ini", "fin", "Nint", "dell", "accu", "norm", "eps", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|ddidiid", kwlist, &Numbers, &ini, &fin, &Nint, &dell, &accu, &norm, &eps))
        return 0;
    
    N = (int) PyList_Size(Numbers);
    lista = PyMem_New(double, N);
    for (i=0; i<N; i++) {
        *(lista+i) = (double) PyFloat_AsDouble(PyList_GetItem(Numbers, i));
    }
    // Decides undefined input
    if (Nint<=0) {
        Nint = 2*sqrt(N);
    }
    if (ini==-999999999999999.) {
        ini = min_list(lista, N);
    }
    if (fin==999999999999999.) {
        fin = max_list(lista, N);
    }
    fin += eps;
    // Histogram!
    xi = PyMem_New(double, Nint);
    ni = PyMem_New(double, Nint);
    histogram(lista, N, ini, fin, Nint, accu, norm, dell, xi, ni);
    
    // transforms in python objs
    XX = PyList_New(Nint);
    YY = PyList_New(Nint);
    for (i=0; i<Nint; i++) {
        PyObject *num = PyFloat_FromDouble(*(xi+i));
        PyList_SetItem(XX, i, num);
        PyObject *num2 = PyFloat_FromDouble(*(ni+i));
        PyList_SetItem(YY, i, num2);
    }
    
    PyMem_Free(lista);
    PyMem_Free(xi);
    PyMem_Free(ni);
    PyObject *toret = Py_BuildValue("[O, O]", XX, YY);
    Py_XDECREF(XX);
    Py_XDECREF(YY);
    return toret;
}










static PyObject *Function_stats(PyObject *self, PyObject *args, PyObject *kwds) {
    int N;
    PyObject *Numbers;
    int i, assy=0;
    double *nums;
    double xb, sb, t;
    
    // Parse the input
    static char *kwlist[] = {"numbers", "assymetry", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|i", kwlist, &Numbers, &assy))
        return 0;

    N = (int) PyList_Size(Numbers);
    nums = PyMem_New(double, N);
    for (i=0;i<N;i++) {
        *(nums + i) = (double) PyFloat_AsDouble(PyList_GetItem(Numbers, i));
    }
        
    if (N==0) {
        printf("Warning: Calculating average of an empty list!\n");
    }
    if (N==1) {
        printf("Warning: Calculating standard deviation of a single element!\n");
    }

    if (assy) {
        stats(nums, N, &xb, &sb, &t, 1);
        PyMem_Free(nums);
        return Py_BuildValue("ddd", xb, sb, t);
    } else {
        stats(nums, N, &xb, &sb, NULL, 0);
        PyMem_Free(nums);
        return Py_BuildValue("dd", xb, sb);
    }
}






static PyObject *Function_stats_his(PyObject *self, PyObject *args, PyObject *kwds) {
    int N;
    PyObject *XX, *YY;
    int i, assy=0, norm=0, NN=0;
    double *xi, *yi;
    double xb, sb, t;
    
    // Parse the input
    static char *kwlist[] = {"xi", "yi", "assymetry", "norm", "N", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO|iii", kwlist, &XX, &YY, &assy, &norm, &NN))
        return 0;

    N = (int) PyList_Size(XX);
    xi = PyMem_New(double, N);
    yi = PyMem_New(double, N);
    for (i=0;i<N;i++) {
        *(xi + i) = (double) PyFloat_AsDouble(PyList_GetItem(XX, i));
        *(yi + i) = (double) PyFloat_AsDouble(PyList_GetItem(YY, i));
    }

    if (N==0) {
        printf("Warning: Calculating average of an empty list!\n");
    }
    if (N==1) {
        printf("Warning: Calculating standard deviation of a single element!\n");
    }

    if (assy) {
        stats_his(xi, yi, N, &xb, &sb, &t, 1, norm, NN);
        PyMem_Free(xi);
        PyMem_Free(yi);
        return Py_BuildValue("ddd", xb, sb, t);
    } else {
        stats_his(xi, yi, N, &xb, &sb, NULL, 0, norm, NN);
        PyMem_Free(xi);
        PyMem_Free(yi);
        return Py_BuildValue("dd", xb, sb);
    }
}








static PyObject *Function_drand_exp(PyObject *self, PyObject *args) {
    double tau;
    
    if (!PyArg_ParseTuple(args, "d", &tau))
        return NULL;
    return Py_BuildValue("d", drand_exp(tau));
}






static PyObject *Function_drand_gauss(PyObject *self, PyObject *args, PyObject *kwds) {   
    double mu=0., sig=1.;
    
    // Parse the input
    static char *kwlist[] = {"mu", "sig", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dd", kwlist, &mu, &sig))
        return 0;
    
    return Py_BuildValue("d", drand_gauss(mu, sig));
}


static PyObject *Function_drand_exg(PyObject *self, PyObject *args) {
    double mu, sig, tau;
    
    if (!PyArg_ParseTuple(args, "ddd", &mu, &sig, &tau))
        return NULL;
    return Py_BuildValue("d", drand_exp(tau)+drand_gauss(mu, sig));
}



static PyObject *Function_gaussian(PyObject *self, PyObject *args) {
    double x, mu, sig;
    
    if (!PyArg_ParseTuple(args, "ddd", &x, &mu, &sig))
        return NULL;
    return Py_BuildValue("d", gaussian(x, mu, sig));
}


static PyObject *Function_exgauss(PyObject *self, PyObject *args) {
    double x, mu, sig, tau;
    
    if (!PyArg_ParseTuple(args, "dddd", &x, &mu, &sig, &tau))
        return NULL;
    return Py_BuildValue("d", exgauss(x, mu, sig, tau));
}


static PyObject *Function_exgauss_lamb(PyObject *self, PyObject *args) {
    double x, lamb;
    
    if (!PyArg_ParseTuple(args, "dd", &x, &lamb))
        return NULL;
    if ((lamb<=0.)|(lamb>=1.)) {
        printf("Warning: lambda must be in the interval (0; 1)!\n");
        return Py_BuildValue("");
    }
    return Py_BuildValue("d", exgauss_lamb(x, lamb));
}

static PyObject *Function_pars_to_stats(PyObject *self, PyObject *args) {
    double mu, sig, tau;
    double M, S, lamb, t;
    
    if (!PyArg_ParseTuple(args, "ddd", &mu, &sig, &tau))
        return NULL;
    pars_to_stats(mu, sig, tau, &M, &S, &lamb);
    t = 2.*lamb*lamb*lamb;
    return Py_BuildValue("ddd", M, S, t);
}

static PyObject *Function_stats_to_pars(PyObject *self, PyObject *args) {
    double mu, sig, tau;
    double M, S, lamb, t;
    
    if (!PyArg_ParseTuple(args, "ddd", &M, &S, &t))
        return NULL;
    if ((t<=0.)|(t>=2.)) {
        printf("Warning: t must be in the interval (0; 2)!\n");
        //return Py_BuildValue("");
    }
    lamb = pow(0.5*t, 1./3);
    stats_to_pars(M, S, lamb, &mu, &sig, &tau);
    return Py_BuildValue("ddd", mu, sig, tau);
}





static PyObject *Function_int_points_gauss(PyObject *self, PyObject *args, PyObject *kwds) {   
    int i, N=20;
    double ini, fin;
    double *xxx;
    PyObject *xxs;

    // Parse the input
    static char *kwlist[] = {"ini", "fin", "N", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "dd|i", kwlist, &ini, &fin, &N))
        return 0;
    
    xxx = PyMem_New(double, 20*N);
    gen_points(ini, fin, N, xxx);
        
    xxs = PyList_New(20*N);
    for (i=0; i<20*N; i++) {
        PyObject *num = PyFloat_FromDouble(*(xxx + i));
        PyList_SetItem(xxs, i, num);
    }
    
    PyMem_Free(xxx);
    PyObject *toret = Py_BuildValue("O", xxs);
    Py_XDECREF(xxs);
    return toret;
}



static PyObject *Function_intsum(PyObject *self, PyObject *args) {
    double ini, fin, integrand;
    int i, N, Np;
    PyObject *yyy;
    double *ys;
    
    if (!PyArg_ParseTuple(args, "ddO", &ini, &fin, &yyy))
        return NULL;
    Np = (int) PyList_Size(yyy);
    if (Np%20) {
        printf("Warning: bad gauss partition!\n");
        return Py_BuildValue("");
    }
    N = Np/20;
    ys = PyMem_New(double, Np);
    for (i=0; i<Np; i++) {
        *(ys + i) = (double) PyFloat_AsDouble(PyList_GetItem(yyy, i));
    }
    integrand = integrate(ini, fin, N, ys);
    
    PyMem_Free(ys);
    return Py_BuildValue("d", integrand);
}




static PyObject *Function_exgLKHD(PyObject *self, PyObject *args) {
    double mu, sig, tau, lnL;
    int i, N;
    PyObject *xxx;
    double *xi, *grad;
    
    if (!PyArg_ParseTuple(args, "Oddd", &xxx, &mu, &sig, &tau))
        return NULL;
    N = (int) PyList_Size(xxx);
    xi = PyMem_New(double, N);
    for (i=0; i<N; i++) {
        *(xi + i) = (double) PyFloat_AsDouble(PyList_GetItem(xxx, i));
    }
    grad = PyMem_New(double, 3);
    exgLKHD(xi, N, mu, sig, tau, &lnL, grad);
    //printf("  lnl=   %10.10f\n", lnL);
    
    PyObject *vec = PyList_New(3);
    for (i=0; i<3; i++) {
        PyObject *num = PyFloat_FromDouble(*(grad+i));
        PyList_SetItem(vec, i, num);
    }
    
    PyMem_Free(xi);
    PyMem_Free(grad);
    return Py_BuildValue("dO", lnL, vec);
}


static PyObject *Function_maxLKHD(PyObject *self, PyObject *args, PyObject *kwds) {   
    double mu=0., sig=0., tau=0.;
    double M, S, t, lambi, lnL;
    double lamb=1., eps=1.e-16;
    PyObject *xxx;
    double *xi, *grad;
    //double lnL;
    int i, N;
    
    // Parse the input
    static char *kwlist[] = {"xi", "mu", "sig", "tau", "lamb", "eps", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|ddddd", kwlist, &xxx, &mu, &sig, &tau, &lamb, &eps))
        return 0;
    N = (int) PyList_Size(xxx);
    xi = PyMem_New(double, N);
    for (i=0; i<N; i++) {
        *(xi + i) = (double) PyFloat_AsDouble(PyList_GetItem(xxx, i));
    }
    // determines the stats and initial parameters
    if (!mu & !sig & !tau) {
        stats(xi, N, &M, &S, &t, 1);
        lambi = pow(t*0.5, 1./3.);
        if (t<0.) {
            //printf("Warning: Distribution has the wrong sign for the skewness!\n");
            lambi = 0.2;
        } else if (lambi>=1.) {
            //printf("Warning: Distribution has lambda > 1!\n");
            lambi = 0.99;
        }
        stats_to_pars(M, S, lambi, &mu, &sig, &tau);
    }
    if (!lamb) {
        grad = PyMem_New(double, 3);
        exgLKHD(xi, N, mu, sig, tau, &lnL, grad);
        lamb = 5./pow(grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2], 0.5);
    } else if (lamb>0) {
        grad = PyMem_New(double, 3);
        exgLKHD(xi, N, mu, sig, tau, &lnL, grad);
        lamb = lamb/pow(grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2], 0.5);
    } else {
        lamb = -lamb;
    }
    // steepest descent
    grad = PyMem_New(double, 3);
    //exgLKHD(xi, N, mu, sig, tau, &lnL, grad);
    //printf("mu=%f   sig=%f    tau=%f\n", mu, sig, tau);
    maxLKHD(xi, N, &mu, &sig, &tau, lamb, eps);
    //printf("mu=%f   sig=%f    tau=%f\n", mu, sig, tau);

    PyMem_Free(xi);
    PyMem_Free(grad);    
    return Py_BuildValue("ddd", mu, sig, tau);
}





static PyObject *Function_exgSQR(PyObject *self, PyObject *args) {
    double mu, sig, tau, lnL;
    int i, N;
    PyObject *xxx, *yyy;
    double *xi, *yi, *grad;
    
    if (!PyArg_ParseTuple(args, "OOddd", &xxx, &yyy, &mu, &sig, &tau))
        return NULL;
    N = (int) PyList_Size(xxx);
    xi = PyMem_New(double, N);
    yi = PyMem_New(double, N);
    for (i=0; i<N; i++) {
        *(xi + i) = (double) PyFloat_AsDouble(PyList_GetItem(xxx, i));
        *(yi + i) = (double) PyFloat_AsDouble(PyList_GetItem(yyy, i));
    }
    grad = PyMem_New(double, 3);
    exgSQR(xi, yi, N, mu, sig, tau, &lnL, grad);
    
    PyObject *vec = PyList_New(3);
    for (i=0; i<3; i++) {
        PyObject *num = PyFloat_FromDouble(*(grad+i));
        PyList_SetItem(vec, i, num);
    }
    
    PyMem_Free(xi);
    PyMem_Free(yi);
    PyMem_Free(grad);
    return Py_BuildValue("dO", lnL, vec);
}


static PyObject *Function_minSQR(PyObject *self, PyObject *args, PyObject *kwds) {   
    double mu=0., sig=0., tau=0.;
    double M, S, t, lambi, lnL, dx, suma=0.;
    double lamb=0., eps=1.e-16;
    PyObject *xxx, *yyy;
    double *xi, *yi, *grad;
    int i, N, NN=0;
    
    // Parse the input
    static char *kwlist[] = {"xi", "yi", "mu", "sig", "tau", "lamb", "eps", "N", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO|dddddi", kwlist, &xxx, &yyy, &mu, &sig, &tau, &lamb, &eps, &NN))
        return 0;
    N = (int) PyList_Size(xxx);
    xi = PyMem_New(double, N);
    yi = PyMem_New(double, N);
    for (i=0; i<N; i++) {
        *(xi + i) = (double) PyFloat_AsDouble(PyList_GetItem(xxx, i));
        *(yi + i) = (double) PyFloat_AsDouble(PyList_GetItem(yyy, i));
        suma += *(yi + i);
    }
    // determines the stats and initial parameters (I think it is ok now.)
    if (!mu & !sig & !tau) {
        if (~NN){
            NN = (int) N*N/4.+ N/4.;
        }
        stats_his(xi, yi, N, &M, &S, &t, 1, 1, NN);
        lambi = pow(t*0.5, 1./3.);
        if (t<0.) {
            lambi = 0.2;
        } else if (lambi>=1.) {
            lambi = 0.99;
        }
        stats_to_pars(M, S, lambi, &mu, &sig, &tau);
    }
    // Normalizes (integral must be one!)
    dx = *(xi+1) - *xi;
    suma *= dx;
    for (i=0; i<N; i++) {
        *(yi + i) /= suma;
    }
    // selects lamb
    if (!lamb) {
        grad = PyMem_New(double, 3);
        exgSQR(xi, yi, N, mu, sig, tau, &lnL, grad);
        lamb = 5./pow(grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2], 0.5);
        //printf("%f %d\n", lamb, NN);
    } else if (lamb>0) {
        grad = PyMem_New(double, 3);
        exgSQR(xi, yi, N, mu, sig, tau, &lnL, grad);
        lamb = lamb/pow(grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2], 0.5);
    } else {
        lamb = -lamb;
    }
    // steepest descent
    grad = PyMem_New(double, 3);
    minSQR(xi, yi, N, &mu, &sig, &tau, lamb, eps);
    
    PyMem_Free(xi);
    PyMem_Free(yi);
    PyMem_Free(grad);
    return Py_BuildValue("ddd", mu, sig, tau);
}













static PyObject *Function_correlation(PyObject *self, PyObject *args) {
    int N, N2;
    PyObject *xs1, *xs2;
    int i;
    double *nums1, *nums2;
    double corr;
    
    if (!PyArg_ParseTuple(args, "OO", &xs1, &xs2))
        return NULL;
    N = (int) PyList_Size(xs1);
    N2 = (int) PyList_Size(xs2);
    if (N!=N2) {
        printf("Warning: Lists do not have the same size!!!\n");
    }
    nums1 = PyMem_New(double, N);
    nums2 = PyMem_New(double, N);
    for (i=0; i<N; i++) {
        *(nums1 + i) = (double) PyFloat_AsDouble(PyList_GetItem(xs1, i));
        *(nums2 + i) = (double) PyFloat_AsDouble(PyList_GetItem(xs2, i));
    }
    corr = correlation(nums1, nums2, N);
    PyMem_Free(nums1);
    PyMem_Free(nums2);
    return Py_BuildValue("d", corr);
}















static PyObject *Function_minsquare(PyObject *self, PyObject *args, PyObject *kwds) {   
    int i, deg=1, N, Nt;
    double *xs, *ys, *sigs, *alphas;
    double chi2;
    PyObject *X, *Y, *err=Py_BuildValue("[]"), *alp;

    // Parse the input
    static char *kwlist[] = {"X", "Y", "err", "deg", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO|Oi", kwlist, &X, &Y, &err, &deg))
        return 0;
    
    N = (int) PyList_Size(X);
    Nt = (int) PyList_Size(err);
    
    xs = PyMem_New(double, N);
    ys = PyMem_New(double, N);
    alphas = PyMem_New(double, deg+1);
    if (Nt!=0) {
        sigs = PyMem_New(double, N);   
        for (i=0; i<N; i++) {
            *(xs+i) = (double) PyFloat_AsDouble(PyList_GetItem(X, i));
            *(ys+i) = (double) PyFloat_AsDouble(PyList_GetItem(Y, i));
            *(sigs+i) = (double) PyFloat_AsDouble(PyList_GetItem(err, i));
        }
        chi2 = minsquare(xs, ys, sigs, N, deg, alphas);
        PyMem_Free(sigs);
    } else {
        for (i=0; i<N; i++) {
            *(xs+i) = (double) PyFloat_AsDouble(PyList_GetItem(X, i));
            *(ys+i) = (double) PyFloat_AsDouble(PyList_GetItem(Y, i));
        }
        chi2 = minsquare(xs, ys, NULL, N, deg, alphas);
    }
    alp = PyList_New(deg+1);
    for (i=0; i<deg+1; i++) {
        PyObject *num = PyFloat_FromDouble(*(alphas+i));
        PyList_SetItem(alp, i, num);
    }
    
    PyMem_Free(xs);
    PyMem_Free(ys);
    PyMem_Free(alphas);
    PyObject *toret = Py_BuildValue("[O, d]", alp, chi2);
    Py_XDECREF(alp);
    return toret;
}



static PyObject *Function_phig(PyObject *self, PyObject *args) {
    double x, mu, sigma;
    double num;
    
    if (!PyArg_ParseTuple(args, "ddd", &x, &mu, &sigma))
        return NULL;
    num = phig(x, mu, sigma);
    return Py_BuildValue("d", num);
}

static PyObject *Function_exgCDF(PyObject *self, PyObject *args) {
    double x, mu, sigma, tau;
    double num;
    
    if (!PyArg_ParseTuple(args, "dddd", &x, &mu, &sigma, &tau))
        return NULL;
    num = exgCDF(x, mu, sigma, tau);
    return Py_BuildValue("d", num);
}



static PyObject *Function_exg_q(PyObject *self, PyObject *args, PyObject *kwds) {   
    double Q, mu, sigma, tau;
    double num;
    double eps=1.e-14, xi=0.;

    // Parse the input
    static char *kwlist[] = {"Q", "mu", "sigma", "tau", "eps", "xi", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "dddd|dd", kwlist, &Q, &mu, &sigma, &tau, &eps, &xi))
        return 0;
    num = exg_q(Q, mu, sigma, tau, eps, xi);
    return Py_BuildValue("d", num);
}


static PyObject *Function_exgCDF_lamb(PyObject *self, PyObject *args) {
    double x, lamb;
    double num;
    
    if (!PyArg_ParseTuple(args, "dd", &x, &lamb))
        return NULL;
    num = exgCDF_lamb(x, lamb);
    return Py_BuildValue("d", num);
}



static PyObject *Function_exg_q_lamb(PyObject *self, PyObject *args, PyObject *kwds) {   
    double Q, lamb;
    double num;
    double eps=1.e-14, xi=0.;

    // Parse the input
    static char *kwlist[] = {"Q", "lamb", "eps", "xi", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "dd|dd", kwlist, &Q, &lamb, &eps, &xi))
        return 0;
    num = exg_q_lamb(Q, lamb, eps, xi);
    return Py_BuildValue("d", num);
}


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

void histogram(double *lista, int N, double ini, double fin, double Nint, int accu, int norm, double dell, double *xi, double *ni); // makes an histogram of elements in lista
void stats(double *nums, int N, double *xb, double *sig, double *t, int assy); // evaluate average and sigma and t for a distribution
void stats_his(double *xx, double *yi, int NN, double *xb, double *sig, double *t, int assy, int norm, int NNN); // evaluate average and sigma and t for an histogram
double gaussian(double x, double mu, double sig); // gaussian distribution
double exgauss(double x, double mu, double sig, double tau); // ex-gaussian function
double exgauss_lamb(double z, double lamb); // ex-gaussian in terms of the assimetry lambda
void pars_to_stats(double mu, double sig, double tau, double *M, double *S, double *lamb); // converts parameters to statisticals
void stats_to_pars(double M, double S, double lamb, double *mu, double *sig, double *tau); // converts statisticals to parameters
void exgLKHD(double *xi, int N, double mu, double sig, double tau, double *lnL, double *grad); // evaluates the likelyhood and its gradient
void maxLKHD(double *xi, int N, double *mu, double *sig, double *tau, double lamb, double eps); // steepest descent to maximize the log likelyhood
void exgSQR(double *xi, double *yi, int N, double mu, double sig, double tau, double *lnL, double *grad); // evaluates the sum of squares and its gradient
void minSQR(double *xi, double *yi, int N, double *mu, double *sig, double *tau, double lamb, double eps); // steepest descent to minimize the sum of squares

double poly(double *a, double x, int N); // Evalueates polynomial a0 + a1*x + a2*x^2 + ...
double correlation(double *xs1, double *xs2, int N); // evaluates the correlation betweeen two lists
double minsquare(double *xs, double *ys, double *err, int N, int deg, double *as); // minimum square method
double min_list(double *lista, int N); // Finds the minimum value of a list
double max_list(double *lista, int N); // Finds the maximum value of a list

double phig(double x, double mu, double sig); // gaussian cumulative distr
double exgCDF(double x, double mu, double sig, double tau); // exgaussian cumulative distr
double exg_q(double Q, double mu, double sig, double tau, double eps, double xi); // finds a quartile
double exgCDF_lamb(double z, double lamb); // exgaussian cumulative distr
double exg_q_lamb(double Q, double lamb, double eps, double xi); // finds a quartile






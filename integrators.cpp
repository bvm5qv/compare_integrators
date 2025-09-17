#include <stdlib.h>
#include <cmath>
#include <stdio.h>
#include "integrators.hpp"
#include "gaussrules.hpp"
#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>
#include <fstream>

using std::cout;
using std::endl;

// Note: this code is not safe wrt the number of intervals
// Make sure the number of intervals used is valid for each integrator


// Integration using trapezoid rule 
// npts : number of points used in calculation (npts>=2)
double trapez (double (*f)(double x), unsigned int npts, double min, double max) 
{
  assert(npts >= 2 && "must have more than two points for trapezoid method");
  const double dx = (max - min) / (npts - 1);

  double sum = 0;
  for(unsigned int i = 1; i <= npts; i++)
  {
    const double x = min + (i - 1)*dx;

    if(i == 1 || i == npts)
      sum += dx/2 * f(x);
    else
      sum += dx * f(x);
  }
  
  return sum;
}      

// Integration using Simpson's rule
// npts : number of points used in calculation (npts odd, and >=3)
double simpson (double (*f)(double x), unsigned npts, double min, double max)
{  
  assert(npts >= 3 && npts % 2 == 1 && "must have odd number of points more than three for Simpson's Rule");
  const double dx = (max - min) / (npts - 1);

  double sum = 0;
  for(unsigned int i = 1; i <= npts; i++)
  {
    const double x = min + (i - 1)*dx;

    if(i == 1 || i == npts)
      sum += dx/3 * f(x);
    else if(i % 2 == 0)
      sum += 4*dx/3 * f(x);
    else
      sum += 2*dx/3 * f(x);
  }
  
  return (sum);
}  

// Integration using Gauss's rule, code is based on the Landau text
// This is not a very good implementation
double gaussint (double (*f)(double x), unsigned npts, double min, double max){
  double result = 0.;
  const unsigned MAXPOINTS=10000;
  double w[MAXPOINTS];    // for points and weights
  double xi[MAXPOINTS];
  
  if (npts>MAXPOINTS) npts=MAXPOINTS;
  gauss (npts, min, max, xi, w);      // returns Legendre polynomials
  // points and weights
  double c1 = (max - min) / 2;
  double c2 = (max + min) / 2;
  for (unsigned n=0; n<npts; n++){
    result += w[n] * f(c1 * xi[n] + c2);   // calculating the integral
  } 
  return result/2;                  
}

// calculate the weights and intervals for the n-point rule
void gauss(unsigned npts, double a, double b, double x[], double w[]){    
  // npts     number of points
  // x, w     output grid interval points and weights.			      

  double  t, t1, pp=0, p1, p2, p3;
  double  eps = 3.e-10;			// limit for accuracy

  // calculating roots of Legendre polynomials and wgt cofficients
  unsigned m = (npts+1)/2;
  for(unsigned i=1; i<=m; i++){  
    t  = cos(M_PI*(i-0.25)/(npts+0.5));
    t1 = 1;
    while(fabs(t-t1)>=eps){ 
      p1 = 1.0;
      p2 = 0.0;
      for(unsigned j=1; j<=npts; j++) {
	p3 = p2;
	p2 = p1;
	p1 = ((2*j-1)*t*p2-(j-1)*p3)/j;
      }
      pp = npts*(t*p1-p2)/(t*t-1);
      t1 = t;
      t  = t1 - p1/pp;
    }   
    x[i-1] = -t;
    x[npts-i] = t;
    w[i-1]    = 2.0/((1-t*t)*pp*pp);
    w[npts-i] = w[i-1];
  } 
}

// modified from rosettacode.org
// this is a MUCH better implemtation of the Gaussian quadrature method
GaussInt::GaussInt(int npoints){
  Init(npoints);
}

#if __GNUC_PREREQ(4,6)
__float128 GaussInt::lege_eval(int n, __float128 x){
  __float128 s = lcoef[n][n];
  for (int i = n; i>0; i--)
    s = s * x + lcoef[n][i - 1];
  return s;
}

__float128 GaussInt::lege_diff(int n, __float128 x){
  return n * (x * lege_eval(n, x) - lege_eval(n - 1, x)) / (x * x - 1);
}
#endif


void GaussInt::Init(int npoints){

  // calculates abscissas and weights to double precision accuracy
  // for n-point quadrature rule
  // Note: In practice you would generally not bother with this step,
  // instead tables of abscissas and weights would be calculated once for
  // a selection of n-point rules and saved for later usage

#if __GNUC_PREREQ(4,6)  
  lroots.assign(npoints,0);
  weight.assign(npoints,0);
  lcoef.assign(npoints+1,vector<__float128>(npoints+1,0));
  
  lcoef[0][0] = lcoef[1][1] = 1;
  for (int n = 2; n <= npoints; n++) {
    lcoef[n][0] = -(n - 1) * lcoef[n - 2][0] / n;
    for (int i = 1; i <= n; i++)
      lcoef[n][i] = ((2 * n - 1) * lcoef[n - 1][i - 1]
		     - (n - 1) * lcoef[n - 2][i] ) / n;
  }

  __float128 x, x1;
  for (int i = 1; i <= npoints; i++) {
    x = cos(M_PI * (i - .25) / (npoints + .5));
    do {
      x1 = x;
      x -= lege_eval(npoints, x) / lege_diff(npoints, x);
    } while (fabs((long double)(x-x1))>1e-16);  // keep double precision 
    lroots[i - 1] = x;
    
    x1 = lege_diff(npoints, x);
    weight[i - 1] = 2 / ((1 - x * x) * x1 * x1);
  }
  // cout << "==== " << npoints << " ====" << endl;
  // cout.precision(20);
  // for (int i = 0; i < npoints; i++)
  //   cout << (double) weight[i] << " " << (double) lroots[i] << endl;
#else
  // if we have an older compiler use the precomputed values
  if (npoints>MAXPOINTRULE){
    npoints=MAXPOINTRULE;
    cout << "Truncating to " << MAXPOINTRULE << "-point rule" << endl;
    cout << "Use gcc 4.6 or higher to enable calculation of higher point rules" << endl;
  }
  double *wa = &GAUSSWA[npoints*npoints-(npoints+2)];
  weight.clear();
  lroots.clear();
  for (int i=0; i<npoints*2; i+=2){
    weight.push_back(wa[i]);
    lroots.push_back(wa[i+1]);
  }
#endif  
}

// notice the efficiency of the integration
// once the constants for the n-point rule are calculated
// the integral is a very simple sum
double GaussInt::Integ(double (*f)(double x), double a, double b)
{
  double sum = 0;

  for(int i = 0; i < this->lroots.size(); i++)
  {
    const double x = (b + a)/2 + (b - a)/2 * this->lroots[i];
    const double w = (b - a)/2 * this->weight[i];

    sum += w*f(x);
  }
  
  return sum;  // integral of function f
}

void GaussInt::PrintWA() const{
  cout << "Weights, abscissas for " << weight.size() << " point rule" << endl;
  for (unsigned i = 0; i < weight.size(); i++){
    cout  << std::setprecision(34) << (long double) weight[i] << "   " << (long double) lroots[i] << endl;
  }
}






double neg_exp(double x) { return std::exp(-x); }
double bad_fcn(double x)
{
  if(std::abs(x) <= 0.5)
    return std::abs(x);
  else
    return -std::abs(x);
}

int main()
{
  std::vector<unsigned int> npts_vec {3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 99, 199, 299, 399, 499, 599, 699, 799, 899, 999};

  //********************************************
  //e^-x integral (0,1) calculation and data recording
  //********************************************
  std::vector<double> trapez_neg_exp_data{};
  std::vector<double> simpson_neg_exp_data{};
  std::vector<double> gauss_neg_exp_data{};
  const double exact_neg_exp_int = 1.0 - neg_exp(1.0);
  for(const auto N : npts_vec)
  {
    const double t_int = trapez(neg_exp, N, 0, 1);
    const double s_int = simpson(neg_exp, N, 0, 1);
    trapez_neg_exp_data.push_back(std::abs((t_int - exact_neg_exp_int)/exact_neg_exp_int));
    simpson_neg_exp_data.push_back(std::abs((s_int - exact_neg_exp_int)/exact_neg_exp_int));

    if(N <= 25)
    {
      GaussInt g(N);
      const double g_int = g.Integ(neg_exp, 0, 1);
      gauss_neg_exp_data.push_back(std::abs((g_int - exact_neg_exp_int)/exact_neg_exp_int));
    }
  }

  std::ofstream neg_exp_data;
  neg_exp_data.open("neg_exp.dat");
  std::cout << "--------------------------------------------\n";
  std::cout << "N        e_t        e_s        e_g    \n";
  std::cout << "--------------------------------------------\n";
  neg_exp_data << "#--------------------------------------------\n";
  neg_exp_data << "#N        e_t        e_s        e_g    \n";
  neg_exp_data << "#--------------------------------------------\n";
  for(int i = 0; i < npts_vec.size(); i++)
  {
    if(i < gauss_neg_exp_data.size())
    {
      std::cout << npts_vec[i] << "    " << trapez_neg_exp_data[i] << "    " << simpson_neg_exp_data[i] << "    " << gauss_neg_exp_data[i] << '\n';
      neg_exp_data << npts_vec[i] << ' ' << trapez_neg_exp_data[i] << ' ' << simpson_neg_exp_data[i] << ' ' << gauss_neg_exp_data[i] << '\n';
    } 
    else
    {
      std::cout << npts_vec[i] << "    " << trapez_neg_exp_data[i] << "    " << simpson_neg_exp_data[i] << "    N/A\n"; 
      neg_exp_data << npts_vec[i] << ' ' << trapez_neg_exp_data[i] << ' ' << simpson_neg_exp_data[i] << '\n';
    } 
  }
  neg_exp_data.close();

  
  //***************************************
  //bad_fcn integral (-1,1) and data recording
  //***************************************
  std::vector<double> bad_trapez_data{};
  std::vector<double> bad_simpson_data{};
  std::vector<double> bad_gauss_data{};
  const double bad_exact_int = -0.5;
  for(const auto N : npts_vec)
  {
    const double t_int = trapez(bad_fcn, N, -1, 1);
    const double s_int = simpson(bad_fcn, N, -1, 1);
    bad_trapez_data.push_back(std::abs((t_int - bad_exact_int)/bad_exact_int));
    bad_simpson_data.push_back(std::abs((s_int - bad_exact_int)/bad_exact_int));
    
    if(N <= 25)
    {
      GaussInt g(N);
      const double g_int = g.Integ(bad_fcn, -1, 1);
      bad_gauss_data.push_back(std::abs((g_int - bad_exact_int)/bad_exact_int));
    }
  }

  std::ofstream bad_x_data;
  bad_x_data.open("bad_x.dat");
  bad_x_data << "#--------------------------------------------\n";
  bad_x_data << "#N        e_t        e_s        e_g    \n";
  bad_x_data << "#--------------------------------------------\n";
  for(int i = 0; i < npts_vec.size(); i++)
  {
    if(i < bad_gauss_data.size())
      bad_x_data << npts_vec[i] << ' ' << bad_trapez_data[i] << ' ' << bad_simpson_data[i] << ' ' << bad_gauss_data[i] << '\n';
    else
      bad_x_data << npts_vec[i] << ' ' << bad_trapez_data[i] << ' ' << bad_simpson_data[i] << '\n'; 
  }
  bad_x_data.close();

  return 0;
}





#include <libcmaes/cmaes.h>
#include <iostream>
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;
using namespace libcmaes;

void transform_to_simplex(const double *vals, double *simplex_point, int len) {
   double cumprod = 1.0;
   for (int i = 0; i < len; i++) {
      simplex_point[i] = cumprod * (1.0 - vals[i]);
      cumprod *= vals[i];
   }
   simplex_point[len] = cumprod;
}

double log_target_dens(const double * vals_inp, const Rcpp::NumericMatrix& lhm) {
   
   const int len = lhm.ncol();
   double ltd = 0.0;
   double vals[len-1];
   
   
   // Logistic transform
   for (int j = 0; j < len - 1; j++) {
      vals[j] = 1.0 / (1.0 + exp(-vals_inp[j]));
   }
   
   // Jacobian term
   for (int j = 0; j < len - 1; j++) {
      ltd += std::log(vals[j] * (1.0 - vals[j]));
   }
   
   // Prior on simplex
   for (int j = 0; j < len - 1; j++) {
      ltd += (len - j - 1) * std::log(vals[j]);
   }
   
   // Map to simplex
   double simplex_point[len+1];
   transform_to_simplex(vals, simplex_point, len);
   
   // Log-likelihood
   for (int i = 0; i < lhm.nrow(); i++) {
      double lh = 0.0;
      for (int j = 0; j < len + 1; j++)
         lh += lhm(i,j) * simplex_point[j];
      ltd += std::log(lh);
   }
   
   return ltd;
}

// [[Rcpp::export]]
double wrapped_log_target_dens(const Rcpp::NumericVector& vals_inp, const Rcpp::NumericMatrix& lhm) {
   return log_target_dens( vals_inp.begin(), lhm );
}
   
libcmaes::FitFunc fsphere = [](const double *x, const int N)
{
   double val = 0.0;
   for (int i=0;i<N;i++)
      val += x[i]*x[i];
   return val;
};

// [[Rcpp::export]]
int test()
{
   int dim = 10; // problem dimensions.
   vector<double> x0(dim,10.0);
   double sigma = 0.1;
   //int lambda = 100; // offsprings at each generation.
   CMAParameters<> cmaparams(x0,sigma);
   //cmaparams._algo = BIPOP_CMAES;
   CMASolutions cmasols = cmaes<>(fsphere,cmaparams);
   cout << "best solution: " << cmasols << endl;
   cout << "optimization took " << cmasols.elapsed_time() / 1000.0 << " seconds\n";
   return cmasols.run_status();
}

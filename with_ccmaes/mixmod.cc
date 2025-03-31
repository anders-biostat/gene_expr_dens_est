#include <Rcpp.h>
#include <armadillo>
#include "c-cmaes/src/cmaes_interface.h"

#include "c-cmaes/src/cmaes.c"

// [[Rcpp::depends(RcppArmadillo)]]

void transform_to_simplex(const arma::vec& vals, arma::vec& simplex_point) {
   const int len = vals.n_elem;
   double cumprod = 1.0;
   for (int i = 0; i < len; i++) {
      simplex_point[i] = cumprod * (1.0 - vals[i]);
      cumprod *= vals[i];
   }
   simplex_point[len] = cumprod;
}

double log_target_dens(const arma::vec& vals_inp, const Rcpp::NumericMatrix& lhm) {
   const int len = vals_inp.n_elem;
   double ltd = 0.0;
   
   // Logistic transform
   arma::vec vals = 1.0 / (1.0 + arma::exp(-vals_inp));
   
   // Jacobian term
   for (int j = 0; j < len - 1; j++) {
      ltd += std::log(vals[j] * (1.0 - vals[j]));
   }
   
   // Prior on simplex
   for (int j = 0; j < len - 1; j++) {
      ltd += (len - j - 1) * std::log(vals[j]);
   }
   
   // Map to simplex
   arma::vec simplex_point(len + 1);
   transform_to_simplex(vals, simplex_point);
   
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
Rcpp::NumericVector mixmod( Rcpp::NumericMatrix lhm ) {
   
   int n = lhm.ncol() - 1;
   double x0[n];
   double sigma0[n];
   double *const *pop;
   double *funvals;
   int iter = 0;
   cmaes_t evo;
   
   std::fill_n( x0, n, 0. );
   std::fill_n( sigma0, n, 1. );
   funvals = cmaes_init( &evo, n, x0, sigma0, 0, 0, NULL );
   evo.sp.stopTolFun = 1e-2;
   
   while( !cmaes_TestForTermination(&evo) ) {
      pop = cmaes_SamplePopulation(&evo);
      
      for( int i = 0; i < evo.sp.lambda; ++i ) {
         arma::vec v( pop[i], n, false, true );
         funvals[i] = -log_target_dens( v, lhm );
         Rcpp::Rcout << funvals[i] << std::endl;
      }
      cmaes_UpdateDistribution( &evo, funvals );  
      iter++;
      
      Rcpp::checkUserInterrupt();
   }
      
   arma::vec v; v.set_size(n);
   for( int i=0; i<n; i++ ) {
      v[i] = 1. / ( 1. + std::exp(-evo.rgxmean[i]) );
   }
   Rcpp::NumericVector ans(n+1);
   arma::vec w( ans.begin(), n+1, false, true);
   transform_to_simplex( v, w );
   for( int i=0; i<n; i++ )
      Rcpp::Rcout << w[i] << " ";
   Rcpp::Rcout << std::endl;
   
   ans.attr("stopping_reason") = cmaes_TestForTermination(&evo);
   ans.attr("population_size") = evo.sp.lambda;
   ans.attr("iterations") = iter;
   
   cmaes_exit( &evo );
   return ans;
}

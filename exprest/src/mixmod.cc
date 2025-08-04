#include <Rcpp.h>
#include "c-cmaes/cmaes_interface.h"

#include "c-cmaes/cmaes.c"

void transform_to_simplex( const double *const vals, double *const simplex_point, int n ) {
   double cumprod = 1.0;
   for (int i = 0; i < n; i++) {
      simplex_point[i] = cumprod * (1.0 - vals[i]);
      cumprod *= vals[i];
   }
   simplex_point[n] = cumprod;
}

double log_target_dens( const double *const vals_inp, 
      const Rcpp::NumericMatrix& lhm, const double penalty_coef=0. ) {
   const int n = lhm.ncol() - 1;
   const int nrow = lhm.nrow();
   double ltd = 0.0;
   
   // Logistic transform
   double vals[n];
   for( int i = 0; i<n; i++ )
      vals[i] = 1.0 / ( 1.0 + std::exp(-vals_inp[i]) );
   
   // Jacobian term
   for (int j = 0; j < n - 1; j++) {
      ltd += std::log(vals[j] * (1.0 - vals[j]));
   }
   
   // Prior on simplex
   for (int j = 0; j < n - 1; j++) {
      ltd += (n - j - 1) * std::log(vals[j]);
   }
   
   // Map to simplex
   double simplex_point[ n+1 ];
   transform_to_simplex(vals, simplex_point, n);
   
   // Log-likelihood
   for (int i = 0; i < nrow; i++) {
      double lh = 0.0;
      for (int j = 0; j < n + 1; j++)
         lh += lhm(i,j) * simplex_point[j];
      ltd += std::log(lh);
   } 
   
   // neighbor-diff penalty
   for (int i = 1; i < n; i++) {
      double ddiff = simplex_point[i+1] - 2*simplex_point[i] + simplex_point[i-1];
      ltd -= penalty_coef*nrow * ddiff*ddiff;
   } 
   return ltd;
} 

// [[Rcpp::export]]
Rcpp::NumericVector mixmod( Rcpp::NumericMatrix lhm, const double penalty_coef=0. ) {
   
   int n = lhm.ncol() - 1;
   double x0[n];
   double sigma0[n];
   double *const *pop;
   double *funvals;
   int iter = 0;
   cmaes_t evo;
   
   if( n <= 1 )
      Rcpp::stop( "matrix must have at least 2 columns" );
   
   Rcpp::Rcout << "FOO" << std::endl;
   std::fill_n( x0, n, 0. );
   std::fill_n( sigma0, n, 1. );
   funvals = cmaes_init( &evo, n, x0, sigma0, 0, 0, NULL );
   evo.sp.stopTolFun = 1e-2;
   evo.sp.flgsupplemented = true;
   
   while( !cmaes_TestForTermination(&evo) ) {
      pop = cmaes_SamplePopulation(&evo);
      
      for( int i = 0; i < evo.sp.lambda; ++i ) {
         funvals[i] = -log_target_dens( pop[i], lhm, penalty_coef );
      }
      cmaes_UpdateDistribution( &evo, funvals );  
      iter++;
      
      Rcpp::checkUserInterrupt();
   }
      
   double v[n];
   for( int i=0; i<n; i++ ) {
      v[i] = 1. / ( 1. + std::exp(-evo.rgxmean[i]) );
   }
   Rcpp::NumericVector ans(n+1);
   transform_to_simplex( v, ans.begin(), n );

   ans.attr("stopping_reason") = cmaes_TestForTermination(&evo);
   ans.attr("population_size") = evo.sp.lambda;
   ans.attr("iterations") = iter;
   
   cmaes_exit( &evo );
   return ans;
}

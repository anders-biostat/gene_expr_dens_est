#define MCMC_ENABLE_ARMA_WRAPPERS
#include "mcmc.hpp"
#include <Rcpp.h>

void transform_to_simplex( const arma::vec& vals, arma::vec& simplex_point )
// Note: simplex_point.n_elem must be vals.n_elem + 1
{
   const int len = vals.n_elem;
   double cumprod = 1.;
   for( int i = 0; i < len; i++ ) {
      simplex_point[i] = cumprod * ( 1 - vals[i] );
      cumprod *= vals[i];
   }
   simplex_point[len] = cumprod;
}

double log_target_dens( const arma::vec& vals_inp, void* ll_data )
{
   const int len = vals_inp.n_elem;
   double ltd = 0.;
   
   // Calculate simplex prior
   for( int i = 0; i < len; i++ ) {
      ltd += (len-i-1) * std::log(vals_inp[i]);
   }

   return ltd;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix run_sampler( )
{
  
   const int ndims = 6;
   arma::vec initial_vals(ndims-1);
   initial_vals.fill(0.5);
  
   mcmc::algo_settings_t settings;
  
   settings.rwmh_settings.n_burnin_draws = 1000;
   settings.rwmh_settings.n_keep_draws = 2000;
   
   settings.vals_bound = true;
   settings.lower_bounds = arma::vec(ndims-1, arma::fill::zeros);
   settings.upper_bounds = arma::vec(ndims-1, arma::fill::ones);

   arma::mat draws_out;
   
   mcmc::rwmh( initial_vals, log_target_dens, draws_out, NULL, settings );

   Rcpp::NumericMatrix simplex_draws_out( draws_out.n_rows, ndims );
   for( int i = 0; i < draws_out.n_rows; i++ ) {
      arma::vec simplex_point(ndims);
      transform_to_simplex( draws_out.row(i).t(), simplex_point );
      for( int j = 0; j < ndims; j++ )
         simplex_draws_out(i,j) = simplex_point[j];
   }
   return simplex_draws_out;
}



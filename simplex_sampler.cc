#define MCMC_ENABLE_ARMA_WRAPPERS
#include "mcmc.hpp"
#include <Rcpp.h>

double log_target_dens( const arma::vec& vals_inp, void* ll_data )
{
   const int dim = vals_inp.n_elem;
   double ltd = 0.;
   for( int i = 0; i < dim; i++ ) {
      const double beta_a = dim-i;
      const double beta_b = 1.0;
      ltd += std::log( std::pow( vals_inp[i], beta_a-1. ) * std::pow( 1. - vals_inp[i], beta_b-1. ) /
         ( std::tgamma( beta_a ) * std::tgamma( beta_b ) / std::tgamma( beta_a+beta_b ) ) );
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

   Rcpp::NumericMatrix draws_out_r( draws_out.n_rows, draws_out.n_cols );
   std::memcpy( draws_out_r.begin(), draws_out.memptr(), draws_out.n_elem * sizeof(double) );
   return draws_out_r;
}



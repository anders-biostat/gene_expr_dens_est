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

double log_target_dens_with_grad( const arma::vec& vals_inp, arma::vec* grad_out, void* data )
{
   const int len = vals_inp.n_elem;
   const arma::mat prob_nb = *reinterpret_cast<arma::mat*>(data);
  
   double ltd = 0.;
   
   // Calculate simplex prior
   for( int j = 0; j < len-1; j++ ) {
      ltd += (len-j-1) * std::log(vals_inp[j]);
   }

   // Get simplex point   
   arma::vec simplex_point(len+1);
   transform_to_simplex( vals_inp, simplex_point );
   
   arma::vec cell_probs = prob_nb * simplex_point;
   
   // Calculate log likelihood
   for( int i = 0; i < prob_nb.n_rows; i++ ) {
      ltd += std::log(cell_probs[i]);
   }

   if( grad_out ) {
      
      //arma::vec simplex_grad = (1./cell_probs).t() * prob_nb;
      arma::vec simplex_grad( len+1, arma::fill::zeros );
      for( int r = 0; r <= len; r++ ) {
         for( int i = 0; i < prob_nb.n_rows; i++ ) {
           simplex_grad[r] += prob_nb(i,r) / cell_probs[i];
         }
      }

      for( int j = 0; j < len; j++ ) {
         double a = (len-j-1) / vals_inp[j];
         a += simplex_grad[j] * simplex_point[j] / (vals_inp[j] - 1);
         for( int r = j+1; r <= len; r++ ) {
            a += simplex_grad[r] * simplex_point[r] / vals_inp[j];
         }
         (*grad_out)[j] = a;
      }

   }
   
   return ltd;
}

double log_target_dens( const arma::vec& vals_inp, void* data ) {
   return log_target_dens_with_grad( vals_inp, NULL, data );
}

// [[Rcpp::export]]
Rcpp::NumericMatrix run_sampler( Rcpp::NumericMatrix prob_nb, double stepsize = .03 )
{
  
   arma::mat prob_nb_arma( prob_nb.begin(), prob_nb.nrow(), prob_nb.ncol(), false, true );
  
   const int ndims = prob_nb.ncol();
   arma::vec initial_vals(ndims-1);
   initial_vals.fill(0.5);
  
   mcmc::algo_settings_t settings;
  
   settings.rwmh_settings.n_burnin_draws = 3000;
   settings.rwmh_settings.n_keep_draws = 20000;
   settings.rwmh_settings.par_scale = stepsize;
   
   settings.vals_bound = true;
   settings.lower_bounds = arma::vec(ndims-1, arma::fill::zeros);
   settings.upper_bounds = arma::vec(ndims-1, arma::fill::ones);

   arma::mat draws_out;
   
   mcmc::rwmh( initial_vals, log_target_dens, draws_out, &prob_nb_arma, settings );

   std::cout << "acceptance rate: " << static_cast<double>(settings.rwmh_settings.n_accept_draws) / 
      settings.rwmh_settings.n_keep_draws << std::endl;   
   
   Rcpp::NumericMatrix simplex_draws_out( draws_out.n_rows, ndims );
   for( int i = 0; i < draws_out.n_rows; i++ ) {
      arma::vec simplex_point(ndims);
      transform_to_simplex( draws_out.row(i).t(), simplex_point );
      for( int j = 0; j < ndims; j++ )
         simplex_draws_out(i,j) = simplex_point[j];
   }
   return simplex_draws_out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix run_sampler2( Rcpp::NumericMatrix prob_nb, Rcpp::NumericVector initial_vals )
{
   
   arma::mat prob_nb_arma( prob_nb.begin(), prob_nb.nrow(), prob_nb.ncol(), false, true );

   const int ndims = prob_nb.ncol();
   arma::vec initial_vals_arma( initial_vals.begin(), initial_vals.size() );

   mcmc::algo_settings_t settings;

   settings.mala_settings.n_burnin_draws = 1000;
   settings.mala_settings.n_keep_draws = 2000;
   
   settings.vals_bound = true;
   settings.lower_bounds = arma::vec(ndims-1, arma::fill::zeros);
   settings.upper_bounds = arma::vec(ndims-1, arma::fill::ones);
   
   arma::mat draws_out;
   
   mcmc::mala( initial_vals_arma, log_target_dens_with_grad, draws_out, &prob_nb_arma, settings );
   
   Rcpp::NumericMatrix simplex_draws_out( draws_out.n_rows, ndims );
   for( int i = 0; i < draws_out.n_rows; i++ ) {
      arma::vec simplex_point(ndims);
      transform_to_simplex( draws_out.row(i).t(), simplex_point );
      for( int j = 0; j < ndims; j++ )
         simplex_draws_out(i,j) = simplex_point[j];
   }
   return simplex_draws_out;
}



// [[Rcpp::export]]
Rcpp::List log_target_dens_and_grad( Rcpp::NumericVector vals, Rcpp::NumericMatrix prob_nb )
{
   arma::mat prob_nb_arma( prob_nb.begin(), prob_nb.nrow(), prob_nb.ncol(), true, true );
   arma::vec vals_arma( vals.begin(), vals.size() );
   arma::vec grad_out( vals.size() );
   double ltd = log_target_dens_with_grad( vals_arma, &grad_out, &prob_nb_arma );
   return Rcpp::List::create( 
     Rcpp::Named("log_density") = ltd,
     Rcpp::Named("gradient") = Rcpp::NumericVector( grad_out.begin(), grad_out.end() ) );
}


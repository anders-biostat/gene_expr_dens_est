#define MCMC_ENABLE_ARMA_WRAPPERS
#include "mcmc.hpp"
#include <Rcpp.h>

struct norm_data_t {
  arma::vec x;
};


double log_target_dens(const arma::vec& vals_inp, arma::vec* grad_out, void* ll_data)
{
  const double pi = arma::datum::pi;
  
  const double mu    = vals_inp(0);
  const double sigma = vals_inp(1);
  
  norm_data_t* dta = reinterpret_cast<norm_data_t*>(ll_data);
  const arma::vec x = dta->x;
  const int n_vals = x.n_rows;
  
  //
  
  const double ret = - n_vals * (0.5 * std::log(2*pi) + std::log(sigma)) - arma::accu( arma::pow(x - mu,2) / (2*sigma*sigma) );
  
  //
  
  if (grad_out) {
    grad_out->set_size(2,1);
    
    //
    
    const double m_1 = arma::accu(x - mu);
    const double m_2 = arma::accu( arma::pow(x - mu,2) );
    
    (*grad_out)(0,0) = m_1 / (sigma*sigma);
    (*grad_out)(1,0) = (m_2 / (sigma*sigma*sigma)) - ((double) n_vals) / sigma;
  }
  
  //
  
  return ret;
}


// [[Rcpp::export]]
Rcpp::NumericVector run_nuts_sampler( Rcpp::NumericVector data )
{
  norm_data_t dta;
  
  dta.x = Rcpp::as<std::vector<double>>(data);
  
  arma::vec initial_val(2);
  initial_val(0) = 3.0; // mu
  initial_val(1) = 1.5; // sigma
  
  mcmc::algo_settings_t settings;

  settings.nuts_settings.n_burnin_draws = 2000;
  settings.nuts_settings.n_keep_draws = 2000;
  
  arma::mat draws_out;
  mcmc::nuts(initial_val, log_target_dens, draws_out, &dta, settings);
  
  //
  
  std::cout << "nuts mean:\n" << arma::mean(draws_out) << std::endl;
  std::cout << "acceptance rate: " << static_cast<double>(settings.nuts_settings.n_accept_draws) / settings.nuts_settings.n_keep_draws << std::endl;
  
  //
  
  auto vecout = arma::conv_to<std::vector<double>>::from( arma::mean(draws_out) );
  return Rcpp::NumericVector( vecout.begin(), vecout.end() );
}



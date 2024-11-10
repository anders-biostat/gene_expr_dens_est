#include <armadillo>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]

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

double log_target_dens( const arma::vec& vals_inp, const Rcpp::NumericMatrix& comp_lh )
{
   const int len = vals_inp.n_elem;

   double ltd = 0.;
   
   // Do logistic transform
   arma::vec vals = 1. / ( 1. + arma::exp(-vals_inp) );

   // Add its derivative
   for( int j = 0; j < len-1; j++ ) {
      ltd += std::log( vals[j] * ( 1. - vals[j] ) );
   }
   
   // Calculate simplex prior
   for( int j = 0; j < len-1; j++ ) {
      ltd += (len-j-1) * std::log(vals[j]);
   }

   // Get simplex point   
   arma::vec simplex_point(len+1);
   transform_to_simplex( vals, simplex_point );

   // Calculate log likelihood
   for( int i = 0; i < comp_lh.nrow(); i++ ) {
      double lh = 0;
      for( int j = 0; j < len+1; j++ )
         lh += comp_lh(i,j) * simplex_point[j];
      ltd += std::log(lh);
   }

   return ltd;
}



// [[Rcpp::export]]
Rcpp::NumericMatrix sample_mixture_weights( Rcpp::NumericMatrix comp_lh, const double stepsize = .03,
      const int n_burnin_draws = 3000, const int n_keep_draws = 20000 )
{
   const int ndims = comp_lh.ncol();

   // matrix to take up result
   Rcpp::NumericMatrix draws( n_keep_draws, ndims );
   
   // initialize RNG
   std::random_device rd;  
   std::mt19937 gen(rd()); 
   std::uniform_real_distribution<> dis(0.0, 1.0);
   // fixme: the step is drawn from another RNG than the acceptance uniform
   
   arma::vec current_draw(ndims-1);

   // set initial vals
   current_draw.fill(.1);

   double current_log_prob = log_target_dens( current_draw, comp_lh );

   int n_accepts = 0;
   for( int i = 0; i < n_burnin_draws+n_keep_draws; i++ ) {
      
      // Get new proposal and its posterior prob
      arma::vec proposal = current_draw + arma::randn<arma::vec>(ndims-1) * stepsize;
      double proposal_log_prob = log_target_dens( proposal, comp_lh );
      
      // Do MCMC step
      double randval = dis(gen);
      if( randval < std::exp( proposal_log_prob - current_log_prob ) ) {
         current_draw = proposal;
         current_log_prob = proposal_log_prob;
         if( i >= n_burnin_draws )
            n_accepts++;
      }
      
      // Save draw, transformed back to simplex position
      if( i >= n_burnin_draws ) {
         arma::vec simplex_point(ndims);
         transform_to_simplex( 1. / ( 1. + arma::exp(-current_draw) ), simplex_point );
         
         for( int j = 0; j < ndims; j++ ) {
            draws( i-n_burnin_draws, j ) = simplex_point[j];
         }
      }
   }
   
   draws.attr("n_accepted_steps") = n_accepts;
   
   return draws;
}

#include <armadillo>
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

double log_target_dens( const arma::vec& vals_inp, void* data )
{
   const int len = vals_inp.n_elem;
   const arma::mat prob_nb = *reinterpret_cast<arma::mat*>(data);
  
   double ltd = 0.;
   
   // Do logistic transform
   const arma::vec vals = 1. / ( 1. + arma::exp(-vals_inp) );
   
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
   
   arma::vec cell_probs = prob_nb * simplex_point;
   
   // Calculate log likelihood
   for( int i = 0; i < prob_nb.n_rows; i++ ) {
      ltd += std::log(cell_probs[i]);
   }

   return ltd;
}



// [[Rcpp::export]]
Rcpp::NumericMatrix run_manual_sampler( Rcpp::NumericMatrix prob_nb, double stepsize = .03 )
{
   const int ndims = prob_nb.ncol();

   const int n_burnin_draws = 3000;
   const int n_keep_draws = 20000;
 
   arma::mat prob_nb_arma( prob_nb.begin(), prob_nb.nrow(), prob_nb.ncol(), false, true );
   
   Rcpp::NumericMatrix draws( n_keep_draws, ndims );
   
   std::random_device rd;  
   std::mt19937 gen(rd()); 
   std::uniform_real_distribution<> dis(0.0, 1.0);
   
   arma::vec current_draw(ndims-1);
   std::cout << ndims << std::endl;
   
   // set initial vals
   current_draw.fill(.1);
   
   double current_log_prob = log_target_dens( current_draw, &prob_nb_arma );

   for( int i = 0; i < n_burnin_draws+n_keep_draws; i++ ) {
      arma::vec proposal = current_draw + arma::randn<arma::vec>(ndims-1) * stepsize;
      double proposal_log_prob = log_target_dens( proposal, &prob_nb_arma );
      double randval = dis(gen);
      if( randval < std::exp( proposal_log_prob - current_log_prob ) ) {
         current_draw = proposal;
         current_log_prob = proposal_log_prob;
         std::cout << "A";
      } else {
         std::cout << "R";
      }
      if( i >= n_burnin_draws ) {
         arma::vec simplex_point(ndims);
         transform_to_simplex( 1. / ( 1. + arma::exp(-current_draw) ), simplex_point );
         for( int j = 0; j < ndims; j++ ) {
            draws( i-n_burnin_draws, j ) = simplex_point[j];
         }
      }
   }
   return draws;
}

#include <armadillo>
#include <Rcpp.h>
#include <thread>
#include <random>
#include <vector>

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

double log_target_dens(const arma::vec& vals_inp, const Rcpp::NumericMatrix& comp_lh) {
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
   for (int i = 0; i < comp_lh.nrow(); i++) {
      double lh = 0.0;
      for (int j = 0; j < len + 1; j++)
         lh += comp_lh(i,j) * simplex_point[j];
      ltd += std::log(lh);
   }
   
   return ltd;
}

Rcpp::NumericMatrix run_one_chain(const Rcpp::NumericMatrix& comp_lh, double stepsize,
                                  int n_burnin_draws, int n_keep_draws, unsigned int seed) {
   const int ndims = comp_lh.ncol();
   Rcpp::NumericMatrix draws(n_keep_draws, ndims);
   
   // RNG per chain
   std::mt19937 gen(seed);
   std::uniform_real_distribution<> uniform_dis(0.0, 1.0);
   std::normal_distribution<double> normal_dis(0.0, 1.0);
   
   arma::vec current_draw(ndims - 1, arma::fill::value(0.1));
   double current_log_prob = log_target_dens(current_draw, comp_lh);
   
   int n_accepts = 0;
   double target_accept_rate = 0.234;
   int adapt_interval = 100;
   int accepted_since_last_adapt = 0;
   
   for (int i = 0; i < n_burnin_draws + n_keep_draws; i++) {
      // Generate proposal
      arma::vec proposal(ndims - 1);
      for (int k = 0; k < ndims - 1; k++) {
         double z = normal_dis(gen);
         proposal[k] = current_draw[k] + z * stepsize;
      }
      
      double proposal_log_prob = log_target_dens(proposal, comp_lh);
      
      double randval = uniform_dis(gen);
      bool accepted = (randval < std::exp(proposal_log_prob - current_log_prob));
      
      if (accepted) {
         current_draw = proposal;
         current_log_prob = proposal_log_prob;
         if (i >= n_burnin_draws)
            n_accepts++;
      }
      
      // Adapt stepsize during burn-in
      if (i < n_burnin_draws) {
         if (accepted) accepted_since_last_adapt++;
         if ((i + 1) % adapt_interval == 0) {
            double current_accept_rate = (double)accepted_since_last_adapt / adapt_interval;
            stepsize *= std::exp((current_accept_rate - target_accept_rate));
            accepted_since_last_adapt = 0;
         }
      }
      
      if (i >= n_burnin_draws) {
         arma::vec simplex_point(ndims);
         transform_to_simplex(1.0 / (1.0 + arma::exp(-current_draw)), simplex_point);
         for (int j = 0; j < ndims; j++) {
            draws(i - n_burnin_draws, j) = simplex_point[j];
         }
      }
      
   }
   
   draws.attr("n_accepted_steps") = n_accepts;
   draws.attr("final_stepsize") = stepsize;
   return draws;
}

// [[Rcpp::export]]
Rcpp::List sample_mixture_weights_threads(Rcpp::NumericMatrix comp_lh, 
                                          double stepsize = 0.03,
                                          int n_burnin_draws = 3000, 
                                          int n_keep_draws = 20000,
                                          unsigned int global_seed = 12345) {
   const int n_chains = 4;
   std::vector<Rcpp::NumericMatrix> chain_results(n_chains);
   
   // Create a single global RNG
   std::mt19937 global_gen(global_seed);
   
   // Draw seeds for each chain from global_gen
   std::vector<unsigned int> chain_seeds(n_chains);
   for (int c = 0; c < n_chains; c++) {
      chain_seeds[c] = global_gen();
   }
   
   std::vector<std::thread> threads;
   threads.reserve(n_chains);
   
   for (int c = 0; c < n_chains; c++) {
      threads.emplace_back([&, c]() {
         Rcpp::NumericMatrix res = run_one_chain(comp_lh, stepsize, n_burnin_draws, n_keep_draws, chain_seeds[c]);
         chain_results[c] = res;
      });
   }
   
   for (auto& t : threads) {
      t.join();
   }
   
   Rcpp::List result(n_chains);
   for (int c = 0; c < n_chains; c++) {
      result[c] = chain_results[c];
   }
   
   return result;
}

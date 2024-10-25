#### Try the MCMC samples

## Compile the C++ code

library(Rcpp)

mcmc_basedir <- "/home/anders/w/gene_expr_dens_est/mcmc"

Sys.setenv( 
  PKG_CXXFLAGS = sprintf( "-I%s/include %s -O0 -UNDEBUG", mcmc_basedir, Sys.getenv("PKG_CXXFLAGS") ),
  PKG_LIBS = sprintf( "-L%s -lmcmc -Wl,-rpath,%s %s", mcmc_basedir, mcmc_basedir, Sys.getenv("PKG_LIBS") ) )

Rcpp::sourceCpp("simplex_sampler.cc", verbose=TRUE, rebuild=TRUE ); 


## Construct example data
# number of cells, read count totals per cell, true fractions, counts
n <- 1500
s <- round( 10^rnorm( n, log10(2500), .5 ) )
truelambda <- 10^ifelse( runif(n)<.7, rnorm( n, -3.3, .4 ), rnorm( n, -1.7, .2 ) )
k <- rpois( n, truelambda*s )

## Construct density basis
stepsize <- log(1.2)  # log(nu)
mu <- exp( seq( log(1e-5), 0, by=stepsize ) )
shape <- 2 / stepsize  # kappa
scale <- mu / shape    # theta

## Pre-calculate NB probs
pmf_nb <- sapply( mu, function(mu_) dnbinom( k, mu=s*mu_, size=shape ) )

# Run sampler
draws <- run_sampler(pmf_nb)

# Get means
spi <- colMeans(draws)

# Plot result
hist( log10(truelambda), freq=FALSE, xlab="log10(lambda)", 30, main="" )
xg <- seq( -6, 0, length.out=1000 )      
lines( xg,
       sapply( 1:length(mu), function(i)
         dgamma( 10^xg, shape=mu[i]/scale[i], scale=scale[i] )*(10^xg)*log(10) ) %*% spi,
       col="red" )

# Do this 20 times
hist( log10(truelambda), freq=FALSE, xlab="log10(lambda)", 30, main="" )
xg <- seq( -6, 0, length.out=1000 )      
for( i in 1:20 ) {
   draws <- run_sampler(pmf_nb)
   spi <- colMeans(draws)
   lines( xg,
       sapply( 1:length(mu), function(i)
         dgamma( 10^xg, shape=mu[i]/scale[i], scale=scale[i] )*(10^xg)*log(10) ) %*% spi,
       col="red" )
}

## Appendix: Show that NB and gamma fit
i <- 12
xg <- seq( -20, 1, length.out=1000 )
ss <- 1234
hist( rpois( 100000, ss*rgamma( 100000, shape=shape, scale=scale[i] ) ), 100, freq=FALSE )
kg <- 0:200
points( kg, dnbinom( kg, mu=ss*mu[i], size=shape ), col="red" )





Rcpp::sourceCpp("simplex_sampler.cc", verbose=TRUE, rebuild=TRUE ); 
x1 <- runif( ncol(pmf_nb)-1 )
r1 <- log_target_dens_and_grad( x1, pmf_nb )
x2 <- x1 
x2[7] <- x2[7] * 1.001
r2 <- log_target_dens_and_grad( x2, pmf_nb )
r1$gradient[7]
( r2$log_density - r1$log_density ) / (x2-x1)[7]

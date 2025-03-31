#### Testbed for the MCMC sampler

## Compile the C++ code

library(Rcpp)

Rcpp::sourceCpp("simplex_sampler_threaded.cc", verbose=TRUE, rebuild=TRUE ); 


## Construct example data
# number of cells, read count totals per cell, true fractions, counts
n <- 1500
s <- round( 10^rnorm( n, log10(2500), .5 ) )
truelambda <- 10^ifelse( runif(n)<.7, rnorm( n, -3.3, .4 ), rnorm( n, -1.7, .2 ) )
#truelambda <- rep( 1e-2, n )
k <- rpois( n, truelambda*s )

## Construct density basis
stepsize <- log(1.2)  # log(nu)
mu <- exp( seq( log(3e-5), 0, by=stepsize ) )
#mu <- sample(mu)
shape <- 10 / stepsize  # kappa
scale <- mu / shape    # theta

## Pre-calculate NB probs
pmf_nb <- sapply( mu, function(mu_) dnbinom( k, mu=s*mu_, size=shape ) )

# Run sampler
draws <- sample_mixture_weights_threads( pmf_nb, n_burnin_draws=0L, n_keep_draws=10000L, 
               temperature_decrease=.95, global_seed = as.integer( runif(1,0,100000) ),  )


plot(draws[[1]][,15], type="l", ylim=c(0,.04), xlab="simulated annealing time", ylab="coef #15")
lines(draws[[2]][,15], col=2)
lines(draws[[3]][,15], col=3)
lines(draws[[4]][,15], col=4)

# Get means
#spi <- sapply( draws, colMeans )
spi <- sapply( draws, function(x) x[ nrow(x), ] )
max( dist(t(spi)) )

# Plot result
hist( log10(truelambda), freq=FALSE, xlab="log10(lambda)", 30, main="", xlim=c(-6,1), ylim=c(0,.8) )
xg <- seq( -6, 0, length.out=1000 )      
for( i in 1:ncol(spi) )
   lines( xg,
       sapply( 1:length(mu), function(i)
         dgamma( 10^xg, shape=mu[i]/scale[i], scale=scale[i] )*(10^xg)*log(10) ) %*% spi[,i],
       col="red" )


###########










a <- t( sapply( 1:length(mu), function(i)
   dgamma( 10^xg, shape=mu[i]/scale[i], scale=scale[i] )*(10^xg)*log(10) ) %*% t(draws) )
lines( xg, colMeans(a) )
lines( xg, colMeans(a)+apply(a,2,sd) )
lines( xg, colMeans(a)-apply(a,2,sd) )

# Do this 20 times
hist( log10(truelambda), freq=FALSE, xlab="log10(lambda)", 30, main="", xlim=c(-6,1) )
xg <- seq( -6, 1, length.out=1000 )      
for( i in 1:10 ) {
   draws <- run_sampler(pmf_nb, .17)
   spi <- colMeans(draws)
   lines( xg,
       sapply( 1:length(mu), function(i)
         dgamma( 10^xg, shape=mu[i]/scale[i], scale=scale[i] )*(10^xg)*log(10) ) %*% spi,
       col="red" )
}
lines( xg,
       sapply( 1:length(mu), function(i)
          dgamma( 10^xg, shape=mu[i]/scale[i], scale=scale[i] )*(10^xg)*log(10) ) %*% rep(1/length(mu),length(mu)),
       col="lightblue" )

## Appendix: Show that NB and gamma fit
i <- 12
xg <- seq( -20, 1, length.out=1000 )
ss <- 1234
hist( rpois( 100000, ss*rgamma( 100000, shape=shape, scale=scale[i] ) ), 100, freq=FALSE )
kg <- 0:200
points( kg, dnbinom( kg, mu=ss*mu[i], size=shape ), col="red" )



## Appendix: Check the gradient
Rcpp::sourceCpp("simplex_sampler.cc", verbose=TRUE, rebuild=TRUE ); 
x1 <- runif( ncol(pmf_nb)-1 )
r1 <- log_target_dens_and_grad( x1, pmf_nb )
sapply( 1:length(x1), function(i) {
   x2 <- x1 
   x2[i] <- x2[i] * 1.001
   r2 <- log_target_dens_and_grad( x2, pmf_nb )
   ( r2$log_density - r1$log_density ) / (x2-x1)[i] } ) -> grad2
grad2
( r1$gradient - grad2 ) / grad2
plot( asinh(r1$gradient), asinh(grad2), asp=1 )


## Traces
matplot( draws, type="l", lty="solid" )
matplot( t(apply( draws, 1, cumsum )), type="l", lty="solid" )


t(sapply( 1:nrow(draws), function(d) {
   a <- round( 300*cumsum( draws[d,] ) )
   x <- rep( -1, 300 )
   for( i in 1:(length(a)-1) )
      x[ a[i] : a[i+1] ] <- i
   x } )) -> a
image(a)




##

spi <- runif( ncol(draws) )


spi <- exp( optim( log(colMeans(draws)), function(x) { x <- exp(x); x <- x/sum(x); -sum( log( pmf_nb %*% x ) ) } )$par )


spi <- exp( optim( rnorm(length(mu)), function(x) { x <- exp(x); x <- x/sum(x); -sum( log( pmf_nb %*% x ) ) } )$par )

spi <- spi/sum(spi)

## Construct test data
n <- 1500
s <- round( 10^rnorm( n, log10(2500), .5 ) )
truelambda <- 10^ifelse( runif(n)<.7, rnorm( n, -3.3, .4 ), rnorm( n, -1.7, .2 ) )
k <- rpois( n, truelambda*s )

hist( log10(truelambda), freq=FALSE, xlab="log10(lambda)", 30, main="", xlim=c(-6,1) )

## Estimator hyperparameters
min_mu <- 3e-5
stepsize <- 1.4
overlap_factor <- 1

## Construct density basis
mu <- exp( seq( log(min_mu), 0, by=log(stepsize) ) )
shape <- 30/overlap_factor / stepsize  # kappa
scale <- mu / shape    # theta
smoothness <- 2

## Calculate NB likelihoods
nbl <- sapply( mu, function(mu_) dnbinom( k, mu=s*mu_, size=shape ) )

## Find likelihood maximum
map_to_simplex <- function(x) {
   v <- 1/(1+exp(-x))
   spi <- c( 1, cumprod(v) ) * c( 1-v, 1 ) }

objective <- function(x) {
   spi <- map_to_simplex( x )
   -( sum( log( nbl %*% spi ) ) + 
         sum( log( v*(1-v) ) ) +
         sum( log(v) * (length(v) - 0:(length(v)-1) - 1) ) -
         sum( diff(spi)^2 ) * length(nbl)*smoothness/100 ) }

o <- optim( rnorm( ncol(nbl)-1 ), objective,
   method="SANN", control=list( maxit=10000, tmax=20) )

o <- optim( rnorm( ncol(nbl)-1 ), objective,
   method="Nelder-Mead", control=list(maxit=10000) )

## Get density estimate
spi <- map_to_simplex( o$par )
dens_est <- data.frame(
   log10lambda = seq( -6, 0, length.out=1000 ) )
dens_est$density <- (
   sapply( 1:length(mu), function(i)
      dgamma( 10^dens_est$log10lambda, shape=mu[i]/scale[i], scale=scale[i] ) * 
         (10^dens_est$log10lambda)*log(10) ) 
   %*% spi )

## Plot result
#hist( log10(truelambda), freq=FALSE, xlab="log10(lambda)", 30, main="", xlim=c(-6,1) )
xg <- seq( -6, 0, length.out=1000 )      
lines( dens_est, col="red" )


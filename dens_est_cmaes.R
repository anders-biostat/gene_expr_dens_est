## Construct test data
n <- 1500
s <- round( 10^rnorm( n, log10(2500), .5 ) )
truelambda <- 10^ifelse( runif(n)<.7, rnorm( n, -3.3, .4 ), rnorm( n, -1.7, .2 ) )
k <- rpois( n, truelambda*s )

hist( log10(truelambda), freq=FALSE, xlab="log10(lambda)", 30, main="", xlim=c(-6,1) )

## Compile C++ cpde
eigen_flags <- system("pkg-config --cflags eigen3", intern = TRUE)
custom_inc  <- normalizePath("libcmaes/include")
Sys.setenv(PKG_CXXFLAGS = paste(eigen_flags, "-I", custom_inc, "-lcmaes"))
lib_path <- normalizePath("libcmaes/lib")
Sys.setenv(PKG_LIBS = paste("-L", lib_path, "-lcmaes", "-Wl,-rpath", lib_path))

Rcpp::sourceCpp("test_cmaes.cc", verbose=TRUE, rebuild=TRUE)
test() 


## Estimator hyperparameters
min_mu <- 3e-5
stepsize <- 1.13
overlap_factor <- 1

## Construct density basis
mu <- exp( seq( log(min_mu), 0, by=log(stepsize) ) )
shape <- 30/overlap_factor / stepsize  # kappa
scale <- mu / shape    # theta
smoothness <- 2

## Calculate NB likelihoods
nbl <- sapply( mu, function(mu_) dnbinom( k, mu=s*mu_, size=shape ) )

x <- do_sann( nbl )

v <- 1/(1+exp(-x))
spi <- c( 1, cumprod(v) ) * c( 1-v, 1 ) 
dens_est <- data.frame(
   log10lambda = seq( -6, 0, length.out=1000 ) )
dens_est$density <- (
   sapply( 1:length(mu), function(i)
      dgamma( 10^dens_est$log10lambda, shape=mu[i]/scale[i], scale=scale[i] ) * 
         (10^dens_est$log10lambda)*log(10) ) 
   %*% spi )

xg <- seq( -6, 0, length.out=1000 )      
lines( dens_est, col="red" )


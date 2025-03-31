est_dens <- function( k, s, min_mu=3e-5, step_factor=1.2, shape_factor=10 ) {

	stepsize <- log(step_factor)
	mu <- exp( seq( log(min_mu), 0, by=stepsize ) )
	shape <- shape_factor / stepsize  # kappa
	scale <- mu / shape    # theta
	lh_nb <- sapply( mu, function(mu_) dnbinom( k, mu=s*mu_, size=shape ) )

	weights <- mixmod(lh_nb)

	dens_est <- data.frame(
	  log10lambda = seq( log10(min_mu)-1, 0, length.out=1000 ) )
	dens_est$density <- (
	  sapply( 1:length(mu), function(i)
	     dgamma( 10^dens_est$log10lambda, shape=mu[i]/scale[i], scale=scale[i] ) * 
	        (10^dens_est$log10lambda)*log(10) ) 
	  %*% weights )

	dens_est

}
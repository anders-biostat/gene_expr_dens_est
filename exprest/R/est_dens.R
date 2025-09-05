est_dens <- function( k, s, subsample=700, min_mu=3e-5, step_factor=1.2, shape_factor=10, pty=.3 ) {

	stopifnot( length(k) == length(s) )
	stopifnot( k == round(k) )
	stopifnot( s == round(s) )
	
	if( length(s) > subsample ) {
		 ss <- sample( 1:length(s), subsample )
		 k <- k[ss]
		 s <- s[ss]
	}

	stepsize <- log(step_factor)
	mu <- exp( seq( log(min_mu), 0, by=stepsize ) )
	shape <- shape_factor / stepsize  # kappa
	scale <- mu / shape    # theta
	lh_nb <- sapply( mu, function(mu_) dnbinom( k, mu=s*mu_, size=shape ) )

	weights <- mixmod( lh_nb, pty )

	dens_est <- data.frame(
	  log10lambda = seq( log10(min_mu)-1, 0, length.out=1000 ) )
	dens_est$density <- drop(
	  sapply( 1:length(mu), function(i)
	     dgamma( 10^dens_est$log10lambda, shape=mu[i]/scale[i], scale=scale[i] ) * 
	        (10^dens_est$log10lambda)*log(10) ) 
	  %*% weights )

	dens_est

}

est_dens_for_groups <- function( k, s, group, progressbar=TRUE, subsample=700, ... ) {
	if( progressbar )
  	pb <- txtProgressBar(max=length(s), style=3)
	a <- tapply( 1:length(s), group, function(ii) {
  	if( length(ii)>1 )
    	a <- exprest::est_dens( k[ii], s[ii] ) 
  	else
    	a <- data.frame( log10lambda=numeric(0), density=numeric(0) )
  	setTxtProgressBar( pb, getTxtProgressBar(pb) + min( subsample, length(ii) ) ) 
  a }, simplify=FALSE )
  a
}

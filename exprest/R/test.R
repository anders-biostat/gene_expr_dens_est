make_test_data <- function( n=1500 ) {
	s <- round( 10^rnorm( n, log10(2500), .5 ) )
	true_lambda <- 10^ifelse( runif(n)<.7, rnorm( n, -3.3, .4 ), rnorm( n, -1.7, .2 ) )
	k <- rpois( n, true_lambda*s )
	data.frame( k, s, true_lambda )
}

test <- function() {
	data <- make_test_data()
	dens <- est_dens( data$k, data$s )
	plot( density ~ log10lambda, dens, type="l" )
}
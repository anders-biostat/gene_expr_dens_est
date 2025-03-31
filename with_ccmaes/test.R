setwd("~/w/gene_expr_dens_est/with_ccmaes")
Rcpp::sourceCpp("mixmod.cc")
#mixmod(matrix())


n <- 1500
s <- round( 10^rnorm( n, log10(2500), .5 ) )
truelambda <- 10^ifelse( runif(n)<.7, rnorm( n, -3.3, .4 ), rnorm( n, -1.7, .2 ) )
k <- rpois( n, truelambda*s )
stepsize <- log(1.2)  # log(nu)
mu <- exp( seq( log(3e-5), 0, by=stepsize ) )
shape <- 10 / stepsize  # kappa
scale <- mu / shape    # theta
pmf_nb <- sapply( mu, function(mu_) dnbinom( k, mu=s*mu_, size=shape ) )

a <- mixmod(pmf_nb)
a
plot(a)

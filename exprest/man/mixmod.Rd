\name{mixmod}
\alias{mixmod}
\title{
Maximum Likelihood for Mixture Models
}
\description{
Performs maximization of the likelihood of a mixture model with 
fixed mixture components using CMA-ES
}
\usage{
mixmod(lhm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lhm}{the matrix of component likelihoods -- each 
    row is an observation, each column a component
}
}
\details{
Assume you have n observations x_i drawn i.i.d. from a mixture model
with probability density function $f = \sum_{j=1}^m w_j f_j(x)$ 
with weights $\sum_j w_j = 1$. You know the $f_j$ but not the $w_j$.
Then, you can compute the "component likelihood matrix" \code{lhm}
(or $L$) with $L_ij = f_j(x_i)$. 

This function takes this matrix and maximizes the likelihood $f$ as
function of the $w_j$ using the CMA-ES optimization algorithm.

Internally, the function uses the code from the c-cmaes reference
implementation by Nikolaus Hansen.
}
\value{
a vector with the maximizing weights $w_i$ -- with attributes 
\code{stopping_reason} with the final status reported by c-cmaes,
\code{iterations} with the number of ES iterations, and 
\code{population_size} with the size of the ES population.
}
\author{
Simon Anders (re-using code by N. Hansen et al.)
}

\seealso{
\code{\link{est_dens}}
}

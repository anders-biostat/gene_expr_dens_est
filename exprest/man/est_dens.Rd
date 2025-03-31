\name{est_dens}
\alias{est_dens}
\title{
Estimate Density underlying Poisson-widened count observations
}
\description{
In RNA-Seq analysis, one often uses the following model: The expression strength $\lambda$
of a gene of interest in each cell (or sample, in case of bulk data) is drawn
from an unknown distribution, which specifies the fraction of the cell's total
mRNA complement that orginates from this gene. The observed read or UMI count $k$
are then modeled as drawn from a Poisson distribution with expectation $s\lambda$,
where $s$ is the total number of reads (summed over all genes) for the cell.

This function estimates the distribution of the expression strengths $\lambda$
by performing maximum likelihood estimation over a space of trial functions
which are linear combinations of Gamma distribution densities with increasing means.
}
\usage{
est_dens(k, s, min_mu = 3e-05, step_factor = 1.2, shape_factor = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k}{a vector of counts for the gene of interest, one value per cell}
  \item{s}{a vector of total counts (summed over all genes), one value per cell}
  \item{min_mu}{the mean of the first (smallest) Gamma basis distribution. This 
  value should be close to the reciprocal of the typical value of $s$}
  \item{step_factor}{The factor by which the mean increases from one Gamma basis 
  distribution to the next. Bigger values result in coarser estimates.}
  \item{shape_factor}{The shape parameter of the basis Gamma distributions
  is log(step_factor) divided by this argument.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a data frame with two columns: 
\item{log10lamba}{A sequence of log10 lambda values, suitable as x values for 
plotting the density on a log scale}
\item{density}{the estimated probability density for the logarithmic expression 
strength distribution, suitable as y values for plotting the density}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Simon Anders
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.

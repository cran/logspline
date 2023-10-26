\name{dlogspline}
\alias{dlogspline}
\alias{plogspline}
\alias{qlogspline}
\alias{rlogspline}
\title{Logspline Density Estimation}
\description{Density (\code{dlogspline}), cumulative probability (\code{plogspline}), quantiles
(\code{qlogspline}), and  random samples (\code{rlogspline}) from 
a logspline density that was fitted using
the 1997 knot addition and deletion algorithm (\code{\link{logspline}}). 
The 1992 algorithm is available using the \code{\link{oldlogspline}} function.
}
\usage{dlogspline(q, fit, log = FALSE) 
plogspline(q, fit) 
qlogspline(p, fit) 
rlogspline(n, fit) }
\arguments{
\item{q}{
vector of quantiles. Missing values (NAs) are allowed. 
}
\item{p}{
vector of probabilities. Missing values (NAs) are allowed. 
}
\item{n}{
sample size. If \code{length(n)} is larger than 1, then 
\code{length(n)} random values are returned.
}
\item{fit}{
\code{logspline} object, typically the result of \code{\link{logspline}}.
}
\item{log}{
should dlogspline return densities (TRUE) or log-densities (FALSE)
}
}
\value{Densities (\code{dlogspline}), probabilities (\code{plogspline}), quantiles (\code{qlogspline}), 
or a random sample (\code{rlogspline}) from a \code{logspline} density that was fitted using
knot addition and deletion.
}
\details{
Elements of \code{q} or \code{p} that are missing will cause the 
corresponding elements of the result to be missing. 
}
\references{Charles Kooperberg and Charles J. Stone.  Logspline density estimation
for censored data (1992). \emph{Journal of Computational and Graphical
Statistics}, \bold{1}, 301--328.

Charles J. Stone, Mark Hansen, Charles Kooperberg, and Young K. Truong. 
The use of polynomial splines and their tensor products in extended
linear modeling (with discussion) (1997).  \emph{Annals of Statistics},
\bold{25}, 1371--1470.}
\author{ Charles Kooperberg \email{clk@fredhutch.org}.}
\seealso{ 
\code{\link{logspline}},      
\code{\link{plot.logspline}},
\code{\link{summary.logspline}},
\code{\link{oldlogspline}}.}
\examples{
x <- rnorm(100)
fit <- logspline(x)
qq <- qlogspline((1:99)/100, fit)
plot(qnorm((1:99)/100), qq)                  # qq plot of the fitted density
pp <- plogspline((-250:250)/100, fit)
plot((-250:250)/100, pp, type = "l")
lines((-250:250)/100,pnorm((-250:250)/100))  # asses the fit of the distribution
dd <- dlogspline((-250:250)/100, fit)
plot((-250:250)/100, dd, type = "l")
lines((-250:250)/100, dnorm((-250:250)/100)) # asses the fit of the density
rr <- rlogspline(100, fit)                   # random sample from fit
}
\keyword{distribution}
\keyword{smooth}

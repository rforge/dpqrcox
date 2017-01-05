\name{qcox}
\alias{qcox}
\title{
Quantile function for Coxian distributions
}
\description{
Quantile  function for the Coxian distribution at probability \code{p}, with index vector equal to \code{probvec} and rate vector equal to \code{lambdvec}.
}
\usage{
qcox(p, probvec, lambdvec, lower.tail = TRUE, log.p = FALSE)
}
\arguments{
  \item{p}{
vector of probabilities.
}
  \item{probvec}{
index vector 
}
  \item{lambdvec}{
rate vector 
}
  \item{lower.tail}{
logical; if TRUE (default), probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}
}
  \item{log.p}{
logical; if TRUE, probabilities \code{p} are given as \code{log(p)}.
}
}
\details{
Assume that \eqn{X_1}, ..., \eqn{X_n} be independent exponential random variables with respective rates \eqn{\lambda_1}, ..., \eqn{\lambda_n}. Assume that \eqn{P(I=i)=p_i}, \eqn{i=1}, ..., \eqn{n} and the  be  random variable \eqn{I} is independent of \eqn{X_1}, ..., \eqn{X_n}. Then the random variable
\deqn{X=X_1+\ldots+X_I} is a Coxian random variable with index vector \code{probvec}=\eqn{(p_1,\ldots,p_n)} and rate vector \code{lambdvec}=\eqn{(\lambda_1,\ldots,\lambda_n)}.

The probability distribution of a Coxian random variable is a mixture of convolutions of exponential distributions. A Coxian distribution can approximate any positive-valued distribution. However, approximations is poor sometimes. [1] indicated that numerical error in computing density and probability functions is one of reasons for the issue. \code{qcox} provides with quantile function of Coxian distributions by using the method of [1].}
\value{
\code{qcox} gives the quantile function.

The length of the result is the length of \code{p}.

Only the first elements of the logical arguments are used.
}
\references{
Qihong Duan and Junrong Liu. (2016), Modelling a Bathtub-Shaped Failure Rate by a Coxian Distribution, IEEE TRANSACTIONS ON RELIABILITY  65(2), pp.878-885.  DOI: 10.1109/TR.2015.2494374 [1]
}
\author{
Junrong Liu
}
\seealso{
\code{Distributions} for other standard distribution
}
\examples{
}
\keyword{ distribution}

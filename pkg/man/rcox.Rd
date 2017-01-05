\name{rcox}
\alias{rcox}
\title{
random generation for Coxian distributions
}
\description{
generates \code{n} random deviates from the Coxian distribution with index vector equal to \code{probvec} and rate vector equal to \code{lambdvec}.
}
\usage{
rcox(n, probvec, lambdvec)
}
\arguments{
  \item{n}{
number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
}
  \item{probvec}{
index vector 
}
  \item{lambdvec}{
rate vector 
}
}
\details{
Assume that \eqn{X_1}, ..., \eqn{X_n} be independent exponential random variables with respective rates \eqn{\lambda_1}, ..., \eqn{\lambda_n}. Assume that \eqn{P(I=i)=p_i}, \eqn{i=1}, ..., \eqn{n} and the  be  random variable \eqn{I} is independent of \eqn{X_1}, ..., \eqn{X_n}. Then the random variable
\deqn{X=X_1+\ldots+X_I} is a Coxian random variable with index vector \code{probvec}=\eqn{(p_1,\ldots,p_n)} and rate vector \code{lambdvec}=\eqn{(\lambda_1,\ldots,\lambda_n)}.

The probability distribution of a Coxian random variable is a mixture of convolutions of exponential distributions. A Coxian distribution can approximate any positive-valued distribution. However, approximations is poor sometimes. [1] indicated that numerical error in computing density and probability functions is one of reasons for the issue. \code{pcox} provides with distribution function of Coxian distributions by using the method of [1].
}
\value{
\code{rcox} generates \code{n} random deviates from the Coxian distribution with index vector equal to \code{probvec} and rate vector equal to \code{lambdvec}. The length of the result is determined by \code{n} for \code{rcox}.
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


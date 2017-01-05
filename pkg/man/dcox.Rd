\name{dcox}
\alias{dcox}
\title{Density function for Coxian distributions}
\description{
Density function for the Coxian distribution with index vector equal to \code{probvec} and rate vector equal to \code{lambdvec}.
}
\usage{
dcox(x, probvec, lambdvec, log = FALSE)
}
\arguments{
  \item{x}{
vector of positions.
}
  \item{probvec}{
index vector 
}
  \item{lambdvec}{
rate vector 
}
  \item{log}{
logical; if TRUE, density d are given as log(d)
}
}
\details{
Assume that \eqn{X_1}, ..., \eqn{X_n} be independent exponential random variables with respective rates \eqn{\lambda_1}, ..., \eqn{\lambda_n}. Assume that \eqn{P(I=i)=p_i}, \eqn{i=1}, ..., \eqn{n} and the  be  random variable \eqn{I} is independent of \eqn{X_1}, ..., \eqn{X_n}. Then the random variable
\deqn{X=X_1+\ldots+X_I} is a Coxian random variable with index vector \code{probvec}=\eqn{(p_1,\ldots,p_n)} and rate vector \code{lambdvec}=\eqn{(\lambda_1,\ldots,\lambda_n)}.

The probability distribution of a Coxian random variable is a mixture of convolutions of exponential distributions. A Coxian distribution can approximate any positive-valued distribution. However, approximations is poor sometimes. [1] indicated that numerical error in computing density and probability functions is one of reasons for the issue. \code{dcox} provides with density function of Coxian distributions by using the method of [1].
}
\value{
\code{dcox} gives the density function.

The length of the result is the length of \code{x}.

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
pv<-c(0.2,0.3,0.1,0.4)
lv<-c(2,1,1,2)
x<-seq(0.1,0.9,0.1)
dcox(x,pv,lv) # density at vector x

##
require(graphics)
pv<-c(0.4,0.3,0.4)
lv<-c(0.1,0.2,0.3)
plot(function(y) dcox(y, pv,lv ), 2, 20, main = "Coxian density")
}
\keyword{ distribution}

\name{QR.ip}
\alias{QR.ip}
\title{Quantile Regression (QR) use Interior Point (ip) Method}
\description{
   The function use the interior point method from Rmosek to solve the quantile regression
}
\usage{
QR.ip(X,y,tau)
}
\arguments{
  \item{X}{the design matrix}
  \item{y}{response variable}
  \item{tau}{quantile level}
 } 
\value{
  a \code{\link{list}} structure is with components
  \item{beta}{the vector of estimated coefficient}
  \item{b}{intercept}
}
\note{
  Need to install Rmosek package from CRAN. There are series of steps need to do before installing the Rmosek package.
}
\references{
 Koenker, Roger. \emph{Quantile Regression}, New York, 2005. Print.

 Kevin A.Mc, Shane Clyde L.Monma and David Shanno.(1989) An Implementation of a Primal-Dual Interior Point Method for Linear Programming. \emph{Journal on Computing, ISSN} 1091--9856

}
\examples{
set.seed(1)
n=100
p=2
a=rnorm(n*p, mean = 1, sd =1)
x=matrix(a,n,p)
beta=rnorm(p,1,1)
beta=matrix(beta,p,1)
y=x\%*\%beta-matrix(rnorm(n,0.1,1),n,1)
# x is 1000*10 matrix, y is 1000*1 vector, beta is 10*1 vector
QR.ip(x,y,0.1)
}

\keyword{IP}

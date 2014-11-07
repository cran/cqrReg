\name{QR.lasso.ip}
\alias{QR.lasso.ip}
\title{Quantile Regression (QR) with Adaptive Lasso Penalty (lasso) use Interior Point (ip) Method}
\description{
    The adaptive lasso parameter base on the estimated coefficient without penalty function.
   The function use the interior point method from Rmosek to solve the quantile regression with adaptive lasso penalty.
}
\usage{
QR.lasso.ip(X,y,tau,lambda)
}
\arguments{
  \item{X}{the design matrix}
  \item{y}{response variable}
  \item{tau}{quantile level}
  \item{lambda}{The constant coefficient of penalty function. (default lambda=1)}
 } 
\value{
  a \code{\link{list}} structure is with components
  \item{beta}{the vector of estimated coefficient}
  \item{b}{intercept}
  \item{lambda}{The constant coefficient of penalty function. (default lambda=1)}
}
\note{
  Need to install Rmosek package from CRAN. There are series of steps need to do before installing the Rmosek package.
}
\references{
 Kevin A.Mc, Shane Clyde L.Monma and David Shanno.(1989) An Implementation of a Primal-Dual Interior Point Method for Linear Programming. \emph{Journal on Computing, ISSN} 1091--9856

Wu, Yichao and Liu, Yufeng (2009). Variable selection in quantile regression. \emph{Statistica Sinica}, \bold{19}, 801--817.
}
\examples{
set.seed(1)
n=100
p=2
a=2*rnorm(n*2*p, mean = 1, sd =1)
x=matrix(a,n,2*p)
beta=2*rnorm(p,1,1)
beta=rbind(matrix(beta,p,1),matrix(0,p,1))
y=x\%*\%beta-matrix(rnorm(n,0.1,1),n,1)
# x is 1000*20 matrix, y is 1000*1 vector, beta is 20*1 vector with last ten zero value elements. 
QR.lasso.ip(x,y,0.1)
}

\keyword{IP}

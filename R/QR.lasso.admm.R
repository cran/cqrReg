QR.lasso.admm = function(X,y,tau,lambda,rho,beta,maxit)
{

if(nargs()<7){
	maxit = 200
}

n=dim(X)[1]
x=cbind(rep(1,n),X)

if(nargs()<6){
	beta=solve(t(x)%*%x,t(x)%*%y)
}

if(nargs()<5){
	rho=0.4
}

if(nargs()<4){
	lambda=1
}


betah=QRPADMMCPP(X,y,beta,maxit,tau,rho,lambda)
beta=matrix(betah[-1])
b=betah[1]
return(list(beta=beta,b=b))
}

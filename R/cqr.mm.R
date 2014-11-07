cqr.mm = function(X, y, tau, beta, maxit, toler)
{

if(nargs()<5){
	toler = 1e-3
	maxit = 200
}

n=dim(X)[1]
k=length(tau)

if(nargs()<4){
	beta=solve(t(X)%*%X,t(X)%*%y)
}

beta=CQRMMCPP(X,y,beta,toler,maxit,tau)
b=quantile(y-X%*%beta,tau)

return(list(beta=beta,b=b))
}

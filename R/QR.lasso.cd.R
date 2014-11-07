QR.lasso.cd = function(X,y,tau,lambda,beta,maxit,toler)
{

if(nargs()<6){
	toler = 1e-3
	maxit = 200
}

if(nargs()<5){
	beta=solve(t(X)%*%X,t(X)%*%y)
}

if(nargs()<4){
	lambda=1
}

beta0=QRCDCPP(X,y,beta,toler,maxit,tau)
betah=QRPCDCPP(X,y,beta,beta0,toler,maxit,tau,lambda)
beta=betah
b=quantile(y-X%*%beta,tau)
return(list(beta=beta,b=b))
}

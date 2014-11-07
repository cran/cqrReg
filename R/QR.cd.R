QR.cd = function(X, y, tau, beta, maxit, toler)
{

if(nargs()<5){
	toler = 1e-3
	maxit = 200
}


if(nargs()<4){
	beta=solve(t(X)%*%X,t(X)%*%y)
}

betah=QRCDCPP(X,y,beta,toler,maxit,tau)
beta=betah
b=quantile(y-X%*%beta,tau)
return(list(beta=beta,b=b))
}

QR.mm = function(X, y, tau, beta, maxit, toler)
{

if(nargs()<5){
	toler = 1e-3
	maxit = 200
}

n=dim(X)[1]
x=cbind(rep(1,n),X)

if(nargs()<4){
	beta=solve(t(x)%*%x,t(x)%*%y)
}

betah=QRMMCPP(X,y,beta,toler,maxit,tau)
beta=matrix(betah[-1])
b=betah[1]
return(list(beta=beta,b=b))
}

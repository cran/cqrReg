QR.admm = function(X, y, tau,rho,beta, maxit, toler)
{

if(nargs()<7){
	toler = 1e-2
	maxit = 200
}

n=dim(X)[1]
x=cbind(rep(1,n),X)

if(nargs()<5){
	beta=solve(t(x)%*%x,t(x)%*%y)
}

if(nargs()<4){
	rho=0.4
}

betah=QRADMMCPP(X,y,beta,toler,maxit,tau,rho)
beta=matrix(betah[-1])
b=betah[1]
return(list(beta=beta,b=b))
}

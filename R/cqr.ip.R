cqr.ip=function(X,y,tau)
{

# variable define
p = dim(X)[2] 
n = dim(X)[1]
k=length(tau)

	
	bmatrix=matrix(c(rep(1,n),rep (0,n*k)),n*k,k+1)[,-(k+1)]
	X= cbind(bmatrix,t(matrix(t (X),p,n*k)))
	y=rep(y,k)
	onematrix=cbind(diag(1,n*k),-diag(1,n*k))
        tauv=c(t(matrix(tau,k,n)))
	tauv=c(tauv,1-tauv)
        subject=cbind(X,onematrix)


lol=list()
lol$sense="min"
lol$c=c(rep(0,k+p),tauv)
lol$A=Matrix(c(subject),nrow=n*k,ncol=(p+k+2*n*k),sparse=TRUE)
lol$bc=t(cbind(y,y))
blx=c(rep(-Inf,p+k),rep(0,2*n*k))
bux=rep(Inf,p+k+2*n*k)
lol$bx=rbind(blx,bux);
betai=mosek(lol)
beta=betai$sol$itr$xx[(k+1):(p+k)]
b=betai$sol$itr$xx[1:k]
return(list(beta=beta,b=b))

}


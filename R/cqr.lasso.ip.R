cqr.lasso.ip=function(X,y,tau,lambda)
{

if(missing(lambda)){
	lambda=1
}

# variable define
p = dim(X)[2] 
n = dim(X)[1]
k=length(tau)
betaold=cqr.ip(X,y,tau)$beta


	
	bmatrix=matrix(c(rep(1,n),rep (0,n*k)),n*k,k+1)[,-(k+1)]
	X= cbind(bmatrix,t(matrix(t (X),p,n*k)))
	y=rep(y,k)
	onematrix=cbind(diag(1,n*k),-diag(1,n*k))
        tauv=c(t(matrix(tau,k,n)))
	tauv=c(tauv,1-tauv)
        subject=cbind(X,onematrix)

lol=list()
lol$sense="min"
lol$c=c(rep(0,k+p),tauv,c(lambda/betaold^2))
a=cbind(subject,matrix(0,n*k,p))
b=cbind(matrix(0,p,k),diag(1,p),matrix(0,p,2*n*k),-diag(1,p))
d=cbind(matrix(0,p,k),-diag(1,p),matrix(0,p,2*n*k),-diag(1,p))
A=rbind(a,b,d)
lol$A=Matrix(c(A),nrow=(n*k+2*p),ncol=(2*p+k+2*n*k),sparse=TRUE)
blc=c(y,rep(-Inf,(2*p)))
buc=c(y,rep(0,(2*p)))
lol$bc=rbind(blc,buc)
blx=c(rep(-Inf,p+k),rep(0,2*n*k+p))
bux=rep(Inf,2*(n*k+p)+k)
lol$bx=rbind(blx,bux);
betai=mosek(lol)
beta=betai$sol$itr$xx[(k+1):(p+k)]
b=betai$sol$itr$xx[1:k]
return(list(beta=beta,b=b))
}

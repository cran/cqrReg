QR.lasso.ip=function(X,y,tau,lambda)
{

if(missing(lambda)){
	lambda=1
}

# variable define
p = dim(X)[2] 
n = dim(X)[1]

betaold=QR.ip(X,y,tau)$beta

X= cbind(rep(1,n),X) 
  
lol=list()
lol$sense="min"
lol$c=c(rep(0,p+1),rep(tau,n),rep(1-tau,n),c(lambda/betaold^2))
a=cbind(X,diag(1,n),-diag(1,n),matrix(0,n,p))
b=cbind(matrix(0,p,1),diag(1,p),matrix(0,p,2*n),-diag(1,p))
d=cbind(matrix(0,p,1),-diag(1,p),matrix(0,p,2*n),-diag(1,p))
A=rbind(a,b,d)
lol$A=Matrix(c(A),nrow=(n+2*p),ncol=(2*p+1+2*n),sparse=TRUE)
blc=c(y,rep(-Inf,(2*p)))
buc=c(y,rep(0,(2*p)))
lol$bc=rbind(blc,buc)
blx=c(rep(-Inf,p+1),rep(0,2*n+p))
bux=rep(Inf,2*(n+p)+1)
lol$bx=rbind(blx,bux);
betai=mosek(lol)
beta=betai$sol$itr$xx[2:(p+1)]
b=betai$sol$itr$xx[1]
return(list(beta=beta,b=b))
}

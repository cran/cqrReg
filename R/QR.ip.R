QR.ip=function(X,y,tau)
{

# variable define
p = dim(X)[2] 
n = dim(X)[1]

X= cbind(rep(1,n),X) 
lol=list()
lol$sense="min"
lol$c=c(rep(0,p+1),rep(tau,n),rep(1-tau,n))
lol$A=Matrix(c(X,diag(1,n),-diag(1,n)),nrow=n,ncol=(p+1+2*n),sparse=T)
lol$bc=t(cbind(y,y))
blx=c(rep(-Inf,p+1),rep(0,2*n))
bux=rep(Inf,p+1+2*n)
lol$bx=rbind(blx,bux);
betai=mosek(lol)
beta=betai$sol$itr$xx[2:(p+1)]
b=betai$sol$itr$xx[1]
return(list(beta=beta,b=b))

}

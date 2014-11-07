#include "RcppArmadillo.h"
// [[Rcpp::export]]

arma::vec QRPCDCPP (arma::mat xr,arma::vec yr,arma::vec betar,arma::vec betaoldr,
double to,int m,double ta,double lamdar) {

double toler = (to);
int maxit=(m);
double tau=(ta), lamda=(lamdar);
arma:: mat x=(xr);
arma:: vec sortz,w,newX,z,signw,uv,r,y=(yr),newz;
arma:: vec betaold,beta=(betar),beta0=(betaoldr);
arma::uvec order, index;

int place,n=x.n_rows;
int p=x.n_cols,flag;
double error=10000;
//toler=0.005;
int iteration=1;
//beta= solve(x, y);;

while (iteration<=maxit&& error>toler)
{
  betaold=beta;
  uv=sort(y-x*beta);

  // u is the quantile of given vector
  double quantile=(n-1)*tau-floor((n-1)*tau);
  double u=quantile*uv(ceil((n-1)*tau))+(1-quantile)*uv(floor((n-1)*tau));
		
  r=y-u-x*beta;
  signw=(1-sign(r))/2*(1-tau)+(sign(r)+1)*tau/2;

	for (int j=0;j<p;j++)
	{          	
		z=(r+beta(j)*x.col(j))/x.col(j);
		z.insert_rows(n,1) ;
		order=sort_index(z);
		sortz=z(order);	
		newX=x.col(j)%signw;
                newX.insert_rows(n,1);
                newX(n)=(lamda/beta0(j)/beta0(j));
		w=abs(newX(order));	
                index=find(cumsum(w)>(sum(w)/2),1,"first");
                place=int(index(0));
		beta(j)=sortz(place);
	}
error=sum(abs(beta-betaold));
iteration++;
}

return beta;
}

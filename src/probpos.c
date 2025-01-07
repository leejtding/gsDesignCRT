#include "R.h"
#include "Rmath.h"
double probpos(double theta,int m,double bk,double *z,double *h,double 
Ikm1,double Ik)
{   int i;
    double xhi,prob,rtIk,rtIkm1,mu,rtdeltak;
    mu=theta*(Ik-Ikm1);
    rtdeltak=sqrt(Ik-Ikm1);
    rtIk=sqrt(Ik);
    rtIkm1=sqrt(Ikm1);
    prob=0.;
    for(i=0;i<=m;i++)
    {   xhi=(z[i]*rtIkm1+mu-bk*rtIk)/rtdeltak;
	    prob += pnorm(xhi,0.,1.,1,0)*h[i];
    }
    return(prob);
}

double probpos2(double theta,int m,double bk,double *z,double *h,double 
Ikm1,double Ik)
{   int i;
    double xhi1,xhi2,prob,rtIk,rtIkm1,mu,rtdeltak;
    mu=theta*(Ik-Ikm1);
    rtdeltak=sqrt(Ik-Ikm1);
    rtIk=sqrt(Ik);
    rtIkm1=sqrt(Ikm1);
    prob=0.;
    for(i=0;i<=m;i++)
    {   xhi1=(z[i]*rtIkm1+mu-bk*rtIk)/rtdeltak;
        xhi2=(z[i]*rtIkm1+mu+bk*rtIk)/rtdeltak;
	    prob += (pnorm(xhi1,0.,1.,1,0)+pnorm(xhi2,0.,1.,0,0))*h[i];
    }
    return(prob);
}

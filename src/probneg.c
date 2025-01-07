#include "R.h"
#include "Rmath.h"
double probneg(double theta,int m,double ak,double *z,double *h,double 
Ikm1,double Ik)
{   int i;
    double xlo,prob,rtIk,rtIkm1,mu,rtdeltak;
    mu=theta*(Ik-Ikm1);
    rtdeltak=sqrt(Ik-Ikm1);
    rtIk=sqrt(Ik);
    rtIkm1=sqrt(Ikm1);
    prob=0.;
    for(i=0;i<=m;i++)
    {   xlo=(z[i]*rtIkm1+mu-ak*rtIk)/rtdeltak;
	    prob += pnorm(xlo,0.,1.,0,0)*h[i];
    }
    return(prob);
}

double probneg2(double theta,int m,double ak,double *z,double *h,double 
Ikm1,double Ik)
{   int i;
    double xlo1,xlo2,prob,rtIk,rtIkm1,mu,rtdeltak;
    mu=theta*(Ik-Ikm1);
    rtdeltak=sqrt(Ik-Ikm1);
    rtIk=sqrt(Ik);
    rtIkm1=sqrt(Ikm1);
    prob=0.;
    for(i=0;i<=m;i++)
    {   xlo1=(z[i]*rtIkm1+mu-ak*rtIk)/rtdeltak;
        xlo2=(z[i]*rtIkm1+mu+ak*rtIk)/rtdeltak;
	    prob += (pnorm(xlo1,0.,1.,0,0)-pnorm(xlo2,0.,1.,0,0))*h[i];
    }
    return(prob);
}

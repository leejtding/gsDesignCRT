#define DEBUG 0
/* note: EXTREMEZ > 3 + log(r) +  Z(1-alpha) + Z(1-beta)
   per bottom of p 349 in Jennison and Turnbull */
#define EXTREMEZ 20
#define MAXR 83
#include "R.h"
#include "Rmath.h"
#include "gsDesignCRT.h"
/* Group sequential probability computation per Jennison & Turnbull
   xnanal- # of possible analyses in the group-sequential designs
            (interims + final)
   I     - statistical information available at each analysis
   a     - lower cutoff points for z statistic at each analysis (output)
   b     - upper cutoff points for z statistic at each analysis (output)
   problo- input vector with probability of rejecting (Z<aj) at
           jth interim analysis, j=1...nanal
   probhi- input vector with probability of rejecting (Z>bj) at
           jth interim analysis, j=1...nanal
   xtol  - relative change between iterations required to stop for 'convergence'
   xr    - determinant of # of grid points for numerical integration
           r=17 will give a max of 201 points which is what they recommend
   retval- error flag returned; 0 if convergence; 1 indicates error
   printerr- 1 if error messages to be printed - other values suppress printing
*/
void gsbounds1(int *xnanal, double *xtheta, double *I, double *a, double *b,
               double *problo, double *probhi, double *xtol, int *xr, int *retval,
               int *printerr) {
    int i,i1,i2,j,m11,m12,m21,m22,r,nanal;
    double plo,phi,dplo,dphi,btem=0.,atem=0.,atem2,btem2,rtdeltak,rtIk,rtIkm1,xlo,xhi,theta,mu1,mu2;
	double adelta,bdelta,tol;
    /* note: should allocat zwk & wwk dynamically...*/
    double zwk11[1000],wwk11[1000],hwk11[1000],zwk12[1000],wwk12[1000],hwk12[1000],
           zwk21[1000],wwk21[1000],hwk21[1000],zwk22[1000],wwk22[1000],hwk22[1000],
           *z11,*z12,*w11,*w12,*h11,*h12,*z21,*z22,*w21,*w22,*h21,*h22,*tem,rt2pi;
    void h1(double,int,double *,double,double *,double *);
    void hupdate(double,double *,int,double,double *,double *,int,double,double *,double *);
    int gridpts(int,double,double,double,double *,double *);
    r=xr[0]; nanal=xnanal[0]; theta=xtheta[0]; tol=xtol[0]; rt2pi=2.506628274631;

    /* compute bounds at 1st interim analysis using inverse normal */
    if (nanal<1 || r<1 || r>MAXR)
    {   retval[0]=1;
        if (*printerr)
        {	Rprintf("gsbounds1 error: illegal argument");
            if (nanal<1) Rprintf("; nanal=%d--must be > 0",nanal);
            if (r<1 || r> MAXR) Rprintf("; r=%d--must be >0 and <84",r);
            Rprintf("\n");
        }
        return;
	}
    rtIk=sqrt(I[0]); mu1=0.; mu2=rtIk*theta;
    if (problo[0] <= 0.) a[0] = -EXTREMEZ;
    else a[0]=qnorm(problo[0],mu2,1.,1,0);
    if (probhi[0] <= 0.) b[0] = EXTREMEZ;
    else b[0]=qnorm(probhi[0],mu1,1.,0,0);
    if (nanal==1) {retval[0]=0; return;}

    /* set up work vectors */
    z11=zwk11; w11=wwk11; h11=hwk11; z12=zwk12; w12=wwk12; h12=hwk12;
    z21=zwk21; w21=wwk21; h21=hwk21; z22=zwk22; w22=wwk22; h22=hwk22;
    
    m11=gridpts(r,mu1,a[0],b[0],z11,w11);
    h1(0.,m11,w11,I[0],z11,h11);
    
    m12=gridpts(r,mu2,a[0],b[0],z12,w12);
    h1(theta,m12,w12,I[0],z12,h12);

    /* use Newton-Raphson to find subsequent interim analysis cutpoints */
    if (*printerr) Rprintf("Start: r=%d mu1=%lf mu2=%lf a[0]=%lf b[0]=%lf\n",r,mu1,mu2,a[0],b[0]);
    retval[0]=0;
    for(i=1;i<nanal;i++)
    {   /* set up constants */
        rtIkm1=rtIk; rtIk=sqrt(I[i]); mu2=rtIk*theta; rtdeltak=sqrt(I[i]-I[i-1]);
        if (rtdeltak < 1) rtdeltak=1;
        if (problo[i]<=0.) atem2= -EXTREMEZ;
        else atem2=qnorm(problo[i],mu2,1.,1,0); 
        if (probhi[i]<=0.) btem2= EXTREMEZ;
        else btem2=qnorm(probhi[i],mu1,1.,0,0);
        
        /* find upper boundary */
        bdelta=1.; j=0;
        if (*printerr) Rprintf("desired probhi=%lf\n", probhi[i]);
        while((bdelta>tol) && j++ < EXTREMEZ)
	    {   phi=0.; dphi=0.; btem=btem2;
            if (*printerr) Rprintf("i=%d m11=%d m12=%d\n",i,m11,m12);
	        
            /* construct upper boundary */
            /* compute probability of crossing upper boundaries & their derivatives under H_0 */
            for(i1=0;i1<=m11;i1++)
            {   xhi=(z11[i1]*rtIkm1-btem*rtIk)/rtdeltak;
                phi+=h11[i1]*pnorm(xhi,0.,1.,1,0);
                dphi-=h11[i1]*exp(-xhi*xhi/2)/rt2pi*rtIk/rtdeltak;
            }

            /* use 1st order Taylor's series to update upper boundary under H_0*/
            /* maximum allowed change is 1 */
            /* maximum value allowed is z1[m1]*rtIk to keep within grid points */
            if (*printerr) Rprintf("i=%2d j=%2d btem=%lf phi=%lf dphi=%lf\n",i,j,btem,phi,dphi);       
            bdelta=probhi[i]-phi;
            if (bdelta<dphi) btem2=btem+1.;
            else if (bdelta > -dphi) btem2=btem-1.;
            else btem2=btem+(probhi[i]-phi)/dphi;
            if (btem2>EXTREMEZ) btem2=EXTREMEZ;
            else if (btem2< -EXTREMEZ) btem2= -EXTREMEZ;
            bdelta=btem2-btem; if (bdelta<0) bdelta= -bdelta;
        }
        b[i]=btem;
        
        /* find lower boundary */
        adelta=1.; j=0;
        if (*printerr) Rprintf("desired problo=%lf\n", problo[i]);
        while((adelta>tol) && j++ < EXTREMEZ)
	    {   plo=0.; dplo=0.; atem=atem2;
            if (*printerr) Rprintf("i=%d m11=%d m12=%d\n",i,m11,m12);

            /* construct lower boundary */
            /* compute probability of crossing lower boundaries & their derivatives under H_1 */
            for(i2=0;i2<=m12;i2++)
            {   xlo=(z12[i2]*rtIkm1-atem*rtIk+theta*(I[i]-I[i-1]))/rtdeltak;
                plo+=h12[i2]*pnorm(xlo,0.,1.,0,0);
                dplo+=h12[i2]*exp(-xlo*xlo/2)/rt2pi*rtIk/rtdeltak;
            }

            /* use 1st order Taylor's series to update upper boundary under H_1 */
            /* maximum allowed change is 1 */
            /* maximum value allowed is z1[m1]*rtIk to keep within grid points */
            if (*printerr) Rprintf("i=%2d j=%2d atem=%lf plo=%lf dplo=%lf\n",i,j,atem,plo,dplo);
            adelta=problo[i]-plo;
            if (adelta>dplo) atem2=atem+1.;
            else if (adelta < -dplo) atem2=atem-1.;
            else atem2=atem+(problo[i]-plo)/dplo;
            if (atem2>EXTREMEZ) atem2=EXTREMEZ;
            else if (atem2 < -EXTREMEZ) atem2= -EXTREMEZ;
            adelta=atem2-atem; if (adelta<0) adelta= -adelta;
        }
        a[i]=atem;
        
        /* if convergence did not occur, set flag for return value */
        if (adelta>tol || bdelta > tol)
        {   if (*printerr) 
            {  Rprintf("gsbound error: No convergence for boundary for interim %d; I=%7.0lf",i+1,I[i]);
	   			if (bdelta>tol) Rprintf("\n last 2 upper boundary values: %lf %lf\n",btem,btem2);
				if (adelta>tol) Rprintf("\n last 2 lower boundary values: %lf %lf\n",atem,atem2);
			}
			retval[0]=1;
		}
        if (i<nanal-1)
        {   m21=gridpts(r,mu1,a[i],b[i],z21,w21);
            m22=gridpts(r,mu2,a[i],b[i],z22,w22);
            hupdate(0.,w21,m11,I[i-1],z11,h11,m21,I[i],z21,h21);
            hupdate(theta,w22,m12,I[i-1],z12,h12,m22,I[i],z22,h22);
            m11=m21; m12=m22;
            tem=z11; z11=z21; z21=tem;
            tem=w11; w11=w21; w21=tem;
            tem=h11; h11=h21; h21=tem;
            tem=z12; z12=z22; z22=tem;
            tem=w12; w12=w22; w22=tem;
            tem=h12; h12=h22; h22=tem;
        }
    }
    retval[0]=0; 
	return;
}


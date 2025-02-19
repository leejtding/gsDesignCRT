#include <math.h>
#include "R.h"
void gridpts1(int r,double mu,double *z)
{  int i,r5,r6;
   double rdbl,r2dbl;
   rdbl=r; r2dbl=2*r;
   for(i=1;i<r;i++) z[i-1]=mu-3-4*log(rdbl/i);
   r5=5*r; r6=6*r;
   for(i=r;i<=r5;i++) z[i-1]=mu-3*(1.-(i-rdbl)/r2dbl);
   for(i=r5+1;i<r6;i++) z[i-1]=mu+3.+4.*log(rdbl/(r6-i));
}
/* returns gridpoints per Jennison & Turnbull, p. 349
   returned value is # of grid pts
*/
int gridpts(int r,double mu,double a, double b, double *z,double *w)
{  int i,r5,r6,j=0,done=0;
   double rdbl,r2dbl,ztem;
   rdbl=r; r2dbl=2*r; r6=6*r; r5=5*r; w[0]=0.;
   ztem=mu-3-4*log(rdbl);
   if (ztem<=a) z[0]=a; 
   else if (ztem>=b) {z[0]=b; w[0]=0.; done=1;}
   else z[0]=ztem;
   for(i=2;i<r6 && done==0;i++) 
   {  if (i<r) ztem=mu-3.-4.*log(rdbl/i);
      else if (i<=r5) ztem= mu+3.*(-1+(i-r)/r2dbl);
      else ztem=mu+3.+4.*log(rdbl/(r6-i));
      if (ztem > a)
      {  j+=2;
         z[j]=ztem;
         if(ztem>=b) {z[j]=b; done=1;}
         z[j-1]=(z[j]+z[j-2])/2.;
   }  }
   if (j>0) 
   {   w[0]=(z[2]-z[0])/6.;
       w[j]=(z[j]-z[j-2])/6.;
       w[j-1]=2.*(z[j]-z[j-2])/3.;
   }
   for(i=1;i<j-1;i+=2)
   {  w[i]= 2.*(z[i+1]-z[i-1])/3.;
      w[i+1]=(z[i+3]-z[i-1])/6.;
   }
   return(j);
}
int gridpts2(int r,double mu,double a, double b, double *z,double *w)
{  int i,r5,r6,j=0,j0=0,done=0;
   double rdbl,r2dbl,ztem;
   rdbl=r; r2dbl=2*r; r6=6*r; r5=5*r; w[0]=0.;
   // (-b, -a)
   ztem=mu-3-4*log(rdbl);
   if (ztem<=-b) z[0]=-b; 
   else if (ztem>=-a) {z[0]=-a; w[0]=0.; done=1;}
   else z[0]=ztem;
   for(i=2;i<r6 && done==0;i++) {
      if (i<r) ztem=mu-3.-4.*log(rdbl/i);
      else if (i<=r5) ztem= mu+3.*(-1+(i-r)/r2dbl);
      else ztem=mu+3.+4.*log(rdbl/(r6-i));
      if (ztem > -b) {
         j+=2;
         z[j]=ztem;
         if(ztem>=-a) {z[j]=-a; done=1;}
         z[j-1]=(z[j]+z[j-2])/2.;
      }
   }
   if (j>0) {
      w[0]=(z[2]-z[0])/6.;
      w[j]=(z[j]-z[j-2])/6.;
      w[j-1]=2.*(z[j]-z[j-2])/3.;
   }
   for(i=1;i<j-1;i+=2) {
      w[i]= 2.*(z[i+1]-z[i-1])/3.;
      w[i+1]=(z[i+3]-z[i-1])/6.;
   }
   // (a, b)
   ztem=mu-3-4*log(rdbl); j0=j+1; j=j0; done=0;
   if (ztem<=a) z[j0]=a; 
   else if (ztem>=b) {z[j0]=b; w[j0]=0.; done=1;}
   else z[j0]=ztem;
   for(i=2;i<r6 && done==0;i++) {
      if (i<r) ztem=mu-3.-4.*log(rdbl/i);
      else if (i<=r5) ztem= mu+3.*(-1+(i-r)/r2dbl);
      else ztem=mu+3.+4.*log(rdbl/(r6-i));
      if (ztem > a) {
         j+=2;
         z[j]=ztem;
         if(ztem>=b) {z[j]=b; done=1;}
         z[j-1]=(z[j]+z[j-2])/2.;
      }
   }
   if (j>0) {
      w[j0]=(z[j0+2]-z[j0])/6.;
      w[j]=(z[j]-z[j-2])/6.;
      w[j-1]=2.*(z[j]-z[j-2])/3.;
   }
   for(i=j0+1;i<j-1;i+=2) {
      w[i]= 2.*(z[i+1]-z[i-1])/3.;
      w[i+1]=(z[i+3]-z[i-1])/6.;
   }
   return(j);
}

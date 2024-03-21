void gsbound(int *xnanal,double *I,double *a,double *b,double *problo,double *probhi,
              double *xtol,int *xr,int *retval,int *printerr);
void gsupper1(int *xnanal,double *xtheta,double *I,double *a,double *b,double *problo,
              double *probhi,double *xtol,int *xr,int *retval,int *printerr);
void gsupper2(int *xnanal,double *xtheta,double *I,double *a,double *b,double *problo,
              double *probhi,double *xtol,int *xr,int *retval,int *printerr);
void gsbounds1(int *xnanal, double *xtheta, double *I, double *a, double *b,
               double *problo, double *probhi, double *xtol, int *xr, int *retval,
               int *printerr);
void gsbounds2(int *xnanal, double *xtheta, double *I, double *a, double *b,
               double *problo, double *probhi, double *xtol, int *xr, int *retval,
               int *printerr);
void gsboundsnb1(int *xnanal, double *xtheta, double *I, double *a, double *b,
                 double *problo, double *probhi, double *xtol, int *xr, int *retval,
                 int *printerr);
void gsboundsnb2(int *xnanal, double *xtheta, double *I, double *a, double *b,
                 double *problo, double *probhi, double *xtol, int *xr, int *retval,
                 int *printerr);
void probrej1(int *xnanal,int *ntheta,double *xtheta,double *I,double *a,double *b,
              double *xproblo,double *xprobhi,int *xr);
void probrej2(int *xnanal,int *ntheta,double *xtheta,double *I,double *a,double *b,
              double *xproblo,double *xprobhi,int *xr);
void gsdensity(double *den, int *xnanal, int *ntheta, double *xtheta,
               double *I, double *a, double *b, double *xz, int *zlen, int *xr);
void stdnorpts(int *r,double *bounds,double *z,double *w);

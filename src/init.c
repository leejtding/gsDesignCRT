#include "R.h"
#include "Rmath.h"
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include "gsDesignCRT.h"

/* This file is to register all .C entry points in gsDesignCRT: 
gsupper; gsbounds1; gsbounds2; probrej; gsdensity; stdnorpts */

/*
void gsbound(int *xnanal,double *I,double *a,double *b,double *problo,double *probhi,
             double *xtol,int *xr,int *retval,int *printerr) */

static R_NativePrimitiveArgType gsbound_t[] = {
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP
};

/*
void gsupper1(int *xnanal,double *xtheta,double *I,double *a,double *b,double *problo,
              double *probhi,double *xtol,int *xr,int *retval,int *printerr) */

static R_NativePrimitiveArgType gsupper1_t[] = {
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP
};

/*
void gsupper2(int *xnanal,double *xtheta,double *I,double *a,double *b,double *problo,
             double *probhi,double *xtol,int *xr,int *retval,int *printerr) */

static R_NativePrimitiveArgType gsupper2_t[] = {
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP
};

/*
void gsbounds1(int *xnanal,double *xtheta,double *I,double *a,double *b,double *problo,
               double *probhi,double *xtol,int *xr,int *retval,int *printerr) */

static R_NativePrimitiveArgType gsbounds1_t[] = {
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP
};

/*
void gsbounds2(int *xnanal,double *xtheta,double *I,double *a,double *b,double *problo,
               double *probhi,double *xtol,int *xr,int *retval,int *printerr) */

static R_NativePrimitiveArgType gsbounds2_t[] = {
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP
};

/*
void gsboundsnb1(int *xnanal,double *xtheta,double *I,double *a,double *b,double *problo,
               double *probhi,double *xtol,int *xr,int *retval,int *printerr) */

static R_NativePrimitiveArgType gsboundsnb1_t[] = {
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP
};

/*
void gsboundsnb2(int *xnanal,double *xtheta,double *I,double *a,double *b,double *problo,
               double *probhi,double *xtol,int *xr,int *retval,int *printerr) */

static R_NativePrimitiveArgType gsboundsnb2_t[] = {
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP
};

/*
void probrej1(int *xnanal,int *ntheta,double *xtheta,double *I,double *a,double *b,
              double *xproblo,double *xprobhi,int *xr) */

static R_NativePrimitiveArgType probrej1_t[] = {
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP
};

/*
void probrej2(int *xnanal,int *ntheta,double *xtheta,double *I,double *a,double *b,
              double *xproblo,double *xprobhi,int *xr) */

static R_NativePrimitiveArgType probrej2_t[] = {
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP
};

/*
void gsdensity(double *den, int *xnanal, int *ntheta, double *xtheta,
               double *I, double *a, double *b, double *xz,
               int *zlen, int *xr) */

static R_NativePrimitiveArgType gsdensity_t[] = {
  REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP
};

/*
void stdnorpts(int *r,double *bounds,double *z,double *w) */

static R_NativePrimitiveArgType stdnorpts_t[] = {
  INTSXP, REALSXP, REALSXP, REALSXP
};

/* define array of all C entry points */
static const R_CMethodDef CEntries[] = {
  {"gsbound", (DL_FUNC) &gsbound, 10, gsbound_t},
  {"gsupper1", (DL_FUNC) &gsupper1, 11, gsupper1_t},
  {"gsupper2", (DL_FUNC) &gsupper2, 11, gsupper2_t},
  {"gsbounds1", (DL_FUNC) &gsbounds1, 11, gsbounds1_t},
  {"gsbounds2", (DL_FUNC) &gsbounds2, 11, gsbounds2_t},
  {"gsboundsnb1", (DL_FUNC) &gsboundsnb1, 11, gsboundsnb1_t},
  {"gsboundsnb2", (DL_FUNC) &gsboundsnb2, 11, gsboundsnb2_t},
  {"probrej1", (DL_FUNC) &probrej1, 9, probrej1_t},
  {"probrej2", (DL_FUNC) &probrej2, 9, probrej2_t},
  {"gsdensity", (DL_FUNC) &gsdensity, 10, gsdensity_t},
  {"stdnorpts", (DL_FUNC) &stdnorpts, 4, stdnorpts_t},
  {NULL, NULL, 0, NULL}
};

/* now register in the init function */
void R_init_gsDesign(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);

  /* the DLL is to be searched for entry points specified by character strings as well */
  R_useDynamicSymbols(dll, TRUE);

  /* allow .C calls by character strings: */
  R_forceSymbols(dll, FALSE);
}

/*  <lalVerbatim file="LALInspiralBCVSpinBankCV">
Author: Tagoshi, H, Van Den Broeck, C, Jones, G.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>
\subsection{Module \texttt{LALInspiralCreateCoarseBank.c}}
\subsubsection*{Prototypes}
\subsubsection*{Description}
\subsubsection*{Algorithm}
\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}
\subsubsection*{Notes}
\clearpage
</lalLaTeX>  */


#include <math.h>
#include <lal/AVFactories.h>
#include <lal/FlatMesh.h>
#include <lal/LALConfig.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALMalloc.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALStdlib.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/MatrixUtils.h>
#include <lal/SeqFactories.h>
#include <lal/LALInspiralBCVSpinBank.h>

#ifndef _HT_BCVSPINMETRIC_H_
#define _HT_BCVSPINMETRIC_H_

/* This is a header for double precision complex routine dcomplex.c */
/* H.Tagoshi */

#ifndef _HT_DCOMPLEX_H_
#define _HT_DCOMPLEX_H_

#ifndef _DCOMPLEX_DECLARE_HT_
typedef struct DCOMPLEX {double r,i;} dcomplex;
#define _DCOMPLEX_DECLARE_HT_
#endif /* _DCOMPLEX_DECLARE_HT_ */

static dcomplex DCadd(dcomplex a, dcomplex b);
static dcomplex DCsub(dcomplex a, dcomplex b);
static dcomplex DCmul(dcomplex a, dcomplex b);
static dcomplex DComplex(double re, double im);
static dcomplex DConjg(dcomplex z);
static dcomplex DCdiv(dcomplex a, dcomplex b);
static double DCabs(dcomplex z);
static dcomplex DCsqrt(dcomplex z);
static dcomplex DRCmul(double x, dcomplex a);

#endif /* _HT_DCOMPLEX_H_ */

/*
#include "nr.h"
#include "nrutil.h"

#include "htnr.h"
#include "htnr.c"
#include "dcomplex.h"
#include "dcomplex.c"

*/
static int 
cos_sin_func(
   /* input */
   int N, double beta,double fmax,
   /* output */
   double *costerm,double *sinterm);

static int 
coef_A(
   /* input */
   int N,double *costerm,double *sinterm,double fmax,
   /* output */
   double *A1, double *A2,double *A3);

static int deriv_A(/* input */
	    int N,double *costerm,double *sinterm,double fmax,
	    /* output */
	    double *dA2,double *dA3);

static int innerR(/* input */
	   int N, double *A, double *B, double *Sn,double fmin,double fmax,
	   /* output */
	   double *result);

static int innerC(/* input */
	   int N, dcomplex *A, dcomplex *B, double *Sn,double fmin,double fmax,
	   /* output */
	   double *result);

static int orthonormalized_A(/* input */
		      int N,double *A1, double *A2,double *A3,double *Sn,
		      double fmin,double fmax,
		      /* output */
		      double *tA2,double *tA3,
		      double *hA1,double *hA2, double *hA3,
		      double *normtA2,double *normtA3);

static int dA2dbeta(/* input */
	     int N,double *Sn,double fmin,double fmax,
	     double *hA1,double *A2,
	     double *tA2,double *hA2,double *dA2,double normtA2,
	     /* output */
	     double *dhA2);

static int dA3dbeta(/* input */
	     int N,double *Sn,double fmin,double fmax,
	     double *hA1,double *hA2,double *dhA2,
	     double *A3,double *tA3,double *dA3,double normtA3,
	     /*  output */
	     double *dhA3);

static int calc_function_G(/* input */
               int N,double *Sn,double fmin,double fmax,
	       double *A1,double *A2,double *A3,
	       double *dhA2db,double *dhA3db,
	       /* output */
	       double funcG[7][7][4][4]);

static int functionG(/* input */
	      int N,double beta,double *Sn,double fmin,double fmax,
	      /* output */
	      double funcG[7][7][4][4]);

static int three_metric(/* input */
	      double funcG[7][7][4][4],double *alpha,
	      /* output */
	      double metric3[4][4]);

static int generate_fit_points(/*input*/ double funcG[7][7][4][4],
			int ndata,
			/*output*/ double fit_point[ndata+1][4]);

static int generate_metric_data(/* input */
			 double funcG[7][7][4][4]);

static void model_func(double xx,double afunc[],int ma);

static int metric_by_fit(/* input */
		  int ndata,
		  /* output */
		  double metric_fit[4][4]);

static int rescale_metric(/*input*/ int ndata, double metric1[4][4],
		   /*output*/ double metric[4][4]);

static double determinant3(double **a);

static int matrix2_determinant_plus(double **matrix,double *eig);

static int matrix3_determinant_plus(double **matrix,double *eig);

static double innerp(int n,double *b,double *c);

static double vector_product(double *a,double *b,double *c);

static int product_matrix(int n,double **A,double **B,double **C);

static int product_mat_vec(int n,double *w,double **A,double *v);

#endif

/* <lalVerbatim file="LALInspiralBCVSpinBankCP"> */
void
LALInspiralBCVSpinBank(
    LALStatus         	 *status,
    SnglInspiralTable   **tiles,
    INT4      		 *ntiles,
    InspiralCoarseBankIn *coarseIn
    )
/* </lalVerbatim> */

{
	double a=2.;
	printf("%e \n", a);
}

/* This is double precision complex computation routines */
/* H.Tagoshi */

#ifndef _HT_DCOMPLEX_C_
#define _HT_DCOMPLEX_C_

#include <math.h>

dcomplex DCadd(dcomplex a, dcomplex b)
{
	dcomplex c;
	c.r=a.r+b.r;
	c.i=a.i+b.i;
	return c;
}

dcomplex DCsub(dcomplex a, dcomplex b)
{
	dcomplex c;
	c.r=a.r-b.r;
	c.i=a.i-b.i;
	return c;
}


dcomplex DCmul(dcomplex a, dcomplex b)
{
	dcomplex c;
	c.r=a.r*b.r-a.i*b.i;
	c.i=a.i*b.r+a.r*b.i;
	return c;
}

dcomplex DComplex(double re, double im)
{
	dcomplex c;
	c.r=re;
	c.i=im;
	return c;
}

dcomplex DConjg(dcomplex z)
{
	dcomplex c;
	c.r=z.r;
	c.i = -z.i;
	return c;
}

dcomplex DCdiv(dcomplex a, dcomplex b)
{
	dcomplex c;
	double r,den;
	if (fabs(b.r) >= fabs(b.i)) {
		r=b.i/b.r;
		den=b.r+r*b.i;
		c.r=(a.r+r*a.i)/den;
		c.i=(a.i-r*a.r)/den;
	} else {
		r=b.r/b.i;
		den=b.i+r*b.r;
		c.r=(a.r*r+a.i)/den;
		c.i=(a.i*r-a.r)/den;
	}
	return c;
}

double DCabs(dcomplex z)
{
	double x,y,ans,temp;
	x=fabs(z.r);
	y=fabs(z.i);
	if (x == 0.0)
		ans=y;
	else if (y == 0.0)
		ans=x;
	else if (x > y) {
		temp=y/x;
		ans=x*sqrt(1.0+temp*temp);
	} else {
		temp=x/y;
		ans=y*sqrt(1.0+temp*temp);
	}
	return ans;
}

dcomplex DCsqrt(dcomplex z)
{
	dcomplex c;
	double x,y,w,r;
	if ((z.r == 0.0) && (z.i == 0.0)) {
		c.r=0.0;
		c.i=0.0;
		return c;
	} else {
		x=fabs(z.r);
		y=fabs(z.i);
		if (x >= y) {
			r=y/x;
			w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
		} else {
			r=x/y;
			w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
		}
		if (z.r >= 0.0) {
			c.r=w;
			c.i=z.i/(2.0*w);
		} else {
			c.i=(z.i >= 0) ? w : -w;
			c.r=z.i/(2.0*c.i);
		}
		return c;
	}
}

dcomplex DRCmul(double x, dcomplex a)
{
	dcomplex c;
	c.r=x*a.r;
	c.i=x*a.i;
	return c;
}


#endif

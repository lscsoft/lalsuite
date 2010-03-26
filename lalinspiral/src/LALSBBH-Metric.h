/* Header file for SBBH-Metric.c
   ver 3.0   2005/3/8 H.Tagoshi
ver 3.1   2005/3/21 H.Tagoshi
ver6.5 2006/2/23  NR functions was replaced with GSL functions H. Takahashi
*/
/* This is a header for double precision complex routine dcomplex.c */
/* H.Tagoshi */

#ifndef _HT_DCOMPLEX_H_
#define _HT_DCOMPLEX_H_

#include        <string.h>

#ifndef _DCOMPLEX_DECLARE_HT_
typedef struct DCOMPLEX {double r,i;} dcomplex;
#define _DCOMPLEX_DECLARE_HT_
#endif /* _DCOMPLEX_DECLARE_HT_ */

dcomplex DCadd(dcomplex a, dcomplex b);
dcomplex DCsub(dcomplex a, dcomplex b);
dcomplex DCmul(dcomplex a, dcomplex b);
dcomplex DComplex(double re, double im);
dcomplex DConjg(dcomplex z);
dcomplex DCdiv(dcomplex a, dcomplex b);
double DCabs(dcomplex z);
dcomplex DCsqrt(dcomplex z);
dcomplex DRCmul(double x, dcomplex a);

#endif /* _HT_DCOMPLEX_H_ */

#ifndef _HT_BCVSPINMETRIC_H_
#define _HT_BCVSPINMETRIC_H_

#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_roots.h>

int cos_sin_func(/* input */
		 int N, double beta,double fmax,
		 /* output */
		 double *costerm,double *sinterm);

int coef_A(/* input */
	   int N,double *costerm,double *sinterm,double fmax,
	   /* output */
	   double *A1, double *A2,double *A3);

int deriv_A(/* input */
	    int N,double *costerm,double *sinterm,double fmax,
	    /* output */
	    double *dA2,double *dA3);

int innerR(/* input */
	   int N, double *A, double *B, double *Sn,double fmin,double fmax,
	   /* output */
	   double *result);

int innerC(/* input */
	   int N, dcomplex *A, dcomplex *B, double *Sn,double fmin,double fmax,
	   /* output */
	   double *result);

int orthonormalized_A(/* input */
		      int N,double *A1, double *A2,double *A3,double *Sn,
		      double fmin,double fmax,
		      /* output */
		      double *tA2,double *tA3,
		      double *hA1,double *hA2, double *hA3,
		      double *normtA2,double *normtA3);

int dA2dbeta(/* input */
	     int N,double *Sn,double fmin,double fmax, double *hA1,
	     double *tA2,double *dA2,double normtA2,
	     /* output */
	     double *dhA2);

int dA3dbeta(/* input */
	     int N,double *Sn,double fmin,double fmax,
	     double *hA1,double *hA2,double *dhA2,
	     double *A3,double *tA3,double *dA3,double normtA3,
	     /*  output */
	     double *dhA3);

int calc_function_G(/* input */
               int N,double *Sn,double fmin,double fmax,
	       double *A1,double *A2,double *A3,
	       double *dhA2db,double *dhA3db,
	       /* output */
	       double funcG[7][7][4][4]);

int functionG(/* input */
	      int N,double beta,double *Sn,double fmin,double fmax,
	      /* output */
	      double funcG[7][7][4][4]);

int three_metric(/* input */
	      double funcG[7][7][4][4],double *alpha,
	      /* output */
	      double metric3[4][4]);

int generate_fit_points(/*input*/double MinMatch, double funcG[7][7][4][4],
			int ndata,
			/*output*/ double **fit_point);

int generate_metric_data(/* input */double MinMatch,
			 double funcG[7][7][4][4]);

void model_func(double xx,double afunc[]);

int metric_by_fit(/* input */
		   double MinMatch, int ndata,
		  /* output */
		  double metric_fit[4][4]);

int rescale_metric(/*input*/ double MinMatch, int ndata, double metric1[4][4],
		   /*output*/ double metric[4][4]);

int BCVspin_metric(/*input*/
		double MinMatch, int N,double *Sn,double fmin,double fmax,double beta,
		/*output*/
		double **bcv2metric,int dbg);

double determinant3(gsl_matrix *matrix);

int matrix3_determinant_plus(gsl_matrix *matrix,gsl_vector *eig);

void svdfit_d_test(double x[], double y[], double sig[], int ndata, gsl_vector *a, int ma,
	gsl_matrix *u, gsl_matrix  *v, gsl_vector *w, double *chisq,
		void (*funcs)(double, double []));

#endif


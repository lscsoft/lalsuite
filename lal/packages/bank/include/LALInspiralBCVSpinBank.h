/*
*  Copyright (C) 2007 Bernd Machenschalk, B.S. Sathyaprakash, Chris Van Den Broeck
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/* Header file for BCVspin_metric_fit.c
   ver 3.0   2005/3/8 H.Tagoshi
ver 3.1   2005/3/21 H.Tagoshi
ver6.5 2006/2/23  NR functions was replaced with GSL functions H. Takahashi
*/

/* <lalVerbatim file="LALInspiralBCVSpinBankHV">

Author: Tagoshi, H., T. Hirotaka, Van Den Broeck, C., Jones, G., B.S. Sathyaprakash
$Id$

</lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{LALInspiralBCVSpinBank.h}}
\label{s:LALInspiralBCVSpinBank.h}

Header file for the template placement codes.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALInspiralBCVSpinBank.h>
\end{verbatim}

\noindent This header file covers routines that are used in template placement
for spinning black hole binaries using phenomenological BCV templates.

</lalLaTeX> */

#include <lal/LIGOMetadataTables.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>

/* This is a header for double precision complex routine dcomplex.c */
/* H.Tagoshi */

#ifndef _HT_DCOMPLEX_H_
#define _HT_DCOMPLEX_H_

#include        <stdio.h>
#include        <stdlib.h>
#include        <string.h>
#include        <math.h>

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

#include <lal/LALRCSID.h>
NRCSID (LALINSPIRALBCVSPINBANKH,"$Id$");

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
			/*output*/ double fit_point[ndata+1][4]);

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
		double bcv2metric[4][4],int dbg);

double determinant3(gsl_matrix *matrix);

int matrix2_determinant_plus(gsl_matrix *matrix,gsl_vector *eig);

int matrix3_determinant_plus(gsl_matrix *matrix,gsl_vector *eig);

double innerp(int n,double *b,double *c);

double vector_product(double *a,double *b,double *c);

int product_matrix(int n,double A[][n+1],double B[][n+1],double C[][n+1]);

int product_mat_vec(int n,double *w,double A[][n+1],double *v);

int BCVspin_spacing(double MinMatch,double metric3[4][4],double a[4][4],double *deltax);

int BCVspin_effmetric(/* input */
		      double MinMatch, double metric3[4][4],double a[4][4],
		      /* output */
		      double effmetric[3][3]);

struct func1_params
  {
	int N;
	double *Sn;
	double fmin;
	double fmax;
	double MinMatch;
};

double func1(double beta,void *params);

struct func2_params
{
	double beta_i;
	int N;
	double *Sn;
	double fmin;
	double fmax;
	double MinMatch;
};

double func2(double beta,void *params);

int BCVspin_beta_placement(double MinMatch,double beta_min,double beta_max, int N,
		double *Sn,double fmin,double fmax, double *beta_list,int *nbeta);

int BCVspin_beta_placement_effmetric(/* input*/ double MinMatch,double beta_min,double beta_max, int N,
		double *Sn,double fmin,double fmax,
		/* output */ double effmetric_list[3][3][1001], double *beta_list,int *nbeta);


void svdfit_d_test(double x[], double y[], double sig[], int ndata, gsl_vector *a, int ma,
	gsl_matrix *u, gsl_matrix  *v, gsl_vector *w, double *chisq,
		void (*funcs)(double, double []));

#endif


/*
*  Copyright (C) 2007 B.S. Sathyaprakash, Thomas Cokelaer, Chris Van Den Broeck
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

/*  <lalVerbatim file="LALInspiralBCVSpinBankCV">
Authors: Tagoshi, H, Takahashi, H, Van Den Broeck, C, Jones, G, Sathyaprakash, BS
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>
\subsection{Module \texttt{LALInspiralBCVSpinBank.c}}
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
#include <lal/SeqFactories.h>
#include <lal/LALInspiralBCVSpinBank.h>
#include <gsl/gsl_linalg.h>

NRCSID(LALINSPIRALBCVSPINBANKC, "$Id$");

/* <lalVerbatim file="LALInspiralBCVSpinBankCP"> */
void
LALInspiralBCVSpinBank(
    LALStatus         	 *status,
    SnglInspiralTable    **tiles,
    INT4      		 *ntiles,
    InspiralCoarseBankIn *coarseIn
    )
/* </lalVerbatim> */

{
  /* Nmax is the maximum number of beta values we are allowed to have */
  INT4 Nmax=1000;
  INT4 N, k, nbeta, totTemps=0;
  REAL8 *Sn, MM, betaMin, betaMax;
  REAL8Vector newSn;
  REAL8 fmin, fmax, beta_list[Nmax+1];
  REAL8 effmetric_list[3][3][Nmax+1];

  /* Let us declare the metric, bank parameters and vector sequences that will
   * hold the psi0-psi3 values at the location of the templates
   */
  static InspiralMetric metric;
  static InspiralBankParams bankParams;
  static CreateVectorSequenceIn in;
  static REAL4VectorSequence *list=NULL;
  static REAL4VectorSequence *totList=NULL;
  SnglInspiralTable *bank=NULL, *tmpBank=NULL;

  INITSTATUS(status, "LALInspiralBCVSpinBank", LALINSPIRALBCVSPINBANKC );
  ATTATCHSTATUSPTR( status );

  ASSERT( coarseIn != NULL, status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL );
  ASSERT( *tiles == NULL, status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL );

  /*
   * Create the structure necessary for vector sequence to
   * temporarily store the values of psi0-psi3-beta.
   * list will store the values of psi0-psi3 at each beta level.
   * totList will store the values of psi0-psi3-beta which will
   * be used to fill the output structure tiles.
   */
  in.length = 1;
  in.vectorLength = 2;

  LALSCreateVectorSequence( status->statusPtr, &list, &in );
  CHECKSTATUSPTR( status );

  LALSCreateVectorSequence( status->statusPtr, &totList, &in );
  CHECKSTATUSPTR( status );

  /*
   * coarseIn structure knows about the range of beta, psi0,
   * and psi3. It also provides the power spectrum of noise,
   * minimal match (in mmCoarse), upper and lower freq cutoffs.
   */

  MM = coarseIn->mmCoarse;
  fmin = coarseIn->fLower;
  fmax = coarseIn->fUpper;
  betaMin = coarseIn->betaMin;
  betaMax = coarseIn->betaMax;
  bankParams.minimalMatch = coarseIn->mmCoarse;
  bankParams.x0Min = coarseIn->psi0Min;
  bankParams.x0Max = coarseIn->psi0Max;
  bankParams.x1Min = coarseIn->psi3Min;
  bankParams.x1Max = coarseIn->psi3Max;

  /* We really don't know the metric; this is just to create a
   * place holder
   */
  bankParams.metric = &metric;

  /*
   * Instead of noisespec(N, Sn, fmin, fmax) let us use the noise
   * PSD provided by coarseIn.
   */

  Sn = coarseIn->shf.data->data;
  N = coarseIn->shf.data->length-1;
  newSn.length=0;

#if 1
  if (N>32768)
  {
    int ratio, i;
    ratio = N/32768;
    N = N / ratio;
    LALInfo(status,  "entering BCVSpin metric computation using the smooth PSD \n");
    newSn.length= N;

    newSn.data  = (REAL8*) LALCalloc(1, sizeof(REAL8) * (newSn.length+1));
    for (i=0; i<N; i++)
    {
      int j;
      for (j=0;j<ratio; j++)
      {
        newSn.data[i]+=  coarseIn->shf.data->data[i*ratio+j];
      }
    newSn.data[i]/=ratio;
    }
    newSn.data[N] = coarseIn->shf.data->data[N*ratio];
    Sn = newSn.data;

  }
else
{
    LALInfo(status,  "entering BCVSpin metric computation using the whole PSD\n");
}
#endif

/*  printf("Num beta steps=%d Sn[0]=%e, Sn[N/2]=%e, MM=%f\n", N, Sn[0], Sn[N/2], MM);*/

  /* Use Tagoshi's code to get values of beta b/w betaMin and betaMax at
   * the given minimal match MM and the effective metric at the corresponding
   * points.
   */
  BCVspin_beta_placement_effmetric(MM, betaMin, betaMax, N, Sn, fmin, fmax,
		  effmetric_list, beta_list, &nbeta);

  /*
  printf("Num beta steps=%d\n", nbeta);
  */
  /*
   * There are nbeta levels of beta.
   * Next, obtain the psi0-psi3 bank at each value of beta
   */

  for(k=1; k<=nbeta; k++)
  {
	  double a, b, c, det, q;
	  int j, ndx;

	  /* a, b, c are the three independent values of the metric
	   */
	  a = metric.G00 = effmetric_list[1][1][k];
	  b = metric.G01 = effmetric_list[1][2][k];
	  c = metric.G11 = effmetric_list[2][2][k];


	  /* Diagonalize the metric and store the diagonalized values in the metric */

	  det = a * c - b * b;
	  q = sqrt( (a-c)*(a-c) + 4. * b*b );
	  metric.g00 = 0.5 * (a + c - q);
	  metric.g11 = 0.5 * (a + c + q);

	/* 	THOMAS UPDATES !!!!!!!!!!! I have the feeling that the bank
		is over efficiency by a factor 3 with respect to classical
		BCV search. Here as a first go, I would use a factor of sqrt(3)
		in the metric component, based on comparison between BCV bank
		and BCVspin bank ffor different MM form 0.8 to 0.97 which shows
		 a constant ration of 3.1

 	metric.g00/=sqrt(3.);
 	metric.g11/=sqrt(3.);
	*/
	  if ( a == c )
	  {
		  metric.theta = LAL_PI/2.;
	  }
	  else
	  {
	  /* metric->theta = 0.5 * atan(2.*b/(a-c)); We want to always
	   * measure the angle from the semi-major axis to the tau0 axis
	   * which is given by the following line as opposed to the line above
	   */
		  metric.theta = atan( b / (metric.g00 - c) );
	  }

	  /* Now we are ready to call Thomas' BCV template bank code
	   * to compute the psi0-psi3 values depending on the value of LowGM.
	   */

	  if (coarseIn->insidePolygon == True)
	  {
		  LALInspiralCreateFlatBankS3S4 (status->statusPtr, list, &bankParams, *coarseIn);
	  }
	  else
	  {
		  LALInspiralCreateFlatBank (status->statusPtr, list, &bankParams);
	  }
	  CHECKSTATUSPTR( status );

	  /* Optionally output the metric on the screen
	   */
	  fprintf(stderr, "k=%d beta=%e g00=%e g11=%e theta=%e Nf=%d\n",
			  k, beta_list[k], metric.g00, metric.g11, metric.theta, list->length);

	  /* Create additional memory to add the psi0-psi3 values at the current
	   * value of beta to the psi0-psi3 values stored from previous values of beta.
	   */
	  ndx = list->length;
	  totTemps += ndx;
	  totList->data = (REAL4 *) LALRealloc( totList->data, (3*totTemps) * sizeof(REAL4) );

	  if ( !totList->data )
	  {
		  ABORT(status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
	  }
	  for (j=0; j<ndx; j++)
	  {
		  int m1 = 2*j;
		  int m2 = 3*j + 3*totTemps-3*ndx;
		  totList->data[m2+0] = list->data[m1+0];
		  totList->data[m2+1] = list->data[m1+1];
		  totList->data[m2+2] = (float) (beta_list[k]);
		  /*
		  printf("%e %e %e %e\n",
		  totList->data[m2], totList->data[m2+1], totList->data[m2+2], x);
		  */
	  }
  }
  /*
  for (k=0; k<totTemps; k++)
  {
	  printf("%e %e %e\n", totList->data[3*k], totList->data[3*k+1], totList->data[3*k+2]);
  }
  */

  /* Convert output data structure. */

  tmpBank = bank = (SnglInspiralTable *) LALCalloc(1, sizeof(SnglInspiralTable));
  if (bank == NULL){
	  ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
  }
  *tiles = bank;
  for( k = 0; k < totTemps; k++ )
  {
	  bank = bank->next = (SnglInspiralTable *) LALCalloc( 1, sizeof( SnglInspiralTable ) );
	  if (bank == NULL)
	  {
		  ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
	  }
	  bank->psi0 = totList->data[3*k+0];
	  bank->psi3 = totList->data[3*k+1];
	  bank->beta = totList->data[3*k+2];
	  bank->f_final = fmax;
  }

  /* Free tiles template, which is blank. */

  *ntiles = totTemps;
  bank = (*tiles)->next;
  LALFree( *tiles );
  *tiles = bank;

  if (newSn.length!=0)
  {
	LALFree(newSn.data);
  }

  LALFree (totList->data);
  LALFree (list->data);
  LALFree (totList);
  LALFree (list);
  DETATCHSTATUSPTR(status);
  RETURN (status);
}



/* BCVspin_metric_fit.c    (by H.Tagoshi and H.Takahashi)

   The program to calculate the metric.
   This routine requre Numerical Recipes library.

   The one sided noise power spectrum density Sn is
   array Sn[0..N]. The frequency interval, df, is df=fmax/N.
   Sn[0] is for f=0 [Hz]. Sn[N] is for f=fmax [Hz].
   Sn[i] must be non-zero from i=fmin/df to i=N.

   The metric is given as bcv2metric[1..3][1..3]
   by BCVspin_metric() routine.
   The index 1 is psi_0, index 2 is psi_3/2,
   and index 3 is beta.

   At ver3.1, this program is working with gcc,
   and PGI C compiler (pgcc)
   (although a lot of warning message appear by pgcc).
   However, this does not work with Intel C Compiler(IA32).
   Thus, this must be compliled by gcc (and probably by PGI C).

   2004/9/22 Revised version for metric part.

   2005/3/7  ver3.0 spacing part added.

   2005/3/21 ver3.1 bug fix

    2006/2/23 ver6.5 NR functions was replaced with GSL functions

*/

#ifndef _HT_BCVSPINMETRIC_C_
#define _HT_BCVSPINMETRIC_C_

#define PI (3.141592653589793238462643383279502)

#define DTRENORM (200)  /* The renormalization factor for the time derivative
                           of the template. This is introduced to have
			   an agreement with Yanbei's mathematica code.
                           This value does not affect the results.
                           This can be set to 1. */
#define N_RANDOM (100)  /* The number of monte carlo to generate
                           the 6 phase factor */
#define JJ (150)       /* The number of direction of contour
                           for which the fitting is done. */

double cont_data[JJ+1][4];

int cos_sin_func(/* input */
		 int N, double beta,double fmax,
		 /* output */
		 double *costerm,double *sinterm)
{
  int i;
  double df,freq,freq2;

  df=fmax/N;

  costerm[0]=0.; sinterm[0]=0.;
  for(i=1;i<=N;i++){
    freq=df*i;
    freq2=pow(freq,-2/3.);
    costerm[i]=cos(beta*freq2);
    sinterm[i]=sin(beta*freq2);
  }

  return 0;
}

int coef_A(/* input */
	   int N,double *costerm,double *sinterm,double fmax,
	   /* output */
	   double *A1, double *A2,double *A3)
{
  int i;
  double df,freq;

  df=fmax/N;

  A1[0]=0.;A2[0]=0.;A3[0]=0.;
  for(i=1;i<=N;i++){
    freq=df*i;
    A1[i]=pow(freq,-7/6.);
    A2[i]=A1[i]*costerm[i];
    A3[i]=A1[i]*sinterm[i];
  }

  return 0;
}


int deriv_A(/* input */
	    int N,double *costerm,double *sinterm,double fmax,
	    /* output */
	    double *dA2,double *dA3)
{
  int i;
  double df,freq,freq2;

  df=(fmax)/N;

  dA2[0]=0.;dA3[0]=0.;
  for(i=1;i<=N;i++){
    freq=df*i;
    freq2=pow(freq,-11/6.);
    dA2[i]=-freq2*sinterm[i];
    dA3[i]=freq2*costerm[i];
  }

  return 0;
}

/* inner product of real functions A(f) and B(f).
   A[0...N], B[0...N]
*/

int innerR(/* input */
	   int N, double *A, double *B, double *Sn,double fmin,double fmax,
	   /* output */
	   double *result)
{

  int i,imin;
  double df,tmp;

  df=(fmax)/N;
  imin=fmin/df;

  tmp=(A[imin]*B[imin]/Sn[imin]+A[N]*B[N]/Sn[N])/2.;
  for(i=imin+1;i<N;i++)
    tmp+=A[i]*B[i]/Sn[i];

  *result=tmp*4*df;

  return 0;
}

/* inner product of complex functions A(f) and B(f).
   A[0...N], B[0...N]
*/

int innerC(/* input */
	   int N, dcomplex *A, dcomplex *B, double *Sn,double fmin,double fmax,
	   /* output */
	   double *result)
{

  int i,imin;
  double df,tmp;

  df=(fmax)/N;
  imin=fmin/df;

  tmp = (A[imin].r * (B[imin].r) + A[imin].i * (B[imin].i))/Sn[imin]/2.;
  tmp+= (A[N].r * B[N].r + A[N].i * B[N].i)/Sn[N]/2.;

  for(i=imin+1;i<N;i++){
    /* here we take the multiplication of conjuguate of A by B and keep only the real part.*/
    tmp+= (A[i].r * B[i].r + A[i].i * B[i].i)/Sn[i];
  }
  *result=tmp*df*4;
  return 0;
}


int orthonormalized_A(/* input */
		      int N,double *A1, double *A2,double *A3,double *Sn,
		      double fmin,double fmax,
		      /* output */
		      double *tA2,double *tA3,
		      double *hA1,double *hA2, double *hA3,
		      double *normtA2,double *normtA3)
{
  int i;
  double tmp,a,b,d;
  double A1A1,A1A2,A1A3,hA2A3;

  innerR(N,A1,A1,Sn,fmin,fmax,&A1A1);
  innerR(N,A1,A2,Sn,fmin,fmax,&A1A2);
  innerR(N,A1,A3,Sn,fmin,fmax,&A1A3);

  a=sqrt(A1A1);
  for(i=0;i<=N;i++)
    hA1[i]=A1[i]/a;

  b=A1A2/a;
  for(i=0;i<=N;i++){
    tA2[i]=A2[i]-b*hA1[i];
  }

  innerR(N,tA2,tA2,Sn,fmin,fmax,&tmp);
  *normtA2=sqrt(tmp);

  for(i=0;i<=N;i++)
    hA2[i]=tA2[i]/(*normtA2);

  innerR(N,hA2,A3,Sn,fmin,fmax,&hA2A3);

  d=A1A3/a;
  for(i=0;i<=N;i++){
    tA3[i]=A3[i]-d*hA1[i]-hA2A3*hA2[i];
  }

  innerR(N,tA3,tA3,Sn,fmin,fmax,&tmp);
  *normtA3=sqrt(tmp);

  for(i=0;i<=N;i++)
    hA3[i]=tA3[i]/(*normtA3);

  return 0;
}

int dA2dbeta(/* input */
	     int N,double *Sn,double fmin,double fmax, double *hA1,
	     double *tA2,double *dA2,double normtA2,
	     /* output */
	     double *dhA2)
{
  int i;
  double inner,inner2,*dtA2, norm;

  dtA2 = (double *) malloc(sizeof(double) *(N+1));
  innerR(N,dA2,hA1,Sn,fmin,fmax,&inner);

  for(i=0;i<=N;i++)
    dtA2[i]=dA2[i]-inner*hA1[i];

  innerR(N,tA2,dtA2,Sn,fmin,fmax,&inner2);

  norm = pow(normtA2, 3.);
  for(i=0;i<=N;i++)
    dhA2[i]=dtA2[i]/normtA2-tA2[i]*inner2/norm;

  free(dtA2);
  return 0;
}

int dA3dbeta(/* input */
	     int N,double *Sn,double fmin,double fmax,
	     double *hA1,double *hA2,double *dhA2,
	     double *A3,double *tA3,double *dA3,double normtA3,
	     /*  output */
	     double *dhA3)
{
  int i;
  double inner,dhA2A3,hA2A3,*dtA3;
  double hA1dA3,hA2dA3, norm;

  dtA3 = (double *) malloc(sizeof(double)*(N+1));

  innerR(N,hA1,dA3,Sn,fmin,fmax,&hA1dA3);

  innerR(N,hA2,dA3,Sn,fmin,fmax,&hA2dA3);

  innerR(N,dhA2,A3,Sn,fmin,fmax,&dhA2A3);

  innerR(N,hA2,A3,Sn,fmin,fmax,&hA2A3);

  for(i=0;i<=N;i++)
    dtA3[i]=dA3[i]-hA1dA3*hA1[i]-hA2dA3*hA2[i]-dhA2A3*hA2[i]-
      hA2A3*dhA2[i];

  innerR(N,tA3,dtA3,Sn,fmin,fmax,&inner);
  norm = pow(normtA3, 3.);
  for(i=0;i<=N;i++)
    dhA3[i]=dtA3[i]/normtA3-tA3[i]*inner/norm;

  free(dtA3);
  return 0;
}

/* calc_function_G */
/* A1,A2,A3 must be orthonormalized */
/* dhA2db=d(A2)/d(beta), d() is pertial derivative */
/* funcG[][][][] :  function G_{ij\alpha\beta} */

int calc_function_G(/* input */
               int N,double *Sn,double fmin,double fmax,
	       double *A1,double *A2,double *A3,
	       double *dhA2db,double *dhA3db,
	       /* output */
	       double funcG[7][7][4][4])
{
  int i,j,k,l,m,imin;
  double freq,df,*freq2,sum;
  double inn1[7][7][4][4],inn2[7][7][4];


  /*
  dcomplex u[7][4][N+1],A[7][N+1],ctmp;
  */
  dcomplex *u[7][4],*A[7],ctmp;

  for(i=0;i<=6;i++)
  for(j=0;j<=3;j++)
  {
     u[i][j] = (dcomplex *) malloc(sizeof(dcomplex)*(N+1));
  }
  for(i=0;i<=6;i++)
  {
     A[i] = (dcomplex *) malloc(sizeof(dcomplex)*(N+1));
  }

  df=(fmax)/N;
  imin=fmin/df;

  for(i=0;i<=N;i++){
    A[1][i]=DComplex(A1[i],0.);
    A[2][i]=DComplex(A2[i],0.);
    A[3][i]=DComplex(A3[i],0.);
    A[4][i]=DComplex(0.,A1[i]);
    A[5][i]=DComplex(0.,A2[i]);
    A[6][i]=DComplex(0.,A3[i]);
  }

  for(j=1;j<=6;j++){
    for(i=0;i<=N;i++){
      freq=df*i;
      ctmp=DComplex(0.,2*PI*freq/DTRENORM);
      u[j][0][i]=DCmul(ctmp,A[j][i]);
    }
  }

  for(i=0;i<=N;i++){
    u[1][3][i]=DComplex(0.,0.);
    u[2][3][i]=DComplex(dhA2db[i],0.);
    u[3][3][i]=DComplex(dhA3db[i],0.);
    u[4][3][i]=DComplex(0.,0.);
    u[5][3][i]=DComplex(0.,dhA2db[i]);
    u[6][3][i]=DComplex(0.,dhA3db[i]);
  }

  freq2 = (double *) malloc(sizeof(double)*(N+1));
  freq2[0]=0.;
  for(i=1;i<=N;i++){
    freq=df*i;
    freq2[i]=pow(freq,-5./3.);
  }

  for(j=1;j<=6;j++){
    for(i=0;i<=N;i++){
      ctmp=DComplex(0.,freq2[i]);
      u[j][1][i]=DCmul(ctmp,A[j][i]);
    }
  }

  freq2[0]=0.;
  for(i=1;i<=N;i++){
    freq=df*i;
    freq2[i]=pow(freq,-2./3.);
  }
  for(j=1;j<=6;j++){
    for(i=0;i<=N;i++){
      ctmp=DComplex(0.,freq2[i]);
      u[j][2][i]=DCmul(ctmp,A[j][i]);
    }
  }
  free(freq2);

  for(i=1;i<=6;i++)
    for(j=1;j<=6;j++)
      for(k=0;k<=3;k++)
	for(l=0;l<=3;l++){
	  innerC(N,u[i][k],u[j][l],Sn,fmin,fmax,&sum);
	  inn1[i][j][k][l]=sum;
	}

  for(i=1;i<=6;i++)
    for(j=1;j<=6;j++)
      for(k=0;k<=3;k++){
	innerC(N,u[i][k],A[j],Sn,fmin,fmax,&sum);
	inn2[i][j][k]=sum;
      }

  for(i=1;i<=6;i++)
    for(j=1;j<=6;j++)
      for(k=0;k<=3;k++)
	for(l=0;l<=3;l++){
	  sum=0;
	  for(m=1;m<=6;m++)
	    sum+=inn2[i][m][k]*inn2[j][m][l];

	  funcG[i][j][k][l]=inn1[i][j][k][l]-sum;
	}

  for(i=0;i<=6;i++)
  for(j=0;j<=3;j++)
  {
     free(u[i][j]);
  }
  for(i=0;i<=6;i++)
  {
     free(A[i]);
  }

  return 0;
}

int functionG(/* input */
	      int N,double beta,double *Sn,double fmin,double fmax,
	      /* output */
	      double funcG[7][7][4][4])
{
  double *A1,*A2,*A3,*tA2,*tA3,*hA1,*hA2,*hA3;
  double *costerm,*sinterm;
  double *dA2,*dA3,*dhA2,*dhA3;
  double normtA2,normtA3;


  A1 = (double *) malloc(sizeof(double)*(N+1));
  A2 = (double *) malloc(sizeof(double)*(N+1));
  A3 = (double *) malloc(sizeof(double)*(N+1));
  dA2 = (double *) malloc(sizeof(double)*(N+1));
  dA3 = (double *) malloc(sizeof(double)*(N+1));
  costerm = (double *) malloc(sizeof(double)*(N+1));
  sinterm = (double *) malloc(sizeof(double)*(N+1));

  cos_sin_func(N,beta,fmax,costerm,sinterm);
  coef_A(N,costerm,sinterm,fmax,A1,A2,A3);
  deriv_A(N,costerm,sinterm,fmax,dA2,dA3);

  free(costerm);
  free(sinterm);

  tA2 = (double *) malloc(sizeof(double)*(N+1));
  tA3 = (double *) malloc(sizeof(double)*(N+1));
  hA1 = (double *) malloc(sizeof(double)*(N+1));
  hA2 = (double *) malloc(sizeof(double)*(N+1));
  hA3 = (double *) malloc(sizeof(double)*(N+1));
  dhA2 = (double *) malloc(sizeof(double)*(N+1));
  dhA3 = (double *) malloc(sizeof(double)*(N+1));

  orthonormalized_A(N,A1,A2,A3,Sn,fmin,fmax,tA2,tA3,hA1,hA2,hA3,
		    &normtA2,&normtA3);
  free(A1);
  free(A2);

  dA2dbeta(N,Sn,fmin,fmax,hA1,tA2,dA2,normtA2,dhA2);
  free(tA2);
  free(dA2);

  dA3dbeta(N,Sn,fmin,fmax,hA1,hA2,dhA2,A3,tA3,dA3,normtA3,dhA3);
  free(A3);
  free(dA3);
  free(tA3);

  calc_function_G(N,Sn,fmin,fmax,hA1,hA2,hA3,dhA2,dhA3,funcG);

  free(hA1);
  free(hA2);
  free(hA3);
  free(dhA2);
  free(dhA3);

  return 0;
}

int three_metric(/* input */
	      double funcG[7][7][4][4],double *alpha,
	      /* output */
	      double metric3[4][4])
{
  int i,j,l,m,a,b;
  double tmp1,tmp2,tmp3;

  tmp3=0;
  for(i=1;i<=6;i++)
    for(j=1;j<=6;j++)
      tmp3+=alpha[i]*alpha[j]*funcG[i][j][0][0];

  for(a=1;a<=3;a++)
    for(b=1;b<=3;b++){

      tmp1=0;
      for(i=1;i<=6;i++)
	for(j=1;j<=6;j++)
	  tmp1+=alpha[i]*alpha[j]*funcG[i][j][a][b];

      tmp2=0;
      for(i=1;i<=6;i++)
	for(j=1;j<=6;j++)
	  for(l=1;l<=6;l++)
	    for(m=1;m<=6;m++)
	      tmp2+=alpha[i]*alpha[j]*funcG[i][j][0][a]
		*alpha[l]*alpha[m]*funcG[l][m][0][b];

      metric3[a][b]=0.5*(tmp1-tmp2/tmp3);
    }

  return 0;
}

int generate_fit_points(/*input*/double MinMatch, double funcG[7][7][4][4],
			int ndata,
			/*output*/ double fit_point[JJ+1][4])
{
  const gsl_rng_type * T;
  gsl_rng * r;
  int i,j,k;
  double x[4],xd[4],sum,alpha[7],metric3[4][4];

  /*  for(i=1;i<=6;i++)
    for(j=1;j<=6;j++)
      for(k=0;k<=3;k++)
	printf("%e %e %e %e\n",funcG[i][j][k][0],funcG[i][j][k][1],
	funcG[i][j][k][2],funcG[i][j][k][3]);*/

	gsl_rng_env_setup();

    T = gsl_rng_default;
     r = gsl_rng_alloc (T);


	gsl_matrix *a = gsl_matrix_alloc (3, 3);
	gsl_vector *eig = gsl_vector_alloc (3);
    gsl_matrix *V = gsl_matrix_alloc (3, 3);

  sum=0;
  for(i=1;i<=6;i++){
		alpha[i]=gsl_rng_uniform_pos (r)-0.5;
		sum+=alpha[i]*alpha[i];
  }
  sum=sqrt(sum);
  for(i=1;i<=6;i++){
    alpha[i]=alpha[i]/sum;
  }
  three_metric(funcG,alpha,metric3);


  for(i=1;i<=3;i++)
    for(j=1;j<=3;j++)
			gsl_matrix_set (a, i-1, j-1,metric3[i][j]);

	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (3);

   gsl_eigen_symmv (a, eig,V,w);

	gsl_eigen_symmv_free (w);

	matrix3_determinant_plus(V,eig);

	for(i=1;i<=ndata;i++){
    x[1]=gsl_rng_uniform_pos (r)-0.5;
    x[2]=gsl_rng_uniform_pos (r)-0.5;
    x[3]=gsl_rng_uniform_pos (r)-0.5;
    sum=x[1]*x[1]+x[2]*x[2]+x[3]*x[3];
		sum=sqrt(sum);

    xd[1]=x[1]/sum*sqrt((1.-MinMatch)/gsl_vector_get (eig, 0));
    xd[2]=x[2]/sum*sqrt((1.-MinMatch)/gsl_vector_get (eig, 1));
	xd[3]=x[3]/sum*sqrt((1.-MinMatch)/gsl_vector_get (eig, 2));

    for(j=1;j<=3;j++){
      fit_point[i][j]=0;
      for(k=1;k<=3;k++)
	fit_point[i][j]+=gsl_matrix_get (V, j-1, k-1) *xd[k];
    }
  }

	  gsl_rng_free (r);
	gsl_vector_free(eig);
	gsl_matrix_free(a);
	gsl_matrix_free(V);

	return 0;
}

int generate_metric_data(/* input */double MinMatch,
			 double funcG[7][7][4][4])
{
  extern double cont_data[JJ+1][4];
  int i,k,i2,j2;
  double x[4],maxr;
  double alpha[7],metric3[N_RANDOM+1][4][4],sum,rr[N_RANDOM+1],distance;
  double fit_point[JJ+1][4], norm;
  const gsl_rng_type * T;
  gsl_rng * r;


  generate_fit_points(MinMatch,funcG,JJ,fit_point);

	gsl_rng_env_setup();

  T = gsl_rng_default;
   r = gsl_rng_alloc (T);

  for(k=1;k<=N_RANDOM;k++){
    sum=0;
    for(i2=1;i2<=6;i2++){
			alpha[i2]=gsl_rng_uniform_pos (r)-0.5;
			sum+=alpha[i2]*alpha[i2];
    }
    sum=sqrt(sum);
    for(i2=1;i2<=6;i2++){
      alpha[i2]=alpha[i2]/sum;
    }
    three_metric(funcG,alpha,metric3[k]);
  }

  norm = sqrt(1.-MinMatch);

  for(i=1;i<=JJ;i++){
    x[1]=fit_point[i][1];
    x[2]=fit_point[i][2];
    x[3]=fit_point[i][3];

    for(k=1;k<=N_RANDOM;k++){
    distance=0;
    for(i2=1;i2<=3;i2++)
      for(j2=1;j2<=3;j2++)
	distance+=metric3[k][i2][j2]*x[i2]*x[j2];
    rr[k]=sqrt(distance);
    }

    maxr=rr[1];
    for(k=2;k<=N_RANDOM;k++)
      if(rr[k]>maxr) maxr=rr[k];

    cont_data[i][1]=x[1]/maxr*norm;
    cont_data[i][2]=x[2]/maxr*norm;
    cont_data[i][3]=x[3]/maxr*norm;
  }

	gsl_rng_free (r);

	return 0;
}

void model_func(double xx,double afunc[])
{
  int i,j,k;
  double x[4];
  extern double cont_data[JJ+1][4];

  x[1]=cont_data[(int)xx][1];
  x[2]=cont_data[(int)xx][2];
  x[3]=cont_data[(int)xx][3];

  k=1;
  for(i=1;i<=3;i++)
    for(j=1;j<=3;j++){
      afunc[k]=(x[i]*x[j]);
      k++;
    }

}

int metric_by_fit(/* input */
		  double MinMatch,int ndata,
		  /* output */
		  double metric_fit[4][4])
{
  int i,j,k,ma=9,ia[ma+1];
  double x[ndata+1],y[ndata+1],sig[ndata+1];
  double chisq;

	gsl_matrix *u = gsl_matrix_alloc (ndata,ma);
	gsl_matrix *v = gsl_matrix_alloc (ma,ma);
	gsl_vector *w = gsl_vector_alloc (ma);
	gsl_vector *a= gsl_vector_alloc (ma);

  for(i=1;i<=ndata;i++){
    x[i]=i;
    y[i]=1.-MinMatch;
    sig[i]=1.;
  }
  for(i=1;i<=ma;i++){
    ia[i]=1.;
  }

  svdfit_d_test(x,y,sig,ndata,a,ma,u,v,w,&chisq,&model_func);

  k=1;
  for(i=1;i<=3;i++)
    for(j=1;j<=3;j++){
      metric_fit[i][j]=gsl_vector_get(a,k-1);
      k++;
    }

	gsl_vector_free(w);
	gsl_vector_free(a);
	gsl_matrix_free(u);
	gsl_matrix_free(v);


  return 0;
}

int rescale_metric(/*input*/double MinMatch, int ndata, double metric1[4][4],
		   /*output*/ double metric[4][4])
{
  int i,j,k;
  double distance,dmin;

  dmin=100;
  for(k=1;k<=ndata;k++){
    distance=0;
    for(i=1;i<=3;i++)
      for(j=1;j<=3;j++)
	distance+=metric1[i][j]*cont_data[k][i]*cont_data[k][j];
    if(distance<dmin) dmin=distance;
  }

  for(i=1;i<=3;i++)
    for(j=1;j<=3;j++)
      metric[i][j]=metric1[i][j]*(1.-MinMatch)/dmin;

  return 0;
}

/* Calculate the BCVspin metric.
   Input:
   MinMatch : The minimal match of the template space
   N            : The number of frequency bin between f=0..fmax .
   Sn[0...N]    : one sided power spectrum density of noise [1/sqrt[Hz]].
   fmin         : minimum frequency cutoff for Sn(f).
   fmax         : maximum frequency for Sn(f).
   beta         : beta value of BCVspin template.
   dbg          : debuging flag. if dbg=1, some internal information
                  will be displayed. Otherwise, nothing is displayed.

   Output:
   bcv2metric[1..3][1..3]
                : The metric of the BCVspin template.
		: The coordinate system is (\psi_0,\psi_{3/2},\beta).
*/
int BCVspin_metric(/*input*/
		double MinMatch, int N,double *Sn,double fmin,double fmax,double beta,
		/*output*/
		double bcv2metric[4][4],int dbg)
{
  extern double cont_data[JJ+1][4];
  int i,j;
  double funcG[7][7][4][4];
  double metric_data3[JJ+1][4];
  double metric_fit[4][4],metric_rescaled[4][4];

  functionG(N,beta,Sn,fmin,fmax,funcG);

  if(dbg==1) fprintf(stderr,"metric data .....................\n");
  generate_metric_data(MinMatch,funcG);

  for(i=1;i<=JJ;i++){
    metric_data3[i][1]=cont_data[i][1];
    metric_data3[i][2]=cont_data[i][2];
    metric_data3[i][3]=cont_data[i][3];
  }

  if(dbg==1) fprintf(stderr,"metric fit .....................\n");

  metric_by_fit(MinMatch,JJ,metric_fit);

  if(dbg==1)
    for(i=1;i<=3;i++)
      fprintf(stderr,"%e %e %e\n",
	      metric_fit[i][1],metric_fit[i][2],metric_fit[i][3]);

  rescale_metric(MinMatch,JJ,metric_fit,metric_rescaled);

  for(i=1;i<=3;i++)
    for(j=i;j<=3;j++)
      bcv2metric[i][j]=metric_rescaled[i][j];
  for(i=2;i<=3;i++)
    for(j=1;j<=i-1;j++)
      bcv2metric[i][j]=metric_rescaled[j][i];

  return 0;
}

/* determinant of 3x3 matrix a[1..3][1..3] */

double determinant3(gsl_matrix *matrix)
{
	int i,j;
  double tmp, a[4][4];

	for (i = 0; i < 3; i++){
		for (j = 0; j < 3; j++){
			a[i+1][j+1]= gsl_matrix_get (matrix, i, j);
		}
	}
  tmp=a[1][1]*a[2][2]*a[3][3];
  tmp+=a[1][2]*a[2][3]*a[3][1];
  tmp+=a[1][3]*a[2][1]*a[3][2];
  tmp-=a[1][3]*a[2][2]*a[3][1];
  tmp-=a[1][2]*a[2][1]*a[3][3];
  tmp-=a[1][1]*a[2][3]*a[3][2];

  return tmp;

}

int matrix2_determinant_plus(gsl_matrix *matrix,gsl_vector *eig)
{
  int i,j;
  double a[3][3],det;

	for (i = 0; i < 2; i++){
		for (j = 0; j < 2; j++){
			a[i+1][j+1]= gsl_matrix_get (matrix, i, j);
		}
	}

  det=a[1][1]*a[2][2]-a[1][2]*a[2][1];

  if(det<0){

		gsl_matrix *V = gsl_matrix_alloc (2, 2);
		gsl_vector *v= gsl_vector_alloc (2);

		for (i = 0; i < 2; i++){
      gsl_vector_set (v, i, gsl_vector_get(eig, i));
		}

	for (i = 0; i < 2; i++){
		for (j = 0; j < 2; j++){
			gsl_matrix_set (V, i, j,gsl_matrix_get (matrix, i, j));
		}
		}


		   gsl_vector_set (eig, 0, gsl_vector_get(v, 1));
    		gsl_vector_set (eig, 1, gsl_vector_get(v, 0));


		for(j=0;j<2;j++){
			gsl_matrix_set (matrix, j, 0,gsl_matrix_get (V, j, 1));
			gsl_matrix_set (matrix, j, 1,gsl_matrix_get (V, j, 0));
		}

    gsl_matrix_free(V);
    gsl_vector_free(v);
	}

  return 0;
}

int matrix3_determinant_plus(gsl_matrix *matrix,gsl_vector *eig)
{
  int i,j;
  double det;


  det=determinant3(matrix);

  if(det<0){
		gsl_matrix *V = gsl_matrix_alloc (3, 3);
		gsl_vector *v= gsl_vector_alloc (3);

		for (i = 0; i < 3; i++){
      gsl_vector_set (v, i, gsl_vector_get(eig, i));
		}

		for (i = 0; i < 3; i++){
		for (j = 0; j < 3; j++){
			gsl_matrix_set (V, i, j,gsl_matrix_get (matrix, i, j));
		}
		}


		    gsl_vector_set (eig, 0, gsl_vector_get(v, 0));
			gsl_vector_set (eig, 1, gsl_vector_get(v, 2));
    		gsl_vector_set (eig, 2, gsl_vector_get(v, 1));

		for(j=0;j<3;j++){
			gsl_matrix_set (matrix, j, 0,gsl_matrix_get (V, j, 0));
			gsl_matrix_set (matrix, j, 1,gsl_matrix_get (V, j, 2));
		    gsl_matrix_set (matrix, j, 2,gsl_matrix_get (V, j, 1));
		}

    gsl_matrix_free(V);
    gsl_vector_free(v);
	}

  return 0;
}


/* Compute the inner product of two vectors */
double innerp(int n,double *b,double *c)
{
  int i;
  double sum;

  sum=0.;
  for(i=1;i<=n;i++)
    sum+=b[i]*c[i];

  return sum;
}

/* Compute the 3-dim vector product */
double vector_product(double *a,double *b,double *c)
{
  a[1]=b[2]*c[3]-b[3]*c[2];
  a[2]=b[3]*c[1]-b[1]*c[3];
  a[3]=b[1]*c[2]-b[2]*c[1];

  return 0;
}


/* compute the product of matrix, B and C, BC. Result is A*/
int product_matrix(int n,double A[][n+1],double B[][n+1],double C[][n+1])
{
  int i,j;
  double b[4][4],c[4][4];

  for(i=1;i<=n;i++)
    for(j=1;j<=n;j++){
      b[i][j]=B[i][j];
      c[i][j]=C[j][i];
    }

  for(i=1;i<=n;i++)
    for(j=1;j<=n;j++)
      A[i][j]=innerp(n,b[i],c[j]);

  return 0;
}

/*Compute the product of a matrix A and a vector v, Av. Result is *w */
int product_mat_vec(int n,double *w,double A[][n+1],double *v)
{
  int i,j;

  for(i=1;i<=n;i++){
    w[i]=0;
    for(j=1;j<=n;j++)
      w[i]+=A[i][j]*v[j];
  }

  return 0;
}

/* BCVspin_spacing()
   Compute the spacing of the template bank.
   Input:
   MinMatch : The minimal match of the template space
   bcvspinmetric[1..3][1..3]: 3-dim metric determined by BCVspin_metric().

   Output:
   deltax[1..3] : (Delta psi_0,Delta psi_3/2,Delta beta)
   a[1..3][1..3]: a[1],a[2],a[3] are vectors which define the inscribed
                  tetrahedron to the equal-match ellipsoid.
*/
int BCVspin_spacing(/* input */
		    double MinMatch, double bcvspinmetric[4][4],
		    /* output */
		    double a[4][4],double *deltax)
{
  int i,j;
  double e[4][4],ad[4][4],f[4][4],P[4][4],Q[4][4],iQP[4][4],PQ[4][4],iP[4][4],iQ[4][4];
  double t1;

  gsl_matrix *mat = gsl_matrix_alloc (3, 3);
  gsl_vector *eig = gsl_vector_alloc (3);
  gsl_matrix *PP = gsl_matrix_alloc (3, 3);

  for(i=1;i<=3;i++)
    for(j=1;j<=3;j++)
      e[i][j]=0;

  for(i=1;i<=3;i++)
    e[i][i]=1;

  for(i=1;i<=3;i++)
    for(j=1;j<=3;j++)
      gsl_matrix_set (mat, i-1, j-1,bcvspinmetric[i][j]);

  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (3);

  gsl_eigen_symmv (mat, eig,PP,w);

  gsl_eigen_symmv_free (w);

  matrix3_determinant_plus(PP,eig);

  for(i=1;i<=3;i++)
  {
    for(j=1;j<=3;j++)
    {
      P[i][j]=gsl_matrix_get (PP, i-1, j-1);
    }
  }


  for(i=1;i<=3;i++)
    for(j=1;j<=3;j++){
      Q[i][j]=0;
      iQ[i][j]=0;
      iP[i][j]=P[j][i];
    }
  for(i=1;i<=3;i++){
    Q[i][i]=1/sqrt(gsl_vector_get (eig, i-1));
    iQ[i][i]=sqrt(gsl_vector_get (eig, i-1));
  }

  product_matrix(3,PQ,P,Q);
  product_matrix(3,iQP,iQ,iP);

  for(i=1;i<=2;i++)
    product_mat_vec(3,f[i],iQP,e[i]);

  vector_product(f[3],f[1],f[2]);

  t1=sqrt(innerp(3,f[1],f[1]));

  for(i=1;i<=3;i++)
    ad[1][i]=f[1][i]/t1;

  t1=sqrt(innerp(3,f[3],f[3]));

  for(i=1;i<=3;i++)
    ad[3][i]=f[3][i]/t1;

  vector_product(ad[2],ad[3],ad[1]);

  for(i=1;i<=3;i++)
    for(j=1;j<=3;j++)
      ad[i][j]*=2./sqrt(3.)*sqrt(1.-MinMatch);

  for(i=1;i<=3;i++)
    product_mat_vec(3,a[i],PQ,ad[i]);

  for(i=1;i<=3;i++)
    deltax[i]=innerp(3,a[i],e[i]);

  gsl_vector_free(eig);
  gsl_matrix_free(PP);
  gsl_matrix_free(mat);

  return 0;
}

/* BCVspin_effmetric()
   Compute the effective 2-dimensional metric of the template bank
   on (\psi_0,\psi_{3/2}) plane.
   This is obtained by rescaling the 3-dim BCVspin metric.

   Input:
   MinMatch : The minimal match of the template space
   bcvspinmetric[1..3][1..3]: 3-dim metric determined by BVC2_metric().
   a[1..3][1..3]: 3x3 matrix determined by BCVspin_spacing().
                  a[1],a[2],a[3] are vectors which define the inscribed
                  tetrahedron to the equal-match ellipsoid.

   Output:
   effmetric[1..2][1..2] : the effective 2-dimensional metric
                           of the template bank on (\psi_0,\psi_{3/2}) plane.
*/
int BCVspin_effmetric(/* input */
		      double MinMatch,double bcvspinmetric[4][4],double a[4][4],
		      /* output */
		      double effmetric[3][3])
{
  int i,j,dim=2;
  double e[3][3],b[3][3],bd[3][3],f[3][3],P[3][3],Q[3][3],iQP[3][3],PQ[3][3],iP[3][3],iQ[3][3];
  double t1;

  gsl_matrix *mat = gsl_matrix_alloc (2, 2);
  gsl_vector *eig = gsl_vector_alloc (2);
  gsl_matrix *PP = gsl_matrix_alloc (2, 2);

  for(i=1;i<=dim;i++)
    for(j=1;j<=dim;j++)
      e[i][j]=0;

  for(i=1;i<=dim;i++)
    e[i][i]=1;

  for(i=1;i<=2;i++)
    for(j=1;j<=2;j++)
      gsl_matrix_set (mat, i-1, j-1,bcvspinmetric[i][j]);


  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (2);

  gsl_eigen_symmv (mat, eig,PP,w);

  gsl_eigen_symmv_free (w);

  matrix2_determinant_plus(PP,eig);

  for(i=1;i<=2;i++){
    for(j=1;j<=2;j++){
	P[i][j]=gsl_matrix_get (PP, i-1, j-1);
			}
	}


  for(i=1;i<=dim;i++)
    for(j=1;j<=dim;j++){
      Q[i][j]=0;
      iQ[i][j]=0;
      iP[i][j]=P[j][i];
    }
  for(i=1;i<=dim;i++){
    Q[i][i]=1/sqrt(gsl_vector_get (eig, i-1));
    iQ[i][i]=sqrt(gsl_vector_get (eig, i-1));
  }

  product_matrix(dim,PQ,P,Q);
  product_matrix(dim,iQP,iQ,iP);

  product_mat_vec(dim,f[1],iQP,e[1]);

  for(i=1;i<=2;i++)
    f[2][i]=e[2][i]-f[1][i]/innerp(2,f[1],f[1])*innerp(2,f[1],e[2]);

  t1=sqrt(innerp(dim,f[1],f[1]));

  for(i=1;i<=dim;i++)
    bd[1][i]=f[1][i]/t1;

  t1=sqrt(innerp(dim,f[2],f[2]));

  for(i=1;i<=dim;i++)
    bd[2][i]=f[2][i]/t1;

  for(i=1;i<=dim;i++)
    for(j=1;j<=dim;j++)
      bd[i][j]*=sqrt(2.)*sqrt(1.-MinMatch);




  for(i=1;i<=dim;i++)
    product_mat_vec(dim, b[i], PQ, bd[i]);

  for(i=1;i<=dim;i++)
    for(j=1;j<=dim;j++)
      effmetric[i][j]=bcvspinmetric[i][j]*innerp(dim,b[1],b[1])
	/(a[1][1]*a[1][1]+a[1][2]*a[1][2]);

  /*
  printf("\n");
  printf("bcvspinmetric\n");
  for(i=1;i<=dim;i++)
    printf("%e %e\n",bcvspinmetric[i][1],bcvspinmetric[i][2]);


  printf("\n");
  printf("effmetric\n");
  for(i=1;i<=dim;i++)
    printf("%e %e\n",effmetric[i][1],effmetric[i][2]);
  */

	gsl_vector_free(eig);
	gsl_matrix_free(PP);
	gsl_matrix_free(mat);


  return 0;
}


double func1(double beta,void *params)
{
	struct func1_params *p
    = (struct func1_params *) params;

	int N = p->N;
	double *Sn =p->Sn;
    double fmin = p ->fmin;
	double fmax = p -> fmax;
	double MinMatch=p->MinMatch;

	 double bcv2metric[4][4],a[4][4],deltax[4];

  BCVspin_metric(MinMatch,N,Sn,fmin,fmax,beta,bcv2metric,0);
  BCVspin_spacing(MinMatch,bcv2metric,a,deltax);

  return beta-deltax[3]/2;
}


double func2(double beta,void *params)
{

		struct func2_params *p
    = (struct func2_params *) params;

	double beta_i=p->beta_i;
	int N = p->N;
	double *Sn =p->Sn;
    double fmin = p ->fmin;
	double fmax = p -> fmax;
	double MinMatch= p -> MinMatch;
  double bcv2metric[4][4],a[4][4],deltax[4], del_beta_i,del_beta;

  BCVspin_metric(MinMatch,N,Sn,fmin,fmax,beta_i,bcv2metric,0);
  BCVspin_spacing(MinMatch,bcv2metric,a,deltax);
  del_beta_i=deltax[3]/2;
  BCVspin_metric(MinMatch,N,Sn,fmin,fmax,beta,bcv2metric,0);
  BCVspin_spacing(MinMatch,bcv2metric,a,deltax);
  del_beta=deltax[3]/2;

  return beta_i+del_beta_i+del_beta-beta;
}


int BCVspin_beta_placement(double MinMatch,double beta_min,double beta_max,
			   int N,double *Sn,double fmin,double fmax,
			   double *beta_list,int *nbeta)
{
  int i,j,Nmax=1000;
  double bcv2metric[4][4],a[4][4],deltax[4],tmplist[Nmax+1],tmp;


  tmp=beta_min;
  if(beta_min==0)
    tmp+=20;

  BCVspin_metric(MinMatch,N,Sn,fmin,fmax,tmp,bcv2metric,0);
  BCVspin_spacing(MinMatch,bcv2metric,a,deltax);

   int status;
   int iter = 0, max_iter = 500;

  const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double r = 0.;
  double x_lo = tmp, x_hi = tmp+deltax[3];
	gsl_function F;


  struct func1_params params = {N,Sn,fmin,fmax,MinMatch};

  F.function = &func1;
  F.params = &params;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);

	do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi,
                                       0., 0.001);

      if (status == GSL_SUCCESS)

		tmplist[1]=r;
    }


	while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free (s);

  for(i=2;i<=Nmax;i++){

    BCVspin_metric(MinMatch,N,Sn,fmin,fmax,tmplist[i-1],bcv2metric,0);
    BCVspin_spacing(MinMatch,bcv2metric,a,deltax);

    if(tmplist[i-1]+deltax[3]/2>beta_max) break;

    r = 0.;
    x_lo = tmplist[i-1], x_hi = tmplist[i-1]+deltax[3]*1.5;

    struct func2_params params = {tmplist[i-1],N,Sn,fmin,fmax,MinMatch};

    F.function = &func2;
     F.params = &params;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);


	do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi,0., 0.001);


      if (status == GSL_SUCCESS)
		tmplist[i]=r;
    }

	while (status == GSL_CONTINUE && iter < max_iter);

		gsl_root_fsolver_free (s);
  }


  *nbeta=i-1;
  if(i>Nmax){
    fprintf(stderr,"Something is wrong in BCVspin_beta_placement.\n");
    fprintf(stderr,"The number of beta points exceeds the limit.\n");
    fprintf(stderr,"i=%d Nmax=%d\n",i,Nmax);
    exit(1);
  }

  for(j=1;j<=*nbeta;j++)
    beta_list[j]=tmplist[j];

  return 0;
}


int BCVspin_beta_placement_effmetric(/* input*/ double MinMatch,double beta_min,double beta_max, int N,
		double *Sn,double fmin,double fmax,
		/* output */ double effmetric_list[3][3][1001], double *beta_list,int *nbeta){
	int i,j ,k ,nbeta2;
	double beta, bcv2metric[4][4],a[4][4],deltax[4],effmetric[3][3];

	BCVspin_beta_placement(MinMatch,beta_min,beta_max,N,Sn,fmin,fmax,
			beta_list,&nbeta2);

	*nbeta=nbeta2;

	for(k=1;k<=nbeta2;k++){
	beta=beta_list[k];
        BCVspin_metric(MinMatch,N,Sn,fmin,fmax,beta,bcv2metric,0);
 	BCVspin_spacing(MinMatch,bcv2metric,a,deltax);
	BCVspin_effmetric(MinMatch,bcv2metric,a,effmetric);

		for(i=1;i<=2;i++){
		for(j=1;j<=2;j++){
				effmetric_list[i][j][k] = effmetric[i][j];
			}
		}
	}

return 0;
}


#define TOL 1.0e-10

void svdfit_d_test(double x[], double y[], double sig[], int ndata, gsl_vector *a, int ma,
	gsl_matrix *u, gsl_matrix *v, gsl_vector *w, double *chisq,
	void (*funcs)(double, double []))
{
	int j,i;
	double wmax,tmp,thresh,sum,afunc[ma+1];

	gsl_vector *b= gsl_vector_alloc (ndata);


	for (i=1;i<=ndata;i++) {
		(*funcs)(x[i],afunc);
		tmp=1.0/sig[i];
		for (j=1;j<=ma;j++) {
			gsl_matrix_set(u,i-1,j-1,afunc[j]*tmp);
		}
		gsl_vector_set(b,i-1,y[i]*tmp);
	}

	gsl_vector *ws= gsl_vector_alloc (ma);
	gsl_linalg_SV_decomp (u, v,w, ws);

	wmax=0.0;
	for (j=1;j<=ma;j++)
		if (gsl_vector_get(w,j-1) > wmax) wmax=gsl_vector_get(w,j-1) ;
	thresh=TOL*wmax;
	for (j=1;j<=ma;j++)
	if (gsl_vector_get(w,j-1)  < thresh) gsl_vector_set(w,j-1,0.0);


	gsl_linalg_SV_solve (u, v,w, b, a);

	*chisq=0.0;
	for (i=1;i<=ndata;i++) {
		(*funcs)(x[i],afunc);
		for (sum=0.0,j=1;j<=ma;j++) sum += gsl_vector_get(a, j-1)*afunc[j];
		*chisq += (tmp=(y[i]-sum)/sig[i],tmp*tmp);
	}
	gsl_vector_free(b);
}
#undef TOL

#undef PI
#undef DTRENORM
#undef N_RANDOM
#undef JJ

#endif

/* example_spacing5.c     2006.3.3  H. Takahashi

   The example program to calculate the value of beta
   of the template space with given beta_min and beta_max .
   This program also calculates the 2dim (\psi_0, \Psi_{3/2}) effective metric
   at the beta value which is calculated above.

   The one sided noise power spectrum density Sn is
   array Sn[0..N]. The frequency interval, df, is df=fmax/N.
   Sn[0] is for f=0 [Hz]. Sn[N] is for f=fmax [Hz].
   Sn[i] must be non-zero from i=fmin/df to i=N.

*/


int noisespec(int N,double *Sn,double fmin,double fmax)
{
  int i,imin;
  double freq,df,f0;

  f0=150.;
  df=fmax/N;
  imin=fmin/df;

  for(i=0;i<=N;i++){
    if(i>=imin){
      freq=df*i;
      Sn[i]=10.0*(pow(4.49*freq/f0,-56.)+0.16*pow(freq/f0,-4.52)+0.52
		  +0.32*pow(freq/f0,2.));
    }
    else
      Sn[i]=0;
  }

  return 0;
}


/* This is double precision complex computation routines */
/* H.Tagoshi */

#ifndef _HT_DCOMPLEX_C_
#define _HT_DCOMPLEX_C_

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

/*  <lalVerbatim file="LALInspiralBCVSpinBankCV">
Authors: Ian Harry, B.S. Sathyaprakash,  Chris Van Den Broeck
$Id$
</lalVerbatim>  */
/*  <lalLaTeX>
\subsection{Module \texttt{LALInspiralBCVSpinRandomBank.c}}
This code uses the random bank developed by I. Harry, B. Allen and B.S. Sathyaprakash (2007, in preparation).
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
#include <gsl/gsl_linalg.h>
#include <time.h>
#include <lal/LALInspiralRandomBank.h>
#include <lal/LALSBBH-Metric.h>


NRCSID(LALINSPIRALBCVSPINRANDOMBANKC, "$Id$");

/* <lalVerbatim file="LALInspiralBCVSpinRandomBankCP"> */
void
LALInspiralBCVSpinRandomBank(
    LALStatus         	 *status,
    SnglInspiralTable    **tiles,
    INT4      		 *ntiles,
    InspiralCoarseBankIn *coarseIn
    )
/* </lalVerbatim> */

{
  INT4 ShMaxSz, N, i, j, k;
  double *Sn;
  REAL8Vector newSn;
  MCBankIn bankIn;
  float *MCbank;
  SnglInspiralTable *bank=NULL;
  MetricFunc *thisMetric;

  INITSTATUS(status, "LALInspiralBCVSpinRandomBank", LALINSPIRALBCVSPINRANDOMBANKC );
  ATTATCHSTATUSPTR( status );

  ASSERT( coarseIn != NULL, status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL );
  ASSERT( *tiles == NULL, status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL );

  /*
   * coarseIn structure knows about the range of beta, psi0,
   * and psi3. It also provides the power spectrum of noise,
   * minimal match (in mmCoarse), upper and lower freq cutoffs.
   */

  thisMetric = &BCVSpinMetric;

  /* The seed should be initialized to a more meaningful random number! */
  bankIn.error = 0.;
  bankIn.seed = coarseIn->iseed;
  bankIn.dim=3;
  bankIn.verbose=0;
  bankIn.nIni=coarseIn->nTIni;
  ShMaxSz = coarseIn->ShMaxSz;
  /*
  EstimateNumberOfTemplates(thisMetric, &bankIn);
  */
  bankIn.nFin=0;

  bankIn.p = (float *) malloc(sizeof(float) * bankIn.dim);
  bankIn.pMax = (float *) malloc(sizeof(float) * bankIn.dim);
  bankIn.pMin = (float *) malloc(sizeof(float) * bankIn.dim);


  /* we will use a one-d array to store the bank. The convention is
   * that the first bankIn.dim-elements of bank store the coordinates of
   * the first template, the bankIn.dim+1 element is reserved for a tag to
   * indicate if the random template is retained or rejected. We
   * therefore use an initial array of size nIni*(bankIn.dim+1).
   */

  MCbank = (float *) malloc(sizeof(float)*bankIn.nIni*(bankIn.dim+1));

  /* We first convert minimal match to maximum mismatch, as
  the random bank uses the latter, and we multiply by a
  factor of 2 because the Tagoshi metric gives an overefficient
  bank. This corresponds to increasing the distance between
  templates by a factor of sqrt(2).
  */

  bankIn.MM=(1-coarseIn->mmCoarse)*2.;
  bankIn.pMin[0] = coarseIn->psi0Min;
  bankIn.pMax[0] = coarseIn->psi0Max;
  bankIn.pMin[1] = coarseIn->psi3Min;
  bankIn.pMax[1] = coarseIn->psi3Max;
  bankIn.pMin[2] = coarseIn->betaMin;
  bankIn.pMax[2] = coarseIn->betaMax;

  Sn = coarseIn->shf.data->data;
  N = coarseIn->shf.data->length-1;
  newSn.length=0;
  newSn.data=NULL;


#if 1
  if (N>ShMaxSz)
  {
    int ratio, m;
    ratio = N/ShMaxSz;
    N = N / ratio;
    LALInfo(status,  "entering BCVSpin metric computation using the smooth PSD \n");
    newSn.length= N;

    newSn.data  = (REAL8*) LALCalloc(1, sizeof(REAL8) * (newSn.length+1));
    for (m=0; m<N; m++)
    {
      int l;
      for (l=0;l<ratio; l++)
      {
        newSn.data[m]+=  coarseIn->shf.data->data[m*ratio+l];
      }
    newSn.data[m]/=ratio;
    }
    newSn.data[N] = coarseIn->shf.data->data[N*ratio];
    Sn = newSn.data;

  }
else
{
    LALInfo(status,  "entering BCVSpin metric computation using the whole PSD\n");
}
#endif

  bankIn.Npsd = N;
  bankIn.Sn = Sn;
  /* Call the routine to generate the bank */
  MonteCarloBank(thisMetric, &bankIn, MCbank);



  if (bankIn.verbose)
  {
    fprintf(stdout, "Numtemplates=%d MaxMismatch=%e\n", bankIn.nFin, bankIn.MM);
    for (i=0; i<bankIn.nIni; i++)
    {
       if ((int)MCbank[i*(bankIn.dim+1)+bankIn.dim])
       {
	 for (j=0; j<bankIn.dim; j++) fprintf(stdout, "%e ", MCbank[i*(bankIn.dim+1)+j]);
         fprintf(stdout, "\n");
       }
     }
   }
  /* Convert output data structure. */

  bank = (SnglInspiralTable *) LALCalloc(1, sizeof(SnglInspiralTable));
  if (bank == NULL){
	  ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
  }

  *tiles = bank;

  for( k = 0; k < bankIn.nIni; k++ )
  {
	  if ((int)MCbank[k*(bankIn.dim+1)+bankIn.dim])
	  {
		  bank = bank->next = (SnglInspiralTable *) LALCalloc( 1, sizeof( SnglInspiralTable ) );
		  if (bank == NULL)
		  {
			  ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
		  }
		  bank->psi0 = MCbank[k*(bankIn.dim+1)];
		  bank->psi3 = MCbank[k*(bankIn.dim+1)+1];
		  bank->beta = MCbank[k*(bankIn.dim+1)+2];
	  }
  }

  /* Free tiles template, which is blank. */

  *ntiles = bankIn.nFin;
  bank = (*tiles)->next;
  LALFree( *tiles );
  free(bankIn.pMax);
  free(bankIn.pMin);
  free(bankIn.p);

  *tiles = bank;

  if (newSn.length!=0)
  {
	LALFree(newSn.data);
  }

  DETATCHSTATUSPTR(status);
  RETURN (status);
}

void MonteCarloBank(
		MetricFunc metric,
		MCBankIn *in,
		float *bank)
{
	int i, j;

	in->error = 1;

	/* Call the function to create a random bank */
	CreateRandomBank(in, bank);
	if (in->verbose)
	{
		for (i=0; i<in->nIni; i++)
		{
			for (j=0; j<in->dim; j++)
				fprintf(stdout, "%e ", bank[i*(in->dim+1)+j]);
			fprintf(stdout, "\n");
		}
		fprintf(stdout, "nIni=%d\n", in->nIni);
	}

	/* Throw away unwanted templates */
	FilterRandomBank(metric, in, bank);

	/* If everything works well return with error=0 */
	in->error = 0;
	return;
}

/*
void EstimateNumberOfTemplates(
		MetricFunc metric,
		MCBankIn  *in)
{
	return;
}
*/


void CreateRandomBank(
		MCBankIn *in,
		float *bank)
{
	int i, j, k;
	float e;

	srand(in->seed);
	for (i=0; i<in->nIni; i++)
	{
		k = i*(in->dim+1);
		for (j=0; j<in->dim; j++)
		{
			e = (float) rand()/(float) RAND_MAX;
			bank[k+j] = in->pMin[j] + e*(in->pMax[j] - in->pMin[j]);
		}
		/* We shall choose the last element in the series for a tag
		 * to tell us whether or not the current template is rejected;
		 * To begin with every template is a member of the bank
		 */
		bank[k+in->dim] = 1;
	}
	return;
}

void FilterRandomBank(
		MetricFunc metric,
		MCBankIn *in,
		float *bank)
{
	int i, j, m, dimp, dimpnIni, diag, sze;
	float dist, *x=NULL, *y=NULL, *gij=NULL;

	dimp = in->dim + 1;
	dimpnIni = dimp * in->nIni;
	in->nFin = in->nIni;

	sze = (int) (in->dim * (in->dim + 1))/2;
	metric(*in, &diag, gij);
	x = (float *) malloc(sizeof(float) * in->dim);
	y = (float *) malloc(sizeof(float) * in->dim);

	if (diag)
		gij = (float *) malloc(sizeof(float) * in->dim);
	else
		gij = (float *) malloc(sizeof(float) * sze);
	/* outer loop over all seed points */
	for (i=0; i<dimpnIni; i+=dimp)
	{
		if (bank[i+in->dim])
		{
		/* coordinates of first point */
		for (j=0; j<in->dim; j++) x[j] = in->p[j] = bank[i+j];
		/* compute the metric at the point (in->p) in question */
		metric(*in, &diag, gij);
		/* inner loop over test points */
		for (m=i+dimp; m<dimpnIni; m+=dimp)
		{
			/* compute dist to current pt only if it is not already rejected */
			if (bank[m+in->dim])
			{
				/* coordinates of m-th point */
				for (j=0; j<in->dim; j++)
				{
					y[j] = bank[m+j];
					if (in->verbose) fprintf(stdout, "%e %e %e ", x[j],y[j],gij[j]);
				}
				MCComputeDistance(in->dim, diag, x, y, gij, &dist);
				if (in->verbose) fprintf(stdout, "%e %d %d\n", dist, i, m);

				if (dist<in->MM)
				{
					bank[m+in->dim] = 0;
					--(in->nFin);
				}
			}
		}
		}
	}
	free(x);
	free(y);
	x = NULL;
	y = NULL;
	if (in->verbose) fprintf(stderr, "NoFinTmplts = %d\n", in->nFin);
	free(gij);
	gij = NULL;
	return;
}

void MCComputeDistance(
		int   dim,
		int   diag,
		float *x,
		float *y,
		float *gij,
		float *dist
		)
{
	int i, j, k;
	float dxi, dxj;
	*dist = 0.;
        k = 0;
	if (diag)
	{
		for (i=0; i<dim; i++)
		{
			dxi = x[i]-y[i];
			*dist += gij[i] * dxi * dxi;
		}
	}
	else
		for (i=0; i<dim; i++)
		{
			dxi = x[i]-y[i];
			*dist += gij[k] * dxi * dxi;
			/* fprintf(stdout,"%d %d \n", i, k); */
			for (j=i+1; j<dim; j++)
			{
		                /* fprintf(stdout,"%d %d %d %d \n", i, j, k, k+j-i); */
				dxj = x[j]-y[j];
				*dist += 2. * gij[k+j-i] * dxi *dxj;
			}
		        k += dim - i;
                }
	return;
}

/*
void BankStats(
		MetricFunc metric,
		MCBankIn *in,
		float *bank)
{
	return;
}
*/

void SchwarzschildMetric(
		MCBankIn in,
		int *diag,
		float *gij)
{
	*diag = 1;
	if (gij==NULL) return;
	gij[0] = 1.-2./in.p[0];
	gij[1] = in.p[0]*in.p[0];
	gij[2] = in.p[0]*in.p[0] * pow(sin(in.p[1]),2.);
	return;
}

void TwoSphereMetric(
		MCBankIn in,
		int *diag,
		float *gij)
{
	*diag = 1;
	if (gij==NULL) return;
	gij[0] = 1;
	gij[1] = pow(sin(in.p[0]),2.);
	return;
}

void FlatSpaceMetric(
		MCBankIn in,
		int *diag,
		float *gij)
{
	int i;
	*diag = 1;
	if (gij==NULL) return;
	for (i=0; i<in.dim; i++) gij[i] = 1;
	return;
}

void TwoDFlatSpacePolarMetric(
		MCBankIn in,
		int *diag,
		float *gij)
{
	*diag = 1;
	if (gij==NULL) return;
	gij[0] = 1;
	gij[1] = in.p[0]*in.p[0];
	return;
}

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

void BCVSpinMetric(
		MCBankIn in,
		int *diag,
		float *gij)
{
	int dbg=0, dim=4, N, i;
	double MinMatch, fmin, fmax, beta, **bcv2metric=NULL;

	*diag = 0;
	if (gij==NULL) return;

	if ( (bcv2metric = (double **) malloc(sizeof(double*)*dim)) == NULL )
	{
		bcv2metric = NULL;
		return;
	}

	for (i=0; i<dim; i++)
	{
		if ( (bcv2metric[i] = (double *) malloc(sizeof(double)*dim)) == NULL )
		{
			bcv2metric = NULL;
			return;
		}
	}

	N = (int)in.Npsd-1;
	beta = in.p[2];
	MinMatch = 1.- (double) in.MM;
	fmax = 1000.;
	fmin = 40.;

	BCVspin_metric(MinMatch, N, in.Sn, fmin, fmax, beta, bcv2metric, dbg);

	gij[0] = (float) bcv2metric[1][1];
	gij[1] = (float) bcv2metric[1][2];
	gij[2] = (float) bcv2metric[1][3];
	gij[3] = (float) bcv2metric[2][2];
	gij[4] = (float) bcv2metric[2][3];
	gij[5] = (float) bcv2metric[3][3];

	for (i=0; i<dim; i++)
	{
		free(bcv2metric[i]);
	}
	free(bcv2metric);
	return;
}

/* determinant of 3x3 matrix a[1..3][1..3] */

static double determinant3BCVSpinRandomBank(gsl_matrix *matrix)
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


int matrix3_determinant_plus(gsl_matrix *matrix,gsl_vector *eig)
{
  int i,j;
  double det;


  det=determinant3BCVSpinRandomBank(matrix);

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
			/*output*/ double **fit_point)
{
  const gsl_rng_type * T;
  gsl_rng * r;
  int i,j,k;
  double x[4],xd[4],sum,alpha[7],metric3[4][4];
  gsl_matrix *a = gsl_matrix_alloc(3, 3);
  gsl_vector *eig = gsl_vector_alloc(3);
  gsl_matrix *V = gsl_matrix_alloc(3, 3);
  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (3);

  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);


  sum=0;
  for(i=1;i<=6;i++)
  {
	  alpha[i]=gsl_rng_uniform_pos (r)-0.5;
	  sum+=alpha[i]*alpha[i];
  }
  sum=sqrt(sum);
  for(i=1;i<=6;i++)
  {
	  alpha[i]=alpha[i]/sum;
  }
  three_metric(funcG,alpha,metric3);


  for(i=1;i<=3;i++)
	  for(j=1;j<=3;j++)
		  gsl_matrix_set (a, i-1, j-1,metric3[i][j]);

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

int generate_metric_data(/* input */double MinMatch,
			 double funcG[7][7][4][4])
{

  extern double cont_data[JJ+1][4];
  int i,k,i2,j2;
  double x[4],maxr;
  double alpha[7],metric3[N_RANDOM+1][4][4],sum,rr[N_RANDOM+1],distance;
  double norm;
  const gsl_rng_type * T;
  gsl_rng * r;
  /*
  double fit_point[JJ+1][4];
  */
  double **fit_point;



  fit_point = (double **) malloc(sizeof(double *)*(JJ+1));
  for (i=0; i<4; i++)
  {
	  fit_point[i] = (double *) malloc(sizeof(double) * 4);
  }

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


#define TOL 1.0e-10

void svdfit_d_test(double x[], double y[], double sig[], int ndata, gsl_vector *a, int ma,
	gsl_matrix *u, gsl_matrix *v, gsl_vector *w, double *chisq,
	void (*funcs)(double, double []))
{
	int j,i;
	double wmax,tmp,thresh,sum,*afunc;
	gsl_vector *b= gsl_vector_alloc (ndata);
	gsl_vector *ws= gsl_vector_alloc (ma);


	if ( (afunc = (double *)malloc(sizeof(double) * (ma+1))) == NULL)
	{
		afunc = NULL;
		gsl_vector_free(b);
		gsl_vector_free(ws);
		return;
	}
	for (i=1;i<=ndata;i++) {
		(*funcs)(x[i],afunc);
		tmp=1.0/sig[i];
		for (j=1;j<=ma;j++) {
			gsl_matrix_set(u,i-1,j-1,afunc[j]*tmp);
		}
		gsl_vector_set(b,i-1,y[i]*tmp);
	}

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
	gsl_vector_free(ws);
	free(afunc);
}
#undef TOL


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
  int i, j, k, ma=9, *ia;
  double *x, *y, *sig;
  double chisq;

  gsl_matrix *u = gsl_matrix_alloc (ndata,ma);
  gsl_matrix *v = gsl_matrix_alloc (ma,ma);
  gsl_vector *w = gsl_vector_alloc (ma);
  gsl_vector *a= gsl_vector_alloc (ma);

  if ( ( ia = (int *)    malloc(sizeof(int)*(ma+1)))  == NULL ||
       (  x = (double *) malloc(sizeof(double)*(ndata+1))) == NULL ||
       (  y = (double *) malloc(sizeof(double)*(ndata+1)))  == NULL ||
       (sig = (double *) malloc(sizeof(double)*(ndata+1)))  == NULL
     )
  {
	  ia = NULL;
	  x = y = sig = NULL;
	  gsl_vector_free(w);
	  gsl_vector_free(a);
	  gsl_matrix_free(u);
	  gsl_matrix_free(v);
	  return 1;
  }

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

  free(ia);
  free(x);
  free(y);
  free(sig);
  gsl_vector_free(w);
  gsl_vector_free(a);
  gsl_matrix_free(u);
  gsl_matrix_free(v);

  return 0;
}
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
		double **bcv2metric,int dbg)
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

#undef PI
#undef DTRENORM
#undef N_RANDOM
#undef JJ

#endif


/* 2006.3.3  H. Takahashi

   The example program to calculate the value of beta
   of the template space with given beta_min and beta_max .
   This program also calculates the 2dim (\psi_0, \Psi_{3/2}) effective metric
   at the beta value which is calculated above.

*/

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

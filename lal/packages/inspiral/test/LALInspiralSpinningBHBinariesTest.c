#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "LALInspiral.h"
#include "LALNoiseModels.h"
#include "RealFFT.h"
#include "AVFactories.h"
#include "utils.h"

void jacobi(double **a, int n, double d[], double **v, int *nrot);
void spinLessWave(LALStatus *status, REAL4Vector *h, double *pars, int nPars);
void LALInspiralWaveWrapper ( LALStatus *status, REAL4Vector *signal, InspiralTemplate params);

INT4 lalDebugLevel=1;
NRCSID (LALINSPIRALSPINMODULATEDWAVEWRAPPERC, "$Id$");

int main()
{

	int nPars=5, length, i, j, k, nrot, label=0;
	static LALStatus status;
	static REAL8Vector psd;
	static REAL4Vector signal;
	static InspiralTemplate template;
	static double *infoMatrix, *covarMatrix, *pars, df, x;
        double eigenValues[nPars], *eigenVectors[nPars], *coMatrix[nPars], tiny =1.e-2;
	static REAL8 dTheta, dPhi;
	char params[nPars];
	FILE *fileout;

	fileout = fopen("pcLIGO.out", "a");
	pars = (double*) malloc(sizeof(double) * (nPars+1));
	infoMatrix = (double*) malloc(sizeof(double) * nPars * nPars);
	covarMatrix = (double*) malloc(sizeof(double) * nPars * nPars);
	for(i=0; i<nPars; i++) (coMatrix[i]) = (double*) malloc(sizeof(double) * nPars);
	for(i=0; i<nPars; i++) (eigenVectors[i]) = (double*) malloc(sizeof(double) * nPars);

	params[0] = 'D';
	params[1] = 'M';
	params[2] = 'h';

	pars[0] = 1.;

	/*
	pars[1] = 11.4;
	pars[2] = 14./(11.4*11.4);
	*/
   
	pars[1] = 3.0;
	pars[2] = 0.2499;
	pars[3] = 0.01;
	pars[4] = 1.5;

	fprintf(stderr, "\\begin{table}\n");
	/*
	pars[8] = LAL_PI/8.L;
	pars[9] = LAL_PI;
	*/

	signal.length = 0;
	spinLessWave(&status, &signal, pars, nPars);
	length = signal.length;

	template.tSampling = 2048.;
	df = template.tSampling/(float) length;
	psd.length = length/2+1;
	psd.data = (double*) malloc(sizeof(double)*psd.length);
	LALNoiseSpectralDensity (&status, &psd, &LALLIGOIPsd, df);
	LALInspiralCovarianceMatrix(&status, infoMatrix, covarMatrix, spinLessWave, psd, pars, nPars);
   
	for (i=0; i<nPars; i++) for (j=0; j<nPars; j++) coMatrix[i][j] = covarMatrix[i*nPars+j];
	for (i=0; i<nPars; i++) {for (j=0; j<nPars; j++) fprintf(stderr, "%e\t", coMatrix[i][j]); fprintf(stderr, "\n");}
	fprintf(stderr, "\n");
	jacobi(coMatrix, nPars, eigenValues, eigenVectors, &nrot);
	for (i=0; i<nPars; i++) {for (j=0; j<nPars; j++) fprintf(stderr, "%e\t",eigenVectors[i][j]); fprintf(stderr, "\n");}
	fprintf(stderr, "Principal components are: $(");
	for(i=0; i<nPars; i++) fprintf(stderr, "%10.4f,\\ ", eigenValues[i]);
	for(i=0; i<nPars; i++) fprintf(fileout, "%10.4f\t", eigenValues[i]);
	fprintf(fileout, "\n}");
	fprintf(stderr, ")$\n}");
	fprintf(stderr,"\\label{label%d}\n", label++);
	fprintf(stderr, "\\begin{tabular}{rrrr}\n");
	fprintf(stderr, "\\hline\n");
	for(i=0; i<nPars; i++)
	{
		for(j=0; j<=i; j++)
	
		{
			/*
			coMatrix[i][j] = coMatrix[j][i] = covarMatrix[i*nPars+j];
			 */
			if (i==j) 
			{
			
				x = sqrt(covarMatrix[i*nPars+j]);
				if (i==0 || i==2) x = x/pars[i];
			}
			else 
			{
				x = covarMatrix[i*nPars+j]/sqrt(covarMatrix[i*nPars+i]*covarMatrix[j*nPars+j]);
			}
			/* if (fabs(x)>0.5) fprintf(stderr, "<%c,%c>=%e\n", params[i], params[j], x); */
			if (fabs(x) > 0.5 && i!=j)
				fprintf(stderr, "& $ \\mathbf{%7.4f} $ ", x);
			else
				fprintf(stderr, "& $ %7.4f $ ", x);
		}
		fprintf(stderr, "\\\\ \n");
	}
			
	fprintf(stderr, "\n");
	fprintf(stderr, "\\hline\n");
	fprintf(stderr, "\\end{tabular}\n");
	fprintf(stderr, "\\end{table}\n");
	/*
	for(i=0; i<nPars; i++)
	{
		for(j=0; j<nPars; j++)
	
		{
			x=0; for(k=0; k<nPars; k++) x+= covarMatrix[i*nPars+k] * infoMatrix[k*nPars+j];
			fprintf(stderr, "%e\t", x);
		}
		fprintf(stderr, "\n");
	}
	*/
	free(psd.data);
	free(infoMatrix);
	free(covarMatrix);
	fclose(fileout);
}

void spinLessWave(LALStatus *status, REAL4Vector *h, double *pars, int nPars)
{

   float preamp;
   InspiralTemplate params;

   INITSTATUS (status, "LALInspiralSpinModulatedWaveWrapper", LALINSPIRALSPINMODULATEDWAVEWRAPPERC);
   ATTATCHSTATUSPTR(status);

   preamp = pars[0]*pars[nPars];

   /*
   params.startTime=0.0; 
   params.startPhase=0.0;
   params.totalMass=pars[1]/pow(pars[2],0.6); 
   params.t0=pars[1];
   params.t3=pars[2];
   */
   params.totalMass= exp(pars[1]);
   params.eta=pars[2];
   params.startTime=pars[3]; 
   params.startPhase=pars[4];
   params.massChoice=t03;
   params.massChoice=totalMassAndEta;

   params.OmegaS = 0.;
   params.Theta = 0.;
   params.ieta=1; 
   params.fLower=40.; 
   params.fCutoff=400.;
   params.tSampling=2048.;
   params.signalAmplitude=1.0;
   params.nStartPad=0;
   params.nEndPad=0;
   params.order=4;
   params.approximant=TaylorT3;
   params.distance = 1.e9 * LAL_PC_SI/LAL_C_SI;
   LALInspiralParameterCalc(status->statusPtr, &params);

   params.orbitTheta0 = 0;
   params.orbitPhi0 = 0;
   params.sourceTheta = 0;
   params.sourcePhi = 0;

   params.spin1[0] = 0;
   params.spin1[1] = 0;
   params.spin1[2] = 0;

   params.spin2[0] = 0;
   params.spin2[1] = 0;
   params.spin2[2] = 0;

   CHECKSTATUSPTR(status);
	   
   if (h->length==0)
   {
   
	   LALInspiralWaveWrapper(status->statusPtr, h, params);
	   CHECKSTATUSPTR(status);
	   fprintf(stderr, "\\caption{\n$n=%d,$\\ $(m_1, m_2)=(%5.2f, %5.2f),$\\ $n_{\\rm cyc}=%5.0f.$ \n", 
			   h->length, params.mass1, params.mass2, 1.6*params.fLower*params.t0);
   }
   else
   {
	   int i;
	   LALInspiralWaveWrapper(status->statusPtr, h, params);
	   CHECKSTATUSPTR(status);
	   for (i=0; i<h->length; i++) h->data[i] *= preamp;
   }
   DETATCHSTATUSPTR(status);
   RETURN (status); 
	
}


void 
LALInspiralWaveWrapper
(
 LALStatus *status, 
 REAL4Vector *signal, 
 InspiralTemplate params)
{

   REAL8 dt, t0=0.0;
   UINT4 n;

   INITSTATUS (status, "LALInspiralSpinModulatedWaveWrapper", LALINSPIRALSPINMODULATEDWAVEWRAPPERC);
   ATTATCHSTATUSPTR(status);

   dt = 1./params.tSampling;
   if (signal->length==0) 
   {
	   params.nEndPad=1000;
	   LALInspiralWaveLength(status->statusPtr, &n, params);
	   CHECKSTATUSPTR(status);
	   signal->length = 2*n;
   }
   else
   {
	   LALInspiralParameterCalc(status->statusPtr, &params);
	   CHECKSTATUSPTR(status);
	   LALInspiralWave(status->statusPtr, signal, &params);
	   CHECKSTATUSPTR(status);
   
   /*
	   printf_timeseries(signal->length, signal->data, dt, t0, params);
	   exit(0);
   */
   }
   DETATCHSTATUSPTR(status);
   RETURN (status); 
}

void LALInspiralCovarianceMatrix
(
 LALStatus *status,
 double *infoMatrix, 
 double *covarMatrix,
 void   (*wavefunc)(LALStatus *status, REAL4Vector *h, double *par, int nPars), 
 REAL8Vector psd,
 double *par,
 int    nPars)
{

	int k,m, LENGTH;
	REAL4Vector *hk=NULL, *hm=NULL, *buff=NULL;
        double tmpAmp, *buffMatrix, c0;
	int *workSpace;
	RealFFTPlan *revp=NULL, *fwdp=NULL;
   
	INITSTATUS (status, "LALInspiralSpinModulatedWaveWrapper", LALINSPIRALSPINMODULATEDWAVEWRAPPERC);
	ATTATCHSTATUSPTR(status);
   
	LENGTH = (psd.length-1)*2;

	buffMatrix = (double*) malloc(sizeof(double) * nPars * nPars);
	workSpace = (int*) malloc(sizeof(int) * nPars);

	LALCreateVector(status->statusPtr, &buff, LENGTH);
	CHECKSTATUSPTR(status);
	LALCreateVector(status->statusPtr, &hm, LENGTH);
	CHECKSTATUSPTR(status);
	LALCreateVector(status->statusPtr, &hk, LENGTH);
	CHECKSTATUSPTR(status);
	LALCreateReverseRealFFTPlan (status->statusPtr, &revp, LENGTH, 0);
	CHECKSTATUSPTR(status);
	LALCreateForwardRealFFTPlan (status->statusPtr, &fwdp, LENGTH, 0);
	CHECKSTATUSPTR(status);
   
	/*
	for (k=0; k<LENGTH; k++)
	{
		printf("%e\n", hm.data[k]);
	}
	printf("&\n");
	*/

	par[nPars] = 1.;
	tmpAmp = par[0];
	par[0] = 1.;
	(*wavefunc)(status->statusPtr, buff, par, nPars);
	CHECKSTATUSPTR(status);
	LALREAL4VectorFFT(status->statusPtr, hm, buff, fwdp);
	CHECKSTATUSPTR(status);
	LALInspiralWaveNormalise(status->statusPtr, hm, &c0, psd);
	CHECKSTATUSPTR(status);
	/*
	fprintf(stderr, "norm=%e\n", c0); 
	*/

	par[nPars] = 1./c0;
	par[0] = tmpAmp;
	(*wavefunc)(status->statusPtr, buff, par, nPars);
	LALREAL4VectorFFT(status->statusPtr, hm, buff, fwdp);
	CHECKSTATUSPTR(status);
	LALInspiralWaveNormalise(status->statusPtr, hm, &c0, psd);
	CHECKSTATUSPTR(status);
	/*
	fprintf(stderr, "norm=%e\n", c0); 
	*/

	for (k=0; k<nPars; k++)
	{
		LALDerivFunc(status->statusPtr, hk, wavefunc, par, k, nPars, fwdp);
		CHECKSTATUSPTR(status);
		LALMakeCopyFloat(hk->data, hm->data, hk->length);
		LALBareOverlap(status->statusPtr, &c0, hk, hm, revp, &psd);
		infoMatrix[k*nPars+k] = c0;
		CHECKSTATUSPTR(status);
		buffMatrix[k*nPars+k] = infoMatrix[k*nPars+k];
		/* fprintf(stderr, "%e ", c0); */
		for (m=k+1; m<nPars; m++)
		{	
			LALDerivFunc(status->statusPtr, hm, wavefunc, par, m, nPars, fwdp);
			CHECKSTATUSPTR(status);
			LALBareOverlap(status->statusPtr, &c0, hk, hm, revp, &psd);
			infoMatrix[k*nPars+m] = c0;
			CHECKSTATUSPTR(status);
			buffMatrix[k*nPars+m] = infoMatrix[k*nPars+m];
			buffMatrix[m*nPars+k] = infoMatrix[m*nPars+k] = infoMatrix[k*nPars+m];
			/* fprintf(stderr, "%e ", c0); */
		}
		/* fprintf(stderr, "\n"); */
	}

	matrixInverse(&nPars, buffMatrix, covarMatrix, workSpace);
	
	LALDestroyVector (status->statusPtr, &hm);
	CHECKSTATUSPTR(status);
	LALDestroyVector (status->statusPtr, &hk);
	CHECKSTATUSPTR(status);
	LALDestroyVector (status->statusPtr, &buff);
	CHECKSTATUSPTR(status);
	LALDestroyRealFFTPlan (status->statusPtr, &fwdp);
	CHECKSTATUSPTR(status);
	LALDestroyRealFFTPlan (status->statusPtr, &revp);
	CHECKSTATUSPTR(status);
	free(buffMatrix);
	free(workSpace);
	DETATCHSTATUSPTR(status);
	RETURN (status); 
}

void LALDerivFunc
(
 LALStatus *status,
 REAL4Vector *h, 
 void        (*wavefunc)(LALStatus *status, REAL4Vector *h, double *par, int nPars),
 double      *par, 
 int         k,
 int         nPars,
 RealFFTPlan *fwdp)
{
	REAL4Vector *hp=NULL, *hm=NULL, *buff=NULL;
        double	epsilon, fac, *tempPar;
	int i, LENGTH;

     	INITSTATUS (status, "LALInspiralSpinModulatedWaveWrapper", LALINSPIRALSPINMODULATEDWAVEWRAPPERC);
	ATTATCHSTATUSPTR(status);
	epsilon = 1.e-6;
	LENGTH = h->length;

	LALCreateVector(status->statusPtr, &buff, LENGTH);
	CHECKSTATUSPTR(status);
	LALCreateVector(status->statusPtr, &hp, LENGTH);
	CHECKSTATUSPTR(status);
	LALCreateVector(status->statusPtr, &hm, LENGTH);
	CHECKSTATUSPTR(status);
	tempPar = (double*) malloc(sizeof(double) * (nPars+1));
	
	LALMakeCopyDouble(par, tempPar, nPars);
	tempPar[nPars]=par[nPars];

	tempPar[k] = par[k] * (1.0L + epsilon);
	(*wavefunc)(status->statusPtr, buff, tempPar, nPars);
	LALREAL4VectorFFT(status->statusPtr, hp, buff, fwdp);
	CHECKSTATUSPTR(status);

	tempPar[k] = par[k] * (1.0L - epsilon);
	(*wavefunc)(status->statusPtr, buff, tempPar, nPars);
	LALREAL4VectorFFT(status->statusPtr, hm, buff, fwdp);
	CHECKSTATUSPTR(status);
		
	fac = 2.0L * epsilon * par[k];
	for (i=0; i<LENGTH; i++)
	{
		h->data[i] = hp->data[i]/fac - hm->data[i]/fac;
	}
	LALDestroyVector(status->statusPtr, &buff);
	CHECKSTATUSPTR(status);
	LALDestroyVector(status->statusPtr, &hp);
	CHECKSTATUSPTR(status);
	LALDestroyVector(status->statusPtr, &hm);
	CHECKSTATUSPTR(status);
	free(tempPar);
	DETATCHSTATUSPTR(status);
	RETURN (status); 
}


void LALMakeCopyFloat(float *a, float *b, int n)
{
	int i;

	for (i=0; i<n; i++)
	{
		b[i] = a[i];
	}
}

void LALMakeCopyDouble(double *a, double *b, int n)
{
	int i;

	for (i=0; i<n; i++)
	{
		b[i] = a[i];
	}
}

void LALBareOverlap 
(
 LALStatus *status,
 REAL8       *x,
 REAL4Vector *s1,  
 REAL4Vector *s2, 
 RealFFTPlan *revp,
 REAL8Vector *psd
) 
{

   REAL4Vector *buff=NULL;
   int i;
   InspiralWaveCorrelateIn corrin;
   /*
    * int i;
    */

   INITSTATUS (status, "LALInspiralSpinModulatedWaveWrapper", LALINSPIRALSPINMODULATEDWAVEWRAPPERC);
   ATTATCHSTATUSPTR(status);

   LALCreateVector(status->statusPtr, &buff, s1->length);
   CHECKSTATUSPTR(status);

   corrin.signal1 = *s1;
   corrin.signal2 = *s2;
   corrin.psd = *psd;
   corrin.revp = revp;

   LALInspiralWaveCorrelate (status->statusPtr, buff, corrin);
   CHECKSTATUSPTR(status);
   /*
   for (i=0; i<buff->length; i++) printf("%e\n", buff->data[i]*buff->length);
   exit(0);
   */
   *x = buff->data[0]*buff->length;

   LALDestroyVector(status->statusPtr, &buff);
   CHECKSTATUSPTR(status);
   DETATCHSTATUSPTR(status);
   RETURN (status); 

}
/* inverse.f -- translated by f2c (version 20000531).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
#include "f2c.h"
*/

/* n is the dimension of the square matrix, a is the input matrix
* (which will be overwritten) b is the inverse and c is workspace
*/
int matrixInverse(int *n, double *a, double *b, int *c__)
{
    /* System generated locals */
    int a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    static double d__;
    static int i__, j;

/* -----------------------------------------------------------------------C */
    /* Parameter adjustments */
    --c__;
    b_dim1 = *n;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;
    a_dim1 = *n;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    b[i__ + j * b_dim1] = 0.;
	}
	b[i__ + i__ * b_dim1] = 1.;
    }
    if (ludcmp(&a[a_offset], n, n, &c__[1], &d__)) return 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	lubksb(&a[a_offset], n, n, &c__[1], &b[j * b_dim1 + 1]);
    }
    return 0;
} /* inverse_ */


/* -----------------------------------------------------------------------C */
/* Subroutine */ int lubksb(double *a, int *n, int *np, int *indx, double *b)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static int i__, j, ii, ll;
    static double sum;

/* -----------------------------------------------------------------------C */
    /* Parameter adjustments */
    --b;
    --indx;
    a_dim1 = *np;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;

    /* Function Body */
    ii = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ll = indx[i__];
	sum = b[ll];
	b[ll] = b[i__];
	if (ii != 0) {
	    i__2 = i__ - 1;
	    for (j = ii; j <= i__2; ++j) {
		sum -= a[i__ + j * a_dim1] * b[j];
/* L11: */
	    }
	} else if (sum != (float)0.) {
	    ii = i__;
	}
	b[i__] = sum;
/* L12: */
    }
    for (i__ = *n; i__ >= 1; --i__) {
	sum = b[i__];
	i__1 = *n;
	for (j = i__ + 1; j <= i__1; ++j) {
	    sum -= a[i__ + j * a_dim1] * b[j];
/* L13: */
	}
	b[i__] = sum / a[i__ + i__ * a_dim1];
/* L14: */
    }
    return 0;
} /* lubksb_ */

/* ludcmp.f -- translated by f2c (version 20000531).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
#include "f2c.h"
*/
/* -----------------------------------------------------------------------C */
/* Subroutine */ int ludcmp(double *a, int *n, int *np, int *indx, double *d__)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2, i__3;
    double d__1, d__2;

    /* Builtin functions */

    /* Local variables */
    static int imax, i__, j, k;
    static double aamax, vv[500], dum, sum;

/* -----------------------------------------------------------------------C */
    /* Parameter adjustments */
    --indx;
    a_dim1 = *np;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;

    /* Function Body */
    *d__ = (float)1.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	aamax = (float)0.;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    if ((d__1 = a[i__ + j * a_dim1], abs(d__1)) > aamax) {
		aamax = (d__2 = a[i__ + j * a_dim1], abs(d__2));
	    }
/* L11: */
	}
	if (aamax == (float)0.) {
	    /*
	    fprintf(stderr, "\n singular matrix in ludcmp \n");
	    */
            return(1);
	}
	vv[i__ - 1] = (float)1. / aamax;
/* L12: */
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    sum = a[i__ + j * a_dim1];
	    i__3 = i__ - 1;
	    for (k = 1; k <= i__3; ++k) {
		sum -= a[i__ + k * a_dim1] * a[k + j * a_dim1];
/* L13: */
	    }
	    a[i__ + j * a_dim1] = sum;
/* L14: */
	}
	aamax = (float)0.;
	i__2 = *n;
	for (i__ = j; i__ <= i__2; ++i__) {
	    sum = a[i__ + j * a_dim1];
	    i__3 = j - 1;
	    for (k = 1; k <= i__3; ++k) {
		sum -= a[i__ + k * a_dim1] * a[k + j * a_dim1];
/* L15: */
	    }
	    a[i__ + j * a_dim1] = sum;
	    dum = vv[i__ - 1] * abs(sum);
	    if (dum >= aamax) {
		imax = i__;
		aamax = dum;
	    }
/* L16: */
	}
	if (j != imax) {
	    i__2 = *n;
	    for (k = 1; k <= i__2; ++k) {
		dum = a[imax + k * a_dim1];
		a[imax + k * a_dim1] = a[j + k * a_dim1];
		a[j + k * a_dim1] = dum;
/* L17: */
	    }
	    *d__ = -(*d__);
	    vv[imax - 1] = vv[j - 1];
	}
	indx[j] = imax;
	if (a[j + j * a_dim1] == (float)0.) {
	    a[j + j * a_dim1] = 1e-20;
	}
	if (j != *n) {
	    dum = (float)1. / a[j + j * a_dim1];
	    i__2 = *n;
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
		a[i__ + j * a_dim1] *= dum;
/* L18: */
	    }
	}
/* L19: */
    }
    return 0;
} /* ludcmp_ */

/*
*/
void printf_timeseries (int n, float *signal, double delta, double t0, InspiralTemplate par) 
{
  int i=0;
  FILE *outfile1;

  outfile1=fopen("wave1.dat","a");
  do 
  {
	  float x=0.;
	  /*
	  x = 0.;
	  if (i>=4096 && i<=4096+2048) x = signal[0] * 2 * cos(2.*LAL_PI*400.*(i-4096)*delta);
	  if (i==8192-10) x = signal[0]*15.;
	  if (i==8192-9) x = signal[0]*15./2.;
	  if (i==8192-11) x = signal[0]*15./2.;
	  */
	  fprintf (outfile1,"%e %e\n", i*delta+t0, *(signal+i) + x);
  }
  while (n-++i); 

  fprintf(outfile1,"&\n");
  fclose(outfile1);
}

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);

void jacobi(double **a, int n, double d[], double **v, int *nrot)
{
	int j,iq,ip,i;
	float tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

	b=(float*) malloc(sizeof(float)*n);
	z=(float*) malloc(sizeof(float)*n);

	for (ip=0;ip<n;ip++) {
		for (iq=0;iq<n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=0;ip<n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	*nrot=0;
	for (i=0;i<50;i++) {
		sm=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0) {
			free(z);
			free(b);
			return;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip])
					&& (float)(fabs(d[iq])+g) == (float)fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((float)(fabs(h)+g) == (float)fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=0;j<ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=1;j<n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}
		for (ip=0;ip<n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	fprintf(stderr, "Too many iterations in routine jacobi\n");
}
#undef ROTATE

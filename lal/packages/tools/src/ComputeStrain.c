/* <lalLaTeX>
\subsection{Module \texttt{ComputeStrain.c}}

\subsubsection*{Description}


\subsubsection*{Algorithm}


\subsubsection*{Uses}
\begin{verbatim}
None
\end{verbatim}

\subsubsection*{Notes}

</lalLaTeX> */

#include <math.h>
#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Calibration.h>
#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/AVFactories.h>
#include <lal/Window.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/RealFFT.h>
#include <lal/AVFactories.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#define MAXALPHAS 10000

#define N_FIR_LP_ALPHAS 8192
#define fhigh_FIRLP_ALPHAS 0.00001220703125

#define N_FIR_LP 100
#define fhigh_FIRLP .6103515625

#define N_FIR_HP 2000
#define fhigh_FIRHP 0.00244140625

NRCSID( COMPUTESTRAINC, "$Id$" );

void LALComputeStrain(
    LALStatus              *status,
    StrainOut              *output,    
    StrainIn               *input
    )
{
/* Inverse sensing, servo, analog actuation, digital x 
actuation  digital y actuation */
static REAL8TimeSeries uphR,hCX,hCY,ALPHAS,upALPHAS;
int p, NWingsR,NWingsC;
REAL8IIRFilter LPFIR;
REAL8IIRFilter HPFIR;
REAL8IIRFilter ALPHASLPFIR;
REAL8IIRFilter *CinvWings=NULL, *AAWings=NULL, *AXWings=NULL, *AYWings=NULL;


 INITSTATUS( status, "LALComputeStrain", COMPUTESTRAINC );
 ATTATCHSTATUSPTR( status );

 LALGetFactors(status->statusPtr, output, input);
 CHECKSTATUSPTR( status );
 
 LALMakeFIRLP(status->statusPtr, &LPFIR, input->CinvUSF);
 CHECKSTATUSPTR( status );

 LALMakeFIRHP(status->statusPtr, &HPFIR);
 CHECKSTATUSPTR( status );

 /* Create vectors that will hold the residual, control and net strain signals */
 LALDCreateVector(status->statusPtr,&hCX.data,input->AS_Q.data->length);
 CHECKSTATUSPTR( status );
 LALDCreateVector(status->statusPtr,&hCY.data,input->AS_Q.data->length);
 CHECKSTATUSPTR( status );

 NWingsR = (int)(input->wings/input->AS_Q.deltaT + 0.5) * input->CinvUSF;
 NWingsC = (int)(input->wings/input->AS_Q.deltaT + 0.5);

 /* copy AS_Q input into residual strain as double */  
 for (p=0; p<(int)output->hR.data->length; p++) 
   {
    output->hR.data->data[p]=input->DARM_ERR.data->data[p];
   }

 /* copy DARM_CTRL input into control strain as double */  
 for (p=0; p<(int)output->hC.data->length; p++) 
   {
     output->hC.data->data[p]=input->DARM.data->data[p];
   }

 /* unit impulse */
 if (input->delta)
   {
     for (p=0; p<(int)output->hR.data->length; p++) 
       {
	 output->hR.data->data[p]=0;
       }
     output->hR.data->data[output->hC.data->length/2]=1.0;
   }

 /* unit impulse */  
 if (input->delta)
   {
     for (p=0; p<(int)output->hC.data->length; p++) 
       {
	 output->hC.data->data[p]=0;
       }
     output->hC.data->data[output->hC.data->length/2]=1.0;
   }
 
 /* ---------- Compute Residual Strain -------------*/

 LALDCreateVector(status->statusPtr,&uphR.data,input->CinvUSF*input->AS_Q.data->length);
 CHECKSTATUSPTR( status );
 uphR.deltaT=input->AS_Q.deltaT/input->CinvUSF;

 /* then we upsample (and smooth it with a low pass filter) */
 if(XLALUpsample(&uphR, &(output->hR), input->CinvUSF)) 
   { 
     ABORT(status,117,"Broke upsampling hR");
   }

 /* apply delay (actually an advance) */ 
 for (p=0; p<(int)uphR.data->length+input->CinvDelay; p++){
   uphR.data->data[p]=uphR.data->data[p-input->CinvDelay];
 }

 /* An odd filtere with N points introduces an (N-1)/2 delay */
 /* apply advance to compensate for FIR delay */
 for (p=0; p<(int)uphR.data->length-(2*N_FIR_LP); p++){
   uphR.data->data[p]=uphR.data->data[p+(2*N_FIR_LP)];
 }
 /* Low pass filter twice to smooth time series */
 LALIIRFilterREAL8Vector(status->statusPtr,uphR.data,&(LPFIR));
 CHECKSTATUSPTR( status );
 LALIIRFilterREAL8Vector(status->statusPtr,uphR.data,&(LPFIR));
 CHECKSTATUSPTR( status );

 /* ===================== */
 /* CAREFUL FILTERING: FILTER UP TO WINGS, THEN COPY FILTERS THEN CONTINUE UNTIL END */
 /* filter only up until the wings */
 for(p = 0; p < input->NCinv; p++){
   int k;
   for(k = NWingsR/2; k < (int)uphR.data->length-3*NWingsR/2; k++){
     uphR.data->data[k]=LALDIIRFilter(uphR.data->data[k], &(input->Cinv[p]));
   }
 }
 /* Here what I need to record filter histories */
 LALCopyFilter(status->statusPtr, &CinvWings, input->Cinv, input->NCinv);
 CHECKSTATUSPTR( status );
 /* then we filter wings as well */
 for(p = 0; p < input->NCinv; p++){
   int k;
   for(k = uphR.data->length-3*NWingsR/2; k < (int)uphR.data->length-NWingsR/2; k++){
     uphR.data->data[k]=LALDIIRFilter(uphR.data->data[k], &CinvWings[p]);
   }
 }
 /* Then we need to destroy the filter */
 LALFreeFilter(status->statusPtr,CinvWings,input->NCinv);
 CHECKSTATUSPTR( status );
 /* ===================== */


 /* apply advance to compensate for Low Pass FIR delay */
 for (p=0; p<(int)uphR.data->length-(2*N_FIR_LP); p++){
   uphR.data->data[p]=uphR.data->data[p+(2*N_FIR_LP)];
 }
 /* Low pass filter twice to smooth time series */
 LALIIRFilterREAL8Vector(status->statusPtr,uphR.data,&(LPFIR));
 CHECKSTATUSPTR( status );
 LALIIRFilterREAL8Vector(status->statusPtr,uphR.data,&(LPFIR));
 CHECKSTATUSPTR( status );
 
 /* then we downsample and voila' */
 for (p=0; p<(int)output->hR.data->length; p++) {
   output->hR.data->data[p]=uphR.data->data[p*input->CinvUSF];
 }

 LALDDestroyVector(status->statusPtr,&uphR.data);
 CHECKSTATUSPTR( status );

 /* Create time series that hold alpha time series and upsampled alpha time-series */
 LALDCreateVector(status->statusPtr,&ALPHAS.data,output->alpha.data->length);
 CHECKSTATUSPTR( status );
 LALDCreateVector(status->statusPtr,&upALPHAS.data,input->AS_Q.data->length);
 CHECKSTATUSPTR( status );
 upALPHAS.deltaT=input->AS_Q.deltaT;
 
 /* copy factors into time series */
 for (p=0; p<(int)ALPHAS.data->length; p++) 
   {
     ALPHAS.data->data[p]=output->alphabeta.data->data[p].re;
   }
 
 /* upsample using a linear interpolation */
 if(XLALUpsampleLinear(&upALPHAS, &ALPHAS, (int) (output->alphabeta.deltaT/input->AS_Q.deltaT+0.5)))
   { 
     ABORT(status,117,"Broke upsampling Alphas");
   }
 
 /* Generate alphas LP filter to smooth linearly interpolated factors time series */
 LALMakeFIRLPALPHAS(status->statusPtr, &ALPHASLPFIR);
 CHECKSTATUSPTR( status );

 /* Note that I do not compensate for delay 
    if factors are computed every second then everything is fine */
 if (input->fftconv)
   {
     LALFIRFilter(status->statusPtr,&upALPHAS,&ALPHASLPFIR);
     CHECKSTATUSPTR( status );
   }else{
     LALIIRFilterREAL8Vector(status->statusPtr,upALPHAS.data,&(ALPHASLPFIR));
     CHECKSTATUSPTR( status );
   }
 
 /* finally we divide by alpha*beta */
 if (input->usefactors)
   {
     if(XLALDivideTimeSeries(&(output->hR), &upALPHAS)) 
       {
	 ABORT(status,116,"Broke at hR/alpha");
       }
   }

 if (input->outalphas)
   {
     /* set output to alphas */
     for (p=0; p<(int)output->hR.data->length; p++) {
       output->hR.data->data[p]=upALPHAS.data->data[p];
     }
   }

 /* destroy low and high pass filters */ 
 LALDDestroyVector(status->statusPtr,&(ALPHASLPFIR.directCoef));
 CHECKSTATUSPTR( status );
 LALDDestroyVector(status->statusPtr,&(ALPHASLPFIR.recursCoef));
 CHECKSTATUSPTR( status );
 LALDDestroyVector(status->statusPtr,&(ALPHASLPFIR.history));
 CHECKSTATUSPTR( status );
     
 /* destroy both alphas time series */
 LALDDestroyVector(status->statusPtr,&upALPHAS.data);
 CHECKSTATUSPTR( status );
 LALDDestroyVector(status->statusPtr,&ALPHAS.data);
 CHECKSTATUSPTR( status );

 /* ---------- Compute Control Strain -------------*/

 /* first implement the time delay filter; start at end to avoid overwriting */
 for (p = (int)output->hC.data->length-1; p >= input->AADelay; p--){
   output->hC.data->data[p]=output->hC.data->data[p-input->AADelay];
 }

 /* then apply advance to compensate for FIR delay */
 for (p=0; p<(int)output->hC.data->length-(2*N_FIR_HP); p++){
   output->hC.data->data[p]=output->hC.data->data[p+(2*N_FIR_HP)];
 }
 /* then high pass filter DARM_CTRL */
 if (input->fftconv)
   {
     LALFIRFilter(status->statusPtr,&(output->hC),&(HPFIR));
     CHECKSTATUSPTR( status );
     LALFIRFilter(status->statusPtr,&(output->hC),&(HPFIR));
     CHECKSTATUSPTR( status );
   }else{
     LALIIRFilterREAL8Vector(status->statusPtr,output->hC.data,&(HPFIR));
     CHECKSTATUSPTR( status );
     LALIIRFilterREAL8Vector(status->statusPtr,output->hC.data,&(HPFIR));
     CHECKSTATUSPTR( status );
   }

 /* ===================== */
 /* CAREFUL FILTERING: FILTER UP TO WINGS, THEN COPY FILTERS THEN CONTINUE UNTIL END */
 /* we filter but only up until the wings */
 for(p = 0; p < input->NAA; p++){
   int k;
   for(k = NWingsC/2; k < (int)output->hC.data->length-3*NWingsC/2; k++){
     output->hC.data->data[k]=LALDIIRFilter(output->hC.data->data[k], &(input->AA[p]));
   }
 }
 /* Here what I need to record filter histories */
 LALCopyFilter(status->statusPtr, &AAWings, input->AA, input->NAA);
 CHECKSTATUSPTR( status );
 /* then we filter wings as well */
 for(p = 0; p < input->NAA; p++){
   int k;
   for(k = output->hC.data->length-3*NWingsC/2; k < (int)output->hC.data->length-NWingsC/2; k++){
     output->hC.data->data[k]=LALDIIRFilter(output->hC.data->data[k], &AAWings[p]);
   }
 }
 /* Then we need to destroy the filter */
 LALFreeFilter(status->statusPtr,AAWings,input->NAA);
 CHECKSTATUSPTR( status );
 /* ===================== */

 /* Copy data into x and y time series for parallel filtering */
 for (p=0; p < (int)output->hC.data->length; p++){
   hCX.data->data[p] = output->hC.data->data[p];
   hCY.data->data[p] = output->hC.data->data[p];
 }

 /* Filter x-arm */
 /* ===================== */
 /* CAREFUL FILTERING: FILTER UP TO WINGS, THEN COPY FILTERS THEN CONTINUE UNTIL END */
 /* we filter only up until the wings */
 for(p = 0; p < input->NAX; p++){
   int k;
   for(k = NWingsC/2; k < (int)hCX.data->length-3*NWingsC/2; k++){
     hCX.data->data[k]=LALDIIRFilter(hCX.data->data[k], &(input->AX[p]));
   }
 }
 /* Here what I need to record filter histories */
 LALCopyFilter(status->statusPtr, &AXWings, input->AX, input->NAX);
 CHECKSTATUSPTR( status );
 /* then we filter wings as well */
 for(p = 0; p < input->NAX; p++){
   int k;
   for(k = hCX.data->length-3*NWingsC/2; k < (int)hCX.data->length-NWingsC/2; k++){
     hCX.data->data[k]=LALDIIRFilter(hCX.data->data[k], &AXWings[p]);
   }
 }
 /* Then we need to destroy the filter */
 LALFreeFilter(status->statusPtr,AXWings,input->NAX);
 CHECKSTATUSPTR( status );
 /* ===================== */


 /* Filter y-arm */
 /* ===================== */
 /* CAREFUL FILTERING: FILTER UP TO WINGS, THEN COPY FILTERS THEN CONTINUE UNTIL END */
 /* we filter only up until the wings */
 for(p = 0; p < input->NAY; p++){
   int k;
   for(k = NWingsC/2; k < (int)hCY.data->length-3*NWingsC/2; k++){
     hCY.data->data[k]=LALDIIRFilter(hCY.data->data[k], &(input->AY[p]));
   }
 }
 /* Here what I need to record filter histories */
 LALCopyFilter(status->statusPtr, &AYWings, input->AY, input->NAY);
 CHECKSTATUSPTR( status );
 /* then we filter wings as well */
 for(p = 0; p < input->NAY; p++){
   int k;
   for(k = hCY.data->length-3*NWingsC/2; k < (int)hCY.data->length-NWingsC/2; k++){
     hCY.data->data[k]=LALDIIRFilter(hCY.data->data[k], &AYWings[p]);
   }
 }
 /* Then we need to destroy the filter */
 LALFreeFilter(status->statusPtr,AYWings,input->NAY);
 CHECKSTATUSPTR( status );
 /* ===================== */

 /* add them together to make the total control signal */
 for (p=0; p < (INT4)output-> h.data->length; p++) {
    output->hC.data->data[p]=(hCX.data->data[p]+hCY.data->data[p])/2; 
 }


 /* ---------- Compute Net Strain -------------*/

 /* add x-arm and y-arm and voila' */
 for (p=0; p< (int)output->h.data->length; p++){
   output->h.data->data[p] = output->hC.data->data[p]+output->hR.data->data[p];
 }

 /* ------------------------------------------*/

 /* destroy vectors that hold the data */
 LALDDestroyVector(status->statusPtr,&hCX.data);
 CHECKSTATUSPTR( status );
 LALDDestroyVector(status->statusPtr,&hCY.data);
 CHECKSTATUSPTR( status );

 /* destroy low and high pass filters */ 
 LALDDestroyVector(status->statusPtr,&(LPFIR.directCoef));
 CHECKSTATUSPTR( status );
 LALDDestroyVector(status->statusPtr,&(LPFIR.recursCoef));
 CHECKSTATUSPTR( status );
 LALDDestroyVector(status->statusPtr,&(LPFIR.history));
 CHECKSTATUSPTR( status );
 LALDDestroyVector(status->statusPtr,&(HPFIR.directCoef));
 CHECKSTATUSPTR( status );
 LALDDestroyVector(status->statusPtr,&(HPFIR.recursCoef));
 CHECKSTATUSPTR( status );
 LALDDestroyVector(status->statusPtr,&(HPFIR.history));
 CHECKSTATUSPTR( status );

 DETATCHSTATUSPTR( status );
 RETURN( status );

}


/*******************************************************************************/
void LALFreeFilter(LALStatus *status, REAL8IIRFilter *F2, int ORDER)
{
  int n;
  INITSTATUS( status, "LALCopyFilter", COMPUTESTRAINC );
  ATTATCHSTATUSPTR( status );

  for(n=0;n<ORDER;n++){
    LALDDestroyVector(status->statusPtr,&(F2[n].directCoef));
    CHECKSTATUSPTR( status );
    LALDDestroyVector(status->statusPtr, &(F2[n].recursCoef));
    CHECKSTATUSPTR( status );
    LALDDestroyVector(status->statusPtr,&(F2[n].history));   
    CHECKSTATUSPTR( status );
  }
  LALFree(F2);

  DETATCHSTATUSPTR( status );
  RETURN( status );
}

/*******************************************************************************/

void LALCopyFilter(LALStatus *status, REAL8IIRFilter **F2, REAL8IIRFilter *F1, int ORDER)
{
  REAL8IIRFilter *F2Array;
  int n;
  INITSTATUS( status, "LALCopyFilter", COMPUTESTRAINC );
  ATTATCHSTATUSPTR( status );

  /* Allocate inverse sensing funtion filters */
  *F2=F2Array=(REAL8IIRFilter *)LALMalloc(sizeof(REAL8IIRFilter) * ORDER); 

  for(n = 0; n < ORDER; n++)
    {
      int l;
      F2Array[n].directCoef=NULL;
      F2Array[n].recursCoef=NULL;
      F2Array[n].history=NULL;

      LALDCreateVector(status->statusPtr,&(F2Array[n].directCoef),F1[n].directCoef->length);
      CHECKSTATUSPTR( status );
      LALDCreateVector(status->statusPtr,&(F2Array[n].recursCoef),F1[n].recursCoef->length);
      CHECKSTATUSPTR( status );
      LALDCreateVector(status->statusPtr,&(F2Array[n].history),F1[n].history->length);
      CHECKSTATUSPTR( status );

      for(l=0;l<(int)F1[n].directCoef->length;l++) 
	F2Array[n].directCoef->data[l] = F1[n].directCoef->data[l];
      for(l=0;l<(int)F1[n].recursCoef->length;l++) 
	F2Array[n].recursCoef->data[l] = F1[n].recursCoef->data[l];
      for(l=0;l<(int)F1[n].history->length;l++) 
	F2Array[n].history->data[l] = F1[n].history->data[l];
    }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}

/*******************************************************************************/
void LALMakeFIRLP(LALStatus *status, REAL8IIRFilter *LPFIR, int USF)
{
  int N=2*N_FIR_LP+1,l;
  int k[2*N_FIR_LP+1];
  REAL8 fN=fhigh_FIRLP/USF;

  INITSTATUS( status, "LALMakeFIRLP", COMPUTESTRAINC );
  ATTATCHSTATUSPTR( status );

  LPFIR->directCoef=NULL;
  LPFIR->recursCoef=NULL;
  LPFIR->history=NULL;

  LALDCreateVector(status->statusPtr,&(LPFIR->directCoef),N);
  CHECKSTATUSPTR( status );
  LALDCreateVector(status->statusPtr,&(LPFIR->recursCoef),N);
  CHECKSTATUSPTR( status );
  LALDCreateVector(status->statusPtr,&(LPFIR->history),N-1);
  CHECKSTATUSPTR( status );

  for(l=0;l<N;l++) LPFIR->recursCoef->data[l]=0.0;
  for(l=0;l<N-1;l++) LPFIR->history->data[l]=0.0;

  for (l=0; l<N;l++)
    {
      k[l]=l-N_FIR_LP;
      if(k[l] != 0)
	{
	  LPFIR->directCoef->data[l]=(sin(LAL_PI*fN*k[l])/(LAL_PI*k[l]))*
	    exp(-0.5*pow(3.0*(double)k[l]/(double)N_FIR_LP,2));
	}else{
	  LPFIR->directCoef->data[l]=fN;
	}
    }


/*   for(l=0;l<N;l++) fprintf(stdout,"%1.16e\n", LPFIR->directCoef->data[l]); */

  DETATCHSTATUSPTR( status );
  RETURN( status );

}

/*******************************************************************************/
void LALMakeFIRLPALPHAS(LALStatus *status, REAL8IIRFilter *LPFIR)
{
  int N=2*N_FIR_LP_ALPHAS+1,l;
  int k[2*N_FIR_LP_ALPHAS+1];
  REAL8 fN=fhigh_FIRLP_ALPHAS;

  INITSTATUS( status, "LALMakeFIRLPALPHAS", COMPUTESTRAINC );
  ATTATCHSTATUSPTR( status );

  LPFIR->directCoef=NULL;
  LPFIR->recursCoef=NULL;
  LPFIR->history=NULL;

  LALDCreateVector(status->statusPtr,&(LPFIR->directCoef),N);
  CHECKSTATUSPTR( status );
  LALDCreateVector(status->statusPtr,&(LPFIR->recursCoef),N);
  CHECKSTATUSPTR( status );
  LALDCreateVector(status->statusPtr,&(LPFIR->history),N-1);
  CHECKSTATUSPTR( status );

  for(l=0;l<N;l++) LPFIR->recursCoef->data[l]=0.0;
  for(l=0;l<N-1;l++) LPFIR->history->data[l]=0.0;

  for (l=0; l<N;l++)
    {
      k[l]=l-N_FIR_LP_ALPHAS;
      if(k[l] != 0)
	{
	  LPFIR->directCoef->data[l]=(sin(LAL_PI*fN*k[l])/(LAL_PI*k[l]))*
	    exp(-0.5*pow(3.0*(double)k[l]/(double)N_FIR_LP_ALPHAS,2))/8.318081379762647e-02;
	}else{
	  LPFIR->directCoef->data[l]=fN;
	}
    }


/*   for(l=0;l<N;l++) fprintf(stdout,"%1.16e\n", LPFIR->directCoef->data[l]); */

  DETATCHSTATUSPTR( status );
  RETURN( status );

}
/*******************************************************************************/

void LALMakeFIRHP(LALStatus *status, REAL8IIRFilter *HPFIR)
{
  int N=2*N_FIR_HP+1,l;
  int k[2*N_FIR_HP+1];
  REAL8 fN=fhigh_FIRHP;

  INITSTATUS( status, "LALMakeFIRHP", COMPUTESTRAINC );
  ATTATCHSTATUSPTR( status );

  HPFIR->directCoef=NULL;
  HPFIR->recursCoef=NULL;
  HPFIR->history=NULL;

  LALDCreateVector(status->statusPtr,&(HPFIR->directCoef),N);
  CHECKSTATUSPTR( status );
  LALDCreateVector(status->statusPtr,&(HPFIR->recursCoef),N);
  CHECKSTATUSPTR( status );
  LALDCreateVector(status->statusPtr,&(HPFIR->history),N-1);
  CHECKSTATUSPTR( status );

  for(l=0;l<N;l++) HPFIR->recursCoef->data[l]=0.0;
  for(l=0;l<N-1;l++) HPFIR->history->data[l]=0.0;

  for (l=0; l<N;l++)
    {
      k[l]=l-N_FIR_HP;
      if(k[l] != 0)
	{
	  HPFIR->directCoef->data[l]=(sin(LAL_PI*k[l])-sin(LAL_PI*fN*k[l]))/(LAL_PI*k[l])*
	    exp(-0.5*pow(3.0*(double)k[l]/(double)N_FIR_HP,2));
	}else{
	  HPFIR->directCoef->data[l]=1.0-fN;
	}
    }


/*   for(l=0;l<N;l++) fprintf(stdout,"%1.16e\n", HPFIR->directCoef->data[l]); */

  DETATCHSTATUSPTR( status );
  RETURN( status );

}

/*******************************************************************************/

int XLALUpsample(REAL8TimeSeries *uphR, REAL8TimeSeries *hR, int up_factor)
{
  int n;

  /* Set all values to 0 */
  for (n=0; n < (int)uphR->data->length; n++) {
    uphR->data->data[n] = 0.0;
  }

  /* Set one in every up_factor to the value of hR x USR */
  for (n=0; n < (int)hR->data->length; n++) {
    uphR->data->data[n * up_factor] = up_factor * hR->data->data[n];
  }

  return 0;
}

/*******************************************************************************/

int XLALUpsampleLinear(REAL8TimeSeries *uphR, REAL8TimeSeries *hR, int up_factor)
{
  int n,m;

  /* Set all values to 0 */
  for (n=0; n < (int)uphR->data->length; n++) {
    uphR->data->data[n] = 0.0;
  }

  /* Set one in every up_factor to the value of hR x USR */
  for (n=0; n < (int)hR->data->length; n++) 
    {
      REAL8 y1=hR->data->data[n],y2;
     
      if(n < hR->data->length-1) y2=hR->data->data[n+1];
      if(n == hR->data->length-1) y2=hR->data->data[n];

      for (m=0; m < (int)up_factor; m++) 
	{
	  uphR->data->data[n*up_factor+m] = y1+m*(y2-y1)/up_factor;
	}
    }
  return 0;
}

/*******************************************************************************/

int XLALDivideTimeSeries(REAL8TimeSeries *hR, REAL8TimeSeries *ALPHAS)
{

  int n;

  if (hR->data->length != ALPHAS->data->length)
    {
      fprintf(stderr,"Length of residual strin time series (%d), not the same as factors time series, (%d)\n"
	      ,hR->data->length, ALPHAS->data->length);
      return 1;
    }

  for (n = 0; n < (int)hR->data->length; n++) 
    {
      hR->data->data[n] /= ALPHAS->data->data[n];
    }

  return 0;
}

/*******************************************************************************/

void LALGetFactors(LALStatus *status, StrainOut *output, StrainIn *input)
{

static REAL4TimeSeries darm;
static REAL4TimeSeries asq;
static REAL4TimeSeries exc;

CalFactors factors;
UpdateFactorsParams params;

REAL4Vector *asqwin=NULL,*excwin=NULL,*darmwin=NULL;  /* windows */
LALWindowParams winparams;
INT4 k,m;
 
REAL4 deltaT=input->AS_Q.deltaT, To=input->To;
INT4 length = input->AS_Q.data->length;
INT4 localtime = input->AS_Q.epoch.gpsSeconds;

 INITSTATUS( status, "LALGetFactors", COMPUTESTRAINC );
 ATTATCHSTATUSPTR( status );

  /* Create local data vectors */
  LALCreateVector(status->statusPtr,&asq.data,(UINT4)(To/input->AS_Q.deltaT +0.5));
  CHECKSTATUSPTR( status );
  LALCreateVector(status->statusPtr,&darm.data,(UINT4)(To/input->DARM.deltaT +0.5));
  CHECKSTATUSPTR( status );
  LALCreateVector(status->statusPtr,&exc.data,(UINT4)(To/input->EXC.deltaT +0.5));
  CHECKSTATUSPTR( status );

  /* Create Window vectors */
  LALCreateVector(status->statusPtr,&asqwin,(UINT4)(To/input->AS_Q.deltaT +0.5));
  CHECKSTATUSPTR( status );
  LALCreateVector(status->statusPtr,&darmwin,(UINT4)(To/input->DARM.deltaT +0.5));
  CHECKSTATUSPTR( status );
  LALCreateVector(status->statusPtr,&excwin,(UINT4)(To/input->EXC.deltaT +0.5));
  CHECKSTATUSPTR( status );

  /* assign time spacing for local time series */
  asq.deltaT=input->AS_Q.deltaT;
  darm.deltaT=input->DARM.deltaT;
  exc.deltaT=input->EXC.deltaT;

  winparams.type=Hann;
   
  /* windows for time domain channels */
  /* asq */
  winparams.length=(INT4)(To/asq.deltaT +0.5);
  LALWindow(status->statusPtr,asqwin,&winparams);
  CHECKSTATUSPTR( status );
  
  /* darm */
  winparams.length=(INT4)(To/darm.deltaT +0.5);
  LALWindow(status->statusPtr,darmwin,&winparams);
  CHECKSTATUSPTR( status );

  /* exc */
  winparams.length=(INT4)(To/exc.deltaT +0.5);
  LALWindow(status->statusPtr,excwin,&winparams);
  CHECKSTATUSPTR( status );

  for(m=0; m < (UINT4)(deltaT*length) / To; m++)
    {
      int facterrflag=0;

      /* assign and window the data */
      for(k=0;k<(INT4)(To/asq.deltaT +0.5);k++)
	{
	  asq.data->data[k] = input->AS_Q.data->data[m * (UINT4)(To/asq.deltaT) + k] * 2.0 * asqwin->data[k];
	}
      for(k=0;k<(INT4)(input->To/darm.deltaT +0.5);k++)
	{
	  darm.data->data[k] = input->DARM.data->data[m *(UINT4)(To/darm.deltaT) + k] * 2.0 * darmwin->data[k];
	}
      for(k=0;k<(INT4)(input->To/exc.deltaT +0.5);k++)
	{
	  exc.data->data[k] = input->EXC.data->data[m * (UINT4)(To/exc.deltaT) + k] * 2.0 * excwin->data[k];
	}

      /* set params to call LALComputeCalibrationFactors */
      params.darmCtrl = &darm;
      params.asQ = &asq;
      params.exc = &exc;

      params.lineFrequency = input->f;
      params.openloop =  input->Go;
      params.digital = input->Do;

      LALComputeCalibrationFactors(status->statusPtr,&factors,&params);
      CHECKSTATUSPTR( status );

      /* check alpha */
      if( (factors.alpha.re < 0.3) ||  (factors.alpha.re > 2.0) ) 
	{
	 factors.alpha.re = 1.0;
	 if (m>0) factors.alpha.re=output->alpha.data->data[m-1].re;
	 facterrflag=1;
	}

      /* check alphabeta */
      if( (factors.alphabeta.re < 0.3) ||  (factors.alphabeta.re > 2.0) ) 
	{
	 factors.alphabeta.re = 1.0;
	 if (m>0) factors.alphabeta.re=output->alphabeta.data->data[m-1].re;
	 facterrflag=1;
	}

      /* check beta */
      if( (factors.beta.re < 0.3) ||  (factors.beta.re > 2.0) ) 
	{
	 factors.beta.re = 1.0;
	 if (m>0) factors.beta.re=output->beta.data->data[m-1].re;
	 facterrflag=1;
	}
      
      localtime=localtime+(int)To;

      if (facterrflag == 1)
	{
	  fprintf(stderr,"BADFACTORS: %d %e %e %e %e %e %e %e %e %e %e %e %e\n",localtime, 
		  factors.alphabeta.re,factors.alphabeta.im,
		  factors.beta.re,factors.beta.im,
		  factors.alpha.re,factors.alpha.im,
		  factors.asq.re*2/To,factors.asq.im*2/To,
		  factors.darm.re*2/To,factors.darm.im*2/To,
		  factors.exc.re*2/To,factors.exc.im*2/To);
	  facterrflag=0;
	}

      output->alpha.data->data[m]= factors.alpha;
      output->beta.data->data[m]= factors.beta;
      output->alphabeta.data->data[m]= factors.alphabeta;

      if(m == MAXALPHAS)
	{
	  fprintf(stderr,"Too many values of the factors, maximum allowed is %d\n",MAXALPHAS);
	  RETURN(status);
	}
    }

  /* Clean up */
  LALDestroyVector(status->statusPtr,&darm.data);
  CHECKSTATUSPTR( status );
  LALDestroyVector(status->statusPtr,&exc.data);
  CHECKSTATUSPTR( status );
  LALDestroyVector(status->statusPtr,&asq.data);
  CHECKSTATUSPTR( status );

  LALDestroyVector(status->statusPtr,&asqwin);
  CHECKSTATUSPTR( status );
  LALDestroyVector(status->statusPtr,&darmwin);
  CHECKSTATUSPTR( status );
  LALDestroyVector(status->statusPtr,&excwin);
  CHECKSTATUSPTR( status );

  DETATCHSTATUSPTR( status );
  RETURN (status);
}



/*******************************************************************************/

void LALFIRFilter(LALStatus *status, REAL8TimeSeries *tseries, REAL8IIRFilter *FIR)
{
  REAL8TimeSeries *tseriesFIR=NULL;
  REAL8TimeSeries *tseriesDATA=NULL;
  COMPLEX16FrequencySeries *vtilde=NULL, *vtildeFIR=NULL;
  REAL8FFTPlan  *fplan = NULL, *rplan = NULL;
  int n;
  int xlerr;

  INITSTATUS( status, "LALFIRFilter", COMPUTESTRAINC );
  ATTATCHSTATUSPTR( status );
  
  /* create time series that will hold FIR filter (with room for zero padding) */
  LALCreateREAL8TimeSeries(status->statusPtr, &tseriesFIR, tseries->name, 
			   tseries->epoch, 0.0, tseries->deltaT, tseries->sampleUnits, 
			   tseries->data->length+FIR->directCoef->length);
  CHECKSTATUSPTR( status );

  /* create time series that will hold data  (with room for zero padding) */
  LALCreateREAL8TimeSeries(status->statusPtr, &tseriesDATA, tseries->name, 
			   tseries->epoch, 0.0, tseries->deltaT, tseries->sampleUnits, 
			   tseries->data->length+FIR->directCoef->length);
  CHECKSTATUSPTR( status );

  /* initialise values in FIR time series */
  for (n = 0; n < (int)tseriesDATA->data->length; n++) 
    {
      tseriesFIR->data->data[n]=0.0;
      tseriesDATA->data->data[n]=0.0;
    }
  /* set first few to values in FIR filter */
  /* Here I want to make sure the FIR filter in in wrap-around order */
/*   for (n = 0; n <= (int)(FIR->directCoef->length-1)/2; n++)  */
/*     { */
/*       tseriesFIR->data->data[n]=FIR->directCoef->data[n+(FIR->directCoef->length-1)/2]; */
/*     } */
/*   for (n = 0; n < (int) (FIR->directCoef->length-1)/2; n++)  */
/*     { */
/*       tseriesFIR->data->data[(tseriesFIR->data->length-1)-n]=FIR->directCoef->data[((FIR->directCoef->length-1)/2-1)-n]; */
/*     } */
  for (n = 0; n < (int)FIR->directCoef->length; n++) 
    {
      tseriesFIR->data->data[n]=FIR->directCoef->data[n];
    }

  /* set first few to values in data series */
  for (n = 0; n < (int)tseries->data->length; n++) 
    {
      tseriesDATA->data->data[n]=tseries->data->data[n];
    }

  /* create frequency series that will hold FT's of both timne series */ 
  LALCreateCOMPLEX16FrequencySeries(status->statusPtr, &vtilde,tseries->name , tseries->epoch, 0.0, 
				   1.0 / (tseries->data->length * tseries->deltaT), 
				   tseries->sampleUnits, tseriesDATA->data->length / 2 + 1);
  CHECKSTATUSPTR( status );
  LALCreateCOMPLEX16FrequencySeries(status->statusPtr, &vtildeFIR,tseries->name , tseries->epoch, 0.0, 
				   1.0 / (tseries->data->length * tseries->deltaT), 
				   tseries->sampleUnits, tseriesDATA->data->length / 2 + 1);
  CHECKSTATUSPTR( status );

  /* make fft plans */
  LALCreateForwardREAL8FFTPlan(status->statusPtr, &fplan, tseriesDATA->data->length, 0 );
  CHECKSTATUSPTR( status );
  LALCreateReverseREAL8FFTPlan(status->statusPtr, &rplan, tseriesDATA->data->length, 0 );
  CHECKSTATUSPTR( status );

  /* fft both series */
  xlerr=XLALREAL8TimeFreqFFT(vtilde, tseriesDATA, fplan);
  if( xlerr < 0 )
    {   
      fprintf(stderr,"Failed creating FT of time series. Errorcode: %d\n",xlerr);
      DETATCHSTATUSPTR( status );
      RETURN(status);
    }
  xlerr=XLALREAL8TimeFreqFFT(vtildeFIR, tseriesFIR, fplan);
  if( xlerr < 0 )
    {   
      fprintf(stderr,"Failed creating FT of FIR filter time series. Errorcode: %d\n",xlerr);
      DETATCHSTATUSPTR( status );
      RETURN(status);
    }
 
  /* multiply both FT's */
  for (n = 0; n < (int)vtilde->data->length; n++) 
    {
      REAL8 re=vtilde->data->data[n].re, im=vtilde->data->data[n].im;
      REAL8 reFIR=vtildeFIR->data->data[n].re, imFIR=vtildeFIR->data->data[n].im;
      vtilde->data->data[n].re=re*reFIR-im*imFIR;
      vtilde->data->data[n].im=re*imFIR+im*reFIR;
    }

  /* reverse FFT into original time series */
  xlerr=XLALREAL8FreqTimeFFT(tseriesDATA, vtilde, rplan);
  if( xlerr < 0 )
    {   
      fprintf(stderr,"Failed creating IFT of FIR and time series. Errorcode: %d\n",xlerr);
      DETATCHSTATUSPTR( status );
      RETURN(status);
    }

  for (n = 0; n < (int)tseries->data->length; n++) 
    {
      tseries->data->data[n] = tseriesDATA->data->data[n]/tseries->deltaT;
    }

  /* Destroy everything */
  LALDestroyREAL8FFTPlan( status->statusPtr, &fplan );
  CHECKSTATUSPTR( status );
  LALDestroyREAL8FFTPlan( status->statusPtr, &rplan );
  CHECKSTATUSPTR( status );

  LALDestroyCOMPLEX16FrequencySeries(status->statusPtr, vtilde);
  CHECKSTATUSPTR( status );
  LALDestroyCOMPLEX16FrequencySeries(status->statusPtr, vtildeFIR);
  CHECKSTATUSPTR( status );
  LALDestroyREAL8TimeSeries(status->statusPtr, tseriesFIR);
  CHECKSTATUSPTR( status );
  LALDestroyREAL8TimeSeries(status->statusPtr, tseriesDATA);
  CHECKSTATUSPTR( status );

  DETATCHSTATUSPTR( status );
  RETURN( status );

}

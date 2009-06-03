/*
*  Copyright (C) 2007 Xavier Siemens
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
#define fhigh_FIRLP .9

#define N_FIR_HP 2000
#define fhigh_FIRHP 0.00244140625

NRCSID( COMPUTESTRAINC, "$Id$" );

void LALComputeStrainDMT(
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
REAL8IIRFilter *CinvWings=NULL, *DWings=NULL, *AAWings=NULL, *AXWings=NULL, *AYWings=NULL;


 INITSTATUS( status, "LALComputeStrainDMT", COMPUTESTRAINC );
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

 /* copy DARM_ERR input into control strain as double */
 for (p=0; p<(int)output->hC.data->length; p++)
   {
     output->hC.data->data[p]=input->DARM_ERR.data->data[p];
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

 /* An odd filter with N points introduces an (N-1)/2 delay */
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
     LALFFTFIRFilter(status->statusPtr,&upALPHAS,&ALPHASLPFIR);
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
 /* then high pass filter DARM_ERR */
 if (input->fftconv)
   {
     LALFFTFIRFilter(status->statusPtr,&(output->hC),&(HPFIR));
     CHECKSTATUSPTR( status );
     LALFFTFIRFilter(status->statusPtr,&(output->hC),&(HPFIR));
     CHECKSTATUSPTR( status );
   }else{
     LALIIRFilterREAL8Vector(status->statusPtr,output->hC.data,&(HPFIR));
     CHECKSTATUSPTR( status );
     LALIIRFilterREAL8Vector(status->statusPtr,output->hC.data,&(HPFIR));
     CHECKSTATUSPTR( status );
   }

 /* Filter through digital servo */
 /* ===================== */
 /* CAREFUL FILTERING: FILTER UP TO WINGS, THEN COPY FILTERS THEN CONTINUE UNTIL END*/
 /* we filter but only up until the wings */
 for(p = 0; p < input->ND; p++){
   int k;
   for(k = NWingsC/2; k < (int)output->hC.data->length-3*NWingsC/2; k++){
     output->hC.data->data[k]=LALDIIRFilter(output->hC.data->data[k], &(input->D[p]));
   }
 }
 /* Here what I need to record filter histories */
 LALCopyFilter(status->statusPtr, &DWings, input->D, input->ND);
 CHECKSTATUSPTR( status );
 /* then we filter wings as well */
 for(p = 0; p < input->ND; p++){
   int k;
   for(k = output->hC.data->length-3*NWingsC/2; k < (int)output->hC.data->length-NWingsC/2; k++){
     output->hC.data->data[k]=LALDIIRFilter(output->hC.data->data[k], &DWings[p]);
   }
 }
 /* Then we need to destroy the filter */
 LALFreeFilter(status->statusPtr,DWings,input->ND);
 CHECKSTATUSPTR( status );
 /* ===================== */

 /* At this point I have something very similar to darm_ctrl in hC */
 /* add the calibration lines */
 if (!input->delta && (input->DARM_ERR.deltaT == input->EXC.deltaT))
   {
     int k;
     for(k = 0; k < (int)output->hC.data->length; k++){
       output->hC.data->data[k] += input->EXC.data->data[k];
     }
   }
 else
   {
     fprintf(stdout, "Warning: Not adding calibration lines to control signal.\n");
   }


 /* Filter through analog actuation */
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
    output->hC.data->data[p]=(hCX.data->data[p]-hCY.data->data[p]);
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


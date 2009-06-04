/*
*  Copyright (C) 2007 Jolien Creighton, Xavier Siemens
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
/* #define fhigh_FIRLP .6103515625 */
#define fhigh_FIRLP .9

#define N_FIR_HP 2000
#define fhigh_FIRHP 0.00244140625

NRCSID( COMPUTESTRAINC, "$Id$" );
RCSID("$Id$");

void LALComputeStrain(
    LALStatus              *status,
    StrainOut              *output,
    StrainIn               *input
    )
{
/* Inverse sensing, servo, analog actuation, digital x
actuation  digital y actuation */
static REAL8TimeSeries uphR,ALPHAS,upALPHAS;
int p;
REAL8IIRFilter LPFIR;
REAL8IIRFilter HPFIR;
REAL8IIRFilter ALPHASLPFIR;

 INITSTATUS( status, "LALComputeStrain", COMPUTESTRAINC );
 ATTATCHSTATUSPTR( status );

 LALGetFactors(status->statusPtr, output, input);
 CHECKSTATUSPTR( status );

 LALMakeFIRLP(status->statusPtr, &LPFIR, input->CinvUSF);
 CHECKSTATUSPTR( status );

 LALMakeFIRHP(status->statusPtr, &HPFIR);
 CHECKSTATUSPTR( status );

 /* copy DARM_ERR input into residual strain as double */
 for (p=0; p<(int)output->hR.data->length; p++)
   {
    output->hR.data->data[p]=input->DARM_ERR.data->data[p];
   }

 if (input->darmctrl)
   {
     /* copy DARM_CTRL input into control strain as double */
     for (p=0; p<(int)output->hC.data->length; p++)
       {
	 output->hC.data->data[p]=input->DARM.data->data[p];
       }
     /* set to 0 the sampleUnits: they are not used anyway, and
      * leaving them uninitialized can cause a mess */
     output->hC.sampleUnits.powerOfTen = 0;
     for (p=0; p < LALNumUnits; p++)
       {
         output->hC.sampleUnits.unitNumerator[p] = 0;
         output->hC.sampleUnits.unitDenominatorMinusOne[p] = 0;
       }
   }
 else
   {
     /* If DARM_ERR only calibration copy DARM_ERR data into control
	signal */
     for (p=0; p<(int)output->hC.data->length; p++)
       {
	 output->hC.data->data[p]=input->DARM_ERR.data->data[p];
       }
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
 XLALFIRFilter(&uphR,&LPFIR);
 XLALFIRFilter(&uphR,&LPFIR);

 /* Filter through inverse of sensing function */
 XLALFIRFilter(&uphR, input->Cinv);

 /* apply advance to compensate for Low Pass FIR delay */
 for (p=0; p<(int)uphR.data->length-(2*N_FIR_LP); p++){
   uphR.data->data[p]=uphR.data->data[p+(2*N_FIR_LP)];
 }

 /* Low pass filter twice to smooth time series (again) */
 XLALFIRFilter(&uphR,&LPFIR);
 XLALFIRFilter(&uphR,&LPFIR);

 /* then we downsample and voila' */
 for (p=0; p<(int)output->hR.data->length; p++) {
   output->hR.data->data[p]=uphR.data->data[p*input->CinvUSF];
 }
 LALDDestroyVector(status->statusPtr,&uphR.data);
 CHECKSTATUSPTR( status );

 /* Create time series that hold alpha time series and upsampled alpha time-series */
 LALDCreateVector(status->statusPtr,&ALPHAS.data,output->alpha.data->length);
 CHECKSTATUSPTR( status );
 LALDCreateVector(status->statusPtr,&upALPHAS.data,input->DARM_ERR.data->length);
 CHECKSTATUSPTR( status );
 upALPHAS.deltaT=input->DARM_ERR.deltaT;

 /* copy factors into time series */
 for (p=0; p<(int)ALPHAS.data->length; p++)
   {
     REAL8 r = output->alphabeta.data->data[p].re;

     /* check alphabeta: If values are outside bounds we replace
      * factors with the last one, or with 1 if it is the first */
     if ( (r < 0.8) ||  (r > 1.2) || (isnan(r)) || isinf(r))
       {
         if (p > 0)
           ALPHAS.data->data[p] = ALPHAS.data->data[p-1];
         else
           ALPHAS.data->data[p] = 1.0;
       }
     else
       {
         ALPHAS.data->data[p] = r;  /* this is the "standard" nice case */
       }
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
     LALFFTFIRFilter(status->statusPtr,&upALPHAS,&(ALPHASLPFIR));
     CHECKSTATUSPTR( status );
   }
 else
   {
     XLALFIRFilter(&upALPHAS,&(ALPHASLPFIR));
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
 /* apply advance to compensate for FIR delay */
 for (p=0; p<(int)output->hC.data->length-(2*N_FIR_HP); p++){
   output->hC.data->data[p]=output->hC.data->data[p+(2*N_FIR_HP)];
 }

 /* then high-pass filter DARM_CTRL */
 if (input->fftconv)
   {
     LALFFTFIRFilter(status->statusPtr,&(output->hC),&(HPFIR));
     CHECKSTATUSPTR( status );
     LALFFTFIRFilter(status->statusPtr,&(output->hC),&(HPFIR));
     CHECKSTATUSPTR( status );
   }
 else
   {
     XLALFIRFilter(&(output->hC),&(HPFIR));
     XLALFIRFilter(&(output->hC),&(HPFIR));
   }

 if (input->darmctrl)
   {
     /* Filter through anti-whitening */
     if (input->fftconv)
       {
	 LALFFTFIRFilter(status->statusPtr,&(output->hC), input->AW);
	 CHECKSTATUSPTR( status );
       }
     else
       {
	 XLALFIRFilter(&(output->hC), input->A);
       }
   }else
   {
     /* Filter through servo */
     if (input->fftconv)
       {
	 LALFFTFIRFilter(status->statusPtr,&(output->hC), input->D);
	 CHECKSTATUSPTR( status );
       }
     else
       {
	 XLALFIRFilter(&(output->hC), input->D);
       }

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
   }

 /* Filter through analog actuation */
 if (input->fftconv)
   {
     LALFFTFIRFilter(status->statusPtr,&(output->hC), input->A);
     CHECKSTATUSPTR( status );
   }
 else
   {
     XLALFIRFilter(&(output->hC), input->A);
   }

 /* ---------- Compute Net Strain -------------*/
 /* add control and residual signals and voila' */
 for (p=0; p< (int)output->h.data->length; p++){
   output->h.data->data[p] = output->hC.data->data[p]+output->hR.data->data[p];
 }
 /* ------------------------------------------*/

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
  UINT4 n,m;

  /* Set all values to 0 */
  for (n=0; n < (UINT4)uphR->data->length; n++) {
    uphR->data->data[n] = 0.0;
  }

  /* Set one in every up_factor to the value of hR x USR */
  for (n=0; n < (UINT4)hR->data->length; n++)
    {
      REAL8 y_1=hR->data->data[n],y_2=0;

      if(n < hR->data->length-1) y_2=hR->data->data[n+1];
      if(n == hR->data->length-1) y_2=hR->data->data[n];

      for (m=0; m < (UINT4)up_factor; m++)
	{
	  uphR->data->data[n*up_factor+m] = y_1+m*(y_2-y_1)/up_factor;
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
      params.whitener = input->Wo;

      LALComputeCalibrationFactors(status->statusPtr,&factors,&params);
      CHECKSTATUSPTR( status );

      if (input->gamma_fudgefactor != 0)
	{
	  factors.alphabeta.re /= input->gamma_fudgefactor;
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

int XLALFIRFilter(REAL8TimeSeries *tseries, REAL8IIRFilter *FIR)
{
  int n,m;
  REAL8 sum;
  int Nfir,Ntseries;
  REAL8 *x, *b;

  x=tseries->data->data;
  b=FIR->directCoef->data;

  Nfir=FIR->directCoef->length;
  Ntseries=tseries->data->length;

  /* initialise values in FIR time series */
  for (n = Ntseries-1; n >= Nfir-1; n--)
    {
      sum = 0;
      for (m = Nfir-1; m >= 0; m--)
	{
	  sum += b[m] * x[n-m];
	}
      x[n]=sum;
    }
  /* set to zero values at the start */
  for (n = 0; n < (int)FIR->directCoef->length-1; n++)
    {
      x[n]=0;
    }
  return 0;
}

/*******************************************************************************/

void LALFFTFIRFilter(LALStatus *status, REAL8TimeSeries *tseries, REAL8IIRFilter *FIR)
{

  REAL8TimeSeries *tseriesFIR=NULL;
  REAL8TimeSeries *tseriesDATA=NULL;
  COMPLEX16FrequencySeries *vtilde=NULL, *vtildeFIR=NULL;
  REAL8FFTPlan  *fplan = NULL, *rplan = NULL;
  int n;
  int xlerr;

  INITSTATUS( status, "LALFIRFilter", COMPUTESTRAINC );
  ATTATCHSTATUSPTR( status );

  /* check that the time series is larger than the length of the filter */
  if ( tseries->data->length < FIR->directCoef->length )
    {
      fprintf(stderr,"ERROR. tseries length (=%d) < filter length (=%d)",
	      tseries->data->length, FIR->directCoef->length);
      ABORT(status,118,"Time series is smaller than filter length.");
    }

  /* create time series that will hold FIR filter (with room for zero padding) */
  tseriesFIR = XLALCreateREAL8TimeSeries(tseries->name,
			   &tseries->epoch, 0.0, tseries->deltaT, &tseries->sampleUnits,
			   2*tseries->data->length);
  if(!tseriesFIR)
    ABORTXLAL( status );

  /* create time series that will hold data  (with room for zero padding) */
  tseriesDATA = XLALCreateREAL8TimeSeries(tseries->name,
			   &tseries->epoch, 0.0, tseries->deltaT, &tseries->sampleUnits,
			   2*tseries->data->length);
  if(!tseriesDATA)
    ABORTXLAL( status );

  /* initialise values in FIR time series */
  for (n = 0; n < (int)tseriesDATA->data->length; n++)
    {
      tseriesFIR->data->data[n]=0.0;
      tseriesDATA->data->data[n]=0.0;
    }
  /* set first few to values in FIR filter */
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
  vtilde = XLALCreateCOMPLEX16FrequencySeries(tseries->name , &tseries->epoch, 0.0,
				   1.0 / (tseries->data->length * tseries->deltaT),
				   &tseries->sampleUnits, tseriesDATA->data->length / 2 + 1);
  if(!vtilde)
    ABORTXLAL( status );
  vtildeFIR = XLALCreateCOMPLEX16FrequencySeries(tseries->name , &tseries->epoch, 0.0,
				   1.0 / (tseries->data->length * tseries->deltaT),
				   &tseries->sampleUnits, tseriesDATA->data->length / 2 + 1);
  if(!vtildeFIR)
    ABORTXLAL( status );

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

  XLALDestroyCOMPLEX16FrequencySeries(vtilde);
  XLALDestroyCOMPLEX16FrequencySeries(vtildeFIR);
  XLALDestroyREAL8TimeSeries(tseriesFIR);
  XLALDestroyREAL8TimeSeries(tseriesDATA);

  DETATCHSTATUSPTR( status );
  RETURN( status );

}


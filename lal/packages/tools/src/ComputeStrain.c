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
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Calibration.h>
#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/AVFactories.h>
#include <lal/Window.h>
#include <lal/BandPassTimeSeries.h>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include <lal/filters-H1-S3.h>

#define UpSamplingFactor 16
#define MAXALPHAS 10000

NRCSID( COMPUTESTRAINC, "$Id$" );

void LALComputeStrain(
    LALStatus              *status,
    StrainOut              *output,    
    StrainIn               *input
    )
{
/* Inverse sensing, servo, analog actuation, digital x 
actuation  digital y actuation */
static REAL8IIRFilter Cinv,G[NGfilt],AA,AX[NAXfilt],AY[NAYfilt];    
static REAL8TimeSeries h,hR,uphR,hC,hCX,hCY;
int p;
PassBandParamStruc highpassfilterpar,lowpassfilterpar;

 INITSTATUS( status, "LALComputeStrain", COMPUTESTRAINC );
 ATTATCHSTATUSPTR( status );

  /* high pass filter parameters */ 
  highpassfilterpar.nMax  = 10;
  highpassfilterpar.f2    = 40.0;
  highpassfilterpar.a2    = 0.5;
  highpassfilterpar.f1    = -1.0;
  highpassfilterpar.a1    = -1.0;
  /* low pass filter parameters */
  lowpassfilterpar.nMax  = 12;
  lowpassfilterpar.f2    = -1.0;
  lowpassfilterpar.a2    = -1.0;
  lowpassfilterpar.f1    = 6000.0;
  lowpassfilterpar.a1    = 0.5;

  XLALMakeFilters(status,&Cinv,G,&AA,AX,AY); 
  CHECKSTATUSPTR( status );

  XLALGetFactors(status, output, input);
  CHECKSTATUSPTR( status );
  
  /* Create vectors that will hold the residual, control and net strain signals */
  LALDCreateVector(status->statusPtr,&h.data,input->AS_Q.data->length);
  CHECKSTATUSPTR( status );
  LALDCreateVector(status->statusPtr,&hR.data,input->AS_Q.data->length);
  CHECKSTATUSPTR( status );
  LALDCreateVector(status->statusPtr,&hC.data,input->AS_Q.data->length);
  CHECKSTATUSPTR( status );
  LALDCreateVector(status->statusPtr,&hCX.data,input->AS_Q.data->length);
  CHECKSTATUSPTR( status );
  LALDCreateVector(status->statusPtr,&hCY.data,input->AS_Q.data->length);
  CHECKSTATUSPTR( status );
  LALDCreateVector(status->statusPtr,&uphR.data,UpSamplingFactor*input->AS_Q.data->length);
  CHECKSTATUSPTR( status );

  h.deltaT=input->AS_Q.deltaT;
  hR.deltaT=input->AS_Q.deltaT;
  hC.deltaT=input->AS_Q.deltaT;
  uphR.deltaT=input->AS_Q.deltaT/UpSamplingFactor;

  /* copy AS_Q input into residual strain as double */  
  for (p=0; p<(int)hR.data->length; p++) {
    hR.data->data[p]=input->AS_Q.data->data[p];
  }
  /* copy AS_Q input into control strain as double */  
  for (p=0; p<(int)hC.data->length; p++) {
    hC.data->data[p]=input->AS_Q.data->data[p];
  }

  /* high pass filter both time series */
  LALButterworthREAL8TimeSeries(status->statusPtr,&hR,&highpassfilterpar);
  CHECKSTATUSPTR( status );
  LALButterworthREAL8TimeSeries(status->statusPtr,&hC,&highpassfilterpar);
  CHECKSTATUSPTR( status );
  
  /* ---------- Compute Residual Strain -------------*/
  /* to get the residual strain we must first divide AS_Q by alpha */
  if(XLALhROverAlpha(&hR, output)) 
    {
      ABORT(status,116,"Broke at hR/alpha") ;
    }
  /* then we upsample (and smooth it with a low pass filter) */
  if(XLALUpsamplehR(&uphR, &hR, UpSamplingFactor)) 
    { 
      ABORT(status,117,"Broke upsampling hR");
    }
  LALButterworthREAL8TimeSeries(status->statusPtr,&uphR,&lowpassfilterpar);
  CHECKSTATUSPTR( status );

  /* then we filter through the inverse of the sensing function */
  LALIIRFilterREAL8Vector(status->statusPtr,uphR.data,&Cinv);
  CHECKSTATUSPTR( status );

  /* Low pass again before downsampling */
  LALButterworthREAL8TimeSeries(status->statusPtr,&uphR,&lowpassfilterpar);
  CHECKSTATUSPTR( status );
 
  /* then we downsample and voila' */
  for (p=0; p<(int)hR.data->length; p++) {
    hR.data->data[p]=uphR.data->data[p*UpSamplingFactor];
  }
  
  /* ---------- Compute Control Strain -------------*/
  /* to get the control strain we first multiply AS_Q by beta */
  if( XLALhCTimesBeta(&hC, output))
    { 
      ABORT(status,120,"Broke in hC x beta");
    }

  /* Now we filter through the servo */
  for(p=NGfilt-1;p>=0;p--){
    LALIIRFilterREAL8Vector(status->statusPtr,hC.data,&G[p]);
    CHECKSTATUSPTR( status );
  }
  /* and adjust to account for servo gain to get darm_ctrl */ 
  for (p=0; p<(int)hC.data->length;p++) {
    hC.data->data[p] *= ServoGain;
  }

  /* Copy data into x and y time series for parallel filtering */
  for (p=0; p < (int)hCX.data->length; p++) {
    hCX.data->data[p]=hC.data->data[p];
  }
  for (p=0; p < (int)hCY.data->length; p++) {
    hCY.data->data[p]=hC.data->data[p];
  }
  
  /* Filter x-arm */
  for(p = NAXfilt-1; p >= 0; p--){
    LALIIRFilterREAL8Vector(status->statusPtr,hCX.data,&AX[p]);
    CHECKSTATUSPTR( status );
  }
  /* Adjust to account for digital gain on x-arm*/ 
  for (p=0; p< (int)hC.data->length;p++) {
    hCX.data->data[p] *= AXGain;
  }
  
  /* Filter y-arm */
  for(p = NAYfilt-1; p >= 0; p--){
    LALIIRFilterREAL8Vector(status->statusPtr,hCY.data,&AY[p]);
    CHECKSTATUSPTR( status );
  }
  /* Adjust to account for digital gain on y-arm*/ 
  for (p = 0; p < (int)hC.data->length; p++) {
    hCY.data->data[p] *= AYGain;
  }

  /* add x-arm and y-arm together */
  for (p=0; p< (int)hC.data->length; p++) {
    hC.data->data[p]=(hCX.data->data[p]+hCY.data->data[p])/2;
  }

  /* filter through analog part of actuation and voila' */ 
  LALIIRFilterREAL8Vector(status->statusPtr,hC.data,&AA);
  CHECKSTATUSPTR( status );

  /* ---------- Compute Net Strain -------------*/

  /* for good measure we high pass filter both residual and control signals again
     before adding them together */
  LALButterworthREAL8TimeSeries(status->statusPtr,&hR,&highpassfilterpar);
  CHECKSTATUSPTR( status );
  LALButterworthREAL8TimeSeries(status->statusPtr,&hC,&highpassfilterpar);
  CHECKSTATUSPTR( status );

  /* now add control and residual signals together and we're done */
  for (p=0; p < (int)h.data->length; p++) {
    h.data->data[p]= hR.data->data[p]+ hC.data->data[p];
    output->h.data->data[p]=h.data->data[p];
  }

  /* destroy vectors that hold the data */
  LALDDestroyVector(status->statusPtr,&h.data);
  CHECKSTATUSPTR( status );
  LALDDestroyVector(status->statusPtr,&hR.data);
  CHECKSTATUSPTR( status );
  LALDDestroyVector(status->statusPtr,&uphR.data);
  CHECKSTATUSPTR( status );
  LALDDestroyVector(status->statusPtr,&hC.data);
  CHECKSTATUSPTR( status );
  LALDDestroyVector(status->statusPtr,&hCX.data);
  CHECKSTATUSPTR( status );
  LALDDestroyVector(status->statusPtr,&hCY.data);
  CHECKSTATUSPTR( status );

  /* Destroy vectors that hold filter coefficients */
  LALDDestroyVector(status->statusPtr,&Cinv.directCoef);
  CHECKSTATUSPTR( status );
  LALDDestroyVector(status->statusPtr,&Cinv.recursCoef);
  CHECKSTATUSPTR( status );
  LALDDestroyVector(status->statusPtr,&Cinv.history);
  CHECKSTATUSPTR( status );

  for(p=0;p<NGfilt;p++){
    LALDDestroyVector(status->statusPtr,&G[p].directCoef);
    CHECKSTATUSPTR( status );
    LALDDestroyVector(status->statusPtr,&G[p].recursCoef);
    CHECKSTATUSPTR( status );
    LALDDestroyVector(status->statusPtr,&G[p].history);   
    CHECKSTATUSPTR( status );
  }

  LALDDestroyVector(status->statusPtr,&AA.directCoef);
  CHECKSTATUSPTR( status );
  LALDDestroyVector(status->statusPtr,&AA.recursCoef);
  CHECKSTATUSPTR( status );
  LALDDestroyVector(status->statusPtr,&AA.history);
  CHECKSTATUSPTR( status );

  for(p=0;p<NAXfilt;p++){
    LALDDestroyVector(status->statusPtr,&AX[p].directCoef);
    CHECKSTATUSPTR( status );
    LALDDestroyVector(status->statusPtr,&AX[p].recursCoef);
    CHECKSTATUSPTR( status );
    LALDDestroyVector(status->statusPtr,&AX[p].history);
    CHECKSTATUSPTR( status );
  }

  for(p=0;p<NAYfilt;p++){
    LALDDestroyVector(status->statusPtr,&AY[p].directCoef);
    CHECKSTATUSPTR( status );
    LALDDestroyVector(status->statusPtr,&AY[p].recursCoef);
    CHECKSTATUSPTR( status );
    LALDDestroyVector(status->statusPtr,&AY[p].history);
    CHECKSTATUSPTR( status );
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}

/*******************************************************************************/


int XLALhCTimesBeta(REAL8TimeSeries *hC, StrainOut *output)
{
  int n,m;
  REAL8 time,InterpolatedBeta;
  static REAL8 beta[MAXALPHAS],tainterp[MAXALPHAS];

  /* copy ouput betas into local array */
  for(m=0; m < (int)output->beta.data->length; m++)
    {
      beta[m]=output->beta.data->data[m].re;
      tainterp[m]= m*output->beta.deltaT;
    }

  time=0.0;    /* time variable */

  {
    gsl_interp_accel *acc_alpha = gsl_interp_accel_alloc();      /* GSL spline interpolation stuff */
    gsl_spline *spline_alpha = gsl_spline_alloc(gsl_interp_cspline,output->beta.data->length);
    gsl_spline_init(spline_alpha,tainterp,beta,output->beta.data->length);

    for (n = 0; n < (int)hC->data->length; n++) 
      {
	InterpolatedBeta=gsl_spline_eval(spline_alpha,time,acc_alpha); 
	hC->data->data[n] *= InterpolatedBeta;
	time=time+hC->deltaT;
      }

    /* clean up GSL spline interpolation stuff */
    gsl_spline_free(spline_alpha);
    gsl_interp_accel_free(acc_alpha);
  }

  return 0;
}


/*******************************************************************************/


int XLALUpsamplehR(REAL8TimeSeries *uphR, REAL8TimeSeries *hR, int up_factor)
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


int XLALhROverAlpha(REAL8TimeSeries *hR, StrainOut *output)
{
  int n,m;
  double time,InterpolatedAlpha;
  double alpha[MAXALPHAS],tainterp[MAXALPHAS];

  /* copy ouput alphas into local array */
  for(m=0; m < (int)output->alpha.data->length; m++)
    {
      alpha[m]=output->alpha.data->data[m].re;
      tainterp[m]= m*output->alpha.deltaT;
    }

  time=0.0;    /* time variable */
  
  {
    gsl_interp_accel *acc_alpha = gsl_interp_accel_alloc();      /* GSL spline interpolation stuff */
    gsl_spline *spline_alpha = gsl_spline_alloc(gsl_interp_cspline,output->alpha.data->length);
    gsl_spline_init(spline_alpha,tainterp,alpha,output->alpha.data->length);

    for (n = 0; n < (int)hR->data->length; n++) 
      {
	InterpolatedAlpha=gsl_spline_eval(spline_alpha,time,acc_alpha);
	
	hR->data->data[n] /= InterpolatedAlpha;
	time=time+hR->deltaT;
      }
  /* clean up GSL spline interpolation stuff */
  gsl_spline_free(spline_alpha);
  gsl_interp_accel_free(acc_alpha);
  }

  return 0;
}

/*******************************************************************************/



void XLALMakeFilters(LALStatus *status, REAL8IIRFilter *Cinv, REAL8IIRFilter *G,REAL8IIRFilter 
  *AA,REAL8IIRFilter *AX,REAL8IIRFilter *AY)
{
  int l,n;
  
  LALDCreateVector(status->statusPtr,&Cinv->directCoef,CinvDirectOrder);
  LALDCreateVector(status->statusPtr,&Cinv->recursCoef,CinvRecursOrder);
  LALDCreateVector(status->statusPtr,&Cinv->history,CinvDirectOrder-1);

  for(l=0;l<CinvDirectOrder;l++) Cinv->directCoef->data[l]=CinvDirectCoefs[l];
  for(l=0;l<CinvDirectOrder;l++) Cinv->recursCoef->data[l]=-CinvRecursCoefs[l];
  for(l=0;l<CinvDirectOrder-1;l++) Cinv->history->data[l]=0.0;

  for(n=0;n<NGfilt;n++){

    LALDCreateVector(status->statusPtr,&G[n].directCoef,G_Dord);
    LALDCreateVector(status->statusPtr,&G[n].recursCoef,G_Dord);
    LALDCreateVector(status->statusPtr,&G[n].history,G_Dord-1);
   
    /* Fill coefficient vectors with coefficients from filters.h */
    for(l=0;l<G_Dord;l++) G[n].directCoef->data[l]=G_D[n][l];
    for(l=0;l<G_Dord;l++) G[n].recursCoef->data[l]=-G_R[n][l];
    for(l=0;l<G_Dord-1;l++) G[n].history->data[l]=0.0;

  }

  LALDCreateVector(status->statusPtr,&AA->directCoef, A_0_Rord);
  LALDCreateVector(status->statusPtr,&AA->recursCoef, A_0_Rord);
  LALDCreateVector(status->statusPtr,&AA->history, A_0_Rord-1);

  /* Fill coefficient vectors with coefficients from filters.h */
  for(l=0;l< A_0_Rord;l++) AA->directCoef->data[l]=A_0_D[l];
  for(l=0;l< A_0_Rord;l++) AA->recursCoef->data[l]=-A_0_R[l];
  for(l=0;l< A_0_Rord-1;l++) AA->history->data[l]=0.0;

  for(n=0;n<NAXfilt;n++){
   
    LALDCreateVector(status->statusPtr,&AX[n].directCoef, A_digital_Rord);
    LALDCreateVector(status->statusPtr,&AX[n].recursCoef, A_digital_Rord);
    LALDCreateVector(status->statusPtr,&AX[n].history, A_digital_Rord-1);

    for(l=0;l< A_digital_Rord;l++) AX[n].directCoef->data[l]=AX_D[n][l];
    for(l=0;l< A_digital_Rord;l++) AX[n].recursCoef->data[l]=-AX_R[n][l];
    for(l=0;l< A_digital_Rord-1;l++) AX[n].history->data[l]=0.0;
    
  }

  for(n=0;n<NAYfilt;n++){
    LALDCreateVector(status->statusPtr,&AY[n].directCoef, A_digital_Rord);
    LALDCreateVector(status->statusPtr,&AY[n].recursCoef, A_digital_Rord);
    LALDCreateVector(status->statusPtr,&AY[n].history, A_digital_Rord-1);

    for(l=0;l< A_digital_Rord;l++) AY[n].directCoef->data[l]=AY_D[n][l];
    for(l=0;l< A_digital_Rord;l++) AY[n].recursCoef->data[l]=-AY_R[n][l];
    for(l=0;l< A_digital_Rord-1;l++) AY[n].history->data[l]=0.0;
  }

  RETURN (status);
}


/*******************************************************************************/

void XLALGetFactors(LALStatus *status, StrainOut *output, StrainIn *input)
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

  /* Create local data vectors */
  LALCreateVector(status->statusPtr,&asq.data,(UINT4)(To/input->AS_Q.deltaT +0.5));
  LALCreateVector(status->statusPtr,&darm.data,(UINT4)(To/input->DARM.deltaT +0.5));
  LALCreateVector(status->statusPtr,&exc.data,(UINT4)(To/input->EXC.deltaT +0.5));

  /* Create Window vectors */
  LALCreateVector(status->statusPtr,&asqwin,(UINT4)(To/input->AS_Q.deltaT +0.5));
  LALCreateVector(status->statusPtr,&darmwin,(UINT4)(To/input->DARM.deltaT +0.5));
  LALCreateVector(status->statusPtr,&excwin,(UINT4)(To/input->EXC.deltaT +0.5));

  /* assign time spacing for local time series */
  asq.deltaT=input->AS_Q.deltaT;
  darm.deltaT=input->DARM.deltaT;
  exc.deltaT=input->EXC.deltaT;

  winparams.type=Hann;
   
  /* windows for time domain channels */
  /* asq */
  winparams.length=(INT4)(To/asq.deltaT +0.5);
  LALWindow(status->statusPtr,asqwin,&winparams);
  
  /* darm */
  winparams.length=(INT4)(To/darm.deltaT +0.5);
  LALWindow(status->statusPtr,darmwin,&winparams);

  /* exc */
  winparams.length=(INT4)(To/exc.deltaT +0.5);
  LALWindow(status->statusPtr,excwin,&winparams);

  for(m=0; m < (UINT4)(deltaT*length) / To; m++)
    {

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

      output->alpha.data->data[m]= factors.alpha;
      output->beta.data->data[m]= factors.beta;

      if(m == MAXALPHAS)
	{
	  fprintf(stderr,"Too many values of the factors, maximum allowed is %d\n",MAXALPHAS);
	  RETURN(status);
	}

    }

  /* Clean up */
  LALDestroyVector(status->statusPtr,&darm.data);
  LALDestroyVector(status->statusPtr,&exc.data);
  LALDestroyVector(status->statusPtr,&asq.data);

  LALDestroyVector(status->statusPtr,&asqwin);
  LALDestroyVector(status->statusPtr,&darmwin);
  LALDestroyVector(status->statusPtr,&excwin);


  RETURN (status);
}

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

NRCSID( UPDATEFACTORSC, "$Id$" );

void LALComputeStrain(
    LALStatus              *status,
    StrainOut              *output,    
    StrainIn               *input
    )
{
/* Inverse sensing, servo, analog actuation, digital x 
actuation  digital y actuation */
static MyIIRFilter Cinv,G[NGfilt],AA,AX[NAXfilt],AY[NAYfilt];    
static REAL8TimeSeries h,hR,uphR,hC,hCX,hCY;
int p;
PassBandParamStruc highpassfilterpar,lowpassfilterpar;

  /* high pass filter parameters */ 
  highpassfilterpar.name  = "Butterworth High Pass";
  highpassfilterpar.nMax  = 10;
  highpassfilterpar.f2    = 40.0;
  highpassfilterpar.a2    = 0.5;
  highpassfilterpar.f1    = -1.0;
  highpassfilterpar.a1    = -1.0;
  /* low pass filter parameters */
  lowpassfilterpar.name  = "Butterworth Low Pass";
  lowpassfilterpar.nMax  = 12;
  lowpassfilterpar.f2    = -1.0;
  lowpassfilterpar.a2    = -1.0;
  lowpassfilterpar.f1    = 6000.0;
  lowpassfilterpar.a1    = 0.5;

  if(XLALMakeFilters(&Cinv,&G,&AA,&AX,&AY)) RETURN( status );
  if(XLALGetFactors(status, output, input)) RETURN( status );
  
  /* Create vectors that will hold the residual, control and net strain signals */
  LALDCreateVector(status,&h.data,input->AS_Q.data->length);
  LALDCreateVector(status,&hR.data,input->AS_Q.data->length);
  LALDCreateVector(status,&hC.data,input->AS_Q.data->length);
  LALDCreateVector(status,&hCX.data,input->AS_Q.data->length);
  LALDCreateVector(status,&hCY.data,input->AS_Q.data->length);
  LALDCreateVector(status,&uphR.data,UpSamplingFactor*input->AS_Q.data->length);
  h.deltaT=input->AS_Q.deltaT;
  hR.deltaT=input->AS_Q.deltaT;
  hC.deltaT=input->AS_Q.deltaT;
  uphR.deltaT=input->AS_Q.deltaT/UpSamplingFactor;

  /* copy AS_Q input into residual strain as double */  
  for (p=0; p<hR.data->length; p++) {
    hR.data->data[p]=input->AS_Q.data->data[p];
  }
  /* copy AS_Q input into control strain as double */  
  for (p=0; p<hC.data->length; p++) {
    hC.data->data[p]=input->AS_Q.data->data[p];
  }

  /* high pass filter both time series */
  LALButterworthREAL8TimeSeries(status,&hR,&highpassfilterpar);
  LALButterworthREAL8TimeSeries(status,&hC,&highpassfilterpar);

  /* ---------- Compute Residual Strain -------------*/
  /* to get the residual strain we must first divide AS_Q by alpha */
  if(XLALhROverAlpha(&hR, output)) RETURN( status );

  /* then we upsample (and smooth it with a low pass filter) */
  if(XLALUpsamplehR(&uphR, &hR, UpSamplingFactor)) RETURN( status );
  LALButterworthREAL8TimeSeries(status,&uphR,&lowpassfilterpar);

  /* then we filter through the inverse of the sensing function */
  if(XLALFilterSeries(&Cinv, &uphR)) RETURN( status );

  /* Low pass again before downsampling */
  LALButterworthREAL8TimeSeries(status,&uphR,&lowpassfilterpar);
 
  /* then we downsample and voila' */
  for (p=0; p<hR.data->length; p++) {
    hR.data->data[p]=uphR.data->data[p*UpSamplingFactor];
  }

  /* ---------- Compute Control Strain -------------*/
  /* to get the control strain we first multiply AS_Q by beta */
  XLALhCTimesBeta(&hC, output);
  
  /* Now we filter through the servo */
  for(p=NGfilt-1;p>=0;p--){
    if (XLALFilterSeries(&G[p],&hC)) RETURN( status );
  }
  /* and adjust to account for servo gain to get darm_ctrl */ 
  for (p=0; p<hC.data->length;p++) {
    hC.data->data[p] *= ServoGain;
  }

  /* Copy data into x and y time series for parallel filtering */
  for (p=0; p < hCX.data->length; p++) {
    hCX.data->data[p]=hC.data->data[p];
  }
  for (p=0; p < hCY.data->length; p++) {
    hCY.data->data[p]=hC.data->data[p];
  }
  
  /* Filter x-arm */
  for(p = NAXfilt-1; p >= 0; p--){
    if (XLALFilterSeries(&AX[p],&hCX)) RETURN( status );
  }
  /* Adjust to account for digital gain on x-arm*/ 
  for (p=0; p< hC.data->length;p++) {
    hCX.data->data[p] *= AXGain;
  }
  
  /* Filter y-arm */
  for(p = NAYfilt-1; p >= 0; p--){
    if (XLALFilterSeries(&AY[p],&hCY)) return 5;
  }
  /* Adjust to account for digital gain on y-arm*/ 
  for (p = 0; p < hC.data->length; p++) {
    hCY.data->data[p] *= AYGain;
  }

  /* add x-arm and y-arm together */
  for (p=0; p< hC.data->length; p++) {
    hC.data->data[p]=(hCX.data->data[p]+hCY.data->data[p])/2;
  }

  /* filter through analog part of actuation and voila' */ 
  if (XLALFilterSeries(&AA,&hC)) RETURN( status );


  /* ---------- Compute Net Strain -------------*/

  /* for good measure we high pass filter both residual and control signals again
     before adding them together */
  LALButterworthREAL8TimeSeries(status,&hR,&highpassfilterpar);
  LALButterworthREAL8TimeSeries(status,&hC,&highpassfilterpar);

  /* now add control and residual signals together and we're done */
  for (p=0; p < h.data->length; p++) {
    h.data->data[p]= hR.data->data[p]+ hC.data->data[p];
    output->h.data->data[p]=h.data->data[p];
  }

  /* destroy vectors that hold the data */
  LALDDestroyVector(status,&h.data);
  LALDDestroyVector(status,&hR.data);
  LALDDestroyVector(status,&uphR.data);
  LALDDestroyVector(status,&hC.data);
  LALDDestroyVector(status,&hCX.data);
  LALDDestroyVector(status,&hCY.data);

  RETURN( status );
}


/*******************************************************************************/


int XLALhCTimesBeta(REAL8TimeSeries *hC, StrainOut *output)
{
  int n,m;
  REAL8 time,InterpolatedBeta;
  static REAL8 beta[MAXALPHAS],tainterp[MAXALPHAS];

  /* copy ouput betas into local array */
  for(m=0; m < output->beta.data->length; m++)
    {
      beta[m]=output->beta.data->data[m].re;
      tainterp[m]= m*output->beta.deltaT;
    }
  
  gsl_interp_accel *acc_alpha = gsl_interp_accel_alloc();      /* GSL spline interpolation stuff */
  gsl_spline *spline_alpha = gsl_spline_alloc(gsl_interp_cspline,output->beta.data->length);
  gsl_spline_init(spline_alpha,tainterp,beta,output->beta.data->length);

  time=0.0;    /* time variable */

  for (n = 0; n < hC->data->length; n++) {

    InterpolatedBeta=gsl_spline_eval(spline_alpha,time,acc_alpha);
    
    hC->data->data[n] *= InterpolatedBeta;
    time=time+hC->deltaT;
  }

  /* clean up GSL spline interpolation stuff */
  gsl_spline_free(spline_alpha);
  gsl_interp_accel_free(acc_alpha);

  return 0;
}



/*******************************************************************************/

int XLALFilterSeries(MyIIRFilter *F, REAL8TimeSeries *TSeries)
{
  int n,r;
  REAL8 yn,xn,xsum,ysum;

  for (n=0; n<TSeries->data->length;n++) 
    {
      xsum=0.0;
      ysum=0.0;

      xn=TSeries->data->data[n];
    
      for(r=0;r<F->xOrder-1;r++)
	{
	  xsum += F->xhist[r]*F->b[r+1];
	}
      xsum=xsum+xn*F->b[0];
    
      for(r=0;r<F->yOrder-1;r++)
	{
	  ysum -= F->yhist[r]*F->a[r+1];
	}
    
      yn=xsum+ysum;

      TSeries->data->data[n]=yn;

      for(r=F->xOrder-2;r>0;r--)
	{
	  F->xhist[r]=F->xhist[r-1];
	}
      for(r=F->yOrder-2;r>0;r--)
	{
	  F->yhist[r]=F->yhist[r-1];
	}

      F->yhist[0]=yn;
      F->xhist[0]=xn;
    }

  return 0;
}

/*******************************************************************************/


int XLALUpsamplehR(REAL8TimeSeries *uphR, REAL8TimeSeries *hR, int up_factor)
{
  int n;

  /* Set all values to 0 */
  for (n=0; n < uphR->data->length; n++) {
    uphR->data->data[n] = 0.0;
  }

  /* Set one in every up_factor to the value of hR x USR */
  for (n=0; n < hR->data->length; n++) {
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
  for(m=0; m < output->alpha.data->length; m++)
    {
      alpha[m]=output->alpha.data->data[m].re;
      tainterp[m]= m*output->alpha.deltaT;
    }
  
  gsl_interp_accel *acc_alpha = gsl_interp_accel_alloc();      /* GSL spline interpolation stuff */
  gsl_spline *spline_alpha = gsl_spline_alloc(gsl_interp_cspline,output->alpha.data->length);
  gsl_spline_init(spline_alpha,tainterp,alpha,output->alpha.data->length);

  time=0.0;    /* time variable */

  for (n = 0; n < hR->data->length; n++) {

    InterpolatedAlpha=gsl_spline_eval(spline_alpha,time,acc_alpha);
    
    hR->data->data[n] /= InterpolatedAlpha;
    time=time+hR->deltaT;

/*     fprintf(stdout,"%e %1.17e\n", time,InterpolatedAlpha); */


  }

  /* clean up GSL spline interpolation stuff */
  gsl_spline_free(spline_alpha);
  gsl_interp_accel_free(acc_alpha);

  return 0;
}

/*******************************************************************************/

int XLALMakeFilters(MyIIRFilter *Cinv,MyIIRFilter *G,MyIIRFilter *AA,MyIIRFilter *AX,MyIIRFilter *AY)
{
  int l,n;
 
  Cinv->yOrder=CinvRecursOrder;
  Cinv->xOrder=CinvDirectOrder;

  /* Fill coefficient vectors with coefficients from filters.h */
  for(l=0;l<Cinv->xOrder;l++) Cinv->b[l]=CinvDirectCoefs[l];
  for(l=0;l<Cinv->yOrder;l++) Cinv->a[l]=CinvRecursCoefs[l];
  for(l=0;l<Cinv->yOrder-1;l++) Cinv->yhist[l]=0.0;
  for(l=0;l<Cinv->xOrder-1;l++) Cinv->xhist[l]=0.0;

  for(n=0;n<NGfilt;n++){
    G[n].xOrder=G_Dord;
    G[n].yOrder=G_Rord;
    
    /* Fill coefficient vectors with coefficients from filters.h */
    for(l=0;l<G[n].xOrder;l++) G[n].b[l]=G_D[n][l];
    for(l=0;l<G[n].yOrder;l++) G[n].a[l]=G_R[n][l];
    for(l=0;l<G[n].yOrder-1;l++) G[n].yhist[l]=0.0;
    for(l=0;l<G[n].xOrder-1;l++) G[n].xhist[l]=0.0;
  }

    AA->yOrder= A_0_Rord;
    AA->xOrder= A_0_Dord;

    /* Fill coefficient vectors with coefficients from filters.h */
    for(l=0;l<AA->xOrder;l++) AA->b[l]=A_0_D[l];
    for(l=0;l<AA->yOrder;l++) AA->a[l]=A_0_R[l];
    for(l=0;l<AA->yOrder-1;l++) AA->yhist[l]=0.0;
    for(l=0;l<AA->xOrder-1;l++) AA->xhist[l]=0.0;


  for(n=0;n<NAXfilt;n++){
    AX[n].yOrder=A_digital_Rord;
    AX[n].xOrder=A_digital_Dord;
    
    /* Fill coefficient vectors with coefficients from filters.h */
    for(l=0;l<AX[n].xOrder;l++) AX[n].b[l]=AX_D[n][l];
    for(l=0;l<AX[n].yOrder;l++) AX[n].a[l]=AX_R[n][l];
    for(l=0;l<AX[n].yOrder-1;l++) AX[n].yhist[l]=0.0;
    for(l=0;l<AX[n].xOrder-1;l++) AX[n].xhist[l]=0.0;
  }

  for(n=0;n<NAYfilt;n++){
    AY[n].yOrder=A_digital_Rord;
    AY[n].xOrder=A_digital_Dord;
    
    /* Fill coefficient vectors with coefficients from filters.h */
    for(l=0;l<AY[n].xOrder;l++) AY[n].b[l]=AY_D[n][l];
    for(l=0;l<AY[n].yOrder;l++) AY[n].a[l]=AY_R[n][l];
    for(l=0;l<AY[n].yOrder-1;l++) AY[n].yhist[l]=0.0;
    for(l=0;l<AY[n].xOrder-1;l++) AY[n].xhist[l]=0.0;
  }

  return 0;
}

/*******************************************************************************/



int XLALGetFactors(LALStatus *status, StrainOut *output, StrainIn *input)
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
  LALCreateVector(status,&asq.data,(UINT4)(To/input->AS_Q.deltaT +0.5));
  LALCreateVector(status,&darm.data,(UINT4)(To/input->DARM.deltaT +0.5));
  LALCreateVector(status,&exc.data,(UINT4)(To/input->EXC.deltaT +0.5));

  /* Create Window vectors */
  LALCreateVector(status,&asqwin,(UINT4)(To/input->AS_Q.deltaT +0.5));
  LALCreateVector(status,&darmwin,(UINT4)(To/input->DARM.deltaT +0.5));
  LALCreateVector(status,&excwin,(UINT4)(To/input->EXC.deltaT +0.5));

  /* assign time spacing for local time series */
  asq.deltaT=input->AS_Q.deltaT;
  darm.deltaT=input->DARM.deltaT;
  exc.deltaT=input->EXC.deltaT;

  winparams.type=Hann;
   
  /* windows for time domain channels */
  /* asq */
  winparams.length=(INT4)(To/asq.deltaT +0.5);
  LALWindow(status,asqwin,&winparams);
  
  /* darm */
  winparams.length=(INT4)(To/darm.deltaT +0.5);
  LALWindow(status,darmwin,&winparams);

  /* exc */
  winparams.length=(INT4)(To/exc.deltaT +0.5);
  LALWindow(status,excwin,&winparams);

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

      LALComputeCalibrationFactors(status,&factors,&params);

      output->alpha.data->data[m]= factors.alpha;
      output->beta.data->data[m]= factors.beta;

/*       fprintf(stdout,"%e %e %e %e\n",output->alpha.data->data[m].re,output->alpha.data->data[m].im, */
/* 	                             output->beta.data->data[m].re,output->beta.data->data[m].im); */



      if(m == MAXALPHAS)
	{
	  fprintf(stderr,"Too many values of the factors, maximum allowed is %d\n",MAXALPHAS);
	  return 1;
	}
    }

  /* Clean up */
  LALDestroyVector(status,&darm.data);
  LALDestroyVector(status,&exc.data);
  LALDestroyVector(status,&asq.data);

  LALDestroyVector(status,&asqwin);
  LALDestroyVector(status,&darmwin);
  LALDestroyVector(status,&excwin);

  return 0;
}

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

/* Butterworth high pass at 40Hz */
const REAL8 b_high[3]={9.892117314337300e-01,    -1.978423462867460e+00,     9.892117314337300e-01};
const REAL8 a_high[3]={1.000000000000000e+00,    -1.978307072742137e+00,     9.785398529927832e-01};

/* Butterworth low pass at 6000Hz */
const REAL8 b_low[3]={4.686543159618028e-03,     9.373086319236057e-03,     4.686543159618028e-03};
const REAL8 a_low[3]={1.000000000000000e+00,    -1.797224507982706e+00,     8.159706806211777e-01};

NRCSID( COMPUTESTRAINC, "$Id$" );

void LALComputeStrain(
    LALStatus              *status,
    StrainOut              *output,    
    StrainIn               *input
    )
{
/* Inverse sensing, servo, analog actuation, digital x 
actuation  digital y actuation */
static MyIIRFilter Cinv,G[NGfilt],AA,AX[NAXfilt],AY[NAYfilt];    
static REAL8IIRFilter CinvLAL,GLAL[NGfilt],AALAL,AXLAL[NAXfilt],AYLAL[NAYfilt];    
static REAL8TimeSeries h,hR,uphR,hC,hCX,hCY;
int p;
PassBandParamStruc highpassfilterpar,lowpassfilterpar;

 INITSTATUS( status, "LALComputeStrain", COMPUTESTRAINC );
 ATTATCHSTATUSPTR( status );


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

  if(XLALMakeFilters(&Cinv,&G,&AA,&AX,&AY)) ABORT(status,111,"Broke making iir filters");
  if(XLALMakeFilters2(status,&CinvLAL,&GLAL,&AALAL,&AXLAL,&AYLAL)) ABORT(status,112,"Broke making LAL iir filters");
  if(XLALGetFactors(status, output, input)) ABORT(status,113,"Broke making factors");
  
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
  for (p=0; p<hR.data->length; p++) {
    hR.data->data[p]=input->AS_Q.data->data[p];
  }
  /* copy AS_Q input into control strain as double */  
  for (p=0; p<hC.data->length; p++) {
    hC.data->data[p]=input->AS_Q.data->data[p];
  }

  /* high pass filter both time series */
  LALButterworthREAL8TimeSeries(status->statusPtr,&hR,&highpassfilterpar);
  CHECKSTATUSPTR( status );
  LALButterworthREAL8TimeSeries(status->statusPtr,&hC,&highpassfilterpar);
  CHECKSTATUSPTR( status );
  
/*   if(XLALHighPass(&hR)) ABORT(status,113,"Broke highpassing hR"); */
/*   if(XLALHighPass(&hC)) ABORT(status,115,"Broke highpassing hC"); */

  /* ---------- Compute Residual Strain -------------*/
  /* to get the residual strain we must first divide AS_Q by alpha */
  if(XLALhROverAlpha(&hR, output)) ABORT(status,116,"Broke at hR/alpha") ;

  /* then we upsample (and smooth it with a low pass filter) */
  if(XLALUpsamplehR(&uphR, &hR, UpSamplingFactor)) ABORT(status,117,"Broke upsampling hR");
  LALButterworthREAL8TimeSeries(status->statusPtr,&uphR,&lowpassfilterpar);
  CHECKSTATUSPTR( status );

/*   if(XLALLowPass(&uphR)) ABORT(status,118,"Broke lowpassing uphR"); */

  /* then we filter through the inverse of the sensing function */
/*   if(XLALFilterSeries(&Cinv, &uphR)) ABORT(status,119,"Broke filtering hR through Cinv"); */

  LALIIRFilterREAL8Vector(status->statusPtr,uphR.data,&CinvLAL);
  CHECKSTATUSPTR( status );

  /* Low pass again before downsampling */
  LALButterworthREAL8TimeSeries(status->statusPtr,&uphR,&lowpassfilterpar);
  CHECKSTATUSPTR( status );
/*   if(XLALLowPass(&uphR)) ABORT(status,120,"Broke lowpassing uphR"); */
 
  /* then we downsample and voila' */
  for (p=0; p<hR.data->length; p++) {
    hR.data->data[p]=uphR.data->data[p*UpSamplingFactor];
  }

  /* ---------- Compute Control Strain -------------*/
  /* to get the control strain we first multiply AS_Q by beta */
  if( XLALhCTimesBeta(&hC, output)) ABORT(status,120,"Broke in hC x beta");
  
  /* Now we filter through the servo */
  for(p=NGfilt-1;p>=0;p--){
/*     if (XLALFilterSeries(&G[p],&hC)) ABORT(status,121,"Broke filtering hC through servo"); */
    LALIIRFilterREAL8Vector(status->statusPtr,hC.data,&GLAL[p]);
    CHECKSTATUSPTR( status );
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
/*     if (XLALFilterSeries(&AX[p],&hCX)) ABORT(status,122,"Broke filtering hC through AX"); */
    LALIIRFilterREAL8Vector(status->statusPtr,hCX.data,&AXLAL[p]);
    CHECKSTATUSPTR( status );
  }
  /* Adjust to account for digital gain on x-arm*/ 
  for (p=0; p< hC.data->length;p++) {
    hCX.data->data[p] *= AXGain;
  }
  
  /* Filter y-arm */
  for(p = NAYfilt-1; p >= 0; p--){
/*     if (XLALFilterSeries(&AY[p],&hCY)) ABORT(status,123,"Broke filtering hC through AY"); */
    LALIIRFilterREAL8Vector(status->statusPtr,hCY.data,&AYLAL[p]);
    CHECKSTATUSPTR( status );
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
/*   if (XLALFilterSeries(&AA,&hC)) ABORT(status,124,"Broke filtering hC through AA"); */
    LALIIRFilterREAL8Vector(status->statusPtr,hC.data,&AALAL);
    CHECKSTATUSPTR( status );


  /* ---------- Compute Net Strain -------------*/

  /* for good measure we high pass filter both residual and control signals again
     before adding them together */
  LALButterworthREAL8TimeSeries(status->statusPtr,&hR,&highpassfilterpar);
  CHECKSTATUSPTR( status );
  LALButterworthREAL8TimeSeries(status->statusPtr,&hC,&highpassfilterpar);
  CHECKSTATUSPTR( status );

/*   if (XLALHighPass(&hR))ABORT(status,125,"Broke in final highpassing of hR"); */
/*   if (XLALHighPass(&hC))ABORT(status,126,"Broke in final highpassing of hC"); */

  /* now add control and residual signals together and we're done */
  for (p=0; p < h.data->length; p++) {
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

  LALDDestroyVector(status->statusPtr,&CinvLAL.directCoef);
  CHECKSTATUSPTR( status );
  LALDDestroyVector(status->statusPtr,&CinvLAL.recursCoef);
  CHECKSTATUSPTR( status );
  LALDDestroyVector(status->statusPtr,&CinvLAL.history);
  CHECKSTATUSPTR( status );

  for(p=0;p<NGfilt;p++){
    LALDDestroyVector(status->statusPtr,&GLAL[p].directCoef);
    CHECKSTATUSPTR( status );
    LALDDestroyVector(status->statusPtr,&GLAL[p].recursCoef);
    CHECKSTATUSPTR( status );
    LALDDestroyVector(status->statusPtr,&GLAL[p].history);   
    CHECKSTATUSPTR( status );
  }

  LALDDestroyVector(status->statusPtr,&AALAL.directCoef);
  CHECKSTATUSPTR( status );
  LALDDestroyVector(status->statusPtr,&AALAL.recursCoef);
  CHECKSTATUSPTR( status );
  LALDDestroyVector(status->statusPtr,&AALAL.history);
  CHECKSTATUSPTR( status );

  for(p=0;p<NAXfilt;p++){
    LALDDestroyVector(status->statusPtr,&AXLAL[p].directCoef);
    CHECKSTATUSPTR( status );
    LALDDestroyVector(status->statusPtr,&AXLAL[p].recursCoef);
    CHECKSTATUSPTR( status );
    LALDDestroyVector(status->statusPtr,&AXLAL[p].history);
    CHECKSTATUSPTR( status );
  }

  for(p=0;p<NAYfilt;p++){
    LALDDestroyVector(status->statusPtr,&AYLAL[p].directCoef);
    CHECKSTATUSPTR( status );
    LALDDestroyVector(status->statusPtr,&AYLAL[p].recursCoef);
    CHECKSTATUSPTR( status );
    LALDDestroyVector(status->statusPtr,&AYLAL[p].history);
    CHECKSTATUSPTR( status );
  }



  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/*******************************************************************************/


int XLALLowPass(REAL8TimeSeries *TSeries)
{
  MyIIRFilter ButterLow;
  int p,l;

  ButterLow.yOrder=3;
  ButterLow.xOrder=3;

  /* Fill coefficient vectors with coefficients */
  for(l=0;l<ButterLow.xOrder;l++) ButterLow.b[l]=b_low[l];
  for(l=0;l<ButterLow.yOrder;l++) ButterLow.a[l]=a_low[l];

  for (p=0;p<3;p++)
    {
      /* Set history to zero */
      for(l=0;l<ButterLow.yOrder-1;l++) ButterLow.yhist[l]=0.0;
      for(l=0;l<ButterLow.xOrder-1;l++) ButterLow.xhist[l]=0.0;

      /* Filter forward */
      XLALFilterSeries(&ButterLow,TSeries);   

      /* Set history to zero */
      for(l=0;l<ButterLow.yOrder-1;l++) ButterLow.yhist[l]=0.0;
      for(l=0;l<ButterLow.xOrder-1;l++) ButterLow.xhist[l]=0.0;

      /* filter backwards */
      XLALFilterSeriesReverse(&ButterLow,TSeries);   

    }
  

  return 0;
}

/*******************************************************************************/


int XLALHighPass(REAL8TimeSeries *TSeries)
{
  MyIIRFilter ButterHigh;
  int p,l;

  ButterHigh.yOrder=3;
  ButterHigh.xOrder=3;

  /* Fill coefficient vectors with coefficients */
  for(l=0;l<ButterHigh.xOrder;l++) ButterHigh.b[l]=b_high[l];
  for(l=0;l<ButterHigh.yOrder;l++) ButterHigh.a[l]=a_high[l];

  for (p=0;p<2;p++)
    {
      /* Set history to zero */
      for(l=0;l<ButterHigh.yOrder-1;l++) ButterHigh.yhist[l]=0.0;
      for(l=0;l<ButterHigh.xOrder-1;l++) ButterHigh.xhist[l]=0.0;

      /* Filter forward */
      XLALFilterSeries(&ButterHigh,TSeries);   

      /* Set history to zero */
      for(l=0;l<ButterHigh.yOrder-1;l++) ButterHigh.yhist[l]=0.0;
      for(l=0;l<ButterHigh.xOrder-1;l++) ButterHigh.xhist[l]=0.0;

      /* filter backwards */
      XLALFilterSeriesReverse(&ButterHigh,TSeries);   

    }
  

  return 0;
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

int XLALFilterSeriesReverse(MyIIRFilter *F, REAL8TimeSeries *TSeries)
{
  int n,r;
  long double yn,xn,xsum,ysum;

  for (n=TSeries->data->length-1; n >= 0; n--) 
    {
      xsum=0.0;
      ysum=0.0;

      xn=TSeries->data->data[n];
    
      for(r=0;r<F->xOrder-1;r++)
	{
	  xsum += (long double) F->xhist[r]*  (long double) F->b[r+1];
	}
      xsum=xsum+xn*F->b[0];
    
      for(r=0;r<F->yOrder-1;r++)
	{
	  ysum -=  (long double) F->yhist[r] *  (long double) F->a[r+1];
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

int XLALFilterSeries(MyIIRFilter *F, REAL8TimeSeries *TSeries)
{
  int n,r;
  long double yn,xn,xsum,ysum;

  for (n=0; n<TSeries->data->length;n++) 
    {
      xsum=0.0;
      ysum=0.0;

      xn=TSeries->data->data[n];
    
      for(r=0;r<F->xOrder-1;r++)
	{
	  xsum += (long double) F->xhist[r]*  (long double) F->b[r+1];
	}
      xsum=xsum+xn*F->b[0];
    
      for(r=0;r<F->yOrder-1;r++)
	{
	  ysum -=  (long double) F->yhist[r] *  (long double) F->a[r+1];
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



int XLALMakeFilters2(LALStatus *status, REAL8IIRFilter *CinvLAL, REAL8IIRFilter *GLAL,REAL8IIRFilter 
  *AALAL,REAL8IIRFilter *AXLAL,REAL8IIRFilter *AYLAL)
{
  int l,n;
  

  LALDCreateVector(status->statusPtr,&CinvLAL->directCoef,CinvDirectOrder);
  LALDCreateVector(status->statusPtr,&CinvLAL->recursCoef,CinvRecursOrder);
  LALDCreateVector(status->statusPtr,&CinvLAL->history,CinvDirectOrder-1);

  for(l=0;l<CinvDirectOrder;l++) CinvLAL->directCoef->data[l]=CinvDirectCoefs[l];
  for(l=0;l<CinvDirectOrder;l++) CinvLAL->recursCoef->data[l]=-CinvRecursCoefs[l];
  for(l=0;l<CinvDirectOrder-1;l++) CinvLAL->history->data[l]=0.0;

  for(n=0;n<NGfilt;n++){

    LALDCreateVector(status->statusPtr,&GLAL[n].directCoef,G_Dord);
    LALDCreateVector(status->statusPtr,&GLAL[n].recursCoef,G_Dord);
    LALDCreateVector(status->statusPtr,&GLAL[n].history,G_Dord-1);
   
    /* Fill coefficient vectors with coefficients from filters.h */
    for(l=0;l<G_Dord;l++) GLAL[n].directCoef->data[l]=G_D[n][l];
    for(l=0;l<G_Dord;l++) GLAL[n].recursCoef->data[l]=-G_R[n][l];
    for(l=0;l<G_Dord-1;l++) GLAL[n].history->data[l]=0.0;

  }


  LALDCreateVector(status->statusPtr,&AALAL->directCoef, A_0_Rord);
  LALDCreateVector(status->statusPtr,&AALAL->recursCoef, A_0_Rord);
  LALDCreateVector(status->statusPtr,&AALAL->history, A_0_Rord-1);

  /* Fill coefficient vectors with coefficients from filters.h */
  for(l=0;l< A_0_Rord;l++) AALAL->directCoef->data[l]=A_0_D[l];
  for(l=0;l< A_0_Rord;l++) AALAL->recursCoef->data[l]=-A_0_R[l];
  for(l=0;l< A_0_Rord-1;l++) AALAL->history->data[l]=0.0;


  for(n=0;n<NAXfilt;n++){
   
    LALDCreateVector(status->statusPtr,&AXLAL[n].directCoef, A_digital_Rord);
    LALDCreateVector(status->statusPtr,&AXLAL[n].recursCoef, A_digital_Rord);
    LALDCreateVector(status->statusPtr,&AXLAL[n].history, A_digital_Rord-1);

    for(l=0;l< A_digital_Rord;l++) AXLAL[n].directCoef->data[l]=AX_D[n][l];
    for(l=0;l< A_digital_Rord;l++) AXLAL[n].recursCoef->data[l]=-AX_R[n][l];
    for(l=0;l< A_digital_Rord-1;l++) AXLAL[n].history->data[l]=0.0;
    
  }

  for(n=0;n<NAYfilt;n++){
    LALDCreateVector(status->statusPtr,&AYLAL[n].directCoef, A_digital_Rord);
    LALDCreateVector(status->statusPtr,&AYLAL[n].recursCoef, A_digital_Rord);
    LALDCreateVector(status->statusPtr,&AYLAL[n].history, A_digital_Rord-1);

    for(l=0;l< A_digital_Rord;l++) AYLAL[n].directCoef->data[l]=AY_D[n][l];
    for(l=0;l< A_digital_Rord;l++) AYLAL[n].recursCoef->data[l]=-AY_R[n][l];
    for(l=0;l< A_digital_Rord-1;l++) AYLAL[n].history->data[l]=0.0;
  }

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

/*       fprintf(stdout,"%e %e %e %e\n",output->alpha.data->data[m].re,output->alpha.data->data[m].im, */
/* 	                             output->beta.data->data[m].re,output->beta.data->data[m].im); */



      if(m == MAXALPHAS)
	{
	  fprintf(stderr,"Too many values of the factors, maximum allowed is %d\n",MAXALPHAS);
	  return 1;
	}
    }

  /* Clean up */
  LALDestroyVector(status->statusPtr,&darm.data);
  LALDestroyVector(status->statusPtr,&exc.data);
  LALDestroyVector(status->statusPtr,&asq.data);

  LALDestroyVector(status->statusPtr,&asqwin);
  LALDestroyVector(status->statusPtr,&darmwin);
  LALDestroyVector(status->statusPtr,&excwin);

  return 0;
}

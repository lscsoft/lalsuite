/******** <lalVerbatim file="FreqSeriesToTFPlaneCV"> ********
Author: Flanagan, E
$Id$
********* </lalVerbatim> ********/


#include <lal/LALRCSID.h>


NRCSID (FREQSERIESTOTFPLANEC, "$Id$");



#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/SeqFactories.h>
#include <lal/RealFFT.h>
#include <lal/ComplexFFT.h>
#include <lal/TFTransform.h>

/******** <lalVerbatim file="FreqSeriesToTFPlaneCP"> ********/
void
LALFreqSeriesToTFPlane (
	       LALStatus                               *status,
	       COMPLEX8TimeFrequencyPlane           *tfp,
	       COMPLEX8FrequencySeries              *freqSeries,
	       HorizontalTFTransformIn              *input
		    )
/******** </lalVerbatim> ********/
{
  COMPLEX8Vector     *tmp  = NULL;
  COMPLEX8Vector     *tmp1 = NULL;
  REAL4              fac;
  INT4               i;
  INT4               j;
  INT4               nt;
  INT4               nf;
  INT4               ntotal;

  INT4               flow1;
  INT4               fseglength;

  INITSTATUS (status, "LALFreqSeriesToTFPlane", FREQSERIESTOTFPLANEC);
  ATTATCHSTATUSPTR (status);

  

  /* make sure that arguments are not NULL */
  ASSERT (freqSeries, status, TFTRANSFORMH_ENULLP, TFTRANSFORMH_MSGENULLP);
  ASSERT (freqSeries->data, status, TFTRANSFORMH_ENULLP, TFTRANSFORMH_MSGENULLP);
  ASSERT (freqSeries->data->data, status, TFTRANSFORMH_ENULLP,
          TFTRANSFORMH_MSGENULLP);

  ASSERT (input, status, TFTRANSFORMH_ENULLP, TFTRANSFORMH_MSGENULLP);
  ASSERT (input->dftParams, status, TFTRANSFORMH_ENULLP, TFTRANSFORMH_MSGENULLP);
  ASSERT (input->dftParams->plan, status, TFTRANSFORMH_ENULLP, 
          TFTRANSFORMH_MSGENULLP);
  ASSERT (input->dftParams->window, status, TFTRANSFORMH_ENULLP, 
          TFTRANSFORMH_MSGENULLP);
  ASSERT (input->dftParams->window->data, status, TFTRANSFORMH_ENULLP, 
          TFTRANSFORMH_MSGENULLP);



  /* make sure that output structure is not NULL */
  ASSERT (tfp, status, TFTRANSFORMH_ENULLP, TFTRANSFORMH_MSGENULLP);
  ASSERT (tfp->params, status, TFTRANSFORMH_ENULLP, TFTRANSFORMH_MSGENULLP);
  ASSERT (tfp->data, status, TFTRANSFORMH_ENULLP, TFTRANSFORMH_MSGENULLP);




  /*
   *
   *
   *  make sure input parameters are reasonable, compatible with
   *  each other, etc.
   *
   *
   */



  nt = tfp->params->timeBins;   /* Number of time bins */
  ASSERT (nt > 0, status, TFTRANSFORMH_EPOSARG, TFTRANSFORMH_MSGEPOSARG);

  nf = tfp->params->freqBins;   /* Number of frequency bins */
  ASSERT (nf > 0, status, TFTRANSFORMH_EPOSARG, TFTRANSFORMH_MSGEPOSARG);

  /* 
   * Next compute fseglength = size of segments in freq domain
   * Round off to nearest integer rather than rounding down.
   */

  ASSERT( freqSeries->deltaF>0.0, status, TFTRANSFORMH_EPOSARG, 
         TFTRANSFORMH_MSGEPOSARG);  
  ASSERT( tfp->params->deltaT>0.0, status, TFTRANSFORMH_EPOSARG, 
         TFTRANSFORMH_MSGEPOSARG);  
  fseglength = (INT4)( 0.5+ 1.0/(tfp->params->deltaT * freqSeries->deltaF));  

  /*
   * Length of segments in freq domain must exceed number of time bins
   * in final TF plane
   */
  ASSERT( fseglength >= nt, status, TFTRANSFORMH_EINCOMP, 
          TFTRANSFORMH_MSGEINCOMP);

  ASSERT( fseglength == (INT4)input->dftParams->window->length, status, 
          TFTRANSFORMH_EINCOMP, TFTRANSFORMH_MSGEINCOMP);


  /* Input lowest frequency must be non-negative */
  ASSERT(freqSeries->f0 >= 0.0, status, TFTRANSFORMH_EPOSARG,
         TFTRANSFORMH_MSGEPOSARG);

  /* Lowest freq of time freq plane >= lowest freq of freq series */
  ASSERT(tfp->params->flow >= freqSeries->f0, status, 
         TFTRANSFORMH_EINCOMP, TFTRANSFORMH_MSGEINCOMP);

  /* low frequency cutoff */
  flow1 = (INT4)( (tfp->params->flow - freqSeries->f0)/ freqSeries->deltaF);

  /* compute total number of data points to be used to construct TF plane */
  ntotal = nf * fseglength;

  /* make sure have enough data points in freq series */
  ASSERT(ntotal + flow1<= (INT4)freqSeries->data->length, status, 
         TFTRANSFORMH_EINCOMP, TFTRANSFORMH_MSGEINCOMP);

  /* sumofsquares parameter must be positive */
  ASSERT( input->dftParams->sumofsquares>0.0, status, TFTRANSFORMH_EPOSARG, 
          TFTRANSFORMH_MSGEPOSARG);


  /* 
   * copy some of the information from the input frequency series
   * to the TFplane structure
   *
   */
  tfp->epoch = freqSeries->epoch;

  /* the parameter input->startT is not used by this function */



  /* This TF plane is a horizontal type */
  tfp->planeType = horizontalPlane;

  /*
   *
   *
   *  Now start the computation of the TF transform.
   *
   *
   */

  /* create temporary vectors */
  LALCCreateVector (status->statusPtr, &tmp, fseglength);
  CHECKSTATUSPTR (status);
  LALCCreateVector (status->statusPtr, &tmp1, fseglength);
  CHECKSTATUSPTR (status);

  fac = 1/sqrt(input->dftParams->sumofsquares);

  /* loop over successive data segments in the freq domain */
  for(i=0; i< nf; i++)
    {
      INT4 offset = flow1 + i * fseglength;

      /* Get the segment of freq domain data and window it */
      for(j=0; j< fseglength; j++)
	{
	  tmp->data[j].re = fac * freqSeries->data->data[offset+j].re * input->dftParams->window->data[j]; 
	  tmp->data[j].im = fac * freqSeries->data->data[offset+j].im * input->dftParams->window->data[j]; 
	}

      LALCOMPLEX8VectorFFT (status->statusPtr, tmp1, tmp, input->dftParams->plan);
      CHECKSTATUSPTR (status);
            
      /* Copy the result into appropriate spot in output structure */
      for(j=0; j< nt; j++)
	{
	  tfp->data[j*nf+i] = tmp1->data[j];
	}
       
    }

  LALCDestroyVector (status->statusPtr, &tmp);
  CHECKSTATUSPTR (status);

  LALCDestroyVector (status->statusPtr, &tmp1);
  CHECKSTATUSPTR (status);

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}

/******** <lalVerbatim file="ModFreqSeriesToTFPlaneCP"> ********/
void
LALModFreqSeriesToTFPlane (
    LALStatus                            *status,
    COMPLEX8TimeFrequencyPlane           *tfp,
    COMPLEX8FrequencySeries              *freqSeries,
    HorizontalTFTransformIn              *input
    )
/******** </lalVerbatim> ********/
{
  REAL4Vector        *filter = NULL;
  REAL4Vector        *snr    = NULL;
  COMPLEX8Vector     *tmp    = NULL;
  COMPLEX8Vector     *fcorr  = NULL;
  REAL4              fac;
  INT4               i;
  INT4               j;
  INT4               numpts=0;
  INT4               nt;
  INT4               nf;
  INT4               ntotal;

  INT4               flow1;
  INT4               fseglength;

  REAL4          delF=0;
  REAL4          delT=0;
  REAL4          dt=0;
  REAL4          twopiOverNumpts=0;
  REAL4          norm = 0.0;
  INT4           filterlen=0;

  RealFFTPlan *pfwd = NULL;   /* FFTW uses a plan to assess best FFT method */
  RealFFTPlan *prev = NULL;   /* This one is for the reverse FFT */

  INITSTATUS (status, "LALModFreqSeriesToTFPlane", FREQSERIESTOTFPLANEC);
  ATTATCHSTATUSPTR (status);

  
  /* make sure that arguments are not NULL */
  ASSERT (freqSeries, status, TFTRANSFORMH_ENULLP, TFTRANSFORMH_MSGENULLP);
  ASSERT (freqSeries->data, status, TFTRANSFORMH_ENULLP, TFTRANSFORMH_MSGENULLP);
  ASSERT (freqSeries->data->data, status, TFTRANSFORMH_ENULLP,
      TFTRANSFORMH_MSGENULLP);
  ASSERT (input, status, TFTRANSFORMH_ENULLP, TFTRANSFORMH_MSGENULLP);
  ASSERT (input->dftParams, status, TFTRANSFORMH_ENULLP, TFTRANSFORMH_MSGENULLP);
  ASSERT (input->dftParams->plan, status, TFTRANSFORMH_ENULLP, 
      TFTRANSFORMH_MSGENULLP);
  ASSERT (input->dftParams->window, status, TFTRANSFORMH_ENULLP, 
      TFTRANSFORMH_MSGENULLP);
  ASSERT (input->dftParams->window->data, status, TFTRANSFORMH_ENULLP, 
      TFTRANSFORMH_MSGENULLP);


  /* make sure that output structure is not NULL */
  ASSERT (tfp, status, TFTRANSFORMH_ENULLP, TFTRANSFORMH_MSGENULLP);
  ASSERT (tfp->params, status, TFTRANSFORMH_ENULLP, TFTRANSFORMH_MSGENULLP);
  ASSERT (tfp->data, status, TFTRANSFORMH_ENULLP, TFTRANSFORMH_MSGENULLP);

  
  /*
   *  make sure input parameters are reasonable, compatible with
   *  each other, etc.
   */
  nt = tfp->params->timeBins;   /* Number of time bins */
  ASSERT (nt > 0, status, TFTRANSFORMH_EPOSARG, TFTRANSFORMH_MSGEPOSARG);

  nf = tfp->params->freqBins;   /* Number of frequency bins */
  ASSERT (nf > 0, status, TFTRANSFORMH_EPOSARG, TFTRANSFORMH_MSGEPOSARG);
  

  /* 
   * Next compute fseglength = size of segments in freq domain
   * Round off to nearest integer rather than rounding down.
   */
  ASSERT( freqSeries->deltaF>0.0, status, TFTRANSFORMH_EPOSARG, 
      TFTRANSFORMH_MSGEPOSARG);  
  ASSERT( tfp->params->deltaT>0.0, status, TFTRANSFORMH_EPOSARG, 
      TFTRANSFORMH_MSGEPOSARG);  

  delF = 1.0 / tfp->params->deltaT;
  delT = 1.0 / delF;
  fseglength = (INT4)( 0.5+ 1.0/(tfp->params->deltaT * freqSeries->deltaF));  

  
  /*
   * Length of segments in freq domain must exceed number of time bins
   * in final TF plane
   */
  ASSERT( fseglength >= nt, status, TFTRANSFORMH_EINCOMP, 
      TFTRANSFORMH_MSGEINCOMP);

  
  /* Input lowest frequency must be non-negative */
  ASSERT(freqSeries->f0 >= 0.0, status, TFTRANSFORMH_EPOSARG,
      TFTRANSFORMH_MSGEPOSARG);

  
  /* Lowest freq of time freq plane >= lowest freq of freq series */
  ASSERT(tfp->params->flow >= freqSeries->f0, status, 
      TFTRANSFORMH_EINCOMP, TFTRANSFORMH_MSGEINCOMP);

  
  /* 
   * low frequency cutoff:  in terms of number of bins of frequency
   * series 
   */
  flow1 = (INT4)( (tfp->params->flow - freqSeries->f0)/ freqSeries->deltaF);

  
  /* compute total number of data points to be used to construct TF plane */
  ntotal = nf * fseglength;

  
  /* make sure have enough data points in freq series */
  ASSERT(ntotal + flow1<= (INT4)freqSeries->data->length, status, 
      TFTRANSFORMH_EINCOMP, TFTRANSFORMH_MSGEINCOMP);

  /* sumofsquares parameter must be positive */
  ASSERT( input->dftParams->sumofsquares>0.0, status, TFTRANSFORMH_EPOSARG, 
      TFTRANSFORMH_MSGEPOSARG);


  /* 
   * copy some of the information from the input frequency series
   * to the TFplane structure
   */
  tfp->epoch = freqSeries->epoch;


  /* This TF plane is a horizontal type */
  tfp->planeType = horizontalPlane;

  /*
   *  Now start the computation of the TF transform.
   */

  /* create temporary vectors */
  numpts = 2*(freqSeries->data->length-1);
  twopiOverNumpts = 2.0 * LAL_PI / (float)numpts;
  LALSCreateVector (status->statusPtr, &filter, numpts);
  CHECKSTATUSPTR (status);
  LALCCreateVector (status->statusPtr, &tmp, freqSeries->data->length);
  CHECKSTATUSPTR (status);
  LALCCreateVector (status->statusPtr, &fcorr, freqSeries->data->length);
  CHECKSTATUSPTR (status);
  LALSCreateVector (status->statusPtr, &snr, numpts);
  CHECKSTATUSPTR (status);

  
  /* sampling rate of time series which gave freqSeries */
  dt = 1.0 / (((REAL4) numpts) * freqSeries->deltaF) ;
  filterlen =  (INT4) (delT/dt);

  /* Create an FFTW plan for forward REAL FFT */
  LALCreateForwardRealFFTPlan( status->statusPtr, &pfwd, numpts, 0);
  CHECKSTATUSPTR (status);

  /* Create an FFTW plan for forward REAL FFT */
  LALCreateReverseRealFFTPlan( status->statusPtr, &prev, numpts, 0); 
  CHECKSTATUSPTR (status);

  /* loop over different basis vectors in the frequency domain */
  for(i=0; i< nf; i++)
  {
    INT4 offset = flow1 + i * fseglength;

    /* PRB - generate the time domain filter function.  This filter should
     * depends on the frequency flow1 and the fseglength.  It should
     * be non-zero for some amount that depends on fseglength and
     * the total bandwidth of the input frequency series.  
     */
    memset(filter->data, 0, numpts * sizeof(REAL4));
    filter->data[0] = twopiOverNumpts * fseglength / (LAL_PI * dt);
    for ( j=1 ; j<filterlen ; j++){
      filter->data[j]=(sin(twopiOverNumpts * j * (offset + fseglength) ) 
          - sin(twopiOverNumpts * j * offset)) / 
        (LAL_PI * 1.0 * ((float)j) * dt);
      filter->data[numpts-j]=(sin(twopiOverNumpts * j * (offset + fseglength) ) 
          - sin(twopiOverNumpts * j * offset )) / (LAL_PI * 1.0 * ((float)j) * dt);
    }

    /* PRB - Fourier transform the filter inot the frequency domain */
    LALForwardRealFFT( status->statusPtr, tmp, filter, pfwd);
    CHECKSTATUSPTR (status);

    /* PRB - Multiply the filter by the data.  Don't forget complex
     * conjugate and any other relevant information */
    for( j=0 ; j<freqSeries->data->length ; j++)
    {
      fcorr->data[j].re = tmp->data[j].re * freqSeries->data->data[j].re
        + tmp->data[j].im * freqSeries->data->data[j].im;
      fcorr->data[j].im = tmp->data[j].re * freqSeries->data->data[j].im
        - tmp->data[j].im * freqSeries->data->data[j].re;
    }

    /* PRB - The normalization constant */
    norm = 0.0;
    for( j=0 ; j<freqSeries->data->length ; j++)
    {
      REAL4 re = tmp->data[j].re;
      REAL4 im = tmp->data[j].im;

      norm += (re*re + im*im);
    }
    norm = sqrt(4.0 * norm);

    /* PRB - Inverse transform the product so that we get a time series at
     * the full sample rate.   Make sure to check that the sample
     * rate agrees with what you had before. */
    LALReverseRealFFT (status->statusPtr, snr, fcorr, prev);
    CHECKSTATUSPTR (status);

    /* PRB - copy the data back into the time series.  In this
     * process,  one only takes the samples corresponding to the
     * uncorupted times.   This means that the data should be longer
     * than the 1/dfmin.   This can mean a change compared to what
     * has been done before.   We'll have to check that carefully.
     * Notice that the originally complex plane is now real.  I still
     * don't understand why the difference arises with the Eanns's
     * original implementation.   We'll have to look at that.  */
    /* Copy the result into appropriate spot in output structure */
    for(j=0; j< nt; j++)
    {
      tfp->data[j*nf+i].re = snr->data[(j+1)*filterlen] / norm;
      tfp->data[j*nf+i].im = 0.0;
    }

    if( nf == 64 ){

      for(j=0; j< numpts; j++)
      {
        snr->data[j] /= norm;
      }
      norm =0;
      for(j=0; j< numpts; j++)
      {
        norm += snr->data[j] * snr->data[j];
      }
      fprintf(stdout,"%e\n",norm / numpts);
      fflush(stdout);
    }
  }

  /* Get rid of the FFT plans */
  LALDestroyRealFFTPlan( status->statusPtr, &pfwd ); 
  LALDestroyRealFFTPlan( status->statusPtr, &prev );

  LALSDestroyVector (status->statusPtr, &filter);
  CHECKSTATUSPTR (status);

  LALCDestroyVector (status->statusPtr, &tmp);
  CHECKSTATUSPTR (status);

  LALCDestroyVector (status->statusPtr, &fcorr);
  CHECKSTATUSPTR (status);

  LALSDestroyVector (status->statusPtr, &snr);
  CHECKSTATUSPTR (status);

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}





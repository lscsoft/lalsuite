/******** <lalVerbatim file="FreqSeriesToTFPlaneCV"> ********
Author: Flanagan, E
$Id$
********* </lalVerbatim> ********/


#include <lal/LALRCSID.h>


NRCSID (FREQSERIESTOTFPLANEC, "$Id$");



#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <lal/ComplexFFT.h>
#include <lal/LALConstants.h>
#include <lal/LALErrno.h>
#include <lal/LALStdlib.h>
#include <lal/RealFFT.h>
#include <lal/SeqFactories.h>
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
  ASSERT(freqSeries, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(freqSeries->data, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(freqSeries->data->data, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(input, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(input->dftParams, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(input->dftParams->plan, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(input->dftParams->window, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(input->dftParams->window->data, status, LAL_NULL_ERR, LAL_NULL_MSG);

  /* make sure that output structure is not NULL */
  ASSERT(tfp, status, LAL_NNULL_ERR, LAL_NNULL_MSG);
  ASSERT(tfp->params, status, LAL_NNULL_ERR, LAL_NNULL_MSG);
  ASSERT(tfp->data, status, LAL_NNULL_ERR, LAL_NNULL_MSG);


  /*
   *  make sure input parameters are reasonable, compatible with
   *  each other, etc.
   */

  nt = tfp->params->timeBins;   /* Number of time bins */
  ASSERT(nt > 0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);

  nf = tfp->params->freqBins;   /* Number of frequency bins */
  ASSERT(nf > 0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);


  /* 
   * Next compute fseglength = size of segments in freq domain
   * Round off to nearest integer rather than rounding down.
   */

  ASSERT(freqSeries->deltaF>0.0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);  
  ASSERT(tfp->params->deltaT>0.0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);  
  fseglength = (INT4)( 0.5+ 1.0/(tfp->params->deltaT * freqSeries->deltaF));  

  /*
   * Length of segments in freq domain must exceed number of time bins
   * in final TF plane
   */
  ASSERT(fseglength >= nt, status, LAL_RANGE_ERR, LAL_RANGE_MSG);
  ASSERT(fseglength == (INT4)input->dftParams->window->length, status, LAL_BADPARM_ERR, LAL_BADPARM_MSG);


  /* Input lowest frequency must be non-negative */
  ASSERT(freqSeries->f0 >= 0.0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);

  /* Lowest freq of time freq plane >= lowest freq of freq series */
  ASSERT(tfp->params->flow >= freqSeries->f0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);

  /* low frequency cutoff */
  flow1 = (INT4)( (tfp->params->flow - freqSeries->f0)/ freqSeries->deltaF);

  /* compute total number of data points to be used to construct TF plane */
  ntotal = nf * fseglength;

  /* make sure have enough data points in freq series */
  ASSERT(ntotal + flow1<= (INT4)freqSeries->data->length, status, LAL_RANGE_ERR, LAL_RANGE_MSG);

  /* sumofsquares parameter must be positive */
  ASSERT(input->dftParams->sumofsquares>0.0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);


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
  INT4               i;
  INT4               j;
  INT4               numpts=0;
  INT4               nt;
  INT4               nf;
  INT4               ntotal;
  INT4               offset;
  UINT4              windowShift;

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
  ASSERT(freqSeries, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(freqSeries->data, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(freqSeries->data->data, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(input, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(input->dftParams, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(input->dftParams->plan, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(input->dftParams->window, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(input->dftParams->window->data, status, LAL_NULL_ERR, LAL_NULL_MSG);


  /* make sure that output structure is not NULL */
  ASSERT(tfp, status, LAL_NNULL_ERR, LAL_NNULL_MSG);
  ASSERT(tfp->params, status, LAL_NNULL_ERR, LAL_NNULL_MSG);
  ASSERT(tfp->data, status, LAL_NNULL_ERR, LAL_NNULL_MSG);

  
  /* check number of time bins is positive */
  nt = tfp->params->timeBins; 
  ASSERT(nt > 0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);

  /* check number of freq bins is positive */
  nf = tfp->params->freqBins;   
  ASSERT(nf > 0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);
  
  windowShift = input->windowShift;
  /* check a couple more input parameters */
  ASSERT(freqSeries->deltaF>0.0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);  
  ASSERT(tfp->params->deltaT>0.0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);  

  /* 
   * delF is the frequency resoltion of the time-frequency plane.  It
   * is directly related to the time resolution so that 
   *    delT * delF = 1
   * Note that the frequency resolution of the frequency series does
   * not have to be the same as delF.
   */
  delT = tfp->params->deltaT;
  delF = 1.0 / delT;

  /* number of bins of the frequency series corresponding to delF */
  fseglength = (INT4)( 0.5 + delF/freqSeries->deltaF );  

  
  /*
   * Length of segments in freq domain must exceed number of time bins
   * in final TF plane
   */
  ASSERT(fseglength >= nt, status, LAL_RANGE_ERR, LAL_RANGE_MSG);

  
  /* Input lowest frequency must be non-negative */
  ASSERT(freqSeries->f0 >= 0.0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);

  
  /* Lowest freq of time freq plane >= lowest freq of freq series */
  ASSERT(tfp->params->flow >= freqSeries->f0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);

  
  /* low frequency cutoff:  in terms of number of bins of frequency series */
  flow1 = (INT4)( (tfp->params->flow - freqSeries->f0)/ freqSeries->deltaF);

  
  /* compute total number of data points to be used to construct TF plane */
  ntotal = nf * fseglength;

  
  /* make sure have enough data points in freq series */
  ASSERT(ntotal + flow1<= (INT4)freqSeries->data->length, status, LAL_RANGE_ERR, LAL_RANGE_MSG);


  /* set the epoch of the TF plane */
  tfp->epoch = freqSeries->epoch;


  /* This TF plane is a horizontal type */
  tfp->planeType = horizontalPlane;


  /* create temporary vectors for filter and correlations */
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

  /* number of points from peak of filter to first zero */
  filterlen =  (INT4) (delT/dt);

  /* Create an FFTW plan for forward REAL FFT */
  LALCreateForwardRealFFTPlan( status->statusPtr, &pfwd, numpts, 0);
  CHECKSTATUSPTR (status);

  /* Create an FFTW plan for forward REAL FFT */
  LALCreateReverseRealFFTPlan( status->statusPtr, &prev, numpts, 0); 
  CHECKSTATUSPTR (status);

    /* PRB - generate the time domain filter function.  This filter should
     * depends on the frequency flow1 and the fseglength.  It should
     * be non-zero for some amount that depends on fseglength and
     * the total bandwidth of the input frequency series.  
     */
    memset(filter->data, 0, numpts * sizeof(REAL4));
    filter->data[0] = twopiOverNumpts * fseglength / (LAL_PI * dt);
    for ( j=1 ; j<filterlen ; j++){
      REAL4 tmpValue = ( sin( twopiOverNumpts * j * (flow1 + fseglength) ) 
                        - sin( twopiOverNumpts * j * flow1)) / 
                        (LAL_PI * 1.0 * ((float)j) * dt);

      filter->data[j] = tmpValue;
      filter->data[numpts-j] = tmpValue;
    }

    /* PRB - Fourier transform the filter into the frequency domain */
    LALForwardRealFFT( status->statusPtr, tmp, filter, pfwd);
    CHECKSTATUSPTR (status);

    /* PRB - The normalization constant */
    norm = 0.0;
    for(j = 0; (unsigned) j < freqSeries->data->length; j++)
    {
      REAL4 re = tmp->data[j].re;
      REAL4 im = tmp->data[j].im;

      norm += (re*re + im*im);
    }
    norm = sqrt(4.0 * norm);

  /* loop over different basis vectors in the frequency domain */
  for(i=0; i< nf; i++)
  {
    offset = flow1 + i * fseglength;

    /* PRB - set all values below i * df to zero */
    for( j=0 ; j < i*fseglength ; j++)
    {
      REAL4 reFilter = 0.0;
      REAL4 imFilter = 0.0;
      REAL4 reData = freqSeries->data->data[j].re;
      REAL4 imData = freqSeries->data->data[j].im;
      
      fcorr->data[j].re = reFilter * reData + imFilter * imData;
      fcorr->data[j].im = reFilter * imData - imFilter * reData;
    }

    /* PRB - Multiply the filter by the data.  Don't forget complex
     * conjugate and any other relevant information */
    for(j = i * fseglength; (unsigned) j < freqSeries->data->length; j++)
    {
      REAL4 reFilter = tmp->data[j-i*fseglength].re;
      REAL4 imFilter = tmp->data[j-i*fseglength].im;
      REAL4 reData = freqSeries->data->data[j].re;
      REAL4 imData = freqSeries->data->data[j].im;
      
      fcorr->data[j].re = reFilter * reData + imFilter * imData;
      fcorr->data[j].im = reFilter * imData - imFilter * reData;
    }

    /* clean up Nyquist */
      fcorr->data[freqSeries->data->length-1].re=0.0;
      fcorr->data[freqSeries->data->length-1].im=0.0;

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
      tfp->data[j*nf+i].re = snr->data[windowShift+(j*filterlen)] / norm;
      tfp->data[j*nf+i].im = 0.0;
    }
  }

  /* Get rid of all temporary memory */
  LALDestroyRealFFTPlan( status->statusPtr, &pfwd ); 
  CHECKSTATUSPTR (status);

  LALDestroyRealFFTPlan( status->statusPtr, &prev );
  CHECKSTATUSPTR (status);

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


void
LALModModFreqSeriesToTFPlane (
    LALStatus                            *status,
    COMPLEX8TimeFrequencyPlane           *tfp,
    COMPLEX8FrequencySeries              *freqSeries,
    HorizontalTFTransformIn              *input,
    REAL4                                *normalisation,
    REAL4FrequencySeries                 *psd
    )
/******** </lalVerbatim> ********/
{
  REAL4Vector        *filter = NULL;
  REAL4Vector        *snr    = NULL;
  COMPLEX8Vector     *tmp    = NULL;
  COMPLEX8Vector     *fcorr  = NULL;
  INT4               i;
  INT4               j;
  INT4               numpts=0;
  INT4               nt;
  INT4               nf;
  INT4               ntotal;
  INT4               offset;
  UINT4              windowShift;

  INT4               flow1;
  INT4               fseglength;

  REAL4          delF=0;
  REAL4          delT=0;
  REAL4          dt=0;
  REAL4          twopiOverNumpts=0;
  REAL4          norm = 0.0;
  INT4           filterlen=0;

  FILE *fp;
  REAL4 zsum = 0;

  RealFFTPlan *pfwd = NULL;   /* FFTW uses a plan to assess best FFT method */
  RealFFTPlan *prev = NULL;   /* This one is for the reverse FFT */

  INITSTATUS (status, "LALModFreqSeriesToTFPlane", FREQSERIESTOTFPLANEC);
  ATTATCHSTATUSPTR (status);

  
  /* make sure that arguments are not NULL */
  ASSERT(freqSeries, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(freqSeries->data, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(freqSeries->data->data, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(input, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(input->dftParams, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(input->dftParams->plan, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(input->dftParams->window, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(input->dftParams->window->data, status, LAL_NULL_ERR, LAL_NULL_MSG);


  /* make sure that output structure is not NULL */
  ASSERT(tfp, status, LAL_NNULL_ERR, LAL_NNULL_MSG);
  ASSERT(tfp->params, status, LAL_NNULL_ERR, LAL_NNULL_MSG);
  ASSERT(tfp->data, status, LAL_NNULL_ERR, LAL_NNULL_MSG);

  
  /* check number of time bins is positive */
  nt = tfp->params->timeBins; 
  ASSERT(nt > 0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);

  /* check number of freq bins is positive */
  nf = tfp->params->freqBins;   
  ASSERT(nf > 0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);
  
  windowShift = input->windowShift;
  /* check a couple more input parameters */
  ASSERT(freqSeries->deltaF>0.0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);  
  ASSERT(tfp->params->deltaT>0.0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);  

  /* 
   * delF is the frequency resoltion of the time-frequency plane.  It
   * is directly related to the time resolution so that 
   *    delT * delF = 1
   * Note that the frequency resolution of the frequency series does
   * not have to be the same as delF.
   */
  delT = tfp->params->deltaT;
  delF = tfp->params->deltaF;

  /*  printf("nt = %d  nf = %d delF = %f delT = %f \n",nt,nf,delF, delT);*/
  /* number of bins of the frequency series corresponding to delF */
  fseglength = (INT4)( 0.5 + delF/freqSeries->deltaF );  
  /* printf("fseglength= %d \n",fseglength);*/
  
  /*
   * Length of segments in freq domain must exceed number of time bins
   * in final TF plane
   */
  /* ASSERT(fseglength >= nt, status, LAL_RANGE_ERR, LAL_RANGE_MSG);*/

  
  /* Input lowest frequency must be non-negative */
  ASSERT(freqSeries->f0 >= 0.0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);

  
  /* Lowest freq of time freq plane >= lowest freq of freq series */
  ASSERT(tfp->params->flow >= freqSeries->f0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);

  
  /* low frequency cutoff:  in terms of number of bins of frequency series */
  flow1 = (INT4)( (tfp->params->flow - freqSeries->f0)/ freqSeries->deltaF);

  
  /* compute total number of data points to be used to construct TF plane */
  /* ntotal = nf * fseglength;*/

  
  /* make sure have enough data points in freq series */
  /*ASSERT(ntotal + flow1<= (INT4)freqSeries->data->length, status, LAL_RANGE_ERR, LAL_RANGE_MSG);*/


  /* set the epoch of the TF plane */
  tfp->epoch = freqSeries->epoch;


  /* This TF plane is a horizontal type */
  tfp->planeType = horizontalPlane;


  /* create temporary vectors for filter and correlations */
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

  /* number of points from peak of filter to first zero */
  /*  filterlen =  (INT4) (delT/dt);*/

  filterlen = (INT4)(numpts/fseglength);

  /*  filterlen = 8192;*/

  /* Create an FFTW plan for forward REAL FFT */
  LALCreateForwardRealFFTPlan( status->statusPtr, &pfwd, numpts, 0);
  CHECKSTATUSPTR (status);

  /* Create an FFTW plan for forward REAL FFT */
  LALCreateReverseRealFFTPlan( status->statusPtr, &prev, numpts, 0); 
  CHECKSTATUSPTR (status);

  /* 
   * PRB - test code to add a delta function to the segment and
   * confirm the output of the code
   */
	/*
	 memset(filter->data, 0, numpts * sizeof(REAL4));
	 filter->data[numpts/2] = 1.0;
	 LALForwardRealFFT( status->statusPtr, freqSeries->data, filter, pfwd);
	 CHECKSTATUSPTR (status);
	 */

    /* PRB - generate the time domain filter function.  This filter should
     * depends on the frequency flow1 and the fseglength.  It should
     * be non-zero for some amount that depends on fseglength and
     * the total bandwidth of the input frequency series.  
     */
    memset(filter->data, 0, numpts * sizeof(REAL4));
    filter->data[0] = twopiOverNumpts * fseglength / (LAL_PI * dt);
    for ( j=1 ; j<filterlen ; j++){
      REAL4 tmpValue = ( sin( twopiOverNumpts * j * (flow1 + fseglength) ) 
                        - sin( twopiOverNumpts * j * flow1)) / 
                        (LAL_PI * 1.0 * ((float)j) * dt);

      filter->data[j] = tmpValue;
      filter->data[numpts-j] = tmpValue;
    }

  
    /* PRB - Fourier transform the filter into the frequency domain */
    LALForwardRealFFT( status->statusPtr, tmp, filter, pfwd);
    CHECKSTATUSPTR (status);

    /* PRB - The normalization constant */
    /*norm = 0.0;
    for(j = 0; (unsigned) j < freqSeries->data->length; j++)
    {
      REAL4 re = tmp->data[j].re;
      REAL4 im = tmp->data[j].im;

      norm += (re*re + im*im);
    }
    norm = sqrt(4.0 * norm);
    */

  /* loop over different basis vectors in the frequency domain */
  for(i=0; i< nf; i++)
  {
    /* PRB - set all values below i * df to zero */
    for( j=0 ; j < i*fseglength ; j++)
    {
      REAL4 reFilter = 0.0;
      REAL4 imFilter = 0.0;
      /*REAL4 reData = freqSeries->data->data[j].re*sqrt(1/psd->data->data[j]);
	REAL4 imData = freqSeries->data->data[j].im*sqrt(1/psd->data->data[j]);*/

      REAL4 reData = freqSeries->data->data[j].re;
      REAL4 imData = freqSeries->data->data[j].im;

      fcorr->data[j].re = reFilter * reData + imFilter * imData;
      fcorr->data[j].im = reFilter * imData - imFilter * reData;

      /*normalisation[i] += (1/psd->data->data[j])*(reFilter*reFilter + imFilter*imFilter);*/
      normalisation[i] += (reFilter*reFilter + imFilter*imFilter);
    }

    /* PRB - Multiply the filter by the data.  Don't forget complex
     * conjugate and any other relevant information */
    for(j = i * fseglength; (unsigned) j < freqSeries->data->length; j++)
    {
      REAL4 reFilter = tmp->data[j-i*fseglength].re;
      REAL4 imFilter = tmp->data[j-i*fseglength].im;
      /*REAL4 reData = freqSeries->data->data[j].re*sqrt(1/psd->data->data[j]);
	REAL4 imData = freqSeries->data->data[j].im*sqrt(1/psd->data->data[j]);*/

      REAL4 reData = freqSeries->data->data[j].re;
      REAL4 imData = freqSeries->data->data[j].im;

      fcorr->data[j].re = reFilter * reData + imFilter * imData;
      fcorr->data[j].im = reFilter * imData - imFilter * reData;

      /*normalisation[i] += (1/psd->data->data[j])*(reFilter*reFilter + imFilter*imFilter);*/
      normalisation[i] += (reFilter*reFilter + imFilter*imFilter);
    }

    normalisation[i] = sqrt(4.0 * normalisation[i]);

    /* clean up Nyquist */
      fcorr->data[freqSeries->data->length-1].re=0.0;
      fcorr->data[freqSeries->data->length-1].im=0.0;

    /* PRB - Inverse transform the product so that we get a time series at
     * the full sample rate.   Make sure to check that the sample
     * rate agrees with what you had before. */
    LALReverseRealFFT (status->statusPtr, snr, fcorr, prev);
    CHECKSTATUSPTR (status);


	/* 
	 * PRB & SRM - to output the time series when the delta function
	 * is inserted instead of other data 
	 */
	/*
	{
		char fname[248];
		sprintf(fname, "z%02d.dat", i);
		fp = fopen(fname,"a");
		for(j=0; j<numpts ; j++)
		{
			fprintf(fp,"%i %f\n",j,snr->data[j]/norm);
		}
		fclose(fp);
	}
	*/

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
	tfp->data[j*nf+i].re = snr->data[(j*(INT4)(delT/dt))+windowShift];
	tfp->data[j*nf+i].im = 0.0;
      }
  }

  /* printf("%f\n", zsum/sqrt(nf));*/

  /* Get rid of all temporary memory */
  LALDestroyRealFFTPlan( status->statusPtr, &pfwd ); 
  CHECKSTATUSPTR (status);

  LALDestroyRealFFTPlan( status->statusPtr, &prev );
  CHECKSTATUSPTR (status);

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





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
  ASSERT (freqSeries, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);
  ASSERT (freqSeries->data, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);
  ASSERT (freqSeries->data->data, status, TFTRANSFORM_ENULLP,
          TFTRANSFORM_MSGENULLP);

  ASSERT (input, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);
  ASSERT (input->dftParams, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);
  ASSERT (input->dftParams->plan, status, TFTRANSFORM_ENULLP, 
          TFTRANSFORM_MSGENULLP);
  ASSERT (input->dftParams->plan->plan, status, TFTRANSFORM_ENULLP, 
          TFTRANSFORM_MSGENULLP);
  ASSERT (input->dftParams->window, status, TFTRANSFORM_ENULLP, 
          TFTRANSFORM_MSGENULLP);
  ASSERT (input->dftParams->window->data, status, TFTRANSFORM_ENULLP, 
          TFTRANSFORM_MSGENULLP);



  /* make sure that output structure is not NULL */
  ASSERT (tfp, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);
  ASSERT (tfp->params, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);
  ASSERT (tfp->data, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);




  /*
   *
   *
   *  make sure input parameters are reasonable, compatible with
   *  each other, etc.
   *
   *
   */



  nt = tfp->params->timeBins;   /* Number of time bins */
  ASSERT (nt > 0, status, TFTRANSFORM_EPOSARG, TFTRANSFORM_MSGEPOSARG);

  nf = tfp->params->freqBins;   /* Number of frequency bins */
  ASSERT (nf > 0, status, TFTRANSFORM_EPOSARG, TFTRANSFORM_MSGEPOSARG);

  /* 
   * Next compute fseglength = size of segments in freq domain
   * Round off to nearest integer rather than rounding down.
   */

  ASSERT( freqSeries->deltaF>0.0, status, TFTRANSFORM_EPOSARG, 
         TFTRANSFORM_MSGEPOSARG);  
  ASSERT( tfp->params->deltaT>0.0, status, TFTRANSFORM_EPOSARG, 
         TFTRANSFORM_MSGEPOSARG);  
  fseglength = (INT4)( 0.5+ 1.0/(tfp->params->deltaT * freqSeries->deltaF));  

  /*
   * Length of segments in freq domain must exceed number of time bins
   * in final TF plane
   */
  ASSERT( fseglength >= nt, status, TFTRANSFORM_EINCOMP, 
          TFTRANSFORM_MSGEINCOMP);

  ASSERT( fseglength == (INT4)input->dftParams->plan->size, status, 
          TFTRANSFORM_EINCOMP, TFTRANSFORM_MSGEINCOMP);
  ASSERT( fseglength == (INT4)input->dftParams->window->length, status, 
          TFTRANSFORM_EINCOMP, TFTRANSFORM_MSGEINCOMP);

  /* Supplied FFT plan must be in inverse direction */
  ASSERT( input->dftParams->plan->sign==-1, status, 
          TFTRANSFORM_EINCOMP, TFTRANSFORM_MSGEINCOMP);

  /* Input lowest frequency must be non-negative */
  ASSERT(freqSeries->f0 >= 0.0, status, TFTRANSFORM_EPOSARG,
         TFTRANSFORM_MSGEPOSARG);

  /* Lowest freq of time freq plane >= lowest freq of freq series */
  ASSERT(tfp->params->flow >= freqSeries->f0, status, 
         TFTRANSFORM_EINCOMP, TFTRANSFORM_MSGEINCOMP);

  /* low frequency cutoff */
  flow1 = (INT4)( (tfp->params->flow - freqSeries->f0)/ freqSeries->deltaF);

  /* compute total number of data points to be used to construct TF plane */
  ntotal = nf * fseglength;

  /* make sure have enough data points in freq series */
  ASSERT(ntotal + flow1<= (INT4)freqSeries->data->length, status, 
         TFTRANSFORM_EINCOMP, TFTRANSFORM_MSGEINCOMP);

  /* sumofsquares parameter must be positive */
  ASSERT( input->dftParams->sumofsquares>0.0, status, TFTRANSFORM_EPOSARG, 
          TFTRANSFORM_MSGEPOSARG);


  /* 
   * copy some of the information from the input frequency series
   * to the TFplane structure
   *
   */
  /*
   * OMITTED
   *
  tfp->sampleUnits = freqSeries->sampleUnits;
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




/******** <lalVerbatim file="CreateTFPlaneCV"> ********
Author: Flanagan, E
$Id$
********* </lalVerbatim> ********/


#include <lal/LALRCSID.h>


NRCSID (CREATETFPLANEC, "$Id$");



#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/SeqFactories.h>
#include <lal/RealFFT.h>
#include <lal/ComplexFFT.h>
#include <lal/TFTransform.h>

/******** <lalVerbatim file="CreateTFPlaneCP"> ********/
void
LALCreateTFPlane (
	       LALStatus                               *status,
	       COMPLEX8TimeFrequencyPlane           **tfp,
	       TFPlaneParams                        *input
	       )
/******** </lalVerbatim> ********/
{
  INITSTATUS (status, "LALCreateTFPlane", CREATETFPLANEC);

  /* Check input structure: report if NULL */
  ASSERT (input, status, TFTRANSFORMH_ENULLP, TFTRANSFORMH_MSGENULLP);
      
  /* Make sure that input parameters are reasonable */
  ASSERT (input->timeBins > 0, status, 
          TFTRANSFORMH_EPOSARG, TFTRANSFORMH_MSGEPOSARG);
  ASSERT (input->freqBins > 0, status,
          TFTRANSFORMH_EPOSARG, TFTRANSFORMH_MSGEPOSARG);
  ASSERT (input->deltaT > 0.0, status,
          TFTRANSFORMH_EPOSARG, TFTRANSFORMH_MSGEPOSARG);

  /* 
   * Check return structure: tfp should point to a valid pointer
   * which should not yet point to anything.
   *
   */

  ASSERT (tfp != NULL, status, TFTRANSFORMH_ENULLP, TFTRANSFORMH_MSGENULLP);
  ASSERT (*tfp == NULL, status, TFTRANSFORMH_EALLOCP, TFTRANSFORMH_MSGEALLOCP);


  /*  Assign memory for *tfp   */
  *tfp = (COMPLEX8TimeFrequencyPlane *) LALMalloc(sizeof(COMPLEX8TimeFrequencyPlane));
  
  /*  Make sure that the allocation was succesful */
  if ( !(*tfp) ){
     ABORT (status, TFTRANSFORMH_EMALLOC, TFTRANSFORMH_MSGEMALLOC);
  }

  /* assign memory for params field */
  (*tfp)->params = NULL;
  (*tfp)->params = (TFPlaneParams *) LALMalloc(sizeof (TFPlaneParams));

  /*  Make sure that the allocation was succesful */
  if ( !((*tfp)->params) ){
    LALFree ( *tfp );
    ABORT (status, TFTRANSFORMH_EMALLOC, TFTRANSFORMH_MSGEMALLOC);
  }


  /* 
   *  Fill some of the fields with nominal values pending the 
   *  allocation of correct values for these fields.
   */
  /* (*tfp)->name = NULL; */
  /* (*tfp)->sampleUnits=NULL; */
  (*tfp)->epoch.gpsSeconds=0;
  (*tfp)->epoch.gpsNanoSeconds=0;
  (*tfp)->data = NULL;    /* until allocated below */
  

  /* 
   * Allocate storage 
   */

  (*tfp)->data = (COMPLEX8 *) 
    LALMalloc (input->timeBins * input->freqBins * sizeof(COMPLEX8));

  if ( !((*tfp)->data) )
  {
    /* Must free storage pointed to by *tfp */
    LALFree ( (*tfp)->params );
    LALFree ( *tfp );
    ABORT (status, TFTRANSFORMH_EMALLOC, TFTRANSFORMH_MSGEMALLOC);
  }
 

  /* 
   *  Set timeBins, freqBins etc.
   *  by copying the values from the input structure 
   */
  *((*tfp)->params) = *input;

  /* Normal exit */
  RETURN (status);
}



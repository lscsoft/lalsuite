/******** <lalVerbatim file="CreateTFPlaneCV"> ********
Author: Flanagan, E
$Id$
********* </lalVerbatim> ********/


#include <lal/LALRCSID.h>


NRCSID (CREATETFPLANEC, "$Id$");



#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <lal/ComplexFFT.h>
#include <lal/LALErrno.h>
#include <lal/LALStdlib.h>
#include <lal/RealFFT.h>
#include <lal/SeqFactories.h>
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
  ASSERT(input, status, LAL_NULL_ERR, LAL_NULL_MSG);
      
  /* Make sure that input parameters are reasonable */
  ASSERT(input->timeBins > 0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);
  ASSERT(input->freqBins > 0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);
  ASSERT(input->deltaT > 0.0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);

  /* 
   * Check return structure: tfp should point to a valid pointer
   * which should not yet point to anything.
   *
   */

  ASSERT(tfp != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(*tfp == NULL, status, LAL_NNULL_ERR, LAL_NNULL_MSG);


  /*  Assign memory for *tfp   */
  *tfp = LALMalloc(sizeof(**tfp));
  ASSERT(*tfp, status, LAL_NOMEM_ERR, LAL_NOMEM_MSG);

  /* assign memory for params field */
  (*tfp)->params = LALMalloc(sizeof(*(*tfp)->params));

  /*  Make sure that the allocation was succesful */
  if ( !((*tfp)->params) ){
    LALFree ( *tfp );
    ABORT(status, LAL_NOMEM_ERR, LAL_NOMEM_MSG);
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

  (*tfp)->data = LALMalloc(input->timeBins * input->freqBins * sizeof(COMPLEX8));

  if ( !((*tfp)->data) )
  {
    /* Must free storage pointed to by *tfp */
    LALFree ( (*tfp)->params );
    LALFree ( *tfp );
    ABORT (status, LAL_NOMEM_ERR, LAL_NOMEM_MSG);
  }
 

  /* 
   *  Set timeBins, freqBins etc.
   *  by copying the values from the input structure 
   */
  *((*tfp)->params) = *input;

  /* Normal exit */
  RETURN (status);
}



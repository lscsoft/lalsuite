/******** <lalVerbatim file="DestroyTFPlaneCV"> ********
Author: Flanagan, E
$Id$
********* </lalVerbatim> ********/


#include <lal/LALRCSID.h>


NRCSID (DESTROYTFPLANEC, "$Id$");



#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/SeqFactories.h>
#include <lal/RealFFT.h>
#include <lal/ComplexFFT.h>
#include <lal/TFTransform.h>


/******** <lalVerbatim file="DestroyTFPlaneCP"> ********/
void
LALDestroyTFPlane (
	       LALStatus                               *status,
	       COMPLEX8TimeFrequencyPlane           **tfp
	       )
/******** </lalVerbatim> ********/
{
  INITSTATUS (status, "LALDestroyTFPlane", DESTROYTFPLANEC);

  /* make sure that arguments are not null */
  ASSERT (tfp, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);
  ASSERT (*tfp, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);

  /* make sure that data pointed to is non-null */
  ASSERT ((*tfp)->data, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);

  /* make sure that pointer to parameter structure is non-null */
  ASSERT ((*tfp)->params, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP); 

  /* Ok, now let's free allocated storage */
  LALFree ( (*tfp)->params );      /* free parameter structure         */
  LALFree ( (*tfp)->data );        /* free allocated data              */
  LALFree ( *tfp );	           /* free TF plane struct itself      */

  *tfp = NULL;	     	      /* make sure we don't point to freed struct */

  /* Normal exit */
  RETURN (status);
}




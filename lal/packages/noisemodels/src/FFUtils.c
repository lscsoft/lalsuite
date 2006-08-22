/* ****************************************************************
 * Author: Babak, S., B. Krishnan
 * ************************************************************** */
#include "lal/FFUtils.h"

/* --- version information --- */
NRCSID( FFUTILSC, "$Id$");
RCSID(  "$Id$");

#define CVS_ID_STRING_C      "$Id$"
#define CVS_REVISION_C      "$Revision$"

/* Compute the length of the longest template */

/*!!!!!!!! Need to change randIn to bankIn */


void GetMaximumTemplateSize(LALStatus *status, RandomInspiralSignalIn  randIn, 
		      UINT4 *checkLength)
{

   INITSTATUS(status, "BEGetMaximumSize", BANKEFFICIENCYNEWC);
   ATTACHSTATUSPTR(status);
   
   InspiralTemplate dummy;

   dummy = randIn.param;
   dummy.massChoice = m1Andm2;
   dummy.mass1 = signalMinMass;
   dummy.mass2 = signalMinMass;
   dummy.approximant = EOB;
   *checkLength = 0;

   LALInspiralWaveLength(status->statusPtr, checkLength, dummy);

   DETACHSTATUSPTR(status);
   RETURN(status);	

}


/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpSBCVTemplate.c
 *
 * Author: Brown D. A., Spinning BCV modifications by Jones, G
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>

NRCSID (FINDCHIRPSBCVTEMPLATEC, "$Id$");

/*documenation later*/
void
LALFindChirpBCVSpinTemplate (
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *tmplt,
    FindChirpSPTmpltParams     *params
		)

{
/*declaration*/
  INITSTATUS( status, "LALFindChirpBCVSpinTemplate", FINDCHIRPSBCVTEMPLATEC );
  ATTATCHSTATUSPTR( status );

  /*code*/

  DETATCHSTATUSPTR( status );
  RETURN( status );
}


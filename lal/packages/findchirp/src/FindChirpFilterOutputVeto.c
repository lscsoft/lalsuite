/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpFilterOutputVeto.c
 *
 * Author: 
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpFilterOutputVeto.h>

double rint(double x);

NRCSID (FINDCHIRPFILTEROUTPUTVETOC, "$Id$");

void LALFindChirpFilterOutputVeto( 
    LALStatus                          *status,
    SnglInspiralTable                 **eventList, 
    COMPLEX8Vector                     *qVec,
    REAL4                               qNorm,
    FindChirpFilterOutputVetoParams    *params
    )
{
  INITSTATUS( status, "LALFindChirpFilterFilterOutputVeto", 
      FINDCHIRPFILTEROUTPUTVETOC );
  ATTATCHSTATUSPTR( status );


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/** \defgroup GenerateInspRing 
 * \ingroup inject 
 * \author S.Fairhurst
 * 
 * \brief Module for pasting a (realistic) ringdown on the end of an inspiral
 *

 *
 */
 
/** \file GenerateInspRing.h
 *  \ingroup GenerateInspRing
 * \date $Date$
 *
 * 
 */


/* includes */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/Units.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

NRCSID( GENERATEINSPRING, "$Id$");

CoherentGW *
XLALGenerateInspRing(
    CoherentGW		    *waveform,
    SimInspiralTable	*thisEvent
    );


#ifdef  __cplusplus
}                /* Close C++ protection */
#endif



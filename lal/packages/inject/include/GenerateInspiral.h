#ifndef _GENERATEINSPIRAL_H
#define _GENERATEINSPIRAL_H

#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/GenerateSpinOrbitCW.h>
#include <lal/GenerateTaylorCW.h>
#include <lal/LIGOMetadataTables.h>


#ifdef  __cplusplus
extern "C" {
#pragma }
#endif


NRCSID( GENERATEINSPIRALH, "$Id$" );

/* A reference number for the method already implemented in the injection package*/
typedef enum {  
   PPN       	= 101,
   SpinOrbitCW	= 102,
   TaylorCW  	= 103
 } Method;



/**/
void
LALInjectInspiral(
		    LALStatus        *status,
		    CoherentGW       *waveform,
		    SimInspiralTable  *params, 
		    REAL4 *buf1, 
		    REAL8 *buf2,
		    UINT4 *buf3
		    );


/**/
void
LALGetApproximantAndOrder(
			  LALStatus *status,
			  CHAR *waveform,
			  UINT4 *order,
			  UINT4 *approximant);
#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _GENERATEINSPIRAL_H */

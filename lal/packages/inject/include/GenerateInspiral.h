#ifndef _GENERATEINSPIRAL_H
#define _GENERATEINSPIRAL_H

#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/GenerateSpinOrbitCW.h>
#include <lal/GenerateTaylorCW.h>


#ifdef  __cplusplus
extern "C" {
#pragma }
#endif


NRCSID( GENERATEINSPIRALH, "$Id$" );

/* A reference number for the method already implemented in the injection package*/
typedef enum {  
   PPN       = 101,
   SpinOrbitCW      = 102,
   TaylorCW  = 103
 } Method;



/* Structure which uses the structures already defined in the inspiral package
 and injection package. */
typedef struct 
tagGeneralInspiralStruc
{
  PPNParamStruc          ppn;
  SpinOrbitCWParamStruc  socw;
  TaylorCWParamStruc     taylorcw;
  InspiralTemplate       inspiral;
  Method method;
} GeneralInspiralStruc;


/**/
void
LALGenerateInspiral(
		    LALStatus        *status,
		    CoherentGW       *waveform,
		    GeneralInspiralStruc  *params
		    );




#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _GENERATEINSPIRAL_H */

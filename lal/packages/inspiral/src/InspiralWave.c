/* 
   Interface routine needed to generate a T- or a P-approximant
   March 08, 00.
*/

#include "Inspiral.h"
#include "LALStdlib.h"

NRCSID (INSPIRALWAVEC, "$Id$");


void InspiralWave(Status *status,
		  REAL8Vector *signal,
		  InspiralTemplate *params)
{

   INITSTATUS(status, "InspiralWave", INSPIRALWAVEC);
   ATTATCHSTATUSPTR(status);

   ASSERT (signal,  status, INSPIRALWAVE_ENULL, INSPIRALWAVE_MSGENULL);
   ASSERT (signal->data,  status, INSPIRALWAVE_ENULL, INSPIRALWAVE_MSGENULL);
   ASSERT (params,  status, INSPIRALWAVE_ENULL, INSPIRALWAVE_MSGENULL);


   switch (params->domain) {
	case time:
		switch (params->method) {

			case one:
			case best:
		        	TimeDomain2(status->statusPtr, signal, params);
   	                	CHECKSTATUSPTR(status);
				break;
			case two:
				TappRpnTdomFreq(status->statusPtr, signal, params);
   	                	CHECKSTATUSPTR(status);
				break;
			case three:
				TappRpnTdomTime(status->statusPtr, signal, params);
   	                	CHECKSTATUSPTR(status);
				break;
			default:
		                fprintf(stderr,"You haven't chosen a method ... exiting\n");
                                exit(0);
				}
	break;

	case frequency:
		fprintf(stderr,"We don't have frequency domain waveforms yet \n");
			}						

   DETATCHSTATUSPTR(status);
   RETURN (status);

}

/* 
   Interface routine needed to generate a T- or a P-approximant
   March 08, 00.
*/

#include <lal/Inspiral.h>
#include <lal/LALStdlib.h>

NRCSID (INSPIRALWAVEC, "$Id$");


void LALInspiralWave(LALStatus *status,
		  REAL8Vector *signal,
		  InspiralTemplate *params)
{

   INITSTATUS(status, "LALInspiralWave", INSPIRALWAVEC);
   ATTATCHSTATUSPTR(status);

   ASSERT (signal,  status, INSPIRALWAVE_ENULL, INSPIRALWAVE_MSGENULL);
   ASSERT (signal->data,  status, INSPIRALWAVE_ENULL, INSPIRALWAVE_MSGENULL);
   ASSERT (params,  status, INSPIRALWAVE_ENULL, INSPIRALWAVE_MSGENULL);


   switch (params->domain) {
	case time:
		switch (params->method) {

			case one:
			case best:
		        	LALTimeDomain2(status->statusPtr, signal, params);
   	                	CHECKSTATUSPTR(status);
				break;
			case two:
				LALTappRpnTdomFreq(status->statusPtr, signal, params);
   	                	CHECKSTATUSPTR(status);
				break;
			case three:
				LALTappRpnTdomTime(status->statusPtr, signal, params);
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

/*------------------------------------------------------------------------------
 *    
 *    File Name: TappRpnTdomTimePhase.c
 *
 *    Created 25/08/9 
 *
 *    Authors: B. S. Sathyaprakash and D. K. Churches
 *
 *    Purpose: To calculate the phase of an inspiralling binary as a
 *             function of time up to second post-Nowtonian order	
 *
 *    Diagnostics:
 *
 *    Inputs:	*status	: pointer to the status structure
 *
 *              *params : pointer to an input parameter structure
 *
 *    Output:   *output : pointer to the output structure which
 *                        will contain the phase
 *
 *    Calls: none
 *
 *    Comments: 
 *
 *    Notes: 
 *
 *    Acknowledgements: 
 *
 *    Revision History:  
 */



#include <lal/LALStdlib.h>
#include <lal/Inspiral.h>

NRCSID (TAPPRPNTDOMTIMEPHASEC, "$Id$");


void LALTappRpnTdomTimePhase0PN (LALStatus *status,
                           InspiralwavePhaseOutput *output,
			   InspiralwavePhaseInput *params) 
{

  INITSTATUS (status, "TappRpnTdomTimePhase", TAPPRPNTDOMTIMEPHASEC);

  ASSERT(output, status, TAPPRPNTDOMTIMEPHASE_ENULL, TAPPRPNTDOMTIMEPHASE_MSGENULL);
  ASSERT(params, status, TAPPRPNTDOMTIMEPHASE_ENULL, TAPPRPNTDOMTIMEPHASE_MSGENULL);
  ASSERT(params->td > 0, status, TAPPRPNTDOMTIMEPHASE_ESIZE, TAPPRPNTDOMTIMEPHASE_MSGESIZE);

  output->phase = -(pow(params->td,0.625))
                / params->etaby2;

  RETURN(status);
}

void LALTappRpnTdomTimePhase1PN (LALStatus *status,
                           InspiralwavePhaseOutput *output,
			   InspiralwavePhaseInput *params) 
{

  INITSTATUS (status, "TappRpnTdomTimePhase", TAPPRPNTDOMTIMEPHASEC);

  ASSERT(output, status, TAPPRPNTDOMTIMEPHASE_ENULL, TAPPRPNTDOMTIMEPHASE_MSGENULL);
  ASSERT(params, status, TAPPRPNTDOMTIMEPHASE_ENULL, TAPPRPNTDOMTIMEPHASE_MSGENULL);
  ASSERT(params->td > 0, status, TAPPRPNTDOMTIMEPHASE_ESIZE, TAPPRPNTDOMTIMEPHASE_MSGESIZE);

  output->phase = -(pow(params->td,0.625)) 
                / params->etaby2;

  RETURN(status);
}

void LALTappRpnTdomTimePhase2PN (LALStatus *status,
                           InspiralwavePhaseOutput *output,
			   InspiralwavePhaseInput *params) 
{

  INITSTATUS (status, "TappRpnTdomTimePhase", TAPPRPNTDOMTIMEPHASEC);

  ASSERT(output, status, TAPPRPNTDOMTIMEPHASE_ENULL, TAPPRPNTDOMTIMEPHASE_MSGENULL);
  ASSERT(params, status, TAPPRPNTDOMTIMEPHASE_ENULL, TAPPRPNTDOMTIMEPHASE_MSGENULL);
  ASSERT(params->td > 0, status, TAPPRPNTDOMTIMEPHASE_ESIZE, TAPPRPNTDOMTIMEPHASE_MSGESIZE);

  output->phase = -(pow(params->td,0.625) 
                + params->a2*pow(params->td,0.375)) 
                / params->etaby2;

  RETURN(status);
}

void LALTappRpnTdomTimePhase3PN (LALStatus *status,
                           InspiralwavePhaseOutput *output,
			   InspiralwavePhaseInput *params) 
{

  INITSTATUS (status, "TappRpnTdomTimePhase", TAPPRPNTDOMTIMEPHASEC);

  ASSERT(output, status, TAPPRPNTDOMTIMEPHASE_ENULL, TAPPRPNTDOMTIMEPHASE_MSGENULL);
  ASSERT(params, status, TAPPRPNTDOMTIMEPHASE_ENULL, TAPPRPNTDOMTIMEPHASE_MSGENULL);
  ASSERT(params->td > 0, status, TAPPRPNTDOMTIMEPHASE_ESIZE, TAPPRPNTDOMTIMEPHASE_MSGESIZE);

  output->phase = -(pow(params->td,0.625) 
                + params->a2*pow(params->td,0.375) 
	        - params->a3*pow(params->td,0.25))
                / params->etaby2;


  RETURN(status);
}

void LALTappRpnTdomTimePhase4PN (LALStatus *status,
                           InspiralwavePhaseOutput *output,
			   InspiralwavePhaseInput *params) 
{

  INITSTATUS (status, "TappRpnTdomTimePhase", TAPPRPNTDOMTIMEPHASEC);

  ASSERT(output, status, TAPPRPNTDOMTIMEPHASE_ENULL, TAPPRPNTDOMTIMEPHASE_MSGENULL);
  ASSERT(params, status, TAPPRPNTDOMTIMEPHASE_ENULL, TAPPRPNTDOMTIMEPHASE_MSGENULL);
  ASSERT(params->td > 0, status, TAPPRPNTDOMTIMEPHASE_ESIZE, TAPPRPNTDOMTIMEPHASE_MSGESIZE);

  output->phase = -(pow(params->td,0.625) 
                + params->a2*pow(params->td,0.375) 
	        - params->a3*pow(params->td,0.25) 
                + params->a4*pow(params->td,0.125))
                / params->etaby2;


  RETURN(status);
}


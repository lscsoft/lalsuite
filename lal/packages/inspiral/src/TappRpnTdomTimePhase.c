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



#include "LALStdlib.h"
#include "Inspiral.h"

NRCSID (TAPPRPNTDOMTIMEPHASEC, "$Id$");


void TappRpnTdomTimePhase0PN (Status *status,
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

void TappRpnTdomTimePhase1PN (Status *status,
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

void TappRpnTdomTimePhase2PN (Status *status,
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

void TappRpnTdomTimePhase3PN (Status *status,
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

void TappRpnTdomTimePhase4PN (Status *status,
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


/*------------------------------------------------------------------------------
 *    
 *    File Name: TappRpnTdomFreqPhase.c
 *
 *    Created 25/08/9 
 *
 *    Authors: B. S. Sathyaprakash and D. K. Churches
 *
 *    Purpose: To calculate the phase of in inspiral waveform as a function
 *	       of frequency up to 2nd post-Newtonian order	
 *
 *    Diagnostics: Null pointer, invalid input size
 *
 *    Inputs:	*status	: pointer to our status structure
 *	        *phases : pointer to the input structure, which is of type
 *			  binarytimes
 *              *params : pointer to an input parameter structure
 *
 *    Output:   *phase  : pointer to the double precision output variable
 *
 *    Calls: none
 *
 *    Comments: See Sathyaprakash, PRD, 50, R7111, 1994, Eq.(5) or document for
 *		Tapp_RPN_tdom_freq.c for further details.
 *
 *    Notes: All frequencies are expressed in units of fs, the frequency of the
 *           waveform when it first enters the detectable part of the 
 *           detector's bandwidth
 *
 *    Acknowledgements: 
 *
 *    Revision History:  
 */


#include "LALStdlib.h"
#include "Inspiral.h"



NRCSID (TAPPRPNTDOMFREQPHASEC, "$Id$");

void LALTappRpnTdomFreqPhase0PN (LALStatus *status,
                              REAL8 *phase, 
                              InspiralPhasesInput *params) 
{
  REAL8 f;

  INITSTATUS (status, "TappRpnTdomFreqPhase", TAPPRPNTDOMFREQPHASEC);

  ASSERT(phase, status, TAPPRPNTDOMFREQPHASE_ENULL, TAPPRPNTDOMFREQPHASE_MSGENULL);
  ASSERT(params, status, TAPPRPNTDOMFREQPHASE_ENULL, TAPPRPNTDOMFREQPHASE_MSGENULL);
  ASSERT(params->f > 0, status, TAPPRPNTDOMFREQPHASE_ESIZE, TAPPRPNTDOMFREQPHASE_MSGESIZE);

  f = params->f;
  
  *phase = - params->p0 * pow(f, fiveby3)
           + params->pc;

  RETURN(status);

}

void LALTappRpnTdomFreqPhase1PN (LALStatus *status,
                              REAL8 *phase, 
                              InspiralPhasesInput *params) 
{
  REAL8 f;

  INITSTATUS (status, "TappRpnTdomFreqPhase", TAPPRPNTDOMFREQPHASEC);

  ASSERT(phase, status, TAPPRPNTDOMFREQPHASE_ENULL, TAPPRPNTDOMFREQPHASE_MSGENULL);
  ASSERT(params, status, TAPPRPNTDOMFREQPHASE_ENULL, TAPPRPNTDOMFREQPHASE_MSGENULL);
  ASSERT(params->f > 0, status, TAPPRPNTDOMFREQPHASE_ESIZE, TAPPRPNTDOMFREQPHASE_MSGESIZE);

  f = params->f;
  
  *phase = - params->p0 / pow(f, fiveby3)
           + params->pc;

  RETURN(status);
}

void LALTappRpnTdomFreqPhase2PN (LALStatus *status,
                              REAL8 *phase, 
                              InspiralPhasesInput *params) 
{
  REAL8 f;

  INITSTATUS (status, "TappRpnTdomFreqPhase", TAPPRPNTDOMFREQPHASEC);

  ASSERT(phase, status, TAPPRPNTDOMFREQPHASE_ENULL, TAPPRPNTDOMFREQPHASE_MSGENULL);
  ASSERT(params, status, TAPPRPNTDOMFREQPHASE_ENULL, TAPPRPNTDOMFREQPHASE_MSGENULL);
  ASSERT(params->f > 0, status, TAPPRPNTDOMFREQPHASE_ESIZE, TAPPRPNTDOMFREQPHASE_MSGESIZE);

  f = params->f;
  
  *phase = - params->p0 / pow(f, fiveby3)
           - params->p2 / f
           + params->pc;

  RETURN(status);
}

void LALTappRpnTdomFreqPhase3PN (LALStatus *status,
                              REAL8 *phase, 
                              InspiralPhasesInput *params) 
{
  REAL8 f;

  INITSTATUS (status, "TappRpnTdomFreqPhase", TAPPRPNTDOMFREQPHASEC);

  ASSERT(phase, status, TAPPRPNTDOMFREQPHASE_ENULL, TAPPRPNTDOMFREQPHASE_MSGENULL);
  ASSERT(params, status, TAPPRPNTDOMFREQPHASE_ENULL, TAPPRPNTDOMFREQPHASE_MSGENULL);
  ASSERT(params->f > 0, status, TAPPRPNTDOMFREQPHASE_ESIZE, TAPPRPNTDOMFREQPHASE_MSGESIZE);

  f = params->f;

  *phase = - params->p0 / pow(f,fiveby3)
           - params->p2 / f
           + params->p3 / pow(f,twoby3)
           + params->pc;


  RETURN(status);

}

void LALTappRpnTdomFreqPhase4PN (LALStatus *status,
                              REAL8 *phase, 
                              InspiralPhasesInput *params) 
{
  REAL8 f;

  INITSTATUS (status, "TappRpnTdomFreqPhase", TAPPRPNTDOMFREQPHASEC);

  ASSERT(phase, status, TAPPRPNTDOMFREQPHASE_ENULL, TAPPRPNTDOMFREQPHASE_MSGENULL);
  ASSERT(params, status, TAPPRPNTDOMFREQPHASE_ENULL, TAPPRPNTDOMFREQPHASE_MSGENULL);
  ASSERT(params->f > 0, status, TAPPRPNTDOMFREQPHASE_ESIZE, TAPPRPNTDOMFREQPHASE_MSGESIZE);

  f = params->f;

  *phase = - params->p0 / pow(f,fiveby3)
           - params->p2 / f
           + params->p3 / pow(f,twoby3)
           - params->p4 / pow(f,oneby3)
           + params->pc;

  RETURN(status);

}

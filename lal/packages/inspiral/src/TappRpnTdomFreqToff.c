/*------------------------------------------------------------------------------
 *
 *   File Name: TappRpnTdomfreqToff.c
 *
 *   Created: 25/08/99
 *
 *   Authors: B.S.Sathyaprakash and D. K. Churches
 *
 *   Purpose: To calculate the LHS of Eq.(4) in Tapp_RPN_tdom_freq_toff.ps
 *
 *   Dependencies: none
 *
 *   Inputs:	params: a pointer to the input parameters, which is of
 *                        type void *
 *		f     : the frequency in units of f_a
 *
 *   Outputs: toff is a pointer to the calculated LHS of Eq.(4) in
 *            TappRpnTdomFreqToff.tex
 *
 *   Comments: See Sathyaprakash, PRD, 50, R7111, 1994 for further details
 *
 *   Notes: All frequencies are expressed in units of fs, the frequency of the 
 *	    waveform when it first enters the detectable part of the detector's 
 *          bandwidth.
 *
 *   Acknowledgements:
 *
 *   Revision History:
 *
 */

#include <lal/LALStdlib.h>
#include <lal/Inspiral.h>



NRCSID (TAPPRPNTDOMFREQTOFFC, "$Id$");

		   
void LALTappRpnTdomFreqTofF0PN (LALStatus *status,
                             REAL8 *toff,
			     REAL8 f,
                             void *params) 
{
  InspiralToffInput *toffIn;


  INITSTATUS (status, "TappRpnTdomFreqToff", TAPPRPNTDOMFREQTOFFC);

  ASSERT(toff, status, TAPPRPNTDOMFREQTOFF_ENULL, TAPPRPNTDOMFREQTOFF_MSGENULL);
  ASSERT(params, status, TAPPRPNTDOMFREQTOFF_ENULL, TAPPRPNTDOMFREQTOFF_MSGENULL);
  ASSERT(f > 0, status, TAPPRPNTDOMFREQTOFF_ESIZE, TAPPRPNTDOMFREQTOFF_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, TAPPRPNTDOMFREQTOFF_ESIZE, TAPPRPNTDOMFREQTOFF_MSGESIZE);


  *toff = toffIn->t 
        + toffIn->t0 / pow(f, (eightby3)) 
        - toffIn->tc;


  RETURN(status);
}
		   
void LALTappRpnTdomFreqTofF1PN (LALStatus *status,
                             REAL8 *toff,
			     REAL8 f,
                             void *params) 
{
  InspiralToffInput *toffIn;


  INITSTATUS (status, "TappRpnTdomFreqToff", TAPPRPNTDOMFREQTOFFC);

  ASSERT(toff, status, TAPPRPNTDOMFREQTOFF_ENULL, TAPPRPNTDOMFREQTOFF_MSGENULL);
  ASSERT(params, status, TAPPRPNTDOMFREQTOFF_ENULL, TAPPRPNTDOMFREQTOFF_MSGENULL);
  ASSERT(f > 0, status, TAPPRPNTDOMFREQTOFF_ESIZE, TAPPRPNTDOMFREQTOFF_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, TAPPRPNTDOMFREQTOFF_ESIZE, TAPPRPNTDOMFREQTOFF_MSGESIZE);


  *toff = toffIn->t 
        + toffIn->t0 / pow(f, (eightby3)) 
        - toffIn->tc;


  RETURN(status);
}
		   
void LALTappRpnTdomFreqTofF2PN (LALStatus *status,
                             REAL8 *toff,
			     REAL8 f,
                             void *params) 
{
  InspiralToffInput *toffIn;


  INITSTATUS (status, "TappRpnTdomFreqToff", TAPPRPNTDOMFREQTOFFC);

  ASSERT(toff, status, TAPPRPNTDOMFREQTOFF_ENULL, TAPPRPNTDOMFREQTOFF_MSGENULL);
  ASSERT(params, status, TAPPRPNTDOMFREQTOFF_ENULL, TAPPRPNTDOMFREQTOFF_MSGENULL);
  ASSERT(f > 0, status, TAPPRPNTDOMFREQTOFF_ESIZE, TAPPRPNTDOMFREQTOFF_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, TAPPRPNTDOMFREQTOFF_ESIZE, TAPPRPNTDOMFREQTOFF_MSGESIZE);


  *toff = toffIn->t 
        + toffIn->t0 / pow(f, (eightby3)) 
        + toffIn->t2 / pow(f, 2.)
        - toffIn->tc;


  RETURN(status);
}
		   
void LALTappRpnTdomFreqTofF3PN (LALStatus *status,
                             REAL8 *toff,
			     REAL8 f,
                             void *params) 
{
  InspiralToffInput *toffIn;


  INITSTATUS (status, "TappRpnTdomFreqToff", TAPPRPNTDOMFREQTOFFC);

  ASSERT(toff, status, TAPPRPNTDOMFREQTOFF_ENULL, TAPPRPNTDOMFREQTOFF_MSGENULL);
  ASSERT(params, status, TAPPRPNTDOMFREQTOFF_ENULL, TAPPRPNTDOMFREQTOFF_MSGENULL);
  ASSERT(f > 0, status, TAPPRPNTDOMFREQTOFF_ESIZE, TAPPRPNTDOMFREQTOFF_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, TAPPRPNTDOMFREQTOFF_ESIZE, TAPPRPNTDOMFREQTOFF_MSGESIZE);


  *toff = toffIn->t 
        + toffIn->t0 / pow(f, (eightby3)) 
        + toffIn->t2 / pow(f, 2.)
        - toffIn->t3 / pow(f, (fiveby3)) 
        - toffIn->tc;


  RETURN(status);
}
		   
void LALTappRpnTdomFreqTofF4PN (LALStatus *status,
                             REAL8 *toff,
			     REAL8 f,
                             void *params) 
{
  InspiralToffInput *toffIn;


  INITSTATUS (status, "TappRpnTdomFreqToff", TAPPRPNTDOMFREQTOFFC);

  ASSERT(toff, status, TAPPRPNTDOMFREQTOFF_ENULL, TAPPRPNTDOMFREQTOFF_MSGENULL);
  ASSERT(params, status, TAPPRPNTDOMFREQTOFF_ENULL, TAPPRPNTDOMFREQTOFF_MSGENULL);
  ASSERT(f > 0, status, TAPPRPNTDOMFREQTOFF_ESIZE, TAPPRPNTDOMFREQTOFF_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, TAPPRPNTDOMFREQTOFF_ESIZE, TAPPRPNTDOMFREQTOFF_MSGESIZE);


  *toff = toffIn->t 
        + toffIn->t0 / pow(f, (eightby3)) 
        + toffIn->t2 / pow(f, 2.)
        - toffIn->t3 / pow(f, (fiveby3)) 
        + toffIn->t4 / pow(f, (fourby3))
        - toffIn->tc;


  RETURN(status);
}

/* -----------------------------------
 *     File Name: FindChirpTDTVeto.c
 *     Author: S. Babak
 *------------------------------------
 */

#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>
#include <lal/FindChirpTDTemplate.h>

NRCSID (FINDCHIRPTDTVETOC, "$Id$");

void LALFindChirpTDTVeto(
   LALStatus			*status,
   REAL4Vector			*chisqVec,
   FindChirpChisqInput		*input,
   FindChirpChisqTDTParams	*params
){

  /*** !!!!!!!!!!!!!!!

       I need to check that that memory is properly aalocated for
       FindChirpChisqTDTParams, which should be the case sinse this structure
       contains old one + InspiralTemplate

    ****/

  INITSTATUS( status, "LALFindChirpTDTVeto", FINDCHIRPCHISQC );
  ATTATCHSTATUSPTR( status );

  /*  ---------------------------------------------------------
  
  1. We need to find chi^2 boundaries. Note that ChisqBinVec must be null
  but memory is allocated. Here I will follow FindChirpSPData

  2. To implement classical chi^2 is straight forward following FindChirpChisq
  with one serious exeption: when we computed crosscorrelation we didn't put to
  zero waveform beyond fFinal. Should we compute boundaries for chi^2 up to
  fFinal or we should go beyond up to Nyquist frequency? Otherwise we can put to
  zero waveform after some frequency, say fFinal + df.
  
   Check UWM implementation, I think they looking for boundaries up to Nyquist, 
   but waveform is chopped off at Flso.

  3. Implement chi^2 (modified) which currently used in GEO++. Think if we can
  maximize over phase at the end of chi^2 (similar to snr computation).

  4. Implement D-test (modified)


  */

    /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );

}

/*----------------------------------------------------------------------- 
 * 
 * File Name: TfrWv.c
 * 
 * Author: Chassande-Mottin, E.
 * 
 * Revision: $Id: 
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * TfrWv
 * 
 * SYNOPSIS 
 *
 *
 * DESCRIPTION 
 * Compute the spectrogram of a given signal
 * Performs foward and inverse real FFTs on vectors and sequences of vectors.
 * 
 * DIAGNOSTICS 
 *
 * CALLS
 * 
 * NOTES
 * 
 * This code has been inspired from the Time-Frequency Toolbox (originally
 * developed by F. Auger, P. Flandrin, P. Goncalves and O. Lemoine. See
 * http://crttsn.univ-nantes.fr/~auger/tftb.html for details) and its translation 
 * in C (written by M. Davy and E. Leroy). 
 *
 *-----------------------------------------------------------------------
 */

#include "TimeFreq.h"

#define MIN(A, B)       ((A) < (B) ? (A) : (B))

NRCSID (TFRWVC, "$Id$");


void TfrWv (Status *stat, REAL4Vector* sig, TimeFreqRep *tfr, TimeFreqParam *param)
{

  INT4    nf;
  INT4    time;
  INT4    column, row;
  INT4    taumax, tau;
  REAL4Vector  *lacf = NULL;      /* local autocorrelation function */
  COMPLEX8Vector  *vtmp = NULL;
  RealFFTPlan  *plan = NULL;

  INITSTATUS (stat, "TfrWv", TFRWVC);
  ATTATCHSTATUSPTR (stat);
  
  /* Make sure the arguments are not NULL: */
  ASSERT (sig, stat, TFR_ENULL, TFR_MSGENULL);
  ASSERT (tfr, stat, TFR_ENULL, TFR_MSGENULL);
  ASSERT (param, stat, TFR_ENULL, TFR_MSGENULL);

  /* Make sure the data pointers are not NULL: */
  ASSERT (sig->data, stat, TFR_ENULL, TFR_MSGENULL);
  ASSERT (tfr->timeInstant, stat, TFR_ENULL, TFR_MSGENULL);
  ASSERT (tfr->freqBin, stat, TFR_ENULL, TFR_MSGENULL);
  ASSERT (tfr->map, stat, TFR_ENULL, TFR_MSGENULL);

  /* Make sure the requested TFR type corresponds to what will be done */
  ASSERT (tfr->type == WignerVille , stat, TFR_ENAME, TFR_MSGENAME);
  ASSERT (param->type == WignerVille , stat, TFR_ENAME, TFR_MSGENAME);

  /* Make sure the number of freq bins is a positive number: */
  nf = tfr->fRow;
  ASSERT (nf > 0 , stat, TFR_EFROW, TFR_MSGEFROW);

  /* Make sure the number of freq bins is a power of 2: */
  while(!(nf & 1))
      nf = nf>>1;

  ASSERT (nf == 1, stat, TFR_EFROW, TFR_MSGEFROW);
 
  /* Make sure the timeInstant indicates existing time instants */
  for (column=0 ; column<tfr->tCol ; column++)
    {
      if ((tfr->timeInstant[column] < 0) || (tfr->timeInstant[column] > (sig->length-1)))
	{
	  ASSERT (tfr->timeInstant[column] > 0, stat, TFR_EBADT, TFR_MSGEBADT);
	  ASSERT (tfr->timeInstant[column] < sig->length, stat, TFR_EBADT, TFR_MSGEBADT);
	}
    }
  
  TRY(EstimateFwdRealFFTPlan(stat->statusPtr, &plan, tfr->fRow), stat);
  TRY(DestroyRealFFTPlan(stat->statusPtr, &plan), stat);
  TRY(MeasureFwdRealFFTPlan(stat->statusPtr, &plan, tfr->fRow), stat);

  TRY(SCreateVector(stat->statusPtr, &lacf, tfr->fRow), stat);
  TRY(CCreateVector(stat->statusPtr, &vtmp, tfr->fRow/2+1), stat);
  
  for (column = 0; column < tfr->tCol; column++)
    {

      for (row = 0; row < tfr->fRow; row++)
	lacf->data[row] = 0.0;
      
      time = tfr->timeInstant[column];
      taumax = MIN (time, (sig->length -1 - time));
      taumax = MIN (taumax, (tfr->fRow / 2 - 1));
      
      for (tau = -taumax; tau <= taumax; tau++)
	{
	  row = (tfr->fRow+tau)%tfr->fRow;
	  lacf->data[row] =   sig->data[time + tau]*sig->data[time - tau];
        }
      
      tau=tfr->fRow/2;
      if ((time<=sig->length-tau-1)&(time>=tau))
	lacf->data[tau] =  sig->data[time+tau]*sig->data[time-tau];
      
      FwdRealFFT(stat->statusPtr, vtmp, lacf, plan);   
      
      for (row = 0; row < (tfr->fRow/2+1); row++)
	tfr->map[column][row] = vtmp->data[row].re;

    }

  for (row = 0; row < (tfr->fRow/2+1) ; row++)
    tfr->freqBin[row] = (REAL4) row / tfr->fRow;
  
  TRY(SDestroyVector(stat->statusPtr, &lacf), stat);
  TRY(CDestroyVector(stat->statusPtr, &vtmp), stat);

  TRY(DestroyRealFFTPlan(stat->statusPtr, &plan), stat);

  DETATCHSTATUSPTR (stat);

  /* normal exit */
  RETURN (stat);
}

#undef MIN

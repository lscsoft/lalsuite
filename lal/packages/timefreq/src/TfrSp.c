/*----------------------------------------------------------------------- 
 * 
 * File Name: TfrSp.c
 * 
 * Author: Chassande-Mottin, E.
 * 
 * Revision: $Id: 
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * TfrSp
 * 
 * SYNOPSIS 
 *
 *
 * DESCRIPTION 
 * Compute the spectrogram of a given signal
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

#include <lal/TimeFreq.h>
#include <lal/RealFFT.h>

#define MIN(A, B)       ((A) < (B) ? (A) : (B))

NRCSID (TFRSPC, "$Id$");


void LALTfrSp (LALStatus *stat, REAL4Vector* sig, TimeFreqRep *tfr, TimeFreqParam *param)
{

  INT4    nf;
  INT4    time;
  INT4    column, row;
  INT4    taumin, taumax, tau;
  INT4    hwl;                    /* half window length */
  REAL4Vector  *windSig = NULL;   /* windowed signal */
  REAL4   win, normh;             /* L^2 norm of the window */
  REAL4Vector  *ptmp = NULL;
  RealFFTPlan  *plan = NULL;

  INITSTATUS (stat, "LALTfrSp", TFRSPC);
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
  ASSERT (param->windowT, stat, TFR_ENULL, TFR_MSGENULL);

  /* Make sure the requested TFR type corresponds to what will be done */
  ASSERT (tfr->type == Spectrogram , stat, TFR_ENAME, TFR_MSGENAME);
  ASSERT (param->type == Spectrogram , stat, TFR_ENAME, TFR_MSGENAME);

  /* Make sure the number of freq bins is a positive number: */
  nf = tfr->fRow;
  ASSERT (nf > 0 , stat, TFR_EFROW, TFR_MSGEFROW);

  /* Make sure the number of freq bins is a power of 2: */
  while(!(nf & 1))
      nf = nf>>1;

  ASSERT (nf == 1, stat, TFR_EFROW, TFR_MSGEFROW);
  
  /* Make sure the window length is a odd number */
  ASSERT (param->windowT->length%2 != 0, stat, TFR_EWSIZ, TFR_MSGEWSIZ);

  /* Make sure the window length is smaller than the number of freq bins: */
  ASSERT ((INT4)param->windowT->length < tfr->fRow, stat, TFR_EWSIZ, TFR_MSGEWSIZ);
 
  /* Make sure the timeInstant indicates existing time instants */
  for (column=0 ; column<tfr->tCol ; column++)
    {
      if ((tfr->timeInstant[column] < 0) || (tfr->timeInstant[column] > (INT4)(sig->length-1)))
	{
	  ASSERT (tfr->timeInstant[column] > 0, stat, TFR_EBADT, TFR_MSGEBADT);
	  ASSERT (tfr->timeInstant[column] < (INT4)sig->length, stat, TFR_EBADT, TFR_MSGEBADT);
	}
    }
  
  /* ??TRY(LALEstimateFwdRealFFTPlan(stat->statusPtr, &plan, tfr->fRow), stat);*/
  /*??*/TRY(LALCreateForwardRealFFTPlan(stat->statusPtr, &plan,(UINT4)tfr->fRow,0),stat);
  /*TRY(LALDestroyRealFFTPlan(stat->statusPtr, &plan), stat);*/
  /* ??TRY(LALMeasureFwdRealFFTPlan(stat->statusPtr, &plan, tfr->fRow), stat);*/

  TRY(LALSCreateVector(stat->statusPtr, &ptmp, tfr->fRow/2 + 1), stat);
  TRY(LALSCreateVector(stat->statusPtr, &windSig, tfr->fRow), stat);

  hwl = (param->windowT->length - 1) / 2.0;
  
  for (column = 0; column < tfr->tCol; column++)
    {
      for (row = 0; row < tfr->fRow; row++)
	windSig->data[row] = 0.0;
      
      time = tfr->timeInstant[column];
      
      taumin = MIN (tfr->fRow / 2, hwl);
      taumin = MIN (taumin, time);
      
      taumax = MIN ((tfr->fRow / 2 - 1), hwl);
      taumax = MIN (taumax, (sig->length - 1.0 - time));
      
      normh = 0.0;
      for (row = -taumin; row <= taumax; row++)
	{
	  win = param->windowT->data[hwl + row];  
	  normh = normh + win * win;
	}
      normh = sqrt (normh);

      for (tau = -taumin; tau <= taumax; tau++)
	{
	  row = (tfr->fRow+tau)%(tfr->fRow);
	  windSig->data[row] = sig->data[time + tau]
	                      * param->windowT->data[hwl + tau]/normh;
	}

      /*       printf("%d\n",column); */
      LALRealPowerSpectrum (stat->statusPtr, ptmp, windSig, plan);
      /* CHANGE BACK TO ORIGINAL NORMALIZATION -- JC */
      {
        REAL4Vector *myvector = ptmp;
        UINT4 mybin;
        for ( mybin = 1; mybin < myvector->length - 1; ++mybin )
          myvector->data[mybin] *= 0.5;
      }


      
      for (row = 0; row < (tfr->fRow/2 +1); ++row)
	tfr->map[column][row] =  ptmp->data[row];
      
      /*       for (row = 0; row < tfr->fRow; row++) */
      /*  	tfr->map[column][row] = windSig->data[row]; */
    }
  
  for (row = 0; row < tfr->fRow/2+1; row++)
    tfr->freqBin[row] = (REAL4) row / tfr->fRow;
  
  TRY(LALSDestroyVector(stat->statusPtr, &ptmp), stat);
  TRY(LALSDestroyVector(stat->statusPtr, &windSig), stat);
  TRY(LALDestroyRealFFTPlan(stat->statusPtr, &plan), stat);

  DETATCHSTATUSPTR (stat);

  /* normal exit */
  RETURN (stat);
}

#undef MIN

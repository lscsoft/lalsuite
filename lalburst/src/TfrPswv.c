 /*
  * Copyright (C) 2004, 2005 Cristina V. Torres
  *                          E. Chassande-Mottin
  *
  *  This program is free software; you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation; either version 2 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU General Public License for more details.
  *
  *  You should have received a copy of the GNU General Public License
  *  along with with program; see the file COPYING. If not, write to the
  *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
  *  MA  02111-1307  USA
  */
/*-----------------------------------------------------------------------
 *
 * File Name: TfrPswv.c
 *
 * Author: Chassande-Mottin, E.
 * Maintainer: Torres, C (Univ TX at Browsville)
 *
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 * TfrPswv
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

#include <lal/TimeFreq.h>

#define MIN(A, B)       ((A) < (B) ? (A) : (B))

void LALTfrPswv (LALStatus *stat, REAL4Vector* sig, TimeFreqRep *tfr, TimeFreqParam *param)
{
  INT4    nf;
  INT4    time;
  INT4    column, row;
  INT4    taumax, tau;
  INT4    mumin, mumax, mu;
  INT4    hwlT,hwlF;                  /* half window length */
  REAL4Vector  *lacf = NULL;          /* local autocorrelation function */
  COMPLEX8Vector  *vtmp = NULL;       /* buffer */
  REAL4   normT, normF;               /* L^2 norm of the window */
  REAL4   R0, R1, R2;

  RealFFTPlan  *plan = NULL;

  INITSTATUS(stat);
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
  ASSERT (param->windowF, stat, TFR_ENULL, TFR_MSGENULL);

  /* Make sure the requested TFR type corresponds to what will be done */
  ASSERT (tfr->type == PSWignerVille , stat, TFR_ENAME, TFR_MSGENAME);
  ASSERT (param->type == PSWignerVille , stat, TFR_ENAME, TFR_MSGENAME);

  /* Make sure the number of freq bins is a positive number: */
  nf = tfr->fRow;
  ASSERT (nf > 0 , stat, TFR_EFROW, TFR_MSGEFROW);

  /* Make sure the number of freq bins is a power of 2: */
  while(!(nf & 1))
      nf = nf>>1;

  ASSERT (nf == 1, stat, TFR_EFROW, TFR_MSGEFROW);

  /* Make sure the window length is a odd number */
  ASSERT (param->windowT->length%2 != 0, stat, TFR_EWSIZ, TFR_MSGEWSIZ);
  ASSERT (param->windowF->length%2 != 0, stat, TFR_EWSIZ, TFR_MSGEWSIZ);

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
  TRY(LALCreateForwardRealFFTPlan(stat->statusPtr, &plan,(UINT4)tfr->fRow,0),stat);

  TRY(LALCCreateVector(stat->statusPtr, &vtmp, tfr->fRow/2 + 1), stat);
  TRY(LALSCreateVector(stat->statusPtr, &lacf, tfr->fRow), stat);

  hwlT = (param->windowT->length - 1) / 2.0;
  hwlF = (param->windowF->length - 1) / 2.0;
  normF = param->windowF->data[hwlF];

  for (column = 0; column < tfr->tCol; column++)
    {

     for (row = 0; row < tfr->fRow; row++)
       lacf->data[row] = 0.0;

      time = (int) tfr->timeInstant[column];

      mumin = MIN(hwlT,(INT4)(sig->length-time-1));
      mumax = MIN(hwlT,time);

      normT = 0.0;
      for(row = -mumin ; row <= mumax ; row++)
	  normT = normT + param->windowT->data[hwlT+row];

      R0 = 0.0;
      for(mu = -mumin ; mu <= mumax ; mu++)
          R0 = R0 + sig->data[time-mu] * sig->data[time-mu]
                          * param->windowT->data[hwlT+mu]/normT;

      lacf->data[0] = R0;

      taumax = MIN((time+hwlT),(INT4)(sig->length-time-1+hwlT));
      taumax = MIN(taumax,(tfr->fRow / 2 - 1));
      taumax = MIN(taumax, hwlF);

      for (tau = 1; tau <= taumax; tau++)
	{

	  mumin = MIN(hwlT,(INT4)(sig->length-time-1-tau));
	  mumax = MIN(hwlT,time-tau);

	  normT = 0.0;
	  for(row = -mumin ; row <= mumax ; row++)
	    normT = normT + param->windowT->data[hwlT+row];

	  R1 = 0.0;
	  R2 = 0.0;
	  for(mu = -mumin ; mu <= mumax ; mu++)
	    {
	      R1 = R1 + sig->data[time+tau-mu] * sig->data[time-tau-mu]
		* param->windowT->data[hwlT+mu]/normT;
	      R2 = R2 + sig->data[time-tau-mu] * sig->data[time+tau-mu]
		* param->windowT->data[hwlT+mu]/normT;
	    }

	  lacf->data[tau] = R1 * param->windowF->data[hwlF+tau]/normF;
	  lacf->data[tfr->fRow-tau] = R2 * param->windowF->data[hwlF-tau]/normF;

        }

     tau = tfr->fRow/2;
     if ((time<=(INT4)sig->length-tau-1)&(time>=tau)&(tau<=hwlF))
       {
	 mumin = MIN(hwlT,(INT4)(sig->length-time-1-tau));
	 mumax = MIN(hwlT,time-tau);

	 normT = 0.0;
	 for(row = -mumin ; row <= mumax ; row++)
	     normT = normT + param->windowT->data[hwlT+row];

	 R1 = 0.0;
	 R2 = 0.0;
	 for(mu = -mumin ; mu <= mumax ; mu++)
	   {
	     R1 = R1 + sig->data[time+tau-mu] * sig->data[time-tau-mu]
	       * param->windowT->data[hwlT+mu]/normT;
	     R2 = R2 + sig->data[time-tau-mu]* sig->data[time+tau-mu]
	       * param->windowT->data[hwlT+mu]/normT;
	   }

	 lacf->data[tau] = 0.5*(R1 * param->windowF->data[hwlF+tau]
			  + R2 * param->windowF->data[hwlF-tau])/normF;
       }

     LALForwardRealFFT (stat->statusPtr, vtmp, lacf, plan);

     for (row = 0; row < tfr->fRow/2+1 ; row++)
       tfr->map[column][row]= crealf(vtmp->data[row]);

    }

  /* Reflecting the effective halving in WV distro */
  for (row = 0; row < tfr->fRow/2+1; row++)
    tfr->freqBin[row] = (REAL4) row / (2 * tfr->fRow);

  TRY(LALCDestroyVector(stat->statusPtr, &vtmp), stat);
  TRY(LALSDestroyVector(stat->statusPtr, &lacf), stat);

  TRY(LALDestroyRealFFTPlan(stat->statusPtr, &plan), stat);

  DETATCHSTATUSPTR (stat);

  /* normal exit */
  RETURN (stat);
}


#undef MIN

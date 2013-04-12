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
 * File Name: TfrRsp.c
 *
 * Author: Chassande-Mottin, E.
 * Maintainer: Torres, C (Univ of TX at Brownsville)
 *
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 * TfrRsp
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

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/TimeFreq.h>

/* well, better macros than these one are welcome! */
#define MIN(A, B)       ((A) < (B) ? (A) : (B))
#define MAX(A, B)       ((A) > (B) ? (A) : (B))
#define ABS(A)          ((A) > 0.0 ? (A) : -(A))
#define SGN(A)          ((A) > 0.0 ? 1.0 : -1.0)
#define ROUND1(x)       (((((x)-(INT4)(x)) >= 0.0) && (((x)-(INT4)(x)) <= 0.50)) ? ((INT4)(x)) : ((INT4)(x+1)))
#define ROUND(x)        ((INT4)(SGN((x))*ROUND1(ABS((x)))))

void LALTfrRsp (LALStatus *stat, REAL4Vector* sig, TimeFreqRep *tfr, TimeFreqParam *param)
{

  INT4    nf;
  INT4    time;
  INT4    stepTime;
  INT4    column, row;
  INT4    taumin, taumax, tau;
  INT4    hwl;                     /* half window length */
  REAL4   win, normH;              /* L^2 norm of the window */
  REAL4Vector  *windowT = NULL;    /* t.window */
  REAL4Vector  *windowD = NULL;    /* d window /dt */
  REAL4Vector  *windSigH = NULL;   /* windowed signal */
  REAL4Vector  *windSigT = NULL;   /* windowed signal */
  REAL4Vector  *windSigD = NULL;   /* windowed signal */
  COMPLEX8Vector  *vtmpH = NULL;
  COMPLEX8Vector  *vtmpT = NULL;
  COMPLEX8Vector  *vtmpD = NULL;
  RealFFTPlan  *plan = NULL;
  REAL4         normhatf;          /* normalization factor */
  REAL4         eps, modulus;
  REAL4         hatt, hatf;        /* reassignment operator */
  INT4          indext,indexf;     /* reassignment index */

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

  /* Make sure the requested TFR type corresponds to what will be done */
  ASSERT (tfr->type == RSpectrogram, stat, TFR_ENAME, TFR_MSGENAME);
  ASSERT (param->type == RSpectrogram , stat, TFR_ENAME, TFR_MSGENAME);

  /* Make sure the window length is a odd number */
  ASSERT (param->windowT->length%2 != 0, stat, TFR_EWSIZ, TFR_MSGEWSIZ);

  /* Make sure the number of freq bins is a positive number: */
  nf = tfr->fRow;
  ASSERT (nf > 0 , stat, TFR_EFROW, TFR_MSGEFROW);

  /* Make sure the number of freq bins is a power of 2: */
  while(!(nf & 1))
      nf = nf>>1;

  ASSERT (nf == 1, stat, TFR_EFROW, TFR_MSGEFROW);

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

  TRY(LALSCreateVector(stat->statusPtr, &windowT, param->windowT->length), stat);
  TRY(LALSCreateVector(stat->statusPtr, &windowD, param->windowT->length), stat);
  TRY(LALCCreateVector(stat->statusPtr, &vtmpH, tfr->fRow/2 + 1), stat);
  TRY(LALCCreateVector(stat->statusPtr, &vtmpT, tfr->fRow/2 + 1), stat);
  TRY(LALCCreateVector(stat->statusPtr, &vtmpD, tfr->fRow/2 + 1), stat);
  TRY(LALSCreateVector(stat->statusPtr, &windSigH, tfr->fRow), stat);
  TRY(LALSCreateVector(stat->statusPtr, &windSigT, tfr->fRow), stat);
  TRY(LALSCreateVector(stat->statusPtr, &windSigD, tfr->fRow), stat);

  stepTime = tfr->timeInstant[1] - tfr->timeInstant[0];

  hwl = (param->windowT->length - 1) / 2.0;

  for (column = 0; column < (INT4)param->windowT->length; column++)
    windowT->data[column] = param->windowT->data[column] * (column - hwl);

  TRY(LALDwindow(stat->statusPtr, param->windowT, windowD), stat);

  for (column = 0; column < tfr->tCol; column++)
    for (row = 0; row < (tfr->fRow/2+1); row++)
	  tfr->map[column][row] = 0.0;

  for (column = 0; column < tfr->tCol; column++)
    {
      for (row = 0; row < tfr->fRow; row++)
	{
	  windSigH->data[row] = 0.0;
	  windSigT->data[row] = 0.0;
	  windSigD->data[row] = 0.0;
	}

      time = tfr->timeInstant[column];

      taumin = MIN (tfr->fRow / 2, hwl);
      taumin = MIN (taumin, time);

      taumax = MIN ((tfr->fRow / 2 - 1), hwl);
      taumax = MIN (taumax, (sig->length - 1.0 - time));

      normH = 0.0;
      for (row = -taumin; row <= taumax; row++)
	{
	  win = param->windowT->data[hwl + row];
	  normH = normH + win * win;
	}

      normH = sqrt (normH);

      for (tau = -taumin; tau <= taumax; tau++)
	{
	  row = (tfr->fRow+tau)%(tfr->fRow);
	  windSigH->data[row] = sig->data[time + tau]
	    * param->windowT->data[hwl + tau]/normH;
	  windSigT->data[row] = sig->data[time + tau]
	    * windowT->data[hwl + tau]/normH;
	  windSigD->data[row] = sig->data[time + tau]
	    * windowD->data[hwl + tau]/normH;
	}

      LALForwardRealFFT(stat->statusPtr, vtmpH, windSigH, plan);
      LALForwardRealFFT(stat->statusPtr, vtmpT, windSigT, plan);
      LALForwardRealFFT(stat->statusPtr, vtmpD, windSigD, plan);

      normhatf = tfr->fRow / (2.0 * LAL_PI);

      eps = 0.0;
      for (row = 0; row < (INT4)sig->length - 1; row++)
	eps = eps + sig->data[row]*sig->data[row];
      eps = 1.0E-6 * eps / sig->length;

      for (row = 0; row < (tfr->fRow/2+1); row++)
	{
	  modulus = crealf(vtmpH->data[row]) * crealf(vtmpH->data[row]) +
	    cimagf(vtmpH->data[row]) * cimagf(vtmpH->data[row]);
	  if  (modulus > eps)
	    {

	      hatt =  crealf(vtmpT->data[row]) * crealf(vtmpH->data[row]) +
		      cimagf(vtmpT->data[row]) * cimagf(vtmpH->data[row]);
 	      hatt = hatt / modulus / stepTime;
              indext = time + ROUND(hatt);
	      indext = MAX(indext,0);
	      indext = MIN(indext,tfr->tCol - 1);
	      hatf = crealf(vtmpD->data[row]) * cimagf(vtmpH->data[row]) -
		cimagf(vtmpD->data[row]) * crealf(vtmpH->data[row]);
	      hatf = hatf / modulus * normhatf;
	      indexf = row - ROUND(hatf);
 	      indexf = ((indexf)%(tfr->fRow)+ tfr->fRow)%(tfr->fRow);
	      indexf = (tfr->fRow)/2 - ABS((tfr->fRow)/2 - indexf);
	      if ((row==0)||(row==tfr->fRow/2))
		tfr->map[indext][indexf] = tfr->map[indext][indexf] + modulus;
	      else
		if ((indexf==0)||(indexf==tfr->fRow/2))
		  tfr->map[indext][indexf] = tfr->map[indext][indexf] + 2.0 * modulus;
		else
		  tfr->map[indext][indexf] = tfr->map[indext][indexf] + modulus;
	    }
	}
    }

  for (row = 0; row < tfr->fRow/2 + 1; row++)
    tfr->freqBin[row] = (REAL4) row / tfr->fRow;

  TRY(LALSDestroyVector(stat->statusPtr, &windowT), stat);
  TRY(LALSDestroyVector(stat->statusPtr, &windowD), stat);
  TRY(LALCDestroyVector(stat->statusPtr, &vtmpH), stat);
  TRY(LALCDestroyVector(stat->statusPtr, &vtmpT), stat);
  TRY(LALCDestroyVector(stat->statusPtr, &vtmpD), stat);
  TRY(LALSDestroyVector(stat->statusPtr, &windSigH), stat);
  TRY(LALSDestroyVector(stat->statusPtr, &windSigT), stat);
  TRY(LALSDestroyVector(stat->statusPtr, &windSigD), stat);
  TRY(LALDestroyRealFFTPlan(stat->statusPtr, &plan), stat);

  DETATCHSTATUSPTR (stat);

  /* normal exit */
  RETURN (stat);
}

#undef MIN
#undef MAX
#undef ABS
#undef SGN
#undef ROUND1
#undef ROUND

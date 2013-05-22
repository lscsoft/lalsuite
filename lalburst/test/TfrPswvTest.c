/*
*  Copyright (C) 2007 Bernd Machenschalk, Cristina Valeria Torres, Jolien Creighton
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
 * File Name: TfrPswvTest.c
 *
 * Maintainer: Torres C, (Univ of TX at Brownsville)
 * Author: Chassande-Mottin, E.
 *
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 *   main()
 *
 * SYNOPSIS
 *
 * DESCRIPTION
 *   Compute the pseudo-smoothed Wigner-Ville Distribution of a test signal
 *   Test of TfrPswv.c
 *
 * DIAGNOSTICS
 *
 * CALLS
 *
 * NOTES
 *
 *-----------------------------------------------------------------------
 */


#include <lal/TimeFreq.h>


int main(void)
{
  const INT4 Nsignal=16;
  const INT4 NwindowT=3;
  const INT4 NwindowF=5;
  const INT4 Nfft=8;

  static LALStatus status;

  REAL4Vector  *signalvec = NULL;
  CreateTimeFreqIn tfrIn;
  TimeFreqRep  *tfr = NULL;
  TimeFreqParam *param = NULL;

  INT4 column;
  INT4 row;


  /*--------------------------------------------------------------------*/

  LALSCreateVector(&status, &signalvec, Nsignal);

  signalvec->data[0]=1.0;
  for (column = 0; column < (INT4)signalvec->length; column++)
    signalvec->data[column]=(rand() % 10) / 2.0;

  /*     signalvec->data[column] = 1.0 - signalvec->data[column-1]; */
  /*     signalvec->data[column] = 1.0; */

  /*--------------------------------------------------------------------*/

  tfrIn.type=PSWignerVille;
  tfrIn.fRow=Nfft;
  tfrIn.tCol=Nsignal;
  tfrIn.wlengthT=NwindowT;
  tfrIn.wlengthF=NwindowF;

  /*--------------------------------------------------------------------*/

  LALCreateTimeFreqRep(&status, &tfr, &tfrIn);

  for (column = 0; column < tfr->tCol; column++)
    tfr->timeInstant[column]=column;

  LALCreateTimeFreqParam(&status, &param, &tfrIn);

  for (column = 0; column < (INT4)param->windowT->length; column++)
    param->windowT->data[column]=1.0;

  for (column = 0; column < (INT4)param->windowF->length; column++)
    param->windowF->data[column]=1.0;

  /*--------------------------------------------------------------------*/

  LALTfrPswv(&status,signalvec,tfr,param);
  REPORTSTATUS(&status);

  /*--------------------------------------------------------------------*/

  printf("Signal:\n");
  for (column= 0; column < (INT4)signalvec->length; column++)
    printf("%1.1f ",signalvec->data[column]);
  printf("\n\n");

  printf("TFR:\n");
  for (row= 0; row < (tfr->fRow/2+1); row++)
    {
    for (column= 0; column < tfr->tCol; column++)
      printf("%2.2f ",tfr->map[column][row]);
    printf("\n");
    }

  /*--------------------------------------------------------------------*/

  LALSDestroyVector(&status,&signalvec);
  LALDestroyTimeFreqRep(&status,&tfr);
  LALDestroyTimeFreqParam(&status,&param);

  LALCheckMemoryLeaks();

  return 0;
}


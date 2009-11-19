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
 * File Name: Dwindow.c
 *
 * Maintainer: Torres, C. (Univ TX at Browsville)
 * Author: Chassande-Mottin, E.
 *
 * Revision: $Id:
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
 * Compute the time derivative of a given window.
 * Used by LALTfrRsp (reassigned spectrogram).
 *
 * DIAGNOSTICS
 *
 * CALLS
 *
 * NOTES
 *
 * This code has been inspired from the Time-Frequency Toolbox (originally
 * developed by F. Auger, P. Flandrin, P. Goncalves and O. Lemoine. See
 * http://www.??.fr/~auger/tftb.html for details) and its translation
 * in C (written by M. Davy and E. Leroy).
 *
 *-----------------------------------------------------------------------
 */

#include <lal/TimeFreq.h>

NRCSID (DWINDOWC, "$Id$");

void LALDwindow (LALStatus *stat, REAL4Vector* window, REAL4Vector* dwindow)
{
  INT4      column, hwl;
  REAL4     step, ramp, dwin1, dwin2;

  INITSTATUS (stat, "LALDwindow", DWINDOWC);
  ATTATCHSTATUSPTR (stat);

  /* Make sure the arguments are not NULL: */
  ASSERT (dwindow, stat, TFR_ENULL, TFR_MSGENULL);
  ASSERT (window, stat, TFR_ENULL, TFR_MSGENULL);

  /* Make sure the arguments are not pointing to the same memory location */
  ASSERT (dwindow != window, stat, TFR_ESAME, TFR_MSGESAME);

  /* Make sure the data pointers are not NULL: */
  ASSERT (dwindow->data, stat, TFR_ENULL, TFR_MSGENULL);
  ASSERT (window->data, stat, TFR_ENULL, TFR_MSGENULL);

  /* Make sure the window length is a odd number */
  ASSERT (window->length%2 != 0, stat, TFR_EWSIZ, TFR_MSGEWSIZ);

  /* Make sure the window lengths are the same */
  ASSERT (dwindow->length == window->length, stat, TFR_EWSIZ, TFR_MSGEWSIZ);

  hwl = (window->length - 1) / 2.0;

  step = (window->data[0] + window->data[window->length - 1]) / 2.0;
  ramp = (window->data[window->length - 1] - window->data[0]) / (window->length - 1);

  for (column = 0 ; column < (INT4)window->length ; column++)
    dwindow->data[column] = window->data[column] - step - ramp * (-hwl + column);

  dwin1 = 0.0;
  for (column = 0 ; column < (INT4)(window->length - 1) ; column++)
    {
      dwin2 = (dwindow->data[column+1] - dwin1) / 2.0 + ramp;
      dwin1 = dwindow->data[column];
      dwindow->data[column] = dwin2;
    }

  dwindow->data[window->length - 1] = (- dwin1) / 2.0 + ramp - step;
  dwindow->data[0] = dwindow->data[0] + step;

  DETATCHSTATUSPTR (stat);

  /* normal exit */
  RETURN (stat);
}

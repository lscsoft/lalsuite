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
 * File Name: CreateTimeFreqParam.c
 *
 * Maintainer: Torres C, (Univ TX at Browsville)
 * Author: Chassande-Mottin, E.
 *
 * Revision:
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 * CreateTimeFreqParam
 *
 * SYNOPSIS
 * void LALCreateTimeFreqParam (LALStatus *, TimeFreqParam **param, CreateTimeFreqIn *in);
 *
 * DESCRIPTION
 * Create a TimeFreqParam object.
 *
 * DIAGNOSTICS
 *
 * CALLS
 * LALMalloc
 *
 * NOTES
 *
 *-----------------------------------------------------------------------
 */


#include <lal/TimeFreq.h>

NRCSID (CREATETIMEFREQPARAMC, "$Id$");

void LALCreateTimeFreqParam (LALStatus *status,
			TimeFreqParam **param,
			CreateTimeFreqIn *in)
{
  /* Initialize status */
  INITSTATUS (status, "LALCreateTimeFreqParam", CREATETIMEFREQPARAMC);

  /* Check input structure: report if NULL */
  ASSERT (in != NULL, status, CREATETFP_ENULL, CREATETFP_MSGENULL);

  /* Check return structure  */
  ASSERT (param != NULL, status, CREATETFP_ENULL, CREATETFP_MSGENULL);
  ASSERT (*param == NULL, status, CREATETFP_ENNUL, CREATETFP_MSGENNUL);

  *param = (TimeFreqParam *) LALMalloc(sizeof(TimeFreqParam));

  /* Check Allocation */
  ASSERT (*param != NULL, status, CREATETFP_EMALL, CREATETFP_MSGEMALL);

  (*param)->type=Undefined; /* Undefined until storage allocated */
  (*param)->windowT = NULL; /* NULL data until allocated */
  (*param)->windowF = NULL; /* NULL data until allocated */

  switch (in->type) {
  case Spectrogram :

    /* Make sure the window length is a odd number */
    ASSERT (in->wlengthT%2 != 0, status, CREATETFP_EWSIZ, CREATETFP_MSGEWSIZ);

    LALSCreateVector(status,&(*param)->windowT,in->wlengthT);
    (*param)->type = Spectrogram;

    break;
  case WignerVille :

    (*param)->type = WignerVille;

    break;
  case PSWignerVille :

    /* Make sure the window length is a odd number */
    ASSERT (in->wlengthT%2 != 0, status, CREATETFP_EWSIZ, CREATETFP_MSGEWSIZ);

    /* Make sure the window length is a odd number */
    ASSERT (in->wlengthF%2 != 0, status, CREATETFP_EWSIZ, CREATETFP_MSGEWSIZ);

    LALSCreateVector(status,&(*param)->windowT,in->wlengthT);
    LALSCreateVector(status,&(*param)->windowF,in->wlengthF);
    (*param)->type = PSWignerVille;

    break;
  case RSpectrogram :

    /* Make sure the window length is a odd number */
    ASSERT (in->wlengthT%2 != 0, status, CREATETFP_EWSIZ, CREATETFP_MSGEWSIZ);

    LALSCreateVector(status,&(*param)->windowT,in->wlengthT);
    (*param)->type = RSpectrogram;

    break;
  default :
    ABORT(status,CREATETFP_ETYPE, CREATETFP_MSGETYPE);

  }

  /* We be done: Normal exit */

  RETURN (status);
}

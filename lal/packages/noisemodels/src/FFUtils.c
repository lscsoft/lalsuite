/*
*  Copyright (C) 2007 Stas Babak
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

/* ****************************************************************
 * Author: Babak, S., B. Krishnan
 * ************************************************************** */
#include "lal/FFUtils.h"

/* --- version information --- */
NRCSID( FFUTILSC, "$Id$");
RCSID(  "$Id$");

#define CVS_ID_STRING_C      "$Id$"
#define CVS_REVISION_C      "$Revision$"

/* Compute the length of the longest template */

/*!!!!!!!! Need to change randIn to bankIn */


void GetMaximumTemplateSize(LALStatus *status, RandomInspiralSignalIn  randIn,
		      UINT4 *checkLength)
{

   INITSTATUS(status, "BEGetMaximumSize", BANKEFFICIENCYNEWC);
   ATTACHSTATUSPTR(status);

   InspiralTemplate dummy;

   dummy = randIn.param;
   dummy.massChoice = m1Andm2;
   dummy.mass1 = signalMinMass;
   dummy.mass2 = signalMinMass;
   dummy.approximant = EOB;
   *checkLength = 0;

   LALInspiralWaveLength(status->statusPtr, checkLength, dummy);

   DETACHSTATUSPTR(status);
   RETURN(status);

}


/*
*  Copyright (C) 2007 Badri Krishnan
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
 * File Name: HoughValidateAM.h
 *
 * Authors: Krishnan, B. 
 *
*-----------------------------------------------------------------------
 */
/*
 *   Protection against double inclusion (include-loop protection)
 *     Note the naming convention!
 */

#ifndef _HOUGHVALIDATEAMH_H
#define _HOUGHVALIDATEAMH_H

#include <time.h>
#include <math.h>

#include <lal/Random.h>
#include <lal/AVFactories.h>
#include <lal/LALStdlib.h>
#include <lal/PulsarDataTypes.h>
#include <lal/SFTfileIO.h>
#include <lal/UserInput.h>
#include <lal/GeneratePulsarSignal.h> 
#include <lal/LALComputeAM.h>
#include <lal/ComputeSky.h>
#include "../src/DriveHoughColor.h" /* proper path*/
#include "../src/MCInjectComputeHough.h"
/******************************************************
 *   Protection against C++ name mangling
 */

#ifdef  __cplusplus
extern "C" {
#endif

/* ***************************************************************
 * Constant Declarations.  Default parameters.
 */



/* ****************************************************** */

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif


#endif     /* Close double-include protection */

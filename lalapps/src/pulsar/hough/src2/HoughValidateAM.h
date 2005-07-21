/*-----------------------------------------------------------------------
 *
 * File Name: HoughValidateAM.h
 *
 * Authors: Krishnan, B. 
 *
 * Revision: $Id$
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

/******************************************************
 *  Assignment of Id string using NRCSID()
 */

NRCSID (HOUGHVALIDATEAMH, "$Id$");



/* ***************************************************************
 * Constant Declarations.  Default parameters.
 */

/*INT4 lalDebugLevel=1; */


/* ****************************************************** */

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif


#endif     /* Close double-include protection */

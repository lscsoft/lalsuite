/*
*  Copyright (C) 2007 Badri Krishnan, Alicia Sintes Olives
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
 * File Name: MCInjectComputeHough.h
 *
 * Authors: Sintes, A.M., Krishnan, B. 
 *
*-----------------------------------------------------------------------
 */
/*
 *   Protection against double inclusion (include-loop protection)
 *     Note the naming convention!
 */

#ifndef _MCINJECTCOMPUTEHOUGH_H
#define _MCINJECTCOMPUTEHOUGH_H

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
#include "./DriveHoughColor.h" /* proper path*/

#ifdef TIMING
#include "./timer/cycle_counter/Intel/GCC/cycle_counter.h" /* proper path*/
#endif

/******************************************************
 *   Protection against C++ name mangling
 */

#ifdef  __cplusplus
extern "C" {
#endif

/* ************************************************************
 * Usage format string. (Similar to DriveHoughColor_velo.c)
 */





/* ***************************************************************
 * Constant Declarations.  Default parameters.
 */


#define EARTHEPHEMERIS "./earth00-04.dat"
#define SUNEPHEMERIS "./sun00-04.dat"

#define ACCURACY 0.00000001 /* of the velocity calculation */
#define MAXFILES 3000 /* maximum number of files to read in a directory */
#define MAXFILENAMELENGTH 256 /* maximum # of characters  of a SFT filename */

#define IFO 2         /*  detector, 1:GEO, 2:LLO, 3:LHO */
#define THRESHOLD 1.6 /* thresold for peak selection, with respect to the
                              the averaged power in the search band */
#define F0 250.0          /*  frequency to build the LUT and start search */
#define FBAND 2.0          /* search frequency band  (in Hz) */
#define ALPHA 0.0		/* center of the sky patch (in radians) */
#define DELTA  (-LAL_PI_2)
#define PATCHSIZEX (LAL_PI*0.99) /* patch size */
#define PATCHSIZEY (LAL_PI*0.99)
#define NFSIZE  21 /* n-freq. span of the cylinder, to account for spin-down
                          search */
#define BLOCKSRNGMED 101 /* Running median window size */
#define NH0 1 /* number of h0 values to be anlyzed */
#define H0MIN 1.0e-23
#define NMCLOOP 10 /* number of Monte-Carlos */

/* #define SFTDIRECTORY "/nfs/morbo/geo600/hannover/sft/S2-LIGO/S2_L1_Funky-v3Calv5DQ30MinSFTs/" */
#define SFTDIRECTORY "/home/badkri/L1sfts/"
#define FILEOUT "./HoughMC"      /* prefix file output */

/* to be removed ? */
#define FILEVELOCITY "./velocity.data"  /* name: file with time-velocity info */
#define FILETIME "./Ts" /* name: file with timestamps */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

/* ******************************************************************
 *  Structure, enum, union, etc., typdefs.
 */


typedef struct tagHoughInjectParams{
  REAL8       h0;
  REAL8       fmin;   /* first_search_frequency_in_Hz */
  REAL8       fSearchBand;  /* search_band_in_Hz */
  REAL8	      deltaF; /* frequency resolution */
  UCHAR       fullSky; /* full sky 1, little patch 0 */
  REAL8       alpha;  /* patch center in equatorial coordinates (in radians) */
  REAL8       delta;
  REAL8       patchSizeAlpha;
  REAL8       patchSizeDelta;
  UINT2       pixelFactor; /* default 2, Width of the thinnest annulus in terms of pixels*/
  REAL8       vTotC;    /* estimate value of v-total/C as VTOT */
  REAL8       timeObs;
  REAL8Vector spnFmax;
} HoughInjectParams;

typedef struct tagHoughPulsarTemplate{
  REAL8        f0;
  REAL8        latitude;   /* of the source in radians */
  REAL8        longitude;  /* of the source in radians */
  REAL8Vector  spindown;   /* SpinOrder and parameters */ 
} HoughPulsarTemplate;

typedef struct tagPulsarData{
  REAL8        f0;
  REAL8        latitude;   /* of the source in radians */
  REAL8        longitude;  /* of the source in radians */
  REAL4        aPlus;
  REAL4        aCross;
  REAL4        psi;
  REAL8        phi0;
  REAL8Vector  spindown;   /* SpinOrder and parameters */
} PulsarData;


/*
 *  Functions Declarations (i.e., prototypes). Not declared in DriveHoughColor.h */
			     
void GenerateInjectTemplateParams(LALStatus   *status,
                        PulsarData           *injectPulsar,
                        HoughPulsarTemplate  *templatePulsar,
                        HoughInjectParams    *params);

/* ****************************************************** */

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif


#endif     /* Close double-include protection */

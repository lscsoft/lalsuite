/*-----------------------------------------------------------------------
 *
 * File Name: MCInjectComputeHough.h
 *
 * Authors: Sintes, A.M., Krishnan, B. 
 *
 * Revision: $Id$
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

/******************************************************
 *  Assignment of Id string using NRCSID()
 */

NRCSID (MCINJECTCOMPUTEHOUGHH, "$Id$");

/* ************************************************************
 * Usage format string. (Similar to DriveHoughColor_velo.c)
 */


#define USAGE "Usage: \n\
 [-d debuglevel] \n\
        Default value:  lalDebugLevel=1\n\
 [-i IFO (1,2,3)] \n\
        Interferometer option\n\
             1:GEO600 (Default), 2:LLO, 3:LHO  \n\
 [-E Earth ephemeris data filename] \n\
         Default value:\n\
	 /afs/aeiw/grawave/Linux/lal/lal/packages/pulsar/test/earth03.dat\n\
 [-S Sun ephemeris data filename] \n\
	 Default value:\n\
	 /afs/aeiw/grawave/Linux/lal/lal/packages/pulsar/test/sun03.dat\n\
 [-D directory] \n\
        Directory where the input SFT data are located. \n\
	If not set, the program will look into ./data1\n\
 [-o outfile-basename] \n\
        This is a string that prefixes some output filenames.\n\
        It might contain a path. Filenames are formed by \n\
        appending _<number>.m and _par\n\
 	If not set, the program will write into ./HoughMC.<number>\n\
 [-V time-velocity data file]\n\
        This is a string of the output time-velocity data file.\n\
        It might contain a path.\n\
        If not set, the program will write on ./velocity.data \n\
 [-f first search frequency (in Hz)] \n\
        Lowest search frequency in Hz. \n\
        Default: 250.0 Hz\n\
 [-b search frequency band (in Hz)] \n\
        Bandwith to be analyzed \n\
        Default: 2.0 Hz\n\
 [-t peak threshold selection ] \n\
        Threshold relative to the PSD for the selection of peak in the\n\
        time-frequency plane \n\
        Default: 1.6 (almost optimal)\n\
 [-w running median window size ] \n\
        To estimate the psd \n\
 	Default: 101\n\
 [-p alpha delta (in radians)] \n\
        Center of the sky patch (in radians) to be analysed.\n\
        Default: alpha = 0.0, delta = - pi/2 \n\
 [-s patchSizeAlpha patchSizeDelta (in radians)]\n\
        if not specified, full sky will be used\n\
 [-H nh0 h0Min h0Max]\n\
	Number of h0 values to be anlyzed and  interval h0Min h0Max\n\
	Default nh0 =1 , h0Min 10^-23\n\
 [-L nMCloop]\n\
	Number of injected signals\n\
	Default: 10  \n\
\n"



/* ***************************************************************
 * Constant Declarations.  Default parameters.
 */

INT4 lalDebugLevel=1;
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

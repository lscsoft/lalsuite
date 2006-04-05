/*-----------------------------------------------------------------------
 *
 * File Name: MCInjectHoughMULTI.h
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

#ifndef _MCINJECTHOUGHMULTI_H
#define _MCINJECTHOUGHMULTI_H

#include <time.h>
#include <math.h>

#include <lal/Random.h>
#include <lal/AVFactories.h>
#include <lal/LALStdlib.h>
#include <lal/PulsarDataTypes.h>
#include <lal/SFTfileIO.h>
#include <lal/UserInput.h>
#include <lal/GeneratePulsarSignal.h> 
#include <lal/SFTClean.h>
#include "./DriveHoughColor.h" /* proper path*/


/******************************************************
 *   Protection against C++ name mangling
 */

#ifdef  __cplusplus
extern "C" {
#endif

/******************************************************
 *  Assignment of Id string using NRCSID()
 */

NRCSID (MCINJECTHOUGHMULTIH, "$Id$");

/* ************************************************************
 * Usage format string. 
 */


/* ***************************************************************
 * Constant Declarations.  Default parameters.
 */


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

typedef struct tagHoughTemplate{
  REAL8        f0;
  REAL8        latitude;   /* of the source in radians */
  REAL8        longitude;  /* of the source in radians */
  REAL8Vector  spindown;   /* SpinOrder and parameters */ 
} HoughTemplate;

typedef struct tagHoughNearTemplates{
  REAL8        f0[2]; /* f0 values */
  REAL8        f1[2]; /* 1st spindown parameters */
  REAL8UnitPolarCoor skytemp[4]; /* near sky template parameters */
} HoughNearTemplates;


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
			     
void GenerateInjectParams(LALStatus  *status,
                        PulsarData           *injectPulsar,
                        HoughTemplate        *templatePulsar,
			HoughNearTemplates   *closeTemplates,
                        HoughInjectParams    *params,
			LineNoiseInfo        *lines  );

void ComputeFoft(LALStatus   *status,
                 REAL8Vector          *foft,
                 HoughTemplate        *pulsarTemplate,
                 REAL8Vector          *timeDiffV,
                 REAL8Cart3CoorVector *velV,
		 REAL8                timeBase);


void PrintLogFile (LALStatus *status, 
                  CHAR *dir, 
		  CHAR *basename, 
		  CHAR *skyfile, 
		  CHAR *executable );

/* ****************************************************** */

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif


#endif     /* Close double-include protection */

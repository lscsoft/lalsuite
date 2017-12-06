/*
*  Copyright (C) 2007 Gregory Mendell, Virginia Re
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

/*********************************************************************************/
/*                                                                               */
/* File: ComputeStackSlideSums.h                                                 */
/* Purpose: Drive StackSlide jobs under Condor and the Grid                      */
/* Author: Mendell, G. and Landry, M.                                            */
/* Started: 03 December 2003                                                     */
/*                                                                               */
/*********************************************************************************/

/* REVISIONS: */
/* 01/23/04 gam; Increase MAXFILES from 40000 to 80000 to allow running on more SFTs */
/*               (Should allocate memory for this in the future.)                    */
/* 05/07/04 gam; add alternative to using glob */
/* 05/11/04 gam; Add code to software inject signals into the SFTs for Monte Carlo simulations */
/* 05/28/04 gam; Use LALUniformDeviate from LAL utilities package (include <lal/Random.h>) to generate random mismatch during Monte Carlo. */
/* 04/12/05 gam; Change default Monte Carlo to use RunStackSlideIsolatedMonteCarloSimulation in StackSlideIsolated.c. */

#ifndef _COMPUTESTACKSLIDESUMS_H
#define _COMPUTESTACKSLIDESUMS_H

/*********************************************/
/*                                           */
/* START SECTION: include header files       */
/*                                           */
/*********************************************/
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
/* #include <glob.h> */
/* 05/07/04 gam; add alternative to using glob; include dirent.h instead of glob.h. */
#include <dirent.h>
#include <time.h>
#include <lal/ComputeSky.h>
#include <lal/LALDatatypes.h>
#include <lal/LALDemod.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/Date.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
/* Next two are used when calibrating BLKs */
/* #include <ExtractSeries.h> */
#include <lal/VectorOps.h>
/* #include <lal/RealFFT.h> */ /* 05/11/04 gam; not needed */
#include <lal/LALConstants.h>
#include <lal/Units.h>
#include <errno.h>
#include "DriveStackSlide.h"
/* 04/12/05 gam; RunStackSlideIsolatedMonteCarloSimulation is in StackSlideIsolated.c. */
#include "StackSlideIsolated.h"
#include "StackSlideBinary.h"
/*********************************************/
/*                                           */
/* END SECTION: include header files         */
/*                                           */
/*********************************************/

/*********************************************/
/*                                           */
/* START SECTION: define eternal global vars */
/*                                           */
/*********************************************/
extern char *optarg;
extern int optind, opterr, optopt;
/*********************************************/
/*                                           */
/* END SECTION: define eternal global vars   */
/*                                           */
/*********************************************/

/*********************************************/
/*                                           */
/* START SECTION: define constants           */
/*                                           */
/*********************************************/
/* #define MAXFILES 40000 */    /* Maximum # of files in a directory  */
#define MAXFILES 80000          /* 01/23/04 gam */
#define MAXFILENAMELENGTH 256   /* Maximum # of characters of a SFT filename */
/*********************************************/
/*                                           */
/* START SECTION: define constants           */
/*                                           */
/*********************************************/

/*********************************************/
/*                                           */
/* START SECTION: define macros              */
/*                                           */
/*********************************************/
#define INTERNAL_CHECKSTATUS_FROMMAIN(status) \
  if (status.statusCode) {                    \
     fprintf(stderr,"Error: statusCode = %i statusDescription = %s \n", status.statusCode, status.statusDescription); \
     return 1; \
  }
/*********************************************/
/*                                           */
/* END SECTION: define macros                */
/*                                           */
/*********************************************/
  
/*********************************************/
/*                                           */
/* START SECTION: define structs             */
/*                                           */
/*********************************************/
typedef struct GlobalVariablesTag {
  double df;
  double tsft;
  INT4 SFTno;
  INT4 nsamples;
  INT4 ifmin;     /* Index of first frequency to get from SFT */
  INT4 ifmax;     /* Index of last frequency to get from SFT */
  INT4 Ti;        /* GPS seconds of first SFT */
  INT4 Tf;        /* GPS seconds of last SFT */
  char filelist[MAXFILES][MAXFILENAMELENGTH];
} GlobalVariables;
  
struct headertag {
    REAL8 endian;
    INT4  gps_sec;
    INT4  gps_nsec;
    REAL8 tbase;
    INT4  firstfreqindex;
    INT4  nsamples;
} header;
/*********************************************/
/*                                           */
/* END SECTION: define structs               */
/*                                           */
/*********************************************/

/*********************************************/
/*                                           */
/* START SECTION: prototype declarations     */
/*                                           */
/*********************************************/
int ReadSFTData(StackSlideSearchParams *params);
int SetGlobalVariables(StackSlideSearchParams *params);
int Freemem(StackSlideSearchParams *params);
/*********************************************/
/*                                           */
/* END SECTION: prototype declarations       */
/*                                           */
/*********************************************/

#endif /* _COMPUTESTACKSLIDESUMS_H */

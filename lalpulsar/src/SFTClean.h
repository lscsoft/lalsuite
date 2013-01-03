/*
 *  Copyright (C) 2005 Badri Krishnan, Alicia Sintes, Greg Mendell
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
#ifndef _SFTCLEAN_H
#define _SFTCLEAN_H

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \defgroup SFTClean_h Header SFTClean.h
 * \ingroup pkg_SFTIO
 * \author Badri Krishnan, Alicia Sintes, Greg Mendell
 *
 * \brief Header file for cleaning routines

Routines for cleaning SFT files using known spectral disturbances.

\heading{Synopsis}

\code
#include <lal/SFTClean.h>
\endcode

Format for list of known spectral disturbances and using
them to clean SFT data

\heading{Error conditions}

 Test program. %%

*/
/*@{*/

/* REVISIONS: */
/* 09/09/05 gam; make RandomParams *randPar a parameter for CleanCOMPLEX8SFT. Thus only need to */
/*               initialze RandomParams *randPar once and avoid repeatly opening /dev/urandom.  */


/*
 * 5. Includes. This header may include others; if so, they go immediately
 *    after include-loop protection. Includes should appear in the following
 *    order:
 *    a. Standard library includes
 *    b. LDAS includes
 *    c. LAL includes
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/SFTfileIO.h>
#include <lal/Random.h>
#include <lal/PulsarDataTypes.h>
#include <lal/UserInput.h>
#include <lal/LUT.h>
#include <lal/RngMedBias.h>

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>

/*
 * 7. Error codes and messages. This must be auto-extracted for
 *    inclusion in the documentation.
 */

/**\name Error Codes */
/*@{*/
#define SFTCLEANH_ENULL 1
#define SFTCLEANH_EFILE 2
#define SFTCLEANH_EHEADER 3
#define SFTCLEANH_EENDIAN 4
#define SFTCLEANH_EVAL 5
#define SFTCLEANH_ELINENAME 6
#define SFTCLEANH_ESEEK 9
#define SFTCLEANH_EREAD 10
#define SFTCLEANH_EWRITE 11

#define SFTCLEANH_MSGENULL "Null pointer"
#define SFTCLEANH_MSGEFILE "Could not open file"
#define SFTCLEANH_MSGEHEADER "Incorrect header in file"
#define SFTCLEANH_MSGEENDIAN "Incorrect endian type"
#define SFTCLEANH_MSGEVAL  "Invalid value"
#define SFTCLEANH_MSGELINENAME  "Invalid linefile name"
#define SFTCLEANH_MSGESEEK "fseek failed"
#define SFTCLEANH_MSGEREAD "fread failed"
#define SFTCLEANH_MSGEWRITE "fwrite failed"
/*@}*/


/* ******************************************************
 * 8. Macros. But, note that macros are deprecated.
 *    They could be moved to the modules where are needed
 */

/* *******************************************************
 * 9. Constant Declarations. (discouraged)
 */



/* **************************************************************
 * 10. Structure, enum, union, etc., typdefs.
 */


/** structure for storing list of spectral lines -- constructed by expanding list of harmonics*/
  typedef struct tagLineNoiseInfo{
    INT4         nLines;     /**< number of lines */
    REAL8        *lineFreq;  /**< central frequency of the line in Hz */
    REAL8        *leftWing;  /**< width to the left from central frequency in Hz */
    REAL8        *rightWing; /**< width to the right in Hz */
  } LineNoiseInfo;

  /** structure for storing the contents of the input list of known
      spectral disturbances */
  typedef struct tagLineHarmonicsInfo{
    INT4         nHarmonicSets; /**< number of sets of harmonics */
    REAL8        *startFreq;    /**< starting frequency of set in Hz */
    REAL8        *gapFreq;      /**< frequency difference between adjacent harmonics in Hz */
    INT4         *numHarmonics; /**< Number of harmonics */
    REAL8        *leftWing;     /**< width to the left of each line in set in Hz */
    REAL8        *rightWing;    /**< width to the right in Hz */
  } LineHarmonicsInfo;

/*
 * 11. Extern Global variables. (discouraged)
 */


/*
 * 12. Functions Declarations (i.e., prototypes).
 */


void LALFindNumberHarmonics (LALStatus           *status,
			  LineHarmonicsInfo   *harmonicInfo,
			  CHAR                *fname
			  );

void  LALReadHarmonicsInfo (LALStatus          *status,
			 LineHarmonicsInfo  *lineInfo,
			 CHAR               *fname
			 );

void  LALHarmonics2Lines (LALStatus          *status,
		       LineNoiseInfo      *lineInfo,
		       LineHarmonicsInfo  *harmonicsInfo
		       );

void LALChooseLines (LALStatus        *status,
		  LineNoiseInfo    *outLine,
		  LineNoiseInfo    *inLine,
		  REAL8            freqMin,
		  REAL8            freqMax
		  );


void LALCheckLines ( LALStatus           *status,
		  INT4                *flag,
		  LineNoiseInfo       *lines,
		  REAL8               freq);


void LALFindNumberLines (LALStatus        *status,
		      LineNoiseInfo    *lineInfo,
		      CHAR             *fname
		      );

void LALReadLineInfo (LALStatus        *status,
		   LineNoiseInfo  *lineInfo,
		   CHAR           *fname
		   );

void LALCleanCOMPLEX8SFT (LALStatus          *status,
		       SFTtype            *sft,
		       INT4               width,
		       INT4               window,
		       LineNoiseInfo      *lineInfo,
		       RandomParams       *randPar);


void LALCleanSFTVector (LALStatus          *status,
			SFTVector          *sftVect,
			INT4               width,
			INT4               window,
			LineNoiseInfo      *lineInfo,
			RandomParams       *randPar);

void LALCleanMultiSFTVect (LALStatus          *status,
			   MultiSFTVector     *multVect,
			   INT4               width,
			   INT4               window,
			   LineNoiseInfo      *lineInfo,
			   RandomParams       *randPar);

void LALRemoveKnownLinesInSFTVect (LALStatus   *status,
				   SFTVector   *sftVect,
				   INT4        width,
				   INT4        window,
				   CHAR        *linefile,
				   RandomParams *randPar);

void LALRemoveKnownLinesInMultiSFTVector (LALStatus        *status,
					  MultiSFTVector   *multiSFTVect,
					  INT4             width,
					  INT4             window,
					  LALStringVector *linefiles,
					  RandomParams     *randPar);


/*@}*/

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif     /* Close double-include protection _SFTCLEAN_H */

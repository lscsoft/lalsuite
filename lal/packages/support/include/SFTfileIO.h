/*
 * Copyright (C) 2004, 2005 R. Prix, B. Machenschalk, A.M. Sintes
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

/** \defgroup SFTfileIO
 * \ingroup support
 * \author R. Prix, B. Machenschalk, A.M. Sintes
 * \date $Date$
 * \brief IO-Module for reading/writing SFTs (Short Fourier transform data-files)
 *
 */
 
/** \file 
 * \ingroup SFTfileIO
 * \date $Date$
 * \brief Header file defining the API for the SFTfileIO modules.
 *
 *
 * Routines for reading and writing SFT binary files.
 * 
 */

#ifndef _SFTFILEIO_H  	/* Double-include protection. */
#define _SFTFILEIO_H

/* includes */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/PulsarDataTypes.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

NRCSID (SFTFILEIOH, "$Id$");

/** \name Error codes */
/*@{*/
#define SFTFILEIOH_ENULL 	1
#define SFTFILEIOH_EFILE 	2
#define SFTFILEIOH_EHEADER 	3
#define SFTFILEIOH_EVERSION 	4
#define SFTFILEIOH_EVAL 	5
#define SFTFILEIOH_EENDIAN 	6
#define SFTFILEIOH_ENONULL 	12
#define SFTFILEIOH_EFREQBAND 	13
#define SFTFILEIOH_EMEM 	14
#define SFTFILEIOH_EGLOB 	15
#define SFTFILEIOH_ECOPYSIZE	16
#define SFTFILEIOH_EDIFFLENGTH  17

#define SFTFILEIOH_MSGENULL 	"Null pointer"
#define SFTFILEIOH_MSGEFILE 	"Error in file-IO"
#define SFTFILEIOH_MSGEHEADER 	"Incorrect header in file"
#define SFTFILEIOH_MSGEVERSION 	"This SFT-version is not currently supported"
#define SFTFILEIOH_MSGEVAL  	"Invalid value"
#define SFTFILEIOH_MSGEENDIAN 	"Wrong endian encoding of SFT (not supported yet"
#define SFTFILEIOH_MSGENONULL  	"Output pointer not NULL"
#define SFTFILEIOH_MSGEFREQBAND "Required frequency-band is not in SFT"
#define SFTFILEIOH_MSGEMEM 	"Out of memory"
#define SFTFILEIOH_MSGEGLOB 	"Failed to get filelist from directory/pattern"
#define SFTFILEIOH_MSGECOPYSIZE	"Target SFT-struct has not enough frequency-bins for copying"
#define SFTFILEIOH_MSGEDIFFLENGTH "Sorry, can only read SFTs of identical length (currently)"
/*@}*/


/** This structure contains the header-info contained in an SFT-file of specification 
 * version v1.0.
 */
typedef struct tagSFTHeader {
  REAL8  version;		/**< SFT version-number (currently only 1.0 allowed )*/
  INT4   gpsSeconds;		/**< gps start-time (seconds)*/
  INT4   gpsNanoSeconds;	/**< gps start-time (nanoseconds) */
  REAL8  timeBase;		/**< length of data-stretch in seconds */
  INT4   fminBinIndex;		/**< first frequency-index contained in SFT */
  INT4   length;  		/**< number of frequency bins */
} SFTHeader;


/*
 * Functions Declarations (i.e., prototypes).
 */
void LALReadSFTheader (LALStatus  *status, SFTHeader *header, const CHAR *fname); 
void LALReadSFTdata (LALStatus  *status, SFTtype *sft, const CHAR *fname, INT4 fminBinIndex);
void LALWriteSFTfile (LALStatus  *status, const SFTtype *sft, const CHAR *outfname);
void LALReadSFTfile (LALStatus *status, SFTtype **sft, REAL8 fMin, REAL8 fMax, const CHAR *fname);

void LALReadSFTfiles (LALStatus *status,
		      SFTVector **sftvect, 
		      REAL8 fMin, 
		      REAL8 fMax, 
		      UINT4 wingBins, 
		      const CHAR *fpattern);

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif     /* Close double-include protection _SFTBIN_H */
 








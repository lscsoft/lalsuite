/*-----------------------------------------------------------------------
 *
 * File Name: SFTfileIO.h
 *
 * Authors: Sintes, A.M., Krishnan, B., Prix, R., Makchenschalk, B.,
 *          inspired from Siemens, X.
 *
 * Revision: $Id$
 *
 * History:   Created by Sintes & Prix May 21, 2003
 *            Modified by Machenschalk Jun 16, 2004
 *
 *-----------------------------------------------------------------------
 */
 
/* *********************************** <lalVerbatim file="SFTfileIOHV">
Author: Sintes, A.M., Prix. R, Machenschalk, B.
$Id$
************************************* </lalVerbatim> */

   	      	   
/* <lalLaTeX>  *********************************************
\section{Header \texttt{SFTfileIO.h}}
\label{s:SFTfileIO.h}

Routines for reading and writing SFT binary files

\subsection*{Synopsis}

\begin{verbatim}
#include <lal/SFTfileIO.h>
\end{verbatim}

*************************************************</lalLaTeX> */

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

/********************************************************** <lalLaTeX>
\subsection*{Error codes}
</lalLaTeX>
***************************************************** <lalErrTable> */
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
/*************************************************** </lalErrTable> */

/****************************************************** <lalLaTeX>
\subsection*{Types}

\subsubsection*{Structure \texttt{SFTHeader}}
\idx[Type]{SFTHeader}

This structure contains the header-info contained in an SFT of specification 
version v1.0.

</lalLaTeX> */
/* <lalVerbatim> */
typedef struct tagSFTHeader {
  REAL8  version;		/* SFT version-number (currently only 1.0 allowed )*/
  INT4   gpsSeconds;		/* gps start-time */
  INT4   gpsNanoSeconds;
  REAL8  timeBase;		/* length of data-stretch in seconds */
  INT4   fminBinIndex;		/* first frequency-index contained in SFT */
  INT4   length;  		/* number of frequency bins */
} SFTHeader;
/* </lalVerbatim> */

/* just for reference: this is the FrequencySeries type:
 * {
 *   CHAR              name[LALNameLength];
 *   LIGOTimeGPS       epoch;
 *   REAL8             f0;
 *   REAL8             deltaF;
 *   LALUnit           sampleUnits
 *   COMPLEX8Sequence *data;
 * }
 * COMPLEX8FrequencySeries;
 */

/********************************************************** <lalLaTeX>
\vfill{\footnotesize\input{SFTfileIOHV}}
\newpage\input{SFTfileIOC}
\newpage\input{SFTfileIOTestC}
******************************************************* </lalLaTeX> */

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

void dump_SFT (FILE *fp, const SFTtype *sft, INT4 format);

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif     /* Close double-include protection _SFTBIN_H */
 








/*-----------------------------------------------------------------------
 *
 * File Name: SFTfileIO.h
 *
 * Authors: Sintes, A.M.,  Krishnan, B.,   Prix, R.  &inspired from Siemens, X.
 *
 * Revision: $Id$
 *
 * History:   Created by Sintes May 21, 2003
 *            Modified...
 *
 *-----------------------------------------------------------------------
 */
 
/* *********************************** <lalVerbatim file="SFTfileIOHV">
Author: Sintes, A.M., Prix. R
$Id$
************************************* </lalVerbatim> */

/*
 * 
 * Documentation:
 * 
 * The Relevant functions are the ones described here.
 * Please, have a look at the structures here defined first.
 * 
 * ReadSFTbinHeader1 (LALStatus  *status,
 *                    SFTHeader1 *header, CHAR  *fname);
 *    Reads the header from a given SFT file.		   
 * 
 * ReadSFTtype (LALStatus  *status,
 *               SFTtype *sft, CHAR *fname, INT4 fminBinIndex );	
 *    Reads from one SFT file  "fname" a COMPLEX8FrequencySeries "sft".
 *    The memory of the sft is assumed to be alocated before, i.e, you know
 *    the numer of elements you want to read, The first element that will
 *    be read will correspond to  fminBinIndex.
 *    
 * ReadCOMPLEX8SFTbinData1 (LALStatus  *status,
 *           COMPLEX8SFTData1    *sft, CHAR   *fname ); 
 *    Reads from one SFT file  "fname" a  COMPLEX8SFTData1 "sft".
 *    The memory for the sft should already be allocated.
 *    You need to specify the sft.length and sft.fminBinIndex before the call.
 *    This function will fill in the rest of the fields in the sft.
 * 	    
 */
   	      	   
/* <lalLaTeX>  *********************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Header \texttt{SFTfileIO.h}}
\label{s:SFTfileIO.h}
Routines for reading SFT binary files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection*{Synopsis}

\begin{verbatim}
#include <lal/SFTfileIO.h>
\end{verbatim}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection*{Error conditions}
\vspace{0.1in}
\input{SFTfileIOHErrorTable}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\vfill{\footnotesize\input{SFTfileIOHV}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newpage\input{SFTfileIOC}
%%%%%%%%%% Test program. %%
\newpage\input{SFTfileIOTestC}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

*************************************************</lalLaTeX> */

/*
 * 4.  Protection against double inclusion (include-loop protection)
 *     Note the naming convention!
 */

#ifndef _SFTFILEIO_H
#define _SFTFILEIO_H

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
#include <lal/PulsarDataTypes.h>

/*
 *   Protection against C++ name mangling
 */

#ifdef  __cplusplus
extern "C" {
#endif

 /*
 * 6. Assignment of Id string using NRCSID()  
 */
  
NRCSID (SFTFILEIOH, "$Id$");
  
/*
 * 7. Error codes and messages. This must be auto-extracted for 
 *    inclusion in the documentation.
 */
  
/* <lalErrTable file="SFTfileIOHErrorTable"> */
  
#define SFTFILEIOH_ENULL 1
#define SFTFILEIOH_EFILE 2
#define SFTFILEIOH_EHEADER 3
#define SFTFILEIOH_EENDIAN 4
#define SFTFILEIOH_EVAL 5
#define SFTFILEIOH_ESEEK 9
#define SFTFILEIOH_EREAD 10
#define SFTFILEIOH_EWRITE 11

#define SFTFILEIOH_ENONULL 12
#define SFTFILEIOH_EFREQBAND 13
#define SFTFILEIOH_EMEM 14

#define SFTFILEIOH_MSGENULL "Null pointer"
#define SFTFILEIOH_MSGEFILE "Could not open file"
#define SFTFILEIOH_MSGEHEADER "Incorrect header in file"
#define SFTFILEIOH_MSGEENDIAN "Incorrect endian type" 
#define SFTFILEIOH_MSGEVAL  "Invalid value"
#define SFTFILEIOH_MSGESEEK "fseek failed"
#define SFTFILEIOH_MSGEREAD "fread failed"
#define SFTFILEIOH_MSGEWRITE "fwrite failed"

#define SFTFILEIOH_MSGENONULL  "Output pointer not NULL"
#define SFTFILEIOH_MSGEFREQBAND "Required frequency-band is not in SFT"
#define SFTFILEIOH_MSGEMEM 	"Out of memory"
/* </lalErrTable>  */


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

  typedef struct tagSFTHeader {
    REAL8  endian;
    INT4   gpsSeconds;
    INT4   gpsNanoSeconds;
    REAL8  timeBase;
    INT4   fminBinIndex;
    INT4   length;  
  } SFTHeader;

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

/*
 * 11. Extern Global variables. (discouraged) 
 */
  

/*
 * 12. Functions Declarations (i.e., prototypes).
 */

void LALReadSFTheader (LALStatus  *status, SFTHeader *header, const CHAR *fname); 
void LALReadSFTtype (LALStatus  *status, SFTtype *sft, const CHAR *fname, INT4 fminBinIndex);
void LALWriteSFTtoFile (LALStatus  *status, const SFTtype *sft, const CHAR *outfname);
void LALReadSFTfile (LALStatus *status, SFTtype **sft, REAL8 fmin, REAL8 fmax, const CHAR *fname);


#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif     /* Close double-include protection _SFTBIN_H */
 








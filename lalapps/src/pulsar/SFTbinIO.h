/*-----------------------------------------------------------------------
 *
 * File Name: SFTbin.h
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
 
/* *********************************** <lalVerbatim file="SFTbinHV">
Author: Sintes, A.M., 
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
\section{Header \texttt{SFTbinIO.h}}
\label{s:SFTbinIO.h}
Routines for reading SFT binary files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection*{Synopsis}

\begin{verbatim}
#include <lal/SFTbinIO.h>
\end{verbatim}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection*{Error conditions}
\vspace{0.1in}
\input{SFTbinHErrorTable}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\vfill{\footnotesize\input{SFTbinHV}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newpage\input{SFTbinIOC}
%%%%%%%%%% Test program. %%
\newpage\input{SFTbinIOTestC}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

*************************************************</lalLaTeX> */

/*
 * 4.  Protection against double inclusion (include-loop protection)
 *     Note the naming convention!
 */

#ifndef _SFTBINIO_H
#define _SFTBINIO_H

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

/*
 *   Protection against C++ name mangling
 */

#ifdef  __cplusplus
extern "C" {
#endif

 /*
 * 6. Assignment of Id string using NRCSID()  
 */
  
NRCSID (SFTBINH, "$Id$");
  
/*
 * 7. Error codes and messages. This must be auto-extracted for 
 *    inclusion in the documentation.
 */
  
/* <lalErrTable file="SFTbinHErrorTable"> */
  
#define SFTBINH_ENULL 1
#define SFTBINH_EFILE 2
#define SFTBINH_EHEADER 3
#define SFTBINH_EENDIAN 4
#define SFTBINH_EVAL 5
#define SFTBINH_ESEEK 9
#define SFTBINH_EREAD 10
#define SFTBINH_EWRITE 11

#define SFTBINH_ENONULL 12
#define SFTBINH_EFREQBAND 13
#define SFTBINH_EMEM 14

#define SFTBINH_MSGENULL "Null pointer"
#define SFTBINH_MSGEFILE "Could not open file"
#define SFTBINH_MSGEHEADER "Incorrect header in file"
#define SFTBINH_MSGEENDIAN "Incorrect endian type" 
#define SFTBINH_MSGEVAL  "Invalid value"
#define SFTBINH_MSGESEEK "fseek failed"
#define SFTBINH_MSGEREAD "fread failed"
#define SFTBINH_MSGEWRITE "fwrite failed"

#define SFTBINH_MSGENONULL  "Output pointer not NULL"
#define SFTBINH_MSGEFREQBAND "Required frequency-band is not in SFT"
#define SFTBINH_MSGEMEM 	"Out of memory"
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

  typedef struct tagSFTHeader1{
    REAL8  endian;
    INT4   gpsSeconds;
    INT4   gpsNanoSeconds;
    REAL8  timeBase;
    INT4   fminBinIndex;
    INT4   length;  
  } SFTHeader1;
  
  typedef struct tagCOMPLEX8SFTData1{  /* simple case */
    LIGOTimeGPS  epoch; /* epoch of first series sample */
    REAL8        timeBase;
    INT4         fminBinIndex;
    INT4         length; /* number of elements in data */ 
    COMPLEX8     *data;  /* pointer to the data */
  } COMPLEX8SFTData1;

  typedef struct tagCOMPLEX8SFTvector1{  
    UINT4                length; /* number of elements  */ 
    COMPLEX8SFTData1     *sft;  /* pointer to the data */
  } COMPLEX8SFTvector1;
  
  typedef struct tagCOMPLEX16SFTData1{  /* simple case */
    LIGOTimeGPS  epoch; /* epoch of first series sample */
    REAL8        timeBase;
    INT4         fminBinIndex;
    INT4         length; /* number of elements in data */ 
    COMPLEX16     *data;  /* pointer to the data */
  } COMPLEX16SFTData1;

typedef COMPLEX8FrequencySeries         SFTtype;

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

void LALReadSFTbinHeader1 (LALStatus  *status,
                   SFTHeader1    *header,
		   CHAR          *fname
		   );

void LALReadSFTtype (LALStatus  *status,
		 SFTtype    *sft,  /* asumed  memory is allocated  */
		 CHAR       *fname,
		 INT4       fminBinIndex
		 );

void ReadCOMPLEX8SFTbinData1 (LALStatus  *status,
		   COMPLEX8SFTData1    *sft,
		   CHAR                *fname
		   );

void ReadCOMPLEX16SFTbinData1 (LALStatus  *status,
		   COMPLEX16SFTData1    *sft,
		   CHAR                 *fname
		   );
void WriteCOMPLEX8SFTbinData1 (LALStatus          *status,
		       COMPLEX8SFTData1   *sft,
		       CHAR               *outfname
		       );

void WriteCOMPLEX16SFTbinData1 (LALStatus          *status,
		       COMPLEX16SFTData1   *sft,
		       CHAR                *outfname
		       );

  /* reinhard's draft */
void ReadSFTfile (LALStatus *status, SFTtype **sft, const CHAR *fname, REAL8 fmin, REAL8 fmax);


#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif     /* Close double-include protection _SFTBIN_H */
 








/*-----------------------------------------------------------------------
 *
 * File Name: SFTClean.h
 *
 * Authors: Krishnan, B.
 *
 * Revision: $Id$
 *
 * History:   Moved from LALAPPS 31/7/05
 *            
 *
 *-----------------------------------------------------------------------
 */
 
/* *********************************** <lalVerbatim file="SFTbinHV">
Author: Krishnan, B 
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>  *********************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Header \texttt{SFTClean.h}}
\label{s:SFTClean.h}

Routines for cleaning SFT files using known spectral disturbances. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection*{Synopsis}

\begin{verbatim}
#include <lal/SFTClean.h>
\end{verbatim}

\noindent Format for list of known spectral disturbances and using 
them to clean SFT data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection*{Error conditions}
\vspace{0.1in}
\input{SFTCleanHErrorTable}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\vfill{\footnotesize\input{SFTCleanHV}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newpage\input{SFTCleanC}
%%%%%%%%%% Test program. %%
\newpage\input{SFTCleanTestC}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

*************************************************</lalLaTeX> */


/**
 * Routines for cleaning SFT files using known spectral disturbances. 
 *
 * Also contains format for list of known spectral disturbances 
 *
 **/

/*
 * 4.  Protection against double inclusion (include-loop protection)
 *     Note the naming convention!
 */

#ifndef _SFTCLEAN_H
#define _SFTCLEAN_H

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
 *   Protection against C++ name mangling
 */

#ifdef  __cplusplus
extern "C" {
#endif

 /*
 * 6. Assignment of Id string using NRCSID()  
 */
  
NRCSID (SFTCLEANH, "$Id$");
  
/*
 * 7. Error codes and messages. This must be auto-extracted for 
 *    inclusion in the documentation.
 */
  
/* <lalErrTable file="SFTbinHErrorTable"> */
  
#define SFTCLEANH_ENULL 1
#define SFTCLEANH_EFILE 2
#define SFTCLEANH_EHEADER 3
#define SFTCLEANH_EENDIAN 4
#define SFTCLEANH_EVAL 5
#define SFTCLEANH_ESEEK 9
#define SFTCLEANH_EREAD 10
#define SFTCLEANH_EWRITE 11

#define SFTCLEANH_MSGENULL "Null pointer"
#define SFTCLEANH_MSGEFILE "Could not open file"
#define SFTCLEANH_MSGEHEADER "Incorrect header in file"
#define SFTCLEANH_MSGEENDIAN "Incorrect endian type" 
#define SFTCLEANH_MSGEVAL  "Invalid value"
#define SFTCLEANH_MSGESEEK "fseek failed"
#define SFTCLEANH_MSGEREAD "fread failed"
#define SFTCLEANH_MSGEWRITE "fwrite failed"


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



  typedef struct tagLineNoiseInfo{
    INT4         nLines; /* number of lines */ 
    REAL8        *lineFreq; /* central frequency of the lines */
    REAL8        *leftWing; /* width to the left from central ferquency */
    REAL8        *rightWing; /* width to the right */
  } LineNoiseInfo; 

  typedef struct tagLineHarmonicsInfo{
    INT4         nHarmonicSets; /* number of sets of harmonics */
    REAL8        *startFreq; /* starting frequency of set */
    REAL8        *gapFreq;  /* frequency difference between adjacent harmonics */
    INT4         *numHarmonics; /* Number of harmonics */  
    REAL8        *leftWing; /* width to the left of each line in set */
    REAL8        *rightWing; /* width to the right */
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
		       LineNoiseInfo      *lineInfo
		       );



#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif     /* Close double-include protection _SFTCLEAN_H */
 








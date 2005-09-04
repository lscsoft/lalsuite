/*-----------------------------------------------------------------------
 *
 * File Name: NormalizeSFTRngMed.h
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
 **/

/*
 * 4.  Protection against double inclusion (include-loop protection)
 *     Note the naming convention!
 */

#ifndef _NORMALIZESFTRNGMED_H
#define _NORMALIZESFTRNGMED_H

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
#include <lal/LALRunningMedian.h>

/*
 *   Protection against C++ name mangling
 */

#ifdef  __cplusplus
extern "C" {
#endif

 /*
 * 6. Assignment of Id string using NRCSID()  
 */
  
NRCSID (NORMALIZESFTRNGMEDH, "$Id$");
  
/*
 * 7. Error codes and messages. This must be auto-extracted for 
 *    inclusion in the documentation.
 */
  
/* <lalErrTable file="SFTbinHErrorTable"> */
  
#define NORMALIZESFTRNGMEDH_ENULL 1
#define NORMALIZESFTRNGMEDH_EVAL 2

#define NORMALIZESFTRNGMEDH_MSGENULL "Null pointer"
#define NORMALIZESFTRNGMEDH_MSGEVAL  "Invalid value"


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





/*
 * 11. Extern Global variables. (discouraged) 
 */
  

/*
 * 12. Functions Declarations (i.e., prototypes).
 */


void LALSFTtoPeriodogram (LALStatus    *status,
			  REAL8FrequencySeries    *periodo,
			  COMPLEX8FrequencySeries *SFT);

void LALPeriodoToPSDRngMed (LALStatus  *status,
			    REAL8FrequencySeries  *psd,
			    REAL8FrequencySeries  *periodo,
			    INT4                  blockSize);

void LALNormalizeSFT (LALStatus  *status,
		      SFTtype  *sft,
		      INT4     blockSize,
		      UCHAR    normSwitch);

void LALNormalizeSFTVect (LALStatus  *status,
			  SFTVector  *sftVect,
			  INT4     blockSize,
			  UCHAR    normSwitch); /* normSwitch == 0 for frequency domain 
						   normalization and == 1 for time domain */

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif     /* Close double-include protection _NORMALIZESFTRNGMED_H */
 








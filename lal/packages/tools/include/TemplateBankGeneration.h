/*
*  Copyright (C) 2007 Bernd Machenschalk, Chad Hanna, Jolien Creighton, Benjamin Owen
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

/*________________________________________________________________________
 *
 * File Name: TemplateBankGeneration.h
 *
 * Author: Hanna C. R.
 *
 * Revision: $Id$
 *
 *________________________________________________________________________
 */

/* ------------------------------ AUTO-DOC ---------------------------- */
/* -------------------------------------------------------------------- */


/*<lalVerbatim file="TemplateBankGenerationHV">
  Author: Hanna, C. R.
  $Id$
  </lalVerbatim> */

/*SECTION - HEADER - "TemplateBankGeneration.h" -------------------- <lalLaTeX>
  \section{Header \texttt{TemplateBankGeneration.h}}
  \label{s:TemplateBankGeneration.h}
  </lalLaTeX> */

  /*SUBSECTION - SYNOPSIS ------------------------------------- <lalLaTeX>
    \subsection*{Synopsis}
    \begin{verbatim}
    #include <lal/TemplateBankGeneration.h>
    \end{verbatim}
    \noindent This header file includes all the necessary types and
    function prototypes for LALNDTemplateBank() and LALMakeTemplateBank().
    NDTemplateBank() provides a general way to tile a parameter space with
    a constant metric ( currently only for less than 12 dimensions).
    MakeTemplateBank() provides a general way for applications to
    generate a template bank with suitable I/O.
  </lalLaTeX>
    END SUBSECTION - SYNOPSIS ----------------------------------------- */

  /*<lalLaTeX>
    \vfill{\footnotesize\input{TemplateBankGenerationHV}}
  </lalLaTeX> */

  /*NEWPAGE - SUBSECTION - ERROR CODES ------------------------ <lalLaTeX>
    \newpage
    \subsection*{Error codes}
    \input{TemplateBankGenerationHE}
    </lalLaTeX>
    END SUBSECTION - ERROR CODES -------------------------------------- */


  /*<lalLaTeX>
    \vfill{\footnotesize\input{TemplateBankGenerationHV}}
  </lalLaTeX> */

  /*NEWPAGE - SUBSECTION - TYPES -------------------------------<lalLaTeX>
    \newpage
    \subsection*{Types}
    \input{TemplateBankGenerationHT-MakeTemplateBankInput}
    \input{TemplateBankGenerationHT-InspiralTmpltBankCInput}
    \input{TemplateBankGenerationHT-PulsarTmpltBankCInput}
    \input{TemplateBankGenerationHT-BurstTmpltBankCInput}
    \input{TemplateBankGenerationHT-TemplateBankType}
    \input{TemplateBankGenerationHT-NDTemplateBankInput}
    \input{TemplateBankGenerationHT-NDTemplateBankOutput}
    \idx[Type]{MakeTemplateBankInput}
    \idx[Type]{InspiralTmpltBankCInput}
    \idx[Type]{PulsarTmpltBankCInput}
    \idx[Type]{BurstTmpltBankCInput}
    \idx[Type]{TemplateBankType}
    \idx[Type]{NDTemplateBankInput}
    \idx[Type]{NDTemplateBankOutput}
    </lalLaTeX>
    END SUBSECTION TYPES ---------------------------------------------- */

  /*SUBSECTION - NOTES -----------------------------------------<lalLaTeX>
    \subsection*{Notes}
    \noindent\begin{itemize}
    \item No notes yet.
    \end{itemize}
    </lalLaTeX>
    END SUBSECTION NOTES ---------------------------------------------- */

  /*<lalLaTeX>
    \vfill{\footnotesize\input{TemplateBankGenerationHV}}
  </lalLaTeX> */

  /*NEWPAGE - MODULES INPUT - ----------------------------------------- */
    /* <lalLaTeX>
    \newpage\input{NDTemplateBankC}
    \newpage\input{MakeTemplateBankC}
    </lalLaTeX>
    END - MODULES INPUT ----------------------------------------------- */

  /*NEWPAGE - TESTS INPUT --------------------------------------------- */
    /* <lalLaTeX>
    </lalLaTeX>
    END - TESTS INPUT ------------------------------------------------- */

/*END SECTION - HEADER - "TemplateBankGeneration.h" ------------------------ */


/* -------------------------------------------------------------------- */
/* ------------------------- END AUTO-DOC ----------------------------- */

#ifndef _TEMPLATEBANKGENERATION_H
#define _TEMPLATEBANKGENERATION_H

#include<lal/LALStdlib.h>
#include<lal/LALStatusMacros.h>
#if 0
#include<lal/LALInspiral.h>
#include<lal/LALInspiralBank.h>
#endif
#include<lal/LALDatatypes.h>
#include<lal/LIGOMetadataTables.h>

NRCSID (TEMPLATEBANKGENERATIONH,"$Id$");

/* <lalErrTable file="TemplateBankGenerationHE"> */
#define TEMPLATEBANKGENERATIONH_ENULL 1
#define TEMPLATEBANKGENERATIONH_MSGENULL "Unexpected NULL pointer to an input type"
/* </lalErrTable> */


#if 0
/* Inspiral Group structure for lalapps tmpltbank.c searches 100-199 */
/* <lalVerbatim file="TemplateBankGenerationHT-InspiralTmpltBankCInput"> */
typedef struct {
  /* InspiralCoarseBankIn Parameters */
     REAL4   mMin;         	/* minimum component mass       */
     REAL4   mMax;         	/* maximum component mass       */
     REAL8   MMax; 		/* Maximum total Mass ???	*/
     REAL4   psi0Min;          	/* minimum value of psi0        */
     REAL4   psi0Max;          	/* maximum value of psi0        */
     REAL4   psi3Min;          	/* minimum value of psi3        */
     REAL4   psi3Max;          	/* maximum value of psi3        */
     REAL4   alpha;            	/* BCV amplitude correction     */
     INT4    numFcutTemplates;     /* num tmplts in fcut direction */
     REAL4   mmCoarse;         	/* minimum requested match      */
     REAL4   fLower;
     REAL4   fUpper;           	/* upper frequency cutoff       */
     Order   order;             /* post-Newtonian order         */
     REAL8FrequencySeries shf;	/* used for PSD stuff		*/
     REAL8   tSampling;  	/* Sampling Rate		*/
     REAL8   etamin;		/* Minimum Eta 			*/
     REAL8   LowGM;		/* ??????????????? 		*/
     REAL8   HighGM;		/* ???????????????		*/
     Approximant approximant;	/* approximation method         */
     CoordinateSpace space;  	/* coordinate space used        */
     InspiralBankMassRange massRange;
     INT4    ntiles;		/*number of templates made 	*/
     REAL4   mmFine;		/* not implemented? */
     INT4    iflso;		/* not implemented? */
  /* InspiralTemplate Parameters */


  /* InspiralSpinBank Additional Parameters */
     REAL4 mTwoMin;		/* minimum mass smaller body 	*/
     REAL4 mTwoMax;		/* maximum mass smaller body 	*/
     REAL4 betaMin; 		/* minimum beta spin parameter 	*/
     REAL4 betaMax;		/* maximum beta spin parameter 	*/
  } InspiralTmpltBankCInput; /* InspiralTmpltBankCInput */

/* </lalVerbatim> */

/* <lalVerbatim file="TemplateBankGenerationHT-PulsarTmpltBankCInput"> */
typedef struct {
     INT4 placeholder;
     } PulsarTmpltBankCInput; /* not used yet */

/* </lalVerbatim> */

/* <lalVerbatim file="TemplateBankGenerationHT-BurstTmpltBankCInput"> */
typedef struct {
     INT4 placeholder;
     } BurstTmpltBankCInput; /* not used yet */

/* </lalVerbatim> */

/* <lalVerbatim file="TemplateBankGenerationHT-MakeTemplateBankInput"> */
typedef union {
     InspiralTmpltBankCInput *InspiralInput; 	/* Searches 100-199 */
     PulsarTmpltBankCInput *PulsarInput;	/* Searches 200-299 */
     BurstTmpltBankCInput *BurstInput;		/* Searches 300-399 */
     } MakeTemplateBankInput;

/* </lalVerbatim> */
#endif

/* <lalVerbatim file="TemplateBankGenerationHT-TemplateBankType"> */
typedef enum {
  /* Binary Inspiral Searches 100-199 */
     BCVType,
     BCVSpinType,
     PrecessingType,
  /* Pulsar Searches 200-299 */
     Pulsar,
  /* Burst Searches 300-399 */
     Burst
  /* Other Searches 400-499 Etc... */
     } TemplateBankType;

/* </lalVerbatim> */

/* <lalVerbatim file="TemplateBankGenerationHT-NDTemplateBankInput"> */
typedef struct {
  REAL4                	minCoordinates[12]; 	/* Psi0, Psi3, etc.      */
  REAL4                	maxCoordinates[12]; 	/* Psi0, Psi3, etc.      */
  REAL4                	minParameters[12];  	/* mass, eta, etc.       */
  REAL4                	maxParameters[12];  	/* mass, eta, etc.       */
  INT2                 	dimension;          	/* 3D?, 4D? -> ND!       */
  REAL4                	mm;                 	/* mismatch              */
  TemplateBankType     	type;               	/* whats the search?     */
  REAL8FrequencySeries *PSD;		   	/* Power Spec. Density   */
  REAL4			f0;			/* Moment scaling freq	 */
  } NDTemplateBankInput;
/* </lalVerbatim> */

/* <lalVerbatim file="TemplateBankGenerationHT-NDTemplateBankOutput"> */
typedef struct NDTemplateBankOutput{
  REAL4				coordinateVals[12];
  REAL4 			parameterVals[12];
  struct NDTemplateBankOutput	*next;
  } NDTemplateBankOutput;
/* </lalVerbatim> */


typedef void (*NDTemplateBankMetricPtr)( LALStatus *, NDTemplateBankInput *, REAL4Array *);
typedef void (*NDTemplateBankTestPtr)(LALStatus *, NDTemplateBankInput *, NDTemplateBankOutput *, INT2 *);

typedef struct {
  NDTemplateBankMetricPtr 	metric; 	/*Ptr to metric function */
  NDTemplateBankTestPtr		test;		/*Ptr boundary test fct  */
  } NDTemplateBankFunctionPtrs;
 /* Function Prototypes */

void
LALNDTemplateBank(
 	LALStatus *,
   	NDTemplateBankInput *,
        NDTemplateBankFunctionPtrs *,
       	NDTemplateBankOutput **);


#if 0
void
LALMakeTemplateBank(
     	LALStatus *,
     	TemplateBankType *,
     	MakeTemplateBankInput *,
     	MetadataTable *);
     /* LALMakeTemplateBank(); */
#endif

#endif

/* <lalVerbatim file="GenerateInspiralHV">
 * Author: Cokelaer, T.
 * $Id$
 * </lalVerbatim>  */


/* <lalLaTeX>
 * 	\section{Header \texttt{GenerateInspiral.h}}
 * 	\label{s:GenerateInspiral.h}
 * 	
 * 	Header file for the inspiral injection interface code.
 * 	
 * 	\subsection*{Synopsis}
 * 	\begin{verbatim}
 * 	#include <lal/GenerateInspiral.h>
 * 	\end{verbatim}
 * </lalLaTeX> */

#ifndef _GENERATEINSPIRAL_H
#define _GENERATEINSPIRAL_H

#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/GenerateSpinOrbitCW.h>
#include <lal/GenerateTaylorCW.h>
#include <lal/LIGOMetadataTables.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif

NRCSID( GENERATEINSPIRALH, "$Id$" );

/*  <lalLaTeX>
 *  	\subsection*{Error codes}
 *  </lalLaTeX>  */

/* <lalErrTable> */
#define GENERATEINSPIRALH_ENORM 	0
#define GENERATEINSPIRALH_ENULL  	1
#define GENERATEINSPIRALH_EAPPROX       2

#define GENERATEINSPIRALH_MSGENORM 	"Normal exit"
#define GENERATEINSPIRALH_MSGENULL  	"Null pointer"
#define GENERATEINSPIRALH_MSGEAPPROX    "non valid Approximant; must be TaylorT1, TaylorT2 \
                                         TaylorT3, EOB, PadeT1 or GeneratePPN; case dependant."

/* </lalErrTable> */


/* EOB a 3PN */
#define GENERATEINSPIRAL_ZETA2       0.
#define GENERATEINSPIRAL_OMEGAS      0.
#define GENERATEINSPIRAL_THETA       0.

/* For the spinning case. might be changed later or include 
   in the injection itself ? */
#define GENERATEINSPIRAL_SOURCETHETA 1.
#define GENERATEINSPIRAL_SOURCEPHI   2.
#define GENERATEINSPIRAL_SPIN1       0.3
#define GENERATEINSPIRAL_SPIN2       0.7
#define GENERATEINSPIRAL_THETA1      0.86
#define GENERATEINSPIRAL_THETA2      0.86

#define GENERATEINSPIRAL_PHI1        1.53
#define GENERATEINSPIRAL_PHI2        1.53

/* idem for CW waveform*/
#define GENERATEINSPIRAL_F0         100.
#define GENERATEINSPIRAL_ARG         0
#define GENERATEINSPIRAL_UDOT        0.5
#define GENERATEINSPIRAL_RP          0.5
#define GENERATEINSPIRAL_E           0.   
#define GENERATEINSPIRAL_ALPHA       0.

/* A reference number for the method already 
 * implemented in the injection package. Should add PPN to 
 * the inspiral strucutre. the two others are useless here
 * */
typedef enum {  
   PPN       	= 101,
   SpinOrbitCW	= 102,
   TaylorCW  	= 103
 } Method;

/*  <lalLaTeX>
 *     \newpage\input{GenerateInspiralC}
 *  </lalLaTeX>  */
void LALGenerateInspiral(LALStatus        *status,
			 CoherentGW       *waveform,
			 SimInspiralTable *params, 
			 PPNParamStruc    *ppnParamsInputOutput );

void ComputeSpin(InspiralTemplate *params);


/* three function to read the order and approximant from a string */
void LALGenerateInspiralGetOrderFromString(LALStatus *status,
					   CHAR *message,
					   UINT4 *result);
     
void LALGenerateInspiralGetApproxFromString(LALStatus *status,
					   CHAR *message,
					   UINT4 *result);

void LALGenerateInspiralGetModelFromString(LALStatus *status,
					   CHAR *message,
					   UINT4 *order,
					   UINT4 *model);

/*  three function to populate the needed structures */
void  LALGenerateInspiralPopulatePPN(LALStatus             *status,
				     PPNParamStruc         *ppnParams,
				     SimInspiralTable      *thisEvent);

void LALGenerateInspiralPopulateInspiral(LALStatus             *status,
					 InspiralTemplate      *inspiralParams,
					 SimInspiralTable      *thisEvent,
					 PPNParamStruc         *ppnParams);

void LALGenerateInspiralPopulateInspiralSpin(LALStatus             *status,
					     InspiralTemplate      *inspiralParams);
#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _GENERATEINSPIRAL_H */

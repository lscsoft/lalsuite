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

#define GENERATEINSPIRALH_MSGENORM 	"Normal exit"
#define GENERATEINSPIRALH_MSGENULL  	"Null pointer"
/* </lalErrTable> */


/* EOB a 3PN */
#define GENERATEINSPIRAL_ZETA2       0.
#define GENERATEINSPIRAL_OMEGAS      0.
#define GENERATEINSPIRAL_THETA       0.

/* For the spinning case. might be changed later or include 
   in the injection itself ? */
#define GENERATEINSPIRAL_SOURCETHETA 1.
#define GENERATEINSPIRAL_SOURCEPHI   2.

#define GENERATEINSPIRAL_SPIN1X       1 
#define GENERATEINSPIRAL_SPIN1Y       1 
#define GENERATEINSPIRAL_SPIN1Z       1 

#define GENERATEINSPIRAL_SPIN2X       1
#define GENERATEINSPIRAL_SPIN2Y       1
#define GENERATEINSPIRAL_SPIN2Z       1

#define GENERATEINSPIRAL_ORBITTHETA0 1.5
#define GENERATEINSPIRAL_ORBITPHI0   0

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

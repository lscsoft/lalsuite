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
#define GENERATEINSPIRALH_EAPPROXIMANT  6
#define GENERATEINSPIRALH_EORDER  	7

#define GENERATEINSPIRALH_MSGENORM 		"Normal exit"
#define GENERATEINSPIRALH_MSGENULL  		"Null pointer"
#define GENERATEINSPIRALH_MSGEAPPROXIMANT  	"Non valid approximant"
#define GENERATEINSPIRALH_MSGEORDER  		"Non valid order"
/* </lalErrTable> */


/* EOB a 3PN */
#define GENERATEINSPIRAL_ZETA2       0.
#define GENERATEINSPIRAL_OMEGAS      0.
#define GENERATEINSPIRAL_THETA       0.

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

void LALGetApproximantAndOrderFromString(LALStatus *status,
					 CHAR *waveform,
					 UINT4 *order,
					 UINT4 *approximant);
#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _GENERATEINSPIRAL_H */

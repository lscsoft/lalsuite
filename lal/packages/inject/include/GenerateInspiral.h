/* <lalVerbatim file="GenerateInspiralHV">
 * Author: Cokelaer, T.
 * $Id$
 * </lalVerbatim>  */


/* <lalLaTeX>
 * 	\section{Header \texttt{GenerateInspiral.h}}
 * 	\label{s:GenerateInspiral.h}
 * 	
 * 	Header file for the inspiral injection interface code. The 
 * 	code contained in GenerateInspiral.c is an interface between the 
 * 	injection package and the inspiral package. More precisely, the 
 * 	function GenerateInspiral.c is used within the FindChirpSimulation.c
 * 	file of the FindChirp package in order to inject waveforms into real
 * 	data. The injection is done through the inject package in order to 
 * 	take into account the interferometer position, binary orientation ...
 * 	
 * 	GenerateInspiral has the capability of injecting both waveform designed 
 * 	within the inspiral package (TaylorT1, T2, T3, PadeT1, EOB, and spinning 
 * 	waveform) and the inject package (so-called PPN waveform).
 *
 * 	There is also a test code as well which allows to check the output of 
 *      the code. It is called InjectionInterfaceTest.c
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
#include <lal/SeqFactories.h>

#include <lal/Units.h>


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
#define GENERATEINSPIRALH_EDFDT 	1


#define GENERATEINSPIRALH_MSGENORM 	"Normal exit"
#define GENERATEINSPIRALH_MSGENULL  	"Null pointer"
#define GENERATEINSPIRALH_MSGEDFDT  	"Waveform sampling interval is too large"

/* </lalErrTable> */


/* parameter for the EOB at 3PN. In principle, the three 
 * following parameter should be set to zero. */
#define GENERATEINSPIRAL_ZETA2       0.
#define GENERATEINSPIRAL_OMEGAS      0.
#define GENERATEINSPIRAL_THETA       0.

/* For the spinning case. might be changed later or include 
   in the injection itself ? */
#define GENERATEINSPIRAL_SOURCETHETA 1.
#define GENERATEINSPIRAL_SOURCEPHI   2.


#define GENERATEINSPIRAL_DEFAULT_FLOWER 40

/* A reference number for the method already 
 * implemented in the injection package. Should add PPN to 
 * the inspiral strucutre. The two others are useless here a
 * */
typedef enum {  
   PPN       	= 101
 } Method;

/*  <lalLaTeX>
 *     \newpage\input{GenerateInspiralC}
 *  </lalLaTeX>  */
void LALGenerateInspiral(LALStatus        *status,
			 CoherentGW       *waveform,
			 SimInspiralTable *params, 
			 PPNParamStruc    *ppnParamsInputOutput );


/* three function to read the order and approximant from a string */
void LALGetOrderFromString(LALStatus *status,
			   CHAR *message,
			   Order *result);
     
void LALGetApproximantFromString(LALStatus *status,
			    CHAR *message,
			    Approximant *result);

/*  three function to populate the needed structures */
void  LALGenerateInspiralPopulatePPN(LALStatus             *status,
				     PPNParamStruc         *ppnParams,
				     SimInspiralTable      *thisEvent);

void LALGenerateInspiralPopulateInspiral(LALStatus             *status,
					 InspiralTemplate      *inspiralParams,
					 SimInspiralTable      *thisEvent,
					 PPNParamStruc         *ppnParams);



#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _GENERATEINSPIRAL_H */

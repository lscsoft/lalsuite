/************************************ <lalVerbatim file="BinaryPulsarTimingCV">
Author: Dupuis, R. D.
$Id$  
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{BinaryPulsarTiming.c}}

The routines in this calculate the time delay dut to the motion of a pulsar in a binary system. 
  
\subsubsection*{Prototypes}
\input{BinaryPulsarTimingCP}
\idx{LALBinaryPulsarDeltaT()}

\subsubsection*{Description}

\noindent The function \texttt{LALBinaryPulsarDeltaT}  ... 
\newline


\subsubsection*{Algorithm}

To be completed.

\subsubsection*{Uses}
\begin{verbatim}
To bo comepleted.
\end{verbatim}


\subsubsection*{Notes}

\vfill{\footnotesize\input{BinaryPulsarTimingCV}}

******************************************************* </lalLaTeX> */ 

/******* INCLUDE STANDARD LIBRARY HEADERS; ************/
/* note LALStdLib.h already includes stdio.h and stdarg.h */

/******* INCLUDE ANY LDAS LIBRARY HEADERS ************/

/******* INCLUDE ANY LAL HEADERS ************/
#include <lal/BinaryPulsarTiming.h>
#include <lal/LALConstants.h>
#include <math.h>
/******* DEFINE RCS ID STRING ************/
NRCSID( BinaryPulsarTimingC, "$Id$" );

/******* DEFINE LOCAL CONSTANTS AND MACROS ************/

/******* DECLARE LOCAL (static) FUNCTIONS ************/
/* (definitions can go here or at the end of the file) */

/******* DEFINE GLOBAL FUNCTIONS ************/

/* <lalVerbatim file="BinaryPulsarTimingCP"> */
void
LALBinaryPulsarDeltaT          ( LALStatus                      *status,
		  	    	 BinaryPulsarTiming		*output,
		   	    	 REAL8				tgps,
				 BinaryPulsarParameters 	*params )
/* </lalVerbatim> */
{
  REAL8 U,du;
  REAL8 phase;
  
  
  INITSTATUS( status, "LALBinaryPulsarDeltaT", BinaryPulsarTimingC );
  ATTATCHSTATUSPTR (status); 

 /******* CHECK VALIDITY OF ARGUMENTS  ************/
  
 /******* EXTRACT INPUTS AND PARAMETERS ************/
  
 

  /******* ALLOCATE MEMORY *************/
 

 /******* DO ANALYSIS ************/

 phase = params->lg + 2.0*LAL_PI * (tgps - params->T0) / params->Pb;

 U = phase + params->e*sin(phase);
 du = phase - (U-params->e*sin(U));
 
 while (du > 1e-15)
 {
   du = phase - (U-params->e*sin(U));
   U += du;
 }
 
 if(params->model==1)
   output->deltaT = params->x*((cos(U) - params->e)*sin(params->w) + sqrt((1 + params->e*params->e))*sin(U)*cos(params->w));
 else
   output->deltaT = 0.0;
   
  /****** CLEAN UP *******/

  DETATCHSTATUSPTR(status); 
  RETURN(status);
}

/*****************************************************************************/

/********************************** <lalVerbatim file="SlopeDetectorFilterCV">
Author: Daw, E. J.
$Id$  
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{SlopeDetectorFilter.c}}

[A one-line description of the function(s) defined in this module.]

\subsubsection*{Prototypes}
\input{SlopeDetectorFilterCP}
\idx{LALSlopeDetectorFilter()}

\subsubsection*{Description}

[A description of the data analysis task performed by this function; 
this is the main place to documeny the module]

\subsubsection*{Algorithm}

[A description of the method used to perform the calculation.]

\subsubsection*{Uses}

% List any external functions called by this function.
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}

[Any relevant notes]

\vfill{\footnotesize\input{SlopeDetectorFilterCV}}

******************************************************* </lalLaTeX> */ 

/******* INCLUDE STANDARD LIBRARY HEADERS; ************/
/* note LALStdLib.h already includes stdio.h and stdarg.h */

/******* INCLUDE ANY LDAS LIBRARY HEADERS ************/

/******* INCLUDE ANY LAL HEADERS ************/
#include <lal/LALStdlib.h>
#include <lal/SlopeDetectorFilter.h>

/******* DEFINE RCS ID STRING ************/
NRCSID( SLOPEDETECTORFILTERC, "$Id$" );

/******* DEFINE LOCAL CONSTANTS AND MACROS ************/

/******* DECLARE LOCAL (static) FUNCTIONS ************/
/* (definitions can go here or at the end of the file) */

/* static REAL4 findMean(); */

/******* DEFINE GLOBAL FUNCTIONS ************/

/* <lalVerbatim file="SlopeDetectorFilterCP"> */
void
LALSlopeDetectorFilter( LALStatus          *status,
			REAL4Vector*       output_data,
			const REAL4Vector* input_data,
			const UINT4        ntaps )
/* </lalVerbatim> */
{
  /******* DECLARE VARIABLES; for example: ************/

  REAL4              meanxt,meant,meanx,meantsquared;
  UINT4              i,j,datalength;
  REAL4*             datain;
  REAL4*             dataout;

  INITSTATUS( status, "LALSlopeDetectorFilter", SLOPEDETECTORFILTERC );
  /* only needed if you call other LAL functions within your code */
  /* ATTATCHSTATUSPTR (status); */

  if ( input_data == NULL ) {
    ABORT( status, SLOPEDETECTORFILTERH_EINPUTNULLP, 
	   SLOPEDETECTORFILTERH_MSGEINPUTNULLP );
  }
  if ( output_data == NULL ) {
    ABORT( status, SLOPEDETECTORFILTERH_EOUTPUTNULLP, 
	   SLOPEDETECTORFILTERH_MSGEOUTPUTNULLP );
  }

  datalength = input_data->length;
  datain = input_data->data;
  dataout = output_data->data;

  /******* CHECK VALIDITY OF ARGUMENTS; for example ************/

  if ( dataout == NULL ) {
    ABORT( status, SLOPEDETECTORFILTERH_EINPUTNULLP, 
	   SLOPEDETECTORFILTERH_MSGEINPUTNULLP );
  }
  if ( datain == NULL ) {
    ABORT( status, SLOPEDETECTORFILTERH_EOUTPUTNULLP, 
	   SLOPEDETECTORFILTERH_MSGEOUTPUTNULLP );
  }
  if ( datalength <= ntaps ) {
    ABORT( status, SLOPEDETECTORFILTERH_EDATATOOSHORT, 
	   SLOPEDETECTORFILTERH_MSGEDATATOOSHORT );
  }

  /* add checks that the input series is valid data, etc */

  /******* DO ANALYSIS ************/

  for(i=0;i<(datalength-ntaps+1);++i) {
    meanx=0; meant=0; meanxt=0; meantsquared=0;
    for(j=0;j<ntaps;++j) {
      meanx += datain[i+j]/ntaps;
      meant += j/ntaps;
      meanxt += j*datain[i+j]/ntaps;
      meantsquared += j*j/ntaps;
    }
    dataout[i] = (meanxt - meanx * meant) /
      ( meantsquared - meant * meant );
  }
  
  /* DETATCHSTATUSPTR(status); */
  RETURN(status);
}


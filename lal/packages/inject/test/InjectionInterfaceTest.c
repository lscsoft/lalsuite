/***************************************** <lalVerbatim file="InterfaceTest">
Author: Cokelaer, T.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Program \texttt{InterfaceTest.c}}
\label{ss:InterfaceTest.c}

Interface to generate any kind of gravitational waves signal.

\subsubsection*{Usage}
\begin{verbatim}
\texttt{InterfaceTest}
\end{verbatim}

\subsubsection*{Description}
Right now this test file read an input xml (injection.xml) file and for 
each valid line it computes the approriate waveform. Those waveforms are 
either produce within the inject package (PPN waveform) or within the 
inspiral package (TaylorT[1,2,3], EOB, SpinTaylor, PadeT1). 

Then the function LALGeneralInspiral create the amplitude, freq and phi 
vectors needed for further injection including the data itself (noise, h(t)) 

Finally, the injections are stored in a vector which is saved in "injection.dat"
file.

\subsubsection*{Uses}
\begin{verbatim}




\end{verbatim}
\subsubsection*{Notes}
\vfill{\footnotesize\input{InterfaceTestCV}}

</lalLaTeX><lalErrTable> */
#define INTERFACETESTC_ENORM 	0
#define INTERFACETESTC_ESUB  	1
#define INTERFACETESTC_EARG  	2
#define INTERFACETESTC_EVAL  	3
#define INTERFACETESTC_EFILE 	4
#define INTERFACETESTC_EMEM  	5
#define INTERFACETESTC_EINJECT  6

#define INTERFACETESTC_MSGENORM 	"Normal exit"
#define INTERFACETESTC_MSGESUB  	"Subroutine failed"
#define INTERFACETESTC_MSGEARG  	"Error parsing arguments"
#define INTERFACETESTC_MSGEVAL  	"Input argument out of valid range"
#define INTERFACETESTC_MSGEFILE 	"Could not open file"
#define INTERFACETESTC_MSGEMEM  	"Out of memory"
#define INTERFACETESTC_MSGEINJECT  	"No valid injection to do ... ? "
/* </lalErrTable><lalLaTeX>*/

/* --- the names of the files to be used --- */
#define INTERFACETEST_INJECTIONXMLFILE    "injection.xml"
#define INTERFACETEST_INJECTIONOUTPUTFILE "injection.dat"

/* --- include files --- */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include <lal/GenerateInspiral.h>
#include <lal/VectorOps.h>
#include <lal/SeqFactories.h>

#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>

#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>

#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/FindChirp.h>
#include <lal/PrintFTSeries.h>

NRCSID( INTERFACETESTC, "$Id$" );

#define ERROR( code, msg, statement )                                \
do                                                                   \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
		 "        %s %s\n", (code), program, __FILE__,       \
		 __LINE__, INTERFACETESTC, statement ? statement :   \
                 "", (msg) );                                        \
}                                                                    \
while (0)

#define SUB( func, statusptr )                                       \
do                                                                   \
if ( (func), (statusptr)->statusCode )                               \
{                                                                    \
  ERROR( INTERFACETESTC_ESUB, INTERFACETESTC_MSGESUB,                \
         "Function call \"" #func "\" failed:" );                    \
  exit( INTERFACETESTC_ESUB );                                       \
}                                                                    \
while (0)

char *program;



/* --- Main program start here --- */
int main(int argc, char **argv)
{
  /* some variables to create the data */
  UINT4 	startTime;
  UINT4		endTime;
  UINT4 	k;
  UINT4 	numPoints 	= 524288  ; /* arbitrary length of the data*/
  UINT4         numInjections   = 0; /* by default no injection. */
  REAL8 	sampling	= 2048.;  
  static  LALStatus 	status;
  /* injection structure */
  SimInspiralTable    *injections = NULL;
  SimInspiralTable    *thisInj    = NULL;
  /* the data */
  REAL4TimeSeries               ts ; /* A time series to store the injection */
  COMPLEX8FrequencySeries       fs ; /* A freq series to store the psd 	*/
  FILE		*output;             /* output result*/
  
   /* --- for debugging --- */
  lalDebugLevel = 0;
  program       = *argv;

  /* --- Start Main part here --- */
  /* First, we test if the injection xml file exist */
  if ((output	= fopen(INTERFACETEST_INJECTIONXMLFILE,"r")) == NULL){
	ERROR(INTERFACETESTC_EFILE,INTERFACETESTC_MSGEFILE, 0);
	exit(0);
  }
  else {
	fclose(output); /*if it exist let's close it */
  }

  /* then let's start to open the output file */
  if ((output 	= fopen(INTERFACETEST_INJECTIONOUTPUTFILE,"w")) == NULL){
      ERROR(INTERFACETESTC_EFILE,INTERFACETESTC_MSGEFILE, 0);
      exit(0);
    }
    
  /* let's allocate memory for to create some dummy data to 
     inject the waveforms */
  memset(&ts, 0, sizeof(REAL4TimeSeries));		
  LALSCreateVector( &status, &(ts.data), numPoints);	
  
  memset( &fs, 0, sizeof(COMPLEX8FrequencySeries));	
  LALCCreateVector( &status, &(fs.data), numPoints / 2 + 1 );

  
  /* --- Let's fix some variables we have to be in agreement with the xml file data --- */
  ts.epoch.gpsSeconds 	= 729273610;		       	/* gps time of the time series		*/
  startTime 		= 729273610;		       	/* start time and end time of ..	*/	
  endTime   		= startTime + 200;	       	/* ..injection; should be in agreement..*/
  						       	/* ..with the xml file			*/ 
  ts.sampleUnits 	= lalADCCountUnit;	       	/*  UNITY ?? 				*/
  ts.deltaT 		= 1./sampling;		       	/* sampling				*/
  ts.name[0]='X'; 					/* L, H, G, V or T for the detector 
							   otherwise it is optimally oriented	*/ 


  fs.deltaF 		= sampling / numPoints;	       /* idem for fs				*/
  fs.sampleUnits 	= lalADCCountUnit;

  /* --- the psd is flat for simplicity --- */		
  for( k = 0 ; k< numPoints/2; k++){
    fs.data->data[k].re = 1;
    fs.data->data[k].im = 0.;
  }
  /* --- and time series is zero --- */
  for( k = 0 ; k< numPoints; k++){
    ts.data->data[k] = 0.;
  }
   
  /* --- read injection  here --- */
  SUB(numInjections = SimInspiralTableFromLIGOLw( &injections, 
						  INTERFACETEST_INJECTIONXMLFILE,
						  startTime,
						  endTime), &status);

  /* any injection to do ? */
  if ( numInjections <= 0 )
    {
      ERROR(INTERFACETESTC_EINJECT, INTERFACETESTC_MSGEINJECT, 0);
      exit( 1 );
    }

  /* --- inject here --- */
  SUB( LALFindChirpInjectSignals( &status,
				  &ts, 
				  injections,
				  &fs ), &status);
  

  /* --- now we save the results --- */
  for (k = 0; k < numPoints; k++){
    fprintf(output,"%15.12f %e\n",  (float)k/sampling, ts.data->data[k]);
  }
  fclose(output);
  
  SUB( LALSDestroyVector( &status, &(ts.data) ), &status );
  SUB( LALCDestroyVector( &status, &(fs.data) ), &status );
  
  /* --- and finally free memory --- */
  while ( injections )
    {
      thisInj 		= injections;
      injections 	= injections->next;
      LALFree(thisInj);
    }

  LALCheckMemoryLeaks();
  return 0;
}

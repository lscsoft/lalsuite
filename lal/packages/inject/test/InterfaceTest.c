/***************************************** <lalVerbatim file="InterfaceTest">
Author: Cokelaer, T.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>


DRAFT DOCUMENTATION
 right now this test file read an input file and for each valid line compute
 the approriate waveform. The waveform are eiter produce inside the inject pacakge 
either inside the inspiral package. 

This file is composed of two main part. One to read the file and parse the parameter and one
to generate the waveform.

Parsing the parameters has required to create a new struct GeneralInspiralStruc which can parse either
"inspiral" data either "inject" data depending of the input data. 

Then the function LALGeneralInspiral create the amplitude, freq and phi vectors needed for further
injection including the data itself (noise, h(t)) 



\subsection{Program \texttt{InterfaceTest.c}}
\label{ss:InterfaceTest.c}

Interface to generate any kind of gravitational waves signal.

\subsubsection*{Usage}
\begin{verbatim}
InterfaceTest [-p parameterfile] 
\end{verbatim}

\subsubsection*{Description}


\subsubsection*{Exit codes}

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
liste des fonctions utilisees
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{InterfaceTestCV}}

******************************************************* </lalLaTeX> */
#define INTERFACETESTC_ENORM 0
#define INTERFACETESTC_ESUB  1
#define INTERFACETESTC_EARG  2
#define INTERFACETESTC_EVAL  3
#define INTERFACETESTC_EFILE 4
#define INTERFACETESTC_EMEM  5

#define INTERFACETESTC_MSGENORM "Normal exit"
#define INTERFACETESTC_MSGESUB  "Subroutine failed"
#define INTERFACETESTC_MSGEARG  "Error parsing arguments"
#define INTERFACETESTC_MSGEVAL  "Input argument out of valid range"
#define INTERFACETESTC_MSGEFILE "Could not open file"
#define INTERFACETESTC_MSGEMEM  "Out of memory"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <FrameL.h>



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
		 __LINE__, INTERFACETESTC, statement ? statement :      \
                 "", (msg) );                                        \
}                                                                    \
while (0)

#define SUB( func, statusptr )                                       \
do                                                                   \
if ( (func), (statusptr)->statusCode )                               \
{                                                                    \
  ERROR( INTERFACETESTC_ESUB, INTERFACETESTC_MSGESUB,                      \
         "Function call \"" #func "\" failed:" );                    \
  exit( INTERFACETESTC_ESUB );                                          \
}                                                                    \
while (0)

char *program;



INT4 lalDebugLevel=1;



int main(int argc, char **argv)
{
  /* input xml data */
  CHAR  	*injectionFile	= NULL;         
  
  /* some variables to create the data */
  UINT4 	startTime;
  UINT4		endTime;
  UINT4 	k;
  UINT4 	numPoints 	= 524288 * 2 / 16;
  UINT4         numInjections   = 0;
  REAL8 	sampling	= 2048.;
  
 static  LALStatus 	status 		;
  
  /* injection structure */
 SimInspiralTable    *injections = NULL;
 SimInspiralTable    *thisInj    = NULL;
 
  /* the data */
  REAL4TimeSeries               ts ;						/* A time series to store the injection */
  COMPLEX8FrequencySeries       fs ;						/* A freq series to store the psd 	*/

  /* output data to check injection */
  FILE		*output;
  PPNParamStruc         ppnParams;
  CoherentGW            waveform;

  
   /* --- MAIN --- */
  lalDebugLevel = 1;
  program = *argv;

  injectionFile = "injection.xml"; 						/* an xml file with the injection	*/
  output 	= fopen("injection.dat","w");					/* an ascii file for results		*/
  
  memset(&ts, 0, sizeof(REAL4TimeSeries));					/* allocate memory 			*/  
  LALSCreateVector( &status, &(ts.data), numPoints);	/* and ts null				*/
  
  memset( &fs, 0, sizeof(COMPLEX8FrequencySeries));				/* idem for fs				*/
  LALCCreateVector( &status, &(fs.data), numPoints / 2 + 1 );
  

  ts.epoch.gpsSeconds 	= 729273610;						/* gps time of the time series		*/
  startTime 		= 729273610;						/* gps start time and end time of ..	*/	
  endTime   		= startTime + 100;					/* ..injection; should be in agreement..*/
  										/* ..with the xml file			*/ 
  ts.sampleUnits 	= lalADCCountUnit;					/*  UNITY ?? 				*/
  ts.deltaT 		= 1./sampling;						/* sampling				*/
  fs.deltaF 		= sampling / numPoints;					/* idem for fs				*/
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
						  injectionFile,
						  startTime,
						  endTime), &status);


  if ( numInjections < 0 )
    {
      fprintf( stderr, "error: cannot read injection file (injection.xml)" );
      exit( 1 );
    }
  /* else we perform the injection */
  /* --- inject here --- */
  SUB( LALFindChirpInjectSignals( &status,
				  &ts, 
				  injections,
				  &fs ), &status);
  
  /*memset( &waveform, 0, sizeof(CoherentGW) );

  
  ppnParams.mTot = 10;
  ppnParams.eta  = .25;
  ppnParams.d    = 1 * 1.0e6 * LAL_PC_SI;
  ppnParams.inc  = 0;
  ppnParams.phi  = 1.;
  
  ppnParams.fStartIn = 40.0;
  ppnParams.fStopIn  = -1.0 / 
    (6.0 * sqrt(6.0) * LAL_PI * ppnParams.mTot * LAL_MTSUN_SI);
  
  ppnParams.position.longitude   = 0;
  ppnParams.position.latitude    = 0;
  ppnParams.position.system      = COORDINATESYSTEM_EQUATORIAL;
  ppnParams.psi                  = 0;
  ppnParams.epoch.gpsSeconds     = 0;
  ppnParams.epoch.gpsNanoSeconds = 0;
  
  SUB( LALGenerateInspiral(&status, &waveform, injections, &ppnParams), &status);
  */

  for (k=0; k<numPoints; k++){
    fprintf(output,"%15.12lf %e\n",  (float)k/sampling, ts.data->data[k]);
  }
  fclose(output);
  
  SUB( LALSDestroyVector( &status, &(ts.data) ), &status );
  SUB( LALCDestroyVector( &status, &(fs.data) ), &status );
  
  /*
    LALSDestroyVectorSequence( &status, &(waveform.a->data) );
    LALSDestroyVector(&status, &(waveform.f->data));
    LALDDestroyVector(&status, &(waveform.phi->data));
    
    LALFree( waveform.a );
    LALFree( waveform.f );
    LALFree( waveform.phi );
  */

  while ( injections )
    {
      thisInj = injections;
      injections = injections->next;
      LALFree( thisInj );
    }



  LALCheckMemoryLeaks();

  return 0;
}


#if 0
    case SpinOrbitCW: 
      params->socw.epoch.gpsSeconds = params->socw.epoch.gpsNanoSeconds = 0;
      params->socw.spinEpoch.gpsSeconds = params->socw.spinEpoch.gpsNanoSeconds = 0;
      params->socw.orbitEpoch.gpsSeconds = params->socw.orbitEpoch.gpsNanoSeconds = 0;
      params->socw.position.latitude = dec* LAL_PI/180.0;
      params->socw.position.longitude = ra* LAL_PI/180.0;
      params->socw.position.system = COORDINATESYSTEM_EQUATORIAL;
      params->socw.psi = psi* LAL_PI/180.0;
      if (deltaT!=0) 
	params->socw.deltaT = deltaT;
      else
	params->socw.deltaT = 1./ sampling;
      params->socw.length = length;
      params->socw.aPlus=a1;
      params->socw.aCross=a2;
      params->socw.phi0=phi0 * LAL_PI/180.0;
      params->socw.f0=f0 ;
      params->socw.omega = arg * LAL_PI/180.0;
      params->socw.rPeriNorm = rp;   /* if =0 -> equivaut a Taylor*/
      params->socw.oneMinusEcc = 1 - e;
      params->socw.angularSpeed = udot ;  
      params->socw.f = NULL;
      
      if ( nspin ) {
	LALDCreateVector( status, &(params->socw.f), nspin );
	ClearStatus(status);
      }
      for (i=0; i<nspin; i++)	
	params->socw.f->data[i] = fdata[i];	  
	
      
      /* vp = rp * udot entre ]0, 1[
	 argument  = omega*/
      
      break;
    case TaylorCW:
      params->taylorcw.epoch.gpsSeconds = params->taylorcw.epoch.gpsNanoSeconds = 0;
      params->taylorcw.epoch.gpsSeconds = params->taylorcw.epoch.gpsNanoSeconds = 0;
      params->taylorcw.position.latitude = dec* LAL_PI/180.0;
      params->taylorcw.position.longitude = ra* LAL_PI/180.0;
      params->taylorcw.position.system = COORDINATESYSTEM_EQUATORIAL;
      params->taylorcw.psi = psi* LAL_PI/180.0;
      if (deltaT!=0) 
	params->taylorcw.deltaT = deltaT;
      else
	params->taylorcw.deltaT = 1./ sampling;
      params->taylorcw.length = length;
      params->taylorcw.aPlus=a1;
      params->taylorcw.aCross=a2;
      params->taylorcw.phi0=phi0 * LAL_PI/180.0;
      params->taylorcw.f0=f0 ;   /* ?? pouirquoi ce terme pi */
          
      params->taylorcw.f = NULL;
      if ( nspin  ) {
	LALDCreateVector( status,&(params->taylorcw.f), nspin );
	ClearStatus(status);
      }
      for (i=0; i<params->taylorcw.f->length; i++)
	params->taylorcw.f->data[i] = fdata[i];	  
      break;
    }
  }


#endif

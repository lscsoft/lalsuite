/*
*  Copyright (C) 2007 David McKechan
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

/********************** <lalVerbatim file="COMPLETEGeneratePPNAmpTruncInspiralTestCV">
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Program \texttt{GeneratePPNAmpTruncInspiralTest.c}}
\label{ss:GeneratePPNAmpTruncInspiralTest.c}

Generates a parametrized post-Newtonian inspiral waveform
using the truncated version.

\subsubsection*{Usage}
\begin{verbatim}
GeneratePPNAmpTruncInspiralTest [-m m1 m2] [-r dist] [-i inc phii] [-f fmin fmax]
                        [-t dt] [-w deltat] [-p order] [-d debuglevel] [-o outfile] [-g fftoutfile]
\end{verbatim}

****************************************** </lalLaTeX><lalErrTable> */
#define GENERATEPPNINSPIRALTESTC_ENORM  0
#define GENERATEPPNINSPIRALTESTC_ESUB   1
#define GENERATEPPNINSPIRALTESTC_EARG   2
#define GENERATEPPNINSPIRALTESTC_EVAL   3
#define GENERATEPPNINSPIRALTESTC_EFILE  4
#define GENERATEPPNINSPIRALTESTC_EPRINT 5
#define DEBUG 1
#define DEBUG2 0
#define BUFFSIZE 1024     /* Number of timesteps buffered */

#define GENERATEPPNINSPIRALTESTC_MSGENORM  "Normal exit"
#define GENERATEPPNINSPIRALTESTC_MSGESUB   "Subroutine failed"
#define GENERATEPPNINSPIRALTESTC_MSGEARG   "Error parsing arguments"
#define GENERATEPPNINSPIRALTESTC_MSGEVAL   "Input argument out of valid range"
#define GENERATEPPNINSPIRALTESTC_MSGEFILE  "Could not open file"
#define GENERATEPPNINSPIRALTESTC_MSGEPRINT "Wrote past end of message string"

#include <math.h>
#include <stdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GeneratePPNInspiral.h>

/* headers from SimulateCoherentGW.c */ 
#include <lal/LALError.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/TimeDelay.h>
#include <lal/LALBarycenter.h>
#include <lal/VectorOps.h>
#include <lal/SkyCoordinates.h>

/* headers from TimeFreqFFTTest.c */
#include <stdio.h>
#include <unistd.h>
#include <lal/Random.h>
#include <lal/PrintFTSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/LALMoment.h>

#include <lal/LALInspiral.h>

NRCSID( GENERATEPPNINSPIRALTESTC, "$Id$" );

/* Default parameter settings. */
int lalDebugLevel = 1;
#define EPOCH (315187200000000000LL) /* about Jan. 1, 1990 */
#define M1    (1.4)
#define M2    (1.4)
#define DIST  (100000)
#define INC   (90.0)
#define PHI   (0.0)
#define FMIN  (40.0)
#define FMAX  (1000.0)
#define DT    (0.0005)
#define ORDER (7)

/* Usage format string. */
#define USAGE "Usage: %s [-g fftoutfile] [-m m1 m2] [-r dist] [-i inc phii]\n\t[-f fmin fmax] [-t dt] [-w deltat] [-p order] [-d debuglevel] [-o outfile]\n"

/* Maximum output message length. */
#define MSGLENGTH (1024)


/* Macros for printing errors and testing subroutines. */
#define ERROR( code, msg, statement )                                \
do                                                                   \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
		 "        %s %s\n", (code), *argv, __FILE__,         \
		 __LINE__, GENERATEPPNINSPIRALTESTC,                 \
		 statement ? statement : "", (msg) );                \
}                                                                    \
while (0)

#define INFO( statement )                                            \
do                                                                   \
if ( lalDebugLevel & LALINFO )                                       \
{                                                                    \
  LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"       \
		 "        %s\n", *argv, __FILE__, __LINE__,          \
		 GENERATEPPNINSPIRALTESTC, (statement) );            \
}                                                                    \
while (0)

#define WARNING( statement )                                         \
do                                                                   \
if ( lalDebugLevel & LALWARNING )                                    \
{                                                                    \
  LALPrintError( "Warning[0]: program %s, file %s, line %d, %s\n"    \
		 "        %s\n", *argv, __FILE__, __LINE__,          \
		 GENERATEPPNINSPIRALTESTC, (statement) );            \
}                                                                    \
while (0)

#define SUB( func, statusptr )                                       \
do                                                                   \
if ( (func), (statusptr)->statusCode )                               \
{                                                                    \
  ERROR( GENERATEPPNINSPIRALTESTC_ESUB,                              \
	 GENERATEPPNINSPIRALTESTC_MSGESUB,                           \
	 "Function call \"" #func "\" failed:" );                    \
  return GENERATEPPNINSPIRALTESTC_ESUB;                              \
}                                                                    \
while (0)

#define CHECKVAL( val, lower, upper )                                \
do                                                                   \
if ( ( (val) < (lower) ) || ( (val) > (upper) ) )                    \
{                                                                    \
  ERROR( GENERATEPPNINSPIRALTESTC_EVAL,                              \
	 GENERATEPPNINSPIRALTESTC_MSGEVAL,                           \
         "Value of " #val " out of range:" );                        \
  LALPrintError( #val " = %f, range = [%f,%f]\n", (REAL8)(val),      \
                 (REAL8)(lower), (REAL8)(upper) );                   \
  return GENERATEPPNINSPIRALTESTC_EVAL;                              \
}                                                                    \
while (0)

/* EXPANSION ********************************************************************************************** */

/* Definition of a data buffer list for storing the waveform. */
typedef struct tagPPNInspiralBuffer {
  REAL4 h[2*BUFFSIZE];               /* polarisation data */
  REAL8 phi[BUFFSIZE];               /* phase data */
  REAL4 f[BUFFSIZE];                 /* frequency data */
  struct tagPPNInspiralBuffer *next; /* next buffer in list */
} PPNInspiralBuffer;

/* Definition of a macro to free the tail of said list, from a given
   node onward. */
#define FREELIST( node )                                             \
do {                                                                 \
  PPNInspiralBuffer *herePtr = (node);                               \
  while ( herePtr ) {                                                \
    PPNInspiralBuffer *lastPtr = herePtr;                            \
    herePtr = herePtr->next;                                         \
    LALFree( lastPtr );                                              \
  }                                                                  \
} while (0)

/* END EXPANSION ****************************************************************************************** */


/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif

/* A function to convert INT8 nanoseconds to LIGOTimeGPS. */
void
I8ToLIGOTimeGPS( LIGOTimeGPS *output, INT8 input );


int
main(int argc, char **argv)
{
  /* Command-line parsing variables. */
  int arg;                      /* command-line argument counter */
  static LALStatus stat;        /* status structure */
  CHAR *outfile = NULL;         /* name of outfile */
  CHAR *fftout  = NULL; 	/* EXPANSION outfile */
  REAL4 m1 = M1, m2 = M2;       /* binary masses */
  REAL4 dist = DIST;            /* binary distance */
  REAL4 inc = 0.0, phii = 0.0;  /* inclination and coalescence phase */
  REAL4 fmin = FMIN, fmax=FMAX; /* start and stop frequencies */
  REAL8 dt = DT;                /* sampling interval */
  REAL8 deltat = 0.0;           /* wave sampling interval */
  INT4 order = ORDER;           /* PN order */
  UINT4 wlength = 0;
  UINT4 flength = 0;

  /* Other variables. */
  UINT4 i;                      /* index */
  CHAR message[MSGLENGTH];      /* signal generation output message */
  PPNParamStruc params;         /* input parameters */
  CoherentGW waveform;          /* output waveform */
  AmpSwitchStruc q;             /* Amplitude switches */
  FILE *fp;                     /* output file pointer */
  REAL4 *hoft;
  static REAL4Vector *htaper1, *htaper2; /* For LALInspiralWaveTaper */ 
  LALDetAMResponseSeries    am_response_series = {NULL,NULL,NULL};
  REAL4TimeSeries           plus_series, cross_series, scalar_series;
  LALTimeIntervalAndNSample time_info;
  LALSource                 pulsar;
  LALDetector               detector;
  LALDetAndSource           det_and_pulsar;
  FILE *fourier;
  static REAL4TimeSeries         ht;
  static COMPLEX8FrequencySeries Hf;
  RealFFTPlan    *fwdRealPlan    = NULL;
  RealFFTPlan    *revRealPlan    = NULL;
  REAL8 f = 0.0; 



  /*******************************************************************
   * ARGUMENT PARSING (arg stores the current position)              *
   *******************************************************************/

  arg = 1;
  while ( arg < argc ) {
    /* Parse mass option. */
    if ( !strcmp( argv[arg], "-m" ) ) {
      if ( argc > arg + 2 ) {
	arg++;
	m1 = atof( argv[arg++] );
	m2 = atof( argv[arg++] );
      }else{
	ERROR( GENERATEPPNINSPIRALTESTC_EARG,
	       GENERATEPPNINSPIRALTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return GENERATEPPNINSPIRALTESTC_EARG;
      }
    }
    /* Parse distance option. */
    else if ( !strcmp( argv[arg], "-r" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	dist = atof( argv[arg++] );
      }else{
	ERROR( GENERATEPPNINSPIRALTESTC_EARG,
	       GENERATEPPNINSPIRALTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return GENERATEPPNINSPIRALTESTC_EARG;
      }
    }
    /* Parse angles option. */
    else if ( !strcmp( argv[arg], "-i" ) ) {
      if ( argc > arg + 2 ) {
	arg++;
	inc = atof( argv[arg++] )*LAL_PI/180.0;
	phii = atof( argv[arg++] )*LAL_PI/180.0;
      }else{
	ERROR( GENERATEPPNINSPIRALTESTC_EARG,
	       GENERATEPPNINSPIRALTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return GENERATEPPNINSPIRALTESTC_EARG;
      }
    }
    /* Parse frequency option. */
    else if ( !strcmp( argv[arg], "-f" ) ) {
      if ( argc > arg + 2 ) {
	arg++;
	fmin = atof( argv[arg++] );
	fmax = atof( argv[arg++] );
      }else{
	ERROR( GENERATEPPNINSPIRALTESTC_EARG,
	       GENERATEPPNINSPIRALTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return GENERATEPPNINSPIRALTESTC_EARG;
      }
    }
    /* Parse sampling time option. */
    else if ( !strcmp( argv[arg], "-t" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	dt = atof( argv[arg++] );
      }else{
	ERROR( GENERATEPPNINSPIRALTESTC_EARG,
	       GENERATEPPNINSPIRALTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return GENERATEPPNINSPIRALTESTC_EARG;
      }
    }
    /* Parse waveform sampling time option. */
    else if ( !strcmp( argv[arg], "-w" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	deltat = atof( argv[arg++] );
      }else{
	ERROR( GENERATEPPNINSPIRALTESTC_EARG,
	       GENERATEPPNINSPIRALTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return GENERATEPPNINSPIRALTESTC_EARG;
      }
    }
    /* Parse PN order option. */
    else if ( !strcmp( argv[arg], "-p" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	order = atoi( argv[arg++] );
      }else{
	ERROR( GENERATEPPNINSPIRALTESTC_EARG,
	       GENERATEPPNINSPIRALTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return GENERATEPPNINSPIRALTESTC_EARG;
      }
    }
    /* Parse output file option. */
    else if ( !strcmp( argv[arg], "-o" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	outfile = argv[arg++];
      }else{
	ERROR( GENERATEPPNINSPIRALTESTC_EARG,
	       GENERATEPPNINSPIRALTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return GENERATEPPNINSPIRALTESTC_EARG;
      }
    }
    /* EXPANSION Parse output file option. */
    else if ( !strcmp( argv[arg], "-g" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	fftout = argv[arg++];
      }else{
	ERROR( GENERATEPPNINSPIRALTESTC_EARG,
	       GENERATEPPNINSPIRALTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return GENERATEPPNINSPIRALTESTC_EARG;
      }
    }
    /* Parse debug level option. */
    else if ( !strcmp( argv[arg], "-d" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	lalDebugLevel = atoi( argv[arg++] );
      }else{
	ERROR( GENERATEPPNINSPIRALTESTC_EARG,
	       GENERATEPPNINSPIRALTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return GENERATEPPNINSPIRALTESTC_EARG;
      }
    }
    /* Check for unrecognized options. */
    else if ( argv[arg][0] == '-' ) {
      ERROR( GENERATEPPNINSPIRALTESTC_EARG,
	     GENERATEPPNINSPIRALTESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return GENERATEPPNINSPIRALTESTC_EARG;
    }
  } /* End of argument parsing loop. */

  /* Make sure that values won't crash the system or anything. */
  CHECKVAL( order, -1, 8 );
  CHECKVAL( dt, LAL_REAL4_MIN, LAL_REAL4_MAX );
  CHECKVAL( deltat, 0.0, LAL_REAL4_MAX );


  /*******************************************************************
   * INPUT SETUP                                                     *
   *******************************************************************/

  /* Fixed parameters. */
  params.position.latitude = params.position.longitude = 0.0;
  params.position.system = COORDINATESYSTEM_EQUATORIAL;
  params.psi = 0.0;
  params.lengthIn = 0;

  /* Variable parameters. */
  I8ToLIGOTimeGPS( &(params.epoch), EPOCH );
  params.deltaT = dt;
  params.mTot = m1 + m2;
  params.eta = m1*m2/( params.mTot*params.mTot );
  params.inc = inc;
  params.phi = 0.0;
  params.d = dist*LAL_PC_SI*1.0e3;
  params.fStartIn = fmin;
  params.fStopIn = fmax;

  /* PPN parameter. */
  params.ppn = NULL;
  SUB( LALSCreateVector( &stat, &(params.ppn), order + 1 ), &stat );
  params.ppn->data[0] = 1.0;
  if ( order > 0 )
    params.ppn->data[1] = 0.0;
  for ( i = 2; i <= (UINT4)( order ); i++ )
    params.ppn->data[i] = 1.0;

  /* Amplitude switches */
  q.q0 = 1;
  q.q1 = 1;
  q.q2 = 1;
  q.q3 = 1;
  q.q4 = 1;
  q.q5 = 1;


  /* Output parameters. */
  memset( &waveform, 0, sizeof(CoherentGW) );

  /*******************************************************************
   * OUTPUT GENERATION                                               *
   *******************************************************************/

  /* Generate waveform. */
  SUB( LALGeneratePPNAmpTruncInspiral( &stat, &waveform, &params ), &stat );
#if DEBUG
  fprintf(stderr,"\n  Left GeneratePPNAmpTruncInspiral  \n\n");
#endif

  /**********************************************
   *                                            *
   *   EXPANSION  to produce h(t), H(f)         *
   *                                            *
   ********************************************************************************************************
   * This Test file now calculates the polar response functions for the detector and sky position defined *
   * below. It also performs the fourier transform to produce H(f).                                       *
   *                                                                                                      *
   * The output now has hPlus, hCross instead of aPlus and aCross. It also has h(t), ReH(f) and ImH(f)    *
   ********************************************************************************************************/
 
  wlength = waveform.h->data->length; 	
  flength = waveform.f->data->length;

#if DEBUG  
  fprintf(stderr," fFinal = %e\n", waveform.f->data->data[flength -1]);
#endif
  /* ************************************************** */
  /* Before we do anything let's taper hplus and hcross */
  /* This is a very inefficient interface for LALInspiralWaveTaper */
  LALCreateVector(&stat, &htaper1, wlength);
  LALCreateVector(&stat, &htaper2, wlength);

  for(i = 0; i < wlength; i++){
    htaper1->data[i] = waveform.h->data->data[2*i];
    htaper2->data[i] = waveform.h->data->data[2*i+1];
  }

  LALInspiralWaveTaper(&stat, htaper1, 3);
  LALInspiralWaveTaper(&stat, htaper2, 3);

  for(i = 0; i < wlength; i++){
    waveform.h->data->data[2*i] = htaper1->data[i];
    waveform.h->data->data[2*i+1] = htaper2->data[i];
  }

  /*************************** h(t)*/
 
  hoft = malloc(wlength*sizeof(REAL4));
 
  /* fake detector */
  detector.location[0] = 0.;
  detector.location[1] = 0.;
  detector.location[2] = LAL_AWGS84_SI;
  detector.response[0][0] = 0.;
  detector.response[1][1] = 0.5;
  detector.response[2][2] = -0.5;
  detector.response[0][1] = detector.response[1][0] = 0.;
  detector.response[0][2] = detector.response[2][0] = 0.;
  detector.response[1][2] = detector.response[2][1] = 0.;
  detector.type = LALDETECTORTYPE_ABSENT;
  
  pulsar.equatorialCoords.longitude = 0.; /* RA */
  pulsar.equatorialCoords.latitude  = 0.; /* Dec */
  pulsar.equatorialCoords.system    = COORDINATESYSTEM_EQUATORIAL;
  pulsar.orientation                = LAL_PI_2; /* orientation */


  det_and_pulsar.pDetector = &detector;
  det_and_pulsar.pSource   = &pulsar;

  plus_series.data = NULL;
  cross_series.data = NULL;
  scalar_series.data = NULL;
  
  am_response_series.pPlus   = &(plus_series);
  am_response_series.pCross  = &(cross_series);
  am_response_series.pScalar = &(scalar_series);

  LALSCreateVector(&stat, &(am_response_series.pPlus->data), 1);
  LALSCreateVector(&stat, &(am_response_series.pCross->data), 1);
  LALSCreateVector(&stat, &(am_response_series.pScalar->data), 1);
  
  time_info.epoch.gpsSeconds     = 61094;
  time_info.epoch.gpsNanoSeconds = 640000000;
  time_info.deltaT               = dt;
  time_info.nSample              = wlength;
  time_info.accuracy             = LALLEAPSEC_STRICT;

  LALComputeDetAMResponseSeries(&stat,
                                &am_response_series,
                                &det_and_pulsar,
                                &time_info);
  
#if DEBUG2
  printf("\n  Done computing AM response vectors\n");
  printf("  am_response_series.pPlus->data->length = %d\n",
          am_response_series.pPlus->data->length);
  printf("  am_response_series.pCross->data->length = %d\n",
          am_response_series.pCross->data->length);
  printf("  am_response_series.pScalar->data->length = %d\n",
          am_response_series.pScalar->data->length);
#endif

#if DEBUG2
  printf("\n Check data length = %d\n\n", wlength);
  

  printf("  TimeSeries data written to files plus_series.txt, ");
  printf("  cross_series.txt, and scalar_series.txt\n");

  LALSPrintTimeSeries(am_response_series.pPlus, "  plus_series.txt");
  LALSPrintTimeSeries(am_response_series.pCross, "  cross_series.txt");
  LALSPrintTimeSeries(am_response_series.pScalar, "  scalar_series.txt");
#endif

  for ( i = 0; i < wlength; i++){
    hoft[i] = waveform.h->data->data[2*i]*am_response_series.pPlus->data->data[i] +
              waveform.h->data->data[2*i+1]*am_response_series.pCross->data->data[i];
#if DEBUG2
    if(i <5){
      printf("\n\n  hplus = %e   hcross = %e    pplus = %e    pcross = %e \n", 
              waveform.h->data->data[2*i],
              waveform.h->data->data[2*i+1],
              am_response_series.pPlus->data->data[i],
              am_response_series.pCross->data->data[i]);		  
      printf("  hoft %e", hoft[i]);
    } 
#endif
  }	

  /*********************** End h(t)*/
  
  /*************************** H(F)*/


  LALSCreateVector( &stat, &ht.data, wlength );
  LALCCreateVector( &stat, &Hf.data, wlength / 2 + 1 );
  
  LALCreateForwardRealFFTPlan( &stat, &fwdRealPlan, wlength, 0 );
  LALCreateReverseRealFFTPlan( &stat, &revRealPlan, wlength, 0 );
  
  ht.f0 = 0;
  ht.deltaT = dt;
  for( i = 0; i < wlength ; i++)
    ht.data->data[i] = hoft[i];
    
  LALTimeFreqRealFFT( &stat, &Hf, &ht, fwdRealPlan );

#if DEBUG
  printf("\n\n h(t)length = %d\n H(F)length = %d\n ", wlength, Hf.data->length);
  printf("\n  Writing FFT data to fourier file...\n\n");
#endif  
 
  /* Write output. */
  if ( ( fourier = fopen( fftout, "w" ) ) == NULL ) 
    fourier = fopen("fftout", "w");
  
  for(i = 0; i < wlength/2+1; i++, f+=Hf.deltaF) 
    fprintf(fourier," %f %10.3e %10.3e\n", f, Hf.data->data[i].re, Hf.data->data[i].im);	  
  fclose(fourier);

  LALFreqTimeRealFFT( &stat, &ht, &Hf, revRealPlan );

  /*********************** End H(f)*/

#if DEBUG
  printf("\n\n Expansion code finished \n\n");
#endif

  /******************************************
   * END EXPANSION                          *
   ******************************************/

  /* Print termination information. */
  LALSnprintf( message, MSGLENGTH, "%d: %s", params.termCode,
	       params.termDescription );
  INFO( message );

  /* Print coalescence phase.
  LALSnprintf( message, MSGLENGTH,
	       "Waveform ends %.3f cycles before coalescence",
	       -waveform.phi->data->data[waveform.phi->data->length-1]
	       / (REAL4)( LAL_TWOPI ) ); */
  {
    INT4 code = sprintf( message,
			 "Waveform ends %.3f cycles before coalescence",
			 -waveform.phi->data->data[waveform.phi->data->length
						  -1]
			 / (REAL4)( LAL_TWOPI ) );
    if ( code >= MSGLENGTH || code < 0 ) {
      ERROR( GENERATEPPNINSPIRALTESTC_EPRINT,
	     GENERATEPPNINSPIRALTESTC_MSGEPRINT, 0 );
      return GENERATEPPNINSPIRALTESTC_EPRINT;
    }
  }
  INFO( message );

  /* Check if sampling interval was too large. */
  if ( params.dfdt > 2.0 ) {
    LALSnprintf( message, MSGLENGTH,
		 "Waveform sampling interval is too large:\n"
		 "\tmaximum df*dt = %f", params.dfdt );
    WARNING( message );
  }

  /* Shift phase. 
  phii -= waveform.phi->data->data[0];
  for ( i = 0; i < waveform.phi->data->length; i++ )
    waveform.phi->data->data[i] += phii;
    */

  /* Write output. */
  if ( outfile ) {
    if ( ( fp = fopen( outfile, "w" ) ) == NULL ) {
      ERROR( GENERATEPPNINSPIRALTESTC_EFILE,
	     GENERATEPPNINSPIRALTESTC_MSGEFILE, outfile );
      return GENERATEPPNINSPIRALTESTC_EFILE;
    }

    /* t phi f h+ hx ht Hfre Hfim ht? */
    if ( deltat == 0.0 ) {
      REAL8 t = 0.0; /* time */
      for ( i = 0; i < waveform.h->data->length; i++, t += dt )
	fprintf( fp, "%f %.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n", t,
	         waveform.phi->data->data[i],
		 waveform.f->data->data[i],
		 waveform.h->data->data[2*i],
		 waveform.h->data->data[2*i+1],
		 hoft[i],
		 ht.data->data[i]);
    }

    /* Waveform: */
    else {
    /*REAL8 t = 0.0;
      REAL8 x = 0.0;
      REAL8 dx = deltat/dt;
      REAL8 xMax = waveform.a->data->length - 1;
      REAL8 *phiData = waveform.phi->data->data;
      REAL4 *fData = waveform.f->data->data;
      REAL4 *hData = waveform.h->data->data;
      for ( ; x < xMax; x += dx, t += deltat ) {
	UINT4 j = floor( x );
	REAL8 frac = x - j;
	REAL8 p = frac*phiData[j+1] + ( 1.0 - frac )*phiData[j];
	REAL8 f = frac*fData[j+1] + ( 1.0 - frac )*fData[j];
	REAL8 ap = frac*aData[2*j+2] + ( 1.0 - frac )*aData[2*j];
	REAL8 ac = frac*aData[2*j+3] + ( 1.0 - frac )*aData[2*j+1];

	fprintf( fp, "%f %.3f %10.3e %10.3e %10.3e\n", t, p, f,
		 ap*cos( p ), ac*sin( p ) );
	fflush( fp );
      }*/
    }

    fclose( fp );
  }

  /*******************************************************************
   * CLEANUP                                                         *
   *******************************************************************/

  SUB( LALSDestroyVector( &stat, &(params.ppn) ), &stat );
  SUB( LALSDestroyVectorSequence( &stat, &(waveform.h->data) ),
       &stat );
  SUB( LALSDestroyVector( &stat, &(waveform.f->data) ), &stat );
  SUB( LALDDestroyVector( &stat, &(waveform.phi->data) ), &stat );
  LALFree( waveform.h );
  LALFree( waveform.f );
  LALFree( waveform.phi );

  /* Housekeeping of the extension */
  free(hoft);
  LALDestroyVector(&stat, &htaper1);
  LALDestroyVector(&stat, &htaper2);
  LALSDestroyVector(&stat, &(am_response_series.pPlus->data));
  LALSDestroyVector(&stat, &(am_response_series.pCross->data));
  LALSDestroyVector(&stat, &(am_response_series.pScalar->data));
  LALDestroyRealFFTPlan( &stat, &fwdRealPlan );
  LALDestroyRealFFTPlan( &stat, &revRealPlan );
  LALCDestroyVector( &stat, &Hf.data );
  LALSDestroyVector( &stat, &ht.data );

  LALCheckMemoryLeaks();
  INFO( GENERATEPPNINSPIRALTESTC_MSGENORM );
  return GENERATEPPNINSPIRALTESTC_ENORM;
}


/* A function to convert INT8 nanoseconds to LIGOTimeGPS. */
void
I8ToLIGOTimeGPS( LIGOTimeGPS *output, INT8 input )
{
  INT8 s = input / 1000000000LL;
  output->gpsSeconds = (INT4)( s );
  output->gpsNanoSeconds = (INT4)( input - 1000000000LL*s );
  return;
}

/*
*  Copyright (C) 2007 David McKechan, Thomas Cokelaer
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

/**
\file

\brief Generates a parametrized post-Newtonian inspiral waveform.
\ingroup GeneratePPNInspiral_h

\heading{Usage}
\code
GeneratePPNAmpCorInspiralTest [-m m1 m2] [-r dist] [-i inc phii psi] [-f f_min f_max]
                        [-t dt] [-p order amp] [-d debuglevel] [-o outfile] [-g FF FFTfile]
                        [-s taper]
\endcode
*/

/** \name Error Codes */
/*@{*/
#define GENERATEPPNINSPIRALTESTC_ENORM  0	/**< Normal exit */
#define GENERATEPPNINSPIRALTESTC_ESUB   1	/**< Subroutine failed */
#define GENERATEPPNINSPIRALTESTC_EARG   2	/**< Error parsing arguments */
#define GENERATEPPNINSPIRALTESTC_EVAL   3	/**< Input argument out of valid range */
#define GENERATEPPNINSPIRALTESTC_EFILE  4	/**< Could not open file */
#define GENERATEPPNINSPIRALTESTC_EPRINT 5	/**< Wrote past end of message string */
/*@}*/

/** \cond DONT_DOXYGEN */
#define GENERATEPPNINSPIRALTESTC_MSGENORM  "Normal exit"
#define GENERATEPPNINSPIRALTESTC_MSGESUB   "Subroutine failed"
#define GENERATEPPNINSPIRALTESTC_MSGEARG   "Error parsing arguments"
#define GENERATEPPNINSPIRALTESTC_MSGEVAL   "Input argument out of valid range"
#define GENERATEPPNINSPIRALTESTC_MSGEFILE  "Could not open file"
#define GENERATEPPNINSPIRALTESTC_MSGEPRINT "Wrote past end of message string"

#define BUFFSIZE 1024     /* Number of timesteps buffered */

#define LAL_USE_OLD_COMPLEX_STRUCTS
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

/* Default parameter settings. */
extern int lalDebugLevel;
#define EPOCH (315187200000000000LL) /* about Jan. 1, 1990 */
#define M1    (1.4)
#define M2    (1.4)
#define DIST  (100000)
#define INC   (90.0)
#define PHI   (0.0)
#define FMIN  (40.0)
#define FMAX  (0.0) /* This means Flso */
#define DT    (0.00048828125) /* 1/2048 */
#define ORDER (7)
#define AMP (5)

/* Usage format string. */
#define USAGE "Usage: %s [-g FFToutfile] [-m m1 m2] [-r dist] [-i inc phii psi]\n\t[-f f_min f_max] [-t dt] [-p order amp] [-d debuglevel] [-o outfile]\n\t [-s taper]"

/* Maximum output message length. */
#define MSGLENGTH (1024)


/* Macros for printing errors and testing subroutines. */
#define ERROR( code, msg, statement )                                \
do                                                                   \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
		 "        %s %s\n", (code), *argv, __FILE__,         \
		 __LINE__, "$Id$",                 \
		 statement ? statement : "", (msg) );                \
}                                                                    \
while (0)

#define INFO( statement )                                            \
do                                                                   \
if ( lalDebugLevel & LALINFO )                                       \
{                                                                    \
  LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"       \
		 "        %s\n", *argv, __FILE__, __LINE__,          \
		 "$Id$", (statement) );            \
}                                                                    \
while (0)

#define WARNING( statement )                                         \
do                                                                   \
if ( lalDebugLevel & LALWARNING )                                    \
{                                                                    \
  LALPrintError( "Warning[0]: program %s, file %s, line %d, %s\n"    \
		 "        %s\n", *argv, __FILE__, __LINE__,          \
		 "$Id$", (statement) );            \
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
  CHAR *fftout  = NULL; 	      /* FFT outfile */
  REAL4 m1 = M1, m2 = M2;       /* binary masses */
  REAL4 dist = DIST;            /* binary distance */
  REAL4 inc = 0.0, psi = LAL_PI_2;  /* inclination, and polarization angle */
  REAL4 f_min = FMIN, f_max=FMAX; /* start and stop frequencies */
  REAL8 dt = DT;                /* sampling interval */
  INT4 order = ORDER;           /* PN order */
  INT4 amp = AMP;               /* Amplitude switches */
  UINT4 taper = 0;		          /* Taper switch (On > 0) */

  /* Other variables. */
  UINT4 i;                      /* index */
  CHAR message[MSGLENGTH];      /* signal generation output message */
  PPNParamStruc params;         /* input parameters */
  CoherentGW waveform;          /* output waveform */
  FILE *fp;                     /* output file pointer */
  static REAL4Vector	      *hoft;
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
  REAL8 t = 0.0; /* time */
  REAL8 f = 0.0;

  lalDebugLevel = 1;

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
      if ( argc > arg + 3 ) {
	arg++;
	inc = atof( argv[arg++] )*LAL_PI/180.0;
        psi = atof(argv[arg++] )*LAL_PI/180.0;
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
	f_min = atof( argv[arg++] );
	f_max = atof( argv[arg++] );
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
    /* Parse PN order option. */
    else if ( !strcmp( argv[arg], "-p" ) ) {
      if ( argc > arg + 2 ) {
	arg++;
	order = atoi( argv[arg++] );
        amp = atoi( argv[arg++] );
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
    /* Parse FFToutput file option. */
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
    /* Parse tapering option. */
    else if ( !strcmp( argv[arg], "-s" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	taper = atoi( argv[arg++] );
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
  CHECKVAL( amp, 0, 5);
  CHECKVAL( dt, LAL_REAL4_MIN, LAL_REAL4_MAX );

  /*******************************************************************
   * INPUT SETUP                                                     *
   *******************************************************************/

  /* Fixed parameters. */
  params.position.latitude = params.position.longitude = 0.0;
  params.position.system = COORDINATESYSTEM_EQUATORIAL;
  params.lengthIn = 0;

  /* Variable parameters. */
  I8ToLIGOTimeGPS( &(params.epoch), EPOCH );
  params.deltaT = dt;
  params.mTot = m1 + m2;
  params.eta = m1*m2/( params.mTot*params.mTot );
  params.inc = inc;
  params.phi = 0.0;
  params.psi = psi;
  params.d = dist*LAL_PC_SI*1.0e3;
  params.fStartIn = f_min;
  params.fStopIn = f_max;

  /* Amplitude switches */
  params.ampOrder = amp;

  /* PPN parameter. */
  params.ppn = NULL;
  SUB( LALSCreateVector( &stat, &(params.ppn), order + 1 ), &stat );
  params.ppn->data[0] = 1.0;
  if ( order > 0 )
    params.ppn->data[1] = 0.0;
  for ( i = 2; i <= (UINT4)( order ); i++ )
    params.ppn->data[i] = 1.0;

  /* Output parameters. */
  memset( &waveform, 0, sizeof(CoherentGW) );

  /*******************************************************************
   * OUTPUT GENERATION                                               *
   *******************************************************************/

  /* Generate waveform. */
  SUB( LALGeneratePPNAmpCorInspiral( &stat, &waveform, &params ), &stat );
   /*                                            *
   ********************************************************************************************************
   * This Test file now calculates the polar response functions for the detector and sky position defined *
   * below. It also performs the fourier transform to produce H(f) if an FFToutfie is specified.          *
   *                                                                                                      *
   ********************************************************************************************************/

  /*************************** h(t)*/
  LALCreateVector( &stat, &hoft, waveform.h->data->length);

  /* fake detector */
  /* This one is overhead */
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
  pulsar.orientation                = params.psi; /* orientation */


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
  time_info.nSample              = waveform.h->data->length;

  LALComputeDetAMResponseSeries(&stat,
                                &am_response_series,
                                &det_and_pulsar,
                                &time_info);

  for ( i = 0; i < waveform.h->data->length; i++)
  {
    hoft->data[i] = waveform.h->data->data[2*i]*am_response_series.pPlus->data->data[i] +
                   waveform.h->data->data[2*i+1]*am_response_series.pCross->data->data[i];
  }

  /* Taper hoft */
  if( taper > 0 )
    XLALSimInspiralREAL4WaveTaper(hoft, LAL_SIM_INSPIRAL_TAPER_STARTEND);

  if( fftout )
  {
    LALSCreateVector( &stat, &ht.data, waveform.h->data->length );
    LALCCreateVector( &stat, &Hf.data, waveform.h->data->length / 2 + 1 );
    LALCreateForwardRealFFTPlan( &stat, &fwdRealPlan, waveform.h->data->length, 0 );

    ht.f0 = 0;
    ht.deltaT = dt;
    for( i = 0; i < waveform.h->data->length ; i++)
      ht.data->data[i] = hoft->data[i];

    LALTimeFreqRealFFT( &stat, &Hf, &ht, fwdRealPlan );

    if( ( fourier = fopen(fftout, "w")) == NULL)
      fourier = fopen("fftout", "w");

    for(i = 0; i < waveform.h->data->length/ 2 + 1; i++, f+=Hf.deltaF)
      fprintf(fourier," %f %1.6e %1.6e\n", f, crealf(Hf.data->data[i]), Hf.data->data[i].im);
    fclose(fourier);

		LALDestroyRealFFTPlan( &stat, &fwdRealPlan );
    LALCDestroyVector( &stat, &Hf.data );
    LALSDestroyVector( &stat, &ht.data );
  }

  /* Print termination information. */
  snprintf( message, MSGLENGTH, "%d: %s", params.termCode,
	       params.termDescription );
  INFO( message );

  /* Print coalescence phase.*/
  snprintf( message, MSGLENGTH,
	       "Waveform ends %.3f cycles before coalescence",
	       -waveform.phi->data->data[waveform.phi->data->length-1]
	       / (REAL4)( LAL_TWOPI ) );
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
    snprintf( message, MSGLENGTH,
		 "Waveform sampling interval is too large:\n"
		 "\tmaximum df*dt = %f", params.dfdt );
    WARNING( message );
  }

  /* Write output. */
  if ( outfile ) {
    if ( ( fp = fopen( outfile, "w" ) ) == NULL ) {
      ERROR( GENERATEPPNINSPIRALTESTC_EFILE,
	     GENERATEPPNINSPIRALTESTC_MSGEFILE, outfile );
      return GENERATEPPNINSPIRALTESTC_EFILE;
    }

    /* t phi f h+ hx ht  */
    for ( i = 0; i < waveform.h->data->length; i++, t += dt )
      fprintf( fp, "%f %.3e %1.6e %1.6e %1.6e %1.6e \n", t,
                  waveform.phi->data->data[i],
		  waveform.f->data->data[i],
		  waveform.h->data->data[2*i],
		  waveform.h->data->data[2*i+1],
		  hoft->data[i]);

    fclose( fp );
  }

  /*******************************************************************
   * CLEANUP                                                         *
   *******************************************************************/

  SUB( LALSDestroyVector( &stat, &(params.ppn) ), &stat );
  SUB( LALSDestroyVectorSequence( &stat, &(waveform.h->data) ),
       &stat );
  SUB( LALSDestroyVectorSequence( &stat, &(waveform.a->data) ),
       &stat );
  SUB( LALSDestroyVector( &stat, &(waveform.f->data) ), &stat );
  SUB( LALDDestroyVector( &stat, &(waveform.phi->data) ), &stat );
  LALFree( waveform.h );
  LALFree( waveform.a );
  LALFree( waveform.f );
  LALFree( waveform.phi );
  LALDestroyVector( &stat, &hoft );
  LALSDestroyVector( &stat, &(am_response_series.pPlus->data) );
  LALSDestroyVector( &stat, &(am_response_series.pCross->data) );
  LALSDestroyVector( &stat, &(am_response_series.pScalar->data) );


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
/** \endcond */

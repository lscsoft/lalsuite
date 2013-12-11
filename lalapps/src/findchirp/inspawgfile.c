/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Patrick Brady
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

/*
 * inspawgfile.c
 * Author: Brown, D. A, and Creighton, T. D.
 */

#include <math.h>
#include <stdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/SeqFactories.h>
#include <lal/Units.h>
#include <lal/Random.h>
#include <lal/VectorOps.h>
#include <lal/Inject.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/StreamInput.h>
#include <lal/TimeFreqFFT.h>

#define INSPAWGFILEC_ENORM  0
#define INSPAWGFILEC_ESUB   1
#define INSPAWGFILEC_EARG   2
#define INSPAWGFILEC_EVAL   3
#define INSPAWGFILEC_EFILE  4
#define INSPAWGFILEC_EINPUT 5
#define INSPAWGFILEC_EMEM   6

#define INSPAWGFILEC_MSGENORM  "Normal exit"
#define INSPAWGFILEC_MSGESUB   "Subroutine failed"
#define INSPAWGFILEC_MSGEARG   "Error parsing arguments"
#define INSPAWGFILEC_MSGEVAL   "Input argument out of valid range"
#define INSPAWGFILEC_MSGEFILE  "Could not open file"
#define INSPAWGFILEC_MSGEINPUT "Error reading file"
#define INSPAWGFILEC_MSGEMEM   "Out of memory"

/* Default parameter settings. */
#define EPOCH (0)
#define DIST  (0.00002*LAL_MRSUN_SI )
#define M1    (1.4)
#define M2    (1.4)
#define INC   (0.0)
#define PHIC  (0.0)
#define SEC   (0)
#define NSEC  (0)
#define NPT   (1048576)
#define DT    (1.0/16384.0)
#define SIGMA (0.0)

/* Global constants. */
#define MSGLEN (256) /* maximum length of warning/info messages */
#define FSTART (40.0) /* initial frequency of waveform */
#define FSTOP  (3000.0) /* termination frequency of waveform */
#define DELTAT (0.00006103515625) /* sampling interval of amplitude and phase */

/* Usage format string. */
#define USAGE "Usage: %s [-s sourcefile] [-r respfile] [-o outfile] [-e seed]\n" \
"                [-i infile | -n sec nsec npt dt sigma] [-fi lowf] [-fe highf] \n" \
"                [-d debuglevel] [-p]\n"

/* Macros for printing errors and testing subroutines. */
#define ERROR( code, msg, statement )                                \
do                                                                   \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
		 "        %s %s\n", (code), *argv, __FILE__,         \
		 __LINE__, "$Id$", statement ? statement : \
                 "", (msg) );                                        \
}                                                                    \
while (0)

#define INFO( statement )                                            \
do                                                                   \
if ( lalDebugLevel & LALINFO )                                       \
{                                                                    \
  LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"       \
		 "        %s\n", *argv, __FILE__, __LINE__,          \
		 "$Id$", (statement) );                    \
}                                                                    \
while (0)

#define WARNING( statement )                                         \
do                                                                   \
if ( lalDebugLevel & LALWARNING )                                    \
{                                                                    \
  LALPrintError( "Warning[0]: program %s, file %s, line %d, %s\n"    \
		 "        %s\n", *argv, __FILE__, __LINE__,          \
		 "$Id$", (statement) );                    \
}                                                                    \
while (0)

#define SUB( func, statusptr )                                       \
do                                                                   \
if ( (func), (statusptr)->statusCode )                               \
{                                                                    \
  ERROR( INSPAWGFILEC_ESUB, INSPAWGFILEC_MSGESUB,            \
         "Function call \"" #func "\" failed:" );                    \
  return INSPAWGFILEC_ESUB;                                      \
}                                                                    \
while (0)

#define CHECKVAL( val, lower, upper )                                \
do                                                                   \
if ( ( (val) <= (lower) ) || ( (val) > (upper) ) )                   \
{                                                                    \
  ERROR( INSPAWGFILEC_EVAL, INSPAWGFILEC_MSGEVAL,            \
         "Value of " #val " out of range:" );                        \
  LALPrintError( #val " = %f, range = (%f,%f]\n", (REAL8)(val),      \
                 (REAL8)(lower), (REAL8)(upper) );                   \
  return INSPAWGFILEC_EVAL;                                      \
}                                                                    \
while (0)

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
  int arg;                   /* command-line argument counter */
  static LALStatus stat;     /* status structure */
  CHAR *sourcefile = NULL;   /* name of sourcefile */
  CHAR *respfile = NULL;     /* name of respfile */
  CHAR *infile = NULL;       /* name of infile */
  CHAR *outfile = NULL;      /* name of outfile */
  CHAR *specfile = NULL;     /* name of spectrum file */
  INT4 seed = 0;             /* random number seed */
  INT4 sec = SEC;            /* ouput epoch.gpsSeconds */
  INT4 nsec = NSEC;          /* ouput epoch.gpsNanoSeconds */
  INT4 npt = NPT;            /* number of output points */
  REAL8 dt = DT;             /* output sampling interval */
  REAL4 sigma = SIGMA;       /* noise amplitude */
  REAL4 fstart = FSTART;     /* start frequency */
  REAL4 fstop  = FSTOP;      /* stop frequency */

  /* File reading variables. */
  FILE *fp = NULL; /* generic file pointer */
  BOOLEAN ok = 1;  /* whether input format is correct */
  UINT4 i;         /* generic index over file lines */
  INT8 epoch;      /* epoch stored as an INT8 */

  /* Other global variables. */
  RandomParams *params = NULL; /* parameters of pseudorandom sequence */
  DetectorResponse detector;   /* the detector in question */
  REAL4TimeSeries output;      /* detector ADC output */


  /*******************************************************************
   * PARSE ARGUMENTS (arg stores the current position)               *
   *******************************************************************/

  /* Exit gracefully if no arguments were given.
  if ( argc <= 1 ) {
    INFO( "No testing done." );
    return 0;
  } */

  arg = 1;
  while ( arg < argc ) {
    /* Parse source file option. */
    if ( !strcmp( argv[arg], "-s" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	sourcefile = argv[arg++];
      }else{
	ERROR( INSPAWGFILEC_EARG, INSPAWGFILEC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return INSPAWGFILEC_EARG;
      }
    }
    /* Parse response file option. */
    else if ( !strcmp( argv[arg], "-r" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	respfile = argv[arg++];
      }else{
	ERROR( INSPAWGFILEC_EARG, INSPAWGFILEC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return INSPAWGFILEC_EARG;
      }
    }
    /* Parse input file option. */
    else if ( !strcmp( argv[arg], "-i" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	infile = argv[arg++];
      }else{
	ERROR( INSPAWGFILEC_EARG, INSPAWGFILEC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return INSPAWGFILEC_EARG;
      }
    }
    /* Parse noise output option. */
    else if ( !strcmp( argv[arg], "-n" ) ) {
      if ( argc > arg + 5 ) {
	arg++;
	sec = atoi( argv[arg++] );
	nsec = atoi( argv[arg++] );
	npt = atoi( argv[arg++] );
	dt = atof( argv[arg++] );
	sigma = atof( argv[arg++] );
      }else{
	ERROR( INSPAWGFILEC_EARG, INSPAWGFILEC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return INSPAWGFILEC_EARG;
      }
    }
    /* Parse output file option. */
    else if ( !strcmp( argv[arg], "-o" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	outfile = argv[arg++];
      }else{
	ERROR( INSPAWGFILEC_EARG, INSPAWGFILEC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return INSPAWGFILEC_EARG;
      }
    }
    /* Parse debug level option. */
    else if ( !strcmp( argv[arg], "-d" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
      }else{
	ERROR( INSPAWGFILEC_EARG, INSPAWGFILEC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return INSPAWGFILEC_EARG;
      }
    }
    /* Parse start frequency */
    else if ( !strcmp( argv[arg], "-fi" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	fstart = atof( argv[arg++] );
      }else{
	ERROR( INSPAWGFILEC_EARG, INSPAWGFILEC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return INSPAWGFILEC_EARG;
      }
    }
    /* Parse stop frequency */
    else if ( !strcmp( argv[arg], "-fe" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	fstop = atof( argv[arg++] );
      }else{
	ERROR( INSPAWGFILEC_EARG, INSPAWGFILEC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return INSPAWGFILEC_EARG;
      }
    }
    /* Parse random seed option. */
    else if ( !strcmp( argv[arg], "-e" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	seed = atoi( argv[arg++] );
      }else{
	ERROR( INSPAWGFILEC_EARG, INSPAWGFILEC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return INSPAWGFILEC_EARG;
      }
    }
    /* Parse print spectrum */
    else if ( !strcmp( argv[arg], "-p" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	specfile = argv[arg++];
      }else{
	ERROR( INSPAWGFILEC_EARG, INSPAWGFILEC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return INSPAWGFILEC_EARG;
      }
    }
    /* Check for unrecognized options. */
    else if ( argv[arg][0] == '-' ) {
      ERROR( INSPAWGFILEC_EARG, INSPAWGFILEC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return INSPAWGFILEC_EARG;
    }
  } /* End of argument parsing loop. */

  /* Check for redundant or bad argument values. */
  CHECKVAL( npt, 0, 2147483647 );
  CHECKVAL( dt, 0, LAL_REAL4_MAX );


  /*******************************************************************
   * SETUP                                                           *
   *******************************************************************/

  /* Set up output, detector, and random parameter structures. */
  output.data = NULL;
  detector.transfer = (COMPLEX8FrequencySeries *)
    LALMalloc( sizeof(COMPLEX8FrequencySeries) );
  if ( !(detector.transfer) ) {
    ERROR( INSPAWGFILEC_EMEM, INSPAWGFILEC_MSGEMEM, 0 );
    return INSPAWGFILEC_EMEM;
  }
  detector.transfer->data = NULL;
  detector.site = NULL;
  SUB( LALCreateRandomParams( &stat, &params, seed ), &stat );

  /* Set up units. */
  {
    RAT4 negOne = { -1, 0 };
    LALUnit unit;
    LALUnitPair pair;
    output.sampleUnits = lalADCCountUnit;
    pair.unitOne = &lalADCCountUnit;
    pair.unitTwo = &lalStrainUnit;
    SUB( LALUnitRaise( &stat, &unit, pair.unitTwo,
		       &negOne ), &stat );
    pair.unitTwo = &unit;
    SUB( LALUnitMultiply( &stat, &(detector.transfer->sampleUnits),
			  &pair ), &stat );
  }

  /* Read response function. */
  if ( respfile ) {
    REAL4VectorSequence *resp = NULL; /* response as vector sequence */
    COMPLEX8Vector *response = NULL;  /* response as complex vector */
    COMPLEX8Vector *unity = NULL;     /* vector of complex 1's */

    if ( ( fp = fopen( respfile, "r" ) ) == NULL ) {
      ERROR( INSPAWGFILEC_EFILE, INSPAWGFILEC_MSGEFILE,
	     respfile );
      return INSPAWGFILEC_EFILE;
    }

    /* Read header. */
    ok &= ( fscanf( fp, "# epoch = %" LAL_INT8_FORMAT "\n", &epoch ) == 1 );
    I8ToLIGOTimeGPS( &( detector.transfer->epoch ), epoch );
    ok &= ( fscanf( fp, "# f0 = %lf\n", &( detector.transfer->f0 ) )
	    == 1 );
    ok &= ( fscanf( fp, "# deltaF = %lf\n",
		    &( detector.transfer->deltaF ) ) == 1 );
    if ( !ok ) {
      ERROR( INSPAWGFILEC_EINPUT, INSPAWGFILEC_MSGEINPUT,
	     respfile );
      return INSPAWGFILEC_EINPUT;
    }

    /* Read and convert body to a COMPLEX8Vector. */
    SUB( LALSReadVectorSequence( &stat, &resp, fp ), &stat );
    fclose( fp );
    if ( resp->vectorLength != 2 ) {
      ERROR( INSPAWGFILEC_EINPUT, INSPAWGFILEC_MSGEINPUT,
	     respfile );
      return INSPAWGFILEC_EINPUT;
    }
    SUB( LALCCreateVector( &stat, &response, resp->length ), &stat );
    memcpy( response->data, resp->data, 2*resp->length*sizeof(REAL4) );
    SUB( LALSDestroyVectorSequence( &stat, &resp ), &stat );

    /* Convert response function to a transfer function. */
    SUB( LALCCreateVector( &stat, &unity, response->length ), &stat );
    for ( i = 0; i < response->length; i++ ) {
      unity->data[i] = 1.0;
    }
    SUB( LALCCreateVector( &stat, &( detector.transfer->data ),
			   response->length ), &stat );
    SUB( LALCCVectorDivide( &stat, detector.transfer->data, unity,
			    response ), &stat );
    SUB( LALCDestroyVector( &stat, &response ), &stat );
    SUB( LALCDestroyVector( &stat, &unity ), &stat );
  }

  /* No response file, so generate a unit response. */
  else {
    I8ToLIGOTimeGPS( &( detector.transfer->epoch ), EPOCH );
    detector.transfer->f0 = 0.0;
    detector.transfer->deltaF = 1.5*fstop;
    SUB( LALCCreateVector( &stat, &( detector.transfer->data ), 2 ),
	 &stat );
    detector.transfer->data->data[0] = 1.0;
    detector.transfer->data->data[1] = 1.0;
  }


  /* Read input data. */
  if ( infile ) {
    REAL4VectorSequence *input = NULL; /* input as vector sequence */
    if ( ( fp = fopen( infile, "r" ) ) == NULL ) {
      ERROR( INSPAWGFILEC_EFILE, INSPAWGFILEC_MSGEFILE,
	     infile );
      return INSPAWGFILEC_EFILE;
    }

    /* Read header. */
    ok &= ( fscanf( fp, "# epoch = %" LAL_INT8_FORMAT "\n", &epoch ) == 1 );
    I8ToLIGOTimeGPS( &( output.epoch ), epoch );
    ok &= ( fscanf( fp, "# deltaT = %lf\n", &( output.deltaT ) )
	    == 1 );
    if ( !ok ) {
      ERROR( INSPAWGFILEC_EINPUT, INSPAWGFILEC_MSGEINPUT,
	     infile );
      return INSPAWGFILEC_EINPUT;
    }

    /* Read and convert body. */
    SUB( LALSReadVectorSequence( &stat, &input, fp ), &stat );
    fclose( fp );
    if ( input->vectorLength != 1 ) {
      ERROR( INSPAWGFILEC_EINPUT, INSPAWGFILEC_MSGEINPUT,
	     infile );
      return INSPAWGFILEC_EINPUT;
    }
    SUB( LALSCreateVector( &stat, &( output.data ), input->length ),
	 &stat );
    for ( i = 0; i < input->length; i++ )
      output.data->data[i] = input->data[i];
    SUB( LALSDestroyVectorSequence( &stat, &input ), &stat );
  }

  /* No input file, so generate one randomly. */
  else {
    output.epoch.gpsSeconds = sec;
    output.epoch.gpsNanoSeconds = nsec;
    output.deltaT = dt;
    SUB( LALSCreateVector( &stat, &( output.data ), npt ), &stat );
    if ( sigma == 0 ) {
      memset( output.data->data, 0, npt*sizeof(REAL4) );
    } else {
      SUB( LALNormalDeviates( &stat, output.data, params ), &stat );
    }
  }


  /*******************************************************************
   * INJECTION                                                       *
   *******************************************************************/

  /* Open sourcefile. */
  if ( sourcefile ) {
    if ( ( fp = fopen( sourcefile, "r" ) ) == NULL ) {
      ERROR( INSPAWGFILEC_EFILE, INSPAWGFILEC_MSGEFILE,
	     sourcefile );
      return INSPAWGFILEC_EFILE;
    }
  }

  /* For each line in the sourcefile... */
  while ( ok ) {
    PPNParamStruc ppnParams;       /* wave generation parameters */
    REAL4 m1, m2, dist, inc, phic; /* unconverted parameters */
    CoherentGW waveform;           /* amplitude and phase structure */
    REAL4TimeSeries signalvec;     /* GW signal */
    REAL8 time;                    /* length of GW signal */
    CHAR timeCode;                 /* code for signal time alignment */
    CHAR message[MSGLEN];          /* warning/info messages */

    /* Read and convert input line. */
    if ( sourcefile ) {
      ok &= ( fscanf( fp, "%c %" LAL_INT8_FORMAT " %f %f %f %f %f\n", &timeCode,
		      &epoch, &m1, &m2, &dist, &inc, &phic ) == 7 );
      fprintf(stderr, "%c %" LAL_INT8_FORMAT " %f %f %f %f %f\n", timeCode,
		      epoch, m1, m2, dist, inc, phic );  fflush(stderr);
      ppnParams.mTot = m1 + m2;
      ppnParams.eta = m1*m2/( ppnParams.mTot*ppnParams.mTot );
      ppnParams.d = dist*LAL_PC_SI*1000.0;
      ppnParams.inc = inc*LAL_PI/180.0;
      ppnParams.phi = phic*LAL_PI/180.0;
    } else {
      timeCode = 'i';
      ppnParams.mTot = M1 + M2;
      ppnParams.eta = M1*M2/( ppnParams.mTot*ppnParams.mTot );
      ppnParams.d = DIST;
      ppnParams.inc = INC;
      ppnParams.phi = PHIC;
      epoch = EPOCH;
    }

    if ( ok ) {
      /* Set up other parameter structures. */
      ppnParams.epoch.gpsSeconds = ppnParams.epoch.gpsNanoSeconds = 0;
      ppnParams.position.latitude = ppnParams.position.longitude = 0.0;
      ppnParams.position.system = COORDINATESYSTEM_EQUATORIAL;
      ppnParams.psi = 0.0;
      ppnParams.fStartIn = fstart;
      ppnParams.fStopIn = fstop;
      ppnParams.lengthIn = 0;
      ppnParams.ppn = NULL;
      ppnParams.deltaT = DELTAT;
      memset( &waveform, 0, sizeof(CoherentGW) );

      /* Generate waveform at zero epoch. */
      SUB( LALGeneratePPNInspiral( &stat, &waveform, &ppnParams ),
	   &stat );
      snprintf( message, MSGLEN, "%d: %s", ppnParams.termCode,
          ppnParams.termDescription );
      INFO( message );
      if ( lalDebugLevel & LALINFO )
      {
        LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"
            "\tWaveform terminared at %e Hz\n", 
            *argv, __FILE__, __LINE__,
            "$Id$", ppnParams.fStop );
      }
      if ( ppnParams.dfdt > 2.0 ) {
        snprintf( message, MSGLEN,
		     "Waveform sampling interval is too large:\n"
		     "\tmaximum df*dt = %f", ppnParams.dfdt );
	WARNING( message );
      }

      /* Compute epoch for waveform. */
      time = waveform.a->data->length*DELTAT;
      if ( timeCode == 'f' )
	epoch -= (INT8)( 1000000000.0*time );
      else if ( timeCode == 'c' )
	epoch -= (INT8)( 1000000000.0*ppnParams.tc );
      I8ToLIGOTimeGPS( &( waveform.a->epoch ), epoch );
      waveform.f->epoch = waveform.phi->epoch = waveform.a->epoch;

      /* Generate and inject signal. */
      signalvec.epoch = waveform.a->epoch;
      signalvec.epoch.gpsSeconds -= 1;
      signalvec.deltaT = output.deltaT/4.0;
      signalvec.f0 = 0;
      signalvec.data = NULL;
      time = ( time + 2.0 )/signalvec.deltaT;
      SUB( LALSCreateVector( &stat, &( signalvec.data ), (UINT4)time ),
	   &stat );
      SUB( LALSimulateCoherentGW( &stat, &signalvec, &waveform,
				  &detector ), &stat );
      SUB( LALSSInjectTimeSeries( &stat, &output, &signalvec ),
	   &stat );
      SUB( LALSDestroyVectorSequence( &stat, &( waveform.a->data ) ),
	   &stat );
      SUB( LALSDestroyVector( &stat, &( waveform.f->data ) ), &stat );
      SUB( LALDDestroyVector( &stat, &( waveform.phi->data ) ), &stat );
      LALFree( waveform.a );
      LALFree( waveform.f );
      LALFree( waveform.phi );
      SUB( LALSDestroyVector( &stat, &( signalvec.data ) ), &stat );
    }

    /* If there is no source file, inject only one source. */
    if ( !sourcefile )
      ok = 0;
  }

  /* Input file is exhausted (or has a badly-formatted line ). */
  if ( sourcefile )
    fclose( fp );

  /*******************************************************************
   * COMPUTE POWER SPECTRUM                                          *
   *******************************************************************/

  /* write diagnostic info to disk */
  if ( specfile ){
    fprintf(stderr, "error, unable to compute power spectrum.\n" );
    fprintf(stderr, "LALRealAverageSpectrum has been removed from LAL.\n" );
    fprintf(stderr, "Fix code to use new LALREAL4AverageSpectrum routine.\n" );
    fprintf(stderr, "Aborting.\n" );
    abort();
#if 0
      REAL4FrequencySeries     *fseries;
      RealDFTParams            *dftparams        = NULL;
      LALWindowParams           winParams;
      INT4 length=2048;

      /* Set up the window parameters */
      winParams.type=0;
      winParams.length=length;

      /* assign temporary memory for the frequency data */
      fseries = (REAL4FrequencySeries *) 
          LALMalloc (sizeof(REAL4FrequencySeries));
      strncpy( fseries->name, "anonymous", LALNameLength );
      fseries->data = NULL;
      SUB( LALCreateVector (&stat, &fseries->data, 
                  length/2 + 1), &stat ); 

      /* create the dft params */
      SUB( LALCreateRealDFTParams(&stat , &dftparams, &winParams, 1),
              &stat);

      /* compute the average spectrum */
      SUB( LALRealAverageSpectrum (&stat, fseries, &output, dftparams,
                  0), &stat );
      LALSPrintFrequencySeries ( fseries, specfile );

      SUB( LALDestroyRealDFTParams(&stat, &dftparams), &stat);
      SUB( LALDestroyVector(&stat, &fseries->data), &stat);
      LALFree( fseries );
#endif
  }


  /*******************************************************************
   * CLEANUP                                                         *
   *******************************************************************/

  /* Print output file. */
  if ( outfile ) {
    if ( ( fp = fopen( outfile, "w" ) ) == NULL ) {
      ERROR( INSPAWGFILEC_EFILE, INSPAWGFILEC_MSGEFILE,
	     outfile );
      return INSPAWGFILEC_EFILE;
    }
    for ( i = 0; i < output.data->length; i++ )
      fprintf( fp, "%e\n", output.data->data[i] );
    fclose( fp );
  }

  /* Destroy remaining memory. */
  SUB( LALDestroyRandomParams( &stat, &params ), &stat );
  SUB( LALSDestroyVector( &stat, &( output.data ) ), &stat );
  SUB( LALCDestroyVector( &stat, &( detector.transfer->data ) ),
       &stat );
  LALFree( detector.transfer );

  /* Done! */
  LALCheckMemoryLeaks();
  INFO( INSPAWGFILEC_MSGENORM );
  return INSPAWGFILEC_ENORM;
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

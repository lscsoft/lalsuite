/*
*  Copyright (C) 2007 Teviet Creighton
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
\author Creighton, T. D.
\file
\ingroup GeneratePPNInspiral_h

\brief Generates a parametrized post-Newtonian inspiral waveform.

\heading{Usage}
\code
GeneratePPNInspiralTest [-m m1 m2] [-r dist] [-i inc phii] [-f f_min f_max]
                        [-t dt] [-w deltat] [-p order] [-d debuglevel] [-o outfile]
\endcode

\heading{Description}

This program generates the amplitude, phase, and frequency of a
post-Newtonian inspiral waveform as functions of time.  The following
option flags are accepted:
<ul>
<li>[<tt>-m</tt>] Sets the binary masses to \c m1 and \c m2
solar massses (default values: \f$1.4M_\odot\f$).</li>
<li>[<tt>-r</tt>] Sets the binary system distance to \c dist kpc
(default value: 8.5kpc).</li>
<li>[<tt>-i</tt>] Sets the inclination and \e initial phase
angles to \c inc and \c phii degrees (default values:
0\ degrees).</li>
<li>[<tt>-f</tt>] Sets the initial and final wave frequencies to
\c f_min and \c f_max Hz (default values: 40Hz and 500Hz).</li>
<li>[<tt>-t</tt>] Sets the waveform sampling interval to \c dt
seconds (default value: 0.01s).</li>
<li>[<tt>-w</tt>] Generates actual waveforms rather than phase and
amplitude functions, sampled at intervals of \c deltat seconds (no
default).</li>
<li>[<tt>-p</tt>] Sets the post\f${}^{n/2}\f$-Newtonian order to
\f$n=\f$\c order (default value: \f$n=4\f$).</li>
<li>[<tt>-d</tt>] Sets the debug level to \c debuglevel (default
value:\ 0).</li>
<li>[<tt>-o</tt>] Sets the output filename to \c outfile (by
default no output is produced).</li>
</ul>

\heading{Algorithm}

This program simply parses the command line, sets the appropriate
fields of a \c PPNParamStruc, and passes it in to
<tt>LALGeneratePPNInspiral()</tt>.  No maximum waveform length is
specified; the function will allocate as much data as necessary.

If the <tt>-w</tt> \e and <tt>-o</tt> options are given, the
amplitude, phase, and frequency are generated as above, but are then
resampled at intervals \c deltat to generate actual wave output.

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

#include <math.h>
#include <stdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GeneratePPNInspiral.h>

/* Default parameter settings. */
#define EPOCH (315187200000000000LL) /* about Jan. 1, 1990 */
#define M1    (1.4)
#define M2    (1.4)
#define DIST  (8.5)
#define INC   (0.0)
#define PHI   (0.0)
#define FMIN  (40.0)
#define FMAX  (500.0)
#define DT    (0.01)
#define ORDER (4)

/* Usage format string. */
#define USAGE "Usage: %s [-m m1 m2] [-r dist] [-i inc phii]\n\t[-f f_min f_max] [-t dt] [-w deltat] [-p order] [-d debuglevel] [-o outfile]\n"

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
  REAL4 m1 = M1, m2 = M2;       /* binary masses */
  REAL4 dist = DIST;            /* binary distance */
  REAL4 inc = 0.0, phii = 0.0;  /* inclination and coalescence phase */
  REAL4 f_min = FMIN, f_max=FMAX; /* start and stop frequencies */
  REAL8 dt = DT;                /* sampling interval */
  REAL8 deltat = 0.0;           /* wave sampling interval */
  INT4 order = ORDER;           /* PN order */

  /* Other variables. */
  UINT4 i;                      /* index */
  CHAR message[MSGLENGTH];      /* signal generation output message */
  PPNParamStruc params;         /* input parameters */
  CoherentGW waveform;          /* output waveform */
  FILE *fp;                     /* output file pointer */


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
    /* Parse debug level option. */
    else if ( !strcmp( argv[arg], "-d" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
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
  CHECKVAL( order, -1, 5 );
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
  params.fStartIn = f_min;
  params.fStopIn = f_max;

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
  SUB( LALGeneratePPNInspiral( &stat, &waveform, &params ), &stat );

  /* Print termination information. */
  snprintf( message, MSGLENGTH, "%d: %s", params.termCode,
	       params.termDescription );
  INFO( message );

  /* Print coalescence phase.
  snprintf( message, MSGLENGTH,
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
    snprintf( message, MSGLENGTH,
		 "Waveform sampling interval is too large:\n"
		 "\tmaximum df*dt = %f", params.dfdt );
    WARNING( message );
  }

  /* Renormalize phase. */
  phii -= waveform.phi->data->data[0];
  for ( i = 0; i < waveform.phi->data->length; i++ )
    waveform.phi->data->data[i] += phii;

  /* Write output. */
  if ( outfile ) {
    if ( ( fp = fopen( outfile, "w" ) ) == NULL ) {
      ERROR( GENERATEPPNINSPIRALTESTC_EFILE,
	     GENERATEPPNINSPIRALTESTC_MSGEFILE, outfile );
      return GENERATEPPNINSPIRALTESTC_EFILE;
    }

    /* Amplitude and phase functions: */
    if ( deltat == 0.0 ) {
      REAL8 t = 0.0; /* time */
      for ( i = 0; i < waveform.a->data->length; i++, t += dt )
	fprintf( fp, "%f %.3f %10.3e %10.3e %10.3e\n", t,
		 waveform.phi->data->data[i],
		 waveform.f->data->data[i],
		 waveform.a->data->data[2*i],
		 waveform.a->data->data[2*i+1] );
    }

    /* Waveform: */
    else {
      REAL8 t = 0.0;
      REAL8 x = 0.0;
      REAL8 dx = deltat/dt;
      REAL8 xMax = waveform.a->data->length - 1;
      REAL8 *phiData = waveform.phi->data->data;
      REAL4 *fData = waveform.f->data->data;
      REAL4 *aData = waveform.a->data->data;
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
      }
    }

    fclose( fp );
  }


  /*******************************************************************
   * CLEANUP                                                         *
   *******************************************************************/

  SUB( LALSDestroyVector( &stat, &(params.ppn) ), &stat );
  SUB( LALSDestroyVectorSequence( &stat, &(waveform.a->data) ),
       &stat );
  SUB( LALSDestroyVector( &stat, &(waveform.f->data) ), &stat );
  SUB( LALDDestroyVector( &stat, &(waveform.phi->data) ), &stat );
  LALFree( waveform.a );
  LALFree( waveform.f );
  LALFree( waveform.phi );

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

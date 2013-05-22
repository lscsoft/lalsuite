/*
*  Copyright (C) 2007 Jolien Creighton, Teviet Creighton
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
   \ingroup StreamInput_h
   \author Creighton, T. D.

   \brief Reads a time or frequency series from a file, and writes it to another file.

\heading{Usage}
\code
StreamSeriesInputTest [-o outfile] [-i infile stype dtype] [-d debuglevel]
\endcode

\heading{Description}

This test program parses data from an input file or from \c stdin,
using the routines in \ref StreamSeriesInput.c, and possibly
generating output using the routines in \ref StreamSeriesOutput.c.
The following option flags are accepted:
<ul>
<li>[<tt>-o</tt>] Writes the output to \c outfile.  If
\c outfile is given as \c stdout, the data is written to
standard output (\e not to a file named \c stdout).  If the
<tt>-o</tt> flag is not given, the input routines are exercised, but no
output is written.</li>
<li>[<tt>-i</tt>] Specifies the input file name \c infile, series
type \c stype, and base datatype \c dtype.  Series type is a
single character: either \c t (time series), \c v (time vector
series), \c a (time array series), or \c f (frequency series).
Base datatype may be \c i2 (\c INT2), \c i4 (\c INT4),
\c i8 (\c INT8), \c u2 (\c UINT2), \c u4
(\c UINT4), \c u8 (\c UINT8), \c s (\c REAL4),
\c d (\c REAL8), \c c (\c COMPLEX8), or \c z
(\c COMPLEX16).  If the <tt>-i</tt> flag is not given,
<tt>-i StreamSeriesInput.dat f s</tt> is assumed (this file is provided
with the distribution so that running the code with no arguments, \'a
la <tt>make check</tt>, will perform a nontrivial test of the
algorithm).</li>
<li>[<tt>-d</tt>] Sets the debug level to \c debuglevel; if
absent, <tt>-d 0</tt> is assumed.</li>
</ul>

See the documentation in \ref StreamSeriesInput_c and
\ref StreamSeriesOutput_c for discussion of the input and output
data file formats.
*/

/** \name Error Codes */ /*@{*/
#define STREAMSERIESINPUTTESTC_ENORM 0  /**< Normal exit */
#define STREAMSERIESINPUTTESTC_ESUB  1  /**< Subroutine failed */
#define STREAMSERIESINPUTTESTC_EARG  2  /**< Error parsing arguments */
#define STREAMSERIESINPUTTESTC_EFILE 3  /**< Could not open file */
/*@}*/

/** \cond DONT_DOXYGEN */
#define STREAMSERIESINPUTTESTC_MSGENORM "Normal exit"
#define STREAMSERIESINPUTTESTC_MSGESUB  "Subroutine failed"
#define STREAMSERIESINPUTTESTC_MSGEARG  "Error parsing arguments"
#define STREAMSERIESINPUTTESTC_MSGEFILE "Could not open file"


#include <stdlib.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/StreamInput.h>
#include <lal/StreamOutput.h>

/* Default parameter settings. */
#define INFILE TEST_DATA_DIR "StreamSeriesInput.data"

/* Usage format string. */
#define USAGE "Usage: %s [-o outfile] [-d debuglevel]\n"        \
"\t[-i infile {t | v | a | f} {i2 | i4 | i8 | u2 | u4 | u8 | s | d | c | z}]\n"

/* Macros for printing errors and testing subroutines. */
#define ERROR( code, msg, statement )                                \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
		 "        %s %s\n", (code), *argv, __FILE__,         \
		 __LINE__, "$Id$",                   \
		 statement ? statement : "", (msg) );                \
}                                                                    \
else (void)(0)

#define INFO( statement )                                            \
if ( lalDebugLevel & LALINFO )                                       \
{                                                                    \
  LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"       \
		 "        %s\n", *argv, __FILE__, __LINE__,          \
		 "$Id$", (statement) );              \
}                                                                    \
else (void)(0)

#define SUB( func, statusptr )                                       \
if ( (func), (statusptr)->statusCode )                               \
{                                                                    \
  ERROR( STREAMSERIESINPUTTESTC_ESUB, STREAMSERIESINPUTTESTC_MSGESUB,\
         "Function call \"" #func "\" failed:" );                    \
  return STREAMSERIESINPUTTESTC_ESUB;                                \
}                                                                    \
else (void)(0)

/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif

int
main(int argc, char **argv)
{
  static LALStatus stat;       /* top-level status structure */
  INT4 arg;                    /* index over command-line options */
  CHAR *outfile = NULL;        /* name of output file */
  const CHAR *infile = INFILE; /* name of input file */
  const CHAR *dtype = "s";     /* data type tag */
  const CHAR *stype = "f";     /* series type tag */
  FILE *fpIn = NULL;           /* input file pointer */
  FILE *fpOut = NULL;          /* output file pointer */


  /* Parse argument list.  arg stores the current position. */
  arg = 1;
  while ( arg < argc ) {
    /* Parse output file option. */
    if ( !strcmp( argv[arg], "-o" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	outfile = argv[arg++];
      } else {
	ERROR( STREAMSERIESINPUTTESTC_EARG,
	       STREAMSERIESINPUTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return STREAMSERIESINPUTTESTC_EARG;
      }
    }
    /* Parse input file option. */
    else if ( !strcmp( argv[arg], "-i" ) ) {
      if ( argc > arg + 3 ) {
	arg++;
	infile = argv[arg++];
	stype = argv[arg++];
	dtype = argv[arg++];
      } else {
	ERROR( STREAMSERIESINPUTTESTC_EARG,
	       STREAMSERIESINPUTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return STREAMSERIESINPUTTESTC_EARG;
      }
    }
    /* Parse debug level option. */
    else if ( !strcmp( argv[arg], "-d" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
      } else {
	ERROR( STREAMSERIESINPUTTESTC_EARG,
	       STREAMSERIESINPUTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return STREAMSERIESINPUTTESTC_EARG;
      }
    }
    /* Check for unrecognized options. */
    else {
      ERROR( STREAMSERIESINPUTTESTC_EARG,
	     STREAMSERIESINPUTTESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return STREAMSERIESINPUTTESTC_EARG;
    }
  } /* End of argument parsing loop. */

  /* Open input and output files. */
  if ( strcmp( infile, "stdin" ) ) {
    if ( !( fpIn = LALOpenDataFile( infile ) ) ) {
      ERROR( STREAMSERIESINPUTTESTC_EFILE, "- "
	     STREAMSERIESINPUTTESTC_MSGEFILE, infile );
      return STREAMSERIESINPUTTESTC_EFILE;
    }
  } else
    fpIn = stdin;
  if ( outfile ) {
    if ( strcmp( outfile, "stdout" ) ) {
      if ( !( fpOut = fopen( outfile, "w" ) ) ) {
	ERROR( STREAMSERIESINPUTTESTC_EFILE, "- "
	       STREAMSERIESINPUTTESTC_MSGEFILE, outfile );
	return STREAMSERIESINPUTTESTC_EFILE;
      }
    } else
      fpOut = stdout;
  }

  /* Read (and possibly write) the data as a time series. */
  if ( !strcmp( stype, "t" ) ) {
    if ( !strcmp( dtype, "i2" ) ) {
      static INT2TimeSeries series;
      SUB( LALI2ReadTSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALI2WriteTSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALI2DestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "i4" ) ) {
      static INT4TimeSeries series;
      SUB( LALI4ReadTSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALI4WriteTSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALI4DestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "i8" ) ) {
      static INT8TimeSeries series;
      SUB( LALI8ReadTSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALI8WriteTSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALI8DestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "u2" ) ) {
      static UINT2TimeSeries series;
      SUB( LALU2ReadTSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALU2WriteTSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALU2DestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "u4" ) ) {
      static UINT4TimeSeries series;
      SUB( LALU4ReadTSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALU4WriteTSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALU4DestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "u8" ) ) {
      static UINT8TimeSeries series;
      SUB( LALU8ReadTSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALU8WriteTSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALU8DestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "s" ) ) {
      static REAL4TimeSeries series;
      SUB( LALSReadTSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALSWriteTSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALSDestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "d" ) ) {
      static REAL8TimeSeries series;
      SUB( LALDReadTSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALDWriteTSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALDDestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "c" ) ) {
      static COMPLEX8TimeSeries series;
      SUB( LALCReadTSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALCWriteTSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALCDestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "z" ) ) {
      static COMPLEX16TimeSeries series;
      SUB( LALZReadTSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALZWriteTSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALZDestroyVector( &stat, &(series.data) ), &stat );
    } else {
      ERROR( STREAMSERIESINPUTTESTC_EARG,
	     STREAMSERIESINPUTTESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return STREAMSERIESINPUTTESTC_EARG;
    }
  }

  /* Read (and possibly write) the data as a time vector series. */
  else if ( !strcmp( stype, "v" ) ) {
    if ( !strcmp( dtype, "i2" ) ) {
      static INT2TimeVectorSeries series;
      SUB( LALI2ReadTVectorSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALI2WriteTVectorSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALI2DestroyVectorSequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "i4" ) ) {
      static INT4TimeVectorSeries series;
      SUB( LALI4ReadTVectorSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALI4WriteTVectorSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALI4DestroyVectorSequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "i8" ) ) {
      static INT8TimeVectorSeries series;
      SUB( LALI8ReadTVectorSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALI8WriteTVectorSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALI8DestroyVectorSequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "u2" ) ) {
      static UINT2TimeVectorSeries series;
      SUB( LALU2ReadTVectorSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALU2WriteTVectorSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALU2DestroyVectorSequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "u4" ) ) {
      static UINT4TimeVectorSeries series;
      SUB( LALU4ReadTVectorSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALU4WriteTVectorSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALU4DestroyVectorSequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "u8" ) ) {
      static UINT8TimeVectorSeries series;
      SUB( LALU8ReadTVectorSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALU8WriteTVectorSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALU8DestroyVectorSequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "s" ) ) {
      static REAL4TimeVectorSeries series;
      SUB( LALSReadTVectorSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALSWriteTVectorSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALSDestroyVectorSequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "d" ) ) {
      static REAL8TimeVectorSeries series;
      SUB( LALDReadTVectorSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALDWriteTVectorSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALDDestroyVectorSequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "c" ) ) {
      static COMPLEX8TimeVectorSeries series;
      SUB( LALCReadTVectorSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALCWriteTVectorSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALCDestroyVectorSequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "z" ) ) {
      static COMPLEX16TimeVectorSeries series;
      SUB( LALZReadTVectorSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALZWriteTVectorSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALZDestroyVectorSequence( &stat, &(series.data) ), &stat );
    } else {
      ERROR( STREAMSERIESINPUTTESTC_EARG,
	     STREAMSERIESINPUTTESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return STREAMSERIESINPUTTESTC_EARG;
    }
  }

  /* Read (and possibly write) the data as a time array series. */
  else if ( !strcmp( stype, "a" ) ) {
    if ( !strcmp( dtype, "i2" ) ) {
      static INT2TimeArraySeries series;
      SUB( LALI2ReadTArraySeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALI2WriteTArraySeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALI2DestroyArraySequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "i4" ) ) {
      static INT4TimeArraySeries series;
      SUB( LALI4ReadTArraySeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALI4WriteTArraySeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALI4DestroyArraySequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "i8" ) ) {
      static INT8TimeArraySeries series;
      SUB( LALI8ReadTArraySeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALI8WriteTArraySeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALI8DestroyArraySequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "u2" ) ) {
      static UINT2TimeArraySeries series;
      SUB( LALU2ReadTArraySeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALU2WriteTArraySeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALU2DestroyArraySequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "u4" ) ) {
      static UINT4TimeArraySeries series;
      SUB( LALU4ReadTArraySeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALU4WriteTArraySeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALU4DestroyArraySequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "u8" ) ) {
      static UINT8TimeArraySeries series;
      SUB( LALU8ReadTArraySeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALU8WriteTArraySeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALU8DestroyArraySequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "s" ) ) {
      static REAL4TimeArraySeries series;
      SUB( LALSReadTArraySeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALSWriteTArraySeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALSDestroyArraySequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "d" ) ) {
      static REAL8TimeArraySeries series;
      SUB( LALDReadTArraySeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALDWriteTArraySeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALDDestroyArraySequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "c" ) ) {
      static COMPLEX8TimeArraySeries series;
      SUB( LALCReadTArraySeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALCWriteTArraySeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALCDestroyArraySequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "z" ) ) {
      static COMPLEX16TimeArraySeries series;
      SUB( LALZReadTArraySeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALZWriteTArraySeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALZDestroyArraySequence( &stat, &(series.data) ), &stat );
    } else {
      ERROR( STREAMSERIESINPUTTESTC_EARG,
	     STREAMSERIESINPUTTESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return STREAMSERIESINPUTTESTC_EARG;
    }
  }

  /* Read (and possibly write) the data as a frequency series. */
  else if ( !strcmp( stype, "f" ) ) {
    if ( !strcmp( dtype, "i2" ) ) {
      static INT2FrequencySeries series;
      SUB( LALI2ReadFSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALI2WriteFSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALI2DestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "i4" ) ) {
      static INT4FrequencySeries series;
      SUB( LALI4ReadFSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALI4WriteFSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALI4DestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "i8" ) ) {
      static INT8FrequencySeries series;
      SUB( LALI8ReadFSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALI8WriteFSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALI8DestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "u2" ) ) {
      static UINT2FrequencySeries series;
      SUB( LALU2ReadFSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALU2WriteFSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALU2DestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "u4" ) ) {
      static UINT4FrequencySeries series;
      SUB( LALU4ReadFSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALU4WriteFSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALU4DestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "u8" ) ) {
      static UINT8FrequencySeries series;
      SUB( LALU8ReadFSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALU8WriteFSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALU8DestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "s" ) ) {
      static REAL4FrequencySeries series;
      SUB( LALSReadFSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALSWriteFSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALSDestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "d" ) ) {
      static REAL8FrequencySeries series;
      SUB( LALDReadFSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALDWriteFSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALDDestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "c" ) ) {
      static COMPLEX8FrequencySeries series;
      SUB( LALCReadFSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALCWriteFSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALCDestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "z" ) ) {
      static COMPLEX16FrequencySeries series;
      SUB( LALZReadFSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALZWriteFSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALZDestroyVector( &stat, &(series.data) ), &stat );
    } else {
      ERROR( STREAMSERIESINPUTTESTC_EARG,
	     STREAMSERIESINPUTTESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return STREAMSERIESINPUTTESTC_EARG;
    }
  }

  /* Unrecognized series type. */
  else {
    ERROR( STREAMSERIESINPUTTESTC_EARG,
	   STREAMSERIESINPUTTESTC_MSGEARG, 0 );
    LALPrintError( USAGE, *argv );
    return STREAMSERIESINPUTTESTC_EARG;
  }

  /* Close files, cleanup, and exit. */
  if ( infile && strcmp( infile, "stdin" ) )
    fclose( fpIn );
  if ( outfile && strcmp( outfile, "stdout" ) )
    fclose( fpOut );
  LALCheckMemoryLeaks();
  INFO( STREAMSERIESINPUTTESTC_MSGENORM );
  return STREAMSERIESINPUTTESTC_ENORM;
}
/** \endcond */

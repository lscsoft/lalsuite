/*
 * Copyright (C) 2008 Karl Wette
 * Copyright (C) 2005 Reinhard Prix
 *
 *  [partially based on the MSG_LOG class in BOINC:
 *  Copyright (C) 2005 University of California]
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

/* Windows version fixed by Bernd Machenschalk */

#include <config.h>

/*---------- INCLUDES ----------*/
#include <stdio.h>
#include <string.h>
#include <strings.h>

#include <lal/XLALError.h>
#include <lal/LALMalloc.h>
#include <lal/LALDebugLevel.h>
#include <lal/Date.h>

#ifdef _MSC_VER
#include <Windows.h>
#include <process.h>
#define getpid _getpid
#else
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#endif

#ifdef HAVE_SYS_RESOURCE_H
#include <sys/resource.h>
#endif

#include <time.h>

#ifndef HAVE_LOCALTIME_R
#define localtime_r(timep, result) memcpy((result), localtime(timep), sizeof(struct tm))
#endif

#include <lal/LogPrintf.h>

/*---------- internal types ----------*/

/*---------- empty initializers ---------- */

/*---------- internal Global variables ----------*/


/*---------- internal prototypes ----------*/
static FILE* LogFile(void);

static const char * LogGetTimestamp (void);
static const char * LogTimeToString(double t);

static const char *LogFormatLevel( LogLevel_t level );

/*==================== FUNCTION DEFINITIONS ====================*/

/**
 * Get log level by examining lalDebugLevel
 */
LogLevel_t LogLevel(void)
{
  if (lalDebugLevel == 0)
    return LOG_NONE;		// If not printing LAL messages, also do not print log messages
  if (lalDebugLevel & LALINFO)
    return LOG_DETAIL;		// Print LOG_DETAIL messages if LAL_DEBUG_LEVEL contains 'info'
  if (lalDebugLevel & LALWARNING)
    return LOG_DEBUG;		// Print LOG_DEBUG messages if LAL_DEBUG_LEVEL contains 'warning'
  return LOG_NORMAL;		// Print LOG_CRITICAL and LOG_NORMAL messages by default
}

/** Decide where to print log messages */
static FILE* LogFile(void)
{
  if (LogLevel() < LOG_NORMAL)
    return stderr;		// Error log messages are printed to standard error
  return stdout;		// All other log messages are printed to standard output
}

/**
 * Verbatim output of the given log-message if LogLevel() >= level
 */
void
LogPrintfVerbatim (LogLevel_t level, const char* format, ...)
{

  if ( LogLevel() < level )
    return;

  va_list va;
  va_start(va, format);

  /* simply print this to output  */
  vfprintf (LogFile(), format, va );
  fflush(LogFile());

  va_end(va);

} /* LogPrintfVerbatim() */


/**
 * Output the given log-message, prefixed by a timestamp and level, if LogLevel() >= level
 */
void
LogPrintf (LogLevel_t level, const char* format, ...)
{

  if ( LogLevel() < level )
    return;

  va_list va;
  va_start(va, format);

  fprintf(LogFile(), "%s (%d) [%s]: ", LogGetTimestamp(), getpid(), LogFormatLevel(level) );
  vfprintf(LogFile(), format, va);
  fflush(LogFile());

  va_end(va);

} /* LogPrintf() */


/* taken from BOINC: return time-string for given unix-time
 */
const char *
LogTimeToString ( double t )
{
  static char buf[100];
  char finer[16];
  time_t x = (time_t)t;
  struct tm tm;
  localtime_r(&x, &tm);

  int hundreds_of_microseconds=(int)(10000*(t-(int)t));

  if (hundreds_of_microseconds == 10000) {
    /* paranoia -- this should never happen! */
    hundreds_of_microseconds=0;
    t+=1.0;
  }

  strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", &tm);
  sprintf(finer, ".%04d", hundreds_of_microseconds);
  strcat(buf, finer);

  return buf;

} /* LogTimeToString() */

///
/// Returns the peak amount of memory (in MB) allocated on the heap so far
/// using either lalMallocTotalPeak if memory-debugging is active or getrusage (if available),
/// otherwise returns -1 (without error) for "dont know"
///
/// \note the reported number always refers to the 'nominal' memory usage *not* including any extra "padding"
/// that would have been added by LALMalloc(), which typically doubles memory usage.
///
REAL8
XLALGetPeakHeapUsageMB ( void )
{
  // first see if lal's memory-debugging can be used
  if ( lalDebugLevel & LALMEMPADBIT ) {
    return lalMallocTotalPeak / (1024.0 * 1024.0);	// lalMallocTotalMax counts bytes, doesn't inlude 'padding'
  }

  // otherwise  try using getrusage
#ifdef HAVE_SYS_RESOURCE_H
  struct rusage usage;
  XLAL_CHECK_REAL8 ( getrusage ( RUSAGE_SELF, &usage ) == 0, XLAL_ESYS, "call to getrusage() failed with errno = %d\n", errno );
  REAL8 peakHeapMB = usage.ru_maxrss / 1024.0;	// maxrss is in KB
  if ( lalDebugLevel & LALMEMPADBIT ) {
    peakHeapMB /= 2.0;	// try to correct for memory-padding added by LALMalloc(), which seems ~factor of 2
  }
  // we're suspicious of a value of '0', which can also indicate that ru_maxrss is unsupported on this platform
  if ( usage.ru_maxrss > 0 ) {
    return peakHeapMB;
  }
#endif

  return -1;	// fallback answer: "dont know"

} // XLALGetMaxHeapUsageMB()

/**
 * Return time of day (seconds since 1970) as a double.
 * Taken from BOINC's dtime():
 *
 */
REAL8
XLALGetTimeOfDay ( void )
{
#ifdef _MSC_VER
  /* Windows version of dtime() from BOINC.
     Compile switch is MS compiler macro,
     because I suspect gettimeofday should be present in MinGW */
#define EPOCHFILETIME_SEC (11644473600.)
#define TEN_MILLION 10000000.

  LARGE_INTEGER time;
  FILETIME sysTime;
  double t;
  GetSystemTimeAsFileTime(&sysTime);
  time.LowPart = sysTime.dwLowDateTime;
  time.HighPart = sysTime.dwHighDateTime;  /* Time is in 100 ns units */
  t = (double)time.QuadPart;    /* Convert to 1 s units */
  t /= TEN_MILLION;                /* In seconds */
  t -= EPOCHFILETIME_SEC;     /* Offset to the Epoch time */
  return t;
#else

  struct timeval tv;
  if ( gettimeofday ( &tv, NULL ) != 0 ) {
    XLALPrintError ("%s: call to gettimeofday() failed with errno = %d\n", __func__, errno );
    XLAL_ERROR_REAL8 ( XLAL_ESYS );
  }

  return tv.tv_sec + (tv.tv_usec/1.e6);

#endif
} /* XLALGetTimeOfDay() */


///
/// High-resolution CPU timer (returns result in seconds), aimed for code-timing purposes.
/// Attempts to provide the highest time resolution available, while adding as little overhead as possible.
///
/// \note uses clock_gettime() with ns precision if available, or falls back to XLALGetTimeOfDay() with
/// mus-precision otherwise.
///
///
REAL8
XLALGetCPUTime ( void )
{
#ifndef HAVE_CLOCK_GETTIME
  return XLALGetTimeOfDay();
#else

  struct timespec ut;
  clockid_t clk_id;
#ifdef CLOCK_THREAD_CPUTIME_ID
  clk_id = CLOCK_THREAD_CPUTIME_ID;	// according to man-page: (since Linux 2.6.12)
#else
  clk_id = CLOCK_REALTIME;	// use this as fallback, guaranteed to exist.
#endif

  clock_gettime ( clk_id, &ut);	// don't bother testing to avoid overheads, and we would notice in timing if unavailable

  return ut.tv_sec + ut.tv_nsec * 1.e-9;

#endif
} // XLALGetCPUTime()


/* returns static timestamps-string for 'now' */
static const char *
LogGetTimestamp (void)
{
  return ( LogTimeToString ( XLALGetTimeOfDay() ) );

} /* LogGetTimestamp() */

/* returns static string characterizing this log-level */
static const char *
LogFormatLevel( LogLevel_t level )
{
  static char buf[100];

  switch ( level )
    {
    case LOG_CRITICAL:
      sprintf (buf, "CRITICAL" );
      break;
    case LOG_NORMAL:
      sprintf (buf, "normal");
      break;
    case LOG_DEBUG:
      sprintf (buf, "debug");
      break;
    case LOG_DETAIL:
      sprintf (buf, "detail");
      break;
    default:
      sprintf (buf, "unknown(%d)", level );
      break;
    } /* switch(level) */

  return buf;

} /* LogFormatLevel() */



/**
 * Output gsl_matrix in octave-format, using the given format for the matrix-entries
 * return -1 on error, 0 if OK.
 */
int
XLALfprintfGSLmatrix ( FILE *fp, const char *fmt, const gsl_matrix *gij )
{
  int cols, rows;
  int i, j;

  /* check user input */
  if ( !fp || !fmt )
    return -1;

  if ( gij == NULL )	/* treat as valid empty matrix */
    {
      fprintf ( fp, "[ ];\n");
      return 0;
    }

  rows = gij->size1;
  cols = gij->size2;

  fprintf (fp, " [ \\\n" );
  for ( i=0; i < rows; i ++ )
    {
      for (j=0; j < cols; j ++ )
	{
	  fprintf (fp, fmt, gsl_matrix_get ( gij, i, j ) );
	  if ( j < cols - 1 )
	    fprintf (fp, ", ");
	  else
	    fprintf (fp, ";\n");
	} /* for j < cols */
    }
  fprintf (fp, " ];\n" );

  return 0;

} /* XLALprintGSLmatrix() */


/**
 * Output gsl_matrix in octave-format, using the given format for the matrix-entries
 * return -1 on error, 0 if OK.
 */
int
XLALfprintfGSLvector ( FILE *fp, const char *fmt, const gsl_vector *vect )
{
  int rows;
  int i;

  /* check user input */
  if ( !vect || !fp || !fmt )
    return -1;

  rows = vect->size;

  fprintf (fp, " [ " );
  for ( i=0; i < rows; i ++ )
    {
      fprintf (fp, fmt, gsl_vector_get ( vect, i ) );
      if ( i < rows - 1 )
	fprintf (fp, ", ");
    } /* for i < rows */

  fprintf (fp, " ];\n" );

  return 0;

} /* XLALprintGSLvector() */

int
XLALfprintfGSLvector_int ( FILE *fp, const char *fmt, const gsl_vector_int *vect )
{
  int rows;
  int i;

  /* check user input */
  if ( !vect || !fp || !fmt )
    return -1;

  rows = vect->size;

  fprintf (fp, " [ " );
  for ( i=0; i < rows; i ++ )
    {
      fprintf (fp, fmt, gsl_vector_int_get ( vect, i ) );
      if ( i < rows - 1 )
	fprintf (fp, ", ");
    } /* for i < rows */

  fprintf (fp, " ];\n" );

  return 0;

} /* XLALprintGSLvector_int() */


/**
 * Returns input string with line-breaks '\n' removed (replaced by space)
 * The original string is unmodified. The returned string is allocated here.
 */
char *
XLALClearLinebreaks ( const char *str )
{
  char *ret, *tmp;

  if ( !str )
    return NULL;

  if ( (ret = LALMalloc( strlen ( str ) + 1) ) == NULL )
    return NULL;
  strcpy ( ret, str );

  tmp = ret;
  while ( (tmp = strchr ( tmp, '\n' ) ) )
    {
      *tmp = ' ';
      tmp ++;
    }

  return ret;

} /* XLALClearLinebreaks() */


/** dump given REAL4 time-series into a text-file */
int
XLALdumpREAL4TimeSeries ( const char *fname, const REAL4TimeSeries *series )
{
  XLAL_CHECK ( fname != NULL, XLAL_EINVAL );
  XLAL_CHECK ( series != NULL, XLAL_EINVAL );

  FILE *fp;
  XLAL_CHECK ( (fp = fopen (fname, "wb")) != NULL, XLAL_ESYS );

  REAL8 dt = series->deltaT;
  UINT4 numSamples = series->data->length;
  char buf[256];
  for ( UINT4 i = 0; i < numSamples; i++ )
  {
    LIGOTimeGPS ti = series->epoch;
    XLALGPSAdd ( &ti, i * dt );
    XLALGPSToStr ( buf, &ti );
    fprintf( fp, "%20s %20.16g\n", buf, series->data->data[i] );
  }
  fclose ( fp );

  return XLAL_SUCCESS;

} // XLALdumpREAL4TimeSeries()

/** dump given REAL8 time-series into a text-file */
int
XLALdumpREAL8TimeSeries ( const char *fname, const REAL8TimeSeries *series )
{
  XLAL_CHECK ( fname != NULL, XLAL_EINVAL );
  XLAL_CHECK ( series != NULL, XLAL_EINVAL );

  FILE *fp;
  XLAL_CHECK ( (fp = fopen (fname, "wb")) != NULL, XLAL_ESYS );

  REAL8 dt = series->deltaT;
  UINT4 numSamples = series->data->length;
  char buf[256];
  for ( UINT4 i = 0; i < numSamples; i++ )
  {
    LIGOTimeGPS ti = series->epoch;
    XLALGPSAdd ( &ti, i * dt );
    XLALGPSToStr ( buf, &ti );
    fprintf( fp, "%20s %20.16g\n", buf, series->data->data[i] );
  }
  fclose ( fp );

  return XLAL_SUCCESS;

} // XLALdumpREAL8TimeSeries()


/** dump given COMPLEX8 time-series into a text-file */
int
XLALdumpCOMPLEX8TimeSeries ( const char *fname, const COMPLEX8TimeSeries *series )
{
  XLAL_CHECK ( fname != NULL, XLAL_EINVAL );
  XLAL_CHECK ( series != NULL, XLAL_EINVAL );

  FILE *fp;
  XLAL_CHECK ( (fp = fopen (fname, "wb")) != NULL, XLAL_ESYS );

  REAL8 dt = series->deltaT;
  UINT4 numSamples = series->data->length;
  char buf[256];
  for ( UINT4 i = 0; i < numSamples; i++ )
  {
    LIGOTimeGPS ti = series->epoch;
    XLALGPSAdd ( &ti, i * dt );
    XLALGPSToStr ( buf, &ti );
    fprintf( fp, "%20s %20.16g %20.16g\n", buf, crealf(series->data->data[i]), cimagf(series->data->data[i]) );
  }
  fclose ( fp );

  return XLAL_SUCCESS;

} // XLALdumpCOMPLEX8TimeSeries()

/*
*  Copyright (C) 2003-9 Bruce Allen, Peter Shawhan
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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

/**
 * \file
 * \ingroup lalpulsar_bin_Tools
 * \author Bruce Allen, Peter Shawhan
 * \brief Multipulsar injection routine
 *
 * How it works
 * ============
 *
 * (adapted from \c README.s3inject by Bruce Allen)
 *
 * There are two executables: \c lalpulsar_Makefakedata_v4 and \c lalpulsar_hwinject.
 *
 * \c lalpulsar_Makefakedata_v4 generates a continuous stream of data for a single
 * pulsar, using a file called \c Pulsar.N to get the pulsar's parameters,
 * and command line arguments to get the starting GPS time and detector
 * name.  It can write either to files or to stdout.
 *
 * \c lalpulsar_hwinject starts up M copies of \c lalpulsar_Makefakedata_v4, each with a
 * different set of parameters, gathers their output, and adds it
 * together, writing it to stdout.
 *
 * Both of these executables support a \c -h command-line argument that
 * will give you a command-line argument summary description.  Try it
 * before reading further:
 * \code
 * lalpulsar_Makefakedata_v4 -h
 * lalpulsar_hwinject -h
 * \endcode
 *
 * In this document we assume you will inject 5 pulsar signals total.
 *
 * Create pulsar parameter files
 * -----------------------------
 *
 * Create a files called \c Pulsar.0 to \c Pulsar.4 containing the parameters
 * of the pulsars that you wish the inject.  These files should be
 * "almost" identical at the different sites.  They only differ because
 * the Actuation function differs between the sites.
 *
 * I have included five files - but note that the parameters to be used
 * "for real" still need to be determined.  People who should help
 * creating these parameter files include: Riles, S. Anderson,
 * G. Mendell, X. Siemens, G. Woan and M. Landry.  These files, once
 * created, should be made read-only, and should only be modified if the
 * instrument Actuation function changes.  Be sure to save a copy of them
 * someplace safe.
 *
 * Note: Each pulsar's parameters are defined at a particular fiducial
 * Solar System Barycenter (SSB) time.  In this document, I call this
 * time \f$ t_S \f$ .  The actual choice of this time is not important,
 * but it should be fixed and invariant.  The team defining pulsar
 * parameters may wish to change the value to something that they find
 * convenient.
 *
 * At the end of this document is a more detailed description of the
 * \c Pulsar.N parameter files.
 *
 * Create command-line argument files
 * ----------------------------------
 *
 * Create files called \c in.0 to \c in.4.  Each file should contain one line.
 * This is used to construct the command-line arguments for
 * \c lalpulsar_Makefakedata_v4.
 *
 * The file \c in.0 should look like this:
 * \code
 * lalpulsar_Makefakedata_v4 @Pulsar.0 -I LHO -S 751680013 -b
 * \endcode
 *
 * and the other \c in.N files should be identical except that they should
 * contain \c Pulsar.N.  The \c in.N files at LLO should contain \c LLO rather
 * than \c LHO.
 *
 * Verify setup
 * ------------
 *
 * To test your setup, do the following:
 * \code
 * lalpulsar_hwinject -n 5 -G 751680013 -s
 * \endcode
 *
 * The \c -s option is a show option.  It makes \c lalpulsar_hwinject read
 * the command line files \c in.N and show you the \b exact commands it would
 * actually run (preceeded by an integer count of the form <tt>[XX]</tt>).  The
 * \c -G command line argument is the GPS time at which to start producing
 * data.
 *
 * Now let's make some output (but just from one pulsar):
 * \code
 * lalpulsar_hwinject -n 1 -G 751680013 -T -X 2> infolog  | head -20
 * \endcode
 *
 * The <tt>2> infolog</tt> redirects stderr into an information log. Have a look.
 *
 * The \c -T option makes \c lalpulsar_hwinject output in human readable text
 * rather than binary
 *
 * The \c -X option makes \c lalpulsar_hwinject output an X axis as well.
 *
 * Notice that the first number output by \c lalpulsar_hwinject is \c always
 * 1234.5, which is a key to use in checking endian ordering and general
 * sanity.
 *
 * Now let's go "whole hog":
 * \code
 * lalpulsar_hwinject -n 5 -G 751680013 2> infolog  | od -w4 -f | more
 * \endcode
 *
 * This shows you the raw binary output in single column format
 *
 * CPU and resource use
 * --------------------
 *
 * On my 1 GHz PIII laptop, this job:
 * \code
 * lalpulsar_hwinject -n 5 -G 751680013 2> infolog  > /dev/null
 * \endcode
 *
 * runs at five time real time speed, has a starup latency of around 1
 * second, and uses 11 MB of memory per pulsar (55 MB total).
 *
 * Defining Pulsar.* files
 * -----------------------
 *
 * Here is a typical file:
 * \code
 * refTime          = 751680013            ## pulsar reference time in SSB frame
 * aPlus            = 4.996981857609109e-26        ## plus-polarization signal amplitude
 * aCross           = 4.868177666869495e-26        ## cross-polarization signal amplitude
 * psi              = 0.770087086          ## polarization angle psi (radians)
 * phi0             = 2.66                 ## phase at reference time
 * Freq             = 265.5771052          ## GW frequency at reference time
 * Delta            = -0.981180225         ## declination (radians)
 * Alpha            = 1.248816734          ## right ascension (radians)
 * f1dot            = -4.15E-12            ## 1st frequency derivative
 * f2dot            = 0.0                  ## 2nd frequency derivative
 * f3dot            = 0.0                  ## 3rd frequency derivative
 * \endcode
 *
 * This structure contains gravitational wave source position (in
 * Equatorial coordinates), and orientation angle.
 * Equatorial coordinates are the standard sky-fixed coordinate
 * system. The z-axis is defined as for geographic coordinates, above;
 * the plane orthogonal to this passing through the Earth's centre is
 * called the equator. The x-axis is defined to be the direction, as
 * viewed from the centre of the Earth, where the Sun appears to cross
 * the equator moving north in spring. This is called the vernal equinox,
 * and is shown in Fig. 16.6. In this coordinate system, the latitude
 * coordinate is called the declination \f$ \delta \f$ and the longitude
 * coordinate is called the right ascension \f$ \alpha \f$ .
 * (See also \ref SkyCoordinates_h.)
 *
 * \c aPlus and \c aCross set the amplitude of the two polarizations.  This is
 * where you can insert an amplitude calibration factor.
 *
 * The only other value to note is \c phi0.  This is where you can insert a
 * phase calibration offset, if needed.  Note: \c phi0 is \b not scaled by a
 * factor of two.  In other words, if you set <tt>phi=PI</tt>, you'll find the
 * output inverted.  If you set <tt>phi=PI/2</tt>, you'll see the output phase
 * retarted.  In other words, the peak of a particular cycle occurs one
 * quarter of a cycle \b earlier.
 *
 * To see this clearly, just do something like:
 * \code
 * lalpulsar_hwinject -n 1 -G 751680013 -T | head -100 > junk1
 * \endcode
 * Then change the value of \c phi0 in \c Pulsar.0, and do it again:
 * \code
 * lalpulsar_hwinject -n 1 -G 751680013 -T | head -100 > junk2
 * \endcode
 * Comparing \c junk1 and \c junk2 should make the sign convention of \c phi0 very
 * clear.
 *
 * Calibration lines
 * -----------------
 *
 * The \c lalpulsar_hwinject executable can inject up to three calibration
 * lines.  Here they are denoted by <tt>L==low</tt>, <tt>M==medium</tt> and
 * <tt>H==high</tt> to indicate the frequency.  They are defined by:
 * \code
 *   DARM = A_L sin(2 pi f_L (GPS-GPS_0)) +
 *          A_M sin(2 pi f_M (GPS-GPS_0)) +
 *          A_H sin(2 pi f_H (GPS-GPS_0))
 * \endcode
 * where <tt>GPS_0 = 751680013</tt>.  In the code, the frequencies are hardwired to:
 * \code
 *   f_L = 52.296875 = 52+1/4+1/32+1/64 Hz
 *   f_M = 166.6875  = 166+1/2+1/8+1/16 Hz
 *   f_H = 927.6875  = 927+1/2+1/8+1/16 Hz
 * \endcode
 * These can be changed, but (1) MUST be exactly represented as IEEE754
 * floats (not doubles) and (2) MUST be positive.
 *
 * The amplitudes of the injected lines are defined by three arguments (\c -L, \c -M and \c -H) to
 * \c lalpulsar_hwinject which set the amplitudes of the three lines.  The arguments are, for example:
 * \code
 * lalpulsar_hwinject -n 0 -T  -G 12345678 -L 17.76 | more
 * \endcode
 * will inject a calibration line at low frequency with an amplitude of
 * 17.76.  You can include any combination of \c -L, \c -M and \c -H.  If one of
 * these arguments is not present, then its assumed amplitude is zero.
 *
 * You can inject five pulsars plus three calibration lines with:
 * \code
 * lalpulsar_hwinject -n 5 -L 0.003 -M 0.0006 -H 0.8 -G 751680013 -T | more
 * \endcode
 *
 * Note: a good check that the calibration line injection code works
 * correctly is to compare the output with GPS times offset by integer
 * multiples of 64 seconds.  Since the smallest fractional part of the
 * frequencies above is 1/64 Hz, the calibration signal should repeat
 * exactly every 64 seconds.
 *
 * The \c -p option to \c lalpulsar_hwinject prints out the built-in
 * calibration line frequencies.
 *
 * Comments
 * --------
 *
 * I've tried to make \c lalpulsar_hwinject fairly robust.  In particular,
 * it catches \c SIGCHLD and if a child has been terminated (rather than
 * just being stopped) it tries to say why.  System related errors print
 * errno and it's interpretation.  Most error messages should come with a
 * PID to help figure out which process is going bad.
 *
 * Note that under Solaris, the pid returned with these error messages
 * appears to be the PID of the shell (started by popen) under which the
 * child was started.  This is usually one less than the PID of the
 * associated \c lalpulsar_Makefakedata_v4 process.
 *
 * If you send SIGUSR1 to \c lalpulsar_hwinject:
 * \code
 * kill -SIGUSR1 PID
 * \endcode
 * then it will report to stderr the amount of simulated data that it has
 * made (days/hours/minutes/seconds). Be careful \b not to send the signal
 * to the entire process group, since the children don't catch this
 * signal and will respond to it by terminating!
 *
 * The \c lalpulsar_hwinject program can be used to inject signals from sources \b other
 * than pulsars.  To use it in this way, your code must do the following:
 *
 *  1. write data to stdout, errors to stderr
 *
 *  2. take a command line argument which is GPS seconds and start its
 *     output at that time: <tt>-G secs</tt>
 *
 *  3. be able to run faster than real time under Solaris.  It should
 *     produce data at a sample rate of 16384 Hz in blocks of an integer
 *     number of seconds (called S below).
 *
 *  4. have the following internal structure.  Here \c S is the number of
 *     seconds, for example, 30, that your code uses internally to compute
 *     the next block of output.
 *     \code
 *     main() {
 *
 *         int length=S*16384;
 *
 *         float magic=1234.5;
 *         fwrite(&magic,  sizeof(float), 1, stdout);
 *
 *         while (1) {
 *
 *             // compute next output data, may be time-consuming ...
 *
 *             fwrite(&data, sizeof(float), length, stdout);
 *
 *             fflush(stdout);
 *         }
 *     }
 *     \endcode
 *     The <tt>fflush(stdout)</tt> is \b very important.  It ensures that your program
 *     will have time to compute its next block of data \b before it is needed
 *     by \c lalpulsar_hwinject.
 *  5. create an \c in.N file for your executable.  See description above.
 *
 * Changelog
 * =========
 *
 * - Multipulsar injection routine, written for E10/S3 by Bruce Allen, 10/2003.
 * - Calls to Signal Injection Library added by Peter Shawhan.
 * - 2005/02: Duncan Brown renamed variable to avoid Condor conflict.
 * - 2005/02: Reinhard Prix removed the actuation scaling options.
 * - 28 May 2009: renamed this code from 's3inject' to 'psinject' in lalsuite-GIT
 * - 19 March - 10 August: Bruce Allen, added ability to output FRAME format
 *   files for VIRGO real-time hardware injections.
 * - 2022/07: Karl Wette added test script, ported FRAME writing to LALFrame,
 *   add '--fmin' and '--Band' options to command line, moved code to LALPulsar,
 *   renamed lalpulsar_hwinject, RPM packaging
 */

#define _GNU_SOURCE   /* for SA_RESTART */

#include "config.h"

#include <stdio.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <errno.h>
#include <stdarg.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>
#include <math.h>
#include <sys/stat.h>

#ifdef ONLINE
#ifdef HAVE_GDS_DTT_SISTR_H
#include <gds/dtt/SIStr.h>
#else
#error SIStr header is not available
#endif
#endif

#include <lal/LALMalloc.h>
#include <lal/XLALError.h>
#include <lal/LALConstants.h>
#include <lal/LALgetopt.h>
#include <lal/LALPulsarVCSInfo.h>

#ifdef HAVE_LIBLALFRAME
#include <lal/Units.h>
#include <lal/TimeSeries.h>
#include <lal/LALFrameIO.h>
#endif

#define MAXPULSARS 64
/* blocksize for gathering results from children */
#define BLOCKSIZE 16384

/* Flag to initiate graceful shutdown, set by sighandler */
volatile int shutdown_pulsar_injection = 0;
int npulsars = 0;
float data[BLOCKSIZE];
float total[BLOCKSIZE];
float testval = 1234.5;
char *directory;
char *channel = NULL;
int  gpstime = 0;
char *programname;
int show = 0;
int do_axis = 0;
int do_text = 0;
int write_frames = 0;
long long count = 0;
long blocks = 0;
const char *ifo_name = NULL;
const char *actuation = NULL;
int sampling_rate = 16384;              /* sample rate of output requested by user (default = 16kHz) */
double starttime_offset_eff = 0;        /* effective starttime offset (in multiple of samples) requested by user*/

/* capacity, occupancy, and pointers to input buffers. Cleared to zero
   on starup because they are initialized data.  Units of buffer size
   and offsets are counted in floats!
*/
int bufflen[MAXPULSARS];
int readfrombuff[MAXPULSARS];
float *buffs[MAXPULSARS];

double calamp[3] = {0.0, 0.0, 0.0};

/*
   Calibration line frequencies. THE CODE ASSUMES THAT THESE ARE
   POSITIVE.  Order is Low, Medium, High frequency.  If you change
   these, you MUST choose a frequency that can be exactly represented
   as an IEEE754 floating point floats. For example
   52+1/4+1/32+1/64=52.296875
*/
double calfreq[3] = {
  52.296875,         /* 52+1/4+1/32+1/64 Hz */
  166.6875,          /* 166+1/2+1/8+1/16 Hz */
  927.6875           /* 927+1/2+1/8+1/16 Hz */
};

/*
  Calibration line fiducial time.
*/
int tfiducial = 751680013;

FILE *fp[MAXPULSARS];

/* forward declaration */
void syserror( int showerrno, const char *fmt, ... );
void sighandler( int sig );
void usage( FILE *filep );
int parseinput( int argc, char **argv );

/* signal handler to monitor the status of children catch SIGCHLD */
void sighandler( int sig )
{
  pid_t pid;
  int status;
  long long seconds = ( ( long long )blocks ) * ( ( long long )BLOCKSIZE ) / ( ( long long )sampling_rate );
  int minutes = seconds / 60;
  int hours = minutes / 60;
  int days = hours / 24;
  int secs = seconds % 60;

  syserror( 0, "made %d days %d hours %d minutes %d seconds of data\n",
            days, hours % 24, minutes % 60, secs );

  syserror( 0, "Received signal %d\n", sig );

  if ( sig == SIGTERM ) {
    syserror( 0, "received SIGTERM, initiating graceful shutdown...\n" );
    shutdown_pulsar_injection = 1;
  }

  if ( sig != SIGCHLD ) {
    return;
  }

  if ( ( pid = waitpid( -1, &status, WNOHANG | WUNTRACED ) ) > 0 ) {
    /* if child stopped, make log entry then return */
    if ( WIFSTOPPED( status ) ) {
      syserror( 0, "Subprocess [PID=%d] stopped because it caught signal %d [%s]\n",
                ( int )pid, WSTOPSIG( status ), strsignal( WSTOPSIG( status ) ) );
      return;
    }

    /* otherwise something more serious is wrong... */
    syserror( 0, "Subprocess [PID=%d] is misbehaving.\n", ( int )pid );
    if ( WIFEXITED( status ) ) {
      syserror( 0, "Subprocess [PID=%d] or shell did exit(%d)\n", ( int )pid, WEXITSTATUS( status ) );
    }

    if ( WIFSIGNALED( status ) )
      syserror( 0, "Subprocess [PID=%d] terminated because it caught signal %d [%s]\n",
                ( int )pid, WTERMSIG( status ), strsignal( WTERMSIG( status ) ) );
    shutdown_pulsar_injection = 1;
  } else {
    syserror( 1, "waitpid() returned -1.  Call Bruce...\n" );
  }
  return;
}

/* Like perror() but takes variable numbers of arguments and includes program name */
void syserror( int showerrno, const char *fmt, ... )
{
  char *thiserror = NULL;
  pid_t pid = getpid();
  time_t t = time( NULL );
  va_list ap;
  /* initialize variable argument list  */
  va_start( ap, fmt );
  /* print a standardized header with time and PID info */
  if ( showerrno && errno && ( thiserror = strerror( errno ) ) ) {
    fprintf( stderr, "%s [PID=%d] %.24s: %s: ", programname, ( int )pid, ctime( &t ), thiserror );
  } else {
    fprintf( stderr, "%s [PID=%d] %.24s: ", programname, ( int )pid, ctime( &t ) );
  }
  vfprintf( stderr, fmt, ap );
  va_end( ap );
  fflush( stderr );
  return;
}

/* usage message */
void usage( FILE *filep )
{
  fprintf( filep,
           "--------------------------------------------------------------------------------\n"
           "%s: \n"
           "--------------------------------------------------------------------------------\n"
           "Options are:\n"
           "-h            THIS help message\n"
           "-v            VCS ID information\n"
           "-n INT        Number of pulsars to simulate: 0, 1, 2, ..., %d\n"
           "-d DIRECTORY  Directory containing command line files in.0, in.1, ...\n"
           "              Default is: . (current working directory)\n"
           "-e CHANNEL    Inject excitation into CHANNEL\n"
           "-D            Turn on Debug output for routines in signal injection library\n"
           "-G INTGPS     GPS time (integer seconds) that will be passed on command line\n"
           "-T            Print human readable text rather than binary to stdout\n"
           "-X            Include X-axis in human-readable text output\n"
           "-s            Print commands that would be fork(2)ed, then exit(0)\n"
           "              ---------- --------------------------------------------------\n"
           "-L DOUBLE     | Inject  calibration lines. Here L,M,H denote Low/Mid/High |\n"
           "-M DOUBLE     | frequency.  The DOUBLE values specifies the corresponding |\n"
           "-H DOUBLE     | amplitude. If NOT given, the amplitude defaults to 0.0    |\n"
           "              ------------- -----------------------------------------------\n"
           "-p            Print the calibration line frequencies in Hz then exit(0)\n"
           "-I STRING     Detector: LHO, LLO, GEO, VIRGO, TAMA, CIT, ROME [REQUIRED]\n"
           "-A STRING     File containing detector actuation-function     [OPTIONAL]\n"
#ifdef HAVE_LIBLALFRAME
           "-F INT        Keep N frame files on disk.  If N==0 write all frames immediately.\n"
           "-S INT        Number of 1-second frames per frame file (default 60).\n"
#endif
           "-r INT        Sampling rate (NOTE: strain-generators must use the same!) (Default:16384)\n"
           "-z DOUBLE     Delay: shift CHANNEL signals by round[offset*samplingRate] samples forward (Default:0)\n"
           "--------------------------------------------------------------------------------\n"
           , programname, MAXPULSARS
         );
  return;
}


int parseinput( int argc, char **argv )
{

  int c;
  const char *optionlist = "hL:M:H:n:d:e:DG:TXspI:A:F:vS:r:z:";
  opterr = 0;

  double starttime_offset_req = 0;      /* requested offset correction to shift start-time of signals into the future */
  double starttime_offset_samples = 0;  /* offset rounded to nearest sample */

  /* set some defaults */
  directory = strdup( "." );

  programname = argv[0];

  while ( -1 != ( c = LALgetopt( argc, argv, optionlist ) ) ) {
    char *end;
    double tempamp;
    switch ( c ) {
    case 'v':
      if ( XLALOutputVCSInfo( stdout, lalPulsarVCSInfoList, 0, "%% " ) != XLAL_SUCCESS ) {
        XLALPrintError( "XLALOutputVCSInfo() failed!\n" );
        exit( 1 );
      }
      exit( 0 );
      break;

    case 'p':
      printf( "The calibration line frequencies are:\n"
              "  L: %18.14f Hz\n"
              "  M: %18.14f Hz\n"
              "  H: %18.14f Hz\n",
              calfreq[0], calfreq[1], calfreq[2] );
      exit( 0 );
      break;
    case 'h':
      /* usage message */
      usage( stdout );
      exit( 0 );
      break;

    case 'n':
      /* number of pulsars to simulate */
      npulsars = atoi( LALoptarg );
      if ( npulsars < 0 || npulsars > MAXPULSARS ) {
        syserror( 0, "%s: Number of pulsars (-n argument = %d) must be non-negative and less than %d\n",
                  argv[0], npulsars, MAXPULSARS + 1 );
        exit( 1 );
      }
      break;
    case 'L':
    case 'M':
    case 'H':
      /* calibration-line amplitude */
      tempamp = strtod( LALoptarg, &end );
      if ( *end ) {
        syserror( 1, "-%c %s is invalid. -%c takes a double-precision amplitude.\n", c, LALoptarg, c );
        exit( 1 );
      }
      /* assign amplitude to the correct line */
      if ( c == 'L' ) {
        calamp[0] = tempamp;
      } else if ( c == 'M' ) {
        calamp[1] = tempamp;
      } else {
        calamp[2] = tempamp;
      }
      break;
    case 'd':
      /* directory path */
      if ( directory ) {
        free( directory );
      }
      if ( !( directory = strdup( LALoptarg ) ) ) {
        syserror( 1, "Out of memory to duplicate -d directory path %s\n", LALoptarg );
        exit( 1 );
      }
      break;
    case 'e':
#ifdef ONLINE
      /* Excitation channel into which to inject */
      if ( !( channel = strdup( LALoptarg ) ) ) {
        syserror( 1, "Out of memory to duplicate -e channel name %s\n", LALoptarg );
        exit( 1 );
      }
#else
      syserror( 0, "The -e option to enable online signal injection requires compilation with -DONLINE\n" );
      exit( 1 );
#endif
      break;
    case 'D':
#ifdef ONLINE
      /* Turn on debugging for SIStr library routines */
      SIStr_debug++;
#else
      syserror( 0, "The -D option to enable SIStr debugging requires compilation with -DONLINE\n" );
      exit( 1 );
#endif
      break;
    case 'G':
      /* GPS time to pass as argument */
      gpstime = atoi( LALoptarg );
      break;
    case 'T':
      /* enable text output */
      do_text = 1;
      break;
    case 'X':
      /* include x axis in text output */
      do_axis = 1;
      break;
    case 's':
      /* show commands rather than executing them */
      show = 1;
      break;
    case 'I':
      ifo_name = LALoptarg;
      break;
    case 'A':
      actuation = LALoptarg;
      break;
    case 'F':
#ifdef HAVE_LIBLALFRAME
    {
      int how_many = atoi( LALoptarg );
      if ( how_many < 0 ) {
        syserror( 0, "%s: fatal error, argument -F %d must be non-negative.\n", argv[0], how_many );
        exit( 1 );
      }
      write_frames = 1 + how_many;
    }
#else
    syserror( 0, "%s: -F specified, but this binary was built without frame support.\n", argv[0] );
    exit( 1 );
#endif
    break;
    case 'r':
      sampling_rate = atoi( LALoptarg );
      if ( sampling_rate <= 0 ) {
        syserror( 0, "%s: need positive sampling rate! %d\n", argv[0], sampling_rate );
        exit( 1 );
      }
      if ( sampling_rate % 2 != 0 ) {
        syserror( 0, "%s: need sampling rate divisible by 2! %d\n", argv[0], sampling_rate );
        exit( 1 );
      }
#ifdef HAVE_LIBLALFRAME
      if ( BLOCKSIZE % sampling_rate != 0 ) {
        syserror( 0, "%s: sampling rate %d must divide evenly into block size %d for frame output\n", argv[0], sampling_rate, BLOCKSIZE );
        exit( 1 );
      }
#endif
      break;

    case 'z':
      starttime_offset_req = strtod( LALoptarg, &end );
      if ( end == LALoptarg ) {
        syserror( 1, "-%c %s is invalid. -%c takes a double-precision amplitude.\n", c, LALoptarg, c );
        exit( 1 );
      }
      break;

    default:
      /* error case -- option not recognized */
      syserror( 0, "%s: Option argument: -%c unrecognized or missing argument.\n"
                "\t\tPlease use the '-h' option to print usage message\n"
                , argv[0], optopt );
      exit( 1 );
      break;

    } /* switch(c) */

  } /* while (LALgetopt) */

  /* sanity checks on command line arguments */
  if ( do_axis && !do_text ) {
    syserror( 0, "The -X (axis) option only works with the -T (human readable) option\n" );
    exit( 1 );
  }
  if ( do_text && channel ) {
    syserror( 0, "Can't use both '-T' and '-e' together\n" );
    exit( 1 );
  }
  if ( channel && ( strlen( channel ) < 3 || channel[2] != ':' ) ) {
    syserror( 0, "Excitation channel %s not of form CC:CCC...\n", channel );
    exit( 1 );
  }
  if ( ifo_name == NULL ) {
    syserror( 0, "You must specify the IFO name (-I)\n" );
    exit( 1 );
  }


  /* "discretize" starttime offset to sampling rate, and
   * provide some debug-info about starttime shifting */
  starttime_offset_samples = floor( starttime_offset_req * sampling_rate + 0.5 );       /* correctly *round*, allowing for negative offsets */
  starttime_offset_eff     = starttime_offset_samples / sampling_rate;
  if ( starttime_offset_req ) {
    syserror( 0, "starttime OFFSET requested = %+.16g s (offset > 0 means a *delay*)\n", starttime_offset_req );
    syserror( 0, "effective OFFSET will be   = %+.16g s (that is: %+.0f samples)\n", starttime_offset_eff, starttime_offset_samples );
  }

#ifndef ONLINE
  if ( channel ) {
    syserror( 0, "Can't do exicitations. Code not compiled with ONLINE defined\n" );
    exit( 1 );
  }
#endif

  return 0;
}

#define MAXLINE 1024
int main( int argc, char *argv[] )
{
  int i, j;
#ifdef ONLINE
  char *cwd;
  int status;
  SIStream sis;
  char info[MAXLINE];
  memset( &sis, '\0', sizeof( SIStream ) );
#endif

#ifdef HAVE_LIBLALFRAME
  int framecounter = 0;
  REAL4TimeSeries *framesim = NULL;
#endif

  parseinput( argc, argv );

  syserror( 0, "Starting up\n" );

  /* install signal handler to catch SIGCHLD. Note that solaris uses
     SysV rather than BSD semantics and doesn't automaticall restart system
     calls like fread and fwrite.  So we use sigaction() rather than
     signal() so that we can set SA_RESTART. */
  {
    struct sigaction sig;
    memset( &sig, '\0', sizeof( sig ) );
    sig.sa_flags = SA_RESTART;
    sig.sa_handler = sighandler;
    if ( sigaction( SIGCHLD, &sig, NULL ) ) {
      syserror( 1, "Unable to install signal handler for messages about troubled children\n" );
    }
    if ( sigaction( SIGTERM, &sig, NULL ) ) {
      syserror( 1, "Unable to install signal handler for logging output rate data and terminating\n" );
    }
    if ( sigaction( SIGUSR1, &sig, NULL ) ) {
      syserror( 1, "Unable to install signal handler for logging output rate data\n" );
    }
  }

  for ( i = 0; i < npulsars; i++ ) {
    char command[MAXLINE];
    char filename[MAXLINE];
    FILE *fpc = NULL;
    int length = 0;
    char *newlineloc = NULL;

    /* construct file name */
    if ( snprintf( filename, MAXLINE, "%s/in.%d", directory, i ) > MAXLINE - 1 ) {
      syserror( 0, "%s: file name %s/in.%d has more than MAXLINE=%d characters\n",
                programname, directory, i, MAXLINE );
      exit( 1 );
    }

    /* open file */
    if ( !( fpc = fopen( filename, "r" ) ) ) {
      syserror( 1, "Can't open file %s for reading\n", filename );
      exit( 1 );
    }

    /* read command line from file */
    if ( !( fgets( command, MAXLINE, fpc ) ) ) {
      syserror( 1, "Command line file %s was empty!\n", filename );
      exit( 1 );
    }
    fclose( fpc );

    /* check that contents are not too large */
    length = strlen( command );
    if ( length >= MAXLINE - 1 ) {
      syserror( 0, "Command line file %s has line >= to MAXLINE=%d characters!\n",
                filename, MAXLINE );
      exit( 1 );
    }

    /* replace first NEWLINE to null terminate string */
    if ( ( newlineloc = strchr( command, '\n' ) ) ) {
      *newlineloc = '\0';
    }

    /* append additional arguments to command line */
    /* GPS starttime */
    length = strlen( command );
    if ( snprintf( command + length, MAXLINE - length, " -G %d", gpstime ) > MAXLINE - length - 1 ) {
      command[length] = '\0';
      syserror( 0, "%s: command line has >= MAXLINE=%d characters\n", programname, MAXLINE );
      exit( 1 );
    }
    /* IFO */
    length = strlen( command );
    if ( snprintf( command + length, MAXLINE - length, " -I %s", ifo_name ) > MAXLINE - length - 1 ) {
      command[length] = '\0';
      syserror( 0, "%s: command line has >= MAXLINE=%d characters\n", programname, MAXLINE );
      exit( 1 );
    }
    /* frequency band; zero to half the sampling rate */
    length = strlen( command );
    if ( snprintf( command + length, MAXLINE - length, " --fmin=0 --Band=%d", sampling_rate / 2 ) > MAXLINE - length - 1 ) {
      command[length] = '\0';
      syserror( 0, "%s: command line has >= MAXLINE=%d characters\n", programname, MAXLINE );
      exit( 1 );
    }
    /* Actuation-function if given */
    if ( actuation ) {
      length = strlen( command );
      if ( snprintf( command + length, MAXLINE - length, " --actuation=%s", actuation ) > MAXLINE - length - 1 ) {
        command[length] = '\0';
        syserror( 0, "%s: command line has >= MAXLINE=%d characters\n", programname, MAXLINE );
        exit( 1 );
      }
    }

    /* now either show the command or execute it */
    if ( show ) {
      printf( "[%02d] %s\n", i, command );
    } else {
      errno = 0;
      if ( !( fp[i] = popen( command, "r" ) ) ) {
        syserror( 1, "Unable to popen(3) %s\n", command );
        exit( 1 );
      }
    }
  }

  /* a useful option for debugging -- show the output */
  if ( show ) {
    if ( !npulsars ) {
      printf( "%s: Warning: n=0 so an infinite-length zero strength signal will be made!\n", argv[0] );
    }
    exit( 0 );
  }

#if 0
  {
    pid_t pid;
    int status;
    /* wait a couple of seconds, then check that all processes are running happily! */
    sleep( 2 );
    pid = waitpid( -1, &status, WNOHANG | WUNTRACED );
    if ( pid ) {
      syserror( 0, "Subprocess with PID=%d is misbehaving.\n", ( int )pid );
      if ( WIFEXITED( status ) ) {
        syserror( 0, "Subprocess or shell did exit(%d)\n", WEXITSTATUS( status ) );
      }

      if ( WIFSIGNALED( status ) )
        syserror( 0, "Subprocess terminated because it caught signal %d [%s]\n",
                  WTERMSIG( status ), strsignal( WTERMSIG( status ) ) );
      exit( 1 );
    }
  }
#endif

  /* processes opened, read data*/
  for ( i = 0; i < npulsars; i++ ) {
    if ( fread( &testval, sizeof( float ), 1, fp[i] ) != 1 ) {
      syserror( 1, "Could not read first float 1234.5 from %d'th signal source\n", i );
      exit( 1 );
    } else if ( testval != 1234.5 ) {
      syserror( 0, "First value (%f) from %d'th signal source was not 1234.5\n", testval,  i );
      exit( 1 );
    } else if ( fread( bufflen + i, sizeof( int ), 1, fp[i] ) != 1 ) {
      syserror( 1, "Could not read buffer size from %d'th signal source\n", i );
      exit( 1 );
    } else if ( bufflen[i] < sampling_rate || bufflen[i] > sampling_rate * 60 ) {
      syserror( 0, "Bad buffer size %d floats from %d'th signal source (expect %d <= size <= %d)\n", bufflen[i], i, sampling_rate, 60 * sampling_rate );
      exit( 1 );
    } else if ( bufflen[i] % BLOCKSIZE ) {
      syserror( 0, "Bad buffer size %d floats from %d'th signal source NOT a multiple of BLOCKSIZE=%d\n", bufflen[i], i, BLOCKSIZE );
      exit( 1 );
    } else if ( !( buffs[i] = ( float * )calloc( bufflen[i], sizeof( float ) ) ) ) {
      syserror( 1, "Can't allocate buffer of %d floats for %d'th signal source\n", bufflen[i], i );
      exit( 1 );
    }
    /* ensure that we read buffers on first pass */
    readfrombuff[i] = bufflen[i];
  }

  /* are we calling the excitation engine directly? */
  if ( channel ) {

#ifdef ONLINE
    /* Report some information about this injection client */
    cwd = getcwd( NULL, 256 );
    if ( cwd ) {
      sprintf( info, "%s %s", argv[0], cwd );
    } else {
      sprintf( info, "%s unknown_directory", argv[0] );
    }
    free( cwd );
    SIStrAppInfo( info );

    /* Open the Signal Injection Stream */
    status = SIStrOpen( &sis, channel, sampling_rate, ( double ) gpstime + starttime_offset_eff );
    if ( SIStr_debug ) {
      syserror( 0, "SIStrOpen() returned %d\n", status );
    }
    if ( status != SIStr_OK ) {
      syserror( 0, "SIStrOpen() error opening SIStream: %s\n", SIStrErrorMsg( status ) );
      exit( 2 );
    }
#endif
  } else if ( do_text ) {
    printf( "1234.5\n" );
  } else {
    /* write out 1234.5 */
    if ( 1 != fwrite( &testval, sizeof( float ), 1, stdout ) ) {
      syserror( 1, "Unable to output key value 1234.5\n" );
      exit( 1 );
    }
  }

  /* now read data blocks unless a SIGTERM has set shutdown */
  while ( !shutdown_pulsar_injection ) {
    int num = 0;
    int line;
    int tdelt = gpstime - tfiducial;

    /* clear block that will contain totals */
    for ( j = 0; j < BLOCKSIZE; j++ ) {
      total[j] = 0.0;
    }

    /* if needed, insert calibration line(s) */
    for ( line = 0; line < 3; line++ ) {
      if ( calamp[line] != 0.0 ) {
        /* normal int and double variables for integer/fractional
           time.  In this and in the code that follows, _fra refers to
           the fractional part [0,1) and _int refers to the integer
           part. */

        double dt_fra;

        /*  Avoid long longs in inner loop as they are slow. */
        long long t_rem = blocks, t_int = BLOCKSIZE;

        /* line amplitude and frequency (integer + fractional parts) */
        double f_fra  = calfreq[line];
        int    f_int  = ( int )f_fra;
        f_fra -= f_int;

        /* integer and fractional time offsets of first sample */
        t_rem   *= t_int;
        t_int    = t_rem;
        t_int   /= sampling_rate;
        t_rem   -= t_int * sampling_rate;
        t_int   += tdelt;

        // unused: int dt_int   = t_int;
        dt_fra   = t_rem;
        dt_fra  /= sampling_rate;

        for ( j = 0; j < BLOCKSIZE; j++ ) {
          double cycles1, cycles2, cycles3;
          double tlocal_fra  = dt_fra + ( double )j / ( ( double )sampling_rate );
          int    tlocal_int  = ( int )tlocal_fra;
          tlocal_fra        -= tlocal_int;
          tlocal_int        += t_int;
          cycles1            = f_fra * tlocal_int;
          cycles1           -= ( int )cycles1;
          cycles2            = tlocal_fra * f_int;
          cycles2           -= ( int )cycles2;
          cycles3            = tlocal_fra * f_fra;
          cycles3           -= ( int )cycles3;

          total[j] = calamp[line] * sin( 2 * LAL_PI * ( cycles1 + cycles2 + cycles3 ) );
        }
      }
    }

    /* loop over the different pulsars */
    for ( i = 0; i < npulsars; i++ ) {
      float *where;

      if ( readfrombuff[i] == bufflen[i] ) {
        /* read data from the i'th signal */
        if ( bufflen[i] != ( num = fread( buffs[i], sizeof( float ), bufflen[i], fp[i] ) ) ) {
          syserror( 1, "Only read %d floats (not %d) from %d'th signal source\n", num, bufflen[i], i );
          exit( 1 );
        }
#ifdef DEBUGTIMING
        syserror( 0, "just read %d seconds of data from %d'th signal\n", num / sampling_rate, i );
#endif
        readfrombuff[i] = 0;
      }

      /* location of signal in buffer */
      where = buffs[i] + readfrombuff[i];

      /* add i'th pulsar to the total signal */
      for ( j = 0; j < BLOCKSIZE; j++ ) {
        total[j] += where[j];
      }

      readfrombuff[i] += BLOCKSIZE;
    }

    /* now output the total signal to frames */

    if ( write_frames ) {
#ifdef HAVE_LIBLALFRAME

      const int secs_per_framefile = BLOCKSIZE / sampling_rate;

      /* allocate time series */
      if ( !framesim ) {
        const char frName[] = "CW_simulated";
        framesim = XLALCreateREAL4TimeSeries( frName, 0, 0, 1.0 / sampling_rate, &lalDimensionlessUnit, BLOCKSIZE );
        if ( !framesim ) {
          syserror( 1, "XLALCreateREAL4TimeSeries() failed" );
          exit( 1 );
        }
      }

      /* set up GPS time, copy data */
      framesim->epoch.gpsSeconds = gpstime + framecounter * secs_per_framefile;
      for ( size_t m = 0; m < BLOCKSIZE; m++ ) {
        framesim->data->data[m] = total[m];
      }

      /* create frame */
      const char framename[] = "CW_Injection";
      LALFrameH *frame = XLALFrameNew( &framesim->epoch, secs_per_framefile, framename, 1, framecounter, 0 );
      if ( !frame ) {
        syserror( 1, "XLALFrameNew() failed" );
        exit( 1 );
      }

      /* add data to frame */
      if ( XLALFrameAddREAL4TimeSeriesSimData( frame, framesim ) != XLAL_SUCCESS ) {
        syserror( 1, "!XLALFrameAddREAL4TimeSeriesSimData() failed" );
        exit( 1 );
      }

      /* write data to framefile */
      char framefilename[256];
      snprintf( framefilename, sizeof( framefilename ), "%s-%d-%d.gwf", framename, framesim->epoch.gpsSeconds, secs_per_framefile );
      if ( XLALFrameWrite( frame, framefilename ) ) {
        syserror( 1, "Error during frame write" );
        exit( 1 );
      }

      /* free memory for frames and for simdata structures */
      XLALFrameFree( frame );

      /* Do we keep a limited set of frames on disk? */
      if ( write_frames > 1 ) {
        char listname[256];
        int watchtime = gpstime + secs_per_framefile * ( framecounter / secs_per_framefile - write_frames + 1 );
        sprintf( listname, "%s-%d-%d.gwf", framename, watchtime, secs_per_framefile );
        /* syserror(0, "Watching for file %s to disappear....\n", listname); */
        struct stat statbuf;
        while ( !stat( listname, &statbuf ) ) {
          /* if enough files already in place, then sleep 0.1 seconds */
          struct timespec rqtp;
          rqtp.tv_sec = 0;
          rqtp.tv_nsec = 100000000;
          nanosleep( &rqtp, NULL );
        }
      }

      /* increment counter for the next second */
      framecounter++;
#else
      syserror( 0, "ERROR: write_frames!=0, but binary was built without Frame support\n" );
      exit( 1 );
#endif
    } /* if (write_frames) */

    /* now output the total signal... */
    else if ( channel ) {
#ifdef ONLINE
      /* ... to the excitation engine ... */
      status = SIStrAppend( &sis, total, BLOCKSIZE, 1.0 );
      if ( SIStr_debug >= 2 ) {
        syserror( 0, "SIStrAppend() returned %d\n", status );
      }
      if ( status != SIStr_OK ) {
        syserror( 0, "SIStrAppend() error adding data to stream: %s\n",
                  SIStrErrorMsg( status ) );
        break;
      }
#endif
    }
    /* ... or as text ... */
    else if ( do_text ) {
      if ( do_axis ) {
        /* ... either as text with both x and y axes ... */
        long long x1 = gpstime;
        long long E9 = 1000000000;
        x1 *= E9;
        x1 += ( long long )( starttime_offset_eff * E9 ); /* account for startime-shift, consistent with CHANNEL injection */

        for ( j = 0; j < BLOCKSIZE; j++ ) {
          long long x2 = count, x3;
          x2 *= E9;
          x2  /= ( sampling_rate );
          x2 += x1;
          x3 =  x2;
          x3 /= E9;
          x2 -= x3 * E9;
          printf( "%lld.%09lld %g\n", x3, x2, total[j] );
          count++;
        }
      } else {
        /* ... or as y-axis text only ... */
        for ( j = 0; j < BLOCKSIZE; j++ ) {
          printf( "%g\n", total[j] );
        }
      }
    } else {
      /* ... or in raw binary form. */
      if ( BLOCKSIZE != ( num = fwrite( total, sizeof( float ), BLOCKSIZE, stdout ) ) ) {
        syserror( 1, "Only wrote %d values (not %d)\n", num, BLOCKSIZE );
        exit( 1 );
      }
#ifdef DEBUGTIMING
      syserror( 0, "Just sent %d seconds of data to system\n", BLOCKSIZE / sampling_rate );
      sleep( BLOCKSIZE / sampling_rate );
#endif
    }

    /* increment counter of blocks sent out */
    blocks++;
  }

  /* We'll be exiting, so clean up */
  if ( channel ) {
#ifdef ONLINE
    /* Close the signal injection stream */
    status = SIStrClose( &sis );
    if ( SIStr_debug ) {
      syserror( 0, "SIStrClose returned %d\n", status );
    }
    if ( status != SIStr_OK ) {
      syserror( 0, "Error while closing SIStream: %s\n", SIStrErrorMsg( status ) );
      exit( 2 );
    }
#endif
  }

#ifdef _LINUX
  /* shut down signal handler for SIGCHLD */
  {
    struct sigaction sig;
    memset( &sig, '\0', sizeof( sig ) );
    sig.sa_flags = SA_RESTART | SA_NOCLDSTOP;
    sig.sa_handler = SIG_IGN;
    if ( sigaction( SIGCHLD, &sig, NULL ) ) {
      syserror( 1, "Unable to install signal handler for exiting\n" );
    }
    if ( sigaction( SIGPIPE, &sig, NULL ) ) {
      syserror( 1, "Unable to install signal handler for exiting\n" );
    }
  }

  for ( i = 0; i < npulsars; i++ ) {
    int status;

    __fpurge( fp[i] );
    status = pclose( fp[i] );
    if ( status != -1 ) {
      syserror( 0, "The %d'th signal generator did exit(%d)\n", i, ( int )WEXITSTATUS( status ) );
    } else {
      syserror( 1, "Trouble shutting down the %d'th signal generator\n", i );
      exit( 1 );
    }
  }
#endif

#ifdef HAVE_LIBLALFRAME
  XLALDestroyREAL4TimeSeries( framesim );
#endif

  LALCheckMemoryLeaks();

  /* and exit cleanly */
  syserror( 0, "Shutdown complete, exiting cleanly\n" );
  exit( 0 );
}

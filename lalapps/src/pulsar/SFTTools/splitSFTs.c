/*
 *  Copyright (C) 2020 Karl Wette
 *  Copyright (C) 2008, 2010 Bernd Machenschalk, Bruce Allen
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
 * \file
 * \ingroup lalapps_pulsar_SFTTools
 * \author Bernd Machenschalk, Bruce Allen
 *
 * \brief This program reads in binary SFTs (v2) and writes out narrow-banded merged SFTs (v2).
 *
 * Writen by Bernd Machenschalk for Einstein\@home 2008
 *
 * Revision splitSFTs.c,v 1.41 2008/10/29 16:54:13 was
 * reviewed by the LSC CW Review Committee Thur Oct 30, 2008
 *
 * Suggested improvements:
 * issue warning at ambiguous frequency values that are so close to the boundary of a bin
 * that rounding might end up giving an unintended bin
 * check for consistency of input SFTs (same timebase, ascending timestamps etc., see spec)
 * and merged SFTs (check last header of a file we are appending to)
 *
 * Other possible improvements not suggested by the committee
 * keep output files open (if there aren't too many)
 * obscure a mystery factor in command-line record even if given with long option --factor
 */

#include <config.h>

#include <math.h>
#include <time.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <glob.h>
#include <sys/stat.h>
#include <lal/LALStdlib.h>
#include <lal/LALHashTbl.h>
#include <lal/Date.h>
#include <lal/SFTfileIO.h>
#include <lal/SFTReferenceLibrary.h>
#include <LALAppsVCSInfo.h>

/**
 * rounding (for positive numbers!)
 * taken from SFTfileIO in LALSupport, should be consistent with that
 */
#define MYROUND(x) ( floor( (x) + 0.5 ) )

#define FALSE 0
#define TRUE (!FALSE)

#define CMT_NONE 0
#define CMT_OLD  1
#define CMT_FULL 2

// Compare two quantities, and return a sort order value if they are unequal
#define COMPARE_BY( x, y ) do { if ( (x) < (y) ) return -1; if ( (x) > (y) ) return +1; } while(0)

typedef struct {
  /* hashed fields */
  char det[2];         /**< detector name */
  long timebase;       /**< SFT timebase */
  long firstbin;       /**< first bin */
  long binwidth;       /**< bin width */
  /* value fields */
  long int nSFT;       /**< number of SFTs */
  long startTime;      /**< GPS start time */
  char *filename;      /**< SFT filename */
} SFT_RECORD;

typedef struct {
  int resource_units;  /**< number of resource units available for consumption */
  int resource_rate;   /**< this many resource units become available in one second */
  time_t last_checked; /**< time we last checked */
} UNIT_SOURCE;

/** throtteling settings */
UNIT_SOURCE read_bandwidth = {0, 0, 0};
UNIT_SOURCE read_open_rate = {0, 0, 0};
UNIT_SOURCE write_bandwidth = {0, 0, 0};
UNIT_SOURCE write_open_rate = {0, 0, 0};

/** request a resurce. Function returns after waiting for throttle time */
static void request_resource( UNIT_SOURCE *us, int units )
{
  time_t now;
  int seconds;

  /* negative rate indicates no throttling */
  if ( us->resource_rate <= 0 ) {
    return;
  }

  us->resource_units -= units;
  while ( us->resource_units <= 0 ) {
    time( &now );
    seconds = now - us->last_checked;
    /* guard against overflow and condor checkpointing */
    if ( seconds < 0 ) {
      seconds = 0;
    }
    if ( seconds > 10 ) {
      seconds = 10;
    }
    if ( seconds == 0 ) {
      sleep( 1 );
    }
    us->resource_units += seconds * us->resource_rate;
    us->last_checked = now;
  }
}

/* determine if the given filepath is an existing directory or not */
static int is_directory( const char *fname )
{
  struct stat stat_buf;

  if ( stat( fname, &stat_buf ) ) {
    return 0;
  }

  if ( ! S_ISDIR( stat_buf.st_mode ) ) {
    return 0;
  } else {
    return 1;
  }

}

/* hash an SFT_RECORD */
static UINT8 hash_SFT_RECORD( const void *x )
{
  const SFT_RECORD *rec = ( const SFT_RECORD * ) x;
  UINT4 hval = 0;
  XLALPearsonHash( &hval, sizeof( hval ), &rec->det, sizeof( rec->det ) );
  XLALPearsonHash( &hval, sizeof( hval ), &rec->timebase, sizeof( rec->timebase ) );
  XLALPearsonHash( &hval, sizeof( hval ), &rec->firstbin, sizeof( rec->firstbin ) );
  XLALPearsonHash( &hval, sizeof( hval ), &rec->binwidth, sizeof( rec->binwidth ) );
  return hval;
}

/* compare an SFT_RECORD */
static int compare_SFT_RECORD( const void *x, const void *y )
{
  const SFT_RECORD *recx = ( const SFT_RECORD * ) x;
  const SFT_RECORD *recy = ( const SFT_RECORD * ) y;
  COMPARE_BY( recx->det[0], recy->det[0] );
  COMPARE_BY( recx->det[1], recy->det[1] );
  COMPARE_BY( recx->timebase, recy->timebase );
  COMPARE_BY( recx->firstbin, recy->firstbin );
  COMPARE_BY( recx->binwidth, recy->binwidth );
  return 0;
}

/* destroy an SFT_RECORD */
static void destroy_SFT_RECORD( void *x )
{
  SFT_RECORD *rec = ( SFT_RECORD * ) x;
  if ( rec != NULL ) {
    XLALFree( rec->filename );
    XLALFree( rec );
  }
}

static int move_to_next_SFT ( FILE *fpin, int nsamples, int comment_length )
{
  int move;
  move = sizeof( struct headertag2 ) + nsamples * 2 * sizeof( float ) + comment_length;
  fseek( fpin, move, SEEK_CUR );
  return XLAL_SUCCESS;
}

/** main program */
int main( int argc, char **argv )
{
  int arg;                        /* current command-line argument */
  unsigned int bin;               /* current bin */
  struct headertag2 hd, lasthd;   /* header of input SFT */
  FILE *fpin;                     /* currently open input filepointer */
  FILE *fpout;                    /* currently open output filepointer */
  char *oldcomment;               /* comment of input SFT */
  char *cmdline = NULL;           /* records command-line to add it to comment */
  char *comment = NULL;           /* comment to be written into output SFT file */
  int swap;                       /* do we need to swap bytes? */
  float *data;                    /* SFT data */
  char *outname;                  /* name of output SFT file */
  char outdir0[] = ".";
  char *outdir = ( char * )outdir0; /* output filename prefix */
  char *detector = NULL;          /* detector name */
  double factor = 1.0;            /* "mystery" factor */
  int firstfile = TRUE;           /* are we processing the first input SFT file? */
  int allcomments = FALSE;        /* write comment into _every_ SFT in the file */
  int add_comment = CMT_FULL;     /* add VCS ID and full command-line to every SFT file */
  unsigned int startBin = 0, endBin = 0;     /* start and end in bins */
  unsigned int width = 0, overlap = 0; /* width and overlap in bins */
  double fMin = -1.0, fMax = -1.0;     /* start and end in Hz */
  double fWidth = -1.0, fOverlap = -1; /* width and overlap in Hz */
  unsigned int nactivesamples;         /* number of bins to actually read in */
  int validate = TRUE;                 /* validate the checksum of each input SFT before using it? */
  LIGOTimeGPS *minStartTime = NULL;    /* pointer for optional constraint */
  LIGOTimeGPS *maxStartTime = NULL;    /* pointer for optional constraint */
  LIGOTimeGPS XLAL_INIT_DECL(userMinStartTime); /* user-given constraint: earliest SFT timestamp to start using */
  LIGOTimeGPS XLAL_INIT_DECL(userMaxStartTime); /* user-given constraint: earliest SFT timestamp to no longer use */
  int assumeSorted = 0;                /* Are SFT input files chronologically sorted? */
  int sfterrno = 0;                    /* SFT error number return from reference library */
  LALHashTbl *nbsfts = NULL;           /* hash table of existing narrow-band SFTs */

  /* initialize throtteling */
  time( &read_bandwidth.last_checked );
  time( &read_open_rate.last_checked );
  time( &write_bandwidth.last_checked );
  time( &write_open_rate.last_checked );

  /* help / usage message */
  if ( ( argv[1] == NULL ) ||
       ( strcmp( argv[1], "-h" ) == 0 ) ||
       ( strcmp( argv[1], "--help" ) == 0 ) ) {
    fprintf( stderr,
             "%s -h\n"
             "\n"
             "  Write this help message\n"
             "\n"
             "%s\n"
             "  [-a|--all-comments]\n"
             "  [-c|--add-comment 0|1|2]\n"
             "  [-v|--no-validation]\n"
             "  [-s|--start-bin <startbin>]\n"
             "  [-e|--end-bin <endbin (exclusively)>]\n"
             "  [-b|--width <sftbins>]\n"
             "  [-x|--overlap <overlap>]\n"
             "  [-fs|--start-frequency <startfrequency>]\n"
             "  [-fe|--end-frequency <endfrequency (exclusively)>]\n"
             "  [-fb|--frequency-bandwidth <frequencywidth>]\n"
             "  [-fx|--frequency-overlap <frequencyoverlap>]\n"
             "  [-ts|--minStartTime <minStartTime>]\n"
             "  [-te|--maxStartTime <maxStartTime (exclusively)>]\n"
	     "  [-as|--assumeSorted 1|0]\n"
             "  [-m|--factor <factor>]\n"
             "  [-d|--detector <detector>]\n"
             "  [-n|--output-directory <outputdirectory>]\n"
             "  [--] <inputfile> ...\n"
             "\n"
             "  This program reads in binary SFTs (v2) and writes out narrow-banded\n"
             "  merged SFTs (v2).\n"
             "\n"
             "  The frequency bands of the ouput SFTs (first frequency bin of first output SFT,\n"
             "  last frequency bin of last output SFT, number of bins in each output SFT)\n"
             "  and a possible overlap of the output files can be specified\n"
             "  either in bins ('-s', '-e', '-b', -'x') or Hz ('-fs', '-fe', '-fb', '-fx')\n"
             "  (or mixed - if both are given, frequency values take precedence).\n"
             "\n"
             "  A 'mystery factor' can be specified with '-m' option.\n"
             "\n"
             "  '-v' skips validation of the CRC checksum of the input files\n"
             "\n"
             "  The name of the output SFT follows the standard SFT filename convention, with\n"
             "  information specific to narrow-band SFTs included in the description:\n"
             "    <outdir>/<obs>-<nSFT>_<det>_<timebase>SFT_NB_F<firstbin>_W<binwidth>-<start>-<span>.sft\n"
             "  where\n"
             "    <outdir> = output directory given by the '-o' option, default is '.'\n"
             "    <nSFT> = number of SFTs in the file\n"
             "    <obs> = 1-character observatory prefix, e.g. 'H', 'L', 'V'\n"
             "    <det> = 2-character detector prefix, e.g. 'H1', 'L1', 'V1'\n"
             "    <timebase> = time base of the SFTs, rounded to the nearest second\n"
             "    <firstbin> = first bin of the SFTs, expressed in the form\n"
             "      <firstbin> = <freq>Hz<remainder>, where\n"
             "        <freq> = <firstbin> / <timebase>, rounded down to nearest Hz\n"
             "        <remainder> = <firstbin> - <freq> * <timebase>, number of remaining bins\n"
             "    <binwidth> = bin width of the SFTs, expressed in the form\n"
             "      <binwidth> = <freq>Hz<remainder>, where\n"
             "        <freq> = <binwidth> / <timebase>, rounded down to nearest Hz\n"
             "        <remainder> = <binwidth> - <freq> * <timebase>, number of remaining bins\n"
             "    <start> = start time of the first SFT as a GPS timestamp\n"
             "    <span> = time spanned by the SFT in seconds\n"
             "  If an output file already exists, the program will append the new SFTs to them, making\n"
             "  it possible to construct the final narrow-band SFTs by running the program multiple\n"
             "  times with different input SFTs. The GPS timestamps of the input SFTs need to be in\n"
             "  ascending order to get valid merged SFT files.\n"
             "\n"
             "  The '-c' options specifies how to deal with comments - 0 means no comment is written\n"
             "  at all, 1 means that the comment is taken unmodified from the input SFTs, 2 (default)\n"
             "  means that the program appends its RCS id and command-line to the comment.\n"
             "  By default a comment is written only to the first SFT of a merged SFT output 'block'\n"
             "  (i.e. call to this program). Adding the option '-a' to the command line specifies that\n"
             "  instead the comment is written into every SFT in the resulting file.\n"
             "\n"
             "  The options 'minStartTime' and 'maxStartTime' constrain input SFTs to be read\n"
             "  according to their time stamps. Following the XLALCWGPSinRange() convention,\n"
             "  these are interpreted as half-open intervals:\n"  
             "  \tminStartTime <= time stamp < maxStartTime.\n"
             "  After encountering a time stamp outside of the given range, the program jumps to\n"
             "  the next input file as specified in the input argument; i.e. time stamp constraints\n"
             "  are imposed file-wise. For a multi-SFT file, headers are checked to ensure proper\n" 
             "  sorting. Mind that all timestamps refer to SFT **start** times.\n" 
             "\n"
             "  If 'assumeSorted' is set to 1, SFT input files will be assumed to be chronologically\n"
             "  sorted, which means the program will stop as soon as an SFT located after the\n" 
             "  specified range is encountered.\n"
             "\n"
             "  After all options (and an optional '--' separator), the input files are given, as many\n"
             "  as you wish (or the OS supports - using xargs should be simple with this command-line\n"
             "  syntax).\n"
             "\n"
             "  The program adds its own VCS ID and command-line to the comment of the written SFTs,\n"
             "  a mystery factor should show up as 'xxx' there.\n",
             argv[0], argv[0] );
    exit( 0 );
  }

  /* record VCS ID and command-line for the comment */
  XLAL_CHECK_MAIN( ( cmdline = ( char * )XLALMalloc( strlen( lalAppsVCSIdentInfo.vcsId ) + strlen( lalAppsVCSIdentInfo.vcsStatus ) + 2 ) ) != NULL, XLAL_ENOMEM, "out of memory allocating cmdline" );
  strcpy( cmdline, lalAppsVCSIdentInfo.vcsId );
  strcat( cmdline, lalAppsVCSIdentInfo.vcsStatus );
  strcat( cmdline, "\n" );
  for ( arg = 0; arg < argc; arg++ ) {
    if ( strcmp( argv[arg], "-m" ) == 0 ) {
      /* obscure the mystery factor */
      XLAL_CHECK_MAIN( ( cmdline = ( char * )XLALRealloc( ( void * )cmdline, strlen( cmdline ) + 8 ) ) != NULL, XLAL_ENOMEM, "out of memory allocating cmdline" );
      strcat( cmdline, "-m xxx " );
      arg++;
    } else {
      XLAL_CHECK_MAIN( ( cmdline = ( char * )XLALRealloc( ( void * )cmdline, strlen( cmdline ) + strlen( argv[arg] ) + 2 ) ) != NULL, XLAL_ENOMEM, "out of memory allocating cmdline" );
      strcat( cmdline, argv[arg] );
      if ( arg == argc - 1 ) {
        strcat( cmdline, "\n" );
      } else {
        strcat( cmdline, " " );
      }
    }
  }

  /* get parameters from command-line */
  for ( arg = 1; arg < argc; arg++ ) {
    if ( ( strcmp( argv[arg], "-d" ) == 0 ) ||
         ( strcmp( argv[arg], "--detector" ) == 0 ) ) {
      detector = argv[++arg];
    } else if ( ( strcmp( argv[arg], "-c" ) == 0 ) ||
                ( strcmp( argv[arg], "--add-comment" ) == 0 ) ) {
      add_comment = atoi( argv[++arg] );
    } else if ( ( strcmp( argv[arg], "-a" ) == 0 ) ||
                ( strcmp( argv[arg], "--all-comments" ) == 0 ) ) {
      allcomments = TRUE;
    } else if ( ( strcmp( argv[arg], "-s" ) == 0 ) ||
                ( strcmp( argv[arg], "--start-bin" ) == 0 ) ) {
      startBin = atoi( argv[++arg] );
    } else if ( ( strcmp( argv[arg], "-e" ) == 0 ) ||
                ( strcmp( argv[arg], "--end-bin" ) == 0 ) ) {
      endBin = atoi( argv[++arg] );
    } else if ( ( strcmp( argv[arg], "-b" ) == 0 ) ||
                ( strcmp( argv[arg], "--width" ) == 0 ) ) {
      width = atoi( argv[++arg] );
    } else if ( ( strcmp( argv[arg], "-x" ) == 0 ) ||
                ( strcmp( argv[arg], "--overlap" ) == 0 ) ) {
      overlap = atoi( argv[++arg] );
    } else if ( ( strcmp( argv[arg], "-fs" ) == 0 ) ||
                ( strcmp( argv[arg], "--start-frequency" ) == 0 ) ) {
      fMin = atof( argv[++arg] );
    } else if ( ( strcmp( argv[arg], "-fe" ) == 0 ) ||
                ( strcmp( argv[arg], "--end-frequency" ) == 0 ) ) {
      fMax = atof( argv[++arg] );
    } else if ( ( strcmp( argv[arg], "-fb" ) == 0 ) ||
                ( strcmp( argv[arg], "--frequency-bandwidth" ) == 0 ) ) {
      fWidth = atof( argv[++arg] );
    } else if ( ( strcmp( argv[arg], "-fx" ) == 0 ) ||
                ( strcmp( argv[arg], "--frequency-overlap" ) == 0 ) ) {
      fOverlap = atof( argv[++arg] );
    } else if ( ( strcmp( argv[arg], "-ts" ) == 0 ) ||
                ( strcmp( argv[arg], "--minStartTime" ) == 0 ) ) {
      XLAL_CHECK_MAIN ( XLALGPSSetREAL8 ( &userMinStartTime, (REAL8)atof( argv[++arg] )) != NULL, XLAL_EDOM );
      minStartTime = &userMinStartTime;
    } else if ( ( strcmp( argv[arg], "-te" ) == 0 ) ||
                ( strcmp( argv[arg], "--maxStartTime" ) == 0 ) ) {
      XLAL_CHECK_MAIN ( XLALGPSSetREAL8 ( &userMaxStartTime, (REAL8)atof( argv[++arg] )) != NULL, XLAL_EDOM );
      maxStartTime = &userMaxStartTime;
    } else if ( ( strcmp( argv[arg], "-as" ) == 0 ) ||
                ( strcmp( argv[arg], "--assumeSorted" ) == 0 ) ) {
      assumeSorted = atoi( argv[++arg] );
    } else if ( ( strcmp( argv[arg], "-m" ) == 0 ) ||
                ( strcmp( argv[arg], "--factor" ) == 0 ) ) {
      factor = atof( argv[++arg] );
    } else if ( ( strcmp( argv[arg], "-v" ) == 0 ) ||
                ( strcmp( argv[arg], "--no-validation" ) == 0 ) ) {
      validate = FALSE;
    } else if ( ( strcmp( argv[arg], "-n" ) == 0 ) ||
                ( strcmp( argv[arg], "--output-directory" ) == 0 ) ) {
      outdir = argv[++arg];
    } else if ( ( strcmp( argv[arg], "-rb" ) == 0 ) ||
                ( strcmp( argv[arg], "--read-bandwidth" ) == 0 ) ) {
      read_bandwidth.resource_rate = atoi( argv[++arg] );
    } else if ( ( strcmp( argv[arg], "-ror" ) == 0 ) ||
                ( strcmp( argv[arg], "--read-open-rate" ) == 0 ) ) {
      read_open_rate.resource_rate = atoi( argv[++arg] );
    } else if ( ( strcmp( argv[arg], "-wb" ) == 0 ) ||
                ( strcmp( argv[arg], "--write-bandwidth" ) == 0 ) ) {
      write_bandwidth.resource_rate = atoi( argv[++arg] );
    } else if ( ( strcmp( argv[arg], "-wor" ) == 0 ) ||
                ( strcmp( argv[arg], "--write-open-rate" ) == 0 ) ) {
      write_open_rate.resource_rate = atoi( argv[++arg] );
    } else if ( strcmp( argv[arg], "--" ) == 0 ) {
      ++arg;
      break;
    } else if ( strncmp( argv[arg], "-", 1 ) == 0 ) {
      fprintf( stderr, "unknown option '%s', try '-h' for help\n", argv[arg] );
      exit( -1 );
    } else {
      break;
    }
  }

  /* build an output string for the user timestamps constraint */
  char *constraint_str;
  const size_t constraint_str_len = 32; /* just needs to be big enough for two INT4 and a comma */
  XLAL_CHECK_MAIN( ( constraint_str = ( char * )XLALMalloc( constraint_str_len ) ) != NULL,
                   XLAL_ENOMEM, "out of memory allocating constraint_str" );
  if ( minStartTime ) {
    XLAL_CHECK_MAIN ( snprintf( constraint_str, constraint_str_len, "%d,", minStartTime->gpsSeconds ) < ( int )constraint_str_len,
                      XLAL_ESYS, "output of snprintf() was truncated" );
  }
  else {
    XLAL_CHECK_MAIN ( snprintf( constraint_str, constraint_str_len, ",") < ( int )constraint_str_len,
                      XLAL_ESYS, "output of snprintf() was truncated" );
  }
  if ( maxStartTime ) {
    XLAL_CHECK_MAIN ( snprintf( constraint_str + strlen(constraint_str), constraint_str_len - strlen(constraint_str), "%d", maxStartTime->gpsSeconds ) < ( int )(constraint_str_len - strlen(constraint_str)),
                      XLAL_ESYS, "output of snprintf() was truncated" );
  }

  /* check if there was an input-file option given at all */
  XLAL_CHECK_MAIN( argv[arg] != NULL, XLAL_EINVAL, "no input files specified" );

  /* check output directory exists */
  XLAL_CHECK_MAIN( is_directory( outdir ), XLAL_ESYS, "output directory does not exist" );

  /* allocate space for output filename */
  const size_t outnamelen = strlen( outdir ) + 256;
  XLAL_CHECK_MAIN( ( outname = ( char * )XLALMalloc( outnamelen ) ) != NULL, XLAL_ENOMEM, "out of memory allocating outname" );

  /* find existing narrow-band SFTs */
  nbsfts = XLALHashTblCreate( destroy_SFT_RECORD, hash_SFT_RECORD, compare_SFT_RECORD );
  XLAL_CHECK_MAIN( nbsfts != NULL, XLAL_EFUNC, "XLALHashTblCreate() failed" );
  {
    XLAL_CHECK_MAIN( snprintf( outname, outnamelen, "%s/*SFT_NB*.sft", outdir ) < ( int )outnamelen,
                     XLAL_ESYS, "output of snprintf() was truncated" );
    glob_t globbuf;
    glob( outname, GLOB_MARK | GLOB_NOSORT, NULL, &globbuf );
    for ( size_t i = 0; i < globbuf.gl_pathc; ++i ) {

      /* find the filename component */
      char *filename = strrchr( globbuf.gl_pathv[i], '/' );
      XLAL_CHECK_MAIN( filename != NULL, XLAL_ESYS, "could not extract filename from '%s'", globbuf.gl_pathv[i] );
      ++filename;
      XLAL_CHECK_MAIN( *filename != '\0', XLAL_ESYS, "'%s' is a directory", globbuf.gl_pathv[i] );

      /* create new SFT record */
      SFT_RECORD *rec = NULL;
      XLAL_CHECK_MAIN( ( rec = XLALCalloc( 1, sizeof( *rec ) ) ) != NULL, XLAL_ENOMEM, "out of memory allocating SFT_RECORD" );

      /* parse SFT filename */
      char obs = 0;
      long int firstbinfreq = 0, firstbinrem = 0, binwidthfreq = 0, binwidthrem = 0, span = 0;
      int sscanf_matched = sscanf( filename, "%c-%ld_%c%c_%ldSFT_NB_F%ldHz%ld_W%ldHz%ld-%ld-%ld.sft",
                                   &obs, &rec->nSFT, &rec->det[0], &rec->det[1], &rec->timebase, &firstbinfreq,
                                   &firstbinrem, &binwidthfreq, &binwidthrem, &rec->startTime, &span );
      XLAL_CHECK_MAIN( sscanf_matched == 11, XLAL_ESYS,
                       "sscanf() matched only %i fields in SFT filename '%s'; partial match = '%c-%ld_%c%c_%ldSFT_NB_F%ldHz%ld_W%ldHz%ld-%ld-%ld.sft'",
                       sscanf_matched, filename,
                       obs, rec->nSFT, rec->det[0], rec->det[1], rec->timebase, firstbinfreq,
                       firstbinrem, binwidthfreq, binwidthrem, rec->startTime, span );
      XLAL_CHECK_MAIN( obs == rec->det[0], XLAL_EINVAL, "inconsistent observatory (%c) vs detector (%c%c) in SFT filename '%s'", obs, rec->det[0], rec->det[1], filename );
      XLAL_CHECK_MAIN( span > 0, XLAL_EINVAL, "nonpositive timespan (%ld) in SFT filename '%s'", span, filename );

      /* complete SFT record */
      rec->firstbin = firstbinfreq * rec->timebase + firstbinrem;
      rec->binwidth = binwidthfreq * rec->timebase + binwidthrem;
      XLAL_CHECK_MAIN( ( rec->filename = XLALStringDuplicate( globbuf.gl_pathv[i] ) ) != NULL, XLAL_ENOMEM, "out of memory allocating filename" );

      /* add SFT record */
      XLALPrintInfo( "Existing SFT: det=%c%c timebase=%ld firstbin=%ld binwidth=%ld filename=%s\n",
                     rec->det[0], rec->det[1], rec->timebase, rec->firstbin, rec->binwidth, rec->filename );
      XLAL_CHECK_MAIN( XLALHashTblAdd( nbsfts, rec ) == XLAL_SUCCESS, XLAL_EFUNC, "XLALHashTblAdd() failed" );

    }
    globfree( &globbuf );
  }

  int stopInputCauseSorted=0;
  /* loop over all input SFT files */
  for ( ; (arg < argc) && !stopInputCauseSorted; arg++ ) {

    /* open input SFT file */
    request_resource( &read_open_rate, 1 );
    XLAL_CHECK_MAIN( ( fpin = fopen( argv[arg], "r" ) ) != NULL, XLAL_EIO, "could not open SFT file for reading" );

    /* loop until end of input SFT file, when ReadSFTHeader() will return SFTENONE */
    int firstread = TRUE;
    int nSFT_this_file = 0;
    while ( 1 ) {

      /* read header, break if ReadSFTHeader() returns SFTENONE
         (unless this is the first read, i.e. fail is SFT file is completely empty) */
      request_resource( &read_bandwidth, 40 );
      sfterrno = ReadSFTHeader( fpin, &hd, &oldcomment, &swap, validate );
      if ( !firstread ) {
        if ( sfterrno == SFTENONE ) {
          break;
        }
      }
      XLAL_CHECK_MAIN( sfterrno == 0, XLAL_EIO, "could not read SFT header: %s", SFTErrorMessage( sfterrno ) );

      /* Check that various bits of header information are consistent.
       * This includes a check for monotonic increasing timestamps.
       */
      XLAL_CHECK_MAIN( firstread || !(sfterrno=CheckSFTHeaderConsistency(&lasthd, &hd)), XLAL_EIO, "Inconsistent SFT headers: %s", SFTErrorMessage(sfterrno));
      lasthd = hd; /* keep copy of header for comparison the next time */

      /* only use SFTs within specified ranges */
      LIGOTimeGPS startTimeGPS = {hd.gps_sec,0};
      int inRange = XLALCWGPSinRange(startTimeGPS, minStartTime, maxStartTime);
      XLAL_CHECK_MAIN ( !( firstread && ( inRange == 1 ) ), XLAL_EIO, "First timestamp %d in file was after user constraint [%s)!", hd.gps_sec, constraint_str );
      firstread = FALSE;
      if ( inRange == -1 ) { /* input timestamp too early, skip */
        XLAL_CHECK_MAIN ( move_to_next_SFT ( fpin, hd.nsamples, hd.comment_length ) == XLAL_SUCCESS, XLAL_EIO );
        continue;
      }
      if ( inRange == 1 ) { /* input timestamp too late, stop with this file */
        stopInputCauseSorted = assumeSorted;
        break;
      }
      nSFT_this_file += 1;

      /* calculate bins from frequency parameters if they were given */
      /* deltaF = 1.0 / tbase; bins = freq / deltaF => bins = freq * tbase */
      if ( fMin >= 0.0 ) {
        startBin = MYROUND( fMin * hd.tbase );
      }
      if ( fMax >= 0.0 ) {
        endBin   = MYROUND( fMax * hd.tbase );
      }
      if ( fWidth >= 0.0 ) {
        width = MYROUND( fWidth * hd.tbase );
      }
      if ( fOverlap >= 0.0 ) {
        overlap = MYROUND( fOverlap * hd.tbase );
      }

      /* allocate space for SFT data */
      XLAL_CHECK_MAIN( ( data = ( float * )XLALCalloc( hd.nsamples, 2 * sizeof( float ) ) ) != NULL, XLAL_ENOMEM, "out of memory allocating data" );

      /* error if desired start bin < hd.firstfreqindex */
      if ( ( int )startBin < hd.firstfreqindex ) {
        fprintf( stderr,
                 "ERROR: start bin (%d) is smaller than first bin in input SFT (%d)\n",
                 startBin, hd.firstfreqindex );
        exit( 9 );
      }

      /* error if desired end bin > hd.firstfreqindex + hd.nsamples - 1 */
      if ( startBin + width > endBin + 1 ) {
        fprintf( stderr,
                 "ERROR: end bin (%d) is larger than last bin in input SFT (%d)\n",
                 endBin, hd.firstfreqindex + hd.nsamples - 1 );
        exit( 10 );
      }

      /* error if overlap is larger than the width */
      if ( overlap >= width ) {
        fprintf( stderr,
                 "ERROR: overlap (%d) is not smaller than the width (%d)\n",
                 overlap, width );
        exit( 11 );
      }

      /* construct comment for output SFTs */
      if ( add_comment > CMT_OLD ) {

        /* allocate space for new comment */
        XLAL_CHECK_MAIN( ( comment = ( char * )XLALMalloc( hd.comment_length + strlen( cmdline ) + 1 ) ) != NULL, XLAL_ENOMEM, "out of memory allocating comment" );

        /* append the commandline of this program to the old comment */
        if ( oldcomment ) {
          strcpy( comment, oldcomment );
        } else {
          *comment = '\0';
        }
        strcat( comment, cmdline );

      } else if ( add_comment == CMT_OLD ) {

        /* only copied existing comment, no additional space needed */
        comment = oldcomment;

      } /* else (add_comment == CMT_NONE) and (comment == NULL) i.e. no comment at all */

      /* get the detector name from SFT header */
      detector = hd.detector;

      /* calculate number of bins to actually read (from width + overlap) */
      /* add width-overlap samples as lon as they are < the total number og bins to write */
      for ( nactivesamples = 0; nactivesamples < endBin - startBin; nactivesamples += width - overlap );
      /* add the last overlap */
      nactivesamples += overlap + 1;
      /* if this number is larger than the bins in the input sft, just use the latter */
      if ( nactivesamples > hd.nsamples + hd.firstfreqindex - startBin ) {
        nactivesamples = hd.nsamples + hd.firstfreqindex - startBin;
      }

      /* read in SFT bins */
      request_resource( &read_bandwidth, nactivesamples * 8 );
      sfterrno = ReadSFTData( fpin, data, startBin, nactivesamples, NULL, NULL );
      XLAL_CHECK_MAIN( sfterrno == 0, XLAL_EIO, "could not read SFT data: %s", SFTErrorMessage( sfterrno ) );

      /* apply mystery factor and possibly normalization factor */
      for ( bin = 0; bin < 2 * nactivesamples; bin++ ) {
        data[bin] *= factor;
      }

      /* loop over start bins for output SFTs */
      for ( bin = startBin; bin < endBin; bin += width - overlap ) {
        /* determine the number of bins actually to write from the desired 'width',
           given that the remaining number of bin may be odd (especially from overlapping)
           and the bins to write need to be present in the input sft
        */
        int last_input_bin   = hd.firstfreqindex + hd.nsamples - 1;
        int last_output_bin  = bin + width - 1;
        int max_input_width  = last_output_bin <= last_input_bin ? width : width - ( last_output_bin - last_input_bin );
        int max_output_width = endBin - bin + 1;
        int this_width       = max_input_width < max_output_width ? max_input_width : max_output_width;

        /* check if this narrow-band SFT exists */
        SFT_RECORD key = {
          .det = { detector[0], detector[1] },
          .timebase = ( int )MYROUND( hd.tbase ),
          .firstbin = bin,
          .binwidth = this_width
        };
        XLALPrintInfo( "Looking for SFT: det=%c%c timebase=%ld firstbin=%ld binwidth=%ld\n",
                       key.det[0], key.det[1], key.timebase, key.firstbin, key.binwidth );
        SFT_RECORD *rec = NULL;
        XLAL_CHECK_MAIN( XLALHashTblExtract( nbsfts, &key, ( void ** ) &rec ) == XLAL_SUCCESS, XLAL_EFUNC, "XLALHashTblExtract() failed" );
        if ( rec == NULL ) {

          /* create new SFT record */
          XLAL_CHECK_MAIN( ( rec = XLALCalloc( 1, sizeof( *rec ) ) ) != NULL, XLAL_ENOMEM, "out of memory allocating SFT_RECORD" );
          *rec = key;
          rec->nSFT = 0;
          rec->startTime = hd.gps_sec;

        }

        /* update number of SFTs */
        rec->nSFT += 1;

        /* construct filename for this output SFT */
        int outfreq = ( int )floor( bin / rec->timebase );
        int outfreqbin = bin - outfreq * rec->timebase;
        int outwidth = ( int )floor( this_width / rec->timebase );
        int outwidthbin = this_width - outwidth * rec->timebase;
        long int span = hd.gps_sec - rec->startTime + rec->timebase;
        if ( hd.gps_nsec > 0 ) {
          span += 1;
        }
        XLAL_CHECK_MAIN( snprintf( outname, outnamelen, "%s/%c-%ld_%c%c_%ldSFT_NB_F%04dHz%d_W%04dHz%d-%ld-%ld.sft",
                                   outdir, detector[0], rec->nSFT, detector[0], detector[1],
                                   rec->timebase, outfreq, outfreqbin, outwidth, outwidthbin,
                                   rec->startTime, span ) < ( int )outnamelen,
                         XLAL_ESYS, "output of snprintf() was truncated" );

        /* update SFT filename */
        if ( rec->filename == NULL ) {
          rec->filename = outname;
          XLAL_CHECK_MAIN( ( outname = ( char * )XLALMalloc( outnamelen ) ) != NULL, XLAL_ENOMEM, "out of memory allocating outname" );
        } else if ( strcmp( rec->filename, outname ) != 0 ) {
          XLALPrintInfo( "Renaming SFT '%s' to '%s'\n", rec->filename, outname );
          rename( rec->filename, outname );
          char *const tmp = rec->filename;
          rec->filename = outname;
          outname = tmp;
        }

        /* store new SFT record */
        XLALPrintInfo( "New/updated SFT: det=%c%c timebase=%ld firstbin=%ld binwidth=%ld filename=%s\n",
                       rec->det[0], rec->det[1], rec->timebase, rec->firstbin, rec->binwidth, rec->filename );
        XLAL_CHECK_MAIN( XLALHashTblAdd( nbsfts, rec ) == XLAL_SUCCESS, XLAL_EFUNC, "XLALHashTblAdd() failed" );

        /* append this SFT */
        request_resource( &write_open_rate, 1 );
        XLAL_CHECK_MAIN( ( fpout = fopen( rec->filename, "a" ) ) != NULL, XLAL_EIO, "could not open SFT for writing" );

        /* write the data */
        /* write the comment only to the first SFT of a "block", i.e. of a call of this program */
        request_resource( &write_bandwidth, 40 + this_width * 8 );
        sfterrno = WriteSFT( fpout, hd.gps_sec, hd.gps_nsec, hd.tbase, bin, this_width, detector, ( firstfile || allcomments ) ? comment : NULL, data + 2 * ( bin - startBin ) );
        XLAL_CHECK_MAIN( sfterrno == 0, XLAL_EIO, "could not write SFT data: %s", SFTErrorMessage( sfterrno ) );

        /* close output SFT file */
        fclose( fpout );

      } /* loop over output SFTs */

      /* cleanup */
      if ( add_comment > CMT_OLD ) {
        XLALFree( comment );
      }
      XLALFree( data );

      /* next file is not the first file anymore */
      firstfile = FALSE;

      /* Move forward to next SFT in merged file */
      XLAL_CHECK_MAIN ( move_to_next_SFT ( fpin, hd.nsamples, hd.comment_length ) == XLAL_SUCCESS, XLAL_EIO );

    } /* end loop over SFTs in this file */

    XLAL_CHECK_MAIN ( nSFT_this_file > 0, XLAL_EIO, "Found no matching SFTs in file." );

    /* close the input SFT file */
    fclose( fpin );

  } /* loop over input SFT files */

  /* cleanup */
  XLALFree( outname );
  XLALFree( constraint_str );
  if ( add_comment > CMT_OLD ) {
    XLALFree( cmdline );
  }
  XLALHashTblDestroy( nbsfts );
  LALCheckMemoryLeaks();

  return ( 0 );
}

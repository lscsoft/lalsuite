/*=============================================================================
ivana - Inspiral Veto ANAlysis
Written June-July 2002 by Peter Shawhan

To compile/link on Solaris:
gcc -g ivana.c -I$LIGOTOOLS/include -L$LIGOTOOLS/lib -ldataflow -lsocket -lnsl\
    -o ivana

Usage example:
ivana ../data/H2_ASQ_triple_c5.xml \
    '../data/H2_MICHCTRL_triple.xml(-0.3,0.3,1.5)' triple_ranges.txt

=============================================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "metaio.h"

/* CLUSWINDOW determines how candidate events are grouped into clusters.  It
   is the maximum time between consecutive candidate events. */
#define CLUSWINDOW 1.0
/* MINLIVETIMEINTERVAL is the time associated with Patrick's "veto clustering"
   criterion.  It is the minimum length of a time interval between veto
   events which can be considered live. */
#define MINLIVEINTERVAL 4.0

#define RINGBUFSIZE 16384
#define MAXTIMERANGES 4096
#define MAXVETOFILES 8

/*====== Structure definitions ==============================================*/
typedef struct Range_struct {
  double t1, t2; /*-- Endpoints of time range (relative to timeBase) --*/
  double secs;   /*-- Length of time range in seconds --*/
  double dead;   /*-- Seconds of deadtime --*/
  int nCand, nCandPass, nCandFail;
  int nClus, nClusPass, nClusFail;
  int nVeto, nVetoUsed;
} Range;

typedef struct VetoFile_struct {
  char filename[256];
  struct MetaioParseEnvironment parseEnv;
  MetaioParseEnv env;
  FILE *datFile;
  double winNeg;
  double winPos;
  float snrRatio;
  int colS, colNS, colDur, colSnr;
  double tLast, tLastNeg, tLastPos;
  float lastSnr;
  int readit;
  int eof;
} VetoFile;


/*====== Function prototypes ================================================*/
void PrintUsage();
void InitRange( Range* range, double t1, double t2 );
void AbortAll( MetaioParseEnv candEnv, VetoFile vetoFile[], int nvetofiles,
	       MetaioParseEnv outEnv );

/*===========================================================================*/
int main( int argc, char **argv )
{
  char candFile[256];
  char vetospec[1024];
  char outFile[256] = "";
  char intext[256], temptext[256];
  char* arg;
  char* nextptr;
  char* chptr;
  char* chptr2;
  char* tokptr;
  FILE* fd;
  int timeBase = 0;
  double tTemp1, tTemp2;

  Range range[MAXTIMERANGES];
  Range *cRange = &range[0];
  Range *vRange = &range[0];
  Range *totals;
  int nRange=0, iRange=0, iRangeV=0, jRange;

  VetoFile vetoFile[MAXVETOFILES];
  VetoFile *vFile;
  int nvetofiles = 0;
  int ivfile, jvfile;

  double tCand, tCandLast=-999.0, tDeadNeg=-2.0e9, tDeadPos=-2.0e9;
  float snrCand;
  double tdead1, tdead2;
  int iCandS, iCandNS, iCandSnr;
  int status, status2, ostatus, candeof=0;
  int iveto, iarg, iposarg=0, ipos, pass, clusPass;
  char ttext[64];
  double secfrac;
  int allPast = 0;
  double pct;

  int debug=0;
  double dur;
  double tLastNeg = 0.0;
  double tUseNeg = 0.0, tUsePos;
  float snrThresh;
  int usevfile;

  struct MetaioParseEnvironment candParseEnv, outParseEnv;
  const MetaioParseEnv candEnv = &candParseEnv;
  const MetaioParseEnv outEnv = &outParseEnv;
  struct MetaioTable *table;
  FILE *candDatFile = NULL;

  /* Ring buffer of veto events */
  int vbufW = 0;
  int vbufR = 0;
  double tVetoNeg[RINGBUFSIZE];
  double tVetoPos[RINGBUFSIZE];
  float cSnrThresh[RINGBUFSIZE];
  Range* rVeto[RINGBUFSIZE];
  int usedVeto[RINGBUFSIZE];

  /*------ Beginning of code ------*/

  /*------ Initialize some things ------*/
  candEnv->fileRec.fp = NULL;
  for ( ivfile=0; ivfile<MAXVETOFILES; ivfile++ ) {
    vetoFile[ivfile].env = &(vetoFile[ivfile].parseEnv);
    vetoFile[ivfile].env->fileRec.fp = NULL;
    vetoFile[ivfile].datFile = NULL;
    vetoFile[ivfile].tLast = 0.0;
    vetoFile[ivfile].tLastNeg = 0.0;
    vetoFile[ivfile].tLastPos = 0.0;
    vetoFile[ivfile].lastSnr = 0.0;
    vetoFile[ivfile].readit = 1;
    vetoFile[ivfile].eof = 0;
  }
  outEnv->fileRec.fp = NULL;

  /*------ Parse command-line arguments ------*/

  for ( iarg=1; iarg<argc; iarg++ ) {
    arg = argv[iarg];

    /*-- First check for a flag (a dash followed by a non-digit) --*/
    if ( *arg == '-' && strchr("0123456789",arg[1])==NULL ) {
      switch (arg[1]) {
      case 'd':
	debug++;
	break;
      default:
	printf( "Invalid flag: %s\n", arg );
	PrintUsage(); return 1;
      }
      continue;
    }

    /*-- If we get here, then we have a positional argument --*/
    iposarg++;
    switch (iposarg) {

    case 1:
      /*-- Candidate event file --*/

      strncpy( candFile, arg, sizeof(candFile) );
      if ( candFile[sizeof(candFile)-1] != '\0' ) {
	printf( "Error: candidate file name is too long\n" ); return 1;
      }

      break;

    case 2:
      /*-- Veto specification (comma-sep list of filename/window pairs) --*/

      /*-- If veto specification is blank, skip it --*/
      if ( *arg == '\0' ) {
	break;
      }

      do {

	vFile = &(vetoFile[nvetofiles]);

	/*-- Split off the window specification (in parentheses) --*/
	chptr = strchr( arg, '(' );
	if ( chptr == NULL || arg[strlen(arg)-1] != ')' ) {
	  printf ( "Each item in the veto specification must consist of"
		   " filename followed, in\nparentheses, by"
		   " (negative,positive) window limits\n" );
	  PrintUsage(); return 1;
	}
	/*-- Null-terminate just the filename part --*/
	*chptr = '\0'; chptr++;

	/*-- Find the end of this item --*/
	chptr2 = strchr( chptr, ')' );
	if ( chptr2 == NULL ) {
	  printf ( "Missing close-paren in veto specification\n" );
	  PrintUsage(); return 1;
	}
	chptr2++;
	/*-- At this point, chptr2 points to the next item (if any) --*/
	switch ( *chptr2 )
	{
	case '\0':
	  /*-- No more veto specification items in list --*/
	  nextptr = NULL;
	  break;
	case ',':
	  /*-- There's another veto specification item in the list --*/
	  nextptr = chptr2 + 1;
	  break;
	default:
	  printf ( "Items in list of veto specifications must be separated"
		   " by a comma\n" );
	  PrintUsage(); return 1;
	}

	/*-- Copy the filename --*/
	strncpy( vFile->filename, arg, sizeof(vFile->filename) );
	if ( vFile->filename[sizeof(vFile->filename)-1] != '\0' ) {
	  printf( "Error: veto file name is too long\n" ); return 1;
	}

	/*-- Set defaults for window --*/
	vFile->winNeg = 0.0;
	vFile->winPos = 0.0;
	vFile->snrRatio = 1000000.0;

	/*-- Copy the window specification string --*/
	strncpy( vetospec, chptr, sizeof(vetospec) );
	if ( vetospec[sizeof(vetospec)-1] != '\0' ) {
	  printf( "Error: veto specification string is too long\n" ); return 1;
	}

	/*-- Loop over tokens in the window specification string, delimited by
	  spaces and/or commas --*/
	chptr = vetospec;
	ipos = 0;
	while ( (tokptr=(char*)strtok(chptr," ,")) != NULL ) {
	  ipos++;
	  switch (ipos) {
	  case 1:
	    sscanf( tokptr, "%lf", &(vFile->winNeg) );
	    break;
	  case 2:
	    sscanf( tokptr, "%lf", &(vFile->winPos) );
	    break;
	  case 3:
	    sscanf( tokptr, "%f", &(vFile->snrRatio) );
	    break;
	  }
	  chptr = NULL;
	}

	nvetofiles++;
	if ( nextptr != NULL ) { arg = nextptr; }

      } while ( nextptr != NULL );

      break;

    case 3:
      /*-- Range or range-file --*/

      /*-- First try to parse as a time range --*/
      chptr = strchr( arg, '-' );
      if ( chptr != NULL ) {
	*chptr = '\0'; chptr++;
	tTemp1 = 0.0; tTemp2 = 2.0e9;
	/*-- Try to parse --*/
	status = 1; status2 = 1;
	if ( strlen(arg) > 0 ) { status = sscanf( arg, "%lf", &tTemp1 ); }
	if ( strlen(chptr) > 0 ) { status2 = sscanf( chptr, "%lf", &tTemp2 ); }
	chptr--; *chptr = '-';
	if ( status && status2 &&
	     (tTemp1 > 1.0e6 || (tTemp2 > 1.0e8 && tTemp2 < 1.9e9) ) ) {
	  /*-- User specified a time range --*/

	  if ( tTemp2 <= tTemp1 ) {
	    printf( "Invalid range specification (t2<=t1): %s\n", arg );
	    return 1;
	  }

	  nRange = 1;
	  InitRange( &range[0], tTemp1, tTemp2 );
	}
      }

      /*-- If nRange is still 0, interpret this as a filename (unless it is
	blank) --*/
      strcpy( temptext, arg );
      if ( nRange == 0 && strtok(temptext," ") != NULL ) {
	fd = fopen( arg, "r" );
	if ( fd == NULL ) {
	  printf( "Error opening range-definition file %s\n", arg );
	  PrintUsage(); return 1;
	}

	/*-- Read lines from the file and parse them --*/
	while ( fgets(intext,sizeof(intext),fd) ) {
	  intext[sizeof(intext)-1] = '\0';

	  /*-- Strip off comments --*/
	  chptr = strchr( intext, '#' );
	  if ( chptr ) { *chptr = '\0'; }

	  /*-- Skip any leading spaces or tabs --*/
	  chptr2 = intext;
	  while ( *chptr2 == ' ' ) { chptr2++; }
	  /*-- If line is blank, go on to the next line --*/
	  if ( strlen(chptr2) == 0 ) { continue; }

	  /*-- Parse as a range --*/

	  chptr = strchr( chptr2, '-' );
	  if ( chptr == NULL ) {
	    /*-- Start and stop times might be separated by a space --*/
	    chptr = strchr( chptr2, ' ' );
	    if ( chptr == NULL ) {
	      printf( "Error parsing time range in %s: %s\n", arg, intext );
	      fclose(fd); return 1;
	    }
	  }
	  *chptr = '\0'; chptr++;
	  tTemp1 = 0.0; tTemp2 = 2.0e9;

	  if ( strlen(chptr2) > 0 ) {
	    if ( sscanf( chptr2, "%lf", &tTemp1 ) < 1 ) {
	      printf( "Error parsing beginning of time range: %s\n", chptr2 );
	      fclose(fd); return 1;
	    }
	  }

	  if ( strlen(chptr) > 0 ) {
	    if ( sscanf( chptr, "%lf", &tTemp2 ) < 1 ) {
	      printf( "Error parsing end of time range: %s\n", chptr );
	      fclose(fd); return 1;
	    }
	  }

	  if ( tTemp2 <= tTemp1 ) {
	    chptr--; *chptr = '-';
	    printf( "Invalid range specification (t2<=t1): %s\n", intext );
	    fclose(fd); return 1;
	  }

	  /*-- Add this to the array of ranges --*/
	  InitRange( &range[nRange], tTemp1, tTemp2 );
	  nRange++;
	}

	fclose( fd );
      }

      break;

    case 4:
      /*-- Output filename --*/

      strncpy( outFile, arg, sizeof(outFile) );
      if ( outFile[sizeof(outFile)-1] != '\0' ) {
	printf( "Error: output file name is too long\n" ); return 1;
      }

      break;

    default:
      printf( "Too many positional arguments\n" );
      PrintUsage(); return 1;

    } /*-- End of switch (iposarg) --*/

  } /*-- End of loop over arguments (iarg) --*/

  if ( iposarg < 2 ) {
    printf( "Not enough arguments\n" );
    PrintUsage(); return 1;
  }

  if (debug >= 1) {
    printf( "nRange=%d\n", nRange );
    for ( jRange=0; jRange<nRange; jRange++ ) {
      printf( "Range %d: %.3lf - %.3lf\n",
	      jRange, range[jRange].t1, range[jRange].t2 );
    }
  }

  /*------ Open the files and read up to the beginning of each Stream ------*/

  /*-- Open file of candidate events --*/
  status = MetaioOpen( candEnv, candFile );
  if ( status != 0 ) {

    /*-- Failed to open this file as a LIGO_LW file --*/
    MetaioAbort( candEnv );
    candEnv->fileRec.fp = NULL;

    /*-- Try opening it as an absGlitch .dat file --*/
    candDatFile = fopen( candFile, "r" );
    if ( candDatFile == NULL ) {
      printf( "Error opening candidate file %s\n", candFile );
      AbortAll( candEnv, vetoFile, nvetofiles, outEnv );
      return 2;
    }

    /*-- Check that this looks like a .dat file --*/
    if ( fgets(intext,sizeof(intext),candDatFile) == NULL ) {
      printf( "Candidate file is empty: %s\n", candFile );
      AbortAll( candEnv, vetoFile, nvetofiles, outEnv );
      return 2;
    }

    if ( strncmp(intext,"### Trigger Statistics:",23) != 0 ) {
      printf( "Candidate file is not a valid LIGO_LW file nor an"
	      " absGlitch .dat file: %s\n", candFile );
      AbortAll( candEnv, vetoFile, nvetofiles, outEnv );
      return 2;
    }
  }

  if ( candDatFile == NULL ) {
    /*-- This is a LIGO_LW file; set up some things --*/

    /*-- Locate the "end_time" and "end_time_ns" columns --*/
    iCandS = MetaioFindColumn( candEnv, "end_time" );
    iCandNS = MetaioFindColumn( candEnv, "end_time_ns" );
    if ( iCandS < 0 || iCandNS < 0 ) {
      printf( "Candidate file %s does not contain end_time or end_time_ns\n",
	      candFile );
      AbortAll( candEnv, vetoFile, nvetofiles, outEnv );
      return 4;
    }

    /*-- It is OK for the "SNR" column to be absent --*/
    iCandSnr = MetaioFindColumn( candEnv, "snr" );

  }

  /*-- Open file(s) of vetoes --*/
  for ( ivfile=0; ivfile<nvetofiles; ivfile++ ) {
    vFile = &(vetoFile[ivfile]);

    /*-- Initialize some structure elements--*/
    vFile->env = &(vFile->parseEnv);

    status = MetaioOpen( vFile->env, vFile->filename );
    if ( status != 0 ) {
      /*-- Failed to open this file as a LIGO_LW file --*/
      MetaioAbort( vFile->env );
      vFile->env->fileRec.fp = NULL;

      /*-- Try opening it as an absGlitch .dat file --*/
      vFile->datFile = fopen( vFile->filename, "r" );
      if ( vFile->datFile == NULL ) {
	printf( "Error opening veto file %s\n", vFile->filename );
	AbortAll( candEnv, vetoFile, nvetofiles, outEnv );
	return 2;
      }

      /*-- Check that this looks like a .dat file --*/
      if ( fgets(intext,sizeof(intext),vFile->datFile) == NULL ) {
	printf( "Veto file is empty: %s\n", vFile->filename );
	AbortAll( candEnv, vetoFile, nvetofiles, outEnv );
	return 2;
      }

      if ( strncmp(intext,"### Trigger Statistics:",23) != 0 ) {
	printf( "Veto file is not a valid LIGO_LW file nor an"
		" absGlitch .dat file: %s\n", vFile->filename );
	AbortAll( candEnv, vetoFile, nvetofiles, outEnv );
	return 2;
      }

      /*-- If we get here, this seems to be a valid absGlitch .dat file --*/
    }

    if ( vFile->datFile == NULL ) {
      /*-- This is a LIGO_LW file; set up some things --*/

      /*-- Locate the relevant columns --*/
      vFile->colS = MetaioFindColumn( vFile->env, "start_time" );
      vFile->colNS = MetaioFindColumn( vFile->env, "start_time_ns" );
      if ( vFile->colS < 0 || vFile->colNS < 0 ) {
	/*-- Try end_time --*/
	vFile->colS = MetaioFindColumn( vFile->env, "end_time" );
	vFile->colNS = MetaioFindColumn( vFile->env, "end_time_ns" );
	if ( vFile->colS < 0 || vFile->colNS < 0 ) {
	  printf( "Veto file %s does not contain (start_time,start_time_ns) or"
		  " (end_time,end_time_ns)\n", vFile->filename );
	  AbortAll( candEnv, vetoFile, nvetofiles, outEnv );
	  return 2;
	}
      }

      /*-- It is OK for the duration column to be absent --*/
      vFile->colDur = MetaioFindColumn( vFile->env, "duration" );

      /*-- It is OK for the SNR column to be absent --*/
      vFile->colSnr = MetaioFindColumn( vFile->env, "snr" );

    }

  } /*-- End loop over veto files --*/

  /*-- Open the output file (if any) --*/
  if ( strlen(outFile) > 0 ) {
    /*-- Make sure input is a LIGO_LW file --*/
    if ( candDatFile ) {
      printf( "Cannot create an output file unless candidate file"
	      " is LIGO_LW\n" );
      AbortAll( candEnv, vetoFile, nvetofiles, outEnv );
      return 2;
    }

    ostatus = MetaioCreate( outEnv, outFile );
    if ( ostatus != 0 ) {
      printf( "Error opening output file %s\n", outFile );
      AbortAll( candEnv, vetoFile, nvetofiles, outEnv );
      return 2;
    }
    MetaioCopyEnv( outEnv, candEnv );
  }


  /*---------- Loop over rows in the candidate file ----------*/

  while ( candeof == 0 ) {

    /*-- Read a row from the candidate file --*/
    if ( candDatFile ) {
      /*-- Reading from an absGlitch .dat file --*/

      while ( 1 ) {

	if ( fgets(intext,sizeof(intext),candDatFile) == NULL ) {
	  /*-- We reached end-of-file --*/
	  candeof = 1;
	  /*-- Set a time extremely far in the future, so that we read all
	    veto events which are "relevant" (i.e. fall within range) --*/
	  tCand = 2.0e9;

	  break;
	}

	/*-- Ignore comment lines --*/
	if ( intext[0] == '#' ) continue;

	/*-- Parse the line; ignore it if unparsable --*/
	if ( sscanf( intext, "%*s %*s %*s %s %lf", ttext, &dur ) < 2 )
	  continue;

	/*-- Parse the time text.  We break it up into integer and
	  fractional parts so that we can subtract the timeBase from the
	  integer part before forming a floating-point number --*/
	chptr = strchr( ttext, '.' );
	if ( chptr == NULL ) {
	  /*-- This is already an integer --*/
	  secfrac = 0.0;
	} else {
	  /*-- Convert fractional part; append 0 in case it's blank --*/
	  strcat( ttext, "0" );
	  secfrac = atof( chptr );
	  /*-- Re-terminate the string to just keep the integer part --*/
	  *chptr = '\0';
	}

	/*-- If this is the first event, set the time base.  This is a trick
	  to retain more numerical precision in time values --*/
	if ( timeBase == 0 ) {
	  timeBase = atoi(ttext);
	  /*-- Modify all time ranges to be relative to the time base --*/
	  for ( jRange=0; jRange<nRange; jRange++ ) {
	    range[jRange].t1 -= (double) timeBase;
	    range[jRange].t2 -= (double) timeBase;
	  }
	}
  
	tCand = (double) ( atoi(ttext) - timeBase ) + secfrac + 0.5*dur;
	/*-- I'm not sure if an "SNR" or "significance" is in the file --*/
	snrCand = 1.0;
	pass = 1;
	if ( debug >= 2 ) printf( "tCand is %.4lf\n", tCand );
	break;

      }  /*-- End of loop over lines in .dat file --*/

    } else {
      /*-- We're reading from a LIGO_LW file --*/

      status = MetaioGetRow(candEnv);
      if ( status == -1 ) {
	printf( "Error while getting row from candidate file %s\n", candFile );
	AbortAll( candEnv, vetoFile, nvetofiles, outEnv );
	return 6;

      } else if ( status == 0 ) {
	/*-- Reached end of file --*/
	candeof = 1;
	/*-- Set a time extremely far in the future, so that we read all
	  veto events which are "relevant" (i.e. fall within range) --*/
	tCand = 2.0e9;

      } else {

	/*-- Get time --*/
	table = &(candEnv->ligo_lw.table);

	/*-- If this is the first event, set the time base.  This is a trick
	  to retain more numerical precision in time values --*/
	if ( timeBase == 0 ) {
	  timeBase = table->elt[iCandS].data.int_4s;
	  /*-- Modify all time ranges to be relative to the time base --*/
	  for ( jRange=0; jRange<nRange; jRange++ ) {
	    range[jRange].t1 -= (double) timeBase;
	    range[jRange].t2 -= (double) timeBase;
	  }
	}

	/*-- Get the time of this event candidate --*/
	tCand = (double) ( table->elt[iCandS].data.int_4s - timeBase )
	  + 1.0e-9 * (double) table->elt[iCandNS].data.int_4s;

	/*-- Also get the SNR (if known) --*/
	if ( iCandSnr >= 0 ) {
	  snrCand = table->elt[iCandSnr].data.real_4;
	} else {
	  snrCand = 1.0;
	}

	pass = 1;
	if ( debug >= 2 ) printf( "tCand is %.4lf\n", tCand );

      }  /*-- End block if we read a row successfully --*/

    }  /*-- End code branch for different kinds of input files --*/


    /*-- If no time range was specified by the user, set up a single range
      beginning at the time of the first candidate event; the end time will be
      set to the time of the last candidate event, later on.  If the user
      specified one or more time ranges but the first time range did not
      include a start time, set it to the time of the first candidate event.
      --*/
    if ( nRange == 0 ) {
      nRange = 1;
      InitRange( &range[0], tCand, 2.0e9 );
    } else if ( range[0].t1 == 0.0 ) {
      range[0].t1 = tCand;
    }

    /*-- Check candidate event time against time range, and increment iRange
      if necessary --*/
    while ( tCand >= cRange->t2 && iRange < nRange-1 ) {
      iRange++;
      cRange = &range[iRange];
      if (debug>=1) printf( "Switching to range %d\n", iRange );
    }
    /*-- Check if candidate event falls between ranges (or beyond end of
      last range) --*/
    if ( tCand < cRange->t1 ) continue;
    /*
    if ( tCand > cRange->t2 ) continue;
    */

    if ( ! candeof ) cRange->nCand++;

    /*-- Discard veto events in ring buffer which are too far in the past to
      be relevant --*/
    while ( vbufR != vbufW && tVetoPos[vbufR] < tCand ) {
      if ( debug >= 2 ) printf( " Discarding an event from the ring buffer\n" );
      vbufR++; if ( vbufR == RINGBUFSIZE ) { vbufR = 0; }
    }

    /*------ Enlarge the veto ring buffer with events from the veto file(s),
      each time choosing the earliest available veto event, until we reach a
      veto event which is too far in the future to be relevant ------*/

    while ( tUseNeg <= tCand && ! allPast ) {

      for ( ivfile=0; ivfile<nvetofiles; ivfile++ ) {
	vFile = &(vetoFile[ivfile]);

	/*-- If flag is set, read a row from this veto file --*/
	if ( vFile->readit && ! vFile->eof ) {

	  vFile->readit = 0;

	  if ( vFile->datFile ) {
	    /*-- Reading from an absGlitch .dat file --*/

	    while ( 1 ) {
	      if ( fgets(intext,sizeof(intext),vFile->datFile) == NULL ) {
		vFile->eof = 1;
		break;
	      }

	      /*-- Ignore comment lines --*/
	      if ( intext[0] == '#' ) continue;

	      /*-- Parse the line; ignore it if unparsable --*/
	      if ( sscanf( intext, "%*s %*s %*s %s %lf", ttext, &dur ) < 2 )
		continue;

	      /*-- Parse the time text.  We break it up into integer and
		fractional parts so that we can subtract the timeBase from the
		integer part before forming a floating-point number --*/
	      chptr = strchr( ttext, '.' );
	      if ( chptr == NULL ) {
		/*-- This is already an integer --*/
		secfrac = 0.0;
	      } else {
		/*-- Convert fractional part; append 0 in case it's blank --*/
		strcat( ttext, "0" );
		secfrac = atof( chptr );
		/*-- Re-terminate the string to just keep the integer part --*/
		*chptr = '\0';
	      }

	      vFile->tLast = (double) ( atoi(ttext) - timeBase ) + secfrac;
	      vFile->tLastNeg = vFile->tLast + vFile->winNeg;
	      vFile->tLastPos = vFile->tLast + dur + vFile->winPos;

	      /*-- I'm not sure whether an SNR or significance is given in a
		.dat file, so just set the SNR to 1 --*/
	      vFile->lastSnr = 1.0;

	      break;
	    }

	  } else {

	    status2 = MetaioGetRow(vFile->env);
	    if ( status2 == -1 ) {
	      printf( "Error while getting row from file %s\n",
		      vFile->filename );
	      AbortAll( candEnv, vetoFile, nvetofiles, outEnv );
	      return 6;
	    } else if ( status2 == 0 ) {
	      /*-- Reached end of file --*/
	      vFile->eof = 1;
	      break;
	    }

	    /*-- Get veto time, etc. --*/
	    table = &(vFile->env->ligo_lw.table);
	    if ( vFile->colDur >= 0 ) {
	      dur = table->elt[vFile->colDur].data.real_8;
	    } else {
	      dur = 0.0;
	    }

	    if ( vFile->colSnr >= 0 ) {
	      vFile->lastSnr = table->elt[vFile->colSnr].data.real_4;
	    } else {
	      vFile->lastSnr = 1.0;
	    }

	    vFile->tLast =
	      (double) ( table->elt[vFile->colS].data.int_4s - timeBase )
	      + 1.0e-9 * (double) table->elt[vFile->colNS].data.int_4s;
	    vFile->tLastNeg = vFile->tLast + vFile->winNeg;
	    vFile->tLastPos = vFile->tLast + dur + vFile->winPos;

	  }  /*-- End code branch for different kinds of input files --*/

	  if ( vFile->eof ) {
	    if ( debug >= 2 ) printf( " vFile %2d is at eof\n", ivfile );
	    continue;
	  }

	  if ( debug >= 3 )
	    printf( " vFile %2d has tLast=%.4lf, duration=%.4lf\n",
		    ivfile, vFile->tLast, dur );

	}  /*-- End block if flag is set to read from this file --*/

      }  /*-- End loop over veto files --*/

      /*---- Now pick the earliest time from among the veto files ----*/
      tUseNeg = 2.0e9;
      usevfile = -1;
      for ( ivfile=0; ivfile<nvetofiles; ivfile++ ) {
	vFile = &(vetoFile[ivfile]);
	if ( vFile->eof ) continue;
	if ( vFile->tLastNeg < tUseNeg ) {
	  tUseNeg = vFile->tLastNeg;
	  /*-- Account for Patrick's "veto clustering" criterion --*/
	  if ( tUseNeg > tUsePos && tUseNeg-tUsePos < MINLIVEINTERVAL ) {
	    /*-- Extend the vetoed interval backward --*/
	    tUseNeg = tUsePos;
	  }
	  tUsePos = vFile->tLastPos;
	  snrThresh = vFile->snrRatio * vFile->lastSnr;
	  usevfile = ivfile;
	}
      }

      /*-- If there are no veto events left, break out of the loop --*/
      if ( usevfile < 0 ) break;

      if ( debug >= 3 ) printf( " Times to use: %.3f, %.3f\n",
				tUseNeg, tUsePos );

      /*-- Mark this veto file to be read again --*/
      vetoFile[usevfile].readit = 1;

      /*-- Check whether veto event is within time range of interest --*/
      while ( tUseNeg >= vRange->t2 && iRangeV < nRange-1 ) {
	iRangeV++;
	vRange = &range[iRangeV];
      }
      /*-- If veto event is fully within a gap between ranges, ignore it --*/
      if ( tUsePos < vRange->t1 ) continue;

      /*-- If beyond end of last range, ignore any remaining events --*/
      if ( iRangeV == nRange-1 && tUseNeg >= vRange->t2 ) {
	allPast = 1;
	break;
      }

      /*-- Now we know this veto event is relevant --*/
      vRange->nVeto++;

      /*-- Figure out what part of veto event is within time range --*/
      tdead1 = ( tUseNeg > vRange->t1 ? tUseNeg : vRange->t1 );
      tdead2 = ( tUsePos < vRange->t2 ? tUsePos : vRange->t2 );

      /*-- Update deadtime total --*/
      if ( tdead1 > tDeadPos ) {
	/*-- There is a live gap between the existing veto cluster and this
	  veto event, so start a new veto cluster --*/
	tDeadNeg = tdead1;
	tDeadPos = tdead2;
	vRange->dead += (tdead2 - tdead1);

      } else if ( tdead2 > tDeadPos ) {
	/*-- This just extends the current veto cluster --*/
	vRange->dead += ( tdead2 - tDeadPos );
	tDeadPos = tdead2;

      } else {
	/*-- This veto event is completely contained in the earlier veto
	  event, so do not change anything --*/
      }

      /*-- If this veto event is too long ago to be relevant, don't add it to
	the ring buffer --*/
      if ( tdead2 < tCand ) continue;

      /*-- Make sure we have more room in the ring buffer for this veto --*/
      if ( vbufW == vbufR-1 || vbufW == vbufR-1+RINGBUFSIZE ) {
	printf( "Ran out of space in ring buffer (size=%d)\n", RINGBUFSIZE );
	AbortAll( candEnv, vetoFile, nvetofiles, outEnv );
	return 5;
      }

      if ( debug >= 2 ) printf( "  adding to ring buffer\n" );

      /*-- Add this veto event to the ring buffer --*/
      tVetoNeg[vbufW] = tUseNeg;
      tVetoPos[vbufW] = tUsePos;
      cSnrThresh[vbufW] = snrThresh;
      rVeto[vbufW] = vRange;
      usedVeto[vbufW] = 0;
      vbufW++; if ( vbufW == RINGBUFSIZE ) { vbufW = 0; }

    }  /*-- End of while loop to fill veto ring buffer --*/

    /*-- If we're at the end of the candidate list, exit the loop --*/
    if ( candeof ) break;

    /*-- See if this is part of the same "cluster" as the last candidate --*/
    if ( tCand-tCandLast < CLUSWINDOW ) {
      /*-- Part of same candidate cluster --*/
    } else {
      /*-- New candidate cluster --*/
      cRange->nClus++;
      clusPass = 0;   /* Will be updated if any event in cluster passes --*/
    }
    tCandLast = tCand;

    /*-- Check for coincidences (within window) with items in ring buffer --*/
    iveto = vbufR;
    while ( iveto != vbufW ) {

      if ( tVetoPos[iveto] < tCand ) {
	/*-- This veto event is no longer relevant; if it is the first item
	  in the ring buffer, drop it --*/
	if ( iveto == vbufR ) {
	  if ( debug >= 2 ) printf( " Dropping an event from the ring buffer\n" );
	  vbufR++; if ( vbufR == RINGBUFSIZE ) { vbufR = 0; }
	}
	
      } else if ( tVetoNeg[iveto] <= tCand ) {
	/*-- This veto event overlaps with the candidate event! --*/

	/*-- Always consider this veto event to be "used", even if the
	  candidate event has such a large SNR that it is not actually vetoed
	  --*/
	if ( usedVeto[iveto] == 0 ) {
	  rVeto[iveto]->nVetoUsed++;
	  usedVeto[iveto] = 1;
	}

	/*-- If candidate SNR is below the threshold, veto the event --*/
	if ( snrCand < cSnrThresh[iveto] ) {
	  pass = 0;
	  if ( debug >= 1 ) printf( "Coinc: cand=%.4lf, veto=%.4lf to %.4lf\n",
				    tCand, tVetoNeg[iveto], tVetoPos[iveto] );
	} else {
	  if ( debug >= 1 ) printf( "Coinc but big SNR: cand=%.4lf (SNR=%.2f),"
				    "veto=%.4lf to %.4lf (SNR thresh=%.2f)\n",
				    tCand, snrCand, tVetoNeg[iveto],
				    tVetoPos[iveto], cSnrThresh[iveto] );
	}

      } else {
	/*-- This veto event is too far in the future, and since the veto
	  events are sorted by time, we can skip all remaining --*/
	break;
      }

      iveto++; if ( iveto == RINGBUFSIZE ) { iveto = 0; }
    }

    /*-- Check whether the candidate event passes --*/
    if (pass) {
      cRange->nCandPass++;
      if ( clusPass == 0 ) {
	cRange->nClusPass++;
	clusPass = 1;
      }

      if ( outEnv->fileRec.fp != NULL ) {
	/*-- Write out the event --*/
	ostatus = MetaioCopyRow( outEnv, candEnv );
	ostatus = MetaioPutRow( outEnv );
      }

    } else {
      cRange->nCandFail++;
    }

  } /*-- End loop over candidate events --*/

  /*------ Close files ------*/

  if ( outEnv->fileRec.fp ) {
    MetaioClose(outEnv);
    outEnv->fileRec.fp = NULL;
  }
  AbortAll( candEnv, vetoFile, nvetofiles, outEnv );

  /*------ Post-calculations ------*/

  /*-- If the upper end of the time range was never specified, set it equal
    to the time of the last candidate event.  Also have to adjust the
    deadtime total if last veto cluster extends beyond end of range. --*/
  if ( cRange->t2 >= 1.9e9 ) {
    cRange->t2 = tCandLast;
    cRange->secs = cRange->t2 - cRange->t1;
    if ( tDeadNeg > cRange->t2 ) {
      vRange->dead -= (tDeadPos - tDeadNeg);
    } else if ( tDeadPos > cRange->t2 ) {
      vRange->dead -= (tDeadPos - cRange->t2 );
    }
  }

  /*-- Fill cluster "fail" fields --*/
  for ( jRange=0; jRange<nRange; jRange++ ) {
    range[jRange].nClusFail = range[jRange].nClus - range[jRange].nClusPass;
  }

  /*-- Initialize an additional range in which to accumulate totals --*/
  InitRange( &range[nRange], 0.0, 0.0 );
  totals = &range[nRange];

  /*-- Accumulate totals --*/
  for ( jRange=0; jRange<nRange; jRange++ ) {
    cRange = &range[jRange];
    totals->secs += cRange->secs;
    totals->dead += cRange->dead;
    totals->nCand += cRange->nCand;
    totals->nCandPass += cRange->nCandPass;
    totals->nCandFail += cRange->nCandFail;
    totals->nClus += cRange->nClus;
    totals->nClusPass += cRange->nClusPass;
    totals->nClusFail += cRange->nClusFail;
    totals->nVeto += cRange->nVeto;
    totals->nVetoUsed += cRange->nVetoUsed;
  }


  /*------ Report results ------*/

  printf( "Candidate clustering window = %.2f seconds\n", CLUSWINDOW );
  printf( "Minimum live interval between vetoes = %.2f seconds\n",
	  MINLIVEINTERVAL );
  for ( ivfile=0; ivfile<nvetofiles; ivfile++ )
    printf( "Veto: %s (%.3f,%.3f,%.3f)\n", vetoFile[ivfile].filename,
	    vetoFile[ivfile].winNeg, vetoFile[ivfile].winPos,
	    vetoFile[ivfile].snrRatio );

/* FORMAT:
___Start__ _Dur_ __Veto__used__used% __Cand___cut___cut% _Clus___cut__cut% dead%
 123456789 12345 123456 123456 12.3% 123456 123456 12.3% 12345 12345 12.3% 12.3%
*/
  printf( "___Start__ _Dur_ __Veto__used__used%% __Cand___cut___cut%%"
	  " _Clus___cut__cut%% dead%%\n" );

  /*-- Loop over all ranges, PLUS the range structure used to hold totals --*/
  for ( jRange=0; jRange<=nRange; jRange++ ) {
    cRange = &range[jRange];

    if ( jRange < nRange ) {
      printf( "%10.0lf %5.0lf", (double)timeBase+cRange->t1, cRange->secs );
    } else {
      printf( "   Total %7.0lf", cRange->secs );
    }

    printf( " %6d %6d", cRange->nVeto, cRange->nVetoUsed );
    if ( cRange->nVetoUsed < cRange->nVeto ) {
      pct = 100.0 * (double) cRange->nVetoUsed / (double) cRange->nVeto;
      if ( pct < 99.95 ) {
	printf( " %4.1f%%", pct );
      } else {
	printf( "  100%%" );
      }
    } else if ( cRange->nVeto == 0 ) {
      printf( "    0 " );
    } else {
      printf( "  100%%" );
    }

    printf( " %6d %6d", cRange->nCand, cRange->nCandFail );
    if ( cRange->nCandFail < cRange->nCand ) {
      pct = 100.0 * (double) cRange->nCandFail / (double) cRange->nCand;
      if ( pct < 99.95 ) {
	printf( " %4.1f%%", pct );
      } else {
	printf( "  100%%" );
      }
    } else if ( cRange->nCand == 0 ) {
      printf( "    0 " );
    } else {
      printf( "  100%%" );
    }

    printf( " %5d %5d", cRange->nClus, cRange->nClusFail );
    if ( cRange->nClusFail < cRange->nClus ) {
      pct = 100.0 * (double) cRange->nClusFail / (double) cRange->nClus;
      if ( pct < 99.95 ) {
	printf( " %4.1f%%", pct );
      } else {
	printf( "  100%%" );
      }
    } else if ( cRange->nClus == 0 ) {
      printf( "    0 " );
    } else {
      printf( "  100%%" );
    }

    pct = 100.0 * cRange->dead / cRange->secs;
    if ( pct < 99.95 ) {
      printf( " %4.1f%%", pct );
    } else {
      printf( "  100%%" );
    }

    printf( "\n" );
  }

  return 0;
}


/*===========================================================================*/
void PrintUsage()
{
  printf( "Usage:  ivana <candidate file> <veto spec> [<range spec>]"
	  " [<output file>]\n" );
  printf( "  <candidate file> is simply the name of a LIGO_LW table file.\n" );
  printf( "  <veto spec> consists of one or more items separated by commas."
                             "  Each item has\n"
          "      the form '<file>(<negWindow>,<posWindow>[,<snrRatio>])'."
                             "  <snrRatio> is\n"
	  "      optional; if specified, then signal candidates will NOT be"
                             " vetoed if the\n"
	  "      signal SNR exceeds the veto SNR by a factor of <snrRatio>"
	                     " or greater.\n"
	  "      NOTE: Enclose <veto spec> in quotes to prevent the shell"
	                     " from\n"
	  "      interpreting the parentheses!\n"
	  "      Examples:  'mich_ctrl5.xml(-0.5,+1,1.5)'\n"
	  "                 'mich_ctrl5.xml(-0.4,0.4),../data/H2_POBQ_snr7.xml"
	                     "(-.3,.3,2)'\n" );
  printf( "  <range spec> can be a GPS time range separated by a dash,\n"
	  "      e.g. '693960000-693961000', or the name of a file containing"
	  " a list of\n      ranges of that form (one per line, optionally"
	  " followed by a comment).\n      If omitted, the range is taken to"
	  " extend from the first candidate event\n      to the last candidate"
	  " event.\n" );
  printf( "  <output file> is the name of the output file to store all"
	  " un-vetoed events in\n      LIGO_LW format.\n" );

  return;
}


/*===========================================================================*/
void InitRange( Range* range, double t1, double t2 )
{
  range->t1 = t1;
  range->t2 = t2;
  range->secs = t2 - t1;
  range->dead = 0.0;
  range->nCand = 0;
  range->nCandPass = 0;
  range->nCandFail = 0;
  range->nClus = 0;
  range->nClusPass = 0;
  range->nClusFail = 0;
  range->nVeto = 0;
  range->nVetoUsed = 0;

  return;
}


/*===========================================================================*/
void AbortAll( MetaioParseEnv candEnv, VetoFile vetoFile[], int nvetofiles,
	       MetaioParseEnv outEnv )
{
  int ivfile;

  if ( candEnv->fileRec.fp ) MetaioAbort( candEnv );

  for ( ivfile=0; ivfile<nvetofiles; ivfile++ ) {
    if ( vetoFile[ivfile].datFile ) {
      fclose( vetoFile[ivfile].datFile );
    } else if ( vetoFile[ivfile].env->fileRec.fp ) {
      MetaioAbort( vetoFile[ivfile].env );
    }
  }

  if ( outEnv->fileRec.fp ) MetaioAbort( outEnv );

  return;
}

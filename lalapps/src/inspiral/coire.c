/*
*  Copyright (C) 2007 Alexander Dietz, Drew Keppel, Gareth Jones, Stephen Fairhurst, Thomas Cokelaer
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

/*----------------------------------------------------------------------- 
 * 
 * File Name: coire.c
 *
 * Author: Fairhurst, S
 * 
 * 
 *-----------------------------------------------------------------------
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <regex.h>
#include <time.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/Date.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataInspiralUtils.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lalapps.h>
#include <processtable.h>
#include <LALAppsVCSInfo.h>

#define PROGRAM_NAME "coire"
#define CVS_ID_STRING "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define CVS_NAME_STRING "$Name$"

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
  this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
calloc( 1, sizeof(ProcessParamsTable) ); \
snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
    PROGRAM_NAME ); \
snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
    long_options[option_index].name ); \
snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );

#define MAX_PATH 4096

/*
 *
 * USAGE
 *
 */


static void print_usage(char *program)
{
  fprintf(stderr,
      "Usage: %s [options] [LIGOLW XML input files]\n"\
      "The following options are recognized.  Options not surrounded in []\n"\
      "are required.\n", program );
  fprintf(stderr,    
      " [--help]                       display this message\n"\
      " [--verbose]                    print progress information\n"\
      " [--version]                    print version information and exit\n"\
      " [--user-tag]          usertag  set the process_params usertag\n"\
      " [--comment]           string   set the process table comment\n"\
      "\n"\
      "  --output             output   write output data to file: output\n"\
      "  --summary-file       summ     write trigger analysis summary to summ\n"\
      "\n"\
      "  --data-type          datatype specify the data type, must be one of\n"\
      "                                (playground_only|exclude_play|all_data)\n"\
      "\n"\
      " [--mass-cut]          masstype keep only triggers in mass range of type\n"\
      "                                (mchirp|mtotal|mcomp)\n"\
      " [--mass-range-low]    lowmass  lower bound on mass range\n"\
      " [--mass-range-high]   highmass upper bound on mass range\n"\
      " [--mass2-range-low]   lowmass  lower bound on mass2 range for mcomp\n"\
      " [--mass2-range-high]  highmass upper bound on mass2 range for mcomp\n"\
      "\n"\
      " [--discard-ifo]       ifo      discard all triggers from ifo\n"\
      " [--coinc-cut]         ifos     only keep triggers from IFOS\n"\
      " [--extract-slide]     slide    only keep triggers from specified slide\n"\
      "\n"\
      " [--num-slides]        slides   number of time slides performed \n"\
      "                                (used in clustering)\n"\
      " [--sort-triggers]              time sort the coincident triggers\n"\
      " [--coinc-stat]        stat     use coinc statistic for cluster/cut\n"\
      "                     [ snrsq | effective_snrsq | s3_snr_chi_stat | bitten_l | new_snrsq]\n"\
      " [--stat-threshold]    thresh   discard all triggers with stat less than thresh\n"\
      " [--rsq-threshold] rsq_thresh   discard all triggers whose rsqveto_duration\n"\
      "                                exceeds rsq_thresh\n"\
      " [--rsq-max-snr]   rsq_max_snr  apply rsq on triggers with snr < rsq_max_snr\n"\
      "                                exceeds rsq_thresh\n"\
      " [--rsq-coeff]     rsq_coeff    apply rsq on triggers with snr > rsq_max_snr\n"\
      "                                exceeds rsq_coeff * snr ^ rsq_power\n"\
      " [--rsq-power]     rsq_power    apply rsq on triggers with snr > rsq_max_snr\n"\
      "                                exceeds rsq_coeff * snr ^ rsq_power\n"\
      " [--eff-snr-denom-fac] number   parameter for clustering effective snr denominator (traditionally 250) \n"\
      " [--chisq-index]       number   parameter for clustering in new-snrsq (traditionally 6) \n"\
      " [--h1-bittenl-a]      bitten   paramater a for clustering\n"\
      " [--h1-bittenl-b]      bitten   paramater b for clustering\n"\
      " [--h2-bittenl-a]      bitten   paramater a for clustering\n"\
      " [--h2-bittenl-b]      bitten   paramater b for clustering\n"\
      " [--l1-bittenl-a]      bitten   paramater a for clustering\n"\
      " [--l1-bittenl-b]      bitten   paramater b for clustering\n"\
      " [--cluster-time]      time     cluster triggers with time ms window\n"\
      "\n"\
      " [--injection-file]    inj_file read injection parameters from inj_file\n"\
      " [--injection-window]  inj_win  trigger and injection coincidence window (ms)\n"\
      " [--missed-injections] missed   write missed injections to file missed\n"\
      "\n");
    fprintf( stderr, "\n");
    fprintf( stderr, "[LIGOLW XML input files] list of the input trigger files.\n");
}

int sortTriggers = 0;
LALPlaygroundDataMask dataType;
extern int vrbflg;

int main( int argc, char *argv[] )
{
  /* lal initialization variables */
  LALStatus status = blank_status ;

  /*  program option variables */
  CHAR *userTag = NULL;
  CHAR comment[LIGOMETA_COMMENT_MAX];
  char *ifos = NULL;
  char *ifo  = NULL;
  char *outputFileName = NULL;
  char *summFileName = NULL;
  CoincInspiralStatistic coincstat = no_stat;
  REAL4 statThreshold = 0;
  INT8 cluster_dt = 0;
  char *injectFileName = NULL;
  INT8 injectWindowNS = -1;
  char *missedFileName = NULL;
  int j;
  FILE *fp = NULL;
  int numInFiles = 0;

  char *massCut = NULL;
  REAL4 massRangeLow = -1;
  REAL4 massRangeHigh = -1;
  REAL4 mass2RangeLow = -1;
  REAL4 mass2RangeHigh = -1;
  REAL4 rsqVetoThresh = -1;
  REAL4 rsqMaxSnr     = -1;
  REAL4 rsqAboveSnrCoeff = -1;
  REAL4 rsqAboveSnrPow     = -1;

  UINT8 triggerInputTimeNS = 0;

  MetadataTable         proctable;
  MetadataTable         procparams;
  ProcessParamsTable   *this_proc_param;

  int                   numSimEvents = 0;
  int                   numSimInData = 0;

  SearchSummvarsTable  *inputFiles = NULL;
  SearchSummvarsTable  *thisInputFile = NULL;
  SearchSummaryTable   *searchSummList = NULL;
  SearchSummaryTable   *thisSearchSumm = NULL;
  SummValueTable       *summValueList = NULL;

  SimInspiralTable     *simEventHead = NULL;
  SimInspiralTable     *thisSimEvent = NULL;
  SimInspiralTable     *missedSimHead = NULL;
  SimInspiralTable     *missedSimCoincHead = NULL;
  SimInspiralTable     *tmpSimEvent = NULL;

  int                   extractSlide = 0;
  int                   numSlides = 0;
  int                   numTriggers = 0;
  int                   numCoincs = 0;
  int                   numEventsInIfos = 0;
  int                   numEventsBelowRsqThresh = 0;
  int                   numEventsPlayTest = 0;
  int                   numEventsInMassRange = 0;
  int                   numEventsAboveThresh = 0;
  int                   numEventsCoinc = 0;
  int                   numClusteredEvents = 0;
  int                   numSnglFound = 0;
  int                   numCoincFound = 0;

  SnglInspiralTable    *inspiralEventList = NULL;
  SnglInspiralTable    *thisSngl = NULL;
  SnglInspiralTable    *missedSnglHead = NULL;
  SnglInspiralTable    *thisInspiralTrigger = NULL;
  SnglInspiralTable    *snglOutput = NULL;

  CoincInspiralTable   *coincHead = NULL;
  CoincInspiralTable   *thisCoinc = NULL;
  CoincInspiralTable   *missedCoincHead = NULL;

  CoincInspiralStatParams    bittenLParams;

  LIGOLwXMLStream       xmlStream;
  MetadataTable         outputTable;
  MetadataTable         searchSummvarsTable;


  /*
   *
   * initialization
   *
   */

  /* set up inital debugging values */
  lal_errhandler = LAL_ERR_EXIT;

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) 
    calloc( 1, sizeof(ProcessTable) );
  XLALGPSTimeNow(&(proctable.processTable->start_time));
  XLALPopulateProcessTable(proctable.processTable, PROGRAM_NAME, lalAppsVCSIdentId,
      lalAppsVCSIdentStatus, lalAppsVCSIdentDate, 0);
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *) 
    calloc( 1, sizeof(ProcessParamsTable) );
  memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );

  memset( &bittenLParams, 0, sizeof(CoincInspiralStatParams   ) );
  /* default value from traditional effective snr formula */
  bittenLParams.eff_snr_denom_fac = 250.0; 
  bittenLParams.chisq_index = 6.0;

  /*
   *
   * parse command line arguments
   *
   */


  while (1)
  {
    /* getopt arguments */
    static struct option long_options[] = 
    {
      {"verbose",                 no_argument,       &vrbflg,              1 },
      {"sort-triggers",           no_argument,  &sortTriggers,             1 },
      {"help",                    no_argument,            0,              'h'},
      {"user-tag",                required_argument,      0,              'Z'},
      {"userTag",                 required_argument,      0,              'Z'},
      {"comment",                 required_argument,      0,              'c'},
      {"version",                 no_argument,            0,              'V'},
      {"data-type",               required_argument,      0,              'k'},
      {"output",                  required_argument,      0,              'o'},
      {"summary-file",            required_argument,      0,              'S'},
      {"extract-slide",           required_argument,      0,              'e'},
      {"num-slides",              required_argument,      0,              'N'},
      {"coinc-stat",              required_argument,      0,              'C'},
      {"stat-threshold",          required_argument,      0,              'E'},
      {"rsq-threshold",           required_argument,      0,              'x'},
      {"rsq-max-snr",             required_argument,      0,              'U'},
      {"rsq-coeff",               required_argument,      0,              'X'},
      {"rsq-power",               required_argument,      0,              'P'},
      {"cluster-time",            required_argument,      0,              't'},
      {"discard-ifo",             required_argument,      0,              'd'},
      {"coinc-cut",               required_argument,      0,              'D'},
      {"injection-file",          required_argument,      0,              'I'},
      {"injection-window",        required_argument,      0,              'T'},
      {"eff-snr-denom-fac",    required_argument,      0,              'A'},
      {"chisq-index",             required_argument,      0,              's'},
      {"h1-bittenl-a",            required_argument,      0,              'a'},
      {"h1-bittenl-b",            required_argument,      0,              'b'},
      {"h2-bittenl-a",            required_argument,      0,              'j'},
      {"h2-bittenl-b",            required_argument,      0,              'n'},
      {"l1-bittenl-a",            required_argument,      0,              'l'},
      {"l1-bittenl-b",            required_argument,      0,              'p'},
      {"missed-injections",       required_argument,      0,              'm'},
      {"mass-cut",                required_argument,      0,              'M'},
      {"mass-range-low",          required_argument,      0,              'q'},
      {"mass-range-high",         required_argument,      0,              'Q'},
      {"mass2-range-low",         required_argument,      0,              'r'},
      {"mass2-range-high",        required_argument,      0,              'R'},

      {0, 0, 0, 0}
    };
    int c;

    /* getopt_long stores the option index here. */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long_only ( argc, argv, "A:a:b:c:d:hj:k:l:m:n:o:p:q:r:s:t:x:"
                                       "C:D:E:I:M:N:P:Q:R:S:T:U:VZ", 
                                       long_options, 
                                       &option_index );

    /* detect the end of the options */
    if ( c == - 1 )
      break;

    switch ( c )
    {
      case 0:
        /* if this option set a flag, do nothing else now */
        if ( long_options[option_index].flag != 0 )
        {
          break;
        }
        else
        {
          fprintf( stderr, "error parsing option %s with argument %s\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        break;

      case 'A':
        bittenLParams.eff_snr_denom_fac = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg);
        break;

      case 'a':
        bittenLParams.param_a[LAL_IFO_H1] = atof(optarg);
      ADD_PROCESS_PARAM( "float", "%s", optarg);
        break;
      
      case 'b':
        bittenLParams.param_b[LAL_IFO_H1] = atof(optarg);
      ADD_PROCESS_PARAM( "float", "%s", optarg);
        break;

      case 'j':
        bittenLParams.param_a[LAL_IFO_H2] = atof(optarg);
      ADD_PROCESS_PARAM( "float", "%s", optarg);
        break;

      case 'n':
        bittenLParams.param_b[LAL_IFO_H2] = atof(optarg);
      ADD_PROCESS_PARAM( "float", "%s", optarg);
        break;

      case 'l':
        bittenLParams.param_a[LAL_IFO_L1] = atof(optarg);
      ADD_PROCESS_PARAM( "float", "%s", optarg);
        break;

      case 'p':
        bittenLParams.param_b[LAL_IFO_L1] = atof(optarg);
      ADD_PROCESS_PARAM( "float", "%s", optarg);
        break;

      case 's':
        bittenLParams.chisq_index = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg);
        break;

      case 'h':
        print_usage(argv[0]);
        exit( 0 );
        break;

      case 'Z':
        /* create storage for the usertag */
        optarg_len = strlen( optarg ) + 1;
        userTag = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( userTag, optarg, optarg_len );

        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", 
            PROGRAM_NAME );
        snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "-userTag" );
        snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%s",
            optarg );
        break;

      case 'c':
        if ( strlen( optarg ) > LIGOMETA_COMMENT_MAX - 1 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "comment must be less than %d characters\n",
              long_options[option_index].name, LIGOMETA_COMMENT_MAX );
          exit( 1 );
        }
        else
        {
          snprintf( comment, LIGOMETA_COMMENT_MAX, "%s", optarg);
        }
        break;

      case 'V':
        fprintf( stdout, "Coincident Inspiral Reader and Injection Analysis\n"
            "Steve Fairhurst\n");
        XLALOutputVersionString(stderr, 0);
        exit( 0 );
        break;

      case 'o':
        /* create storage for the output file name */
        optarg_len = strlen( optarg ) + 1;
        outputFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( outputFileName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'e':
        /* store the number of slides */
        extractSlide = atoi( optarg );
        if ( extractSlide == 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "extractSlide must be non-zero: "
              "(%d specified)\n",
              long_options[option_index].name, extractSlide );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", extractSlide );
        break;
        
      case 'N':
        /* store the number of slides */
        numSlides = atoi( optarg );
        if ( numSlides < 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "numSlides >= 0: "
              "(%d specified)\n",
              long_options[option_index].name, numSlides );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", numSlides );
        break;
        
      case 'S':
        /* create storage for the summ file name */
        optarg_len = strlen( optarg ) + 1;
        summFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( summFileName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'k':
        /* type of data to analyze */
        if ( ! strcmp( "playground_only", optarg ) )
        {
          dataType = playground_only;
        }
        else if ( ! strcmp( "exclude_play", optarg ) )
        {
          dataType = exclude_play;
        }
        else if ( ! strcmp( "all_data", optarg ) )
        {
          dataType = all_data;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown data type, %s, specified: "
              "(must be playground_only, exclude_play or all_data)\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'C':
        /* choose the coinc statistic */
        {        
          if ( ! strcmp( "snrsq", optarg ) )
          {
            coincstat = snrsq;
          }
          else if ( ! strcmp( "bitten_l", optarg ) )
          {
            coincstat = bitten_l;
          }
          else if ( ! strcmp( "bitten_lsq", optarg ) )
          {
            coincstat = bitten_lsq;
          }

          else if ( ! strcmp( "s3_snr_chi_stat", optarg) )
          {
            coincstat = s3_snr_chi_stat;
          }
          else if ( ! strcmp( "effective_snrsq", optarg) )
          {
            coincstat = effective_snrsq;
          }
          else if ( ! strcmp( "new_snrsq", optarg) )
          {
            coincstat = new_snrsq;
          }
          else
          {
            fprintf( stderr, "invalid argument to  --%s:\n"
                "unknown coinc statistic:\n "
                "%s (must be one of:\n"
                "snrsq, effective_snrsq, bitten_l, s3_snr_chi_stat, new_snrsq)\n",
                long_options[option_index].name, optarg);
            exit( 1 );
          }
          ADD_PROCESS_PARAM( "string", "%s", optarg );
        }
        break;

      case 'E':
        /* store the stat threshold for a cut */
        statThreshold = atof( optarg );
        if ( statThreshold <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "statThreshold must be positive: (%f specified)\n",
              long_options[option_index].name, statThreshold );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%f", statThreshold );
        break;
        
      case 'x':
        rsqVetoThresh = (REAL4) atof( optarg );
        if ( rsqVetoThresh < 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "threshold must be >= 0: "
              "(%f specified)\n",
              long_options[option_index].name, rsqVetoThresh );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", rsqVetoThresh );
        break;

      case 'U':
        rsqMaxSnr = (REAL4) atof( optarg );
        if ( rsqMaxSnr < 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "threshold must be >= 0: "
              "(%f specified)\n",
              long_options[option_index].name, rsqMaxSnr );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", rsqMaxSnr );
        break;

      case 'X':
        rsqAboveSnrCoeff = (REAL4) atof( optarg );
        if ( rsqAboveSnrCoeff < 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "coefficient must be >= 0: "
              "(%f specified)\n",
              long_options[option_index].name, rsqAboveSnrCoeff );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", rsqAboveSnrCoeff );
        break;

      case 'P':
        rsqAboveSnrPow = (REAL4) atof( optarg );
        if ( rsqAboveSnrPow < 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "power must be >= 0: "
              "(%f specified)\n",
              long_options[option_index].name, rsqAboveSnrPow );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", rsqAboveSnrPow );
        break;


      case 't':
        /* cluster time is specified on command line in ms */
        cluster_dt = (INT8) atoi( optarg );
        if ( cluster_dt <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "custer window must be > 0: "
              "(%" LAL_INT8_FORMAT " specified)\n",
              long_options[option_index].name, cluster_dt );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%" LAL_INT8_FORMAT, cluster_dt );
        /* convert cluster time from ms to ns */
        cluster_dt *= 1000000LL;
        break;

      case 'I':
        /* create storage for the injection file name */
        optarg_len = strlen( optarg ) + 1;
        injectFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( injectFileName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'd':
        /* discard all triggers from ifo */
        optarg_len = strlen( optarg ) + 1;
        ifo = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( ifo, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'D':
        /* keep only coincs found in ifos */
        optarg_len = strlen( optarg ) + 1;
        ifos = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( ifos, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'T':
        /* injection coincidence time is specified on command line in ms */
        injectWindowNS = (INT8) atoi( optarg );
        if ( injectWindowNS < 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "injection coincidence window must be >= 0: "
              "(%" LAL_INT8_FORMAT " specified)\n",
              long_options[option_index].name, injectWindowNS );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%" LAL_INT8_FORMAT, injectWindowNS );
        /* convert inject time from ms to ns */
        injectWindowNS *= LAL_INT8_C(1000000);
        break;

      case 'm':
        /* create storage for the missed injection file name */
        optarg_len = strlen( optarg ) + 1;
        missedFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( missedFileName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'M':
        /* create storage for the missed injection file name */
        optarg_len = strlen( optarg ) + 1;
        massCut = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( massCut, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'q':
        massRangeLow = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg);
        break;

      case 'Q':
        massRangeHigh = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg);
        break;

      case 'r':
        mass2RangeLow = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg);
        break;

      case 'R':
        mass2RangeHigh = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg);
        break;

      case '?':
        exit( 1 );
        break;

      default:
        fprintf( stderr, "unknown error while parsing options\n" );
        exit( 1 );
    }   
  }


  /*
   *
   * can use LALCalloc() / LALMalloc() from here
   *
   */


  /* don't buffer stdout if we are in verbose mode */
  if ( vrbflg ) setvbuf( stdout, NULL, _IONBF, 0 );

  /* fill the comment, if a user has specified it, or leave it blank */
  if ( ! *comment )
  {
    snprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX, " " );
  }
  else
  {
    snprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX,
        "%s", comment );
  }

  /* check that the output file name has been specified */
  if ( ! outputFileName )
  {
    fprintf( stderr, "--output must be specified\n" );
    exit( 1 );
  }

  /* check that Data Type has been specified */
  if ( dataType == unspecified_data_type )
  {
    fprintf( stderr, "Error: --data-type must be specified\n");
    exit(1);
  }


  /* check that if clustering is being done that we have all the options */
  if ( cluster_dt && (coincstat == no_stat) )
  {
    fprintf( stderr, 
        "--coinc-stat must be specified if --cluster-time is given\n" );
    exit( 1 );
  }

  /* check that if the rsq veto is being preformed,
                         we have the required options */
  if ( ( (rsqVetoThresh > 0) || (rsqMaxSnr > 0) ) && ( (rsqVetoThresh < 0)
    || (rsqMaxSnr < 0) ) )
  {
    fprintf( stderr, "--rsq-threshold and --rsq-max-snr and must be "
      "specified together" );
    exit( 1 );
  }
  else if ( (rsqAboveSnrCoeff > 0) && ( (rsqMaxSnr < 0) || (rsqVetoThresh < 0)
    || (rsqAboveSnrPow < 0) ) )
  {
    fprintf( stderr, "--rsq-max-snr --rsq-threshold and --rsq-power "
      "must be specified if --rsq-coeff is given\n" );
    exit( 1 );
  }
  else if ( (rsqAboveSnrPow > 0) && ( (rsqMaxSnr < 0) || (rsqVetoThresh < 0)
    || (rsqAboveSnrCoeff < 0) ) )
  {
    fprintf( stderr, "--rsq-max-snr --rsq-threshold and --rsq-coeff "
      "must be specified if --rsq-power is given\n" );
    exit( 1 );
  }

  /* check that if stat cut is being done that we have all the options */
  if ( statThreshold && (coincstat == no_stat) )
  {
    fprintf( stderr, 
        "--coinc-stat must be specified if --stat-threshold is given\n" );
    exit( 1 );
  }
  
  /* check that we have all the options to do injections */
  if ( injectFileName && injectWindowNS < 0 )
  {
    fprintf( stderr, "--injection-window must be specified if "
        "--injection-file is given\n" );
    exit( 1 );
  }
  else if ( ! injectFileName && injectWindowNS >= 0 )
  {
    fprintf( stderr, "--injection-file must be specified if "
        "--injection-window is given\n" );
    exit( 1 );
  }

  if ( numSlides && extractSlide )
  {
    fprintf( stderr, "--num-slides and --extract-slide both specified\n"
        "this doesn't make sense\n" );
    exit( 1 );
  }
 
  if ( ( massCut || massRangeLow >= 0 || massRangeHigh >= 0 ) && 
       ! ( massCut && massRangeLow >= 0 && massRangeHigh >= 0 ) )
  {
    fprintf( stderr, "--mass-cut, --mass-range-low, and --mass-rang-high "
        "must all be used together\n" );
    exit( 1 );
  }

  if ( massCut && ( massRangeLow >= massRangeHigh ) )
  {
    fprintf( stderr, "--mass-range-low must be less than "
        "--mass-range-high\n" );
    exit( 1 );
  }

  if ( massCut && ( ! strcmp( "mcomp", massCut ) ) &&
       ( mass2RangeLow < 0 || mass2RangeHigh <0 ) )
  {
    fprintf( stderr, "--mass2-range-low and --mass2-rang-high \n"
        "must be specified if using --mass-cut mcomp\n" );
    exit( 1 );
  }

  if ( massCut && ( ! strcmp( "mcomp", massCut ) ) &&
       mass2RangeLow >= mass2RangeHigh )
  {
    fprintf( stderr, "--mass2-range-low must be less than "
        "--mass2-range-high\n" );
    exit( 1 );
  }

  if ( massCut && strcmp( "mchirp", massCut ) &&
       strcmp("mtotal", massCut) && strcmp( "mcomp", massCut ) )
  {
    fprintf( stderr, "--mass-cut must be either mchirp, mtotal, or mcomp\n" );
    exit( 1 );
  }
 
  /* save the sort triggers flag */
  if ( sortTriggers )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *) 
      calloc( 1, sizeof(ProcessParamsTable) ); 
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s",
        PROGRAM_NAME ); 
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--sort-triggers" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" ); 
    snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, " " );
  }


  /*
   *
   * read in the input triggers from the xml files
   *
   */


  /* if we have run out of arguments on the command line, throw an error */
  if ( ! (optind < argc) )
  {
    fprintf( stderr, "Error: No input trigger files specified.\n" );
    exit( 1 );
  }

  /* read in the triggers */
  for( j = optind; j < argc; ++j )
  {
    INT4 numFileTriggers = 0;
    INT4 numFileCoincs   = 0;
    SnglInspiralTable   *inspiralFileList = NULL;
    SnglInspiralTable   *thisFileTrigger  = NULL;
    CoincInspiralTable  *coincFileHead    = NULL;

    numInFiles++;
    
    numFileTriggers = XLALReadInspiralTriggerFile( &inspiralFileList,
        &thisFileTrigger, &searchSummList, &inputFiles, argv[j] );
    if (numFileTriggers < 0)
    {
      fprintf(stderr, "Error reading triggers from file %s\n",
          argv[j]);
      exit( 1 );
    }
    else
    {
      if ( vrbflg )
      {
        fprintf(stdout, "Read %d reading triggers from file %s\n",
            numFileTriggers, argv[j]);
      }
    }

    /* read the summ value table as well. */
    XLALReadSummValueFile(&summValueList, argv[j]);

    /* Discard triggers from requested ifo */
    if ( ifo )
    {
      SnglInspiralTable *ifoTrigList = NULL;

      ifoTrigList = XLALIfoCutSingleInspiral( &inspiralFileList, ifo );

      /* discard events from ifo */
      while ( ifoTrigList )
      {
        thisSngl = ifoTrigList;
        ifoTrigList = ifoTrigList->next;
        XLALFreeSnglInspiral( &thisSngl );
      }
      numFileTriggers = XLALCountSnglInspiral( inspiralFileList );
      if ( vrbflg ) fprintf( stdout, 
          "Have %d triggers after discarding those from ifo %s\n",
          numFileTriggers, ifo );
    }

    /*  Do rsq cut */
    if ( rsqVetoThresh > 0 )
    {
      inspiralFileList = XLALRsqCutSingleInspiral( inspiralFileList, 
          rsqVetoThresh, rsqMaxSnr, rsqAboveSnrCoeff, rsqAboveSnrPow );
      /* count the triggers  */
      numFileTriggers = XLALCountSnglInspiral( inspiralFileList );

      if ( vrbflg ) fprintf( stdout, "Have %d triggers after rsq cut\n",
          numFileTriggers );
      numEventsBelowRsqThresh += numFileTriggers;
    }

    /* reconstruct the coincs */
    numFileCoincs = XLALRecreateCoincFromSngls( &coincFileHead, 
        &inspiralFileList );
    if( numFileCoincs < 0 )
    {
      fprintf(stderr, 
          "Unable to reconstruct coincs from single ifo triggers");
      exit( 1 );
    }
   
    /* If there are any remaining triggers ... */
    if ( inspiralFileList )
    {
      /* add inspirals to list */
      if ( thisInspiralTrigger )
      {
        thisInspiralTrigger->next = inspiralFileList;
      }
      else
      {
        inspiralEventList = thisInspiralTrigger = inspiralFileList;
      }
      for( ; thisInspiralTrigger->next; 
          thisInspiralTrigger = thisInspiralTrigger->next);
      numTriggers += numFileTriggers;
    }
 
    if ( vrbflg )
    {
      fprintf( stdout,
          "Recreated %d coincs from the %d triggers in file %s\n", 
          numFileCoincs, numFileTriggers, argv[j] );
    }
    numCoincs += numFileCoincs;

    /* keep only the requested coincs */
    if( ifos )
    {
      numFileCoincs = XLALCoincInspiralIfosCut( &coincFileHead, ifos );
      if ( vrbflg ) fprintf( stdout,
          "Kept %d coincs from %s instruments\n", numFileCoincs, ifos );
      numEventsInIfos += numFileCoincs;
    }

    /* Do playground_only or exclude_play cut */
    if ( dataType != all_data )
    {
      coincFileHead = XLALPlayTestCoincInspiral( coincFileHead, 
          &dataType );
      /* count the triggers, scroll to end of list */
      numFileCoincs = XLALCountCoincInspiral( coincFileHead );

      if ( dataType == playground_only && vrbflg ) fprintf( stdout, 
        "Have %d playground triggers\n", numFileCoincs );
      else if ( dataType == exclude_play && vrbflg ) fprintf( stdout, 
        "Have %d non-playground triggers\n", numFileCoincs );
      numEventsPlayTest += numFileCoincs;
    }

    /* Do mean mass cut */
    if ( massCut )
    {
      numFileCoincs = XLALCountCoincInspiral( coincFileHead );

      coincFileHead = XLALMeanMassCut( coincFileHead, massCut,
          massRangeLow, massRangeHigh, mass2RangeLow, mass2RangeHigh);
      /* count the triggers, scroll to end of list */
      numFileCoincs = XLALCountCoincInspiral( coincFileHead );

      if ( vrbflg ) fprintf( stdout,
          "Kept %d coincs in mass range %f to %f\n", numFileCoincs, 
            massRangeLow, massRangeHigh );
      numEventsInMassRange += numFileCoincs;
    }


    /* perform the statistic cut */
    if( statThreshold )
    {
      coincFileHead = XLALStatCutCoincInspiral ( coincFileHead, coincstat , 
          &bittenLParams, statThreshold);
      numFileCoincs = XLALCountCoincInspiral( coincFileHead );
      if ( vrbflg ) fprintf( stdout,
          "Kept %d coincs above threshold of %6.2f\n", numFileCoincs, 
          statThreshold );
      numEventsAboveThresh += numFileCoincs;
    }
    
    /* add coincs to list */
    if( numFileCoincs )
    {
      if ( thisCoinc )
      {
        thisCoinc->next = coincFileHead;
      }
      else
      {
        coincHead = thisCoinc = coincFileHead;
      }
      for ( ; thisCoinc->next; thisCoinc = thisCoinc->next );
    }
    
  }
        
  if ( vrbflg )
  {
    fprintf( stdout, "Read in %d triggers\n", numTriggers );
    fprintf( stdout, "Recreated %d coincs\n", numCoincs );
    if ( ifos )
    {
      fprintf( stdout, "Have %d coincs from %s\n", numEventsInIfos,
          ifos );
    }
    if ( dataType != all_data ) 
    {
      fprintf( stdout, 
          "Have %d coincs after play test\n", numEventsPlayTest);
    }
    if ( statThreshold ) 
    {
      fprintf( stdout, 
          "Have %d coincs above threshold of %6.2f\n", numEventsAboveThresh, 
          statThreshold);
    }
  }

  for ( thisSearchSumm = searchSummList; thisSearchSumm; 
      thisSearchSumm = thisSearchSumm->next )
  {
    UINT8 outPlayNS, outStartNS, outEndNS, triggerTimeNS;
    LIGOTimeGPS inPlay, outPlay;
    outStartNS = XLALGPSToINT8NS( &(thisSearchSumm->out_start_time) );
    outEndNS = XLALGPSToINT8NS( &(thisSearchSumm->out_end_time) );
    triggerTimeNS = outEndNS - outStartNS;

    /* check for events and playground */
    if ( dataType != all_data )
    {
      XLALPlaygroundInSearchSummary( thisSearchSumm, &inPlay, &outPlay );
      outPlayNS = XLALGPSToINT8NS( &outPlay );

      if ( dataType == playground_only )
      {
        /* increment the total trigger time by the amount of playground */
        triggerInputTimeNS += outPlayNS;
      }
      else if ( dataType == exclude_play )
      {
        /* increment the total trigger time by the out time minus */
        /* the time that is in the playground                     */
        triggerInputTimeNS += triggerTimeNS - outPlayNS;
      }
    }
    else
    {
      /* increment the total trigger time by the out time minus */
      triggerInputTimeNS += triggerTimeNS;
    }
  }




  /*
   *
   * sort the inspiral events by time
   *
   */


  if ( sortTriggers || cluster_dt )
  {
    if ( vrbflg ) fprintf( stdout, "sorting coinc inspiral trigger list..." );
    coincHead = XLALSortCoincInspiral( coincHead, 
        *XLALCompareCoincInspiralByTime );
    if ( vrbflg ) fprintf( stdout, "done\n" );
  }


  /*
   *
   * read in the injection XML file, if we are doing an injection analysis
   *
   */

  if ( injectFileName )
  {
    if ( vrbflg ) 
      fprintf( stdout, "reading injections from %s... ", injectFileName );

    numSimEvents = SimInspiralTableFromLIGOLw( &simEventHead, 
        injectFileName, 0, 0 );

    if ( vrbflg ) fprintf( stdout, "got %d injections\n", numSimEvents );

    if ( numSimEvents < 0 )
    {
      fprintf( stderr, "error: unable to read sim_inspiral table from %s\n", 
          injectFileName );
      exit( 1 );
    }


    /* keep play/non-play/all injections */
    if ( dataType == playground_only && vrbflg ) fprintf( stdout, 
        "Keeping only playground injections\n" );
    else if ( dataType == exclude_play && vrbflg ) fprintf( stdout, 
        "Keeping only non-playground injections\n" );
    else if ( dataType == all_data && vrbflg ) fprintf( stdout, 
        "Keeping all injections\n" );
    XLALPlayTestSimInspiral( &simEventHead, &dataType );

    /* keep only injections in times analyzed */
    numSimInData = XLALSimInspiralInSearchedData( &simEventHead, 
        &searchSummList ); 

    if ( vrbflg ) fprintf( stdout, "%d injections in analyzed data\n", 
        numSimInData );


    /* check for events that are coincident with injections */

    
    if ( vrbflg ) fprintf( stdout, 
        "Sorting single inspiral triggers before injection coinc test\n" );
    inspiralEventList = XLALSortSnglInspiral( inspiralEventList, 
        *LALCompareSnglInspiralByTime );
    
    /* first find singles which are coincident with injections */
    numSnglFound = XLALSnglSimInspiralTest( &simEventHead, 
        &inspiralEventList, &missedSimHead, &missedSnglHead, injectWindowNS );

    if ( vrbflg ) fprintf( stdout, "%d injections found in single ifo\n", 
        numSnglFound );

    /* then check for coincs coincident with injections */
    numCoincFound = XLALCoincSimInspiralTest ( &simEventHead,  &coincHead, 
        &missedSimCoincHead, &missedCoincHead );

    if ( vrbflg ) fprintf( stdout, "%d injections found in coincidence\n", 
        numCoincFound );

    if ( numCoincFound )
    {
      for ( thisCoinc = coincHead; thisCoinc; thisCoinc = thisCoinc->next,
          numEventsCoinc++ );
      if ( vrbflg ) fprintf( stdout, "%d coincs found at times of injection\n",
          numEventsCoinc );
    }
    
    if ( missedSimCoincHead )
    {
      /* add these to the list of missed Sim's */
      if ( missedSimHead )
      {
        for (thisSimEvent = missedSimHead; thisSimEvent->next; 
            thisSimEvent = thisSimEvent->next );
        thisSimEvent->next = missedSimCoincHead;
      }
      else
      {
        missedSimHead = missedSimCoincHead;
      }
    }

    /* free the missed singles and coincs */
    while ( missedCoincHead )
    {
      thisCoinc = missedCoincHead;
      missedCoincHead = missedCoincHead->next;
      XLALFreeCoincInspiral( &thisCoinc );
    }

    while ( missedSnglHead )
    {
      thisSngl = missedSnglHead;
      missedSnglHead = missedSnglHead->next;
      XLALFreeSnglInspiral( &thisSngl );
    }

  } 


  /*
   *
   * extract specified slide
   *
   */

  if ( extractSlide )
  {
    CoincInspiralTable *slideCoinc = NULL;
    slideCoinc = XLALCoincInspiralSlideCut( &coincHead, extractSlide );
    /* free events from other slides */
    while ( coincHead )
    {
      thisCoinc = coincHead;
      coincHead = coincHead->next;
      XLALFreeCoincInspiral( &thisCoinc );
    }

    /* move events to coincHead */
    coincHead = slideCoinc;
    slideCoinc = NULL;
  }
    

  
  /*
   *
   * cluster the remaining events
   *
   */


  if ( coincHead && cluster_dt )
  {
    if ( vrbflg ) fprintf( stdout, "clustering remaining triggers... " );

    if ( !numSlides )
    {
      numClusteredEvents = XLALClusterCoincInspiralTable( &coincHead, 
          cluster_dt, coincstat , &bittenLParams);
    }
    else
    { 
      int slide = 0;
      int numClusteredSlide = 0;
      CoincInspiralTable *slideCoinc = NULL;
      CoincInspiralTable *slideClust = NULL;
      
      if ( vrbflg ) fprintf( stdout, "splitting events by slide\n" );

      for( slide = -numSlides; slide < (numSlides + 1); slide++)
      {
        if ( vrbflg ) fprintf( stdout, "slide number %d; ", slide );
        /* extract the slide */
        slideCoinc = XLALCoincInspiralSlideCut( &coincHead, slide );
        /* run clustering */
        numClusteredSlide = XLALClusterCoincInspiralTable( &slideCoinc, 
          cluster_dt, coincstat, &bittenLParams);
        
        if ( vrbflg ) fprintf( stdout, "%d clustered events \n", 
          numClusteredSlide );
        numClusteredEvents += numClusteredSlide;

        /* add clustered triggers */
        if( slideCoinc )
        {
          if( slideClust )
          {
            thisCoinc = thisCoinc->next = slideCoinc;
          }
          else
          {
            slideClust = thisCoinc = slideCoinc;
          }
          /* scroll to end of list */
          for( ; thisCoinc->next; thisCoinc = thisCoinc->next);
        }
      }

      /* free coincHead -- although we expect it to be empty */
      while ( coincHead )
      {
        thisCoinc = coincHead;
        coincHead = coincHead->next;
        XLALFreeCoincInspiral( &thisCoinc );
      }

      /* move events to coincHead */
      coincHead = slideClust;
      slideClust = NULL;
    }
    
    if ( vrbflg ) fprintf( stdout, "done\n" );
    if ( vrbflg ) fprintf( stdout, "%d clustered events \n", 
        numClusteredEvents );
  }


  /*
   *
   * update search_summary->nevents with an authoritative count of coincs
   *
   */

  searchSummList->nevents = 0;
  thisCoinc = coincHead;
  while (thisCoinc) {
    searchSummList->nevents += 1;
    thisCoinc = thisCoinc->next;
  }

  /*
   *
   * write output data
   *
   */


  /* write out all coincs as singles with event IDs */
  snglOutput = XLALExtractSnglInspiralFromCoinc( coincHead, 
      NULL, 0);


  /* write the main output file containing found injections */
  if ( vrbflg ) fprintf( stdout, "writing output xml files... " );
  memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlStream, outputFileName ), 
      &status );

  /* write out the process and process params tables */
  if ( vrbflg ) fprintf( stdout, "process... " );
  XLALGPSTimeNow(&(proctable.processTable->end_time));
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, process_table ), 
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, proctable, 
        process_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );

  /* erase the first empty process params entry */
  {
    ProcessParamsTable *emptyPPtable = procparams.processParamsTable;
    procparams.processParamsTable = procparams.processParamsTable->next;
    free( emptyPPtable );
  }

  /* write the process params table */
  if ( vrbflg ) fprintf( stdout, "process_params... " );
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, 
        process_params_table ), &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, procparams, 
        process_params_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );

  /* write search_summary table */
  if ( vrbflg ) fprintf( stdout, "search_summary... " );
  outputTable.searchSummaryTable = searchSummList;
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, 
        search_summary_table ), &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, outputTable, 
        search_summary_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );

  /* write the search_summvars table */
  if ( vrbflg ) fprintf( stdout, "search_summvars... " );
  LAL_CALL( LALBeginLIGOLwXMLTable( &status ,&xmlStream,
        search_summvars_table), &status );
  searchSummvarsTable.searchSummvarsTable = inputFiles;
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, searchSummvarsTable,
        search_summvars_table), &status );
  LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream), &status );

  /* write summ_value table */
  if ( summValueList )
  {
    if ( vrbflg ) fprintf( stdout, "search_summary... " );
    outputTable.summValueTable = summValueList;
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, 
          summ_value_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, outputTable, 
          summ_value_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );
  }

  /* Write the found injections to the sim table */
  if ( simEventHead )
  {
    if ( vrbflg ) fprintf( stdout, "sim_inspiral... " );
    outputTable.simInspiralTable = simEventHead;
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, 
          sim_inspiral_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, outputTable, 
          sim_inspiral_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream ), &status );
  }

  /* Write the results to the inspiral table */
  if ( snglOutput )
  {
    if ( vrbflg ) fprintf( stdout, "sngl_inspiral... " );
    outputTable.snglInspiralTable = snglOutput;
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, 
          sngl_inspiral_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, outputTable, 
          sngl_inspiral_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream ), &status);
  }

  /* close the output file */
  LAL_CALL( LALCloseLIGOLwXMLFile(&status, &xmlStream), &status);
  if ( vrbflg ) fprintf( stdout, "done\n" );

  if ( missedFileName )
  {
    /* open the missed injections file and write the missed injections to it */
    if ( vrbflg ) fprintf( stdout, "writing missed injections... " );
    memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
    LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlStream, missedFileName ), 
        &status );

    /* write out the process and process params tables */
    if ( vrbflg ) fprintf( stdout, "process... " );
    XLALGPSTimeNow(&(proctable.processTable->end_time));
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, process_table ),
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, proctable,
          process_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );

    /* write the process params table */
    if ( vrbflg ) fprintf( stdout, "process_params... " );
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream,
         process_params_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, procparams,
          process_params_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );

    /* write search_summary table */
    if ( vrbflg ) fprintf( stdout, "search_summary... " );
    outputTable.searchSummaryTable = searchSummList;
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream,
          search_summary_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, outputTable,
          search_summary_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );

    if ( missedSimHead )
    {
      outputTable.simInspiralTable = missedSimHead;
      LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, 
            sim_inspiral_table ), &status );
      LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, outputTable, 
            sim_inspiral_table ), &status );
      LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream ), &status );
    }

    LAL_CALL( LALCloseLIGOLwXMLFile( &status, &xmlStream ), &status );
    if ( vrbflg ) fprintf( stdout, "done\n" );
  }
  
  /* free process proctable data after found and missed files have been written */  
  free( proctable.processTable );

  if ( summFileName )
  {
    LIGOTimeGPS triggerTime;

    /* write out a summary file */
    fp = fopen( summFileName, "w" );

    switch ( dataType )
    {
      case playground_only:
        fprintf( fp, "using data from playground times only\n" );
        break;
      case exclude_play:
        fprintf( fp, "excluding all triggers in playground times\n" );
        break;
      case all_data:
        fprintf( fp, "using all input data\n" );
        break;
      default:
        fprintf( stderr, "data set not defined\n" );
        exit( 1 );
    }

    fprintf( fp, "read triggers from %d files\n", numInFiles );
    fprintf( fp, "number of triggers in input files: %d \n", numTriggers );

    if ( rsqVetoThresh > 0 )
    {
      fprintf( fp, "performed R-squared veto on triggers with snr < %f\n",
          rsqMaxSnr);
      fprintf( fp, "with rsqveto_duration below %f\n",
          rsqVetoThresh);
      if ( (rsqAboveSnrCoeff > 0) && (rsqAboveSnrPow > 0) )
      {
        fprintf( fp, "and on triggers with snr > %f\n",
            rsqMaxSnr);
        fprintf( fp, "with rsqveto_duration above %f * snr ^ %f\n",
            rsqAboveSnrCoeff, rsqAboveSnrPow );
      }
      fprintf( fp, "the number of triggers below the R-squared veto are: %d \n",
          numEventsBelowRsqThresh);
    }

    fprintf( fp, "number of reconstructed coincidences: %d \n", numCoincs );
    if ( ifos )
    {
      fprintf( fp, "number of triggers from %s ifos: %d \n", ifos, 
          numEventsInIfos );
    }
    if ( dataType != all_data )
    {
      fprintf( fp, "number of triggers after playground test: %d \n", 
          numEventsPlayTest);
    }
    if ( statThreshold )
    {
      fprintf( fp, "number of triggers with statistic above %6.2f is: %d \n", 
          statThreshold, numEventsAboveThresh);
    }
    XLALINT8NSToGPS( &triggerTime, triggerInputTimeNS );
    fprintf( fp, "amount of time analysed for triggers %d sec %d ns\n", 
        triggerTime.gpsSeconds, triggerTime.gpsNanoSeconds );

    if ( injectFileName )
    {
      fprintf( fp, "read %d injections from file %s\n", 
          numSimEvents, injectFileName );

      fprintf( fp, "number of injections in input data: %d\n", numSimInData );
      fprintf( fp, "number of injections found in input data: %d\n", 
          numCoincFound );
      fprintf( fp, 
          "number of triggers found within %lld msec of injection: %d\n",
          (injectWindowNS / 1000000LL), numEventsCoinc );

      fprintf( fp, "efficiency: %f \n", 
          (REAL4) numCoincFound / (REAL4) numSimInData );
    }

    if ( extractSlide )
    {
      fprintf( fp, "kept only triggers from slide %d\n", extractSlide );
    }

    if ( cluster_dt )
    {
      if ( numSlides )
      {
        fprintf( fp, "clustering triggers from %d slides separately\n",
            numSlides );
      }
      fprintf( fp, "number of event clusters with %lld msec window: %d\n",
          cluster_dt/ 1000000LL, numClusteredEvents ); 
    }

    fclose( fp ); 
  }


  /*
   *
   * free memory and exit
   *
   */


  /* free the coinc inspirals */
  while ( coincHead )
  {
    thisCoinc = coincHead;
    coincHead = coincHead->next;
    XLALFreeCoincInspiral( &thisCoinc );
  }

  /* free the inspiral events we saved */
  while ( inspiralEventList )
  {
    thisSngl = inspiralEventList;
    inspiralEventList = inspiralEventList->next;
    XLALFreeSnglInspiral( &thisSngl );
  }

  while ( snglOutput )
  {
    thisInspiralTrigger = snglOutput;
    snglOutput = snglOutput->next;
    XLALFreeSnglInspiral( &thisInspiralTrigger );
  }

  /* free the process params */
  while( procparams.processParamsTable )
  {
    this_proc_param = procparams.processParamsTable;
    procparams.processParamsTable = this_proc_param->next;
    free( this_proc_param );
  }

  while ( summValueList )
  {
    SummValueTable *thisSummValue;
    thisSummValue = summValueList;
    summValueList = summValueList->next;
    LALFree( thisSummValue );
  } 

  /* free the found injections */
  while ( simEventHead )
  {
    thisSimEvent = simEventHead;
    simEventHead = simEventHead->next;
    LALFree( thisSimEvent );
  }

  /* free the temporary memory containing the missed injections */
  while ( missedSimHead )
  {
    tmpSimEvent = missedSimHead;
    missedSimHead = missedSimHead->next;
    LALFree( tmpSimEvent );
  }

  /* free input files list */
  while ( inputFiles )
  {
    thisInputFile = inputFiles;
    inputFiles = thisInputFile->next;
    LALFree( thisInputFile );
  }

  /* free search summaries read in */
  while ( searchSummList )
  {
    thisSearchSumm = searchSummList;
    searchSummList = searchSummList->next;
    LALFree( thisSearchSumm );
  }

  if ( vrbflg ) fprintf( stdout, "checking memory leaks and exiting\n" );
  LALCheckMemoryLeaks();
  exit( 0 );
}

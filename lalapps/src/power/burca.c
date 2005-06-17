#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOMetadataUtils.h>
#include <lalapps.h>

int snprintf(char *str, size_t size, const char *format, ...);
long long int llabs(long long int j);

RCSID("$Id$");

#define INCA_EARG   1
#define INCA_EROW   2
#define INCA_EFILE  3

#define INCA_MSGEARG   "Error parsing arguments"
#define INCA_MSGROW    "Error reading row from XML table"
#define INCA_MSGEFILE  "Could not open file"
#define INCA_TIMEWARNING "only triggers before 00:00 Dec 31, \
2010 will be used"
#define BURCA_SEARCHSUMWARNING "The Search summary info of \
of the two ifos do not match "

#define TRUE     1
#define FALSE    0 
#define MAXIFO   16
#define MAXFILES 128
#define MSEC   (1000000LL)

/* cluster options */
enum { undefined, clusterbypeaktimeandfreq, clusterbytimeandfreq } clusterchoice = undefined;

/* Usage format string. */
#define USAGE "Usage: %s --ifo-a trigfile.a --ifo-b trigfile.b \
--start-time startCoincidence --stop-time stopCoincidence \
--drho-plus dRhoPlus --drho-minus dRhoMinus  --slide-time sec.s to slide \
--slide-time-ns nsec.s to slide --number-slides no. of slides \
--dt deltat [--noncoincident] [--noplayground] [--help] [--user-tag]\n"

/************************************************************************************
 * The main program
 ***********************************************************************************/
int main(int argc, char **argv)
{
  static LALStatus       stat;
  INT4                   usePlayground = TRUE;
  INT4                   verbose = FALSE;
  INT4                   cluster = FALSE;
  INT4                   noncoincident = FALSE;
  INT4                   ignorePlayground = FALSE;
  INT4                   ignoreTFcomparison = FALSE;
  INT4                   ignoreTcomparison = FALSE;

  CHAR                   *comment = "";

  INT4                   isPlay = 0;
  INT4                   deltaT=0;
  INT8                   triggerTime = 0;
  INT8                   startCoincidenceNS=0;
  INT8                   endCoincidenceNS=0;
  LIGOTimeGPS            startCoincidence={0,0};
  LIGOTimeGPS            endCoincidence={977788813,0};

  LIGOTimeGPS            slideData={0,0};
  INT8                   slideDataNS=0;
  INT8                   slideStep = 0;
  INT8                   numSlides = 0;
  INT8                   slideNum = 0;
  INT8                   dtSlide = 0;

  CHAR                 **trigFile;
  INT4                   nTrigFile[MAXIFO];

  CHAR                   outfileName[512];
  CHAR                   outnoncfileName[512];

  SnglBurstTable         *currentEvent=NULL,*tmpEvent=NULL;
  SnglBurstTable         *prevEvent=NULL,*prevNonCEvent=NULL,*outEvent=NULL;
  SnglBurstTable         *coincidentEvents=NULL;
  SnglBurstTable         *noncoincidentEvents = NULL;
  SnglBurstTable        **burstEventList=NULL;
  SnglBurstTable        **currentTrigger=NULL;
 
  SnglBurstAccuracy       accParams;
  SearchSummaryTable      **searchSummary=NULL;
  MetadataTable           myTable;
  LIGOLwXMLStream         xmlStream;
  INT4                    i, j;

  CHAR                   *outputdir = "./";

  /* getopt arguments */
  struct option long_options[] =
    {
      /* these options set a flag */
      {"verbose",                 no_argument,       &verbose,           1 },
      {"noplayground",            no_argument,       &usePlayground,     0 },
      {"noncoincident",           no_argument,       &noncoincident,     1 },
      {"ignore-playground",       no_argument,       &ignorePlayground,  1 },
      {"ignore-tfcomparison",     no_argument,       &ignoreTFcomparison,  1 },
      {"ignore-tcomparison",      no_argument,       &ignoreTcomparison, 1 },
      /* parameters used to generate calibrated power spectrum */
      {"ifo-a",                   required_argument, 0,                'a'},
      {"ifo-b",                   required_argument, 0,                'b'},
      {"drhoplus",                required_argument, 0,                'c'},
      {"clustertype",             required_argument, 0,                'f'},
      {"drhominus",               required_argument, 0,                'd'},
      {"dt",                      required_argument, 0,                't'},
      {"start-time",              required_argument, 0,                'r'},
      {"stop-time",               required_argument, 0,                's'},
      {"slide-time",              required_argument, 0,                'X'},
      {"slide-time-ns",           required_argument, 0,                'Y'},
      {"output-dir",              required_argument, 0,                'o'},
      {"number-slides",           required_argument, 0,                'k'},
      {"user-tag",                required_argument, 0,                'i'},
      {"help",                    no_argument,       0,                'h'}, 
      {0, 0, 0, 0}
    };
  int c;

  /*******************************************************************
   * initialize things
   *******************************************************************/
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "2" );
  memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );

  trigFile = (CHAR **) LALCalloc( MAXIFO * MAXFILES, sizeof(CHAR*) );
  burstEventList = (SnglBurstTable **) LALCalloc( 2, 
						  sizeof(SnglBurstTable) );
  currentTrigger    = (SnglBurstTable **) LALCalloc( 2, 
						     sizeof(SnglBurstTable) );
 
  searchSummary = (SearchSummaryTable **) LALCalloc( 2,
						    sizeof(SearchSummaryTable) );

  memset( &accParams, 0, sizeof(SnglBurstAccuracy) );
  memset( &nTrigFile, 0, MAXIFO * sizeof(INT4) );

  /*******************************************************************
   * BEGIN PARSE ARGUMENTS (inarg stores the current position)        *
   *******************************************************************/
  /* parse the arguments */
  while ( 1 )
    {
      /* getopt_long stores long option here */
      int option_index = 0;

      c = getopt_long_only( argc, argv, 
			    "a:b:c:d:f:t:r:s:i:j:k:X:Y:o:h", long_options, &option_index );

      /* detect the end of the options */
      if ( c == - 1 )
        {
	  break;
        }

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

                
	case 'a':
	  /* trigger files for interferometer A */
	  {
	    size_t len = strlen(optarg) + 1;
	    trigFile[nTrigFile[0]] = 
	      (CHAR *) LALCalloc( len, sizeof(CHAR) );
	    memcpy( trigFile[nTrigFile[0]], optarg, len);
	    nTrigFile[0]++;
	  }
	  break;

	case 'b':
	  /* trigger files for interferometer B */
	  {
	    size_t len = strlen(optarg) + 1;
	    *(trigFile+MAXFILES+nTrigFile[1]) = 
	      (CHAR *) LALCalloc( len, sizeof(CHAR) );
	    memcpy( *(trigFile+MAXFILES+nTrigFile[1]), optarg, len);
	    nTrigFile[1]++;
	  }
	  break;

	case 'c':
	  /* SNR error upward */
	  accParams.dRhoPlus = atof(optarg);
	  break;

	case 'd':
	  /* SNR error downward */
	  accParams.dRhoMinus = atof(optarg);
	  break;

	case 't':
	  /* time coincidence window */
	  accParams.dtime = atof(optarg) * MSEC;
	  deltaT = accParams.dtime;
	  break;
           
	case 'r':
	  /* time coincidence window */
	  startCoincidence.gpsSeconds = atoi(optarg);
	  LAL_CALL( LALGPStoINT8( &stat, &startCoincidenceNS, &startCoincidence ), &stat );
	  break;

	case 's':
	  /* time coincidence window */
	  endCoincidence.gpsSeconds = atoi(optarg);
	  LAL_CALL( LALGPStoINT8( &stat, &endCoincidenceNS, &endCoincidence ), &stat );
	  break;

	case 'o':
	  /* output directory for burca files */
	  outputdir = optarg;
	  break;

	case 'X':
	  slideData.gpsSeconds = (INT4) atoi( optarg );
	  break;

	case 'Y':
	  slideData.gpsNanoSeconds = (INT4) atoi( optarg );
	  break;

	case 'f':
	  /*
	   * set the cluster option
	   */			
	  {
	    if ( ! strcmp( "clusterbypeaktimeandfreq", optarg ) )
	      {
		cluster = TRUE;
		clusterchoice = clusterbypeaktimeandfreq;
	      }
	    else if ( ! strcmp( "clusterbytimeandfreq", optarg ) )
	      {
		cluster = TRUE;
		clusterchoice = clusterbytimeandfreq;
	      }
	    else
	      {
		fprintf( stderr, "invalid argument to --clustertype\n"
			 "unknown clusterchoice specified;\n"
			 "(must be one of:clusterbypeaktimeandfreq ,clusterbytimeandfreq )\n");
	      }
	  }
	  break;		    
		    
	case 'h':
	  /* help message */
	  fprintf( stderr, USAGE , argv[0]);
	  exit( 1 );
	  break;

	case 'i':
	  comment = optarg;
	  break;

	case 'k':
	  numSlides = atoi( optarg );
	  if( numSlides < 0 ){
	    fprintf( stderr, "--number-slides must be > 0\n");
	    exit ( 1 );
	  }
	  break;

	case '?':
	  fprintf( stderr, USAGE , argv[0]);
	  exit( 1 );
	  break;

	default:
	  fprintf( stderr, "unknown error while parsing options\n" );
	  fprintf( stderr, USAGE, argv[0] );
	  exit( 1 );
        }
    }

  if ( *trigFile == NULL || (*trigFile+MAXFILES) == NULL ){
    LALPrintError( "Must supply two input files\n" );
    return INCA_EARG;
  }

  if ( deltaT == 0 ){
    fprintf(stderr,"Warning: coincidence window has been set to 0 ms \n");
  }

  /*******************************************************************
   * END PARSE ARGUMENTS                                              *
   *******************************************************************/
  if (endCoincidence.gpsSeconds == 977788813)
    fprintf(stderr,"Warning: %s\n", INCA_TIMEWARNING);

  /****************************************************************
   * READ IN THE SEARCH SUMMARY INFO                              *
   ***************************************************************/
  for(j=0 ; j<2 ; j++)
    {
      SearchSummaryTableFromLIGOLw( &searchSummary[j], *(trigFile+j*MAXFILES));
    }

  /*****************************************************************
   * loop over input files for both ifos
   *****************************************************************/
  if (verbose) fprintf(stdout,"Looping over files\n");
  for(j=0 ; j<2 ; j++)
    {

      currentEvent = tmpEvent = burstEventList[j] = NULL;
      for(i=0; i<nTrigFile[j] ; i++)
        {
	  LAL_CALL( LALSnglBurstTableFromLIGOLw (&stat, &tmpEvent, 
						 *(trigFile+j*MAXFILES+i)), &stat);

	  /* connect results to linked list */
	  if (currentEvent == NULL)
            {
	      burstEventList[j] = currentEvent = tmpEvent;
            }
	  else
            {
	      currentEvent->next = tmpEvent;
            }

	  /* move to the end of the linked list for next input file */
	  if (currentEvent){
	    while (currentEvent->next != NULL)
	      {
		currentEvent = currentEvent->next;
	      }
	    tmpEvent = currentEvent->next;
	  }
	}

      /* sort the triggers in time */
      if ( burstEventList[j] )
        {
	  LAL_CALL( LALSortSnglBurst(&stat, &(burstEventList[j]), 
				     XLALCompareSnglBurstByStartTime ), &stat);
        }
    }

  /*****************************************************************
   * slide the data:  always slide burstEventList[1]
   * and put it on a ring with start and end time given by 
   * startCoincidence and endCoincidence ........
   *****************************************************************/
  LAL_CALL( LALGPStoINT8( &stat, &slideStep, &slideData ), &stat );

  for( slideNum = -numSlides; slideNum <= numSlides; slideNum++)
    {      
      if (slideNum == -numSlides)
	{
	  slideDataNS = (slideNum * slideStep);
	  dtSlide = slideDataNS;
	}
      else
	{
	  slideDataNS = slideStep;
	  dtSlide += slideStep;
	  /* if doing slides and reached 0-lag then skip it */
	  if ( numSlides && slideNum == 0 )
	    {
	      slideNum++;
	      slideDataNS += slideStep;
	      dtSlide += slideStep;
	    }	 
	}

      /* slide the data only if not 0 lag analysis*/    
      if (numSlides)
	XLALTimeSlideSnglBurst(burstEventList[1],startCoincidenceNS,endCoincidenceNS, slideDataNS);
	
      /*****************************************************************
       * find the first trigger after coincidence start time
       *****************************************************************/
      if (verbose) fprintf(stdout,"Moving to first trigger in window\n");
      for(j=0 ; j<2 ; j++)
	{
	  currentTrigger[j] = burstEventList[j];
	
	  while ( (currentTrigger[j] != NULL) && 
		  (currentTrigger[j]->start_time.gpsSeconds < startCoincidence.gpsSeconds) )
	    {
	      currentTrigger[j] = currentTrigger[j]->next;
	    }
	}

      /*****************************************************************
       * outer loop over triggers from interferometer A
       ****************************************************************/
      if (verbose) fprintf(stdout,"Start loop over ifo A\n");
      while ( (currentTrigger[0] != NULL) && 
	      (currentTrigger[0]->start_time.gpsSeconds < endCoincidence.gpsSeconds) )
	{
	  INT8 ta, tb;
	  INT8 coin = 0;
	
	  LAL_CALL( LALGPStoINT8(&stat, &ta, &(currentTrigger[0]->peak_time)), &stat);
	
	  /*catch up triggers from ifo B */
	  while (currentTrigger[1] != NULL)
	    {
	      LAL_CALL( LALGPStoINT8(&stat, &tb, &(currentTrigger[1]->peak_time)), &stat);
	      if (tb > ta-deltaT)
		{
		  break;
		}
	      currentTrigger[1] = currentTrigger[1]->next;
	    }
	
	  if (verbose) fprintf(stdout,"Start loop over ifo B\n");
	
	  /* look for coincident events in B within the time window */
	  LAL_CALL( LALGPStoINT8( &stat, &triggerTime, 
				  &(currentTrigger[0]->start_time) ), &stat );
	
	  LAL_CALL( LALINT8NanoSecIsPlayground( &stat, &isPlay, 
						&triggerTime ), &stat );
	
	  if ( verbose )
	    {
	      if ( isPlay )
		{
		  fprintf( stdout, "trigger is playground\n" );
		} 
	      else
		{
		  fprintf( stdout, "trigger is not playground\n" );
		}
	    }

	  /* if we are playground only and the trigger is in playground or
	   * we are not using playground and the trigger is not in the 
	   * playground  ... */
	
	  if ( ( usePlayground && isPlay ) || ( ! usePlayground && ! isPlay) || ignorePlayground )
	    {
	      tmpEvent = currentTrigger[1];
	      while (tmpEvent != NULL)
		{
		  LAL_CALL( LALGPStoINT8(&stat, &tb, &(tmpEvent->peak_time)), &stat);
		  if (tb > ta+deltaT)
		    break;

		  /* this is a LAL function which compares events */
		  if( ignoreTcomparison )
		    accParams.difference = XLALCompareSnglBurstByFreq((const SnglBurstTable * const *)&currentTrigger[0], (const SnglBurstTable * const *)&tmpEvent);
		  else if ( !ignoreTFcomparison )
		    accParams.difference = XLALCompareSnglBurst((const SnglBurstTable * const *)&currentTrigger[0], (const SnglBurstTable * const *)&tmpEvent);
		
		  if (!accParams.difference || ignoreTFcomparison)
		    {
		      coin = 1;
		      if (coincidentEvents == NULL)
			{
			  outEvent = coincidentEvents = (SnglBurstTable *)
			    LALCalloc(1, sizeof(SnglBurstTable) );
			  prevEvent = outEvent;
			}
		      else 
			{
			  outEvent = (SnglBurstTable *)
			    LALCalloc(1, sizeof(SnglBurstTable) );
			  prevEvent->next = outEvent;
			}
		      memcpy( outEvent, currentTrigger[0], sizeof(SnglBurstTable));
		      prevEvent = outEvent;
		      outEvent = outEvent->next = NULL;
		    }
		  tmpEvent = tmpEvent->next;
		}
	    }
	
	  if (noncoincident && coin == 0)
	    {
	      if (noncoincidentEvents == NULL)
		{
		  outEvent = noncoincidentEvents = (SnglBurstTable *)
		    LALCalloc(1, sizeof(SnglBurstTable) );
		  prevNonCEvent = outEvent;
		}
	      else 
		{
		  outEvent = (SnglBurstTable *)
		    LALCalloc(1, sizeof(SnglBurstTable) );
		  prevNonCEvent->next = outEvent;
		}
	      memcpy( outEvent, currentTrigger[0], sizeof(SnglBurstTable));
	      prevNonCEvent = outEvent;
	      outEvent = outEvent->next = NULL;
	    }
	
	  currentTrigger[0] = currentTrigger[0]->next;
	}



      /* cluster the events if asked to */
      if(cluster && clusterchoice == clusterbypeaktimeandfreq)
	{
	  if (noncoincident && noncoincidentEvents)
	    XLALClusterSnglBurstTable(&noncoincidentEvents, XLALCompareSnglBurstByPeakTime, XLALCompareSnglBurstByPeakTimeAndFreq, XLALSnglBurstCluster);
	  if(coincidentEvents)
	    XLALClusterSnglBurstTable(&coincidentEvents, XLALCompareSnglBurstByPeakTime, XLALCompareSnglBurstByPeakTimeAndFreq, XLALSnglBurstCluster);
	}
      else if (cluster && clusterchoice == clusterbytimeandfreq)
	{
	  if (noncoincident && noncoincidentEvents)
	    XLALClusterSnglBurstTable(&noncoincidentEvents, NULL, XLALCompareSnglBurst, XLALSnglBurstCluster);
	  if(coincidentEvents)
	    XLALClusterSnglBurstTable(&coincidentEvents,  NULL, XLALCompareSnglBurst, XLALSnglBurstCluster);
	}

      /*****************************************************************
       * open output xml file to write coincident triggers
       *****************************************************************/
      if ( dtSlide < 0 )
	  snprintf(outfileName, sizeof(outfileName)-1, "%s/%s-BURCA_%s%s_M_%lld_%s-%d-%d.xml", outputdir,searchSummary[0]->ifos,searchSummary[0]->ifos,searchSummary[1]->ifos,-dtSlide,comment,searchSummary[0]->in_start_time.gpsSeconds,searchSummary[0]->in_end_time.gpsSeconds - searchSummary[0]->in_start_time.gpsSeconds);	  
      else
	  snprintf(outfileName, sizeof(outfileName)-1, "%s/%s-BURCA_%s%s_P_%lld_%s-%d-%d.xml", outputdir,searchSummary[0]->ifos,searchSummary[0]->ifos,searchSummary[1]->ifos,dtSlide,comment,searchSummary[0]->in_start_time.gpsSeconds,searchSummary[0]->in_end_time.gpsSeconds - searchSummary[0]->in_start_time.gpsSeconds);

      outfileName[sizeof(outfileName)-1] = '\0';
    
      LAL_CALL( LALOpenLIGOLwXMLFile(&stat, &xmlStream, outfileName), &stat);
    
      /* write the search cummary table[Note currently we only keep events from Ifo-a */
      snprintf(searchSummary[0]->comment, LIGOMETA_COMMENT_MAX, "%s", comment); 
      searchSummary[0]->nevents = XLALCountSnglBurst(coincidentEvents);   
      LAL_CALL( LALBeginLIGOLwXMLTable( &stat, &xmlStream, search_summary_table ), 
		&stat );
      myTable.searchSummaryTable = searchSummary[0];
      LAL_CALL( LALWriteLIGOLwXMLTable( &stat, &xmlStream, myTable, 
					search_summary_table ), &stat );
      LAL_CALL( LALEndLIGOLwXMLTable ( &stat, &xmlStream ), &stat );
    
      /*write the triggers */
      LAL_CALL( LALBeginLIGOLwXMLTable (&stat, &xmlStream, sngl_burst_table), &stat);
      myTable.snglBurstTable = coincidentEvents;
      LAL_CALL( LALWriteLIGOLwXMLTable (&stat, &xmlStream, myTable,
					sngl_burst_table), &stat);
      LAL_CALL( LALEndLIGOLwXMLTable (&stat, &xmlStream), &stat);
      LAL_CALL( LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);

      /*****************************************************************
       * open output xml file to write non coincident triggers
       *****************************************************************/
      if (noncoincident)
	{ 
	  if ( dtSlide < 0 )
	    snprintf(outnoncfileName, sizeof(outnoncfileName)-1, "%s/%s-NEGBURCA_%s%s_M_%lld_%s-%d-%d.xml", outputdir,searchSummary[0]->ifos,searchSummary[0]->ifos,searchSummary[1]->ifos,-dtSlide,comment,searchSummary[0]->in_start_time.gpsSeconds,searchSummary[0]->in_end_time.gpsSeconds - searchSummary[0]->in_start_time.gpsSeconds);	  
	  else
	    snprintf(outnoncfileName, sizeof(outnoncfileName)-1, "%s/%s-NEGBURCA_%s%s_P_%lld_%s-%d-%d.xml", outputdir,searchSummary[0]->ifos,searchSummary[0]->ifos,searchSummary[1]->ifos,dtSlide,comment,searchSummary[0]->in_start_time.gpsSeconds,searchSummary[0]->in_end_time.gpsSeconds - searchSummary[0]->in_start_time.gpsSeconds );
	  
	  outnoncfileName[sizeof(outnoncfileName)-1] = '\0';
	
	  LAL_CALL( LALOpenLIGOLwXMLFile(&stat, &xmlStream, outnoncfileName), &stat);

	  /* write the search cummary table */
	  snprintf(searchSummary[0]->comment, LIGOMETA_COMMENT_MAX, "%s", comment); 
	  searchSummary[0]->nevents = XLALCountSnglBurst(noncoincidentEvents);    
	  LAL_CALL( LALBeginLIGOLwXMLTable( &stat, &xmlStream, search_summary_table ), 
		    &stat );
	  myTable.searchSummaryTable = searchSummary[0];
	  LAL_CALL( LALWriteLIGOLwXMLTable( &stat, &xmlStream, myTable, 
					    search_summary_table ), &stat );
	  LAL_CALL( LALEndLIGOLwXMLTable ( &stat, &xmlStream ), &stat );      

	
	  /*write the triggers */
      
	  LAL_CALL( LALBeginLIGOLwXMLTable (&stat, &xmlStream, sngl_burst_table), &stat);
	  myTable.snglBurstTable = noncoincidentEvents;
	  LAL_CALL( LALWriteLIGOLwXMLTable (&stat, &xmlStream, myTable,
					    sngl_burst_table), &stat);
	  LAL_CALL( LALEndLIGOLwXMLTable (&stat, &xmlStream), &stat);
	  LAL_CALL( LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);
	}
 
      /*****************************************************************
       * clean up the memory that has been allocated 
       *****************************************************************/
      while ( coincidentEvents != NULL)
	{
	  prevEvent = coincidentEvents;
	  coincidentEvents = coincidentEvents->next;
	  LALFree( prevEvent );
	}
    
      if (noncoincident)
	{
	  while ( noncoincidentEvents != NULL)
	    {
	      prevEvent = noncoincidentEvents;
	      noncoincidentEvents = noncoincidentEvents->next;
	      LALFree( prevEvent );
	    } 
	}
    }

  for(j=0; j<2 ; j++)
    {
      tmpEvent = burstEventList[j];
      while (tmpEvent)
        {
	  prevEvent = tmpEvent;
	  tmpEvent = tmpEvent->next;
	  LALFree( prevEvent );
        }
      for (i=0; i<nTrigFile[j] ; i++){
	LALFree( trigFile[j*MAXFILES +i] );
      }

      while (searchSummary[j])
	{
	  SearchSummaryTable *thisEvent;
	  thisEvent = searchSummary[j];
	  searchSummary[j] = searchSummary[j]->next;
	  LALFree( thisEvent );
	}
    }

  LALFree( currentTrigger );
  LALFree( burstEventList );
  LALFree( trigFile );
  LALFree( searchSummary );
 
  LALCheckMemoryLeaks();

  return 0;
}

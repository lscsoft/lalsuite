/*----------------------------------------------------------------------- 
 * 
 * File Name: inspinj_find.c
 *
 * Author: Fairhurst, S
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/Date.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOLwXMLRead.h>
#include <lalapps.h>
#include <processtable.h>
#include "event_utils.h"


RCSID("$Id$");

#define MAXSTR 2048

/* Usage format string. */
#define USAGE \
"%s [options]\n"\
"  --help                       display this message\n"\
"  --input file                 read list of input XML files from file\n"\
"  --inject file                read list of simulated injections from file\n"\
"  --outfile file               write found injections to file\n"\
"  --missedinjections file      write missed injections to file\n"\
"  --snrstar thresh             discard triggers with snr less than thresh\n"\
"  --noplayground               use all data, not just playground\n"\
"  --hardware gpstime           assume hardware injections starting at gpstime\n"\
"  --deltat t                   trigger and injection time coinc window (ms)\n"\
"  --cluster t                  cluster coincident triggers with t ms window\n"\
"  --clusteralgorithm alg       use cluster algorithm alg\n"\
"\n"

#define INSPINJFIND_EARG   1
#define INSPINJFIND_EROW   2
#define INSPINJFIND_EFILE  3

#define INSPINJFIND_MSGEARG   "Error parsing arguments"
#define INSPINJFIND_MSGROW    "Error reading row from XML table"
#define INSPINJFIND_MSGEFILE  "Could not open file"

#define TRUE  1
#define FALSE 0

#define PROGRAM_NAME "inspinj_find"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
  this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
  calloc( 1, sizeof(ProcessParamsTable) ); \
  LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
  PROGRAM_NAME ); \
  LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
  long_options[option_index].name ); \
  LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
  LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );


static int getline(char *line, int max, FILE *fpin)
{
  int i;
  CHAR tmpline[MAXSTR];

  for (i=0;i<MAXSTR;i++) { *(line+i) = '\0' ; }
  if (fgets(tmpline, max, fpin) == NULL){
    return 0;
  }
  else{
    strncpy(line,tmpline,strlen(tmpline)-1);
    return strlen(line);
  }
}


/****************************************************************************
 * 
 * FUNCTION TESTS IF THE FILE CONTAINS ANY PLAYGROUND DATA
 * 
 ***************************************************************************/
static int isPlayground(INT4 gpsStart, INT4 gpsEnd){
  INT4 runStart=729273613;
  INT4 playInterval=6370;
  INT4 playLength=600;
  INT4 segStart,segEnd,segMiddle;

  segStart = (gpsStart - runStart)%playInterval;
  segEnd   = (gpsEnd - runStart)%playInterval;
  segMiddle = gpsStart + (INT4) (0.5 * (gpsEnd - gpsStart));
  segMiddle = (segMiddle - runStart)%playInterval;

  if (segStart < playLength || segEnd < playLength || segMiddle < playLength)
  {
    return TRUE;
  }

  return FALSE;
}

/*****************************************************************************
 *
 * The main program
 *
 *****************************************************************************/

int main(int argc, char **argv)
{
  static LALStatus      stat;
  LALLeapSecAccuracy	accuracy = LALLEAPSEC_LOOSE;
  INT4                  retVal=0;
  FILE                 *fpin=NULL;
  INT4                  fileCounter=0;
  BOOLEAN               playground=TRUE;
  INT4                  cluster=FALSE;
  Clusterchoice	        clusterchoice;
  INT4		        saveMissedInjections=FALSE;
  INT4		        coincidence = FALSE;
  INT4		        hardware = FALSE;

  INT8		        simTime;
  INT8		        inspiralTime;
  INT8		        starttime;

  REAL4                 snrstar=0.0;
  INT4                  clust=0;
  INT4		        dt = 20;

  CHAR                  inputfile[MAXSTR],line[MAXSTR];

  CHAR                 *outfileName=NULL;
  CHAR		       *injectfile=NULL;
  CHAR		       *missedInjectionsFile=NULL;

  const CHAR           *searchSummaryName="search_summary";
 
  SearchSummaryIndex    searchSummaryIndex;
  SearchSummaryTable    searchSummaryTable;

  SnglInspiralTable    *inspiralEventList=NULL;
  SnglInspiralTable    *currentInspiralEvent=NULL,*prevInspiralEvent=NULL;
  SnglInspiralTable   **inspiralHandle=NULL;
      
  SimInspiralTable     *currentSimEvent=NULL,*prevSimEvent=NULL;
  SimInspiralTable     *simEventList=NULL;
  SimInspiralTable     *simEventMissed=NULL;
  SimInspiralTable     *prevSimEventMissed=NULL;
  SimInspiralTable    **simHandle=NULL;

  MetadataTable         myTable;
  MetadataTable	        simFoundTable;
  MetadataTable	        simMissedTable;
  ProcessParamsTable   *this_proc_param;
  MetadataTable	        proctable;
  MetadataTable	        procparams;
  LIGOLwXMLStream       xmlStream;

  INT4		        numEvents = 0;
  INT4			currentNumEvents = 0;
  INT4			numKeptEvents = 0;
  INT4		        numInjects = 0; 
  INT4			numPlayInjects = 0;
  INT4		        numKeptInjects = 0; 
  INT4		        numFound = 0;

  static int		verbose_flag;

  struct MetaioParseEnvironment triggerParseEnv;
  const MetaioParseEnv triggerEnv = &triggerParseEnv;

  int			c;

  /********************************************************************
   * BEGIN PARSE ARGUMENTS						*
   ********************************************************************/

  /* set up inital debugging values */
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "33" );

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) 
    calloc( 1, sizeof(ProcessTable) );
  LAL_CALL( LALGPSTimeNow ( &stat, &(proctable.processTable->start_time),
        &accuracy ), &stat );
  LAL_CALL( populate_process_table( &stat, proctable.processTable, 
        PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE ), &stat );
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *) 
    calloc( 1, sizeof(ProcessParamsTable) );


  while (1)
  {
    /* getopt arguments */
    static struct option long_options[] = 
    {
      /* these options set a flag */
      {"verbose",               no_argument,	&verbose_flag, 1 },
      {"input",                 required_argument,  0,	'a'},
      {"inject",                required_argument,  0,	'b'},
      {"outfile",               required_argument,  0,	'c'},
      {"snrstar",               required_argument,  0,	'd'},
      {"noplayground",          no_argument,	    0,	'f'},
      {"deltat",                required_argument,  0,  'g'},
      {"help",                  no_argument,	    0,	'h'},    
      {"cluster",               required_argument,  0,	'j'},
      {"clusteralgorithm",      required_argument,  0,	'k'},
      {"missedinjections",      required_argument,  0,	'l'},
      {"hardware",              required_argument,  0,	'm'},
      {0, 0, 0, 0}
    };

    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long ( argc, argv, "a:b:c:d:e:f:g:h:j:k:l:", 
        long_options, &option_index );

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

      case 'a':
        /* file containing list of xml files to use */
        {
          sprintf(inputfile,"%s",optarg);
          ADD_PROCESS_PARAM( "string", "%s", inputfile );
        }
        break;	

      case 'b':
        /* give name of file containing the injections */
        {
          injectfile = optarg;
          ADD_PROCESS_PARAM( "string", "%s", injectfile );
        }
        break;

      case 'c':
        /* output file name */
        {
          outfileName = optarg;
          ADD_PROCESS_PARAM( "string", "%s", outfileName );
        }
        break;

      case 'd':
        /* a threshold cut on signal to noise, only triggers with
         * SNR > snrstar are kept */
        {
          snrstar = atof( optarg );
          ADD_PROCESS_PARAM( "float", "%e", snrstar );
        }
        break;

      case 'f':
        /* don't restrict to the playground data */
        {
          playground = FALSE;
          ADD_PROCESS_PARAM( "string", "%s", " " );
        }
        break;

      case 'g':
        /* provide the time discrepancy allowed when checking for 
         * coincidence of injections and inspirals 
         * If not specified, then default is 20ms */ 
        {
          dt = atoi( optarg );
          ADD_PROCESS_PARAM( "int", "%d", dt );
        }
        break;

      case 'h':
        /* print help */
        {
          LALPrintError( USAGE, *argv );
          return INSPINJFIND_EARG;
        }
        break;

      case 'j':
        /* cluster the events */ 
        {
          clust = atoi( optarg );
          cluster = TRUE;
          ADD_PROCESS_PARAM( "int", "%d", clust );
        }
        break;

      case 'k':
        /* choose the clustering algorithm, default is 0
         * options are detailed in event_utils.c */
        {	
          if ( ! strcmp( "snr_and_chisq", optarg ) )
          {
            clusterchoice = snr_and_chisq;
          }
          else if ( ! strcmp( "snrsq_over_chisq", optarg) )
          {
            clusterchoice = snrsq_over_chisq;
          }	    
          else
          {
            fprintf( stderr, "invalid argument to  --%s:\n"
                "unknown clustering specified:\n "
                "%s (must be one of: snr_and_chisq, \n"
                "   snrsq_over_chisq)\n",
                long_options[option_index].name, optarg);
            exit( 1 );
          }
          ADD_PROCESS_PARAM( "string", "%s", optarg );
        }
        break;

      case 'l':
        /* save a copy of those injected events which were missed, 
         * and specify the filename where those events should be stored */
        {
          missedInjectionsFile = optarg;
          saveMissedInjections = TRUE;
          ADD_PROCESS_PARAM( "string", "%s", missedInjectionsFile );
        }
        break;

      default:
        {
          return INSPINJFIND_EARG;
        }

      case 'm':
        /* if they are hardware injections, then need to pass in the start
         * time of the injections.  The inject file is assumed to contain a 
         * list of the hardware injections, where the time is the time after
         * starting injection that the coalescence occured */
        {
          starttime = atol( optarg );
          hardware = TRUE;
        }
        break;

    }   
  }

  if ( optind < argc )
  {
    fprintf( stderr, "extraneous command line arguments:\n" );
    while ( optind < argc )
    {
      fprintf ( stderr, "%s\n", argv[optind++] );
    }
    exit( 1 );
  }	  



  /* check validity of required arguments */

  /* check validity of input file */
  if ( inputfile == NULL)
  {
    LALPrintError( "Input file must be specified\n" );
    exit( 1 );
  }

  /* check outfile has been given */
  if ( ! outfileName )
  {
    LALPrintError( "Outfile name must be specified\n" );
    exit( 1 );
  }

  /*******************************************************************
   * END PARSE ARGUMENTS                                              *
   *******************************************************************/



  /*******************************************************************
   * initialize things
   *******************************************************************/
    memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );


  /*****************************************************************
   * OPEN THE FILE WITH THE INJECTIONS (if appropriate)
   ******************************************************************/ 

  if ( ! (injectfile == NULL))
  {	

    /****************************************************************
     *  open xml file at injection table 
     ****************************************************************/
    simHandle = &simEventList;
    numInjects = SimInspiralTableFromLIGOLw( simHandle, injectfile, 0, 0 );

    fprintf( stdout, "number of injections %d \n",numInjects);

    if ( numInjects < 0 )
    {
      fprintf( stderr, "error: unable to read sim_inspiral table from %s\n", 
	injectfile );
      exit( 1 );
    }

/*    if ( numInjects )
    {
      currentSimEvent = *simHandle;
      while ( currentSimEvent->next )
      {
        currentSimEvent = currentSimEvent->next;
      }
    }
*/    
    numPlayInjects = numInjects;
    currentSimEvent = simEventList;

    if (  playground || hardware )
    {
      numPlayInjects = 0;
      while( currentSimEvent->next )
      {
	/* if the injections are hardware, then add the starttime time to the 
	* time contained in sim_inspiral table */
	if ( hardware )
	{
	  currentSimEvent->geocent_end_time.gpsSeconds = starttime +
	    currentSimEvent->geocent_end_time.gpsSeconds;
	}

	if ( (playground) && 
	  (currentSimEvent->geocent_end_time.gpsSeconds-729273613)%6370 > 600 ) 
	{
	  /* injection not in playground so ditch */
	  currentSimEvent = currentSimEvent->next;
	}
	else
	{
	  /* keep injection */
	  ++numPlayInjects;
	  if (prevSimEvent == NULL)
	  {
	    simEventList = prevSimEvent = currentSimEvent;	
	    currentSimEvent = currentSimEvent->next;
	  }
	  else
	  {
	    prevSimEvent = prevSimEvent->next = currentSimEvent;
	    currentSimEvent = currentSimEvent->next;
	  }
	}
      }
      currentSimEvent = simEventList;   
      prevSimEvent = NULL;
    }
  }



  /*****************************************************************
   * OPEN FILE WITH LIST OF XML FILES (one file per line)
   ****************************************************************/
  if ( !(fpin = fopen(inputfile,"r")) )
  {
    LALPrintError("Could not open input file\n");
    exit (1);
  }


  /*****************************************************************
   * loop over the xml files
   *****************************************************************/
  inspiralHandle = &inspiralEventList;

  while ( getline(line, MAXSTR, fpin) )
  {
    fileCounter++;

    /**************************************************************
     *  open xml file at search summary table 
     **************************************************************/
    if ( (retVal = MetaioOpenTable( triggerEnv, line, 
            searchSummaryName)) !=0 )
    {
      fprintf(stderr, "Error opening injection file %s for %s\n", 
          line, searchSummaryName );
      if (fileCounter>1)
      {
        MetaioAbort( triggerEnv ); 
        exit(INSPINJFIND_EFILE);
      } 
      else 
      {
        fprintf(stderr,"Proceeding anyway\n");
      }
    }
    if ( !retVal )
    {
      /* get the search summary table */
      LAL_CALL( buildSearchSummaryIndex(&stat, triggerEnv, 
            &searchSummaryIndex), &stat);
      LAL_CALL( getSearchSummaryTable(&stat, triggerEnv, 
            &searchSummaryTable, &searchSummaryIndex), &stat);

      /* check for events and playground */
      if ( ( playground && 
            !(isPlayground(searchSummaryTable.out_start_time.gpsSeconds,
                searchSummaryTable.out_end_time.gpsSeconds ))) )
      {
        fprintf(stdout,"File %i %s: not in playground, continuing\n",
            fileCounter,line);
        continue;
      }
      else 
      {

        /* if there are injections, keep only those which occur
         * during the times searched in the input files */
        if( numPlayInjects != 0)
        {
          while (!( currentSimEvent == NULL) &&
              currentSimEvent->geocent_end_time.gpsSeconds < 
              searchSummaryTable.out_end_time.gpsSeconds)
          {
            /*  check if injection is before file start time */
            if (currentSimEvent->geocent_end_time.gpsSeconds < 
                searchSummaryTable.out_start_time.gpsSeconds)
            {
              /* discard the current injection */    
              if (prevSimEvent != NULL)
              {
                prevSimEvent->next = currentSimEvent->next;
                LALFree(currentSimEvent);
                currentSimEvent = prevSimEvent->next;
              }
              else
              {
                simEventList = simEventList->next;
                LALFree(currentSimEvent);
                currentSimEvent = simEventList;
              }
            }
            else 
            {
              /* keep the current injection */    
              numKeptInjects++;
              prevSimEvent = currentSimEvent;
              currentSimEvent = currentSimEvent->next;
            }
          }		
        }    
        /* check to see whether there are any events in the file */
        if ( searchSummaryTable.nevents == 0 )
        {
          fprintf(stdout,"File %i %s: no events, continuing\n",
              fileCounter,line);
          continue;
        }
        else 
        {
          fprintf(stdout,"File %i %s: processing\n", 
              fileCounter,line);
        }



      }
      /* close the stream */
      MetaioAbort( triggerEnv );
    }

    /**************************************************************
     *  open xml file at inspiral table 
     **************************************************************/
    currentNumEvents = 
      SnglInspiralTableFromLIGOLw( inspiralHandle, line, 0, -1 );

    if ( currentNumEvents < 0 )
    {
      fprintf( stderr, "error: unable to read sngl_inspiral table from %s\n", 
	line );
      exit( 1 );
    }

    if ( currentNumEvents )
    {
      numEvents += currentNumEvents;

      currentInspiralEvent = *inspiralHandle;
      while ( currentInspiralEvent->next )
      {
        currentInspiralEvent = currentInspiralEvent->next;
      }
      inspiralHandle = &(currentInspiralEvent->next);
    }
  }

  currentInspiralEvent = inspiralEventList;
  prevInspiralEvent = NULL;
  
  /* check for events and playground, and event satisfying SNR */
  while( currentInspiralEvent->next ) 
  {  
    if ( (( playground ) && 
          (currentInspiralEvent->end_time.gpsSeconds-729273613)%6370 > 600 )
	    || ( currentInspiralEvent->snr < snrstar) )
    {
      /* conditions not satisfied so ditch trigger */
      currentInspiralEvent = currentInspiralEvent->next;
    }
    else
    {
      /* keep trigger */
      ++numKeptEvents;
      if (prevInspiralEvent == NULL)
      {
	inspiralEventList = prevInspiralEvent = currentInspiralEvent;	
	currentInspiralEvent = currentInspiralEvent->next;
      }
      else
      {
	prevInspiralEvent = prevInspiralEvent->next = currentInspiralEvent;
	currentInspiralEvent = currentInspiralEvent->next;
      }
    }
  }
  
  /* sort the events */
  LAL_CALL( LALSortSnglInspiral( &stat, &inspiralEventList, 
    *LALCompareSnglInspiralByTime), &stat);

  

  /* discard the remaining injections which occured after the last file */
  while (currentSimEvent != NULL) 
  {
    if (prevSimEvent != NULL)
    {
      prevSimEvent->next = currentSimEvent->next;
      LALFree(currentSimEvent);
      currentSimEvent = prevSimEvent->next;
    }
    else
    {
      simEventList = simEventList->next;
      LALFree(currentSimEvent);
      currentSimEvent = simEventList;
    }
  }

  fprintf( stdout, "number of triggers %d \n", numEvents);
  fprintf( stdout, "number of triggers in playground with good snr %d \n",
      numKeptEvents);
  fprintf( stdout, "number of injections %d \n",numInjects);
  fprintf( stdout, "number of injections in playground %d \n", 
      numPlayInjects);	
  fprintf( stdout, "number of relevant injections %d \n", numKeptInjects);
  
  /*************************************************************
   * CHECK FOR COINCIDENCE OF ARRIVAL TIME BETWEEN INJECTIONS
   * AND FOUND INSPIRALS 
   **************************************************************/
  if ( ! (injectfile == NULL))
  {
    /* allowed delay time given by dt */

    numFound = 0;
    coincidence = FALSE;
    /*  Note: we are assuming that both the inspiral and injection events
     *  are time sorted */

    currentSimEvent = simEventList;
    prevSimEvent = NULL;
    currentInspiralEvent = inspiralEventList;
    prevInspiralEvent = NULL;

    /* loop over all injections */
    while ( currentInspiralEvent != NULL && currentSimEvent != NULL )
    {

      /* compute the time in nanosec for the injection */
      LALGPStoINT8(&stat, &simTime, 
          &(currentSimEvent->geocent_end_time));

      /* compute the time in nanosec for the inspiral */
      LALGPStoINT8(&stat, &inspiralTime, 
          &(currentInspiralEvent->end_time));
      /*  check if inspiral is before injection 
       *  (up to allowed difference -- dt) */
      if ( ((inspiralTime) < (simTime + 1000000LL * dt)) )
      {
        /*  check if inspiral is coincident with injection
         *  (up to allowed difference -- dt) */
        if ( ((inspiralTime) > (simTime - 1000000LL * dt)) )
        {
          /* record the fact that there's an inspiral coincident 
           * with the injection */
          if (coincidence == FALSE)
          {
            ++numFound;
            coincidence = TRUE;
          }
          /* keep the inspiral event */
          prevInspiralEvent = currentInspiralEvent;
          currentInspiralEvent = currentInspiralEvent->next;
        }
        else
          /* discard the inspiral event and check the next one,
           * do nothing to the injection event */
        {
          if (prevInspiralEvent != NULL)
          {
            prevInspiralEvent->next = currentInspiralEvent->next;
            LALFree(currentInspiralEvent);
            currentInspiralEvent = prevInspiralEvent->next;
          }

          else
          {
            inspiralEventList = inspiralEventList->next;
            LALFree(currentInspiralEvent);
            currentInspiralEvent = inspiralEventList;	
          }
        }
      }
      else if (coincidence == TRUE)
        /* keep the injection event */
      {
        prevSimEvent = currentSimEvent;
        currentSimEvent = currentSimEvent->next;
        coincidence = FALSE;
      }
      else if (coincidence == FALSE)
        /* Put the injection event into the simEventMissed list and 
         * check the next one, do nothing to the inspiral event. */
      {
        if (prevSimEvent != NULL)
        {
          prevSimEvent->next = currentSimEvent->next;
          currentSimEvent->next = NULL;
          if(prevSimEventMissed != NULL)
          {    
            prevSimEventMissed->next = currentSimEvent;
            prevSimEventMissed = prevSimEventMissed->next;
          }
          else
          {
            prevSimEventMissed = simEventMissed = currentSimEvent;
          }	
          currentSimEvent = prevSimEvent->next;		    
        }
        else
        {
          simEventList = simEventList->next;
          currentSimEvent->next = NULL;
          if(prevSimEventMissed != NULL)
          {    
            prevSimEventMissed->next = currentSimEvent;
            prevSimEventMissed = prevSimEventMissed->next;
          }
          else
          {
            prevSimEventMissed = simEventMissed = currentSimEvent;
          } 
          currentSimEvent = simEventList;
        }
      }
    }

    fprintf( stdout, "done sorting, number of coincidences found is %d\n", 
        numFound);

    if( currentInspiralEvent == NULL )
    {
      /* Remove all remaining injection events since they occur 
       * after the last recorded inspiral -- put them in the missed table 
       */

      while ( currentSimEvent != NULL )
      {
        if (prevSimEvent != NULL)
        {
          prevSimEvent->next = currentSimEvent->next;
          currentSimEvent->next = NULL;
          if(prevSimEventMissed != NULL)
          {    
            prevSimEventMissed->next = currentSimEvent;
            prevSimEventMissed = prevSimEventMissed->next;
          }
          else
          {
            prevSimEventMissed = simEventMissed = currentSimEvent;
          }	
          currentSimEvent = prevSimEvent->next;		    
        }
        else
        {
          simEventList = simEventList->next;
          currentSimEvent->next = NULL;
          if(prevSimEventMissed != NULL)
          {    
            prevSimEventMissed->next = currentSimEvent;
            prevSimEventMissed = prevSimEventMissed->next;
          }
          else
          {
            prevSimEventMissed = simEventMissed = currentSimEvent;
          } 
          currentSimEvent = simEventList;
        }

      }
    }

    if( currentSimEvent == NULL )
    {
      /* Remove all remaining inspiral events since they occur after 
       * the last recorded injection */
      while (currentInspiralEvent != NULL )	
      {
        if (prevInspiralEvent != NULL)
        {
          prevInspiralEvent->next = currentInspiralEvent->next;
          LALFree(currentInspiralEvent);
          currentInspiralEvent = prevInspiralEvent->next;
        }
        else
        {
          inspiralEventList = inspiralEventList->next;
          LALFree(currentInspiralEvent);
          currentInspiralEvent = inspiralEventList;	
        }
      }
    }
  }
  /*********************************************************************/


  /******************************************************************
   * Now, cluster the events which were found to be coincident
   * (within the allowed window) with the injections.  Note, here it 
   * is natural to cluster with a time window  
   *		    clust = 2 * dt
   ******************************************************************/

  /* cluster the events */
  if ( cluster && inspiralEventList)
  {
    LAL_CALL( LALClusterSnglInspiralTable( &stat, inspiralEventList,
          clust, clusterchoice), &stat);
  }


  /*****************************************************************
   * open output xml file
   *****************************************************************/
  LAL_CALL( LALOpenLIGOLwXMLFile(&stat, &xmlStream, outfileName), &stat);


  /* write out the process and process params tables */
  LAL_CALL( LALGPSTimeNow ( &stat, &(proctable.processTable->end_time),
        &accuracy ), &stat );
  LAL_CALL( LALBeginLIGOLwXMLTable( &stat, &xmlStream, process_table ), 
      &stat );
  LAL_CALL( LALWriteLIGOLwXMLTable( &stat, &xmlStream, proctable, 
        process_table ), &stat );
  LAL_CALL( LALEndLIGOLwXMLTable ( &stat, &xmlStream ), &stat );
  LALFree( proctable.processTable );

  /* erase the first empty process params entry */
  {
    ProcessParamsTable *emptyPPtable = procparams.processParamsTable;
    procparams.processParamsTable = procparams.processParamsTable->next;
    LALFree( emptyPPtable );
  }

  /* write the process params table */
  LAL_CALL( LALBeginLIGOLwXMLTable( &stat, &xmlStream, 
        process_params_table ),	&stat );
  LAL_CALL( LALWriteLIGOLwXMLTable( &stat, &xmlStream, procparams, 
        process_params_table ), &stat );
  LAL_CALL( LALEndLIGOLwXMLTable ( &stat, &xmlStream ), &stat );
  while( procparams.processParamsTable )
  {
    this_proc_param = procparams.processParamsTable;
    procparams.processParamsTable = this_proc_param->next;
    LALFree( this_proc_param );
  }
  
  /* Write the found injections to the sim table */
  if ( simEventList )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable (&stat, &xmlStream, 
          sim_inspiral_table),&stat);
    simFoundTable.simInspiralTable = simEventList;
    LAL_CALL( LALWriteLIGOLwXMLTable (&stat, &xmlStream, simFoundTable, 
          sim_inspiral_table), &stat);
    LAL_CALL( LALEndLIGOLwXMLTable (&stat, &xmlStream), &stat);
  
     /* free the temporary memory containing the events */
    while (simEventList)
    {
      SimInspiralTable *currentSimEvent;
      currentSimEvent = simEventList;
      simEventList = simEventList->next;
      LALFree( currentSimEvent );
    }
  }
  
  /* Write the results to the inspiral table */
  if ( inspiralEventList )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable (&stat, &xmlStream, 
          sngl_inspiral_table),&stat);
    myTable.snglInspiralTable = inspiralEventList;
    LAL_CALL( LALWriteLIGOLwXMLTable (&stat, &xmlStream, myTable, 
          sngl_inspiral_table), &stat);
    LAL_CALL( LALEndLIGOLwXMLTable (&stat, &xmlStream), &stat);
  }

  /* free the temporary memory containing the events */
  while (inspiralEventList)
  {
    SnglInspiralTable *currentInspiralEvent;
    currentInspiralEvent = inspiralEventList;
    inspiralEventList = inspiralEventList->next;
    LALFree( currentInspiralEvent );
  }

  /* close the output file */
  LAL_CALL( LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);


  /* Open the missedInjectionsFile and copy the missed injections into it */
  if(saveMissedInjections == TRUE)
  {    
    /* Write the missed injections to the sim table */
    LAL_CALL( LALOpenLIGOLwXMLFile(&stat, &xmlStream, 
          missedInjectionsFile), &stat);
    if ( simEventMissed )
    {
      LAL_CALL( LALBeginLIGOLwXMLTable (&stat, &xmlStream, 
            sim_inspiral_table),&stat);
      simMissedTable.simInspiralTable = simEventMissed;
      LAL_CALL( LALWriteLIGOLwXMLTable (&stat, &xmlStream, simMissedTable, 
            sim_inspiral_table), &stat);
      LAL_CALL( LALEndLIGOLwXMLTable (&stat, &xmlStream), &stat);
    }
    LAL_CALL( LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);
  }

  /* free the temporary memory containing the missed injections */
  while (simEventMissed)
  {
    SimInspiralTable *currentSimEventMissed;
    currentSimEventMissed = simEventMissed;
    simEventMissed = simEventMissed->next;
    LALFree( currentSimEventMissed );
  }

  LALCheckMemoryLeaks();
  return 0;
}

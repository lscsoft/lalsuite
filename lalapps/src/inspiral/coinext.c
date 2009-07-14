/*
*  Copyright (C) 2007 Alexander Dietz
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


/******************************
include
*******************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> 
#include <errno.h>
#include <signal.h>
#include <getopt.h>
#include <time.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/Date.h>
#include <lal/lalGitID.h>
#include <lalappsGitID.h>

#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "coinExt"


#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
  this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
  calloc( 1, sizeof(ProcessParamsTable) ); \
  LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
      PROGRAM_NAME ); \
  LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
      long_options[option_index].name ); \
  LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
  LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );

/*******************************
  structures
********************************/

typedef struct S_InspiralList {
  char  filename[256];
  long  timeStart;
  long  timeLength;
  int   coin;
  int   dagNumber;

  LIGOTimeGPS startTime;
  LIGOTimeGPS endTime;
} InspiralList;

typedef struct S_ExternalList {
  long             gps;
  long             gpsNano;
  int              status[3];
  int              dagStatus[3];

  ExtTriggerTable* table;
} ExternalList;

int checkCGNtrigger( void );
int readCGN( void );
void writeCGN( void );	
void checkInspiralTrigger( void );
int coinInspiralCGN( void );
int createInspiralTable( int lifo, int dag );
int writeXML( int ifo, int cgnIndex, int dag );
void startAnalysisJob(ExternalList* eList, int ifo);
char* getIFOname(int index);
int checkGPS( long gps, long gpsNano );
int getPeriod( long gps );
char* getCurrentTime( void );
void generateDirName( void );
int checkProxy( void );
int ce_startup( void );
int ce_shutdown( void );
int printOut(int mode, char* text);
int ce_wait(int minutes);
void ce_restart(void);
char* basename1( char* name);
void sigint_handler(int sig);
int arg_parse_check( int argc, char *argv[]);

/******************************
   variables
*******************************/
  
LALLeapSecAccuracy accuracy=LALLEAPSEC_LOOSE;

int flagTest=0;    /* flag indicating a test of the program */
int flagRestart;   /* flag indicating to restart the complete DAQ */
int flagRecalc=0;  /* flag indicating to redo the inspiral analysis */
int flagTriggerFile=0; /* flag indicating to use special trigger file */
int flagVerbose=0; /* flag for debugging the code */
int flagDirData=0; /* flag indicating that the data directory is explicit specified by the user */

/*list of current available Inspiral XML trigger files
  index  0: L1 |  1: H1 | 2: H2 */
int numberInspiralList[3];
InspiralList inspiralList[3][1000];

/* list of all CGN stored locally */
int numberStoredCGN=0;
ExternalList extList[2000];

/* new CGN read from data */
int numberNewCGN=0;
ExternalList extNewList[2000];

/* directories and filenames */
static char nameLogfile[256]="";
static char nameSaveInspiralList[256] = "saveInspiral.cel";
static char nameInputList[256]="";
static char nameSaveExternalList[256] = "";
static char dirData[256] ="";
static char dirJobs[256] ="";
static char dirInspiral[256]          = "/data2/dietz/S4/CoinExt/cluster";
static char destCopyURL[]          = "gsiftp://dietz.phys.lsu.edu/home/alex/Work2/E_ExtTrigger/Online/Data/";
static char serverName[50]           = "ldas.ligo-la.caltech.edu";
static char calibrationVersion[10]="V1";

/*  web page with avaliable triggers */
static char webPage[256] = "http://www.uoregon.edu/~ileonor/ligo/s5/grb/online/currentgrbs_html_S5.txt";

/* where does the prog run? LLO or LHO */
static char ifo[] = "LLO";
static char instrName[3][3]={"L1","H1", "H2"};
static char ifoName[3][2]={"L","H","H"};

/* times and name of periods (extentable!) */
int numberPeriods=3;
long* startPeriod;
long* endPeriod;
static char namePeriod[3][3]={"S4","A4","S5"};
static char run[10]="";

static char dagName[3][2]={"A","B"};
int startDAG=0;
int endDAG=1;

/* file for logfile */
FILE* logfile;

/* some other things... */
double fDiffPlus=0;
double fDiffMinus=0;
int nNumberPlus=0;
int nNumberMinus=0;
int nNumberTotal=0;
double fTimeTotal=0;

/* some boundary parameters */
double snrCut=1.0;
float condorTresholdSNR=6.0;
int condorMaxJobs = 40;
int condorNumBanks = 80;
long timeWindow=3600;

/* standard waiting time in minutes */
int waitingTime=180;

/* loop boundaries */
int startIFO=0;
int endIFO=0;

/* output stream */
LIGOLwXMLStream  xmlStream;

/* some global tables */
SnglInspiralTable*    list = NULL;
ExtTriggerTable*      headExt  = NULL;
ExtTriggerTable*      thisExt  = NULL;
ExtTriggerTable*      thisNewExt[1000]; /* vector for 100 new externals maximum */

MetadataTable    proctable;
MetadataTable    processParamsTable;
ProcessParamsTable* this_proc_param=NULL;  

/* tables to use for summary (times of each used xml file) */
SummValueTable* listSumm=NULL;  
MetadataTable   listSummTable;

char message[256];

/* LALStatus */
static LALStatus status;

/*------------------------------
  main:
  main program part: reading program agruments
  and start the 'endless' loop 
  ------------------------------ */
int main(int argc, char* argv[])
{
  int num;
  char* env;

  /* create the process and process params tables */  
  proctable.processTable = (ProcessTable *) calloc( 1, sizeof(ProcessTable) );
  this_proc_param = processParamsTable.processParamsTable =  
    (ProcessParamsTable *) calloc( 1, sizeof(ProcessParamsTable) );

  /* call the argument parse and check function */
  arg_parse_check( argc, argv);

  /* setting all the paths... */
  env=getenv("PWD");
  sprintf(message, "%s/../Data/%s/", env,run);
  strcpy(dirData, message);
  sprintf( nameLogfile, "coinExt%s.log",run);
  sprintf( nameSaveExternalList, "externalTrigger%s.cel",run);
  if ( !flagTriggerFile) {
    sprintf( nameInputList, "triggers%s.list", run);
  }
  sprintf( destCopyURL, "gsiftp://dietz.phys.lsu.edu/home/alex/Work2/E_ExtTrigger/Online/Data/%s/",run );
  if ( !strcmp(run,"A4") ) {
    sprintf( webPage, "http://www.uoregon.edu/~ileonor/ligo/a4/grb/online/currentgrbs_html_A4.txt");
  }

  /*  startup */
  ce_startup();
  if (signal(SIGINT, sigint_handler) == SIG_ERR) {
    perror("signal");
    exit(1);
  }	
 
  /* check if to restart all */
  if (flagRestart) {
    sprintf(message,"RESTARTING program. Deleting all previous data files...");
    printOut(3, message);
    ce_restart();
  }

  if (!flagRestart) {
    /* reading information from CGN trigger that are already available */
    num=readCGN();
    sprintf(message, "Reading %d CGN triggers from file '%s'",num, nameSaveExternalList);
    printOut(3, message);
  }

  /* setting some loop parameters */
  if (strcmp (ifo,"LLO") == 0) {
    startIFO=0;
    endIFO=0;
  } else {
    startIFO=1;
    endIFO=2;
  }

  if (flagTest) {
    printf("*** TEST mode activated ***\n");
    sprintf(nameInputList,"trigger.list");
    waitingTime=1;
  }
  /* printing some user parameters */  
  printf("--- nameLogFile: %s\n",nameLogfile );
  printf("--- waitingTime: %d minutes \n",waitingTime);
  printf("--- dirInspiral: %s\n",dirInspiral);
  printf("--- dirData:       %s\n",dirData);
  printf("--- nameInputList: %s\n",nameInputList);
  printf("--- nameSaveList : %s\n",nameSaveExternalList);
  printf("--- running at site %s\n",ifo);
  printf("--- run %s\n",run);
  printf("--- calibration version:%s\n",calibrationVersion);
  printf("--- data-server used %s\n",serverName);
  printf("--- webPage %s\n",webPage);
  printf("--- desCopyURL %s\n", destCopyURL);
  printf("--- timeWindow: %ld seconds\n",timeWindow);
  if (flagTriggerFile) {
    printf("--- using trigger file %s\n",nameInputList);
  }
  if (flagRecalc) {
    printf("--- starting analysis dag for every GRB\n");
  } 
  printf("\n");

  /* start 'endless' loop */
  do {
    /* download latest list of external trigger */
    if (checkCGNtrigger()) {
      
      /* store all CGN triggers to file */
      writeCGN();
    }
     
    /* check for new XML files */
    checkInspiralTrigger();
    
    /* check coincidents with external trigger */
    coinInspiralCGN();

    /* waiting for some time... */
    ce_wait(waitingTime);

  } while(1);

}

/******************************
checkCGNtrigger
  downloading the latest CGN trigger from Isabels list
*******************************/
int checkCGNtrigger( void )
{
  /* specify used variables here ... */
  int i,nifo;
  int counter=0;
  int numEvents=0;
  int index=0;
  int newTrigger=0;
  char command[256];
  char input[256];
  char inputGPS[256];
  char input1[256];
  char input2[256];
  LIGOTimeGPS* gps;
  struct tm date;
  FILE* fileTrigger;

  /* initialize date */
  gps=(LIGOTimeGPS*)LALCalloc(1, sizeof(LIGOTimeGPS));

  /* check if to get some trigger file or not */
  if (!flagTriggerFile) {
    
    /* copy test file if it is a test */
    if (flagTest) {
      system( "cp Test/trigger.list ." );
      
      /* downloading file */
    } else {
      sprintf(command, "wget %s -O %s > download.log 2>&1",
	      webPage, nameInputList);
      system(command);
    }
  }

  /* delete something, (only neccessary for C code) */
  while ( headExt )  {
    thisExt = headExt;
    headExt = headExt->next;
    LALFree( thisExt );
  }
  
  /* open the trigger file */
  fileTrigger=fopen(nameInputList, "r");
  
  /* reading GCN number */
  fscanf(fileTrigger, "%s", input);

  do {
    
    /* allocate new memory for next table */
    if ( ! headExt ) {
      thisExt = headExt = (ExtTriggerTable *) LALCalloc( 1, sizeof(ExtTriggerTable) );
    } else {
      thisExt = thisExt->next = (ExtTriggerTable *) LALCalloc( 1, sizeof(ExtTriggerTable) );
    }
    
    thisExt->event_number_gcn=atoi(input);

    /* reading GRB number  */
    fscanf(fileTrigger, "%s", input);
    strncpy(thisExt->event_number_grb, input,8);

    /* reading UT time -- not stored */
    fscanf(fileTrigger, "%s", input);    

    /* reading  GPS time  */
    fscanf(fileTrigger, "%s", inputGPS);
    strncpy(input1, inputGPS, 9);
    strncpy(input2, inputGPS+9, strlen(inputGPS)-10);
    thisExt->start_time=atol(input1);    
    thisExt->start_time_ns=(long)(atof(input2)*1.0e+9);

    /*  reading RA */
    fscanf(fileTrigger, "%s", input);
    thisExt->event_ra=atof(input);

    /* reading DEC */
    fscanf(fileTrigger, "%s", input);
    thisExt->event_dec=atof(input);

    /* reading fave LHO */
    fscanf(fileTrigger, "%s", input);
    thisExt->ligo_fave_lho=atof(input);

    /* reading  fave LLO */
    fscanf(fileTrigger, "%s", input);
    thisExt->ligo_fave_llo=atof(input);   

    /* reading    */
    fscanf(fileTrigger, "%s", input);
    thisExt->ligo_delay=atof(input);

    /* reading GCN number */
    fscanf(fileTrigger, "%s", input);

    /* check end of file */
    numEvents++;
  } while (!feof(fileTrigger));

  if (numEvents<0) {
    sprintf(message, "Error: Unable to read data from file %s",nameInputList);
  } else {
    
    /* loop over all entries */
    for ( thisExt=headExt; thisExt; thisExt=thisExt->next ) {

      /* check if there is a new CGN */
      index=checkGPS( thisExt->start_time, thisExt->start_time_ns );
      counter++;

      if (index==-1) {

	/* store new trigger in external triggers structure */
	extList[numberStoredCGN].gps=thisExt->start_time;
	extList[numberStoredCGN].gpsNano=thisExt->start_time_ns;
	extList[numberStoredCGN].status[0]=0;
	extList[numberStoredCGN].status[1]=0;
	extList[numberStoredCGN].status[2]=0;
	extList[numberStoredCGN].dagStatus[0]=0;
	extList[numberStoredCGN].dagStatus[1]=0;
	extList[numberStoredCGN].dagStatus[2]=0;
	extList[numberStoredCGN].table=(ExtTriggerTable*)LALCalloc( 1, sizeof(ExtTriggerTable) );
	memcpy( extList[numberStoredCGN].table, thisExt, sizeof(ExtTriggerTable));

	/* getting user friendly date and time */
	gps->gpsSeconds=thisExt->start_time;
	gps->gpsNanoSeconds=thisExt->start_time_ns;
	XLALGPSToUTC(&date, &gps);

	/* generate output saying that there is a new CGN trigger */
	sprintf(message,"New CGN trigger occured at %d-%02d-%02d %d:%02d:%02d (%d)",
		date.tm_year+1900, date.tm_mon+1,
		date.tm_mday,date.tm_hour,
		date.tm_min, date.tm_sec, 
		gps->gpsSeconds);
	printOut(3, message);

	/* start special analysis job if data should be recalculated */
	if (flagRecalc) {
	  for ( nifo=startIFO;nifo<=endIFO;nifo++) {
	    startAnalysisJob( &extList[numberStoredCGN], nifo );
	  }
	}

	/* increase number of stored triggers from the Global GRB Network */
        numberStoredCGN++;

      }  else {

	/* old trigger found -> set pointer to correct table */
	extList[index].table=thisExt;		
      }
    }
  }
  
  /* check if there is a new E.T. */
  newTrigger=0;
  for (nifo=startIFO;nifo<=endIFO;nifo++) {
    for (i=0;i<numberStoredCGN;i++) {
      if (extList[i].status[nifo]==0) {
	
	/* set flag and exit loop */
	newTrigger=1;
	continue;
      }
    }
  }
  
  /* return flag if a new E.T. is recorded */
  return newTrigger;
 
}

/******************************
readCGN:
  reading data for stored CGN (only once at start of program)
  executed in the 'startup' method 'externalTriggers.cel'
*******************************/
int readCGN( void )
{ 
  FILE* file;
  int c=0;
  long gps;
  long gpsNano;
  int status0, status1, status2;
  int dagStatus0, dagStatus1, dagStatus2;


  /* open file for reading data */
  file=fopen( nameSaveExternalList , "r" );
  
  /* check file */
  if (!file) {
    sprintf(message,"File %s does not exist",nameSaveExternalList);
    printOut(1, message);
    
    numberStoredCGN=0;
    return numberStoredCGN;
  }

  /* loop over all entries in the file */
  do {

    /* read next line and store data */
    fscanf(file, "%ld,%ld,%d,%d,%d, %d, %d, %d", &gps, &gpsNano, &status0, &status1, &status2,
	   &dagStatus0, &dagStatus1, &dagStatus2); 
    extList[c].gps      = gps;
    extList[c].gpsNano  = gpsNano;
    extList[c].status[0]= status0;
    extList[c].status[1]= status1;
    extList[c].status[2]= status2;
    extList[c].dagStatus[0]= dagStatus0;
    extList[c].dagStatus[1]= dagStatus1;
    extList[c].dagStatus[2]= dagStatus2;
    extList[c].table    = NULL;
    
    /* increase counter */
    c++;

  } while ( !feof(file));
  fclose(file);

  /* store counter and return value */
  numberStoredCGN=c-1;
  return numberStoredCGN;
}


/********************************
writeCGN:
  writing complete list of CGN triggers to a file
  This file is just for lookup, it is not needed in this code
  REMARK: Not very elegant, but easy. All previos data in the file are overwritten
          by the same data and eventually some more entries are added
*********************************/
void writeCGN( void )
{	
  FILE* file;
  int i;

  /* open output file for overwriting all data! */
  file=fopen(nameSaveExternalList,"w");
  for ( i=0; i<numberStoredCGN; i++ ) {    
    fprintf( file, "%ld,%ld,%d,%d,%d, %d, %d, %d\n",
	     extList[i].gps, extList[i].gpsNano, 
	     extList[i].status[0], extList[i].status[1], extList[i].status[2],
	     extList[i].dagStatus[0], extList[i].dagStatus[1], extList[i].dagStatus[2]);
  }
  fclose(file);

}

/*******************************
checkInspiralTrigger:
  create a list of all trigger-xml files that has been recopied from the finished DAG
  If new data are present, a 1 is returned, 0 else.
  REMARK: Again not very nifty. Always reads the previous data
********************************/
void checkInspiralTrigger( void )
{
  int i,j,dag, c=0;
  char dummy[256];
  char command[256];
  char text[256];
  char* filename;
  char savename[256];
  FILE* file;
  FILE* out;

  /* loop over IFO's to take into account */
  for ( i=startIFO;i<=endIFO;i++) {

    c=0;
    for ( dag=startDAG; dag<=endDAG; dag++) {

      /* create new file containing all present inspiral trigger lists */
      if (flagRecalc) {

	/* use directory for recalculated trigger */
	sprintf( command, "ls -1 ../Trigger/%s/%s/%s/*.xml > namesInspiral.cel", run, instrName[i], dagName[dag]);
      } else{

	/* use directory for online created triggers */
	sprintf( command,"ls -1 %s%s/*.xml > namesInspiral.cel", dirInspiral, instrName[i] );
      }
      system( command );
    
      /* read list and extract time data */
      file=fopen( "namesInspiral.cel","r" );
      do {      

	/* read next filename and remove directory path from it */
	fscanf(file, "%s", text); 
	filename = basename1( (char*)text );

	/* extracting relevant parts of the filename
	   assuming the length are always the same! */
	strcpy(inspiralList[i][c].filename, filename);
	strncpy(dummy, filename+8, 9);dummy[9]='\0';
	inspiralList[i][c].timeStart=atoi(dummy);
	strncpy(dummy, filename+18, 4);dummy[4]='\0';
	inspiralList[i][c].timeLength=atoi(dummy);
	inspiralList[i][c].coin=0;
	inspiralList[i][c].dagNumber=dag;

	if (flagVerbose) {
	  sprintf( message,"Found inspiral-trigger, entry %d: %s| GPS: %ld | dag: %d", 
		   c, filename, inspiralList[i][c].timeStart, dag );
	  printOut( 7, message );
	}
	
	/* increasing counter */
	c++;
	
      } while ( !feof(file) );
      fclose(file);
    }
    
    /* storing total number of inspiral xml files */
    numberInspiralList[i]=c-1;
    
    /* open file and store all entries of the inspiral list */
    sprintf(savename, "%s%s",ifo, nameSaveInspiralList);    
    if (flagVerbose) {
      sprintf( message,"Storing files-infos into file %s", savename);
      printOut( 7, message );
    }
    out=fopen( savename, "w" );
    for (j=0;j<numberInspiralList[i]; j++) {
      fprintf(  out, "%s %ld %ld\n",inspiralList[i][j].filename,
		inspiralList[i][j].timeStart,inspiralList[i][j].timeLength);
    }
    fclose(out);
    
  } 

}

/********************************
coinInspiralCGN:
  search for coincidences of inspiral and CGN triggers
  and store data to new xml file
*********************************/
int coinInspiralCGN( void )
{
  /* variables */
  int cgn, lifo, insp, dag;  
  int numberFound=0;
  int copyStatus;
  double timeCGN;
  SnglInspiralTable* element;
  ExtTriggerTable* table;
  long timeStart, timeEnd;

  /* loop over ALL CGN triggers  */
  for ( cgn=0; cgn<numberStoredCGN; cgn++ ) {
    
    /* loop over IFO's */
    for (lifo=startIFO;lifo<=endIFO;lifo++) { 

      /* check status of the current CGN triggers */
      if (extList[cgn].status[lifo]<4) {

	/* get pointer to corresponding E.T. structure */
	table=extList[cgn].table;
	if ( !table ) {
	  sprintf(message,"WARNING in coinInspiralCGN: A ExtTriggerTable does not exist for cgn %d. " 
		  "Check program!\n",cgn);
	  printOut( 1, message );
	  continue;
	}
	
	/* get time of CGN trigger */
	timeCGN=(double)table->start_time + (double)table->start_time_ns/1.0e+9;

	/* loop over the different DAG's */
	for ( dag=startDAG; dag<=endDAG; dag++ ){
			
	  /* loop over all available trigger-XML-files  */
	  numberFound=0;
	  for (insp=0; insp<numberInspiralList[lifo]; insp++ ) {
	    
	    inspiralList[lifo][insp].coin=0;
	    
	    /* get starting end ending time of the time period
	       the inspiral triggers are from this file */ 
	    timeStart=inspiralList[lifo][insp].timeStart;
	    timeEnd=timeStart+inspiralList[lifo][insp].timeLength;
	    /* */
	    if ( (timeStart>=timeCGN-timeWindow && timeStart<=timeCGN+timeWindow) || 
		 (timeEnd>=timeCGN-timeWindow && timeEnd<=timeCGN+timeWindow)) {
	      printf("overlap, dag:%d  dagNumber:%d\n",dag,inspiralList[lifo][insp].dagNumber );
	    }
	    /* */

	    /* check if time overlaps */ 
	    if ( ( (timeStart>=timeCGN-timeWindow && timeStart<=timeCGN+timeWindow) || 
		   (timeEnd>=timeCGN-timeWindow && timeEnd<=timeCGN+timeWindow) ) &&
		 (inspiralList[lifo][insp].dagNumber==dag) ) {
	    
	      /* mark corresponding inspiral file */
	      numberFound++;
	      inspiralList[lifo][insp].coin=1;
	    }
	  }	
	
	  /* happy: coincidence trigger found!!! 
	     Create xml file and copy data for end-analysis !! */
	  if (numberFound) {
	    
	    /* output */
	    sprintf(message, "ExtTrigger found in coincidence at time %f", timeCGN);
	    printOut(4, message);

	    /* creating final xml file and copy it to somewhere */
	    createInspiralTable( lifo, dag );
	    copyStatus=writeXML( lifo, cgn, dag );

	    /* set new status */
	    if (copyStatus==0) {
	      extList[cgn].dagStatus[lifo]+=pow(3, dag);
	      sprintf(message, "coinInspiralCGN: DAG status for ifo %d (CGN=%d) set to %d\n", 
		      lifo, cgn, extList[cgn].dagStatus[lifo]);
		if (flagVerbose) printOut(7, message);
	    } else {
	      sprintf(message, "Unable to transfer data to %s. Trying later again...", destCopyURL);
	      printOut(2,message);
	    }

	  }

	} /* {dag} */

	/* set new status */
	if (extList[cgn].dagStatus[lifo]>=pow(3, endDAG+1)-1) {
	  extList[cgn].status[lifo]=5;
	  sprintf(message, "Status for CGN %d set to 5 (all DAG's ready and processed)\n", cgn);
	  if (flagVerbose) printOut(7,message);
	}
	
	/* delete old inspiral triggers */
	while ( list )  {
	  element = list;
	  list = list->next;
	  LALFree( element );
	}  
      }
      
    }
  }
  
  return 0;
}

/********************************
createInspiralTable:
  create table with all needed inspiral triggers in it
*********************************/
int createInspiralTable( int lifo, int dag )
{
  /* variables */
  int insp,numTrig, numSum,flag;
  int numberDiscard=0;
  char filename[256]; 
  char comment[256]; 
  SnglInspiralTable* head;
  SnglInspiralTable* tail;
  SnglInspiralTable* mainPointer;
  SnglInspiralTable* prevElement=NULL;
  SnglInspiralTable* element;
  SummValueTable* summary;
  SummValueTable* tailSumm=NULL;
  tail=NULL;

  /* loop over all inspiral files */
  flag=0;
  for (insp=0; insp<numberInspiralList[lifo]; insp++ ) {

    /* check if this entry if marked */
    if (inspiralList[lifo][insp].coin) {

      /* read data from the XML file using LAL functions 
	 into new allocated memory, pointed by 'head' */
      head=NULL; 

      if (flagRecalc) {
	/* use directory for recalculated trigger */
	sprintf( filename,  "../Trigger/%s/%s/%s/%s", run, getIFOname(lifo), dagName[dag], 
		 inspiralList[lifo][insp].filename);
      } else{
	/* use directory for online created triggers */
	sprintf( filename, " %s%s/%s",
		 dirInspiral, getIFOname(lifo), inspiralList[lifo][insp].filename);
      }

      if (flagVerbose) {
	sprintf(comment, "Trying opening file %s, lifo: %d",filename, lifo);
	printOut(7, comment);
      }

      numTrig=LALSnglInspiralTableFromLIGOLw(&head, filename, 1,100000);  

      /* read search summary table from same file 
	 and store values into inspiral-list */
      summary=NULL;
      numSum=SummValueTableFromLIGOLw(&summary, filename);
      inspiralList[lifo][insp].startTime=summary->start_time;
      inspiralList[lifo][insp].endTime=summary->end_time;

      /* apply SNR cut */
      if (numTrig>0) {
	
	/* some basic pointer settings */
	numberDiscard=0;
	mainPointer=head;
	prevElement=NULL;
	while ( mainPointer )  {	 

	  /* check for SNR cut */
	  if (mainPointer->snr < snrCut) {
	    /* discard one element */
	    if (prevElement) {
	      prevElement->next=mainPointer->next; 
	    } else {
	      head=head->next;  /* the head element is changing */
	    }
	    element=mainPointer;  /* store in dummy struct*/
	    mainPointer=mainPointer->next;  /* go one up */
	    LALFree( element );	  /* free dummy here */
	    numberDiscard++;

	  } else {
	    /* do not discard */
	    mainPointer=mainPointer->next; /* go one up */
	    if (prevElement) {
	      prevElement=prevElement->next; /* go one up */
	    } else {
	      prevElement=head; /* or specify */
	    }

	  }
	}

      }

      /* store pointer to data only when there are some data */
      if (head) {

	if (tail)
	  tail->next=head;
	else
	  list=head;

	/* store last element */
	tail=prevElement;

	if (!prevElement) {
	  printf("Hopefully this line will never be printed to screen....");
	  exit(2);
	}
      }

      /* store pointer to FIRST table only - in list */
      if (!flag) {
	flag=1;

	/* set summ values in table */
	listSumm=summary;
	tailSumm=listSumm;
      } else {

	/* set summ values in table */
	tailSumm->next=summary;
	tailSumm=tailSumm->next;
      }

      /* setting values in the summ-value table */
      LALSnprintf( tailSumm->name, LIGOMETA_SUMMVALUE_NAME_MAX, "length" ); 
      tailSumm->value=summary->end_time.gpsSeconds - summary->start_time.gpsSeconds;
      LALSnprintf( tailSumm->comment, LIGOMETA_SUMMVALUE_COMM_MAX, "%s",inspiralList[lifo][insp].filename ); 

      if (flagVerbose) {
	sprintf(comment, "+++ Adding %d-%d triggers from file %s ", 
		numTrig, numberDiscard,inspiralList[lifo][insp].filename);  
	printOut( 7, comment);
      }
    }
  }

  return 0;
}


/********************************
writeXML:
  writing XML file with all needed informations
*********************************/
int writeXML( int nifo, int cgnIndex, int dag )
{
  /* variables */
  char command[256];
  char filename[256];
  char filenameOut[256];
  MetadataTable outputTable; 
  SnglInspiralTable*    thisInsp = NULL;
  int exitStatus;
 
  /* create filename for output file */
  if (flagVerbose) {
    sprintf(message, "writeXML: GPS time: %9ld ", extList[cgnIndex].gps);
    printOut( 7,message);
  }
  sprintf( filename, "externalTrig%9ld%s-%s.xml", extList[cgnIndex].gps, dagName[dag], getIFOname(nifo));
  sprintf( filenameOut, "%s%s", dirData, filename);
  sprintf( message, "Creating output file %s", filenameOut );
  printOut(5, message);

  /* allocate memory for output stream and open file */
  memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
  LALOpenLIGOLwXMLFile( &status, &xmlStream, filenameOut );
	    
  /* write process table */
  LALGPSTimeNow ( &status, &(proctable.processTable->end_time), &accuracy );
  LALBeginLIGOLwXMLTable( &status, &xmlStream, process_table );
  LALWriteLIGOLwXMLTable( &status, &xmlStream, proctable, process_table );
  LALEndLIGOLwXMLTable ( &status, &xmlStream );
  
  /* write process_params table */
  LALBeginLIGOLwXMLTable( &status, &xmlStream,process_params_table );
  LALWriteLIGOLwXMLTable( &status, &xmlStream, processParamsTable, process_params_table );
  LALEndLIGOLwXMLTable ( &status, &xmlStream );
  
  /* write summ_value table */
  listSummTable.summValueTable=listSumm;
  LALBeginLIGOLwXMLTable( &status, &xmlStream, summ_value_table );
  LALWriteLIGOLwXMLTable( &status, &xmlStream, listSummTable, summ_value_table );
  LALEndLIGOLwXMLTable ( &status, &xmlStream );
  
  /* write external_triggers table */
  /* set pointer to the correct table */
  outputTable.extTriggerTable = extList[cgnIndex].table;
  /* there is no next table (exactly ONE exttrigger entry in the table */
  outputTable.extTriggerTable->next=NULL;

  LALBeginLIGOLwXMLTable( &status, &xmlStream, ext_triggers_table );
  LALWriteLIGOLwXMLTable( &status, &xmlStream, outputTable, ext_triggers_table );
  LALEndLIGOLwXMLTable( &status, &xmlStream );	

  /* writing inspiral table to xml file */
  outputTable.snglInspiralTable = list;
  LALBeginLIGOLwXMLTable( &status, &xmlStream, sngl_inspiral_table );
  LALWriteLIGOLwXMLTable( &status, &xmlStream, outputTable, sngl_inspiral_table );
  LALEndLIGOLwXMLTable( &status, &xmlStream );

  /* closing XML file */
  LALCloseLIGOLwXMLFile( &status, &xmlStream );

  /* remove the 'list'-table */
  while ( list ) {
    thisInsp = list;
    list = list->next;
    LALFree( thisInsp );
  }
  LALCheckMemoryLeaks();
  list=NULL;

  /* transfer the files to the master computer 'adietz' */
  sprintf(command ,"globus-url-copy file://%s %s%s", filenameOut, destCopyURL, filename);
  if (flagVerbose) {
    sprintf( message,"Executing command: %s",command);
    printOut( 7, message);
  }
  exitStatus=system(command);
  if (flagVerbose) {
    sprintf( message,"exitStatus: %d",exitStatus);
    printOut(7,message);
  }

  return exitStatus;
}

/****************************** 
startAnalysisJob
   creates and starts a special analysis job
********************************/
void startAnalysisJob(ExternalList* eList, int nifo)
{
  int indexPeriod;
  LIGOTimeGPS gpsnow;
  int numberOfSegments=0;
  int minSegmentLength=512;
  int pad_data=8.0;
  int dag;
  char dirJob[256];
  /*char sedcal[256];*/
  char command[256];  
  FILE* fileSegment;
  FILE* fileSED; 
  double d=0, sumTime=0, maxSegment=0, minSegment=-1.0;
  char dummy[256];
  time_t rawtime;

  /* discard any times before S4 */
  if (eList->gps<793130413) {
    sprintf(message,"GPS time %9ld lies before run S4! Discarded...",eList->gps);
    printOut(1, message);
    return;
  }

  /* check user proxy */
  if ( !checkProxy() ) {
    sprintf(message,"No valid user grid-proxy found. Unable to proceed...");
    printOut(2, message);
    return;
  }

  /* check if enough time to create the frames (40 minutes) */
  LALGPSTimeNow( &status, &gpsnow, LALLEAPSEC_LOOSE );

  
  time ( &rawtime );
  gpsnow.gpsSeconds=rawtime-315964787;
  if ( gpsnow.gpsSeconds - eList->gps < 2400 ) {
    sprintf(message,"GRB trigger occured only %ld seconds before, analysing event in next round...", 
	    gpsnow.gpsSeconds - eList->gps);
    printOut(3, message);
    return;
  }

  /* get corresponding time period  */
  indexPeriod=getPeriod( eList->gps );
  if (indexPeriod==-1) return;

  /* find possible segments that can be used */
  if (indexPeriod==1) {
    /* user EVERY available data */
    sprintf(command, "LSCdataFind --server %s --observatory %s --gps-start-time %9ld --gps-end-time %9ld --type RDS_R_L3 2> .notfound.cel | xargs FrFileRanges 2> .notfound.cel | gawk '{s++; print s, $1, $2, $2-$1;}' > segment.txt 2> .notfound.cel ",
	  serverName, ifoName[nifo], eList->gps-timeWindow, eList->gps+timeWindow );
  } else {
    /* use only Science segments and comparable types */
    sprintf(command, "LSCsegFind --server=%s --interferometer %s --type Science --output-format segwizard --gps-start-time %9ld --gps-end-time %9ld > segment.txt 2> .notfound.cel ",
	    serverName, instrName[nifo], eList->gps-timeWindow, eList->gps+timeWindow );
  }
  /*sprintf(command, "LSCdataFind --server %s --observatory %s --gps-start-time %9ld --gps-end-time %9ld --type RDS_R_L3 | xargs FrFileRanges | gawk '{s++; print s, $1, $2, $2-$1;}' > segment.txt 2> .notfound.cel ",
    serverName, instrName[nifo], eList->gps-timeWindow, eList->gps+timeWindow );*/

  if (flagVerbose) printOut(7, command);
  system (command);

  /* check file segment.txt */
  fileSegment=fopen("segment.txt", "r");
  do {    
    fscanf(fileSegment, "%s",dummy);
    fscanf(fileSegment, "%s",dummy);
    fscanf(fileSegment, "%s",dummy);
    fscanf(fileSegment, "%s",dummy);
    d=atof(dummy);
    
    /* get sum of times and maximum segment (has to be > 2048 seconds) */
    if (!feof(fileSegment)) {
      if (d>minSegmentLength+2*pad_data) {
	if (d>maxSegment) maxSegment=d;
	if (minSegment<0 || d<minSegment) minSegment=d;
	sumTime=sumTime+d;
      } 

    }
  } while (!feof(fileSegment));


  /* use a minimum segment length of 256 seconds! Corresponds to 3 segments. */
  if (minSegment>=minSegmentLength) {

    if (flagVerbose) {
      sprintf( message, "minimum segment: %f | maximum segment: %f | totaltime: %f",minSegment, maxSegment, sumTime);
      printOut(7, message);
    }

    /* calculate some basic parameters */
    /* calculate the number of segments */
    numberOfSegments = (int)((minSegment-2*pad_data)/128.0) -1;
    if (numberOfSegments>15)  numberOfSegments=15;
    
    sprintf(message, "Longest segment found is %.0f s, using %d number of segments in lalapps_inspiral.",
	    maxSegment, numberOfSegments);
    printOut( 3, message);
    
    /* EXTRA FOR BEGINNING OF S5 */
    /*sprintf( sedcal, "A4");
    printf("##########################################\n");
    printf("#### USING CALIBRATION A4 INSTEAD S5 !! ##\n");
    printf("##########################################\n");*/
    /* END EXTRA */      

    /* loop over the different DAG's (low-mass, high-mass */
    for ( dag=startDAG; dag<=endDAG; dag++ ) {

      /* create name of new directory for job execution */
      sprintf(dirJob,"../OnlineAnalysis/Jobs/%s/%s%s-%9ld%s", dirJobs, instrName[nifo], 
	      namePeriod[indexPeriod], eList->gps, dagName[dag]);
    
      /* create directory */
      sprintf(command,"mkdir -p %s",dirJob);
      if (flagVerbose) printOut(7, command);
      system(command);
    
      /* copy the segment file to that directory */
      sprintf(command, "cp segment.txt %s",dirJob);
      if (flagVerbose) printOut(7, command);
      system("cat segment.txt");
      system(command);
      /*system("rm -f segment.txt");*/

      /* symbolic links to all needed basic files */
      sprintf(command,"cd %s; ln -s ../../../ProgInspiral/Executables/* .",dirJob);
      if (flagVerbose) printOut(7, command);
      system(command);           
      
      /* create sed-file */
      fileSED=fopen("sed.file", "w");
      fprintf(fileSED, "s/sedifo/%s/\n",instrName[nifo]);
      fprintf(fileSED, "s/sedrun/%s/\n", namePeriod[indexPeriod]);
      fprintf(fileSED, "s/sedcal/%s/\n", run);
      fprintf(fileSED, "s/sedcav/%s/\n", calibrationVersion);
      fprintf(fileSED, "s/sednumbanks/%d/\n", condorNumBanks);
      fprintf(fileSED, "s/sedtreshsnr/%f/\n", condorTresholdSNR);
      fprintf(fileSED, "s/sedchannel/LSC-DARM_ERR/\n");
      fprintf(fileSED, "s/sednumberofsegments/%d/\n",(int)numberOfSegments);
      fprintf(fileSED, "s/sedservername/%s/\n",serverName);
      fprintf(fileSED, "s/seddag/%s/\n", dagName[dag]);
      fclose(fileSED);

      /* create ini file for job execution */
      sprintf(command,"sed -f sed.file ../OnlineAnalysis/ProgInspiral/online%s.ini > %s/online.ini", 
	      dagName[dag], dirJob);
      if (flagVerbose) printOut(7, command);
      system(command);    

   
      /* create the DAG */ 
      sprintf(command,"cd %s; ./lalapps_inspiral_online_pipe --config-file online.ini --log-path /usr1/dietz",
	      dirJob );
      if (flagVerbose) printOut(7, command);
      system(command);
      
      sprintf(command,"Starting DAG for analysis in directory %s",dirJob);
      printOut(5, command);
      sprintf(command, "cd %s; condor_submit_dag --maxjobs %d online.dag",dirJob, condorMaxJobs);
      if (flagVerbose) printOut(7, command);
      system (command);

      eList->dagStatus[nifo]+=pow(3, dag);
      sprintf(message, "startAnalysisJob: DAG status for ifo %d set to %d\n", nifo, eList->dagStatus[nifo]);
      if (flagVerbose) printOut(7, message);
    }
      
    /* set status */
    eList->status[nifo]=2;

  } else {

    /* set status */
    eList->status[nifo]=4;

    /* output */
    if (d==0) {
      sprintf(message, "No segments found for times %ld to %ld on %s",
	     eList->gps-timeWindow, eList->gps+timeWindow, instrName[nifo] );      
      printOut(3, message);
    } else {
      sprintf(message, 
	      "Max. segment for times %ld to %ld on %s is only %.0f seconds long. No data analyzed.",
	     eList->gps-timeWindow, eList->gps+timeWindow, instrName[nifo], maxSegment );
      printOut(3, message);
    }
  }

}

/****************************** 
getIFOname:
return the IFO name corresponding to 'index'
********************************/
char* getIFOname(int index)
{
  static char ifoname[3];

  switch (index) {
  case 0:
    sprintf(ifoname, "L1");
    break;
  case 1:
    sprintf(ifoname, "H1");
    break;
  case 2:
    sprintf(ifoname, "H2");
  }
  /*printf("ifoname: %s,  index: %d\n",ifoname, index);
    fflush(stdout);*/
  return ifoname;
}

/******************************
checkGPS:
   checking the given GPS time and compare 
   if a CGN trigger at this time exist
*******************************/
int checkGPS( long gps, long gpsNano )
{
  int i;

  /* loop over all stored CGN triggers */
  for ( i=0; i<numberStoredCGN; i++ ) {

    /* compare entries for same time */
    if (gps==extList[i].gps && gpsNano==extList[i].gpsNano) {
      return i; /* returning index */
    }
  }
  
  return -1;
}


/******************************
getPeriod:
   get the index of the period the 
   time lies in
*******************************/
int getPeriod( long gps )
{
  int index=-1;
  int i; 

  /* search for corresponding time period */
  for (i=0;i<numberPeriods;i++) {
    if (gps>startPeriod[i] && gps<endPeriod[i]) {
      index=i;
    }
  }

  /* check if a correct period was found */
  if (index==-1) {
    sprintf(message,"GPS time %9ld lies outside of senseful range!", gps);
    printOut(2, message);
  }

  return index;
}

/*******************************
checkProxy:
   check if a valid proxy is set-up
*******************************/
int checkProxy( void )
{
  char command[256];
  struct stat file;

  sprintf(command, "grid-proxy-info > .proxy 2>.err");
  system(command);
    
  if(!stat(".proxy",&file)) {
    if (file.st_size==0) {
      return 0;
    }
  }

  return 1;
}

/*******************************
getCurrectTime:
   returns the current time as char string
*******************************/
char* getCurrentTime( void )
{
  time_t rawtime;
  struct tm * timeinfo;
  static char currentTime[256];
  
  /* get current time */
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );

  sprintf(currentTime, "%d-%02d-%02d %d:%02d:%02d", timeinfo->tm_year+1900,  timeinfo->tm_mon+1, timeinfo->tm_mday,
	  timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);

  return currentTime;	     
}

/*******************************
generateDirName
   generate unique directory name string from date
*******************************/
void generateDirName( void )
{
  time_t rawtime;
  struct tm * timeinfo;

  /* get current time */
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  sprintf(dirJobs, "Run%d-%02d-%02d_%d-%02d-%02d", timeinfo->tm_year+1900,  timeinfo->tm_mon+1, timeinfo->tm_mday,
          timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);

}


/******************************
ce_startup:
   startup function: looks for existence of the file .coinExt.lock
   if that file exist -> abort. Program running.
   if not exist -> create this file with command stopping this program 
   (sending kill signal 2 to this program, which is fetched by the function sigint_handler)
*******************************/
int ce_startup( void )
{ 
  FILE* out;

  /* checking/creating lock-file */
  out=fopen( ".coinext.lock" ,"r");
  if (out) {
    if (!flagTest) {
      printf("ERROR: Program already running!\n");
      exit(0);
    }
  } else {
    out=fopen( ".coinext.lock", "w" );
    fprintf( out, "kill -2 %d\n", getpid());
    fprintf( out, "rm -f .coinext.lock\n" );
    fclose(out);
  }

  /* open log file for appending */
  logfile=fopen( nameLogfile , "a" );

  /* create startup message */
  sprintf(message, "-------------------------------------------------\n--- Starting coinext at %s",
	  getCurrentTime());
  printOut(6, message);

  startPeriod=(long*)LALCalloc( numberPeriods, sizeof(long) );
  endPeriod=(long*)LALCalloc( numberPeriods, sizeof(long) );
  startPeriod[0]=793130413; /* S4 */
  endPeriod[0]=795679213;
  startPeriod[1]=795813760; /* A4 */
  endPeriod[1]=815961612;
  startPeriod[2]=815155213; /* start S5 in LHO, Nov 4, 2005 9:00 PST = 16:00 UTC */
  endPeriod[2]=999999999;

  /* generate unique directory name for this session*/
  generateDirName();

  /* set correct server name */
  if (strcmp (ifo,"LHO") == 0 || strcmp (ifo,"lho") == 0) {
    sprintf(serverName,"ldas.ligo-wa.caltech.edu");
  } 

  return 1;
}


/******************************
ce_shutdown:
   shutdown function
*******************************/
int ce_shutdown( void )
{

  /* store CGN triggers for next time (status is important!!!) */
  writeCGN();

  /* creating shutdown message */
  sprintf(message, "---"
	  "\n--- Stopping coinext at %s "
	  "\n-------------------------------------------------",
	  getCurrentTime() );
  printOut(6, message);  
  
  /* closing log file */
  fclose(logfile);

  /* deleting lock file */
  sprintf(message,"rm -f .coinext.lock");
  system( message );
  
  /* ... and exit program (good style?) */
  exit(0);

  /* the answer never returned... */
  return 42;
}

/******************************
printOut:
   function to print messages on screen AND maybe into some file
*******************************/
int printOut( int mode, char* text )
{

  switch (mode) {
  case 0:
    printf("%s [TEST] %s\n",getCurrentTime(), text);
    break;
  case 1:
     printf("%s [WARNING] %s\n",getCurrentTime(), text);
    break;
  case 2:
     printf("%s [ERROR] %s\n",getCurrentTime(), text);
    break;
  case 3:
     printf("%s [INFO] %s\n",getCurrentTime(), text);
    break;
  case 4:
     printf("%s [COIN] %s\n",getCurrentTime(), text);
    break;
  case 5:
     printf("%s [FILE] %s\n",getCurrentTime(), text);
    break;
  case 6:
     printf("%s\n",text);
    break;
  case 7:
    printf("%s [DEBUG] %s\n",getCurrentTime(), text);
    break;
  }

  fflush(stdout);
  return 1;
}

/* In general, we can't use the builtin `basename' function if available,
   since it has different meanings in different environments.
   In some environments the builtin `basename' modifies its argument.  */
char* basename1( char* name )
{
  char* base = name;

  for (; *name; name++)
    if ( *name == '/' )
      base = name + 1;

  return (char *) base;
}

/******************************
ce_wait:
   waiting for some time: executing sleep command
*******************************/
int ce_wait(int minutes)
{
  /* printing message */
  sprintf(message, "Waiting for %d minutes.",minutes);
  printOut(3, message);

  /* good night... */
  sleep(60*minutes);

  return 1;
}


/*******************************
ce_restart:
   restarting program, beginning from scratch 
   -> deleting all 'cel' files
********************************/
void ce_restart(void)
{
  char command[256];
  sprintf( command, "rm -f *.cel");
  system( command );
}


/******************************
sigint_handler:
   handles the SIGINT event send to this program
   resulting in shutdown
*******************************/
void sigint_handler(int sig)
{
  /* the next two lines just to avoid warning .... */
  int dummy;
  dummy=sig;  

  /* shutdown */
  ce_shutdown();
}


/******************************
print_usage:
   prints the usage of the code
*******************************/
static void print_usage(char *program)
{
  fprintf(stderr,
	  "Usage:  %s [options] \n" \
	  "The following options are recognized.  Options are required unless \n"\
	  "surrounded in []\n", program);
  fprintf(stderr, 
	  "  [--help]                      display this message\n"\
	  "  [--version]                   print version information and exit\n"\
	  "  [--dirInspiral]               specifies the directory where to find the inspiral XML\n"\
	  "                                Either this option must be specified or --recalc\n"\
	  "  [--dirData]                   Specifies the directoet where to put the data\n"\
	  "  [--refresh]                   Specifies the refresh cycle in minutes. Default: 30 minutes\n");
  fprintf(stderr,
	  "  [--timeWindow]                Specifies the time window to analyze around a external trigger\n"\
	  "  [--snrCut]                    Specifies a general SNR cut\n"\
	  "  [--restart]                   restarts the program, delets any intermediate data\n"\
	  "  [--recalc]                    starts a seperate DAG for each found external trigger\n");
  fprintf(stderr,
          "  [--ifo]                       need to be specified if running on LHO. Choices: (LLO, LHO)\n"\
	  "  [--trigger]                   Besides downloading the latest trigger file, this option can be used\n"\
	  "                                to choose a specific trigger file\n\n");
}



/*------------------------------
  arg_parse_check:
  parses the arguments given to the variables 
  ------------------------------ */
int arg_parse_check( int argc, char *argv[])
{
  int flagRun=0;

  /* getopt arguments */
  struct option long_options[] =
    {
      /* these options set a flag */
      {"restart",                 no_argument,       &flagRestart,  1 },
      {"test",                    no_argument,       &flagTest,     1 },
      {"recalc",                  no_argument,       &flagRecalc,   1 },
      {"verbose",                 no_argument,       &flagVerbose,  1 },
      /* these options don't set a flag */
      {"dirInspiral",          required_argument, 0,                'a'},
      {"dirData",              required_argument, 0,                'A'},
      {"refresh",              required_argument, 0,                'b'},
      {"timeWindow",           required_argument, 0,                'B'},
      {"snrCut",               required_argument, 0,                'd'},
      {"ifo",                  required_argument, 0,                'i'},
      {"trigger",              required_argument, 0,                't'},
      {"calibration",          required_argument, 0,                'c'},
      {"run",                  required_argument, 0,                'r'},
      {"help",                 no_argument,       0,                'h'}, 
      {"version",              no_argument,       0,                'V'},
      {0, 0, 0, 0}
    };
  
  int c;
  
  
  /*
   *
   * parse command line arguments
   *
   */
  
  
  while ( 1 ) {
    /* getopt_long stores long option here */
    int option_index = 0;
    size_t optarg_len;
    
    c = getopt_long_only( argc, argv, 
			  "A:B:C:D:E:F:G:H:I:J:K:L:M:N:O:P:Q:R:S:T:U:V:W:X:Y:Z:"
			  "a:b:c:d:e:f:g:h:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:",
			  long_options, &option_index );
    
    /* detect the end of the options */
    if ( c == - 1 ) {
      break;
    }
    
    switch ( c ) {
    case 0:
      /* if this option set a flag, do nothing else now */
      if ( long_options[option_index].flag != 0 ) {
	break;
      } else {
	fprintf( stderr, "error parsing option %s with argument %s\n",
		 long_options[option_index].name, optarg );
	exit( 1 );
      }
      break;
      
    case 'a': 
      
      /* get name for inspiral data path */
      optarg_len = strlen( optarg ) + 1;
      memcpy( dirInspiral, optarg, optarg_len );
      ADD_PROCESS_PARAM( "string" , "%s", optarg);	      
      break;
      
    case 'A': { 
      
      /* get name for directory to put the data in */
      optarg_len = strlen( optarg ) + 1;
      memcpy( dirData, optarg, optarg_len );	      
      ADD_PROCESS_PARAM( "string" , "%s", optarg);
      flagDirData=1;
      break;
    }
      
    case 'b': {
      waitingTime = atol( optarg );
      if ( waitingTime < 0 ) {
	fprintf( stderr, "invalid argument to --%s:\n"
		 "waiting time can'y be negative\n",
		 long_options[option_index].name );
	exit( 1 );
      }  
      ADD_PROCESS_PARAM( "int" , "%ld", waitingTime);
    }
      break;
      
    case 'B': {
      timeWindow = atoi( optarg );
      if ( timeWindow < 0 ) {
	fprintf( stderr, "invalid argument to --%s:\n"
		 "window time must be positive number\n",
		 long_options[option_index].name );
	exit( 1 );
      } 
      ADD_PROCESS_PARAM( "float" , "%s", optarg);	      
    }
      break;
      
      
    case 'd': {
      snrCut = atof( optarg );
      if ( snrCut < 0 ) {
	fprintf( stderr, "invalid argument to --%s:\n"
		 "snrCut must be a positive number\n",
		 long_options[option_index].name );
	exit( 1 );
      } 
      ADD_PROCESS_PARAM( "float" , "%s", snrCut);	      
    }
      break;
      
    case 'i': {
      /* read IFO site */
      optarg_len = strlen( optarg ) + 1;
      memcpy( ifo, optarg, optarg_len );	
      ADD_PROCESS_PARAM( "string" , "%s", optarg);	      	      
      break;
    }
      
    case 't': {
      /* read triger file */
      optarg_len = strlen( optarg ) + 1;
      memcpy( nameInputList, optarg, optarg_len );	
      ADD_PROCESS_PARAM( "string" , "%s", optarg);
      flagTriggerFile=1;	      	      
      break;
    }

    case 'c': {
      /* read calibration vesrion to use */
      optarg_len = strlen( optarg ) + 1;
      memcpy( calibrationVersion, optarg, optarg_len );	
      ADD_PROCESS_PARAM( "string" , "%s", optarg);
      break;
    }

    case 'r': {
      /* read run file */
      optarg_len = strlen( optarg ) + 1;
      memcpy( run, optarg, optarg_len );	
      ADD_PROCESS_PARAM( "string" , "%s", optarg);
      flagRun=1;
      break;
    }
      
    case 'h': {
      /* help message */
      print_usage(argv[0]);
      exit( 1 );
      break;
    }
      
    case 'V': {
      /* print version information and exit */
      fprintf( stdout, "COINcidences with EXTernal triggers\n" 
	       "Alexander Dietz\n"
	       "CVS Version: " CVS_ID_STRING "\n"
	       "CVS Tag: " CVS_NAME_STRING "\n" );
      fprintf( stdout, lalappsGitID );
      exit( 0 );
      break;
    }
    
    }
  }
  
  if ( optind < argc ) {
    fprintf( stderr, "extraneous command line arguments:\n" );
    while ( optind < argc ) {
      fprintf ( stderr, "%s\n", argv[optind++] );
    }
    exit( 1 );
  }

  /* check for needed input arguments */
  if ( !flagRun ) {
    fprintf( stderr, "Must specify '--run [RUN]', e.g. A4 or S5. Aborted.\n");
    exit( 0 );
  }

  /* setting flags into table */
  
  /* restart-flag */
  if ( flagRestart ) {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
		 "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
		 "--restart" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }
  
  /* test-flag */
  if ( flagTest ) {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
		 "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--test" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }
  
  /* recalc-flag */
  if ( flagRecalc ) {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
		 "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
		 "--recalc" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }

  /* verbose-flag */
  if ( flagRecalc ) {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX,
                 "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
                 "--verbose" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }

  
  return 0;
}






// ****************************************************************************
// File:        coinExt.c
// Created:     2004/10/22
// Rev:         2005/07/25 16:13
// Version:     1.0.0
// Author:      Alexander Dietz
// Description: C(!) program for automatic search for coincidences
//              between Inspiral trigger and External triggers (e.g. GRB)
// ****************************************************************************

// ******************************
// include
// ******************************
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

  LIGOTimeGPS startTime;
  LIGOTimeGPS endTime;
} InspiralList;

typedef struct S_ExternalList {
  long             gps;
  long             gpsNano;
  int              status;
  ExtTriggerTable* table;
} ExternalList;

char* basename1( char* name);
void sigint_handler(int sig); 
void writeCGN( void );
void checkInspiralTrigger( void );
char* getIFOname(int index);
void restart( void );
void startAnalysisJob(ExternalList* eList, int ifo);


/******************************
   variables
*******************************/
  
LALLeapSecAccuracy accuracy=LALLEAPSEC_LOOSE;

int flagTest=0;    /* flag indicating a test of the program */
int flagRestart;   /* flag indicating to restart the complete DAQ */
int flagRecalc=0;  /* flag indicating to redo the inspiral analysis */
int flagTriggerFile=0; /* flag indicating to use special trigger file */

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
char* nameLogfile          = "coinExt.log";
char* nameSaveInspiralList = "saveInspiral.cel";
char* nameInputList        = "triggers.xml";  
char* nameInputList2       = "triggersIsabel.list";  
char* nameSaveExternalList = "externalTriggers.cel";
//char* dirData              = "this text will never ne used, just for memory allocation";
char* dirData=NULL;
char* dirInspiral          = "/data2/dietz/S4/CoinExt/cluster";
char* destCopyURL          = "gsiftp://adietz.phys.lsu.edu/home/alex/Work/E_ExtTrigger/Online/Data/";
char* serverName           = "ldas.ligo-la.caltech.edu";

/*  web page with avaliable triggers */
char* webPage = " http://www.ligo.caltech.edu/~wk_sz/triggers/triggers.xml";
char* webPage2 = "http://www.uoregon.edu/~ileonor/ligo/a4/grb/online/currentgrbs_html_A4.txt";

/* where does the prog run? LLO or LHO */
char* ifo = "LLO";
char* ifoName[]={"L1","H1", "H2"};

/* times and name of periods (extentable!) */
int numberPeriods=3;
long* startPeriod;
long* endPeriod;
char* namePeriod[]={"S4","A4","S5"};

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
double snrCut=7.0;
long timeWindow=1800;

/* standard waiting time in minutes */
int waitingTime=180; // 3 hours

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
ProcessParamsTable* this_proc_param=NULL;  // usage??? 

SummValueTable* listSumm=NULL;  // tables to use for summary (times of each used xml file)
MetadataTable   listSummTable;

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
  char message[256];
  char* env;

  /* create the process and process params tables */  
  LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;
  proctable.processTable = (ProcessTable *) calloc( 1, sizeof(ProcessTable) );
  LALGPSTimeNow ( &status, &(proctable.processTable->start_time), &accuracy );
  populate_process_table( &status, proctable.processTable, 
			  PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE );
  this_proc_param = processParamsTable.processParamsTable = 
    (ProcessParamsTable *) calloc( 1, sizeof(ProcessParamsTable) ); // ????

  /* first of all set default dirData */
  env=getenv("PWD");
  sprintf(message, "%s/../Data/", env);
  dirData=(char*)calloc( strlen(message+1), sizeof(char));
  strcpy(dirData, message);

  /* call the argument parse and check function */
  arg_parse_check( argc, argv);

  /*  startup */
  startup();
  if (signal(SIGINT, sigint_handler) == SIG_ERR) {
    perror("signal");
    exit(1);
  }	

  /* check if to restart all */
  if (flagRestart) {
    printOut("--- RESTARTING program. Deleting all data files...");
    restart();
  }

  if (!flagRestart) {
    /* reading information from CGN trigger that are already available */
    num=readCGN();
    sprintf(message, "--- Reading %d CGN triggers from file '%s'",num, nameSaveExternalList);
    printOut(message);
  }

  /* setting some loop parameters */
  if (strcmp (ifo,"LLO") == 0) {
    startIFO=0;
    endIFO=1;
  } else {
    startIFO=1;
    endIFO=2;
  }

  /* printing some user parameters */  
  printf("--- waitingTime: %d minutes \n",waitingTime);
  printf("--- dirInspiral: %s\n",dirInspiral);
  printf("--- dirData:      %s\n",dirData);
  printf("--- running at site %s\n",ifo);
  if (flagTriggerFile) {
    printf("--- using trigger file %s\n",nameInputList2);
  }
  if (flagRecalc) {
    printf("--- starting analysis dag for every GRB\n");
  }
  if (flagTest) {
    printOut("--- TEST mode activated!");
    waitingTime=1;
  }

  printf("\n");


  /* start 'endless' loop */
  do {
    
    /* check latest GRB... */
    //downloadCGN();
    /* and check if there are some new CGN */
    //int number=analyzeNewCGN();  
    //printf("--- Number of unanalyzed CGN triggers: %d \n",number);

    if (checkCGNtrigger2()) {

      /* store all CGN triggers to file */
      writeCGN();
    }
     
    /* check for new XML files */
    checkInspiralTrigger();
    
    /* check coincidents with external trigger */
    coinInspiralCGN();


    /* waiting for some time... */
    wait(waitingTime);

  } while(1);

}


/******************************
checkCGNtrigger
  downloading the latest CGN trigger from the web page
  and compare it with the stored list to check for new CGN triggers
*******************************/
int checkCGNtrigger( void )
{
  /* specify used variables here ... */
  int i,ifo;
  int counter=0;
  int numEvents=0;
  int index=0;
  int triggersUnanalyzed=0;
  char command[256];
  char message[256];
  char dirJob[256];
  LIGOTimeGPS* gps;
  LALDate* date;

  date=(LALDate*)LALCalloc(1, sizeof(LALDate));
  gps=(LIGOTimeGPS*)LALCalloc(1, sizeof(LIGOTimeGPS));


  if (flagTest) {

    /* copy test file if it is a test */
    system( "cp Test/triggers.xml ." );

  } else {

    /* downloading file */
    sprintf(command, "wget %s -O triggersUncheck.xml > download.log 2>&1",
	    webPage);
    system(command);

    /* recheck file (taking out corrupt lines) [added 05/03/01] */
    sprintf(command,"gawk -f recheckTriggers.awk triggersUncheck.xml > "
	    "triggers.xml");
    system(command);
  }

  /* delete something, (only neccessary for C code) */
  while ( headExt )  {
    thisExt = headExt;
    headExt = headExt->next;
    LALFree( thisExt );
  }

    
  /* read data from xml file */
  numEvents=LALExtTriggerTableFromLIGOLw( &headExt, nameInputList, 0, -1); 
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
	//extList[numberStoredCGN].table=(ExtTriggerTable*)LALCalloc( 1, sizeof(ExtTriggerTable) );
	//memcpy( extList[numberStoredCGN].table, thisExt, sizeof(ExtTriggerTable));

	extList[numberStoredCGN].gps=thisExt->start_time;
	extList[numberStoredCGN].gpsNano=thisExt->start_time_ns;
	extList[numberStoredCGN].status=0;
	/* extList[numberStoredCGN].table=thisExt;*/  /* just setting the pointer */	
	extList[numberStoredCGN].table=(ExtTriggerTable*)LALCalloc( 1, sizeof(ExtTriggerTable) );
	memcpy( extList[numberStoredCGN].table, thisExt, sizeof(ExtTriggerTable));

	//printf("CGN gps: %d\n",thisExt->start_time);
	/*
	sprintf(message,"creating trigger table for CGN %d",numberStoredCGN);
	printOut(message);
	if (extList[numberStoredCGN].table) {
	  sprintf(message,"table allocated...");
	} else {
	  sprintf(message,"table NOT allocated...");
	}
	  printOut(message);
	*/

	/* increase number of stored triggers from the C... GRB Network */
	//	numberStoredCGN++;
	/* number increased at the end... */
	/* getting user friendly date and time */
	gps->gpsSeconds=thisExt->start_time;
	gps->gpsNanoSeconds=thisExt->start_time_ns;
	LALLeapSecAccuracy accuracy= LALLEAPSEC_STRICT;
	LALGPStoUTC(&status, date, gps, &accuracy); 

	/* generate output saying that there is a new CGN trigger */
	sprintf(message,"--- New CGN trigger occured at %d/%d/%d %d:%d:%d (%d)",
		date->unixDate.tm_year+1900, date->unixDate.tm_mon+1,
		date->unixDate.tm_mday,date->unixDate.tm_hour,
		date->unixDate.tm_min, date->unixDate.tm_sec, 
		gps->gpsSeconds);
	printOut(message);

	// start special analysis job if data should be recalculated
	if (flagRecalc) {
	  for ( ifo=startIFO;ifo<endIFO;ifo++) {
	    startAnalysisJob( &extList[numberStoredCGN], ifo );
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

  /* get the number of unanalyzed triggers */
  triggersUnanalyzed=0;
  for (i=0;i<numberStoredCGN;i++) {
    if (extList[i].status==0) triggersUnanalyzed++;
  }

  /* return number of new E.T. */
  return triggersUnanalyzed;
}


/******************************
checkCGNtrigger2
  downloading the latest CGN trigger from Isabels list
*******************************/
int checkCGNtrigger2( void )
{
  /* specify used variables here ... */
  int i,ifo;
  int counter=0;
  int numEvents=0;
  int index=0;
  int triggersUnanalyzed=0;

  char command[256];
  char message[256];
  char input[256];
  char inputGPS[256];
  char input1[256];
  char input2[256];
  LIGOTimeGPS* gps;
  LALDate* date;
  FILE* fileTrigger;

  /* initialize date */
  date=(LALDate*)LALCalloc(1, sizeof(LALDate));
  gps=(LIGOTimeGPS*)LALCalloc(1, sizeof(LIGOTimeGPS));

  if (flagTest) {

    /* copy test file if it is a test */
    system( "cp Test/triggersIsabel.list ." );

  } else {

    /* downloading file */
    sprintf(command, "wget %s -O triggersIsabel.list > download2.log 2>&1",
	    webPage2);
    system(command);
  }

  /* delete something, (only neccessary for C code) */
  while ( headExt )  {
    thisExt = headExt;
    headExt = headExt->next;
    LALFree( thisExt );
  }
  
  /* open the trigger file */
  fileTrigger=fopen(nameInputList2, "r");
  
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

    // reading GRB number 
    fscanf(fileTrigger, "%s", input);
    strncpy(thisExt->event_number_grb, input,8);

    // reading UT time -- not stored
    fscanf(fileTrigger, "%s", input);    

    // reading  GPS time 
    fscanf(fileTrigger, "%s", inputGPS);
    strncpy(input1, inputGPS, 9);
    strncpy(input2, inputGPS+9, strlen(inputGPS)-10);
    thisExt->start_time=atol(input1);    
    thisExt->start_time_ns=(long)(atof(input2)*1.0e+9);

    // reading RA
    fscanf(fileTrigger, "%s", input);
    thisExt->event_ra=atof(input);

    // reading DEC
    fscanf(fileTrigger, "%s", input);
    thisExt->event_dec=atof(input);

    // reading fave LHO
    fscanf(fileTrigger, "%s", input);
    thisExt->ligo_fave_lho=atof(input);

    // reading  fave LLO
    fscanf(fileTrigger, "%s", input);
    thisExt->ligo_fave_llo=atof(input);   

    // reading   
    fscanf(fileTrigger, "%s", input);
    thisExt->ligo_delay=atof(input);

    /* reading GCN number */
    fscanf(fileTrigger, "%5s", input);

    /* check end of file */
    numEvents++;
  } while (!feof(fileTrigger));

  /* TODO: delete */
  /* test output */
  if (1==2) {
    thisExt=headExt;
    while ( thisExt )  {
      printf("Reading data from structure:\n");
      printf("GCN: %d | GRB: %s \n",thisExt->event_number_gcn, thisExt->event_number_grb);
      printf("GPS: %d %d  | ra: %f | dec: %f | lho: %f | llo: %f | delay: %f\n\n",thisExt->start_time, thisExt->start_time_ns, thisExt->event_ra,thisExt->event_dec, thisExt->ligo_fave_lho, thisExt->ligo_fave_llo, thisExt->ligo_delay);    
      thisExt = thisExt->next;
    }
  }

  if (numEvents<0) {
    sprintf(message, "Error: Unable to read data from file %s",nameInputList2);
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
	extList[numberStoredCGN].status=0;
	extList[numberStoredCGN].table=(ExtTriggerTable*)LALCalloc( 1, sizeof(ExtTriggerTable) );
	memcpy( extList[numberStoredCGN].table, thisExt, sizeof(ExtTriggerTable));

	//printf("CGN gps: %d\n",thisExt->start_time);
	/*
	sprintf(message,"creating trigger table for CGN %d",numberStoredCGN);
	printOut(message);
	if (extList[numberStoredCGN].table) {
	  sprintf(message,"table allocated...");
	} else {
	  sprintf(message,"table NOT allocated...");
	}
	  printOut(message);
	*/

	/* increase number of stored triggers from the C... GRB Network */
	//	numberStoredCGN++;
	/* number increased at the end... */
	/* getting user friendly date and time */
	gps->gpsSeconds=thisExt->start_time;
	gps->gpsNanoSeconds=thisExt->start_time_ns;
	LALLeapSecAccuracy accuracy= LALLEAPSEC_STRICT;
	LALGPStoUTC(&status, date, gps, &accuracy); 

	/* generate output saying that there is a new CGN trigger */
	sprintf(message,"--- New CGN trigger occured at %d/%d/%d %d:%d:%d (%d)",
		date->unixDate.tm_year+1900, date->unixDate.tm_mon+1,
		date->unixDate.tm_mday,date->unixDate.tm_hour,
		date->unixDate.tm_min, date->unixDate.tm_sec, 
		gps->gpsSeconds);
	printOut(message);

	// start special analysis job if data should be recalculated
	if (flagRecalc) {
	  for ( ifo=startIFO;ifo<endIFO;ifo++) {
	    startAnalysisJob( &extList[numberStoredCGN], ifo );
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

  /* get the number of unanalyzed triggers */
  triggersUnanalyzed=0;
  for (i=0;i<numberStoredCGN;i++) {
    if (extList[i].status==0) triggersUnanalyzed++;
  }


  /* return number of new E.T. */
  return triggersUnanalyzed;
}

/******************************
readCGN:
  reading data for stored CGN (only once at start of program)
  executed in the 'startup' method 'externalTriggers.cel'
*******************************/
int readCGN( void )
{ 
  FILE* file;
  char message[256];
  int c=0;
  long gps;
  long gpsNano;
  int status;
  int i=0;

  /* open file for reading data */
  file=fopen( nameSaveExternalList , "r" );
  
  /* check file */
  if (!file) {
    sprintf(message,">>> File %s does not exist",nameSaveExternalList);
    printOut(message);
    
    numberStoredCGN=0;
    return numberStoredCGN;
  }

  /* loop over all entries in the file */
  do {

    /* read next line and store data */
    fscanf(file, "%ld,%ld,%d", &gps, &gpsNano, &status); 
    extList[c].gps      = gps;
    extList[c].gpsNano  = gpsNano;
    extList[c].status   = status;
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
  char message[256];
  int i;

  /* open output file for overwriting all data! */
  file=fopen(nameSaveExternalList,"w");
  for ( i=0; i<numberStoredCGN; i++ ) {    
    fprintf( file, "%ld,%ld,%d\n",
	     extList[i].gps, extList[i].gpsNano, extList[i].status);
  }
  fclose(file);

}

/*******************************
checkInspiralTrigger:
  checking all available xml-files and store corresponding times in a special file.
  If new data are present, a 1 is returned, 0 else.
  checking database...
  REMARK: Again not very nifty. Always reads the previous data
********************************/
void checkInspiralTrigger( void )
{
  int i,j,c=0;
  char* ifo;
  char dummy[256];
  char command[256];
  char text[256];
  char* filename;
  char savename[256];
  FILE* file;
  FILE* out;

  /* loop over IFO's to take into account */
  for ( i=startIFO;i<endIFO;i++) {

    /* get name of IFO */
    ifo=ifoName[i];
        
    /* create new file containing all present inspiral trigger lists */
    if (flagRecalc) {
      /* use directory for recalculated trigger */
      sprintf( command, "ls -1 ../Trigger/%s/*.xml > namesInspiral.cel",ifo);
    } else{
      /* use directory for online created triggers */
      sprintf( command,"ls -1 %s%s/*.xml > namesInspiral.cel",dirInspiral,ifo );
    }
    //printOut(command);
    system( command );
    
    /* read list and extract time data */
    c=0;
    file=fopen( "namesInspiral.cel","r" );
    do {      

      /* read next filename and remove directory path from it */
      fscanf(file, "%s", &text); 
      filename = basename1( (char*)text );
      //printf("before: %s | after: %s \n",text, filename);

      /* extracting relevant parts of the filename
         assuming the length are always the same! */
      strcpy(inspiralList[i][c].filename, filename);
      strncpy(dummy, filename+8, 9);dummy[9]='\0';
      inspiralList[i][c].timeStart=atoi(dummy);
      strncpy(dummy, filename+18, 4);dummy[4]='\0';
      inspiralList[i][c].timeLength=atoi(dummy);
      inspiralList[i][c].coin=0;

      //printf("filename: %s  | start: %ld  length: %ld\n",filename, inspiralList[i][c].timeStart,  inspiralList[i][c].timeLength);

      /* increasing counter */
      c++;
      
    } while ( !feof(file) );
    fclose(file);
    
    /* storing total number of inspiral xml files */
    numberInspiralList[i]=c-1;
    
    /* open file and store all entries of the inspiral list */
    sprintf(savename, "%s%s",ifo, nameSaveInspiralList);    
    out=fopen( savename, "w" );
    for (j=0;j<numberInspiralList[i]; j++) {
      fprintf(  out, "%s %d %d\n",inspiralList[i][j].filename,
		inspiralList[i][j].timeStart,inspiralList[i][j].timeLength);
    }
    fclose(out);

  } 

  /* removing unused file */
  //system( "rm -f namesInspiral.cel" );
}

/********************************
coinInspiralCGN:
  search for coincidences of inspiral and CGN triggers
  and store data to new xml file
*********************************/
int coinInspiralCGN( void )
{
  /* variables */
  int cgn, lifo, insp;  
  int numberFound=0;
  double timeCGN;
  char message[256];
  SnglInspiralTable* element;
  ExtTriggerTable* table;
  long timeStart, timeEnd;

  /* loop over ALL CGN triggers  */
  for ( cgn=0; cgn<numberStoredCGN; cgn++ ) {
    
    /* check status of the current CGN triggers */
    if (extList[cgn].status<4) {
      
      /* get pointer to corresponding E.T. structure */
      table=extList[cgn].table;
      if ( !table ) {
	printf("WARNING in coinInspiralCGN: A ExtTriggerTable does not exist for cgn %d. " 
	       "Check program!\n",cgn);
	continue;
      }

      /* get time of CGN trigger */
      timeCGN=(double)table->start_time + (double)table->start_time_ns/1.0e+9;
      
      /* loop over IFO's */
      for (lifo=startIFO;lifo<endIFO;lifo++) { 

	/* loop over accessible inspiral-trigger XML files */
	numberFound=0;
	for (insp=0; insp<numberInspiralList[lifo]; insp++ ) {

	  inspiralList[lifo][insp].coin=0;

	  /* get starting end ending time of the time period
	     the inspiral triggers are from this file */ 
	  timeStart=inspiralList[lifo][insp].timeStart;
	  timeEnd=timeStart+inspiralList[lifo][insp].timeLength;


	  /* check if time overlaps */ 
	  if ( (timeStart>=timeCGN-timeWindow && timeStart<=timeCGN+timeWindow) || 
	       (timeEnd>=timeCGN-timeWindow && timeEnd<=timeCGN+timeWindow) ) {
	    
	    //printf("timeStart: %d  timeEnd: %d   timeCGN: %f\n",timeStart, timeEnd, timeCGN);
	    /* mark corresponding Inspiral file */
	    numberFound++;
	    inspiralList[lifo][insp].coin=1;
	  }
	}
	
	// set status: inspiral data found
	//	extList[cgn].status=4;

	if (numberFound) {

	  /* output */
	  sprintf(message, "### ExtTrigger found in coincidence at time %f", timeCGN);
	  printOut(message);

	  /* creating final xml file and copy it */
	  createInspiralTable( lifo );
	  writeXML( lifo, cgn );

	  /* set new status */
	  extList[cgn].status=5;

	  /* delete old inspiral triggers */
	  while ( list )  {
	    element = list;
	    list = list->next;
	    LALFree( element );
	  }
  
	}
      }

    }
  }

}

/********************************
createInspiralTable:
  create table with all needed inspiral triggers in it
*********************************/
int createInspiralTable( int lifo )
{
  /* variables */
  int insp,numTrig, numSum,flag;
  int numberDiscard;
  char filename[256];  
  SnglInspiralTable* head;
  SnglInspiralTable* tail;
  SnglInspiralTable* mainPointer;
  SnglInspiralTable* prevElement;
  SnglInspiralTable* element;
  SummValueTable* summary;
  SummValueTable* tailSumm;
  tail=NULL;

  /* loop over all inspiral files */
  flag=0;
  for (insp=0; insp<numberInspiralList[lifo]; insp++ ) {

    /* check if file contain useful data */
    if (inspiralList[lifo][insp].coin) {

      /* read data from the XML file using LAL functions 
	 into new allocated memory, pointed by 'head' */
      head=NULL; 

      if (flagRecalc) {
	/* use directory for recalculated trigger */
	sprintf( filename,  "../Trigger/%s/%s", getIFOname(lifo), inspiralList[lifo][insp].filename);
      } else{
	/* use directory for online created triggers */
	sprintf( filename, " %s%s/$s",
		 dirInspiral, getIFOname(lifo), inspiralList[lifo][insp].filename);
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

	/* set inspiral data into table ONLY IF there are some inspiral data */	
	//if (numTrig-numberDiscard>0)
	//list=head;

	/* set summ values in table */
	listSumm=summary;
	tailSumm=listSumm;
      } else {

	/* set inpiral data into table */
	//if (numTrig-numberDiscard>0)
	//tail->next=head;

	/* set summ values in table */
	tailSumm->next=summary;
	tailSumm=tailSumm->next;
      }


      /* store last element in the total list 
	 (to set the 'next' parameter to next part of list */
      //if (prevElement) 
      //tail=prevElement;


      /* searching last entry in (inspiral) table */
      //if (numTrig-numberDiscard>0) {
      //tail=head;
      //while (tail->next) { tail=tail->next; }
      //}

      /* setting values in the summ-value table */
      LALSnprintf( tailSumm->name, LIGOMETA_SUMMVALUE_NAME_MAX, "length" ); 
      tailSumm->value=summary->end_time.gpsSeconds - summary->start_time.gpsSeconds;
      LALSnprintf( tailSumm->comment, LIGOMETA_SUMMVALUE_COMM_MAX, "%s",inspiralList[lifo][insp].filename ); 

      printf("+++ Adding %d-%d triggers from file %s \n", numTrig, numberDiscard,inspiralList[lifo][insp].filename);  
   
    }
  }
}


/********************************
writeXML:
  writing XML file with all needed informations
*********************************/
int writeXML( int ifo, int cgnIndex )
{
  /* variables */
  char message[256];
  char command[256];
  char filename[256];
  char filenameOut[256];
  MetadataTable outputTable; 
  ExtTriggerTable* exTrig;
  ExtTriggerTable* table;
  SnglInspiralTable*    thisInsp = NULL;
 
  /* create filename for output file */
  sprintf(filename, "externalTrig%d-%s.xml", extList[cgnIndex].gps, getIFOname(ifo));
  sprintf( filenameOut, "%s%s", dirData, filename);
  sprintf( message, "+++ Creating output file %s", filenameOut );
  printOut(message);

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
  printf("%s\n",command);
  system(command);

}

/****************************** 
startAnalysisJob
   creates and starts a special analysis job
********************************/
void startAnalysisJob(ExternalList* eList, int ifo)
{
  int i,indexPeriod;
  char dirJob[256];
  char command[256];
  char message[256];
  int flagPrint=1;
  FILE* fileSegment;
  double a,b,c,d,sumTime=0,maxSegment=0;
  char dummy[256];

  if (ifo>0) {
    serverName="ldas.ligo-wa.caltech.edu";
  }

  // discard any times before S4
  if (eList->gps<793130413) 
    return;

  /* search for corresponding time period */
  indexPeriod=-1;
  for (i=0;i<numberPeriods;i++) {
    if (eList->gps>startPeriod[i] && eList->gps<endPeriod[i]) {
      indexPeriod=i;
    }
  }

  /* check if a correct period was found */
  if (indexPeriod==-1) {
    sprintf(message,"ERROR: GPS time %9d lies outside of senseful range!",eList->gps);
    printOut(message);
    return;
  }

  /* find possible segments that can be used */
  sprintf(command, "LSCsegFind --server=%s --interferometer %s --type Science --output-format segwizard --gps-start-time %9d --gps-end-time %9d > segment.txt 2> .notfound.cel ", 
	  serverName, ifoName[ifo], eList->gps-timeWindow, eList->gps+timeWindow );
  if (flagPrint) printOut(command);
  system (command);

  /* check file segment.txt */
  fileSegment=fopen("segment.txt", "r");
  do {    
    fscanf(fileSegment, "%s",dummy);
    fscanf(fileSegment, "%s",dummy);
    fscanf(fileSegment, "%s",dummy);
    fscanf(fileSegment, "%s",dummy);

    //    fscanf(fileSegment, "%f %f %f %f",a,b,c,d);
    d=atof(dummy);
    
    // get sum of times and maximum segment (has to be > 2048 seconds)
    if (!feof(fileSegment)) {
      if (d>maxSegment) maxSegment=d;

      sumTime=sumTime+d;
    }
  } while (!feof(fileSegment));

  printf("sum: %f\n",sumTime);

  if (d>=2048) {

    /* create name of new directory for job execution */
    sprintf(dirJob,"../OnlineAnalysis/Jobs/%s%s-%9d", ifoName[ifo], namePeriod[indexPeriod], eList->gps);
    
    /* create directory */
    sprintf(command,"mkdir %s",dirJob);
    if (flagPrint) printOut(command);
    system(command);
    
    /* copy the segment file to that directory */
    sprintf(command, "cp segment.txt %s",dirJob);
    if (flagPrint) printOut(command);
    system(command);

    /* symbolic links to all needed basic files */
    sprintf(command,"cd %s; ln -s ../../ProgInspiral/Copy/* .",dirJob);
    if (flagPrint) printOut(command);
    system(command);
    
    /* create ini file for job execution */
    sprintf(command,"cp ../OnlineAnalysis/ProgInspiral/online%s%s.ini %s", 
	    ifoName[ifo], namePeriod[indexPeriod], dirJob);
    if (flagPrint) printOut(command);
    system(command);
    
    
    /* copy the segment to the directory */
    sprintf(command, "cp segment.txt %s", dirJob);
    if (flagPrint) printOut(command);
    system (command);

    /* find the segments that can be used */
    //sprintf(command, "cd %s; LSCsegFind --server=%s --interferometer %s --type Science --output-format segwizard --gps-start-time %9d --gps-end-time %9d > segment.txt", dirJob, serverName, ifoName[ifo], eList->gps-timeWindow, eList->gps+timeWindow );
    //if (flagPrint) printOut(command);
    //system (command);
    
    sprintf(command,"cd %s; ./lalapps_inspiral_online_pipe --config-file online%s%s.ini --log-path /usr1/dietz",
	    dirJob,ifoName[ifo], namePeriod[indexPeriod] );
    if (flagPrint) printOut(command);
    system(command);
    
    sprintf(command,"+++ Starting DAG for analysis in directory %s",dirJob);
    printOut(command);
    sprintf(command, "cd %s; condor_submit_dag online%s%s.dag",
	    dirJob,ifoName[ifo], namePeriod[indexPeriod]);
    if (flagPrint) printOut(command);
    system (command);
    
    // set status 
    eList->status=2;

  } else {

    /* set status */
    eList->status=4;

    /* output */
    if (d==0) {
      printf("--- No segments found for times %d to %d on %s\n",
	     eList->gps-timeWindow, eList->gps+timeWindow, ifoName[ifo] );      
    } else {
      printf("--- Maximum segment for times %d to %d on %s is only %f seconds long\n",
	     eList->gps-timeWindow, eList->gps+timeWindow, ifoName[ifo], maxSegment );
      
    }
  }


  /* TODO idea: remove executables in the directories with the post-script */
}

/****************************** 
getIFOname:
return the IFO name corresponding to 'index'
********************************/
char* getIFOname(int index)
{
  switch (index) {
  case 0:
    return "L1";
  case 1:
    return "H1";
  case 2:
    return "H2";
  }
  return "";
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
startup:
   startup function: looks for existence of the file .coinExt.lock
   if that file exist -> abort. Program running.
   if not exist -> create this file with command stopping this program 
   (sending kill signal 2 to this program, which is fetched by the function sigint_handler)
*******************************/
int startup( void )
{ 
  FILE* out;
  int num;
  char message[256];
  time_t rawtime;
  struct tm * timeinfo;

  /* checking/creating lock-file */
  out=fopen( ".coinExt.lock" ,"r");
  if (out) {
    if (!flagTest) {
      printf("ERROR: Program already running!\n");
      exit(0);
    }
  } else {
    out=fopen( ".coinExt.lock", "w" );
    fprintf( out, "kill -2 %d\n", getpid());
    fprintf( out, "rm -f .coinExt.lock\n" );
    fclose(out);
  }

  /* open log file for appending */
  logfile=fopen( nameLogfile , "a" );

  /* get curreent time */
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );

  /* create startup message */
  sprintf(message, "-------------------------------------------------\n--- Starting coinExt at %s",
	   asctime (timeinfo));
  printOut(message);

  /* reading already stored CGN triggers */
  //num=readCGN();
  //sprintf(message, "--- Reading %d events from file %s",num, nameSaveExternalList);
  //printOut(message);

  startPeriod=(long*)LALCalloc( numberPeriods, sizeof(long) );
  endPeriod=(long*)LALCalloc( numberPeriods, sizeof(long) );
  startPeriod[0]=793130413; // S4
  endPeriod[0]=795679213;
  startPeriod[1]=795813760; // A4
  endPeriod[1]=811711310;
  startPeriod[2]=811711311; //S5 (preliminary)
  endPeriod[2]=999999999;
  //long startPeriod=
  //long endPeriod[3]={800000000, 900000000, 1000000000};
  //{700000001, 800000001, 900000001};  // TODO

  return 1;
}

/*******************************
shutdown:
   shutdown function
*******************************/
int shutdown( void )
{
  char message[256];
  time_t rawtime;
  struct tm * timeinfo;

  /* store CGN triggers for next time (status is important!!!) */
  writeCGN();

  /* get curreent time */
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );

  /* creating shutdown message */
  sprintf(message, "--- Stopping coinExt at %s "
	  "\n-------------------------------------------------",
	  asctime (timeinfo));
  printOut(message);  
  
  /* closing log file */
  fclose(logfile);

  /* deleting lock file */
  sprintf(message,"rm -f .coinExt.lock");
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
int printOut(char* text)
{
  printf("%s\n",text);
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
wait:
   waiting for some time: executing sleep command
*******************************/
int wait(int minutes)
{
  /* printing message */
  char message[256];
  sprintf(message, "--- Waiting for %d minutes.",minutes);
  printOut(message);

  /* good night... */
  sleep(60*minutes);

  return 1;
}


/*******************************
restart:
   restarting program, beginning from scratch 
   -> deleting all 'cel' files
********************************/
void restart(void)
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
  shutdown();
}


/*------------------------------
  arg_parse_check:
  parses the arguments given to the variables 
  ------------------------------ */
int arg_parse_check( int argc, char *argv[])
{
  // char* dummy;  

  /* getopt arguments */
  struct option long_options[] =
    {
      /* these options set a flag */
      {"restart",                 no_argument,       &flagRestart,  1 },
      {"test",                    no_argument,       &flagTest,     1 },
      {"recalc",                  no_argument,       &flagRecalc,   1 },
      /* these options don't set a flag */
      {"dirInspiral",          required_argument, 0,                'a'},
      {"dirData",              required_argument, 0,                'A'},
      {"refresh",              required_argument, 0,                'b'},
      {"timeWindow",           required_argument, 0,                'B'},
      {"snrCut",               required_argument, 0,                'd'},
      {"ifo",                  required_argument, 0,                'C'},
      {"trigger",              required_argument, 0,                'c'},
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
      dirInspiral = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
      memcpy( dirInspiral, optarg, optarg_len );
      ADD_PROCESS_PARAM( "string" , "%s", optarg);	      
      break;
      
    case 'A': { 
      
      /* get name for directory to put the data in */
      optarg_len = strlen( optarg ) + 1;
      //if (dirData) free(dirData);
      dirData = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
      memcpy( dirData, optarg, optarg_len );	      
      ADD_PROCESS_PARAM( "string" , "%s", optarg);
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
	
    case 'C': {
      /* get name for data to put inh */
      optarg_len = strlen( optarg ) + 1;
      ifo = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
      memcpy( ifo, optarg, optarg_len );	
      ADD_PROCESS_PARAM( "string" , "%s", optarg);	      	      
      break;
    }

    case 'c': {
      /* get name for data to put inh */
      optarg_len = strlen( optarg ) + 1;
      nameInputList2 = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
      memcpy( nameInputList2, optarg, optarg_len );	
      ADD_PROCESS_PARAM( "string" , "%s", optarg);
      flagTriggerFile=1;	      	      
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
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
		 "--test" );
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
}





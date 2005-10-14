
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
  int              status[3];
  ExtTriggerTable* table;
} ExternalList;



int checkCGNtrigger( void );
int checkCGNtrigger2( void );
int readCGN( void );
void writeCGN( void );	
void checkInspiralTrigger( void );
int coinInspiralCGN( void );
int createInspiralTable( int lifo );
int writeXML( int ifo, int cgnIndex );
void startAnalysisJob(ExternalList* eList, int ifo);
char* getIFOname(int index);
int checkGPS( long gps, long gpsNano );
int getPeriod( long gps );
char* getCurrentTime( void );
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
static char nameLogfile[]          = "coinext.log";
static char nameSaveInspiralList[] = "saveInspiral.cel";
static char nameInputList[]        = "triggers.xml";  
static char nameInputList2[]       = "triggersIsabel.list";  
static char nameSaveExternalList[] = "externalTriggers.cel";
static char dirData[] ="";
static char dirInspiral[]          = "/data2/dietz/S4/CoinExt/cluster";
static char destCopyURL[]          = "gsiftp://dietz.phys.lsu.edu/home/alex/Work2/E_ExtTrigger/Online/Data/";
static char serverName[]           = "ldas.ligo-la.caltech.edu";

/*  web page with avaliable triggers */
static char webPage[] = " http://www.ligo.caltech.edu/~wk_sz/triggers/triggers.xml";
static char webPage2[] = "http://www.uoregon.edu/~ileonor/ligo/a4/grb/online/currentgrbs_html_A4.txt";

/* where does the prog run? LLO or LHO */
static char ifo[] = "LLO";
static char instrName[3][3]={"L1","H1", "H2"};
static char ifoName[3][2]={"L","H","H"};

/* times and name of periods (extentable!) */
int numberPeriods=3;
long* startPeriod;
long* endPeriod;
static char namePeriod[3][3]={"S4","A4","S5"};

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

  proctable.processTable = (ProcessTable *) calloc( 1, sizeof(ProcessTable) );
  /* LALGPSTimeNow ( &status, &(proctable.processTable->start_time), &accuracy ); */
  /*populate_process_table( &status, proctable.processTable,   PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE ); */
  this_proc_param = processParamsTable.processParamsTable =  
    (ProcessParamsTable *) calloc( 1, sizeof(ProcessParamsTable) );
  

  /* first of all set default dirData */
  env=getenv("PWD");
  sprintf(message, "%s/../Data/", env);
  /*dirData=(char*)calloc( strlen(message+1), sizeof(char));*/
  strcpy(dirData, message);

  /* call the argument parse and check function */
  arg_parse_check( argc, argv);


  /*  startup */
  ce_startup();
  if (signal(SIGINT, sigint_handler) == SIG_ERR) {
    perror("signal");
    exit(1);
  }	


  /* check if to restart all */
  if (flagRestart) {
    printf(message,"RESTARTING program. Deleting all data files...");
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
    printf("--- TEST mode activated!");
    waitingTime=1;
  }
  printf("\n");


  /* start 'endless' loop */
  do {
    
    /* download latest list of external trigger */
    if (checkCGNtrigger2()) {

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
  downloading the latest CGN trigger from the web page
  and compare it with the stored list to check for new CGN triggers
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
  char message[256];
  /* char dirJob[256];*/
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
	extList[numberStoredCGN].gps=thisExt->start_time;
	extList[numberStoredCGN].gpsNano=thisExt->start_time_ns;
	extList[numberStoredCGN].status[0]=0;
	extList[numberStoredCGN].status[1]=0;
	extList[numberStoredCGN].status[2]=0;
	extList[numberStoredCGN].table=(ExtTriggerTable*)LALCalloc( 1, sizeof(ExtTriggerTable) );
	memcpy( extList[numberStoredCGN].table, thisExt, sizeof(ExtTriggerTable));

	/* getting user friendly date and time */
	gps->gpsSeconds=thisExt->start_time;
	gps->gpsNanoSeconds=thisExt->start_time_ns;
	LALGPStoUTC(&status, date, gps, &accuracy); 

	/* generate output saying that there is a new CGN trigger */
	sprintf(message,"New CGN trigger occured at %d/%d/%d %d:%d:%d (%d)",
		date->unixDate.tm_year+1900, date->unixDate.tm_mon+1,
		date->unixDate.tm_mday,date->unixDate.tm_hour,
		date->unixDate.tm_min, date->unixDate.tm_sec, 
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
checkCGNtrigger2
  downloading the latest CGN trigger from Isabels list
*******************************/
int checkCGNtrigger2( void )
{
  /* specify used variables here ... */
  int i,nifo;
  int counter=0;
  int numEvents=0;
  int index=0;
  int newTrigger=0;

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

  /* check if to get some trigger file or not */
  if (!flagTriggerFile) {
    
    /* copy test file if it is a test */
    if (flagTest) {
      system( "cp Test/triggersIsabel.list ." );
      
      /* downloading file */
    } else {
      sprintf(command, "wget %s -O triggersIsabel.list > download2.log 2>&1",
	      webPage2);
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
	extList[numberStoredCGN].status[0]=0;
	extList[numberStoredCGN].status[1]=0;
	extList[numberStoredCGN].status[2]=0;
	extList[numberStoredCGN].table=(ExtTriggerTable*)LALCalloc( 1, sizeof(ExtTriggerTable) );
	memcpy( extList[numberStoredCGN].table, thisExt, sizeof(ExtTriggerTable));

	/* getting user friendly date and time */
	gps->gpsSeconds=thisExt->start_time;
	gps->gpsNanoSeconds=thisExt->start_time_ns;
	LALGPStoUTC(&status, date, gps, &accuracy); 

	/* generate output saying that there is a new CGN trigger */
	sprintf(message,"New CGN trigger occured at %d/%d/%d %d:%d:%d (%d)",
		date->unixDate.tm_year+1900, date->unixDate.tm_mon+1,
		date->unixDate.tm_mday,date->unixDate.tm_hour,
		date->unixDate.tm_min, date->unixDate.tm_sec, 
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
  char message[256];
  int c=0;
  long gps;
  long gpsNano;
  int status0, status1, status2;

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
    fscanf(file, "%ld,%ld,%d,%d,%d", &gps, &gpsNano, &status0, &status1, &status2); 
    extList[c].gps      = gps;
    extList[c].gpsNano  = gpsNano;
    extList[c].status[0]= status0;
    extList[c].status[1]= status1;
    extList[c].status[2]= status2;
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
    fprintf( file, "%ld,%ld,%d,%d,%d\n",
	     extList[i].gps, extList[i].gpsNano, 
	     extList[i].status[0], extList[i].status[1], extList[i].status[2] );
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
  char dummy[256];
  char command[256];
  char text[256];
  char* filename;
  char savename[256];
  FILE* file;
  FILE* out;

  /* loop over IFO's to take into account */
  for ( i=startIFO;i<=endIFO;i++) {

    /* create new file containing all present inspiral trigger lists */
    if (flagRecalc) {
      /* use directory for recalculated trigger */
      sprintf( command, "ls -1 ../Trigger/%s/*.xml > namesInspiral.cel", instrName[i]);
    } else{
      /* use directory for online created triggers */
      sprintf( command,"ls -1 %s%s/*.xml > namesInspiral.cel", dirInspiral, instrName[i] );
    }
    system( command );
    
    /* read list and extract time data */
    c=0;
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
      fprintf(  out, "%s %ld %ld\n",inspiralList[i][j].filename,
		inspiralList[i][j].timeStart,inspiralList[i][j].timeLength);
    }
    fclose(out);

  } 

  /* TODO remove outcomment command some time... */
  /* removing unused file */
  /* system( "rm -f namesInspiral.cel" ); */
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
    
    /* loop over IFO's */
    for (lifo=startIFO;lifo<=endIFO;lifo++) { 

      /* check status of the current CGN triggers */
      if (extList[cgn].status[lifo]<4) {
	
	/* get pointer to corresponding E.T. structure */
	table=extList[cgn].table;
	if ( !table ) {
	  printf("WARNING in coinInspiralCGN: A ExtTriggerTable does not exist for cgn %d. " 
		 "Check program!\n",cgn);
	  continue;
	}
	
	/* get time of CGN trigger */
	timeCGN=(double)table->start_time + (double)table->start_time_ns/1.0e+9;
		
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
	    
	    /* mark corresponding Inspiral file */
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

	  /* creating final xml file and copy it */
	  createInspiralTable( lifo );
	  writeXML( lifo, cgn );

	  /* set new status */
	  extList[cgn].status[lifo]=5;

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

  return 0;
}

/********************************
createInspiralTable:
  create table with all needed inspiral triggers in it
*********************************/
int createInspiralTable( int lifo )
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
	sprintf( filename, " %s%s/%s",
		 dirInspiral, getIFOname(lifo), inspiralList[lifo][insp].filename);
      }
      
      sprintf(comment, "Trying opening file %s, lifo: %d\n",filename, lifo);
      printOut(0, comment);

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

      
      sprintf(comment, "+++ Adding %d-%d triggers from file %s \n", 
	      numTrig, numberDiscard,inspiralList[lifo][insp].filename);  
      printOut(0, comment);
    }
  }

  return 0;
}


/********************************
writeXML:
  writing XML file with all needed informations
*********************************/
int writeXML( int nifo, int cgnIndex )
{
  /* variables */
  char message[256];
  char command[256];
  char filename[256];
  char filenameOut[256];
  MetadataTable outputTable; 
  SnglInspiralTable*    thisInsp = NULL;
 
  /* create filename for output file */
  sprintf(filename, "externalTrig%9ld-%s.xml", extList[cgnIndex].gps, getIFOname(nifo));
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
  printf("%s\n",command);
  system(command);

  return 0;
}

/****************************** 
startAnalysisJob
   creates and starts a special analysis job
********************************/
void startAnalysisJob(ExternalList* eList, int nifo)
{
  int indexPeriod, code;
  int dirNumber=0;
  int numberOfSegments=0;
  char dirJob[256];
  char command[256];
  char message[256];
  int flagPrint=1;
  FILE* fileSegment;
  FILE* fileSED; 
  struct stat file;
  double d=0, sumTime=0, maxSegment=0;
  char dummy[256];

  if (nifo>0) {
    printf(serverName,"ldas.ligo-wa.caltech.edu");
  }

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
    sprintf(command, "LSCsegFind --server=%s --interferometer %s --type Science,Science_Up,Commissioning_Up --output-format segwizard --gps-start-time %9ld --gps-end-time %9ld > segment.txt 2> .notfound.cel ",
	    serverName, instrName[nifo], eList->gps-timeWindow, eList->gps+timeWindow );
  }
  /*sprintf(command, "LSCdataFind --server %s --observatory %s --gps-start-time %9ld --gps-end-time %9ld --type RDS_R_L3 | xargs FrFileRanges | gawk '{s++; print s, $1, $2, $2-$1;}' > segment.txt 2> .notfound.cel ",
    serverName, instrName[nifo], eList->gps-timeWindow, eList->gps+timeWindow );*/

  if (flagPrint) printOut(0, command);
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
      if (d>maxSegment) maxSegment=d;
      sumTime=sumTime+d;
    }
  } while (!feof(fileSegment));


  /* use a minimum segment length of 256 seconds! Corresponds to 3 segments. */
  if (maxSegment>=512) {

    /* create name of new directory for job execution */
    sprintf(dirJob,"../OnlineAnalysis/Jobs/%s%s-%9ld", instrName[nifo], namePeriod[indexPeriod], eList->gps);
    
    /* check for directory already */
    code=stat(dirJob, & file);
    if (code==0) {

      /* loop to find a directory that is NOT already existing */
      do {
	dirNumber++;
	sprintf(dirJob,"../OnlineAnalysis/Jobs/%s%s-%9ld_%d", 
		instrName[nifo], namePeriod[indexPeriod], eList->gps,dirNumber);
	code=stat( dirJob, &file );
      } while (code==0);
      
    }

    /* create directory */
    sprintf(command,"mkdir %s",dirJob);
    if (flagPrint) printOut(0, command);
    system(command);
    
    /* copy the segment file to that directory */
    sprintf(command, "mv segment.txt %s",dirJob);
    if (flagPrint) printOut(0, command);
    system("cat segment.txt");
    system(command);
    /*system("rm -f segment.txt");*/

    /* symbolic links to all needed basic files */
    sprintf(command,"cd %s; ln -s ../../ProgInspiral/Executables/* .",dirJob);
    if (flagPrint) printOut(0, command);
    system(command);

    /* calculate the number of segments */
    numberOfSegments = (int)(maxSegment/128.0) -1;
    if (numberOfSegments>15)  numberOfSegments=15;
    sprintf(message, "Longest segment found is %.0f s, using %d number of segments in lalapps_inspiral.\n",
	    maxSegment, numberOfSegments);

    /* create sed-file */
    fileSED=fopen("sed.file", "w");
    fprintf(fileSED, "s/sedifo/%s/\n",instrName[nifo]);
    fprintf(fileSED, "s/sedrun/%s/\n", namePeriod[indexPeriod]);
    fprintf(fileSED, "s/sedchannel/LSC-DARM_ERR/\n");
    fprintf(fileSED, "s/sednumberofsegments/%d/\n",(int)numberOfSegments);
    fprintf(fileSED, "s/sedservername/%s/\n",serverName);
    fclose(fileSED);

    /* create ini file for job execution */
    sprintf(command,"sed -f sed.file ../OnlineAnalysis/ProgInspiral/online.ini > %s/online.ini", 
	    dirJob);
    if (flagPrint) printOut(0, command);
    system(command);    

   
    /* create the DAG */ 
    sprintf(command,"cd %s; ./lalapps_inspiral_online_pipe --config-file online.ini --log-path /usr1/dietz",
	    dirJob );
    if (flagPrint) printOut(0, command);
    system(command);
    
    sprintf(command,"Starting DAG for analysis in directory %s",dirJob);
    printOut(5, command);
    sprintf(command, "cd %s; condor_submit_dag online.dag",dirJob);
    if (flagPrint) printOut(0, command);
    system (command);
    
    /* set status */
    eList->status[nifo]=2;

  } else {

    /* set status */
    eList->status[nifo]=4;

    /* output */
    if (d==0) {
      sprintf(message, "No segments found for times %ld to %ld on %s\n",
	     eList->gps-timeWindow, eList->gps+timeWindow, instrName[nifo] );      
      printOut(3, message);
    } else {
      sprintf(message, 
	      "Maximum segment for times %ld to %ld on %s is only %f seconds long. No data analyzed.\n",
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

  printf("ifoname: %s,  index: %d\n",ifoname, index);
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
  printf("ifoname: %s,  index: %d\n",ifoname, index);
  fflush(stdout);
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
  char message[256];

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
  
  /* get curreent time */
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );

  sprintf(currentTime, "%d-%d-%d %d:%d:%d", timeinfo->tm_year+1900,  timeinfo->tm_mon, timeinfo->tm_mday,
	  timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);

  return currentTime;	     
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
  char message[256];

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
  startPeriod[2]=815961613; /* S5 (preliminary) */
  endPeriod[2]=999999999;

  return 1;
}


/*******************************
ce_shutdown:
   shutdown function
*******************************/
int ce_shutdown( void )
{
  char message[256];

  /* store CGN triggers for next time (status is important!!!) */
  writeCGN();

  /* creating shutdown message */
  sprintf(message, "--- Stopping coinext at %s "
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
  char message[256];
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
      {"help",                no_argument,       0,                    'h'}, 
      {"version",             no_argument,       0,                    'V'},

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
      /*dirInspiral = (CHAR *) calloc( optarg_len, sizeof(CHAR) );*/
      memcpy( dirInspiral, optarg, optarg_len );
      ADD_PROCESS_PARAM( "string" , "%s", optarg);	      
      break;
      
    case 'A': { 
      
      /* get name for directory to put the data in */
      optarg_len = strlen( optarg ) + 1;
      /*dirData = (CHAR *) calloc( optarg_len, sizeof(CHAR) );*/
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
      /* read IFO site */
      optarg_len = strlen( optarg ) + 1;
      /*ifo = (CHAR *) calloc( optarg_len, sizeof(CHAR) );*/
      memcpy( ifo, optarg, optarg_len );	
      ADD_PROCESS_PARAM( "string" , "%s", optarg);	      	      
      break;
    }
      
    case 'c': {
      /* read triger file */
      optarg_len = strlen( optarg ) + 1;
      /*nameInputList2 = (CHAR *) calloc( optarg_len, sizeof(CHAR) );*/
      memcpy( nameInputList2, optarg, optarg_len );	
      ADD_PROCESS_PARAM( "string" , "%s", optarg);
      flagTriggerFile=1;	      	      
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
  
  return 0;
}






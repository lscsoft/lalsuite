#include <lal/LALStdlib.h>
#include "trump.h"

INT4 lalDebugLevel = LALMSGLVL3;

static int getline(char *line, int max, FILE *fpin)
{
  int i;
  for (i=0;i<MAXSTR;i++) { *(line+i) = '\0' ; }
  if (fgets(line, max, fpin) == NULL)
    return 0;
  else 
    return strlen(line);
}

/*
 *
 * The main program
 *
 */

int main(int argc, char **argv)
{
    int    arg=1;                          /* counters                             */
    int    countsamples=0;                 /* number of vetoes applied             */
    char  *vetofile = NULL;                /* name of ascii veto metadata file     */
    char  *outfile = NULL;                 /* name of ascii outfile                */
    char  *trigfile = NULL;                /* name of candidate xml file           */
    char  *timefile = NULL;                /* name of candidate xml file           */
    FILE  *fp = NULL;                      /* generic file pointer                 */
    FILE  *fpout = NULL;                      /* generic file pointer                 */
    char   vetoline[1024];
    int    i,pass,coincidence_flag=0;
    vetoParams *veto_entry=NULL, *thisentry=NULL, *preventry=NULL;
    candParams candidates;
    char  *vcol[128],*eventfile[32];
    timeWindow *thiswindow=NULL,*vwindows=NULL, *awindows=NULL;
    double timeVetoed;
    candEvent *eventhead=NULL, *thisCEvent=NULL;
    candEvent *Ieventhead=NULL, *thisIEvent=NULL;
    float time_analyzed, coincidence_window=0.01, safety=0.0;
    int **triggerHistogram, **injectHistogram, numbins=40;

    /*******************************************************************
    * PARSE ARGUMENTS (arg stores the current position)               *
    *******************************************************************/

    if (argc <= 1){
        fprintf(stderr,  USAGE, *argv );
        return 0;
    }

    while ( arg < argc ) {
        /*********************************************************
        * Output filename. 
        *********************************************************/
        if ( !strcmp( argv[arg], "-o" ) ) {
            if ( argc > arg + 1 ) {
                arg++;
                outfile = argv[arg++];
            }else{
                fprintf(stderr,  USAGE, *argv );
                return RESPONSEC_EARG;
            }
        }
        /*********************************************************
        * File containing veto information metadata 
        *********************************************************/
        else if ( !strcmp( argv[arg], "--help" ) ) {
                fprintf(stderr,  USAGE, *argv );
                return 0;
        }
        /*********************************************************
        * File containing veto information metadata 
        *********************************************************/
        else if ( !strcmp( argv[arg], "--veto" ) ) {
            if ( argc > arg + 1 ) {
                arg++;
                vetofile = argv[arg++];
            }else{
                fprintf(stderr,  USAGE, *argv );
                return RESPONSEC_EARG;
            }
        }
        /*********************************************************
        * Name of input file with trigger parameters
        *********************************************************/
        else if ( !strcmp( argv[arg], "--trigger" ) ) {
            if ( argc > arg + 1 ) {
                arg++;
                trigfile = argv[arg++];
            }else{
                fprintf(stderr,  USAGE, *argv );
                return RESPONSEC_EARG;
            }
        }
        /*********************************************************
        * Name of file with time intervals
        *********************************************************/
        else if ( !strcmp( argv[arg], "--times" ) ) {
            if ( argc > arg + 1 ) {
                arg++;
                timefile = argv[arg++];  /* file with times analyzed */
                safety = atof(argv[arg++]);    /* secs ignored per chunk   */
            }else{
                fprintf(stderr,  USAGE, *argv );
                return RESPONSEC_EARG;
            }
        }
        /*********************************************************
        * Set-up for coincidences
        *********************************************************/
        else if ( !strcmp( argv[arg], "--coincidence" ) ) {
            if ( argc > arg + 2 ) {
                arg++;
                eventfile[0] = argv[arg++];
                eventfile[1] = argv[arg++];
                coincidence_window = atof(argv[arg++]);
                coincidence_flag=1;
                fprintf(stderr,"Doing coincidence, all other options ignored\n");
            }else{
                fprintf(stderr,  USAGE, *argv );
                return RESPONSEC_EARG;
            }
        }
        /* Check for unrecognized options. */
        else if ( argv[arg][0] == '-' ) {
            fprintf(stderr,  USAGE, *argv );
            return RESPONSEC_EARG;
        }
    } /* End of argument parsing loop. */


    /****************************************************************
     * Open the log file
     ****************************************************************/
    if ( ( fpout = fopen( "trump.log", "w" ) ) == NULL ) {
        fprintf(stderr,  "Error opening log file\n");
        return RESPONSEC_EFILE;
    }


    /****************************************************************
     * If we're doing a coincidence analysis,  this is the code that
     * does it ....... it uses the output of running the code to
     * determine which events survive veto,  but only looks at events
     * which are present ......
     ****************************************************************/
    if (coincidence_flag){
        FILE *fpevent[32];
        candEvent *cevents[32], myevent;
        int dummyI;

        fprintf(fpout,"Doing coincidence analysis\n"); fflush(fpout);

        /************************************************************
        * First we open the files containing the events from each
        * instrument.   These files have a format int, time, snr,
        * chisq ........
        ************************************************************/
        for(i=0;i<2;i++){
            if ( ( fpevent[i] = fopen( eventfile[i], "r" ) ) == NULL ) {
                fprintf(stderr,  "Error opening file %s\n", eventfile[i] );
                return RESPONSEC_EFILE;
            }
        }

        /************************************************************
        * Get the events from the first file and construct a list of
        * events .....
        ************************************************************/
        cevents[0]=NULL;
        while ( getline(vetoline, MAXSTR, fpevent[0]) ){
            /* skip lines containing "#" */
            if ( !(strstr(vetoline, "#")) ) {
                if ( cevents[0] == NULL){
                    /* Allocate head of list first time */
                    thisCEvent = cevents[0] = (candEvent *)malloc(sizeof(candEvent));
                    (*thisCEvent).next_event = NULL;
                }
                else 
                {
                    /* Otherwise allocate next node */
                    (*thisCEvent).next_event = 
                        (candEvent *)malloc(sizeof(candEvent));
                    thisCEvent = (*thisCEvent).next_event;
                    (*thisCEvent).next_event = NULL;
                }
                /**************************************************
                * Note the formating of the data ............
                **************************************************/
                sscanf(vetoline,"%i %lf %f %f\n",
                        &dummyI,
                        &((*thisCEvent).time),
                        &((*thisCEvent).snr),
                        &((*thisCEvent).chisq));
            }
        }


        /************************************************************
         * Get the events from the second file and find coincidences
         * with the original list .....
         ************************************************************/
        while ( getline(vetoline, MAXSTR, fpevent[1]) ){
                 sscanf(vetoline,"%i %lf %f %f\n",
                         &dummyI,
                         &(myevent.time),
                         &(myevent.snr),
                         &(myevent.chisq));
                 thisCEvent = cevents[0];
                 while ( thisCEvent != NULL ){
                     /***********************************************
                      * If there is a coincidence,  then we print out
                      * the event information ......
                      ***********************************************/
                     if ( myevent.time > (*thisCEvent).time - coincidence_window &&
                             myevent.time < (*thisCEvent).time + coincidence_window){
                         fprintf(stdout,"%f %f %f %f %f %f\n",
                                (*thisCEvent).time, (*thisCEvent).snr, 
                                (*thisCEvent).chisq, myevent.time, 
                                myevent.snr, myevent.chisq);
                     }
                     thisCEvent = (*thisCEvent).next_event;
                 }
        }
      
        /***********************************************************
         * Close the files
         ***********************************************************/
        for(i=0;i<2;i++){
            fclose(fpevent[i]);
        }
        fclose(fpout);

        /*************************************************************
         * This code is not tested as of 28 oct 2002
         ************************************************************/
        thisCEvent=cevents[0];
        while (thisCEvent != NULL){
           candEvent *nextEvent;

           nextEvent=(*thisCEvent).next_event;
           free(thisCEvent);
           thisCEvent=nextEvent;
        }
        exit(0);
    }
    
    /****************************************************************
     * First read determine the total time analyzed
     ****************************************************************/
    if ( ( fp = fopen( timefile, "r" ) ) == NULL ) {
                fprintf(stderr,  "Error opening timefile\n" );
        return RESPONSEC_EFILE;
    }

    fprintf(fpout,"# Determining time analyzed\n"); fflush(fpout);
    time_analyzed=0.0;
    while ( getline(vetoline, MAXSTR, fp) ){
        double dummy;
        float interval;
        if ( ! (awindows) ){
            thiswindow = awindows = (timeWindow *)malloc(sizeof(timeWindow));
            (*thiswindow).next_time_window = NULL;
        }
        else
        {
            (*thiswindow).next_time_window = (timeWindow *)malloc(sizeof(timeWindow));
            thiswindow = (*thiswindow).next_time_window;
            (*thiswindow).next_time_window = NULL;
        }
        sscanf(vetoline, "%lf %f\n", &dummy, &interval);
        (*thiswindow).start_time = dummy;
        (*thiswindow).end_time = (dummy + interval);
        time_analyzed+=interval;
    }

    fclose(fp);
        
    
    /****************************************************************
     * Open the veto configuration file, read it, count no. of vetoes
    ****************************************************************/
    if ( vetofile != NULL ){

        if ( ( fp = fopen( vetofile, "r" ) ) == NULL ) {
            fprintf(stderr,  "Error opening vetofile\n" );
            return RESPONSEC_EFILE;
        }

        countsamples=0;
        while ( getline(vetoline, MAXSTR, fp) ){
            /* skip lines containing "#" */
            if ( !(strstr(vetoline, "#")) ) {
                countsamples++;

                if (veto_entry == NULL){
                    /* Allocate head of list first time */
                    thisentry = veto_entry = (vetoParams *)malloc(sizeof(vetoParams));
                    (*thisentry).next_veto = NULL;
                }
                else 
                {
                    /* Otherwise allocate next node */
                    preventry=thisentry;
                    (*thisentry).next_veto = (vetoParams *)malloc(sizeof(vetoParams));
                    thisentry = (*thisentry).next_veto;
                    (*thisentry).next_veto = NULL;
                }

                /* copy the veto information into the structure */
                vcol[0]=strtok(vetoline,";");
                for( i=1; i<5 && (vcol[i]=strtok(NULL,";")) ; i++) { }

                /***********************************************************
                *
                * The information about the vetoes is:
                *          name:  any string
                *          filename: xml file containing veto triggers
                *          table_column: column to use
                *          threshold: number above which to make veto
                *          dtime: how much time to veto
                *          
                ***********************************************************/
                strcpy(((*thisentry).name), vcol[0]);
                strcpy(((*thisentry).filename), vcol[1]);
                strcpy(((*thisentry).table_column), vcol[2]);
                (*thisentry).threshold = atof(vcol[3]);
                (*thisentry).dtime = atof(vcol[4]);
            }
        }

        fclose(fp);


        /******************************************************************** 
        * Construct a list of vetoed times
        ********************************************************************/
        fprintf(fpout,"# Determining vetoes\n"); fflush(fpout);
        thisentry = veto_entry;
        while ( thisentry != NULL ){
            if ( (buildVetoTimes( thisentry )) != 0 ){
                fprintf(stderr,"Error building veto time list\n");
                exit(1);
            }

            /******************************************************************** 
            * Determine vetoed times
            ********************************************************************/
            thiswindow = (*thisentry).vwindows;
            timeVetoed = 0.0;
            while ( thiswindow != NULL ){
                timeVetoed += ((*thiswindow).end_time - (*thiswindow).start_time);
                thiswindow = (*thiswindow).next_time_window;
            }


            /******************************************************************** 
            * Print out the information that was supplied in the veto files
            ********************************************************************/
            fprintf(fpout,"Veto filename = %s\n",(*thisentry).filename);
            fprintf(fpout,"%s > %f, dtime = +/-%f\n",
                    (*thisentry).table_column,
                    (*thisentry).threshold,
                    (*thisentry).dtime);
            fprintf(fpout,"Time vetoed %f, %f%% \n----\n", timeVetoed, 100.0 * 
                    timeVetoed/time_analyzed);
            fflush(fpout);
            thisentry = (*thisentry).next_veto;
        }

        /******************************************************** 
        * Resolve all vetoed times into a single list
        ********************************************************/
        resolveVetoTimes( &vwindows, veto_entry);
        thiswindow = vwindows;
        timeVetoed = 0.0;
        while ( thiswindow != NULL ){
            timeVetoed += ((*thiswindow).end_time - (*thiswindow).start_time);
            thiswindow = (*thiswindow).next_time_window;
        }
        fprintf(fpout,"Total time vetoed %f, %f%%\n----\n", timeVetoed, 100.0 *
                timeVetoed/time_analyzed);

    }

    /******************************************************** 
    * Read information about the candidate event list
    ********************************************************/
    if ( ( fp = fopen( trigfile, "r" ) ) == NULL ) {
                fprintf(stderr,  "Error opening candidate file\n" );
        return RESPONSEC_EFILE;
    }


    countsamples=0;
    while ( getline(vetoline, MAXSTR, fp) ){
        /* skip lines containing "#" */
        if ( !(strstr(vetoline, "#")) ) {
            countsamples++;

            /* copy the candidate information into the structure */
            vcol[0]=strtok(vetoline,";");
            for( i=1; i<7 && (vcol[i]=strtok(NULL,";")) ; i++) { }

             /***********************************************************
             *
             * The information about the triggers is:
             *          name:  any string
             *          injectfile: xml file containing injections
             *          triggerfile: xml file containing triggers
             *          snr_threshold: threshold on SNR
             *          chi_threshold: threshold on CHISQ
             *          dtime: ????
             *          
             ***********************************************************/
            strcpy((candidates.name), vcol[0]);
            strcpy((candidates.triggerfile), vcol[1]);
            strcpy((candidates.injectfile), vcol[2]);
            candidates.snr_threshold = atof(vcol[3]);
            candidates.chi_threshold = atof(vcol[4]);
            candidates.dtime = atof(vcol[5]);
            strcpy((candidates.injepochs), vcol[6]);
        }
    }

    fclose(fp);

    
    /******************************************************************** 
    * Print out the information that was supplied in the cadidate files
    ********************************************************************/
    fprintf(fpout,"Trigger filename = %s\n",(candidates).triggerfile);
    fprintf(fpout,"Injection filename = %s\n",(candidates).injectfile);
    fprintf(fpout,"SNR > %f, CHISQ < %f, dtime = +/-%f\n----\n",
            (candidates).snr_threshold,
            (candidates).chi_threshold,
            (candidates).dtime);
    fflush(fpout);

    
    /******************************************************************** 
    * Read in list of candidate events and apply vetoes to them.
    ********************************************************************/
    buildEventList( &eventhead, vwindows, candidates, TRIGGERS);
    triggerHistogram = calloc( numbins, sizeof(int*) );
    for ( i = 0; i < numbins; ++i )
    {
        triggerHistogram[i] = calloc( numbins, sizeof(int) );
    }
    build2DHistogram( eventhead , "2d-event.dat", triggerHistogram, numbins);

    
    /******************************************************************** 
    * Print out the list of events that were found
    ********************************************************************/
    if ( ( fp = fopen( "triggers.dat", "w" ) ) == NULL ) {
                fprintf(stderr,  USAGE, *argv );
        return RESPONSEC_EFILE;
    }
    i=0;
    thisCEvent = eventhead;
    while ((thisCEvent) != NULL){
        i++;
        fprintf(fp,"%i %f %f %f\n",i,(*thisCEvent).time,
                (*thisCEvent).snr,(*thisCEvent).chisq);
        thisCEvent = (*thisCEvent).next_event;
    }
    fclose(fp);


    /******************************************************************** 
    * Read in list of injection events and apply vetoes to them.
    ********************************************************************/
    buildEventList( &Ieventhead, vwindows, candidates, INJECTIONS);
    injectHistogram = calloc( numbins, sizeof(int*) );
    for ( i = 0; i < numbins; ++i )
    {
        injectHistogram[i] = calloc( numbins, sizeof(int) );
    }   
    build2DHistogram( Ieventhead , "2d-inject.dat", injectHistogram, numbins);

    
    /******************************************************************** 
    * Print out the list of events that were found
    ********************************************************************/
    if ( ( fp = fopen( "injections.dat", "w" ) ) == NULL ) {
                fprintf(stderr,  USAGE, *argv );
        return RESPONSEC_EFILE;
    }
    i=0;
    thisIEvent = Ieventhead;
    while ((thisIEvent) != NULL){
        i++;
        fprintf(fp,"%i %f %f %f\n",i,(*thisIEvent).time,
                (*thisIEvent).snr,(*thisIEvent).chisq);
        thisIEvent = (*thisIEvent).next_event;
    }
    fclose(fp);

    /******************************************************** 
    * Read information about the injections
    ********************************************************/
 /*   if ( ( fp = fopen( candidates.injepochs, "r" ) ) == NULL ) {
                fprintf(stderr,  USAGE, *argv );
        return RESPONSEC_EFILE;
    }
*/
    /******************************************************** 
    * Count the number of injections that survive
    ********************************************************/
  /*  countsamples=0;
    while ( getline(vetoline, MAXSTR, fp) ){
        int injtimeS,injtimeNS;
  */      /* skip lines containing "#" */
 /*       if ( !(strstr(vetoline, "#")) ) {

            sscanf(vetoline,"%i %i\n",&injtimeS,&injtimeNS);
            thiswindow=awindows;
            while ( thiswindow != NULL ){
                if ((double)injtimeS > (*thiswindow).start_time &&
                        (double)injtimeS < (*thiswindow).end_time){
                    pass = 1;
                    break;
                }
                thiswindow = (*thiswindow).next_time_window;
            }
            thiswindow=vwindows;
            while ( thiswindow != NULL ){
                if ((double)injtimeS > (*thiswindow).start_time &&
                        (double)injtimeS < (*thiswindow).end_time){
                    pass = 0;
                    break;
                }
                thiswindow = (*thiswindow).next_time_window;
            }
        }
        countsamples+=pass;
    }
    fclose(fp);
   */

    fprintf(fpout,"ninject:=%i;\n", countsamples);
    fprintf(fpout,"T:=%1.2f/3600\n", (time_analyzed - timeVetoed));
    fflush(fpout);
    
    fclose(fpout);

    return 0;
}

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/Thresholds.h>
#include <lal/LALStdlib.h>
#include <lalapps.h>
#include "metaio.h"
#include "event_utils.h"
#include "donald.h"
RCSID("$Id$");

static int getline(char *line, int max, FILE *fpin)
{
    int i;
    for (i=0;i<MAXSTR;i++) { *(line+i) = '\0' ; }
    if (fgets(line, max, fpin) == NULL)
        return 0;
    else 
        return strlen(line);
}


/**********************************************************************************
 *
 * Read in the information which defines an single interferometer
 * search
 *
 *********************************************************************************/
int readconfig(snglIFO *ifo, FILE *fp, FILE *logfile)
{
    char   vetoline[1024],*vcol[128];
    int    first=1,i;
    const char *ifoname;
    vetoParams *thisentry=NULL;
    timeWindow *thiswindow=NULL;

    while ( getline(vetoline, MAXSTR, fp) ){
        /* skip lines containing "#" */
        if ( !(strstr(vetoline, "#")) ) {

            /* copy the veto information into the structure */
            vcol[0]=strtok(vetoline,";");
            if (first){
                ifoname = vcol[0];
                strncpy(ifo->ifoname,ifoname,2);
                first =0;
            }else if ( (strcmp(vcol[0],ifoname)) ){
                fprintf(stderr,"Configuration file has incorrect IFO\n");
                return -1;
            }

            /* find out what information this is */
            vcol[1]=strtok(NULL,";");

            /********************************************************************
             *
             * Determine the veto information 
             *
             *******************************************************************/
            if ( !strcmp(vcol[1],"veto") )
            {
                fprintf(logfile,"# Determining veto information for %s\n",vcol[0]);
                fflush(logfile);
                if (ifo->veto_entry == NULL){
                    /* Allocate head of list first time */
                    thisentry = ifo->veto_entry = (vetoParams *)malloc(sizeof(vetoParams));
                    (*thisentry).next_veto = NULL;
                    (*thisentry).vwindows = NULL;
                }
                else 
                {
                    /* Otherwise allocate next node */
                    (*thisentry).next_veto = (vetoParams *)malloc(sizeof(vetoParams));
                    thisentry = (*thisentry).next_veto;
                    (*thisentry).next_veto = NULL;
                    (*thisentry).vwindows = NULL;
                }

                /* copy the veto information into the structure */
                for( i=2; i<8 && (vcol[i]=strtok(NULL,";")) ; i++) { }

                /***********************************************************
                 *
                 * The information about the vetoes is:
                 *          name:  any string
                 *          filename: xml file containing veto triggers
                 *          table_column: column to use
                 *          threshold: number above which to make veto
                 *          minusdtime: how much time to veto before
                 *          plusdtime: how much time to veto after
                 *
                 ***********************************************************/
                strcpy(((*thisentry).name), vcol[0]);
                strcpy(((*thisentry).filename), vcol[2]);
                strcpy(((*thisentry).table_column), vcol[3]);
                (*thisentry).threshold = atof(vcol[4]);
                (*thisentry).minusdtime = atof(vcol[5]);
                (*thisentry).plusdtime = atof(vcol[6]);
                (*thisentry).ratio = atof(vcol[7]);
            }
            /********************************************************************
             *
             * Determine the trigger information 
             *
             *******************************************************************/
            else if ( !strcmp(vcol[1],"triggers") )
            {
                fprintf(logfile,"# Determining trigger information fopr %s\n",vcol[0]);
                fflush(logfile);

                /* copy the candidate information into the structure */
                for( i=2; i<8 && (vcol[i]=strtok(NULL,";")) ; i++) { }

                /***********************************************************
                 *
                 * The information about the triggers is:
                 *          name:  any string
                 *          injectfile: xml file containing injections
                 *          triggerfile: xml file containing triggers
                 *          snr_threshold: threshold on SNR
                 *          chi_threshold: threshold on CHISQ
                 *          dtime: cluster window
                 *          
                 ***********************************************************/
                strcpy((ifo->candidates.name), vcol[0]);
                strcpy((ifo->candidates.triggerfile), vcol[2]);
                strcpy((ifo->candidates.injectfile), vcol[3]);
                ifo->candidates.snr_threshold = atof(vcol[4]);
                ifo->candidates.chi_threshold = atof(vcol[5]);
                ifo->candidates.dtime = atof(vcol[6]);
                strcpy((ifo->candidates.injepochs), vcol[7]);
            }
            /********************************************************************
             *
             * Determine the trigger information 
             *
             *******************************************************************/
            else if ( !strcmp(vcol[1],"dqbit") )
            {
                fprintf(logfile,"# Determining data quality bit for %s\n",vcol[0]);
                fflush(logfile);

                /* copy the candidate information into the structure */
                for( i=2; i<3 && (vcol[i]=strtok(NULL,";")) ; i++) { }

                ifo->dqbit = (int)(pow(2.0,atof(vcol[2])));
            }
            /********************************************************************
             *
             * Determine the times analyzed 
             *
             *******************************************************************/
            else if ( !strcmp(vcol[1],"times") )
            {
                FILE *fpintmp;

                fprintf(logfile,"# Determining time analyzed\n"); fflush(logfile);

                /* Open the file with the segment information */
                for( i=2; i<4 && (vcol[i]=strtok(NULL,";")) ; i++) { }

                if ( ( fpintmp = fopen( vcol[2], "r" ) ) == NULL ) {
                    fprintf(stderr,  "Error opening timefile\n" );
                    return RESPONSEC_EFILE;
                }

                /* determine the time analyzed */
                ifo->safety=atof(vcol[3]);
                ifo->time_analyzed=0.0;
                while ( getline(vetoline, MAXSTR, fpintmp) ){
                    double dummys,dummye,ovrlap=ifo->safety;
                    int    doneflag, segnum;

                    sscanf(vetoline, "%i %lf %lf %i\n", 
                            &doneflag, &dummys, &dummye, &segnum);
                    if (doneflag == 1) {
                        if ( ! (ifo->awindows) ){
                            thiswindow = ifo->awindows = 
                                (timeWindow *)malloc(sizeof(timeWindow));
                            (*thiswindow).next_time_window = NULL;
                        }
                        else
                        {
                            (*thiswindow).next_time_window = 
                                (timeWindow *)malloc(sizeof(timeWindow));
                            thiswindow = (*thiswindow).next_time_window;
                            (*thiswindow).next_time_window = NULL;
                        }

                        (*thiswindow).start_time = dummys+ovrlap/2.0;
                        (*thiswindow).end_time = dummye-ovrlap/2.0;
                        ifo->end_time = (*thiswindow).end_time;
                        ifo->time_analyzed+=(dummye-dummys-ovrlap);
                    }
                }
                ifo->start_time = ifo->awindows->start_time;
                fclose(fpintmp);
            }
        }
    }
    return 0;
}



/************************************************************************************
 *
 * The main program
 *
 ***********************************************************************************/

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
    char   vetoline[1024],filename[1024];
    int    i,pass=0,numIFO=1, dummyMask, j, printRaw=0;
    vetoParams *veto_entry=NULL, *thisentry=NULL;
    candParams candidates;
    timeWindow *thiswindow=NULL,*vwindows=NULL, *awindows=NULL;
    double timeVetoed=0.0;
    candEvent *thisCEvent=NULL, *myevent=NULL, *thisIEvent=NULL;
    float time_analyzed, delm=DELM;
    int **triggerHistogram, **injectHistogram, numbins=40, *coincident_times;
    snglIFO ifo[2];
    double dummyStart=0,dummyEnd=0,coincidence_window=COINWINDOW;
    float distance=H1L1DISTANCE;
    multiInspiral *multInspEv;
    static LALStatus        status;
    lal_errhandler = LAL_ERR_EXIT;

  set_debug_level( "3" );

    /*******************************************************************
     * PARSE ARGUMENTS (arg stores the current position)               *
     *******************************************************************/

    if (argc <= 1){
        fprintf(stderr,  USAGE, *argv );
        return 0;
    }

    while ( arg < argc ) {
        /*********************************************************
         * File containing veto information metadata 
         *********************************************************/
        if ( !strcmp( argv[arg], "--help" ) ) {
            fprintf(stderr,  USAGE, *argv );
            return 0;
        }
        /*********************************************************
         * Name of file with time intervals
         *********************************************************/
        else if ( !strcmp( argv[arg], "--ifofiles" ) ) {
            if ( argc > arg + 1 ) {
                arg++;
                numIFO=0;
                while (argv[arg] && (argv[arg][0] != '-'))
                {
                    ifo[numIFO].cfgfile = argv[arg++];
                    numIFO++;
                }
            }else{
                fprintf(stderr,  USAGE, *argv );
                return RESPONSEC_EARG;
            }
        }
        else if ( !strcmp( argv[arg], "--delm" ) ) {
            if ( argc > arg + 1 ) {
                arg++;
                delm = atof(argv[arg++]);
            }else{
                fprintf(stderr,  USAGE, *argv );
                return RESPONSEC_EARG;
            }
        }
        else if ( !strcmp( argv[arg], "--distance" ) ) {
            if ( argc > arg + 1 ) {
                arg++;
                distance = atof(argv[arg++]);
            }else{
                fprintf(stderr,  USAGE, *argv );
                return RESPONSEC_EARG;
            }
        }
        else if ( !strcmp( argv[arg], "--dtime" ) ) {
            if ( argc > arg + 1 ) {
                arg++;
                coincidence_window = atof(argv[arg++]);
            }else{
                fprintf(stderr,  USAGE, *argv );
                return RESPONSEC_EARG;
            }
        }
        else if ( !strcmp( argv[arg], "--raw" ) ) {
            if ( argc > arg + 1 ) {
                arg++;
                printRaw = 1;
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
     * 
     * Read in the information about each of the interferometers
     * 
     ****************************************************************/
    for ( i=0 ; i<numIFO ; i++ ){

        float calfudge=1.0;
        
        if ( ( fp = fopen( ifo[i].cfgfile, "r" ) ) == NULL ) {
            fprintf(stderr,  "Error opening IFO file %s\n", ifo[i].cfgfile );
            return RESPONSEC_EFILE;
        }
        ifo[i].veto_entry=NULL;
        ifo[i].awindows=NULL;
        ifo[i].vwindows=NULL;
        ifo[i].eventhead=NULL;
        readconfig(&ifo[i], fp, fpout);
        fclose(fp);

        fprintf(fpout,"# Working on IFO %s out of %i\n",ifo[i].ifoname,numIFO);
        fflush(fpout);

        /******************************************************************** 
         * Construct a list of vetoed times
         ********************************************************************/
        fprintf(fpout,"# Determining vetoes\n"); fflush(fpout);
        thisentry = ifo[i].veto_entry;
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
            fprintf(fpout,"%s > %f, dtime = -%f, +%f\n",
                    (*thisentry).table_column,
                    (*thisentry).threshold,
                    (*thisentry).minusdtime,
                    (*thisentry).plusdtime);
            fprintf(fpout,"Time vetoed %f of %f, %f%% \n----\n", timeVetoed,
                    ifo[i].time_analyzed , 100.0 * timeVetoed/ifo[i].time_analyzed);
            fflush(fpout);
            thisentry = (*thisentry).next_veto;
        }

        /******************************************************** 
         * Resolve all vetoed times into a single list
         ********************************************************/
        resolveVetoTimes( &(ifo[i].vwindows), ifo[i].veto_entry);
        thiswindow = ifo[i].vwindows;
        timeVetoed = 0.0;
        while ( thiswindow != NULL ){
            timeVetoed += ((*thiswindow).end_time - (*thiswindow).start_time);
            thiswindow = (*thiswindow).next_time_window;
        }
        fprintf(fpout,"Total time vetoed %f, %f%%\n----\n", timeVetoed, 100.0 *
                timeVetoed/ifo[i].time_analyzed);

        /******************************************************************** 
         * Print out the information that was supplied in the trigger files
         ********************************************************************/
        fprintf(fpout,"Trigger filename = %s\n",ifo[i].candidates.triggerfile);
        fprintf(fpout,"Injection filename = %s\n",ifo[i].candidates.injectfile);
        fprintf(fpout,"SNR > %f, CHISQ < %f, dtime = %f\n----\n",
                ifo[i].candidates.snr_threshold,
                ifo[i].candidates.chi_threshold,
                ifo[i].candidates.dtime);
        fflush(fpout);

        /******************************************************************** 
         * Read in list of candidate events and apply vetoes to them.
         ********************************************************************/
        fprintf(fpout,"# Building a list of triggers\n"); fflush(fpout);
        buildEventList( &(ifo[i].eventhead), ifo[i].vwindows, 
                ifo[i].candidates, TRIGGERS, 0, 1.0);

        /******************************************************************** 
         * Print out the list of events that were found
         ********************************************************************/
        if(printRaw){
            sprintf(filename,"triggers-%s.dat",ifo[i].ifoname);
            if ( ( fp = fopen( filename, "w" ) ) == NULL ) {
                fprintf(stderr,  USAGE, *argv );
                return RESPONSEC_EFILE;
            }
            j=0;
            thisCEvent = ifo[i].eventhead;
            while ((thisCEvent) != NULL){
                j++;
                fprintf(fp,"%i %f %f %f %f %f\n",j,(*thisCEvent).time,
                        (*thisCEvent).snr,(*thisCEvent).chisq,
                        (*thisCEvent).eff_distance,(*thisCEvent).mchirp);
                thisCEvent = (*thisCEvent).next_event;
            }
            fclose(fp);
        }

        /******************************************************************** 
         * Read in list of injection events and apply vetoes to them.
         ********************************************************************/
        fprintf(fpout,"# Building a list of injection triggers\n"); fflush(fpout);
        if (!(strcmp(ifo[i].ifoname,"H1"))){
            calfudge = 1.11/0.87;
        }
        fprintf(fpout,"# Using calibration fudge factor %f for %s\n",
                calfudge, ifo[i].ifoname);
        buildEventList( &(ifo[i].Ieventhead), ifo[i].vwindows, 
                ifo[i].candidates, INJECTIONS, 0, calfudge);

        /******************************************************************** 
         * Print out the list of events that were found
         ********************************************************************/
        if(printRaw){
            sprintf(filename,"injections-%s.dat",ifo[i].ifoname);
            if ( ( fp = fopen( filename, "w" ) ) == NULL ) {
                fprintf(stderr,  USAGE, *argv );
                return RESPONSEC_EFILE;
            }
            j=0;
            thisIEvent = ifo[i].Ieventhead;
            while ((thisIEvent) != NULL){
                j++;
                fprintf(fp,"%i %f %f %f %f %f\n",j,(*thisIEvent).time,
                        (*thisIEvent).snr,(*thisIEvent).chisq,
                        (*thisIEvent).eff_distance,(*thisIEvent).mchirp);
                thisIEvent = (*thisIEvent).next_event;
            }
            fclose(fp);
        }
    }

    /*************************************************************************
     *
     * Determine coincident times
     *
     *************************************************************************/
    fprintf(fpout,"# Determining coincident times\n"); fflush(fpout);
    buildDataQaulity(&coincident_times, ifo, numIFO, &dummyStart, &dummyEnd);


    /*************************************************************************
     *
     * Search for events in multiple instruments
     *
     *************************************************************************/
    fprintf(fpout,"# Looking for coincident triggers\n"); fflush(fpout);
    LALSortTriggers(&status,ifo,numIFO);
    fprintf(fpout,"# Delta M is %e \n",delm); fflush(fpout);
    buildMultiInspiralEvents(&multInspEv, coincident_times, ifo, numIFO, 
            TRIGGERS, dummyStart, delm, distance, coincidence_window);

    fp = fopen("triggers.dat","w");
    printInspiralEvents(fp,&ifo[0], 3, TRIGGERS);
    /* Identify singles in each of the insterferometers */
    dummyMask=1;
    for ( i=0 ; i<numIFO ; i++ ){
        printInspiralEvents(fp,&ifo[i], dummyMask, TRIGGERS);
        dummyMask *= 2;
    }
    fclose(fp);

    fprintf(fpout,"# Looking for coincident injections\n"); fflush(fpout);
    buildMultiInspiralEvents(&multInspEv, coincident_times, ifo, numIFO, 
            INJECTIONS, dummyStart, delm, distance, coincidence_window);

    fp = fopen("injections.dat","w");
    printInspiralEvents(fp,&ifo[0], 3, INJECTIONS);
    /* Identify singles in each of the insterferometers */
    dummyMask=1;
    for ( i=0 ; i<numIFO ; i++ ){
        printInspiralEvents(fp,&ifo[i], dummyMask, INJECTIONS);
        dummyMask *= 2;
    }
    fclose(fp);

    /******************************************************** 
     * Read information about the injections
     ********************************************************/
    if ( ( fp = fopen( ifo[0].candidates.injepochs, "r" ) ) == NULL ) {
        fprintf(stderr,  USAGE, *argv );
        return RESPONSEC_EFILE;
    }

    /******************************************************** 
     * Count the number of injections that survive
     ********************************************************/
    fprintf(fpout,"# Determining the number of injections that were made\n"); fflush(fpout);
    countsamples=0;
    while ( getline(vetoline, MAXSTR, fp) ){
        double injtimeS,injtimeNS;
        pass = 0; 
        /* skip lines containing "#" */
        if ( !(strstr(vetoline, "#")) ) {

            sscanf(vetoline,"%lf",&injtimeS);
            if ( injtimeS >= dummyStart && injtimeS <= dummyEnd ){
                if ( coincident_times[(int)( 10*(injtimeS - dummyStart))] > 0 ){
                    pass = 1;
                }
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


    /**********************************************************
     * Determine coincident times and times analyzed
     *********************************************************/
    {
        int   count, numPts= 10 * (int)(dummyEnd-dummyStart), prev_value;
        FILE *tmpPtr;
        double time_done=0.0;

        tmpPtr = fopen("coin.times","w");
        prev_value = coincident_times[0];
        fprintf(tmpPtr,"%lf %i %lf\n",dummyStart,coincident_times[0],time_done);
        for(i=1;i<numPts;i++){
            if ( coincident_times[i] != prev_value ){
                fprintf(tmpPtr,"%lf %i %lf\n",dummyStart+0.1*(double)(i-1),
                        prev_value,time_done);
                fprintf(tmpPtr,"%lf %i %lf\n",dummyStart+0.1*(double)i,
                        coincident_times[i],time_done);
            }
            if ( coincident_times[i] > 0 ){
                count++;
            }
            prev_value = coincident_times[i];
        }
        fprintf(tmpPtr,"%lf %i %lf\n",dummyStart+0.1*(double)(i-1),coincident_times[i-1],
                time_done);
        fclose(tmpPtr);
        time_analyzed = 0.1*count;
    }

    fprintf(fpout,"ninject:=%i;\n", countsamples);
    fprintf(fpout,"T = %1.2f sec = %1.2f hr\n", 
            (time_analyzed),
            (time_analyzed)/3600.0);
    fprintf(fpout,"# Finished\n");
    fflush(fpout);


    /**********************************************************
     * Compute the upper limit
     *********************************************************/
    {
        FILE  *ulout,*compout;
        float *snrbin, snrdiv, minsnr=8.0,maxsnr=minsnr+10.0,*efficiency,*numevents;
        double mu;
        int    numbins=800,k,first;
        Chi2ThresholdIn  thresholdIn;

        ulout = fopen("ul.dat","w");

        /* set up the bin boundaries for the histogram */
        snrbin     = (float *) calloc( numbins+1, sizeof(float));
        efficiency = (float *) calloc (numbins+1, sizeof(float));
        numevents  = (float *) calloc (numbins+1, sizeof(float));
        snrdiv     = ( maxsnr - minsnr ) / (float) numbins;
        for ( i = 0; i <= numbins; ++i )
        {
            snrbin[i] = minsnr + (float) i * snrdiv;
        }

        for (i=0 ; i<numIFO ; i++){
            for ( thisCEvent = ifo[i].eventhead; thisCEvent; 
                    thisCEvent = thisCEvent->next_event )
            {
                if ( thisCEvent->candidate ){
                    for ( k = 0; k < numbins; ++k )
                    {
                        if ( thisCEvent->snr >= snrbin[k] ) numevents[k]+=1.0;
                    }
                }
            }
        }

        for (i=0 ; i<numIFO ; i++){
            for ( thisCEvent = ifo[i].Ieventhead; thisCEvent; 
                    thisCEvent = thisCEvent->next_event )
            {
                if ( thisCEvent->candidate ){
                    for ( k = 0; k < numbins; ++k )
                    {
                        if ( thisCEvent->snr >= snrbin[k] ) efficiency[k]+=1.0;
                    }
                }
            }
        }
        fprintf(ulout,"# snr numevents efficiency\n");

        compout = fopen("comparisons.txt","a");
        first=1;
        for ( k = 0; k < numbins; ++k )
        {
            thresholdIn.falseAlarm=0.1;
            thresholdIn.dof=2.0*(numevents[k]+1);
            LALChi2Threshold(&status, &mu, &thresholdIn);
            mu /= 2.0;
            fprintf( ulout, "%f %f %f %e %e %e\n", snrbin[k],numevents[k],
                    efficiency[k]/(float)countsamples,
                    (3600.0*countsamples*(float)mu)/
                    (efficiency[k]*(time_analyzed)), mu,
                    thresholdIn.dof);
            if ( numevents[k] == 0 && first ){
                fprintf(compout,"%f\t%e\t%f\t%f\n",coincidence_window,
                        delm,distance,efficiency[k]/(float)countsamples);
                first = 0;
            }
        }

        fclose(compout);
        fclose(ulout);
        free(snrbin);
        free(efficiency);
        free(numevents);
        return 0;
    }

    fclose(fpout);

    return 0;
}

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include "metaio.h"
#include "event_utils.h"
#include "trump.h"

#define OVRLAP 32.0

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
    int         arg=1;                  /* counters                             */
    char       *outfile = NULL;         /* name of ascii outfile                */
    FILE       *fpout = NULL;           /* generic file pointer                 */
    char        vetoline[1024];
    int         i;
    char       *eventfile[32];
    candEvent  *thisCEvent=NULL;
    float       coincidence_window=0.01;
    FILE       *fpevent[32];
    candEvent  *cevents[32], myevent;
    int         dummyI;

    /*******************************************************************
     * PARSE ARGUMENTS (arg stores the current position)               *
     *******************************************************************/

    if (argc <= 1){
        fprintf(stderr,  ICOINUSAGE, *argv );
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
                fprintf(stderr,  ICOINUSAGE, *argv );
                return RESPONSEC_EARG;
            }
        }
        /*********************************************************
         * File containing veto information metadata 
         *********************************************************/
        else if ( !strcmp( argv[arg], "--help" ) ) {
            fprintf(stderr,  ICOINUSAGE, *argv );
            return 0;
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
                fprintf(stderr,"Doing coincidence, all other options ignored\n");
            }else{
                fprintf(stderr,  ICOINUSAGE, *argv );
                return RESPONSEC_EARG;
            }
        }
        /* Check for unrecognized options. */
        else if ( argv[arg][0] == '-' ) {
            fprintf(stderr,  ICOINUSAGE, *argv );
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
            sscanf(vetoline,"%i %lf %f %f %f %f\n",
                    &dummyI,
                    &((*thisCEvent).time),
                    &((*thisCEvent).snr),
                    &((*thisCEvent).chisq),
                    &((*thisCEvent).eff_distance),
                    &((*thisCEvent).mchirp)
                  );
        }
    }


    /************************************************************
     * Get the events from the second file and find coincidences
     * with the original list .....
     ************************************************************/
    while ( getline(vetoline, MAXSTR, fpevent[1]) ){
        sscanf(vetoline,"%i %lf %f %f %f %f\n",
                &dummyI,
                &(myevent.time),
                &(myevent.snr),
                &(myevent.chisq),
                &(myevent.eff_distance),
                &(myevent.mchirp)
              );
        thisCEvent = cevents[0];
        while ( thisCEvent != NULL ){
            /***********************************************
             * If there is a coincidence,  then we print out
             * the event information ......
             ***********************************************/
            if ( myevent.time > (*thisCEvent).time - coincidence_window &&
                    myevent.time < (*thisCEvent).time + coincidence_window){
                fprintf(stdout,"%f %f %f %f %f %f %f %f\n",
                        (*thisCEvent).time-myevent.time,
                        (*thisCEvent).mchirp-myevent.mchirp,
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
    return(0);

}



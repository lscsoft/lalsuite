#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/Thresholds.h>
#include <lal/LALStdlib.h>
#include <lalapps.h>
#include "metaio.h"
RCSID("$Id$");

#define MAXSTR 2048
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

/* Usage format string. */
#define USAGE "Usage: %s --input infile [--threshold threshold] \
    --freq fstart fstop df [--help]\n"

#define POWERC_ENORM  0
#define POWERC_ESUB   1
#define POWERC_EARG   2
#define POWERC_EVAL   3
#define POWERC_EFILE  4
#define POWERC_EINPUT 5
#define POWERC_EMEM   6

#define POWERC_MSGENORM  "Normal exit"
#define POWERC_MSGESUB   "Subroutine failed"
#define POWERC_MSGEARG   "Error parsing arguments"
#define POWERC_MSGEVAL   "Input argument out of valid range"
#define POWERC_MSGEFILE  "Could not open file"
#define POWERC_MSGEINPUT "Error reading file"
#define POWERC_MSGEMEM   "Out of memory"

#define TRUE  1
#define FALSE 0

/************************************************************************************
 *
 * The main program
 *
 ***********************************************************************************/

int main(int argc, char **argv)
{
    static LALStatus         stat;
    INT4                     retVal=0;
    INT4                     inarg=1;
    BOOLEAN                  freqHistogram=TRUE;
    BOOLEAN                  tfHistogram=FALSE;
    REAL4                    fstart=0.0,fstop=1000.0,df=0.0;
    REAL4                    threshold=1.0,logThreshold=0.0;
    FILE                     *fpin=NULL,*fpTFHist=NULL;
    INT4                     fileCounter=0;

    CHAR                     inputfile[MAXSTR],line[MAXSTR];

    INT4                      iTriggerConfidence=0;
    INT4                      iTriggerFreq=0;
    INT4                      iTriggerBand=0;
    INT4                      i;
    REAL4                     tmpConfidence=0.0;
    REAL4                     tmpFreq=0.0;
    REAL4                     bandwidth=0.0;

    INT4                      nFreqBins=0;
    INT8                      *freqHistData=NULL;
    INT8                      *tmpHistData=NULL;
    size_t len;

    struct MetaioParseEnvironment triggerParseEnv;
    const MetaioParseEnv triggerEnv = &triggerParseEnv;

    /*******************************************************************
     * BEGIN PARSE ARGUMENTS (inarg stores the current position)        *
     *******************************************************************/
    if (argc <= 1){
        LALPrintError( USAGE, *argv );
        return 0;
    }

    while ( inarg < argc ) {
        /* Parse output file option. */
        if ( !strcmp( argv[inarg], "--freq" ) ) {
            freqHistogram = TRUE;
            inarg++;
            if (argv[inarg] && (argv[inarg][0] != '-'))
            {
                fstart = atof(argv[inarg++]);
                if (argv[inarg] && (argv[inarg][0] != '-'))
                {
                    fstop = atof(argv[inarg++]);
                    if (argv[inarg] && (argv[inarg][0] != '-'))
                    {
                        df = atof(argv[inarg++]);
                    }
                }
            }
        }
        /* threshold (in probability) to apply */
        else if ( !strcmp( argv[inarg], "--input" ) )
        {
            inarg++;
            sprintf(inputfile,"%s",argv[inarg++]);
        }
        /* time-frequency histogram */
        else if ( !strcmp( argv[inarg], "--tfhist" ) )
        {
            inarg++;
            fpTFHist=fopen(argv[inarg++],"w");
            tfHistogram=TRUE;
        }
        /* threshold (in probability) to apply */
        else if ( !strcmp( argv[inarg], "--threshold" ) )
        {
            inarg++;
            threshold = atof(argv[inarg++]);
            logThreshold = log(threshold);
        }
        /* print a help message */
        else if ( !strcmp( argv[inarg], "--help" ) ) {
            LALPrintError( USAGE, *argv );
            return POWERC_EARG;
        }
        /* Check for unrecognized options. */
        else if ( argv[inarg][0] == '-' ) {
            LALPrintError( "Unrecognized options\n" );
            return POWERC_EARG;
        }
        /* Default case for unknown arguments */
        else
        {
            LALPrintError( "Unknown arguments\n" );
            return POWERC_EARG;
        }
    }
    if (inputfile == NULL){
        LALPrintError( "Must supply an xml file to parse\n" );
        return POWERC_EARG;
    }

    /*******************************************************************
    * END PARSE ARGUMENTS                                              *
    *******************************************************************/


    /*******************************************************************
    * INITIALIZE ERROR HANDLER
    *******************************************************************/
    lal_errhandler = LAL_ERR_DFLT;
    set_debug_level( "3" );

    /*****************************************************************
     * INITIALIZE THE HISTOGRAM
     *****************************************************************/
    if ( freqHistogram ){
        nFreqBins = (INT4)( (fstop-fstart)/df );
        freqHistData = (INT8 *) LALCalloc(nFreqBins , sizeof(INT8));
        tmpHistData = (INT8 *) LALCalloc(nFreqBins , sizeof(INT8));
    }

    if (tfHistogram){
        fprintf(fpTFHist,"# fstart=%f, fstop=%f, df=%f\n",fstart,fstop,df);
    }

    /*****************************************************************
     * OPEN FILE WITH LIST OF XML FILES (one file per line)
     ****************************************************************/
    if ( !(fpin = fopen(inputfile,"r")) ){
        LALPrintError("Could not open input file\n");
    }

    while ( getline(line, MAXSTR, fpin) ){
        fileCounter++;
        fprintf(stderr,"Processing file %i %s\n",fileCounter,line);
    
        /******************************************************************
         * OPEN XML FILE AND READ IT
         *****************************************************************/
        if ( (retVal = MetaioOpen( triggerEnv, line)) !=0 ){
            fprintf(stderr, "Error opening injection file %s\n", line );
            MetaioAbort( triggerEnv ); 
            exit(2);
        }

        /* Locate the relevant columns */
        iTriggerConfidence = MetaioFindColumn( triggerEnv, "CONFIDENCE");
        iTriggerFreq       = MetaioFindColumn( triggerEnv, "CENTRAL_FREQ" );
        iTriggerBand       = MetaioFindColumn( triggerEnv, "BANDWIDTH" );

        /* Loop over the triggers */
        memset(tmpHistData, 0, nFreqBins*sizeof(INT8));
        while (1) {

            /* assume candidate is an event unless proven otherwise */
            retVal = MetaioGetRow(triggerEnv);
            if ( retVal == -1 ) {
                printf( "Error while getting row from injection or trigger file\n");
                MetaioAbort( triggerEnv ); 
                return 6;
            } else if ( retVal == 0 ) {
                /*-- Reached end of file --*/
                break;
            }

            /* get the confidence from the table */
            tmpConfidence = 
                triggerEnv->ligo_lw.table.elt[iTriggerConfidence].data.real_4;  
            tmpFreq       = 
                triggerEnv->ligo_lw.table.elt[iTriggerFreq].data.real_4;  
            bandwidth       = 
                triggerEnv->ligo_lw.table.elt[iTriggerBand].data.real_4;  

            /* require the confidence to exceed threshold */
            if ( tmpConfidence < logThreshold ){
                if (freqHistogram){
                    i=(INT4)((tmpFreq-fstart-0.5*bandwidth)/df);
                    while (i<(INT4)((tmpFreq-fstart+0.5*bandwidth)/df) &&
                            i < nFreqBins){
                        freqHistData[i] += 1;
                        tmpHistData[i] += 1;
                        i++;
                    }
                }
            }
        }
        /********************************************************************
         * SPIT OUT THE RESULTS
         *******************************************************************/
        if (freqHistogram){
            FILE *fpFreqHist=NULL;
            fpFreqHist = fopen("freq-hist.txt","w");
            for(i=0;i<nFreqBins;i++){
                fprintf(fpFreqHist,"%f %i\n",fstart+i*df+0.5*df,
                        freqHistData[i]);
            }
            fclose(fpFreqHist);
        }

        /********************************************************************
         * histograms as function of "time"
         *******************************************************************/
        if (tfHistogram){
            INT4 i;
            for(i=0;i<nFreqBins;i++){
                fprintf(fpTFHist,"%i %f %i\n",fileCounter,
                        fstart+ i*df + 0.5*df, tmpHistData[i]);
            }
            fprintf(fpTFHist,"\n");
        }


        MetaioAbort(triggerEnv);
    }

    LALFree( freqHistData );
    LALFree( tmpHistData );
    if (tfHistogram){
        fclose(fpTFHist);
    }

    LALCheckMemoryLeaks();

    return 0;
}

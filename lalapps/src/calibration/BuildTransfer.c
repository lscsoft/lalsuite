#include <math.h>
#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/Calibration.h>
#include <lal/AVFactories.h>
#include <lal/LALDatatypes.h>
#include <lal/LALConstants.h>
#include <lal/Interpolate.h>

NRCSID( RESPONSEC, "$Id$" );


#define MAXSTR 256
#define PPOLE  0.76      /* Pendulum frequency in Hz                         */
#define FMIN   20.0     /* Lowest frequency filled with correct calibration */

#define RESPONSEC_ENORM  0
#define RESPONSEC_ESUB   1
#define RESPONSEC_EARG   2
#define RESPONSEC_EVAL   3
#define RESPONSEC_EFILE  4
#define RESPONSEC_EINPUT 5
#define RESPONSEC_EMEM   6

#define RESPONSEC_MSGENORM  "Normal exit"
#define RESPONSEC_MSGESUB   "Subroutine failed"
#define RESPONSEC_MSGEARG   "Error parsing arguments"
#define RESPONSEC_MSGEVAL   "Input argument out of valid range"
#define RESPONSEC_MSGEFILE  "Could not open file"
#define RESPONSEC_MSGEINPUT "Error reading file"
#define RESPONSEC_MSGEMEM   "Out of memory"

/* Usage format string. */
#define USAGE "Usage: %s [-f npoints fmin fmax] [-o outfile]\n" \
                        "[-i ilwdfile] [-r respfile] [-a amplitude] [-t] \n"

/* Macros for printing errors and testing subroutines. */
#define ERROR( code, msg, statement )                                \
do                                                                   \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
		 "        %s %s\n", (code), *argv, __FILE__,         \
		 __LINE__, RESPONSEC, statement ? statement : \
                 "", (msg) );                                        \
}                                                                    \
while (0)

typedef struct 
tagfresponse
{
  REAL4                freq;     /* frequency          */
  REAL4                a;        /* amplitude          */
  REAL4                da;       /* error in amplitude */
  REAL4                phi;      /* phase              */
  REAL4                dphi;     /* error in phase     */
  INT4                 ntrials;  /* number of averages */
  struct tagfresponse *nextfresponse; 
}
fresponse;

static UINT4 mylocate(REAL4 *farray, REAL4 target, UINT4 n)
{
  UINT4 i=0;
  UINT4 jmin;
  REAL4 dummy;

  jmin = 0;
  dummy = abs(farray[jmin]-target);
  for (i=1 ; i<n ; i++){
    if ( abs(farray[i]-target) < dummy ){
      dummy = abs(farray[i]-target);
      jmin = i;
    }
  }
      
  return jmin;
}

INT4 lalDebugLevel = LALMSGLVL3;

static const CHAR *ilwdTop="<?ilwd?>\n";

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
    static LALStatus stat;
    INT4   arg=1;                          /* counters                             */
    INT4   countsamples=0,npoints=10000;   /* descriptors                          */
    CalibrationRecord calrec;              /* place-holder for calibration record  */
    CHAR  *infile = NULL;                  /* name of ascii input file             */
    CHAR  *outfile = NULL;                 /* name of ascii outfile                */
    CHAR  *ilwdfile = NULL;                /* name of ilwd outfile                 */
    FILE  *fp = NULL;                      /* generic file pointer                 */
    REAL4  fmax=1024.0;                    /* max freq                             */
    REAL4  fmin=0,freq;
    REAL4  amplitude=0.0;
    REAL4 *farray = NULL, *amparray = NULL, *phasearray = NULL;
    INT4   response=1;                     /* compute response if 1, transfer if 0 */
    INT4   i,gpsTime=0;
    UINT4  location;
    CHAR   respline[MAXSTR];
    fresponse *fresponse_entry = NULL, *thisentry = NULL, *preventry = NULL;
    SInterpolateOut  intOutput;
    SInterpolatePar  intParams;

    /*  Create memory for the transfer function */
    calrec.transfer = (COMPLEX8FrequencySeries *)LALMalloc( 
            sizeof(COMPLEX8FrequencySeries));

    /*******************************************************************
    * PARSE ARGUMENTS (arg stores the current position)               *
    *******************************************************************/

    if (argc <= 1){
        LALPrintError( USAGE, *argv );
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
                ERROR( RESPONSEC_EARG, RESPONSEC_MSGEARG, 0 );
                LALPrintError( USAGE, *argv );
                return RESPONSEC_EARG;
            }
        }
        /*********************************************************
         * ILWD filename. 
         *********************************************************/
        else if ( !strcmp( argv[arg], "-i" ) ) {
            if ( argc > arg + 1 ) {
                arg++;
                ilwdfile = argv[arg++];
            }else{
                ERROR( RESPONSEC_EARG, RESPONSEC_MSGEARG, 0 );
                LALPrintError( USAGE, *argv );
                return RESPONSEC_EARG;
            }
        }
        /*********************************************************
         * Name of input file with reponse function
         *********************************************************/
        else if ( !strcmp( argv[arg], "-r" ) ) {
            if ( argc > arg + 1 ) {
                arg++;
                infile = argv[arg++];
            }else{
                ERROR( RESPONSEC_EARG, RESPONSEC_MSGEARG, 0 );
                LALPrintError( USAGE, *argv );
                return RESPONSEC_EARG;
            }
        }
        /*********************************************************
         * Output a transfer function instead of response
         *********************************************************/
        else if ( !strcmp( argv[arg], "-t" ) ) {
            arg++;
            response = 0;
        }
        /*********************************************************
         * Parse the frequency range of the output
         *********************************************************/
        else if ( !strcmp( argv[arg], "-f" ) ) {
            if ( argc > arg + 2 ) {
                arg++;
                npoints = atoi( argv[arg++] );
                calrec.transfer->data=NULL;
                LALCCreateVector(&stat, &(calrec.transfer->data), npoints);
                fmin = atof( argv[arg++] );
                fmax = atof( argv[arg++] );
                calrec.transfer->deltaF = (fmax-fmin) / ( (REAL4)(npoints-1) );
                calrec.transfer->f0 = fmin;
            }else{
                LALPrintError( USAGE, *argv );
                return RESPONSEC_EARG;
            }
        }
        /*********************************************************
         * What is the amplitude normalization
         *********************************************************/
        else if (!strcmp( argv[arg], "-a") ) {
            arg++;
            /* Factor of two to account for light travel path */
            amplitude = atof( argv[arg++] ); 
        }
        /* Check for unrecognized options. */
        else if ( argv[arg][0] == '-' ) {
            LALPrintError( USAGE, *argv );
            return RESPONSEC_EARG;
        }
    } /* End of argument parsing loop. */


    /****************************************************************
     * Tell the operator about some of the hardwired numbers
     ****************************************************************/
    fprintf(stderr,"\nThe following are hardwired into the BuildTransfer program:\n");
    fprintf(stderr,"\t Pendulum frequency = %f\n",PPOLE);
    fprintf(stderr,"\t Below %f Hz,  response set to unity\n\n",FMIN);


    /****************************************************************
     * Open the transfer function file 
     ****************************************************************/
    if ( ( fp = fopen( infile, "r" ) ) == NULL ) {
        ERROR( RESPONSEC_EFILE, RESPONSEC_MSGEFILE, infile );
        return RESPONSEC_EFILE;
    }


    /****************************************************************
     * read in the file and count the number of frequency samples 
     ****************************************************************/
    countsamples=0;
    while ( getline(respline, MAXSTR, fp) ){
        /* skip lines containing "#" */
        if ( !(strstr(respline, "#")) ) {
            countsamples++;

            if (fresponse_entry == NULL){
                /* Allocate head of list first time */
                thisentry = fresponse_entry = (fresponse *)malloc(sizeof(fresponse));
                (*thisentry).nextfresponse = NULL;
            }
            else 
            {
                /* Otherwise allocate next node */
                preventry=thisentry;
                (*thisentry).nextfresponse = (fresponse *)malloc(sizeof(fresponse));
                thisentry = (*thisentry).nextfresponse;
                (*thisentry).nextfresponse = NULL;
            }

            /* copy the calibration information */
            sscanf(respline," %f %f %f %f %f %i", &((*thisentry).freq),
                    &((*thisentry).a),
                    &((*thisentry).da),
                    &((*thisentry).phi),
                    &((*thisentry).dphi),
                    &((*thisentry).ntrials));

            /* Do a check for duplicate frequencies */
            if( preventry != NULL ){
                if ( (*preventry).freq == (*thisentry).freq ){
                    /* keep entry with smallest measurement error */
                    if ( pow((*preventry).dphi/(*preventry).phi,2.0) >
                            pow((*thisentry).dphi/(*thisentry).phi,2.0)){
                        (*preventry) =  (*thisentry);
                    }
                    /* free the node that was not needed */
                    free(thisentry);
                    thisentry = preventry;
                    (*thisentry).nextfresponse = NULL;
                    countsamples--;
                }
            }

        }
        else if ( (strstr(respline, "GPS")) ){
            sscanf(respline,"#     ETMX  GPS: %i",&gpsTime);
        } 
    }


    /******************************************************************** 
     * Allocate and fill arrays for interpolation 
     ********************************************************************/
    farray     = (REAL4 *)LALMalloc( countsamples*sizeof(REAL4) );
    amparray   = (REAL4 *)LALMalloc( countsamples*sizeof(REAL4) );
    phasearray = (REAL4 *)LALMalloc( countsamples*sizeof(REAL4) );
    thisentry = fresponse_entry;
    i=0;
    for (i=0 ; i<countsamples ; i++){
        farray[i]     = (*thisentry).freq;
        amparray[i]   = (*thisentry).a;
        phasearray[i] = (*thisentry).phi;
        thisentry     = (*thisentry).nextfresponse;
    } 


    /******************************************************************** 
     * Do the interpolation
     ********************************************************************/
    intParams.n=4;
    location = 0;
    for (i=0 ; i<npoints ; i++){
        REAL4 amp;
        freq = fmin + ((float) i)*calrec.transfer->deltaF; 
        if ( freq >= FMIN )
        {
            location=mylocate(farray,freq,countsamples);
            if (location > (countsamples-intParams.n)){
                location = (countsamples-intParams.n);
            } else {
                location-=(intParams.n/2);
            }
            intParams.x = &farray[location];

            /******************************************************** 
             * interpolate the AMPLITUDE and fold in the actuator 
             * response and absolute calibration 
             ********************************************************/
            intParams.y = &amparray[location];
            LALSPolynomialInterpolation(&stat, &intOutput, freq, &intParams); 
            amp = PPOLE*PPOLE * amplitude / ((-freq*freq + PPOLE*PPOLE) * intOutput.y);

            /********************************************************
             * interpolate the PHASE
             ********************************************************/
            intParams.y = &phasearray[location];
            LALSPolynomialInterpolation(&stat, &intOutput, freq, &intParams); 

            /******************************************************** 
             * fill the response function array 
             ********************************************************/
            calrec.transfer->data->data[i].re = amp * cos(intOutput.y);
            calrec.transfer->data->data[i].im = - amp * sin(intOutput.y);
        }
        else {
            /********************************************************
             * Below FMIN,  we fill the response function array with 
             * unity 
             ********************************************************/
            calrec.transfer->data->data[i].re = 1.0 ;
            calrec.transfer->data->data[i].im = 0.0 ;
        }
    }

    /* Compute the transfer function */
    /* LALComputeTransfer(&stat, &calrec); */

    /******************************************************** 
    * Output to an ascii file in the format:  freq Re[T] Im[T]
    ********************************************************/
    if ( outfile ) {
        if ( ( fp = fopen( outfile, "w" ) ) == NULL ) {
            ERROR( RESPONSEC_EFILE, RESPONSEC_MSGEFILE, outfile );
            return RESPONSEC_EFILE;
        }
        for(i=0;i<npoints;i++){
            REAL4 re, im;
            re = calrec.transfer->data->data[i].re;
            im = calrec.transfer->data->data[i].im;
            fprintf(fp,"%e %e %e\n",fmin+ calrec.transfer->deltaF * i, re, im);
        }
        fclose( fp );
    }
    /********************************************************* 
    * otherwise print it to an ilwd file
    ********************************************************/
    else if (ilwdfile ) {
        if ( ( fp = fopen( ilwdfile, "w" ) ) == NULL ) {
            ERROR( RESPONSEC_EFILE, RESPONSEC_MSGEFILE, ilwdfile );
            return RESPONSEC_EFILE;
        }
        fprintf(fp,ilwdTop );
        fprintf(fp," <ilwd name='response::sequence' size='9'>\n");
        fprintf(fp,"  <lstring name='complex:domain' size='4'>FREQ</lstring>\n");
        fprintf(fp,"  <int_4u name='gps_sec:start_time' units='sec'>%i</int_4u>\n",gpsTime);
        fprintf(fp,"  <int_4u name='gps_nan:start_time' units='nanosec'>%i</int_4u>\n",0);
        fprintf(fp,"  <int_4u name='gps_sec:stop_time' units='sec'>%i</int_4u>\n",gpsTime);
        fprintf(fp,"  <int_4u name='gps_nan:stop_time' units='nanosec'>%i</int_4u>\n",0);
        fprintf(fp,"  <real_8 name='start_freq' units='hz'>%16.16e</real_8>\n", fmin);
        fprintf(fp,"  <real_8 name='stop_freq' units='hz'>%16.16e</real_8>\n", fmax);
        fprintf(fp,"  <real_8 name='freq:step_size' units='hz'>%16.16e</real_8>\n",
                calrec.transfer->deltaF);
        if ( response ){
            fprintf(fp,"  <complex_8 dims='%i' name='data' units='strain/count'>",npoints);
        }
        else {
            fprintf(fp,"  <complex_8 dims='%i' name='data' units='counts/attostrain'>",
                    npoints);
        }

        for(i=0;i<npoints;i++){
            REAL4 re, im, tmpAmplitude;
            re = calrec.transfer->data->data[i].re;
            im = calrec.transfer->data->data[i].im;
            if ( response ){
                fprintf(fp,"%e %e ", re, im);
            }
            else {
                re *= 1.0e18;
                im *= 1.0e18;
                tmpAmplitude = sqrt(re*re + im*im);
                fprintf(fp,"%e %e ", re/tmpAmplitude, -im/tmpAmplitude);
            }
        }
        fprintf(fp,"  </complex_8>\n");
        fprintf(fp,"</ilwd>");
        fclose(fp);
    }


    return 0;
}

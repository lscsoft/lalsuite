#include <math.h>
#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/Calibration.h>
#include <lal/AVFactories.h>
#include <lal/LALDatatypes.h>
#include <lal/LALConstants.h>

NRCSID( RESPONSEC, "$Id$" );

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
#define USAGE "Usage: %s [-z nz z_1 ... z_nz] [-p np p_1 .... p_np] \n" \
                        "[-g gain] [-f npoints fmin fmax] [-o outfile]\n" \
                        "[-i ilwdfile] [-t] \n"

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

INT4 lalDebugLevel = LALMSGLVL3;

static const CHAR *ilwdTop="<?ilwd?>\n"
"    <ilwd name='response::sequence' size='9'>\n"
"        <lstring name='complex:domain' size='4'>FREQ</lstring>\n"
"        <int_4u name='gps_sec:start_time' units='sec'>60000000</int_4u>\n"
"        <int_4u name='gps_nan:start_time' units='nanosec'>0</int_4u>\n"
"        <int_4u name='gps_sec:stop_time' units='sec'>60000128</int_4u>\n"
"        <int_4u name='gps_nan:stop_time' units='nanosec'>0</int_4u>\n";

int main(int argc, char **argv)
{
  static LALStatus stat;
  INT4 arg=1,i;                   /* counters                            */
  INT4 nz=4,np=5,npoints=10000;   /* descriptors                         */
  CalibrationRecord calrec;       /* place-holder for calibration record */
  CHAR *outfile = NULL;           /* name of ascii outfile */
  CHAR *ilwdfile = NULL;          /* name of ilwd outfile */
  FILE *fp = NULL;                /* generic file pointer */
  REAL4 fmax=1024.0;              /* max freq */
  REAL4 fmin;
  INT4  response=1;               /* compute response if 1, transfer if 0 */

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
    /* Parse output file option. */
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
    else if ( !strcmp( argv[arg], "-t" ) ) {
	arg++;
	response = 0;
    }
    /* Parse zeros option. */
    else if ( !strcmp( argv[arg], "-z" ) ) {
      arg++;
      nz = atoi( argv[arg] );
      if ( nz == 0 ){
        calrec.zeros = (REAL8Vector *) LALMalloc( sizeof(REAL8Vector) );
        calrec.zeros->length=nz;
        arg++;
      }
      else if ( argc > arg + nz ) {
        calrec.zeros=NULL;
        LALDCreateVector(&stat, &calrec.zeros, nz);
        for (i=0 ; i<nz ; i++)
        {
          arg++;
          calrec.zeros->data[i] = atof( argv[arg] );
        }
        arg++;
      }else{
        LALPrintError( USAGE, *argv );
        return RESPONSEC_EARG;
      }
    }
    /* Parse poles option. */
    else if ( !strcmp( argv[arg], "-p" ) ) {
      arg++;
      np = atoi( argv[arg] );
      if ( np == 0 ){
        calrec.poles = (REAL8Vector *) LALMalloc( sizeof(REAL8Vector) );
        calrec.poles->length=np;
        arg++;
      }
      else if ( argc > arg + np ) {
        calrec.poles=NULL;
        LALDCreateVector(&stat, &calrec.poles, np);
        for (i=0 ; i<np ; i++)
        {
          arg++;
          calrec.poles->data[i] = atof( argv[arg] );
        }
        arg++;
      }else{
        LALPrintError( USAGE, *argv );
        return RESPONSEC_EARG;
      }
    }
    /* Parse gain option. */
    else if ( !strcmp( argv[arg], "-g" ) ) {
      arg++;
      calrec.gain=atof( argv[arg++] );
    }
    /* Parse output option. */
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
    /* Check for unrecognized options. */
    else if ( argv[arg][0] == '-' ) {
      LALPrintError( USAGE, *argv );
      return RESPONSEC_EARG;
    }
  } /* End of argument parsing loop. */
  /* Compute the transfer function */
  LALComputeTransfer(&stat, &calrec);
 
  /* print to an ascii text file.  Format:  freq Re[T] Im[T]  */
  if ( outfile ) {
    if ( ( fp = fopen( outfile, "w" ) ) == NULL ) {
      ERROR( RESPONSEC_EFILE, RESPONSEC_MSGEFILE, outfile );
      return RESPONSEC_EFILE;
    }
    for(i=0;i<npoints;i++){
      REAL4 re, im;
      re = calrec.transfer->data->data[i].re;
      im = calrec.transfer->data->data[i].im;
      fprintf(fp,"%e %e %e\n",fmin+calrec.transfer->deltaF * i, re, im);
    }
    fclose( fp );
  }
  /* otherwise print it to an ilwd file */
  else if (ilwdfile ) {
    if ( ( fp = fopen( ilwdfile, "w" ) ) == NULL ) {
      ERROR( RESPONSEC_EFILE, RESPONSEC_MSGEFILE, ilwdfile );
      return RESPONSEC_EFILE;
    }
    fprintf( fp, ilwdTop );
    fprintf( fp, "        <real_8 name='start_freq' units='hz'>%16.16e</real_8>\n",
        0.0);
    fprintf( fp, "        <real_8 name='stop_freq' units='hz'>%16.16e</real_8>\n",
        fmax);
    fprintf( fp, "        <real_8 name='freq:step_size' units='hz'>%16.16e</real_8>\n",
        calrec.transfer->deltaF);
    if ( response ){
      fprintf( fp, "        <complex_8 dims='%i' name='data' units='strain/count'>",
          npoints);
    }
    else {
      fprintf( fp, "        <complex_8 dims='%i' name='data' units='counts/attostrain'>",
          npoints);
    }

    for(i=0;i<npoints;i++){
      REAL4 re, im;
      re = calrec.transfer->data->data[i].re;
      im = calrec.transfer->data->data[i].im;
      fprintf(fp,"%e %e ", re, im);
    }
    fprintf( fp, "</complex_8>\n");
    fprintf( fp, "    </ilwd>");
    fclose( fp );
  }

         
  return 0;
}

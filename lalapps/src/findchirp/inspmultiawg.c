/*
 * inspmultiawg.c
 * Author:  Creighton, T and Fairhurst, S
 * $Id$
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <lal/Date.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/SeqFactories.h>
#include <lal/Units.h>
#include <lal/Random.h>
#include <lal/VectorOps.h>
#include <lal/Inject.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/StreamInput.h>
#include <lal/TimeFreqFFT.h>
#include <lalapps.h>
#include <processtable.h>

RCSID( "$Id$" );

#define PROGRAM_NAME "inspmultiawg"
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

#define INSPAWGFILEC_ENORM  0
#define INSPAWGFILEC_ESUB   1
#define INSPAWGFILEC_EARG   2
#define INSPAWGFILEC_EVAL   3
#define INSPAWGFILEC_EFILE  4
#define INSPAWGFILEC_EINPUT 5
#define INSPAWGFILEC_EMEM   6

#define INSPAWGFILEC_MSGENORM  "Normal exit"
#define INSPAWGFILEC_MSGESUB   "Subroutine failed"
#define INSPAWGFILEC_MSGEARG   "Error parsing arguments"
#define INSPAWGFILEC_MSGEVAL   "Input argument out of valid range"
#define INSPAWGFILEC_MSGEFILE  "Could not open file"
#define INSPAWGFILEC_MSGEINPUT "Error reading file"
#define INSPAWGFILEC_MSGEMEM   "Out of memory"

/* Default parameter settings. */
#define EPOCH (0)
#define DIST  (1.0)
#define M1    (1.4)
#define M2    (1.4)
#define INC   (0.0)
#define PHIC  (0.0)
#define LENGTH (64)
#define NPT   (1048576)
#define FREQ  (16384)
#define DT    (1.0/16384.0)

#define TRUE  1
#define FALSE 0

/* Global constants. */
#define MSGLEN (256) /* maximum length of warning/info messages */
#define FSTART (40.0) /* initial frequency of waveform */
#define FSTOP (3000.0) /* final frequency of waveform */
#define DELTAT (0.00006103515625) /* sampling interval of amplitude and phase */

/* Usage format string. */
#define USAGE \
  "%s [options]\n"\
"  --help               display this message\n"\
"  --source sfile       source file containing details of injection\n"\
"  --actuation actfile  file containing the actuation function\n"\
"  --darm2inj dcfactor  calibration between darm_ctrl and injection point\n"\
"  --summary sumfile    write details of injections to file\n"\
"  --ifo ifo            name of interferomter (optional)\n"\
"  --flow fstart        start frequency of injection (default 40 Hz)\n"\
"  --fhigh fstop        end frequency of injection (default: end at ISCO)\n"\
"  --smooth Qfac        ringdown the end of the injection with Q factor Qfac\n"\
"  --length length      length of the data (default 64 seconds)\n"\
"  --samplerate	freq    rate at which data is sampled (default 16384Hz)\n"\
"  --debug-level debug  the lal debug level\n"\
"  --user-tag tag       user-tag added to output file names\n"\
"\n"


/* Macros for printing errors and info subroutines. */
#define ERROR( code, msg, statement )                                \
  do                                                                   \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
      "        %s %s\n", (code), *argv, __FILE__,         \
      __LINE__, rcsid, statement ? statement : \
      "", (msg) );                                        \
}                                                                    \
while (0)

#define INFO( statement )                                            \
  do                                                                   \
if ( lalDebugLevel & LALINFO )                                       \
{                                                                    \
  LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"       \
      "        %s\n", *argv, __FILE__, __LINE__,          \
      rcsid, (statement) );                    \
}                                                                    \
while (0)

/*****************************************************************************
 *
 * The main program
 *
 *****************************************************************************/

  int
main(int argc, char **argv)
{
  /* Command-line parsing variables. */
  int c;		      /* command-line argument counter */
  static int verbose_flag;
  static LALStatus stat;      /* status structure */
  CHAR	*sourcefile = NULL;   /* name of sourcefile */
  CHAR	*actfile = NULL;      /* name of respfile */
  CHAR	*summfile = NULL;      /* name of summfile */
  CHAR  *tag = NULL;	      /* user-tag */
  CHAR  ifo[3]="";	      /* name of IFO */
  CHAR  fname[256];	      /* name of outfile */
  INT4	length = LENGTH;      /* length of data segment */
  INT4	freq = FREQ;	      /* sampling frequency */
  INT4	npt = NPT;            /* number of output points */
  REAL8	dt = DT;	      /* output sampling interval */
  REAL4 fstart = FSTART;      /* start frequency */
  INT4  fstopset = FALSE;     /* check whether end frequency specified */
  REAL4 fstop  = FSTOP;	      /* stop frequency */
  INT4	xmloutput = FALSE;  
  LIGOTimeGPS inj_length; /* length of the injection */  

  /* File reading variables. */
  FILE		       *fp = NULL,*fq = NULL;  /* generic file pointer */
  BOOLEAN		ok = 1;	    /* whether input format is correct */
  UINT4			i;	    /* generic index over file lines */
  INT8			epoch;      /* epoch stored as an INT8 */
  LALLeapSecAccuracy	accuracy = LALLEAPSEC_LOOSE;
  ProcessParamsTable   *this_proc_param;
  MetadataTable	        proctable;
  MetadataTable	        procparams;
  LIGOLwXMLStream       xmlStream;
  SimInspiralTable     *currentSimEvent=NULL,*prevSimEvent=NULL;
  SimInspiralTable     *simEventList=NULL;
  MetadataTable	        simTable;


  /* Other global variables. */
  DetectorResponse detector;   /* the detector in question */
  REAL4TimeSeries output;      /* detector ADC output */
  INT4	  numinjects=0;
  REAL4	  Qfac = 0; /* Q factor for the "ringdown" */
  INT4	  smoothEnd = FALSE; /* do we include the "ringdown" */
  REAL4	  dcfactor = 1; /* calibration factor between darm and inj */

  /*******************************************************************
   * BEGIN PARSE ARGUMENTS					     *
   *******************************************************************/

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
      {"verbose",	      no_argument,  &verbose_flag, 1},
      {"source",              required_argument,  0,  'a'},
      {"actuation",           required_argument,  0,  'b'},
      {"summary",	      required_argument,  0,  'c'},
      {"ifo",		      required_argument,  0,  'd'},
      {"darm2inj",	      required_argument,  0,  'e'},
      {"length",	      required_argument,  0,  'f'},
      {"samplerate",	      required_argument,  0,  'g'},
      {"help",                no_argument,	  0,  'h'},
      {"smooth",	      required_argument,  0,  'i'},    
      {"flow",		      required_argument,  0,  'j'},
      {"fhigh",		      required_argument,  0,  'k'},
      {"debug-level",	      required_argument,  0,  'l'},
      {"user-tag",	      required_argument,  0,  'm'},
      {0, 0, 0, 0}
    };

    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long ( argc, argv, "a:b:c:d:f:g:h:j:k:l:p", 
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
        /* Parse source file option. */
        {
          sourcefile = optarg;
          ADD_PROCESS_PARAM( "string", "%s", sourcefile );
        }
        break;

      case 'b':
        /* Parse actuation file option. */
        {
          actfile = optarg;
          ADD_PROCESS_PARAM( "string", "%s", actfile );
        }
        break;

      case 'c':
        /* Parse summfile option */
        {
          summfile = optarg;
          xmloutput = TRUE;
          ADD_PROCESS_PARAM( "string", "%s", summfile );
        }
        break;

      case 'd':
        /* IFO */
        {
          LALSnprintf( ifo, sizeof(ifo), optarg);
          ADD_PROCESS_PARAM( "ifo", "%s", ifo );
        }
        break;

      case 'e':
        /* Calibration between darm and injection */
        {
          dcfactor = atof( optarg );
          ADD_PROCESS_PARAM( "float", "%e", dcfactor );
        }
        break;

      case 'f':
        /* Specify length of data segment */
        {
          length = atoi( optarg );
          ADD_PROCESS_PARAM( "int", "%d", length );
        }
        break;

      case 'g':
        /* Specify the sampling rate of the channel */
        {
          freq = atoi( optarg );
          ADD_PROCESS_PARAM( "int", "%d", freq );
        }
        break;

      case 'h':
        /* print help */
        {
          LALPrintError( USAGE, *argv );
          return INSPAWGFILEC_EARG;
        }
        break;

      case 'i':
        /* end the injection with a "ringdown", Q factor = Q */
        {
          Qfac = atof( optarg );
          ADD_PROCESS_PARAM( "float", "%e", Qfac );
          smoothEnd = TRUE;
        }
        break;	

      case 'j':
        /* Parse start frequency */
        {
          fstart = atof( optarg );
          ADD_PROCESS_PARAM( "float", "%e", fstart );
        }
        break;

      case 'k':
        /* Parse stop frequency */
        {
          fstopset = TRUE;
          fstop = atof( optarg );
          ADD_PROCESS_PARAM( "float", "%e", fstop );
        }
        break;

      case 'l':
        /* Parse debug level option. */
        {
          set_debug_level( optarg );
          ADD_PROCESS_PARAM( "string", "%e", optarg );
        }
        break;

      case 'm':
        /* Parse user-tag. */
        {
          tag = optarg;
          ADD_PROCESS_PARAM( "string", "%s", tag );
        }
        break;


      default:
        {
          fprintf( stderr, "unknown error while parsing options\n" );	 
          return INSPAWGFILEC_EARG;
        }

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


  /*******************************************************************
   * initialize things
   *******************************************************************/
  memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );

  /*******************************************************************
   * SETUP                                                           *
   *******************************************************************/

  /* Set up output and detector structures. */
  output.data = NULL;
  detector.transfer = (COMPLEX8FrequencySeries *)
    LALMalloc( sizeof(COMPLEX8FrequencySeries) );
  if ( !(detector.transfer) ) 
  {
    ERROR( INSPAWGFILEC_EMEM, INSPAWGFILEC_MSGEMEM, 0 );
    return INSPAWGFILEC_EMEM;
  }
  detector.transfer->data = NULL;
  detector.site = NULL;

  /* Set up units. */
  {
    RAT4 negOne = { -1, 0 };
    LALUnit unit;
    LALUnitPair pair;
    output.sampleUnits = lalADCCountUnit;
    pair.unitOne = &lalADCCountUnit;
    pair.unitTwo = &lalStrainUnit;
    LAL_CALL( LALUnitRaise( &stat, &unit, pair.unitTwo, &negOne ), &stat );
    pair.unitTwo = &unit;
    LAL_CALL( LALUnitMultiply( &stat, &(detector.transfer->sampleUnits),
          &pair ), &stat );
  }

  /* Read actuation function -- which is stored as a 3 column vector.
     the first column is the frequency, second is amplitude, third is
     phase. */

  if ( actfile ) 
  {
    REAL4VectorSequence *actuation = NULL; /* actuation as vector sequence */
    COMPLEX8Vector *response = NULL;  /* response as complex vector */
    COMPLEX8Vector *unity = NULL;     /* complex vector of 1's */

    if ( ( fp = fopen( actfile, "r" ) ) == NULL )
    {
      ERROR( INSPAWGFILEC_EFILE, INSPAWGFILEC_MSGEFILE, actfile );
      return INSPAWGFILEC_EFILE;
    }

    LAL_CALL( LALSReadVectorSequence( &stat, &actuation, fp ), &stat );
    fclose( fp );
    if ( actuation->vectorLength != 3 ) 
    {
      ERROR( INSPAWGFILEC_EINPUT, INSPAWGFILEC_MSGEINPUT, actfile );
      return INSPAWGFILEC_EINPUT;
    }

    /* extract f0 and deltaF from the actuation function */
    detector.transfer->f0 = actuation->data[0];
    detector.transfer->deltaF = actuation->data[3] - actuation->data[0];

    /* Read and convert body to a COMPLEX8Vector. */

    LAL_CALL( LALCCreateVector( &stat, &response, actuation->length ), &stat );
    LAL_CALL( LALCCreateVector( &stat, &unity, response->length ), &stat );
    for ( i = 0; i < response->length; i++ ) 
    {
      response->data[i].re = actuation->data[3*i-2] * 
        cos( actuation->data[3*i-1] ) / dcfactor;
      response->data[i].im = actuation->data[3*i-2] * 
        sin( actuation->data[3*i-1] ) / dcfactor;
      unity->data[i].re = 1.0;
      unity->data[i].im = 0.0;
    }

    LAL_CALL( LALSDestroyVectorSequence( &stat, &actuation ), &stat );

    /* convert response function to transfer function */    
    LAL_CALL( LALCCreateVector( &stat, &( detector.transfer->data ),
          response->length ), &stat );
    LAL_CALL( LALCCVectorDivide( &stat, detector.transfer->data, unity,
          response ), &stat );
    LAL_CALL( LALCDestroyVector( &stat, &response ), &stat );
    LAL_CALL( LALCDestroyVector( &stat, &unity ), &stat );
  }

  /* No response file, so generate a unit response. */
  else 
  {
    epoch = EPOCH;
    LALINT8toGPS( &stat, &( detector.transfer->epoch ), &epoch );
    detector.transfer->f0 = 0.0;
    detector.transfer->deltaF = 1.5*fstop;
    LAL_CALL( LALCCreateVector( &stat, &( detector.transfer->data ), 2 ),
        &stat );
    detector.transfer->data->data[0].re = 1.0;
    detector.transfer->data->data[1].re = 1.0;
    detector.transfer->data->data[0].im = 0.0;
    detector.transfer->data->data[1].im = 0.0;
  }




  /*******************************************************************
   * INJECTION                                                       *
   *******************************************************************/

  /* Open sourcefile. */
  if ( sourcefile ) 
  {
    if ( ( fp = fopen( sourcefile, "r" ) ) == NULL )
    {
      ERROR( INSPAWGFILEC_EFILE, INSPAWGFILEC_MSGEFILE,
          sourcefile );
      return INSPAWGFILEC_EFILE;
    }
  }

  /* For each line in the sourcefile... */
  while ( ok ) 
  {
    PPNParamStruc ppnParams;       /* wave generation parameters */
    REAL4 m1, m2, dist, inc, phic; /* unconverted parameters */
    CoherentGW waveform,*wf;       /* amplitude and phase structure */
    REAL4TimeSeries signal;        /* GW signal */
    REAL8 time;                    /* length of GW signal */
    CHAR timeCode;                 /* code for signal time alignment */
    CHAR message[MSGLEN];          /* warning/info messages */

    /* Read and convert input line. */
    if ( sourcefile ) 
    {
      ok &= ( fscanf( fp, "%c %lli %f %f %f %f %f\n", &timeCode,
            &epoch, &m1, &m2, &dist, &inc, &phic ) == 7 );
      if ( ok )
      {
        fprintf(stderr, "%c %lli %f %f %f %f %f\n", timeCode,
            epoch, m1, m2, dist, inc, phic );  fflush(stderr);
        ppnParams.mTot = m1 + m2;
        ppnParams.eta = m1*m2/( ppnParams.mTot*ppnParams.mTot );
        ppnParams.d = dist*LAL_PC_SI*1000000.0;
        ppnParams.inc = inc*LAL_PI/180.0;
        ppnParams.phi = phic*LAL_PI/180.0;
      }
    } 
    else 
    {
      timeCode = 'i';
      ppnParams.mTot = M1 + M2;
      ppnParams.eta = M1*M2/( ppnParams.mTot*ppnParams.mTot );
      ppnParams.d = DIST*LAL_PC_SI*1000000.0;
      ppnParams.inc = INC;
      ppnParams.phi = PHIC;
      epoch = EPOCH;
    }

    if ( ok ) 
    {
      /* Set up other parameter structures. */
      ppnParams.epoch.gpsSeconds = ppnParams.epoch.gpsNanoSeconds = 0;
      ppnParams.position.latitude = ppnParams.position.longitude = 0.0;
      ppnParams.position.system = COORDINATESYSTEM_EQUATORIAL;
      ppnParams.psi = 0.0;
      ppnParams.fStartIn = fstart;
      /* set end frequency to fstop, if given, otherwise continue to ISCO */
      if ( fstopset == TRUE )
      {
        ppnParams.fStopIn = fstop;
      }
      else
      {
        ppnParams.fStopIn =  -1.0 / 
          (6.0 * sqrt(6.0) * LAL_PI * ppnParams.mTot * LAL_MTSUN_SI);
      }
      ppnParams.lengthIn = 0;
      ppnParams.ppn = NULL;
      ppnParams.deltaT = DELTAT;
      memset( &waveform, 0, sizeof(CoherentGW) );

      /* Generate empty output file to which injection will be added. */
      output.epoch.gpsSeconds = 0;
      output.epoch.gpsNanoSeconds = 0;
      dt = (DT*FREQ)/freq;
      npt = length*freq;
      output.deltaT = dt; 
      LAL_CALL( LALSCreateVector( &stat, &( output.data ), npt ), &stat );
      memset( output.data->data, 0, npt*sizeof(REAL4) );

      /* Generate waveform at zero epoch. */
      LAL_CALL( LALGeneratePPNInspiral( &stat, &waveform, &ppnParams ),
          &stat );
      LALSnprintf( message, MSGLEN, "%d: %s", ppnParams.termCode,
          ppnParams.termDescription );
      INFO( message );
      if ( lalDebugLevel & LALINFO )
      {
        LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"
            "\tWaveform terminated at %e Hz\n", 
            *argv, __FILE__, __LINE__,
            rcsid, ppnParams.fStop );
      }
      if ( ppnParams.dfdt > 2.0 ) 
      {
        LALPrintError( "Waveform sampling interval is too large:\n"
            "\tmaximum df*dt = %f", ppnParams.dfdt );
      }

      /* add the smooth ending */
      if( smoothEnd == TRUE )
      {
        wf = &waveform;
        LAL_CALL( LALGenerateInspiralSmooth( &stat, &wf, &ppnParams, 
              &Qfac ), &stat );
        waveform = *wf;
      }

      /* Compute epoch for waveform. */
      time = waveform.a->data->length*DELTAT;
      if ( timeCode == 'f' )
        epoch -= (INT8)( 1000000000.0*time );
      else if ( timeCode == 'c' )
        epoch -= (INT8)( 1000000000.0*ppnParams.tc );
      LALINT8toGPS( &stat, &( waveform.a->epoch ), &epoch );
      waveform.f->epoch = waveform.phi->epoch = waveform.a->epoch;

      /* Generate and inject signal. */
      signal.epoch = waveform.a->epoch;
      signal.epoch.gpsSeconds -= 1;
      signal.deltaT = output.deltaT/4.0;
      signal.f0 = 0;
      signal.data = NULL;
      time = ( time + 2.0 )/signal.deltaT;
      LAL_CALL( LALSCreateVector( &stat, &( signal.data ), (UINT4)time ),
          &stat );
      LAL_CALL( LALSimulateCoherentGW( &stat, &signal, &waveform,
            &detector ), &stat );
      LAL_CALL( LALSSInjectTimeSeries( &stat, &output, &signal ),
          &stat );
      LAL_CALL( LALSDestroyVectorSequence( &stat, &( waveform.a->data ) ),
          &stat );
      LAL_CALL( LALSDestroyVector( &stat, &( waveform.f->data ) ), &stat );
      LAL_CALL( LALDDestroyVector( &stat, &( waveform.phi->data ) ), &stat );
      LALFree( waveform.a );
      LALFree( waveform.f );
      LALFree( waveform.phi );
      LAL_CALL( LALSDestroyVector( &stat, &( signal.data ) ), &stat );

      /* add 1 to number of injections */
      numinjects += 1;

      /* allocate memory for storing the injected events */
      if (prevSimEvent == NULL)
      {
        currentSimEvent = simEventList = (SimInspiralTable *)            
          LALCalloc( 1 , sizeof(SimInspiralTable) );
      }
      else
      {
        currentSimEvent = (SimInspiralTable *)
          LALCalloc( 1 , sizeof(SimInspiralTable) );
        prevSimEvent->next = currentSimEvent;
      }
      /* add information about current event */
      LALSnprintf(currentSimEvent->waveform, sizeof(currentSimEvent->waveform),
          "GeneratePPNtwoPN");

      LALFloatToGPS( &stat, &inj_length, &(ppnParams.tc));
      currentSimEvent->geocent_end_time.gpsSeconds = 
        epoch + inj_length.gpsSeconds;
      currentSimEvent->geocent_end_time.gpsNanoSeconds = 
        inj_length.gpsNanoSeconds;
      currentSimEvent->h_end_time = currentSimEvent->l_end_time = 
        currentSimEvent->geocent_end_time;
      currentSimEvent->end_time_gmst = 0;
      LALSnprintf(currentSimEvent->source, sizeof(currentSimEvent->waveform),
          "HW");
      currentSimEvent->mass1 = m1;
      currentSimEvent->mass2 = m2;
      currentSimEvent->eta = ppnParams.eta;
      currentSimEvent->distance = dist; 
      currentSimEvent->longitude = 0;
      currentSimEvent->latitude = 0;
      currentSimEvent->inclination = inc;
      currentSimEvent->coa_phase = phic;
      currentSimEvent->eff_dist_h = dist;
      currentSimEvent->eff_dist_l = dist;

      prevSimEvent = currentSimEvent;

      /* Print output file. */
      if ( tag ) 
      {
        if ( !strcmp(ifo,"") )
        {
          LALSnprintf( fname, sizeof(fname), "%s_inspiral_%d.out", tag, 
              numinjects);
        }
        else
        {
          LALSnprintf( fname, sizeof(fname), "%s_inspiral_%d_%s.out",
              tag, numinjects, ifo);
        }
      }	
      else
      {
        if ( !strcmp(ifo,"") )
        {
          LALSnprintf( fname, sizeof(fname), "inspiral_%d.out", numinjects);
        }
        else
        {
          LALSnprintf( fname, sizeof(fname), "inspiral_%d_%s.out",
              numinjects, ifo);
        }
      }

      if ( ( fq = fopen( fname, "w" ) ) == NULL ) 
      {
        ERROR( INSPAWGFILEC_EFILE, INSPAWGFILEC_MSGEFILE,fname );
        return INSPAWGFILEC_EFILE;
      }
      for ( i = 0; i < output.data->length; i++ )
        fprintf( fq, "%e\n", output.data->data[i] );
      fclose( fq ); 

      /* destroy current injection */
      LAL_CALL( LALSDestroyVector( &stat, &( output.data ) ), &stat );    
    }

    /* If there is no source file, inject only one source. */
    if ( !sourcefile )
    {
      ok = 0;
    }

  }

  /* Input file is exhausted (or has a badly-formatted line ). */
  if ( sourcefile )
  {  
    fclose( fp );
  }
  /*****************************************************************
   * open summary xml file
   *****************************************************************/
  if ( summfile )
  {
    LAL_CALL( LALOpenLIGOLwXMLFile(&stat, &xmlStream, summfile), &stat);


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

    /* Write the list of injections to the sim table */
    if ( simEventList )
    {
      LAL_CALL( LALBeginLIGOLwXMLTable (&stat, &xmlStream, sim_inspiral_table),
          &stat);
      simTable.simInspiralTable = simEventList;
      LAL_CALL( LALWriteLIGOLwXMLTable (&stat, &xmlStream, simTable, 
            sim_inspiral_table), &stat);
      LAL_CALL( LALEndLIGOLwXMLTable (&stat, &xmlStream), &stat);

      /* free the temporary memory containing the events */
      while (simEventList)
      {
        currentSimEvent = simEventList;
        simEventList = simEventList->next;
        LALFree( currentSimEvent );
      }
    }

    /* close the summary file */
    LAL_CALL( LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);
  }

  /*******************************************************************
   * CLEANUP                                                         *
   *******************************************************************/


  /* Destroy remaining memory. */
  LAL_CALL( LALCDestroyVector( &stat, &( detector.transfer->data ) ), &stat );
  LALFree( detector.transfer );

  /* Done! */
  LALCheckMemoryLeaks();
  INFO( INSPAWGFILEC_MSGENORM );
  return INSPAWGFILEC_ENORM;
}


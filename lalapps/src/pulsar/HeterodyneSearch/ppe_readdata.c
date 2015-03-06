/*
*  Copyright (C) 2014 Matthew Pitkin
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


/******************************************************************************/
/*                        DATA READING FUNCTIONS                              */
/******************************************************************************/

#include "ppe_readdata.h"

/**
 * \brief Reads in the input gravitational wave data files, or creates fake data sets.
 *
 * The function will check whether data files are being input of fake data is to be generated. If using real data the \c
 * detectors command line input must list the names of the detectors from which each of the data sets comes, with names
 * separated by commas - allowed detector names are H1, H2, L1, V1, G1 or T1, for LIGO Hanford Observatory 1, LIGO
 * Hanford Observatory 2, LIGO Livingston Observatory, Virgo, GEO600/HF and TAMA300 respectively. If requiring fake data
 * to be generated then the \c fake-data command line argument must list the detectors for which fake data is required
 * (again separated by commas) - these can be the same names as above, although if an 'A' is appended to the LIGO of
 * Virgo detector name then the Advanced detector is assumed (for use if generating data from designed sensitivity
 * curves).
 *
 * If generating fake data the noise floor can either be created using models of the detector design sensitivities (if
 * just \c fake-data is set), or with noise levels defined using the \c fake-psd command line argument. If using \c
 * fake-psd it should list the signle-sided power spectral density required for each detector listed in \c fake-data
 * (again separated by commas). If using the design sensitivity models the \c par-file argument will be used to find the
 * noise at the correct frequency, which is here assumed to be twice the rotation frequency. The start time (in GPS
 * seconds), data length (in seconds) and sample interval (in seconds) can be specified for each fake data set by
 * supplying comma separated lists to the \c fake-starts, \c fake-lengths and \c fake-dt command line arguments. By
 * default these values are GPS 900000000 (13th July 2008 at 15:59:46 UTC), 86400 seconds (1 solar day) and 60 seconds
 * (1/60 Hz sample rate) respectively. The fake data is drawn from a normal distribution using \c XLALNormalDeviates and
 * scaled with the appropriate PSD.
 *
 * The number of data files required to be read in, or number of fake data sets generated depends on the pulsar model
 * type, which is specified by the number of frequency harmonics given by the \c harmonics command line argument. This
 * should be a list of comma separated values giving the frequency of the signal harmonics to be included. E.g.
 * <ol>
 * <li>For a model with \f$l=m=2\f$ (i.e. a triaxial star with a signal defined in e.g. \cite DupuisWoan2005), which
 * purely emits at twice the rotation frequency, the \c harmonics input would just be \c 2. This requires that for each
 * specified detector there is <b>one</b> input file containing data heterodyned at twice the rotation frequency of the
 * pulsar.</li>
 * <li>For a model including the two harmonics \f$l=2\f$, \f$m=1,2\f$ (see e.g. \cite Jones2010), which produces
 * emission at both the rotation frequency <i>and</i> twice the rotation frequency, the \c harmonics input would be \c
 * 1,2. This requires that for each specified detector there are two input files containing data heterodyned and the
 * rotation frequency <i>and</i> twice the rotation frequency (these must be in the same order as the harmonics
 * are listed).</li>
 * </ol>
 * The default model, if none is specified, is the triaxial source model with emission at just twice the rotation
 * frequency. At the moment only these two models above can be used, although this could be expanded in the future.
 *
 * If creating fake data for a specific model then the number of PSDs, start time, lengths and sample intervals
 * specified must be equivalent to the number of input files that would have been required e.g. if using the pinned
 * superfluid model and requiring data for H1 and L1 then four fake PSDs would be required (the first pair at the
 * pulsars rotation frequency and twice that in H1, and the seconds pair at the pulsars rotation frequency and twice
 * that in L1). These most be specified in the same order as the detectors.
 *
 * Any new models added can require and arbitrary number of inputs for a given detector, however the heterodyned
 * frequencies of each input must be hardcoded into the function.
 *
 * If using real data the files must be specified in the \c input-files command line argument - these should be comma
 * separated for multiple files and be in the same order at the associated detector from which they came given by the
 * \c detectors command.
 *
 * The function also checks that valid Earth and Sun ephemeris files (from the lalpulsar suite) are set with the \c
 * ephem-earth and \c ephem-sun arguments, and that a valid output file for the nested samples is set via the \c outfile
 * argument.
 *
 * The function will by default also call \c chop_n_merge() for each data set, which will split the data into chunks
 * during which it can be considered Gaussian and stationary. The command line arguments \c chunk-min and \c chunk-max
 * can be used to specify hardwire minimum and maximum lengths of chunks that are allowable. By default the maximum
 * chunk length is 0, which corresponds to no maximum value being set. If the \c --oldChunks flag is set then data will
 * be split as in the older version of the parameter estimation code, in which the chunk length is fixed, except for the
 * possibility of shorter segments at the end of science segments.
 *
 * The \c log_factorial array is also filled in with values of the log of the factorial of all numbers up to the maximum
 * length of the data. This is used in likelihood calculations.
 *
 * \param runState [in] A pointer to the LALInferenceRunState
 */
void read_pulsar_data( LALInferenceRunState *runState ){
  ProcessParamsTable *ppt = NULL, *ppt2 = NULL;
  ProcessParamsTable *commandLine = runState->commandLine;

  CHAR *detectors = NULL;
  CHAR *inputfile = NULL;

  CHAR *filestr = NULL;

  CHAR *efile = NULL, *sfile = NULL, *tfile = NULL;
  TimeCorrectionType ttype;

  CHAR *tempdets = NULL;
  CHAR *tempdet = NULL;

  REAL8 *fpsds = NULL;
  CHAR *fakestarts = NULL, *fakelengths = NULL, *fakedt = NULL;
  REAL8 *fstarts = NULL, *flengths = NULL, *fdt = NULL;

  CHAR dets[MAXDETS][256];
  INT4 numDets = 0, i = 0, numPsds = 0, numLengths = 0, numStarts = 0;
  INT4 numDt = 0, count = 0;
  UINT4 maxlen = 0;

  LALInferenceIFOData *ifodata = NULL;
  LALInferenceIFOData *prev = NULL;

  LALInferenceIFOModel *ifomodel = NULL;
  LALInferenceIFOModel *prevmodel = NULL;

  UINT4 seed = 0; /* seed for data generation */
  RandomParams *randomParams = NULL;

  CHAR *harmonics = NULL;
  REAL8Vector *modelFreqFactors = NULL;
  INT4 ml = 1;

  CHAR *parFile = NULL;
  BinaryPulsarParams pulsar;

  runState->data = NULL;

  /* Initialize the model, as it will hold IFO params and signal buffers */
  runState->model = XLALMalloc(sizeof(LALInferenceModel));

  /* timing values */
  struct timeval time1, time2;
  REAL8 tottime;

  if ( LALInferenceCheckVariable( runState->algorithmParams, "timefile" ) ){ gettimeofday(&time1, NULL); }

  /* check pulsar model required by getting the frequency harmonics */
  ppt = LALInferenceGetProcParamVal( commandLine, "--harmonics" );
  if ( ppt ) { harmonics = XLALStringDuplicate( ppt->value ); }
  else { harmonics = XLALStringDuplicate( "2" ); } /* default to using harmonic at twice the rotation rate */

  CHAR *tmpharms = NULL, *tmpharm = NULL, harmval[256];

  ml = count_csv( harmonics );
  modelFreqFactors = XLALCreateREAL8Vector( ml );

  tmpharms = XLALStringDuplicate( harmonics );
  for( i = 0; i < ml; i++ ){
    tmpharm = XLALStringToken( &tmpharms, "," , 0);
    XLALStringCopy( harmval, tmpharm, strlen(tmpharm)+1 );

    modelFreqFactors->data[i] = atof(harmval);
  }

  ppt = LALInferenceGetProcParamVal( runState->commandLine, "--par-file" );
  if( ppt == NULL ) { fprintf(stderr,"Must specify --par-file!\n"); exit(1); }
  parFile = XLALStringDuplicate( ppt->value );

  /* get the pulsar parameters to give a value of f */
  XLALReadTEMPOParFile( &pulsar, parFile );

  /* get the detectors - must */
  ppt = LALInferenceGetProcParamVal( commandLine, "--detectors" );
  ppt2 = LALInferenceGetProcParamVal( commandLine, "--fake-data" );
  if( ppt && !ppt2 ){
    detectors = XLALStringDuplicate( ppt->value );

    /* count the number of detectors from command line argument of comma separated vales and set their names */
    tempdets = XLALStringDuplicate( detectors );

    if ( (numDets = count_csv( detectors )) > MAXDETS ){
      fprintf(stderr, "Error... too many detectors specified. Increase MAXDETS to be greater than %d if necessary.\n", MAXDETS);
      exit(0);
    }

    for( i = 0; i < numDets; i++ ){
      tempdet = XLALStringToken( &tempdets, "," , 0);
      XLALStringCopy( dets[i], tempdet, strlen(tempdet)+1 );
    }
  }
  /*Get psd values for generating fake data.*/
  /*=========================================================================*/
  /*if using fake data and no detectors are specified*/
  else if( ppt2 && !ppt ){
    detectors = XLALStringDuplicate( ppt2->value );

    fpsds = XLALCalloc( MAXDETS * ml, sizeof(REAL8) );

    if ( (numDets = count_csv( detectors )) > MAXDETS ){
      fprintf(stderr, "Error... too many detectors specified. Increase MAXDETS to be greater than %d if necessary.\n", MAXDETS);
      exit(0);
    }
    /*------------------------------------------------------------------------*/
    /* get noise psds if specified */
    ppt = LALInferenceGetProcParamVal(commandLine,"--fake-psd");
    if( ppt ){
      CHAR *psds = NULL, *tmppsds = NULL, *tmppsd = NULL, psdval[256];

      psds = XLALStringDuplicate( ppt->value );

      tmppsds = XLALStringDuplicate( psds );
      tempdets = XLALStringDuplicate( detectors );

      /* count the number of PSDs (comma seperated values) to compare to number of detectors */
      if( (numPsds = count_csv( psds )) != ml*numDets ){
        fprintf(stderr, "Error... for %d harmonics the number of PSDs for fake data must be %d times the number of \
detectors specified (no. dets = %d)\n", ml, ml, numDets);
        exit(0);
      }

      for( i = 0; i < ml*numDets; i++ ){
        CHAR *tmpstr = NULL;

        tmppsd = XLALStringToken( &tmppsds, "," , 0);
        XLALStringCopy( psdval, tmppsd, strlen(tmppsd)+1 );
        fpsds[i] = atof(psdval);

        /* set detector */
        if ( i%ml == 0 ){
          tempdet = XLALStringToken( &tempdets, "," , 0);

          if( (tmpstr = strstr(tempdet, "A")) != NULL ){ /* have advanced */
            XLALStringCopy( dets[FACTOR(i,ml)], tmpstr+1, strlen(tmpstr)+1 );
          }
          else { XLALStringCopy( dets[FACTOR(i,ml)], tempdet, strlen(tempdet)+1 ); }
        }
      }
    }
    /*------------------------------------------------------------------------*/
    else{ /* get PSDs from model functions and set detectors */
      REAL8 pfreq = 0.;

      /* putting in pulsar frequency at f here */
      pfreq = pulsar.f0;

      tempdets = XLALStringDuplicate( detectors );

      for( i = 0; i < ml*numDets; i++ ){
        CHAR *tmpstr = NULL;
        REAL8 psdvalf = 0.;

        numPsds++;

        if( i%ml == 0 ) { tempdet = XLALStringToken( &tempdets, "," , 0); }

        if( (tmpstr = strstr(tempdet, "A")) != NULL ){ /* have Advanced */
          XLALStringCopy( dets[FACTOR(i,ml)], tmpstr+1, strlen(tmpstr)+1 );

          if( !strcmp(dets[FACTOR(i,ml)], "H1") || !strcmp(dets[FACTOR(i,ml)], "L1") || !strcmp(dets[FACTOR(i,ml)], "H2") ){ /* ALIGO */
            psdvalf = XLALSimNoisePSDaLIGOZeroDetHighPower( pfreq*modelFreqFactors->data[(INT4)(fmod(i,ml))] );
          }
          else if( !strcmp(dets[FACTOR(i,ml)], "V1") ){ /* AVirgo */
            psdvalf = XLALSimNoisePSDAdvVirgo( pfreq*modelFreqFactors->data[(INT4)(fmod(i,ml))] );
          }
          else{
            fprintf(stderr, "Error... trying to use Advanced detector that is not available!\n");
            exit(0);
          }
        }
        else{ /* initial detector */
          XLALStringCopy( dets[FACTOR(i,ml)], tempdet, strlen(tempdet)+1 );

          if( !strcmp(dets[FACTOR(i,ml)], "H1") || !strcmp(dets[FACTOR(i,ml)], "L1") || !strcmp(dets[FACTOR(i,ml)], "H2") ){ /* Initial LIGO */
            psdvalf = XLALSimNoisePSDiLIGOSRD( pfreq*modelFreqFactors->data[(INT4)(fmod(i,ml))] );

            /* divide H2 psds by 2 */
            if( !strcmp(dets[FACTOR(i,ml)], "H2") ) { psdvalf /= 2.; }
          }
          else if( !strcmp(dets[FACTOR(i,ml)], "V1") ){ /* Initial Virgo */
            psdvalf = XLALSimNoisePSDVirgo( pfreq*modelFreqFactors->data[(INT4)(fmod(i,ml))] );
          }
          else if( !strcmp(dets[FACTOR(i,ml)], "G1") ){ /* GEO 600 */
            psdvalf = XLALSimNoisePSDGEO( pfreq*modelFreqFactors->data[(INT4)(fmod(i,ml))] );
          }
          else if( !strcmp(dets[FACTOR(i,ml)], "T1") ){ /* TAMA300 */
            psdvalf = XLALSimNoisePSDTAMA( pfreq*modelFreqFactors->data[(INT4)(fmod(i,ml))] );
          }
          else{
            fprintf(stderr, "Error... trying to use detector that is not available!\n");
            exit(0);
          }
        }

        fpsds[i] = psdvalf;
      }
    }
    /*generate the fake data timestamps.*/
    /*====================================================================*/

    fstarts = XLALCalloc( MAXDETS*ml, sizeof(REAL8) );
    ppt = LALInferenceGetProcParamVal(commandLine,"--fake-starts");
    if( ppt ){
      CHAR *tmpstarts = NULL, *tmpstart = NULL, startval[256];
      fakestarts = XLALStringDuplicate( ppt->value );

      if( (numStarts = count_csv( fakestarts )) != numDets*ml ){
        fprintf(stderr, "Error... for %d harmonics the number of start times for fake data must be %d times the number \
of detectors specified (no. dets = %d)\n", ml, ml, numDets);
        exit(0);
      }

      tmpstarts = XLALStringDuplicate( fakestarts );
      for( i = 0; i < ml*numDets; i++ ){
        tmpstart = XLALStringToken( &tmpstarts, "," , 0);
        XLALStringCopy( startval, tmpstart, strlen(tmpstart)+1 );

        fstarts[i] = atof(startval);
      }
    }
    else{ /* set default GPS 900000000 - 13th July 2008 at 15:59:46 */
      for(i = 0; i < ml*numDets; i++ ){ fstarts[i] = 900000000.; }
    }

    flengths = XLALCalloc( MAXDETS*ml, sizeof(REAL8) );
    ppt = LALInferenceGetProcParamVal( commandLine, "--fake-lengths" );
    if( ppt ){
      CHAR *tmplengths = NULL, *tmplength = NULL, lengthval[256];
      fakelengths = XLALStringDuplicate( ppt->value );

      if( (numLengths = count_csv( fakelengths )) != numDets*ml ){
        fprintf(stderr, "Error... for %d harmonics the number of data lengths for fake data must be %d times the \
number of detectors specified (no. dets = %d)\n", ml, ml, numDets);
        exit(0);
      }

      tmplengths = XLALStringDuplicate( fakelengths );
      for( i = 0; i < ml*numDets; i++ ){
        tmplength = XLALStringToken( &tmplengths, "," , 0);
        XLALStringCopy( lengthval, tmplength, strlen(tmplength)+1 );
        flengths[i] = atof(lengthval);
      }
    }
    else{ /* set default (86400 seconds or 1 day) */
      for(i = 0; i < ml*numDets; i++ ) { flengths[i] = 86400.; }
    }

    fdt = XLALCalloc( MAXDETS*ml, sizeof(REAL8) );
    ppt = LALInferenceGetProcParamVal( commandLine, "--fake-dt" );
    if( ppt ){
      CHAR *tmpdts = NULL, *tmpdt = NULL, dtval[256];
      fakedt = XLALStringDuplicate( ppt->value );

      if( (numDt = count_csv( fakedt )) != ml*numDets ){
        fprintf(stderr, "Error... for %d harmonics the number of sample time steps for fake data must be %d times the \
number of detectors specified (no. dets =%d)\n", ml, ml, numDets);
        exit(0);
      }

      tmpdts = XLALStringDuplicate( fakedt );

      for( i = 0; i < ml*numDets; i++ ){
        tmpdt = XLALStringToken( &tmpdts, "," , 0);
        XLALStringCopy( dtval, tmpdt, strlen(tmpdt)+1 );
        fdt[i] = atof(dtval);
      }
    }
    else{ /* set default (60 sesonds) */
      for(i = 0; i < ml*numDets; i++) { fdt[i] = 60.; }
    }

  }
  /*psds set and timestamps set.*/
  /*====================================================================*/
  else{
    fprintf(stderr, "Error... --detectors OR --fake-data needs to be set.\n");
    fprintf(stderr, USAGE, commandLine->program);
    exit(0);
  }

  runState->model->ifo_loglikelihoods = XLALMalloc( sizeof(REAL8)*ml*numDets );
  runState->model->ifo_SNRs = XLALMalloc( sizeof(REAL8)*ml*numDets );

  UINT4 nstreams = ml*numDets;
  LALInferenceAddVariable( runState->algorithmParams, "numstreams", &nstreams, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED );

  ppt = LALInferenceGetProcParamVal( commandLine,"--input-files" );
  if( ppt ) { inputfile = XLALStringDuplicate( ppt->value ); }

  if ( inputfile == NULL && !ppt2 ){
    fprintf(stderr, "Error... an input file or fake data needs to be set.\n");
    fprintf(stderr, USAGE, commandLine->program);
    exit(0);
  }

  /* get the output directory */
  ppt = LALInferenceGetProcParamVal( commandLine, "--outfile" );
  if( !ppt ){
    fprintf(stderr, "Error... --outfile needs to be set.\n");
    fprintf(stderr, USAGE, commandLine->program);
    exit(0);
  }

  /* count the number of input files (by counting commas) and check it's equal to twice the number of detectors */
  if ( !ppt2 ){ /* if using real data */
    count = count_csv( inputfile );

    if ( count != ml*numDets ){
      fprintf(stderr, "Error... for %d harmonics the number of input files given must be %d times the number of \
detectors specified (no. dets =%d)\n", ml, ml, numDets);
      exit(0);
    }
  }

  /* set random number generator in case when that fake data is used */
  ppt = LALInferenceGetProcParamVal( commandLine, "--randomseed" );
  if ( ppt != NULL ) { seed = atoi( ppt->value ); }
  else { seed = 0; } /* will be set from system clock */

  /* reset filestr if using real data (i.e. not fake) */
  if ( !ppt2 ) { filestr = XLALStringDuplicate( inputfile ); }

  /* set flags for whether noise sigma is input or required to be calculated */
  INT4 inputsigma = 0, computesigma = 0, gaussianLike = 0;

  if ( LALInferenceGetProcParamVal( commandLine, "--gaussian-like" ) ) { gaussianLike = 1; }

  /* set reheterodyne frequency (0 if no reheterodyne is required) */
  REAL8 rehetfreq = 0;
  if ( LALInferenceGetProcParamVal( commandLine, "--reheterodyne" ) ) {
    if ( *LALInferenceGetProcParamVal( commandLine, "--reheterodyne")->value == '\0'){
    fprintf(stderr, "Error... --reheterodyne needs frequency as argument.\n");
    fprintf(stderr, "Provide argument or remove flag\n.");
    exit(0);
    }
    else {
    rehetfreq = atof( LALInferenceGetProcParamVal( commandLine, "--reheterodyne" )->value );
    }
  }

  /* read in data, needs to read in two sets of data for each ifo for pinsf model */
  for( i = 0, prev=NULL, prevmodel=NULL ; i < ml*numDets ; i++, prev=ifodata, prevmodel=ifomodel ){
    CHAR *datafile = NULL;
    REAL8 times = 0;
    LIGOTimeGPS gpstime;
    REAL8 dataValsRe = 0., dataValsIm = 0., sigmaVals = 0.;
    REAL8Vector *temptimes = NULL;
    INT4 j = 0, k = 0, datalength = 0;
    ProcessParamsTable *ppte = NULL, *ppts = NULL, *pptt = NULL;

    CHAR *filebuf = NULL;

    count = 0;

    /* initialise random number generator */
    /* Moved into det loop so same random seed can be used with */
    /* different detector combos and still get same noise realisation */
    randomParams = XLALCreateRandomParams( seed+i );

    ifodata = XLALCalloc( 1, sizeof(LALInferenceIFOData) );
    ifodata->likeli_counter = 0;
    ifodata->templa_counter = 0;
    ifodata->next = NULL;

    ifomodel = XLALMalloc(sizeof(LALInferenceIFOModel));
    ifomodel->params = XLALCalloc(1, sizeof(LALInferenceVariables) );
    ifomodel->next = NULL;

    /* add frequency factors variable */
    LALInferenceAddVariable( ifomodel->params, "freqfactors", &modelFreqFactors, LALINFERENCE_REAL8Vector_t,
                             LALINFERENCE_PARAM_FIXED );

    /* check if using non-GR model */
    if ( LALInferenceGetProcParamVal( commandLine, "--nonGR" ) ){
      UINT4 nonGR = 1;
      LALInferenceAddVariable( ifomodel->params, "nonGR", &nonGR, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED );
      /* check if nonGR argument has an associated input variable */
      if ( *LALInferenceGetProcParamVal( commandLine, "--nonGR")->value != '\0'){
        CHAR *nonGRmodel = NULL;
        nonGRmodel =  XLALStringDuplicate( LALInferenceGetProcParamVal( commandLine, "--nonGR" )->value);
        LALInferenceAddVariable( ifomodel->params, "nonGRmodel", &nonGRmodel, LALINFERENCE_string_t, LALINFERENCE_PARAM_FIXED );
      }
    }

    if( i == 0 ) {
        runState->data = ifodata;
        runState->model->ifo = ifomodel;
    }
    if( i > 0 ) {
        prev->next = ifodata;
        prevmodel->next = ifomodel;
    }

    /* set detector */
    ifodata->detector = XLALGetSiteInfo( dets[FACTOR(i,ml)] );
    ifomodel->detector = XLALGetSiteInfo( dets[FACTOR(i,ml)] );
    snprintf(ifodata->name, sizeof(char)*DETNAMELEN, "%s", dets[FACTOR(i,ml)]);

    /* set dummy initial time */
    gpstime.gpsSeconds = 0;
    gpstime.gpsNanoSeconds = 0;

    /* allocate time domain data - will be dynamically allocated as data read*/
    ifodata->compTimeData = NULL;
    ifodata->compTimeData = XLALCreateCOMPLEX16TimeSeries( "", &gpstime, 0., 1., &lalSecondUnit, 1 );

    /* always allocate memory for variances (as these will always be needed in calculating SNRs) */
    ifodata->varTimeData = NULL;
    ifodata->varTimeData = XLALCreateREAL8TimeSeries( "", &gpstime, 0., 1., &lalSecondUnit, 1 );

    /* allocate time domain model */
    ifomodel->compTimeSignal = NULL;
    ifomodel->compTimeSignal = XLALCreateCOMPLEX16TimeSeries( "", &gpstime, 0., 1., &lalSecondUnit, 1 );

    /*============================ GET DATA ==================================*/
    /* get i'th filename from the comma separated list */
    if ( !ppt2 ){ /* if using real data read in from the file */
      datafile = XLALStringToken(&filestr, ",", 0);

      j=0;

      /* read in data */
      temptimes = XLALCreateREAL8Vector( 1 );

      /* read in all the data then ignore lines starting with # or % */
      filebuf = XLALFileLoad( datafile );

      /* separate data into lines */
      TokenList *tlist = NULL;
      if ( XLALCreateTokenList( &tlist, filebuf, "\n" ) != XLAL_SUCCESS ){
        fprintf(stderr, "Error... could not convert data into separate lines.\n");
        exit(3);
      }

      INT4 nvals = 0; /* number of values in a line */

      for ( k = 0; k < (INT4)tlist->nTokens; k++ ){
        /* search for a comment character in the string */
        if ( strchr(tlist->tokens[k], '#') || strchr(tlist->tokens[k], '%') ){ continue; }
        else{ /* read in data from string */
          if ( j == 0 ){
            /* check the number of values in the line by counting the number of value separated by whitespace  */
            TokenList *tline = NULL;
            XLALCreateTokenList( &tline, tlist->tokens[k], " \t" );
            nvals = (INT4)tline->nTokens;
            XLALDestroyTokenList( tline );

            /* set whether we need to compute variance, or whether the standard deviation is in the input file */
            if ( nvals == 3 && gaussianLike ){ computesigma = 1; }
            if ( nvals == 4 ){ inputsigma = 1; }
          }

          if ( nvals == 3 ){ /* no sigma value is in the input file */
            int rc = sscanf( tlist->tokens[k], "%lf%lf%lf", &times, &dataValsRe, &dataValsIm );
            if ( rc != nvals ){ continue; } /* ignore the line */
          }
          else if( nvals == 4 ){ /* sigma is in the input file */
            int rc = sscanf( tlist->tokens[k], "%lf%lf%lf%lf", &times, &dataValsRe, &dataValsIm, &sigmaVals );
            if ( rc != nvals ){ continue; } /* ignore the line */
          }
          else{
            fprintf(stderr, "Error... unrecognised number of values in first line of data file %s.\n", datafile);
            exit(3);
          }
        }
        j++;

        /* dynamically allocate more memory */
        ifodata->compTimeData = XLALResizeCOMPLEX16TimeSeries( ifodata->compTimeData, 0, j );
        ifomodel->compTimeSignal = XLALResizeCOMPLEX16TimeSeries( ifomodel->compTimeSignal, 0, j );
        ifodata->varTimeData = XLALResizeREAL8TimeSeries( ifodata->varTimeData, 0, j );

        temptimes = XLALResizeREAL8Vector( temptimes, j );

        // Note: j-1 because we added to j above (553)
        temptimes->data[j-1] = times;

        /* reheterodyne data if required */
        if ( rehetfreq != 0 ) {
          /* create template */
          REAL8 deltaPhase = 2. * LAL_PI * rehetfreq * times;
          COMPLEX16 dataTemp = ( dataValsRe*cos(-deltaPhase) - dataValsIm*sin(-deltaPhase) ) +
            I * ( dataValsRe*sin(-deltaPhase) + dataValsIm*cos(-deltaPhase) );
          ifodata->compTimeData->data->data[j-1] = dataTemp;
        }
        else {
          ifodata->compTimeData->data->data[j-1] = dataValsRe + I*dataValsIm;
        }

        if ( inputsigma ){ ifodata->varTimeData->data->data[j-1] = SQUARE( sigmaVals ); }
      }

      if ( j == 0 ){
        fprintf(stderr, "Error... nothing read in from data file %s.\n", datafile);
        exit(3);
      }

      XLALDestroyTokenList( tlist );
      XLALFree( filebuf );

      datalength = j;

      /* allocate data time stamps */
      ifomodel->times = NULL;
      ifomodel->times = XLALCreateTimestampVector( datalength );

      UINT4 epochint = 0; /* index of the earliest time in the data */

      /* fill in time stamps as LIGO Time GPS Vector */
      REAL8 sampledt = INFINITY; /* sample interval */
      for ( k = 0; k < datalength; k++ ) {
        XLALGPSSetREAL8( &ifomodel->times->data[k], temptimes->data[k] );

        if ( k > 0 ){
          /* get sample interval from the minimum time difference in the data */
          if ( temptimes->data[k] - temptimes->data[k-1] < sampledt ) {
            sampledt = temptimes->data[k] - temptimes->data[k-1];
          }
        }

        if ( temptimes->data[k] < temptimes->data[epochint] ){ epochint = k; }
      }

      ifodata->compTimeData->epoch = ifomodel->times->data[epochint];
      ifomodel->compTimeSignal->epoch = ifomodel->times->data[epochint];
      ifodata->varTimeData->epoch = ifomodel->times->data[epochint];

      /* check whether to randomise the data by shuffling the time stamps (this will preserve the order of
       * the data for working out stationary chunk, but randomise the signal) */
      if ( LALInferenceGetProcParamVal( commandLine, "--randomise" ) ){
        gsl_ran_shuffle( runState->GSLrandom, &ifomodel->times->data[0], (size_t)datalength, sizeof(LIGOTimeGPS) );
      }

      /* add data sample interval */
      ppt = LALInferenceGetProcParamVal( commandLine, "--sample-interval" );
      if( ppt ){ sampledt = atof( ppt->value ); }
      LALInferenceAddVariable( ifomodel->params, "dt", &sampledt, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );

      XLALDestroyREAL8Vector( temptimes );
    }
    else{ /* set up fake data */
      /* if a Gaussian likelihood is required compute sigma from the data (to mimic real life) */
      if ( gaussianLike ) { computesigma = 1; }

      datalength = flengths[i] / fdt[i];

      /* temporary real and imaginary data vectors */
      REAL4Vector *realdata = NULL;
      REAL4Vector *imagdata = NULL;

      REAL8 psdscale = 0.;

      /* allocate data time stamps */
      ifomodel->times = NULL;
      ifomodel->times = XLALCreateTimestampVector( (UINT4)datalength );

      /* add data sample interval */
      LALInferenceAddVariable( ifomodel->params, "dt", &fdt[0], LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );

      /* resize the data and model times series */
      ifodata->compTimeData = XLALResizeCOMPLEX16TimeSeries( ifodata->compTimeData, 0, datalength );
      ifomodel->compTimeSignal = XLALResizeCOMPLEX16TimeSeries( ifomodel->compTimeSignal, 0, datalength );
      ifodata->varTimeData = XLALResizeREAL8TimeSeries( ifodata->varTimeData, 0, datalength );

      /* create data drawn from normal distribution with zero mean and unit variance */
      realdata = XLALCreateREAL4Vector( datalength );
      imagdata = XLALCreateREAL4Vector( datalength );

      XLALNormalDeviates( realdata, randomParams );
      XLALNormalDeviates( imagdata, randomParams );

      /* converts single sided psd into double sided psd, and then into a time domain noise standard deviation */
      psdscale = sqrt( ( fpsds[i] / 2.) / ( 2. * fdt[i] ) );

      /* create time stamps and scale data with the PSD */
      for( k = 0; k < datalength; k++ ){
        /* set time stamp */
        XLALGPSSetREAL8( &ifomodel->times->data[k], fstarts[i] + fdt[i] * (REAL8)k );

        ifodata->compTimeData->data->data[k] = (REAL8)realdata->data[k] * psdscale + I * (REAL8)imagdata->data[k] * psdscale;
      }

      ifodata->compTimeData->epoch = ifomodel->times->data[0];
      ifomodel->compTimeSignal->epoch = ifomodel->times->data[0];
      ifodata->varTimeData->epoch = ifomodel->times->data[0];

      XLALDestroyREAL4Vector( realdata );
      XLALDestroyREAL4Vector( imagdata );
    }

    /* set ephemeris data */
    ifomodel->ephem = XLALMalloc( sizeof(EphemerisData) );

    /* get ephemeris files */
    ppte = LALInferenceGetProcParamVal( commandLine, "--ephem-earth" );
    ppts = LALInferenceGetProcParamVal( commandLine, "--ephem-sun" );
    pptt = LALInferenceGetProcParamVal( commandLine, "--ephem-timecorr" );
    if( ppte && ppts ){
      efile = XLALStringDuplicate( ppte->value );
      sfile = XLALStringDuplicate( ppts->value );

      if ( pptt ){
        tfile = XLALStringDuplicate( pptt->value );

        if ( pulsar.units != NULL ){
          if( !strcmp(pulsar.units, "TDB") ) { ttype = TIMECORRECTION_TDB; }
          else { ttype = TIMECORRECTION_TCB; } /* default to TCB otherwise */
        }
        else { ttype = TIMECORRECTION_TCB; }
      }
      else{
        tfile = NULL;
        ttype = TIMECORRECTION_ORIGINAL;
      }
    }
    else{ /* try getting files automatically */
      CHAR efiletmp[1024], sfiletmp[1024], tfiletmp[1024];

      if( !( ttype = XLALAutoSetEphemerisFiles( efiletmp, sfiletmp, tfiletmp, pulsar,
            ifomodel->times->data[0].gpsSeconds, ifomodel->times->data[datalength-1].gpsSeconds ) ) ){
        fprintf(stderr, "Error... not been able to set ephemeris files!\n");
        exit(3);
      }

      efile = XLALStringDuplicate(efiletmp);
      sfile = XLALStringDuplicate(sfiletmp);
      tfile = XLALStringDuplicate(tfiletmp);
    }

    /* check ephemeris files exist and if not output an error message */
    if( fopen(sfile, "r") == NULL || fopen(efile, "r") == NULL ){
      fprintf(stderr, "Error... ephemeris files not, or incorrectly, defined!\n");
      exit(3);
    }

    /* set up ephemeris information */
    XLAL_CHECK_VOID( ( ifomodel->ephem = XLALInitBarycenter( efile, sfile ) ) != NULL, XLAL_EFUNC );
    if( tfile ){ XLAL_CHECK_VOID( (ifomodel->tdat = XLALInitTimeCorrections( tfile ) ) != NULL, XLAL_EFUNC ); }
    else { ifomodel->tdat = NULL; }
    ifomodel->ttype = ttype;

    XLALDestroyRandomParams( randomParams );

    /* get maximum data length */
    if ( ifodata->compTimeData->data->length > maxlen ) { maxlen = ifodata->compTimeData->data->length; }
  }

  /* chop the data into stationary chunks and also calculate the noise variance if required
   * (note that if there is going to be a signal injected then this variance will be recalculated
   * after the injection has been made to make the analysis most similar to a real case). */
  INT4 chunkMin, chunkMax;

  /* Get chunk min and chunk max */
  ppt = LALInferenceGetProcParamVal( commandLine, "--chunk-min" );
  if( ppt ) { chunkMin = atoi( ppt->value ); }
  else { chunkMin = CHUNKMIN; } /* default minimum chunk length */

  ppt = LALInferenceGetProcParamVal( commandLine, "--chunk-max" );
  if( ppt ) { chunkMax = atoi( ppt->value ); }
  else { chunkMax = CHUNKMAX; } /* default maximum chunk length */

  LALInferenceIFOData *datatmp = runState->data;
  LALInferenceIFOModel *modeltmp = runState->model->ifo;
  while ( modeltmp ){
    UINT4Vector *chunkLength = NULL;

    LALInferenceAddVariable( modeltmp->params, "chunkMin", &chunkMin, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED );
    LALInferenceAddVariable( modeltmp->params, "chunkMax", &chunkMax, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED );

    ppt = LALInferenceGetProcParamVal( commandLine, "--oldChunks" );
    if ( ppt ){ /* use old style quasi-fixed data chunk lengths */
      /* if a chunk max wasn't set use 30 mins by default */
      if ( !LALInferenceGetProcParamVal( commandLine, "--chunk-max" ) ){
        chunkMax = 30;
        LALInferenceSetVariable( modeltmp->params, "chunkMax", &chunkMax );
      }

      /* if sigma's have been input then there is just one chunk with a length of the full dataset */
      if ( inputsigma ){
        chunkLength = XLALCreateUINT4Vector( 1 );
        chunkLength->data[0] = datatmp->varTimeData->data->length;
      }
      else{ chunkLength = get_chunk_lengths( modeltmp, chunkMax ); }
    }
    /* use new change points analysis to get chunks */
    else {
      /* if sigma's have been input then there is just one chunk with a length of the full dataset */
      if ( inputsigma ){
        chunkLength = XLALCreateUINT4Vector( 1 );
        chunkLength->data[0] = datatmp->varTimeData->data->length;
      }
      else{ chunkLength = chop_n_merge( datatmp, chunkMin, chunkMax ); }
    }

    LALInferenceAddVariable( modeltmp->params, "chunkLength", &chunkLength, LALINFERENCE_UINT4Vector_t,
                             LALINFERENCE_PARAM_FIXED );

    /* compute the variance */
    if( !inputsigma ){ compute_variance( datatmp, modeltmp ); }

    /* if we have a Gaussian likelihood and have computed sigma now we can reset the
     * number of chunks to be just one chunk (so we only need to sum over one data chunk) */
    if ( computesigma ){
      UINT4Vector *chunkLengthNew = NULL;
      LALInferenceRemoveVariable( modeltmp->params, "chunkLength" );

      chunkLengthNew = XLALCreateUINT4Vector( 1 );
      chunkLengthNew->data[0] = datatmp->varTimeData->data->length;

      LALInferenceAddVariable( modeltmp->params, "chunkLength", &chunkLengthNew, LALINFERENCE_UINT4Vector_t,
                               LALINFERENCE_PARAM_FIXED );
    }

    /* set whether using Gaussian likelihood */
    if ( gaussianLike ){
      LALInferenceAddVariable( modeltmp->params, "gaussianLikelihood", &gaussianLike, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED );
    }

    datatmp = datatmp->next;
    modeltmp = modeltmp->next;
  }

  /* free memory */
  XLALFree( fdt );
  XLALFree( flengths );
  XLALFree( fstarts );
  XLALFree( fpsds );

  if ( LALInferenceCheckVariable( runState->algorithmParams, "timefile" ) ){
    gettimeofday(&time2, NULL);

    FILE *timefile = *(FILE **)LALInferenceGetVariable( runState->algorithmParams, "timefile" );
    UINT4 timenum = *(UINT4 *)LALInferenceGetVariable( runState->algorithmParams, "timenum" );
    tottime = (REAL8)((time2.tv_sec + time2.tv_usec*1.e-6) - (time1.tv_sec + time1.tv_usec*1.e-6));
    fprintf(timefile, "[%d] %s: %.9le secs\n", timenum, __func__, tottime);
    timenum++;
    check_and_add_fixed_variable( runState->algorithmParams, "timenum", &timenum, LALINFERENCE_UINT4_t );
  }
}


/**
 * \brief Reads in the parameters of the pulsar being searched for
 *
 * This function reads in a pulsars parameters from the specified TEMPO-style .par file given by \c par-file using \c
 * XLALReadTEMPOParFile. This file must be specified and should contain at least the pulsars frequency, right
 * ascension and declination (any value not included will be zero by default). The file should contain the parameters
 * with which the detector data was heterodyned, as these are used to produce a signal phase template based on this
 * assumption.
 *
 * A example .par file may look like
 * \code
 * RA 12:31:56.17643
 * DEC 43:21:35.2531
 * F0 100.78634 1 0.00005
 * F1 2.34e-15
 * PEPOCH 54323.785634
 * \endcode
 * which shows several parameters mostly defined by the parameter name and a parameter value. However, the \c F0 value
 * contains 4 items. If a parameter has a \c 1 as the third entry then it means that this was a parameter that was fit
 * by TEMPO with the entry after the \c 1 being the 1 standard deviation error on that parameter. For parameters where
 * an error is present the code will attempt to search over that parameter using a Gaussian prior defined by the
 * 1\f$\sigma\f$ error value. Other parameters will be set as fixed by default. These can be overridden by the prior
 * file values described in \c initialise_prior().
 *
 * Based on the defined sky position defined in the par file a lookup table of the detector antenna response over time
 * and polarisation will be set by \c setup_lookup_tables().
 *
 * The function \c add_initial_variables() is used to pass the parameter values from the .par file to the algorithm.
 *
 * Using the parameters from the .par file the phase template, including the solar system and binary system barycentring
 * time delays will be setup. These define the phase template used to perform the initial heterodyne, which is used as
 * the reference in cases when phase parameters (other than the initial phase) are being searched over.
 *
 * Values used for scaling the parameters (to avoid dynamic range issues) are initialised although will be set as
 * default values.
 *
 * \param runState [in] A pointer to the LALInferenceRunState
 *
 * \sa setup_lookup_tables
 * \sa add_initial_variables
 * \sa get_phase_model
 * \sa add_correlation_matrix
 */
void setup_from_par_file( LALInferenceRunState *runState )
/* Read the PAR file of pulsar parameters and setup the code using them */
/* Generates lookup tables also */
{
  LALSource psr;
  PulsarParameters *pulsar;
  REAL8Vector *phase_vector = NULL;
  LALInferenceIFOData *data = runState->data;
  ProcessParamsTable *ppt = NULL;
  REAL8 DeltaT = 0.; /* maximum data time span */

  ppt = LALInferenceGetProcParamVal( runState->commandLine, "--par-file" );
  if( ppt == NULL ) { fprintf(stderr,"Must specify --par-file!\n"); exit(1); }
  CHAR *parFile = ppt->value;

  /* get the pulsar parameters */
  pulsar = XLALReadTEMPOParFileNew( parFile );

  REAL8 ra = 0.;
  if ( PulsarCheckParam( pulsar, "RA" ) ) { ra = PulsarGetREAL8Param( pulsar, "RA" ); }
  else if ( PulsarCheckParam( pulsar, "RAJ" ) ) { ra = PulsarGetREAL8Param( pulsar, "RAJ" ); }
  else {
    XLALPrintError ("%s: No source right ascension specified!", __func__ );
    XLAL_ERROR_VOID( XLAL_EINVAL );
  }
  REAL8 dec = 0.;
  if ( PulsarCheckParam( pulsar, "DEC" ) ) { dec = PulsarGetREAL8Param( pulsar, "DEC" ); }
  else if ( PulsarCheckParam( pulsar, "DECJ" ) ) { dec = PulsarGetREAL8Param( pulsar, "DECJ" ); }
  else {
    XLALPrintError ("%s: No source declination specified!", __func__ );
    XLAL_ERROR_VOID( XLAL_EINVAL );
  }
  psr.equatorialCoords.longitude = ra;
  psr.equatorialCoords.latitude = dec;
  psr.equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;

  /* Setup lookup tables for amplitudes */
  setup_lookup_tables( runState, &psr );

  runState->model->params = XLALCalloc( 1, sizeof(LALInferenceVariables) );
  runState->model->domain = LAL_SIM_DOMAIN_TIME;

  runState->currentParams = XLALCalloc( 1, sizeof(LALInferenceVariables) );

  /* Add initial (unchanging) variables for the model from the par file */
  add_initial_variables( runState->currentParams, pulsar );

  /* check for binary model */
  CHAR *binarymodel = NULL;
  if ( LALInferenceCheckVariable( runState->currentParams, "BINARY") ){
    binarymodel = XLALStringDuplicate(*(CHAR**)LALInferenceGetVariable( runState->currentParams, "BINARY" ));

    /* now remove from runState->params (as it conflict with calls to LALInferenceCompareVariables in the proposal) */
    LALInferenceRemoveVariable( runState->currentParams, "BINARY" );
  }

  /* Setup initial phase, and barycentring delays */
  LALInferenceIFOModel *ifo_model = runState->model->ifo;
  while( data ){
    REAL8Vector *freqFactors = NULL;
    UINT4 j = 0;
    REAL8 dt = XLALGPSGetREAL8( &ifo_model->times->data[ifo_model->times->length-1] ) - XLALGPSGetREAL8( &ifo_model->times->data[0] );

    if ( dt > DeltaT ){ DeltaT = dt; }

    freqFactors = *(REAL8Vector **)LALInferenceGetVariable( ifo_model->params, "freqfactors" );

    for( j = 0; j < freqFactors->length; j++ ){
      UINT4 i = 0;
      REAL8Vector *dts = NULL, *bdts = NULL;

      /* check whether using original Jones (2010) signal model or a biaxial model (in the amplitude/phase parameterisation) */
      if ( freqFactors->length == 2 ){
        UINT4 dummyvar = 1;

        if ( LALInferenceGetProcParamVal( runState->commandLine, "--jones-model" ) ){
          LALInferenceAddVariable( ifo_model->params, "jones-model", &dummyvar, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED );
        }
        else if ( LALInferenceGetProcParamVal( runState->commandLine, "--biaxial" ) ){
          LALInferenceAddVariable( ifo_model->params, "biaxial", &dummyvar, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED );
        }
      }

      /* add binary model to the general parameters */
      if ( binarymodel != NULL ){
        LALInferenceAddVariable( ifo_model->params, "BINARY", &binarymodel, LALINFERENCE_string_t, LALINFERENCE_PARAM_FIXED );
      }

      dts = get_ssb_delay( pulsar, ifo_model->times, ifo_model->ephem, ifo_model->tdat, ifo_model->ttype, data->detector, 0. );
      LALInferenceAddVariable( ifo_model->params, "ssb_delays", &dts, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );

      bdts = get_bsb_delay( pulsar, ifo_model->times, dts, ifo_model->ephem );
      LALInferenceAddVariable( ifo_model->params, "bsb_delays", &bdts, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );

      phase_vector = get_phase_model( pulsar, ifo_model, freqFactors->data[j] );

      ifo_model->timeData = NULL;
      ifo_model->timeData = XLALCreateREAL8TimeSeries( "", &ifo_model->times->data[0], 0., 1., &lalSecondUnit, phase_vector->length );

      for ( i=0; i<phase_vector->length; i++ ) { ifo_model->timeData->data->data[i] = phase_vector->data[i]; }

      data = data->next;
      ifo_model = ifo_model->next;
    }
  }

  /* set frequency bin step from longest data time span */
  REAL8 df = 1./(2.*DeltaT);
  LALInferenceAddVariable( runState->currentParams, "df", &df, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );

  return;
}


/**
 * \brief Read in an ascii text file of nested samples, convert to posterior samples and create k-d tree
 *
 * This function reads in an ascii file of nested samples, converted them into posterior samples and them add them to a
 * k-d tree. The file name containing the samples must have been given as the command line argument \c sample-file and
 * there must be an accompanying file with the names of each column with the same file name with _params.txt appended.
 *
 * It is assumed that the samples are in ascending log likelihood order. It is also assumed that variable values in the
 * file (and are not likelihood-like values) are consistent with those given that have prior ranges defined in the prior
 * file/par file (as these ranges will be used as bounds in a k-d tree created from this data).
 *
 * As it is assumed that the points read in are from a previous nested sampling run the number of live points used for
 * that run are also required to be given with the \c sample-nlive argument. This will be used during the conversion to
 * posterior samples.
 *
 * If given the k-d tree cell size for using the posterior as a prior can be set with the \c prior-cell argument, if not
 * set this defaults to 32.
 *
 * In the future this will be altered so as to also read in an XML file of samples.
 *
 * NOTE: add the ability to read in multiple files and combine the posterior samples
 *
 * \param runState [in] A pointer to the LALInferenceRunState
 */
void samples_prior( LALInferenceRunState *runState ){
  ProcessParamsTable *ppt = NULL;

  UINT4 Ncell = 8; /* default prior cell size */

  UINT4 i = 0, k = 0, nsamps = 0, nnlive = 0, n = 0;
  UINT4Vector *nlive = NULL, *Nsamps = NULL;
  CHAR *nlivevals = NULL, *templives = NULL, *templive = NULL;

  CHAR *sampfile = NULL;
  CHAR *tempsamps = NULL, *tempsamp = NULL;

  LALStringVector *sampfilenames = NULL;

  LALInferenceVariables ***params = NULL;

  FILE *fp = NULL;

  const CHAR *fn = __func__;

  /* get names of nested sample file columns */
  ppt = LALInferenceGetProcParamVal( runState->commandLine, "--sample-files" );
  if ( ppt != NULL ){
    sampfile = XLALStringDuplicate( ppt->value );

    /* count the number of sample files from the comma separated vales and set their names */
    tempsamps = XLALStringDuplicate( sampfile );

    nsamps = count_csv( tempsamps );

    for( i = 0; i < nsamps; i++ ){
      tempsamp = XLALStringToken( &tempsamps, "," , 0);
      sampfilenames = XLALAppendString2Vector( sampfilenames, tempsamp );
    }
  }
  else return; /* no file so we don't use this function */

  ppt = LALInferenceGetProcParamVal( runState->commandLine, "--sample-nlives" );
  if ( ppt != NULL ){
    nlivevals = XLALStringDuplicate( ppt->value );

    templives = XLALStringDuplicate( nlivevals );

    nnlive = count_csv( templives );

    if( nnlive != nsamps ){
      XLALPrintError("%s: Number of live points not equal to number of posterior files!\n", fn);
      XLAL_ERROR_VOID(XLAL_EIO);
    }

    for( i = 0; i < nnlive; i++ ){
      templive = XLALStringToken( &templives, "," , 0);
      nlive = XLALResizeUINT4Vector( nlive, i+1 );
      nlive->data[i] = atoi( templive );
    }

    LALInferenceAddVariable( runState->algorithmParams, "numberlive", &nlive, LALINFERENCE_UINT4Vector_t,
                             LALINFERENCE_PARAM_FIXED );
  }
  else{
    fprintf(stderr, "Must set the number of live points used in the input nested samples file.\n\n");
    fprintf(stderr, USAGE, runState->commandLine->program);
    exit(0);
  }

  /* allocate memory for nested samples */
  params = XLALCalloc( nsamps, sizeof(LALInferenceVariables**) );

  /* loop over files, convert to posterior samples and combine them */
  for ( n = 0; n < nsamps; n++ ){
    CHAR *namefile = NULL, name[256];
    LALStringVector *paramNames = NULL;

    i = 0;

    /* initialise array as NULL */
    params[n] = NULL;

    namefile = XLALStringDuplicate( sampfilenames->data[n] );
    namefile = XLALStringAppend( namefile, "_params.txt" );

    /* check file exists */
    if ( fopen(namefile, "r") == NULL || fopen(sampfilenames->data[n], "r") == NULL ){
      XLALPrintError("%s: Cannot access either %s or %s!\n", fn, namefile, sampfilenames->data[n]);
      XLAL_ERROR_VOID(XLAL_EIO);
    }

    /* read in parameter names */
    fp = fopen( namefile, "r" );
    while( fscanf(fp, "%s", name) != EOF ){ paramNames = XLALAppendString2Vector( paramNames, name ); }

    fclose(fp);

    /* read in parameter values */
    fp = fopen( sampfilenames->data[n], "r" );
    while( !feof(fp) ){
      REAL8 ps[paramNames->length];
      UINT4 j = 0;

      for( j=0; j<paramNames->length; j++ ) { if( fscanf(fp, "%lf", &ps[j]) == EOF ) { break; } }

      if( feof(fp) ) { break; }

      /* dynamically allocate memory */
      params[n] = XLALRealloc( params[n], (i+1)*sizeof(LALInferenceVariables*) );
      params[n][i] = NULL;
      params[n][i] = XLALCalloc( 1, sizeof(LALInferenceVariables) );

      /* add variables */
      for( j=0; j<paramNames->length; j++ ){
        /* use vary type of this analyses parameters i.e. those set by the prior
           and par file, otherwise set the parameter to fixed */
        LALInferenceParamVaryType vary;

        if ( LALInferenceCheckVariable( runState->currentParams, paramNames->data[j] ) ){
          vary = LALInferenceGetVariableVaryType( runState->currentParams, paramNames->data[j] );
        }
        else { vary = LALINFERENCE_PARAM_FIXED; }

        LALInferenceAddVariable( params[n][i], paramNames->data[j], &ps[j], LALINFERENCE_REAL8_t, vary );
      }

      i++;
    }

    /* check that non-fixed, or output parameters actually do vary, otherwise
      complain */
    LALInferenceVariableItem *item1 = params[n][0]->head;

    while ( item1 ){
      UINT4 allsame = 0;

      for ( k=1; k<i; k++ ){
        LALInferenceVariableItem *item2 = LALInferenceGetItem( params[n][k], item1->name );

        if( item1->vary != LALINFERENCE_PARAM_FIXED && item1->vary != LALINFERENCE_PARAM_OUTPUT ){
          if ( *(REAL8*)item1->value != *(REAL8*)item2->value ) { allsame++; }
        }
      }

      if( ( item1->vary != LALINFERENCE_PARAM_FIXED && item1->vary != LALINFERENCE_PARAM_OUTPUT ) && allsame == 0 ){
        XLALPrintError("%s: Apparently variable parameter %s does not vary!\n", fn, item1->name );
        XLAL_ERROR_VOID(XLAL_EFUNC);
      }

      item1 = item1->next;
    }

    Nsamps = XLALResizeUINT4Vector( Nsamps, n+1 );
    Nsamps->data[n] = i;
  }

  LALInferenceAddVariable( runState->algorithmParams, "nestedsamples", &params, LALINFERENCE_void_ptr_t,
                           LALINFERENCE_PARAM_FIXED );
  LALInferenceAddVariable( runState->algorithmParams, "Nsamps", &Nsamps, LALINFERENCE_UINT4Vector_t,
                           LALINFERENCE_PARAM_FIXED );

  /* get cell size */
  ppt = LALInferenceGetProcParamVal( runState->commandLine, "--prior-cell" );
  if ( ppt != NULL ) { Ncell = atoi( ppt->value ); }

  LALInferenceAddVariable( runState->priorArgs, "kDTreePriorNcell", &Ncell, LALINFERENCE_UINT4_t,
                           LALINFERENCE_PARAM_FIXED );

  /* convert samples to posterior */
  ns_to_posterior( runState );

  /* create k-d tree of the samples for use as a prior */
  create_kdtree_prior( runState );
}

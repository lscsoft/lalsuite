/*
*  Copyright (C) 2007 Matt Pitkin
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

/*******************************************************************************
  Matt Pitkin - 08/08/07

  pulsar_parameter_estimation.c

  This code will be used to produce a posterior probability distribution for
  the four unknown pulsar parameters:  h0, cosiota, psi and phi0. The code can
  output either the full posterior or the marginalised posteriors. This can be
  for multiple detector or a single detector. The log(evidence) or "marginalised
  posterior" can also be output.

  This code can also use a Markov Chain Monte Carlo to produce the posterior.

*******************************************************************************/

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include "pulsar_parameter_estimation.h"

/* global variable */
INT4 verbose=0;

static char joint_string[] = "Joint";
static char uniform_string[] = "uniform";
const OutputParams empty_OutputParams;

/* Usage format string */
static char USAGE1[] = \
" --help              display this message\n"\
" --verbose           display all error messages\n"\
" --detectors         all IFOs with data to be analysed e.g. H1,H2\n\
                     (delimited by commas)\n"\
" --pulsar            name of pulsar e.g. J0534+2200\n"\
" --par-file          pulsar parameter (.par) file (full path) \n"\
" --input-dir         directory containing the input data files\n"\
" --output-dir        directory for output data files\n"\
" --dob-ul            (REAL8) percentage degree-of-belief for upper limit\n\
                     - if 0 (default no upper limit is produced)\n"\
" --output-post       output the full log(posterior)\n"\
" --chunk-min         (INT4) minimum stationary length of data to be used in\n\
                     the likelihood e.g. 5 mins\n"\
" --chunk-max         (INT4) maximum stationary length of data to be used in\n\
                     the likelihood e.g. 30 mins\n"\
"\n"\
" Parameter space grid values:-\n"\
" --minh0             (REAL8) minimum of the h0 grid\n"\
" --maxh0             (REAL8) maximum of the h0 grid, if maxh0=0 then\n\
                     calculate range from the data\n"\
" --h0steps           (INT4) number of steps in the h0 grid\n"\
" --minphi0           (REAL8) minimum of the phi0 grid\n"\
" --maxphi0           (REAL8) maximum of the phi0 grid\n"\
" --phi0steps         (INT4) number of steps in the phi0 grid\n"\
" --minpsi            (REAL8) minimum of the psi grid\n"\
" --maxpsi            (REAL8) maximum of the psi grid\n"\
" --psisteps          (INT4) number of steps in the psi grid\n"\
" --minci             (REAL8) minimum of the cos(iota) grid\n"\
" --maxci             (REAL8) maximum of the cos(iota) grid\n"\
" --cisteps           (INT4) number of steps in the cos(iota) grid\n"\
" --psi-bins          (INT4) no. of psi bins in the time-psi lookup table\n"\
" --time-bins         (INT4) no. of time bins in the time-psi lookup table\n"\
"\n"\
" Prior values:-\n"\
" --use-priors        set to use priors in the posterior calculation\n"\
" --h0prior           type of prior on h0 - uniform, jeffreys or gaussian\n"\
" --h0mean            (REAL8) mean of a gaussian prior on h0\n"\
" --h0sig             (REAL8) standard deviation of a gaussian prior on h0\n"\
" --priorfile         (string) a binary file containing a h0 x cos(iota)\n\
                     prior probability distribution. The file should contain\n\
                     a header of six doubles with:\n\
                       - the minimum value for h0\n\
                       - dh0 - the step size in h0\n\
                       - N - number of h0 values\n\
                       - the minimum value for cos(iota)\n\
                       - dci - the step size in cos(iota)\n\
                       - M - number of cos(iota) values\n\
                     followed by an NxM double array of posterior values.\n"\
" --phi0prior         type of prior on phi0 - uniform or gaussian\n"\
" --phi0mean          (REAL8) mean of a gaussian prior on phi0\n"\
" --phi0sig           (REAL8) std. dev. of a gaussian prior on phi0\n"\
" --psiprior          type of prior on psi - uniform or gaussian\n"\
" --psimean           (REAL8) mean of a gaussian prior on psi\n"\
" --psisig            (REAL8) std. dev. of a gaussian prior on psi\n"\
" --iotaprior         type of prior on iota - uniform or gaussian\n"\
" --iotamean          (REAL8) mean of a gaussian prior on iota\n"\
" --iotasig           (REAL8) std. dev. of a gaussian prior on iota\n"\
"\n";

static char USAGE2[] = \
" MCMC parameters:-\n"\
" --mcmc              set to perform an MCMC\n"\
" --iterations        (INT4) the number of iteraction in the MCMC chain\n"\
" --burn-in           (INT4) the number of burn in iterations\n"\
" --temperature       (REAL8) the temperatue to start of the simulated\n\
                     annealing in the burn in stage\n"\
" --h0-width          (REAL8) width of the h0 proposal distribution (if set\n\
                     to 0 this will be worked out in the code)\n"\
" --h0-scale          (REAL8) scale factor for h0 proposal width\n"\
" --psi-width         (REAL8) width of the psi proposal distribution\n"\
" --phi0-width        (REAL8) width of the phi proposal distribution\n"\
" --ci-width          (REAL8) width of the cos(iota) proposal distribution\n"\
" --output-rate       (INT4) rate at which to output chain e.g. 10 means\n\
                     output every tenth sample\n"\
" --nglitch           (INT4) number of glitches\n"\
" --glitch-times      (CHAR) a string of pulsar glitch times given in MJD\n\
                     e.g. 45623.872,52839.243,53992.091 - at each\n\
                     glitch an additional MCMC phase parameter will be used\n"\
" --glitch-cut        (REAL8) the number of seconds of data to ignore\n\
                     before and after a glitch\n"\
" --earth-ephem       Earth ephemeris file\n"\
" --sun-ephem         Sun ephemeris file\n"\
" --use-cov           if this flag is set and/or a covariance file is\n\
                     specified (with --covariance) then that will be used, if\n\
                     just this flag is set then a covariance matrix\n\
                     constructed from standard deviations in the par file\n\
                     (with zero off diagonal elements) will be used\n"\
" --covariance        pulsar parameter covariance matrix file (.mat)\n"\
" --only-joint        set this to only produce the joint MCMC when given \n\
                     muliple detectors (MCMC only)\n"\
" --output-burn-in    set this to also output the burn in stage to the chain\n"\
"\n";




INT4 main(INT4 argc, CHAR *argv[]){
  REAL8 ****singleLike=NULL;
  REAL8 ****jointLike=NULL;

  InputParams inputs;
  BinaryPulsarParams pulsar;
  LALDetAndSource detAndSource;

  INT4 numDets=0; /* number of detectors */
  CHAR dets[5][3]; /* we'll have a max of five detectors */
  LALDetector detPos[5];

  INT4 i=0, j=0, k=0, n=0;

  DataStructure *data=NULL;
  REAL8 times=0.;
  COMPLEX16 dataVals;
  REAL8 stdh0=0.;       /* approximate h0 limit from data */

  FILE *fp=NULL;
  CHAR dataFile[256];
  CHAR outputFile[256];

  OutputParams output = empty_OutputParams;
  REAL8 logNoiseEv[5]; /* log evidence for noise only (no signal) */
  Results results;
  REAL8 h0ul=0.;
  REAL8 maxPost=0.;
  
  CHAR params[][10]={"h0", "phi", "psi", "ciota"};

  EphemerisData *edat=NULL;

  /*===================== GET AND SET THE INPUT PARAMETER ====================*/
  get_input_args(&inputs, argc, argv);

  /* if we want to output in verbose mode set global variable */
  if(inputs.verbose) verbose = 1;

  /* get the pulsar parameters */
  XLALReadTEMPOParFile(&pulsar, inputs.parFile);
  inputs.psr.equatorialCoords.longitude = pulsar.ra;
  inputs.psr.equatorialCoords.latitude = pulsar.dec;
  inputs.psr.equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;

  /* find the number of detectors being used */
  if( strstr(inputs.detectors, "H1") != NULL ){
    sprintf(dets[numDets], "H1");
    detPos[numDets] = *XLALGetSiteInfo( dets[numDets] );
    numDets++;
  }
  if( strstr(inputs.detectors, "H2") != NULL ){
     sprintf(dets[numDets], "H2");
     detPos[numDets] = *XLALGetSiteInfo( dets[numDets] );
     numDets++;
  }
  if( strstr(inputs.detectors, "L1") != NULL ){
     sprintf(dets[numDets], "L1");
     detPos[numDets] = *XLALGetSiteInfo( dets[numDets] );
     numDets++;
  }
  if( strstr(inputs.detectors, "G1") != NULL ){
     sprintf(dets[numDets], "G1");
     detPos[numDets] = *XLALGetSiteInfo( dets[numDets] );
     numDets++;
  }
  if( strstr(inputs.detectors, "V1") != NULL ){
     sprintf(dets[numDets], "V1");
     detPos[numDets] = *XLALGetSiteInfo( dets[numDets] );
     numDets++;
  }

  if( verbose ){
    fprintf(stderr, "Analysing data from %d detector(s):\n  ", numDets);
    for( i = 0 ; i < numDets ; i++ )
      fprintf(stderr, "%s ", dets[i]);
    fprintf(stderr, "\n");
  }
  /*==========================================================================*/

  /*====================== SET OUTPUT PARAMETERS =============================*/
  output.outputDir = inputs.outputDir;
  output.psr = inputs.pulsar;
  output.dob = inputs.dob; /* set degree of belief for UL */
  output.outPost = inputs.outputPost;
  sprintf(outputFile, "%s/evidence_%s", output.outputDir, output.psr);
  /*==========================================================================*/

  if( inputs.mcmc.doMCMC == 0 ){
    /* allocate likelihood memory */
    singleLike = allocate_likelihood_memory(inputs.mesh);

    /* if more than one detector create a joint likelihood */
    if( numDets > 1 )
      jointLike = allocate_likelihood_memory(inputs.mesh);

    if( verbose )
      fprintf(stderr, "I've allocated the memory for the likelihood.\n");

    /* if we're doing the grid search we only need to store one data set at a
       time */
    data = XLALMalloc(sizeof(DataStructure));
  }
  else{
    /* if we're doing the MCMC we need to store all the Bks */
    data = XLALMalloc(numDets*sizeof(DataStructure));

    /* if there's a covariance matrix file then set up the earth and sun
       ephemeris */
    if( inputs.matrixFile != NULL || inputs.usecov ){
      /* check files exist and if not output an error message */
      if( fopen(inputs.earthfile, "r") == NULL || 
          fopen(inputs.sunfile, "r") == NULL ){
        fprintf(stderr, "Error... ephemeris files not, or incorrectly, \
defined!\n");
        return 0;
      }

      XLAL_CHECK( ( edat = XLALInitBarycenter( inputs.earthfile, 
                    inputs.sunfile ) ) != NULL, XLAL_EFUNC );
    }
    else
      edat = NULL;
  }

  k = -1;

  /* read in data for each detector in turn an compute the likelihood */
  for( i = 0 ; i < numDets ; i++ ){
    /*============================ GET DATA ==================================*/
    /* get detector B_ks data file in form finehet_JPSR_DET */
    sprintf(dataFile, "%s/data%s/finehet_%s_%s", inputs.inputDir, dets[i],
      inputs.pulsar, dets[i]);

    /* open data file */
    if((fp = fopen(dataFile, "r"))==NULL){
      fprintf(stderr, "Error... can't open data file %s!\n", dataFile);
      exit(0);
    }

    j=0;

    /* only store one set of data for grid search (saves memory), but store
       them all for the MCMC */
    if( inputs.mcmc.doMCMC == 0 ) k = 0;
    else k++;

    /* read in data */
    data[k].data = NULL;
    data[k].data = XLALCreateCOMPLEX16Vector( MAXLENGTH );
    data[k].times = NULL;
    data[k].times = XLALCreateREAL8Vector( MAXLENGTH );

    /* set the minimum and maximum data segments lengths */
    data[k].chunkMin = inputs.chunkMin;
    data[k].chunkMax = inputs.chunkMax;

    stdh0 = 0.;
    /* read in data */
    while(fscanf(fp, "%lf%lf%lf", &times, &dataVals.re, &dataVals.im) != EOF){
      /* check that size of data file is not to large */
      if( j == MAXLENGTH ){
        fprintf(stderr, "Error... size of MAXLENGTH not large enough.\n");
        exit(0);
      }

      /* exclude values smaller than 1e-28 as most are spurious points caused
         during a the heterodyne stage (e.g. when frame files were missing in
         the original S5 analysis) */
      /* if( fabs(dataVals.re) > 1.e-28 && fabs(dataVals.im) > 1.e-28 ){ */
        data[k].times->data[j] = times;
        data[k].data->data[j] = dataVals;

        /* get the power from the time series */
        stdh0 += dataVals.re*dataVals.re + dataVals.im*dataVals.im;

        j++;
      /*}*/
    }

    fclose(fp);

    if( verbose )
      fprintf(stderr, "I've read in the data for %s.\n", dets[i]);

    data[k].data = XLALResizeCOMPLEX16Vector(data[k].data, j);
    data[k].times = XLALResizeREAL8Vector(data[k].times, j);

    /* if there is no input range for h0 then estimate it from the data */
    /* only do this once if performing grid search, but do for each seperate
       data set if doing MCMC */
    if( ( inputs.mesh.maxVals.h0 == 0 || inputs.mcmc.sigmas.h0 == 0 ) && 
        ( inputs.mcmc.doMCMC == 1 || i == 0 ) ){
      if( verbose ) fprintf(stderr, "Calculating h0 UL estimate: ");

      /* get the power spectral density power/bandwidth (1/60 Hz) */
      stdh0 = stdh0/((REAL8)j*(1./60.));

      /* upper limit estimate comes from ~ h0 = 10.8*sqrt(Sn/T) */
      stdh0 = 10.8*sqrt(stdh0/((REAL8)j*60.));

      /* set the MCMC h0 proposal step size at stdh0*scalefac */
      if( inputs.mcmc.doMCMC == 1 ){
        inputs.mcmc.sigmas.h0 = stdh0*inputs.mcmc.h0scale;
        
        if( inputs.mesh.maxVals.h0 == 0 )
          inputs.mesh.maxVals.h0 = stdh0;
      }

      /* set h0 max value for the grid at 5 times the expected ul */
      if( inputs.mesh.maxVals.h0 == 0 ){
        inputs.mesh.maxVals.h0 = 5.*stdh0;
        inputs.mesh.delta.h0 = (inputs.mesh.maxVals.h0 -
          inputs.mesh.minVals.h0)/(REAL8)(inputs.mesh.h0Steps - 1.);
      }

      if( verbose ) fprintf(stderr, "%le\n", stdh0);
    }

    /*========================================================================*/

    output.det = dets[i];

    /* create lookup table */
    data[k].lookupTable = NULL;
    data[k].lookupTable = XLALCalloc(1, sizeof(DetRespLookupTable));
    data[k].lookupTable->lookupTable=NULL;
    detAndSource.pSource = &inputs.psr;
    detAndSource.pDetector = &detPos[i];

    /* create memory for the lookup table */
    data[k].lookupTable->lookupTable = XLALCalloc(inputs.mesh.psiRangeSteps, 
      sizeof(LALDetAMResponse *));

    for( j = 0 ; j < inputs.mesh.psiRangeSteps ; j++ ){
      data[k].lookupTable->lookupTable[j] =
        XLALCalloc(inputs.mesh.timeRangeSteps, sizeof(LALDetAMResponse));
    }

    data[k].lookupTable->psiSteps = inputs.mesh.psiRangeSteps;
    data[k].lookupTable->timeSteps = inputs.mesh.timeRangeSteps;

    /* create lookup table */
    response_lookup_table(data[k].times->data[0], detAndSource,
      data[k].lookupTable);

    if( verbose ) fprintf(stderr, "Created look-up table.\n");

    /*======================== CALCULATE LIKELIHOOD ==========================*/
    if( inputs.mcmc.doMCMC == 0 ){
      logNoiseEv[i] = create_likelihood_grid(data[k], singleLike, inputs.mesh);

      if( verbose )
        fprintf(stderr, "I've calculated the likelihood for %s.\n", dets[i]);

      /* if there's more than one detector calculate the joint likelihood */
      if( numDets > 1 ){
        /* add the single detector log likelihood onto the joint likelihood */
        combine_likelihoods(singleLike, jointLike, inputs.mesh);
        output.outPost = 0; /* don't output individual detector posteriors */
      }
      /*======================================================================*/

      /*========== CREATE THE SINGLE DETECTOR POSTERIORS =====================*/
      maxPost = log_posterior(singleLike, inputs.priors, inputs.mesh, output);

      if( isinf(maxPost) ){
        fprintf(stderr, "Error... posterior is infinite!\n");
        return 0;
      }
      
      /* marginalise over each parameter and output the data */
      for( n = 0 ; n < 4 ; n++ ){
        output.margParam = params[n];
        if( verbose )
          fprintf(stderr, "Marginalising over %s.\n", output.margParam);

        results = marginalise_posterior(singleLike, inputs.mesh, output);

        if( output.dob != 0. && strcmp( output.margParam, "h0" ) == 0)
          h0ul = results.h0UpperLimit;
      }

      /* open file to output the evidence */
      if((fp = fopen(outputFile, "a"))==NULL){
        fprintf(stderr, "Error... can't open file %s!\n", outputFile);
        return 0;
      }

      /* output the log odds ratio and an UL if requested */
      fprintf(fp, "%s\t%le\n", output.det, results.evidence-logNoiseEv[i]);
      if( output.dob != 0. )
        fprintf(fp, "%.1lf%% h0 upper limit = %le\n", output.dob, h0ul);

      fclose(fp);

      XLALDestroyCOMPLEX16Vector(data[k].data);
      XLALDestroyREAL8Vector(data[k].times);
    }
    /*========================================================================*/

    /*================== PERFORM THE MCMC ====================================*/
    if( inputs.mcmc.doMCMC != 0 && inputs.onlyjoint == 0 )
      perform_mcmc(&data[k], inputs, 1, dets[i], &detPos[i], edat);

    /*========================================================================*/
  }

  /*=================== CREATE THE JOINT POSTERIOR IF REQUIRED ===============*/
  if( numDets > 1 ){
    output.det = joint_string;

    if( inputs.mcmc.doMCMC == 0 ){
      REAL8 totLogNoiseEv=0.;

      for( n = 0 ; n < numDets ; n++ ) totLogNoiseEv += logNoiseEv[n];

      output.outPost = inputs.outputPost; /* set for whether we want to output
                                            the full posterior */

      maxPost = log_posterior(jointLike, inputs.priors, inputs.mesh, output);
      
      if( isinf(maxPost) ){
        fprintf(stderr, "Error... posterior is infinite!\n");
        return 0;
      }
      
      if( verbose )
        fprintf(stderr, "I've calculated the joint posterior.\n");

      /* marginalise over each parameter and output the data */
      for( n = 0 ; n < 4 ; n++ ){
        output.margParam = params[n];
        if( verbose )
          fprintf(stderr, "Marginalising over %s.\n", output.margParam);

        results = marginalise_posterior(jointLike, inputs.mesh, output);
      }

      /* open file to output the evidence */
      if((fp = fopen(outputFile, "a"))==NULL){
        fprintf(stderr, "Error... can't open file %s!\n", outputFile);
        return 0;
      }

      /* output the evidence */
      fprintf(fp, "%s\t%le\n", output.det, results.evidence - totLogNoiseEv);
      if( output.dob != 0. ){
        fprintf(fp, "%.1lf%% h0 upper limit = %le\n", output.dob,
          results.h0UpperLimit);
      }
      fclose(fp);
    }
    /*========================================================================*/

    /*======================= PERFORM JOINT MCMC =============================*/
    if( inputs.mcmc.doMCMC == 1 ){
      perform_mcmc(data, inputs, numDets, output.det, detPos, edat);

      /* destroy data */
      for( i = 0 ; i < numDets ; i++ ){
        XLALDestroyCOMPLEX16Vector(data[i].data);
        XLALDestroyREAL8Vector(data[i].times);
      }
    }
  }

  /*====================== FREE THE LIKELIHOOD MEMORY ========================*/
  if( inputs.mcmc.doMCMC == 0 ){
    for( i = 0 ; i < inputs.mesh.phiSteps ; i++ ){
      for( j = 0 ; j < inputs.mesh.ciotaSteps ; j++ ){
        if( numDets > 1 )
          XLALFree(jointLike[i][j]);

        XLALFree(singleLike[i][j]);
      }
      if( numDets > 1 )
        XLALFree(jointLike[i]);

      XLALFree(singleLike[i]);
    }

    if( numDets > 1 )
      XLALFree(jointLike);

    XLALFree(singleLike);
  }
  /*=========================================================================*/ 

  /* free Ephemeris data */
  if ( edat != NULL ) XLALDestroyEphemerisData( edat );
  
  return 0;
}



/* function to get the input arguments from the command line */
void get_input_args(InputParams *inputParams, INT4 argc, CHAR *argv[]){
  struct option long_options[] =
  {
    { "help",           no_argument,       0, 'h' },
    { "verbose",        no_argument,    NULL, 'R' },
    { "detectors",      required_argument, 0, 'D' },
    { "pulsar",         required_argument, 0, 'p' },
    { "par-file",       required_argument, 0, 'P' },
    { "input-dir",      required_argument, 0, 'i' },
    { "output-dir",     required_argument, 0, 'o' },
    { "minh0",          required_argument, 0, 'a' },
    { "maxh0",          required_argument, 0, 'A' },
    { "h0steps",        required_argument, 0, 'j' },
    { "minphi0",        required_argument, 0, 'b' },
    { "maxphi0",        required_argument, 0, 'B' },
    { "phi0steps",      required_argument, 0, 'k' },
    { "minpsi",         required_argument, 0, 's' },
    { "maxpsi",         required_argument, 0, 'S' },
    { "psisteps",       required_argument, 0, 'm' },
    { "minci",          required_argument, 0, 'c' },
    { "maxci",          required_argument, 0, 'C' },
    { "cisteps",        required_argument, 0, 'n' },
    { "psi-bins",       required_argument, 0, 'l' },
    { "time-bins",      required_argument, 0, 'L' },
    { "h0prior",        required_argument, 0, 'q' },
    { "priorfile",      required_argument, 0, ']' },
    { "phi0prior",      required_argument, 0, 'Q' },
    { "psiprior",       required_argument, 0, 'U' },
    { "iotaprior",      required_argument, 0, 'u' },
    { "h0mean",         required_argument, 0, 'Y' },
    { "h0sig",          required_argument, 0, 'T' },
    { "phi0mean",       required_argument, 0, 'v' },
    { "phi0sig",        required_argument, 0, 'V' },
    { "psimean",        required_argument, 0, 'z' },
    { "psisig",         required_argument, 0, 'Z' },
    { "iotamean",       required_argument, 0, 'e' },
    { "iotasig",        required_argument, 0, 'E' },
    { "output-post",    no_argument,    NULL, 'f' },
    { "dob-ul",         required_argument, 0, 'd' },
    { "mcmc",           no_argument,    NULL, 'F' },
    { "iterations",     required_argument, 0, 'I' },
    { "burn-in",        required_argument, 0, 'x' },
    { "temperature",    required_argument, 0, 't' },
    { "h0-width",       required_argument, 0, 'H' },
    { "psi-width",      required_argument, 0, 'w' },
    { "phi0-width",     required_argument, 0, 'W' },
    { "ci-width",       required_argument, 0, 'y' },
    { "glitch-times",   required_argument, 0, 'g' },
    { "glitch-cut",     required_argument, 0, 'G' },
    { "chunk-min",      required_argument, 0, 'K' },
    { "chunk-max",      required_argument, 0, 'N' },
    { "output-rate",    required_argument, 0, 'X' },
    { "nglitch",        required_argument, 0, 'O' },
    { "earth-ephem",    required_argument, 0, 'J' },
    { "sun-ephem",      required_argument, 0, 'M' },
    { "use-cov",        no_argument,       0, '(' },
    { "covariance",     required_argument, 0, 'r' },
    { "use-priors",     no_argument,    NULL, '>' },
    { "only-joint",     no_argument,    NULL, '<' },
    { "output-burn-in", no_argument,    NULL, ')' },
    { "h0-scale",       required_argument, 0, '[' },
    { 0, 0, 0, 0 }
  };

  CHAR args[] =
"hD:p:P:i:o:a:A:j:b:B:k:s:S:m:c:C:n:l:L:q:]:Q:U:u:Y:T:v:V:z:Z:e:E:d:I:x:t:H:w:\
W:y:g:G:K:N:X:O:J:M:(r:fFR><)[:" ;
  CHAR *program = argv[0];

  /* set defaults */
  inputParams->mcmc.doMCMC = 0;/* by default don't perform an MCMC */
  inputParams->verbose = 0;    /* by default don't do verbose */
  inputParams->outputPost = 0; /* by default don't output the full posterior */
  inputParams->dob = 0.;       /* by default don't calculate an upper limit */

  inputParams->chunkMin = 5;   /* default to 5 minute minimum chunk length */
  inputParams->chunkMax = 30;  /* default to 30 minute maximum chunk length */

  /* default grid 50x50x50x50 = 6250000 points ~ 50 Mbytes */
  inputParams->mesh.minVals.h0 = 0.;     /* 0 <= h0 <= estimate */
  inputParams->mesh.maxVals.h0 = 0.;     /* estimate h0 from data by default */
  inputParams->mesh.h0Steps = 50;

  inputParams->mesh.minVals.phi0 = 0.;        /* 0 <= phi0 <= 2pi */
  inputParams->mesh.maxVals.phi0 = LAL_TWOPI;
  inputParams->mesh.phiSteps = 50;

  inputParams->mesh.minVals.psi = -LAL_PI/4.; /* -pi/4 <= psi < pi/4 */
  inputParams->mesh.maxVals.psi = LAL_PI/4.;
  inputParams->mesh.psiSteps = 50;

  inputParams->mesh.minVals.ci = -1.;         /* -1 <= cos(iota) < 1 */
  inputParams->mesh.maxVals.ci = 1.;
  inputParams->mesh.ciotaSteps = 50;

  /* default psi-time lookup table grid size 50x1440 (one per minute) */
  inputParams->mesh.psiRangeSteps = 50;
  inputParams->mesh.timeRangeSteps = 1440;

  inputParams->usepriors = 0;   /* by default don't use priors on h0, phi,
                                   iota, psi */

  /* default priors are uniform */
  inputParams->priors.h0Prior = uniform_string;
  inputParams->priors.phiPrior = uniform_string;
  inputParams->priors.psiPrior = uniform_string;
  inputParams->priors.iotaPrior = uniform_string;
  inputParams->priors.priorFile = NULL;
  inputParams->priors.h0vals = NULL;
  inputParams->priors.civals = NULL;
  inputParams->priors.h0cipdf = NULL;
  
  /* default MCMC parameters */
  inputParams->mcmc.sigmas.h0 = 0.;           /* estimate from data */
  inputParams->mcmc.sigmas.phi0 = LAL_PI_2/2.;   /* eighth of phi range */
  inputParams->mcmc.sigmas.psi = LAL_PI/16.;   /* eighth of psi range */
  inputParams->mcmc.sigmas.ci = 0.25;          /* eighth of cosi range */
  
  inputParams->mcmc.h0scale = 0.5;            /* half of default value */

  inputParams->mcmc.outputRate = 1;           /* output every sample */

  inputParams->mcmc.iterations = 10000;       /* default 10000 points */
  inputParams->mcmc.temperature = 1.;         /* default annealing */
  inputParams->mcmc.burnIn = 1000;            /* default burn in time */
  inputParams->mcmc.outputBI = 0;             /* output the burn in chain - default to no */

  inputParams->mcmc.nGlitches = 0;            /* no glitches is default */
  
  inputParams->usecov = 0;
  inputParams->matrixFile = NULL;             /* no covariance file */

  inputParams->onlyjoint = 0;       /* by default output all posteriors */

  /* parse input arguments */
  while( 1 ){
    INT4 option_index = 0;
    INT4 c;

    c = getopt_long_only( argc, argv, args, long_options, &option_index );
    if( c == -1 ) /* end of options */
      break;

    switch( c ){
      case 0:
        if( long_options[option_index].flag )
          break;
        else
          fprintf(stderr, "Error passing option %s with argument %s\n",
            long_options[option_index].name, optarg);
      case 'h': /* help message */
        printf("Usage: %s [options]\n\n%s%s\n", program, USAGE1, USAGE2);
        exit(0);
      case 'R': /* verbose */
        inputParams->verbose = 1;
        break;
      case 'f': /* output posterior distribution */
        inputParams->outputPost = 1;
        break;
      case 'F': /* perform MCMC */
        inputParams->mcmc.doMCMC = 1;
        break;
      case 'D': /* detectors */
        sprintf(inputParams->detectors, "%s", optarg);
        break;
      case 'p': /* pulsar name */
        sprintf(inputParams->pulsar, "%s", optarg);
        break;
      case 'P': /* pulsar parameter file */
        sprintf(inputParams->parFile, "%s", optarg);
        break;
      case 'i': /* input data file directory */
        sprintf(inputParams->inputDir, "%s", optarg);
        break;
      case 'o': /* output directory */
        sprintf(inputParams->outputDir, "%s", optarg);
        break;
      case 'a': /* minimum of h0 range */
        inputParams->mesh.minVals.h0 = atof(optarg);
        break;
      case 'A': /* maxmimum of h0 range */
        inputParams->mesh.maxVals.h0 = atof(optarg);
        break;
      case 'j': /* number of grid steps in h0 space */
        inputParams->mesh.h0Steps = atoi(optarg);
        break;
      case 'b': /* minimum of phi0 range */
        inputParams->mesh.minVals.phi0 = atof(optarg);
        break;
      case 'B': /* maxmimum of phi0 range */
        inputParams->mesh.maxVals.phi0 = atof(optarg);
        break;
      case 'k': /* number of grid steps in phi0 space */
        inputParams->mesh.phiSteps = atoi(optarg);
        break;
      case 's': /* minimum of psi range */
        inputParams->mesh.minVals.psi = atof(optarg);
        break;
      case 'S': /* maximum of psi range */
        inputParams->mesh.maxVals.psi = atof(optarg);
        break;
      case 'm': /* number of grid steps in psi space */
        inputParams->mesh.psiSteps = atoi(optarg);
        break;
      case 'c': /* minimum of cos(iota) range */
        inputParams->mesh.minVals.ci = atof(optarg);
        break;
      case 'C': /* maximum of cos(iota) range */
        inputParams->mesh.maxVals.ci = atof(optarg);
        break;
      case 'n': /* number of grid steps in cos(iota) space */
        inputParams->mesh.ciotaSteps = atoi(optarg);
        break;
      case 'l': /* number of bins in psi for the lookup table */
        inputParams->mesh.psiRangeSteps = atoi(optarg);
        break;
      case 'L': /* number of bins in time for the lookup table */
        inputParams->mesh.timeRangeSteps = atoi(optarg);
        break;
      case 'q': /* prior on h0 */
        inputParams->priors.h0Prior = optarg;
        break;
      case ']': /* prior file for h0 */
        inputParams->priors.priorFile = XLALStringDuplicate( optarg );
        inputParams->priors.h0Prior = optarg;
        break;
      case 'Q': /* prior on phi0 */
        inputParams->priors.phiPrior = optarg;
        break;
      case 'U': /* prior on psi */
        inputParams->priors.psiPrior = optarg;
        break;
      case 'u': /* prior on iota */
        inputParams->priors.iotaPrior = optarg;
        break;
      case 'Y': /* mean of h0 Gaussian prior */
        inputParams->priors.meanh0 = atof(optarg);
        break;
      case 'T': /* standard deviation of Gaussian h0 prior */
        inputParams->priors.stdh0 = atof(optarg);
        break;
      case 'v': /* mean of phi0 Gaussian prior */
        inputParams->priors.meanphi = atof(optarg);
        break;
      case 'V': /* standard deviation of Gaussian phi0 prior */
        inputParams->priors.stdphi = atof(optarg);
        break;
      case 'z': /* mean of psi Gaussian prior */
        inputParams->priors.meanpsi = atof(optarg);
        break;
      case 'Z': /* standard deviation of Gaussian psi prior */
        inputParams->priors.stdpsi = atof(optarg);
        break;
      case 'e': /* mean of iota Gaussian prior */
        inputParams->priors.meaniota = atof(optarg);
        break;
      case 'E': /* standard deviation of Gaussian iota prior */
        inputParams->priors.stdiota = atof(optarg);
        break;
      case 'd': /* percentage degree-of-belief at which to get UL */
        inputParams->dob = atof(optarg);
        break;
      case 'I': /* number of iterations in the MCMC chain */
        inputParams->mcmc.iterations = atoi(optarg);
        break;
      case 'x': /* number of point in the burn in stage */
        inputParams->mcmc.burnIn = atoi(optarg);
        break;
      case 't': /* temperature of the simulated annealing during burn in */
        inputParams->mcmc.temperature = atof(optarg);
        break;
      case 'H': /* standard deviation of h0 proposal distribution */
        inputParams->mcmc.sigmas.h0 = atof(optarg);
        break;
      case 'w': /* standard deviation of psi proposal distribution */
        inputParams->mcmc.sigmas.psi = atof(optarg);
        break;
      case 'W': /* standard deviation of phi0 proposal distribution */
        inputParams->mcmc.sigmas.phi0 = atof(optarg);
        break;
      case 'y': /* standard deviation of cos(iota) proposal distribution */
        inputParams->mcmc.sigmas.ci = atof(optarg);
        break;
      case 'g':
        sprintf(inputParams->mcmc.glitchTimes, "%s", optarg);
        break;
      case 'G':
        inputParams->mcmc.glitchCut = atof(optarg);
        break;
      case 'K':
        inputParams->chunkMin = atoi(optarg);
        break;
      case 'N':
        inputParams->chunkMax = atoi(optarg);
        break;
      case 'X':
        inputParams->mcmc.outputRate = atoi(optarg);
        break;
      case 'O':
        inputParams->mcmc.nGlitches = atoi(optarg);
        break;
      case 'J':
        sprintf(inputParams->earthfile, "%s", optarg);
        break;
      case 'M':
        sprintf(inputParams->sunfile, "%s", optarg);
        break;
      case '(':
        inputParams->usecov = 1;
        break;
      case 'r':
        inputParams->matrixFile = optarg;
        break;
      case '>': /* use priors on parameters */
        inputParams->usepriors = 1;
        break;
      case '<': /* only calculate/output the joint posterior for MCMC */
        inputParams->onlyjoint = 1;
        break;
      case ')': /* output the burn in chain as well as the full chain */
        inputParams->mcmc.outputBI = 1;
        break;
      case '[': /* scale factor for stdev of h0 proposal distribution */
        inputParams->mcmc.h0scale = atof(optarg);
        break;
      case '?':
        fprintf(stderr, "Unknown error while parsing options\n");
      default:
        fprintf(stderr, "Unknown error while parsing options\n");
    }
  }

  /* check parameters for weird values */
  if( inputParams->mesh.minVals.h0 < 0. || inputParams->mesh.maxVals.h0 < 0. ||
      inputParams->mesh.maxVals.h0 < inputParams->mesh.minVals.h0 ){
    fprintf(stderr, "Error... h0 grid range is wrong!\n");
    exit(0);
  }

  if( inputParams->mesh.minVals.phi0 < 0. || inputParams->mesh.maxVals.phi0 >
      LAL_TWOPI || inputParams->mesh.maxVals.phi0 <
      inputParams->mesh.minVals.phi0 ){
    fprintf(stderr, "Error... phi0 grid range is wrong!\n");
    exit(0);
  }

  if( inputParams->mesh.minVals.psi < -LAL_PI/4. ||
      inputParams->mesh.maxVals.psi >  LAL_PI/4. ||
      inputParams->mesh.maxVals.psi < inputParams->mesh.minVals.psi ){
    fprintf(stderr, "Error... psi grid range is wrong!\n");
    exit(0);
  }

  if( inputParams->mesh.minVals.ci < -1. || inputParams->mesh.maxVals.ci > 1. ||
      inputParams->mesh.maxVals.ci < inputParams->mesh.minVals.ci ){
    fprintf(stderr, "Error... cos(iota) grid range is wrong!\n");
    exit(0);
  }

  if( inputParams->dob != 0. && inputParams->mcmc.doMCMC == 1 ){
    fprintf(stderr, "Error... can't output an upper limit in MCMC mode!\n");
    fprintf(stderr, "\tNon-fatal error so continue anyway...\n");
    inputParams->dob = 0.;
  }

  if( inputParams->outputPost == 1 && inputParams->mcmc.doMCMC == 1 ){
    fprintf(stderr, "Error... can't output full posterior in MCMC mode!\n");
    fprintf(stderr, "\tNon-fatal error so continue anyway...\n");
    inputParams->outputPost = 0;
  }

  /* set mesh step sizes */
  if( inputParams->mesh.h0Steps > 1 ){
    inputParams->mesh.delta.h0 = (inputParams->mesh.maxVals.h0 -
      inputParams->mesh.minVals.h0)/(REAL8)(inputParams->mesh.h0Steps - 1.);
  }
  else
    inputParams->mesh.delta.h0 = 1.;

  if( inputParams->mesh.phiSteps > 1 ){
    inputParams->mesh.delta.phi0 = (inputParams->mesh.maxVals.phi0 -
      inputParams->mesh.minVals.phi0)/(REAL8)(inputParams->mesh.phiSteps - 1.);
  }
  else
    inputParams->mesh.delta.phi0 = 1.;

  if( inputParams->mesh.psiSteps > 1 ){
    inputParams->mesh.delta.psi = (inputParams->mesh.maxVals.psi -
      inputParams->mesh.minVals.psi)/(REAL8)(inputParams->mesh.psiSteps - 1.);
  }
  else
    inputParams->mesh.delta.psi = 1.;

  if( inputParams->mesh.ciotaSteps > 1 ){
    inputParams->mesh.delta.ci = (inputParams->mesh.maxVals.ci -
      inputParams->mesh.minVals.ci)/(REAL8)(inputParams->mesh.ciotaSteps - 1.);
  }
  else
      inputParams->mesh.delta.ci = 1.;

  if( inputParams->chunkMin < 1 || inputParams->chunkMax <
      inputParams->chunkMin ){
    fprintf(stderr, "Error... data chunk lengths are wrong!\n");
    exit(0);
  }
  
  /* read in h0 prior file if required */
  if( !inputParams->usepriors && inputParams->priors.priorFile ){
    fprintf(stderr, "Error... if h0 prior file is given then priors should\
 be used!\n");
    exit(0);
  }
  else if( inputParams->priors.priorFile ){ /* read in h0 prior file */
    FILE *fp = NULL;
    UINT4 i = 0;
    
    /* set iotaPrior to priorfile */
    inputParams->priors.iotaPrior = inputParams->priors.priorFile;
    
    if( (fp = fopen(inputParams->priors.priorFile, "rb")) == NULL ){
      fprintf(stderr, "Error... could not open prior file %s\n",
              inputParams->priors.priorFile);
      exit(0);
    }
    
    /* file should contain a header of six doubles with:
     *  - h0 minimum
     *  - dh0 - the step size in h0
     *  - N - number of h0 values
     *  - cos(iota) minimum
     *  - dci - the step size in cos(iota)
     *  - M - number of cos(iota) values
     * followed by an NxM double array of posterior values. */
   
    double header[6];
    
    /* read in header data */
    if ( !fread(header, sizeof(double), 6, fp) ){
      fprintf(stderr, "Error... could not read prior file header %s\n",
              inputParams->priors.priorFile);
      exit(0);
    }
    
    /* allocate h0 and cos(iota) vectors */
    inputParams->priors.h0vals = XLALCreateREAL8Vector( (UINT4)header[2] );
    for( i = 0; i < (UINT4)header[2]; i++ )
      inputParams->priors.h0vals->data[i] = header[0] + (REAL8)i*header[1];
    
    inputParams->priors.civals = XLALCreateREAL8Vector( (UINT4)header[5] );
    for( i = 0; i < (UINT4)header[5]; i++ )
      inputParams->priors.civals->data[i] = header[3] + (REAL8)i*header[4];
    
    /* read in prior NxM array */
    inputParams->priors.h0cipdf = XLALMalloc((UINT4)header[2]*sizeof(double*));
    for( i = 0; i < (UINT4)header[2]; i++ ){
      inputParams->priors.h0cipdf[i] = XLALMalloc( (UINT4)header[5] *
                                                   sizeof(double) );
   
      if ( !fread(inputParams->priors.h0cipdf[i], sizeof(double),
                 (UINT4)(header[5]), fp ) ){
        fprintf(stderr, "Error... could not read prior file array %s\n",
                inputParams->priors.priorFile);
        exit(0);
      }
    }
    
    fclose(fp);
  } 
}



/* function to allocate memory for a likelihood array - it will be index as
   logLike[phi][cosiota][psi][h0] */
REAL8 ****allocate_likelihood_memory(MeshGrid mesh){
  INT4 i=0, j=0, k=0;
  REAL8 ****logLike=NULL;

  /* allocate the h0 positions using calloc (i.e. array will be initialise to 
     zero */
  logLike = XLALCalloc(mesh.phiSteps, sizeof(REAL8 ***));

  for( i = 0 ; i < mesh.phiSteps ; i++ ){
    logLike[i] = XLALCalloc(mesh.ciotaSteps, sizeof(REAL8 **));

    for( j = 0 ; j < mesh.ciotaSteps ; j++ ){
      logLike[i][j] = XLALCalloc(mesh.psiSteps, sizeof(REAL8 *));

      for( k = 0 ; k < mesh.psiSteps ; k++ ){
        logLike[i][j][k] = XLALCalloc(mesh.h0Steps, sizeof(REAL8));
      }
    }
  }

  return logLike;
}



/* function to create a log likelihood array over the parameter grid */ 
REAL8 create_likelihood_grid(DataStructure data, REAL8 ****logLike, 
  MeshGrid mesh){
  IntrinsicPulsarVariables vars;

  INT4 i=0, j=0, k=0;

  REAL4 cosphi=0., sinphi=0.;

  REAL8 noiseEvidence=0.;

  /* create vector of data segment length */
  data.chunkLengths = NULL;
  data.chunkLengths = XLALCreateINT4Vector(
    (INT4)ceil(data.data->length/data.chunkMin) );

  get_chunk_lengths(data);

  /* allocate memory for summations */
  data.sumData = NULL;
  data.sumData = XLALCreateREAL8Vector(data.chunkLengths->length);

  /* get the sum over the data */
  sum_data(data);
  
  /* calculate likelihood array */
  for( i = 0 ; i < mesh.phiSteps ; i++ ){
    if( verbose )
      fprintf(stderr, "In phi0 loop %d of %d.\n", i+1, mesh.phiSteps);

    vars.phi0 = mesh.minVals.phi0 + (REAL8)i*mesh.delta.phi0;
    
    sin_cos_LUT( &sinphi, &cosphi, vars.phi0 );

    for( j = 0 ; j < mesh.ciotaSteps ; j++ ){
      vars.ci = mesh.minVals.ci + (REAL8)j*mesh.delta.ci;
      
      vars.Xplus = 0.5*(1.+vars.ci*vars.ci);
      vars.Xcross = vars.ci;
      vars.Xpsinphi_2 = 0.5*vars.Xplus*sinphi;
      vars.Xcsinphi_2 = 0.5*vars.Xcross*sinphi;
      vars.Xpcosphi_2 = 0.5*vars.Xplus*cosphi;
      vars.Xccosphi_2 = 0.5*vars.Xcross*cosphi;

      for( k = 0 ; k < mesh.psiSteps ; k++ ){
        vars.psi = mesh.minVals.psi + (REAL8)k*mesh.delta.psi;
        
        /* perform final loop over h0 within log_likelihood function */
        noiseEvidence = log_likelihood(logLike[i][j][k], data, vars, mesh,
          NULL);
      }
    }
  }

  /* free memory */
  XLALDestroyREAL8Vector(data.sumData);

  return noiseEvidence;
}



/* a function to sum over the data */
void sum_data(DataStructure data){
  INT4 chunkLength=0, length=0, i=0, j=0, count=0;
  COMPLEX16 B;

  length = data.data->length + 1 -
           data.chunkLengths->data[data.chunkLengths->length-1];

  for( i = 0 ; i < length ; i+= chunkLength ){
    chunkLength = data.chunkLengths->data[count];
    data.sumData->data[count] = 0.;

    for( j = i ; j < i + chunkLength ; j++){
      B.re = data.data->data[j].re;
      B.im = data.data->data[j].im;

      /* sum up the data */
      data.sumData->data[count] += (B.re*B.re + B.im*B.im);
    }

    count++;
  }
}



/* function to calculate the log(likelihood) given some data and a set of
   particular pulsar parameters */
REAL8 log_likelihood( REAL8 *likeArray, DataStructure data,
  IntrinsicPulsarVariables vars, MeshGrid mesh, REAL8Vector *dphi ){
  INT4 i=0, j=0, count=0, k=0;
  INT4 length=0;
  REAL8 chunkLength=0.;

  REAL8 tstart=0., T=0.;

  COMPLEX16 model;
  INT4 psibin=0, timebin=0;

  REAL8 plus=0., cross=0.;
  REAL8 sumModel=0., sumDataModel=0.;
  REAL8 chiSquare=0.;
  COMPLEX16 B;

  REAL8 exclamation[data.chunkMax+1]; /* all factorials up to chunkMax */
  REAL8 logOf2=log(2.);

  REAL8 noiseEvidence=0.; /* the log evidence that the data is just noise */

  INT4 first=0, through=0;

  REAL8 psteps = (REAL8)data.lookupTable->psiSteps;
  REAL8 tsteps = (REAL8)data.lookupTable->timeSteps;
  
  /* to save time get all log factorials up to chunkMax */
  for( i = 0 ; i < data.chunkMax+1 ; i++ )
    exclamation[i] = log_factorial(i);

  /* set the psi bin for the lookup table */
  psibin = (INT4)ROUND( ( vars.psi + LAL_PI/4. ) * ( psteps-1. )/LAL_PI_2 );

  length = (INT4)data.data->length + 1 -
           data.chunkLengths->data[(INT4)data.chunkLengths->length-1];

  tstart = data.times->data[0]; /* time of first B_k */
  
  for( i = 0 ; i < length ; i += chunkLength ){
    chunkLength = (REAL8)data.chunkLengths->data[count];

    if( chunkLength < data.chunkMin ){
      count++;

      if( through == 0 ) first = 0;

      continue;
    }

    through = 1;

    sumModel = 0.;
    sumDataModel = 0.;

    for( j = i ; j < i + chunkLength ; j++){
      /* set the time bin for the lookup table */
      /* sidereal day in secs*/
      T = fmod(data.times->data[j] - tstart, LAL_DAYSID_SI);
      timebin = (INT4)fmod( ROUND(T*tsteps/LAL_DAYSID_SI), tsteps );

      plus = data.lookupTable->lookupTable[psibin][timebin].plus;
      cross = data.lookupTable->lookupTable[psibin][timebin].cross;

      B.re = data.data->data[j].re;
      B.im = data.data->data[j].im;

      /*********************************************************/
      /* stuff for phase offset due to parameter uncertainties - MCMC only */
      if( dphi != NULL ){
        REAL4 cphi=0., sphi=0.;

        /* create the signal model */
        sin_cos_2PI_LUT( &sphi, &cphi, -dphi->data[j] );

        model.re = (plus*vars.Xpcosphi_2 + cross*vars.Xcsinphi_2)*cphi +
                 (cross*vars.Xccosphi_2 - plus*vars.Xpsinphi_2)*sphi;
        model.im = (plus*vars.Xpsinphi_2 - cross*vars.Xccosphi_2)*cphi +
                 (cross*vars.Xcsinphi_2 + plus*vars.Xpcosphi_2)*sphi;
      }
      /*********************************************************/
      else{
        /* create the signal model */
        model.re = plus*vars.Xpcosphi_2 + cross*vars.Xcsinphi_2;
        model.im = plus*vars.Xpsinphi_2 - cross*vars.Xccosphi_2;
      }
      
      /* sum over the model */
      sumModel += model.re*model.re + model.im*model.im;
      
      /* sum over that data and model */
      sumDataModel += B.re*model.re + B.im*model.im;
    }

    for( k = 0 ; k < mesh.h0Steps ; k++ ){
      vars.h0 = mesh.minVals.h0 + (REAL8)k*mesh.delta.h0;
      
      chiSquare = data.sumData->data[count];
      chiSquare -= 2.*vars.h0*sumDataModel;
      chiSquare += vars.h0*vars.h0*sumModel;
      
      /* log(likelihood)
         logL = (m-1)log(2) + log(m!) - m*log(sum((Bk - yk)^2)) */

      /* reset array if first time in loop - otherwise joint likelihoods
         will be wrong */
      if( first == 0 )
        likeArray[k] = (chunkLength - 1.)*logOf2;
      else likeArray[k] += (chunkLength - 1.)*logOf2;

      likeArray[k] += exclamation[(INT4)chunkLength];
      likeArray[k] -= chunkLength*log(chiSquare);
      
      /*** SET LIKELIHOOD TO CONSTANT TO CHECK PRIOR IS RETURNED PROPERLY ****/
      //likeArray[k] = 0.;
      /***********************************************************************/
    }
    
    /* get the log evidence for the data not containing a signal */
    noiseEvidence += (chunkLength - 1.)*logOf2;
    noiseEvidence += exclamation[(INT4)chunkLength];
    noiseEvidence -= chunkLength*log(data.sumData->data[count]);

    first++;
    count++;
  }

  return noiseEvidence;
}



/* function to combine log likelihoods to give a joint likelihood 
   log(p(data_joint|a) = log(p(data1|a)) + log(p(data2|a))        */
void combine_likelihoods(REAL8 ****logLike1, REAL8 ****logLike2, 
  MeshGrid mesh){
  INT4 i=0, j=0, k=0, n=0;

  for( i = 0 ; i < mesh.phiSteps ; i++ ){
    for( j = 0 ; j < mesh.ciotaSteps ; j++ ){
      for( k = 0 ; k < mesh.psiSteps ; k++ ){
        for( n = 0 ; n < mesh.h0Steps ; n++ ){
          logLike2[i][j][k][n] += logLike1[i][j][k][n];
        }
      }
    }
  }
}



/* function to calculate the log prior */
REAL8 log_prior(PriorVals prior, MeshGrid mesh){
  REAL8 pri=0.;

  if ( prior.priorFile == NULL ){
    if(strcmp(prior.h0Prior, "uniform") == 0){
      pri = 0.; /* set h0 prior to be one for all values if uniform */
      /* pri *= 1./(mesh.maxVals.h0 - mesh.minVals.h0); */
    }
    else if(strcmp(prior.h0Prior, "jeffreys") == 0) pri += -log(prior.vars.h0);
    else if(strcmp(prior.h0Prior, "gaussian") == 0){
      pri += -log(prior.stdh0*sqrt(LAL_TWOPI)) + (-(prior.vars.h0 -
        prior.meanh0)*(prior.vars.h0 -
        prior.meanh0)/(2.*prior.stdh0*prior.stdh0));
    }
  
    if(strcmp(prior.iotaPrior, "uniform") == 0){
      pri += -log(fabs(acos(mesh.maxVals.ci) - acos(mesh.minVals.ci)));
    }
    /* wrap around at zero and pi */
    else if(strcmp(prior.iotaPrior, "gaussian") == 0){
      REAL8 iota = acos(prior.vars.ci);
      if( iota < prior.meaniota - LAL_PI_2 ) iota += LAL_PI;
      else if( prior.meaniota + LAL_PI_2 < iota ) prior.vars.ci -= LAL_PI;

      pri += -log(prior.stdiota*sqrt(LAL_TWOPI) ) + ( -( iota -
              prior.meaniota ) * ( iota - prior.meaniota ) / 
             ( 2.*prior.stdiota*prior.stdiota ) );
    }
  }
  else{
    UINT4 i = 0, j = 0;
    
    /* use h0 prior read in from a file */
    REAL8 h0low = 0., h0high = 0., pdflow = 0., pdfhigh = 0.;
    REAL8 cilow = 0., cihigh = 0.;
    REAL8 pdf00 = 0., pdf01 = 0., pdf10 = 0., pdf11 = 0.;
    REAL8 grad = 0.;
    
    /* find the nearest h0 point */
    if ( prior.vars.h0 <= prior.h0vals->data[0] ) i = 0;
    else if ( prior.vars.h0 >= prior.h0vals->data[prior.h0vals->length-1] )
      i = prior.h0vals->length;
    else{
      for( i = 1; i < prior.h0vals->length; i++ ){
        if( prior.h0vals->data[i-1] < prior.vars.h0 && 
            prior.vars.h0 <= prior.h0vals->data[i] ) break;
      }
    }
    
    /* find the nearest cos(iota) point */
    if ( prior.vars.ci <= prior.civals->data[0] ) j = 0;
    else if ( prior.vars.ci >= prior.civals->data[prior.civals->length-1] )
      j = prior.civals->length;
    else{
      for( j = 1; j < prior.civals->length; j++ ){
         if( prior.civals->data[j-1] < prior.vars.ci && 
            prior.vars.ci <= prior.civals->data[j] ) break;
      }
    }
    
    if( i == 0 || i == prior.h0vals->length ){
      /* if the point is less than or greater than the edge of h0 range then
       * just linearly interpolate in cos(iota) */
      if( i == 0 && j == 0 ) pri += log( prior.h0cipdf[i][j] );
      else if( i == 0 && j == prior.civals->length )
        pri += log( prior.h0cipdf[i][j-1] );
      else if( i == prior.h0vals->length && j == 0 )
        pri += log( prior.h0cipdf[i-1][j] );
      else if( i == prior.h0vals->length && j == prior.civals->length )
        pri += log( prior.h0cipdf[i-1][j-1] );
      else{
        cilow = prior.civals->data[j-1];
        cihigh = prior.civals->data[j];
        
        if ( i == 0 ){
          pdflow = prior.h0cipdf[i][j-1];
          pdfhigh = prior.h0cipdf[i][j];
        }
        else{
          pdflow = prior.h0cipdf[i-1][j-1];
          pdfhigh = prior.h0cipdf[i-1][j];
        }
        
        grad = (pdfhigh-pdflow)/(cihigh-cilow);
        pri += log( pdflow + grad*(prior.vars.ci - cilow) );
      }
    }
    else if ( j == 0 || j == prior.civals->length ){
      /* if the point is less than or greater than the edge of ci range then
       * just linearly interpolate in h0 */
      if( j == 0 && i == 0 ) pri += log( prior.h0cipdf[i][j] );
      else if( j == 0 && i == prior.h0vals->length ) 
        pri += log( prior.h0cipdf[i-1][j] );
      else if( j == prior.civals->length && i == 0 )
        pri += log( prior.h0cipdf[i][j-1] );
      else if( j == prior.civals->length && i == prior.h0vals->length )
        pri += log( prior.h0cipdf[i-1][j-1] );
      else{
        h0low = prior.h0vals->data[i-1];
        h0high = prior.h0vals->data[i];
        
        if ( j == 0 ){
          pdflow = prior.h0cipdf[i-1][j];
          pdfhigh = prior.h0cipdf[i][j];
        }
        else{
          pdflow = prior.h0cipdf[i-1][j-1];
          pdfhigh = prior.h0cipdf[i][j-1];
        }
        
        grad = (pdfhigh-pdflow)/(h0high-h0low);
        pri += log( pdflow + grad*(prior.vars.h0 - h0low) );
      }
    }
    else{ /* bilinearly interpolate */
      h0low = prior.h0vals->data[i-1], h0high = prior.h0vals->data[i];
      cilow = prior.civals->data[j-1], cihigh = prior.civals->data[j];
      
      REAL8 h0scaled = (prior.vars.h0 - h0low)/(h0high - h0low);
      REAL8 ciscaled = (prior.vars.ci - cilow)/(cihigh - cilow);
      
      pdf00 = prior.h0cipdf[i-1][j-1];
      pdf01 = prior.h0cipdf[i-1][j];
      pdf10 = prior.h0cipdf[i][j-1];
      pdf11 = prior.h0cipdf[i][j];
            
      pri += log( pdf00*(1. - h0scaled)*(1. - ciscaled) + 
              pdf10*h0scaled*(1. - ciscaled) + pdf01*(1. - h0scaled)*ciscaled + 
              pdf11*h0scaled*ciscaled );
    }
  }

  if(strcmp(prior.phiPrior, "uniform") == 0){
    pri += -log(mesh.maxVals.phi0 - mesh.minVals.phi0);
  }
  /* wrap around Gaussian priors so that the mean is always at the centre */
  else if(strcmp(prior.phiPrior, "gaussian") == 0){
    if( prior.vars.phi0 < prior.meanphi - LAL_PI )
      prior.vars.phi0 += LAL_TWOPI;
    else if( prior.meanphi + LAL_PI < prior.vars.phi0 )
      prior.vars.phi0 -= LAL_TWOPI;

    pri += -log(prior.stdphi*sqrt(LAL_TWOPI) ) + ( -( prior.vars.phi0 -
            prior.meanphi ) * ( prior.vars.phi0 - prior.meanphi ) /
           ( 2.*prior.stdphi*prior.stdphi ) );
  }

  if(strcmp(prior.psiPrior, "uniform") == 0){
    pri += -log(mesh.maxVals.psi - mesh.minVals.psi);
  }
  else if(strcmp(prior.psiPrior, "gaussian") == 0){
    if( prior.vars.psi < prior.meanpsi - LAL_PI_4 )
      prior.vars.psi += LAL_PI_2;
    else if( prior.meanpsi + LAL_PI_4 < prior.vars.psi )
      prior.vars.psi -= LAL_PI_2;

    pri += -log(prior.stdpsi*sqrt(LAL_TWOPI) ) + ( -( prior.vars.psi -
            prior.meanpsi ) * ( prior.vars.psi - prior.meanpsi ) / 
           ( 2.*prior.stdpsi*prior.stdpsi ) );
  }

  return pri;
}



/* function to calculate the unnormalised log posterior and output the max value 
   - print out the log posterior if requested */
REAL8 log_posterior(REAL8 ****logLike, PriorVals prior, MeshGrid mesh,
  OutputParams output){
  REAL8 maxPost=-INFINITY, mP=0.;
  REAL8 logPi=0.;

  INT4 i=0, j=0, k=0, n=0;

  CHAR outputFile[256];
  FILE *fp=NULL;

  /* if outputPost != NULL then print out the log Posterior to a file */
  if( output.outPost == 1 ){
    sprintf(outputFile, "%s/log_posterior.%s.%s", output.outputDir, output.psr,
      output.det);

    if( ( fp = fopen(outputFile, "w") ) == NULL ){
      fprintf(stderr, "Error... Can't output posterior file!");
    }
  }

  /* create posterior */
  for( i = 0 ; i < mesh.phiSteps ; i++ ){
    prior.vars.phi0 = mesh.minVals.phi0 + (REAL8)i*mesh.delta.phi0;

    for( j = 0 ; j < mesh.ciotaSteps ; j++ ){
      prior.vars.ci = mesh.minVals.ci + (REAL8)j*mesh.delta.ci;

      for( k = 0 ; k < mesh.psiSteps ; k++ ){
        prior.vars.psi = mesh.minVals.psi + (REAL8)k*mesh.delta.psi;

        for( n = 0 ; n < mesh.h0Steps ; n++ ){
          prior.vars.h0 = mesh.minVals.h0 + (REAL8)n*mesh.delta.h0;

          logPi = log_prior(prior, mesh);

          /* posterior \propto log(likelihood) + log(prior) */
          logLike[i][j][k][n] += logPi;
          mP = logLike[i][j][k][n];

          if( mP > maxPost )
            maxPost = mP;

          if( output.outPost == 1 && fp != NULL )
            fprintf(fp, "%le\n", mP);
        }
      }
    }
  }

  return maxPost;
}



/* function to marginalise posterior over requested parameter and output the log 
   evidence if requested */
Results marginalise_posterior(REAL8 ****logPost, MeshGrid mesh, 
  OutputParams output){
  REAL8 dval1=0., dval2=0., dval3=0., dval4=0.;

  REAL8 post1=0., post2=0.;

  REAL8 ***evSum1=NULL;      /* first integral */
  REAL8 **evSum2=NULL;       /* second integral */
  REAL8 *evSum3=NULL;        /* third integral */
  REAL8 *cumsum=NULL;        /* cumulative probability */
  REAL8 evSum4=-INFINITY;    /* fouth integral */
  REAL8 maxPost=-INFINITY;
  REAL8 evVal=0., sumVal=0.;

  Results results={0.,0.,0.,0.}; /* the results structure */

  INT4 numSteps1=0, numSteps2=0, numSteps3=0, numSteps4=0;
  REAL8 step=0., minVal=0.;

  INT4 i=0, j=0, k=0, n=0;

  CHAR outputFile[256];  /* filname to output the marginalised pdf */
  FILE *fp=NULL;

  /* check which parameter we will not be integrating over for the output */
  if( strcmp( output.margParam, "h0" ) == 0 ){
    numSteps1 = mesh.h0Steps;
    numSteps2 = mesh.phiSteps;
    numSteps3 = mesh.ciotaSteps;
    numSteps4 = mesh.psiSteps;

    dval1 = mesh.delta.h0; dval2 = mesh.delta.phi0; dval3 = mesh.delta.ci;
    dval4 = mesh.delta.psi;
    minVal = mesh.minVals.h0;
  }
  else if( strcmp( output.margParam, "phi" ) == 0 ){
    numSteps1 = mesh.phiSteps;
    numSteps2 = mesh.ciotaSteps;
    numSteps3 = mesh.psiSteps;
    numSteps4 = mesh.h0Steps;

    dval1 = mesh.delta.phi0; dval2 = mesh.delta.ci; dval3 = mesh.delta.psi;
    dval4 = mesh.delta.h0;
    minVal = mesh.minVals.phi0;
  }
  else if( strcmp( output.margParam, "psi" ) == 0 ){
    numSteps1 = mesh.psiSteps;
    numSteps2 = mesh.phiSteps;
    numSteps3 = mesh.ciotaSteps;
    numSteps4 = mesh.h0Steps;

    dval1 = mesh.delta.psi; dval2 = mesh.delta.phi0; dval3 = mesh.delta.ci;
    dval4 = mesh.delta.h0;
    minVal = mesh.minVals.psi;
  }
  else if( strcmp( output.margParam, "ciota" ) == 0 ){
    numSteps1 = mesh.ciotaSteps;
    numSteps2 = mesh.phiSteps;
    numSteps3 = mesh.psiSteps;
    numSteps4 = mesh.h0Steps;

    dval1 = mesh.delta.ci; dval2 = mesh.delta.phi0; dval3 = mesh.delta.psi;
    dval4 = mesh.delta.h0;
    minVal = mesh.minVals.ci;
  }
  else{
    fprintf(stderr, "Error... No parameter selected for posterior output!\n");
    exit(0);
  }

  /* allocate memory for integration */
  evSum1 = XLALCalloc(numSteps1, sizeof(REAL8 **));

  /* perform first integral */
  for( i = 0 ; i < numSteps1 ; i++ ){
    evSum1[i] = XLALCalloc(numSteps2, sizeof(REAL8 *));
    for( j = 0 ; j < numSteps2 ; j++ ){
      evSum1[i][j] = XLALCalloc(numSteps3, sizeof(REAL8));
      for( k = 0 ; k < numSteps3 ; k++ ){
        evSum1[i][j][k] = -INFINITY; /* initialise */

        /* if we only have one point in the parameter space */
        if( numSteps4 == 1 ){
          if( strcmp( output.margParam, "h0" ) == 0 ) 
            evSum1[i][j][k] = logPost[j][k][0][i];
          else if( strcmp( output.margParam, "phi" ) == 0 )
            evSum1[i][j][k] = logPost[i][j][k][0];
          else if( strcmp( output.margParam, "psi" ) == 0 )
            evSum1[i][j][k] = logPost[j][k][i][0];
          else if( strcmp( output.margParam, "ciota" ) == 0 )
            evSum1[i][j][k] = logPost[j][i][k][0];
        }
        else{
          for( n = 0 ; n < numSteps4 - 1 ; n++ ){
            if( strcmp( output.margParam, "h0" ) == 0 ){
              post1 = logPost[j][k][n][i];
              post2 = logPost[j][k][n+1][i];

              /* use trapezium rule for intergration */
              evVal = evSum1[i][j][k];
              sumVal = LOG_TRAPEZIUM(post1, post2, dval4);
              evSum1[i][j][k] = PLUS(sumVal, evVal);
            }
            else if( strcmp( output.margParam, "phi" ) == 0 ){
              post1 = logPost[i][j][k][n];
              post2 = logPost[i][j][k][n+1];

              evVal = evSum1[i][j][k];
              sumVal = LOG_TRAPEZIUM(post1, post2, dval4);
              evSum1[i][j][k] = PLUS(sumVal, evVal);
            }
            else if( strcmp( output.margParam, "psi" ) == 0 ){
              post1 = logPost[j][k][i][n];
              post2 = logPost[j][k][i][n+1];

              evVal = evSum1[i][j][k];
              sumVal = LOG_TRAPEZIUM(post1, post2, dval4);
              evSum1[i][j][k] = PLUS(sumVal, evVal);
            }
            else if( strcmp( output.margParam, "ciota" ) == 0 ){
              post1 = logPost[j][i][k][n];
              post2 = logPost[j][i][k][n+1];

              evVal = evSum1[i][j][k];
              sumVal = LOG_TRAPEZIUM(post1, post2, dval4);
              evSum1[i][j][k] = PLUS(sumVal, evVal);
            }
          }
/*           for( n = 0 ; n < numSteps4 - 1 ; n++ ){
             if( strcmp( output.margParam, "h0" ) == 0 ){
              post1 = logPost[j][k][n][i];

              evVal = evSum1[i][j][k];
              sumVal = post1 + log(dval4);
              evSum1[i][j][k] = PLUS(sumVal, evVal);
            }
            else if( strcmp( output.margParam, "phi" ) == 0 ){
              post1 = logPost[i][j][k][n];

              evVal = evSum1[i][j][k];
              sumVal = post1 + log(dval4);
              evSum1[i][j][k] = PLUS(sumVal, evVal);
            }
            else if( strcmp( output.margParam, "psi" ) == 0 ){
              post1 = logPost[j][k][i][n];

              evVal = evSum1[i][j][k];
              sumVal = post1 + log(dval4);
              evSum1[i][j][k] = PLUS(sumVal, evVal);
            }
            else if( strcmp( output.margParam, "ciota" ) == 0 ){
              post1 = logPost[j][i][k][n];

              evVal = evSum1[i][j][k];
              sumVal = post1 + log(dval4);
              evSum1[i][j][k] = PLUS(sumVal, evVal);
            }
          } */
        }
      }
    }
  }

  /* allocate memory for second integral */
  evSum2 = XLALCalloc(numSteps1, sizeof(REAL8 *));

  /* perform the second integration */
  for( i = 0 ; i < numSteps1 ; i++ ){
    evSum2[i] = XLALCalloc(numSteps2, sizeof(REAL8));
    for( j = 0 ; j < numSteps2 ; j++ ){
      evSum2[i][j] = -INFINITY;

      if( numSteps3 == 1 ) evSum2[i][j] = evSum1[i][j][0];
      else{  
        for( k = 0 ; k <numSteps3 - 1 ; k++ ){
          evVal = evSum2[i][j];
          sumVal = LOG_TRAPEZIUM(evSum1[i][j][k], evSum1[i][j][k+1], dval3);
          evSum2[i][j] = PLUS(sumVal, evVal);
        }
/*        for( k = 0 ; k < numSteps3 - 1 ; k++ ){
              evVal = evSum2[i][j];
              sumVal = evSum1[i][j][k] + log(dval3);
              evSum2[i][j] = PLUS(sumVal, evVal);
          } */
      }
    }
  }

  /* allocate memory for third integration */
  evSum3 = XLALCalloc(numSteps1, sizeof(REAL8));

  /* perform third integration */
  for( i = 0 ; i < numSteps1 ; i++ ){
    evSum3[i] = -INFINITY;

    if( numSteps2 == 1 ) evSum3[i] = evSum2[i][0];
    else{
      for( j = 0 ; j < numSteps2 - 1 ; j++ ){
        evVal = evSum3[i];
        sumVal = LOG_TRAPEZIUM(evSum2[i][j], evSum2[i][j+1], dval2);
        evSum3[i] = PLUS(sumVal, evVal);
      }
/*       for( j = 0 ; j < numSteps2 - 1 ; j++ ){
           evVal = evSum3[i];
           sumVal = evSum2[i][j] + log(dval2);
           evSum3[i] = PLUS(sumVal, evVal);
         } */
    }

    if( evSum3[i] > maxPost )
      maxPost = evSum3[i];
  }

  if( numSteps1 == 1 ) evSum4 = evSum3[0];
  else{
    /* perform final integration to get the log evidence */
    for( i = 0 ; i < numSteps1 - 1 ; i++ ){
      evVal = evSum4;
      sumVal = LOG_TRAPEZIUM(evSum3[i], evSum3[i+1], dval1);
      evSum4 = PLUS(sumVal, evVal);
    }
/*     for( i = 0 ; i < numSteps1 - 1 ; i++ ){
         evVal = evSum4;
         sumVal = evSum3[i] + log(dval1);
         evSum4 = PLUS(sumVal, evVal);
       } */
  }

  /* create the file to output the marginalise posterior to */
  sprintf(outputFile, "%s/pdf_%s.%s.%s", output.outputDir, output.margParam,
    output.psr, output.det);
  if( ( fp = fopen(outputFile, "w") ) == NULL ){
    fprintf(stderr, "Error... Can't open file to output pdf for %s.\n!",
      output.margParam);
    exit(0);
  }

  /* allocate memory for cumulative sum (needed when we want to produce an UL */
  cumsum = XLALCalloc(numSteps1, sizeof(REAL8));

  /* normalise the marginalised pdf using the evidence and free memory */
  for( i = 0 ; i < numSteps1 ; i++ ){
    for( j = 0 ; j < numSteps2 ; j++ ){
      XLALFree(evSum1[i][j]);
    }

    /* the parameter value at a given point */
    step = minVal + (REAL8)i*dval1;

    /* calculate the cumulative probability */
    if( numSteps1 == 1 ) cumsum[i] = exp(evSum3[0] - evSum4);
    else{
      if( i == 0 ){
        cumsum[i] = 0.;
      }
      else{
        cumsum[i] = cumsum[i-1] + exp(LOG_TRAPEZIUM(evSum3[i-1], evSum3[i],
dval1) - evSum4);
/*         cumsum[i] = cumsum[i-1] + exp(evSum3[i-1] - evSum4) * dval1; */
      }
    }
    /* print out marginalised posterior */
    fprintf(fp, "%le\t%le\t%le\n", step, exp(evSum3[i] - evSum4), cumsum[i]);

    XLALFree(evSum1[i]);
    XLALFree(evSum2[i]);
  }

  results.evidence = evSum4;

  /* get the h0 upper limit if required */
  if( strcmp( output.margParam, "h0" ) == 0 && output.dob != 0 )
    results.h0UpperLimit = get_upper_limit(cumsum, output.dob, mesh); 

  fclose(fp);

  XLALFree(evSum1);
  XLALFree(evSum2);
  XLALFree(evSum3);
  XLALFree(cumsum);

  return results; /* return the log evidence */
}



/* detector response lookup table function  - this function will output a lookup
table of points in time and psi, covering a sidereal day from the start time
(t0) and from -pi/4 to pi/4 in psi */
void response_lookup_table(REAL8 t0, LALDetAndSource detAndSource,
  DetRespLookupTable *lookupTable){ 
  LIGOTimeGPS gps;
  REAL8 T=0;

  REAL8 fplus=0., fcross=0.;
  REAL8 psteps = (REAL8)lookupTable->psiSteps;
  REAL8 tsteps = (REAL8)lookupTable->timeSteps;

  INT4 i=0, j=0;

  for( i = 0 ; i < lookupTable->psiSteps ; i++ ){
    detAndSource.pSource->orientation = -(LAL_PI/4.) +
        (REAL8)i*(LAL_PI/2.) / ( psteps - 1. );

    for( j = 0 ; j < lookupTable->timeSteps ; j++ ){
      T = t0 + (REAL8)j*LAL_DAYSID_SI / tsteps;

      gps.gpsSeconds = (INT4)floor(T);
      gps.gpsNanoSeconds = (INT4)floor((fmod(T,1.0)*1.e9));

      XLALComputeDetAMResponse(&fplus, &fcross,
        detAndSource.pDetector->response,
        detAndSource.pSource->equatorialCoords.longitude,
        detAndSource.pSource->equatorialCoords.latitude,
        detAndSource.pSource->orientation,
        XLALGreenwichMeanSiderealTime(&gps));

      lookupTable->lookupTable[i][j].plus = fplus;
      lookupTable->lookupTable[i][j].cross = fcross;
    }
  }
}



/* function to return the (REAL8) log factorial of an integer */
REAL8 log_factorial(INT4 num){
  INT4 i=0;
  REAL8 logFac=0.;

  for( i=2 ; i <= num ; i++ ) logFac += log((REAL8)i);

  return logFac;
}



/* function to calculate the upper limit - use quadratic spline interpolation 
  between points around the upper limit */
REAL8 get_upper_limit(REAL8 *cumsum, REAL8 limit, MeshGrid mesh){
  REAL8 ul1=0., ul2=0.;
  INT4 i=0;

  INT4 point1=0, point2=0;
  REAL8 vals[3]={0., 0., 0.};
  REAL8 h0s[3]={0., 0., 0.};  
  REAL8 x[3]={0., 0., 0.}, y[3]={0., 0., 0.};

  REAL8 det=0., c[3]={0., 0., 0.};

  /* find the point above and below the required upper limit */
  for( i = 0 ; i < mesh.h0Steps ; i++ ){
    if( cumsum[i] < limit/100. ){
      vals[0] = cumsum[i];
      point1 = i;
      continue;
    }
    else{
      vals[1] = cumsum[i];
      point2 = i;
      break;
    }
  }

  /* couldn't calculate an upper limit! */
  if( point2 == 0 || point1 < 1 || point2 > (mesh.h0Steps - 2) ){
    fprintf(stderr, "Couldn't calculate a %.1lf%% upper limit!\n", limit);
  }

  h0s[0] = mesh.minVals.h0 + (REAL8)(point1*mesh.delta.h0);
  h0s[1] = mesh.minVals.h0 + (REAL8)(point2*mesh.delta.h0);

  /* find out if point one or point two is closer to UL */
  if( fabs(cumsum[point1] - limit/100.) < fabs(cumsum[point2] - limit/100.) ){
    vals[2] = cumsum[point1-1];
    h0s[2] = h0s[0] - mesh.delta.h0;
    x[0] = h0s[2]; x[1] = h0s[0]; x[2] = h0s[1];
    y[0] = vals[2]; y[1] = vals[0]; y[2] = vals[1];
  }
  else{
    vals[2] = cumsum[point2+1];
    h0s[2] = h0s[1] + mesh.delta.h0;
    x[0] = h0s[0]; x[1] = h0s[1]; x[2] = h0s[2];
    y[0] = vals[0]; y[1] = vals[1]; y[2] = vals[2];
  }

  /* calculate determinant and quadratic coefficients */
  det = x[0]*x[0]*(x[1]-x[2])-x[0]*(x[1]*x[1]-x[2]*x[2])+(x[1]*x[1]*x[2]
        - x[1]*x[2]*x[2]);
  c[0] = ((x[1]-x[2])*y[0]-(x[0]-x[2])*y[1]+(x[0]-x[1])*y[2])/det;
  c[1] = (-(x[1]*x[1]-x[2]*x[2])*y[0] + (x[0]*x[0]-x[2]*x[2])*y[1] -
           (x[0]*x[0]-x[1]*x[1])*y[2])/det;
  c[2] = ( (x[1]*x[1]*x[2]-x[2]*x[2]*x[1])*y[0] -
           (x[0]*x[0]*x[2]-x[0]*x[2]*x[2])*y[1] +
           (x[0]*x[0]*x[1]-x[0]*x[1]*x[1])*y[2] )/det;

  /* put equation in the form ax^2 + bx + c = 0 */
  c[2] -= limit/100.;

  /* use quadratic formula to calculate upper limit */
  ul1 = (-c[1] + sqrt(c[1]*c[1] - 4.*c[0]*c[2]))/(2.*c[0]);
  ul2 = (-c[1] - sqrt(c[1]*c[1] - 4.*c[0]*c[2]))/(2.*c[0]);

  /* work out derivative of the gradient to see which solution to use */
  if( (y[2]-y[1]) < (y[1]-y[0]) ){ /* quadratic is negative /^\*/
    /* return the smaller of the two solutions */
    if( ul1 < ul2 ) return ul1;
    else return ul2;
  }
  else{ /* quadratic is positive  \_/ */
    /* return the larger of the two solutions */
    if( ul1 > ul2 ) return ul1;
    else return ul2;
  }
}



/* function to perform the MCMC parameter estimation */
void perform_mcmc(DataStructure *data, InputParams input, INT4 numDets, 
  CHAR *det, LALDetector *detPos, EphemerisData *edat ){
  static LALStatus status;

  IntrinsicPulsarVariables vars, varsNew;

  /* extra h0 and phase parameter required if a glitch is obsevered */
  IntrinsicPulsarVariables *extraVars=NULL, *extraVarsNew=NULL;
  BinaryPulsarParams pulsarParamsFixed;
  REAL8 *glitchTimes=NULL;
  CHAR *gtimestr=NULL;
  INT4 nGlitches=input.mcmc.nGlitches;
  INT4 **g1=NULL, **g2=NULL; /* start - end positions of each glitch segment */
  DataStructure **glitchData=NULL;

  /* INT4 below0=0; */

  REAL4Vector *randNum=NULL; /* LAL random variable params */
  UINT4 seed=0;              /* set to get seed from clock time */
  RandomParams *randomParams=NULL;

  CHAR *pos1=NULL, *pos2=NULL;
  INT4 i=0, j=0, k=0, n=0;

  REAL8 like1[1], like2[1];
  REAL8 logL1=0., logL2=0.; /* log likelihoods */
  REAL8 prior=0., priorNew=0., priorTmp=0., priorTmpNew=0.;
  REAL8 ratio=0.;           /* logL2 - logL1 = log(L2/L1) */
  REAL8 log_invtemp=0.;

  FILE *fp=NULL;
  CHAR outFile[256];

  /* variables for pulsar parameters */
  INT4 matTrue=0;
  BinaryPulsarParams pulsarParams, pulsarParamsNew;

  REAL8Array *chol=NULL, *invmat=NULL;
  REAL8 determinant;

  ParamData *paramData=NULL, *randVals=NULL, *vals=NULL;
  INT4Vector *matPos=NULL; /* position of parameters in ParamData */

  BarycenterInput baryinput = empty_BarycenterInput;
  REAL8Vector *phi1[numDets], *phi2=NULL;

  INT4 iterations = input.mcmc.iterations + input.mcmc.burnIn;  
  INT4 burnInLength = input.mcmc.burnIn; /* length of burn in */

  INT4 acc=0, rej=0; /* count acceptance and rejection of new point */

  INT4 onlyonce=0; /* use this variable only once */

  /* read the TEMPO par file for the pulsar */
  XLALReadTEMPOParFile( &pulsarParamsFixed, input.parFile );

  if( verbose ){
    fprintf(stderr, "Performing an MCMC for %s with %d iterations.\n",
      det, input.mcmc.iterations);
  }

  /* set up random parameters */
  randomParams = XLALCreateRandomParams( seed );

  /* work out how many glitches have been input */
  if( input.mcmc.nGlitches > 0 ){
    extraVars = XLALCalloc(nGlitches, sizeof(IntrinsicPulsarVariables));
    extraVarsNew = XLALCalloc(nGlitches, sizeof(IntrinsicPulsarVariables));
    glitchTimes = XLALCalloc(nGlitches, sizeof(REAL8));

    /* start and end point of the data before each glitch */
    g1 = XLALCalloc(numDets, sizeof(INT4 *));
    g2 = XLALCalloc(numDets, sizeof(INT4 *));

    for( j = 0 ; j < numDets ; j++ ){
      g1[j] = XLALCalloc(nGlitches + 1, sizeof(INT4));
      g2[j] = XLALCalloc(nGlitches + 1, sizeof(INT4));
    }

    /* check that the number of glitch times entered as the same as the number
       of glitches specified (count commas in the glitch time string */
    pos1 = input.mcmc.glitchTimes;
    j=0;
    while( (pos2 = strchr(pos1, ',')) != NULL ){
      pos1 = pos2 + 1;
      j++;
    }
    if( j != input.mcmc.nGlitches-1 ){
      fprintf(stderr, "Error... number of glitches specified is not the same as\
 the number of glitch times given!\n");
      exit(0);
    }
    pos1 = NULL;
    pos2 = NULL;

    /* glitch times are seperated by commas so count them up */
    for( i = 0 ; i < nGlitches ; i++ ){
      if( nGlitches == 1 )
        glitchTimes[i] = LALTDBMJDtoGPS(atof(input.mcmc.glitchTimes));
      else{
        if( i == 0 )
          pos1 = input.mcmc.glitchTimes;/*string starts "*/
        else
          pos1 = pos2+1; /* first position equals the previous value */

        gtimestr = NULL;
        gtimestr = XLALMalloc(256*sizeof(CHAR));

        /* values are delimited by commas , */
        if( i < nGlitches - 1 ){
          pos2 = strchr(pos1, ',');
          XLALStringCopy(gtimestr, pos1, (pos2-pos1)+1);/*copy time of glitch*/
        }
        else
          gtimestr = XLALStringDuplicate(pos1);

        /* crude check that the time is valid */
        if( atof(gtimestr) == 0. ){
          fprintf(stderr, "Error... Glitch time is not valid!\n");
          exit(0);
        }

        glitchTimes[i] = LALTDBMJDtoGPS(atof(gtimestr)); /* convert to GPS */

        XLALFree(gtimestr);
      }

      /* set initial values for extra params (all start at same point) */
      if( i==0 ){
        /* extraVars[i].h0 = input.mesh.minVals.h0 +
          (REAL8)XLALUniformDeviate(randomParams) *
          (input.mesh.maxVals.h0 - input.mesh.minVals.h0); */
        extraVars[i].phi0 = input.mesh.minVals.phi0 +
          (REAL8)XLALUniformDeviate(randomParams) *
          (input.mesh.maxVals.phi0 - input.mesh.minVals.phi0);
      }
      else{
        /* extraVars[i].h0 = extraVars[i-1].h0; */
        extraVars[i].phi0 = extraVars[i-1].phi0;
      }

      /* find the position before and after the glitch +/- the cut-off time/2 */
      for( j = 0 ; j < numDets ; j++ ){
        if( i == 0 )
          g1[j][i] = 0;

        for( k = 0 ; k < (INT4)data[j].times->length ; k++ ){
          /* first point after the last glitch */
          if( i > 0 && data[j].times->data[k] < glitchTimes[i-1] +
              input.mcmc.glitchCut/2. ){
              g1[j][i] = k+2;
           }

          /* end point of segment before current glitch */
          if( data[j].times->data[k] > glitchTimes[i] -
            input.mcmc.glitchCut/2. ){
            if( g2[j][i] == 0. ) /* check if it's been set already */
              g2[j][i] = k;

            if(i != nGlitches - 1)
              break;
          }

          /* point from final glitch to end of data */
          if( data[j].times->data[k] > glitchTimes[i] +
             input.mcmc.glitchCut/2. && i == nGlitches - 1 ){
            g1[j][i+1] = k+1;
            g2[j][i+1] = data[j].times->length-1;

            break;
          }
        }
      }
    }

    /* put each section of data into a new structure */
    glitchData = XLALCalloc(numDets, sizeof(DataStructure *));

    for( i = 0 ; i < numDets ; i++ ){
      glitchData[i] = XLALCalloc(nGlitches + 1, sizeof(DataStructure));

      for( j = 0 ; j < nGlitches + 1 ; j++ ){
        glitchData[i][j].data = NULL;
        glitchData[i][j].data = XLALCreateCOMPLEX16Vector(g2[i][j]-g1[i][j]+1);

        glitchData[i][j].times = NULL;
        glitchData[i][j].times = XLALCreateREAL8Vector(g2[i][j]-g1[i][j]+1);

        glitchData[i][j].chunkLengths = NULL;
        glitchData[i][j].chunkLengths =
          XLALCreateINT4Vector(g2[i][j]-g1[i][j]+1);
        glitchData[i][j].chunkMin = data[i].chunkMin;
        glitchData[i][j].chunkMax = data[i].chunkMax;

        n=0;
        for( k = g1[i][j] ; k <= g2[i][j] ; k++ ){
          glitchData[i][j].data->data[n] = data[i].data->data[k];
          glitchData[i][j].times->data[n] = data[i].times->data[k];
          n++;
        }

        /* get chunk lengths and data sums */
        get_chunk_lengths(glitchData[i][j]);

        glitchData[i][j].sumData = NULL;
        glitchData[i][j].sumData =
          XLALCreateREAL8Vector(glitchData[i][j].chunkLengths->length);

        sum_data(glitchData[i][j]);

        /* set lookup tables */
        glitchData[i][j].lookupTable = data[i].lookupTable;
      }
    }
  }

  /* if there's an input covarince matrix file then set up parameters */
  if( input.matrixFile != NULL || input.usecov )
    matTrue = 1;

  if( matTrue ){
    REAL8Array *cormat=NULL, *covmat=NULL, *posdef=NULL, *tempinvmat=NULL;

    paramData = XLALMalloc( MAXPARAMS*sizeof(ParamData) );

    /* get correlation matrix */
    if( ( cormat = read_correlation_matrix( input.matrixFile, pulsarParamsFixed,
paramData ) ) == NULL ){
      fprintf(stderr, "Error... couldn't read in correlation matrix file!\n");
      exit(0);
    }

    /* allocate memory for inverse matrix */
    tempinvmat = XLALCreateREAL8Array( cormat->dimLength );

    /* check if correlation matrix is positive definite if necessary */
    if( (posdef = check_positive_definite( cormat )) == NULL ){
      /* set covariance matrix */
      covmat = create_covariance_matrix( paramData, cormat, 0 );

      /* cholesky decompose matrix */
      chol = cholesky_decomp(covmat, "lower");
    }
    else{
      /* set covariance matrix */
      covmat = create_covariance_matrix( paramData, posdef, 0 );

      /* cholesky decompose covariance matrix */
      chol = cholesky_decomp(covmat, "lower");

      XLALDestroyREAL8Array( posdef );
    }

    /* create inverse of correlation matrix */
    LAL_CALL( LALDMatrixInverse( &status, &determinant, cormat, tempinvmat),
      &status );

    /* set inverse of covariance matrix */
    invmat = create_covariance_matrix( paramData, tempinvmat, 1 );

    XLALDestroyREAL8Array( tempinvmat );
    XLALDestroyREAL8Array( cormat );
    XLALDestroyREAL8Array( covmat );

    /* generate random parameters */
    randVals = multivariate_normal_deviates( chol, paramData, randomParams );

    vals = XLALMalloc( MAXPARAMS*sizeof(ParamData) );
    memcpy( vals, randVals, MAXPARAMS*sizeof(ParamData) );

    /* get the posistions of the parameters within paramData */
    i=0;
    j=0;
    matPos = XLALCreateINT4Vector( MAXPARAMS );
    while(i < MAXPARAMS){
      if( paramData[i].matPos != 0 ){
        matPos->data[paramData[i].matPos-1] = i;
        j++;
      }
      i++;
    }
    matPos = XLALResizeINT4Vector(matPos, j);

    /* set up initial pulsar parameters */
    memcpy( &pulsarParams, &pulsarParamsFixed, sizeof(pulsarParamsFixed) );

    set_mcmc_pulsar_params( &pulsarParams, randVals );

    /* check eccentricities so that 0 <= e < 1 i.e. circular or elliptical
       orbits */
    if( pulsarParams.e < 0. )
      pulsarParams.e = fabs(pulsarParams.e);
    else if( pulsarParams.e >= 1. )
      pulsarParams.e = 1. - fmod(pulsarParams.e, 1.);

    if( pulsarParams.e2 < 0. )
      pulsarParams.e2 = fabs(pulsarParams.e2);
    else if( pulsarParams.e2 >= 1. )
      pulsarParams.e2 = 1. - fmod(pulsarParams.e2, 1.);

    if( pulsarParams.e3 < 0. )
      pulsarParams.e3 = fabs(pulsarParams.e3);
    else if( pulsarParams.e3 >= 1. )
      pulsarParams.e3 = 1. - fmod(pulsarParams.e3, 1.);

    /* check sini parameter and make sure it's between -1 and 1 */
    if( pulsarParams.s > 1. )
      pulsarParams.s = 1. - fmod(pulsarParams.s, 1.);
    else if( pulsarParams.s < -1. )
      pulsarParams.s = -1. - fmod(pulsarParams.s, 1.);
  }
  else
    edat = NULL; /* make sure edat is NULL as this is how later functions
                    decide whether or not to perform MCMC over all
                    parameters */

  /* set initial chain parameters */
  vars.h0 = input.mesh.minVals.h0 + (REAL8)XLALUniformDeviate(randomParams) *
             (input.mesh.maxVals.h0 - input.mesh.minVals.h0);

  /* if( input.mcmc.nGlitches > 0 )
    vars.h0 = extraVars[0].h0; */

  vars.phi0 = input.mesh.minVals.phi0 + (REAL8)XLALUniformDeviate(randomParams)
             * (input.mesh.maxVals.phi0 - input.mesh.minVals.phi0);
  vars.psi = input.mesh.minVals.psi + (REAL8)XLALUniformDeviate(randomParams) *
              (input.mesh.maxVals.psi - input.mesh.minVals.psi);
  vars.ci = input.mesh.minVals.ci + (REAL8)XLALUniformDeviate(randomParams) *
            (input.mesh.maxVals.ci - input.mesh.minVals.ci);

  /* fprintf(stderr, "Give me a start h0 value for the chain:\n");
  fscanf(stdin, "%lf", &vars.h0);

  fprintf(stderr, "Give me a start phi0 value for the chain:\n");
  fscanf(stdin, "%lf", &vars.phi0);

  fprintf(stderr, "Give me a start psi value for the chain:\n");
  fscanf(stdin, "%lf", &vars.psi);

  fprintf(stderr, "Give me a start cos(iota) value for the chain:\n");
  fscanf(stdin, "%lf", &vars.ci); */
  
  vars.Xplus = 0.5*(1.+vars.ci*vars.ci);
  vars.Xcross = vars.ci;
  vars.Xpsinphi_2 = 0.5*vars.Xplus*sin(vars.phi0);
  vars.Xcsinphi_2 = 0.5*vars.Xcross*sin(vars.phi0);
  vars.Xpcosphi_2 = 0.5*vars.Xplus*cos(vars.phi0);
  vars.Xccosphi_2 = 0.5*vars.Xcross*cos(vars.phi0);

  for( i = 0 ; i < nGlitches ; i++ ){
    extraVars[i].h0 = vars.h0;
    extraVars[i].Xpsinphi_2 = 0.5*vars.Xplus * sin(extraVars[i].phi0);
    extraVars[i].Xcsinphi_2 = 0.5*vars.Xcross * sin(extraVars[i].phi0);
    extraVars[i].Xpcosphi_2 = 0.5*vars.Xplus * cos(extraVars[i].phi0);
    extraVars[i].Xccosphi_2 = 0.5*vars.Xcross * cos(extraVars[i].phi0);
    extraVars[i].psi = vars.psi;
  }

  input.mesh.h0Steps = 1;
  input.mesh.phiSteps = 1;
  input.mesh.ciotaSteps = 1;
  input.mesh.psiSteps = 1;

  /* get data chunk lengths - only if no glitches here */
  if( nGlitches == 0 ){
    for( k = 0 ; k < numDets ; k++ ){
      /* create vector of data segment length */
      data[k].chunkLengths = NULL;
      data[k].chunkLengths = XLALCreateINT4Vector(
      (INT4)ceil(data[k].data->length/data[k].chunkMin) );

      get_chunk_lengths(data[k]);

      /* allocate memory for summations */
      data[k].sumData = NULL;
      data[k].sumData = XLALCreateREAL8Vector(data[k].chunkLengths->length);

      sum_data(data[k]);
    }
  }

  /* open output file */
  if( input.mcmc.outputBI == 0 )  
    sprintf(outFile, "%s/MCMCchain_%s_%s", input.outputDir, input.pulsar, det);
  else{ /* append number of burn in steps to the file name */
    sprintf(outFile, "%s/MCMCchain_%s_%s_burn_in_%d", input.outputDir, 
      input.pulsar, det, input.mcmc.burnIn );
  }

  if( (fp=fopen(outFile, "w")) == NULL ){
    fprintf(stderr, "Error... Can't open MCMC output chain file!\n");
    exit(0);
  }

  /* write MCMC chain header info */
  fprintf(fp, "%% MCMC for %s with %s data using %d iterations\n",
    input.pulsar, det, input.mcmc.iterations);
  fprintf(fp, "%% log(L)    \th0          \tphi0 (rads) \tcos(iota)\tpsi \
(rads)");
  for( j = 0 ; j < nGlitches ; j++ )
    fprintf(fp, "\tphi0(%d) ", j+2);
    /* fprintf(fp, "\th0(%d)       \tphi0(%d) ", j+2, j+2); */
  if( matTrue ){
    for( j = 0 ; j < (INT4)matPos->length ; j++ )
        fprintf(fp, "\t%s", paramData[matPos->data[j]].name);
  }
  fprintf(fp, "\n");

  /* create vector for random Gaussian numbers */
  /* randNum = XLALCreateREAL4Vector(4+2*nGlitches); */
  randNum = XLALCreateREAL4Vector(4+nGlitches);

  if( matTrue && numDets == 1 ){
    /* set up detector location - if not doing joint analysis */
    baryinput.site = *detPos;
    baryinput.site.location[0] /= LAL_C_SI;
    baryinput.site.location[1] /= LAL_C_SI;
    baryinput.site.location[2] /= LAL_C_SI;

    /* set phase of initial heterodyne */
    if( (phi1[0] = get_phi( data[0], pulsarParamsFixed, baryinput, edat )) ==
      NULL ){
      fprintf(stderr, "Error... Phase generation produces NULL!\n");
      exit(0);
    }
  }

  /* if a annealing temperature is set then set the value of log(1/temp) */
  if( input.mcmc.temperature > 0. )
    log_invtemp = log( 1. / input.mcmc.temperature );

  /* if we also want to output the burn in chain set to -1 */
  if( input.mcmc.outputBI == 1 )
    burnInLength = -1;

  fprintf(stderr, "Entering MCMC stage\n");

  /*=================== MCMC LOOP =====================*/
  for( i = 0 ; i < iterations ; i++ ){
    REAL4 sp=0., cp=0.; /* sin and cos values */
    INT4 nege=0;        /* eccentricity check */

    if( verbose ){
      fprintf(stderr, "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
      fprintf(stderr, "%06.2lf%% complete", ((REAL8)(i)+1.)*100. /
        (REAL8)(input.mcmc.iterations + input.mcmc.burnIn));
    }

    /* get new values of parameters using Gaussian proposal distributions */
    XLALNormalDeviates(randNum, randomParams);

    varsNew.h0 = vars.h0 + input.mcmc.sigmas.h0*randNum->data[0];

    /* below0 = 0; */
    for( j = 0 ; j < nGlitches ; j++ ){
      /* extraVarsNew[j].h0 = extraVars[j].h0 +
        input.mcmc.sigmas.h0*randNum->data[4+2*j];
      extraVarsNew[j].phi0 = extraVars[j].phi0 +
        input.mcmc.sigmas.phi0*randNum->data[4+2*j+1]; */
      extraVarsNew[j].phi0 = extraVars[j].phi0 +
        input.mcmc.sigmas.phi0*randNum->data[4+j];

      if( extraVarsNew[j].phi0 < 0. ){
        extraVarsNew[j].phi0 = fmod(extraVarsNew[j].phi0, LAL_TWOPI);
        extraVarsNew[j].phi0 += LAL_TWOPI;
      }
      else if( extraVarsNew[j].phi0 > LAL_TWOPI )
        extraVarsNew[j].phi0 = fmod(extraVarsNew[j].phi0, LAL_TWOPI);

      /* if any of the h0 jump below zero then write out data below */
      /* if( extraVarsNew[j].h0 < 0. )
        below0 = 1; */

      /* make it so that values of h0 after a glitch must be smaller than
         before the glitch i.e. the pulsar always relaxes to a lower energy
         state - if this happens set below0 = 1 to output the current state */
      /* if( j==0 ){
        if( extraVarsNew[j].h0 > varsNew.h0 )
          below0 = 1;
      }
      else{
        if( extraVarsNew[j].h0 > extraVarsNew[j-1].h0 )
          below0 = 1;
      } */
    }

    if( matTrue ){
      /* generate new pulsars parameters from the pulsar parameter covariance */
      randVals = multivariate_normal_deviates( chol, vals, randomParams );
      memcpy(&pulsarParamsNew, &pulsarParams, sizeof(BinaryPulsarParams));
      set_mcmc_pulsar_params( &pulsarParamsNew, randVals );

      /* I've discovered that TEMPO can produce negative eccentricities (when
         they're close to zero), so will allow this to happen here, but still
         won't allow it to go over 1 */
      if( /* pulsarParamsNew.e < 0.  || */ pulsarParamsNew.e >= 1. ||
          /* pulsarParamsNew.e2 < 0. || */ pulsarParamsNew.e2 >= 1. ||
          /* pulsarParamsNew.e3 < 0. || */ pulsarParamsNew.e3 >= 1. ||
          pulsarParamsNew.s > 1. || pulsarParamsNew.s < -1. )
        nege = 1; 
    }

    /* if h0 jumps negative, or eccentricity is negative or greater than 1 then
       this is equivalent to having jumped outside our prior range so the
       likelihood is always zero and this move always rejected - therefore it's
       quickest just to output the only step now and move on to the next step */
    /* if( ( varsNew.h0 < 0. || below0 == 1 || nege == 1 ) && i > 0 ){ */
    if( ( varsNew.h0 < 0. || nege == 1 ) && i > 0 ){
      if( fmod(i, input.mcmc.outputRate) == 0. && i >= burnInLength ){
        fprintf(fp, "%le\t%le\t%lf\t%lf\t%lf", logL1, vars.h0, vars.phi0,
          vars.ci, vars.psi);

        for( j = 0 ; j < nGlitches ; j++ )
          fprintf(fp, "\t%lf", extraVars[j].phi0);
          /* fprintf(fp, "\t%le\t%lf", extraVars[j].h0, extraVars[j].phi0); */

        if( matTrue ){
          /* print out pulsar parameters */
          for( j = 0 ; j < (INT4)matPos->length ; j++ ){
            fprintf(fp, "\t%.17le", vals[matPos->data[j]].val -
              paramData[matPos->data[j]].val);
          }
        }
        fprintf(fp, "\n");
      }

      if( i > input.mcmc.burnIn - 1 )
        rej++; /* count rejected steps */

      continue;
    }
    /* else if( ( varsNew.h0 < 0. || below0 == 1 || nege == 1 ) && i == 0 ){ */
    else if( ( varsNew.h0 < 0. || nege == 1 ) && i == 0 ){
      onlyonce = 1; /* if h0 goes below zero on the first step then we still 
                       have to calculate logL1, so continue but make sure 
                       logL2 gets set to -Inf (or close to!) later on */
      /* set values of h0 so that they aren't negative, as this could screw
         up other functions */
      varsNew.h0 = 1e-30;
      pulsarParamsNew.e = 0.;
      pulsarParamsNew.e2 = 0.;
      pulsarParamsNew.e3 = 0.;

      for( j = 0 ; j < nGlitches ; j++ ) extraVarsNew[j].h0 = 1e-30;
    }

    varsNew.phi0 = vars.phi0 + input.mcmc.sigmas.phi0*randNum->data[1];
    varsNew.psi = vars.psi + input.mcmc.sigmas.psi*randNum->data[2];
    varsNew.ci = vars.ci + input.mcmc.sigmas.ci*randNum->data[3];

    /* wrap parameters around or bounce */
    if( varsNew.ci > 1.0 )
      varsNew.ci = 1. - fmod(varsNew.ci, 1.);
    else if( varsNew.ci < -1.0 )
      varsNew.ci = -1. - fmod(varsNew.ci, 1.);

    if( varsNew.phi0 < 0. ){
      varsNew.phi0 = fmod(varsNew.phi0, LAL_TWOPI);
      varsNew.phi0 += LAL_TWOPI;
    }
    else if( varsNew.phi0 > LAL_TWOPI )
      varsNew.phi0 = fmod(varsNew.phi0, LAL_TWOPI);

    if( varsNew.psi > LAL_PI/4. ){
      varsNew.psi = fmod(varsNew.psi, LAL_PI/4.) - LAL_PI/4.;
      varsNew.phi0 = fmod(varsNew.phi0 + LAL_PI, LAL_TWOPI);

      for( j = 0 ; j < nGlitches ; j++ )
        extraVarsNew[j].phi0 = fmod(extraVarsNew[j].phi0 + LAL_PI, LAL_TWOPI);
    }
    else if(varsNew.psi < -LAL_PI/4.){
      varsNew.psi = fmod(varsNew.psi, LAL_PI/4.) + LAL_PI/4.;
      varsNew.phi0 = fmod(varsNew.phi0 + LAL_PI, LAL_TWOPI);

      for( j = 0 ; j < nGlitches ; j++ )
        extraVarsNew[j].phi0 = fmod(extraVarsNew[j].phi0 + LAL_PI, LAL_TWOPI);
    }

    /* set combined parameters */
    varsNew.Xplus = 0.5*(1.+varsNew.ci*varsNew.ci);
    varsNew.Xcross = varsNew.ci;
    sin_cos_LUT( &sp, &cp, varsNew.phi0 );
    varsNew.Xpsinphi_2 = 0.5*varsNew.Xplus*sp;
    varsNew.Xcsinphi_2 = 0.5*varsNew.Xcross*sp;
    varsNew.Xpcosphi_2 = 0.5*varsNew.Xplus*cp;
    varsNew.Xccosphi_2 = 0.5*varsNew.Xcross*cp;

    for( j = 0 ; j < nGlitches ; j++ ){
      sin_cos_LUT( &sp, &cp, extraVarsNew[j].phi0 );
      extraVarsNew[j].h0 = varsNew.h0;
      extraVarsNew[j].Xpsinphi_2 = 0.5*varsNew.Xplus * sp;
      extraVarsNew[j].Xcsinphi_2 = 0.5*varsNew.Xcross * sp;
      extraVarsNew[j].Xpcosphi_2 = 0.5*varsNew.Xplus * cp;
      extraVarsNew[j].Xccosphi_2 = 0.5*varsNew.Xcross * cp;
      extraVarsNew[j].psi = varsNew.psi;
    }

    /* calculate log likelihoods */
    /* only calculate the likelhood twice in the first loop */
    if( i == 0 ){
      for( k = 0 ; k < numDets ; k++ ){
        if( matTrue && numDets > 1 ){
          /* set up detector location for joint likelihood */
          baryinput.site = detPos[k];
          baryinput.site.location[0] /= LAL_C_SI;
          baryinput.site.location[1] /= LAL_C_SI;
          baryinput.site.location[2] /= LAL_C_SI;

          if( (phi1[k] = get_phi( data[k], pulsarParamsFixed, baryinput, edat )) == NULL ){
            fprintf(stderr, "Error... Phase generation produces NULL!");
            exit(0);
          }
        }

        /* get new phase */
        /* set up the deltaphi, and put it into the phi2 variable */
        if( matTrue ){
          if( (phi2 = get_phi( data[k], pulsarParams, baryinput, edat )) == NULL ){
            fprintf(stderr, "Error... Phase generation produces NULL!");
            exit(0);
          }

          for( j=0; j<(INT4)data[k].times->length; j++ )
            phi2->data[j] = fmod(phi2->data[j] - phi1[k]->data[j], 1.);
        }
        else
          phi2 = NULL; /* if there's no cov file log_likelihood requires phi to be NULL */

        /* first likelihood */
        if( nGlitches == 0 ){
          input.mesh.minVals.h0 = vars.h0;
          log_likelihood(like1, data[k], vars, input.mesh, phi2);
          logL1 += *like1;
        }
        else{
          for( j = 0 ; j < nGlitches + 1 ; j++ ){
            if( j == 0 ){
              input.mesh.minVals.h0 = vars.h0;

              log_likelihood(like1, glitchData[k][j], vars, input.mesh, phi2);
            }
            else{
              /* input.mesh.minVals.h0 = extraVars[j-1].h0; */
              log_likelihood(like1, glitchData[k][j], extraVars[j-1],
                input.mesh, phi2);
            }
            logL1 += *like1;
          }
        }
        XLALDestroyREAL8Vector( phi2 );
        phi2 = NULL;
      }
    }

    for( k = 0 ; k < numDets ; k++ ){
      /* second likelihood - for new position in parameter space */
      if( matTrue && numDets > 1 ){
        /* set up detector location for joint likelihood */
        baryinput.site = detPos[k];
        baryinput.site.location[0] /= LAL_C_SI;
        baryinput.site.location[1] /= LAL_C_SI;
        baryinput.site.location[2] /= LAL_C_SI;
      }

      /* get new phase */
      /* set up the deltaphi, and put it into the phi2 variable */
      if( matTrue ){
        if( (phi2 = get_phi( data[k], pulsarParamsNew, baryinput, edat )) == NULL ){
          fprintf(stderr, "Error... Phase generation produces NULL!");
          exit(0); 
        }

        for( j=0; j<(INT4)data[k].times->length; j++ )
          phi2->data[j] = fmod(phi2->data[j] - phi1[k]->data[j], 1.);
      }
      else
        phi2 = NULL;

      if( nGlitches == 0 ){
        input.mesh.minVals.h0 = varsNew.h0;
        log_likelihood(like2, data[k], varsNew, input.mesh, phi2);
        logL2 += *like2;
      }
      else{
        for( j = 0 ; j < nGlitches + 1 ; j++ ){
          if( j == 0 ){
            input.mesh.minVals.h0 = varsNew.h0;
            log_likelihood(like2, glitchData[k][j], varsNew, input.mesh, phi2);
          }
          else{
            /* input.mesh.minVals.h0 = extraVarsNew[j-1].h0; */
            log_likelihood(like2, glitchData[k][j], extraVarsNew[j-1],
              input.mesh, phi2);
          }
          logL2 += *like2;
        }
      }
      XLALDestroyREAL8Vector( phi2 );
      phi2 = NULL;
    }

    /* most importantly use the covariance matrix as a prior on the pulsar
       parameters */
    priorNew = 0.;
    if( matTrue ){
      for( j=0; j<(INT4)matPos->length ; j++){
        if(i==0)
          priorTmp = 0.;

        priorTmpNew = 0.;

        for( k=0; k<(INT4)matPos->length ; k++){
          if(i==0){
            priorTmp += (vals[matPos->data[k]].val -
              paramData[matPos->data[k]].val) *
              invmat->data[j * (INT4)invmat->dimLength->data[0]+ k];
          }

          priorTmpNew += (randVals[matPos->data[k]].val -
              paramData[matPos->data[k]].val) *
              invmat->data[j * (INT4)invmat->dimLength->data[0] + k];
        }

        if(i==0){
          prior += priorTmp*(vals[matPos->data[j]].val -
            paramData[matPos->data[j]].val);
        }

        priorNew += priorTmpNew*(randVals[matPos->data[j]].val -
          paramData[matPos->data[j]].val);
      }

      if(i==0){
        prior *= -0.5;
        logL1 += prior;
      }

      priorNew *= -0.5;

      logL2 += priorNew;
    }

    /* set and apply priors on h0, psi, cos(iota) and phi if requested */
    if( input.usepriors ){
      if(i==0){
        input.priors.vars.h0 = vars.h0;
        input.priors.vars.phi0 = vars.phi0;
        input.priors.vars.ci = vars.ci;
        input.priors.vars.psi = vars.psi;

        logL1 += log_prior(input.priors, input.mesh);

        if( nGlitches > 0 ){
          for( j = 0 ; j < nGlitches ; j++ ){
            /* input.priors.vars.h0 = extraVars[j].h0; */
            input.priors.vars.phi0 = extraVars[j].phi0;

            logL1 += log_prior(input.priors, input.mesh);
          }
        }
      }

      input.priors.vars.h0 = varsNew.h0;
      input.priors.vars.phi0 = varsNew.phi0;
      input.priors.vars.ci = varsNew.ci;
      input.priors.vars.psi = varsNew.psi;

      logL2 += log_prior(input.priors, input.mesh);

      if( nGlitches > 0 ){
        for( j = 0 ; j < nGlitches ; j++ ){
          /* input.priors.vars.h0 = extraVarsNew[j].h0; */
          input.priors.vars.phi0 = extraVarsNew[j].phi0;

          logL2 += log_prior(input.priors, input.mesh);
        }
      }
    }

    /* if this is the first time in the chain an h0 was negative then set
       the value of logL2 to something close to -Inf */
    if( onlyonce == 1 ){
      logL2 = -INFINITY;
      onlyonce = 0; /* reset value */
    }

    /* accept new values if Lnew/Lold >=1 (or log(Lnew/Lold) >= 0) */
    /* include simulated annealing factor */
    ratio = logL2 - logL1;
    if( i < input.mcmc.burnIn ){
      ratio *= input.mcmc.temperature * exp( log_invtemp *
        (double)i/(double)input.mcmc.burnIn );
    }

    /* accept new step (Metropolis-Hastings algorithm) if Lnew/Lold > 1)
       or with a probability Lnew/Lold if < 1 */
    if( ratio - log(XLALUniformDeviate(randomParams)) >= 0. ){
      vars = varsNew;

      /* update vals structure */
      if( matTrue ){
        memcpy(vals, randVals, MAXPARAMS*sizeof(ParamData));
        memcpy(&pulsarParams, &pulsarParamsNew, sizeof(BinaryPulsarParams));
      }

      for( j = 0 ; j < nGlitches ; j++)
        extraVars[j] = extraVarsNew[j];

      logL1 = logL2;

      if( i > input.mcmc.burnIn - 1 )
        acc++; /* count acceptance number */
    }
    else{
      if( i > input.mcmc.burnIn - 1 )
        rej++; /* count rejection number */
    }

    /* printf out chains */
    if( fmod(i, input.mcmc.outputRate) == 0. && i >= burnInLength ){
      fprintf(fp, "%le\t%le\t%lf\t%lf\t%lf", logL1, vars.h0, vars.phi0, vars.ci,
        vars.psi);

      for( j = 0 ; j < nGlitches ; j++ )
        fprintf(fp, "\t%lf", extraVars[j].phi0);
        /* fprintf(fp, "\t%le\t%lf", extraVars[j].h0, extraVars[j].phi0); */

      if( matTrue ){
        /* print out pulsar parameters */
        for( j = 0 ; j < (INT4)matPos->length ; j++ ){
          /* output difference in parameter from heterodyne parameter */
          fprintf(fp, "\t%.17le", vals[matPos->data[j]].val -
            paramData[matPos->data[j]].val);
        }
      }

      fprintf(fp, "\n");
    }

    logL2 = 0.;
  }
  /*===============================================*/

  if( verbose ) fprintf(stderr, "\n");

  /* print out acceptance ratio */
  if( verbose ){
    if( rej > 0 )
      fprintf(stderr, "Acceptance ratio %.3lf\n", (double)acc/(double)rej);
    else
      fprintf(stderr, "No step rejected!!\n");
  }

  /* output acceptance ratio to end of the chain */
  if( rej > 0 )
    fprintf(fp, "%% Acceptance ratio = %.5lf\n", (double)acc/(double)rej);
  else
    fprintf(fp, "%% Acceptance ratio - no step rejected!\n");

  /* destroy vectors */
  XLALDestroyREAL4Vector(randNum);

  for( j = 0 ; j < nGlitches ; j++ ){
    for( i = 0 ; i < numDets ; i++ ){
      XLALDestroyCOMPLEX16Vector(glitchData[i][j].data);
      XLALDestroyREAL8Vector(glitchData[i][j].times);
      XLALDestroyINT4Vector(glitchData[i][j].chunkLengths);
      XLALDestroyREAL8Vector(glitchData[i][j].sumData);
    }
  }

  if( matTrue ){
    XLALFree( paramData );
    XLALFree( vals );
    XLALFree( randVals );

    XLALDestroyREAL8Array( invmat );
    XLALDestroyREAL8Array( chol );

    XLALDestroyINT4Vector( matPos );

    for( i=0; i < numDets; i++ )
      XLALDestroyREAL8Vector( phi1[i] );
  }

  fclose(fp);
}



/* function to return a vector of the pulsar phase for each data point */
REAL8Vector *get_phi( DataStructure data, BinaryPulsarParams params,
  BarycenterInput bary, EphemerisData *edat ){
  INT4 i=0;

  REAL8 T0=0., DT=0., DTplus=0., deltat=0., deltat2=0.;
  REAL8 interptime = 1800.; /* calulate every 30 mins (1800 secs) */

  EarthState earth, earth2;
  EmissionTime emit, emit2;
  REAL8 emitdt=0.;

  BinaryPulsarInput binput;
  BinaryPulsarOutput boutput;

  REAL8Vector *phis=NULL;

  /* if edat is NULL then return a NULL poniter */
  if( edat == NULL )
    return NULL;

   /* set the position and frequency epochs if not already set */
  if(params.pepoch == 0. && params.posepoch != 0.)
    params.pepoch = params.posepoch;
  else if(params.posepoch == 0. && params.pepoch != 0.)
    params.posepoch = params.pepoch;

  /* allocate memory for phases */
  phis = XLALCreateREAL8Vector( data.times->length );

  /* set 1/distance if parallax or distance value is given (1/sec) */
  if( params.px != 0. )
    bary.dInv = params.px*1e-3*LAL_C_SI/LAL_PC_SI;
  else if( params.dist != 0. )
    bary.dInv = LAL_C_SI/(params.dist*1e3*LAL_PC_SI);
  else
    bary.dInv = 0.;

  for( i=0; i<(INT4)data.times->length; i++){
    T0 = params.pepoch;

    DT = data.times->data[i] - T0;

    /* only do call the barycentring routines every 30 minutes, otherwise just
       linearly interpolate between them */
    if( i==0 || DT > DTplus ){
      bary.tgps.gpsSeconds = (UINT8)floor(data.times->data[i]);
      bary.tgps.gpsNanoSeconds =
(UINT8)floor((fmod(data.times->data[i],1.)*1e9));

      bary.delta = params.dec +
        (data.times->data[i]-params.posepoch) * params.pmdec;
      bary.alpha = params.ra + (data.times->data[i]-params.posepoch) *
         params.pmra/cos(bary.delta);

      /* call barycentring routines */
      XLAL_CHECK_NULL( XLALBarycenterEarth( &earth, &bary.tgps, edat ) ==
                       XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_NULL( XLALBarycenter( &emit, &bary, &earth ) ==
                       XLAL_SUCCESS, XLAL_EFUNC );

      /* add interptime to the time */
      DTplus = DT + interptime;
      bary.tgps.gpsSeconds = (UINT8)floor(data.times->data[i]+interptime);
      bary.tgps.gpsNanoSeconds =
(UINT8)floor((fmod(data.times->data[i]+interptime,1.)*1e9));

      /* No point in updating the positions as difference will be tiny */
      XLAL_CHECK_NULL( XLALBarycenterEarth( &earth2, &bary.tgps, edat ) ==
                       XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_NULL( XLALBarycenter( &emit2, &bary, &earth2 ) ==
                       XLAL_SUCCESS, XLAL_EFUNC );
    }

    /* linearly interpolate to get emitdt */
    emitdt = emit.deltaT + (DT - (DTplus - interptime)) *
      (emit2.deltaT - emit.deltaT)/interptime;

    /* check if need to perform binary barycentring */
    if( params.model != NULL ){
      binput.tb = data.times->data[i] + emitdt;
      binput.earth = earth;
      
      XLALBinaryPulsarDeltaT( &boutput, &binput, &params );

      deltat = DT + emitdt + boutput.deltaT;
    }
    else
      deltat = DT + emitdt;

    /* work out phase */
    deltat2 = deltat*deltat;
    phis->data[i] = 2.*deltat*(params.f0 + 0.5*params.f1*deltat +
      SIXTH*params.f2*deltat2 + TWENTYFOURTH*params.f3*deltat*deltat2 
      + (1./120.)*params.f4*deltat2*deltat2 
      + (1./720.)*params.f5*deltat2*deltat2*deltat);
  }

  return phis;
}



/* function to get the lengths of consecutive chunks of data */
void get_chunk_lengths(DataStructure data){
  INT4 i=0, j=0, count=0;

  /* create vector of data segment length */
  while( 1 ){
    count++; /* counter */

    /* break clause */
    if( i > (INT4)data.data->length - 2 ){
      /* set final value of chunkLength */
      data.chunkLengths->data[j] = count;
      j++;
      break;
    }

    i++;

    /* if consecutive points are within 180 seconds of each other count as in
       the same chunk */
    if( data.times->data[i] - data.times->data[i-1] > 180. || count ==
      data.chunkMax ){
      data.chunkLengths->data[j] = count;
      count = 0; /* reset counter */

      j++;
    }
  }

  data.chunkLengths = XLALResizeINT4Vector(data.chunkLengths, j);
}



/* this functions set the value of a particular parameter from the structure
   ParamData into a BinaryPulsarParams structure - as this is what is taken
   in by the get_phi function */
void set_mcmc_pulsar_params( BinaryPulsarParams *pulsarParams, ParamData *data){
  /* go through params and set the next step in the MCMC */
  pulsarParams->f0 = data[0].val;
  pulsarParams->f1 = data[1].val;
  pulsarParams->f2 = data[2].val;
  pulsarParams->dec = data[3].val;
  pulsarParams->ra = data[4].val;
  pulsarParams->pmdec = data[5].val;
  pulsarParams->pmra = data[6].val;
  pulsarParams->x = data[7].val;
  pulsarParams->e = data[8].val;
  pulsarParams->T0 = data[9].val;
  pulsarParams->Pb = data[10].val;
  pulsarParams->w0 = data[11].val;
  pulsarParams->wdot = data[12].val;
  pulsarParams->gamma = data[13].val;
  pulsarParams->Pbdot = data[14].val;
  pulsarParams->s = data[15].val;
  pulsarParams->M = data[16].val;
  pulsarParams->m2 = data[17].val;
  pulsarParams->dth = data[18].val;
  pulsarParams->xdot = data[19].val;
  pulsarParams->edot = data[20].val;
  pulsarParams->x2 = data[21].val;
  pulsarParams->e2 = data[22].val;
  pulsarParams->T02 = data[23].val;
  pulsarParams->Pb2 = data[24].val;
  pulsarParams->w02 = data[25].val;
  pulsarParams->x3 = data[26].val;
  pulsarParams->e3 = data[27].val;
  pulsarParams->T03 = data[28].val;
  pulsarParams->Pb3 = data[29].val;
  pulsarParams->w03 = data[30].val;
  pulsarParams->xpbdot = data[31].val;
  pulsarParams->f3 = data[32].val;
  pulsarParams->f4 = data[33].val;
  pulsarParams->f5 = data[34].val;

  if( pulsarParams->model != NULL && !strcmp(pulsarParams->model, "ELL1") ){
    pulsarParams->eps1 = data[8].val;
    pulsarParams->eps2 = data[11].val;
    pulsarParams->e = 0.;
    pulsarParams->Tasc = data[9].val;
    pulsarParams->eps1dot = data[12].val;
    pulsarParams->eps2dot = data[20].val;
  }
}



/* this function performs Cholesky decomposition on M and outputs the
   lower, or upper triagular matrix depending if uOrl is set to "upper" or
   "lower" - if nothing is specified then the default is lower
   This is pretty much copied from the GSL function gsl_linalg_cholesky_decomp
   although this works with matrices where the values are close to zero
*/
REAL8Array *cholesky_decomp(REAL8Array *M, const CHAR* uOrl){
  INT4 i=0, j=0, k=0;

  REAL8 A_00=0., L_00=0.;

  REAL8Array *A=NULL;

  INT4 length=0;

  if( verbose ) fprintf(stderr, "Performing Cholesky decomposition of \
matrix\n");

  /* check dimensions are equal */
  if( M->dimLength->data[0] != M->dimLength->data[1] ){
    fprintf(stderr, "Error... input matrix has unequal dimensions!\n");
    exit(0);
  }

  length = M->dimLength->data[0];

  /* allocate memory */
  A = XLALCreateREAL8Array( M->dimLength );

  if(M == NULL || A == NULL){
    fprintf(stderr, "Error... input or output matrix is NULL!\n");
    exit(0);
  }

  /* initialise L be same as input matrix M */
  for(i=0; i < length; i++)
    for(j=0; j < length; j++)
      A->data[i*length + j] = M->data[i*length + j];

  A_00 = A->data[0];
  L_00 = sqrt(A_00);

  if( A_00 <= 0 )
    fprintf(stderr, "Error... matrix must be positive definite!\n");

  A->data[0] = L_00;

  if( length > 1 ){
    REAL8 A_10 = A->data[1*length + 0];
    REAL8 A_11 = A->data[1*length + 1];

    REAL8 L_10 = A_10/L_00;
    REAL8 diag = A_11 - L_10*L_10;
    REAL8 L_11 = sqrt(diag);

    if( diag <= 0 ){
      fprintf(stderr, "Error... input matrix is not pos def!\n");
      exit(0);
    }

    A->data[1*length + 0] = L_10;
    A->data[1*length + 1] = L_11;
  }

  for( k=2; k<length; k++ ){
    REAL8 A_kk = A->data[k*length + k];

    for( i=0; i<k; i++ ){
      REAL8 sum = 0.;

      REAL8 A_ki = A->data[k*length + i];
      REAL8 A_ii = A->data[i*length + i];

      REAL8 ci[length];
      REAL8 ck[length];

      for( j=0; j<length; j++ ){
        ci[j] = A->data[i*length + j];
        ck[j] = A->data[k*length + j];
      }

      if( i>0 ){
        for( j=0; j<i; j++ )
          sum += ci[j]*ck[j];
      }

      A_ki = (A_ki - sum) / A_ii;
      A->data[k*length + i] = A_ki;
    }

    {
      REAL8 sum = 0.;
      REAL8 diag = 0.;

      for( j = 0 ; j < k ; j++ )
        sum += A->data[k*length + j]*A->data[k*length + j];

      diag = A_kk - sum;

      /* check if the diag value is negative, but also not close to the minimum
         resolvable difference between to REAL8 numbers - if it is negative
         and also close to this value set it to +LAL_REAL8_EPS (see
         LALConstants.h), as it could be anywhere inbetween LAL_REAL8_MIN and
         LAL_REAL8_EPS */
      /* these problems are caused by the fact the when computing the eigen
         values/vectors to determine if the matrix is positive definite the
         process uses iterative methods to check on the convergence of values
         and these will only be good down to the precision of REAL8s */
      if( diag <= 0. && fabs(diag) >= LAL_REAL8_EPS && k != length-1 ){
        fprintf(stderr, "Error... input matrix is not pos def!\n");
        exit(0);
      }
      else if( diag < 0. && fabs(diag) <= LAL_REAL8_EPS ){
        /* diag = LAL_REAL8_EPS; */
        /* diag = 0.; */ /* set to zero as setting it to LAL_REAL8_EPS sometimes
                      is a value that's far larger than it should be */
        diag = fabs(diag);
      }
      else if( diag <= 0. && fabs(diag) >= LAL_REAL8_EPS && k == length-1 ){
        /* this is a kludge as a lot of the matricies seem to have entries
           where m(length, length) diagonal value is small but less than zero,
           so I'll just set it to zero manually */
        diag = 0.;
      }

      A->data[k*length + k] = sqrt(diag);

    }
  }

  /* set upper triangular matrix to zeros - for lower value */
  for(i=0; i<length; i++)
    for(j=i+1; j<length; j++)
      A->data[i*length + j] = 0.;

  /* check if the upper triangle is wanted - if so perform transpose */
  if(strstr(uOrl, "upper")!=NULL){
    REAL8 tempdata = 0.;

    /* perform transpose */
    for(j=0; j<length-1; j++){
      for(i=j; i<length; i++){
        tempdata = A->data[j*length + i];
        A->data[j*length + i] = A->data[i*length + j];
        A->data[i*length + j] = tempdata;
      }
    }
  }

  /* if( verbose ){
    fprintf(stderr, "\nCholesky decomposed matrix:\n");
    for(i=0; i<length; i++){
       for(j=0; j<length; j++)
         fprintf(stderr, "%.2e  ", A->data[i*length + j]);
       fprintf(stderr, "\n");
    }
  } */

  return A;
}



/* this function will draw a set of random numbers from a multivariate Gaussian
distribution, with a cholesky decomposed covariance matrix given by cholmat and
parameter mean values */
ParamData *multivariate_normal_deviates( REAL8Array *cholmat, ParamData *data,
  RandomParams *randomParams ){
  REAL4Vector *randNum=NULL;

  ParamData *deviates=NULL;

  INT4 dim=cholmat->dimLength->data[0]; /* covariance matrix dimensions */

  INT4 i=0, j=0;
  REAL8Vector *Z=NULL;

  /* check dimensions of covariance matrix and mean vector length are equal */
  if( cholmat->dimLength->data[0] != cholmat->dimLength->data[1] ){
    fprintf(stderr, "Error... wrong number of dimensions in input matrix!\n");
    exit(0);
  }

  deviates = XLALMalloc(MAXPARAMS*sizeof(ParamData));

  /* create a vector of random Gaussian numbers */
  randNum = XLALCreateREAL4Vector( dim );

  XLALNormalDeviates( randNum, randomParams );

  /* multiply cholsky decomposed matrix by randNum */
  Z = XLALCreateREAL8Vector( dim );
  for(i=0;i<dim;i++)
    Z->data[i] = 0.;

  for(i=0;i<dim;i++)
    for(j=0;j<dim;j++)
      Z->data[i] += cholmat->data[i*dim + j]*randNum->data[j];

  /* get the output random deviates by doing the mean plus Z */
  for(i=0;i<MAXPARAMS;i++){
    deviates[i].name = data[i].name;
    deviates[i].sigma = data[i].sigma;
    deviates[i].matPos = data[i].matPos;
    if( data[i].matPos != 0 )
      deviates[i].val = data[i].val + Z->data[data[i].matPos-1];
    else
      deviates[i].val = data[i].val;
  }

  XLALDestroyREAL4Vector( randNum );
  XLALDestroyREAL8Vector( Z );

  return deviates;
}



/* I need to define a standard set of positions in which various pulsar
   parameters will sit within the internal correlation matrix - this will as
   far as possible following the standard in the matrix files I have 
*/
/* function to read in the correlation matrix */
REAL8Array *read_correlation_matrix( CHAR *matrixFile, 
  BinaryPulsarParams params, ParamData *data ){
  FILE *fp=NULL;

  CHAR matrixParams[MAXPARAMS][6]; /* parameters in the correlation matrix */
  CHAR paramTmp[256];

  INT4 numParams=0, i=0, j=0, k=0, n=0;

  CHAR tmpStr[256], tmpStr2[256];

  INT4 arraySize=0;

  REAL8Array *corMat=NULL;
  UINT4Vector *matdims=NULL;

  INT4 DMpos=0, DM1pos=0; /* position of dispersion measure in matrix */
  INT4 numDM=0;
  REAL8 corTemp=0., junk=0.;

  int rc;

  ParamData paramData[]=
  {
    { "f0",  0., 0., 0 },{ "f1",  0., 0., 0 },{ "f2",  0., 0., 0 },
    { "Dec", 0., 0., 0 },{ "RA",  0., 0., 0 },{ "pmdc",0., 0., 0 },
    { "pmra",0., 0., 0 },{ "x",   0., 0., 0 },{ "e",   0., 0., 0 },
    { "T0",  0., 0., 0 },{ "Pb",  0., 0., 0 },{ "Om",  0., 0., 0 },
    { "Omdt",0., 0., 0 },{ "gamma",0., 0., 0 },{ "Pbdt",0.,0., 0 },
    { "s",   0., 0., 0 },{ "M",   0., 0., 0 },{ "m2",  0., 0., 0 },
    { "dth", 0., 0., 0 },{ "xdot",0., 0., 0 },{ "edot",0., 0., 0 },
    { "x2",  0., 0., 0 },{ "e2",  0., 0., 0 },{ "T02", 0., 0., 0 },
    { "Pb2", 0., 0., 0 },{ "Om2", 0., 0., 0 },{ "x3",  0., 0., 0 },
    { "e3",  0., 0., 0 },{ "T03", 0., 0., 0 },{ "Pb3", 0., 0., 0 },
    { "Om3", 0., 0., 0 },{ "Xpbd",0., 0., 0 },{ "f3",  0., 0., 0 },
    { "f4",  0., 0., 0 },{ "f5",  0., 0., 0 }
  };

  /* set the values */
  paramData[0].val = params.f0;      paramData[0].sigma = params.f0Err;
  paramData[1].val = params.f1;      paramData[1].sigma = params.f1Err;
  paramData[2].val = params.f2;      paramData[2].sigma = params.f2Err;
  paramData[3].val = params.dec;     paramData[3].sigma = params.decErr;
  paramData[4].val = params.ra;      paramData[4].sigma = params.raErr;
  paramData[5].val = params.pmdec;   paramData[5].sigma = params.pmdecErr;
  paramData[6].val = params.pmra;    paramData[6].sigma = params.pmraErr;
  paramData[7].val = params.x;       paramData[7].sigma = params.xErr;
  paramData[8].val = params.e;       paramData[8].sigma = params.eErr;
  paramData[9].val = params.T0;      paramData[9].sigma = params.T0Err;
  paramData[10].val = params.Pb;     paramData[10].sigma = params.PbErr;
  paramData[11].val = params.w0;     paramData[11].sigma = params.w0Err;
  paramData[12].val = params.wdot;   paramData[12].sigma = params.wdotErr;
  paramData[13].val = params.gamma;  paramData[13].sigma = params.gammaErr;
  paramData[14].val = params.Pbdot;  paramData[14].sigma = params.PbdotErr;
  paramData[15].val = params.s;      paramData[15].sigma = params.sErr;
  paramData[16].val = params.M;      paramData[16].sigma = params.MErr;
  paramData[17].val = params.m2;     paramData[17].sigma = params.m2Err;
  paramData[18].val = params.dth;    paramData[18].sigma = params.dthErr;
  paramData[19].val = params.xdot;   paramData[19].sigma = params.xdotErr;
  paramData[20].val = params.edot;   paramData[20].sigma = params.edotErr;
  paramData[21].val = params.x2;     paramData[21].sigma = params.x2Err;
  paramData[22].val = params.e2;     paramData[22].sigma = params.e2Err;
  paramData[23].val = params.T02;    paramData[23].sigma = params.T02Err;
  paramData[24].val = params.Pb2;    paramData[24].sigma = params.Pb2Err;
  paramData[25].val = params.w02;    paramData[25].sigma = params.w02Err;
  paramData[26].val = params.x3;     paramData[26].sigma = params.x3Err;
  paramData[27].val = params.e3;     paramData[27].sigma = params.e3Err;
  paramData[28].val = params.T03;    paramData[28].sigma = params.T03Err;
  paramData[29].val = params.Pb3;    paramData[29].sigma = params.Pb3Err;
  paramData[30].val = params.w03;    paramData[30].sigma = params.w03Err;
  paramData[31].val = params.xpbdot; paramData[31].sigma = params.xpbdotErr;
  paramData[32].val = params.f3;     paramData[32].sigma = params.f3Err;
  paramData[33].val = params.f4;     paramData[33].sigma = params.f4Err;
  paramData[34].val = params.f5;     paramData[34].sigma = params.f5Err;

  arraySize = MAXPARAMS;

  if( params.model != NULL && !strcmp(params.model, "ELL1") ){
    paramData[8].name = "eps1";
    paramData[8].val = params.eps1;
    paramData[8].sigma = params.eps1Err;

    paramData[11].name = "eps2";
    paramData[11].val = params.eps2;
    paramData[11].sigma = params.eps2Err;

    paramData[9].name = "Tasc";
    paramData[9].val = params.Tasc;
    paramData[9].sigma = params.TascErr;

    paramData[12].name = "e1dt";
    paramData[12].val = params.eps1dot;
    paramData[12].sigma = params.eps1dotErr;

    paramData[20].name = "e2dt";
    paramData[20].val = params.eps2dot;
    paramData[20].sigma = params.eps2dotErr;
  }

  /* if we have a correlation matrix file then read it in */
  if( matrixFile != NULL ){
    /* read in data from correlation matrix file */
    if((fp = fopen(matrixFile, "r")) == NULL){
      fprintf(stderr, "Error...No correlation matrix file!\n" );
      return NULL;
    }

    /* read in the first line of the matrix file */
    while(fscanf(fp, "%s", paramTmp)){
      if(strchr(paramTmp, '-') != NULL)
        break;

      if(feof(fp)){
        fprintf(stderr, "Error... I've reached the end of the file without \
reading any correlation data!");
        fclose(fp);
        exit(0);
      }

      sprintf(matrixParams[numParams], "%s", paramTmp);
      numParams++;

      /*check if parameter is actually for a dispersion measure (ignore if so)*/
      if(!strcmp(paramTmp, "DM")){
        numParams--;
        DMpos = i;
        numDM++;
      }
      if(!strcmp(paramTmp, "DM1")){
        numParams--;
        DM1pos = i;
        numDM++;
      }

      i++;
    };

    if(numParams > arraySize){
      fprintf(stderr, "More parameters in matrix file than there should be!\n");
      exit(0);
    }

    matdims = XLALCreateUINT4Vector( 2 );
    matdims->data[0] = numParams;
    matdims->data[1] = numParams;

    corMat = XLALCreateREAL8Array( matdims );

    /* find positions of each parameter */
    /* the strings that represent parameters in a matrix are given in the param
       variable in the tempo code mxprt.f */
    /* read in matrix */
    k=0;
    for(i=0;i<numParams+numDM;i++){
      n=0;
      rc = fscanf(fp, "%s%s", tmpStr, tmpStr2);

      /* if its a dispersion measure then just skip the line */
      if( (DMpos != 0 && i == DMpos) || (DM1pos != 0 && i == DM1pos) ){
        rc = fscanf(fp, "%*[^\n]");
        k--;
        continue;
      }

      for(j=0;j<i+1;j++){
        if( (DMpos != 0 && j == DMpos) || (DM1pos != 0 && j == DM1pos) ){
          rc = fscanf(fp, "%lf", &junk);
          n--;
          continue;
        }

        rc = fscanf(fp, "%lf", &corTemp);

        /* if covariance equals 1 set as 0.9999999, because values of 1
           can cause problems of giving singular matrices */
        if( (n != k) && (corTemp == 1.) )
          corTemp = 0.9999999;
        else if( (n != k) && (corTemp == -1.) )
          corTemp = -0.9999999;

        corMat->data[k*corMat->dimLength->data[0] + n] = corTemp;

        if(n != k)
          corMat->data[n*corMat->dimLength->data[0] + k] = corTemp;

        n++;
      }

      /* send an error if we hit the end of the file */
      if( feof(fp) || rc == EOF ){
        fprintf(stderr, "Error reading in matrix - hit end of file!\n");
        exit(0);
      }

      k++;
    }

    fclose(fp);

    /* give the correlation matrix positions of the parameters */
    for(i=1;i<numParams+1;i++){
      for(j=0;j<arraySize;j++){
        if(!strcmp(matrixParams[i-1], paramData[j].name)){
          paramData[j].matPos = i;
          break;
        }
      }
    }
  }
  else{ /* create files with just variances from par file and no correlations */
    j = 0;
    for( i = 0; i < MAXPARAMS; i++ ){
      if( paramData[i].sigma != 0.0 ){
        j++;
        paramData[i].matPos = j;
      }
    }
    
    /* create array */
    matdims = XLALCreateUINT4Vector( 2 );
    matdims->data[0] = j;
    matdims->data[1] = j;

    corMat = XLALCreateREAL8Array( matdims );
    
    /* set diagonal elements to one - they'll be converted to variances later */
    for( i = 0; i < j; i++ ){
      for ( k = 0; k < j; k++){
        if ( i == k )
          corMat->data[i*corMat->dimLength->data[0]+k] = 1.;
        else
          corMat->data[i*corMat->dimLength->data[0]+k] = 0.;
      }
    }
  }
  
  /* pass the parameter data to be output */
  memcpy(data, paramData, sizeof(paramData));

  XLALDestroyUINT4Vector( matdims );

  return corMat;
}



/* function to turn the input /correlation/ matrix into a covariance matrix */
REAL8Array *create_covariance_matrix( ParamData *data, REAL8Array *corMat, 
  INT4 isinv ){
  REAL8Array *covMat=NULL;
  INT4 i=0, j=0;

  covMat = XLALCreateREAL8Array( corMat->dimLength );

  /* convert correlation matrix into a covariance matrix */
  for(i=0;i<MAXPARAMS;i++){
    if( data[i].matPos != 0 ){
      for(j=0;j<MAXPARAMS;j++){
        if( data[j].matPos != 0 ){
          /* if covariance matrix then isinv = 0 */
          if( isinv == 0 ){
            covMat->data[(data[i].matPos-1)*covMat->dimLength->data[0] +
              data[j].matPos-1] =
              corMat->data[(data[i].matPos-1)*corMat->dimLength->data[0] +
              data[j].matPos-1] * data[i].sigma * data[j].sigma;
          }
          else if( isinv == 1 ){ /* doing matrix inverse */
            covMat->data[(data[i].matPos-1)*covMat->dimLength->data[0] +
              data[j].matPos-1] =
              corMat->data[(data[i].matPos-1)*corMat->dimLength->data[0] + 
              data[j].matPos-1] / ( data[i].sigma * data[j].sigma );
          }
          else{
            fprintf(stderr, "Error... in setting covariance matrix isinv must \
be 0 or 1\n");
            exit(0);
          } 
        }
      }
    }
  }

  return covMat;
}



/* function to check that a 2D matrix is positive definite - if not 
   positive definite it will be converted so that it is */
REAL8Array *check_positive_definite( REAL8Array *matrix ){
  static LALStatus status;

  REAL8Vector *eigenval=NULL;
  REAL8Array *eigenvec=NULL;

  REAL8Array *posdef=NULL;

  INT4 i=0, j=0;

  /* copy input array into eigenvec as this gets converted by function */
  eigenvec = XLALCreateREAL8Array( matrix->dimLength );

  for( i=0; i<(INT4)eigenvec->dimLength->data[0]; i++ ){
    for( j=0; j<(INT4)eigenvec->dimLength->data[1]; j++ ){
      eigenvec->data[i*eigenvec->dimLength->data[0] + j] =
        matrix->data[i*matrix->dimLength->data[0] + j];
    }
  }

  eigenval = XLALCreateREAL8Vector( matrix->dimLength->data[0] );

  /* calculate the eigen values and vectors */
  LAL_CALL( LALDSymmetricEigenVectors( &status, eigenval, eigenvec ), &status );

  for( i=0; i<(INT4)matrix->dimLength->data[0]; i++ ){
    /* first check if any eigen values are zero and if so convert to positive
       definite matrix */
    if( eigenval->data[i] < 0. && fabs(eigenval->data[i]) > 10.*LAL_REAL8_EPS ){
      fprintf(stderr, "Eigenvalue is negative. Non-postive definite matrix!\n");
      posdef = convert_to_positive_definite( matrix );
      break;
    }
  }

  /* if matrix is positive definite return NULL it i.e. posdef hasn't been set */
  if( posdef == NULL ){
    XLALDestroyREAL8Array( eigenvec );
    XLALDestroyREAL8Vector( eigenval );
    return NULL;
  }

  /* re-check new matrix for positive definiteness, but be aware of values
     close to the precision of REAL8 numbers */
  for( i=0; i<(INT4)eigenvec->dimLength->data[0]; i++ ){
    for( j=0; j<(INT4)eigenvec->dimLength->data[1]; j++ ){
      eigenvec->data[i*eigenvec->dimLength->data[0] + j] =
        posdef->data[i*posdef->dimLength->data[0] + j];
    }
    eigenval->data[i] = 0.;
  }

  LAL_CALL( LALDSymmetricEigenVectors( &status, eigenval, eigenvec ), &status );

  for( i=0; i<(INT4)matrix->dimLength->data[0]; i++ ){
    if( eigenval->data[i] < 0. && fabs(eigenval->data[i]) > 10.*LAL_REAL8_EPS){
      fprintf(stderr, "ABORT! Eigenvalue is negative. Non-postive definite \
matrix!\n");
      exit(0);
    }
  }

  XLALDestroyREAL8Array( eigenvec );
  XLALDestroyREAL8Vector( eigenval );

  return posdef;
}



/* this function takes a matrix that isn't positive definite and converts it
into a positive definite matrix using the method (number 2) of Rebonato and
Jackel (see their paper at
http://www.riccardorebonato.co.uk/papers/ValCorMat.pdf) */
REAL8Array *convert_to_positive_definite( REAL8Array *nonposdef ){
  static LALStatus status;

  REAL8Vector *eigenval=NULL;
  REAL8Array *eigenvec=NULL;

  REAL8Array *posdef = NULL; /* output positive definite matrix */

  REAL8Array *Lprime=NULL;
  REAL8Array *B=NULL, *Bprime=NULL, *Btrans=NULL;
  REAL8Array *T=NULL;

  REAL8 Tval=0.;

  INT4 i=0, j=0, length=0;

  if( verbose ) fprintf(stderr, "Converting to positive definite matrix\n");

  /* check that matrix is square */
  if( nonposdef->dimLength->data[0] != nonposdef->dimLength->data[1] ){
    fprintf(stderr, "Error... matrix must be square!\n");
    exit(0);
  }
  else
    length = nonposdef->dimLength->data[0];

  Lprime = XLALCreateREAL8Array( nonposdef->dimLength );
  T = XLALCreateREAL8Array( nonposdef->dimLength );

  /* copy input array into eigenvec as this gets converted by function */
  eigenvec = XLALCreateREAL8Array( nonposdef->dimLength );

  for( i=0; i<length; i++ ){
    for( j=0; j<length; j++ ){
      eigenvec->data[i*length + j] = nonposdef->data[i*length + j];

      /* initialise Lprime and T to zeros */
      Lprime->data[i*length + j] = 0.;
      T->data[i*length + j] =  0.;
    }
  }

  eigenval = XLALCreateREAL8Vector( length );

  /* calculate the eigen values and vectors */
  LAL_CALL( LALDSymmetricEigenVectors( &status, eigenval, eigenvec ), &status );

  /* if eigen value is > 0 set Lprime to that value i.e. have eigen values of 
     zero if eigen value is negative */
  for( i=0; i<length; i++ )
    if( eigenval->data[i] > 0. )
      Lprime->data[i*length + i] = eigenval->data[i];

  /* compute scaling matrix T */
  for( i=0; i<length; i++ ){
    Tval = 0.;
    for( j=0; j<length; j++ ){
      Tval += eigenvec->data[i*length + j] * eigenvec->data[i*length + j] *
         Lprime->data[j*length + j];
    }

    Tval = 1./Tval;

    /* really we just want the sqrt of T */
    T->data[i*length + i] = sqrt(Tval);
  }

  /* convert Lprime to sqrt(lambdaprime) */
  for( i=0; i<length; i++ ){
    REAL8 tempL = Lprime->data[i*length + i];

    Lprime->data[i*length + i] = sqrt(tempL);
  }

  /* Bprime = S*sqrt(lambdaprime); */
  Bprime = XLALCreateREAL8Array( nonposdef->dimLength );
  LAL_CALL( LALDMatrixMultiply( &status, Bprime, eigenvec, Lprime ), &status );

  /* B = sqrt(T)*Bprime */
  B = XLALCreateREAL8Array( nonposdef->dimLength );
  LAL_CALL( LALDMatrixMultiply( &status, B, T, Bprime ), &status );

  /* transpose(B) */
  Btrans = XLALCreateREAL8Array( nonposdef->dimLength );
  LAL_CALL( LALDMatrixTranspose( &status, Btrans, B ), &status );

  /* posdef = B*transpose(B) this is our new positive definite matrix */
  posdef = XLALCreateREAL8Array( nonposdef->dimLength );
  LAL_CALL( LALDMatrixMultiply( &status, posdef, B, Btrans ), &status );

  XLALDestroyREAL8Array( eigenvec );
  XLALDestroyREAL8Vector( eigenval );
  XLALDestroyREAL8Array( Lprime );
  XLALDestroyREAL8Array( T );
  XLALDestroyREAL8Array( Bprime );
  XLALDestroyREAL8Array( B );
  XLALDestroyREAL8Array( Btrans );

  /* check that the difference between new and old values are greater than the
     maximum precision between REAL8 values (LAL_REAL8_EPS) - if not use
     original value */
  for( i=0; i<length; i++ ){
    for( j=0; j<length; j++ ){
      if( fabs(posdef->data[i*length + j]-nonposdef->data[i*length + j]) <=
        LAL_REAL8_EPS)
        posdef->data[i*length + j] = nonposdef->data[i*length + j];
    }
  }

  return posdef;
}

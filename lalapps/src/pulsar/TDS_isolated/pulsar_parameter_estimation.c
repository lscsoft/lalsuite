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

#include "pulsar_parameter_estimation.h"

RCSID("$Id$");

/* global variable */
INT4 verbose=0;

INT4 main(INT4 argc, CHAR *argv[]){
  static LALStatus status;
  
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
  
  OutputParams output;
  REAL8 maxPost=0.;
  REAL8 logNoiseEv[5]; /* log evidence for noise only (no signal) */
  Results results;
  REAL8 h0ul=0.;  

  CHAR *params[]={"h0", "phi", "psi", "ciota"};
  
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
    if( inputs.matrixFile != NULL ){
      edat = XLALMalloc(sizeof(*edat));

      (*edat).ephiles.earthEphemeris = inputs.earthfile;
      (*edat).ephiles.sunEphemeris = inputs.sunfile;
      LAL_CALL( LALInitBarycenter(&status, edat), &status );
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
      /* exclude values smaller than 1e-28 as most are spurious points caused
         during a the heterodyne stage (e.g. when frame files were missing in
         the original S5 analysis) */
      if( fabs(dataVals.re) > 1.e-28 && fabs(dataVals.im) > 1.e-28 ){
        data[k].times->data[j] = times;
        data[k].data->data[j] = dataVals;
        
        /* get the power from the time series */
        stdh0 += dataVals.re*dataVals.re + dataVals.im*dataVals.im;
        
        j++;
      }
    }
    
    fclose(fp);
    
    if( verbose )
      fprintf(stderr, "I've read in the data for %s.\n", dets[i]);
    
    data[k].data = XLALResizeCOMPLEX16Vector(data[k].data, j);
    data[k].times = XLALResizeREAL8Vector(data[k].times, j);
    
    /* if there is no input range for h0 then estimate it from the data */
    /* only do this once if performing grid search, but do for each seperate
       data set if doing MCMC */
    if( ( inputs.mesh.maxVals.h0 == 0 || inputs.mcmc.sigmas.h0 == 0 ) && (
          inputs.mcmc.doMCMC != 0 || i == 0 ) ){      
      if( verbose ) fprintf(stderr, "Calculating h0 UL estimate: ");
      
      /* get the power spectral density power/bandwidth (1/60 Hz) */
      stdh0 = stdh0/((REAL8)j*(1./60.));
     
      /* upper limit estimate comes from ~ h0 = 10.8*sqrt(Sn/T) */
      stdh0 = 10.8*sqrt(stdh0/((REAL8)j*60.));
      
      /* set h0 max value for the grid at 5 times the expected ul */
      if( inputs.mesh.maxVals.h0 == 0 ){
        inputs.mesh.maxVals.h0 = 5.*stdh0;
        inputs.mesh.delta.h0 = (inputs.mesh.maxVals.h0 -
          inputs.mesh.minVals.h0)/(REAL8)(inputs.mesh.h0Steps - 1.);
      }
      
      /* set the MCMC h0 proposal step size at stdh0/20 */
      if( inputs.mcmc.doMCMC != 0 )
        inputs.mcmc.sigmas.h0 = stdh0/20.;
      
      if( verbose ) fprintf(stderr, "%le\n", stdh0);
    }

    /*========================================================================*/
    
    output.det = dets[i];
    
    fprintf(stderr, "allocating memory for look-up table\n");
    /* create lookup table */
    data[k].lookupTable = NULL;
    data[k].lookupTable = XLALCalloc(1, sizeof(DetRespLookupTable));
    data[k].lookupTable->lookupTable=NULL;
    detAndSource.pSource = &inputs.psr;
    detAndSource.pDetector = &detPos[i];
    
    /* create memory for the lookup table */
    data[k].lookupTable->lookupTable = XLALCalloc(inputs.mesh.psiRangeSteps, 
      sizeof(LALDetAMResponse *));
    fprintf(stderr, "still allocating memory for look-up table\n");
    for( j = 0 ; j < inputs.mesh.psiRangeSteps ; j++ ){
      data[k].lookupTable->lookupTable[j] =
        XLALCalloc(inputs.mesh.timeRangeSteps, sizeof(LALDetAMResponse));
    }

    data[k].lookupTable->psiSteps = inputs.mesh.psiRangeSteps;
    data[k].lookupTable->timeSteps = inputs.mesh.timeRangeSteps;
    
    fprintf(stderr, "making look-up table\n");
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
    if( inputs.mcmc.doMCMC != 0 )
      perform_mcmc(&data[k], inputs, 1, dets[i], &detPos[i], edat);
              
    /*========================================================================*/
  }

  /*=================== CREATE THE JOINT POSTERIOR IF REQUIRED ===============*/
  if( numDets > 1 ){
    output.det = "Joint";
    
    if( inputs.mcmc.doMCMC == 0 ){
      REAL8 totLogNoiseEv=0.;

      for( n = 0 ; n < numDets ; n++ ) totLogNoiseEv += logNoiseEv[n];

      output.outPost = inputs.outputPost; /* set for whether we want to output
                                            the full posterior */
   
      maxPost = log_posterior(jointLike, inputs.priors, inputs.mesh, output);
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
        fprintf(fp, "%.1lf%% h0 upper limit = %lf\n", output.dob,
          results.h0UpperLimit);
      }
      fclose(fp);
    }
    /*========================================================================*/
  
    /*======================= PERFORM JOINT MCMC =============================*/
    if( inputs.mcmc.doMCMC != 0 ){
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
          free(jointLike[i][j]);
        
        free(singleLike[i][j]);
      }
      if( numDets > 1 )
        free(jointLike[i]);
      
      free(singleLike[i]);
    }
  
    if( numDets > 1 )
      free(jointLike);
    
    free(singleLike);
  }
  /*=========================================================================*/ 
                                            
  return 0;
}


/* function to get the input arguments from the command line */
void get_input_args(InputParams *inputParams, INT4 argc, CHAR *argv[]){
  struct option long_options[] = 
  {
    { "help",         no_argument,       0, 'h' },
    { "verbose",      no_argument, &inputParams->verbose, 1 },
    { "detectors",    required_argument, 0, 'D' },
    { "pulsar",       required_argument, 0, 'p' },
    { "par-file",     required_argument, 0, 'P' },
    { "input-dir",    required_argument, 0, 'i' },
    { "output-dir",   required_argument, 0, 'o' },
    { "minh0",        required_argument, 0, 'a' },
    { "maxh0",        required_argument, 0, 'A' },
    { "h0steps",      required_argument, 0, 'j' },
    { "minphi0",      required_argument, 0, 'b' },
    { "maxphi0",      required_argument, 0, 'B' },
    { "phi0steps",    required_argument, 0, 'k' },
    { "minpsi",       required_argument, 0, 's' },
    { "maxpsi",       required_argument, 0, 'S' },
    { "psisteps",     required_argument, 0, 'm' },
    { "minci",        required_argument, 0, 'c' },
    { "maxci",        required_argument, 0, 'C' },
    { "cisteps",      required_argument, 0, 'n' },
    { "psi-bins",     required_argument, 0, 'l' },
    { "time-bins",    required_argument, 0, 'L' },
    { "h0prior",      required_argument, 0, 'q' },
    { "phi0prior",    required_argument, 0, 'Q' },
    { "psiprior",     required_argument, 0, 'U' },
    { "ciprior",      required_argument, 0, 'u' },
    { "h0mean",       required_argument, 0, 'Y' },
    { "h0sig",        required_argument, 0, 'T' },
    { "phi0mean",     required_argument, 0, 'v' },
    { "phi0sig",      required_argument, 0, 'V' },
    { "psimean",      required_argument, 0, 'z' },
    { "psisig",       required_argument, 0, 'Z' },
    { "cimean",       required_argument, 0, 'e' },
    { "cisig",        required_argument, 0, 'E' },
    { "output-post",  no_argument, &inputParams->outputPost, 1 },
    { "dob-ul",       required_argument, 0, 'd' },
    { "mcmc",         no_argument, &inputParams->mcmc.doMCMC, 1 },
    { "iterations",   required_argument, 0, 'I' },
    { "burn-in",      required_argument, 0, 'x' },
    { "temperature",  required_argument, 0, 't' },
    { "h0-width",     required_argument, 0, 'H' },
    { "psi-width",    required_argument, 0, 'w' },
    { "phi0-width",   required_argument, 0, 'W' },
    { "ci-width",     required_argument, 0, 'y' },
    { "glitch-times", required_argument, 0, 'g' },
    { "glitch-cut",   required_argument, 0, 'G' },
    { "chunk-min",    required_argument, 0, 'K' },
    { "chunk-max",    required_argument, 0, 'N' },
    { "output-rate",  required_argument, 0, 'X' },
    { "nglitch",      required_argument, 0, 'O' },
    { "earth-ephem",  required_argument, 0, 'J' },
    { "sun-ephem",    required_argument, 0, 'M' },
    { "covariance",   required_argument, 0, 'r' },
    { "use-priors",   no_argument, &inputParams->usepriors, 1 },
    { 0, 0, 0, 0 }
  };
  
  CHAR args[] =
"hD:p:P:i:o:a:A:j:b:B:k:s:S:m:c:C:n:l:L:q:Q:U:u:Y:T:v:V:z:Z:e:E:d:I:x:t:H:w:W:\
y:g:G:K:N:X:O:J:M:r:" ;
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
  inputParams->priors.h0Prior = "uniform";
  inputParams->priors.phiPrior = "uniform";
  inputParams->priors.psiPrior = "uniform";
  inputParams->priors.ciotaPrior = "uniform";
  
  /* default MCMC parameters */
  inputParams->mcmc.sigmas.h0 = 0.;           /* estimate from data */
  inputParams->mcmc.sigmas.phi0 = LAL_PI/10.; /* 20th of phi range */
  inputParams->mcmc.sigmas.psi = LAL_PI/40.;  /* 20th of psi range */
  inputParams->mcmc.sigmas.ci = 0.1;          /* 20th of cosi range */
  inputParams->mcmc.outputRate = 1;           /* output every sample */
  
  inputParams->mcmc.iterations = 10000;       /* default 10000 points */
  inputParams->mcmc.temperature = 1;          /* default annealing */
  inputParams->mcmc.burnIn = 1000;            /* default burn in time */
  
  inputParams->mcmc.nGlitches = 0;            /* no glitches is default */
  
  inputParams->matrixFile = NULL;             /* no covriance file */

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
        fprintf(stderr, USAGE, program);
        exit(0);
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
      case 'Q': /* prior on phi0 */
        inputParams->priors.phiPrior = optarg;
        break;
      case 'U': /* prior on psi */
        inputParams->priors.psiPrior = optarg;
        break;
      case 'u': /* prior on cos(iota) */
        inputParams->priors.ciotaPrior = optarg;
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
      case 'e': /* mean of cos(iota) Gaussian prior */
        inputParams->priors.meanciota = atof(optarg);
        break;
      case 'E': /* standard deviation of Gaussian cos(iota) prior */
        inputParams->priors.stdciota = atof(optarg);
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
      case 'r':
        inputParams->matrixFile = optarg;
        break;
      case '?':
        fprintf(stderr, "Unknown error while parsing options\n");
      default:
        fprintf(stderr, "Unknown error while parsing options\n");
    }
  }
  
  /* check parameters for wierd values */
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
}


/* function to allocate memory for a likelihood array - it will be index as
   logLike[phi][cosiota][psi][h0] */
REAL8 ****allocate_likelihood_memory(MeshGrid mesh){
  INT4 i=0, j=0, k=0;
  REAL8 ****logLike=NULL;
  
  /* allocate the h0 positions using calloc (i.e. array will be initialise to 
     zero */
  logLike = calloc(mesh.phiSteps, sizeof(REAL8 ***));

  for( i = 0 ; i < mesh.phiSteps ; i++ ){
    logLike[i] = calloc(mesh.ciotaSteps, sizeof(REAL8 **));
    
    for( j = 0 ; j < mesh.ciotaSteps ; j++ ){
      logLike[i][j] = calloc(mesh.psiSteps, sizeof(REAL8 *));
      
      for( k = 0 ; k < mesh.psiSteps ; k++ ){
        logLike[i][j][k] = calloc(mesh.h0Steps, sizeof(REAL8));
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
  
  REAL8 cosphi=0., sinphi=0.;

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
    cosphi = cos(vars.phi0);
    sinphi = sin(vars.phi0);
    
    for( j = 0 ; j < mesh.ciotaSteps ; j++ ){
      vars.ci = mesh.minVals.ci + (REAL8)j*mesh.delta.ci;
      vars.Xplus = 0.5*(1.+vars.ci*vars.ci);
      vars.Xcross = vars.ci;
      vars.Xpsinphi_2 = 0.5*vars.Xplus*sinphi;
      vars.Xcsinphi_2 = 0.5*vars.Xcross*sinphi;
      vars.Xpcosphi_2 = 0.5*vars.Xplus*cosphi;
      vars.Xccosphi_2 = 0.5*vars.Xcross*cosphi;

      for( k = 0 ; k < mesh.psiSteps ; k++ ){
        /* set up these parameters as they need to be passed to log_likelihood
           function - although tey aren't used in the grid-based search */
        BinaryPulsarParams params[2]; 
        BarycenterInput barys[2]; 
        EphemerisData *edat=NULL;

        vars.psi = mesh.minVals.psi + (REAL8)k*mesh.delta.psi;
        
        /* perform final loop over h0 within log_likelihood function */
        noiseEvidence = log_likelihood(logLike[i][j][k], data, vars, mesh,
params, barys, edat);
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

  length = (INT4)data.data->length + 1 - 
           data.chunkLengths->data[(INT4)data.chunkLengths->length-1];
           
  for( i = 0 ; i < length ; i+= chunkLength ){
    chunkLength = data.chunkLengths->data[count];
    data.sumData->data[count] = 0.;

    for( j = i ; j < i + chunkLength ; j++){
      B.re = data.data->data[j].re;
      B.im = data.data->data[j].im;

      /* sum up the data */
      data.sumData->data[count] += (REAL8)(B.re*B.re + B.im*B.im);
    }
    
    count++;
  }
}


REAL8 log_likelihood( REAL8 *likeArray, DataStructure data,
  IntrinsicPulsarVariables vars, MeshGrid mesh,
  BinaryPulsarParams params[2], BarycenterInput barys[2], EphemerisData *edat ){
  static LALStatus status;

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
  REAL8 log2=log(2.);
  
  REAL8 noiseEvidence=0.; /* the log evidence that the data is just noise */

  INT4 first=0, through=0;
  
  /* variables for delta phi caused by parameter uncertainties */
  EmissionTime emit1, emit2;
  EarthState earth;
  REAL8 phi1=0., phi2=0., dphi=0.;

  BinaryPulsarInput binInput1, binInput2;
  BinaryPulsarOutput binOutput1, binOutput2;
  REAL8 T01=0., T02=0., dt1=0., dt2=0.;

  /* to save time get all log factorials up to chunkMax */
  for( i = 0 ; i < data.chunkMax+1 ; i++ )
    exclamation[i] = log_factorial(i);
  
  /* set the psi bin for the lookup table */
  psibin = (INT4)ROUND( ( vars.psi + LAL_PI/4. )
           *(REAL8)(data.lookupTable->psiSteps-1.)/LAL_PI_2 );
           
  length = (INT4)data.data->length + 1 - 
           data.chunkLengths->data[(INT4)data.chunkLengths->length-1];
           
  tstart = data.times->data[0]; /* time of first B_k */

  /* set the position and frequency epochs if not already set */
  if(params[0].pepoch == 0. && params[0].posepoch != 0.)
    params[0].pepoch = params[0].posepoch;
  else if(params[0].posepoch == 0. && params[0].pepoch != 0.)
    params[0].posepoch = params[0].pepoch;

  if(params[1].pepoch == 0. && params[1].posepoch != 0.)
    params[1].pepoch = params[1].posepoch;
  else if(params[1].posepoch == 0. && params[1].pepoch != 0.)
    params[1].posepoch = params[1].pepoch;

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
      T = fmod(data.times->data[j] - tstart, 86400.);
      timebin = (INT4)fmod(ROUND(T*(REAL8)data.lookupTable->timeSteps/86400.),
        data.lookupTable->timeSteps);
      
      plus = data.lookupTable->lookupTable[psibin][timebin].plus;
      cross = data.lookupTable->lookupTable[psibin][timebin].cross;
      
      B.re = data.data->data[j].re;
      B.im = data.data->data[j].im;

      /*********************************************************/
      /* stuff for phase offset due to parameter uncertainties - MCMC only */
      if( edat != NULL ){
        if(data.times->data[j] <= 820108813)
          (*edat).leap = 13;
        else
          (*edat).leap = 14;

        T01 = params[0].pepoch;
        T02 = params[1].pepoch;

        barys[0].tgps.gpsSeconds = (UINT8)floor(data.times->data[j]);
        barys[0].tgps.gpsNanoSeconds =
          (UINT8)floor((fmod(data.times->data[j],1.0)*1.e9));
      
        barys[0].delta = params[0].dec + (data.times->data[j] -
          T01)*params[0].pmdec;
        barys[0].alpha = params[0].ra + (data.times->data[j] -
          T01)*params[0].pmra/cos(barys[0].delta);

        LAL_CALL( LALBarycenterEarth(&status, &earth, &barys[0].tgps, edat),
          &status ); 
        LAL_CALL( LALBarycenter(&status, &emit1, &barys[0], &earth), &status );

        /* if sky positions aren't the same then work out new one */
        if( params[0].ra != params[1].ra || params[0].dec != params[1].dec ||
            params[0].pmra != params[1].pmra || params[0].pmdec !=
            params[1].pmdec ){
          barys[1].tgps.gpsSeconds = (UINT8)floor(data.times->data[j]);
          barys[1].tgps.gpsNanoSeconds =
            (UINT8)floor((fmod(data.times->data[j],1.0)*1.e9));

          barys[1].delta = params[1].dec + (data.times->data[j] -
            T01)*params[1].pmdec;
          barys[1].alpha = params[1].ra + (data.times->data[j] -
            T01)*params[1].pmra/cos(barys[1].delta);

          LAL_CALL( LALBarycenterEarth(&status, &earth, &barys[1].tgps, edat),
            &status ); 
          LAL_CALL(LALBarycenter(&status, &emit2, &barys[1], &earth), &status);
        }
        else
          emit2.deltaT = emit1.deltaT;
      
        if(params[0].model!=NULL){
          binInput1.tb = data.times->data[j] + emit1.deltaT;
          binInput2.tb = data.times->data[j] + emit2.deltaT;

          XLALBinaryPulsarDeltaT(&binOutput1, &binInput1, &params[0]);
          XLALBinaryPulsarDeltaT(&binOutput2, &binInput2, &params[1]);
          
          dt1 = (data.times->data[j] - T01) + emit1.deltaT + binOutput1.deltaT;
          dt2 = (data.times->data[j] - T02) + emit2.deltaT + binOutput2.deltaT;
        }
        else{
          dt1 = (data.times->data[j] - T01) + emit1.deltaT;
          dt2 = (data.times->data[j] - T02) + emit2.deltaT;
        }

        /* work out phase */
        phi1 = 2.*(params[0].f0*dt1 + 0.5*params[0].f1*dt1*dt1 +
  (1./6.)*params[0].f2*dt1*dt1*dt1 + (1./24.)*params[0].f3*dt1*dt1*dt1*dt1);

        phi2 = 2.*(params[1].f0*dt2 + 0.5*params[1].f1*dt2*dt2 +
  (1./6.)*params[1].f2*dt2*dt2*dt2 + (1./24.)*params[1].f3*dt2*dt2*dt2*dt2);
      
        dphi = phi2 - phi1;

        /* create the signal model */
        model.re = (plus*vars.Xpcosphi_2 + cross*vars.Xcsinphi_2)*cos(-dphi) +
                 (cross*vars.Xccosphi_2 - plus*vars.Xpsinphi_2)*sin(-dphi);
        model.im = (plus*vars.Xpsinphi_2 - cross*vars.Xccosphi_2)*cos(-dphi) +
                 (cross*vars.Xcsinphi_2 + plus*vars.Xpcosphi_2)*sin(-dphi);
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
         logL = (m-1)log(2) + m! - m*sum((Bk - yk)^2) */
         
      /* reset array if first time in loop - else joint likelihoods
         will be wrong */
      if( first == 0 && through == 1 ) 
        likeArray[k] = (chunkLength - 1.)*log2;
      else likeArray[k] += (chunkLength - 1.)*log2;
        
      likeArray[k] += exclamation[(INT4)chunkLength];
      likeArray[k] -= chunkLength*log(chiSquare);

      /* get the log evidence for the data not containing a signal */
      if( k == 0 ){
        noiseEvidence += (chunkLength - 1.)*log2;
        noiseEvidence += exclamation[(INT4)chunkLength];
        noiseEvidence -= chunkLength*log(data.sumData->data[count]);
      }
    }
    
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


/* calculate the log prior */
REAL8 log_prior(PriorVals prior, MeshGrid mesh){
  REAL8 pri=1.;
  
  /* FIXME: Add ability to read in a old pdf file to use as a prior */
  
  if(strcmp(prior.h0Prior, "uniform") == 0){
    pri *= 1.; /* set h0 prior to be one for all values if uniform */
    /* pri *= 1./(mesh.maxVals.h0 - mesh.minVals.h0); */
  }
  else if(strcmp(prior.h0Prior, "jeffreys") == 0)
    pri *= 1./prior.vars.h0;
  else if(strcmp(prior.h0Prior, "gaussian") == 0){
    pri *= (1./(prior.stdh0*sqrt(2.*LAL_PI)))*exp(-(prior.vars.h0 -
prior.meanh0)*(prior.vars.h0 -
    prior.meanh0)/(2.*prior.stdh0*prior.stdh0));
  }
  
  if(strcmp(prior.phiPrior, "uniform") == 0){
    pri *= 1./(mesh.maxVals.phi0 - mesh.minVals.phi0);
  }
  /* wrap around Gaussian priors so that the mean is always at the centre */
  else if(strcmp(prior.phiPrior, "gaussian") == 0){
    if( prior.vars.phi0 < prior.meanphi - LAL_PI )
      prior.vars.phi0 += LAL_TWOPI;
    else if( prior.meanphi + LAL_PI < prior.vars.phi0 )
      prior.vars.phi0 -= LAL_TWOPI;
    
    pri *= ( 1./ (prior.stdphi*sqrt(2.*LAL_PI) ) )*exp( -( prior.vars.phi0 -
            prior.meanphi ) * ( prior.vars.phi0 - prior.meanphi ) /
           ( 2.*prior.stdphi*prior.stdphi ) );
  }
  
  if(strcmp(prior.psiPrior, "uniform") == 0){
    pri *= 1./(mesh.maxVals.psi - mesh.minVals.psi);
  }
  else if(strcmp(prior.psiPrior, "gaussian") == 0){
    if( prior.vars.psi < prior.meanpsi - LAL_PI/4. )
      prior.vars.psi += LAL_PI/2.;
    else if( prior.meanpsi + LAL_PI/4. < prior.vars.psi )
      prior.vars.psi -= LAL_PI/2.;
    
    pri *= ( 1./ (prior.stdpsi*sqrt(2.*LAL_PI) ) )*exp( -( prior.vars.psi -
            prior.meanpsi ) * ( prior.vars.psi - prior.meanpsi ) / 
           ( 2.*prior.stdpsi*prior.stdpsi ) );
  }
  
  if(strcmp(prior.ciotaPrior, "uniform") == 0){
    pri *= 1./(mesh.maxVals.ci - mesh.minVals.ci);
  }
  else if(strcmp(prior.ciotaPrior, "gaussian") == 0){
    if( prior.vars.ci < prior.meanciota - 1. )
      prior.vars.ci += 2.;
    else if( prior.meanciota + 1. < prior.vars.ci )
      prior.vars.ci -= 2.;
    
    pri *= ( 1./ (prior.stdciota*sqrt(2.*LAL_PI) ) ) * exp( -( prior.vars.ci -
            prior.meanciota ) * ( prior.vars.ci - prior.meanciota ) / 
           ( 2.*prior.stdciota*prior.stdciota ) );
  }

  return log(pri);
}


/* calculate the unnormalised log posterior and output the max value - print
out the log posterior if requested */
REAL8 log_posterior(REAL8 ****logLike, PriorVals prior, MeshGrid mesh,
  OutputParams output){
  REAL8 maxPost=-1.e200, mP=0.;
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


/* marginalise posterior over requested parameter and output the log evidence if
requested */
Results marginalise_posterior(REAL8 ****logPost, MeshGrid mesh, 
  OutputParams output){
  REAL8 dval1=0., dval2=0., dval3=0., dval4=0.;
 
  REAL8 post1=0., post2=0.;
  
  REAL8 ***evSum1=NULL;      /* first integral */
  REAL8 **evSum2=NULL;       /* second integral */
  REAL8 *evSum3=NULL;        /* third integral */
  REAL8 *cumsum=NULL;        /* cumulative probability */
  REAL8 evSum4=-1.e200;      /* fouth integral */
  REAL8 maxPost=-1.e200;
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
  evSum1 = calloc(numSteps1, sizeof(REAL8 **));
  
  /* perform first integral */
  for( i = 0 ; i < numSteps1 ; i++ ){
    evSum1[i] = calloc(numSteps2, sizeof(REAL8 *));
    for( j = 0 ; j < numSteps2 ; j++ ){
      evSum1[i][j] = calloc(numSteps3, sizeof(REAL8));
      for( k = 0 ; k < numSteps3 ; k++ ){
        evSum1[i][j][k] = -1.e200; /* initialise */
        
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
              sumVal = log_trapezium(post1, post2, dval4);
              evSum1[i][j][k] = PLUS(sumVal, evVal);
            }
            else if( strcmp( output.margParam, "phi" ) == 0 ){
              post1 = logPost[i][j][k][n];
              post2 = logPost[i][j][k][n+1];
  
              evVal = evSum1[i][j][k];
              sumVal = log_trapezium(post1, post2, dval4);
              evSum1[i][j][k] = PLUS(sumVal, evVal);
            }
            else if( strcmp( output.margParam, "psi" ) == 0 ){
              post1 = logPost[j][k][i][n];
              post2 = logPost[j][k][i][n+1];
  
              evVal = evSum1[i][j][k];
              sumVal = log_trapezium(post1, post2, dval4);
              evSum1[i][j][k] = PLUS(sumVal, evVal);
            }
            else if( strcmp( output.margParam, "ciota" ) == 0 ){
              post1 = logPost[j][i][k][n];
              post2 = logPost[j][i][k][n+1];
  
              evVal = evSum1[i][j][k];
              sumVal = log_trapezium(post1, post2, dval4);
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
  evSum2 = calloc(numSteps1, sizeof(REAL8 *));
  
  /* perform the second integration */
  for( i = 0 ; i < numSteps1 ; i++ ){
    evSum2[i] = calloc(numSteps2, sizeof(REAL8));
    for( j = 0 ; j < numSteps2 ; j++ ){
      evSum2[i][j] = -1.e200;
      
      if( numSteps3 == 1 ) evSum2[i][j] = evSum1[i][j][0];
      else{  
        for( k = 0 ; k <numSteps3 - 1 ; k++ ){
          evVal = evSum2[i][j];
          sumVal = log_trapezium(evSum1[i][j][k], evSum1[i][j][k+1], dval3);
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
  evSum3 = calloc(numSteps1, sizeof(REAL8));
  
  /* perform third integration */
  for( i = 0 ; i < numSteps1 ; i++ ){
    evSum3[i] = -1.e200;
    
    if( numSteps2 == 1 ) evSum3[i] = evSum2[i][0];
    else{
      for( j = 0 ; j < numSteps2 - 1 ; j++ ){
        evVal = evSum3[i];
        sumVal = log_trapezium(evSum2[i][j], evSum2[i][j+1], dval2);
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
      sumVal = log_trapezium(evSum3[i], evSum3[i+1], dval1);
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
  cumsum = calloc(numSteps1, sizeof(REAL8));
  
  /* normalise the marginalised pdf using the evidence and free memory */
  for( i = 0 ; i < numSteps1 ; i++ ){
    for( j = 0 ; j < numSteps2 ; j++ ){
      free(evSum1[i][j]);
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
        cumsum[i] = cumsum[i-1] + exp(log_trapezium(evSum3[i-1], evSum3[i],
dval1) - evSum4); 
/*         cumsum[i] = cumsum[i-1] + exp(evSum3[i-1] - evSum4) * dval1; */
      }
    }
    /* print out marginalised posterior */
    fprintf(fp, "%le\t%le\t%le\n", step, exp(evSum3[i] - evSum4), cumsum[i]);
    
    free(evSum1[i]);
    free(evSum2[i]);
  }
  
  results.evidence = evSum4;
  
  /* get the h0 upper limit if required */
  if( strcmp( output.margParam, "h0" ) == 0 && output.dob != 0 )
    results.h0UpperLimit = get_upper_limit(cumsum, output.dob, mesh); 
  
  fclose(fp);

  free(evSum1);
  free(evSum2);
  free(evSum3);
  free(cumsum);

  return results; /* return the log evidence */
}


/* function to do the trapezium rule for integration on logs  */
REAL8 log_trapezium(REAL8 logHeight1, REAL8 logHeight2, REAL8 width){
  REAL8 area=0.;
  
  /* area = 0.5*(height1 + height2)*width 
     logarea = log(0.5) + log(width) + log(exp(logHeight1) + exp(logHeight2)) */
  area = log(0.5) + log(width) + PLUS(logHeight1, logHeight2);
  
  return area;
}


/* detector response lookup table function  - this function will output a lookup
table of points in time and psi, covering a day from the start time (t0) and
from -pi/4 to pi/4 in psi */
void response_lookup_table(REAL8 t0, LALDetAndSource detAndSource,
  DetRespLookupTable *lookupTable){ 
  LIGOTimeGPS gps;
  REAL8 T=0;
  
  REAL8 fplus=0., fcross=0.;

  INT4 i=0, j=0;

  for( i = 0 ; i < lookupTable->psiSteps ; i++ ){
    detAndSource.pSource->orientation = -(LAL_PI/4.) +
        (REAL8)i*(LAL_PI/2.) / ( (REAL8)lookupTable->psiSteps - 1. );

    for( j = 0 ; j < lookupTable->timeSteps ; j++ ){
      /* one day is 86400 seconds */
      T = t0 + (REAL8)j*86400./(REAL8)lookupTable->timeSteps;
    
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


/* function to return the log factorial of an integer */
REAL8 log_factorial(INT4 num){
  INT4 i=0;
  REAL8 logFac=0.;
  
  for( i=2 ; i <= num ; i++ ) logFac += (REAL8)i;
  
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
  if( point2 == 0. ){
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

  /* return smaller value as cumulative prob gradient will always be positive */
  if( ul1 < ul2 )
    return ul1;
  else
    return ul2;
}


/* function to perform the MCMC parameter estimation */
void perform_mcmc(DataStructure *data, InputParams input, INT4 numDets, 
  CHAR *det, LALDetector *detPos, EphemerisData *edat ){
  IntrinsicPulsarVariables vars, varsNew;
  
  /* extra h0 and phase parameter required if a glitch is obsevered */
  IntrinsicPulsarVariables *extraVars=NULL, *extraVarsNew=NULL;
  BinaryPulsarParams pulsarParamsMCMC[2];
  REAL8 *glitchTimes=NULL;
  CHAR *gtimestr=NULL;
  INT4 nGlitches=input.mcmc.nGlitches;
  INT4 **g1=NULL, **g2=NULL; /* start - end positions of each glitch segment */
  DataStructure **glitchData=NULL;
  
  INT4 below0=0;
  
  REAL4Vector *randNum=NULL; /* LAL random variable params */
  UINT4 seed=0;              /* set to get seed from clock time */
  RandomParams *randomParams=NULL;
  
  CHAR *pos1=NULL, *pos2=NULL;
  INT4 i=0, j=0, k=0, n=0, count=0, m=0;
  
  REAL8 like1[1], like2[1];
  REAL8 logL1=0., logL2=0.; /* log likelihoods */
  REAL8 ratio;              /* logL2 - logL1 = log(L2/L1) */
  
  FILE *fp=NULL;
  CHAR outFile[256];
  
  /* variables for pulsar parameters */
  INT4 matTrue=0;
  BinaryPulsarParams pulsarParams, pulsarParamsNew;
  REAL8Array *covmat=NULL, *posdef=NULL, *chol=NULL;
  
  ParamData *paramData=NULL, *randVals=NULL, *vals=NULL;

  BarycenterInput baryinput[2];

  REAL8 prior=0., priorNew=0.;

  /* read the TEMPO par file for the pulsar */
  XLALReadTEMPOParFile( &pulsarParams, input.parFile );

  if( verbose ){
    fprintf(stderr, "Performing an MCMC for %s with %d iterations.\n",
      det, input.mcmc.iterations);
  }  

  /* set up random parameters */
  randomParams = XLALCreateRandomParams( seed );
 
  /* work out how many glitches have been input */
  if( input.mcmc.nGlitches > 0 ){
    extraVars = calloc(nGlitches, sizeof(IntrinsicPulsarVariables));
    extraVarsNew = calloc(nGlitches, sizeof(IntrinsicPulsarVariables));
    glitchTimes = calloc(nGlitches, sizeof(REAL8));
    
    /* start and end point of the data before each glitch */
    g1 = XLALCalloc(numDets, sizeof(INT4 *));
    g2 = XLALCalloc(numDets, sizeof(INT4 *));
    
    for( j = 0 ; j < numDets ; j++ ){
        g1[j] = XLALCalloc(nGlitches + 1, sizeof(INT4));
        g2[j] = XLALCalloc(nGlitches + 1, sizeof(INT4));
    }

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
        gtimestr = malloc(256*sizeof(CHAR));
          
        /* values are delimited by commas , */        
        if( i < nGlitches - 1 ){
          pos2 = strchr(pos1, ',');
          XLALStringCopy(gtimestr, pos1, (pos2-pos1)+1);/*copy time of glitch*/
        }
        else
          gtimestr = XLALStringDuplicate(pos1);
         
        glitchTimes[i] = LALTDBMJDtoGPS(atof(gtimestr)); /* convert to GPS */
        
        XLALFree(gtimestr);
      }
      
      /* set initial values for extra params (all start at same point) */
      if( i==0 ){
        extraVars[i].h0 = input.mesh.minVals.h0 +
          (REAL8)XLALUniformDeviate(randomParams) *
          (input.mesh.maxVals.h0 - input.mesh.minVals.h0);
        extraVars[i].phi0 = input.mesh.minVals.phi0 +
          (REAL8)XLALUniformDeviate(randomParams) *
          (input.mesh.maxVals.phi0 - input.mesh.minVals.phi0);
      }
      else{
        extraVars[i].h0 = extraVars[i-1].h0;
        extraVars[i].phi0 = extraVars[i-1].phi0;
      }

      /* find the position before and after the glitch +/- the cut-off time/2 */
      for( j = 0 ; j < numDets ; j++ ){        
        if( i == 0 )
          g1[j][i] = 0;
         
        for( k = 0 ; k < (INT4)data[j].times->length ; k++ ){
          /* first point after the last glitch */
          if( i > 0 && data[j].times->data[k] < glitchTimes[i-1] + 
              input.mcmc.glitchCut/2.){
              g1[j][i] = k;
           }
          
          /* first point before current glitch */
          if( data[j].times->data[k] > glitchTimes[i] - 
            input.mcmc.glitchCut/2. ){
            g2[j][i] = k;
            
            if(i != nGlitches - 1)
              break;           
          }
          
          /* point from final glitch to end of data */
          if( data[j].times->data[k] > glitchTimes[i] +
             input.mcmc.glitchCut/2. && i == nGlitches - 1 ){
            g1[j][i+1] = k;
            g2[j][i+1] = data[j].times->length;
            
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
        for( k = g1[i][j] ; k < g2[i][j] ; k++ ){
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
  if( input.matrixFile != NULL )
    matTrue = 1;

  if( matTrue ){
    paramData = XLALMalloc(MAXPARAMS*sizeof(ParamData));

    /* get covariance matrix */
    covmat = CreateCovarianceMatrix(input.matrixFile, pulsarParams, paramData);

    /* get covariance matrix, convert it to positive definite if necessary */
    if( (posdef = XLALCheckPositiveDefinite( covmat )) == NULL ){
      /* perform the cholesky decomposition of the matrix */
      chol = CholeskyDecomp(covmat, "lower");
    }
    else
      chol = CholeskyDecomp(posdef, "lower");

    /* generate random parameters */
    randVals = MultivariateNormalDeviates( chol, paramData, randomParams );
  
    vals = XLALMalloc(MAXPARAMS*sizeof(ParamData));
    memcpy(vals, randVals, MAXPARAMS*sizeof(ParamData) );

    /* set up initial pulsar parameters */
    memcpy( &pulsarParamsMCMC[0], &pulsarParams, sizeof(pulsarParams) );
    memcpy( &pulsarParamsMCMC[1], &pulsarParams, sizeof(pulsarParams) );

    SetMCMCPulsarParams( &pulsarParamsMCMC[1], randVals );

    /* check eccentricities so that 0 <= e < 1 i.e. circular of elliptical
       orbits */
    if( pulsarParamsMCMC[1].e < 0. )
      pulsarParamsMCMC[1].e = fabs(pulsarParamsMCMC[1].e);
    else if( pulsarParamsMCMC[1].e >= 1. )
      pulsarParamsMCMC[1].e = 1. - fmod(pulsarParamsMCMC[1].e, 1.);

    if( pulsarParamsMCMC[1].e2 < 0. )
      pulsarParamsMCMC[1].e2 = fabs(pulsarParamsMCMC[1].e2);
    else if( pulsarParamsMCMC[1].e2 >= 1. )
      pulsarParamsMCMC[1].e2 = 1. - fmod(pulsarParamsMCMC[1].e2, 1.);

    if( pulsarParamsMCMC[1].e3 < 0. )
      pulsarParamsMCMC[1].e3 = fabs(pulsarParamsMCMC[1].e3);
    else if( pulsarParamsMCMC[1].e3 >= 1. )
      pulsarParamsMCMC[1].e3 = 1. - fmod(pulsarParamsMCMC[1].e3, 1.);
  }
  else{
    /* make sure edat is NULL as this is what is used by the log_likelihood
       function to determine if a covariance matrix has been used or not */
    edat = NULL;
  }

  /* set initial chain parameters */
  vars.h0 = input.mesh.minVals.h0 + (REAL8)XLALUniformDeviate(randomParams) *
             (input.mesh.maxVals.h0 - input.mesh.minVals.h0);
  if( input.mcmc.nGlitches > 0 )
    vars.h0 = extraVars[0].h0;
  vars.phi0 = input.mesh.minVals.phi0 + (REAL8)XLALUniformDeviate(randomParams)
             * (input.mesh.maxVals.phi0 - input.mesh.minVals.phi0);
  vars.psi = input.mesh.minVals.psi + (REAL8)XLALUniformDeviate(randomParams) *
              (input.mesh.maxVals.psi - input.mesh.minVals.psi);
  vars.ci = input.mesh.minVals.ci + (REAL8)XLALUniformDeviate(randomParams) *
            (input.mesh.maxVals.ci - input.mesh.minVals.ci);

  vars.Xplus = 0.5*(1.+vars.ci*vars.ci);
  vars.Xcross = vars.ci;
  vars.Xpsinphi_2 = 0.5*vars.Xplus*sin(vars.phi0);
  vars.Xcsinphi_2 = 0.5*vars.Xcross*sin(vars.phi0);
  vars.Xpcosphi_2 = 0.5*vars.Xplus*cos(vars.phi0);
  vars.Xccosphi_2 = 0.5*vars.Xcross*cos(vars.phi0);
  
  for( i = 0 ; i < nGlitches ; i++ ){
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
  sprintf(outFile, "%s/MCMCchain_%s_%s", input.outputDir, input.pulsar, det);
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
    fprintf(fp, "\th0(%d)       \tphi0(%d) ", j+2, j+2);
  if( matTrue ){
    for( j = 0 ; j < MAXPARAMS ; j++ ){
      if( paramData[j].matPos != 0 )
        fprintf(fp, "\t%s", paramData[j].name);
    }
  }
  fprintf(fp, "\n");

  /* create vector for random Gaussian numbers */
  randNum = XLALCreateREAL4Vector(4+2*nGlitches);
  
  count = 0;

  if( matTrue && numDets == 1 ){
    /* set up detector location - if not doing joint analysis */
    baryinput[0].site = *detPos;
    baryinput[0].site.location[0] /= LAL_C_SI;
    baryinput[0].site.location[1] /= LAL_C_SI;
    baryinput[0].site.location[2] /= LAL_C_SI;

    baryinput[1].site = *detPos;
    baryinput[1].site.location[0] /= LAL_C_SI;
    baryinput[1].site.location[1] /= LAL_C_SI;
    baryinput[1].site.location[2] /= LAL_C_SI;
  }

  /*====================== MCMC LOOP =========================================*/
  for( i = 0 ; i < input.mcmc.iterations + input.mcmc.burnIn ; i++ ){
    if( verbose ){ 
      fprintf(stderr, "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
      fprintf(stderr, "%06.2lf%% complete", ((double)(i)+1.)*100. /
        (double)(input.mcmc.iterations + input.mcmc.burnIn));
    }

    /* get new values of parameters using Gaussian proposal distributions */
    XLALNormalDeviates(randNum, randomParams);
    varsNew.h0 = vars.h0 + input.mcmc.sigmas.h0*randNum->data[0];
    
    below0 = 0;
    for( j = 0 ; j < nGlitches ; j++ ){
      extraVarsNew[j].h0 = extraVars[j].h0 +
        input.mcmc.sigmas.h0*randNum->data[4+2*j];
      extraVarsNew[j].phi0 = extraVars[j].phi0 +
        input.mcmc.sigmas.phi0*randNum->data[4+2*j+1];
      
      if( extraVarsNew[j].phi0 < 0. ){
        extraVarsNew[j].phi0 = fmod(extraVarsNew[j].phi0, LAL_TWOPI);
        extraVarsNew[j].phi0 += LAL_TWOPI;
      }
      else if( extraVarsNew[j].phi0 > LAL_TWOPI )
        extraVarsNew[j].phi0 = fmod(extraVarsNew[j].phi0, LAL_TWOPI);
      
      /* if any of the h0 jump below zero then write out data below */
      if( extraVarsNew[j].h0 < 0. )
        below0 = 1;

      /* make it so that values of h0 after a glitch must be smaller than
         before the glitch i.e. the pulsar always relaxes to a lower energy
         state - if this happens set below0 = 1 to output the current state */
      if( j==0 ){
        if( extraVarsNew[j].h0 > varsNew.h0 )
          below0 = 1;
      }
      else{
        if( extraVarsNew[j].h0 > extraVarsNew[j-1].h0 )
          below0 = 1;
      }
    }
    
    /* if h0 jumps negative then this is equivalent to having jumped outside
       our prior range so the likelihood is always zero and this move always
       rejected */
    if( varsNew.h0 < 0. || below0 == 1 ){
      if( fmod(count, input.mcmc.outputRate) == 0. && i > input.mcmc.burnIn-1 ){
        fprintf(fp, "%le\t%le\t%lf\t%lf\t%lf", logL1, vars.h0, vars.phi0,
          vars.ci, vars.psi);
         
        for( j = 0 ; j < nGlitches ; j++ )
          fprintf(fp, "\t%le\t%lf", extraVars[j].h0, extraVars[j].phi0);
        
        if( matTrue ){
          /* print out pulsar parameters */
          for( j = 0 ; j < MAXPARAMS ; j++ ){
            if( paramData[j].matPos != 0 )
              fprintf(fp, "\t%.17le", vals[j].val-paramData[j].val);
          }
        }
        fprintf(fp, "\n");
      }
      
      count++;
      continue;
    }

    if( matTrue ){
      /* generate new pulsars parameters from the pulsar parameter covariance */
      randVals = MultivariateNormalDeviates( chol, paramData, randomParams );
      memcpy(&pulsarParamsNew, &pulsarParams, sizeof(BinaryPulsarParams));
      SetMCMCPulsarParams( &pulsarParamsNew, randVals );

      /* check eccentricities so that 0 <= e < 1 i.e. circular of elliptical
       orbits */
      if( pulsarParamsNew.e < 0. )
        pulsarParamsNew.e = fabs(pulsarParamsNew.e);
      else if( pulsarParamsNew.e >= 1. )
        pulsarParamsNew.e = 1. - fmod(pulsarParamsNew.e, 1.);

      if( pulsarParamsNew.e2 < 0. )
        pulsarParamsNew.e2 = fabs(pulsarParamsNew.e2);
      else if( pulsarParamsNew.e2 >= 1. )
        pulsarParamsNew.e2 = 1. - fmod(pulsarParamsNew.e2, 1.);

      if( pulsarParamsNew.e3 < 0. )
        pulsarParamsNew.e3 = fabs(pulsarParamsNew.e3);
      else if( pulsarParamsNew.e3 >= 1. )
        pulsarParamsNew.e3 = 1. - fmod(pulsarParamsNew.e3, 1.);
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
    varsNew.Xpsinphi_2 = 0.5*varsNew.Xplus*sin(varsNew.phi0);
    varsNew.Xcsinphi_2 = 0.5*varsNew.Xcross*sin(varsNew.phi0);
    varsNew.Xpcosphi_2 = 0.5*varsNew.Xplus*cos(varsNew.phi0);
    varsNew.Xccosphi_2 = 0.5*varsNew.Xcross*cos(varsNew.phi0);
    
    for( j = 0 ; j < nGlitches ; j++ ){
      extraVarsNew[j].Xpsinphi_2 = 0.5*varsNew.Xplus *
                                   sin(extraVarsNew[j].phi0);
      extraVarsNew[j].Xcsinphi_2 = 0.5*varsNew.Xcross *
                                   sin(extraVarsNew[j].phi0);
      extraVarsNew[j].Xpcosphi_2 = 0.5*varsNew.Xplus *
                                   cos(extraVarsNew[j].phi0);
      extraVarsNew[j].Xccosphi_2 = 0.5*varsNew.Xcross *
                                   cos(extraVarsNew[j].phi0);
      extraVarsNew[j].psi = varsNew.psi;
    }
    
    /* set single h0 value for log_likelihood function */
    input.mesh.h0Steps = 1;
    
    /* calculate log likelihoods */
    /* only calculate the likelhood twice in the first loop */
    if( i == 0 ){
      for( k = 0 ; k < numDets ; k++ ){
        if( matTrue && numDets > 1 ){
          /* set up detector location for joint likelihood */
          baryinput[0].site = detPos[k];
          baryinput[0].site.location[0] /= LAL_C_SI;
          baryinput[0].site.location[1] /= LAL_C_SI;
          baryinput[0].site.location[2] /= LAL_C_SI;

          baryinput[1].site = detPos[k];
          baryinput[1].site.location[0] /= LAL_C_SI;
          baryinput[1].site.location[1] /= LAL_C_SI;
          baryinput[1].site.location[2] /= LAL_C_SI;
        }

        /* first likelihood */
        if( nGlitches == 0 ){
          input.mesh.minVals.h0 = vars.h0;
          log_likelihood(like1, data[k], vars, input.mesh, pulsarParamsMCMC,
            baryinput, edat);
          logL1 += *like1;
        }
        else{
          for( j = 0 ; j < nGlitches + 1 ; j++ ){
            if( j == 0 ){
              input.mesh.minVals.h0 = vars.h0;
              
              log_likelihood(like1, glitchData[k][j], vars, input.mesh,
                pulsarParamsMCMC, baryinput, edat);
            }
            else{
              input.mesh.minVals.h0 = extraVars[j-1].h0;
              log_likelihood(like1, glitchData[k][j], extraVars[j-1],
                input.mesh, pulsarParamsMCMC, baryinput, edat);
            }
            logL1 += *like1;
          }
        }
      }
    }
      
    if( matTrue ){
      /* set new pulsar parameters */
      /* store originals */
      memcpy(&pulsarParams, &pulsarParamsMCMC[1], sizeof(BinaryPulsarParams));
      memcpy(&pulsarParamsMCMC[1], &pulsarParamsNew,
        sizeof(BinaryPulsarParams));
    }  

    for( k = 0 ; k < numDets ; k++ ){
      /* second likelihood - for new position in parameter space */
      if( matTrue && numDets > 1 ){
        /* set up detector location for joint likelihood */
        baryinput[0].site = detPos[k];
        baryinput[0].site.location[0] /= LAL_C_SI;
        baryinput[0].site.location[1] /= LAL_C_SI;
        baryinput[0].site.location[2] /= LAL_C_SI;

        baryinput[1].site = detPos[k];
        baryinput[1].site.location[0] /= LAL_C_SI;
        baryinput[1].site.location[1] /= LAL_C_SI;
        baryinput[1].site.location[2] /= LAL_C_SI;
      }

      if( nGlitches == 0 ){
        input.mesh.minVals.h0 = varsNew.h0;
        log_likelihood(like2, data[k], varsNew, input.mesh, pulsarParamsMCMC,
          baryinput, edat);
        logL2 += *like2;
      }
      else{
        for( j = 0 ; j < nGlitches + 1 ; j++ ){
          if( j == 0 ){
            input.mesh.minVals.h0 = varsNew.h0;
            log_likelihood(like2, glitchData[k][j], varsNew, input.mesh,
              pulsarParamsMCMC, baryinput, edat);
          }
          else{
            input.mesh.minVals.h0 = extraVarsNew[j-1].h0;
            log_likelihood(like2, glitchData[k][j], extraVarsNew[j-1],
              input.mesh, pulsarParamsMCMC, baryinput, edat);
          }  
          logL2 += *like2;          
        }
      }
    }
    
    /* most importantly use the covariance matrix as a prior on the pulsar
       parameters */
    prior = 0.;
    priorNew = 0.;
    if( matTrue ){
      n=0;
      m=0;
      for( j=0; j<MAXPARAMS ; j++){
        if( paramData[j].matPos != 0 ){
          for( k=0; k<MAXPARAMS ; k++){
            if( paramData[k].matPos != 0 ){
              /* log prior */
              if( i == 0 ){ /* first likelihood */
                prior -= (vals[k].val - paramData[k].val) *
                  XLALGetREAL8MatrixValue(covmat, n, m) * (vals[j].val - 
                  paramData[j].val);
              }
          
              priorNew -= (randVals[k].val - paramData[k].val) *
                XLALGetREAL8MatrixValue(covmat, n, m) * (randVals[j].val -
                paramData[j].val);

              m++;
            }
          }
          n++;
        }
      }
    }

    logL1 += prior;
    logL2 += priorNew;

    /* set and apply priors on h0, psi, cos(iota) and phi is requested */
    if( input.usepriors ){
      input.priors.vars.h0 = vars.h0;
      input.priors.vars.phi0 = vars.phi0;
      input.priors.vars.ci = vars.ci;
      input.priors.vars.psi = vars.psi;

      logL1 += log_prior(input.priors, input.mesh);
      
      if( nGlitches > 0 ){
        for( j = 0 ; j < nGlitches ; j++ ){
          input.priors.vars.h0 = extraVars[j].h0;

          logL1 += log_prior(input.priors, input.mesh);
        }
      }

      input.priors.vars.h0 = varsNew.h0;
      input.priors.vars.phi0 = varsNew.phi0;
      input.priors.vars.ci = varsNew.ci;
      input.priors.vars.psi = varsNew.psi;

      logL2 += log_prior(input.priors, input.mesh);

      if( nGlitches > 0 ){
        for( j = 0 ; j < nGlitches ; j++ ){
          input.priors.vars.h0 = extraVarsNew[j].h0;

          logL2 += log_prior(input.priors, input.mesh);
        }
      }
    }

    /* accept new values if Lnew/Lold >=1 (or log(Lnew/Lold) >= 0) */
    /* include simulated annealing factor */
    if( i < input.mcmc.burnIn-1 ){
      ratio = input.mcmc.temperature * exp( log(1./input.mcmc.temperature) *
              (double)i/(double)input.mcmc.burnIn) * (logL2 - logL1);
    }
    else
      ratio = logL2 - logL1;
      
    if( ratio >= 0. ){ /* always accept new value */
      vars = varsNew;
      /* update vals structure */ 
      if( matTrue )
        memcpy(vals, randVals, MAXPARAMS*sizeof(ParamData));

      for( j = 0 ; j < nGlitches ; j++)
        extraVars[j] = extraVarsNew[j];
        
      logL1 = logL2;     
    }
    /* otherwise accept with a certain probability */
    else if( log(XLALUniformDeviate(randomParams)) < ratio ){
      vars = varsNew;
      
      if( matTrue )
        memcpy(vals, randVals, MAXPARAMS*sizeof(ParamData));   

      for( j = 0 ; j < nGlitches ; j++)
        extraVars[j] = extraVarsNew[j];
        
      logL1 = logL2;
    }

    /* printf out chains */
    if( fmod(count, input.mcmc.outputRate) == 0. && i > input.mcmc.burnIn-1 ){
      fprintf(fp, "%le\t%le\t%lf\t%lf\t%lf", logL1, vars.h0, vars.phi0, vars.ci,
        vars.psi);
         
      for( j = 0 ; j < nGlitches ; j++ )
        fprintf(fp, "\t%le\t%lf", extraVars[j].h0, extraVars[j].phi0);
      
      if( matTrue ){
        /* print out pulsar parameters */
        for( j = 0 ; j < MAXPARAMS ; j++ ){
          /* output difference in parameter from heterodyne parameter */
          if( paramData[j].matPos != 0 )
            fprintf(fp, "\t%.17le", vals[j].val-paramData[j].val);
        }
      }
      
      fprintf(fp, "\n");
    }
    
    count++;
    logL2 = 0.;
  }
  /*==========================================================================*/

  if( verbose ) fprintf(stderr, "\n");
  
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

    XLALDestroyREAL8Array( covmat );
    XLALDestroyREAL8Array( chol );

    if( posdef != NULL )
      XLALDestroyREAL8Array( posdef );
  }

  fclose(fp);
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

void SetMCMCPulsarParams( BinaryPulsarParams *pulsarParams, ParamData *data ){
  INT4 i=0;

  /* go through params and set the next step in the MCMC */
  for( i=0 ; i<MAXPARAMS; i++ ){
    if( !strcmp(data[i].name, "f0") && data[i].matPos != 0 )
      pulsarParams->f0 = data[i].val;
    else if( !strcmp(data[i].name, "f1") && data[i].matPos != 0 )
      pulsarParams->f1 = data[i].val;
    else if( !strcmp(data[i].name, "f2") && data[i].matPos != 0 )
      pulsarParams->f2 = data[i].val;
    else if( !strcmp(data[i].name, "Dec") && data[i].matPos != 0 )
      pulsarParams->dec = data[i].val;
    else if( !strcmp(data[i].name, "RA") && data[i].matPos != 0 )
      pulsarParams->ra = data[i].val;
    else if( !strcmp(data[i].name, "pmdc") && data[i].matPos != 0 )
      pulsarParams->pmdec = data[i].val;
    else if( !strcmp(data[i].name, "pmra") && data[i].matPos != 0 )
      pulsarParams->pmra = data[i].val;
    else if( !strcmp(data[i].name, "x") && data[i].matPos != 0 )
      pulsarParams->x = data[i].val;
    else if( !strcmp(data[i].name, "e") && data[i].matPos != 0 )
      pulsarParams->e = data[i].val;
    else if( !strcmp(data[i].name, "T0") && data[i].matPos != 0 )
      pulsarParams->T0 = data[i].val;
    else if( !strcmp(data[i].name, "Pb") && data[i].matPos != 0 )
      pulsarParams->Pb = data[i].val;
    else if( !strcmp(data[i].name, "Om") && data[i].matPos != 0 )
      pulsarParams->w0 = data[i].val;
    else if( !strcmp(data[i].name, "Omdt") && data[i].matPos != 0 )
      pulsarParams->wdot = data[i].val;
    else if( !strcmp(data[i].name, "gamma") && data[i].matPos != 0 )
      pulsarParams->gamma = data[i].val;
    else if( !strcmp(data[i].name, "Pbdt") && data[i].matPos != 0 )
      pulsarParams->Pbdot = data[i].val;
    else if( !strcmp(data[i].name, "s") && data[i].matPos != 0 )
      pulsarParams->s = data[i].val;
    else if( !strcmp(data[i].name, "M") && data[i].matPos != 0 )
      pulsarParams->M = data[i].val;
    else if( !strcmp(data[i].name, "m2") && data[i].matPos != 0 )
      pulsarParams->m2 = data[i].val;
    else if( !strcmp(data[i].name, "dth") && data[i].matPos != 0 )
      pulsarParams->dth = data[i].val;
    else if( !strcmp(data[i].name, "xdot") && data[i].matPos != 0 )
      pulsarParams->xdot = data[i].val;
    else if( !strcmp(data[i].name, "edot") && data[i].matPos != 0 )
      pulsarParams->edot = data[i].val;
    else if( !strcmp(data[i].name, "x2") && data[i].matPos != 0 )
      pulsarParams->x2 = data[i].val;
    else if( !strcmp(data[i].name, "e2") && data[i].matPos != 0 )
      pulsarParams->e2 = data[i].val;
    else if( !strcmp(data[i].name, "T02") && data[i].matPos != 0 )
      pulsarParams->T02 = data[i].val;
    else if( !strcmp(data[i].name, "Pb2") && data[i].matPos != 0 )
      pulsarParams->Pb2 = data[i].val;
    else if( !strcmp(data[i].name, "Om2") && data[i].matPos != 0 )
      pulsarParams->w02 = data[i].val;
    else if( !strcmp(data[i].name, "x3") && data[i].matPos != 0 )
      pulsarParams->x3 = data[i].val;
    else if( !strcmp(data[i].name, "e3") && data[i].matPos != 0 )
      pulsarParams->e3 = data[i].val;
    else if( !strcmp(data[i].name, "T03") && data[i].matPos != 0 )
      pulsarParams->T03 = data[i].val;
    else if( !strcmp(data[i].name, "Pb3") && data[i].matPos != 0 )
      pulsarParams->Pb3 = data[i].val;
    else if( !strcmp(data[i].name, "Om3") && data[i].matPos != 0 )
      pulsarParams->w03 = data[i].val;
    else if( !strcmp(data[i].name, "Xpbd") && data[i].matPos != 0 )
      pulsarParams->xpbdot = data[i].val;
    else if( !strcmp(data[i].name, "eps1") && data[i].matPos != 0 )
      pulsarParams->eps1 = data[i].val;
    else if( !strcmp(data[i].name, "eps2") && data[i].matPos != 0 )
      pulsarParams->eps2 = data[i].val;
    else if( !strcmp(data[i].name, "e1dt") && data[i].matPos != 0 )
      pulsarParams->eps1dot = data[i].val;
    else if( !strcmp(data[i].name, "e2dt") && data[i].matPos != 0 )
      pulsarParams->eps2dot = data[i].val;
  }
}

/* this function perform Cholesky decomposition on M and outputs the
   lower, or upper triagular matrix depending if uOrl is set to "upper" or
   "lower" - if nothing is specified then the default is lower 
   This is pretty much copied from the GSL function gsl_linalg_cholesky_decomp
   although this works with floats rather than doubles
*/ 
REAL8Array *CholeskyDecomp(REAL8Array *M, CHAR* uOrl){
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
      XLALSetREAL8MatrixValue( A, i, j, XLALGetREAL8MatrixValue(M, i, j) );

  A_00 = XLALGetREAL8MatrixValue( A, 0, 0 );
  L_00 = sqrt(A_00);

  if( A_00 <= 0 )
    fprintf(stderr, "Error... matrix must be positive definite!\n");

  XLALSetREAL8MatrixValue( A, 0, 0, L_00 );

  if( length > 1 ){
    REAL8 A_10 = XLALGetREAL8MatrixValue( A, 1, 0 );
    REAL8 A_11 = XLALGetREAL8MatrixValue( A, 1, 1 );

    REAL8 L_10 = A_10/L_00;
    REAL8 diag = A_11 - L_10*L_10;
    REAL8 L_11 = sqrt(diag);

    if( diag <= 0 ){
      fprintf(stderr, "Error... input matrix is not pos def!\n");
      exit(0);
    }

    XLALSetREAL8MatrixValue( A, 1, 0, L_10 );
    XLALSetREAL8MatrixValue( A, 1, 1, L_11 );
  }

  for( k=2; k<length; k++ ){
    REAL8 A_kk = XLALGetREAL8MatrixValue( A, k, k );

    for( i=0; i<k; i++ ){
      REAL8 sum = 0.;
      
      REAL8 A_ki = XLALGetREAL8MatrixValue( A, k, i );
      REAL8 A_ii = XLALGetREAL8MatrixValue( A, i, i );
    
      REAL8 ci[length];
      REAL8 ck[length];

      for( j=0; j<length; j++ ){
        ci[j] = XLALGetREAL8MatrixValue( A, i, j );
        ck[j] = XLALGetREAL8MatrixValue( A, k, j );
      }

      if( i>0 ){
        for( j=0; j<i; j++ )
          sum += ci[j]*ck[j];
      }

      A_ki = (A_ki - sum) / A_ii;
      XLALSetREAL8MatrixValue( A, k, i, A_ki );
    }

    {
      REAL8 sum = 0.;
      REAL8 diag = 0.;
      
      for( j=0; j<k; j++ ){
        sum += XLALGetREAL8MatrixValue( A, k, j )*
          XLALGetREAL8MatrixValue( A, k, j );
      }

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
      else if( diag <= 0. && fabs(diag) <= LAL_REAL8_EPS ){
        diag = LAL_REAL8_EPS;
      }
      else if( diag <= 0. && fabs(diag) >= LAL_REAL8_EPS && k == length-1 ){
        /* this is a kludge as a lot of the matricies seem to have entries
           there m(length, length) diagonal value as small but less than zero,
           so I'll just set it to zero manually */
        diag = 0.;
      }
      
      XLALSetREAL8MatrixValue( A, k, k, sqrt(diag) );
      
    }
  }

  /* set upper triangular matrix to zeros - for lower value */
  for(i=0; i<length; i++)
    for(j=i+1; j<length; j++)
      XLALSetREAL8MatrixValue( A, i, j, 0. );

  /* check if the upper triangle is wanted - if so perform transpose */
  if(strstr(uOrl, "upper")!=NULL){
    REAL8 tempdata = 0.;

    /* perform transpose */
    for(j=0; j<length-1; j++){
      for(i=j; i<length; i++){
        tempdata = XLALGetREAL8MatrixValue( A, j, i );
        XLALSetREAL8MatrixValue( A, j, i, XLALGetREAL8MatrixValue(A, i, j) );
        XLALSetREAL8MatrixValue( A, i, j, tempdata );
      }
    }
  }

  /* print out matrix */
  /* if( verbose ){
    fprintf(stderr, "\nCholesky decomposed matrix:\n");
    for(i=0; i<length; i++){
      for(j=0; j<length; j++)
        fprintf(stderr, "%.2le  ", XLALGetREAL8MatrixValue( A, i, j ));
      fprintf(stderr, "\n");
    }
  } */

  return A;
}

/* this function will draw a set of random numbers from a multivariate Gaussian
distribution, with a cholesky decomposed covariance matrix given by cholmat and
parameter mean
values */
ParamData *MultivariateNormalDeviates( REAL8Array *cholmat, ParamData *data,
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

  /* multiply L by randNum */
  Z = XLALCreateREAL8Vector( dim );
  for(i=0;i<dim;i++)
    Z->data[i] = 0.;

  for(i=0;i<dim;i++){
    for(j=0;j<dim;j++)
      Z->data[i] += XLALGetREAL8MatrixValue(cholmat,i,j)*randNum->data[j];
    /* fprintf(stderr, "Z->data[%d] = %le\n", i, Z->data[i]); */
  }

  /* get the output random deviates by doing the mean plus Z */
  j=0;
  for(i=0;i<MAXPARAMS;i++){
    deviates[i].name = data[i].name;
    deviates[i].sigma = data[i].sigma;
    deviates[i].matPos = data[i].matPos;
    if( data[i].matPos != 0 ){
      deviates[i].val = data[i].val + Z->data[j];
      
      j++;
    }
    else
      deviates[i].val = data[i].val;
  }

  XLALDestroyREAL4Vector( randNum );
  XLALDestroyREAL8Vector( Z );

  return deviates;
}

/* I need to define a standard set of positions in which various pulsar
   parameters will sit within the internal correlation matrix - this will as
   far as possible following the standard in the matrix files I have from
   Michael Kramer
*/
REAL8Array *CreateCovarianceMatrix( CHAR *matrixFile, 
BinaryPulsarParams params, ParamData *data ){
  FILE *fp=NULL;

  CHAR matrixParams[MAXPARAMS][6]; /* parameters in the correlation matrix */
  CHAR paramTmp[256];

  INT4 numParams=0, i=0, j=0, k=0, n=0;

  CHAR tmpStr[256], tmpStr2[256];

  INT4 arraySize=0;

  REAL8Array *corMat=NULL, *covMat=NULL;
  UINT4Vector *matdims=NULL;

  INT4 DMpos=0, DM1pos=0; /* position of dispersion measure in matrix */
  INT4 numDM=0;
  REAL8 corTemp=0., junk=0.;

  ParamData paramData[]=
  {
    { "f0",    params.f0,     params.f0Err,      0 },
    { "f1",    params.f1,     params.f1Err,      0 },
    { "f2",    params.f2,     params.f2Err,      0 },
    { "Dec",   params.dec,    params.decErr,     0 },
    { "RA",    params.ra,     params.decErr,     0 },
    { "pmdc",  params.pmdec,  params.pmdecErr,   0 },
    { "pmra",  params.pmra,   params.pmraErr,    0 },
    { "x",     params.x,      params.xErr,       0 },
    { "e",     params.e,      params.eErr,       0 },
    { "T0",    params.T0,     params.T0Err,      0 },
    { "Pb",    params.Pb,     params.PbErr,      0 },
    { "Om",    params.w0,     params.w0Err,      0 },
    { "Omdt",  params.wdot,   params.wdotErr,    0 },
    { "gamma", params.gamma,  params.gammaErr,   0 },
    { "Pbdt",  params.Pbdot,  params.PbdotErr,   0 },
    { "s",     params.s,      params.sErr,       0 },
    { "M",     params.M,      params.MErr,       0 },
    { "m2",    params.m2,     params.m2Err,      0 },
    { "dth",   params.dth,    params.dthErr,     0 },
    { "xdot",  params.xdot,   params.xdotErr,    0 },
    { "edot",  params.edot,   params.edotErr,    0 },
    { "x2",    params.x2,     params.x2Err,      0 },
    { "e2",    params.e2,     params.e2Err,      0 },
    { "T02",   params.T02,    params.T02Err,     0 },
    { "Pb2",   params.Pb2,    params.Pb2Err,     0 },
    { "Om2",   params.w02,    params.w02Err,     0 },
    { "x3",    params.x3,     params.x3Err,      0 },
    { "e3",    params.e3,     params.e3Err,      0 },
    { "T03",   params.T03,    params.T03Err,     0 },
    { "Pb3",   params.Pb3,    params.Pb3Err,     0 },
    { "Om3",   params.w03,    params.w03Err,     0 },
    { "Xpbd",  params.xpbdot, params.xpbdotErr,  0 }
  };

  arraySize = MAXPARAMS;

  if( params.model != NULL && !strcmp(params.model, "ELL1") ){
    paramData[8].name = "eps1";
    paramData[8].val = params.eps1;
    paramData[8].sigma = params.eps1Err;

    paramData[11].name = "eps2";
    paramData[11].val = params.eps2;
    paramData[11].sigma = params.eps2Err;

    paramData[10].val = params.Tasc;
    paramData[10].sigma = params.TascErr;

    paramData[12].name = "e1dt";
    paramData[12].val = params.eps1dot;
    paramData[12].sigma = params.eps1dotErr;

    paramData[20].name = "e2dt";
    paramData[20].val = params.eps2dot;
    paramData[20].sigma = params.eps2dotErr;
  }

  /* read in data from correlation matrix file */
  if((fp = fopen(matrixFile, "r")) == NULL){
    fprintf(stderr, "No correlation matrix file" );
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

    /* check if parameter is actually for a dispersion measure (ignore if so) */
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

  covMat = XLALCreateREAL8Array( matdims );
  corMat = XLALCreateREAL8Array( matdims );

  /* find positions of each parameter */
  /* the strings that represent parameters in a matrix are given in the paraem
     variable in the tempo code mxprt.f */  
  /* read in matrix */
  k=0;
  for(i=0;i<numParams+numDM;i++){
    n=0; 
    fscanf(fp, "%s%s", tmpStr, tmpStr2);
    
    /* if its a dispersion measure then just skip the line */
    if( (DMpos != 0 && i == DMpos) || (DM1pos != 0 && i == DM1pos) ){
      fscanf(fp, "%*[^\n]");
      k--;
      continue;
    }

    for(j=0;j<i+1;j++){
      if( (DMpos != 0 && j == DMpos) || (DM1pos != 0 && j == DM1pos) ){
        fscanf(fp, "%lf", &junk);
        n--;
        continue;
      }
      
      fscanf(fp, "%lf", &corTemp);
      XLALSetREAL8MatrixValue( corMat, k, n, corTemp );
      if(n != k)
        XLALSetREAL8MatrixValue( corMat, n, k, corTemp );
      n++;
    }

    /* send an error if we hit the end of the file */
    if(feof(fp)){
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

  /* pass the parameter data to be output */
  memcpy(data, paramData, sizeof(paramData));

  /* if( verbose ){
    fprintf(stderr, "\nCorrelation matrix:\n");
    for(i=0;i<numParams;i++){
      for(j=0;j<numParams;j++)
        fprintf(stderr, "%.2lf  ", XLALGetREAL8MatrixValue(corMat,i,j));
      fprintf(stderr, "\n");
    }
  } */

  /* convert correlation matrix into a covariance matrix */
  for(i=0;i<arraySize;i++){
    if( paramData[i].matPos != 0 ){
      for(j=0;j<arraySize;j++){
        if( paramData[j].matPos != 0 ){
          /* get (co)variances */
          XLALSetREAL8MatrixValue( covMat, paramData[i].matPos-1,
paramData[j].matPos-1, XLALGetREAL8MatrixValue( corMat, paramData[i].matPos-1,
paramData[j].matPos-1) * paramData[i].sigma * paramData[j].sigma );
        }
      }
    }
  }

  /* if( verbose ){
    fprintf(stderr, "\nCovariance matrix:\n");
    for(i=0;i<numParams;i++){
      for(j=0;j<numParams;j++)
        fprintf(stderr, "%.2e  ", XLALGetREAL8MatrixValue(covMat,i,j));
      fprintf(stderr, "\n");
    }
  } */

  XLALDestroyUINT4Vector( matdims );
  XLALDestroyREAL8Array( corMat );

  return covMat;
}

REAL8Array *XLALCheckPositiveDefinite( REAL8Array *matrix ){
  static LALStatus status;

  REAL8Vector *eigenval=NULL;
  REAL8Array *eigenvec=NULL;

  REAL8Array *posdef=NULL;

  INT4 i=0, j=0;

  /* copy input array into eigenvec as this gets converted by function */
  eigenvec = XLALCreateREAL8Array( matrix->dimLength );

  for( i=0; i<(INT4)eigenvec->dimLength->data[0]; i++ ){
    for( j=0; j<(INT4)eigenvec->dimLength->data[1]; j++ ){
      XLALSetREAL8MatrixValue( eigenvec, i, j, 
        XLALGetREAL8MatrixValue(matrix, i, j) );
    }
  }

  eigenval = XLALCreateREAL8Vector( matrix->dimLength->data[0] );

  /* calculate the eigen values and vectors */
  LAL_CALL( LALDSymmetricEigenVectors( &status, eigenval, eigenvec ), &status );

  for( i=0; i<(INT4)matrix->dimLength->data[0]; i++ ){
    /* first check if any eigen values are zero and if so convert to positive
       definite matrix */
    if( eigenval->data[i] < 0. ){
      fprintf(stderr, "Eigenvalue is negative. Non-postive definite matrix!\n");
      posdef = XLALConvertToPositiveDefinite( matrix );
      break;
    }
  }

  /* if matrix is positive definite return it i.e. posdef hasn't been set */
  if( posdef == NULL ){
    XLALDestroyREAL8Array( eigenvec );
    XLALDestroyREAL8Vector( eigenval );
    return NULL;
  }

  /* re-check new matrix for positive definiteness, but be aware of values
     close to the precision of REAL8 numbers */
  for( i=0; i<(INT4)eigenvec->dimLength->data[0]; i++ ){
    for( j=0; j<(INT4)eigenvec->dimLength->data[1]; j++ ){
      XLALSetREAL8MatrixValue( eigenvec, i, j, 
        XLALGetREAL8MatrixValue( posdef, i, j) );
    }
    eigenval->data[i] = 0.;
  }
  
  LAL_CALL( LALDSymmetricEigenVectors( &status, eigenval, eigenvec ), &status );

  for( i=0; i<(INT4)matrix->dimLength->data[0]; i++ ){
    if( eigenval->data[i] < 0. && fabs(eigenval->data[i]) > LAL_REAL8_EPS){
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
into a postive definite matrix using the method (number 2) of Rebonato and
Jackel (see their paper at
http://www.riccardorebonato.co.uk/papers/ValCorMat.pdf */
REAL8Array *XLALConvertToPositiveDefinite( REAL8Array *nonposdef ){
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
  if( nonposdef->dimLength->data[0] != nonposdef->dimLength->data[0] ){
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
      XLALSetREAL8MatrixValue( eigenvec, i, j, 
        XLALGetREAL8MatrixValue(nonposdef, i, j) );

      /* initialise Lprime and T to zeros */
      XLALSetREAL8MatrixValue( Lprime, i, j, 0. );
      XLALSetREAL8MatrixValue( T, i, j, 0. );
    }
  }

  eigenval = XLALCreateREAL8Vector( length );

  /* calculate the eigen values and vectors */
  LAL_CALL( LALDSymmetricEigenVectors( &status, eigenval, eigenvec ), &status );
  /* LALDSymmetricEigenVectors( &status, eigenval, eigenvec ); */

  /* if eigen value is > 0 set Lprime to that value i.e. have eigen values of 
     zero if eigen value is negative */
  for( i=0; i<length; i++ )
    if( eigenval->data[i] > 0. )
      XLALSetREAL8MatrixValue( Lprime, i, i, eigenval->data[i] );

  /* compute scaling matrix T */
  for( i=0; i<length; i++ ){
    Tval = 0.;
    for( j=0; j<length; j++ ){
      Tval += XLALGetREAL8MatrixValue( eigenvec, i, j ) * 
              XLALGetREAL8MatrixValue( eigenvec, i, j ) *
              XLALGetREAL8MatrixValue( Lprime, j, j );
    }

    Tval = 1./Tval;

    /* really we just want the sqrt of T */
    XLALSetREAL8MatrixValue( T, i, i, sqrt(Tval) );
  }

  /* convert Lprime to sqrt(lambdaprime) */
  for( i=0; i<length; i++ ){
    REAL8 tempL = XLALGetREAL8MatrixValue(Lprime, i, i);

    XLALSetREAL8MatrixValue( Lprime, i, i, sqrt(tempL) );
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
      if( fabs(XLALGetREAL8MatrixValue( posdef, i, j ) -
            XLALGetREAL8MatrixValue( nonposdef, i, j )) <= LAL_REAL8_EPS ){
        XLALSetREAL8MatrixValue( posdef, i, j, 
          XLALGetREAL8MatrixValue( nonposdef, i, j ) );
      }
    }
  }

  return posdef;
}

REAL8 XLALGetREAL8MatrixValue( REAL8Array *matrix, INT4 i, INT4 j ){
  REAL8 val=0.;

  /* check that matrix is not NULL */
  if( matrix == NULL ){
    fprintf(stderr, "Error... matrix is NULL!\n");
    exit(0);
  }

  /* check that matrix dimemsions are not NULL */
  if( matrix->dimLength == NULL ){
    fprintf(stderr, "Error... matrix dimensions not defined!\n");
    exit(0);
  }

  /* check that matrix is two dimensional */
  if( matrix->dimLength->length != 2 ){
    fprintf(stderr, "Error... matrix is not 2D!\n");
    exit(0);
  }

  /* check that i and j are in the range */
  if( (i >= (INT4)matrix->dimLength->data[0] || 
       j >= (INT4)matrix->dimLength->data[1]) && (i < 0 || j < 0) ){
    fprintf(stderr, "Error... i or j not in matrix range!\n");
    exit(0);
  }

  val = matrix->data[i*matrix->dimLength->data[0] + j];

  return val;
}

void XLALSetREAL8MatrixValue( REAL8Array *matrix, INT4 i, INT4 j, REAL8 val ){
  /* check that matrix is not NULL */
  if( matrix == NULL ){
    fprintf(stderr, "Error... matrix is NULL!\n");
    exit(0);
  }

  /* check that matrix dimemsions are not NULL */
  if( matrix->dimLength == NULL ){
    fprintf(stderr, "Error... matrix dimensions not defined!\n");
    exit(0);
  }

  /* check that matrix is two dimensional */
  if( matrix->dimLength->length != 2 ){
    fprintf(stderr, "Error... matrix is not 2D!\n");
    exit(0);
  }

  /* check that i and j are in the range */
  if( (i >= (INT4)matrix->dimLength->data[0] || 
       j >= (INT4)matrix->dimLength->data[1]) && (i < 0 || j < 0) ){
    fprintf(stderr, "Error... i or j not in matrix range!\n");
    exit(0);
  }

  matrix->data[i*matrix->dimLength->data[0] + j] = val;
}

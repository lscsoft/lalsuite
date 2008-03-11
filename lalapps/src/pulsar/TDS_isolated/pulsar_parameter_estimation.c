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

/* 
  Author:
  $Id$
*/

#include "pulsar_parameter_estimation.h"

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
  
  /*===================== GET AND SET THE INPUT PARAMETER ====================*/
  get_input_args(&inputs, argc, argv);
  
  /* if we want to output in verbose mode set global variable */ 
  if(inputs.verbose) verbose = 1;
  
  /* get the pulsar parameters */
  LALReadTEMPOParFile(&status, &pulsar, inputs.parFile);
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
    data = malloc(sizeof(DataStructure));
  }
  else{
    /* if we're doing the MCMC we need to store all the Bks */
    data = malloc(numDets*sizeof(DataStructure));
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
    
    /* create lookup table */
    if( verbose ) fprintf(stderr, "Creating look up table.\n");
    
    data[k].lookupTable = NULL;
    data[k].lookupTable = calloc(1, sizeof(DetRespLookupTable));
    data[k].lookupTable->lookupTable=NULL;
    detAndSource.pSource = &inputs.psr;
    detAndSource.pDetector = &detPos[i];
    
    /* create memory for the lookup table */
    data[k].lookupTable->lookupTable=NULL;
    data[k].lookupTable->lookupTable = calloc(inputs.mesh.psiRangeSteps, 
      sizeof(LALDetAMResponse *));
  
    for( j = 0 ; j < inputs.mesh.psiRangeSteps ; j++ ){
      data[k].lookupTable->lookupTable[j] = calloc(inputs.mesh.timeRangeSteps,
        sizeof(LALDetAMResponse));
    }
    
    data[k].lookupTable->psiSteps = inputs.mesh.psiRangeSteps;
    data[k].lookupTable->timeSteps = inputs.mesh.timeRangeSteps;
    
    /* create lookup table */  
    response_lookup_table(data[k].times->data[0], detAndSource,
      data[k].lookupTable);
    
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
    
        results = marginalise_posterior(singleLike, inputs.priors, inputs.mesh,
          output);

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
      perform_mcmc(&data[k], inputs, 1, dets[i]);
              
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
      
        results = marginalise_posterior(jointLike, inputs.priors, inputs.mesh,
          output);
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
      perform_mcmc(data, inputs, numDets, output.det);
    
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
    { 0, 0, 0, 0 }
  };
  
  CHAR args[] =
"hD:p:P:i:o:a:A:j:b:B:k:s:S:m:c:C:n:l:L:q:Q:U:u:Y:T:v:V:z:Z:e:E:d:I:x:t:H:w:W:\
y:g:K:N:X:O:" ;
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
REAL8 **** allocate_likelihood_memory(MeshGrid mesh){
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
        vars.psi = mesh.minVals.psi + (REAL8)k*mesh.delta.psi;
        
        /* perform final loop over h0 within log_likelihood function */
        noiseEvidence = log_likelihood(logLike[i][j][k], data, vars, mesh);
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


REAL8 log_likelihood(REAL8 *likeArray, DataStructure data,
  IntrinsicPulsarVariables vars, MeshGrid mesh){
  INT4 i=0, j=0, count=0, k=0;
  INT4 length=0, chunkLength=0;

  REAL8 tstart=0., T=0.;

  COMPLEX16 model;
  INT4 psibin=0, timebin=0;

  REAL8 plus=0., cross=0.;
  REAL8 sumModel=0., sumDataModel=0.;
  REAL8 chiSquare=0.;
  COMPLEX16 B;

  REAL8 exclamation[data.chunkMax+1]; /* all factorials up to chunkMax */
  REAL8 log2=0.;
  
  REAL8 noiseEvidence=0.; /* the log evidence that the data is just noise */

  INT4 first=0, through=0;
  
  /* to save time get all log factorials up to chunkMax */
  for( i = 0 ; i < data.chunkMax+1 ; i++ )
    exclamation[i] = log_factorial(i);
  
  log2 = log(2.);
  
  /* set the psi bin for the lookup table */
  psibin = (INT4)ROUND( ( vars.psi + LAL_PI/4. )
           *(REAL8)(data.lookupTable->psiSteps-1.)/LAL_PI_2 );
           
  length = (INT4)data.data->length + 1 - 
           data.chunkLengths->data[(INT4)data.chunkLengths->length-1];
           
  tstart = data.times->data[0]; /* time of first B_k */

  for( i = 0 ; i < length ; i += chunkLength ){
    chunkLength = data.chunkLengths->data[count];

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

      /* create the signal model */
      model.re = plus*vars.Xpcosphi_2 + cross*vars.Xcsinphi_2;
      model.im = plus*vars.Xpsinphi_2 - cross*vars.Xccosphi_2;

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
        likeArray[k] = ((REAL8)chunkLength - 1.)*log2;
      else likeArray[k] += ((REAL8)chunkLength - 1.)*log2;
        
      likeArray[k] += exclamation[chunkLength];
      likeArray[k] -= (REAL8)chunkLength*log(chiSquare);

      /* get the log evidence for the data not containing a signal */
      if( k == 0 ){
        noiseEvidence += ((REAL8)chunkLength - 1.)*log2;
        noiseEvidence += exclamation[chunkLength];
        noiseEvidence -= (REAL8)chunkLength*log(data.sumData->data[count]);
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
Results marginalise_posterior(REAL8 ****logPost, PriorVals prior, 
  MeshGrid mesh, OutputParams output){
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
  static LALStatus status;
  
  LALDetAMResponse response;
  
  LALGPSandAcc tgps;
  REAL8 T=0;
  
  INT4 i=0, j=0;
  
  tgps.accuracy = 1.0;
  
  for( i = 0 ; i < lookupTable->psiSteps ; i++ ){
    detAndSource.pSource->orientation = -(LAL_PI/4.) +
        (REAL8)i*(LAL_PI/2.) / ( (REAL8)lookupTable->psiSteps - 1. );
        
    for( j = 0 ; j < lookupTable->timeSteps ; j++ ){
      /* one day is 86400 seconds */
      T = t0 + (REAL8)j*86400./(REAL8)lookupTable->timeSteps;
    
      tgps.gps.gpsSeconds = (INT4)floor(T);
      tgps.gps.gpsNanoSeconds = (INT4)floor((fmod(T,1.0)*1.e9));
        
      LALComputeDetAMResponse(&status, &response, &detAndSource, &tgps);
      
      lookupTable->lookupTable[i][j].plus = response.plus;
      lookupTable->lookupTable[i][j].cross = response.cross;
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
  CHAR *det){
  IntrinsicPulsarVariables vars, varsNew;
  
  /* extra h0 and phase parameter required if a glitch is obsevered */
  IntrinsicPulsarVariables *extraVars=NULL, *extraVarsNew=NULL;
  REAL8 *glitchTimes=NULL;
  CHAR *gtimestr=NULL;
  INT4 nGlitches=input.mcmc.nGlitches;
  INT4 **g1=NULL, **g2=NULL; /* start - end positions of each glitch segment */
  DataStructure **glitchData=NULL;
  
  INT4 below0=0;
  
  REAL4Vector *randNum=NULL; /* LAL random variable params */
  UINT4 seed=0;              /* set to get seed from clock time */
  RandomParams *randomParams;
  
  CHAR *pos1=NULL, *pos2=NULL;
  INT4 i=0, j=0, k=0, n=0, count=0;
  
  REAL8 like1[1], like2[1];
  REAL8 logL1=0., logL2=0.; /* log likelihoods */
  REAL8 ratio;              /* logL2 - logL1 = log(L2/L1) */
  
  FILE *fp=NULL;
  CHAR outFile[256];
  
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
    g1 = calloc(numDets, sizeof(INT4 *));
    g2 = calloc(numDets, sizeof(INT4 *));
    
    for( j = 0 ; j < numDets ; j++ ){
        g1[j] = calloc(nGlitches + 1, sizeof(INT4));
        g2[j] = calloc(nGlitches + 1, sizeof(INT4));
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
        
        free(gtimestr);
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
         
        for( k = 0 ; k < data[j].times->length ; k++ ){
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
    glitchData = calloc(numDets, sizeof(DataStructure *));
    
    for( i = 0 ; i < numDets ; i++ ){
      glitchData[i] = calloc(nGlitches + 1, sizeof(DataStructure));
      
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
  fprintf(fp, "\n");
  
  /* create vector for random Gaussian numbers */
  randNum = XLALCreateREAL4Vector(4+2*nGlitches);
  
  count = 0;
  
  fprintf(stderr, "var.h0 = %le, extravars.h0 = %le.\n",
vars.h0, extraVars[0].h0);

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
        
        fprintf(fp, "\n");
      }
      
      count++;
      continue;
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
    for( k = 0 ; k < numDets ; k++ ){
      /* only calculate the likelhood twice in the first loop */
      if( i == 0 ){
        /* first likelihood */
        if( nGlitches == 0 ){
          input.mesh.minVals.h0 = vars.h0;
          
          log_likelihood(like1, data[k], vars, input.mesh);
          
          logL1 += *like1;
        }
        else{
          for( j = 0 ; j < nGlitches + 1 ; j++ ){
            if( j == 0){
              input.mesh.minVals.h0 = vars.h0;
              
              log_likelihood(like1, glitchData[k][j], vars, input.mesh);
            }
            else{
              input.mesh.minVals.h0 = extraVars[j-1].h0;
              log_likelihood(like1, glitchData[k][j], extraVars[j-1],
                input.mesh);
            }
            logL1 += *like1;
          }
        }
      }
        
      /* second likelihood - for new position in parameter space */
      if( nGlitches == 0 ){
        input.mesh.minVals.h0 = varsNew.h0;
        log_likelihood(like2, data[k], varsNew, input.mesh);
        logL2 += *like2;
      }
      else{
        for( j = 0 ; j < nGlitches + 1 ; j++ ){
          if( j == 0){
            input.mesh.minVals.h0 = varsNew.h0;
            log_likelihood(like2, glitchData[k][j], varsNew, input.mesh);
          }
          else{
            input.mesh.minVals.h0 = extraVarsNew[j-1].h0;
            log_likelihood(like2, glitchData[k][j], extraVarsNew[j-1],
              input.mesh);
          }  
          logL2 += *like2;            
        }
      }
    }
    
    /** FOR NOW JUST ALLOW UNIFORM PRIORS - SO WE JUST COMPARE LIKELIHOODS */
    /** DON'T CALCULATE THE POSTERIOR */
    
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
      
      for( j = 0 ; j < nGlitches ; j++)
        extraVars[j] = extraVarsNew[j];
        
      logL1 = logL2;     
    }
    /* otherwise accept with a certain probability */
    else if( log(XLALUniformDeviate(randomParams)) < ratio ){
      vars = varsNew;
      
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
  
  fclose(fp);
}


/* function to get the lengths of consecutive chunks of data */
void get_chunk_lengths(DataStructure data){
  INT4 i=0, j=0, count=0;

  /* create vector of data segment length */
  while( 1 ){
    count++; /* counter */
    
    /* break clause */
    if( i > data.data->length - 2 ){
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

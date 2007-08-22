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

#include <time.h>
#include <sys/timeb.h>

#include "pulsar_parameter_estimation.h"

int main(int argc, char *argv[]){
  static LALStatus status;
  
  double ****singleLike=NULL;
  double ****jointLike=NULL;
  
  InputParams inputs;
  BinaryPulsarParams pulsar;
  
  int numDets=0; /* number of detectors */
  char dets[5][3]; /* we'll have a max of five detectors */
  LALDetector detPos[5];
  
  int i=0, j=0, k=0;
  
  DataStructure *data=NULL;
  double times=0.;
  COMPLEX16 dataVals;
  
  FILE *fp=NULL;
  char dataFile[256];
  char outputFile[256];
  
  OutputParams output;
  double maxPost=0., evidence=0.;
  
  /*===================== GET AND SET THE INPUT PARAMETER ====================*/
  get_input_args(&inputs, argc, argv);
  
  /* get the pulsar parameters */
  LALReadTEMPOParFile(&status, &pulsar, inputs.parFile);
  inputs.mesh.pulsar.equatorialCoords.longitude = pulsar.ra;
  inputs.mesh.pulsar.equatorialCoords.latitude = pulsar.dec;	      
  inputs.mesh.pulsar.equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
  
  /* find the number of detectors being used */
  if( strstr(inputs.detectors, "H1") != NULL ){
    sprintf(dets[numDets], "H1");
    detPos[numDets] = lalCachedDetectors[LALDetectorIndexLHODIFF];
    numDets++;
  }
  if( strstr(inputs.detectors, "H2") != NULL ){
     sprintf(dets[numDets], "H2");
     detPos[numDets] = lalCachedDetectors[LALDetectorIndexLHODIFF];
     numDets++;
  }
  if( strstr(inputs.detectors, "L1") != NULL ){
     sprintf(dets[numDets], "L1");
     detPos[numDets] = lalCachedDetectors[LALDetectorIndexLLODIFF];
     numDets++;
  }
  if( strstr(inputs.detectors, "G1") != NULL ){
     sprintf(dets[numDets], "G1");
     detPos[numDets] = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
     numDets++;
  }
  if( strstr(inputs.detectors, "V1") != NULL ){
     sprintf(dets[numDets], "V1");
     detPos[numDets] = lalCachedDetectors[LALDetectorIndexVIRGODIFF];
     numDets++;
  }
  
  if( inputs.verbose ){
    fprintf(stderr, "Analysing data from %d detector(s):\n", numDets);
    for( i = 0 ; i < numDets ; i++ )
      fprintf(stderr, "%s\n", dets[i]);
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
    
    if( inputs.verbose ) 
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
    
    /* read in data */
    while(fscanf(fp, "%lf%lf%lf", &times, &dataVals.re, &dataVals.im) != EOF){
      /* exclude values smaller than 1e-28 as most are spurious points caused
         during a the heterodyne stage (e.g. when frame files were missing in
         the original S5 analysis) */
      if( fabs(dataVals.re) > 1.e-28 && fabs(dataVals.im) > 1.e-28 ){
        data[k].times->data[j] = times;
        data[k].data->data[j] = dataVals;
        j++;
      }
    }
    
    fclose(fp);
    
    if( inputs.verbose )
      fprintf(stderr, "I've read in the data for %s.\n", dets[i]);
    
    data[k].data = XLALResizeCOMPLEX16Vector(data[k].data, j);
    data[k].times = XLALResizeREAL8Vector(data[k].times, j);
    /*========================================================================*/
    
    output.det = dets[i];
    
    inputs.mesh.detector = detPos[i];
       
    /*======================== CALCULATE LIKELIHOOD ==========================*/
    create_likelihood_grid(data[k], singleLike, inputs.mesh);
    
    if( inputs.verbose )
      fprintf(stderr, "I've calculated the likelihood for %s.\n", dets[i]);
    
    /* if there's more than one detector calculate the joint likelihood */
    if( numDets > 1 ){
      /* add the single detector log likelihood onto the joint likelihood */
      combine_likelihoods(singleLike, jointLike, inputs.mesh);
      output.outPost = 0; /* don't output individual detector full posteriors */
    }
    /*========================================================================*/
    
    /*========== CREATE THE SINGLE DETECTOR POSTERIORS =======================*/
    maxPost = log_posterior(singleLike, inputs.priors, inputs.mesh, output);
    
    /* marginalise over each parameter and output the data */
    output.margParam = "h0"; /* create posterior for h0 */
    evidence = marginalise_posterior(singleLike, inputs.priors, inputs.mesh,
      output);
    
    output.margParam = "phi"; /* create posterior for phi0 */
    evidence = marginalise_posterior(singleLike, inputs.priors, inputs.mesh,
      output);
      
    output.margParam = "psi"; /* create posterior for psi */
    evidence = marginalise_posterior(singleLike, inputs.priors, inputs.mesh,
      output);
      
    output.margParam = "ciota"; /* create posterior for psi */
    evidence = marginalise_posterior(singleLike, inputs.priors, inputs.mesh,
      output);
    
    /* open file to output the evidence */
    if((fp = fopen(outputFile, "a"))==NULL){
      fprintf(stderr, "Error... can't open file %s!\n", outputFile);
      return 0;
    }
    
    /* output the evidence - FIXME: in future also output an UL if requested */
    fprintf(fp, "%s\t%le\n", output.det, evidence);
    fclose(fp);
    /*========================================================================*/
      
    if( inputs.mcmc.doMCMC == 0 ){
      XLALDestroyCOMPLEX16Vector(data[k].data);
      XLALDestroyREAL8Vector(data[k].times);
    }
  }
  
  /*=================== CREATE THE JOINT POSTERIOR IF REQUIRED ===============*/
  if( numDets > 1 ){
    output.outPost = inputs.outputPost; /* set for whether we want to output
                                            the full posterior */
    output.det = "Joint";
    
    maxPost = log_posterior(jointLike, inputs.priors, inputs.mesh, output);
    
     /* marginalise over each parameter and output the data */
    output.margParam = "h0"; /* create posterior for h0 */
    evidence = marginalise_posterior(singleLike, inputs.priors, inputs.mesh,
      output);
    
    output.margParam = "phi"; /* create posterior for phi0 */
    evidence = marginalise_posterior(singleLike, inputs.priors, inputs.mesh,
      output);
      
    output.margParam = "psi"; /* create posterior for psi */
    evidence = marginalise_posterior(singleLike, inputs.priors, inputs.mesh,
      output);
      
    output.margParam = "ciota"; /* create posterior for psi */
    evidence = marginalise_posterior(singleLike, inputs.priors, inputs.mesh,
      output);
    
    /* open file to output the evidence */
    if((fp = fopen(outputFile, "a"))==NULL){
      fprintf(stderr, "Error... can't open file %s!\n", outputFile);
      return 0;
    }
    
    /* output the evidence - FIXME: in future also output an UL if requested */
    fprintf(fp, "%s\t%le\n", output.det, evidence);
    fclose(fp);
  }
  /*==========================================================================*/
  
  /*====================== FREE THE LIKELIHOOD MEMORY ========================*/
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
  /*=========================================================================*/ 
                                            
  return 0;
}


/* function to combine log likelihoods to give a joint likelihood 
   log(p(data_joint|a) = log(p(data1|a)) + log(p(data2|a))        */
void combine_likelihoods(double ****logLike1, double ****logLike2, 
  MeshGrid mesh){
  int i=0, j=0, k=0, n=0;
  
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


/* function to get the input arguments from the command line */
void get_input_args(InputParams *inputParams, int argc, char *argv[]){
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
    { 0, 0, 0, 0 }
  };
  
  char args[] =
"hD:p:P:i:o:a:A:j:b:B:k:s:S:m:c:C:n:l:L:q:Q:U:u:Y:T:v:V:z:Z:e:E:d:I:x:t:H:w:W:\
y:g:K:N:" ;
  char *program = argv[0];
  
  /* set defaults */
  inputParams->mcmc.doMCMC = 0;/* by default don't perform an MCMC */
  inputParams->verbose = 0;    /* by default don't do verbose */
  inputParams->outputPost = 0; /* by default don't output the full posterior */
  inputParams->dob = 0.;       /* by default don't calculate an upper limit */
  
  inputParams->chunkMin = 5;   /* default to 5 minute minimum chunk length */
  inputParams->chunkMax = 30;  /* default to 30 minute maximum chunk length */
  
  /* default grid 50x50x50x50 = 6250000 points ~ 50 Mbytes */
  inputParams->mesh.minVals.h0 = 0.;          /* 0 <= h0 <= 1e-21; */
  inputParams->mesh.maxVals.h0 = 1.e-21;
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
  
  /* parse input arguments */
  while( 1 ){
    int option_index = 0;
    int c;
    
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
      inputParams->mesh.minVals.h0)/(double)(inputParams->mesh.h0Steps - 1.);
  }
  else
    inputParams->mesh.delta.h0 = 1.;
  
  if( inputParams->mesh.phiSteps > 1 ){
    inputParams->mesh.delta.phi0 = (inputParams->mesh.maxVals.phi0 -
      inputParams->mesh.minVals.phi0)/(double)(inputParams->mesh.phiSteps - 1.);
  }
  else
    inputParams->mesh.delta.phi0 = 1.;
  
  if( inputParams->mesh.psiSteps > 1 ){
    inputParams->mesh.delta.psi = (inputParams->mesh.maxVals.psi -
      inputParams->mesh.minVals.psi)/(double)(inputParams->mesh.psiSteps - 1.);
  }
  else
    inputParams->mesh.delta.psi = 1.;
  
  if( inputParams->mesh.ciotaSteps > 1 ){
    inputParams->mesh.delta.ci = (inputParams->mesh.maxVals.ci -
      inputParams->mesh.minVals.ci)/(double)(inputParams->mesh.ciotaSteps - 1.);
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
double **** allocate_likelihood_memory(MeshGrid mesh){
  int i=0, j=0, k=0;
  double ****logLike=NULL;
  
  /* allocate the h0 positions using calloc (i.e. array will be initialise to 
     zero */
  logLike = calloc(mesh.phiSteps, sizeof(double ***));

  for( i = 0 ; i < mesh.phiSteps ; i++ ){
    logLike[i] = calloc(mesh.ciotaSteps, sizeof(double **));
    
    for( j = 0 ; j < mesh.ciotaSteps ; j++ ){
      logLike[i][j] = calloc(mesh.psiSteps, sizeof(double *));
      
      for( k = 0 ; k < mesh.psiSteps ; k++ ){
        logLike[i][j][k] = calloc(mesh.h0Steps, sizeof(double));
      }
    }
  }
  
  return logLike;
}


/* function to create a log likelihood array over the parameter grid */ 
void create_likelihood_grid(DataStructure data, double ****logLike, 
  MeshGrid mesh){
  IntrinsicPulsarVariables vars;
  DetRespLookupTable lookupTable;
  
  LALDetAndSource detAndSource;
  int i=0, j=0, k=0, n=0;
  
  double cosphi=0., sinphi=0.;

  struct timeb t1, t2;
  
  detAndSource.pSource = &mesh.pulsar;
  detAndSource.pDetector = &mesh.detector;
  
  /* create memory for the lookup table */
  lookupTable.lookupTable = calloc(mesh.psiRangeSteps, 
    sizeof(LALDetAMResponse *));
  
  for( i=0 ; i < mesh.psiRangeSteps ; i++ ){
    lookupTable.lookupTable[i] = calloc(mesh.timeRangeSteps,
      sizeof(LALDetAMResponse));
  }
  
  lookupTable.psiSteps = mesh.psiRangeSteps;
  lookupTable.timeSteps = mesh.timeRangeSteps;
   
  /* create lookup table */  
  response_lookup_table(data.times->data[0], detAndSource, lookupTable);

  /* create vector of data segment length */
  data.chunkLengths = NULL;
  data.chunkLengths = XLALCreateINT4Vector(
    (int)ceil(data.data->length/data.chunkMin) );
  
  i=0;
  
  while( 1 ){
    n++; /* counter */ 

    /* break clause */
    if( i > (int)data.data->length - 2 ){
      /* set final value of chunkLength */
      data.chunkLengths->data[j] = n;
      j++;
      break;
    }

    i++;
    
    /* if consecutive points are within 180 seconds of each other count as in
       the same chunk */
    if( data.times->data[i] - data.times->data[i-1] > 180. || n ==
      data.chunkMax ){
      data.chunkLengths->data[j] = n;
      n = 0; /* reset counter */
           
      j++;
    }
  }

  data.chunkLengths = XLALResizeINT4Vector(data.chunkLengths, j);

  /* allocate memory for summations */
  data.sumData = NULL;
  data.sumData = XLALCreateREAL8Vector(j);
  
  /* get the sum over the data */
  sum_data(data);
  
  fprintf(stderr, "phiSteps = %d, ciSteps = %d, psiSteps = %d, h0Steps = \
%d.\n", mesh.phiSteps, mesh.ciotaSteps, mesh.psiSteps , mesh.h0Steps);
  
  /* calculate likelihood array */
  for( i = 0 ; i < mesh.phiSteps ; i++ ){
    vars.phi0 = mesh.minVals.phi0 + (double)i*mesh.delta.phi0;
    cosphi = cos(vars.phi0);
    sinphi = sin(vars.phi0);

    ftime(&t1);
    
    for( j = 0 ; j < mesh.ciotaSteps ; j++ ){
      vars.ci = mesh.minVals.ci + (double)j*mesh.delta.ci;
      vars.Xplus = 0.5*(1.+vars.ci*vars.ci);
      vars.Xcross = vars.ci;
      vars.Xpsinphi_2 = 0.5*vars.Xplus*sinphi;
      vars.Xcsinphi_2 = 0.5*vars.Xcross*sinphi;
      vars.Xpcosphi_2 = 0.5*vars.Xplus*cosphi;
      vars.Xccosphi_2 = 0.5*vars.Xcross*cosphi;

      for( k = 0 ; k < mesh.psiSteps ; k++ ){
        vars.psi = mesh.minVals.psi + (double)k*mesh.delta.psi;
        
        /* perform final loop over h0 within log_likelihood function */
        log_likelihood(logLike[i][j][k], data, vars, mesh, lookupTable);
      }
    }

    ftime(&t2);
    fprintf(stderr, "time = %ld\n", 1000*(t2.time - t1.time) + (t2.millitm -
t1.millitm));
    fprintf(stderr, "In h0 loop %d of %d.\n", i+1, mesh.h0Steps);
  }

  /* free memory */ 
  XLALDestroyREAL8Vector(data.sumData);
}


/* a function to sum over the data */
void sum_data(DataStructure data){
  int chunkLength=0, length=0, i=0, j=0, count=0;
  COMPLEX16 B;

  length = (int)data.data->length + 1 - 
           data.chunkLengths->data[(int)data.chunkLengths->length-1];
  
  for( i = 0 ; i < length ; i+= chunkLength ){
    chunkLength = (int)data.chunkLengths->data[count];
    data.sumData->data[count] = 0.;

    for( j = i ; j < i + chunkLength ; j++){
      B.re = data.data->data[j].re;
      B.im = data.data->data[j].im;

      /* sum up the data */
      data.sumData->data[count] += (double)(B.re*B.re + B.im*B.im);
    }
    
    count++;
  }
}


void log_likelihood(double *likeArray, DataStructure data,
  IntrinsicPulsarVariables vars, MeshGrid mesh, DetRespLookupTable lookupTable){
  int i=0, j=0, count=0, k=0;
  int length=0, chunkLength=0;

  double tstart=0., T=0.;

  COMPLEX16 model;
  int psibin=0, timebin=0;

  double plus=0., cross=0.;
  double sumModel=0., sumDataModel=0.;
  double chiSquare=0.;
  COMPLEX16 B;

  double exclamation[data.chunkMax+1]; /* all factorials up to chunkMax */
  double log2=0.;
  
  /* to save time get all log factorials up to chunkMax */
  for( i = 0 ; i < data.chunkMax+1 ; i++ )
    exclamation[i] = log((double)factorial(i));
  
  log2 = log(2.);
    
  /* set the psi bin for the lookup table */
  psibin = (int)ROUND( ( vars.psi + LAL_PI/4. )
           *(double)(lookupTable.psiSteps-1.)/LAL_PI_2 );
           
  length = (int)data.data->length + 1 - 
           data.chunkLengths->data[(int)data.chunkLengths->length-1];

  tstart = data.times->data[0]; /* time of first B_k */

  for( i = 0 ; i < length ; i += chunkLength ){
    chunkLength = data.chunkLengths->data[count];

    if( chunkLength < data.chunkMin ){
      count++;
      continue;
    }
    
    sumModel = 0.;
    sumDataModel = 0.;
    
    for( j = i ; j < i + chunkLength ; j++){
      /* set the time bin for the lookup table */
      T = fmod(data.times->data[j] - tstart, 86400.);
      timebin = (int)fmod(ROUND(T*(double)lookupTable.timeSteps/86400.),
        lookupTable.timeSteps);

      plus = lookupTable.lookupTable[psibin][timebin].plus;
      cross = lookupTable.lookupTable[psibin][timebin].cross;
      
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
      vars.h0 = mesh.minVals.h0 + (double)k*mesh.delta.h0;
      
      chiSquare = data.sumData->data[count];
      chiSquare -= 2.*vars.h0*sumDataModel;
      chiSquare += vars.h0*vars.h0*sumModel;
      
      likeArray[k] += ((double)chunkLength - 1.)*log2;
      likeArray[k] += exclamation[chunkLength];
      likeArray[k] -= (double)chunkLength*log(chiSquare);
    }
    
    count++;
  }
}


/* calculate the log prior */
double log_prior(PriorVals prior, MeshGrid mesh){
  double pri=1.;
  
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
double log_posterior(double ****logLike, PriorVals prior, MeshGrid mesh,
  OutputParams output){
  double maxPost=-1.e200, mP=0.;
  double logPi=0.;
  
  int i=0, j=0, k=0, n=0;
  
  char outputFile[256];
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
    prior.vars.phi0 = mesh.minVals.phi0 + (double)i*mesh.delta.phi0;
    
    for( j = 0 ; j < mesh.ciotaSteps ; j++ ){
      prior.vars.ci = mesh.minVals.ci + (double)j*mesh.delta.ci;
      
      for( k = 0 ; k < mesh.psiSteps ; k++ ){
        prior.vars.psi = mesh.minVals.psi + (double)k*mesh.delta.psi;
        
        for( n = 0 ; n < mesh.h0Steps ; n++ ){
          prior.vars.h0 = mesh.minVals.h0 + (double)n*mesh.delta.h0;
          
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
double marginalise_posterior(double ****logPost, PriorVals prior, 
  MeshGrid mesh, OutputParams output){
  double dval1=0., dval2=0., dval3=0., dval4=0.;
 
  double post1=0., post2=0.;
  
  double ***evSum1=NULL; /* first integral */
  double **evSum2=NULL;  /* second integral */
  double *evSum3=NULL;   /* third integral */
  double evSum4=-1.e200;      /* fouth integral */
  double maxPost=-1.e200;
  double evVal=0., sumVal=0.;
  
  int numSteps1=0, numSteps2=0, numSteps3=0, numSteps4=0;
  double step=0., minVal=0., cumsum=0.;
  
  int i=0, j=0, k=0, n=0;
  
  char outputFile[256];  /* filname to output the marginalised pdf */
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
  evSum1 = calloc(numSteps1, sizeof(double **));
  
  /* perform first integral */
  for( i = 0 ; i < numSteps1 ; i++ ){
    evSum1[i] = calloc(numSteps2, sizeof(double *));
    for( j = 0 ; j < numSteps2 ; j++ ){
      evSum1[i][j] = calloc(numSteps3, sizeof(double));
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
//           for( n = 0 ; n < numSteps4 ; n++ ){
//              if( strcmp( output.margParam, "h0" ) == 0 ){
//               post1 = logPost[j][k][n][i];
//   
//               /* use trapezium rule for intergration */
//               evVal = evSum1[i][j][k];
//               sumVal = post1 + log(dval4);
//               evSum1[i][j][k] = PLUS(sumVal, evVal);
//             }
//             else if( strcmp( output.margParam, "phi" ) == 0 ){
//               post1 = logPost[i][j][k][n];
//   
//               evVal = evSum1[i][j][k];
//               sumVal = post1 + log(dval4);
//               evSum1[i][j][k] = PLUS(sumVal, evVal);
//             }
//             else if( strcmp( output.margParam, "psi" ) == 0 ){
//               post1 = logPost[j][k][i][n];
//   
//               evVal = evSum1[i][j][k];
//               sumVal = post1 + log(dval4);
//               evSum1[i][j][k] = PLUS(sumVal, evVal);
//             }
//             else if( strcmp( output.margParam, "ciota" ) == 0 ){
//               post1 = logPost[j][i][k][n];
//   
//               evVal = evSum1[i][j][k];
//               sumVal = post1 + log(dval4);
//               evSum1[i][j][k] = PLUS(sumVal, evVal);
//             }
//           }
        }
      }
    }
  }
  
  /* allocate memory for second integral */
  evSum2 = calloc(numSteps1, sizeof(double *));
  
  /* perform the second integration */
  for( i = 0 ; i < numSteps1 ; i++ ){
    evSum2[i] = calloc(numSteps2, sizeof(double));
    for( j = 0 ; j < numSteps2 ; j++ ){
      evSum2[i][j] = -1.e200;
      
      if( numSteps3 == 1 ) evSum2[i][j] = evSum1[i][j][0];
      else{  
        for( k = 0 ; k <numSteps3 - 1 ; k++ ){
          evVal = evSum2[i][j];
          sumVal = log_trapezium(evSum1[i][j][k], evSum1[i][j][k+1], dval3);
          evSum2[i][j] = PLUS(sumVal, evVal);
        }
//         for( k = 0 ; k <numSteps3 ; k++ ){
//           evVal = evSum2[i][j];
//           sumVal = evSum1[i][j][k] + log(dval3);
//           evSum2[i][j] = PLUS(sumVal, evVal);
//         }
      }
    }
  }
  
  /* allocate memory for third integration */
  evSum3 = calloc(numSteps1, sizeof(double));
  
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
//       for( j = 0 ; j < numSteps2 ; j++ ){
//         evVal = evSum3[i];
//         sumVal = evSum2[i][j] + log(dval2);
//         evSum3[i] = PLUS(sumVal, evVal);
//       }
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
//     for( i = 0 ; i < numSteps1 ; i++ ){
//       evVal = evSum4;
//       sumVal = evSum3[i] + log(dval1);
//       evSum4 = PLUS(sumVal, evVal);
//     }
  }
  
  /* create the file to output the marginalise posterior to */
  sprintf(outputFile, "%s/pdf_%s.%s.%s", output.outputDir, output.margParam,
    output.psr, output.det);
  if( ( fp = fopen(outputFile, "w") ) == NULL ){
    fprintf(stderr, "Error... Can't open file to output pdf for %s.\n!",
      output.margParam);
    exit(0);
  }
  
  /* normalise the marginalised pdf using the evidence and free memory */
  for( i = 0 ; i < numSteps1 ; i++ ){
    for( j = 0 ; j < numSteps2 ; j++ ){
      free(evSum1[i][j]);
    }
    
    /* the parameter value at a given point */
    step = minVal + (double)i*dval1;
  
    /* calculate the cumulative probability */
    if( numSteps1 == 1 ) cumsum = exp(evSum3[0] - evSum4);
    else{
      if( i == 0 ){
        cumsum = 0.;
      }
      else{
        cumsum += exp(log_trapezium(evSum3[i-1], evSum3[i], dval1) - evSum4);
      }
    }
    /* print out marginalised posterior */
    fprintf(fp, "%le\t%le\t%le\n", step, exp(evSum3[i] - evSum4), cumsum);
    
    free(evSum1[i]);
    free(evSum2[i]);
  }
  
  fclose(fp);

  free(evSum1);
  free(evSum2);
  free(evSum3);

  return evSum4; /* return the log evidence */
}


/* function to do the trapezium rule for integration on logs  */
double log_trapezium(double logHeight1, double logHeight2, double width){
  double area=0.;
  
  /* area = 0.5*(height1 + height2)*width 
     logarea = log(0.5) + log(width) + log(exp(logHeight1) + exp(logHeight2)) */
  area = log(0.5) + log(width) + PLUS(logHeight1, logHeight2);
  
  return area;
}


/* detector response lookup table function  - this function will output a lookup
table of points in time and psi, covering a day from the start time (t0) and
from -pi/4 to pi/4 in psi */
void response_lookup_table(double t0, LALDetAndSource detAndSource,
  DetRespLookupTable lookupTable){
  static LALStatus status;
  
  LALDetAMResponse response;
  
  LALGPSandAcc tgps;
  double T=0;
  
  int i=0, j=0;
  
  tgps.accuracy = 1.0;
  
  for( i = 0 ; i < lookupTable.psiSteps ; i++ ){
    detAndSource.pSource->orientation = -(LAL_PI/4.) +
        (double)i*(LAL_PI/2.) / ( (double)lookupTable.psiSteps - 1. );
        
    for( j = 0 ; j < lookupTable.timeSteps ; j++ ){
      /* one day is 86400 seconds */
      T = t0 + (double)j*86400./(double)lookupTable.timeSteps;
    
      tgps.gps.gpsSeconds = (int)floor(T);
      tgps.gps.gpsNanoSeconds = (int)floor((fmod(T,1.0)*1.e9));
        
      LALComputeDetAMResponse(&status, &response, &detAndSource, &tgps);
      
      lookupTable.lookupTable[i][j].plus = response.plus;
      lookupTable.lookupTable[i][j].cross = response.cross;
    }
  }
}

/* function to return to factorial of an integer */
int factorial(int num){
  int fac=1, i=0;
  
  for( i=2 ; i <= num ; i++ ) fac *= i;
  
  return fac;
}

/** need to add code to calculate the upper limit for h0 if requested - this
will be called within marginalise posterior */

/** FIXME: I need to add the MCMC code */

/** The MCMC code needs to be able to include additional phase parameters for
potential phase jumps at the times of glitches */

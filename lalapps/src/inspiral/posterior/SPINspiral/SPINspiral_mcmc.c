/* 
   
   SPINspiral:                parameter estimation on binary inspirals detected by LIGO, including spins of the binary members
   SPINspiral_mcmc.c:         routines that form the MCMC core of the code
   
   
   Copyright 2007, 2008, 2009 Christian Roever, Marc van der Sluys, Vivien Raymond, Ilya Mandel
   
   
   This file is part of SPINspiral.
   
   SPINspiral is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   SPINspiral is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with SPINspiral.  If not, see <http://www.gnu.org/licenses/>.
   
*/



#include "SPINspiral.h"


/**
 * \file SPINspiral_mcmc.c
 * \brief Contains MCMC routines
 */



// ****************************************************************************************************************************************************  
/**
 * \brief MCMC routine - forms the MCMC core of the program
 *
 * Initialise and build a Markov chain
 */
// ****************************************************************************************************************************************************  
void MCMC(struct runPar run, struct interferometer *ifo[])
{
  
  struct parSet state;                        // MCMC/template parameter set struct
  
  // *** MCMC struct ***
  struct MCMCvariables mcmc;                  // MCMC variables struct
  copyRun2MCMC(run, &mcmc);                   // Copy elements from run struct to mcmc struct
  
  mcmc.tempOverlap = 1.1;                     // Set the overlap factor for the overlap between adjacent sinusoidal temperatures, e.g. 1.1.  1.0 means extremes touch. Keep <~ 1.3
  mcmc.decreaseSigma = 0.5;                   // Factor with which to decrease the jump-size sigma when a proposal is rejected
  mcmc.increaseSigma = exp(log(mcmc.decreaseSigma)* -(1.0/mcmc.acceptRateTarget - 1.0));  // (0.5)^(-(1/0.25 -1)) gives 8 -> jump up 1x needs 3 jumps down (0.5^3=1/8), so ~every 1:4 jumps gets accepted
  
  
  
  if(mcmc.beVerbose >= 1) {
    printf("\n");
    printf("   GPS base time:%12d,    target acceptance rate:%7.3f,    minimum logL:%9.1e\n\n",(int)mcmc.baseTime,mcmc.acceptRateTarget,mcmc.minlogL);
  }
  
  if(mcmc.parallelTempering==0) mcmc.nTemps=1;
  mcmc.ran = gsl_rng_alloc(gsl_rng_mt19937);  // GSL random-number seed
  gsl_rng_set(mcmc.ran, mcmc.seed);           // Set seed for this run
  
  
  char outfileName[99];
  char outfilePath[512];
  mcmc.fouts = (FILE**)calloc(mcmc.nTemps,sizeof(FILE*));
  for(mcmc.iTemp=0;mcmc.iTemp<mcmc.nTemps;mcmc.iTemp++) {
    if(mcmc.iTemp==0 || mcmc.saveHotChains>0) {
      if(run.outputPath) {
        strcpy(outfilePath,run.outputPath);
      } else {
        sprintf(outfilePath,"./");                       // In current dir, allows for multiple copies to run
      }
      sprintf(outfileName,"SPINspiral.output.%6.6d.%2.2d",mcmc.seed,mcmc.iTemp);
      strcat(outfilePath,outfileName);
      mcmc.fouts[mcmc.iTemp] = fopen(outfilePath,"w");
      if(mcmc.fouts[mcmc.iTemp] == NULL) {
        fprintf(stderr, "\n\n   ERROR:  could not open/create output file %s. Check that output directory %s exists.\n   Aborting...",outfilePath,run.outputPath); 
        exit(1);
      }
    }
  }
  
  if(mcmc.annealNburn0>=mcmc.annealNburn) {
    //fprintf(stderr, "\n ***  Warning: mcmc.annealNburn0 > mcmc.annealNburn, setting mcmc.annealNburn0 = mcmc.annealNburn*0.9 ***\n\n");
    mcmc.annealNburn0 = (int)(0.9*(double)mcmc.annealNburn);
  }
  
  
  
  // *** MEMORY ALLOCATION ********************************************************************************************************************************************************
  
  int i=0,j=0,j1=0,j2=0,injectionWF=0;
  
  //Allocate memory for (most of) the MCMCvariables struct
  allocateMCMCvariables(&mcmc);
  
  double **tempcovar;
  tempcovar = (double**)calloc(mcmc.nMCMCpar,sizeof(double*)); // A temp Cholesky-decomposed matrix
  for(i=0;i<mcmc.nMCMCpar;i++) tempcovar[i] = (double*)calloc(mcmc.nMCMCpar,sizeof(double));
  
  double ***covar1;
  covar1  = (double***)calloc(mcmc.nTemps,sizeof(double**)); // The actual covariance matrix
  for(i=0;i<mcmc.nTemps;i++) {
    covar1[i]  = (double**)calloc(mcmc.nMCMCpar,sizeof(double*));
    for(j=0;j<mcmc.nMCMCpar;j++) covar1[i][j]  = (double*)calloc(mcmc.nMCMCpar,sizeof(double));
  }
  
  
  
  
  
  // *** INITIALISE PARALLEL TEMPERING ********************************************************************************************************************************************
  
  // *** Set up temperature ladder ***
  if(mcmc.nTemps == 1) {
    mcmc.tempLadder[0] = 1.0;
  } else {
    setTemperatureLadderOld(&mcmc);
    //setTemperatureLadder(&mcmc);
  }
  mcmc.iTemp = 0;  //MUST be zero
  
  
  
  
  
  // *** WRITE RUN 'HEADER' TO SCREEN AND FILE ************************************************************************************************************************************
  
  writeMCMCheader(ifo, mcmc, run);
  mcmc.iTemp = 0;  //MUST be zero
  
  
  
  
  
  // *** INITIALISE MARKOV CHAIN **************************************************************************************************************************************************
  
  // *** Get the injection/best-guess values for signal ***
  if(mcmc.injectSignal >= 1) { // If a software injection was done:
    getInjectionParameters(&state, mcmc.nInjectPar, mcmc.injParVal);
    allocParset(&state, mcmc.networkSize);
    
    // *** Write injection/best-guess values to screen and file ***
    par2arr(state, mcmc.param, mcmc);  //Put the variables in their array
    injectionWF = 1;                                                 // Call localPar, netLogLikelihood with an injection waveform
    localPar(&state, ifo, mcmc.networkSize, injectionWF, run);
    mcmc.logL[mcmc.iTemp] = netLogLikelihood(&state, mcmc.networkSize, ifo, mcmc.injectionWaveform, injectionWF, run);  //Calculate the likelihood using the injection waveform
    freeParset(&state);
    
    // Store Injection parameters in temp array nParam[][]:
    for(i=0;i<mcmc.nInjectPar;i++) mcmc.nParam[mcmc.iTemp][i] = mcmc.param[mcmc.iTemp][i];
    
    // Copy Injection to MCMC parameters, as far as possible; otherwise use MCMC BestValue
    int iInj=0, nDiffPar=0;
    for(i=0;i<mcmc.nMCMCpar;i++) {
      iInj = mcmc.injRevID[mcmc.parID[i]];  //Get the index of this parameter in the injection set.  -1 if not available.
      if(mcmc.injParUse[mcmc.parID[i]] == 1) { // If an MCMC parameter was used for the injection
        mcmc.param[mcmc.iTemp][i] = mcmc.nParam[mcmc.iTemp][iInj];  // Set the MCMC parameter to the corresponding injection parameter
	
      } else { // If an MCMC parameter was not used for the injection, try to translate:
        if(mcmc.parID[i]==21 && mcmc.injID[i]==22) {
          mcmc.param[mcmc.iTemp][i] = exp(3.0*mcmc.nParam[mcmc.iTemp][i]);  // Injection uses log(d), MCMC uses d^3
          if(mcmc.beVerbose>=1) printf("   I translated  log(d_L/Mpc) = %lf  to  d_L^3 = %lf Mpc^3\n",mcmc.nParam[mcmc.iTemp][i],mcmc.param[mcmc.iTemp][i]);
        } else if(mcmc.parID[i]==22 && mcmc.injID[i]==21) {
          mcmc.param[mcmc.iTemp][i] = log(mcmc.nParam[mcmc.iTemp][i])/3.0;  // Injection uses d^3, MCMC uses log(d)
          if(mcmc.beVerbose>=1) printf("   I translated  d_L^3 = %lf  to  log(d_L/Mpc) = %lf Mpc^3\n",mcmc.nParam[mcmc.iTemp][i],mcmc.param[mcmc.iTemp][i]);
        } else {
          mcmc.param[mcmc.iTemp][i] = mcmc.parBestVal[i];        // Set the MCMC parameter to BestValue - this should only happen if the injection waveform has different parameters than the MCMC waveform
          nDiffPar += 1;
        }
	
      }
      
    } //for(i...)
    
    // Safety check:
    if(mcmc.mcmcWaveform == mcmc.injectionWaveform && nDiffPar != 0) {
      if(nDiffPar==1) {
        fprintf(stderr, "\n ***  Warning:  The injection and MCMC waveform are identical, but 1 parameter was found to be different ***\n\n");
      } else {
        fprintf(stderr, "\n ***  Warning:  The injection and MCMC waveform are identical, but %i parameters were found to be different ***\n\n",nDiffPar);
      }
    }
    
  } else { // If no software injection was done (injectSignal<=0):
    for(i=0;i<mcmc.nMCMCpar;i++) mcmc.param[mcmc.iTemp][i] = mcmc.parBestVal[i];  // Set the MCMC parameter to BestValue
  }
  
  
  // Print/save injection parameters as MCMC output, line -1:
  mcmc.iIter = -1;
  for(mcmc.iTemp=0;mcmc.iTemp<mcmc.nTemps;mcmc.iTemp++) {
    mcmc.iTemp = mcmc.iTemp;
    for(j1=0;j1<mcmc.nMCMCpar;j1++) mcmc.param[mcmc.iTemp][j1] = mcmc.param[0][j1];
    mcmc.logL[mcmc.iTemp] = mcmc.logL[0];
    writeMCMCoutput(mcmc, ifo);  //Write output line with injection parameters to screen and/or file (iteration -1)
  }
  mcmc.iTemp = 0;  //MUST be zero
  
  
  //Determine the number of parameters that is actually fitted/varied (i.e. not kept fixed at the true values)
  for(i=0;i<mcmc.nMCMCpar;i++) if(mcmc.parFix[i]==0) mcmc.nParFit += 1;
  
  
  // *** Initialise covariance matrix (initially diagonal), to do updates in the first block ***
  mcmc.corrUpdate[0] = 0;
  if(mcmc.correlatedUpdates>0) {
    for(j1=0;j1<mcmc.nMCMCpar;j1++) mcmc.covar[mcmc.iTemp][j1][j1] = mcmc.parSigma[j1];
    mcmc.corrUpdate[0] = 1; //Use the matrix above and don't change it
    if(mcmc.correlatedUpdates==2) mcmc.corrUpdate[0] = 2; //Use the matrix above and update it every nCorr iterations
  }
  
  
  
  
  // ***  GET (OFFSET) STARTING VALUES  ***********************************************************************************************************************************************
  
  // Get the best-guess values for the chain:
  getStartParameters(&state, run);
  allocParset(&state, mcmc.networkSize);
  
  par2arr(state, mcmc.param, mcmc);  //Put the variables in their array
  startMCMCOffset(&state, &mcmc, ifo, run);  // Start MCMC offset if and where wanted
  
  
  // *** Set the NEW array, sigma and scale ***
  for(i=0;i<mcmc.nMCMCpar;i++) {
    mcmc.nParam[mcmc.iTemp][i] = mcmc.param[mcmc.iTemp][i];
    mcmc.adaptSigma[mcmc.iTemp][i]   = 0.1  * mcmc.parSigma[i];
    if(mcmc.adaptiveMCMC==1) mcmc.adaptSigma[mcmc.iTemp][i] = mcmc.parSigma[i]; //Don't use adaptation (?)
    mcmc.adaptScale[mcmc.iTemp][i] = 10.0 * mcmc.parSigma[i];
    //mcmc.adaptScale[mcmc.iTemp][i] = 0.0 * mcmc.parSigma[i]; //No adaptation
    sigmaPeriodicBoundaries(mcmc.adaptSigma[mcmc.iTemp][i], i, mcmc);
  }
  
  
  
  
  
  // *** WRITE STARTING STATE TO SCREEN AND FILE **********************************************************************************************************************************
  
  arr2par(mcmc.param, &state, mcmc);                         //Get the parameters from their array
  injectionWF = 0;                                                 // Call netLogLikelihood with an MCMC waveform
  localPar(&state, ifo, mcmc.networkSize, injectionWF, run);
  mcmc.logL[mcmc.iTemp] = netLogLikelihood(&state, mcmc.networkSize, ifo, mcmc.mcmcWaveform, injectionWF, run);  //Calculate the likelihood
  
  // *** Write output line to screen and/or file
  printf("\n");
  mcmc.iIter = 0;
  for(mcmc.iTemp=0;mcmc.iTemp<mcmc.nTemps;mcmc.iTemp++) {
    mcmc.iTemp = mcmc.iTemp;
    for(j1=0;j1<mcmc.nMCMCpar;j1++) mcmc.param[mcmc.iTemp][j1] = mcmc.param[0][j1];
    mcmc.logL[mcmc.iTemp] = mcmc.logL[0];
    writeMCMCoutput(mcmc, ifo);  //Write output line to screen and/or file
  }
  mcmc.iTemp = 0;  //MUST be zero
  
  
  
  
  
  // *** INITIALISE PARALLEL TEMPERING ********************************************************************************************************************************************
  
  // *** Put the initial values of the parameters, sigmas etc in the different temperature chains ***
  if(mcmc.nTemps>1) {
    for(mcmc.iTemp=1;mcmc.iTemp<mcmc.nTemps;mcmc.iTemp++) {
      for(j=0;j<mcmc.nMCMCpar;j++) {
        mcmc.param[mcmc.iTemp][j] = mcmc.param[0][j];
        mcmc.nParam[mcmc.iTemp][j] = mcmc.nParam[0][j];
        mcmc.adaptSigma[mcmc.iTemp][j] = mcmc.adaptSigma[0][j];
        mcmc.adaptScale[mcmc.iTemp][j] = mcmc.adaptScale[0][j];
        mcmc.logL[mcmc.iTemp] = mcmc.logL[0];
        mcmc.nlogL[mcmc.iTemp] = mcmc.nlogL[0];
	
        for(j1=0;j1<mcmc.nMCMCpar;j1++) {
          for(j2=0;j2<=j1;j2++) mcmc.covar[mcmc.iTemp][j1][j2] = mcmc.covar[0][j1][j2];
        }
      }
      mcmc.corrUpdate[mcmc.iTemp] = mcmc.corrUpdate[0];
      //mcmc.corrUpdate[mcmc.iTemp] = 0; //Correlated update proposals only for T=1 chain?
      //mcmc.corrUpdate[mcmc.nTemps-1] = 0; //Correlated update proposals not for hottest chain
    }
  }
  
  
  
  
  
  
  
  
  // ********************************************************************************************************************************************************************************
  // ***  CREATE MARKOV CHAIN   *****************************************************************************************************************************************************
  // ********************************************************************************************************************************************************************************
  
  mcmc.iIter = 1;
  while(mcmc.iIter<=mcmc.nIter) {  // loop over Markov-chain states 
    
    for(mcmc.iTemp=0;mcmc.iTemp<mcmc.nTemps;mcmc.iTemp++) {  // loop over temperature ladder
      mcmc.iTemp = mcmc.iTemp;
      
      //Set temperature
      if(mcmc.parallelTempering==1 || mcmc.parallelTempering==3) { //Chains at fixed T
        mcmc.chTemp = mcmc.tempLadder[mcmc.iTemp];
      }
      if(mcmc.parallelTempering==2 || mcmc.parallelTempering==4) { //Chains with sinusoid T
        if(mcmc.iTemp==0) {
          mcmc.chTemp = 1.0;
        } else {
          mcmc.chTemp = mcmc.tempLadder[mcmc.iTemp]  +  mcmc.tempAmpl[mcmc.iTemp] * pow((-1.0),mcmc.iTemp) * sin(tpi*(double)mcmc.iIter/(5.0*(double)mcmc.nCorr));  //Sinusoid around the temperature T_i with amplitude tempAmpl and period 5.0 * nCorr
          mcmc.chTemp = max(mcmc.chTemp, 1.0);  // Make sure T>=1
        }
      }
      
      
      
      // *** UPDATE MARKOV CHAIN STATE **************************************************************************************************************************************************
      
      // *** Uncorrelated update *************************************************************************************************
      if(gsl_rng_uniform(mcmc.ran) > mcmc.corrFrac) {                                               //Do correlated updates from the beginning (quicker, but less efficient start); this saves ~4-5h for 2D, nCorr=1e4, nTemps=5
        if(gsl_rng_uniform(mcmc.ran) < mcmc.blockFrac){   
          uncorrelatedMCMCblockUpdate(ifo, &state, &mcmc, run);                                          //Block update for the current temperature chain
        } else {                                         
          uncorrelatedMCMCsingleUpdate(ifo, &state, &mcmc, run);                                         //Componentwise update for the current temperature chain (e.g. 90% of the time)
        }
	
        // *** Correlated update ****************************************************************************************************
      } else {
        correlatedMCMCupdate(ifo, &state, &mcmc, run);
      }
      
      
      // Update the dlogL = logL - logLo, and remember the parameter values where it has a maximum
      mcmc.dlogL[mcmc.iTemp] = mcmc.logL[mcmc.iTemp];
      if(mcmc.dlogL[mcmc.iTemp]>mcmc.maxdlogL[mcmc.iTemp]) {
        mcmc.maxdlogL[mcmc.iTemp] = mcmc.dlogL[mcmc.iTemp];
        for(i=0;i<mcmc.nMCMCpar;i++) mcmc.maxLparam[mcmc.iTemp][i] = mcmc.param[mcmc.iTemp][i];
      }
      
      
      // *** ACCEPT THE PROPOSED UPDATE *************************************************************************************************************************************************
      
      if(mcmc.acceptPrior[0]==1) { //Then write output and take care of the correlation matrix
	
	
        // *** WRITE STATE TO SCREEN AND FILE *******************************************************************************************************************************************
	
        writeMCMCoutput(mcmc, ifo);  //Write output line to screen and/or file
	
	
	
        // *** CORRELATION MATRIX *******************************************************************************************************************************************************
	
        //if(mcmc.corrUpdate[mcmc.iTemp]==2) {  //Calculate correlations only once
        if(mcmc.corrUpdate[mcmc.iTemp]>=2) { //Calculate correlations multiple times
	  
          // *** Save state to calculate correlations ***
          if(mcmc.iHist[mcmc.iTemp]<mcmc.nCorr) {
            for(j1=0;j1<mcmc.nMCMCpar;j1++) mcmc.hist[mcmc.iTemp][j1][mcmc.iHist[mcmc.iTemp]] = mcmc.param[mcmc.iTemp][j1];
            mcmc.iHist[mcmc.iTemp] += 1;
          }
	  
	  
	  
          // ***  Update covariance matrix  and  print parallel-tempering info  *************************************************************
          if(mcmc.iHist[mcmc.iTemp]>=mcmc.nCorr) {
	    
            updateCovarianceMatrix(&mcmc);  // Calculate the new covariance matrix for the current temperature chain and determine whether the matrix should be updated
	    
	    
            if(mcmc.parallelTempering>=1 && mcmc.prParTempInfo>0) writeChainInfo(mcmc);  //Print info on the current (temperature) chain(s) to screen
	    
	    
            mcmc.iHist[mcmc.iTemp] = 0;        //Reset history counter for the covariance matrix  This is also used for parallel tempering, that should perhaps get its own counter
	    
            /* 
               if( (mcmc.prMatrixInfo>0 || mcmc.prParTempInfo>0) && mcmc.iTemp==mcmc.nTemps-1 ) {
               printf("\n\n");
               //printf("\n%10s  %15s  %8s  %8s  %16s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s\n","cycle","logL","Mc","eta","tc","logdL","spin","kappa","RA","sindec","phase","sinthJ0","phiJ0","alpha");
               printf("\n%9s %10s  %7s %7s %8s %6s %6s %6s %6s %6s %6s %6s %6s %6s\n","cycle","logL","Mc","eta","tc","logdL","spin","kappa","RA","sindec","phase","snthJ0","phiJ0","alpha");
               }
            */
	    
          } //if(mcmc.iHist[mcmc.iTemp]>=mcmc.nCorr)
        } //if(mcmc.corrUpdate[mcmc.iTemp]>=2)
        // *** END CORRELATION MATRIX *************************************************************
	
      } //if(mcmc.acceptPrior[mcmc.iTemp]==1)
      
    } // for(mcmc.iTemp=0;mcmc.iTemp<mcmc.nTemps;mcmc.iTemp++) {  //loop over temperature ladder
    
    
    
    
    
    
    // *** ANNEALING ****************************************************************************************************************************************************************
    
    //Doesn't work with parallel tempering.  Use only when not using parallel tempering (and of course, temp0>1)
    if(mcmc.parallelTempering==0 && mcmc.annealTemp0>1.0) mcmc.chTemp = annealTemperature(mcmc.annealTemp0, mcmc.annealNburn, mcmc.annealNburn0, mcmc.iIter);
    
    
    // *** A test with adaptive parallel tempering was here.   
    // The commented-out, incomplete subroutine adaptive_parallel_tempering() was removed in rev.170
    
    
    
    
    // *** PARALLEL TEMPERING:  Swap states between T-chains *************************************************************************
    
    if(mcmc.acceptPrior[0]==1 && mcmc.parallelTempering>=1 && mcmc.nTemps>1) swapChains(&mcmc);
    
    
    if(mcmc.acceptPrior[0]==1) mcmc.iIter++;
  } // while(iIter<=mcmc.nIter) {  //loop over markov chain states 
  
  
  // ********************************************************************************************************************************************************************************
  // ***  END CREATE MARKOV CHAIN   *************************************************************************************************************************************************
  // ********************************************************************************************************************************************************************************
  
  
  
  
  
  for(mcmc.iTemp=0;mcmc.iTemp<mcmc.nTemps;mcmc.iTemp++) if(mcmc.iTemp==0 || mcmc.saveHotChains>0) fclose(mcmc.fouts[mcmc.iTemp]);
  free(mcmc.fouts);
  
  
  // *** FREE MEMORY **************************************************************************************************************************************************************
  
  printf("\n");
  freeMCMCvariables(&mcmc);
  
  
  for(i=0;i<mcmc.nMCMCpar;i++) free(tempcovar[i]);
  free(tempcovar);
  
  for(i=0;i<mcmc.nTemps;i++) {
    for(j=0;j<mcmc.nMCMCpar;j++) free(covar1[i][j]);
    free(covar1[i]);
  }
  free(covar1);
  
  freeParset(&state);
  
} // End MCMC()
// ****************************************************************************************************************************************************  















// ****************************************************************************************************************************************************  
/**
 * \brief Put the MCMC parameters in an array
 *
 */
// ****************************************************************************************************************************************************  
void par2arr(struct parSet par, double **param, struct MCMCvariables mcmc)
{
  int i=0;
  for(i=0;i<mcmc.nMCMCpar;i++) param[mcmc.iTemp][i] = par.par[i];
} // End par2arr
// ****************************************************************************************************************************************************  



// ****************************************************************************************************************************************************  
/**
 * \brief Get the MCMC parameters out of their array
 *
 */
// ****************************************************************************************************************************************************  
void arr2par(double **param, struct parSet *par, struct MCMCvariables mcmc)
{
  int i=0;
  for(i=0;i<mcmc.nMCMCpar;i++) par->par[i] = param[mcmc.iTemp][i];
} // End arr2par
// ****************************************************************************************************************************************************  













// ****************************************************************************************************************************************************  
/**
 * \brief Compute the prior for the given parameter set
 *
 * Contains boundary conditions and prior information for the MCMC.  Try to avoid returning 0, to increase jump sizes
 */
// ****************************************************************************************************************************************************  
double prior(double *par, int p, struct MCMCvariables mcmc)
{
  double priorValue = 1.0;
  
  if(mcmc.priorType[p]==21) {                                               // Periodic boundary condition to bring the variable between 0 and 2pi
    *par = fmod(*par+mtpi,tpi); 
  } else if(mcmc.priorType[p]==22) {                                        // Periodic boundary condition to bring the variable between 0 and pi
    *par = fmod(*par+mtpi,pi); 
  } else {                                                                  // Bounce back from the wall
    if(*par < mcmc.priorBoundLow[p] || *par > mcmc.priorBoundUp[p]) {                                        // Do only one bounce
      if(*par < mcmc.priorBoundLow[p]) {
        *par = mcmc.priorBoundLow[p] + fabs(*par - mcmc.priorBoundLow[p]);
      } else {
        *par = mcmc.priorBoundUp[p] - fabs(*par - mcmc.priorBoundUp[p]);
      }
      if(*par<mcmc.priorBoundLow[p] || *par>mcmc.priorBoundUp[p]) priorValue = 0.0;                             // If, after bouncing once, still outside the range, reject
    }
  }
  
  return priorValue;
} // End prior
// ****************************************************************************************************************************************************  



// ****************************************************************************************************************************************************  
/**
 * \brief Bring the adaptation sigma between its periodic boundaries
 *
 */
// ****************************************************************************************************************************************************  
double sigmaPeriodicBoundaries(double sigma, int p, struct MCMCvariables mcmc)
{
  if(mcmc.priorType[p] == 21) {
    return min(tpi,sigma);                                     //Bring sigma between 0 and 2pi;
  } else if(mcmc.priorType[p] == 22) {
    return min(pi,sigma);                                      //Bring sigma between 0 and pi;
  } else {
    return sigma;                                              //Don't do anything
  }
} // End sigmaPeriodicBoundaries()
// ****************************************************************************************************************************************************  








// ****************************************************************************************************************************************************  
/**
 * \brief Do a correlated block MCMC update
 *
 * Do an update for all non-fixed MCMC parameters. Use the covariance matrix to take into account correlations. 
 * The covariance matrix has been constructed from previous iterations.
 * Use adaptation to scale the whole matrix.
 * Some experiments with larger jumps using the 'hotter' covariance matrix.
 */
// ****************************************************************************************************************************************************  
void correlatedMCMCupdate(struct interferometer *ifo[], struct parSet *state, struct MCMCvariables *mcmc, struct runPar run)
// ****************************************************************************************************************************************************  
{
  int p1=0, p2=0, tempi=mcmc->iTemp, tempj=0;
  double temparr[mcmc->nMCMCpar], dparam=0.0;
  double ran=0.0, largejump1=0.0, largejumpall=0.0;
  
  //Prepare the proposal by creating a vector of univariate gaussian random numbers
  largejumpall = 1.0;
  ran = gsl_rng_uniform(mcmc->ran);
  if(ran < 1.0e-3) {
    largejumpall = 1.0e1;    //Every 1e3 iterations, take a 10x larger jump in all parameters
    if(ran < 1.0e-4) largejumpall = 1.0e2;    //Every 1e4 iterations, take a 100x larger jump in all parameters
  }
  
  for(p1=0;p1<mcmc->nMCMCpar;p1++) {
    largejump1 = 1.0;
    ran = gsl_rng_uniform(mcmc->ran);
    if(ran < 1.0e-2) {
      largejump1 = 1.0e1;    //Every 1e2 iterations, take a 10x larger jump in this parameter
      if(ran < 1.0e-3) largejump1 = 1.0e2;    //Every 1e3 iterations, take a 100x larger jump in this parameter
    }
    tempj = tempi;
    /*
      if(largejump1*largejumpall > 1.01) { //When making a larger jump, use the 'hotter' covariance matrix
      tempj = min(tempj+1,mcmc->nTemps);
      if(largejump1*largejumpall > 10.01) tempj = min(tempj+1,mcmc->nTemps);
      }
    */
    temparr[p1] = gsl_ran_gaussian(mcmc->ran,1.0) * mcmc->corrSig[tempj] * largejump1 * largejumpall;   //Univariate gaussian random numbers, with sigma=1, times the adaptable sigma_correlation times the large-jump factors
  }
  
  //Do the proposal
  mcmc->acceptPrior[tempi] = 1;
  for(p1=0;p1<mcmc->nMCMCpar;p1++){
    if(mcmc->parFix[p1]==0) {
      dparam = 0.0;
      for(p2=0;p2<=p1;p2++) dparam += mcmc->covar[tempi][p1][p2]*temparr[p2];     //Temparr is now a univariate gaussian random vector
      mcmc->nParam[tempi][p1] = mcmc->param[tempi][p1] + dparam;                  //Jump from the previous parameter value
      mcmc->adaptSigmaOut[tempi][p1] = fabs(dparam);                              //This isn't really sigma, but the proposed jump size
      mcmc->acceptPrior[tempi] *= (int)prior(&mcmc->nParam[tempi][p1],p1,*mcmc);
    }
  }
  
  
  
  /*
  //Testing with sky position/orientation updates
  if(gsl_rng_uniform(mcmc->ran) < 0.33) mcmc->nParam[tempi][6]  = fmod(mcmc->nParam[tempi][6]+pi,tpi);  //Move RA over 12h
  if(gsl_rng_uniform(mcmc->ran) < 0.33) mcmc->nParam[tempi][7]  *= -1.0;                                //Flip declination
  if(gsl_rng_uniform(mcmc->ran) < 0.33) mcmc->nParam[tempi][9]  *= -1.0;                                //Flip theta_Jo
  if(gsl_rng_uniform(mcmc->ran) < 0.33) mcmc->nParam[tempi][10] = fmod(mcmc->nParam[tempi][10]+pi,tpi); //Move phi_Jo over 12h
  */
  
  
  
  
  //Decide whether to accept
  if(mcmc->acceptPrior[tempi]==1) {                                            //Then calculate the likelihood
    arr2par(mcmc->nParam, state, *mcmc);                                       //Get the parameters from their array
    int injectionWF = 0;                                                 // Call netLogLikelihood with an MCMC waveform
    localPar(state, ifo, mcmc->networkSize, injectionWF, run);
    mcmc->nlogL[tempi] = netLogLikelihood(state, mcmc->networkSize, ifo, mcmc->mcmcWaveform, injectionWF, run); //Calculate the likelihood
    par2arr(*state, mcmc->nParam, *mcmc);                                      //Put the variables back in their array
    
    if(exp(max(-30.0,min(0.0,mcmc->nlogL[tempi]-mcmc->logL[tempi]))) > pow(gsl_rng_uniform(mcmc->ran),mcmc->chTemp) && mcmc->nlogL[tempi] > mcmc->minlogL) {  // Accept proposal
      for(p1=0;p1<mcmc->nMCMCpar;p1++) {
        if(mcmc->parFix[p1]==0) {
          mcmc->param[tempi][p1] = mcmc->nParam[tempi][p1];
          mcmc->accepted[tempi][p1] += 1;
        }
      }
      mcmc->logL[tempi] = mcmc->nlogL[tempi];
      if(mcmc->adaptiveMCMC==1){ 
        mcmc->corrSig[tempi] *= mcmc->increaseSigma;                    // Increase sigma
      }
    } else {                                                      // Reject proposal because of low likelihood
      if(mcmc->adaptiveMCMC==1){ 
        mcmc->corrSig[tempi] *= mcmc->decreaseSigma;                    // Decrease sigma
      }
    }
  } else {                                                        // Reject proposal because of boundary conditions.  Perhaps one should increase the step size, or at least not decrease it?
    /*
      if(mcmc->adaptiveMCMC==1){ 
      mcmc->corrSig[tempi] *= 0.8;
      }
    */
  }
} // End correlatedMCMCupdate
// ****************************************************************************************************************************************************  












// ****************************************************************************************************************************************************  
/**
 * \brief Do an uncorrelated, per-parameter MCMC update
 *
 * Do an update for all non-fixed MCMC parameters in the current T chain. 
 * Propose a jump and decide whether to accept or not on a per-parameter basis.
 * Use adaptation.
 */
// ****************************************************************************************************************************************************  
void uncorrelatedMCMCsingleUpdate(struct interferometer *ifo[], struct parSet *state, struct MCMCvariables *mcmc, struct runPar run)
// ****************************************************************************************************************************************************  
{
  int p=0, tempi=mcmc->iTemp;
  double gamma=0.0;
  double ran=0.0, largejump1=0.0, largejumpall=0.0;
  
  largejumpall = 1.0;
  ran = gsl_rng_uniform(mcmc->ran);
  if(ran < 1.0e-3) largejumpall = 1.0e1;    //Every 1e3 iterations, take a 10x larger jump in all parameters
  if(ran < 1.0e-4) largejumpall = 1.0e2;    //Every 1e4 iterations, take a 100x larger jump in all parameters
  
  for(p=0;p<mcmc->nMCMCpar;p++) if(mcmc->parFix[p]==0) mcmc->nParam[tempi][p] = mcmc->param[tempi][p];
  for(p=0;p<mcmc->nMCMCpar;p++){
    if(mcmc->parFix[p]==0) {
      largejump1 = 1.0;
      ran = gsl_rng_uniform(mcmc->ran);
      if(ran < 1.0e-2) largejump1 = 1.0e1;    //Every 1e2 iterations, take a 10x larger jump in this parameter
      if(ran < 1.0e-3) largejump1 = 1.0e2;    //Every 1e3 iterations, take a 100x larger jump in this parameter
      
      mcmc->nParam[tempi][p] = mcmc->param[tempi][p] + gsl_ran_gaussian(mcmc->ran,mcmc->adaptSigma[tempi][p]) * largejump1 * largejumpall;
      
      /*
      //Testing with sky position/orientation updates
      if(p==6  && gsl_rng_uniform(mcmc->ran) < 0.3) mcmc->nParam[tempi][6]  = fmod(mcmc->nParam[tempi][6]+pi,tpi);  //Move RA over 12h
      if(p==7  && gsl_rng_uniform(mcmc->ran) < 0.3) mcmc->nParam[tempi][7]  *= -1.0;                                //Flip declination
      if(p==9  && gsl_rng_uniform(mcmc->ran) < 0.3) mcmc->nParam[tempi][9]  *= -1.0;                                //Flip theta_Jo
      if(p==10 && gsl_rng_uniform(mcmc->ran) < 0.3) mcmc->nParam[tempi][10] = fmod(mcmc->nParam[tempi][10]+pi,tpi); //Move phi_Jo over 12h
      */
      
      mcmc->acceptPrior[tempi] = (int)prior(&mcmc->nParam[tempi][p],p,*mcmc);
      
      if(mcmc->acceptPrior[tempi]==1) {
        arr2par(mcmc->nParam, state, *mcmc);                                            //Get the parameters from their array
        int injectionWF = 0;                                                 // Call netLogLikelihood with an MCMC waveform
        localPar(state, ifo, mcmc->networkSize, injectionWF, run);
        mcmc->nlogL[tempi] = netLogLikelihood(state, mcmc->networkSize, ifo, mcmc->mcmcWaveform, injectionWF, run);   //Calculate the likelihood
        par2arr(*state, mcmc->nParam, *mcmc);                                           //Put the variables back in their array
	
        if(exp(max(-30.0,min(0.0,mcmc->nlogL[tempi]-mcmc->logL[tempi]))) > pow(gsl_rng_uniform(mcmc->ran),mcmc->chTemp) && mcmc->nlogL[tempi] > mcmc->minlogL) {  //Accept proposal
          mcmc->param[tempi][p] = mcmc->nParam[tempi][p];
          mcmc->logL[tempi] = mcmc->nlogL[tempi];
          if(mcmc->adaptiveMCMC==1){
            gamma = mcmc->adaptScale[tempi][p]*pow(1.0/((double)(mcmc->iIter+1)),1.0/6.0);
            mcmc->adaptSigma[tempi][p] = max(0.0,mcmc->adaptSigma[tempi][p] + gamma*(1.0 - mcmc->acceptRateTarget)); //Accept - increase sigma
            sigmaPeriodicBoundaries(mcmc->adaptSigma[tempi][p], p, *mcmc);              //Bring the sigma between 0 and 2pi
          }
          mcmc->accepted[tempi][p] += 1;
        } else {                                                                        //Reject proposal
          mcmc->nParam[tempi][p] = mcmc->param[tempi][p];
          if(mcmc->adaptiveMCMC==1){
            gamma = mcmc->adaptScale[tempi][p]*pow(1.0/((double)(mcmc->iIter+1)),1.0/6.0);
            mcmc->adaptSigma[tempi][p] = max(0.0,mcmc->adaptSigma[tempi][p] - gamma*mcmc->acceptRateTarget); //Reject - decrease sigma
            sigmaPeriodicBoundaries(mcmc->adaptSigma[tempi][p], p, *mcmc);              //Bring the sigma between 0 and 2pi
            //mcmc->adaptSigma[tempi][p] = max(0.01*mcmc->adaptSigma[tempi][p], mcmc->adaptSigma[tempi][p] - gamma*mcmc->acceptRateTarget);
          }
        }
      } else {  //If new state not within boundaries
        mcmc->nParam[tempi][p] = mcmc->param[tempi][p];
        if(mcmc->adaptiveMCMC==1) {
          gamma = mcmc->adaptScale[tempi][p]*pow(1.0/((double)(mcmc->iIter+1)),1.0/6.0);
          mcmc->adaptSigma[tempi][p] = max(0.0,mcmc->adaptSigma[tempi][p] - gamma*mcmc->acceptRateTarget);   //Reject - decrease sigma
          sigmaPeriodicBoundaries(mcmc->adaptSigma[tempi][p], p, *mcmc);                                     //Bring the sigma between 0 and 2pi
        }
      } //if(mcmc->acceptPrior[tempi]==1)
    } //if(mcmc->parFix[p]==0)
    mcmc->adaptSigmaOut[tempi][p] = mcmc->adaptSigma[tempi][p]; //Save sigma for output
  } //p
} // End uncorrelatedMCMCsingleUpdate
// ****************************************************************************************************************************************************  






// ****************************************************************************************************************************************************  
/**
 * \brief Do an uncorrelated block update
 *
 * Do an update for all non-fixed MCMC parameters in the current T chain. 
 * Propose a jump and decide whether to accept or not for all parameters at once.
 * No adaptation here, some experimenting with larger jumps every now and then.
 */
// ****************************************************************************************************************************************************  
void uncorrelatedMCMCblockUpdate(struct interferometer *ifo[], struct parSet *state, struct MCMCvariables *mcmc, struct runPar run)
// ****************************************************************************************************************************************************  
{
  int p=0;
  double ran=0.0, largejump1=0.0, largejumpall=0.0;
  
  largejumpall = 1.0;
  ran = gsl_rng_uniform(mcmc->ran);
  if(ran < 1.0e-3) largejumpall = 1.0e1;    //Every 1e3 iterations, take a 10x larger jump in all parameters
  if(ran < 1.0e-4) largejumpall = 1.0e2;    //Every 1e4 iterations, take a 100x larger jump in all parameters
  
  mcmc->acceptPrior[mcmc->iTemp] = 1;
  for(p=0;p<mcmc->nMCMCpar;p++){
    if(mcmc->parFix[p]==0) {
      largejump1 = 1.0;
      ran = gsl_rng_uniform(mcmc->ran);
      if(ran < 1.0e-2) largejump1 = 1.0e1;    //Every 1e2 iterations, take a 10x larger jump in this parameter
      if(ran < 1.0e-3) largejump1 = 1.0e2;    //Every 1e3 iterations, take a 100x larger jump in this parameter
      
      mcmc->nParam[mcmc->iTemp][p] = mcmc->param[mcmc->iTemp][p] + gsl_ran_gaussian(mcmc->ran,mcmc->adaptSigma[mcmc->iTemp][p]) * largejump1 * largejumpall;
      mcmc->acceptPrior[mcmc->iTemp] *= (int)prior(&mcmc->nParam[mcmc->iTemp][p],p,*mcmc);
    }
  }
  
  if(mcmc->acceptPrior[mcmc->iTemp]==1) {
    arr2par(mcmc->nParam, state, *mcmc);                                      //Get the parameters from their array
    int injectionWF = 0;                                                 // Call netLogLikelihood with an MCMC waveform
    localPar(state, ifo, mcmc->networkSize, injectionWF, run);                               //Calculate local variables
    mcmc->nlogL[mcmc->iTemp] = netLogLikelihood(state, mcmc->networkSize, ifo, mcmc->mcmcWaveform, injectionWF, run);  //Calculate the likelihood
    par2arr(*state, mcmc->nParam, *mcmc);                                     //Put the variables back in their array
    
    if(exp(max(-30.0,min(0.0,mcmc->nlogL[mcmc->iTemp]-mcmc->logL[mcmc->iTemp]))) > pow(gsl_rng_uniform(mcmc->ran),mcmc->chTemp) && mcmc->nlogL[mcmc->iTemp] > mcmc->minlogL){  //Accept proposal if L>Lo
      for(p=0;p<mcmc->nMCMCpar;p++){
        if(mcmc->parFix[p]==0) {
          mcmc->param[mcmc->iTemp][p] = mcmc->nParam[mcmc->iTemp][p];
          mcmc->accepted[mcmc->iTemp][p] += 1;
        }
      }
      mcmc->logL[mcmc->iTemp] = mcmc->nlogL[mcmc->iTemp];
    }
  }
} // End uncorrelatedMCMCblockUpdate
// ****************************************************************************************************************************************************  




















// ****************************************************************************************************************************************************  
/**
 * \brief Write MCMC header to screen and file
 *
 */
// ****************************************************************************************************************************************************  
void writeMCMCheader(struct interferometer *ifo[], struct MCMCvariables mcmc, struct runPar run)
// ****************************************************************************************************************************************************  
{
  int i=0, tempi=0;
  // *** Print run parameters to screen ***
  if(mcmc.offsetMCMC==0) printf("   Starting MCMC from the true initial parameters\n\n");
  if(mcmc.offsetMCMC>=1) printf("   Starting MCMC from offset initial parameters\n\n");
  
  // *** Open the output file and write run parameters in the header ***
  for(tempi=0;tempi<mcmc.nTemps;tempi++) {
    if(tempi==0 || mcmc.saveHotChains>0) {
      fprintf(mcmc.fouts[tempi], "  SPINspiral version:%8.2f\n\n",1.0);
      fprintf(mcmc.fouts[tempi], "%10s  %10s  %6s  %20s  %6s %8s   %6s  %8s  %10s  %12s  %9s  %9s  %8s\n",
              "nIter","Nburn","seed","null likelihood","Ndet","nCorr","nTemps","Tmax","Tchain","Network SNR","Waveform","pN order","Npar");
      fprintf(mcmc.fouts[tempi], "%10d  %10d  %6d  %20.10lf  %6d %8d   %6d%10d%12.1f%14.6f  %9i  %9.1f  %8i\n",
              mcmc.nIter,mcmc.annealNburn,mcmc.seed,0.0,run.networkSize,mcmc.nCorr,mcmc.nTemps,(int)mcmc.maxTemp,mcmc.tempLadder[tempi],run.netsnr,run.mcmcWaveform,run.mcmcPNorder,run.nMCMCpar);
      fprintf(mcmc.fouts[tempi], "\n%16s  %16s  %10s  %10s  %10s  %10s  %20s  %15s  %12s  %12s  %12s\n",
              "Detector","SNR","f_low","f_high","before tc","after tc","Sample start (GPS)","Sample length","Sample rate","Sample size","FT size");
      for(i=0;i<run.networkSize;i++) {
        fprintf(mcmc.fouts[tempi], "%16s  %16.8lf  %10.2lf  %10.2lf  %10.2lf  %10.2lf  %20.8lf  %15.7lf  %12d  %12d  %12d\n",
                ifo[i]->name,ifo[i]->snr,ifo[i]->lowCut,ifo[i]->highCut,ifo[i]->before_tc,ifo[i]->after_tc,
                ifo[i]->FTstart,ifo[i]->deltaFT,ifo[i]->samplerate,ifo[i]->samplesize,ifo[i]->FTsize);
      }
      
      //Parameter numbers:
      fprintf(mcmc.fouts[tempi], "\n\n%31s","");
      for(i=0;i<mcmc.nMCMCpar;i++) {
        if(mcmc.parID[i]>=11 && mcmc.parID[i]<=19) {  //GPS time
          fprintf(mcmc.fouts[tempi], " %17i",mcmc.parID[i]);
        } else {
          fprintf(mcmc.fouts[tempi], " %9i",mcmc.parID[i]);
        }
      }
      fprintf(mcmc.fouts[tempi],"\n");
      
      //Parameter symbols:
      fprintf(mcmc.fouts[tempi], "%8s %12s %9s","Cycle","log_Post.","Prior");
      for(i=0;i<mcmc.nMCMCpar;i++) {
        if(mcmc.parID[i]>=11 && mcmc.parID[i]<=19) {  //GPS time
          fprintf(mcmc.fouts[tempi], " %17s",mcmc.parAbrev[mcmc.parID[i]]);
        } else {
          fprintf(mcmc.fouts[tempi], " %9s",mcmc.parAbrev[mcmc.parID[i]]);
        }
      }
      fprintf(mcmc.fouts[tempi],"\n");
      
      fflush(mcmc.fouts[tempi]);
    }
  }
} // End writeMCMCheader
// ****************************************************************************************************************************************************  






// ****************************************************************************************************************************************************  
/**
 * \brief Write an MCMC iteration as an output line to screen and/or file
 *
 */
// ****************************************************************************************************************************************************  
void writeMCMCoutput(struct MCMCvariables mcmc, struct interferometer *ifo[])
// ****************************************************************************************************************************************************  
{
  int p=0, tempi=mcmc.iTemp, iIter=mcmc.iIter, i=0;
  
  //printf("%d  %d",tempi,iIter);
  
  // *** Write output to screen ***
  if(tempi==0) { //Only for the T=1 chain
    /*ILYA*/
    // if((iIter % (50*thinScreenOutput))==0 || iIter<0) printf("Previous iteration has match of %10g with true signal\n\n", 
    //  matchBetweenParameterArrayAndTrueParameters(mcmc.param[tempi], ifo, mcmc); //CHECK need support for two different waveforms
    // While the above is commented out, get rid of 'not used' warnings for the ifo struct:
    ifo[0]->index = ifo[0]->index;
    
    if((iIter % (50*mcmc.thinScreenOutput))==0 || iIter<0) {
      printf("\n%9s%10s","Cycle","logL");
      for(i=0;i<mcmc.nMCMCpar;i++) {
        if(mcmc.parID[i]>=11 && mcmc.parID[i]<=19) {  //GPS time
          printf(" %18s",mcmc.parAbrev[mcmc.parID[i]]);
        } else {
          printf(" %9s",mcmc.parAbrev[mcmc.parID[i]]);
        }
      }
      printf("\n");
    }
    
    
    if((iIter % mcmc.thinScreenOutput)==0 || iIter<0){  printf("%9d%10.3lf",iIter,mcmc.logL[tempi]);
      for(i=0;i<mcmc.nMCMCpar;i++) {
        if(mcmc.parID[i]>=11 && mcmc.parID[i]<=19) {  //GPS time
          printf(" %18.4f",mcmc.param[tempi][i]);
        } else {
          printf(" %9.4f",mcmc.param[tempi][i]);
        }
      }
      printf("\n");}
    
    
  }
  
  
  double *accrat;
  accrat = (double*)calloc(mcmc.nMCMCpar,sizeof(double));
  for(p=0;p<mcmc.nMCMCpar;p++) accrat[p] = 0.0;
  if(iIter > 0) {
    for(p=0;p<mcmc.nMCMCpar;p++) accrat[p] = mcmc.accepted[tempi][p]/(double)iIter;
  }
  
  // *** Write output to file ***
  if(tempi==0 || mcmc.saveHotChains>0) { //For all T-chains if desired, otherwise the T=1 chain only
    if((iIter % mcmc.thinOutput)==0 || iIter<=0){
      if(iIter<=0 || tempi==0 || (iIter % (mcmc.thinOutput*mcmc.saveHotChains))==0) { //Save every mcmc.thinOutput-th line for the T=1 chain, but every (mcmc.thinOutput*mcmc.saveHotChains)-th line for the T>1 ones
        fprintf(mcmc.fouts[tempi], "%8d %12.5lf %9.6lf", iIter,mcmc.logL[tempi],1.0);
        for(i=0;i<mcmc.nMCMCpar;i++) {
          if(mcmc.parID[i]>=11 && mcmc.parID[i]<=19) {  //GPS time
            fprintf(mcmc.fouts[tempi]," %17.6f",mcmc.param[tempi][i]);
          } else {
            fprintf(mcmc.fouts[tempi]," %9.5f",mcmc.param[tempi][i]);
          }
        }
        fprintf(mcmc.fouts[tempi],"\n");
	
	
        fflush(mcmc.fouts[tempi]); //Make sure any 'snapshot' you take halfway is complete
	
      } //if(tempi==0 || (iIter % (mcmc.thinOutput*mcmc.saveHotChains))==0)
    } //if((iIter % mcmc.thinOutput)==0 || iIter<0)
  } //if(tempi==0)
  
  free(accrat);
} // End writeMCMCoutput
// ****************************************************************************************************************************************************  









// ****************************************************************************************************************************************************  
/**
 * \brief Allocate memory for the MCMCvariables struct.
 *
 * Allocate memory for the MCMCvariables struct.  Don't forget to deallocate whatever you put here in freeMCMCvariables()
 */
// ****************************************************************************************************************************************************  
void allocateMCMCvariables(struct MCMCvariables *mcmc)
// ****************************************************************************************************************************************************  
{
  int i=0, j=0, nPar=0;
  
  nPar = max(mcmc->nMCMCpar, mcmc->nInjectPar);  // param[] and nParam[] may be used for injection parameters at initialisation...
  
  mcmc->nParFit=0;
  
  mcmc->histMean  = (double*)calloc(mcmc->nMCMCpar,sizeof(double));     // Mean of hist block of iterations, used to get the covariance matrix
  mcmc->histDev = (double*)calloc(mcmc->nMCMCpar,sizeof(double));       // Standard deviation of hist block of iterations, used to get the covariance matrix
  for(i=0;i<mcmc->nMCMCpar;i++) {
    mcmc->histMean[i] = 0;
    mcmc->histDev[i] = 0;
  }
  
  mcmc->corrUpdate = (int*)calloc(mcmc->nTemps,sizeof(int));
  mcmc->acceptElems = (int*)calloc(mcmc->nTemps,sizeof(int));
  for(i=0;i<mcmc->nTemps;i++) {
    mcmc->corrUpdate[i] = 0;
    mcmc->acceptElems[i] = 0;
  }
  
  mcmc->tempAmpl = (double*)calloc(mcmc->nTemps,sizeof(double));     // Temperature amplitudes for sinusoid T in parallel tempering                                                     
  mcmc->logL = (double*)calloc(mcmc->nTemps,sizeof(double));         // Current log(L)                                                                                          
  mcmc->nlogL = (double*)calloc(mcmc->nTemps,sizeof(double));        // New log(L)                                                                                                      
  mcmc->dlogL = (double*)calloc(mcmc->nTemps,sizeof(double));        // log(L)-log(Lo)                                                                                          
  mcmc->maxdlogL = (double*)calloc(mcmc->nTemps,sizeof(double));     // Remember the maximum dlog(L)                                                                                    
  for(i=0;i<mcmc->nTemps;i++) {
    mcmc->tempAmpl[i] = 0.0;
    mcmc->logL[i] = 0.0;
    mcmc->nlogL[i] = 0.0;
    mcmc->dlogL[i] = 0.0;
    mcmc->maxdlogL[i] = -1.e30;
  }
  
  mcmc->corrSig = (double*)calloc(mcmc->nTemps,sizeof(double));      // Sigma for correlated update proposals                                                
  mcmc->swapTs1 = (int*)calloc(mcmc->nTemps,sizeof(int));                  // Totals for the columns in the chain-swap matrix                                        
  mcmc->swapTs2 = (int*)calloc(mcmc->nTemps,sizeof(int));                  // Totals for the rows in the chain-swap matrix                                              
  mcmc->acceptPrior = (int*)calloc(mcmc->nTemps,sizeof(int));        // Check boundary conditions and choose to accept (1) or not(0)                                 
  mcmc->iHist = (int*)calloc(mcmc->nTemps,sizeof(int));            // Count the iteration number in the current history block to calculate the covar matrix from
  for(i=0;i<mcmc->nTemps;i++) {
    mcmc->corrSig[i] = 1.0;
    mcmc->swapTs1[i] = 0;
    mcmc->swapTs2[i] = 0;
    mcmc->acceptPrior[i] = 1;
    mcmc->iHist[i] = 0;
  }
  
  mcmc->accepted = (int**)calloc(mcmc->nTemps,sizeof(int*));         // Count accepted proposals
  mcmc->swapTss = (int**)calloc(mcmc->nTemps,sizeof(int*));          // Count swaps between chains
  mcmc->param = (double**)calloc(mcmc->nTemps,sizeof(double*));      // The old parameters for all chains
  mcmc->nParam = (double**)calloc(mcmc->nTemps,sizeof(double*));     // The new parameters for all chains
  mcmc->maxLparam = (double**)calloc(mcmc->nTemps,sizeof(double*));   // The best parameters for all chains (max logL)
  mcmc->adaptSigma = (double**)calloc(mcmc->nTemps,sizeof(double*));        // The standard deviation of the gaussian to draw the jump size from
  mcmc->adaptSigmaOut = (double**)calloc(mcmc->nTemps,sizeof(double*));     // The sigma that gets written to output
  mcmc->adaptScale = (double**)calloc(mcmc->nTemps,sizeof(double*));      // The rate of adaptation
  for(i=0;i<mcmc->nTemps;i++) {
    mcmc->accepted[i] = (int*)calloc(mcmc->nMCMCpar,sizeof(int));
    mcmc->swapTss[i] = (int*)calloc(mcmc->nTemps,sizeof(int));
    mcmc->param[i] = (double*)calloc(nPar,sizeof(double));
    mcmc->nParam[i] = (double*)calloc(nPar,sizeof(double));
    mcmc->maxLparam[i] = (double*)calloc(mcmc->nMCMCpar,sizeof(double));
    mcmc->adaptSigma[i] = (double*)calloc(mcmc->nMCMCpar,sizeof(double));
    mcmc->adaptSigmaOut[i] = (double*)calloc(mcmc->nMCMCpar,sizeof(double));
    mcmc->adaptScale[i] = (double*)calloc(mcmc->nMCMCpar,sizeof(double));
  }
  
  mcmc->hist    = (double***)calloc(mcmc->nTemps,sizeof(double**));  // Store a block of iterations, to calculate the covariances
  mcmc->covar  = (double***)calloc(mcmc->nTemps,sizeof(double**));   // The Cholesky-decomposed covariance matrix
  for(i=0;i<mcmc->nTemps;i++) {
    mcmc->hist[i]    = (double**)calloc(mcmc->nMCMCpar,sizeof(double*));
    mcmc->covar[i]  = (double**)calloc(mcmc->nMCMCpar,sizeof(double*));
    for(j=0;j<mcmc->nMCMCpar;j++) {
      mcmc->hist[i][j]    = (double*)calloc(mcmc->nCorr,sizeof(double));
      mcmc->covar[i][j]  = (double*)calloc(mcmc->nMCMCpar,sizeof(double));
    }
  }  
} // End allocateMCMCvariables
// ****************************************************************************************************************************************************  




// ****************************************************************************************************************************************************  
/**
 * \brief Deallocate memory for the MCMCvariables struct
 *
 */
// ****************************************************************************************************************************************************  
void freeMCMCvariables(struct MCMCvariables *mcmc)
// ****************************************************************************************************************************************************  
{
  int i=0, j=0;
  gsl_rng_free(mcmc->ran);
  
  free(mcmc->histMean);
  free(mcmc->histDev);
  
  free(mcmc->corrUpdate);
  free(mcmc->acceptElems);
  
  free(mcmc->tempAmpl);
  free(mcmc->logL);
  free(mcmc->nlogL);
  free(mcmc->dlogL);
  free(mcmc->maxdlogL);
  
  free(mcmc->corrSig);
  free(mcmc->acceptPrior);
  free(mcmc->iHist);
  free(mcmc->swapTs1);
  free(mcmc->swapTs2);
  
  for(i=0;i<mcmc->nTemps;i++) {
    free(mcmc->accepted[i]);
    free(mcmc->swapTss[i]);
    free(mcmc->param[i]);
    free(mcmc->nParam[i]);
    free(mcmc->maxLparam[i]);
    free(mcmc->adaptSigma[i]);
    free(mcmc->adaptSigmaOut[i]);
    free(mcmc->adaptScale[i]);
  }
  free(mcmc->accepted);
  free(mcmc->swapTss);
  free(mcmc->param);
  free(mcmc->nParam);
  free(mcmc->maxLparam);
  free(mcmc->adaptSigma);
  free(mcmc->adaptSigmaOut);
  free(mcmc->adaptScale);
  
  for(i=0;i<mcmc->nTemps;i++) {
    for(j=0;j<mcmc->nMCMCpar;j++) {
      free(mcmc->hist[i][j]);
      free(mcmc->covar[i][j]);
    }
    free(mcmc->hist[i]);
    free(mcmc->covar[i]);
  }
  free(mcmc->hist);
  free(mcmc->covar);
} // End freeMCMCvariables
// ****************************************************************************************************************************************************  












// ****************************************************************************************************************************************************  
/**
 * \brief Calculate the new covariance matrix and determine whether the new matrix should be accepted
 *
 * Compute the new covariance matrix for the current T chain. Determine by the 'improvement' of the new matrix whether it should be accepted.
 */
// ****************************************************************************************************************************************************  
void updateCovarianceMatrix(struct MCMCvariables *mcmc)
// ****************************************************************************************************************************************************  
{
  int i=0, j=0, i1=0, p1=0, p2=0;
  double **tempcovar;
  tempcovar = (double**)calloc(mcmc->nMCMCpar,sizeof(double*));           // A temp Cholesky-decomposed matrix
  for(i=0;i<mcmc->nMCMCpar;i++) tempcovar[i] = (double*)calloc(mcmc->nMCMCpar,sizeof(double));
  double ***covar1;
  covar1  = (double***)calloc(mcmc->nTemps,sizeof(double**));   // The actual covariance matrix
  for(i=0;i<mcmc->nTemps;i++) {
    covar1[i]  = (double**)calloc(mcmc->nMCMCpar,sizeof(double*));
    for(j=0;j<mcmc->nMCMCpar;j++) covar1[i][j]  = (double*)calloc(mcmc->nMCMCpar,sizeof(double));
  }
  
  
  //Calculate the mean
  for(p1=0;p1<mcmc->nMCMCpar;p1++){
    mcmc->histMean[p1]=0.0;
    for(i1=0;i1<mcmc->nCorr;i1++) mcmc->histMean[p1]+=mcmc->hist[mcmc->iTemp][p1][i1];
    mcmc->histMean[p1]/=((double)mcmc->nCorr);
  }
  
  //Calculate the standard deviation. Only for printing, not used in the code
  for(p1=0;p1<mcmc->nMCMCpar;p1++){
    mcmc->histDev[p1]=0.0;
    for(i1=0;i1<mcmc->nCorr;i1++) mcmc->histDev[p1] += (mcmc->hist[mcmc->iTemp][p1][i1]-mcmc->histMean[p1])*(mcmc->hist[mcmc->iTemp][p1][i1]-mcmc->histMean[p1]);
    mcmc->histDev[p1] = sqrt(mcmc->histDev[p1]/(double)(mcmc->nCorr-1));
  }
  
  //Calculate the covariances and save them in covar1
  for(p1=0;p1<mcmc->nMCMCpar;p1++){
    for(p2=0;p2<=p1;p2++){
      covar1[mcmc->iTemp][p1][p2]=0.0;
      for(i1=0;i1<mcmc->nCorr;i1++) covar1[mcmc->iTemp][p1][p2] += (mcmc->hist[mcmc->iTemp][p1][i1] - mcmc->histMean[p1]) * (mcmc->hist[mcmc->iTemp][p2][i1] - mcmc->histMean[p2]);
      covar1[mcmc->iTemp][p1][p2] /= (double)(mcmc->nCorr-1);
    }
  }
  
  //Store the covariance matrix in a temporary variable tempcovar, for Cholesky decomposition
  for(p1=0;p1<mcmc->nMCMCpar;p1++){
    for(p2=0;p2<=p1;p2++) tempcovar[p1][p2] = covar1[mcmc->iTemp][p1][p2];
  }
  
  if(mcmc->iTemp==0 && mcmc->prMatrixInfo>0) printf("\n\n");
  
  //Do Cholesky decomposition
  CholeskyDecompose(tempcovar,mcmc);
  
  //Get conditions to decide whether to accept the new matrix or not
  mcmc->acceptElems[mcmc->iTemp] = 0;
  for(p1=0;p1<mcmc->nMCMCpar;p1++) {
    if(mcmc->parFix[p1]==0) {
      if(tempcovar[p1][p1] < mcmc->covar[mcmc->iTemp][p1][p1]) mcmc->acceptElems[mcmc->iTemp] += 1; //Smaller diagonal element is better; count for how many this is the case
      if(tempcovar[p1][p1]<=0.0 || isnan(tempcovar[p1][p1])!=0 || isinf(tempcovar[p1][p1])!=0) mcmc->acceptElems[mcmc->iTemp] -= 9999;  //If diagonal element is <0, NaN or Inf
    }
  }
  mcmc->acceptElems[mcmc->iTemp] = max(mcmc->acceptElems[mcmc->iTemp],-1); //Now -1 means there is a diagonal element that is 0, NaN or Inf
  
  //Print matrix information  
  if(mcmc->prMatrixInfo==2 && mcmc->iTemp==0){
    printf("\n  Update for the covariance matrix proposed at iteration:  %10d\n",mcmc->iIter);
    printf("\n    AcceptElems: %d\n",mcmc->acceptElems[mcmc->iTemp]);
    printf("\n    Covariance matrix:\n");
    for(p1=0;p1<mcmc->nMCMCpar;p1++){
      for(p2=0;p2<=p1;p2++) printf("    %10.3g",covar1[mcmc->iTemp][p1][p2]);
      printf("\n");
    }
    printf("\n    Old Cholesky-decomposed matrix:\n");
    for(p1=0;p1<mcmc->nMCMCpar;p1++){
      for(p2=0;p2<=p1;p2++) printf("    %10.3g",mcmc->covar[mcmc->iTemp][p1][p2]);
      printf("\n");
    }
    printf("\n    New Cholesky-decomposed matrix:\n");
    for(p1=0;p1<mcmc->nMCMCpar;p1++){
      for(p2=0;p2<=p1;p2++) printf("    %10.3g",tempcovar[p1][p2]);
      //for(p2=0;p2<mcmc->nMCMCpar;p2++) printf("    %10.3g",tempcovar[p1][p2]);
      printf("\n");
    }
  }
  
  // Copy the new covariance matrix from tempcovar into mcmc->covar
  if(mcmc->acceptElems[mcmc->iTemp]>=0) { //Accept new matrix only if no Infs, NaNs, 0s occur.
    if(mcmc->corrUpdate[mcmc->iTemp]<=2 || (double)mcmc->acceptElems[mcmc->iTemp] >= (double)mcmc->nParFit*mcmc->matAccFr) { //Always accept the new matrix on the first update, otherwise only if the fraction matAccFr of diagonal elements are better (smaller) than before
      for(p1=0;p1<mcmc->nMCMCpar;p1++){
        for(p2=0;p2<=p1;p2++) mcmc->covar[mcmc->iTemp][p1][p2] = tempcovar[p1][p2]; 
      }
      mcmc->corrUpdate[mcmc->iTemp] += 1;
      if(mcmc->prMatrixInfo>0 && mcmc->iTemp==0) printf("  Proposed covariance-matrix update at iteration %d accepted.  AcceptElems: %d.  Accepted matrices: %d/%d \n", mcmc->iIter, mcmc->acceptElems[mcmc->iTemp], mcmc->corrUpdate[mcmc->iTemp]-2, (int)((double)mcmc->iIter/(double)mcmc->nCorr));  // -2 since you start with 2
    } else {
      if(mcmc->prMatrixInfo>0 && mcmc->iTemp==0) printf("  Proposed covariance-matrix update at iteration %d rejected.  AcceptElems: %d.  Accepted matrices: %d/%d \n", mcmc->iIter, mcmc->acceptElems[mcmc->iTemp], mcmc->corrUpdate[mcmc->iTemp]-2, (int)((double)mcmc->iIter/(double)mcmc->nCorr));  // -2 since you start with 2
    }
  } else {
    if(mcmc->prMatrixInfo>0 && mcmc->iTemp==0) printf("  Proposed covariance-matrix update at iteration %d rejected.  AcceptElems: %d.  Accepted matrices: %d/%d \n", mcmc->iIter, mcmc->acceptElems[mcmc->iTemp], mcmc->corrUpdate[mcmc->iTemp]-2, (int)((double)mcmc->iIter/(double)mcmc->nCorr));  // -2 since you start with 2
  }
  
  
  
  //Deallocate memory
  for(i=0;i<mcmc->nMCMCpar;i++) free(tempcovar[i]);
  free(tempcovar);
  
  for(i=0;i<mcmc->nTemps;i++) {
    for(j=0;j<mcmc->nMCMCpar;j++) free(covar1[i][j]);
    free(covar1[i]);
  }
  free(covar1);
  
} // End updateCovarianceMatrix
// ****************************************************************************************************************************************************  







// ****************************************************************************************************************************************************  
/**
 * \brief Compute Cholesky decomposition matrix
 *
 * Performs Cholesky decompositon of matrix A and returns result in the same matrix - adapted from PJG Fortran function
 * If matrix is not positive definite, return zeroes
 */
// ****************************************************************************************************************************************************  
void CholeskyDecompose(double **A, struct MCMCvariables *mcmc)
{
  int j1=0,j2=0,j3=0,notposdef=0;
  int n=mcmc->nMCMCpar;
  double sum=0.0;
  
  for(j1=0;j1<n;j1++){
    if(mcmc->parFix[j1]==0) {
      sum = A[j1][j1];
      for(j2=0;j2<j1;j2++){
        if(mcmc->parFix[j2]==0) {
          sum -= A[j1][j2]*A[j1][j2];
        }
      }
      if(sum<0.0) {
        notposdef=1;
      } else {
        A[j1][j1]=sqrt(sum);
        for(j2=j1+1;j2<n;j2++){
          if(mcmc->parFix[j2]==0){
            sum = A[j2][j1];
            for(j3=0;j3<j1;j3++){
              if(mcmc->parFix[j3]==0){
                sum -= A[j2][j3]*A[j1][j3];
              }
            }
          }
          A[j2][j1] = sum/A[j1][j1];
        }
      }
    }
  }
  if(notposdef==1) {
    //printf("  CholeskyDecompose():  Matrix %i is not positive definite\n",mcmc->iTemp);
    for(j1=0;j1<n;j1++){
      for(j2=0;j2<n;j2++){
        A[j1][j2] = 0.0;
      }
    }
  }
} // End CholeskyDecompose()
// ****************************************************************************************************************************************************  










// ****************************************************************************************************************************************************  
/**
 * \brief Annealing: set the temperature according to the iteration number and burnin
 */
// ****************************************************************************************************************************************************  
double annealTemperature(double temp0, int nburn, int nburn0, int iIter)
// ****************************************************************************************************************************************************  
{
  double temperature = 1.0;
  //temperature = temp0*pow(1.0/((double)(max(iIter,mcmc->nCorr))),1.0/10.0);
  //temperature = temp0*pow(((double)(iIter)),-0.1);
  //temperature = temp0*pow(((double)(iIter)),-0.25);
  //temperature = temp0 * pow(((double)(iIter)), (-log(temp0)/log((double)nburn)) );  //Temp drops from temp0 to 1.0 in nburn iterations
  //printf("%f\n",-log(temp0)/log((double)nburn));
  //temperature = min( max( temp0*pow( max((double)iIter-0.1*(double)nburn,1.0) ,-log10(temp0)/log10(0.9*(double)nburn)) , 1.0) , temp0);  //Temp stays at temp0 for 10% of burnin and then drops to 1.0 at the end of the burnin (too drastic for short burnins/high percentages?)
  //temperature = temp0*pow( max((double)iIter-1000.0,1.0) ,-log10(temp0)/log10((double)nburn-1000.0));  //Temp stays at temp0 for the first 1000 iterations and then drops to 1.0 at the end of the burnin (too drastic for short burnins/high percentages?)
  temperature = exp(log(temp0) * ((double)(nburn-iIter)/((double)(nburn-nburn0))));
  
  //if(gsl_rng_uniform(mcmc.ran)<1.e-2 && iIter > nburn && mcmc.dlogL[tempi] < 0.95*mcmc.maxdlogL[tempi]){
  //if(iIter > nburn && mcmc.dlogL[tempi] < 0.85*mcmc.maxdlogL[tempi]){
  //  temp=pow(temp0,0.5);
  //  mcmc.corrSig[tempi] *= 10.0;
  //}
  
  temperature = min( max(temperature,1.0) , temp0);  //Just in case...
  return temperature;
} // End annealTemperature
// ****************************************************************************************************************************************************  









// ****************************************************************************************************************************************************  
/**
 * \brief Parallel tempering: Swap states between two chains
 */
// ****************************************************************************************************************************************************  
void swapChains(struct MCMCvariables *mcmc)
// ****************************************************************************************************************************************************  
{
  int i=0, tempi=0, tempj=0;
  double tmpdbl = 0.0;
  
  //Swap parameters and likelihood between any two chains
  for(tempi=0;tempi<mcmc->nTemps-1;tempi++) {
    for(tempj=tempi+1;tempj<mcmc->nTemps;tempj++) {
      
      if(exp(max(-30.0,min(0.0, (1.0/mcmc->tempLadder[tempi]-1.0/mcmc->tempLadder[tempj]) * (mcmc->logL[tempj]-mcmc->logL[tempi]) ))) > gsl_rng_uniform(mcmc->ran)) { //Then swap...
        for(i=0;i<mcmc->nMCMCpar;i++) {
          tmpdbl = mcmc->param[tempj][i]; //Temp var
          mcmc->param[tempj][i] = mcmc->param[tempi][i];
          mcmc->param[tempi][i] = tmpdbl;
        }
        tmpdbl = mcmc->logL[tempj];
        mcmc->logL[tempj] = mcmc->logL[tempi];
        mcmc->logL[tempi] = tmpdbl;
        mcmc->swapTss[tempi][tempj] += 1;
        mcmc->swapTs1[tempi] += 1;
        mcmc->swapTs2[tempj] += 1;
      }
    } //tempj
  } //tempi
} // End swapChains
// ****************************************************************************************************************************************************  








// ****************************************************************************************************************************************************  
/**
 * \brief Parallel tempering: Print chain and swap info to screen
 */
// ****************************************************************************************************************************************************  
void writeChainInfo(struct MCMCvariables mcmc)
// ****************************************************************************************************************************************************  
{
  int tempi=mcmc.iTemp, p=0, t1=0, t2=0;
  double tmpdbl = 0.0;
  
  if(tempi==0) {
    printf("\n\n      Chain  log(T)   dlog(L)  AccEls AccMat     Swap  AccRat lgStdv:");
    for(p=0;p<mcmc.nMCMCpar;p++) printf(" %6s",mcmc.parAbrv[mcmc.parID[p]]);
    printf("\n");
  }
  
  printf("        %3d   %5.3f     %3d    %3d   %6.4f  %6.4f        ",
         tempi,log10(mcmc.chTemp),mcmc.acceptElems[tempi],mcmc.corrUpdate[tempi]-2,  (double)mcmc.swapTs1[tempi]/(double)mcmc.iIter,(double)mcmc.accepted[tempi][0]/(double)mcmc.iIter);
  for(p=0;p<mcmc.nMCMCpar;p++) {
    tmpdbl = log10(mcmc.histDev[p]+1.e-30);
    if(tmpdbl<-9.99) tmpdbl = 0.0;
    printf(" %6.3f",tmpdbl);
  }
  printf("\n");
  
  
  if(mcmc.prParTempInfo==2 && tempi==mcmc.nTemps-1) { //Print swap-rate matrix for parallel tempering
    printf("\n   Chain swap: ");
    for(t1=0;t1<mcmc.nTemps-1;t1++) printf("  %6d",t1);
    printf("   total\n");
    for(t1=1;t1<mcmc.nTemps;t1++) {
      printf("            %3d",t1);
      for(t2=0;t2<mcmc.nTemps-1;t2++) {
        if(t2<t1) {
          printf("  %6.4f",(double)mcmc.swapTss[t2][t1]/(double)mcmc.iIter);
        } else {
          printf("        ");
        }
      }
      printf("  %6.4f\n",(double)mcmc.swapTs2[t1]/(double)mcmc.iIter);
    }
    printf("          total");
    for(t1=0;t1<mcmc.nTemps-1;t1++) printf("  %6.4f",(double)mcmc.swapTs1[t1]/(double)mcmc.iIter);
    printf("\n\n");
  } //if(mcmc.prParTempInfo==2 && tempi==mcmc.nTemps-1)
} // End writeChainInfo
// ****************************************************************************************************************************************************  







// ****************************************************************************************************************************************************  
/**
 * \brief Choose and print offset starting values for the Markov chain
 * 
 * Set MCMC parameters to either the best-guess values or the injection values (where possible).
 * Then, depending on the detailed per-parameter settings, add a random offset from this value.
 * Require a good value for logL (determined by mcmc->minlogL) in order to accept the proposed starting values, unless MCMC parameters do not fully match the injection parameters.
 * This happens independently of whether parameters are fixed or not!
 * Finally, print the selected starting values to screen.
 */
// ****************************************************************************************************************************************************  
void startMCMCOffset(struct parSet *par, struct MCMCvariables *mcmc, struct interferometer *ifo[], struct runPar run)
{
  int i=0, iInj=0, nStart=0, nDiffPar=0;
  double db = 0.0;
  
  printf("\n");
  mcmc->logL[mcmc->iTemp] = -9999.999;
  
  
  // *** Set each MCMC parameter to either the best-guess value or the injection value, depending on the per-parameter settings - nothing random about this bit
  // At the start, nParam[][] contains the injection values
  if(mcmc->injectSignal >= 1) {
    for(i=0;i<mcmc->nMCMCpar;i++) {
      
      //Start at or around BestValue:
      if(mcmc->parStartMCMC[i]==1 || mcmc->parStartMCMC[i]==2) mcmc->param[mcmc->iTemp][i] = mcmc->parBestVal[i];
      
      //Start at or around the injection value where possible:
      if(mcmc->offsetMCMC == 0 || mcmc->parStartMCMC[i]==3 || mcmc->parStartMCMC[i]==4) {
        iInj = mcmc->injRevID[mcmc->parID[i]];  //Get the index of this parameter in the injection set.  -1 if not available.
	
        if(mcmc->injParUse[mcmc->parID[i]] == 1) {
          mcmc->param[mcmc->iTemp][i] = mcmc->injParVal[iInj];  //Start at or around the injection value
        } else {
          if(mcmc->parID[i]==21 && mcmc->injID[i]==22) {
            mcmc->param[mcmc->iTemp][i] = exp(3.0*mcmc->nParam[mcmc->iTemp][i]);  // Injection uses log(d), MCMC uses d^3
            if(mcmc->beVerbose>=1) printf("   I translated  log(d_L/Mpc) = %lf  to  d_L^3 = %lf Mpc^3\n",mcmc->nParam[mcmc->iTemp][i],mcmc->param[mcmc->iTemp][i]);
          } else if(mcmc->parID[i]==22 && mcmc->injID[i]==21) {
            mcmc->param[mcmc->iTemp][i] = log(mcmc->nParam[mcmc->iTemp][i])/3.0;  // Injection uses d^3, MCMC uses log(d)
            if(mcmc->beVerbose>=1) printf("   I translated  d_L^3 = %lf  to  log(d_L/Mpc) = %lf Mpc^3\n",mcmc->nParam[mcmc->iTemp][i],mcmc->param[mcmc->iTemp][i]);
          } else {
            mcmc->param[mcmc->iTemp][i] = mcmc->parBestVal[i];        // Set the MCMC parameter to BestValue - this should only happen if the injection waveform has different parameters than the MCMC waveform
            nDiffPar += 1;
          }
        }
      }
      mcmc->nParam[mcmc->iTemp][i] = mcmc->param[mcmc->iTemp][i];
    } // for i...
    
    
    // Safety check:
    if(mcmc->mcmcWaveform == mcmc->injectionWaveform && nDiffPar != 0) {
      if(nDiffPar==1) {
        fprintf(stderr, "\n ***  Warning:  The injection and MCMC waveform are identical, but 1 parameter was found to be different ***\n\n");
      } else {
        fprintf(stderr, "\n ***  Warning:  The injection and MCMC waveform are identical, but %i parameters were found to be different ***\n\n",nDiffPar);
      }
    }
    
  } else {  // if(mcmc->injectSignal <= 0), i.e., no signal injected; always use bestValue
    for(i=0;i<mcmc->nMCMCpar;i++) {
      mcmc->nParam[mcmc->iTemp][i] = mcmc->parBestVal[i];
      mcmc->param[mcmc->iTemp][i] = mcmc->nParam[mcmc->iTemp][i];
    }
  }
  
  
  // Now, param[][] == nParam[][]
  
  
  // *** Add a random offset to the MCMC starting parameters:
  if(mcmc->offsetMCMC != 0) {
    //while(mcmc->logL[mcmc->iTemp] < mcmc->minlogL+1.0) { // Accept only good starting values *** don't do this for starting values - only for updates later on ***
    while(mcmc->logL[mcmc->iTemp] < 0.1) {   // Accept only good starting values
      mcmc->acceptPrior[mcmc->iTemp] = 1;
      
      for(i=0;i<mcmc->nMCMCpar;i++) {  //For each MCMC parameter
        if(mcmc->parStartMCMC[i]==2 || mcmc->parStartMCMC[i]==4 || mcmc->parStartMCMC[i]==5) {  //Then find random offset parameters
	  
          if(mcmc->parStartMCMC[i]==2 || mcmc->parStartMCMC[i]==4) {
            mcmc->param[mcmc->iTemp][i] = mcmc->nParam[mcmc->iTemp][i] + mcmc->offsetX * gsl_ran_gaussian(mcmc->ran, mcmc->parSigma[i]);  //Gaussian with width parSigma around either Injection or BestValue
          } else if(mcmc->parStartMCMC[i]==5) {
            db = mcmc->priorBoundUp[i]-mcmc->priorBoundLow[i];                                     // Width of range
            mcmc->param[mcmc->iTemp][i] = mcmc->priorBoundLow[i] + gsl_rng_uniform(mcmc->ran)*db;        // Draw random number uniform on range with width db
          }
          mcmc->acceptPrior[mcmc->iTemp] *= (int)prior(&mcmc->param[mcmc->iTemp][i],i,*mcmc);
	  
        } // if(mcmc->parStartMCMC[i]==2 || mcmc->parStartMCMC[i]==4 || mcmc->parStartMCMC[i]==5) {  //Then find random offset parameters
      } //i
      
      if(mcmc->acceptPrior[mcmc->iTemp]==1) {                     //Check the value of the likelihood for this draw
        arr2par(mcmc->param, par, *mcmc);                             //Get the parameters from their array
        int injectionWF = 0;                                                 // Call netLogLikelihood with an MCMC waveform
        localPar(par, ifo, mcmc->networkSize, injectionWF, run);
        mcmc->logL[mcmc->iTemp] = netLogLikelihood(par, mcmc->networkSize, ifo, mcmc->mcmcWaveform, injectionWF, run);  //Calculate the likelihood
      }
      nStart = nStart + 1;
      
      
      // Print trial starting values to screen:
      if(mcmc->beVerbose>=1 && (nStart % mcmc->thinScreenOutput)==0) {
        printf("%9d%10.3lf",nStart,max(mcmc->logL[mcmc->iTemp],-99999.999));
        for(i=0;i<mcmc->nMCMCpar;i++) {
          if(mcmc->parID[i]>=11 && mcmc->parID[i]<=19) {  //GPS time
            printf(" %18.4f",mcmc->param[mcmc->iTemp][i]);
          } else {
            printf(" %9.4f",mcmc->param[mcmc->iTemp][i]);
          }
        }
        printf("\n");
      }
      
      if(mcmc->mcmcWaveform != mcmc->injectionWaveform && nDiffPar != 0 && nStart > 1e4) break;  //Don't require good logL if not all parameters match between injection and MCMC waveforms  *** This gives starting values much farther from the signal! *** -> do it after 10^4 trials
      
    }  //while(mcmc->logL[mcmc->iTemp] < mcmc->minlogL+1.0) // Accept only good starting values
  } //if(mcmc->offsetMCMC != 0)
  
  
  
  // *** Print selected starting parameters to screen:
  //Print parameter names:
  printf("\n%9s%10s", "nDraws","logL");
  for(i=0;i<mcmc->nMCMCpar;i++) {
    if(mcmc->parID[i]>=11 && mcmc->parID[i]<=19) {  //GPS time
      printf(" %18s",mcmc->parAbrev[mcmc->parID[i]]);
    } else {
      printf(" %9s",mcmc->parAbrev[mcmc->parID[i]]);
    }
  }
  printf("\n");
  
  //Print parameter values:
  printf("%9d%10.3lf",nStart,mcmc->logL[mcmc->iTemp]);
  for(i=0;i<mcmc->nMCMCpar;i++) {
    if(mcmc->parID[i]>=11 && mcmc->parID[i]<=19) {  //GPS time
      printf(" %18.4f",mcmc->param[mcmc->iTemp][i]);
    } else {
      printf(" %9.4f",mcmc->param[mcmc->iTemp][i]);
    }
  }
  printf("\n");
  
} // End void startMCMCOffset()
// ****************************************************************************************************************************************************  





// ****************************************************************************************************************************************************  
/**
 * \brief Set up the temperature ladder for parallel tempering
 */
// ****************************************************************************************************************************************************  
void setTemperatureLadder(struct MCMCvariables *mcmc)
{
  int tempi=0;
  double tempratio = 0.0;
  
  if(mcmc->nTemps < 3) { 
    exp(log(mcmc->maxTemp)/(double)(mcmc->nTemps-1));
  } else {
    tempratio = exp(log(mcmc->maxTemp)/(double)(mcmc->nTemps-2));
  }
  
  if(mcmc->prParTempInfo>0 && mcmc->beVerbose>=1) {
    printf("   Temperature ladder:\n     Number of chains:%3d,  Tmax:%7.2lf, Ti/Ti-1:%7.3lf, Overlap:%5.2lf\n",mcmc->nTemps,mcmc->maxTemp,tempratio,mcmc->tempOverlap);
    if(mcmc->parallelTempering==1) printf("     Using fixed temperatures for the chains\n");
    if(mcmc->parallelTempering==2) printf("     Using sinusoid temperatures for the chains\n");
    if(mcmc->parallelTempering==3) printf("     Using a manual temperature ladder with fixed temperatures for the chains\n");
    if(mcmc->parallelTempering==4) printf("     Using a manual temperature ladder with sinusoid temperatures for the chains\n");
    printf("     Chain     To     Ampl.    Tmin     Tmax\n");
  }
  
  for(tempi=0;tempi<mcmc->nTemps;tempi++) {
    // Set (zero-amplitude) temperatures:
    if(mcmc->nTemps < 3) {
      mcmc->tempLadder[tempi] = exp(log(mcmc->maxTemp)/(double)(mcmc->nTemps-1)*(double)tempi);
    } else {
      
      if(tempi == 0) {
        mcmc->tempLadder[tempi] = 1.0;
      } else if(tempi == 1) {
        mcmc->tempLadder[tempi+1] = exp(log(mcmc->maxTemp)/(double)(mcmc->nTemps-2)*(double)tempi);
        mcmc->tempLadder[tempi] = exp(log((mcmc->tempLadder[tempi-1]) + log(mcmc->tempLadder[tempi+1]))/2.0);
      } else if(tempi < mcmc->nTemps-1) {
        mcmc->tempLadder[tempi+1] = exp(log(mcmc->maxTemp)/(double)(mcmc->nTemps-2)*(double)tempi);
      }
      
    }
    if(mcmc->parallelTempering==3 || mcmc->parallelTempering==4) mcmc->tempLadder[tempi] = mcmc->tempLadder[tempi];  //Set manual ladder
  } //  for(tempi=0;tempi<mcmc->nTemps;tempi++) {
  
  for(tempi=mcmc->nTemps-1;tempi>=0;tempi--) {
    if(tempi==0 || mcmc->parallelTempering<=1 || mcmc->parallelTempering==3) {
      mcmc->tempAmpl[tempi] = 0.0;
    } else {
      if(mcmc->nTemps < 3) {
        mcmc->tempAmpl[tempi] = mcmc->tempOverlap*(mcmc->tempLadder[tempi] - mcmc->tempLadder[tempi-1])/(tempratio+1.0)*tempratio;  //Temperatures of adjacent chains overlap a lot at extrema (since in antiphase), make sure Ti,min>=T-i1,0
      } else {
	
        if(tempi == 1) {
          if(mcmc->tempOverlap < 1.01) {
            mcmc->tempLadder[tempi] = (mcmc->tempLadder[tempi-1] + (mcmc->tempLadder[tempi+1] - mcmc->tempAmpl[tempi+1]))/2.0;
            mcmc->tempAmpl[tempi] = fabs(mcmc->tempLadder[tempi-1] - (mcmc->tempLadder[tempi+1] - mcmc->tempAmpl[tempi+1]))/2.0;  //Temperatures of adjacent chains overlap a lot at extrema (since in antiphase), make sure Ti,min>=T-i1,0
          } else {
            //mcmc->tempLadder[tempi] = (mcmc->tempLadder[tempi-1] + (mcmc->tempLadder[tempi+1] - mcmc->tempAmpl[tempi+1]/mcmc->tempOverlap))/2.0;
            //mcmc->tempAmpl[tempi] = mcmc->tempOverlap * fabs(mcmc->tempLadder[tempi-1] - (mcmc->tempLadder[tempi+1] - mcmc->tempAmpl[tempi+1]/mcmc->tempOverlap))/2.0;  //Temperatures of adjacent chains overlap a lot at extrema (since in antiphase), make sure Ti,min>=T-i1,0
            mcmc->tempAmpl[tempi] = mcmc->tempAmpl[tempi+1]/tempratio;
            mcmc->tempLadder[tempi] = exp((log(mcmc->tempLadder[tempi-1]) + log(mcmc->tempLadder[tempi+1] - mcmc->tempAmpl[tempi+1]/mcmc->tempOverlap))/2.0);
          }
        } else if(tempi == 2) {
          mcmc->tempAmpl[tempi] = mcmc->tempOverlap*(mcmc->tempLadder[tempi] - mcmc->tempLadder[tempi-2])/(tempratio+1.0)*tempratio;  //Temperatures of adjacent chains overlap a lot at extrema (since in antiphase), make sure Ti,min>=T-i1,0
        } else {
          mcmc->tempAmpl[tempi] = mcmc->tempOverlap*(mcmc->tempLadder[tempi] - mcmc->tempLadder[tempi-1])/(tempratio+1.0)*tempratio;  //Temperatures of adjacent chains overlap a lot at extrema (since in antiphase), make sure Ti,min>=T-i1,0
        }
	
      }
      
    }
  }
  
  if(mcmc->prParTempInfo>0 && mcmc->beVerbose>=1) {
    for(tempi=0;tempi<mcmc->nTemps;tempi++) {
      if(mcmc->tempLadder[tempi]-mcmc->tempAmpl[tempi] < 1.0) {
        printf("     %3d  %7.2lf  %7.2lf  %7.2lf* %7.2lf    * I will use max( T, 1.0 )",tempi,mcmc->tempLadder[tempi],mcmc->tempAmpl[tempi],mcmc->tempLadder[tempi]-mcmc->tempAmpl[tempi],mcmc->tempLadder[tempi]+mcmc->tempAmpl[tempi]);
      } else {
        printf("     %3d  %7.2lf  %7.2lf  %7.2lf  %7.2lf",tempi,mcmc->tempLadder[tempi],mcmc->tempAmpl[tempi],mcmc->tempLadder[tempi]-mcmc->tempAmpl[tempi],mcmc->tempLadder[tempi]+mcmc->tempAmpl[tempi]);
      }
      printf("\n");
    }
    printf("\n\n");
  }
} // End setTemperatureLadder()
// ****************************************************************************************************************************************************  





// ****************************************************************************************************************************************************  
/**
 * \brief Set up the temperature ladder for parallel tempering - old routine
 */
// ****************************************************************************************************************************************************  
void setTemperatureLadderOld(struct MCMCvariables *mcmc)
{
  int tempi=0;
  double tempratio = exp(log(mcmc->maxTemp)/(double)(mcmc->nTemps-1));
  if(mcmc->prParTempInfo>0 && mcmc->beVerbose>=1) {
    printf("   Temperature ladder:  (old routine)\n     Number of chains:%3d,  Tmax:%7.2lf, Ti/Ti-1:%7.3lf\n",mcmc->nTemps,mcmc->maxTemp,tempratio);
    if(mcmc->parallelTempering==1) printf("     Using fixed temperatures for the chains\n");
    if(mcmc->parallelTempering==2) printf("     Using sinusoid temperatures for the chains\n");
    if(mcmc->parallelTempering==3) printf("     Using a manual temperature ladder with fixed temperatures for the chains\n");
    if(mcmc->parallelTempering==4) printf("     Using a manual temperature ladder with sinusoid temperatures for the chains\n");
    printf("     Chain     To     Ampl.    Tmin     Tmax\n");
  }
  for(tempi=0;tempi<mcmc->nTemps;tempi++) {
    mcmc->tempLadder[tempi] = pow(10.0,log10(mcmc->maxTemp)/(double)(mcmc->nTemps-1)*(double)tempi);
    if(mcmc->parallelTempering==3 || mcmc->parallelTempering==4) mcmc->tempLadder[tempi] = mcmc->tempLadder[tempi];  //Set manual ladder
    
    if(tempi==0 || mcmc->parallelTempering<=1 || mcmc->parallelTempering==3) {
      mcmc->tempAmpl[tempi] = 0.0;
    } else {
      //mcmc->tempAmpl[tempi] = (mcmc->tempLadder[tempi] - mcmc->tempLadder[tempi-1])/(tempratio+1.0)*tempratio;  //Temperatures of adjacent chains just touch at extrema (since in antiphase)
      //mcmc->tempAmpl[tempi] = 1.5*(mcmc->tempLadder[tempi] - mcmc->tempLadder[tempi-1])/(tempratio+1.0)*tempratio;  //Temperatures of adjacent chains overlap somewhat at extrema (since in antiphase)
      mcmc->tempAmpl[tempi] = min(3.0*(mcmc->tempLadder[tempi] - mcmc->tempLadder[tempi-1])/(tempratio+1.0)*tempratio , fabs(mcmc->tempLadder[tempi]-mcmc->tempLadder[tempi-1]));  //Temperatures of adjacent chains overlap a lot at extrema (since in antiphase), make sure Ti,min>=T-i1,0
      //if(mcmc->nTemps>10)  mcmc->tempAmpl[tempi] = min(3.0*(mcmc->tempLadder[tempi] - mcmc->tempLadder[tempi-1])/(tempratio+1.0)*tempratio , fabs(mcmc->tempLadder[tempi]-mcmc->tempLadder[tempi-2]));  //Temperatures of adjacent chains overlap a lot at extrema (since in antiphase), make sure Ti,min>=T-i1,0
      //mcmc->tempAmpl[tempi] = fabs(mcmc->tempLadder[tempi]-mcmc->tempLadder[tempi-1]);  //Temperatures of adjacent chains overlap: Amplitude = (T(i) - T(i-1))  (may be a bit smallish for large nTemps)
    }
    if(mcmc->prParTempInfo>0 && mcmc->beVerbose>=1) printf("     %3d  %7.2lf  %7.2lf  %7.2lf  %7.2lf\n",tempi,mcmc->tempLadder[tempi],mcmc->tempAmpl[tempi],mcmc->tempLadder[tempi]-mcmc->tempAmpl[tempi],mcmc->tempLadder[tempi]+mcmc->tempAmpl[tempi]);
  }
  if(mcmc->prParTempInfo>0 && mcmc->beVerbose>=1) printf("\n\n");
} // End of setTemperatureLadderOld
// ****************************************************************************************************************************************************  

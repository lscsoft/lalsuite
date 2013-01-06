/* 
   
   SPINspiral:                parameter estimation on binary inspirals detected by LIGO, including spins of the binary members
   SPINspiral_main.c:         main routine
   
   
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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <SPINspiral.h>

double Ms,Mpc,G,c,Mpcs,pi,tpi,mtpi;

/**
 * \file
 * \brief SPINspiral documentation
 * SPINspiral is a parameter-estimation code designed to extract the physical parameters of compact-binary coalescences (CBCs) from observed gravitational-wave signals.
 *
 * More information about SPINspiral can be found on its web page: 
 * <a href="http://www.astro.northwestern.edu/~sluys/index.php?title=SPINspiral">http://www.astro.northwestern.edu/~sluys/index.php?title=SPINspiral</a>.
 *
 * The pages in this documentation provide information on the SPINspiral code.
 *
 */


// Main program:
int main(int argc, char* argv[])
{
  printf("\n\n   Starting SPINspiral...\n");
  printf("   Compiled from source code version $Id$ \n");
  
  clock_t time0 = clock();
  int ifonr=0,i=0,injectionWF=0,mcmcWF=0;
  double snr=0.0;
  
  
  //Initialise stuff for the run:
  struct runPar run;
  run.maxnPar = 20;                        //The maximum number of allowed MCMC/injection parameters (this number is hardcoded in many places in SPINspiral.h)
  run.parDBn = 200;                        //The size of the hardcoded parameter database (this number is hardcoded in many places in SPINspiral.h)
  for(i=0;i<run.parDBn;i++) {
    run.parDef[i]     = 0;
    run.parRevID[i]   = -1;
    run.injRevID[i]   = -1;
    run.mcmcParUse[i] = 0;
    run.injParUse[i]  = 0;
  }
  sprintf(run.executable,"%s",argv[0]);
  run.lowFrequencyCut = 0.0;
  run.injXMLfilename = NULL;
  run.injXMLnr = -1;
  for(i=0;i<99;i++) run.commandSettingsFlag[i] = 0;
  setConstants();                          //Set the global constants (which are variable in C)
  setParameterNames(&run);                 //Set the names of the parameters in the hardcoded parameter database
  
  sprintf(run.mainFilename,"SPINspiral.input");  //Default input filename
  readCommandLineOptions(argc,argv,&run);  //Read the command-line options
  
  
  readMainInputfile(&run);                 //Read main input data file for this run from input.mcmc
  readMCMCinputfile(&run);                 //Read the input data on how to do MCMC 
  if(run.MCMCseed==0) {
    setSeed(&run.MCMCseed);                  //Set MCMCseed if 0, otherwise keep the current value
    if(run.beVerbose>=1) printf("   Picking seed from the system clock to start Markov chains from randomly offset values: %d\n", run.MCMCseed);
  }
  readInjectionInputfile(&run);            //Read the input data on whether and how to do a software injection
  if(run.injXMLfilename != NULL && run.injXMLnr >= 0) readInjectionXML(&run);    //Read injection XML file if specified:
  readParameterInputfile(&run);            //Read the input data on how to handle MCMC parameters
  readSystemInputfile(&run);               //Read system-dependent data, e.g. path to data files
  
  
  
  //Set up the data for the IFOs in an IFO database you may want to use (H1,L1 + VIRGO by default)
  run.maxIFOdbaseSize = 3;  //The maximum number of IFOs to read the properties in for from the data input file (SPINspiral.input.data or equivalent)
  struct interferometer database[run.maxIFOdbaseSize];
  setIFOdata(&run, database);
  
  
  
  //Define interferometer network with IFOs.  The first run.networkSize are actually used
  //struct interferometer *network[3] = {&database[0], &database[1], &database[2]}; //H1L1V
  struct interferometer *network[run.networkSize];
  for(ifonr=0;ifonr<run.networkSize;ifonr++) network[ifonr] = &database[run.selectifos[ifonr]-1];
  int networkSize = run.networkSize;
  
  
  
  //Initialise interferometers, read and prepare data, inject signal (takes some time)
  if(networkSize == 1) {
    printf("   Initialising 1 IFO: %s, reading noise and data...\n",database[run.selectifos[0]-1].name);
  } else {
    printf("   Initialising %d IFOs: ",networkSize);
    for(ifonr=0;ifonr<run.networkSize;ifonr++) printf(" %s,",database[run.selectifos[ifonr]-1].name);
    printf(" reading noise and data files...\n");
  }
  IFOinit(network, networkSize, run); //Do the actual initialisation
  
  
  // Get a injection parameter set to calculate SNR or write the wavefrom to disc:
  struct parSet injParSet;
  if(run.injectSignal >= 1) {
    if(run.injectionSNR < 0.001) printf("   A signal with the injection parameter values was injected into the data.\n");
    
    getInjectionParameters(&injParSet, run.nInjectPar, run.injParVal);
    allocParset(&injParSet, networkSize);
    injectionWF = 1;
    localPar(&injParSet, network, networkSize, injectionWF, run);
    
    
    // Calculate SNR:
    run.netsnr = 0.0;
    if(run.doSNR==1) {
      for(ifonr=0; ifonr<networkSize; ++ifonr) {
        injectionWF = 1;                           //Call signalToNoiseRatio with the injection waveform
        snr = signalToNoiseRatio(&injParSet, network, ifonr, run.injectionWaveform, injectionWF, run);
        network[ifonr]->snr = snr;
        run.netsnr += snr*snr;
      }
      run.netsnr = sqrt(run.netsnr);
    }
    
    
    // Get the desired SNR by scaling the distance:
    if(run.injectionSNR > 0.001 && run.injectSignal>=1) {
      run.injParVal[3] += log(run.netsnr/run.injectionSNR);  //Use total network SNR
      printf("   Setting distance to %lf Mpc (log(d/Mpc)=%lf) to get a network SNR of %lf.\n",exp(run.injParVal[3]),run.injParVal[3],run.injectionSNR);
      freeParset(&injParSet);
      getInjectionParameters(&injParSet, run.nMCMCpar, run.injParVal);
      allocParset(&injParSet, networkSize);
      injectionWF = 1;
      localPar(&injParSet, network, networkSize, injectionWF, run);
      
      // Recalculate SNR:
      run.netsnr = 0.0;
      if(run.doSNR==1) {
        for(ifonr=0; ifonr<networkSize; ++ifonr) {
          injectionWF = 1;                           //Call signalToNoiseRatio with the injection waveform
          snr = signalToNoiseRatio(&injParSet, network, ifonr, run.injectionWaveform, injectionWF, run);
          network[ifonr]->snr = snr;
          run.netsnr += snr*snr;
        }
        run.netsnr = sqrt(run.netsnr);
      }
      
      for(ifonr=0; ifonr<networkSize; ++ifonr) IFOdispose(network[ifonr], run);
      //Reinitialise interferometers, read and prepare data, inject signal (takes some time)
      if(networkSize == 1) {
        printf("   Reinitialising 1 IFO, reading data...\n");
      } else {
        printf("   Reinitialising %d IFOs, reading datafiles...\n",networkSize);
      }
      IFOinit(network, networkSize, run);
      printf("   A signal with the 'true' parameter values was injected.\n");
    }
    
  } else { //if(run.injectSignal <= 0)
    for(ifonr=0; ifonr<networkSize; ++ifonr) network[ifonr]->snr = 0.0;
    run.netsnr = 0.0;
    printf("   No signal was injected.\n");
  }
  
  
  
  
  //Write some data parameters to screen:
  if(run.beVerbose >= 1) {
    printf("\n\n");
    printf("%10s  %11s  %10s  %10s  %9s  %17s  %14s  %12s  %12s  %12s  %12s\n",
           "Detector","f_low","f_high","before tc","after tc","Sample start (GPS)","Sample length","Downsample","Sample rate","Sample size","FT size");
    for(ifonr=0;ifonr<run.networkSize;ifonr++) {
      printf("%10s  %8.2lf Hz  %7.2lf Hz  %8.2lf s  %7.2lf s  %16.5lf s  %12.4lf s  %10d x  %9d Hz  %9d pt  %9d pt\n",
             network[ifonr]->name,network[ifonr]->lowCut,network[ifonr]->highCut,network[ifonr]->before_tc,network[ifonr]->after_tc,
             network[ifonr]->FTstart,network[ifonr]->deltaFT,run.downsampleFactor,network[ifonr]->samplerate,network[ifonr]->samplesize,network[ifonr]->FTsize);
    }
  }
  
  
  //Write the data and its FFT, the signal and its FFT, and the noise ASD to disc
  if(run.writeSignal) {
    writeDataToFiles(network, networkSize, run);
    writeNoiseToFiles(network, networkSize, run);
    writeSignalsToFiles(network, networkSize, run);
  }  
  
  //Write some injection parameters to screen:  CHECK: needs fix
  //printf("\n");
  //printf("   Global     :    Source position:  RA: %5.2lfh, dec: %6.2lfd;  J0 points to:  RA: %5.2lfh, dec: %6.2lfd;   inclination J0: %5.2lfd  \n",  injParSet.par[run.parRevID[31]]*r2h, asin(injParSet.par[run.parRevID[32]])*r2d,  injParSet.par[run.parRevID[54]]*r2h,  asin(injParSet.par[run.parRevID[53]])*r2d, (pi/2.0-acos(injParSet.NdJ))*r2d );
  //printf("   Global     :    Source position:  RA: %5.2lfh, dec: %6.2lfd\n",  injParSet.par[run.parRevID[31]]*r2h, asin(injParSet.par[run.parRevID[32]])*r2d);
  //for(ifonr=0;ifonr<networkSize;ifonr++) printf("   %-11s:    theta: %5.1lfd,  phi: %5.1lfd;   azimuth: %5.1lfd,  altitude: %5.1lfd\n",network[ifonr]->name,injParSet.localti[ifonr]*r2d,injParSet.locazi[ifonr]*r2d,fmod(pi-(injParSet.locazi[ifonr]+network[ifonr]->rightArm)+mtpi,tpi)*r2d,(pi/2.0-injParSet.localti[ifonr])*r2d);
  
  if(run.beVerbose >= 1) {
    printf("\n  %10s  %10s  %6s  %6s  ","nIter","nBurn","seed","nDet");
    for(ifonr=0;ifonr<networkSize;ifonr++) printf("%18s%4s     ",network[ifonr]->name,"SNR");
    printf("%18s  ","Network SNR");
    
    double maxSNR=0.0, relSNR[99];
    for(ifonr=0;ifonr<networkSize;ifonr++) maxSNR = max(maxSNR,network[ifonr]->snr);
    for(ifonr=0;ifonr<networkSize;ifonr++) relSNR[ifonr] = network[ifonr]->snr/maxSNR;
    
    printf("\n  %10d  %10d  %6d  %6d  ",run.nIter,run.annealNburn,run.MCMCseed,networkSize);
    for(ifonr=0;ifonr<networkSize;ifonr++) printf("%18.8lf (%4.2lf)  ",network[ifonr]->snr,relSNR[ifonr]);
    printf("%18.8lf\n\n",run.netsnr);
  }
  
  //Print actual injection parameters to screen:
  if(run.injectSignal >= 1 && run.doMCMC==0 && run.beVerbose>=1) {
    printf("   Injection param:");
    for(i=0;i<run.nMCMCpar;i++) {
      if(run.parID[i]>=11 && run.parID[i]<=19) {  //GPS time
        printf(" %18s",run.parAbrev[run.parID[i]]);
      } else {
        printf(" %9s",run.parAbrev[run.parID[i]]);
      }
    }
    printf("\n");
    
    printf("                   ");
    for(i=0;i<run.nMCMCpar;i++) {
      if(run.parID[i]>=11 && run.parID[i]<=19) {  //GPS time
        printf(" %18.4f",injParSet.par[i]);
      } else {
        printf(" %9.4f",injParSet.par[i]);
      }
    }
    printf("\n\n");
  }
  
  
  //Do MCMC
  clock_t time1 = clock();
  if(run.doMCMC==1) {
    MCMC(run, network);
  }
  clock_t time2 = clock();
  
  
  
  
  
  //Calculate matches between two signals
  if(run.doMatch==1) {
    /*
      if(1==2) {
      printf("\n\n");
      getInjectionParameters(&injParSet, run.nMCMCpar, run.injParVal);
      allocParset(&injParSet, networkSize);
      
      FILE *fout;
      fout = fopen("tc.dat","w");
      double fac=0.0;
      double matchpar = injParSet.tc,matchres=0.0;
      for(fac=-0.002;fac<0.002;fac+=0.00005) {
      injParSet.tc = matchpar+fac;
      for(ifonr=0;ifonr<networkSize;ifonr++) {
      localPar(&injParSet, network, networkSize, injectionWF, run);
      matchres = match(&injParSet,network,ifonr,networkSize);
      printf("%10.6f  %10.6f\n",fac,matchres);
      fprintf(fout,"%10.6f  %10.6f\n",fac,matchres);
      }
      }
      fclose(fout);
      }
    */
    
    
    //Compute match between waveforms with parameter sets injctPar and startPar
    if(1==1) {
      printf("\n\n");
      
      injectionWF = 1;  // When calling localPar, parMatch and parOverlap with the injection waveform (1)
      mcmcWF = 0;       // When calling localPar, parMatch and parOverlap with the MCMC waveform (0)
      
      struct parSet injctPar;
      allocParset(&injctPar, networkSize);
      
      struct parSet startPar;
      allocParset(&startPar, networkSize);
      
      getInjectionParameters(&injctPar, run.nInjectPar, run.injParVal);  // injctPar contains injection parameters
      //getStartParameters(&injctPar, run);                                // injctPar contains MCMC starting parameters
      //getInjectionParameters(&startPar, run.nInjectPar, run.injParVal);  // startPar contains injection parameters
      getStartParameters(&startPar, run);                                // startPar contains MCMC starting parameters
      
      localPar(&injctPar, network, networkSize, injectionWF, run);
      localPar(&startPar, network, networkSize, mcmcWF, run);
      
      //run.injectionWaveform = 1;
      double matchres=0.0, overlap=0.0;
      printf("\n\n  Match:\n");
      matchres = parMatch(&injctPar, run.injectionWaveform, injectionWF, &startPar, run.mcmcWaveform, mcmcWF, network, networkSize, run);
      printf("\n\n  Overlap:\n");
      overlap = parOverlap(&injctPar, run.injectionWaveform, injectionWF, &startPar, run.mcmcWaveform, mcmcWF, network, 0, run); //For IFO 0
      
      printf("\n\n   Match: %10.5lf,  overlap IFO 0: %g \n",matchres,overlap);
      printf("\n");
      
      freeParset(&injctPar);
      freeParset(&startPar);
    }
    
    
    //Compute Fisher matrix for parameter set par
    if(1==2) {
      printf("\n\n  Computing Fisher matrix...\n\n");
      i = 0;
      int j=0;
      struct parSet par;
      double **matrix  = (double**)calloc(run.nMCMCpar,sizeof(double*));
      for(i=0;i<run.nMCMCpar;i++) matrix[i]  = (double*)calloc(run.nMCMCpar,sizeof(double));
      allocParset(&par, networkSize);
      //getparameterset(&par, 3.0,0.11,700009012.346140,3.0, 0.5,0.9,3.0,0.5, 1.0,0.1,2.0,3.0);
      
      //computeFisherMatrixIFO(par,run.nMCMCpar,network,networkSize,0,matrix);
      //computeFisherMatrix(&par,run.nMCMCpar,network,networkSize,matrix);
      
      for(i=0;i<run.nMCMCpar;i++) {
        for(j=0;j<run.nMCMCpar;j++) {
          printf("  %12.4g",matrix[i][j]/(double)networkSize);
        }
        printf("\n\n");
      }
      printf("\n\n");
      
      
      
      freeParset(&par);
      for(i=0;i<run.nMCMCpar;i++) {
        free(matrix[i]);
      }
      free(matrix);
    }
  } //if(run.doMatch=1)
  
  
  
  
  
  //Get rid of allocated memory and quit
  for(ifonr=0; ifonr<networkSize; ++ifonr) IFOdispose(network[ifonr], run);
  if(run.injectSignal >= 1) freeParset(&injParSet);
  
  
  clock_t time3 = clock();
  if(run.beVerbose>=1) { 
    printf("   Timing:\n");
    if(run.doMCMC>=1) {
      printf("     initialisation:%10.2lfs\n", ((double)time1 - (double)time0)*1.e-6 );
      printf("     MCMC:          %10.2lfs\n", ((double)time2 - (double)time1)*1.e-6 );
    }
    printf("     total time:    %10.2lfs\n", (double)time3*1.e-6);
  }
  
  printf("\n   SPINspiral done.\n\n");
  if(run.doMCMC>=1) printf("\n");
  return 0;
}



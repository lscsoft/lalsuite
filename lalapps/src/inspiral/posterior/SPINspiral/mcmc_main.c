/* 
   
   SPINspiral:                parameter estimation on binary inspirals detected by LIGO, including spins of the binary members
   mcmc_main.c:               main routine
   
   
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


#include <mcmc.h>


// Main program:
int main(int argc, char * argv[])
{
  if(doMCMC>=1) printf("\n");
  printf("\n   Starting SPINspiral...\n");
  printf("   Produced with source code version $Id: mcmc_main.c 149 2009-05-07 23:55:57Z sluys $ \n");
  
  clock_t time0 = clock();
  int ifonr=0,i=0;
  double snr=0.0;
  
  
  
  //Initialise stuff for the run
  struct runPar run;
  sprintf(run.executable,argv[0]);
  /*
  if(system( NULL )) {  //Is incompatible with condor_compile
    char shellCommand[199];
    sprintf(shellCommand,"echo -n '   Executable:  %s,  compiled:  ';  ls -l --time-style=\"+%%a %%e %%b %%Y, %%T %%Z (UTC%%z)\" %s | gawk '{print $6,$7,$8,$9,$10,$11,$12}'",run.executable,run.executable);
    system(shellCommand);
    
    sprintf(shellCommand,"echo '   Run path:    '`uname -n`':'`pwd`");
    system(shellCommand);
    printf("\n");
  }
  */
  
  run.maxnPar = 20;                      //The maximum number of allowed MCMC/injection parameters (this number is hardcoded in many places in mcmc.h)
  setconstants();                        //Set the global constants (which are variable in C)
  setParameterNames(&run);               //Set the names of the parameters in the hardcoded parameter database
  sprintf(run.mainFilename,"mcmc.input");  //Default input filename
  if(argc > 1) sprintf(run.mainFilename,argv[1]);
  
  readLocalInputfile();                  //Read system-dependent data, e.g. path to data files
  readMainInputfile(&run);               //Read main input data file for this run from input.mcmc
  readMCMCinputfile(&run);               //Read the input data on how to do MCMC 
  setseed(&run.MCMCseed);                //Set MCMCseed if 0, otherwise keep the current value
  readInjectionInputfile(&run);          //Read the input data on whether and how to do a software injection
  //setRandomInjectionParameters(&run);    //Randomise the injection parameters where wanted (do this in readInjectionInputfile
  readParameterInputfile(&run);          //Read the input data on how to handle MCMC parameters
  //writeInputfile(&run);                //Write run data to nicely formatted input.mcmc.<MCMCseed>
  
  
  
  //Set up the data for the IFOs in an IFO database you may want to use (H1,L1 + VIRGO by default)
  run.maxIFOdbaseSize = 4;  //The maximum number of IFOs to read the properties in for from the data input file (mcmc.data or equivalent)
  struct interferometer database[run.maxIFOdbaseSize];
  setIFOdata(&run, database);
  
  
  
  //Define interferometer network with IFOs.  The first run.networksize are actually used
  //struct interferometer *network[3] = {&database[0], &database[1], &database[2]}; //H1L1V
  struct interferometer *network[run.networksize];
  for(ifonr=0;ifonr<run.networksize;ifonr++) network[ifonr] = &database[run.selectifos[ifonr]-1];
  int networksize = run.networksize;
  
  
  
  //Initialise interferometers, read and prepare data, inject signal (takes some time)
  if(networksize == 1) {
    printf("   Initialising 1 IFO: %s, reading noise and data...\n",database[run.selectifos[0]-1].name);
  } else {
    printf("   Initialising %d IFOs: ",networksize);
    for(ifonr=0;ifonr<run.networksize;ifonr++) printf(" %s,",database[run.selectifos[ifonr]-1].name);
    printf(" reading noise and data files...\n");
  }
  ifoinit(network, networksize, run); //Do the actual initialisation
  if(run.injectSignal>=1) {
    if(run.injectionSNR < 0.001) printf("   A signal with the injection parameter values was injected into the data.\n");
  } else {
    printf("   No signal was injected.\n");
  }
  
  
  
  //Get a parameter set to calculate SNR or write the wavefrom to disc
  struct parset dummypar;
  getInjectionParameters(&dummypar, run.nMCMCpar, run.injParVal);
  dummypar.loctc    = (double*)calloc(networksize,sizeof(double));
  dummypar.localti  = (double*)calloc(networksize,sizeof(double));
  dummypar.locazi   = (double*)calloc(networksize,sizeof(double));
  dummypar.locpolar = (double*)calloc(networksize,sizeof(double));
  localpar(&dummypar, network, networksize);
  
  
  
  //Calculate SNR
  run.netsnr = 0.0;
  if(doSNR==1) {
    for(ifonr=0; ifonr<networksize; ++ifonr) {
      snr = signaltonoiseratio(&dummypar, network, ifonr, run.injectionWaveform);
      network[ifonr]->snr = snr;
      run.netsnr += snr*snr;
    }
    run.netsnr = sqrt(run.netsnr);
  }
  
  
  
  //Get the desired SNR by scaling the distance
  if(run.injectionSNR > 0.001 && run.injectSignal>=1) {
    run.injParVal[3] += log(run.netsnr/run.injectionSNR);  //Use total network SNR
    printf("   Setting distance to %lf Mpc (log(d/Mpc)=%lf) to get a network SNR of %lf.\n",exp(run.injParVal[3]),run.injParVal[3],run.injectionSNR);
    getInjectionParameters(&dummypar, run.nMCMCpar, run.injParVal);
    dummypar.loctc    = (double*)calloc(networksize,sizeof(double));
    dummypar.localti  = (double*)calloc(networksize,sizeof(double));
    dummypar.locazi   = (double*)calloc(networksize,sizeof(double));
    dummypar.locpolar = (double*)calloc(networksize,sizeof(double));
    localpar(&dummypar, network, networksize);
    
    //Recalculate SNR
    run.netsnr = 0.0;
    if(doSNR==1) {
      for(ifonr=0; ifonr<networksize; ++ifonr) {
	snr = signaltonoiseratio(&dummypar, network, ifonr, run.injectionWaveform);
	network[ifonr]->snr = snr;
	run.netsnr += snr*snr;
      }
      run.netsnr = sqrt(run.netsnr);
    }
    
    for(ifonr=0; ifonr<networksize; ++ifonr) ifodispose(network[ifonr]);
    //Reinitialise interferometers, read and prepare data, inject signal (takes some time)
    if(networksize == 1) {
      printf("   Reinitialising 1 IFO, reading data...\n");
    } else {
      printf("   Reinitialising %d IFOs, reading datafiles...\n",networksize);
    }
    ifoinit(network, networksize, run);
    printf("   A signal with the 'true' parameter values was injected.\n");
  }
  
  
  
  
  
  //Write some data parameters to screen:
  printf("\n\n");
  printf("%10s  %11s  %10s  %10s  %9s  %17s  %14s  %12s  %12s  %12s  %12s\n",
	 "Detector","f_low","f_high","before tc","after tc","Sample start (GPS)","Sample length","Downsample","Sample rate","Sample size","FT size");
  for(ifonr=0;ifonr<run.networksize;ifonr++) {
    printf("%10s  %8.2lf Hz  %7.2lf Hz  %8.2lf s  %7.2lf s  %16.5lf s  %12.4lf s  %10d x  %9d Hz  %9d pt  %9d pt\n",
	   network[ifonr]->name,network[ifonr]->lowCut,network[ifonr]->highCut,network[ifonr]->before_tc,network[ifonr]->after_tc,
	   network[ifonr]->FTstart,network[ifonr]->deltaFT,downsamplefactor,network[ifonr]->samplerate,network[ifonr]->samplesize,network[ifonr]->FTsize);
  }
  
  
  //Write the data and its FFT, the signal and its FFT, and the noise ASD to disc
  if(writeSignal)
  {
     writeDataToFiles(network, networksize, run.MCMCseed);
     writeNoiseToFiles(network, networksize, run.MCMCseed);
     writeSignalsToFiles(network, networksize, run);
  }  
  
  //Write some injection parameters to screen:
  printf("\n");
  //printf("   Global     :    Source position:  RA: %5.2lfh, dec: %6.2lfd;  J0 points to:  RA: %5.2lfh, dec: %6.2lfd;   inclination J0: %5.2lfd  \n",  dummypar.par[run.parRevID[31]]*r2h, asin(dummypar.par[run.parRevID[32]])*r2d,  dummypar.par[run.parRevID[54]]*r2h,  asin(dummypar.par[run.parRevID[53]])*r2d, (pi/2.0-acos(dummypar.NdJ))*r2d );
  printf("   Global     :    Source position:  RA: %5.2lfh, dec: %6.2lfd\n",  dummypar.par[run.parRevID[31]]*r2h, asin(dummypar.par[run.parRevID[32]])*r2d);
  for(ifonr=0;ifonr<networksize;ifonr++) printf("   %-11s:    theta: %5.1lfd,  phi: %5.1lfd;   azimuth: %5.1lfd,  altitude: %5.1lfd\n",network[ifonr]->name,dummypar.localti[ifonr]*r2d,dummypar.locazi[ifonr]*r2d,fmod(pi-(dummypar.locazi[ifonr]+network[ifonr]->rightarm)+mtpi,tpi)*r2d,(pi/2.0-dummypar.localti[ifonr])*r2d);
  
  printf("\n  %10s  %10s  %6s  %6s  ","niter","nburn","seed","ndet");
  for(ifonr=0;ifonr<networksize;ifonr++) printf("%16s%4s  ",network[ifonr]->name,"SNR");
  printf("%20s  ","Network SNR");
  printf("\n  %10d  %10d  %6d  %6d  ",iter,nburn,run.MCMCseed,networksize);
  for(ifonr=0;ifonr<networksize;ifonr++) printf("%20.10lf  ",network[ifonr]->snr);
  printf("%20.10lf\n\n",run.netsnr);
  
  
  //Print actual injection parameters to screen:
  if(doMCMC==0) {
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
	printf(" %18.4f",dummypar.par[i]);
      } else {
	printf(" %9.4f",dummypar.par[i]);
      }
    }
    printf("\n\n");
  }
  
  
  //Do MCMC
  clock_t time1 = clock();
  if(doMCMC==1) {
    //printMuch=1;
    mcmc(run, network);
    //printMuch=0;
  }
  clock_t time2 = clock();
  
  
  
  
  
  //Calculate matches between two signals
  if(doMatch==1) {
    /*
    if(1==2) {
      printf("\n\n");
      getInjectionParameters(&dummypar, run.nMCMCpar, run.injParVal);
      dummypar.loctc    = (double*)calloc(networksize,sizeof(double));
      dummypar.localti  = (double*)calloc(networksize,sizeof(double));
      dummypar.locazi   = (double*)calloc(networksize,sizeof(double));
      dummypar.locpolar = (double*)calloc(networksize,sizeof(double));
      
      FILE *fout;
      fout = fopen("tc.dat","w");
      double fac=0.0;
      double matchpar = dummypar.tc,matchres=0.0;
      for(fac=-0.002;fac<0.002;fac+=0.00005) {
	dummypar.tc = matchpar+fac;
	for(ifonr=0;ifonr<networksize;ifonr++) {
	  localpar(&dummypar, network, networksize);
	  matchres = match(&dummypar,network,ifonr,networksize);
	  printf("%10.6f  %10.6f\n",fac,matchres);
	  fprintf(fout,"%10.6f  %10.6f\n",fac,matchres);
	}
      }
      fclose(fout);
    }
    */
    
    //Compute match between waveforms with parameter sets par1 and par2
    if(1==2) {
      printf("\n\n");
      struct parset par1;
      allocparset(&par1, networksize);
      struct parset par2;
      allocparset(&par2, networksize);
      
      double eta=0.0;
      //for(eta=0.01;eta<0.25001;eta+=0.001) {
      //for(eta=0.1;eta<0.12001;eta+=0.001) {
      for(eta=0.111;eta<0.1111;eta+=0.001) {
	//getparameterset(&par1, 3.0,0.11,700009012.346140,3.0, 0.5,0.9,3.0,0.5, 1.0,0.1,2.0,3.0);
	
	//getparameterset(&par2, 3.0,eta,700009012.346140,3.0, 0.5,0.9,3.0,0.5, 1.0,0.1,2.0,3.0);
	
	double matchres = parmatch(&par1,&par2,network, networksize, run.injectionWaveform);
	double overlap = paroverlap(&par1,&par2,network,0, run.injectionWaveform);
	
	printf("   Eta: %6.4f,  match: %10.5lf,  overlap: %g \n",eta,matchres,overlap);
      }
      printf("\n");
      
      freeparset(&par1);
      freeparset(&par2);
    }
  
    //Compute Fisher matrix for parameter set par
    if(1==2) {
      printf("\n\n  Computing Fisher matrix...\n\n");
      i = 0;
      int j=0;
      struct parset par;
      double **matrix  = (double**)calloc(run.nMCMCpar,sizeof(double*));
      for(i=0;i<run.nMCMCpar;i++) matrix[i]  = (double*)calloc(run.nMCMCpar,sizeof(double));
      allocparset(&par, networksize);
      //getparameterset(&par, 3.0,0.11,700009012.346140,3.0, 0.5,0.9,3.0,0.5, 1.0,0.1,2.0,3.0);
      
      //computeFishermatrixIFO(par,run.nMCMCpar,network,networksize,0,matrix);
      //computeFishermatrix(&par,run.nMCMCpar,network,networksize,matrix);
      
      for(i=0;i<run.nMCMCpar;i++) {
	for(j=0;j<run.nMCMCpar;j++) {
	  printf("  %12.4g",matrix[i][j]/(double)networksize);
	}
	printf("\n\n");
      }
      printf("\n\n");
      
      
      
      freeparset(&par);
      for(i=0;i<run.nMCMCpar;i++) {
	free(matrix[i]);
      }
      free(matrix);
    }
  } //if(doMatch=1)

  
  
    
  
  //Get rid of allocated memory and quit
  for(ifonr=0; ifonr<networksize; ++ifonr) ifodispose(network[ifonr]);
  freeparset(&dummypar);
  
  
  clock_t time3 = clock();
  printf("   Timimg:\n");
  if(1==2) { 
    if(doMCMC>=1) {
      printf("     initialisation:%10.2lfs\n", ((double)time1 - (double)time0)*1.e-6 );
      printf("     MCMC:          %10.2lfs\n", ((double)time2 - (double)time1)*1.e-6 );
    }
    printf("     total time:    %10.2lfs\n", (double)time3*1.e-6);
  }
  
  printf("\n   SPINspiral done.\n\n");
  if(doMCMC>=1) printf("\n");
  return 0;
}



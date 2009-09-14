/* 
   
   SPINspiral:                parameter estimation on binary inspirals detected by LIGO, including spins of the binary members
   SPINspiral_signal.c:       routines to calculate likelihood, SNR, match, etc.
   
   
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



#include <SPINspiral.h>


/**
 * \file SPINspiral_signal.c
 * \brief Contains routines compute signal-related information, like likelihood, SNR, match, etc.
 */



// ****************************************************************************************************************************************************  
/**
 * \brief Compute the log(Likelihood) for a network of IFOs
 */
// ****************************************************************************************************************************************************  
double netLogLikelihood(struct parSet *par, int networkSize, struct interferometer *ifo[], int waveformVersion, int injectionWF, struct runPar run)
{
  double result = 0.0;
  int i;
  for (i=0; i<networkSize; ++i){
    result += IFOlogLikelihood(par, ifo, i, waveformVersion, injectionWF, run);
  }
  return result;
} // End of netLogLikelihood()
// ****************************************************************************************************************************************************  



// ****************************************************************************************************************************************************  
/**
 * \brief Compute the log(Likelihood) for a single IFO
 */
// ****************************************************************************************************************************************************  
double IFOlogLikelihood(struct parSet *par, struct interferometer *ifo[], int ifonr, int waveformVersion, int injectionWF, struct runPar run)
{
  int j=0;
  
  // Fill ifo[ifonr]->FTin with time-domain template:
  waveformTemplate(par, ifo, ifonr, waveformVersion, injectionWF, run);
  
  // Window template, FTwindow is a Tukey window:
  for(j=0; j<ifo[ifonr]->samplesize; ++j) 
    ifo[ifonr]->FTin[j] *= ifo[ifonr]->FTwindow[j];
  
  // Execute Fourier transform of signal template:
  fftw_execute(ifo[ifonr]->FTplan);
  
  // Compute the overlap between waveform and data:
  double overlaphd = vecOverlap(ifo[ifonr]->raw_dataTrafo, 
				ifo[ifonr]->FTout, ifo[ifonr]->noisePSD,
				ifo[ifonr]->lowIndex, ifo[ifonr]->highIndex, ifo[ifonr]->deltaFT);
  //correct FFT for sampling rate of waveform
  overlaphd/=((double)ifo[ifonr]->samplerate);  
  
  // Compute the overlap between waveform and itself:
  double overlaphh = vecOverlap(ifo[ifonr]->FTout,
				ifo[ifonr]->FTout, ifo[ifonr]->noisePSD,
				ifo[ifonr]->lowIndex, ifo[ifonr]->highIndex, ifo[ifonr]->deltaFT);
  //correct FFT for sampling rate of waveform
  overlaphh/=((double)ifo[ifonr]->samplerate);
  overlaphh/=((double)ifo[ifonr]->samplerate);
  
  return (overlaphd-0.5*overlaphh);
  
  /*
  //Alternative: about 8% slower because of extra copies of FFT output
  //   and division of each FFT point by ((double)ifo[ifonr]->samplerate)
  fftw_complex *FFTwaveform =
  fftw_malloc(sizeof(fftw_complex) * (ifo[ifonr]->FTsize));
  signalFFT(FFTwaveform, par, ifo, ifonr, waveformVersion, injectionWF, run);                             
  // Compute the overlap between waveform and data:
  double overlaphd = vecOverlap(ifo[ifonr]->raw_dataTrafo,
  FFTwaveform, ifo[ifonr]->noisePSD,
  ifo[ifonr]->lowIndex, ifo[ifonr]->highIndex, ifo[ifonr]->deltaFT);
  // Compute the overlap between waveform and itself:
  double overlaphh = vecOverlap(
  FFTwaveform, FFTwaveform, ifo[ifonr]->noisePSD,
  ifo[ifonr]->lowIndex, ifo[ifonr]->highIndex, ifo[ifonr]->deltaFT);
  return (overlaphd-0.5*overlaphh);
  */
  
} // End IFOlogLikelihood()
// ****************************************************************************************************************************************************  





// ****************************************************************************************************************************************************  
/**
 * \brief Compute the SNR of the waveform with a given parameter set for a single IFO
 */
// ****************************************************************************************************************************************************  
double signalToNoiseRatio(struct parSet *par, struct interferometer *ifo[], int ifonr, int waveformVersion, int injectionWF, struct runPar run)
// SNR of signal corresponding to parameter set, w.r.t. i-th interferometer's noise.
// (see SNR definition in Christensen/Meyer/Libson (2004), p.323)
{
  int j=0;
  
  // Fill ifo[ifonr]->FTin with time-domain template:
  waveformTemplate(par, ifo, ifonr, waveformVersion, injectionWF, run);
  
  // Window template, FTwindow is a Tukey window:
  for(j=0; j<ifo[ifonr]->samplesize; ++j)
    ifo[ifonr]->FTin[j] *= ifo[ifonr]->FTwindow[j];
  
  // Execute Fourier transform of signal template:
  fftw_execute(ifo[ifonr]->FTplan);
  
  // Compute the overlap between waveform and itself:
  double overlaphh = vecOverlap(ifo[ifonr]->FTout,
				ifo[ifonr]->FTout, ifo[ifonr]->noisePSD,
				ifo[ifonr]->lowIndex, ifo[ifonr]->highIndex, ifo[ifonr]->deltaFT);
  
  // Correct FFT for sampling rate of waveform
  overlaphh /= ((double)ifo[ifonr]->samplerate);
  overlaphh /= ((double)ifo[ifonr]->samplerate);
  
  //printf("\n\n\n DELTAFT: %10.4f  %10d %10d  %10d\n\n\n" ,ifo[ifonr]->deltaFT,ifo[ifonr]->samplesize,ifo[ifonr]->lowIndex, ifo[ifonr]->highIndex);
  
  return sqrt(overlaphh);
  
} // End signalToNoiseRatio()
// ****************************************************************************************************************************************************  






// ****************************************************************************************************************************************************  
/**
 * \brief Compute match between waveforms with parameter sets par1 and par2
 */
// ****************************************************************************************************************************************************  
double parMatch(struct parSet* par1, int waveformVersion1, int injectionWF1, struct parSet* par2, int waveformVersion2, int injectionWF2, struct interferometer *ifo[], int networkSize, struct runPar run)
{
  double overlap11=0.0, overlap12=0.0, overlap22=0.0;
  int ifonr;
  fftw_complex *FFT1, *FFT2; 
  
  
  /*
    printf("  parMatch:   %i  %i    %i  %i    %i\n", waveformVersion1, injectionWF1, waveformVersion2, injectionWF2, networkSize);
    int i=0;
    for(i=0;i<12;i++) {
    printf("  %i %f",i,par1->par[i]);
    }
    printf("\n");
    for(i=0;i<12;i++) {
    printf("  %i %f",i,par2->par[i]);
    }
    printf("\n");
  */
  
  
  for(ifonr=0; ifonr<networkSize; ifonr++){
    FFT1 = fftw_malloc(sizeof(fftw_complex) * (ifo[ifonr]->FTsize));
    FFT2 = fftw_malloc(sizeof(fftw_complex) * (ifo[ifonr]->FTsize));
    
    signalFFT(FFT1, par1, ifo, ifonr, waveformVersion1, injectionWF1, run);
    signalFFT(FFT2, par2, ifo, ifonr, waveformVersion2, injectionWF2, run);
    
    overlap11 += vecOverlap(FFT1, FFT1, ifo[ifonr]->noisePSD, ifo[ifonr]->lowIndex, ifo[ifonr]->highIndex, ifo[ifonr]->deltaFT);
    overlap12 += vecOverlap(FFT1, FFT2, ifo[ifonr]->noisePSD, ifo[ifonr]->lowIndex, ifo[ifonr]->highIndex, ifo[ifonr]->deltaFT);
    overlap22 += vecOverlap(FFT2, FFT2, ifo[ifonr]->noisePSD, ifo[ifonr]->lowIndex, ifo[ifonr]->highIndex, ifo[ifonr]->deltaFT);
  }
  
  double match = overlap12/sqrt(overlap11*overlap22);
  
  free(FFT1);
  free(FFT2);
  
  return match;
  
} // End of parMatch()
// ****************************************************************************************************************************************************  






// ****************************************************************************************************************************************************  
/**
 * \brief Compute frequency-domain overlap of waveform of given parameters with raw data
 */
// ****************************************************************************************************************************************************  
double overlapWithData(struct parSet *par, struct interferometer *ifo[], int ifonr, int waveformVersion, int injectionWF, struct runPar run)
{
  fftw_complex *FFTwaveform = fftw_malloc(sizeof(fftw_complex) * (ifo[ifonr]->FTsize));
  signalFFT(FFTwaveform, par, ifo, ifonr, waveformVersion, injectionWF, run);
  
  double overlap = vecOverlap(ifo[ifonr]->raw_dataTrafo, FFTwaveform, ifo[ifonr]->noisePSD, 
                              ifo[ifonr]->lowIndex, ifo[ifonr]->highIndex, ifo[ifonr]->deltaFT);
  
  free(FFTwaveform);
  return overlap;
} // End of overlapWithData()
// ****************************************************************************************************************************************************  





// ****************************************************************************************************************************************************  
/**
 * \brief Compute the overlap in the frequency domain between two waveforms with parameter sets par1 and par2
 */
// ****************************************************************************************************************************************************  
double parOverlap(struct parSet* par1, int waveformVersion1, int injectionWF1, struct parSet* par2, int waveformVersion2, int injectionWF2, struct interferometer* ifo[], int ifonr, struct runPar run)
{
  double overlap = 0.0;
  fftw_complex *FFT1 = fftw_malloc(sizeof(fftw_complex) * (ifo[ifonr]->FTsize));
  fftw_complex *FFT2 = fftw_malloc(sizeof(fftw_complex) * (ifo[ifonr]->FTsize));
  
  // Get waveforms, FFT them and store them in FFT1,2
  signalFFT(FFT1, par1, ifo, ifonr, waveformVersion1, injectionWF1, run);
  signalFFT(FFT2, par2, ifo, ifonr, waveformVersion2, injectionWF2, run);
  
  // Compute the overlap between the vectors FFT1,2, between index i1 and i2:
  overlap = vecOverlap(FFT1, FFT2, ifo[ifonr]->noisePSD, ifo[ifonr]->lowIndex, ifo[ifonr]->highIndex, ifo[ifonr]->deltaFT);
  
  fftw_free(FFT1);
  fftw_free(FFT2);
  
  return overlap;
  
} // End of parOverlap()
// ****************************************************************************************************************************************************  





// ****************************************************************************************************************************************************  
/**
 * \brief Compute the overlap in the frequency domain between two waveforms with parameter sets par1 and par2
 */
// ****************************************************************************************************************************************************  
double vecOverlap(fftw_complex *vec1, fftw_complex *vec2, double * noise, int j1, int j2, double deltaFT)
//Compute the overlap of vectors vec1 and vec2, between indices j1 and j2
{
  int j=0;
  double overlap = 0.0;
  
  // Sum over the Fourier frequencies within operational range (some 40-1500 Hz or similar):
  for (j=j1; j<=j2; ++j){
    overlap += 4.0*creal( vec1[j]*conj(vec2[j]) / noise[j-j1] ) / deltaFT;   
  }
  return overlap;
  
}
// ****************************************************************************************************************************************************  





// ****************************************************************************************************************************************************  
/**
 * \brief Compute the FFT of a waveform with given parameter set
 */
// ****************************************************************************************************************************************************  
void signalFFT(fftw_complex* FFTout, struct parSet* par, struct interferometer* ifo[], int ifonr, int waveformVersion, int injectionWF, struct runPar run)
{
  int j=0;
  if(FFTout==NULL) {
    fprintf(stderr,"\n\n   ERROR: memory should be allocated for FFTout vector before call to signalFFT()\n   Aborting...\n\n");
    exit(1);
  }
  
  //printf("  signalFFT:  %i  %i\n", waveformVersion, injectionWF);
  
  // Fill ifo[i]->FTin with time-domain template:
  waveformTemplate(par, ifo, ifonr, waveformVersion, injectionWF, run);
  
  // Window template, FTwindow is a Tukey window:
  for(j=0; j<ifo[ifonr]->samplesize; ++j) ifo[ifonr]->FTin[j] *= ifo[ifonr]->FTwindow[j];
  
  // Execute Fourier transform of signal template:
  fftw_execute(ifo[ifonr]->FTplan);
  
  for(j=0; j<ifo[ifonr]->FTsize; j++) FFTout[j] = ifo[ifonr]->FTout[j]/((double)ifo[ifonr]->samplerate);
  
} // End of SignalFFT()
// ****************************************************************************************************************************************************  



// ****************************************************************************************************************************************************  
/**
 * \brief Compute the match between a given parameter array and the injection parameters
 * \todo Add support for two different waveform families
 */
// ****************************************************************************************************************************************************  
double matchBetweenParameterArrayAndTrueParameters(double * pararray, struct interferometer *ifo[], struct MCMCvariables mcmc, struct runPar run)
{
  struct parSet par, injectPar;
  int i=0;
  int injectionWF = 1;                                     // Call localPar or parMatch with the injection template
  int mcmcWF = 0;                                          // Call localPar or parMatch with the MCMC template
  for(i=0;i<mcmc.nMCMCpar;i++) {
    par.par[i] = pararray[i];
  }
  allocParset(&par,mcmc.networkSize);
  localPar(&par, ifo, mcmc.networkSize, mcmcWF, run);  // Call localPar for an MCMC template
  
  //Get the injection parameters:
  getInjectionParameters(&injectPar, mcmc.nInjectPar, mcmc.injParVal);
  allocParset(&injectPar,mcmc.networkSize);
  localPar(&injectPar, ifo, mcmc.networkSize, injectionWF, run);  // Call localPar for an injection template
  
  freeParset(&par);
  freeParset(&injectPar);
  
  return parMatch(&injectPar, mcmc.injectionWaveform, injectionWF, &par, mcmc.mcmcWaveform, mcmcWF, ifo, mcmc.networkSize, run);
} // End of matchBetweenParameterArrayAndTrueParameters()
// ****************************************************************************************************************************************************  





// ****************************************************************************************************************************************************  
/**
 * \brief Compute Fisher matrix for parameter set par for a given IFO
 */
// ****************************************************************************************************************************************************  
/*
  void computeFisherMatrixIFO(struct parSet *par, int nParameters, struct interferometer *ifo[], int networkSize, int ifonr, double **matrix)
  {
  int ip=0,jp=0,j=0,j1=0,j2=0,nFT=0;
  struct parSet par1;
  allocParset(&par1, networkSize);
  double pars[nParameters];
  
  nFT = ifo[ifonr]->FTsize;
  double dx[nParameters], noise[nFT];
  double _Complex FFT0[nFT], FFT1[nFT], dFFTs[nParameters][nFT];
  
  j1 = ifo[ifonr]->lowIndex;
  j2 = ifo[ifonr]->highIndex;
  
  // Compute the FFTed signal for the default parameter set FFT0
  signalFFT(FFT0, par, ifo, ifonr, waveformVersion, injectionWF, run);
  
  for(ip=0;ip<nParameters;ip++) {
  // Change parameter ip with dx
  dx[ip] = 1.e-5;       // Should this be the same for each parameter?
  par2arr(par,pars);    // Stick the default parameter set par into the array pars
  pars[ip] += dx[ip];   // Change parameter ip with dx
  arr2par(pars,&par1);  // Put the changed parameter set into struct par1
  
  // Compute the FFTed signal for this parameter set FFT1
  signalFFT(par1, ifo, networkSize, ifonr, FFT1, waveformVersion, injectionWF, run);
  
  // Compute the partial derivative to parameter ip
  for(j=j1;j<=j2;j++) {
  dFFTs[ip][j] = (FFT1[j]-FFT0[j])/dx[ip];
  }
  }
  
  
  // Compute the actual Fisher matrix (diagonal + lower triangle)
  for(ip=0;ip<nParameters;ip++) {
  for(jp=0;jp<=ip;jp++) {
  matrix[ip][jp] = vecOverlap(dFFTs[ip], dFFTs[jp], 
  ifo[ifonr]->noisePSD, j1, j2, ifo[ifonr]->deltaFT);
  }
  }
  
  // Copy the lower to the upper triangle to get a complete matrix
  for(ip=0;ip<nParameters;ip++) {
  for(jp=ip;jp<nParameters;jp++) {
  matrix[ip][jp] = matrix[jp][ip];
  }
  }
  
  freeParset(&par1);
  } // End of computeFisherMatrixIFO()
  // ****************************************************************************************************************************************************  
  */

// ****************************************************************************************************************************************************  
/**
 * \brief Compute the Fisher matrix for a network of IFOs.
 * 
 * Use computeFisherMatrixIFO() to compute the elements per IFO
 */
// ****************************************************************************************************************************************************  
/*
  void computeFisherMatrix(struct parSet *par, int nParameters, struct interferometer *ifo[], int networkSize, double **matrix)
  {
  int ip=0,jp=0,ifonr=0;
  double **dmatrix  = (double**)calloc(nParameters,sizeof(double*));
  for(ip=0;ip<nParameters;ip++) dmatrix[ip]  = (double*)calloc(nParameters,sizeof(double));
  
  for(ifonr=0;ifonr<networkSize;ifonr++) {
  computeFisherMatrixIFO(par,nParameters,ifo,networkSize,ifonr,dmatrix);
  
  for(ip=0;ip<nParameters;ip++) {
  for(jp=0;jp<nParameters;jp++) {
  matrix[ip][jp] += dmatrix[ip][jp];
  }
  }
  }
  } // End of computeFisherMatrix()
  // ****************************************************************************************************************************************************  
  */

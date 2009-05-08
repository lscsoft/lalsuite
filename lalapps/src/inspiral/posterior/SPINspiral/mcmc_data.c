/* 
   
   SPINspiral:                parameter estimation on binary inspirals detected by LIGO, including spins of the binary members
   mcmc_data.c:               data-handling routines: data I/O, noise PSD, filtering, windowing, downsampling, FFTing
   
   
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


// *** Routines that handle IFOs ***
void setIFOdata(struct runPar *run, struct interferometer ifo[])
// Set all the data for all IFOs that may be used
{
  int i=0,numberofdatasets = 6;
  // Description of the data sets below
  char datadescriptions[10][99];
  sprintf(datadescriptions[0],"File:  %s",run->dataFilename);
  sprintf(datadescriptions[1],"Gaussian, stationary noise (GPS ~894377000)");
  sprintf(datadescriptions[2],"clean S5 data (GPS ~846226044)");
  sprintf(datadescriptions[3],"playground trigger data (GPS ~845348295)");
  sprintf(datadescriptions[4],"glitchy data (GPS ~846471090)");
  sprintf(datadescriptions[5],"NINJA");
  sprintf(datadescriptions[6],"Original Gaussian noise (GPS ~700006000)");
  
  run->selectdata = 0;  //Now, you *have* to read from a file
  
  //run->selectdata = max(min(run->selectdata,numberofdatasets),1);
  if(run->selectdata < 0 || run->selectdata > numberofdatasets) {
    printf("\n\n   Unknown dataset %d selected.  Please set the parameter  selectdata  in  %s  to a value between 0 and %d: \n",run->selectdata, run->mainFilename,numberofdatasets);
    for(i=0;i<=numberofdatasets;i++) {
      printf("     %3d:  %s\n",i,datadescriptions[i]);
    }
    printf("\n\n");
    exit(1);
  }
  
  
  // Print selected type of data/noise:
  //printf("   Using data set %d: %s.\n",run->selectdata,datadescriptions[run->selectdata]);
  
  
  // *** Read detector position, orientation, file and channel names for detector strains from file (e.g. mcmc.data):
  if(run->selectdata == 0) readDataInputfile(run,ifo);        //Read data on the noise files to use
  
  // Print selected type of data/noise:
  printf("   Data set used: %s\n",run->datasetName);
  
  // WGS-84 data:
  for(i=0;i<run->maxIFOdbaseSize;i++) {
    ifo[i].radius_eqt  = 6378137.0;
    ifo[i].radius_pole = 6356752.314;
  }
  
  
  // *** Data time and frequency cutoffs per detector (same for all detectors for now):
  for(i=0;i<run->maxIFOdbaseSize;i++) {
    ifo[i].lowCut    = run->lowfrequencycut;   // Define lower and upper limits of overlap integral
    ifo[i].highCut   = run->highfrequencycut;
    ifo[i].before_tc = run->databeforetc;      // Define data segment: [t_c-before_tc, t_c+after_tc]
    ifo[i].after_tc  = run->dataaftertc;
  }
  
  // ***  Initialise hard-coded data sets ***
  // ***  Removed at revision 148
}
// End of set_ifo_data()





void ifoinit(struct interferometer **ifo, int networksize, struct runPar run)
// Determines interferometer arm (unit-) vectors (and others)                  
// given position (lat/long) and arm angles.                                   
// Vectors refer to the (right-handed) earth                                    
// coordinate system spanned by the three vectors:                             
//   x) from geocenter to intersection of greenwich meridian with equator plane
//   y) from geocenter to intersection of 90E meridian with equator plane      
//   z) from geocenter to north pole                                           
{
  double merinormal[3];  // Normal vector of meridian plane
  char latchar[2];
  char longchar[2];
  int ifonr,j;
  double f;
  double flattening, eccentricitySQ, curvatureradius;
  for(ifonr=0; ifonr<networksize; ++ifonr){
    ifo[ifonr]->index = ifonr;
    // Some text output...
    if(intscrout==1) printf(" | Interferometer %d: `%s'", ifo[ifonr]->index+1, ifo[ifonr]->name);
    if((ifo[ifonr]->lati)  < 0.0) sprintf(latchar, "S");
    else sprintf(latchar, "N");
    if((ifo[ifonr]->longi) < 0.0) sprintf(longchar, "W");
    else sprintf(longchar, "E");
    if(intscrout==1) printf(" at  %1.0f*%2.1f'%s  %1.0f*%2.1f'%s  (%3.0f/%3.0f)\n",
           floor(fabs(ifo[ifonr]->lati*r2d)),  (fabs(ifo[ifonr]->lati*r2d)-floor(fabs(ifo[ifonr]->lati*r2d)))*60.0, latchar,
           floor(fabs(ifo[ifonr]->longi*r2d)), (fabs(ifo[ifonr]->longi*r2d)-floor(fabs(ifo[ifonr]->longi*r2d)))*60.0, longchar,
           360.0 - ifo[ifonr]->rightarm*r2d, 360.0 - ifo[ifonr]->leftarm*r2d);
    if(ifo[ifonr]->ch1doubleprecision && intscrout==1) printf(" | frame file precision: double (64 bit)\n"); 
    else if(intscrout==1) printf(" | frame file precision: float (32 bit)\n"); 
    if(intscrout==1) printf(" | frequency range: %.0f to %.0f Hz.\n", ifo[ifonr]->lowCut, ifo[ifonr]->highCut);
    if(intscrout==1) printf(" | initialising vectors etc...");
    
    
    //Longitude: East is positive
    
    // Place arms on equator plane, so that its designated N-S-direction is aligned with its meridian plane:
    ifo[ifonr]->rightvec[0]  = -cos(ifo[ifonr]->longi + ifo[ifonr]->rightarm);
    ifo[ifonr]->rightvec[1]  = -sin(ifo[ifonr]->longi + ifo[ifonr]->rightarm);
    ifo[ifonr]->rightvec[2]  = 0.0;
    
    ifo[ifonr]->leftvec[0]   = -cos(ifo[ifonr]->longi + ifo[ifonr]->leftarm);
    ifo[ifonr]->leftvec[1]   = -sin(ifo[ifonr]->longi + ifo[ifonr]->leftarm);
    ifo[ifonr]->leftvec[2]   = 0.0;
    
    ifo[ifonr]->normalvec[0] = 0.0;
    ifo[ifonr]->normalvec[1] = 0.0;
    ifo[ifonr]->normalvec[2] = 1.0;
    
    // The following vector is rightarm + 90deg and usually (but not necessarily) identical to the left arm (leftvec):
    ifo[ifonr]->orthoarm[0]  = -cos(ifo[ifonr]->longi + ifo[ifonr]->rightarm + 0.5*pi);
    ifo[ifonr]->orthoarm[1]  = -sin(ifo[ifonr]->longi + ifo[ifonr]->rightarm + 0.5*pi);
    ifo[ifonr]->orthoarm[2]  = 0.0;
    
    // Determine normal vector of meridian plane (i.e., the vector that points E(?) when standing at the equator at longi):
    merinormal[0] = cos(ifo[ifonr]->longi - 0.5*pi);
    merinormal[1] = sin(ifo[ifonr]->longi - 0.5*pi);
    merinormal[2] = 0.0;
    
    // The three vectors:                                                     
    //   x) from geocenter to intersection of ifo meridian with equator plane
    //   y) from geocenter to north pole
    //   z) the above normal vector (merinormal(?))
    // again form another (orthonormal) right-handed system.
    // Now turn all arms clockwise around the normal vector of meridian plane, to account for the latitude of the detectors:
    rotate(ifo[ifonr]->rightvec,  pi/2.0 - ifo[ifonr]->lati, merinormal);
    rotate(ifo[ifonr]->leftvec,   pi/2.0 - ifo[ifonr]->lati, merinormal);
    rotate(ifo[ifonr]->orthoarm,  pi/2.0 - ifo[ifonr]->lati, merinormal);
    rotate(ifo[ifonr]->normalvec, pi/2.0 - ifo[ifonr]->lati, merinormal);
    
    // Initialise the ifo position (!NOT! unit-) vector:
    coord2vec(sin(ifo[ifonr]->lati), ifo[ifonr]->longi, ifo[ifonr]->positionvec);
    if(ifo[ifonr]->radius_eqt < ifo[ifonr]->radius_pole) printf("  CHECK EARTH MODEL RADII !!  ");
    flattening      = (ifo[ifonr]->radius_eqt - ifo[ifonr]->radius_pole) / ifo[ifonr]->radius_eqt;
    eccentricitySQ  = flattening*(2.0-flattening);  // (squared eccentricity)
    curvatureradius = ifo[ifonr]->radius_eqt / sqrt(1.0-eccentricitySQ*pow(sin(ifo[ifonr]->lati),2.0));
    ifo[ifonr]->positionvec[0] *= curvatureradius;
    ifo[ifonr]->positionvec[1] *= curvatureradius;
    ifo[ifonr]->positionvec[2] *= curvatureradius*(1.0-eccentricitySQ);
    if(intscrout==1) printf(" ok.\n");
    // printf("== normalvec (%1.2f, %1.2f, %1.2f)\n", ifo[ifonr]->normalvec[0],ifo[ifonr]->normalvec[1],ifo[ifonr]->normalvec[2]);
    // if(intscrout==1) printf(" : f=%f  e^2=%f  v=%f \n", flattening, eccentricitySQ, curvatureradius);
    
    
    // Read 'detector' noise and estimate PSD
    if(printMuch>=1) printf("   Reading noise for the detector in %s...\n",ifo[ifonr]->name);
    noisePSDestimate(ifo[ifonr],run);
  
  
    // Read 'detector' data for injection
    if(printMuch>=1) printf("   Reading data for the detector in %s...\n",ifo[ifonr]->name);
    dataFT(ifo,ifonr,networksize,run);
    
    
    
    // Initialise array of different powers of Fourier frequencies corresponding to the elements of 'ifo[ifonr]->dataTrafo':       
    // First loop to determine index bounds & range:               
    ifo[ifonr]->lowIndex = 0; ifo[ifonr]->highIndex=0;
    for(j=1; j<ifo[ifonr]->FTsize; ++j){
      f = (((double)j)/((double)ifo[ifonr]->deltaFT));
      if((ifo[ifonr]->lowIndex==0)  && (f>=ifo[ifonr]->lowCut)) ifo[ifonr]->lowIndex = j;
      if((ifo[ifonr]->highIndex==0) && (f>ifo[ifonr]->highCut)) ifo[ifonr]->highIndex = j-1;
      // ...so 'lowIndex' and 'highIndex' are the extreme indexes WITHIN frequency band
    }
    ifo[ifonr]->indexRange = ifo[ifonr]->highIndex - (ifo[ifonr]->lowIndex - 1);
    // 'lowIndex' and 'highIndex' are the indices analogous to 'lowCut' and 'highCut'
    // but can be used to access the respective elements of 'raw_dataTrafo'.
    
    ifo[ifonr]->noisePSD   = ((double*) malloc(sizeof(double) * ifo[ifonr]->indexRange));
    ifo[ifonr]->dataTrafo  = ((fftw_complex*) malloc(sizeof(fftw_complex) * ifo[ifonr]->indexRange));
    for(j=0; j<ifo[ifonr]->indexRange; ++j){
      f = (((double)(j+ifo[ifonr]->lowIndex))/((double)ifo[ifonr]->deltaFT));
      ifo[ifonr]->noisePSD[j]      = interpol_log_noisePSD(f,ifo[ifonr]);
      
      // Although smoothing was done for log noise, we store real noise on output
      ifo[ifonr]->noisePSD[j] = exp(ifo[ifonr]->noisePSD[j]);
      ifo[ifonr]->dataTrafo[j]  = ifo[ifonr]->raw_dataTrafo[j+ifo[ifonr]->lowIndex];
    }
    if(intscrout==1) printf(" | %d Fourier frequencies within operational range %.0f--%.0f Hz.\n", ifo[ifonr]->indexRange, ifo[ifonr]->lowCut, ifo[ifonr]->highCut);
    if(ifonr<networksize-1 && intscrout==1) printf(" | --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --\n");
  } //for(ifonr=0; ifonr<networksize; ++ifonr)
}
// End of ifoinit()




void ifodispose(struct interferometer *ifo)
{
  if(intscrout==1) printf(" | Interferometer %d `%s' is taken offline.\n", ifo->index, ifo->name);
  free(ifo->raw_noisePSD);       ifo->raw_noisePSD = NULL;
  fftw_free(ifo->raw_dataTrafo); ifo->raw_dataTrafo = NULL;
  free(ifo->noisePSD);           ifo->noisePSD = NULL;
  free(ifo->dataTrafo);          ifo->dataTrafo = NULL;
  fftw_destroy_plan(ifo->FTplan);
  fftw_free(ifo->FTin);          ifo->FTin = NULL;
  fftw_free(ifo->rawDownsampledWindowedData); ifo->rawDownsampledWindowedData = NULL;  
  fftw_free(ifo->FTout);         ifo->FTout = NULL;
  free(ifo->FTwindow);           ifo->FTwindow = NULL;
}






// *** Routines that do data I/O and data handling ***

double *filter(int *order, int samplerate, double upperlimit)
// Computes FIR filter coefficients.                               
// The filter has a total of N=(order+order-1) coefficients        
// that are symmetric, i.e. coef[k] = coef[N-k].                   
// For details & filter properties see downsampling function below.
{
  // Specify filter characteristics:
  //int ncoef      = 129;  // number of filter coefficients... 129 should be sufficient
  int ncoef      = 129;  // number of filter coefficients... 129 should be sufficient
  int totalcoef  = ncoef+ncoef-1;   // total number of coefficients                  
  double desired[2] = {1.0, 0.0};      // desired gain                               
  double weights[2] = {1.0, 1.0};      // weight for `loss' in pass- & stopband      
  //double transitionbandwidth=0.0125;    //0.0125 was suggested by Christian Roever via 07/30/08 e-mail
  double transitionbandwidth=0.025;    
  // Place transition bandwidth half-way between upper edge of pass band, which is
  // (upperlimit/samplerate) in relative units, and new Nyquist frequency, which is
  // 0.5/downsamplefactor in relative units.
  // band vector contains 0, end of pass band, end of transition band, 0.5 (old Nyquist)
  // if there's not enough room for a transition band, we have a problem!
  
  if(0.5/(double)downsamplefactor - upperlimit/((double)samplerate) < transitionbandwidth){
    printf(" Problem in filter() function while downsampling!\n");
    printf(" New Nyquist frequency after downsampling by %d is %g\n",
	   downsamplefactor, ((double)samplerate)/(double)downsamplefactor/2.0);
    printf(" Desired upper limit is %g\n", upperlimit);
    printf(" This doesn't leave enough room for a transition band of relative width %g\n",
	   transitionbandwidth);
    printf(" Maximum upper limit: %g Hz.\n", (double)samplerate * 
	   (0.5/(double)downsamplefactor - transitionbandwidth));
    printf(" Aborting!\n");
    exit(1);
  }
  
  double endpassband= upperlimit/((double)samplerate) + (0.5/(double)downsamplefactor - upperlimit/((double)samplerate) - transitionbandwidth)/2.0;
  double endtransitionband=0.5/(double)downsamplefactor -  (0.5/(double)downsamplefactor - upperlimit/((double)samplerate)  - transitionbandwidth)/2.0;
  double bands[4]   = {0.0, endpassband, endtransitionband, 0.5};
  double *coef;                        // Vector of coefficients (symmetric)
  coef = (double*) malloc(sizeof(double)*totalcoef);
  // Determine filter coefficients:
  remez(coef, totalcoef, 2, bands, desired, weights, BANDPASS);
  *order = ncoef;
  return coef;
}





double *downsample(double data[], int *datalength, double filtercoef[], int ncoef)
// Downsamples a time series by factor downsamplefactor by first low-pass
// filtering it using a finite-impulse-response (FIR) filter
// and then thinning the data.                              
// Filter coefficients are determined using the             
// `Parks-McClellan' or `Remez exchange' algorithm.         
// Resulting data vector is shorter than original           
// Returned vector is allocated using `fftw_malloc()' and   
// thus must be freed again using `fftw_free()'.            
{
  int flength = *datalength-2*(ncoef-1);
  int tlength = (int)ceil(((double)flength)/(double)downsamplefactor);
  double *thinned;      // Vector of filtered & thinned data
  int i,j,k;
  
  // Filter & thin data:
  thinned = (double*) fftw_malloc(sizeof(double) * tlength);
  k = 0;
  for(i=ncoef-1; i<*datalength-ncoef; i+=downsamplefactor) {
    thinned[k] = filtercoef[ncoef-1]*data[i];
    for(j=1; j<ncoef; ++j)
      thinned[k] += filtercoef[(ncoef-1)+j]*(data[i-j]+data[i+j]);
    ++k;
  }
  
  // Return results:
  *datalength = tlength;
  return thinned;
}



double hann(int j, int N)
// `hann window' for windowing data 
// j = 0, ..., N-1                  
{
  return 0.5*(1.0-cos(((double)j/(double)N)*2.0*pi));
}

double tukey(int j, int N, double r)
// Tukey window... for r=0 equal to rectangular window, for r=1 equal to Hann window.
// ( 0 < r < 1 denotes the fraction of the window in which it behaves sinusoidal)    
// j = 0, ..., N-1                                                                   
{
  double win = 1.0;
  if(((double)j) > (((double)N)/2.0)) j = N-j;
  if(((double)j) < (r*(((double)N)/2.0)))
    win = 0.5*(1.0-cos(((2.0*pi)/r)*(((double)j)/((double)N))));
  return win;
}




void dataFT(struct interferometer *ifo[], int ifonr, int networksize, struct runPar run)
// Read the data and do a software injection
// Computes the Fourier Transform for the specified range of the specified Frame (".gwf") file,
// after adding up the two (signal & noise) channels, or injecting a waveform template into the noise.
// Also takes care of preparing FT stuff  (ifo[ifonr]->FTplan, ->FTin, ->FTout, ...).
{
  struct FrFile *iFile=NULL;                // Frame file(s)
  struct FrVect *svect=NULL, *nvect=NULL;   // data vectors (signal & noise)
  int           N;                          // size of input
  double        *raw;                       // downsampling input
  int           j, ncoef;
  double        *filtercoef;
  long          filestart;
  char          filenames[1000]="";
  int           filecount = 0;
  double        from, to, delta;
  double        *injection;
  
  
  // `from' and `to' are determined so that the range specified by `before_tc' and `after_tc'
  // falls into the flat part of the (Tukey-) window:                                        
  from  = floor(prior_tc_mean - ifo[ifonr]->before_tc - (ifo[ifonr]->before_tc+ifo[ifonr]->after_tc) * 0.5 * (run.tukeywin/(1.0-run.tukeywin)));
  to    =  ceil(prior_tc_mean + ifo[ifonr]->after_tc  + (ifo[ifonr]->before_tc+ifo[ifonr]->after_tc) * 0.5 * (run.tukeywin/(1.0-run.tukeywin)));
  delta = (to) - (from);
  if(intscrout==1) printf(" | investigated time range : from %.1f to %.1f (%.1f seconds)\n", from, to, delta);
  
  
  // Starting time of first(!) Frame file to be read:
  filestart = (((((long)(from))-ifo[ifonr]->ch1fileoffset) / ifo[ifonr]->ch1filesize) * ifo[ifonr]->ch1filesize) + ifo[ifonr]->ch1fileoffset;
  
  
  // Assemble the filename character string:
  while (((double)filestart) < to){
    if(filecount == 0)
      sprintf(filenames, "%s/%s%ld%s", ifo[ifonr]->ch1filepath, ifo[ifonr]->ch1fileprefix, (long)filestart, ifo[ifonr]->ch1filesuffix);  // Fill in filename etc. for first file
    else
      sprintf(filenames, "%s %s/%s%ld%s", filenames, ifo[ifonr]->ch1filepath, ifo[ifonr]->ch1fileprefix, (long)filestart, ifo[ifonr]->ch1filesuffix);  // Append filename etc. for following files
    //if(intscrout==1) printf(" |   %s%ld%s\n", ifo[ifonr]->framefileprefix, (long)filestart, ifo[ifonr]->framefilesuffix);
    filestart += ifo[ifonr]->ch1filesize;
    filecount += 1;
  }
  
  
  
  // Open frame file(s):
  if(intscrout==1) printf(" | opening %d chirp data file(s)... \n", filecount);
  iFile = FrFileINew(filenames);
  if(iFile == NULL) {
    printf("\n\n   ERROR opening data file: %s (channel 1), aborting.\n\n\n",filenames);
    exit(1);
  }
  if(intscrout==1) printf(" | %s\n",filenames);
  
  // Read 1st channel (noise or noise+signal):
  if(ifo[ifonr]->ch1doubleprecision) //Seems precision is double
    nvect = FrFileIGetVectD(iFile, ifo[ifonr]->ch1name, from, delta);
  else
    nvect = FrFileIGetVectF(iFile, ifo[ifonr]->ch1name, from, delta);
  if(nvect == NULL) {
    printf("\n\n   ERROR reading data file: %s (channel 1), aborting.\n\n\n",filenames);
    exit(1);
  }
  
  FrFileIEnd(iFile);
  
  N = nvect->nData;
  ifo[ifonr]->samplerate = (int)(1.0 / (nvect->dx[0]) + 0.5);  // Add 0.5 for correct truncation/rounding
  if(intscrout==1) printf(" | original sampling rate: %d Hz\n", ifo[ifonr]->samplerate);
  
  // Inject the signal into the noise
  if(run.injectSignal >= 1) {
    if(intscrout==1) printf(" :  injecting signal:\n");
    
    // Define injection parameters:
    struct parset injectpar;
    getInjectionParameters(&injectpar, run.nInjectPar, run.injParVal);
    double m1=0.0,m2=0.0;
    mceta2masses(injectpar.mc, injectpar.eta, &m1, &m2);
    injectpar.loctc    = (double*)calloc(networksize,sizeof(double));
    injectpar.localti  = (double*)calloc(networksize,sizeof(double));
    injectpar.locazi   = (double*)calloc(networksize,sizeof(double));
    injectpar.locpolar = (double*)calloc(networksize,sizeof(double));
    
    if(intscrout==1) {
      printf(" :   m1 = %.1f Mo,  m2 = %.1f Mo  (Mc = %.3f Mo,  eta = %.4f)\n", m1, m2, injectpar.mc, injectpar.eta);
      printf(" :   tc = %.4f s,  dist = %.1f Mpc\n", injectpar.tc, exp(injectpar.logdl));
      printf(" :   ra = %.2f h,  dec = %.2f deg  (GMST = %.2f h)\n",(rightAscension(injectpar.longi,GMST(injectpar.tc))/pi)*12.0, (asin(injectpar.sinlati)/pi)*180.0, (GMST(injectpar.tc)/pi)*12.0);
      printf(" :   phase = %.2f rad\n", injectpar.phase);
    }
    ifo[ifonr]->FTstart = from; // Temporary setting so that localpar() works properly
    
    localpar(&injectpar, ifo, networksize);
    
    if(intscrout==1) {
      printf(" :   local parameters:\n");
      printf(" :   tc           = %.5f s\n",injectpar.loctc[ifonr]+from);
      printf(" :   altitude     = %.2f rad\n",injectpar.localti[ifonr]);
      printf(" :   azimuth      = %.2f rad\n",injectpar.locazi[ifonr]);
      printf(" :   polarisation = %.2f rad\n",injectpar.locpolar[ifonr]);
    }
    
    
    // Generate injection waveform template:
    injection = malloc(sizeof(double) * N);
    double *tempinj = ifo[ifonr]->FTin;
    double tempfrom = ifo[ifonr]->FTstart;
    int tempN = ifo[ifonr]->samplesize;
    ifo[ifonr]->FTin = injection;
    ifo[ifonr]->FTstart = from;
    ifo[ifonr]->samplesize = N;
    template(&injectpar,ifo,ifonr, run.injectionWaveform);
    ifo[ifonr]->FTin = tempinj;
    ifo[ifonr]->FTstart = tempfrom;
    ifo[ifonr]->samplesize = tempN;
    
    freeparset(&injectpar);
    
    
  } // if(run.injectSignal >= 1)
  else if(ifo[ifonr]->add2channels) { // Read 2nd channel (signal only)  (if not doing a software injection)
    
    filestart = (((((long)(from))-ifo[ifonr]->ch2fileoffset) / ifo[ifonr]->ch2filesize) * ifo[ifonr]->ch2filesize) + ifo[ifonr]->ch2fileoffset;
    
    // Assemble the filename character string:
    sprintf(filenames, " "); filecount = 0;
    while (((double)filestart) < to){
      if(filecount == 0)
        sprintf(filenames, "%s/%s%ld%s", ifo[ifonr]->ch2filepath, ifo[ifonr]->ch2fileprefix, (long)filestart, ifo[ifonr]->ch2filesuffix);  // Fill in filename etc. for first file
      else
        sprintf(filenames, "%s %s/%s%ld%s", filenames, ifo[ifonr]->ch2filepath, ifo[ifonr]->ch2fileprefix, (long)filestart, ifo[ifonr]->ch2filesuffix);  // Append filename etc. for following files:
      filestart += ifo[ifonr]->ch2filesize;
      filecount += 1;
    }
    
    // Open file:
    iFile = FrFileINew(filenames);
    if(iFile == NULL) {
      printf("\n\n   ERROR opening data file: %s (channel 2), aborting.\n\n\n",filenames);
      exit(1);
    }
    
    if(ifo[ifonr]->ch2doubleprecision) {
      svect = FrFileIGetVectD(iFile, ifo[ifonr]->ch2name, from, delta);
    } else {
      svect = FrFileIGetVectF(iFile, ifo[ifonr]->ch2name, from, delta);
    }
    if(svect == NULL) {
      printf("\n\n   ERROR reading data file: %s (channel 2), aborting.\n\n\n",filenames);
      exit(1);
    }
    FrFileIEnd(iFile);
  } //End if not doing a software injection
  
  
  // Allocate memory for transform input:
  raw  = malloc(sizeof(double) * N);
  // Fill in values:
  for(j=0; j<N; ++j) raw[j] = nvect->dataF[j];
  
  // Add channels (noise plus signal):
  if(run.injectSignal >= 1){
    for(j=0; j<N; ++j) raw[j] += injection[j];
  } else if(ifo[ifonr]->add2channels) {
    for(j=0; j<N; ++j) raw[j] += svect->dataF[j];
  }
  
  
  // Release the FrVect objects:
  FrVectFree(svect);
  FrVectFree(nvect);
  if(run.injectSignal >= 1) free(injection);
  
  
  int screwcount = 0;
  for(j=0; j<N; ++j)
    if(!(raw[j]<HUGE_VAL)) ++screwcount;
  if(screwcount>0){
    printf(" : %d missing data points in DATA file(s) !!\n",screwcount);
    printf(" : (maybe the precision is incorrect)\n");
  }
 
  // Downsample (by factor downsamplefactor):    *** changes value of N ***
  if(downsamplefactor!=1){
    if(intscrout==1) printf(" | downsampling... \n");
    filtercoef = filter(&ncoef, ifo[ifonr]->samplerate, ifo[ifonr]->highCut);
    
    /*
    //Print the filter coefficients:
    for(j=1; j<ncoef; ++j) {
      printf("  %d  %lf\n",j,filtercoef[ncoef-1+j]);
    }
    printf("\n\n");
    */
    
    ifo[ifonr]->FTin = downsample(raw, &N, filtercoef, ncoef);
    ifo[ifonr]->FTstart = from + ((double)(ncoef-1))/((double)(ifo[ifonr]->samplerate));
    ifo[ifonr]->deltaFT = delta - ((double)((ncoef-1)*2))/((double)(ifo[ifonr]->samplerate));
    ifo[ifonr]->samplesize = N;
    ifo[ifonr]->FTsize = (N/2)+1;
    ifo[ifonr]->samplerate = (int)((double)ifo[ifonr]->samplerate/(double)downsamplefactor);
    free(raw);
    free(filtercoef);
  }else{
    ifo[ifonr]->FTin = (double*) fftw_malloc(sizeof(double)*N);
    for (j=0; j<N; ++j)
      ifo[ifonr]->FTin[j] = raw[j];
    ifo[ifonr]->FTstart = from;
    ifo[ifonr]->deltaFT = delta;
    ifo[ifonr]->samplesize = N;
    ifo[ifonr]->FTsize = (N/2)+1;  
    free(raw);
  }
  
  // Window input data with a Tukey window:
  ifo[ifonr]->FTwindow = malloc(sizeof(double) * N);
  ifo[ifonr]->rawDownsampledWindowedData = (double*) fftw_malloc(sizeof(double)*N);
  for(j=0; j<N; ++j){
    ifo[ifonr]->FTwindow[j] =  tukey(j, N, run.tukeywin);
    ifo[ifonr]->FTin[j] *= ifo[ifonr]->FTwindow[j];
    ifo[ifonr]->rawDownsampledWindowedData[j]=ifo[ifonr]->FTin[j];    
  }
  
  
  
  // Allocate memory for Fourier-transform output:
  ifo[ifonr]->FTout = fftw_malloc(sizeof(fftw_complex) * (ifo[ifonr]->FTsize));  
  
  // Contruct a FFTW plan:
  ifo[ifonr]->FTplan = fftw_plan_dft_r2c_1d(N, ifo[ifonr]->FTin, ifo[ifonr]->FTout, FFTW_ESTIMATE);
  //ifo[ifonr]->FTplan = fftw_plan_dft_r2c_1d(N, ifo[ifonr]->FTin, ifo[ifonr]->FTout, FFTW_MEASURE);  //This must be done before initialisation of FTin and could optimise the FFT
  
  // Compute the FFT:
  if(intscrout==1) printf(" | performing data Fourier transform (%.1f s at %d Hz)... ",delta, ifo[ifonr]->samplerate);
  fftw_execute(ifo[ifonr]->FTplan);
  if(ifo[ifonr]->FTout == NULL){
    printf("\n\n   ERROR performing Fourier transform: %s, aborting.\n\n\n",filenames);
    exit(1);
  }
  else if(intscrout==1) printf("ok.\n");
  
  
  // Normalise transform (divide by sampling rate, see Mark's code, line 184):
  for(j=0; j<ifo[ifonr]->FTsize; ++j) ifo[ifonr]->FTout[j] /= (double)ifo[ifonr]->samplerate;
  
  // Copy to 'raw_dataTrafo':
  ifo[ifonr]->raw_dataTrafo = fftw_malloc(sizeof(fftw_complex) * (ifo[ifonr]->FTsize));  
  for(j=0; j<ifo[ifonr]->FTsize; ++j) ifo[ifonr]->raw_dataTrafo[j] = ifo[ifonr]->FTout[j];
  
}
// End dataFT()









void noisePSDestimate(struct interferometer *ifo, struct runPar run)
// Returns a (smoothed) estimate of the log- Power Spectral Density. 
// Data is split into K segments of M seconds,
// and K-1 overlapping segments of length 2M are eventually
// windowed and transformed.
  
{
  struct FrFile  *iFile=NULL;  // Frame File
  struct FrVect   *vect=NULL;
  double           *raw=NULL;  // FFTW vectors etc.               
  double            *in=NULL;
  fftw_complex     *out=NULL;
  fftw_plan           FTplan;
  double           *PSD=NULL;  // vector containing PSD
  double          *sPSD=NULL;  // vector containing smoothed PSD
  double           *win=NULL;  // window
  
  //Use K segments of Mseconds, hence a total of Nseconds of data for the PSD estimation:
  int K = run.PSDsegmentNumber;            // Number of data segments
  double Mseconds = run.PSDsegmentLength;  // Number of seconds in each data segment
  double Nseconds = Mseconds * (double)K;  // Total number of seconds to estimate the PSD for.  Default is 256 seconds
  if(Nseconds < 127.0) printf(" ***  Warning: you're using only %5.1f seconds of data for the PSD estimation;  ~256s is recommended  *** \n",Nseconds);
  if(Nseconds > 1025.0) printf(" ***  Warning: you're using %6.1f seconds of data for the PSD estimation;  ~256s is recommended  *** \n",Nseconds);
  
  
  double wss=0.0,log2=log(2.0);  // squared & summed window coefficients  etc.
  int     i, j, M, N, dummyN;
  int             samplerate;
  int           lower, upper;  // indices of lower & upper frequency bounds in FT vector
  double             nyquist;  // the critical nyquist frequency
  int          smoothrange=0;  // number of samples to left & right that are averaged (=0!!)
  int               PSDrange;
  double                 sum;
  int                 FTsize;
  double         *filtercoef;
  int                  ncoef; 
  char    filenames[2000];
  long             filestart;
  int            filecount=0;
  
  
  // Starting time of first(!) frame file to be read:
  filestart = (((ifo->noiseGPSstart-ifo->noisefileoffset) / ifo->noisefilesize) * ifo->noisefilesize) + ifo->noisefileoffset;
  // Assemble the filename character string:
  while (((double)filestart) < (((double)ifo->noiseGPSstart)+Nseconds)){
    if(filecount == 0) // Fill in filename for first file:
      sprintf(filenames,"%s/%s%ld%s",ifo->noisefilepath,ifo->noisefileprefix,(long)filestart,ifo->noisefilesuffix);
    else // Append filename for following files:
      sprintf(filenames,"%s %s/%s%ld%s",filenames,ifo->noisefilepath,ifo->noisefileprefix,(long)filestart,ifo->noisefilesuffix);
    filestart += ifo->noisefilesize;
    filecount += 1;
  }
  
  // Open Frame file:
  if(intscrout==1) printf(" | opening %d noise data file(s)... \n",filecount);
  if(intscrout==1) printf(" | %s\n",filenames);
  iFile = FrFileINew(filenames);
  if(iFile == NULL) {
    printf("\n\n   ERROR opening noise data file: %s, aborting.\n\n\n",filenames);
    exit(1);
  }
  //else if(intscrout==1) printf("ok.\n");
  if(intscrout==1) printf(" | estimating noise PSD... ");
  // Read first two bits (2M seconds)
  // Access (noise) channel:
  if(ifo->noisedoubleprecision)
    vect = FrFileIGetVectD(iFile, ifo->noisechannel, ((double)ifo->noiseGPSstart), Mseconds*2);
  else
    vect = FrFileIGetVectF(iFile, ifo->noisechannel, ((double)ifo->noiseGPSstart), Mseconds*2);
  if(vect == NULL) {
    printf("\n\n   ERROR reading noise data file: %s, aborting.\n\n\n",filenames);
    exit(1);
  }
  
  N = vect->nData; // Length of filtered & downsampled data (not yet!)
  M = (int)(N/2.0);
  samplerate = (int)(1.0 / (vect->dx[0]) + 0.5);
  // Add 0.5 for correct truncation/rounding
  
  
  // Copy data to vector `raw'
  raw = (double*)malloc(sizeof(double)*N);
  for(i=0; i<N; ++i) 
    raw[i] = vect->dataF[i];
  int screwcount = 0;
  for(i=0; i<N; ++i)
    if(!(raw[i]<HUGE_VAL))
      ++screwcount;
  
  // Downsample (by factor downsamplefactor)  !! changes value of `N' as well !!
  // NINJA:
  if(downsamplefactor!=1){
    filtercoef = filter(&ncoef, samplerate, ifo->highCut);
    in = downsample(raw, &N, filtercoef, ncoef);
    FTsize = (N/2)+1;
    samplerate = (int)((double)samplerate/(double)downsamplefactor);
  } else{
    in = (double*) fftw_malloc(sizeof(double)*N);
    for (i=0; i<N; ++i)
      in[i] = raw[i];
    FTsize = (N/2)+1;    
  }
  nyquist      = ((double)samplerate)/2.0;
  lower        = (int)(floor((ifo->lowCut/nyquist)*(FTsize-1)));
  upper        = (int)(ceil((ifo->highCut/nyquist)*(FTsize-1)));
  ifo->PSDsize = upper-lower;
  PSDrange     = ifo->PSDsize+2*smoothrange;
  
  
  
  // *** Window the data
  // Allocate memory for window:
  win = (double *) malloc(sizeof(double) * N);
  
  // Compute Hahn windowing coefficients:
  for(i=0; i<N; ++i) {
    win[i] = hann(i, N);
    wss += win[i]*win[i];
  }
  // Normalise window coefs as in Mark's/Nelson's code (see line 695):
  for(i=0; i<N; ++i) {
    win[i] /= sqrt(wss * ((double)(K-1)) * ((double)samplerate));
  }
  
  // Apply window to data:
  for(i=0; i<N; ++i) {
    in[i] *= win[i];
  }
  
  
  
  // Put last M values in first M places for following iterations:
  //  ('vect' vectors in later iterations are only half as long)
  for(i=0; i<M; ++i)
    vect->dataF[i] = vect->dataF[i+M];
  
  // Allocate memory for transform output
  out = fftw_malloc(sizeof(fftw_complex) * FTsize);  
  
  // Contruct a transform plan:
  FTplan = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);
  // ('FFTW_MEASURE' option not appropriate here.)
  
  // Execute Fourier transform:
  fftw_execute(FTplan);
  if(out == NULL){
    printf("\n\n   ERROR performing noise Fourier transform: %s, aborting.\n\n\n",filenames);
    exit(1);
  }
  fftw_free(in); 
  
  
  // Allocate & initialize PSD vector:
  PSD = (double *) malloc(PSDrange*sizeof(double));
  for(i=0; i<PSDrange; ++i)
    PSD[i] = pow(cabs(out[(lower-smoothrange)+i]), 2.0);
  
  
  // Read segments 3 to K:
  for(j=3; j<=K; ++j){ 
    // Copy first half of data from previous iteration (M seconds):
    for(i=0; i<M; ++i) 
      raw[i] = vect->dataF[i];
    // Read 2nd half of data (again, M seconds):
    FrVectFree(vect);
    if(ifo->noisedoubleprecision) {
      vect = FrFileIGetVectD(iFile, ifo->noisechannel, ((double)ifo->noiseGPSstart)+((double)(j-1))*Mseconds, Mseconds);
    } else {
      vect = FrFileIGetVectF(iFile, ifo->noisechannel, ((double)ifo->noiseGPSstart)+((double)(j-1))*Mseconds, Mseconds);
    }
    if(vect == NULL) printf("\n : error accessing noise channel!\n");
    
    // Copy 2nd half of data:
    for(i=0; i<M; ++i) 
      raw[i+M] = vect->dataF[i];
    
    for(i=0; i<(2*M); ++i)
      if(!(raw[i]<HUGE_VAL))
        ++screwcount;
    
    // Downsample:
    if(downsamplefactor!=1){
    	dummyN = 2*M;
    	in = downsample(raw, &dummyN, filtercoef, ncoef);
    }else{
    	in = (double*) fftw_malloc(sizeof(double)*N);
    	for (i=0; i<N; ++i)
      		in[i] = raw[i];
    }
    
    // Window data:
    for(i=0; i<N; ++i)
      in[i] *= win[i];
    
    // Execute FT:
    fftw_destroy_plan(FTplan); // Previous `in'-vector was freed in the meantime
    FTplan = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);
    fftw_execute(FTplan);
    if(out == NULL){
      printf("\n : error performing (noise) Fourier transform!\n");
    }
    fftw_free(in); 
    
    // Add to PSD vector
    for(i=0; i<PSDrange; ++i) PSD[i] += pow(cabs(out[(lower-smoothrange)+i]),2.0);
  }
  
  FrVectFree(vect);
  FrFileIEnd(iFile);
  fftw_destroy_plan(FTplan);
  fftw_free(out);
  free(win);
  if(downsamplefactor!=1)  free(filtercoef);
  
  
  // 'PSD' now contains the squared FT magnitudes, summed over (K-1) segments
  
  // Normalise & log the summed PSD's:
  for(i=0; i<PSDrange; ++i)
    PSD[i] = log(PSD[i]) + log2;
  
  if(intscrout==1) printf("ok.\n");
  //if(intscrout==1) printf(" | averaged over %d overlapping segments of %1.0fs each (%.0f s total).\n", K-1, Mseconds*2, Nseconds);
  if(screwcount>0){
    printf(" : %d missing data points in NOISE file(s) !!\n",screwcount);
    printf(" : (maybe the precision is incorrect)\n");
  }
  
  // Smooth PSD:
  //sPSD = (double *) malloc((ifo->PSDsize)*sizeof(double));
  //Since ifo->raw_noisePSD = sPSD, ifo->raw_noisePSD[highindex] goes up to ifo->PSDsize+2 (perhaps more?) in interpol_log_noisePSD()
  sPSD = (double *) malloc((ifo->PSDsize+10)*sizeof(double));
  for(i=0;i<ifo->PSDsize+10;i++) sPSD[i] = 0.0;
  for(i=smoothrange; i<(PSDrange-smoothrange); ++i) {
    sum = 0.0;
    for(j=-smoothrange; j<=smoothrange; ++j)
      sum += PSD[i+j];
    sPSD[i-smoothrange] = sum / (2.0*smoothrange+1.0);
  }
  if(smoothrange>0 && intscrout==1) printf(" | and a range of +/- %0.2f Hz.\n", (((double)smoothrange)/((double)FTsize))*nyquist);

  
  // PSD estimation finished
  free(PSD);
  ifo->raw_noisePSD = sPSD;
  free(raw);
}
// End of void noisePSDestimate







double interpol_log_noisePSD(double f, struct interferometer *ifo)
// Returns linearly interpolated (log-) noise PSD.
{
  double dblindex = (((f-ifo->lowCut)/(ifo->highCut-ifo->lowCut)) * (double)ifo->PSDsize);
  int lowindex    = (int)dblindex;   // (truncated!)
  int highindex   = lowindex + 1;
  double weight1  = ((double)highindex) - dblindex;
  double weight2  = dblindex - ((double)lowindex);
  double result   = weight1*ifo->raw_noisePSD[lowindex] + weight2*ifo->raw_noisePSD[highindex];
  return result;
}






void writeDataToFiles(struct interferometer *ifo[], int networksize, int MCMCseed){
  int i, j;
  for(i=0; i<networksize; i++){ 
    // Write windowed, time-domain data (signal + noise) to disc
    char filename[1000]="";
    sprintf(filename, "%s-data.dat.%6.6d", ifo[i]->name, MCMCseed);  // Write in current dir
    FILE *dump = fopen(filename,"w");
    
    // Get true signal parameters and write them to the header
    if(writeSignal==2){
	printParameterHeaderToFile(dump);
    	fprintf(dump, "       GPS time (s)         H(t)\n");
    }
    for(j=0; j<ifo[i]->samplesize; ++j)
      fprintf(dump, "%9.9f %13.6e\n", 
	ifo[i]->FTstart+(((double)j)/((double) (ifo[i]->samplerate))), 
	ifo[i]->rawDownsampledWindowedData[j]);
    fclose(dump);
    if(intscrout) printf(" : (data written to file)\n");
  
    // Write data FFT to disc (i.e., FFT or amplitude spectrum of signal+noise)
    double f;
    sprintf(filename, "%s-dataFFT.dat.%6.6d", ifo[i]->name, MCMCseed);  // Write in current dir
    FILE *dump1 = fopen(filename,"w");
    
    // Get true signal parameters and write them to the header
    if(writeSignal==2){
    	printParameterHeaderToFile(dump1);
    	fprintf(dump1, "       f (Hz)    real(H(f))    imag(H(f))\n");
    }
    // Loop over the Fourier frequencies 
    for(j=0; j<ifo[i]->indexRange; ++j){
      f = ((double)(j+ifo[i]->lowIndex))/((double) ifo[i]->deltaFT);
      //if(f>0.9*ifo[i]->lowCut) 
      fprintf(dump1, "%13.6e %13.6e %13.6e\n", 
	f, creal(ifo[i]->dataTrafo[j]), cimag(ifo[i]->dataTrafo[j]) );  
	//Save the real and imaginary parts of the data FFT
	//Note that data FFT is already properly normalized in dataFT()
    }
    fclose(dump1);
    if(intscrout) printf(" : (data FFT written to file)\n");
  }
}





void writeNoiseToFiles(struct interferometer *ifo[], int networksize, int MCMCseed){
  int i, j;
  for(i=0; i<networksize; i++){ 
    //Write noise ASD  (Square root of the estimated noise PSD (i.e., no injected$
    char filename[1000]="";
    sprintf(filename, "%s-noiseASD.dat.%6.6d", ifo[i]->name, MCMCseed);  // Write in current dir
    FILE *dump = fopen(filename,"w");

    if(writeSignal==2){
	printParameterHeaderToFile(dump);
    	fprintf(dump, "       f (Hz)          H(f)\n");
    }

    // Loop over the Fourier frequencies within operational range (some 40-1500 H$
    double f=0.0;
    for(j=0; j<ifo[i]->indexRange; ++j){
      f = ((double)(j+ifo[i]->lowIndex))/((double) ifo[i]->deltaFT);
      fprintf(dump, "%13.6f %13.6e %13.6e\n",f, sqrt(ifo[i]->noisePSD[j]), 0.0);
    }
    fclose(dump);
    if(intscrout) printf(" : (noise ASD written to file)\n");
  }
}






void writeSignalsToFiles(struct interferometer *ifo[], int networksize, struct runPar run){
  int i, j;
  //Set local values in parameter struct (needed for template computation)
  struct parset par;
  getInjectionParameters(&par, run.nInjectPar, run.injParVal);
  par.loctc    = (double*)calloc(networksize,sizeof(double));
  par.localti  = (double*)calloc(networksize,sizeof(double));
  par.locazi   = (double*)calloc(networksize,sizeof(double));
  par.locpolar = (double*)calloc(networksize,sizeof(double));
  localpar(&par, ifo, networksize);
  
  
  for(i=0; i<networksize; i++){
    double f;
    double complex FFTout;
    
    // Fill `ifo[i]->FTin' with time-domain template:
    template(&par, ifo, i, run.injectionWaveform);
    // And FFT it
    fftw_execute(ifo[i]->FTplan);
    
    // Write signal in time domain:
    char filename[1000]="";
    sprintf(filename, "%s-signal.dat.%6.6d", ifo[i]->name, run.MCMCseed);  // Write in current dir
    FILE *dump = fopen(filename,"w");
    if(writeSignal==2){
      printParameterHeaderToFile(dump);
      fprintf(dump, "       GPS time (s)         H(t)\n");
    }
    for(j=0; j<ifo[i]->samplesize; ++j)
      fprintf(dump, "%9.9f %13.6e\n", 
	      ifo[i]->FTstart+((double)j)/((double)ifo[i]->samplerate),
	      ifo[i]->FTin[j]);
    fclose(dump);             
    if(intscrout) printf(" : (signal written to file)\n");
    
    // Write signal FFT to disc (i.e., amplitude spectrum of signal w/o noise):
    sprintf(filename, "%s-signalFFT.dat.%6.6d", ifo[i]->name, run.MCMCseed);  // Write in current dir
    FILE *dump2 = fopen(filename,"w");
    if(writeSignal==2){
      printParameterHeaderToFile(dump2);
      fprintf(dump2, "       f (Hz)    real(H(f))    imag(H(f))\n");
    }
    // Loop over the Fourier frequencies within operational range
    //(some 40-1500 Hz or similar):
    for(j=ifo[i]->lowIndex; j<=ifo[i]->highIndex; ++j){
      f=((double) j)/((double) ifo[i]->deltaFT);
      FFTout=ifo[i]->FTout[j]/((double)ifo[i]->samplerate);
      fprintf(dump2, "%13.6e %13.6e %13.6e\n",f, creal(FFTout), cimag(FFTout) ); //Save the real and imaginary parts of the signal FFT
    }
    fclose(dump2);
    if(intscrout) printf(" : (signal FFT written to file)\n");
  }
}




void printParameterHeaderToFile(FILE * dump)
{
  struct parset par;
  //getInjectionParameters(&par, 12, run.injParVal);  //CHECK - get nMCMCpar/nInjectPar and run here somehow
  fprintf(dump,"%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n","m1","m2","mc","eta","tc","dl","lat","lon","phase","spin","kappa","thJ0","phJ0","alpha");
  fprintf(dump,"%12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
	  //par.m1,par.m2,par.mc,par.eta,par.tc,exp(par.logdl),asin(par.sinlati)*r2d,par.longi*r2d,par.phase,par.spin,par.kappa,par.sinthJ0,par.phiJ0,par.alpha);
	  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);  //CHECK: replace with par.par[]
  freeparset(&par);
}


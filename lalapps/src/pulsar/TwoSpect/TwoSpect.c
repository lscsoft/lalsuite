
#include <stdio.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/Sequence.h>
#include <lal/RealFFT.h>
#include <lal/Window.h>
#include <lal/LALMalloc.h>
#include <lal/LALRunningMedian.h>
#include <lal/RngMedBias.h>
#include <lal/DetResponse.h>
#include <lal/Date.h>
#include <lal/SFTfileIO.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "cmdline.h"
#include "IHS.h"
#include "candidates.h"
#include "templates.h"
#include "antenna.h"
#include "TwoSpect.h"


//Global variables
candidate *ihsCandidates[10000];
candidate *gaussCandidates1[10000];
candidate *gaussCandidates2[10000];
candidate *gaussCandidates3[10000];
candidate *gaussCandidates4[10000];
candidate *exactCandidates1[10000];
candidate *exactCandidates2[10000];
//signalParamsStruct *params;
inputParamsStruct *inputParams;
REAL4FFTPlan *secondFFTplan;

FILE *LOG;//, *CAND1, *CAND2;//, *DATOUT;

CHAR *earth_ephemeris = NULL, *sun_ephemeris = NULL;



//Main program
int main(int argc, char *argv[])
{

   INT4 ii, jj, kk, ll, cols, numofcandidates, numofcandidates2;
   REAL4 fsig, per0, B, ihsfarthresh;
   char s[20000];

   struct gengetopt_args_info args_info;
   if ( cmdline_parser(argc, argv, &args_info) ) exit(-1);
   
   //Create directory
   if (args_info.outdirectory_given) {
      mkdir(args_info.outdirectory_arg, 0777);
      snprintf(s, 20000, "%s/logfile.txt", args_info.outdirectory_arg);
   } else {
      mkdir("output",0777);
      snprintf(s, 20000, "%s/logfile.txt", "output");
   }
   
   LOG = fopen(s,"w");
   fprintf(LOG,"Starting TwoSpect analysis...\n");
   fprintf(stderr,"Starting TwoSpect analysis...\n");
   
   //params = new_params();
   inputParams = new_inputParams();
   
   /* if (args_info.A_given) params->amplitude = args_info.A_arg;
   else params->amplitude = 0.015;
   if (args_info.f0_given) params->f0 = args_info.f0_arg;
   else params->f0 = 1000;
   if (args_info.P_given) params->period0 = args_info.P_arg;
   else params->period0 = 16*3600;
   if (args_info.df_given) params->moddepth0 = args_info.df_arg;
   else params->moddepth0 = 0.01; */
   if (args_info.Tobs_given) inputParams->Tobs = args_info.Tobs_arg;
   else inputParams->Tobs = 168*3600;
   //if (args_info.fs_given) params->samplerate = 1/args_info.fs_arg;
   //else params->samplerate = 0.05;
   if (args_info.fmin_given) inputParams->fmin = args_info.fmin_arg;
   else inputParams->fmin = 99.99;
   if (args_info.fspan_given) inputParams->fspan = args_info.fspan_arg;
   else inputParams->fspan = 0.02;
   if (args_info.cols_given) cols = args_info.cols_arg;
   else cols = 20;
   //if (args_info.P0_given) params->P0 = args_info.P0_arg;
   //else params->P0 = 0.1;
   if (args_info.t0_given) inputParams->searchstarttime = args_info.t0_arg;
   else inputParams->searchstarttime = 900000000.0;
   if (args_info.ra_given) inputParams->ra = args_info.ra_arg;
   else inputParams->ra = 0.0;
   if (args_info.dec_given) inputParams->dec = args_info.dec_arg;
   else inputParams->dec = 0.0;
   
   //Defaults given or option passed
   inputParams->Tcoh = args_info.Tcoh_arg;
   ihsfarthresh = args_info.ihsfar_arg;
   //params->Pslope = args_info.Pslope_arg;
   //params->Pquad = args_info.Pquad_arg;
   inputParams->blksize = args_info.blksize_arg;
   //params->varySFTs = args_info.varySFTs_arg;
   earth_ephemeris = (CHAR*)XLALMalloc(strlen(args_info.ephemDir_arg)+20);
   sun_ephemeris = (CHAR*)XLALMalloc(strlen(args_info.ephemDir_arg)+20);
   sprintf(earth_ephemeris,"%s/earth05-09.dat",args_info.ephemDir_arg);
   sprintf(sun_ephemeris,"%s/sun05-09.dat",args_info.ephemDir_arg);
   
   //Blocksize should be an odd number
   if (inputParams->blksize % 2 != 1) inputParams->blksize += 1;
   
   //Adjust maximum columns, if necessary
   if (cols > maxModDepth(inputParams->Tobs*0.2, inputParams->Tcoh)*inputParams->Tcoh) {
      cols = (INT4)floor(maxModDepth(inputParams->Tobs*0.2, inputParams->Tcoh)*inputParams->Tcoh);
      fprintf(LOG,"WARNING! Adjusting number of columns due to maximum modulation depth allowed\n");
      fprintf(stderr,"WARNING! Adjusting number of columns due to maximum modulation depth allowed\n");
   }
   if (cols > (INT4)roundf(inputParams->fspan*inputParams->Tcoh)+1) {
      cols = (INT4)roundf(inputParams->fspan*inputParams->Tcoh)+1;
      fprintf(LOG,"WARNING! Adjusting number of columns due to frequency span of band\n");
      fprintf(stderr,"WARNING! Adjusting number of columns due to frequency span of band\n");
   }
   
   //Set existance of SFTs for fake data
   /* INT4Vector *sftexistlist = NULL;
   if (args_info.removeSFTs_arg==1) {
      gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
      sftexistlist = XLALCreateINT4Vector((INT4)floor(2*(params->Tobs/params->Tcoh)-1));
      for (ii=0; ii<sftexistlist->length; ii++) {
         REAL4 randval = gaussRandNum(1.0, rng) + 1.0;
         if (randval < 0.0) sftexistlist->data[ii] = 0;
         else sftexistlist->data[ii] = 1;
      }
      params->SFTexistlist = sftexistlist;
      gsl_rng_free(rng);
   } else {
      sftexistlist = XLALCreateINT4Vector((INT4)floor(2*(params->Tobs/params->Tcoh)-1));
      for (ii=0; ii<sftexistlist->length; ii++) sftexistlist->data[ii] = 1;
      params->SFTexistlist = sftexistlist;
   } */
   
   //Set input parameters this is almost the same as the fake data parameters
   //inputParams->fmin = params->fmin;
   //inputParams->fspan = params->fspan;
   //inputParams->Tobs = params->Tobs;
   //inputParams->Tcoh = params->Tcoh;
   //inputParams->samplerate = params->samplerate;
   //inputParams->searchstarttime = params->searchstarttime;
   //inputParams->ra = params->ra;
   //inputParams->dec = params->dec;
   //inputParams->blksize = params->blksize;
   //inputParams->SFTexistlist = sftexistlist;
   inputParams->dopplerMultiplier = args_info.dopplerMultiplier_arg;
   
   //fprintf(LOG,"Amplitude = %f\n",params->amplitude);
   //fprintf(LOG,"f0 = %f Hz\n",params->f0);
   //fprintf(LOG,"Period = %f sec\n",params->period0);
   //fprintf(LOG,"Mod. depth = %f Hz\n",params->moddepth0);
   fprintf(LOG,"Tobs = %f sec\n",inputParams->Tobs);
   fprintf(LOG,"Tcoh = %f sec\n",inputParams->Tcoh);
   //fprintf(LOG,"dt = %f sec\n",params->samplerate);
   fprintf(LOG,"fmin = %f Hz\n",inputParams->fmin);
   fprintf(LOG,"fspan = %f Hz\n",inputParams->fspan);
   //fprintf(LOG,"Freq. independent noise spectral density = %f\n",params->P0);
   //fprintf(LOG,"Noise spectral density slope coefficient = %f\n",params->Pslope);
   //fprintf(LOG,"Noise spectral density quadratic coefficient = %f\n",params->Pquad);
   fprintf(LOG,"Running median blocksize = %d\n",inputParams->blksize);
   //fprintf(stderr,"Amplitude = %f\n",params->amplitude);
   //fprintf(stderr,"f0 = %f Hz\n",params->f0);
   //fprintf(stderr,"Period = %f sec\n",params->period0);
   //fprintf(stderr,"Mod. depth = %f Hz\n",params->moddepth0);
   fprintf(stderr,"Tobs = %f sec\n",inputParams->Tobs);
   fprintf(stderr,"Tcoh = %f sec\n",inputParams->Tcoh);
   //fprintf(stderr,"dt = %f sec\n",params->samplerate);
   fprintf(stderr,"fmin = %f Hz\n",inputParams->fmin);
   fprintf(stderr,"fspan = %f Hz\n",inputParams->fspan);
   //fprintf(stderr,"Freq. independent noise spectral density = %f\n",params->P0);
   //fprintf(stderr,"Noise spectral density slope coefficient = %f\n",params->Pslope);
   //fprintf(stderr,"Noise spectral density quadratic coefficient = %f\n",params->Pquad);
   fprintf(stderr,"Running median blocksize = %d\n",inputParams->blksize);
   
   //Allocate memory for ffdata structure
   ffdataStruct *ffdata = new_ffdata(inputParams, 0);
   
   //Second fft plan, only need to make this once for all the exact templates
   secondFFTplan = XLALCreateForwardREAL4FFTPlan((UINT4)floor(2*(inputParams->Tobs/inputParams->Tcoh)-1), 0);
   
   //Create fake data
   fprintf(LOG,"Input SFT data... ");
   fprintf(stderr,"Input SFT data... ");
   //makeFakeBinarySig(ffdata, params);
   makeFakeBinarySigFF(ffdata, inputParams);
   fprintf(LOG,"done\n");
   fprintf(stderr,"done\n");
   
   //Calculate the mean value of the 2nd FFT plane
   REAL4 meanval = calcMean(ffdata->ffdata);
   REAL4 stddev = calcStddev(ffdata->ffdata);
   fprintf(stderr,"2nd PSD mean value = %g\n",meanval);
   fprintf(stderr,"2nd PSD std dev value = %g\n",stddev);
   fprintf(LOG,"2nd PSD mean value = %g\n",meanval);
   fprintf(LOG,"2nd PSD std dev value = %g\n",stddev);
   meanval = calcRms(ffdata->antweights);
   fprintf(stderr,"RMS value of antenna weights = %g\n",meanval);
   
   
   //Average noise floor of FF plane for each 1st FFT frequency bin
   REAL4Vector *aveNoise = ffPlaneNoise(inputParams, ffdata->backgrnd, ffdata->antweights);
   
   //Find the FAR of IHS sum
   ihsfarStruct *ihsfarstruct = new_ihsfarStruct(cols);
   genIhsFar(ihsfarstruct, ffdata, cols, ihsfarthresh);
   fprintf(LOG,"Maximum column width to be searched = %d, at threshold = %f\n",cols,ihsfarthresh);
   fprintf(stderr,"Maximum column width to be searched = %d, at threshold = %f\n",cols,ihsfarthresh);
   
   //IHS FOM FAR (allows a relative offset of +/- 1 bin between maximum values
   REAL4 ihsfomfar = 6.0;
   fprintf(LOG,"IHS FOM FAR = %f\n",ihsfomfar);
   fprintf(stderr,"IHS FOM FAR = %f\n",ihsfomfar);
   
   
/****** Incoherent Harmonic Sum step starts here *****/
   //Allocate memory for the new ihsMaximaStruct
   ihsMaximaStruct *ihsmaxima = new_ihsMaxima(ffdata, cols);
   
   //Run the IHS algorithm on the data
   fprintf(LOG,"Running IHS search on data... ");
   fprintf(stderr,"Running IHS search on data... ");
   runIHS(ihsmaxima, ffdata, cols);
   fprintf(LOG,"done\n");
   fprintf(stderr,"done\n");
   
   numofcandidates = 0;
   REAL4Vector *ihss, *noiseinrange, *ihsexpect;
   INT4Vector *locs;
   ffdataStruct *ffdata_temp;
   ihsMaximaStruct *ihsmaxima_temp;
   fprintf(LOG,"Checking IHS values for candidates... ");
   fprintf(stderr,"Checking IHS values for candidates... ");
   //Check the IHS values against the FAR, checking for >1 column widths
   for (ii=1; ii<ihsfarstruct->ihsfar->length; ii++) {
      ihss = XLALCreateREAL4Vector((UINT4)(ii+1));
      locs = XLALCreateINT4Vector((UINT4)(ii+1));
      noiseinrange = XLALCreateREAL4Vector((UINT4)(ii+1));
      ihsexpect = XLALCreateREAL4Vector((UINT4)(ii+1));
      for (jj=0; jj<ffdata->f->length-ii; jj++) {
      
         //Noise in the range of the columns, mean and rms values for IHS
         for (kk=0; kk<=ii; kk++) {
            noiseinrange->data[kk] = aveNoise->data[jj + kk];
            ihsexpect->data[kk] = aveNoise->data[jj + kk]*ihsfarstruct->ihsdistMean->data[ii];
         }
         REAL4 meanNoise = calcMean(noiseinrange);
         REAL4 rmsNoise = calcRms(noiseinrange);
         
         //Check the IHS sum against the FAR (scaling FAR with mean of the noise in the range of columns)
         if (ihsmaxima->maxima->data[ii*ihsmaxima->locations->length + jj] > 
            ihsfarstruct->ihsfar->data[ii]*meanNoise) {
         
            //Load temporary vectors for determining the FOM
            for (kk=0; kk<=ii; kk++) {
               ihss->data[kk] = ihsmaxima->maxima->data[jj + kk];
               locs->data[kk] = ihsmaxima->locations->data[jj + kk];
            }
            
            //Compute the IHS FOM
            REAL4 fom = ihsFOM(ihss, locs, ihsexpect);
         
            //Check the IHS FOM against the FAR
            if  (fom<=ihsfomfar) {
               
               //Candidate frequency
               fsig = inputParams->fmin + (0.5*ii + (jj))/inputParams->Tcoh;
               
               //Allocate memory for temporary ffdata structure
               ffdata_temp = new_ffdata(inputParams, 1);
            
               //Create vector to hold ffdata with second frequency >5/Tobs to make sure
               //candidates have >5 periods in the observation time
               for (kk=0; kk<ffdata_temp->f->length; kk++) {
                  for (ll=5; ll<ffdata->fpr->length; ll++) {
                     ffdata_temp->ffdata->data[kk*(ffdata->fpr->length-5) + (ll-5)] = 
                        ffdata->ffdata->data[kk*ffdata->fpr->length + ll];
                  }
               }
            
               //Run the IHS maxima algorithm again, so we need memory allocated
               ihsmaxima_temp = new_ihsMaxima(ffdata_temp, ii);
               runIHS(ihsmaxima_temp, ffdata_temp, ii);
               
               //Determine locations. Load temporary vectors of IHS values and locations.
               //Then compute the best location
               REAL4 loc = 0.0;
               for (kk=0; kk<=ii; kk++) {
                  ihss->data[kk] = ihsmaxima_temp->maxima->data[jj + kk];
                  locs->data[kk] = ihsmaxima_temp->locations->data[jj + kk] + 5;
               }
               loc = ihsLoc(ihss, locs, ihsexpect);
               
               //Candidate modulation depth
               B = 0.5*ii/inputParams->Tcoh;
               
               //If the location is not zero, then let's load it as a candidate!
               if (loc!=0.0) {
                  per0 = inputParams->Tobs/loc;
                  
                  REAL4 ihs_sum = ihsmaxima->maxima->data[ii*ihsmaxima->locations->length + jj];
                  //REAL4 ihsSnr = ihs_sum/sqrt(ii+1)/ihsStd;
                  REAL4 ihsSnr = (ihs_sum - meanNoise*ihsfarstruct->ihsdistMean->data[ii])/(rmsNoise*
                     ihsfarstruct->ihsdistSigma->data[ii]);
                  ihsCandidates[numofcandidates] = new_candidate();
                  loadCandidateData(ihsCandidates[numofcandidates], fsig, per0, B, 
                     inputParams->Tobs, inputParams->Tcoh, inputParams->fmin, inputParams->fspan, ihs_sum, ihsSnr);
                  numofcandidates++;
               }
               
               //Destroy unneeded variables
               free_ffdata(ffdata_temp);
               free_ihsMaxima(ihsmaxima_temp);
               ihsmaxima_temp = NULL;
               ffdata_temp = NULL;
            }
         }
      }
      XLALDestroyREAL4Vector(ihss);
      XLALDestroyREAL4Vector(ihsexpect);
      XLALDestroyREAL4Vector(noiseinrange);
      XLALDestroyINT4Vector(locs);
      ihss = NULL;
      locs = NULL;
      ihsexpect = NULL;
      noiseinrange = NULL;
   }
   fprintf(LOG,"done\n");
   fprintf(stderr,"done\n");
   fprintf(LOG,"Candidates found in IHS step = %d\n",numofcandidates);
   fprintf(stderr,"Candidates found in IHS step = %d\n",numofcandidates);
/***** End of Incoherent Harmonic Sum step *****/
   
   
   
/***** Run Gaussian template test on IHS candidates *****/
   farStruct *farval = NULL;
   fprintf(LOG,"Starting Gaussian template search...\n");
   fprintf(stderr,"Starting Gaussian template search...\n");
   numofcandidates2 = 0;
   for (ii=0; ii<numofcandidates; ii++) {
      if (ihsCandidates[ii]->fsig-ihsCandidates[ii]->moddepth-2/inputParams->Tcoh > inputParams->fmin && 
            ihsCandidates[ii]->fsig+ihsCandidates[ii]->moddepth+2/inputParams->Tcoh < inputParams->fmin+inputParams->fspan && 
            ihsCandidates[ii]->moddepth < maxModDepth(ihsCandidates[ii]->period,inputParams->Tcoh)) {
         
         //Allocate memory for template
         ffdataStruct *template = new_ffdata(inputParams, 0);
         
         //Make a Gaussian train template
         makeTemplateGaussians(template, ihsCandidates[ii]);
         
         //Find the top 50 bins in the template
         topbinsStruct *topbinsstruct = new_topbinsStruct(50);
         topbins(topbinsstruct, template, 50);
         
         //Estimate the FAR for these bin weights
         farval = new_farStruct();
         estimateFAR(farval, template->ffdata, topbinsstruct, 0.001, aveNoise);
         
         //Caclulate R
         REAL4 R = calculateR(ffdata->ffdata, template->ffdata, topbinsstruct, aveNoise);
         
         //Destroy unneeded vectors
         free_topbinsStruct(topbinsstruct);
         topbinsstruct = NULL;
         free_ffdata(template);
         template = NULL;
         
         //Log the candidate if R exceeds the FAR or check other possibilities
         //of different periods. Use same farval.far because the FAR is ~independent
         //on period.
         REAL4 Rfirst = R;
         REAL4 Rtemp = R;
         REAL4 bestPeriod = ihsCandidates[ii]->period;
         if (Rfirst > farval->far) {
            gaussCandidates1[numofcandidates2] = new_candidate();
            loadCandidateData(gaussCandidates1[numofcandidates2], ihsCandidates[ii]->fsig, 
               ihsCandidates[ii]->period, ihsCandidates[ii]->moddepth, inputParams->Tobs, inputParams->Tcoh, 
               inputParams->fmin, inputParams->fspan, Rfirst, (Rfirst-farval->distMean)/farval->distSigma);
            numofcandidates2++;
         } else {
            if (ihsCandidates[ii]->period*0.5 > minPeriod(ihsCandidates[ii]->moddepth, 
               inputParams->Tcoh)) {
               ihsCandidates[ii]->period *= 0.5;
               template = new_ffdata(inputParams, 0);
               makeTemplateGaussians(template, ihsCandidates[ii]);
               topbinsstruct = new_topbinsStruct(50);
               topbins(topbinsstruct, template, 50);
               R = calculateR(ffdata->ffdata, template->ffdata, topbinsstruct, aveNoise);
               if (R > Rtemp) {
                  bestPeriod = ihsCandidates[ii]->period;
                  Rtemp = R;
               }
               free_topbinsStruct(topbinsstruct);
               topbinsstruct = NULL;
               free_ffdata(template);
               template = NULL;
               ihsCandidates[ii]->period *= 2.0;
            }
            if (ihsCandidates[ii]->period*2.0/3.0 > minPeriod(ihsCandidates[ii]->moddepth, 
               inputParams->Tcoh)) {
               ihsCandidates[ii]->period *= 2.0/3.0;
               template = new_ffdata(inputParams, 0);
               makeTemplateGaussians(template, ihsCandidates[ii]);
               topbinsstruct = new_topbinsStruct(50);
               topbins(topbinsstruct, template, 50);
               R = calculateR(ffdata->ffdata, template->ffdata, topbinsstruct, aveNoise);
               if (R > Rtemp) {
                  bestPeriod = ihsCandidates[ii]->period;
                  Rtemp = R;
               }
               free_topbinsStruct(topbinsstruct);
               topbinsstruct = NULL;
               free_ffdata(template);
               template = NULL;
               ihsCandidates[ii]->period *= 1.5;
            }
            if (ihsCandidates[ii]->period*0.75 > minPeriod(ihsCandidates[ii]->moddepth, 
               inputParams->Tcoh)) {
               ihsCandidates[ii]->period *= 0.75;
               template = new_ffdata(inputParams, 0);
               makeTemplateGaussians(template, ihsCandidates[ii]);
               topbinsstruct = new_topbinsStruct(50);
               topbins(topbinsstruct, template, 50);
               R = calculateR(ffdata->ffdata, template->ffdata, topbinsstruct, aveNoise);
               if (R > Rtemp) {
                  bestPeriod = ihsCandidates[ii]->period;
                  Rtemp = R;
               }
               free_topbinsStruct(topbinsstruct);
               topbinsstruct = NULL;
               free_ffdata(template);
               template = NULL;
               ihsCandidates[ii]->period /= 0.75;
            }
            if (ihsCandidates[ii]->Tobs/ihsCandidates[ii]->period*0.75 > 5) {
               ihsCandidates[ii]->period *= 4.0/3.0;
               template = new_ffdata(inputParams, 0);
               makeTemplateGaussians(template, ihsCandidates[ii]);
               topbinsstruct = new_topbinsStruct(50);
               topbins(topbinsstruct, template, 50);
               R = calculateR(ffdata->ffdata, template->ffdata, topbinsstruct, aveNoise);
               if (R > Rtemp) {
                  bestPeriod = ihsCandidates[ii]->period;
                  Rtemp = R;
               }
               free_topbinsStruct(topbinsstruct);
               topbinsstruct = NULL;
               free_ffdata(template);
               template = NULL;
               ihsCandidates[ii]->period *= 0.75;
            }
            if (ihsCandidates[ii]->Tobs/ihsCandidates[ii]->period/1.5 > 5) {
               ihsCandidates[ii]->period *= 1.5;
               template = new_ffdata(inputParams, 0);
               makeTemplateGaussians(template, ihsCandidates[ii]);
               topbinsstruct = new_topbinsStruct(50);
               topbins(topbinsstruct, template, 50);
               R = calculateR(ffdata->ffdata, template->ffdata, topbinsstruct, aveNoise);
               if (R > Rtemp) {
                  bestPeriod = ihsCandidates[ii]->period;
                  Rtemp = R;
               }
               free_topbinsStruct(topbinsstruct);
               topbinsstruct = NULL;
               free_ffdata(template);
               template = NULL;
               ihsCandidates[ii]->period /= 1.5;
            }
            if (ihsCandidates[ii]->Tobs/ihsCandidates[ii]->period*0.5 > 5) {
               ihsCandidates[ii]->period *= 2.0;
               template = new_ffdata(inputParams, 0);
               makeTemplateGaussians(template, ihsCandidates[ii]);
               topbinsstruct = new_topbinsStruct(50);
               topbins(topbinsstruct, template, 50);
               R = calculateR(ffdata->ffdata, template->ffdata, topbinsstruct, aveNoise);
               if (R > Rtemp) {
                  bestPeriod = ihsCandidates[ii]->period;
                  Rtemp = R;
               }
               free_topbinsStruct(topbinsstruct);
               topbinsstruct = NULL;
               free_ffdata(template);
               template = NULL;
               ihsCandidates[ii]->period *= 0.5;
            }
            R = Rtemp;
         }
         
         if (R > farval->far) {
            for (jj=0; jj<10; jj++) {
               //REAL4 periodfact = (jj*2+1)*ihsCandidates[ii]->period/(jj+1)*0.5;
               REAL4 periodfact = ihsCandidates[ii]->period/(powf(2,jj+1));
               if ( (ihsCandidates[ii]->period - periodfact) > 
                  minPeriod(ihsCandidates[ii]->moddepth, inputParams->Tcoh)) {
                  ihsCandidates[ii]->period -= periodfact;
                  template = new_ffdata(inputParams, 0);
                  makeTemplateGaussians(template, ihsCandidates[ii]);
                  topbinsstruct = new_topbinsStruct(50);
                  topbins(topbinsstruct, template, 50);
                  R = calculateR(ffdata->ffdata, template->ffdata, topbinsstruct, aveNoise);
                  if (R > Rtemp) {
                     bestPeriod = ihsCandidates[ii]->period;
                     Rtemp = R;
                  }
                  free_topbinsStruct(topbinsstruct);
                  topbinsstruct = NULL;
                  free_ffdata(template);
                  template = NULL;
                  ihsCandidates[ii]->period += periodfact;
               }
               periodfact = ihsCandidates[ii]->period*0.5;
               if ( inputParams->Tobs/(ihsCandidates[ii]->period + (jj+1)*periodfact) > 5) {
                  ihsCandidates[ii]->period += (jj+1)*periodfact;
                  template = new_ffdata(inputParams, 0);
                  makeTemplateGaussians(template, ihsCandidates[ii]);
                  topbinsstruct = new_topbinsStruct(50);
                  topbins(topbinsstruct, template, 50);
                  R = calculateR(ffdata->ffdata, template->ffdata, topbinsstruct, aveNoise);
                  if (R > Rtemp) {
                     bestPeriod = ihsCandidates[ii]->period;
                     Rtemp = R;
                  }
                  free_topbinsStruct(topbinsstruct);
                  topbinsstruct = NULL;
                  free_ffdata(template);
                  template = NULL;
                  ihsCandidates[ii]->period -= (jj+1)*periodfact;
               }
            }
            R = Rtemp;
            
            if (Rfirst<R && Rfirst>farval->far) {
               free_candidate(gaussCandidates1[numofcandidates2-1]);
               gaussCandidates1[numofcandidates2-1] = NULL;
               gaussCandidates1[numofcandidates2-1] = new_candidate();
               loadCandidateData(gaussCandidates1[numofcandidates2-1], ihsCandidates[ii]->fsig, 
                  bestPeriod, ihsCandidates[ii]->moddepth, inputParams->Tobs, inputParams->Tcoh, 
                  inputParams->fmin, inputParams->fspan, R, (R-farval->distMean)/farval->distSigma);
            } else if (Rfirst<R && Rfirst<=farval->far) {
               gaussCandidates1[numofcandidates2] = new_candidate();
               loadCandidateData(gaussCandidates1[numofcandidates2], ihsCandidates[ii]->fsig, 
                  bestPeriod, ihsCandidates[ii]->moddepth, inputParams->Tobs, inputParams->Tcoh, 
                  inputParams->fmin, inputParams->fspan, R, (R-farval->distMean)/farval->distSigma);
               numofcandidates2++;
            }
         }
         
         free_farStruct(farval);
         farval = NULL;
      }
   }
   numofcandidates = numofcandidates2;
   numofcandidates2 = 0;
   fprintf(LOG,"Initial stage done with candidates = %d\n",numofcandidates);
   fprintf(stderr,"Initial stage done with candidates = %d\n",numofcandidates);
   
   //Cluster candidates algorithm
   fprintf(LOG,"Starting to cluster...\n");
   fprintf(stderr,"Starting to cluster...\n");
   clusterCandidates(gaussCandidates2, gaussCandidates1, ffdata, inputParams, aveNoise, numofcandidates, 0);
   numofcandidates = 0;
   for (ii=0; ii<10000; ii++) {
      if (gaussCandidates2[ii]!=NULL) numofcandidates++;
      else ii = 9999;
   }
   fprintf(LOG,"Clustering done with candidates = %d\n",numofcandidates);
   fprintf(stderr,"Clustering done with candidates = %d\n",numofcandidates);
   
   //Perform detailed search over the clustered candidates
   fprintf(LOG,"Starting detailed search using Gaussian train templates... ");
   fprintf(stderr,"Starting detailed search using Gaussian train templates... ");
   REAL4 quadparam = 2.4e-3;
   REAL4 linparam = 4.1e-3;
   for (ii=0; ii<numofcandidates; ii++) {
   
      REAL4Vector *trialf, *trialb, *trialp;
      REAL4 minf, maxf, minb, maxb;
      UINT4 numf, numb, nump;
      
      //Set up parameters of modulation depth search
      minb = gaussCandidates2[ii]->moddepth-1/inputParams->Tcoh;
      maxb = gaussCandidates2[ii]->moddepth+1/inputParams->Tcoh;
      if (minb<0.25/inputParams->Tcoh) minb = 0.5/inputParams->Tcoh;
      numb = (UINT4)roundf(2*(maxb-minb)*inputParams->Tcoh)+1;
      trialb = XLALCreateREAL4Vector(numb);
      for (jj=0; jj<(INT4)numb; jj++) trialb->data[jj] = minb + 0.5*jj/inputParams->Tcoh;
      
      //Set up parameters of signal frequency search
      minf = gaussCandidates2[ii]->fsig-2/inputParams->Tcoh;
      maxf = gaussCandidates2[ii]->fsig+2/inputParams->Tcoh;
      if (minf<inputParams->fmin) minf = inputParams->fmin;
      if (maxf>inputParams->fmin+inputParams->fspan) maxf = inputParams->fmin+inputParams->fspan;
      numf = (UINT4)roundf(2*(maxf-minf)*inputParams->Tcoh)+1;
      trialf = XLALCreateREAL4Vector(numf);
      for (jj=0; jj<(INT4)numf; jj++) trialf->data[jj] = minf + 0.5*jj/inputParams->Tcoh;
      
      //Search over 5 different periods
      nump = 5;
      trialp = XLALCreateREAL4Vector(nump);
      
      //Now search over the parameter space. Frequency, then modulation depth, then period
      REAL4 bestSNR = 0;
      REAL4 bestf, bestp, bestdf, bestR;
      for (jj=0; jj<trialf->length; jj++) {
         for (kk=0; kk<trialb->length; kk++) {
            
            //Start with period of the first guess, then determine nearest neighbor from the
            //modulation depth amplitude to find the other period guesses. These parameters 
            //are determined from simulation to scale the N.N. distance w.r.t. mod. depth with
            //20% mismatch parameter
            trialp->data[2] = gaussCandidates2[ii]->period;
            for (ll=0; ll<2; ll++) {
               REAL4 nnp = (quadparam*(trialp->data[2+ll]/3600)*(trialp->data[2+ll]/3600)+
                  linparam*(trialp->data[2+ll]/3600))*sqrt(3e-3/trialb->data[kk])/
                  powf(2.0,(inputParams->Tobs-7.0*24.0*3600.0)/(7.0*24.0*3600.0))*3600.0;
               trialp->data[2+ll+1] = trialp->data[2+ll] + nnp;
               nnp = (quadparam*(trialp->data[2-ll]/3600)*(trialp->data[2-ll]/3600)+
                  linparam*(trialp->data[2-ll]/3600))*sqrt(3e-3/trialb->data[kk])/
                  powf(2.0,(inputParams->Tobs-7.0*24.0*3600.0)/(7.0*24.0*3600.0))*3600;
               trialp->data[2-(ll+1)] = trialp->data[2-ll] - nnp;
            }
            
            //Take the mean period and compute a template/FAR pair.
            REAL4 tempP = calcMean(trialp);
            candidate *cand = new_candidate();
            loadCandidateData(cand, trialf->data[jj], tempP, trialb->data[kk], 
               inputParams->Tobs, inputParams->Tcoh, inputParams->fmin, inputParams->fspan, 0, 0);
            ffdataStruct *template = new_ffdata(inputParams, 0);
            makeTemplateGaussians(template, cand);
            //INT4Vector *top = topbins(template->ffdata, 50);
            topbinsStruct *topbinsstruct = new_topbinsStruct(50);
            topbins(topbinsstruct, template, 50);
            farStruct *farval = new_farStruct();
            estimateFAR(farval, template->ffdata, topbinsstruct, 0.001, aveNoise);
            //farStruct farval = estimateFAR(template->ffdata, top, 0.01);
            free_topbinsStruct(topbinsstruct);
            topbinsstruct = NULL;
            free_candidate(cand);
            cand = NULL;
            free_ffdata(template);
            template = NULL;
            for (ll=0; ll<trialp->length; ll++) {
               if ( trialf->data[jj]-trialb->data[kk]-2/inputParams->Tcoh > inputParams->fmin && 
                  trialf->data[jj]+trialb->data[kk]+2/inputParams->Tcoh < inputParams->fmin+inputParams->fspan &&
                  trialb->data[kk] < maxModDepth(trialp->data[ll], inputParams->Tcoh) && 
                  trialp->data[ll] > minPeriod(trialb->data[kk], inputParams->Tcoh) && 
                  inputParams->Tobs/trialp->data[ll] > 5.0 ) {
                  cand = new_candidate();
                  loadCandidateData(cand, trialf->data[jj], trialp->data[ll], trialb->data[kk], 
                     inputParams->Tobs, inputParams->Tcoh, inputParams->fmin, inputParams->fspan, 0, 0);
                  template = new_ffdata(inputParams, 0);
                  makeTemplateGaussians(template, cand);
                  topbinsstruct = new_topbinsStruct(50);
                  topbins(topbinsstruct, template, 50);
                  REAL4 R = calculateR(ffdata->ffdata, template->ffdata, topbinsstruct, aveNoise);
                  REAL4 snr = (R - farval->distMean)/farval->distSigma;
                  if (R > farval->far && snr > bestSNR) {
                     bestf = trialf->data[jj];
                     bestp = trialp->data[ll];
                     bestdf = trialb->data[kk];
                     bestR = R;
                     bestSNR = snr;
                  }
                  free_topbinsStruct(topbinsstruct);
                  topbinsstruct = NULL;
                  free_candidate(cand);
                  cand = NULL;
                  free_ffdata(template);
                  template = NULL;
               }
            }
            
            free_farStruct(farval);
            farval = NULL;
         }
      }
      
      gaussCandidates3[numofcandidates2] = new_candidate();
      loadCandidateData(gaussCandidates3[numofcandidates2], bestf, bestp, bestdf, inputParams->Tobs, 
         inputParams->Tcoh, inputParams->fmin, inputParams->fspan, bestR, bestSNR);
      numofcandidates2++;
      
      //free_candidate(gaussCandidates2[ii]);
      //gaussCandidates2[ii] = NULL;
      
      XLALDestroyREAL4Vector(trialf);
      XLALDestroyREAL4Vector(trialb);
      XLALDestroyREAL4Vector(trialp);
      trialf = NULL;
      trialb = NULL;
      trialp = NULL;
   }
   fprintf(LOG,"done\n");
   fprintf(stderr,"done\n");
   
   //Cluster candidates algorithm, second time
   fprintf(LOG,"Starting the second round of clustering...\n");
   fprintf(stderr,"Starting the second round of clustering...\n");
   clusterCandidates(gaussCandidates4, gaussCandidates3, ffdata, inputParams, aveNoise, numofcandidates2, 0);
   numofcandidates = 0;
   numofcandidates2 = 0;
   for (ii=0; ii<10000; ii++) {
      if (gaussCandidates4[ii]!=NULL) numofcandidates++;
      else ii = 9999;
   }
   fprintf(LOG,"Clustering done with candidates = %d\n",numofcandidates);
   fprintf(stderr,"Clustering done with candidates = %d\n",numofcandidates);
/***** End of Gaussian template step *****/
   
   
/***** Start of exact template step *****/
   //Start with a test of each candidate using the exact template
   fprintf(LOG,"Starting exact template search...\n");
   fprintf(stderr,"Starting exact template search...\n");
   for (ii=0; ii<numofcandidates; ii++) {
      ffdataStruct *template = new_ffdata(inputParams, 0);
      makeTemplate(template, gaussCandidates4[ii], secondFFTplan);
      //INT4Vector *top = topbins(template->ffdata, 50);
      topbinsStruct *topbinsstruct = new_topbinsStruct(50);
      topbins(topbinsstruct, template, 50);
      farStruct *farval = new_farStruct();
      estimateFAR(farval, template->ffdata, topbinsstruct, 0.001, aveNoise);
      //farStruct farval = estimateFAR(template->ffdata, top, 0.01);
      REAL4 R = calculateR(ffdata->ffdata, template->ffdata, topbinsstruct, aveNoise);
      REAL4 SNR = (R - farval->distMean)/farval->distSigma;
      if (R > farval->far) {
         exactCandidates1[numofcandidates2] = new_candidate();
         loadCandidateData(exactCandidates1[numofcandidates2], gaussCandidates4[ii]->fsig, 
            gaussCandidates4[ii]->period, gaussCandidates4[ii]->moddepth, gaussCandidates4[ii]->Tobs, 
            gaussCandidates4[ii]->Tcoh, gaussCandidates4[ii]->fmin, gaussCandidates4[ii]->fspan, 
            R, SNR);
         numofcandidates2++;
      }
      
      free_topbinsStruct(topbinsstruct);
      topbinsstruct = NULL;
      free_ffdata(template);
      template = NULL;
      free_farStruct(farval);
      farval = NULL;
   }
   numofcandidates = numofcandidates2;
   numofcandidates2 = 0;
   fprintf(LOG,"Number of candidates confirmed with exact templates = %d\n",numofcandidates);
   fprintf(stderr,"Number of candidates confirmed with exact templates = %d\n",numofcandidates);
   
   //Perform detailed search using exact templates
   fprintf(LOG,"Starting detailed search with exact templates...\n");
   fprintf(stderr,"Starting detailed search with exact templates...\n");
   for (ii=0; ii<numofcandidates; ii++) {
   
      REAL4Vector *trialf, *trialb, *trialp;
      REAL4 minf, maxf, minb, maxb;
      UINT4 numf, numb, nump;
      
      minb = exactCandidates1[ii]->moddepth-1/inputParams->Tcoh;
      maxb = exactCandidates1[ii]->moddepth+1/inputParams->Tcoh;
      if (minb<0.25/inputParams->Tcoh) minb = 0.5/inputParams->Tcoh;
      numb = (UINT4)roundf(2*(maxb-minb)*inputParams->Tcoh)+1;
      trialb = XLALCreateREAL4Vector(numb);
      for (jj=0; jj<(INT4)numb; jj++) trialb->data[jj] = minb + 0.5*jj/inputParams->Tcoh;

      minf = exactCandidates1[ii]->fsig-2/inputParams->Tcoh;
      maxf = exactCandidates1[ii]->fsig+2/inputParams->Tcoh;
      if (minf<inputParams->fmin) minf = inputParams->fmin;
      if (maxf>inputParams->fmin+inputParams->fspan) maxf = inputParams->fmin+inputParams->fspan;
      numf = (UINT4)roundf(2*(maxf-minf)*inputParams->Tcoh)+1;
      trialf = XLALCreateREAL4Vector(numf);
      for (jj=0; jj<(INT4)numf; jj++) trialf->data[jj] = minf + 0.5*jj/inputParams->Tcoh;
      
      nump = 5;
      trialp = XLALCreateREAL4Vector(nump);
      
      REAL4 bestSNR = 0;
      REAL4 bestf, bestp, bestdf, bestR;
      for (jj=0; jj<trialf->length; jj++) {
         for (kk=0; kk<trialb->length; kk++) {
            trialp->data[2] = exactCandidates1[ii]->period;
            for (ll=0; ll<2; ll++) {
               REAL4 nnp = (quadparam*(trialp->data[2+ll]/3600)*(trialp->data[2+ll]/3600)+
                  linparam*(trialp->data[2+ll]/3600))*sqrt(3e-3/trialb->data[kk])/
                  powf(2.0,(inputParams->Tobs-7.0*24.0*3600.0)/(7.0*24.0*3600.0))*3600;
               trialp->data[2+ll+1] = trialp->data[2+ll] + nnp;
               nnp = (quadparam*(trialp->data[2-ll]/3600)*(trialp->data[2-ll]/3600)+
                  linparam*(trialp->data[2-ll]/3600))*sqrt(3e-3/trialb->data[kk])/
                  powf(2.0,(inputParams->Tobs-7.0*24.0*3600.0)/(7.0*24.0*3600.0))*3600;
               trialp->data[2-(ll+1)] = trialp->data[2-ll] - nnp;
            }
            
            REAL4 tempP = calcMean(trialp);
            candidate *cand = new_candidate();
            loadCandidateData(cand, trialf->data[jj], tempP, trialb->data[kk], 
               inputParams->Tobs, inputParams->Tcoh, inputParams->fmin, inputParams->fspan, 0, 0);
            ffdataStruct *template = new_ffdata(inputParams, 0);
            makeTemplate(template, cand, secondFFTplan);
            //INT4Vector *top = topbins(template->ffdata, 50);
            topbinsStruct *topbinsstruct = new_topbinsStruct(50);
            topbins(topbinsstruct, template, 50);
            //farStruct farval = estimateFAR(template->ffdata, top, 0.01);
            farStruct *farval = new_farStruct();
            estimateFAR(farval, template->ffdata, topbinsstruct, 0.001, aveNoise);
            free_topbinsStruct(topbinsstruct);
            topbinsstruct = NULL;
            free_candidate(cand);
            cand = NULL;
            free_ffdata(template);
            template = NULL;
            for (ll=0; ll<trialp->length; ll++) {
               if ( trialf->data[jj]-trialb->data[kk]-2/inputParams->Tcoh > inputParams->fmin && 
                  trialf->data[jj]+trialb->data[kk]+2/inputParams->Tcoh < inputParams->fmin+inputParams->fspan &&
                  trialb->data[kk]<maxModDepth(trialp->data[ll], inputParams->Tcoh) && 
                  trialp->data[ll] > minPeriod(trialb->data[kk], inputParams->Tcoh) && 
                  inputParams->Tobs/trialp->data[ll]>5 ) {
                  cand = new_candidate();
                  loadCandidateData(cand, trialf->data[jj], trialp->data[ll], trialb->data[kk], 
                     inputParams->Tobs, inputParams->Tcoh, inputParams->fmin, inputParams->fspan, 0, 0);
                  template = new_ffdata(inputParams, 0);
                  makeTemplate(template, cand, secondFFTplan);
                  topbinsstruct = new_topbinsStruct(50);
                  topbins(topbinsstruct, template, 50);
                  REAL4 R = calculateR(ffdata->ffdata, template->ffdata, topbinsstruct, aveNoise);
                  REAL4 SNR = (R - farval->distMean)/farval->distSigma;
                  if (R > farval->far && SNR > bestSNR) {
                     bestf = trialf->data[jj];
                     bestp = trialp->data[ll];
                     bestdf = trialb->data[kk];
                     bestR = R;
                     bestSNR = SNR;
                  }
                  free_topbinsStruct(topbinsstruct);
                  topbinsstruct = NULL;
                  free_candidate(cand);
                  cand = NULL;
                  free_ffdata(template);
                  template = NULL;
               }
            }
            
            free_farStruct(farval);
            farval = NULL;
         }
      }
      
      exactCandidates2[numofcandidates2] = new_candidate();
      loadCandidateData(exactCandidates2[numofcandidates2], bestf, bestp, bestdf, inputParams->Tobs, 
         inputParams->Tcoh, inputParams->fmin, inputParams->fspan, bestR, bestSNR);
      numofcandidates2++;
      
      XLALDestroyREAL4Vector(trialf);
      XLALDestroyREAL4Vector(trialb);
      XLALDestroyREAL4Vector(trialp);
      trialf = NULL;
      trialb = NULL;
      trialp = NULL;
   }
   fprintf(LOG,"Exact step is done with the number of candidates = %d\n",numofcandidates2);
   fprintf(stderr,"Exact step is done with the number of candidates = %d\n",numofcandidates2);
   
   
/***** End of exact template step *****/
   
   
   /* if (args_info.outdirectory_given) {
      snprintf(s, 20000, "%s/candidate_confirmed.txt", args_info.outdirectory_arg);
      CAND1 = fopen(s,"a");
      snprintf(s, 20000, "%s/candidate_possible.txt", args_info.outdirectory_arg);
      CAND2 = fopen(s,"a");
   } else {
      snprintf(s, 20000, "output/candidate_confirmed.txt");
      CAND1 = fopen(s,"a");
      snprintf(s, 20000, "output/candidate_possible.txt");
      CAND2 = fopen(s,"a");
   }
   INT4 num1 = 0;
   INT4 num2 = 0; */
   
   if (numofcandidates2!=0) {
      fprintf(LOG,"\n**Report of candidates:**\n");
      fprintf(stderr,"\n**Report of candidates:**\n");
      
      for (ii=0; ii<numofcandidates2; ii++) {
         fprintf(LOG,"fsig = %f, period = %f, df = %f, R = %f, SNR = %f\n",
            exactCandidates2[ii]->fsig, exactCandidates2[ii]->period, exactCandidates2[ii]->moddepth, 
            exactCandidates2[ii]->stat, exactCandidates2[ii]->snr);
         fprintf(stderr,"fsig = %f, period = %f, df = %f, R = %f, SNR = %f\n",
            exactCandidates2[ii]->fsig, exactCandidates2[ii]->period, exactCandidates2[ii]->moddepth, 
            exactCandidates2[ii]->stat, exactCandidates2[ii]->snr);
         
         
         /* if ( fabs(exactCandidates2[ii]->fsig-params->f0)*params->Tcoh <= 2. && 
            fabs(params->Tobs/exactCandidates2[ii]->period-params->Tobs/params->period0) <= 2. &&  
            fabs(exactCandidates2[ii]->moddepth-params->moddepth0)*params->Tcoh <= 2. ) {
            num2++;
         }
         if ( fabs(exactCandidates2[ii]->fsig-params->f0)*params->Tcoh <= 1. && 
            fabs(params->Tobs/exactCandidates2[ii]->period-params->Tobs/params->period0) <= 1. &&  
            fabs(exactCandidates2[ii]->moddepth-params->moddepth0)*params->Tcoh <= 1. ) {
            num1++;
         } */
      }
      /* fprintf(CAND1,"%d\n",num1);
      fclose(CAND1);
      fprintf(CAND2,"%d\n",num2);
      fclose(CAND2); */

      fprintf(LOG,"\n");
      fprintf(stderr,"\n");
   } else {
      /* fprintf(CAND1,"%d\n",num1);
      fclose(CAND1);
      fprintf(CAND2,"%d\n",num2);
      fclose(CAND2); */
   }
   
   
   
   //Destroy varaibles
   free_ffdata(ffdata);
   free_ihsMaxima(ihsmaxima);
   free_ihsfarStruct(ihsfarstruct);
   //XLALDestroyINT4Vector(sftexistlist);
   //free_params(params);
   free_inputParams(inputParams);
   XLALDestroyREAL4Vector(aveNoise);
   XLALDestroyREAL4FFTPlan(secondFFTplan);
   XLALFree((CHAR*)earth_ephemeris);
   XLALFree((CHAR*)sun_ephemeris);
   cmdline_parser_free(&args_info);
   for (ii=0; ii<10000; ii++) {
      if (ihsCandidates[ii]!=NULL) {
         free_candidate(ihsCandidates[ii]);
         ihsCandidates[ii] = NULL;
      } else {
         ii = 9999;
      }
   }
   for (ii=0; ii<10000; ii++) {
      if (gaussCandidates1[ii]!=NULL) {
         free_candidate(gaussCandidates1[ii]);
         gaussCandidates1[ii] = NULL;
      } else {
         ii = 9999;
      }
   }
   for (ii=0; ii<10000; ii++) {
      if (gaussCandidates2[ii]!=NULL) {
         free_candidate(gaussCandidates2[ii]);
         gaussCandidates2[ii] = NULL;
      } else {
         ii = 9999;
      }
   }
   for (ii=0; ii<10000; ii++) {
      if (gaussCandidates3[ii]!=NULL) {
         free_candidate(gaussCandidates3[ii]);
         gaussCandidates3[ii] = NULL;
      } else {
         ii = 9999;
      }
   }
   for (ii=0; ii<10000; ii++) {
      if (gaussCandidates4[ii]!=NULL) {
         free_candidate(gaussCandidates4[ii]);
         gaussCandidates4[ii] = NULL;
      } else {
         ii = 9999;
      }
   }
   for (ii=0; ii<10000; ii++) {
      if (exactCandidates1[ii]!=NULL) {
         free_candidate(exactCandidates1[ii]);
         exactCandidates1[ii] = NULL;
      } else {
         ii = 9999;
      }
   }
   
   for (ii=0; ii<numofcandidates2; ii++) {
      free_candidate(exactCandidates2[ii]);
      exactCandidates2[ii] = NULL;
   }
   
   fclose(LOG);
   
   return 0;

}
/*** End of main() ***/


//////////////////////////////////////////////////////////////
// Allocate memory for signal parameters struct  -- done
/* signalParamsStruct * new_params(void)
{
   
   signalParamsStruct *param = (signalParamsStruct*)XLALMalloc(sizeof(signalParamsStruct));
   
   return param;
   
} */


//////////////////////////////////////////////////////////////
// Free signal parameters struct memory  -- done
/* void free_params(signalParamsStruct *param)
{

   XLALFree((signalParamsStruct*)param);

} */


inputParamsStruct * new_inputParams(void)
{

   inputParamsStruct *input = (inputParamsStruct*)XLALMalloc(sizeof(inputParamsStruct));
   
   return input;

}



void free_inputParams(inputParamsStruct *input)
{

   XLALFree((inputParamsStruct*)input);

}


//////////////////////////////////////////////////////////////
// Create Gaussian distributed noise vector  -- done
REAL4Vector * gaussRandNumVector(REAL4 sigma, UINT4 length, gsl_rng *ptrToGenerator)
{

   INT4 ii;
   
   REAL4Vector *noise = XLALCreateREAL4Vector(length);
   
   //Create the exponentially distributed noise
   for (ii=0; ii<noise->length; ii++) noise->data[ii] = gaussRandNum(sigma, ptrToGenerator);
   
   return noise;

}


//////////////////////////////////////////////////////////////
// Create a Gaussian distributed noise value  -- done
REAL4 gaussRandNum(REAL4 sigma, gsl_rng *ptrToGenerator)
{

   REAL4 noise = gsl_ran_gaussian(ptrToGenerator, sigma);
   
   return noise;

}


//////////////////////////////////////////////////////////////
// Create exponentially distributed noise vector  -- done
REAL4Vector * expRandNumVector(REAL4 mu, UINT4 length, gsl_rng *ptrToGenerator)
{

   INT4 ii;
   
   REAL4Vector *noise = XLALCreateREAL4Vector(length);
   
   //Create the exponentially distributed noise
   for (ii=0; ii<noise->length; ii++) noise->data[ii] = expRandNum(mu, ptrToGenerator);
   
   return noise;

}


//////////////////////////////////////////////////////////////
// Create a exponentially distributed noise value  -- done
REAL4 expRandNum(REAL4 mu, gsl_rng *ptrToGenerator)
{

   REAL4 noise = gsl_ran_exponential(ptrToGenerator, mu);
   
   return noise;

}



//////////////////////////////////////////////////////////////
// Allocate ffdataStruct vectors  -- done
ffdataStruct * new_ffdata(inputParamsStruct *param, INT4 mode)
{
   
   ffdataStruct *ffdata;
   ffdata = (ffdataStruct*)XLALMalloc(sizeof(ffdataStruct));
   
   UINT4 numfbins = (UINT4)(roundf(param->fspan*param->Tcoh)+1);
   UINT4 numffts = (UINT4)floor(2*(param->Tobs/param->Tcoh)-1);
   
   ffdata->f = XLALCreateREAL4Vector(numfbins);
   if (mode==1) ffdata->fpr = XLALCreateREAL4Vector((UINT4)(floor(numffts*0.5)+1)-5);
   else ffdata->fpr = XLALCreateREAL4Vector((UINT4)floor(numffts*0.5)+1);
   ffdata->ffdata = XLALCreateREAL4Vector(ffdata->f->length * ffdata->fpr->length);
   ffdata->backgrnd = XLALCreateREAL4Vector(ffdata->f->length * numffts);
   ffdata->antweights = XLALCreateREAL4Vector(numffts);
   
   return ffdata;
   
}


//////////////////////////////////////////////////////////////
// Destroy ffdataStruct vectors  -- done
void free_ffdata(ffdataStruct *data)
{

   XLALDestroyREAL4Vector(data->f);
   XLALDestroyREAL4Vector(data->fpr);
   XLALDestroyREAL4Vector(data->ffdata);
   XLALDestroyREAL4Vector(data->backgrnd);
   XLALFree((ffdataStruct*)data);

} 


//////////////////////////////////////////////////////////////
// Read in SFT data to produce a TF vector
REAL4Vector * readInSFTs(inputParamsStruct *input)
{

   INT4 ii, jj;
   LALStatus status;
   status.statusPtr = NULL;
   SFTCatalog *catalog = NULL;
   SFTConstraints *constraints = NULL;
   const CHAR *filenames = "*.sft";
   SFTVector *sfts = NULL;
   
   //Find SFT files
   LALSFTdataFind(&status, &catalog, filenames, constraints);
   
   //Determine band size (remember to get extra bins because of the running median and the bin shifts due to detector velocity)
   REAL8 minfbin = (REAL8)round(input->fmin*input->Tcoh - 0.5*(input->blksize-1) - (input->maxbinshift))/input->Tcoh;
   REAL8 maxfbin = (REAL8)round((input->fmin + input->fspan)*input->Tcoh + 0.5*(input->blksize-1) + (input->maxbinshift))/input->Tcoh;
   
   //Now extract the data
   LALLoadSFTs(&status, &sfts, catalog, minfbin, maxfbin);
   
   //Now put the power data into the TF plane, looping through each SFT
   //If an SFT doesn't exit, fill the TF pixels of the SFT with zeros
   INT4 numffts = (INT4)floor(2*(input->Tobs/input->Tcoh)-1);
   INT4 sftlength = sfts->data->data->length;
   INT4 nonexistantsft = 0;
   REAL4Vector *tfdata = XLALCreateREAL4Vector((UINT4)(numffts*sftlength));
   for (ii=0; ii<numffts; ii++) {
      SFTDescriptor *sftdescription = &(catalog->data[ii - nonexistantsft]); //catalog->data + (UINT4)(ii - nonexistantsft);
      SFTtype *sft = &(sfts->data[ii - nonexistantsft]);
      if (sftdescription->header.epoch.gpsSeconds == (INT4)(ii*0.5*input->Tcoh+input->searchstarttime)) {
         for (jj=0; jj<sftlength; jj++) {
            COMPLEX8 sftcoeff = sft->data->data[jj];
            tfdata->data[ii*sftlength + jj] = (REAL4)(2.0*(sftcoeff.re*sftcoeff.re + sftcoeff.im*sftcoeff.im)); //TODO: check this for consistancy
         }
      } else {
         for (jj=0; jj<sftlength; jj++) {
            tfdata->data[ii*sftlength + jj] = 0.0;
         }
         nonexistantsft++;
      }
   }
   
   LALDestroySFTCatalog(&status, &catalog);
   XLALDestroySFTVector(sfts);
   
   return tfdata;

}


//////////////////////////////////////////////////////////////
// Slide SFT TF data  -- done
REAL4Vector * slideTFdata(inputParamsStruct *input, REAL4Vector *tfdata, INT4Vector *binshifts)
{
   
   INT4 ii, jj;
   INT4 numffts = (INT4)floor(2*(input->Tobs/input->Tcoh)-1);
   REAL4 tempfspan = input->fspan + (input->blksize-1)/input->Tcoh;
   INT4 tempnumfbins = (INT4)(roundf(tempfspan*input->Tcoh)+1);
   
   REAL4Vector *outtfdata = XLALCreateREAL4Vector((UINT4)(numffts*tempnumfbins));
   for (ii=0; ii<numffts; ii++) {
      for (jj=0; jj<tempnumfbins; jj++) outtfdata->data[ii*tempnumfbins + jj] = tfdata->data[ii*(tempnumfbins+2*input->maxbinshift) + jj + input->maxbinshift + binshifts->data[ii]];
   }
   
   return outtfdata;
   
}



//////////////////////////////////////////////////////////////
// Load in SFTs from Makefakedata  -- should make this just a general loading data into FF vector
void makeFakeBinarySigFF(ffdataStruct *ffdata, inputParamsStruct *input)
{
   
   INT4 ii, jj;
   
   INT4 numffts = (INT4)floor(2*(input->Tobs/input->Tcoh)-1);
   INT4 numfbins = (INT4)(roundf(input->fspan*input->Tcoh)+1);
   REAL4 tempfspan = input->fspan + (input->blksize-1)/input->Tcoh;
   INT4 tempnumfbins = (INT4)(roundf(tempfspan*input->Tcoh)+1);
   
   LALDetector det = lalCachedDetectors[LALDetectorIndexLHODIFF]; //H1
   
   //Do ephemeris setup
   EphemerisData *edat = new_Ephemeris(earth_ephemeris, sun_ephemeris);
   
   //Calculate detector velocity w.r.t. a sky location for the SFTs
   REAL4Vector *detVelocity = CompAntennaVelocity(input->ra, input->dec, input->searchstarttime, input->Tcoh, input->Tobs, det, edat);
   INT4Vector *binshifts = CompBinShifts(input->fmin+input->fspan*0.5, detVelocity, input->Tcoh, input->dopplerMultiplier);
   
   //Maximum detector velocity during observation time
   INT4 maxShift = 0;
   for (ii=0; ii<binshifts->length; ii++) if (abs(binshifts->data[ii]) > maxShift) maxShift = abs(binshifts->data[ii]);
   input->maxbinshift = maxShift;
   
   //Compute antenna pattern weights
   REAL4Vector *antenna = CompAntennaPatternWeights(input->ra, input->dec, input->searchstarttime, input->Tcoh, input->Tobs, det);
   memcpy(ffdata->antweights->data, antenna->data, antenna->length*sizeof(*antenna->data));
   //Setting all weights of antenna pattern to 1 -- Remove this before running for real!
   //for (ii=0; ii<antptrnweights->length; ii++) ffdata->antweights->data[ii] = 1.0;
   
   //Create initial SFTs (TF data has more frequency bins because of need to calculate running means)
   REAL4Vector *inSFTs = readInSFTs(input);
   //REAL4Vector *initialTFdata = makeFakeBinarySigTF(in);
   
   //Slide SFTs here
   REAL4Vector *initialTFdata = slideTFdata(input, inSFTs, binshifts);
   REAL4 meanval = calcMean(initialTFdata);
   REAL4 sigval = calcStddev(initialTFdata);
   fprintf(stderr,"Mean value = %g, standard deviation %g\n",meanval,sigval);
   
   //Calculate the running mean values of the SFTs (output here is smaller than initialTFdata). Here,
   //numfbins needs to be the bins you expect to come out of the running means -- the band you are going
   //to search!
   REAL4Vector *background = tfRngMeans(initialTFdata, numffts, numfbins, input->blksize);
   memcpy(ffdata->backgrnd->data, background->data, background->length*sizeof(*background->data));
   
   //This is necessary for the output of data to an ASCII file. Go ahead and remove later
   //char s[20000];
   //snprintf(s, 20000, "output/aftersliding.dat");
   //DATOUT = fopen(s,"w");
   //can be removed
   
   //Need to reduce the original TF data so the weighted TF data can be calculated
   REAL4Vector *usableTFdata = XLALCreateREAL4Vector(ffdata->backgrnd->length);
   for (ii=0; ii<numffts; ii++) {
      for (jj=0; jj<numfbins; jj++) {
         usableTFdata->data[ii*numfbins + jj] = initialTFdata->data[ii*tempnumfbins + jj + (INT4)roundf((input->blksize-1)*0.5)];
         //fprintf(DATOUT,"%g\n",usableTFdata->data[ii*numfbins + jj]);   //remove this
      }
   }
   //fclose(DATOUT);   //remove this
   meanval = calcMean(ffdata->backgrnd);
   sigval = calcStddev(ffdata->backgrnd);
   fprintf(stderr,"Mean value = %g, standard deviation %g\n",meanval,sigval);
   
   //Compute the weighted TF data
   REAL4Vector *weightedTFdata = tfWeightMeanSubtract(usableTFdata, ffdata->backgrnd, ffdata->antweights, input);
   
   //Do the second FFT
   REAL4Vector *secFFTdata = makeSecondFFT(weightedTFdata, input);
   memcpy(ffdata->ffdata->data, secFFTdata->data, secFFTdata->length*sizeof(*secFFTdata->data));
   //for(ii=0; ii<secFFTdata->length; ii++) fprintf(DATOUT,"%g\n",secFFTdata->data[ii]); //remove this
   //fclose(DATOUT); //remove this
   
   //Make the first FFT frequencies
   for (ii=0; ii<ffdata->f->length; ii++) ffdata->f->data[ii] = input->fmin + ii/input->Tcoh;
   
   //Make the second FFT frequencies
   for (ii=0; ii<ffdata->fpr->length; ii++) ffdata->fpr->data[ii] = ii/input->Tobs;
   
   XLALDestroyREAL4Vector(inSFTs);
   XLALDestroyREAL4Vector(antenna);
   XLALDestroyREAL4Vector(initialTFdata);
   XLALDestroyREAL4Vector(background);
   XLALDestroyREAL4Vector(usableTFdata);
   XLALDestroyREAL4Vector(weightedTFdata);
   XLALDestroyREAL4Vector(secFFTdata);
   free_Ephemeris(edat);
   
}


//////////////////////////////////////////////////////////////
// Create fake TF data  -- done
/* REAL4Vector * makeFakeBinarySigTF(signalParamsStruct *in)
{
   
   INT4 numffts, numfbins, tempnumfbins, ii, jj, block;
   REAL4 A, f0, period0, B, Tobs, dT, dt, fmin, fspan, noiseA, periodf, tempfmin, tempfspan;
   
   A = in->amplitude;         //Amplitude of signal
   f0 = in->f0;               //Frequency of signal
   period0 = in->period0;     //Period of signal
   B = in->moddepth0;         //Modulation depth of signal
   Tobs = in->Tobs;           //Observation time
   dT = in->Tcoh;             //Coherence time
   dt = in->samplerate;       //Sampling rate
   fmin = in->fmin;           //Minimum frequency
   fspan = in->fspan;         //Frequency span
   noiseA = sqrt(in->P0*0.5/dt);
   periodf = 1.0/period0;
   block = in->blksize;
   
   numffts = (INT4)floor(2*(Tobs/dT)-1);    //Number of FFTs
   numfbins = (INT4)(roundf(fspan*dT)+1);    //Number of frequency bins
   tempfmin = fmin - roundf((block-1)*0.5)/dT;
   tempfspan = fspan + (block-1)/dT;
   tempnumfbins = (INT4)(roundf(tempfspan*dT)+1);
   printf("Number of FFTs = %d, number of frequency bins = %d\n",numffts,numfbins);
   fprintf(LOG,"Number of FFTs = %d, number of frequency bins = %d\n",numffts,numfbins);
   
   REAL4Vector *tfdata = XLALCreateREAL4Vector((UINT4)(tempnumfbins*numffts));
   
   //Create time series
   REAL4Vector *t = XLALCreateREAL4Vector((UINT4)floor(dT/dt)+1);
   for (ii=0; ii<t->length; ii++) t->data[ii] = ii*dt;
   
   //Create Hann window
   REAL4Window *win = XLALCreateHannREAL4Window(t->length);
   REAL4 winFactor = 8.0/3.0;
   
   //Create FFT plan
   REAL4FFTPlan *plan = XLALCreateForwardREAL4FFTPlan(t->length, 0 );
   
   //Initialize GSL random number generator
   gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
   srand(time(NULL));
   UINT8 randseed = rand();
   gsl_rng_set(rng, randseed);
   
   //Calculate noise across band
   REAL4Vector *Pvals = XLALCreateREAL4Vector((UINT4)(numfbins + block-1));
   INT4 flag = 0;
   for (ii=0; ii<Pvals->length; ii++) {
      Pvals->data[ii] = in->P0*(1.0 + in->Pslope*(fmin+(ii-(INT4)roundf((block-1)*0.5))/dT - 
         (fmin+0.5*fspan))/fspan + in->Pquad*((fmin+(ii-(INT4)roundf((block-1)*0.5))/dT - 
         (fmin+0.5*fspan))/fspan)*((fmin+(ii-(INT4)roundf((block-1)*0.5))/dT - (fmin+0.5*fspan))/fspan));
      if (Pvals->data[ii] < 0.0 && flag!=1) flag = 1; 
   }
   if (flag==1) {
      printf("Noise spectrum values cannot be calculated for values less than zero!!!!\n");
   }
   
   //Create values in the time series, compute PSD
   REAL4Vector *x = XLALCreateREAL4Vector(t->length);
   REAL4 windowAndMatlabFactor = winFactor/x->length*dt;
   REAL4Vector *psd = XLALCreateREAL4Vector((UINT4)floor(t->length/2)+1);
   for (ii=0; ii<numffts; ii++) {
      if (params->SFTexistlist->data[ii] > 0) {
         
         //Time series for each PSD, including noise, and window
         for (jj=0; jj<t->length; jj++) x->data[jj] = (A*sin(LAL_TWOPI*(f0 + 
            B*sin(LAL_TWOPI*periodf*(t->data[jj] + ii*0.5*dT)))*t->data[jj]))*win->data->data[jj];
      
         //Do the PSD
         INT4 check = XLALREAL4PowerSpectrum(psd,x,plan);
         if (check != 0) {
            printf("Something wrong with PSD...\n");
            fprintf(LOG,"Something wrong with PSD...\n");
         }
         
         REAL4 noisecoeff = 1.0;
         if (in->varySFTs==1) noisecoeff *= (fabsf(gaussRandNum(0.2, rng) + 2.0)*0.5);
         
         //Scale PSD by window factor, 1/N, and 1/fs to give same results as Matlab and add noise
         for (jj=0; jj<tempnumfbins; jj++) {
            REAL4 noise = 0.0;
            if (Pvals->data[jj]*noisecoeff > 0.0) noise = expRandNum(Pvals->data[jj]*noisecoeff, rng);
            tfdata->data[ii*tempnumfbins + jj] = psd->data[jj+(INT4)roundf(tempfmin*dT)]*windowAndMatlabFactor + noise;
         }
      } else {
         for (jj=0; jj<tempnumfbins; jj++) tfdata->data[ii*tempnumfbins + jj] = 0.0;
      }
      
   }
   
   XLALDestroyREAL4Vector(t);
   XLALDestroyREAL4Vector(x);
   XLALDestroyREAL4Window(win);
   XLALDestroyREAL4FFTPlan(plan);
   XLALDestroyREAL4Vector(psd);
   XLALDestroyREAL4Vector(Pvals);
   gsl_rng_free(rng);
   
   
   return tfdata;
   
} */


//////////////////////////////////////////////////////////////
// Determine the TF running mean of each SFT  -- done
REAL4Vector * tfRngMeans(REAL4Vector *tfdata, INT4 numffts, INT4 numfbins, INT4 blksize)
{

   REAL4Vector *rngMeans = XLALCreateREAL4Vector((UINT4)(numffts*numfbins));
   
   LALStatus status;
   status.statusPtr = NULL;
   REAL8 bias;
   INT4 ii, jj;
   
   //Blocksize of running median. This needs to be >=501 bins because of accuracy in
   //determination of background is essential for 2nd PSD computation.
   LALRunningMedianPar block = {blksize};
   
   //Running median bias calculation
   LALRngMedBias(&status, &bias, blksize);
   REAL4 invbias = 1.0/(REAL4)bias;
   
   REAL4Sequence *inpsd = XLALCreateREAL4Sequence((UINT4)(numfbins+blksize-1));
   REAL4Sequence *mediansout = XLALCreateREAL4Sequence((UINT4)numfbins);
   for (ii=0; ii<numffts; ii++) {
      //Determine running median value, convert to mean value
      for (jj=0; jj<numfbins+blksize-1; jj++) inpsd->data[jj] = tfdata->data[ii*(numfbins+blksize-1) + jj];
      LALSRunningMedian2(&status, mediansout, inpsd, block);
      //Now make the output medians into means by multiplying by 1/bias
      for (jj=0; jj<mediansout->length; jj++) rngMeans->data[ii*numfbins + jj] = mediansout->data[jj]*invbias;
      //for (jj=0; jj<mediansout->length; jj++) fprintf(stderr,"%g\n",rngMeans->data[ii*numfbins + jj]);
   }
   
   XLALDestroyREAL4Sequence(inpsd);
   XLALDestroyREAL4Sequence(mediansout);
   
   return rngMeans;

}



//////////////////////////////////////////////////////////////
// Do the weighting by noise variance (from tfRngMeans), mean subtraction, and antenna pattern weights  -- done
REAL4Vector * tfWeightMeanSubtract(REAL4Vector *tfdata, REAL4Vector *rngMeans, REAL4Vector *antPatternWeights, inputParamsStruct *params)
{

   REAL4Vector *out = XLALCreateREAL4Vector(tfdata->length);
   
   INT4 ii, jj;
   
   INT4 numffts = (INT4)floor(2*(params->Tobs/params->Tcoh)-1);         //Number of FFTs
   INT4 numfbins = (INT4)(roundf(params->fspan*params->Tcoh)+1);     //Number of frequency bins
   
   for (ii=0; ii<numfbins; ii++) {
      
      //Get sum of antenna pattern weight/variances for each frequency bin as a function of time
      REAL4 sumofweights = 0.0;
      for (jj=0; jj<numffts; jj++) {
         if (rngMeans->data[ii+jj*numfbins] != 0.0) sumofweights += (antPatternWeights->data[jj]*antPatternWeights->data[jj])/(rngMeans->data[ii+jj*numfbins]*rngMeans->data[ii+jj*numfbins]);
      }
      REAL4 invsumofweights = 1.0/sumofweights;
      
      //Now do mean subtraction, noise weighting, antenna pattern weighting
      for (jj=0; jj<numffts; jj++) {
         if (rngMeans->data[ii+jj*numfbins] != 0.0) out->data[ii+jj*numfbins] = invsumofweights*antPatternWeights->data[jj]*(tfdata->data[ii+jj*numfbins]/rngMeans->data[ii+jj*numfbins] - 1.0)/rngMeans->data[ii+jj*numfbins];
      }
   }
   
   return out;

}


//////////////////////////////////////////////////////////////
// Make the second FFT powers
REAL4Vector * makeSecondFFT(REAL4Vector *tfdata, inputParamsStruct *params)
{
   
   INT4 ii, jj;
   
   INT4 numffts = (INT4)floor(2*(params->Tobs/params->Tcoh)-1);    //Number of FFTs
   INT4 numfbins = (INT4)(roundf(params->fspan*params->Tcoh)+1);    //Number of frequency bins
   INT4 numfprbins = (INT4)floor(numffts*0.5)+1;
   
   REAL4Vector *ffdata = XLALCreateREAL4Vector((UINT4)(numfbins*numfprbins));
   
   //Do the second FFT
   REAL4Vector *x = XLALCreateREAL4Vector((UINT4)numffts);
   printf("length of 2nd PSD TS = %d\n",numffts);
   fprintf(LOG,"length of 2nd PSD TS = %d\n",numffts);
   REAL4Window *win = XLALCreateHannREAL4Window(x->length);
   REAL4 winFactor = 8.0/3.0;
   REAL4FFTPlan *plan = XLALCreateForwardREAL4FFTPlan(x->length, 0 );
   REAL4Vector *psd = XLALCreateREAL4Vector((UINT4)floor(x->length*0.5)+1);
   printf("length of 2nd PSD = %d\n",psd->length);
   fprintf(LOG,"length of 2nd PSD = %d\n",psd->length);
   //First loop over frequencies
   for (ii=0; ii<numfbins; ii++) {
   
      //Next, loop over times and pick the right frequency bin for each FFT and window
      for (jj=0; jj<x->length; jj++) x->data[jj] = tfdata->data[ii + jj*numfbins]*win->data->data[jj];
      
      //Make the FFT
      INT4 check = XLALREAL4PowerSpectrum(psd,x,plan);
      if (check != 0) {
         printf("Something wrong with second PSD...\n");
         fprintf(LOG,"Something wrong with second PSD...\n");
      }
      
      //Scale the data points by 1/N and window factor and (1/fs)
      //Order of vector is by second frequency then first frequency
      //for (jj=0; jj<psd->length; jj++) out->ffdata->data[psd->length*ii + jj] = 
      //   psd->data[jj]*winFactor/x->length*0.5*dT;
      for (jj=0; jj<psd->length; jj++) ffdata->data[psd->length*ii + jj] = psd->data[jj]*winFactor/x->length*0.5*params->Tcoh;
      
   }
   
   XLALDestroyREAL4Vector(x);
   XLALDestroyREAL4Vector(psd);
   XLALDestroyREAL4Window(win);
   XLALDestroyREAL4FFTPlan(plan);
   
   
   return ffdata;
   
}



//////////////////////////////////////////////////////////////
// Measure of the average noise power in each 1st FFT frequency bin  -- done
REAL4Vector * ffPlaneNoise(inputParamsStruct *param, REAL4Vector *rngMeans, REAL4Vector *antPatternWeights)
{

   INT4 ii, jj, numfbins, numffts;
   
   numfbins = (INT4)(roundf(param->fspan*param->Tcoh)+1);   //Number of frequency bins
   numffts = (INT4)floor(2*(param->Tobs/param->Tcoh)-1);     //Number of FFTs
   
   //Initialize the random number generator
   gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
   srand(time(NULL));
   UINT8 randseed = rand();
   gsl_rng_set(rng, randseed);
   
   //Mean value of F_n^4
   REAL4Vector *sqAntWeights = XLALCreateREAL4Vector(antPatternWeights->length);
   for (ii=0; ii<antPatternWeights->length; ii++) sqAntWeights->data[ii] = antPatternWeights->data[ii]*antPatternWeights->data[ii];
   REAL4 sqAntWeightsMean = calcMean(sqAntWeights);
   
   //REAL4 bandMean = calcMean(rngMeans);
   REAL4Vector *aveNoise = XLALCreateREAL4Vector((UINT4)numfbins);
   REAL4Vector *rngMeansInFreqBin = XLALCreateREAL4Vector((UINT4)numffts);
   for (ii=0; ii<numfbins; ii++) {
   
      REAL4 sumofinvvariances = 0.0;
      for (jj=0; jj<numffts; jj++) if (rngMeans->data[jj*numfbins + ii] != 0.0) sumofinvvariances += (antPatternWeights->data[jj]*antPatternWeights->data[jj])/(rngMeans->data[jj*numfbins + ii]*rngMeans->data[jj*numfbins + ii]);
      REAL4 invsumofinvvariances = 1.0/sumofinvvariances;
   
      /* for (jj=0; jj<numffts; jj++) {
         if (param->SFTexistlist->data[jj] != 0) {
            rngMeansInFreqBin->data[jj] = (expRandNum(rngMeans->data[jj*numfbins + ii], rng)/
               rngMeans->data[jj*numfbins + ii] - 1.0)/rngMeans->data[jj*numfbins + ii]*invsumofinvvariances;
            rngMeansInFreqBin->data[jj] *= rngMeansInFreqBin->data[jj];
         } else {
            rngMeansInFreqBin->data[jj] = 0.0;
         }
      } */
      
      //TODO: See if there is a more efficient way of calculating the expected background in 2nd FFT
      //This code is here because the above commented out code gives poor results. Need to work on this
      /*for (jj=0; jj<numffts; jj++) {
         //if (param->SFTexistlist->data[jj] != 0) {
         if (rngMeans->data[jj] != 0.0) {
            REAL4 valforbin = 0.0;
            for (kk=0; kk<1000; kk++) {
               REAL4 val = (expRandNum(rngMeans->data[jj*numfbins + ii], rng)/rngMeans->data[jj*numfbins + ii] - 
                  1.0)/rngMeans->data[jj*numfbins + ii]*invsumofinvvariances*antPatternWeights->data[jj];
               val *= val;
               valforbin += val;
            }
            rngMeansInFreqBin->data[jj] = valforbin/kk; //Taking the average
         } else {
            rngMeansInFreqBin->data[jj] = 0.0;
         }
      }
      aveNoise->data[ii] = calcMean(rngMeansInFreqBin)*param->Tcoh; */
      
      //I think this is the efficient code I was trying to write in the first place
      for (jj=0; jj<numffts; jj++) {
         if (rngMeans->data[jj*numfbins + ii] != 0.0) rngMeansInFreqBin->data[jj] = 1.0/rngMeans->data[jj*numfbins + ii]/rngMeans->data[jj*numfbins + ii];
         else rngMeansInFreqBin->data[jj] = 0.0;
      }
      aveNoise->data[ii] = calcMean(rngMeansInFreqBin)*invsumofinvvariances*invsumofinvvariances*sqAntWeightsMean*param->Tcoh;
      //fprintf(stderr,"Plane noise in bin %d = %g\n",ii,aveNoise->data[ii]);
   }
   
   XLALDestroyREAL4Vector(rngMeansInFreqBin);
   XLALDestroyREAL4Vector(sqAntWeights);
   gsl_rng_free(rng);
   
   return aveNoise;

}





//////////////////////////////////////////////////////////////
// Calculate the R statistic
REAL4 calculateR(REAL4Vector *ffdata, REAL4Vector *template, topbinsStruct *topbinsstruct, REAL4Vector *noise)
{
   
   INT4 ii;
   
   REAL4 sumofsqweights = 0.0;
   for (ii=0; ii<topbinsstruct->topbins->length; ii++) {
      sumofsqweights += (template->data[ topbinsstruct->topbins->data[ii] ]*
         template->data[ topbinsstruct->topbins->data[ii] ]);
   }
   
   REAL4 R = 0.0;
   for (ii=0; ii<topbinsstruct->topbins->length; ii++) {
      R += ((ffdata->data[ topbinsstruct->topbins->data[ii] ] - 
         noise->data[ topbinsstruct->freqoftopbins->data[ii] ])*template->data[ topbinsstruct->topbins->data[ii] ]);
   }
   R /= sumofsqweights;
   
   return R;
   
}


//////////////////////////////////////////////////////////////
// Allocate topbinsStruct vectors  -- done
topbinsStruct * new_topbinsStruct(INT4 number)
{

   topbinsStruct *topbinsstruct = (topbinsStruct*)XLALMalloc(sizeof(topbinsStruct));
   topbinsstruct->topbins = XLALCreateINT4Vector((UINT4)number);
   topbinsstruct->freqoftopbins = XLALCreateINT4Vector((UINT4)number);

   return topbinsstruct;

}



//////////////////////////////////////////////////////////////
// Free topbinsStruct vectors  -- done
void free_topbinsStruct(topbinsStruct *topbinsstruct)
{

   XLALDestroyINT4Vector(topbinsstruct->topbins);
   XLALDestroyINT4Vector(topbinsstruct->freqoftopbins);
   XLALFree((topbinsStruct*)topbinsstruct);

}


//////////////////////////////////////////////////////////////
// Find the locations of the top weighted bins
void topbins(topbinsStruct *out, ffdataStruct *in, INT4 number)
{
   
   INT4 ii, jj;
   REAL4Vector *backup = XLALCreateREAL4Vector((UINT4)number);
   
   //Find the maximum powers in the template
   for (ii=0; ii<number; ii++) {
      REAL4 max = 0.0;
      INT4 maxloc = 0;
      for (jj=0; jj<in->ffdata->length; jj++) {
         if (in->ffdata->data[jj] > max) {
            max = in->ffdata->data[jj];
            maxloc = jj;
         }
      }
      out->topbins->data[ii] = maxloc;
      out->freqoftopbins->data[ii] = (INT4)(floor(maxloc/in->fpr->length));
      backup->data[ii] = in->ffdata->data[maxloc];
      in->ffdata->data[maxloc] = 0.0;
   }
   
   //Making sure that the original data still exists
   for (ii=0; ii<backup->length; ii++) {
      in->ffdata->data[out->topbins->data[ii]] = backup->data[ii];
   }
   XLALDestroyREAL4Vector(backup);
   

}




//////////////////////////////////////////////////////////////
// Calculates maximum modulation depth
REAL4 maxModDepth(REAL4 period, REAL4 cohtime)
{

   REAL4 maxB = 0.25*period/cohtime/cohtime;
   
   return maxB;

}


//////////////////////////////////////////////////////////////
// Calculates minimum period allowable for modulation depth and Tcoh
REAL4 minPeriod(REAL4 moddepth, REAL4 cohtime)
{

   REAL4 minP = 4.0*moddepth*cohtime*cohtime;
   
   return minP;

}


//////////////////////////////////////////////////////////////
// Calculates y = sin(pi*x)/(pi*x)/(x^2-1)
REAL4 sincxoverxsqminusone(REAL4 x)
{
   
   REAL4 val;
   
   if (x==1.0 || x==-1.0) val = -0.5;
   else val = sinc(x)/(x*x-1);
   
   return val;
   
}


//////////////////////////////////////////////////////////////
// Calculates y = sin(pi*x)/(pi*x)
REAL4 sinc(REAL4 x)
{

   REAL4 val;
   
   if (x==0.0) val = 1.0;
   else val = sin(LAL_PI*x)/(LAL_PI*x);
   
   return val;

}



//////////////////////////////////////////////////////////////
// Compute the mean value of a vector of values
REAL4 calcMean(REAL4Vector *vector)
{

   INT4 ii;
   REAL4 meanval = 0;
   for (ii=0; ii<vector->length; ii++) meanval += vector->data[ii];
   meanval = meanval / vector->length;
   
   return meanval;

}


//////////////////////////////////////////////////////////////
// Compute the standard deviation of a vector of values
REAL4 calcStddev(REAL4Vector *vector)
{

   INT4 ii;
   REAL4 meanval = calcMean(vector);
   REAL4 stddev = 0;
   for (ii=0; ii<vector->length; ii++) stddev += (vector->data[ii] - meanval)*(vector->data[ii] - meanval);
   stddev = sqrt( stddev / (vector->length - 1) );
   
   return stddev;

}



//////////////////////////////////////////////////////////////
// Compute the RMS of a vector of values
REAL4 calcRms(REAL4Vector *vector)
{

   INT4 ii;
   REAL4 rms = 0;
   REAL4Vector *sqvector = XLALCreateREAL4Vector(vector->length);
   for (ii=0; ii<vector->length; ii++) sqvector->data[ii] = (vector->data[ii]*vector->data[ii]);
   rms = calcMean(sqvector);
   rms = sqrt(rms);
   
   XLALDestroyREAL4Vector(sqvector);
   
   return rms;

}











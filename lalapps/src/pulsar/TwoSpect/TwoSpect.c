/*
*  Copyright (C) 2010 Evan Goetz
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
#include <lal/DopplerScan.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "cmdline.h"
#include "IHS.h"
#include "candidates.h"
#include "templates.h"
#include "antenna.h"
#include "TwoSpect.h"


//Global variables
candidate *ihsCandidates[10000], *gaussCandidates1[10000], *gaussCandidates2[10000], *gaussCandidates3[10000], *gaussCandidates4[10000], *exactCandidates1[10000], *exactCandidates2[10000];
inputParamsStruct *inputParams;
REAL4FFTPlan *secondFFTplan;

FILE *LOG, *TFDATA, *FFDATA;

CHAR *earth_ephemeris = NULL, *sun_ephemeris = NULL;



//Main program
int main(int argc, char *argv[])
{

   INT4 ii, jj, kk, ll, cols, numofcandidates, numofcandidates2, numofcandidatesadded;
   REAL4 ihsfarthresh, templatefarthresh;
   LALStatus status;
   status.statusPtr = NULL;
   char s[20000], t[20000], u[20000];

   struct gengetopt_args_info args_info;
   if ( cmdline_parser(argc, argv, &args_info) ) exit(-1);
   
   //Create directory
   if (args_info.outdirectory_given) {
      mkdir(args_info.outdirectory_arg, 0777);
      snprintf(s, 20000, "%s/logfile.txt", args_info.outdirectory_arg);
      snprintf(t, 20000, "%s/tfdata.dat", args_info.outdirectory_arg);
      snprintf(u, 20000, "%s/ffdata.dat", args_info.outdirectory_arg);
   } else {
      mkdir("output",0777);
      snprintf(s, 20000, "%s/logfile.txt", "output");
      snprintf(t, 20000, "%s/tfdata.dat", "output");
      snprintf(u, 20000, "%s/ffdata.dat", "output");
   }
   
   LOG = fopen(s,"w");
   fprintf(LOG,"Starting TwoSpect analysis...\n");
   fprintf(stderr,"Starting TwoSpect analysis...\n");
   
   inputParams = new_inputParams();
   
   if (args_info.Tobs_given) inputParams->Tobs = args_info.Tobs_arg;
   else inputParams->Tobs = 3*168*3600;
   if (args_info.fmin_given) inputParams->fmin = args_info.fmin_arg;
   else inputParams->fmin = 99.9;
   if (args_info.fspan_given) inputParams->fspan = args_info.fspan_arg;
   else inputParams->fspan = 0.2;
   if (args_info.cols_given) cols = args_info.cols_arg;
   else cols = 20;
   if (args_info.t0_given) inputParams->searchstarttime = args_info.t0_arg;
   else inputParams->searchstarttime = 900000000.0;
   
   //Defaults given or option passed
   inputParams->Tcoh = args_info.Tcoh_arg;
   ihsfarthresh = args_info.ihsfar_arg;
   templatefarthresh = args_info.tmplfar_arg;
   inputParams->blksize = args_info.blksize_arg;
   earth_ephemeris = (CHAR*)XLALMalloc(strlen(args_info.ephemDir_arg)+20);
   sun_ephemeris = (CHAR*)XLALMalloc(strlen(args_info.ephemDir_arg)+20);
   sprintf(earth_ephemeris,"%s/earth05-09.dat",args_info.ephemDir_arg);
   sprintf(sun_ephemeris,"%s/sun05-09.dat",args_info.ephemDir_arg);
   inputParams->dopplerMultiplier = args_info.dopplerMultiplier_arg;
   inputParams->templatelength = args_info.templateLength_arg;
   
   //Blocksize should be an odd number
   if (inputParams->blksize % 2 != 1) inputParams->blksize += 1;
   
   //Adjust maximum columns, if necessary
   if (cols > 2.0*maxModDepth(inputParams->Tobs*0.2, inputParams->Tcoh)*inputParams->Tcoh) {
      cols = (INT4)floorf(2.0*maxModDepth(inputParams->Tobs*0.2, inputParams->Tcoh)*inputParams->Tcoh);
      fprintf(LOG,"WARNING! Adjusting number of columns due to maximum modulation depth allowed\n");
      fprintf(stderr,"WARNING! Adjusting number of columns due to maximum modulation depth allowed\n");
   }
   if (cols > (INT4)roundf(inputParams->fspan*inputParams->Tcoh)+1) {
      cols = (INT4)roundf(inputParams->fspan*inputParams->Tcoh)+1;
      fprintf(LOG,"WARNING! Adjusting number of columns due to frequency span of band\n");
      fprintf(stderr,"WARNING! Adjusting number of columns due to frequency span of band\n");
   }
   
   //Parameters for the sky-grid
   CHAR *sky = (CHAR*)XLALMalloc(strlen(args_info.skyRegion_arg));
   sprintf(sky,"%s",args_info.skyRegion_arg);
   DopplerSkyScanInit scanInit = empty_DopplerSkyScanInit;
   DopplerSkyScanState scan = empty_DopplerSkyScanState;
   PulsarDopplerParams dopplerpos;
   scanInit.gridType = 1;     //Default value for an approximate-isotropic grid
   scanInit.skyRegionString = sky;      //"allsky" = Default value for all-sky search
   scanInit.numSkyPartitions = 1;   //Default value so sky is not broken into chunks
   scanInit.Freq = (REAL8)args_info.fmin_arg;  
   
   fprintf(LOG,"Tobs = %f sec\n",inputParams->Tobs);
   fprintf(LOG,"Tcoh = %f sec\n",inputParams->Tcoh);
   fprintf(LOG,"fmin = %f Hz\n",inputParams->fmin);
   fprintf(LOG,"fspan = %f Hz\n",inputParams->fspan);
   fprintf(LOG,"Running median blocksize = %d\n",inputParams->blksize);
   fprintf(stderr,"Tobs = %f sec\n",inputParams->Tobs);
   fprintf(stderr,"Tcoh = %f sec\n",inputParams->Tcoh);
   fprintf(stderr,"fmin = %f Hz\n",inputParams->fmin);
   fprintf(stderr,"fspan = %f Hz\n",inputParams->fspan);
   fprintf(stderr,"Running median blocksize = %d\n",inputParams->blksize);
   
   //Basic units
   INT4 numffts = (INT4)floor(2*(inputParams->Tobs/inputParams->Tcoh)-1);
   INT4 numfbins = (INT4)(roundf(inputParams->fspan*inputParams->Tcoh)+1);
   REAL4 tempfspan = inputParams->fspan + (inputParams->blksize-1)/inputParams->Tcoh;
   INT4 tempnumfbins = (INT4)(roundf(tempfspan*inputParams->Tcoh)+1);
   
   //Allocate memory for ffdata structure
   ffdataStruct *ffdata = new_ffdata(inputParams, 0);
   
   //Second fft plan, only need to make this once for all the exact templates
   secondFFTplan = XLALCreateForwardREAL4FFTPlan((UINT4)floor(2*(inputParams->Tobs/inputParams->Tcoh)-1), 0);
   
   //Detector velocity during SFTs
   LALDetector det = lalCachedDetectors[LALDetectorIndexLHODIFF]; //H1
   inputParams->det = &det;
   EphemerisData *edat = new_Ephemeris(earth_ephemeris, sun_ephemeris);
   REAL4 detectorVmax = 9.93e-5;   //Average earth speed in units of c
   
   //Initialize the sky-grid
   scanInit.dAlpha = (REAL8)0.5/(inputParams->fmin * inputParams->Tcoh * detectorVmax);
   scanInit.dDelta = scanInit.dAlpha;
   InitDopplerSkyScan(&status, &scan, &scanInit);
   XLALNextDopplerSkyPos(&dopplerpos, &scan); //Start at first location
   
   //Find the FAR of IHS sum
   fprintf(stderr,"Determining IHS FAR values\n");
   ihsfarStruct *ihsfarstruct = new_ihsfarStruct(cols);
   genIhsFar(ihsfarstruct, ffdata, cols, ihsfarthresh);
   fprintf(LOG,"Maximum column width to be searched = %d, at threshold = %f\n",cols,ihsfarthresh);
   fprintf(stderr,"Maximum column width to be searched = %d, at threshold = %f\n",cols,ihsfarthresh);
   //IHS FOM FAR (allows a relative offset of +/- 1 bin between maximum values
   REAL4 ihsfomfar = 6.0;
   fprintf(LOG,"IHS FOM FAR = %f\n",ihsfomfar);
   fprintf(stderr,"IHS FOM FAR = %f\n",ihsfomfar);
   
   //Assume maximum bin shift possible
   inputParams->maxbinshift = detectorVmax * inputParams->fmin * inputParams->Tcoh; //TODO: better way to do this?
   
   //Read in the T-F data
   fprintf(LOG,"Loading in SFTs... ");
   fprintf(stderr,"Loading in SFTs... ");
   REAL4Vector *tfdata = readInSFTs(inputParams);
   fprintf(LOG,"done\n");
   fprintf(stderr,"done\n");
   
   //Calculate the running mean values of the SFTs (output here is smaller than initialTFdata). Here,
   //numfbins needs to be the bins you expect to come out of the running means -- the band you are going
   //to search!
   fprintf(LOG,"Assessing background... ");
   fprintf(stderr,"Assessing background... ");
   REAL4Vector *background = tfRngMeans(tfdata, numffts, numfbins + 2*inputParams->maxbinshift, inputParams->blksize);
   fprintf(LOG,"done\n");
   fprintf(stderr,"done\n");
   
   //Need to reduce the original TF data to remove the excess bins used for running median calculation
   REAL4Vector *usableTFdata = XLALCreateREAL4Vector(background->length);
   for (ii=0; ii<numffts; ii++) {
      for (jj=0; jj<numfbins+2*inputParams->maxbinshift; jj++) usableTFdata->data[ii*(numfbins+2*inputParams->maxbinshift) + jj] = tfdata->data[ii*(tempnumfbins+2*inputParams->maxbinshift) + jj + (INT4)roundf(0.5*(inputParams->blksize-1))];
   }
   //At this point the TF plane and the running median calculation are the same size=numffts*(numfbins + 2*maxbinshift)
   
   numofcandidates = numofcandidates2 = numofcandidatesadded = 0;
   
   //Search over the sky
   while (scan.state != STATE_FINISHED) {
      fprintf(LOG,"Sky location: RA = %g, DEC = %g\n", dopplerpos.Alpha, dopplerpos.Delta);
      fprintf(stderr,"Sky location: RA = %g, DEC = %g\n", dopplerpos.Alpha, dopplerpos.Delta);
      
      numofcandidates = numofcandidates2 = 0;
      
      //Determine detector velocity w.r.t. a sky location for each SFT
      REAL4Vector *detectorVelocities = CompAntennaVelocity((REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, inputParams->searchstarttime, inputParams->Tcoh, inputParams->Tobs, det, edat);
      
      //Compute the bin shifts for each SFT
      INT4Vector *binshifts = CompBinShifts(inputParams->fmin+inputParams->fspan*0.5, detectorVelocities, inputParams->Tcoh, inputParams->dopplerMultiplier);
      
      //Compute antenna pattern weights
      REAL4Vector *antenna = CompAntennaPatternWeights((REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, inputParams->searchstarttime, inputParams->Tcoh, inputParams->Tobs, det);
      memcpy(ffdata->antweights->data, antenna->data, antenna->length*sizeof(*antenna->data));
      
      //Slide SFTs here -- need to slide the data and the estimated background
      REAL4Vector *initialTFdata = slideTFdata(inputParams, usableTFdata, binshifts);
      REAL4Vector *backgroundslide = slideTFdata(inputParams, background, binshifts);
      memcpy(ffdata->backgrnd->data, backgroundslide->data, backgroundslide->length*sizeof(*backgroundslide->data));
      
      //Compute the weighted TF data
      REAL4Vector *weightedTFdata = tfWeightMeanSubtract(initialTFdata, ffdata->backgrnd, ffdata->antweights, inputParams);
      
      //Do the second FFT
      REAL4Vector *secFFTdata = makeSecondFFT(weightedTFdata, inputParams);
      memcpy(ffdata->ffdata->data, secFFTdata->data, secFFTdata->length*sizeof(*secFFTdata->data));
      fprintf(stderr,"Average FF = %g\n",calcMean(ffdata->ffdata));
      fprintf(stderr,"Std dev. FF = %g\n",calcStddev(ffdata->ffdata));
      
      //Average noise floor of FF plane for each 1st FFT frequency bin
      REAL4Vector *aveNoise = ffPlaneNoise(inputParams, ffdata->backgrnd, ffdata->antweights);
      fprintf(stderr,"Average expected noise = %g\n",calcMean(aveNoise));
      fprintf(stderr,"Std dev. expected noise = %g\n",calcStddev(aveNoise));
      
////////Start of the IHS step!
      ihsMaximaStruct *ihsmaxima = new_ihsMaxima(ffdata, cols);
      
      //Run the IHS algorithm on the data
      runIHS(ihsmaxima, ffdata, cols);
      
      //Find any IHS candidates
      fprintf(LOG,"Checking IHS values for candidates...\n");
      fprintf(stderr,"Checking IHS values for candidates...\n");
      findIHScandidates(ihsCandidates, &numofcandidates, ihsfarstruct, aveNoise, inputParams, ffdata, ihsmaxima, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta);
      fprintf(LOG,"done\n");
      fprintf(stderr,"done\n");
      fprintf(LOG,"Candidates found in IHS step = %d\n",numofcandidates);
      fprintf(stderr,"Candidates found in IHS step = %d\n",numofcandidates);
////////End of the IHS step
      
////////Start of the Gaussian template search!
      farStruct *farval = NULL;
      fprintf(LOG,"Starting Gaussian template search...\n");
      fprintf(stderr,"Starting Gaussian template search...\n");
      for (ii=0; ii<numofcandidates; ii++) {
         if (ihsCandidates[ii]->fsig-ihsCandidates[ii]->moddepth-6/inputParams->Tcoh > inputParams->fmin && ihsCandidates[ii]->fsig+ihsCandidates[ii]->moddepth+6/inputParams->Tcoh < inputParams->fmin+inputParams->fspan && ihsCandidates[ii]->moddepth < maxModDepth(ihsCandidates[ii]->period,inputParams->Tcoh) && ihsCandidates[ii]->period >= 2.0*3600.0) {
            
            /* TFDATA = fopen(t,"w");
            FFDATA = fopen(u,"w");
            for (jj=0; jj<(INT4)weightedTFdata->length; jj++) fprintf(TFDATA,"%g\n",weightedTFdata->data[jj]);
            for (jj=0; jj<(INT4)ffdata->ffdata->length; jj++) fprintf(FFDATA,"%g\n",ffdata->ffdata->data[jj]);
            fclose(TFDATA);
            fclose(FFDATA); */
            
            //Allocate memory for template
            templateStruct *template = new_templateStruct(inputParams->templatelength);
            
            //Make a Gaussian train template
            makeTemplateGaussians(template, ihsCandidates[ii], inputParams);
            
            //Estimate the FAR for these bin weights
            farval = new_farStruct();
            estimateFAR(farval, template, (INT4)roundf(10000*.01/templatefarthresh), templatefarthresh, aveNoise);
            
            //Caclulate R
            REAL4 R = calculateR(ffdata->ffdata, template, aveNoise);
            
            //Destroy unneeded things
            free_templateStruct(template);
            template = NULL;
            
            /* Log the candidate if R exceeds the FAR or check other possibilities of different 
            periods. Use same farval.far because the FAR is ~independent on period. */
            REAL4 Rfirst = R;
            REAL4 Rtemp = R;
            REAL4 bestPeriod = ihsCandidates[ii]->period;
            if (Rfirst > farval->far) {
               gaussCandidates1[numofcandidates2] = new_candidate();
               loadCandidateData(gaussCandidates1[numofcandidates2], ihsCandidates[ii]->fsig, ihsCandidates[ii]->period, ihsCandidates[ii]->moddepth, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, Rfirst, (Rfirst-farval->distMean)/farval->distSigma, 0.0);
               numofcandidates2++;
            } else {
               if (ihsCandidates[ii]->period*0.5 > minPeriod(ihsCandidates[ii]->moddepth, inputParams->Tcoh) && ihsCandidates[ii]->period*0.5 >= 2.0*3600.0) {
                  ihsCandidates[ii]->period *= 0.5;
                  template = new_templateStruct(inputParams->templatelength);
                  makeTemplateGaussians(template, ihsCandidates[ii], inputParams);
                  R = calculateR(ffdata->ffdata, template, aveNoise);
                  if (R > Rtemp) {
                     bestPeriod = ihsCandidates[ii]->period;
                     Rtemp = R;
                  }
                  free_templateStruct(template);
                  template = NULL;
                  ihsCandidates[ii]->period *= 2.0;
               }
               if (ihsCandidates[ii]->period/3.0 > minPeriod(ihsCandidates[ii]->moddepth, inputParams->Tcoh) && ihsCandidates[ii]->period/3.0 >= 2.0*3600.0) {
                  ihsCandidates[ii]->period /= 3.0;
                  template = new_templateStruct(inputParams->templatelength);
                  makeTemplateGaussians(template, ihsCandidates[ii], inputParams);
                  R = calculateR(ffdata->ffdata, template, aveNoise);
                  if (R > Rtemp) {
                     bestPeriod = ihsCandidates[ii]->period;
                     Rtemp = R;
                  }
                  free_templateStruct(template);
                  template = NULL;
                  ihsCandidates[ii]->period *= 3.0;
               }
               if (ihsCandidates[ii]->period*0.25 > minPeriod(ihsCandidates[ii]->moddepth, inputParams->Tcoh) && ihsCandidates[ii]->period*0.25 >= 2.0*3600.0) {
                  ihsCandidates[ii]->period *= 0.25;
                  template = new_templateStruct(inputParams->templatelength);
                  makeTemplateGaussians(template, ihsCandidates[ii], inputParams);
                  R = calculateR(ffdata->ffdata, template, aveNoise);
                  if (R > Rtemp) {
                     bestPeriod = ihsCandidates[ii]->period;
                     Rtemp = R;
                  }
                  free_templateStruct(template);
                  template = NULL;
                  ihsCandidates[ii]->period *= 4.0;
               }
               if (ihsCandidates[ii]->period*0.2 > minPeriod(ihsCandidates[ii]->moddepth, inputParams->Tcoh) && ihsCandidates[ii]->period*0.2 >= 2.0*3600.0) {
                  ihsCandidates[ii]->period *= 0.2;
                  template = new_templateStruct(inputParams->templatelength);
                  makeTemplateGaussians(template, ihsCandidates[ii], inputParams);
                  R = calculateR(ffdata->ffdata, template, aveNoise);
                  if (R > Rtemp) {
                     bestPeriod = ihsCandidates[ii]->period;
                     Rtemp = R;
                  }
                  free_templateStruct(template);
                  template = NULL;
                  ihsCandidates[ii]->period *= 5.0;
               }
               /* if (ihsCandidates[ii]->period*0.5 > minPeriod(ihsCandidates[ii]->moddepth, inputParams->Tcoh) && ihsCandidates[ii]->period*0.5 >= 2.0*3600.0) {
                  ihsCandidates[ii]->period *= 0.5;
                  template = new_templateStruct(inputParams->templatelength);
                  makeTemplateGaussians(template, ihsCandidates[ii], inputParams);
                  R = calculateR(ffdata->ffdata, template, aveNoise);
                  if (R > Rtemp) {
                     bestPeriod = ihsCandidates[ii]->period;
                     Rtemp = R;
                  }
                  free_templateStruct(template);
                  template = NULL;
                  ihsCandidates[ii]->period *= 2.0;
               }
               if (ihsCandidates[ii]->period*2.0/3.0 > minPeriod(ihsCandidates[ii]->moddepth, inputParams->Tcoh) && ihsCandidates[ii]->period*2.0/3.0 >= 2.0*3600.0) {
                  ihsCandidates[ii]->period *= 2.0/3.0;
                  template = new_templateStruct(inputParams->templatelength);
                  makeTemplateGaussians(template, ihsCandidates[ii], inputParams);
                  R = calculateR(ffdata->ffdata, template, aveNoise);
                  if (R > Rtemp) {
                     bestPeriod = ihsCandidates[ii]->period;
                     Rtemp = R;
                  }
                  free_templateStruct(template);
                  template = NULL;
                  ihsCandidates[ii]->period *= 1.5;
               }
               if (ihsCandidates[ii]->period*0.75 > minPeriod(ihsCandidates[ii]->moddepth, inputParams->Tcoh) && ihsCandidates[ii]->period*0.75 >= 2.0*3600.0) {
                  ihsCandidates[ii]->period *= 0.75;
                  template = new_templateStruct(inputParams->templatelength);
                  makeTemplateGaussians(template, ihsCandidates[ii], inputParams);
                  R = calculateR(ffdata->ffdata, template, aveNoise);
                  if (R > Rtemp) {
                     bestPeriod = ihsCandidates[ii]->period;
                     Rtemp = R;
                  }
                  free_templateStruct(template);
                  template = NULL;
                  ihsCandidates[ii]->period /= 0.75;
               }
               if (inputParams->Tobs/ihsCandidates[ii]->period*0.75 > 5) {
                  ihsCandidates[ii]->period *= 4.0/3.0;
                  template = new_templateStruct(inputParams->templatelength);
                  makeTemplateGaussians(template, ihsCandidates[ii], inputParams);
                  R = calculateR(ffdata->ffdata, template, aveNoise);
                  if (R > Rtemp) {
                     bestPeriod = ihsCandidates[ii]->period;
                     Rtemp = R;
                  }
                  free_templateStruct(template);
                  template = NULL;
                  ihsCandidates[ii]->period *= 0.75;
               }
               if (inputParams->Tobs/ihsCandidates[ii]->period/1.5 > 5) {
                  ihsCandidates[ii]->period *= 1.5;
                  template = new_templateStruct(inputParams->templatelength);
                  makeTemplateGaussians(template, ihsCandidates[ii], inputParams);
                  R = calculateR(ffdata->ffdata, template, aveNoise);
                  if (R > Rtemp) {
                     bestPeriod = ihsCandidates[ii]->period;
                     Rtemp = R;
                  }
                  free_templateStruct(template);
                  template = NULL;
                  ihsCandidates[ii]->period /= 1.5;
               }
               if (inputParams->Tobs/ihsCandidates[ii]->period*0.5 > 5) {
                  ihsCandidates[ii]->period *= 2.0;
                  template = new_templateStruct(inputParams->templatelength);
                  makeTemplateGaussians(template, ihsCandidates[ii], inputParams);
                  R = calculateR(ffdata->ffdata, template, aveNoise);
                  if (R > Rtemp) {
                     bestPeriod = ihsCandidates[ii]->period;
                     Rtemp = R;
                  }
                  free_templateStruct(template);
                  template = NULL;
                  ihsCandidates[ii]->period *= 0.5;
               } */
               R = Rtemp;
            }
            
            if (R > farval->far) {
               ihsCandidates[ii]->period = bestPeriod;
               //for (jj=0; jj<10; jj++) {
               for (jj=0; jj<5; jj++) {
                  //REAL4 periodfact = (jj*2+1)*ihsCandidates[ii]->period/(jj+1)*0.5;
                  REAL4 periodfact = (jj+1)*ihsCandidates[ii]->period/(jj+2);
                  if ( (ihsCandidates[ii]->period - periodfact) > minPeriod(ihsCandidates[ii]->moddepth, inputParams->Tcoh) && ihsCandidates[ii]->period-periodfact>2.0*3600.0) {
                     ihsCandidates[ii]->period -= periodfact;
                     template = new_templateStruct(inputParams->templatelength);
                     makeTemplateGaussians(template, ihsCandidates[ii], inputParams);
                     R = calculateR(ffdata->ffdata, template, aveNoise);
                     if (R > Rtemp) {
                        bestPeriod = ihsCandidates[ii]->period;
                        Rtemp = R;
                     }
                     free_templateStruct(template);
                     template = NULL;
                     ihsCandidates[ii]->period += periodfact;
                  }
                  /* periodfact = ihsCandidates[ii]->period*0.5;
                  if ( inputParams->Tobs/(ihsCandidates[ii]->period + (jj+1)*periodfact) > 5) { */
                  if ( inputParams->Tobs/(ihsCandidates[ii]->period + periodfact) > 5) {
                     ihsCandidates[ii]->period += periodfact;
                     template = new_templateStruct(inputParams->templatelength);
                     makeTemplateGaussians(template, ihsCandidates[ii], inputParams);
                     R = calculateR(ffdata->ffdata, template, aveNoise);
                     if (R > Rtemp) {
                        bestPeriod = ihsCandidates[ii]->period;
                        Rtemp = R;
                     }
                     free_templateStruct(template);
                     template = NULL;
                     ihsCandidates[ii]->period -= periodfact;
                  }
               }
               R = Rtemp;
               
               if (Rfirst<R && Rfirst>farval->far) {
                  free_candidate(gaussCandidates1[numofcandidates2-1]);
                  gaussCandidates1[numofcandidates2-1] = NULL;
                  gaussCandidates1[numofcandidates2-1] = new_candidate();
                  loadCandidateData(gaussCandidates1[numofcandidates2-1], ihsCandidates[ii]->fsig, bestPeriod, ihsCandidates[ii]->moddepth, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, R, (R-farval->distMean)/farval->distSigma, 0.0);
               } else if (Rfirst<R && Rfirst<=farval->far) {
                  gaussCandidates1[numofcandidates2] = new_candidate();
                  loadCandidateData(gaussCandidates1[numofcandidates2], ihsCandidates[ii]->fsig, bestPeriod, ihsCandidates[ii]->moddepth, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, R, (R-farval->distMean)/farval->distSigma, 0.0);
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
      
      ii = 0;
      while (gaussCandidates1[ii]!=NULL) {
         fprintf(stderr,"Candidate %d: f0=%g, P=%g, df=%g\n",ii,gaussCandidates1[ii]->fsig,gaussCandidates1[ii]->period,gaussCandidates1[ii]->moddepth);
         ii++;
      }
////////End of the Gaussian template search

      //Destroy IHS candidates
      for (ii=0; ii<10000; ii++) {
         if (ihsCandidates[ii]!=NULL) {
            free_candidate(ihsCandidates[ii]);
            ihsCandidates[ii] = NULL;
         } else {
            ii = 9999;
         }
      }

////////Start clustering!
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
      
      ii = 0;
      while (gaussCandidates2[ii]!=NULL) {
         fprintf(stderr,"Candidate %d: f0=%g, P=%g, df=%g\n",ii,gaussCandidates2[ii]->fsig,gaussCandidates2[ii]->period,gaussCandidates2[ii]->moddepth);
         ii++;
      }
////////End clustering

      //Destroy first set of Gaussian template candidates
      for (ii=0; ii<10000; ii++) {
         if (gaussCandidates1[ii]!=NULL) {
            free_candidate(gaussCandidates1[ii]);
            gaussCandidates1[ii] = NULL;
         } else {
            ii = 9999;
         }
      }

////////Start detailed Gaussian template search!
      fprintf(LOG,"Starting detailed search using Gaussian train templates... ");
      fprintf(stderr,"Starting detailed search using Gaussian train templates... ");
      //REAL4 quadparam = 2.4e-3;
      //REAL4 linparam = 4.1e-3;
      REAL4 tcohfactor = 1.49e-3*inputParams->Tcoh + 1.76;
      for (ii=0; ii<numofcandidates; ii++) {
      
         REAL4Vector *trialf, *trialb, *trialp;
         REAL4 minf, maxf, minb, maxb;
         UINT4 numf, numb, nump;
         
         //Set up parameters of modulation depth search
         minb = gaussCandidates2[ii]->moddepth-1.0/inputParams->Tcoh;
         maxb = gaussCandidates2[ii]->moddepth+1.0/inputParams->Tcoh;
         if (minb<0.25/inputParams->Tcoh) minb = 0.5/inputParams->Tcoh;
         numb = (UINT4)roundf(2*(maxb-minb)*inputParams->Tcoh)+1;
         trialb = XLALCreateREAL4Vector(numb);
         for (jj=0; jj<(INT4)numb; jj++) trialb->data[jj] = minb + 0.5*jj/inputParams->Tcoh;
         
         //Set up parameters of signal frequency search
         minf = gaussCandidates2[ii]->fsig-1.0/inputParams->Tcoh;
         maxf = gaussCandidates2[ii]->fsig+1.0/inputParams->Tcoh;
         if (minf<inputParams->fmin) minf = inputParams->fmin;
         if (maxf>inputParams->fmin+inputParams->fspan) maxf = inputParams->fmin+inputParams->fspan;
         numf = (UINT4)roundf(2*(maxf-minf)*inputParams->Tcoh)+1;
         trialf = XLALCreateREAL4Vector(numf);
         for (jj=0; jj<(INT4)numf; jj++) trialf->data[jj] = minf + 0.5*jj/inputParams->Tcoh;
         
         //Search over 5 different periods
         nump = 5;
         trialp = XLALCreateREAL4Vector(nump);
         
         //Now search over the parameter space. Frequency, then modulation depth, then period
         REAL4 bestf, bestp, bestdf, bestR, bestSNR;
         bestf = bestp = bestdf = bestR = bestSNR = 0.0;
         for (jj=0; jj<(INT4)trialf->length; jj++) {
            for (kk=0; kk<(INT4)trialb->length; kk++) {
               //Start with period of the first guess, then determine nearest neighbor from the
               //modulation depth amplitude to find the other period guesses. These parameters 
               //are determined from simulation to scale the N.N. distance w.r.t. mod. depth with
               //20% mismatch parameter
               INT4 midposition = (INT4)((nump-1)*0.5);
               trialp->data[midposition] = gaussCandidates2[ii]->period;
               for (ll=0; ll<midposition; ll++) {
                  REAL4 nnp = trialp->data[midposition+ll]*trialp->data[midposition+ll]*(1+trialp->data[midposition+ll]/tcohfactor/inputParams->Tobs)/tcohfactor/inputParams->Tobs*sqrt(3.6e-3/trialb->data[kk]);
                  trialp->data[midposition+(ll+1)] = trialp->data[midposition+ll] + nnp;
                  nnp = trialp->data[midposition-ll]*trialp->data[midposition-ll]*(1+trialp->data[midposition-ll]/tcohfactor/inputParams->Tobs)/tcohfactor/inputParams->Tobs*sqrt(3.6e-3/trialb->data[kk]);
                  trialp->data[midposition-(ll+1)] = trialp->data[midposition-ll] - nnp;
               }
               
               //Take the mean period and compute a template/FAR pair.
               REAL4 tempP = calcMean(trialp);
               candidate *cand = new_candidate();
               loadCandidateData(cand, trialf->data[jj], tempP, trialb->data[kk], (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, 0, 0, 0.0);
               templateStruct *template = new_templateStruct(inputParams->templatelength);
               makeTemplateGaussians(template, cand, inputParams);
               farval = new_farStruct();
               estimateFAR(farval, template, (INT4)roundf(10000*.01/templatefarthresh), templatefarthresh, aveNoise);
               free_candidate(cand);
               cand = NULL;
               free_templateStruct(template);
               template = NULL;
               for (ll=0; ll<(INT4)trialp->length; ll++) {
                  if ( trialf->data[jj]-trialb->data[kk]-6/inputParams->Tcoh > inputParams->fmin && trialf->data[jj]+trialb->data[kk]+6/inputParams->Tcoh < inputParams->fmin+inputParams->fspan && trialb->data[kk] < maxModDepth(trialp->data[ll], inputParams->Tcoh) && trialp->data[ll] > minPeriod(trialb->data[kk], inputParams->Tcoh) && inputParams->Tobs/trialp->data[ll] > 5.0 && trialp->data[ll] >= 2.0*3600.0) {
                     cand = new_candidate();
                     loadCandidateData(cand, trialf->data[jj], trialp->data[ll], trialb->data[kk], (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, 0, 0, 0.0);
                     template = new_templateStruct(inputParams->templatelength);
                     makeTemplateGaussians(template, cand, inputParams);
                     REAL4 R = calculateR(ffdata->ffdata, template, aveNoise);
                     REAL4 snr = (R - farval->distMean)/farval->distSigma;
                     if (R > farval->far && snr > bestSNR) {
                        bestf = trialf->data[jj];
                        bestp = trialp->data[ll];
                        bestdf = trialb->data[kk];
                        bestR = R;
                        bestSNR = snr;
                     }
                     free_candidate(cand);
                     cand = NULL;
                     free_templateStruct(template);
                     template = NULL;
                  }
               }
               
               free_farStruct(farval);
               farval = NULL;
            }
         }
         
         if (bestf!=0.0) {
            gaussCandidates3[numofcandidates2] = new_candidate();
            loadCandidateData(gaussCandidates3[numofcandidates2], bestf, bestp, bestdf, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, bestR, bestSNR, 0.0);
            numofcandidates2++;
         } else {
            fprintf(stderr,"WTF?!\n");
         }
         
         XLALDestroyREAL4Vector(trialf);
         XLALDestroyREAL4Vector(trialb);
         XLALDestroyREAL4Vector(trialp);
         trialf = NULL;
         trialb = NULL;
         trialp = NULL;
      }
      fprintf(LOG,"done\n");
      fprintf(stderr,"done\n");
      
      ii = 0;
      while (gaussCandidates3[ii]!=NULL) {
         fprintf(stderr,"Candidate %d: f0=%g, P=%g, df=%g\n",ii,gaussCandidates3[ii]->fsig,gaussCandidates3[ii]->period,gaussCandidates3[ii]->moddepth);
         ii++;
      }
////////End detailed Gaussian template search

      //Destroy 2nd round of Gaussian template candidates
      for (ii=0; ii<10000; ii++) {
         if (gaussCandidates2[ii]!=NULL) {
            free_candidate(gaussCandidates2[ii]);
            gaussCandidates2[ii] = NULL;
         } else {
            ii = 9999;
         }
      }

////////Start clustering!
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
      
      ii = 0;
      while (gaussCandidates4[ii]!=NULL) {
         fprintf(stderr,"Candidate %d: f0=%g, P=%g, df=%g\n",ii,gaussCandidates4[ii]->fsig,gaussCandidates4[ii]->period,gaussCandidates4[ii]->moddepth);
         ii++;
      }
////////End clustering

      //Destroy 3rd set of Gaussian template candidates
      for (ii=0; ii<10000; ii++) {
         if (gaussCandidates3[ii]!=NULL) {
            free_candidate(gaussCandidates3[ii]);
            gaussCandidates3[ii] = NULL;
         } else {
            ii = 9999;
         }
      }

////////Initial check using "exact" template
      fprintf(LOG,"Starting exact template search...\n");
      fprintf(stderr,"Starting exact template search...\n");
      for (ii=0; ii<numofcandidates; ii++) {
         templateStruct *template = new_templateStruct(inputParams->templatelength);
         makeTemplate(template, gaussCandidates4[ii], inputParams, secondFFTplan);
         farval = new_farStruct();
         estimateFAR(farval, template, (INT4)roundf(10000*.01/templatefarthresh), templatefarthresh, aveNoise);
         REAL4 R = calculateR(ffdata->ffdata, template, aveNoise);
         REAL4 SNR = (R - farval->distMean)/farval->distSigma;
         if (R > farval->far) {
            exactCandidates1[numofcandidates2] = new_candidate();
            loadCandidateData(exactCandidates1[numofcandidates2], gaussCandidates4[ii]->fsig, gaussCandidates4[ii]->period, gaussCandidates4[ii]->moddepth, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, R, SNR, 0.0);
            numofcandidates2++;
         }
         
         free_templateStruct(template);
         template = NULL;
         free_farStruct(farval);
         farval = NULL;
      }
      numofcandidates = numofcandidates2;
      numofcandidates2 = 0;
      fprintf(LOG,"Number of candidates confirmed with exact templates = %d\n",numofcandidates);
      fprintf(stderr,"Number of candidates confirmed with exact templates = %d\n",numofcandidates);
////////Done with initial check

      //Destroy 4th set of Gaussian template candidates
      for (ii=0; ii<10000; ii++) {
         if (gaussCandidates4[ii]!=NULL) {
            free_candidate(gaussCandidates4[ii]);
            gaussCandidates4[ii] = NULL;
         } else {
            ii = 9999;
         }
      }

////////Start detailed "exact" template search!
      fprintf(LOG,"Starting detailed search with exact templates...\n");
      fprintf(stderr,"Starting detailed search with exact templates...\n");
      for (ii=0; ii<numofcandidates; ii++) {
      
         REAL4Vector *trialf, *trialb, *trialp;
         REAL4 minf, maxf, minb, maxb;
         UINT4 numf, numb, nump;
         
         minb = exactCandidates1[ii]->moddepth-1.0/inputParams->Tcoh;
         maxb = exactCandidates1[ii]->moddepth+1.0/inputParams->Tcoh;
         if (minb<0.25/inputParams->Tcoh) minb = 0.5/inputParams->Tcoh;
         numb = (UINT4)roundf(2*(maxb-minb)*inputParams->Tcoh)+1;
         trialb = XLALCreateREAL4Vector(numb);
         for (jj=0; jj<(INT4)numb; jj++) trialb->data[jj] = minb + 0.5*jj/inputParams->Tcoh;

         minf = exactCandidates1[ii]->fsig-1.0/inputParams->Tcoh;
         maxf = exactCandidates1[ii]->fsig+1.0/inputParams->Tcoh;
         if (minf<inputParams->fmin) minf = inputParams->fmin;
         if (maxf>inputParams->fmin+inputParams->fspan) maxf = inputParams->fmin+inputParams->fspan;
         numf = (UINT4)roundf(2*(maxf-minf)*inputParams->Tcoh)+1;
         trialf = XLALCreateREAL4Vector(numf);
         for (jj=0; jj<(INT4)numf; jj++) trialf->data[jj] = minf + 0.5*jj/inputParams->Tcoh;
         
         nump = 5;
         trialp = XLALCreateREAL4Vector(nump);
         
         REAL4 bestf, bestp, bestdf, bestR, bestSNR, bestProb;
         bestf = bestp = bestdf = bestR = bestSNR = bestProb = 0.0;
         for (jj=0; jj<(INT4)trialf->length; jj++) {
            for (kk=0; kk<(INT4)trialb->length; kk++) {
               INT4 midposition = (INT4)((nump-1)*0.5);
               trialp->data[midposition] = exactCandidates1[ii]->period;
               for (ll=0; ll<midposition; ll++) {
                  REAL4 nnp = trialp->data[midposition+ll]*trialp->data[midposition+ll]*(1+trialp->data[midposition+ll]/tcohfactor/inputParams->Tobs)/tcohfactor/inputParams->Tobs*sqrt(3.6e-3/trialb->data[kk]);
                  trialp->data[midposition+(ll+1)] = trialp->data[midposition+ll] + nnp;
                  nnp = trialp->data[midposition-ll]*trialp->data[midposition-ll]*(1+trialp->data[midposition-ll]/tcohfactor/inputParams->Tobs)/tcohfactor/inputParams->Tobs*sqrt(3.6e-3/trialb->data[kk]);
                  trialp->data[midposition-(ll+1)] = trialp->data[midposition-ll] - nnp;
               }
               
               REAL4 tempP = calcMean(trialp);
               candidate *cand = new_candidate();
               loadCandidateData(cand, trialf->data[jj], tempP, trialb->data[kk], (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, 0, 0, 0.0);
               templateStruct *template = new_templateStruct(inputParams->templatelength);
               makeTemplate(template, cand, inputParams, secondFFTplan);
               farval = new_farStruct();
               estimateFAR(farval, template, (INT4)roundf(10000*.01/templatefarthresh), templatefarthresh, aveNoise);
               free_candidate(cand);
               cand = NULL;
               free_templateStruct(template);
               template = NULL;
               for (ll=0; ll<(INT4)trialp->length; ll++) {
                  if ( trialf->data[jj]-trialb->data[kk]-6/inputParams->Tcoh > inputParams->fmin && trialf->data[jj]+trialb->data[kk]+6/inputParams->Tcoh < inputParams->fmin+inputParams->fspan && trialb->data[kk]<maxModDepth(trialp->data[ll], inputParams->Tcoh) && trialp->data[ll] > minPeriod(trialb->data[kk], inputParams->Tcoh) && inputParams->Tobs/trialp->data[ll]>5 && trialp->data[ll] >= 2.0*3600.0) {
                     
                     cand = new_candidate();
                     loadCandidateData(cand, trialf->data[jj], trialp->data[ll], trialb->data[kk], (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, 0, 0, 0.0);
                     template = new_templateStruct(inputParams->templatelength);
                     makeTemplate(template, cand, inputParams, secondFFTplan);
                     
                     REAL4 R = calculateR(ffdata->ffdata, template, aveNoise);
                     REAL4 SNR = (R - farval->distMean)/farval->distSigma;
                     
                     if (R > farval->far && SNR > bestSNR) {
                        bestf = trialf->data[jj];
                        bestp = trialp->data[ll];
                        bestdf = trialb->data[kk];
                        bestR = R;
                        bestSNR = SNR;
                     }
                     free_candidate(cand);
                     cand = NULL;
                     free_templateStruct(template);
                     template = NULL;
                  }
               }
               
               free_farStruct(farval);
               farval = NULL;
            }
         }
         
         //Load candidate
         exactCandidates2[numofcandidates2+numofcandidatesadded] = new_candidate();
         loadCandidateData(exactCandidates2[numofcandidates2+numofcandidatesadded], bestf, bestp, bestdf, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, bestR, bestSNR, 0.0);
         
         //Best template
         templateStruct *template = new_templateStruct(inputParams->templatelength);
         makeTemplate(template, exactCandidates2[numofcandidates2+numofcandidatesadded], inputParams, secondFFTplan);
         farval = new_farStruct();
         estimateFAR(farval, template, (INT4)roundf(1000000*.01/templatefarthresh), templatefarthresh, aveNoise);
         
         //Determine log-likelihood
         REAL4 prob = 1.0;
         INT4 locinlist = (INT4)farval->topRvalues->length-1;
         if (bestR > farval->far) {
            while (locinlist > 0 && bestR > farval->topRvalues->data[locinlist-1]) locinlist--;
         }
         prob = log10((locinlist+1)/(0.01/templatefarthresh*1000000));
         
         //Load the likelihood value
         loadCandidateData(exactCandidates2[numofcandidates2+numofcandidatesadded], bestf, bestp, bestdf, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, bestR, bestSNR, prob);
         numofcandidates2++;
         
         //Clean up after log-likelihood measure
         free_templateStruct(template);
         free_farStruct(farval);
         template = NULL;
         farval = NULL;
         
         //Destroy parameter space values
         XLALDestroyREAL4Vector(trialf);
         XLALDestroyREAL4Vector(trialb);
         XLALDestroyREAL4Vector(trialp);
         trialf = NULL;
         trialb = NULL;
         trialp = NULL;
      }
      fprintf(LOG,"Exact step is done with the number of candidates = %d\n",numofcandidates2);
      fprintf(stderr,"Exact step is done with the number of candidates = %d\n",numofcandidates2);
////////End of detailed search

      //Destroy first round of exact template candidates
      for (ii=0; ii<10000; ii++) {
         if (exactCandidates1[ii]!=NULL) {
            free_candidate(exactCandidates1[ii]);
            exactCandidates1[ii] = NULL;
         } else {
            ii = 9999;
         }
      }
      
      //Add number of new candidates to the total number of candidates added
      numofcandidatesadded += numofcandidates2;
      
      //Destroy stuff
      XLALDestroyREAL4Vector(detectorVelocities);
      XLALDestroyINT4Vector(binshifts);
      XLALDestroyREAL4Vector(antenna);
      XLALDestroyREAL4Vector(initialTFdata);
      XLALDestroyREAL4Vector(backgroundslide);
      XLALDestroyREAL4Vector(weightedTFdata);
      XLALDestroyREAL4Vector(secFFTdata);
      XLALDestroyREAL4Vector(aveNoise);
      free_ihsMaxima(ihsmaxima);
      
      //Iterate to next sky location
      XLALNextDopplerSkyPos(&dopplerpos, &scan);
      
   }
   
   if (numofcandidatesadded!=0) {
      fprintf(LOG,"\n**Report of candidates:**\n");
      fprintf(stderr,"\n**Report of candidates:**\n");
      
      for (ii=0; ii<numofcandidatesadded; ii++) {
         fprintf(LOG,"fsig = %g, period = %g, df = %g, RA = %g, DEC = %g, R = %g, SNR = %g, Prob = %g\n", exactCandidates2[ii]->fsig, exactCandidates2[ii]->period, exactCandidates2[ii]->moddepth, exactCandidates2[ii]->ra, exactCandidates2[ii]->dec, exactCandidates2[ii]->stat, exactCandidates2[ii]->snr, exactCandidates2[ii]->prob);
         fprintf(stderr,"fsig = %g, period = %g, df = %g, RA = %g, DEC = %g, R = %g, SNR = %g, Prob = %g\n", exactCandidates2[ii]->fsig, exactCandidates2[ii]->period, exactCandidates2[ii]->moddepth, exactCandidates2[ii]->ra, exactCandidates2[ii]->dec, exactCandidates2[ii]->stat, exactCandidates2[ii]->snr, exactCandidates2[ii]->prob);
      }
   }
   
   //Destroy varaibles
   XLALDestroyREAL4Vector(tfdata);
   XLALDestroyREAL4Vector(background);
   XLALDestroyREAL4Vector(usableTFdata);
   free_ffdata(ffdata);
   free_ihsfarStruct(ihsfarstruct);
   free_inputParams(inputParams);
   XLALDestroyREAL4FFTPlan(secondFFTplan);
   //XLALFree((CHAR*)earth_ephemeris);
   //XLALFree((CHAR*)sun_ephemeris);
   free_Ephemeris(edat);
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
   
   for (ii=0; ii<numofcandidatesadded; ii++) {
      free_candidate(exactCandidates2[ii]);
      exactCandidates2[ii] = NULL;
   }
   
   fclose(LOG);
   
   return 0;

}
/*** End of main() ***/




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
   for (ii=0; ii<(INT4)noise->length; ii++) noise->data[ii] = gaussRandNum(sigma, ptrToGenerator);
   
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
   for (ii=0; ii<(INT4)noise->length; ii++) noise->data[ii] = expRandNum(mu, ptrToGenerator);
   
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
   XLALDestroyREAL4Vector(data->antweights);
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
   if (status.statusCode != 0) exit(-1);
   
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
            tfdata->data[ii*sftlength + jj] = (REAL4)(2.0*(sftcoeff.re*sftcoeff.re + sftcoeff.im*sftcoeff.im)); //TODO: check this for consistancy. Doing --noiseSqh=1/sqrt(1800) in MFD_v4, I need to do 2*abs(z)^2 to recover 1.
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
   
   fprintf(stderr,"TF before weighting, mean subtraction = %g\n",calcMean(tfdata));
   
   return tfdata;

}


//////////////////////////////////////////////////////////////
// Slide SFT TF data  -- 
REAL4Vector * slideTFdata(inputParamsStruct *input, REAL4Vector *tfdata, INT4Vector *binshifts)
{
   
   INT4 ii, jj;
   INT4 numffts = (INT4)floor(2*(input->Tobs/input->Tcoh)-1);
   INT4 numfbins = (INT4)(roundf(inputParams->fspan*inputParams->Tcoh)+1);
   //REAL4 tempfspan = input->fspan + (input->blksize-1)/input->Tcoh;
   //INT4 tempnumfbins = (INT4)(roundf(tempfspan*input->Tcoh)+1);
   
   REAL4Vector *outtfdata = XLALCreateREAL4Vector((UINT4)(numffts*numfbins));
   for (ii=0; ii<numffts; ii++) {
      for (jj=0; jj<numfbins; jj++) outtfdata->data[ii*numfbins + jj] = tfdata->data[ii*(numfbins+2*input->maxbinshift) + jj + input->maxbinshift + binshifts->data[ii]];
   }
   
   return outtfdata;
   
}


//////////////////////////////////////////////////////////////
// Slide SFT background TF data  -- done
REAL4Vector * slideBackgroundData(inputParamsStruct *input, REAL4Vector *background, INT4Vector *binshifts)
{
   
   INT4 ii, jj;
   INT4 numffts = (INT4)floor(2*(input->Tobs/input->Tcoh)-1);
   INT4 numfbins = (INT4)(roundf(input->fspan*input->Tcoh)+1);
   
   REAL4Vector *outtfdata = XLALCreateREAL4Vector((UINT4)(numffts*numfbins));
   for (ii=0; ii<numffts; ii++) {
      for (jj=0; jj<numfbins; jj++) outtfdata->data[ii*numfbins + jj] = background->data[ii*(numfbins+2*input->maxbinshift) + jj + input->maxbinshift + binshifts->data[ii]];
   }
   
   return outtfdata;
   
}



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
      for (jj=0; jj<(INT4)mediansout->length; jj++) rngMeans->data[ii*numfbins + jj] = mediansout->data[jj]*invbias;
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
         if (rngMeans->data[jj*numfbins + ii] != 0.0) sumofweights += (antPatternWeights->data[jj]*antPatternWeights->data[jj])/(rngMeans->data[jj*numfbins+ii]*rngMeans->data[jj*numfbins+ii]);
      }
      REAL4 invsumofweights = 1.0/sumofweights;
      
      //Now do mean subtraction, noise weighting, antenna pattern weighting
      for (jj=0; jj<numffts; jj++) {
         if (rngMeans->data[jj*numfbins+ii] != 0.0) out->data[jj*numfbins+ii] = invsumofweights*antPatternWeights->data[jj]*(tfdata->data[jj*numfbins+ii]/rngMeans->data[jj*numfbins+ii] - 1.0)/rngMeans->data[jj*numfbins+ii];
      }
   }
   
   fprintf(stderr,"TF after weighting, mean subtraction = %g\n",calcMean(out));
   
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
   //printf("length of 2nd PSD TS = %d\n",numffts);
   //fprintf(LOG,"length of 2nd PSD TS = %d\n",numffts);
   REAL4Window *win = XLALCreateHannREAL4Window(x->length);
   REAL4 winFactor = 8.0/3.0;
   REAL4FFTPlan *plan = XLALCreateForwardREAL4FFTPlan(x->length, 0 );
   REAL4Vector *psd = XLALCreateREAL4Vector((UINT4)floor(x->length*0.5)+1);
   //printf("length of 2nd PSD = %d\n",psd->length);
   //fprintf(LOG,"length of 2nd PSD = %d\n",psd->length);
   //First loop over frequencies
   for (ii=0; ii<numfbins; ii++) {
   
      //Next, loop over times and pick the right frequency bin for each FFT and window
      for (jj=0; jj<(INT4)x->length; jj++) x->data[jj] = tfdata->data[ii + jj*numfbins]*win->data->data[jj];
      
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
      for (jj=0; jj<(INT4)psd->length; jj++) ffdata->data[psd->length*ii + jj] = psd->data[jj]*winFactor/x->length*0.5*params->Tcoh;
      
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
   for (ii=0; ii<(INT4)antPatternWeights->length; ii++) sqAntWeights->data[ii] = antPatternWeights->data[ii]*antPatternWeights->data[ii];
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
//REAL4 calculateR(REAL4Vector *ffdata, REAL4Vector *template, topbinsStruct *topbinsstruct, REAL4Vector *noise)
REAL4 calculateR(REAL4Vector *ffdata, templateStruct *templatestruct, REAL4Vector *noise)
{
   
   INT4 ii;
   
   REAL4 sumofsqweights = 0.0;
   for (ii=0; ii<(INT4)templatestruct->templatedata->length; ii++) sumofsqweights += (templatestruct->templatedata->data[ii]*templatestruct->templatedata->data[ii]);
   REAL4 sumofsqweightsinv = 1.0/sumofsqweights;
   
   REAL4 R = 0.0;
   for (ii=0; ii<(INT4)templatestruct->templatedata->length; ii++) R += (ffdata->data[ templatestruct->pixellocations->data[ii] ] - noise->data[ templatestruct->firstfftfrequenciesofpixels->data[ii] ])*templatestruct->templatedata->data[ii];
   
   R *= sumofsqweightsinv;
   
   return R;
   
}



//////////////////////////////////////////////////////////////
// Calculates maximum modulation depth
REAL4 maxModDepth(REAL4 period, REAL4 cohtime)
{

   REAL4 maxB = 0.5*period/cohtime/cohtime;
   
   return maxB;

}


//////////////////////////////////////////////////////////////
// Calculates minimum period allowable for modulation depth and Tcoh
REAL4 minPeriod(REAL4 moddepth, REAL4 cohtime)
{

   REAL4 minP = 2.0*moddepth*cohtime*cohtime;
   
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
   REAL4 total = 0;
   for (ii=0; ii<(INT4)vector->length; ii++) total += vector->data[ii];
   REAL4 meanval = total / vector->length;
   
   return meanval;

}


//////////////////////////////////////////////////////////////
// Compute the standard deviation of a vector of values
REAL4 calcStddev(REAL4Vector *vector)
{

   INT4 ii;
   REAL4 meanval = calcMean(vector);
   REAL4 values = 0;
   for (ii=0; ii<(INT4)vector->length; ii++) values += (vector->data[ii] - meanval)*(vector->data[ii] - meanval);
   REAL4 stddev = sqrt( values / (vector->length - 1) );
   
   return stddev;

}



//////////////////////////////////////////////////////////////
// Compute the RMS of a vector of values
REAL4 calcRms(REAL4Vector *vector)
{

   INT4 ii;
   REAL4 rms = 0;
   REAL4Vector *sqvector = XLALCreateREAL4Vector(vector->length);
   for (ii=0; ii<(INT4)vector->length; ii++) sqvector->data[ii] = (vector->data[ii]*vector->data[ii]);
   rms = calcMean(sqvector);
   rms = sqrt(rms);
   
   XLALDestroyREAL4Vector(sqvector);
   
   return rms;

}











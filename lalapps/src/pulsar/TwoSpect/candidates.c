/*
*  Copyright (C) 2010, 2011, 2013, 2014 Evan Goetz
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

#include <math.h>
#include <lal/LALMalloc.h>
#include "candidates.h"
#include "templates.h"


//Allocate a candidateVector
candidateVector * new_candidateVector(UINT4 length)
{

   candidateVector *vector = NULL;
   XLAL_CHECK_NULL( (vector = XLALMalloc(sizeof(*vector))) != NULL, XLAL_ENOMEM );

   vector->length = length;
   vector->numofcandidates = 0;
   if (length==0) vector->data = NULL;
   else {
      XLAL_CHECK_NULL( (vector->data = XLALMalloc( length*sizeof(*vector->data) )) != NULL, XLAL_ENOMEM );
   }

   return vector;

} /* new_candidateVector() */


//Resize a candidateVector
candidateVector * resize_candidateVector(candidateVector *vector, UINT4 length)
{

   if (vector==NULL) return new_candidateVector(length);
   if (length==0) {
      free_candidateVector(vector);
      return NULL;
   }

   XLAL_CHECK_NULL( (vector->data = XLALRealloc(vector->data, length*sizeof(*vector->data))) != NULL, XLAL_ENOMEM );

   vector->length = length;

   return vector;

} /* resize_candidateVector() */


//Free a candidateVector
void free_candidateVector(candidateVector *vector)
{

   if (vector==NULL) return;
   if ((!vector->length || !vector->data) && (vector->length || vector->data)) XLAL_ERROR_VOID(XLAL_EINVAL);
   if (vector->data) XLALFree((candidate*)vector->data);
   vector->data = NULL;
   XLALFree((candidateVector*)vector);
   return;

} /* delete_candidateVector() */



//////////////////////////////////////////////////////////////
// Load candidate data
void loadCandidateData(candidate* output, REAL8 fsig, REAL8 period, REAL8 moddepth, REAL4 ra, REAL4 dec, REAL8 statval, REAL8 h0, REAL8 prob, INT4 proberrcode, REAL8 normalization)
{

   output->fsig = fsig;
   output->period = period;
   output->moddepth = moddepth;
   output->ra = ra;
   output->dec = dec;
   output->stat = statval;
   output->h0 = h0;
   output->prob = prob;
   output->proberrcode = proberrcode;
   output->normalization = normalization;

} /* loadCandidateData() */


//////////////////////////////////////////////////////////////
// Cluster candidates by frequency and period using templates:
// option = 0 uses Gaussian templates (default)
// option = 1 uses exact templates
INT4 clusterCandidates(candidateVector **output, candidateVector *input, ffdataStruct *ffdata, inputParamsStruct *params, REAL4Vector *ffplanenoise, REAL4Vector *fbinaveratios, INT4Vector *sftexist, INT4 option)
{

   INT4 ii, jj, kk, loc, loc2, numcandoutlist;
   REAL8 avefsig, aveperiod, mindf, maxdf;

   //Allocate int vectors for storage
   INT4Vector *locs = NULL, *locs2 = NULL, *usedcandidate = NULL;
   XLAL_CHECK( (locs = XLALCreateINT4Vector(input->numofcandidates)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (locs2 = XLALCreateINT4Vector(input->numofcandidates)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (usedcandidate = XLALCreateINT4Vector(input->numofcandidates)) != NULL, XLAL_EFUNC );

   //Initialize arrays
   for (ii=0; ii<(INT4)input->numofcandidates; ii++) {
      locs->data[ii] = -1;
      locs2->data[ii] = -1;
      usedcandidate->data[ii] = 0;
   }

   //Set default if bad option given
   if (option!=0 && option!=1) option = 0;

   //Make FFT plan if option 1 is given
   REAL4FFTPlan *plan = NULL;
   if (option==1) {
      XLAL_CHECK( (plan = XLALCreateForwardREAL4FFTPlan(ffdata->numffts, 1)) != NULL, XLAL_EFUNC );
   }

   numcandoutlist = 0;
   for (ii=0; ii<(INT4)input->numofcandidates; ii++) {

      //Make note of first candidate available
      locs->data[0] = ii;
      loc = 1;

      INT4 foundany = 0;   //Switch to determine if any other candidates in the group. 1 if true
      INT4 iter = 1;
      //Find any in the list that are within +1/2 bin in first FFT frequency
      for (jj=ii+1; jj<(INT4)input->numofcandidates; jj++) {
         if ( usedcandidate->data[jj] == 0 && (input->data[jj].fsig-input->data[locs->data[0]].fsig <= 0.5*iter/params->Tcoh+1.0e-6 && input->data[jj].fsig-input->data[locs->data[0]].fsig >= -0.25*iter/params->Tcoh) ) {
            locs->data[loc] = jj;
            loc++;
            if (foundany==0) foundany = 1;
         }
      } /* for jj < input->numofcandidates */
      //Keep checking as long as there are more connected frequencies going higher in frequency
      while (foundany==1) {
         foundany = 0;
         iter++;
         for (jj=ii+1; jj<(INT4)input->numofcandidates; jj++) {
            if ( usedcandidate->data[jj] == 0 && (input->data[jj].fsig-input->data[locs->data[0]].fsig-0.25/params->Tcoh <= 0.5*iter/params->Tcoh && input->data[jj].fsig-input->data[locs->data[0]].fsig+0.25/params->Tcoh >= 0.5*iter/params->Tcoh) ) {
               locs->data[loc] = jj;
               loc++;
               if (foundany==0) foundany = 1;
            }
         } /* for jj < input->numofcandidates */
      } /* while foundany==1 */
      //Now check frequencies 1/2 bin below and keep going as long as there are more connected frequencies
      foundany = 1;
      iter = 0;
      while (foundany==1) {
         foundany = 0;
         iter++;
         for (jj=ii+1; jj<(INT4)input->numofcandidates; jj++) {
            if ( usedcandidate->data[jj] == 0 && (input->data[locs->data[0]].fsig-input->data[jj].fsig-0.25/params->Tcoh <= 0.5*iter/params->Tcoh && input->data[locs->data[0]].fsig-input->data[jj].fsig+0.25/params->Tcoh >= 0.5*iter/params->Tcoh) ) {
               locs->data[loc] = jj;
               loc++;
               if (foundany==0) foundany = 1;
            }
         } /* for jj < input->numofcandidates */
      } /* while foundany==1 */

      //Using the list of locations, find the subset that have periods within 1 bin
      //of the second FFT frequencies
      INT4 subsetloc = 0;
      INT4 nextsubsetloc = 0;
      INT4 subsetlocset = 0;
      loc2 = 0;
      for (jj=subsetloc; jj<loc; jj++) {
         if ( usedcandidate->data[locs->data[jj]] == 0 && fabs(params->Tobs/input->data[locs->data[jj]].period - params->Tobs/input->data[locs->data[subsetloc]].period) <= 1.0 ) {
            locs2->data[loc2] = locs->data[jj];
            loc2++;
         } else if (usedcandidate->data[locs->data[jj]] == 0 && subsetlocset == 0) {
            subsetlocset = 1;
            nextsubsetloc = jj;
         }

         if (jj+1 == loc) {
            if (subsetlocset==1) {
               subsetloc = nextsubsetloc;   //Reset subsetloc and jj to the next candidate period
               jj = subsetloc-1;
               subsetlocset = 0;    //Reset the logic of whether there are any more periods to go
            }

            //find best candidate moddepth
            fprintf(stderr,"Finding best modulation depth with number to try %d\n",loc2);
            avefsig = 0.0;
            aveperiod = 0.0;
            mindf = 0.0;
            maxdf = 0.0;
            REAL8 weight = 0.0;
            REAL8 bestmoddepth = 0.0;
            REAL8 bestR = 0.0;
            REAL8 besth0 = 0.0;
            REAL8 bestProb = 0.0;
            INT4 bestproberrcode = 0;
            for (kk=0; kk<loc2; kk++) {
               /* avefsig += input->data[locs2->data[kk]].fsig*(-input->data[locs2->data[kk]].prob);
               aveperiod += input->data[locs2->data[kk]].period*(-input->data[locs2->data[kk]].prob);
               weight += -input->data[locs2->data[kk]].prob; */
               avefsig += input->data[locs2->data[kk]].fsig*(input->data[locs2->data[kk]].prob*input->data[locs2->data[kk]].prob);
               aveperiod += input->data[locs2->data[kk]].period*(input->data[locs2->data[kk]].prob*input->data[locs2->data[kk]].prob);
               weight += input->data[locs2->data[kk]].prob*input->data[locs2->data[kk]].prob;
               if (mindf > input->data[locs2->data[kk]].moddepth || mindf == 0.0) mindf = input->data[locs2->data[kk]].moddepth;
               if (maxdf < input->data[locs2->data[kk]].moddepth) maxdf = input->data[locs2->data[kk]].moddepth;

               if (loc2==1 && input->data[locs2->data[kk]].fsig>=params->fmin && input->data[locs2->data[kk]].fsig<(params->fmin+params->fspan) && input->data[locs2->data[kk]].period>=params->Pmin && input->data[locs2->data[kk]].period<=params->Pmax) {
                  besth0 = input->data[locs2->data[kk]].h0;
                  bestmoddepth = input->data[locs2->data[kk]].moddepth;
                  bestR = input->data[locs2->data[kk]].stat;
                  bestProb = input->data[locs2->data[kk]].prob;
                  bestproberrcode = input->data[locs2->data[kk]].proberrcode;
               }

               usedcandidate->data[locs2->data[kk]] = 1;
            } /* for kk < loc2 */
            avefsig = avefsig/weight;
            aveperiod = aveperiod/weight;

            INT4 proberrcode = 0;

            if (loc2 > 1 && aveperiod >= params->Pmin && aveperiod <= params->Pmax) {
               INT4 numofmoddepths = (INT4)floorf(2*(maxdf-mindf)*params->Tcoh)+1;
               candidate cand;
               templateStruct *template = NULL;
               XLAL_CHECK( (template = new_templateStruct(params->maxtemplatelength)) != NULL, XLAL_EFUNC );

               for (kk=0; kk<numofmoddepths; kk++) {
                  if ((mindf+kk*0.5/params->Tcoh)>=params->dfmin && (mindf+kk*0.5/params->Tcoh)<=params->dfmax) {

                     loadCandidateData(&cand, avefsig, aveperiod, mindf + kk*0.5/params->Tcoh, input->data[0].ra, input->data[0].dec, 0, 0, 0.0, 0, 0.0);

                     if (option==1) {
                        XLAL_CHECK( makeTemplate(template, cand, params, sftexist, plan) == XLAL_SUCCESS, XLAL_EFUNC );
                     } else {
                        XLAL_CHECK( makeTemplateGaussians(template, cand, params, ffdata->numfbins, ffdata->numfprbins) == XLAL_SUCCESS, XLAL_EFUNC );
                     }

                     REAL8 R = calculateR(ffdata->ffdata, template, ffplanenoise, fbinaveratios);
                     XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
                     REAL8 prob = probR(template, ffplanenoise, fbinaveratios, R, params, &proberrcode);
                     XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
                     //REAL8 h0 = 2.7426*pow(R/(params->Tcoh*params->Tobs),0.25);

                     if (prob < bestProb) {
                        //besth0 = h0;
                        bestmoddepth = mindf + kk*0.5/params->Tcoh;
                        bestR = R;
                        bestProb = prob;
                        bestproberrcode = proberrcode;
                     }
                  } /* if test moddepth is within user specified range */
               } /* for kk < numofmoddepths */

               free_templateStruct(template);
               template = NULL;
            } /* if loc2 > 1 ... */

            if (bestProb != 0.0) {
               if (bestR > 0.0) besth0 = 2.7426*pow(bestR/(params->Tcoh*params->Tobs),0.25);
               else besth0 = 0.0;

               if ((*output)->numofcandidates == (*output)->length-1) {
                  XLAL_CHECK( (*output = resize_candidateVector(*output, 2*(*output)->length)) != NULL, XLAL_EFUNC );
               }
               loadCandidateData(&((*output)->data[(*output)->numofcandidates]), avefsig, aveperiod, bestmoddepth, input->data[0].ra, input->data[0].dec, bestR, besth0, bestProb, bestproberrcode, input->data[0].normalization);
               numcandoutlist++;
               (*output)->numofcandidates++;
            }

            loc2 = 0;
         } /* if jj+1 == loc */
      } /* for jj < loc */

      //Find location of first entry to be searched next time or finish the cluster search
      for (jj=ii; jj<(INT4)input->numofcandidates; jj++) {
         if (usedcandidate->data[jj]==0) {
            ii = jj - 1;
            jj = (INT4)input->numofcandidates - 1;
         } else if (jj==(INT4)input->numofcandidates-1) {
            ii = (INT4)input->numofcandidates - 1;
         }
      }

      //Reinitialize values, just in case
      for (jj=0; jj<(INT4)locs->length; jj++) {
         locs->data[jj] = -1;
         locs2->data[jj] = -1;
      }
   } /* for ii < numofcandidates */

   //Destroy stuff
   XLALDestroyINT4Vector(locs);
   XLALDestroyINT4Vector(locs2);
   XLALDestroyINT4Vector(usedcandidate);
   if (option==1) XLALDestroyREAL4FFTPlan(plan);

   return 0;

} /* clusterCandidates() */



//Big function to test the IHS candidates against Gaussian templates
INT4 testIHScandidates(candidateVector **output, candidateVector *ihsCandidates, ffdataStruct *ffdata, REAL4Vector *aveNoise, REAL4Vector *aveTFnoisePerFbinRatio, REAL4 alpha, REAL4 delta, inputParamsStruct *inputParams)
{

   //R probability calculator errorcode
   INT4 proberrcode = 0;

   INT4 ii, jj, kk;

   //Allocate memory for FAR struct
   farStruct *farval = NULL;
   XLAL_CHECK( (farval = new_farStruct()) != NULL, XLAL_EFUNC );

   //Allocate memory for template
   templateStruct *template = NULL;
   XLAL_CHECK( (template = new_templateStruct(inputParams->maxtemplatelength)) != NULL, XLAL_EFUNC );

   INT4 candidatesoutsideofmainULrange = 0;
   REAL8 log10templatefar = inputParams->log10templatefar;

   for (ii=0; ii<(INT4)ihsCandidates->numofcandidates; ii++) {
      //Assess the IHS candidate if the signal is away from the band edges, the modulation depth is greater or equal to minimum allowed and less than or equal to the maximum allowed, and if the period/modulation depth combo is within allowable limits for a template to be made. We will cut the period space in the next step.
      if ( ihsCandidates->data[ii].fsig>=inputParams->fmin && ihsCandidates->data[ii].fsig<(inputParams->fmin+inputParams->fspan) ) {
         if ( inputParams->followUpOutsideULrange || (ihsCandidates->data[ii].fsig>=inputParams->ULfmin && ihsCandidates->data[ii].fsig<=(inputParams->ULfmin + inputParams->ULfspan) && ihsCandidates->data[ii].moddepth>=inputParams->ULmindf && ihsCandidates->data[ii].moddepth<=inputParams->ULmaxdf) ) {

            resetTemplateStruct(template);

            REAL8 R, prob, bestPeriod = 0.0, bestR = 0.0, bestProb = 0.0;
            INT4 bestproberrcode = 0;

            if (ihsCandidates->data[ii].period>=fmax(2.0*3600.0, minPeriod(ihsCandidates->data[ii].moddepth, inputParams->Tcoh)) && ihsCandidates->data[ii].period<=(0.2*inputParams->Tobs)) {
               //Make a Gaussian train template
               XLAL_CHECK( makeTemplateGaussians(template, ihsCandidates->data[ii], inputParams, ffdata->numfbins, ffdata->numfprbins) == XLAL_SUCCESS, XLAL_EFUNC );

               //remove this--for testing purposes only
               //for (jj=0; jj<(INT4)template->templatedata->length; jj++) fprintf(stderr, "%g %d %d %d %g\n", template->templatedata->data[jj], template->pixellocations->data[jj], template->firstfftfrequenciesofpixels->data[jj], template->secondfftfrequencies->data[jj], aveNoise->data[template->secondfftfrequencies->data[jj]]*aveTFnoisePerFbinRatio->data[template->firstfftfrequenciesofpixels->data[jj]]);
               /* for (jj=0; jj<50; jj++) {
                  REAL8 probval = probR(template, aveNoise, aveTFnoisePerFbinRatio, 0.3*jj-2.0, inputParams, &proberrcode);
                  fprintf(stderr, "%f %g\n", 0.3*jj-2.0, pow(10.0, probval));
               } */
               /* resetTemplateStruct(template);
               REAL4FFTPlan *FFTplan = XLALCreateForwardREAL4FFTPlan(ffdata->numffts, inputParams->FFTplanFlag);
               INT4Vector *sftexist = XLALCreateINT4Vector(ffdata->numffts);
               for (jj=0; jj<(INT4)ffdata->numffts; jj++) sftexist->data[jj] = 1;
               makeTemplate(template, ihsCandidates->data[ii], inputParams, sftexist, FFTplan);
               if (xlalErrno!=0) {
                  fprintf(stderr,"%s: makeTemplate() failed.\n", __func__);
                  XLAL_ERROR(XLAL_EFUNC);
               }
               fprintf(stderr, "\n"); */
               //for (jj=0; jj<(INT4)template->templatedata->length; jj++) fprintf(stderr, "%g %d %d %d %g\n", template->templatedata->data[jj], template->pixellocations->data[jj], template->firstfftfrequenciesofpixels->data[jj], template->secondfftfrequencies->data[jj], aveNoise->data[template->secondfftfrequencies->data[jj]]*aveTFnoisePerFbinRatio->data[template->firstfftfrequenciesofpixels->data[jj]]);
               /* for (jj=0; jj<50; jj++) {
                  REAL8 probval = probR(template, aveNoise, aveTFnoisePerFbinRatio, 0.75*jj-8.0, inputParams, &proberrcode);
                  fprintf(stderr, "%f %g\n", 0.75*jj-8.0, pow(10.0, probval));
               } */
               /* XLALDestroyREAL4FFTPlan(FFTplan);
               XLALDestroyINT4Vector(sftexist); */

               //Estimate the FAR for these bin weights if the option was given
               if (inputParams->calcRthreshold) {
                  XLAL_CHECK( numericFAR(farval, template, inputParams->templatefar, aveNoise, aveTFnoisePerFbinRatio, inputParams, inputParams->rootFindingMethod) == XLAL_SUCCESS, XLAL_EFUNC );
               }

               //Caclulate R, probability noise caused the candidate, and estimate of h0
               R = calculateR(ffdata->ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
               XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
               prob = probR(template, aveNoise, aveTFnoisePerFbinRatio, R, inputParams, &proberrcode);
               XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );

               /* Note the candidate if R exceeds the FAR or check other possibilities of different
                periods */
               if ((!inputParams->calcRthreshold && prob<log10templatefar) || (inputParams->calcRthreshold && R>farval->far)) {
                  bestR = R;
                  bestProb = prob;
                  bestPeriod = ihsCandidates->data[ii].period;
               } /* if prob<log10templatefar || R > farval->far */
            } // if within moddepth/period range

            // longer or shorter
            REAL8 periodfact = 0.0;
            for (jj=0; jj<=1; jj++) {
               //Shift by harmonics
               for (kk=1; kk<=inputParams->periodHarmToCheck; kk++) {
                  if (jj==0) periodfact = 1.0/(REAL8)kk;
                  else periodfact = (REAL8)kk;
                  if (ihsCandidates->data[ii].period*periodfact>=fmax(inputParams->Pmin, minPeriod(ihsCandidates->data[ii].moddepth, inputParams->Tcoh)) && ihsCandidates->data[ii].period*periodfact<=fmin(inputParams->Pmax, inputParams->Tobs*0.2)) {
                     ihsCandidates->data[ii].period *= periodfact;
                     XLAL_CHECK( makeTemplateGaussians(template, ihsCandidates->data[ii], inputParams, ffdata->numfbins, ffdata->numfprbins) == XLAL_SUCCESS, XLAL_EFUNC );
                     R = calculateR(ffdata->ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
                     XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
                     prob = probR(template, aveNoise, aveTFnoisePerFbinRatio, R, inputParams, &proberrcode);
                     XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
                     if (inputParams->calcRthreshold && bestProb==0.0) {
                        XLAL_CHECK( numericFAR(farval, template, inputParams->templatefar, aveNoise, aveTFnoisePerFbinRatio, inputParams, inputParams->rootFindingMethod) == XLAL_SUCCESS, XLAL_EFUNC );
                     }
                     if ((bestProb!=0.0 && prob<bestProb) || (bestProb==0.0 && !inputParams->calcRthreshold && prob<log10templatefar) || (bestProb==0.0 && inputParams->calcRthreshold && R>farval->far)) {
                        bestPeriod = ihsCandidates->data[ii].period;
                        bestR = R;
                        bestProb = prob;
                        bestproberrcode = proberrcode;
                     }
                     ihsCandidates->data[ii].period /= periodfact;  //reset the period back to the original value
                  } // in range?
               } // shift by harmonics for kk <= inputParams->periodHarmToCheck (harmonics)

               //shift by fractions
               for (kk=1; kk<=inputParams->periodFracToCheck; kk++) {
                  if (jj==0) periodfact = (kk+1.0)/(kk+2.0);
                  else periodfact = (kk+2.0)/(kk+1.0);
                  if (ihsCandidates->data[ii].period*periodfact>=fmax(inputParams->Pmin, minPeriod(ihsCandidates->data[ii].moddepth, inputParams->Tcoh)) && ihsCandidates->data[ii].period*periodfact<=fmin(inputParams->Pmax, inputParams->Tobs*0.2)) {
                     ihsCandidates->data[ii].period *= periodfact;
                     XLAL_CHECK( makeTemplateGaussians(template, ihsCandidates->data[ii], inputParams, ffdata->numfbins, ffdata->numfprbins) == XLAL_SUCCESS, XLAL_EFUNC );
                     R = calculateR(ffdata->ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
                     XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
                     prob = probR(template, aveNoise, aveTFnoisePerFbinRatio, R, inputParams, &proberrcode);
                     XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
                     if (inputParams->calcRthreshold && bestProb==0.0) {
                        XLAL_CHECK( numericFAR(farval, template, inputParams->templatefar, aveNoise, aveTFnoisePerFbinRatio, inputParams, inputParams->rootFindingMethod) == XLAL_SUCCESS, XLAL_EFUNC );
                     }
                     if ((bestProb!=0.0 && prob<bestProb) || (bestProb==0.0 && !inputParams->calcRthreshold && prob<log10templatefar) || (bestProb==0.0 && inputParams->calcRthreshold && R>farval->far)) {
                        bestPeriod = ihsCandidates->data[ii].period;
                        bestR = R;
                        bestProb = prob;
                        bestproberrcode = proberrcode;
                     }
                     ihsCandidates->data[ii].period /= periodfact;  //reset the period back to the original value
                  } // in range?
               } // shift by fractions kk <= inputParams->periodFracToCheck
            } // longer or shorter

            if (bestProb != 0.0) {
               REAL8 h0 = 0.0;
               if (bestR > 0.0) h0 = 2.7426*sqrt(sqrt(bestR/(inputParams->Tcoh*inputParams->Tobs)));  //Now compute the h0 value

               if ((*output)->numofcandidates == (*output)->length-1) {
                  XLAL_CHECK( (*output = resize_candidateVector(*output, 2*(*output)->length)) != NULL, XLAL_EFUNC );
               }
               loadCandidateData(&((*output)->data[(*output)->numofcandidates]), ihsCandidates->data[ii].fsig, bestPeriod, ihsCandidates->data[ii].moddepth, alpha, delta, bestR, h0, bestProb, bestproberrcode, ihsCandidates->data[ii].normalization);
               (*output)->numofcandidates++;

            } /* if bestR != 0.0, add candidate or replace if something better is found */
         } /* if within UL boundaries */
         else {
            candidatesoutsideofmainULrange++;
         }
      } /* if within outer boundaries */
   } /* for ii < numofcandidates */

   fprintf(stderr, "%d remaining candidate(s) inside UL range.\n", ihsCandidates->numofcandidates-candidatesoutsideofmainULrange);

   //Destroy allocated memory
   free_templateStruct(template);
   template = NULL;
   free_farStruct(farval);
   farval = NULL;

   return 0;

} /* testIHScandidates() */


//Keep the most significant candidates, potentially reducing the number of candidates if there are more than allowed
candidateVector * keepMostSignificantCandidates(candidateVector *input, inputParamsStruct *params)
{

   INT4 ii, jj;
   candidateVector *output = NULL;

   //If the number to keep is > 0 and the number of candidates is less than the number to keep,
   //just move the input vector to the output vector
   if (params->keepOnlyTopNumIHS>0 && (INT4)input->numofcandidates<=params->keepOnlyTopNumIHS) {
      XLAL_CHECK_NULL( (output = new_candidateVector(input->numofcandidates)) != NULL, XLAL_EFUNC );

      for (ii=0; ii<(INT4)input->numofcandidates; ii++) {
         loadCandidateData(&(output->data[ii]), input->data[ii].fsig, input->data[ii].period, input->data[ii].moddepth, input->data[ii].ra, input->data[ii].dec, input->data[ii].stat, input->data[ii].h0, input->data[ii].prob, input->data[ii].proberrcode, input->data[ii].normalization);
      }
      output->numofcandidates = input->numofcandidates;

   } else if (params->keepOnlyTopNumIHS>0 && (INT4)input->numofcandidates>params->keepOnlyTopNumIHS) {
      //If keep is > 0 and the number of candidates is > the number to keep,
      //we sort through the list and find the most significant candidates to keep
      XLAL_CHECK_NULL( (output = new_candidateVector(params->keepOnlyTopNumIHS)) != NULL, XLAL_EFUNC );

      for (ii=0; ii<(INT4)output->length; ii++) {
         REAL8 highestsignificance = 0.0;
         INT4 candidateWithHighestSignificance = 0;
         for (jj=0; jj<(INT4)input->numofcandidates; jj++) {
            if (input->data[jj].prob>highestsignificance) {
               highestsignificance = input->data[jj].prob;
               candidateWithHighestSignificance = jj;
            }
         }

         loadCandidateData(&(output->data[ii]), input->data[candidateWithHighestSignificance].fsig, input->data[candidateWithHighestSignificance].period, input->data[candidateWithHighestSignificance].moddepth, input->data[candidateWithHighestSignificance].ra, input->data[candidateWithHighestSignificance].dec, input->data[candidateWithHighestSignificance].stat, input->data[candidateWithHighestSignificance].h0, input->data[candidateWithHighestSignificance].prob, input->data[candidateWithHighestSignificance].proberrcode, input->data[candidateWithHighestSignificance].normalization);

         input->data[candidateWithHighestSignificance].prob = 0.0;

      }
      output->numofcandidates = params->keepOnlyTopNumIHS;

   } else {
      //Otherwise, we need to fail
      fprintf(stderr, "%s: keepOnlyTopNumIHS given (%d) is not greater than 0, but it should be to use this function.\n", __func__, params->keepOnlyTopNumIHS);
      XLAL_ERROR_NULL(XLAL_EINVAL);
   }

   return output;

} /* keepMostSignificantCandidates() */


//////////////////////////////////////////////////////////////
// Calculate the R statistic
REAL8 calculateR(REAL4Vector *ffdata, templateStruct *templatestruct, REAL4Vector *noise, REAL4Vector *fbinaveratios)
{

   INT4 ii;

   REAL8 sumofsqweights = 0.0;
   for (ii=0; ii<(INT4)templatestruct->templatedata->length; ii++) if (templatestruct->templatedata->data[ii]!=0.0) sumofsqweights += (templatestruct->templatedata->data[ii]*templatestruct->templatedata->data[ii]);
   XLAL_CHECK_REAL8( sumofsqweights != 0.0, XLAL_EFPDIV0 );
   REAL8 sumofsqweightsinv = 1.0/sumofsqweights;

   REAL8 R = 0.0;
   for (ii=0; ii<(INT4)templatestruct->templatedata->length; ii++) {
      if (templatestruct->templatedata->data[ii]!=0.0) {
         R += (ffdata->data[ templatestruct->pixellocations->data[ii] ] - noise->data[ templatestruct->secondfftfrequencies->data[ii] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[ii] ])*templatestruct->templatedata->data[ii]*sumofsqweightsinv;
      }
   }

   return R;

} /* calculateR() */


//////////////////////////////////////////////////////////////
// Calculates maximum modulation depth allowed
REAL8 maxModDepth(REAL8 period, REAL8 cohtime)
{

   REAL8 maxB = 0.5*period/(cohtime*cohtime);

   return maxB;

} /* maxModDepth() */


//////////////////////////////////////////////////////////////
// Calculates minimum period allowable for modulation depth and Tcoh
REAL8 minPeriod(REAL8 moddepth, REAL8 cohtime)
{

   REAL8 minP = 2.0*moddepth*cohtime*cohtime;

   return minP;

} /* minPeriod() */

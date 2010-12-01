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

#include <math.h>
#include <lal/LALMalloc.h>
#include "candidates.h"
#include "templates.h"


candidateVector * new_candidateVector(UINT4 length)
{
   
   const CHAR *fn = __func__;
   
   candidateVector *vector = XLALMalloc(sizeof(*vector));
   if (vector==NULL) {
      fprintf(stderr,"%s: XLALMalloc(%lu) failed.\n", fn, sizeof(*vector));
      XLAL_ERROR_NULL(fn, XLAL_ENOMEM);
   }
      
   vector->length = length;
   vector->numofcandidates = 0;
   if (length==0) vector->data = NULL;
   else {
      vector->data = XLALMalloc( length*sizeof(*vector->data) );
      if (vector->data==NULL) {
         XLALFree((candidateVector*)vector);
         fprintf(stderr,"%s: XLALMalloc(%lu) failed.\n", fn, length*sizeof(*vector->data));
         XLAL_ERROR_NULL(fn, XLAL_ENOMEM);
      }
   }
   
   return vector;
   
} /* new_candidateVector() */



candidateVector * resize_candidateVector(candidateVector *vector, UINT4 length)
{
   
   const CHAR *fn = __func__;
   
   if (vector==NULL) return new_candidateVector(length);
   if (length==0) {
      free_candidateVector(vector);
      return NULL;
   }
   
   vector->data = XLALRealloc(vector->data, length*sizeof(*vector->data));
   if (vector->data==NULL) {
      vector->length = 0;
      fprintf(stderr,"%s: XLALRealloc(%lu) failed.\n", fn, length*sizeof(*vector->data));
      XLAL_ERROR_NULL(fn, XLAL_ENOMEM);
   }
   vector->length = length;
   
   return vector;
   
} /* resize_candidateVector() */



void free_candidateVector(candidateVector *vector)
{
   
   const CHAR *fn = __func__;
   
   if (vector==NULL) return;
   if ((!vector->length || !vector->data) && (vector->length || vector->data)) XLAL_ERROR_VOID(fn, XLAL_EINVAL);
   if (vector->data) XLALFree((candidate*)vector->data);
   vector->data = NULL;
   XLALFree((candidateVector*)vector);
   return;
   
} /* delete_candidateVector() */



//////////////////////////////////////////////////////////////
// Load candidate data
void loadCandidateData(candidate* output, REAL8 fsig, REAL8 period, REAL8 moddepth, REAL4 ra, REAL4 dec, REAL8 stat, REAL8 h0, REAL8 prob, INT4 proberrcode, REAL8 normalization)
{
   
   output->fsig = fsig;
   output->period = period;
   output->moddepth = moddepth;
   output->ra = ra;
   output->dec = dec;
   output->stat = stat;
   output->h0 = h0;
   output->prob = prob;
   output->proberrcode = proberrcode;
   output->normalization = normalization;

} /* loadCandidateData() */


//////////////////////////////////////////////////////////////
// Cluster candidates by frequency and period using templates:
// option = 0 uses Gaussian templates (default)
// option = 1 uses exact templates
void clusterCandidates(candidateVector *output, candidateVector *input, ffdataStruct *ffdata, inputParamsStruct *params, REAL4Vector *ffplanenoise, REAL4Vector *fbinaveratios, INT4 option)
{
   
   const char *fn = __func__;
   
   INT4 ii, jj, kk, loc, loc2, numcandoutlist;
   REAL8 avefsig, aveperiod, mindf, maxdf;
   
   //Allocate int vectors for storage
   INT4Vector *locs = XLALCreateINT4Vector((UINT4)input->numofcandidates);
   INT4Vector *locs2 = XLALCreateINT4Vector((UINT4)input->numofcandidates);
   INT4Vector *usedcandidate = XLALCreateINT4Vector((UINT4)input->numofcandidates);
   if (locs==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", fn, input->numofcandidates);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (locs2==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", fn, input->numofcandidates);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (usedcandidate==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", fn, input->numofcandidates);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   
   //Initialize arrays
   for (ii=0; ii<(INT4)input->numofcandidates; ii++) {
      locs->data[ii] = -1;
      locs2->data[ii] = -1;
      usedcandidate->data[ii] = 0;
   }
   
   //If not equal, make output vector the same length as the input vector
   if (output->length < input->length) {
      output = resize_candidateVector(output, input->length);
      if (output->length < input->length) {
         fprintf(stderr,"%s: resize_candidateVector(%d) failed.\n", fn, input->length);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
   }
   
   //Set default if bad option given
   if (option!=0 || option!=1) option = 0;
   
   //Make FFT plan if option 1 is given
   REAL4FFTPlan *plan = NULL;
   if (option==1) {
      plan = XLALCreateForwardREAL4FFTPlan((UINT4)floor(2*(params->Tobs/params->Tcoh)-1), 0);
      if (plan==NULL) {
         fprintf(stderr,"%s: XLALCreateForwardREAL4FFTPlan(%d, 0) failed.\n", fn, (INT4)floor(2*(params->Tobs/params->Tcoh)-1));
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
   }
   
   numcandoutlist = 0;
   for (ii=0; ii<(INT4)input->numofcandidates; ii++) {
      
      if (input->data[ii].period==0.0) exit(-1);
      
      //Make note of first candidate available
      locs->data[0] = ii;
      loc = 1;
      
      INT4 foundany = 0;   //Switch to determing if any other candidates in the group. 1 if true
      INT4 iter = 1;
      //Find any in the list that are within +1/2 bin in first FFT frequency
      for (jj=ii+1; jj<(INT4)input->numofcandidates; jj++) {
         if ( usedcandidate->data[jj] == 0 && (input->data[jj].fsig-input->data[locs->data[0]].fsig <= 0.5*iter/params->Tcoh+1e-6 && input->data[jj].fsig-input->data[locs->data[0]].fsig >= -0.25*iter/params->Tcoh) ) {
            locs->data[loc] = jj;
            loc++;
            if (foundany==0) foundany = 1;
         }
      }
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
         }
      }
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
         }
      }
      
      
      
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
            REAL8 bestProb = 1.0;
            INT4 bestproberrcode = 0;
            for (kk=0; kk<loc2; kk++) {
               avefsig += input->data[locs2->data[kk]].fsig*input->data[locs2->data[kk]].stat;
               aveperiod += input->data[locs2->data[kk]].period*input->data[locs2->data[kk]].stat;
               weight += input->data[locs2->data[kk]].stat;
               if (mindf > input->data[locs2->data[kk]].moddepth || mindf == 0.0) mindf = input->data[locs2->data[kk]].moddepth;
               if (maxdf < input->data[locs2->data[kk]].moddepth) maxdf = input->data[locs2->data[kk]].moddepth;
               
               if (loc2==1 && aveperiod/weight >= params->Pmin && aveperiod/weight <= params->Pmax) {
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
               for (kk=0; kk<numofmoddepths; kk++) {
                  
                  candidate cand;
                  loadCandidateData(&cand, avefsig, aveperiod, mindf + kk*0.5/params->Tcoh, input->data[0].ra, input->data[0].dec, 0, 0, 0.0, 0, 0.0);
                  
                  templateStruct *template = new_templateStruct(params->templatelength);
                  if (template==NULL) {
                     fprintf(stderr,"%s: new_templateStruct(%d) failed.\n", fn, params->templatelength);
                     XLAL_ERROR_VOID(fn, XLAL_EFUNC);
                  }
                  if (option==1) {
                     makeTemplate(template, cand, params, plan);
                     if (xlalErrno!=0) {
                        fprintf(stderr,"%s: makeTemplate() failed.\n", fn);
                        XLAL_ERROR_VOID(fn, XLAL_EFUNC);
                     }
                  } else {
                     makeTemplateGaussians(template, cand, params);
                     if (xlalErrno!=0) {
                        fprintf(stderr,"%s: makeTemplateGaussians() failed.\n", fn);
                        XLAL_ERROR_VOID(fn, XLAL_EFUNC);
                     }
                  }
                  
                  farStruct *farval = new_farStruct();
                  if (farval==NULL) {
                     fprintf(stderr,"%s: new_farStruct() failed.\n", fn);
                     XLAL_ERROR_VOID(fn, XLAL_EFUNC);
                  }
                  numericFAR(farval, template, 0.01, ffplanenoise, fbinaveratios);
                  if (xlalErrno!=0) {
                     fprintf(stderr,"%s: numericFAR() failed.\n", fn);
                     XLAL_ERROR_VOID(fn, XLAL_EFUNC);
                  }
                  
                  REAL8 R = calculateR(ffdata->ffdata, template, ffplanenoise, fbinaveratios);
                  REAL8 prob = probR(template, ffplanenoise, fbinaveratios, R, &proberrcode);
                  if (XLAL_IS_REAL8_FAIL_NAN(R)) {
                     fprintf(stderr,"%s: calculateR() failed.\n", fn);
                     XLAL_ERROR_VOID(fn, XLAL_EFUNC);
                  } else if (XLAL_IS_REAL8_FAIL_NAN(prob)) {
                     fprintf(stderr,"%s: probR() failed.\n", fn);
                     XLAL_ERROR_VOID(fn, XLAL_EFUNC);
                  }
                  REAL8 h0 = 2.9569*pow(R/(params->Tcoh*params->Tobs),0.25);
                  
                  if (R > farval->far && prob < bestProb) {
                     besth0 = h0;
                     bestmoddepth = mindf + kk*0.5/params->Tcoh;
                     bestR = R;
                     bestProb = prob;
                     bestproberrcode = proberrcode;
                  }
                  
                  free_templateStruct(template);
                  template = NULL;
                  free_farStruct(farval);
                  farval = NULL;
               } /* for kk < numofmoddepths */
            }
            
            if (bestR != 0.0) {
               loadCandidateData(&output->data[output->numofcandidates], avefsig, aveperiod, bestmoddepth, input->data[0].ra, input->data[0].dec, bestR, besth0, bestProb, bestproberrcode, input->data[0].normalization);
               numcandoutlist++;
               (output->numofcandidates)++;
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
   
   XLALDestroyINT4Vector(locs);
   XLALDestroyINT4Vector(locs2);
   XLALDestroyINT4Vector(usedcandidate);
   if (option==1) XLALDestroyREAL4FFTPlan(plan);
   
} /* clusterCandidates() */







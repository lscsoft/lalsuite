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


//////////////////////////////////////////////////////////////
// Create memory for new candidate
candidate * new_candidate(void)
{

   candidate *cand = (candidate*)XLALMalloc(sizeof(candidate));
   
   cand->fsig = 0.0;
   cand->period = 0.0;
   cand->moddepth = 0.0;
   cand->ra = 0.0;
   cand->dec = 0.0;
   cand->stat = 0.0;
   cand->h0 = 0.0;
   cand->prob = 0.0;
   cand->proberrcode = 0;
   cand->normalization = 0.0;
   
   return cand;

}


//////////////////////////////////////////////////////////////
// Free memory of a candidate
void free_candidate(candidate *cand)
{
   
   cand->fsig = 0.0;
   cand->period = 0.0;
   cand->moddepth = 0.0;
   cand->ra = 0.0;
   cand->dec = 0.0;
   cand->stat = 0.0;
   cand->h0 = 0.0;
   cand->prob = 0.0;
   cand->proberrcode = 0;
   cand->normalization = 0.0;
   
   XLALFree((candidate*)cand);

}


//////////////////////////////////////////////////////////////
// Load candidate data
void loadCandidateData(candidate *output, REAL8 fsig, REAL8 period, REAL8 moddepth, REAL4 ra, REAL4 dec, REAL8 stat, REAL8 h0, REAL8 prob, INT4 proberrcode, REAL8 normalization)
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

}


//////////////////////////////////////////////////////////////
// Cluster candidates by frequency and period using templates:
// option = 0 uses Gaussian templates (default)
// option = 1 uses exact templates
void clusterCandidates(candidate *output[], candidate *input[], ffdataStruct *ffdata, inputParamsStruct *params, REAL4Vector *ffplanenoise, REAL4Vector *fbinaveratios, INT4 numofcandidates, INT4 option)
{

   INT4 ii, jj, kk, loc, loc2, numcandoutlist;
   REAL8 avefsig, aveperiod, mindf, maxdf;
   
   INT4Vector *locs = XLALCreateINT4Vector((UINT4)numofcandidates);
   INT4Vector *locs2 = XLALCreateINT4Vector((UINT4)numofcandidates);
   INT4Vector *usedcandidate = XLALCreateINT4Vector((UINT4)numofcandidates);
   for (ii=0; ii<numofcandidates; ii++) {
      locs->data[ii] = -1;
      locs2->data[ii] = -1;
      usedcandidate->data[ii] = 0;
   }
   
   //Set default if bad option given
   if (option!=0 || option!=1) option = 0;
   
   //Make FFT plan if option 1 is given
   REAL4FFTPlan *plan = NULL;
   if (option==1) plan = XLALCreateForwardREAL4FFTPlan((UINT4)floor(2*(params->Tobs/params->Tcoh)-1), 0);
   
   numcandoutlist = 0;
   for (ii=0; ii<numofcandidates; ii++) {
      
      if (input[ii]->period==0.0) exit(-1);
      
      //Make note of first candidate available
      locs->data[0] = ii;
      loc = 1;
      
      INT4 foundany = 0;   //Switch to determing if any other candidates in the group. 1 if true
      INT4 iter = 1;
      //Find any in the list that are within +1/2 bin in first FFT frequency
      for (jj=ii+1; jj<numofcandidates; jj++) {
         if ( usedcandidate->data[jj] == 0 && (input[jj]->fsig-input[locs->data[0]]->fsig <= 0.5*iter/params->Tcoh+1e-6 && input[jj]->fsig-input[locs->data[0]]->fsig >= -0.25*iter/params->Tcoh) ) {
            locs->data[loc] = jj;
            loc++;
            if (foundany==0) foundany = 1;
         }
      }
      //Keep checking as long as there are more connected frequencies going higher in frequency
      while (foundany==1) {
         foundany = 0;
         iter++;
         for (jj=ii+1; jj<numofcandidates; jj++) {
            if ( usedcandidate->data[jj] == 0 && (input[jj]->fsig-input[locs->data[0]]->fsig-0.25/params->Tcoh <= 0.5*iter/params->Tcoh && input[jj]->fsig-input[locs->data[0]]->fsig+0.25/params->Tcoh >= 0.5*iter/params->Tcoh) ) {
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
         for (jj=ii+1; jj<numofcandidates; jj++) {
            if ( usedcandidate->data[jj] == 0 && (input[locs->data[0]]->fsig-input[jj]->fsig-0.25/params->Tcoh <= 0.5*iter/params->Tcoh && input[locs->data[0]]->fsig-input[jj]->fsig+0.25/params->Tcoh >= 0.5*iter/params->Tcoh) ) {
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
         if ( usedcandidate->data[locs->data[jj]] == 0 && fabs(params->Tobs/input[locs->data[jj]]->period - params->Tobs/input[locs->data[subsetloc]]->period) <= 1.0 ) {
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
               avefsig += input[locs2->data[kk]]->fsig*input[locs2->data[kk]]->stat;
               aveperiod += input[locs2->data[kk]]->period*input[locs2->data[kk]]->stat;
               weight += input[locs2->data[kk]]->stat;
               if (mindf > input[locs2->data[kk]]->moddepth || mindf == 0.0) mindf = input[locs2->data[kk]]->moddepth;
               if (maxdf < input[locs2->data[kk]]->moddepth) maxdf = input[locs2->data[kk]]->moddepth;
               
               if (loc2==1 && aveperiod/weight >= params->Pmin && aveperiod/weight <= params->Pmax) {
                  besth0 = input[locs2->data[kk]]->h0;
                  bestmoddepth = input[locs2->data[kk]]->moddepth;
                  bestR = input[locs2->data[kk]]->stat;
                  bestProb = input[locs2->data[kk]]->prob;
                  bestproberrcode = input[locs2->data[kk]]->proberrcode;
               }
               
               usedcandidate->data[locs2->data[kk]] = 1;
            }
            avefsig = avefsig/weight;
            aveperiod = aveperiod/weight;
            
            INT4 proberrcode = 0;
            
            if (loc2 > 1 && aveperiod >= params->Pmin && aveperiod <= params->Pmax) {
               INT4 numofmoddepths = (INT4)floorf(2*(maxdf-mindf)*params->Tcoh)+1;
               for (kk=0; kk<numofmoddepths; kk++) {
                  
                  candidate *cand = new_candidate();
                  loadCandidateData(cand, avefsig, aveperiod, mindf + kk*0.5/params->Tcoh, input[0]->ra, input[0]->dec, 0, 0, 0.0, 0, 0.0);
                  templateStruct *template = new_templateStruct(params->templatelength);
                  if (option==1) makeTemplate(template, cand, params, plan);
                  else makeTemplateGaussians(template, cand, params);
                  farStruct *farval = new_farStruct();
                  numericFAR(farval, template, 0.01, ffplanenoise, fbinaveratios);
                  REAL8 R = calculateR(ffdata->ffdata, template, ffplanenoise, fbinaveratios);
                  REAL8 prob = probR(template, ffplanenoise, fbinaveratios, R, &proberrcode);
                  REAL8 h0 = 2.9569*pow(R/(params->Tcoh*params->Tobs),0.25);
                  if (R > farval->far && prob < bestProb) {
                     besth0 = h0;
                     bestmoddepth = mindf + kk*0.5/params->Tcoh;
                     bestR = R;
                     bestProb = prob;
                     bestproberrcode = proberrcode;
                  }
                  free_candidate(cand);
                  cand = NULL;
                  free_templateStruct(template);
                  template = NULL;
                  free_farStruct(farval);
                  farval = NULL;
               }
            }
            
            if (bestR != 0.0) {
               output[numcandoutlist] = new_candidate();
               loadCandidateData(output[numcandoutlist], avefsig, aveperiod, bestmoddepth, input[0]->ra, input[0]->dec, bestR, besth0, bestProb, bestproberrcode, input[0]->normalization);
               numcandoutlist++;
            }
            
            loc2 = 0;
         }
      }
      
      //Find location of first entry to be searched next time or finish the cluster search
      for (jj=ii; jj<numofcandidates; jj++) {
         if (usedcandidate->data[jj]==0) {
            ii = jj - 1;
            jj = numofcandidates - 1;
         } else if (jj==numofcandidates-1) {
            ii = numofcandidates - 1;
         }
      }
      
      //Reinitialize values, just in case
      for (jj=0; jj<(INT4)locs->length; jj++) {
         locs->data[jj] = -1;
         locs2->data[jj] = -1;
      }
   }
   
   XLALDestroyINT4Vector(locs);
   XLALDestroyINT4Vector(locs2);
   XLALDestroyINT4Vector(usedcandidate);
   if (option==1) XLALDestroyREAL4FFTPlan(plan);
   
}







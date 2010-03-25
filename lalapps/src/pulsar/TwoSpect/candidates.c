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
#include <lal/RealFFT.h>
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
   cand->snr = 0.0;
   cand->prob = 0.0;
   
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
   cand->snr = 0.0;
   cand->prob = 0.0;
   
   XLALFree((candidate*)cand);

}


//////////////////////////////////////////////////////////////
// Load candidate data
void loadCandidateData(candidate *out, REAL4 fsig, REAL4 period, REAL4 moddepth, REAL4 ra, REAL4 dec, REAL4 stat, REAL4 snr, REAL4 prob)
{

   out->fsig = fsig;
   out->period = period;
   out->moddepth = moddepth;
   out->ra = ra;
   out->dec = dec;
   out->stat = stat;
   out->snr = snr;
   out->prob = prob;

}


//////////////////////////////////////////////////////////////
// Cluster candidates by frequency and period using templates:
// option = 0 uses Gaussian templates (default)
// option = 1 uses exact templates
void clusterCandidates(candidate *out[], candidate *in[], ffdataStruct *ffdata, inputParamsStruct *params, REAL4Vector *ffplanenoise, INT4 numofcandidates, INT4 option)
{

   INT4 ii, jj, kk, loc, loc2, numcandoutlist;
   REAL4 avefsig, aveperiod, mindf, maxdf;
   
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
      
      if (in[ii]->period==0.0) exit(-1);
      
      //Make note of first candidate available
      locs->data[0] = ii;
      loc = 1;
      
      INT4 foundany = 0;   //Switch to determing if any other candidates in the group. 1 if true
      INT4 iter = 1;
      //Find any in the list that are within +1/2 bin in first FFT frequency
      for (jj=ii+1; jj<numofcandidates; jj++) {
         if ( usedcandidate->data[jj] == 0 && (in[jj]->fsig-in[locs->data[0]]->fsig <= 0.5*iter/params->Tcoh+1e-6 && in[jj]->fsig-in[locs->data[0]]->fsig >= -0.25*iter/params->Tcoh) ) {
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
            if ( usedcandidate->data[jj] == 0 && (in[jj]->fsig-in[locs->data[0]]->fsig-0.25/params->Tcoh <= 0.5*iter/params->Tcoh && in[jj]->fsig-in[locs->data[0]]->fsig+0.25/params->Tcoh >= 0.5*iter/params->Tcoh) ) {
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
            if ( usedcandidate->data[jj] == 0 && (in[locs->data[0]]->fsig-in[jj]->fsig-0.25/params->Tcoh <= 0.5*iter/params->Tcoh && in[locs->data[0]]->fsig-in[jj]->fsig+0.25/params->Tcoh >= 0.5*iter/params->Tcoh) ) {
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
         if ( usedcandidate->data[locs->data[jj]] == 0 && fabs(params->Tobs/in[locs->data[jj]]->period - params->Tobs/in[locs->data[subsetloc]]->period) <= 1.0 ) {
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
            REAL4 weight = 0;
            REAL4 bestmoddepth = 0;
            REAL4 bestR = 0;
            REAL4 bestSNR = 0;
            for (kk=0; kk<loc2; kk++) {
               avefsig += in[locs2->data[kk]]->fsig*in[locs2->data[kk]]->snr;
               aveperiod += in[locs2->data[kk]]->period*in[locs2->data[kk]]->snr;
               weight += in[locs2->data[kk]]->snr;
               if (mindf > in[locs2->data[kk]]->moddepth || mindf == 0.0) mindf = in[locs2->data[kk]]->moddepth;
               if (maxdf < in[locs2->data[kk]]->moddepth) maxdf = in[locs2->data[kk]]->moddepth;
               
               if (loc2==1 && aveperiod/weight >= params->Pmin && aveperiod/weight <= params->Pmax) {
                  bestSNR = in[locs2->data[kk]]->snr;
                  bestmoddepth = in[locs2->data[kk]]->moddepth;
                  bestR = in[locs2->data[kk]]->stat;
               }
               
               usedcandidate->data[locs2->data[kk]] = 1;
            }
            avefsig = avefsig/weight;
            aveperiod = aveperiod/weight;
            
            if (loc2 > 1 && aveperiod >= params->Pmin && aveperiod <= params->Pmax) {
               INT4 numofmoddepths = (INT4)floorf(2*(maxdf-mindf)*params->Tcoh)+1;
               for (kk=0; kk<numofmoddepths; kk++) {
                  
                  candidate *cand = new_candidate();
                  loadCandidateData(cand, avefsig, aveperiod, mindf + kk*0.5/params->Tcoh, in[0]->ra, in[0]->dec, 0, 0, 0.0);
                  templateStruct *template = new_templateStruct(params->templatelength);
                  if (option==1) makeTemplate(template, cand, params, plan);
                  else makeTemplateGaussians(template, cand, params);
                  farStruct *farval = new_farStruct();
                  estimateFAR(farval, template, 10000, 0.01, ffplanenoise);
                  REAL4 R = calculateR(ffdata->ffdata, template, ffplanenoise);
                  REAL4 snr = (R - farval->distMean)/farval->distSigma;
                  if (snr > bestSNR) {
                     bestSNR = snr;
                     bestmoddepth = mindf + kk*0.5/params->Tcoh;
                     bestR = R;
                  }
                  free_candidate(cand);
                  cand = NULL;
                  free_templateStruct(template);
                  template = NULL;
                  free_farStruct(farval);
                  farval = NULL;
               }
            }
            
            if (bestSNR != 0.0) {
               out[numcandoutlist] = new_candidate();
               loadCandidateData(out[numcandoutlist], avefsig, aveperiod, bestmoddepth, in[0]->ra, in[0]->dec, bestR, bestSNR, 0.0);
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







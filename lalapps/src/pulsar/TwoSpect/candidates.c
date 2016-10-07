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

#include <sys/stat.h>
#include <lal/UserInput.h>
#include "TwoSpect.h"
#include "candidates.h"
#include "falsealarm.h"
#include "templates.h"

/**
 * Allocate a candidateVector
 * \param [in] length Length of the candidateVector
 * \return Pointer to the allocated candidateVector
 */
candidateVector * createcandidateVector(const UINT4 length)
{

   candidateVector *vector = NULL;
   XLAL_CHECK_NULL( (vector = XLALMalloc(sizeof(*vector))) != NULL, XLAL_ENOMEM );

   vector->length = length;
   vector->numofcandidates = 0;
   if (length==0) vector->data = NULL;
   else XLAL_CHECK_NULL( (vector->data = XLALMalloc( length*sizeof(*vector->data) )) != NULL, XLAL_ENOMEM );

   for (UINT4 ii=0; ii<vector->length; ii++) vector->data[ii].prob = 0.0;

   return vector;

} /* createcandidateVector() */


/**
 * Resize a candidateVector
 * \param [in] vector Pointer of vector to resize
 * \param [in] length New length of candidateVector
 * \return Pointer to resized vector
 */
candidateVector * resizecandidateVector(candidateVector *vector, const UINT4 length)
{

   if (vector==NULL) return createcandidateVector(length);
   if (length==0) {
      destroycandidateVector(vector);
      return NULL;
   }

   XLAL_CHECK_NULL( (vector->data = XLALRealloc(vector->data, length*sizeof(*vector->data))) != NULL, XLAL_ENOMEM );

   if (length > vector->length) for (UINT4 ii=vector->length; ii<length; ii++) vector->data[ii].prob = 0.0;

   vector->length = length;

   return vector;

} /* resizecandidateVector() */


/**
 * Free a candidateVector
 * \param [in] vector Pointer of candidateVector to be freed
 */
void destroycandidateVector(candidateVector *vector)
{
   if (vector==NULL) return;
   if ((!vector->length || !vector->data) && (vector->length || vector->data)) XLAL_ERROR_VOID(XLAL_EINVAL);
   if (vector->data) XLALFree((candidate*)vector->data);
   vector->data = NULL;
   XLALFree((candidateVector*)vector);
   return;
} /* destroycandidateVector() */


/**
 * Load candidate data
 * \param [out] output              Pointer to candidate
 * \param [in]  fsig                Frequency of candidate
 * \param [in]  period              Orbital period of candidate
 * \param [in]  moddepth            Modulation depth of candidate
 * \param [in]  ra                  Right ascension of candidate
 * \param [in]  dec                 Declination of candidate
 * \param [in]  statval             Detection statistic
 * \param [in]  h0                  Estimated strain amplitude
 * \param [in]  prob                False alarm probability
 * \param [in]  proberrcode         Davies' method error code
 * \param [in]  normalization       Time-frequency normalization
 * \param [in]  templateVectorIndex Index value of the template in a templateVector (can be -1 if not from a vector)
 * \param [in]  lineContamination   Boolean flag to indicate 0 = no contamination from lines or 1 = likely contaminated by one or more lines
 */
void loadCandidateData(candidate *output, const REAL8 fsig, const REAL8 period, const REAL8 moddepth, const REAL4 ra, const REAL4 dec, const REAL8 statval, const REAL8 h0, const REAL8 prob, const INT4 proberrcode, const REAL8 normalization, const INT4 templateVectorIndex, const BOOLEAN lineContamination)
{
   XLAL_CHECK_VOID( output != NULL, XLAL_EINVAL );
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
   output->templateVectorIndex = templateVectorIndex;
   output->lineContamination = lineContamination;
} // loadCandidateData()


/**
 * Analyze a single template
 * \param [out] output                 Pointer to candidate structure
 * \param [in]  input                  Pointer to candidate structure
 * \param [in]  ffdata                 Pointer to ffdataStruct
 * \param [in]  aveNoise               Pointer to REAL4VectorAligned of expected 2nd FFT background powers
 * \param [in]  aveTFnoisePerFbinRatio Pointer to REAL4VectorAligned of normalized power across the frequency band
 * \param [in]  params                 Pointer to UserInput_t
 * \param [in]  plan                   Pointer to REAL4FFTPlan
 * \param [in]  rng                    Pointer to gsl_rng
 * \param [in]  exactflag              Boolean value to indicate using exact templates
 * \return Status value
 */
INT4 analyzeOneTemplate(candidate *output, const candidate *input, const ffdataStruct *ffdata, const REAL4VectorAligned *aveNoise, const REAL4VectorAligned *aveTFnoisePerFbinRatio, const UserInput_t *params, const REAL4FFTPlan *plan, const gsl_rng *rng, const BOOLEAN exactflag)
{

   XLAL_CHECK( output!=NULL && input!=NULL && ffdata!=NULL && aveNoise!=NULL && aveTFnoisePerFbinRatio!=NULL && params!=NULL && plan!=NULL && rng!=NULL, XLAL_EINVAL );

   INT4 proberrcode = 0;

   //Allocate and make the template
   TwoSpectTemplate *template = NULL;
   XLAL_CHECK( (template = createTwoSpectTemplate(params->maxTemplateLength)) != NULL, XLAL_EFUNC );
   resetTwoSpectTemplate(template);
   if (exactflag) XLAL_CHECK( makeTemplate(template, *input, params, plan) == XLAL_SUCCESS, XLAL_EFUNC );
   else XLAL_CHECK( makeTemplateGaussians(template, *input, params) == XLAL_SUCCESS, XLAL_EFUNC );

   //Calculate R from the template and the data
   REAL8 R = calculateR(ffdata->ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
   XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );

   //Calculate FAP
   REAL8 prob = probR(template, aveNoise, aveTFnoisePerFbinRatio, R, params, rng, &proberrcode);
   XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );

   //Estimate the h0 if R>0.0
   REAL8 h0 = 0.0;
   if ( R > 0.0 ) h0 = 2.7426*pow(R/(params->Tsft*params->Tobs),0.25);

   loadCandidateData(output, input->fsig, input->period, input->moddepth, input->ra, input->dec, R, h0, prob, proberrcode, 1.0, -1, 0);

   destroyTwoSpectTemplate(template);

   return XLAL_SUCCESS;
}

INT4 analyzeCandidatesTemplateFromVector(candidateVector *output, const candidateVector *input, const TwoSpectTemplateVector *vector, const ffdataStruct *ffdata, const REAL4VectorAligned *aveNoise, const REAL4VectorAligned *aveTFnoisePerFbinRatio, const UserInput_t *params, const gsl_rng *rng, const UINT4 templateLen)
{

   XLAL_CHECK( output!=NULL && input!=NULL && vector!=NULL && ffdata!=NULL && aveNoise!=NULL && aveTFnoisePerFbinRatio!=NULL && params!=NULL && rng!=NULL, XLAL_EINVAL );

   TwoSpectTemplate *template = NULL;
   XLAL_CHECK( (template = createTwoSpectTemplate(templateLen)) != NULL, XLAL_EFUNC );

   for (UINT4 ii=0; ii<input->length; ii++) {
      INT4 proberrcode = 0;

      //First convert the template to the right pixels
      resetTwoSpectTemplate(template);
      XLAL_CHECK( convertTemplateForSpecificFbin(template, vector->data[input->data[ii].templateVectorIndex], input->data[ii].fsig, params) == XLAL_SUCCESS, XLAL_EFUNC );

      //Calculate R from the template and the data
      REAL8 R = calculateR(ffdata->ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
      XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );

      //Calculate FAP
      REAL8 prob = probR(template, aveNoise, aveTFnoisePerFbinRatio, R, params, rng, &proberrcode);
      XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );

      //Estimate the h0 if R>0.0
      REAL8 h0 = 0.0;
      if ( R > 0.0 ) h0 = 2.7426*pow(R/(params->Tsft*params->Tobs),0.25);

      if (prob < output->data[output->length-1].prob) {
         UINT4 insertionPoint = output->length - 1;
         while(insertionPoint>0 && prob<output->data[insertionPoint - 1].prob) insertionPoint--;
         for (INT4 kk=(INT4)output->length-2; kk>=(INT4)insertionPoint; kk--) loadCandidateData(&(output->data[kk+1]), output->data[kk].fsig, output->data[kk].period, output->data[kk].moddepth, output->data[kk].ra, output->data[kk].dec, output->data[kk].stat, output->data[kk].h0, output->data[kk].prob, output->data[kk].proberrcode, output->data[kk].normalization, output->data[kk].templateVectorIndex, 0);
         loadCandidateData(&(output->data[insertionPoint]), template->f0, template->period, template->moddepth, input->data[ii].ra, input->data[ii].dec, R, h0, prob, proberrcode, ffdata->tfnormalization, input->data[ii].templateVectorIndex, 0);
         if (output->numofcandidates<output->length) output->numofcandidates++;
      }
   }

   destroyTwoSpectTemplate(template);

   return XLAL_SUCCESS;
}

/**
 * A brute force template search to find the most significant template around a candidate
 * \param [out] output                 Pointer to candidate structure
 * \param [in]  input                  Input candidate structure
 * \param [in]  paramspace             Pointer to TwoSpectParamSpaceSearchVals containing the parameter space to be searched
 * \param [in]  params                 Pointer to UserInput_t
 * \param [in]  ffdata                 Pointer to ffdataStruct
 * \param [in]  aveNoise               Pointer to REAL4VectorAligned of 2nd FFT background powers
 * \param [in]  aveTFnoisePerFbinRatio Pointer to REAL4VectorAligned of normalized power across the frequency band
 * \param [in]  secondFFTplan          Pointer to REAL4FFTPlan
 * \param [in]  rng                    Pointer to gsl_rng
 * \param [in]  useExactTemplates      Boolean of 0 (use Gaussian templates) or 1 (use exact templates)
 * \return Status value
 */
INT4 bruteForceTemplateSearch(candidate *output, const candidate input, const TwoSpectParamSpaceSearchVals *paramspace, const UserInput_t *params, const REAL4VectorAligned *ffdata, const REAL4VectorAligned *aveNoise, const REAL4VectorAligned *aveTFnoisePerFbinRatio, const REAL4FFTPlan *secondFFTplan, const gsl_rng *rng, const BOOLEAN useExactTemplates)
{

   XLAL_CHECK( output != NULL && params != NULL && aveNoise != NULL && aveTFnoisePerFbinRatio != NULL && secondFFTplan != NULL, XLAL_EINVAL );

   fprintf(stderr, "Performing brute force template search... ");

   REAL8Vector *trialf, *trialb, *trialp;
   REAL8 fstepsize, dfstepsize;
   //REAL4 tcohfactor = 1.49e-3*params->Tsft + 1.76;    //From in-text equation after Eq. 23 of E.G. and K.R. 2011
   REAL4 alpha0 = 45.0*(params->Tsft/1800.0)+30.0;
   REAL8 log10templatefar = log10(params->tmplfar);

   TwoSpectParamSpaceSearchVals search = {paramspace->fminimum, paramspace->fmaximum, paramspace->numfsteps, paramspace->numperiodslonger, paramspace->numperiodsshorter, paramspace->periodSpacingFactor, paramspace->dfmin, paramspace->dfmax, paramspace->numdfsteps};

   //Set up parameters of modulation depth search
   if (search.dfmin<params->dfmin) search.dfmin = params->dfmin;
   if (search.dfmax>params->dfmax) search.dfmax = params->dfmax;
   XLAL_CHECK( (trialb = XLALCreateREAL8Vector(search.numdfsteps)) != NULL, XLAL_EFUNC );
   if (search.numdfsteps>1) {
      dfstepsize = (search.dfmax-search.dfmin)/(REAL8)(search.numdfsteps-1);
      for (UINT4 ii=0; ii<search.numdfsteps; ii++) trialb->data[ii] = search.dfmin + dfstepsize*ii;
   } else {
      trialb->data[0] = 0.5*(search.dfmin+search.dfmax);
   }

   //Set up parameters of signal frequency search
   if (search.fminimum<params->fmin) search.fminimum = params->fmin;
   if (search.fmaximum>params->fmin+params->fspan) search.fmaximum = params->fmin+params->fspan;
   XLAL_CHECK( (trialf = XLALCreateREAL8Vector(search.numfsteps)) != NULL, XLAL_EFUNC );
   if (search.numfsteps>1) {
      fstepsize = (search.fmaximum-search.fminimum)/(REAL8)(search.numfsteps-1);
      for (UINT4 ii=0; ii<search.numfsteps; ii++) trialf->data[ii] = search.fminimum + fstepsize*ii;
   } else {
      trialf->data[0] = 0.5*(search.fminimum+search.fmaximum);
   }

   //Search over numperiods different periods
   XLAL_CHECK( (trialp = XLALCreateREAL8Vector(search.numperiodslonger+search.numperiodsshorter+1)) != NULL, XLAL_EFUNC );

   //Now search over the parameter space. Frequency, then modulation depth, then period
   //Initialze best values as the initial point we are searching around
   INT4 bestproberrcode = 0;
   REAL8 bestf = 0.0, bestp = 0.0, bestdf = 0.0, bestR = 0.0, besth0 = 0.0, bestProb = 0.0;
   candidate cand;
   TwoSpectTemplate *template = NULL;
   XLAL_CHECK( (template = createTwoSpectTemplate(params->maxTemplateLength)) != NULL, XLAL_EFUNC );
   farStruct *farval = NULL;
   if (params->calcRthreshold) {
      XLAL_CHECK( (farval = createfarStruct()) != NULL, XLAL_EFUNC );
   }

   INT4 startposition = search.numperiodsshorter, proberrcode = 0;
   //Search over frequency
   for (UINT4 ii=0; ii<trialf->length; ii++) {
      //Search over modulation depth
      for (UINT4 jj=0; jj<trialb->length; jj++) {
         //Start with period of the first guess, then determine nearest neighbor from the
         //modulation depth amplitude to find the other period guesses. These parameters
         //are determined from simulation to scale the N.N. distance w.r.t. mod. depth with
         //20% mismatch parameter
         trialp->data[startposition] = input.period;
         for (UINT4 kk=0; kk<search.numperiodsshorter; kk++) {
            REAL8 nnp = search.periodSpacingFactor*trialp->data[startposition-kk]*trialp->data[startposition-kk]*(1+trialp->data[startposition-kk]/(alpha0*sqrt(trialb->data[jj])*params->Tobs))/(alpha0*params->Tobs*sqrt(trialb->data[jj]));
            trialp->data[startposition-(kk+1)] = trialp->data[startposition-kk] - nnp;
         }
         for (UINT4 kk=0; kk<search.numperiodslonger; kk++) {
            REAL8 nnp = search.periodSpacingFactor*trialp->data[startposition+kk]*trialp->data[startposition+kk]*(1+trialp->data[startposition+kk]/(alpha0*sqrt(trialb->data[jj])*params->Tobs))/(alpha0*params->Tobs*sqrt(trialb->data[jj]));
            trialp->data[startposition+(kk+1)] = trialp->data[startposition+kk] + nnp;
         }

         //Search over period
         for (UINT4 kk=0; kk<trialp->length; kk++) {
            //Within boundaries?
            if ( trialf->data[ii]>=params->fmin &&
                trialf->data[ii]<(params->fmin+params->fspan) &&
                trialb->data[jj]<maxModDepth(trialp->data[kk], params->Tsft) &&
                trialp->data[kk]>minPeriod(trialb->data[jj], params->Tsft) &&
                trialp->data[kk]<=(0.2*params->Tobs) &&
                trialp->data[kk]>=(4.0*params->Tsft) &&
                trialb->data[jj]>=params->dfmin &&
                trialb->data[jj]<=params->dfmax &&
                trialp->data[kk]<=params->Pmax &&
                trialp->data[kk]>=params->Pmin ) {

               loadCandidateData(&cand, trialf->data[ii], trialp->data[kk], trialb->data[jj], input.ra, input.dec, 0, 0, 0.0, 0, 0.0, -1, 0);

               resetTwoSpectTemplate(template);

               if (useExactTemplates) XLAL_CHECK( makeTemplate(template, cand, params, secondFFTplan) == XLAL_SUCCESS, XLAL_EFUNC );
               else XLAL_CHECK( makeTemplateGaussians(template, cand, params) == XLAL_SUCCESS, XLAL_EFUNC );

               if (params->calcRthreshold && bestProb==0.0) XLAL_CHECK( numericFAR(farval, template, params->tmplfar, aveNoise, aveTFnoisePerFbinRatio, params, rng, params->BrentsMethod) == XLAL_SUCCESS, XLAL_EFUNC );

               REAL8 R = calculateR(ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
               XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
               REAL8 prob = probR(template, aveNoise, aveTFnoisePerFbinRatio, R, params, rng, &proberrcode);
               XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
               REAL8 h0 = 0.0;
               if ( R > 0.0 ) h0 = 2.7426*pow(R/(params->Tsft*params->Tobs),0.25);

               //fprintf(stderr, "%.8g %.9g %g %.14g\n", trialf->data[ii], trialp->data[kk], trialb->data[jj], R);

               if ( (bestProb!=0.0 && prob < bestProb) || (bestProb==0.0 && !params->calcRthreshold && prob<log10templatefar) || (bestProb==0.0 && params->calcRthreshold && R > farval->far) ) {
                  bestf = trialf->data[ii];
                  bestp = trialp->data[kk];
                  bestdf = trialb->data[jj];
                  bestR = R;
                  besth0 = h0;
                  bestProb = prob;
                  bestproberrcode = proberrcode;
               }

            } /* if within boundaries */
         } /* for kk < trialp */
      } /* for jj < trialb */
   } /* for ii < trialf */
   destroyTwoSpectTemplate(template);
   template = NULL;
   if (params->calcRthreshold) {
      destroyfarStruct(farval);
      farval = NULL;
   }
   XLALDestroyREAL8Vector(trialf);
   XLALDestroyREAL8Vector(trialb);
   XLALDestroyREAL8Vector(trialp);
   trialf = NULL;
   trialb = NULL;
   trialp = NULL;

   if (bestProb==0.0) loadCandidateData(output, input.fsig, input.period, input.moddepth, input.ra, input.dec, input.stat, input.h0, input.prob, input.proberrcode, input.normalization, input.templateVectorIndex, input.lineContamination);
   else loadCandidateData(output, bestf, bestp, bestdf, input.ra, input.dec, bestR, besth0, bestProb, bestproberrcode, input.normalization, input.templateVectorIndex, input.lineContamination);

   fprintf(stderr, "done\n");

   return XLAL_SUCCESS;

}

/**
 * A brute force template search to test templates around a candidate
 * \param [out] output                 Pointer to a pointer of a candidateVector
 * \param [in]  input                  Input candidate structure
 * \param [in]  paramspace             Pointer to TwoSpectParamSpaceSearchVals containing the parameter space to be searched
 * \param [in]  params                 Pointer to UserInput_t
 * \param [in]  ffdata                 Pointer to ffdataStruct
 * \param [in]  aveNoise               Pointer to REAL4VectorAligned of 2nd FFT background powers
 * \param [in]  aveTFnoisePerFbinRatio Pointer to REAL4VectorAligned of normalized power across the frequency band
 * \param [in]  secondFFTplan          Pointer to REAL4FFTPlan
 * \param [in]  rng                    Pointer to gsl_rng
 * \param [in]  useExactTemplates      Boolean of 0 (use Gaussian templates) or 1 (use exact templates)
 * \return Status value
 */
INT4 bruteForceTemplateTest(candidateVector **output, const candidate input, const TwoSpectParamSpaceSearchVals *paramspace, const UserInput_t *params, const REAL4VectorAligned *ffdata, const REAL4VectorAligned *aveNoise, const REAL4VectorAligned *aveTFnoisePerFbinRatio, const REAL4FFTPlan *secondFFTplan, const gsl_rng *rng, const BOOLEAN useExactTemplates)
{

   XLAL_CHECK( *output != NULL && paramspace!=NULL && params != NULL && ffdata!=NULL && aveNoise != NULL && aveTFnoisePerFbinRatio != NULL && secondFFTplan != NULL && rng!=NULL, XLAL_EINVAL );

   REAL8Vector *trialf, *trialb, *trialp;
   REAL8 fstepsize, dfstepsize;
   REAL4 tcohfactor = 1.49e-3*params->Tsft + 1.76;    //From in-text equation after Eq. 23 of E.G. and K.R. 2011

   TwoSpectParamSpaceSearchVals search = {paramspace->fminimum, paramspace->fmaximum, paramspace->numfsteps, paramspace->numperiodslonger, paramspace->numperiodsshorter, paramspace->periodSpacingFactor, paramspace->dfmin, paramspace->dfmax, paramspace->numdfsteps};

   //Set up parameters of modulation depth search
   if (search.dfmin<params->dfmin) search.dfmin = params->dfmin;
   if (search.dfmax>params->dfmax) search.dfmax = params->dfmax;
   XLAL_CHECK( (trialb = XLALCreateREAL8Vector(search.numdfsteps)) != NULL, XLAL_EFUNC );
   if (search.numdfsteps>1) {
      dfstepsize = (search.dfmax-search.dfmin)/(REAL8)(search.numdfsteps-1);
      for (UINT4 ii=0; ii<search.numdfsteps; ii++) trialb->data[ii] = search.dfmin + dfstepsize*ii;
   } else {
      trialb->data[0] = 0.5*(search.dfmin+search.dfmax);
   }

   //Set up parameters of signal frequency search
   if (search.fminimum<params->fmin) search.fminimum = params->fmin;
   if (search.fmaximum>params->fmin+params->fspan) search.fmaximum = params->fmin+params->fspan;
   XLAL_CHECK( (trialf = XLALCreateREAL8Vector(search.numfsteps)) != NULL, XLAL_EFUNC );
   if (search.numfsteps>1) {
      fstepsize = (search.fmaximum-search.fminimum)/(REAL8)(search.numfsteps-1);
      for (UINT4 ii=0; ii<search.numfsteps; ii++) trialf->data[ii] = search.fminimum + fstepsize*ii;
   } else {
      trialf->data[0] = 0.5*(search.fminimum+search.fmaximum);
   }

   //Search over numperiods different periods
   XLAL_CHECK( (trialp = XLALCreateREAL8Vector(search.numperiodslonger+search.numperiodsshorter+1)) != NULL, XLAL_EFUNC );

   //Now search over the parameter space. Frequency, then modulation depth, then period
   candidate cand;
   TwoSpectTemplate *template = NULL;
   XLAL_CHECK( (template = createTwoSpectTemplate(params->maxTemplateLength)) != NULL, XLAL_EFUNC );

   INT4 startposition = search.numperiodsshorter, proberrcode = 0;
   //Search over frequency
   for (UINT4 ii=0; ii<trialf->length; ii++) {
      //Search over modulation depth
      for (UINT4 jj=0; jj<trialb->length; jj++) {
         //Start with period of the first guess, then determine nearest neighbor from the
         //modulation depth amplitude to find the other period guesses. These parameters
         //are determined from simulation to scale the N.N. distance w.r.t. mod. depth with
         //20% mismatch parameter
         trialp->data[startposition] = input.period;
         for (UINT4 kk=0; kk<search.numperiodsshorter; kk++) {
            REAL8 nnp = search.periodSpacingFactor*trialp->data[startposition-kk]*trialp->data[startposition-kk]*(1+trialp->data[startposition-kk]/tcohfactor/params->Tobs)/tcohfactor/params->Tobs*sqrt(3.6e-3/trialb->data[jj]);
            trialp->data[startposition-(kk+1)] = trialp->data[startposition-kk] - nnp;
         }
         for (UINT4 kk=0; kk<search.numperiodslonger; kk++) {
            REAL8 nnp = search.periodSpacingFactor*trialp->data[startposition+kk]*trialp->data[startposition+kk]*(1+trialp->data[startposition+kk]/tcohfactor/params->Tobs)/tcohfactor/params->Tobs*sqrt(3.6e-3/trialb->data[jj]);
            trialp->data[startposition+(kk+1)] = trialp->data[startposition+kk] + nnp;
         }

         //Search over period
         for (UINT4 kk=0; kk<trialp->length; kk++) {
            //Within boundaries?
            if ( trialf->data[ii]>=params->fmin &&
                trialf->data[ii]<(params->fmin+params->fspan) &&
                trialb->data[jj]<maxModDepth(trialp->data[kk], params->Tsft) &&
                trialp->data[kk]>minPeriod(trialb->data[jj], params->Tsft) &&
                trialp->data[kk]<=(0.2*params->Tobs) &&
                trialp->data[kk]>=(4.0*params->Tsft) &&
                trialb->data[jj]>=params->dfmin &&
                trialb->data[jj]<=params->dfmax &&
                trialp->data[kk]<=params->Pmax &&
                trialp->data[kk]>=params->Pmin ) {

               loadCandidateData(&cand, trialf->data[ii], trialp->data[kk], trialb->data[jj], input.ra, input.dec, 0, 0, 0.0, 0, 0.0, -1, 0);

               resetTwoSpectTemplate(template);

               if (useExactTemplates) XLAL_CHECK( makeTemplate(template, cand, params, secondFFTplan) == XLAL_SUCCESS, XLAL_EFUNC );
               else XLAL_CHECK( makeTemplateGaussians(template, cand, params) == XLAL_SUCCESS, XLAL_EFUNC );

               REAL8 R = calculateR(ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
               XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
               REAL8 prob = probR(template, aveNoise, aveTFnoisePerFbinRatio, R, params, rng, &proberrcode);
               XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
               REAL8 h0 = 0.0;
               if ( R > 0.0 ) h0 = 2.7426*pow(R/(params->Tsft*params->Tobs),0.25);

               //Resize the output candidate vector if necessary
               if ((*output)->numofcandidates == (*output)->length-1) XLAL_CHECK( (*output = resizecandidateVector(*output, 2*(*output)->length)) != NULL, XLAL_EFUNC );

               loadCandidateData(&((*output)->data[(*output)->numofcandidates]), trialf->data[ii], trialp->data[kk], trialb->data[jj], input.ra, input.dec, R, h0, prob, proberrcode, input.normalization, input.templateVectorIndex, input.lineContamination);
               (*output)->numofcandidates++;

            } /* if within boundaries */
         } /* for kk < trialp */
      } /* for jj < trialb */
   } /* for ii < trialf */
   destroyTwoSpectTemplate(template);
   XLALDestroyREAL8Vector(trialf);
   XLALDestroyREAL8Vector(trialb);
   XLALDestroyREAL8Vector(trialp);

   return XLAL_SUCCESS;

}


/**
 * A brute force template search to find the most significant template around a putative source whose parameters are somewhat constrained
 * \param [out] output                 Pointer to a pointer of a candidateVector
 * \param [in]  fminimum               Lower frequency bound to search (inclusive)
 * \param [in]  fspan                  Span of the frequency band (inclusive of endpoint)
 * \param [in]  period                 Specific orbital period (measured in seconds)
 * \param [in]  asini                  Specific projected semi-major axis (measured in light seconds)
 * \param [in]  asinisigma             Uncertainty on the specific asini value (measured in light seconds)
 * \param [in]  skypos                 SkyPosition struct of the sky position (in RA and DEC) being searched
 * \param [in]  params                 Pointer to UserInput_t
 * \param [in]  ffdata                 Pointer to ffdataStruct
 * \param [in]  aveNoise               Pointer to REAL4VectorAligned of 2nd FFT background powers
 * \param [in]  aveTFnoisePerFbinRatio Pointer to REAL4VectorAligned of normalized power across the frequency band
 * \param [in]  trackedlines           Pointer to REAL4VectorSequence of lines (allowed to be NULL if no lines)
 * \param [in]  secondFFTplan          Pointer to REAL4FFTPlan
 * \param [in]  rng                    Pointer to gsl_rng
 * \param [in]  useExactTemplates      Boolean of 0 (use Gaussian templates) or 1 (use exact templates)
 * \return Status value
 */
INT4 templateSearch_scox1Style(candidateVector **output, const REAL8 fminimum, const REAL8 fspan, const REAL8 period, const REAL8 asini, const REAL8 asinisigma, const SkyPosition skypos, const UserInput_t *params, const REAL4VectorAligned *ffdata, const REAL4VectorAligned *aveNoise, const REAL4VectorAligned *aveTFnoisePerFbinRatio, const REAL4VectorSequence *trackedlines, const REAL4FFTPlan *secondFFTplan, const gsl_rng *rng, BOOLEAN useExactTemplates)
{

   XLAL_CHECK( *output != NULL && params != NULL && ffdata != NULL && aveNoise != NULL && aveTFnoisePerFbinRatio != NULL && secondFFTplan != NULL && rng != NULL, XLAL_EINVAL );

   REAL8Vector *trialf;
   REAL8Vector *trialdf;
   REAL8 fstepsize;
   REAL8 dfstepsize;
   
   //Set up parameters of signal frequency search
   INT4 numfsteps = (INT4)round(2.0*fspan*params->Tsft)+1;
   XLAL_CHECK( (trialf = XLALCreateREAL8Vector(numfsteps)) != NULL, XLAL_EFUNC );
   fstepsize = fspan/(REAL8)(numfsteps-1);
   for (INT4 ii=0; ii<numfsteps; ii++) trialf->data[ii] = fminimum + fstepsize*ii;

   //Now search over the frequencies
   INT4 proberrcode = 0;
   candidate cand;
   TwoSpectTemplate *template = NULL;
   XLAL_CHECK( (template = createTwoSpectTemplate(params->maxTemplateLength)) != NULL, XLAL_EFUNC );

   //Search over frequency
   for (UINT4 ii=0; ii<trialf->length; ii++) {

      //Set up parameters of signal modulation depth search
      /* Modulation depth is 2*pi*f*asini*period, or rearranged
         0.8727*(f/1000.0)*(7200.0/period)*asini
         Assuming sigma = 0.18 uncertainty in an asini of 1.44 for
         Scorpius X-1 and giving +/- 3*sigma leeway, the conservative 
         number of df steps should cover
         0.8727*(fmax/1000.0)*(7200.0/period)*6*0.18 
         with, as empirical testing has found necessary, a spacing of
         4*Tsft */
      /* While this initialization could be moved inside the loop
         that searches over frequency, it is slightly faster not to have to 
         recalculate these variables every time,
         and it gives us a bit of extra data*/
      //REAL8 asinisigma = 0.18; typical for Scorpius X-1 with 2014 data
      REAL8 moddepth = 0.8727*(trialf->data[ii]/1000.0)*(7200.0/period)*asini;
      //fprintf(stderr,"Making the first computation involving asinisigma, for moddepthmin\n");
      //Note, 6*asinsigma for a span of plus/minus 3*asinisigma
      REAL8 moddepthspan = 0.8727*(trialf->data[numfsteps-1]/1000.0)*(7200.0/period)*6*asinisigma;
      //fprintf(stderr,"intended moddepthspan: %f \n", moddepthspan);
      //fprintf(stderr,"Done with moddepthspan, making moddepthmin\n");
      REAL8 moddepthmin = moddepth - 0.5*moddepthspan;
      //fprintf(stderr,"intended moddepthmin: %f \n", moddepthmin);
      INT4 numdfsteps = (INT4)round(4.0*moddepthspan*params->Tsft) + 1;
      //fprintf(stderr,"intended numdfsteps: %d \n", numdfsteps);
      trialdf = XLALCreateREAL8Vector(numdfsteps);
      XLAL_CHECK( trialdf != NULL, XLAL_EFUNC);
      if (numdfsteps > 1) {
         dfstepsize = moddepthspan/(REAL8)(numdfsteps-1);
      } else {
         dfstepsize = 0;
      }
      for (INT4 jj=0; jj<numdfsteps; jj++) trialdf->data[jj] = moddepthmin + dfstepsize*jj;

      for (UINT4 jj=0; jj<trialdf->length; jj++) {
         //Determine modulation depth
         //REAL8 moddepth = 0.8727*(trialf->data[ii]/1000.0)*(7200.0/period)*asini;

         //load candidate
         //printf(stderr,"Loading candidate. Remember to get the RA and dec from outside in production run\n");
         loadCandidateData(&cand, trialf->data[ii], period, trialdf->data[jj], skypos.longitude, skypos.latitude, 0, 0, 0.0, 0, 0.0, -1, 0);

         //Make the template
         resetTwoSpectTemplate(template);
         if (useExactTemplates!=0) XLAL_CHECK( makeTemplate(template, cand, params, secondFFTplan) == XLAL_SUCCESS, XLAL_EFUNC );
         else XLAL_CHECK( makeTemplateGaussians(template, cand, params) == XLAL_SUCCESS, XLAL_EFUNC );

         REAL8 R = calculateR(ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
         XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC);
         REAL8 prob = probR(template, aveNoise, aveTFnoisePerFbinRatio, R, params, rng, &proberrcode);
         XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC);
         REAL8 h0 = 0.0;
         if ( R > 0.0 ) h0 = 2.7426*pow(R/(params->Tsft*params->Tobs),0.25);

         //Line contamination?
         BOOLEAN lineContamination = 0;
         if (trackedlines!=NULL) {
            UINT4 kk = 0;
            while (kk<trackedlines->length && lineContamination==0) {
               if (2.0*trialdf->data[jj]>=(trackedlines->data[kk*3+2]-trackedlines->data[kk*3+1])) {
                  if ((trackedlines->data[kk*3+2]>=(REAL4)(trialf->data[ii]-trialdf->data[jj]) && trackedlines->data[kk*3+2]<=(REAL4)(trialf->data[ii]+trialdf->data[jj])) ||
                      (trackedlines->data[kk*3+1]>=(REAL4)(trialf->data[ii]-trialdf->data[jj]) && trackedlines->data[kk*3+1]<=(REAL4)(trialf->data[ii]+trialdf->data[jj]))) {
                     lineContamination = 1;
                  }
               } // if the band spanned by the line is smaller than the band spanned by the signal
               else {
                  if (((REAL4)(trialf->data[ii]+trialdf->data[jj])>=trackedlines->data[kk*3+1] && (REAL4)(trialf->data[ii]+trialdf->data[jj])<=trackedlines->data[kk*3+2]) ||
                      ((REAL4)(trialf->data[ii]-trialdf->data[jj])>=trackedlines->data[kk*3+1] && (REAL4)(trialf->data[ii]-trialdf->data[jj])<=trackedlines->data[kk*3+2])) {
                     lineContamination = 1;
                  }
               } // instead if the band spanned by the line is larger than the band spanned by the signal
               kk++;
            } // while kk < trackedlines->length && lineContamination==0
         } // if trackedlines != NULL

         //Resize the output candidate vector if necessary
         if ((*output)->numofcandidates == (*output)->length-1) {
            *output = resizecandidateVector(*output, 2*((*output)->length));
            XLAL_CHECK( *output != NULL, XLAL_EFUNC);
         }

         loadCandidateData(&((*output)->data[(*output)->numofcandidates]), trialf->data[ii], period, trialdf->data[jj], skypos.longitude, skypos.latitude, R, h0, prob, proberrcode, 0.0, -1, lineContamination);
         (*output)->numofcandidates++;
      } /* for jj < trialdf */   
      XLALDestroyREAL8Vector(trialdf);
      trialdf = NULL;
   } /* for ii < trialf */
   destroyTwoSpectTemplate(template);
   template = NULL;
   XLALDestroyREAL8Vector(trialf);
   trialf = NULL;

   return XLAL_SUCCESS;

}

/**
 * A brute force template search to find the most significant template at a fixed modulation depth around a putative source whose parameters are somewhat constrained
 * \param [out] output                 Pointer to a pointer of a candidateVector
 * \param [in]  dffixed                Modulation depth (fixed by user)
 * \param [in]  fminimum               Lower frequency bound to search (inclusive)
 * \param [in]  fspan                  Span of the frequency band (inclusive of endpoint)
 * \param [in]  period                 Specific orbital period (measured in seconds)
 * \param [in]  skypos                 SkyPosition struct of the sky position (in RA and DEC) being searched
 * \param [in]  params                 Pointer to UserInput_t
 * \param [in]  ffdata                 Pointer to ffdataStruct
 * \param [in]  aveNoise               Pointer to REAL4VectorAligned of 2nd FFT background powers
 * \param [in]  aveTFnoisePerFbinRatio Pointer to REAL4VectorAligned of normalized power across the frequency band
* \param [in]  trackedlines           Pointer to REAL4VectorSequence of lines (allowed to be NULL if no lines)
 * \param [in]  secondFFTplan          Pointer to REAL4FFTPlan
 * \param [in]  rng                    Pointer to gsl_rng
 * \param [in]  useExactTemplates      Boolean of 0 (use Gaussian templates) or 1 (use exact templates)
 * \return Status value
 */


INT4 templateSearch_fixedDf(candidateVector **output, const LALStringVector *dffixed, const REAL8 fminimum, const REAL8 fspan, const REAL8 period, const SkyPosition skypos, const UserInput_t *params, const REAL4VectorAligned *ffdata, const REAL4VectorAligned *aveNoise, const REAL4VectorAligned *aveTFnoisePerFbinRatio, const REAL4VectorSequence *trackedlines, const REAL4FFTPlan *secondFFTplan, const gsl_rng *rng, BOOLEAN useExactTemplates)
{

   XLAL_CHECK( *output != NULL && dffixed !=NULL && params != NULL && ffdata != NULL && aveNoise != NULL && aveTFnoisePerFbinRatio != NULL && secondFFTplan != NULL && rng != NULL, XLAL_EINVAL );

   REAL8Vector *trialf;
   REAL8Vector *trialdf;
   REAL8 fstepsize;

   // Create the vector trialdf
   XLAL_CHECK( (trialdf = XLALCreateREAL8Vector((dffixed->length))) != NULL, XLAL_EFUNC );
   for (UINT4 ii=0;ii<dffixed->length;ii++) {
      XLAL_CHECK( XLALParseStringValueAsREAL8(&(trialdf->data[ii]), dffixed->data[ii])== XLAL_SUCCESS, XLAL_EFUNC );

      // Check that the specified df is ok; if not, return an error and end the program.
      XLAL_CHECK ( (trialdf->data[ii] < maxModDepth(period, params->Tsft)) && (trialdf->data[ii] > 0.5/params->Tsft), XLAL_EFAILED, "ERROR: Modulation depth must be between %.5f and %.5f.\n",0.5/params->Tsft,maxModDepth(period,params->Tsft) );
   }

   //Set up parameters of signal frequency search
   INT4 numfsteps = (INT4)round(2.0*fspan*params->Tsft)+1;
   XLAL_CHECK( (trialf = XLALCreateREAL8Vector(numfsteps)) != NULL, XLAL_EFUNC );
   fstepsize = fspan/(REAL8)(numfsteps-1);
   for (INT4 ii=0; ii<numfsteps; ii++) trialf->data[ii] = fminimum + fstepsize*ii;

   //Now search over the frequencies
   INT4 proberrcode = 0;
   candidate cand;
   TwoSpectTemplate *template = NULL;
   XLAL_CHECK( (template = createTwoSpectTemplate(params->maxTemplateLength)) != NULL, XLAL_EFUNC );

   // loop over dfs
   for (UINT4 jj=0; jj<trialdf->length; jj++){

   //Search over frequency
   for (UINT4 ii=0; ii<trialf->length; ii++) {

      loadCandidateData(&cand, trialf->data[ii], period, trialdf->data[jj], skypos.longitude, skypos.latitude, 0, 0, 0.0, 0, 0.0, -1, 0);

      //Make the template
      resetTwoSpectTemplate(template);
      if (useExactTemplates!=0) XLAL_CHECK( makeTemplate(template, cand, params, secondFFTplan) == XLAL_SUCCESS, XLAL_EFUNC );
      else XLAL_CHECK( makeTemplateGaussians(template, cand, params) == XLAL_SUCCESS, XLAL_EFUNC );

      REAL8 R = calculateR(ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
      XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC);
      REAL8 prob = probR(template, aveNoise, aveTFnoisePerFbinRatio, R, params, rng, &proberrcode);
      XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC);
      REAL8 h0 = 0.0;
      if ( R > 0.0 ) h0 = 2.7426*pow(R/(params->Tsft*params->Tobs),0.25);

         //Line contamination?
         BOOLEAN lineContamination = 0;
         if (trackedlines!=NULL) {
            UINT4 kk = 0;
            while (kk<trackedlines->length && lineContamination==0) {
               if (2.0*trialdf->data[jj]>=(trackedlines->data[kk*3+2]-trackedlines->data[kk*3+1])) {
                  if ((trackedlines->data[kk*3+2]>=(REAL4)(trialf->data[ii]-trialdf->data[jj]) && trackedlines->data[kk*3+2]<=(REAL4)(trialf->data[ii]+trialdf->data[jj])) ||
                      (trackedlines->data[kk*3+1]>=(REAL4)(trialf->data[ii]-trialdf->data[jj]) && trackedlines->data[kk*3+1]<=(REAL4)(trialf->data[ii]+trialdf->data[jj]))) {
                     lineContamination = 1;
                  }
               } // if the band spanned by the line is smaller than the band spanned by the signal
               else {
                  if (((REAL4)(trialf->data[ii]+trialdf->data[jj])>=trackedlines->data[kk*3+1] && (REAL4)(trialf->data[ii]+trialdf->data[jj])<=trackedlines->data[kk*3+2]) ||
                      ((REAL4)(trialf->data[ii]-trialdf->data[jj])>=trackedlines->data[kk*3+1] && (REAL4)(trialf->data[ii]-trialdf->data[jj])<=trackedlines->data[kk*3+2])) {
                     lineContamination = 1;
                  }
               } // instead if the band spanned by the line is larger than the band spanned by the signal
               kk++;
            } // while kk < trackedlines->length && lineContamination==0
         } // if trackedlines != NULL

         //Resize the output candidate vector if necessary
         if ((*output)->numofcandidates == (*output)->length-1) {
            *output = resizecandidateVector(*output, 2*((*output)->length));
            XLAL_CHECK( *output != NULL, XLAL_EFUNC);
         }

         loadCandidateData(&((*output)->data[(*output)->numofcandidates]), trialf->data[ii], period, trialdf->data[jj], skypos.longitude, skypos.latitude, R, h0, prob, proberrcode, 0.0, -1, lineContamination);
         (*output)->numofcandidates++;

         } /* for ii < trialf */
      } /* for jj < trialdf */

      XLALDestroyREAL8Vector(trialdf);
      trialdf = NULL;
      destroyTwoSpectTemplate(template);
      template = NULL;
      XLALDestroyREAL8Vector(trialf);
      trialf = NULL;

   return XLAL_SUCCESS;

}

/**
 * Cluster candidates by frequency, period, and modulation depth using templates
 * \param [out] output        Pointer to pointer of a candidateVector
 * \param [in]  input         Pointer to a candidateVector
 * \param [in]  ffdata        Pointer to ffdataStruct
 * \param [in]  params        Pointer to UserInput_t
 * \param [in]  ffplanenoise  Pointer to REAL4VectorAligned of 2nd FFT background powers
 * \param [in]  fbinaveratios Pointer to REAL4VectorAligned of normalized SFT background
 * \param [in]  rng           Pointer to gsl_rng
 * \param [in]  exactflag     Flag to use Gaussian templates (0) or exact templates (1)
 * \return Status value
 */
INT4 clusterCandidates(candidateVector **output, const candidateVector *input, const ffdataStruct *ffdata, const UserInput_t *params, const REAL4VectorAligned *ffplanenoise, const REAL4VectorAligned *fbinaveratios, const gsl_rng *rng, const BOOLEAN exactflag)
{

   XLAL_CHECK( *output != NULL && input != NULL && ffdata != NULL && params != NULL && ffplanenoise != NULL && fbinaveratios != NULL && rng != NULL, XLAL_EINVAL );

   UINT4 loc, loc2, numcandoutlist;
   REAL8 avefsig, aveperiod, mindf, maxdf;

   //Allocate int vectors for storage
   INT4Vector *locs = NULL, *locs2 = NULL, *usedcandidate = NULL;
   XLAL_CHECK( (locs = XLALCreateINT4Vector(input->numofcandidates)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (locs2 = XLALCreateINT4Vector(input->numofcandidates)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (usedcandidate = XLALCreateINT4Vector(input->numofcandidates)) != NULL, XLAL_EFUNC );

   //Initialize arrays
   for (UINT4 ii=0; ii<input->numofcandidates; ii++) {
      locs->data[ii] = -1;
      locs2->data[ii] = -1;
      usedcandidate->data[ii] = 0;
   }

   //Make FFT plan if exactflag is given
   REAL4FFTPlan *plan = NULL;
   if (exactflag==1) XLAL_CHECK( (plan = XLALCreateForwardREAL4FFTPlan(ffdata->numffts, 1)) != NULL, XLAL_EFUNC );

   numcandoutlist = 0;
   for (UINT4 ii=0; ii<input->numofcandidates; ii++) {

      //Make note of first candidate available
      locs->data[0] = ii;
      loc = 1;

      INT4 foundany = 0;   //Switch to determine if any other candidates in the group. 1 if true
      INT4 iter = 1;
      //Find any in the list that are within +1/2 bin in first FFT frequency
      for (UINT4 jj=ii+1; jj<input->numofcandidates; jj++) {
         if ( usedcandidate->data[jj] == 0 && (input->data[jj].fsig-input->data[locs->data[0]].fsig <= 0.5*iter/params->Tsft+1.0e-6 && input->data[jj].fsig-input->data[locs->data[0]].fsig >= -0.25*iter/params->Tsft) ) {
            locs->data[loc] = jj;
            loc++;
            if (foundany==0) foundany = 1;
         }
      } /* for jj < input->numofcandidates */
      //Keep checking as long as there are more connected frequencies going higher in frequency
      while (foundany==1) {
         foundany = 0;
         iter++;
         for (UINT4 jj=ii+1; jj<input->numofcandidates; jj++) {
            if ( usedcandidate->data[jj] == 0 && (input->data[jj].fsig-input->data[locs->data[0]].fsig-0.25/params->Tsft <= 0.5*iter/params->Tsft && input->data[jj].fsig-input->data[locs->data[0]].fsig+0.25/params->Tsft >= 0.5*iter/params->Tsft) ) {
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
         for (UINT4 jj=ii+1; jj<input->numofcandidates; jj++) {
            if ( usedcandidate->data[jj] == 0 && (input->data[locs->data[0]].fsig-input->data[jj].fsig-0.25/params->Tsft <= 0.5*iter/params->Tsft && input->data[locs->data[0]].fsig-input->data[jj].fsig+0.25/params->Tsft >= 0.5*iter/params->Tsft) ) {
               locs->data[loc] = jj;
               loc++;
               if (foundany==0) foundany = 1;
            }
         } /* for jj < input->numofcandidates */
      } /* while foundany==1 */

      //Using the list of locations, find the subset that have periods within 1 bin
      //of the second FFT frequencies
      UINT4 subsetloc = 0, nextsubsetloc = 0, subsetlocset = 0;
      loc2 = 0;
      for (UINT4 jj=subsetloc; jj<loc; jj++) {
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
            avefsig = 0.0, aveperiod = 0.0, mindf = 0.0, maxdf = 0.0;
            REAL8 weight = 0.0, bestmoddepth = 0.0, bestR = 0.0, besth0 = 0.0, bestProb = 0.0;
            INT4 bestproberrcode = 0;
            for (UINT4 kk=0; kk<loc2; kk++) {
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

            if (loc2 > 1 && aveperiod >= params->Pmin-1.0 && aveperiod <= params->Pmax+1.0) {
               UINT4 numofmoddepths = (UINT4)floorf(2*(maxdf-mindf)*params->Tsft)+1;
               candidate cand;
               TwoSpectTemplate *template = NULL;
               XLAL_CHECK( (template = createTwoSpectTemplate(params->maxTemplateLength)) != NULL, XLAL_EFUNC );

               for (UINT4 kk=0; kk<numofmoddepths; kk++) {
                  if ((mindf+kk*0.5/params->Tsft)>=params->dfmin && (mindf+kk*0.5/params->Tsft)<=params->dfmax) {

                     loadCandidateData(&cand, avefsig, aveperiod, mindf + kk*0.5/params->Tsft, input->data[0].ra, input->data[0].dec, 0, 0, 0.0, 0, 0.0, -1, 0);

                     if (exactflag==1) XLAL_CHECK( makeTemplate(template, cand, params, plan) == XLAL_SUCCESS, XLAL_EFUNC );
                     else XLAL_CHECK( makeTemplateGaussians(template, cand, params) == XLAL_SUCCESS, XLAL_EFUNC );

                     REAL8 R = calculateR(ffdata->ffdata, template, ffplanenoise, fbinaveratios);
                     XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
                     REAL8 prob = probR(template, ffplanenoise, fbinaveratios, R, params, rng, &proberrcode);
                     XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );

                     if (prob < bestProb) {
                        bestmoddepth = mindf + kk*0.5/params->Tsft;
                        bestR = R;
                        bestProb = prob;
                        bestproberrcode = proberrcode;
                     }
                  } /* if test moddepth is within user specified range */
               } /* for kk < numofmoddepths */

               destroyTwoSpectTemplate(template);
               template = NULL;
            } /* if loc2 > 1 ... */

            if (bestProb != 0.0) {
               if (bestR > 0.0) besth0 = 2.7426*pow(bestR/(params->Tsft*params->Tobs),0.25);
               else besth0 = 0.0;

               if ((*output)->numofcandidates == (*output)->length-1) XLAL_CHECK( (*output = resizecandidateVector(*output, 2*(*output)->length)) != NULL, XLAL_EFUNC );
               loadCandidateData(&((*output)->data[(*output)->numofcandidates]), avefsig, aveperiod, bestmoddepth, input->data[0].ra, input->data[0].dec, bestR, besth0, bestProb, bestproberrcode, input->data[0].normalization, -1, 0);
               numcandoutlist++;
               (*output)->numofcandidates++;
            }

            loc2 = 0;
         } /* if jj+1 == loc */
      } /* for jj < loc */

      //Find location of first entry to be searched next time or finish the cluster search
      for (UINT4 jj=ii; jj<input->numofcandidates; jj++) {
         if (usedcandidate->data[jj]==0) {
            ii = jj - 1;
            jj = input->numofcandidates - 1;
         } else if (jj==input->numofcandidates-1) {
            ii = input->numofcandidates - 1;
         }
      }

      //Reinitialize values, just in case
      for (UINT4 jj=0; jj<locs->length; jj++) {
         locs->data[jj] = -1;
         locs2->data[jj] = -1;
      }
   } /* for ii < numofcandidates */

   //Destroy stuff
   XLALDestroyINT4Vector(locs);
   XLALDestroyINT4Vector(locs2);
   XLALDestroyINT4Vector(usedcandidate);
   if (exactflag==1) XLALDestroyREAL4FFTPlan(plan);

   fprintf(stderr, "Clustering done with candidates = %d\n", (*output)->numofcandidates);
   fprintf(LOG, "Clustering done with candidates = %d\n", (*output)->numofcandidates);

   return XLAL_SUCCESS;

} /* clusterCandidates() */


/**
 * Function to test the IHS candidates against Gaussian templates
 * \param [out] output                 Pointer to pointer of a candidateVector
 * \param [in]  ihsCandidates          Pointer to candidateVector of IHS candidates
 * \param [in]  ffdata                 Pointer to ffdataStruct
 * \param [in]  aveNoise               Pointer to REAL4VectorAligned of 2nd FFT background powers
 * \param [in]  aveTFnoisePerFbinRatio Pointer to REAL4VectorAligned of normalized SFT background spectra
 * \param [in]  pos                    The current sky position
 * \param [in]  params                 Pointer to UserInput_t
 * \param [in]  rng                    Pointer to gsl_rng
 * \return Status value
 */
INT4 testIHScandidates(candidateVector **output, const candidateVector *ihsCandidates, const ffdataStruct *ffdata, const REAL4VectorAligned *aveNoise, const REAL4VectorAligned *aveTFnoisePerFbinRatio, const SkyPosition pos, const UserInput_t *params, const gsl_rng *rng)
{

   XLAL_CHECK( *output != NULL && ihsCandidates != NULL && ffdata != NULL && aveNoise != NULL && aveTFnoisePerFbinRatio != NULL && params != NULL && rng != NULL, XLAL_EINVAL );

   //R probability calculator errorcode
   INT4 proberrcode = 0;

   //Allocate memory for FAR struct
   farStruct *farval = NULL;
   XLAL_CHECK( (farval = createfarStruct()) != NULL, XLAL_EFUNC );

   //Allocate memory for template
   TwoSpectTemplate *template = NULL;
   XLAL_CHECK( (template = createTwoSpectTemplate(params->maxTemplateLength)) != NULL, XLAL_EFUNC );

   INT4 candidatesoutsideofmainULrange = 0;
   REAL8 log10templatefar = log10(params->tmplfar);

   for (UINT4 ii=0; ii<ihsCandidates->numofcandidates; ii++) {
      //Assess the IHS candidate if the signal is away from the band edges, the modulation depth is greater or equal to minimum allowed and less than or equal to the maximum allowed, and if the period/modulation depth combo is within allowable limits for a template to be made. We will cut the period space in the next step.
      if ( ihsCandidates->data[ii].fsig>=params->fmin && ihsCandidates->data[ii].fsig<(params->fmin+params->fspan) ) {
         if ( params->followUpOutsideULrange || (ihsCandidates->data[ii].fsig>=params->ULfmin && ihsCandidates->data[ii].fsig<=(params->ULfmin + params->ULfspan) && ihsCandidates->data[ii].moddepth>=params->ULminimumDeltaf && ihsCandidates->data[ii].moddepth<=params->ULmaximumDeltaf) ) {

            resetTwoSpectTemplate(template);

            REAL8 R, prob, bestPeriod = 0.0, bestR = 0.0, bestProb = 0.0;
            INT4 bestproberrcode = 0;

            if (ihsCandidates->data[ii].period>=fmax(4.0*params->Tsft, minPeriod(ihsCandidates->data[ii].moddepth, params->Tsft)) && ihsCandidates->data[ii].period<=(0.2*params->Tobs)) {
               //Make a Gaussian train template
               XLAL_CHECK( makeTemplateGaussians(template, ihsCandidates->data[ii], params) == XLAL_SUCCESS, XLAL_EFUNC );

               //Estimate the FAR for these bin weights if the option was given
               if (params->calcRthreshold) XLAL_CHECK( numericFAR(farval, template, params->tmplfar, aveNoise, aveTFnoisePerFbinRatio, params, rng, params->BrentsMethod) == XLAL_SUCCESS, XLAL_EFUNC );

               //Caclulate R and probability noise caused the candidate
               R = calculateR(ffdata->ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
               XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
               prob = probR(template, aveNoise, aveTFnoisePerFbinRatio, R, params, rng, &proberrcode);
               XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );

               /* Note the candidate if R exceeds the FAR or check other possibilities of different
                periods */
               if ((!params->calcRthreshold && prob<log10templatefar) || (params->calcRthreshold && R>farval->far)) {
                  bestR = R;
                  bestProb = prob;
                  bestPeriod = ihsCandidates->data[ii].period;
               } /* if prob<log10templatefar || R > farval->far */
            } // if within moddepth/period range

            // longer or shorter
            REAL8 periodfact = 0.0;
            for (UINT4 jj=0; jj<=1; jj++) {
               //Shift by harmonics
               for (INT4 kk=2; kk<=params->periodHarmToCheck; kk++) {
                  if (jj==0) periodfact = 1.0/(REAL8)kk;
                  else periodfact = (REAL8)kk;
                  if (ihsCandidates->data[ii].period*periodfact>=fmax(params->Pmin, minPeriod(ihsCandidates->data[ii].moddepth, params->Tsft)) && ihsCandidates->data[ii].period*periodfact<=fmin(params->Pmax, params->Tobs*0.2)) {
                     ihsCandidates->data[ii].period *= periodfact;
                     XLAL_CHECK( makeTemplateGaussians(template, ihsCandidates->data[ii], params) == XLAL_SUCCESS, XLAL_EFUNC );
                     R = calculateR(ffdata->ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
                     XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
                     prob = probR(template, aveNoise, aveTFnoisePerFbinRatio, R, params, rng, &proberrcode);
                     XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
                     if (params->calcRthreshold && bestProb==0.0) XLAL_CHECK( numericFAR(farval, template, params->tmplfar, aveNoise, aveTFnoisePerFbinRatio, params, rng, params->BrentsMethod) == XLAL_SUCCESS, XLAL_EFUNC );
                     if ((bestProb!=0.0 && prob<bestProb) || (bestProb==0.0 && !params->calcRthreshold && prob<log10templatefar) || (bestProb==0.0 && params->calcRthreshold && R>farval->far)) {
                        bestPeriod = ihsCandidates->data[ii].period;
                        bestR = R;
                        bestProb = prob;
                        bestproberrcode = proberrcode;
                     }
                     ihsCandidates->data[ii].period /= periodfact;  //reset the period back to the original value
                  } // in range?
               } // shift by harmonics for kk <= inputParams->periodHarmToCheck (harmonics)

               //shift by fractions
               for (INT4 kk=1; kk<=params->periodFracToCheck; kk++) {
                  if (jj==0) periodfact = (kk+1.0)/(kk+2.0);
                  else periodfact = (kk+2.0)/(kk+1.0);
                  if (ihsCandidates->data[ii].period*periodfact>=fmax(params->Pmin, minPeriod(ihsCandidates->data[ii].moddepth, params->Tsft)) && ihsCandidates->data[ii].period*periodfact<=fmin(params->Pmax, params->Tobs*0.2)) {
                     ihsCandidates->data[ii].period *= periodfact;
                     XLAL_CHECK( makeTemplateGaussians(template, ihsCandidates->data[ii], params) == XLAL_SUCCESS, XLAL_EFUNC );
                     R = calculateR(ffdata->ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
                     XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
                     prob = probR(template, aveNoise, aveTFnoisePerFbinRatio, R, params, rng, &proberrcode);
                     XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
                     if (params->calcRthreshold && bestProb==0.0) XLAL_CHECK( numericFAR(farval, template, params->tmplfar, aveNoise, aveTFnoisePerFbinRatio, params, rng, params->BrentsMethod) == XLAL_SUCCESS, XLAL_EFUNC );
                     if ((bestProb!=0.0 && prob<bestProb) || (bestProb==0.0 && !params->calcRthreshold && prob<log10templatefar) || (bestProb==0.0 && params->calcRthreshold && R>farval->far)) {
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
               if (bestR > 0.0) h0 = 2.7426*sqrt(sqrt(bestR/(params->Tsft*params->Tobs)));  //Now compute the h0 value

               if ((*output)->numofcandidates == (*output)->length-1) XLAL_CHECK( (*output = resizecandidateVector(*output, 2*(*output)->length)) != NULL, XLAL_EFUNC );
               loadCandidateData(&((*output)->data[(*output)->numofcandidates]), ihsCandidates->data[ii].fsig, bestPeriod, ihsCandidates->data[ii].moddepth, pos.longitude, pos.latitude, bestR, h0, bestProb, bestproberrcode, ihsCandidates->data[ii].normalization, -1, ihsCandidates->data[ii].lineContamination);
               (*output)->numofcandidates++;

            } /* if bestR != 0.0, add candidate or replace if something better is found */
         } /* if within UL boundaries */
         else {
            candidatesoutsideofmainULrange++;
         }
      } /* if within outer boundaries */
   } /* for ii < numofcandidates */

   //Destroy allocated memory
   destroyTwoSpectTemplate(template);
   template = NULL;
   destroyfarStruct(farval);
   farval = NULL;

   fprintf(stderr, "%d remaining candidate(s) inside UL range.\n", ihsCandidates->numofcandidates-candidatesoutsideofmainULrange);
   fprintf(stderr,"Initial stage done with candidates = %d\n", (*output)->numofcandidates);
   fprintf(LOG,"Initial stage done with candidates = %d\n", (*output)->numofcandidates);

   return XLAL_SUCCESS;

} /* testIHScandidates() */


/**
 * Test each of the templates in a TwoSpectTemplateVector and keep the top 10
 * This will not check the false alarm probability of any R value less than 0.
 * \param [out] output                 Pointer to pointer of a candidateVector storing a list of all candidates
 * \param [in]  templateVec            Pointer to a TwoSpectTemplateVector containing all the templates to be searched
 * \param [in]  ffdata                 Pointer to ffdataStruct
 * \param [in]  aveNoise               Pointer to REAL4VectorAligned of 2nd FFT background powers
 * \param [in]  aveTFnoisePerFbinRatio Pointer to REAL4VectorAligned of normalized SFT background spectra
 * \param [in]  skypos                 The current sky position
 * \param [in]  params                 Pointer to UserInput_t
 * \param [in]  rng                    Pointer to gsl_rng
 * \param [in]  templateLen            Maximum length of a template
 * \return Status value
 */
INT4 testTwoSpectTemplateVector(candidateVector *output, const TwoSpectTemplateVector *templateVec, const ffdataStruct *ffdata, const REAL4VectorAligned *aveNoise, const REAL4VectorAligned *aveTFnoisePerFbinRatio, const SkyPosition skypos, const UserInput_t *params, const gsl_rng *rng, const UINT4 templateLen)
{

   XLAL_CHECK( output!=NULL && templateVec!=NULL && ffdata!=NULL && aveNoise!=NULL && aveTFnoisePerFbinRatio!=NULL && params!=NULL && rng!=NULL, XLAL_EINVAL );

   fprintf(stderr, "Testing TwoSpectTemplateVector... ");
   
   TwoSpectTemplate *template = NULL;
   XLAL_CHECK( (template = createTwoSpectTemplate(templateLen)) != NULL, XLAL_EFUNC );

   INT4 proberrcode = 0;

   FILE *RVALS = NULL;
   if (XLALUserVarWasSet(&params->saveRvalues)) XLAL_CHECK( (RVALS = fopen(params->saveRvalues, "w")) != NULL, XLAL_EIO, "Couldn't open %s for writing", params->saveRvalues );
   
   UINT4 numfbins = (UINT4)round(params->fspan*params->Tsft);
   for (UINT4 ii=0; ii<numfbins; ii++) {
      REAL8 freq = params->fmin + ii/params->Tsft;
      for (UINT4 jj=0; jj<templateVec->length; jj++) {
         if (templateVec->data[jj]->templatedata->data[0] == 0.0) break;
         
         XLAL_CHECK( convertTemplateForSpecificFbin(template, templateVec->data[jj], freq, params) == XLAL_SUCCESS, XLAL_EFUNC );
         
         REAL8 R = calculateR(ffdata->ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
         XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
         REAL8 prob = 0.0, h0 = 0.0;
         if ( R > 0.0 ) {
            prob = probR(template, aveNoise, aveTFnoisePerFbinRatio, R, params, rng, &proberrcode);
            XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
            h0 = 2.7426*pow(R/(params->Tsft*params->Tobs),0.25);
         }

         if (XLALUserVarWasSet(&params->saveRvalues)) fprintf(RVALS, "%g\n", R);

         if (prob < output->data[output->length-1].prob) {
            UINT4 insertionPoint = output->length - 1;
            while(insertionPoint>0 && prob<output->data[insertionPoint - 1].prob) insertionPoint--;
            for (INT4 kk=(INT4)output->length-2; kk>=(INT4)insertionPoint; kk--) loadCandidateData(&(output->data[kk+1]), output->data[kk].fsig, output->data[kk].period, output->data[kk].moddepth, output->data[kk].ra, output->data[kk].dec, output->data[kk].stat, output->data[kk].h0, output->data[kk].prob, output->data[kk].proberrcode, output->data[kk].normalization, output->data[kk].templateVectorIndex, output->data[kk].lineContamination);
            loadCandidateData(&(output->data[insertionPoint]), template->f0, template->period, template->moddepth, skypos.longitude, skypos.latitude, R, h0, prob, proberrcode, ffdata->tfnormalization, jj, 0);
            if (output->numofcandidates<output->length) output->numofcandidates++;
         }
      }
   }

   destroyTwoSpectTemplate(template);

   if (XLALUserVarWasSet(&params->saveRvalues)) fclose(RVALS);

   fprintf(stderr, "done\n");

   return XLAL_SUCCESS;

}


/**
 * Keep the most significant candidates, potentially reducing the number of candidates if there are more than allowed
 * \param [in] input  Pointer to input candidateVector
 * \param [in] params Pointer to UserInput_t
 * \return Pointer to newly allocated candidateVector containing reduced number of candidates
 */
candidateVector * keepMostSignificantCandidates(const candidateVector *input, const UserInput_t *params)
{

   XLAL_CHECK_NULL( input != NULL && params != NULL, XLAL_EINVAL );

   fprintf(stderr, "Reducing total number of IHS candidates %d to user input %d\n", input->numofcandidates, params->keepOnlyTopNumIHS);
   fprintf(LOG, "Reducing total number of IHS candidates %d to user input %d\n", input->numofcandidates, params->keepOnlyTopNumIHS);

   candidateVector *output = NULL;

   //If the number to keep is > 0 and the number of candidates is less than the number to keep,
   //just move the input vector to the output vector
   if (params->keepOnlyTopNumIHS>0 && (INT4)input->numofcandidates<=params->keepOnlyTopNumIHS) {
      XLAL_CHECK_NULL( (output = createcandidateVector(input->numofcandidates)) != NULL, XLAL_EFUNC );

      for (UINT4 ii=0; ii<input->numofcandidates; ii++) {
         loadCandidateData(&(output->data[ii]), input->data[ii].fsig, input->data[ii].period, input->data[ii].moddepth, input->data[ii].ra, input->data[ii].dec, input->data[ii].stat, input->data[ii].h0, input->data[ii].prob, input->data[ii].proberrcode, input->data[ii].normalization, input->data[ii].templateVectorIndex, input->data[ii].lineContamination);
      }
      output->numofcandidates = input->numofcandidates;

   } else if (params->keepOnlyTopNumIHS>0 && (INT4)input->numofcandidates>params->keepOnlyTopNumIHS) {
      //If keep is > 0 and the number of candidates is > the number to keep,
      //we sort through the list and find the most significant candidates to keep
      XLAL_CHECK_NULL( (output = createcandidateVector(params->keepOnlyTopNumIHS)) != NULL, XLAL_EFUNC );

      for (UINT4 ii=0; ii<output->length; ii++) {
         REAL8 highestsignificance = 0.0;
         INT4 candidateWithHighestSignificance = 0;
         for (UINT4 jj=0; jj<input->numofcandidates; jj++) {
            if (input->data[jj].prob>highestsignificance) {
               highestsignificance = input->data[jj].prob;
               candidateWithHighestSignificance = jj;
            }
         }

         loadCandidateData(&(output->data[ii]), input->data[candidateWithHighestSignificance].fsig, input->data[candidateWithHighestSignificance].period, input->data[candidateWithHighestSignificance].moddepth, input->data[candidateWithHighestSignificance].ra, input->data[candidateWithHighestSignificance].dec, input->data[candidateWithHighestSignificance].stat, input->data[candidateWithHighestSignificance].h0, input->data[candidateWithHighestSignificance].prob, input->data[candidateWithHighestSignificance].proberrcode, input->data[candidateWithHighestSignificance].normalization, input->data[candidateWithHighestSignificance].templateVectorIndex, input->data[candidateWithHighestSignificance].lineContamination);

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


/**
 * Calculate the R statistic from equation 13 of E. Goetz and K. Riles (2011)
 * \param [in] ffdata Pointer to REAL4VectorAligned of the 2nd FFT data
 * \param [in] template Pointer to the template
 * \param [in] noise Pointer to the REAL4VectorAligned containing the background 2nd FFT powers
 * \param [in] fbinaveratios Pointer to the REAL4VectorAligned of normalized SFT background powers
 * \return Value of the R statistic
 */
REAL8 calculateR(const REAL4VectorAligned *ffdata, const TwoSpectTemplate *template, const REAL4VectorAligned *noise, const REAL4VectorAligned *fbinaveratios)
{

   XLAL_CHECK_REAL8( ffdata != NULL && template != NULL && noise != NULL && fbinaveratios != NULL, XLAL_EINVAL );

   REAL8 sumofsqweights = 0.0;
   for (UINT4 ii=0; ii<template->templatedata->length; ii++) if (template->templatedata->data[ii]!=0.0) sumofsqweights += (template->templatedata->data[ii]*template->templatedata->data[ii]);
   XLAL_CHECK_REAL8( sumofsqweights != 0.0, XLAL_EFPDIV0 );
   REAL8 sumofsqweightsinv = 1.0/sumofsqweights;

   UINT4 numfprbins = noise->length;

   REAL8 R = 0.0;
   for (UINT4 ii=0; ii<template->templatedata->length; ii++) {
      if (template->templatedata->data[ii]!=0.0) {
         UINT4 firstfreqbin = template->pixellocations->data[ii]/numfprbins;
         UINT4 secfreqbin = template->pixellocations->data[ii] - firstfreqbin*numfprbins;
         R += (ffdata->data[ template->pixellocations->data[ii] ] - noise->data[secfreqbin]*fbinaveratios->data[firstfreqbin])*template->templatedata->data[ii]*sumofsqweightsinv;
      }
   }

   return R;

} /* calculateR() */

INT4 writeCandidateVector2File(const CHAR *outputfile, const candidateVector *input)
{
   XLAL_CHECK( outputfile != NULL && input != NULL, XLAL_EINVAL );

   FILE *CANDFILE = NULL;
   struct stat XLAL_INIT_DECL(buf);
   if ( stat(outputfile, &buf) == -1 && errno == ENOENT ) {
      XLAL_CHECK( (CANDFILE = fopen(outputfile, "w")) != NULL, XLAL_EIO, "Couldn't fopen file %s to output candidates\n", outputfile );
      fprintf(CANDFILE, "# TwoSpect candidates output file\n");
      fprintf(CANDFILE, "# Freq     Period       Mod. depth  RA    DEC  Statistic val.  Est. h0  False alarm prob.  TF normalization  Template vec. num.  Line contam.\n");
   } else {
      XLAL_CHECK( (CANDFILE = fopen(outputfile, "a")) != NULL, XLAL_EIO, "Couldn't fopen file %s to output candidates\n", outputfile );
   }

   for (UINT4 ii=0; ii<input->numofcandidates; ii++) {
      fprintf(CANDFILE, "%.6f %.6f %.7f %.4f %.4f %.4f %g %.4f %g %d %d\n", input->data[ii].fsig, input->data[ii].period, input->data[ii].moddepth, input->data[ii].ra, input->data[ii].dec, input->data[ii].stat, input->data[ii].h0, input->data[ii].prob, input->data[ii].normalization, input->data[ii].templateVectorIndex, input->data[ii].lineContamination);
   }

   fclose(CANDFILE);

   return XLAL_SUCCESS;
}

/**
 * Calculates maximum modulation depth allowed, equation 6 of E. Goetz and K. Riles (2011)
 * \param [in] period  Orbital period value
 * \param [in] cohtime SFT coherence length
 * \return Maximum modulation depth allowed
 */
REAL8 maxModDepth(const REAL8 period, const REAL8 cohtime)
{
   REAL8 maxB = 0.5*period/(cohtime*cohtime);
   return maxB;
} /* maxModDepth() */


/**
 * Calculates minimum period allowed, equation 6 of E. Goetz and K. Riles (2011)
 * \param [in] moddepth Modulation depth value
 * \param [in] cohtime  SFT coherence length
 * \return Maximum modulation depth allowed
 */
REAL8 minPeriod(const REAL8 moddepth, const REAL8 cohtime)
{
   REAL8 minP = 2.0*moddepth*cohtime*cohtime;
   return minP;
} /* minPeriod() */


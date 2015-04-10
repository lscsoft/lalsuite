/*
*  Copyright (C) 2010 -- 2014 Evan Goetz
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

#include <lal/Window.h>
#include <lal/VectorOps.h>
#include <lal/SinCosLUT.h>

#include "templates.h"
#include "vectormath.h"
#include "statistics.h"
#include "TwoSpectSpecFunc.h"


/**
 * Allocate a new TwoSpectTemplate
 * \param [in] length Length of the template
 * \return Pointer to TwoSpectTemplate
 */
TwoSpectTemplate * new_TwoSpectTemplate(INT4 length)
{

   TwoSpectTemplate *template = NULL;
   XLAL_CHECK_NULL( (template = XLALMalloc(sizeof(*template))) != NULL, XLAL_ENOMEM );

   XLAL_CHECK_NULL( (template->templatedata = XLALCreateREAL4VectorAligned(length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (template->pixellocations = XLALCreateINT4Vector(length)) != NULL, XLAL_EFUNC );

   memset(template->templatedata->data, 0, sizeof(REAL4)*length);
   memset(template->pixellocations->data, 0, sizeof(INT4)*length);

   template->f0 = 0.0;
   template->period = 0.0;
   template->moddepth = 0.0;

   return template;

} /* new_TwoSpectTemplate() */


/**
 * Reset the values in a TwoSpectTemplate
 * \param [in] template Pointer to a TwoSpectTemplate
 */
void resetTwoSpectTemplate(TwoSpectTemplate *template)
{

   INT4 length = (INT4)template->templatedata->length;
   memset(template->templatedata->data, 0, sizeof(REAL4)*length);
   memset(template->pixellocations->data, 0, sizeof(INT4)*length);

   template->f0 = 0.0;
   template->period = 0.0;
   template->moddepth = 0.0;

} /* resetTwoSpectTemplate() */


/**
 * Free a TwoSpectTemplate
 * \param [in] template Pointer to a TwoSpectTemplate
 */
void free_TwoSpectTemplate(TwoSpectTemplate *template)
{
   XLALDestroyREAL4VectorAligned(template->templatedata);
   XLALDestroyINT4Vector(template->pixellocations);
   XLALFree((TwoSpectTemplate*)template);
} /* free_TwoSpectTemplate() */


/**
 * Create a TwoSpectTemplateVector
 * \param [in] length         The number of templates in the vector
 * \param [in] templateLength The maximum number of pixels in a template
 * \return Pointer to a TwoSpectTemplateVector
 */
TwoSpectTemplateVector * new_TwoSpectTemplateVector(UINT4 length, UINT4 templateLength)
{
   TwoSpectTemplateVector *vector = NULL;
   XLAL_CHECK_NULL( (vector = XLALMalloc(sizeof(*vector))) != NULL, XLAL_ENOMEM );

   vector->length = length;
   vector->data = NULL;
   if (vector->length>0) {
      XLAL_CHECK_NULL( (vector->data = XLALMalloc(vector->length*sizeof(*vector->data))) != NULL, XLAL_ENOMEM );
      for (UINT4 ii=0; ii<vector->length; ii++) XLAL_CHECK_NULL( (vector->data[ii] = new_TwoSpectTemplate(templateLength)) != NULL, XLAL_EFUNC );
   }

   vector->Tsft = 0.0;
   vector->SFToverlap = 0.0;
   vector->Tobs = 0.0;

   return vector;
}


/**
 * Free a TwoSpectTemplateVector
 * \param [in] vector Pointer to the TwoSpectTemplateVector
 */
void free_TwoSpectTemplateVector(TwoSpectTemplateVector *vector)
{
   if (vector==NULL) return;
   if ((!vector->length || !vector->data) && (vector->length || vector->data)) XLAL_ERROR_VOID(XLAL_EINVAL);
   if (vector->data) {
      for (UINT4 ii=0; ii<vector->length; ii++) free_TwoSpectTemplate(vector->data[ii]);
      XLALFree((TwoSpectTemplate*)vector->data);
   }
   vector->data = NULL;
   XLALFree((TwoSpectTemplateVector*)vector);
   return;
}


/**
 * Generate a TwoSpectTemplateVector containing the template data
 * \param [in] Pmin              Minimum orbital period (s)
 * \param [in] Pmax              Maximum orbital period (s)
 * \param [in] dfmin             Minimum modulation depth (Hz)
 * \param [in] dfmax             Maximum modulation depth (Hz)
 * \param [in] Tsft              Timespan of SFTs (s)
 * \param [in] SFToverlap        Overlap of SFTs (s)
 * \param [in] Tobs              Observation time (s)
 * \param [in] maxvectorlength   Limit the number of templates generated to be the maximum value specified here
 * \param [in] minTemplateLength The minimum number of pixels in a template
 * \param [in] maxTemplateLength The maximum number of pixels in a template
 * \param [in] vectormathflag    Flag specifying what type of vector math to use: 0 = none, 1 = SSE, 2 = AVX (for SSE and AVX, must compile using correct flags and instructions allowed on the CPU)
 * \param [in] exactflag         Flag specifying to use exact templates: 1 = enabled, 0 = disabled
 * \return Pointer to a TwoSpectTemplateVector
 */
TwoSpectTemplateVector * generateTwoSpectTemplateVector(REAL8 Pmin, REAL8 Pmax, REAL8 dfmin, REAL8 dfmax, REAL8 Tsft, REAL8 SFToverlap, REAL8 Tobs, UINT4 maxvectorlength, UINT4 minTemplateLength, UINT4 maxTemplateLength, UINT4 vectormathflag, BOOLEAN exactflag)
{
   TwoSpectTemplateVector *vector = NULL;
   XLAL_CHECK_NULL( (vector = new_TwoSpectTemplateVector(maxvectorlength, maxTemplateLength)) != NULL, XLAL_EFUNC );

   vector->Tsft = Tsft;
   vector->SFToverlap = SFToverlap;
   vector->Tobs = Tobs;

   UINT4 numdfvals = (UINT4)(floor(2.0*Tsft*(dfmax-dfmin)))+1, numtemplatesgenerated = 0;
   REAL8 alpha0 = 45.0*(Tsft/1800.0)+30.0;

   UINT4 numffts = (UINT4)floor(Tobs/(Tsft-SFToverlap)-1);
   REAL4FFTPlan *plan = NULL;
   XLAL_CHECK_NULL( (plan = XLALCreateForwardREAL4FFTPlan(numffts, 1)) != NULL, XLAL_EFUNC );

   REAL8Vector *dfvals = NULL;
   XLAL_CHECK_NULL( (dfvals = XLALCreateREAL8Vector(numdfvals)) != NULL, XLAL_EFUNC );
   for (UINT4 ii=0; ii<numdfvals; ii++) dfvals->data[ii] = dfmin + 0.5*ii/Tsft;
   for (UINT4 ii=0; ii<numdfvals; ii++) {
      REAL8 P = Pmax;
      while (P>=Pmin && P>=2.0*dfvals->data[ii]*Tsft*Tsft && numtemplatesgenerated<vector->length) {
         if (!exactflag) XLAL_CHECK_NULL( makeTemplateGaussians2(vector->data[numtemplatesgenerated], 0.0, P, dfvals->data[ii], Tsft, SFToverlap, Tobs, minTemplateLength, vectormathflag) == XLAL_SUCCESS, XLAL_EFUNC );
         else XLAL_CHECK_NULL( makeTemplate2(vector->data[numtemplatesgenerated], 0.0, P, dfvals->data[ii], Tsft, SFToverlap, Tobs, minTemplateLength, vectormathflag, plan) == XLAL_SUCCESS, XLAL_EFUNC );
         numtemplatesgenerated++;
         if (numtemplatesgenerated == vector->length) break;
         
         if (!exactflag) XLAL_CHECK_NULL( makeTemplateGaussians2(vector->data[numtemplatesgenerated], 0.5, P, dfvals->data[ii], Tsft, SFToverlap, Tobs, minTemplateLength, vectormathflag) == XLAL_SUCCESS, XLAL_EFUNC );
         else XLAL_CHECK_NULL( makeTemplate2(vector->data[numtemplatesgenerated], 0.5, P, dfvals->data[ii], Tsft, SFToverlap, Tobs, minTemplateLength, vectormathflag, plan) == XLAL_SUCCESS, XLAL_EFUNC );
         numtemplatesgenerated++;
         if (numtemplatesgenerated == vector->length) break;

         REAL8 dP = P*P/(alpha0*Tobs*sqrt(dfvals->data[ii]))*(1+P/(alpha0*Tobs));
         P -= dP;
      }
      if (numtemplatesgenerated == vector->length) break;
   }

   XLALDestroyREAL8Vector(dfvals);
   XLALDestroyREAL4FFTPlan(plan);

   fprintf(stderr, "Templates generated = %d\n", numtemplatesgenerated);
   return vector;
}


/**
 * Write a TwoSpectTemplateVector to binary file
 * \param [in] vector   Pointer to the TwoSpectTemplateVector
 * \param [in] filename String of the filename
 * \return Status value
 */
INT4 writeTwoSpectTemplateVector(TwoSpectTemplateVector *vector, CHAR *filename)
{
   FILE *fp = NULL;
   XLAL_CHECK( (fp = fopen(filename, "wb")) != NULL, XLAL_EIO, "Couldn't fopen %s for writing\n", filename );
   fwrite(&(vector->length), sizeof(UINT4), 1, fp);
   fwrite(&(vector->data[0]->templatedata->length), sizeof(UINT4), 1, fp);
   fwrite(&(vector->Tsft), sizeof(REAL8), 1, fp);
   fwrite(&(vector->SFToverlap), sizeof(REAL8), 1, fp);
   fwrite(&(vector->Tobs), sizeof(REAL8), 1, fp);
   for (UINT4 ii=0; ii<vector->length; ii++) {
      fwrite(vector->data[ii]->templatedata->data, sizeof(REAL4), vector->data[ii]->templatedata->length, fp);
      fwrite(vector->data[ii]->pixellocations->data, sizeof(INT4), vector->data[ii]->pixellocations->length, fp);
      fwrite(&(vector->data[ii]->f0), sizeof(REAL8), 1, fp);
      fwrite(&(vector->data[ii]->period), sizeof(REAL8), 1, fp);
      fwrite(&(vector->data[ii]->moddepth), sizeof(REAL8), 1, fp);
   }
   fclose(fp);

   return XLAL_SUCCESS;
}


/**
 * Read a TwoSpectTemplateVector from a binary file
 * \param [in] filename String of the filename
 * \return Pointer to a new TwoSpectTemplateVector
 */
TwoSpectTemplateVector * readTwoSpectTemplateVector(CHAR *filename)
{
   TwoSpectTemplateVector *vector = NULL;
   UINT4 vectorlength, templatelength;
   FILE *fp = NULL;
   XLAL_CHECK_NULL( (fp = fopen(filename, "rb")) != NULL, XLAL_EIO, "Couldn't fopen %s for reading\n", filename );
   
   fread(&vectorlength, sizeof(UINT4), 1, fp);
   fread(&templatelength, sizeof(UINT4), 1, fp);
   XLAL_CHECK_NULL( (vector = new_TwoSpectTemplateVector(vectorlength, templatelength)) != NULL, XLAL_EFUNC );
   fread(&(vector->Tsft), sizeof(REAL8), 1, fp);
   fread(&(vector->SFToverlap), sizeof(REAL8), 1, fp);
   fread(&(vector->Tobs), sizeof(REAL8), 1, fp);
   for (UINT4 ii=0; ii<vectorlength; ii++) {
      fread(vector->data[ii]->templatedata->data, sizeof(REAL4), templatelength, fp);
      fread(vector->data[ii]->pixellocations->data, sizeof(INT4), templatelength, fp);
      fread(&(vector->data[ii]->f0), sizeof(REAL8), 1, fp);
      fread(&(vector->data[ii]->period), sizeof(REAL8), 1, fp);
      fread(&(vector->data[ii]->moddepth), sizeof(REAL8), 1, fp);
   }
   fclose(fp);
   
   return vector;
}


/**
 * \brief Make an estimated template based on FFT of train of Gaussians
 *
 * This is eq. 18 of E. Goetz and K. Riles (2011). This handles spillages of power into neighboring frequency bins a little
 * more gracefully. Numerical stability issues mean that we need to compute exp(log(eq. 18)) = eq. 18
 * exp(log(eq. 18)) = exp(log(4*pi*sigma^2) - sigma^2*omegapr^2) * (1+cos(delta*omegapr)) * exp(log(1-cos(N*P*omegapr))-log(P*omegapr))
 * \param [out] output     Pointer to TwoSpectTemplate
 * \param [in]  input      An input candidate structure
 * \param [in]  params     Pointer to UserInput_t
 * \return Status value
 */
INT4 makeTemplateGaussians(TwoSpectTemplate *output, candidate input, UserInput_t *params)
{
   XLAL_CHECK( output != NULL && params != NULL, XLAL_EINVAL );
   XLAL_CHECK( input.period != 0.0 && input.moddepth != 0.0, XLAL_EINVAL, "Invalid input (%f, %f, %f)\n", input.fsig, input.period, input.moddepth );
   
   REAL8 freqbin = input.fsig*params->Tsft;
   INT4 roundedbinval = (INT4)round(freqbin);
   REAL8 offset = freqbin - roundedbinval;
   
   INT4 ffdatabin0 = (INT4)round((params->fmin-params->dfmax)*params->Tsft) - 6;
   INT4 sigbin0 = roundedbinval - ffdatabin0;

   UINT4 numffts = (UINT4)floor(params->Tobs/(params->Tsft-params->SFToverlap)-1);
   UINT4 numfprbins = (UINT4)floorf(0.5*numffts) + 1;
   
   XLAL_CHECK( makeTemplateGaussians2(output, offset, input.period, input.moddepth, params->Tsft, params->SFToverlap, params->Tobs, params->minTemplateLength, (UINT4)params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );

   for (UINT4 ii=0; ii<output->pixellocations->length; ii++) {
      output->pixellocations->data[ii] += sigbin0*numfprbins;
   }
   
   //makeTemplateGaussians2 sets moddepth and P but not fsig, so we do that here
   output->f0 = input.fsig;

   return XLAL_SUCCESS;
   
}


/**
 * \brief Convert an arbitrary frequency bin template into a template for a specific frequency bin.
 *
 * When using this function, the input template should be generated first using mateTemplate*2 functions. 
 * Then call this function with a specific bin center frequency. 
 * The input template can contain an offset up to half a frequency bin in order to generate any frequency of a template.
 * \param [in,out] output Pointer to a TwoSpectTemplate for a specific frequency bin
 * \param [in]     input  Pointer to a TwoSpectTemplate for an arbitrary frequency bin
 * \param [in]     freq   Frequency of a bin centered signal
 * \param [in]     params Pointer to a UserInput_t struct
 */
INT4 convertTemplateForSpecificFbin(TwoSpectTemplate *output, TwoSpectTemplate *input, REAL8 freq, UserInput_t *params)
{
   XLAL_CHECK( output != NULL && input != NULL && params != NULL && freq>=params->fmin && freq<params->fmin+params->fspan, XLAL_EINVAL );

   resetTwoSpectTemplate(output);
   
   UINT4 numffts = (UINT4)floor(params->Tobs/(params->Tsft-params->SFToverlap)-1);
   UINT4 numfprbins = (UINT4)floorf(0.5*numffts) + 1;

   memcpy(output->templatedata->data, input->templatedata->data, sizeof(REAL4)*input->templatedata->length);
   memcpy(output->pixellocations->data, input->pixellocations->data, sizeof(INT4)*input->pixellocations->length);
   output->f0 = freq + input->f0/params->Tsft;
   output->period = input->period;
   output->moddepth = input->moddepth;
   
   REAL8 freqbin = freq*params->Tsft;
   INT4 roundedbinval = (INT4)round(freqbin);
   INT4 ffdatabin0 = (INT4)round((params->fmin-params->dfmax)*params->Tsft) - 6;
   INT4 sigbin0 = roundedbinval - ffdatabin0;
   for (UINT4 ii=0; ii<input->pixellocations->length; ii++) output->pixellocations->data[ii] += sigbin0*numfprbins;

   return XLAL_SUCCESS;
}


/**
 * \brief Make an estimated template based on FFT of train of Gaussians
 *
 * This is eq. 18 of E. Goetz and K. Riles (2011). This handles spillages of power into neighboring frequency bins a little
 * more gracefully. Numerical stability issues mean that we need to compute exp(log(eq. 18)) = eq. 18
 * exp(log(eq. 18)) = exp(log(4*pi*sigma^2) - sigma^2*omegapr^2) * (1+cos(delta*omegapr)) * exp(log(1-cos(N*P*omegapr))-log(P*omegapr))
 * \param [in,out] output            Pointer to TwoSpectTemplate
 * \param [in]     offset            Amount of offset from bin centered signal (-0.5 <= offset <= 0.5)
 * \param [in]     P                 Orbital period (seconds)
 * \param [in]     deltaf            Modulation depth of the signal (Hz)
 * \param [in]     Tsft              Length of an SFT (s)
 * \param [in]     SFToverlap        SFT overlap (s)
 * \param [in]     Tobs              Observation time (s)
 * \param [in]     minTemplateLength Minimum number of pixels in a template
 * \param [in]     vectormathflag    Flag indicating to use vector math: 0 = none, 1 = SSE, 2 = AVX (must compile for vector math appropriately and CPU must have those instructions)
 * \return Status value
 */
INT4 makeTemplateGaussians2(TwoSpectTemplate *output, REAL8 offset, REAL8 P, REAL8 deltaf, REAL8 Tsft, REAL8 SFToverlap, REAL8 Tobs, UINT4 minTemplateLength, UINT4 vectormathflag)
{

   XLAL_CHECK( output != NULL, XLAL_EINVAL );
   XLAL_CHECK( offset <= 0.5 && offset >= -0.5 && P != 0.0 && deltaf != 0.0, XLAL_EINVAL, "Invalid input (%f, %f, %f)\n", offset, P, deltaf );

   UINT4 numffts = (UINT4)floor(Tobs/(Tsft-SFToverlap)-1);
   UINT4 numfprbins = (UINT4)floorf(0.5*numffts) + 1;
   
   //Set data for output template
   output->f0 = offset;
   output->period = P;
   output->moddepth = deltaf;

   //Reset the data values to zero, just in case
   memset(output->templatedata->data, 0, sizeof(REAL4)*output->templatedata->length);

   INT4 N = (INT4)floor(Tobs/P);     //Number of Gaussians = observation time / period
   REAL8 periodf = 1.0/P;

   //Create second FFT frequencies and other useful values
   REAL4VectorAligned *fpr = NULL;
   XLAL_CHECK( (fpr = XLALCreateREAL4VectorAligned(numfprbins, 32)) != NULL, XLAL_EFUNC );
   for (UINT4 ii=0; ii<fpr->length; ii++) fpr->data[ii] = (REAL4)ii*(1.0/Tobs);

   //For speed, we will precompute a number of useful vectors described by their names
   //This part is the allocation
   REAL4VectorAligned *omegapr = NULL, *omegapr_squared = NULL, *cos_ratio = NULL;
   XLAL_CHECK( (omegapr = XLALCreateREAL4VectorAligned(fpr->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (omegapr_squared = XLALCreateREAL4VectorAligned(fpr->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (cos_ratio = XLALCreateREAL4VectorAligned(fpr->length, 32)) != NULL, XLAL_EFUNC );

   //Doing the precomputation of the useful values
   if (vectormathflag==1) {
      XLAL_CHECK( sseScaleREAL4Vector(omegapr, fpr, (REAL4)LAL_TWOPI) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( sseSSVectorMultiply(omegapr_squared, omegapr, omegapr) == XLAL_SUCCESS, XLAL_EFUNC );
   } else if (vectormathflag==2) {
      XLAL_CHECK( avxScaleREAL4Vector(omegapr, fpr, (REAL4)LAL_TWOPI) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( avxSSVectorMultiply(omegapr_squared, omegapr, omegapr) == XLAL_SUCCESS, XLAL_EFUNC );
   } else {
      for (UINT4 ii=0; ii<fpr->length; ii++) {
         omegapr->data[ii] = (REAL4)LAL_TWOPI*fpr->data[ii];
         omegapr_squared->data[ii] = omegapr->data[ii]*omegapr->data[ii];
      }
   }
   if (N*P*omegapr->data[omegapr->length-1]<2.147483647e9) {
      for (UINT4 ii=0; ii<fpr->length; ii++) {
         REAL4 tempSinValue = 0.0, cos_omegapr_times_period = 0.0, cos_N_times_omegapr_times_period = 0.0;
         XLAL_CHECK( XLALSinCosLUT(&tempSinValue, &cos_omegapr_times_period, P*omegapr->data[ii]) == XLAL_SUCCESS, XLAL_EFUNC );
         XLAL_CHECK( XLALSinCosLUT(&tempSinValue, &cos_N_times_omegapr_times_period, N*P*omegapr->data[ii]) == XLAL_SUCCESS, XLAL_EFUNC );
         if (cos_N_times_omegapr_times_period>1.0) cos_N_times_omegapr_times_period = 1.0;
         if (cos_omegapr_times_period>1.0) cos_omegapr_times_period = 1.0;
         if (cos_N_times_omegapr_times_period<-1.0) cos_N_times_omegapr_times_period = -1.0;
         if (cos_omegapr_times_period<-1.0) cos_omegapr_times_period = -1.0;
         if (cos_N_times_omegapr_times_period<=(1.0-100.0*LAL_REAL4_EPS) && cos_omegapr_times_period<=(1.0-100.0*LAL_REAL4_EPS)) cos_ratio->data[ii] = (1.0 - cos_N_times_omegapr_times_period)/(1.0 - cos_omegapr_times_period);
         else if (cos_N_times_omegapr_times_period<=(1.0-100.0*LAL_REAL4_EPS)) {
            REAL8 fmodval = fmod(P*omegapr->data[ii], LAL_TWOPI);
            cos_ratio->data[ii] = 2.0*(1.0 - cos_N_times_omegapr_times_period)/(fmodval*fmodval);
         } else if (cos_omegapr_times_period<=(1.0-100.0*LAL_REAL4_EPS)) {
            REAL8 fmodval = fmod(N*P*omegapr->data[ii], LAL_TWOPI);
            cos_ratio->data[ii] = 0.5*fmodval*fmodval/(1.0 - cos_omegapr_times_period);
         } else cos_ratio->data[ii] = (REAL4)(N*N);
      }
   } else if (P*omegapr->data[omegapr->length-1]<2.147483647e9) {
      for (UINT4 ii=0; ii<fpr->length; ii++) {
         REAL4 tempSinValue = 0.0, cos_omegapr_times_period = 0.0;
         XLAL_CHECK( XLALSinCosLUT(&tempSinValue, &cos_omegapr_times_period, P*omegapr->data[ii]) == XLAL_SUCCESS, XLAL_EFUNC );
         if (cos_omegapr_times_period>1.0) cos_omegapr_times_period = 1.0;
         if (cos_omegapr_times_period<-1.0) cos_omegapr_times_period = -1.0;
         REAL4 cos_N_times_omegapr_times_period = cosf((REAL4)(N*P*omegapr->data[ii]));
         if (cos_N_times_omegapr_times_period<=(1.0-100.0*LAL_REAL4_EPS) && cos_omegapr_times_period<=(1.0-100.0*LAL_REAL4_EPS)) cos_ratio->data[ii] = (1.0 - cos_N_times_omegapr_times_period)/(1.0 - cos_omegapr_times_period);
         else if (cos_N_times_omegapr_times_period<=(1.0-100.0*LAL_REAL4_EPS)) {
            REAL8 fmodval = fmod(P*omegapr->data[ii], LAL_TWOPI);
            cos_ratio->data[ii] = 2.0*(1.0 - cos_N_times_omegapr_times_period)/(fmodval*fmodval);
         } else if (cos_omegapr_times_period<=(1.0-100.0*LAL_REAL4_EPS)) {
            REAL8 fmodval = fmod(N*P*omegapr->data[ii], LAL_TWOPI);
            cos_ratio->data[ii] = 0.5*fmodval*fmodval/(1.0 - cos_omegapr_times_period);
         } else cos_ratio->data[ii] = (REAL4)(N*N);
      }
   } else {
      for (UINT4 ii=0; ii<fpr->length; ii++) {
         REAL4 cos_omegapr_times_period = cosf((REAL4)(P*omegapr->data[ii]));
         REAL4 cos_N_times_omegapr_times_period = cosf((REAL4)(N*P*omegapr->data[ii]));
         if (cos_N_times_omegapr_times_period<=(1.0-100.0*LAL_REAL4_EPS) && cos_omegapr_times_period<=(1.0-100.0*LAL_REAL4_EPS)) cos_ratio->data[ii] = (1.0 - cos_N_times_omegapr_times_period)/(1.0 - cos_omegapr_times_period);
         else if (cos_N_times_omegapr_times_period<=(1.0-100.0*LAL_REAL4_EPS)) {
            REAL8 fmodval = fmod(P*omegapr->data[ii], LAL_TWOPI);
            cos_ratio->data[ii] = 2.0*(1.0 - cos_N_times_omegapr_times_period)/(fmodval*fmodval);
         } else if (cos_omegapr_times_period<=(1.0-100.0*LAL_REAL4_EPS)) {
            REAL8 fmodval = fmod(N*P*omegapr->data[ii], LAL_TWOPI);
            cos_ratio->data[ii] = 0.5*fmodval*fmodval/(1.0 - cos_omegapr_times_period);
         } else cos_ratio->data[ii] = (REAL4)(N*N);
      }
   }

   //Determine span of the template
   REAL8 binamplitude = deltaf*Tsft;
   REAL8 binmin = -binamplitude + offset;
   REAL8 binmax = binamplitude + offset;
   INT4 templatemin = (INT4)round(binmin), templatemax = (INT4)round(binmax);
   if (templatemin > binmin) templatemin--;
   if (templatemin - binmin >= -0.5) templatemin--;
   if (templatemax < binmax) templatemax++;
   if (templatemax - binmax <= 0.5) templatemax++;
   UINT4 templatespan = (UINT4)(templatemax - templatemin) + 1;
   REAL8 disttominbin = binmin - templatemin, disttomaxbin = templatemax - binmax, disttomidbin = offset - templatemin;

   //Determine the weighting scale for the leakage bins and the distance between Gaussians := phi_actual
   REAL4Vector *scale = NULL, *phi_actual = NULL;
   XLAL_CHECK( (scale = XLALCreateREAL4Vector(templatespan)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (phi_actual = XLALCreateREAL4Vector(templatespan)) != NULL, XLAL_EFUNC );
   memset(scale->data, 0, templatespan*sizeof(REAL4));
   memset(phi_actual->data, 0, templatespan*sizeof(REAL4));
   for (UINT4 ii=0; ii<scale->length; ii++) {
      if (ii!=0 && ii!=scale->length-1) {
         scale->data[ii] = 1.0;
      } else if (ii==0) {
         scale->data[ii] = sqsincxoverxsqminusone(disttominbin);
         XLAL_CHECK( xlalErrno==0, XLAL_EFUNC );
      } else {
         scale->data[ii] = sqsincxoverxsqminusone(disttomaxbin);
         XLAL_CHECK( xlalErrno==0, XLAL_EFUNC );
      }

      if ( fabs(ii-disttomidbin)/(deltaf*Tsft) <= 1.0 ) phi_actual->data[ii] = 0.5*P - asin(fabs(ii-disttomidbin)/(deltaf*Tsft))*LAL_1_PI*P;
   }

   //Make sigmas for each frequency
   //First, allocate vectors
   REAL4Vector *sigmas = NULL, *wvals = NULL, *allsigmas = NULL, *weightvals = NULL;
   XLAL_CHECK( (sigmas = XLALCreateREAL4Vector(templatespan)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (wvals = XLALCreateREAL4Vector((UINT4)floor(30.0*P/Tsft))) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (allsigmas = XLALCreateREAL4Vector(wvals->length * sigmas->length)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (weightvals = XLALCreateREAL4Vector(wvals->length * sigmas->length)) != NULL, XLAL_EFUNC );

   //Here is where the sigmas are computed. It is a weighted average. t = (ii+1)*in->Tsft*0.5
   REAL4 sin2pix = 0.0, cos2pix = 0.0;
   for (UINT4 ii=0; ii<wvals->length; ii++) {
      //calculate sin and cos of 2*pi*t/P and then the bin the signal is in and the signal velocity
      XLAL_CHECK( XLALSinCos2PiLUT(&sin2pix, &cos2pix, periodf*((ii+1)*Tsft*0.5)) == XLAL_SUCCESS, XLAL_EFUNC );
      if (cos2pix>1.0) cos2pix = 1.0;
      else if (cos2pix<-1.0) cos2pix = -1.0;
      REAL4 sigbin = (deltaf*cos2pix)*Tsft + offset;
      REAL4 sigbinvelocity = fabs(-deltaf*sin2pix*Tsft*Tsft*LAL_PI*periodf);

      //if the velocity approaches zero, the sigma calculation will diverge (which it should do) but this is bad numerically, so we cap it
      if (sigbinvelocity<1.0e-4) sigbinvelocity = 1.0e-4;

      //REAL4 sigma = 0.5*Tsft * (0.5346 * powf(sigbinvelocity, -1.0213f));   //Derived fit from simulation
      REAL4 sigma = 0.5*Tsft * (0.5979 / (sigbinvelocity - 3.2895e-5));  //Could think about using this fit in the future

      for (UINT4 jj=0; jj<sigmas->length; jj++) {
         REAL4 bindiff = sigbin-((REAL4)jj+templatemin);
         if (fabs(bindiff)<1.75) {
            weightvals->data[ii*sigmas->length + jj] = sqsincxoverxsqminusone(bindiff);
            XLAL_CHECK( xlalErrno==0, XLAL_EFUNC );
         } else weightvals->data[ii*sigmas->length + jj] = 0.0;
         allsigmas->data[ii*sigmas->length + jj] = weightvals->data[ii*sigmas->length + jj]*sigma;
      }

   } /* for ii < wvals->length */
   for (UINT4 ii=0; ii<sigmas->length; ii++) {
      REAL8 wavesigma = 0.0;
      REAL8 totalw = 0.0;
      for (UINT4 jj=0; jj<wvals->length; jj++) {
         wavesigma += allsigmas->data[ii + jj*sigmas->length];
         totalw += weightvals->data[ii + jj*sigmas->length];
      }
      sigmas->data[ii] = (REAL4)(wavesigma/totalw);
   } /* for ii < sigmas->length */

   //Allocate more useful data vectors. These get computed for each different first FFT frequency bin in the F-F plane
   REAL4VectorAligned *exp_neg_sigma_sq_times_omega_pr_sq = NULL, *sin_phi_times_omega_pr = NULL, *cos_phi_times_omega_pr = NULL, *phi_times_fpr = NULL, *datavector = NULL;
   XLAL_CHECK( (exp_neg_sigma_sq_times_omega_pr_sq = XLALCreateREAL4VectorAligned(omegapr_squared->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (sin_phi_times_omega_pr = XLALCreateREAL4VectorAligned(omegapr->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (cos_phi_times_omega_pr = XLALCreateREAL4VectorAligned(omegapr->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (phi_times_fpr = XLALCreateREAL4VectorAligned(fpr->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (datavector = XLALCreateREAL4VectorAligned(fpr->length, 32)) != NULL, XLAL_EFUNC );

   //Create template. We are going to do exp(log(Eq. 18))
   REAL8 sum = 0.0;
   REAL4 log4pi = 2.53102424697f;
   for (UINT4 ii=0; ii<sigmas->length; ii++) {
      INT4 bins2middlebin = ii + templatemin;  //Number of frequency bins away from the "central frequency bin" of the template

      memset(datavector->data, 0, sizeof(REAL4)*datavector->length);

      //Scaling factor for leakage
      REAL4 scale1 = sqrtf((REAL4)(1.0/(1.0+expf((REAL4)(-phi_actual->data[ii]*phi_actual->data[ii]*0.5/(sigmas->data[ii]*sigmas->data[ii]))))));

      //pre-factor
      //REAL8 prefact0 = scale1 * 2.0 * LAL_TWOPI * sigmas->data[ii] * sigmas->data[ii];
      //REAL4 prefact0 = log(scale1 * 2.0 * LAL_TWOPI * sigmas->data[ii] * sigmas->data[ii]);     //We are going to do exp(log(Eq. 18))
      REAL4 prefact0 = log4pi + 2.0*logf((REAL4)(scale1*sigmas->data[ii]));

      if (vectormathflag==1 || vectormathflag==2) {
         //Compute exp(log(4*pi*s*s*exp(-s*s*omegapr_squared))) = exp(log(4*pi*s*s)-s*s*omegapr_squared)
         if (vectormathflag==1) {
            XLAL_CHECK( sseScaleREAL4Vector(exp_neg_sigma_sq_times_omega_pr_sq, omegapr_squared, -sigmas->data[ii]*sigmas->data[ii]) == XLAL_SUCCESS, XLAL_EFUNC );
            XLAL_CHECK( sseAddScalarToREAL4Vector(exp_neg_sigma_sq_times_omega_pr_sq, exp_neg_sigma_sq_times_omega_pr_sq, prefact0) == XLAL_SUCCESS, XLAL_EFUNC );
         } else {
            XLAL_CHECK( avxScaleREAL4Vector(exp_neg_sigma_sq_times_omega_pr_sq, omegapr_squared, -sigmas->data[ii]*sigmas->data[ii]) == XLAL_SUCCESS, XLAL_EFUNC );
            XLAL_CHECK( avxAddScalarToREAL4Vector(exp_neg_sigma_sq_times_omega_pr_sq, exp_neg_sigma_sq_times_omega_pr_sq, prefact0) == XLAL_SUCCESS, XLAL_EFUNC );
         }
         UINT4 truncationLength = 1;
         while (truncationLength<omegapr_squared->length && exp_neg_sigma_sq_times_omega_pr_sq->data[truncationLength]>-88.0) truncationLength++;
         XLAL_CHECK( XLALVectorExpREAL4(exp_neg_sigma_sq_times_omega_pr_sq->data, exp_neg_sigma_sq_times_omega_pr_sq->data, truncationLength) == XLAL_SUCCESS, XLAL_EFUNC );

         //Compute phi_actual*fpr
         if (vectormathflag==1) XLAL_CHECK( truncated_sseScaleREAL4Vector(phi_times_fpr, fpr, phi_actual->data[ii], truncationLength) == XLAL_SUCCESS, XLAL_EFUNC );
         else XLAL_CHECK( truncated_avxScaleREAL4Vector(phi_times_fpr, fpr, phi_actual->data[ii], truncationLength) == XLAL_SUCCESS, XLAL_EFUNC );

         //Start computing the datavector values
         INT4 maxindex = max_index_in_range((REAL4Vector*)phi_times_fpr, 0, truncationLength-1);
         if (phi_times_fpr->data[maxindex]<=2.147483647e9) {
            //Compute cos(2*pi*phi_actual*fpr) using LUT and SSE
            XLAL_CHECK( XLALVectorSinCos2PiREAL4(sin_phi_times_omega_pr->data, cos_phi_times_omega_pr->data, phi_times_fpr->data, truncationLength) == XLAL_SUCCESS, XLAL_EFUNC );
         } else {
            //Compute cos(2*pi*phi_actual*fpr) without using LUT
            for (UINT4 jj=0; jj<truncationLength; jj++) cos_phi_times_omega_pr->data[jj] = cosf((REAL4)LAL_TWOPI*phi_times_fpr->data[jj]);
         }
         //datavector = cos(phi_actual*omega_pr) + 1.0
         if (vectormathflag==1) XLAL_CHECK( truncated_sseAddScalarToREAL4Vector(datavector, cos_phi_times_omega_pr, 1.0, truncationLength) == XLAL_SUCCESS, XLAL_EFUNC );
         else XLAL_CHECK( truncated_avxAddScalarToREAL4Vector(datavector, cos_phi_times_omega_pr, 1.0, truncationLength) == XLAL_SUCCESS, XLAL_EFUNC );
         //datavector = prefact0 * exp(-s*s*omega_pr*omega_pr) * [cos(phi_actual*omega_pr) + 1.0]
         if (vectormathflag==1) XLAL_CHECK( truncated_sseSSVectorMultiply(datavector, datavector, exp_neg_sigma_sq_times_omega_pr_sq, truncationLength) == XLAL_SUCCESS, XLAL_EFUNC );
         else XLAL_CHECK( truncated_avxSSVectorMultiply(datavector, datavector, exp_neg_sigma_sq_times_omega_pr_sq, truncationLength) == XLAL_SUCCESS, XLAL_EFUNC );
         //datavector = scale * exp(-s*s*omega_pr*omega_pr) * [cos(phi_actual*omega_pr) + 1.0]
         if (vectormathflag==1) XLAL_CHECK( truncated_sseScaleREAL4Vector(datavector, datavector, scale->data[ii], truncationLength) == XLAL_SUCCESS, XLAL_EFUNC );
         else XLAL_CHECK( truncated_avxScaleREAL4Vector(datavector, datavector, scale->data[ii], truncationLength) == XLAL_SUCCESS, XLAL_EFUNC );
         //datavector *= cos_ratio
         if (vectormathflag==1) XLAL_CHECK( truncated_sseSSVectorMultiply(datavector, datavector, cos_ratio, truncationLength) == XLAL_SUCCESS, XLAL_EFUNC );
         else XLAL_CHECK( truncated_avxSSVectorMultiply(datavector, datavector, cos_ratio, truncationLength) == XLAL_SUCCESS, XLAL_EFUNC );
      } else {
         for (UINT4 jj=0; jj<omegapr_squared->length; jj++) {
            //Do all or nothing if the exponential is too negative
            if ((prefact0-sigmas->data[ii]*sigmas->data[ii]*omegapr_squared->data[jj])>-88.0) {
               exp_neg_sigma_sq_times_omega_pr_sq->data[jj] = expf((REAL4)(prefact0-sigmas->data[ii]*sigmas->data[ii]*omegapr_squared->data[jj]));
               XLAL_CHECK( XLALSinCos2PiLUT(&sin2pix, &cos2pix, phi_actual->data[ii]*fpr->data[jj]) == XLAL_SUCCESS, XLAL_EFUNC );
               if (cos2pix>1.0) cos2pix = 1.0;
               else if (cos2pix<-1.0) cos2pix = -1.0;
               cos_phi_times_omega_pr->data[jj] = (REAL4)cos2pix;
               datavector->data[jj] = scale->data[ii]*exp_neg_sigma_sq_times_omega_pr_sq->data[jj]*(cos_phi_times_omega_pr->data[jj]+1.0)*cos_ratio->data[jj];
            }
            /* Don't need to compute the else because the datavector was set to zero in the outer for loop */
         } /* for jj = 0 --> omegapr_squared->length */
      } /* use SSE or not */

      //Now loop through the second FFT frequencies, starting with index 4
      for (UINT4 jj=4; jj<omegapr->length; jj++) {
         //Sum up the weights in total
         sum += (REAL8)(datavector->data[jj]);

         //Compare with weakest top bins and if larger, launch a search to find insertion spot (insertion sort)
         if (datavector->data[jj] > output->templatedata->data[output->templatedata->length-1]) {
            //insertionSort_template(output, datavector->data[jj], bins2middlebin*fpr->length+jj, bins2middlebin, jj);
            insertionSort_template(output, datavector->data[jj], bins2middlebin*fpr->length+jj);
         }
      } /* for jj < omegapr->length */
   } /* for ii < sigmas->length */

   //Normalize
   REAL4 invsum = (REAL4)(1.0/sum);
   if (vectormathflag==1) XLAL_CHECK( sseScaleREAL4Vector(output->templatedata, output->templatedata, invsum) == XLAL_SUCCESS, XLAL_EFUNC );
   else if (vectormathflag==2) XLAL_CHECK( avxScaleREAL4Vector(output->templatedata, output->templatedata, invsum) == XLAL_SUCCESS, XLAL_EFUNC );
   else for (UINT4 ii=0; ii<output->templatedata->length; ii++) if (output->templatedata->data[ii]!=0.0) output->templatedata->data[ii] *= invsum;

   //Truncate weights when they don't add much to the total sum of weights
   sum = 0.0;
   for (UINT4 ii=0; ii<minTemplateLength; ii++) sum += (REAL8)output->templatedata->data[ii];
   UINT4 counter = minTemplateLength;
   while (counter<output->templatedata->length && output->templatedata->data[counter]>=epsval_float((REAL4)sum)) {
      sum += (REAL8)output->templatedata->data[counter];
      counter++;
   }
   for (/* last counter val */; counter<output->templatedata->length; counter++) output->templatedata->data[counter] = 0.0;

   //Destroy variables
   XLALDestroyREAL4Vector(phi_actual);
   XLALDestroyREAL4Vector(scale);
   XLALDestroyREAL4Vector(sigmas);
   XLALDestroyREAL4Vector(allsigmas);
   XLALDestroyREAL4Vector(weightvals);
   XLALDestroyREAL4Vector(wvals);
   XLALDestroyREAL4VectorAligned(fpr);
   XLALDestroyREAL4VectorAligned(omegapr);
   XLALDestroyREAL4VectorAligned(omegapr_squared);
   XLALDestroyREAL4VectorAligned(cos_ratio);
   XLALDestroyREAL4VectorAligned(exp_neg_sigma_sq_times_omega_pr_sq);
   XLALDestroyREAL4VectorAligned(phi_times_fpr);
   XLALDestroyREAL4VectorAligned(sin_phi_times_omega_pr);
   XLALDestroyREAL4VectorAligned(cos_phi_times_omega_pr);
   XLALDestroyREAL4VectorAligned(datavector);

   return XLAL_SUCCESS;

} /* makeTemplateGaussians() */


INT4 truncated_sseScaleREAL4Vector(REAL4VectorAligned *output, REAL4VectorAligned *input, REAL4 scale, UINT4 length)
{

   XLAL_CHECK( output != NULL && input != NULL && length > 0 && length<=input->length , XLAL_EINVAL );

   REAL4VectorAligned *truncatedVector = NULL;
   XLAL_CHECK( (truncatedVector = XLALCreateREAL4VectorAligned(length, 16)) != NULL, XLAL_EFUNC );
   memcpy(truncatedVector->data, input->data, sizeof(REAL4)*length);
   XLAL_CHECK( sseScaleREAL4Vector(output, truncatedVector, scale) == XLAL_SUCCESS, XLAL_EFUNC );
   XLALDestroyREAL4VectorAligned(truncatedVector);

   return XLAL_SUCCESS;

}
INT4 truncated_avxScaleREAL4Vector(REAL4VectorAligned *output, REAL4VectorAligned *input, REAL4 scale, UINT4 length)
{

   XLAL_CHECK( output != NULL && input != NULL && length > 0 && length<=input->length , XLAL_EINVAL );

   REAL4VectorAligned *truncatedVector = NULL;
   XLAL_CHECK( (truncatedVector = XLALCreateREAL4VectorAligned(length, 32)) != NULL, XLAL_EFUNC );
   memcpy(truncatedVector->data, input->data, sizeof(REAL4)*length);
   XLAL_CHECK( avxScaleREAL4Vector(output, truncatedVector, scale) == XLAL_SUCCESS, XLAL_EFUNC );
   XLALDestroyREAL4VectorAligned(truncatedVector);

   return XLAL_SUCCESS;

}
INT4 truncated_sseAddScalarToREAL4Vector(REAL4VectorAligned *output, REAL4VectorAligned *input, REAL4 scalar, UINT4 length)
{

   XLAL_CHECK( output != NULL && input != NULL && length > 0 && length<=input->length , XLAL_EINVAL );

   REAL4VectorAligned *truncatedVector = NULL;
   XLAL_CHECK( (truncatedVector = XLALCreateREAL4VectorAligned(length, 16)) != NULL, XLAL_EFUNC );
   memcpy(truncatedVector->data, input->data, sizeof(REAL4)*length);
   XLAL_CHECK( sseAddScalarToREAL4Vector(output, truncatedVector, scalar) == XLAL_SUCCESS, XLAL_EFUNC );
   XLALDestroyREAL4VectorAligned(truncatedVector);

   return XLAL_SUCCESS;

}
INT4 truncated_avxAddScalarToREAL4Vector(REAL4VectorAligned *output, REAL4VectorAligned *input, REAL4 scalar, UINT4 length)
{

   XLAL_CHECK( output != NULL && input != NULL && length > 0 && length<=input->length , XLAL_EINVAL );

   REAL4VectorAligned *truncatedVector = NULL;
   XLAL_CHECK( (truncatedVector = XLALCreateREAL4VectorAligned(length, 32)) != NULL, XLAL_EFUNC );
   memcpy(truncatedVector->data, input->data, sizeof(REAL4)*length);
   XLAL_CHECK( avxAddScalarToREAL4Vector(output, truncatedVector, scalar) == XLAL_SUCCESS, XLAL_EFUNC );
   XLALDestroyREAL4VectorAligned(truncatedVector);

   return XLAL_SUCCESS;

}
INT4 truncated_sseSSVectorMultiply(REAL4VectorAligned *output, REAL4VectorAligned *input1, REAL4VectorAligned *input2, UINT4 length)
{

   XLAL_CHECK( output != NULL && input1 != NULL && input2 != NULL && length > 0 && length<=input1->length && input1->length==input2->length , XLAL_EINVAL );

   REAL4VectorAligned *truncatedVector1 = NULL, *truncatedVector2 = NULL;
   XLAL_CHECK( (truncatedVector1 = XLALCreateREAL4VectorAligned(length, 16)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (truncatedVector2 = XLALCreateREAL4VectorAligned(length, 16)) != NULL, XLAL_EFUNC );
   memcpy(truncatedVector1->data, input1->data, sizeof(REAL4)*length);
   memcpy(truncatedVector2->data, input2->data, sizeof(REAL4)*length);
   XLAL_CHECK( sseSSVectorMultiply(output, truncatedVector1, truncatedVector2) == XLAL_SUCCESS, XLAL_EFUNC );
   XLALDestroyREAL4VectorAligned(truncatedVector1);
   XLALDestroyREAL4VectorAligned(truncatedVector2);

   return XLAL_SUCCESS;

}
INT4 truncated_avxSSVectorMultiply(REAL4VectorAligned *output, REAL4VectorAligned *input1, REAL4VectorAligned *input2, UINT4 length)
{

   XLAL_CHECK( output != NULL && input1 != NULL && input2 != NULL && length > 0 && length<=input1->length && input1->length==input2->length , XLAL_EINVAL );

   REAL4VectorAligned *truncatedVector1 = NULL, *truncatedVector2 = NULL;
   XLAL_CHECK( (truncatedVector1 = XLALCreateREAL4VectorAligned(length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (truncatedVector2 = XLALCreateREAL4VectorAligned(length, 32)) != NULL, XLAL_EFUNC );
   memcpy(truncatedVector1->data, input1->data, sizeof(REAL4)*length);
   memcpy(truncatedVector2->data, input2->data, sizeof(REAL4)*length);
   XLAL_CHECK( avxSSVectorMultiply(output, truncatedVector1, truncatedVector2) == XLAL_SUCCESS, XLAL_EFUNC );
   XLALDestroyREAL4VectorAligned(truncatedVector1);
   XLALDestroyREAL4VectorAligned(truncatedVector2);

   return XLAL_SUCCESS;

}


/**
 * \brief Make an template based on FFT of sinc squared functions
 *
 * This is eq. 20 of E. Goetz and K. Riles (2011)
 * \param [out] output   Pointer to TwoSpectTemplate
 * \param [in]  input    An input candidate structure
 * \param [in]  params   Pointer to UserInput_t
 * \param [in]  sftexist Pointer to INT4Vector of existing SFTs
 * \param [in]  plan     Pointer to REAL4FFTPlan
 * \return Status value
 */
INT4 makeTemplate(TwoSpectTemplate *output, candidate input, UserInput_t *params, INT4Vector *sftexist, REAL4FFTPlan *plan)
{

   XLAL_CHECK( output != NULL && params != NULL && sftexist != NULL && plan != NULL, XLAL_EINVAL );

   REAL8 freqbin = input.fsig*params->Tsft;
   INT4 roundedbinval = (INT4)round(freqbin);
   REAL8 offset = freqbin - roundedbinval;

   INT4 ffdatabin0 = (INT4)round((params->fmin-params->dfmax)*params->Tsft) - 6;
   INT4 sigbin0 = roundedbinval - ffdatabin0;

   UINT4 numffts = (UINT4)floor(params->Tobs/(params->Tsft-params->SFToverlap)-1);
   UINT4 numfprbins = (UINT4)floorf(0.5*numffts) + 1;

   XLAL_CHECK( makeTemplate2(output, offset, input.period, input.moddepth, params->Tsft, params->SFToverlap, params->Tobs, params->minTemplateLength, (UINT4)params->vectorMath, plan) == XLAL_SUCCESS, XLAL_EFUNC );

   for (UINT4 ii=0; ii<output->pixellocations->length; ii++) {
      output->pixellocations->data[ii] += sigbin0*numfprbins;
   }

   output->f0 = input.fsig;

   return XLAL_SUCCESS;

}


/**
 * \brief Make an template based on FFT of sinc squared functions
 *
 * This is eq. 20 of E. Goetz and K. Riles (2011)
 * \param [in,out] output            Pointer to TwoSpectTemplate
 * \param [in]     offset            Amount of offset from bin centered signal (-0.5 <= offset <= 0.5)
 * \param [in]     P                 Orbital period (seconds)
 * \param [in]     deltaf            Modulation depth of the signal (Hz)
 * \param [in]     Tsft              Length of an SFT (s)
 * \param [in]     SFToverlap        SFT overlap (s)
 * \param [in]     Tobs              Observation time (s)
 * \param [in]     minTemplateLength Minimum number of pixels in a template
 * \param [in]     vectormathflag    Flag indicating to use vector math: 0 = none, 1 = SSE, 2 = AVX (must compile for vector math appropriately and CPU must have those instructions)
 * \param [in]     plan              Pointer to REAL4FFTPlan
 * \return Status value
 */
INT4 makeTemplate2(TwoSpectTemplate *output, REAL8 offset, REAL8 P, REAL8 deltaf, REAL8 Tsft, REAL8 SFToverlap, REAL8 Tobs, UINT4 minTemplateLength, UINT4 vectormathflag, REAL4FFTPlan *plan)
{
   output->f0 = offset;
   output->period = P;
   output->moddepth = deltaf;

   //Reset to zero, just in case
   memset(output->templatedata->data, 0, sizeof(REAL4)*output->templatedata->length);

   UINT4 numffts = (UINT4)floor(Tobs/(Tsft-SFToverlap)-1);
   UINT4 numfprbins = (UINT4)floorf(0.5*numffts) + 1;

   //Determine span of the template
   REAL8 binamplitude = deltaf*Tsft;
   REAL8 binmin = -binamplitude + offset;
   REAL8 binmax = binamplitude + offset;
   INT4 templatemin = (INT4)round(binmin), templatemax = (INT4)round(binmax);
   if (templatemin > binmin) templatemin--;
   if (templatemin - binmin >= -0.5) templatemin--;
   if (templatemax < binmax) templatemax++;
   if (templatemax - binmax <= 0.5) templatemax++;
   UINT4 templatespan = (UINT4)(templatemax - templatemin) + 1;

   REAL8 periodf = 1.0/P;

   REAL4Vector *psd1 = NULL;
   alignedREAL8Vector *freqbins = NULL, *bindiffs = NULL;
   XLAL_CHECK( (psd1 = XLALCreateREAL4Vector(templatespan*numffts)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (freqbins = createAlignedREAL8Vector(templatespan, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (bindiffs = createAlignedREAL8Vector(templatespan, 32)) != NULL, XLAL_EFUNC );
   memset(psd1->data, 0, sizeof(REAL4)*psd1->length);

   //Bin numbers of the frequencies
   for (UINT4 ii=0; ii<templatespan; ii++) freqbins->data[ii] = (REAL8)(templatemin + (INT4)ii);

   //Determine the signal modulation in bins with time at center of coherence time and create
   //Hann windowed PSDs
   REAL8 PSDprefact = 2.0/3.0;
   REAL4VectorAligned *t = NULL, *sigbin_sin2PiPeriodfT = NULL, *cos2PiPeriodfT = NULL;
   XLAL_CHECK( (t = XLALCreateREAL4VectorAligned(numffts, 32)) != NULL, XLAL_EFUNC  );
   XLAL_CHECK( (sigbin_sin2PiPeriodfT = XLALCreateREAL4VectorAligned(numffts, 32)) != NULL, XLAL_EFUNC  );
   XLAL_CHECK( (cos2PiPeriodfT = XLALCreateREAL4VectorAligned(numffts, 32)) != NULL, XLAL_EFUNC  );
   for (UINT4 ii=0; ii<numffts; ii++) t->data[ii] = 0.5*Tsft*(ii+1); //Assumed 50% overlapping SFTs
   if (vectormathflag==1) {
      XLAL_CHECK( sseScaleREAL4Vector(t, t, periodf) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( XLALVectorSinCos2PiREAL4(sigbin_sin2PiPeriodfT->data, cos2PiPeriodfT->data, t->data, t->length) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( sseScaleREAL4Vector(sigbin_sin2PiPeriodfT, sigbin_sin2PiPeriodfT, binamplitude) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( sseAddScalarToREAL4Vector(sigbin_sin2PiPeriodfT, sigbin_sin2PiPeriodfT, offset) == XLAL_SUCCESS, XLAL_EFUNC );
   } else if (vectormathflag==2) {
      XLAL_CHECK( avxScaleREAL4Vector(t, t, periodf) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( XLALVectorSinCos2PiREAL4(sigbin_sin2PiPeriodfT->data, cos2PiPeriodfT->data, t->data, t->length) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( avxScaleREAL4Vector(sigbin_sin2PiPeriodfT, sigbin_sin2PiPeriodfT, binamplitude) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( avxAddScalarToREAL4Vector(sigbin_sin2PiPeriodfT, sigbin_sin2PiPeriodfT, offset) == XLAL_SUCCESS, XLAL_EFUNC );
   } else {
      for (UINT4 ii=0; ii<numffts; ii++) t->data[ii] *= periodf;
      XLAL_CHECK( XLALVectorSinCos2PiREAL4(sigbin_sin2PiPeriodfT->data, cos2PiPeriodfT->data, t->data, t->length) == XLAL_SUCCESS, XLAL_EFUNC );
      for (UINT4 ii=0; ii<numffts; ii++) {
         sigbin_sin2PiPeriodfT->data[ii] *= binamplitude;
         sigbin_sin2PiPeriodfT->data[ii] += offset;
      }
   }
   XLALDestroyREAL4VectorAligned(t);
   XLALDestroyREAL4VectorAligned(cos2PiPeriodfT);
   for (UINT4 ii=0; ii<numffts; ii++) {
      REAL8 sigbin = sigbin_sin2PiPeriodfT->data[ii];
      if (vectormathflag==1) XLAL_CHECK( sseAddScalarToREAL8Vector(bindiffs, freqbins, -sigbin) == XLAL_SUCCESS, XLAL_EFUNC );
      else if (vectormathflag==2) XLAL_CHECK( avxAddScalarToREAL8Vector(bindiffs, freqbins, -sigbin) == XLAL_SUCCESS, XLAL_EFUNC );
      else for (UINT4 jj=0; jj<templatespan; jj++) bindiffs->data[jj] = freqbins->data[jj] - sigbin;
      for (UINT4 jj=0; jj<templatespan; jj++) {
         //Create PSD values organized by f0 => psd1->data[0...numffts-1], sft1 => psd1->data[numffts...2*numffts-1]
         //Restricting to +/- 1.75 bins means >99.9% of the total power is included in the template calculation
         if ( fabs(bindiffs->data[jj]) <= 1.75 ) {
            psd1->data[ii + jj*numffts] = sqsincxoverxsqminusone(bindiffs->data[jj])*PSDprefact;
            XLAL_CHECK( xlalErrno==0, XLAL_EFUNC );
         }
      } /* for jj < numfbins */
   } /* for ii < numffts */
   XLALDestroyREAL4VectorAligned(sigbin_sin2PiPeriodfT);

   //Do the second FFT
   REAL4VectorAligned *x = NULL, *psd = NULL, *windowdata = NULL;
   REAL4Window *win = NULL;
   XLAL_CHECK( (x = XLALCreateREAL4VectorAligned(numffts, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (psd = XLALCreateREAL4VectorAligned(numfprbins, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (windowdata = XLALCreateREAL4VectorAligned(x->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (win = XLALCreateHannREAL4Window(x->length)) != NULL, XLAL_EFUNC );
   memcpy(windowdata->data, win->data->data, x->length*sizeof(REAL4));
   REAL8 winFactor = 8.0/3.0, secPSDfactor = winFactor/x->length*0.5*Tsft, sum = 0.0;
   BOOLEAN doSecondFFT;
   //First loop over frequencies
   for (UINT4 ii=0; ii<templatespan; ii++) {
      //Set doSecondFFT check flag to 0. Value becomes 1 if we are to do the second FFT
      doSecondFFT = 0;

      INT4 bins2middlebin = ii + templatemin;  //Number of frequency bins away from the "central frequency bin" of the template

      //Next, loop over times and check to see if we need to do second FFT
      //Sum up the power in the row and see if it exceeds 5.0*(sinc(3.0)/(3.0^2-1))^2
      REAL4 rowpowersum = 0.0;
      for (UINT4 jj=0; jj<x->length; jj++) rowpowersum += psd1->data[ii*numffts+jj];
      if (rowpowersum > 1.187167e-34) doSecondFFT = 1;

      //If we are to do the second FFT then do it!
      if (doSecondFFT) {
         //Obtain and window the time series
         memcpy(x->data, &(psd1->data[ii*numffts]), sizeof(REAL4)*x->length);
         if (vectormathflag==1) XLAL_CHECK( sseSSVectorMultiply(x, x, windowdata) == XLAL_SUCCESS, XLAL_EFUNC );
         else if (vectormathflag==2) XLAL_CHECK( avxSSVectorMultiply(x, x, windowdata) == XLAL_SUCCESS, XLAL_EFUNC );
         else {
            XLAL_CHECK( XLALSSVectorMultiply((REAL4Vector*)x, (REAL4Vector*)x, win->data) != NULL, XLAL_EFUNC );
            XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
         }

         //Do the FFT
         XLAL_CHECK( XLALREAL4PowerSpectrum((REAL4Vector*)psd, (REAL4Vector*)x, plan) == XLAL_SUCCESS, XLAL_EFUNC );

         //Scale the data points by 1/N and window factor and (1/fs)
         //Order of vector is by second frequency then first frequency
         if (vectormathflag==1) XLAL_CHECK( sseScaleREAL4Vector(psd, psd, secPSDfactor) == XLAL_SUCCESS, XLAL_EFUNC );
         else if (vectormathflag==2) XLAL_CHECK( avxScaleREAL4Vector(psd, psd, secPSDfactor) == XLAL_SUCCESS, XLAL_EFUNC );
         else for (UINT4 jj=0; jj<psd->length; jj++) psd->data[jj] *= secPSDfactor;

         //Ignore the DC to 3rd frequency bins in sum
         for (UINT4 jj=4; jj<psd->length; jj++) {
            sum += (REAL8)psd->data[jj];     //sum up the total weight

            //Sort the weights, insertion sort technique
            if (psd->data[jj] > output->templatedata->data[output->templatedata->length-1]) insertionSort_template(output, psd->data[jj], bins2middlebin*psd->length+jj);
         } /* for jj < psd->length */
      } /* if doSecondFFT */
   } /* if ii < numfbins */

   //Normalize
   REAL4 invsum = (REAL4)(1.0/sum);
   if (vectormathflag==1) XLAL_CHECK( sseScaleREAL4Vector(output->templatedata, output->templatedata, invsum) == XLAL_SUCCESS, XLAL_EFUNC );
   else if (vectormathflag==2) XLAL_CHECK( avxScaleREAL4Vector(output->templatedata, output->templatedata, invsum) == XLAL_SUCCESS, XLAL_EFUNC );
   else for (UINT4 ii=0; ii<output->templatedata->length; ii++) if (output->templatedata->data[ii]!=0.0) output->templatedata->data[ii] *= invsum;

   //Truncate weights if they don't contribute much to the sum
   sum = 0.0;
   for (UINT4 ii=0; ii<minTemplateLength; ii++) sum += (REAL8)output->templatedata->data[ii];
   UINT4 counter = minTemplateLength;
   while (counter<output->templatedata->length && output->templatedata->data[counter]>=epsval_float((REAL4)sum)) {
      sum += (REAL8)output->templatedata->data[counter];
      counter++;
   }
   for (/* last counter val */; counter<output->templatedata->length; counter++) output->templatedata->data[counter] = 0.0;

   //Destroy stuff
   XLALDestroyREAL4Vector(psd1);
   destroyAlignedREAL8Vector(freqbins);
   destroyAlignedREAL8Vector(bindiffs);
   XLALDestroyREAL4Window(win);
   XLALDestroyREAL4VectorAligned(x);
   XLALDestroyREAL4VectorAligned(psd);
   XLALDestroyREAL4VectorAligned(windowdata);

   return XLAL_SUCCESS;

}


/**
 * Insertion sort for the template weights
 * \param [out] output        Pointer to TwoSpectTemplate
 * \param [in]  weight        Pixel weight
 * \param [in]  pixelloc      Index of the pixel in the REAL4Vector of the frequency-frequency plane
 */
void insertionSort_template(TwoSpectTemplate *output, REAL4 weight, INT4 pixelloc)
{

   INT4 insertionpoint = output->templatedata->length-1;
   INT4 numbertomove = 0;
   while (insertionpoint > 0 && weight > output->templatedata->data[insertionpoint-1]) {
      insertionpoint--;
      numbertomove++;
   }

   if (insertionpoint<(INT4)output->templatedata->length-1) {
      memmove(&(output->templatedata->data[insertionpoint+1]), &(output->templatedata->data[insertionpoint]), sizeof(REAL4)*numbertomove);
      memmove(&(output->pixellocations->data[insertionpoint+1]), &(output->pixellocations->data[insertionpoint]), sizeof(INT4)*numbertomove);
   }

   output->templatedata->data[insertionpoint] = weight;
   output->pixellocations->data[insertionpoint] = pixelloc;

} /* insertionSort_template() */


/**
 * Calculate sin(pi*x)/(pi*x)/(x^2-1)
 * \param x Value from which to compute
 */
REAL8 sincxoverxsqminusone(REAL8 x)
{
   if (fabs(x*x-1.0)<1.0e-8) return -0.5;
   if (fabs(x)<1.0e-8) return -1.0;
   REAL8 pix = LAL_PI*x;
   //return sin(pix)/(pix*(x*x-1.0));
   REAL4 sinpix = 0.0, cospix = 0.0;
   XLAL_CHECK_REAL8( XLALSinCosLUT(&sinpix, &cospix, pix) == XLAL_SUCCESS, XLAL_EFUNC );
   if (sinpix>1.0) sinpix = 1.0;
   else if (sinpix<-1.0) sinpix = -1.0;
   return sinpix/(pix*(x*x-1.0));
} /* sincxoverxsqminusone() */


/**
 * Calculate [sin(pi*x)/(pi*x)/(x^2-1)]^2
 * \param x Value from which to compute
 */
REAL8 sqsincxoverxsqminusone(REAL8 x)
{
   REAL8 val = sincxoverxsqminusone(x);
   XLAL_CHECK_REAL8( xlalErrno == 0, XLAL_EFUNC );
   return val*val;
} /* sqsincxoverxsqminusone() */



/* 
 *  LALInferenceLikelihood.c:  Bayesian Followup likelihood functions
 *
 *  Copyright (C) 2009 Ilya Mandel, Vivien Raymond, Christian Roever, Marc van der Sluys and John Veitch
 *
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

#include <lal/LALInferenceLikelihood.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/Sequence.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeFreqFFT.h>


#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/* ============ Likelihood computations: ========== */

REAL8 ZeroLogLikelihood(LALInferenceVariables UNUSED *currentParams, LALInferenceIFOData UNUSED *data, LALInferenceTemplateFunction UNUSED *template) {
  return 0.0;
}

REAL8 UndecomposedFreqDomainLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData * data, 
                              LALInferenceTemplateFunction *template)
/***************************************************************/
/* (log-) likelihood function.                                 */
/* Returns the non-normalised logarithmic likelihood.          */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Required (`currentParams') parameters are:                  */
/*   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       */
/*   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  */
/*   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        */
/*   - "distance"        (REAL8, Mpc, >0)                      */
/*   - "time"            (REAL8, GPS sec.)                     */
/***************************************************************/
{
  static int timeDomainWarning = 0;
  double Fplus, Fcross;
  double FplusScaled, FcrossScaled;
  double diffRe, diffIm, diffSquared;
  double dataReal, dataImag;
  REAL8 loglikeli;
  REAL8 plainTemplateReal, plainTemplateImag;
  REAL8 templateReal, templateImag;
  int i, lower, upper;
  LALInferenceIFOData *dataPtr;
  double ra, dec, psi, distMpc, gmst;
  double GPSdouble;
  LIGOTimeGPS GPSlal;
  double chisquared;
  double timedelay;  /* time delay b/w iterferometer & geocenter w.r.t. sky location */
  double timeshift;  /* time shift (not necessarily same as above)                   */
  double deltaT, TwoDeltaToverN, deltaF, twopit, f, re, im;
  double timeTmp;
  int different;
  LALStatus status;
  memset(&status,0,sizeof(status));
  LALInferenceVariables intrinsicParams;

  /* determine source's sky location & orientation parameters: */
  ra        = *(REAL8*) LALInferenceGetVariable(currentParams, "rightascension"); /* radian      */
  dec       = *(REAL8*) LALInferenceGetVariable(currentParams, "declination");    /* radian      */
  psi       = *(REAL8*) LALInferenceGetVariable(currentParams, "polarisation");   /* radian      */
  GPSdouble = *(REAL8*) LALInferenceGetVariable(currentParams, "time");           /* GPS seconds */
  distMpc   = *(REAL8*) LALInferenceGetVariable(currentParams, "distance");       /* Mpc         */

  /* figure out GMST: */
  //XLALINT8NSToGPS(&GPSlal, floor(1e9 * GPSdouble + 0.5));
  XLALGPSSetREAL8(&GPSlal, GPSdouble);
  //UandA.units    = MST_RAD;
  //UandA.accuracy = LALLEAPSEC_LOOSE;
  //LALGPStoGMST1(&status, &gmst, &GPSlal, &UandA);
  gmst=XLALGreenwichMeanSiderealTime(&GPSlal);
  intrinsicParams.head      = NULL;
  intrinsicParams.dimension = 0;
  LALInferenceCopyVariables(currentParams, &intrinsicParams);
  LALInferenceRemoveVariable(&intrinsicParams, "rightascension");
  LALInferenceRemoveVariable(&intrinsicParams, "declination");
  LALInferenceRemoveVariable(&intrinsicParams, "polarisation");
  LALInferenceRemoveVariable(&intrinsicParams, "time");
  LALInferenceRemoveVariable(&intrinsicParams, "distance");
  // TODO: add pointer to template function here.
  // (otherwise same parameters but different template will lead to no re-computation!!)

  chisquared = 0.0;
  /* loop over data (different interferometers): */
  dataPtr = data;

  while (dataPtr != NULL) {
    /* The parameters the Likelihood function can handle by itself   */
    /* (and which shouldn't affect the template function) are        */
    /* sky location (ra, dec), polarisation and signal arrival time. */
    /* Note that the template function shifts the waveform to so that*/
	/* t_c corresponds to the "time" parameter in                    */
	/* IFOdata->modelParams (set, e.g., from the trigger value).     */
    
    /* Reset log-likelihood */
    dataPtr->loglikelihood = 0.0;

    /* Compare parameter values with parameter values corresponding  */
    /* to currently stored template; ignore "time" variable:         */
    if (LALInferenceCheckVariable(dataPtr->modelParams, "time")) {
      timeTmp = *(REAL8 *) LALInferenceGetVariable(dataPtr->modelParams, "time");
      LALInferenceRemoveVariable(dataPtr->modelParams, "time");
    }
    else timeTmp = GPSdouble;
    different = LALInferenceCompareVariables(dataPtr->modelParams, &intrinsicParams);
    /* "different" now may also mean that "dataPtr->modelParams" */
    /* wasn't allocated yet (as in the very 1st iteration).      */

    if (different) { /* template needs to be re-computed: */
      LALInferenceCopyVariables(&intrinsicParams, dataPtr->modelParams);
      LALInferenceAddVariable(dataPtr->modelParams, "time", &timeTmp, REAL8_t,PARAM_LINEAR);
      template(dataPtr);
      if (dataPtr->modelDomain == timeDomain) {
	if (!timeDomainWarning) {
	  timeDomainWarning = 1;
	  fprintf(stderr, "WARNING: using time domain template with frequency domain likelihood (in %s, line %d)\n", __FILE__, __LINE__);
	}
        LALInferenceExecuteFT(dataPtr);
        /* note that the dataPtr->modelParams "time" element may have changed here!! */
        /* (during "template()" computation)  */
      }
    }
    else { /* no re-computation necessary. Return back "time" value, do nothing else: */
      LALInferenceAddVariable(dataPtr->modelParams, "time", &timeTmp, REAL8_t,PARAM_LINEAR);
    }

    /*-- Template is now in dataPtr->freqModelhPlus and dataPtr->freqModelhCross. --*/
    /*-- (Either freshly computed or inherited.)                            --*/

    /* determine beam pattern response (F_plus and F_cross) for given Ifo: */
    XLALComputeDetAMResponse(&Fplus, &Fcross,
                             dataPtr->detector->response,
			     ra, dec, psi, gmst);
    /* signal arrival time (relative to geocenter); */
    timedelay = XLALTimeDelayFromEarthCenter(dataPtr->detector->location,
                                             ra, dec, &GPSlal);
    /* (negative timedelay means signal arrives earlier at Ifo than at geocenter, etc.) */
    /* amount by which to time-shift template (not necessarily same as above "timedelay"): */
    timeshift =  (GPSdouble - (*(REAL8*) LALInferenceGetVariable(dataPtr->modelParams, "time"))) + timedelay;
    twopit    = LAL_TWOPI * timeshift;

    /* include distance (overall amplitude) effect in Fplus/Fcross: */
    FplusScaled  = Fplus  / distMpc;
    FcrossScaled = Fcross / distMpc;

    if (LALInferenceCheckVariable(currentParams, "crazyInjectionHLSign") &&
        *((INT4 *)LALInferenceGetVariable(currentParams, "crazyInjectionHLSign"))) {
      if (strstr(dataPtr->name, "H") || strstr(dataPtr->name, "L")) {
        FplusScaled *= -1.0;
        FcrossScaled *= -1.0;
      }
    }

    dataPtr->fPlus = FplusScaled;
    dataPtr->fCross = FcrossScaled;
    dataPtr->timeshift = timeshift;

 //FILE *testout=fopen("test_likeliLAL.txt","w");
 //fprintf(testout, "f PSD dataRe dataIm signalRe signalIm\n");
    /* determine frequency range & loop over frequency bins: */
    deltaT = dataPtr->timeData->deltaT;
    deltaF = 1.0 / (((double)dataPtr->timeData->data->length) * deltaT);
    // printf("deltaF %g, Nt %d, deltaT %g\n", deltaF, dataPtr->timeData->data->length, dataPtr->timeData->deltaT);
    lower = ceil(dataPtr->fLow / deltaF);
    upper = floor(dataPtr->fHigh / deltaF);
    TwoDeltaToverN = 2.0 * deltaT / ((double) dataPtr->timeData->data->length);
    for (i=lower; i<=upper; ++i){
      /* derive template (involving location/orientation parameters) from given plus/cross waveforms: */
      plainTemplateReal = FplusScaled * dataPtr->freqModelhPlus->data->data[i].re  
                          +  FcrossScaled * dataPtr->freqModelhCross->data->data[i].re;
      plainTemplateImag = FplusScaled * dataPtr->freqModelhPlus->data->data[i].im  
                          +  FcrossScaled * dataPtr->freqModelhCross->data->data[i].im;

      /* do time-shifting...             */
      /* (also un-do 1/deltaT scaling): */
      f = ((double) i) * deltaF;
      /* real & imag parts of  exp(-2*pi*i*f*deltaT): */
      re = cos(twopit * f);
      im = - sin(twopit * f);
      templateReal = (plainTemplateReal*re - plainTemplateImag*im) / deltaT;
      templateImag = (plainTemplateReal*im + plainTemplateImag*re) / deltaT;
      dataReal     = dataPtr->freqData->data->data[i].re / deltaT;
      dataImag     = dataPtr->freqData->data->data[i].im / deltaT;
      /* compute squared difference & 'chi-squared': */
      diffRe       = dataReal - templateReal;         // Difference in real parts...
      diffIm       = dataImag - templateImag;         // ...and imaginary parts, and...
      diffSquared  = diffRe*diffRe + diffIm*diffIm ;  // ...squared difference of the 2 complex figures.
      REAL8 temp = ((TwoDeltaToverN * diffSquared) / dataPtr->oneSidedNoisePowerSpectrum->data->data[i]);
      chisquared  += temp;
      dataPtr->loglikelihood -= temp;
 //fprintf(testout, "%e %e %e %e %e %e\n",
 //        f, dataPtr->oneSidedNoisePowerSpectrum->data->data[i], 
 //        dataPtr->freqData->data->data[i].re, dataPtr->freqData->data->data[i].im,
 //        templateReal, templateImag);
    }
    dataPtr = dataPtr->next;
 //fclose(testout);
  }
  loglikeli = -1.0 * chisquared; // note (again): the log-likelihood is unnormalised!
  LALInferenceDestroyVariables(&intrinsicParams);
  return(loglikeli);
}


REAL8 FreqDomainStudentTLogLikelihood(LALVariables *currentParams, LALIFOData *data, 
                                      LALTemplateFunction *template,
                                      LALVariables *df)
/***************************************************************/
/* Student-t (log-) likelihood function                        */
/* as described in Roever/Meyer/Christensen (2011):            */
/*   "Modelling coloured residual noise                        */
/*   in gravitational-wave signal processing."                 */
/*   Classical and Quantum Gravity, 28(1):015010.              */
/*   http://dx.doi.org/10.1088/0264-9381/28/1/015010           */
/*   http://arxiv.org/abs/0804.3853                            */
/* Returns the non-normalised logarithmic likelihood.          */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Required (`currentParams') parameters are:                  */
/*   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       */
/*   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  */
/*   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        */
/*   - "distance"        (REAL8, Mpc, > 0)                     */
/*   - "time"            (REAL8, GPS sec.)                     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* This function is essentially the same as the                */
/* "UndecomposedFreqDomainLogLikelihood()" function.           */
/* The additional parameter to be supplied is the (REAL8)      */
/* degrees-of-freedom parameter (nu) for each Ifo.             */
/* The additional "df" argument gives the corresponding        */
/* d.f. parameter for each element of the "*data" list.        */
/* The names of "df" must match the "->name" slot of           */
/* the elements of "data".                                     */
/*                                                             */
/* (TODO: allow for d.f. parameter to vary with frequency,     */
/*        i.e., to be a set of vectors corresponding to        */
/*        frequencies)                                         */
/***************************************************************/
{
  static int timeDomainWarning = 0;
  double Fplus, Fcross;
  double FplusScaled, FcrossScaled;
  double diffRe, diffIm, diffSquared;
  double dataReal, dataImag;
  REAL8 loglikeli;
  REAL8 plainTemplateReal, plainTemplateImag;
  REAL8 templateReal, templateImag;
  int i, lower, upper;
  LALIFOData *dataPtr;
  double ra, dec, psi, distMpc, gmst;
  double GPSdouble;
  LIGOTimeGPS GPSlal;
  double chisquared;
  double timedelay;  /* time delay b/w iterferometer & geocenter w.r.t. sky location */
  double timeshift;  /* time shift (not necessarily same as above)                   */
  double deltaT, FourDeltaToverN, deltaF, twopit, f, re, im, singleFreqBinTerm;
  double degreesOfFreedom, nu;
  double timeTmp;
  int different;
  LALStatus status;
  memset(&status,0,sizeof(status));
  LALVariables intrinsicParams;

  /* determine source's sky location & orientation parameters: */
  ra        = *(REAL8*) getVariable(currentParams, "rightascension"); /* radian      */
  dec       = *(REAL8*) getVariable(currentParams, "declination");    /* radian      */
  psi       = *(REAL8*) getVariable(currentParams, "polarisation");   /* radian      */
  GPSdouble = *(REAL8*) getVariable(currentParams, "time");           /* GPS seconds */
  distMpc   = *(REAL8*) getVariable(currentParams, "distance");       /* Mpc         */

  /* figure out GMST: */
  /* XLALINT8NSToGPS(&GPSlal, floor(1e9 * GPSdouble + 0.5)); */
  XLALGPSSetREAL8(&GPSlal, GPSdouble);
  /* UandA.units    = MST_RAD;                       */
  /* UandA.accuracy = LALLEAPSEC_LOOSE;              */
  /* LALGPStoGMST1(&status, &gmst, &GPSlal, &UandA); */
  gmst=XLALGreenwichMeanSiderealTime(&GPSlal);
  intrinsicParams.head      = NULL;
  intrinsicParams.dimension = 0;
  copyVariables(currentParams, &intrinsicParams);
  removeVariable(&intrinsicParams, "rightascension");
  removeVariable(&intrinsicParams, "declination");
  removeVariable(&intrinsicParams, "polarisation");
  removeVariable(&intrinsicParams, "time");
  removeVariable(&intrinsicParams, "distance");
  /*  TODO: add pointer to template function here.                                         */
  /*  (otherwise same parameters but different template will lead to no re-computation!!)  */

  chisquared = 0.0;
  /* loop over data (different interferometers): */
  dataPtr = data;

  while (dataPtr != NULL) {
    /* The parameters the Likelihood function can handle by itself    */
    /* (and which shouldn't affect the template function) are         */
    /* sky location (ra, dec), polarisation and signal arrival time.  */
    /* Note that the template function shifts the waveform to so that */
    /* t_c corresponds to the "time" parameter in                     */
    /* IFOdata->modelParams (set, e.g., from the trigger value).      */
    
    /* Reset log-likelihood */
    dataPtr->loglikelihood = 0.0;

    /* Compare parameter values with parameter values corresponding */
    /* to currently stored template; ignore "time" variable:        */
    if (checkVariable(dataPtr->modelParams, "time")) {
      timeTmp = *(REAL8 *) getVariable(dataPtr->modelParams, "time");
      removeVariable(dataPtr->modelParams, "time");
    }
    else timeTmp = GPSdouble;
    different = compareVariables(dataPtr->modelParams, &intrinsicParams);
    /* "different" now may also mean that "dataPtr->modelParams" */
    /* wasn't allocated yet (as in the very 1st iteration).      */

    if (different) { /* template needs to be re-computed: */
      copyVariables(&intrinsicParams, dataPtr->modelParams);
      addVariable(dataPtr->modelParams, "time", &timeTmp, REAL8_t,PARAM_LINEAR);
      template(dataPtr);
      if (dataPtr->modelDomain == timeDomain) {
	if (!timeDomainWarning) {
	  timeDomainWarning = 1;
	  fprintf(stderr, "WARNING: using time domain template with frequency domain likelihood (in %s, line %d)\n", __FILE__, __LINE__);
	}
        executeFT(dataPtr);
        /* note that the dataPtr->modelParams "time" element may have changed here!! */
        /* (during "template()" computation)  */
      }
    }
    else { /* no re-computation necessary. Return back "time" value, do nothing else: */
      addVariable(dataPtr->modelParams, "time", &timeTmp, REAL8_t,PARAM_LINEAR);
    }

    /*-- Template is now in dataPtr->freqModelhPlus and dataPtr->freqModelhCross. --*/
    /*-- (Either freshly computed or inherited.)                                  --*/

    /* determine beam pattern response (F_plus and F_cross) for given Ifo: */
    XLALComputeDetAMResponse(&Fplus, &Fcross,
                             dataPtr->detector->response,
			     ra, dec, psi, gmst);
    /* signal arrival time (relative to geocenter); */
    timedelay = XLALTimeDelayFromEarthCenter(dataPtr->detector->location,
                                             ra, dec, &GPSlal);
    /* (negative timedelay means signal arrives earlier at Ifo than at geocenter, etc.)    */
    /* amount by which to time-shift template (not necessarily same as above "timedelay"): */
    timeshift =  (GPSdouble - (*(REAL8*) getVariable(dataPtr->modelParams, "time"))) + timedelay;
    twopit    = LAL_TWOPI * timeshift;

    /* include distance (overall amplitude) effect in Fplus/Fcross: */
    FplusScaled  = Fplus  / distMpc;
    FcrossScaled = Fcross / distMpc;

    if (checkVariable(currentParams, "crazyInjectionHLSign") &&
        *((INT4 *)getVariable(currentParams, "crazyInjectionHLSign"))) {
      if (strstr(dataPtr->name, "H") || strstr(dataPtr->name, "L")) {
        FplusScaled *= -1.0;
        FcrossScaled *= -1.0;
      }
    }

    dataPtr->fPlus = FplusScaled;
    dataPtr->fCross = FcrossScaled;
    dataPtr->timeshift = timeshift;

    /* extract the element from the "df" vector that carries the current Ifo's name: */
    degreesOfFreedom = *(REAL8*) getVariable(df, dataPtr->name);
    if (!(degreesOfFreedom>0)) die(" ERROR in StudentTLogLikelihood(): degrees-of-freedom parameter (\"df\") must be positive.\n");

    /* determine frequency range & loop over frequency bins: */
    deltaT = dataPtr->timeData->deltaT;
    deltaF = 1.0 / (((double)dataPtr->timeData->data->length) * deltaT);
    lower = ceil(dataPtr->fLow / deltaF);
    upper = floor(dataPtr->fHigh / deltaF);
    FourDeltaToverN = 4.0 * deltaT / ((double) dataPtr->timeData->data->length);
    for (i=lower; i<=upper; ++i){
      /* degrees-of-freedom parameter (nu_j) for this particular frequency bin: */
      nu = degreesOfFreedom;
      /* (for now constant across frequencies)                                  */
      /* derive template (involving location/orientation parameters) from given plus/cross waveforms: */
      plainTemplateReal = FplusScaled * dataPtr->freqModelhPlus->data->data[i].re  
                          +  FcrossScaled * dataPtr->freqModelhCross->data->data[i].re;
      plainTemplateImag = FplusScaled * dataPtr->freqModelhPlus->data->data[i].im  
                          +  FcrossScaled * dataPtr->freqModelhCross->data->data[i].im;

      /* do time-shifting...            */
      /* (also un-do 1/deltaT scaling): */
      f = ((double) i) * deltaF;
      /* real & imag parts of  exp(-2*pi*i*f*deltaT): */
      re = cos(twopit * f);
      im = - sin(twopit * f);
      templateReal = (plainTemplateReal*re - plainTemplateImag*im) / deltaT;
      templateImag = (plainTemplateReal*im + plainTemplateImag*re) / deltaT;
      dataReal     = dataPtr->freqData->data->data[i].re / deltaT;
      dataImag     = dataPtr->freqData->data->data[i].im / deltaT;
      /* compute squared difference & 'chi-squared': */
      diffRe       = dataReal - templateReal;         /* Difference in real parts...                     */
      diffIm       = dataImag - templateImag;         /* ...and imaginary parts, and...                  */
      diffSquared  = diffRe*diffRe + diffIm*diffIm ;  /* ...squared difference of the 2 complex figures. */
      singleFreqBinTerm = ((nu+2.0)/2.0) * log((FourDeltaToverN * diffSquared) / (nu * dataPtr->oneSidedNoisePowerSpectrum->data->data[i]));
      chisquared  += singleFreqBinTerm;   /* (This is a sum-of-squares, or chi^2, term in the Gaussian case, not so much in the Student-t case...)  */
      dataPtr->loglikelihood -= singleFreqBinTerm;
    }
    dataPtr = dataPtr->next;
  }
  loglikeli = -1.0 * chisquared; /* note (again): the log-likelihood is unnormalised! */
  destroyVariables(&intrinsicParams);  
  return(loglikeli);
}



REAL8 FreqDomainLogLikelihood(LALVariables *currentParams, LALIFOData * data, 
                              LALTemplateFunction *template)
/***************************************************************/
/* (log-) likelihood function.                                 */
/* Returns the non-normalised logarithmic likelihood.          */
/* Slightly slower but cleaner than							   */
/* UndecomposedFreqDomainLogLikelihood().          `		   */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Required (`currentParams') parameters are:                  */
/*   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       */
/*   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  */
/*   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        */
/*   - "distance"        (REAL8, Mpc, >0)                      */
/*   - "time"            (REAL8, GPS sec.)                     */
/***************************************************************/
{
  REAL8 loglikeli, totalChiSquared=0.0;
  LALInferenceIFOData *ifoPtr=data;
  COMPLEX16Vector *freqModelResponse=NULL;

  /* loop over data (different interferometers): */
  while (ifoPtr != NULL) {
    ifoPtr->loglikelihood = 0.0;

	if(freqModelResponse==NULL)
		freqModelResponse= XLALCreateCOMPLEX16Vector(ifoPtr->freqData->data->length);
	else
		freqModelResponse= XLALResizeCOMPLEX16Vector(freqModelResponse, ifoPtr->freqData->data->length);
	/*compute the response*/
	ComputeFreqDomainResponse(currentParams, ifoPtr, template, freqModelResponse);
	/*if(residual==NULL)
		residual=XLALCreateCOMPLEX16Vector(ifoPtr->freqData->data->length);
	else
		residual=XLALResizeCOMPLEX16Vector(residual, ifoPtr->freqData->data->length);
	
	COMPLEX16VectorSubtract(residual, ifoPtr->freqData->data, freqModelResponse);
	totalChiSquared+=ComputeFrequencyDomainOverlap(ifoPtr, residual, residual); 
	*/
        REAL8 temp = ComputeFrequencyDomainOverlap(ifoPtr, ifoPtr->freqData->data, ifoPtr->freqData->data)
          -2.0*ComputeFrequencyDomainOverlap(ifoPtr, ifoPtr->freqData->data, freqModelResponse)
          +ComputeFrequencyDomainOverlap(ifoPtr, freqModelResponse, freqModelResponse);
	totalChiSquared+=temp;
        ifoPtr->loglikelihood -= 0.5*temp;

    ifoPtr = ifoPtr->next;
  }
  loglikeli = -0.5 * totalChiSquared; // note (again): the log-likelihood is unnormalised!
  XLALDestroyCOMPLEX16Vector(freqModelResponse);
  return(loglikeli);
}

REAL8 ChiSquareTest(LALInferenceVariables *currentParams, LALInferenceIFOData * data, LALInferenceTemplateFunction *template)
/***************************************************************/
/* Chi-Square function.                                        */
/* Returns the chi square of a template:                       */
/* chisq= p * sum_i (dx_i)^2, with dx_i  =  <s,h>_i  - <s,h>/p */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Required (`currentParams') parameters are:                  */
/*   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       */
/*   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  */
/*   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        */
/*   - "distance"        (REAL8, Mpc, >0)                      */
/*   - "time"            (REAL8, GPS sec.)                     */
/***************************************************************/
{
  REAL8 ChiSquared=0.0, dxp, xp, x, norm, binPower, nextBin;
  REAL8 lowerF, upperF, deltaT, deltaF;
  REAL8 *segnorm;
  INT4  i, chisqPt, imax, kmin, kmax, numBins=0;
  INT4  *chisqBin;
  LALInferenceIFOData *ifoPtr=data;
  COMPLEX16Vector *freqModelResponse=NULL;
  //COMPLEX16Vector *segmentFreqModelResponse-NULL;

  /* Allocate memory for local pointers */
  segnorm=malloc(sizeof(REAL8) * ifoPtr->freqData->data->length);
  chisqBin=malloc(sizeof(INT4) * (numBins + 1));

  /* loop over data (different interferometers): */
  while (ifoPtr != NULL) {
    if(freqModelResponse==NULL)
      freqModelResponse= XLALCreateCOMPLEX16Vector(ifoPtr->freqData->data->length);
    else
      freqModelResponse= XLALResizeCOMPLEX16Vector(freqModelResponse, ifoPtr->freqData->data->length);
    /*compute the response*/
    ComputeFreqDomainResponse(currentParams, ifoPtr, template, freqModelResponse);

    deltaT = ifoPtr->timeData->deltaT;
    deltaF = 1.0 / (((REAL8)ifoPtr->timeData->data->length) * deltaT);
    
    /* Store values of fLow and fHigh to use later */
    lowerF = ifoPtr->fLow;
    upperF = ifoPtr->fHigh;
   
    /* Generate bin boundaries */
    numBins = *(INT4*) LALInferenceGetVariable(currentParams, "numbins");
    kmin = ceil(ifoPtr->fLow / deltaF);
    kmax = floor(ifoPtr->fHigh / deltaF);
    imax = kmax > (INT4) ifoPtr->freqData->data->length-1 ? (INT4) ifoPtr->freqData->data->length-1 : kmax;
    
    memset(segnorm,0,sizeof(REAL8) * ifoPtr->freqData->data->length);
    norm = 0.0;
    
    for (i=1; i < imax; ++i){  	  	  
      norm += ((4.0 * deltaF * (freqModelResponse->data[i].re*freqModelResponse->data[i].re
              +freqModelResponse->data[i].im*freqModelResponse->data[i].im)) 
              / ifoPtr->oneSidedNoisePowerSpectrum->data->data[i]);
      segnorm[i] = norm;
    }


    memset(chisqBin,0,sizeof(INT4) * (numBins +1));

    binPower = norm / (REAL8) numBins;
    nextBin   = binPower;
    chisqPt   = 0;
    chisqBin[chisqPt++] = 0;

    for ( i = 1; i < imax; ++i )
    {
      if ( segnorm[i] >= nextBin )
      {
        chisqBin[chisqPt++] = i;
        nextBin += binPower;
        if ( chisqPt == numBins ) break;
      }
    }
    chisqBin[16]=imax;
    /* check that we have sucessfully allocated all the bins */
    if ( i == (INT4) ifoPtr->freqData->data->length && chisqPt != numBins )
    {
      /* if we have reaced the end of the template power vec and not
       * */
      /* allocated all the bin boundaries then there is a problem
       * */
      fprintf(stderr,"Error constructing frequency bins\n"); 
    }

    /* the last bin boundary is at can be at Nyquist since   */
    /* qtilde is zero above the ISCO of the current template */
    // chisqBin[numBins] = ifoPtr->freqData->data->length;

    /* end */
    
    x = ComputeFrequencyDomainOverlap(ifoPtr, ifoPtr->freqData->data, freqModelResponse)/(sqrt(norm));
    
    ChiSquared=0.0;
    
    for (i=0; i < numBins; ++i){
      
      ifoPtr->fLow = chisqBin[i] * deltaF;
      ifoPtr->fHigh = chisqBin[i+1] * deltaF;
      
      xp = ComputeFrequencyDomainOverlap(ifoPtr, ifoPtr->freqData->data, freqModelResponse)/(sqrt(norm));
      dxp = ((REAL8) numBins) * xp - x;
      ChiSquared += (dxp * dxp);
      
    }
    ChiSquared = ChiSquared / (REAL8) numBins;
    ifoPtr->fLow = lowerF;
    ifoPtr->fHigh = upperF;
    
    printf("Chi-Square for %s\t=\t%f\n",ifoPtr->detector->frDetector.name,ChiSquared);
    
    ifoPtr = ifoPtr->next;
  }
  free(chisqBin);
  free(segnorm);
  return(ChiSquared);
}

REAL8 TimeDomainLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData * data, 
                              LALInferenceTemplateFunction *template)
/***************************************************************/
/* (log-) likelihood function.                                 */
/* Returns the non-normalised logarithmic likelihood.          */
/* Slightly slower but cleaner than							   */
/* UndecomposedFreqDomainLogLikelihood().          `		   */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Required (`currentParams') parameters are:                  */
/*   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       */
/*   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  */
/*   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        */
/*   - "distance"        (REAL8, Mpc, >0)                      */
/*   - "time"            (REAL8, GPS sec.)                     */
/***************************************************************/
{
  REAL8 loglikeli, totalChiSquared=0.0;
  LALInferenceIFOData *ifoPtr=data;
  REAL8TimeSeries *timeModelResponse=NULL;

  /* loop over data (different interferometers): */
  while (ifoPtr != NULL) {
    ifoPtr->loglikelihood = 0.0;

    if(timeModelResponse==NULL) {
      timeModelResponse = 
	XLALCreateREAL8TimeSeries("time detector response", &(ifoPtr->timeData->epoch), 
				  0.0, ifoPtr->timeData->deltaT,
				  &lalDimensionlessUnit,
				  ifoPtr->timeData->data->length);
    } else if (timeModelResponse->data->length != ifoPtr->timeData->data->length) {
      /* Cannot resize *up* a time series, so just dealloc and reallocate it. */
      XLALDestroyREAL8TimeSeries(timeModelResponse);
      timeModelResponse =                   
	XLALCreateREAL8TimeSeries("time detector response", &(ifoPtr->timeData->epoch), 
				  0.0, ifoPtr->timeData->deltaT,
				  &lalDimensionlessUnit,
				  ifoPtr->timeData->data->length);
    }
     
    /*compute the response*/
    ComputeTimeDomainResponse(currentParams, ifoPtr, template, timeModelResponse);
    REAL8 temp = (WhitenedTimeDomainOverlap(ifoPtr->whiteTimeData, ifoPtr->windowedTimeData)
                  -2.0*WhitenedTimeDomainOverlap(ifoPtr->whiteTimeData, timeModelResponse)
                  +timeDomainOverlap(ifoPtr->timeDomainNoiseWeights, timeModelResponse, timeModelResponse));
    totalChiSquared+=temp;
    ifoPtr->loglikelihood -= 0.5*temp;
    
    ifoPtr = ifoPtr->next;
  }
  loglikeli = -0.5*totalChiSquared; 
  XLALDestroyREAL8TimeSeries(timeModelResponse);
  return(loglikeli);
}

void ComputeFreqDomainResponse(LALInferenceVariables *currentParams, LALInferenceIFOData * dataPtr, 
                              LALInferenceTemplateFunction *template, COMPLEX16Vector *freqWaveform)
/***************************************************************/
/* Frequency-domain single-IFO response computation.           */
/* Computes response for a given template.                     */
/* Will re-compute template only if necessary                  */
/* (i.e., if previous, as stored in data->freqModelhCross,     */
/* was based on different parameters or template function).    */
/* Carries out timeshifting for a given detector               */
/* and projection onto this detector.                          */
/* Result stored in freqResponse, assumed to be correctly      */
/* initialized												   */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Required (`currentParams') parameters are:                  */
/*   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       */
/*   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  */
/*   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        */
/*   - "distance"        (REAL8, Mpc, >0)                      */
/*   - "time"            (REAL8, GPS sec.)                     */
/***************************************************************/							  
{
  static int timeDomainWarning = 0;

	double ra, dec, psi, distMpc, gmst;
	
	double GPSdouble;
	double timeTmp;
	LIGOTimeGPS GPSlal;
	double timedelay;  /* time delay b/w iterferometer & geocenter w.r.t. sky location */
	double timeshift;  /* time shift (not necessarily same as above)                   */
	double deltaT, deltaF, twopit, f, re, im;

	int different;
	LALInferenceVariables intrinsicParams;
	LALStatus status;
	memset(&status,0,sizeof(status));

	double Fplus, Fcross;
	double FplusScaled, FcrossScaled;
	REAL8 plainTemplateReal, plainTemplateImag;
	UINT4 i;
	REAL8 mc;
	
	/* Fill in derived parameters if necessary */
	if(LALInferenceCheckVariable(currentParams,"logdistance")){
		distMpc=exp(*(REAL8 *) LALInferenceGetVariable(currentParams,"logdistance"));
		LALInferenceAddVariable(currentParams,"distance",&distMpc,REAL8_t,PARAM_OUTPUT);
	}

	if(LALInferenceCheckVariable(currentParams,"logmc")){
		mc=exp(*(REAL8 *)LALInferenceGetVariable(currentParams,"logmc"));
		LALInferenceAddVariable(currentParams,"chirpmass",&mc,REAL8_t,PARAM_OUTPUT);
	}
		
	
	/* determine source's sky location & orientation parameters: */
	ra        = *(REAL8*) LALInferenceGetVariable(currentParams, "rightascension"); /* radian      */
	dec       = *(REAL8*) LALInferenceGetVariable(currentParams, "declination");    /* radian      */
	psi       = *(REAL8*) LALInferenceGetVariable(currentParams, "polarisation");   /* radian      */
	GPSdouble = *(REAL8*) LALInferenceGetVariable(currentParams, "time");           /* GPS seconds */
	distMpc   = *(REAL8*) LALInferenceGetVariable(currentParams, "distance");       /* Mpc         */
		
	/* figure out GMST: */
	//XLALINT8NSToGPS(&GPSlal, floor(1e9 * GPSdouble + 0.5));
	XLALGPSSetREAL8(&GPSlal, GPSdouble);
	//UandA.units    = MST_RAD;
	//UandA.accuracy = LALLEAPSEC_LOOSE;
	//LALGPStoGMST1(&status, &gmst, &GPSlal, &UandA);
	gmst=XLALGreenwichMeanSiderealTime(&GPSlal);
	intrinsicParams.head      = NULL;
	intrinsicParams.dimension = 0;
	LALInferenceCopyVariables(currentParams, &intrinsicParams);
	LALInferenceRemoveVariable(&intrinsicParams, "rightascension");
	LALInferenceRemoveVariable(&intrinsicParams, "declination");
	LALInferenceRemoveVariable(&intrinsicParams, "polarisation");
	LALInferenceRemoveVariable(&intrinsicParams, "time");
	LALInferenceRemoveVariable(&intrinsicParams, "distance");
	// TODO: add pointer to template function here.
	// (otherwise same parameters but different template will lead to no re-computation!!)
      
	/* The parameters the response function can handle by itself     */
    /* (and which shouldn't affect the template function) are        */
    /* sky location (ra, dec), polarisation and signal arrival time. */
    /* Note that the template function shifts the waveform to so that*/
	/* t_c corresponds to the "time" parameter in                    */
	/* IFOdata->modelParams (set, e.g., from the trigger value).     */
    
    /* Compare parameter values with parameter values corresponding  */
    /* to currently stored template; ignore "time" variable:         */
    if (LALInferenceCheckVariable(dataPtr->modelParams, "time")) {
      timeTmp = *(REAL8 *) LALInferenceGetVariable(dataPtr->modelParams, "time");
      LALInferenceRemoveVariable(dataPtr->modelParams, "time");
    }
    else timeTmp = GPSdouble;
    different = LALInferenceCompareVariables(dataPtr->modelParams, &intrinsicParams);
    /* "different" now may also mean that "dataPtr->modelParams" */
    /* wasn't allocated yet (as in the very 1st iteration).      */

    if (different) { /* template needs to be re-computed: */
      LALInferenceCopyVariables(&intrinsicParams, dataPtr->modelParams);
      LALInferenceAddVariable(dataPtr->modelParams, "time", &timeTmp, REAL8_t,PARAM_LINEAR);
      template(dataPtr);
      if (dataPtr->modelDomain == timeDomain) {
	if (!timeDomainWarning) {
	  timeDomainWarning = 1;
	  fprintf(stderr, "WARNING: using time domain template with frequency domain likelihood (in %s, line %d)\n", __FILE__, __LINE__);
	}
	LALInferenceExecuteFT(dataPtr);
      /* note that the dataPtr->modelParams "time" element may have changed here!! */
      /* (during "template()" computation)                                      */
      }
    }
    else { /* no re-computation necessary. Return back "time" value, do nothing else: */
      LALInferenceAddVariable(dataPtr->modelParams, "time", &timeTmp, REAL8_t,PARAM_LINEAR);
    }

    /*-- Template is now in dataPtr->freqModelhPlus and dataPtr->freqModelhCross. --*/
    /*-- (Either freshly computed or inherited.)                            --*/

    /* determine beam pattern response (F_plus and F_cross) for given Ifo: */
    XLALComputeDetAMResponse(&Fplus, &Fcross, dataPtr->detector->response,
			     ra, dec, psi, gmst);
		 
    /* signal arrival time (relative to geocenter); */
    timedelay = XLALTimeDelayFromEarthCenter(dataPtr->detector->location,
                                             ra, dec, &GPSlal);
    /* (negative timedelay means signal arrives earlier at Ifo than at geocenter, etc.) */

    /* amount by which to time-shift template (not necessarily same as above "timedelay"): */
    timeshift =  (GPSdouble - (*(REAL8*) LALInferenceGetVariable(dataPtr->modelParams, "time"))) + timedelay;
    twopit    = LAL_TWOPI * timeshift;

    /* include distance (overall amplitude) effect in Fplus/Fcross: */
    FplusScaled  = Fplus  / distMpc;
    FcrossScaled = Fcross / distMpc;

    
    if (LALInferenceCheckVariable(currentParams, "crazyInjectionHLSign") &&
        *((INT4 *)LALInferenceGetVariable(currentParams, "crazyInjectionHLSign"))) {
      if (strstr(dataPtr->name, "H") || strstr(dataPtr->name, "L")) {
        FplusScaled *= -1.0;
        FcrossScaled *= -1.0;
      }
    }

	if(freqWaveform->length!=dataPtr->freqModelhPlus->data->length){
		printf("fW%d data%d\n", freqWaveform->length, dataPtr->freqModelhPlus->data->length);
		printf("Error!  Frequency data vector must be same length as original data!\n");
		exit(1);
	}
	
	deltaT = dataPtr->timeData->deltaT;
    deltaF = 1.0 / (((double)dataPtr->timeData->data->length) * deltaT);

#ifdef DEBUG
FILE* file=fopen("TempSignal.dat", "w");	
#endif
	for(i=0; i<freqWaveform->length; i++){
		/* derive template (involving location/orientation parameters) from given plus/cross waveforms: */
		plainTemplateReal = FplusScaled * dataPtr->freqModelhPlus->data->data[i].re  
                          +  FcrossScaled * dataPtr->freqModelhCross->data->data[i].re;
		plainTemplateImag = FplusScaled * dataPtr->freqModelhPlus->data->data[i].im  
                          +  FcrossScaled * dataPtr->freqModelhCross->data->data[i].im;

		/* do time-shifting...             */
		/* (also un-do 1/deltaT scaling): */
		f = ((double) i) * deltaF;
		/* real & imag parts of  exp(-2*pi*i*f*deltaT): */
		re = cos(twopit * f);
		im = - sin(twopit * f);

		freqWaveform->data[i].re= (plainTemplateReal*re - plainTemplateImag*im);
		freqWaveform->data[i].im= (plainTemplateReal*im + plainTemplateImag*re);		
#ifdef DEBUG
		fprintf(file, "%lg %lg \t %lg\n", f, freqWaveform->data[i].re, freqWaveform->data[i].im);
#endif
	}
#ifdef DEBUG
fclose(file);
#endif
	LALInferenceDestroyVariables(&intrinsicParams);
}

void ComputeTimeDomainResponse(LALInferenceVariables *currentParams, LALInferenceIFOData * dataPtr, 
                               LALInferenceTemplateFunction *template, REAL8TimeSeries *timeWaveform)
/***************************************************************/
/* Based on ComputeFreqDomainResponse above.                   */
/* Time-domain single-IFO response computation.                */
/* Computes response for a given template.                     */
/* Will re-compute template only if necessary                  */
/* (i.e., if previous, as stored in data->timeModelhCross,     */
/* was based on different parameters or template function).    */
/* Carries out timeshifting for a given detector               */
/* and projection onto this detector.                          */
/* Result stored in timeResponse, assumed to be correctly      */
/* initialized												   */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Required (`currentParams') parameters are:                  */
/*   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       */
/*   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  */
/*   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        */
/*   - "distance"        (REAL8, Mpc, >0)                      */
/*   - "time"            (REAL8, GPS sec.)                     */
/***************************************************************/							  
{
  static int freqDomainWarning = 0;
	double ra, dec, psi, distMpc, gmst;
	
	double GPSdouble;
	double timeTmp;
	LIGOTimeGPS GPSlal;
	double timedelay;  /* time delay b/w iterferometer & geocenter w.r.t. sky location */
	double timeshift;  /* time shift (not necessarily same as above)                   */

	int different;
	LALInferenceVariables intrinsicParams;
	LALStatus status;
	memset(&status,0,sizeof(status));

	double Fplus, Fcross;
	double FplusScaled, FcrossScaled;
	UINT4 i;
	REAL8 mc;

	/* Fill in derived parameters if necessary */
	if(LALInferenceCheckVariable(currentParams,"logdistance")){
		distMpc=exp(*(REAL8 *) LALInferenceGetVariable(currentParams,"logdistance"));
		LALInferenceAddVariable(currentParams,"distance",&distMpc,REAL8_t,PARAM_OUTPUT);
	}

	if(LALInferenceCheckVariable(currentParams,"logmc")){
		mc=exp(*(REAL8 *)LALInferenceGetVariable(currentParams,"logmc"));
		LALInferenceAddVariable(currentParams,"chirpmass",&mc,REAL8_t,PARAM_OUTPUT);
	}
		
	
	/* determine source's sky location & orientation parameters: */
	ra        = *(REAL8*) LALInferenceGetVariable(currentParams, "rightascension"); /* radian      */
	dec       = *(REAL8*) LALInferenceGetVariable(currentParams, "declination");    /* radian      */
	psi       = *(REAL8*) LALInferenceGetVariable(currentParams, "polarisation");   /* radian      */
	GPSdouble = *(REAL8*) LALInferenceGetVariable(currentParams, "time");           /* GPS seconds */
	distMpc   = *(REAL8*) LALInferenceGetVariable(currentParams, "distance");       /* Mpc         */
		
	/* figure out GMST: */
	XLALGPSSetREAL8(&GPSlal, GPSdouble);
	//UandA.units    = MST_RAD;
	//UandA.accuracy = LALLEAPSEC_LOOSE;
	//LALGPStoGMST1(&status, &gmst, &GPSlal, &UandA);
	gmst=XLALGreenwichMeanSiderealTime(&GPSlal);
	intrinsicParams.head      = NULL;
	intrinsicParams.dimension = 0;
	LALInferenceCopyVariables(currentParams, &intrinsicParams);
	LALInferenceRemoveVariable(&intrinsicParams, "rightascension");
	LALInferenceRemoveVariable(&intrinsicParams, "declination");
	LALInferenceRemoveVariable(&intrinsicParams, "polarisation");
	LALInferenceRemoveVariable(&intrinsicParams, "time");
	LALInferenceRemoveVariable(&intrinsicParams, "distance");
	// TODO: add pointer to template function here.
	// (otherwise same parameters but different template will lead to no re-computation!!)
      
	/* The parameters the response function can handle by itself     */
    /* (and which shouldn't affect the template function) are        */
    /* sky location (ra, dec), polarisation and signal arrival time. */
    /* Note that the template function shifts the waveform to so that*/
	/* t_c corresponds to the "time" parameter in                    */
	/* IFOdata->modelParams (set, e.g., from the trigger value).     */
    
    /* Compare parameter values with parameter values corresponding  */
    /* to currently stored template; ignore "time" variable:         */
    if (LALInferenceCheckVariable(dataPtr->modelParams, "time")) {
      timeTmp = *(REAL8 *) LALInferenceGetVariable(dataPtr->modelParams, "time");
      LALInferenceRemoveVariable(dataPtr->modelParams, "time");
    }
    else timeTmp = GPSdouble;
    different = LALInferenceCompareVariables(dataPtr->modelParams, &intrinsicParams);
    /* "different" now may also mean that "dataPtr->modelParams" */
    /* wasn't allocated yet (as in the very 1st iteration).      */

    if (different) { /* template needs to be re-computed: */
      LALInferenceCopyVariables(&intrinsicParams, dataPtr->modelParams);
      LALInferenceAddVariable(dataPtr->modelParams, "time", &timeTmp, REAL8_t,PARAM_LINEAR);
      template(dataPtr);
      if (dataPtr->modelDomain == frequencyDomain) {
	if (!freqDomainWarning) {
	  freqDomainWarning = 1;
	  fprintf(stderr, "WARNING: frequency domain template used with time domain calculation (in %s, line %d)\n", __FILE__, __LINE__);
	}
        LALInferenceExecuteInvFT(dataPtr);
	/* note that the dataPtr->modelParams "time" element may have changed here!! */
	/* (during "template()" computation)  */
      }
    }
    else { /* no re-computation necessary. Return back "time" value, do nothing else: */
      LALInferenceAddVariable(dataPtr->modelParams, "time", &timeTmp, REAL8_t,PARAM_LINEAR);
    }

    /*-- Template is now in dataPtr->timeModelhPlus and dataPtr->timeModelhCross. --*/
    /*-- (Either freshly computed or inherited.)                            --*/

    /* determine beam pattern response (F_plus and F_cross) for given Ifo: */
    XLALComputeDetAMResponse(&Fplus, &Fcross, dataPtr->detector->response,
			     ra, dec, psi, gmst);
		 
    /* signal arrival time (relative to geocenter); */
    timedelay = XLALTimeDelayFromEarthCenter(dataPtr->detector->location,
                                             ra, dec, &GPSlal);
    /* (negative timedelay means signal arrives earlier at Ifo than at geocenter, etc.) */

    /* amount by which to time-shift template (not necessarily same as above "timedelay"): */
    timeshift =  (GPSdouble - (*(REAL8*) LALInferenceGetVariable(dataPtr->modelParams, "time"))) + timedelay;

    /* include distance (overall amplitude) effect in Fplus/Fcross: */
    FplusScaled  = Fplus  / distMpc;
    FcrossScaled = Fcross / distMpc;

    if (LALInferenceCheckVariable(currentParams, "crazyInjectionHLSign") &&
        *((INT4 *)LALInferenceGetVariable(currentParams, "crazyInjectionHLSign"))) {
      if (strstr(dataPtr->name, "H") || strstr(dataPtr->name, "L")) {
        FplusScaled *= -1.0;
        FcrossScaled *= -1.0;
      }
    }

	if(timeWaveform->data->length!=dataPtr->timeModelhPlus->data->length){
		printf("fW%d data%d\n", timeWaveform->data->length, dataPtr->freqModelhPlus->data->length);
		printf("Error!  Time data vector must be same length as original data!\n");
		exit(1);
	}
	
	timeWaveform->deltaT = dataPtr->timeModelhPlus->deltaT;
        
        /* Shift to correct start time. */
        timeWaveform->epoch = dataPtr->timeModelhPlus->epoch;
        XLALGPSAdd(&(timeWaveform->epoch), timeshift);

#ifdef DEBUG
FILE* file=fopen("TempSignal.dat", "w");	
#endif
	for(i=0; i<timeWaveform->data->length; i++){
#ifdef DEBUG
          double t = timeWaveform->deltaT*i + XLALGPSGetREAL8(&(timeWaveform->epoch));
#endif
                timeWaveform->data->data[i] = 
                  FplusScaled*dataPtr->timeModelhPlus->data->data[i] + 
                  FcrossScaled*dataPtr->timeModelhCross->data->data[i];
#ifdef DEBUG
		fprintf(file, "%lg %lg\n", t, timeWaveform->data[i]);
#endif
	}
#ifdef DEBUG
fclose(file);
#endif
	LALInferenceDestroyVariables(&intrinsicParams);
}	
							  						  
REAL8 ComputeFrequencyDomainOverlap(LALInferenceIFOData * dataPtr,
                                    COMPLEX16Vector * freqData1, 
                                    COMPLEX16Vector * freqData2)
{
  int lower, upper, i;
  double deltaT, deltaF;
  
  double overlap=0.0;
  
  /* determine frequency range & loop over frequency bins: */
  deltaT = dataPtr->timeData->deltaT;
  deltaF = 1.0 / (((double)dataPtr->timeData->data->length) * deltaT);
  lower = ceil(dataPtr->fLow / deltaF);
  upper = floor(dataPtr->fHigh / deltaF);
	
  for (i=lower; i<=upper; ++i){  	  	  
    overlap  += ((4.0*deltaF*(freqData1->data[i].re*freqData2->data[i].re+freqData1->data[i].im*freqData2->data[i].im)) 
                 / dataPtr->oneSidedNoisePowerSpectrum->data->data[i]);
  }

  return overlap;
}

REAL8 NullLogLikelihood(LALInferenceIFOData *data)
/*Idential to FreqDomainNullLogLikelihood                        */
{
	REAL8 loglikeli, totalChiSquared=0.0;
	LALInferenceIFOData *ifoPtr=data;
	
	/* loop over data (different interferometers): */
	while (ifoPtr != NULL) {
          ifoPtr->nullloglikelihood = 0.0;
          REAL8 temp = ComputeFrequencyDomainOverlap(ifoPtr, ifoPtr->freqData->data, ifoPtr->freqData->data);
          totalChiSquared+=temp;
          ifoPtr->nullloglikelihood -= 0.5*temp;
		ifoPtr = ifoPtr->next;
	}
	loglikeli = -0.5 * totalChiSquared; // note (again): the log-likelihood is unnormalised!
	return(loglikeli);
}

REAL8 WhitenedTimeDomainOverlap(const REAL8TimeSeries *whitenedData, const REAL8TimeSeries *data) {
  return 2.0*integrateSeriesProduct(whitenedData, data);
}

REAL8 TimeDomainNullLogLikelihood(LALInferenceIFOData *data) {
  REAL8 logL = 0.0;
  LALInferenceIFOData *ifoPtr = data;
  /* UINT4 ifoIndex = 0; */
  /* UINT4 i; */
  /* char fileName[256]; */
  /* FILE *out; */
  
  while (ifoPtr != NULL) {
    ifoPtr->nullloglikelihood = 0.0;
    REAL8 temp = WhitenedTimeDomainOverlap(ifoPtr->whiteTimeData, ifoPtr->windowedTimeData);
    logL += temp;
    ifoPtr->nullloglikelihood -= 0.5*temp;
 
    ifoPtr = ifoPtr->next;
    /* ifoIndex++; */
  }

  logL *= -0.5;
  
  return logL;
}


/* The time-domain weight corresponding to the (two-sided) noise power
   spectrum S(f) is defined to be:

   s(tau) == \int_{-infty}^\infty df \frac{\exp(2 \pi i f tau)}{S(f)}

*/
void PSDToTDW(REAL8TimeSeries *TDW, const REAL8FrequencySeries *PSD, const REAL8FFTPlan *plan,
              const REAL8 fMin, const REAL8 fMax) {
  COMPLEX16FrequencySeries *CPSD = NULL;
  UINT4 i;
  UINT4 PSDLength = TDW->data->length/2 + 1;
  
  if (PSD->data->length != PSDLength) {
    fprintf(stderr, "PSDToTDW: lengths of PSD and TDW do not match (in %s, line %d)", 
            __FILE__, __LINE__);
    exit(1);
  }

  CPSD = 
    XLALCreateCOMPLEX16FrequencySeries(PSD->name, &(PSD->epoch), PSD->f0, PSD->deltaF, &(PSD->sampleUnits), PSD->data->length);

  for (i = 0; i < PSD->data->length; i++) {
    REAL8 f = PSD->f0 + i*PSD->deltaF;

    if (fMin <= f && f <= fMax) {
      CPSD->data->data[i].re = 1.0 / (2.0*PSD->data->data[i]);
      CPSD->data->data[i].im = 0.0;
    } else {
      CPSD->data->data[i].re = 0.0;
      CPSD->data->data[i].im = 0.0;
    }
  }

  XLALREAL8FreqTimeFFT(TDW, CPSD, plan);

  /* FILE *PSDf = fopen("PSD.dat", "w"); */
  /* for (i = 0; i < PSD->data->length; i++) { */
  /*   fprintf(PSDf, "%g %g\n", i*PSD->deltaF, PSD->data->data[i]); */
  /* } */
  /* fclose(PSDf); */

  /* FILE *TDWf = fopen("TDW.dat", "w"); */
  /* for (i = 0; i < TDW->data->length; i++) { */
  /*   fprintf(TDWf, "%g %g\n", i*TDW->deltaT, TDW->data->data[i]); */
  /* } */
  /* fclose(TDWf); */
}

UINT4 nextPowerOfTwo(const UINT4 n) {
  UINT4 np2 = 1;
  
  for (np2 = 1; np2 < n; np2 *= 2) ; /* Keep doubling until >= n */

  return np2;
}

void padREAL8Sequence(REAL8Sequence *padded, const REAL8Sequence *data) {
  if (padded->length < data->length) {
    fprintf(stderr, "padREAL8Sequence: padded sequence too short (in %s, line %d)", __FILE__, __LINE__);
    exit(1);
  }

  memset(padded->data, 0, padded->length*sizeof(padded->data[0]));
  memcpy(padded->data, data->data, data->length*sizeof(data->data[0]));
}

void padWrappedREAL8Sequence(REAL8Sequence *padded, const REAL8Sequence *data) {
  UINT4 i;
  UINT4 np = padded->length;
  UINT4 nd = data->length;

  if (np < nd) {
    fprintf(stderr, "padWrappedREAL8Sequence: padded sequence too short (in %s, line %d)", __FILE__, __LINE__);
    exit(1);
  }

  memset(padded->data, 0, np*sizeof(padded->data[0]));

  padded->data[0] = data->data[0];
  for (i = 1; i <= (nd-1)/2; i++) {
    padded->data[i] = data->data[i]; /* Positive times/frequencies. */
    padded->data[np-i] = data->data[nd-i]; /* Wrapped, negative times/frequencies. */
  }
  if (nd % 2 == 0) { /* If even, take care of singleton positive frequency. */
    padded->data[nd/2] = data->data[nd/2];
  }
}

UINT4 LIGOTimeGPSToNearestIndex(const LIGOTimeGPS *tm, const REAL8TimeSeries *series) {
  REAL8 dT = XLALGPSDiff(tm, &(series->epoch));

  return (UINT4) (round(dT/series->deltaT));
}

REAL8 integrateSeriesProduct(const REAL8TimeSeries *s1, const REAL8TimeSeries *s2) {
  LIGOTimeGPS start, stop;
  LIGOTimeGPS stopS1, stopS2;
  UINT4 i1, i2;
  REAL8 sum = 0.0;
  REAL8 t1Start, t2Start, t, tStop;

  /* Compute stop times. */
  stopS1 = s1->epoch;
  stopS2 = s2->epoch;
  XLALGPSAdd(&stopS1, (s1->data->length-1)*s1->deltaT);
  XLALGPSAdd(&stopS2, (s2->data->length-1)*s2->deltaT);

  /* The start time is the max of the two start times, the stop time
     is the min of the two stop times */
  start = (XLALGPSCmp(&(s1->epoch), &(s2->epoch)) <= 0 ? s2->epoch : s1->epoch); /* Start at max start time. */
  stop = (XLALGPSCmp(&stopS1, &stopS2) <= 0 ? stopS1 : stopS2); /* Stop at min end time. */

  t = 0;
  tStop = XLALGPSDiff(&stop, &start);
  t1Start = XLALGPSDiff(&(s1->epoch), &start);
  t2Start = XLALGPSDiff(&(s2->epoch), &start);
  i1 = LIGOTimeGPSToNearestIndex(&start, s1);
  i2 = LIGOTimeGPSToNearestIndex(&start, s2);
  
  REAL8 *data1 = s1->data->data;
  REAL8 *data2 = s2->data->data;

  do {
    REAL8 nextTime1, nextTime2, nextTime;
    REAL8 dt;
    nextTime1 = t1Start + (i1+0.5)*s1->deltaT;
    nextTime2 = t2Start + (i2+0.5)*s2->deltaT;

    /* Whichever series needs updating first gives us the next time. */
    nextTime = (nextTime1 <= nextTime2 ? nextTime1 : nextTime2);

    /* Ensure we don't go past the stop time. */
    nextTime = (tStop < nextTime ? tStop : nextTime);

    dt = nextTime - t;

    sum += dt*data1[i1]*data2[i2];

    if (nextTime1 == nextTime) {
      i1++;
    }
    if (nextTime2 == nextTime) {
      i2++;
    }

    t = nextTime;    
  } while (t < tStop);

  return sum;
}

void convolveTimeSeries(REAL8TimeSeries *conv, const REAL8TimeSeries *data, const REAL8TimeSeries *response) {
  UINT4 responseSpan = (response->data->length + 1)/2;
  UINT4 paddedLength = nextPowerOfTwo(data->data->length + responseSpan);
  REAL8FFTPlan *fwdPlan = XLALCreateForwardREAL8FFTPlan(paddedLength, 1); /* Actually measure---rely on FFTW to store the best plan for a given length. */
  REAL8FFTPlan *revPlan = XLALCreateReverseREAL8FFTPlan(paddedLength, 1); /* Same. */
  REAL8Sequence *paddedData, *paddedResponse, *paddedConv;
  COMPLEX16Sequence *dataFFT, *responseFFT;
  UINT4 i;

  if (data->deltaT != response->deltaT) {
    fprintf(stderr, "convolveTimeSeries: sample spacings differ: %g vs %g (in %s, line %d)\n", 
            data->deltaT, response->deltaT, __FILE__, __LINE__);
    exit(1);
  }

  if (conv->data->length < data->data->length) {
    fprintf(stderr, "convolveTimeSeries: output length smaller than input length (in %s, line %d)", __FILE__, __LINE__);
    exit(1);
  }

  paddedData = XLALCreateREAL8Sequence(paddedLength);
  paddedResponse = XLALCreateREAL8Sequence(paddedLength);
  paddedConv = XLALCreateREAL8Sequence(paddedLength);

  dataFFT = XLALCreateCOMPLEX16Sequence(paddedLength/2 + 1); /* Exploit R -> C symmetry. */
  responseFFT = XLALCreateCOMPLEX16Sequence(paddedLength/2 + 1);

  padREAL8Sequence(paddedData, data->data);
  padWrappedREAL8Sequence(paddedResponse, response->data);

  XLALREAL8ForwardFFT(dataFFT, paddedData, fwdPlan);
  XLALREAL8ForwardFFT(responseFFT, paddedResponse, fwdPlan);

  for (i = 0; i < paddedLength/2 + 1; i++) {
    /* Store product in dataFFT. */
    double dataRe, dataIm, resRe, resIm;
    dataRe = dataFFT->data[i].re;
    dataIm = dataFFT->data[i].im;
    resRe = responseFFT->data[i].re;
    resIm = responseFFT->data[i].im;

    dataFFT->data[i].re = dataRe*resRe - dataIm*resIm;
    dataFFT->data[i].im = dataRe*resIm + resRe*dataIm;
  }

  XLALREAL8ReverseFFT(paddedConv, dataFFT, revPlan);

  memset(conv->data->data, 0, conv->data->length*sizeof(conv->data->data[0]));
  for (i = 0; i < data->data->length; i++) {
    conv->data->data[i] = conv->deltaT*paddedConv->data[i]/paddedConv->length; /* Normalize */
  }

  strncpy(conv->name, "convolved", LALNameLength);
  conv->epoch = data->epoch;
  conv->deltaT = data->deltaT;
  conv->f0 = data->f0;
  conv->sampleUnits = data->sampleUnits;  

  XLALDestroyREAL8FFTPlan(fwdPlan);
  XLALDestroyREAL8FFTPlan(revPlan);

  XLALDestroyREAL8Sequence(paddedData);
  XLALDestroyREAL8Sequence(paddedResponse);
  XLALDestroyREAL8Sequence(paddedConv);

  XLALDestroyCOMPLEX16Sequence(dataFFT);
  XLALDestroyCOMPLEX16Sequence(responseFFT);
}

void wrappedTimeSeriesToLinearTimeSeries(REAL8TimeSeries *linear, const REAL8TimeSeries *wrapped) {
  UINT4 NNeg, NPos, N, i;

  if (linear->data->length != wrapped->data->length) {
    fprintf(stderr, "wrappedTimeSeriesToLinearTimeSeries: lengths differ (in %s, line %d)",
            __FILE__, __LINE__);
    exit(1);
  }

  N = linear->data->length;
  NNeg = (N-1)/2;
  NPos = N-NNeg-1; /* 1 for the zero component. */

  for (i = 0; i < NNeg; i++) {
    linear->data->data[i] = wrapped->data->data[N-i-1];
  }
  linear->data->data[NNeg] = wrapped->data->data[0];
  for (i = 1; i <= NPos; i++) {
    linear->data->data[NNeg+i] = wrapped->data->data[i];
  }

  linear->epoch = wrapped->epoch;
  linear->deltaT = wrapped->deltaT;
  linear->f0 = wrapped->f0;
  linear->sampleUnits = wrapped->sampleUnits;
  
  /* Adjust start time for linear to account for negative components. */
  XLALGPSAdd(&linear->epoch, -(NNeg*linear->deltaT));
}

void linearTimeSeriesToWrappedTimeSeries(REAL8TimeSeries *wrapped, const REAL8TimeSeries *linear) {
  UINT4 NNeg, NPos, N, i;

  if (wrapped->data->length != linear->data->length) {
    fprintf(stderr, "linearTimeSeriesToWrappedTimeSeries: lengths differ (in %s, line %d)",
            __FILE__, __LINE__);
    exit(1);
  }

  N = wrapped->data->length;
  NNeg = (N-1)/2;
  NPos = N-NNeg-1;

  wrapped->data->data[0] = linear->data->data[NNeg];
  for (i = 1; i <= NPos; i++) {
    wrapped->data->data[i] = linear->data->data[NNeg+i];
  }
  for (i = 0; i < NNeg; i++) {
    wrapped->data->data[N-i-1] = linear->data->data[i];
  }

  wrapped->epoch = linear->epoch;
  wrapped->deltaT = linear->deltaT;
  wrapped->f0 = linear->f0;
  wrapped->sampleUnits = linear->sampleUnits;
  
  /* Adjust start time. */
  XLALGPSAdd(&wrapped->epoch, NNeg*wrapped->deltaT);
}

REAL8 timeDomainOverlap(const REAL8TimeSeries *TDW, const REAL8TimeSeries *A, const REAL8TimeSeries *B) {
  REAL8TimeSeries *Bconv;
  REAL8 overlap;

  Bconv = XLALCreateREAL8TimeSeries(B->name, &(B->epoch), 0.0, B->deltaT, &(B->sampleUnits), B->data->length);

  convolveTimeSeries(Bconv, B, TDW);

  overlap = integrateSeriesProduct(A, Bconv);

  XLALDestroyREAL8TimeSeries(Bconv);

  return 4.0*overlap; /* This is the overlap definition. */
}

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

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInference.h>
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


/** Internal functions which are not exported */
//static UINT4 LIGOTimeGPSToNearestIndex(const LIGOTimeGPS *tm, const REAL8TimeSeries *series);
//static UINT4 LIGOTimeGPSToNearestIndex(const LIGOTimeGPS *tm, const REAL8TimeSeries *series) {
//  REAL8 dT = XLALGPSDiff(tm, &(series->epoch));
//
//  return (UINT4) (round(dT/series->deltaT));
//}

/* The time-domain weight corresponding to the (two-sided) noise power
   spectrum S(f) is defined to be:

   s(tau) == \int_{-infty}^\infty df \frac{\exp(2 \pi i f tau)}{S(f)}

*/

//static UINT4 nextPowerOfTwo(const UINT4 n);
//static UINT4 nextPowerOfTwo(const UINT4 n) {
//  UINT4 np2 = 1;
//  
//  for (np2 = 1; np2 < n; np2 *= 2) ; /* Keep doubling until >= n */
//
//  return np2;
//}

//static void padREAL8Sequence(REAL8Sequence *padded, const REAL8Sequence *data);
//static void padREAL8Sequence(REAL8Sequence *padded, const REAL8Sequence *data) {
//  if (padded->length < data->length) {
//    fprintf(stderr, "padREAL8Sequence: padded sequence too short (in %s, line %d)", __FILE__, __LINE__);
//    exit(1);
//  }
//
//  memset(padded->data, 0, padded->length*sizeof(padded->data[0]));
//  memcpy(padded->data, data->data, data->length*sizeof(data->data[0]));
//}

//static void padWrappedREAL8Sequence(REAL8Sequence *padded, const REAL8Sequence *data);
//static void padWrappedREAL8Sequence(REAL8Sequence *padded, const REAL8Sequence *data) {
//  UINT4 i;
//  UINT4 np = padded->length;
//  UINT4 nd = data->length;
//
//  if (np < nd) {
//    fprintf(stderr, "padWrappedREAL8Sequence: padded sequence too short (in %s, line %d)", __FILE__, __LINE__);
//    exit(1);
//  }
//
//  memset(padded->data, 0, np*sizeof(padded->data[0]));
//
//  padded->data[0] = data->data[0];
//  for (i = 1; i <= (nd-1)/2; i++) {
//    padded->data[i] = data->data[i]; /* Positive times/frequencies. */
//    padded->data[np-i] = data->data[nd-i]; /* Wrapped, negative times/frequencies. */
//  }
//  if (nd % 2 == 0) { /* If even, take care of singleton positive frequency. */
//    padded->data[nd/2] = data->data[nd/2];
//  }
//}

/* Covariance matrix for use in analytic likelihoods */
  static const REAL8 CM[15][15] = {{0.045991865933182365,-0.005489748382557155,-0.0016101718296658227,0.000631074138709115,-0.0005128470284006775,-0.0010749983983025836,-0.001167836861249167,-0.07694897715679858,0.0005260905282822459,0.00013607957548231718,0.0001970785895702776,0.0010537612689396098,-0.0008023132316053279,0.0013831032271139259,-0.0010523947356712027},
    {-0.005489748382557152,0.05478640427684032,-0.0007518149961062337,-0.0024914078235453523,-0.00009338553047165586,0.0015263321896359199,-0.0036655537192248114,0.03169780190169035,0.0006761345004654851,-0.00037599761532668477,0.0005571796842520161,-0.0011168093514619126,-0.0007033670229857137,-0.0035345099781093734,0.00234527085387442},
    {-0.0016101718296658227,-0.0007518149961062335,0.0010936682441487272,-0.00005217477897378074,-0.00024391971119620107,-0.0002384555935682541,-0.00022582577511239625,0.003943831557506022,-0.00015907316035613088,0.00018120630281353927,0.00009068293162194053,0.00013044059594397145,-0.00013678975289569217,0.00022094200221693018,-0.0003061634969958406},
    {0.000631074138709115,-0.0024914078235453523,-0.0000521747789737807,0.002947196332915054,-0.00038508616134763673,-0.00011317420216125936,0.0009796309417027076,-0.0024078472818094164,-0.0001809367534662932,-8.560285041708448e-6,0.00012206188815326092,0.000011035343589960223,0.00032617765310252674,-0.00032497871992205406,-0.000579663407072727},
    {-0.0005128470284006775,-0.00009338553047165613,-0.00024391971119620107,-0.0003850861613476367,0.0013284571608980058,-0.0003674776827218034,-0.0001244386086190235,0.007067905265651232,0.00015065561412672175,1.3528871245427267e-6,0.00005319182621224511,-0.00006417063490265885,0.000015326214000894772,-0.0003258872829973067,-0.0004560360202322925},
    {-0.0010749983983025836,0.0015263321896359199,-0.00023845559356825414,-0.0001131742021612593,-0.0003674776827218035,0.004317612325601109,0.0008845843191354933,-0.034870585314942824,-0.00002079809771946314,-0.00003955745826391885,-0.0002167693575226643,0.00009737526755558733,-0.00012068930992937891,-0.0006370611332578868,0.0003286931361162724},
    {-0.001167836861249168,-0.0036655537192248123,-0.0002258257751123962,0.0009796309417027076,-0.0001244386086190235,0.0008845843191354933,0.003616798210497918,-0.029751915470875086,-0.00019521045260144132,0.00022364294241159524,-0.000021178223746799266,-0.00029625467953365886,0.0002822034243865473,-0.000051234056572972575,-0.000835616647661762},
    {-0.07694897715679858,0.031697801901690345,0.003943831557506019,-0.002407847281809419,0.007067905265651228,-0.03487058531494284,-0.029751915470875107,5.7734267068088,0.005521731225009533,-0.017008048805405164,0.006749693090695894,-0.009972137823002502,-0.01237668867616082,-0.016718782760481683,0.03495582063898905},
    {0.0005260905282822459,0.0006761345004654851,-0.00015907316035613088,-0.00018093675346629324,0.00015065561412672175,-0.000020798097719463143,-0.00019521045260144135,0.005521731225009533,0.0004610670018969681,-0.00010427010812879567,-9.861561285861989e-6,-0.00013973836161050322,-0.00005910518389959231,0.000010588566237447082,-0.00009970008161638036},
    {0.00013607957548231745,-0.00037599761532668477,0.0001812063028135393,-8.560285041708332e-6,1.352887124542701e-6,-0.00003955745826391886,0.00022364294241159524,-0.01700804880540517,-0.00010427010812879569,0.0005909125052583999,0.000021925458163952995,-0.00003232183992063706,-0.00007542207331389262,-0.00044187524915708764,-0.00017652647061896536},
    {0.00019707858957027763,0.0005571796842520161,0.00009068293162194053,0.00012206188815326099,0.0000531918262122451,-0.00021676935752266437,-0.000021178223746799252,0.006749693090695894,-9.861561285862006e-6,0.00002192545816395299,0.00024417715762416555,-0.00004770765187366059,-0.00017551566664082666,-0.000025739468872420018,-0.00022421586669614112},
    {0.0010537612689396102,-0.0011168093514619126,0.0001304405959439715,0.000011035343589960283,-0.00006417063490265878,0.00009737526755558737,-0.0002962546795336587,-0.009972137823002498,-0.00013973836161050322,-0.00003232183992063706,-0.000047707651873660604,0.0011750671719392885,0.00002175988100809145,-0.00003583994582735929,-0.0003198947382559093},
    {-0.0008023132316053279,-0.0007033670229857146,-0.00013678975289569214,0.00032617765310252674,0.00001532621400089473,-0.00012068930992937891,0.0002822034243865474,-0.012376688676160815,-0.000059105183899592275,-0.00007542207331389265,-0.00017551566664082666,0.000021759881008091417,0.0010520988043566145,0.000050310054858351346,0.00016715495303381306},
    {0.001383103227113926,-0.0035345099781093747,0.0002209420022169302,-0.00032497871992205406,-0.0003258872829973067,-0.0006370611332578866,-0.00005123405657297253,-0.016718782760481683,0.000010588566237447109,-0.00044187524915708786,-0.000025739468872420018,-0.00003583994582735929,0.00005031005485835143,0.0051765089984675185,-0.00002532055516308609},
    {-0.0010523947356712027,0.00234527085387442,-0.00030616349699584056,-0.0005796634070727272,-0.0004560360202322925,0.00032869313611627257,-0.000835616647661762,0.03495582063898905,-0.0000997000816163803,-0.0001765264706189654,-0.00022421586669614112,-0.0003198947382559093,0.000167154953033813,-0.00002532055516308609,0.0053108216113592396}};




/* ============ Likelihood computations: ========== */

/** For testing purposes (for instance sampling the prior), likelihood that returns 0.0 = log(1) every
 time.  Activated with the --zeroLogLike command flag. */
REAL8 LALInferenceZeroLogLikelihood(LALInferenceVariables UNUSED *currentParams, LALInferenceIFOData UNUSED *data, LALInferenceTemplateFunction UNUSED *template) {
  return 0.0;
}

REAL8 LALInferenceUndecomposedFreqDomainLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData * data, 
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
  //static int timeDomainWarning = 0;
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
	double mc;
  int different;
	UINT4 logDistFlag=0;
  LALStatus status;
  memset(&status,0,sizeof(status));
  LALInferenceVariables intrinsicParams;

  logDistFlag=LALInferenceCheckVariable(currentParams, "logdistance");
  if(LALInferenceCheckVariable(currentParams,"logmc")){
    mc=exp(*(REAL8 *)LALInferenceGetVariable(currentParams,"logmc"));
    LALInferenceAddVariable(currentParams,"chirpmass",&mc,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
  }

  /* determine source's sky location & orientation parameters: */
  ra        = *(REAL8*) LALInferenceGetVariable(currentParams, "rightascension"); /* radian      */
  dec       = *(REAL8*) LALInferenceGetVariable(currentParams, "declination");    /* radian      */
  psi       = *(REAL8*) LALInferenceGetVariable(currentParams, "polarisation");   /* radian      */
  GPSdouble = *(REAL8*) LALInferenceGetVariable(currentParams, "time");           /* GPS seconds */
	if(logDistFlag)
		 distMpc = exp(*(REAL8*)LALInferenceGetVariable(currentParams,"logdistance"));
	else
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
	if(logDistFlag)
			LALInferenceRemoveVariable(&intrinsicParams, "logdistance");
	else
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
      LALInferenceAddVariable(dataPtr->modelParams, "time", &timeTmp, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
      template(dataPtr);
      if(XLALGetBaseErrno()==XLAL_FAILURE) /* Template generation failed in a known way, set -Inf likelihood */
          return(-DBL_MAX);

      if (dataPtr->modelDomain == LALINFERENCE_DOMAIN_TIME) {
//	if (!timeDomainWarning) {
//	  timeDomainWarning = 1;
//	  fprintf(stderr, "WARNING: using time domain template with frequency domain likelihood (in %s, line %d)\n", __FILE__, __LINE__);
//	}
        LALInferenceExecuteFT(dataPtr);
        /* note that the dataPtr->modelParams "time" element may have changed here!! */
        /* (during "template()" computation)  */
      }
    }
    else { /* no re-computation necessary. Return back "time" value, do nothing else: */
      LALInferenceAddVariable(dataPtr->modelParams, "time", &timeTmp, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
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
    lower = (UINT4)ceil(dataPtr->fLow / deltaF);
    upper = (UINT4)floor(dataPtr->fHigh / deltaF);
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

REAL8 LALInferenceFreqDomainStudentTLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData *data, 
                                      LALInferenceTemplateFunction *template)
{
  //static int timeDomainWarning = 0;
  double Fplus, Fcross;
  double FplusScaled, FcrossScaled;
  double diffRe, diffIm, diffSquared;
  double dataReal, dataImag;
  REAL8 loglikeli;
  REAL8 plainTemplateReal, plainTemplateImag;
  REAL8 templateReal, templateImag;
  int i, lower, upper;
  LALInferenceIFOData *dataPtr;
  double ra, dec, psi, distMpc, gmst,mc;
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
  LALInferenceVariables intrinsicParams;
  
  /* Fill in derived parameters if necessary */
  if(LALInferenceCheckVariable(currentParams,"logdistance")){
    distMpc=exp(*(REAL8 *) LALInferenceGetVariable(currentParams,"logdistance"));
    LALInferenceAddVariable(currentParams,"distance",&distMpc,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
  }

  if(LALInferenceCheckVariable(currentParams,"logmc")){
    mc=exp(*(REAL8 *)LALInferenceGetVariable(currentParams,"logmc"));
    LALInferenceAddVariable(currentParams,"chirpmass",&mc,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
  }
  /* determine source's sky location & orientation parameters: */
  ra        = *(REAL8*) LALInferenceGetVariable(currentParams, "rightascension"); /* radian      */
  dec       = *(REAL8*) LALInferenceGetVariable(currentParams, "declination");    /* radian      */
  psi       = *(REAL8*) LALInferenceGetVariable(currentParams, "polarisation");   /* radian      */
  GPSdouble = *(REAL8*) LALInferenceGetVariable(currentParams, "time");           /* GPS seconds */
  distMpc   = *(REAL8*) LALInferenceGetVariable(currentParams, "distance");       /* Mpc         */

  /* figure out GMST: */
  /* XLALINT8NSToGPS(&GPSlal, floor(1e9 * GPSdouble + 0.5)); */
  XLALGPSSetREAL8(&GPSlal, GPSdouble);
  /* UandA.units    = MST_RAD;                       */
  /* UandA.accuracy = LALLEAPSEC_LOOSE;              */
  /* LALGPStoGMST1(&status, &gmst, &GPSlal, &UandA); */
  gmst=XLALGreenwichMeanSiderealTime(&GPSlal);
  intrinsicParams.head      = NULL;
  intrinsicParams.dimension = 0;
  LALInferenceCopyVariables(currentParams, &intrinsicParams);
  LALInferenceRemoveVariable(&intrinsicParams, "rightascension");
  LALInferenceRemoveVariable(&intrinsicParams, "declination");
  LALInferenceRemoveVariable(&intrinsicParams, "polarisation");
  LALInferenceRemoveVariable(&intrinsicParams, "time");
  LALInferenceRemoveVariable(&intrinsicParams, "distance");
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
      LALInferenceAddVariable(dataPtr->modelParams, "time", &timeTmp, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
      template(dataPtr);
      if (dataPtr->modelDomain == LALINFERENCE_DOMAIN_TIME) {
//	if (!timeDomainWarning) {
//	  timeDomainWarning = 1;
//	  fprintf(stderr, "WARNING: using time domain template with frequency domain likelihood (in %s, line %d)\n", __FILE__, __LINE__);
//	}
        LALInferenceExecuteFT(dataPtr);
        /* note that the dataPtr->modelParams "time" element may have changed here!! */
        /* (during "template()" computation)  */
      }
    }
    else { /* no re-computation necessary. Return back "time" value, do nothing else: */
      LALInferenceAddVariable(dataPtr->modelParams, "time", &timeTmp, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
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

    /* extract the element from the "df" vector that carries the current Ifo's name: */
    CHAR df_variable_name[64];
    sprintf(df_variable_name,"df_%s",dataPtr->name);
    if(LALInferenceCheckVariable(currentParams,df_variable_name)){
      degreesOfFreedom = *(REAL8*) LALInferenceGetVariable(currentParams,df_variable_name);
    }
    else {
      degreesOfFreedom = dataPtr->STDOF;
    }
    if (!(degreesOfFreedom>0)) {
      XLALPrintError(" ERROR in StudentTLogLikelihood(): degrees-of-freedom parameter must be positive.\n");
      XLAL_ERROR_REAL8(XLAL_EDOM);
    }

    /* determine frequency range & loop over frequency bins: */
    deltaT = dataPtr->timeData->deltaT;
    deltaF = 1.0 / (((double)dataPtr->timeData->data->length) * deltaT);
    lower = (UINT4)ceil(dataPtr->fLow / deltaF);
    upper = (UINT4)floor(dataPtr->fHigh / deltaF);
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
      singleFreqBinTerm = ((nu+2.0)/2.0) * log(1.0 + (FourDeltaToverN * diffSquared) / (nu * dataPtr->oneSidedNoisePowerSpectrum->data->data[i]));
      chisquared  += singleFreqBinTerm;   /* (This is a sum-of-squares, or chi^2, term in the Gaussian case, not so much in the Student-t case...)  */
      dataPtr->loglikelihood -= singleFreqBinTerm;
    }
    dataPtr = dataPtr->next;
  }
  loglikeli = -1.0 * chisquared; /* note (again): the log-likelihood is unnormalised! */
  LALInferenceDestroyVariables(&intrinsicParams);  
  return(loglikeli);
}



REAL8 LALInferenceFreqDomainLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData * data, 
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
  COMPLEX16Vector *freqModelResponse=NULL;

  /* loop over data (different interferometers): */
  while (ifoPtr != NULL) {
    ifoPtr->loglikelihood = 0.0;

	if(freqModelResponse==NULL)
		freqModelResponse= XLALCreateCOMPLEX16Vector(ifoPtr->freqData->data->length);
	else
		freqModelResponse= XLALResizeCOMPLEX16Vector(freqModelResponse, ifoPtr->freqData->data->length);
	/*compute the response*/
	LALInferenceComputeFreqDomainResponse(currentParams, ifoPtr, template, freqModelResponse);
	/*if(residual==NULL)
		residual=XLALCreateCOMPLEX16Vector(ifoPtr->freqData->data->length);
	else
		residual=XLALResizeCOMPLEX16Vector(residual, ifoPtr->freqData->data->length);
	
	COMPLEX16VectorSubtract(residual, ifoPtr->freqData->data, freqModelResponse);
	totalChiSquared+=ComputeFrequencyDomainOverlap(ifoPtr, residual, residual); 
	*/
        REAL8 temp = LALInferenceComputeFrequencyDomainOverlap(ifoPtr, ifoPtr->freqData->data, ifoPtr->freqData->data)
          -2.0*LALInferenceComputeFrequencyDomainOverlap(ifoPtr, ifoPtr->freqData->data, freqModelResponse)
          +LALInferenceComputeFrequencyDomainOverlap(ifoPtr, freqModelResponse, freqModelResponse);
	totalChiSquared+=temp;
        ifoPtr->loglikelihood -= 0.5*temp;

    ifoPtr = ifoPtr->next;
  }
  loglikeli = -0.5 * totalChiSquared; // note (again): the log-likelihood is unnormalised!
  XLALDestroyCOMPLEX16Vector(freqModelResponse);
  return(loglikeli);
}

REAL8 LALInferenceChiSquareTest(LALInferenceVariables *currentParams, LALInferenceIFOData * data, LALInferenceTemplateFunction *template)
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
  INT4  i, chisqPt, imax,  kmax, numBins=0;//kmin - set but not used
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
    LALInferenceComputeFreqDomainResponse(currentParams, ifoPtr, template, freqModelResponse);

    deltaT = ifoPtr->timeData->deltaT;
    deltaF = 1.0 / (((REAL8)ifoPtr->timeData->data->length) * deltaT);
    
    /* Store values of fLow and fHigh to use later */
    lowerF = ifoPtr->fLow;
    upperF = ifoPtr->fHigh;
   
    /* Generate bin boundaries */
    numBins = *(INT4*) LALInferenceGetVariable(currentParams, "numbins");
    //kmin = ceil(ifoPtr->fLow / deltaF); - set but not used
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
    
    x = LALInferenceComputeFrequencyDomainOverlap(ifoPtr, ifoPtr->freqData->data, freqModelResponse)/(sqrt(norm));
    
    ChiSquared=0.0;
    
    for (i=0; i < numBins; ++i){
      
      ifoPtr->fLow = chisqBin[i] * deltaF;
      ifoPtr->fHigh = chisqBin[i+1] * deltaF;
      
      xp = LALInferenceComputeFrequencyDomainOverlap(ifoPtr, ifoPtr->freqData->data, freqModelResponse)/(sqrt(norm));
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

//REAL8 LALInferenceTimeDomainLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData * data, 
//                              LALInferenceTemplateFunction *template)
///***************************************************************/
///* Time domain (log-) likelihood function.                     */
///* Returns the non-normalised logarithmic likelihood.          */
///* Mathematically equivalent to Frequency domain likelihood.   */
///* Time domain version of LALInferenceFreqDomainLogLikelihood()*/
///* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
///* Required (`currentParams') parameters are:                  */
///*   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       */
///*   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  */
///*   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        */
///*   - "distance"        (REAL8, Mpc, >0)                      */
///*   - "time"            (REAL8, GPS sec.)                     */
///***************************************************************/
//{
//  REAL8 loglikeli, totalChiSquared=0.0;
//  LALInferenceIFOData *ifoPtr=data;
//  REAL8TimeSeries *timeModelResponse=NULL;
//
//  /* loop over data (different interferometers): */
//  while (ifoPtr != NULL) {
//    ifoPtr->loglikelihood = 0.0;
//
//    if(timeModelResponse==NULL) {
//      timeModelResponse = 
//	XLALCreateREAL8TimeSeries("time detector response", &(ifoPtr->timeData->epoch), 
//				  0.0, ifoPtr->timeData->deltaT,
//				  &lalDimensionlessUnit,
//				  ifoPtr->timeData->data->length);
//    } else if (timeModelResponse->data->length != ifoPtr->timeData->data->length) {
//      /* Cannot resize *up* a time series, so just dealloc and reallocate it. */
//      XLALDestroyREAL8TimeSeries(timeModelResponse);
//      timeModelResponse =                   
//	XLALCreateREAL8TimeSeries("time detector response", &(ifoPtr->timeData->epoch), 
//				  0.0, ifoPtr->timeData->deltaT,
//				  &lalDimensionlessUnit,
//				  ifoPtr->timeData->data->length);
//    }
//     
//    /*compute the response*/
//    LALInferenceComputeTimeDomainResponse(currentParams, ifoPtr, template, timeModelResponse);
//    REAL8 temp = (LALInferenceWhitenedTimeDomainOverlap(ifoPtr->whiteTimeData, ifoPtr->windowedTimeData)
//                  -2.0*LALInferenceWhitenedTimeDomainOverlap(ifoPtr->whiteTimeData, timeModelResponse)
//                  +LALInferenceTimeDomainOverlap(ifoPtr->timeDomainNoiseWeights, timeModelResponse, timeModelResponse));
//    totalChiSquared+=temp;
//    ifoPtr->loglikelihood -= 0.5*temp;
//    
//    ifoPtr = ifoPtr->next;
//  }
//  loglikeli = -0.5*totalChiSquared; 
//  XLALDestroyREAL8TimeSeries(timeModelResponse);
//  return(loglikeli);
//}

void LALInferenceComputeFreqDomainResponse(LALInferenceVariables *currentParams, LALInferenceIFOData * dataPtr, 
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
  //static int timeDomainWarning = 0;

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
		LALInferenceAddVariable(currentParams,"distance",&distMpc,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
	}

	if(LALInferenceCheckVariable(currentParams,"logmc")){
		mc=exp(*(REAL8 *)LALInferenceGetVariable(currentParams,"logmc"));
		LALInferenceAddVariable(currentParams,"chirpmass",&mc,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
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
      LALInferenceAddVariable(dataPtr->modelParams, "time", &timeTmp, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
      template(dataPtr);
      if (dataPtr->modelDomain == LALINFERENCE_DOMAIN_TIME) {
//		  if (!timeDomainWarning) {
//			  timeDomainWarning = 1;
//			  fprintf(stderr, "WARNING: using time domain template with frequency domain likelihood (in %s, line %d)\n", __FILE__, __LINE__);
//		  }
		  LALInferenceExecuteFT(dataPtr);
		  /* note that the dataPtr->modelParams "time" element may have changed here!! */
		  /* (during "template()" computation)                                      */
		}
    }
    else { /* no re-computation necessary. Return back "time" value, do nothing else: */
      LALInferenceAddVariable(dataPtr->modelParams, "time", &timeTmp, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
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

//void LALInferenceComputeTimeDomainResponse(LALInferenceVariables *currentParams, LALInferenceIFOData * dataPtr, 
//                               LALInferenceTemplateFunction *template, REAL8TimeSeries *timeWaveform)
///***************************************************************/
///* Based on ComputeFreqDomainResponse above.                   */
///* Time-domain single-IFO response computation.                */
///* Computes response for a given template.                     */
///* Will re-compute template only if necessary                  */
///* (i.e., if previous, as stored in data->timeModelhCross,     */
///* was based on different parameters or template function).    */
///* Carries out timeshifting for a given detector               */
///* and projection onto this detector.                          */
///* Result stored in timeResponse, assumed to be correctly      */
///* initialized												   */
///* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
///* Required (`currentParams') parameters are:                  */
///*   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       */
///*   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  */
///*   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        */
///*   - "distance"        (REAL8, Mpc, >0)                      */
///*   - "time"            (REAL8, GPS sec.)                     */
///***************************************************************/							  
//{
//  static int freqDomainWarning = 0;
//	double ra, dec, psi, distMpc, gmst;
//	
//	double GPSdouble;
//	double timeTmp;
//	LIGOTimeGPS GPSlal;
//	double timedelay;  /* time delay b/w iterferometer & geocenter w.r.t. sky location */
//	double timeshift;  /* time shift (not necessarily same as above)                   */
//
//	int different;
//	LALInferenceVariables intrinsicParams;
//	LALStatus status;
//	memset(&status,0,sizeof(status));
//
//	double Fplus, Fcross;
//	double FplusScaled, FcrossScaled;
//	UINT4 i;
//	REAL8 mc;
//
//	/* Fill in derived parameters if necessary */
//	if(LALInferenceCheckVariable(currentParams,"logdistance")){
//		distMpc=exp(*(REAL8 *) LALInferenceGetVariable(currentParams,"logdistance"));
//		LALInferenceAddVariable(currentParams,"distance",&distMpc,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
//	}
//
//	if(LALInferenceCheckVariable(currentParams,"logmc")){
//		mc=exp(*(REAL8 *)LALInferenceGetVariable(currentParams,"logmc"));
//		LALInferenceAddVariable(currentParams,"chirpmass",&mc,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
//	}
//		
//	
//	/* determine source's sky location & orientation parameters: */
//	ra        = *(REAL8*) LALInferenceGetVariable(currentParams, "rightascension"); /* radian      */
//	dec       = *(REAL8*) LALInferenceGetVariable(currentParams, "declination");    /* radian      */
//	psi       = *(REAL8*) LALInferenceGetVariable(currentParams, "polarisation");   /* radian      */
//	GPSdouble = *(REAL8*) LALInferenceGetVariable(currentParams, "time");           /* GPS seconds */
//	distMpc   = *(REAL8*) LALInferenceGetVariable(currentParams, "distance");       /* Mpc         */
//		
//	/* figure out GMST: */
//	XLALGPSSetREAL8(&GPSlal, GPSdouble);
//	//UandA.units    = MST_RAD;
//	//UandA.accuracy = LALLEAPSEC_LOOSE;
//	//LALGPStoGMST1(&status, &gmst, &GPSlal, &UandA);
//	gmst=XLALGreenwichMeanSiderealTime(&GPSlal);
//	intrinsicParams.head      = NULL;
//	intrinsicParams.dimension = 0;
//	LALInferenceCopyVariables(currentParams, &intrinsicParams);
//	LALInferenceRemoveVariable(&intrinsicParams, "rightascension");
//	LALInferenceRemoveVariable(&intrinsicParams, "declination");
//	LALInferenceRemoveVariable(&intrinsicParams, "polarisation");
//	LALInferenceRemoveVariable(&intrinsicParams, "time");
//	LALInferenceRemoveVariable(&intrinsicParams, "distance");
//	// TODO: add pointer to template function here.
//	// (otherwise same parameters but different template will lead to no re-computation!!)
//      
//	/* The parameters the response function can handle by itself     */
//    /* (and which shouldn't affect the template function) are        */
//    /* sky location (ra, dec), polarisation and signal arrival time. */
//    /* Note that the template function shifts the waveform to so that*/
//	/* t_c corresponds to the "time" parameter in                    */
//	/* IFOdata->modelParams (set, e.g., from the trigger value).     */
//    
//    /* Compare parameter values with parameter values corresponding  */
//    /* to currently stored template; ignore "time" variable:         */
//    if (LALInferenceCheckVariable(dataPtr->modelParams, "time")) {
//      timeTmp = *(REAL8 *) LALInferenceGetVariable(dataPtr->modelParams, "time");
//      LALInferenceRemoveVariable(dataPtr->modelParams, "time");
//    }
//    else timeTmp = GPSdouble;
//    different = LALInferenceCompareVariables(dataPtr->modelParams, &intrinsicParams);
//    /* "different" now may also mean that "dataPtr->modelParams" */
//    /* wasn't allocated yet (as in the very 1st iteration).      */
//
//    if (different) { /* template needs to be re-computed: */
//      LALInferenceCopyVariables(&intrinsicParams, dataPtr->modelParams);
//      LALInferenceAddVariable(dataPtr->modelParams, "time", &timeTmp, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
//      template(dataPtr);
//      if (dataPtr->modelDomain == LALINFERENCE_DOMAIN_FREQUENCY) {
//	if (!freqDomainWarning) {
//	  freqDomainWarning = 1;
//	  fprintf(stderr, "WARNING: frequency domain template used with time domain calculation (in %s, line %d)\n", __FILE__, __LINE__);
//	}
//        LALInferenceExecuteInvFT(dataPtr);
//	/* note that the dataPtr->modelParams "time" element may have changed here!! */
//	/* (during "template()" computation)  */
//      }
//    }
//    else { /* no re-computation necessary. Return back "time" value, do nothing else: */
//      LALInferenceAddVariable(dataPtr->modelParams, "time", &timeTmp, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
//    }
//
//    /*-- Template is now in dataPtr->timeModelhPlus and dataPtr->timeModelhCross. --*/
//    /*-- (Either freshly computed or inherited.)                            --*/
//
//    /* determine beam pattern response (F_plus and F_cross) for given Ifo: */
//    XLALComputeDetAMResponse(&Fplus, &Fcross, dataPtr->detector->response,
//			     ra, dec, psi, gmst);
//		 
//    /* signal arrival time (relative to geocenter); */
//    timedelay = XLALTimeDelayFromEarthCenter(dataPtr->detector->location,
//                                             ra, dec, &GPSlal);
//    /* (negative timedelay means signal arrives earlier at Ifo than at geocenter, etc.) */
//
//    /* amount by which to time-shift template (not necessarily same as above "timedelay"): */
//    timeshift =  (GPSdouble - (*(REAL8*) LALInferenceGetVariable(dataPtr->modelParams, "time"))) + timedelay;
//
//    /* include distance (overall amplitude) effect in Fplus/Fcross: */
//    FplusScaled  = Fplus  / distMpc;
//    FcrossScaled = Fcross / distMpc;
//
//    if (LALInferenceCheckVariable(currentParams, "crazyInjectionHLSign") &&
//        *((INT4 *)LALInferenceGetVariable(currentParams, "crazyInjectionHLSign"))) {
//      if (strstr(dataPtr->name, "H") || strstr(dataPtr->name, "L")) {
//        FplusScaled *= -1.0;
//        FcrossScaled *= -1.0;
//      }
//    }
//
//	if(timeWaveform->data->length!=dataPtr->timeModelhPlus->data->length){
//		printf("fW%d data%d\n", timeWaveform->data->length, dataPtr->freqModelhPlus->data->length);
//		printf("Error!  Time data vector must be same length as original data!\n");
//		exit(1);
//	}
//	
//	timeWaveform->deltaT = dataPtr->timeModelhPlus->deltaT;
//        
//        /* Shift to correct start time. */
//        timeWaveform->epoch = dataPtr->timeModelhPlus->epoch;
//        XLALGPSAdd(&(timeWaveform->epoch), timeshift);
//
//#ifdef DEBUG
//FILE* file=fopen("TempSignal.dat", "w");	
//#endif
//	for(i=0; i<timeWaveform->data->length; i++){
//#ifdef DEBUG
//          double t = timeWaveform->deltaT*i + XLALGPSGetREAL8(&(timeWaveform->epoch));
//#endif
//                timeWaveform->data->data[i] = 
//                  FplusScaled*dataPtr->timeModelhPlus->data->data[i] + 
//                  FcrossScaled*dataPtr->timeModelhCross->data->data[i];
//#ifdef DEBUG
//		fprintf(file, "%lg %lg\n", t, timeWaveform->data[i]);
//#endif
//	}
//#ifdef DEBUG
//fclose(file);
//#endif
//	LALInferenceDestroyVariables(&intrinsicParams);
//}	
							  						  
REAL8 LALInferenceComputeFrequencyDomainOverlap(LALInferenceIFOData * dataPtr,
                                    COMPLEX16Vector * freqData1, 
                                    COMPLEX16Vector * freqData2)
{
  if (dataPtr==NULL || freqData1 ==NULL || freqData2==NULL){
  	XLAL_ERROR_REAL8(XLAL_EFAULT); 
  	}
  	
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

REAL8 LALInferenceNullLogLikelihood(LALInferenceIFOData *data)
/*Identical to FreqDomainNullLogLikelihood                        */
{
	REAL8 loglikeli, totalChiSquared=0.0;
	LALInferenceIFOData *ifoPtr=data;
	
	/* loop over data (different interferometers): */
	while (ifoPtr != NULL) {
          ifoPtr->nullloglikelihood = 0.0;
          REAL8 temp = LALInferenceComputeFrequencyDomainOverlap(ifoPtr, ifoPtr->freqData->data, ifoPtr->freqData->data);
          totalChiSquared+=temp;
          ifoPtr->nullloglikelihood -= 0.5*temp;
		ifoPtr = ifoPtr->next;
	}
	loglikeli = -0.5 * totalChiSquared; // note (again): the log-likelihood is unnormalised!
	return(loglikeli);
}

//REAL8 LALInferenceWhitenedTimeDomainOverlap(const REAL8TimeSeries *whitenedData, const REAL8TimeSeries *data) {
//  return 2.0*LALInferenceIntegrateSeriesProduct(whitenedData, data);
//}

//REAL8 LALInferenceTimeDomainNullLogLikelihood(LALInferenceIFOData *data) {
//  REAL8 logL = 0.0;
//  LALInferenceIFOData *ifoPtr = data;
//  /* UINT4 ifoIndex = 0; */
//  /* UINT4 i; */
//  /* char fileName[256]; */
//  /* FILE *out; */
//  
//  while (ifoPtr != NULL) {
//    ifoPtr->nullloglikelihood = 0.0;
//    REAL8 temp = LALInferenceWhitenedTimeDomainOverlap(ifoPtr->whiteTimeData, ifoPtr->windowedTimeData);
//    logL += temp;
//    ifoPtr->nullloglikelihood -= 0.5*temp;
// 
//    ifoPtr = ifoPtr->next;
//    /* ifoIndex++; */
//  }
//
//  logL *= -0.5;
//  
//  return logL;
//}



//REAL8 LALInferenceIntegrateSeriesProduct(const REAL8TimeSeries *s1, const REAL8TimeSeries *s2) {
//  if (s1 == NULL || s2 == NULL)
//      XLAL_ERROR_REAL8(XLAL_EFAULT, "Null arguments received.");
//
//  LIGOTimeGPS start, stop;
//  LIGOTimeGPS stopS1, stopS2;
//  UINT4 i1, i2;
//  REAL8 sum = 0.0;
//  REAL8 t1Start, t2Start, t, tStop;
//
//  /* Compute stop times. */
//  stopS1 = s1->epoch;
//  stopS2 = s2->epoch;
//  XLALGPSAdd(&stopS1, (s1->data->length-1)*s1->deltaT);
//  XLALGPSAdd(&stopS2, (s2->data->length-1)*s2->deltaT);
//
//  /* The start time is the max of the two start times, the stop time
//     is the min of the two stop times */
//  start = (XLALGPSCmp(&(s1->epoch), &(s2->epoch)) <= 0 ? s2->epoch : s1->epoch); /* Start at max start time. */
//  stop = (XLALGPSCmp(&stopS1, &stopS2) <= 0 ? stopS1 : stopS2); /* Stop at min end time. */
//
//  t = 0;
//  tStop = XLALGPSDiff(&stop, &start);
//  t1Start = XLALGPSDiff(&(s1->epoch), &start);
//  t2Start = XLALGPSDiff(&(s2->epoch), &start);
//  i1 = LIGOTimeGPSToNearestIndex(&start, s1);
//  i2 = LIGOTimeGPSToNearestIndex(&start, s2);
//  
//  REAL8 *data1 = s1->data->data;
//  REAL8 *data2 = s2->data->data;
//
//  do {
//    REAL8 nextTime1, nextTime2, nextTime;
//    REAL8 dt;
//    nextTime1 = t1Start + (i1+0.5)*s1->deltaT;
//    nextTime2 = t2Start + (i2+0.5)*s2->deltaT;
//
//    /* Whichever series needs updating first gives us the next time. */
//    nextTime = (nextTime1 <= nextTime2 ? nextTime1 : nextTime2);
//
//    /* Ensure we don't go past the stop time. */
//    nextTime = (tStop < nextTime ? tStop : nextTime);
//
//    dt = nextTime - t;
//
//    sum += dt*data1[i1]*data2[i2];
//
//    if (nextTime1 == nextTime) {
//      i1++;
//    }
//    if (nextTime2 == nextTime) {
//      i2++;
//    }
//
//    t = nextTime;    
//  } while (t < tStop);
//
//  return sum;
//}

//void LALInferenceConvolveTimeSeries(REAL8TimeSeries *conv, const REAL8TimeSeries *data, const REAL8TimeSeries *response) {
//  if (conv == NULL || data == NULL || response == NULL)
//    XLAL_ERROR_VOID(XLAL_EFAULT, "Null arguments received.");
//
//  UINT4 responseSpan = (response->data->length + 1)/2;
//  UINT4 paddedLength = nextPowerOfTwo(data->data->length + responseSpan);
//  REAL8FFTPlan *fwdPlan = XLALCreateForwardREAL8FFTPlan(paddedLength, 1); /* Actually measure---rely on FFTW to store the best plan for a given length. */
//  REAL8FFTPlan *revPlan = XLALCreateReverseREAL8FFTPlan(paddedLength, 1); /* Same. */
//  REAL8Sequence *paddedData, *paddedResponse, *paddedConv;
//  COMPLEX16Sequence *dataFFT, *responseFFT;
//  UINT4 i;
//
//  if (data->deltaT != response->deltaT) {
//    fprintf(stderr, "convolveTimeSeries: sample spacings differ: %g vs %g (in %s, line %d)\n", 
//            data->deltaT, response->deltaT, __FILE__, __LINE__);
//    exit(1);
//  }
//
//  if (conv->data->length < data->data->length) {
//    fprintf(stderr, "convolveTimeSeries: output length smaller than input length (in %s, line %d)", __FILE__, __LINE__);
//    exit(1);
//  }
//
//  paddedData = XLALCreateREAL8Sequence(paddedLength);
//  paddedResponse = XLALCreateREAL8Sequence(paddedLength);
//  paddedConv = XLALCreateREAL8Sequence(paddedLength);
//
//  dataFFT = XLALCreateCOMPLEX16Sequence(paddedLength/2 + 1); /* Exploit R -> C symmetry. */
//  responseFFT = XLALCreateCOMPLEX16Sequence(paddedLength/2 + 1);
//
//  padREAL8Sequence(paddedData, data->data);
//  padWrappedREAL8Sequence(paddedResponse, response->data);
//
//  XLALREAL8ForwardFFT(dataFFT, paddedData, fwdPlan);
//  XLALREAL8ForwardFFT(responseFFT, paddedResponse, fwdPlan);
//
//  for (i = 0; i < paddedLength/2 + 1; i++) {
//    /* Store product in dataFFT. */
//    double dataRe, dataIm, resRe, resIm;
//    dataRe = dataFFT->data[i].re;
//    dataIm = dataFFT->data[i].im;
//    resRe = responseFFT->data[i].re;
//    resIm = responseFFT->data[i].im;
//
//    dataFFT->data[i].re = dataRe*resRe - dataIm*resIm;
//    dataFFT->data[i].im = dataRe*resIm + resRe*dataIm;
//  }
//
//  XLALREAL8ReverseFFT(paddedConv, dataFFT, revPlan);
//
//  memset(conv->data->data, 0, conv->data->length*sizeof(conv->data->data[0]));
//  for (i = 0; i < data->data->length; i++) {
//    conv->data->data[i] = conv->deltaT*paddedConv->data[i]/paddedConv->length; /* Normalize */
//  }
//
//  strncpy(conv->name, "convolved", LALNameLength);
//  conv->epoch = data->epoch;
//  conv->deltaT = data->deltaT;
//  conv->f0 = data->f0;
//  conv->sampleUnits = data->sampleUnits;  
//
//  XLALDestroyREAL8FFTPlan(fwdPlan);
//  XLALDestroyREAL8FFTPlan(revPlan);
//
//  XLALDestroyREAL8Sequence(paddedData);
//  XLALDestroyREAL8Sequence(paddedResponse);
//  XLALDestroyREAL8Sequence(paddedConv);
//
//  XLALDestroyCOMPLEX16Sequence(dataFFT);
//  XLALDestroyCOMPLEX16Sequence(responseFFT);
//}

//void LALInferenceWrappedTimeSeriesToLinearTimeSeries(REAL8TimeSeries *linear, const REAL8TimeSeries *wrapped) {
//  UINT4 NNeg, NPos, N, i;
//
//  if (linear->data->length != wrapped->data->length) {
//    fprintf(stderr, "wrappedTimeSeriesToLinearTimeSeries: lengths differ (in %s, line %d)",
//            __FILE__, __LINE__);
//    exit(1);
//  }
//
//  N = linear->data->length;
//  NNeg = (N-1)/2;
//  NPos = N-NNeg-1; /* 1 for the zero component. */
//
//  for (i = 0; i < NNeg; i++) {
//    linear->data->data[i] = wrapped->data->data[N-i-1];
//  }
//  linear->data->data[NNeg] = wrapped->data->data[0];
//  for (i = 1; i <= NPos; i++) {
//    linear->data->data[NNeg+i] = wrapped->data->data[i];
//  }
//
//  linear->epoch = wrapped->epoch;
//  linear->deltaT = wrapped->deltaT;
//  linear->f0 = wrapped->f0;
//  linear->sampleUnits = wrapped->sampleUnits;
//  
//  /* Adjust start time for linear to account for negative components. */
//  XLALGPSAdd(&linear->epoch, -(NNeg*linear->deltaT));
//}

//void LALInferenceLinearTimeSeriesToWrappedTimeSeries(REAL8TimeSeries *wrapped, const REAL8TimeSeries *linear) {
//  UINT4 NNeg, NPos, N, i;
//
//  if (wrapped->data->length != linear->data->length) {
//    fprintf(stderr, "linearTimeSeriesToWrappedTimeSeries: lengths differ (in %s, line %d)",
//            __FILE__, __LINE__);
//    exit(1);
//  }
//
//  N = wrapped->data->length;
//  NNeg = (N-1)/2;
//  NPos = N-NNeg-1;
//
//  wrapped->data->data[0] = linear->data->data[NNeg];
//  for (i = 1; i <= NPos; i++) {
//    wrapped->data->data[i] = linear->data->data[NNeg+i];
//  }
//  for (i = 0; i < NNeg; i++) {
//    wrapped->data->data[N-i-1] = linear->data->data[i];
//  }
//
//  wrapped->epoch = linear->epoch;
//  wrapped->deltaT = linear->deltaT;
//  wrapped->f0 = linear->f0;
//  wrapped->sampleUnits = linear->sampleUnits;
//  
//  /* Adjust start time. */
//  XLALGPSAdd(&wrapped->epoch, NNeg*wrapped->deltaT);
//}

//REAL8 LALInferenceTimeDomainOverlap(const REAL8TimeSeries *TDW, const REAL8TimeSeries *A, const REAL8TimeSeries *B) {
//  REAL8TimeSeries *Bconv;
//  REAL8 overlap;
//
//  Bconv = XLALCreateREAL8TimeSeries(B->name, &(B->epoch), 0.0, B->deltaT, &(B->sampleUnits), B->data->length);
//
//  LALInferenceConvolveTimeSeries(Bconv, B, TDW);
//
//  overlap = LALInferenceIntegrateSeriesProduct(A, Bconv);
//
//  XLALDestroyREAL8TimeSeries(Bconv);
//
//  return 4.0*overlap; /* This is the overlap definition. */
//}

static void extractDimensionlessVariableVector(LALInferenceVariables *currentParams, REAL8 *x, INT4 mode) {
  REAL8 m1, m2, d, iota, phi, psi, ra, dec, t, a1, a2, theta1, theta2, phi1, phi2;
  
  REAL8 mean[15];
  REAL8 Mc;
  
  memset(x, 0, 15*sizeof(REAL8));
  memset(mean, 0, 15*sizeof(REAL8));

  if (mode==0) {
    mean[0] = 16.0;
    mean[1] = 7.0;
    mean[2] = M_PI/2.0;
    mean[3] = M_PI;
    mean[4] = M_PI/2.0;
    mean[5] = M_PI;
    mean[6] = 0.0;
    mean[7] = 50.0;
    mean[8] = 0.0;
    mean[9] =0.5;
    mean[10]=0.5;
    mean[11]=M_PI/2.0;
    mean[12]=M_PI/2.0;
    mean[13]=M_PI;
    mean[14]=M_PI;
  } else if (mode==1) {
    mean[0] = 16.0;
    mean[1] = 7.0;
    mean[2] = 1.0*M_PI/4.0;
    mean[3] = 1.0*M_PI/2.0;
    mean[4] = 1.0*M_PI/4.0;
    mean[5] = 1.0*M_PI/2.0;
    mean[6] = -M_PI/4.0;
    mean[7] = 25.0;
    mean[8] = -0.03;
    mean[9] =0.2;
    mean[10]=0.2;
    mean[11]=1.0*M_PI/4.0;
    mean[12]=1.0*M_PI/4.0;
    mean[13]=1.0*M_PI/2.0;
    mean[14]=1.0*M_PI/2.0;
  } else if (mode==2) {
    /* set means of second mode to be 8 sigma from first mode */
    mean[0] = 16.0 + 8.*sqrt(CM[0][0]);
    mean[1] = 7.0 + 8.*sqrt(CM[1][1]);
    mean[2] = 1.0*M_PI/4.0 + 8.*sqrt(CM[2][2]);
    mean[3] = 1.0*M_PI/2.0 + 8.*sqrt(CM[3][3]);
    mean[4] = 1.0*M_PI/4.0 + 8.*sqrt(CM[4][4]);
    mean[5] = 1.0*M_PI/2.0 + 8.*sqrt(CM[5][5]);
    mean[6] = -M_PI/4.0 + 8.*sqrt(CM[6][6]);
    mean[7] = 25.0 + 8.*sqrt(CM[7][7]);
    mean[8] = -0.03 + 8.*sqrt(CM[8][8]);
    mean[9] =0.2 + 8.*sqrt(CM[9][9]);
    mean[10]=0.2 + 8.*sqrt(CM[10][10]);
    mean[11]=1.0*M_PI/4.0 + 8.*sqrt(CM[11][11]);
    mean[12]=1.0*M_PI/4.0 + 8.*sqrt(CM[12][12]);
    mean[13]=1.0*M_PI/2.0 + 8.*sqrt(CM[13][13]);
    mean[14]=1.0*M_PI/2.0 + 8.*sqrt(CM[14][14]);
    /* mean[0] = 16.0;
    mean[1] = 7.0;
    mean[2] = 3.0*M_PI/4.0;
    mean[3] = 3.0*M_PI/2.0;
    mean[4] = 3.0*M_PI/4.0;
    mean[5] = 3.0*M_PI/2.0;
    mean[6] = M_PI/4.0;
    mean[7] = 75.0;
    mean[8] = 0.03;
    mean[9] =0.8;
    mean[10]=0.8;
    mean[11]=3.0*M_PI/4.0;
    mean[12]=3.0*M_PI/4.0;
    mean[13]=3.0*M_PI/2.0;
    mean[14]=3.0*M_PI/2.0; */
  } else {
    printf("Error!  Unrecognized mode in analytic likelihood!\n");
    exit(1);
  }


  if (LALInferenceCheckVariable(currentParams,"m1")&&LALInferenceCheckVariable(currentParams,"m2"))
  {
    m1=*(REAL8 *)LALInferenceGetVariable(currentParams,"m1");
    m2=*(REAL8 *)LALInferenceGetVariable(currentParams,"m2");
  }
  else
  {
  	if (LALInferenceCheckVariable(currentParams, "chirpmass")) {
    	Mc = *(REAL8 *)LALInferenceGetVariable(currentParams, "chirpmass");
  	} else if (LALInferenceCheckVariable(currentParams, "logmc")) {
    	Mc = exp(*(REAL8 *)LALInferenceGetVariable(currentParams, "logmc"));
  	} else {
    	fprintf(stderr, "Could not find chirpmass or logmc in LALInferenceCorrelatedAnalyticLogLikelihood (in %s, line %d)\n", 
        	    __FILE__, __LINE__);
    	exit(1);
  	}

  	if (LALInferenceCheckVariable(currentParams, "massratio")) {
    	REAL8 eta = *(REAL8 *)LALInferenceGetVariable(currentParams, "massratio");
    	LALInferenceMcEta2Masses(Mc, eta, &m1, &m2);
  	} else if (LALInferenceCheckVariable(currentParams, "asym_massratio")) {
    	REAL8 q = *(REAL8 *)LALInferenceGetVariable(currentParams, "asym_massratio");
    	LALInferenceMcQ2Masses(Mc, q, &m1, &m2);
  	} else {
    	fprintf(stderr, "Could not find eta or q in LALInferenceCorrelatedAnalyticLogLikelihood (in %s, line %d)\n",
        	    __FILE__, __LINE__);
    	exit(1);
  	}
  }

  if (LALInferenceCheckVariable(currentParams, "distance")) {
    d = *(REAL8 *)LALInferenceGetVariable(currentParams, "distance");
  } else if (LALInferenceCheckVariable(currentParams, "logdistance")) {
    d = exp(*(REAL8 *)LALInferenceGetVariable(currentParams, "logdistance"));
  } else {
    fprintf(stderr, "Could not find distance or log(d) in LALInferenceCorrelatedAnalyticLogLikelihood (in %s, line %d)\n",
            __FILE__, __LINE__);
    exit(1);
  }

  iota = *(REAL8 *)LALInferenceGetVariable(currentParams, "inclination");
  psi = *(REAL8 *)LALInferenceGetVariable(currentParams, "polarisation");
  phi = *(REAL8 *)LALInferenceGetVariable(currentParams, "phase");
  ra = *(REAL8 *)LALInferenceGetVariable(currentParams, "rightascension");
  dec = *(REAL8 *)LALInferenceGetVariable(currentParams, "declination");
  t = *(REAL8 *)LALInferenceGetVariable(currentParams, "time");
  
  if (LALInferenceCheckVariable(currentParams, "a_spin1")) {
    a1 = *(REAL8 *)LALInferenceGetVariable(currentParams, "a_spin1");
  } else {
    a1 = 0.0;
  }

  if (LALInferenceCheckVariable(currentParams, "a_spin2")) {
    a2 = *(REAL8 *)LALInferenceGetVariable(currentParams, "a_spin2");
  } else {
    a2 = 0.0;
  }

  if (LALInferenceCheckVariable(currentParams, "phi_spin1")) {
    phi1 = *(REAL8 *)LALInferenceGetVariable(currentParams, "phi_spin1");
  } else {
    phi1 = 0.0;
  }
  
  if (LALInferenceCheckVariable(currentParams, "phi_spin2")) {
    phi2 = *(REAL8 *)LALInferenceGetVariable(currentParams, "phi_spin2");
  } else {
    phi2 = 0.0;
  }
  
  if (LALInferenceCheckVariable(currentParams, "theta_spin1")) {
    theta1 = *(REAL8 *)LALInferenceGetVariable(currentParams, "theta_spin1");
  } else {
    theta1 = 0.0;
  }

  if (LALInferenceCheckVariable(currentParams, "theta_spin2")) {
    theta2 = *(REAL8 *)LALInferenceGetVariable(currentParams, "theta_spin2");
  } else {
    theta2 = 0.0;
  }

  x[0] = (m1    - mean[0]);
  x[1] = (m2    - mean[1]);
  x[2] = (iota  - mean[2])/(M_PI/20.0);
  x[3] = (phi   - mean[3])/(M_PI/10.0);
  x[4] = (psi   - mean[4])/(M_PI/20.0);
  x[5] = (ra    - mean[5])/(M_PI/10.0);
  x[6] = (dec   - mean[6])/(M_PI/10.0);
  x[7] = (d     - mean[7])/10.0;
  x[8] = (t     - mean[8])/0.1;
  x[9] =(a1     - mean[9])/0.1;
  x[10]=(a2     - mean[10])/0.1;
  x[11]=(theta1 - mean[11])/(M_PI/20.0);
  x[12]=(theta2 - mean[12])/(M_PI/20.0);
  x[13]=(phi1   - mean[13])/(M_PI/10.0);
  x[14]=(phi2   - mean[14])/(M_PI/10.0);
}

REAL8 LALInferenceCorrelatedAnalyticLogLikelihood(LALInferenceVariables *currentParams, 
                                                  LALInferenceIFOData UNUSED *data, 
                                                  LALInferenceTemplateFunction UNUSED *template) {
  const INT4 DIM = 15;
  static gsl_matrix *LUCM = NULL;
  static gsl_permutation *LUCMPerm = NULL;
  INT4 mode = 0;
  
  REAL8 x[DIM];
  REAL8 xOrig[DIM];

  gsl_vector_view xView = gsl_vector_view_array(x, DIM);

  extractDimensionlessVariableVector(currentParams, x, mode);

  memcpy(xOrig, x, DIM*sizeof(REAL8));

  if (LUCM==NULL) {
    gsl_matrix_const_view CMView = gsl_matrix_const_view_array(&(CM[0][0]), DIM, DIM);
    int signum;

    LUCM = gsl_matrix_alloc(DIM, DIM);
    LUCMPerm = gsl_permutation_alloc(DIM);

    gsl_matrix_memcpy(LUCM, &(CMView.matrix));

    gsl_linalg_LU_decomp(LUCM, LUCMPerm, &signum);
  }

  gsl_linalg_LU_svx(LUCM, LUCMPerm, &(xView.vector));

  INT4 i;
  REAL8 sum = 0.0;
  for (i = 0; i < DIM; i++) {
    sum += xOrig[i]*x[i];
  }
  //gsl_matrix_free(LUCM);
  //gsl_permutation_free(LUCMPerm);
  return -sum/2.0;
}

REAL8 LALInferenceBimodalCorrelatedAnalyticLogLikelihood(LALInferenceVariables *currentParams,
                                                  LALInferenceIFOData UNUSED *data,
                                                  LALInferenceTemplateFunction UNUSED *template) {
  const INT4 DIM = 15;
  const INT4 MODES = 2;
  INT4 i, mode;
  REAL8 sum = 0.0;
  REAL8 a, b;
  static gsl_matrix *LUCM = NULL;
  static gsl_permutation *LUCMPerm = NULL;
  gsl_vector_view xView;

  REAL8 x[DIM];
  REAL8 xOrig[DIM];
  REAL8 exps[MODES];

  if (LUCM==NULL) {
    gsl_matrix_const_view CMView = gsl_matrix_const_view_array(&(CM[0][0]), DIM, DIM);
    int signum;

    LUCM = gsl_matrix_alloc(DIM, DIM);
    LUCMPerm = gsl_permutation_alloc(DIM);

    gsl_matrix_memcpy(LUCM, &(CMView.matrix));

    gsl_linalg_LU_decomp(LUCM, LUCMPerm, &signum);
  }

  for(mode = 1; mode < 3; mode++) {
    xView = gsl_vector_view_array(x, DIM);

    extractDimensionlessVariableVector(currentParams, x, mode);

    memcpy(xOrig, x, DIM*sizeof(REAL8));

    gsl_linalg_LU_svx(LUCM, LUCMPerm, &(xView.vector));

    sum = 0.0;
    for (i = 0; i < DIM; i++) {
      sum += xOrig[i]*x[i];
    }
    exps[mode-1] = -sum/2.0;
  }

  /* Assumes only two modes used from here on out */
  if (exps[0] > exps[1]) {
    a = exps[0];
    b = exps[1];
  } else {
    a = exps[1];
    b = exps[0];
  }
  //gsl_matrix_free(LUCM);
  //gsl_permutation_free(LUCMPerm);

  /* attempt to keep returned values finite */
  return a + log1p(exp(b-a));
}

REAL8 LALInferenceRosenbrockLogLikelihood(LALInferenceVariables *currentParams,
                                          LALInferenceIFOData UNUSED *data,
                                          LALInferenceTemplateFunction UNUSED *template) {
  const INT4 DIM = 15;
  REAL8 x[DIM];

  REAL8 sum = 0;
  INT4 mode = 0;
  INT4 i;

  extractDimensionlessVariableVector(currentParams, x, mode);

  for (i = 0; i < DIM-1; i++) {
    REAL8 oneMX = 1.0 - x[i];
    REAL8 dx = x[i+1] - x[i]*x[i];

    sum += oneMX*oneMX + 100.0*dx*dx;
  }

  return -sum;
}

//void LALInferencePSDToTDW(REAL8TimeSeries *TDW, const REAL8FrequencySeries *PSD, const REAL8FFTPlan *plan,
//              const REAL8 fMin, const REAL8 fMax) {
//  COMPLEX16FrequencySeries *CPSD = NULL;
//  UINT4 i;
//  UINT4 PSDLength = TDW->data->length/2 + 1;
//  
//  if (PSD->data->length != PSDLength) {
//    fprintf(stderr, "PSDToTDW: lengths of PSD and TDW do not match (in %s, line %d)", 
//            __FILE__, __LINE__);
//    exit(1);
//  }
//
//  CPSD = 
//    XLALCreateCOMPLEX16FrequencySeries(PSD->name, &(PSD->epoch), PSD->f0, PSD->deltaF, &(PSD->sampleUnits), PSD->data->length);
//
//  for (i = 0; i < PSD->data->length; i++) {
//    REAL8 f = PSD->f0 + i*PSD->deltaF;
//
//    if (fMin <= f && f <= fMax) {
//      CPSD->data->data[i].re = 1.0 / (2.0*PSD->data->data[i]);
//      CPSD->data->data[i].im = 0.0;
//    } else {
//      CPSD->data->data[i].re = 0.0;
//      CPSD->data->data[i].im = 0.0;
//    }
//  }
//
//  XLALREAL8FreqTimeFFT(TDW, CPSD, plan);
//
//  /* FILE *PSDf = fopen("PSD.dat", "w"); */
//  /* for (i = 0; i < PSD->data->length; i++) { */
//  /*   fprintf(PSDf, "%g %g\n", i*PSD->deltaF, PSD->data->data[i]); */
//  /* } */
//  /* fclose(PSDf); */
//
//  /* FILE *TDWf = fopen("TDW.dat", "w"); */
//  /* for (i = 0; i < TDW->data->length; i++) { */
//  /*   fprintf(TDWf, "%g %g\n", i*TDW->deltaT, TDW->data->data[i]); */
//  /* } */
//  /* fclose(TDWf); */
//}

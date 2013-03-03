/* 
 *  LALInferenceLikelihood.c:  Bayesian Followup likelihood functions
 *
 *  Copyright (C) 2009 Ilya Mandel, Vivien Raymond, Christian Roever,
 *  Marc van der Sluys and John Veitch, Will M. Farr
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

/* Scaling used for the analytic likelihood parameters */
  static const REAL8 scaling[15] = {
    1.0,
    1.0,
    20.0/M_PI,
    10.0/M_PI,
    20.0/M_PI,
    10.0/M_PI,
    10.0/M_PI,
    0.1,
    10.0,
    10.0,
    10.0,
    20.0/M_PI,
    20.0/M_PI,
    10.0/M_PI,
    10.0/M_PI};

/* Covariance matrix for use in analytic likelihoods */
  static const REAL8 CM[15][15] = {{0.045991865933182365, -0.005489748382557155, -0.01025067223674548, 0.0020087713726603213, -0.0032648855847982987, -0.0034218261781145264, -0.0037173401838545774, -0.007694897715679858, 0.005260905282822458, 0.0013607957548231718, 0.001970785895702776, 0.006708452591621081, -0.005107684668720825, 0.004402554308030673, -0.00334987648531921},
                              {-0.005489748382557152, 0.05478640427684032, -0.004786202916836846, -0.007930397407501268, -0.0005945107515129139, 0.004858466255616657, -0.011667819871670204, 0.003169780190169035, 0.006761345004654851, -0.0037599761532668475, 0.005571796842520162, -0.0071098291510566895, -0.004477773540640284, -0.011250694688474672, 0.007465228985669282},
                              {-0.01025067223674548, -0.004786202916836844, 0.044324704403674524, -0.0010572820723801645, -0.009885693540838514, -0.0048321205972943464, -0.004576186966267275, 0.0025107211483955676, -0.010126911913571181, 0.01153595152487264, 0.005773054728678472, 0.005286558230422045, -0.0055438798694137734, 0.0044772210361854765, -0.00620416958073918},
                              {0.0020087713726603213, -0.007930397407501268, -0.0010572820723801636, 0.029861342087731065, -0.007803477134405363, -0.0011466944120756021, 0.009925736654597632, -0.0007664415942207051, -0.0057593957402320385, -0.00027248233573270216, 0.003885350572544307, 0.00022362281488693097, 0.006609741178005571, -0.003292722856107429, -0.005873218251875897},
                              {-0.0032648855847982987, -0.0005945107515129156, -0.009885693540838514, -0.007803477134405362, 0.0538403407841302, -0.007446654755103316, -0.0025216534232170153, 0.004499568241334517, 0.009591034277125177, 0.00008612746932654223, 0.003386296829505543, -0.002600737873367083, 0.000621148057330571, -0.006603857049454523, -0.009241221870695395},
                              {-0.0034218261781145264, 0.004858466255616657, -0.004832120597294347, -0.0011466944120756015, -0.007446654755103318, 0.043746559133865104, 0.008962713024625965, -0.011099652042761613, -0.0006620240117921668, -0.0012591530037708058, -0.006899982952117269, 0.0019732354732442878, -0.002445676747004324, -0.006454778807421816, 0.0033303577606412765},
                              {-0.00371734018385458, -0.011667819871670206, -0.004576186966267273, 0.009925736654597632, -0.0025216534232170153, 0.008962713024625965, 0.03664582756831382, -0.009470328827284009, -0.006213741694945105, 0.007118775954484294, -0.0006741237990418526, -0.006003374957986355, 0.005718636997353189, -0.0005191095254772077, -0.008466566781233205},
                              {-0.007694897715679857, 0.0031697801901690347, 0.002510721148395566, -0.0007664415942207059, 0.004499568241334515, -0.011099652042761617, -0.009470328827284016, 0.057734267068088, 0.005521731225009532, -0.017008048805405164, 0.006749693090695894, -0.006348460110898, -0.007879244727681924, -0.005321753837620446, 0.011126783289057604},
                              {0.005260905282822458, 0.0067613450046548505, -0.010126911913571181, -0.00575939574023204, 0.009591034277125177, -0.0006620240117921668, -0.006213741694945106, 0.005521731225009532, 0.04610670018969681, -0.010427010812879566, -0.0009861561285861987, -0.008896020395949732, -0.0037627528719902485, 0.00033704453138913093, -0.003173552163182467},
                              {0.0013607957548231744, -0.0037599761532668475, 0.01153595152487264, -0.0002724823357326985, 0.0000861274693265406, -0.0012591530037708062, 0.007118775954484294, -0.01700804880540517, -0.010427010812879568, 0.05909125052583998, 0.002192545816395299, -0.002057672237277737, -0.004801518314458135, -0.014065326026662672, -0.005619012077114913},
                              {0.0019707858957027763, 0.005571796842520162, 0.005773054728678472, 0.003885350572544309, 0.003386296829505542, -0.006899982952117272, -0.0006741237990418522, 0.006749693090695893, -0.0009861561285862005, 0.0021925458163952988, 0.024417715762416557, -0.003037163447600162, -0.011173674374382736, -0.0008193127407211239, -0.007137012700864866},
                              {0.006708452591621083, -0.0071098291510566895, 0.005286558230422046, 0.00022362281488693216, -0.0026007378733670806, 0.0019732354732442886, -0.006003374957986352, -0.006348460110897999, -0.008896020395949732, -0.002057672237277737, -0.003037163447600163, 0.04762367868805726, 0.0008818947598625008, -0.0007262691465810616, -0.006482422704208912},
                              {-0.005107684668720825, -0.0044777735406402895, -0.005543879869413772, 0.006609741178005571, 0.0006211480573305693, -0.002445676747004324, 0.0057186369973531905, -0.00787924472768192, -0.003762752871990247, -0.004801518314458137, -0.011173674374382736, 0.0008818947598624995, 0.042639958466440225, 0.0010194948614718209, 0.0033872675386130637},
                              {0.004402554308030674, -0.011250694688474675, 0.004477221036185477, -0.003292722856107429, -0.006603857049454523, -0.006454778807421815, -0.0005191095254772072, -0.005321753837620446, 0.0003370445313891318, -0.014065326026662679, -0.0008193127407211239, -0.0007262691465810616, 0.0010194948614718226, 0.05244900188599414, -0.000256550861960499},
                              {-0.00334987648531921, 0.007465228985669282, -0.006204169580739178, -0.005873218251875899, -0.009241221870695395, 0.003330357760641278, -0.008466566781233205, 0.011126783289057604, -0.0031735521631824654, -0.005619012077114915, -0.007137012700864866, -0.006482422704208912, 0.0033872675386130632, -0.000256550861960499, 0.05380987317762257}};




/* ============ Likelihood computations: ========== */

/** For testing purposes (for instance sampling the prior), likelihood that returns 0.0 = log(1) every
 time.  Activated with the --zeroLogLike command flag. */
REAL8 LALInferenceZeroLogLikelihood(LALInferenceVariables UNUSED *currentParams, LALInferenceIFOData UNUSED *data, LALInferenceTemplateFunction UNUSED template) {
  return 0.0;
}

REAL8 LALInferenceNoiseOnlyLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData *data, LALInferenceTemplateFunction UNUSED template)
/***************************************************************/
/* (log-) likelihood function.                                 */
/* Returns the non-normalised logarithmic likelihood           */
/* for noise-only models of the data                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Required (`currentParams') parameters are:                  */
/*   - "psdscale"  (gslMatrix)                                 */
/***************************************************************/
{
  double diffRe, diffIm, diffSquared;
  double dataReal, dataImag;
  REAL8 loglikeli;
  int i, j, lower, upper, ifo, n;
  LALInferenceIFOData *dataPtr;
  double chisquared;
  double deltaT, TwoDeltaToverN, deltaF;

  //noise model meta parameters
  gsl_matrix *lines   = NULL;//pointer to matrix holding line centroids
  gsl_matrix *widths  = NULL;//pointer to matrix holding line widths
  gsl_matrix *nparams = NULL;//pointer to matrix holding noise parameters
  double dflog=1.0;        //logarithmic spacing of psd parameters

  int Nblock = 1;            //number of frequency blocks per IFO
  int Nlines = 1;            //number of lines to be removed
  int psdFlag;               //flag for including psd fitting
  int lineFlag;              //flag for excluding lines from integration
  int lineimin = 0;
  int lineimax = 0;
  int lineSwitch;          //switch inside integration to exclude bins

  //line removal parameters
  lineFlag = *((INT4 *)LALInferenceGetVariable(currentParams, "removeLinesFlag"));
  if(lineFlag)
  {
    //Add line matrices to variable lists
    lines  = *(gsl_matrix **)LALInferenceGetVariable(currentParams, "line_center");
    widths = *(gsl_matrix **)LALInferenceGetVariable(currentParams, "line_width");
    Nlines = (int)lines->size2;
  }
  int lines_array[Nlines];
  int widths_array[Nlines];

  //check if psd parameters are included in the model
  psdFlag = *((INT4 *)LALInferenceGetVariable(currentParams, "psdScaleFlag"));
  if(psdFlag)
  {
    //if so, store current noise parameters in easily accessible matrix
    nparams = *((gsl_matrix **)LALInferenceGetVariable(currentParams, "psdscale"));
    Nblock = (int)nparams->size2;

    dflog = *(REAL8*) LALInferenceGetVariable(currentParams, "logdeltaf");
  }
  double alpha[Nblock];
  double lnalpha[Nblock];


  chisquared = 0.0;
  /* loop over data (different interferometers): */
  dataPtr = data;
  ifo=0;

  while (dataPtr != NULL) {
    /* The parameters the Likelihood function can handle by itself   */
    /* (and which shouldn't affect the template function) are        */
    /* sky location (ra, dec), polarisation and signal arrival time. */
    /* Note that the template function shifts the waveform to so that*/
    /* t_c corresponds to the "time" parameter in                    */
    /* IFOdata->modelParams (set, e.g., from the trigger value).     */

    /* Reset log-likelihood */
    dataPtr->loglikelihood = 0.0;

    /* determine frequency range & loop over frequency bins: */
    deltaT = dataPtr->timeData->deltaT;
    deltaF = 1.0 / (((double)dataPtr->timeData->data->length) * deltaT);

    lower = (UINT4)ceil(dataPtr->fLow / deltaF);
    upper = (UINT4)floor(dataPtr->fHigh / deltaF);
    TwoDeltaToverN = 2.0 * deltaT / ((double) dataPtr->timeData->data->length);

    //Set up noise PSD meta parameters
    for(i=0; i<Nblock; i++)
    {
      if(psdFlag)
      {
        alpha[i]   = gsl_matrix_get(nparams,ifo,i);
        lnalpha[i] = log(alpha[i]);
      }
      else
      {
        alpha[i]=1.0;
        lnalpha[i]=0.0;
      }
    }

    //Set up psd line arrays
    for(j=0;j<Nlines;j++)
    {
      if(lineFlag)
      {

        //find range of fourier fourier bins which are excluded from integration
        lines_array[j]  = (int)gsl_matrix_get(lines,ifo,j);//lineimin = (int)(gsl_matrix_get(lines,ifo,j) - gsl_matrix_get(widths,ifo,j));
        widths_array[j] = (int)gsl_matrix_get(widths,ifo,j);//lineimax = (int)(gsl_matrix_get(lines,ifo,j) + gsl_matrix_get(widths,ifo,j));
      }
      else
      {
        lines_array[j]=0;
        widths_array[j]=0;
      }
    }

    for (i=lower; i<=upper; ++i)
    {

      dataReal     = dataPtr->freqData->data->data[i].re / deltaT;
      dataImag     = dataPtr->freqData->data->data[i].im / deltaT;
      
      /* compute squared difference & 'chi-squared': */
      diffRe       = dataReal;         // Difference in real parts...
      diffIm       = dataImag;         // ...and imaginary parts, and...
      diffSquared  = diffRe*diffRe + diffIm*diffIm ;  // ...squared difference of the 2 complex figures.
      
      REAL8 temp = ((TwoDeltaToverN * diffSquared) / dataPtr->oneSidedNoisePowerSpectrum->data->data[i]);

      /* Add noise PSD parameters to the model */
      if(psdFlag)
      {
        n = (int)( log( (double)i/(double)lower )/dflog );
        temp  /= alpha[n];
        temp  += lnalpha[n];
      }

      /* Remove lines from model */
      lineSwitch=1;

      if(lineFlag)
      {
        for(j=0;j<Nlines;j++)
        {
          //find range of fourier fourier bins which are excluded from integration
          lineimin = lines_array[j] - widths_array[j];//gsl_matrix_get(widths,ifo,j));
          lineimax = lines_array[j] + widths_array[j];//(int)(gsl_matrix_get(lines,ifo,j) + gsl_matrix_get(widths,ifo,j));

          //if the current bin is inside the exluded region, set the switch to 0
          if(i>lineimin && i<lineimax) lineSwitch=0;
        }
      }

      /*only sum over bins which are outside of excluded regions */
      if(lineSwitch)
      {
        chisquared  += temp;
        dataPtr->loglikelihood -= temp;
      }
      
     }
    ifo++; //increment IFO counter for noise parameters
    dataPtr = dataPtr->next;
  }

  loglikeli = -1.0 * chisquared; // note (again): the log-likelihood is unnormalised!
  return(loglikeli);
}

REAL8 LALInferenceUndecomposedFreqDomainLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData * data, 
                              LALInferenceTemplateFunction templt)
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
  int i, j, lower, upper, ifo, n;
  LALInferenceIFOData *dataPtr;
  double ra, dec, psi, distMpc, gmst;
  double GPSdouble;
  LIGOTimeGPS GPSlal;
  double chisquared;
  double timedelay;  /* time delay b/w iterferometer & geocenter w.r.t. sky location */
  double timeshift;  /* time shift (not necessarily same as above)                   */
  double deltaT, TwoDeltaToverN, deltaF, twopit, re, im, dre, dim, newRe, newIm;
  double timeTmp;
	double mc;
  int different;
	UINT4 logDistFlag=0;
  LALStatus status;
  memset(&status,0,sizeof(status));
  LALInferenceVariables intrinsicParams;

  //noise model meta parameters
  gsl_matrix *lines   = NULL;//pointer to matrix holding line centroids
  gsl_matrix *widths  = NULL;//pointer to matrix holding line widths
  gsl_matrix *nparams = NULL;//pointer to matrix holding noise parameters
  double dflog = 1.0;        //logarithmic spacing of psd parameters
  
  int Nblock = 1;            //number of frequency blocks per IFO
  int Nlines = 1;            //number of lines to be removed
  int psdFlag = 0;           //flag for including psd fitting
  int lineFlag = 0;          //flag for excluding lines from integration
  int lineimin = 0;
  int lineimax = 0;
  int lineSwitch;          //switch inside integration to exclude bins

  //line removal parameters
  if(LALInferenceCheckVariable(currentParams, "removeLinesFlag"))
    lineFlag = *((INT4 *)LALInferenceGetVariable(currentParams, "removeLinesFlag"));
  if(lineFlag)
  {
    //Add line matrices to variable lists
    lines  = *(gsl_matrix **)LALInferenceGetVariable(currentParams, "line_center");
    widths = *(gsl_matrix **)LALInferenceGetVariable(currentParams, "line_width");
    Nlines = (int)lines->size2;
  }
  int lines_array[Nlines];
  int widths_array[Nlines];

  //check if psd parameters are included in the model
  if(LALInferenceCheckVariable(currentParams, "psdScaleFlag"))
    psdFlag = *((INT4 *)LALInferenceGetVariable(currentParams, "psdScaleFlag"));
  if(psdFlag)
  {
    //if so, store current noise parameters in easily accessible matrix
    nparams = *((gsl_matrix **)LALInferenceGetVariable(currentParams, "psdscale"));
    Nblock = (int)nparams->size2;

    dflog = *(REAL8*) LALInferenceGetVariable(currentParams, "logdeltaf");
  }
  double alpha[Nblock];
  double lnalpha[Nblock];

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
  ifo=0;

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
      templt(dataPtr);
      if(XLALGetBaseErrno()==XLAL_FAILURE) /* Template generation failed in a known way, set -Inf likelihood */
          return(-DBL_MAX);

      if (dataPtr->modelDomain == LAL_SIM_DOMAIN_TIME) {
	/* TD --> FD. */
	LALInferenceExecuteFT(dataPtr);
      }
    }
    else { /* no re-computation necessary. Return back "time" value, do nothing else: */
      LALInferenceAddVariable(dataPtr->modelParams, "time", &timeTmp, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
    }

    /* Template is now in dataPtr->timeFreqModelhPlus and hCross */

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

    /* Employ a trick here for avoiding cos(...) and sin(...) in time
       shifting.  We need to multiply each template frequency bin by
       exp(-J*twopit*deltaF*i) = exp(-J*twopit*deltaF*(i-1)) +
       exp(-J*twopit*deltaF*(i-1))*(exp(-J*twopit*deltaF) - 1) .  This
       recurrance relation has the advantage that the error growth is
       O(sqrt(N)) for N repetitions. */
    
    /* Values for the first iteration: */
    re = cos(twopit*deltaF*lower);
    im = -sin(twopit*deltaF*lower);

    /* Incremental values, using cos(theta) - 1 = -2*sin(theta/2)^2 */
    dim = -sin(twopit*deltaF);
    dre = -2.0*sin(0.5*twopit*deltaF)*sin(0.5*twopit*deltaF);

    //Set up noise PSD meta parameters
    for(i=0; i<Nblock; i++)
    {
      if(psdFlag)
      {
        alpha[i]   = gsl_matrix_get(nparams,ifo,i);
        lnalpha[i] = log(alpha[i]);
      }
      else
      {
        alpha[i]=1.0;
        lnalpha[i]=0.0;
      }
    }

    //Set up psd line arrays
    for(j=0;j<Nlines;j++)
    {
      if(lineFlag)
      {

        //find range of fourier fourier bins which are excluded from integration
        lines_array[j]  = (int)gsl_matrix_get(lines,ifo,j);//lineimin = (int)(gsl_matrix_get(lines,ifo,j) - gsl_matrix_get(widths,ifo,j));
        widths_array[j] = (int)gsl_matrix_get(widths,ifo,j);//lineimax = (int)(gsl_matrix_get(lines,ifo,j) + gsl_matrix_get(widths,ifo,j));
      }
      else
      {
        lines_array[j]=0;
        widths_array[j]=0;
      }
    }


    for (i=lower; i<=upper; ++i){
      /* derive template (involving location/orientation parameters) from given plus/cross waveforms: */
      plainTemplateReal = FplusScaled * dataPtr->freqModelhPlus->data->data[i].re  
                          +  FcrossScaled * dataPtr->freqModelhCross->data->data[i].re;
      plainTemplateImag = FplusScaled * dataPtr->freqModelhPlus->data->data[i].im  
                          +  FcrossScaled * dataPtr->freqModelhCross->data->data[i].im;

      /* do time-shifting...             */
      /* (also un-do 1/deltaT scaling): */
      templateReal = (plainTemplateReal*re - plainTemplateImag*im) / deltaT;
      templateImag = (plainTemplateReal*im + plainTemplateImag*re) / deltaT;
      dataReal     = dataPtr->freqData->data->data[i].re / deltaT;
      dataImag     = dataPtr->freqData->data->data[i].im / deltaT;
      /* compute squared difference & 'chi-squared': */
      diffRe       = dataReal - templateReal;         // Difference in real parts...
      diffIm       = dataImag - templateImag;         // ...and imaginary parts, and...
      diffSquared  = diffRe*diffRe + diffIm*diffIm ;  // ...squared difference of the 2 complex figures.
      REAL8 temp = ((TwoDeltaToverN * diffSquared) / dataPtr->oneSidedNoisePowerSpectrum->data->data[i]);

      /* Add noise PSD parameters to the model */
      if(psdFlag)
      {
        n = (int)( log( (double)i/(double)lower )/dflog );
        temp  /= alpha[n];
        temp  += lnalpha[n];
      }

      /* Remove lines from model */
      lineSwitch=1;
      
      if(lineFlag)
      {
        for(j=0;j<Nlines;j++)
        {
          //find range of fourier fourier bins which are excluded from integration
          lineimin = lines_array[j] - widths_array[j];//gsl_matrix_get(widths,ifo,j));
          lineimax = lines_array[j] + widths_array[j];//(int)(gsl_matrix_get(lines,ifo,j) + gsl_matrix_get(widths,ifo,j));

          //if the current bin is inside the exluded region, set the switch to 0
          if(i>lineimin && i<lineimax) lineSwitch=0;
        }
      }

      /*only sum over bins which are outside of excluded regions */
      if(lineSwitch)
      {
        chisquared  += temp;
        dataPtr->loglikelihood -= temp;
      }
 
      /* Now update re and im for the next iteration. */
      newRe = re + re*dre - im*dim;
      newIm = im + re*dim + im*dre;

      re = newRe;
      im = newIm;
    }
    ifo++; //increment IFO counter for noise parameters
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
                                      LALInferenceTemplateFunction templt)
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
  double deltaT, FourDeltaToverN, deltaF, twopit, re, im, singleFreqBinTerm, dre, dim, newRe, newIm;
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
      templt(dataPtr);

      if (dataPtr->modelDomain == LAL_SIM_DOMAIN_TIME) {
	/* TD --> FD. */
	LALInferenceExecuteFT(dataPtr);
      }
    }
    else { /* no re-computation necessary. Return back "time" value, do nothing else: */
      LALInferenceAddVariable(dataPtr->modelParams, "time", &timeTmp, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
    }

    /* Template is now in dataPtr->freqModelhPlus hCross. */

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

    /* Employ a trick here for avoiding cos(...) and sin(...) in time
       shifting.  We need to multiply each template frequency bin by
       exp(-J*twopit*deltaF*i) = exp(-J*twopit*deltaF*(i-1)) +
       exp(-J*twopit*deltaF*(i-1))*(exp(-J*twopit*deltaF) - 1) .  This
       recurrance relation has the advantage that the error growth is
       O(sqrt(N)) for N repetitions. */
    
    /* Values for the first iteration: */
    re = cos(twopit*deltaF*lower);
    im = -sin(twopit*deltaF*lower);

    /* Incremental values, using cos(theta) - 1 = -2*sin(theta/2)^2 */
    dim = -sin(twopit*deltaF);
    dre = -2.0*sin(0.5*twopit*deltaF)*sin(0.5*twopit*deltaF);

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

      /* Now update re and im for the next iteration. */
      newRe = re + re*dre - im*dim;
      newIm = im + re*dim + im*dre;

      re = newRe;
      im = newIm;
    }
    dataPtr = dataPtr->next;
  }
  loglikeli = -1.0 * chisquared; /* note (again): the log-likelihood is unnormalised! */
  LALInferenceDestroyVariables(&intrinsicParams);  
  return(loglikeli);
}



REAL8 LALInferenceFreqDomainLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData * data, 
                              LALInferenceTemplateFunction templt)
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
	LALInferenceComputeFreqDomainResponse(currentParams, ifoPtr, templt, freqModelResponse);
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

REAL8 LALInferenceChiSquareTest(LALInferenceVariables *currentParams, LALInferenceIFOData * data, LALInferenceTemplateFunction templt)
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
    LALInferenceComputeFreqDomainResponse(currentParams, ifoPtr, templt, freqModelResponse);

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
//                              LALInferenceTemplateFunction templt)
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
                              LALInferenceTemplateFunction templt, COMPLEX16Vector *freqWaveform)
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
	double deltaT, deltaF, twopit, re, im, dre, dim, newRe, newIm;

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
      templt(dataPtr);

      if (dataPtr->modelDomain == LAL_SIM_DOMAIN_TIME) {
	/* TD --> FD. */
	LALInferenceExecuteFT(dataPtr);
      }
    }
    else { /* no re-computation necessary. Return back "time" value, do nothing else: */
      LALInferenceAddVariable(dataPtr->modelParams, "time", &timeTmp, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
    }
    /* Template is now in dataPtr->freqModelhPlus and
       dataPtr->freqModelhCross */

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
/* Employ a trick here for avoiding cos(...) and sin(...) in time
   shifting.  We need to multiply each template frequency bin by
   exp(-J*twopit*deltaF*i) = exp(-J*twopit*deltaF*(i-1)) +
   exp(-J*twopit*deltaF*(i-1))*(exp(-J*twopit*deltaF) - 1) .  This
   recurrance relation has the advantage that the error growth is
   O(sqrt(N)) for N repetitions. */
    
/* Values for the first iteration: */
 re = 1.0;
 im = 0.0;

 /* Incremental values, using cos(theta) - 1 = -2*sin(theta/2)^2 */
 dim = -sin(twopit*deltaF);
 dre = -2.0*sin(0.5*twopit*deltaF)*sin(0.5*twopit*deltaF);

	for(i=0; i<freqWaveform->length; i++){
		/* derive template (involving location/orientation parameters) from given plus/cross waveforms: */
		plainTemplateReal = FplusScaled * dataPtr->freqModelhPlus->data->data[i].re  
                          +  FcrossScaled * dataPtr->freqModelhCross->data->data[i].re;
		plainTemplateImag = FplusScaled * dataPtr->freqModelhPlus->data->data[i].im  
                          +  FcrossScaled * dataPtr->freqModelhCross->data->data[i].im;

		/* do time-shifting...             */
		freqWaveform->data[i].re= (plainTemplateReal*re - plainTemplateImag*im);
		freqWaveform->data[i].im= (plainTemplateReal*im + plainTemplateImag*re);		
#ifdef DEBUG
		fprintf(file, "%lg %lg \t %lg\n", f, freqWaveform->data[i].re, freqWaveform->data[i].im);
#endif
		/* Now update re and im for the next iteration. */
		newRe = re + re*dre - im*dim;
		newIm = im + re*dim + im*dre;

		re = newRe;
		im = newIm;
	}
#ifdef DEBUG
fclose(file);
#endif
	LALInferenceDestroyVariables(&intrinsicParams);
}

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
    mean[0] = 16.0 + 8./scaling[0]*sqrt(CM[0][0]);
    mean[1] = 7.0 + 8./scaling[1]*sqrt(CM[1][1]);
    mean[2] = 1.0*M_PI/4.0 + 8./scaling[2]*sqrt(CM[2][2]);
    mean[3] = 1.0*M_PI/2.0 + 8./scaling[3]*sqrt(CM[3][3]);
    mean[4] = 1.0*M_PI/4.0 + 8./scaling[4]*sqrt(CM[4][4]);
    mean[5] = 1.0*M_PI/2.0 + 8./scaling[5]*sqrt(CM[5][5]);
    mean[6] = -M_PI/4.0 + 8./scaling[6]*sqrt(CM[6][6]);
    mean[7] = 25.0 + 8./scaling[7]*sqrt(CM[7][7]);
    mean[8] = -0.03 + 8./scaling[8]*sqrt(CM[8][8]);
    mean[9] =0.2 + 8./scaling[9]*sqrt(CM[9][9]);
    mean[10]=0.2 + 8./scaling[10]*sqrt(CM[10][10]);
    mean[11]=1.0*M_PI/4.0 + 8./scaling[11]*sqrt(CM[11][11]);
    mean[12]=1.0*M_PI/4.0 + 8./scaling[12]*sqrt(CM[12][12]);
    mean[13]=1.0*M_PI/2.0 + 8./scaling[13]*sqrt(CM[13][13]);
    mean[14]=1.0*M_PI/2.0 + 8./scaling[14]*sqrt(CM[14][14]);
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

  x[0] = scaling[0] * (m1    - mean[0]);
  x[1] = scaling[1] * (m2    - mean[1]);
  x[2] = scaling[2] * (iota  - mean[2]);
  x[3] = scaling[3] * (phi   - mean[3]);
  x[4] = scaling[4] * (psi   - mean[4]);
  x[5] = scaling[5] * (ra    - mean[5]);
  x[6] = scaling[6] * (dec   - mean[6]);
  x[7] = scaling[7] * (d     - mean[7]);
  x[8] = scaling[8] * (t     - mean[8]);
  x[9] = scaling[9] * (a1     - mean[9]);
  x[10]= scaling[10] * (a2     - mean[10]);
  x[11]= scaling[11] * (theta1 - mean[11]);
  x[12]= scaling[12] * (theta2 - mean[12]);
  x[13]= scaling[13] * (phi1   - mean[13]);
  x[14]= scaling[14] * (phi2   - mean[14]);
}

REAL8 LALInferenceCorrelatedAnalyticLogLikelihood(LALInferenceVariables *currentParams, 
                                                  LALInferenceIFOData UNUSED *data, 
                                                  LALInferenceTemplateFunction UNUSED template) {
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
                                                  LALInferenceTemplateFunction UNUSED template) {
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
                                          LALInferenceTemplateFunction UNUSED template) {
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

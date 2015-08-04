/*
 *  LALInferenceCalibrationErrors.c:  Bayesian Followup Calibration Errors routines.
 *
 *  Copyright (C) 2014 Salvatore Vitale, Ryan Lynch
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
#include <lal/LALStdlib.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/StringInput.h>
#include <lal/TimeSeries.h>
#include <lal/LALInference.h>
#include <lal/LALDatatypes.h>
#include <math.h>
#include <lal/LALInferenceCalibrationErrors.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <lal/LALInferenceLikelihood.h>

/* Hard coded fmin and fmax for CE calculation. Those should be large enough to accomodate any realistic CBC WF */
REAL8 freq_min=1.0;
REAL8 freq_max=4096.01;
/* Number of points to sample. Need to be an odd number */
INT4 Npoints= 13;
/* Order of the fit */
INT4 FitOrder = 7;

static REAL8  ConvertRandTransitionSlopeToFunction(REAL8 *coeff, REAL8 f);
static void fill_IFO_Amp_vars_from_IFOname(REAL8 * stddev,REAL8* fbin, char* ifoname);
static void fill_IFO_Pha_vars_from_IFOname(REAL8 * stddev,REAL8* fbin, char* ifoname);
static void CreateRandomAmplitudeCalibrationErrors(REAL8 * ampcoeffs, int calib_seed_ampli, char* ifoname);
static void CreateRandomPhaseCalibrationErrors(REAL8 * phacoeffs, int calib_seed_pha, char* ifoname);
static void FitErrorRealisation(INT4	R,	INT4 N,	REAL8    *y,REAL8 dlogf,REAL8	*D);
static void InvertMatrixSVD (gsl_matrix *A,gsl_matrix	*InvA,	int	N);
static REAL8 ConvertCoefficientsToFunction(REAL8 *coeff, REAL8 f);
static void ApplySquaredAmplitudeErrors(REAL8FrequencySeries * Spectrum,REAL8 * Acoeffs);
static void ApplyBothPhaseAmplitudeErrors(COMPLEX16FrequencySeries *doff,REAL8 * Acoeffs,REAL8 * Pcoeffs);
static void ApplyAmplitudeCalibrationErrors(COMPLEX16FrequencySeries *doff,REAL8 * Acoeffs);
static void ApplyPhaseCalibrationErrors(COMPLEX16FrequencySeries *doff,REAL8 * Pcoeffs);
static void PrintCEtoFile(REAL8* Acoeffs,REAL8* Pcoeffs,LALInferenceIFOData* IFOdata, ProcessParamsTable *commandLine);
static void  fill_IFO_Amp_vars_from_IFOname(REAL8 * stddev,REAL8* fbin, char* ifoname){

    if (!strcmp(ifoname,"H1")){
        stddev[0]=0.104;
        stddev[1]=0.154;
        stddev[2]=0.242;
        fbin[1]=2000.0;
        fbin[2]=4000.0;
    }
    else if(!strcmp(ifoname,"L1")){
        stddev[0]=0.144;
        stddev[1]=0.139;
        stddev[2]=0.138; 
        fbin[1]=2000.0;
        fbin[2]=4000.0;
    }
    else if(!strcmp(ifoname,"V1")){
        stddev[0]=0.10;
        stddev[1]=0.10;
        stddev[2]=0.20;
        fbin[1]=2000.0;
        fbin[2]=4000.0;
    }
    else{
        fprintf(stderr,"Unknown IFO in fill_IFO_vars_from_IFOname! Valid codes are H1, L1, V1. Aborting\n");
        exit(-1);
        }
}

void  fill_IFO_Pha_vars_from_IFOname(REAL8 * stddev,REAL8* fbin, char* ifoname){

    /* Errors here are in degrees. Will convert to rads in CreatePhaseCalibrationErrors */
    if (!strcmp(ifoname,"H1")){
        stddev[0]=4.5;  // 1-500Hz
        stddev[1]=4.5;  // 500-1000Hz
        stddev[2]=4.5;   // 1k-2k
        stddev[3]=4.9;   // 2k -2.8k
        stddev[4]=4.9;   //2.8k-4k
        stddev[5]=5.8;
        fbin[1]=500.0;
        fbin[2]=1000.0;
        fbin[3]=2000.0;
        fbin[4]=2800.0;
        fbin[5]=4000.0;
    }
    else if(!strcmp(ifoname,"L1")){
        stddev[0]=4.2;  // 1-500
        stddev[1]=4.2;  // 500-1000
        stddev[2]=4.2;   // 1k-2k
        stddev[3]=3.6;   // 2k -2.8k
        stddev[4]=3.6;   //2.8k-4k
        stddev[5]=3.3;   // >4k 
        fbin[1]=500.0;
        fbin[2]=1000.0;
        fbin[3]=2000.0;
        fbin[4]=2800.0;
        fbin[5]=4000.0;
    }
    else if(!strcmp(ifoname,"V1")){
        stddev[0]=2.2918;  // 1-500
        stddev[1]=0.5729;  // 500-1000
        stddev[2]=6.87;   // 1k-2k
        stddev[3]=6.87;   // 2k -2.8k
        stddev[4]=360.0*7e-06;   //2.8k-4k
        stddev[5]=360.0*7e-06;   // >4k 
        fbin[1]=500.0;
        fbin[2]=1000.0;
        fbin[3]=2000.0;
        fbin[4]=2800.0;
        fbin[5]=4000.0;
    }
    else{
        fprintf(stderr,"Unknown IFO in fill_IFO_Pha_vars_from_IFOname! Valid codes are H1, L1, V1. Aborting\n");
        exit(-1);
        }
    }

void CreateRandomAmplitudeCalibrationErrors(REAL8 * ampcoeffs, int calib_seed_ampli, char* ifoname){

    /* GSL generator */
    const gsl_rng_type *type;           // RNG type
  	gsl_rng *p;                         // Generator
  	gsl_rng_env_setup();                // Setup environment
  	gsl_rng_default_seed = calib_seed_ampli;    // vary generation sequence
  	type = gsl_rng_default;             // set RNG type to default
  	p = gsl_rng_alloc (type);  
    int i,j;
    REAL8 ampErr[Npoints];

    /* 3 as Amplitude CE are given in three frequency bins */
    REAL8 amp_stdev[3]={0.0};

    /* The stdevs are the errors in the bins [freq_min,freq[bin][0]],[freq[bin][1],freq[bin][2]],.., [freq[bin][last],freq_max]. In Hertz */
    /* plus 2 to accomodate freq_min and freq_max in the first/last position */
    REAL8 amp_freq_bin[2+2]={0.0};
    amp_freq_bin[0]=freq_min;
    amp_freq_bin[3]=freq_max;

    fill_IFO_Amp_vars_from_IFOname(amp_stdev,amp_freq_bin,ifoname);

    /* Space frequencies uniformly in logF */
    REAL8 logF[Npoints];
    REAL8 deltalogf=(log10(freq_max)-log10(freq_min))/(REAL8)(Npoints-1);
    for (i=0; i<Npoints; i++) {
        logF[i]=log10(freq_min)+deltalogf*i;
    }
    int nbins=3;
    /* draw amplitude errors */
    for (i=0; i<Npoints; i++) {
        for(j=0;j<nbins;j++){
            if (logF[i]>=log10(amp_freq_bin[j]) && logF[i]<=log10(amp_freq_bin[j+1])) {
                ampErr[i]=gsl_ran_gaussian(p, amp_stdev[j]);
            }   
        }
    }

    FitErrorRealisation(FitOrder,Npoints,ampErr,deltalogf,ampcoeffs);    
    free(p);
}

void CreateRandomPhaseCalibrationErrors(REAL8 * phacoeffs, int calib_seed_pha, char* ifoname){

    /* GSL generator */
    const gsl_rng_type *type;           // RNG type
  	gsl_rng *p;                         // Generator
  	gsl_rng_env_setup();                // Setup environment
  	gsl_rng_default_seed = calib_seed_pha;    // vary generation sequence
  	type = gsl_rng_default;             // set RNG type to default
  	p = gsl_rng_alloc (type);  
    int i,j;
    REAL8 phaErr[Npoints];
    /* 6 as Phase CE are given in six frequency bins */
    REAL8 pha_stdev[6]={0.0};
    /* The stdevs are the errors in the bins [freq_min,freq[bin][0]],[freq[bin][1],freq[bin][2]],.., [freq[bin][last],freq_max]. In Hertz */
    /* plus 2 to accomodate freq_min and freq_max in the first/last position */
    REAL8 pha_freq_bin[5+2]={0.0};

    pha_freq_bin[0]=freq_min;
    pha_freq_bin[6]=freq_max;

    fill_IFO_Pha_vars_from_IFOname(pha_stdev,pha_freq_bin,ifoname);
    /* Space frequencies uniformly in logF */
    REAL8 logF[Npoints];
    REAL8 deltalogf=(log10(freq_max)-log10(freq_min))/(REAL8)(Npoints-1);
    for (i=0; i<Npoints; i++) {
        logF[i]=log10(freq_min)+deltalogf*i;
    }

    int nbins=6;
    /* draw phase errors */
    for (i=0; i<Npoints; i++) {
        for(j=0;j<nbins;j++){
            if (logF[i]>=log10(pha_freq_bin[j]) && logF[i]<=log10(pha_freq_bin[j+1])) {
                phaErr[i]=gsl_ran_gaussian(p, pha_stdev[j]);
                phaErr[i]*=LAL_PI/180.0;
            }   
        }
    }

    FitErrorRealisation(FitOrder,Npoints,phaErr,deltalogf,phacoeffs);    
    free(p);
}

void ApplyPhaseCalibrationErrors(COMPLEX16FrequencySeries *doff,REAL8 * Pcoeffs){
    /* Take a Complex16FrequencySeries d(f) and a set of phase error coefficients.
     * Return d(f)*exp(I calpha(f))*/

    REAL8 f=0.0;
    /* Amplitude and phase of the WF (in polar coordinates) */
    REAL8 ampli,phase=0.0;
    COMPLEX16 datum;
    REAL8 df=doff->deltaF;
    UINT4 ui;
    for (ui=0;ui<doff->data->length;ui++){
      f=ui*df;
      datum=doff->data->data[ui];
      ampli=sqrt(creal(datum)*creal(datum)+cimag(datum)*cimag(datum));
      phase=atan2(cimag(datum),creal(datum));

      if (Pcoeffs[Npoints+3]==-300.)
        phase+=ConvertRandTransitionSlopeToFunction(Pcoeffs,f);
      else //catch all random errors
        phase+=ConvertCoefficientsToFunction(Pcoeffs,f);
      doff->data->data[ui]=crect(ampli*cos(phase),ampli*sin(phase));
      /* Note: I (salvo) checked that this way of introducing the phase errors does not introduced significant numerical differences w.r.t. expanding both exp(i cE) and datum in their real and imag part, as done in the commented line below and in the likelihood. */
      //doff->data->data[ui]=crect(creal(datum)*cos(tmp)-cimag(datum)*sin(tmp),cos(tmp)*cimag(datum)+sin(tmp)*creal(datum));    
    }
    
}

void ApplyAmplitudeCalibrationErrors(COMPLEX16FrequencySeries *doff,REAL8 * Acoeffs){
    /* Take a Complex16FrequencySeries d(f) and a set of amplitude error coefficients.
     * Return d(f)*calamp(f)  */    
    REAL8 f=0.0;
    /* Amplitude and phase of the WF (in polar coordinates) */
    REAL8 ampli;
    COMPLEX16 datum;
    REAL8 df=doff->deltaF;
    UINT4 ui;
    for (ui=0;ui<doff->data->length;ui++){
        f=ui*df;
        datum=doff->data->data[ui];
        
        if (Acoeffs[Npoints+3]==-300.)
          ampli= 1. + ConvertRandTransitionSlopeToFunction(Acoeffs,f);
        else
          ampli= 1. + ConvertCoefficientsToFunction(Acoeffs,f);

        doff->data->data[ui]=crect(creal(datum)*ampli,cimag(datum)*ampli);     
    }
}

void ApplyBothPhaseAmplitudeErrors(COMPLEX16FrequencySeries *doff,REAL8 * Acoeffs,REAL8 * Pcoeffs){
    /* Apply both phase and amplitude errors to complex 16 stream */
    ApplyAmplitudeCalibrationErrors(doff,Acoeffs);
    ApplyPhaseCalibrationErrors(doff,Pcoeffs);
}

void ApplySquaredAmplitudeErrors(REAL8FrequencySeries * Spectrum,REAL8 * Acoeffs){
    /* Take a REAL8Frequency series S(f) and a set of amplitude error coefficients.
     * Return S(f)*calamp(f)*calamp(f)  */
    REAL8 f=0.0;
    REAL8 ampli;
    REAL8 df=Spectrum->deltaF;
    UINT4 ui;
    for (ui=0;ui<Spectrum->data->length;ui++){
      f=ui*df;
      if (Acoeffs[Npoints+3]==-300.)
        ampli= 1.+ConvertRandTransitionSlopeToFunction(Acoeffs,f);
      else
        ampli= 1. + ConvertCoefficientsToFunction(Acoeffs,f);

      Spectrum->data->data[ui]*=(ampli*ampli);
    }
}

void LALInferenceApplyCalibrationErrors(LALInferenceRunState *state, ProcessParamsTable *commandLine ){
    /*
     * This function takes a pointer to a LALInferenceRunState and applies calibration errors to state->data->freqData (i.e. the frequency domain stream), state->data->oneSidedNoisePowerSpectrum (i.e. the PSD) and state->data->freqData. These arrays must already have been filled by LALInferenceReadData()  and (if injtable is used) by LALInferenceInjectInspiralSignal().
     * CE are either randomly generated or constant.
     * 
     * */
     
     char help[]="\
\n\
------------------------------------------------------------------------------------------------------------------\n\
--- Calibration Errors Handling Arguments ------------------------------------------------------------------------\n\
------------------------------------------------------------------------------------------------------------------\n\
(--AddCalibrationErrors) Adds calibration errors into the f domain datastream (that includes both noise and signal)\n\
(--RandomCE) Add a random realization of phase and amplitude CE, using the S6/VSR2-3 error budget as an indication of the 1-sigma errors\n\
(--ConstantCE) Assumes calibration errors are constant over the bandwidth (requires ConstantCalAmp and ConstantCalPha)\n\
(--ConstantCalAmp [IFO1err,IFO2err,IFO3err]) List of constant amplitude CE. 0.0 means no error, 0.1 means 10 percent\n\
(--ConstantCalPha [IFO1err,IFO2err,IFO3err]) List of constant phase CE. 0.0 means no error, 5 means a  5 degree shift \n\
(--RandomLinearCE ) Assumes CE are given by a contant plateau plus a random jittering of a few percent.\n\t\t After a given frequency f CE increase linearly with a given slope (requires RandomLinearCalAmp and RandomLinearCalPha)\n\
(--RandomLinearCalAmp [IF01_c,IFO1_f,IFO1_slope, ...] ) Add on the i-th IFO's stream errors on the form (IFOi_c + jitter) for f<IFOi_f and (IFOi_c-f)*IFOi_slope for f>IFOi_f\n\
(--RandomLinearCalPha [IF01_c, IFO1_f,IFO1_slope, ...] ) Add on the i-th IFO's stream errors on the form (IFOi_c + jitter) for f<IFOi_f and (IFOi_c-f)*IFOi_slope for f>IFOi_f\n\
 * Constant Calibration Model \n\
  (--MarginalizeConstantCalAmp ) If given, will add a constant value of Amplitude CE per each IFO on the top of the CBC parameters.\n\
  (--MarginalizeConstantCalPha ) If given, will add a constant value of Phase CE per each IFO on the top of the CBC parameters.\n\
  (--constcal_ampsigma ) If given, will use gaussian prior on the constant amplitude error with this sigma (e.g. 0.05=5%) .\n\
  (--constcal_phasigma ) If given, will use gaussian prior on the constant phase error with this sigma (e.g. 5=5degs).\n\
 * Spline Calibration Model \n\
  (--enable-spline-calibration)            Enable cubic-spline calibration error model.\n\
  (--spcal-nodes N)           Set the number of spline nodes per detector (default 5)\n\
  (--spcal-amp-uncertainty X) Set the prior on relative amplitude uncertainty (default 0.1)\n\
  (--spcal-phase-uncertainty X) Set the prior on phase uncertanity in degrees (default 5)\n\n\n";

    static LALStatus   status;
      /* Print command line arguments if state was not allocated */
    if(state==NULL)
    {
      fprintf(stdout,"%s",help);
      return ;
    }
    /* Print command line arguments if help requested */
    if(LALInferenceGetProcParamVal(state->commandLine,"--help"))
    {
      fprintf(stdout,"%s",help);
      return;
    }

    ProcessParamsTable *ppt=NULL;

    if(!LALInferenceGetProcParamVal(commandLine,"--AddCalibrationErrors")) 
    {
      fprintf(stdout,"No --AddCalibrationErrors option give. Not applying calibration errors in injection...\n"); 
      return;
    }

    /* Set calibration seed for random errors */
    if(!LALInferenceGetProcParamVal(commandLine,"--dataseed")){
      fprintf(stdout,"--dataseed is required when running with --AddCalibrationErrors\n");
      exit(1);
    }
    int dataseed=atoi(LALInferenceGetProcParamVal(commandLine,"--dataseed")->value);
    RandomParams *datarandparam=XLALCreateRandomParams(dataseed);
    
    int calib_seed_ampli=0.0;
    int calib_seed_phase=0.0;
    REAL4 tmpampli,tmpphase;
    LALUniformDeviate(&status,&tmpampli,datarandparam);
    LALUniformDeviate(&status,&tmpphase,datarandparam);
    calib_seed_ampli=floor(1E6*tmpampli);
    calib_seed_phase=floor(1E6*tmpphase);
    fprintf(stdout,"Using calibseedAmp %d and calibseedPha %d\n",calib_seed_ampli,calib_seed_phase);
  
    LALInferenceIFOData * tmpdata=state->data;
    int num_ifos=0;
    int i;
    UINT4 unum_ifos=(UINT4) num_ifos;

    while (tmpdata!=NULL) {
      num_ifos++;
      tmpdata=tmpdata->next;
    }
    tmpdata=state->data;
    int this_ifo=0;  
    REAL8 phaseCoeffs[num_ifos][Npoints*2];
    REAL8 ampCoeffs[num_ifos][Npoints*2];
    while (tmpdata!=NULL) {
      memset(ampCoeffs[this_ifo],0.0,2*Npoints*sizeof(REAL8));
      memset(phaseCoeffs[this_ifo],0.0,2*Npoints*sizeof(REAL8));
      this_ifo++;
      tmpdata=tmpdata->next;
    }
    if(LALInferenceGetProcParamVal(commandLine,"--RandomCE")){
        /* Random phase and amplitude calibration errors. Use S6 Budget. That is an overkill, but may be good to test worst case scenarios. */
        fprintf(stdout,"Applying random phase and amplitude errors. \n");
        tmpdata=state->data;
        this_ifo=0;
        /* For each IFO create a CE realization and write in amp/phaseCoeffs the coefficients of the polynomial expansion */
        while (tmpdata!=NULL){
          CreateRandomAmplitudeCalibrationErrors(ampCoeffs[this_ifo],calib_seed_ampli,tmpdata->name);
          CreateRandomPhaseCalibrationErrors(phaseCoeffs[this_ifo],calib_seed_phase,tmpdata->name);
          LALUniformDeviate(&status,&tmpampli,datarandparam);
          LALUniformDeviate(&status,&tmpphase,datarandparam);
          calib_seed_ampli+=floor(1E6*tmpampli);
          calib_seed_phase+=floor(1E6*tmpphase);
          this_ifo++;
          tmpdata=tmpdata->next;
        }
    }
    else if (LALInferenceGetProcParamVal(commandLine,"--ConstantCE")){
      /* NOTE: If we want to apply constant CE we simply have to set ampCoeffs and phaseCoeffs in such a way that the 0th element is non null, while the others are zero.     */ 
        ppt=LALInferenceGetProcParamVal(commandLine,"--ConstantCalAmp");
        if(!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--constantcalamp");
        if (!ppt){ fprintf(stderr,"Must provide a list of constant amplitude calibration errors. E.g: --constantcalamp [1.1,1.2,0.9]. Exiting... \n"); exit(1);}
        else
          fprintf(stdout,"Applying constant amplitude calibration errors. \n");

        char** calamps;
        LALInferenceParseCharacterOptionString(ppt->value,&calamps,&unum_ifos);
        ppt=LALInferenceGetProcParamVal(commandLine,"--ConstantCalPha");
        if(!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--constantcalpha");
        if (!ppt){ fprintf(stderr,"Must provide a list of constant phase calibration errors [Degs]. E.g: --constantcalpha [4.0,2.2,0.1]. Exiting... \n"); exit(1);}
        else
          fprintf(stdout,"Applying constant phase calibration errors. \n");
          
        char** calphases;               
        LALInferenceParseCharacterOptionString(ppt->value,&calphases,&unum_ifos);

        for(i=0;i<num_ifos;i++){
            (ampCoeffs[i])[0]=atof(calamps[i]);
            (phaseCoeffs[i])[0]=atof(calphases[i])*LAL_PI/180.0;
         }
    }
    else if (LALInferenceGetProcParamVal(commandLine,"--RandomLinearCE")){
	
      ppt=LALInferenceGetProcParamVal(commandLine,"--RandomLinearCalAmp");
      if (!ppt){ 
        fprintf(stderr,"Must provide a list of coefficients with random linear amplitude calibration errors. E.g: --RandomLinearCalamp [1.1,1000,0.5, 0.95, 1100,0.5,...]. Exiting... \n"); 
        exit(1);
      }
      else
        fprintf(stdout,"Applying quasi constant amplitude calibration errors. \n");

      char** calamps;
      UINT4 threeTimesNIFO=3*unum_ifos;
      LALInferenceParseCharacterOptionString(ppt->value,&calamps,&threeTimesNIFO);
      i=0;
      tmpdata=state->data;
      while (tmpdata!=NULL){
        /* Store variables for random jitter in the first (Npoints-1) positions of ampCoeffs. 
         * Store constant plateau, knee position and slope in the next 3 positions.
         * Store -300 in the last position to make the code recognize our choice  */
         
        /* Fill random part. Will take 10% of it later on */
        CreateRandomAmplitudeCalibrationErrors(ampCoeffs[i],calib_seed_ampli,tmpdata->name);
        LALUniformDeviate(&status,&tmpampli,datarandparam);
        calib_seed_ampli+=floor(1E6*tmpampli);
        /* Consant plateau, knee, slope*/
        (ampCoeffs[i])[Npoints]=atof(calamps[i*3]);
        (ampCoeffs[i])[Npoints+1]=atof(calamps[i*3+1]);
        (ampCoeffs[i])[Npoints+2]=atof(calamps[i*3+2]);
        (ampCoeffs[i])[Npoints+3]=-300;
        i++;
        tmpdata=tmpdata->next;
       }
      ppt=LALInferenceGetProcParamVal(commandLine,"--RandomLinearCalPha");
      if (!ppt){ 
        fprintf(stderr,"Must provide a list of coefficients with random linear phase calibration errors. E.g: --RandomLinearCalpha [5,1000,0.1, 3, 1100,0.01,...]. Exiting... \n"); 
        exit(1);
      }
      else
        fprintf(stdout,"Applying quasi constant phase calibration errors. \n");
        
      char** calphas;
      LALInferenceParseCharacterOptionString(ppt->value,&calphas,&threeTimesNIFO);
      i=0;
      tmpdata=state->data;
      while (tmpdata!=NULL){
        /* Store variables for random jitter in the first (Npoints-1) positions of phaCoeffs. 
         * Store constant plateau, knee position and slope in the next 3 positions.
         * Store -300 in the last position to make the code recognize our choice  */
         
        /* Fill random part. Will take 10% of it later on */
        CreateRandomPhaseCalibrationErrors(phaseCoeffs[i],calib_seed_phase,tmpdata->name);
        
        LALUniformDeviate(&status,&tmpphase,datarandparam);
        calib_seed_phase+=floor(1E6*tmpphase);
        /* Consant plateau, knee, slope*/
        // user gave degrees, convert to radiands
        (phaseCoeffs[i])[Npoints]=LAL_PI/180.*atof(calphas[i*3]);
        // transition frequency (Hz)
        (phaseCoeffs[i])[Npoints+1]=atof(calphas[i*3+1]);
        // slope. After the transition freq the errors go like f^slope
        (phaseCoeffs[i])[Npoints+2]=atof(calphas[i*3+2]);
        (phaseCoeffs[i])[Npoints+3]=-300;
        i++;
        tmpdata=tmpdata->next;
       }
    }
    else{
      fprintf(stderr, "Must provide a calibration error flag together with --AddCalibrationErrors\n");
      exit(1);
    }
    /* Now apply CE to various quantities */
    tmpdata=state->data;
    this_ifo=0;
    while (tmpdata!=NULL){
      PrintCEtoFile(ampCoeffs[this_ifo],phaseCoeffs[this_ifo],tmpdata, commandLine);
      ApplyBothPhaseAmplitudeErrors(tmpdata->freqData,ampCoeffs[this_ifo],phaseCoeffs[this_ifo]);
      ApplyBothPhaseAmplitudeErrors(tmpdata->whiteFreqData,ampCoeffs[this_ifo],phaseCoeffs[this_ifo]);
      ApplySquaredAmplitudeErrors(tmpdata->oneSidedNoisePowerSpectrum,ampCoeffs[this_ifo]);
      this_ifo++;
      tmpdata=tmpdata->next;
    }
    XLALDestroyRandomParams(datarandparam);
}

void PrintCEtoFile(REAL8* Acoeffs,REAL8* Pcoeffs,LALInferenceIFOData* IFOdata, ProcessParamsTable *commandLine ){

    REAL8 f=0.0;
    UINT4 ui;
    REAL8 df=IFOdata->freqData->deltaF;
    UINT4 f_low_idx=ceil( IFOdata->fLow/df);
    UINT4 f_high_idx=floor(IFOdata->fHigh/df);
    if (df*f_low_idx<freq_min || df*f_high_idx >freq_max) {
        fprintf(stderr,"The min and max frequency in LALInspiralCalibrationErrors.c are inside the range [flow,fmin] of the integral overlap. Exiting...\n");
        exit(1);
    }
    
    ProcessParamsTable *ppt_order=NULL;
    ppt_order=LALInferenceGetProcParamVal(commandLine, "--outfile");
    char *outfile=ppt_order->value;
    
    FILE *calibout;
    char caliboutname[100];
    sprintf(caliboutname,"%s_CE_%s.dat",outfile, IFOdata->name);
    calibout=fopen(caliboutname,"w");
    
    for(ui=f_low_idx;ui<f_high_idx;ui++){
      f=ui*df;
      if (Acoeffs[3+Npoints]==-300.)
        fprintf(calibout,"%lf \t%10.10e \t%10.10e\n",f,ConvertRandTransitionSlopeToFunction(Acoeffs,f),ConvertRandTransitionSlopeToFunction(Pcoeffs,f));
      else 
        fprintf(calibout,"%lf \t%10.10e \t%10.10e\n",f,ConvertCoefficientsToFunction(Acoeffs,f),ConvertCoefficientsToFunction(Pcoeffs,f));
    }
    fclose(calibout);
}

void FitErrorRealisation(INT4	R,	INT4 N,	REAL8    *y,REAL8 dlogf,REAL8	*D)
{
	
  int i=0; int j=0;
  /*********************************************************************
   *
   * Savitsky-Golay Filter
   * - Based on least square fitting of a polynomial of Rth order
   * - Smoothens function by extrapolating from m neighbouring points
   * - See Abraham Savitsky and Marcel J. E. Golay 
   * 		"Smoothing and differentiation of data by simplified least squares procedures"
   * 
   ********************************************************************/ 
	
	/*********************************************************************
	 * 
	 * LAL error handling
	 * 
	 *********************************************************************/	
  
  //INITSTATUS( status, "LALSavitskyGolayFilter", LALCALIBRATIONERRORSC);
  //ATTATCHSTATUSPTR(status);	 
	
	/*********************************************************************
	 * 
	 * Read input
	 * 
	 *********************************************************************/	   
	
  /* LOAD IN TIMESERIES PARAMETERS */
  INT4 M = (N-1)/2;
  //printf("Input parameters: R = %d | M = %d | N = %d | dt = %e \n", R, M, N, dt);	
	
  /*********************************************************************
   *
   * Create Temporary Variables
   * 
   ********************************************************************/  
  
  //printf("Initialising Variables \n");
  
  /* COUNTERS */
  int k;
  
  /* factorial of D (used for derivatives) */
  INT4 factorial = 1;
  
  /* MATRICES AND VECTORS */
	gsl_matrix *m     	= gsl_matrix_calloc (R+1, 2*M+1);   /* m_ij = j^i */
	gsl_matrix *U				= gsl_matrix_calloc (R+1, R+1);		/* U_ii = deltaT^i */
	gsl_vector *a				= gsl_vector_calloc (R+1);			/* a_j, for y(t_i) = Sum_{j=0}^{R} a_j t_i^j */
	gsl_matrix *c				= gsl_matrix_calloc (R+1, 2*M+1);		/* c_ij = U_-1 (m m^T)^-1 m, in a_j = c_ji y_i */ 
	gsl_vector *ym			= gsl_vector_calloc (2*M+1);	/* y_m = [y_-M, ... , y_M] */
	gsl_matrix *tmr			= gsl_matrix_calloc (R+1, 2*M+1);		/* t_m^r = U*m */
	
	/* COMBINED MATRICES AND VECTORS */
	gsl_matrix *mT			= gsl_matrix_calloc (2*M+1, R+1);		/* m^T */
	gsl_matrix *mmT			= gsl_matrix_calloc (R+1, R+1);		/* mm^T */
	gsl_matrix *InvmmT	= gsl_matrix_calloc (R+1, R+1);		/* (mm^T)^-1 */
	gsl_matrix *InvmmTm	= gsl_matrix_calloc (R+1, 2*M+1);		/* (mm^T)^-1 m */
	gsl_matrix *InvU		= gsl_matrix_calloc (R+1, R+1);		/* U^-1 */
	
  /*********************************************************************
   *
   * Filling matrices
   * 
   ********************************************************************/ 
	//printf("Filling parameters \n");
	
  /* m_ij = j^i */
  //printf("Filling parameters -  m %dx%d\n", m->size1, m->size2);
  for(i=0;i<(R+1);i++)
  {
		for(j=0;j<(2*M+1);j++)
		{
			gsl_matrix_set(m, i, j, pow((j-M), i));
		}
	}
	
  //printf("m %dx%d\n", m->size1, m->size2);
	//for(i=0;i<((R+1)*(2*M+1));i++)
	//{
	//printf("%e", gsl_matrix_get(m, i/(2*M+1), i%(2*M+1)));
	//if(i%(2*M+1)==(2*M)) printf("\n");
	//else printf("\t");
	//}  
	//printf("\n");	 
	
  /* U_ii = deltaT^i */
  //printf("Filling parameters -  U %dx%d \n", U->size1, U->size2);
  for(i=0;i<(R+1);i++)
  {
		gsl_matrix_set(U, i, i, pow(dlogf,i));
	} 
	
  //printf("U %dx%d\n", U->size1, U->size2);
	//for(i=0;i<((R+1)*(R+1));i++)
	//{
	//printf("%e", gsl_matrix_get(U, i/(R+1), i%(R+1)));
	//if(i%(R+1)==(R)) printf("\n");
	//else printf("\t");
	//}  
	//printf("\n");	 	
	
	/* m^T */
	//printf("Filling parameters -  mT %dx%d\n", mT->size1, mT->size2);
	for(i=0;i<(R+1); i++)
	{
		for(j=0;j<(2*M+1);j++)
		{
			gsl_matrix_set(mT, j, i, gsl_matrix_get(m, i, j));
		}
	}
	
  //printf("mT %dx%d\n", mT->size1, mT->size2);
	//for(i=0;i<((2*M+1)*(R+1));i++)
	//{
	//printf("%e", gsl_matrix_get(mT, i/(R+1), i%(R+1)));
	//if(i%(R+1)==(R)) printf("\n");
	//else printf("\t");
	//}  
	//printf("\n");	 	
	
	/* mm^T */
	//printf("Filling parameters -  mmT %dx%d\n", mmT->size1, mmT->size2);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
									1.0, m, mT,
									0.0, mmT);	
	
  //printf("mmT %dx%d\n", mmT->size1, mmT->size2);
	//for(i=0;i<((R+1)*(R+1));i++)
	//{
	//printf("%e", gsl_matrix_get(mmT, i/(R+1), i%(R+1)));
	//if(i%(R+1)==(R)) printf("\n");
	//else printf("\t");
	//}  
	//printf("\n");							
	
	/* (mm^T)^-1 */
	//printf("Filling parameters -  InvmmT %dx%d\n", InvmmT->size1, InvmmT->size2);
	InvertMatrixSVD(mmT, InvmmT, R+1);
	
	/* U^-1*/
	//printf("Filling parameters -  InvU %dx%d\n", InvU->size1, InvU->size2);
	//InvertMatrixSVD(U, InvU, R+1);
	
	for(i=0;i<(R+1);i++)
	{
		//printf("%e | %e \n", 1.0/gsl_matrix_get(U, i, i), gsl_matrix_get(InvU, i, i));
		gsl_matrix_set(InvU, i, i, 1.0/gsl_matrix_get(U, i, i));
	}
	
	/* (mm^T)^-1 m */
	//printf("Filling parameters -  InvmmTm %dx%d \n", InvmmTm->size1, InvmmTm->size2 );
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
									1.0, InvmmT, m,
									0.0, InvmmTm);	
	
  //printf("InvmmTm \n");
	//for(i=0;i<((R+1)*(2*M+1));i++)
	//{
	//printf("%e", gsl_matrix_get(InvmmTm, i/(2*M+1), i%(2*M+1)));
	//if(i%(2*M+1)==(2*M)) printf("\n");
	//else printf("\t");
	//}  
	//printf("\n");	
	
	/* c_ij = U_-1 (m m^T)^-1 m */
	//printf("Filling parameters -  c \n");
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
									1.0, InvU, InvmmTm,
									0.0, c);
	
  //printf("c \n");
	//for(i=0;i<((c->size1)*(c->size2));i++)
	//{
	//printf("%e", gsl_matrix_get(c, i/(c->size2), i%(c->size2)));
	//if(i%(c->size2)==(c->size2-1)) printf("\n");
	//else printf("\t");
	//}  
	//printf("\n");	
	
	/* t_m^r = U*m */
	//printf("%dx%d -> (%dx%d)x(%dx%d)\n", tmr->size1, tmr->size2, U->size1, U->size2, m->size1, m->size2);
	//printf("Filling parameters -  tmr \n");
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
									1.0, U, m,
									0.0, tmr);	
	
	
	/*********************************************************************
   *
   * Set polynomial prefactors and smooth function
   * 
   ********************************************************************/ 
  
  //printf("Smoothing \n");
	/* READ DATA POINTS INTO VECTOR */
	for(j=0;j<2*M+1;j++)
	{
		gsl_vector_set(ym, j, y[j]);		
	}
	
	/* a = c*y */
	gsl_blas_dgemv( CblasNoTrans, 
								 1.0, c, ym, 
								 0.0, a );
	
	for(k=0; k<R; k++)
	{
		D[k] = factorial*gsl_vector_get(a, k);
	}
	
	gsl_vector_set_zero (ym);
	gsl_vector_set_zero (a);
	
  /*********************************************************************
   *
   * Output to file
   * 
   ********************************************************************/  
  //printf("Write to file \n");
  //FILE *smoothOut;
  //smoothOut = fopen("smoothOut.dat", "w");
  
  //for(k=0;k<N;k++)
  //{
	//fprintf(smoothOut, "%e\t%e\n", k*dt, Output[k]);
	//}
	//fclose(smoothOut);
	
  /*********************************************************************
   *
   * Clean Up
   *  
   ********************************************************************/  
	
	//printf("Cleaning up \n");
	gsl_matrix_free(m);
	gsl_matrix_free(U);
	gsl_vector_free(a);
	gsl_matrix_free(c);  
	gsl_vector_free(ym);
	gsl_matrix_free(tmr);
	gsl_matrix_free(mT);
	gsl_matrix_free(mmT);
	gsl_matrix_free(InvmmT);
	gsl_matrix_free(InvmmTm);
	gsl_matrix_free(InvU);
	
	//DETATCHSTATUSPTR(status);
	//RETURN(status);
}

void InvertMatrixSVD (gsl_matrix *A,gsl_matrix	*InvA,	int	N)
{ 
	
  /*********************************************************************
   *
   *  CREATING TEMPORARY VARIABLES
   * 
   ********************************************************************/  	
  int i=0;
  // Initialise matrices A, U, S and V
  gsl_matrix *InvS  = gsl_matrix_calloc (N, N); // inverse S
  gsl_matrix *V     = gsl_matrix_calloc (N, N); // V
  gsl_matrix *U     = gsl_matrix_calloc (N, N); // U
  gsl_matrix *C     = gsl_matrix_calloc (N, N); // temporary storage
  gsl_vector *s     = gsl_vector_alloc (N);     // eigenvalues AA^T
  gsl_matrix *II= gsl_matrix_calloc (N,N); // testing idenity
	
  //printf("INPUT \n");
	//for(i=0;i<(N*N);i++)
	//{
	//printf("%e", gsl_matrix_get(A, i/N, i%N));
	//if(i%N==(N-1)) printf("\n");
	//else printf("\t");
	//}  
	//printf("\n");
  
  /*********************************************************************
   *
   *  COMPUTING INVERSE
   * 		- PERFORM SVD
   * 		- CALCULATE INVERSE
   * 
   ********************************************************************/ 
	
	// Prepare U for SVD
	gsl_matrix_memcpy(U, A);
	
	// Perform SVD
	gsl_linalg_SV_decomp_jacobi(U, V, s);  
	
	// Compute Inverse S
	for (i = 0; i<N; i++)
	{
		gsl_vector_set( s, i, 1./gsl_vector_get( s, i) );
		gsl_matrix_set( InvS, i, i, gsl_vector_get( s, i) );
	}
	
	//printf("EIGENVECTORS \n");
	//for(i=0;i<N;i++)
	//{
	//printf("%e", gsl_vector_get(s,i));
	//if(i==(N-1)) printf("\n");
	//else printf("\t");
	//}
	//printf("\n");
	
	// Tranpose U
	gsl_matrix_transpose(U);
	
	// Multiply V and InvS
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
									1.0, V, InvS,
									0.0, C);
	
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
									1.0, C, U,
									0.0, InvA);                             
  
  //printf("INVERSE \n");
	//for(i=0;i<(N*N);i++)
	//{
	//printf("%e", gsl_matrix_get(InvA, i/N, i%N));
	//if(i%N==(N-1)) printf("\n");
	//else printf("\t");
	//}  
	//printf("\n");  
  
  /*********************************************************************
   *
   *  TESTING ACCURACY
   * 		- A * INVA = 1
   * 
   ********************************************************************/  
  
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
									1.0, A, InvA,
									0.0, II);      
	
	//printf("UNIT\n");
	//for(i=0;i<(N*N);i++)
	//{
	//printf("%e", gsl_matrix_get(I, i/N, i%N));
	//if(i%N==(N-1)) printf("\n");
	//else printf("\t");
	//}
	//printf("\n");
  
  /*********************************************************************
   *
   *  CLEANING UP
   * 
   ********************************************************************/  
	
	/* MATRICES */
  //gsl_matrix_free(A);
  gsl_matrix_free(U);
  gsl_matrix_free(InvS);
  //gsl_matrix_free(InvA);
  gsl_matrix_free(V);
  gsl_matrix_free(C);
  gsl_matrix_free(II);
  gsl_vector_free(s);
  
  /*********************************************************************
   *
   *  Detach Error handling
   * 
   ********************************************************************/
  
  return;
}

REAL8 ConvertCoefficientsToFunction(REAL8 *coeff, REAL8 f)
{
  /* Takes a set of N coefficients c_i and return 
   * c[0]+ c[1]*logf ..+c[N-1]*logf^(N-1).
   * To simulate constant CE, call with only c[0] non zero */ 

  REAL8 output = 0.0;
  int i;
  /* Space frequencies uniformly in logF */
  REAL8 logFreqs[Npoints];
  REAL8 deltalogf=(log10(freq_max)-log10(freq_min))/(REAL8)(Npoints-1);
  for (i=0; i<Npoints; i++) {
    logFreqs[i]=log10(freq_min)+deltalogf*i;
  } 
  if (Npoints%2==0) fprintf(stderr,"The number of points used for the fit must be odd, in ConvertCoefficientsToFunction. Exiting\n");

  REAL8 cen = logFreqs[(Npoints-1)/2]; 
  REAL8 logF = log10(f)-cen; // FIT USED CEN AS CENTRAL POINT!	

  for(i=0;i<FitOrder;i++)
  {
    output += coeff[i]*pow(logF, (REAL8) i);
  }

  return output;
}

static REAL8 ConvertRandTransitionSlopeToFunction(REAL8 *coeff,REAL8 f)
{
  /* Takes an array of 3 numbers and return a calibraation error realization which is:
   * 
   * flat (coeff[0]) + fluctuation (given by ~10% of the S6 1-sigma) for freq< coeff[1]
   * increasing with constant slope coeff[2] for freq>coeff[1]
   *  */
  REAL8 output=coeff[Npoints];
  REAL8 heavi=coeff[Npoints+1];
  REAL8 slope=coeff[Npoints+2];
  if (f<heavi)
    output+=(ConvertCoefficientsToFunction(coeff,f)/10.);
  if (f>= heavi) 
    output = slope*(f-heavi)+(output+(ConvertCoefficientsToFunction(coeff,heavi)/10.));
	return output;
}

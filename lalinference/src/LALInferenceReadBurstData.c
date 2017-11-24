/* 
 *  LALInferenceReadData.c:  Bayesian Followup functions
 *
 *  Copyright (C) 2013 Salvatore Vitale
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

#include <stdio.h>
#include <stdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/LALCache.h>
#include <lal/LALFrStream.h>
#include <lal/TimeFreqFFT.h>
#include <lal/LALDetectors.h>
#include <lal/AVFactories.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/Date.h>
#include <lal/StringInput.h>
#include <lal/VectorOps.h>
#include <lal/Random.h>
#include <lal/LALNoiseModels.h>
#include <lal/XLALError.h>
#include <lal/GenerateInspiral.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/SeqFactories.h>
#include <lal/DetectorSite.h>
#include <lal/GenerateInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/Inject.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOMetadataInspiralUtils.h>
#include <lal/LIGOMetadataRingdownUtils.h>
#include <lal/LALInspiralBank.h>
#include <lal/FindChirp.h>
#include <lal/LALInspiralBank.h>
#include <lal/GenerateInspiral.h>
#include <lal/NRWaveInject.h>
#include <lal/GenerateInspRing.h>
#include <math.h>
#include <lal/LALInspiral.h>
#include <lal/LALSimulation.h>
#include <lal/LALInference.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceTemplate.h>
#include <lal/LIGOLwXMLBurstRead.h>
#include <lal/GenerateBurst.h>
#include <lal/LALInferenceBurstRoutines.h>
#include <lal/LALInferenceReadBurstData.h>
#include <lal/LALSimNoise.h>

#define LALINFERENCE_DEFAULT_FLOW "40.0"
//typedef void (NoiseFunc)(LALStatus *statusPtr,REAL8 *psd,REAL8 f);
static void PrintSNRsToFile(LALInferenceIFOData *IFOdata , char SNRpath[] );
void InjectBurstFD(LALInferenceIFOData *IFOdata, SimBurst *inj_table, ProcessParamsTable *commandLine);
//typedef void (NoiseFunc)(LALStatus *statusPtr,REAL8 *psd,REAL8 f);

void LALInferenceInjectBurstSignal(LALInferenceIFOData *IFOdata, ProcessParamsTable *commandLine)
{
	LALStatus status;
	memset(&status,0,sizeof(status));
	SimBurst *injTable=NULL;
  SimBurst *injEvent=NULL;
  TimeSlide *tslide=NULL;
	INT4 Ninj=0;
	INT4 event=0;
	UINT4 i=0,j=0;
	int si=0;
  LIGOTimeGPS injstart;
	REAL8 SNR=0.0,NetworkSNR=0.0;// previous_snr=0.0; 
	memset(&injstart,0,sizeof(LIGOTimeGPS));
	COMPLEX16FrequencySeries *injF=NULL;
 // FILE *rawWaveform=NULL;
	ProcessParamsTable *ppt=NULL;
	REAL8 bufferLength = 2048.0; /* Default length of buffer for injections (seconds) */
	LIGOTimeGPS bufferStart;

	LALInferenceIFOData *thisData=IFOdata->next;
	REAL8 minFlow=IFOdata->fLow;
	REAL8 MindeltaT=IFOdata->timeData->deltaT;
    char SNRpath[FILENAME_MAX+10]="";
	while(thisData){
          minFlow   = minFlow>thisData->fLow ? thisData->fLow : minFlow;
          MindeltaT = MindeltaT>thisData->timeData->deltaT ? thisData->timeData->deltaT : MindeltaT;
          thisData  = thisData->next;
	}
  
	thisData=IFOdata;
	
	if(!LALInferenceGetProcParamVal(commandLine,"--binj")) {fprintf(stdout,"No injection file specified, not injecting\n"); return;}
	if(LALInferenceGetProcParamVal(commandLine,"--event")){
	    event= atoi(LALInferenceGetProcParamVal(commandLine,"--event")->value);
	    fprintf(stdout,"Injecting event %d\n",event);
	}
	else
	    fprintf(stdout,"WARNING: you did not give --event. Injecting event 0 of the xml table, which may not be what you want!\n");

  ppt = LALInferenceGetProcParamVal(commandLine,"--outfile");
	if (ppt)
	    snprintf(SNRpath, sizeof(SNRpath), "%s_snr.txt", ppt->value);
	else
		snprintf(SNRpath, sizeof(SNRpath), "snr.txt");

	injTable=XLALSimBurstTableFromLIGOLw(LALInferenceGetProcParamVal(commandLine,"--binj")->value,0,0);
	REPORTSTATUS(&status);
  Ninj=-1;
  while(injTable){Ninj++;injTable=injTable->next;}
	if(Ninj < event){ 
	    fprintf(stderr,"Error reading event %d from %s\n",event,LALInferenceGetProcParamVal(commandLine,"--binj")->value);
	    exit(1);
  }
	injTable=XLALSimBurstTableFromLIGOLw(LALInferenceGetProcParamVal(commandLine,"--binj")->value,0,0);
	while(si<event) {si++; injTable = injTable->next;} /* Select event */
	injEvent = injTable;
	injEvent->next = NULL;
  tslide=XLALTimeSlideTableFromLIGOLw(LALInferenceGetProcParamVal(commandLine,"--binj")->value);
	REPORTSTATUS(&status);
    
  /* If it is the case, inject burst in the FreqDomain */
  int FDinj=0;
  if (injTable)
    if(XLALSimBurstImplementedFDApproximants(XLALGetBurstApproximantFromString(injEvent->waveform))) FDinj=1;
    
    
	if (FDinj)
    {
         InjectBurstFD(thisData, injEvent, commandLine);
         return;
    }
    
	/* Begin loop over interferometers */
	while(thisData){
   
		memcpy(&bufferStart,&thisData->timeData->epoch,sizeof(LIGOTimeGPS));
		XLALGPSAdd(&bufferStart,(REAL8) thisData->timeData->data->length * thisData->timeData->deltaT);
		XLALGPSAdd(&bufferStart,-bufferLength);
		char series_name[320];
    sprintf(series_name,"%s:injection",thisData->name);
    REAL8TimeSeries *inj8Wave=(REAL8TimeSeries *)XLALCreateREAL8TimeSeries(series_name,
                                                                           &thisData->timeData->epoch,
                                                                           0.0,
                                                                           thisData->timeData->deltaT,
                                                                           &lalStrainUnit,
                                                                           thisData->timeData->data->length);
		if(!inj8Wave) XLAL_ERROR_VOID(XLAL_EFUNC);
		for(i=0;i<inj8Wave->data->length;i++) inj8Wave->data->data[i]=0.0;
    
    REAL8 Q, centre_frequency;
    Q=injEvent->q;
    centre_frequency=injEvent->frequency;
    /* Check that 2*width_gauss_envelope is inside frequency range */
    if ((centre_frequency + 3.0*centre_frequency/Q)>=  1.0/(2.0*thisData->timeData->deltaT)){
      XLALPrintWarning("WARNING: Your sample rate is too low to ensure a good analysis for a SG centered at f0=%lf and with Q=%lf. Consider increasing it to more than %lf. Exiting...\n",centre_frequency,Q,2.0*(centre_frequency + 3.0*centre_frequency/Q));
    }
    if ((centre_frequency -3.0*centre_frequency/Q)<=  thisData->fLow){
      XLALPrintWarning(
      "WARNING: The low frenquency tail of your SG centered at f0=%lf and with Q=%lf will lie below the low frequency cutoff. Whit your current settings and parameters the minimum f0 you can analyze without cuts is %lf.\n Continuing... \n",centre_frequency,Q,centre_frequency -3.0*centre_frequency/Q);
    }
    XLALBurstInjectSignals(inj8Wave,injEvent,tslide,NULL);
    XLALResampleREAL8TimeSeries(inj8Wave,thisData->timeData->deltaT);

    injF=(COMPLEX16FrequencySeries *)XLALCreateCOMPLEX16FrequencySeries("injF",
										&thisData->timeData->epoch,
										0.0,
										thisData->freqData->deltaF,
										&lalDimensionlessUnit,
										thisData->freqData->data->length);
    if(!injF) {
      XLALPrintError("Unable to allocate memory for injection buffer\n");
      XLAL_ERROR_VOID(XLAL_EFUNC);
    }

    /* Window the data */
    //REAL4 WinNorm = sqrt(thisData->window->sumofsquares/thisData->window->data->length);
    for(j=0;j<inj8Wave->data->length;j++){
      inj8Wave->data->data[j]*=thisData->window->data->data[j]; /* /WinNorm; */ /* Window normalisation applied only in freq domain */
    }

    XLALREAL8TimeFreqFFT(injF,inj8Wave,thisData->timeToFreqFFTPlan);

    if(thisData->oneSidedNoisePowerSpectrum){
        UINT4 upper=thisData->fHigh/injF->deltaF;
        for(SNR=0.0,j=thisData->fLow/injF->deltaF;j<upper;j++){
          SNR+=pow(creal(injF->data->data[j]),2.0)/thisData->oneSidedNoisePowerSpectrum->data->data[j];
          SNR+=pow(cimag(injF->data->data[j]),2.0)/thisData->oneSidedNoisePowerSpectrum->data->data[j];
        }
        SNR*=4.0*injF->deltaF;
    }
    thisData->SNR=sqrt(SNR);
    NetworkSNR+=SNR;

    /* Actually inject the waveform */
    for(j=0;j<inj8Wave->data->length;j++) thisData->timeData->data->data[j]+=inj8Wave->data->data[j];
    fprintf(stdout,"Injected SNR in detector %s = %.1f\n",thisData->name,thisData->SNR);
    char filename[320];
    sprintf(filename,"%s_timeInjection.dat",thisData->name);
    FILE* file=fopen(filename, "w");
    for(j=0;j<inj8Wave->data->length;j++){   
      fprintf(file, "%.6f\t%lg\n", XLALGPSGetREAL8(&thisData->timeData->epoch) + thisData->timeData->deltaT*j, inj8Wave->data->data[j]);
    }
    fclose(file);
    sprintf(filename,"%s_freqInjection.dat",thisData->name);
    file=fopen(filename, "w");
    /* NOTE: Here I (salvo) got rid of the division by WinNorm that was done in the CBC version of this routine. This is because burst signals are short and centered in the middle of the segment. Thus, the actual signal will be in the flat part of the Tuckey window. */
    for(j=0;j<injF->data->length;j++){   
      thisData->freqData->data->data[j]+=(injF->data->data[j]);
      fprintf(file, "%lg %lg \t %lg\n", thisData->freqData->deltaF*j, creal(injF->data->data[j]), cimag(injF->data->data[j]));    
    }
    fclose(file);

    XLALDestroyREAL8TimeSeries(inj8Wave);
    XLALDestroyCOMPLEX16FrequencySeries(injF);
    thisData=thisData->next;
    }

    ppt=LALInferenceGetProcParamVal(commandLine,"--dont-dump-extras");
    if (!ppt){
      PrintSNRsToFile(IFOdata , SNRpath);
    }
    
    NetworkSNR=sqrt(NetworkSNR);
    fprintf(stdout,"Network SNR of event %d = %.1e\n",event,NetworkSNR);
    thisData=IFOdata;

    return;
}

/** Fill the variables passed in vars with the parameters of the injection passed in event
    will over-write and destroy any existing parameters. Param vary type will be fixed */
void LALInferenceBurstInjectionToVariables(SimBurst *theEventTable, LALInferenceVariables *vars)
{
    if(!vars) {
	XLALPrintError("Encountered NULL variables pointer");
   	XLAL_ERROR_VOID(XLAL_EINVAL);
	}
    /* Destroy existing parameters */
    if(vars->head!=NULL) LALInferenceClearVariables(vars);
    REAL8 q = theEventTable->q;
    REAL8 psi = theEventTable->psi;
    REAL8 injGPSTime = XLALGPSGetREAL8(&(theEventTable->time_geocent_gps));
    REAL8 hrss = theEventTable->hrss;
    REAL8 loghrss=log(hrss);
    REAL8 f0 = theEventTable->frequency;
    REAL8 pol_angle = theEventTable->pol_ellipse_angle;
    REAL8 eccentricity = theEventTable->pol_ellipse_e;
    REAL8 duration=theEventTable->duration;
    REAL8 dec = theEventTable->dec;
    REAL8 ra = theEventTable->ra;
    
    LALInferenceAddVariable(vars, "quality", &q, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars, "frequency", &f0, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars, "time", &injGPSTime, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars, "hrss", &hrss, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars, "polarisation", &(psi), LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars, "declination", &dec, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars, "rightascension", &ra, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars, "loghrss", &loghrss, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars,"duration",&duration,LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars, "polar_angle", &pol_angle, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars, "polar_eccentricity", &eccentricity, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
}

static void PrintSNRsToFile(LALInferenceIFOData *IFOdata , char SNRpath[] ){
  REAL8 NetSNR=0.0;
  LALInferenceIFOData *thisData=IFOdata;
  int nIFO=0;

  while(thisData){
    thisData=thisData->next;
    nIFO++;
  }
  FILE * snrout = fopen(SNRpath,"w");
  if(!snrout){
    fprintf(stderr,"Unable to open the path %s for writing SNR files\n",SNRpath);
    fprintf(stderr,"Error code %i: %s\n",errno,strerror(errno));
    exit(errno);
  }
  thisData=IFOdata;
  while(thisData){
    fprintf(snrout,"%s:\t %4.2f\n",thisData->name,thisData->SNR);
    nIFO++;
    NetSNR+=(thisData->SNR*thisData->SNR);
    thisData=thisData->next;
  }
  if (nIFO>1){ 
    fprintf(snrout,"Network:\t");
    fprintf(snrout,"%4.2f\n",sqrt(NetSNR));
  }
  fclose(snrout);
}

void InjectBurstFD(LALInferenceIFOData *IFOdata, SimBurst *inj_table, ProcessParamsTable *commandLine)
///*-------------- Inject in Frequency domain -----------------*/
{
  /* Inject a gravitational wave into the data in the frequency domain */
  LALStatus status;
  memset(&status,0,sizeof(LALStatus));
  INT4 errnum;
  char SNRpath[FILENAME_MAX+16];
  ProcessParamsTable *ppt=NULL;
  ppt = NULL; 
  ppt = LALInferenceGetProcParamVal(commandLine,"--outfile");
  if (ppt)
    snprintf(SNRpath,sizeof(SNRpath), "%s_snr.txt", ppt->value);
  else
    snprintf(SNRpath,sizeof(SNRpath), "snr.txt");
  //REAL8 WinNorm = sqrt(IFOdata->window->sumofsquares/IFOdata->window->data->length);
  BurstApproximant approx = XLALGetBurstApproximantFromString(inj_table->waveform);
  
  if( (int) approx == XLAL_FAILURE)
    ABORTXLAL(&status);

  REAL8 injtime=0.0;
  injtime=inj_table->time_geocent_gps.gpsSeconds + 1e-9*inj_table->time_geocent_gps.gpsNanoSeconds;
  REAL8 hrss_one=1.0;
  REAL8 deltaT = IFOdata->timeData->deltaT;
  REAL8 deltaF = IFOdata->freqData->deltaF;

  REAL8 f_min = IFOdata->fLow;
  REAL8 f_max = 0.0;

  LALSimBurstExtraParam *extraParams = NULL;

  /* Print a line with information about approximant, amp_order, phaseorder, tide order and spin order */
  fprintf(stdout,"\n\n---\t\t ---\n");
  fprintf(stdout,"Injection will run using Approximant %d (%s) in the frequency domain.\n",approx,XLALGetStringFromBurstApproximant(approx));
  fprintf(stdout,"---\t\t ---\n\n");

  COMPLEX16FrequencySeries *hptilde=NULL, *hctilde=NULL;

  XLALSimBurstChooseFDWaveform(&hptilde, &hctilde, deltaF,deltaT,inj_table->frequency,inj_table->q,inj_table->duration,f_min,f_max,hrss_one,inj_table->pol_ellipse_angle,inj_table->pol_ellipse_e,extraParams,approx);

  /* Fail if injection waveform generation was not successful */
  errnum = *XLALGetErrnoPtr();
  if (errnum != XLAL_SUCCESS) {
    XLALPrintError(" ERROR in InjectFD(): error encountered when injecting waveform. errnum=%d\n",errnum);
    exit(1);
  }
  LALInferenceIFOData *dataPtr;
  REAL8 Fplus, Fcross;
  REAL8 plainTemplateReal, plainTemplateImag;
  REAL8 templateReal, templateImag;
  LIGOTimeGPS GPSlal;
  REAL8 gmst;
  REAL8 chisquared;
  REAL8 timedelay;  /* time delay b/w iterferometer & geocenter w.r.t. sky location */
  REAL8 timeshift;  /* time shift (not necessarily same as above)                   */
  REAL8 twopit, re, im, dre, dim, newRe, newIm;
  UINT4 i, lower, upper;

  REAL8 temp=0.0;
  REAL8 NetSNR=0.0;

  /* figure out GMST: */
  XLALGPSSetREAL8(&GPSlal, injtime);
  gmst=XLALGreenwichMeanSiderealTime(&GPSlal);

  /* loop over data (different interferometers): */
  dataPtr = IFOdata;

  while (dataPtr != NULL) {
    /*-- WF to inject is now in hptilde and hctilde. --*/
    /* determine beam pattern response (Fplus and Fcross) for given Ifo: */
    XLALComputeDetAMResponse(&Fplus, &Fcross,
                                (const REAL4(*)[3])dataPtr->detector->response,
                                inj_table->ra, inj_table->dec,
                                inj_table->psi, gmst);

    /* signal arrival time (relative to geocenter); */
    timedelay = XLALTimeDelayFromEarthCenter(dataPtr->detector->location,
                                                inj_table->ra, inj_table->dec,
                                                &GPSlal);

    /* (negative timedelay means signal arrives earlier at Ifo than at geocenter, etc.) */
    /* amount by which to time-shift template (not necessarily same as above "timedelay"): */
    REAL8 instant = dataPtr->timeData->epoch.gpsSeconds + 1e-9*dataPtr->timeData->epoch.gpsNanoSeconds;

    timeshift = (injtime - instant) + timedelay;
    twopit    = LAL_TWOPI * (timeshift);
    /* Restore hrss (template has been calculated for hrss=1) effect in Fplus/Fcross: */
    Fplus*=inj_table->hrss;
    Fcross*=inj_table->hrss;
    dataPtr->fPlus = Fplus;
    dataPtr->fCross = Fcross;
    dataPtr->timeshift = timeshift;

    char InjFileName[320];
    sprintf(InjFileName,"injection_%s.dat",dataPtr->name);
    FILE *outInj=fopen(InjFileName,"w");

     /* determine frequency range & loop over frequency bins: */
    lower = (UINT4)ceil(dataPtr->fLow / deltaF);
    upper = (UINT4)floor(dataPtr->fHigh / deltaF);
    chisquared = 0.0;

    re = cos(twopit * deltaF * lower);
    im = -sin(twopit * deltaF * lower);
    for (i=lower; i<=upper; ++i){
      /* derive template (involving location/orientation parameters) from given plus/cross waveforms: */
      if (i < hptilde->data->length) {
          plainTemplateReal = Fplus * creal(hptilde->data->data[i])
                              +  Fcross * creal(hctilde->data->data[i]);
          plainTemplateImag = Fplus * cimag(hptilde->data->data[i])
                              +  Fcross * cimag(hctilde->data->data[i]);
      } else {
          plainTemplateReal = 0.0;
          plainTemplateImag = 0.0;
      }

      /* do time-shifting...             */
      /* (also un-do 1/deltaT scaling): */
      /* real & imag parts of  exp(-2*pi*i*f*deltaT): */
      templateReal = (plainTemplateReal*re - plainTemplateImag*im);
      templateImag = (plainTemplateReal*im + plainTemplateImag*re);

      /* Incremental values, using cos(theta) - 1 = -2*sin(theta/2)^2 */
      dim = -sin(twopit*deltaF);
      dre = -2.0*sin(0.5*twopit*deltaF)*sin(0.5*twopit*deltaF);
      newRe = re + re*dre - im * dim;
      newIm = im + re*dim + im*dre;
      re = newRe;
      im = newIm;

      fprintf(outInj,"%lf %e %e %e\n",i*deltaF ,templateReal,templateImag,1.0/dataPtr->oneSidedNoisePowerSpectrum->data->data[i]);
      dataPtr->freqData->data->data[i] += crect( templateReal, templateImag );

      temp = ((2.0/( deltaT*(double) dataPtr->timeData->data->length) * (templateReal*templateReal+templateImag*templateImag)) / dataPtr->oneSidedNoisePowerSpectrum->data->data[i]);
      chisquared  += temp;
    }
    printf("injected SNR %.1f in IFO %s\n",sqrt(2.0*chisquared),dataPtr->name);
    NetSNR+=2.0*chisquared;
    dataPtr->SNR=sqrt(2.0*chisquared);
    dataPtr = dataPtr->next;

    fclose(outInj);
  }
  printf("injected Network SNR %.1f \n",sqrt(NetSNR));

  ppt=LALInferenceGetProcParamVal(commandLine,"--dont-dump-extras");
  if (!ppt){
    PrintSNRsToFile(IFOdata , SNRpath);
  }

  XLALDestroyCOMPLEX16FrequencySeries(hctilde);
  XLALDestroyCOMPLEX16FrequencySeries(hptilde);
}

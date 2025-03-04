/**********************************************************************
*        Copyright (C) 2009 John Veitch, Stephen Fairhurst
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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
**********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <config.h>
#include <math.h>
#include <string.h>

#include <LALAppsVCSInfo.h>
#include <lal/LALNoiseModels.h>
#include <lal/LALStdio.h>
#include <lal/LALgetopt.h>
#include <lal/LALStdlib.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/Date.h>
#include <lal/FindChirp.h>
#include <lal/Units.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/InspiralInjectionParams.h>
#include <lal/VectorOps.h>
#include <lal/LALFrStream.h>
#include <lal/LALDetectors.h>
#include <lal/LALFrameIO.h>

#define PROGRAM_NAME "coinj"

#define USAGE \
"lalpps_coinj [options]\n\
--help                       display this message \n\
--input <injection.xml>      Specify input SimInspiralTable xml file\n\
--output-path OUTPUTPATH directory path where frames are going to be written to\n\
--skip-ascii-output           skip generation of  ASCII txt output file\n\
--response-type TYPE         TYPE of injection, [ strain | etmx | etmy ]\n\
--frames                     Create h(t) frame files\n\n\
[--maxSNR snrhigh --minSNR snrlow                     Adjust injections to have combined SNR between snrlow and snrhigh in the H1H2L1V1 network]\n\
[--SNR snr      adjust distance to get precisely this snr]\n\
[--GPSstart A --GPSend B     Only generate waveforms for injection between GPS seconds A and B (int)]\n\
[--max-chirp-dist DIST	     Set the maximum chirp distance in H, L or V to DIST]\n\
lalapps_coinj: create coherent injection files for LIGO and VIRGO\n"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

extern int vrbflg;

typedef enum
{
  noResponse,
  unityResponse,
  design,
  actuationX,
  actuationY
} ResponseType;

typedef struct actuationparameters
{
  REAL4        ETMXcal;
  REAL4        ETMYcal;
  REAL4        pendFX;
  REAL4        pendFY;
  REAL4        pendQX;
  REAL4        pendQY;
  REAL4        length;
} ActuationParametersType;

typedef void (NoiseFunc)(LALStatus *status,REAL8 *psd,REAL8 f);

double chirpDist(SimInspiralTable *inj, char ifo);
double chirpDist(SimInspiralTable *inj, char ifo){
  /* eff_dist(IFO)*(2.8*0.25^(3/5)/mchirp)^(5/6) */
  double eff_dist=0.0;
  if(!inj) {printf("Null pointer passed to chirpDist!\n"); exit(1);}
  switch (ifo)
  {
	case 'H': {eff_dist=inj->eff_dist_h; break;}
	case 'L': {eff_dist=inj->eff_dist_l; break;}
	case 'G': {eff_dist=inj->eff_dist_g; break;}
	case 'T': {eff_dist=inj->eff_dist_t; break;}
	case 'V': {eff_dist=inj->eff_dist_v; break;}
        default: {printf("ERROR: Unknown IFO code %c\n",ifo);}
  }
  return( eff_dist*pow(1.218/inj->mchirp , (5./6.)));
}

int main(int argc, char *argv[])
{
  LALStatus        XLAL_INIT_DECL(status);
  CHAR                inputfile[FILENAME_MAX];
  CHAR                outputpath[1000];
  CHAR                injtype[30];
  CHAR                det_name[10];
  INT8                detectorFlags;
  LIGOTimeGPS inj_epoch;
  REAL8                deltaT= 1.0/16384.0;
  REAL8                injLength=100.0; /* Ten seconds at end */
  REAL8                LeadupTime=95.0;
  REAL8                dynRange=1.0/3.0e-23;
  int                  result;

  UINT4                Nsamples,det_idx,i,inj_num=0;
  ActuationParametersType actuationParams[LAL_NUM_IFO];
  ActuationParametersType actData = {0,0,0,0,0,0,0};
  ResponseType injectionResponse=noResponse;
  FILE *                outfile;
  LIGOLwXMLStream                *xmlfp=NULL;
  CHAR                outfilename[FILENAME_MAX];
  CHAR                fname[FILENAME_MAX];

  const LALUnit strainPerCount={0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};
  const LALUnit countPerStrain={0,{0,0,0,0,0,-1,1},{0,0,0,0,0,0,0}};
  NoiseFunc *PSD = NULL;
  REAL8 PSDscale=1.0;
  int c;
  int repeatLoop=0;
  int hitTarget=0;
  SimInspiralTable *injTable=NULL;
  SimInspiralTable this_injection;
  SimInspiralTable *headTable=NULL;
  REAL4TimeSeries *TimeSeries;
  double max_chirp_dist=0.0;

  REAL4TimeSeries *actuationTimeSeries;
  COMPLEX8FrequencySeries *resp;
  COMPLEX8FrequencySeries *actuationResp;
  COMPLEX8FrequencySeries *transfer;
  COMPLEX8Vector *unity;

  FrOutPar        UNUSED VirgoOutPars;
  CHAR                VirgoParsSource[100];
  CHAR                VirgoParsInfo[100];
  REAL8                NetworkSNR=0.0;
  INT4  makeFrames=0;
  INT4 outputRaw=0;
  INT4 skipASCIIoutput=0;
  COMPLEX8FrequencySeries *fftData;
  REAL8 mySNRsq,mySNR;
  REAL4FFTPlan *fwd_plan;
  REAL8 minSNR=0.0,maxSNR=0.0;
  REAL8 UNUSED maxRatio=1.0, UNUSED minRatio=1.0;
  REAL8 targetSNR=0.0;
  INT4 GPSstart=0,GPSend=2147483647;
  int UNUSED SNROK=1;
  int rewriteXML=0;
  LALFrameH *frame;

  struct LALoption long_options[]=
    {
      {"help", no_argument, 0, 'h'},
      {"input",required_argument,0, 'i'},
      {"output-path",required_argument,0,'o'},
      {"skip-ascii-output",no_argument,&skipASCIIoutput,'A'},
      {"response-type",required_argument,0,'r'},
      {"frames",no_argument,&makeFrames,'F'},
      {"rawstrain",no_argument,&outputRaw,'s'},
      {"verbose",no_argument,&vrbflg,1},
      {"minSNR",required_argument,0,2},
      {"maxSNR",required_argument,0,3},
      {"SNR",required_argument,0,6},
      {"GPSstart",required_argument,0,4},
      {"GPSend",required_argument,0,5},
      {"max-chirp-dist",required_argument,0,'d'},
      {0,0,0,0}
    };

  /*taken from Calibration CVS file:
   * calibration/frequencydomain/runs/S5/H1/model/V3/H1DARMparams_849677446.m */
  actuationParams[LAL_IFO_H1].ETMXcal = -0.795e-9;
  actuationParams[LAL_IFO_H1].pendFX  = 0.767;
  actuationParams[LAL_IFO_H1].pendQX  = 10.0;
  actuationParams[LAL_IFO_H1].ETMYcal = -0.827e-9;
  actuationParams[LAL_IFO_H1].pendFY  = 0.761;
  actuationParams[LAL_IFO_H1].pendQY  = 10.0;
  actuationParams[LAL_IFO_H1].length  = 4000.0;

  /*taken from Calibration CVS file:
   * calibration/frequencydomain/runs/S5/H2/model/V3/H2DARMparams_849678155.m */
  actuationParams[LAL_IFO_H2].ETMXcal = -0.876e-9;
  actuationParams[LAL_IFO_H2].pendFX  = 0.749;
  actuationParams[LAL_IFO_H2].pendQX  = 10.0;
  actuationParams[LAL_IFO_H2].ETMYcal = -0.912e-9;
  actuationParams[LAL_IFO_H2].pendFY  = 0.764;
  actuationParams[LAL_IFO_H2].pendQY  = 10.0;
  actuationParams[LAL_IFO_H2].length  = 2000.0;

  /*taken from Calibration CVS file:
   * calibration/frequencydomain/runs/S5/L1/model/V3/L1DARMparams_841930071.m */
  actuationParams[LAL_IFO_L1].ETMXcal = -0.447e-9;
  actuationParams[LAL_IFO_L1].pendFX  = 0.766;
  actuationParams[LAL_IFO_L1].pendQX  = 100.0;
  actuationParams[LAL_IFO_L1].ETMYcal = -0.438e-9;
  actuationParams[LAL_IFO_L1].pendFY  = 0.756;
  actuationParams[LAL_IFO_L1].pendQY  = 100.0;
  actuationParams[LAL_IFO_L1].length  = 4000.0;


/*******************************************************************************/

/* Process input arguments */
  while(1)
    {
      int option_idx=0;
      c=LALgetopt_long_only(argc,argv,"hFi:d:",long_options,&option_idx);
      if(c==-1) break;
      switch(c)
        {
        case 'h':
          fprintf(stderr,USAGE);
          exit(0);
          break;
        case 'i':
          strncpy(inputfile,LALoptarg,FILENAME_MAX-1);
          break;
        case 'o':
          strncpy(outputpath,LALoptarg,999);
          break;
        case 'r':
          if(!strcmp("strain",LALoptarg))          injectionResponse = unityResponse;
          else if(!strcmp("etmx",LALoptarg))   injectionResponse = actuationX;
          else if(!strcmp("etmy",LALoptarg))   injectionResponse = actuationY;
          else {fprintf(stderr,"Invalid argument to response-type: %s\nResponse type must be strain, etmy or etmx\n", \
                        LALoptarg); exit(1);}
          break;
        case 2:
          minSNR=atof(LALoptarg);
          fprintf(stderr,"Using minimum SNR of %f\n",minSNR);
          break;
        case 3:
          maxSNR=atof(LALoptarg);
          fprintf(stderr,"Using maximum SNR of %f\n",maxSNR);
          break;
        case 4:
          GPSstart=atoi(LALoptarg);
          break;
        case 5:
          GPSend=atoi(LALoptarg);
          break;
        case 6:
          targetSNR=atof(LALoptarg);
          fprintf(stderr,"Target SNR = %lf\n",targetSNR);
          break;
        case 'd':
	  max_chirp_dist=atof(LALoptarg);
          fprintf(stderr,"Using maximum chirp distance of %lf\n",max_chirp_dist);
	  break;
	}
    }

  if(minSNR!=0 && maxSNR!=0 && (maxSNR<minSNR)){
    fprintf(stderr,"Error: minSNR must be less than maxSNR\n");
    exit(1);
    if(targetSNR!=0.0 && (targetSNR<minSNR || targetSNR>maxSNR)){
      fprintf(stderr,"Target SNR %lf is not in range %lf to %lf, ignoring min and max\n",targetSNR,minSNR,maxSNR);
    }
  }

  memset(&status,0,sizeof(status));

  /* Read in the input XML */
  injTable = XLALSimInspiralTableFromLIGOLw(inputfile);
  headTable=injTable;
  Nsamples = (UINT4)injLength/deltaT;

  do{
    memcpy(&this_injection,injTable,sizeof(SimInspiralTable));
    this_injection.next=NULL;
    NetworkSNR=0.0;
    /* Set epoch */
    memcpy(&inj_epoch,&(this_injection.geocent_end_time),sizeof(LIGOTimeGPS));
    inj_epoch = this_injection.geocent_end_time;
    XLALGPSAdd(&inj_epoch, -LeadupTime);
    inj_epoch.gpsNanoSeconds=0;
    SNROK=0; /* Reset this to 0 = OK */
    minRatio=2.0;
    maxRatio=0.0;
    repeatLoop=0;
    /* Loop over detectors */
    for(det_idx=0;det_idx<LAL_NUM_IFO;det_idx++){
      /* Only generate within chosen bounds, if specified */
      if((this_injection.geocent_end_time.gpsSeconds-(int)LeadupTime )<GPSstart || (this_injection.geocent_end_time.gpsSeconds-(int)LeadupTime)>GPSend) continue;

      if(det_idx==LAL_IFO_T1||det_idx==LAL_IFO_G1||det_idx==LAL_IFO_H2) continue; /* Don't generate for GEO or TAMA */

      switch(det_idx)
        {
        case LAL_IFO_H1: sprintf(det_name,"H1"); PSD=&LALLIGOIPsd; PSDscale=9E-46; detectorFlags = LAL_LHO_4K_DETECTOR_BIT; break;
        case LAL_IFO_H2: sprintf(det_name,"H2"); PSD=&LALLIGOIPsd; PSDscale=9E-46; detectorFlags = LAL_LHO_2K_DETECTOR_BIT; break;
        case LAL_IFO_L1: sprintf(det_name,"L1"); PSD=&LALLIGOIPsd; PSDscale=9E-46; detectorFlags = LAL_LLO_4K_DETECTOR_BIT; break;
        case LAL_IFO_V1: sprintf(det_name,"V1"); PSD=&LALVIRGOPsd; PSDscale=1.0; detectorFlags = LAL_VIRGO_DETECTOR_BIT;  break;
        case LAL_IFO_G1: sprintf(det_name,"G1"); PSD=&LALGEOPsd; PSDscale=1E-46; detectorFlags = LAL_GEO_600_DETECTOR_BIT;  break;
        case LAL_IFO_T1: sprintf(det_name,"T1"); PSD=&LALTAMAPsd; PSDscale=75E-46; detectorFlags = LAL_TAMA_300_DETECTOR_BIT; break;
        default: fprintf(stderr,"Unknown IFO\n"); exit(1); break;
        }

      TimeSeries=XLALCreateREAL4TimeSeries(det_name,&inj_epoch,0.0,deltaT,&lalADCCountUnit,(size_t)Nsamples);
      for(i=0;i<Nsamples;i++) TimeSeries->data->data[i]=0.0;
      resp = XLALCreateCOMPLEX8FrequencySeries("response",&inj_epoch,0.0,1.0/injLength,&strainPerCount,(size_t)Nsamples/2+1);
      for(i=0;i<resp->data->length;i++) {resp->data->data[i] = (REAL4)1.0/dynRange;}

      /* Create h(t) time series for this detector */
      LAL_CALL( LALFindChirpInjectSignals(&status,TimeSeries,&this_injection,resp) , &status);

      XLALDestroyCOMPLEX8FrequencySeries(resp);

      if(det_idx==LAL_IFO_V1 && injectionResponse!=unityResponse){
        actuationTimeSeries=NULL;
        goto calcSNRandwriteFrames;
      }

      /* -=-=-=-=-=-=- Prepare actuations -=-=-=-=-=-=- */

      if(injectionResponse==actuationX || injectionResponse==actuationY) actData=actuationParams[det_idx];
      actuationResp = XLALCreateCOMPLEX8FrequencySeries("actuationResponse",&inj_epoch,0.0,1.0/(2.0*injLength),&strainPerCount,(size_t)Nsamples/2+1);
      /* Create actuation response */
      switch(injectionResponse){
      case unityResponse:
        sprintf(injtype,"STRAIN");
        for(i=0;i<actuationResp->data->length;i++){actuationResp->data->data[i] = 1.0;}
        break;
      case actuationX:
        sprintf(injtype,"ETMX");
        actuationResp=generateActuation(actuationResp,actData.ETMXcal/actData.length,actData.pendFX,actData.pendQX);
        break;
      case actuationY:
        sprintf(injtype,"ETMY");
        actuationResp=generateActuation(actuationResp,actData.ETMYcal/actData.length,actData.pendFY,actData.pendQY);
        break;
      default:
        fprintf(stderr,"Must specify response function: strain, etmy or etmx\n"); exit(1);
        break;
      }



      if(injectionResponse!=unityResponse) {
        actuationTimeSeries=XLALCreateREAL4TimeSeries(det_name,&inj_epoch,0.0,deltaT,&lalADCCountUnit,(size_t)Nsamples);
        unity = XLALCreateCOMPLEX8Vector(actuationResp->data->length);
        transfer = XLALCreateCOMPLEX8FrequencySeries("transfer",&inj_epoch,0.0,1.0/(2.0*injLength),&countPerStrain,(size_t)Nsamples/2+1);
        for(i=0;i<unity->length;i++) {unity->data[i] = 1.0;}
        XLALCCVectorDivide(transfer->data,unity,actuationResp->data);
        for(i=0;i<Nsamples;i++) actuationTimeSeries->data->data[i]=TimeSeries->data->data[i];
        actuationTimeSeries = XLALRespFilt(actuationTimeSeries,transfer);
        XLALDestroyCOMPLEX8FrequencySeries(transfer);
        XLALDestroyCOMPLEX8Vector(unity);
        for(i=0;i<actuationTimeSeries->data->length;i++) actuationTimeSeries->data->data[i]/=dynRange;
      }
      else actuationTimeSeries=TimeSeries;

      XLALDestroyCOMPLEX8FrequencySeries(actuationResp);

      /* Output the actuation time series */
      /*sprintf(outfilename,"HWINJ_%i_%s_%i_%s.out",inj_num,injtype,inj_epoch.gpsSeconds,det_name);*/
      char massBin[5];
      if(this_injection.mass1<2.0 && this_injection.mass2<2.0) sprintf(massBin,"BNS");
      else if(this_injection.mass1<2.0 || this_injection.mass2<2.0) sprintf(massBin,"NSBH");
      else sprintf(massBin,"BBH");

      if (!(skipASCIIoutput)){
        result = snprintf(outfilename,sizeof(outfilename),"%s%s%i_CBC_%s_%i_%s_%s.txt",outputpath,"/",inj_epoch.gpsSeconds,massBin,inj_num,injtype,det_name);
        if (result < 0 || (size_t)result > sizeof(outfilename) - 1)
        {
          fprintf(stderr, "ERROR: file name too long: '%s'\n", outfilename);
          exit(1);
        }
        outfile=fopen(outfilename,"w");
        fprintf(stdout,"Injected signal %i for %s into file %s\n",inj_num,det_name,outfilename);
        for(i=0;i<actuationTimeSeries->data->length;i++) fprintf(outfile,"%10.10e\n",actuationTimeSeries->data->data[i]);
        fclose(outfile);
      };

    calcSNRandwriteFrames:

      /* Calculate SNR for this injection */
      fwd_plan = XLALCreateForwardREAL4FFTPlan( TimeSeries->data->length, 0 );
      fftData = XLALCreateCOMPLEX8FrequencySeries(TimeSeries->name,&(TimeSeries->epoch),0,1.0/TimeSeries->deltaT,&lalDimensionlessUnit,TimeSeries->data->length/2 +1);
      XLALREAL4TimeFreqFFT(fftData,TimeSeries,fwd_plan);
      XLALDestroyREAL4FFTPlan(fwd_plan);

      mySNRsq = 0.0;
      mySNR=0.0;
      for(i=1;i<fftData->data->length;i++){
        REAL8 freq;
        REAL8 sim_psd_value=0;
        freq = fftData->deltaF * i;
        PSD( &status, &sim_psd_value, freq );
        mySNRsq += crealf(fftData->data->data[i]) * crealf(fftData->data->data[i]) /
          (sim_psd_value*PSDscale);
        mySNRsq += cimagf(fftData->data->data[i]) * cimagf(fftData->data->data[i]) /
          (sim_psd_value*PSDscale);
      }
      mySNRsq *= 4.0*fftData->deltaF;
      XLALDestroyCOMPLEX8FrequencySeries( fftData );
      if(det_idx==LAL_IFO_H2) mySNRsq/=4.0;
      mySNR = sqrt(mySNRsq)/dynRange;
      fprintf(stdout,"SNR in design %s of injection %i = %lf\n",det_name,inj_num,mySNR);

      for(i=0;i<TimeSeries->data->length;i++) {
        TimeSeries->data->data[i]=TimeSeries->data->data[i]/dynRange +0.0;
      }
      NetworkSNR+=mySNR*mySNR;

      if(makeFrames){ /* Also output frames for Virgo */
        sprintf(VirgoParsSource,"%s-INSP%i",det_name,inj_num);
        VirgoOutPars.source=VirgoParsSource;
        sprintf(VirgoParsInfo,"HWINJ-STRAIN");
        VirgoOutPars.description=VirgoParsInfo;
        VirgoOutPars.type=ProcDataChannel;
        VirgoOutPars.nframes=(UINT4)injLength;
        VirgoOutPars.frame=0;
        VirgoOutPars.run=2;
        fprintf(stdout,"Generating frame file for %s-%s-%i\n",VirgoParsSource,VirgoParsInfo,TimeSeries->epoch.gpsSeconds);
       /* LALFrWriteREAL4TimeSeries(&status,TimeSeries,&VirgoOutPars); */
        /* define frame */
        frame = XLALFrameNew( &TimeSeries->epoch, injLength, "HWINJ-STRAIN", 0, 1, detectorFlags );
        /* add time series as a channel to the frame */
        XLALFrameAddREAL4TimeSeriesSimData( frame, TimeSeries );
        /* write frame */
        result = snprintf(fname,sizeof(fname),"%s%s%s-INSP%i_HWINJ_STRAIN-%i-%i.gwf",outputpath,"/",det_name, inj_num, inj_epoch.gpsSeconds, (UINT4)injLength);
        if (result < 0 || (size_t)result > sizeof(fname) - 1)
        {
          fprintf( stderr, "ERROR: file name too long: '%s'\n", fname );
          exit( 1 );
        }
        /*sprintf(fname, "%s%s%s",outputpath, "/", fname);*/
        if (XLALFrameWrite( frame, fname) != 0)
        {
          fprintf( stderr, "ERROR: Cannot save frame file: '%s'\n", fname );
          exit( 1 );
        }

        /* clear frame */
        XLALFrameFree( frame );


      }

      if(TimeSeries==actuationTimeSeries) XLALDestroyREAL4TimeSeries(TimeSeries);
      else {
        if(injectionResponse) XLALDestroyREAL4TimeSeries(actuationTimeSeries);
        XLALDestroyREAL4TimeSeries(TimeSeries);
      }
} /* End loop over detectors */

    /*fprintf(stdout,"Finished injecting signal %i, network SNR %f\n",inj_num,sqrt(NetworkSNR));*/
    NetworkSNR=sqrt(NetworkSNR);
    if(NetworkSNR!=0.0){ /* Check the SNR if we did this injection */
      if(targetSNR!=0.0 && repeatLoop==0 && hitTarget==0) {
        injTable->distance*=(REAL4)(NetworkSNR/targetSNR);
        injTable->eff_dist_h*=(REAL4)(NetworkSNR/targetSNR);
        injTable->eff_dist_l*=(REAL4)(NetworkSNR/targetSNR);
        injTable->eff_dist_v*=(REAL4)(NetworkSNR/targetSNR);
        injTable->eff_dist_g*=(REAL4)(NetworkSNR/targetSNR);
        injTable->eff_dist_t*=(REAL4)(NetworkSNR/targetSNR);
        rewriteXML=1; repeatLoop=1; hitTarget=1;}
      else {repeatLoop=0; hitTarget=0;}
      if(targetSNR==0.0 && minSNR>NetworkSNR) {
        injTable->distance*=(REAL4)(0.99*NetworkSNR/minSNR);
        injTable->eff_dist_h*=(REAL4)(0.99*NetworkSNR/minSNR);
        injTable->eff_dist_l*=(REAL4)(0.99*NetworkSNR/minSNR);
        injTable->eff_dist_v*=(REAL4)(0.99*NetworkSNR/minSNR);
        injTable->eff_dist_t*=(REAL4)(0.99*NetworkSNR/minSNR);
        injTable->eff_dist_g*=(REAL4)(0.99*NetworkSNR/minSNR);
        rewriteXML=1; repeatLoop=1;}
      else {
        if(targetSNR==0.0 && maxSNR!=0.0 && maxSNR<NetworkSNR) {
          injTable->distance*=(1.01*NetworkSNR/maxSNR);
          injTable->eff_dist_h*=(1.01*NetworkSNR/maxSNR);
          injTable->eff_dist_l*=(1.01*NetworkSNR/maxSNR);
          injTable->eff_dist_v*=(1.01*NetworkSNR/maxSNR);
          injTable->eff_dist_t*=(1.01*NetworkSNR/maxSNR);
          injTable->eff_dist_g*=(1.01*NetworkSNR/maxSNR);
          /*injTable->distance+=0.01;*/
          rewriteXML=1;
          repeatLoop=1;
          fprintf(stderr,"Multiplying by %lf to get from %lf to target\n",1.01*(NetworkSNR/maxSNR),NetworkSNR);}
      }
    }
    if(max_chirp_dist!=0.0 && (maxSNR==0 || (maxSNR!=0 && (maxSNR>NetworkSNR) ) )  ){
	double this_max=0.0;
	double second_max=0.0;
	char ifostr[]="HLV";
	for(int cidx=0;cidx<3;cidx++)
	{
		if(this_max<chirpDist(injTable,ifostr[cidx])) this_max=chirpDist(injTable,ifostr[cidx]);
	}
	for(int cidx=0;cidx<3;cidx++)
		if(second_max<chirpDist(injTable,ifostr[cidx]) && chirpDist(injTable,ifostr[cidx])<this_max) second_max=chirpDist(injTable,ifostr[cidx]);
	if(second_max>max_chirp_dist){
	injTable->distance*=0.95*(max_chirp_dist/second_max);
	injTable->eff_dist_h*=0.95*max_chirp_dist/second_max;
	injTable->eff_dist_l*=0.95*max_chirp_dist/second_max;
	injTable->eff_dist_v*=0.95*max_chirp_dist/second_max;
	injTable->eff_dist_t*=0.95*max_chirp_dist/second_max;
	injTable->eff_dist_g*=0.95*max_chirp_dist/second_max;
	rewriteXML=1; repeatLoop=1;
	fprintf(stderr,"MCD: Multiplying distance by %lf to get from %lf to target\n",0.95*max_chirp_dist/second_max,second_max);
    	}
    }
    if(repeatLoop==1) fprintf(stderr,"Reinjecting with new distance %f for desired SNR\n\n",injTable->distance);

    if(repeatLoop==0){
      fprintf(stderr,"\nNetwork SNR of %i = %lf\n",inj_num,NetworkSNR);
      injTable=injTable->next;
      inj_num++;
    }


  }while(injTable!=NULL);

  /* If the distances were adjusted, re-write the SimInspiral table */
  if(rewriteXML){
    fprintf(stderr,"Overwriting %s with adjusted distances\n",inputfile);
    inputfile[strlen(inputfile)-4] = 0; /* Cut off the .xml */
    strncat(inputfile,"_adj.xml",FILENAME_MAX-1);
    xmlfp=XLALOpenLIGOLwXMLFile(inputfile);
    if(xmlfp==NULL) fprintf(stderr,"Error! Cannot open %s for writing\n",inputfile);
    XLALWriteLIGOLwXMLSimInspiralTable( xmlfp, headTable );
    XLALCloseLIGOLwXMLFile( xmlfp );
  }

  return(0);
}

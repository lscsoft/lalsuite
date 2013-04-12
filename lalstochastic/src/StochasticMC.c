/*
*  Copyright (C) 2007 Jolien Creighton, Kaice T. Reilly, Robert Adam Mercer, Tania Regimbau
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


/**
\author Tania Regimbau, Sukanta Bose, Jeff Noel
\addtogroup StochasticMC_h

\heading{Description}

This routine simulates time-domain signal in a pair
of detectors using Sukanta Bose's code SimulateSB.c, whitened with the adequate response function that can be used in LALwrapper.

In this version, long time-series are constructed by concatenating short segments of simulated data
whitened with the adequate response function. Segments are sinusoidally spliced into consecutive
segments using Jeff Noel's function SinusoidalSplice, in order to avoid discontinuities in the final
time serie.

\heading{Algorithm}

The following program shows how to use the routines LALStochasticMCDso and LALStochasticMCDsoSplice

\code
#include <math.h>
#include <string.h>
#include <stdio.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif
#include <FrameL.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/PrintFTSeries.h>
#include <lal/FrameStream.h>
#include <lal/FrameCalibration.h>
#include <lal/Calibration.h>
#include <lal/LALConstants.h>
#include <lal/LALStatusMacros.h>
#include <lal/StochasticCrossCorrelation.h>
#include <lal/RealFFT.h>
#include <lal/ComplexFFT.h>
#include <lal/Units.h>
#include <lal/StreamInput.h>
#include <lal/PrintVector.h>
#include <lal/Random.h>
#include <lal/SimulateSB.h>
#include <lal/StochasticMC.h>

int main( ){

  static LALStatus status;
  SSSimStochBGOutput MCoutput;
  StochasticMCParams MCparams;
  StochasticMCSInput MCinput;

  //output structure
  REAL4TimeSeries SimStochBG1;
  REAL4TimeSeries SimStochBG2;

  //input structure
  CalibrationFunctions calfuncs1,calfuncs2;
  CalibrationUpdateParams  calfacts1,calfacts2;
  COMPLEX8FrequencySeries responseFunction1,responseFunction2,sensingFunction1,sensingFunction2;
  LIGOTimeGPS epoch;
  COMPLEX8TimeSeries oloopfactor1,sfactor1,oloopfactor2,sfactor2;
  COMPLEX8Vector *response[2]={NULL,NULL};
  COMPLEX8Vector *sensing[2]={NULL,NULL};
  COMPLEX8Vector *oloopfactor[2]={NULL,NULL};
  COMPLEX8Vector *sfactor[2]={NULL,NULL};

  //parameters
  UINT4 timeref, starttime, caltime;
  UINT4 freqlen,length,lengthseg,caloffset;
  REAL8 sRate,deltaT, deltaF;

  INT4 i, j;
  LALUnit countPerStrain = {0,{0,0,0,0,0,-1,1},{0,0,0,0,0,0,0}};

  CHAR Outfile1[LALNameLength];
  CHAR Outfile2[LALNameLength];
  FILE *pfone, *pftwo;

  status.statusPtr = NULL;
  lengthseg = 61440;
  sRate = 1024.;
  starttime = 729331170;

  //write parameters
  MCparams.lengthseg = lengthseg; MCparams.numseg = 3;
  MCparams.sRate = sRate;
  MCparams.starttime = starttime;
  MCparams.seed = 173;
  MCparams.fRef = 100.; MCparams.f0 = 0.;
  MCparams.omegaRef = 1.; MCparams.alpha = 0;
  MCparams.site1 = 0; MCparams.site2 = 0;

  //set other parameters
  length = MCparams.numseg * lengthseg;
  deltaT = 1./sRate;
  freqlen = lengthseg / 2 + 1;
  caloffset =  lengthseg / (2 * sRate);
  deltaF = 1.0/(deltaT*lengthseg);
  timeref = 0;
  caltime = starttime + caloffset;
  for (i=0;i<2;i++)
  {
    LALCCreateVector(&status,&response[i],freqlen);
    LALCCreateVector(&status,&sensing[i],freqlen);
  }

 for (i=0;i<2;i++)
  {
    LALCCreateVector(&status,&oloopfactor[i],length);
    LALCCreateVector(&status,&sfactor[i],length);
  }

 for (i=0;i<2;i++)
  {
    for (j=0; j < freqlen; j++)
     {
      response[i]->data[j].re = 1.;
      response[i]->data[j].im = 1.;
      sensing[i]->data[j].re = 1.;
      sensing[i]->data[j].im = 1.;
     }
  }

  for (i=0;i<2;i++)
  {
    for (j=0; j < length; j++)
     {
      oloopfactor[i]->data[j].re = 1.;
      oloopfactor[i]->data[j].im = 1.;
      sfactor[i]->data[j].re = 1.;
      sfactor[i]->data[j].im = 1.;
     }
  }

  memset(&responseFunction1,0,sizeof(COMPLEX8FrequencySeries));
  memset(&responseFunction2,0,sizeof(COMPLEX8FrequencySeries));
  memset(&sensingFunction1,0,sizeof(COMPLEX8FrequencySeries));
  memset(&sensingFunction2,0,sizeof(COMPLEX8FrequencySeries));

  responseFunction1.epoch.gpsSeconds = timeref;
  responseFunction1.epoch.gpsNanoSeconds = 0;
  responseFunction1.deltaF = deltaF;
  responseFunction1.f0 = 0;
  responseFunction1.sampleUnits = countPerStrain;
  responseFunction1.data = response[0];

  responseFunction2.epoch.gpsSeconds = timeref;
  responseFunction2.epoch.gpsNanoSeconds = 0;
  responseFunction2.deltaF = deltaF;
  responseFunction2.f0 = 0;
  responseFunction2.sampleUnits = countPerStrain;
  responseFunction2.data = response[1];

  sensingFunction1.epoch.gpsSeconds = timeref;
  sensingFunction1.epoch.gpsNanoSeconds = 0;
  sensingFunction1.deltaF = deltaF;
  sensingFunction1.f0 = 0;
  sensingFunction1.sampleUnits = countPerStrain;
  sensingFunction1.data = sensing[0];

  sensingFunction2.epoch.gpsSeconds = timeref;
  sensingFunction2.epoch.gpsNanoSeconds = 0;
  sensingFunction2.deltaF = deltaF;
  sensingFunction2.f0 = 0;
  sensingFunction2.sampleUnits = countPerStrain;
  sensingFunction2.data = sensing[1];

  memset(&oloopfactor1,0,sizeof(COMPLEX8TimeSeries));
  memset(&oloopfactor2,0,sizeof(COMPLEX8TimeSeries));
  memset(&sfactor1,0,sizeof(COMPLEX8TimeSeries));
  memset(&sfactor2,0,sizeof(COMPLEX8TimeSeries));

  oloopfactor1.epoch.gpsSeconds = timeref;
  oloopfactor1.epoch.gpsNanoSeconds = 0;
  oloopfactor1.deltaT = 1. / sRate;
  oloopfactor1.f0 = 0;
  oloopfactor1.sampleUnits = lalADCCountUnit;
  oloopfactor1.data = oloopfactor[0];

  oloopfactor2.epoch.gpsSeconds = timeref;
  oloopfactor2.epoch.gpsNanoSeconds = 0;
  oloopfactor2.deltaT = 1. / sRate;
  oloopfactor2.f0 = 0;
  oloopfactor2.sampleUnits = lalADCCountUnit;
  oloopfactor2.data = oloopfactor[1];

  sfactor1.epoch.gpsSeconds = timeref;
  sfactor1.epoch.gpsNanoSeconds = 0;
  sfactor1.deltaT = 1. / sRate;
  sfactor1.f0 = 0;
  sfactor1.sampleUnits = lalADCCountUnit;
  sfactor1.data = sfactor[0];

  sfactor2.epoch.gpsSeconds = timeref;
  sfactor2.epoch.gpsNanoSeconds = 0;
  sfactor2.deltaT = 1. / sRate;
  sfactor2.f0 = 0;
  sfactor2.sampleUnits = lalADCCountUnit;
  sfactor2.data = sfactor[1];

  calfuncs1.responseFunction = &responseFunction1;
  calfuncs1.sensingFunction = &sensingFunction1;
  calfuncs2.responseFunction = &responseFunction2;
  calfuncs2.sensingFunction = &sensingFunction2;

  calfacts1.openLoopFactor = &oloopfactor1;
  calfacts1.sensingFactor = &sfactor1;
  calfacts2.openLoopFactor = &oloopfactor2;
  calfacts2.sensingFactor = &sfactor2;

  epoch.gpsSeconds = caltime;
  epoch.gpsNanoSeconds = 0;
  calfacts1.epoch = epoch;calfacts2.epoch = epoch;

  //write input parameters
  MCinput.calfuncs1 = calfuncs1;
  MCinput.calfuncs2 = calfuncs2;
  MCinput.calfacts1 = calfacts1;
  MCinput.calfacts2 = calfacts2;

  SimStochBG1.data = NULL;
  LALSCreateVector(&status, &(SimStochBG1.data), length);
  SimStochBG2.data = NULL;
  LALSCreateVector(&status, &(SimStochBG2.data), length);

  MCoutput.SSimStochBG1 = &SimStochBG1;
  MCoutput.SSimStochBG2 = &SimStochBG2;

  LALStochasticMCDso (&status,&MCoutput,&MCinput,&MCparams);

  //create output files

  if (print){
  snprintf( Outfile1, LALNameLength * sizeof(CHAR),
            "%s_1_%d.ilwd",_IFO1,starttime);
  snprintf( Outfile2, LALNameLength * sizeof(CHAR),
            "%s_2_%d.ilwd",_IFO2,starttime);

  //print output to ilwds

   if (print)
    {
     pfone=LALFopen(Outfile1,"w");pftwo=LALFopen(Outfile2,"w");
     fprintf(pfone,"<?ilwd?>\n");fprintf(pftwo,"<?ilwd?>\n");
     fprintf(pfone,"<ilwd comment='%s'",Outfile1);
     fprintf(pfone," name='%s' size='1'>\n",Outfile1);
     fprintf(pfone," <real_8 dims='%d' name='%s'>",length,Outfile1);
     fprintf(pftwo,"<ilwd comment='%s'",Outfile2);
     fprintf(pftwo," name='%s' size='1'>\n",Outfile2);
     fprintf(pftwo," <real_8 dims='%d' name='%s'>",length,Outfile2);

     for(i=0;(UINT4)i<length;i++)
       {
        fprintf(pfone,"% e",SimStochBG1.data->data[i]);
        fprintf(pftwo,"% e",SimStochBG2.data->data[i]);
       }
     fprintf(pfone,"</real_8>");fprintf(pftwo,"</real_8>");
     fprintf(pfone," </ilwd>");fprintf(pftwo," </ilwd>");
     LALFclose(pfone);LALFclose(pftwo);
    }
  }
}
\endcode

\heading{Uses}
\code
LALSSSimStochBGTimeSeries()
LALUpdateCalibration()
LALResponseConvert()
\endcode
*/


#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <math.h>
#include <string.h>
#include <stdio.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/Calibration.h>
#include <lal/LALConstants.h>
#include <lal/LALStatusMacros.h>
#include <lal/StochasticCrossCorrelation.h>
#include <lal/RealFFT.h>
#include <lal/ComplexFFT.h>
#include <lal/Units.h>
#include <lal/Random.h>
#include <lal/SimulateSB.h>
#include <lal/StochasticMC.h>

static int verbose = 0;

static void SinusoidalSplice(REAL4Vector **longData, REAL4Vector **shortData, REAL4Vector *output, UINT4 nSpliceSegs, UINT4 offset);

void LALStochasticMCDso (LALStatus *status,
		         SSSimStochBGOutput *MCoutput,
	  	         StochasticMCInput  *MCinput,
		         StochasticMCParams *MCparams)
{

  /* This stores the parameters of functions used */
  StochasticOmegaGWParameters        parametersOmega;

  /* This stores the structures of SimulateSB */
  SSSimStochBGParams                 SBParams;
  SSSimStochBGInput                  SBInput;
  SSSimStochBGOutput                 SBOutput;
  REAL4TimeSeries                    whitenedSSimStochBG1;
  REAL4TimeSeries                    whitenedSSimStochBG2;

  REAL4FrequencySeries               omegaGW;
  COMPLEX8FrequencySeries            wFilter1;
  COMPLEX8FrequencySeries            wFilter2;


  /* vector to store response functions of a pair of detectors  */
  COMPLEX8Vector                    *response[2]={NULL,NULL};

  /* Counters & output */
  INT4                               i,loop;
  INT4                               seed;
  REAL8                              totnorm;
  REAL8                              totnorm2;

  /* various times, frequencies, sample rates */
	INT4                               caltime;
  UINT4                              starttime;
	INT4                               lengthseg;
  INT4                               numseg;
  UINT4                              caliboffset;
  UINT4                              freqlen;
  REAL8                              deltaF, deltaT, sRate;
  INT4 site1, site2;
  REAL8 omegaRef, f0, fRef, alpha;

  /*calibration information*/
  CalibrationFunctions          calfuncs1,calfuncs2 ;
  CalibrationUpdateParams       calfacts1,calfacts2 ;

  const LALUnit countPerStrain = {0,{0,0,0,0,0,-1,1},{0,0,0,0,0,0,0}};

  /* initialize status pointer */
   INITSTATUS(status);
   ATTATCHSTATUSPTR(status);

    /* ERROR CHECKING */

 /* **** check input/output structures exist *****/

 /* output structure */

 ASSERT(MCoutput !=NULL,status,
        STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
 ASSERT(MCoutput->SSimStochBG1->data !=NULL,status,
        STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
 ASSERT(MCoutput->SSimStochBG2->data !=NULL,status,
        STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
 ASSERT(MCoutput->SSimStochBG1->data->data !=NULL,status,
        STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
 ASSERT(MCoutput->SSimStochBG2->data->data !=NULL,status,
 STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);

 /* input structure */
 ASSERT(MCinput != NULL, status,
        STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);

 ASSERT(MCinput->calfuncs1.responseFunction->data !=NULL,status,
         STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
 ASSERT(MCinput->calfuncs2.responseFunction->data !=NULL,status,
         STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
 ASSERT(MCinput->calfuncs1.responseFunction->data->data !=NULL,
         status,STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
 ASSERT(MCinput->calfuncs2.responseFunction->data->data !=NULL,
        status,STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);

 ASSERT(MCinput->calfuncs1.sensingFunction->data !=NULL,status,
         STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
 ASSERT(MCinput->calfuncs2.sensingFunction->data !=NULL,status,
         STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
 ASSERT(MCinput->calfuncs1.sensingFunction->data->data !=NULL,
        status, STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
 ASSERT(MCinput->calfuncs2.sensingFunction->data->data !=NULL,
        status,STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);

  ASSERT(MCinput->calfacts1.openLoopFactor->data !=NULL,status,
         STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
  ASSERT(MCinput->calfacts2.openLoopFactor->data !=NULL,status,
         STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
  ASSERT(MCinput->calfacts1.openLoopFactor->data->data !=NULL,
         status,STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
  ASSERT(MCinput->calfacts2.openLoopFactor->data->data !=NULL,
         status,STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);

  ASSERT(MCinput->calfacts1.sensingFactor->data !=NULL,status,
         STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
  ASSERT(MCinput->calfacts2.sensingFactor->data !=NULL,status,
         STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
  ASSERT(MCinput->calfacts1.sensingFactor->data->data !=NULL,
         status,STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
  ASSERT(MCinput->calfacts2.sensingFactor->data->data !=NULL,
         status,STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);


  /* ************ check parameter structures ***********/

  /* lengthseg is non-zero */
  ASSERT(MCparams->lengthseg > 0, status,
         STOCHASTICMCH_ENULLLEN,STOCHASTICMCH_MSGENULLLEN);

  /* numseg is non-zero */
  ASSERT(MCparams->numseg > 0, status,
         STOCHASTICMCH_ENULLSEG,STOCHASTICMCH_MSGENULLSEG);

  /* sampling rate is non zero */
  ASSERT(MCparams->sRate > 0, status,
         STOCHASTICMCH_ENULLSRATE,STOCHASTICMCH_MSGENULLSRATE);


  /* ************ done with null pointers *****************/


  /* *** check for legality ****/

  /* start frequency must not be negative */
  if (MCparams->f0 < 0)
    {
      ABORT(status,STOCHASTICMCH_ENEGFMIN,STOCHASTICMCH_MSGENEGFMIN );
    }


  /* read input stracture for calibration info*/

  calfuncs1 = MCinput->calfuncs1;
  calfuncs2 = MCinput->calfuncs2;
  calfacts1 = MCinput->calfacts1;
  calfacts2 = MCinput->calfacts2;

  /*read parameters*/

  fRef = MCparams->fRef;f0 = MCparams->f0;
  omegaRef = MCparams->omegaRef; alpha = MCparams->alpha;
  site1 = MCparams->site1; site2 = MCparams->site2;
  starttime = MCparams->starttime;
  lengthseg = MCparams->lengthseg; numseg = MCparams->numseg;
  sRate =  MCparams->sRate;
  seed = MCparams->seed;

  /*derive other parameters*/

  freqlen = lengthseg / 2 + 1;
  caliboffset =  lengthseg / (2 * sRate);
  deltaT = 1.0 / sRate;
  deltaF = 1.0/(deltaT*lengthseg);
  caltime = starttime + caliboffset;


/* * check for mismatches **/
  /* epoch */

  if (calfacts1.epoch.gpsSeconds != caltime)
    {
      ABORT(status,STOCHASTICMCH_EMMEPOCH,STOCHASTICMCH_MSGEMMEPOCH);
    }
  if (calfacts2.epoch.gpsSeconds != caltime)
    {
      ABORT(status,STOCHASTICMCH_EMMEPOCH,STOCHASTICMCH_MSGEMMEPOCH);
    }
   if (calfacts1.epoch.gpsNanoSeconds != 0)
    {
      ABORT(status,STOCHASTICMCH_EMMEPOCH,STOCHASTICMCH_MSGEMMEPOCH);
    }
  if (calfacts2.epoch.gpsNanoSeconds != 0)
    {
      ABORT(status,STOCHASTICMCH_EMMEPOCH,STOCHASTICMCH_MSGEMMEPOCH);
    }


  /* Create vectors */

  omegaGW.data = NULL;
  LALSCreateVector(status->statusPtr, &(omegaGW.data), freqlen);

  whitenedSSimStochBG1.data = NULL;
  LALSCreateVector(status->statusPtr, &(whitenedSSimStochBG1.data), lengthseg);
  whitenedSSimStochBG2.data = NULL;
  LALSCreateVector(status->statusPtr, &(whitenedSSimStochBG2.data), lengthseg);


  /* define SimulateSBParams */
  SBParams.length         = lengthseg;
  SBParams.deltaT         = 1.0/sRate;
  SBParams.detectorOne    = lalCachedDetectors[site1];
  SBParams.detectorTwo    = lalCachedDetectors[site2];
  SBParams.SSimStochBGTimeSeries1Unit = lalADCCountUnit;
  SBParams.SSimStochBGTimeSeries2Unit = lalADCCountUnit;


  /* find omegaGW */
  parametersOmega.length   = freqlen;
  parametersOmega.f0       = f0;
  parametersOmega.deltaF   = deltaF;
  parametersOmega.alpha    = alpha;
  parametersOmega.fRef     = fRef;
  parametersOmega.omegaRef = omegaRef;
  LALStochasticOmegaGW(status->statusPtr, &omegaGW, &parametersOmega);

  memset( &wFilter1, 0, sizeof(COMPLEX8FrequencySeries) );
  memset( &wFilter2, 0, sizeof(COMPLEX8FrequencySeries) );

  for (i=0;i<2;i++)
    {
      LALCCreateVector(status->statusPtr, &response[i],freqlen);
    }
  /*wFilter1.epoch.gpsSeconds = caltime;*/
  wFilter1.epoch.gpsNanoSeconds = 0;
  wFilter1.deltaF = deltaF ;
  wFilter1.f0 = f0;
  wFilter1.sampleUnits = countPerStrain;
  wFilter1.data = response[0];

  /*wFilter2.epoch.gpsSeconds = caltime;*/
  wFilter2.epoch.gpsNanoSeconds = 0;
  wFilter2.deltaF = deltaF  ;
  wFilter2.f0 = f0;
  wFilter2.sampleUnits = countPerStrain;
  wFilter2.data = response[1];


  for (loop = 0; loop < numseg; loop ++)

   {
    seed = seed + loop*2;
    SBParams.seed = seed;

    wFilter1.epoch.gpsSeconds = caltime;
    wFilter2.epoch.gpsSeconds = caltime;

    calfacts1.epoch.gpsSeconds = caltime;
    calfacts2.epoch.gpsSeconds = caltime;

    /* compute response function from calibration parameters */
    LALUpdateCalibration(status->statusPtr,&calfuncs1,&calfuncs1,&calfacts1);
    LALUpdateCalibration(status->statusPtr,&calfuncs2,&calfuncs2,&calfacts2);
    LALResponseConvert(status->statusPtr,&wFilter1,calfuncs1.responseFunction);
    LALResponseConvert(status->statusPtr,&wFilter2,calfuncs2.responseFunction);


    /* force DC to be 0 and nyquist to be real */
    response[0]->data[0].realf_FIXME = 0.0;
    response[0]->data[0].im = 0.0;
    response[1]->data[0].realf_FIXME = 0.0;
    response[1]->data[0].im = 0.0;
    response[0]->data[freqlen-1].im = 0;
    response[1]->data[freqlen-1].im = 0;


    /* define SSSimStochBGInput */
    SBInput.omegaGW                  = &omegaGW;
    SBInput.whiteningFilter1         = &wFilter1;
    SBInput.whiteningFilter2         = &wFilter2;

    SBOutput.SSimStochBG1 = &whitenedSSimStochBG1;
    SBOutput.SSimStochBG2 = &whitenedSSimStochBG2;

    /* generate whitened simulated SB data */
    LALSSSimStochBGTimeSeries(status->statusPtr, &SBOutput, &SBInput, &SBParams);
     if ( verbose )
       {

        /* Mean square */
        totnorm2=0.0;
        for (i=0;i<lengthseg;i++)
        totnorm2+=((whitenedSSimStochBG1.data->data[i])*(whitenedSSimStochBG1.data->data[i]));

        totnorm2/=lengthseg;
        printf("Mean square of whitened output is: %e\n",totnorm2);

        /* check normalizations */
        totnorm=0.0;
        for (i=1;(UINT4)i<freqlen;i++){
         REAL8 freq=i*sRate/lengthseg;
         REAL8 resp = crealf(response[0]->data[i]);
         totnorm+=resp*resp*(omegaGW.data->data[i])/(freq*freq*freq);
         }
        totnorm*=0.3*LAL_H0FAC_SI*LAL_H0FAC_SI*sRate/(LAL_PI*LAL_PI*lengthseg);
        printf("Mean square of whitened output should be: %e.  Ratio is %e\n",totnorm,totnorm/totnorm2);
        }

      for (i = 0; i < lengthseg; i ++)
	{
          MCoutput->SSimStochBG1->data->data[i + loop * lengthseg] = whitenedSSimStochBG1.data->data[i];
	  MCoutput->SSimStochBG2->data->data[i + loop * lengthseg] = whitenedSSimStochBG2.data->data[i];
        }
     caltime = caltime + caliboffset;
   }

     /* assign parameters and data to output */

   MCoutput->SSimStochBG1->f0 = f0;
   MCoutput->SSimStochBG1->deltaT = deltaT;
   MCoutput->SSimStochBG1->epoch.gpsSeconds = starttime;
   MCoutput->SSimStochBG1->epoch.gpsNanoSeconds = 0;
   MCoutput->SSimStochBG1->sampleUnits = lalADCCountUnit;
   strncpy( MCoutput->SSimStochBG1->name,
	       "Whitened-SimulatedSBOne", LALNameLength );

   MCoutput->SSimStochBG2->f0 = f0;
   MCoutput->SSimStochBG2->deltaT = deltaT;
   MCoutput->SSimStochBG2->epoch.gpsSeconds = starttime;
   MCoutput->SSimStochBG2->epoch.gpsNanoSeconds = 0;
   MCoutput->SSimStochBG2->sampleUnits =lalADCCountUnit;
   strncpy( MCoutput->SSimStochBG2->name,
	       "Whitened-SimulatedSBTwo", LALNameLength );


  /* clean up, and exit */

  /* deallocate memory for output vector */

  LALSDestroyVector(status->statusPtr, &(omegaGW.data));
  for (i=0;i<2;i++){
    LALCDestroyVector(status->statusPtr, &response[i]);
  }
  /* clean up valid data */
  LALSDestroyVector(status->statusPtr, &(whitenedSSimStochBG1.data));
  LALSDestroyVector(status->statusPtr, &(whitenedSSimStochBG2.data));


  DETATCHSTATUSPTR(status);
  RETURN(status);
}



void LALStochasticMCDsoSplice (LALStatus *status,
		         SSSimStochBGOutput *MCoutput,
	  	         StochasticMCInput  *MCinput,
		         StochasticMCParams *MCparams)
{

  /* This stores the parameters of functions used */
  StochasticOmegaGWParameters        parametersOmega;

  /* This stores the structures of SimulateSB */
  SSSimStochBGParams                 SBParams;
  SSSimStochBGInput                  SBInput;
  SSSimStochBGOutput                 SBOutput;
  REAL4TimeSeries                    whitenedSSimStochBG1;
  REAL4TimeSeries                    whitenedSSimStochBG2;

  REAL4FrequencySeries               omegaGW;
  COMPLEX8FrequencySeries            wFilter1;
  COMPLEX8FrequencySeries            wFilter2;

  REAL4Vector **longTrain1 = NULL, **shortTrain1 = NULL, **longTrain2 = NULL, **shortTrain2 = NULL;
  /* vector to store response functions of a pair of detectors  */
  COMPLEX8Vector                    *response[2]={NULL,NULL};

  /* Counters & output */
  INT4                               i,m,k,loop;
  INT4                               seed;

  /* various times, frequencies, sample rates */
	INT4                               caltime;
  UINT4                              starttime;
	INT4                               lengthseg;
	INT4                               numseg;
	INT4                               numsegsplice;
	INT4                               numsegtot;
  UINT4                              spliceoffset;
  UINT4                              freqlen;
  REAL8                              deltaF, deltaT, sRate;
  INT4 site1, site2;
  REAL8 omegaRef, f0, fRef, alpha;

  /*calibration information*/
  CalibrationFunctions          calfuncs1,calfuncs2 ;
  CalibrationUpdateParams       calfacts1,calfacts2 ;

  const LALUnit countPerStrain = {0,{0,0,0,0,0,-1,1},{0,0,0,0,0,0,0}};

  /* initialize status pointer */
   INITSTATUS(status);
   ATTATCHSTATUSPTR(status);

    /* ERROR CHECKING */

 /* **** check input/output structures exist *****/

 /* output structure */

 ASSERT(MCoutput !=NULL,status,
        STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
 ASSERT(MCoutput->SSimStochBG1->data !=NULL,status,
        STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
 ASSERT(MCoutput->SSimStochBG2->data !=NULL,status,
        STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
 ASSERT(MCoutput->SSimStochBG1->data->data !=NULL,status,
        STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
 ASSERT(MCoutput->SSimStochBG2->data->data !=NULL,status,
 STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);

 /* input structure */
 ASSERT(MCinput != NULL, status,
        STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);

 ASSERT(MCinput->calfuncs1.responseFunction->data !=NULL,status,
         STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
 ASSERT(MCinput->calfuncs2.responseFunction->data !=NULL,status,
         STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
 ASSERT(MCinput->calfuncs1.responseFunction->data->data !=NULL,
         status,STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
 ASSERT(MCinput->calfuncs2.responseFunction->data->data !=NULL,
        status,STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);

 ASSERT(MCinput->calfuncs1.sensingFunction->data !=NULL,status,
         STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
 ASSERT(MCinput->calfuncs2.sensingFunction->data !=NULL,status,
         STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
 ASSERT(MCinput->calfuncs1.sensingFunction->data->data !=NULL,
        status, STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
 ASSERT(MCinput->calfuncs2.sensingFunction->data->data !=NULL,
        status,STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);

  ASSERT(MCinput->calfacts1.openLoopFactor->data !=NULL,status,
         STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
  ASSERT(MCinput->calfacts2.openLoopFactor->data !=NULL,status,
         STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
  ASSERT(MCinput->calfacts1.openLoopFactor->data->data !=NULL,
         status,STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
  ASSERT(MCinput->calfacts2.openLoopFactor->data->data !=NULL,
         status,STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);

  ASSERT(MCinput->calfacts1.sensingFactor->data !=NULL,status,
         STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
  ASSERT(MCinput->calfacts2.sensingFactor->data !=NULL,status,
         STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
  ASSERT(MCinput->calfacts1.sensingFactor->data->data !=NULL,
         status,STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);
  ASSERT(MCinput->calfacts2.sensingFactor->data->data !=NULL,
         status,STOCHASTICMCH_ENULLP,STOCHASTICMCH_MSGENULLP);


  /* ************ check parameter structures ***********/

  /* lengthseg is non-zero */
  ASSERT(MCparams->lengthseg > 0, status,
         STOCHASTICMCH_ENULLLEN,STOCHASTICMCH_MSGENULLLEN);

  /* numseg is non-zero */
  ASSERT(MCparams->numseg > 0, status,
         STOCHASTICMCH_ENULLSEG,STOCHASTICMCH_MSGENULLSEG);

  /* sampling rate is non zero */
  ASSERT(MCparams->sRate > 0, status,
         STOCHASTICMCH_ENULLSRATE,STOCHASTICMCH_MSGENULLSRATE);


  /* ************ done with null pointers *****************/


  /* *** check for legality ****/

  /* start frequency must not be negative */
  if (MCparams->f0 < 0)
    {
      ABORT(status,STOCHASTICMCH_ENEGFMIN,STOCHASTICMCH_MSGENEGFMIN );
    }


  /* read input stracture for calibration info*/

  calfuncs1 = MCinput->calfuncs1;
  calfuncs2 = MCinput->calfuncs2;
  calfacts1 = MCinput->calfacts1;
  calfacts2 = MCinput->calfacts2;

  /*read parameters*/

  fRef = MCparams->fRef;f0 = MCparams->f0;
  omegaRef = MCparams->omegaRef; alpha = MCparams->alpha;
  site1 = MCparams->site1; site2 = MCparams->site2;
  starttime = MCparams->starttime;
  lengthseg = MCparams->lengthseg; numseg = MCparams->numseg;
  sRate =  MCparams->sRate;
  seed = MCparams->seed;

  /*derive other parameters*/

  numsegsplice = numseg - 1;
  numsegtot = numseg + numsegsplice;
  freqlen = lengthseg / 2 + 1;
  spliceoffset =  lengthseg / (2 * sRate);
  deltaT = 1.0 / sRate;
  deltaF = 1.0/(deltaT*lengthseg);
  caltime = starttime + spliceoffset;


/* * check for mismatches **/
  /* epoch */

  if (calfacts1.epoch.gpsSeconds != caltime)
    {
      ABORT(status,STOCHASTICMCH_EMMEPOCH,STOCHASTICMCH_MSGEMMEPOCH);
    }
  if (calfacts2.epoch.gpsSeconds != caltime)
    {
      ABORT(status,STOCHASTICMCH_EMMEPOCH,STOCHASTICMCH_MSGEMMEPOCH);
    }
   if (calfacts1.epoch.gpsNanoSeconds != 0)
    {
      ABORT(status,STOCHASTICMCH_EMMEPOCH,STOCHASTICMCH_MSGEMMEPOCH);
    }
  if (calfacts2.epoch.gpsNanoSeconds != 0)
    {
      ABORT(status,STOCHASTICMCH_EMMEPOCH,STOCHASTICMCH_MSGEMMEPOCH);
    }


  /* Create vectors */

  /* allocate memory for longer data train */
  longTrain1 = (REAL4Vector**) LALMalloc( sizeof(REAL4Vector*) * numseg);
  for(i=0; i<numseg; i++) {
   longTrain1[i] = NULL;
   LALSCreateVector(status->statusPtr, &longTrain1[i], lengthseg);
   for (k=0;k<lengthseg;k++) {longTrain1[i]->data[k]=0;}
  }

  longTrain2 = (REAL4Vector**) LALMalloc( sizeof(REAL4Vector*) * numseg);
  for(i=0; i<numseg; i++) {
   longTrain2[i] = NULL;
   LALSCreateVector(status->statusPtr, &longTrain2[i], lengthseg);
   for (k=0;k<lengthseg;k++) {longTrain2[i]->data[k]=0;}
   longTrain2[i] = NULL;
  }

  /* allocate memory for shorter data train */
  shortTrain1 = (REAL4Vector**) LALMalloc( sizeof(REAL4Vector*) * (numsegsplice));
  for(i=0; i<numsegsplice; i++) {
  shortTrain1[i] = NULL;
  LALSCreateVector(status->statusPtr, &shortTrain1[i], lengthseg);
  for (k=0;k<lengthseg;k++) {shortTrain1[i]->data[k]=0;}
  }
  shortTrain2 = (REAL4Vector**) LALMalloc( sizeof(REAL4Vector*) * (numsegsplice));
  for(i=0; i<numsegsplice; i++) {
  shortTrain2[i] = NULL;
  LALSCreateVector(status->statusPtr, &shortTrain2[i], lengthseg);
  for (k=0;k<lengthseg;k++) {shortTrain2[i]->data[k]=0;}
  }



  omegaGW.data = NULL;
  LALSCreateVector(status->statusPtr, &(omegaGW.data), freqlen);

  whitenedSSimStochBG1.data = NULL;
  LALSCreateVector(status->statusPtr, &(whitenedSSimStochBG1.data), lengthseg);
  whitenedSSimStochBG2.data = NULL;
  LALSCreateVector(status->statusPtr, &(whitenedSSimStochBG2.data), lengthseg);


  /* define SimulateSBParams */
  SBParams.length         = lengthseg;
  SBParams.deltaT         = 1.0/sRate;
  SBParams.detectorOne    = lalCachedDetectors[site1];
  SBParams.detectorTwo    = lalCachedDetectors[site2];
  SBParams.SSimStochBGTimeSeries1Unit = lalADCCountUnit;
  SBParams.SSimStochBGTimeSeries2Unit = lalADCCountUnit;


  /* find omegaGW */
  parametersOmega.length   = freqlen;
  parametersOmega.f0       = f0;
  parametersOmega.deltaF   = deltaF;
  parametersOmega.alpha    = alpha;
  parametersOmega.fRef     = fRef;
  parametersOmega.omegaRef = omegaRef;
  LALStochasticOmegaGW(status->statusPtr, &omegaGW, &parametersOmega);

  memset( &wFilter1, 0, sizeof(COMPLEX8FrequencySeries) );
  memset( &wFilter2, 0, sizeof(COMPLEX8FrequencySeries) );

  for (i=0;i<2;i++)
    {
      LALCCreateVector(status->statusPtr, &response[i],freqlen);
    }
  /*wFilter1.epoch.gpsSeconds = caltime;*/
  wFilter1.epoch.gpsNanoSeconds = 0;
  wFilter1.deltaF = deltaF ;
  wFilter1.f0 = f0;
  wFilter1.sampleUnits = countPerStrain;
  wFilter1.data = response[0];

  /*wFilter2.epoch.gpsSeconds = caltime;*/
  wFilter2.epoch.gpsNanoSeconds = 0;
  wFilter2.deltaF = deltaF  ;
  wFilter2.f0 = f0;
  wFilter2.sampleUnits = countPerStrain;
  wFilter2.data = response[1];


  for (loop = 0; loop < numsegtot; loop ++)
	{
		seed = seed + loop*2;
    SBParams.seed = seed;

    wFilter1.epoch.gpsSeconds = caltime;
    wFilter2.epoch.gpsSeconds = caltime;

    calfacts1.epoch.gpsSeconds = caltime;
    calfacts2.epoch.gpsSeconds = caltime;

    /* compute response function from calibration parameters */
    LALUpdateCalibration(status->statusPtr,&calfuncs1,&calfuncs1,&calfacts1);
    LALUpdateCalibration(status->statusPtr,&calfuncs2,&calfuncs2,&calfacts2);
    LALResponseConvert(status->statusPtr,&wFilter1,calfuncs1.responseFunction);
    LALResponseConvert(status->statusPtr,&wFilter2,calfuncs2.responseFunction);


    /* force DC to be 0 and nyquist to be real */
    response[0]->data[0].realf_FIXME = 0.0;
    response[0]->data[0].im = 0.0;
    response[1]->data[0].realf_FIXME = 0.0;
    response[1]->data[0].im = 0.0;
    response[0]->data[freqlen-1].im = 0;
    response[1]->data[freqlen-1].im = 0;


    /* define SSSimStochBGInput */
    SBInput.omegaGW                  = &omegaGW;
    SBInput.whiteningFilter1         = &wFilter1;
    SBInput.whiteningFilter2         = &wFilter2;

    SBOutput.SSimStochBG1 = &whitenedSSimStochBG1;
    SBOutput.SSimStochBG2 = &whitenedSSimStochBG2;

    /* generate whitened simulated SB data */
    LALSSSimStochBGTimeSeries(status->statusPtr, &SBOutput, &SBInput, &SBParams);

		/* initialise m */
		m = 0;

		if (loop%2==0)
		{
			m = (UINT4)(loop/2);
			longTrain1[m] =  whitenedSSimStochBG1.data;
			longTrain2[m] =  whitenedSSimStochBG2.data;
		}
		else
		{
			shortTrain1[m] = whitenedSSimStochBG1.data;
		  shortTrain2[m] = whitenedSSimStochBG2.data;
		}
		caltime = caltime + spliceoffset;
	}

    /* splice the long and short data trains into the output vector */
    SinusoidalSplice(longTrain1, shortTrain1,MCoutput->SSimStochBG1->data , numsegsplice, spliceoffset);
    SinusoidalSplice(longTrain2, shortTrain2,MCoutput->SSimStochBG2->data, numsegsplice, spliceoffset);
    /* deallocate memory for longer data train */
    /* for(i=0; i<numseg; i++){
    LALSDestroyVector(status->statusPtr, &longTrain1[i]);
    LALSDestroyVector(status->statusPtr, &longTrain2[i]);}
    LALFree(longTrain1); LALFree(longTrain2);*/


    /* deallocate memory for shorter data train */
    /*for(i=0; i<numsegsplice; i++){
    LALSDestroyVector(status->statusPtr, &shortTrain1[i]);
    LALSDestroyVector(status->statusPtr, &shortTrain2[i]);}
    LALFree(shortTrain1); LALFree(shortTrain2);*/

     /* assign parameters and data to output */

   MCoutput->SSimStochBG1->f0 = f0;
   MCoutput->SSimStochBG1->deltaT = deltaT;
   MCoutput->SSimStochBG1->epoch.gpsSeconds = starttime;
   MCoutput->SSimStochBG1->epoch.gpsNanoSeconds = 0;
   MCoutput->SSimStochBG1->sampleUnits = lalADCCountUnit;
   strncpy( MCoutput->SSimStochBG1->name,
	       "Whitened-SimulatedSBOne", LALNameLength );

   MCoutput->SSimStochBG2->f0 = f0;
   MCoutput->SSimStochBG2->deltaT = deltaT;
   MCoutput->SSimStochBG2->epoch.gpsSeconds = starttime;
   MCoutput->SSimStochBG2->epoch.gpsNanoSeconds = 0;
   MCoutput->SSimStochBG2->sampleUnits =lalADCCountUnit;
   strncpy( MCoutput->SSimStochBG2->name,
	       "Whitened-SimulatedSBTwo", LALNameLength );

  /* clean up, and exit */

  /* deallocate memory for output vector */

  LALSDestroyVector(status->statusPtr, &(omegaGW.data));
  for (i=0;i<2;i++){
    LALCDestroyVector(status->statusPtr, &response[i]);
  }
  /* clean up valid data */
  LALSDestroyVector(status->statusPtr, &(whitenedSSimStochBG1.data));
  LALSDestroyVector(status->statusPtr, &(whitenedSSimStochBG2.data));


  DETATCHSTATUSPTR(status);
  RETURN(status);
}



void SinusoidalSplice(REAL4Vector **longData, REAL4Vector **shortData, REAL4Vector *output,
		      UINT4 nSpliceSegs, UINT4 offset) {

  UINT4           segLen, nOutputPts, nSplicePts, lastSpliceIndexPlusOne;
  UINT4           leftOverlap, rightOverlap, i, j, iMod, jMod;
  REAL4           leftScale, rightScale, phase;

  segLen          = longData[0]->length;
  nOutputPts      = output->length;
  nSplicePts      = nSpliceSegs * segLen;
  lastSpliceIndexPlusOne = nSplicePts + offset;

  leftOverlap     = offset % segLen;
  rightOverlap    = segLen - leftOverlap;

  /* leftOverlap and rightOverlap guaranteed to be non-zero by
     CheckArguments */
  leftScale       = LAL_PI / (2.0 * leftOverlap);
  rightScale      = LAL_PI / (2.0 * rightOverlap);

  for (i=0; i<offset; i++)
    output->data[i] = longData[i/segLen]->data[i%segLen];

  /* check to make sure this is correct later */
  for(i=offset, j=0; i<lastSpliceIndexPlusOne; i++, j++) {
    iMod = i % segLen;
    jMod = j % segLen;

    /* splice from start to middle of shortData segment */
    if (iMod < leftOverlap) {
      phase = iMod * leftScale;

      output->data[i] = longData[i/segLen]->data[iMod] * sin( phase )
	+ shortData[j/segLen]->data[jMod] * cos( phase );
    }
    /* splice from middle to end of shortData segment */
    else {
      phase = (iMod - leftOverlap) * rightScale;

      output->data[i] = longData[i/segLen]->data[iMod] * cos( phase )
	+ shortData[j/segLen]->data[jMod] * sin( phase );
    }
  }

  for (i=lastSpliceIndexPlusOne; i<nOutputPts; i++)
     output->data[i] = longData[i/segLen]->data[i%segLen];

}

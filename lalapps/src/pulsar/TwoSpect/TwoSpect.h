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

#ifndef __TWOSPECT_H__
#define __TWOSPECT_H__

#include <lal/LALStdlib.h>
#include <lal/DetResponse.h>
#include <lal/AVFactories.h>
#include <gsl/gsl_rng.h>

typedef struct
{
   REAL4Vector *f;         //First PSD frequencies
   REAL4Vector *fpr;       //Second PSD frequencies
   REAL4Vector *ffdata;    //Doubly Fourier transformed data
   REAL4Vector *backgrnd;  //Noise background
   REAL4Vector *antweights;   //Antenna pattern weights (F**2)
} ffdataStruct;

/* typedef struct 
{
   REAL4 amplitude;     //Amplitude of signal
   REAL4 f0;            //Frequency of signal
   REAL4 period0;       //Period of signal
   REAL4 moddepth0;     //Modulation depth of signal
   REAL8 Tobs;          //Observation time
   REAL4 Tcoh;          //Coherence time
   REAL4 samplerate;    //Sampling rate
   REAL4 fmin;          //Minimum frequency
   REAL4 fspan;         //Frequency spanned
   REAL4 P0;            //Frequency independent noise spectral density
   REAL4 Pslope;
   REAL4 Pquad;
   REAL8 searchstarttime;
   REAL4 ra;
   REAL4 dec;
   INT4 blksize;
   INT4Vector *SFTexistlist;
   INT4 varySFTs;
} signalParamsStruct; */

typedef struct
{
   REAL4 fmin;
   REAL4 fspan;
   REAL8 Tobs;
   REAL4 Tcoh;
   //REAL4 samplerate;
   REAL8 searchstarttime;
   REAL4 ra;
   REAL4 dec;
   REAL4 dopplerMultiplier;
   INT4 blksize;
   INT4 maxbinshift;
   //INT4Vector *SFTexistlist;
   LALDetector *det;
} inputParamsStruct;

typedef struct
{
   INT4Vector *topbins;
   INT4Vector *freqoftopbins;
} topbinsStruct;


//signalParamsStruct * new_params(void);
//void free_params(signalParamsStruct *params);

inputParamsStruct * new_inputParams(void);
void free_inputParams(inputParamsStruct *input);

ffdataStruct * new_ffdata(inputParamsStruct *param, INT4 mode);
void free_ffdata(ffdataStruct *data);
void makeFakeBinarySigFF(ffdataStruct *ffdata, inputParamsStruct *input);

topbinsStruct * new_topbinsStruct(INT4 number);
void free_topbinsStruct(topbinsStruct *topbinsstruct);
void topbins(topbinsStruct *out, ffdataStruct *in, INT4 number);

REAL4Vector * gaussRandNumVector(REAL4 sigma, UINT4 length, gsl_rng *ptrToGenerator);
REAL4Vector * expRandNumVector(REAL4 mu, UINT4 length, gsl_rng *ptrToGenerator);
REAL4Vector * readInSFTs(inputParamsStruct *input);
REAL4Vector * slideTFdata(inputParamsStruct *input, REAL4Vector *tfdata, INT4Vector *binshifts);
REAL4Vector * ffPlaneNoise(inputParamsStruct *param, REAL4Vector *rngMeans, REAL4Vector *antPatternWeights);
REAL4Vector * tfWeightMeanSubtract(REAL4Vector *tfdata, REAL4Vector *rngMeans, REAL4Vector *antPatternWeights, inputParamsStruct *params);
REAL4Vector * tfRngMeans(REAL4Vector *tfdata, INT4 numffts, INT4 numfbins, INT4 blksize);
//REAL4Vector * makeFakeBinarySigTF(signalParamsStruct *in);
REAL4Vector * makeSecondFFT(REAL4Vector *tfdata, inputParamsStruct *params);

REAL4 calculateR(REAL4Vector *ffdata, REAL4Vector *template, topbinsStruct *topbinsstruct, REAL4Vector *noise);
REAL4 gaussRandNum(REAL4 sigma, gsl_rng *ptrToGenerator);
REAL4 expRandNum(REAL4 mu, gsl_rng *ptrToGenerator);
REAL4 maxModDepth(REAL4 period, REAL4 cohtime);
REAL4 minPeriod(REAL4 moddepth, REAL4 cohtime);
REAL4 sincxoverxsqminusone(REAL4 overage);
REAL4 sinc(REAL4 x);
REAL4 calcMean(REAL4Vector *vector);
REAL4 calcStddev(REAL4Vector *vector);
REAL4 calcRms(REAL4Vector *vector);



#endif




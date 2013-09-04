/*
*  Copyright (C) 2010, 2011 Evan Goetz
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

#ifndef __TWOSPECTTYPES_H__
#define __TWOSPECTTYPES_H__

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/LALDetectors.h>

#include <gsl/gsl_rng.h>

typedef struct
{
   REAL4Vector *ffdata;    //Doubly Fourier transformed data
   REAL8 tfnormalization;
   REAL8 ffnormalization;
   INT4 numffts;
   INT4 numfbins;
   INT4 numfprbins;
} ffdataStruct;

typedef struct
{
   REAL8 fmin;
   REAL8 fspan;
   REAL8 Tobs;
   REAL8 Tcoh;
   REAL8 SFToverlap;
   REAL8 searchstarttime;
   REAL8 Pmin;
   REAL8 Pmax;
   REAL8 dfmin;
   REAL8 dfmax;
   REAL4 dopplerMultiplier;
   REAL8 ihsfar;
   REAL4 ihsfom;
   REAL8 ihsfomfar;
   REAL8 templatefar;
   REAL8 log10templatefar;
   REAL8 ULfmin;
   REAL8 ULfspan;
   REAL8 ULmindf;
   REAL8 ULmaxdf;
   REAL8 lineDetection;
   INT4 ihsfactor;
   INT4 harmonicNumToSearch;
   INT4 keepOnlyTopNumIHS;
   INT4 blksize;
   INT4 maxbinshift;
   INT4 mintemplatelength;
   INT4 maxtemplatelength;
   INT4 rootFindingMethod;
   INT4 numofIFOs;
   LALDetector *det;
   CHAR* sftType;
   INT4 markBadSFTs;
   INT4 FFTplanFlag;
   INT4 calcRthreshold;
   INT4 noNotchHarmonics;
   INT4 antennaOff;
   INT4 noiseWeightOff;
   INT4 printAllULvalues;
   INT4 fastchisqinv;
   INT4 useSSE;
   INT4 followUpOutsideULrange;
   INT4 randSeed;
   INT4 ULsolver;
   INT4 signalOnly;
   gsl_rng *rng;
} inputParamsStruct;

typedef struct
{
   REAL8 fsig; /* 0 value means candidate not valid */
   REAL8 period;
   REAL8 moddepth;
   REAL4 ra;
   REAL4 dec;
   REAL8 stat;
   REAL8 h0;
   REAL8 prob;
   INT4 proberrcode;
   REAL8 normalization;
} candidate;

typedef struct
{
   candidate *data;
   UINT4 length;
   UINT4 numofcandidates;
} candidateVector;

typedef struct
{
   REAL4 alpha;
   REAL4 delta;
   REAL8Vector *fsig;
   REAL8Vector *period;
   REAL8Vector *moddepth;
   REAL8Vector *ULval;
   REAL8Vector *effSNRval;
   REAL8 normalization;
} UpperLimit;

typedef struct
{
   UpperLimit *data;
   UINT4 length;
} UpperLimitVector;

typedef struct
{
   REAL4Vector *maxima;
   INT4Vector *locations;
   REAL4Vector *foms;
   REAL4Vector *maximaForEachFbin;
   INT4Vector *locationsForEachFbin;
   INT4 rows;
} ihsMaximaStruct;

typedef struct
{
   REAL4 ihs;
   INT4 loc;
} ihsVals;

typedef struct
{
   REAL4Vector *ihsfar;
   REAL4Vector *ihsdistMean;
   REAL4Vector *ihsdistSigma;
   REAL4Vector *fomfarthresh;
   REAL4Vector *ihsfomdistMean;
   REAL4Vector *ihsfomdistSigma;
   REAL4Vector *expectedIHSVector;
} ihsfarStruct;

typedef struct
{
   REAL4 far;
   REAL4 distMean;
   REAL4 distSigma;
   REAL4Vector *topRvalues;
   INT4 farerrcode;
} farStruct;

typedef struct
{
   REAL4Vector *templatedata;       //weights
   INT4Vector *pixellocations;      //pixel locations
   INT4Vector *firstfftfrequenciesofpixels;  //pixel first frequency values
   INT4Vector *secondfftfrequencies;   //pixel second frequency values
   REAL8 f0;
   REAL8 period;
   REAL8 moddepth;
} templateStruct;


#endif


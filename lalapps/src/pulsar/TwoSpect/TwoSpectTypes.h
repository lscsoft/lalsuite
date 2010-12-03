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

#ifndef __TWOSPECTTYPES_H__
#define __TWOSPECTTYPES_H__

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/LALDetectors.h>

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
   INT4 blksize;
   INT4 maxbinshift;
   INT4 templatelength;
   LALDetector *det;
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
   REAL4Vector *maxima;
   INT4Vector *locations;
   INT4 columns;
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
} templateStruct;


#endif


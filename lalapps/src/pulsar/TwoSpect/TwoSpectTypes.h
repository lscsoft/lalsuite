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
#include <lal/LALDetectors.h>

typedef struct
{
   REAL8Vector *f;         //First PSD frequencies
   REAL8Vector *fpr;       //Second PSD frequencies
   REAL8Vector *ffdata;    //Doubly Fourier transformed data
   REAL8Vector *backgrnd;  //TF Noise background
   REAL8Vector *antweights;   //Antenna pattern weights (F**2)
} ffdataStruct;

typedef struct
{
   REAL8 fmin;
   REAL8 fspan;
   REAL8 Tobs;
   REAL8 Tcoh;
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
   REAL8 snr;
   REAL8 prob;
} candidate;

typedef struct
{
   REAL8Vector *maxima;
   INT4Vector *locations;
   INT4 columns;
} ihsMaximaStruct;

typedef struct
{
   REAL8 ihs;
   INT4 loc;
} ihsVals;

typedef struct
{
   REAL8Vector *ihsfar;
   REAL8Vector *ihsdistMean;
   REAL8Vector *ihsdistSigma;
} ihsfarStruct;

typedef struct
{
   REAL8 far;
   REAL8 distMean;
   REAL8 distSigma;
   REAL8Vector *topRvalues;
} farStruct;

typedef struct
{
   REAL8Vector *templatedata;       //weights
   INT4Vector *pixellocations;      //pixel locations
   INT4Vector *firstfftfrequenciesofpixels;  //pixel frequency values
} templateStruct;


#endif


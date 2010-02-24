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
   REAL4Vector *f;         //First PSD frequencies
   REAL4Vector *fpr;       //Second PSD frequencies
   REAL4Vector *ffdata;    //Doubly Fourier transformed data
   REAL4Vector *backgrnd;  //TF Noise background
   REAL4Vector *antweights;   //Antenna pattern weights (F**2)
} ffdataStruct;

typedef struct
{
   REAL4 fmin;
   REAL4 fspan;
   REAL8 Tobs;
   REAL4 Tcoh;
   REAL8 searchstarttime;
   REAL4 ra;
   REAL4 dec;
   REAL4 dopplerMultiplier;
   INT4 blksize;
   INT4 maxbinshift;
   INT4 templatelength;
   LALDetector *det;
} inputParamsStruct;

typedef struct
{
   REAL4 fsig; /* 0 value means candidate not valid */
   REAL4 period;
   REAL4 moddepth;
   REAL4 ra;
   REAL4 dec;
   //REAL4 Tobs;
   //REAL4 Tcoh;
   //REAL4 fmin;
   //REAL4 fspan;
   REAL4 stat;
   REAL4 snr;
} candidate;

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
} farStruct;

typedef struct
{
   REAL4Vector *templatedata;       //weights
   INT4Vector *pixellocations;      //pixel locations
   INT4Vector *firstfftfrequenciesofpixels;  //pixel frequency values
} templateStruct;


#endif


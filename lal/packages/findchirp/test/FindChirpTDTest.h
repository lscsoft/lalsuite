/*
*  Copyright (C) 2007 Bernd Machenschalk, Duncan Brown, Jolien Creighton
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

#include <lal/DataBuffer.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>

#include <lal/LALRCSID.h>
NRCSID (FINDCHIRPTDTESTH,"$Id$");

int Start(
    DataSegmentVector      **dataSegVec,
    FindChirpFilterInput   **filterInput,
    FindChirpFilterParams  **filterParams,
    FindChirpSegmentVector **fcSegVec,
    FindChirpInitParams     *initParams
    );

int Stop(
    DataSegmentVector      **dataSegVec,
    FindChirpFilterInput   **filterInput,
    FindChirpFilterParams  **filterParams,
    FindChirpSegmentVector **fcSegVec,
    UINT4 numChisqBins
    );

int SPInit(
    FindChirpTmpltParams **spTmpltParams,
    FindChirpDataParams  **spDataParams,
    FindChirpInitParams     *initParams,
    REAL4 srate,
    REAL4 f_min,
    REAL4 dynRange,
    UINT4 trunc
    );

int SPFini(
    FindChirpTmpltParams **spTmpltParams,
    FindChirpDataParams  **spDataParams
    );

int TDInit(
    FindChirpDataParams **tdDataParams,
    FindChirpInitParams    *initParams,
    REAL4 srate,
    REAL4 f_min,
    REAL4 dynRange,
    UINT4 trunc
    );

int TDFini(
    FindChirpDataParams  **tdDataParams
    );

int Init(
    FindChirpTmpltParams **TmpltParams,
    FindChirpDataParams  **DataParams,
    FindChirpInitParams   *initParams,
    REAL4 srate,
    REAL4 f_min,
    REAL4 dynRange,
    UINT4 trunc
    );

int Fini(
    FindChirpTmpltParams **TmpltParams,
    FindChirpDataParams  **DataParams
    );

int SPFilter(
    DataSegmentVector *dataSegVec,
    REAL4 mass1,
    REAL4 mass2,
    FindChirpFilterInput *filterInput,
    FindChirpFilterParams *filterParams,
    FindChirpSegmentVector *fcSegVec,
    FindChirpTmpltParams *tmpltParams,
    FindChirpDataParams *dataParams
    );

int TDFilter(
    DataSegmentVector *dataSegVec,
    REAL4 mass1,
    REAL4 mass2,
    REAL4 f_max,
    FindChirpFilterInput *filterInput,
    FindChirpFilterParams *filterParams,
    FindChirpSegmentVector *fcSegVec,
    FindChirpDataParams *dataParams
    );

int MakeData(
    DataSegmentVector *dataSegVec,
    REAL4 mass1,
    REAL4 mass2,
    REAL4 srate,
    REAL4 f_min,
    REAL4 f_max
    );

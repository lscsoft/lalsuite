/*
 * Copyright (C) 2013 Reinhard Prix
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

#ifndef _CWMAKEFAKEDATA_H  // Double-include protection
#define _CWMAKEFAKEDATA_H

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \defgroup CWMakeFakeData_h Header CWMakeFakeData.h
 * \ingroup pkg_pulsarCommon
 * \author Reinhard Prix, Karl Wette
 *
 * \brief Module for generating 'fake' data containing CW signals and/or Gaussian noise.
 * This basically presents a high-level wrapper API to the lower-level CW signal-generation
 * functions in lalsuite.
 */
/*@{*/

// ---------- exported INCLUDES ----------
#include <math.h>

// gsl includes

// lal includes
#include <lal/LALDatatypes.h>
#include <lal/SFTfileIO.h>
#include <lal/PulsarDataTypes.h>
#include <lal/ConfigFile.h>
#include <lal/DetectorStates.h>

// ---------- exported TYPES ----------
/**
 * Straightforward vector type of N PulsarParams structs
 */
typedef struct tagPulsarParamsVector
{
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(ARRAY_1D(PulsarParamsVector, PulsarParams, data, UINT4, length));
#endif /* SWIG */
  UINT4 length;		//!< number of pulsar-signals
  PulsarParams *data;	//!< array of pulsar-signal parameters
} PulsarParamsVector;

/**
 * Struct controlling all the aspects of the fake data (time-series + SFTs)
 * to be produced by XLALCWMakeFakeData() and XLALCWMakeFakeMultiData()
 */
#ifdef SWIG /* SWIG interface directives */
SWIGLAL(IMMUTABLE_MEMBERS(tagCWMFDataParams, SFTWindowType));
#endif /* SWIG */
typedef struct tagCWMFDataParams
{
  REAL8 fMin;					//!< smallest frequency guaranteed to be generated [returned fMin can be smaller]
  REAL8 Band;					//!< smallest frequency band guaranteed to be generated [returned Band can be larger]
  MultiLALDetector multiIFO;			//!< detectors to generate data for
  MultiNoiseFloor multiNoiseFloor;		//!< ... and corresponding noise-floors to generate Gaussian white noise for
  MultiLIGOTimeGPSVector multiTimestamps;	//!< timestamps to generate SFTs for
  const char *SFTWindowType;			//!< window to apply to the SFT timeseries
  REAL8 SFTWindowBeta;				//!< 'beta' parameter required for *some* windows [otherwise must be 0]
  UINT4 randSeed;				//!< seed value for random-number generator
} CWMFDataParams;

// ---------- Global variables ----------

// ---------- exported prototypes [API] ----------

#ifdef SWIG // SWIG interface directives
SWIGLAL(INOUT_STRUCTS(MultiSFTVector**, multiSFTs));
SWIGLAL(INOUT_STRUCTS(MultiREAL4TimeSeries**, multiTseries));
SWIGLAL(INOUT_STRUCTS(SFTVector**, SFTVect));
SWIGLAL(INOUT_STRUCTS(REAL4TimeSeries**, Tseries));
#endif

int XLALFindSmallestValidSamplingRate ( UINT4 *n1, UINT4 n0, const LIGOTimeGPSVector *timestamps );
int XLALCWMakeFakeMultiData ( MultiSFTVector **multiSFTs, MultiREAL4TimeSeries **multiTseries,
                              const PulsarParamsVector *injectionSources, const CWMFDataParams *dataParams, const EphemerisData *edat );
int XLALCWMakeFakeData ( SFTVector **SFTVect, REAL4TimeSeries **Tseries,
                         const PulsarParamsVector *injectionSources, const CWMFDataParams *dataParams, const EphemerisData *edat );

REAL4TimeSeries *
XLALGenerateCWSignalTS ( const PulsarParams *pulsarParams, const LALDetector *site, LIGOTimeGPS startTime, REAL8 duration, REAL8 fSamp, REAL8 fHet, const EphemerisData *edat );


int XLALReadPulsarParams ( PulsarParams *pulsarParams, const LALParsedDataFile *cfgdata, const CHAR *secName );
PulsarParamsVector *XLALPulsarParamsFromFile ( const char *fname );
PulsarParamsVector *XLALPulsarParamsFromUserInput ( const char *UserInput );

PulsarParamsVector *XLALCreatePulsarParamsVector ( UINT4 numPulsars );

void XLALDestroyPulsarParamsVector ( PulsarParamsVector *ppvect );
void XLALDestroyPulsarParams ( PulsarParams *params );


/*@}*/

#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */

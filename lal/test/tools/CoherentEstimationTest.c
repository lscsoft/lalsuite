/*
*  Copyright (C) 2007 Jolien Creighton
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

#include <complex.h>
#include <lal/CoherentEstimation.h>
#include <lal/DetectorSite.h>
#include <lal/LALNoiseModels.h>
#include <lal/AVFactories.h>
#include <lal/Random.h>
#include <lal/SkyCoordinates.h>


int
main(void)
{

  static LALStatus stat;
  UINT4 i/*, j*/;
  UINT4 Ndetectors = 2; /* number of detectors */
  REAL4TimeSeries output;
  REAL4Vector *outputData;
  CoherentEstimation params;
  DetectorsData in;

  COMPLEX16ZPGFilter filterInput;
  UINT4 filLength = 3;
  COMPLEX16 buffer[3];
  COMPLEX16Vector *Z, *P;

  SkyPosition Pos;

  static RandomParams *noiseParams;
  INT4 seed = 0;
  UINT4 length = 360 * 2048;

  REAL8IIRFilter *filt = NULL;


  /* initialize CoherentEstimation structure */
  params.preProcessed = 0; /* for pre-processing to happen */

  /* detectors info */
  params.Ndetectors = Ndetectors;
  params.detectors = (LALDetector *)LALMalloc(Ndetectors * sizeof(LALDetector));
  params.detectors[0] = lalCachedDetectors[LALDetectorIndexLHODIFF];
  params.detectors[1] = lalCachedDetectors[LALDetectorIndexLLODIFF];

  /* pre-processing filters */
  filterInput.deltaT = 4.8828125E-4; /* 2048 Hz sampling */
  Z = P = NULL;
  LALZCreateVector(&stat, &Z, filLength);
  LALZCreateVector(&stat, &P, filLength);
  filterInput.zeros = Z;
  filterInput.poles = P;
  filterInput.gain = 31.55634918610406103312016;

  P->data[0] = P->data[1] = P->data[2] = -0.8;

  Z->data[0] = 0.522407749927482845109239;
  Z->data[0] += 0.4524183825710683670706658 * I;
  Z->data[1] = 0.522407749927482845109239;
  Z->data[1] += -0.4524183825710683670706658 * I;
  Z->data[2] = 0.4142135623730950899634706;

  params.filters = (REAL8IIRFilter **)LALMalloc(Ndetectors * sizeof(REAL8IIRFilter *));
  for(i=0; i<Ndetectors; i++) {
    params.filters[i] = NULL;
    LALCreateREAL8IIRFilter(&stat, params.filters + i, &filterInput);
  }

  /* sky position */
  params.position = &Pos;
  Pos.longitude = 0.524; /* right ascension, radians */
  Pos.latitude = 0.524;  /* declination, radians */
  Pos.system = COORDINATESYSTEM_EQUATORIAL;

  /* polarization */
  params.polAngle = 0.2; /* radians */

  /* create direct filter */
  filterInput.gain = 0.031689343849711;
  memcpy(buffer, P, filLength * sizeof(COMPLEX16));
  memcpy(P, Z, filLength * sizeof(COMPLEX16));
  memcpy(Z, buffer, filLength * sizeof(COMPLEX16));
  LALCreateREAL8IIRFilter(&stat, &filt, &filterInput);


  /* Initialize input */
  in.Ndetectors = Ndetectors;
  in.data = (REAL4TimeSeries *)LALMalloc(Ndetectors * sizeof(REAL4TimeSeries));
  LALCreateRandomParams( &stat, &noiseParams, seed );

  for(i=0;i<Ndetectors;i++) {

    REAL4Vector *tmpVect = NULL;

    in.data[i].epoch.gpsSeconds = 700000000;
    in.data[i].epoch.gpsNanoSeconds = 0;
    in.data[i].deltaT = 4.8828125E-4; /* 2048 Hz sampling */
    in.data[i].f0 = 0.0;

    LALCreateVector(&stat, &tmpVect, length);
    in.data[i].data = tmpVect;
    LALNormalDeviates(&stat, in.data[i].data, noiseParams);

    /*
    LALIIRFilterREAL4Vector(&stat, in.data[i].data, params.filters[0]);
    for(j=0;j<length;j++)
      printf("%g\n",in.data[i].data->data[j]);
    return 0;
    */

    LALDIIRFilterREAL4Vector(&stat, in.data[i].data, filt);
  }

  /* allocate output */
  outputData = NULL;
  LALSCreateVector(&stat, &outputData, length);
  output.data = outputData;

  /*
  LALIIRFilterREAL4Vector(&stat, in.data[0].data, params.filters[0]);
  for(i=0;i<length;i++)
    printf("%g\n",in.data[0].data->data[i]);
  */

  /* run code */
  LALDoCoherentEstimation(&stat, &output, &params, &in);

  for(i=0;i<length;i++)
    printf("%g\n",output.data->data[i]);


  return 0;
}

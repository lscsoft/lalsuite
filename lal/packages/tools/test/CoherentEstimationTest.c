#include <lal/CoherentEstimation.h>
#include <lal/DetectorSite.h>
#include <lal/LALNoiseModels.h>
#include <lal/AVFactories.h>
#include <lal/Random.h>
#include <lal/SkyCoordinates.h>

NRCSID( COHERENTESTIMATIONTESTC, "$Id$" );

int lalDebugLevel = 3;

int
main(int argc, char **argv)
{

  static LALStatus stat;
  UINT4 i, j;
  UINT4 Ndetectors = 2; /* number of detectors */
  COMPLEX8TimeSeries output;
  COMPLEX8Vector *outputData;
  CoherentEstimation params;
  DetectorsData in;

  COMPLEX8ZPGFilter filterInput;
  UINT4 filLength = 3;
  COMPLEX8 buffer[3];
  COMPLEX8Vector *Z, *P;

  SkyPosition Pos;

  static RandomParams *noiseParams;
  INT4 seed = 0;
  UINT4 length = 360 * 2048;

  REAL4IIRFilter *filt = NULL;

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
  LALCCreateVector(&stat, &Z, filLength);
  LALCCreateVector(&stat, &P, filLength);
  filterInput.zeros = Z;
  filterInput.poles = P;
  filterInput.gain.re = 31.55634918610406103312016;
  filterInput.gain.im = 0.0;

  P->data[0].re = P->data[1].re = P->data[2].re = -0.8;
  P->data[0].im = P->data[1].im = P->data[2].im = 0;

  Z->data[0].re = 0.522407749927482845109239;
  Z->data[0].im = 0.4524183825710683670706658;
  Z->data[1].re = 0.522407749927482845109239;
  Z->data[1].im = -0.4524183825710683670706658;
  Z->data[2].re = 0.4142135623730950899634706;
  Z->data[2].im = 0.0;

  params.filters = (REAL4IIRFilter **)LALMalloc(Ndetectors * sizeof(REAL4IIRFilter *));
  for(i=0; i<Ndetectors; i++) {
    params.filters[i] = NULL;
    LALCreateREAL4IIRFilter(&stat, params.filters + i, &filterInput);
  }

  /* sky position */
  params.position = &Pos;
  Pos.longitude = 0.524; /* right ascension, radians */
  Pos.latitude = 0.524;  /* declination, radians */
  Pos.system = COORDINATESYSTEM_EQUATORIAL;

  /* polarization */
  params.polAngle = 0.2; /* radians */
  
  /* create direct filter */
  filterInput.gain.re = 0.031689343849711;
  filterInput.gain.im = 0.0;
  memcpy(buffer, P, filLength * sizeof(COMPLEX8));
  memcpy(P, Z, filLength * sizeof(COMPLEX8));
  memcpy(Z, buffer, filLength * sizeof(COMPLEX8));
  LALCreateREAL4IIRFilter(&stat, &filt, &filterInput);
  

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

    LALIIRFilterREAL4Vector(&stat, in.data[i].data, filt);
  }
  
  /* allocate output */
  outputData = NULL;
  LALCCreateVector(&stat, &outputData, length);
  output.data = outputData;

  /*
  LALIIRFilterREAL4Vector(&stat, in.data[0].data, params.filters[0]);
  for(i=0;i<length;i++)
    printf("%g\n",in.data[0].data->data[i]);
  */

  /* run code */
  LALCoherentEstimation(&stat, &output, &params, &in);

  for(i=0;i<length;i++)
    printf("%g\t%g\n",output.data->data[i].re, output.data->data[i].im);


  return 0;
}

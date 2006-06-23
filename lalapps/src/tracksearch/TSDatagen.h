
#ifndef _TSDATAGEN_H
#define _TSDATAGEN_H

#include <tracksearch.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif


  /* Error Codes */
#define TSDATAGENC_NARGS 8
  /* End N-E Defines */

  /* Define Error Codes */
#define TSDATAGENC_ENORM              0
#define TSDATAGENC_ESUB               1
#define TSDATAGENC_EARGS              2
#define TSDATAGENC_EVAL               3
#define TSDATAGENC_EFILE              4
#define TSDATAGENC_EREAD              5
#define TSDATAGENC_EMEM               6
#define TSDATAGENC_EMISC              7

#define TSDATAGENC_MSGENORM          "Normal Exit"
#define TSDATAGENC_MSGESUB           "Subroutine Fail"
#define TSDATAGENC_MSGEARGS          "Arguement Parse Error"
#define TSDATAGENC_MSGEVAL           "Invalid Argument(s)"
#define TSDATAGENC_MSGEFILE          "Error Opening File"
#define TSDATAGENC_MSGEREAD          "Error Reading File"
#define TSDATAGENC_MSGEMEM           "Memory Error"
#define TSDATAGENC_MSGEMISC          "Unknown Error"

  typedef struct
  tagTSDataGenParams
  {
    REAL8          ampFinal;
    REAL8          ampInitial;
    REAL8          freqInitial;
    REAL8          freqFinal;
    REAL8          sampleFreq;
    INT4           numSamplePoints;
    REAL8          noiseAmp;
    INT4           seed;
    CHAR          *name;
    BOOLEAN        quiet;
    BOOLEAN        noiseOnly;
    REAL8          SNR;
    BOOLEAN        externalSignal;
    CHAR          *signalFileName;
    INT4           multipleInjects;
    INT4           multipleInjectSpacing;
    REAL8          noisePSDDeltaF;
    REAL8          noisePSDMaxF;
    INT4          gpsSeconds;
    INT4          gpsNanoSeconds;
  }TSDataGenParams;

#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif

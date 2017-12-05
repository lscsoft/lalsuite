/*
*  Copyright (C) 2007 Cristina Valeria Torres
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


#ifndef _TSDATAGEN_H
#define _TSDATAGEN_H

#include "tracksearch.h"

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
 

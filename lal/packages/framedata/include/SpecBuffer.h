/*----------------------------------------------------------------------- 
 * 
 * File Name: SpecBuffer.h
 *
 * Author: Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _SPECBUFFER_H
#define _SPECBUFFER_H

#ifndef _LALDATATYPES_H
#include "LALDatatypes.h"
#ifndef _LALDATATYPES_H
#define _LALDATATYPES_H
#endif
#endif

#ifndef _REALFFT_H
#include "RealFFT.h"
#ifndef _REALFFT_H
#define _REALFFT_H
#endif
#endif

#ifndef _WINDOW_H
#include "Window.h"
#ifndef _WINDOW_H
#define _WINDOW_H
#endif
#endif

NRCSID (SPECBUFFERH, "$Id$");

#define SPECBUFFER_ENULL 1
#define SPECBUFFER_ENNUL 2
#define SPECBUFFER_ESIZE 4
#define SPECBUFFER_ESZMM 8
#define SPECBUFFER_EZERO 16
#define SPECBUFFER_ENONE 32

#define SPECBUFFER_MSGENULL "Null pointer"
#define SPECBUFFER_MSGENNUL "Non-null pointer"
#define SPECBUFFER_MSGESIZE "Invalid input size"
#define SPECBUFFER_MSGESZMM "Size mismatch"
#define SPECBUFFER_MSGEZERO "Zero divide"
#define SPECBUFFER_MSGENONE "No stored spectra"

typedef struct
tagComputeSpectrumPar
{
  WindowType   windowType;
  REAL4Vector *window;
  REAL4        wss;
  RealFFTPlan *plan;
}
ComputeSpectrumPar;

typedef struct
tagSpectrumBuffer
{
  INT4                  numSpec;
  INT4                  specFilled;
  REAL4FrequencySeries *specBuffer;
  ComputeSpectrumPar   *specParams;
}
SpectrumBuffer;

typedef struct
tagSpectrumBufferPar
{
  INT4         numSpec;
  INT4         numPoints;
  WindowType   windowType;
  RealFFTPlan *plan;
}
SpectrumBufferPar;

void
ComputeSpectrum (
    Status               *status,
    REAL4FrequencySeries *spectrum,
    INT2TimeSeries       *timeSeries,
    ComputeSpectrumPar   *parameters
    );

void
CreateSpectrumBuffer (
    Status             *status,
    SpectrumBuffer    **buffer,
    SpectrumBufferPar  *params
    );

void
DestroySpectrumBuffer (
    Status          *status,
    SpectrumBuffer **buffer
    );

void
AddSpectrum (
    Status         *status,
    SpectrumBuffer *specBuffer,
    INT2TimeSeries *timeSeries
    );

void
AverageSpectrum (
    Status               *status,
    REAL4FrequencySeries *spectrum,
    SpectrumBuffer       *buffer
    );

#endif

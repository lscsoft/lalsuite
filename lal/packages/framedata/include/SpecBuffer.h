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

#include <lal/LALDatatypes.h>
#include <lal/RealFFT.h>
#include <lal/Window.h>

#ifdef  __cplusplus
extern "C" {
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
LALComputeSpectrum (
    LALStatus               *status,
    REAL4FrequencySeries *spectrum,
    INT2TimeSeries       *timeSeries,
    ComputeSpectrumPar   *parameters
    );

void
LALCreateSpectrumBuffer (
    LALStatus             *status,
    SpectrumBuffer    **buffer,
    SpectrumBufferPar  *params
    );

void
LALDestroySpectrumBuffer (
    LALStatus          *status,
    SpectrumBuffer **buffer
    );

void
LALAddSpectrum (
    LALStatus         *status,
    SpectrumBuffer *specBuffer,
    INT2TimeSeries *timeSeries
    );

void
LALAverageSpectrum (
    LALStatus               *status,
    REAL4FrequencySeries *spectrum,
    SpectrumBuffer       *buffer
    );


#ifdef  __cplusplus
}
#endif

#endif /* _SPECBUFFER_H */

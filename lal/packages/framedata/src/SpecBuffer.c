/*----------------------------------------------------------------------- 
 * 
 * File Name: SpecBuffer.c
 *
 * Author: Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _LALSTDLIB_H
#include "LALStdlib.h"
#ifndef _LALSTDLIB_H
#define _LALSTDLIB_H
#endif
#endif

#ifndef _AVFACTORIES_H
#include "AVFactories.h"
#ifndef _AVFACTORIES_H
#define _AVFACTORIES_H
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

#ifndef _VECTOROPS_H
#include "VectorOps.h"
#ifndef _VECTOROPS_H
#define _VECTOROPS_H
#endif
#endif

#ifndef _SPECBUFFER_H
#include "SpecBuffer.h"
#ifndef _SPECBUFFER_H
#define _SPECBUFFER_H
#endif
#endif

NRCSID (SPECBUFFERC, "$Id$");

void
ComputeSpectrum (
    Status               *status,
    REAL4FrequencySeries *spectrum,
    INT2TimeSeries       *timeSeries,
    ComputeSpectrumPar   *parameters
    )
{
  REAL4Vector *tmp = NULL;
  REAL4        fac;
  INT4         i;
  INT4         n;

  INITSTATUS (status, SPECBUFFERC);
  ATTATCHSTATUSPTR (status);

  /* make sure that arguments are not NULL */
  ASSERT (spectrum, status, SPECBUFFER_ENULL, SPECBUFFER_MSGENULL);
  ASSERT (timeSeries, status, SPECBUFFER_ENULL, SPECBUFFER_MSGENULL);
  ASSERT (parameters, status, SPECBUFFER_ENULL, SPECBUFFER_MSGENULL);

  /* make sure sizes are reasonable and agree */
  n = parameters->plan->size;
  ASSERT (n > 0, status, SPECBUFFER_ESIZE, SPECBUFFER_MSGESIZE);
  ASSERT (timeSeries->data->length == n, status,
          SPECBUFFER_ESZMM, SPECBUFFER_MSGESZMM);
  ASSERT (spectrum->data->length == n/2 + 1, status,
          SPECBUFFER_ESZMM, SPECBUFFER_MSGESZMM);
  ASSERT (parameters->window->length == n, status,
          SPECBUFFER_ESZMM, SPECBUFFER_MSGESZMM);

  /* create temporary vector */
  CreateVector (status->statusPtr, &tmp, timeSeries->data->length);
  CHECKSTATUSPTR (status);

  spectrum->epoch = timeSeries->epoch;
  spectrum->f0    = timeSeries->f0;

  /* make sure that these do not divide by zero... */
  ASSERT (timeSeries->deltaT, status, SPECBUFFER_EZERO, SPECBUFFER_MSGEZERO);
  ASSERT (parameters->wss, status, SPECBUFFER_EZERO, SPECBUFFER_MSGEZERO);
  spectrum->deltaF = 1/timeSeries->deltaT;
  fac              = timeSeries->deltaT/parameters->wss;

  for (i = 0; i < n; ++i)
  {
    tmp->data[i] = fac*timeSeries->data->data[i]*parameters->window->data[i];
  }

  RealPowerSpectrum (status->statusPtr, spectrum->data, tmp, parameters->plan);
  CHECKSTATUSPTR (status);

  DestroyVector (status->statusPtr, &tmp);
  CHECKSTATUSPTR (status);

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}


void
CreateSpectrumBuffer (
    Status             *status,
    SpectrumBuffer    **buffer,
    SpectrumBufferPar  *params
    )
{
  WindowParams winparams;
  INT4 specSize;
  INT4 spec;

  INITSTATUS (status, SPECBUFFERC);
  ATTATCHSTATUSPTR (status);

  /* make sure that arguments are not NULL */
  ASSERT (buffer, status, SPECBUFFER_ENULL, SPECBUFFER_MSGENULL);
  ASSERT (params, status, SPECBUFFER_ENULL, SPECBUFFER_MSGENULL);

  /* make sure that buffer does not exist */
  ASSERT (*buffer == NULL, status, SPECBUFFER_ENNUL, SPECBUFFER_MSGENNUL);
  
  /* make sure that parameters are reasonable */
  ASSERT (params->numSpec > 0, status, SPECBUFFER_ESIZE, SPECBUFFER_MSGESIZE);
  ASSERT (params->numPoints > 0, status, SPECBUFFER_ESIZE, SPECBUFFER_MSGESIZE);
  ASSERT (params->plan, status, SPECBUFFER_ENULL, SPECBUFFER_MSGENULL);

  specSize = params->numPoints/2 + 1;

  /* assign memory for buffer */
  *buffer = (SpectrumBuffer *) LALMalloc (sizeof (SpectrumBuffer));

  /* make sure that the allocation was successful */
  ASSERT (*buffer, status, SPECBUFFER_ENULL, SPECBUFFER_MSGENULL);

  /* assign memory for buffer fields */
  (*buffer)->numSpec    = params->numSpec;
  (*buffer)->specBuffer = (REAL4FrequencySeries *)
    LALMalloc (params->numSpec*sizeof(REAL4FrequencySeries));

  ASSERT ((*buffer)->specBuffer, status, SPECBUFFER_ENULL, SPECBUFFER_MSGENULL);

  (*buffer)->specParams = (ComputeSpectrumPar *)
    LALMalloc (sizeof (ComputeSpectrumPar));

  ASSERT ((*buffer)->specBuffer, status, SPECBUFFER_ENULL, SPECBUFFER_MSGENULL);

  /* allocate memory for spectra */
  for (spec = 0; spec < params->numSpec; ++spec)
  {
    REAL4FrequencySeries *thisSpec = (*buffer)->specBuffer + spec;

    thisSpec->data = NULL;
    CreateVector (status->statusPtr, &thisSpec->data, specSize);
    CHECKSTATUSPTR (status);
  }

  /* create window and assign specParams fields */
  (*buffer)->specParams->window = NULL;
  CreateVector (
      status->statusPtr,
      &(*buffer)->specParams->window,
      params->numPoints
      );
  CHECKSTATUSPTR (status);

  winparams.length = params->numPoints;
  winparams.type   = params->windowType;
  Window (status->statusPtr, (*buffer)->specParams->window, &winparams);
  CHECKSTATUSPTR (status);
  (*buffer)->specParams->windowType = winparams.type;
  (*buffer)->specParams->wss        = winparams.sumofsquares;
  (*buffer)->specParams->plan       = params->plan;

  /* no spectra have been added yet */
  (*buffer)->specFilled = 0;

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}


void
DestroySpectrumBuffer (
    Status          *status,
    SpectrumBuffer **buffer
    )
{
  INT4 spec;

  INITSTATUS (status, SPECBUFFERC);
  ATTATCHSTATUSPTR (status);

  /* make sure that arguments are not NULL */
  ASSERT (buffer, status, SPECBUFFER_ENULL, SPECBUFFER_MSGENULL);
  ASSERT (*buffer, status, SPECBUFFER_ENULL, SPECBUFFER_MSGENULL);

  /* destroy spectra */
  for (spec = 0; spec < (*buffer)->numSpec; ++spec)
  {
    REAL4FrequencySeries *thisSpec = (*buffer)->specBuffer + spec;
    DestroyVector (status->statusPtr, &thisSpec->data);
    CHECKSTATUSPTR (status);
  }

  LALFree ((*buffer)->specBuffer);

  /* destroy window */
  DestroyVector (status->statusPtr, &(*buffer)->specParams->window);
  CHECKSTATUSPTR (status);
  
  LALFree ((*buffer)->specParams);

  LALFree (*buffer);
  *buffer = NULL;

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}


void
AddSpectrum (
    Status         *status,
    SpectrumBuffer *specBuffer,
    INT2TimeSeries *timeSeries
    )
{
  INT4 whichSpec;

  INITSTATUS (status, SPECBUFFERC);
  ATTATCHSTATUSPTR (status);

  /* make sure that arguments are not NULL */
  ASSERT (specBuffer, status, SPECBUFFER_ENULL, SPECBUFFER_MSGENULL);
  ASSERT (timeSeries, status, SPECBUFFER_ENULL, SPECBUFFER_MSGENULL);
 
  /* the next spectrum to fill */
  whichSpec = specBuffer->specFilled % specBuffer->numSpec;

  /* fill it */
  ComputeSpectrum (
      status->statusPtr,
      specBuffer->specBuffer + whichSpec,
      timeSeries,
      specBuffer->specParams
      );
  CHECKSTATUSPTR (status);

  ++(specBuffer->specFilled);

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}


void
AverageSpectrum (
    Status               *status,
    REAL4FrequencySeries *spectrum,
    SpectrumBuffer       *buffer
    )
{
  INT4  i;
  INT4  spec;
  INT4  nspec;
  REAL4 fac;

  INITSTATUS (status, SPECBUFFERC);

  /* make sure that arguments are not NULL */
  ASSERT (spectrum, status, SPECBUFFER_ENULL, SPECBUFFER_MSGENULL);
  ASSERT (buffer, status, SPECBUFFER_ENULL, SPECBUFFER_MSGENULL);

  /* make sure that fields are not NULL */
  ASSERT (spectrum->data, status, SPECBUFFER_ENULL, SPECBUFFER_MSGENULL);
  ASSERT (spectrum->data->data, status, SPECBUFFER_ENULL, SPECBUFFER_MSGENULL);
  ASSERT (buffer->specBuffer, status, SPECBUFFER_ENULL, SPECBUFFER_MSGENULL);
  ASSERT (buffer->specParams, status, SPECBUFFER_ENULL, SPECBUFFER_MSGENULL);

  if (buffer->specFilled > buffer->numSpec)
    nspec = buffer->numSpec;
  else
    nspec = buffer->specFilled;

  ASSERT (nspec, status, SPECBUFFER_ENONE, SPECBUFFER_MSGENONE);
  ASSERT (spectrum->data->length > 0, status,
          SPECBUFFER_ESIZE, SPECBUFFER_MSGESIZE);

  fac = 1/(REAL4)nspec;
  for (spec = 0; spec < nspec; ++spec)
  {
    REAL4FrequencySeries *thisSpec = buffer->specBuffer + spec;
    /* should check sample rates, f0, and average epochs ... */
    ASSERT (thisSpec->data, status, SPECBUFFER_ENULL, SPECBUFFER_MSGENULL);
    ASSERT (thisSpec->data->data, status,
            SPECBUFFER_ENULL, SPECBUFFER_MSGENULL);
    ASSERT (thisSpec->data->length == spectrum->data->length, status,
            SPECBUFFER_ESZMM, SPECBUFFER_MSGESZMM);
    for (i = 0; i < spectrum->data->length; ++i)
      if (spec == 0)
        spectrum->data->data[i]  = thisSpec->data->data[i];
      else
        spectrum->data->data[i] += thisSpec->data->data[i];
  }

  for (i = 0; i < spectrum->data->length; ++i)
    spectrum->data->data[i] *= fac;

  spectrum->deltaF = buffer->specBuffer->deltaF;
  spectrum->f0     = buffer->specBuffer->f0;

  /* normal exit */
  RETURN (status);
}

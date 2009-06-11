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

#if 0  /* autodoc block */

<lalVerbatim file="SpecBufferCV">
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{SpecBuffer.c}}
\label{ss:SpecBuffer.c}

Functions for root finding.

\subsubsection*{Prototypes}
\vspace{0.1in}
\idx{LALComputeSpectrum()}
\idx{LALCreateSpectrumBuffer()}
\idx{LALDestroySpectrumBuffer()}
\idx{LALAddSpectrum()}
\idx{LALAverageSpectrum()}

\input{SpecBufferCP}

\subsubsection*{Description}

The routine \texttt{LALComputeSpectrum()} computes the spectrum of a
\texttt{INT2TimeSeries} to produce an \texttt{REAL4FrequencySeries} power
spectrum.  It requires parameters including the forward real FFT plan and an
appropriate window.  This routine is used by the routine
\texttt{LALAddSpectrum()}.

The routine \texttt{LALCreateSpectrumBuffer()} creates a spectrum buffer.  This
requires a parameter input that specifies the number of points in each time
series that will be used to construct the spectra, the number of spectra to
store, the window type to be used when constructing the spectra, and the real
FFT plan to be used.  The spectrum buffer is destroyed using the routine
\texttt{LALDestroySpectrumBuffer()}.

The routine \texttt{LALAddSpectrum()} adds the spectrum of time series data to a
spectrum buffer.  The routine \texttt{LALAverageSpectrum()} outputs the power
spectrum frequency series consisting of the average of the (most recent)
spectra that have been stored so far.

\subsubsection*{Operating Instructions}

\begin{verbatim}
const  UINT4                 numSpec   = 8;
const  UINT4                 numPoints = 1024;
static LALStatus             status;
static SpectrumBufferPar     buffPar;
static SpectrumBuffer        buff;
static INT2TimeSeries        data;
static REAL4FrequencySeries  spec;

UINT4 i;

buffPar.numSpec    = numSpec;
buffPar.numPoints  = numPoints;
buffPar.windowType = Welch;
LALCreateForwardRealFFTPlan( &status, &buffPar.plan, numPoints, 0 );
LALCreateSpectrumBuffer( &status, &buff, &buffPar );
LALI2CreateVector( &status, &data.data, numPoints );
LALSCreateVector( &status, &spec.data, numPoints/2 + 1 );

for ( i = 0; i < numSpec; ++i )
{

  /* get data here; break out if end of data */

  LALAddSpectrum (&status, buff, &data);
}

LALAverageSpectrum( &status, &spec, buff );

/* do something with the spectrum */

LALSDestroyVector( &status, &spec.data );
LALI2DestroyVector( &status, &data.data );
LALDestroySpectrumBuffer( &status, &buff );
LALDestroyRealFFTPlan( &status, &buffPar.plan );
\end{verbatim}

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}
\vfill{\footnotesize\input{SpecBufferCV}}

</lalLaTeX>

#endif /* autodoc block */


#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/RealFFT.h>
#include <lal/Window.h>
#include <lal/VectorOps.h>
#include <lal/SpecBuffer.h>

NRCSID (SPECBUFFERC, "$Id$");

/* <lalVerbatim file="SpecBufferCP"> */
void
LALComputeSpectrum (
    LALStatus            *status,
    REAL4FrequencySeries *spectrum,
    INT2TimeSeries       *timeSeries,
    ComputeSpectrumPar   *parameters
    )
{ /* </lalVerbatim> */
  REAL4Vector *tmp = NULL;
  REAL4        fac;
  UINT4        i;
  UINT4        n;

  INITSTATUS (status, "LALComputeSpectrum", SPECBUFFERC);
  ATTATCHSTATUSPTR (status);

  /* make sure that arguments are not NULL */
  ASSERT (spectrum, status, SPECBUFFERH_ENULL, SPECBUFFERH_MSGENULL);
  ASSERT (timeSeries, status, SPECBUFFERH_ENULL, SPECBUFFERH_MSGENULL);
  ASSERT (parameters, status, SPECBUFFERH_ENULL, SPECBUFFERH_MSGENULL);

  /* make sure sizes are reasonable and agree */
  n = timeSeries->data->length;
  ASSERT (n > 0, status, SPECBUFFERH_ESIZE, SPECBUFFERH_MSGESIZE);
  ASSERT (spectrum->data->length == n/2 + 1, status,
          SPECBUFFERH_ESZMM, SPECBUFFERH_MSGESZMM);
  ASSERT (parameters->window->length == n, status,
          SPECBUFFERH_ESZMM, SPECBUFFERH_MSGESZMM);

  /* create temporary vector */
  LALCreateVector (status->statusPtr, &tmp, timeSeries->data->length);
  CHECKSTATUSPTR (status);

  spectrum->epoch = timeSeries->epoch;
  spectrum->f0    = timeSeries->f0;

  /* make sure that these do not divide by zero... */
  ASSERT (timeSeries->deltaT, status, SPECBUFFERH_EZERO, SPECBUFFERH_MSGEZERO);
  ASSERT (parameters->wss, status, SPECBUFFERH_EZERO, SPECBUFFERH_MSGEZERO);
  spectrum->deltaF = 1/(timeSeries->deltaT*timeSeries->data->length);
  fac              = sqrt( timeSeries->deltaT/parameters->wss );

  for (i = 0; i < n; ++i)
  {
    tmp->data[i] = fac*timeSeries->data->data[i]*parameters->window->data[i];
  }

  LALRealPowerSpectrum (status->statusPtr, spectrum->data, tmp, parameters->plan);
  CHECKSTATUSPTR (status);

  LALDestroyVector (status->statusPtr, &tmp);
  CHECKSTATUSPTR (status);

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/* <lalVerbatim file="SpecBufferCP"> */
void
LALCreateSpectrumBuffer (
    LALStatus          *status,
    SpectrumBuffer    **buffer,
    SpectrumBufferPar  *params
    )
{ /* </lalVerbatim> */
  LALWindowParams winparams;
  INT4 specSize;
  INT4 spec;

  INITSTATUS (status, "LALCreateSpectrumBuffer", SPECBUFFERC);
  ATTATCHSTATUSPTR (status);

  /* make sure that arguments are not NULL */
  ASSERT (buffer, status, SPECBUFFERH_ENULL, SPECBUFFERH_MSGENULL);
  ASSERT (params, status, SPECBUFFERH_ENULL, SPECBUFFERH_MSGENULL);

  /* make sure that buffer does not exist */
  ASSERT (*buffer == NULL, status, SPECBUFFERH_ENNUL, SPECBUFFERH_MSGENNUL);

  /* make sure that parameters are reasonable */
  ASSERT (params->numSpec > 0, status, SPECBUFFERH_ESIZE, SPECBUFFERH_MSGESIZE);
  ASSERT (params->numPoints > 0, status, SPECBUFFERH_ESIZE, SPECBUFFERH_MSGESIZE);
  ASSERT (params->plan, status, SPECBUFFERH_ENULL, SPECBUFFERH_MSGENULL);

  specSize = params->numPoints/2 + 1;

  /* assign memory for buffer */
  *buffer = (SpectrumBuffer *) LALMalloc (sizeof (SpectrumBuffer));

  /* make sure that the allocation was successful */
  ASSERT (*buffer, status, SPECBUFFERH_ENULL, SPECBUFFERH_MSGENULL);

  /* assign memory for buffer fields */
  (*buffer)->numSpec    = params->numSpec;
  (*buffer)->specBuffer = (REAL4FrequencySeries *)
    LALMalloc (params->numSpec*sizeof(REAL4FrequencySeries));

  ASSERT ((*buffer)->specBuffer, status, SPECBUFFERH_ENULL, SPECBUFFERH_MSGENULL);

  (*buffer)->specParams = (ComputeSpectrumPar *)
    LALMalloc (sizeof (ComputeSpectrumPar));

  ASSERT ((*buffer)->specBuffer, status, SPECBUFFERH_ENULL, SPECBUFFERH_MSGENULL);

  /* allocate memory for spectra */
  for (spec = 0; spec < params->numSpec; ++spec)
  {
    REAL4FrequencySeries *thisSpec = (*buffer)->specBuffer + spec;

    thisSpec->data = NULL;
    LALCreateVector (status->statusPtr, &thisSpec->data, specSize);
    CHECKSTATUSPTR (status);
  }

  /* create window and assign specParams fields */
  (*buffer)->specParams->window = NULL;
  LALCreateVector (
      status->statusPtr,
      &(*buffer)->specParams->window,
      params->numPoints
      );
  CHECKSTATUSPTR (status);

  winparams.length = params->numPoints;
  winparams.type   = params->windowType;
  LALWindow (status->statusPtr, (*buffer)->specParams->window, &winparams);
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


/* <lalVerbatim file="SpecBufferCP"> */
void
LALDestroySpectrumBuffer (
    LALStatus       *status,
    SpectrumBuffer **buffer
    )
{ /* </lalVerbatim> */
  INT4 spec;

  INITSTATUS (status, "LALDestroySpectrumBuffer", SPECBUFFERC);
  ATTATCHSTATUSPTR (status);

  /* make sure that arguments are not NULL */
  ASSERT (buffer, status, SPECBUFFERH_ENULL, SPECBUFFERH_MSGENULL);
  ASSERT (*buffer, status, SPECBUFFERH_ENULL, SPECBUFFERH_MSGENULL);

  /* destroy spectra */
  for (spec = 0; spec < (*buffer)->numSpec; ++spec)
  {
    REAL4FrequencySeries *thisSpec = (*buffer)->specBuffer + spec;
    LALDestroyVector (status->statusPtr, &thisSpec->data);
    CHECKSTATUSPTR (status);
  }

  LALFree ((*buffer)->specBuffer);

  /* destroy window */
  LALDestroyVector (status->statusPtr, &(*buffer)->specParams->window);
  CHECKSTATUSPTR (status);

  LALFree ((*buffer)->specParams);

  LALFree (*buffer);
  *buffer = NULL;

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/* <lalVerbatim file="SpecBufferCP"> */
void
LALAddSpectrum (
    LALStatus      *status,
    SpectrumBuffer *specBuffer,
    INT2TimeSeries *timeSeries
    )
{ /* </lalVerbatim> */
  INT4 whichSpec;

  INITSTATUS (status, "LALAddSpectrum", SPECBUFFERC);
  ATTATCHSTATUSPTR (status);

  /* make sure that arguments are not NULL */
  ASSERT (specBuffer, status, SPECBUFFERH_ENULL, SPECBUFFERH_MSGENULL);
  ASSERT (timeSeries, status, SPECBUFFERH_ENULL, SPECBUFFERH_MSGENULL);

  /* the next spectrum to fill */
  whichSpec = specBuffer->specFilled % specBuffer->numSpec;

  /* fill it */
  LALComputeSpectrum (
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


/* <lalVerbatim file="SpecBufferCP"> */
void
LALAverageSpectrum (
    LALStatus            *status,
    REAL4FrequencySeries *spectrum,
    SpectrumBuffer       *buffer
    )
{ /* </lalVerbatim> */
  UINT4 i;
  INT4  spec;
  INT4  nspec;
  REAL4 fac;

  INITSTATUS (status, "LALAverageSpectrum", SPECBUFFERC);

  /* make sure that arguments are not NULL */
  ASSERT (spectrum, status, SPECBUFFERH_ENULL, SPECBUFFERH_MSGENULL);
  ASSERT (buffer, status, SPECBUFFERH_ENULL, SPECBUFFERH_MSGENULL);

  /* make sure that fields are not NULL */
  ASSERT (spectrum->data, status, SPECBUFFERH_ENULL, SPECBUFFERH_MSGENULL);
  ASSERT (spectrum->data->data, status, SPECBUFFERH_ENULL, SPECBUFFERH_MSGENULL);
  ASSERT (buffer->specBuffer, status, SPECBUFFERH_ENULL, SPECBUFFERH_MSGENULL);
  ASSERT (buffer->specParams, status, SPECBUFFERH_ENULL, SPECBUFFERH_MSGENULL);

  if (buffer->specFilled > buffer->numSpec)
    nspec = buffer->numSpec;
  else
    nspec = buffer->specFilled;

  ASSERT (nspec, status, SPECBUFFERH_ENONE, SPECBUFFERH_MSGENONE);
  ASSERT (spectrum->data->length > 0, status,
          SPECBUFFERH_ESIZE, SPECBUFFERH_MSGESIZE);

  fac = 1/(REAL4)nspec;
  for (spec = 0; spec < nspec; ++spec)
  {
    REAL4FrequencySeries *thisSpec = buffer->specBuffer + spec;
    /* should check sample rates, f0, and average epochs ... */
    ASSERT (thisSpec->data, status, SPECBUFFERH_ENULL, SPECBUFFERH_MSGENULL);
    ASSERT (thisSpec->data->data, status,
            SPECBUFFERH_ENULL, SPECBUFFERH_MSGENULL);
    ASSERT (thisSpec->data->length == spectrum->data->length, status,
            SPECBUFFERH_ESZMM, SPECBUFFERH_MSGESZMM);
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

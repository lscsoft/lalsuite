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

<lalVerbatim file="FrameDataCV">
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{FrameData.c}}
\label{ss:FrameData.c}

Functions for reading frame data.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{FrameDataCP}
\idx{LALInitializeFrameData()}
\idx{LALFinalizeFrameData()}
\idx{LALGetFrameData()}
\idx{LALGetFrameDataResponse()}

\subsubsection*{Description}

The routine \texttt{LALInitializeFrameData()} searches for frame files in a
specified directory and performs the necessary preparation for reading these
files.  When the user is finished reading frame data from the path, the
routine \texttt{LALFinalizeFrameData()} should be called.

The routine \texttt{LALGetFrameData()} gets the next frame IFO\_DMRO data while
the routine \texttt{LALGetFrameDataResponse()} gets the current response
function.  The routine \texttt{LALGetFrameDataResponse()} does a spline-fit to
the sweptsine response of the instrument in order to get the response
function.

If the output time series has a \texttt{NULL} data field, then
\texttt{LALGetFrameData()} enters seek mode in which no data is returned but
the frame data is advanced the required amount.

\subsubsection*{Operating Instructions}

\begin{verbatim}
const  UINT4                    numPoints = 1024;
const  CHAR                    *framePath = "/data/frames";
static Status                   status;
static FrameData                frameData;
static INT2TimeSeries           dmro;
static COMPLEX8FrequencySeries  resp;

LALI2CreateVector( &status, &dmro.data, numPoints );
LALCCreateVector( &status, &resp.data, numPoints/2 + 1 );
LALInitializeFrameData( &status, &frameData, framePath );

/* infinite loop reading frame data */
while ( 1 )
{
  LALGetFrameData( &status, &dmro, frameData );

  /* break out of loop if end of data */
  if ( frameData->endOfData )
  {
    break;
  }

  /* get response function if new calibration info */
  if ( frameData->newCalibration )
  {
    LALGetFrameDataResponse( &status, &resp, frameData );
  }

  /* seek 3 minutes into each new locked section */
  if ( frameData->newLock )
  {
    INT2TimeSeries seek;
    INT2Vector     svec;

    svec.length = 180/dmro.deltaT; /* 3 minutes */
    svec.data   = NULL;            /* seek mode */
    seek.data   = &svec;
    LALGetFrameData( &status, &seek, frameData );

    /* break out of loop if end of data */
    if ( frameData->endOfData )
    {
      break;
    }

    /* get response function if new calibration info */
    if ( frameData->newCalibration )
    {
      LALGetFrameDataResponse( &status, &resp, frameData );
    }

    /* go back to the beginning of the infinite loop */
    continue;
  }

  /* do something with the data here */

}

LALFinalizeFrameData( &status, &frameData );
LALCDestroyVector( &status, &resp.data );
LALI2DestroyVector( &status, &dmro.data );
\end{verbatim}

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}
\vfill{\footnotesize\input{FrameDataCV}}

</lalLaTeX>

#endif /* autodoc block */


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/SeqFactories.h>
#include <lal/FrameData.h>
#include <lal/LALFrameL.h>

NRCSID (FRAMEDATAC, "$Id$");

/* <lalVerbatim file="FrameDataCP"> */
void
LALInitializeFrameData (
    LALStatus  *status,
    FrameData **frameData,
    CHAR       *framePath
    )
{ /* </lalVerbatim> */
  const CHAR *headNames[]       = {"C1-*.F", "H-*.F", "H-*.T", "L-*.F",
                                   "L-*.T", "C1-*[0-9]"};
  const INT4  numHeadNames      = 6;
  const INT4  maxNumFiles       = 2048;
  const INT4  maxFileNameLength = 256;

  CreateVectorSequenceIn frameFileNameIn;
  CHAR                   command[1024];
  INT4                   nameType;

  INITSTATUS (status, "LALInitializeFrameData", FRAMEDATAC);
  ATTATCHSTATUSPTR (status);

  /* make sure arguments are reasonable */
  ASSERT (framePath, status, FRAMEDATAH_ENULL, FRAMEDATAH_MSGENULL);
  ASSERT (frameData, status, FRAMEDATAH_ENULL, FRAMEDATAH_MSGENULL);
  ASSERT (!(*frameData), status, FRAMEDATAH_ENNUL, FRAMEDATAH_MSGENNUL);

  /* allocate memory */
  *frameData = LALCalloc (1, sizeof(FrameData));
  ASSERT (*frameData, status, FRAMEDATAH_ENULL, FRAMEDATAH_MSGENULL);

  /* debuglevel zero: don't report errors */
  FrLibIni (NULL, stderr, 0);

  /* set some flags and modes */
  (*frameData)->dataBreak = 1;
  (*frameData)->newLock   = 1;
  (*frameData)->inLock    = 1;

  /* set up frame files */

  frameFileNameIn.length       = maxNumFiles;
  frameFileNameIn.vectorLength = maxFileNameLength;
  LALCHARCreateVectorSequence (
      status->statusPtr,
      &(*frameData)->frameFileNames,
      &frameFileNameIn
      );
  CHECKSTATUSPTR (status);

  /* search directory for frame files */
  for (nameType = 0; nameType < numHeadNames; ++nameType)
  {
    FILE *fp;
    INT4  nbytes;
    INT4  numFiles;

    /* command to list frame files of current name type */
    nbytes = sprintf (command, "ls %s/%s 2>/dev/null",
                      framePath, headNames[nameType]);
    ASSERT (nbytes > 0, status, FRAMEDATAH_EREAD, FRAMEDATAH_MSGEREAD);
    ASSERT (nbytes < (INT4)sizeof(command), status,
            FRAMEDATAH_EREAD, FRAMEDATAH_MSGEREAD);

    /* fp is a stream containing the filenames */
    fp = popen (command, "r");
    ASSERT (fp, status, FRAMEDATAH_EREAD, FRAMEDATAH_MSGEREAD);

    /* read the filenames into the stored filename list */
    numFiles = (*frameData)->numFiles;
    while (numFiles < maxNumFiles)
    {
      CHAR *fileName;
      fileName  = (*frameData)->frameFileNames->data;
      fileName += numFiles*maxFileNameLength;
      if (EOF == fscanf (fp, "%s\n", fileName))
        break;
      ASSERT ((INT4)strlen(fileName) < maxFileNameLength, status,
              FRAMEDATAH_EREAD, FRAMEDATAH_MSGEREAD);
      ++numFiles;
    }
    (*frameData)->numFiles = numFiles;

    pclose (fp);
  }
  ASSERT ((*frameData)->numFiles > 0, status,
          FRAMEDATAH_EREAD, FRAMEDATAH_MSGEREAD);
  ASSERT ((*frameData)->numFiles < maxNumFiles, status,
          FRAMEDATAH_EREAD, FRAMEDATAH_MSGEREAD);

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/* <lalVerbatim file="FrameDataCP"> */
void
LALFinalizeFrameData (
    LALStatus  *status,
    FrameData **frameData
    )
{ /* </lalVerbatim> */
  INITSTATUS (status, "LALFinalizeFrameData", FRAMEDATAC);
  ATTATCHSTATUSPTR (status);

  /* make sure argument is reasonable */
  ASSERT (frameData, status, FRAMEDATAH_ENULL, FRAMEDATAH_MSGENULL);
  ASSERT (*frameData, status, FRAMEDATAH_ENULL, FRAMEDATAH_MSGENULL);

  /* free an existing frame */
  if ((*frameData)->frame)
  {
    FrameFree ((struct FrameH *)((*frameData)->frame));
  }

  /* close an open file */
  if ((*frameData)->fileOpen)
  {
    FrFileOEnd ((*frameData)->frameFile);
  }

  /* destroy file name vector sequence */
  LALCHARDestroyVectorSequence (status->statusPtr, &(*frameData)->frameFileNames);
  CHECKSTATUSPTR (status);

  /* free memory */
  LALFree (*frameData);
  *frameData = NULL;

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


static void
GetNewFrame (
    LALStatus *status,
    FrameData *frameData
    )
{
  const REAL4 resolution = 3.2e-3;

  INITSTATUS (status, "GetNewFrame", FRAMEDATAC);

  /* make sure argument is not NULL */
  ASSERT (frameData, status, FRAMEDATAH_ENULL, FRAMEDATAH_MSGENULL);

  /* free an existing frame */
  if (frameData->frame)
  {
    FrameFree ((struct FrameH *)(frameData->frame));
  }

  /* get a new frame from the frame file */
  frameData->frame = FrameRead ((struct FrFile *)(frameData->frameFile));

  if (frameData->frame)
  {
    struct FrameH    *frame = frameData->frame;
    struct FrAdcData *dmro;
    struct FrAdcData *lock;
    REAL8 startTime;
    REAL8 finishTime;
    INT4  oldCalibTime;
    INT4  newCalibTime;

    /* set old calibration time */
    oldCalibTime = frameData->calibrationTime.gpsSeconds;

    /* get calibration info */
    {
      char name[]="sweptsine"; /* hack to get non-const string */
      frameData->calibration =
        FrStatDataFind (frame->detectProc, name, frame->GTimeS);
    }
    ASSERT (frameData->calibration, status,
            FRAMEDATAH_ENOSS, FRAMEDATAH_MSGENOSS);

    /* get new calibration time */
    newCalibTime = ((struct FrStatData *)(frameData->calibration))->timeStart;

    /* indicate if calibration data is new */
    if (newCalibTime != oldCalibTime)
    {
      frameData->calibrationTime.gpsSeconds     = newCalibTime;
      frameData->calibrationTime.gpsNanoSeconds = 0;
      frameData->newCalibration                 = 1;
    }

    /* get IFO_DMRO */
    {
      char name[]="IFO_DMRO"; /* hack to get non-const string */
      frameData->dmro = dmro = FrAdcDataFind (frame, name);
    }
    ASSERT (dmro, status, FRAMEDATAH_EDMRO, FRAMEDATAH_MSGEDMRO);

    frameData->numDmro = dmro->data->nData;
    frameData->curDmro = 0;

    /* get IFO_Lock */
    {
      char name[]="IFO_Lock"; /* hack to get non-const string */
      frameData->lock = lock = FrAdcDataFind (frame, name);
    }
    ASSERT (lock, status, FRAMEDATAH_ELOCK, FRAMEDATAH_MSGELOCK);

    frameData->numLock = lock->data->nData;
    frameData->curLock = 0;

    /* ratio is number of IFO_DMRO samples to each IFO_Lock sample */
    frameData->ratio = frameData->numDmro/frameData->numLock;

    /* get lock-low and lock-high values */
    {
      char name[] = "locklo/lockhi"; /* hack to get non-const string */
      frameData->lockLowHigh =
        FrStatDataFind (frame->detectProc, name, frame->GTimeS);
    }
    ASSERT (frameData->lockLowHigh, status,
            FRAMEDATAH_ELOHI, FRAMEDATAH_MSGELOHI);

    frameData->lockLow =
      ((struct FrStatData *)(frameData->lockLowHigh))->data->dataS[0];
    frameData->lockHigh =
      ((struct FrStatData *)(frameData->lockLowHigh))->data->dataS[1];

    /* get sample rate of IFO_DMRO */
    frameData->sampleRate = dmro->sampleRate;

    /* set start time for this frame */
    frameData->frameStartTime.gpsSeconds     = frame->GTimeS;
    frameData->frameStartTime.gpsNanoSeconds = frame->GTimeN;

    /* make sure frame data immediately follows previous frame's data */
    startTime  = frameData->frameStartTime.gpsSeconds
      + 1e-9*frameData->frameStartTime.gpsNanoSeconds;
    finishTime = frameData->frameFinishTime.gpsSeconds
      + 1e-9*frameData->frameFinishTime.gpsNanoSeconds;
    if (!frameData->dataBreak && fabs(startTime - finishTime) > resolution)
    {
      frameData->dataBreak = 1;
    }
    else
    {
      frameData->dataBreak = 0;
    }

    /* set finish time for this frame */
    finishTime = startTime + frame->dt;
    frameData->frameFinishTime.gpsSeconds     = (INT4)finishTime;
    frameData->frameFinishTime.gpsNanoSeconds = (INT4)fmod(1e9*finishTime, 1e9);
  }
  else
  {
    /* no more frames in file: close the current frame file */
    FrFileOEnd (frameData->frameFile);
    frameData->fileOpen = 0;
  }

  RETURN (status);
}


/* <lalVerbatim file="FrameDataCP"> */
void
LALGetFrameData (
    LALStatus      *status,
    INT2TimeSeries *data,
    FrameData      *frameData
    )
{ /* </lalVerbatim> */
  INT4 seek       = 0;
  INT4 brokenLock = 0;
  INT4 numPoints;
  INT4 needed;

  INITSTATUS (status, "LALGetFrameData", FRAMEDATAC);
  ATTATCHSTATUSPTR (status);

  /* make sure arguments are reasonable */
  ASSERT (frameData, status, FRAMEDATAH_ENULL, FRAMEDATAH_MSGENULL);
  ASSERT (data, status, FRAMEDATAH_ENULL, FRAMEDATAH_MSGENULL);
  ASSERT (data->data, status, FRAMEDATAH_ENULL, FRAMEDATAH_MSGENULL);
  numPoints = data->data->length;
  ASSERT (numPoints, status, FRAMEDATAH_ESIZE, FRAMEDATAH_MSGESIZE);

  /* if data->data->data is NULL enter seek mode */
  if (data->data->data == NULL)
  {
    seek = 1;
  }

  /* initialize */
  needed = numPoints;
  frameData->newCalibration = 0;

  /* loop until all requested points are obtained */
  while (needed > 0)
  {

    /* do we need to open a new frame file? */
    if (!frameData->fileOpen)
    {
      CHAR *fileName;

      /* check to see if we have read all the files in list */
      if (frameData->fileNum >= frameData->numFiles)
      {
        /* we are done */
        frameData->endOfData = 1;
        DETATCHSTATUSPTR (status);
        RETURN (status);
      }

      /* open next file in list */
      fileName = frameData->frameFileNames->data;
      fileName += frameData->fileNum*frameData->frameFileNames->vectorLength;
      frameData->frameFile = FrFileINew (fileName);
      ASSERT (frameData->frameFile, status,
              FRAMEDATAH_EOPEN, FRAMEDATAH_MSGEOPEN);
      frameData->fileOpen = 1;
      ++frameData->fileNum;
    }

    /* are there points left in current frame? */
    if (frameData->numDmro > frameData->curDmro)
    {
      short *dmroData;
      short *lockData;
      INT4   top;
      INT4   imax;
      INT4   i;

      dmroData = ((struct FrAdcData *)(frameData->dmro))->data->dataS;
      lockData = ((struct FrAdcData *)(frameData->lock))->data->dataS;

      /* find next point in lock stream covering needed range */
      top = (frameData->curDmro + needed)/frameData->ratio;
      if (top%frameData->ratio != 0)
        ++top;
      imax = (top < frameData->numLock ? top : frameData->numLock);

      /* scan to see if we drop out of lock during the interval */
      i = frameData->curLock;
      while (frameData->inLock && i < imax)
      {
        if (lockData[i] < frameData->lockLow)
          break;
        if (lockData[i] > frameData->lockHigh)
          break;
        ++i;
      }

      /* did we drop out of lock (or do we care)? */
      if (frameData->inLock && i < imax) /* lost lock: look for return */
      {
        imax = frameData->numLock;
        while (i < imax)
        {
          if (lockData[i] >= frameData->lockLow
              && lockData[i] <= frameData->lockHigh)
            break;
          ++i;
        }
        frameData->curLock = i;
        frameData->curDmro = i*frameData->ratio;
        brokenLock         = 1;
        needed             = numPoints;
      }
      else /* in lock or don't care: get data */
      {
        INT4  remaining = frameData->numDmro - frameData->curDmro;
        INT4  ncopy     = remaining < needed ? remaining : needed;
        INT4  wordSize  = ((struct FrAdcData *)(frameData->dmro))->data->wSize;
        INT4  wordRatio = wordSize > (INT4)sizeof(INT2) ? wordSize/sizeof(INT2) : 1;

        /* copy the data if desired */
        if (!seek)
        {
          INT2 *location = data->data->data + (numPoints - needed)*wordRatio;
          memcpy (location, dmroData + frameData->curDmro, ncopy*wordSize);
        }

        /* if starting buffer: mark start time and time step */
        if (needed == numPoints)
        {
          REAL8 starttime;
          starttime  = frameData->frameStartTime.gpsSeconds;
          starttime += 1e-9*frameData->frameStartTime.gpsNanoSeconds;
          starttime += frameData->curDmro/frameData->sampleRate;
          data->epoch.gpsSeconds     = (INT4)starttime;
          data->epoch.gpsNanoSeconds = (INT4)fmod(1e9*starttime, 1e9);
          data->deltaT               = 1/frameData->sampleRate;
          data->f0                   = 0;
        }

        /* bookkeeping */
        needed  -= ncopy;
        frameData->curDmro += ncopy;
        frameData->curLock  = frameData->curDmro/frameData->ratio;
      }

    }
    else /* no more points left in current frame */
    {
      GetNewFrame (status->statusPtr, frameData);
      CHECKSTATUSPTR (status);

      /* check to see if there is a gap in the data */
      if (frameData->dataBreak)
      {
        brokenLock = 1;
        needed     = numPoints;
      }
    }

  }

  /* if lock was broken or there is a data gap, set newLock flag */
  if (brokenLock == 1)
  {
    frameData->newLock = 1;
  }
  else
  {
    frameData->newLock = 0;
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}



static void
SplineFit (
    LALStatus   *status,
    REAL4Vector *yout,
    REAL4Vector *yinp,
    REAL4Vector *xinp
    )
{
  REAL4Vector *yppvec = NULL;
  REAL4       *ypp;
  UINT4        n;

  INITSTATUS (status, "SplineFit", FRAMEDATAC);
  ATTATCHSTATUSPTR (status);

  /* make sure arguments are reasonable */

  ASSERT (yout, status, FRAMEDATAH_ENULL, FRAMEDATAH_MSGENULL);
  ASSERT (yout->data, status, FRAMEDATAH_ENULL, FRAMEDATAH_MSGENULL);
  ASSERT (yout->length > 2, status, FRAMEDATAH_ESIZE, FRAMEDATAH_MSGESIZE);

  ASSERT (yinp, status, FRAMEDATAH_ENULL, FRAMEDATAH_MSGENULL);
  ASSERT (yinp->data, status, FRAMEDATAH_ENULL, FRAMEDATAH_MSGENULL);
  n = yinp->length;
  ASSERT (n > 2, status, FRAMEDATAH_ESIZE, FRAMEDATAH_MSGESIZE);

  ASSERT (xinp, status, FRAMEDATAH_ENULL, FRAMEDATAH_MSGENULL);
  ASSERT (xinp->data, status, FRAMEDATAH_ENULL, FRAMEDATAH_MSGENULL);
  ASSERT (xinp->length == n, status, FRAMEDATAH_ESIZE, FRAMEDATAH_MSGESIZE);

  /* create temporary vector */
  LALCreateVector (status->statusPtr, &yppvec, n);
  CHECKSTATUSPTR (status);
  ypp = yppvec->data;

  { /* setup second derivative array */

    REAL4Vector *uvec = NULL;
    REAL4       *u;

    REAL4 *x = xinp->data;
    REAL4 *y = yinp->data;

    UINT4  i;
    INT4   j;

    LALCreateVector (status->statusPtr, &uvec, n);
    CHECKSTATUSPTR (status);
    u = uvec->data;

    /* decomposition loop */
    ypp[0] = u[0] = 0;
    for (i = 1; i < n - 1; ++i)
    {
      REAL4 dx0   = x[i]   - x[i-1];
      REAL4 dx1   = x[i+1] - x[i];
      REAL4 dx2   = x[i+1] - x[i-1];
      REAL4 dy0   = y[i]   - y[i-1];
      REAL4 dy1   = y[i+1] - y[i];
      REAL4 ddydx = dy1/dx1 - dy0/dx0;
      REAL4 sigma = dx0/dx2;
      REAL4 fac   = 1/(sigma*ypp[i-1] + 2);

      ypp[i] = fac*(sigma - 1);
      u[i]   = fac*(6*ddydx/dx2 - sigma*u[i-1]);
    }

    /* backsubstitution loop */
    ypp[n-1] = 0;
    for (j = n - 2; j >= 0; --j)
    {
      ypp[j] = ypp[j]*ypp[j+1] + u[j];
    }

    LALDestroyVector (status->statusPtr, &uvec);
    CHECKSTATUSPTR (status);

  }

  { /* do fit */

    REAL4 *x  = xinp->data;
    REAL4 *y  = yinp->data;
    REAL4  dx = (x[n - 1] - x[0])/(yout->length - 1);

    UINT4 i;
    UINT4 j = 0;

    for (i = 0; i < yout->length; ++i)
    {
      REAL4 xx = x[0] + i*dx;
      REAL4 a;
      REAL4 b;
      REAL4 c;
      REAL4 d;
      REAL4 h;
      REAL4 dx0;
      REAL4 dx1;

      while (j < n - 1 && x[j] <= xx)
      {
        ++j;
      }

      h   = x[j] - x[j-1];
      dx0 = xx - x[j-1];
      dx1 = x[j] - xx;
      a   = dx1/h;
      b   = dx0/h;
      c   = dx1*(dx1*a - h)/6;
      d   = dx0*(dx0*b - h)/6;

      yout->data[i] = a*y[j-1] + b*y[j] + c*ypp[j-1] + d*ypp[j];
    }

  }

  LALDestroyVector (status->statusPtr, &yppvec);
  CHECKSTATUSPTR (status);

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/* <lalVerbatim file="FrameDataCP"> */
void
LALGetFrameDataResponse (
    LALStatus               *status,
    COMPLEX8FrequencySeries *response,
    FrameData               *frameData
    )
{ /* </lalVerbatim> */
  REAL4Vector *re = NULL;
  REAL4Vector *im = NULL;
  REAL4Vector  x;
  REAL4Vector  y;

  REAL4  factor;

  REAL4 *fri;
  REAL4 *ssf;
  REAL4 *ssr;
  REAL4 *ssi;
  INT4   ssn;
  INT4   n;
  UINT4  i;

  INITSTATUS (status, "LALGetFrameDataResponse", FRAMEDATAC);
  ATTATCHSTATUSPTR (status);

  ASSERT (frameData, status, FRAMEDATAH_ENULL, FRAMEDATAH_MSGENULL);
  ASSERT (response, status, FRAMEDATAH_ENULL, FRAMEDATAH_MSGENULL);
  ASSERT (response->data, status, FRAMEDATAH_ENULL, FRAMEDATAH_MSGENULL);
  ASSERT (response->data->data, status, FRAMEDATAH_ENULL, FRAMEDATAH_MSGENULL);
  ASSERT (response->data->length > 2, status,
          FRAMEDATAH_ESIZE, FRAMEDATAH_MSGESIZE);

  fri = ((struct FrStatData *)(frameData->calibration))->data->dataF;
  n   = ((struct FrStatData *)(frameData->calibration))->data->nData;
  ssn = n/3;
  ssf = fri;
  ssr = ssf + ssn;
  ssi = ssr + ssn;
  ASSERT (ssn > 2, status, FRAMEDATAH_ESSSZ, FRAMEDATAH_MSGESSSZ);

  LALCreateVector (status->statusPtr, &re, response->data->length);
  CHECKSTATUSPTR (status);
  LALCreateVector (status->statusPtr, &im, response->data->length);
  CHECKSTATUSPTR (status);

  x.data   = ssf;
  x.length = ssn;
  y.length = ssn;

  /* spline fit to real part */
  y.data = ssr;
  SplineFit (status->statusPtr, re, &y, &x);
  CHECKSTATUSPTR (status);

  /* spline fit to imagingary part */
  y.data = ssi;
  SplineFit (status->statusPtr, im, &y, &x);
  CHECKSTATUSPTR (status);

  factor  = sqrt(9.35);        /* square root bandwidth     */
  factor /= 21399;             /* Spero's constant k        */
  factor *= 10/(REAL4)2048;    /* volts IFO / ADC bit       */
  factor /= -4*LAL_PI*LAL_PI;  /* from two time derivatives */

  /* demand DC is zero */
  response->data->data[0].re = response->data->data[0].im = 0;

  /* other components (including possible Nyquist) */
  for (i = 1; i < response->data->length; ++i)
  {
    REAL4 f     = 0.5*i*frameData->sampleRate/response->data->length;
    REAL4 fsq   = f*f;
    REAL4 ssre  = re->data[i];
    REAL4 ssim  = im->data[i];
    REAL4 ssmod = ssre*ssre + ssim*ssim;
    REAL4 fac   = factor/(fsq*ssmod);

    response->data->data[i].re = fac*ssre;
    response->data->data[i].im = fac*ssim;
  }

  response->epoch  = frameData->calibrationTime;
  response->f0     = 0;
  response->deltaF = 0.5*frameData->sampleRate/response->data->length;

  LALDestroyVector (status->statusPtr, &re);
  CHECKSTATUSPTR (status);
  LALDestroyVector (status->statusPtr, &im);
  CHECKSTATUSPTR (status);

  DETATCHSTATUSPTR (status);
  RETURN (status);
}

/*----------------------------------------------------------------------- 
 * 
 * File Name: FrameData.c
 *
 * Author: Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "FrameL.h"
#include "LALStdlib.h"
#include "LALConstants.h"
#include "SeqFactories.h"
#include "FrameData.h"

NRCSID (FRAMEDATAC, "$Id$");

void
InitializeFrameData (
    Status     *status,
    FrameData **frameData,
    CHAR       *framePath
    )
{
  const CHAR *headNames[]       = {"C1-*.F", "H-*.F", "H-*.T", "L-*.F",
                                   "L-*.T", "C1-*[0-9]"};
  const INT4  numHeadNames      = 6;
  const INT4  maxNumFiles       = 2048;
  const INT4  maxFileNameLength = 256;

  CreateVectorSequenceIn frameFileNameIn;
  CHAR                   command[1024];
  INT4                   nameType;

  INITSTATUS (status, "InitializeFrameData", FRAMEDATAC);
  ATTATCHSTATUSPTR (status);

  /* make sure arguments are reasonable */
  ASSERT (framePath, status, FRAMEDATA_ENULL, FRAMEDATA_MSGENULL);
  ASSERT (frameData, status, FRAMEDATA_ENULL, FRAMEDATA_MSGENULL);
  ASSERT (!(*frameData), status, FRAMEDATA_ENNUL, FRAMEDATA_MSGENNUL);

  /* allocate memory */
  *frameData = LALCalloc (1, sizeof(FrameData));
  ASSERT (*frameData, status, FRAMEDATA_ENULL, FRAMEDATA_MSGENULL);

  /* debuglevel zero: don't report errors */
  FrLibIni (NULL, stderr, 0);

  /* set some flags and modes */
  (*frameData)->dataBreak = 1;
  (*frameData)->newLock   = 1;
  (*frameData)->inLock    = 1;

  /* set up frame files */

  frameFileNameIn.length       = maxNumFiles;
  frameFileNameIn.vectorLength = maxFileNameLength;
  CHARCreateVectorSequence (
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
    ASSERT (nbytes > 0, status, FRAMEDATA_EREAD, FRAMEDATA_MSGEREAD);
    ASSERT (nbytes < sizeof(command), status,
            FRAMEDATA_EREAD, FRAMEDATA_MSGEREAD);

    /* fp is a stream containing the filenames */
    fp = popen (command, "r");
    ASSERT (fp, status, FRAMEDATA_EREAD, FRAMEDATA_MSGEREAD);

    /* read the filenames into the stored filename list */
    numFiles = (*frameData)->numFiles;
    while (numFiles < maxNumFiles)
    {
      CHAR *fileName;
      fileName  = (*frameData)->frameFileNames->data;
      fileName += numFiles*maxFileNameLength;
      if (EOF == fscanf (fp, "%s\n", fileName))
        break;
      ASSERT (strlen(fileName) < maxFileNameLength, status,
              FRAMEDATA_EREAD, FRAMEDATA_MSGEREAD);
      ++numFiles;
    }
    (*frameData)->numFiles = numFiles;

    pclose (fp);
  }
  ASSERT ((*frameData)->numFiles > 0, status,
          FRAMEDATA_EREAD, FRAMEDATA_MSGEREAD);
  ASSERT ((*frameData)->numFiles < maxNumFiles, status,
          FRAMEDATA_EREAD, FRAMEDATA_MSGEREAD);

  DETATCHSTATUSPTR (status);
  RETURN (status);  
}


void
FinalizeFrameData (
    Status     *status,
    FrameData **frameData
    )
{
  INITSTATUS (status, "FinalizeFrameData", FRAMEDATAC);
  ATTATCHSTATUSPTR (status);

  /* make sure argument is reasonable */
  ASSERT (frameData, status, FRAMEDATA_ENULL, FRAMEDATA_MSGENULL);
  ASSERT (*frameData, status, FRAMEDATA_ENULL, FRAMEDATA_MSGENULL);

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
  CHARDestroyVectorSequence (status->statusPtr, &(*frameData)->frameFileNames);
  CHECKSTATUSPTR (status);

  /* free memory */
  LALFree (*frameData);
  *frameData = NULL;

  DETATCHSTATUSPTR (status);
  RETURN (status);  
}


static void
GetNewFrame (
    Status    *status,
    FrameData *frameData
    )
{
  const REAL4 resolution = 3.2e-3;

  INITSTATUS (status, "GetNewFrame", FRAMEDATAC);

  /* make sure argument is not NULL */
  ASSERT (frameData, status, FRAMEDATA_ENULL, FRAMEDATA_MSGENULL);

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
    frameData->calibration =
      FrStatDataFind (frame->detectProc, "sweptsine", frame->GTimeS);
    ASSERT (frameData->calibration, status,
            FRAMEDATA_ENOSS, FRAMEDATA_MSGENOSS);

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
    frameData->dmro = dmro = FrAdcDataFind (frame, "IFO_DMRO");
    ASSERT (dmro, status, FRAMEDATA_EDMRO, FRAMEDATA_MSGEDMRO);

    frameData->numDmro = dmro->data->nData;
    frameData->curDmro = 0;

    /* get IFO_Lock */
    frameData->lock = lock = FrAdcDataFind (frame, "IFO_Lock");
    ASSERT (lock, status, FRAMEDATA_ELOCK, FRAMEDATA_MSGELOCK);

    frameData->numLock = lock->data->nData;
    frameData->curLock = 0;

    /* ratio is number of IFO_DMRO samples to each IFO_Lock sample */
    frameData->ratio = frameData->numDmro/frameData->numLock;

    /* get lock-low and lock-high values */
    frameData->lockLowHigh =
      FrStatDataFind (frame->detectProc, "locklo/lockhi", frame->GTimeS);
    ASSERT (frameData->lockLowHigh, status,
            FRAMEDATA_ELOHI, FRAMEDATA_MSGELOHI);

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


void
GetFrameData (
    Status         *status,
    INT2TimeSeries *data,
    FrameData      *frameData
    )
{
  INT4 seek       = 0;
  INT4 brokenLock = 0;
  INT4 numPoints;
  INT4 needed;

  INITSTATUS (status, "GetFrameData", FRAMEDATAC);
  ATTATCHSTATUSPTR (status);

  /* make sure arguments are reasonable */
  ASSERT (frameData, status, FRAMEDATA_ENULL, FRAMEDATA_MSGENULL);
  ASSERT (data, status, FRAMEDATA_ENULL, FRAMEDATA_MSGENULL);
  ASSERT (data->data, status, FRAMEDATA_ENULL, FRAMEDATA_MSGENULL);
  numPoints = data->data->length;
  ASSERT (numPoints, status, FRAMEDATA_ESIZE, FRAMEDATA_MSGESIZE);

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
              FRAMEDATA_EOPEN, FRAMEDATA_MSGEOPEN);
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
        INT4  wordRatio = wordSize > sizeof(INT2) ? wordSize/sizeof(INT2) : 1;

        /* copy the data if desired */
        if (!seek)
        {
          INT2 *location = data->data->data + (numPoints - needed)*wordRatio;
          memcpy (location, dmroData + frameData->curDmro, ncopy*wordSize);
        }

        /* if starting buffer: mark start time and time step */
        if (needed == numPoints)
        {
          REAL8 time;
          time  = frameData->frameStartTime.gpsSeconds;
          time += 1e-9*frameData->frameStartTime.gpsNanoSeconds;
          time += frameData->curDmro/frameData->sampleRate;
          data->epoch.gpsSeconds     = (INT4)time;
          data->epoch.gpsNanoSeconds = (INT4)fmod(1e9*time, 1e9);
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
    Status      *status,
    REAL4Vector *yout,
    REAL4Vector *yinp,
    REAL4Vector *xinp
    )
{
  REAL4Vector *yppvec = NULL;
  REAL4       *ypp;
  INT4         n;

  INITSTATUS (status, "SplineFit", FRAMEDATAC);
  ATTATCHSTATUSPTR (status);

  /* make sure arguments are reasonable */

  ASSERT (yout, status, FRAMEDATA_ENULL, FRAMEDATA_MSGENULL);
  ASSERT (yout->data, status, FRAMEDATA_ENULL, FRAMEDATA_MSGENULL);
  ASSERT (yout->length > 2, status, FRAMEDATA_ESIZE, FRAMEDATA_MSGESIZE);

  ASSERT (yinp, status, FRAMEDATA_ENULL, FRAMEDATA_MSGENULL);
  ASSERT (yinp->data, status, FRAMEDATA_ENULL, FRAMEDATA_MSGENULL);
  n = yinp->length;
  ASSERT (n > 2, status, FRAMEDATA_ESIZE, FRAMEDATA_MSGESIZE);

  ASSERT (xinp, status, FRAMEDATA_ENULL, FRAMEDATA_MSGENULL);
  ASSERT (xinp->data, status, FRAMEDATA_ENULL, FRAMEDATA_MSGENULL);
  ASSERT (xinp->length == n, status, FRAMEDATA_ESIZE, FRAMEDATA_MSGESIZE);

  /* create temporary vector */
  CreateVector (status->statusPtr, &yppvec, n);
  CHECKSTATUSPTR (status);
  ypp = yppvec->data;

  { /* setup second derivative array */
    
    REAL4Vector *uvec = NULL;
    REAL4       *u;

    REAL4 *x = xinp->data;
    REAL4 *y = yinp->data;

    INT4  i;

    CreateVector (status->statusPtr, &uvec, n);
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
    for (i = n - 2; i >= 0; --i)
    {
      ypp[i] = ypp[i]*ypp[i+1] + u[i];
    }

    DestroyVector (status->statusPtr, &uvec);
    CHECKSTATUSPTR (status);

  }

  { /* do fit */

    REAL4 *x  = xinp->data;
    REAL4 *y  = yinp->data;
    REAL4  dx = (x[n - 1] - x[0])/(yout->length - 1);

    INT4 i;
    INT4 j = 0;

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

  DestroyVector (status->statusPtr, &yppvec);
  CHECKSTATUSPTR (status);

  DETATCHSTATUSPTR (status);
  RETURN (status);  
}


void
GetFrameDataResponse (
    Status                  *status,
    COMPLEX8FrequencySeries *response,
    FrameData               *frameData
    )
{
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
  INT4   i;

  INITSTATUS (status, "GetFrameDataResponse", FRAMEDATAC);
  ATTATCHSTATUSPTR (status);

  ASSERT (frameData, status, FRAMEDATA_ENULL, FRAMEDATA_MSGENULL);
  ASSERT (response, status, FRAMEDATA_ENULL, FRAMEDATA_MSGENULL);
  ASSERT (response->data, status, FRAMEDATA_ENULL, FRAMEDATA_MSGENULL);
  ASSERT (response->data->data, status, FRAMEDATA_ENULL, FRAMEDATA_MSGENULL);
  ASSERT (response->data->length > 2, status,
          FRAMEDATA_ESIZE, FRAMEDATA_MSGESIZE);

  fri = ((struct FrStatData *)(frameData->calibration))->data->dataF;
  n   = ((struct FrStatData *)(frameData->calibration))->data->nData;
  ssn = n/3;
  ssf = fri;
  ssr = ssf + ssn;
  ssi = ssr + ssn;
  ASSERT (ssn > 2, status, FRAMEDATA_ESSSZ, FRAMEDATA_MSGESSSZ);

  CreateVector (status->statusPtr, &re, response->data->length);
  CHECKSTATUSPTR (status);
  CreateVector (status->statusPtr, &im, response->data->length);
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
  response->deltaF = frameData->sampleRate;

  DestroyVector (status->statusPtr, &re);
  CHECKSTATUSPTR (status);
  DestroyVector (status->statusPtr, &im);
  CHECKSTATUSPTR (status);

  DETATCHSTATUSPTR (status);
  RETURN (status);  
}

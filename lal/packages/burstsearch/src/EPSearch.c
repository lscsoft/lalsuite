/******** <lalVerbatim file="EPSearchCV"> ********
Author: Brady, P
$Id$
********* </lalVerbatim> ********/

#include <lal/LALRCSID.h>

NRCSID (EPSEARCHC, "$Id$");


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/SeqFactories.h>
#include <lal/RealFFT.h>
#include <lal/Thresholds.h>
#include <lal/ExcessPower.h>
#include <lal/BurstSearch.h>
#include <lal/EPData.h>
#include <lal/Random.h>


#define TRUE 1
#define FALSE 0

extern INT4 lalDebugLevel;

/******** <lalVerbatim file="EPSearchCP"> ********/
void
EPSearch (
		   LALStatus               *status,
                   EPSearchParams          *params,
                   BurstEvent             **burstEvent,
                   UINT4                    tmpDutyCycle
		   )
/******** </lalVerbatim> ********/
{ 
    INT4                      i,j;
    REAL4                     redummy, imdummy;
    EPDataSegment            *dummySegment     = NULL;
    BurstEvent               *currentEvent     = NULL;
    BurstEvent               *prevEvent        = NULL;
    COMPLEX8FrequencySeries  *fseries;
    RealDFTParams            *dftparams        = NULL;
    LALWindowParams           winParams;
    REAL4                    *dummySpec        = NULL;

    INITSTATUS (status, "EPSearch", EPSEARCHC);
    ATTATCHSTATUSPTR (status);

    /* make sure that arguments are not NULL */
    ASSERT (params, status, EXCESSPOWERH_ENULLP, EXCESSPOWERH_MSGENULLP);
    ASSERT (burstEvent, status, EXCESSPOWERH_ENULLP, EXCESSPOWERH_MSGENULLP);

    /* Set up the window parameters */
    winParams.type=params->winParams.type;
    winParams.length=params->ntotT;

    /* assign temporary memory for the frequency data */
    fseries = (COMPLEX8FrequencySeries *) LALMalloc (sizeof(COMPLEX8FrequencySeries));
    strncpy( fseries->name, "anonymous", LALNameLength * sizeof(CHAR) );
    fseries->data = NULL;
    LALCCreateVector (status->statusPtr, &fseries->data, 
        params->initParams->numPoints/2 + 1);
    CHECKSTATUSPTR (status);

    /* create the dft params */
    LALCreateRealDFTParams(status->statusPtr , &dftparams, &winParams, 1);
    CHECKSTATUSPTR (status);

    /* point to the start of event list */
    params->numEvents=0;

    /* allocate temporary memory for spectrum */
    dummySpec = (REAL4 *)LALMalloc(fseries->data->length * sizeof(REAL4));
    for (j=0 ; j<(INT4)fseries->data->length ; j++)
    {
      dummySpec[j] = 0.0;
    }

    if (params->initParams->method == useMedian)
    /* run new code using median to obtain power spectrum */
    { /* loop over data computing spectrum */
      /* allocate two dimensional array to hold power at each freq in each time slice */
      INT4 flength = fseries->data->length;
      INT4 numSegs = (INT4)tmpDutyCycle;
      REAL4 *ptr = (REAL4 *)LALMalloc(flength * numSegs * sizeof(REAL4));
      ASSERT(ptr, status, EPSEARCHH_EMALLOC, EPSEARCHH_MSGEMALLOC);

      /* zero out memory array */
      for (j = 0; j < flength; j++)
      for (i = 0; i < numSegs; i++)
        {ptr[i * flength + j] = 0.0;}

      /* obtain power spectrum in each time slice */
      for (i = 0; i < numSegs; i++)
      {
        /* point dummySegment to the segment to analyze */
        dummySegment = params->epSegVec->data + params->currentSegment + i;

        /* compute the DFT of input time series */
        LALComputeFrequencySeries (status->statusPtr, fseries, 
            dummySegment->data, dftparams);
        CHECKSTATUSPTR (status);

        /* copy modulus of DFT output into two dimensional array */
        for (j=0 ; j < flength ; j++)
        {
          redummy = fseries->data->data[j].re ;
          imdummy = fseries->data->data[j].im ;
          ptr[i * flength + j] = redummy*redummy + imdummy*imdummy;
        }
      }  /* done computing spectrum */

      /* find median value over time slices for each frequency */
      for (j = 0; j < flength; j++)
      {
        dummySpec[j] = EPMedian(ptr, j, flength, numSegs, status);
        dummySpec[j] *= sqrt(2.0); /* scale to match mean method */
      }

      LALFree(ptr);     
    }  /* end of new code using median */
    
    else if (params->initParams->method == useMean)
    /* run original code using mean to obtain power spectrum */
    {  /* loop over data computing spectrum */
      for ( i=0 ; i<(INT4)tmpDutyCycle ; i++)
      {
        /* point dummySegment to the segment to analyze */
        dummySegment = params->epSegVec->data + params->currentSegment + i;

        /* compute the DFT of input time series */
        LALComputeFrequencySeries (status->statusPtr, fseries, 
            dummySegment->data, dftparams);
        CHECKSTATUSPTR (status);

        /* normalize the data stream so that rms of Re or Im is 1 */
        redummy=imdummy=0.0;
        for (j=0 ; j<(INT4)fseries->data->length ; j++)
        {
          redummy = fseries->data->data[j].re ;
          imdummy = fseries->data->data[j].im ;
          dummySpec[j] += redummy*redummy + imdummy*imdummy;
        }
      }
      for (j=0 ; j<(INT4)fseries->data->length ; j++)
      {
        dummySpec[j] /= ((REAL4) tmpDutyCycle);
      }
    }  /* end of original code using mean */
    
    else if (params->initParams->method == useUnity)
    /* force power spectrum to unity */
    {
      for (j=0 ; j<(INT4)fseries->data->length ; j++)
        {dummySpec[j] = 1.0;}
    }

    /* default case for unknown method */
    else
    {
        ABORT (status, EPSEARCHH_EINCOMP, EPSEARCHH_MSGEINCOMP );
    }
    
    /* write diagnostic info to disk */
    { /* output normalized power spectrum MSW 7/19/02 */
      FILE *fp;
      fp = fopen("./freqseries.dat","w");
      for (j=0 ; j<(INT4)fseries->data->length ; j++)
      {
        fprintf(fp, "%i\t%g\n", j, dummySpec[j]);
      }    
      fclose(fp);
    }
    
    /* loop over data applying excess power method */
    for ( i=0 ; i<(INT4)tmpDutyCycle ; i++)
    {

      /*
       * 
       * determine the type of run we're doing:  
       *                                         0. Analyze data;
       *                                         1. Gaussian Sim;  
       *                                         2. Injection;
       *                                         
       *
       */
      
      if ( params->simType == 1 ) 
      {
        /* point dummySegment to the segment to analyze */
        dummySegment = params->epSegVec->data;
      }
      else
      {
        /* point dummySegment to the segment to analyze */
        dummySegment = params->epSegVec->data + params->currentSegment + i;

        /* if we're doing simulated injections */
        if ( params->simType == 2 ) {

        }
      }

      /* compute the DFT of input time series */
      LALComputeFrequencySeries (status->statusPtr, fseries, 
          dummySegment->data, dftparams);
      CHECKSTATUSPTR (status);

      /* check that deltaF agrees with that of response */
      if ( fabs( dummySegment->spec->deltaF - fseries->deltaF ) > 0.000001 )
      {
        ABORT (status, EXCESSPOWERH_EDELF, EXCESSPOWERH_MSGEDELF );
      }

      /* normalize the data stream so that rms of Re or Im is 1 */
      redummy=imdummy=0.0;
      for (j=0 ; j<(INT4)fseries->data->length ; j++)
      {
        REAL4 tmpVar = sqrt( 4 * dummySegment->data->deltaT / 
            dummySegment->spec->data->data[j] );
        tmpVar = sqrt( 2.0 / dummySpec[j] );
        fseries->data->data[j].re *= tmpVar;
        fseries->data->data[j].im *= tmpVar;
        redummy += fseries->data->data[j].re * fseries->data->data[j].re;
        imdummy += fseries->data->data[j].im * fseries->data->data[j].im;
      }
      
#ifdef PRINTAVERAGE
      {
        FILE *fp;
        fp = fopen("./power.dat","a");
        fprintf (fp , "Real average = %f, Imag average = %f\n",redummy/((REAL4)
              fseries->data->length), imdummy/((REAL4) fseries->data->length));
        fclose(fp);
      }
#endif


        

      /* create time-frequency tiling of plane.  */
      if ( params->tfTiling == NULL ){
        params->tfTilingInput->deltaF=fseries->deltaF;
        LALCreateTFTiling (status->statusPtr, &(params->tfTiling), params->tfTilingInput);
        CHECKSTATUSPTR (status);
      }

      /* compute the TFplanes for the data segment */
      LALComputeTFPlanes (status->statusPtr, params->tfTiling, fseries);
      CHECKSTATUSPTR (status);

      /* search these planes */
      LALComputeExcessPower (status->statusPtr, params->tfTiling, params->compEPInput);
      CHECKSTATUSPTR (status);

      /* compute the likelihood for slightly better detection method */
      /*
       * LALComputeLikelihood  (status->statusPtr, &(params->lambda), params->tfTiling);
       * CHECKSTATUSPTR (status);
       */

      /* sort the results. */
      LALSortTFTiling (status->statusPtr, params->tfTiling);
      CHECKSTATUSPTR (status);


      /* count the number of events  */
      /* change alphaThreshold to match with confidence */
      /*
       * LALCountEPEvents(status->statusPtr, &(params->numEvents), 
       * params->tfTiling, params->alphaThreshold); 
       * CHECKSTATUSPTR (status);
       */

      {
        TFTile *thisTile = params->tfTiling->firstTile;
        INT4 tileCount   = 0;

        while ( (thisTile != NULL) && (thisTile->alpha <= params->alphaThreshold) 
            && (tileCount < params->events2Master) )
        {
          INT8 tstartNS = 0;
          
          /* increment local and global counter */
          tileCount++;
          (params->numEvents)++;

          /* convert epoch to GPS nanoseconds */
          tstartNS  = 1000000000L * 
            (INT8) dummySegment->data->epoch.gpsSeconds;
          tstartNS += (INT8) dummySegment->data->epoch.gpsNanoSeconds;

          /* allocate memory for the burst event */
          if ( (*burstEvent) == NULL )
          {
            currentEvent=(*burstEvent)=(BurstEvent *) LALMalloc( sizeof(BurstEvent) );
          }
          else 
          {
            currentEvent = (BurstEvent *) LALMalloc( sizeof(BurstEvent) );
          }
          
          /* build a burst event from TFTile */
          LALTFTileToBurstEvent(status->statusPtr, currentEvent, thisTile,
              tstartNS, params->tfTilingInput->flow); 
          CHECKSTATUSPTR (status);

          /* point to the next event */
          currentEvent->nextEvent = NULL;
          if (prevEvent != NULL) prevEvent->nextEvent = currentEvent;
          prevEvent = currentEvent;
          currentEvent = currentEvent->nextEvent;
          thisTile = thisTile->nextTile;
        }
      }

      /* reset the flags on the tftiles */
      params->tfTiling->planesComputed=FALSE;
      params->tfTiling->excessPowerComputed=FALSE;
      params->tfTiling->tilesSorted=FALSE;

    }

    /* destroy time-frequency tiling of planes */
    LALDestroyTFTiling (status->statusPtr, &(params->tfTiling));
    CHECKSTATUSPTR (status);

    /* destroy temporary spectrum */
    LALFree(dummySpec);

    /* destroy the dftparams for computing frequency series */
    LALDestroyRealDFTParams (status->statusPtr, &dftparams);
    CHECKSTATUSPTR (status);

    /* destroy temporary memory for the frequency data */
    LALCDestroyVector (status->statusPtr, &fseries->data);
    CHECKSTATUSPTR (status);
    LALFree(fseries);

    /* normal exit */
    DETATCHSTATUSPTR (status);
    RETURN (status);
}

/*************************************************************/
/* calculation of median for power spectrum data 7/17/02 MSW */
/* implements an insert sort using temporary array           */
/* assumes inputs are positive, sorts in descending order    */
/*************************************************************/

REAL4 EPMedian(REAL4 *p, INT4 j, INT4 flength, INT4 numSegs, LALStatus *status)
{
    /* p points to array of power spectra data over time slices */
    /* j is desired frequency offset into power spectra array   */
    /* flength is size of frequency series obtained from DFT    */
    /* numSegs is the number of time slices to be evaluated     */
    /* status points to LAL status struct passed into main      */
    /* returns the median value, over time slice at given freq. */

    /* TODO should check for valid or reasonable inputs */
    INT4  outer  = 0;       /* local loop counter */
    INT4  middle = 0;       /* local loop counter */
    INT4  inner  = 0;       /* local loop counter */
    REAL4 returnVal = 0.0;  /* holder for return value */

    /* allocate memory array for insert sort, test for success */
    REAL4 *s = (REAL4 *)LALMalloc(numSegs * sizeof(REAL4));
    ASSERT(s, status, EPSEARCHH_EMALLOC, EPSEARCHH_MSGEMALLOC);

    /* zero out the sort array */
    for (outer = 0; outer < numSegs; outer++)
    {
      s[outer] = 0.0;
    }

    /* scan time slices for a given frequency */
    for (outer = 0; outer < numSegs; outer++)
    {
      /* insert power value into sort array */
      REAL4 tmp = p[outer * flength + j]; /* obtain value to insert */
      for (middle = 0; middle < numSegs; middle++)
      {
        if (tmp > s[middle])
        {
          /* insert taking place of s[middle] */
          for (inner = numSegs - 1; inner > middle; inner--)
          {
            s[inner] = s [inner - 1];  /* move old values */
          }
          s[middle] = tmp;   /* insert new value */
          break;  /* terminate loop */
        }
      }
    }  /* done inserting into sort array */

    /* check for odd or even number of segments */
    if (numSegs % 2)
    {
      /* if odd number of segments, return median */
      returnVal = s[numSegs / 2];
    }
    else
    {
      /* if even number of segments, return average of two medians */
      returnVal = 0.5 * (s[numSegs/2] + s[(numSegs/2) - 1]);
    }

    /* free memory used for sort array */
    LALFree(s);

    return returnVal;
}


/******** <lalVerbatim file="EPSearchCV"> ********
Author: Brady, P
$Id$
********* </lalVerbatim> ********/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/SeqFactories.h>
#include <lal/RealFFT.h>
#include <lal/Thresholds.h>
#include <lal/ExcessPower.h>
#include <lal/BurstSearch.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/EPData.h>
#include <lal/Random.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LALRCSID.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/PrintFTSeries.h>

NRCSID (EPSEARCHC, "$Id$");

#define TRUE 1
#define FALSE 0

extern INT4 lalDebugLevel;

static void EPMedian(REAL4 *ans, REAL4 *p, INT4 j, INT4 flength, 
    INT4 numSegs, LALStatus *status);
static void EPClean(INT4 freqcl,  COMPLEX8FrequencySeries *freqSeries);

/****************************************************************
 *
 * Weights tiles according to the number present of a given
 * time-frequency volume.
 *
 ***************************************************************/

/******** <lalVerbatim file="LALWeighTFTileListCP"> ********/
void LALWeighTFTileList (
        LALStatus         *status,
        TFTiling          *tfTiling,
        INT4               maxDOF
        )
/******** </lalVerbatim> ********/
{
  TFTile *thisTile;
  INT4   *tmparray;
  INT4    i=0;

  INITSTATUS (status, "LALPrintTFTileList", EPSEARCHC);
  ATTATCHSTATUSPTR (status);

  ASSERT(tfTiling, status, EXCESSPOWERH_ENULLP, EXCESSPOWERH_MSGENULLP); 
  ASSERT(tfTiling->firstTile, status, EXCESSPOWERH_ENULLP, \
         EXCESSPOWERH_MSGENULLP); 

  
  tmparray = (INT4 *) LALMalloc (2*maxDOF*sizeof(INT4));
  for (i=0 ; i < 2*maxDOF ; i++){
        tmparray[i]=0;
  }

  thisTile = tfTiling->firstTile;
  while ( (thisTile != NULL) )
  {
      INT4 t1 = thisTile->tstart;
      INT4 t2 = thisTile->tend;
      INT4 f1 = thisTile->fstart;
      INT4 f2 = thisTile->fend;
      INT4 dof = 2*(t2-t1+1)*(f2-f1+1);
      tmparray[dof]+=1;

      thisTile = thisTile->nextTile;
  }
  thisTile = tfTiling->firstTile;
  while ( (thisTile != NULL) )
  {
      INT4 t1 = thisTile->tstart;
      INT4 t2 = thisTile->tend;
      INT4 f1 = thisTile->fstart;
      INT4 f2 = thisTile->fend;
      INT4 dof = 2*(t2-t1+1)*(f2-f1+1);

      thisTile->weight=(REAL4)tmparray[dof];
      thisTile = thisTile->nextTile;
  }

  LALFree(tmparray);

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}

static REAL4 square_modulus(COMPLEX8 z)
{
	return(z.re*z.re + z.im*z.im);
}

static void ComputeAverageSpectrum(
	LALStatus *status,
	REAL4 *AverageSpec,
        EPSearchParams *params,
	INT4 numSegs,
	EPDataSegment *segment
)
{
    int i, j;
    COMPLEX8FrequencySeries fseries;
    RealDFTParams *dftparams = NULL;
    LALWindowParams winParams;
    INT4 numPoints = params->initParams->numPoints;
    AvgSpecMethod method = params->initParams->method;

    INITSTATUS (status, "ComputeAverageSpectrum", EPSEARCHC);
    ATTATCHSTATUSPTR (status);

    /* create the dft params */
    winParams.type = params->winParams.type;
    winParams.length = 2 * params->ntotT;
    LALCreateRealDFTParams(status->statusPtr , &dftparams, &winParams, 1);
    CHECKSTATUSPTR(status);

    /* Initialize the frequency series container */
    strncpy(fseries.name, "anonymous", LALNameLength * sizeof(CHAR));
    fseries.data = NULL;
    LALCCreateVector(status->statusPtr, &fseries.data, numPoints + 1);
    CHECKSTATUSPTR(status);

    if (useMedian)
    /* run new code using median to obtain power spectrum */
    { /* loop over data computing spectrum */
      /* allocate two dimensional array to hold power at each freq in each time slice */
      INT4 flength = fseries.data->length;
      REAL4 *ptr = LALMalloc(flength * numSegs * sizeof(*ptr));
      ASSERT(ptr, status, EPSEARCHH_EMALLOC, EPSEARCHH_MSGEMALLOC);

      /* obtain power spectrum in each time slice */
      for (i = 0; i < numSegs; i++)
      {
        /* compute the DFT of input time series */
        LALComputeFrequencySeries (status->statusPtr, &fseries, (segment + i)->data, dftparams);
        CHECKSTATUSPTR (status);

        /* copy modulus of DFT output into two dimensional array */
        for (j=0 ; j < flength ; j++)
          ptr[i * flength + j] = square_modulus(fseries.data->data[j]);
      }  /* done computing spectrum */

      /* find median value over time slices for each frequency */
      for (j = 0; j < flength; j++)
      {
        EPMedian(&AverageSpec[j], ptr, j, flength, numSegs, status);
        AverageSpec[j] /= LAL_LN2; /* scale to match mean method */
      }

      LALFree(ptr);     
    }  /* end of new code using median */
    
    else if (method == useMean)
    /* run original code using mean to obtain power spectrum */
    {  /* loop over data computing spectrum */
      for ( i=0 ; i < numSegs; i++)
      {
        /* compute the DFT of input time series */
        LALComputeFrequencySeries (status->statusPtr, &fseries, (segment + i)->data, dftparams);
        CHECKSTATUSPTR (status);

        /* normalize the data stream so that rms of Re or Im is 1 */
        for (j=0 ; j < (INT4)fseries.data->length ; j++)
          AverageSpec[j] += square_modulus(fseries.data->data[j]);
      }
      for (j=0 ; j < (INT4)fseries.data->length ; j++)
      {
        AverageSpec[j] /= numSegs;
      }
    }  /* end of original code using mean */
    
    else if (method == useUnity)
    /* force power spectrum to unity */
    {
      for (j=0 ; j < (INT4)fseries.data->length ; j++)
        AverageSpec[j] = 1.0;
    }

    /* default case for unknown method */
    else
      ABORT (status, EPSEARCHH_EINCOMP, EPSEARCHH_MSGEINCOMP );

    LALCDestroyVector(status->statusPtr, &fseries.data);
    CHECKSTATUSPTR (status);
    LALDestroyRealDFTParams (status->statusPtr, &dftparams);
    CHECKSTATUSPTR (status);

    DETATCHSTATUSPTR( status );
    RETURN( status );
}


/******** <lalVerbatim file="EPSearchCP"> ********/
void
EPSearch (
        LALStatus               *status,
        REAL4TimeSeries         *tseries,
        EPSearchParams          *params,
        SnglBurstTable         **burstEvent,
        UINT4                    tmpDutyCycle
        )
/******** </lalVerbatim> ********/
{ 
    INT4                      i,j,freqcl;
    EPDataSegment            *segment;
    SnglBurstTable           *currentEvent     = NULL;
    SnglBurstTable           *prevEvent        = NULL;
    COMPLEX8FrequencySeries  *fseries;
    RealDFTParams            *dftparams        = NULL;
    LALWindowParams           winParams;
    AverageSpectrumParams     AverageSpecParams;
    REAL4FrequencySeries     *AverageSpec;
    INT4                      nevents, dumevents;

    INITSTATUS (status, "EPSearch", EPSEARCHC);
    ATTATCHSTATUSPTR (status);

    /* make sure that arguments are not NULL */
    ASSERT (params, status, EXCESSPOWERH_ENULLP, EXCESSPOWERH_MSGENULLP);
    ASSERT (burstEvent, status, EXCESSPOWERH_ENULLP, EXCESSPOWERH_MSGENULLP);

    /* assign temporary memory for the frequency data */
    fseries = (COMPLEX8FrequencySeries *) LALMalloc (sizeof(COMPLEX8FrequencySeries));
    strncpy( fseries->name, "anonymous", LALNameLength * sizeof(CHAR) );
    fseries->data = NULL;
    LALCCreateVector (status->statusPtr, &fseries->data, params->initParams->numPoints + 1);
    CHECKSTATUSPTR (status);

    /* Set up the window parameters */
    winParams.type=params->winParams.type;
    winParams.length=2 * params->ntotT;

    /* create the dft params */
    LALCreateRealDFTParams(status->statusPtr , &dftparams, &winParams, 1);
    CHECKSTATUSPTR (status);

    /* point to the start of event list */
    params->numEvents=0;

    /* Compute the average spectrum */
    LALNewREAL4FrequencySeries(status->statusPtr, &AverageSpec, "", LIGOTIMEGPSINITIALIZER, 0, 0, LALUNITINITIALIZER, fseries->data->length);
    CHECKSTATUSPTR(status);

    ComputeAverageSpectrum(status->statusPtr, AverageSpec->data->data, params, tmpDutyCycle, params->epSegVec->data + params->currentSegment);
    CHECKSTATUSPTR(status);

    /* write diagnostic info to disk */
    if ( params->printSpectrum == TRUE )
    { 
        FILE *fp;
        fp = fopen("./freqseries.dat","w");
        for (j=0 ; j<(INT4)fseries->data->length ; j++)
        {
            fprintf(fp, "%f\t%g\n", j*fseries->deltaF, AverageSpec->data->data[j]);
        }    
        fclose(fp);
    }

    /* loop over data applying excess power method */
    for ( i=0 ; i<(INT4)tmpDutyCycle ; i++)
    {

      /* point segment to the segment to analyze */
      segment = params->epSegVec->data + params->currentSegment + i;

      /* compute the DFT of input time series */
      LALInfo(status->statusPtr, "Computing the frequency series");
      CHECKSTATUSPTR (status);
      LALComputeFrequencySeries (status->statusPtr, fseries, 
          segment->data, dftparams);
      CHECKSTATUSPTR (status);

      /* check that deltaF agrees with that of response */
      if ( fabs( segment->spec->deltaF - fseries->deltaF ) > 0.000001 )
      {
        ABORT (status, EXCESSPOWERH_EDELF, EXCESSPOWERH_MSGEDELF );
      }

      /* normalize the data stream so that rms of Re or Im is 1 */
      for (j=0 ; j<(INT4)fseries->data->length ; j++)
      {
        REAL4 tmpVar;
        tmpVar = sqrt( 2.0 / AverageSpec->data->data[j] );
        fseries->data->data[j].re *= tmpVar;
        fseries->data->data[j].im *= tmpVar;
      }

      /* write diagnostic info to disk */
      if ( params->printSpectrum == TRUE )
      { 
        FILE *fp;
        fp = fopen("./dummy.dat","w");
        for (j=0 ; j<(INT4)fseries->data->length ; j++)
        {
          fprintf(fp, "%f\t%g\n", j*fseries->deltaF, 
              sqrt(fseries->data->data[j].re * fseries->data->data[j].re
                + fseries->data->data[j].im * fseries->data->data[j].im));
        }    
        fclose(fp);
      }

      /* create time-frequency tiling of plane.  */
      if ( params->tfTiling == NULL ){
        /* the factor of 2 here is to account for the overlapping */
        params->tfTilingInput->deltaF=2.0*fseries->deltaF;
        LALCreateTFTiling (status->statusPtr, &(params->tfTiling), 
            params->tfTilingInput);
        CHECKSTATUSPTR (status);
      }

      /* compute the TFplanes for the data segment */
      LALInfo(status->statusPtr, "Computing the TFPlanes");
      CHECKSTATUSPTR (status);
      LALComputeTFPlanes (status->statusPtr, params->tfTiling, fseries);
      CHECKSTATUSPTR (status);

      /* search these planes */
      LALInfo(status->statusPtr, "Computing the excess power");
      CHECKSTATUSPTR (status);
      LALComputeExcessPower (status->statusPtr, params->tfTiling, 
          params->compEPInput);
      CHECKSTATUSPTR (status);

      /* compute the likelihood for slightly better detection method */
      /*
       * LALComputeLikelihood  (status->statusPtr, &(params->lambda), params->tfTiling);
       * CHECKSTATUSPTR (status);
       */

      /* sort the results. */
      LALInfo(status->statusPtr, "Sorting TFTiling");
      CHECKSTATUSPTR (status);
      LALSortTFTiling (status->statusPtr, params->tfTiling);
      CHECKSTATUSPTR (status);

      /* determine the weighting for each tile */
      LALWeighTFTileList ( status->statusPtr, params->tfTiling, 10000);
      CHECKSTATUSPTR (status);

      /* convert the TFTiles into sngl_burst events for output */
      LALInfo(status->statusPtr, "Converting times into sngl_burst events");
      CHECKSTATUSPTR (status);
      {
        TFTile *thisTile = params->tfTiling->firstTile;
        INT4 tileCount   = 0;

        while ( (thisTile != NULL) 
                && (thisTile->alpha <= params->alphaThreshold/thisTile->weight) 
                && (tileCount < params->events2Master) )
        {
          INT8 tstartNS = 0;
          
          /* increment local and global counter */
          tileCount++;
          (params->numEvents)++;

          /* convert epoch to GPS nanoseconds */
          tstartNS  = 1000000000L * 
            (INT8) segment->data->epoch.gpsSeconds;
          tstartNS += (INT8) segment->data->epoch.gpsNanoSeconds;

          /* allocate memory for the burst event */
          if ( (*burstEvent) == NULL )
          {
            currentEvent=(*burstEvent)=(SnglBurstTable *) 
              LALMalloc( sizeof(SnglBurstTable) );
          }
          else 
          {
            currentEvent = (SnglBurstTable *) 
              LALMalloc( sizeof(SnglBurstTable) );
          }

          /* build a burst event from TFTile */
          LALTFTileToBurstEvent(status->statusPtr, currentEvent, thisTile,
              tstartNS, params); 
          CHECKSTATUSPTR (status);

          /* point to the next event */
          currentEvent->next = NULL;
          if (prevEvent != NULL) prevEvent->next = currentEvent;
          prevEvent = currentEvent;
          currentEvent = currentEvent->next;
          thisTile = thisTile->nextTile;
        }
      }

      /* reset the flags on the tftiles */
      params->tfTiling->planesComputed=FALSE;
      params->tfTiling->excessPowerComputed=FALSE;
      params->tfTiling->tilesSorted=FALSE;
    }
  
    nevents = 0;
    dumevents = 1;
    j = 0;
    /* cluster the events if requested */
    if (params->cluster && (*burstEvent) != NULL)
      {
	while ( (dumevents != nevents) && j < 500) 
	{
	  dumevents = nevents;
	  LALSortSnglBurst(status->statusPtr, burstEvent, 
              LALCompareSnglBurstByTimeAndFreq);
	  CHECKSTATUSPTR (status);
	  LALClusterSnglBurstTable(status->statusPtr, *burstEvent, 
              &nevents);
	  CHECKSTATUSPTR (status);
	  j++;
	}
      }
    
    /* destroy time-frequency tiling of planes */
    LALDestroyTFTiling (status->statusPtr, &(params->tfTiling));
    CHECKSTATUSPTR (status);

    /* destroy temporary spectrum */
    LALDestroyREAL4FrequencySeries(status->statusPtr, AverageSpec);

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

static void EPMedian(REAL4 *ans, REAL4 *p, INT4 j, INT4 flength, 
    INT4 numSegs, LALStatus *status)
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
    *ans = 0;
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

    *ans = returnVal;
    return;
}


/***************************************************************
 *
 * Initialize parameters and memory for the EP search
 *
 ***************************************************************/

/* <lalVerbatim file="EPInitSearchCP"> */
void EPInitSearch(
        LALStatus             *status,
        void                 **searchParams,
        CHAR                  *argv[],
        INT4                   argc
        )
/* </lalVerbatim> */
{

    EPSearchParams *params;


    INITSTATUS (status, "LALInitSearchPRB", EPSEARCHC);
    ATTATCHSTATUSPTR (status);

    if ( !(argv) ){
        ABORT( status, EPSEARCHH_ENULLP, EPSEARCHH_MSGENULLP);
    }
    if ( argc < 19 ){
        ABORT( status, EPSEARCHH_EARGS, EPSEARCHH_MSGEARGS ); 
    }

    if (strcmp( argv[0], "-filterparams" )){
        ABORT( status, EPSEARCHH_EARGS, EPSEARCHH_MSGEARGS);
    }
    if ( atoi( argv[1] ) <= 0 ){
        ABORT( status, EPSEARCHH_ENUMZ, EPSEARCHH_MSGENUMZ);
    }
    if ( atoi( argv[2] ) <= 0 ){
        ABORT( status, EPSEARCHH_ESEGZ, EPSEARCHH_MSGESEGZ);
    }
    if ( atoi( argv[3] ) < 0 ){
        ABORT( status, EPSEARCHH_EOVLP, EPSEARCHH_MSGEOVLP);
    }
    if ( atoi( argv[4] ) < 0 ){
        ABORT( status, EPSEARCHH_EOVPF, EPSEARCHH_MSGEOVPF);
    }
    if ( atoi( argv[5] ) <= 0 ){
        ABORT( status, EPSEARCHH_EMFBZ, EPSEARCHH_MSGEMFBZ);
    }
    if ( atoi( argv[6] ) <= 0 ){
        ABORT( status, EPSEARCHH_EMTBZ, EPSEARCHH_MSGEMTBZ);
    }
    if ( atof( argv[7] ) <= 0.0 ){
        ABORT( status, EPSEARCHH_EFLOW, EPSEARCHH_MSGEFLOW);
    }
    if ( atof( argv[8] ) <= 0.0 ){
        ABORT( status, EPSEARCHH_EDELF, EPSEARCHH_MSGEDELF);
    }
    if ( atoi( argv[9] ) <= 0 ){
        ABORT( status, EPSEARCHH_ELTFZ, EPSEARCHH_MSGELTFZ);
    }
    if ( atof( argv[10] ) <= 1.0 ){
        ABORT( status, EPSEARCHH_ESIGM, EPSEARCHH_MSGESIGM);
    }
    if ( atof( argv[11] ) <= 0.0 || 
            atof( argv[11] ) >= 1.0 ){
        ABORT( status, EPSEARCHH_EALPH, EPSEARCHH_MSGEALPH);
    }
    if ( atoi( argv[12] ) < 1 ){
        ABORT( status, EPSEARCHH_EDUTY, EPSEARCHH_MSGEDUTY);
    }
    if ( atof( argv[13] ) < 0 ){
        ABORT( status, EPSEARCHH_EAMAX, EPSEARCHH_MSGEAMAX);
    }
    /* EK - modified, argv[14] is the number of communicated events (integer) */
    if ( atoi( argv[14] ) < 1 ||
            atoi( argv[14] ) > 999){
        ABORT( status, EPSEARCHH_EE2MS, EPSEARCHH_MSGEE2MS);
    }
    /* EK - modified, argv[15] is the string representing the channel to apply */
    if ( strlen( argv[15] ) == 0 ){
        ABORT( status, EPSEARCHH_ECHNL, EPSEARCHH_MSGECHNL);
    }
    /* Simulation type */
    if ( atoi( argv[16] ) < 0 ||
            atoi( argv[16] ) > 3 ){
        ABORT( status, EPSEARCHH_ESIM, EPSEARCHH_MSGESIM);
    }
    /* Spectrum estimator to use when whitening the data */
    if ( strlen( argv[17] ) ==0 ){
        ABORT( status, EPSEARCHH_ESPEC, EPSEARCHH_MSGESPEC);
    }
    /* Type of window to use on the data */
    if ( atoi( argv[18] ) < 0 ||
            atoi( argv[18] ) > 6 ){
        ABORT( status, EPSEARCHH_EWIN, EPSEARCHH_MSGEWIN);
    }

    /*
     *
     * allocate memory 
     *
     */

    *searchParams = params = LALMalloc (sizeof( EPSearchParams )); 
    if ( !params )
    {
        ABORT (status, EPSEARCHH_EALOC, EPSEARCHH_MSGEALOC);
    }

    params->tfTilingInput = 
        (CreateTFTilingIn *) LALMalloc (sizeof(CreateTFTilingIn));
    if ( !params->tfTilingInput )
    {
        LALFree (*searchParams); *searchParams = NULL;
        ABORT (status, EPSEARCHH_EALOC, EPSEARCHH_MSGEALOC);
    }

    params->initParams = (EPInitParams *) LALMalloc (sizeof(EPInitParams));
    if ( !params-> initParams)
    {
        LALFree (params->tfTilingInput); params->tfTilingInput = NULL;
        LALFree (*searchParams); *searchParams = NULL;
        ABORT (status, EPSEARCHH_EALOC, EPSEARCHH_MSGEALOC);
    }

    params->compEPInput = 
        (ComputeExcessPowerIn *) LALMalloc (sizeof(ComputeExcessPowerIn));
    if ( !params-> compEPInput)
    {
        LALFree (params->initParams); params->initParams = NULL;
        LALFree (params->tfTilingInput); params->tfTilingInput = NULL;
        LALFree (*searchParams); *searchParams = NULL;
        ABORT (status, EPSEARCHH_EALOC, EPSEARCHH_MSGEALOC);
    }

    /* EK - Channel name to process (e.g. "ifodmro") */
    params->channelName = (CHAR *) LALMalloc(strlen( argv[15] )+1 );
    if ( !params->channelName )
    {
        LALFree (params->compEPInput); params->compEPInput = NULL;
        LALFree (params->initParams); params->initParams = NULL;
        LALFree (params->tfTilingInput); params->tfTilingInput = NULL;
        ABORT (status, EPSEARCHH_EALOC, EPSEARCHH_MSGEALOC);
    }     

    /*
     * 
     * assign the initial parameter values
     *
     */

    /* Number of data points in a segment */
    params->initParams->numPoints         = atoi( argv[1] );
    /* Number of overlapping data segments */
    params->initParams->numSegments       = atoi( argv[2] );
    /* Number of segments sent to slave */
    params->initParams->segDutyCycle      = atoi( argv[12] ); 
    /* Overlap betweeen segments (# of points) */
    params->ovrlap                        = atoi( argv[3] );  
    /* Identify events with alpha less that this value */
    params->alphaThreshold                = atof( argv[13] ); 
    /* Number of data points in a segment */
    params->ntotT                         = params->initParams->numPoints;      

    /* Amount of overlap between neighboring TF tiles */
    params->tfTilingInput->overlapFactor  = atoi( argv[4] );  
    /* Smallest extent in freq of TF tiles to search */
    params->tfTilingInput->minFreqBins    = atoi( argv[5] );  
    /* Smallest extent in time of TF tiles to search */
    params->tfTilingInput->minTimeBins    = atoi( argv[6] );  
    /* Lowest frequency in Hz to be searched */
    params->tfTilingInput->flow           = atof( argv[7] );  
    /* Frequency resolution of first TF plane */
    params->tfTilingInput->deltaF         = atof( argv[8] );  
    /* Length (N_F) of first TF plane (with N_T = 1) */
    params->tfTilingInput->length         = atoi( argv[9] );  
    /* threshold number of sigma */
    params->compEPInput->numSigmaMin      = atof( argv[10] ); 
    /* default alpha value for tiles with sigma < numSigmaMin */
    params->compEPInput->alphaDefault     = atof( argv[11] ); 
    /* EK - Max. number of events to communicate to master */
    params->events2Master                 = atoi( argv[14] );
    /* EK - Channel name to process (e.g. "ifodmro") */
    strcpy( params->channelName, argv[15] );
    /* Simulation type:  currently ignored. */
    params->simType                       =  atoi( argv[16] );
    params->simType = 0;
    /* Spectrum method to use */
    if ( !strcmp( argv[17], "useMean" ) ) {
        params->initParams->method            = useMean;
    } 
    else if ( !strcmp( argv[17], "useMedian" ) ) {
        params->initParams->method            = useMedian;
    }
    else {
        ABORT( status, EPSEARCHH_ESPEC, EPSEARCHH_MSGESPEC);
    }
    /* Window to use on the data */
    params->winParams.type                = atoi( argv[18] );

    /* initialize parameter structures */
    params->tfTiling     = NULL;
    params->epSegVec     = NULL;
    params->numSlaves    = NULL;

    /* allocate memory for the conditioned segments */
    LALCreateEPDataSegmentVector (status->statusPtr, &(params->epSegVec), 
            params->initParams);  
    CHECKSTATUSPTR (status);

    /* initialize parameters */
    params->haveData        = 0;
    params->currentSegment  = 0;
    params->numEvents       = 0;
    params->searchMaster    = 0;
    params->tfTilingInput->maxTileBand = 128.0;


    DETATCHSTATUSPTR (status);
    RETURN (status);
}


/********************************************************************
 *
 * This function breaks the time series into segments to be analyzed
 * by the power code
 *
 ********************************************************************/
/* <lalVerbatim file="EPConditionDataCP"> */
void EPConditionData(
        LALStatus             *status,
        REAL4TimeSeries       *series,
	REAL4                  flow,
	REAL8                  resampledeltaT,
	ResampleTSFilter       resampFiltType,
        void                  *searchParams
        )
/* </lalVerbatim> */
{
    UINT4                         i,j;
    INT8                          dataTimeNS  = 0;
    REAL4                        *dummyData    = NULL;  
    REAL4                         fsafety=0;
    EPDataSegmentVector          *dataSegVec = NULL;
    EPSearchParams               *params;
    PassBandParamStruc            highpassParam;

    /* resample parameters */
    UINT4                  resampleChan = 0;
    ResampleTSParams       resampleParams;
    const REAL8            epsilon = 1.0e-8;

    INITSTATUS (status, "LALConditionData", EPSEARCHC);
    ATTATCHSTATUSPTR (status);

    params = (EPSearchParams *) searchParams;

    /****************************************************************
     *
     * identify multiDimData sequences by their names 
     *
     ***************************************************************/

    ASSERT (searchParams, status, EPSEARCHH_ENULLP, EPSEARCHH_MSGENULLP);

    /****************************************************************
     * 
     * check that there is enough data to construct the required segments
     * :TODO:  This assumes that overlap is half the number of points
     * in a segment.  
     *
     * Also,  ovrlap points of data are ignored at the beginning and
     * end because of the Butterworth filter corruption.  As a result,
     * we should have 2 * ovrlap more points than one expects
     *
     * This changes with new time-frequency plane construction
     *
     ****************************************************************/

    if ( series->data->length < 
            params->initParams->numPoints * ( params->initParams->numSegments - 1 ) / 2 +
            4 * params->ovrlap )
    {
        ABORT (status, EPSEARCHH_EDATZ, EPSEARCHH_MSGEDATZ); 
    }


    /* **************************************************************
     *
     * Resample the data if necessary
     *
     * *************************************************************/

    /* set the time series parameters of the input data and resample params */
    memset( &resampleParams, 0, sizeof(ResampleTSParams) );
    resampleParams.deltaT = resampledeltaT;

    if( ! ( fabs( resampleParams.deltaT - series->deltaT ) < epsilon ) )
      {
	resampleChan = 1;
	
	if ( resampFiltType == 0 )
	  {
	    resampleParams.filterType = LDASfirLP;
	  }
	else if ( resampFiltType == 1 )
	  {
	    resampleParams.filterType = defaultButterworth;
	  }
      }
    
    /*resample if required to */
    if ( resampleChan )
      {
	LALResampleREAL4TimeSeries(status->statusPtr, series, &resampleParams);
	CHECKSTATUSPTR (status);
      }

    /* ***************************************************************
     *
     * translate InPut to DataSegment and allocate memory for response
     * function data
     *
     ****************************************************************/
    dataSegVec = params->epSegVec;

    /* Set up for a highpass filter */
    highpassParam.nMax = 4;
    fsafety = params->tfTilingInput->flow - 10.0;
    highpassParam.f2 = fsafety > 150.0 ? 150.0 : fsafety;
    highpassParam.f1 = -1.0;
    highpassParam.a2 = 0.1;
    highpassParam.a1 = -1.0;
    LALButterworthREAL4TimeSeries(status->statusPtr, series, &highpassParam);
    CHECKSTATUSPTR (status);
            
    /* Point dummyData ovrlap into series to avoid corruption */
    dummyData = series->data->data;
    /* no need for the next line with segments twice as long 
     * dummyData += (params->ovrlap);
     */
    dataTimeNS  = 1000000000L * (INT8) series->epoch.gpsSeconds;
    dataTimeNS += (INT8) series->epoch.gpsNanoSeconds;
    /* no need for the next line with segments twice as long 
     * dataTimeNS += (INT8) (1e9 * params->ovrlap * series->deltaT);
     */

    /****************************************************************
     *
     * I think memory allocation could be substantially reduced once the 
     * INT2 -> REAL4 since we dont need to allocate this memory,  rather we
     * just use the wrapperAPI's data.  This should be ok.  I need to check
     * that the wrapperAPI does not destroy the memory until it cleans up at
     * the end of an analysis. -- PRB 08/14/01
     *
     ****************************************************************/

    for ( i = 0; i < dataSegVec->length; ++i )
    {
        /* point to current segment */
        EPDataSegment *dummySegment = dataSegVec->data + i;
        INT4 ptsPerSeg = 2 * params->initParams->numPoints;

        /* 
         * Would only have relevance in a standalone code.  For wrapperAPI,  
         * we know how much data we are getting
         */
        dummySegment->endOfData = 0;
        /* this should be set to a unique number for each segment   */
        dummySegment->number = i;

        /* copy the ifodmro */
        for ( j = 0; j < ptsPerSeg ; ++j)
        {
            dummySegment->data->data->data[j] = (REAL4) dummyData[j];
        }
        dummySegment->data->data->length = ptsPerSeg;
        dummyData += (params->initParams->numPoints - params->ovrlap);


        strncpy( dummySegment->data->name, INPUTNAME_CHANNEL, 
                LALNameLength * sizeof(CHAR) );
        dummySegment->data->deltaT = series->deltaT;

        {
            INT8 dummyNS = 0;

            dummyNS = dataTimeNS + (INT8) (1e9 * 
                    (params->initParams->numPoints - params->ovrlap) 
                    * i * dummySegment->data->deltaT);

            dummySegment->data->epoch.gpsSeconds = 
                (INT4) (dummyNS/1000000000L);

            dummySegment->data->epoch.gpsNanoSeconds = 
                (INT4) (dummyNS%1000000000L);
        }

        dummySegment->data->f0 = 0.0;

      /* copy the spectrum:  this does nothing at the moment */
      dummySegment->spec->data->length = (ptsPerSeg/2 + 1);
      strncpy( dummySegment->spec->name, INPUTNAME_SPECTRUM,
          LALNameLength * sizeof(CHAR) );
      dummySegment->spec->deltaF = 1.0 / 
          ( series->deltaT * (REAL4) ptsPerSeg );

      dummySegment->spec->epoch.gpsSeconds =
          dummySegment->data->epoch.gpsSeconds;
      dummySegment->spec->epoch.gpsNanoSeconds = 
          dummySegment->data->epoch.gpsNanoSeconds;

      dummySegment->spec->f0 = 0.0;
    }

    /****************************************************************
     * 
     * clean up and return
     *
     ****************************************************************/
    DETATCHSTATUSPTR (status);
    RETURN (status);

}


/********************************************************************
 *
 * This function cleans up all the search stuff
 *
 ********************************************************************/

/* <lalVerbatim file="EPFinalizeSearchCP"> */
void
EPFinalizeSearch(
        LALStatus             *status,
        void                 **searchParams
        )
/* </lalVerbatim> */
{
  EPSearchParams *params;

  INITSTATUS (status, "LALFinalizeSearch", EPSEARCHC);
  ATTATCHSTATUSPTR (status);

  /* check searchParams exists and cast approriately */
  if ( !(*searchParams) ){
    ABORT( status, EPSEARCHH_ENULLP, EPSEARCHH_MSGENULLP );
  }
  params = (EPSearchParams *) *searchParams;

  /* destroy memory for conditioned segments [ Now moved inside the power code]*/
  /* LALDestroyEPDataSegmentVector(status->statusPtr, &(params->epSegVec));*/
  
  if ( params->searchMaster )
  {
    /* free numSlaves */
    LALFree (params->numSlaves);
  }
  
  /* free compEPInput */
  LALFree (params->compEPInput);

  /* free EPInitParams */
  LALFree (params->initParams); 

  /* free TFTiling initialisation parameters */
  LALFree (params->tfTilingInput);
  
  /* EK - modified; added channel name */
  LALFree(params->channelName);

  /* free searchParams */
  LALFree (*searchParams);
  *searchParams = NULL;

  DETATCHSTATUSPTR (status);
  RETURN (status);
}



/********************************** <lalVerbatim file="PowerConditionDataCV">
Author: Anderson, W G, Brady, P and Brown, D
$Id$
**************************************************** </lalVerbatim> */

#include <lal/LALStdlib.h>
#include <lal/ExcessPower.h>
#include <lal/TFTransform.h>
#include <lal/EPData.h>
#include "ConditionData.h"

NRCSID( CONDITIONDATAC , "power $Id$" );

void
LALConditionData(
    LALStatus             *status,
    LALSearchInput        *inout,
    void                  *searchParams,
    LALMPIParams          *mpiParams
    )
{
  INT4                          i,j;
  INT4                          dmroIndex  = 0;
  INT4                          specIndex  = 0;
  INT4                          respIndex  = 0;
  INT8                          dataTimeNS  = 0;
  REAL4                         *dummyData    = NULL;  /* Should be REAL4 */
  REAL4                        *dummyR4    = NULL;
  COMPLEX8                     *dummyC8    = NULL;
  EPDataSegmentVector          *dataSegVec = NULL;
  EPSearchParams               *params;

  INITSTATUS (status, "LALConditionData", CONDITIONDATAC);
  ATTATCHSTATUSPTR (status);

  params = (EPSearchParams *) searchParams;

  /*
   *
   * do not condition the data on the search master
   *
   */

  if ( params->searchMaster )
  {
    /* return immediately */
    DETATCHSTATUSPTR (status);
    RETURN (status);
  }
  else
  {

    /*
     *
     * identify multiDimData sequences by their names 
     *
     */

    ASSERT (searchParams, status, CONDITIONDATAH_ENULL, CONDITIONDATAH_MSGENULL);

    for ( i = 0; i < inout->numberSequences; ++i )
    {
        if ( strstr( inout->sequences[i].name, INPUTNAME_CHANNEL))  
        {
        dmroIndex = i;
        dummyData   = inout->sequences[i].data.real4;
        dataTimeNS  = 1000000000L * 
            (INT8) inout->sequences[i].range.dTime.startSec;
        dataTimeNS += (INT8) inout->sequences[i].range.dTime.startNan;
      }
      else if ( strstr( inout->sequences[i].name, INPUTNAME_SPECTRUM ) )
      {
        specIndex = i;
        dummyR4   = inout->sequences[i].data.real4;
      }
      else if ( strstr( inout->sequences[i].name, INPUTNAME_RESPONSE ) )
      {
        respIndex = i;
        dummyC8   = inout->sequences[i].data.complex8;
      }
      else
        ABORT (status, CONDITIONDATAH_EINPUT, CONDITIONDATAH_MSGEINPUT); 
    }


    /*
     * 
     * check that there is enough data to construct the required segments
     *
     */

    if ( inout->sequences[dmroIndex].dimensions[0] < 
        params->initParams->numPoints * ( params->initParams->numSegments + 1 ) / 2 )
    {
      ABORT (status, CONDITIONDATAH_EDATZ, CONDITIONDATAH_MSGEDATZ); 
    }
    

    /* 
     *
     * translate InPut to DataSegment 
     *
     */

    dataSegVec = params->epSegVec;

    /*
     *
     * I think memory allocation could be substantially reduced once the 
     * INT2 -> REAL4 since we dont need to allocate this memory,  rather we
     * just use the wrapperAPI's data.  This should be ok.  I need to check
     * that the wrapperAPI does not destroy the memory until it cleans up at
     * the end of an analysis. -- PRB 08/14/01
     *
     */
    
    for ( i = 0; i < dataSegVec->length; ++i )
    {
      /* point to current segment */
      EPDataSegment *dummySegment = dataSegVec->data + i;

      /* 
       * Would only have relevance in a standalone code.  For wrapperAPI,  
       * we know how much data we are getting
       */
      dummySegment->endOfData = 0;
      /* this should be set to a unique number for each segment   */
      dummySegment->number = i;

      /* copy the ifodmro */
      for ( j = 0; j < params->initParams->numPoints ; ++j)
      {
        dummySegment->data->data->data[j] = (REAL4) dummyData[j];
      }
      dummySegment->data->data->length = params->initParams->numPoints;
#ifdef POWER_SO_DBGLVL3
        {
          FILE *fpout;
          char fname[100];
          INT4 tmpi;

          sprintf(fname,"data%d.dat",i);
          fpout = fopen(fname,"w");
          for (tmpi =0 ; tmpi < 1026 ; tmpi++)
          { 
            fprintf( fpout, "%e \n", dummyData[tmpi]); 
          }
          fclose(fpout);
        }
#endif
      dummyData += (params->initParams->numPoints - params->ovrlap);


      strncpy( dummySegment->data->name, INPUTNAME_CHANNEL, 
          LALNameLength * sizeof(CHAR) );
      dummySegment->data->deltaT = inout->sequences[dmroIndex].range.dTime.timeStepSize;

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

      /* copy the spectrum */
       if ( inout->sequences[specIndex].range.dFreq.startFreq < 0 )
       {
         /* ...if the spectrum has numPoints samples... */
         if ( inout->sequences[specIndex].range.dFreq.numberSamples
             == params->initParams->numPoints )
         {
           /* copy the two sided spectrum */
           for ( j = 0; j < params->initParams->numPoints/2 ; ++j )
           {
            dummySegment->spec->data->data[j] = 
               2.0 * dummyR4[j + params->initParams->numPoints/2];
           }
           dummySegment->spec->data->data[params->initParams->numPoints/2] 
             = 2.0 * dummyR4[0];
         }
         else /* ...otherwise abort */
         {
           ABORT( status, CONDITIONDATAH_ESPEC, CONDITIONDATAH_MSGESPEC );
         }
       }
       else /* it _could_ be a one sided psd... */
       {
         /* ...if the spectrum has numPoints / 2 + 1 samples... */
         if ( inout->sequences[specIndex].range.dFreq.numberSamples
             == (params->initParams->numPoints/2 + 1) )
         {
           memcpy( dummySegment->spec->data->data, dummyR4, 
                 (params->initParams->numPoints/2 +1)*sizeof(REAL4) ); 
         }
         else /* ...otherwise abort */
        {
           ABORT( status, CONDITIONDATAH_ESPEC, CONDITIONDATAH_MSGESPEC );
         }
       }

      dummySegment->spec->data->length = (params->initParams->numPoints/2 + 1);
      strncpy( dummySegment->spec->name, INPUTNAME_SPECTRUM,
          LALNameLength * sizeof(CHAR) );
      dummySegment->spec->deltaF = inout->sequences[specIndex].range.dFreq.freqStepSize;

      dummySegment->spec->epoch.gpsSeconds = 
        ( inout->sequences[specIndex].range.dTime.startSec );

      dummySegment->spec->epoch.gpsNanoSeconds = 
        ( inout->sequences[specIndex].range.dTime.startNan );

      dummySegment->spec->f0 = 0.0;

      /* copy the response */
      memcpy( dummySegment->resp->data->data, dummyC8, 
          (params->initParams->numPoints/2 +1)*sizeof(COMPLEX8) ); 

      dummySegment->resp->data->length=(params->initParams->numPoints/2 + 1);
      strncpy( dummySegment->resp->name, INPUTNAME_RESPONSE, 
          LALNameLength * sizeof(CHAR) );

      dummySegment->resp->deltaF = 
        inout->sequences[respIndex].range.dFreq.freqStepSize;

      dummySegment->resp->epoch.gpsSeconds = 
        ( inout->sequences[respIndex].range.dTime.startSec );

      dummySegment->resp->epoch.gpsNanoSeconds = 
        ( inout->sequences[respIndex].range.dTime.startNan );  
      dummySegment->resp->f0 = 0.0;

    }
  }
  /*
   * 
   * clean up and return
   *
   */

  DETATCHSTATUSPTR (status);
  RETURN (status);
}



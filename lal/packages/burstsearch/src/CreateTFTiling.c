/******** <lalVerbatim file="CreateTFTilingCV"> ********
Author: Eanna Flanagan
$Id$
********* </lalVerbatim> **********/


#include <lal/LALRCSID.h>


NRCSID (CREATETFTILINGC, "$Id$");


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <lal/ExcessPower.h>
#include <lal/LALConstants.h>
#include <lal/LALErrno.h>
#include <lal/LALStdlib.h>
#include <lal/Random.h>
#include <lal/RealFFT.h>
#include <lal/SeqFactories.h>
#include <lal/Thresholds.h>


#define TRUE 1
#define FALSE 0

extern INT4 lalDebugLevel;

static INT4 pow1(INT4 a, INT4 b)
{
  /* returns a^b */
  INT4 t=1;
  INT4 i;
  for(i=0;i<b;i++) t*=a;
  return(t);
}

/******** <lalVerbatim file="CreateTFTilingCP"> ********/
void
LALCreateTFTiling (
		 LALStatus              *status,
		 TFTiling               **tfTiling,
		 CreateTFTilingIn       *input,
                 TFPlaneParams          *planeParams
		 )
/******** </lalVerbatim> *********/
{
  INT4          numPlanes;
  INT4          tileCount=0;

  REAL8         fhigh;        /* max frequency of the TF plane */
  REAL8         flow;         /* min frequency of the TF plane */
  REAL8         timeDuration; /* time duration of the plane    */

  TFTile        **currentTile;

  INITSTATUS (status, "LALCreateTFTiling", CREATETFTILINGC);
  ATTATCHSTATUSPTR (status);

  /* Check input structure: report if NULL */
  ASSERT(input, status, LAL_NULL_ERR, LAL_NULL_MSG);
      
  /* 
   * Check return structure: tfTiling should point to a valid pointer
   * which should not yet point to anything.
   *
   */
  ASSERT(tfTiling != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(*tfTiling == NULL, status, LAL_NNULL_ERR, LAL_NNULL_MSG);

  /*
   *
   *  compute some parameters 
   *
   */

  /* lowest frequency to be used in the plane. */

  flow = planeParams->flow;
  
  /* highest frequency to be used in the plane  */ 

  fhigh = planeParams->fhigh;

  /* time duration to be searched for, this will determine the no. of time
   * bins in the plane. 
   */

  timeDuration = planeParams->timeDuration;  

  /* number of time frequency planes to be constructed 
   *  THERE IS ONLY ONE PLANE
   */
  numPlanes = 1;

  /*  Assign memory for *tfTiling   */
  *tfTiling = LALMalloc(sizeof(**tfTiling));
  ASSERT(*tfTiling, status, LAL_NOMEM_ERR, LAL_NOMEM_MSG);

  /* set some tile parameters */
  (*tfTiling)->numPlanes = numPlanes;
  (*tfTiling)->planesComputed = FALSE;
  (*tfTiling)->excessPowerComputed = FALSE;
  (*tfTiling)->tilesSorted = FALSE;
  (*tfTiling)->tfp = NULL;
  (*tfTiling)->dftParams = NULL;

  /* set things up for recursive generation of linked list below */
  currentTile = &((*tfTiling)->firstTile);
  *currentTile = NULL;

  /* allocate memory for vector of pointers to TF planes 
   * WILL CHANGE TO REAL8TFPLANE : saikat
   */
  (*tfTiling)->tfp = (COMPLEX8TimeFrequencyPlane **) 
       LALMalloc (numPlanes*sizeof(COMPLEX8TimeFrequencyPlane *));

  /*  Make sure that the allocation was succesful */
  if ( !((*tfTiling)->tfp) ){
    LALFree( *tfTiling );
    ABORT (status, LAL_NOMEM_ERR, LAL_NOMEM_MSG);
  }

  /* allocate memory for vector of pointers to DFTParams */
  (*tfTiling)->dftParams = (ComplexDFTParams **) 
       LALMalloc (numPlanes*sizeof(ComplexDFTParams *));

  /*  Make sure that the allocation was succesful */
  if ( !((*tfTiling)->dftParams) ){
    LALFree( (*tfTiling)->tfp );
    LALFree( *tfTiling );
    ABORT (status, LAL_NOMEM_ERR, LAL_NOMEM_MSG);
  }

  /* 
   *  
   *  create the time frequency plane and DFTParams
   *  
   */

    {
      COMPLEX8TimeFrequencyPlane    **thisPlane;
      ComplexDFTParams              **thisDftParams;

      /* setup parameter structure for creating TF plane 
       * 
       */

      planeParams->timeBins = timeDuration / planeParams->deltaT;
      planeParams->freqBins = (fhigh - flow) / planeParams->deltaF;

      /* Create TF plane structure */
      thisPlane = (*tfTiling)->tfp;
      *thisPlane=NULL;
      LALCreateTFPlane( status->statusPtr, thisPlane, planeParams);
      CHECKSTATUSPTR (status);

      /* create the DFTParams structure */
      {
	LALWindowParams winParams;
	winParams.type=Rectangular;
	winParams.length=planeParams->timeBins;
	thisDftParams  = (*tfTiling)->dftParams;
	*thisDftParams = NULL;
	LALCreateComplexDFTParams( status->statusPtr, thisDftParams, 
                             &winParams, -1);
	CHECKSTATUSPTR (status);
	/* Its an inverse transform instead of a forward transform */
      }

    }

  /* 
   *  
   *  compute the linked list of Time Frequency Tiles
   *
   */

    {
      /* coordinates of a given TF tile */
      INT4                          fstart;
      INT4                          deltaf;
      INT4                          tstart;
      INT4                          deltat;
      INT4                          incrementT;
      INT4                          incrementF;
      INT4                          timeBins ;
      INT4                          freqBins ;

      COMPLEX8TimeFrequencyPlane    **thisPlane;

      /* since only one plane now*/      
      thisPlane = (*tfTiling)->tfp ;      

      /* copy the no. of timebins and freqbins into the variables */
      timeBins = (*thisPlane)->params->timeBins;
      freqBins = (*thisPlane)->params->freqBins;

    /* deltat should correspond to no. of time bins for a particular tile;   
       * and NOTE the plane must be constructed with resolutions
       * TWO times finer than the tile resolution
       */
      deltat = 2;
      while (deltat <= timeBins && 
	     deltat*(*thisPlane)->params->deltaT <= input->maxTileDuration)
	{
	  incrementT = deltat/input->overlapFactor;
	  tstart=0;
	  while (tstart <= timeBins - deltat)
	    {
	      /* deltaf is set by the deltat, requiring that the 
               * TF vol is >=1 
               */
	      deltaf = 1/(deltat*(*thisPlane)->params->deltaT*(*thisPlane)->params->deltaF);
	      while (deltaf <= freqBins && 
                  deltaf*(*thisPlane)->params->deltaF <= input->maxTileBand)
		{
		  incrementF = deltaf/input->overlapFactor;
		  fstart=0;
		  while (fstart <= freqBins - deltaf)
		    {
		      /* 
		       * 
		       *  add new Tile to linked list 
		       *
		       */

		      /*  Assign memory for tile */
		      *currentTile = (TFTile *) LALMalloc(sizeof(TFTile));

		      /*  Make sure that the allocation was succesful */
                      if ( ! (*currentTile) ){
                        ABORT(status, LAL_NOMEM_ERR, LAL_NOMEM_MSG);
                      }

		      /* assign the various fields */
		      (*currentTile)->fstart=fstart;
		      (*currentTile)->fend=fstart+deltaf;
		      (*currentTile)->tstart=tstart;
		      (*currentTile)->tend=tstart+deltat;
		      (*currentTile)->whichPlane=0;
		      (*currentTile)->nextTile=NULL;
                      (*currentTile)->deltaT=(*thisPlane)->params->deltaT;
                      (*currentTile)->deltaF=(*thisPlane)->params->deltaF;
		      (*currentTile)->excessPower=0.0;
		      (*currentTile)->alpha=0.0;
		      (*currentTile)->weight=1.0;
		      (*currentTile)->firstCutFlag=FALSE;


		      /* keep track of how many tiles */
		      tileCount++;

		      /* 
		       *  update currentTile to point at location of 
		       *  pointer to next tile
		       *
		       */
		      currentTile = &((*currentTile)->nextTile);

  		      fstart += incrementF;
		    }
		  deltaf = 2*deltaf;
		};  /* while (deltaf .. */
	      tstart += incrementT;
	    }
	  deltat = 2*deltat;
	};  /* while (deltat .. */
    } 

    /* printf("no. of tile in the plane = %d\n",tileCount);*/

  (*tfTiling)->numTiles=tileCount;

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}

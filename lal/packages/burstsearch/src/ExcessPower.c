/*----------------------------------------------------------------------- 
 * 
 * File Name: ExcessPower.c 
 * 
 * Author: Eanna Flanagan
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * ExcessPower
 * 
 * SYNOPSIS 
 * 
 * DESCRIPTION 
 * 
 * 
 * 
 * DIAGNOSTICS 
 *
 * CALLS
 * 
 * NOTES
 * 
 *-----------------------------------------------------------------------
 */


#include "LALRCSID.h"


NRCSID (EXCESSPOWERC, "$Id$");


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "LALStdlib.h"
#include "LALConstants.h"
#include "SeqFactories.h"
#include "RealFFT.h"
#include "Thresholds.h"
#include "ExcessPower.h"
#include "Random.h"


#define TRUE 1
#define FALSE 0


extern INT4 LALDebugLevel;

static INT4 pow1(INT4 a, INT4 b)
{
  /* returns a^b */
  INT4 t=1;
  INT4 i;
  for(i=0;i<b;i++) t*=a;
  return(t);
}


static int TileCompare( TFTile **tiles1, TFTile **tiles2 )
{
  if ( (*tiles1)->alpha > (*tiles2)->alpha )
    return 1;
  if ( (*tiles1)->alpha < (*tiles2)->alpha )
    return -1;
  return 0;
}

static void 
DestroyTFTile (LALStatus *status, TFTile *tfTile)
{
  /*
   *  this function destroys linked list of tiles
   *  The easiest coding method would be to use a recursive function,
   *  but that would use a lot of memory on the stack, and memory
   *  overflows on stack are hard to trap.  So, code it by copying out the
   *  linked list of pointers into an array of pointers.
   */

  INT4 tileCount;
  INT4 i;
  TFTile *thisTile;
  TFTile **tiles;

  INITSTATUS (status, "DestroyTFTile", EXCESSPOWERC);
  ATTATCHSTATUSPTR (status);

  ASSERT(tfTile, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);

  tileCount=0;
  thisTile = tfTile;
  while (thisTile != NULL)
    {
      tileCount++;
      thisTile = thisTile->nextTile;
    }


  /* allocate memory for array of pointers to tiles */
  tiles = (TFTile **) LALMalloc (tileCount * sizeof(TFTile *));
  
  /*  Make sure that the allocation was succesful */
  ASSERT (tiles, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);

  tileCount=0;
  thisTile = tfTile;

  while (thisTile != NULL)
    {
      tileCount++;
      *(tiles + tileCount-1) = thisTile;
      thisTile = thisTile->nextTile;
    }
  
  for(i=tileCount-1; i>=0; i--)
    {
      thisTile = *(tiles + i);
      LALFree (thisTile);
    }

  LALFree (tiles);

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}





static void 
PrintTFTile (
	     LALStatus                                 *status,
	     FILE                                   *fp,
	     TFTile                                 *tfTile,
	     TFTiling                               *tfTiling
	     )
{
  /*
   *  this function prints out information contained in one time-frequency
   *  tile.
   */

  INT4 t1;
  INT4 t2;
  INT4 f1;
  INT4 f2;
  INT4 dof;
  COMPLEX8TimeFrequencyPlane *thisPlane=NULL;
  REAL8 flow;
  REAL8 epoch;
  REAL8 deltaT;
  

  INITSTATUS (status, "PrintTFTile", EXCESSPOWERC);
  ATTATCHSTATUSPTR (status);

  ASSERT(tfTile, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);
  ASSERT(tfTiling, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);
  ASSERT(tfTiling->tfp, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);
  ASSERT(fp, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);

  t1 = tfTile->tstart;
  t2 = tfTile->tend;
  f1 = tfTile->fstart;
  f2 = tfTile->fend;
  dof = 2*(t2-t1+1)*(f2-f1+1);
  thisPlane = *(tfTiling->tfp + tfTile->whichPlane);
  flow = thisPlane->params->flow;
  epoch = (REAL8)(thisPlane->epoch.gpsSeconds) + (REAL8)(thisPlane->epoch.gpsNanoSeconds)/1000000000.0;
  deltaT = thisPlane->params->deltaT;

  fprintf(fp,"\n");
  fprintf(fp," Time frequency tile in TF Plane number %d\n", 
          tfTile->whichPlane);
  fprintf(fp," Frequency interval: %f Hz to %f Hz, i.e.,  bins %d to %d of %d\n",
	 flow + (REAL8)(f1)/deltaT, flow + (REAL8)(f2+1)/deltaT, f1, f2,
	 thisPlane->params->freqBins);
  fprintf(fp," Time interval    :  %f s to %f s, i.e.,  bins %d to %d of %d\n",
	 epoch + (REAL8)(t1)*deltaT, epoch + (REAL8)(t2+1)*deltaT, t1, t2,
	 thisPlane->params->timeBins);
  fprintf(fp," Total number of degrees of freedom:  %d\n",dof);
  fprintf(fp," Excess power:  %f,   1 / alpha    :  %f\n", tfTile->excessPower,
         1/tfTile->alpha);
  {
    /* print out effective number of sigma */
    REAL8 sigma2;
    Chi2ThresholdIn input;

    input.dof = 1.0;
    input.falseAlarm = tfTile->alpha;
    LALChi2Threshold (status->statusPtr, &sigma2, &input);
    CHECKSTATUSPTR (status);
    fprintf(fp," Effective number of sigma:  %f\n", sqrt(sigma2));
  }

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}










void
LALAddWhiteNoise (
	       LALStatus                               *status,
	       COMPLEX8Vector                       *v,
	       REAL8                                noiseLevel
	       )
{
  /*
   *
   *  Add white noise to complex vector
   *
   */

  RandomParams           *params=NULL;
  REAL4Vector            *vr=NULL;
  REAL4Vector            *vi=NULL;
  INT4                   i;

  INITSTATUS (status, "LALAddWhiteNoise", EXCESSPOWERC);
  ATTATCHSTATUSPTR (status);


  /* make sure that arguments are not NULL */
  ASSERT (v, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);
  ASSERT (v->data, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);

  /* make sure length of series is nonzero */
  ASSERT (v->length>0, status, EXCESSPOWER_EINCOMP, 
          EXCESSPOWER_MSGEINCOMP);

  
  /* Seed Random Number Generator with current time for seed */
  LALCreateRandomParams (status->statusPtr, &params, 0);  
  CHECKSTATUSPTR (status);

  /* create temporary vectors */
  LALSCreateVector (status->statusPtr, &vr, v->length);
  CHECKSTATUSPTR (status);
  LALSCreateVector (status->statusPtr, &vi, v->length);
  CHECKSTATUSPTR (status);
  
  /* Fill temporary vectors with Gaussian deviates */
  LALNormalDeviates (status->statusPtr, vr, params);
  CHECKSTATUSPTR (status);
  LALNormalDeviates (status->statusPtr, vi, params);
  CHECKSTATUSPTR (status);

  for(i=0;i<v->length;i++) 
    {
      v->data[i].re += noiseLevel * vr->data[i];
      v->data[i].im += noiseLevel * vi->data[i];
    }
    
  LALSDestroyVector (status->statusPtr, &vr);
  CHECKSTATUSPTR (status);
  
  LALSDestroyVector (status->statusPtr, &vi);
  CHECKSTATUSPTR (status);
  
  LALDestroyRandomParams (status->statusPtr, &params);
  CHECKSTATUSPTR (status);
  
  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}    












void
LALCreateTFTiling (
		 LALStatus                             *status,
		 TFTiling                           **tfTiling,
		 CreateTFTilingIn                   *input
		 )
{
  INT4                          numPlanes;
  INT4                          nf;
  INT4                          i;
  INT4                          tileCount=0;

  REAL8                         fhigh;
  REAL8                         flow;

  TFTile                        **currentTile;

  INITSTATUS (status, "LALCreateTFTiling", EXCESSPOWERC);
  ATTATCHSTATUSPTR (status);

  /* Check input structure: report if NULL */
  ASSERT (input, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);
      

  /* 
   * Check return structure: tfTiling should point to a valid pointer
   * which should not yet point to anything.
   *
   */
  ASSERT (tfTiling != NULL, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);
  ASSERT (*tfTiling == NULL, status, EXCESSPOWER_ENONNULL, 
           EXCESSPOWER_MSGENONNULL);


  /* 
   *
   *  Make sure that input parameters are reasonable, compatible etc.
   *
   */


  ASSERT (input->overlapFactor > 1, status, 
          EXCESSPOWER_EPOSARG, EXCESSPOWER_MSGEPOSARG);
  ASSERT (input->length > 0, status,
          EXCESSPOWER_EPOSARG, EXCESSPOWER_MSGEPOSARG);
  ASSERT (input->deltaF > 0.0, status,
          EXCESSPOWER_EPOSARG, EXCESSPOWER_MSGEPOSARG);
  ASSERT (input->flow >= 0.0, status,
          EXCESSPOWER_EPOSARG, EXCESSPOWER_MSGEPOSARG);
  ASSERT (input->minFreqBins > 0, status, 
          EXCESSPOWER_EPOSARG, EXCESSPOWER_MSGEPOSARG);
  ASSERT (input->minTimeBins > 0, status, 
          EXCESSPOWER_EPOSARG, EXCESSPOWER_MSGEPOSARG);


  /*
   *
   *  compute some parameters 
   *
   */


  /* lowest and highest frequency to be used */
  flow = input->flow;
  fhigh = input->flow + input->deltaF*(REAL8)(input->length);

  /* number of frequency bins to be used in computation */
  nf = input->length;

  /* 
   * minimum sizes of TF tiles to be searched over in time
   * and freq directions must not exceed length of data used
   *
   */
  ASSERT(input->minFreqBins < nf, status, EXCESSPOWER_EINCOMP, 
          EXCESSPOWER_MSGEINCOMP);
  ASSERT(input->minTimeBins < nf, status, EXCESSPOWER_EINCOMP, 
          EXCESSPOWER_MSGEINCOMP);

  /* number of time frequency planes to be constructed */
  numPlanes = 1+(INT4)(0.5+log( (REAL8)(nf)) / log(2.0));
  /*printf("nf: %d,  numPlanes: %d\n", nf, numPlanes);*/

  /* check that length of data to be used is a power of 2 */
  ASSERT( nf == pow1(2, numPlanes-1), status, EXCESSPOWER_EPOW2, 
          EXCESSPOWER_MSGEPOW2);

  /*  Assign memory for *tfTiling   */
  *tfTiling = (TFTiling *) LALMalloc(sizeof(TFTiling));
  
  /*  Make sure that the allocation was succesful */
  ASSERT (*tfTiling, status, EXCESSPOWER_EMALLOC, EXCESSPOWER_MSGEMALLOC);

  /* set some parameters */
  (*tfTiling)->numPlanes = numPlanes;
  (*tfTiling)->planesComputed = FALSE;
  (*tfTiling)->excessPowerComputed = FALSE;
  (*tfTiling)->tilesSorted = FALSE;

  /* set things up for recursive generation of linked list below */
  currentTile = &((*tfTiling)->firstTile);
  *currentTile = NULL;

  /* allocate memory for vector of pointers to TF planes */
  (*tfTiling)->tfp = (COMPLEX8TimeFrequencyPlane **) 
       LALMalloc (numPlanes*sizeof(COMPLEX8TimeFrequencyPlane *));

  /*  Make sure that the allocation was succesful */
  ASSERT ((*tfTiling)->tfp, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);


  /* allocate memory for vector of pointers to DFTParams */
  (*tfTiling)->dftParams = (ComplexDFTParams **) 
       LALMalloc (numPlanes*sizeof(ComplexDFTParams *));

  /*  Make sure that the allocation was succesful */
  ASSERT ((*tfTiling)->dftParams, status, EXCESSPOWER_ENULLP, 
          EXCESSPOWER_MSGENULLP);




  /* 
   *  
   *  create the set of time frequency planes and DFTParams
   *
   */


  for(i=0;i<numPlanes;i++)
    {
      TFPlaneParams                 params;
      COMPLEX8TimeFrequencyPlane    **thisPlane;
      ComplexDFTParams              **thisDftParams;

      /* setup parameter structure for creating TF plane */
      params.timeBins = pow1(2,i);
      params.freqBins = nf / params.timeBins;

      params.deltaT = 1.0 / ( (REAL8)(params.timeBins) * input->deltaF);
      params.flow = flow;

      /* Create TF plane structure */
      thisPlane = (*tfTiling)->tfp + i;
      *thisPlane=NULL;
      LALCreateTFPlane( status->statusPtr, thisPlane, &params);
      CHECKSTATUSPTR (status);

      /* create the DFTParams structure */
      {
	LALWindowParams winParams;
	winParams.type=Rectangular;
	winParams.length=params.timeBins;

	thisDftParams  = (*tfTiling)->dftParams + i;
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


  /* loop over time-frequency Planes */
  for(i=0; i<numPlanes; i++)
    {
      /* coordinates of a given TF tile */
      INT4                          fstart;
      INT4                          deltaf;
      INT4                          tstart;
      INT4                          deltat;
      INT4                          incrementT;
      INT4                          incrementF;
  
      INT4                          timeBins = pow1(2,i);
      INT4                          freqBins = nf/timeBins; 

      

      deltat=input->minTimeBins;
      while (deltat <= timeBins)      
	{
	  incrementT = 1+deltat/input->overlapFactor;
	  tstart=0;
	  while (tstart <= timeBins - deltat)
	    {
	      deltaf=input->minFreqBins;
	      while (deltaf <= freqBins)      
		{
		  incrementF = 1+deltaf/input->overlapFactor;
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
		      ASSERT (*currentTile, status, EXCESSPOWER_EMALLOC, 
			      EXCESSPOWER_MSGEMALLOC); 

		      /* assign the various fields */
		      (*currentTile)->fstart=fstart;
		      (*currentTile)->fend=fstart+deltaf-1;
		      (*currentTile)->tstart=tstart;
		      (*currentTile)->tend=tstart+deltat-1;
		      (*currentTile)->whichPlane=i;
		      (*currentTile)->nextTile=NULL;
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
		  
		  deltaf += incrementF;
		};  /* while (deltaf .. */
	      tstart += incrementT;
	    }
	  deltat += incrementT;
	};  /* while (deltat .. */
    } /* for(i=.. */


  (*tfTiling)->numTiles=tileCount;

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}






void
LALDestroyTFTiling (
		 LALStatus                             *status,
		 TFTiling                           **tfTiling
		 )
{
  INT4                         i;

  INITSTATUS (status, "LALDestroyTFTiling", EXCESSPOWERC);
  ATTATCHSTATUSPTR (status);

  /* make sure that arguments are not null */
  ASSERT (tfTiling, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);
  ASSERT (*tfTiling, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);
  ASSERT ((*tfTiling)->tfp, status, EXCESSPOWER_ENULLP, 
	  EXCESSPOWER_MSGENULLP); 
  ASSERT ((*tfTiling)->firstTile, status, EXCESSPOWER_ENULLP, 
      EXCESSPOWER_MSGENULLP);

  /* make sure that number of TF planes is positive */
  ASSERT ( (*tfTiling)->numPlanes>0, status, EXCESSPOWER_EPOSARG,
           EXCESSPOWER_MSGEPOSARG);



  /* destroy the set of time frequency planes and DFTParams */
  for(i=0;i<(*tfTiling)->numPlanes;i++)
    {
      COMPLEX8TimeFrequencyPlane    **thisPlane;
      ComplexDFTParams              **thisDftParams;

      thisPlane = (*tfTiling)->tfp + i;
      LALDestroyTFPlane( status->statusPtr, thisPlane );
      CHECKSTATUSPTR (status);

      thisDftParams = (*tfTiling)->dftParams + i;
      LALDestroyComplexDFTParams( status->statusPtr, thisDftParams );
      CHECKSTATUSPTR (status);
    }

  /* free the vector of pointers to TF planes */
  LALFree ( (*tfTiling)->tfp );

  /* free the vector of pointers to DFTParams */
  LALFree ( (*tfTiling)->dftParams );

  /* destroy the linked list of TF Tiles */
  DestroyTFTile (status->statusPtr, (*tfTiling)->firstTile);
  CHECKSTATUSPTR (status);

  /* free tfTiling struct itself */
  LALFree ( *tfTiling );           

  *tfTiling = NULL;	    /* make sure we don't point to freed struct */ 

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}






void
LALComputeTFPlanes (
		 LALStatus                             *status,
		 TFTiling                           *tfTiling,
		 COMPLEX8FrequencySeries            *freqSeries
		 )
{
  INT4               i;

  INITSTATUS (status, "LALComputeTFPlanes", EXCESSPOWERC);
  ATTATCHSTATUSPTR (status);


    /* make sure that arguments are not NULL */
  ASSERT (freqSeries, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);
  ASSERT (freqSeries->data, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);
  ASSERT (freqSeries->data->data, status, EXCESSPOWER_ENULLP,
          EXCESSPOWER_MSGENULLP);

  ASSERT (tfTiling, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);
  ASSERT (tfTiling->firstTile, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);
  ASSERT (tfTiling->numPlanes>0, status, EXCESSPOWER_EPOSARG,
          EXCESSPOWER_MSGEPOSARG);

  /* compute the ith TFPlane from input frequency series */
  for(i=0; i<tfTiling->numPlanes; i++)
    {
      COMPLEX8TimeFrequencyPlane    **thisPlane;
      ComplexDFTParams              **thisDftParams;
      HorizontalTFTransformIn       transformparams;

      thisPlane = tfTiling->tfp + i;
      ASSERT(thisPlane, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);
      ASSERT(*thisPlane, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);

      thisDftParams = tfTiling->dftParams + i;
      ASSERT(thisDftParams, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);
      ASSERT(*thisDftParams, status, EXCESSPOWER_ENULLP, 
              EXCESSPOWER_MSGENULLP);


      /* setup input structure for computing TF transform */
      transformparams.startT=0;  /* not used for horizontal transforms */
      transformparams.dftParams=*thisDftParams;

      /* Compute TF transform */
      LALFreqSeriesToTFPlane( status->statusPtr, *thisPlane, freqSeries, 
                           &transformparams); 
      CHECKSTATUSPTR (status);
    }

  /* set flags saying TF planes have been computed, but not EP or sorted */
  tfTiling->planesComputed=TRUE;
  tfTiling->excessPowerComputed=FALSE;
  tfTiling->tilesSorted=FALSE;

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}








void
LALComputeExcessPower (
		   LALStatus                             *status,
		   TFTiling                           *tfTiling,
		   ComputeExcessPowerIn               *input
		   )

{
  INT4               i;
  TFTile             *thisTile;


  INITSTATUS (status, "LALComputeExcessPower", EXCESSPOWERC);
  ATTATCHSTATUSPTR (status);


    /* make sure that arguments are not NULL */
  ASSERT (tfTiling, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);
  ASSERT (tfTiling->tfp, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);
  ASSERT (tfTiling->dftParams, status, EXCESSPOWER_ENULLP, 
          EXCESSPOWER_MSGENULLP);
  ASSERT (tfTiling->firstTile, status, EXCESSPOWER_ENULLP, 
          EXCESSPOWER_MSGENULLP);
  ASSERT (input, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);
  ASSERT (input->numSigmaMin >= 1.0, status, EXCESSPOWER_EINCOMP,
          EXCESSPOWER_MSGEINCOMP);
  ASSERT ((input->alphaDefault >0.0) && (input->alphaDefault < 1.0), 
          status, EXCESSPOWER_EINCOMP, EXCESSPOWER_MSGEINCOMP);
  

  ASSERT (tfTiling->numPlanes>0, status, EXCESSPOWER_EPOSARG,
          EXCESSPOWER_MSGEPOSARG);

  /* make sure TF planes have already been computed */
  ASSERT (tfTiling->planesComputed, status, EXCESSPOWER_EORDER,
          EXCESSPOWER_MSGEORDER);

  for(i=0; i<tfTiling->numPlanes; i++)
    {
      COMPLEX8TimeFrequencyPlane **thisPlane = tfTiling->tfp + i;
      ASSERT(thisPlane, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);
      ASSERT(*thisPlane, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);
      ASSERT((*thisPlane)->data, status, EXCESSPOWER_ENULLP, 
             EXCESSPOWER_MSGENULLP);
      ASSERT((*thisPlane)->params, status, EXCESSPOWER_ENULLP, 
             EXCESSPOWER_MSGENULLP);
    }
      

  thisTile = tfTiling->firstTile;
  
  while (thisTile != NULL)
    {
      COMPLEX8TimeFrequencyPlane *tfPlane;
      REAL8 sum;
      REAL8 dof;
      REAL8 rho2;
      REAL8 numsigma;
      

      INT4 j;
      INT4 i;
      INT4 nf;
      INT4 nt;
      INT4 t1;
      INT4 t2;
      INT4 f1;
      INT4 f2;


      /* check plane index is in required range */
      ASSERT( (thisTile->whichPlane >=0) && (thisTile->whichPlane <tfTiling->numPlanes), status, EXCESSPOWER_EINCOMP, EXCESSPOWER_MSGEINCOMP); 

      tfPlane = *(tfTiling->tfp + thisTile->whichPlane);
      nf = tfPlane->params->freqBins;   
      ASSERT( nf>0, status, EXCESSPOWER_EPOSARG, EXCESSPOWER_MSGEPOSARG);
      nt = tfPlane->params->timeBins;   
      ASSERT( nt>0, status, EXCESSPOWER_EPOSARG, EXCESSPOWER_MSGEPOSARG);

      t1=thisTile->tstart;
      t2=thisTile->tend;
      f1=thisTile->fstart;
      f2=thisTile->fend;

      ASSERT( (t1>=0) && (t1<=t2) && (t2<nt), status, EXCESSPOWER_EINCOMP, 
              EXCESSPOWER_MSGEINCOMP);
      ASSERT( (f1>=0) && (f1<=f2) && (f2<nf), status, EXCESSPOWER_EINCOMP, 
              EXCESSPOWER_MSGEINCOMP);

      dof = (REAL8)(2*(t2-t1+1)*(f2-f1+1));


      sum=0.0;
      for(j=t1; j<=t2; j++)
	{
	  for(i=f1; i<=f2; i++)
	    {
	      INT4 offset = j*nf;
	      COMPLEX8 z;
	      z = tfPlane->data[offset+i];
	      sum += z.re*z.re + z.im*z.im;
	    }
	}

      rho2 = sum - dof;
      thisTile->excessPower = rho2;
      numsigma = rho2 / sqrt(2*dof);
      thisTile->weight = 1.0;


      /*
       *  need to compute an accurate value of likelihood only if
       *  excess power is greater than a few sigma
       *
       */

      thisTile->alpha =  input->alphaDefault;         /* default value */

      if(numsigma > input->numSigmaMin)
	{
	  ChisqCdfIn input;
	  REAL8 alpha; /* false alarm probability */

	  thisTile->firstCutFlag=TRUE;

	  /* compute alpha value */
	  input.chi2= sum;
	  input.dof = dof;
	  /* input->nonCentral not used by LALChisqCdf() */

	  LALOneMinusChisqCdf( status->statusPtr, &alpha, &input);

	  /* 
           *  trap error where alpha=0.0.
           *  If alpha=0 we replace alpha with a small number,
           *  otherwise code will frequently crash while testing etc.
           *
           */

	  if( ((alpha==0.0) || (1.0/alpha > LAL_REAL8_MAX)) 
               && (status->statusPtr->statusCode==THRESHOLDS_ERANGE) )
	    {
	      status->statusPtr->statusCode=0;
	      alpha = exp(-700.0);
	    }

	  /* check for other possible errors from LALOneMinusChisqCdf() */
	  CHECKSTATUSPTR (status);

	  thisTile->alpha = alpha;
	}
      
      /* go onto next tile */
      thisTile = thisTile->nextTile;
    }
  

  /* set flag saying alpha for each tile has been computed */
  tfTiling->excessPowerComputed=TRUE;

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}












void
LALSortTFTiling (
	      LALStatus                                *status,
	      TFTiling                              *tfTiling
	      )
{
  INT4               tileCount;
  INT4               numTiles;
  TFTile             *thisTile;
  TFTile             **tiles;

  INITSTATUS (status, "LALSortTFTiling", EXCESSPOWERC);
  ATTATCHSTATUSPTR (status);


    /* make sure that arguments are not NULL */
  ASSERT (tfTiling, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);
  ASSERT (tfTiling->tfp, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);
  ASSERT (tfTiling->dftParams, status, EXCESSPOWER_ENULLP, 
          EXCESSPOWER_MSGENULLP);
  ASSERT (tfTiling->firstTile, status, EXCESSPOWER_ENULLP, 
          EXCESSPOWER_MSGENULLP);
  ASSERT (tfTiling->numPlanes>0, status, EXCESSPOWER_EPOSARG,
          EXCESSPOWER_MSGEPOSARG);


  /* make sure excess power has already been computed */
  ASSERT (tfTiling->excessPowerComputed, status, EXCESSPOWER_EORDER,
          EXCESSPOWER_MSGEORDER);

  /* compute number of tiles */
  thisTile = tfTiling->firstTile;
  tileCount=0;
  while (thisTile != NULL)
    {
      tileCount++;
      thisTile = thisTile->nextTile;
    }
  numTiles = tileCount;

  /* 
   *
   *  Make an array of pointers to be used to sort the tiles.
   *
   */

  /* allocate memory for array of pointers to tiles */
  tiles = (TFTile **) LALMalloc (numTiles * sizeof(TFTile *));
  
  /*  Make sure that the allocation was succesful */
  ASSERT (tiles, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);


  /* copy out pointers into array */
  tileCount=0;
  thisTile = tfTiling->firstTile;
  while (thisTile != NULL)
    {
      tileCount++;
      *(tiles + tileCount-1) = thisTile;
      thisTile = thisTile->nextTile;
    }
  
  qsort( tiles, numTiles, sizeof( TFTile * ), TileCompare );

  /* copy sorted array back into linked list */
  { 
    TFTile **currentTile = &(tfTiling->firstTile);

    tileCount=0;
    while (tileCount < numTiles)
      {
	*currentTile = *(tiles + tileCount);
	tileCount++;
	currentTile = &((*currentTile)->nextTile);
      }

    /* correctly terminate the linked list */
    *currentTile = NULL;
  }

  LALFree (tiles);
  
  /* set flag saying tiles have been sorted */
  tfTiling->tilesSorted=TRUE;

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}





void
LALCountEPEvents (
               LALStatus                               *status,
               INT4                                 *numEvents,
               TFTiling                             *tfTiling,
               REAL8                                alphaThreshold
               )
{
  INT4               tileCount;
  TFTile             *thisTile;

  INITSTATUS (status, "LALCountEPEvents", EXCESSPOWERC);
  ATTATCHSTATUSPTR (status);


  /* make sure that arguments are not NULL */
  ASSERT (tfTiling, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);
  ASSERT (tfTiling->firstTile, status, EXCESSPOWER_ENULLP, 
          EXCESSPOWER_MSGENULLP);
  ASSERT (numEvents, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);
  ASSERT (alphaThreshold > 0.0, status, EXCESSPOWER_EPOSARG,
          EXCESSPOWER_MSGEPOSARG);

  /* 
   *  Should call this routine after calling SortTFTiling, so tiles
   *  are already in order 
   *
   */

  /* check that already sorted */
  ASSERT( tfTiling->tilesSorted, status, EXCESSPOWER_EORDER,
          EXCESSPOWER_MSGEORDER);

  thisTile = tfTiling->firstTile;
  tileCount=0;

  while ( (thisTile != NULL) && (thisTile->alpha <= alphaThreshold))
    {
      tileCount++;
      thisTile = thisTile->nextTile;
    }
  
  /* return number of events */
  *numEvents = tileCount;

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}








void
LALComputeLikelihood (
		   LALStatus                             *status,
		   REAL8                              *lambda,
		   TFTiling                           *tfTiling
		   )
{
  REAL8              avglambda=0.0;
  TFTile             *thisTile;


  INITSTATUS (status, "LALComputeLikelihood", EXCESSPOWERC);
  ATTATCHSTATUSPTR (status);


    /* make sure that arguments are not NULL */
  ASSERT (tfTiling, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);
  ASSERT (tfTiling->firstTile, status, EXCESSPOWER_ENULLP, 
          EXCESSPOWER_MSGENULLP);
  ASSERT (lambda, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);

  /* make sure LALComputeExcessPower() has been already called */
  ASSERT ( tfTiling->excessPowerComputed, status, EXCESSPOWER_EORDER, 
           EXCESSPOWER_MSGEORDER);

  thisTile = tfTiling->firstTile;
  while (thisTile != NULL)
    {
      if(thisTile->firstCutFlag)
	{
	  INT4 t1;
	  INT4 t2;
	  INT4 f1;
	  INT4 f2;
	  REAL8 dof;
	  REAL8 lambda;
	  REAL8 rho4;

	  t1=thisTile->tstart;
	  t2=thisTile->tend;
	  f1=thisTile->fstart;
	  f2=thisTile->fend;
	  dof = (REAL8)(2*(t2-t1+1)*(f2-f1+1));
	  rho4 = thisTile->excessPower * thisTile->excessPower;

	  lambda = dof / (rho4 * thisTile->alpha);

	  avglambda += lambda * thisTile->weight;
	}

      /* go onto next tile */
      thisTile = thisTile->nextTile;
    }
  
  /* compute the likelihood averaged over TF tiles */
  avglambda /= (REAL8)(tfTiling->numTiles);

  /* return value of statistic */
  *lambda=avglambda;
  

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}









void 
LALPrintTFTileList (
		 LALStatus                                 *status,
		 FILE                                   *fp,
		 TFTiling                               *tfTiling,
		 INT4                                   maxTiles
		 )
{
  TFTile *thisTile;
  INT4   tileCount=0;

  INITSTATUS (status, "LALPrintTFTileList", EXCESSPOWERC);
  ATTATCHSTATUSPTR (status);

  ASSERT(fp, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP);
  ASSERT(tfTiling, status, EXCESSPOWER_ENULLP, EXCESSPOWER_MSGENULLP); 
  ASSERT(tfTiling->firstTile, status, EXCESSPOWER_ENULLP, \
         EXCESSPOWER_MSGENULLP); 

  ASSERT(tfTiling->excessPowerComputed, status, EXCESSPOWER_EORDER,
         EXCESSPOWER_MSGEORDER);


  thisTile = tfTiling->firstTile;
  while ( (thisTile != NULL) && (tileCount < maxTiles))
    {
      PrintTFTile (status->statusPtr, fp, thisTile, tfTiling);
      tileCount++;
      thisTile = thisTile->nextTile;
    }

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}





















/******** <lalVerbatim file="ComputeExcessPowerCV"> ********
Author: Flanagan, E
$Id$
********* </lalVerbatim> ********/

#include <lal/LALRCSID.h>

NRCSID (COMPUTEEXCESSPOWERC, "$Id$");


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
#include <lal/Random.h>


#define TRUE 1
#define FALSE 0


extern INT4 lalDebugLevel;

/******** <lalVerbatim file="ComputeExcessPowerCP"> ********/
void
LALComputeExcessPower (
		   LALStatus               *status,
		   TFTiling                *tfTiling,
		   ComputeExcessPowerIn    *input
		   )
/******** </lalVerbatim> ********/
{
  INT4               i;
  TFTile             *thisTile;


  INITSTATUS (status, "LALComputeExcessPower", COMPUTEEXCESSPOWERC);
  ATTATCHSTATUSPTR (status);


    /* make sure that arguments are not NULL */
  ASSERT (tfTiling, status, EXCESSPOWERH_ENULLP, EXCESSPOWERH_MSGENULLP);
  ASSERT (tfTiling->tfp, status, EXCESSPOWERH_ENULLP, EXCESSPOWERH_MSGENULLP);
  ASSERT (tfTiling->dftParams, status, EXCESSPOWERH_ENULLP, 
          EXCESSPOWERH_MSGENULLP);
  ASSERT (tfTiling->firstTile, status, EXCESSPOWERH_ENULLP, 
          EXCESSPOWERH_MSGENULLP);
  ASSERT (input, status, EXCESSPOWERH_ENULLP, EXCESSPOWERH_MSGENULLP);
  ASSERT (input->numSigmaMin >= 1.0, status, EXCESSPOWERH_EINCOMP,
          EXCESSPOWERH_MSGEINCOMP);
  ASSERT ((input->alphaDefault >0.0) && (input->alphaDefault < 1.0), 
          status, EXCESSPOWERH_EINCOMP, EXCESSPOWERH_MSGEINCOMP);
  

  ASSERT (tfTiling->numPlanes>0, status, EXCESSPOWERH_EPOSARG,
          EXCESSPOWERH_MSGEPOSARG);

  /* make sure TF planes have already been computed */
  ASSERT (tfTiling->planesComputed, status, EXCESSPOWERH_EORDER,
          EXCESSPOWERH_MSGEORDER);

  for(i=0; i<tfTiling->numPlanes; i++)
    {
      COMPLEX8TimeFrequencyPlane **thisPlane = tfTiling->tfp + i;
      ASSERT(thisPlane, status, EXCESSPOWERH_ENULLP, EXCESSPOWERH_MSGENULLP);
      ASSERT(*thisPlane, status, EXCESSPOWERH_ENULLP, EXCESSPOWERH_MSGENULLP);
      ASSERT((*thisPlane)->data, status, EXCESSPOWERH_ENULLP, 
             EXCESSPOWERH_MSGENULLP);
      ASSERT((*thisPlane)->params, status, EXCESSPOWERH_ENULLP, 
             EXCESSPOWERH_MSGENULLP);
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
      INT4 ii;
      INT4 nf;
      INT4 nt;
      INT4 t1;
      INT4 t2;
      INT4 f1;
      INT4 f2;


      /* check plane index is in required range */
      ASSERT( (thisTile->whichPlane >=0) && (thisTile->whichPlane <tfTiling->numPlanes), status, EXCESSPOWERH_EINCOMP, EXCESSPOWERH_MSGEINCOMP); 

      tfPlane = *(tfTiling->tfp + thisTile->whichPlane);
      nf = tfPlane->params->freqBins;   
      ASSERT( nf>0, status, EXCESSPOWERH_EPOSARG, EXCESSPOWERH_MSGEPOSARG);
      nt = tfPlane->params->timeBins;   
      ASSERT( nt>0, status, EXCESSPOWERH_EPOSARG, EXCESSPOWERH_MSGEPOSARG);

      t1=thisTile->tstart;
      t2=thisTile->tend;
      f1=thisTile->fstart;
      f2=thisTile->fend;

      ASSERT( (t1>=0) && (t1<=t2) && (t2<nt), status, EXCESSPOWERH_EINCOMP, 
              EXCESSPOWERH_MSGEINCOMP);
      ASSERT( (f1>=0) && (f1<=f2) && (f2<nf), status, EXCESSPOWERH_EINCOMP, 
              EXCESSPOWERH_MSGEINCOMP);

      dof = (REAL8)(2*(t2-t1+1)*(f2-f1+1));


      sum=0.0;
      for(j=t1; j<=t2; j++)
	{
	  for(ii=f1; ii<=f2; ii++)
	    {
	      INT4 offset = j*nf;
	      COMPLEX8 z;
	      z = tfPlane->data[offset+ii];
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
	  ChisqCdfIn locinput;
	  REAL8 alpha; /* false alarm probability */

	  thisTile->firstCutFlag=TRUE;

	  /* compute alpha value */
	  locinput.chi2= sum;
	  locinput.dof = dof;
	  /* locinput->nonCentral not used by LALChisqCdf() */

	  LALOneMinusChisqCdf( status->statusPtr, &alpha, &locinput);

	  /* 
           *  trap error where alpha=0.0.
           *  If alpha=0 we replace alpha with a small number,
           *  otherwise code will frequently crash while testing etc.
           *
           */

	  if( ((alpha==0.0) || (1.0/alpha > LAL_REAL8_MAX)) 
               && (status->statusPtr->statusCode==THRESHOLDSH_ERANGE) )
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



/* -----------------------------------
 *     File Name: FindChirpTDtemplate.c
 *     Author: S. Babak
 *------------------------------------
 */

#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/AVFactories.h>
#include <lal/FindChirpTDTemplate.h>

NRCSID (FINDCHIRPTDTEMPLATEC, "$Id$");




void LALFindChirpTDTemplate(
   LALStatus			*status,
   FindChirpTemplate		*fcTmplt,
   InspiralTemplate		*tmplt,
   FindChirpDataParams		*params  /* note that params are different from SP*/
){

  UINT4		n,i;
  UINT4		numPoints = 0;
  REAL4Vector   *waveform;
/*  REAL4		fMax; */
  UINT4 	cutLow;
  UINT4		cutHigh;
  REAL4		norm;
  REAL4 	deltaF;
  
  INITSTATUS( status, "LALFindChirpTDTemplateInit", FINDCHIRPTDTEMPLATEC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */
   
 /* check that the output structures exist */
  ASSERT( fcTmplt, status, FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
  ASSERT( fcTmplt->data, status, FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
  ASSERT( fcTmplt->data->data, status, FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
 
  /* make sure that the parameter structure exists */
  ASSERT( params, status, FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
  ASSERT( tmplt, status, FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL);
 
  /* all the checks on params are done in FindChirpTDTData */
 
  /* check that we have a valid approximant */ 
  ASSERT( params->approximant == TaylorT1 || params->approximant == TaylorT3
        || params->approximant == PadeT1 || params->approximant == EOB,
      status, FINDCHIRPTDT_ETAPX, FINDCHIRPTDT_MSGETAPX );  


    /* check that the fft plans exist */
  ASSERT( params->fwdPlan, status, FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );

  
  /* check that the parameter values are reasonable */
  ASSERT( params->deltaT > 0, status,
      FINDCHIRPSPH_EDELT, FINDCHIRPSPH_MSGEDELT );
  ASSERT( params->fLow >= 0, status,
      FINDCHIRPSPH_EFLOW, FINDCHIRPSPH_MSGEFLOW );

/*  I SHOULD CHECK THAT PN ORDER IS 2 OR 3 !!!!*/


  
  /* compute the Waveform in time domain */

  numPoints = fcTmplt->data->length; /* !!! this must be equal to data
  					segment size (time series). Hope it is!!! */
					
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     Here I assume that all parameters for 
     in InspiralTemplate are set up properly
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

   LALInspiralParameterCalc(status->statusPtr, tmplt); /* to compute tau's
     							   otherwise useless*/
   CHECKSTATUSPTR(status);
   n = 0;

  /* compute the approximate duration of template */

   LALInspiralWaveLength(status->statusPtr, &n, *tmplt);
 
   /* length of template must be at leat 4 times less then data segment length */ 

   ASSERT( n <= numPoints/4, status, FINDCHIRPTDTEMPLATE_TLEN, \
   				FINDCHIRPTDTEMPLATE_MSGETLEN );
   CHECKSTATUSPTR(status);
   
   /* Allocate memory for waveform and create it */

   LALCreateVector(status->statusPtr, &waveform, numPoints);
   CHECKSTATUSPTR(status);

   memset(waveform, 0.0, numPoints*sizeof(REAL4));
   
   LALInspiralWave(status->statusPtr, waveform, tmplt);
   CHECKSTATUSPTR(status);

   /* perform fft of a waveform */

   LALForwardRealFFT( status->statusPtr, fcTmplt->data, waveform,
   params->fwdPlan);
   CHECKSTATUSPTR(status);
   
   /*   ***************************************************
        compute template normalisation constant (segNorm !) 
        ***************************************************            */ 
  
   deltaF = 1.0/((REAL4)numPoints * params->deltaT);

   /* check that fLow is consistent in different structures */

   ASSERT(params->fLow == tmplt->fLower, status, FINDCHIRPTDTEMPLATE_EFLOW,
         FINDCHIRPTDTEMPLATE_MSGEFLOW);

   cutLow = params->fLow/deltaF > 1 ? params->fLow/deltaF : 1;
   
   /*  cutHigh is defined by the last frequency */

   cutHigh = (tmplt->fFinal/deltaF < fcTmplt->data->length-1 ? tmplt->fFinal/deltaF : fcTmplt->data->length-1);


   for (i=cutLow; i<cutHigh; ++i){
      REAL4 re = fcTmplt->data->data[i].re;
      REAL4 im = fcTmplt->data->data[i].im;
      norm += (re*re + im*im)*params->wtildeVec->data[i].re;
   }
   
/* store the normalization (segNorm) as tmpltNorm !!!! */

   fcTmplt->tmpltNorm = norm;

  /* fcTmplt->tmpltNorm = sqrt(2.0*norm*deltaT/(REAL4)numPoints);*/

   
   LALDestroyVector(status->statusPtr, &waveform);
   CHECKSTATUSPTR(status);

    /* normal exit */
   DETATCHSTATUSPTR( status );
   RETURN( status );
   
}



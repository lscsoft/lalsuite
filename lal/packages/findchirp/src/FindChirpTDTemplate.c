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
#include <lal/LALInstpiral.h>
#include <lal/AVFactories.h>
#include <lal/FindChirpTDTemplate.h>

NRCSID (FINDCHIRPTDTEMPLATEC, "$Id$");




void LALFindChirpTDTemplate(
   LALStatus			*status,
   FindChirpTemplate		*fcTmplt,
   InspiralTemplate		*tmplt,
   FindChirpSPDataParams	*params  /* note that params are different from SP*/
){

  UINT4		n,i;
  UINT4		numPoints = 0;
  REAL4Vector   *waveform;
  REAL4		fMax;
  UINT4 	cutLow;
  UINT4		cutHigh;
  REAL4		norm;
  
  
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
  
  /* check that we have a valid approximant */ 
  if ( params->approximant != TaylorT1 && params->approximant != TaylorT2 
  	&& params->approximant != TaylorT3 && params->approximant != PadeT1 
	&& params->approximant != EOB ) {
    ABORT( status, FINDCHIRPSPH_EUAPX, FINDCHIRPSPH_MSGEUAPX );
  }

    /* check that the fft plans exist */
  ASSERT( params->fwdPlan, status, FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );

  
  /* check that the parameter values are reasonable */
  ASSERT( params->deltaT > 0, status,
      FINDCHIRPSPH_EDELT, FINDCHIRPSPH_MSGEDELT );
  ASSERT( params->fLow >= 0, status,
      FINDCHIRPSPH_EFLOW, FINDCHIRPSPH_MSGEFLOW );

  
  /* compute the Waveform in time domain */

  numPoints = 2*(fcTmplt->data->length-1); /* !!! this must be equal to data
  					segment size (time series) !!! */
					
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     Here I assume that all parameters for 
     in InspiralTemplate are set up properly
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

   LALInspiralParameterCalc(status->statusPtr, tmplt); /* to compute tau's
     							   otherwise useless*/
   CHECKSTATUSPTR(status);
   n = 0;
   LALInspiralWaveLength(status-statusPtr, &n, tmplt); 
   /* length of template must be at leat 4 times less then data segment length */ 
   ASSERT( n <= numPoints/4, status, FINDCHIRPTDTEMPLATE_TLEN, \
   				FINDCHIRPTDTEMPLATE_MSGETLEN );
   CHECKSTATUSPTR(status);
   
   /* Allocate memory for waveform and create it */
   LALCreateVector(status->statusPtr, waveform, numPoints);
   CHECKSTATUSPTR(status);
   
   memset(waveform, 0, waveform);
   
   LALInspiralWave(status->statusPtr, waveform, tmplt);
   CHECKSTATUSPTR(status);

   /* perform fft of a waveform */

   LALForwardRealFFT( status->statusPtr, fcTmplt->data->data, waveform,
   params->fwdPlan);
   CHECKSTATUSPTR(status);
   
   /* compute template normalisation constant */ 
  
   REAL4 deltaF = 1.0/((REAL4)numPoints * params->deltaT);
   cutLow = params->fLow/deltaF > 1 ? params->fLow/deltaF : 1;
   /*  cutHigh is put so far to the Nyquist frequency */
   cutHigh = fcTmplt->data->length-1;
   for (i=cutLow; i<cutHigh; ++i){
      REAL4 re = fcTmplt->data->data[i].re;
      REAL4 im = fcTmplt->data->data[i].im;
      norm += (re*re + im*im)*wtilde[i].re;
   }
    
   fcTmplt->tmpltNorm = sqrt(2.0*norm*deltaT/(REAL4)numPoints);

   
   LALDestroyVector(status->statusPtr, waveform);
   CHECKSTATUSPTR(status);

    /* normal exit */
   DETATCHSTATUSPTR( status );
   RETURN( status );
   
}



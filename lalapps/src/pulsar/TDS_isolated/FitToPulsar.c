/*
*  Copyright (C) 2007  Rejean Dupuis
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

/**
 * \author Dupuis, R. J.
 *
 * \heading{Module \ref FitToPulsar.c}
 *
 * \heading{Algorithm}
 * To be completed.
 *
 * \heading{Uses}
 * \code
 * LALSCreateVector()
 * LALSDestroyVector()
 * LALComputeDetAMResponse()
 * \endcode
 */

#include <math.h>
#include <lal/LALStdlib.h>
#include "FitToPulsar.h" /* use local version until up-to-date file commited to lal */

#define INICHISQU 1.e10



void
LALCoarseFitToPulsar	( 	LALStatus            *status,
		    		CoarseFitOutput      *output,
		    		CoarseFitInput       *input,
		    		CoarseFitParams      *params )
		

{
  /******* DECLARE VARIABLES ************/

  UINT4 			n,i,j;		/* integer indices */
  REAL8 			psi;
  REAL8 			h0;
  REAL8 			cosIota;
  REAL8			phase;
  REAL8   		chiSquare;
  LALDetector 		detector;
  LALSource 		pulsar;
  LALDetAndSource 	detAndSource;
  LALDetAMResponse 	computedResponse;
  REAL4Vector 		*Fp, *Fc;	/* pointers to amplitude responses of the detector */
  REAL8Vector		*var;
  REAL8 			cos2phase,sin2phase;
  COMPLEX16		Xp, Xc;
  REAL8 			Y, cosIota2;
  COMPLEX16		A,B;
  REAL8			sumAA, sumBB, sumAB;	
  REAL8			h0BestFit=0,phaseBestFit=0, psiBestFit=0, cosIotaBestFit=0;
  REAL8			minChiSquare;
  UINT4			iH0, iCosIota, iPhase, iPsi, arg;
  LALGPSandAcc		pGPSandAcc;
	
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  /******* CHECK VALIDITY OF ARGUMENTS  ************/	 
  ASSERT(output != (CoarseFitOutput *)NULL, status,
         FITTOPULSARH_ENULLOUTPUT, FITTOPULSARH_MSGENULLOUTPUT);	 
  
  ASSERT(params != (CoarseFitParams *)NULL, status,
         FITTOPULSARH_ENULLPARAMS, FITTOPULSARH_MSGENULLPARAMS);
	 
  ASSERT(input != (CoarseFitInput *)NULL, status,
         FITTOPULSARH_ENULLINPUT, FITTOPULSARH_MSGENULLINPUT);

  ASSERT(input->B->length == input->var->length, status, 
  FITTOPULSARH_EVECSIZE , FITTOPULSARH_MSGEVECSIZE );

  /******* EXTRACT INPUTS AND PARAMETERS ************/
  
  n = input->B->length;

  for (i = 0; i < n; i++)
     ASSERT(input->var->data[i].re > 0.0 && input->var->data[i].im > 0.0, status,
         FITTOPULSARH_EVAR, FITTOPULSARH_MSGEVAR);
  
  detector = params->detector;
  detAndSource.pDetector = &detector;
  
  pulsar = params->pulsarSrc;

  /******* ALLOCATE MEMORY TO VECTORS **************/
  
  Fp = NULL;
  LALSCreateVector(status->statusPtr, &Fp, n);
  Fp->length = n;
 
  Fc = NULL;
  LALSCreateVector(status->statusPtr, &Fc, n);
  Fc->length = n;
  
  var = NULL;
  LALDCreateVector(status->statusPtr, &var, n);
  var->length = n; 
  /******* INITIALIZE VARIABLES *******************/
  
  minChiSquare = INICHISQU;
  
  /******* DO ANALYSIS ************/
  for (iPsi = 0; iPsi < params->meshPsi[2]; iPsi++)
  {
    fprintf(stderr,"%d out of %f iPsi\n", iPsi+1, params->meshPsi[2]);
    fflush(stderr);
    psi = params->meshPsi[0] + (REAL8)iPsi*params->meshPsi[1];
    pulsar.orientation = psi;
    detAndSource.pSource = &pulsar;
   
   for (i=0;i<n;i++)
   {    
      Fp->data[i] = 0.0;
      Fc->data[i] = 0.0;
      for(j=0;j<input->N;j++)
      {
        pGPSandAcc.gps.gpsNanoSeconds =  input->t[i].gpsNanoSeconds;
        pGPSandAcc.gps.gpsSeconds =  input->t[i].gpsSeconds + (double)j*60.0; 
        pGPSandAcc.accuracy = 1.0; 
        LALComputeDetAMResponse(status->statusPtr, &computedResponse,&detAndSource, &pGPSandAcc);        
	Fp->data[i] += (REAL8)computedResponse.plus;
        Fc->data[i] += (REAL8)computedResponse.cross;	  
      }	  
      Fp->data[i] /= (double)input->N;
      Fc->data[i] /= (double)input->N;
   }
   for (iPhase = 0; iPhase < params->meshPhase[2]; iPhase++)
   {   
     phase = params->meshPhase[0] + (REAL8)iPhase*params->meshPhase[1];
     cos2phase = cos(2.0*phase);
     sin2phase = sin(2.0*phase);
     for (iCosIota = 0; iCosIota < params->meshCosIota[2]; iCosIota++)
     {
       cosIota = params->meshCosIota[0] + (REAL8)iCosIota*params->meshCosIota[1];
       cosIota2 = 1.0 + cosIota*cosIota; 
       Xp.re = 0.25*cosIota2 * cos2phase;
       Xp.im = 0.25*cosIota2 * sin2phase;

       Y = 0.5*cosIota;
		
       Xc.re = Y*sin2phase;
       Xc.im = -Y*cos2phase; 
	  
       sumAB = 0.0;
       sumAA = 0.0;
       sumBB = 0.0;
       
       for (i = 0; i < n; i++)
       {
	 B.re = input->B->data[i].re;
	 B.im = input->B->data[i].im;
	    
	 A.re = Fp->data[i]*Xp.re + Fc->data[i]*Xc.re;
	 A.im = Fp->data[i]*Xp.im + Fc->data[i]*Xc.im;	
	
	 sumBB += B.re*B.re/input->var->data[i].re +
	 B.im*B.im/input->var->data[i].im;
	 sumAA += A.re*A.re/input->var->data[i].re +
	 A.im*A.im/input->var->data[i].im;
	 sumAB += B.re*A.re/input->var->data[i].re +
	 B.im*A.im/input->var->data[i].im;
       }

	for (iH0 = 0; iH0 < params->meshH0[2]; iH0++)
        {
	  h0 = params->meshH0[0] + (float)iH0*params->meshH0[1];
	  chiSquare = sumBB - 2.0*h0*sumAB + h0*h0*sumAA;
	
	  if (chiSquare<minChiSquare)
	  {
	    minChiSquare = chiSquare;
	    h0BestFit = h0;
	    cosIotaBestFit = cosIota;
	    psiBestFit = psi;
	    phaseBestFit = phase;
	  }
	  
	  arg = iH0 + params->meshH0[2]*(iCosIota +  params->meshCosIota[2]*(iPhase + params->meshPhase[2]*iPsi));
	  output->mChiSquare->data[arg] = chiSquare;
        }
      }  
    }
  }

  /******* CONSTRUCT OUTPUT ************/

  output->h0 = h0BestFit;
  output->cosIota = cosIotaBestFit;
  output->phase = phaseBestFit;
  output->psi = psiBestFit;
  output->chiSquare = minChiSquare;

  ASSERT(minChiSquare < INICHISQU,  status,
  FITTOPULSARH_EMAXCHI , FITTOPULSARH_MSGEMAXCHI); 
  
  LALSDestroyVector(status->statusPtr, &Fp);
  LALSDestroyVector(status->statusPtr, &Fc);
  LALDDestroyVector(status->statusPtr, &var);
  DETATCHSTATUSPTR(status);
  RETURN(status);  
}


		
void
LALFitToPulsarStudentT	( 	LALStatus            *status,
		    		CoarseFitOutput      *output,
		    		FitInputStudentT       *input,
		    		CoarseFitParams      *params )		

{
  /******* DECLARE VARIABLES ************/

	UINT4 			n,i,j,k;		
	REAL8 			psi;
	REAL8 			h0;
	REAL8 			cosIota;
	REAL8			phase;
	REAL8   		chiSquare;
	LALDetector 		detector;
	LALSource 		pulsar;
	LALDetAndSource 	detAndSource;
 	LALDetAMResponse 	computedResponse;
	REAL4Vector 		*Fp, *Fc;	/* pointers to amplitude responses of the detector */
	REAL8 			cos2phase,sin2phase;
 	COMPLEX16		Xp, Xc;
	REAL8 			Y, cosIota2;
	COMPLEX16		A,B;
  	REAL8			sumAA, sumBB, sumAB;	
	REAL8			h0BestFit=0,phaseBestFit=0, psiBestFit=0, cosIotaBestFit=0;
  	REAL8			minChiSquare;
	UINT4			iH0, iCosIota, iPhase, iPsi, arg;
        LALGPSandAcc		pGPSandAcc;

	
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  /******* CHECK VALIDITY OF ARGUMENTS  ************/	 
	 
  
  ASSERT(params != (CoarseFitParams *)NULL, status,
         FITTOPULSARH_ENULLPARAMS, FITTOPULSARH_MSGENULLPARAMS);
	 
  ASSERT(input != (FitInputStudentT *)NULL, status,
         FITTOPULSARH_ENULLINPUT, FITTOPULSARH_MSGENULLINPUT);

  /******* EXTRACT INPUTS AND PARAMETERS ************/
  
  n = input->B->length;

  detector = params->detector;
  detAndSource.pDetector = &detector;
  
  pulsar = params->pulsarSrc;

  /******* ALLOCATE MEMORY TO VECTORS **************/
  
  Fp = NULL;
  LALSCreateVector(status->statusPtr, &Fp, n);
  Fp->length = n;
 
  Fc = NULL;
  LALSCreateVector(status->statusPtr, &Fc, n);
  Fc->length = n;

  /******* INITIALIZE VARIABLES *******************/
  for (i=0;i<params->meshH0[2]*params->meshCosIota[2]*params->meshPhase[2]*params->meshPsi[2];i++) 
     output->mChiSquare->data[i] = 0.0;

  
  minChiSquare = INICHISQU;
  
 /******* DO ANALYSIS ************/
  for (iPsi = 0; iPsi < params->meshPsi[2]; iPsi++)
  {
    fprintf(stderr,"%d out of %f iPsi\n", iPsi+1, params->meshPsi[2]);
    fflush(stderr);
    psi = params->meshPsi[0] + (REAL8)iPsi*params->meshPsi[1];
    pulsar.orientation = psi;
    detAndSource.pSource = &pulsar;
   
   /* create vectors containing amplitude response of detector for specified times */ 
   for (i=0;i<n;i++)
   { 
     pGPSandAcc.gps.gpsNanoSeconds =  input->t[i].gpsNanoSeconds;
     pGPSandAcc.gps.gpsSeconds =  input->t[i].gpsSeconds;
     pGPSandAcc.accuracy = 1.0; 
     LALComputeDetAMResponse(status->statusPtr, &computedResponse,&detAndSource, &pGPSandAcc);          
     Fp->data[i] = (REAL8)computedResponse.plus;
     Fc->data[i] = (REAL8)computedResponse.cross;	
   }
   for (iPhase = 0; iPhase < params->meshPhase[2]; iPhase++)
   {   
     phase = params->meshPhase[0] + (REAL8)iPhase*params->meshPhase[1];
     cos2phase = cos(2.0*phase);
     sin2phase = sin(2.0*phase);
     for (iCosIota = 0; iCosIota < params->meshCosIota[2]; iCosIota++)
     {
       cosIota = params->meshCosIota[0] + (REAL8)iCosIota*params->meshCosIota[1];
       cosIota2 = 1.0 + cosIota*cosIota; 
       Xp.re = 0.25*cosIota2 * cos2phase;
       Xp.im = 0.25*cosIota2 * sin2phase;

       Y = 0.5*cosIota;
		
       Xc.re = Y*sin2phase;
       Xc.im = -Y*cos2phase; 
	  
       sumAB = 0.0;
       sumAA = 0.0;
       sumBB = 0.0;

        k = input->N;     
       for (i = 0; i < n-k; i+=k)
       { 
       
	sumBB = 0.0;
	sumAA = 0.0;
	sumAB = 0.0;
	
         for (j=i;j<i+k;j++)
	 {
	 B.re = input->B->data[j].re;
	 B.im = input->B->data[j].im;
	    
	 A.re = Fp->data[j]*Xp.re + Fc->data[j]*Xc.re;
	 A.im = Fp->data[j]*Xp.im + Fc->data[j]*Xc.im;	
	
	 sumBB += B.re*B.re + B.im*B.im;
	 sumAA += A.re*A.re + A.im*A.im;
	 sumAB += B.re*A.re + B.im*A.im;
         }

	for (iH0 = 0; iH0 < params->meshH0[2]; iH0++)
        {
	  h0 = params->meshH0[0] + (float)iH0*params->meshH0[1]; 
	  chiSquare = sumBB - 2.0*h0*sumAB + h0*h0*sumAA;
	  arg = iH0 + params->meshH0[2]*(iCosIota +  params->meshCosIota[2]*(iPhase + params->meshPhase[2]*iPsi));
	  output->mChiSquare->data[arg] += 2.0*(k)*log(chiSquare); 
        } /* iH0 */

       } /* n*/
       /* find best fit */
       	for (iH0 = 0; iH0 < params->meshH0[2]; iH0++)
        {
	  arg = iH0 + params->meshH0[2]*(iCosIota +  params->meshCosIota[2]*(iPhase + params->meshPhase[2]*iPsi));
	  if ((output->mChiSquare->data[arg])<minChiSquare)
	  {
	     minChiSquare = output->mChiSquare->data[arg];
	     h0 = params->meshH0[0] + (float)iH0*params->meshH0[1]; 
	     h0BestFit = h0;
	     cosIotaBestFit = cosIota;
	     psiBestFit = psi;
	     phaseBestFit = phase;	 
	  }	  
        }  /* end find best fit */
      }  /* iCosIota */
    } /* iPhase */
  } /* iPsi */

  
  /******* CONSTRUCT OUTPUT ************/
  output->h0 = h0BestFit;
  output->cosIota = cosIotaBestFit;
  output->phase = phaseBestFit;
  output->psi = psiBestFit;
  output->chiSquare = minChiSquare;

  ASSERT(minChiSquare < INICHISQU,  status,
  FITTOPULSARH_EMAXCHI , FITTOPULSARH_MSGEMAXCHI); 
  
  LALSDestroyVector(status->statusPtr, &Fp);
  LALSDestroyVector(status->statusPtr, &Fc);
  DETATCHSTATUSPTR(status);
  RETURN(status);
 
  
}

/* ******************************* */

void
LALPulsarMarginalize   ( 	LALStatus              *status,
		    		PulsarPdfs             *output,
		    		CoarseFitOutput        *input,
		    		CoarseFitParams        *params )

{
  /******* DECLARE VARIABLES ************/
   UINT4 arg,i; 
   UINT4 iH0, iPsi, iPhase, iCosIota; 
   REAL8 area;
  
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  /* initialize outputs */
  for (i=0;i<params->meshH0[2];i++)
  {
    output->pdf->data[i] = 0.0;
    output->cdf->data[i] = 0.0;  
  }
  for (iPhase = 0; iPhase < params->meshPhase[2]; iPhase++)
  {
    output->pdfPhase->data[iPhase] = 0.0;
    output->cdfPhase->data[iPhase] = 0.0; 
  }
  for (iPsi = 0; iPsi < params->meshPsi[2]; iPsi++)
  {
    output->pdfPsi->data[iPsi] = 0.0;
    output->cdfPsi->data[iPsi] = 0.0; 
  }
  for (iCosIota = 0; iCosIota < params->meshCosIota[2];iCosIota++)   
  {
    output->pdfCosIota->data[iCosIota] = 0.0;
    output->cdfCosIota->data[iCosIota] = 0.0; 
  }

 /* marginalize over psi, phase, cosIota for p(h_0|B_k) */
 for (iPsi = 0; iPsi < params->meshPsi[2]; iPsi++)
   for (iPhase = 0; iPhase < params->meshPhase[2]; iPhase++)
     for (iCosIota = 0; iCosIota < params->meshCosIota[2];iCosIota++)
       for (iH0 = 0; iH0 < params->meshH0[2]; iH0++)
       {
         arg = iH0 + params->meshH0[2]*(iCosIota +  params->meshCosIota[2]*(iPhase + params->meshPhase[2]*iPsi)); 	 
         output->pdf->data[iH0] += exp((input->chiSquare - input->mChiSquare->data[arg])/2.0); 	      
       }
 
  area = 0.0;
  for (iH0 = 0; iH0 < params->meshH0[2]; iH0++)
    area += output->pdf->data[iH0]*params->meshH0[1];
  
  for (iH0 = 0; iH0 < params->meshH0[2]; iH0++)
    output->pdf->data[iH0] = output->pdf->data[iH0]/area;
   
  output->cdf->data[0] = output->pdf->data[0]*params->meshH0[1];  
  for (iH0 = 1; iH0 < params->meshH0[2]; iH0++)
    output->cdf->data[iH0] = output->pdf->data[iH0]*params->meshH0[1] + output->cdf->data[iH0-1];
  
        
 /* marginalize over h0, phase, cosIota for p(psi|B_k) */
 for (iH0 = 0; iH0 < params->meshH0[2]; iH0++)  
   for (iPhase = 0; iPhase < params->meshPhase[2]; iPhase++)   
     for (iCosIota = 0; iCosIota < params->meshCosIota[2];iCosIota++)     
       for (iPsi = 0; iPsi < params->meshPsi[2]; iPsi++)
       {
         arg = iH0 + params->meshH0[2]*(iCosIota +  params->meshCosIota[2]*(iPhase + params->meshPhase[2]*iPsi)); 	 
         output->pdfPsi->data[iPsi] += exp((input->chiSquare - input->mChiSquare->data[arg])/2.0); 	      
       }  
  
  area = 0.0;
  for (iPsi = 0; iPsi < params->meshPsi[2]; iPsi++)
    area += output->pdfPsi->data[iPsi]*params->meshPsi[1];
 
  for (iPsi = 0; iPsi < params->meshPsi[2]; iPsi++)
    output->pdfPsi->data[iPsi] = output->pdfPsi->data[iPsi]/area;
  
  output->cdfPsi->data[0] = output->pdfPsi->data[0]*params->meshPsi[1];
  for (iPsi = 1; iPsi < params->meshPsi[2]; iPsi++)  
    output->cdfPsi->data[iPsi] = output->pdfPsi->data[iPsi]*params->meshPsi[1] + output->cdfPsi->data[iPsi-1];
 
/* marginalize over psi, h0, cosIota for p(phase|B_k) */
 for (iPsi = 0; iPsi < params->meshPsi[2]; iPsi++)
   for (iH0 = 0; iH0 < params->meshH0[2]; iH0++)
     for (iCosIota = 0; iCosIota < params->meshCosIota[2];iCosIota++)
       for (iPhase = 0; iPhase < params->meshPhase[2]; iPhase++)
       {
         arg = iH0 + params->meshH0[2]*(iCosIota +  params->meshCosIota[2]*(iPhase + params->meshPhase[2]*iPsi)); 	 
         output->pdfPhase->data[iPhase] += exp((input->chiSquare - input->mChiSquare->data[arg])/2.0); 	      
       }
 
  area = 0.0;
  for (iPhase = 0; iPhase < params->meshPhase[2]; iPhase++)
    area += output->pdfPhase->data[iPhase]*params->meshPhase[1];
  
  for (iPhase = 0; iPhase < params->meshPhase[2]; iPhase++)
    output->pdfPhase->data[iPhase] = output->pdfPhase->data[iPhase]/area;
   
  output->cdfPhase->data[0] = output->pdfPhase->data[0]*params->meshPhase[1];  
  for (iPhase = 1; iPhase < params->meshPhase[2]; iPhase++)
    output->cdfPhase->data[iPhase] = output->pdfPhase->data[iPhase]*params->meshPhase[1] + output->cdfPhase->data[iPhase-1];

/* marginalize over h0, phase, psi for p(cosIota|B_k) */
 for (iPsi = 0; iPsi < params->meshPsi[2]; iPsi++)
   for (iPhase = 0; iPhase < params->meshPhase[2]; iPhase++)
     for (iH0 = 0; iH0 < params->meshH0[2];iH0++)
       for (iCosIota = 0; iCosIota < params->meshCosIota[2]; iCosIota++)
       {
         arg = iH0 + params->meshH0[2]*(iCosIota +  params->meshCosIota[2]*(iPhase + params->meshPhase[2]*iPsi)); 	 
         output->pdfCosIota->data[iCosIota] += exp((input->chiSquare - input->mChiSquare->data[arg])/2.0); 	      
       }
 
  area = 0.0;
  for (iCosIota = 0; iCosIota < params->meshCosIota[2]; iCosIota++)
    area += output->pdfCosIota->data[iCosIota]*params->meshCosIota[1];
  
  for (iCosIota = 0; iCosIota < params->meshCosIota[2]; iCosIota++)
    output->pdfCosIota->data[iCosIota] = output->pdfCosIota->data[iCosIota]/area;
   
  output->cdfCosIota->data[0] = output->pdfCosIota->data[0]*params->meshCosIota[1];  
  for (iCosIota = 1; iCosIota < params->meshCosIota[2]; iCosIota++)
    output->cdfCosIota->data[iCosIota] = output->pdfCosIota->data[iCosIota]*params->meshCosIota[1] +
    output->cdfCosIota->data[iCosIota-1];
  

  DETATCHSTATUSPTR(status);
  RETURN(status);

}				



/*
*  Copyright (C) 2007 Chad Hanna, Duncan Brown
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

/*-----------------------------------------------------------------------
 *
 * File Name: FindChirpBankVeto.c
 *
 * Author: Brown D. A., and Hanna, C.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpBankVetoCV">
Author: Brown, D. A., and Hanna, C.
$Id$
</lalVerbatim>

<lalLaTeX>
\vfill{\footnotesize\input{FindChirpBankVetoCV}}
</lalLaTeX>
#endif

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <time.h>


NRCSID (FINDCHIRPBANKVETOC, "$Id$");

/* macro to "use" unused function parameters */
#define UNUSED(expr) do { (void)(expr); } while(0)

/* Some static function prototypes */

static int compareTemplateByLevel (const void * a, const void * b);
static int compareTemplateByChirpMass (const void * a, const void * b);
static int compareTemplateByMass (const void * a, const void * b);
static int compareTemplateByEta (const void * a, const void * b);
static int  breakUpRegions(InspiralTemplate *bankHead);
static InspiralTemplate *
XLALFindChirpSortTemplatesByChirpMass( InspiralTemplate *bankHead, UINT4 num );
static InspiralTemplate *
XLALFindChirpSortTemplatesByMass( InspiralTemplate *bankHead, UINT4 num );
static InspiralTemplate *
XLALFindChirpSortTemplatesByLevel( InspiralTemplate *bankHead, UINT4 num);

static int writeCCmatToFile(COMPLEX8Vector *ccmat,UINT4 maxSubBankSize,FILE *fdOut)
{
    int row;
    int col;

    fprintf(stderr,"in");

    for ( row = 0; row< maxSubBankSize; row++ )
    {
	for ( col = 0; col < maxSubBankSize; col++ )
	{
	    fprintf(fdOut,"%.8g+i%.8g\t",ccmat->data[row*maxSubBankSize+col].re,ccmat->data[row*maxSubBankSize+col].im);
	}

	fprintf(fdOut,"\n");

    }

    fprintf(stderr,"out");

    return 0;
}



static double chirp_time (double m1, double m2, double fLower, int order)
	{

	/* variables used to compute chirp time */
	double c0T, c2T, c3T, c4T, c5T, c6T, c6LogT, c7T;
	double xT, x2T, x3T, x4T, x5T, x6T, x7T, x8T;
	double m = m1 + m2;
	double eta = m1 * m2 / m / m;

	c0T = c2T = c3T = c4T = c5T = c6T = c6LogT = c7T = 0.;

	/* Switch on PN order, set the chirp time coeffs for that order */
	switch (order)
		{
		case 8:
		case 7:
			c7T = LAL_PI * (14809.0 * eta * eta - 75703.0 * eta / 756.0 - 15419335.0 / 127008.0);
		case 6:
			c6T = LAL_GAMMA * 6848.0 / 105.0 - 10052469856691.0 / 23471078400.0 + LAL_PI * LAL_PI * 128.0 / 3.0 + eta * (3147553127.0 / 3048192.0 - LAL_PI * LAL_PI * 451.0 / 12.0) - eta * eta * 15211.0 / 1728.0 + eta * eta * eta * 25565.0 / 1296.0 + log (4.0) * 6848.0 / 105.0;
     			c6LogT = 6848.0 / 105.0;
		case 5:
			c5T = 13.0 * LAL_PI * eta / 3.0 - 7729.0 / 252.0;
		case 4:
			c4T = 3058673.0 / 508032.0 + eta * (5429.0 / 504.0 + eta * 617.0 / 72.0);
			c3T = -32.0 * LAL_PI / 5.0;
			c2T = 743.0 / 252.0 + eta * 11.0 / 3.0;
			c0T = 5.0 * m * LAL_MTSUN_SI / (256.0 * eta);	
			break;
		default:
			fprintf (stderr, "ERROR!!!\n");
			break;
		}

	/* This is the PN parameter v evaluated at the lower freq. cutoff */
	xT = pow (LAL_PI * m * LAL_MTSUN_SI * fLower, 1.0 / 3.0);
	x2T = xT * xT;
	x3T = xT * x2T;
	x4T = x2T * x2T;
	x5T = x2T * x3T;
	x6T = x3T * x3T;
	x7T = x3T * x4T;
	x8T = x4T * x4T;

	/* Computes the chirp time as tC = t(v_low)    */
	/* tC = t(v_low) - t(v_upper) would be more    */
	/* correct, but the difference is negligble.   */

	/* This formula works for any PN order, because */
	/* higher order coeffs will be set to zero.     */

	return c0T * (1 + c2T * x2T + c3T * x3T + c4T * x4T + c5T * x5T + (c6T + c6LogT * log (xT)) * x6T + c7T * x7T) / x8T;
}


/* A convenience function to compute the time between two frequencies */
static double chirp_time_between_f1_and_f2(double m1, double m2, double fLower, double fUpper, int order)
	{
	return chirp_time(m1,m2,fLower,order) - chirp_time(m1,m2,fUpper,order);
	}


void XLALInitBankVetoData(FindChirpBankVetoData *bankVetoData)
{
  bankVetoData->acorr = NULL;
  bankVetoData->workspace = NULL;
  bankVetoData->acorrMat = NULL;
  bankVetoData->revplan = NULL;
  bankVetoData->tchirp = NULL;
}

void XLALDestroyBankVetoData(FindChirpBankVetoData *bvdata)
{
  if (bvdata->acorr) XLALDestroyREAL4Vector(bvdata->acorr);
  if (bvdata->workspace) XLALDestroyCOMPLEX8Vector(bvdata->workspace);
  if (bvdata->acorrMat) XLALDestroyREAL4Vector(bvdata->acorrMat);
  if (bvdata->revplan) XLALDestroyREAL4FFTPlan(bvdata->revplan);
  if (bvdata->tchirp) XLALDestroyREAL4Vector(bvdata->tchirp);
  return; 
}


/* <lalVerbatim file="FindChirpBankVetoCP"> */

/*
|
|  XLALFindChirpCreateSubBanks takes a linked list of InspiralTemplates
|  and breaks the templates up into groups of roughly the same size.
|  Each group is called a subbank.  The subbanks are disjoint so that 
|  each template belongs to exactly one subbank.  The function returns
|  a pointer to the first element
|
|
*/

FindChirpSubBank*
XLALFindChirpCreateSubBanks(
    UINT4                      *maxSubBankSize,
    UINT4                       subBankSize,
    UINT4                       bankSize,
    InspiralTemplate           *bankHead
    )
/* </lalVerbatim> */
{
  UINT4                 i;
  UINT4                 numSubBanks = 0;
  UINT4                 subBankRemainder = 0;
  UINT4                *bankSizes = NULL;
  FindChirpSubBank     *subBankHead = NULL;
  FindChirpSubBank     *thisSubBank = NULL;
  InspiralTemplate     *thisTmplt = NULL;
  InspiralTemplate     *nextTmplt = NULL;

  numSubBanks = bankSize / subBankSize;
  subBankRemainder = bankSize % subBankSize;

  /* if there are no templates return NULL */
  if ( ! bankSize )
    return NULL;

  if ( ! numSubBanks )
  {
    /* the bank is smaller than the subbank size, so return the entire */
    /* template bank as the subbank                                    */
    subBankHead = (FindChirpSubBank *) LALCalloc( 1, sizeof(FindChirpSubBank) );
    subBankHead->bankHead = bankHead;
    subBankHead->subBankSize = bankSize;
    *maxSubBankSize = bankSize;
    return subBankHead;
  }

  /* determine whether all subBanks are the same size */
  *maxSubBankSize = subBankSize;
  if ( subBankRemainder )
  {
    *maxSubBankSize = subBankSize + 1; 
  }

  /* create an array of subbank sizes with the minimum size */
  bankSizes = (UINT4 *) LALCalloc( numSubBanks, sizeof(UINT4) );

  for ( i = 0; i < numSubBanks; ++i )
  {
    bankSizes[i] = subBankSize;
  }

  /* disperse the remainder through the subbanks */
  while( subBankRemainder )
  {
    for ( i = 0; i < numSubBanks; ++i )
    {
      if ( ! subBankRemainder )
      {
        break;
      }
      ++bankSizes[i];
      --subBankRemainder;
    }
  }

  /* allocate storage for the subbanks */
  for ( i = 0; i < numSubBanks; ++i )
  {
    if ( ! subBankHead )
    {
      thisSubBank = subBankHead =
        (FindChirpSubBank *) LALCalloc( 1, sizeof(FindChirpSubBank) );
    }
    else
    {
      thisSubBank = thisSubBank->next =
        (FindChirpSubBank *) LALCalloc( 1, sizeof(FindChirpSubBank) );
    }

    thisSubBank->subBankSize = bankSizes[i];

  }

  /* chop up the template bank into subbanks */
  for ( thisSubBank = subBankHead, nextTmplt = bankHead; thisSubBank;
      thisSubBank = thisSubBank->next )
  {
    thisTmplt = nextTmplt;
    thisSubBank->bankHead = thisTmplt;
    for ( i = 0; i < thisSubBank->subBankSize - 1; ++i )
    {
      thisTmplt = thisTmplt->next;
    }
    nextTmplt = thisTmplt->next;
    thisTmplt->next = NULL;
  }

  LALFree( bankSizes );

  return subBankHead;
}


REAL4
XLALComputeFullChisq(
    FindChirpBankVetoData      *bankVetoData,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params,
    COMPLEX8                   *q,
    UINT4 			i,
    UINT4                       snrIX,
    UINT4                      *dof,
    REAL4                       norm
)

{
  UINT4 k;
  REAL4 chisqnorm, tmp, C, chisq, angle, snri, snrk;
  /* input, params and norm are unused */
  UNUSED(input);
  UNUSED(params);
  /* test isn't being done */

  if (!bankVetoData->acorrMat) return 0.0;
  chisq = 0;
  chisqnorm = 0;

  angle = atan2(q[snrIX].im, q[snrIX].re);
  snri = q[snrIX].re * cos(angle) + q[snrIX].im * sin(angle);
  *dof = 0;

  for (k = 0; k < bankVetoData->acorrMatSize; k++)
  {
    /* we skipped saving the first sample of the autocorrelation which is why */
    /* we need the -1 in the snrk here */
    C = bankVetoData->acorrMat->data[i * bankVetoData->acorrMatSize + k];
    chisqnorm = 1.0 - C * C;
    snrk = q[snrIX-k-1].re * cos(angle) + q[snrIX-k-1].im * sin(angle);
    (*dof) += 1;
    tmp = snrk - C * snri;
    chisq += tmp * tmp * norm / chisqnorm;
  }

  /*chisq /= chisqnorm;*/

  return chisq;  
}

void
XLALBankVetoCCMat ( FindChirpBankVetoData *bankVetoData,
                    FindChirpSubBank *vetoBank,
                    FindChirpDataParams *params,
                    REAL4 dynRange,
                    REAL4 fLow,
                    REAL4 deltaF)
{
    /* Temporary variables that hold the real and imaginary parts
     * of the two templates, which will be multiplied together
     * and added to the cross-correlation integral */
    REAL8 tmpltAreal = 0;  
    REAL8 tmpltAimag = 0;    
    REAL8 tmpltBreal = 0;    
    REAL8 tmpltBimag = 0;  

    /* The cumulative real and imaginary parts of the cross-correlation
     * integral int( A*(f)B(f)/Sn(f) )df */
    REAL8 crossCorrABreal; 
    REAL8 crossCorrABimag;   

    /* Indices used to iterate of CC matrix row and columns */
    UINT4 row; 
    UINT4 col;   

    /* Index used to integrate over frequencies */
    UINT4 sample;           
    
    /* Indices which dictate the start and end of integration */
    UINT4 startIndex = floor( fLow / deltaF ); 
    UINT4 sampleMaxTmpltA;  // Maximum index available in template A
    UINT4 sampleMaxTmpltB;  // Maximum index available in template B
    UINT4 templateLength = bankVetoData->fcInputArray[0]->fcTmplt->data->length;

    /* Parameters used to time shift templates */
    REAL4 phaseOffset;        // phase shift for template A
    double fBucket = 150.0;   // approximate frequency of the PSD bucket
    REAL4 deltaChirp;


    /* FIXME what do these do? */
    REAL8 sqResp, normfac; 

    /* FIXME this should be a command line argument */
    bankVetoData->acorrMatSize = 200; /*200 points of autocorrelation function stored */

    /* Allocate memory for workspace and the autocorrelation if necessary */
    if ( !bankVetoData->acorr ) bankVetoData->acorr = XLALCreateREAL4Vector((templateLength-1) * 2);
    if ( !bankVetoData->workspace) bankVetoData->workspace = XLALCreateCOMPLEX8Vector(templateLength);
    if ( !bankVetoData->acorrMat) bankVetoData->acorrMat = XLALCreateREAL4Vector(bankVetoData->acorrMatSize * bankVetoData->length);
    if ( !bankVetoData->revplan) bankVetoData->revplan = XLALCreateReverseREAL4FFTPlan((templateLength-1) * 2 , 0);
    if ( !bankVetoData->tchirp) bankVetoData->tchirp = XLALCreateREAL4Vector(bankVetoData->length);

    /* Compute chirp times for time shifting */
    for (row = 0; row < vetoBank->subBankSize; row++ )
    {
	bankVetoData->tchirp->data[row] = (REAL4) ( chirp_time_between_f1_and_f2(bankVetoData->fcInputArray[row]->fcTmplt->tmplt.mass1, 
									   bankVetoData->fcInputArray[row]->fcTmplt->tmplt.mass2,
									   fBucket,
									   bankVetoData->fcInputArray[row]->fcTmplt->tmplt.fFinal, 7) );
    }		
	

    /* Compute the cross-correlation between template row and template col */
    for ( row = 0; row < bankVetoData->length; row++ )
    {
	
	for ( col = row; col < bankVetoData->length; col++ )
	{          
	    if ( row == col  )
	    {
		memset(bankVetoData->workspace->data, 0, bankVetoData->workspace->length * sizeof(COMPLEX8));
	    }


	    fprintf(stderr,"i=%d,j=%d\n",row,col);
	    /* initialize the cross-correlation integral */
	    crossCorrABreal = crossCorrABimag = 0;
	    	    
	    if ( (row < vetoBank->subBankSize) && (col < vetoBank->subBankSize) )
	    {
		/* if the templates are the same then do an autocorrelation function */
		
		/* find the ending sample point of the templates */
		sampleMaxTmpltA = bankVetoData->fcInputArray[row]->fcTmplt->tmplt.fFinal / deltaF;
		sampleMaxTmpltB = bankVetoData->fcInputArray[col]->fcTmplt->tmplt.fFinal / deltaF;
		
		/* Begin computation of the cross-correlation. */	  
		for ( sample = startIndex; sample < templateLength-1; sample++ )
		{
		    if ( (sample > sampleMaxTmpltA) || (sample > sampleMaxTmpltB) ) break;
		    sqResp = ( bankVetoData->resp->data[sample].re *
			       bankVetoData->resp->data[sample].re +
			       bankVetoData->resp->data[sample].im *
			       bankVetoData->resp->data[sample].im ) * dynRange * dynRange;
		    
		    if ( bankVetoData->spec->data[sample] != 0 && sqResp != 0 )
		    {

			deltaChirp = bankVetoData->tchirp->data[col] - bankVetoData->tchirp->data[row];
			phaseOffset = 2*LAL_PI*deltaF*sample*deltaChirp;

			
			/*  
			 *  tmpltA* x tmpltB = (tmpltAreal-i*tmpltAimag) x (tmpltBreal+i*tmpltBimag) =
			 *     (tmpltAreal*tmpltBreal + tmpltAimag*tmpltBimag) +
			 *   i*(tmpltAreal*tmpltBimag - tmpltAimag*tmpltBreal)
			 *
			 */

			tmpltAreal = bankVetoData->fcInputArray[col]->fcTmplt->data->data[sample].re;
			tmpltAimag = bankVetoData->fcInputArray[col]->fcTmplt->data->data[sample].im;
			
			/* Time shift tmplt A by deltaChirp */
			tmpltAreal = tmpltAreal*cos(phaseOffset) - tmpltAimag*sin(phaseOffset);
			tmpltAimag = tmpltAimag*cos(phaseOffset) + tmpltAreal*sin(phaseOffset);
			
			tmpltBreal = bankVetoData->fcInputArray[row]->fcTmplt->data->data[sample].re /
			  ( bankVetoData->spec->data[sample] * sqResp ) *
			  params->ampVec->data[sample] * params->ampVec->data[sample];
			
			tmpltBimag = bankVetoData->fcInputArray[row]->fcTmplt->data->data[sample].im /
			  ( bankVetoData->spec->data[sample] * sqResp ) *
			  params->ampVec->data[sample] * params->ampVec->data[sample];
			
			crossCorrABreal +=    (tmpltAreal*tmpltBreal + tmpltAimag*tmpltBimag);
			
			crossCorrABimag +=    (tmpltAreal*tmpltBimag - tmpltAimag*tmpltBreal);
			
		    }
		    
		    if ( row == col )
		    {	      
			bankVetoData->workspace->data[sample].re = (REAL4) tmpltBreal * bankVetoData->fcInputArray[row]->fcTmplt->data->data[sample].re 
			  + tmpltBimag * bankVetoData->fcInputArray[row]->fcTmplt->data->data[sample].im;
			bankVetoData->workspace->data[sample].im = 0.0; /* must be 0 */
		    }
		    
		}
		
		
		if ( row == col )
		{ 
		    XLALREAL4ReverseFFT( bankVetoData->acorr, bankVetoData->workspace, bankVetoData->revplan );
		    /* since the normalized first sample of the auto correlation is always 1, lets skip it */
		    for (sample=1; sample < bankVetoData->acorrMatSize + 1; sample++)
		    {
			/* Save the last acorrMatSize samples */
			bankVetoData->acorrMat->data[row*bankVetoData->acorrMatSize + sample-1] = bankVetoData->acorr->data[sample] / bankVetoData->acorr->data[0];
		    }
		}

		
	    }
	    
	    /* Fill in the cross-correlation matrix.
	       The matrix is hermitian so the imaginary flips sign. */
	    bankVetoData->ccMat->data[row*bankVetoData->length + col].re = (REAL4) crossCorrABreal;
	    bankVetoData->ccMat->data[col*bankVetoData->length + row].re = (REAL4) crossCorrABreal;
	    bankVetoData->ccMat->data[row*bankVetoData->length + col].im = (REAL4) crossCorrABimag;
	    bankVetoData->ccMat->data[col*bankVetoData->length + row].im = (REAL4) -crossCorrABimag;
	    
	    /* what is normMat ? FIXME */
	    bankVetoData->normMat->data[row*bankVetoData->length + col] = (REAL4) crossCorrABreal;
	    bankVetoData->normMat->data[col*bankVetoData->length + row] = (REAL4) crossCorrABreal;
	    
	}
    } /* end loop over templates */
    

    // HUH?
    for ( row = 0; row < vetoBank->subBankSize; row++ ) 
    {
	for ( col = row; col < vetoBank->subBankSize; col++ ) 
	{
	    normfac = (REAL8) 1.0 / sqrt(  bankVetoData->normMat->data[row*bankVetoData->length + row] * bankVetoData->normMat->data[col*bankVetoData->length + col]);
	    bankVetoData->ccMat->data[row*bankVetoData->length + col].re *= (REAL4) normfac;
	    if (col != row) bankVetoData->ccMat->data[col*bankVetoData->length + row].re *= (REAL4) normfac;
	    bankVetoData->ccMat->data[row*bankVetoData->length + col].im *= (REAL4) normfac;
	    if (col != row) bankVetoData->ccMat->data[col*bankVetoData->length + row].im *= (REAL4) normfac;
	}
    } /* end norm loop */



    if ( 0 )
    {

    /* Output data */
    CHAR ccFileName[100];
    sprintf(ccFileName,"ccmat_%d.txt",time(NULL)-1257471720);
    FILE *fdOut = fopen(ccFileName,"w");

    writeCCmatToFile(bankVetoData->ccMat,bankVetoData->length,fdOut);

    fclose(fdOut);
    /* End output data */
    
    }

}


REAL4
XLALComputeBankVeto( FindChirpBankVetoData *bankVetoData,
                     UINT4 row,
                     UINT4 snrIndexRow,
		     REAL4 deltaT,
                     UINT4 *dof)
{

    REAL8 tmpltA_SNR_real = 0;
    REAL8 tmpltA_SNR_imag = 0;
    REAL8 tmpltA_SNR_mag = 0;
    REAL8 tmpltA_SNR_phi = 0;
    
    REAL8 tmpltB_SNR_real = 0;
    REAL8 tmpltB_SNR_imag = 0;
    
    UINT8 col;
    UINT8 snrIndexCol;

    COMPLEX16 chisq;
    REAL4 chisq_mag = 0;
    REAL4 bankNorm;

    COMPLEX16 crossCorrAB;
    REAL8 crossCorrAB_mag;

    /* Handle trivial case */
    if (bankVetoData->length == 1) return 0;
   
    /* Lookup SNR for template A at time given by snrIndexRow */
    tmpltA_SNR_real = bankVetoData->qVecArray[row]->data[snrIndexRow].re
      * sqrt(bankVetoData->fcInputArray[row]->fcTmplt->norm);
    tmpltA_SNR_imag = bankVetoData->qVecArray[row]->data[snrIndexRow].im
      * sqrt(bankVetoData->fcInputArray[row]->fcTmplt->norm);
    
    tmpltA_SNR_mag = sqrt(tmpltA_SNR_real * tmpltA_SNR_real + tmpltA_SNR_imag * tmpltA_SNR_imag);
    tmpltA_SNR_phi = atan2(tmpltA_SNR_imag, tmpltA_SNR_real);

    chisq_mag = 0;
    *dof = 0;

    for (col = 0; col < bankVetoData->length; col++)
    {
	if (row == col) continue;
	if (bankVetoData->ccMat->data[row*bankVetoData->length + col].re == 0 && bankVetoData->ccMat->data[row*bankVetoData->length + col].im == 0) continue;
	

	tmpltA_SNR_phi = 0; //FIXME
	/* Time shift the cross correlation of template B with template B by tmpltA_SNR_phi */ 
	crossCorrAB.re = bankVetoData->ccMat->data[row*bankVetoData->length + col].re * cos(tmpltA_SNR_phi)
	                 - bankVetoData->ccMat->data[row*bankVetoData->length + col].im * sin(tmpltA_SNR_phi);
	crossCorrAB.im = bankVetoData->ccMat->data[row*bankVetoData->length + col].im * cos(tmpltA_SNR_phi) 
	                 + bankVetoData->ccMat->data[row*bankVetoData->length + col].re * sin(tmpltA_SNR_phi);

	crossCorrAB_mag = crossCorrAB.re * crossCorrAB.re + crossCorrAB.im * crossCorrAB.im;

	bankNorm = 2.0 - crossCorrAB_mag;

	/* FIXME: warning -- this could go off the edge of the SNR time series */
	snrIndexCol = (UINT4) ( (double) snrIndexRow - round( ( bankVetoData->tchirp->data[col] -   bankVetoData->tchirp->data[row]) / deltaT ) );

	tmpltB_SNR_real = bankVetoData->qVecArray[col]->data[snrIndexCol].re
	  * sqrt(bankVetoData->fcInputArray[col]->fcTmplt->norm);
	tmpltB_SNR_imag = bankVetoData->qVecArray[col]->data[snrIndexCol].im
	  * sqrt(bankVetoData->fcInputArray[col]->fcTmplt->norm);
      
	chisq.re = tmpltB_SNR_real - crossCorrAB.re * (tmpltA_SNR_real);
	chisq.im = tmpltB_SNR_imag - crossCorrAB.im * (tmpltA_SNR_imag);
	chisq_mag += (chisq.re*chisq.re + chisq.im*chisq.im) / bankNorm;
	(*dof)++;

    }

    /*FIXME Normalization now not necessarily chisquare distributed */
    return chisq_mag;
    /*return (REAL4) chisq / bankNorm;*/
}


InspiralTemplate *
XLALFindChirpSortTemplates( InspiralTemplate *bankHead, UINT4 num)
  {
  bankHead = XLALFindChirpSortTemplatesByChirpMass(bankHead,num);
  breakUpRegions(bankHead);
  return XLALFindChirpSortTemplatesByLevel(bankHead,num);
  }


static InspiralTemplate *
XLALFindChirpSortTemplatesByChirpMass( InspiralTemplate *bankHead, UINT4 num )
  {
  InspiralTemplate **bankArray = NULL;
  InspiralTemplate *bankFirst = NULL;
  UINT4 i = 0;
  bankFirst = bankHead;
  bankArray = (InspiralTemplate **) LALCalloc(num, sizeof(InspiralTemplate *));

  for (i = 0; (i < num); bankHead = bankHead->next, i++)
  {
    bankArray[i] = bankHead; /* populate pointer array */
  }

  qsort(bankArray, num, sizeof(InspiralTemplate *), compareTemplateByChirpMass);

  bankFirst = bankHead = bankArray[0];
  /* repopulate linked list */
  for (i=1; i < num; i++)
  {
    bankHead = bankHead->next = bankArray[i];
  }
  bankHead->next = NULL;
  LALFree(bankArray);
  return bankFirst;
}

static InspiralTemplate *
XLALFindChirpSortTemplatesByMass( InspiralTemplate *bankHead, UINT4 num )
  {
  InspiralTemplate **bankArray = NULL;
  InspiralTemplate *bankFirst = NULL;
  UINT4 i = 0;
  bankFirst = bankHead;
  bankArray = (InspiralTemplate **) LALCalloc(num, sizeof(InspiralTemplate *));

  for (i = 0; (i < num); bankHead = bankHead->next, i++)
  {
    bankArray[i] = bankHead; /* populate pointer array */
  }

  qsort(bankArray, num, sizeof(InspiralTemplate *), compareTemplateByMass);

  bankFirst = bankHead = bankArray[0];
  /* repopulate linked list */
  for (i=1; i < num; i++)
  {
    bankHead = bankHead->next = bankArray[i];
  }
  bankHead->next = NULL;
  LALFree(bankArray);
  return bankFirst;
}

static InspiralTemplate *
XLALFindChirpSortTemplatesByEta( InspiralTemplate *bankHead, UINT4 num )
  {
  InspiralTemplate **bankArray = NULL;
  InspiralTemplate *bankFirst = NULL;
  UINT4 i = 0;
  bankFirst = bankHead;
  bankArray = (InspiralTemplate **) LALCalloc(num, sizeof(InspiralTemplate *));

  for (i = 0; (i < num); bankHead = bankHead->next, i++)
  {
    bankArray[i] = bankHead; /* populate pointer array */
  }

  qsort(bankArray, num, sizeof(InspiralTemplate *), compareTemplateByEta);

  bankFirst = bankHead = bankArray[0];
  /* repopulate linked list */
  for (i=1; i < num; i++)
  {
    bankHead = bankHead->next = bankArray[i];
  }
  bankHead->next = NULL;
  LALFree(bankArray);
  return bankFirst;
}


static int  breakUpRegions(InspiralTemplate *bankHead)
  {
  gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
  gsl_rng_set (r, 1);
  while (bankHead)
    {
    bankHead->level = (UINT4) gsl_rng_uniform_int (r, 1000000000);
    bankHead = bankHead->next;
    }
  gsl_rng_free(r);
  return 0;
  }


static InspiralTemplate *
XLALFindChirpSortTemplatesByLevel( InspiralTemplate *bankHead, UINT4 num)
  {
  InspiralTemplate **bankArray = NULL;
  InspiralTemplate *bankFirst = NULL;
  UINT4 i = 0;
  bankFirst = bankHead;
  bankArray = (InspiralTemplate **) LALCalloc(num, sizeof(InspiralTemplate *));

  i = 0;
/*  for (i = 0; (i < num); bankHead = bankHead->next, i++)*/
  while(bankHead)
  {
    bankArray[i] = bankHead; /* populate pointer array */
    bankHead=bankHead->next;
    i++;
  }

  qsort(bankArray, num, sizeof(InspiralTemplate *), compareTemplateByLevel);

  bankFirst = bankHead = bankArray[0];
  /* repopulate linked list */
  for (i=1; i < num; i++)
  {
    bankHead = bankHead->next = bankArray[i];
  }
  bankHead->next = NULL;
  LALFree(bankArray);
  return bankFirst;
}


static int compareTemplateByLevel (const void * a, const void * b)
{
  REAL4 mVal1 = (*(InspiralTemplate * const *)a)->level;
  REAL4 mVal2 = (*(InspiralTemplate * const *)b)->level;

  if ( mVal1 > mVal2 ) return 1;
  if ( mVal1 == mVal2 ) return 0;
  if ( mVal1 < mVal2 ) return -1;
  return 0;
}


static int compareTemplateByChirpMass (const void * a, const void * b)
{
  REAL4 mVal1 = (*(InspiralTemplate * const *)a)->chirpMass;
  REAL4 mVal2 = (*(InspiralTemplate * const *)b)->chirpMass;

  if ( mVal1 > mVal2 ) return 1;
  if ( mVal1 == mVal2 ) return 0;
  if ( mVal1 < mVal2 ) return -1;
  return 0;
}

static int compareTemplateByEta (const void * a, const void * b)
{
  REAL4 m1Val1 = (*(InspiralTemplate * const *)a)->mass1;
  REAL4 m2Val1 = (*(InspiralTemplate * const *)a)->mass2;
  REAL4 m1Val2 = (*(InspiralTemplate * const *)b)->mass1;
  REAL4 m2Val2 = (*(InspiralTemplate * const *)b)->mass2;

  REAL4 eta1 = m1Val1 * m2Val1 / (m1Val1 + m2Val1);
  REAL4 eta2 = m1Val2 * m2Val2 / (m1Val2 + m2Val2);

  if ( eta1 > eta2 ) return 1;
  if ( eta1 == eta2 ) return 0;
  if ( eta1 < eta2 ) return -1;
  return 0;
}


static int compareTemplateByMass (const void * a, const void * b)
{
  REAL4 mVal1 = (*(InspiralTemplate * const *)a)->mass1 + (*(InspiralTemplate * const *)a)->mass2;
  REAL4 mVal2 = (*(InspiralTemplate * const *)b)->mass1 + (*(InspiralTemplate * const *)b)->mass2;

  if ( mVal1 > mVal2 ) return 1;
  if ( mVal1 == mVal2 ) return 0;
  if ( mVal1 < mVal2 ) return -1;
  return 0;
}


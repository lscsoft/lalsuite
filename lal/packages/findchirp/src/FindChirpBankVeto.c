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
  bankVetoData->timeshift= NULL;
}

void XLALDestroyBankVetoData(FindChirpBankVetoData *bvdata)
{
  if (bvdata->acorr) XLALDestroyREAL4Vector(bvdata->acorr);
  if (bvdata->workspace) XLALDestroyCOMPLEX8Vector(bvdata->workspace);
  if (bvdata->acorrMat) XLALDestroyREAL4Vector(bvdata->acorrMat);
  if (bvdata->revplan) XLALDestroyREAL4FFTPlan(bvdata->revplan);
  if (bvdata->timeshift) XLALDestroyREAL4Vector(bvdata->timeshift);
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
                    REAL4Vector *ampVec,
                    UINT4 subBankSize,
                    REAL4 dynRange,
                    REAL4 fLow,
                    REAL4 deltaF,
		    REAL4 deltaT)
{
    /* 
     * Temporary variables that hold the real and imaginary parts
     * of the two templates, which will be multiplied together
     * and added to the cross-correlation integral.  One should remember
     * that the ROW variable gets the complex conjugate in the integral.
     */
    REAL8 tmpltRowReal = 0;  
    REAL8 tmpltRowImag = 0;    
    REAL8 tmpltColReal = 0;    
    REAL8 tmpltColImag = 0;  
    REAL8 tmpltRowRealShifted = 0;  
    REAL8 tmpltRowImagShifted = 0;    
    REAL8 tmpltColRealShifted = 0;    
    REAL8 tmpltColImagShifted = 0;  

    /* 
     * The cumulative real and imaginary parts of the cross-correlation
     * integral int( A*(f)B(f)/Sn(f) )df 
     */
    REAL8 crossCorrReal; 
    REAL8 crossCorrImag;   

    /* Indices used to iterate over the templates in the subbank */
    UINT4 row; 
    UINT4 col;   

    /* Index used to integrate over frequencies */
    UINT4 sample;           
    
    /* Indices which dictate the start and end of integration */
    UINT4 startIndex = floor( fLow / deltaF ); 
    UINT4 sampleMaxTmpltRow;  // Maximum index available in template A
    UINT4 sampleMaxTmpltCol;  // Maximum index available in template B
    UINT4 templateLength = bankVetoData->fcInputArray[0]->fcTmplt->data->length;

   /* Absolute value of calibration response function*/
    REAL8 sqResp; 
    REAL8 spectralDensity;
 
   /* Parameters used to time shift templates */
     REAL8 phaseOffsetRow;
     REAL8 phaseOffsetCol;

    /* PSD bucket frequency used for chirp time offset */
    double fBucket = 150.0;   

 
    /* FIXME this should be a command line argument */
     /*200 points of autocorrelation function stored */
    bankVetoData->acorrMatSize = 200;

    /* Allocate memory for workspace and the autocorrelation if necessary */
    if ( !bankVetoData->acorr ) 
      bankVetoData->acorr = XLALCreateREAL4Vector((templateLength-1) * 2);
    if ( !bankVetoData->workspace) 
      bankVetoData->workspace = XLALCreateCOMPLEX8Vector(templateLength);
    if ( !bankVetoData->acorrMat) 
      bankVetoData->acorrMat = XLALCreateREAL4Vector(bankVetoData->acorrMatSize * bankVetoData->length);
    if ( !bankVetoData->revplan) 
      bankVetoData->revplan = XLALCreateReverseREAL4FFTPlan((templateLength-1) * 2 , 0);

    /* Decide how to time shift each template if necessary */
    if ( !bankVetoData->timeshift) 
    {
	bankVetoData->timeshift = XLALCreateREAL4Vector(bankVetoData->length);
	
	for (row = 0; row < subBankSize; row++ )
	{
	    bankVetoData->timeshift->data[row] = 0;

	    /*(REAL4) ( chirp_time_between_f1_and_f2(bankVetoData->fcInputArray[row]->fcTmplt->tmplt.mass1,
						     bankVetoData->fcInputArray[row]->fcTmplt->tmplt.mass2,
						     fBucket,
						     bankVetoData->fcInputArray[row]->fcTmplt->tmplt.fFinal, 
						     7) );*/
	
	}		
    }

    /* Compute the cross-correlation between template row and template col */
    for ( row = 0; row < bankVetoData->length; row++ )
    {
	for ( col = row; col < bankVetoData->length; col++ )
	{          
	    /* zero out the workspace vector that will contain A*(f) A(f) */
	    if ( row == col  )
	      memset(bankVetoData->workspace->data, 0, bankVetoData->workspace->length * sizeof(COMPLEX8));
	   
	    /* initialize the cross-correlation integral */
	    crossCorrReal = crossCorrImag = 0;
	    	    
	    if ( (row < subBankSize) && (col < subBankSize) )
	    {
		/* find the ending sample point of the templates */
		sampleMaxTmpltRow = 
		  bankVetoData->fcInputArray[row]->fcTmplt->tmplt.fFinal / deltaF;
		sampleMaxTmpltCol = 
		  bankVetoData->fcInputArray[col]->fcTmplt->tmplt.fFinal / deltaF;

		/* Begin computation of the cross-correlation. */	  
		for ( sample = startIndex; sample < templateLength; sample++ )
		{
		    //if ( (sample > sampleMaxTmpltRow) || (sample > sampleMaxTmpltCol) ) break;

		    sqResp = ( bankVetoData->resp->data[sample].re *
			       bankVetoData->resp->data[sample].re +
			       bankVetoData->resp->data[sample].im *
			       bankVetoData->resp->data[sample].im ) * dynRange * dynRange;

		    spectralDensity = bankVetoData->spec->data[sample]*sqResp;
		    
		    if ( spectralDensity != 0 )
		    {
			tmpltRowReal = (bankVetoData->fcInputArray[row]->fcTmplt->data->data[sample].re)*ampVec->data[sample];
			tmpltRowImag = (bankVetoData->fcInputArray[row]->fcTmplt->data->data[sample].im)*ampVec->data[sample];

			tmpltColReal = (bankVetoData->fcInputArray[col]->fcTmplt->data->data[sample].re)*ampVec->data[sample];			
			tmpltColImag = (bankVetoData->fcInputArray[col]->fcTmplt->data->data[sample].im)*ampVec->data[sample];

			/* 
			 * Time shift tmpltRow forward by timeshift[row] and 
			 * tmpltCol forward by timeshift[col].
			 */
			phaseOffsetRow = 2.0*LAL_PI*((REAL8)deltaF)*((REAL8)sample)*((REAL8)bankVetoData->timeshift->data[row]);
			phaseOffsetCol = 2.0*LAL_PI*((REAL8)deltaF)*((REAL8)sample)*((REAL8)bankVetoData->timeshift->data[col]);			

			tmpltRowRealShifted = tmpltRowReal*cos(phaseOffsetRow) - tmpltRowImag*sin(phaseOffsetRow);
			tmpltRowImagShifted = tmpltRowImag*cos(phaseOffsetRow) + tmpltRowReal*sin(phaseOffsetRow);
			
			tmpltColRealShifted = tmpltColReal*cos(phaseOffsetCol) - tmpltColImag*sin(phaseOffsetCol);
			tmpltColImagShifted = tmpltColImag*cos(phaseOffsetCol) + tmpltColReal*sin(phaseOffsetCol);

			
			/*  
			 *  tmpltRow* x tmpltCol/Sn = (tmpltRowReal-i*tmpltRowImag) x (tmpltColReal+i*tmpltColImag) / Sn =
			 *     (tmpltRowReal*tmpltColReal + tmpltRowImag*tmpltColImag) / Sn +
			 *   i*(tmpltRowReal*tmpltColImag - tmpltRowImag*tmpltColReal) / Sn
			 */
			crossCorrReal += (REAL8) (tmpltRowRealShifted*tmpltColRealShifted + tmpltRowImagShifted*tmpltColImagShifted)/spectralDensity;
			crossCorrImag += (REAL8) (tmpltRowRealShifted*tmpltColImagShifted - tmpltRowImagShifted*tmpltColRealShifted)/spectralDensity;			
		    }

		    if ( row == col )
		    {	      
			/* set the workspace to contain the frequency series A*(f) x A(f) / Sn(f) */
			bankVetoData->workspace->data[sample].re = (REAL4)  (tmpltRowReal*tmpltColReal + tmpltRowImag*tmpltColImag)/spectralDensity;
			/* imaginary part must be 0, so we force it to be here */
			bankVetoData->workspace->data[sample].im = 0.0;
		    }		    
		}//end cross-correlation integration
	
			       
		/* if the templates are the same then do an autocorrelation function */
		if ( row == col )
		{ 
		    /* 
		     * The autocorrelation is integral( A*(f) x A(f) exp(2*pi*i*f*t) )df, i.e.,
		     * the reverse FFT of the workspace vector. 
		     */
		    XLALREAL4ReverseFFT( bankVetoData->acorr, bankVetoData->workspace, bankVetoData->revplan );
		    /* since the normalized first sample of the auto correlation is always 1, lets skip it */
		    for (sample=1; sample < bankVetoData->acorrMatSize + 1; sample++)
		    {
			/* Save the last acorrMatSize samples */
			bankVetoData->acorrMat->data[row*bankVetoData->acorrMatSize + sample-1] = bankVetoData->acorr->data[sample] / bankVetoData->acorr->data[0];
		    }

		}//end autocorrelation computation
	       
	    }	    
	    /* 
	     * Fill in the cross-correlation matrix with computed values.
	     * The matrix is hermitian so the imaginary flips sign. 
	     */
	    bankVetoData->ccMat->data[row*bankVetoData->length + col].re = (REAL4) crossCorrReal;
	    bankVetoData->ccMat->data[row*bankVetoData->length + col].im = (REAL4) crossCorrImag;

	    bankVetoData->ccMat->data[col*bankVetoData->length + row].re = (REAL4) crossCorrReal;
	    bankVetoData->ccMat->data[col*bankVetoData->length + row].im = (REAL4) -crossCorrImag;
	    
	    /* Store normalization factors normMat(i) = sqrt(<hi|hi>) (see below) */
	    if ( row == col ) 
	    {
		bankVetoData->normMat->data[row] = (REAL4) sqrt(crossCorrReal);
	    }
	}
    } /* end cross-correlation matrix computation */
    

    /* 
     * Normalize the templates in the template bank such that <hi|hi> = 1 for all hi 
     * in the template bank.  This amounts to dividing each hi by sqrt(<hi|hi>), or
     * equivalently (and as implemented here) dividing ccMat(i,j) = <hi|hj> by 
     * normMat(i)*normMat(j) = sqrt(<hi|hi>)sqrt(<hj|hj>).
     */
    for ( row = 0; row < subBankSize; row++ ) 
    {
	for ( col = 0; col < subBankSize; col++ ) 
	{
	    bankVetoData->ccMat->data[row*bankVetoData->length + col].re /= 
	      (bankVetoData->normMat->data[row] *bankVetoData->normMat->data[col]);

	    bankVetoData->ccMat->data[row*bankVetoData->length + col].im /=
	      (bankVetoData->normMat->data[row] *bankVetoData->normMat->data[col]);
	}
    } /* end norm loop */

    return;

}


REAL4
XLALComputeBankVeto( FindChirpBankVetoData *bankVetoData,
                     UINT4 row,
                     UINT4 snrIndexRow,
		     REAL4 deltaT,
                     UINT4 *dof)
{

    /* SNR (and derived quantities) for reference (input) template */
    REAL8 rowSNR_real = 0;
    REAL8 rowSNR_imag = 0;
    REAL8 rowSNR_mag = 0;
    REAL8 rowSNR_phi = 0;
    
    /* SNR (and derived quantities) for comparison templates */
    REAL8 colSNR_real = 0;
    REAL8 colSNR_imag = 0;
    
    /* index to loop over comparison templates */
    UINT8 col;

    /* Time index at which to get the SNR of the comparison template.  This
     * will be related to the SNR index of the input template by the time shift
     * introduced in the templates in the cross-correlation computation. */
    UINT8 snrIndexCol;

    /* The time shift between the reference and comparison template */
    REAL4 deltaTimeShift;

    /* temporary storage for cross-corr data */
    COMPLEX16 crossCorr;
    REAL8 crossCorr_mag;

    /* complex cross correlation */
    COMPLEX16 chisq;
    
    /* normalization factor to make chi-sq distributed */
    REAL4 bankNorm;

    /* Output of function */
    REAL4 chisq_mag = 0;


    /* Handle trivial case */
    if (bankVetoData->length == 1) return 0;
   
    /* Lookup SNR for template A at time given by snrIndexRow */
    rowSNR_real = bankVetoData->qVecArray[row]->data[snrIndexRow].re
      * sqrt(bankVetoData->fcInputArray[row]->fcTmplt->norm);
    rowSNR_imag = bankVetoData->qVecArray[row]->data[snrIndexRow].im
      * sqrt(bankVetoData->fcInputArray[row]->fcTmplt->norm);
    
    rowSNR_mag = sqrt(rowSNR_real * rowSNR_real + rowSNR_imag * rowSNR_imag);
    rowSNR_phi = atan2(rowSNR_imag, rowSNR_real);

    chisq_mag = 0;
    *dof = 0;

    /* 
     * For each comparison template, compare its SNR to the SNR it would have if the signal were
     * precisely equal to the reference template. Deviations from this expected SNR will be chi-sq
     * distributed in Guassian noise. 
     */ 
    for (col = 0; col < bankVetoData->length; col++)
    {
	if (row == col) continue;
	if (bankVetoData->ccMat->data[row*bankVetoData->length + col].re == 0 && bankVetoData->ccMat->data[row*bankVetoData->length + col].im == 0) continue;
	
	/* Phase shift the cross correlation of template B with template B by rowSNR_phi */ 
	crossCorr.re = bankVetoData->ccMat->data[row*bankVetoData->length + col].re * cos(rowSNR_phi)
	                 - bankVetoData->ccMat->data[row*bankVetoData->length + col].im * sin(rowSNR_phi);
	crossCorr.im = bankVetoData->ccMat->data[row*bankVetoData->length + col].im * cos(rowSNR_phi) 
	                 + bankVetoData->ccMat->data[row*bankVetoData->length + col].re * sin(rowSNR_phi);

	crossCorr_mag = crossCorr.re * crossCorr.re + crossCorr.im * crossCorr.im;
	bankNorm = 2.0 - crossCorr_mag;

	/* FIXME: warning -- this could go off the edge of the SNR time series */
	/* Look at the col SNR deltaTimeShift/deltaT samples in the past.  The choice of sign
	 * depends on the sign convention choice made in computing ccmat. */
	deltaTimeShift = bankVetoData->timeshift->data[row] - bankVetoData->timeshift->data[col];
	snrIndexCol = snrIndexRow - (UINT4) round(deltaTimeShift/deltaT);

	colSNR_real = bankVetoData->qVecArray[col]->data[snrIndexCol].re
	  * sqrt(bankVetoData->fcInputArray[col]->fcTmplt->norm);
	colSNR_imag = bankVetoData->qVecArray[col]->data[snrIndexCol].im
	  * sqrt(bankVetoData->fcInputArray[col]->fcTmplt->norm);
      
	chisq.re = colSNR_real - crossCorr.re * (rowSNR_mag);
	chisq.im = colSNR_imag - crossCorr.im * (rowSNR_mag);
	chisq_mag += (chisq.re*chisq.re + chisq.im*chisq.im) / bankNorm;
	(*dof)++;

    }

    /*FIXME Normalization now not necessarily chisquare distributed */
    return chisq_mag;
}


InspiralTemplate *
XLALFindChirpSortTemplates( InspiralTemplate *bankHead, UINT4 num)
  {
  /*return XLALFindChirpSortTemplatesByEta(bankHead, num);*/
  /*bankHead = XLALFindChirpSortTemplatesByChirpMass(bankHead,num);*/
  breakUpRegions(bankHead);
  return XLALFindChirpSortTemplatesByLevel(bankHead,num);
  //return XLALFindChirpSortTemplatesByChirpMass(bankHead,num);
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


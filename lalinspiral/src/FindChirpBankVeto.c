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
 *-----------------------------------------------------------------------
 */

/**

\author Brown, D. A., and Hanna, C.
\file

*/

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <time.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/* Some static function prototypes */
static int compareTemplateByLevel (const void * a, const void * b);
/*static int compareTemplateByChirpMass (const void * a, const void * b);*/
/*static int compareTemplateByMass (const void * a, const void * b);*/
static int compareTemplateByFfinal (const void * a, const void * b);
static int  breakUpRegions(InspiralTemplate *bankHead);
static int  breakUpSubBanks(InspiralTemplate *bankHead, UINT4 num, UINT4 size);
static InspiralTemplate * XLALFindChirpSortTemplatesByLevel( InspiralTemplate *bankHead, UINT4 num);
/*static InspiralTemplate * XLALFindChirpSortTemplatesByChirpMass( InspiralTemplate *bankHead, UINT4 num );*/
static InspiralTemplate * XLALFindChirpSortTemplatesByFfinal( InspiralTemplate *bankHead, UINT4 num );
/*static InspiralTemplate * XLALFindChirpSortTemplatesByMass( InspiralTemplate *bankHead, UINT4 num );*/

/* helper functions */

static COMPLEX8 rotate(REAL4 a, REAL4 phi)
	{
	COMPLEX8 out;
        out = crectf( a * cos(phi), a * sin(phi) );
	return out;
	}
#if 0
static COMPLEX16 rotate_complex16(COMPLEX16 a, REAL8 phi)
	{
	COMPLEX16 out;
        out.re = a.re * cos(phi) - a.im * sin(phi);
	out.im = a.re * sin(phi) + a.im * cos(phi);
	return out;
	}
#endif

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

{
  UINT4                 i;
  UINT4                 numSubBanks = 0;
  UINT4                 subBankRemainder = 0;
  UINT4                *bankSizes = NULL;
  FindChirpSubBank     *subBankHead = NULL;
  FindChirpSubBank     *thisSubBank = NULL;
  InspiralTemplate     *thisTmplt = NULL;
  InspiralTemplate     *nextTmplt = NULL;

  /* if there are no templates return NULL */
  if ( ! bankSize )
    return NULL;

  numSubBanks = bankSize / subBankSize;
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

  /* create an array of subbank sizes with the minimum size */
  bankSizes = (UINT4 *) LALCalloc( numSubBanks, sizeof(UINT4) );
  for ( i = 0; i < numSubBanks; ++i )
  {
	bankSizes[i] = subBankSize;
  }

  /* disperse the remainder through the subbanks */
  subBankRemainder = bankSize % subBankSize;
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

  /* find the largest subbank */
  *maxSubBankSize = 0;
  for ( i = 0; i < numSubBanks; i++)
  {
      if ( *maxSubBankSize < bankSizes[i] )
      {
	  *maxSubBankSize = bankSizes[i];
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
    FindChirpFilterInput       UNUSED *input,
    FindChirpFilterParams      UNUSED *params,
    COMPLEX8                   *q,
    UINT4 			i,
    UINT4                       snrIX,
    UINT4                      *dof,
    REAL4                      norm
)

{
  UINT4 k;
  REAL4 chisqnorm, tmp, C, chisq, angle, snri, snrk;
  UINT4 numsamps = bankVetoData->acorrMatSize;

  /* test isn't being done */
  if (!bankVetoData->acorrMat) return 0.0;
  chisq = 0;
  chisqnorm = 0;

  angle = atan2(cimagf(q[snrIX]), crealf(q[snrIX]));
  snri = crealf(q[snrIX]) * cos(angle) + cimagf(q[snrIX]) * sin(angle);
  /* have the dof vary with snr^1/2 up until snr = 1000 */
  //if (snri < 1000) numsamps = (UINT4) ceil((double) numsamps * sqrt(snri) / 31.622776601683793);
  *dof = 0;

  for (k = 0; k < numsamps; k+=bankVetoData->autochisqStride)
  {
    /* we skipped saving the first sample of the autocorrelation which is why */
    /* we need the -1 in the snrk here */
    C = bankVetoData->acorrMat->data[i * bankVetoData->acorrMatSize + k];
    chisqnorm = 1.0 - C * C;
    snrk = crealf(q[snrIX-k-1]) * cos(angle) + cimagf(q[snrIX-k-1]) * sin(angle);
    (*dof) += 1;
    tmp = snrk - C * snri;
    chisq += tmp * tmp * norm / chisqnorm;
    /* Other side */
    if (bankVetoData->two_sided_auto_chisq)
      {
      snrk = crealf(q[snrIX+k+1]) * cos(angle) + cimagf(q[snrIX+k+1]) * sin(angle);
      (*dof) += 1;
      tmp = snrk - C * snri;
      chisq += tmp * tmp * norm / chisqnorm;
      }
  }

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
    COMPLEX16 tmp = 0.0;
    /* FIXME someday this could track time offset phase */
    /*REAL8 tphase = 0.0; */

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

    /* Allocate memory for workspace and the autocorrelation if necessary */
    if ( !bankVetoData->acorr )
      bankVetoData->acorr = XLALCreateREAL4Vector((templateLength-1) * 2);
    if ( !bankVetoData->workspace)
      bankVetoData->workspace = XLALCreateCOMPLEX8Vector(templateLength);
    if ( !bankVetoData->acorrMat)
      bankVetoData->acorrMat = XLALCreateREAL4Vector(bankVetoData->acorrMatSize * bankVetoData->length);
    if ( !bankVetoData->revplan)
      bankVetoData->revplan = XLALCreateReverseREAL4FFTPlan((templateLength-1) * 2 , 1);

    /* FIXME no time shifts are done yet, we have to linearly shift the phase during
     * the inner product too
     */
    /* Decide how to time shift each template if necessary */
    if (!bankVetoData->timeshift)
    {
	bankVetoData->timeshift = XLALCreateREAL4Vector(bankVetoData->length);

	for (row = 0; row < subBankSize; row++ )
	{
	/* FIXME someday maybe this could use time offsets */
	if (bankVetoData->time_freq_bank_veto) bankVetoData->timeshift->data[row] = (REAL8) 0.0 * deltaT;
	else bankVetoData->timeshift->data[row] = 0.0;
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
		    if ( (sample > sampleMaxTmpltRow) || (sample > sampleMaxTmpltCol) ) break;

		    sqResp = ( crealf(bankVetoData->resp->data[sample]) *
			       crealf(bankVetoData->resp->data[sample]) +
			       cimagf(bankVetoData->resp->data[sample]) *
			       cimagf(bankVetoData->resp->data[sample]) ) * dynRange * dynRange;

		    spectralDensity = bankVetoData->spec->data[sample]*sqResp;
		    /*FIXME time shifting stuff */
		    /*tphase = (REAL8) (col - row) * 2.0 * LAL_PI * bankVetoData->timeshift->data[col] * sample * deltaF;*/

		    if ( spectralDensity != 0 )
		    {
			tmpltRowReal = crealf(bankVetoData->fcInputArray[row]->fcTmplt->data->data[sample])*ampVec->data[sample];
			tmpltRowImag = cimagf(bankVetoData->fcInputArray[row]->fcTmplt->data->data[sample])*ampVec->data[sample];

			tmpltColReal = crealf(bankVetoData->fcInputArray[col]->fcTmplt->data->data[sample])*ampVec->data[sample];
			tmpltColImag = cimagf(bankVetoData->fcInputArray[col]->fcTmplt->data->data[sample])*ampVec->data[sample];

			tmp = crect(
                          (REAL8) (tmpltRowReal*tmpltColReal + tmpltRowImag*tmpltColImag)/spectralDensity,
                          (REAL8) (tmpltRowReal*tmpltColImag - tmpltRowImag*tmpltColReal)/spectralDensity
                          );
			/* FIXME Apply time shift */
			/*tmp = rotate_complex16(tmp, tphase);*/

			crossCorrReal += (REAL8) creal(tmp);
			crossCorrImag += (REAL8) cimag(tmp);
		    }

		    if ( row == col )
		    {
			/* set the workspace to contain the frequency series A*(f) x A(f) / Sn(f) */
			/* imaginary part must be 0, so we force it to be here */
                        bankVetoData->workspace->data[sample] = (REAL4)((tmpltRowReal*tmpltColReal + tmpltRowImag*tmpltColImag)/spectralDensity);
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

	    bankVetoData->ccMat->data[row*bankVetoData->length + col] = crectf( (REAL4) crossCorrReal, (REAL4) crossCorrImag );

	    bankVetoData->ccMat->data[col*bankVetoData->length + row] = crectf( (REAL4) crossCorrReal, (REAL4) -crossCorrImag );

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
	    bankVetoData->ccMat->data[row*bankVetoData->length + col] /=
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
    COMPLEX8 rowSNR = 0.0;

    /* SNR (and derived quantities) for comparison templates */
    COMPLEX8 colSNR = 0.0;
    COMPLEX8 expSNR = 0.0;

    /* index to loop over comparison templates */
    UINT4 col;

    /* Time index at which to get the SNR of the comparison template.  This
     * will be related to the SNR index of the input template by the time shift
     * introduced in the templates in the cross-correlation computation. */
    UINT4 snrIndexCol;

    /* The time shift between the reference and comparison template */
    REAL4 deltaTimeShift;

    /* temporary storage for cross-corr data */
    COMPLEX8 crossCorr;

    /* complex cross correlation */
    COMPLEX8 chisq;

    /* Output of function */
    REAL4 chisq_mag = 0;

    /* Handle trivial case */
    if (bankVetoData->length == 1) return 0;

    /* Lookup SNR for template A at time given by snrIndexRow */
    rowSNR = bankVetoData->qVecArray[row]->data[snrIndexRow] *
             sqrt(bankVetoData->fcInputArray[row]->fcTmplt->norm);

    *dof = 0;

    /*
     * For each comparison template, compare its SNR to the SNR it would have if the signal were
     * precisely equal to the reference template. Deviations from this expected SNR will be chi-sq
     * distributed in Guassian noise.
     */
    for (col = 0; col < bankVetoData->length; col++)
    {

        /* don't do the diagonal case as it is trivial */
	if (row == col) continue;

	crossCorr = bankVetoData->ccMat->data[row*bankVetoData->length + col];

        /* having the crossCorr == 0 means that we are outside the valid range for this subbank */
	if (crossCorr == 0) continue;

	/* FIXME: warning -- this could go off the edge of the SNR time series */
	/* Look at the col SNR deltaTimeShift/deltaT samples in the past.  The choice of sign
	 * depends on the sign convention choice made in computing ccmat. */
	/* FIXME there are not any time shifts yet, so this is a no-op for now */
	/* i.e. SNRIndexCol = SNRIndexRow */
	deltaTimeShift = ((REAL4) col - (REAL4) row) * bankVetoData->timeshift->data[col];
	if (deltaTimeShift > 0) snrIndexCol = snrIndexRow - (UINT4) round(deltaTimeShift/deltaT);
	else snrIndexCol = snrIndexRow - (UINT4) round(deltaTimeShift/deltaT);
	//fprintf(stderr, "snrIndexRow %d snrIndexCol %d\n", snrIndexRow, snrIndexCol);

	/* Get the "column" snr in the same way as the row */
	colSNR = bankVetoData->qVecArray[col]->data[snrIndexCol] *
		                      sqrt(bankVetoData->fcInputArray[col]->fcTmplt->norm);

	expSNR = rotate(cabsf(rowSNR)*cabsf(crossCorr), cargf(rowSNR) - cargf(crossCorr));

	/* Subtract the expected SNR from the measured one */
	chisq = colSNR - expSNR;

	/* Square the result */
	chisq_mag += cabsf(chisq) * cabsf(chisq);

	//fprintf(stderr, "%.2f %.2f expSNR %.2f colSNR %.2f\n", sqrt(rowSNR.re*rowSNR.re + rowSNR.im*rowSNR.im), chisq_mag, cabsf(expSNR), cabsf(colSNR));

	(*dof)+=2;
    }
    return chisq_mag;
}


InspiralTemplate *
XLALFindChirpSortTemplates( InspiralTemplate *bankHead, FindChirpBankVetoData *bvdata, UINT4 num, UINT4 max_subbank_size)
  {
  /* Either do something clever */
  if (bvdata->time_freq_bank_veto)
  {
    bankHead = XLALFindChirpSortTemplatesByFfinal(bankHead, num);
    breakUpSubBanks(bankHead, num, max_subbank_size);
    return XLALFindChirpSortTemplatesByLevel(bankHead,num);
  }
  /* Or just randomly choose templates */
  breakUpRegions(bankHead);
  return XLALFindChirpSortTemplatesByLevel(bankHead,num);
  }


#if 0
InspiralTemplate *
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
#endif

#if 0
InspiralTemplate *
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
#endif

InspiralTemplate *
XLALFindChirpSortTemplatesByFfinal( InspiralTemplate *bankHead, UINT4 num )
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

  qsort(bankArray, num, sizeof(InspiralTemplate *), compareTemplateByFfinal);

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

static int  breakUpSubBanks(InspiralTemplate *bankHead, UINT4 num, UINT4 size)
  {
  UINT4 cnt = 0;
  while (bankHead)
    {
    if (num / size > 0) bankHead->level = (UINT4) cnt % (num / size);
    else bankHead->level = 0;
    cnt++;
    bankHead = bankHead->next;
    }
  return 0;
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

#if 0
static int compareTemplateByChirpMass (const void * a, const void * b)
{
  REAL4 mVal1 = (*(InspiralTemplate * const *)a)->chirpMass;
  REAL4 mVal2 = (*(InspiralTemplate * const *)b)->chirpMass;

  if ( mVal1 > mVal2 ) return 1;
  if ( mVal1 == mVal2 ) return 0;
  if ( mVal1 < mVal2 ) return -1;
  return 0;
}
#endif

static REAL4 imr_ring(double m1, double m2){
        double chi = 0.0;
	double totalMass = m1+m2;
        double eta = m1*m2/pow(totalMass,2.);
        double piM = totalMass*LAL_PI*LAL_MTSUN_SI;
        double fRing = (1. - 0.63*pow(1.-chi,0.3))/2. +
                1.4690e-01*eta + -1.2281e-01*eta*chi + -2.6091e-02*eta*pow(chi,2.) +
                -2.4900e-02*pow(eta,2.) + 1.7013e-01*pow(eta,2.)*chi +
                2.3252e+00*pow(eta,3.);
        fRing /= piM;
        return (REAL4) fRing;
        }

static int compareTemplateByFfinal (const void * a, const void * b)
{
  REAL4 m1Val1 = (*(InspiralTemplate * const *)a)->mass1;
  REAL4 m2Val1 = (*(InspiralTemplate * const *)a)->mass2;
  REAL4 m1Val2 = (*(InspiralTemplate * const *)b)->mass1;
  REAL4 m2Val2 = (*(InspiralTemplate * const *)b)->mass2;

  REAL4 ffinal1 = imr_ring(m1Val1, m2Val1);
  REAL4 ffinal2 = imr_ring(m1Val2, m2Val2);

  if ( ffinal1 > ffinal2 ) return 1;
  if ( ffinal1 == ffinal2 ) return 0;
  if ( ffinal1 < ffinal2 ) return -1;
  return 0;
}

#if 0
static int compareTemplateByMass (const void * a, const void * b)
{
  REAL4 mVal1 = (*(InspiralTemplate * const *)a)->mass1 + (*(InspiralTemplate * const *)a)->mass2;
  REAL4 mVal2 = (*(InspiralTemplate * const *)b)->mass1 + (*(InspiralTemplate * const *)b)->mass2;

  if ( mVal1 > mVal2 ) return 1;
  if ( mVal1 == mVal2 ) return 0;
  if ( mVal1 < mVal2 ) return -1;
  return 0;
}
#endif

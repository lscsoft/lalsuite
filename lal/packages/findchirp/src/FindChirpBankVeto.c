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
}

void XLALDestroyBankVetoData(FindChirpBankVetoData *bvdata)
{
  if (bvdata->acorr) XLALDestroyREAL4Vector(bvdata->acorr);
  if (bvdata->workspace) XLALDestroyCOMPLEX8Vector(bvdata->workspace);
  if (bvdata->acorrMat) XLALDestroyREAL4Vector(bvdata->acorrMat);
  if (bvdata->revplan) XLALDestroyREAL4FFTPlan(bvdata->revplan);
  return; 
}


/* <lalVerbatim file="FindChirpBankVetoCP"> */
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

  subBankRemainder = bankSize % subBankSize;

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
  *maxSubBankSize = 0;
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

    /* store the size of the biggest bank */
    if ( bankSizes[i] > *maxSubBankSize )
    {
      *maxSubBankSize = bankSizes[i];
    }
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
                    REAL4 deltaF,
                    REAL4 deltaT )

{
  UINT4 i,j,k,iMax,jMax,correctFlag;
  UINT4 iSize = bankVetoData->length;
  UINT4 tmpLen = bankVetoData->fcInputArray[0]->fcTmplt->data->length;
  REAL8 ABr, ABi, Br, Bi, sqResp, normfac;
  UINT4 stIX = floor( fLow / deltaF );
  double deltaChirp;

  /* FIXME this should be a command line argument */
  bankVetoData->acorrMatSize = 200; /*200 points of autocorrelation function stored */

  /* if there isn't alread memory allocated for workspace and the autocorrelation 
   * then do it
   */
  if ( !bankVetoData->acorr ) bankVetoData->acorr = XLALCreateREAL4Vector((tmpLen-1) * 2);
  if ( !bankVetoData->workspace) bankVetoData->workspace = XLALCreateCOMPLEX8Vector(tmpLen);
  if ( !bankVetoData->acorrMat) bankVetoData->acorrMat = XLALCreateREAL4Vector(bankVetoData->acorrMatSize * iSize);
  if ( !bankVetoData->revplan) bankVetoData->revplan = XLALCreateReverseREAL4FFTPlan((tmpLen-1) * 2 , 0);
  if ( !bankVetoData->tchirp) bankVetoData->tchirp = (double *) calloc(iSize*sizeof(double));

  /* deltaT is unused in this function */
  UNUSED(deltaT);

  for ( i = 0; i < iSize; i++ )
    {
      /* FIXME find actual lal function */
      bankVetoData->tchirp[i] = chirp_time_between_f1_and_f2(bankVetoData->fcInputArray[i]->fcTmplt->tmplt.mass1, 
							     bankVetoData->fcInputArray[i]->fcTmplt->tmplt.mass2,
							     fLow,
							     bankVetoData->fcInputArray[i]->fcTmplt->tmplt.fFinal,
							     7);
    }


  for ( i = 0; i < iSize; i++ )
  {
    correctFlag = 1;
    for ( j = i; j < iSize; j++ )
    {
      /* FIXME find an XLAL Function to do this?*/
      if (i == j) memset(bankVetoData->workspace->data, 0, bankVetoData->workspace->length * sizeof(COMPLEX8));
      if ( (i >= vetoBank->subBankSize) || (j >= vetoBank->subBankSize) )
      {
        bankVetoData->ccMat->data[i*iSize + j].re = 0.0;
        bankVetoData->ccMat->data[j*iSize + i].re = 0.0;
        bankVetoData->ccMat->data[i*iSize + j].im = 0.0;
        bankVetoData->ccMat->data[j*iSize + i].im = 0.0;
        bankVetoData->normMat->data[i*iSize + j] = 0.0;
        bankVetoData->normMat->data[j*iSize + i] = 0.0;
      }
      else
      {
        /* find the ending sample point of the templates */
        iMax = bankVetoData->fcInputArray[i]->fcTmplt->tmplt.fFinal / deltaF;
        jMax = bankVetoData->fcInputArray[j]->fcTmplt->tmplt.fFinal / deltaF;
        ABr = ABi = 0;
        for ( k = stIX; k < tmpLen-1; k++ )
        {
          if ( (k > iMax) || (k > jMax) ) break;
          sqResp = ( bankVetoData->resp->data[k].re *
                     bankVetoData->resp->data[k].re +
                     bankVetoData->resp->data[k].im *
                     bankVetoData->resp->data[k].im ) * dynRange * dynRange;

          if ( bankVetoData->spec->data[k] != 0 && sqResp != 0 )
          {
            /* Convention is that i is 'in the data' */
            Br = bankVetoData->fcInputArray[i]->fcTmplt->data->data[k].re /
                 ( bankVetoData->spec->data[k] * sqResp ) *
                 params->ampVec->data[k] * params->ampVec->data[k];

            Bi = bankVetoData->fcInputArray[i]->fcTmplt->data->data[k].im /
                 ( bankVetoData->spec->data[k] * sqResp ) *
                 params->ampVec->data[k] * params->ampVec->data[k];

	    deltaChirp = bankVetoData->tchirp[i] - bankVetoData->tchirp[j];

            ABr +=    cos(2*LAL_PI*k*deltaF*deltaChirp)*(bankVetoData->fcInputArray[j]->fcTmplt->data->data[k].re*Br
		            +   bankVetoData->fcInputArray[j]->fcTmplt->data->data[k].im*Bi)
	            - sin(2*LAL_PI*k*deltaF*deltaChirp)*(Br*bankVetoData->fcInputArray[j]->fcTmplt->data->data[k].im
		            +   Bi*bankVetoData->fcInputArray[j]->fcTmplt->data->data[k].re);

            ABi +=    cos(2*LAL_PI*k*deltaF*deltaChirp)*(Br*bankVetoData->fcInputArray[j]->fcTmplt->data->data[k].im
		            + Bi*bankVetoData->fcInputArray[j]->fcTmplt->data->data[k].re)
	            + sin(2*LAL_PI*k*deltaF*deltaChirp)*(bankVetoData->fcInputArray[j]->fcTmplt->data->data[k].re*Br
            		    + bankVetoData->fcInputArray[j]->fcTmplt->data->data[k].im*Bi);
          }
          /* if the templates are the same then do an autocorrelation function */
          if (i==j) 
          { 
            bankVetoData->workspace->data[k].re = (REAL4) Br * bankVetoData->fcInputArray[i]->fcTmplt->data->data[k].re 
                                                + Bi * bankVetoData->fcInputArray[i]->fcTmplt->data->data[k].im;
            bankVetoData->workspace->data[k].im = 0.0; /* must be 0 */
          }
        }

        if (i==j)
        {          XLALREAL4ReverseFFT( bankVetoData->acorr, bankVetoData->workspace, bankVetoData->revplan );
          /* since the normalized first sample of the auto correlation is always 1, lets skip it */
	  for (k=1; k < bankVetoData->acorrMatSize + 1; k++)
          {
            /* Save the last acorrMatSize samples */
            bankVetoData->acorrMat->data[i*bankVetoData->acorrMatSize + k-1] = bankVetoData->acorr->data[k] / bankVetoData->acorr->data[0];
          }
        }
        bankVetoData->ccMat->data[i*iSize + j].re = (REAL4) ABr;
        bankVetoData->ccMat->data[j*iSize + i].re = (REAL4) ABr;
        bankVetoData->ccMat->data[i*iSize + j].im = (REAL4) ABi;
        bankVetoData->ccMat->data[j*iSize + i].im = (REAL4) -ABi;

        /* This is hermitian so the imaginary flips sign */
        bankVetoData->normMat->data[i*iSize + j] = (REAL4) ABr;
        /* the convention is that i is in the data */
        bankVetoData->normMat->data[j*iSize + i] = (REAL4) ABr;
      }
    }
  } /* end loop over templates */

  for ( i = 0; i < iSize; i++ ) {
    if (i >= vetoBank->subBankSize) break;
    for ( j = i; j < iSize; j++ ) {
      if (j >= vetoBank->subBankSize) break;
      normfac = (REAL8) 1.0 / sqrt(  bankVetoData->normMat->data[i*iSize + i] * bankVetoData->normMat->data[j*iSize + j]);
      bankVetoData->ccMat->data[i*iSize + j].re *= (REAL4) normfac;
      if (j != i) bankVetoData->ccMat->data[j*iSize + i].re *= (REAL4) normfac;
      bankVetoData->ccMat->data[i*iSize + j].im *= (REAL4) normfac;
      if (j != i) bankVetoData->ccMat->data[j*iSize + i].im *= (REAL4) normfac;
    }
  } /* end norm loop */

}


REAL4
XLALComputeBankVeto( FindChirpBankVetoData *bankVetoData,
                     UINT4 i,
                     UINT4 snrIX,
                     UINT4 *dof)
{

  
  REAL8 iSNR_r = 0;
  REAL8 iSNR_i = 0;
  REAL8 iSNR = 0;
  REAL8 jSNR_r = 0;
  REAL8 jSNR_i = 0;
  REAL8 chisq = 0;
  REAL8 phi = 0;
  UINT8 j = 0;
  COMPLEX16 ij;
  REAL8 ijsq;

  REAL8 denomFac, bankNorm, chi_i, chi_r;
  UINT4 iSize = bankVetoData->length;
  iSNR_r = bankVetoData->qVecArray[i]->data[snrIX].re
         * sqrt(bankVetoData->fcInputArray[i]->fcTmplt->norm);
  iSNR_i = bankVetoData->qVecArray[i]->data[snrIX].im
         * sqrt(bankVetoData->fcInputArray[i]->fcTmplt->norm);
  iSNR = sqrt(iSNR_r * iSNR_r + iSNR_i * iSNR_i);
  phi = atan2(iSNR_i, iSNR_r);
  bankNorm = 0.0;
  denomFac = 0.0;
  *dof = 0;
  if (iSize == 1) return 0;

  for (j = 0; j < iSize; j++)
  {
    if (i == j) continue;
    if (bankVetoData->ccMat->data[i*iSize + j].re == 0 && bankVetoData->ccMat->data[i*iSize + j].im == 0) continue;

    ij.re = bankVetoData->ccMat->data[i*iSize + j].re * cos(phi) - bankVetoData->ccMat->data[i*iSize + j].im * sin(phi);
    ij.im = bankVetoData->ccMat->data[i*iSize + j].im * cos(phi) + bankVetoData->ccMat->data[i*iSize + j].re * sin(phi);
    ijsq = ij.re * ij.re + ij.im * ij.im;

    bankNorm = 2.0 - ijsq;

    jSNR_r = bankVetoData->qVecArray[j]->data[snrIX].re
           * sqrt(bankVetoData->fcInputArray[j]->fcTmplt->norm);
    jSNR_i = bankVetoData->qVecArray[j]->data[snrIX].im
           * sqrt(bankVetoData->fcInputArray[j]->fcTmplt->norm);
      
    chi_r = jSNR_r - ij.re * (iSNR);
    chi_i = jSNR_i - ij.im * (iSNR);
    chisq += (chi_r*chi_r + chi_i*chi_i) / bankNorm;
    (*dof)++;
  }

  /*FIXME Normalization now not necessarily chisquare distributed */
  return (REAL4) chisq;
  /*return (REAL4) chisq / bankNorm;*/
}


InspiralTemplate *
XLALFindChirpSortTemplates( InspiralTemplate *bankHead, UINT4 num, UINT4 subbanksize)
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


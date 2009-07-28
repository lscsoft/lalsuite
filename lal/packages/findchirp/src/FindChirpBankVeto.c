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

NRCSID (FINDCHIRPBANKVETOC, "$Id$");

/* macro to "use" unused function parameters */
#define UNUSED(expr) do { (void)(expr); } while(0)

/* Some static function prototypes */

static int compareTemplateByLevel (const void * a, const void * b);
static int compareTemplateByChirpMass (const void * a, const void * b);
static int compareTemplateByEta (const void * a, const void * b);
static int  breakUpRegions(InspiralTemplate *bankHead, UINT4 num, UINT4 subbanksize );
static InspiralTemplate *
XLALFindChirpSortTemplatesByChirpMass( InspiralTemplate *bankHead, UINT4 num );
static InspiralTemplate *
XLALFindChirpSortTemplatesByEta( InspiralTemplate *bankHead, UINT4 num );
static InspiralTemplate *
XLALFindChirpSortTemplatesByLevel( InspiralTemplate *bankHead, UINT4 num);

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
  REAL4 chisqnorm, tmp, C, chisq, angle, snr, snrk;
  /* input, params and norm are unused */
  UNUSED(input);
  UNUSED(params);
  /* test isn't being done */
  if (!bankVetoData->acorrMat) return 0.0;
  chisq = 0;
  chisqnorm = 0;

  angle = atan2(q[snrIX].im, q[snrIX].re);
  snr = q[snrIX].re * cos(angle) + q[snrIX].im * sin(angle);

  /*Start at k=1 because k=0 is Gauranteed to be zero by construction */
  for (k = 1; k < bankVetoData->acorrMatSize; k++)
  {
    snrk = q[snrIX-k].re * cos(angle) + q[snrIX-k].im * sin(angle);
    C = bankVetoData->acorrMat->data[i * bankVetoData->acorrMatSize + k];
    *dof++;
    /*FIXME what? This should never happen, and it probably doesn't but if it did it would screw the whole value */
    if (C >=1 ) continue;
    tmp = snrk - C * snr;
    chisq += tmp * tmp * norm;
    chisqnorm += 1.0 / (1.0 - C * C);
  }
  chisq /= chisqnorm;
  if (chisq <= 0 )
  {
    fprintf(stderr,"chisq < 0");
    exit(1);
  }
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
  REAL4 ABr, ABi, Br, Bi, sqResp;
  UINT4 stIX = floor( fLow / deltaF );
  /*char fname[30];
  FILE *FP = NULL;*/

  /* FIXME this should be a command line argument */
  bankVetoData->acorrMatSize = 30; /*30 points of autocorrelation function stored */

  /* if there isn't alread memory allocated for workspace and the autocorrelation 
   * then do it
   */
  if ( !bankVetoData->acorr ) bankVetoData->acorr = XLALCreateREAL4Vector((tmpLen-1) * 2);
  if ( !bankVetoData->workspace) bankVetoData->workspace = XLALCreateCOMPLEX8Vector(tmpLen);
  if ( !bankVetoData->acorrMat) bankVetoData->acorrMat = XLALCreateREAL4Vector(bankVetoData->acorrMatSize * iSize);
  if ( !bankVetoData->revplan) bankVetoData->revplan = XLALCreateReverseREAL4FFTPlan((tmpLen-1) * 2 , 0);

  /* deltaT is unused in this function */
  UNUSED(deltaT);

  for ( i = 0; i < iSize; i++ )
  {
    correctFlag = 1;
    for ( j = i; j < iSize; j++ )
    {
      /* FIXME find an XLAL Function to do this?*/
      if (i == j) memset(bankVetoData->workspace->data, 0, bankVetoData->workspace->length * sizeof(COMPLEX8));
      if ( i >= vetoBank->subBankSize )
      {
        bankVetoData->ccMat->data[i*iSize + j] = 0;
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
            /* FIXME set spec to 1 for testing */
             /*bankVetoData->spec->data[k] = 1.0; */
            /* Convention is that i is 'in the data' */
            Br = bankVetoData->fcInputArray[i]->fcTmplt->data->data[k].re /
                 ( bankVetoData->spec->data[k] * sqResp ) *
                 params->ampVec->data[k] * params->ampVec->data[k];

            Bi = bankVetoData->fcInputArray[i]->fcTmplt->data->data[k].im /
                 ( bankVetoData->spec->data[k] * sqResp ) *
                 params->ampVec->data[k] * params->ampVec->data[k];

            ABr += bankVetoData->fcInputArray[j]->fcTmplt->data->data[k].re*Br+
                   bankVetoData->fcInputArray[j]->fcTmplt->data->data[k].im*Bi;

            ABi+=0.0-Br*bankVetoData->fcInputArray[j]->fcTmplt->data->data[k].im
                + Bi*bankVetoData->fcInputArray[j]->fcTmplt->data->data[k].re;
          }
          /* if the templates are the same then do an autocorrelation function */
          if (i==j) 
          { 
            bankVetoData->workspace->data[k].re = Br * bankVetoData->fcInputArray[i]->fcTmplt->data->data[k].re 
                                                + Bi * bankVetoData->fcInputArray[i]->fcTmplt->data->data[k].im;
            bankVetoData->workspace->data[k].im = 0.0; /* must be 0 */
          }
        }

        if (i==j)
        {
          XLALREAL4ReverseFFT( bankVetoData->acorr, bankVetoData->workspace, bankVetoData->revplan );
          /*fprintf(stderr, "\n storing autocorrelation \n");
          sprintf(fname, "chisqtest%d.txt", i);
          FP = fopen(fname, "w");
          for (k=0; k < bankVetoData->acorr->length; k++)
          {
            fprintf(FP,"%d %e\n",k, bankVetoData->acorr->data[k] / bankVetoData->acorr->data[0]);
          }
          fclose(FP);*/
	  for (k=0; k < bankVetoData->acorrMatSize; k++)
          {
            /* Save the last acorrMatSize samples */
            bankVetoData->acorrMat->data[i*bankVetoData->acorrMatSize + k] = bankVetoData->acorr->data[k] / bankVetoData->acorr->data[0];
            if (bankVetoData->acorrMat->data[i*bankVetoData->acorrMatSize + k] > 1.0) 
            {
              fprintf(stderr, "Autocorrelation > 1.0\n");
              exit(1);
            }
           
          }
        }
        bankVetoData->ccMat->data[i*iSize + j] = sqrt(ABr*ABr+ABi*ABi);
        bankVetoData->ccMat->data[j*iSize + i] = sqrt(ABr*ABr+ABi*ABi);
        /* This is hermitian so the imaginary flips sign */
        bankVetoData->normMat->data[i*iSize + j] = atan2(-ABi, ABr);
        /* the convention is that i is in the data */
        bankVetoData->normMat->data[j*iSize + i] = atan2(ABi, ABr);
      }
    }
  }
}


REAL4
XLALComputeBankVeto( FindChirpBankVetoData *bankVetoData,
                     UINT4 i,
                     UINT4 snrIX,
                     UINT4 *dof)
{


  REAL4 iSNR_r = 0;
  REAL4 iSNR_i = 0;
  REAL4 jSNR_r = 0;
  REAL4 jSNR_i = 0;
  REAL4 iSNR = 0;
  REAL4 jSNR = 0;
  REAL4 chisq = 0;
  UINT4 j = 0;
  REAL4 ii,jj,ji,denomFac,phi_i, phi_j, chi_i, chi_r;
  UINT4 iSize = bankVetoData->length;
  iSNR_r = bankVetoData->qVecArray[i]->data[snrIX].re
         * sqrt(bankVetoData->fcInputArray[i]->fcTmplt->norm);
  iSNR_i = bankVetoData->qVecArray[i]->data[snrIX].im
         * sqrt(bankVetoData->fcInputArray[i]->fcTmplt->norm);
  iSNR = sqrt(iSNR_r*iSNR_r + iSNR_i*iSNR_i);

  phi_i =  atan2(iSNR_i,  iSNR_r);

  *dof = 0;
  if (iSize == 1) return 0;

  for (j = 0; j < iSize; j++)
  {
    ji = bankVetoData->ccMat->data[j*iSize + i];
    ii = bankVetoData->ccMat->data[i*iSize + i];
    jj = bankVetoData->ccMat->data[j*iSize + j];
    phi_j = bankVetoData->normMat->data[j*iSize + i];
    denomFac =  1.0 - (ji * ji / ii / jj) ;

    if ( denomFac != 0.0 && ji != 0.0 && ii != 0 && jj != 0 )
    {
      jSNR_r = bankVetoData->qVecArray[j]->data[snrIX].re
             * sqrt(bankVetoData->fcInputArray[j]->fcTmplt->norm);
      jSNR_i = bankVetoData->qVecArray[j]->data[snrIX].im
             * sqrt(bankVetoData->fcInputArray[j]->fcTmplt->norm);
      jSNR = sqrt(jSNR_r*jSNR_r + jSNR_i*jSNR_i);

      chi_r = (jSNR_r*cos(0.0-phi_j) - jSNR_i*sin(0.0-phi_j)) -
              ji / sqrt(ii*jj) * (iSNR_r);

      chi_i = (jSNR_r*sin(0.0-phi_j) + jSNR_i*cos(0.0-phi_j)) -
              ji / sqrt(ii*jj) * (iSNR_i);
      /* Previously I used the modulus only */
      /*  chisq += (jSNR - ji / sqrt(ii*jj) * iSNR)                        */
      /*         * (jSNR - ji / sqrt(ii*jj) * iSNR) / denomFac;            */

      /* Now I use the real and imaginary parts seperately */
      chisq += chi_r*chi_r/denomFac + chi_i*chi_i/denomFac;

      (*dof)++;
    }
  }
  return chisq;
}


InspiralTemplate *
XLALFindChirpSortTemplates( InspiralTemplate *bankHead, UINT4 num, UINT4 subbanksize)
  {
  bankHead = XLALFindChirpSortTemplatesByEta(bankHead,num);
  breakUpRegions(bankHead, num, subbanksize );
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


static int  breakUpRegions(InspiralTemplate *bankHead, UINT4 num, UINT4 subbanksize )
  {
  UINT4 i = 0;
  UINT4 subbanknum = (UINT4) floor(num / (subbanksize+1) );
  if (subbanknum == 0) subbanknum = 1;

  for (i = 0; i < num; bankHead = bankHead->next, i++)
    {
    bankHead->level = (UINT4) (i % subbanknum);
    }
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
  REAL4 mVal1 = (*(InspiralTemplate * const *) a)->level;
  REAL4 mVal2 = (*(InspiralTemplate * const *) b)->level;

  if ( mVal1 > mVal2 ) return 1;
  if ( mVal1 == mVal2 ) return 0;
  if ( mVal1 < mVal2 ) return -1;
  return 0;
}


static int compareTemplateByChirpMass (const void * a, const void * b)
{
  REAL4 mVal1 = (*(InspiralTemplate * const *) a)->chirpMass;
  REAL4 mVal2 = (*(InspiralTemplate * const *) b)->chirpMass;

  if ( mVal1 > mVal2 ) return 1;
  if ( mVal1 == mVal2 ) return 0;
  if ( mVal1 < mVal2 ) return -1;
  return 0;
}

static int compareTemplateByEta (const void * a, const void * b)
{
  const REAL4 mVal1 = ((const InspiralTemplate*)a)->eta;
  const REAL4 mVal2 = ((const InspiralTemplate*)b)->eta;

  if ( mVal1 > mVal2 ) return 1;
  if ( mVal1 == mVal2 ) return 0;
  if ( mVal1 < mVal2 ) return -1;
  return 0;
}


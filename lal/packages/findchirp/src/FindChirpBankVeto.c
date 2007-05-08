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
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>

NRCSID (FINDCHIRPBANKVETOC, "$Id$");

static int compareTemplate (const void * a, const void * b);

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
  static const char *func = "XLALFindChirpCreateSubBanks";
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

void 
XLALBankVetoCCMat ( FindChirpBankVetoData *bankVetoData, 
                    FindChirpSubBank *vetoBank,
                    FindChirpDataParams *params,
                    REAL4 dynRange, 
                    REAL4 fLow,
                    REAL4 deltaF,
                    REAL4 deltaT )

{

  UINT4 i,j,k,iMax,jMax;
  UINT4 iSize = bankVetoData->length;
  UINT4 tmpLen = bankVetoData->fcInputArray[0]->fcTmplt->data->length;
  REAL4 ABr, ABi, Br, Bi, sqResp;
  UINT4 stIX = floor( fLow / deltaF );
  
  for ( i = 0; i < iSize; i++ )
  {
    for ( j = i; j < iSize; j++ )
    {
      if ( i >= vetoBank->subBankSize ) 
      {
        bankVetoData->ccMat->data[i*iSize + j] = 0;
        /* bankVetoData->normMat->data[j] = 0; */
      }
      else
      {
        iMax = bankVetoData->fcInputArray[i]->fcTmplt->tmplt.fFinal / deltaF;
        jMax = bankVetoData->fcInputArray[j]->fcTmplt->tmplt.fFinal / deltaF;
        ABr = ABi = 0;
        for ( k = stIX; k < tmpLen; k++ )
        {
          /* This stops the integration */
          if ( (k > iMax) || (k > jMax) ) break;
          /* remove it possibly */
          sqResp = ( bankVetoData->resp->data[k].re * 
                     bankVetoData->resp->data[k].re +
                     bankVetoData->resp->data[k].im *
                     bankVetoData->resp->data[k].im ) * dynRange * dynRange;

          if ( bankVetoData->spec->data[k] != 0 && sqResp != 0 )
          {
            Br = bankVetoData->fcInputArray[j]->fcTmplt->data->data[k].re /
                 ( bankVetoData->spec->data[k] * sqResp ) *
                 params->ampVec->data[k] * params->ampVec->data[k];
            Bi = bankVetoData->fcInputArray[j]->fcTmplt->data->data[k].im /
                 ( bankVetoData->spec->data[k] * sqResp ) *
                 params->ampVec->data[k] * params->ampVec->data[k];
            ABr += bankVetoData->fcInputArray[i]->fcTmplt->data->data[k].re*Br+
                   bankVetoData->fcInputArray[i]->fcTmplt->data->data[k].im*Bi;
            ABi += bankVetoData->fcInputArray[i]->fcTmplt->data->data[k].re*Bi-
                   bankVetoData->fcInputArray[i]->fcTmplt->data->data[k].im*Br; 
          }
        }  

        bankVetoData->ccMat->data[i*iSize + j] = sqrt(ABr*ABr+ ABi*ABi); 
        bankVetoData->ccMat->data[j*iSize + i] = sqrt(ABr*ABr + ABi*ABi);  
        
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

  REAL4 iSNR = 0;
  REAL4 jSNR = 0;
  REAL4 chisq = 0;
  UINT4 j = 0;
  REAL4 ii,jj,ij,denomFac;
  UINT4 iSize = bankVetoData->length;
  
  iSNR = sqrt( ( bankVetoData->qVecArray[i]->data[snrIX].re *
                 bankVetoData->qVecArray[i]->data[snrIX].re +
                 bankVetoData->qVecArray[i]->data[snrIX].im *
                 bankVetoData->qVecArray[i]->data[snrIX].im ) 
               * bankVetoData->fcInputArray[i]->fcTmplt->norm ); 
  *dof = 0;
  if (iSize == 1) return 0;

  for (j = 0; j < iSize; j++)
  {
    ij = bankVetoData->ccMat->data[i*iSize + j];
    ii = bankVetoData->ccMat->data[i*iSize + i];
    jj = bankVetoData->ccMat->data[j*iSize + j];
    
    denomFac =  ij/sqrt(ii*jj) - sqrt(ii*jj)/ij ;
               /* ij/sqrt(ii*jj) - sqrt(ii*jj)/ij */
               /* ( (ij / ii) - (jj / ij) ) ; */
               /* bankVetoData->fcInputArray[j]->fcTmplt->norm /
               bankVetoData->fcInputArray[i]->fcTmplt->norm; */
    if ( denomFac != 0.0 && ij != 0.0 && ii != 0 && jj != 0 )
    {
      jSNR = sqrt( ( bankVetoData->qVecArray[j]->data[snrIX].re *
                     bankVetoData->qVecArray[j]->data[snrIX].re +
                     bankVetoData->qVecArray[j]->data[snrIX].im *
                     bankVetoData->qVecArray[j]->data[snrIX].im ) 
                   * bankVetoData->fcInputArray[j]->fcTmplt->norm );

      chisq += (iSNR - sqrt(ii*jj) / ij * jSNR) / denomFac *
               (iSNR - sqrt(ii*jj) / ij * jSNR) / denomFac;

      (*dof)++;
    }
  }
  return chisq;
}


InspiralTemplate * 
XLALFindChirpSortTemplates( InspiralTemplate *bankHead, UINT4 num ){

  InspiralTemplate **bankArray = NULL;
  InspiralTemplate *bankFirst = NULL;
  UINT4 i = 0;
  bankFirst = bankHead;
  bankArray = (InspiralTemplate **) LALCalloc(num, sizeof(InspiralTemplate *));


  for (i = 0; (i < num); bankHead = bankHead->next, i++)
  {
    bankArray[i] = bankHead; /* populate pointer array */
  }

  qsort(bankArray, num, sizeof(InspiralTemplate *), compareTemplate);

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

static int compareTemplate (const void * a, const void * b)
{
  REAL4 mVal1 =   (*((InspiralTemplate**)a))->mass1 +
                  (*((InspiralTemplate**)a))->mass2 ;
  REAL4 mVal2 =   (*((InspiralTemplate**)b))->mass1 +
                  (*((InspiralTemplate**)b))->mass2 ;

  if ( mVal1 > mVal2 ) return 1;
  if ( mVal1 == mVal2 ) return 0;
  if ( mVal1 < mVal2 ) return -1;
}


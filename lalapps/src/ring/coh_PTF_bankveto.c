#include "coh_PTF.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>


void initialise_sub_bank(
InspiralTemplate        *PTFBankTemplates,
FindChirpTemplate       *bankFcTmplts,
UINT4                    subBankSize,
UINT4                    numPoints,
UINT4                    spinBank)
{
  UINT4 i;
  unsigned int iseed = (unsigned int)time(NULL);
  srand(iseed);

  REAL4 maxmass1 = 30.;
  REAL4 minmass1 = 1.;
  REAL4 minmass2 = 1.;
  REAL4 maxchi = 1.;
  REAL4 minchi = 0.;
  REAL4 maxkappa = 1.;
  REAL4 minkappa = -1.;
 
  for ( i=0 ; i < subBankSize ; i++ )
  {
    if (spinBank)
      bankFcTmplts[i].PTFQtilde =
        XLALCreateCOMPLEX8VectorSequence( 5, numPoints / 2 + 1 );
    else
      bankFcTmplts[i].PTFQtilde =
        XLALCreateCOMPLEX8VectorSequence( 2, numPoints / 2 + 1 );
    PTFBankTemplates[i].approximant = FindChirpPTF;
    PTFBankTemplates[i].order = LAL_PNORDER_TWO;
    PTFBankTemplates[i].mass1 = pow(rand()/(float)RAND_MAX,2)*(maxmass1-minmass1)+minmass1;
    PTFBankTemplates[i].mass2 = rand()/(float)RAND_MAX*(PTFBankTemplates[i].mass1 - minmass2) + minmass2;
    PTFBankTemplates[i].chi = rand()/(float)RAND_MAX*(maxchi-minchi)+minchi;
    PTFBankTemplates[i].kappa = rand()/(float)RAND_MAX*(maxkappa-minkappa)+minkappa;
    PTFBankTemplates[i].fLower = 38.;
  }
}

REAL4 calculate_bank_veto(
UINT4           numPoints,
UINT4           position,
UINT4           subBankSize,
UINT4           vecLength,
REAL4           a[LAL_NUM_IFO],
REAL4           b[LAL_NUM_IFO],
REAL4           SNR,
REAL8Array      *PTFM[LAL_NUM_IFO+1],
struct coh_PTF_params      *params,
struct bankTemplateOverlaps *bankOverlaps,
struct bankTemplateOverlaps *bankNormOverlaps,
struct bankDataOverlaps *dataOverlaps,
REAL4TimeSeries         *pValues[10],
REAL4TimeSeries         *gammaBeta[2] )
{
  UINT4 ui,uj,uk,ifoNumber,halfNumPoints;
  REAL4 overlapCount,PTFMcomp,pVals1,pVals2;
  REAL4 bankVeto,bankVetoTemp,TjwithS,alpha;
  COMPLEX8 Qoverlap;
  REAL4 overlapQwithQ[subBankSize+1];
  REAL8Array *PTFMtest;
  REAL4 gammaBetaMag,cosPhase,sinPhase;

  gammaBetaMag = pow(gammaBeta[0]->data->data[position - numPoints/4],2);
  gammaBetaMag += pow(gammaBeta[1]->data->data[position - numPoints/4],2);
  gammaBetaMag = pow(gammaBetaMag,0.5);
  cosPhase = gammaBeta[0]->data->data[position - numPoints/4]/gammaBetaMag;
  sinPhase = gammaBeta[1]->data->data[position - numPoints/4]/gammaBetaMag;

  /* We need to calculate normalization factors of (Q_i,Q_i)*/
  halfNumPoints = 3*numPoints/4 - numPoints/4;

  for ( ui = 0 ; ui < subBankSize + 1 ; ui++ )
  {
    overlapCount = 0;
    for ( uj = 0; uj < vecLength; uj++ )
    {
      for ( uk = 0; uk < vecLength ; uk++ )
      {
        for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO ; ifoNumber++ )
        {
          if ( params->haveTrig[ifoNumber] )
          {
            if ( ui == subBankSize )
            {
              PTFMcomp = PTFM[ifoNumber]->data[5 * uj + uk];
            }
            else
            {
              PTFMtest = bankNormOverlaps[ui].PTFM[ifoNumber];
              PTFMcomp = PTFMtest->data[vecLength * uj + uk];
            }
            pVals1 = a[ifoNumber]+pValues[uj]->data->data[position-numPoints/4];
            pVals1 += b[ifoNumber]+pValues[uj+vecLength]->data->data[position-numPoints/4];
            pVals2 = a[ifoNumber]+pValues[uk]->data->data[position-numPoints/4];
            pVals2 += b[ifoNumber]+pValues[uk+vecLength]->data->data[position-numPoints/4];
            overlapCount += pVals1*pVals2*PTFMcomp;
//            fprintf(stderr,"Overlap contributions: %e %e %e %d %d\n",pVals1,pVals2,PTFMcomp,ui,subBankSize);
          }
        }
      }
    } 
    overlapQwithQ[ui] = pow(overlapCount,0.5);
//    fprintf(stderr,"Norm factors: %e \n", overlapQwithQ[ui]);
  }    

  /* Now we can calculate the bank veto itself */
  bankVeto = 0;

  for ( ui = 0 ; ui < subBankSize ; ui++ )
  {
    /* Calculate the SNR between the bank template and the noise */
    TjwithS = 0;
    for ( uj = 0; uj < vecLength; uj++ )
    {
      for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO ; ifoNumber++ )
      {
        if ( params->haveTrig[ifoNumber] )
        {
          pVals1 = a[ifoNumber]+pValues[uj]->data->data[position-numPoints/4];
          pVals1 += b[ifoNumber]+pValues[uj+vecLength]->data->data[position-numPoints/4];
          if (vecLength == 2 && uj == 1 )
          {
            /* Special case for non spinning: to save memory Q2 is not made*/
            Qoverlap = dataOverlaps[ui].PTFqVec[ifoNumber]->data[0*halfNumPoints + (position - numPoints/4)];
            TjwithS += pVals1 * (Qoverlap.im*cosPhase - Qoverlap.re * sinPhase);
          }
          else
          {
            Qoverlap = dataOverlaps[ui].PTFqVec[ifoNumber]->data[uj*halfNumPoints + (position - numPoints/4)];
            TjwithS += pVals1 * (Qoverlap.re*cosPhase + Qoverlap.im * sinPhase);
          }
//          fprintf(stderr,"TjwithS cont: %e %e %e %e %e\n,",Qoverlap.re,Qoverlap.im,pVals1,cosPhase,sinPhase);
        }
      }
    }
    TjwithS = TjwithS / overlapQwithQ[ui];
 
    /* Calculate the overlap between the bank template and the interesting one*/
    alpha = 0;
    for ( uj = 0; uj < vecLength; uj++ )
    {
      for ( uk = 0; uk < vecLength ; uk++ )
      {
        for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO ; ifoNumber++ )
        {
          if ( params->haveTrig[ifoNumber] )
          {
            PTFMtest = bankOverlaps[ui].PTFM[ifoNumber];
            PTFMcomp = PTFMtest->data[vecLength * uj + uk];
            pVals1 = a[ifoNumber]+pValues[uj]->data->data[position-numPoints/4];
            pVals1 += b[ifoNumber]+pValues[uj+vecLength]->data->data[position-numPoints/4];
            pVals2 = a[ifoNumber]+pValues[uk]->data->data[position-numPoints/4];
            pVals2 += b[ifoNumber]+pValues[uk+vecLength]->data->data[position-numPoints/4];
            alpha += pVals1*pVals2 * PTFMcomp;
//            fprintf(stderr,"alpha: %e %e %e \n",PTFMcomp,pVals1,pVals2);
          }
        }
      }
    }
    alpha = alpha / (overlapQwithQ[ui]*overlapQwithQ[subBankSize]);
      
    bankVetoTemp = TjwithS - alpha * SNR;
    bankVetoTemp = pow(bankVetoTemp,2)/pow(1-alpha,2);
    bankVeto += bankVetoTemp;
//    fprintf(stderr,"BV comps: %e %e %e \n",TjwithS,alpha,SNR);
  }

//  fprintf(stderr,"Bank Veto %e \n", bankVeto);  

  return bankVeto;

}

void free_bank_veto_memory(
  struct bankTemplateOverlaps *bankNormOverlaps,
  InspiralTemplate        *PTFBankTemplates,
  FindChirpTemplate       *bankFcTmplts,
  UINT4 subBankSize,
  struct bankTemplateOverlaps *bankOverlaps,
  struct bankDataOverlaps *dataOverlaps)
{
  UINT4 ui,ifoNumber;
 
  if ( bankNormOverlaps )
  {
    for ( ui = 0 ; ui < subBankSize ; ui++ )
    {
      for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        if ( bankNormOverlaps[ui].PTFM[ifoNumber] )
        {
          XLALDestroyREAL8Array(bankNormOverlaps[ui].PTFM[ifoNumber]);
        }
      }
    }
    LALFree(bankNormOverlaps);
  }

  if ( PTFBankTemplates )
    LALFree( PTFBankTemplates);

  if ( bankFcTmplts )
  {
    for ( ui = 0 ; ui < subBankSize ; ui++ )
    {
      if ( bankFcTmplts[ui].PTFQtilde )
        XLALDestroyCOMPLEX8VectorSequence( bankFcTmplts[ui].PTFQtilde );
    }
    LALFree( bankFcTmplts );
  } 

  if (dataOverlaps)
    LALFree(dataOverlaps);
  if (bankOverlaps)
    LALFree(bankOverlaps);
}

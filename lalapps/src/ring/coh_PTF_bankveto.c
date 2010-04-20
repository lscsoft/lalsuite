#include "coh_PTF.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>


void initialise_sub_bank(
struct coh_PTF_params   *params,
InspiralTemplate        *PTFBankTemplates,
FindChirpTemplate       *bankFcTmplts,
UINT4                    subBankSize,
UINT4                    numPoints,
UINT4                    spinBank)
{
  UINT4 i;
  srand(params->randomSeed);

  REAL4 maxmass1 = 30.;
  REAL4 minmass1 = 28.;
  REAL4 minmass2 = 28.;
/*  REAL4 maxmass1 = 30.;
  REAL4 minmass1 = 29.;
  REAL4 minmass2 = 28.;*/
  REAL4 maxchi = 0.1;
  REAL4 minchi = 0.;
  REAL4 maxkappa = 1.;
  REAL4 minkappa = 0.9;
 
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
REAL4TimeSeries         *gammaBeta[2],
COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO+1],
INT4            timeOffsetPoints[LAL_NUM_IFO] )
{
  UINT4 ui,uj,uk,ifoNumber,halfNumPoints,bankVecLength,calTimeOffset;
  REAL4 overlapCount,PTFMcomp,pVals1,pVals2;
  REAL4 bankVeto,bankVetoTemp,TjwithS,TjwithS2,alpha;
  COMPLEX8 Qoverlap,cSNR;
  REAL4 overlapQwithQ[subBankSize+1];
  REAL8Array *PTFMtest;
  REAL4 gammaBetaMag,cosPhase,sinPhase;

  if ( params->spinBank )
    bankVecLength = 5;
  else
    bankVecLength = 2;

  gammaBetaMag = pow(gammaBeta[0]->data->data[position - numPoints/4],2);
  gammaBetaMag += pow(gammaBeta[1]->data->data[position - numPoints/4],2);
  gammaBetaMag = pow(gammaBetaMag,0.5);
  cosPhase = gammaBeta[0]->data->data[position - numPoints/4]/gammaBetaMag;
  sinPhase = gammaBeta[1]->data->data[position - numPoints/4]/gammaBetaMag;

  /* We need to calculate normalization factors of (Q_i,Q_i)*/
  halfNumPoints = 3*numPoints/4 - numPoints/4 + 10000;

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
              PTFMcomp = PTFMtest->data[bankVecLength * uj + uk];
            }
            pVals1 = a[ifoNumber]*pValues[uj]->data->data[position-numPoints/4];
            pVals1 += b[ifoNumber]*pValues[uj+vecLength]->data->data[position-numPoints/4];
            pVals2 = a[ifoNumber]*pValues[uk]->data->data[position-numPoints/4];
            pVals2 += b[ifoNumber]*pValues[uk+vecLength]->data->data[position-numPoints/4];
            overlapCount += pVals1*pVals2*PTFMcomp;
//            fprintf(stderr,"Overlap contributions: %e %e %e %d %d %e %e\n",pVals1,pVals2,PTFMcomp,ui,subBankSize,pValues[uj]->data->data[position-numPoints/4],pValues[uj]->data->data[position-numPoints/4]);
          }
        }
      }
    } 
    overlapQwithQ[ui] = pow(overlapCount,0.5);
//    fprintf(stderr,"Norm factors: %e %e \n", overlapCount,overlapQwithQ[ui]);
  }    

  /* Now we can calculate the bank veto itself */
  bankVeto = 0;

  for ( ui = 0 ; ui < subBankSize ; ui++ )
  {
    /* Calculate the SNR between the bank template and the noise */
    TjwithS = 0;
    TjwithS2 = 0;
    for ( uj = 0; uj < vecLength; uj++ )
    {
      for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO ; ifoNumber++ )
      {
        if ( params->haveTrig[ifoNumber] )
        {
          pVals1 = a[ifoNumber]*pValues[uj]->data->data[position-numPoints/4];
          pVals1 += b[ifoNumber]*pValues[uj+vecLength]->data->data[position-numPoints/4];
          if ( position+timeOffsetPoints[ifoNumber] >= 3*numPoints/4 +5000)
          {
            calTimeOffset = 3*numPoints/4 -1;
            fprintf(stderr,"Overflow occured in time shifting in bank veto\n");
          }
          else if  ( position+timeOffsetPoints[ifoNumber] < numPoints/4 - 5000)
          {
            calTimeOffset = numPoints/4;
            fprintf(stderr,"Overflow occured in time shifting in bank veto\n");
          }
          else
            calTimeOffset =position+timeOffsetPoints[ifoNumber];

          if (vecLength == 2 && uj == 1 )
          {
            /* Special case for non spinning: to save memory Q2 is not made*/
            /* Note that (s|H) is calculated not (H|s) */
            Qoverlap = PTFqVec[ifoNumber]->data[1*numPoints + calTimeOffset];
            TjwithS2 += pVals1 * (Qoverlap.re*cosPhase + Qoverlap.im * sinPhase);
            Qoverlap = dataOverlaps[ui].PTFqVec[ifoNumber]->data[(calTimeOffset - numPoints/4+5000)];
            TjwithS += pVals1 * (-Qoverlap.im*cosPhase + Qoverlap.re * sinPhase);
          }
          else
          {
            Qoverlap = PTFqVec[ifoNumber]->data[uj*numPoints + calTimeOffset];
            TjwithS2 += pVals1 * (Qoverlap.re*cosPhase + Qoverlap.im * sinPhase);
            Qoverlap = dataOverlaps[ui].PTFqVec[ifoNumber]->data[uj*halfNumPoints + (calTimeOffset - numPoints/4+5000)];
            TjwithS += pVals1 * (Qoverlap.re*cosPhase + Qoverlap.im * sinPhase);
          }
//          fprintf(stderr,"TjwithS cont: %e %e %e %e %e %e %e %e %e %d\n,",Qoverlap.re,Qoverlap.im,PTFqVec[ifoNumber]->data[uj*numPoints + calTimeOffset].re,PTFqVec[ifoNumber]->data[uj*numPoints + calTimeOffset].im,pVals1,cosPhase,sinPhase,TjwithS,TjwithS2,calTimeOffset);
        }
      }
    }
    TjwithS = TjwithS / overlapQwithQ[ui];
    TjwithS2 = TjwithS2 / overlapQwithQ[subBankSize];
 
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
            PTFMcomp = PTFMtest->data[bankVecLength * uj + uk];
            pVals1 = a[ifoNumber]*pValues[uj]->data->data[position-numPoints/4];
            pVals1 += b[ifoNumber]*pValues[uj+vecLength]->data->data[position-numPoints/4];
            pVals2 = a[ifoNumber]*pValues[uk]->data->data[position-numPoints/4];
            pVals2 += b[ifoNumber]*pValues[uk+vecLength]->data->data[position-numPoints/4];
            alpha += pVals1*pVals2 * PTFMcomp;
//            fprintf(stderr,"alpha: %e %e %e \n",PTFMcomp,pVals1,pVals2);
          }
        }
      }
    }
    alpha = alpha / (overlapQwithQ[ui]*overlapQwithQ[subBankSize]);
      
//    fprintf(stderr,"BankVeto : %e \n", bankVeto);
    bankVetoTemp = TjwithS - alpha * SNR;
//    bankVetoTemp = pow(bankVetoTemp,2)/(1-pow(alpha,2));

    bankVetoTemp = bankVetoTemp/pow(1-pow(alpha,2),0.5);

    bankVeto += bankVetoTemp;

    fprintf(stderr,"BV comps: %e %e %e %e\n",TjwithS,alpha,SNR,TjwithS2);
  }

  fprintf(stderr,"Bank Veto %e \n", bankVeto);  

  return bankVeto;

}

REAL4 calculate_bank_veto_max_phase(
UINT4           numPoints,
UINT4           position,
UINT4           subBankSize,
UINT4           vecLength,
REAL4           a[LAL_NUM_IFO],
REAL4           b[LAL_NUM_IFO],
REAL4           SNR,
REAL8Array      *PTFM[LAL_NUM_IFO+1],
struct coh_PTF_params      *params,
struct bankComplexTemplateOverlaps *bankOverlaps,
struct bankTemplateOverlaps *bankNormOverlaps,
struct bankDataOverlaps *dataOverlaps,
REAL4TimeSeries         *pValues[10],
REAL4TimeSeries         *gammaBeta[2],
INT4            timeOffsetPoints[LAL_NUM_IFO],
UINT4           singleDetector )
{
//  fprintf(stderr,"Entering bank veto calculator\n");
  UINT4 ui,uj,uk,ifoNumber,halfNumPoints,bankVecLength,calTimeOffset;
  REAL4 overlapCount,PTFMcomp,pVals1,pVals2;
  REAL4 bankVeto;
  COMPLEX8 Qoverlap,TjwithS,PTFMComplexcomp,alpha,bankVetoTemp,cSNR;
  REAL4 overlapQwithQ[subBankSize+1];
  REAL8Array *PTFMtest;
  COMPLEX8Array *PTFMComptest;
  REAL4 gammaBetaMag,cosPhase,sinPhase;

  if ( params->spinBank )
    bankVecLength = 5;
  else
    bankVecLength = 2;

  gammaBetaMag = pow(gammaBeta[0]->data->data[position - numPoints/4],2);
  gammaBetaMag += pow(gammaBeta[1]->data->data[position - numPoints/4],2);
  gammaBetaMag = pow(gammaBetaMag,0.5);
  cosPhase = gammaBeta[0]->data->data[position - numPoints/4]/gammaBetaMag;
  sinPhase = gammaBeta[1]->data->data[position - numPoints/4]/gammaBetaMag;
  cSNR.re = SNR * cosPhase;
  cSNR.im = SNR * sinPhase;  

  /* We need to calculate normalization factors of (Q_i,Q_i)*/
  halfNumPoints = 3*numPoints/4 - numPoints/4 + 10000;

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
              PTFMcomp = PTFMtest->data[bankVecLength * uj + uk];
            }
            pVals1 = a[ifoNumber]*pValues[uj]->data->data[position-numPoints/4];
            if (! singleDetector)
              pVals1 += b[ifoNumber]*pValues[uj+vecLength]->data->data[position-numPoints/4];
            pVals2 = a[ifoNumber]*pValues[uk]->data->data[position-numPoints/4];
            if (! singleDetector)
              pVals2 += b[ifoNumber]*pValues[uk+vecLength]->data->data[position-numPoints/4];
            overlapCount += pVals1*pVals2*PTFMcomp;
            fprintf(stderr,"Overlap contributions: %e %e %e %d %d\n",pVals1,pVals2,PTFMcomp,ui,subBankSize);
          }
        }
      }
    }
    overlapQwithQ[ui] = pow(overlapCount,0.5);
    fprintf(stderr,"Norm factors: %e %e \n", overlapCount,overlapQwithQ[ui]);
 }

  /* Now we can calculate the bank veto itself */
  bankVeto = 0;

  for ( ui = 0 ; ui < subBankSize ; ui++ )
  {
    /* Calculate the SNR between the bank template and the noise */
    TjwithS.re = 0;
    TjwithS.im = 0;
    for ( uj = 0; uj < vecLength; uj++ )
    {
      for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO ; ifoNumber++ )
      {
        if ( params->haveTrig[ifoNumber] )
        {
          pVals1 = a[ifoNumber]*pValues[uj]->data->data[position-numPoints/4];
          if (! singleDetector)
            pVals1 += b[ifoNumber]*pValues[uj+vecLength]->data->data[position-numPoints/4];
          if ( position+timeOffsetPoints[ifoNumber] >= 3*numPoints/4 +5000)
          {
            calTimeOffset = 3*numPoints/4 +4999;
            fprintf(stderr,"Overflow occured in time shifting in bank veto\n");
          }
          else if  ( position+timeOffsetPoints[ifoNumber] < numPoints/4 - 5000)
          {
            calTimeOffset = numPoints/4 - 4999;
            fprintf(stderr,"Overflow occured in time shifting in bank veto\n");
          }
          else
            calTimeOffset =position+timeOffsetPoints[ifoNumber];

          if (vecLength == 2 && uj == 1 )
          {
            /* Special case for non spinning: to save memory Q2 is not made*/
            Qoverlap = dataOverlaps[ui].PTFqVec[ifoNumber]->data[ (calTimeOffset - numPoints/4 + 5000)];
            TjwithS.re += pVals1 * Qoverlap.im; 
            TjwithS.im += -pVals1 * Qoverlap.re;
          }
          else
          {
            Qoverlap = dataOverlaps[ui].PTFqVec[ifoNumber]->data[uj*halfNumPoints + (calTimeOffset - numPoints/4 + 5000)];
            TjwithS.re += pVals1 * Qoverlap.re;
            TjwithS.im += pVals1 * Qoverlap.im;
          }
          fprintf(stderr,"TjwithS cont: %e %e %e %e %e\n,",Qoverlap.re,Qoverlap.im,pVals1,cosPhase,sinPhase);
        }
      }
    }
    TjwithS.re = TjwithS.re / overlapQwithQ[ui];
    TjwithS.im = TjwithS.im / overlapQwithQ[ui];

    /* Calculate the overlap between the bank template and the interesting one*/
    alpha.re = 0;
    alpha.im = 0;
    for ( uj = 0; uj < vecLength; uj++ )
    {
      for ( uk = 0; uk < vecLength ; uk++ )
      {
        for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO ; ifoNumber++ )
        {
          if ( params->haveTrig[ifoNumber] )
          {
            PTFMComptest = bankOverlaps[ui].PTFM[ifoNumber];
            PTFMComplexcomp = PTFMComptest->data[bankVecLength * uj + uk];
            pVals1 = a[ifoNumber]*pValues[uj]->data->data[position-numPoints/4];
            if (! singleDetector)
              pVals1 += b[ifoNumber]*pValues[uj+vecLength]->data->data[position-numPoints/4];
            pVals2 = a[ifoNumber]*pValues[uk]->data->data[position-numPoints/4];
            if (! singleDetector)
              pVals2 += b[ifoNumber]*pValues[uk+vecLength]->data->data[position-numPoints/4];
            alpha.re += pVals1*pVals2 * PTFMComplexcomp.re;
            alpha.im += pVals1*pVals2 * PTFMComplexcomp.im;
            fprintf(stderr,"alpha: %e %e %e %e \n",PTFMComplexcomp.re,PTFMComplexcomp.im,pVals1,pVals2);
          }
        }
      }
    }
    alpha.re = alpha.re / (overlapQwithQ[ui]*overlapQwithQ[subBankSize]);
    alpha.im = alpha.im / (overlapQwithQ[ui]*overlapQwithQ[subBankSize]);

    bankVetoTemp.re = TjwithS.re - alpha.re * cSNR.re + alpha.im * cSNR.im;
    bankVetoTemp.im = TjwithS.im - alpha.im * cSNR.re - alpha.re * cSNR.im;
    bankVetoTemp.re = bankVetoTemp.re / pow( 1 - pow(alpha.re,2) - pow(alpha.im,2),0.5);
    bankVetoTemp.im = bankVetoTemp.im / pow( 1 - pow(alpha.re,2) - pow(alpha.im,2),0.5);
    bankVeto += pow(bankVetoTemp.re,2) + pow(bankVetoTemp.im,2);
    fprintf(stderr,"BV comps: %e %e %e %e %e %e\n",TjwithS.re,TjwithS.im,alpha.re,alpha.im,cSNR.re,cSNR.im);
    fprintf(stderr,"BankVeto : %e \n", bankVeto);
  }

  fprintf(stderr,"Bank Veto %e \n", bankVeto);  

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

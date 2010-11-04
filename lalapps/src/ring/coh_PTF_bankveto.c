#include "coh_PTF.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

UINT4 read_sub_bank(
struct coh_PTF_params   *params,
InspiralTemplate        **PTFBankTemplates)
{
  UINT4 i,numTemplates;
  InspiralTemplate *bankTemplate;

  numTemplates = InspiralTmpltBankFromLIGOLw( PTFBankTemplates,
      params->bankVetoBankName,-1, -1 );

  bankTemplate = *PTFBankTemplates;

  for (i=0; (i < numTemplates); bankTemplate = bankTemplate->next,i++)
  {
    bankTemplate->fLower = 38.;
    bankTemplate->approximant = FindChirpPTF;
    bankTemplate->order = LAL_PNORDER_TWO;
  }
  return numTemplates;
}

/* FIXME: parameter spinBank is unused */
void initialise_sub_bank(
struct coh_PTF_params   *params,
InspiralTemplate        *PTFBankTemplates,
FindChirpTemplate       *bankFcTmplts,
UINT4                    subBankSize,
UINT4                    numPoints,
UINT4                    UNUSED spinBank)
{
  spinBank = 0;
  UINT4 i;
  srand(params->randomSeed);

  REAL4 maxmass1 = 30.;
  REAL4 minmass1 = 1.;
  REAL4 minmass2 = 1.;
/*  REAL4 maxmass1 = 30.;
  REAL4 minmass1 = 29.;
  REAL4 minmass2 = 28.;*/
  REAL4 maxchi = 0.;
  REAL4 minchi = 0.;
  REAL4 maxkappa = 1.;
  REAL4 minkappa = 0.9;
 
  for ( i=0 ; i < subBankSize ; i++ )
  {
/*    if (spinBank)
      bankFcTmplts[i].PTFQtilde =
        XLALCreateCOMPLEX8VectorSequence( 5, numPoints / 2 + 1 );
    else*/
    bankFcTmplts[i].PTFQtilde =
      XLALCreateCOMPLEX8VectorSequence( 1, numPoints / 2 + 1 );
    PTFBankTemplates[i].approximant = FindChirpPTF;
    PTFBankTemplates[i].order = LAL_PNORDER_TWO;
    PTFBankTemplates[i].mass1 = pow(rand()/(float)RAND_MAX,2)*(maxmass1-minmass1)+minmass1;
    PTFBankTemplates[i].mass2 = rand()/(float)RAND_MAX*(PTFBankTemplates[i].mass1 - minmass2) + minmass2;
    PTFBankTemplates[i].chi = rand()/(float)RAND_MAX*(maxchi-minchi)+minchi;
    PTFBankTemplates[i].kappa = rand()/(float)RAND_MAX*(maxkappa-minkappa)+minkappa;
    PTFBankTemplates[i].fLower = 38.;
    fprintf(stderr,"Masses %e %e \n",PTFBankTemplates[i].mass1,PTFBankTemplates[i].mass2);
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
struct bankComplexTemplateOverlaps *bankOverlaps,
struct bankTemplateOverlaps *bankNormOverlaps,
struct bankDataOverlaps *dataOverlaps,
REAL4TimeSeries         *pValues[10],
REAL4TimeSeries         *gammaBeta[2],
COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO+1],
INT4            timeOffsetPoints[LAL_NUM_IFO],
UINT4           singleDetector )
{
  UINT4 ui,uj,uk,ifoNumber,halfNumPoints,bankVecLength,calTimeOffset;
  REAL4 overlapCount,PTFMcomp,pVals1,pVals2;
  REAL4 bankVeto,bankVetoTemp,TjwithS,TjwithS2;
  COMPLEX8 Qoverlap,alpha,PTFMComplexcomp;
  REAL4 overlapQwithQ[subBankSize+1];
  REAL8Array *PTFMtest;
  COMPLEX8Array *PTFMCompTest;
  REAL4 gammaBetaMag,cosPhase,sinPhase;

 /* if ( params->spinBank )
    bankVecLength = 5;
  else*/
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
    for ( uj = 0; uj < bankVecLength; uj++ )
    {
      for ( uk = 0; uk < bankVecLength ; uk++ )
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
              if (uj != uk)
                PTFMcomp = 0;
              else
              {
                PTFMtest = bankNormOverlaps[ui].PTFM[ifoNumber];
                PTFMcomp = PTFMtest->data[0];
              }
            }
            pVals1 = a[ifoNumber]*pValues[uj]->data->data[position-numPoints/4];
            if (! singleDetector)
              pVals1 += b[ifoNumber]*pValues[uj+vecLength]->data->data[position-numPoints/4];
            pVals2 = a[ifoNumber]*pValues[uk]->data->data[position-numPoints/4];
            if (! singleDetector)
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
    for ( uj = 0; uj < bankVecLength; uj++ )
    {
      for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO ; ifoNumber++ )
      {
        if ( params->haveTrig[ifoNumber] )
        {
          pVals1 = a[ifoNumber]*pValues[uj]->data->data[position-numPoints/4];
          if (! singleDetector )
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

          if (bankVecLength == 2 && uj == 1 )
          {
            /* Special case for non spinning: to save memory Q2 is not made*/
            /* Note that (s|H) is calculated not (H|s) */
            Qoverlap = PTFqVec[ifoNumber]->data[0*numPoints + calTimeOffset];
            TjwithS2 += pVals1 * (-Qoverlap.im*cosPhase + Qoverlap.re * sinPhase);
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
    alpha.re = 0;
    alpha.im = 0;
    for ( uj = 0; uj < bankVecLength; uj++ )
    {
      for ( uk = 0; uk < bankVecLength ; uk++ )
      {
        for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO ; ifoNumber++ )
        {
          if ( params->haveTrig[ifoNumber] )
          {
            PTFMCompTest = bankOverlaps[ui].PTFM[ifoNumber];
            if (uj == uk)           
              PTFMComplexcomp = PTFMCompTest->data[0];
            else 
            {
              PTFMComplexcomp.im = PTFMCompTest->data[0].re;
              PTFMComplexcomp.re = PTFMCompTest->data[0].im;
            }
            if (uk > uj)
            {
              PTFMComplexcomp.re= - PTFMComplexcomp.re;
              PTFMComplexcomp.im= - PTFMComplexcomp.im;
            }
            pVals1 = a[ifoNumber]*pValues[uj]->data->data[position-numPoints/4];
            if (! singleDetector)
              pVals1 += b[ifoNumber]*pValues[uj+vecLength]->data->data[position-numPoints/4];
            pVals2 = a[ifoNumber]*pValues[uk]->data->data[position-numPoints/4];
            if (! singleDetector)
              pVals2 += b[ifoNumber]*pValues[uk+vecLength]->data->data[position-numPoints/4];
            alpha.re += pVals1*pVals2 * PTFMComplexcomp.re;
            alpha.im += pVals1*pVals2 * PTFMComplexcomp.im;
//            fprintf(stderr,"alpha: %e %e %e \n",PTFMComplexcomp.re,pVals1,pVals2);
          }
        }
      }
    }
    alpha.re = alpha.re / (overlapQwithQ[ui]*overlapQwithQ[subBankSize]);
    alpha.im = alpha.im / (overlapQwithQ[ui]*overlapQwithQ[subBankSize]);
      
//    fprintf(stderr,"BankVeto : %e \n", bankVeto);
    bankVetoTemp = TjwithS - alpha.re * SNR;
    bankVetoTemp = pow(bankVetoTemp,2)/(1-pow(alpha.re,2) - pow(alpha.im,2));

//    bankVetoTemp = bankVetoTemp/pow(1-pow(alpha.re,2) - pow(alpha.im,2),0.5);

    bankVeto += bankVetoTemp;

//    fprintf(stderr,"BV comps: %e %e %e %e %e\n",TjwithS,alpha.re,alpha.im,SNR,TjwithS2);
  }

//  fprintf(stderr,"Bank Veto %e \n", bankVeto);  

  return bankVeto;

}

REAL4 calculate_bank_veto_max_phase(
UINT4           numPoints,
UINT4           position,
UINT4           subBankSize,
REAL8Array      *PTFM[LAL_NUM_IFO+1],
struct coh_PTF_params      *params,
struct bankComplexTemplateOverlaps *bankOverlaps,
struct bankTemplateOverlaps *bankNormOverlaps,
struct bankDataOverlaps *dataOverlaps,
COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO+1],
INT4            timeOffsetPoints[LAL_NUM_IFO]
)
{
//  fprintf(stderr,"Entering bank veto calculator\n");
  UINT4 ui, ifoNumber, UNUSED halfNumPoints, UNUSED bankVecLength, calTimeOffset;
  REAL4 overlapCount,PTFMcomp;
  REAL4 bankVeto;
  COMPLEX8 Qoverlap,TjwithS,PTFMComplexcomp,alpha,bankVetoTemp,cSNR;
  REAL4 overlapQwithQ[subBankSize+1];
  REAL8Array *PTFMtest;
  COMPLEX8Array *PTFMComptest;

  if ( params->spinBank )
    bankVecLength = 5;
  else
    bankVecLength = 2;

  /* We need to seperate SNR into a complex number */

  cSNR.re = 0;
  cSNR.im = 0;

  for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO ; ifoNumber++ )
  {
    if ( params->haveTrig[ifoNumber] )
    {
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
       
      Qoverlap = PTFqVec[ifoNumber]->data[calTimeOffset];
      cSNR.re +=  Qoverlap.re;
      cSNR.im +=  Qoverlap.im; 
    }
  }


  /* We need to calculate normalization factors of (Q_i,Q_i)*/
  halfNumPoints = 3*numPoints/4 - numPoints/4 + 10000;

  for ( ui = 0 ; ui < subBankSize + 1 ; ui++ )
  {
    overlapCount = 0;
    for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO ; ifoNumber++ )
    {
      if ( params->haveTrig[ifoNumber] )
      {
        if ( ui == subBankSize )
        {
          PTFMcomp = PTFM[ifoNumber]->data[0];
        }
        else
        {
          PTFMtest = bankNormOverlaps[ui].PTFM[ifoNumber];
          PTFMcomp = PTFMtest->data[0];
        }
        overlapCount += PTFMcomp;
//      fprintf(stderr,"Overlap contributions: %e %e %e %d %d\n",pVals1,pVals2,PTFMcomp,ui,subBankSize);
      }
    }
    overlapQwithQ[ui] = pow(overlapCount,0.5);
//  fprintf(stderr,"Norm factors: %e %e \n", overlapCount,overlapQwithQ[ui]);
 }

  /* Now we can calculate the bank veto itself */
  bankVeto = 0;

  cSNR.re = cSNR.re / overlapQwithQ[subBankSize];
  cSNR.im = cSNR.im / overlapQwithQ[subBankSize];

  for ( ui = 0 ; ui < subBankSize ; ui++ )
  {
    /* Calculate the SNR between the bank template and the noise */
    TjwithS.re = 0;
    TjwithS.im = 0;
    for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO ; ifoNumber++ )
    {
      if ( params->haveTrig[ifoNumber] )
      {
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

        Qoverlap = dataOverlaps[ui].PTFqVec[ifoNumber]->data[(calTimeOffset - numPoints/4 + 5000)];
        TjwithS.re += Qoverlap.re;
        TjwithS.im += Qoverlap.im;
//      fprintf(stderr,"TjwithS cont: %e %e %e %e %e\n,",Qoverlap.re,Qoverlap.im,pVals1,cosPhase,sinPhase);
      }
    }
    TjwithS.re = TjwithS.re / overlapQwithQ[ui];
    TjwithS.im = TjwithS.im / overlapQwithQ[ui];

    /* Calculate the overlap between the bank template and the interesting one*/
    alpha.re = 0;
    alpha.im = 0;
    for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO ; ifoNumber++ )
    {
      if ( params->haveTrig[ifoNumber] )
      {
        PTFMComptest = bankOverlaps[ui].PTFM[ifoNumber];
        PTFMComplexcomp = PTFMComptest->data[0];
        alpha.re += PTFMComplexcomp.re;
        alpha.im += PTFMComplexcomp.im;
//      fprintf(stderr,"alpha: %e %e %e %e \n",PTFMComplexcomp.re,PTFMComplexcomp.im,pVals1,pVals2);
      }
    }
    alpha.re = alpha.re / (overlapQwithQ[ui]*overlapQwithQ[subBankSize]);
    alpha.im = alpha.im / (overlapQwithQ[ui]*overlapQwithQ[subBankSize]);

    bankVetoTemp.re = TjwithS.re - alpha.re * cSNR.re - alpha.im * cSNR.im;
    bankVetoTemp.im = TjwithS.im + alpha.im * cSNR.re - alpha.re * cSNR.im;
    bankVetoTemp.re = bankVetoTemp.re / pow( 1 - pow(alpha.re,2) - pow(alpha.im,2),0.5);
    bankVetoTemp.im = bankVetoTemp.im / pow( 1 - pow(alpha.re,2) - pow(alpha.im,2),0.5);
    bankVeto += pow(bankVetoTemp.re,2) + pow(bankVetoTemp.im,2);
//    fprintf(stderr,"BV comps: %e %e %e %e %e %e\n",TjwithS.re,TjwithS.im,alpha.re,alpha.im,cSNR.re,cSNR.im);
//    fprintf(stderr,"BankVeto : %e \n", bankVeto);
//    fprintf(stderr,"%e %e %e %e %e\n", bankVetoTemp.re,bankVetoTemp.im,pValues[0]->data->data[position-numPoints/4],pValues[1]->data->data[position-numPoints/4],SNR);
  }

//  fprintf(stderr,"Bank Veto %e \n", bankVeto);  

  return bankVeto;

}

REAL4 calculate_bank_veto_max_phase_coherent(
UINT4           numPoints,
UINT4           position,
UINT4           subBankSize,
REAL4           a[LAL_NUM_IFO],
REAL4           b[LAL_NUM_IFO],
struct coh_PTF_params      *params,
struct bankCohTemplateOverlaps *cohBankOverlaps,
struct bankDataOverlaps *dataOverlaps,
COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO+1],
INT4            timeOffsetPoints[LAL_NUM_IFO],
gsl_matrix *Bankeigenvecs[50],
gsl_vector *Bankeigenvals[50]
)
{
  UINT4 ui,uj,uk,ifoNumber,halfNumPoints;
  INT4 calTimeOffsetPoints[LAL_NUM_IFO];
  gsl_matrix *rotReOverlaps;
  gsl_matrix *rotImOverlaps;

  REAL4 *SNRu1,*SNRu2;
  REAL4 BankVetoTemp[4];
  REAL4 BankVeto=0;
  REAL4 normFac;
  REAL4 *TjwithS1,*TjwithS2;
  SNRu1 = LALCalloc(2,sizeof(REAL4));
  SNRu2 = LALCalloc(2,sizeof(REAL4));
  TjwithS1 = LALCalloc(2,sizeof(REAL4));
  TjwithS2 = LALCalloc(2,sizeof(REAL4));

  /* Ensure that the time sliding doesn't push into data points that don't exist */
  for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO ; ifoNumber++ )
  {
    if ( params->haveTrig[ifoNumber] )
    {
      if ( position+timeOffsetPoints[ifoNumber] >= 3*numPoints/4 +5000)
      {
        calTimeOffsetPoints[ifoNumber] = 3*numPoints/4 +4999 - position;
        fprintf(stderr,"Overflow occured in time shifting in bank veto\n");
      }
      else if  ( position+timeOffsetPoints[ifoNumber] < numPoints/4 - 5000)
      {
        calTimeOffsetPoints[ifoNumber] = numPoints/4 - 4999 - position;
        fprintf(stderr,"Overflow occured in time shifting in bank veto\n");
      }
      else
        calTimeOffsetPoints[ifoNumber]=timeOffsetPoints[ifoNumber];
    }
  }

  /* Begin by calculating the components of the SNR */
  calculate_rotated_vectors(params,PTFqVec,SNRu1,SNRu2,a,b,
      timeOffsetPoints,Bankeigenvecs[subBankSize],Bankeigenvals[subBankSize],
      numPoints,position,1,2);
  
  /* The normalization factors are already calculated, they are the eigenvalues*/
  halfNumPoints = 3*numPoints/4 - numPoints/4 + 10000;

  for ( ui = 0 ; ui < subBankSize ; ui++ )
  {
    rotReOverlaps = cohBankOverlaps[ui].rotReOverlaps;
    rotImOverlaps = cohBankOverlaps[ui].rotImOverlaps;
    /* Calculate the components of subBank template with data */

    calculate_rotated_vectors(params,dataOverlaps[ui].PTFqVec,TjwithS1,
        TjwithS2,a,b,calTimeOffsetPoints,Bankeigenvecs[ui],
        Bankeigenvals[ui],halfNumPoints,position-numPoints/4+5000,1,2);
    for (uj = 0; uj < 4; uj++)
    {
      normFac = 0;
      if (uj < 2)
        BankVetoTemp[uj] = TjwithS1[uj];
      else
        BankVetoTemp[uj] = TjwithS2[uj-2];
      for (uk = 0; uk < 4; uk++)
      {
        if (uj < 2 && uk < 2)
        {
          BankVetoTemp[uj] -= gsl_matrix_get(rotReOverlaps,uk,uj)*SNRu1[uk];
          normFac += pow(gsl_matrix_get(rotReOverlaps,uk,uj),2);
        }
        if (uj < 2 && uk > 1 )
        {
          BankVetoTemp[uj] -= gsl_matrix_get(rotImOverlaps,uk-2,uj)*SNRu2[uk-2];
          normFac += pow(gsl_matrix_get(rotImOverlaps,uk-2,uj),2);
        }
        if (uj > 1 && uk < 2 )
        {
          BankVetoTemp[uj] += gsl_matrix_get(rotImOverlaps,uk,uj-2)*SNRu1[uk];
          normFac += pow(gsl_matrix_get(rotImOverlaps,uk,uj-2),2);
        }
        if (uj > 1 && uk > 1 )
        {
          BankVetoTemp[uj]-=gsl_matrix_get(rotReOverlaps,uk-2,uj-2)*SNRu2[uk-2];
          normFac += pow(gsl_matrix_get(rotReOverlaps,uk-2,uj-2),2);
        }
      }
      BankVeto+=pow(BankVetoTemp[uj],2)/(1-normFac);
     
//      fprintf(stderr,"BV comps: %e %e\n",BankVetoTemp[uj],normFac);
//      fprintf(stderr,"Bank Veto: %e \n",BankVeto);
    }
//    fprintf(stderr,"SNR comps: %e %e %e %e \n",SNRu1[0],SNRu1[1],SNRu2[0],SNRu2[1]);
//    fprintf(stderr,"Bank SNR comps: %e %e %e %e \n",TjwithS1[0],TjwithS1[1],TjwithS2[0],TjwithS2[1]);
//    fprintf(stderr,"Overlaps: %e %e %e %e %e %e %e %e",gsl_matrix_get(rotReOverlaps,0,0),gsl_matrix_get(rotReOverlaps,0,1),gsl_matrix_get(rotReOverlaps,1,0),gsl_matrix_get(rotReOverlaps,1,1),gsl_matrix_get(rotImOverlaps,0,0),gsl_matrix_get(rotImOverlaps,0,1),gsl_matrix_get(rotImOverlaps,1,0),gsl_matrix_get(rotImOverlaps,1,1));
  }

  LALFree(SNRu1);
  LALFree(SNRu2);
  LALFree(TjwithS1);
  LALFree(TjwithS2);

  return BankVeto;
}
 
REAL4 calculate_auto_veto_max_phase_coherent(
UINT4           numPoints,
UINT4           position,
REAL4           a[LAL_NUM_IFO],
REAL4           b[LAL_NUM_IFO],
struct coh_PTF_params      *params,
struct bankCohTemplateOverlaps *cohAutoOverlaps,
COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO+1],
INT4            timeOffsetPoints[LAL_NUM_IFO],
gsl_matrix *Autoeigenvecs,
gsl_vector *Autoeigenvals
)
{
  UINT4 ui,uj,uk;
  gsl_matrix *rotReOverlaps;
  gsl_matrix *rotImOverlaps;
  UINT4 timeStepPoints = 0;
  timeStepPoints = params->autoVetoTimeStep*params->sampleRate;

  REAL4 *SNRu1,*SNRu2;
  REAL4 AutoVetoTemp[4];
  REAL4 AutoVeto=0;
  REAL4 normFac;
  REAL4 *TjwithS1,*TjwithS2;
  SNRu1 = LALCalloc(4,sizeof(REAL4));
  SNRu2 = LALCalloc(4,sizeof(REAL4));
  TjwithS1 = LALCalloc(4,sizeof(REAL4));
  TjwithS2 = LALCalloc(4,sizeof(REAL4));

  /* Begin by calculating the components of the SNR */
  calculate_rotated_vectors(params,PTFqVec,SNRu1,SNRu2,a,b,
      timeOffsetPoints,Autoeigenvecs,Autoeigenvals,
      numPoints,position,1,2);

  for ( ui = 0 ; ui < params->numAutoPoints ; ui++ )
  {
    rotReOverlaps = cohAutoOverlaps[ui].rotReOverlaps;
    rotImOverlaps = cohAutoOverlaps[ui].rotImOverlaps;
    /* Calculate the components of subBank template with data */

    calculate_rotated_vectors(params,PTFqVec,TjwithS1,
        TjwithS2,a,b,timeOffsetPoints,Autoeigenvecs,
        Autoeigenvals,numPoints,position-((ui+1) * timeStepPoints),1,2);
    for (uj = 0; uj < 4; uj++)
    {
      normFac = 0;
      if (uj < 2)
        AutoVetoTemp[uj] = TjwithS1[uj];
      else
        AutoVetoTemp[uj] = TjwithS2[uj-2];
//      fprintf(stderr,"Initial %e \n",AutoVetoTemp[uj]);
      for (uk = 0; uk < 4; uk++)
      {
        if (uj < 2 && uk < 2)
        {
          AutoVetoTemp[uj] -= gsl_matrix_get(rotReOverlaps,uk,uj)*SNRu1[uk];
          normFac += pow(gsl_matrix_get(rotReOverlaps,uk,uj),2);
//          fprintf(stderr,"1 %e %e \n",gsl_matrix_get(rotReOverlaps,uk,uj),SNRu1[uk]);
        }
        if (uj < 2 && uk > 1 )
        {
          AutoVetoTemp[uj] += gsl_matrix_get(rotImOverlaps,uk-2,uj)*SNRu2[uk-2];
          normFac += pow(gsl_matrix_get(rotImOverlaps,uk-2,uj),2);
//          fprintf(stderr,"2 %e %e \n",gsl_matrix_get(rotImOverlaps,uk-2,uj),SNRu2[uk-2]);
        }
        if (uj > 1 && uk < 2 )
        {
          AutoVetoTemp[uj] -= gsl_matrix_get(rotImOverlaps,uk,uj-2)*SNRu1[uk];
          normFac += pow(gsl_matrix_get(rotImOverlaps,uk,uj-2),2);
//          fprintf(stderr,"3 Plus %e %e \n",gsl_matrix_get(rotImOverlaps,uk,uj-2),SNRu1[uk]);
        }
        if (uj > 1 && uk > 1 )
        {
          AutoVetoTemp[uj]-=gsl_matrix_get(rotReOverlaps,uk-2,uj-2)*SNRu2[uk-2];
          normFac += pow(gsl_matrix_get(rotReOverlaps,uk-2,uj-2),2);
//          fprintf(stderr,"4 %e %e \n",gsl_matrix_get(rotReOverlaps,uk-2,uj-2),SNRu2[uk-2]);
        }
//        fprintf(stderr,"Auto veto comps: %e %e %e %e %e %e\n",SNRu1[0],SNRu2[0],SNRu1[1],SNRu2[1],AutoVetoTemp[uj],normFac);
//        fprintf(stderr,"Overlap comps: %e %e %e %e\n\n",gsl_matrix_get(rotReOverlaps,0,0),gsl_matrix_get(rotReOverlaps,1,1),gsl_matrix_get(rotImOverlaps,0,0),gsl_matrix_get(rotImOverlaps,1,1));
        
      }
      AutoVeto+=pow(AutoVetoTemp[uj],2)/(1-normFac);
    }
  }

  LALFree(SNRu1);
  LALFree(SNRu2);
  LALFree(TjwithS1);
  LALFree(TjwithS2);

  return AutoVeto;
}


void free_bank_veto_memory(
  struct bankTemplateOverlaps *bankNormOverlaps,
  InspiralTemplate        *PTFBankTemplates,
  FindChirpTemplate       *bankFcTmplts,
  UINT4 subBankSize,
  struct bankComplexTemplateOverlaps *bankOverlaps,
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


void calculate_coherent_bank_overlaps(
  struct coh_PTF_params   *params,
  struct bankComplexTemplateOverlaps bankOverlaps,
  struct bankCohTemplateOverlaps cohBankOverlaps,
  REAL4           a[LAL_NUM_IFO],
  REAL4           b[LAL_NUM_IFO],
  gsl_matrix *eigenvecs,
  gsl_vector *eigenvals,
  gsl_matrix *Bankeigenvecs,
  gsl_vector *Bankeigenvals
)
{
  UINT4 uk,uj;
  gsl_matrix *rotReOverlaps = cohBankOverlaps.rotReOverlaps;
  gsl_matrix *rotImOverlaps = cohBankOverlaps.rotImOverlaps;  

  gsl_matrix *reOverlaps = gsl_matrix_alloc(2,2);
  gsl_matrix *imOverlaps = gsl_matrix_alloc(2,2);
  gsl_matrix *tempM = gsl_matrix_alloc(2,2);
  REAL4 reOverlapsA[4],imOverlapsA[4];

  /* First step would be to create a matrix of overlaps in non rotated frame*/
  for (uk = 0 ; uk < 4; uk++)
  {
    reOverlapsA[uk] = 0;
    imOverlapsA[uk] = 0;
  }
  for (uk = 0; uk < LAL_NUM_IFO; uk++)
  {
    if ( params->haveTrig[uk] )
    {
      reOverlapsA[0] += a[uk] * a[uk] *  bankOverlaps.PTFM[uk]->data[0].re;
      reOverlapsA[1] += a[uk] * b[uk] *  bankOverlaps.PTFM[uk]->data[0].re;
      reOverlapsA[2] += b[uk] * a[uk] *  bankOverlaps.PTFM[uk]->data[0].re;
      reOverlapsA[3] += b[uk] * b[uk] *  bankOverlaps.PTFM[uk]->data[0].re;
      imOverlapsA[0] += a[uk] * a[uk] *  bankOverlaps.PTFM[uk]->data[0].im;
      imOverlapsA[1] += a[uk] * b[uk] *  bankOverlaps.PTFM[uk]->data[0].im;
      imOverlapsA[2] += b[uk] * a[uk] *  bankOverlaps.PTFM[uk]->data[0].im;
      imOverlapsA[3] += b[uk] * b[uk] *  bankOverlaps.PTFM[uk]->data[0].im;
    }
  }
  for (uj = 0;uj < 2;uj++)
  {
    for (uk=0; uk < 2; uk++)
    {
      gsl_matrix_set(reOverlaps,uj,uk,reOverlapsA[uj*2+uk]);
      gsl_matrix_set(imOverlaps,uj,uk,imOverlapsA[uj*2+uk]);
    }
  }
 
  /* And rotate by both set of eigenvectors */

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1,reOverlaps,eigenvecs,0.,tempM);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1,Bankeigenvecs,tempM,0.,rotReOverlaps);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1,imOverlaps,eigenvecs,0.,tempM);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1,Bankeigenvecs,tempM,0.,rotImOverlaps);

  for (uj=0; uj<2; uj++)
  {
    for (uk=0; uk<2 ;uk++ )
    {
      gsl_matrix_set(rotReOverlaps,uj,uk,gsl_matrix_get(rotReOverlaps,uj,uk)/(sqrt(gsl_vector_get(eigenvals,uj))*sqrt(gsl_vector_get(Bankeigenvals,uk))));
      gsl_matrix_set(rotImOverlaps,uj,uk,gsl_matrix_get(rotImOverlaps,uj,uk)/(sqrt(gsl_vector_get(eigenvals,uj))*sqrt(gsl_vector_get(Bankeigenvals,uk))));
    }
  }
//  fprintf(stderr,"Initial overlaps: %e %e\n",reOverlapsA[0],imOverlapsA[0]);
//  fprintf(stderr,"Rotated overlaps: %e %e %e %e\n",gsl_matrix_get(rotImOverlaps,0,0),gsl_matrix_get(rotImOverlaps,1,1),gsl_matrix_get(rotImOverlaps,1,0),gsl_matrix_get(rotImOverlaps,0,1));

  gsl_matrix_free(reOverlaps);
  gsl_matrix_free(imOverlaps);
  gsl_matrix_free(tempM);

}

void calculate_standard_chisq_freq_ranges(
    struct coh_PTF_params   *params,
    FindChirpTemplate       *fcTmplt,
    REAL4FrequencySeries    *invspec[LAL_NUM_IFO+1],
    REAL8Array              *PTFM[LAL_NUM_IFO+1],
    REAL4 a[LAL_NUM_IFO],
    REAL4 b[LAL_NUM_IFO],
    REAL4 *frequencyRangesPlus,
    REAL4 *frequencyRangesCross,
    gsl_matrix *eigenvecs
)
{
  UINT4 i,k,kmin,kmax,len,freqBinPlus,freqBinCross,numFreqBins;
  REAL4 v1,v2,v3,u1,u2,overlapCont,SNRtempPlus,SNRtempCross,SNRmaxPlus,SNRmaxCross;
  REAL8         f_min, deltaF, fFinal;
  COMPLEX8     *PTFQtilde   = NULL;
  REAL4 a2[LAL_NUM_IFO];
  REAL4 b2[LAL_NUM_IFO];
//  FILE *outfile;

  PTFQtilde = fcTmplt->PTFQtilde->data;  
  len = 0;
  deltaF = 0;
  for ( k = 0; k < LAL_NUM_IFO; k++)
  {
    if ( params->haveTrig[k] )
    {
      len       = invspec[k]->data->length;
      deltaF    = invspec[k]->deltaF;
      break;
    }
  }
  /* This is explicit as I want f_min of template lower than f_min of filter*/
  /* Note that these frequencies are not just hardcoded here, if you change*/
  /* these values you will need to change them in other places as well */
  f_min     = 40;
  kmin      = f_min / deltaF > 1 ?  f_min / deltaF : 1;
  fFinal    = 1000;
  kmax      = fFinal / deltaF < (len - 1) ? fFinal / deltaF : (len - 1);

  numFreqBins = params->numChiSquareBins;

  v1 = 0;
  v2 = 0;
  v3 = 0;
  for( k = 0; k < LAL_NUM_IFO; k++)
  {
    if ( params->haveTrig[k] )
    {
      a2[k] = a[k]*gsl_matrix_get(eigenvecs,0,0) + b[k]*gsl_matrix_get(eigenvecs,1,0);
      b2[k] = a[k]*gsl_matrix_get(eigenvecs,0,1) + b[k]*gsl_matrix_get(eigenvecs,1,1);
      v1 += a2[k]*a2[k]*PTFM[k]->data[0];
      v2 += b2[k]*b2[k]*PTFM[k]->data[0];
      v3 += a2[k]*b2[k]*PTFM[k]->data[0];
    }
  }
//  fprintf(stderr,"%e %e %e\n",v1,v2,v3);

  u1 = v1;
  u2 = v2;
  SNRmaxPlus = u1;
  if (SNRmaxPlus < 0) SNRmaxPlus = -SNRmaxPlus;
  SNRmaxCross = u2;
  if (SNRmaxCross < 0) SNRmaxCross = -SNRmaxCross;

//  fprintf(stderr,"%e %e \n",SNRmaxPlus,SNRmaxCross);

  v1 = 0;
  v2 = 0;

  freqBinPlus = 1;
  freqBinCross = 1;
  SNRtempPlus = 0;
  SNRtempCross = 0;
//  outfile = fopen("temp.dat","w");
  for ( i = kmin; i < kmax ; ++i )
  {
    for( k = 0; k < LAL_NUM_IFO; k++)
    {
      if ( params->haveTrig[k] )
      {
        overlapCont = (PTFQtilde[i].re * PTFQtilde[i].re +
               PTFQtilde[i].im * PTFQtilde[i].im )* invspec[k]->data->data[i] ;
        v1 += a2[k] * a2[k] * overlapCont * 4 * deltaF;
        v2 += b2[k] * b2[k] * overlapCont * 4 * deltaF;
      }
    }
    /* Calculate SNR */
//    SNRtemp = v1*v1 + v2*v2;
    u1 = v1;
    u2 = v2;
    SNRtempPlus = u1;
    if (SNRtempPlus < 0) SNRtempPlus = -SNRtempPlus;
    SNRtempCross = u2;
    if (SNRtempCross < 0) SNRtempCross = -SNRtempCross;
    /* Compare to max SNR */
//    fprintf(outfile,"%e %e %e %e\n",SNRtempPlus,SNRtempCross,v1,v2);
    if (SNRtempPlus > SNRmaxPlus * ((REAL4)freqBinPlus/(REAL4)numFreqBins))
    {
      if (freqBinPlus < numFreqBins)
      {
        /* Record the frequency */
        frequencyRangesPlus[freqBinPlus-1] = i*deltaF;
//        fprintf(stderr,"Frequency bin:%e \n",frequencyRangesPlus[freqBinPlus-1]);
        freqBinPlus+=1;
      }
    }
    if (SNRtempCross > SNRmaxCross * ((REAL4)freqBinCross/(REAL4)numFreqBins))
    {
      if (freqBinCross < numFreqBins)
      {
        /* Record the frequency */
        frequencyRangesCross[freqBinCross-1] = i*deltaF;
//        fprintf(stderr,"Frequency bin:%e \n",frequencyRangesCross[freqBinCross-1]);
        freqBinCross+=1;
      }
    }
  }
//  fclose(outfile);
//  fprintf(stderr,"%e %e \n", SNRtemp,SNRmax);
}

void calculate_standard_chisq_power_bins(
    struct coh_PTF_params   *params,
    FindChirpTemplate       *fcTmplt,
    REAL4FrequencySeries    *invspec[LAL_NUM_IFO+1],
    REAL8Array              *PTFM[LAL_NUM_IFO+1],
    REAL4 a[LAL_NUM_IFO],
    REAL4 b[LAL_NUM_IFO],
    REAL4 *frequencyRanges,
    REAL4 *powerBinsPlus,
    REAL4 *powerBinsCross,
    gsl_matrix *eigenvecs
)
{
  UINT4 i,k,kmin,kmax,len,freqBin,numFreqBins;
  REAL4 v1,v2,v3,overlapCont,SNRtempPlus,SNRtempCross,SNRmaxPlus,SNRmaxCross;
  REAL4 SNRplusLast,SNRcrossLast;
  REAL8         f_min, deltaF, fFinal;
  COMPLEX8     *PTFQtilde   = NULL;
  REAL4 a2[LAL_NUM_IFO];
  REAL4 b2[LAL_NUM_IFO];
//  FILE *outfile;

  PTFQtilde = fcTmplt->PTFQtilde->data;
  len = 0;
  deltaF = 0;
  for ( k = 0; k < LAL_NUM_IFO; k++)
  {
    if ( params->haveTrig[k] )
    {
      len       = invspec[k]->data->length;
      deltaF    = invspec[k]->deltaF;
      break;
    }
  }
  /* This is explicit as I want f_min of template lower than f_min of filter*/
  /* Note that these frequencies are not just hardcoded here, if you change*/
  /* these values you will need to change them in other places as well */
  f_min     = 40;
  kmin      = f_min / deltaF > 1 ?  f_min / deltaF : 1;
  fFinal    = 1000;
  kmax      = fFinal / deltaF < (len - 1) ? fFinal / deltaF : (len - 1);

  numFreqBins = params->numChiSquareBins;

  // NOTE: v3 is calculated for verification. It should = 0.
  v1 = 0;
  v2 = 0;
  v3 = 0;
  for( k = 0; k < LAL_NUM_IFO; k++)
  {
    if ( params->haveTrig[k] )
    {
      a2[k] = a[k]*gsl_matrix_get(eigenvecs,0,0) + b[k]*gsl_matrix_get(eigenvecs,1,0);
      b2[k] = a[k]*gsl_matrix_get(eigenvecs,0,1) + b[k]*gsl_matrix_get(eigenvecs,1,1);
      v1 += a2[k]*a2[k]*PTFM[k]->data[0];
      v2 += b2[k]*b2[k]*PTFM[k]->data[0];
      v3 += a2[k]*b2[k]*PTFM[k]->data[0];
    }
  }
//  fprintf(stderr,"%e %e %e\n",v1,v2,v3);

  SNRmaxPlus = v1;
  SNRmaxCross = v2;
  if (SNRmaxPlus < 0) SNRmaxPlus = -SNRmaxPlus;
  if (SNRmaxCross < 0) SNRmaxCross = -SNRmaxCross;

  v1 = 0;
  v2 = 0;

  freqBin = 0;
  SNRtempPlus = 0;
  SNRtempCross = 0;
  SNRplusLast = 0.;
  SNRcrossLast = 0.;

  for ( i = kmin; i < kmax ; ++i )
  {
    for( k = 0; k < LAL_NUM_IFO; k++)
    {
      if ( params->haveTrig[k] )
      {
        overlapCont = (PTFQtilde[i].re * PTFQtilde[i].re +
               PTFQtilde[i].im * PTFQtilde[i].im )* invspec[k]->data->data[i] ;
        v1 += a2[k] * a2[k] * overlapCont * 4 * deltaF;
        v2 += b2[k] * b2[k] * overlapCont * 4 * deltaF;
      }
    }
    SNRtempPlus = v1;
    if (SNRtempPlus < 0) SNRtempPlus = -SNRtempPlus;
    SNRtempCross = v2;
    if (SNRtempCross < 0) SNRtempCross = -SNRtempCross;

    if (i * deltaF > frequencyRanges[freqBin] && freqBin < (numFreqBins-1))
    {
      powerBinsPlus[freqBin] = SNRtempPlus/SNRmaxPlus - SNRplusLast;
      powerBinsCross[freqBin] = SNRtempCross/SNRmaxCross - SNRcrossLast;
//      fprintf(stderr,"Constructing plus: %e %e %e \n",SNRtempPlus,SNRmaxPlus,SNRplusLast);
//      fprintf(stderr,"Constructing cross: %e %e %e \n",SNRtempCross,SNRmaxCross,SNRcrossLast);
      SNRplusLast = SNRtempPlus/SNRmaxPlus;
      SNRcrossLast = SNRtempCross/SNRmaxCross;
      freqBin++;
    }
  }
  if (freqBin == (numFreqBins-1))  
  {
    powerBinsPlus[freqBin] = SNRtempPlus/SNRmaxPlus - SNRplusLast;
    powerBinsCross[freqBin] = SNRtempCross/SNRmaxCross - SNRcrossLast;
  }   

  /* Ensure that the power Bins add to 1. This should already be true but
  numerical counting errors can have an effect here.*/
  SNRplusLast = 0;
  SNRcrossLast = 0;
  for ( i = 0 ; i < (numFreqBins); i++)
  {
//    fprintf(stderr,"Power bins: %e %e \n",powerBinsPlus[i],powerBinsCross[i]);
    SNRplusLast += powerBinsPlus[i];
    SNRcrossLast += powerBinsCross[i];
//    fprintf(stderr,"Totals: %e %e \n", SNRplusLast,SNRcrossLast);
  }
  for ( i = 0 ; i < (numFreqBins); i++)
  {
    powerBinsPlus[freqBin] = powerBinsPlus[freqBin]/SNRplusLast;
    powerBinsCross[freqBin] = powerBinsCross[freqBin]/SNRcrossLast;
//    fprintf(stderr,"Power bins: %e %e \n",powerBinsPlus[i],powerBinsCross[i]);

  }
}

        

REAL4 calculate_chi_square(
struct coh_PTF_params   *params,
UINT4           numPoints,
UINT4           position,
struct bankDataOverlaps *chisqOverlaps,    
COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO+1],
REAL4           a[LAL_NUM_IFO],
REAL4           b[LAL_NUM_IFO],
INT4            timeOffsetPoints[LAL_NUM_IFO],
gsl_matrix *eigenvecs,
gsl_vector *eigenvals,   
REAL4 *powerBinsPlus,
REAL4 *powerBinsCross
)
{
  UINT4 i,halfNumPoints;
  REAL4 *v1,*v2,*v1full,*v2full;
  REAL4 chiSq,SNRtemp,SNRexp;
  UINT4 numChiSquareBins = params->numChiSquareBins;

  v1 = LALCalloc(2,sizeof(REAL4));
  v2 = LALCalloc(2,sizeof(REAL4));
  v1full = LALCalloc(2,sizeof(REAL4));
  v2full = LALCalloc(2,sizeof(REAL4));

  halfNumPoints = 3*numPoints/4 - numPoints/4 + 10000;

  chiSq = 0;

  calculate_rotated_vectors(params,PTFqVec,v1full,v2full,a,b,
        timeOffsetPoints,eigenvecs,eigenvals,numPoints,
        position,1,2);

  SNRexp = v1full[0]*v1full[0] + v1full[1]*v1full[1];
  SNRexp += v2full[0]*v2full[0] + v2full[1]*v2full[1];
  SNRexp = pow(SNRexp,0.5);

  for (i = 0; i < numChiSquareBins; i++ )
  {
    /* calculate SNR in this frequency bin */
    calculate_rotated_vectors(params,chisqOverlaps[i].PTFqVec,v1,v2,a,b,
        timeOffsetPoints,eigenvecs,eigenvals,halfNumPoints,
        position-numPoints/4+5000,1,2);

    SNRtemp = pow((v1[0] - v1full[0]*powerBinsPlus[i]),2)/powerBinsPlus[i];
//    fprintf(stderr,"Comp 1: %e %e \n",v1[0],v1full[0]);
    SNRtemp += pow((v1[1] - v1full[1]*powerBinsCross[i]),2)/powerBinsCross[i];
//    fprintf(stderr,"Comp 2: %e %e \n",v1[1],v1full[1]);
    SNRtemp += pow((v2[0] - v2full[0]*powerBinsPlus[i]),2)/powerBinsPlus[i];
//    fprintf(stderr,"Comp 3: %e %e \n",v2[0],v2full[0]);
    SNRtemp += pow((v2[1] - v2full[1]*powerBinsCross[i]),2)/powerBinsCross[i];
//    fprintf(stderr,"Comp 4: %e %e \n",v2[1],v2full[1]);
    chiSq += SNRtemp;
  }
  LALFree(v1);
  LALFree(v2);
  LALFree(v1full);
  LALFree(v2full);
  return chiSq;
}

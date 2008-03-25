/*
*  Copyright (C) 2007 Chad Hanna, Duncan Brown, Benjamin Owen, B.S. Sathyaprakash, Anand Sengupta, Thomas Cokelaer
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



#include<lal/LALStdlib.h>
#include<lal/LALStatusMacros.h>
#include<lal/LALInspiral.h>
#include<lal/LALInspiralBank.h>
#include<lal/LIGOMetadataTables.h>


NRCSID(INSPIRALBANKGENERATIONC, "$Id$");

void
LALInspiralBankGeneration(
     LALStatus *status,
     InspiralCoarseBankIn *input,
     SnglInspiralTable **first,
     INT4 *ntiles )
{
  InspiralTemplateList *coarseList = NULL;
  SnglInspiralTable *bank;
  InspiralMomentsEtc moments;
  INT4 cnt        = 0;
  REAL8 fFinal    = 0;
  REAL8 MinfFinal = 0;
  REAL8 MaxfFinal = 0;
  REAL8 q         = 0;
  INT4  chicnt    = 0;
  INT4  kappacnt  = 0;
  INT4  numTmplts = 0;
  INT4  i;
  REAL4 chi[3], kappa[4];
  
  INITSTATUS(status, "LALInspiralBankGeneration", INSPIRALBANKGENERATIONC);
  ATTATCHSTATUSPTR(status);
    
  ASSERT( input != NULL, status, LALINSPIRALBANKH_ENULL,
          LALINSPIRALBANKH_MSGENULL );
  ASSERT( *first == NULL, status, LALINSPIRALBANKH_ENULL,
          LALINSPIRALBANKH_MSGENULL );

  
  
  /* For nonspinning approximants, call LALInspiralCreateCoarseBank(). */
  switch( input->approximant )
  {
  case BCV:
  case EOB:
  case EOBNR:
  case PadeT1:
  case PadeF1:
  case TaylorF1:
  case TaylorF2:
  case TaylorT1:
  case TaylorT2:
  case TaylorT3:
  case AmpCorPPN:
  case Eccentricity:

    /* Use LALInspiralCreateCoarseBank(). */
    TRY( LALInspiralCreateCoarseBank( status->statusPtr, &coarseList, ntiles,
         *input ), status );
    /* */ 
    /* Convert output data structure. */
    bank = (SnglInspiralTable *) LALCalloc(1, sizeof(SnglInspiralTable));
    if (bank == NULL){
      ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
    }
    *first = bank;
    for( cnt = 0; cnt < *ntiles; cnt++ )
    {
      /* Set the Min and Max fFinals using the appropriate formula*/
      if( input->MaxFreqCut == SchwarzISCO )
	{
	  MaxfFinal = 1.0 / (6.0 * sqrt(6.0)*LAL_PI*coarseList[cnt].params.totalMass*LAL_MTSUN_SI);
	}
      else if( input->MaxFreqCut == BKLISCO )
	{
	  if( coarseList[cnt].params.mass1 > coarseList[cnt].params.mass2 )
	    {
	      q = coarseList[cnt].params.mass2 / coarseList[cnt].params.mass1;
	    }
	  else
	      q = coarseList[cnt].params.mass1 / coarseList[cnt].params.mass2;
	  MaxfFinal = 1.0 / (6.0 * sqrt(6.0)*LAL_PI*coarseList[cnt].params.totalMass*LAL_MTSUN_SI) * ( 1 + 2.8*q - 2.6*q*q + 0.8*q*q*q );
	}
      else if( input->MaxFreqCut == ERD )
	{
	  MaxfFinal = 1.07*0.5326/(2*LAL_PI*0.955*coarseList[cnt].params.totalMass*LAL_MTSUN_SI);
	}
      else
	ABORT( status, LALINSPIRALBANKH_EFCUT, LALINSPIRALBANKH_MSGEFCUT );

      if( input->MinFreqCut == SchwarzISCO )
	{
	  MinfFinal = 1.0 / (6.0 * sqrt(6.0)*LAL_PI*coarseList[cnt].params.totalMass*LAL_MTSUN_SI);
	}
      else if( input->MinFreqCut == BKLISCO )
	{
	  if( coarseList[cnt].params.mass1 > coarseList[cnt].params.mass2 )
	    {
	      q = coarseList[cnt].params.mass2 / coarseList[cnt].params.mass1;
	    }
	  else
	      q = coarseList[cnt].params.mass1 / coarseList[cnt].params.mass2;
	  MinfFinal = 1.0 / (6.0 * sqrt(6.0)*LAL_PI*coarseList[cnt].params.totalMass*LAL_MTSUN_SI) * ( 1 + 2.8*q - 2.6*q*q + 0.8*q*q*q );
	}
      else if( input->MinFreqCut == ERD )
	{
	  MinfFinal = 1.07*0.5326/(2*LAL_PI*0.955*coarseList[cnt].params.totalMass*LAL_MTSUN_SI);
	}
      else
	ABORT( status, LALINSPIRALBANKH_EFCUT, LALINSPIRALBANKH_MSGEFCUT );

      /* For 1 upper frequency cutoff, fill the bank as usual with the
       * specified fFinal (checked MinFreqCut = MaxFreqCut in tmpltbank.c)
       */
      if( input->NumFreqCut == 1 )
	{
	  bank = bank->next = (SnglInspiralTable *) LALCalloc( 1, sizeof(
             SnglInspiralTable ) );
	  if (bank == NULL)
	    {
	      ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
	    }
	  bank->mass1 = coarseList[cnt].params.mass1;
	  bank->mass2 = coarseList[cnt].params.mass2;
	  bank->mchirp = coarseList[cnt].params.chirpMass;
	  bank->mtotal = coarseList[cnt].params.totalMass;
	  bank->eta = coarseList[cnt].params.eta;
	  bank->tau0 = coarseList[cnt].params.t0;
	  bank->tau2 = coarseList[cnt].params.t2;
	  bank->tau3 = coarseList[cnt].params.t3;
	  bank->tau4 = coarseList[cnt].params.t4;
	  bank->tau5 = coarseList[cnt].params.t5;
	  bank->ttotal = coarseList[cnt].params.tC;
	  bank->psi0 = coarseList[cnt].params.psi0;
	  bank->psi3 = coarseList[cnt].params.psi3;
      
	  /* This calucation is only valid for the PN case. For EOB, we 
	   * should use the correct value of v (close to lightring). What 
	   * about the amplitude corrected one ? */

	  fFinal = MinfFinal;
	  if (fFinal > input->fUpper)
	    {
	      fFinal = input->fUpper;
	    }
	  coarseList[cnt].params.fFinal = fFinal;
    
	  /* Update the Gamma parameter if requested, using the proper cut-off 
	   * frequency */
	  if ( input->computeMoments )
	    {
	      coarseList[cnt].params.fCutoff = coarseList[cnt].params.fFinal;
	      LALGetInspiralMoments( status->statusPtr, &moments, &(input->shf), 
				     &(coarseList[cnt].params) );

	      LALInspiralComputeMetric(status->statusPtr, &(coarseList[cnt].metric),
				       &(coarseList[cnt].params), &moments);
	    }


	  bank->f_final = coarseList[cnt].params.fFinal;
	  bank->eta = coarseList[cnt].params.eta;
	  bank->beta = coarseList[cnt].params.beta;
      
	  /* Copy the 10 metric co-efficients ... */
	  memcpy (bank->Gamma, coarseList[cnt].metric.Gamma, 10*sizeof(REAL4));
	}

      /* If we have multiple frequency cutoffs, create duplicate
       * templates evenly incremented between MinfFinal and MaxfFinal
       */
      else for( i = 0; i < input->NumFreqCut; i++ )
      {
      bank = bank->next = (SnglInspiralTable *) LALCalloc( 1, sizeof(
             SnglInspiralTable ) );
      if (bank == NULL)
      {
        ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
      }
      bank->mass1 = coarseList[cnt].params.mass1;
      bank->mass2 = coarseList[cnt].params.mass2;
      bank->mchirp = coarseList[cnt].params.chirpMass;
      bank->mtotal = coarseList[cnt].params.totalMass;
      bank->eta = coarseList[cnt].params.eta;
      bank->tau0 = coarseList[cnt].params.t0;
      bank->tau2 = coarseList[cnt].params.t2;
      bank->tau3 = coarseList[cnt].params.t3;
      bank->tau4 = coarseList[cnt].params.t4;
      bank->tau5 = coarseList[cnt].params.t5;
      bank->ttotal = coarseList[cnt].params.tC;
      bank->psi0 = coarseList[cnt].params.psi0;
      bank->psi3 = coarseList[cnt].params.psi3;
      
      /* This calucation is only valid for the PN case. For EOB, we 
       * should use the correct value of v (close to lightring). What 
       * about the amplitude corrected one ? */

      fFinal = MinfFinal + i *
	(MaxfFinal - MinfFinal)/(input->NumFreqCut - 1);
      if (fFinal > input->fUpper)
      {
        fFinal = input->fUpper;
      }
      coarseList[cnt].params.fFinal = fFinal;
    
      /* Update the Gamma parameter if requested, using the proper cut-off 
       * frequency */
      if ( input->computeMoments )
      {
        coarseList[cnt].params.fCutoff = coarseList[cnt].params.fFinal;
        LALGetInspiralMoments( status->statusPtr, &moments, &(input->shf), 
            &(coarseList[cnt].params) );

        LALInspiralComputeMetric(status->statusPtr, &(coarseList[cnt].metric),
            &(coarseList[cnt].params), &moments);
      }


      bank->f_final = coarseList[cnt].params.fFinal;
      bank->eta = coarseList[cnt].params.eta;
      bank->beta = coarseList[cnt].params.beta;
      
      /* Copy the 10 metric co-efficients ... */
      memcpy (bank->Gamma, coarseList[cnt].metric.Gamma, 10*sizeof(REAL4));
      }
      
    }
    /* Free first template, which is blank. */
    bank = (*first)->next;
    LALFree( *first );
    *first = bank;
    /* free the coarse list returned by create coarse bank */
    LALFree( coarseList );
    break;
  
  case FindChirpPTF:
    
    for (i=0; i<5; i++)
    {
      if ( i < 3 ) chi[i]     = 0.1 + i * 0.4 ;
      if ( i < 2 ) kappa[i]   = -0.9 + i * 0.4 ;
      if ( i > 2 ) kappa[i-1] = 0.5 + (i-3) * 0.4;
    }
    
    /* Use LALInspiralCreateCoarseBank(). */
    TRY( LALInspiralCreateCoarseBank( status->statusPtr, &coarseList, ntiles,
         *input ), status );
 
    /* Convert output data structure. */
    bank = (SnglInspiralTable *) LALCalloc(1, sizeof(SnglInspiralTable));
    if (bank == NULL){
      ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
    }
    *first = bank;

    for ( chicnt = 0; chicnt < 3; chicnt++ )
    {
      for( kappacnt = 0; kappacnt < 4; kappacnt++ )
      {
        for( cnt = 0; cnt < *ntiles; cnt++ )
        {
          /* restrict the bank boundaries to the region of validity of PTF */
          if ( coarseList[cnt].params.mass1 < 6.0 || 
               coarseList[cnt].params.mass2 > 3.0 ) continue;
          bank = bank->next = (SnglInspiralTable *) LALCalloc( 1, sizeof(
                SnglInspiralTable ) );
          if (bank == NULL)
          {
            ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
          }
          numTmplts     = numTmplts + 1 ;
          bank->mass1   = coarseList[cnt].params.mass1;
          bank->mass2   = coarseList[cnt].params.mass2;
          bank->mchirp  = coarseList[cnt].params.chirpMass;
          bank->mtotal  = coarseList[cnt].params.totalMass;
          bank->eta     = coarseList[cnt].params.eta;
          bank->kappa   = kappa[kappacnt];
          bank->chi     = chi[chicnt];
          bank->tau0    = coarseList[cnt].params.t0;
          bank->tau2    = coarseList[cnt].params.t2;
          bank->tau3    = coarseList[cnt].params.t3;
          bank->tau4    = coarseList[cnt].params.t4;
          bank->tau5    = coarseList[cnt].params.t5;
          bank->ttotal  = coarseList[cnt].params.tC;
          bank->psi0    = coarseList[cnt].params.psi0;
          bank->psi3    = coarseList[cnt].params.psi3;
          bank->f_final = coarseList[cnt].params.fFinal;
          bank->eta     = coarseList[cnt].params.eta;
          bank->beta    = coarseList[cnt].params.beta;
          

          /* Copy the 10 metric co-efficients ... */
          memcpy (bank->Gamma, coarseList[cnt].metric.Gamma, 10*sizeof(REAL4));

        }
      }
    }

    /* Free first template, which is blank. */
    bank = (*first)->next;
    LALFree( *first );
    *first = bank;
    /* free the coarse list returned by create coarse bank */
    LALFree( coarseList );
    *ntiles = numTmplts;
    break;

  case BCVSpin:
    if (input->spinBank==0)
    {
    /* Use LALInspiralSpinBank(); no need to convert output. */
    TRY( LALInspiralSpinBank( status->statusPtr, first, ntiles, input ),
         status );   
    }
    else if (input->spinBank==1)
    {
    /* For extended bank use LALInspiralBCVSpinBank() */
    /*
    TRY( LALInspiralBCVSpinBank( status->statusPtr, first, ntiles, input ),
         status );   
    */
    }
    else if (input->spinBank==2)
    {
    /* For extended bank use LALInspiralBCVSpinBank() */
    
    /*
    TRY( LALInspiralBCVSpinRandomBank( status->statusPtr, first, ntiles, input ),
         status );   
     */
     
    }
    else
    {
      ABORT( status, LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE );
    }

    if (*ntiles < 1){       
      ABORT( status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL );
    }
    break;

  default:
    ABORT( status, LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE );

  }

  DETATCHSTATUSPTR(status);
  RETURN(status); 
}



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
  INT4 cnt = 0;
  
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
  case PadeT1:
  case PadeF1:
  case TaylorF1:
  case TaylorF2:
  case TaylorT1:
  case TaylorT2:
  case TaylorT3:

    /* Use LALInspiralCreateCoarseBank(). */
    TRY( LALInspiralCreateCoarseBank( status->statusPtr, &coarseList, ntiles,
         *input ), status );
 
    /* Convert output data structure. */
    bank = (SnglInspiralTable *) LALCalloc(1, sizeof(SnglInspiralTable));
    if (bank == NULL){
      ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
    }
    *first = bank;
    for( cnt = 0; cnt < *ntiles; cnt++ )
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
      bank->eta = coarseList[cnt].params.eta;
      bank->tau0 = coarseList[cnt].params.t0;
      bank->tau2 = coarseList[cnt].params.t2;
      bank->tau3 = coarseList[cnt].params.t3;
      bank->tau4 = coarseList[cnt].params.t4;
      bank->tau5 = coarseList[cnt].params.t5;
      bank->ttotal = coarseList[cnt].params.tC;
      bank->psi0 = coarseList[cnt].params.psi0;
      bank->psi3 = coarseList[cnt].params.psi3;
      bank->f_final = coarseList[cnt].params.fFinal;
      bank->eta = coarseList[cnt].params.eta;
      bank->beta = coarseList[cnt].params.beta;
    }
    /* Free first template, which is blank. */
    bank = (*first)->next;
    LALFree( *first );
    *first = bank;
    /* free the coarse list returned by create coarse bank */
    LALFree( coarseList );
    break;

  case BCVSpin:
    /* Use LALInspiralSpinBank(); no need to convert output. */
    TRY( LALInspiralSpinBank( status->statusPtr, first, ntiles, input ),
         status );   
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

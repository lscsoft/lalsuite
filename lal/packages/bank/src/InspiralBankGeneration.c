

#include<lal/LALStdlib.h>
#include<lal/LALStatusMacros.h>
#include<lal/LALInspiral.h>
#include<lal/LALInspiralBank.h>
#include<lal/LIGOMetadataTables.h>


NRCSID(INSPIRALBANKGENERATIONC, "$Id$");

static void CalculateEOBFFinal( LALStatus *status, InspiralTemplateList *input, REAL8 *fFinal );

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
  INT4 cnt = 0;
  REAL8 fFinal =0;
  
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
    /* */ 
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
      if (input->approximant == EOB)
      {
        TRY( CalculateEOBFFinal( status->statusPtr, 
              &coarseList[cnt], &fFinal ), status) ;
        printf("EOB fFinal is %e\n", fFinal );
      }
      else
      {
        fFinal = pow(sqrt(1./6.),3.)
        /(LAL_PI * (coarseList[cnt].params.mass1+coarseList[cnt].params.mass2) 
        * LAL_MTSUN_SI );
      }
      if (fFinal > input->fUpper)
      {
        fFinal = input->fUpper;
      }
      if (fFinal <= input->fLower)
      {
        ABORT( status, LALINSPIRALBANKH_EFHIGH, LALINSPIRALBANKH_MSGEFHIGH );
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
    /* Free first template, which is blank. */
    bank = (*first)->next;
    LALFree( *first );
    *first = bank;
    /* free the coarse list returned by create coarse bank */
    LALFree( coarseList );
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


static void CalculateEOBFFinal( LALStatus *status, InspiralTemplateList *input, REAL8 *fFinal )
{
  InspiralTemplate params;
  InspiralInit     paramsInit;
  REAL4Vector      *signal = NULL;
  REAL4            dt;


  INITSTATUS(status, "CalculateEOBFFinal", INSPIRALBANKGENERATIONC);
  ATTATCHSTATUSPTR(status);

  /* Fill the template with appropriate values */
  memcpy( &params, &(input->params), sizeof( InspiralTemplate ) );

  TRY ( LALInspiralInit( status->statusPtr, &params, &paramsInit ), status );

  dt    = 1./params.tSampling;
  printf( "tSampling = %e\n", params.tSampling );

  TRY( LALSCreateVector( status->statusPtr, &signal, paramsInit.nbins ), status );

  /* Assuming everything is okay, we create the waveform and get the fFinal */

  LALInspiralWave( status->statusPtr, signal, &params );
  BEGINFAIL( status )
  {
    LALSDestroyVector( status->statusPtr, &signal );
    ABORT( status, LALINSPIRALH_ENOWAVEFORM, LALINSPIRALH_MSGENOWAVEFORM );
  }
  ENDFAIL( status );

  *fFinal = params.fFinal;

  /* Free memory */
  LALSDestroyVector( status->statusPtr, &signal );
  CHECKSTATUSPTR( status );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}

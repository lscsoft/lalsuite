                                                                                                                                                
#include<lal/LALStdlib.h>
#include<lal/LALStatusMacros.h>
#include<lal/InspiralBankGeneration.h>   
#include<lal/LALInspiral.h>
#include<lal/LALInspiralBank.h>
#include<lal/LIGOMetadataTables.h>
                                                                                                                                         
NRCSID(INSPIRALBANKGENERATIONC, "$Id$");

void
LALInspiralBankGeneration(LALStatus *status,
                          TemplateBankType *type,
                          InspiralTmpltBankCInput *input,
                          SnglInspiralTable *first)
{
  InspiralCoarseBankIn CoarseIn;
  InspiralTemplateList *CoarseList;
  INT4 ntiles = 0;
  SnglInspiralTable *bank;
  INT4 cnt = 0;
  
  INITSTATUS(status, "LALInspiralBankGeneration", INSPIRALBANKGENERATIONC);
  ATTATCHSTATUSPTR(status);
    
  /* Check the inputs 
  if (type == NULL){
    ABORT(status, INSPIRALBANKGENERATIONH_ENULL, INSPIRALBANKGENERATIONH_MSGENULL);
    }
  if (input == NULL){
    ABORT(status, INSPIRALBANKGENERATIONH_ENULL, INSPIRALBANKGENERATIONH_MSGENULL);
    }
  if (first != NULL){
    ABORT(status, INSPIRALBANKGENERATIONH_ENONNULL, INSPIRALBANKGENERATIONH_MSGENONNULL);
    }
  */
  /* For InspiralSpinBank() */                               
  if (*type == 102){
    CoarseIn.mMin = (REAL8) input->mMin;
    CoarseIn.mMax = (REAL8) input->mMax;
    CoarseIn.MMax = (REAL8) input->MMax;
    CoarseIn.psi0Min = (REAL8) input->psi0Min;
    CoarseIn.psi0Max = (REAL8) input->psi0Max;
    CoarseIn.psi3Min = (REAL8) input->psi3Min;
    CoarseIn.psi3Max = (REAL8) input->psi3Max;
    CoarseIn.alpha = (REAL8) input->alpha;
    CoarseIn.numFcutTemplates = input->numFcutTemplates;
    CoarseIn.mmCoarse = (REAL8) input->mmCoarse;
    CoarseIn.fLower = (REAL8) input->fLower;
    CoarseIn.fUpper = (REAL8) input->fUpper;
    CoarseIn.order = input->order;
    CoarseIn.shf = input->shf;
    CoarseIn.tSampling = input->tSampling;
    CoarseIn.etamin = (REAL8) input->etamin;
    CoarseIn.LowGM = input->LowGM;
    CoarseIn.HighGM = input->HighGM;
    CoarseIn.approximant = input->approximant;
    CoarseIn.space = input->space;
    CoarseIn.massRange = MinMaxComponentMass;
    
    /* call LALInspiralSpinBank() */
    TRY(LALInspiralSpinBank(status->statusPtr, &CoarseList, &ntiles, CoarseIn), status);   
    
    if (!ntiles){       
      ABORT(status, INSPIRALBANKGENERATIONH_ENULL, INSPIRALBANKGENERATIONH_MSGENULL);
      }
 
    input->ntiles = ntiles;   
    first = (SnglInspiralTable *) LALCalloc(1, sizeof(SnglInspiralTable));
    bank = first;

    if (bank == NULL){
      ABORT(status, INSPIRALBANKGENERATIONH_EMEM, INSPIRALBANKGENERATIONH_MSGEMEM);
      }

    for(cnt = 0; cnt < input->ntiles; cnt++){
      bank->mass1 = CoarseList[cnt].params.mass1;
      bank->mass2 = CoarseList[cnt].params.mass2;
      bank->psi0 = CoarseList[cnt].params.psi0;
      bank->psi3 = CoarseList[cnt].params.psi3;
      bank->eta = CoarseList[cnt].params.eta;
      bank->mchirp = CoarseList[cnt].params.chirpMass; 
      /* CoarseList[cnt].params.beta = bank->beta  THIS NEEDS TO BE ADDED TO SnglInspiralTable */
      bank->next = (SnglInspiralTable *) LALCalloc(1, sizeof(SnglInspiralTable));
      bank = bank->next;
      if (bank == NULL){
        ABORT(status, INSPIRALBANKGENERATIONH_EMEM, INSPIRALBANKGENERATIONH_MSGEMEM);
        }
      }
    bank->next == NULL; 
  }
     
  DETATCHSTATUSPTR(status);
  RETURN(status); 
}
                                                                                                                                                


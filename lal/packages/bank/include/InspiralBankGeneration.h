

#include<lal/LALStatusMacros.h>
#include<lal/LALStdlib.h>
#include<lal/LALDatatypes.h>
#include<lal/LIGOMetadataTables.h>
#include<lal/TemplateBankGeneration.h>
#include<lal/LALInspiralBank.h>

#define INSPIRALBANKGENERATIONH_ENULL 1
#define INSPIRALBANKGENERATIONH_MSGENULL "Unexpected NULL pointer to an input type"
#define INSPIRALBANKGENERATIONH_ENONNULL 2
#define INSPIRALBANKGENERATIONH_MSGENONNULL "Argument should be passed as a NULL POINTER"
#define INSPIRALBANKGENERATIONH_EMEM 3
#define INSPIRALBANKGENERATIONH_MSGEMEM "Memory allocation error"
        

void
LALInspiralBankGeneration(
     LALStatus *,
     TemplateBankType *,
     InspiralTmpltBankCInput *,
     SnglInspiralTable *);
     /* LALInspiralBankGeneration(); */

void
LALInspiralSpinBankNew(
    LALStatus            *,
    InspiralTemplateList **,
    INT4                 *,
    InspiralCoarseBankIn
    );
                                                                                                                                                
                                                                                                                                                
void
LALInspiralSpinBankBoundary(
   LALStatus            *,
   NDTemplateBankInput  *,
   NDTemplateBankOutput *,
   INT2                 *
   );
                                                                                                                                                
                                                                                                                                                
void
LALInspiralSpinBankMetric(
   LALStatus            *,
   NDTemplateBankInput  *,
   REAL4Array           *
   ); /* LALInspiralSpinBankMetric() Prototype */
                                                                                                                                                


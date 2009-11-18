#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>

#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOMetadataTables.h>

#include <math.h>

#define SUBBANKSIZE 5
#define TEMPLATE_LENGTH 1000
#define FLOW 0.0
#define DELTAF 1.0
#define DYNRANGE 1.0


SnglInspiralTable * convertTemplateBankToInspiralTable(InspiralTemplate *bankHead,SnglInspiralTable *rowHead);

int writeCCmatToFile(COMPLEX8Vector *ccmat,UINT4 maxSubBankSize,FILE *fdOut);

/* Program to test the essential bank chisquare functions.
 *
 * 1) This program tests the bank sorting algorithm XLALFindChirpSortTemplates,
 *  printing the sorted template bank out to file.
 * 2) This program tests the subbank creation algorithm XLALFindChirpCreateSubBanks,
 *  printing a few of the subbanks out to file.
 * 3) This program tests the cross-correlation matrix calculation in
 *  XLALBankVetoCCMat by computing the CC matrix for the first subbank.
 *
 */
int main(int argc, char ** argv)
{
    /* Input params to XLALBankVetoCCMat */
    FindChirpBankVetoData *bankVetoData = 
      (FindChirpBankVetoData *) calloc(1,sizeof(FindChirpBankVetoData));
    UINT4 subBankSize = SUBBANKSIZE;
    REAL4Vector *ampVec = XLALCreateREAL4Vector(TEMPLATE_LENGTH);
    REAL4 dynRange = DYNRANGE;
    REAL4 fLow = FLOW;
    REAL4 deltaF = DELTAF;
    REAL4 deltaT = DELTAF;

    /* Params needs for main */
    UINT4 sampleIndex;
    UINT4 templateIndex;
    REAL8 phase;

    /* initiate ccmat */
    bankVetoData->ccMat = XLALCreateCOMPLEX8Vector(subBankSize*subBankSize);
    bankVetoData->normMat = XLALCreateREAL4Vector(subBankSize);

    /* Create unity response function, spectrum and ampVec */
    bankVetoData->resp = XLALCreateCOMPLEX8Vector(TEMPLATE_LENGTH);
    bankVetoData->spec = XLALCreateREAL4Vector(TEMPLATE_LENGTH);
    for (sampleIndex = 0; sampleIndex < TEMPLATE_LENGTH; sampleIndex++)
    {
	bankVetoData->resp->data[sampleIndex].re = 1;
	bankVetoData->resp->data[sampleIndex].im = 0;
	bankVetoData->spec->data[sampleIndex] = 1;
	ampVec->data[sampleIndex] = 1;
    }


    /* Fill in templates */
    bankVetoData->length = subBankSize;
    bankVetoData->fcInputArray = 
      (FindChirpFilterInput **) calloc(subBankSize,sizeof(FindChirpFilterInput *));
    for ( templateIndex = 0; templateIndex < subBankSize; templateIndex++)
    {
	bankVetoData->fcInputArray[templateIndex] = 
	  (FindChirpFilterInput *) calloc(1,sizeof(FindChirpFilterInput));
	bankVetoData->fcInputArray[templateIndex]->fcTmplt = 
	  (FindChirpTemplate *) calloc(1,sizeof(FindChirpTemplate));
	bankVetoData->fcInputArray[templateIndex]->fcTmplt->data =
	  XLALCreateCOMPLEX8Vector(TEMPLATE_LENGTH);

	for (sampleIndex = 0; sampleIndex < TEMPLATE_LENGTH; sampleIndex++)
	{
	    phase = 2*LAL_PI*( (REAL8) (templateIndex*sampleIndex) ) / ( (REAL8) TEMPLATE_LENGTH );
	    bankVetoData->fcInputArray[templateIndex]->fcTmplt->data->data[sampleIndex].re =
	      cos(phase);
	    bankVetoData->fcInputArray[templateIndex]->fcTmplt->data->data[sampleIndex].im =
	      sin(phase);
	}

	bankVetoData->fcInputArray[templateIndex]->fcTmplt->tmplt.fFinal = 
	  (REAL8) (DELTAF*TEMPLATE_LENGTH-FLOW);
    }


    UINT4 trial;
    CHAR ccFileName[100];
    FILE *fdOut;
    bankVetoData->timeshift = XLALCreateREAL4Vector(subBankSize);
    for (trial = 0; trial < subBankSize; trial++)
    {
	/* set the time shift for all but one template to be 1 unit */
	memset(bankVetoData->timeshift->data, 0, subBankSize * sizeof(REAL4));
	bankVetoData->timeshift->data[trial] = 1.0/ (REAL4) TEMPLATE_LENGTH;

	/* Pass to CCMat function */
	XLALBankVetoCCMat (bankVetoData,
			   ampVec,
			   subBankSize,
			   dynRange,
			   fLow,
			   deltaF,
			   deltaT);


	/* write out results */
	sprintf(ccFileName,"ccmat_%d.txt",trial);
	fdOut = fopen(ccFileName,"w");
	writeCCmatToFile(bankVetoData->ccMat,bankVetoData->length,fdOut);
	fclose(fdOut);
    }

    return 0;

}


SnglInspiralTable *convertTemplateBankToInspiralTable(InspiralTemplate *bankHead,
				       SnglInspiralTable *rowHead)
{
    // Return with error if bank is empty
    if (!bankHead)
      return NULL;

    // Allocate space for inspiral table if necessary
    if (!rowHead)
      rowHead = (SnglInspiralTable *) calloc(1,sizeof(SnglInspiralTable));

    // Make copies of the pointers
    InspiralTemplate *bank = bankHead;
    SnglInspiralTable *snglInspiralRow  = rowHead;

    // Write data in bank to inspiral rows
    for ( bank = bankHead; bank; bank = bank->next)
    {

	snglInspiralRow->mass1 = bank->mass1;
	snglInspiralRow->mass2 = bank->mass2;
	snglInspiralRow->snr = 3.7;
	snglInspiralRow->chisq = 3.8;
	snglInspiralRow->chisq_dof = 100;
	snglInspiralRow->bank_chisq = 3.9;
	snglInspiralRow->bank_chisq_dof = 100;
	snglInspiralRow->cont_chisq = 4.0;
	snglInspiralRow->cont_chisq_dof = 100;

	if (!snglInspiralRow->next) 
	  snglInspiralRow->next = (SnglInspiralTable *) calloc(1,sizeof(SnglInspiralTable));

	snglInspiralRow = snglInspiralRow->next;
    
    }

    // Null terminate the table (to tell where the table ends)
    snglInspiralRow->next = NULL;




    return 0;

}



int writeCCmatToFile(COMPLEX8Vector *ccmat,UINT4 maxSubBankSize,FILE *fdOut)
{
    int row;
    int col;

    // print a matrix of the real parts
        for ( row = 0; row< maxSubBankSize; row++ )
    {
	for ( col = 0; col < maxSubBankSize; col++ )
	{
	    fprintf(fdOut,"%.4g\t",ccmat->data[row*maxSubBankSize+col].re);
	}
	fprintf(fdOut,"\n");
    }

    // print a matrix of the imaginary parts
    for ( row = 0; row< maxSubBankSize; row++ )
    {
	for ( col = 0; col < maxSubBankSize; col++ )
	{
	    fprintf(fdOut,"%.4g\t",ccmat->data[row*maxSubBankSize+col].im);
	}
	fprintf(fdOut,"\n");
    }




    return 0;
}

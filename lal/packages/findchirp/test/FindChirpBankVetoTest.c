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

#define SUBBANKSIZE 16
#define STARTTEMPLATE -1
#define ENDTEMPLATE -1
#define FLOW 40
#define DELTAF 0.04
#define DYNRANGE 1
#define NUMBER_OF_SUBBANKS_TO_WRITE_TO_FILE 10


FindChirpDataParams * initCCmatParams(FindChirpDataParams *params);

SnglInspiralTable * convertTemplateBankToInspiralTable(InspiralTemplate *bankHead,SnglInspiralTable *rowHead);


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
    /* Input and output file names */
    CHAR tmpltBankFileName[100];
    CHAR sortedBankFileName[100];
    CHAR subBankFileName[100];
    CHAR ccOutFileName[100];

    /* Construct needed for writing xml output */
    LIGOLwXMLStream *xmlStreamOut;

    /* Input to sub bank functions */
    InspiralTemplate *bankHead = NULL;
    FindChirpSubBank *subBankHead = NULL;
    FindChirpSubBank *subBank = NULL;
    FindChirpBankVetoData *bankVetoData = NULL;
    SnglInspiralTable *snglInspiralHead  = NULL;
    FindChirpDataParams *params = NULL;

    UINT4 numberTemplates;
    UINT4 subBankSize = SUBBANKSIZE;
    UINT4 maxSubBankSize = SUBBANKSIZE;
    REAL4 fLow = FLOW;
    REAL4 deltaF = DELTAF;
    REAL4 dynRange = DYNRANGE;

    UINT4 subbankIndex=0;

    sprintf(tmpltBankFileName,"%s","in.xml");
    sprintf(sortedBankFileName,"%s","sorted.xml");


    /* Read in the template bank from file */
    numberTemplates = InspiralTmpltBankFromLIGOLw( &bankHead, tmpltBankFileName,
						   STARTTEMPLATE,ENDTEMPLATE);

    /* Sort templates, convert to inspiral table and 
       write out table to file */
    bankHead = XLALFindChirpSortTemplates( bankHead, numberTemplates );
    
    snglInspiralHead = convertTemplateBankToInspiralTable( bankHead, snglInspiralHead);
   
    xmlStreamOut = XLALOpenLIGOLwXMLFile(sortedBankFileName);
    XLALWriteLIGOLwXMLSnglInspiralTable(xmlStreamOut,snglInspiralHead);
    XLALCloseLIGOLwXMLFile(xmlStreamOut);    

    /* Create subbanks from template bank and write a few of the
       subbanks out to file */
    subBankHead = XLALFindChirpCreateSubBanks(&maxSubBankSize,
					      subBankSize,
					      numberTemplates,
					      bankHead);

    
    /* Add the minimal information needed for CC mat to FCparams and bankVetoData */
    params = initCCmatParams(params);

    /* Write a few subbanks to file and compute a few CC mats */
    for (subBank=subBankHead; subBank; subBank = subBank->next)
    {
	if (subbankIndex == NUMBER_OF_SUBBANKS_TO_WRITE_TO_FILE) 
	  break;
	subbankIndex++;

	convertTemplateBankToInspiralTable(subBankHead->bankHead,snglInspiralHead);

	sprintf(subBankFileName,"subbank%d.xml",subbankIndex);
	xmlStreamOut = XLALOpenLIGOLwXMLFile(subBankFileName);
	XLALWriteLIGOLwXMLSnglInspiralTable(xmlStreamOut,snglInspiralHead);
	XLALCloseLIGOLwXMLFile(xmlStreamOut);	


	//convertSubbankTemplateToFreqSeries(subbank,bankVetoData);

	//XLALBankVetoCCMat(bankVetoData,subBank,params,
	//		  dynRange,fLow, deltaF);
	
	//sprintf(ccOutFileName,"subBank%d-ccOut.txt",subbankIndex);
	
	//writeCCMatToFile(ccOutFileName,bankVetoData);


    }





    

    
    /* Generate templates from first bank and compute CC */
    //    bankVetoData->length = maxSubBankSize;
    // bankVetoData->fcInputArray = (FindChirpFilterInput **)
    //calloc( bankVetoData->length, sizeof(FindChirpFilterInput*) );

    /*
    LALStatus status;
    FindChirpTmpltParams *tmpltParams = NULL;
    FindChirpInitParams *initParams;
    initParams->approximant = FindChirpSP;
    initParams->numPoints = 1000;

    tmpltParams->deltaT = 1;
    tmpltParams->fLow = fLow;
    tmpltParams->reverseChirpBank = 0;
    tmpltParams->order = LAL_PNORDER_THREE_POINT_FIVE;
    tmpltParams->bandPassTmplt = 0;
    
    int j;
    int k;
    for (j=0, bank = subBankHead->bankHead; bank; bank = bank->next,j++)
    {

      /* initialize the template functions 
      LALFindChirpTemplateInit( &status, &tmpltParams,initParams );

      LALFindChirpSPTemplate( &status, bankVetoData->fcInputArray[j]->fcTmplt,
			      bank,tmpltParams );

      for (k=0; (UINT4) k < bankVetoData->fcInputArray[j]->fcTmplt->data->length; k++)
	fprintf(stderr,"%g+i%g\n",
		bankVetoData->fcInputArray[j]->fcTmplt->data->data[k].re,
		bankVetoData->fcInputArray[j]->fcTmplt->data->data[k].im);

    }

    */


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

    return rowHead;

}



FindChirpDataParams * initCCmatParams(FindChirpDataParams *params)
{

    FindChirpInitParams *initParams = (FindChirpInitParams *) calloc(1,sizeof(FindChirpInitParams));


    // FIXME status pointer not working correctly
    LALStatus *status = (LALStatus *) calloc(1,sizeof(LALStatus));

    initParams->numSegments = 15;
    initParams->numPoints = 1048576;
    initParams->ovrlap = 524288;
    initParams->numChisqBins = 0;
    initParams->createRhosqVec = 0;
    initParams->createCVec = 0;
    initParams->approximant = FindChirpSP;
    initParams->order = LAL_PNORDER_THREE;

    if (!params)
      params = (FindChirpDataParams *) calloc(1,sizeof(FindChirpDataParams));

    LALFindChirpDataInit( status, &params, initParams );

    fprintf(stderr,"%d\n",params->fLow);

    fprintf(stderr,"in\n");    
    return params;

}

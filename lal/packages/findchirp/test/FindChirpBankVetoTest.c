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

#define SUBBANKSIZE 10
#define TEMPLATE_LENGTH 10000
#define FLOW 0.0
#define DELTAT 1.0
#define DYNRANGE 1.0


static void initBankVetoData( FindChirpBankVetoData *bankVetoData,REAL4Vector *ampVec,UINT4 subBankSize,UINT4 templateLength);
static void makeOrthogonalTemplates(FindChirpBankVetoData *bankVetoData,UINT4 subBankSize,UINT4 templateLength,REAL4 deltaF);
static void timeshiftTemplates(REAL4Vector * timeshift,UINT4 trial);
static int writeCCmatToFile(COMPLEX8Vector *ccmat,UINT4 subBankSize,CHAR *ccFileName);

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
    REAL4Vector *ampVec = XLALCreateREAL4Vector(TEMPLATE_LENGTH);
    REAL4 dynRange = DYNRANGE;
    REAL4 fLow = FLOW;
    REAL4 deltaT = DELTAT;
    UINT4 subBankSize = SUBBANKSIZE;
    REAL4 templateLength = TEMPLATE_LENGTH;
    REAL4 deltaF = 1.0/(deltaT*templateLength);

    /* Params needs for main */
    UINT4 trial;
    CHAR ccFileName[100];

    /* Initiate/fill basic structures needed by XLALBankVetoCCMat */
    initBankVetoData(bankVetoData,ampVec,subBankSize,templateLength);

    /* Create mock templates which are orthogonal */
    makeOrthogonalTemplates(bankVetoData,subBankSize,templateLength,deltaF);

    /* 
     * Make repeated calls to XLALBankVetoCCMat with time shifts 
     * in some of the templates. 
     */
    for (trial = 0; trial < subBankSize; trial++)
    {
	/* set the time shift for all but one template to be 1 unit */
	timeshiftTemplates(bankVetoData->timeshift,trial);

	/* Pass to CCMat function */
	XLALBankVetoCCMat (bankVetoData, ampVec, subBankSize,
			   dynRange, fLow, deltaF, deltaT);

	/* write out results */
	sprintf(ccFileName,"ccmat_%d.txt",trial);
	writeCCmatToFile(bankVetoData->ccMat,subBankSize,ccFileName);
    }

    return 0;

}



static void initBankVetoData( FindChirpBankVetoData *bankVetoData,REAL4Vector *ampVec,UINT4 subBankSize,UINT4 templateLength)
{
    UINT4 sampleIndex;

    /* initiate ccmat */
    bankVetoData->ccMat = XLALCreateCOMPLEX8Vector(subBankSize*subBankSize);
    
    /* initiate normMat */
    bankVetoData->normMat = XLALCreateREAL4Vector(subBankSize);

    /* Create unity response function, spectrum and ampVec */
    bankVetoData->resp = XLALCreateCOMPLEX8Vector(templateLength);
    bankVetoData->spec = XLALCreateREAL4Vector(templateLength);
    for (sampleIndex = 0; sampleIndex < templateLength; sampleIndex++)
    {
	bankVetoData->resp->data[sampleIndex].re = 1;
	bankVetoData->resp->data[sampleIndex].im = 0;
	bankVetoData->spec->data[sampleIndex] = 1;
	ampVec->data[sampleIndex] = 1;
    }

    /* initiate timeshift */
    bankVetoData->timeshift = XLALCreateREAL4Vector(subBankSize);

    return;

}


static void makeOrthogonalTemplates(FindChirpBankVetoData *bankVetoData,UINT4 subBankSize,UINT4 templateLength,REAL4 deltaF)
{
    UINT4 templateIndex;
    UINT4 sampleIndex;
    REAL4 phase;

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
	  XLALCreateCOMPLEX8Vector(templateLength);

	for (sampleIndex = 0; sampleIndex < templateLength; sampleIndex++)
	{
	    phase = 2*LAL_PI*( (REAL4) (templateIndex*sampleIndex) ) / ( (REAL4) templateLength );
	    bankVetoData->fcInputArray[templateIndex]->fcTmplt->data->data[sampleIndex].re =
	      cos(phase);
	    bankVetoData->fcInputArray[templateIndex]->fcTmplt->data->data[sampleIndex].im =
	      sin(phase);
	}
	
	bankVetoData->fcInputArray[templateIndex]->fcTmplt->tmplt.fFinal = 
	  (REAL4) (deltaF*templateLength-FLOW);
	
    }

    return;

}


static void timeshiftTemplates(REAL4Vector * timeshift,UINT4 trial)
{
    /* shift only one of the templates by only one time unit */
    memset(timeshift->data, 0, (timeshift->length) * sizeof(REAL4));
    timeshift->data[trial] = 2.0*DELTAT;

    return;

}



static int writeCCmatToFile(COMPLEX8Vector *ccmat,UINT4 subBankSize,CHAR *ccFileName)
{
    FILE *fdOut = fopen(ccFileName,"w");
    int row;
    int col;

    // print a matrix of the real parts
    fprintf(fdOut,"#Real part:\n");
    for ( row = 0; row< subBankSize; row++ )
    {
	for ( col = 0; col < subBankSize; col++ )
	{
	    fprintf(fdOut,"%.8f\t",ccmat->data[row*subBankSize+col].re);
	}
	fprintf(fdOut,"\n");
    }

    // print a matrix of the imaginary parts
    fprintf(fdOut,"#Imaginary part:\n");
    for ( row = 0; row< subBankSize; row++ )
    {
	for ( col = 0; col < subBankSize; col++ )
	{
	    fprintf(fdOut,"%.8f\t",ccmat->data[row*subBankSize+col].im);
	}
	fprintf(fdOut,"\n");
    }

    fclose(fdOut);

    return 0;
}

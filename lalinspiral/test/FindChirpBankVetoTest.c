#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>

#include <math.h>

/* TEST FUNCTION WILL (LIKELY) BREAK IF THESE PARAMETERS ARE CHANGED*/
#define FLOW 0.0
#define TOLERANCE 0.000001

/* THESE PARAMETERS CAN TAKE ANY VALUE */
#define SUBBANKSIZE 15
#define TEMPLATE_LENGTH 1024
#define DELTAT 0.25
#define DYNRANGE 2.0
#define SHIFTUNITS 0.0

static void initBankVetoData( FindChirpBankVetoData *bankVetoData,REAL4Vector *ampVec,UINT4 subBankSize,UINT4 templateLength);
static void makeDeltaFunctionTemplates(FindChirpBankVetoData *bankVetoData,UINT4 subBankSize,UINT4 templateLength,REAL4 deltaF);
static void timeshiftTemplates(REAL4Vector *timeshift,REAL4 deltaT,UINT4 trial);
static int checkCCmatForErrors(COMPLEX8Vector *ccMat,REAL4Vector *timeshift, REAL4 deltaT);
static int writeCCmatToFile(COMPLEX8Vector *ccmat,UINT4 subBankSize,CHAR *ccFileName);

/*
 * This program is intended to check the computation of the cross-correlation
 * matrix in XLALBankVetoCCMat.  The program generates time domain delta functions
 * and computes the CC matrix for these templates using XLALBankVetoCCMat.
 * The program then inspects the output to make sure that it agrees with the
 * expected result.  Any matrix that fails the test is printed to file.  The program
 * exits with exit value equal to the number of matrices that failed the test.
 *
 * NOTE: This is really useful for testing time shifting in the bank veto code
 * which is currently disabled.  However when it is enabled you can change the
 * hash define SHIFTUNITS to 1
 *
 *
 * author: Stephen Privitera
 * contact: sprivite@ligo.caltech.edu
 * date: 11/20/2009
 *
 */
int main(int argc, char *argv[])
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
    UINT4 error_flag;
    UINT4 trial_max = 10*subBankSize;
    int exit_status = 0;

    /* verbosity level (warning: this program can be very verbose) */
    int verbose = ( (argc == 2) && (strcmp(argv[1],"dump") == 0)
		    ? 1
		    : 0 );

    /* Fill basic structures needed by XLALBankVetoCCMat */
    initBankVetoData(bankVetoData,ampVec,subBankSize,templateLength);

    /* Create mock templates which are orthogonal */
    makeDeltaFunctionTemplates(bankVetoData,subBankSize,templateLength,deltaF);

    /*
     * Make repeated calls to XLALBankVetoCCMat with time shifts
     * in some of the templates.
     */
    for (trial = 0; trial < trial_max; trial++)
    {
	/* set the time shift for all but one template to be 1 unit */
	timeshiftTemplates(bankVetoData->timeshift,deltaT,trial);

	/* Pass to CCMat function */
	XLALBankVetoCCMat (bankVetoData, ampVec, subBankSize,
			   dynRange, fLow, deltaF, deltaT);

	/* check whether output is as expected */
	error_flag = checkCCmatForErrors(bankVetoData->ccMat,bankVetoData->timeshift,deltaT);

	if ( error_flag )
	{
	    /* write out the offending matrix */
	    if ( verbose )
	    {
		sprintf(ccFileName,"cross_corr_mat_%d.txt",trial);
		writeCCmatToFile(bankVetoData->ccMat,subBankSize,ccFileName);
	    }

	    /* raise the exit status */
	    exit_status++;
	}
    }

    return exit_status;
}



/*
 * This function fills in the basic structures in the bankVetoData
 * structure (making most of them unity or zero).
 */
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
	bankVetoData->resp->data[sampleIndex].realf_FIXME = 1;
	bankVetoData->resp->data[sampleIndex].im = 0;
	bankVetoData->spec->data[sampleIndex] = 1;
	ampVec->data[sampleIndex] = 1;
    }

    /* initiate timeshift */
    bankVetoData->timeshift = XLALCreateREAL4Vector(subBankSize);

    return;

}



/*
 * This function fills the bankVetoData templates with the frequency domain
 * representation of time domain delta functions.  These function are orthogonal
 * and remain orthogonal under an integral number of deltaT time shifts.  The
 * cross-correlation matrix that results from this template bank and time shifts
 * of this bank will consist of 1's and 0's only.
 */
static void makeDeltaFunctionTemplates(FindChirpBankVetoData *bankVetoData,UINT4 subBankSize,UINT4 templateLength,REAL4 deltaF)
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
	    bankVetoData->fcInputArray[templateIndex]->fcTmplt->data->data[sampleIndex].realf_FIXME =
	      cos(phase);
	    bankVetoData->fcInputArray[templateIndex]->fcTmplt->data->data[sampleIndex].im =
	      sin(phase);
	}
	bankVetoData->fcInputArray[templateIndex]->fcTmplt->tmplt.fFinal =
	  (REAL4) (deltaF*templateLength-FLOW);
    }

    return;
}



/*
 * This function fills in the time shift array for the templates.
 * The timeshifts are forced to be integral multiples of deltaT so that
 * the templates (which are originally orthogonal) will remain orthogonal
 * under time shifting.
 */
static void timeshiftTemplates(REAL4Vector *timeshift,REAL4 deltaT,UINT4 trial)
{
    /* Force all time shifts to be integral multiples of deltaT
     * FIXME right now the bank veto does not do timeshifts, so these
     * are set to 0 because SHIFTUNITS = 0 as defined at the top,
     * change it when the bank veto supports time shifting.
     */

    UINT4 thisShift;
    UINT4 templateIndex;

    UINT4 shiftAmount = floor(trial/timeshift->length) + 1;
    UINT4 whoToShift = trial % timeshift->length;

    /* Choose a time shift for each template */
    for ( templateIndex = 0; templateIndex < timeshift->length; templateIndex++)
    {
	thisShift = ( (templateIndex == whoToShift )
		      ? shiftAmount
		      : 0 );

	timeshift->data[templateIndex] = SHIFTUNITS * ((REAL4)thisShift)*deltaT;
    }
    return;
}



/*
 * This function checks whether the cross-correlation matrix returned
 * by XLALBankVetoCCMat is consistent with the time shift array.
 */
static int checkCCmatForErrors(COMPLEX8Vector *ccMat,REAL4Vector *timeshift,REAL4 deltaT)
{
    REAL4 tolerance = TOLERANCE;
    REAL4 expectedValue;
    INT4 row;
    INT4 col;

    for ( row = 0; (UINT4)row < timeshift->length; row++)
    {
	for ( col = 0; (UINT4)col < timeshift->length; col++)
	{
	    /* The expected value depends on the specific way we have chosen to set up the template
	       bank, but does not depend on how we choose the time shifts, except that the time shifts
	       must be integral multiples of deltaT. */
	    expectedValue = (REAL4) ( ( (INT4)( (timeshift->data[row]-timeshift->data[col])/deltaT ) == (col - row) )
				      ? 1.0
				      : 0.0 );

	    /* real part will be 1 or 0 */
	    if ( fabs(expectedValue - crealf(ccMat->data[row*timeshift->length+col])) > tolerance)
		return 1;

	    /* imaginary part should always be zero */
	    if ( fabs( ccMat->data[row*timeshift->length+col].im) > tolerance)
	      return 1;

	}
    }
    return 0;
}



/*
 * This function can be used to write out cross-corr matrices.
 */
static int writeCCmatToFile(COMPLEX8Vector *ccmat,UINT4 subBankSize,CHAR *ccFileName)
{
    FILE *fdOut = fopen(ccFileName,"w");
    UINT4 row;
    UINT4 col;

    // print a matrix of the real parts
    fprintf(fdOut,"#Real part:\n");
    for ( row = 0; row< subBankSize; row++ )
    {
	for ( col = 0; col < subBankSize; col++ )
	{
	    fprintf(fdOut,"%.8f\t",crealf(ccmat->data[row*subBankSize+col]));
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

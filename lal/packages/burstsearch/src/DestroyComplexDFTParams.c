/******** <lalVerbatim file="DestroyComplexDFTParamsCV"> ********
Author: Flanagan, E
$Id$
********* </lalVerbatim> ********/


#include <lal/LALRCSID.h>

NRCSID (DESTROYCOMPLEXDFTPARAMSC, "$Id$");

#include <lal/ComplexFFT.h>
#include <lal/LALStdlib.h>
#include <lal/Sequence.h>
#include <lal/TFTransform.h>

/******** <lalVerbatim file="DestroyComplexDFTParamsCP"> ********/
void
XLALDestroyComplexDFTParams(
	ComplexDFTParams *dftParams
)
/******** </lalVerbatim> ********/
{
	if(dftParams) {
		XLALDestroyREAL4Sequence(dftParams->window);
		XLALDestroyCOMPLEX8FFTPlan(dftParams->plan);
		LALFree(dftParams);
	}
}

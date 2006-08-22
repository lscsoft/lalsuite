#include <stdio.h>

#include <lal/LALInspiral.h>
#include <lal/LALInspiralBank.h>
#include <lal/RealFFT.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>

#include <lal/LALConfig.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/Window.h>
#include <lal/TimeFreqFFT.h>
#include <lal/IIRFilter.h>
#include <lal/Date.h>
#include <lal/LIGOMetadataTables.h>




#define CVS_ID_STRING      "$Id$"
#define CVS_NAME_STRING    "$Name$"
#define CVS_REVISION       "$Revision$"
#define CVS_SOURCE         "$Source$"
#define CVS_DATE           "$Date$"
#define PROGRAM_NAME       "FFUtils"




/*  !!!!!  Structures to define   !!!!  */

 

RandomSignalIn




/*    !!!!      functions to define    !!!!  */ 

/* this function should compute the max length of template in the bank 
  and padd this size to the nearest power of 2 */
GetMaximumTemplateSize(&status, &bankIn , &checkLength); 

/* should compute the expected sensitivity (psd) using f-ns in the noise models */
GetNoisePSD(&status, bankIn.shf, ifo, df);

/* computes fourier trandform of the signal and template */
CreateFSignal(&status, &fwdp, signal, &fSignal);
CreateFTemplate(&status, &fwdp, bankCurrent, &fTemplate);

/** this function suppose to compute overlap */
ComputeOverlap(&status, &revp, &fSignal, &fTemplate, bankIn.shf, &match);


/**** !!!!!!!! The rest has to be modified. ********/


/* --- Some Error messages --- */
#define BANKEFFICIENCYNEW_ENORM  0
#define BANKEFFICIENCYNEW_ESUB   1
#define BANKEFFICIENCYNEW_EARG   2
#define BANKEFFICIENCYNEW_EVAL   3
#define BANKEFFICIENCYNEW_EFILE  4
#define BANKEFFICIENCYNEW_EINPUT 5
#define BANKEFFICIENCYNEW_EMEM   6

#define BANKEFFICIENCYNEW_MSGENORM  "Normal exit"
#define BANKEFFICIENCYNEW_MSGESUB   "Subroutine failed"
#define BANKEFFICIENCYNEW_MSGEARG   "Error parsing arguments"
#define BANKEFFICIENCYNEW_MSGEVAL   "Input argument out of valid range"
#define BANKEFFICIENCYNEW_MSGEFILE  "Could not open file"
#define BANKEFFICIENCYNEW_MSGEINPUT "Error reading file"
#define BANKEFFICIENCYNEW_MSGEMEM   "Out of memory"
#define BANKEFFICIENCYNEW_MSGPARSER "Missing arguments ??  "










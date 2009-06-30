/*
*  Copyright (C) 2007 Jolien Creighton
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/* <lalVerbatim file="LALfctCV">
 * Author: Edlund, Jeffrey A.
 * $Id$
 * </lalVerbatim>
 */
/* <lalLaTeX>
   \subsection{Module \texttt{LALfct.c}}
   \label{ss:LALfct.c}
     In this document I will need to refer to the different axes in the
     N-dimensional space.  I will use K[p] to the denote the different axes.
     K[0] is the axis usually refered to as the time or frequency axis.

     This package currently makes direct calls to FFTW to calculate the FFTs
     along the K[0] axis.

   \subsubsection*{Prototypes}
   \vspace{0.1in}
   </lalLaTeX> */

#include <config.h>
#if defined LAL_FFTW3_ENABLED
/* fftw3 not yet supported */
#else /* fftw2 implementation */

#ifdef HAVE_SFFTW_H
#include <sfftw.h>
#elif HAVE_FFTW_H
#include <fftw.h>
#else
#error "don't have either sfftw.h or fftw.h"
#endif

#include <lal/LALStdlib.h>  /* Include any required headers */
#include <lal/LALConstants.h>
#include <math.h>
#include <lal/LALfct.h>
#include <lal/FFTWMutex.h>
#include <lal/AVFactories.h>

/* Define RCS ID string */
NRCSID( LALFCTC, "$Id$" );

/* tell FFTW to use LALMalloc and LALFree */
#define FFTWHOOKS \
  do { fftw_malloc_hook = LALMallocShort; fftw_free_hook = LALFree; } while(0)


/* Define local constants and macros: */
/* This is the LALFCTPlan structure. */
struct LALFCTPlan {
  UINT2   numOfDims;           /* Number of terms and phase functions. */
  UINT4   numOfDataPoints;
  UINT4  *lengthOfDim;         /* An array of the lengths of the dimensions */
  UINT4  *maxDeltaOfDim;       /* An array of the largest Delta in each
				  dimension */
  LALFCT_FP  *phaseFuncForDim; /* An array of phase Functions */
  UINT4     **deltaNForDim;    /* An array of arrays telling when each dimension
				  should change */
  LALFCTCOMP  **prePhaseArray;    /* An array of arrays that contain
				     the Prephase data. */
  fftw_plan     planForFFT;    /* These fftw things will probably be changed */
  fftw_complex *work;          /* later. */
};


/* Declare local (static) functions (definitions can go here or at the
   end of the file) */
static INT2 fctMakeDeltaArray( LALFCTPlan *fctPlan, UINT2 dim );
static INT2 fctGenPrePhaseArray( LALFCTPlan *fctPlan, UINT2 dim );
/*
static INT2 fctGenRowIndex( LALFCTPlan *fctPlan, \
		     LALfctGenRowIndexParams *fctGenRowIndexParams, \
		     LALfctCalcOutput *fctCalcOutput );
 */
static INT2 fctGenRowData( LALFCTPlan *fctPlan, LALFCTCOMPVector *inputDataVector, \
		    LALfctCalcOutput *fctCalcOutput );
static INT4 CheckStatus( LALStatus *status);

/* Define global functions.  Function prototypes must be
   auto-extracted for inclusion in the documentation */
/* <lalVerbatim> */
void LALfctInitialize(LALStatus *status,
		      LALFCTPlan **fctPlan,
		      LALfctInitParams *fctInitParams
		      ){ /* </lalVerbatim> */
  /* <lalLaTeX>
     This function initializes the fctPlan structure. In needs to called before
     any of the other functions.
     \vspace{0.1in}


     </lalLaTeX>*/


  /* Initialize status structure to indicate nominal execution by default. */
  INITSTATUS( status, "LALfctInitialize", LALFCTC );

  FFTWHOOKS;

  /* Check to make sure that the pointer that we were passed will work. */
  ASSERT( fctPlan != 0, status, LALFCTH_EPOINTERIS0, LALFCTH_MSGEPOINTERIS0 );
  ASSERT( *fctPlan == 0, status, LALFCTH_EPOINTERNOT0, LALFCTH_MSGEPOINTERNOT0 );

  /* Allocate the memory for the fctPlan structure. */
  *fctPlan = (LALFCTPlan *) LALMalloc( sizeof( LALFCTPlan ) );

  /* Check to make sure that the memory was allocated. */
  if ( *fctPlan == NULL ){
    ABORT(status, LALFCTH_EALLOC_MEM, LALFCTH_MSGEALLOC_MEM);
  }

  /* Allocate memory for the arrays in the fctPlan. */
  /* Also Calculate the plan for the FFT. */
  (*fctPlan)->lengthOfDim = (UINT4 *) LALCalloc( fctInitParams->numOfDims, \
					   sizeof( UINT4 ) );
  (*fctPlan)->maxDeltaOfDim = (UINT4 *) LALCalloc( fctInitParams->numOfDims, \
					     sizeof( UINT4 ) );
  (*fctPlan)->phaseFuncForDim = \
    (LALFCT_FP *) LALCalloc( fctInitParams->numOfDims, \
			     sizeof( LALFCT_FP ) );
  LAL_FFTW_PTHREAD_MUTEX_LOCK;
  (*fctPlan)->planForFFT = fftw_create_plan( fctInitParams->numOfDataPoints, \
					     FFTW_FORWARD, FFTW_THREADSAFE | \
					     FFTW_ESTIMATE | FFTW_IN_PLACE);
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;

  (*fctPlan)->work = \
    (fftw_complex *) LALCalloc( fctInitParams->numOfDataPoints, \
				sizeof( fftw_complex ) );
  (*fctPlan)->deltaNForDim = (UINT4 **) LALCalloc( fctInitParams->numOfDims,\
						   sizeof(UINT4 *));
  (*fctPlan)->prePhaseArray = \
    (LALFCTCOMP **) LALCalloc(fctInitParams->numOfDims, \
			      sizeof(LALFCTCOMP *));

  /* Check to make sure that the memory was allocated and the that the
     FFT plan was calculated. */
  if ( (*fctPlan)->lengthOfDim && (*fctPlan)->maxDeltaOfDim && \
       (*fctPlan)->phaseFuncForDim && (*fctPlan)->planForFFT && \
       (*fctPlan)->work && (*fctPlan)->deltaNForDim && \
       (*fctPlan)->prePhaseArray == 0 ) {

    if ( (*fctPlan)->lengthOfDim != 0 ) {
     LALFree( (*fctPlan)->lengthOfDim );
    }
    if ( (*fctPlan)->maxDeltaOfDim != 0 ) {
      LALFree( (*fctPlan)->maxDeltaOfDim );
    }
    if ( (*fctPlan)->phaseFuncForDim != 0 ) {
      LALFree( (*fctPlan)->phaseFuncForDim );
    }
    if ( (*fctPlan)->planForFFT != 0 ) {
      LAL_FFTW_PTHREAD_MUTEX_LOCK;
      fftw_destroy_plan( (*fctPlan)->planForFFT );
      LAL_FFTW_PTHREAD_MUTEX_UNLOCK;
    }
    if ( (*fctPlan)->work != 0 ) {
      LALFree( (*fctPlan)->work );
    }
    if ( (*fctPlan)->deltaNForDim != 0 ) {
      LALFree( (*fctPlan)->deltaNForDim );
    }
    if ( (*fctPlan)->prePhaseArray != 0 ) {
      LALFree( (*fctPlan)->prePhaseArray );
    }
    LALFree( *fctPlan );
    ABORT(status, LALFCTH_EALLOC_MEM, LALFCTH_MSGEALLOC_MEM);
  }

  /* Initialize a few variables inside the fctPlan. */
  (*fctPlan)->numOfDims = fctInitParams->numOfDims;
  (*fctPlan)->numOfDataPoints = fctInitParams->numOfDataPoints;
  (*fctPlan)->lengthOfDim[0] = fctInitParams->numOfDataPoints;

  RETURN( status );
}

/* <lalVerbatim> */
void LALfctAddPhaseFunc(LALStatus *status,
			LALFCTPlan *fctPlan,
			LALfctAddPhaseFuncParams *fctAddPhaseFuncParams
			){ /* </lalVerbatim> */
  /* <lalLaTeX>
     This function puts the length of the dimension and the phaseFunc for that
     dimension into the fctPlan.

     This function also calls the internal functions fctMakeDeltaArray and fctGenPrePhaseArray.
     \vspace{0.1in}

     </lalLaTeX> */
  UINT2 dim;
  UINT4 length;
  INT2 check;

  /* Initialize status to indicate nominal execution by default. */
  INITSTATUS( status, "LALfctAddPhaseFunc", LALFCTC );

  /* Check our parameters. */
  ASSERT( fctAddPhaseFuncParams != 0, status, LALFCTH_EPOINTERIS0, \
	  LALFCTH_MSGEPOINTERIS0 );

  /* Get dim from fctAddPhaseFuncParams and check it to make sure that it's
     in the proper range. */
  dim = fctAddPhaseFuncParams->dim;
  if ((dim >= fctPlan->numOfDims)||(dim == 0)) {
    ABORT( status, LALFCTH_EDIM_OUT_OF_RANGE, LALFCTH_MSGEDIM_OUT_OF_RANGE );
  }

  /* Get length from fctAddPhaseFuncParams */
  length = fctAddPhaseFuncParams->lengthOfDim;

  /* Check the fctPlan */
  if (fctPlan == 0) {
    ABORT(status, LALFCTH_ECALL_INIT_FIRST,  LALFCTH_MSGECALL_INIT_FIRST);
  } else if (fctPlan->lengthOfDim && fctPlan->phaseFuncForDim == 0) {
    ABORT(status, LALFCTH_ECALL_INIT_FIRST,  LALFCTH_MSGECALL_INIT_FIRST);
  }

  /* Set the length of the dimension. */
  fctPlan->lengthOfDim[dim] = length;

  /* Set the phase function for the dimension. */
  fctPlan->phaseFuncForDim[dim] = fctAddPhaseFuncParams->phaseFuncForDim;

  /* Call fctMakeDeltaArray and check to make sure that it succeeded. */
  check = fctMakeDeltaArray(fctPlan, dim);
  if (check <= 0) {
    ABORT(status, LALFCTH_EMAKE_DELTA_ARRAY, \
	  LALFCTH_MSGEMAKE_DELTA_ARRAY);
  }

  /* Call fctGenPrePhaseArray and check to make sure that it succeeded. */
  check = fctGenPrePhaseArray(fctPlan, dim);
  if (check <= 0) {
    ABORT(status, LALFCTH_EGEN_PREPHASE_ARRAY, \
	  LALFCTH_MSGEGEN_PREPHASE_ARRAY);
  }

  RETURN( status );
}

/* <lalVerbatim> */
void LALfctCalc( LALStatus *status,
		 LALfctCalcOutput *fctCalcOutput,
		 LALFCTCOMPVector *inputDataVector,
		 LALfctCalcParams *fctCalcParams){ /* </lalVerbatim> */
  /* <lalLaTeX>
     This function Calculates the N-dimensional FCT.  It requires that
     LALfctInitialize and LALfctAddPhaseFunc were called to set up fctPlan
     first.

     It calls fctGenRowIndex to make the list of rows parallel to K[0] axis.
     This list of rows is actually stored in the columns of
     $fctCalcOutput \rightarrow rowIndex \rightarrow data$.  This is because LAL
     and FFTW use row-major order and I think that it's easier to deal with
     the different rows if they are in contiguous arrays.

     Then it calls fctGenRowData to create $fctCalcOutput \rightarrow outputData$.
     $fctCalcOutput \rightarrow OutputData$ is a 2 dimensional LALFCTCOMPArray that stores
     the the data corresponding to the list of indices in $fctCalcOutput /rightarrow rowIndex$
     in the columns of $fctCalcOutput \rightarrow outputData \rightarrow data$.

     FFTW is then called to calculate the FFT of all of the columns of
     $fctCalcOutput \rightarrow outputData \rightarrow data$.

     Remember: The columns of $fctCalcOutput \rightarrow outputData \rightarrow data$ and
     $fctCalcOutput \rightarrow rowIndex \rightarrow data$ correspond to the rows
     parallel to the K[0] axis.
     \vspace{0.1in}

     </lalLaTeX>
  */

  INT2 check;
  LALfctGenRowIndexOutput fctGenRowIndexOutput;

  /* Initialize status to indicate nominal execution by default. */
  INITSTATUS( status, "LALfctCalc", LALFCTC );

  FFTWHOOKS;

  /* Check the various parameters. */
  ASSERT(fctCalcOutput != 0, status, LALFCTH_EPOINTERIS0, \
	 LALFCTH_MSGEPOINTERIS0 );
  ASSERT(inputDataVector != 0, status, LALFCTH_EPOINTERIS0, \
	 LALFCTH_MSGEPOINTERIS0 );
  ASSERT(fctCalcParams != 0, status, LALFCTH_EPOINTERIS0, \
	 LALFCTH_MSGEPOINTERIS0 );
  ASSERT(fctCalcParams->fctPlan != 0, status, LALFCTH_EPOINTERIS0, \
	 LALFCTH_MSGEPOINTERIS0 );
  ASSERT(fctCalcParams->fctGenRowIndexParams != 0, status, \
	 LALFCTH_EPOINTERIS0, LALFCTH_MSGEPOINTERIS0 );
  ASSERT(fctCalcParams->fctGenRowIndexParams->fctPlan != 0, status, \
	 LALFCTH_EPOINTERIS0, LALFCTH_MSGEPOINTERIS0 );
  ASSERT(fctCalcOutput->rowIndex == 0, status, LALFCTH_EPOINTERNOT0, \
	 LALFCTH_MSGEPOINTERNOT0);
  ASSERT(fctCalcOutput->outputData == 0, status, LALFCTH_EPOINTERNOT0, \
	 LALFCTH_MSGEPOINTERNOT0);

  /*Fill in the fctGenRowIndexOutput structure */
  fctGenRowIndexOutput.fctCalcOutput = fctCalcOutput;

  /* Call fctGenRowIndex and check to make sure it succeeded. */
  if (fctCalcParams->fctGenRowIndexFunc) {
    fctCalcParams->fctGenRowIndexFunc(status, &fctGenRowIndexOutput, \
				      fctCalcParams->fctGenRowIndexParams);
  } else {
    LALfctGenRowIndex(status, &fctGenRowIndexOutput, \
		      fctCalcParams->fctGenRowIndexParams);
  }
  if (CheckStatus( status )) {
    RETURN( status );
  }


  /* Call fctGenRowData and check to make sure it succeeded. */
  check = fctGenRowData(fctCalcParams->fctPlan, inputDataVector, fctCalcOutput);
  if (check <= 0) {
    ABORT(status, LALFCTH_EGEN_ROW_DATA, LALFCTH_MSGEGEN_ROW_DATA);
  }

  /* Call fftw to calculate the fft of all the columns of
     fctCalcOutput->outputData->data.  Remember: these columns correspond
     to the rows parallel to the K[0] axis. */
  fftw( fctCalcParams->fctPlan->planForFFT, fctGenRowIndexOutput.numOfRows, \
	(fftw_complex *) fctCalcOutput->outputData->data, 1, \
	fctCalcParams->fctPlan->lengthOfDim[0], \
	fctCalcParams->fctPlan->work, 1, 0 );

  RETURN( status );
}

/* <lalVerbatim> */
void LALfctGenRowIndex(LALStatus *status, \
		       LALfctGenRowIndexOutput *fctGenRowIndexOutput, \
		       LALfctGenRowIndexParams *fctGenRowIndexParams) {
/* </lalVerbatim> */

  /* <lalLaTeX>
     This function generates a list of rows to calculate.

     This is the default fctGenRowIndex function. It just generates a list of
     all the rows in the entire space.

     If you want a smarter function you should write your own and pass a pointer
     to your new GenRowIndex function in the fctCalcParams structure.
     \vspace{0.1in}

     </lalLaTeX>
 */

  /* I want this function easy to copy and change for different wave types and
     dependencies.  To that end I'm making a few things more general. */

  UINT4 J, location;
  UINT4Vector *dimlength=0;
  /*LALStatus status={0};*/
  UINT4 *dimCursors=0;
  UINT4 totalRows=0;
  UINT4 skipCursor;
  UINT4 rowCursor=0;
  UINT2 numOfDims;
  UINT4 *lengthOfDim;
  BOOLEAN continueRowLooping;
  BOOLEAN continueIndexGenLooping;
  BOOLEAN numOfRowsKnown=0;

  /* Initialize status indicate nominal execution by default. */
  INITSTATUS( status, "LALfct: fctGenRowIndex", LALFCTC );

  LALfctAccessNumOfDims( status, fctGenRowIndexParams->fctPlan, &numOfDims );
  LALfctAccessLengthOfDims( status, fctGenRowIndexParams->fctPlan, &lengthOfDim );

  /* Create an array of cursors to use while generating the rowIndex. */
  dimCursors = (UINT4 *) LALCalloc(numOfDims - 1, sizeof(UINT4));
  if (dimCursors == 0) {
    ABORT( status, LALFCTH_EALLOC_MEM, LALFCTH_MSGEALLOC_MEM );
  }


  /* Determine the number of rows that will be calculated */
  if (!fctGenRowIndexParams->goToEndOfRows) {
    totalRows = fctGenRowIndexParams->numOfRows;
    numOfRowsKnown = 1;
  }

  continueIndexGenLooping = 1;
  do {
    if (numOfRowsKnown && fctGenRowIndexParams->createIndex) {
      /* Create a vector to use to create the array.
	 :TODO: This vector probably doesn't need to be created dynamically.*/
      LALU4CreateVector ( status, &dimlength, 2);
      if (CheckStatus(status)) {
	RETURN( status );
      }
      dimlength->data[0] = totalRows; /* LAL uses Row major order,
					 we should too.*/
      dimlength->data[1] = numOfDims - 1; /* I want the data in a
					     contiguous piece. */

      /* Create the Array to hold the index. */
      LALU4CreateArray ( status, &(fctGenRowIndexOutput->fctCalcOutput->rowIndex), \
		      dimlength );
      if (CheckStatus(status)) {
	LALU4DestroyVector ( status, &dimlength );
	RETURN( status );
      }

      /* Destroy the vector that was used to create the array. */
      LALU4DestroyVector ( status, &dimlength );
      if (CheckStatus(status)) {
	LALU4DestroyArray ( status, \
			 &(fctGenRowIndexOutput->fctCalcOutput->rowIndex));
	RETURN( status );
      }
    }


    continueRowLooping = 1;
    rowCursor = 0;
    skipCursor = 0;
    /* This loops until we've gone through all the possible rows. */
    do {

      /*      if (1) {  */ /* This will have a test to see if the row is one
		    we wish to calculate. */
      if (skipCursor >= fctGenRowIndexParams->skipRows) {
	if (fctGenRowIndexParams->createIndex&&numOfRowsKnown){
	  /* Add the row Indices to the list of rows to calculate. */
	  /* Find the location of the beginning of the index column. */
	  location =  rowCursor * (numOfDims - 1);

	  /* Add the indices to the list.  Note that J = 0 is dimension K[1]. */
	  for (J = 0; ((INT4)J) < numOfDims - 1; J++) {
	    fctGenRowIndexOutput->fctCalcOutput->rowIndex->data[location+J] = \
	      dimCursors[J];
	  }
	}
	rowCursor++;
      } else {
	skipCursor++;
      }

      /* } This is the end of the test block. */

      /* Increase the cursors for the K[1] dimension by 1. */
      dimCursors[0]++;

      /* loop through the dimensions. If any of them have reached their
	 maximum values, then set that dimension cursor to 0 and increase
	 the next dimension cursor.  If the last dimension has reached it's
	 maximum value than end. */
      for(J = 0; ((INT4)J) < numOfDims - 1; J++) {
	if(dimCursors[J] == lengthOfDim[J+1]) {
	  dimCursors[J] = 0;
	  if (((INT4)J) == numOfDims - 2) {
	    continueRowLooping = 0;
	  } else {
	    ++dimCursors[J+1];
	  }
	}
      }
      if (!fctGenRowIndexParams->goToEndOfRows) {
	if (rowCursor == totalRows) {
	  continueRowLooping = 0;
	}
      }

    } while (continueRowLooping);

    if (numOfRowsKnown) {
      if (rowCursor != totalRows) {
	printf("rowCursor = %d, totalRows = %d\n", rowCursor, totalRows);
	ABORT(status, LALFCTH_EGEN_ROW_INDEX_NUMROW_MISMATCH,
	      LALFCTH_MSGEGEN_ROW_INDEX_NUMROW_MISMATCH);
      }
      continueIndexGenLooping = 0;
    } else {
      numOfRowsKnown = 1;
      totalRows = rowCursor;
    }

    if (!fctGenRowIndexParams->createIndex) {
      continueIndexGenLooping = 0;
    }

  } while (continueIndexGenLooping);
  fctGenRowIndexOutput->numOfRows = totalRows;

  /* This is debug code that will print out the list of rows. */
  /*for (rowCursor = 0; rowCursor < TotalRows; rowCursor++) {
    printf("rC = %d :", rowCursor);
    location =  rowCursor * (fctPlan->numOfDims - 1);
    for (J = 0; J < fctPlan->numOfDims - 1; J++) {
    printf("%d ", fctCalcOutput->rowIndex->data[location+J]);
    }
    printf("\n");
    }*/


  LALFree(dimCursors);

  RETURN( status );
}


/* <lalVerbatim> */
void LALfctAccessNumOfDims( LALStatus *status, LALFCTPlan *fctPlan, \
			    UINT2 *numOfDims ) {
/* </lalVerbatim> */
  /* <lalLaTeX>
     Use this function to access the number of dimensions in an fctPlan.
     \vspace{0.1in}


     </lalLaTeX> */
/* Initialize status indicate nominal execution by default. */
  INITSTATUS( status, "LALfctAccessNumOfDims", LALFCTC );

  *numOfDims = fctPlan->numOfDims;

  RETURN( status );
}

/* <lalVerbatim> */
void LALfctAccessLengthOfDims( LALStatus *status, LALFCTPlan *fctPlan, \
			      UINT4 **lengthOfDim ) {
/* </lalVerbatim> */
 /* <lalLaTeX>
     Use this function to access the length of the dimensions in an fctPlan.
     \vspace{0.1in}

     </lalLaTeX> */

/* Initialize status indicate nominal execution by default. */
  INITSTATUS( status, "LALfctAccessNumOfDims", LALFCTC );

  *lengthOfDim = fctPlan->lengthOfDim;

  RETURN( status );
}

/* <lalVerbatim> */
void LALfctDestroyPlan( LALStatus *status,
			LALFCTPlan **fctPlan ){ /* </lalVerbatim> */
  /* <lalLaTeX>
     This function destroys the fctPlan.
     \vspace{0.1in}

     </lalLaTeX> */

  UINT2 J;

  /* Initialize status indicate nominal execution by default. */
  INITSTATUS( status, "LALfctDestroyPlan", LALFCTC );

  FFTWHOOKS;

  /* Systematically destroy the fctPlan. */
  if ( (*fctPlan) != 0 ) {
    if ( (*fctPlan)->work != 0 ) {
      LALFree( (*fctPlan)->work );
    }
    if ( (*fctPlan)->planForFFT != 0 ) {
      LAL_FFTW_PTHREAD_MUTEX_LOCK;
      fftw_destroy_plan( (*fctPlan)->planForFFT );
      LAL_FFTW_PTHREAD_MUTEX_UNLOCK;
    }
    if ((*fctPlan)->prePhaseArray != 0) {
      for (J = 1; J < (*fctPlan)->numOfDims; J++){
	if ((*fctPlan)->prePhaseArray[J] != 0) {
	  LALFree((*fctPlan)->prePhaseArray[J]);
	}
      }
      LALFree((*fctPlan)->prePhaseArray);
    }
    if ((*fctPlan)->deltaNForDim != 0) {
      for (J = 1; J < (*fctPlan)->numOfDims; J++){
	if ((*fctPlan)->deltaNForDim[J] != 0) {
	  LALFree((*fctPlan)->deltaNForDim[J]);
	}
      }
      LALFree((*fctPlan)->deltaNForDim);
    }
    if ( (*fctPlan)->phaseFuncForDim != 0 ) {
      LALFree( (*fctPlan)->phaseFuncForDim );
    }
    if ( (*fctPlan)->maxDeltaOfDim != 0 ) {
      LALFree( (*fctPlan)->maxDeltaOfDim );
    }
    if ( (*fctPlan)->lengthOfDim != 0 ) {
      LALFree( (*fctPlan)->lengthOfDim );
    }
    LALFree( (*fctPlan) );
    *fctPlan = 0;
  }

  RETURN( status );
}


static INT2 fctMakeDeltaArray( LALFCTPlan *fctPlan, UINT2 dim ) {
  /* This function generates the delta array for dimension dim.
     It is called by LALfctAddPhaseFunc.  */
  /* Much of this code was taken from Fredrick A. Jenet's Original code.*/
  UINT4 j;
  UINT4 jMinCurrentIndex;
  UINT4 jMinPreviousIndex;
  UINT4 sum;
  LALFCTREAL jMinUpperValue, jMinLowerValue, jMinMidValue;
  LALFCTREAL diff;
  LALFCTREAL value;
  LALFCTREAL phaseFuncAt0, phaseFuncAt1;

  /* If for whatever reason LALfctAddPhaseFunc was called for the same dimension
     twice, free the old array and make a new one. */
  if ( fctPlan->deltaNForDim[dim] != 0 ) {
    LALFree( fctPlan->deltaNForDim[dim] );
  }

  /* Allocate Memory for the delta array for this dimension and check to
     make sure it worked. */
  fctPlan->deltaNForDim[dim] = (UINT4 *) LALCalloc(fctPlan->lengthOfDim[dim], \
						sizeof(UINT4));
  if (fctPlan->deltaNForDim[dim] == 0) {
    printf("Problem Allocating Memory\n");
    return(0);
  }

  /* Store the values of the function at 0 and 1. */
  phaseFuncAt0 = fctPlan->phaseFuncForDim[dim](0);
  phaseFuncAt1 = fctPlan->phaseFuncForDim[dim](1);

  /* Initialize some of the variables. */
  jMinPreviousIndex = 0;
  fctPlan->maxDeltaOfDim[dim] = 0;
  sum = 0;

  for(j = 1; j < fctPlan->lengthOfDim[dim]; j++){

    value = (fftw_real)(j)/fctPlan->lengthOfDim[dim];
    jMinUpperValue = 1;
    jMinLowerValue = 0;
    diff = jMinUpperValue - jMinLowerValue;

    while(diff > .001) {

      jMinMidValue = jMinLowerValue + diff/2;

      if( ( fctPlan->phaseFuncForDim[dim](jMinMidValue) - phaseFuncAt0 ) / \
	 ( phaseFuncAt1 - phaseFuncAt0 ) > value ) {
	jMinUpperValue = jMinMidValue;
      } else {
	jMinLowerValue = jMinMidValue;
      }

      diff = jMinUpperValue - jMinLowerValue;

    }

    jMinCurrentIndex  = (int)floor(fctPlan->numOfDataPoints * jMinLowerValue);
    fctPlan->deltaNForDim[dim][j - 1] = jMinCurrentIndex - jMinPreviousIndex;
    jMinPreviousIndex = jMinCurrentIndex;

    if(fctPlan->deltaNForDim[dim][j - 1] > fctPlan->maxDeltaOfDim[dim]){
      fctPlan->maxDeltaOfDim[dim] = fctPlan->deltaNForDim[dim][j - 1];
    }
    sum += fctPlan->deltaNForDim[dim][j - 1];
  }

  fctPlan->deltaNForDim[dim][fctPlan->lengthOfDim[dim] - 1] = \
    fctPlan->numOfDataPoints - jMinPreviousIndex;

  if(fctPlan->deltaNForDim[dim][fctPlan->lengthOfDim[dim] - 1] > \
     fctPlan->maxDeltaOfDim[dim]){
    fctPlan->maxDeltaOfDim[dim] = \
      fctPlan->deltaNForDim[dim][fctPlan->lengthOfDim[dim] - 1];
  }
  sum += fctPlan->deltaNForDim[dim][fctPlan->lengthOfDim[dim] - 1];

  /*Check for consistency*/
  if(sum != fctPlan->lengthOfDim[0]){
    printf("ERROR:The sum of the delta_N[i] array does not\n");
    printf("      equal Length_K[0].\n");
    printf("      This is either an internal error or a problem\n");
    printf("      with the phase function specified.\n");
    return(0);
  }
  if(fctPlan->lengthOfDim[0] < fctPlan->maxDeltaOfDim[dim]){
    printf("ERROR:Length_K[0] is less then the maximum\n");
    printf("      interval length (= %d).\n",fctPlan->maxDeltaOfDim[dim]);
    printf("      Increase Length_K[0] to at least %d\n", \
	   fctPlan->maxDeltaOfDim[dim]);
    return(0);
  }
  return(1);
}

static INT2 fctGenPrePhaseArray( LALFCTPlan *fctPlan, UINT2 dim ) {
  /* This function generates the prePhaseArray for dimension dim.
     The prePhaseArray is used later in fctGenRowData to multiply the proper
     phase factors to each element of each row (parallel to K[0]). The result
     is the same as what you would get by taking the FFT of every dimension
     except K[0].  */


  UINT4 J;
  LALFCTREAL x = 0;
  LALFCTREAL realValue;
  LALFCTREAL imagValue;

  /* If LALfctAddPhaseFunc is called twice on the same dimension, free the old
     prePhaseArray and recalculate it. */
  if (fctPlan->prePhaseArray[dim] != 0) {
    LALFree(fctPlan->prePhaseArray[dim]);
  }

  /* Allocate the prePhase Array for this dimension and check to make sure
     it worked.*/
  fctPlan->prePhaseArray[dim] = \
    (LALFCTCOMP *) LALCalloc(fctPlan->lengthOfDim[dim]/2 + 1, \
			     sizeof(LALFCTCOMP));
  if (fctPlan->prePhaseArray[dim] == 0) {
    printf("Problem Allocating Memory\n");
    return(0);
  }

  /* It turns out that all of the phase factors that we need can be easily
     derived from fctPlan->lengthOfDim[dim]/2 + 1 complex numbers.  This is
     because exp(-2*PI*i*J/M) is periodic and also that exp(-2*PI*i*8/9) is
     the same as exp(2*PI*i*2/9). */
  for(J = 0; J <= fctPlan->lengthOfDim[dim]/2; J++) {

    x = 2*LAL_PI*(J)/fctPlan->lengthOfDim[dim];

    realValue = cos(x);
    imagValue = -sin(x); /* this is negative because I've simplified from
			    exp(-2*PI*i*J/fctPlan->lengthOfDim[dim]) */

    fctPlan->prePhaseArray[dim][J].re = realValue;
    fctPlan->prePhaseArray[dim][J].im = imagValue;

  }
  return(1);
}


static const LALStatus statusInit;

static INT2 fctGenRowData( LALFCTPlan *fctPlan, LALFCTCOMPVector *inputDataVector, \
		    LALfctCalcOutput *fctCalcOutput ) {
  /* This function generates the data for the rows (parallel to the K[0] axis)
     that are specified in fctCalcOutput->rowIndex */

  UINT4 J;
  UINT4 rowCursor=0;
  UINT4 *dataRowCursors;
  UINT4 *dataRowMaxs;
  UINT4 TotalRows;
  UINT4 indexLocation;
  UINT4 preLocation;
  UINT4 tempRowIndex;
  UINT4 rowLocation;
  LALFCTCOMP *outputDataPointer;
  LALFCTCOMP *inputDataPointer;
  UINT4 indexLocationRowLength;
  UINT4 rowLocationRowLength;

  UINT4Vector *dimlength=NULL;
  LALStatus status = statusInit;

  BOOLEAN CachedFlag;
  LALFCTCOMP *TempPhaseArray;
  LALFCTCOMP prePhaseData;
  LALFCTCOMP TempPhaseStorage;

  /* Get the Total number of rows. */
  TotalRows = fctCalcOutput->rowIndex->dimLength->data[0];

  /* Create the vector to create the array to store the data in. */
  LALU4CreateVector ( &status, &dimlength, 2);
  if (CheckStatus(&status)) {
    return(0);
  }
  dimlength->data[0] = TotalRows; /* LAL uses Row major order, we should too.*/
  dimlength->data[1] = fctPlan->lengthOfDim[0];

  /* Create the output array */
  LALFCTCOMPCreateArray ( &status, &(fctCalcOutput->outputData), dimlength );
  if (CheckStatus(&status)) {
    LALU4DestroyVector ( &status, &dimlength );
    return(0);
  }

  /* Destroy the vector that was used to create the array. */
  LALU4DestroyVector ( &status, &dimlength );
  if (CheckStatus(&status)) {
    LALFCTCOMPDestroyArray ( &status, &(fctCalcOutput->outputData));
    return(0);
  }

  /* Allocate the Temporary Phase Array.  (We use this as a cache
     so that we don't have to multiply all of the phase factors together
     each time.) */
  TempPhaseArray = (LALFCTCOMP *) LALCalloc(TotalRows, sizeof(LALFCTCOMP));

  /* Allocate the dataRowCursors.  We use these to keep track of where the the
     data would be in the N-dimensional space. */
  dataRowCursors = (UINT4 *) LALCalloc(fctPlan->numOfDims, sizeof(UINT4));

  /* Allocate the dataRowMaxs.  These keep track of when the dataRowCursors need
     to change. */
  dataRowMaxs = (UINT4 *) LALCalloc(fctPlan->numOfDims, sizeof(UINT4));

  /* Check to make sure that they were allocated. */
  if (TempPhaseArray && dataRowCursors && dataRowMaxs == 0) {
    if (TempPhaseArray != 0) {
      LALFree(TempPhaseArray);
    }
    if (dataRowCursors != 0) {
      LALFree(dataRowCursors);
    }
    if (dataRowMaxs != 0) {
      LALFree(dataRowMaxs);
    }
    return(0);
  }

  /* Initialize the dataRow Maxs */
  dataRowMaxs[0] = fctPlan->lengthOfDim[0];
  for(J = 1; J < fctPlan->numOfDims; J++) {
    dataRowMaxs[J] =  fctPlan->deltaNForDim[J][dataRowCursors[J]];
  }

  /* These things make it a bit faster. */
  outputDataPointer = fctCalcOutput->outputData->data;
  inputDataPointer = inputDataVector->data;
  indexLocationRowLength = fctCalcOutput->rowIndex->dimLength->data[1];
  rowLocationRowLength = fctCalcOutput->outputData->dimLength->data[1];

  /* this part takes longer than any other part of the program.  (Other than the
     fft call.) It can still be optimized. */

  /* Make sure that the Cache will be calculated the first time through.*/
  CachedFlag = 0;

  /* Move along the K[0] axis */
  for (dataRowCursors[0] = 0; dataRowCursors[0] < fctPlan->lengthOfDim[0]; \
	 dataRowCursors[0]++) {

    /* As we move along the K[0] axis, check to see if we need to increase any
       of the dataRowCursors. */
    for(J = 1; J < fctPlan->numOfDims; J++) {
      if(dataRowCursors[0] == dataRowMaxs[J]) {
	dataRowCursors[J]++;
	dataRowMaxs[J] += fctPlan->deltaNForDim[J][dataRowCursors[J]];
	/* Because of this change we need to recalculate the Cache. */
	CachedFlag = 0;
      }
    }

    /*Reset these locations. */
    indexLocation = 0;
    rowLocation = 0;


    if (!CachedFlag) { /* We need to recalculate the TempPhaseArray. */

      /* Loop over the rows.*/
      for (rowCursor = 0; rowCursor < TotalRows; rowCursor++) {

	TempPhaseArray[rowCursor].re = 1;
	TempPhaseArray[rowCursor].im = 0;

	for(J = 1; J < fctPlan->numOfDims; J++){
	  tempRowIndex = fctCalcOutput->rowIndex->data[indexLocation+J-1];
	  preLocation = (dataRowCursors[J] * tempRowIndex) % \
	    fctPlan->lengthOfDim[J];
	  if (preLocation <= fctPlan->lengthOfDim[J]/2) {
	    prePhaseData = fctPlan->prePhaseArray[J][preLocation];
	  } else {
	    prePhaseData = \
	      fctPlan->prePhaseArray[J][fctPlan->lengthOfDim[J]/2 - \
				       preLocation % \
				       (fctPlan->lengthOfDim[J]/2)];
	    prePhaseData.im = -prePhaseData.im;
	  }
	  TempPhaseStorage = TempPhaseArray[rowCursor];
	  TempPhaseArray[rowCursor].re = \
	    TempPhaseStorage.re * prePhaseData.re - \
	    TempPhaseStorage.im * prePhaseData.im;
	  TempPhaseArray[rowCursor].im = \
	    TempPhaseStorage.re * prePhaseData.im + \
	    TempPhaseStorage.im * prePhaseData.re;
	}

	outputDataPointer[rowLocation+dataRowCursors[0]].re = \
	  TempPhaseArray[rowCursor].re * \
	  inputDataPointer[dataRowCursors[0]].re - \
	  TempPhaseArray[rowCursor].im *\
	  inputDataPointer[dataRowCursors[0]].im;

	outputDataPointer[rowLocation+dataRowCursors[0]].im = \
	  TempPhaseArray[rowCursor].re * \
	  inputDataPointer[dataRowCursors[0]].im + \
	  TempPhaseArray[rowCursor].im *\
	  inputDataPointer[dataRowCursors[0]].re;

	indexLocation += indexLocationRowLength;
	rowLocation += rowLocationRowLength;
      }
      CachedFlag = 1;
    } else {
      for (rowCursor = 0; rowCursor < TotalRows; rowCursor++) {
	outputDataPointer[rowLocation+dataRowCursors[0]].re = \
	  TempPhaseArray[rowCursor].re * \
	  inputDataPointer[dataRowCursors[0]].re - \
	  TempPhaseArray[rowCursor].im *\
	  inputDataPointer[dataRowCursors[0]].im;

	outputDataPointer[rowLocation+dataRowCursors[0]].im = \
	  TempPhaseArray[rowCursor].re * \
	  inputDataPointer[dataRowCursors[0]].im + \
	  TempPhaseArray[rowCursor].im *\
	  inputDataPointer[dataRowCursors[0]].re;

	indexLocation += indexLocationRowLength;
	rowLocation += rowLocationRowLength;
      }
    }

  }
  LALFree(TempPhaseArray);
  LALFree(dataRowCursors);
  LALFree(dataRowMaxs);
  return(1);
}

static INT4 CheckStatus( LALStatus *status) {
  /* This function reports nonzero values of status. */
  if ( status->statusCode != 0) {
    REPORTSTATUS( status );
    return(1);
  }
  return(0);
}

#endif /* fftw2 implementation */

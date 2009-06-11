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

/*  <lalVerbatim file="LALfctHV">
    Author: Edlund, Jeffrey A.
    $Id$
    </lalVerbatim>
*/

/* <lalLaTeX>

\section{Header \texttt{LALfct.h}}
\label{s:LALfct.h}

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALfct.h>
\end{verbatim}

</lalLaTeX>

 */

#ifndef _LALFCT_H    /* Protect against double-inclusion */
#define _LALFCT_H

#include <lal/LALStdlib.h>  /* Include any other headers */

#ifdef __cplusplus
extern "C" {            /* Protect against C++ name mangling */
#endif

/* Define the RCS ID string */
  NRCSID(LALFCTH,"$Id$");

/*
<lalLaTeX>
%\subsection*{Error conditions}
%\input{LALfctHErrTab}
</lalLaTeX>
*/

/* Define error codes and messages.  These must be auto-extracted for
   inclusion into the documentation */

/* <lalErrTable file="LALfctHErrTab"> */

#define LALFCTH_EPOINTERNOT0 1
#define LALFCTH_EPOINTERIS0 2
#define LALFCTH_EALLOC_MEM 3
#define LALFCTH_ECALL_INIT_FIRST 4
#define LALFCTH_EINTERNAL_SUB 5
#define LALFCTH_EDIM_OUT_OF_RANGE 6
#define LALFCTH_EMAKE_DELTA_ARRAY 7
#define LALFCTH_EGEN_PREPHASE_ARRAY 8
#define LALFCTH_EGEN_ROW_INDEX 9
#define LALFCTH_EGEN_ROW_DATA 10
#define LALFCTH_EGEN_ROW_INDEX_NUMROW_MISMATCH 11

#define LALFCTH_MSGEPOINTERNOT0 "A pointer that was supposed to be 0 was not."
#define LALFCTH_MSGEPOINTERIS0  "Pointer = 0"
#define LALFCTH_MSGEALLOC_MEM "Problem Allocating Memory."
#define LALFCTH_MSGECALL_INIT_FIRST "Make sure that you call LALfctInitialize first."
#define LALFCTH_MSGEINTERNAL_SUB "Problem in internal subroutine."
#define LALFCTH_MSGEDIM_OUT_OF_RANGE "Dim parameter is out of range."
#define LALFCTH_MSGEMAKE_DELTA_ARRAY "Problem in fctMakeDeltaArray subroutine."
#define LALFCTH_MSGEGEN_PREPHASE_ARRAY "Problem in fctGenPrePhaseArray subroutine."
#define LALFCTH_MSGEGEN_ROW_INDEX "Problem in fctGenRowIndex subroutine."
#define LALFCTH_MSGEGEN_ROW_DATA "Problem in fctGenRowData subroutine."
#define LALFCTH_MSGEGEN_ROW_INDEX_NUMROW_MISMATCH "The number of rows requested does not match the number of rows found."

/* </lalErrTable> */

/* Define other global constants or macros */

/* Define new structures and types */

/*<lalLaTeX>
  \subsection*{Structures and Datatypes}
  The following LALFCT data types exist so that we can easily change from
  float to double precision.
  </lalLaTeX>
  <lalVerbatim> */
#define LALFCTREAL REAL4
#define LALFCTCOMP COMPLEX8
#define LALFCTCOMPArray COMPLEX8Array
#define LALFCTCOMPVector COMPLEX8Vector
/* </lalVerbatim>
<lalLaTeX>
The following function names are replaced with LAL functions that
correspond to the data types used above.
</lalLaTeX>
<lalVerbatim> */
#define LALFCTCOMPCreateVector LALCCreateVector
#define LALFCTCOMPDestroyVector LALCDestroyVector
#define LALFCTCOMPCreateArray LALCCreateArray
#define LALFCTCOMPDestroyArray LALCDestroyArray
/* </lalVerbatim>
<lalLaTeX>
   LALFCT\_FP is a pointer to a function that returns a LALFCTREAL and takes a
   LALFCTREAL as an argument.  Your phase functions must match this so that
   they can be put in an array inside the LALFCTPlan.
</lalLaTeX>

<lalVerbatim>
*/
typedef LALFCTREAL (*LALFCT_FP)(LALFCTREAL);

/*
</lalVerbatim>

<lalLaTeX>

LALFCTPlan is used internally to store information that will
be used in the calculation of the FCT.

</lalLaTeX> */

/*The structure of LALFCTPlan is defined in LALfct.c but other programs
  should know that it exists. */
typedef struct LALFCTPlan LALFCTPlan;

/* <lalLaTeX>
 LALfctInit takes the following structure as one of  it's arguments.
 </lalLaTeX>
 <lalVerbatim>*/
typedef struct tagLALfctInitParams {
  UINT4         numOfDataPoints;
  UINT2         numOfDims;
} LALfctInitParams;
/* </lalVerbatim>
   <lalLaTeX>
Where numofDataPoints is the size of the input data array and numOfDims
is the number of phase functions that will be used in the FCT.

LALfctAddPhaseFunc takes the following structure as one of it's
arguments.
</lalLaTeX>
<lalVerbatim>*/
typedef struct tagLALfctAddPhaseFuncParams {
  UINT2         dim;
  UINT4         lengthOfDim;
  LALFCT_FP     phaseFuncForDim;
} LALfctAddPhaseFuncParams;
/* </lalVerbatim>
   <lalLaTeX>
Where dim is the dimension of the N dimensional
output array that you would like to set the phase function for.
</lalLaTeX>

<lalLaTeX>
LALfctGenRowIndex takes the following structure as one of it's arguments.
</lalLaTeX>
 <lalVerbatim> */
typedef struct tagLALfctGenRowIndexParams {
  LALFCTPlan *fctPlan;
  UINT4 skipRows;
  UINT4 numOfRows; /* This only works if goToEndOfRows is FALSE */
  BOOLEAN goToEndOfRows; /* This means that the function should ignore
			    numOfRows. */
  BOOLEAN createIndex;
} LALfctGenRowIndexParams;
  /* </lalVerbatim>
 <lalLaTeX>
     Where skipRows is the number of rows in the N-dim space that you want
     to skip over and numOfRows is the number of rows you want the function
     to output indices for. If goToEndOfRows is true, then the function
     ignores numOfRows.  If createIndex is FALSE then the function will not
     create the index, it will only the total number of rows after skipRows.

     The following is the LALfctCalcOutput structure.
     </lalLaTeX>
<lalVerbatim> */
typedef struct tagLALfctCalcOutput {
  UINT4Array *rowIndex;
  LALFCTCOMPArray *outputData;
} LALfctCalcOutput;
  /* </lalVerbatim> */

  /* <lalLaTeX>
     The following structure is returned from the LALfctGenRowIndexOutput.
     </lalLaTeX>
<lalVerbatim>*/
typedef struct tagLALfctGenRowIndexOutput {
  LALfctCalcOutput *fctCalcOutput;
  UINT4 numOfRows;
}LALfctGenRowIndexOutput;
  /* </lalVerbatim>*/

  /* <lalLaTeX>
     It is possible to replace the standard LALFCTGenRowIndex with your
     own function.  This is useful if you want more control over the selection
     of the rows that you wish calculated. The following is the typedef
     for that function.
</lalLaTeX>
 <lalVerbatim> */

typedef void (*LALFCTGenRowIndex_FP)(LALStatus *status, \
		       LALfctGenRowIndexOutput *fctGenRowIndexOutput, \
		       LALfctGenRowIndexParams *fctGenRowIndexParams);
  /* </lalVerbatim> */

/* <lalLaTeX>
The following is the parameter structure for LALfctCalc.
</lalLaTeX>
 <lalVerbatim>*/
typedef struct tagLALfctCalcParams {
  LALFCTPlan    *fctPlan;
  LALfctGenRowIndexParams *fctGenRowIndexParams;
  LALFCTGenRowIndex_FP fctGenRowIndexFunc;
} LALfctCalcParams;
  /* </lalVerbatim> */




/* Include external global variables */

       /* Declare global function prototypes */
void LALfctInitialize( LALStatus *status,
		       LALFCTPlan **fctPlan,
		       LALfctInitParams *fctInitParams );
void LALfctAddPhaseFunc( LALStatus *status,
			 LALFCTPlan *fctPlan,
			 LALfctAddPhaseFuncParams *fctAddPhaseFuncParams );
void LALfctCalc( LALStatus *status,
		 LALfctCalcOutput *fctCalcOutput,
		 LALFCTCOMPVector *data,
		 LALfctCalcParams *fctCalcParams );
void LALfctDestroyPlan( LALStatus *status,
			LALFCTPlan **fctPlan );
void LALfctGenRowIndex(LALStatus *status, \
		       LALfctGenRowIndexOutput *fctGenRowIndexOutput, \
		       LALfctGenRowIndexParams *fctGenRowIndexParams);
void LALfctAccessNumOfDims( LALStatus *status, LALFCTPlan *fctplan, \
			    UINT2 *numOfDims );
void LALfctAccessLengthOfDims( LALStatus *status, LALFCTPlan *fctplan, \
			       UINT4 **lengthOfDim );

/* <lalLaTeX>
\newpage\input{LALfctC}
</lalLaTeX> */


#ifdef __cplusplus
}                   /* Close C++ protection */
#endif

#endif              /* Close double-include protection */









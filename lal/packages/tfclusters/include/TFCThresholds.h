/*-----------------------------------------------------------------------
 *
 * File Name: TFCThresholds.h
 *
 * Author: Julien Sylvestre
 *
 * Revision: $Id$ 
 *
 -----------------------------------------------------------------------*/


/************************************ <lalVerbatim file="TFCTHRESHOLDSH">
Author: Sylvestre, J.
$Id$
************************************* </lalVerbatim> */

#ifndef _TFCTHRESHOLDS_H
#define _TFCTHRESHOLDS_H

#include <lal/LALDatatypes.h>
#include <lal/LALStatusMacros.h>
#include "lal/LALRCSID.h"

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

NRCSID (TFCTHRESHOLDSH, "$Id$");


 /******************************** <lalErrTable file="TFCThresholdsHErrTab"> */
#define TFCTHRESHOLDSH_ENULLP       1
#define TFCTHRESHOLDSH_E01 128
#define TFCTHRESHOLDSH_EEGOAL 129
#define TFCTHRESHOLDSH_ENERR 130

#define TFCTHRESHOLDSH_MSGENULLP "Null pointer"
#define TFCTHRESHOLDSH_MSGE01 "Argument must be in [0,1]"
#define TFCTHRESHOLDSH_MSGEEGOAL "Error goal smaller than numerical precision"
#define TFCTHRESHOLDSH_MSGENERR "A numerical error occured"
/*************************************************** </lalErrTable> */



/*************************************<lalLaTeX file="TFCThresholdsStructs">
\subsubsection*{struct \texttt{RiceThreshold}}
\noindent A container for the parameters used in the computation of the thresholds.
\begin{description}
\item[\texttt{nFreq}] Number of frequency bins.
\item[\texttt{meanRe}] Mean values of the real part of the Fourier transform of the data, as a function of frequency (size = \texttt{nFreq}).
\item[\texttt{varRe}] Variances of the real part of the Fourier transform of the data, as a function of frequency (size = \texttt{nFreq}).
\item[\texttt{meanIm, varIm}] As above for imaginary part.
\item[\texttt{bpp}] The target black pixel probability.
\item[\texttt{eGoal}] The absolute error goal on the power thresholds.
\end{description}
******************************************************* </lalLaTeX> */
typedef struct tagRiceThresholdParams {

  UINT4 nFreq;

  REAL4* meanRe;

  REAL4* meanIm;

  REAL4* varRe;
  
  REAL4* varIm;

  REAL4 bpp;

  REAL4 eGoal;

} RiceThresholdParams;


void
LALTFCRiceThreshold ( LALStatus *status,
		      REAL4* rho,
		      RiceThresholdParams* thr
		      );

#ifdef  __cplusplus
}
#endif

#endif

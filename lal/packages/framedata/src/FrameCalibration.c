/**** <lalVerbatim file="FrameCalibrationCV">
 * Author: Brown, D. A.
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \subsection{Module \texttt{FrameCalibration.c}}
 * \label{ss:FrameCalibration.c}
 *
 * This module contains code used to extract calibration information contained
 * in frame files, and to construct a response (or transfer) function.  This
 * is supposed to provide a high-level interface for search authors to obtain a
 * response function in the desired form.
 *
 * \subsection*{Prototypes}
 * \input{FrameCalibrationCP}
 *
 * The routine \texttt{LALFrameExtractResponse()} extracts the necessary
 * calibration information from the frames. The frames used to construct
 * the calibration are located using the specified catalog. The function
 * constructs a response function (as a frequency series) from this
 * information.  If the fourth argument is non-\texttt{NULL} then this
 * string specifies the detector (H1, H2, L1, etc.) for which the calibration
 * is required.
 *
 * Certain fields of the output should be set before this routine is
 * called.  In particular: 1. The epoch field of the frequency series
 * should be set to the correct epoch so that the routine can generate
 * a response function tailored to that time (accounting for calibration
 * drifts).  2. The units of the response function should be set to
 * be either strain-per-count (for a response function) or count-per-strain
 * (for a transfer function); the routine will then return either the response
 * function or its inverse depending on the specified units.  Furthermore,
 * the power-of-ten field of the units is examined to scale the response
 * function accordingly.  3. The data vector should be allocated to the
 * required length and the frequency step size should be set to the required
 * value so that the routine can interpolate the response function to the
 * required frequencies.
 *
 * \vfill{\footnotesize\input{FrameCalibrationCV}}
 *
 **** </lalLaTeX> */
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/Calibration.h>
#include <lal/FrameCalibration.h>

NRCSID( FRAMECALIBRATIONC, "$Id$" );

/* <lalVerbatim file="FrameCalibrationCP"> */
void
LALExtractFrameResponse(
    LALStatus               *status,
    COMPLEX8FrequencySeries *output,
    const CHAR              *catalog,
    const CHAR              *ifo
    )
{ /* </lalVerbatim> */
  UINT4 k;
  const LALUnit strainPerCount = {0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};

  INITSTATUS( status, "LALFrameExtractResponse", FRAMECALIBRATIONC );
  ATTATCHSTATUSPTR( status );

  ASSERT( output, status, 
      FRAMECALIBRATIONH_ENULL, FRAMECALIBRATIONH_MSGENULL );
  ASSERT( output->data, status, 
      FRAMECALIBRATIONH_ENULL, FRAMECALIBRATIONH_MSGENULL );
  ASSERT( output->data->data, status, 
      FRAMECALIBRATIONH_ENULL, FRAMECALIBRATIONH_MSGENULL );
  ASSERT( catalog, status, 
      FRAMECALIBRATIONH_ENULL, FRAMECALIBRATIONH_MSGENULL );
  ASSERT( ifo, status, 
      FRAMECALIBRATIONH_ENULL, FRAMECALIBRATIONH_MSGENULL );

  output->sampleUnits = strainPerCount;

  for ( k = 0; k < output->data->length; ++k)
  {
    output->data->data[k].re = 1.0;
    output->data->data[k].im = 0;
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}

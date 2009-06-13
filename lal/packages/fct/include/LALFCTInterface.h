/*
*  Copyright (C) 2007 Philip Charlton
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

#ifndef _LALFCTINTERFACE_H
#define _LALFCTINTERFACE_H

#include <lal/LALDatatypes.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif

NRCSID(LALFCTINTERFACEH, "$Id$");

/*
<lalVerbatim file="LALFCTInterfaceHV">
Author: Charlton, P. R.
$Id$
</lalVerbatim>
*/

/*
<lalLaTeX>

\section{Header \texttt{LALFCTInterface.h}}
\label{s:LALFCTInterface.h}

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALFCTInterface.h>
\end{verbatim}

This header provides a LAL front-end to the ``generalised''
Fast Chirp Transform (FCT) of a
data set using an arbitrary number of phase functions. The package
{\tt fct} currently provided in LAL is restricted to phase functions
which are monotonic. The generalised FCT is valid for any phase function.

</lalLaTeX>
*/

/*
<lalLaTeX>
\subsection*{Error conditions}
\input{LALFCTInterfaceHErrTab}
</lalLaTeX>
*/

/* <lalErrTable file="LALFCTInterfaceHErrTab"> */

#define LALFCTINTERFACEH_ENULL 1
#define LALFCTINTERFACEH_ENNUL 2
#define LALFCTINTERFACEH_ESIZE 4
#define LALFCTINTERFACEH_ESZMM 8
#define LALFCTINTERFACEH_ENEGA 32
#define LALFCTINTERFACEH_ENPOS 64
#define LALFCTINTERFACEH_ECUBE 128
#define LALFCTINTERFACEH_EDATATYPE 256
#define LALFCTINTERFACEH_EINTERNAL 512

#define LALFCTINTERFACEH_MSGENULL "Null pointer"
#define LALFCTINTERFACEH_MSGENNUL "Non-null pointer"
#define LALFCTINTERFACEH_MSGESIZE "Invalid input size"
#define LALFCTINTERFACEH_MSGESZMM "Size mismatch"
#define LALFCTINTERFACEH_MSGENEGA "Negative number"
#define LALFCTINTERFACEH_MSGENPOS "Non-positive number"
#define LALFCTINTERFACEH_MSGECUBE "Invalid data cube bounds"
#define LALFCTINTERFACEH_MSGEDATATYPE "Internal FCT data type incompatible with REAL4"
#define LALFCTINTERFACEH_MSGEINTERNAL "Internal FCT error"

/* </lalErrTable> */

/*
<lalLaTeX>

\subsection*{Structures and Datatypes}
The following LALFCT data types exist so that we can easily change from
float to double precision.

</lalLaTeX>
*/

/* <lalVerbatim> */

/* Forward declarations */
struct tagLALFCTPlan;
typedef struct tagLALFCTPlan LALFCTPlan;

/*
  Typedef for phase functions.
  All phase functions must have the signature

  REAL4 phi(REAL4)
*/
typedef REAL4 (*LALFCTPhaseFn)(REAL4);

/* Maximum number of dimensions allowed in this implementation */
#define LALFCT_MAX_DIMS 8

/*
  Structure to contain parameters for initialising the FCT
*/
typedef struct
tagLALCreateFCTPlanInput
{
    /* Number of samples in input data. Must be positive */
    INT4 data_length;
    /* Number of dimensions. Must be > 1 */
    INT4 number_of_dimensions;
    /* Step through the input data with this stride. Must be positive */
    INT4 dimension_0_stride;
    /* An array of pointers to phase functions */
    LALFCTPhaseFn phase_fn[LALFCT_MAX_DIMS];
}
LALCreateFCTPlanInput;

/*
  Input structure for LALFCTSetOversamplingFactor

*/
typedef struct
tagLALFCTSetOversamplingFactorInput
{
    INT4 ofac;
    INT4 dim;
}
LALFCTSetOversamplingFactorInput;

/*
  Structure to hold information about the limits of a "data cube",
  that is, a parallelpiped (not necessarily cubic!) in N dimensions.
  The bounds for dimension i are given by start_locations[i] and
  end_locations[i].
*/
typedef struct
tagLALFCTDataCube
{
    /* List of start locations indexed by dimension. Must be non-negative */
    INT4 start_locations[LALFCT_MAX_DIMS];
    /*
      List of end locations indexed by dimension. Must be larger than
      the corresponding start_location
    */
    INT4 end_locations[LALFCT_MAX_DIMS];
    /* Stride for each dimension. Must be positive */
    INT4 stride[LALFCT_MAX_DIMS];
}
LALFCTDataCube;

/*
  Input structure for LALFCTSetDataCubes
*/
typedef struct
tagLALFCTSetDataCubesInput
{
    /* The address of an array of data cubes */
    LALFCTDataCube* data_cube;
    /* The number of cubes in the array */
    UINT4 num_data_cubes;
}
LALFCTSetDataCubesInput;

/*
  Input structure for LALFCTSetUnits

  The actual values at which the phase functions will be evaluated are

  phi[i] = phi(offset + delta*i);
*/
typedef struct
tagLALFCTSetUnitsInput
{
    REAL4 offset;
    REAL4 delta;
}
LALFCTSetUnitsInput;

/* </lalVerbatim> */

/*
<lalLaTeX>
\vfill{\footnotesize\input{LALFCTInterfaceHV}}
%\newpage\input{LALSampleC}
</lalLaTeX>
*/

/*
<lalLaTeX>
\subsection*{Prototypes}
</lalLaTeX>
*/

/* <lalVerbatim> */

/*
  Create the LALFCTPlan

  On input, *plan_ptr must have the value 0. On exit, it will be a
  pointer to a valid FCT plan.

  This function must be called before the plan can be used in any other
  function.
*/
extern
void LALCreateFCTPlan(LALStatus* const status,
		      LALFCTPlan** plan_ptr,
		      const LALCreateFCTPlanInput* const in);

/*
  Destroy the given plan.

  All memory associated with the plan is freed and the plan pointer
  is set to 0
*/
extern
void LALDestroyFCTPlan(LALStatus* const status,
		       LALFCTPlan** plan);

/*
  Set the oversampling factor for each phase function

  If this function is not called, the default oversampling factor is 1.
*/
extern
void LALFCTSetOversamplingFactor(LALStatus* const status,
			      const LALFCTSetOversamplingFactorInput* const in,
			      LALFCTPlan* const plan);

/*
  Set the units for the plan.

  If this function is not called, the default values are offset = 0
  and delta = 1.
*/
extern
void LALFCTSetUnits(LALStatus* const status,
		    const LALFCTSetUnitsInput* const in,
		    LALFCTPlan* const plan);

/*
  Set the maximum number of FFT's which will be done at one time
  by the FCT. Setting this to a value greater than 1 can improve
  the efficiency of the FCT.

  If this function is not called, the default value is 1.
*/
extern
void LALFCTSetMaxSegments(LALStatus* const status,
			  const INT4 max_segments,
			  LALFCTPlan* const plan);

/*
  Set the regions (data cube) for which the FCT is calculated.

  The input structure provides a list of data cubes. The FCT will
  only be calculated for the specific regions specified in the data
  cubes.

  At least one data cube must be set before the FCT can be calculated.
*/
extern
void LALFCTSetDataCubes(LALStatus* const status,
		       const LALFCTSetDataCubesInput* const in,
		       LALFCTPlan* const plan);

/*
  Obtain the size of the vector needed to hold the output.

  On return, output_data_size is the total length of the complex vector
  needed to hold the result of the FCT.
  Since the output size of the FCT depends on the data cube(s) being
  calculated, this must be obtained after calling LALFCTSetDataCubes.
*/
extern
void LALFCTOutputDataSize(LALStatus* const status,
			  UINT8* const output_data_size,
			  LALFCTPlan* const plan);

/*
  Calculate the FCT

  The input must be a vector with length data_length,
  as set in the LALCreateFCTPlanInput structure. The output must be
  a vector with the length given by LALFCTOutputDataSize.
*/
extern
void LALFCTCalculate(LALStatus* const status,
		     COMPLEX8Vector* const out,
		     const COMPLEX8Vector* const in,
		     LALFCTPlan* const plan);

/* </lalVerbatim> */

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _LALFCTINTERFACE_H */

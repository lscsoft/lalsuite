/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton
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

/**** <lalVerbatim file="ComplexFFTHV">
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \section{Header \texttt{ComplexFFT.h}}
 * \label{s:ComplexFFT.h}
 *
 * Performs complex-to-complex FFTs.
 *
 * \subsection*{Synopsis}
 * \begin{verbatim}
 * #include <lal/ComplexFFT.h>
 * \end{verbatim}
 *
 * Perform complex-to-complex fast Fourier transforms of vectors using the
 * package FFTW~\cite{fj:1998}.
 *
 **** </lalLaTeX> */

#ifndef _COMPLEXFFT_H
#define _COMPLEXFFT_H

#include <lal/LALDatatypes.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif

NRCSID( COMPLEXFFTH, "$Id$" );

/**** <lalLaTeX>
 * \subsection*{Error conditions}
 **** </lalLaTeX> */
/**** <lalErrTable> */

#define COMPLEXFFTH_ENULL 1
#define COMPLEXFFTH_ENNUL 2
#define COMPLEXFFTH_ESIZE 4
#define COMPLEXFFTH_ESZMM 8
#define COMPLEXFFTH_ESLEN 16
#define COMPLEXFFTH_ESAME 32
#define COMPLEXFFTH_EALOC 64
#define COMPLEXFFTH_EFFTW 128
#define COMPLEXFFTH_ESNGL 256
#define COMPLEXFFTH_EINTL 512
#define COMPLEXFFTH_ESIGN 1024

#define COMPLEXFFTH_MSGENULL "Null pointer"
#define COMPLEXFFTH_MSGENNUL "Non-null pointer"
#define COMPLEXFFTH_MSGESIZE "Invalid input size"
#define COMPLEXFFTH_MSGESZMM "Size mismatch"
#define COMPLEXFFTH_MSGESLEN "Invalid/mismatched sequence lengths"
#define COMPLEXFFTH_MSGESAME "Input/Output data vectors are the same"
#define COMPLEXFFTH_MSGEFFTW "Error in FFTW"
#define COMPLEXFFTH_MSGEALOC "Memory allocation failed"
#define COMPLEXFFTH_MSGESNGL "FFTW library is not single-precision"
#define COMPLEXFFTH_MSGEINTL "Error in Intel FFT library"
#define COMPLEXFFTH_MSGESIGN "Unknown sign of transform in plan"

/**** </lalErrTable> */
/**** <lalLaTeX>
 *
 * \subsection*{Structures}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct tagCOMPLEX8FFTPlan COMPLEX8FFTPlan;
typedef struct tagCOMPLEX16FFTPlan COMPLEX16FFTPlan;
#define tagComplexFFTPlan tagCOMPLEX8FFTPlan
#define ComplexFFTPlan COMPLEX8FFTPlan
/**** </lalVerbatim> */
/**** <lalLaTeX>
 *
 * This structure contains the parameters necessary for performing an FFT of a
 * given size and direction.  The contents should not be manually adjusted.
 *
 * \newpage\input{ComplexFFTC}
 * \newpage\input{ComplexFFTTestC}
 **** </lalLaTeX> */

/*
 *
 * XLAL COMPLEX8 functions
 *
 */

COMPLEX8FFTPlan * XLALCreateCOMPLEX8FFTPlan( UINT4 size, int fwdflg, int measurelvl );
COMPLEX8FFTPlan * XLALCreateForwardCOMPLEX8FFTPlan( UINT4 size, int measurelvl );
COMPLEX8FFTPlan * XLALCreateReverseCOMPLEX8FFTPlan( UINT4 size, int measurelvl );
void XLALDestroyCOMPLEX8FFTPlan( COMPLEX8FFTPlan *plan );
int XLALCOMPLEX8VectorFFT( COMPLEX8Vector *output, COMPLEX8Vector *input,
    const COMPLEX8FFTPlan *plan );

/*
 *
 * XLAL COMPLEX16 functions
 *
 */

COMPLEX16FFTPlan * XLALCreateCOMPLEX16FFTPlan( UINT4 size, int fwdflg, int measurelvl );
COMPLEX16FFTPlan * XLALCreateForwardCOMPLEX16FFTPlan( UINT4 size, int measurelvl );
COMPLEX16FFTPlan * XLALCreateReverseCOMPLEX16FFTPlan( UINT4 size, int measurelvl );
void XLALDestroyCOMPLEX16FFTPlan( COMPLEX16FFTPlan *plan );
int XLALCOMPLEX16VectorFFT( COMPLEX16Vector *output, COMPLEX16Vector *input,
    const COMPLEX16FFTPlan *plan );

/*
 *
 * LAL COMPLEX8 functions
 *
 */

void
LALCreateForwardCOMPLEX8FFTPlan(
    LALStatus       *status,
    COMPLEX8FFTPlan **plan,
    UINT4            size,
    INT4             measure
    );
#define LALCreateForwardComplexFFTPlan LALCreateForwardCOMPLEX8FFTPlan

void
LALCreateReverseCOMPLEX8FFTPlan(
    LALStatus       *status,
    COMPLEX8FFTPlan **plan,
    UINT4            size,
    INT4             measure
    );
#define LALCreateReverseComplexFFTPlan LALCreateReverseCOMPLEX8FFTPlan

void
LALDestroyCOMPLEX8FFTPlan (
    LALStatus       *status,
    COMPLEX8FFTPlan **plan
    );
#define LALDestroyComplexFFTPlan LALDestroyCOMPLEX8FFTPlan

void
LALCOMPLEX8VectorFFT (
    LALStatus      *status,
    COMPLEX8Vector *output,
    COMPLEX8Vector *input,
    COMPLEX8FFTPlan *plan
    );

/*
 *
 * LAL COMPLEX8 functions
 *
 */

void
LALCreateForwardCOMPLEX16FFTPlan(
    LALStatus       *status,
    COMPLEX16FFTPlan **plan,
    UINT4            size,
    INT4             measure
    );

void
LALCreateReverseCOMPLEX16FFTPlan(
    LALStatus       *status,
    COMPLEX16FFTPlan **plan,
    UINT4            size,
    INT4             measure
    );

void
LALDestroyCOMPLEX16FFTPlan (
    LALStatus       *status,
    COMPLEX16FFTPlan **plan
    );

void
LALCOMPLEX16VectorFFT (
    LALStatus      *status,
    COMPLEX16Vector *output,
    COMPLEX16Vector *input,
    COMPLEX16FFTPlan *plan
    );


#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _COMPLEXFFT_H */

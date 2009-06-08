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

#if 0 /* autodoc block */

<lalVerbatim file="InterpolateHV">
$Id$
</lalVerbatim>

<lalLaTeX>

\section{Header \texttt{Interpolate.h}}
\label{s:Interpolate.h}

Generates random numbers.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/Interpolate.h>
\end{verbatim}

\noindent This header covers the routines for interpolation.

</lalLaTeX>

#endif /* autodoc block */

#ifndef _INTERPOLATE_H
#define _INTERPOLATE_H

#include <lal/LALDatatypes.h>

#ifdef __cplusplus
extern "C" {
#endif

NRCSID (INTERPOLATEH, "$Id: Interpolate.h");

#if 0 /* autodoc block */

<lalLaTeX>
\subsection*{Error conditions}
\input{InterpolateHErrTab}
</lalLaTeX>

<lalErrTable file="InterpolateHErrTab">

#endif /* autodoc block */

#define INTERPOLATEH_ENULL 1
#define INTERPOLATEH_ESIZE 2
#define INTERPOLATEH_EZERO 4

#define INTERPOLATEH_MSGENULL "Null pointer"
#define INTERPOLATEH_MSGESIZE "Invalid size"
#define INTERPOLATEH_MSGEZERO "Zero divide"

#if 0 /* autodoc block */

</lalErrTable>

<lalLaTeX>

\subsection*{Structures}

\begin{verbatim}
typedef struct
tagSInterpolateOut
{
  REAL4  y;
  REAL4 dy;
}
SInterpolateOut;

typedef struct
tagDInterpolateOut
{
  REAL8  y;
  REAL8 dy;
}
DInterpolateOut;
\end{verbatim}

These structures contain the output of the interpolation.  The two fields are:
\begin{description}
\item[\texttt{y}] The interpolated value.
\item[\texttt{dy}] The estimated error in the interpolated value.
\end{description}

\begin{verbatim}
typedef struct
tagSInterpolatePar
{
  UINT4  n;
  REAL4 *x;
  REAL4 *y;
}
SInterpolatePar;

typedef struct
tagDInterpolatePar
{
  UINT4  n;
  REAL8 *x;
  REAL8 *y;
}
DInterpolatePar;
\end{verbatim}

These structures contain the interpolation parameters.  These are the arrays
of \verb+n+ domain values \verb+x[0]+\ldots\verb+x[n-1]+ and their
corresponding values \verb+y[0]+\ldots\verb+y[n-1]+.  The fields are:
\begin{description}
\item[\texttt{n}] The number of points in the arrays to use in the
  interpolation.
\item[\texttt{x}] The array of domain values.
\item[\texttt{y}] The array of values to interpolate.
\end{description}

</lalLaTeX>

#endif /* autodoc block */

typedef struct
tagSInterpolateOut
{
  REAL4  y;
  REAL4 dy;
}
SInterpolateOut;

typedef struct
tagDInterpolateOut
{
  REAL8  y;
  REAL8 dy;
}
DInterpolateOut;

typedef struct
tagSInterpolatePar
{
  UINT4  n;
  REAL4 *x;
  REAL4 *y;
}
SInterpolatePar;

typedef struct
tagDInterpolatePar
{
  UINT4  n;
  REAL8 *x;
  REAL8 *y;
}
DInterpolatePar;

#if 0 /* autodoc block */

<lalLaTeX>
\newpage\input{InterpolateC}
</lalLaTeX>

#endif /* autodoc block */

void
LALSPolynomialInterpolation (
    LALStatus          *status,
    SInterpolateOut *output,
    REAL4            target,
    SInterpolatePar *params
    );

void
LALDPolynomialInterpolation (
    LALStatus          *status,
    DInterpolateOut *output,
    REAL8            target,
    DInterpolatePar *params
    );

#if 0 /* autodoc block */

<lalLaTeX>
\newpage\input{InterpolateTestC}
</lalLaTeX>

#endif /* autodoc block */

#ifdef __cplusplus
}
#endif

#endif /* _INTERPOLATE_H */

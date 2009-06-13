/*
 *  Copyright (C) 2005 Badri Krishnan, Alicia Sintes
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

/** \file
 *  \ingroup modukeHoughPulsar
 *  \author Sintes, A. M.
 *  \date $Date$
 *  \brief Conversion from peaks in a spectrum into a partial Hough map derivative
 *
 *  $Id$
 *
 * History:   Created by Sintes June 21, 2001
 *            Modified...
 *
 *

\par Description

The Hough map is an histogram, thus additive. It can be seen as the sum of several
partial Hough maps constructed using just one periodogram, or equivalently, as
the sum of partial Hough map derivatives phmd and then integrating the
result.

A phmd can be represented by a set of borders, here called \it
left and right. They indicate the beginning and the end of the annuli.
The position of the so-called left borders should be marked with $+1$,
and
the position of the right borders should be marked with $-1$ in the {\sc phmd}.
To obtain a partial Hough map, one needs to integrate each row of the {\sc phmd}
from left to right.

The representation of a {\sc phmd} is simplified by considering
pointers to the borders in a pre-calculated look-up-table, plus some
extra information about
their character and edge effects when clipping on a finite patch.

 */


/************************************ <lalVerbatim file="PHMDHV">
Author: Sintes, A.M., Krishnan, B.
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Header \texttt{PHMD.h}}
\label{s:PHMD.h}

Conversion from peaks in a spectrum into a partial Hough map derivative.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection*{Synopsis}

\begin{verbatim}
#include <lal/PHMD.h>
\end{verbatim}

 The Hough map is an histogram, thus additive. It can be seen as the sum of several
 partial Hough maps constructed using just one periodogram, or equivalently, as
 the sum of partial Hough map derivatives ({\sc phmd}) and then integrating the
 result.\\

  A  {\sc phmd} can be represented by a set of borders, here called {\it
 left} and {\it right}. They indicate the beginning and the end of the annuli.
 The position of the so-called left borders should be marked with $+1$,
  and
the position of the right borders should be marked with $-1$ in the {\sc phmd}.
To obtain a partial Hough map, one needs to integrate each row of the {\sc phmd}
from left to right.\\

 The representation of a  {\sc phmd} is simplified by considering
 pointers to the borders in a pre-calculated look-up-table, plus some extra information about
 their character and edge effects when clipping on a finite patch.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection*{Error conditions}
\vspace{0.1in}
\input{PHMDHErrorTable}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection*{Structures}

\begin{verbatim}
struct HOUGHPeakGram
\end{verbatim}
\index{\texttt{HOUGHPeakGram}}

\noindent This structure stores the \lq\lq peak-gram".  The fields are:

\begin{description}
\item[\texttt{INT2    timeIndex}]  The time index of the peak-gram.
\item[\texttt{REAL8   deltaF}]  Frequency resolution: \texttt{df=1/TCOH}.
\item[\texttt{UINT8   fBinIni }]  Frequency index of the first element of the
spectrum covered by this peak-gram. It can be seen as an offset
\item[\texttt{UINT8   fBinFin}]  Frequency index of the last element of the
spectrum covered by this peak-gram.
\item[\texttt{UINT4   length}]   Number of peaks present in the peak-gram.
\item[\texttt{INT4    *peak}]  The peak indices relative to \texttt{fBinIni},
i.e., the zero peak  corresponds to \texttt{fBinIni}.
\end{description}

\begin{verbatim}
struct HOUGHphmd
\end{verbatim}
\index{\texttt{HOUGHphmd}}

\noindent This structure stores a partial Hough map derivative.  The fields are:

\begin{description}
\item[\texttt{UINT8          fBin }]   Frequency bin of this partial map
derivative
\item[\texttt{UINT2          lengthLeft }]  Exact number of {\it Left} borders.
\item[\texttt{UINT2          lengthRight }]  Exact number of {\it Right} borders.
\item[\texttt{UINT2          maxNBorders }]  Maximun number of borders of each type (for
     memory allocation purposes), i.e. length of \verb@*leftBorderP@ and \verb@*rightBorderP@.
\item[\texttt{HOUGHBorder   **leftBorderP }]  Pointers to borders.
\item[\texttt{HOUGHBorder   **rightBorderP }]   Pointers to borders.
\item[\texttt{UINT2          ySide }]  Number of elements of \verb@firstColumn@.
\item[\texttt{UCHAR         *firstColumn }]  First column border,
containing the edge effects  when clipping on a finite patch.
\end{description}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\vfill{\footnotesize\input{PHMDHV}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newpage\input{Peak2PHMDC}

%%%%%%%%%%Test program. %%
\newpage\input{TestPeak2PHMDC}
\newpage\input{TestNDPeak2PHMDC}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

</lalLaTeX> */



/*
 * 4.  Protection against double inclusion (include-loop protection)
 *     Note the naming convention!
 */

#ifndef _PHMD_H
#define _PHMD_H

/*
 * 5. Includes. This header may include others; if so, they go immediately
 *    after include-loop protection. Includes should appear in the following
 *    order:
 *    a. Standard library includes
 *    b. LDAS includes
 *    c. LAL includes
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
# include <stdlib.h>
# include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>

#include <lal/LUT.h>

/*
 *  #include "LALRCSID.h"
 *  not needed, it is included in "LALConstants.h"
 */



/*
 *   Protection against C++ name mangling
 */

#ifdef  __cplusplus
extern "C" {
#endif

/*
 * 6. Assignment of Id string using NRCSID()
 */

NRCSID (PHMDH, "$Id$");

/*
 * 7. Error codes and messages. This must be auto-extracted for
 *    inclusion in the documentation.
 */

/* <lalErrTable file="PHMDHErrorTable"> */

#define PHMDH_ENULL 1
#define PHMDH_ESIZE 2
#define PHMDH_ESZMM 4
#define PHMDH_EINT  6
#define PHMDH_ESAME 8
#define PHMDH_EFREQ 10
#define PHMDH_EVAL 12

#define PHMDH_MSGENULL "Null pointer"
#define PHMDH_MSGESIZE "Invalid input size"
#define PHMDH_MSGESZMM "Size mismatch"
#define PHMDH_MSGEINT  "Invalid interval"
#define PHMDH_MSGESAME "Input/Output data vectors are the same"
#define PHMDH_MSGEFREQ "Invalid frequency"
#define PHMDH_MSGEVAL  "Invalid value"

/* </lalErrTable>  */


/* ******************************************************
 * 8. Macros. But, note that macros are deprecated.
 *    They could be moved to the modules where are needed
 */


/* *******************************************************
 * 9. Constant Declarations. (discouraged)
 */



/* **************************************************************
 * 10. Structure, enum, union, etc., typdefs.
 */

/*  Hough Map derivative pixel type */
/* typedef CHAR  HoughDT; */
/*  typedef INT2  HoughDT; */
  typedef REAL8 HoughDT; /* for weighted hough maps */

/** structure to store peakgrams */
  typedef struct tagHOUGHPeakGram{
    INT2    timeIndex;  /**< time index of the Peakgram */
    REAL8   deltaF;     /**< df=1/TCOH */
    UINT8   fBinIni;    /**< frequency bin of the zero peak (initial offset) */
    UINT8   fBinFin;    /**< maximum frequency bin of the peakgram */
    UINT4   length;     /**< number of peaks present in the peakgram  */
    INT4    *peak;      /**< the peak indexes relative to fBinIni*/
  } HOUGHPeakGram;

  /** partial hough map derivative structure */
typedef struct tagHOUGHphmd{
  UINT8          fBin;  /**< frequency bin of this partial map derivative */
  UINT2          lengthLeft; /**< exact number of Left borders */
  UINT2          lengthRight;/**< exact number of Right borders */
  UINT2          maxNBorders; /**< maximun number of borders of each type (for memory allocation purposes) */
  HOUGHBorder    **leftBorderP; /**< pointer the borders[x] +1 */
  HOUGHBorder    **rightBorderP; /**< pointer the borders[x] -1 */
  UINT2          ySide;  /**< number of elements of firstColumn */
  UCHAR          *firstColumn; /**< border corrections on 1st column */
  HoughDT        weight; /**< weighting factor of phmd -- for improving sensitivity */
} HOUGHphmd;


/*
 * 11. Extern Global variables. (discouraged)
 */


/*
 * 12. Functions Declarations (i.e., prototypes).
 */

void LALHOUGHPeak2PHMD (LALStatus    *status,
			HOUGHphmd    *phmd,
			HOUGHptfLUT  *lut,
			HOUGHPeakGram *pg
			);


#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif     /* Close double-include protection _PHMD_H */

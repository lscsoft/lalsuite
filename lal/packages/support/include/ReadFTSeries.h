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

/************************************ <lalVerbatim file="ReadFTSeriesHV">
Author: Torres, C. W.
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{ReadFTSeries.h}}
\label{s:ReadFTSeries.h}

This is a simple utility to Read time and frequency series into a
file.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/ReadFTSeries.h>
\end{verbatim}

\noindent Provides prototype information for the routines in
\verb+ReadTimeSeries.c+ and \verb+ReadFrequencySeries.c+.

\vfill{\footnotesize\input{ReadFTSeriesHV}}
\newpage\input{ReadFrequencySeriesC}
\newpage\input{ReadTimeSeriesC}
\newpage\input{ReadFTSeriesTestC}

</lalLaTeX> */

#ifndef _READFTSERIES_H
#define _READFTSERIES_H

#include <lal/LALRCSID.h>
#include <lal/LALStdlib.h>

#ifdef  __cplusplus
extern "C" {
#endif

enum { LALSupportUnitTextSize = sizeof("10^-32768 m^-32768/32767 "
				       "kg^-32768/32767 "
				       "s^-32768/32767 A^-32768/32767 "
				       "K^-32768/32767 strain^-32768/32767 "
				       "count^-32768/32767") };

enum { MaxLineLength = LALSupportUnitTextSize + sizeof("Units are ()\n") };

NRCSID( READFTSERIESH, "$Id$" );

void LALReadTimeSeries(LALStatus* status,  REAL4TimeSeries *series , const CHAR *filename );
void LALSReadTimeSeries(LALStatus* status,  REAL4TimeSeries *series , const CHAR *filename );
void LALDReadTimeSeries(LALStatus* status,  REAL8TimeSeries *series , const CHAR *filename );
void LALCReadTimeSeries(LALStatus* status,  COMPLEX8TimeSeries *series , const CHAR *filename );
void LALZReadTimeSeries(LALStatus* status,  COMPLEX16TimeSeries *series , const CHAR *filename );

void LALReadFrequencySeries(LALStatus* status,  REAL4FrequencySeries *series , const CHAR *filename );
void LALSReadFrequencySeries(LALStatus* status,  REAL4FrequencySeries *series , const CHAR *filename );
void LALDReadFrequencySeries(LALStatus* status,  REAL8FrequencySeries *series , const CHAR *filename );
void LALCReadFrequencySeries(LALStatus* status, COMPLEX8FrequencySeries *series , const CHAR *filename );
void LALZReadFrequencySeries(LALStatus* status,  COMPLEX16FrequencySeries *series , const CHAR *filename );


#ifdef  __cplusplus
}
#endif

#endif /* _READSERIES_H */



/****************************** <lalErrTable file="ReadFTSeriesHE"> */

#define  READFTSERIESH_EFILENOTFOUND       1
#define  READFTSERIESH_EPARSE              2

#define  READFTSERIESH_MSGEFILENOTFOUND    "Invalid Filename or File Not Found"
#define  READFTSERIESH_MSGEPARSE           "Error Parsing File"

/****************************** </lalErrTable> */

/*
*  Copyright (C) 2007 Jolien Creighton, John Whelan
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

#ifndef _PRINTFTSERIES_H
#define _PRINTFTSERIES_H

#include <lal/LALStdlib.h>

#ifdef  __cplusplus
extern "C" {
#endif


/**
   \addtogroup PrintFTSeries_h
\author J. T. Whelan <jtwhelan@loyno.edu>

   \brief This is a simple utility to print time and frequency series into a file.

\heading{Synopsis}
\code
#include <lal/PrintFTSeries.h>
\endcode

   Provides prototype information for the routines in \ref PrintTimeSeries_c and \ref PrintFrequencySeries_c.

*/
/*@{*/

/** \defgroup PrintTimeSeries_c Module PrintTimeSeries.c

    Print a \<datatype\>TimeSeries object into a
    file.  For use in non-production and test code only.

    \heading{Description}

    Each member of this family of functions prints the elements of
    \f$\langle\mbox{datatype}\rangle\f$\c TimeSeries into a file.  Note:
    the file name is specified using a character string.  This function is
    for debugging use only: its arguments do not conform to LAL standards
    so it should not be used in any real analysis codes.

    \heading{Notes}

    This function's arguments do not conform to the LAL spec.  For this
    reason it should only be used for debugging purposes in test
    functions, not in any production code.

    Additionally, since printf cannot handle INT8 as integers, the
    functions <tt>LALI8PrintTimeSeries()</tt> and
    <tt>LALU8PrintTimeSeries()</tt> use a typecast to REAL8 and are thus
    only valid for numbers between around \f$-10^{15}\f$ and \f$10^{15}\f$.

    The first five lines of the file are a header containing:
    <ol>
    <li> the name of the series</li>
    <li> the starting epoch </li>
    <li> the units expressed in terms of the basic SI units</li>
    <li> column labels</li>
    </ol>
    after which come the data, one per line.

    The output format is two or three tab-separated columns: the first
    column is the time corresponding to the row in question, in seconds
    after the series' starting epoch; for real and integer time
    series, the second column is the value of the series; for complex time
    series, the second column is the real part and the third the imaginary
    part of the value.
*/
/*@{*/
void LALI2PrintTimeSeries( INT2TimeSeries *series , const CHAR *filename );
void LALI4PrintTimeSeries( INT4TimeSeries *series , const CHAR *filename );
void LALI8PrintTimeSeries( INT8TimeSeries *series , const CHAR *filename );
void LALU2PrintTimeSeries( UINT2TimeSeries *series , const CHAR *filename );
void LALU4PrintTimeSeries( UINT4TimeSeries *series , const CHAR *filename );
void LALU8PrintTimeSeries( UINT8TimeSeries *series , const CHAR *filename );
void LALPrintTimeSeries( REAL4TimeSeries *series , const CHAR *filename );
void LALSPrintTimeSeries( REAL4TimeSeries *series , const CHAR *filename );
void LALDPrintTimeSeries( REAL8TimeSeries *series , const CHAR *filename );
void LALCPrintTimeSeries( COMPLEX8TimeSeries *series , const CHAR *filename );
void LALZPrintTimeSeries( COMPLEX16TimeSeries *series , const CHAR *filename );
/*@}*/


/** \defgroup PrintFrequencySeries_c Module PrintFrequencySeries.c

Print a \<datatype\>FrequencySeries object into a
file.  For use in non-production and test code only.

\heading{Description}

Each member of this family of functions prints the elements of
\<datatype\>FrequencySeries into a file.
Note: the file name is specified using a character string.  This
function is for debugging use only: its arguments do not conform to
LAL standards so it should not be used in any real analysis codes.

\heading{Notes}

This function's arguments do not conform to the LAL spec.  For this
reason it should only be used for debugging purposes in test
functions, not in any production code.

Additionally, since printf cannot handle INT8 as integers, the
functions <tt>LALI8PrintFrequencySeries()</tt> and
<tt>LALU8PrintFrequencySeries()</tt> use a typecast to REAL8 and are
thus only valid for numbers between around \f$-10^{15}\f$ and \f$10^{15}\f$.

The first four lines of the file are a header containing:
<ol>
<li> the name of the series</li>
<li> heterodyning information, if any</li>
<li> the starting epoch, relative to the GPS reference epoch (1980 January 6)</li>
<li> the units expressed in terms of the basic SI units</li>
<li> column labels</li>
</ol>
after which come the data, one per line.

The output format is two or three tab-separated columns: the first
column is the frequency in hertz corresponding to the row in question;
for real and integer frequency series, the second column is the value
of the series; for complex frequency series, the second column is the
real part and the third the imaginary part of the value.

Note that the frequency given is the physical frequency.  In the case
of a heterodyned frequency series, this is the heterodyning frequency
plus the frequency offset.  A frequency series of length \f$[N]\f$ is
assumed to be packed so that the 0th element corresponds to zero
frequency offset, elements 1 through \f$[N/2]\f$ to positive frequency
offsets (in ascending order), and elements \f$N-[N/2]\f$ to \f$N-1\f$ to
negative frequency offsets (also in ascending order, so that the
frequency corresponding to the \f$N-1\f$st element is just below that of
the 0th element).  If \f$N\f$ is even, the element in position \f$N/2\f$ is
assumed to correspond both the maximum poitive and negative frequency
offset.
*/
/*@{*/
void LALI2PrintFrequencySeries( INT2FrequencySeries *series , const CHAR *filename );
void LALI4PrintFrequencySeries( INT4FrequencySeries *series , const CHAR *filename );
void LALI8PrintFrequencySeries( INT8FrequencySeries *series , const CHAR *filename );
void LALU2PrintFrequencySeries( UINT2FrequencySeries *series , const CHAR *filename );
void LALU4PrintFrequencySeries( UINT4FrequencySeries *series , const CHAR *filename );
void LALU8PrintFrequencySeries( UINT8FrequencySeries *series , const CHAR *filename );
void LALPrintFrequencySeries( REAL4FrequencySeries *series , const CHAR *filename );
void LALSPrintFrequencySeries( REAL4FrequencySeries *series , const CHAR *filename );
void LALDPrintFrequencySeries( REAL8FrequencySeries *series , const CHAR *filename );
void LALCPrintFrequencySeries( COMPLEX8FrequencySeries *series , const CHAR *filename );
void LALZPrintFrequencySeries( COMPLEX16FrequencySeries *series , const CHAR *filename );
/*@}*/


/*@}*/

#ifdef  __cplusplus
}
#endif

#endif /* _PRINTFTSERIES_H */

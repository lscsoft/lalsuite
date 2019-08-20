/*
*  Copyright (C) 2007 Alexander Dietz, Drew Keppel, Duncan Brown, Eirini Messaritaki, Gareth Jones, George Birthisel, Jolien Creighton, Kipp Cannon, Lisa M. Goggin, Patrick Brady, Peter Shawhan, Saikat Ray-Majumder, Stephen Fairhurst, Xavier Siemens, Craig Robinson , Thomas Cokelaer
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

/*-----------------------------------------------------------------------
 *
 * File Name: LIGOMetadataUtils.h
 *
 * Author: Brown, D. A. and Fairhurst, S.
 *
 *-----------------------------------------------------------------------
 */

/**
 * \author Brown, D. A. and Fairhurst, S.
 * \file
 * \ingroup lalmetaio_general
 * \brief Provides functions for manipulating the LAL structures that correspond
 * to the LIGO metadata database tables defined in \ref LIGOMetadataTables.h.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/LIGOMetadataUtils.h>
 * \endcode
 *
 * This header provides prototypes for routines that perform processing
 * on the LAL structures that correspond to the LIGO metadata database tables
 * defined in \ref LIGOMetadataTables.h, such as sorting and eliminating
 * duplictaes. The functions specific to a particular metadata table (e.g.
 * \c sngl_inspiral, \c sngl_burst, etc.) are all prototyped in
 * this header.
 *
 * ### Types ###
 *
 * None.
 *
 */

#ifndef _LIGOMETADATAUTILS_H
#define _LIGOMETADATAUTILS_H

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

#include <lal/LIGOMetadataTables.h>
#include <lal/LALDetectors.h>

/**
 * The \c LALPlaygroundDataMask contains an enum type for describing the
 * subset of data to be used, \c playground_only, \c exclude_play and \c all_data.
 */
typedef enum tagLALPlaygroundDataMask
{
  unspecified_data_type,
  playground_only,
  exclude_play,
  all_data
}
LALPlaygroundDataMask;


/*
 *
 * general manipulation functions
 *
 */

ProcessTable *XLALCreateProcessTableRow(void);
#ifndef SWIG   // exclude from SWIG interface
void XLALDestroyProcessTableRow(ProcessTable *);
#endif   // SWIG
void XLALDestroyProcessTable(ProcessTable *);
long XLALProcessTableGetNextID(ProcessTable *);
int XLALPopulateProcessTable(
	ProcessTable *ptable,
	const char *program_name,
	const char *cvs_revision,
	const char *cvs_source,
	const char *cvs_date,
	long process_id
);

ProcessParamsTable *XLALCreateProcessParamsTableRow(const ProcessTable *);
#ifndef SWIG   // exclude from SWIG interface
void XLALDestroyProcessParamsTableRow(ProcessParamsTable *);
#endif   // SWIG
void XLALDestroyProcessParamsTable(ProcessParamsTable *);

TimeSlide *XLALCreateTimeSlide(void);
#ifndef SWIG   // exclude from SWIG interface
void XLALDestroyTimeSlide(TimeSlide *);
#endif   // SWIG
void XLALDestroyTimeSlideTable(TimeSlide *);
const TimeSlide *XLALTimeSlideConstGetByIDAndInstrument(const TimeSlide *, long, const char *);
TimeSlide *XLALTimeSlideGetByIDAndInstrument(TimeSlide *, long, const char *);

SearchSummaryTable *XLALCreateSearchSummaryTableRow(const ProcessTable *);
#ifndef SWIG   // exclude from SWIG interface
void XLALDestroySearchSummaryTableRow(SearchSummaryTable *);
#endif   // SWIG
void XLALDestroySearchSummaryTable(SearchSummaryTable *);

SegmentTable *XLALCreateSegmentTableRow(const ProcessTable *process);
#ifndef SWIG   // exclude from SWIG interface
void XLALDestroySegmentTableRow(SegmentTable *row);
#endif   // SWIG
void XLALDestroySegmentTable(SegmentTable *head);

int
XLALCountProcessTable(
    ProcessTable *head
    );

int
XLALCountProcessParamsTable(
    ProcessParamsTable *head
    );

int
XLALIFONumber(
    const char *ifo
    );

void
XLALReturnIFO(
    char                *ifo,
    InterferometerNumber IFONumber
    );

void
XLALReturnDetector(
    LALDetector           *det,
    InterferometerNumber   IFONumber
    );


int
XLALCompareSearchSummaryByOutTime (
    const void *a,
    const void *b
    );

int
XLALTimeSortSearchSummary(
    SearchSummaryTable  **summHead,
    int(*comparfunc)    (const void *, const void *)
    );

SearchSummaryTable *
XLALIfoScanSearchSummary(
    SearchSummaryTable         *input,
    CHAR                       *ifos
    );

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LIGOMETADATAUTILS_H */

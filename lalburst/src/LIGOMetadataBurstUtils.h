/*
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

#ifndef _LIGOMETADATABURSTUTILS_H
#define _LIGOMETADATABURSTUTILS_H

#ifdef  __cplusplus
extern "C" {
#endif

#include <lal/LIGOMetadataTables.h>

/*
 *
 * burst specific functions
 *
 */


void
XLALDestroySnglBurst(
	SnglBurst *event
);

void
XLALDestroySnglBurstTable(
	SnglBurst *head
);

int
XLALSnglBurstTableLength(
	SnglBurst *head
);

long
XLALSnglBurstAssignIDs(
	SnglBurst *head,
	long process_id,
	long event_id
);

SnglBurst **
XLALSortSnglBurst(
	SnglBurst **head,
	int (*comparefunc)(const SnglBurst * const *, const SnglBurst * const *)
);

int
XLALCompareSnglBurstByStartTime(
	const SnglBurst * const *a,
	const SnglBurst * const *b
);

int
XLALCompareSnglBurstByExactPeakTime(
	const SnglBurst * const *a,
	const SnglBurst * const *b
);

int
XLALCompareSnglBurstBySNR(
	const SnglBurst * const *a,
	const SnglBurst * const *b
);

int
XLALCompareSnglBurstByPeakTimeAndSNR(
	const SnglBurst * const *a,
	const SnglBurst * const *b
);

int
XLALCompareSnglBurstByLowFreq(
	const SnglBurst * const *a,
	const SnglBurst * const *b
);

int
XLALCompareSnglBurstByStartTimeAndLowFreq(
	const SnglBurst * const *a,
	const SnglBurst * const *b
);

void
XLALStringBurstCluster(
	SnglBurst *a,
	const SnglBurst *b
);

void
XLALClusterSnglBurstTable(
	SnglBurst  **list,
	int (*bailoutfunc)(const SnglBurst * const *, const SnglBurst * const *),
	int (*testfunc)(const SnglBurst * const *, const SnglBurst * const *),
	void (*clusterfunc)(SnglBurst *, const SnglBurst *)
);

SnglBurst *
XLALCreateSnglBurst(
	void
);

SimBurst *
XLALCreateSimBurst(
	void
);

void
XLALDestroySimBurst(
	SimBurst *row
);

void
XLALDestroySimBurstTable(
	SimBurst *head
);

int XLALCompareSimBurstByGeocentTimeGPS(
	const SimBurst * const *a,
	const SimBurst * const *b
);

int XLALSimBurstTableLength(
	SimBurst *head
);

SimBurst **XLALSortSimBurst(
	SimBurst **head,
	int (*comparefunc)(const SimBurst * const *, const SimBurst * const *)
);

long
XLALSimBurstAssignIDs(
	SimBurst *head,
	long process_id,
	long event_id
);

#ifdef  __cplusplus
}
#endif

#endif /* _LIGOMETADATABURSTUTILS_H */

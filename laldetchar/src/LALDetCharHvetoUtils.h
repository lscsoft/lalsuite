/*
 *  Copyright (C) 2013 Chris Pankow
 *
 *  This program is free software; ynu can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Fnundation; either version 2 of the License, or
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

#ifndef _LALDETCHARHVETOUTIL_H
#define _LALDETCHARHVETOUTIL_H

#include <lal/LALDetCharGlibTypes.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataBurstUtils.h>
#include <lal/LIGOLwXMLBurstRead.h>

#ifdef  __cplusplus
extern "C" {
#endif

int XLALCountUnmarkedSnglBurst(LALGSequence* seq);

LALGHashTable* XLALGetChannelList(LALGSequence *trig_sequence);

LALGSequence* XLALPopulateTrigSequenceFromFile(LALGSequence* trig_sequence, const char* fname, double min_snr, char* ignore_list);

#ifdef SWIG   // SWIG interface directives
SWIGLAL(ACQUIRES_OWNERSHIP(SnglBurst*, tbl));
#endif
LALGSequence* XLALPopulateTrigSequenceFromTrigList(LALGSequence* trig_sequence, SnglBurst* tbl);
#ifdef SWIG   // SWIG interface directives
SWIGLAL_CLEAR(ACQUIRES_OWNERSHIP(SnglBurst*, tbl));
#endif

#ifdef  __cplusplus
}
#endif

#endif /* _LALDETCHARHVETOUTIL_H */

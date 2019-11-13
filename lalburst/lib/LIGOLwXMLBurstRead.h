/*
 * LIGOLwXMLBurstRead.h
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with with program; see the file COPYING. If not, write to the Free
 * Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 */

#ifndef _LIGOLWXMLBURSTREAD_H_
#define _LIGOLWXMLBURSTREAD_H_

#include <lal/Date.h>
#include <lal/LIGOMetadataTables.h>

#ifdef  __cplusplus
extern "C" {
#endif

/*
 * function prototypes
 */

SnglBurst *XLALSnglBurstTableFromLIGOLw(
    const char *filename
);

SimBurst *XLALSimBurstTableFromLIGOLw(
    const char *filename,
    const LIGOTimeGPS *start,
    const LIGOTimeGPS *end
);

#ifdef  __cplusplus
}
#endif

#endif /* _LIGOLWXMLBURSTREAD_H_ */

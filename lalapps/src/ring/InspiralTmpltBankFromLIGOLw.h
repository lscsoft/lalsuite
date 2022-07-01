/*
 * LIGOLwXMLInspiralRead.h
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
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301  USA
 */

#ifndef _INSPIRALTMPLTBANKFROMLIGOLW_H
#define _INSPIRALTMPLTBANKFROMLIGOLW_H

#include <lal/LALDatatypes.h>
#include <lal/LALInspiral.h>
#include <lal/LIGOMetadataTables.h>

#ifdef  __cplusplus
extern "C" {
#endif

/*
 * function prototypes
 */

int
InspiralTmpltBankFromLIGOLw (
    InspiralTemplate   **bankHead,
    const CHAR         *fileName,
    INT4                startTmplt,
    INT4                stopTmplt
    );

#ifdef  __cplusplus
}
#endif

#endif /* _INSPIRALTMPLTBANKFROMLIGOLW_H */

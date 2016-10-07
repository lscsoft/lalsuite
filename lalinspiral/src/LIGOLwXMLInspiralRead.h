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
 * Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 */

#ifndef _LIGOLWXMLINSPIRALREAD_H
#define _LIGOLWXMLINSPIRALREAD_H

#include <lal/LIGOLwXMLRead.h>
#include <lal/LALInspiral.h>
#include <lal/LALInspiralBank.h>

#ifdef  __cplusplus
extern "C" {
#endif

/*
 * function prototypes
 */

MultiInspiralTable* XLALMultiInspiralTableFromLIGOLw (
    CHAR               *fileName
    );

int
LALSnglInspiralTableFromLIGOLw (
    SnglInspiralTable **eventHead,
    const CHAR         *fileName,
    INT4                startEvent,
    INT4                stopEvent
    );

int
InspiralTmpltBankFromLIGOLw (
    InspiralTemplate   **bankHead,
    const CHAR         *fileName,
    INT4                startTmplt,
    INT4                stopTmplt
    );

int
SimInspiralTableFromLIGOLw (
    SimInspiralTable   **simHead,
    const CHAR          *fileName,
    INT4                 startTime,
    INT4                 endTime
    );

int
SummValueTableFromLIGOLw (
    SummValueTable **sumHead,
    CHAR           *fileName
    );

int
LALExtTriggerTableFromLIGOLw (
    ExtTriggerTable   **eventHead,
    CHAR               *fileName,
    INT4                startEvent,
    INT4                stopEvent
    );

int
XLALReadSummValueFile (
    SummValueTable **summValueList,
    CHAR                  *fileName
    );
int
XLALReadInspiralTriggerFile (
    SnglInspiralTable    **inspiralEventList,
    SnglInspiralTable    **lastTrigger,
    SearchSummaryTable   **searchSummList,
    SearchSummvarsTable  **inputFileList,
    CHAR                  *fileName
    );

void
XLALCleanSummValueTable (
    SummValueTable **summValueList
    );

int
LALMultiInspiralTableFromLIGOLw (
    MultiInspiralTable **eventHead,
    CHAR                *fileName
    );

int
XLALReadMultiInspiralTriggerFile (
    MultiInspiralTable    **inspiralEventList,
    MultiInspiralTable    **lastTrigger,
    SearchSummaryTable   **searchSummList,
    SearchSummvarsTable  **inputFileList,
    CHAR                  *fileName
    );

#ifdef  __cplusplus
}
#endif

#endif /* _LIGOLWXMLINSPIRALREAD_H */

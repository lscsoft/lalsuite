/*----------------------------------------------------------------------- 
 * 
 * File Name: ligolw_tmpltbank.h
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#ifndef LIGOLW_TMPLTBANK_H_
#define LIGOLW_TMPLTBANK_H_

#include <stdio.h>
#include <stdlib.h>
#include <metaio.h>
#include <lal/LALInspiral.h>
#include <lal/LALInspiralBank.h>
#include <lal/LIGOLwXMLHeaders.h>

int
InspiralTmpltBankFromLIGOLw (
    InspiralTemplate   **bankHead,
    CHAR                *fileName,
    INT4                startTmplt,
    INT4                stopTmplt
    );

int
InspiralTmpltBankToLIGOLw (
    InspiralTemplate   **bankHead,
    CHAR                *fileName,
    INT4                startTmplt,
    INT4                stopTmplt
    );

#endif

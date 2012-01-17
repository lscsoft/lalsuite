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

/*-----------------------------------------------------------------------
 *
 * File Name: olapredfcn.h
 *
 * Author: John T. Whelan
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <lalapps.h>
#include <lal/LALStdlib.h>
#include <lal/StochasticCrossCorrelation.h>
#include <lal/DetectorSite.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/PrintFTSeries.h>

#include "lalapps.h"

#define OLAPREDFCNH_TRUE     1
#define OLAPREDFCNH_FALSE    0
#define OLAPREDFCNH_OOR    1e+5

/**\name Error Codes */ /*@{*/
#define OLAPREDFCNH_ENOM 0
#define OLAPREDFCNH_EARG 1
#define OLAPREDFCNH_ESUB 2
#define OLAPREDFCNH_MSGENOM "Nominal exit"
#define OLAPREDFCNH_MSGEARG "Invalid command-line arguments"
#define OLAPREDFCNH_MSGECHK "LAL Subroutine Returned Error"
/*@}*/

void
olapredfcn_usage (
    const char *program,
    int         exitcode
    );

void
olapredfcn_parse_options (
    int         argc,
    char       *argv[]
    );

void
olapredfcn_print_options (
    void
    );

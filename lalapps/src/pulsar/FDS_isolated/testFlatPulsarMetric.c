/*
 * Copyright (C) 2005 Reinhard Prix
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

/**
 * \author Reinhard Prix
 * \date 2005
 * \file 
 * \ingroup flatPulsarMetric
 * \brief test-functions for the module FlatPulsarMetric.
 *
 * $Id$
 *
 */


/*---------- INCLUDES ----------*/
#include <math.h>

#include "FlatPulsarMetric.h"

#include <lalapps.h>

RCSID ("$Id$");

/*---------- DEFINES ----------*/
#define TRUE (1==1)
#define FALSE (1==0)


/* prototypes */
int main(void);

/*--------------------------------------------------*/
/* Test function(s) */
/*--------------------------------------------------*/
int main (void)
{
  LALStatus status = blank_status;
  REAL8Vector *metric = NULL;
  LIGOTimeGPS startTime = {714180733, 0};
  REAL8 duration = 180000;
  LALDetector site = lalCachedDetectors[LALDetectorIndexGEO600DIFF];

  lalDebugLevel = 1;

  LAL_CALL ( LALFlatPulsarMetric( &status, &metric, startTime, duration, &site ), &status);

  return 0;
}

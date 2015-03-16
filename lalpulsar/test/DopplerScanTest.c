/*
 * Copyright (C) 2015 Reinhard Prix
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

#include <math.h>

#include <lal/LALMalloc.h>
#include <lal/DopplerScan.h>
#include <lal/DopplerFullScan.h>

/**
 * \author Reinhard Prix
 * \file
 * \brief Tests for exported functions in Doppler[Full]Scan
 *
 */

// ---------- defines --------------------

// ---------- Macros --------------------
// ---------- global variables --------------------

// ---------- local prototypes
static int test_XLALParseSkyRegionString ( void );
static int compareSkyPosition ( SkyPosition *pos1, SkyPosition *pos2, REAL8 tol );
// ---------- function definitions --------------------

/**
 * MAIN function: calls a number of unit-tests
 */
int main( void )
{

  XLAL_CHECK ( test_XLALParseSkyRegionString() == XLAL_SUCCESS, XLAL_EFUNC );

  /* check for memory leaks */
  LALCheckMemoryLeaks();

  /* all tests passed */
  return XLAL_SUCCESS;

} /* main() */


/**
 * Unit test for XLALParseSkyRegionString()
 */
static int
test_XLALParseSkyRegionString ( void )
{

  const char *skyRegionString1 = "(10:25:1.1, -5:30:44.123), (1:2:3.45, 2:33:15), (5:0:0, 0:0:0 )";
  const char *skyRegionString2 = "(2.72715695049852,-0.0962070232002041),(0.270776925139095,0.0445786179780218),(1.30899693899575,0)";

  SkyRegion XLAL_INIT_DECL(skyRegion1);
  SkyRegion XLAL_INIT_DECL(skyRegion2);

  XLAL_CHECK ( XLALParseSkyRegionString ( &skyRegion1, skyRegionString1 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALParseSkyRegionString ( &skyRegion2, skyRegionString2 ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK ( skyRegion1.numVertices == skyRegion2.numVertices, XLAL_EFAILED );

  REAL8 tol = 5e-15;
  XLAL_CHECK ( compareSkyPosition ( &skyRegion1.lowerLeft,  &skyRegion2.lowerLeft, tol ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( compareSkyPosition ( &skyRegion1.upperRight, &skyRegion2.upperRight, tol ) == XLAL_SUCCESS, XLAL_EFUNC );
  for ( UINT4 i = 0; i < skyRegion1.numVertices; i ++ )
    {
      XLAL_CHECK ( compareSkyPosition ( &skyRegion1.vertices[i], &skyRegion2.vertices[i], tol ) == XLAL_SUCCESS, XLAL_EFUNC );
    }

  XLALFree ( skyRegion1.vertices );
  XLALFree ( skyRegion2.vertices );

  return XLAL_SUCCESS;

} // test_XLALParseSkyRegionString()

static int
compareSkyPosition ( SkyPosition *pos1, SkyPosition *pos2, REAL8 tol )
{
  XLAL_CHECK ( (pos1 != NULL) && (pos2 != NULL), XLAL_EINVAL );

  REAL8 diff;
  XLAL_CHECK ( pos1->system == pos2->system, XLAL_EINVAL );
  XLAL_CHECK ( (diff=fabs ( pos1->longitude - pos2->longitude )) < tol, XLAL_ETOL, "Longitude %.16g differs from %.16g by %.2g exceeding tolerance %.2g\n",
               pos1->longitude, pos2->longitude, diff, tol );
  XLAL_CHECK ( (diff=fabs ( pos1->latitude  - pos2->latitude )) < tol, XLAL_ETOL,  "Latitude %.16g differs from %.16g by %.2g exceeding tolerance %.2g\n",
               pos1->latitude, pos2->latitude, diff, tol );

  return XLAL_SUCCESS;
} // compareSkyPosition()

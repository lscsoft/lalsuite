/*
 *  Copyright (C) 2017 Karl Wette
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

#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <lal/LALStdio.h>
#include <lal/LALBitset.h>

int main( void )
{

  /* Create hash table */
  LALBitset *bs = XLALBitsetCreate();
  XLAL_CHECK_MAIN( bs != NULL, XLAL_EFUNC );

  /* Create some random bits */
  BOOLEAN XLAL_INIT_DECL( bits, [4096] );
  gsl_rng *r = gsl_rng_alloc( gsl_rng_mt19937 );
  XLAL_CHECK_MAIN( r != NULL, XLAL_ESYS );
  int nbits = 0;
  for ( size_t n = 0; n < XLAL_NUM_ELEM( bits ); ++n ) {
    bits[n] = ( gsl_rng_uniform( r ) > 0.44 );
    nbits += bits[n] ? 1 : 0;
  }

  /* Create random index offset into bitset */
  const UINT8 n0 = gsl_rng_get( r ) % 65536;

  /* Print information */
  printf("nbits = %i, n0 = %"LAL_UINT8_FORMAT"\n", nbits, n0);

  /* Set bits */
  for ( size_t n = 0; n < XLAL_NUM_ELEM( bits ); ++n ) {
    XLAL_CHECK_MAIN( XLALBitsetSet( bs, n0 + n, bits[n] ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  /* Get bits */
  for ( size_t n = 0; n < XLAL_NUM_ELEM( bits ); ++n ) {
    BOOLEAN is_set = 0;
    XLAL_CHECK_MAIN( XLALBitsetGet( bs, n0 + n, &is_set ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK_MAIN( !is_set == !bits[n], XLAL_EFAILED, "Inconsistent bit at index %"LAL_UINT8_FORMAT": LALBitset=%i, reference=%i", n0 + n, is_set, bits[n] );
  }

  /* Clear bitset */
  XLAL_CHECK_MAIN( XLALBitsetClear( bs ) == XLAL_SUCCESS, XLAL_EFUNC );
  for ( size_t n = 0; n < XLAL_NUM_ELEM( bits ); ++n ) {
    BOOLEAN is_set = 0;
    XLAL_CHECK_MAIN( XLALBitsetGet( bs, n0 + n, &is_set ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK_MAIN( !is_set, XLAL_EFAILED, "Bit still set at index %zu", n0 + n );
  }

  /* Cleanup */
  gsl_rng_free( r );
  XLALBitsetDestroy( bs );

  /* Check for memory leaks */
  LALCheckMemoryLeaks();

  return EXIT_SUCCESS;

}

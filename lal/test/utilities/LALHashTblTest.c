/*
 *  Copyright (C) 2016 Karl Wette
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
#include <stdint.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <lal/LALHashTbl.h>

typedef struct {
  int key;
  int value;
} elem;

static void *new_elem( int key, int value )
{
  elem e = { .key = key, .value = value };
  return memcpy( XLALMalloc( sizeof( e ) ), &e, sizeof( e ) );
}

static UINT8 hash_elem( const void *x )
{
  const elem *ex = ( const elem * ) x;
  UINT2 hval = 0;
  XLALPearsonHash( &hval, sizeof( hval ), &ex->key, sizeof( ex->key ) );
  return hval;
}

static int cmp_elem( const void *x, const void *y )
{
  const elem *ex = ( const elem * ) x;
  const elem *ey = ( const elem * ) y;
  return ex->key - ey->key;
}

int main( void )
{

  /* Create hash table */
  LALHashTbl *ht = XLALHashTblCreate( XLALFree, hash_elem, cmp_elem );
  XLAL_CHECK_MAIN( ht != NULL, XLAL_EFUNC );
  XLAL_CHECK_MAIN( XLALHashTblSize( ht ) == 0, XLAL_EFAILED );

  /* Repeat hash table test a few times */
  for ( int n = 0; n < 4; ++n ) {

    /* Add 100 elements with keys in 100*n + [0,99] to table in a random order */
    {
      gsl_rng *r = gsl_rng_alloc( gsl_rng_mt19937 );
      XLAL_CHECK_MAIN( r != NULL, XLAL_ESYS );
      gsl_permutation *p = gsl_permutation_calloc( 100 );
      XLAL_CHECK_MAIN( p != NULL, XLAL_ESYS );
      gsl_ran_shuffle( r, p->data, 100, sizeof( size_t ) );
      for ( int i = 0; i < 100; ++i ) {
        int key = 100*n + gsl_permutation_get( p, i );
        XLAL_CHECK_MAIN( XLALHashTblAdd( ht, new_elem( key, 3*key - n ) ) == XLAL_SUCCESS, XLAL_EFUNC );
        XLAL_CHECK_MAIN( XLALHashTblSize( ht ) == 100*n + i + 1, XLAL_EFAILED );
      }
      gsl_rng_free( r );
      gsl_permutation_free( p );
    }

    /* Try finding all 100 elements by key */
    for ( int i = 0; i < 100; ++i ) {
      elem x = { .key = 100*n + i };
      const elem *y;
      XLAL_CHECK_MAIN( XLALHashTblFind( ht, &x, ( const void ** ) &y ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( y != NULL, XLAL_EFAILED );
      XLAL_CHECK_MAIN( y->value == 3*y->key - n, XLAL_EFAILED );
    }

    /* Try extracting all 100 elements, then adding them back */
    for ( int i = 0; i < 100; ++i ) {
      elem x = { .key = 100*n + i };
      elem *y;
      XLAL_CHECK_MAIN( XLALHashTblExtract( ht, &x, ( void ** ) &y ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( y != NULL, XLAL_EFAILED );
      XLAL_CHECK_MAIN( y->value == 3*y->key - n, XLAL_EFAILED );
      const elem *z;
      XLAL_CHECK_MAIN( XLALHashTblFind( ht, &x, ( const void ** ) &z ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( z == NULL, XLAL_EFAILED );
      XLAL_CHECK_MAIN( XLALHashTblAdd( ht, y ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( XLALHashTblFind( ht, &x, ( const void ** ) &z ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( z != NULL, XLAL_EFAILED );
      XLAL_CHECK_MAIN( z->value == 3*z->key - n, XLAL_EFAILED );
    }

  }
  XLAL_CHECK_MAIN( XLALHashTblSize( ht ) == 400, XLAL_EFAILED );

  /* Try removing some elements */
  for ( int i = 0; i < 250; ++i ) {
    elem x = { .key = i };
    XLAL_CHECK_MAIN( XLALHashTblRemove( ht, &x ) == XLAL_SUCCESS, XLAL_EFAILED );
    XLAL_CHECK( XLALHashTblSize( ht ) == 400 - i - 1, XLAL_EFAILED );
  }

  /* Try finding the rest of the elements */
  for ( int i = 250; i < 400; ++i ) {
    elem x = { .key = i };
    const elem *y;
    XLAL_CHECK_MAIN( XLALHashTblFind( ht, &x, ( const void ** ) &y ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK_MAIN( y != NULL, XLAL_EFAILED );
    XLAL_CHECK_MAIN( y->value == 3*y->key - ( i / 100 ), XLAL_EFAILED );
  }

  /* Cleanup */
  XLALHashTblDestroy( ht );

  /* Check for memory leaks */
  LALCheckMemoryLeaks();

  return EXIT_SUCCESS;

}

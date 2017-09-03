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

#include <limits.h>

#include <lal/LALBitset.h>
#include <lal/LALHashTbl.h>

#define BITS_PER_ELEM (sizeof(UINT8) * CHAR_BIT)

struct tagLALBitset {
  LALHashTbl *ht;               /* Hash table which stores bits in UINT8s */
};

typedef struct {
  UINT8 key;
  UINT8 bits;
} elem;

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

LALBitset *XLALBitsetCreate(
  void
  )
{

  /* Allocate memory for bitset struct */
  LALBitset *bs = XLALCalloc( 1, sizeof( *bs ) );
  XLAL_CHECK_NULL( bs != NULL, XLAL_ENOMEM );

  /* Create hash table */
  bs->ht = XLALHashTblCreate( XLALFree, hash_elem, cmp_elem );
  XLAL_CHECK_NULL( bs->ht != NULL, XLAL_EFUNC );

  return bs;

}

void XLALBitsetDestroy(
  LALBitset *bs
  )
{
  if ( bs != NULL ) {
    XLALHashTblDestroy( bs->ht );
    XLALFree( bs );
  }
}

int XLALBitsetClear(
  LALBitset *bs
  )
{

  /* Check input */
  XLAL_CHECK( bs != NULL, XLAL_EFAULT );

  /* Clear hash table */
  XLAL_CHECK( XLALHashTblClear( bs->ht ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;

}

int XLALBitsetSet(
  LALBitset *bs,
  const UINT8 idx,
  const BOOLEAN is_set
  )
{

  /* Check input */
  XLAL_CHECK( bs != NULL, XLAL_EFAULT );

  /* Compute key and bit index */
  const UINT8 key = idx / BITS_PER_ELEM;
  const UINT8 bitidx = idx % BITS_PER_ELEM;

  /* Extract element corresponding to key, or else create new element */
  const elem x = { .key = key };
  elem *y = NULL;
  XLAL_CHECK( XLALHashTblExtract( bs->ht, &x, ( void ** ) &y ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( y == NULL ) {
    y = XLALCalloc( 1, sizeof( *y ) );
    XLAL_CHECK( y != NULL, XLAL_ENOMEM );
    y->key = key;
  }

  /* Set/unset bit in element */
  if ( is_set ) {
    y->bits |=  ( ( (UINT8) 1 ) << bitidx );
  } else {
    y->bits &= ~( ( (UINT8) 1 ) << bitidx );
  }

  /* Add element back into hash */
  XLAL_CHECK( XLALHashTblAdd( bs->ht, y ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;

}

int XLALBitsetGet(
  const LALBitset *bs,
  const UINT8 idx,
  BOOLEAN *is_set
  )
{

  /* Check input */
  XLAL_CHECK( bs != NULL, XLAL_EFAULT );
  XLAL_CHECK( is_set != NULL, XLAL_EFAULT );

  /* Compute key and bit index */
  const UINT8 key = idx / BITS_PER_ELEM;
  const UINT8 bitidx = idx % BITS_PER_ELEM;

  /* Find element corresponding to key */
  const elem x = { .key = key };
  const elem *y = NULL;
  XLAL_CHECK( XLALHashTblFind( bs->ht, &x, ( const void ** ) &y ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( y == NULL ) {
    *is_set = 0;
  } else {

    /* Check if bit is set in element */
    *is_set = ( y->bits & ( ( (UINT8) 1 ) << bitidx ) ) ? 1 : 0;

  }

  return XLAL_SUCCESS;

}

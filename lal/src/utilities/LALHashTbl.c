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

#include <lal/LALHashTbl.h>

/* Special hash table element value to indicate elements that have been deleted */
static const void *hash_del = 0;
#define DEL   ((void*) &hash_del)

/* Evaluates to the hash value of x, restricted to the length of the hash table */
#define HASHIDX(ht, x)   ((int)((ht)->hash((ht)->hash_param, (x)) % (ht)->data_len))

/* Increment the next hash index, restricted to the length of the hash table */
#define INCRIDX(ht, i)   do { if (++(i) == (ht)->data_len) { (i) = 0; } } while(0)

/* Evaluates true if the elements x and y are equal, according to the hash table comparison function */
#define EQUAL(ht, x, y)   ((ht)->cmp((ht)->cmp_param, (x), (y)) == 0)

struct tagLALHashTbl {
  void **data;                  /* Hash table with open addressing and linear probing */
  int data_len;                 /* Size of the memory block 'data', in number of elements */
  int n;                        /* Number of valid elements in the hash */
  int q;                        /* Number of non-NULL elements in the hash */
  LALHashTblDtorFcn dtor;       /* Function to free memory of elements of hash, if required */
  LALHashTblHashParamFcn hash;  /* Parameterised hash function for hash table elements */
  void *hash_param;             /* Parameter to pass to hash function */
  LALHashTblCmpParamFcn cmp;    /* Parameterised hash table element comparison function */
  void *cmp_param;              /* Parameter to pass to comparison function */
};

/* Call a non-parameterised hash function, which is passed in 'param' */
static UINT8 hashtbl_no_param_hash( void *param, const void *x )
{
  LALHashTblHashFcn hash = ( LALHashTblHashFcn ) param;
  return hash( x );
}

/* Call a non-parameterised compare function, which is passed in 'param' */
static int hashtbl_no_param_cmp( void *param, const void *x, const void *y )
{
  LALHashTblCmpFcn cmp = ( LALHashTblCmpFcn ) param;
  return cmp( x, y );
}

/* Resize and rebuild the hash table */
static int hashtbl_resize( LALHashTbl *ht )
{
  void **old_data = ht->data;
  int old_data_len = ht->data_len;
  ht->data_len = 2;
  while ( ht->data_len < 3*ht->n ) {
    ht->data_len *= 2;
  }
  ht->data = XLALCalloc( ht->data_len, sizeof( ht->data[0] ) );
  XLAL_CHECK( ht->data != NULL, XLAL_ENOMEM );
  ht->q = ht->n;
  for ( int k = 0; k < old_data_len; ++k ) {
    if ( old_data[k] != NULL && old_data[k] != DEL ) {
      int i = HASHIDX( ht, old_data[k] );
      while ( ht->data[i] != NULL ) {
        INCRIDX( ht, i );
      }
      ht->data[i] = old_data[k];
    }
  }
  XLALFree( old_data );
  return XLAL_SUCCESS;
}

LALHashTbl *XLALHashTblCreate(
  LALHashTblDtorFcn dtor,
  LALHashTblHashFcn hash,
  LALHashTblCmpFcn cmp
  )
{

  /* Create a hash table using hashtbl_no_param_hash/cmp as the hash/comparison functions */
  LALHashTbl *ht = XLALHashTblCreate2( dtor, hashtbl_no_param_hash, hash, hashtbl_no_param_cmp, cmp );
  XLAL_CHECK_NULL( ht != NULL, XLAL_EFUNC );

  return ht;

}

LALHashTbl *XLALHashTblCreate2(
  LALHashTblDtorFcn dtor,
  LALHashTblHashParamFcn hash,
  void *hash_param,
  LALHashTblCmpParamFcn cmp,
  void *cmp_param
  )
{

  /* Check input */
  XLAL_CHECK_NULL( hash != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( cmp != NULL, XLAL_EFAULT );

  /* Allocate memory for hash table struct */
  LALHashTbl *ht = XLALCalloc( 1, sizeof( *ht ) );
  XLAL_CHECK_NULL( ht != NULL, XLAL_ENOMEM );

  /* Set hash table struct parameters */
  ht->dtor = dtor;
  ht->hash = hash;
  ht->hash_param = hash_param;
  ht->cmp = cmp;
  ht->cmp_param = cmp_param;

  return ht;

}

void XLALHashTblDestroy(
  LALHashTbl *ht
  )
{
  if ( ht != NULL ) {
    if ( ht->data != NULL ) {
      if ( ht->dtor != NULL ) {
        for ( int i = 0; i < ht->data_len; ++i ) {
          if ( ht->data[i] != NULL && ht->data[i] != DEL ) {
            ht->dtor( ht->data[i] );
          }
        }
      }
      XLALFree( ht->data );
    }
    XLALFree( ht );
  }
}

int XLALHashTblSize(
  const LALHashTbl *ht
  )
{
  XLAL_CHECK( ht != NULL, XLAL_EFAULT );
  return ht->n;
}

int XLALHashTblFind(
  const LALHashTbl *ht,
  const void *x,
  const void **y
  )
{

  /* Check input */
  XLAL_CHECK( ht != NULL, XLAL_EFAULT );
  XLAL_CHECK( x != NULL && x != DEL, XLAL_EINVAL );
  XLAL_CHECK( y != NULL, XLAL_EFAULT );

  /* Try to find element matching 'x' in hash table, if found return in 'y' */
  if ( ht->data_len > 0 ) {
    int i = HASHIDX( ht, x );
    while ( ht->data[i] != NULL ) {
      *y = ht->data[i];
      if ( *y != DEL && EQUAL( ht, x, *y ) ) {
        return XLAL_SUCCESS;
      }
      INCRIDX( ht, i );
    }
  }

  /* No element matches 'x' */
  *y = NULL;
  return XLAL_SUCCESS;

}

int XLALHashTblAdd(
  LALHashTbl *ht,
  void *x
  )
{

  /* Check input */
  XLAL_CHECK( ht != NULL, XLAL_EFAULT );
  XLAL_CHECK( x != NULL && x != DEL, XLAL_EINVAL );

  /* Check that no element matching 'x' exists in the hash table */
  {
    const void *y;
    XLAL_CHECK( XLALHashTblFind( ht, x, &y ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( y == NULL, XLAL_EFAILED, "Hash table already contains given element" );
  }

  /* Resize hash table to preserve maximum 50% occupancy */
  if ( 2*( ht->q + 1 ) > ht->data_len ) {
    XLAL_CHECK( hashtbl_resize( ht ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  /* Add 'x' to the hash table */
  int i = HASHIDX( ht, x );
  while ( ht->data[i] != NULL && ht->data[i] != DEL ) {
    INCRIDX( ht, i );
  }
  if ( ht->data[i] == NULL ) {
    ++ht->q;
  }
  ++ht->n;
  ht->data[i] = x;

  return XLAL_SUCCESS;

}

int XLALHashTblExtract(
  LALHashTbl *ht,
  const void *x,
  void **y
  )
{

  /* Check input */
  XLAL_CHECK( ht != NULL, XLAL_EFAULT );
  XLAL_CHECK( x != NULL && x != DEL, XLAL_EINVAL );
  XLAL_CHECK( y != NULL, XLAL_EFAULT );

  /* Try to find element matching 'x' in hash table, if found remove it from table and return in 'y' */
  if ( ht->data_len > 0 ) {
    int i = HASHIDX( ht, x );
    while ( ht->data[i] != NULL ) {
      *y = ht->data[i];
      if ( *y != DEL && EQUAL( ht, x, *y ) ) {
        ht->data[i] = DEL;
        --ht->n;
        if ( 8*ht->n < ht->data_len ) { /* Resize hash table to preserve minimum 50% occupancy */
          XLAL_CHECK( hashtbl_resize( ht ) == XLAL_SUCCESS, XLAL_EFUNC );
        }
        return XLAL_SUCCESS;
      }
      INCRIDX( ht, i );
    }
  }

  /* No element matches 'x' */
  *y = NULL;
  return XLAL_SUCCESS;

}

int XLALHashTblRemove(
  LALHashTbl *ht,
  const void *x
  )
{

  /* Check input */
  XLAL_CHECK( ht != NULL, XLAL_EFAULT );
  XLAL_CHECK( x != NULL && x != DEL, XLAL_EINVAL );

  /* Remove element matching 'x' from hash table, if it exists */
  void *y;
  XLAL_CHECK( XLALHashTblExtract( ht, x, &y ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* Free memory associated with element, if required */
  if ( ht->dtor != NULL ) {
    ht->dtor( y );
  }

  return XLAL_SUCCESS;

}

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
static const void* hash_del = 0;
#define DEL   ((void*) &hash_del)

/* Evaluates to the hash value of x, restricted to the length of the hash table */
#define HASHIDX(ht, x)   ((int)((ht)->hash((ht)->hash_param, (x)) % (ht)->data_len))

/* Increment the next hash index, restricted to the length of the hash table */
#define INCRIDX(ht, i)   do { if (++(i) == (ht)->data_len) { (i) = 0; } } while(0)

/* Evaluates true if the elements x and y are equal, according to the hash table comparison function */
#define EQUAL(ht, x, y)   ((ht)->cmp((ht)->cmp_param, (x), (y)) == 0)

struct tagLALHashTbl {
  void **data;			/* Hash table with open addressing and linear probing */
  int data_len;			/* Size of the memory block 'data', in number of elements */
  int n;			/* Number of valid elements in the hash */
  int q;			/* Number of non-NULL elements in the hash */
  LALHashTblDtorFcn dtor;	/* Function to free memory of elements of hash, if required */
  LALHashTblHashParamFcn hash;	/* Parameterised hash function for hash table elements */
  void *hash_param;		/* Parameter to pass to hash function */
  LALHashTblCmpParamFcn cmp;	/* Parameterised hash table element comparison function */
  void *cmp_param;		/* Parameter to pass to comparison function */
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

UINT8 XLALComputeHash(
  const void *data,
  const size_t len
  )
{
  if ( data == NULL || len == 0 ) {
    return 0;
  }

  /* Compute the hash value using Pearson hashing, following https://en.wikipedia.org/wiki/Pearson_hashing */
  const unsigned char T[256] = {
    171,38,41,159,93,52,71,169,246,73,191,232,196,33,13,34,
    31,143,228,59,96,63,91,23,250,36,221,126,214,57,90,94,
    157,58,88,213,29,4,83,179,176,121,77,0,226,25,72,15,
    102,76,153,10,137,134,254,20,128,65,198,89,229,47,160,12,
    42,189,203,141,100,40,53,78,182,35,130,68,197,212,55,111,
    18,194,131,252,22,125,147,124,39,11,21,174,249,97,209,144,
    218,27,82,5,220,67,129,150,92,193,48,45,139,216,255,110,
    140,156,192,87,215,69,185,104,175,181,109,75,231,138,180,105,
    6,80,54,190,135,227,50,164,146,8,148,115,123,32,206,7,
    238,74,204,223,177,248,149,188,230,239,170,51,61,30,24,98,
    28,1,37,200,46,184,244,113,116,154,84,86,199,112,133,107,
    145,207,241,166,162,101,172,208,205,236,211,106,19,132,225,56,
    108,99,186,9,251,183,210,114,245,242,187,62,165,60,127,17,
    161,64,3,234,237,16,168,152,120,95,167,155,49,219,81,202,
    122,26,136,44,70,224,79,117,201,253,222,103,233,2,243,178,
    235,240,66,119,173,195,163,85,247,151,158,142,43,14,217,118,
  };
  const char *x = ( const char * ) data;
  union {
    UINT8 val;
    unsigned char hh[8];
  } out;
  out.val = 0;
  for ( size_t j = 0; j < 8; ++j ) {
    unsigned char h = T[( x[0] + j ) % 256];
    for ( size_t i = 1; i < len; ++i ) {
      h = T[h ^ x[i]];
    }
    out.hh[j] = h;
  }

  return out.val;

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

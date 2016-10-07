/*
*  Copyright (C) 2007 Jolien Creighton
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
 * File Name: HeapSort.c
 *
 * Author: Creighton, T. D.
 *
 *
 *-----------------------------------------------------------------------*/

/* ---------- see Sort.h for doxygen documentation ---------- */

#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/Sort.h>

/* helpful macros for generic routines */
/* copy element j of array y to element i of array x; elements have size s */
#define COPY(x,i,y,j,s) (memcpy((char*)(x)+(i)*(s),(char*)(y)+(j)*(s),(s)))
/* compare element j of array y with element i of array x; elements have size s
 * use compare function c with params p */
#define CMP(x,i,y,j,s,p,c) ((c)((p),(char*)(x)+(i)*(s),(char*)(y)+(j)*(s)))


int XLALHeapSort( void *base, UINT4 nobj, UINT4 size, void *params,
    int (*compar)(void *, const void *, const void *) )
{
  INT4 i, j, k, n = nobj;
  void *temp;

  if ( (INT4)size <= 0 )
    XLAL_ERROR( XLAL_EINVAL );
  if ( ! base || ! compar )
    XLAL_ERROR( XLAL_EFAULT );

  /* 0 or 1 objects are already sorted. */
  if (n<2)
    return 0;

  temp = LALMalloc( size );
  if ( ! temp )
    XLAL_ERROR( XLAL_ENOMEM );

  /* Here is the heapsort algorithm. */
  j=n-1;
  n >>= 1;
  while(1){
    if(n>0)
      COPY(temp,0,base,--n,size);
    else{
      COPY(temp,0,base,j,size);
      COPY(base,j,base,0,size);
      if(--j==0){
        COPY(base,0,temp,0,size);
        break;
      }
    }
    i=n;
    k=(n << 1)+1;
    while(k<=j){
      if(k<j && CMP(base,k,base,k+1,size,params,compar) < 0)
        k++;
      if(CMP(temp,0,base,k,size,params,compar) < 0){
        COPY(base,i,base,k,size);
        i=k;
        k<<=1;
        k++;
      }else
        k=j+1;
    }
    COPY(base,i,temp,0,size);
  }

  LALFree( temp );
  return 0;
}


int XLALHeapIndex( INT4 *indx, void *base, UINT4 nobj, UINT4 size, void *params,
    int (*compar)(void *, const void *, const void *) )
{
  INT4 i, j, k, n = nobj;
  INT4 itemp;

  if ( (INT4)size <= 0 )
    XLAL_ERROR( XLAL_EINVAL );
  if ( ! indx || ! base || ! compar )
    XLAL_ERROR( XLAL_EFAULT );

  /* Initialize the indx vector */
  for (i=0;i<n;++i)
    indx[i]=i;

  /* 0 or 1 objects are already sorted. */
  if (n<2)
    return 0;

  /* Here is the heapsort algorithm. */
  j=n-1;
  n >>= 1;

  while(1){
    if(n>0)
      itemp=indx[--n];
    else{
      itemp=indx[j];
      indx[j]=indx[0];
      if(--j==0){
        indx[0]=itemp;
        break;
      }
    }
    i=n;
    k=(n << 1)+1;
    while(k<=j){
      if(k<j && CMP(base,indx[k],base,indx[k+1],size,params,compar)<0)
        k++;
      if(CMP(base,itemp,base,indx[k],size,params,compar)<0){
        indx[i]=indx[k];
        i=k;
        k<<=1;
        k++;
      }else
        k=j+1;
    }
    indx[i]=itemp;
  }

  return 0;
}


int XLALHeapRank( INT4 *rank, void *base, UINT4 nobj, UINT4 size, void *params,
    int (*compar)(void *, const void *, const void *) )
{
  INT4 i, n = nobj;
  INT4 *indx;

  if ( (INT4)size <= 0 )
    XLAL_ERROR( XLAL_EINVAL );
  if ( ! rank || ! base || ! compar )
    XLAL_ERROR( XLAL_EFAULT );

  indx = LALMalloc( nobj * sizeof( *rank ) );
  if ( ! indx )
    XLAL_ERROR( XLAL_ENOMEM );

  if ( XLALHeapIndex( indx, base, nobj, size, params, compar ) < 0 )
    XLAL_ERROR( XLAL_EFUNC );

  for(i=0;i<n;++i)
    rank[indx[i]]=i;

  LALFree( indx );
  return 0;
}

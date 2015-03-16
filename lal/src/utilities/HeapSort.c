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

#undef COPY
#undef CMP


void LALSHeapSort(LALStatus      *stat,
	       REAL4Vector *vector)
{
  INT4 i;
  INT4 j;
  INT4 k;
  INT4 n;
  REAL4 temp;
  REAL4 *data;

  INITSTATUS(stat);

  /* Make sure all pointers are valid. */
  ASSERT(vector,stat,SORTH_ENUL,SORTH_MSGENUL);
  ASSERT(vector->data,stat,SORTH_ENUL,SORTH_MSGENUL);

  n=vector->length;
  data=vector->data;

  /* A vector of length 0 or 1 is already sorted. */
  if(n<2)
  {
    RETURN(stat);
  }

  /* Here is the heapsort algorithm. */
  j=n-1;
  n >>= 1;

  while(1){
    if(n>0)
      temp=data[--n];
    else{
      temp=data[j];
      data[j]=data[0];
      if(--j==0){
	data[0]=temp;
	break;
      }
    }
    i=n;
    k=(n << 1)+1;
    while(k<=j){
      if(k<j && data[k]<data[k+1])
	k++;
      if(temp<data[k]){
	data[i]=data[k];
	i=k;
	k<<=1;
        k++;
      }else
	k=j+1;
    }
    data[i]=temp;
  }
}



void LALSHeapIndex(LALStatus      *stat,
		INT4Vector  *idx,
		REAL4Vector *vector)
{
  INT4 i;
  INT4 j;
  INT4 k;
  INT4 n;
  INT4 itemp;
  INT4 *indx;
  REAL4 temp;
  REAL4 *data;

  INITSTATUS(stat);

  /* Make sure all pointers are valid, and the idx vector is of the
     same length as the data vector. */
  ASSERT(vector,stat,SORTH_ENUL,SORTH_MSGENUL);
  ASSERT(idx,stat,SORTH_ENUL,SORTH_MSGENUL);
  ASSERT(vector->data,stat,SORTH_ENUL,SORTH_MSGENUL);
  ASSERT(idx->data,stat,SORTH_ENUL,SORTH_MSGENUL);
  ASSERT(vector->length==idx->length,stat,SORTH_ELEN,SORTH_MSGELEN);

  n=vector->length;
  data=vector->data;

  /* Initialize the idx vector. */
  for(i=0,indx=idx->data;i<n;i++,indx++)
    *indx=i;
  indx=idx->data;

  /* A vector of length 0 or 1 is already sorted. */
  if(n<2)
  {
    RETURN(stat);
  }

  /* Here is the heapsort algorithm. */
  j=n-1;
  n >>= 1;

  while(1){
    if(n>0){
      itemp=indx[--n];
      temp=data[itemp];
    }
    else{
      itemp=indx[j];
      temp=data[itemp];
      indx[j]=indx[0];
      if(--j==0){
	indx[0]=itemp;
	break;
      }
    }
    i=n;
    k=(n << 1)+1;
    while(k<=j){
      if(k<j && data[indx[k]]<data[indx[k+1]])
	k++;
      if(temp<data[indx[k]]){
	indx[i]=indx[k];
	i=k;
	k<<=1;
        k++;
      }else
	k=j+1;
    }
    indx[i]=itemp;
  }
}



void LALSHeapRank(LALStatus      *stat,
	       INT4Vector  *rank,
	       REAL4Vector *vector)
{
  INT4       i;
  INT4       *indx;
  INT4       *rnk;
  INT4Vector *idx=NULL;

  INITSTATUS(stat);
  ATTATCHSTATUSPTR(stat);

  /* Make sure all pointers are valid, and the rank vector is of the
     same length as the data vector. */
  ASSERT(vector,stat,SORTH_ENUL,SORTH_MSGENUL);
  ASSERT(rank,stat,SORTH_ENUL,SORTH_MSGENUL);
  ASSERT(vector->data,stat,SORTH_ENUL,SORTH_MSGENUL);
  ASSERT(rank->data,stat,SORTH_ENUL,SORTH_MSGENUL);
  ASSERT(vector->length==rank->length,stat,SORTH_ELEN,SORTH_MSGELEN);

  /* Make the temporary idx vector. */
  TRY(LALI4CreateVector(stat->statusPtr,&idx,vector->length),stat);
  TRY(LALSHeapIndex(stat->statusPtr,idx,vector),stat);

  /* Invert to get the rank vector. */
  indx=idx->data;
  rnk=rank->data;
  for(i=0;i<(int)vector->length;i++)
  {
    rnk[indx[i]]=i;
  }

  /* Free memory and exit. */
  TRY(LALI4DestroyVector(stat->statusPtr,&idx),stat);
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}



void LALDHeapSort(LALStatus      *stat,
	       REAL8Vector *vector)
{
  INT4 i;
  INT4 j;
  INT4 k;
  INT4 n;
  REAL8 temp;
  REAL8 *data;

  INITSTATUS(stat);

  /* Make sure all pointers are valid. */
  ASSERT(vector,stat,SORTH_ENUL,SORTH_MSGENUL);
  ASSERT(vector->data,stat,SORTH_ENUL,SORTH_MSGENUL);

  n=vector->length;
  data=vector->data;

  /* A vector of length 0 or 1 is already sorted. */
  if(n<2)
  {
    RETURN(stat);
  }

  /* Here is the heapsort algorithm. */
  j=n-1;
  n >>= 1;

  while(1){
    if(n>0)
      temp=data[--n];
    else{
      temp=data[j];
      data[j]=data[0];
      if(--j==0){
	data[0]=temp;
	break;
      }
    }
    i=n;
    k=(n << 1)+1;
    while(k<=j){
      if(k<j && data[k]<data[k+1])
	k++;
      if(temp<data[k]){
	data[i]=data[k];
	i=k;
	k<<=1;
        k++;
      }else
	k=j+1;
    }
    data[i]=temp;
  }
}



void LALDHeapIndex(LALStatus      *stat,
		INT4Vector  *idx,
		REAL8Vector *vector)
{
  INT4 i;
  INT4 j;
  INT4 k;
  INT4 n;
  INT4 itemp;
  INT4 *indx;
  REAL8 temp;
  REAL8 *data;

  INITSTATUS(stat);

  /* Make sure all pointers are valid, and the idx vector is of the
     same length as the data vector. */
  ASSERT(vector,stat,SORTH_ENUL,SORTH_MSGENUL);
  ASSERT(idx,stat,SORTH_ENUL,SORTH_MSGENUL);
  ASSERT(vector->data,stat,SORTH_ENUL,SORTH_MSGENUL);
  ASSERT(idx->data,stat,SORTH_ENUL,SORTH_MSGENUL);
  ASSERT(vector->length==idx->length,stat,SORTH_ELEN,SORTH_MSGELEN);

  n=vector->length;
  data=vector->data;

  /* Initialize the idx vector. */
  for(i=0,indx=idx->data;i<n;i++,indx++)
  {
    *indx=i;
  }
  indx=idx->data;

  /* A vector of length 0 or 1 is already sorted. */
  if(n<2)
  {
    RETURN(stat);
  }

  /* Here is the heapsort algorithm. */
  j=n-1;
  n >>= 1;

  while(1){
    if(n>0){
      itemp=indx[--n];
      temp=data[itemp];
    }
    else{
      itemp=indx[j];
      temp=data[itemp];
      indx[j]=indx[0];
      if(--j==0){
	indx[0]=itemp;
	break;
      }
    }
    i=n;
    k=(n << 1)+1;
    while(k<=j){
      if(k<j && data[indx[k]]<data[indx[k+1]])
	k++;
      if(temp<data[indx[k]]){
	indx[i]=indx[k];
	i=k;
	k<<=1;
        k++;
      }else
	k=j+1;
    }
    indx[i]=itemp;
  }
}



void LALDHeapRank(LALStatus      *stat,
	       INT4Vector  *rank,
	       REAL8Vector *vector)
{
  INT4       i;
  INT4       *indx;
  INT4       *rnk;
  INT4Vector *idx=NULL;

  INITSTATUS(stat);
  ATTATCHSTATUSPTR(stat);

  /* Make sure all pointers are valid, and the rank vector is of the
     same length as the data vector. */
  ASSERT(vector,stat,SORTH_ENUL,SORTH_MSGENUL);
  ASSERT(rank,stat,SORTH_ENUL,SORTH_MSGENUL);
  ASSERT(vector->data,stat,SORTH_ENUL,SORTH_MSGENUL);
  ASSERT(rank->data,stat,SORTH_ENUL,SORTH_MSGENUL);
  ASSERT(vector->length==rank->length,stat,SORTH_ELEN,SORTH_MSGELEN);

  /* Make the temporary idx vector. */
  TRY(LALI4CreateVector(stat->statusPtr,&idx,vector->length),stat);
  TRY(LALDHeapIndex(stat->statusPtr,idx,vector),stat);

  /* Invert to get the rank vector. */
  indx=idx->data;
  rnk=rank->data;
  for(i=0;i<(int)vector->length;i++)
    rnk[indx[i]]=i;

  /* Free memory and exit. */
  TRY(LALI4DestroyVector(stat->statusPtr,&idx),stat);
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}

/*----------------------------------------------------------------------- 
 * 
 * File Name: HeapSort.c
 * 
 * Author: Creighton, T. D.
 * 
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------*/

/*
 
<lalVerbatim file="HeapSortCV">
$Id$
</lalVerbatim>

<lalLaTeX>

\subsection{Module \texttt{HeapSort.c}}

Sorts, indexes, or ranks vector elements using the heap sort
algorithm.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{HeapSortCP}
\idx{XLALHeapSort()}
\idx{XLALHeapIndex()}
\idx{XLALHeapRank()}
\idx{LALSHeapSort()}
\idx{LALSHeapIndex()}
\idx{LALSHeapRank()}
\idx{LALDHeapSort()}
\idx{LALDHeapIndex()}
\idx{LALDHeapRank()}

\subsubsection*{Description}

These routines sort a vector \verb@*data@ (of type \verb@REAL4Vector@
or \verb@REAL8Vector@) into ascending order using the in-place
heapsort algorithm, or construct an index vector \verb@*index@ that
indexes \verb@*data@ in increasing order (leaving \verb@*data@
unchanged), or construct a rank vector \verb@*rank@ that gives the
rank order of the corresponding \verb@*data@ element.

The relationship between sorting, indexing, and ranking can be a bit
confusing.  One way of looking at it is that the original array is
ordered by index, while the sorted array is ordered by rank.  The
index array gives the index as a function of rank; i.e.\ if you're
looking for a given rank (say the 0th, or smallest element), the index
array tells you where to look it up in the unsorted array:
\begin{verbatim}
unsorted_array[index[i]] = sorted_array[i]
\end{verbatim}
The rank array gives the rank as a function of index; i.e.\ it tells
you where a given element in the unsorted array will appear in the
sorted array:
\begin{verbatim}
unsorted_array[j] = sorted_array[rank[j]]
\end{verbatim}
Clearly these imply the following relationships, which can be used to
construct the index array from the rank array or vice-versa:
\begin{verbatim}
index[rank[j]] = j
rank[index[i]] = i
\end{verbatim}

The XLAL versions of these routines, \verb@XLALHeapSort@, \verb@XLALHeapIndex@,
and \verb@XLALHeapRank@, perform the same operations but on arrays of
\verb@nobj@ generic objects of size \verb@size@ pointed to by \verb@base@ and
using the comparison function \verb@compar@.  The function \verb@compar@ has
the prototype
\begin{verbatim}
int compar( void *p, const void *x, const void *y )
\end{verbatim}
and returns $-1$ if ${\mathtt{x}}<{\mathtt{y}}$,
$0$ if ${\mathtt{x}}={\mathtt{y}}$,
and $+1$ if ${\mathtt{x}}>{\mathtt{y}}$.  Here \verb@p@ (which may be NULL)
is a pointer to additional data that may be used in the comparison function.
This pointer is passed to the comparison function unchanged from the argument
\verb@params@ of \verb@XLALHeapSort@, \verb@XLALHeapIndex@, and
\verb@XLALHeapRank@.


\subsubsection*{Algorithm}

These routines use the standard heap sort algorithm described in
Sec.~8.3 of Ref.~\cite{ptvf:1992}.

The \verb@LALSHeapSort()@ and \verb@LALDHeapSort()@ routines are entirely
in-place, with no auxiliary storage vector.  The \verb@LALSHeapIndex()@
and \verb@LALDHeapIndex()@ routines are also technically in-place, but
they require two input vectors (the data vector and the index vector),
and leave the data vector unchanged.  The \verb@LALSHeapRank()@ and
\verb@LALDHeapRank()@ routines require two input vectors (the data and
rank vectors), and also allocate a temporary index vector internally;
these routines are therefore the most memory-intensive.  All of these
algorithms are $N\log_2(N)$ algorithms, regardless of the ordering of
the initial dataset.

Note: if you can use \verb@qsort@, you should.

\subsubsection*{Uses}
\begin{verbatim}
LALI4CreateVector()
LALI4DestroyVector()
\end{verbatim}

\subsubsection*{Notes}
\vfill{\footnotesize\input{HeapSortCV}}

</lalLaTeX> */

#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/Sort.h>

NRCSID(HEAPSORTC,"$Id$");

/* helpful macros for generic routines */
/* copy element j of array y to element i of array x; elements have size s */
#define COPY(x,i,y,j,s) (memcpy((char*)(x)+(i)*(s),(char*)(y)+(j)*(s),(s)))
/* compare element j of array y with element i of array x; elements have size s
 * use compare function c with params p */
#define CMP(x,i,y,j,s,p,c) ((c)((p),(char*)(x)+(i)*(s),(char*)(y)+(j)*(s)))

/* <lalVerbatim file="HeapSortCP"> */
int XLALHeapSort( void *base, UINT4 nobj, UINT4 size, void *params,
    int (*compar)(void *, const void *, const void *) )
{ /* </lalVerbatim> */
  static const char *func = "XLALHeapSort";
  INT4 i, j, k, n = nobj;
  void *temp;

  if ( (INT4)size <= 0 )
    XLAL_ERROR( func, XLAL_EINVAL );
  if ( ! base || ! compar )
    XLAL_ERROR( func, XLAL_EFAULT );

  /* 0 or 1 objects are already sorted. */
  if (n<2)
    return 0;

  temp = LALMalloc( size );
  if ( ! temp )
    XLAL_ERROR( func, XLAL_ENOMEM );

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

/* <lalVerbatim file="HeapSortCP"> */
int XLALHeapIndex( INT4 *indx, void *base, UINT4 nobj, UINT4 size, void *params,
    int (*compar)(void *, const void *, const void *) )
{ /* </lalVerbatim> */
  static const char *func = "XLALHeapIndex";
  INT4 i, j, k, n = nobj;
  INT4 itemp;

  if ( (INT4)size <= 0 )
    XLAL_ERROR( func, XLAL_EINVAL );
  if ( ! indx || ! base || ! compar )
    XLAL_ERROR( func, XLAL_EFAULT );

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

/* <lalVerbatim file="HeapSortCP"> */
int XLALHeapRank( INT4 *rank, void *base, UINT4 nobj, UINT4 size, void *params,
    int (*compar)(void *, const void *, const void *) )
{ /* </lalVerbatim> */
  static const char *func = "XLALHeapRank";
  INT4 i, n = nobj;
  INT4 *indx;

  if ( (INT4)size <= 0 )
    XLAL_ERROR( func, XLAL_EINVAL );
  if ( ! rank || ! base || ! compar )
    XLAL_ERROR( func, XLAL_EFAULT );

  indx = LALMalloc( nobj * sizeof( *rank ) );
  if ( ! indx )
    XLAL_ERROR( func, XLAL_ENOMEM );

  if ( XLALHeapIndex( indx, base, nobj, size, params, compar ) < 0 )
    XLAL_ERROR( func, XLAL_EFUNC );

  for(i=0;i<n;++i)
    rank[indx[i]]=i;

  LALFree( indx );
  return 0;
}

#undef COPY
#undef CMP

/* <lalVerbatim file="HeapSortCP"> */
void LALSHeapSort(LALStatus      *stat,
	       REAL4Vector *vector)
{ /* </lalVerbatim> */
  INT4 i;
  INT4 j;
  INT4 k;
  INT4 n;
  REAL4 temp;
  REAL4 *data;

  INITSTATUS(stat,"LALSHeapSort",HEAPSORTC);

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


/* <lalVerbatim file="HeapSortCP"> */
void LALSHeapIndex(LALStatus      *stat,
		INT4Vector  *idx,
		REAL4Vector *vector)
{ /* </lalVerbatim> */
  INT4 i;
  INT4 j;
  INT4 k;
  INT4 n;
  INT4 itemp;
  INT4 *indx;
  REAL4 temp;
  REAL4 *data;

  INITSTATUS(stat,"LALSHeapIndex",HEAPSORTC);

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


/* <lalVerbatim file="HeapSortCP"> */
void LALSHeapRank(LALStatus      *stat,
	       INT4Vector  *rank,
	       REAL4Vector *vector)
{ /* </lalVerbatim> */
  INT4       i;
  INT4       *indx;
  INT4       *rnk;
  INT4Vector *idx=NULL;

  INITSTATUS(stat,"LALSHeapRank",HEAPSORTC);
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


/* <lalVerbatim file="HeapSortCP"> */
void LALDHeapSort(LALStatus      *stat,
	       REAL8Vector *vector)
{ /* </lalVerbatim> */
  INT4 i;
  INT4 j;
  INT4 k;
  INT4 n;
  REAL8 temp;
  REAL8 *data;

  INITSTATUS(stat,"LALDHeapSort",HEAPSORTC);

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


/* <lalVerbatim file="HeapSortCP"> */
void LALDHeapIndex(LALStatus      *stat,
		INT4Vector  *idx,
		REAL8Vector *vector)
{ /* </lalVerbatim> */
  INT4 i;
  INT4 j;
  INT4 k;
  INT4 n;
  INT4 itemp;
  INT4 *indx;
  REAL8 temp;
  REAL8 *data;

  INITSTATUS(stat,"LALDHeapIndex",HEAPSORTC);

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


/* <lalVerbatim file="HeapSortCP"> */
void LALDHeapRank(LALStatus      *stat,
	       INT4Vector  *rank,
	       REAL8Vector *vector)
{ /* </lalVerbatim> */
  INT4       i;
  INT4       *indx;
  INT4       *rnk;
  INT4Vector *idx=NULL;

  INITSTATUS(stat,"LALDHeapRank",HEAPSORTC);
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

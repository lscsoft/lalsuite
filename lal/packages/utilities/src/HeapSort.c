/*----------------------------------------------------------------------- 
 * 
 * File Name: HeapSort.c
 * 
 * Author: Creighton, T. D.
 * 
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------*/

/* <lalLaTeX>

\subsection{Module \texttt{HeapSort.c}}

Sorts, indexes, or ranks vector elements using the heap sort
algorithm.

\subsubsection{Prototypes}
\vspace{0.1in}
\input{HeapSortD}

\subsubsection{Description}

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

\subsubsection{Algorithm}

These routines use the standard heap sort algorithm described in
Sec.~8.3 of Ref.~\cite{ptvf:1992}.

The \verb@SHeapSort()@ and \verb@DHeapSort()@ routines are entirely
in-place, with no auxiliary storage vector.  The \verb@SHeapIndex()@
and \verb@DHeapIndex()@ routines are also technically in-place, but
they require two input vectors (the data vector and the index vector),
and leave the data vector unchanged.  The \verb@SHeapRank()@ and
\verb@DHeapRank()@ routines require two input vectors (the data and
rank vectors), and also allocate a temporary index vector internally;
these routines are therefore the most memory-intensive.  All of these
algorithms are $N\log_2(N)$ algorithms, regardless of the ordering of
the initial dataset.


\subsubsection{Uses}
\begin{verbatim}
I4CreateVector()
I4DestroyVector()
\end{verbatim}

\subsubsection{Notes}

</lalLaTeX> */

#ifndef _LALSTDLIB_H
#include "LALStdlib.h"
#ifndef _LALSTDLIB_H
#define _LALSTDLIB_H
#endif
#endif

#ifndef _AVFACTORIES_H
#include "AVFactories.h"
#ifndef _AVFACTORIES_H
#define _AVFACTORIES_H
#endif
#endif

#ifndef _SORT_H
#include "Sort.h"
#ifndef _SORT_H
#define _SORT_H
#endif
#endif

NRCSID(HEAPSORTC,"$Id$");

/* <lalVerbatim file="HeapSortD"> */
void SHeapSort(Status      *stat,
	       REAL4Vector *vector)
{ /* </lalVerbatim> */
  INT4 i;
  INT4 j;
  INT4 k;
  INT4 n;
  REAL4 temp;
  REAL4 *data;

  INITSTATUS(stat,HEAPSORTC);

  /* Make sure all pointers are valid. */
  ASSERT(vector,stat,SORT_ENUL,SORT_MSGENUL);
  ASSERT(vector->data,stat,SORT_ENUL,SORT_MSGENUL);

  n=vector->length;
  data=vector->data;

  /* A vector of length 0 or 1 is already sorted. */
  if(n<2)
    RETURN(stat);

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
      }else
	k=j+1;
    }
    data[i]=temp;
  }
}


/* <lalVerbatim file="HeapSortD"> */
void SHeapIndex(Status      *stat,
		INT4Vector  *index,
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

  INITSTATUS(stat,HEAPSORTC);

  /* Make sure all pointers are valid, and the index vector is of the
     same length as the data vector. */
  ASSERT(vector,stat,SORT_ENUL,SORT_MSGENUL);
  ASSERT(index,stat,SORT_ENUL,SORT_MSGENUL);
  ASSERT(vector->data,stat,SORT_ENUL,SORT_MSGENUL);
  ASSERT(index->data,stat,SORT_ENUL,SORT_MSGENUL);
  ASSERT(vector->length==index->length,stat,SORT_ELEN,SORT_MSGELEN);

  n=vector->length;
  data=vector->data;

  /* Initialize the index vector. */
  for(i=0,indx=index->data;i<n;i++,indx++)
    *indx=i;
  indx=index->data;

  /* A vector of length 0 or 1 is already sorted. */
  if(n<2)
    RETURN(stat);

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
      }else
	k=j+1;
    }
    indx[i]=itemp;
  }
}


/* <lalVerbatim file="HeapSortD"> */
void SHeapRank(Status      *stat,
	       INT4Vector  *rank,
	       REAL4Vector *vector)
{ /* </lalVerbatim> */
  INT4       i;
  INT4       *indx;
  INT4       *rnk;
  INT4Vector *index=NULL;

  INITSTATUS(stat,HEAPSORTC);
  ATTATCHSTATUSPTR(stat);

  /* Make sure all pointers are valid, and the rank vector is of the
     same length as the data vector. */
  ASSERT(vector,stat,SORT_ENUL,SORT_MSGENUL);
  ASSERT(rank,stat,SORT_ENUL,SORT_MSGENUL);
  ASSERT(vector->data,stat,SORT_ENUL,SORT_MSGENUL);
  ASSERT(rank->data,stat,SORT_ENUL,SORT_MSGENUL);
  ASSERT(vector->length==rank->length,stat,SORT_ELEN,SORT_MSGELEN);

  /* Make the temporary index vector. */
  TRY(I4CreateVector(stat->statusPtr,&index,vector->length),stat);
  TRY(SHeapIndex(stat->statusPtr,index,vector),stat);

  /* Invert to get the rank vector. */
  indx=index->data;
  rnk=rank->data;
  for(i=0;i<vector->length;i++)
    rnk[indx[i]]=i;

  /* Free memory and exit. */
  TRY(I4DestroyVector(stat->statusPtr,&index),stat);
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}


/* <lalVerbatim file="HeapSortD"> */
void DHeapSort(Status      *stat,
	       REAL8Vector *vector)
{ /* </lalVerbatim> */
  INT4 i;
  INT4 j;
  INT4 k;
  INT4 n;
  REAL8 temp;
  REAL8 *data;

  INITSTATUS(stat,HEAPSORTC);

  /* Make sure all pointers are valid. */
  ASSERT(vector,stat,SORT_ENUL,SORT_MSGENUL);
  ASSERT(vector->data,stat,SORT_ENUL,SORT_MSGENUL);

  n=vector->length;
  data=vector->data;

  /* A vector of length 0 or 1 is already sorted. */
  if(n<2)
    RETURN(stat);

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
      }else
	k=j+1;
    }
    data[i]=temp;
  }
}


/* <lalVerbatim file="HeapSortD"> */
void DHeapIndex(Status      *stat,
		INT4Vector  *index,
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

  INITSTATUS(stat,HEAPSORTC);

  /* Make sure all pointers are valid, and the index vector is of the
     same length as the data vector. */
  ASSERT(vector,stat,SORT_ENUL,SORT_MSGENUL);
  ASSERT(index,stat,SORT_ENUL,SORT_MSGENUL);
  ASSERT(vector->data,stat,SORT_ENUL,SORT_MSGENUL);
  ASSERT(index->data,stat,SORT_ENUL,SORT_MSGENUL);
  ASSERT(vector->length==index->length,stat,SORT_ELEN,SORT_MSGELEN);

  n=vector->length;
  data=vector->data;

  /* Initialize the index vector. */
  for(i=0,indx=index->data;i<n;i++,indx++)
    *indx=i;
  indx=index->data;

  /* A vector of length 0 or 1 is already sorted. */
  if(n<2)
    RETURN(stat);

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
      }else
	k=j+1;
    }
    indx[i]=itemp;
  }
}


/* <lalVerbatim file="HeapSortD"> */
void DHeapRank(Status      *stat,
	       INT4Vector  *rank,
	       REAL8Vector *vector)
{ /* </lalVerbatim> */
  INT4       i;
  INT4       *indx;
  INT4       *rnk;
  INT4Vector *index=NULL;

  INITSTATUS(stat,HEAPSORTC);
  ATTATCHSTATUSPTR(stat);

  /* Make sure all pointers are valid, and the rank vector is of the
     same length as the data vector. */
  ASSERT(vector,stat,SORT_ENUL,SORT_MSGENUL);
  ASSERT(rank,stat,SORT_ENUL,SORT_MSGENUL);
  ASSERT(vector->data,stat,SORT_ENUL,SORT_MSGENUL);
  ASSERT(rank->data,stat,SORT_ENUL,SORT_MSGENUL);
  ASSERT(vector->length==rank->length,stat,SORT_ELEN,SORT_MSGELEN);

  /* Make the temporary index vector. */
  TRY(I4CreateVector(stat->statusPtr,&index,vector->length),stat);
  TRY(DHeapIndex(stat->statusPtr,index,vector),stat);

  /* Invert to get the rank vector. */
  indx=index->data;
  rnk=rank->data;
  for(i=0;i<vector->length;i++)
    rnk[indx[i]]=i;

  /* Free memory and exit. */
  TRY(I4DestroyVector(stat->statusPtr,&index),stat);
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}

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


\subsubsection*{Uses}
\begin{verbatim}
LALI4CreateVector()
LALI4DestroyVector()
\end{verbatim}

\subsubsection*{Notes}
\vfill{\footnotesize\input{HeapSortCV}}

</lalLaTeX> */

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/Sort.h>

NRCSID(HEAPSORTC,"$Id$");

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

/************************************ <lalVerbatim file="LALRunningMedianCV">
Author: Somya D. Mohanty, B. Machenschalk
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{LALRunningMedian.c}}
\label{ss:LALRunningMedian.c}

Functions to efficiently calculate running medians

\subsubsection*{Prototypes}
\vspace{0.1in}
\idx[Function]{LALSRunningMedian}
\idx[Function]{LALDRunningMedian}
\input{LALRunningMedianCP}

\idx{LALRunningMedian()}

\subsubsection*{Description}

The routine \verb+LALDRunningMedian()+ calculates the running medians of a
REAL8Sequence. The routine \verb+LALSRunningMedian()+ does the same for a REAL4Sequence.
\verb+input+ ist a REAL4/REAL8Sequence containing the input array, \verb+blocksize+
is the length of the block the medians are calculated of.
With n being the lenght of the input array and b being the blocksize,
the medians array must be a REAL4/REAL8 sequence of length (n-b+1).
\verb+LALDRunningMedian2()+ is a different implentation of the same algorithm. It 
should behave exactly like \verb+LALDRunningMedian()+, but has proven to be a
little faster. Check if it works for you.
\subsubsection*{Algorithm}

For a detailed description of the algorithm see the
LIGO document T-030168-00-D, Somya D. Mohanty:
Efficient Algorithm for computing a Running Median

\subsubsection*{Uses}
\begin{verbatim}
LALCalloc()
LALFree()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALRunningMedianCV}}

</lalLaTeX> */


#include <stdio.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/LALRunningMedian.h>
/* #include "LALRunningMedian.h" */

/* <lalVerbatim file="LALRunningMedianNRCSID"> */
NRCSID( LALRUNNINGMEDIANC, "$Id$" );
/* </lalVerbatim> */


/*----------------------------------
  A structure to store values and indices
  of elements in an array
  -----------------------------------*/
struct rngmed_val_index8 {
  REAL8 data;
  UINT4 index;
};
struct rngmed_val_index4 {
  REAL4 data;
  UINT4 index;
};


static int rngmed_sortindex8(const void *elem1, const void *elem2){
  /*Used in running qsort*/
  
  const struct rngmed_val_index8 *A = elem1;
  const struct rngmed_val_index8 *B = elem2;
  REAL8 data1, data2;
  
  data1=A->data;
  data2=B->data;
  if (data1 < data2)
    return -1;
  else if (data1==data2)
    return 0;
  else
    return 1;
}

static int rngmed_sortindex4(const void *elem1, const void *elem2){
  /*Used in running qsort*/
  
  const struct rngmed_val_index4 *A = elem1;
  const struct rngmed_val_index4 *B = elem2;
  REAL4 data1, data2;
  
  data1=A->data;
  data2=B->data;
  if (data1 < data2)
    return -1;
  else if (data1==data2)
    return 0;
  else
    return 1;
}


/* Used in running qsort, RunningMedian2 version */
static int rngmed_qsortindex8(const void *elem1, const void *elem2){
  struct qsnode{
    REAL8 value;
    UINT4 index;
  };
  
  const struct qsnode *A = elem1;
  const struct qsnode *B = elem2;

  if (B->value > A->value)
    return -1;
  else if (A->value > B->value)
    return 1;
  else if (A->index > B->index)
    return -1;
  else
    return 1;
}

/* <lalVerbatim file="LALRunningMedianCP"> */
void LALDRunningMedian( LALStatus *status,
			REAL8Sequence *medians,
			const REAL8Sequence *input,
			LALRunningMedianPar param)
/* </lalVerbatim> */
{
  /*----------------------------------
    Two types of pointers:
    (a)next_sorted: point to the next node in sorted list
    (b)next_sequence: point to the next node in sequential list
    ------------------------------------*/
  struct node{
    REAL8 data;
    struct node *next_sorted, *next_sequence, *prev_sorted;
    int rank; /*Used for constructing optional output*/
  };
  
  /*----------------------------------
    checks: Array to hold pointers to Checkpoint nodes.
    first_sequence: Pointer to first node of sequential list
    ------------------------------------*/
  struct node **checks = NULL;
  struct node **node_addresses = NULL;    
  struct node *first_sequence = NULL;
  struct node *last_sequence = NULL;
  struct node *currentnode = NULL;
  struct node *previousnode = NULL; 
  struct node *leftnode = NULL;
  struct node *rightnode = NULL;
  struct node *reuse_next_sorted = NULL;
  struct node *reuse_prev_sorted = NULL;
  struct node *dummy_node = NULL;
  struct node *dummy_node1 = NULL;
  struct node *dummy_node2 = NULL;
  UINT4 ncheckpts,stepchkpts;
  UINT4 nextchkptindx,*checks4shift;
  UINT4 nearestchk,midpoint,offset,numberoffsets;
  UINT4 samplecount,counter_chkpt,chkcount=0, shiftcounter=0;
  INT8 k;
  REAL8 nextsample,deletesample,dummy;
  INT4 shift,dummy_int;
  /* for initial qsort */
  REAL8 *sorted_indices;
  struct rngmed_val_index8 *index_block;

  INITSTATUS( status, "LALDRunningMedian", LALRUNNINGMEDIANC );
  
  /* check input parameters */
  /* input must not be NULL */
  ASSERT(input,status,LALRUNNINGMEDIANH_ENULL,LALRUNNINGMEDIANH_MSGENULL);
  /* param.blocksize must be >2 */
  ASSERT(param.blocksize>2,status,LALRUNNINGMEDIANH_EZERO,LALRUNNINGMEDIANH_MSGEZERO);
  /* blocksize must not be larger than input size */
  ASSERT(param.blocksize <= input->length,status,LALRUNNINGMEDIANH_ELARGE,LALRUNNINGMEDIANH_MSGELARGE);
  /* medians must point to a valid sequence of correct size */
  ASSERT(medians,status,LALRUNNINGMEDIANH_EIMED,LALRUNNINGMEDIANH_MSGEIMED);
  ASSERT(medians->length == (input->length - param.blocksize + 1),status,LALRUNNINGMEDIANH_EIMED,LALRUNNINGMEDIANH_MSGEIMED);

  ATTATCHSTATUSPTR( status );

  /*-----------------------------------
    Sort the first block of param.blocksize samples
    ------------------------------------*/
  index_block =(struct rngmed_val_index8 *)LALCalloc(param.blocksize, sizeof(struct rngmed_val_index8));
  if(!index_block) {
    ABORT(status,LALRUNNINGMEDIANH_EMALOC1,LALRUNNINGMEDIANH_MSGEMALOC1);
  }
  for(k=0;k<param.blocksize;k++){
    index_block[k].data=input->data[k];
    index_block[k].index=k;
  }

  qsort(index_block, param.blocksize, sizeof(struct rngmed_val_index8),rngmed_sortindex8);
  
  sorted_indices=(REAL8 *)LALCalloc(param.blocksize,sizeof(REAL8));
  if(!sorted_indices) {
    LALFree(index_block);
    ABORT(status,LALRUNNINGMEDIANH_EMALOC1,LALRUNNINGMEDIANH_MSGEMALOC1);
  }

  for(k=0;k<param.blocksize;k++){
    sorted_indices[k]=index_block[k].index; 
  }

  LALFree(index_block);

  /*----------------------------------
    Indices of checkpoint nodes.
    Number of nodes per checkpoint=floor(sqrt(param.blocksize))
    ------------------------------------*/
  stepchkpts = sqrt(param.blocksize);
  ncheckpts = param.blocksize/stepchkpts;
  checks = (struct node **)LALCalloc(ncheckpts,sizeof(struct node*));
  if(!checks){
    LALFree(sorted_indices);
    ABORT(status,LALRUNNINGMEDIANH_EMALOC2,LALRUNNINGMEDIANH_MSGEMALOC2);
  }
  checks4shift = (INT4*)LALCalloc(ncheckpts,sizeof(INT4));
  if(!checks4shift){
    LALFree(sorted_indices);
    LALFree(checks);
    ABORT(status,LALRUNNINGMEDIANH_EMALOC3,LALRUNNINGMEDIANH_MSGEMALOC3);
  }

  /*---------------------------------
    Offsets for getting median from nearest
    checkpoint: For param.blocksize even, 
    (node(offset(1))+node(offset(2)))/2;
    for param.blocksize odd,
    (node(offset(1))+node(offset(1)))/2;
    ----------------------------------*/
  if((int)fmod(param.blocksize,2.0)){
    /*Odd*/
    midpoint=(param.blocksize+1)/2-1;
    numberoffsets=1;
  }
  else{
    /*Even*/
    midpoint=param.blocksize/2-1;
    numberoffsets=2;   
  }
  nearestchk=floor(midpoint/stepchkpts);
  offset=midpoint-nearestchk*stepchkpts;

  /*----------------------------------
    Build up linked list using first nblock points
    in sequential order
    ------------------------------------*/
  node_addresses=(struct node **)LALCalloc(param.blocksize,sizeof(struct node *));
  if(!node_addresses){
    LALFree(sorted_indices);
    LALFree(checks4shift);
    LALFree(checks);
    ABORT(status,LALRUNNINGMEDIANH_EMALOC4,LALRUNNINGMEDIANH_MSGEMALOC4);
  }
  first_sequence=(struct node *)LALCalloc(1,sizeof(struct node));
  if(!first_sequence){
    LALFree(node_addresses);
    LALFree(sorted_indices);
    LALFree(checks4shift);
    LALFree(checks);
    ABORT(status,LALRUNNINGMEDIANH_EMALOC5,LALRUNNINGMEDIANH_MSGEMALOC5);
  }
  node_addresses[0]=first_sequence;
  first_sequence->next_sequence=NULL;
  first_sequence->next_sorted=NULL;
  first_sequence->prev_sorted=NULL;
  first_sequence->data=input->data[0];
  previousnode=first_sequence;
  for(samplecount=1;samplecount<param.blocksize;samplecount++){
    currentnode=(struct node *)LALCalloc(1,sizeof(struct node));
    if(!currentnode){
      LALFree(first_sequence);
      LALFree(sorted_indices);
      LALFree(node_addresses);
      LALFree(checks4shift);
      LALFree(checks);
      ABORT(status,LALRUNNINGMEDIANH_EMALOC6,LALRUNNINGMEDIANH_MSGEMALOC6);
    }
    node_addresses[samplecount]=currentnode;
    previousnode->next_sequence=currentnode;
    currentnode->next_sequence=NULL;
    currentnode->prev_sorted=NULL;
    currentnode->next_sorted=NULL;
    currentnode->data=input->data[samplecount];
    previousnode=currentnode;
  }
  last_sequence=currentnode;

  /*------------------------------------
    Set the sorted sequence pointers and
    the pointers to checkpoint nodes
    -------------------------------------*/
  currentnode=node_addresses[(int)sorted_indices[0]];
  previousnode=NULL;
  checks[0]=currentnode;
  nextchkptindx=stepchkpts;
  counter_chkpt=1;
  for(samplecount=1;samplecount<param.blocksize;samplecount++){
    dummy_node=node_addresses[(int)sorted_indices[samplecount]];
    currentnode->next_sorted=dummy_node;
    currentnode->prev_sorted=previousnode;
    previousnode=currentnode;
    currentnode=dummy_node;
    if(samplecount==nextchkptindx && counter_chkpt<ncheckpts){
      checks[counter_chkpt]=currentnode;
      nextchkptindx+=stepchkpts;
      counter_chkpt++;
    }
  }
  currentnode->prev_sorted=previousnode;
  currentnode->next_sorted=NULL;
  LALFree(sorted_indices);

  /*------------------------------
    get the first output element
    -------------------------------*/
  currentnode=checks[nearestchk];
  for(k=1;k<=offset;k++){
    currentnode=currentnode->next_sorted;
  }
  dummy=0;
  for(k=1;k<=numberoffsets;k++){
    dummy+=currentnode->data;
    currentnode=currentnode->next_sorted;
  }
  medians->data[0]=dummy/numberoffsets;

  /*---------------------------------
    This is the main part
    ----------------------------------*/
  for(samplecount=param.blocksize;samplecount<input->length;samplecount++){
    nextsample=input->data[samplecount];
    if(nextsample>=checks[0]->data){ /* corrected */
      for(chkcount=1;chkcount<ncheckpts;chkcount++){
	if(nextsample>checks[chkcount]->data){
	}
	else{
	  break;
	}
      }
      chkcount-=1;
      rightnode=checks[chkcount];
      leftnode=NULL;  /* corrected */
      while(rightnode){
	if(nextsample<rightnode->data){ /* corrected */
	  break;
	}
	leftnode=rightnode;
	rightnode=rightnode->next_sorted;
      }
    } else {
      /* corrected */
      /* new case */
      if(nextsample<checks[0]->data){  
	chkcount=0;
	/* dummy_node=checks[0]; */
	rightnode=checks[0];
	leftnode=NULL;  
      }
    }

    /* corrected */
    /* from here ... */
    dummy_node=NULL;
    if(rightnode==first_sequence){
      dummy_node=rightnode;
    }
    else if(leftnode==first_sequence){
      dummy_node=leftnode;
    }
    if(dummy_node) {
      dummy_node->data=nextsample;
      first_sequence=first_sequence->next_sequence;
      dummy_node->next_sequence=NULL;
      last_sequence->next_sequence=dummy_node;
      last_sequence=dummy_node;
      shift=0;                
    }
    else{
      reuse_next_sorted=rightnode;
      reuse_prev_sorted=leftnode;
      shift=1; /*shift maybe required*/
    }
    /* to here */

    /*-----------------------------------
      Getting check points to be shifted
      -----------------------------------*/
    if(shift){
      deletesample=first_sequence->data;
      if(deletesample>nextsample){
	shiftcounter=0;
	for(k=chkcount;k<ncheckpts;k++){
	  dummy=checks[k]->data;
	  if(dummy>nextsample){ /* corrected */
	    if(dummy<=deletesample){
	      checks4shift[shiftcounter]=k;
	      shiftcounter++;
	    }
	    else{
	      break;
	    }
	  }
	}
	shift=-1; /*Left shift*/
      }
      else 
	if(deletesample<nextsample){
	  shiftcounter=0;
	  for(k=chkcount;k>=0;k--){
	    dummy=checks[k]->data;
	    if(dummy>=deletesample){
	      checks4shift[shiftcounter]=k;
	      shiftcounter++;
	    }
	    else{
	      break;
	    }
	  }
	  shift=1; /*Shift Right*/
	}
      /* corrected: else case deleted */
    }
    
    /*------------------------------
      Delete and Insert
      --------------------------------*/
    if(shift){
      dummy_node=first_sequence;
      first_sequence=dummy_node->next_sequence;
      dummy_node->next_sequence=NULL;
      last_sequence->next_sequence=dummy_node;
      last_sequence=dummy_node;
      dummy_node->data=nextsample;
      dummy_node1=dummy_node->prev_sorted;
      dummy_node2=dummy_node->next_sorted;
      
      /*-----------------------
	Repair deletion point
	------------------------*/
      if(!dummy_node1){
	dummy_node2->prev_sorted=dummy_node1;
      }
      else {
	if(!dummy_node2){
	  dummy_node1->next_sorted=dummy_node2;
	}
	else{
	  dummy_node1->next_sorted=dummy_node2;
	  dummy_node2->prev_sorted=dummy_node1;
	}
      }  
      
      /*------------------------
	Set pointers from neighbours to new node at insertion point
	-------------------------*/
      if(!rightnode){
	leftnode->next_sorted=dummy_node;
      }
      else {
	if(!leftnode){
	  rightnode->prev_sorted=dummy_node;
	}
	else{
	  leftnode->next_sorted=dummy_node;
	  rightnode->prev_sorted=dummy_node;
	}
      }
      
      /*-------------------------------
	Shift check points before resetting sorted list
	--------------------------------*/
      if(shift==-1){
	for(k=0;k<shiftcounter;k++){
	  dummy_int=checks4shift[k];
	  checks[dummy_int]=checks[dummy_int]->prev_sorted;
	}
      }
      else
	if(shift==1){
	  for(k=0;k<shiftcounter;k++){
	    dummy_int=checks4shift[k];
	    checks[dummy_int]=checks[dummy_int]->next_sorted;
	  } 
	}
      
      /*--------------------------------
	insert node
	--------------------------------*/
      dummy_node->next_sorted=reuse_next_sorted;
      dummy_node->prev_sorted=reuse_prev_sorted;
    }

    /*--------------------------------
      Get the median
      ---------------------------------*/
    currentnode=checks[nearestchk];
         for(k=1;k<=offset;k++){
	   currentnode=currentnode->next_sorted;
         }
         dummy=0;
         for(k=1;k<=numberoffsets;k++){
	   dummy+=currentnode->data;
	   currentnode=currentnode->next_sorted;
         }
         medians->data[samplecount-param.blocksize+1]=dummy/numberoffsets;
  }/*Outer For Loop*/
  

  /*--------------------------------
    Clean Up
    ---------------------------------*/
  LALFree(node_addresses);
  currentnode=first_sequence;
  while(currentnode){
    previousnode=currentnode;
    currentnode=currentnode->next_sequence;
    LALFree(previousnode);
  }
  LALFree(checks4shift);
  LALFree(checks);
  
  DETATCHSTATUSPTR( status );
  RETURN( status );  
}

/* <lalVerbatim file="LALRunningMedianCP"> */
void LALSRunningMedian( LALStatus *status,
			REAL4Sequence *medians,
			const REAL4Sequence *input,
			LALRunningMedianPar param)
/* </lalVerbatim> */
{
  /*----------------------------------
    Two types of pointers:
    (a)next_sorted: point to the next node in sorted list
    (b)next_sequence: point to the next node in sequential list
    ------------------------------------*/
  struct node{
    REAL4 data;
    struct node *next_sorted, *next_sequence, *prev_sorted;
    int rank; /*Used for constructing optional output*/
  };
  
  /*----------------------------------
    checks: Array to hold pointers to Checkpoint nodes.
    first_sequence: Pointer to first node of sequential list
    ------------------------------------*/
  struct node **checks = NULL;
  struct node **node_addresses = NULL;    
  struct node *first_sequence = NULL;
  struct node *last_sequence = NULL;
  struct node *currentnode = NULL;
  struct node *previousnode = NULL; 
  struct node *leftnode = NULL;
  struct node *rightnode = NULL;
  struct node *reuse_next_sorted = NULL;
  struct node *reuse_prev_sorted = NULL;
  struct node *dummy_node = NULL;
  struct node *dummy_node1 = NULL;
  struct node *dummy_node2 = NULL;
  UINT4 ncheckpts,stepchkpts;
  UINT4 nextchkptindx,*checks4shift;
  UINT4 nearestchk,midpoint,offset,numberoffsets;
  UINT4 samplecount,counter_chkpt,chkcount=0,shiftcounter=0;
  INT8 k;
  REAL4 nextsample,deletesample,dummy;
  INT4 shift,dummy_int;
  /* for initial qsort */
  REAL4 *sorted_indices;
  struct rngmed_val_index4 *index_block;


  INITSTATUS( status, "LALSRunningMedian", LALRUNNINGMEDIANC );
  
  /* check input parameters */
  /* input must not be NULL */
  ASSERT(input,status,LALRUNNINGMEDIANH_ENULL,LALRUNNINGMEDIANH_MSGENULL);
  /* param.blocksize must be >2 */
  ASSERT(param.blocksize>2,status,LALRUNNINGMEDIANH_EZERO,LALRUNNINGMEDIANH_MSGEZERO);
  /* param.blocksize must not be larger than input size */
  ASSERT(param.blocksize <= input->length,status,LALRUNNINGMEDIANH_ELARGE,LALRUNNINGMEDIANH_MSGELARGE);
  /* medians must point to a valid sequence of correct size */
  ASSERT(medians,status,LALRUNNINGMEDIANH_EIMED,LALRUNNINGMEDIANH_MSGEIMED);
  ASSERT(medians->length == (input->length - param.blocksize + 1),status,LALRUNNINGMEDIANH_EIMED,LALRUNNINGMEDIANH_MSGEIMED);

  ATTATCHSTATUSPTR( status );

  /*-----------------------------------
    Sort the first block of param.blocksize samples
    ------------------------------------*/
  index_block =(struct rngmed_val_index4 *)LALCalloc(param.blocksize, sizeof(struct rngmed_val_index4));
  if(!index_block) {
    ABORT(status,LALRUNNINGMEDIANH_EMALOC1,LALRUNNINGMEDIANH_MSGEMALOC1);
  }
  for(k=0;k<param.blocksize;k++){
    index_block[k].data=input->data[k];
    index_block[k].index=k;
  }

  qsort(index_block, param.blocksize, sizeof(struct rngmed_val_index4),rngmed_sortindex4);
  
  sorted_indices=(REAL4 *)LALCalloc(param.blocksize,sizeof(REAL4));
  if(!sorted_indices) {
    LALFree(index_block);
    ABORT(status,LALRUNNINGMEDIANH_EMALOC1,LALRUNNINGMEDIANH_MSGEMALOC1);
  }

  for(k=0;k<param.blocksize;k++){
    sorted_indices[k]=index_block[k].index; 
  }

  LALFree(index_block);

  /*----------------------------------
    Indices of checkpoint nodes.
    Number of nodes per checkpoint=floor(sqrt(param.blocksize))
    ------------------------------------*/
  stepchkpts = sqrt(param.blocksize);
  ncheckpts = param.blocksize/stepchkpts;
  checks = (struct node **)LALCalloc(ncheckpts,sizeof(struct node*));
  if(!checks){
    LALFree(sorted_indices);
    ABORT(status,LALRUNNINGMEDIANH_EMALOC2,LALRUNNINGMEDIANH_MSGEMALOC2);
  }
  checks4shift = (INT4*)LALCalloc(ncheckpts,sizeof(INT4));
  if(!checks4shift){
    LALFree(sorted_indices);
    LALFree(checks);
    ABORT(status,LALRUNNINGMEDIANH_EMALOC3,LALRUNNINGMEDIANH_MSGEMALOC3);
  }

  /*---------------------------------
    Offsets for getting median from nearest
    checkpoint: For param.blocksize even, 
    (node(offset(1))+node(offset(2)))/2;
    for param.blocksize odd,
    (node(offset(1))+node(offset(1)))/2;
    ----------------------------------*/
  if((int)fmod(param.blocksize,2.0)){
    /*Odd*/
    midpoint=(param.blocksize+1)/2-1;
    numberoffsets=1;
  }
  else{
    /*Even*/
    midpoint=param.blocksize/2-1;
    numberoffsets=2;   
  }
  nearestchk=floor(midpoint/stepchkpts);
  offset=midpoint-nearestchk*stepchkpts;

  /*----------------------------------
    Build up linked list using first nblock points
    in sequential order
    ------------------------------------*/
  node_addresses=(struct node **)LALCalloc(param.blocksize,sizeof(struct node *));
  if(!node_addresses){
    LALFree(sorted_indices);
    LALFree(checks4shift);
    LALFree(checks);
    ABORT(status,LALRUNNINGMEDIANH_EMALOC4,LALRUNNINGMEDIANH_MSGEMALOC4);
  }
  first_sequence=(struct node *)LALCalloc(1,sizeof(struct node));
  if(!first_sequence){
    LALFree(node_addresses);
    LALFree(sorted_indices);
    LALFree(checks4shift);
    LALFree(checks);
    ABORT(status,LALRUNNINGMEDIANH_EMALOC5,LALRUNNINGMEDIANH_MSGEMALOC5);
  }
  node_addresses[0]=first_sequence;
  first_sequence->next_sequence=NULL;
  first_sequence->next_sorted=NULL;
  first_sequence->prev_sorted=NULL;
  first_sequence->data=input->data[0];
  previousnode=first_sequence;
  for(samplecount=1;samplecount<param.blocksize;samplecount++){
    currentnode=(struct node *)LALCalloc(1,sizeof(struct node));
    if(!currentnode){
      LALFree(first_sequence);
      LALFree(sorted_indices);
      LALFree(node_addresses);
      LALFree(checks4shift);
      LALFree(checks);
      ABORT(status,LALRUNNINGMEDIANH_EMALOC6,LALRUNNINGMEDIANH_MSGEMALOC6);
    }
    node_addresses[samplecount]=currentnode;
    previousnode->next_sequence=currentnode;
    currentnode->next_sequence=NULL;
    currentnode->prev_sorted=NULL;
    currentnode->next_sorted=NULL;
    currentnode->data=input->data[samplecount];
    previousnode=currentnode;
  }
  last_sequence=currentnode;

  /*------------------------------------
    Set the sorted sequence pointers and
    the pointers to checkpoint nodes
    -------------------------------------*/
  currentnode=node_addresses[(int)sorted_indices[0]];
  previousnode=NULL;
  checks[0]=currentnode;
  nextchkptindx=stepchkpts;
  counter_chkpt=1;
  for(samplecount=1;samplecount<param.blocksize;samplecount++){
    dummy_node=node_addresses[(int)sorted_indices[samplecount]];
    currentnode->next_sorted=dummy_node;
    currentnode->prev_sorted=previousnode;
    previousnode=currentnode;
    currentnode=dummy_node;
    if(samplecount==nextchkptindx && counter_chkpt<ncheckpts){
      checks[counter_chkpt]=currentnode;
      nextchkptindx+=stepchkpts;
      counter_chkpt++;
    }
  }
  currentnode->prev_sorted=previousnode;
  currentnode->next_sorted=NULL;
  LALFree(sorted_indices);

  /*------------------------------
    get the first output element
    -------------------------------*/
  currentnode=checks[nearestchk];
  for(k=1;k<=offset;k++){
    currentnode=currentnode->next_sorted;
  }
  dummy=0;
  for(k=1;k<=numberoffsets;k++){
    dummy+=currentnode->data;
    currentnode=currentnode->next_sorted;
  }
  medians->data[0]=dummy/numberoffsets;

  /*---------------------------------
    This is the main part
    ----------------------------------*/
  for(samplecount=param.blocksize;samplecount<input->length;samplecount++){
    nextsample=input->data[samplecount];
    if(nextsample>=checks[0]->data){ /* corrected */
      for(chkcount=1;chkcount<ncheckpts;chkcount++){
	if(nextsample>checks[chkcount]->data){
	}
	else{
	  break;
	}
      }
      chkcount-=1;
      rightnode=checks[chkcount];
      leftnode=NULL;  /* corrected */
      while(rightnode){
	if(nextsample<rightnode->data){ /* corrected */
	  break;
	}
	leftnode=rightnode;
	rightnode=rightnode->next_sorted;
      }
    } else {
      /* corrected */
      /* new case */
      if(nextsample<checks[0]->data){  
	chkcount=0;
	/* dummy_node=checks[0]; */
	rightnode=checks[0];
	leftnode=NULL;  
      }
    }

    /* corrected */
    /* from here ... */
    dummy_node=NULL;
    if(rightnode==first_sequence){
      dummy_node=rightnode;
    }
    else if(leftnode==first_sequence){
      dummy_node=leftnode;
    }
    if(dummy_node) {
      dummy_node->data=nextsample;
      first_sequence=first_sequence->next_sequence;
      dummy_node->next_sequence=NULL;
      last_sequence->next_sequence=dummy_node;
      last_sequence=dummy_node;
      shift=0;                
    }
    else{
      reuse_next_sorted=rightnode;
      reuse_prev_sorted=leftnode;
      shift=1; /*shift maybe required*/
    }
    /* to here */

    /*-----------------------------------
      Getting check points to be shifted
      -----------------------------------*/
    if(shift){
      deletesample=first_sequence->data;
      if(deletesample>nextsample){
	shiftcounter=0;
	for(k=chkcount;k<ncheckpts;k++){
	  dummy=checks[k]->data;
	  if(dummy>nextsample){ /* corrected */
	    if(dummy<=deletesample){
	      checks4shift[shiftcounter]=k;
	      shiftcounter++;
	    }
	    else{
	      break;
	    }
	  }
	}
	shift=-1; /*Left shift*/
      }
      else 
	if(deletesample<nextsample){
	  shiftcounter=0;
	  for(k=chkcount;k>=0;k--){
	    dummy=checks[k]->data;
	    if(dummy>=deletesample){
	      checks4shift[shiftcounter]=k;
	      shiftcounter++;
	    }
	    else{
	      break;
	    }
	  }
	  shift=1; /*Shift Right*/
	}
      /* corrected: else case deleted */
    }
    
    /*------------------------------
      Delete and Insert
      --------------------------------*/
    if(shift){
      dummy_node=first_sequence;
      first_sequence=dummy_node->next_sequence;
      dummy_node->next_sequence=NULL;
      last_sequence->next_sequence=dummy_node;
      last_sequence=dummy_node;
      dummy_node->data=nextsample;
      dummy_node1=dummy_node->prev_sorted;
      dummy_node2=dummy_node->next_sorted;
      
      /*-----------------------
	Repair deletion point
	------------------------*/
      if(!dummy_node1){
	dummy_node2->prev_sorted=dummy_node1;
      }
      else {
	if(!dummy_node2){
	  dummy_node1->next_sorted=dummy_node2;
	}
	else{
	  dummy_node1->next_sorted=dummy_node2;
	  dummy_node2->prev_sorted=dummy_node1;
	}
      }  
      
      /*------------------------
	Set pointers from neighbours to new node at insertion point
	-------------------------*/
      if(!rightnode){
	leftnode->next_sorted=dummy_node;
      }
      else {
	if(!leftnode){
	  rightnode->prev_sorted=dummy_node;
	}
	else{
	  leftnode->next_sorted=dummy_node;
	  rightnode->prev_sorted=dummy_node;
	}
      }
      
      /*-------------------------------
	Shift check points before resetting sorted list
	--------------------------------*/
      if(shift==-1){
	for(k=0;k<shiftcounter;k++){
	  dummy_int=checks4shift[k];
	  checks[dummy_int]=checks[dummy_int]->prev_sorted;
	}
      }
      else
	if(shift==1){
	  for(k=0;k<shiftcounter;k++){
	    dummy_int=checks4shift[k];
	    checks[dummy_int]=checks[dummy_int]->next_sorted;
	  } 
	}
      
      /*--------------------------------
	insert node
	--------------------------------*/
      dummy_node->next_sorted=reuse_next_sorted;
      dummy_node->prev_sorted=reuse_prev_sorted;
    }

    /*--------------------------------
      Get the median
      ---------------------------------*/
    currentnode=checks[nearestchk];
         for(k=1;k<=offset;k++){
	   currentnode=currentnode->next_sorted;
         }
         dummy=0;
         for(k=1;k<=numberoffsets;k++){
	   dummy+=currentnode->data;
	   currentnode=currentnode->next_sorted;
         }
         medians->data[samplecount-param.blocksize+1]=dummy/numberoffsets;
  }/*Outer For Loop*/
  

  /*--------------------------------
    Clean Up
    ---------------------------------*/
  LALFree(node_addresses);
  currentnode=first_sequence;
  while(currentnode){
    previousnode=currentnode;
    currentnode=currentnode->next_sequence;
    LALFree(previousnode);
  }
  LALFree(checks4shift);
  LALFree(checks);
  
  DETATCHSTATUSPTR( status );
  RETURN( status );  
}

/* <lalVerbatim file="LALRunningMedianCP"> */
void LALDRunningMedian2( LALStatus *status,
			 REAL8Sequence *medians,
			 const REAL8Sequence *input,
			 LALRunningMedianPar param)
/* </lalVerbatim> */
{
  /* a single "node"
   lesser  points to the next node with less or equal value
   greater points to the next node with greater or equal value
   an index == blocksize is an end marker 
  */
  struct node{
    REAL8 value;
    UINT4 lesser;
    UINT4 greater;
  };

  /* a node of the quicksort array */
  struct qsnode{
    REAL8 value;
    UINT4 index;
  };

  const UINT4 bsize = param.blocksize; /* just an abbrevation */
  const UINT4 nil = bsize;       /* invalid index used as end marker */
  const BOOLEAN isodd = bsize&1; /* bsize is odd = median is a single element */

  struct node* nodes;           /* array of nodes, will be of size blocksize */
  struct qsnode* qsnodes;       /* array of indices for initial qsort */
  UINT4* checkpts;              /* array of checkpoints */
  UINT4  ncheckpts,stepchkpts;  /* checkpoints: number and distance between */
  UINT4  oldestnode;            /* index of "oldest" node */
  UINT4  i;                     /* loop counter (up to input length) */
  INT4   j;                     /* loop counter (might get negative) */
  UINT4  nmedian;               /* current median, outer loop counter */
  UINT4  midpoint;              /* index of middle node in sorting order */
  UINT4  mdnnearest;            /* checkpoint "nearest" to the median */
  UINT4  nextnode;              /* node after an insertion point,
			           also used to find a median */
  UINT4  prevnode;              /* node before an insertion point */
  UINT4  rightcheckpt;          /* checkpoint 'right' of an insertion point */
  REAL8 oldvalue,newvalue;      /* old + new value of the node being replaced */
  UINT4 oldlesser,oldgreater;   /* remember the pointers of the replaced node */

  INITSTATUS( status, "LALDRunningMedian", LALRUNNINGMEDIANC );

  /* check input parameters */
  /* input must not be NULL */
  ASSERT(input,status,LALRUNNINGMEDIANH_ENULL,LALRUNNINGMEDIANH_MSGENULL);
  /* param.blocksize must be >2 */
  ASSERT(param.blocksize>2,
	 status,LALRUNNINGMEDIANH_EZERO,LALRUNNINGMEDIANH_MSGEZERO);
  /* blocksize must not be larger than input size */
  ASSERT(param.blocksize <= input->length,
	 status,LALRUNNINGMEDIANH_ELARGE,LALRUNNINGMEDIANH_MSGELARGE);
  /* medians must point to a valid sequence of correct size */
  ASSERT(medians,status,LALRUNNINGMEDIANH_EIMED,LALRUNNINGMEDIANH_MSGEIMED);
  ASSERT(medians->length == (input->length - param.blocksize + 1),
	 status,LALRUNNINGMEDIANH_EIMED,LALRUNNINGMEDIANH_MSGEIMED);

  ATTATCHSTATUSPTR( status );

  /* create nodes array */
  nodes = (struct node*)LALCalloc(bsize, sizeof(struct node));  

  /* determine checkpoint positions */
  stepchkpts = sqrt(bsize);
  /* the old form
     ncheckpts = bsize/stepchkpts;
     caused too less checkpoints at the end, leading to break the
     cost calculation */
  ncheckpts = ceil((REAL4)bsize/(REAL4)stepchkpts);

  /* set checkpoint nearest to the median and offset of the median to it */
  midpoint = (bsize+(bsize&1)) / 2 - 1;
  /* this becomes the median checkpoint */
  mdnnearest = ceil((REAL4)midpoint / (REAL4)stepchkpts);

  /* add a checkpoint for the median if necessary */
  if (ceil((REAL4)midpoint / (REAL4)stepchkpts) != (REAL4)midpoint / (REAL4)stepchkpts)
    ncheckpts++;

  /* create checkpoints array */
  checkpts = (UINT4*)LALCalloc(ncheckpts,sizeof(UINT4));

  /* create array for qsort */
  qsnodes = (struct qsnode*)LALCalloc(bsize, sizeof(struct qsnode));  

  /* init qsort array
   the nodes get their values from the input,
   the indices are only identities qi[0]=0,qi[1]=1,... */
  for(i=0;i<bsize;i++) {
    qsnodes[i].value = input->data[i];
    qsnodes[i].index = i;
  }

  /* sort qsnodes by value and index(!) */
  qsort(qsnodes, bsize, sizeof(struct qsnode),rngmed_qsortindex8);

  /* init nodes array */
  for(i=0;i<bsize;i++)
    nodes[i].value = input->data[i];
  for(i=1;i<bsize-1;i++) {
    nodes[qsnodes[i-1].index].greater = qsnodes[i].index;
    nodes[qsnodes[i+1].index].lesser  = qsnodes[i].index;
  }
  nodes[qsnodes[0].index].lesser      = nil; /* end marker */
  nodes[qsnodes[1].index].lesser      = qsnodes[0].index;
  nodes[qsnodes[bsize-2].index].greater = qsnodes[bsize-1].index;
  nodes[qsnodes[bsize-1].index].greater = nil; /* end marker */

  /* setup checkpoints */
  /* j is the current checkpoint
     i is the stepping
     they are out of sync after a median checkpoint has been added */
  for(i=0,j=0; (UINT4)j<ncheckpts; i++,j++) {
    if ((UINT4)j == mdnnearest) {
      checkpts[j] = qsnodes[midpoint].index;
      if (i*stepchkpts != midpoint)
	j++;
    }
    checkpts[j] = qsnodes[i*stepchkpts].index;
  }

  /* don't need the qsnodes anymore */
  LALFree(qsnodes);

  /* find first median */
  nextnode = checkpts[mdnnearest];
  if(isodd)
    medians->data[0] = nodes[nextnode].value;
  else
    medians->data[0] = (nodes[nextnode].value
			+ nodes[nodes[nextnode].greater].value) / 2.0;

  /* the "oldest" node (first in sequence) is the one with index 0 */
  oldestnode = 0;

  /* outer loop: find a median with each iteration */
  for(nmedian=1; nmedian < medians->length; nmedian++) {

    /* remember value of sample to be deleted */
    oldvalue = nodes[oldestnode].value;

    /* get next value to be inserted from input */
    newvalue = input->data[nmedian+bsize-1];

    /** find point of insertion: **/

    /* find checkpoint greater or equal newvalue */
    /* possible optimisation: use bisectional search instaed of linear */
    for(rightcheckpt=0; rightcheckpt<ncheckpts; rightcheckpt++)
      if(newvalue <= nodes[checkpts[rightcheckpt]].value)
	break;
      
    /* assume we are inserting at the beginning: */
    prevnode = nil;
    if (rightcheckpt == 0)
      /* yes, we are */
      nextnode = checkpts[0];
    else {
      /* we're beyond the first checkpoint, find the node we're inserting at: */
      nextnode = checkpts[rightcheckpt-1]; /* this also works if we found no
					      checkpoint > newvalue, as
					      then rightcheckpt == ncheckpts */ 
      /* the following loop is always ran at least once, as
	 nodes[checkpts[rightcheckpt-1]].value < newvalue
         after 'find checkpoint' loop */ 
      while((nextnode != nil) && (newvalue > nodes[nextnode].value)) {
	prevnode = nextnode;
	nextnode = nodes[nextnode].greater;
      }
    }
    /* ok, we have:
       - case 1: insert at beginning: prevnode == nil, nextnode == smallest node
       - case 2: insert at end: nextnode == nil (terminated loop),
                 prevnode == last node
       - case 3: ordinary insert: insert between prevnode and nextnode
    */

    /* insertion deletion and shifting are unnecessary if we are replacing
       at the same pos */
    if ((oldestnode != prevnode) && (oldestnode != nextnode)) {

      /* delete oldest node from list */
      if (nodes[oldestnode].lesser == nil) {
	/* case 1: at beginning */
	nodes[nodes[oldestnode].greater].lesser = nil;
	/* this shouldn't be necessary, but doesn't harm */
	checkpts[0] = nodes[oldestnode].greater;
      } else if (nodes[oldestnode].greater == nil)
	/* case 2: at end */
	nodes[nodes[oldestnode].lesser].greater = nil;
      else {
	/* case 3: anywhere else */
	nodes[nodes[oldestnode].lesser].greater = nodes[oldestnode].greater;
	nodes[nodes[oldestnode].greater].lesser = nodes[oldestnode].lesser;
      }
      /* remember the old links for special case in shifting below */
      oldgreater = nodes[oldestnode].greater;
      oldlesser = nodes[oldestnode].lesser;


      /* insert new node - actually we reuse the oldest one */
      /* the value is set outside the outer "if" */
      nodes[oldestnode].lesser = prevnode;
      nodes[oldestnode].greater = nextnode;
      if (prevnode != nil)
	nodes[prevnode].greater = oldestnode;
      if (nextnode != nil)
	nodes[nextnode].lesser = oldestnode;


      /* shift checkpoints */

      /* if there is a sequence of identical values, new values are inserted
	 always at the left end. Thus, the oldest value has to be the rightmost
	 of such a sequence. This requires proper init.
	 
	 This makes shifting of the checkpoints rather easy:
	 if (oldvalue < newvalue), all checkpoints with
	     oldvalue <(=) chkptvalue < newvalue are shifted,
	 if (newvalue <= oldvalue), all checkpoints with
	     newvalue <= chkptvalue <= oldvalue are shifted.
	 <(=) means that only a checkpoint at the deleted node must be
	     shifted, no other accidently pointing to the same value.

	 Care is needed if a checkpoint to shift is the node we just deleted

	 We start at the checkpoint we know to be closest to the new node
	 satifying the above condition:
	 rightcheckpt-1 if (oldvalue < newvalue)
	 rightcheckpt othewise
	 and proceed in the direction towards the deleted node
      */

      if (oldvalue < newvalue) {
	/* we shift towards larger values */
	for(j=rightcheckpt-1; (j>0)&&(nodes[checkpts[j]].value >= oldvalue);j--)
	  if (nodes[checkpts[j]].value > oldvalue)
	    checkpts[j] = nodes[checkpts[j]].greater;
	  else if (checkpts[j] == oldestnode)
	    checkpts[j] = oldgreater;
      } else /* newvalue <= oldvalue */
	/* we shift towards smaller values */
	for(i=rightcheckpt;
	    (i<ncheckpts) && (nodes[checkpts[i]].value <= oldvalue); i++)
	  if (checkpts[i] == oldestnode)
	    checkpts[i] = oldlesser;
	  else
	    checkpts[i] = nodes[checkpts[i]].lesser;

    } /* if ((oldestnode != prevnode) && (oldestnode != nextnode)) */

    /* in any case set new value */
    nodes[oldestnode].value = newvalue;


    /* find median */
    if (newvalue == oldvalue)
      medians->data[nmedian] = medians->data[nmedian-1];
    else {
      nextnode = checkpts[mdnnearest];
      if(isodd)
	medians->data[nmedian] = nodes[nextnode].value;
      else
	medians->data[nmedian] = (nodes[nextnode].value
				  + nodes[nodes[nextnode].greater].value) / 2.0;
    }
    
    /* next oldest node */
    oldestnode = (oldestnode + 1) % bsize; /* wrap around */
    
  } /* for (nmedian...) */
  
  /* cleanup */
  LALFree(checkpts);
  LALFree(nodes);

  DETATCHSTATUSPTR( status );
  RETURN( status );  
}


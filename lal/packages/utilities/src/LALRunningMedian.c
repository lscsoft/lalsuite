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
  
  struct rngmed_val_index8 *A, *B;
  REAL8 data1, data2;
  
  A=(struct rngmed_val_index8 *)elem1;
  B=(struct rngmed_val_index8 *)elem2;
  
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
  
  struct rngmed_val_index4 *A, *B;
  REAL4 data1, data2;
  
  A=(struct rngmed_val_index4 *)elem1;
  B=(struct rngmed_val_index4 *)elem2;
  
  data1=A->data;
  data2=B->data;
  if (data1 < data2)
    return -1;
  else if (data1==data2)
    return 0;
  else
    return 1;
  
}


/* <lalVerbatim file="LALRunningMedianCP"> */
void LALDRunningMedian( LALStatus *status,
			REAL8Sequence *medians,
			REAL8Sequence *input,
			LALRunningMedianPar param)
/* </lalVerbatim> */
{
  REAL8 *sorted_indices;

  struct rngmed_val_index8 *index_block;

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
  UINT4 samplecount,counter_chkpt,chkcount,shiftcounter = 0;
  INT8 k;
  REAL8 nextsample,deletesample,dummy;
  UINT4 shift,dummy_int;

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
    ABORT(status,LALRUNNINGMEDIANH_EMALOC2,LALRUNNINGMEDIANH_MSGEMALOC2);
  }
  checks4shift = (INT4*)LALCalloc(ncheckpts,sizeof(INT4));
  if(!checks4shift){
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
    LALFree(checks4shift);
    LALFree(checks);
    ABORT(status,LALRUNNINGMEDIANH_EMALOC4,LALRUNNINGMEDIANH_MSGEMALOC4);
  }
  first_sequence=(struct node *)LALCalloc(1,sizeof(struct node));
  if(!first_sequence){
    LALFree(node_addresses);
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
    if(nextsample>checks[0]->data){
      for(chkcount=1;chkcount<ncheckpts;chkcount++){
	if(nextsample>checks[chkcount]->data){
	}
	else{
	  break;
	}
      }
      chkcount-=1;
      rightnode=checks[chkcount];
      while(rightnode){
	if(nextsample<=rightnode->data){
	  break;
	}
	leftnode=rightnode;
	rightnode=rightnode->next_sorted;
      }
      /*-------------------------
	Guard against case when node to
	be deleted is currentnode, otherwise
	the inserted node will point to 
	itself
	---------------------------*/
      if(rightnode==first_sequence){
	rightnode->data=nextsample;
	first_sequence=first_sequence->next_sequence;
	rightnode->next_sequence=NULL;
	last_sequence->next_sequence=rightnode;
	last_sequence=rightnode;
	shift=0;
      }
      else{ 
	if(leftnode==first_sequence){
	  leftnode->data=nextsample;
	  first_sequence=first_sequence->next_sequence;
	  leftnode->next_sequence=NULL;
	  last_sequence->next_sequence=leftnode;
	  last_sequence=leftnode;
	  shift=0; 
	}
	else {
	  reuse_next_sorted=rightnode;
	  reuse_prev_sorted=leftnode;
	  shift=1; /*shift maybe required*/
	}
      }
    }
    else{  
      chkcount=0;
      dummy_node=checks[0];
      if(dummy_node==first_sequence){
	dummy_node->data=nextsample;
	first_sequence=first_sequence->next_sequence;
	dummy_node->next_sequence=NULL;
	last_sequence->next_sequence=dummy_node;
	last_sequence=dummy_node;
	shift=0;
      }
      else{
	reuse_next_sorted=checks[0];
	reuse_prev_sorted=NULL;
	shift=1;  /*shift maybe required*/
      }
      rightnode=checks[0];
      leftnode=NULL;
    }

    /*-----------------------------------
      Getting check points to be shifted
      -----------------------------------*/
    if(shift){
      deletesample=first_sequence->data;
      if(deletesample>nextsample){
	shiftcounter=0;
	for(k=chkcount;k<ncheckpts;k++){
	  dummy_node=checks[k];
	  dummy=dummy_node->data;
	  if(dummy>=nextsample){
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
	    dummy_node=checks[k];
	    dummy=dummy_node->data;
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
	else{
	  /* nextsample=deletesample but they are separated in
	     the ordered list. This implies that all values in 
	     between are equal. So, just change the sequential list.
	  */
	  dummy_node=first_sequence;
	  dummy_node->data=nextsample;
	  last_sequence->next_sequence=dummy_node;
	  first_sequence=dummy_node->next_sequence;
	  dummy_node->next_sequence=0;
	  last_sequence=dummy_node;
	  shift=0; 
	}
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
			REAL4Sequence *input,
			LALRunningMedianPar param)
/* </lalVerbatim> */
{
  REAL4 *sorted_indices;

  struct rngmed_val_index4 *index_block;

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
  UINT4 samplecount,counter_chkpt,chkcount,shiftcounter = 0;
  INT8 k;
  REAL4 nextsample,deletesample,dummy;
  UINT4 shift,dummy_int;

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
    ABORT(status,LALRUNNINGMEDIANH_EMALOC2,LALRUNNINGMEDIANH_MSGEMALOC2);
  }
  checks4shift = (INT4*)LALCalloc(ncheckpts,sizeof(INT4));
  if(!checks4shift){
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
    LALFree(checks4shift);
    LALFree(checks);
    ABORT(status,LALRUNNINGMEDIANH_EMALOC4,LALRUNNINGMEDIANH_MSGEMALOC4);
  }
  first_sequence=(struct node *)LALCalloc(1,sizeof(struct node));
  if(!first_sequence){
    LALFree(node_addresses);
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
    if(nextsample>checks[0]->data){
      for(chkcount=1;chkcount<ncheckpts;chkcount++){
	if(nextsample>checks[chkcount]->data){
	}
	else{
	  break;
	}
      }
      chkcount-=1;
      rightnode=checks[chkcount];
      while(rightnode){
	if(nextsample<=rightnode->data){
	  break;
	}
	leftnode=rightnode;
	rightnode=rightnode->next_sorted;
      }
      /*-------------------------
	Guard against case when node to
	be deleted is currentnode, otherwise
	the inserted node will point to 
	itself
	---------------------------*/
      if(rightnode==first_sequence){
	rightnode->data=nextsample;
	first_sequence=first_sequence->next_sequence;
	rightnode->next_sequence=NULL;
	last_sequence->next_sequence=rightnode;
	last_sequence=rightnode;
	shift=0;
      }
      else{ 
	if(leftnode==first_sequence){
	  leftnode->data=nextsample;
	  first_sequence=first_sequence->next_sequence;
	  leftnode->next_sequence=NULL;
	  last_sequence->next_sequence=leftnode;
	  last_sequence=leftnode;
	  shift=0; 
	}
	else {
	  reuse_next_sorted=rightnode;
	  reuse_prev_sorted=leftnode;
	  shift=1; /*shift maybe required*/
	}
      }
    }
    else{  
      chkcount=0;
      dummy_node=checks[0];
      if(dummy_node==first_sequence){
	dummy_node->data=nextsample;
	first_sequence=first_sequence->next_sequence;
	dummy_node->next_sequence=NULL;
	last_sequence->next_sequence=dummy_node;
	last_sequence=dummy_node;
	shift=0;
      }
      else{
	reuse_next_sorted=checks[0];
	reuse_prev_sorted=NULL;
	shift=1;  /*shift maybe required*/
      }
      rightnode=checks[0];
      leftnode=NULL;
    }

    /*-----------------------------------
      Getting check points to be shifted
      -----------------------------------*/
    if(shift){
      deletesample=first_sequence->data;
      if(deletesample>nextsample){
	shiftcounter=0;
	for(k=chkcount;k<ncheckpts;k++){
	  dummy_node=checks[k];
	  dummy=dummy_node->data;
	  if(dummy>=nextsample){
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
	    dummy_node=checks[k];
	    dummy=dummy_node->data;
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
	else{
	  /* nextsample=deletesample but they are separated in
	     the ordered list. This implies that all values in 
	     between are equal. So, just change the sequential list.
	  */
	  dummy_node=first_sequence;
	  dummy_node->data=nextsample;
	  last_sequence->next_sequence=dummy_node;
	  first_sequence=dummy_node->next_sequence;
	  dummy_node->next_sequence=0;
	  last_sequence=dummy_node;
	  shift=0; 
	}
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


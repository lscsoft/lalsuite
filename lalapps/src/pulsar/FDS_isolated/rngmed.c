/** \file rngmed.c
Implementation file for running median.
\author Soumya D. Mohanty
\date May 2001
*/


#include <stdlib.h>
#include <math.h>

#include <lal/LALStdlib.h>

#include "rngmed.h"

/** 
Computes running median in an efficient manner.

      \param data Pointer to input data array

      \param lendata Length of input data array

      \param nblocks Block size for computing running median
    
      \param medians Pointer to output array. Number of elements is 
                      lendata - nblocks+1. Must be 
                       allocated outside this function.

*/
void rngmed(const double *data, unsigned int lendata, unsigned int nblocks, double *medians){


/*--------------------------
Structure to store samples in one block
----------------------------*/
     struct node{
           double data;
           struct node *next_sorted; 
           struct node *next_sequence;
           struct node *prev_sorted;
     };


/**\page page1 Variables
\code
     double *sorted_indices;
     struct rngmed_val_index *index_block;
     struct node **checks,**node_addresses;    
     struct node *first_sequence,*last_sequence;
     struct node *currentnode,*previousnode; 
     struct node *leftnode, *rightnode;
     struct node *reuse_next_sorted,*reuse_prev_sorted;
     struct node *dummy_node,*dummy_node1,*dummy_node2;
     int ncheckpts,stepchkpts;
     int nextchkptindx,*checks4shift;
     int nearestchk,midpoint,offset,numberoffsets;
     int samplecount,k,counter_chkpt,chkcount,shiftcounter;
     double nextsample,deletesample,dummy;
     int shift,dummy_int;
\endcode
*/

/*-------------------------------
checks: points to checkpoints 
---------------------------------*/
     double *sorted_indices;
     struct rngmed_val_index *index_block;
     struct node **checks,**node_addresses;    
     struct node *first_sequence,*last_sequence;
     struct node *currentnode,*previousnode; 
     struct node *leftnode, *rightnode;
     struct node *reuse_next_sorted,*reuse_prev_sorted;
     struct node *dummy_node,*dummy_node1,*dummy_node2;
     int ncheckpts,stepchkpts;
     int nextchkptindx,*checks4shift;
     int nearestchk,midpoint,offset,numberoffsets;
     int samplecount,k,counter_chkpt,chkcount,shiftcounter;
     double nextsample,deletesample,dummy;
     int shift,dummy_int;


/**\page page2 Sort the first block
 Allocate storage
for an array of \b rngmed_val_index. 
\code
     index_block =(struct rngmed_val_index *)LALCalloc(nblocks, sizeof(struct rngmed_val_index));
\endcode
Store the samples in the \b data member of each array element. Store the index of the
sample in the \b index member of each array element.
\code
     for(k=0;k<nblocks;k++){
          index_block[k].data=data[k];
          index_block[k].index=k;
     }
\endcode
Pass the array to \e qsort along with pointer to function \e rngmed_sortindex.
\code
     qsort(index_block, nblocks, sizeof(struct rngmed_val_index),rngmed_sortindex);
\endcode
Get the original indices of the samples in the sorted list. This list of indices
is used at the start of the core part of the code.
\code
     sorted_indices=(double *)LALCalloc(nblocks,sizeof(double));
     for(k=0;k<nblocks;k++){
         sorted_indices[k]=index_block[k].index; 
     }
\endcode
*/     

/*-----------------------------------
  Sort the first block of nblocks samples
  using the MATLAB sort function
------------------------------------*/
     index_block =(struct rngmed_val_index *)LALCalloc(nblocks, sizeof(struct rngmed_val_index));
     if(!index_block){
          printf("Could not allocate memory for index_block\n");
          return;
     }
     for(k=0;k<nblocks;k++){
          index_block[k].data=data[k];
          index_block[k].index=k;
     }

     qsort(index_block, nblocks, sizeof(struct rngmed_val_index),rngmed_sortindex);

     sorted_indices=(double *)LALCalloc(nblocks,sizeof(double));
     if (!sorted_indices){
           printf("Could not allocate memory for sorted_indices\n");
           return;
     }
     for(k=0;k<nblocks;k++){
         sorted_indices[k]=index_block[k].index; 
     }

     LALFree(index_block);

/** \page page_setupchecks Set up checkpoints
\e Checkpoints are special nodes in the linked list
containing the block of \b nblocks samples. A new sample
is first compared against the values stored in these
special nodes.

Get the number of nodes between consecutive checkpoints.
\code
     stepchkpts=sqrt(nblocks);
\endcode
Get the number of checkpoints to use.
\code
     ncheckpts=nblocks/stepchkpts;
\endcode
Allocate array to hold the pointers to checkpoints.
\code
     checks=(struct node **)LALCalloc(ncheckpts,sizeof(struct node*));
\endcode
The insertion of a new sample and the deletion of an old sample
from the linked list containing the block of \b nblocks samples leads
to the shifting of the checkpoints. The indices of elements in \b checks
that need to be changed is stored in \b checks4shift.
\code
     if(!(checks4shift=(int*)LALCalloc(ncheckpts,sizeof(int)))){
           printf("Could not allocate storage for checks4shift\n");
           return;
     }
\endcode
*/
/*----------------------------------
Indices of checkpoint nodes.
Number of nodes per checkpoint=floor(sqrt(nblocks))
------------------------------------*/
     stepchkpts=sqrt(nblocks);
     ncheckpts=nblocks/stepchkpts;
     checks=(struct node **)LALCalloc(ncheckpts,sizeof(struct node*));
     if(!checks){
           printf("Could not allocate storage for checks\n");
           return;
     }
     if(!(checks4shift=(int*)LALCalloc(ncheckpts,sizeof(int)))){
           printf("Could not allocate storage for checks4shift\n");
           return;
     }
/** \page page_offset Get the nearest checkpoint 
Get the nearest checkpoint to the node(s) containing the median
of the block of \b nblocks samples. For \e odd \b nblocks, the 
median is the datum in the node that lies at the
\b midpoint of the sorted list. For \e even \b nblocks, the median
is the average of the data in the two nodes at the middle. The two cases are
distinguished by the flag \b numberoffsets. 
\code
     if((int)fmod(nblocks,2.0)){
           midpoint=(nblocks+1)/2-1;
           numberoffsets=1;
     }
     else{
           midpoint=nblocks/2-1;
           numberoffsets=2;   
     }
\endcode
Get the nearest checkpoint to the median.
This is used for fast access to the node(s) containing the 
median value. The usual method of access in an elementary
 linked list is
\e sequential, which is slower.
\code
     nearestchk=floor(midpoint/stepchkpts);
     offset=midpoint-nearestchk*stepchkpts;
\endcode
*/

/*---------------------------------
  Offsets for getting median from nearest
  checkpoint: For nblocks even, 
  (node(offset(1))+node(offset(2)))/2;
  for nblocks odd,
  (node(offset(1))+node(offset(1)))/2;
 ----------------------------------*/
     if((int)fmod(nblocks,2.0)){
    /*Odd*/
           midpoint=(nblocks+1)/2-1;
           numberoffsets=1;
     }
     else{
    /*Even*/
           midpoint=nblocks/2-1;
           numberoffsets=2;   
     }
     nearestchk=floor(midpoint/stepchkpts);
     offset=midpoint-nearestchk*stepchkpts;

/** \page page_setuplist Set up the linked list
The linked list containing \b nblocks samples is set up and initialized.
This list has three types of links from one node to another. See the
documentation for \b struct \b node. First, the sequential ordering is
between nodes is set up. But the addresses of the nodes
are needed to set up later the links representing the sorted order. 
The linked list has 
bidirectional links going in both the \e ascending and \e descending order.

Allocate storage to
hold the node addresses.
\code
     node_addresses=(struct node **)LALCalloc(nblocks,sizeof(struct node *));
\endcode
Create a node. This stores the sequentially first sample in the 
block of \b nblocks samples.
\code
     first_sequence=(struct node *)LALCalloc(1,sizeof(struct node));
\endcode
Store its address.
\code
     node_addresses[0]=first_sequence;
\endcode
Initialize this node. 
\code
     first_sequence->next_sequence=NULL;
     first_sequence->next_sorted=NULL;
     first_sequence->prev_sorted=NULL;
     first_sequence->data=data[0];
\endcode
Start loop to setup links between nodes and load the sample values.
Only the links representing the \e sequential order are setup. Links representing
the sorted order are setup later.

BEGIN FOR LOOP.
\code
     previousnode=first_sequence;
     for(samplecount=1;samplecount<nblocks;samplecount++){
           currentnode=(struct node *)LALCalloc(1,sizeof(struct node));
           if(!currentnode){
                printf("Could not create node ");
                return;
           }
\endcode
Store the address of the node.
\code
           node_addresses[samplecount]=currentnode;
\endcode
Link from already allocated node \b previousnode to the new one, \b currentnode.
Load current data sample into \b currentnode.
\code
           previousnode->next_sequence=currentnode;
           currentnode->next_sequence=NULL;
           currentnode->prev_sorted=NULL;
           currentnode->next_sorted=NULL;
           currentnode->data=data[samplecount];
\endcode
Set \b currentnode as \b previousnode for next iteration of the loop.
\code
           previousnode=currentnode;
     }
\endcode
END FOR LOOP.

This stores the sequentially last sample of the block of \b nblocks samples.
\code
     last_sequence=currentnode;
\endcode
*/
/*----------------------------------
Build up linked list using first nblock points
in sequential order
------------------------------------*/
     node_addresses=(struct node **)LALCalloc(nblocks,sizeof(struct node *));
     if(!node_addresses){
           printf("Could not allocate storage for node_addresses\n");
           return;
     }
     first_sequence=(struct node *)LALCalloc(1,sizeof(struct node));
     if(!first_sequence){
           printf("Could not create first node\n");
           return;
     }
     node_addresses[0]=first_sequence;
     first_sequence->next_sequence=NULL;
     first_sequence->next_sorted=NULL;
     first_sequence->prev_sorted=NULL;
     first_sequence->data=data[0];
     previousnode=first_sequence;
     for(samplecount=1;samplecount<nblocks;samplecount++){
           currentnode=(struct node *)LALCalloc(1,sizeof(struct node));
           if(!currentnode){
                printf("Could not create node ");
                return;
           }
           node_addresses[samplecount]=currentnode;
           previousnode->next_sequence=currentnode;
           currentnode->next_sequence=NULL;
           currentnode->prev_sorted=NULL;
           currentnode->next_sorted=NULL;
           currentnode->data=data[samplecount];
           previousnode=currentnode;
     }
     last_sequence=currentnode;

/** \page page_setuplist_sorted Set up sorted list
The links between nodes representing sorted order are
set up.
For this the previously stored indices, after sorting using
\e qsort, are used.

Node containing the lowest value.
\code
     currentnode=node_addresses[(int)sorted_indices[0]];
\endcode
This is also the first checkpoint.
\code
     checks[0]=currentnode;
\endcode
BEGIN FOR LOOP.
\code
     previousnode=NULL;
     nextchkptindx=stepchkpts;
     counter_chkpt=1;
     for(samplecount=1;samplecount<nblocks;samplecount++){
\endcode
Get the address of the node containing the next highest value.
\code
          dummy_node=node_addresses[(int)sorted_indices[samplecount]];
\endcode
Make a link from \b current_node to this address. A second link is
made to the previous lower value also. Thus, the linked list has 
bidirectional links going in both the ascending and descending order.
\code
          currentnode->next_sorted=dummy_node;
          currentnode->prev_sorted=previousnode;
          previousnode=currentnode;
          currentnode=dummy_node;
\endcode
If the node is also a checkpoint then record its address in the \b checks array.
\code
          if(samplecount==nextchkptindx && counter_chkpt<ncheckpts){
                checks[counter_chkpt]=currentnode;
                nextchkptindx+=stepchkpts;
                counter_chkpt++;
          }
     }
\endcode
END FOR LOOP.

Set up links for the last node in sorted order.
\code
     currentnode->prev_sorted=previousnode;
     currentnode->next_sorted=NULL;
\endcode

*/
/*------------------------------------
Set the sorted sequence pointers and
the pointers to checkpoint nodes
-------------------------------------*/
     currentnode=node_addresses[(int)sorted_indices[0]];
     previousnode=NULL;
     checks[0]=currentnode;
     nextchkptindx=stepchkpts;
     counter_chkpt=1;
     for(samplecount=1;samplecount<nblocks;samplecount++){
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


/** \page page_core The core part
This is the core engine of the code.

The output is stored in the array \b medians which should be
allocated outside the code.

First, get the median of the first block of samples.
Go to the nearest checkpoint.
\code
     currentnode=checks[nearestchk];
\endcode
Follow link from this node to the node containing the median. 
\code
     for(k=1;k<=offset;k++){
           currentnode=currentnode->next_sorted;
     }
\endcode
Depending on odd or even \b nblocks, get the value from the \b data
member of \b currentnode or calculate the average of this value and
the value stored in the nest node.
\code
     dummy=0;
     for(k=1;k<=numberoffsets;k++){
           dummy+=currentnode->data;
           currentnode=currentnode->next_sorted;
     }
     medians[0]=dummy/numberoffsets;
\endcode
Now move to calculation of medians for successive blocks.

BEGIN FOR LOOP.
\code
     for(samplecount=nblocks;samplecount<lendata;samplecount++){
          nextsample=data[samplecount];
\endcode

\section sec_main_part Locate point of insertion

First the point where the nextsample must be inserted in the
linked list constructed above must be located.

Compare the new sample, \b nextsample, with checkpoints. 
There are
two cases to be considered. 
\subsection case1 Case 1
\code
          if(nextsample>checks[0]->data){
\endcode
 Find a checkpoint that is greater than 
the new sample. 
\code
                  for(chkcount=1;chkcount<ncheckpts;chkcount++){
                          if(nextsample>checks[chkcount]->data){
                          }
                          else{
                               break;
                          }
                  }
\endcode
Back up to previous checkpoint.
\code
                  chkcount-=1;
\endcode
\b rightnode is the node that
lies immediately to the right of the new sample 
in ascending order 
and \b leftnode
is the node immediately on the left.

Follow the link in ascending order, starting from the checkpoint to
the left of \b nextsample, until the bracketing nodes are found.
\code
                  rightnode=checks[chkcount];
                  while(rightnode){
                          if(nextsample<=rightnode->data){
                                break;
                          }
                          leftnode=rightnode;
                          rightnode=rightnode->next_sorted;
                  }
\endcode
The new sample must be inserted as a node between \b leftnode and \b rightnode.

The node containing the 
sequentially first sample of the block, \b first_sequence, must be removed
from the list and the same node is reused to store the new sample. This node
then becomes the sequentially last, i.e., \b last_sequence.

Special care is needed if the node to be removed
also happens to be either \b rightnode or \b leftnode. In this case, the checkpoints
need not be shifted (\b shift = 0).  
\code
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
\endcode
Otherwise the checkpoints may need to be shifted (\b shift = 1). 
\code
                          else {
                                 reuse_next_sorted=rightnode;
                                 reuse_prev_sorted=leftnode;
                                 shift=1; 
                          }
                  }
\endcode
The nodes to the right and left of the node that will be recycled and inserted
are \b rightnode and \b leftnode. So, store the addresses of these nodes in
\b reuse_next(prev)_sorted. The links between \b rightnode and \b leftnode 
must be broken and pointed to the new node. Conversely the new node must form 
links to \b rightnode and \b leftnode.

\subsection case2 Case 2
This is the case where \b nextsample \f$\leq\f$ \b checks[0]->data. 
Recall that the first checkpoint \b checks[0] also holds the lowest 
value in sorted order. But this need not be the sequentially first value.
Hence, distinguish the two cases below.
\code
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
                        shift=1;  
                  }
                  rightnode=checks[0];
                  leftnode=NULL;
          }
\endcode
\section sec3 Find checkpoints to shift

If the sequentially first sample is not immediately 
adjacent to \b nextsample then deleting the node
containing the first sample and reinserting the node
elsewhere in the list (as located in the above code) 
requires that the intermediate checkpoints be shifted
by one node. 
\code
          if(shift){
                  deletesample=first_sequence->data;

\endcode
The direction of the shift depends on whether
the sequentially first sample, \b deletesample,
is less than, greater than or equal to \b nextsample. 
\subsection case_shift_1 Sequentially first greater than sequentially last 
The checkpoints must be shifted to point to the nodes with
immediately lower values in the ordered list.
\code
                  if(deletesample>nextsample){
                       shiftcounter=0;
                       for(k=chkcount;k<ncheckpts;k++){
                             dummy_node=checks[k];
                             dummy=dummy_node->data;
                             if(dummy>=nextsample){
                                   if(dummy<=deletesample){
\endcode
This checkpoint falls between \b deletesample and \b nextsample. So, it 
must be shifted. Store its index in \b checks4shift.
\code
                                         checks4shift[shiftcounter]=k;
                                         shiftcounter++;
                                   }
                                   else{
                                         break;
                                   }
                             }
                        }
                        shift=-1; 
                  }
\endcode
\subsection case_shift_2 Sequentially first less than sequentially last
The checkpoints must be shifted to point to the nodes with
immediately higher values in the ordered list.
\code
                  else 
                        if(deletesample<nextsample){
                              shiftcounter=0;
                              for(k=chkcount;k>=0;k--){
                                    dummy_node=checks[k];
                                    dummy=dummy_node->data;
\endcode
This checkpoint falls between \b deletesample and \b nextsample. So, it 
must be shifted. Store its index in \b checks4shift.
\code
                                    if(dummy>=deletesample){
                                         checks4shift[shiftcounter]=k;
                                         shiftcounter++;
                                    }
                                    else{
                                         break;
                                    }
                              }
                              shift=1; 
                       }
\endcode
\subsection case_shift_3 Sequentially first equal to sequentially last

The sequentially first and last samples are equal
 but they are separated in
the ordered list. This implies that all values in 
between in the ordered list must be
 equal. So, only the sequential list links need be changed.
\code
                       else{
                              dummy_node=first_sequence;
                              dummy_node->data=nextsample;
                              last_sequence->next_sequence=dummy_node;
                              first_sequence=dummy_node->next_sequence;
                              dummy_node->next_sequence=0;
                              last_sequence=dummy_node;
                              shift=0; 
                       }
          }
\endcode

\section sec_implement Implementing the link changes

Now the node containing the sequentially first sample, \b first_sequence,
must be recycled. This means its old links (both 
sequential and ordered) must be severed. The nodes
linking into \b first_sequence must be relinked. The links between 
the nodes immediately adjacent to the insertion point must be broken
and relinked to the recycled \b first_sequence node.
\code
         if(shift){
\endcode
Some cases did not require any shift and in those cases \b first_sequence 
has already been recycled. Do the following only if \b shift \f$\neq 0\f$.

First reset the sequential links and load \b nextsample into \b first_sequence.
The \b first_sequence node becomes \b dummy_node.
\code
                  dummy_node=first_sequence;
                  first_sequence=dummy_node->next_sequence;
                  dummy_node->next_sequence=NULL;
                  last_sequence->next_sequence=dummy_node;
                  last_sequence=dummy_node;
                  dummy_node->data=nextsample;
\endcode
Store the sorted order links.
\code
                  dummy_node1=dummy_node->prev_sorted;
                  dummy_node2=dummy_node->next_sorted;
\endcode

\subsection sec_implement_1 Repair the deletion point
The nodes adjacent to \b first_sequence in the sorted list are relinked.
Two cases must be considered. 

\paragraph sec_implement_1_1 Case 1
If \b first_sequence was also the first
checkpoint then its \b prev_sorted link (=\b dummy_node1 above) must be NULL.
This may be a redundant check since this case has been addressed earlier in the code
 (\b shift = 0).
But there is no harm in repeating it.
\code
                 if(!dummy_node1){
                        dummy_node2->prev_sorted=dummy_node1;
                 }
\endcode
\paragraph sec_implement_1_2 Case 2
If \b first_sequence was the last node in the ascending order list then
\b next_sorted must be NULL.
\code
                 else {
                        if(!dummy_node2){
                               dummy_node1->next_sorted=dummy_node2;
                        }
\endcode
\paragraph sec_implement_1_3 Case 3
The normal case where there are nodes on either side.
\code
                        else{
                               dummy_node1->next_sorted=dummy_node2;
                               dummy_node2->prev_sorted=dummy_node1;
                        }
                 } 
\endcode

\subsection sec_implement_2 Relink nodes at insertion point 
 
Recall that the adjacent nodes were \b rightnode and \b leftnode.
Link them into the recycled node.

If the insertion point is at the end of the 
ascending order sorted list, \b rightnode must
be NULL.
\code
               if(!rightnode){
                       leftnode->next_sorted=dummy_node;
                }
\endcode
If the insertion point is at the beginning of the ascending order
sorted list, \b leftnode must be NULL.
\code
                else {
                       if(!leftnode){
                               rightnode->prev_sorted=dummy_node;
                       }
\endcode
Else both \b leftnode and \b rightnode exist.
\code
                       else{
                               leftnode->next_sorted=dummy_node;
                               rightnode->prev_sorted=dummy_node;
                       }
                }
\endcode

\section sec_implement_shifts Implementing checkpoint shifts
Two cases here: shift to left or right.
\code
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
\endcode
Link the recycled node to \b leftnode and \b rightnode.
\code
              dummy_node->next_sorted=reuse_next_sorted;
              dummy_node->prev_sorted=reuse_prev_sorted;
         }
\endcode

\section sec_median Get the median
Go to the nearest check point.
\code
         currentnode=checks[nearestchk];
\endcode
Follow link through \b offset number of nodes.
\code
         for(k=1;k<=offset;k++){
              currentnode=currentnode->next_sorted;
         }
         dummy=0;
\endcode
Even and odd \b nblocks cases.
\code
         for(k=1;k<=numberoffsets;k++){
              dummy+=currentnode->data;
              currentnode=currentnode->next_sorted;
         }
         medians[samplecount-nblocks+1]=dummy/numberoffsets;
\endcode
END FOR LOOP.
\code
     }
\endcode
*/

/*------------------------------
  Get the first output element
-------------------------------*/
     if(!medians){
           printf("Null output array");
           return;
     }
     currentnode=checks[nearestchk];
     for(k=1;k<=offset;k++){
           currentnode=currentnode->next_sorted;
     }
     dummy=0;
     for(k=1;k<=numberoffsets;k++){
           dummy+=currentnode->data;
           currentnode=currentnode->next_sorted;
     }
     medians[0]=dummy/numberoffsets;

/*---------------------------------
This is the main part
----------------------------------*/
     for(samplecount=nblocks;samplecount<lendata;samplecount++){
          nextsample=data[samplecount];
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
         medians[samplecount-nblocks+1]=dummy/numberoffsets;
     }/*Outer For Loop*/

/** \page page_final Clean up
Deallocate memory that was dynamically allocated. 
Not shown in this document are some sections of the code
where additional cleaning up is done during the course of 
execution. 
\code
     LALFree(node_addresses);
     currentnode=first_sequence;
     while(currentnode){
             previousnode=currentnode;
             currentnode=currentnode->next_sequence;
             LALFree(previousnode);
     }
     LALFree(checks);
\endcode
*/
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
     LALFree(checks);
     LALFree(checks4shift);
}

/*--------------------------
This function is passed to qsort
----------------------------*/
/**  This function is passed to \em qsort.
\param elem1 element of a rngmed_val_index array
\param elem2 another element of same rngmed_val_index array 
*/
int rngmed_sortindex(const void *elem1, const void *elem2){
    /*Used in running qsort*/
    
    struct rngmed_val_index *A, *B;
    double data1, data2;

    A=(struct rngmed_val_index *)elem1;
    B=(struct rngmed_val_index *)elem2;

    data1=A->data;
    data2=B->data;
    if (data1 < data2)
         return -1;
    else if (data1==data2)
         return 0;
    else
         return 1;

}


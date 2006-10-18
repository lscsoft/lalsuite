/** \file rngmed.c
Implementation file for running median.
\author Soumya D. Mohanty
\date May 2001
$Id$


Revision 1.2  2004/07/26 14:36:14  bema
New version by Soumya Mohanty to fix the "same values" bug

Revision 1.6  2004/07/23 16:28:27  mohanty
Removed print commands used in debugging. Changed some comments.

Revision 1.5  2004/07/23 16:14:16  mohanty
*** empty log message ***

Revision 1.4  2004/07/23 16:07:22  mohanty
Major changes in response to bug: segmentation violation
when the running block comes upon a string of identical values that are
guaranteed to be the smallest values in the data. (In
this case, the first check point is guranteed to contain this smallest value.)
The bug was caused by the fact that the next identical value was placed in the
sorted list (ascending) *before* an old value.
This caused the check point to move asynchronously with the pointer
to the node containing the sequentially first sample. The segmentation
violation actually happened much later and would not be detected 
if the total data length were not long enough. This explains why the bug
was hard to catch.
Also took this opportunity to reduce over-engineering in the code:
(1) Combined the nextsample>checks[0]->data and nextsample==checks[0]->data
cases. Same code works for both except nextsample > checks[k]->data
needs to be replaced with nextsample >= checks[k]->data
(2) Combined the deletesample<nextsample and deletesample==nextsample
cases into one. Same code works for both.
NOTE: Need to check for this effect in the first block too because
qsort does not preserve sequence when sorting by value.

Revision 1.3  2004/07/22 17:25:08  mohanty

*/


#include <lal/LALStdlib.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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
int rngmed(const double *data, unsigned int lendata, unsigned int nblocks, double *medians){


/*--------------------------
Structure to store samples in one block
----------------------------*/
     struct node{
           double data;
           struct node *next_sorted; 
           struct node *next_sequence;
           struct node *prev_sorted;
     };



/*-------------------------------
checks: pointers to subset of nodes to use
as checkpoints 
---------------------------------*/
     double *sorted_indices;
     struct rngmed_val_index *index_block;
     struct node **checks,**node_addresses;    
     struct node *first_sequence,*last_sequence;
     struct node *currentnode=NULL,*previousnode; 
     struct node *leftnode=NULL, *rightnode=NULL;
     struct node *reuse_next_sorted=NULL,*reuse_prev_sorted=NULL;
     struct node *dummy_node,*dummy_node1,*dummy_node2;
     int ncheckpts,stepchkpts;
     int nextchkptindx,*checks4shift;
     int nearestchk,midpoint,offset,numberoffsets;
     int samplecount,k,counter_chkpt,chkcount=0,shiftcounter=0;
     double nextsample,deletesample,dummy;
     int shift,dummy_int;

/*-----------------------------------
  Sort the first block of nblocks samples
  using the qsort function
------------------------------------*/
     index_block =(struct rngmed_val_index *)LALCalloc(nblocks, sizeof(struct rngmed_val_index));
     if(!index_block){
          printf("Could not allocate memory for index_block\n");
          return 1;
     }
     for(k=0;k<(int)nblocks;k++){
          index_block[k].data=data[k];
          index_block[k].index=k;
     }

     qsort(index_block, nblocks, sizeof(struct rngmed_val_index),rngmed_sortindex);

     sorted_indices=(double *)LALCalloc(nblocks,sizeof(double));
     if (!sorted_indices){
           printf("Could not allocate memory for sorted_indices\n");
           return 1;
     }
     for(k=0;k<(int)nblocks;k++){
         sorted_indices[k]=index_block[k].index; 
     }

     LALFree(index_block);

/*----------------------------------
Indices of checkpoint nodes.
Number of nodes per checkpoint=floor(sqrt(nblocks))
------------------------------------*/
     stepchkpts=sqrt(nblocks);
     ncheckpts=nblocks/stepchkpts;
     checks=(struct node **)LALCalloc(ncheckpts,sizeof(struct node*));
     if(!checks){
           printf("Could not allocate storage for checks\n");
           return 1;
     }
     if(!(checks4shift=(int*)LALCalloc(ncheckpts,sizeof(int)))){
           printf("Could not allocate storage for checks4shift\n");
           return 1;
     }

/*---------------------------------
  Offsets for getting median from nearest
  checkpoint: For nblocks even, 
  (node(offset(1))+node(offset(2)))/2;
  for nblocks odd,
  (node(offset(1))+node(offset(1)))/2;
  THIS CAN BE OPTIMISED.
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


/*----------------------------------
Build up linked list using first nblock points
in sequential order
------------------------------------*/
     node_addresses=(struct node **)LALCalloc(nblocks,sizeof(struct node *));
     if(!node_addresses){
           printf("Could not allocate storage for node_addresses\n");
           return 1;
     }
     first_sequence=(struct node *)LALCalloc(1,sizeof(struct node));
     if(!first_sequence){
           printf("Could not create first node\n");
           return 1;
     }
     node_addresses[0]=first_sequence;
     first_sequence->next_sequence=NULL;
     first_sequence->next_sorted=NULL;
     first_sequence->prev_sorted=NULL;
     first_sequence->data=data[0];
     previousnode=first_sequence;
     for(samplecount=1;samplecount<(int)nblocks;samplecount++){
           currentnode=(struct node *)LALCalloc(1,sizeof(struct node));
           if(!currentnode){
                printf("Could not create node \n");
                return 1;
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

/*------------------------------------
Set the sorted sequence pointers and
the pointers to checkpoint nodes
-------------------------------------*/
     currentnode=node_addresses[(int)sorted_indices[0]];
     previousnode=NULL;
     checks[0]=currentnode;
     nextchkptindx=stepchkpts;
     counter_chkpt=1;
     for(samplecount=1;samplecount<(int)nblocks;samplecount++){
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
  Get the first output element
-------------------------------*/
     if(!medians){
           printf("Null output array");
           return 1;
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
This is the main part.
Find the nodes whose values
form the smallest closed interval
around the new incoming value.
The right limit is always >
the new value.
----------------------------------*/
     for(samplecount=nblocks;samplecount<(int)lendata;samplecount++){         
          nextsample=data[samplecount];
          if(nextsample>=checks[0]->data){
                  for(chkcount=1;chkcount<ncheckpts;chkcount++){
                          if(nextsample>=checks[chkcount]->data){
                          }
                          else{
                               break;
                          }
                  }
                  chkcount-=1;
                  rightnode=checks[chkcount];
                  leftnode=NULL; /*NEW*/
                  while(rightnode){
                          if(nextsample<rightnode->data){                                
                                break;
                          }
                          leftnode=rightnode;
                          rightnode=rightnode->next_sorted;
                  }
   
          }
          else {
                  if(nextsample<checks[0]->data){  
                      chkcount=0;
                      /* dummy_node=checks[0]; */
                      rightnode=checks[0];
                      leftnode=NULL;  
                  }
          }

/*-------------------------
     Determine if checkpoints need to be 
     shifted or not.
   ---------------------------*/
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

/*-----------------------------------
   Getting check points to be shifted
 -----------------------------------*/
          if(shift){
                  deletesample=first_sequence->data;
                  if(deletesample>nextsample){
                       shiftcounter=0;
                       for(k=chkcount;k<ncheckpts;k++){
                             dummy=checks[k]->data;
                             if(dummy>nextsample){
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
                        if(deletesample<=nextsample){
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
          }

 /*------------------------------
  Recycle the node with the 
  oldest value. 
 --------------------------------*/
          if(shift){
   /*---------------------
    Reset sequential links
    ---------------------*/
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

     return 0;
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
    
    const struct rngmed_val_index *A = (struct rngmed_val_index *)elem1;
    const struct rngmed_val_index *B = (struct rngmed_val_index *)elem2;

    double data1, data2;

    data1=A->data;
    data2=B->data;
    if (data1 < data2)
         return -1;
    else if (data1==data2)
         return 0;
    else
         return 1;

}


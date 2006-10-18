#ifndef RNGMED_H
#define RNGMED_H

int rngmed(const double *, unsigned int, unsigned int, double *);
int rngmed_sortindex(const void *, const void *); 

/*----------------------------
Structure for storing elements of
an array along with position
of the element in the array
---------------------------*/
    struct rngmed_val_index{
        double data;
        unsigned int index;
    };

#else
#endif

/**
\struct node
This structure is used to make a linked list. The list holds the samples
in one block in the running median algorithm.

         \param data Holds a single number.
         \param next_sorted Points to the next node in the sorted list.
         \param prev_sorted Points to the previous node in the sorted list.
         \param next_sequence point to the next node in the sequential list.
*/

/** 
\struct rngmed_val_index
A structure to store values and indices
of elements in an array

\param data Stores a single number 
\param index Stores the original position of the number

This structure is used to track the indices
of elements after sorting by qsort. An array of
rngmed_val_index is passed to qsort which 
rearranges the array according to the values in
the data member. The
indices of these elements in the original unsorted 
array can then be read off from the index member.
*/

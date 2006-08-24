#include <stdio.h>
#include "HeapToplist.h"

static int smaller(const void*a,const void*b) { 
  if(*((double*)a) > *((double*)b))
    return(1);
  else if(*((double*)a) > *((double*)b))
    return(-1);
  else
    return(0);
}

static void print_elem(const void*e){
  printf(" %d",*((double*)e));
}

main(int argc, char**argv){
  double elem;
  int i;
  toplist_t *l;
  create_toplist(&l,10,sizeof(double),smaller);
  for(i=0;i<10;i++){
    elem = i;
    insert_into_toplist(l,&elem);
  }
  go_through_toplist(l,print_elem);
  free_toplist(&l);
}

  

/*
*  Copyright (C) 2007 Bernd Machenschalk
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

/**
 * \file
 * \ingroup pulsarApps
 * \author Bernd Machenschalk
 */

#include <stdio.h>
#include <stdlib.h>
#include "HeapToplist.h"

#define elem_t float

static int smaller(const void*a,const void*b) { 
  if(*((elem_t*)a) < *((elem_t*)b))
    return(1);
  else if(*((elem_t*)a) > *((elem_t*)b))
    return(-1);
  else
    return(0);
}

static void print_elem(const void*e){
  printf("%f ",*((elem_t*)e));
}

main(int argc, char**argv){
  elem_t elem;
  int i,n,m;
  toplist_t *l;

  if (argc < 2)
    n=20;
  else
    n=atoi(argv[1]);

  if (argc < 3)
    m=10;
  else
    m=atoi(argv[2]);

  create_toplist(&l,m,sizeof(elem_t),smaller);
  for(i=0;i<n;i++){
    elem = rand()/(double)RAND_MAX;
    insert_into_toplist(l,&elem);
    go_through_toplist(l,print_elem); printf("\n");
  }
  qsort_toplist_r(l,smaller);
  go_through_toplist(l,print_elem); printf("\n");
  qsort_toplist(l,smaller);
  go_through_toplist(l,print_elem); printf("\n");
  free_toplist(&l);
}

  

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
 */

#include <stdlib.h>
#include "HeapToplist.h"
#include "FstatToplist.h"

static int fstat_smaller(const void*a, const void*b) {
  if     (((FstatOutputEntry*)a)->Fstat < ((FstatOutputEntry*)b)->Fstat)
    return(1);
  else if(((FstatOutputEntry*)a)->Fstat > ((FstatOutputEntry*)b)->Fstat)
    return(-1);
  else if(((FstatOutputEntry*)a)->Freq  < ((FstatOutputEntry*)b)->Freq)
    return(1);
  else if(((FstatOutputEntry*)a)->Freq  > ((FstatOutputEntry*)b)->Freq)
    return(-1);
  else
    return(0);
}

static void print_fstat(const void*e){
  fprintf(stderr,"%e %e %e %e %e\n",
	  ((FstatOutputEntry*)e)->Freq,
	  ((FstatOutputEntry*)e)->f1dot,
	  ((FstatOutputEntry*)e)->Alpha,
	  ((FstatOutputEntry*)e)->Delta,
	  ((FstatOutputEntry*)e)->Fstat);
}

/* short test for the basic functions of this lib. */
int test_fstat_toplist(UINT8 n, UINT8 m, char*filename) {
    REAL8 epsilon=1e-5;
    FILE*fp;
    toplist_t *tl=NULL, *tl2=NULL;
    FstatOutputEntry FstatLine;
    FstatOutputEntry *ll;
    UINT8 i,ins=0;
    UINT4 checksum;

    ll=(FstatOutputEntry*)malloc(m*sizeof(FstatOutputEntry));
    if(ll == NULL) {
	XLALPrintError("Couldn't create list ll\n");
	return(-1);
    }

    fprintf(stderr,"creating first toplist...\n");
    if(create_fstat_toplist(&tl,n)) {
	XLALPrintError("Couldn't create toplist tl\n");
	return(-1);
    }
    fprintf(stderr,"creating second toplist...\n");
    if(create_fstat_toplist(&tl2,n)) {
	XLALPrintError("Couldn't create toplist tl2\n");
	free_fstat_toplist(&tl);
	return(-1);
    }

    fprintf(stderr,"open file %s for writing...\n",filename);
    fp=fopen(filename,"w");
    if(!fp) {
	XLALPrintError("Couldn't open file %s for writing\n",filename);
	free_fstat_toplist(&tl);
	free_fstat_toplist(&tl2);
	return(-2);
    }

    checksum = 0;

    fprintf(stderr,"filling toplist...\n");
    for(i=0;i<m;i++) {
	FstatLine.Freq = (double)rand() / (double)RAND_MAX;
	FstatLine.f1dot = (double)rand() / (double)RAND_MAX;
	FstatLine.Alpha = (double)rand() / (double)RAND_MAX;
	FstatLine.Delta = (double)rand() / (double)RAND_MAX;
	FstatLine.Fstat = (double)rand() / (double)RAND_MAX;

	ll[i]= FstatLine;
	
	if(insert_into_fstat_toplist(tl, FstatLine)) {
	    ins++;
	    write_fstat_toplist_item_to_fp(FstatLine,fp,&checksum);
	}

	if(i==n)
	    fprintf(stderr,"array filled now\n");
    }
    fprintf(stderr,"%d inserts actually done of %d\n",(int)ins,m);

    fclose(fp);

    fprintf(stderr,"file checksum: %d\n",checksum);

    fprintf(stderr,"\nLet's see if we really kept the top elements:\n");

    fprintf(stderr,"sort the lists completely\n");
    qsort_toplist(tl,fstat_smaller);
    qsort(ll,m,sizeof(FstatOutputEntry),fstat_smaller);

    /*
      for(i=0;i<m;i++)
      print_fstat((void*)&(ll[i]));
      fprintf(stderr,"...\n");
      go_through_toplist(tl,print_fstat);
    */

    fprintf(stderr,"comparing...\n");
    for(i=0;i<n;i++)
      if((((FstatOutputEntry*)toplist_elem(tl,i))->Freq  -
	  ll[i].Freq > epsilon) ||
	 (((FstatOutputEntry*)toplist_elem(tl,i))->f1dot -
	  ll[i].f1dot > epsilon) ||
	 (((FstatOutputEntry*)toplist_elem(tl,i))->Alpha -
	  ll[i].Alpha > epsilon) ||
	 (((FstatOutputEntry*)toplist_elem(tl,i))->Delta -
	  ll[i].Delta > epsilon) ||
	 (((FstatOutputEntry*)toplist_elem(tl,i))->Fstat -
	  ll[i].Fstat > epsilon)) {
	
	XLALPrintError("line %d differs\n",i);
	    fprintf(stderr,"%e %e %e %e %e\n",
		    ((FstatOutputEntry*)toplist_elem(tl,i))->Freq,
		    ((FstatOutputEntry*)toplist_elem(tl,i))->f1dot,
		    ((FstatOutputEntry*)toplist_elem(tl,i))->Alpha,
		    ((FstatOutputEntry*)toplist_elem(tl,i))->Delta,
		    ((FstatOutputEntry*)toplist_elem(tl,i))->Fstat);
	    fprintf(stderr,"%e %e %e %e %e\n",
		    ll[i].Freq,
		    ll[i].f1dot,
		    ll[i].Alpha,
		    ll[i].Delta,
		    ll[i].Fstat);
	}

    fprintf(stderr,"\nNow see if we can rebuild the toplist from the written file:\n");

    fprintf(stderr,"open file %s for reading...\n",filename);
    fp=fopen(filename,"r");
    if(!fp) {
	XLALPrintError("Couldn't open file %s for reading\n",filename);
	free_fstat_toplist(&tl);
	free_fstat_toplist(&tl2);
	return(-2);
    }

    fprintf(stderr,"reading...\n");
    read_fstat_toplist_from_fp(tl2, fp, &checksum, 0);

    fclose(fp);

    fprintf(stderr,"checksum: %d\n",checksum);

    sort_fstat_toplist(tl);

    fprintf(stderr,"open file %s for writing...\n",filename);
    fp=fopen(filename,"w");
    if(!fp) {
	XLALPrintError("Couldn't open file %s for writing\n",filename);
	free_fstat_toplist(&tl);
	free_fstat_toplist(&tl2);
	return(-2);
    }
    
    fprintf(stderr,"writing...\n");
    if(write_fstat_toplist_to_fp(tl, fp, &checksum)<0) {
	XLALPrintError("Couldn't write toplist\n",filename);
	fclose(fp);
	free_fstat_toplist(&tl);
	free_fstat_toplist(&tl2);
	return(-2);
    }

    fprintf(stderr,"checksum: %d\n",checksum);

    fclose(fp);

    sort_fstat_toplist(tl2);

    fprintf(stderr,"comparing...\n");
    for(i=0;i<n;i++)
      if((((FstatOutputEntry*)toplist_elem(tl,i))->Freq  -
	  ((FstatOutputEntry*)toplist_elem(tl2,i))->Freq > epsilon) ||
	 (((FstatOutputEntry*)toplist_elem(tl,i))->f1dot -
	  ((FstatOutputEntry*)toplist_elem(tl2,i))->f1dot > epsilon) ||
	 (((FstatOutputEntry*)toplist_elem(tl,i))->Alpha -
	  ((FstatOutputEntry*)toplist_elem(tl2,i))->Alpha > epsilon) ||
	 (((FstatOutputEntry*)toplist_elem(tl,i))->Delta -
	  ((FstatOutputEntry*)toplist_elem(tl2,i))->Delta > epsilon) ||
	 (((FstatOutputEntry*)toplist_elem(tl,i))->Fstat -
	  ((FstatOutputEntry*)toplist_elem(tl2,i))->Fstat > epsilon)) {
	
	XLALPrintError("line %d differs\n",i);
	    fprintf(stderr,"%e %e %e %e %e\n",
		    ((FstatOutputEntry*)toplist_elem(tl,i))->Freq,
		    ((FstatOutputEntry*)toplist_elem(tl,i))->f1dot,
		    ((FstatOutputEntry*)toplist_elem(tl,i))->Alpha,
		    ((FstatOutputEntry*)toplist_elem(tl,i))->Delta,
		    ((FstatOutputEntry*)toplist_elem(tl,i))->Fstat);
	    fprintf(stderr,"%e %e %e %e %e\n",
		    ((FstatOutputEntry*)toplist_elem(tl2,i))->Freq,
		    ((FstatOutputEntry*)toplist_elem(tl2,i))->f1dot,
		    ((FstatOutputEntry*)toplist_elem(tl2,i))->Alpha,
		    ((FstatOutputEntry*)toplist_elem(tl2,i))->Delta,
		    ((FstatOutputEntry*)toplist_elem(tl2,i))->Fstat);
	}

    fprintf(stderr,"cleanup...\n");
    free_fstat_toplist(&tl);
    free_fstat_toplist(&tl2);
    free(ll);

    return(0);
}

int main(int argc,char**argv) {
    char*name="Toplist.test";
    int n=10000;
    int m=0;
    if(argc>1)
      n=atoi(argv[1]);
    if(argc>2)
      m=atoi(argv[2]);
    if(argc>3)
      name=argv[3];
    if(m<n)
      m=2*n;
    return(test_fstat_toplist(n,m,name));
}

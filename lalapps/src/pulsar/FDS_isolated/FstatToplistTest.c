#include <stdlib.h>
#include "HeapToplist.h"
#include "FstatToplist.h"

/* short test for the basic functions of this lib. */
int test_fstat_toplist(UINT8 n, char*filename) {
    REAL8 epsilon=1e-5;
    FILE*fp;
    toplist_t *tl=NULL, *tl2=NULL;
    FstatOutputEntry FstatLine;
    UINT8 i,ins=0;
    UINT4 checksum;

    fprintf(stderr,"creating first toplist...\n");
    if(create_fstat_toplist(&tl,n)) {
	LALPrintError("Couldn't create toplist tl\n");
	return(-1);
    }
    fprintf(stderr,"creating second toplist...\n");
    if(create_fstat_toplist(&tl2,n)) {
	LALPrintError("Couldn't create toplist tl2\n");
	free_fstat_toplist(&tl);
	return(-1);
    }

    fprintf(stderr,"open file %s for writing...\n",filename);
    fp=fopen(filename,"w");
    if(!fp) {
	LALPrintError("Couldn't open file %s for writing\n",filename);
	free_fstat_toplist(&tl);
	free_fstat_toplist(&tl2);
	return(-2);
    }

    checksum = 0;

    fprintf(stderr,"filling first toplist...\n");
    for(i=0;i<n*2;i++) {
	FstatLine.Freq = (double)rand() / (double)RAND_MAX;
	FstatLine.f1dot = (double)rand() / (double)RAND_MAX;
	FstatLine.Alpha = (double)rand() / (double)RAND_MAX;
	FstatLine.Delta = (double)rand() / (double)RAND_MAX;
	FstatLine.Fstat = (double)rand() / (double)RAND_MAX;

	if(insert_into_fstat_toplist(tl, FstatLine)) {
	    ins++;
	    write_fstat_toplist_item_to_fp(FstatLine,fp,&checksum);
	}

	if(i==n)
	    fprintf(stderr,"array filled now\n");
    }
    fprintf(stderr,"%d inserts actually done of %d\n",(int)ins,2*(int)n);

    fclose(fp);

    fprintf(stderr,"checksum: %d\n",checksum);

    fprintf(stderr,"open file %s for reading...\n",filename);
    fp=fopen(filename,"r");
    if(!fp) {
	LALPrintError("Couldn't open file %s for reading\n",filename);
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
	LALPrintError("Couldn't open file %s for writing\n",filename);
	free_fstat_toplist(&tl);
	free_fstat_toplist(&tl2);
	return(-2);
    }
    
    fprintf(stderr,"writing...\n");
    if(write_fstat_toplist_to_fp(tl, fp, &checksum)<0) {
	LALPrintError("Couldn't write toplist\n",filename);
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
      if((((FstatOutputEntry*)(tl->heap[i]))->Freq  -
	  ((FstatOutputEntry*)(tl2->heap[i]))->Freq > epsilon) ||
	 (((FstatOutputEntry*)(tl->heap[i]))->f1dot -
	  ((FstatOutputEntry*)(tl2->heap[i]))->f1dot > epsilon) ||
	 (((FstatOutputEntry*)(tl->heap[i]))->Alpha -
	  ((FstatOutputEntry*)(tl2->heap[i]))->Alpha > epsilon) ||
	 (((FstatOutputEntry*)(tl->heap[i]))->Delta -
	  ((FstatOutputEntry*)(tl2->heap[i]))->Delta > epsilon) ||
	 (((FstatOutputEntry*)(tl->heap[i]))->Fstat -
	  ((FstatOutputEntry*)(tl2->heap[i]))->Fstat > epsilon)) {
	
	LALPrintError("line %d differs\n",i);
	    fprintf(stderr,"%e %e %e %e %e\n",
		    ((FstatOutputEntry*)(tl->heap[i]))->Freq,
		    ((FstatOutputEntry*)(tl->heap[i]))->f1dot,
		    ((FstatOutputEntry*)(tl->heap[i]))->Alpha,
		    ((FstatOutputEntry*)(tl->heap[i]))->Delta,
		    ((FstatOutputEntry*)(tl->heap[i]))->Fstat);
	    fprintf(stderr,"%e %e %e %e %e\n",
		    ((FstatOutputEntry*)(tl2->heap[i]))->Freq,
		    ((FstatOutputEntry*)(tl2->heap[i]))->f1dot,
		    ((FstatOutputEntry*)(tl2->heap[i]))->Alpha,
		    ((FstatOutputEntry*)(tl2->heap[i]))->Delta,
		    ((FstatOutputEntry*)(tl2->heap[i]))->Fstat);
	}

    fprintf(stderr,"cleanup...\n");
    free_fstat_toplist(&tl);
    free_fstat_toplist(&tl2);

    return(0);
}

int main(int argc,char**argv) {
    char*name="Toplist.test";
    int n=10000;
    if(argc>1)
	n=atoi(argv[1]);
    if(argc>2)
	name=argv[2];
    return(test_fstat_toplist(n,name));
}

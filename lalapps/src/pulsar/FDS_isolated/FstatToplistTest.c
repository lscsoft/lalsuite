#include <stdlib.h>
#include "FstatToplist.h"

/* short test for the basic functions of this lib. */
int test_toplist(UINT8 n, char*filename) {
    REAL8 epsilon=1e-5;
    FILE*fp;
    toplist *tl=NULL, *tl2=NULL;
    FstatsClusterOutput FstatLine;
    UINT8 i,ins=0;

    fprintf(stderr,"creating first toplist...\n");
    if(create_toplist(&tl,n)) {
	LALPrintError("Couldn't create toplist tl\n");
	return(-1);
    }
    fprintf(stderr,"creating second toplist...\n");
    if(create_toplist(&tl2,n)) {
	LALPrintError("Couldn't create toplist tl2\n");
	free_toplist(&tl);
	return(-1);
    }

    fprintf(stderr,"filling first toplist...\n");
    for(i=0;i<2*n;i++) {
	FstatLine.Freq = (double)rand() / (double)RAND_MAX;
	FstatLine.f1dot = (double)rand() / (double)RAND_MAX;
	FstatLine.Alpha = (double)rand() / (double)RAND_MAX;
	FstatLine.Delta = (double)rand() / (double)RAND_MAX;
	FstatLine.Nbins = rand();
	FstatLine.mean = (double)rand() / (double)RAND_MAX;
	FstatLine.std = (double)rand() / (double)RAND_MAX;
	FstatLine.max = (double)rand() / (double)RAND_MAX;

	if(insert_into_toplist(tl, FstatLine))
	    ins++;
    }
    fprintf(stderr,"%d inserts actually done of %d\n",(int)ins,2*(int)n);
    
    fprintf(stderr,"open file %s for writing...\n",filename);
    fp=fopen(filename,"w");
    if(!fp) {
	LALPrintError("Couldn't open file %s for writing\n",filename);
	free_toplist(&tl);
	free_toplist(&tl2);
	return(-2);
    }

    fprintf(stderr,"writing...\n");
    if(write_toplist_to_fp(tl,fp)<0) {
	LALPrintError("Couldn't write toplist\n",filename);
	fclose(fp);
	free_toplist(&tl);
	free_toplist(&tl2);
	return(-2);
    }

    fclose(fp);

    fprintf(stderr,"open file %s for reading...\n",filename);
    fp=fopen(filename,"r");
    if(!fp) {
	LALPrintError("Couldn't open file %s for reading\n",filename);
	free_toplist(&tl);
	free_toplist(&tl2);
	return(-2);
    }

    fprintf(stderr,"reading...\n");
    read_toplist_from_fp(tl2,fp);

    fprintf(stderr,"comparing...\n");
    for(i=0;i<n;i++)
	if((tl->data[i].Freq  - tl2->data[i].Freq > epsilon) ||
	   (tl->data[i].f1dot - tl2->data[i].f1dot > epsilon) ||
	   (tl->data[i].Alpha - tl2->data[i].Alpha > epsilon) ||
	   (tl->data[i].Delta - tl2->data[i].Delta > epsilon) ||
	   (tl->data[i].Nbins - tl2->data[i].Nbins > epsilon) ||
	   (tl->data[i].mean  - tl2->data[i].mean > epsilon) ||
	   (tl->data[i].std   - tl2->data[i].std > epsilon) ||
	   (tl->data[i].max   - tl2->data[i].max > epsilon))
	    LALPrintError("line %d differs\n",i);

    fprintf(stderr,"cleanup...\n");
    free_toplist(&tl);
    free_toplist(&tl2);

    return(0);
}

int main(int argc,char**argv) {
    char*name="Toplist.test";
    int n=10000;
    if(argc>1)
	n=atoi(argv[1]);
    if(argc>2)
	name=argv[2];
    return(test_toplist(n,name));
}

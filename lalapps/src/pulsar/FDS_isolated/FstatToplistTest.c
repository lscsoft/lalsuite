#include <stdlib.h>
#include "FstatToplist.h"

/* short test for the basic functions of this lib. */
int test_toplist(UINT8 n, char*filename) {
    REAL8 epsilon=1e-5;
    FILE*fp;
    toplist_t *tl=NULL, *tl2=NULL;
    FstatsClusterOutput FstatLine;
    UINT8 i,ins=0;
    UINT4 checksum;

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

    fprintf(stderr,"open file %s for writing...\n",filename);
    fp=fopen(filename,"w");
    if(!fp) {
	LALPrintError("Couldn't open file %s for writing\n",filename);
	free_toplist(&tl);
	free_toplist(&tl2);
	return(-2);
    }

    checksum = 0;

    fprintf(stderr,"filling first toplist...\n");
    for(i=0;i<n*2;i++) {
	FstatLine.Freq = (double)rand() / (double)RAND_MAX;
	FstatLine.f1dot = (double)rand() / (double)RAND_MAX;
	FstatLine.Alpha = (double)rand() / (double)RAND_MAX;
	FstatLine.Delta = (double)rand() / (double)RAND_MAX;
	FstatLine.Nbins = rand();
	FstatLine.mean = (double)rand() / (double)RAND_MAX;
	FstatLine.std = (double)rand() / (double)RAND_MAX;
	FstatLine.max = (double)rand() / (double)RAND_MAX;

	if(insert_into_toplist(tl, FstatLine)) {
	    ins++;
	    write_toplist_item_to_fp(FstatLine,fp,&checksum);
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
	free_toplist(&tl);
	free_toplist(&tl2);
	return(-2);
    }

    fprintf(stderr,"reading...\n");
    read_toplist_from_fp(tl2, fp, &checksum, 0);

    fclose(fp);

    fprintf(stderr,"checksum: %d\n",checksum);

    sort_toplist(tl);

    fprintf(stderr,"open file %s for writing...\n",filename);
    fp=fopen(filename,"w");
    if(!fp) {
	LALPrintError("Couldn't open file %s for writing\n",filename);
	free_toplist(&tl);
	free_toplist(&tl2);
	return(-2);
    }
    
    fprintf(stderr,"writing...\n");
    if(write_toplist_to_fp(tl, fp, &checksum)<0) {
	LALPrintError("Couldn't write toplist\n",filename);
	fclose(fp);
	free_toplist(&tl);
	free_toplist(&tl2);
	return(-2);
    }

    fprintf(stderr,"checksum: %d\n",checksum);

    fclose(fp);

    sort_toplist(tl2);

    fprintf(stderr,"comparing...\n");
    for(i=0;i<n;i++)
	if((tl->sorted[i]->Freq  - tl2->sorted[i]->Freq > epsilon) ||
	   (tl->sorted[i]->f1dot - tl2->sorted[i]->f1dot > epsilon) ||
	   (tl->sorted[i]->Alpha - tl2->sorted[i]->Alpha > epsilon) ||
	   (tl->sorted[i]->Delta - tl2->sorted[i]->Delta > epsilon) ||
	   (tl->sorted[i]->Nbins - tl2->sorted[i]->Nbins > epsilon) ||
	   (tl->sorted[i]->mean  - tl2->sorted[i]->mean > epsilon) ||
	   (tl->sorted[i]->std   - tl2->sorted[i]->std > epsilon) ||
	   (tl->sorted[i]->max   - tl2->sorted[i]->max > epsilon)) {

	    LALPrintError("line %d differs\n",i);
	    fprintf(stderr,"%e %e %e %e %d %e %e %e\n",
		    tl->sorted[i]->Freq,
		    tl->sorted[i]->f1dot,
		    tl->sorted[i]->Alpha,
		    tl->sorted[i]->Delta,
		    tl->sorted[i]->Nbins,
		    tl->sorted[i]->mean,
		    tl->sorted[i]->std,
		    tl->sorted[i]->max);
	    fprintf(stderr,"%e %e %e %e %d %e %e %e\n",
		    tl2->sorted[i]->Freq,
		    tl2->sorted[i]->f1dot,
		    tl2->sorted[i]->Alpha,
		    tl2->sorted[i]->Delta,
		    tl2->sorted[i]->Nbins,
		    tl2->sorted[i]->mean,
		    tl2->sorted[i]->std,
		    tl2->sorted[i]->max);
	}

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

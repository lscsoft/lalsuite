#include "bayespputils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

char * SPIN_PAR_LIST="chain log(post)     prior    mchirp       eta              time    log(dist)     RA  sin(dec) cos(iota)   phi0   psi  a1 cos(theta1) phi1 a2 cos(theta2) phi2";
char * NOSPIN_PAR_LIST="chain log(post)     prior    mchirp       eta              time    log(dist)     RA  sin(dec) cos(iota)   phi0   psi";
char * LALINFERENCE_NOSPIN_PAR_LIST="chain logl     prior    iota       psi       dec      RA      dist   phi0       time       eta        mchirp	temperature	MPIrank";
char * LALINFERENCE_SPIN_PAR_LIST="chain logl     prior    phi2   theta2   a2   phi1   theta1   a1    iota       psi       dec      RA      dist   phi0       time       eta        mchirp	temperature	MPIrank";
char * NOSPIN_COMMON_PAR_LIST="chain logl     prior    mchirp       eta              time    dist     RA  sin(dec) cos(iota)   phi0   psi";
char * SPIN_COMMON_PAR_LIST="chain logl     prior    mchirp       eta              time    dist     RA  sin(dec) cos(iota)   phi0   psi  a1 cos(theta1) phi1 a2 cos(theta2) phi2";

int spin;
double deltaLogL;
int thin;

double findMaxLogL(FILE * input, double maxLogL)
{
    int i;
    char line[MAX_REC_LEN];
    char* dummy;
    double blah, logL;
    for(i=1; i<=14; i++){
        dummy=fgets(line, MAX_REC_LEN, input);
    }
    while((dummy=fgets(line, MAX_REC_LEN, input))){
        sscanf(line, "%lf %lf", &blah, &logL);
        if(logL>maxLogL) maxLogL=logL;
    }
    return maxLogL;
}

void printBurnIn(FILE * input, FILE * output, double maxLogL, int chain, double * bayes, int * numpoints)
{
    int i, n;
    int numentries;
	if(spin==0)
		numentries=NOSPIN_PAR_NUM+3;
	if(spin==1)
		numentries=SPIN_PAR_NUM+3;
	if(spin==2)
		numentries=LALINFERENCE_NOSPIN_PAR_NUM+3;
	if(spin==3)
		numentries=LALINFERENCE_SPIN_PAR_NUM+3;
	if(spin==4)
		numentries=NOSPIN_COMMON_PAR_NUM+3;
	if(spin==5)
		numentries=SPIN_COMMON_PAR_NUM+3;
	
    int burnedIn=0;
    char buffer[MAX_REC_LEN];
    char * line=buffer;
    double values[LALINFERENCE_SPIN_PAR_NUM+3];

    char* dummy;

    for(i=1; i<=14; i++){
        dummy=fgets(buffer, MAX_REC_LEN, input);
    }

    while(fgets(buffer, MAX_REC_LEN, input)){
        line=buffer;
        //fprintf(stdout, "|%s|\n", line); fflush(stdout);
        for(i=0; i<numentries; i++){
            sscanf(line, "%lf%n", &(values[i]),&n);
            line+=n;
            //fprintf(stdout, "%p %lg ", line, values[i]); fflush(stdout);
        }
        //fprintf(stdout, "%p\n", line); fflush(stdout);

        if(!burnedIn && values[1]<(maxLogL-deltaLogL)) continue;
        if(!burnedIn) burnedIn=1;
		
		if(spin==4 || spin==5){
			values[3]=pow(values[3],6);
			values[6]=pow(values[6],1.0/3.0);
		}
		
        (*bayes)+=1.0/exp(values[1]);
        (*numpoints)++;
        //printf("%d: %g to %g\n", *numpoints, values[1], 1.0/exp(values[1]));
        fprintf(output, "%d ", chain);
        for(i=1; i<numentries; i++){
            fprintf(output, "%.13lf ", values[i]);
            //fprintf(stdout, "%lg ", values[i]);fflush(stdout);
        }
        fprintf(output,"\n");
		for(i=0;i<thin;i++) fgets(buffer, MAX_REC_LEN, input);
    }
}

int BurnIn(BurnInInput* inputs,BurnInOutput* outputs){

    FILE * output = fopen(inputs->output_file_name, "w");

    int i;

    double maxLogL=0;
    double bayes=0;
    int numpoints=0;
    FILE * input;

    for(i=0; i<inputs->nfiles; i++){
    	input=fopen(inputs->files[i], "r");
        maxLogL=findMaxLogL(input, maxLogL);
        fclose(input);
    }
    printf("maxLogL: %lg\n", maxLogL);

	if(spin==0)
		fprintf(output, "%s\n", NOSPIN_PAR_LIST);
	if(spin==1)
		fprintf(output, "%s\n", SPIN_PAR_LIST);
	if(spin==2)
		fprintf(output, "%s\n", LALINFERENCE_NOSPIN_PAR_LIST);
	if(spin==3)
		fprintf(output, "%s\n", LALINFERENCE_SPIN_PAR_LIST);
	if(spin==4)
		fprintf(output, "%s\n", NOSPIN_COMMON_PAR_LIST);
	if(spin==5)
		fprintf(output, "%s\n", SPIN_COMMON_PAR_LIST);
	
    for(i=0; i<inputs->nfiles; i++){
        input=fopen(inputs->files[i], "r");
        printf("Parsing file %s\n", inputs->files[i]);
        printBurnIn(input,output, maxLogL, i, &bayes, &numpoints);
        printf("numpoints: %d \n",numpoints);
        fclose(input);
    }
    fclose(output);

    bayes=1.0/bayes*numpoints;
    printf("Bayes factor according to harmonic mean of %d points is %g [log(B)=%g]\n",
    numpoints, bayes, log(bayes));

    outputs->bayesFactorFromHarmonicMean=log(bayes);

    return 0;
}


int main(int argc, char * argv [])
{
	if(argc!=6){
		printf("Usage: burnin <list file>\n");
		printf("              <0 for SPINspiral original non-spinning, 1 for SPINspiral original spinning,\n");
		printf("               2 for LALInference non-spinning, 3 for LALInference spinning,\n");
		printf("               4 for SPINspiral common input non-spinning, 5 for SPINspiral common input spinning>\n");
		printf("              <delta log(L) for burnin>\n"); 
		printf("              <output file>\n"); 
		printf("              <thinning>\n");
		
		exit(1);
    }
    FILE * list = fopen(argv[1],"r");	// File produced via "ls SPINspiral.output.*.00 > list"
    spin = atoi(argv[2]);
    deltaLogL = atof(argv[3]);
	thin = atoi(argv[5]);

    char line[MAX_REC_LEN];

    int i;
    char** files = calloc(MAX_IN_FILES, sizeof(char*));
    for (i = 0; i < MAX_IN_FILES; ++i) files[i] = malloc(MAX_REC_LEN);

    int nfiles=0;
    while(fgets(line, MAX_REC_LEN, list) && nfiles<MAX_IN_FILES){
        sscanf(line, "%s", files[nfiles]);
        //printf("%s\n", files[nfiles]);
        nfiles++;
    }

    fclose(list);

    BurnInInput* input=(BurnInInput*)malloc(sizeof(BurnInInput));
    input->files=files;
    input->nfiles=nfiles;
    input->spin=spin;
    input->deltaLogL=deltaLogL;
	input->thin=thin;
    input->output_file_name=argv[4];

    BurnInOutput* output=(BurnInOutput*)malloc(sizeof(BurnInOutput));

    int burnin_exit = BurnIn(input,output);

    return burnin_exit;

}

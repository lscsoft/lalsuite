#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> 

#define EOF_MARKER 26    /* Decimal code of DOS end-of-file marker */
#define MAX_REC_LEN 2048 /* Maximum size of input buffer */
#define MAX_IN_FILES 128

#define SPIN 1
#define NOSPIN 0

#define SPIN_PAR_NUM 15 //+cycles,posterior,prior
#define NOSPIN_PAR_NUM 9  //+cycles, posterior, prior

char * SPIN_PAR_LIST="chain log(post)     prior    mchirp       eta              time    log(dist)     RA  sin(dec) cos(iota)   phi0   psi  a1 cos(theta1) phi1 a2 cos(theta2) phi2";
char * NOSPIN_PAR_LIST="chain log(post)     prior    mchirp       eta              time    log(dist)     RA  sin(dec) cos(iota)   phi0   psi";

int spin;
double deltaLogL;
double findMaxLogL(FILE * input, double	maxLogL);
void printBurnIn(FILE * input, FILE * output, double maxLogL, int chain);

int main(int argc, char * argv [])
{
  if(argc!=5){
	printf("Usage: burnin <list file> <0 for non-spinning, 1 for spinning> <delta log(L) for burnin> <output file>\n");
  	exit(1);
  }
  FILE * list = fopen(argv[1],"r");	// File produced via "ls SPINspiral.output.*.00 > list"
  spin = atoi(argv[2]);
  deltaLogL = atof(argv[3]);
  FILE * output = fopen(argv[4], "w");

  char line[MAX_REC_LEN];
  char files[MAX_IN_FILES][MAX_REC_LEN];
  int i, nfiles=0;
  while(fgets(line, MAX_REC_LEN, list) && nfiles<MAX_IN_FILES){
	sscanf(line, "%s", files[nfiles]);
	//printf("%s\n", files[nfiles]);
	nfiles++;
  }
  printf("nfiles: %d\n", nfiles);
  fclose(list);

  double maxLogL=0;
  FILE * input;
  for(i=0; i<nfiles; i++){
	input=fopen(files[i], "r");
        maxLogL=findMaxLogL(input, maxLogL);
	fclose(input);
  }
  printf("maxLogL: %lg\n", maxLogL);

  if(spin)
  	fprintf(output, "%s\n", SPIN_PAR_LIST);
  else
  	fprintf(output, "%s\n", NOSPIN_PAR_LIST);
  for(i=0; i<nfiles; i++){
        input=fopen(files[i], "r");
	printf("Parsing file %s\n", files[i]);
	printBurnIn(input, output, maxLogL, i);
        fclose(input);
  }
  fclose(output);

  return 0;
}        
  
double findMaxLogL(FILE * input, double maxLogL)
{ 	
  int i;
  char line[MAX_REC_LEN];
  double blah, logL;
  for(i=1; i<=14; i++){
	fgets(line, MAX_REC_LEN, input);
  }
  while(fgets(line, MAX_REC_LEN, input)){
	sscanf(line, "%lf %lf", &blah, &logL);
	if(logL>maxLogL) maxLogL=logL;
  }
  return maxLogL;
}

void printBurnIn(FILE * input, FILE * output, double maxLogL, int chain)
{ 
  int i, n;
  int numentries=spin?SPIN_PAR_NUM+3:NOSPIN_PAR_NUM+3;
  int burnedIn=0;
  char buffer[MAX_REC_LEN];
  char * line=buffer;
  double values[SPIN_PAR_NUM+3];
  
  for(i=1; i<=14; i++){
	fgets(buffer, MAX_REC_LEN, input);
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

	fprintf(output, "%d ", chain);
	for(i=1; i<numentries; i++){
                fprintf(output, "%.13lf ", values[i]); 
//fprintf(stdout, "%lg ", values[i]);fflush(stdout);
        }
	fprintf(output,"\n");	
  }
}

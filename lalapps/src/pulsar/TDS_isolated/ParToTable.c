/* Matt Pitkin 22/06/04 */
/*
$Id$
*/
/* ParToTable.c
   Preprocessor code to take inputs from a TEMPO .par file and output the values
   into a table with a specific format like that used for S2 (don't know if I
	 think this is a good idea to do this for the actual analysis when we could 
	 just have all the .par files in a directory to read) */
	 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <dirent.h>

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>

#include "BinaryPulsarTiming.h"

int main(int argc, char *argv[]){
	static LALStatus status;
	char pAndP[256];
		
	struct dirent **namelist;
	int n;
	int i;
	
	char pardir[256]; /* .par file directory */
	
	FILE *fp; /* output file */
	
	char *loc;
	
	/* input .par file directory as first command line arg */
	sprintf(pardir, "%s", argv[1]);
	
	/* scan in names of all files in dir */
	n = scandir(pardir, &namelist, 0, alphasort);
	if(n<0)
		perror("scandir");
	
	fp = fopen("pulsar.list", "w");
	
	/* first 2 dir are . and .. */
	free(namelist[0]);
	free(namelist[1]);
	
	for(i=2;i<n;i++){
		int length;
		BinaryPulsarParams params;
		char psrname[15]="";
		char psr[15]="";
		
		/* read in data from .par files and output */
		sprintf(pAndP, "%s/%s", pardir, namelist[i]->d_name);
		
		/* check that file is a .par file */
		if(!strstr(namelist[i]->d_name, ".par")){
			fprintf(stderr, "%s is not a .par file, try the next one!.\n",
			namelist[i]->d_name);
			continue;
		}
		
		loc = strchr(namelist[i]->d_name, '.');
		
		length = loc - namelist[i]->d_name;
		
		strncpy(psr, namelist[i]->d_name, length);
		sprintf(psrname, "J%s", psr);
		
		free(namelist[i]);
		
		LALReadTEMPOParFile(&status, &params, pAndP);
		
		if(params.model == NULL)
			params.model = "*";
		else{
			/* print out max deltaf due to binary orbit */
			double df;
			
			df = LAL_TWOPI*params.f0*params.x/(params.Pb*86400.0);
			
			fprintf(stderr, "Max df for %s is %f Hz.\n", psrname, df);
		}
		
		fprintf(fp,
	"%d\t%s\t%.11lf\t%.11lf\t%.4e\t%.4e\t%.6lf\t%.16lf\t%.13e\t%.13e\t%.6lf\t%s\t%.9lf\t%.9lf\t%.12lf\t%.5e\t%.10lf\t%.5e\t%.7lf\t%.5e\t%.8e\t%.6e\t%.5e\t%.6e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5f\t%.5f\t%.5f\t%.5e\t%.5e\t%.5e\t%.5e\n", 
		i-1, psrname, params.ra, params.dec, params.pmra, params.pmdec, 
		params.posepoch, params.f0,params.f1, params.f2, params.pepoch, 
		params.model, params.T0, params.Tasc, params.Pb, params.Pbdot, params.x, 
		params.xdot, params.w0, params.wdot, params.e, params.edot, params.eps1, 
		params.eps1dot, params.eps2, params.eps2dot, params.xpbdot, params.gamma, 
		params.s, params.M, params.m2, params.dr, params.dth, params.a0, params.b0);		
	}
	
	free(namelist);
	
	fclose(fp);
	
	return 0;
}

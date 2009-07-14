#include <stdio.h>
#include "LALInference.h"


LALVariables variables;
REAL4 number,five;
ProcessParamsTable *ppt, *ptr;
int i;


int main(int argc, char *argv[]){
	number = 10.0;
	five=5.0;
	variables.head=NULL;
	variables.dimension=0;
	addVariable(&variables,"number",&number,REAL4_t);
	number=*(REAL4 *)getVariable(&variables,"number");
	fprintf(stdout,"Got %lf\n",number);
	setVariable(&variables,"number",&five);
	number=*(REAL4 *)getVariable(&variables,"number");
	fprintf(stdout,"Got %lf\n",number);
	fprintf(stdout,"Checkvariable?: %i\n",checkVariable(&variables,"number"));	
	removeVariable(&variables,"number");
	fprintf(stdout,"Removed, Checkvariable?: %i\n",checkVariable(&variables,"number"));
	destroyVariables(&variables);
  ppt = (ProcessParamsTable*) parseCommandLine(argc,argv);
  printf("parsed command line arguments:\n");
  ptr = ppt;
  i=1;
  while (ptr != NULL){
    printf(" (%d)  %s  %s  %s  \"%s\"\n", i, ptr->program, ptr->param, ptr->type, ptr->value);
    ptr = ptr->next;
    ++i;
  }

	
	LALVariables param;
	param.head=NULL;
	param.dimension=5;
	REAL4 m1 = 1.4;
	addVariable(&variables,"m1",&m1,REAL4_t);
	REAL4 m2 = 1.4;
	addVariable(&variables,"m2",&m2,REAL4_t);
	REAL4 dist = 8.5;
	addVariable(&variables,"dist",&m1,REAL4_t);
	REAL4 inc = 0.0;
	addVariable(&variables,"inc",&m1,REAL4_t);
	REAL4 phii = 0.0;
	addVariable(&variables,"phii",&m1,REAL4_t);
	
	
	LALIFOData ifo;
	ifo.fLow=40.0;
	ifo.fHigh=500.0;
	//REAL8TimeSeries timeData;
	//REAL8TimeSeries timeModelhPlus;
	//REAL8TimeSeries timeModelhCross;
	
	LALTemplateWrapper(&variables, &ifo);

	
}

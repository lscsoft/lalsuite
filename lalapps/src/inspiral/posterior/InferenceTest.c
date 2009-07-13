#include <stdio.h>
#include "LALInference.h"


LALVariables variables;
REAL8 number,five;

int main(){
	number = 10.0;
	five=5.0;
	char numberstr[]="number";
	variables.head=NULL;
	variables.dimension=0;
	addVariable(&variables,numberstr,&number,REAL8_t);
	number=*(REAL8 *)getVariable(&variables,numberstr);
	fprintf(stdout,"Got %lf\n",number);
	setVariable(&variables,numberstr,&five);
	number=*(REAL8 *)getVariable(&variables,numberstr);
	fprintf(stdout,"Got %lf\n",number);
	
}
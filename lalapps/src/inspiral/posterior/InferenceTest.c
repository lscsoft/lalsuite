/* 
 *  InferenceTest.c:  Bayesian Followup function testing site
 *
 *  Copyright (C) 2009 Ilya Mandel, Vivien Raymond, Christian Roever, Marc van der Sluys and John Veitch
 *
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

#include <stdio.h>
#include "LALInference.h"


LALVariables variables;
LALVariables variables2;

REAL4 number,five;
ProcessParamsTable *ppt, *ptr;
LALInferenceRunState *runstate=NULL;
int i;


int main(int argc, char *argv[]){
  /* test "LALVariables" stuff: */
  //LALIFOData *IFOdata=NULL;
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
  printVariables(&variables);
  copyVariables(&variables, &variables2);
  printVariables(&variables2);
  fprintf(stdout,"compareVariables?: %i\n",
          compareVariables(&variables,&variables2));

  removeVariable(&variables,"number");
  fprintf(stdout,"Removed, Checkvariable?: %i\n",checkVariable(&variables,"number"));
  
  fprintf(stdout,"compareVariables?: %i\n",
          compareVariables(&variables,&variables2));
  destroyVariables(&variables);
  destroyVariables(&variables2);
  printVariables(&variables2);
  
  /* test "parseCommandLine()" function: */
  ppt = (ProcessParamsTable*) parseCommandLine(argc,argv);
  printf("parsed command line arguments:\n");
  ptr = ppt;
  i=1;
  while (ptr != NULL){
    printf(" (%d)  %s  %s  %s  \"%s\"\n", i, ptr->program, ptr->param, ptr->type, ptr->value);
    ptr = ptr->next;
    ++i;
  }
	
  /* Test the data setup */
  runstate = initialize(ppt);

  //IFOdata=ReadData(ppt);

  if(runstate->data) fprintf(stdout,"Successfully read in the data!\n");

 // IFOdata->modelParams=calloc(1, sizeof(LALVariables));
	
 // IFOdata->modelParams->head=NULL;
 // IFOdata->modelParams->dimension=0;
 // REAL4 m1 = 1.4;
 // addVariable(IFOdata->modelParams,"m1",&m1,REAL4_t);
 // REAL4 m2 = 1.4;
 // addVariable(IFOdata->modelParams,"m2",&m2,REAL4_t);
 // REAL4 inc = 0.0;
 // addVariable(IFOdata->modelParams,"inc",&inc,REAL4_t);
 // REAL4 phii = 0.0;
 // addVariable(IFOdata->modelParams,"phii",&phii,REAL4_t);
	
  //REAL8TimeSeries timeModelhPlus;
  //REAL8TimeSeries timeModelhCross;
  
	  LALIFOData *ifo=NULL;
	
  LALTemplateWrapper(ifo);
  return 0;
}

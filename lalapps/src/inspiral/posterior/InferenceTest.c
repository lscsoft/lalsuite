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

/* example command line: */
/* ./InferenceTest --IFO [H1] --cache [/Users/john/data/triple/H1/frames.cache] --PSDstart 864162143.0 --PSDlength 1000 --srate 1024 --seglen 10 --trigtime 864162943.0  */

#include <stdio.h>
#include <lal/Date.h>
#include "LALInference.h"


LALVariables variables;
LALVariables variables2;

REAL4 number,five;
REAL8 numberR8;
INT4 numberI4;
INT8 numberI8;
COMPLEX8 numberC8;
COMPLEX16 numberC16;

ProcessParamsTable *ppt, *ptr;
LALInferenceRunState *runstate=NULL;
int i;


int main(int argc, char *argv[]){
  /* test "LALVariables" stuff: */
  number = 10.0;
  five=5.0;
  variables.head=NULL;
  variables.dimension=0;
  addVariable(&variables, "number", &number, REAL4_t);
  numberR8 = 7.0;
  addVariable(&variables, "seven", &numberR8, REAL8_t);
  numberR8 = LAL_PI;
  addVariable(&variables, "pi", &numberR8, REAL8_t);
  numberI4 = 123;
  addVariable(&variables, "small", &numberI4, INT4_t);
  numberI8 = 256*256*256*64;
  addVariable(&variables, "large", &numberI8, INT8_t);
  numberC8.re = 2.0;  numberC8.im = 3.0;
  addVariable(&variables, "complex1", &numberC8, COMPLEX8_t);
  numberC16.re = 1.23;  numberC16.im = -3.45;
  addVariable(&variables, "complex2", &numberC16, COMPLEX16_t);

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
  numberC16.im = 4.56;
  setVariable(&variables2,"complex2",&numberC16);
  fprintf(stdout,"compareVariables?: %i\n",
          compareVariables(&variables,&variables2));
  numberC16.im = -3.45;
  setVariable(&variables2,"complex2",&numberC16);
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

	
  /* Test the data initialisation &c. */
  runstate = initialize(ppt);

  
  if(runstate->data) {
    fprintf(stdout," data found --> trying some template computations etc.\n");

    REAL4 m1 = 1.4;
    addVariable(runstate->data->modelParams,"m1",&m1,REAL4_t);
    REAL4 m2 = 1.4;
    addVariable(runstate->data->modelParams,"m2",&m2,REAL4_t);
    REAL4 inc = 0.0;
    addVariable(runstate->data->modelParams,"inc",&inc,REAL4_t);
    REAL4 phii = 0.0;
    addVariable(runstate->data->modelParams,"phii",&phii,REAL4_t);
    LALTemplateGeneratePPN(runstate->data);

    /* templateStatPhase() test: */
    REAL8 mc   = 2.0;
    REAL8 eta  = 0.24;
    REAL8 iota = 1.0;
    REAL8 phi  = 2.0;
    REAL8 tc   = XLALGPSGetREAL8(&(runstate->data->timeData->epoch)) + (runstate->data->timeData->data->length / runstate->data->timeData->deltaT) - 1.0;
    destroyVariables(runstate->data->modelParams);
    addVariable(runstate->data->modelParams, "chirpmass",   &mc,   REAL8_t);
    addVariable(runstate->data->modelParams, "massratio",   &eta,  REAL8_t);
    addVariable(runstate->data->modelParams, "inclination", &iota, REAL8_t);
    addVariable(runstate->data->modelParams, "phase",       &phi,  REAL8_t);
    addVariable(runstate->data->modelParams, "time",        &tc,   REAL8_t);
    printVariables(runstate->data->modelParams);
    templateStatPhase(runstate->data);
  }


  printf(" done.\n");
  return 0;
}

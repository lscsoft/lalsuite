/* 

LALInference.c:	Bayesian Followup functions

Copyright 2009 Ilya Mandel, Vivien Raymond, Christian Roever, Marc van der Sluys and John Veitch

*/

#include <stdio.h>
#include <stdlib.h>
#include <LAL/LALInspiral.h>
#include "LALInference.h"

size_t typeSize[]={sizeof(REAL8),sizeof(REAL4),sizeof(gsl_matrix *)};


void die(char *message){
 fprintf(stderr,message);
 exit(1);
}

/* ============ Accessor functions for the Variable structure ========== */

LALVariableItem *getItem(LALVariables *vars,const char *name){
LALVariableItem *this=vars->head;
while(this!=NULL) 
 { if(!strcmp(this->name,name)) break;
   else this=this->next;
 }
return(this);
}

void *getVariable(LALVariables * vars,const char * name){
/* Return the value of variable name from the vars structure by walking the list*/
LALVariableItem *item;
item=getItem(vars,name);
if(!item) die("Error: variable not found in getVariable.\n");
return(item->value);
}

void setVariable(LALVariables * vars,const char * name, void *value){
/* Set the value of variable name in the vars structure to value */
LALVariableItem *item;
item=getItem(vars,name);
if(!item) die("Error: variable not found in setVariable.\n");
memcpy(item->value,value,typeSize[item->type]);

return;
}

void addVariable(LALVariables * vars,const char * name, void *value, VariableType type){
/* Add the variable name with type type and value value to vars */
 /* Check the name doesn't already exist -TODO */

 LALVariableItem *new=malloc(sizeof(LALVariableItem));
 memset(new,sizeof(LALVariableItem),0);
 new->value = (void *)malloc(typeSize[type]);
 if(new==NULL||new->value==NULL) die("Unable to allocate memory for list item\n");
 
 memcpy(new->name,name,VARNAME_MAX);
 new->type = type;
 memcpy(new->value,value,typeSize[type]);
 new->next = vars->head;
 vars->head = new;
 vars->dimension++;
 return;
}
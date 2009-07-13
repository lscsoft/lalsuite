/* 

LALInference.c:	Bayesian Followup functions

Copyright 2009 Ilya Mandel, Vivien Raymond, Christian Roever, Marc van der Sluys and John Veitch

*/

#include <stdio.h>
#include <stdlib.h>
#include <LAL/LALInspiral.h>
#include "LALInference.h"

void die(char *message){
 fprintf(stderr,message);
 exit(1);
}

/* ============ Accessor functions for the Variable structure ========== */

LALVariableItem *getItem(LALVariables *vars,char *name){
LALVariableItem *this=vars->head;
while(!strcmp(this->name,name) && this!=NULL) this=this->next;
return(this);
}

void *getVariable(LALVariables * vars, char * name){
/* Return the value of variable name from the vars structure by walking the list*/
LALVariableItem *item;
item=getItem(vars,name);
if(!item) die("Error: variable not found in getVariable.\n");
return(item->value);
}

void setVariable(LALVariables * vars, char * name, void value){
/* Set the value of variable name in the vars structure to value */
LALVariableItem *item;
item=getItem(vars,name);
if(!item) die("Error: variable not found in setVariable.\n");
*(item->value)=value;
return;
}

void addVariable(LALVariables * vars, char * name, void value, VariableType type){
/* Add the variable name with type type and value value to vars */
 /* Check the name doesn't already exist -TODO */

 LALVariableItem *new=malloc(sizeof(LALVariableItem));
 new.value = malloc(typeSize[type]);
 if(new==NULL||new.value==NULL) die("Unable to allocate memory for list item\n");
 
 memcpy(new.name,name,NAME_MAX);
 new.type = type;
 *new.value = value;
 new.next = vars->head;
 vars->head = &new;
 vars->dimension++;
 return;
}
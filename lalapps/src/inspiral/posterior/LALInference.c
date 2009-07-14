/* 

LALInference.c:	Bayesian Followup functions

Copyright 2009 Ilya Mandel, Vivien Raymond, Christian Roever, Marc van der Sluys and John Veitch

*/

#include <stdio.h>
#include <stdlib.h>
/*#include <LAL/LALInspiral.h>*/
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
 /* Check the name doesn't already exist */
 if(checkVariable(vars,name)) {fprintf(stderr,"addVariable: Cannot re-add %s\n",name); exit(1);}

 LALVariableItem *new=malloc(sizeof(LALVariableItem));
 memset(new,0,sizeof(LALVariableItem));
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

void removeVariable(LALVariables *vars,const char *name)
{
	LALVariableItem *this=vars->head;
	LALVariableItem *parent=NULL;
	
	while(this){
		if(!strcmp(this->name,name)) break;
		else {parent=this; this=this->next;}
	}
	if(!this) {fprintf(stderr,"removeVariable: warning, %s not found to remove\n",name); return;}
	if(!parent) vars->head=this->next;
	else parent->next=this->next;
	free(this->value);
	free(this);
	vars->dimension--;
	return;
}

int checkVariable(LALVariables *vars,const char *name){
/* Check for existance of name */
if(getItem(vars,name)) return 1;
else return 0;
}

/* Free the entire structure */
void destroyVariables(LALVariables *vars){
 LALVariableItem *this,*next;
 if(!vars) return;
 this=vars->head;
 if(this) next=this->next;
 while(this){
  free(this->value);
  free(this);
  this=next;
  if(this) next=this->next;
 }
 vars->head=NULL;
 vars->dimension=0;
return;
}


ProcessParamsTable *parseCommandLine(int argc, char *argv[])
/* parse command line and set up & fill in 'ProcessParamsTable' linked list.          */
/* If no command line arguments are supplied, the 'ProcessParamsTable' still contains */
/* one empty entry.                                                                   */
{
  int i, state=1;
  int dbldash;
  ProcessParamsTable *head, *ptr=NULL;
  // always (even for argc==1, i.e. no arguments) put one element in list:
  head = (ProcessParamsTable*) calloc(1, sizeof(ProcessParamsTable));
  strcpy(head->program, argv[0]);
  ptr = head;

  i=1;
  while ((i<argc) & (state<=3)) {
    // check for a double-dash at beginning of argument #i:
    dbldash = ((argv[i][0]=='-') && (argv[i][1]=='-'));
    // react depending on current state:
    if (state==1){
      if (dbldash) {
        strcpy(head->param, argv[i]);
        state = 2;
      }
      else {
        fprintf(stderr, " WARNING (1): orphaned command line argument '%s' in parseCommandLine().\n", argv[i]);
        state = 4;
      }
    } 
    else if (state==2) {
      if (dbldash) {
        ptr->next = (ProcessParamsTable*) calloc(1, sizeof(ProcessParamsTable));
        ptr = ptr->next;
        strcpy(ptr->program, argv[0]);
        strcpy(ptr->param, argv[i]);
      }
      else {
        state = 3;
        strcpy(ptr->value, argv[i]);          
        strcpy(ptr->type, "string");
      }
    }
    else if (state==3) {
      if (dbldash) {
        ptr->next = (ProcessParamsTable*) calloc(1, sizeof(ProcessParamsTable));
        ptr = ptr->next;
        strcpy(ptr->program, argv[0]);
        strcpy(ptr->param, argv[i]);
        state = 2;
      }
      else {
        fprintf(stderr, " WARNING (2): orphaned command line argument '%s' in parseCommandLine().\n", argv[i]);
        state = 4;
      }     
    }
    ++i;
  }
  if (0) { // check results:
    printf("-----\n");
    ptr = head;
    i=1;
    while (ptr != NULL){
      printf(" (%d)  %s  %s  %s  \"%s\"\n", i, ptr->program, ptr->param, ptr->type, ptr->value);
      ptr = ptr->next;
      ++i;
    }
  }
  return(head);
}

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

ProcessParamsTable *getProcParamVal(ProcessParamsTable *procparams,const char *name)
{
ProcessParamsTable *this=procparams;
while(this!=NULL) 
 { if(!strcmp(this->param,name)) break;
   else this=this->next;
 }
return(this);
}
 
void parseCharacterOptionString(char *input, char **strings[], int *n)
/* parses a character string (passed as one of the options) and decomposes   */
/* it into individual parameter character strings. Input is of the form      */
/*   input   :  "[one,two,three]"                                            */
/* and the resulting output is                                               */
/*   strings :  {"one", "two", "three"}                                      */
/* length of parameter names is by now limited to 512 characters.            */
/* (should 'theoretically' (untested) be able to digest white space as well. */
/* Irrelevant for command line options, though.) */
{
  int i,j,k,l;
  /* perform a very basic well-formedness-check and count number of parameters: */
  i=0; j=0;
  *n = 0;
  while (input[i] != '\0') {
    if ((j==0) & (input[i]=='[')) j=1;
    if ((j==1) & (input[i]==',')) ++*n;
    if ((j==1) & (input[i]==']')) {++*n; j=2;}
    ++i;
  }
  if (j!=2) printf(" : ERROR: argument vector '%s' not well-formed!\n", input);
  /* now allocate memory for results: */
  *strings  = (char**)  malloc(sizeof(char*) * (*n));
  for (i=0; i<(*n); ++i) (*strings)[i] = (char*) malloc(sizeof(char)*512);
  i=0; j=0; 
  k=0; /* string counter    */
  l=0; /* character counter */
  while ((input[i] != '\0') & (j<3)) {
    /* state transitions: */
    if ((j==0) & ((input[i]!='[') & (input[i]!=' '))) j=1;
    if (((j==1)|(j==2)) & (input[i]==',')) {(*strings)[k][l]='\0'; j=2; ++k; l=0;}
    if ((j==1) & (input[i]==' ')) j=2;
    if ((j==1) & (input[i]==']')) {(*strings)[k][l]='\0'; j=3;}
    if ((j==2) & (input[i]==']')) {(*strings)[k][l]='\0'; j=3;}
    if ((j==2) & ((input[i]!=']') & (input[i]!=',') & (input[i]!=' '))) j=1;
    /* actual copying: */
    if (j==1) {
      if (l>=511) {
        printf(" : WARNING: character argument too long!\n");
        printf(" : \"%s\"\n",(*strings)[k]);
        (*strings)[k][l] = input[i];
        ++l;
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
        fprintf(stderr, " WARNING: orphaned 1st command line argument \"%s\" in parseCommandLine().\n", argv[i]);
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
        fprintf(stderr, " WARNING: orphaned command line argument \"%s\" in parseCommandLine().\n", argv[i]);
        state = 4;
      }     
    }
    ++i;
  }
  if (state==4) die(" ERROR: failed parsing command line options.\n");
  return(head);
}

/* 

LALInference.c:	Bayesian Followup functions

Copyright 2009 Ilya Mandel, Vivien Raymond, Christian Roever, Marc van der Sluys and John Veitch

*/

#include <stdio.h>
#include <stdlib.h>
/*#include <LAL/LALInspiral.h>*/
#include "LALInference.h"
#include <lal/DetResponse.h>


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
/* Irrelevant for command line options, though.)                             */
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
      }
      else {
        (*strings)[k][l] = input[i];
        ++l;
      }
    }
    ++i;
  } 
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
    if (state==1){ // ('state 1' means handling very 1st argument)
      if (dbldash) {
        strcpy(head->param, argv[i]);
        strcpy(ptr->type, "string");
        state = 2;
      }
      else { // (very 1st argument needs to start with "--...")
        fprintf(stderr, " WARNING: orphaned first command line argument \"%s\" in parseCommandLine().\n", argv[i]);
        state = 4;
      }
    } 
    else if (state==2) { // ('state 2' means last entry was a parameter starting with "--")
      if (dbldash) {
        ptr->next = (ProcessParamsTable*) calloc(1, sizeof(ProcessParamsTable));
        ptr = ptr->next;
        strcpy(ptr->program, argv[0]);
        strcpy(ptr->param, argv[i]);
        strcpy(ptr->type, "string");
      }
      else {
        state = 3;
        strcpy(ptr->value, argv[i]);          
      }
    }
    else if (state==3) { // ('state 3' means last entry was a value)
      if (dbldash) {
        ptr->next = (ProcessParamsTable*) calloc(1, sizeof(ProcessParamsTable));
        ptr = ptr->next;
        strcpy(ptr->program, argv[0]);
        strcpy(ptr->param, argv[i]);
        strcpy(ptr->type, "string");
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




REAL8 FreqDomainLogLikelihood(LALVariables *currentParams, LALIFOData * data, 
                              LALTemplateFunction *template)
// (log-) likelihood function
{
  double Fplus, Fcross;
  double FplusScaled, FcrossScaled;
  double TwoDeltaToverN;
  double diff2;
  REAL8 likeli;
  REAL8 templateReal, templateImag;
  REAL8 dataReal, dataImag;
  int i, lower, upper;
  LALIFOData *dataPtr;
  double ra, dec, psi, distMpc, gmst;
  double GPSdouble;
  LIGOTimeGPS GPSlal;
  LALMSTUnitsAndAcc UandA;
  double chisquared;
  double timeshift;
  double deltaT=0.0;
  LALStatus status;

  // determine source's sky location & orientation parameters:
  ra        = *(REAL8*) getVariable(currentParams, "rightascension");
  dec       = *(REAL8*) getVariable(currentParams, "declination");
  psi       = *(REAL8*) getVariable(currentParams, "polarisation");
  GPSdouble = *(REAL8*) getVariable(currentParams, "time");
  distMpc   = *(REAL8*) getVariable(currentParams, "distance");

  GPSlal.gpsSeconds     = ((INT4) floor(GPSdouble));
  GPSlal.gpsNanoSeconds = 0; //((INT4) round(fmod(GPSdouble,1.0)*1e9));
  UandA.units = MST_RAD;
  UandA.accuracy = LALLEAPSEC_LOOSE;
  LALGPStoGMST1(&status, &gmst, &GPSlal, &UandA);

  chisquared = 0.0;

  // loop over data (different interferometers):
  dataPtr = data;
  while (dataPtr != NULL) {
    // compute template (deposited in elements of `data'):
    template(currentParams, data);

    // determine beam pattern response (F_plus and F_cross) for given Ifo:
    XLALComputeDetAMResponse(&Fplus, &Fcross,
                             dataPtr->detector->response,
			     ra, dec, psi, gmst);
    // signal arrival time (relative to geocentre);
    timeshift = XLALTimeDelayFromEarthCenter(dataPtr->detector->location,
                                             ra, dec, GPSlal);
    // (negative timeshift means signal arrives earlier than at geocenter etc.)

    // include distance (overall amplitude) effect in Fplus/Fcross:
    FplusScaled  = Fplus  / distMpc;
    FcrossScaled = Fcross / distMpc;

    // determine frequency range & loop over frequency bins:
    lower = ceil(dataPtr->fLow * dataPtr->timeData->data->length * dataPtr->timeData->deltaT);
    upper = floor(dataPtr->fHigh * dataPtr->timeData->data->length * dataPtr->timeData->deltaT);
    TwoDeltaToverN = 2.0 * dataPtr->timeData->deltaT / ((double)(upper-lower+1));
    for (i=lower; i<=upper; ++i){
      // derive template (involving location/orientation parameters) from given plus/cross waveforms:
      templateReal = FplusScaled * data->freqModelhPlus->data->data[i].re  +  FcrossScaled * data->freqModelhCross->data->data[i].re;
      templateImag = FplusScaled * data->freqModelhPlus->data->data[i].im  +  FcrossScaled * data->freqModelhCross->data->data[i].im;
      // (still need to do time-shifting!)

      // squared difference & `chi-squared':
      diff2 = TwoDeltaToverN * (pow(data->freqData->data->data[i].re - templateReal, 2.0) 
                                + pow(data->freqData->data->data[i].im - templateImag, 2.0));
      chisquared += (diff2 / data->oneSidedNoisePowerSpectrum->data->data[i]);
    }
    dataPtr = dataPtr->next;
  }
  likeli = -1.0 * chisquared;
  return(likeli);
}

/* 
 *  LALInference.c:  Bayesian Followup functions
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
#include <stdlib.h>
#include "LALInference.h"
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>


size_t typeSize[]={sizeof(REAL8),sizeof(REAL4),sizeof(gsl_matrix *)};


void die(char *message)
{
  fprintf(stderr, message);
  exit(1);
}



/* ============ Accessor functions for the Variable structure ========== */



LALVariableItem *getItem(LALVariables *vars,const char *name)
{
  LALVariableItem *this=vars->head;
  while(this!=NULL) { 
    if(!strcmp(this->name,name)) break;
    else this=this->next;
  }
  return(this);
}



void *getVariable(LALVariables * vars,const char * name)
/* Return the value of variable name from the vars structure by walking the list*/
{
  LALVariableItem *item;
  item=getItem(vars,name);
  if(!item) {
    fprintf(stderr, " ERROR in getVariable(): entry \"%s\" not found.\n", name);
    exit(1);
  }
  return(item->value);
}



void setVariable(LALVariables * vars,const char * name, void *value)
/* Set the value of variable name in the vars structure to value */
{
  LALVariableItem *item;
  item=getItem(vars,name);
  if(!item) {
    fprintf(stderr, " ERROR in setVariable(): entry \"%s\" not found.\n", name);
    exit(1);
  }
  memcpy(item->value,value,typeSize[item->type]);
  return;
}



void addVariable(LALVariables * vars,const char * name, void *value, VariableType type)
/* Add the variable name with type type and value value to vars */
{
  /* Check the name doesn't already exist */
  if(checkVariable(vars,name)) {fprintf(stderr," ERROR in addVariable(): Cannot re-add \"%s\".\n",name); exit(1);}

  LALVariableItem *new=malloc(sizeof(LALVariableItem));
  memset(new,0,sizeof(LALVariableItem));
  new->value = (void *)malloc(typeSize[type]);
  if(new==NULL||new->value==NULL) die(" ERROR in addVariable(): unable to allocate memory for list item.\n");
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
  if(!this) {fprintf(stderr," WARNING in removeVariable(): entry \"%s\" not found.\n",name); return;}
  if(!parent) vars->head=this->next;
  else parent->next=this->next;
  free(this->value);
  free(this);
  vars->dimension--;
  return;
}



int checkVariable(LALVariables *vars,const char *name)
/* Check for existance of name */
{
  if(getItem(vars,name)) return 1;
  else return 0;
}



void destroyVariables(LALVariables *vars)
/* Free the entire structure */
{
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



void copyVariables(LALVariables *origin, LALVariables *target)
/*  copy contents of "origin" over to "target"  */
{
  LALVariableItem *ptr;
  /* first dispose contents of "target": */
  destroyVariables(target);
  /* then copy over elements of "origin": */
  ptr = origin->head;
  while (ptr != NULL) {
    addVariable(target, ptr->name, ptr->value, ptr->type);
    ptr = ptr->next;
  }
  return;
}



void printVariables(LALVariables *var)
/* output contents of a 'LALVariables' structure       */
/* (by now only prints names and types, but no values) */
{
  LALVariableItem *ptr = var->head;
  fprintf(stderr, "LALVariables:\n");
  if (ptr==NULL) fprintf(stderr, "  <empty>\n");
  else {
    while (ptr != NULL) {
      fprintf(stderr, "  \"%s\" (type #%d)\n", ptr->name, ((int) ptr->type));
      ptr = ptr->next;
    }  
  }
  return;
}



int compareVariables(LALVariables *var1, LALVariables *var2)
/*  compare contents of "var1" and "var2".                       */
/*  returns zero for equal entries, and one if difference found. */
{
  int result = 0;
  LALVariableItem *ptr1 = var1->head;
  LALVariableItem *ptr2 = NULL;
  if (var1->dimension != var2->dimension) result = 1;  // differing dimension
  while ((ptr1 != NULL) & (result == 0)) {
    ptr2 = getItem(var2, ptr1->name);
    if (ptr2 != NULL) {  // corrsesponding entry exists; now compare type, then value:
      if (ptr2->type == ptr1->type) {  // entry type identical
        switch (ptr1->type) {  // do value comparison depending on type:
          case REAL8_t:
            result = ((*(REAL8 *) ptr2->value) != (*(REAL8 *) ptr1->value));
            break;
          case REAL4_t: 
            result = ((*(REAL4 *) ptr2->value) != (*(REAL4 *) ptr1->value));
            break;
          case gslMatrix_t: 
            fprintf(stderr, " WARNING: compareVariables() cannot yet compare \"gslMatrix\" type entries.\n");
            fprintf(stderr, "          (entry: \"%s\").\n", ptr1->name);
            fprintf(stderr, "          For now entries are by default assumed different.\n");
            result = 1;
            break;
          default:
            fprintf(stderr, " ERROR: encountered unknown LALVariables type in compareVariables()\n");
            fprintf(stderr, "        (entry: \"%s\").\n", ptr1->name);
            exit(1);
        }
      }
      else result = 1;  // same name but differing type
    }
    else result = 1;  // entry of given name doesn't exist in var2
    ptr1 = ptr1->next;
  }
  return(result);
}



ProcessParamsTable *getProcParamVal(ProcessParamsTable *procparams,const char *name)
{
  ProcessParamsTable *this=procparams;
  while(this!=NULL) { 
    if(!strcmp(this->param,name)) break;
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
  if (j!=2) printf(" ERROR: argument vector \"%s\" not well-formed!\n", input);
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
  /* always (even for argc==1, i.e. no arguments) put one element in list: */
  head = (ProcessParamsTable*) calloc(1, sizeof(ProcessParamsTable));
  strcpy(head->program, argv[0]);
  ptr = head;
  i=1;
  while ((i<argc) & (state<=3)) {
    /* check for a double-dash at beginning of argument #i: */
    dbldash = ((argv[i][0]=='-') && (argv[i][1]=='-'));
    /* react depending on current state: */
    if (state==1){ /* ('state 1' means handling very 1st argument) */
      if (dbldash) {
        strcpy(head->param, argv[i]);
        strcpy(ptr->type, "string");
        state = 2;
      }
      else { /* (very 1st argument needs to start with "--...") */
        fprintf(stderr, " WARNING: orphaned first command line argument \"%s\" in parseCommandLine().\n", argv[i]);
        state = 4;
      }
    } 
    else if (state==2) { /* ('state 2' means last entry was a parameter starting with "--") */
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
    else if (state==3) { /* ('state 3' means last entry was a value) */
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
/* (log-) likelihood function.                        */
/* Returns the non-normalised logarithmic likelihood. */
{
  double Fplus, Fcross;
  double FplusScaled, FcrossScaled;
  double TwoDeltaToverN;
  double diff2;
  REAL8 loglikeli;
  REAL8 plainTemplateReal, plainTemplateImag;
  REAL8 templateReal, templateImag;
  int i, lower, upper;
  LALIFOData *dataPtr;
  double ra, dec, psi, distMpc, gmst;
  double GPSdouble;
  LIGOTimeGPS GPSlal;
  LALMSTUnitsAndAcc UandA;
  double chisquared;
  double timedelay;  /* time delay b/w iterferometer & geocenter w.r.t. sky location */
  double timeshift;  /* time shift (not necessarily same as above)                   */
  double NDeltaT, twopit, f, re, im;
  double timeTmp;
  int different;
  LALStatus status;
  LALVariables intrinsicParams;

  /* determine source's sky location & orientation parameters: */
  ra        = *(REAL8*) getVariable(currentParams, "rightascension"); /* radian      */
  dec       = *(REAL8*) getVariable(currentParams, "declination");    /* radian      */
  psi       = *(REAL8*) getVariable(currentParams, "polarisation");   /* radian      */
  GPSdouble = *(REAL8*) getVariable(currentParams, "time");           /* GPS seconds */
  distMpc   = *(REAL8*) getVariable(currentParams, "distance");       /* Mpc         */

  GPSlal.gpsSeconds     = ((INT4) floor(GPSdouble));
  GPSlal.gpsNanoSeconds = 0; /*((INT4) round(fmod(GPSdouble,1.0)*1e9)); */
  UandA.units = MST_RAD;
  UandA.accuracy = LALLEAPSEC_LOOSE;
  LALGPStoGMST1(&status, &gmst, &GPSlal, &UandA);

  intrinsicParams.head = NULL;
  intrinsicParams.dimension = 0;
  copyVariables(currentParams, &intrinsicParams);
  removeVariable(&intrinsicParams, "rightascension");
  removeVariable(&intrinsicParams, "declination");
  removeVariable(&intrinsicParams, "polarisation");
  removeVariable(&intrinsicParams, "time");
  removeVariable(&intrinsicParams, "distance");

  chisquared = 0.0;
  /* loop over data (different interferometers): */
  dataPtr = data;
  while (dataPtr != NULL) {
    /* The parameters the Likelihood function can handle by itself   */
    /* (and which shouldn't affect the template function) are        */
    /* sky location (ra, dec), polarisation and signal arrival time. */
    /* Note that the template function still needs _some_ reasonable */
    /* arrival time parameter value (e.g. something like the trigger */
    /* value).                                                       */
    
    /* Compare parameter values with parameter values corresponding  */
    /* to currently stored template; ignore "time" variable:         */
    if (checkVariable(data->modelParams, "time")) {
      timeTmp = *(REAL8 *) getVariable(data->modelParams, "time");
      removeVariable(data->modelParams, "time");
    }
    else timeTmp = GPSdouble;
    different = compareVariables(data->modelParams, &intrinsicParams);
    /* "different" now may also mean that "data->modelParams" */
    /* wasn't allocated yet (as in the very 1st iteration).   */

    if (different) { /* template needs to be re-computed: */
      copyVariables(&intrinsicParams, data->modelParams);
      addVariable(data->modelParams, "time", &timeTmp, REAL8_t);
      template(data);
    }
    else { /* no re-computation necessary. Return back "time" value, do nothing else: */
      addVariable(data->modelParams, "time", &timeTmp, REAL8_t);
    }

    /*-- Template is now in data->freqModelhPlus and data->freqModelhCross. --*/
    /*-- (Either freshly computed or inherited.)                            --*/

    /* determine beam pattern response (F_plus and F_cross) for given Ifo: */
    XLALComputeDetAMResponse(&Fplus, &Fcross,
                             dataPtr->detector->response,
			     ra, dec, psi, gmst);
    /* signal arrival time (relative to geocenter); */
    timedelay = XLALTimeDelayFromEarthCenter(dataPtr->detector->location,
                                             ra, dec, &GPSlal);
    /* (negative timedelay means signal arrives earlier than at geocenter etc.) */

    /* amount by which to time-shift template (not necessarily same as above "timedelay"): */
    timeshift =  (GPSdouble - (*(REAL8*) getVariable(data->modelParams, "time"))) + timedelay;
    twopit    = LAL_TWOPI * timeshift;

    /* include distance (overall amplitude) effect in Fplus/Fcross: */
    FplusScaled  = Fplus  / distMpc;
    FcrossScaled = Fcross / distMpc;

    /* determine frequency range & loop over frequency bins: */
    NDeltaT = dataPtr->timeData->data->length * dataPtr->timeData->deltaT;
    lower = ceil(dataPtr->fLow * NDeltaT);
    upper = floor(dataPtr->fHigh * NDeltaT);
    TwoDeltaToverN = 2.0 * dataPtr->timeData->deltaT / ((double)(upper-lower+1));
    for (i=lower; i<=upper; ++i){
      /* derive template (involving location/orientation parameters) from given plus/cross waveforms: */
      plainTemplateReal = FplusScaled * data->freqModelhPlus->data->data[i].re  
                          +  FcrossScaled * data->freqModelhCross->data->data[i].re;
      plainTemplateImag = FplusScaled * data->freqModelhPlus->data->data[i].im  
                          +  FcrossScaled * data->freqModelhCross->data->data[i].im;

      /* do time-shifting: */
      f = ((double) i) / NDeltaT;
      /* real & imag parts of  exp(-2*pi*i*f*deltaT): */
      re = cos(twopit * f);
      im = - sin(twopit * f);
      templateReal = plainTemplateReal*re - plainTemplateImag*im;
      templateImag = plainTemplateReal*im + plainTemplateImag*re;

      /* squared difference & 'chi-squared': */
      diff2 = TwoDeltaToverN * (pow(data->freqData->data->data[i].re - templateReal, 2.0) 
                                + pow(data->freqData->data->data[i].im - templateImag, 2.0));
      chisquared += (diff2 / data->oneSidedNoisePowerSpectrum->data->data[i]);
    }
    dataPtr = dataPtr->next;
  }
  loglikeli = -1.0 * chisquared;
  return(loglikeli);
}

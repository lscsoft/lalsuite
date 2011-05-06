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
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/TimeFreqFFT.h>
#include <lal/VectorOps.h>
#include <lal/Date.h>
#include <lal/Sequence.h>

size_t typeSize[] = {sizeof(INT4), 
                     sizeof(INT8),
                     sizeof(UINT4),
                     sizeof(REAL4), 
                     sizeof(REAL8), 
                     sizeof(COMPLEX8), 
                     sizeof(COMPLEX16), 
                     sizeof(gsl_matrix *),
                     sizeof(REAL8Vector *),
                     sizeof(UINT4Vector *),
                     sizeof(CHAR *)};


void die(const char message[])
{
  fprintf(stderr, "%s", message);
  exit(1);
}



/* ============ Accessor functions for the Variable structure: ========== */



LALVariableItem *getItem(LALVariables *vars,const char *name)
/* (this function is only to be used internally) */
/* Returns pointer to item for given item name.  */
{
  if(vars==NULL) return NULL;
  LALVariableItem *this = vars->head;
  while (this != NULL) { 
    if (!strcmp(this->name,name)) break;
    else this = this->next;
  }
  return(this);
}


LALVariableItem *getItemNr(LALVariables *vars, int idx)
/* (this function is only to be used internally)  */
/* Returns pointer to item for given item number. */
{
  int i=1;
  if (idx < i) die(" Error in getItemNr(): requesting zero or negative idx entry.\n");
  LALVariableItem *this=vars->head;
  while (this != NULL) { 
    if (i == idx) break;
    else {
      this = this->next;
      ++i;
    }
  }
  return(this);
}

ParamVaryType getVariableVaryType(LALVariables *vars, const char *name)
{
	return (getItem(vars,name)->vary);
}


void *getVariable(LALVariables * vars,const char * name)
/* Return the value of variable name from the vars structure by walking the list */
{
  LALVariableItem *item;
  item=getItem(vars,name);
  if(!item) {
    fprintf(stderr, " ERROR in getVariable(): entry \"%s\" not found.\n", name);
    exit(1);
  }
  return(item->value);
}


INT4 getVariableDimension(LALVariables *vars)
{
  return(vars->dimension);
}


INT4 getVariableDimensionNonFixed(LALVariables *vars)
{
	INT4 count=0;
	LALVariableItem *ptr = vars->head;
	if (ptr==NULL) return count;
	else {
		/* loop over entries: */
		while (ptr != NULL) {
			/* print name: */
			if (ptr->vary != PARAM_FIXED) ++count;
			ptr = ptr->next;
		}  
	}
	return count;
}


VariableType getVariableType(LALVariables *vars, const char *name)
{
	return getItem(vars,name)->type;
}

VariableType getVariableTypeByIndex(LALVariables *vars, int idx)
/* Returns type of the i-th entry, */
/* where  1 <= idx <= dimension. */
{
  LALVariableItem *item;
  if ((idx < 1) | (idx > vars->dimension)){
    fprintf(stderr, " ERROR in getVariableTypeByIndex(...,idx=%d): idx needs to be 1 <= idx <= dimension = %d.\n", 
            idx, vars->dimension);
    exit(1);
  }
  item = getItemNr(vars, idx);
  return(item->type);
}


char *getVariableName(LALVariables *vars, int idx)
/* Returns (pointer to) the name of the i-th entry, */
/* where  1 <= idx <= dimension.                  */
{
  LALVariableItem *item;
  if ((idx < 1) | (idx > vars->dimension)){
    fprintf(stderr, " ERROR in getVariableName(...,idx=%d): idx needs to be 1 <= idx <= dimension = %d.\n", 
            idx, vars->dimension);
    exit(1);
  }
  item = getItemNr(vars, idx);
  return(item->name);
}


void setVariable(LALVariables * vars, const char * name, void *value)
/* Set the value of variable name in the vars structure to value */
{
  LALVariableItem *item;
  item=getItem(vars,name);
  if(!item) {
    fprintf(stderr, " ERROR in setVariable(): entry \"%s\" not found.\n", name);
    exit(1);
  }
  if (item->vary==PARAM_FIXED) return;
  memcpy(item->value,value,typeSize[item->type]);
  return;
}



void addVariable(LALVariables * vars, const char * name, void *value, VariableType type, ParamVaryType vary)
/* Add the variable name with type type and value value to vars */
/* If variable already exists, it will over-write the current value if type compatible*/
{
  LALVariableItem *old=NULL;
  /* Check the name doesn't already exist */
  if(checkVariable(vars,name)) {
	  old=getItem(vars,name);
	  if(old->type != type)
	  {fprintf(stderr," ERROR in addVariable(): Cannot re-add \"%s\" as previous definition has wrong type.\n",name); exit(1);}
	  setVariable(vars,name,value);
	  return;
  }
	
  if(!value) {fprintf(stderr,"Unable to access value through null pointer in addVariable, trying to add %s\n",name); exit(1);}

  LALVariableItem *new=malloc(sizeof(LALVariableItem));

  memset(new,0,sizeof(LALVariableItem));
	if(new) {
		new->value = (void *)malloc(typeSize[type]);
	}
  if(new==NULL||new->value==NULL) die(" ERROR in addVariable(): unable to allocate memory for list item.\n");
  memcpy(new->name,name,VARNAME_MAX);
  new->type = type;
  new->vary = vary;
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
  if(!origin)
  {
	  fprintf(stderr,"Unable to access origin pointer in copyVariables\n");
	  exit(1);
  }

  /* Make sure the structure is initialised */
	if(!target) fprintf(stderr,"ERROR: Unable to copy to uninitialised LALVariables structure\n");

	
  /* first dispose contents of "target" (if any): */
  destroyVariables(target);
  
	
  /* then copy over elements of "origin": */
  ptr = origin->head;
  if(!ptr)
  {
	  fprintf(stderr,"Bad LALVariable structure found while trying to copy\n");
	  exit(1);
  }
  while (ptr != NULL) {
	  if(!ptr->value || !ptr->name){
		  fprintf(stderr,"Badly formed LALVariableItem structure found in copyVariables!\n");
		  exit(1);
	  }
    addVariable(target, ptr->name, ptr->value, ptr->type, ptr->vary);
    ptr = ptr->next;
  }
  return;
}

void printVariables(LALVariables *var)
/* output contents of a 'LALVariables' structure       */
/* (by now only prints names and types, but no values) */
{
  LALVariableItem *ptr = var->head;
  fprintf(stdout, "LALVariables:\n");
  if (ptr==NULL) fprintf(stdout, "  <empty>\n");
  else {
    /* loop over entries: */
    while (ptr != NULL) {
      /* print name: */
      fprintf(stdout, "  \"%s\"", ptr->name); 
      /* print type: */
      fprintf(stdout, "  (type #%d, ", ((int) ptr->type));
      switch (ptr->type) {
        case INT4_t:
          fprintf(stdout, "'INT4'");
          break;
        case INT8_t:
          fprintf(stdout, "'INT8'");
          break;
		case UINT4_t:
		  fprintf(stdout, "'UINT4'");
		  break;			  
        case REAL4_t:
          fprintf(stdout, "'REAL4'");
          break;
        case REAL8_t:
          fprintf(stdout, "'REAL8'");
          break;
        case COMPLEX8_t:
          fprintf(stdout, "'COMPLEX8'");
          break;
        case COMPLEX16_t:
          fprintf(stdout, "'COMPLEX16'");
          break;
        case gslMatrix_t:
          fprintf(stdout, "'gslMatrix'");
          break;
        default:
          fprintf(stdout, "<unknown type>");
      }
      fprintf(stdout, ")  ");
      /* print value: */
      switch (ptr->type) {
        case INT4_t:
          fprintf(stdout, "%d", *(INT4 *) ptr->value);
          break;
        case INT8_t:
          fprintf(stdout, "%lld", *(INT8 *) ptr->value);
          break;
		case UINT4_t:
		  fprintf(stdout, "%ud", *(UINT4 *) ptr->value);
		  break;			  
        case REAL4_t:
          fprintf(stdout, "%.15lf", *(REAL4 *) ptr->value);
          break;
        case REAL8_t:
          fprintf(stdout, "%.15lf", *(REAL8 *) ptr->value);
          break;
        case COMPLEX8_t:
          fprintf(stdout, "%e + i*%e", 
                  (REAL4) ((COMPLEX8 *) ptr->value)->re, (REAL4) ((COMPLEX8 *) ptr->value)->im);
          break;
        case COMPLEX16_t:
          fprintf(stdout, "%e + i*%e", 
                  (REAL8) ((COMPLEX16 *) ptr->value)->re, (REAL8) ((COMPLEX16 *) ptr->value)->im);
          break;
        case gslMatrix_t:
          fprintf(stdout, "<can't print matrix>");          
          break;
        default:
          fprintf(stdout, "<can't print>");          
      }
      fprintf(stdout, "\n");
      ptr = ptr->next;
    }  
  }
  return;
}

void fprintSample(FILE *fp,LALVariables *sample){
	if(sample==NULL) return;
	LALVariableItem *ptr=sample->head;
	if(fp==NULL) return;
	while(ptr!=NULL) {
		switch (ptr->type) {
			case INT4_t:
				fprintf(fp, "%d", *(INT4 *) ptr->value);
				break;
			case INT8_t:
				fprintf(fp, "%lld", *(INT8 *) ptr->value);
				break;
			case UINT4_t:
				fprintf(fp, "%ud", *(UINT4 *) ptr->value);
				break;
			case REAL4_t:
				fprintf(fp, "%9.5f", *(REAL4 *) ptr->value);
				break;
			case REAL8_t:
				fprintf(fp, "%9.5lf", *(REAL8 *) ptr->value);
				break;
			case COMPLEX8_t:
				fprintf(fp, "%e + i*%e",
						(REAL4) ((COMPLEX8 *) ptr->value)->re, (REAL4) ((COMPLEX8 *) ptr->value)->im);
				break;
			case COMPLEX16_t:
				fprintf(fp, "%e + i*%e",
						(REAL8) ((COMPLEX16 *) ptr->value)->re, (REAL8) ((COMPLEX16 *) ptr->value)->im);
				break;
			case gslMatrix_t:
				fprintf(stdout, "<can't print matrix>");
				break;
			default:
				fprintf(stdout, "<can't print>");
			}
	
	fprintf(fp,"\t");
	ptr=ptr->next;
	}
	return;
}

void fprintSampleNonFixed(FILE *fp,LALVariables *sample){
	if(sample==NULL) return;
	LALVariableItem *ptr=sample->head;
	if(fp==NULL) return;
	while(ptr!=NULL) {
		if (ptr->vary != PARAM_FIXED) {
			switch (ptr->type) {
				case INT4_t:
					fprintf(fp, "%d", *(INT4 *) ptr->value);
					break;
				case INT8_t:
					fprintf(fp, "%lld", *(INT8 *) ptr->value);
					break;
				case UINT4_t:
					fprintf(fp, "%ud", *(UINT4 *) ptr->value);
					break;
				case REAL4_t:
					fprintf(fp, "%9.5f", *(REAL4 *) ptr->value);
					break;
				case REAL8_t:
					fprintf(fp, "%9.5f", *(REAL8 *) ptr->value);
					break;
				case COMPLEX8_t:
					fprintf(fp, "%e + i*%e",
							(REAL4) ((COMPLEX8 *) ptr->value)->re, (REAL4) ((COMPLEX8 *) ptr->value)->im);
					break;
				case COMPLEX16_t:
					fprintf(fp, "%e + i*%e",
							(REAL8) ((COMPLEX16 *) ptr->value)->re, (REAL8) ((COMPLEX16 *) ptr->value)->im);
					break;
				case gslMatrix_t:
					fprintf(stdout, "<can't print matrix>");
					break;
				default:
					fprintf(stdout, "<can't print>");
			}
		fprintf(fp,"\t");
		}
		ptr=ptr->next;
	}
	return;
}

const char *translateInternalToExternalParamName(const char *inName) {
  if (!strcmp(inName, "a_spin1")) {
    return "a1";
  } else if (!strcmp(inName, "a_spin2")) {
    return "a2";
  } else if (!strcmp(inName, "phi_spin1")) {
    return "phi1";
  } else if (!strcmp(inName, "phi_spin2")) {
    return "phi2";
  } else if (!strcmp(inName, "theta_spin1")) {
    return "theta1";
  } else if (!strcmp(inName, "theta_spin2")) {
    return "theta2";
  } else if (!strcmp(inName, "chirpmass")) {
    return "mc";
  } else if (!strcmp(inName, "massratio")) {
    return "eta";
  } else if (!strcmp(inName, "rightascension")) {
    return "ra";
  } else if (!strcmp(inName, "declination")) {
    return "dec";
  } else if (!strcmp(inName, "phase")) {
    return "phi_orb";
  } else if (!strcmp(inName, "polarisation")) {
    return "psi";
  } else if (!strcmp(inName, "inclination")) {
    return "iota";
  } else if (!strcmp(inName, "distance")) {
    return "dist";
  } else {
    return inName;
  }
}

int fprintParameterNonFixedHeaders(FILE *out, LALVariables *params) {
  LALVariableItem *head = params->head;

  while (head != NULL) {
    if (head->vary != PARAM_FIXED) {
      fprintf(out, "%s\t", translateInternalToExternalParamName(head->name));
    }
    head = head->next;
  }

  return 0;
}

int compareVariables(LALVariables *var1, LALVariables *var2)
/*  Compare contents of "var1" and "var2".                       */
/*  Returns zero for equal entries, and one if difference found. */
/*  Make sure to only call this function when all entries are    */
/*  actually comparable. For example, "gslMatrix" type entries   */
/*  cannot (yet?) be checked for equality.                       */
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
          case INT4_t: 
            result = ((*(INT4 *) ptr2->value) != (*(INT4 *) ptr1->value));
            break;
          case INT8_t: 
            result = ((*(INT8 *) ptr2->value) != (*(INT8 *) ptr1->value));
            break;
		  case UINT4_t: 
			result = ((*(UINT4 *) ptr2->value) != (*(UINT4 *) ptr1->value));
			break;
          case REAL4_t: 
            result = ((*(REAL4 *) ptr2->value) != (*(REAL4 *) ptr1->value));
            break;
          case REAL8_t:
            result = ((*(REAL8 *) ptr2->value) != (*(REAL8 *) ptr1->value));
            break;
          case COMPLEX8_t: 
            result = (((REAL4) ((COMPLEX8 *) ptr2->value)->re != (REAL4) ((COMPLEX8 *) ptr1->value)->re)
                      || ((REAL4) ((COMPLEX8 *) ptr2->value)->im != (REAL4) ((COMPLEX8 *) ptr1->value)->im));
            break;
          case COMPLEX16_t: 
            result = (((REAL8) ((COMPLEX16 *) ptr2->value)->re != (REAL8) ((COMPLEX16 *) ptr1->value)->re)
                      || ((REAL8) ((COMPLEX16 *) ptr2->value)->im != (REAL8) ((COMPLEX16 *) ptr1->value)->im));
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



/* ============ Command line parsing functions etc.: ========== */



ProcessParamsTable *getProcParamVal(ProcessParamsTable *procparams,const char *name)
/* Returns pointer to element "name" of the ProcessParamsTable, */
/* if present, and NULL otherwise.                              */
{
  ProcessParamsTable *this=procparams;
  while (this!=NULL) { 
    if (!strcmp(this->param, name)) break;
    else this=this->next;
  }
  return(this);
}



void parseCharacterOptionString(char *input, char **strings[], UINT4 *n)
/* parses a character string (passed as one of the options) and decomposes   */
/* it into individual parameter character strings. Input is of the form      */
/*   input   :  "[one,two,three]"                                            */
/* and the resulting output is                                               */
/*   strings :  {"one", "two", "three"}                                      */
/* length of parameter names is for now limited to 512 characters.           */
/* (should 'theoretically' (untested) be able to digest white space as well. */
/* Irrelevant for command line options, though.)                             */
{
  UINT4 i,j,k,l;
  /* perform a very basic well-formedness-check and count number of parameters: */
  i=0; j=0;
  *n = 0;
  while (input[i] != '\0') {
    if ((j==0) & (input[i]=='[')) j=1;
    if ((j==1) & (input[i]==',')) ++*n;
    if ((j==1) & (input[i]==']')) {++*n; j=2;}
    ++i;
  }
  if (j!=2) fprintf(stderr, " ERROR: argument vector \"%s\" not well-formed!\n", input);
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
        fprintf(stderr, " WARNING: character argument too long!\n");
        fprintf(stderr, " \"%s\"\n",(*strings)[k]);
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
  if (state==4) die(" ERROR in parseCommandLine(): failed parsing command line options.\n");
  return(head);
}


void printCommandLine(ProcessParamsTable *procparams, char *str)
{
	ProcessParamsTable *this=procparams;
	strcpy (str,"Command line: ");
	//strcat (str,this->program);
	while (this!=NULL) { 
		strcat (str," ");
		strcat (str,this->param);
		strcat (str," ");
		strcat (str,this->value);
		this=this->next;
	}
}


/* ============ Likelihood computations: ========== */

REAL8 ZeroLogLikelihood(LALVariables *currentParams, LALIFOData *data, LALTemplateFunction *template) {
  (void) currentParams; /* avoid warning about unused parameter */
  (void) data; /* avoid warning about unused parameter */
  (void) template; /* avoid warning about unused parameter */
  return 0.0;
}

REAL8 UndecomposedFreqDomainLogLikelihood(LALVariables *currentParams, LALIFOData * data, 
                              LALTemplateFunction *template)
/***************************************************************/
/* (log-) likelihood function.                                 */
/* Returns the non-normalised logarithmic likelihood.          */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Required (`currentParams') parameters are:                  */
/*   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       */
/*   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  */
/*   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        */
/*   - "distance"        (REAL8, Mpc, >0)                      */
/*   - "time"            (REAL8, GPS sec.)                     */
/***************************************************************/
{
  static int timeDomainWarning = 0;
  double Fplus, Fcross;
  double FplusScaled, FcrossScaled;
  double diffRe, diffIm, diffSquared;
  double dataReal, dataImag;
  REAL8 loglikeli;
  REAL8 plainTemplateReal, plainTemplateImag;
  REAL8 templateReal, templateImag;
  int i, lower, upper;
  LALIFOData *dataPtr;
  double ra, dec, psi, distMpc, gmst;
  double GPSdouble;
  LIGOTimeGPS GPSlal;
  double chisquared;
  double timedelay;  /* time delay b/w iterferometer & geocenter w.r.t. sky location */
  double timeshift;  /* time shift (not necessarily same as above)                   */
  double deltaT, TwoDeltaToverN, deltaF, twopit, f, re, im;
  double timeTmp;
  int different;
  LALStatus status;
  memset(&status,0,sizeof(status));
  LALVariables intrinsicParams;

  /* determine source's sky location & orientation parameters: */
  ra        = *(REAL8*) getVariable(currentParams, "rightascension"); /* radian      */
  dec       = *(REAL8*) getVariable(currentParams, "declination");    /* radian      */
  psi       = *(REAL8*) getVariable(currentParams, "polarisation");   /* radian      */
  GPSdouble = *(REAL8*) getVariable(currentParams, "time");           /* GPS seconds */
  distMpc   = *(REAL8*) getVariable(currentParams, "distance");       /* Mpc         */

  /* figure out GMST: */
  //XLALINT8NSToGPS(&GPSlal, floor(1e9 * GPSdouble + 0.5));
  XLALGPSSetREAL8(&GPSlal, GPSdouble);
  //UandA.units    = MST_RAD;
  //UandA.accuracy = LALLEAPSEC_LOOSE;
  //LALGPStoGMST1(&status, &gmst, &GPSlal, &UandA);
  gmst=XLALGreenwichMeanSiderealTime(&GPSlal);
  intrinsicParams.head      = NULL;
  intrinsicParams.dimension = 0;
  copyVariables(currentParams, &intrinsicParams);
  removeVariable(&intrinsicParams, "rightascension");
  removeVariable(&intrinsicParams, "declination");
  removeVariable(&intrinsicParams, "polarisation");
  removeVariable(&intrinsicParams, "time");
  removeVariable(&intrinsicParams, "distance");
  // TODO: add pointer to template function here.
  // (otherwise same parameters but different template will lead to no re-computation!!)

  chisquared = 0.0;
  /* loop over data (different interferometers): */
  dataPtr = data;

  while (dataPtr != NULL) {
    /* The parameters the Likelihood function can handle by itself   */
    /* (and which shouldn't affect the template function) are        */
    /* sky location (ra, dec), polarisation and signal arrival time. */
    /* Note that the template function shifts the waveform to so that*/
	/* t_c corresponds to the "time" parameter in                    */
	/* IFOdata->modelParams (set, e.g., from the trigger value).     */
    
    /* Reset log-likelihood */
    dataPtr->loglikelihood = 0.0;

    /* Compare parameter values with parameter values corresponding  */
    /* to currently stored template; ignore "time" variable:         */
    if (checkVariable(dataPtr->modelParams, "time")) {
      timeTmp = *(REAL8 *) getVariable(dataPtr->modelParams, "time");
      removeVariable(dataPtr->modelParams, "time");
    }
    else timeTmp = GPSdouble;
    different = compareVariables(dataPtr->modelParams, &intrinsicParams);
    /* "different" now may also mean that "dataPtr->modelParams" */
    /* wasn't allocated yet (as in the very 1st iteration).      */

    if (different) { /* template needs to be re-computed: */
      copyVariables(&intrinsicParams, dataPtr->modelParams);
      addVariable(dataPtr->modelParams, "time", &timeTmp, REAL8_t,PARAM_LINEAR);
      template(dataPtr);
      if (dataPtr->modelDomain == timeDomain) {
	if (!timeDomainWarning) {
	  timeDomainWarning = 1;
	  fprintf(stderr, "WARNING: using time domain template with frequency domain likelihood (in %s, line %d)\n", __FILE__, __LINE__);
	}
        executeFT(dataPtr);
        /* note that the dataPtr->modelParams "time" element may have changed here!! */
        /* (during "template()" computation)  */
      }
    }
    else { /* no re-computation necessary. Return back "time" value, do nothing else: */
      addVariable(dataPtr->modelParams, "time", &timeTmp, REAL8_t,PARAM_LINEAR);
    }

    /*-- Template is now in dataPtr->freqModelhPlus and dataPtr->freqModelhCross. --*/
    /*-- (Either freshly computed or inherited.)                            --*/

    /* determine beam pattern response (F_plus and F_cross) for given Ifo: */
    XLALComputeDetAMResponse(&Fplus, &Fcross,
                             dataPtr->detector->response,
			     ra, dec, psi, gmst);
    /* signal arrival time (relative to geocenter); */
    timedelay = XLALTimeDelayFromEarthCenter(dataPtr->detector->location,
                                             ra, dec, &GPSlal);
    /* (negative timedelay means signal arrives earlier at Ifo than at geocenter, etc.) */
    /* amount by which to time-shift template (not necessarily same as above "timedelay"): */
    timeshift =  (GPSdouble - (*(REAL8*) getVariable(dataPtr->modelParams, "time"))) + timedelay;
    twopit    = LAL_TWOPI * timeshift;

    /* include distance (overall amplitude) effect in Fplus/Fcross: */
    FplusScaled  = Fplus  / distMpc;
    FcrossScaled = Fcross / distMpc;

    if (checkVariable(currentParams, "crazyInjectionHLSign") &&
        *((INT4 *)getVariable(currentParams, "crazyInjectionHLSign"))) {
      if (strstr(dataPtr->name, "H") || strstr(dataPtr->name, "L")) {
        FplusScaled *= -1.0;
        FcrossScaled *= -1.0;
      }
    }

    dataPtr->fPlus = FplusScaled;
    dataPtr->fCross = FcrossScaled;
    dataPtr->timeshift = timeshift;

 //FILE *testout=fopen("test_likeliLAL.txt","w");
 //fprintf(testout, "f PSD dataRe dataIm signalRe signalIm\n");
    /* determine frequency range & loop over frequency bins: */
    deltaT = dataPtr->timeData->deltaT;
    deltaF = 1.0 / (((double)dataPtr->timeData->data->length) * deltaT);
    // printf("deltaF %g, Nt %d, deltaT %g\n", deltaF, dataPtr->timeData->data->length, dataPtr->timeData->deltaT);
    lower = ceil(dataPtr->fLow / deltaF);
    upper = floor(dataPtr->fHigh / deltaF);
    TwoDeltaToverN = 2.0 * deltaT / ((double) dataPtr->timeData->data->length);
    for (i=lower; i<=upper; ++i){
      /* derive template (involving location/orientation parameters) from given plus/cross waveforms: */
      plainTemplateReal = FplusScaled * dataPtr->freqModelhPlus->data->data[i].re  
                          +  FcrossScaled * dataPtr->freqModelhCross->data->data[i].re;
      plainTemplateImag = FplusScaled * dataPtr->freqModelhPlus->data->data[i].im  
                          +  FcrossScaled * dataPtr->freqModelhCross->data->data[i].im;

      /* do time-shifting...             */
      /* (also un-do 1/deltaT scaling): */
      f = ((double) i) * deltaF;
      /* real & imag parts of  exp(-2*pi*i*f*deltaT): */
      re = cos(twopit * f);
      im = - sin(twopit * f);
      templateReal = (plainTemplateReal*re - plainTemplateImag*im) / deltaT;
      templateImag = (plainTemplateReal*im + plainTemplateImag*re) / deltaT;
      dataReal     = dataPtr->freqData->data->data[i].re / deltaT;
      dataImag     = dataPtr->freqData->data->data[i].im / deltaT;
      /* compute squared difference & 'chi-squared': */
      diffRe       = dataReal - templateReal;         // Difference in real parts...
      diffIm       = dataImag - templateImag;         // ...and imaginary parts, and...
      diffSquared  = diffRe*diffRe + diffIm*diffIm ;  // ...squared difference of the 2 complex figures.
      REAL8 temp = ((TwoDeltaToverN * diffSquared) / dataPtr->oneSidedNoisePowerSpectrum->data->data[i]);
      chisquared  += temp;
      dataPtr->loglikelihood -= temp;
 //fprintf(testout, "%e %e %e %e %e %e\n",
 //        f, dataPtr->oneSidedNoisePowerSpectrum->data->data[i], 
 //        dataPtr->freqData->data->data[i].re, dataPtr->freqData->data->data[i].im,
 //        templateReal, templateImag);
    }
    dataPtr = dataPtr->next;
 //fclose(testout);
  }
  loglikeli = -1.0 * chisquared; // note (again): the log-likelihood is unnormalised!
  destroyVariables(&intrinsicParams);
  return(loglikeli);
}

REAL8 FreqDomainLogLikelihood(LALVariables *currentParams, LALIFOData * data, 
                              LALTemplateFunction *template)
/***************************************************************/
/* (log-) likelihood function.                                 */
/* Returns the non-normalised logarithmic likelihood.          */
/* Slightly slower but cleaner than							   */
/* UndecomposedFreqDomainLogLikelihood().          `		   */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Required (`currentParams') parameters are:                  */
/*   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       */
/*   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  */
/*   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        */
/*   - "distance"        (REAL8, Mpc, >0)                      */
/*   - "time"            (REAL8, GPS sec.)                     */
/***************************************************************/
{
  REAL8 loglikeli, totalChiSquared=0.0;
  LALIFOData *ifoPtr=data;
  COMPLEX16Vector *freqModelResponse=NULL;

  /* loop over data (different interferometers): */
  while (ifoPtr != NULL) {
    ifoPtr->loglikelihood = 0.0;

	if(freqModelResponse==NULL)
		freqModelResponse= XLALCreateCOMPLEX16Vector(ifoPtr->freqData->data->length);
	else
		freqModelResponse= XLALResizeCOMPLEX16Vector(freqModelResponse, ifoPtr->freqData->data->length);
	/*compute the response*/
	ComputeFreqDomainResponse(currentParams, ifoPtr, template, freqModelResponse);
	/*if(residual==NULL)
		residual=XLALCreateCOMPLEX16Vector(ifoPtr->freqData->data->length);
	else
		residual=XLALResizeCOMPLEX16Vector(residual, ifoPtr->freqData->data->length);
	
	COMPLEX16VectorSubtract(residual, ifoPtr->freqData->data, freqModelResponse);
	totalChiSquared+=ComputeFrequencyDomainOverlap(ifoPtr, residual, residual); 
	*/
        REAL8 temp = ComputeFrequencyDomainOverlap(ifoPtr, ifoPtr->freqData->data, ifoPtr->freqData->data)
          -2.0*ComputeFrequencyDomainOverlap(ifoPtr, ifoPtr->freqData->data, freqModelResponse)
          +ComputeFrequencyDomainOverlap(ifoPtr, freqModelResponse, freqModelResponse);
	totalChiSquared+=temp;
        ifoPtr->loglikelihood -= 0.5*temp;

    ifoPtr = ifoPtr->next;
  }
  loglikeli = -0.5 * totalChiSquared; // note (again): the log-likelihood is unnormalised!
  XLALDestroyCOMPLEX16Vector(freqModelResponse);
  return(loglikeli);
}

REAL8 ChiSquareTest(LALVariables *currentParams, LALIFOData * data, LALTemplateFunction *template)
/***************************************************************/
/* Chi-Square function.                                        */
/* Returns the chi square of a template:                       */
/* chisq= p * sum_i (dx_i)^2, with dx_i  =  <s,h>_i  - <s,h>/p */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Required (`currentParams') parameters are:                  */
/*   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       */
/*   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  */
/*   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        */
/*   - "distance"        (REAL8, Mpc, >0)                      */
/*   - "time"            (REAL8, GPS sec.)                     */
/***************************************************************/
{
  REAL8 ChiSquared=0.0, dxp, xp, x, norm, binPower, nextBin;
  REAL8 lowerF, upperF, deltaT, deltaF;
  REAL8 *segnorm;
  INT4  i, chisqPt, imax, kmin, kmax, numBins=0;
  INT4  *chisqBin;
  LALIFOData *ifoPtr=data;
  COMPLEX16Vector *freqModelResponse=NULL;
  //COMPLEX16Vector *segmentFreqModelResponse-NULL;

  /* Allocate memory for local pointers */
  segnorm=malloc(sizeof(REAL8) * ifoPtr->freqData->data->length);
  chisqBin=malloc(sizeof(INT4) * (numBins + 1));

  /* loop over data (different interferometers): */
  while (ifoPtr != NULL) {
    if(freqModelResponse==NULL)
      freqModelResponse= XLALCreateCOMPLEX16Vector(ifoPtr->freqData->data->length);
    else
      freqModelResponse= XLALResizeCOMPLEX16Vector(freqModelResponse, ifoPtr->freqData->data->length);
    /*compute the response*/
    ComputeFreqDomainResponse(currentParams, ifoPtr, template, freqModelResponse);

    deltaT = ifoPtr->timeData->deltaT;
    deltaF = 1.0 / (((REAL8)ifoPtr->timeData->data->length) * deltaT);
    
    /* Store values of fLow and fHigh to use later */
    lowerF = ifoPtr->fLow;
    upperF = ifoPtr->fHigh;
   
    /* Generate bin boundaries */
    numBins = *(INT4*) getVariable(currentParams, "numbins");
    kmin = ceil(ifoPtr->fLow / deltaF);
    kmax = floor(ifoPtr->fHigh / deltaF);
    imax = kmax > (INT4) ifoPtr->freqData->data->length-1 ? (INT4) ifoPtr->freqData->data->length-1 : kmax;
    
    memset(segnorm,0,sizeof(REAL8) * ifoPtr->freqData->data->length);
    norm = 0.0;
    
    for (i=1; i < imax; ++i){  	  	  
      norm += ((4.0 * deltaF * (freqModelResponse->data[i].re*freqModelResponse->data[i].re
              +freqModelResponse->data[i].im*freqModelResponse->data[i].im)) 
              / ifoPtr->oneSidedNoisePowerSpectrum->data->data[i]);
      segnorm[i] = norm;
    }


    memset(chisqBin,0,sizeof(INT4) * (numBins +1));

    binPower = norm / (REAL8) numBins;
    nextBin   = binPower;
    chisqPt   = 0;
    chisqBin[chisqPt++] = 0;

    for ( i = 1; i < imax; ++i )
    {
      if ( segnorm[i] >= nextBin )
      {
        chisqBin[chisqPt++] = i;
        nextBin += binPower;
        if ( chisqPt == numBins ) break;
      }
    }
    chisqBin[16]=imax;
    /* check that we have sucessfully allocated all the bins */
    if ( i == (INT4) ifoPtr->freqData->data->length && chisqPt != numBins )
    {
      /* if we have reaced the end of the template power vec and not
       * */
      /* allocated all the bin boundaries then there is a problem
       * */
      fprintf(stderr,"Error constructing frequency bins\n"); 
    }

    /* the last bin boundary is at can be at Nyquist since   */
    /* qtilde is zero above the ISCO of the current template */
    // chisqBin[numBins] = ifoPtr->freqData->data->length;

    /* end */
    
    x = ComputeFrequencyDomainOverlap(ifoPtr, ifoPtr->freqData->data, freqModelResponse)/(sqrt(norm));
    
    ChiSquared=0.0;
    
    for (i=0; i < numBins; ++i){
      
      ifoPtr->fLow = chisqBin[i] * deltaF;
      ifoPtr->fHigh = chisqBin[i+1] * deltaF;
      
      xp = ComputeFrequencyDomainOverlap(ifoPtr, ifoPtr->freqData->data, freqModelResponse)/(sqrt(norm));
      dxp = ((REAL8) numBins) * xp - x;
      ChiSquared += (dxp * dxp);
      
    }
    ChiSquared = ChiSquared / (REAL8) numBins;
    ifoPtr->fLow = lowerF;
    ifoPtr->fHigh = upperF;
    
    printf("Chi-Square for %s\t=\t%f\n",ifoPtr->detector->frDetector.name,ChiSquared);
    
    ifoPtr = ifoPtr->next;
  }
  free(chisqBin);
  free(segnorm);
  return(ChiSquared);
}

REAL8 TimeDomainLogLikelihood(LALVariables *currentParams, LALIFOData * data, 
                              LALTemplateFunction *template)
/***************************************************************/
/* (log-) likelihood function.                                 */
/* Returns the non-normalised logarithmic likelihood.          */
/* Slightly slower but cleaner than							   */
/* UndecomposedFreqDomainLogLikelihood().          `		   */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Required (`currentParams') parameters are:                  */
/*   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       */
/*   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  */
/*   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        */
/*   - "distance"        (REAL8, Mpc, >0)                      */
/*   - "time"            (REAL8, GPS sec.)                     */
/***************************************************************/
{
  REAL8 loglikeli, totalChiSquared=0.0;
  LALIFOData *ifoPtr=data;
  REAL8TimeSeries *timeModelResponse=NULL;

  /* loop over data (different interferometers): */
  while (ifoPtr != NULL) {
    ifoPtr->loglikelihood = 0.0;

    if(timeModelResponse==NULL) {
      timeModelResponse = 
	XLALCreateREAL8TimeSeries("time detector response", &(ifoPtr->timeData->epoch), 
				  0.0, ifoPtr->timeData->deltaT,
				  &lalDimensionlessUnit,
				  ifoPtr->timeData->data->length);
    } else if (timeModelResponse->data->length != ifoPtr->timeData->data->length) {
      /* Cannot resize *up* a time series, so just dealloc and reallocate it. */
      XLALDestroyREAL8TimeSeries(timeModelResponse);
      timeModelResponse =                   
	XLALCreateREAL8TimeSeries("time detector response", &(ifoPtr->timeData->epoch), 
				  0.0, ifoPtr->timeData->deltaT,
				  &lalDimensionlessUnit,
				  ifoPtr->timeData->data->length);
    }
     
    /*compute the response*/
    ComputeTimeDomainResponse(currentParams, ifoPtr, template, timeModelResponse);
    REAL8 temp = (WhitenedTimeDomainOverlap(ifoPtr->whiteTimeData, ifoPtr->windowedTimeData)
                  -2.0*WhitenedTimeDomainOverlap(ifoPtr->whiteTimeData, timeModelResponse)
                  +timeDomainOverlap(ifoPtr->timeDomainNoiseWeights, timeModelResponse, timeModelResponse));
    totalChiSquared+=temp;
    ifoPtr->loglikelihood -= 0.5*temp;
    
    ifoPtr = ifoPtr->next;
  }
  loglikeli = -0.5*totalChiSquared; 
  XLALDestroyREAL8TimeSeries(timeModelResponse);
  return(loglikeli);
}

void ComputeFreqDomainResponse(LALVariables *currentParams, LALIFOData * dataPtr, 
                              LALTemplateFunction *template, COMPLEX16Vector *freqWaveform)
/***************************************************************/
/* Frequency-domain single-IFO response computation.           */
/* Computes response for a given template.                     */
/* Will re-compute template only if necessary                  */
/* (i.e., if previous, as stored in data->freqModelhCross,     */
/* was based on different parameters or template function).    */
/* Carries out timeshifting for a given detector               */
/* and projection onto this detector.                          */
/* Result stored in freqResponse, assumed to be correctly      */
/* initialized												   */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Required (`currentParams') parameters are:                  */
/*   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       */
/*   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  */
/*   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        */
/*   - "distance"        (REAL8, Mpc, >0)                      */
/*   - "time"            (REAL8, GPS sec.)                     */
/***************************************************************/							  
{
  static int timeDomainWarning = 0;

	double ra, dec, psi, distMpc, gmst;
	
	double GPSdouble;
	double timeTmp;
	LIGOTimeGPS GPSlal;
	double timedelay;  /* time delay b/w iterferometer & geocenter w.r.t. sky location */
	double timeshift;  /* time shift (not necessarily same as above)                   */
	double deltaT, deltaF, twopit, f, re, im;

	int different;
	LALVariables intrinsicParams;
	LALStatus status;
	memset(&status,0,sizeof(status));

	double Fplus, Fcross;
	double FplusScaled, FcrossScaled;
	REAL8 plainTemplateReal, plainTemplateImag;
	UINT4 i;
	REAL8 mc;
	
	/* Fill in derived parameters if necessary */
	if(checkVariable(currentParams,"logdistance")){
		distMpc=exp(*(REAL8 *) getVariable(currentParams,"logdistance"));
		addVariable(currentParams,"distance",&distMpc,REAL8_t,PARAM_OUTPUT);
	}

	if(checkVariable(currentParams,"logmc")){
		mc=exp(*(REAL8 *)getVariable(currentParams,"logmc"));
		addVariable(currentParams,"chirpmass",&mc,REAL8_t,PARAM_OUTPUT);
	}
		
	
	/* determine source's sky location & orientation parameters: */
	ra        = *(REAL8*) getVariable(currentParams, "rightascension"); /* radian      */
	dec       = *(REAL8*) getVariable(currentParams, "declination");    /* radian      */
	psi       = *(REAL8*) getVariable(currentParams, "polarisation");   /* radian      */
	GPSdouble = *(REAL8*) getVariable(currentParams, "time");           /* GPS seconds */
	distMpc   = *(REAL8*) getVariable(currentParams, "distance");       /* Mpc         */
		
	/* figure out GMST: */
	//XLALINT8NSToGPS(&GPSlal, floor(1e9 * GPSdouble + 0.5));
	XLALGPSSetREAL8(&GPSlal, GPSdouble);
	//UandA.units    = MST_RAD;
	//UandA.accuracy = LALLEAPSEC_LOOSE;
	//LALGPStoGMST1(&status, &gmst, &GPSlal, &UandA);
	gmst=XLALGreenwichMeanSiderealTime(&GPSlal);
	intrinsicParams.head      = NULL;
	intrinsicParams.dimension = 0;
	copyVariables(currentParams, &intrinsicParams);
	removeVariable(&intrinsicParams, "rightascension");
	removeVariable(&intrinsicParams, "declination");
	removeVariable(&intrinsicParams, "polarisation");
	removeVariable(&intrinsicParams, "time");
	removeVariable(&intrinsicParams, "distance");
	// TODO: add pointer to template function here.
	// (otherwise same parameters but different template will lead to no re-computation!!)
      
	/* The parameters the response function can handle by itself     */
    /* (and which shouldn't affect the template function) are        */
    /* sky location (ra, dec), polarisation and signal arrival time. */
    /* Note that the template function shifts the waveform to so that*/
	/* t_c corresponds to the "time" parameter in                    */
	/* IFOdata->modelParams (set, e.g., from the trigger value).     */
    
    /* Compare parameter values with parameter values corresponding  */
    /* to currently stored template; ignore "time" variable:         */
    if (checkVariable(dataPtr->modelParams, "time")) {
      timeTmp = *(REAL8 *) getVariable(dataPtr->modelParams, "time");
      removeVariable(dataPtr->modelParams, "time");
    }
    else timeTmp = GPSdouble;
    different = compareVariables(dataPtr->modelParams, &intrinsicParams);
    /* "different" now may also mean that "dataPtr->modelParams" */
    /* wasn't allocated yet (as in the very 1st iteration).      */

    if (different) { /* template needs to be re-computed: */
      copyVariables(&intrinsicParams, dataPtr->modelParams);
      addVariable(dataPtr->modelParams, "time", &timeTmp, REAL8_t,PARAM_LINEAR);
      template(dataPtr);
      if (dataPtr->modelDomain == timeDomain) {
	if (!timeDomainWarning) {
	  timeDomainWarning = 1;
	  fprintf(stderr, "WARNING: using time domain template with frequency domain likelihood (in %s, line %d)\n", __FILE__, __LINE__);
	}
	executeFT(dataPtr);
      /* note that the dataPtr->modelParams "time" element may have changed here!! */
      /* (during "template()" computation)                                      */
      }
    }
    else { /* no re-computation necessary. Return back "time" value, do nothing else: */
      addVariable(dataPtr->modelParams, "time", &timeTmp, REAL8_t,PARAM_LINEAR);
    }

    /*-- Template is now in dataPtr->freqModelhPlus and dataPtr->freqModelhCross. --*/
    /*-- (Either freshly computed or inherited.)                            --*/

    /* determine beam pattern response (F_plus and F_cross) for given Ifo: */
    XLALComputeDetAMResponse(&Fplus, &Fcross, dataPtr->detector->response,
			     ra, dec, psi, gmst);
		 
    /* signal arrival time (relative to geocenter); */
    timedelay = XLALTimeDelayFromEarthCenter(dataPtr->detector->location,
                                             ra, dec, &GPSlal);
    /* (negative timedelay means signal arrives earlier at Ifo than at geocenter, etc.) */

    /* amount by which to time-shift template (not necessarily same as above "timedelay"): */
    timeshift =  (GPSdouble - (*(REAL8*) getVariable(dataPtr->modelParams, "time"))) + timedelay;
    twopit    = LAL_TWOPI * timeshift;

    /* include distance (overall amplitude) effect in Fplus/Fcross: */
    FplusScaled  = Fplus  / distMpc;
    FcrossScaled = Fcross / distMpc;

    
    if (checkVariable(currentParams, "crazyInjectionHLSign") &&
        *((INT4 *)getVariable(currentParams, "crazyInjectionHLSign"))) {
      if (strstr(dataPtr->name, "H") || strstr(dataPtr->name, "L")) {
        FplusScaled *= -1.0;
        FcrossScaled *= -1.0;
      }
    }

	if(freqWaveform->length!=dataPtr->freqModelhPlus->data->length){
		printf("fW%d data%d\n", freqWaveform->length, dataPtr->freqModelhPlus->data->length);
		printf("Error!  Frequency data vector must be same length as original data!\n");
		exit(1);
	}
	
	deltaT = dataPtr->timeData->deltaT;
    deltaF = 1.0 / (((double)dataPtr->timeData->data->length) * deltaT);

#ifdef DEBUG
FILE* file=fopen("TempSignal.dat", "w");	
#endif
	for(i=0; i<freqWaveform->length; i++){
		/* derive template (involving location/orientation parameters) from given plus/cross waveforms: */
		plainTemplateReal = FplusScaled * dataPtr->freqModelhPlus->data->data[i].re  
                          +  FcrossScaled * dataPtr->freqModelhCross->data->data[i].re;
		plainTemplateImag = FplusScaled * dataPtr->freqModelhPlus->data->data[i].im  
                          +  FcrossScaled * dataPtr->freqModelhCross->data->data[i].im;

		/* do time-shifting...             */
		/* (also un-do 1/deltaT scaling): */
		f = ((double) i) * deltaF;
		/* real & imag parts of  exp(-2*pi*i*f*deltaT): */
		re = cos(twopit * f);
		im = - sin(twopit * f);

		freqWaveform->data[i].re= (plainTemplateReal*re - plainTemplateImag*im);
		freqWaveform->data[i].im= (plainTemplateReal*im + plainTemplateImag*re);		
#ifdef DEBUG
		fprintf(file, "%lg %lg \t %lg\n", f, freqWaveform->data[i].re, freqWaveform->data[i].im);
#endif
	}
#ifdef DEBUG
fclose(file);
#endif
	destroyVariables(&intrinsicParams);
}

void ComputeTimeDomainResponse(LALVariables *currentParams, LALIFOData * dataPtr, 
                               LALTemplateFunction *template, REAL8TimeSeries *timeWaveform)
/***************************************************************/
/* Based on ComputeFreqDomainResponse above.                   */
/* Time-domain single-IFO response computation.                */
/* Computes response for a given template.                     */
/* Will re-compute template only if necessary                  */
/* (i.e., if previous, as stored in data->timeModelhCross,     */
/* was based on different parameters or template function).    */
/* Carries out timeshifting for a given detector               */
/* and projection onto this detector.                          */
/* Result stored in timeResponse, assumed to be correctly      */
/* initialized												   */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Required (`currentParams') parameters are:                  */
/*   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       */
/*   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  */
/*   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        */
/*   - "distance"        (REAL8, Mpc, >0)                      */
/*   - "time"            (REAL8, GPS sec.)                     */
/***************************************************************/							  
{
  static int freqDomainWarning = 0;
	double ra, dec, psi, distMpc, gmst;
	
	double GPSdouble;
	double timeTmp;
	LIGOTimeGPS GPSlal;
	double timedelay;  /* time delay b/w iterferometer & geocenter w.r.t. sky location */
	double timeshift;  /* time shift (not necessarily same as above)                   */

	int different;
	LALVariables intrinsicParams;
	LALStatus status;
	memset(&status,0,sizeof(status));

	double Fplus, Fcross;
	double FplusScaled, FcrossScaled;
	UINT4 i;
	REAL8 mc;

	/* Fill in derived parameters if necessary */
	if(checkVariable(currentParams,"logdistance")){
		distMpc=exp(*(REAL8 *) getVariable(currentParams,"logdistance"));
		addVariable(currentParams,"distance",&distMpc,REAL8_t,PARAM_OUTPUT);
	}

	if(checkVariable(currentParams,"logmc")){
		mc=exp(*(REAL8 *)getVariable(currentParams,"logmc"));
		addVariable(currentParams,"chirpmass",&mc,REAL8_t,PARAM_OUTPUT);
	}
		
	
	/* determine source's sky location & orientation parameters: */
	ra        = *(REAL8*) getVariable(currentParams, "rightascension"); /* radian      */
	dec       = *(REAL8*) getVariable(currentParams, "declination");    /* radian      */
	psi       = *(REAL8*) getVariable(currentParams, "polarisation");   /* radian      */
	GPSdouble = *(REAL8*) getVariable(currentParams, "time");           /* GPS seconds */
	distMpc   = *(REAL8*) getVariable(currentParams, "distance");       /* Mpc         */
		
	/* figure out GMST: */
	XLALGPSSetREAL8(&GPSlal, GPSdouble);
	//UandA.units    = MST_RAD;
	//UandA.accuracy = LALLEAPSEC_LOOSE;
	//LALGPStoGMST1(&status, &gmst, &GPSlal, &UandA);
	gmst=XLALGreenwichMeanSiderealTime(&GPSlal);
	intrinsicParams.head      = NULL;
	intrinsicParams.dimension = 0;
	copyVariables(currentParams, &intrinsicParams);
	removeVariable(&intrinsicParams, "rightascension");
	removeVariable(&intrinsicParams, "declination");
	removeVariable(&intrinsicParams, "polarisation");
	removeVariable(&intrinsicParams, "time");
	removeVariable(&intrinsicParams, "distance");
	// TODO: add pointer to template function here.
	// (otherwise same parameters but different template will lead to no re-computation!!)
      
	/* The parameters the response function can handle by itself     */
    /* (and which shouldn't affect the template function) are        */
    /* sky location (ra, dec), polarisation and signal arrival time. */
    /* Note that the template function shifts the waveform to so that*/
	/* t_c corresponds to the "time" parameter in                    */
	/* IFOdata->modelParams (set, e.g., from the trigger value).     */
    
    /* Compare parameter values with parameter values corresponding  */
    /* to currently stored template; ignore "time" variable:         */
    if (checkVariable(dataPtr->modelParams, "time")) {
      timeTmp = *(REAL8 *) getVariable(dataPtr->modelParams, "time");
      removeVariable(dataPtr->modelParams, "time");
    }
    else timeTmp = GPSdouble;
    different = compareVariables(dataPtr->modelParams, &intrinsicParams);
    /* "different" now may also mean that "dataPtr->modelParams" */
    /* wasn't allocated yet (as in the very 1st iteration).      */

    if (different) { /* template needs to be re-computed: */
      copyVariables(&intrinsicParams, dataPtr->modelParams);
      addVariable(dataPtr->modelParams, "time", &timeTmp, REAL8_t,PARAM_LINEAR);
      template(dataPtr);
      if (dataPtr->modelDomain == frequencyDomain) {
	if (!freqDomainWarning) {
	  freqDomainWarning = 1;
	  fprintf(stderr, "WARNING: frequency domain template used with time domain calculation (in %s, line %d)\n", __FILE__, __LINE__);
	}
        executeInvFT(dataPtr);
	/* note that the dataPtr->modelParams "time" element may have changed here!! */
	/* (during "template()" computation)  */
      }
    }
    else { /* no re-computation necessary. Return back "time" value, do nothing else: */
      addVariable(dataPtr->modelParams, "time", &timeTmp, REAL8_t,PARAM_LINEAR);
    }

    /*-- Template is now in dataPtr->timeModelhPlus and dataPtr->timeModelhCross. --*/
    /*-- (Either freshly computed or inherited.)                            --*/

    /* determine beam pattern response (F_plus and F_cross) for given Ifo: */
    XLALComputeDetAMResponse(&Fplus, &Fcross, dataPtr->detector->response,
			     ra, dec, psi, gmst);
		 
    /* signal arrival time (relative to geocenter); */
    timedelay = XLALTimeDelayFromEarthCenter(dataPtr->detector->location,
                                             ra, dec, &GPSlal);
    /* (negative timedelay means signal arrives earlier at Ifo than at geocenter, etc.) */

    /* amount by which to time-shift template (not necessarily same as above "timedelay"): */
    timeshift =  (GPSdouble - (*(REAL8*) getVariable(dataPtr->modelParams, "time"))) + timedelay;

    /* include distance (overall amplitude) effect in Fplus/Fcross: */
    FplusScaled  = Fplus  / distMpc;
    FcrossScaled = Fcross / distMpc;

    if (checkVariable(currentParams, "crazyInjectionHLSign") &&
        *((INT4 *)getVariable(currentParams, "crazyInjectionHLSign"))) {
      if (strstr(dataPtr->name, "H") || strstr(dataPtr->name, "L")) {
        FplusScaled *= -1.0;
        FcrossScaled *= -1.0;
      }
    }

	if(timeWaveform->data->length!=dataPtr->timeModelhPlus->data->length){
		printf("fW%d data%d\n", timeWaveform->data->length, dataPtr->freqModelhPlus->data->length);
		printf("Error!  Time data vector must be same length as original data!\n");
		exit(1);
	}
	
	timeWaveform->deltaT = dataPtr->timeModelhPlus->deltaT;
        
        /* Shift to correct start time. */
        timeWaveform->epoch = dataPtr->timeModelhPlus->epoch;
        XLALGPSAdd(&(timeWaveform->epoch), timeshift);

#ifdef DEBUG
FILE* file=fopen("TempSignal.dat", "w");	
#endif
	for(i=0; i<timeWaveform->data->length; i++){
#ifdef DEBUG
          double t = timeWaveform->deltaT*i + XLALGPSGetREAL8(&(timeWaveform->epoch));
#endif
                timeWaveform->data->data[i] = 
                  FplusScaled*dataPtr->timeModelhPlus->data->data[i] + 
                  FcrossScaled*dataPtr->timeModelhCross->data->data[i];
#ifdef DEBUG
		fprintf(file, "%lg %lg\n", t, timeWaveform->data[i]);
#endif
	}
#ifdef DEBUG
fclose(file);
#endif
	destroyVariables(&intrinsicParams);
}	
							  						  
REAL8 ComputeFrequencyDomainOverlap(LALIFOData * dataPtr,
                                    COMPLEX16Vector * freqData1, 
                                    COMPLEX16Vector * freqData2)
{
  int lower, upper, i;
  double deltaT, deltaF;
  
  double overlap=0.0;
  
  /* determine frequency range & loop over frequency bins: */
  deltaT = dataPtr->timeData->deltaT;
  deltaF = 1.0 / (((double)dataPtr->timeData->data->length) * deltaT);
  lower = ceil(dataPtr->fLow / deltaF);
  upper = floor(dataPtr->fHigh / deltaF);
	
  for (i=lower; i<=upper; ++i){
    overlap  += ((4.0*deltaF*(freqData1->data[i].re*freqData2->data[i].re+freqData1->data[i].im*freqData2->data[i].im)) 
                 / dataPtr->oneSidedNoisePowerSpectrum->data->data[i]);
  }

  return overlap;
}

void dumptemplateFreqDomain(LALVariables *currentParams, LALIFOData * data, 
                            LALTemplateFunction *template, const char *filename)
/* de-bugging function writing (frequency-domain) template to a CSV file */
/* File contains real & imaginary parts of plus & cross components.      */
/* Template amplitude is scaled to 1Mpc distance.                        */
{
  FILE *outfile=NULL; 
  LALIFOData *dataPtr;
  double deltaT, deltaF, f;
  UINT4 i;

  copyVariables(currentParams, data->modelParams);
  dataPtr = data;
  while (dataPtr != NULL) { /* this loop actually does nothing (yet) here. */
    template(data);
    if (data->modelDomain == timeDomain)
      executeFT(data);

    outfile = fopen(filename, "w");
    /*fprintf(outfile, "f PSD dataRe dataIm signalPlusRe signalPlusIm signalCrossRe signalCrossIm\n");*/
    fprintf(outfile, "\"f\",\"PSD\",\"signalPlusRe\",\"signalPlusIm\",\"signalCrossRe\",\"signalCrossIm\"\n");
    deltaT = dataPtr->timeData->deltaT;
    deltaF = 1.0 / (((double)dataPtr->timeData->data->length) * deltaT);
    for (i=0; i<data->freqModelhPlus->data->length; ++i){
      f = ((double) i) * deltaF;
      fprintf(outfile, "%f,%e,%e,%e,%e,%e\n",
              f, data->oneSidedNoisePowerSpectrum->data->data[i],
              /*data->freqData->data->data[i].re, data->freqData->data->data[i].im,*/
              data->freqModelhPlus->data->data[i].re,
              data->freqModelhPlus->data->data[i].im,
              data->freqModelhCross->data->data[i].re,
              data->freqModelhCross->data->data[i].im);
    }
    fclose(outfile);
    dataPtr = NULL;
  }
  fprintf(stdout, " wrote (frequency-domain) template to CSV file \"%s\".\n", filename);
}


void dumptemplateTimeDomain(LALVariables *currentParams, LALIFOData * data, 
                            LALTemplateFunction *template, const char *filename)
/* de-bugging function writing (frequency-domain) template to a CSV file */
/* File contains real & imaginary parts of plus & cross components.      */
/* Template amplitude is scaled to 1Mpc distance.                        */
{
  FILE *outfile=NULL; 
  LALIFOData *dataPtr;
  double deltaT, deltaF, t, epoch;
  UINT4 i;

  copyVariables(currentParams, data->modelParams);
  dataPtr = data;
  while (dataPtr != NULL) { /* this loop actually does nothing (yet) here. */
    template(data);
    if (data->modelDomain == frequencyDomain)
      executeInvFT(data);

    outfile = fopen(filename, "w");
    fprintf(outfile, "\"t\",\"signalPlus\",\"signalCross\"\n");
    deltaT = dataPtr->timeData->deltaT;
    deltaF = 1.0 / (((double)dataPtr->timeData->data->length) * deltaT);
    epoch = XLALGPSGetREAL8(&data->timeData->epoch);
    for (i=0; i<data->timeModelhPlus->data->length; ++i){
      t =  epoch + ((double) i) * deltaT;
      fprintf(outfile, "%f,%e,%e\n",
              t,
              data->timeModelhPlus->data->data[i],
              data->timeModelhCross->data->data[i]);
    }
    fclose(outfile);
    dataPtr = NULL;
  }
  fprintf(stdout, " wrote (time-domain) template to CSV file \"%s\".\n", filename);
}



void executeFT(LALIFOData *IFOdata)
/* Execute (forward, time-to-freq) Fourier transform.         */
/* Contents of IFOdata->timeModelh... are windowed and FT'ed, */
/* results go into IFOdata->freqModelh...                     */
/*  CHECK: keep or drop normalisation step here ?!?  */
{
  UINT4 i;
  double norm;
  
  for(;IFOdata;IFOdata=IFOdata->next){
    /* h+ */
    if(!IFOdata->freqModelhPlus)
      IFOdata->freqModelhPlus=(COMPLEX16FrequencySeries *)XLALCreateCOMPLEX16FrequencySeries("freqData",&(IFOdata->timeData->epoch),0.0,IFOdata->freqData->deltaF,&lalDimensionlessUnit,IFOdata->freqData->data->length);
    
    XLALDDVectorMultiply(IFOdata->timeModelhPlus->data,IFOdata->timeModelhPlus->data,IFOdata->window->data);
    XLALREAL8TimeFreqFFT(IFOdata->freqModelhPlus,IFOdata->timeModelhPlus,IFOdata->timeToFreqFFTPlan);
    
    /* hx */
    if(!IFOdata->freqModelhCross)
      IFOdata->freqModelhCross=(COMPLEX16FrequencySeries *)XLALCreateCOMPLEX16FrequencySeries("freqData",&(IFOdata->timeData->epoch),0.0,IFOdata->freqData->deltaF,&lalDimensionlessUnit,IFOdata->freqData->data->length);
    
    XLALDDVectorMultiply(IFOdata->timeModelhCross->data,IFOdata->timeModelhCross->data,IFOdata->window->data);
    XLALREAL8TimeFreqFFT(IFOdata->freqModelhCross,IFOdata->timeModelhCross,IFOdata->timeToFreqFFTPlan);
    
    norm=sqrt(IFOdata->window->sumofsquares/IFOdata->window->data->length);
    
    for(i=0;i<IFOdata->freqModelhPlus->data->length;i++){
      IFOdata->freqModelhPlus->data->data[i].re*=norm;
      IFOdata->freqModelhPlus->data->data[i].im*=norm;
      IFOdata->freqModelhCross->data->data[i].re*=norm;
      IFOdata->freqModelhCross->data->data[i].im*=norm;
    }
  }
}



void executeInvFT(LALIFOData *IFOdata)
/* Execute inverse (freq-to-time) Fourier transform. */
/* Results go into 'IFOdata->timeModelh...'          */
{
  while (IFOdata != NULL) {
    if (IFOdata->freqToTimeFFTPlan==NULL) die(" ERROR in executeInvFT(): encountered unallocated 'freqToTimeFFTPlan'.\n");

    /*  h+ :  */
    if (IFOdata->timeModelhPlus==NULL) die(" ERROR in executeInvFT(): encountered unallocated 'timeModelhPlus'.\n");
    if (IFOdata->freqModelhPlus==NULL) die(" ERROR in executeInvFT(): encountered unallocated 'freqModelhPlus'.\n");
    
    XLALREAL8FreqTimeFFT(IFOdata->timeModelhPlus, IFOdata->freqModelhPlus, IFOdata->freqToTimeFFTPlan);

    if (xlalErrno) {
      fprintf(stderr, "XLAL Error: %s (in %s, line %d)\n",
              XLALErrorString(xlalErrno), __FILE__, __LINE__);
      exit(1);
    }
    
    /*  hx :  */
    if (IFOdata->timeModelhCross==NULL) die(" ERROR in executeInvFT(): encountered unallocated 'timeModelhCross'.\n");
    if (IFOdata->freqModelhCross==NULL) die(" ERROR in executeInvFT(): encountered unallocated 'freqModelhCross'.\n");
    
    XLALREAL8FreqTimeFFT(IFOdata->timeModelhCross, IFOdata->freqModelhCross, IFOdata->freqToTimeFFTPlan);

    if (xlalErrno) {
      fprintf(stderr, "XLAL Error: %s (in %s, line %d)\n",
              XLALErrorString(xlalErrno), __FILE__, __LINE__);
      exit(1);
    }
    
    IFOdata=IFOdata->next;
  }
}


/* Function to add the min and max values for the prior onto the priorArgs */
void addMinMaxPrior(LALVariables *priorArgs, const char *name, void *min, void *max, VariableType type){
  char minName[VARNAME_MAX];
  char maxName[VARNAME_MAX];
  
  sprintf(minName,"%s_min",name);
  sprintf(maxName,"%s_max",name);
  
  addVariable(priorArgs,minName,min,type,PARAM_FIXED);
  addVariable(priorArgs,maxName,max,type,PARAM_FIXED);    
  return;
}

/* Get the min and max values of the prior from the priorArgs list, given a name */
void getMinMaxPrior(LALVariables *priorArgs, const char *name, void *min, void *max)
{
		char minName[VARNAME_MAX];
		char maxName[VARNAME_MAX];
		
		sprintf(minName,"%s_min",name);
		sprintf(maxName,"%s_max",name);
    
		*(REAL8 *)min=*(REAL8 *)getVariable(priorArgs,minName);
		*(REAL8 *)max=*(REAL8 *)getVariable(priorArgs,maxName);
		return;
		
}


REAL8 NullLogLikelihood(LALIFOData *data)
/*Idential to FreqDomainNullLogLikelihood                        */
{
	REAL8 loglikeli, totalChiSquared=0.0;
	LALIFOData *ifoPtr=data;
	
	/* loop over data (different interferometers): */
	while (ifoPtr != NULL) {
          ifoPtr->nullloglikelihood = 0.0;
          REAL8 temp = ComputeFrequencyDomainOverlap(ifoPtr, ifoPtr->freqData->data, ifoPtr->freqData->data);
          totalChiSquared+=temp;
          ifoPtr->nullloglikelihood -= 0.5*temp;
		ifoPtr = ifoPtr->next;
	}
	loglikeli = -0.5 * totalChiSquared; // note (again): the log-likelihood is unnormalised!
	return(loglikeli);
}

REAL8 WhitenedTimeDomainOverlap(const REAL8TimeSeries *whitenedData, const REAL8TimeSeries *data) {
  return 2.0*integrateSeriesProduct(whitenedData, data);
}

REAL8 TimeDomainNullLogLikelihood(LALIFOData *data) {
  REAL8 logL = 0.0;
  LALIFOData *ifoPtr = data;
  /* UINT4 ifoIndex = 0; */
  /* UINT4 i; */
  /* char fileName[256]; */
  /* FILE *out; */
  
  while (ifoPtr != NULL) {
    ifoPtr->nullloglikelihood = 0.0;
    REAL8 temp = WhitenedTimeDomainOverlap(ifoPtr->whiteTimeData, ifoPtr->windowedTimeData);
    logL += temp;
    ifoPtr->nullloglikelihood -= 0.5*temp;
 
    ifoPtr = ifoPtr->next;
    /* ifoIndex++; */
  }

  logL *= -0.5;
  
  return logL;
}

/* The time-domain weight corresponding to the (two-sided) noise power
   spectrum S(f) is defined to be:

   s(tau) == \int_{-infty}^\infty df \frac{\exp(2 \pi i f tau)}{S(f)}

*/
void PSDToTDW(REAL8TimeSeries *TDW, const REAL8FrequencySeries *PSD, const REAL8FFTPlan *plan,
              const REAL8 fMin, const REAL8 fMax) {
  COMPLEX16FrequencySeries *CPSD = NULL;
  UINT4 i;
  UINT4 PSDLength = TDW->data->length/2 + 1;
  
  if (PSD->data->length != PSDLength) {
    fprintf(stderr, "PSDToTDW: lengths of PSD and TDW do not match (in %s, line %d)", 
            __FILE__, __LINE__);
    exit(1);
  }

  CPSD = 
    XLALCreateCOMPLEX16FrequencySeries(PSD->name, &(PSD->epoch), PSD->f0, PSD->deltaF, &(PSD->sampleUnits), PSD->data->length);

  for (i = 0; i < PSD->data->length; i++) {
    REAL8 f = PSD->f0 + i*PSD->deltaF;

    if (fMin <= f && f <= fMax) {
      CPSD->data->data[i].re = 1.0 / (2.0*PSD->data->data[i]);
      CPSD->data->data[i].im = 0.0;
    } else {
      CPSD->data->data[i].re = 0.0;
      CPSD->data->data[i].im = 0.0;
    }
  }

  XLALREAL8FreqTimeFFT(TDW, CPSD, plan);

  /* FILE *PSDf = fopen("PSD.dat", "w"); */
  /* for (i = 0; i < PSD->data->length; i++) { */
  /*   fprintf(PSDf, "%g %g\n", i*PSD->deltaF, PSD->data->data[i]); */
  /* } */
  /* fclose(PSDf); */

  /* FILE *TDWf = fopen("TDW.dat", "w"); */
  /* for (i = 0; i < TDW->data->length; i++) { */
  /*   fprintf(TDWf, "%g %g\n", i*TDW->deltaT, TDW->data->data[i]); */
  /* } */
  /* fclose(TDWf); */
}

UINT4 nextPowerOfTwo(const UINT4 n) {
  UINT4 np2 = 1;
  
  for (np2 = 1; np2 < n; np2 *= 2) ; /* Keep doubling until >= n */

  return np2;
}

void padREAL8Sequence(REAL8Sequence *padded, const REAL8Sequence *data) {
  if (padded->length < data->length) {
    fprintf(stderr, "padREAL8Sequence: padded sequence too short (in %s, line %d)", __FILE__, __LINE__);
    exit(1);
  }

  memset(padded->data, 0, padded->length*sizeof(padded->data[0]));
  memcpy(padded->data, data->data, data->length*sizeof(data->data[0]));
}

void padWrappedREAL8Sequence(REAL8Sequence *padded, const REAL8Sequence *data) {
  UINT4 i;
  UINT4 np = padded->length;
  UINT4 nd = data->length;

  if (np < nd) {
    fprintf(stderr, "padWrappedREAL8Sequence: padded sequence too short (in %s, line %d)", __FILE__, __LINE__);
    exit(1);
  }

  memset(padded->data, 0, np*sizeof(padded->data[0]));

  padded->data[0] = data->data[0];
  for (i = 1; i <= (nd-1)/2; i++) {
    padded->data[i] = data->data[i]; /* Positive times/frequencies. */
    padded->data[np-i] = data->data[nd-i]; /* Wrapped, negative times/frequencies. */
  }
  if (nd % 2 == 0) { /* If even, take care of singleton positive frequency. */
    padded->data[nd/2] = data->data[nd/2];
  }
}

UINT4 LIGOTimeGPSToNearestIndex(const LIGOTimeGPS *tm, const REAL8TimeSeries *series) {
  REAL8 dT = XLALGPSDiff(tm, &(series->epoch));

  return (UINT4) (round(dT/series->deltaT));
}

REAL8 integrateSeriesProduct(const REAL8TimeSeries *s1, const REAL8TimeSeries *s2) {
  LIGOTimeGPS start, stop;
  LIGOTimeGPS stopS1, stopS2;
  UINT4 i1, i2;
  REAL8 sum = 0.0;
  REAL8 t1Start, t2Start, t, tStop;

  /* Compute stop times. */
  stopS1 = s1->epoch;
  stopS2 = s2->epoch;
  XLALGPSAdd(&stopS1, (s1->data->length-1)*s1->deltaT);
  XLALGPSAdd(&stopS2, (s2->data->length-1)*s2->deltaT);

  /* The start time is the max of the two start times, the stop time
     is the min of the two stop times */
  start = (XLALGPSCmp(&(s1->epoch), &(s2->epoch)) <= 0 ? s2->epoch : s1->epoch); /* Start at max start time. */
  stop = (XLALGPSCmp(&stopS1, &stopS2) <= 0 ? stopS1 : stopS2); /* Stop at min end time. */

  t = 0;
  tStop = XLALGPSDiff(&stop, &start);
  t1Start = XLALGPSDiff(&(s1->epoch), &start);
  t2Start = XLALGPSDiff(&(s2->epoch), &start);
  i1 = LIGOTimeGPSToNearestIndex(&start, s1);
  i2 = LIGOTimeGPSToNearestIndex(&start, s2);
  
  REAL8 *data1 = s1->data->data;
  REAL8 *data2 = s2->data->data;

  do {
    REAL8 nextTime1, nextTime2, nextTime;
    REAL8 dt;
    nextTime1 = t1Start + (i1+0.5)*s1->deltaT;
    nextTime2 = t2Start + (i2+0.5)*s2->deltaT;

    /* Whichever series needs updating first gives us the next time. */
    nextTime = (nextTime1 <= nextTime2 ? nextTime1 : nextTime2);

    /* Ensure we don't go past the stop time. */
    nextTime = (tStop < nextTime ? tStop : nextTime);

    dt = nextTime - t;

    sum += dt*data1[i1]*data2[i2];

    if (nextTime1 == nextTime) {
      i1++;
    }
    if (nextTime2 == nextTime) {
      i2++;
    }

    t = nextTime;    
  } while (t < tStop);

  return sum;
}

void convolveTimeSeries(REAL8TimeSeries *conv, const REAL8TimeSeries *data, const REAL8TimeSeries *response) {
  UINT4 responseSpan = (response->data->length + 1)/2;
  UINT4 paddedLength = nextPowerOfTwo(data->data->length + responseSpan);
  REAL8FFTPlan *fwdPlan = XLALCreateForwardREAL8FFTPlan(paddedLength, 1); /* Actually measure---rely on FFTW to store the best plan for a given length. */
  REAL8FFTPlan *revPlan = XLALCreateReverseREAL8FFTPlan(paddedLength, 1); /* Same. */
  REAL8Sequence *paddedData, *paddedResponse, *paddedConv;
  COMPLEX16Sequence *dataFFT, *responseFFT;
  UINT4 i;

  if (data->deltaT != response->deltaT) {
    fprintf(stderr, "convolveTimeSeries: sample spacings differ: %g vs %g (in %s, line %d)\n", 
            data->deltaT, response->deltaT, __FILE__, __LINE__);
    exit(1);
  }

  if (conv->data->length < data->data->length) {
    fprintf(stderr, "convolveTimeSeries: output length smaller than input length (in %s, line %d)", __FILE__, __LINE__);
    exit(1);
  }

  paddedData = XLALCreateREAL8Sequence(paddedLength);
  paddedResponse = XLALCreateREAL8Sequence(paddedLength);
  paddedConv = XLALCreateREAL8Sequence(paddedLength);

  dataFFT = XLALCreateCOMPLEX16Sequence(paddedLength/2 + 1); /* Exploit R -> C symmetry. */
  responseFFT = XLALCreateCOMPLEX16Sequence(paddedLength/2 + 1);

  padREAL8Sequence(paddedData, data->data);
  padWrappedREAL8Sequence(paddedResponse, response->data);

  XLALREAL8ForwardFFT(dataFFT, paddedData, fwdPlan);
  XLALREAL8ForwardFFT(responseFFT, paddedResponse, fwdPlan);

  for (i = 0; i < paddedLength/2 + 1; i++) {
    /* Store product in dataFFT. */
    double dataRe, dataIm, resRe, resIm;
    dataRe = dataFFT->data[i].re;
    dataIm = dataFFT->data[i].im;
    resRe = responseFFT->data[i].re;
    resIm = responseFFT->data[i].im;

    dataFFT->data[i].re = dataRe*resRe - dataIm*resIm;
    dataFFT->data[i].im = dataRe*resIm + resRe*dataIm;
  }

  XLALREAL8ReverseFFT(paddedConv, dataFFT, revPlan);

  memset(conv->data->data, 0, conv->data->length*sizeof(conv->data->data[0]));
  for (i = 0; i < data->data->length; i++) {
    conv->data->data[i] = conv->deltaT*paddedConv->data[i]/paddedConv->length; /* Normalize */
  }

  strncpy(conv->name, "convolved", LALNameLength);
  conv->epoch = data->epoch;
  conv->deltaT = data->deltaT;
  conv->f0 = data->f0;
  conv->sampleUnits = data->sampleUnits;  

  XLALDestroyREAL8FFTPlan(fwdPlan);
  XLALDestroyREAL8FFTPlan(revPlan);

  XLALDestroyREAL8Sequence(paddedData);
  XLALDestroyREAL8Sequence(paddedResponse);
  XLALDestroyREAL8Sequence(paddedConv);

  XLALDestroyCOMPLEX16Sequence(dataFFT);
  XLALDestroyCOMPLEX16Sequence(responseFFT);
}

void wrappedTimeSeriesToLinearTimeSeries(REAL8TimeSeries *linear, const REAL8TimeSeries *wrapped) {
  UINT4 NNeg, NPos, N, i;

  if (linear->data->length != wrapped->data->length) {
    fprintf(stderr, "wrappedTimeSeriesToLinearTimeSeries: lengths differ (in %s, line %d)",
            __FILE__, __LINE__);
    exit(1);
  }

  N = linear->data->length;
  NNeg = (N-1)/2;
  NPos = N-NNeg-1; /* 1 for the zero component. */

  for (i = 0; i < NNeg; i++) {
    linear->data->data[i] = wrapped->data->data[N-i-1];
  }
  linear->data->data[NNeg] = wrapped->data->data[0];
  for (i = 1; i <= NPos; i++) {
    linear->data->data[NNeg+i] = wrapped->data->data[i];
  }

  linear->epoch = wrapped->epoch;
  linear->deltaT = wrapped->deltaT;
  linear->f0 = wrapped->f0;
  linear->sampleUnits = wrapped->sampleUnits;
  
  /* Adjust start time for linear to account for negative components. */
  XLALGPSAdd(&linear->epoch, -(NNeg*linear->deltaT));
}

void linearTimeSeriesToWrappedTimeSeries(REAL8TimeSeries *wrapped, const REAL8TimeSeries *linear) {
  UINT4 NNeg, NPos, N, i;

  if (wrapped->data->length != linear->data->length) {
    fprintf(stderr, "linearTimeSeriesToWrappedTimeSeries: lengths differ (in %s, line %d)",
            __FILE__, __LINE__);
    exit(1);
  }

  N = wrapped->data->length;
  NNeg = (N-1)/2;
  NPos = N-NNeg-1;

  wrapped->data->data[0] = linear->data->data[NNeg];
  for (i = 1; i <= NPos; i++) {
    wrapped->data->data[i] = linear->data->data[NNeg+i];
  }
  for (i = 0; i < NNeg; i++) {
    wrapped->data->data[N-i-1] = linear->data->data[i];
  }

  wrapped->epoch = linear->epoch;
  wrapped->deltaT = linear->deltaT;
  wrapped->f0 = linear->f0;
  wrapped->sampleUnits = linear->sampleUnits;
  
  /* Adjust start time. */
  XLALGPSAdd(&wrapped->epoch, NNeg*wrapped->deltaT);
}

REAL8 timeDomainOverlap(const REAL8TimeSeries *TDW, const REAL8TimeSeries *A, const REAL8TimeSeries *B) {
  REAL8TimeSeries *Bconv;
  REAL8 overlap;

  Bconv = XLALCreateREAL8TimeSeries(B->name, &(B->epoch), 0.0, B->deltaT, &(B->sampleUnits), B->data->length);

  convolveTimeSeries(Bconv, B, TDW);

  overlap = integrateSeriesProduct(A, Bconv);

  XLALDestroyREAL8TimeSeries(Bconv);

  return 4.0*overlap; /* This is the overlap definition. */
}

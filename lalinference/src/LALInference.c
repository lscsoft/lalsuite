/*
 *  LALInference.c:  Bayesian Followup functions
 *
 *  Copyright (C) 2009, 2012 Ilya Mandel, Vivien Raymond, Christian
 *  Roever, Marc van der Sluys, John Veitch, and Will M. Farr
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
#include <math.h>
#include <lal/LALInference.h>
#include <lal/Units.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/VectorOps.h>
#include <lal/Date.h>
#include <lal/XLALError.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#define COL_MAX 128
#define STR_MAX 2048

size_t LALInferenceTypeSize[] = {sizeof(INT4),
                                   sizeof(INT8),
                                   sizeof(UINT4),
                                   sizeof(REAL4),
                                   sizeof(REAL8),
                                   sizeof(COMPLEX8),
                                   sizeof(COMPLEX16),
                                   sizeof(gsl_matrix *),
                                   sizeof(REAL8Vector *),
                                   sizeof(UINT4Vector *),
                                   sizeof(CHAR *),
                                   sizeof(LALInferenceMCMCRunPhase *),
                                   sizeof(void *)
};


/* ============ Accessor functions for the Variable structure: ========== */

static char *colNameToParamName(const char *colName);
/* Helper functions for sanity check */
static INT4 checkREAL8TimeSeries(REAL8TimeSeries *series);
static INT4 checkREAL8FrequencySeries(REAL8FrequencySeries *series);
static INT4 checkCOMPLEX16FrequencySeries(COMPLEX16FrequencySeries *series);
static INT4 matrix_equal(gsl_matrix *a, gsl_matrix *b);

/* This replaces gsl_matrix_equal which is only available with gsl 1.15+ */
/* Return 1 if matrices are equal, 0 otherwise */
static INT4 matrix_equal(gsl_matrix *a, gsl_matrix *b)
{
    if(!a||!b) return 0;
    if(a->size1!=b->size1 || a->size2!=b->size2) return 0;
    UINT4 i,j;
    for(i=0;i<a->size1;i++)
        for(j=0;j<a->size2;j++)
            if(gsl_matrix_get(a,i,j)!=gsl_matrix_get(b,i,j))
                return 0;

    return 1;
}

LALInferenceVariableItem *LALInferenceGetItem(const LALInferenceVariables *vars,const char *name)
/* (this function is only to be used internally) */
/* Returns pointer to item for given item name.  */
{
  if(vars==NULL) return NULL;
  LALInferenceVariableItem *this = vars->head;
  while (this != NULL) {
    if (!strcmp(this->name,name)) break;
    else this = this->next;
  }
  return(this);
}


LALInferenceVariableItem *LALInferenceGetItemNr(LALInferenceVariables *vars, int idx)
/* (this function is only to be used internally)  */
/* Returns pointer to item for given item number. */
{
  int i=1;
  if (idx < i) {
    XLAL_ERROR_NULL(XLAL_EINVAL, "Requesting zero or negative idx entry.");
  }
  LALInferenceVariableItem *this=vars->head;
  while (this != NULL) {
    if (i == idx) break;
    else {
      this = this->next;
      ++i;
    }
  }
  return(this);
}

LALInferenceParamVaryType LALInferenceGetVariableVaryType(LALInferenceVariables *vars, const char *name)
{
  return (LALInferenceGetItem(vars,name)->vary);
}

void LALInferenceSetParamVaryType(LALInferenceVariables *vars, const char *name, LALInferenceParamVaryType vary)
{
  LALInferenceVariableItem *item = LALInferenceGetItem(vars,name);
  item->vary = vary;
  return;
}

void *LALInferenceGetVariable(const LALInferenceVariables * vars,const char * name)
/* Return the value of variable name from the vars structure by walking the list */
{
  LALInferenceVariableItem *item;
  item=LALInferenceGetItem(vars,name);
  if(!item) {
    XLAL_ERROR_NULL(XLAL_EFAILED, "Entry \"%s\" not found.", name);
  }
  return(item->value);
}


INT4 LALInferenceGetVariableDimension(LALInferenceVariables *vars)
{
  return(vars->dimension);
}


INT4 LALInferenceGetVariableDimensionNonFixed(LALInferenceVariables *vars)
{
  INT4 count=0;
  gsl_matrix *m=NULL;
  UINT4Vector *v=NULL;
  LALInferenceVariableItem *ptr = vars->head;
  if (ptr==NULL) return count;
  else {
    /* loop over entries: */
    while (ptr != NULL) {
      /* print name: */
      //TBL: LALInferenceGetVariableDimensionNonFixed had to be modified for noise-parameters, which are stored in a gsl_matrix
      if (LALInferenceCheckVariableNonFixed(vars,ptr->name))
      {
        //Generalize to allow for other data types
        if(ptr->type == LALINFERENCE_gslMatrix_t)
        {
          m = *((gsl_matrix **)ptr->value);
          count += (int)( (m->size1)*(m->size2) );
        }
        else if(ptr->type == LALINFERENCE_UINT4Vector_t)
        {
          v = *((UINT4Vector **)ptr->value);
          count += (int)( v->length );
        }
        else count++;
      }
      ptr = ptr->next;
    }
  }
  return count;
}

LALInferenceVariableType LALInferenceGetVariableType(const LALInferenceVariables *vars, const char *name)
{
  return LALInferenceGetItem(vars,name)->type;
}

INT4 LALInferenceGetVariableTypeByIndex(LALInferenceVariables *vars, int idx)
/* Returns type of the i-th entry, */
/* where  1 <= idx <= dimension. */
{
  LALInferenceVariableItem *item;
  if ((idx < 1) || (idx > vars->dimension)){
    XLAL_ERROR(XLAL_EINVAL, "idx = %d, but needs to be 1 <= idx <= dimension = %d.", idx, vars->dimension);
  }
  item = LALInferenceGetItemNr(vars, idx);
  return(item->type);
}


char *LALInferenceGetVariableName(LALInferenceVariables *vars, int idx)
/* Returns (pointer to) the name of the i-th entry, */
/* where  1 <= idx <= dimension.                  */
{
  LALInferenceVariableItem *item;
  if ((idx < 1) || (idx > vars->dimension)){
    XLAL_ERROR_NULL(XLAL_EINVAL, "idx = %d, but needs to be 1 <= idx <= dimension = %d.", idx, vars->dimension);
  }
  item = LALInferenceGetItemNr(vars, idx);
  return(item->name);
}


void LALInferenceSetVariable(LALInferenceVariables * vars, const char * name, void *value)
/* Set the value of variable name in the vars structure to value */
{
  LALInferenceVariableItem *item;
  item=LALInferenceGetItem(vars,name);
  if(!item) {
    XLAL_ERROR_VOID(XLAL_EINVAL, "Entry \"%s\" not found.", name);
  }
  if (item->vary==LALINFERENCE_PARAM_FIXED) return;
  memcpy(item->value,value,LALInferenceTypeSize[item->type]);
  return;
}



void LALInferenceAddVariable(LALInferenceVariables * vars, const char * name, void *value, LALInferenceVariableType type, LALInferenceParamVaryType vary)
/* Add the variable name with type type and value value to vars */
/* If variable already exists, it will over-write the current value if type compatible*/
{
  LALInferenceVariableItem *old=NULL;
  /* Check the name doesn't already exist */
  if(LALInferenceCheckVariable(vars,name)) {
    old=LALInferenceGetItem(vars,name);
    if(old->type != type)
    {
      XLAL_ERROR_VOID(XLAL_EINVAL, "Cannot re-add \"%s\" as previous definition has wrong type.", name);
    }
    LALInferenceSetVariable(vars,name,value);
    return;
  }

  if(!value) {
    XLAL_ERROR_VOID(XLAL_EFAULT, "Unable to access value through null pointer; trying to add \"%s\".", name);
  }

  LALInferenceVariableItem *new=XLALMalloc(sizeof(LALInferenceVariableItem));

  memset(new,0,sizeof(LALInferenceVariableItem));
  if(new) {
    new->value = (void *)XLALMalloc(LALInferenceTypeSize[type]);
  }
  if(new==NULL||new->value==NULL) {
    XLAL_ERROR_VOID(XLAL_ENOMEM, "Unable to allocate memory for list item.");
  }
  memcpy(new->name,name,VARNAME_MAX);
  new->type = type;
  new->vary = vary;
  memcpy(new->value,value,LALInferenceTypeSize[type]);
  new->next = vars->head;
  vars->head = new;
  vars->dimension++;
  return;
}

void LALInferenceRemoveVariable(LALInferenceVariables *vars,const char *name)
{
  LALInferenceVariableItem *this;
  if(!vars)
    XLAL_ERROR_VOID(XLAL_EFAULT);
  this=vars->head;
  LALInferenceVariableItem *parent=NULL;
  while(this){
    if(!strcmp(this->name,name)) break;
    else {parent=this; this=this->next;}
  }
  if(!this){
    XLAL_PRINT_WARNING("Entry \"%s\" not found.", name);
    return;
  }
  if(!parent) vars->head=this->next;
  else parent->next=this->next;
  XLALFree(this->value);
  this->value=NULL;
  XLALFree(this);
  this=NULL;
  vars->dimension--;
  return;
}

int LALInferenceCheckVariableNonFixed(LALInferenceVariables *vars, const char *name)
/* Checks for a writeable variable */
{
  LALInferenceParamVaryType type;
  if(!LALInferenceCheckVariable(vars,name)) return 0;
  type=LALInferenceGetVariableVaryType(vars,name);
  if(type==LALINFERENCE_PARAM_CIRCULAR||type==LALINFERENCE_PARAM_LINEAR) return 1;
  else return 0;

}

int LALInferenceCheckVariable(LALInferenceVariables *vars,const char *name)
/* Check for existance of name */
{
  if(LALInferenceGetItem(vars,name)) return 1;
  else return 0;
}

void LALInferenceClearVariables(LALInferenceVariables *vars)
/* Free all variables inside the linked list, leaving only the head struct */
{
  LALInferenceVariableItem *this,*next;
  if(!vars) return;
  this=vars->head;
  if(this) next=this->next;
  while(this){
    if(this->type==LALINFERENCE_gslMatrix_t) gsl_matrix_free(*(gsl_matrix **)this->value);
    if(this->type==LALINFERENCE_UINT4Vector_t) XLALDestroyUINT4Vector(*(UINT4Vector **)this->value);
    if(this->type==LALINFERENCE_REAL8Vector_t) XLALDestroyREAL8Vector(*(REAL8Vector **)this->value);
    XLALFree(this->value);
    XLALFree(this);
    this=next;
    if(this) next=this->next;
  }
  vars->head=NULL;
  vars->dimension=0;
  return;
}

void LALInferenceCopyVariables(LALInferenceVariables *origin, LALInferenceVariables *target)
/*  copy contents of "origin" over to "target"  */
{
  int dims = 0, i = 0;

  /* Check that the source and origin differ */
  if(origin==target) return;

  LALInferenceVariableItem *ptr;
  if(!origin)
  {
    XLAL_ERROR_VOID(XLAL_EFAULT, "Unable to access origin pointer.");
  }

  /* Make sure the structure is initialised */
  if(!target) XLAL_ERROR_VOID(XLAL_EFAULT, "Unable to copy to uninitialised LALInferenceVariables structure.");
  /* first dispose contents of "target" (if any): */
  LALInferenceClearVariables(target);

  /* get the number of elements in origin */
  dims = LALInferenceGetVariableDimension( origin );

  if ( !dims ){
    XLAL_ERROR_VOID(XLAL_EFAULT, "Origin variables has zero dimensions!");
  }

  /* then copy over elements of "origin" - due to how elements are added by
     LALInferenceAddVariable this has to be done in reverse order to preserve
     the ordering of "origin"  */
  for ( i = dims; i > 0; i-- ){
    ptr = LALInferenceGetItemNr(origin, i);

    if(!ptr)
    {
      XLAL_ERROR_VOID(XLAL_EFAULT, "Bad LALInferenceVariable structure found while trying to copy.");
    }
    else
    {
      if(!ptr->value || !ptr->name){
        XLAL_ERROR_VOID(XLAL_EFAULT, "Badly formed LALInferenceVariableItem structure!");
      }
      /* Deep copy matrix and vector types */
      switch (ptr->type)
      {
          case LALINFERENCE_gslMatrix_t:
          {
            gsl_matrix *old=*(gsl_matrix **)ptr->value;
            gsl_matrix *new=gsl_matrix_alloc(old->size1,old->size2);
            if(!new) XLAL_ERROR_VOID(XLAL_ENOMEM,"Unable to create %ix%i matrix\n",old->size1,old->size2);
            gsl_matrix_memcpy(new,old);
            LALInferenceAddVariable(target,ptr->name,(void *)&new,ptr->type,ptr->vary);
            break;
          }
          case LALINFERENCE_UINT4Vector_t:
          {
            UINT4Vector *old=*(UINT4Vector **)ptr->value;
            UINT4Vector *new=XLALCreateUINT4Vector(old->length);
            if(new) memcpy(new->data,old->data,new->length*sizeof(new->data[0]));
            else XLAL_ERROR_VOID(XLAL_ENOMEM,"Unable to copy vector!\n");
            LALInferenceAddVariable(target,ptr->name,(void *)&new,ptr->type,ptr->vary);
            break;
          }
          case LALINFERENCE_REAL8Vector_t:
          {
            REAL8Vector *old=*(REAL8Vector **)ptr->value;
            REAL8Vector *new=XLALCreateREAL8Vector(old->length);
            if(new) memcpy(new->data,old->data,new->length);
            else XLAL_ERROR_VOID(XLAL_ENOMEM,"Unable to copy vector!\n");
            LALInferenceAddVariable(target,ptr->name,(void *)&new,ptr->type,ptr->vary);
            break;
          }
          default:
          { /* Just memcpy */
            LALInferenceAddVariable(target, ptr->name, ptr->value, ptr->type,
                                    ptr->vary);
            break;
          }
      }
    }
  }

  return;
}

/** Prints a variable item to a string (must be pre-allocated!) */
void LALInferencePrintVariableItem(char *out, LALInferenceVariableItem *ptr)
{
  if(ptr==NULL) {
    XLAL_ERROR_VOID(XLAL_EFAULT, "Null LALInferenceVariableItem pointer.");
  }
  if(out==NULL) {
    XLAL_ERROR_VOID(XLAL_EFAULT, "Null output string pointer.");
  }
  switch (ptr->type) {
        case LALINFERENCE_INT4_t:
          sprintf(out, "%d", *(INT4 *) ptr->value);
          break;
        case LALINFERENCE_INT8_t:
          sprintf(out, "%" LAL_INT8_FORMAT, *(INT8 *) ptr->value);
          break;
        case LALINFERENCE_UINT4_t:
          sprintf(out, "%ud", *(UINT4 *) ptr->value);
          break;
        case LALINFERENCE_REAL4_t:
          sprintf(out, "%.15lf", *(REAL4 *) ptr->value);
          break;
        case LALINFERENCE_REAL8_t:
          sprintf(out, "%.15lf", *(REAL8 *) ptr->value);
          break;
        case LALINFERENCE_COMPLEX8_t:
          sprintf(out, "%e + i*%e",
                 (REAL4) crealf(*(COMPLEX8 *) ptr->value), (REAL4) cimagf(*(COMPLEX8 *) ptr->value));
          break;
        case LALINFERENCE_COMPLEX16_t:
          sprintf(out, "%e + i*%e",
                 (REAL8) creal(*(COMPLEX16 *) ptr->value), (REAL8) cimag(*(COMPLEX16 *) ptr->value));
          break;
        case LALINFERENCE_UINT4Vector_t:
          sprintf(out, "<can't print vector>");
        break;
        case LALINFERENCE_gslMatrix_t:
          sprintf(out, "<can't print matrix>");
          break;
        default:
          sprintf(out, "<can't print>");
      }
  return;
}

void LALInferencePrintVariables(LALInferenceVariables *var)
/**
 * output contents of a 'LALInferenceVariables' structure * /
 * / * (by now only prints names and types, but no values)
 */
{
  LALInferenceVariableItem *ptr = var->head;
  fprintf(stdout, "LALInferenceVariables:\n");
  if (ptr==NULL) fprintf(stdout, "  <empty>\n");
  else {
    /* loop over entries: */
    while (ptr != NULL) {
      /* print name: */
      fprintf(stdout, "  \"%s\"", ptr->name);
      /* print type: */
      fprintf(stdout, "  (type #%d, ", ((int) ptr->type));
      switch (ptr->type) {
        case LALINFERENCE_INT4_t:
          fprintf(stdout, "'INT4'");
          break;
        case LALINFERENCE_INT8_t:
          fprintf(stdout, "'INT8'");
          break;
        case LALINFERENCE_UINT4_t:
          fprintf(stdout, "'UINT4'");
          break;
        case LALINFERENCE_REAL4_t:
          fprintf(stdout, "'REAL4'");
          break;
        case LALINFERENCE_REAL8_t:
          fprintf(stdout, "'REAL8'");
          break;
        case LALINFERENCE_COMPLEX8_t:
          fprintf(stdout, "'COMPLEX8'");
          break;
        case LALINFERENCE_COMPLEX16_t:
          fprintf(stdout, "'COMPLEX16'");
          break;
        case LALINFERENCE_UINT4Vector_t:
          fprintf(stdout, "'UINT4Vector'");
          break;
        case LALINFERENCE_gslMatrix_t:
          fprintf(stdout, "'gslMatrix'");
          break;
        default:
          fprintf(stdout, "<unknown type>");
      }
      fprintf(stdout, ")  ");
      /* print value: */
      gsl_matrix *matrix = NULL;
      switch (ptr->type) {
        case LALINFERENCE_INT4_t:
          fprintf(stdout, "%d", *(INT4 *) ptr->value);
          break;
        case LALINFERENCE_INT8_t:
          fprintf(stdout, "%" LAL_INT8_FORMAT, *(INT8 *) ptr->value);
          break;
        case LALINFERENCE_UINT4_t:
          fprintf(stdout, "%u", *(UINT4 *) ptr->value);
          break;
        case LALINFERENCE_REAL4_t:
          fprintf(stdout, "%.15lf", *(REAL4 *) ptr->value);
          break;
        case LALINFERENCE_REAL8_t:
          fprintf(stdout, "%.15lf", *(REAL8 *) ptr->value);
          break;
        case LALINFERENCE_COMPLEX8_t:
          fprintf(stdout, "%e + i*%e",
                 (REAL4) crealf(*(COMPLEX8 *) ptr->value), (REAL4) cimagf(*(COMPLEX8 *) ptr->value));
          break;
        case LALINFERENCE_COMPLEX16_t:
          fprintf(stdout, "%e + i*%e",
                 (REAL8) creal(*(COMPLEX16 *) ptr->value), (REAL8) cimag(*(COMPLEX16 *) ptr->value));
          break;
        case LALINFERENCE_UINT4Vector_t:
          //fprintf(stdout,"%iD matrix", (int)((UINT4Vector **)ptr->value)->size);
          fprintf(stdout,"[");
          fprintf(stdout,"%i]",(int)(*(UINT4Vector **)ptr->value)->length);
          break;
        case LALINFERENCE_gslMatrix_t:
          fprintf(stdout,"[");
          matrix = *((gsl_matrix **)ptr->value);
          fprintf(stdout,"%ix%i]",(int)(matrix->size1),(int)(matrix->size2));
          /*
          for(i=0; i<(int)( matrix->size1 ); i++)
          {
            for(j=0;j<(int)( matrix->size2 );j++)
            {
              fprintf(stdout,"%.2g",gsl_matrix_get(matrix, i, j));
              if(j<(int)( matrix->size2 )-1)fprintf(stdout,",");
            }
            if(i<(int)( matrix->size1 )-1)fprintf(stdout,"; ");
          }
          fprintf(stdout,"]");
           */
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

void LALInferencePrintSample(FILE *fp,LALInferenceVariables *sample){
  UINT4 i;//,j;
  //gsl_matrix *m=NULL;
  UINT4Vector *v=NULL;
  if(sample==NULL) return;
  LALInferenceVariableItem *ptr=sample->head;
  if(fp==NULL) return;
  while(ptr!=NULL) {
    switch (ptr->type) {
      case LALINFERENCE_INT4_t:
        fprintf(fp, "%"LAL_INT4_FORMAT, *(INT4 *) ptr->value);
        break;
      case LALINFERENCE_INT8_t:
        fprintf(fp, "%"LAL_INT8_FORMAT , *(INT8 *) ptr->value);
        break;
      case LALINFERENCE_UINT4_t:
        fprintf(fp, "%"LAL_UINT4_FORMAT , *(UINT4 *) ptr->value);
        break;
      case LALINFERENCE_REAL4_t:
        fprintf(fp, "%9.20e", *(REAL4 *) ptr->value);
        break;
      case LALINFERENCE_REAL8_t:
        fprintf(fp, "%9.20le", *(REAL8 *) ptr->value);
        break;
      case LALINFERENCE_COMPLEX8_t:
        fprintf(fp, "%e + i*%e",
            (REAL4) crealf(*(COMPLEX8 *) ptr->value), (REAL4) cimagf(*(COMPLEX8 *) ptr->value));
        break;
      case LALINFERENCE_COMPLEX16_t:
        fprintf(fp, "%e + i*%e",
            (REAL8) creal(*(COMPLEX16 *) ptr->value), (REAL8) cimag(*(COMPLEX16 *) ptr->value));
        break;
      case LALINFERENCE_string_t:
        fprintf(fp, "%s", *((CHAR **)ptr->value));
        break;
      case LALINFERENCE_UINT4Vector_t:
        v=*((UINT4Vector **)ptr->value);
        for(i=0;i<v->length;i++) fprintf(fp,"%"LAL_UINT4_FORMAT" ",v->data[i]);
        break;
	  case LALINFERENCE_gslMatrix_t:
        /*
        m = *((gsl_matrix **)ptr->value);
        for(i=0; i<(int)( m->size1 ); i++)
        {
          for(j=0; j<(int)( m->size2); j++)
          {
            fprintf(fp,"%11.7f",gsl_matrix_get(m, i, j));
            if(i<(int)( m->size1 )-1 && j<(int)( m->size2)-1) fprintf(fp,"\t");
          }
        }
         */
        break;
      default:
        XLALPrintWarning("<can't print>");
      }

  fprintf(fp,"\t");
  ptr=ptr->next;
  }
  return;
}

void LALInferencePrintSampleNonFixed(FILE *fp,LALInferenceVariables *sample){
  UINT4 i;//,j;
  //gsl_matrix *m=NULL;
  UINT4Vector *v=NULL;
	if(sample==NULL) return;
	LALInferenceVariableItem *ptr=sample->head;
	if(fp==NULL) return;
	while(ptr!=NULL) {
		if (ptr->vary != LALINFERENCE_PARAM_FIXED && ptr->type != LALINFERENCE_gslMatrix_t ) {
			switch (ptr->type) {
				case LALINFERENCE_INT4_t:
					fprintf(fp, "%"LAL_INT4_FORMAT, *(INT4 *) ptr->value);
					break;
				case LALINFERENCE_INT8_t:
					fprintf(fp, "%"LAL_INT8_FORMAT, *(INT8 *) ptr->value);
					break;
				case LALINFERENCE_UINT4_t:
					fprintf(fp, "%"LAL_UINT4_FORMAT, *(UINT4 *) ptr->value);
					break;
				case LALINFERENCE_REAL4_t:
					fprintf(fp, "%9.20e", *(REAL4 *) ptr->value);
					break;
				case LALINFERENCE_REAL8_t:
					fprintf(fp, "%9.20le", *(REAL8 *) ptr->value);
					break;
				case LALINFERENCE_COMPLEX8_t:
					fprintf(fp, "%e + i*%e",
							(REAL4) crealf(*(COMPLEX8 *) ptr->value), (REAL4) cimagf(*(COMPLEX8 *) ptr->value));
					break;
				case LALINFERENCE_COMPLEX16_t:
					fprintf(fp, "%e + i*%e",
							(REAL8) creal(*(COMPLEX16 *) ptr->value), (REAL8) cimag(*(COMPLEX16 *) ptr->value));
					break;
        case LALINFERENCE_UINT4Vector_t:
          v = *((UINT4Vector **)ptr->value);
          for(i=0;i<v->length;i++)
          {
            fprintf(fp,"%11.7f",(REAL8)v->data[i]);
            if( i!=(UINT4)(v->length-1) )fprintf(fp,"\t");
          }
          break;
          /*
				case LALINFERENCE_gslMatrix_t:
                    m = *((gsl_matrix **)ptr->value);
                    for(i=0; i<(int)( m->size1 ); i++)
                    {
                        for(j=0; j<(int)( m->size2); j++)
                        {
                            fprintf(fp,"%11.7f",gsl_matrix_get(m, i, j));
                            if(i<(int)( m->size1 )-1 && j<(int)( m->size2)-1) fprintf(fp,"\t");
                        }
                    }
					break;
           */
				default:
					fprintf(stdout, "<can't print>");
			}
		fprintf(fp,"\t");
		}
		ptr=ptr->next;
	}
	return;
}

void LALInferenceReadSampleNonFixed(FILE *fp, LALInferenceVariables *p) {
  if (p == NULL || fp == NULL) return;
  LALInferenceVariableItem *item = p->head;
  while (item != NULL) {
    if (item->vary != LALINFERENCE_PARAM_FIXED) {
      switch (item->type) {
      case LALINFERENCE_INT4_t:
	fscanf(fp, "%"LAL_INT4_FORMAT, (INT4 *)item->value);
	break;
      case LALINFERENCE_INT8_t:
	fscanf(fp, "%"LAL_INT8_FORMAT, (INT8 *)item->value);
	break;
      case LALINFERENCE_UINT4_t:
	fscanf(fp, "%"LAL_UINT4_FORMAT, (UINT4 *)item->value);
	break;
      case LALINFERENCE_REAL4_t:
	fscanf(fp, "%"LAL_REAL4_FORMAT, (REAL4 *)item->value);
	break;
      case LALINFERENCE_REAL8_t:
	fscanf(fp, "%"LAL_REAL8_FORMAT, (REAL8 *)item->value);
	break;
      default:
	/* Pass on reading */
	XLAL_ERROR_VOID(XLAL_EINVAL, "cannot read data type into LALINferenceVariables");
	break;
      }
    }

    item = item->next;
  }
}

/**
 * Utility for readling in delimited ASCII files.
 *
 * Reads in an ASCII (delimited) file, and returns the results in a REAL8 array.
 * @param[in]  input       Input stream to be parsed.
 * @param[in]  nCols       Number of total columns in the input stream.
 * @param[in]  nWantedCols Number of columns to be saved.
 * @param[in]  wantedCols  Array of 0/1 flags (should be \a nCols long), indicating desired columns.
 * @param[out] nLines      Total number of lines read.
 * @return A REAL8 array containing the parsed data.
 */
REAL8 *LALInferenceParseDelimitedAscii(FILE *input, UINT4 nCols, UINT4 nWantedCols, UINT4 *wantedCols, UINT4 *nLines) {
    UINT4 nread;
    UINT4 i=0, par=0, col=0;
    REAL8 val=0;

    // Determine number of samples to be read
    unsigned long startPostBurnin = ftell(input);
    UINT4 nSamples=0;

    INT4 ch;
    while ( (ch = getc(input)) != EOF) {
        if (ch=='\n')
            ++nSamples;
    }
    fseek(input,startPostBurnin,SEEK_SET);

    // Read in samples
    REAL8 *sampleArray;
    sampleArray = (REAL8*) XLALMalloc(nSamples*nWantedCols*sizeof(REAL8));

    for (i = 0; i < nSamples; i++) {
        par=0;
        for (col = 0; col < nCols; col++) {
            nread = fscanf(input, "%lg", &val);
            if (nread != 1) {
                fprintf(stderr, "Cannot read sample from file (in %s, line %d)\n",
                __FILE__, __LINE__);
                exit(1);
            }

            if (wantedCols[col]) {
                sampleArray[i*nWantedCols + par] = val;
                par++;
            }
        }
    }

    *nLines = nSamples;
    return sampleArray;
}


/**
 * Parse a single line of delimited ASCII.
 *
 * Splits a line using \a delim into an array of strings.
 * @param[in]  record Line to be parsed.
 * @param[in]  delim  Delimiter to split on.
 * @param[out] arr    Parsed string.
 * @param[out] cnt    Total number of fields read.
 */
void parseLine(char *record, const char *delim, char arr[][VARNAME_MAX], UINT4 *cnt) {
    // Discard newline character at end of line to avoid including it in a field
    char *newline = strchr(record, '\n');
    if (newline)
        *newline=0;

    INT4 i=0;
    char *p = strtok(record, delim);
    while (p) {
        strcpy(arr[i], p);
        i++;
        p = strtok(NULL, delim);
    }
    *cnt = i;
}


/**
 * Discard the standard header of a PTMCMC chain file.
 *
 * Reads a PTMCMC input stream, moving foward until it finds a line beginning
 * with "cycle", which is typically the line containing the column names of a
 * PTMCMC output file.  The final stream points to the line containing the
 * column names.
 * @param filestream The PTMCMC input stream to discard the header of.
 */
void LALInferenceDiscardPTMCMCHeader(FILE *filestream) {
    char str[STR_MAX];
    char row[COL_MAX][VARNAME_MAX];
    const char *delimiters = " \t";
    UINT4 nCols;
    INT4 last_line;

    fgets(str, sizeof(str), filestream);
    parseLine(str, delimiters, row, &nCols);
    while (strcmp(row[0], "cycle") && str != NULL) {
        last_line = ftell(filestream);
        fgets(str, sizeof(str), filestream);
        parseLine(str, delimiters, row, &nCols);
    }

    if (str == NULL) {
        fprintf(stderr, "Couldn't find column headers in PTMCMC file.\n");
        exit(1);
    }

    fseek(filestream, last_line, SEEK_SET);
    return;
}


/**
 * Burn-in a PTMCMC output file.
 *
 * Read through a PTMCMC output file until the desired cycle is reached.
 * @param     filestream  The PTMCMC input stream to be burned in.
 * @param[in] burninCycle The cycle for \a filestream to point to.
 */
void LALInferenceBurninPTMCMC(FILE *filestream, UINT4 burninCycle) {
    char str[STR_MAX];
    char row[COL_MAX][VARNAME_MAX];
    const char *delimiters = " \t";
    UINT4 nCols;

    fgets(str, sizeof(str), filestream);
    parseLine(str, delimiters, row, &nCols);
    while ((UINT4)atoi(row[0]) <= burninCycle && str != NULL) {
        fgets(str, sizeof(str), filestream);
        parseLine(str, delimiters, row, &nCols);
    }

    if (str == NULL) {
        fprintf(stderr, "Error burning in PTMCMC file.\n");
        exit(1);
    }
}


/**
 * Burn-in a generic ASCII stream.
 *
 * Reads past the desired number of lines of a filestream.
 * @param     filestream The filestream to be burned in.
 * @param[in] burnin     Number of lines to read past in \a filestream.
 */
void LALInferenceBurninStream(FILE *filestream, UINT4 burnin) {
    char str[STR_MAX];
    UINT4 i=0;

    printf("Burnin: %i.\n", burnin);
    for (i=0; i<burnin; i++)
        fgets(str, sizeof(str), filestream);

    if (str == NULL && burnin > 0) {
        fprintf(stderr, "Error burning in file.\n");
        exit(1);
    }
}

/**
 * Read desired column names from an ASCII file.
 *
 * Reads the column names for an output file, and determines which columns are
 * parameter names by looking for standard column outputs that aren't parameter
 * names (e.g. cycle, logpost).  The internal list should be updated if other
 * non-parameter columns are added in the future.
 * @param[in]  input      The input files stream to parse.
 * @param[out] params     The column names of found to be valid parameters.
 * @param[out] nTotalCols Total number of columns in \a input.
 * @param[out] nValidCols Number of columns that are parameters.
 * @param[out] validCols  Array of 0/1 flags indicating which columns are parameters.
 */
void LALInferenceReadAsciiHeader(FILE *input, char params[][VARNAME_MAX], UINT4 *nTotalCols, UINT4 *nValidCols, UINT4 **validCols) {
    char str[STR_MAX];
    char row[COL_MAX][VARNAME_MAX];
    const char *delimiters = " \n\t";
    UINT4 nCols=0, nPar=0, par=0;
    UINT4 i, j;

    const char *non_params[] = {"cycle","timestamp","logpost","logprior","logl","loglH1","loglL1","loglV1","",NULL};

    fgets(str, sizeof(str), input);
    parseLine(str, delimiters, row, &nCols);

    UINT4 is_param[COL_MAX];
    for (i=0; i<nCols; i++) {
        nPar++;
        is_param[i] = 1;
        j=0;
        while (non_params[j] != NULL) {
            if (!strcmp(non_params[j], row[i])) {
                is_param[i] = 0;
                nPar--;
                break;
            }
            j++;
        }
    }

    for (i=0; i<nCols; i++) {
        if (is_param[i]) {
            strcpy(params[par], row[i]);
            par++;
        }
    }

    *nTotalCols = nCols;
    *validCols = is_param;
    *nValidCols = nPar;
}


/**
 * Utility for selecting columns from an array, in the specified order.
 *
 * Selects a subset of columns from an existing array and create a new array with them.
 * @param[in]  inarray  The array to select columns from.
 * @param[in]  nRows    Number of rows in \a inarray.
 * @param[in]  nCols    Number of columns in \a inarray.
 * @param[in]  nSelCols Number of columns being extracted.
 * @param[in]  selCols  Array of column numbers to be extracted from \a inarray.
 * @returns An array containing the requested columns in the order specified in \a selCols.
 */
REAL8 **LALInferenceSelectColsFromArray(REAL8 **inarray, UINT4 nRows, UINT4 nCols, UINT4 nSelCols, UINT4 *selCols) {
    UINT4 i,j;

    REAL8 **array = (REAL8**) XLALMalloc(nRows * sizeof(REAL8*));
    for (i = 0; i < nRows; i++) {
        array[i] = XLALMalloc(nSelCols * sizeof(REAL8));

        for (j = 0; j < nSelCols; j++) {
            if (selCols[j] > nCols) {
                XLAL_ERROR_NULL(XLAL_FAILURE, "Requesting a column number that is out of bounds.");
            }
            array[i][j] = inarray[i][selCols[j]];
        }
    }
    return array;
}

int LALInferencePrintProposalStatsHeader(FILE *fp,LALInferenceVariables *propStats) {
  LALInferenceVariableItem *head = propStats->head;
  while (head != NULL) {
    fprintf(fp, "%s\t", head->name);
    head = head->next;
  }
  fprintf(fp, "\n");
  return 0;
}

void LALInferencePrintProposalStats(FILE *fp,LALInferenceVariables *propStats){
  REAL4 accepted = 0;
  REAL4 proposed = 0;
  REAL4 acceptanceRate = 0;

  if(propStats==NULL || fp==NULL) return;
  LALInferenceVariableItem *ptr=propStats->head;
  while(ptr!=NULL) {
    accepted = (REAL4) (*(LALInferenceProposalStatistics *) ptr->value).accepted;
    proposed = (REAL4) (*(LALInferenceProposalStatistics *) ptr->value).proposed;
    acceptanceRate = accepted/(proposed==0 ? 1.0 : proposed);
    fprintf(fp, "%9.5f\t", acceptanceRate);
    ptr=ptr->next;
  }
  fprintf(fp, "\n");
  return;
}

const char *LALInferenceTranslateInternalToExternalParamName(const char *inName) {
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
  } else if (!strcmp(inName, "tilt_spin1")) {
    return "tilt1";
  } else if (!strcmp(inName, "tilt_spin2")) {
    return "tilt2";
  } else if (!strcmp(inName, "chirpmass")) {
    return "mc";
  } else if (!strcmp(inName, "massratio")) {
    return "eta";
  } else if (!strcmp(inName, "asym_massratio")) {
    return "q";
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
  } else if (!strcmp(inName, "fRef")) {
    return "f_ref";
  } else {
    return inName;
  }
}

void LALInferenceTranslateExternalToInternalParamName(char *outName, const char *inName) {
  if (!strcmp(inName, "a1")) {
    strcpy(outName, "a_spin1");
  } else if (!strcmp(inName, "a2")) {
    strcpy(outName, "a_spin2");
  } else if (!strcmp(inName, "phi1")) {
    strcpy(outName, "phi_spin1");
  } else if (!strcmp(inName, "phi2")) {
    strcpy(outName, "phi_spin2");
  } else if (!strcmp(inName, "theta1")) {
    strcpy(outName, "theta_spin1");
  } else if (!strcmp(inName, "theta2")) {
    strcpy(outName, "theta_spin2");
  } else if (!strcmp(inName, "tilt1")) {
    strcpy(outName, "tilt_spin1");
  } else if (!strcmp(inName, "tilt2")) {
    strcpy(outName, "tilt_spin2");
  } else if (!strcmp(inName, "mc")) {
    strcpy(outName, "chirpmass");
  } else if (!strcmp(inName, "eta")) {
    strcpy(outName, "massratio");
  } else if (!strcmp(inName, "q")) {
    strcpy(outName, "asym_massratio");
  } else if (!strcmp(inName, "ra")) {
    strcpy(outName, "rightascension");
  } else if (!strcmp(inName, "dec")) {
    strcpy(outName, "declination");
  } else if (!strcmp(inName, "phi_orb")) {
    strcpy(outName, "phase");
  } else if (!strcmp(inName, "psi")) {
    strcpy(outName, "polarisation");
  } else if (!strcmp(inName, "iota")) {
    strcpy(outName, "inclination");
  } else if (!strcmp(inName, "dist")) {
    strcpy(outName, "distance");
  } else if (!strcmp(inName, "f_ref")) {
    strcpy(outName, "fRef");
  } else {
    strcpy(outName, inName);
  }
}


int LALInferenceFprintParameterHeaders(FILE *out, LALInferenceVariables *params) {
  LALInferenceVariableItem *head = params->head;
  int i,j;
  gsl_matrix *matrix = NULL;
  UINT4Vector *vector = NULL;

  while (head != NULL) {
      if(head->type==LALINFERENCE_gslMatrix_t)
      {
          matrix = *((gsl_matrix **)head->value);
          for(i=0; i<(int)matrix->size1; i++)
          {
              for(j=0; j<(int)matrix->size2; j++)
              {
                  fprintf(out, "%s%i%i\t", LALInferenceTranslateInternalToExternalParamName(head->name),i,j);
              }
          }
      }
      else if(head->type==LALINFERENCE_UINT4Vector_t)
      {
        vector = *((UINT4Vector **)head->value);
        for(i=0; i<(int)vector->length; i++) fprintf(out, "%s%i\t", LALInferenceTranslateInternalToExternalParamName(head->name),i);
      }
      else fprintf(out, "%s\t", LALInferenceTranslateInternalToExternalParamName(head->name));
      head = head->next;
  }
  return 0;
}

int LALInferenceFprintParameterNonFixedHeaders(FILE *out, LALInferenceVariables *params) {
  LALInferenceVariableItem *head = params->head;

  int i;//,j;
  //gsl_matrix *matrix = NULL;
  UINT4Vector *vector = NULL;
  while (head != NULL) {
    if (head->vary != LALINFERENCE_PARAM_FIXED) {
      if(head->type==LALINFERENCE_gslMatrix_t)
      {
        /*
        fprintf(stdout,"\n");
        fprintf(stdout,"Skipping noise/glitch amplitudes in output files\n");
        fprintf(stdout,"   edit LALInferenceFprintParameterNonFixedHeaders()\n");
        fprintf(stdout,"   and LALInferencePrintSampleNonFixed() to modify\n");
         */
        /*
        matrix = *((gsl_matrix **)head->value);
        for(i=0; i<(int)matrix->size1; i++)
        {
          for(j=0; j<(int)matrix->size2; j++)
          {
            fprintf(out, "%s%i%i\t", LALInferenceTranslateInternalToExternalParamName(head->name),i,j);
          }
        }
        */
      }
      else if(head->type==LALINFERENCE_UINT4Vector_t)
      {
        vector = *((UINT4Vector **)head->value);
        for(i=0; i<(int)vector->length; i++) fprintf(out, "%s%i\t", LALInferenceTranslateInternalToExternalParamName(head->name),i);
      }
      else fprintf(out, "%s\t", LALInferenceTranslateInternalToExternalParamName(head->name));
    }
    head = head->next;
  }

  return 0;
}

INT4 LALInferenceFprintParameterNonFixedHeadersWithSuffix(FILE *out, LALInferenceVariables *params, const char *suffix) {
  LALInferenceVariableItem *head = params->head;

  INT4 i,j;
  gsl_matrix *matrix = NULL;

  while (head != NULL) {
    if (head->vary != LALINFERENCE_PARAM_FIXED) {
      if(head->type==LALINFERENCE_gslMatrix_t)
      {
        matrix = *((gsl_matrix **)head->value);
        for(i=0; i<(int)matrix->size1; i++)
        {
          for(j=0; j<(int)matrix->size2; j++)
          {
            fprintf(out, "%s_%s%i%i\t", LALInferenceTranslateInternalToExternalParamName(head->name),suffix,i,j);
          }
        }
      }
      else fprintf(out, "%s_%s\t", LALInferenceTranslateInternalToExternalParamName(head->name), suffix);
    }
    head = head->next;
  }

  return 0;
}

int LALInferenceCompareVariables(LALInferenceVariables *var1, LALInferenceVariables *var2)
/*  Compare contents of "var1" and "var2".                       */
/*  Returns zero for equal entries, and one if difference found. */
/*  Make sure to only call this function when all entries are    */
/*  actually comparable. For example, "gslMatrix" type entries   */
/*  cannot (yet?) be checked for equality.                       */
{
  int result = 0;
  LALInferenceVariableItem *ptr1 = var1->head;
  LALInferenceVariableItem *ptr2 = NULL;
  if (var1->dimension != var2->dimension) result = 1;  // differing dimension
  while ((ptr1 != NULL) && (result == 0)) {
    ptr2 = LALInferenceGetItem(var2, ptr1->name);
    if (ptr2 != NULL) {  // corrsesponding entry exists; now compare type, then value:
      if (ptr2->type == ptr1->type) {  // entry type identical
        switch (ptr1->type) {  // do value comparison depending on type:
          case LALINFERENCE_INT4_t:
            result = ((*(INT4 *) ptr2->value) != (*(INT4 *) ptr1->value));
            break;
          case LALINFERENCE_INT8_t:
            result = ((*(INT8 *) ptr2->value) != (*(INT8 *) ptr1->value));
            break;
          case LALINFERENCE_UINT4_t:
            result = ((*(UINT4 *) ptr2->value) != (*(UINT4 *) ptr1->value));
            break;
          case LALINFERENCE_REAL4_t:
            result = ((*(REAL4 *) ptr2->value) != (*(REAL4 *) ptr1->value));
            break;
          case LALINFERENCE_REAL8_t:
            result = ((*(REAL8 *) ptr2->value) != (*(REAL8 *) ptr1->value));
            break;
          case LALINFERENCE_COMPLEX8_t:
            result = (((REAL4) crealf(*(COMPLEX8 *) ptr2->value) != (REAL4) crealf(*(COMPLEX8 *) ptr1->value))
                      || ((REAL4) cimagf(*(COMPLEX8 *) ptr2->value) != (REAL4) cimagf(*(COMPLEX8 *) ptr1->value)));
            break;
          case LALINFERENCE_COMPLEX16_t:
            result = (((REAL8) creal(*(COMPLEX16 *) ptr2->value) != (REAL8) creal(*(COMPLEX16 *) ptr1->value))
                      || ((REAL8) cimag(*(COMPLEX16 *) ptr2->value) != (REAL8) cimag(*(COMPLEX16 *) ptr1->value)));
            break;
          case LALINFERENCE_gslMatrix_t:
            if( matrix_equal(*(gsl_matrix **)ptr1->value,*(gsl_matrix **)ptr2->value) )
                result = 0;
            else
                result = 1;
            break;
          default:
            XLAL_ERROR(XLAL_EFAILED, "Encountered unknown LALInferenceVariables type (entry: \"%s\").", ptr1->name);
        }
      }
      else result = 1;  // same name but differing type
    }
    else result = 1;  // entry of given name doesn't exist in var2
    ptr1 = ptr1->next;
  }
  return(result);
}

INT4 LALInferenceBufferToArray(LALInferenceRunState *state, INT4 startCycle, INT4 endCycle, REAL8** DEarray) {
  LALInferenceVariableItem *ptr;
  INT4 i=0,p=0;

  INT4 Nskip = *(INT4*) LALInferenceGetVariable(state->algorithmParams, "Nskip");
  INT4 totalPoints = state->differentialPointsLength;
  INT4 start = (INT4)ceil((REAL8)startCycle/(REAL8)Nskip);
  INT4 end = (INT4)floor((REAL8)endCycle/(REAL8)Nskip);
  /* Include last point */
  if (end > totalPoints-1)
    end = totalPoints-1;

  for (i = start; i <= end; i++) {
    ptr=state->differentialPoints[i]->head;
    p=0;
    while(ptr!=NULL) {
      if (LALInferenceCheckVariableNonFixed(state->differentialPoints[i], ptr->name) && ptr->type == LALINFERENCE_REAL8_t) {
        DEarray[i-start][p]=*(REAL8 *)ptr->value;
        p++;
      }
      ptr=ptr->next;
    }
  }
  return end-start+1;
}

void LALInferenceArrayToBuffer(LALInferenceRunState *runState, REAL8** DEarray) {
  LALInferenceVariableItem *ptr;
  UINT4 i=0,p=0;
  UINT4 nPoints = sizeof(DEarray) / sizeof(REAL8*);

  /* Save last LALInferenceVariables item from buffer to keep fixed params consistent for chain */
  LALInferenceVariables templateParamSet;
  LALInferenceCopyVariables(runState->differentialPoints[runState->differentialPointsLength-1], &templateParamSet);

  /* Free old DE buffer */
  XLALFree(runState->differentialPoints);

  /* Expand DE buffer */
  size_t newSize = runState->differentialPointsSize;
  while (nPoints > newSize) {
    newSize = newSize*2;
  }

  runState->differentialPoints = XLALCalloc(newSize, sizeof(LALInferenceVariables *));
  runState->differentialPointsLength = nPoints;
  runState->differentialPointsSize = newSize;

  for (i=0; i<nPoints; i++) {
    runState->differentialPoints[i] = XLALCalloc(1, sizeof(LALInferenceVariables));
    LALInferenceCopyVariables(&templateParamSet, runState->differentialPoints[i]);
    ptr = runState->differentialPoints[i]->head;
    while(ptr!=NULL) {
      if (LALInferenceCheckVariableNonFixed(runState->differentialPoints[i], ptr->name) && ptr->type == LALINFERENCE_REAL8_t) {
        *((REAL8 *)ptr->value) = (REAL8)DEarray[i][p];
        p++;
      }
      ptr=ptr->next;
    }
  }
}


REAL8Vector *LALInferenceCopyVariablesToArray(LALInferenceVariables *origin) {
  INT4 nPar = LALInferenceGetVariableDimensionNonFixed(origin);
  REAL8Vector * parameters = NULL;
  gsl_matrix *m = NULL; //for dealing with noise parameters
  UINT4Vector *v = NULL; //for dealing with dimension parameters
  UINT4 j,k;

  parameters = XLALCreateREAL8Vector(nPar);

  LALInferenceVariableItem *ptr=origin->head;
  INT4 p=0;
  while(ptr!=NULL) {
    if (LALInferenceCheckVariableNonFixed(origin, ptr->name)) {
      //Generalized to allow for parameters stored in gsl_matrix or UINT4Vector
      if(ptr->type == LALINFERENCE_gslMatrix_t)
      {
        m = *((gsl_matrix **)ptr->value);
        for(j=0; j<m->size1; j++)
        {
          for(k=0; k<m->size2; k++)
          {
            parameters->data[p]=gsl_matrix_get(m,j,k);
            p++;
          }
        }
      }
      else if(ptr->type == LALINFERENCE_UINT4Vector_t)
      {
        v = *(UINT4Vector **)ptr->value;
        for(j=0; j<v->length; j++)
        {
          parameters->data[p]=v->data[j];
          p++;
        }
      }
      else
      {
        parameters->data[p]=*(REAL8 *)ptr->value;
        p++;
      }
    }
    ptr=ptr->next;
  }

  return parameters;
}

void LALInferenceCopyArrayToVariables(REAL8Vector *origin, LALInferenceVariables *target) {
  gsl_matrix *m = NULL; //for dealing with noise parameters
  UINT4Vector *v = NULL; //for dealing with dimension parameters
  UINT4 j,k;

  LALInferenceVariableItem *ptr = target->head;
  INT4 p=0;
  while(ptr!=NULL) {
    if (LALInferenceCheckVariableNonFixed(target, ptr->name))
    {
      //Generalized to allow for parameters stored in gsl_matrix
      if(ptr->type == LALINFERENCE_gslMatrix_t)
      {
        m = *((gsl_matrix **)ptr->value);
        for(j=0; j<m->size1; j++)
        {
          for(k=0; k<m->size2; k++)
          {
            gsl_matrix_set(m,j,k,origin->data[p]);
            p++;
          }
        }
      }
      else if(ptr->type == LALINFERENCE_UINT4Vector_t)
      {
        v = *(UINT4Vector **)ptr->value;
        for(j=0; j<v->length; j++)
        {
          v->data[j] = origin->data[p];
          p++;
        }
      }
      else
      {
        memcpy(ptr->value,&(origin->data[p]),LALInferenceTypeSize[ptr->type]);
        p++;
      }
    }
    ptr=ptr->next;
  }
  return;
}

/* ============ Command line parsing functions etc.: ========== */



ProcessParamsTable *LALInferenceGetProcParamVal(ProcessParamsTable *procparams,const char *name)
/* Returns pointer to element "name" of the ProcessParamsTable, */
/* if present, and NULL otherwise.                              */
{
  ProcessParamsTable *this=procparams;

  if (this==NULL) {
    fprintf(stderr, " Warning:  ProcessParamsTable is a NULL pointer\n");
    exit(1);
  }

  while (this!=NULL) {
    if (!strcmp(this->param, name)) break;
    else this=this->next;
  }

  return(this);
}



void LALInferenceParseCharacterOptionString(char *input, char **strings[], UINT4 *n)
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
  if (j!=2) XLAL_ERROR_VOID(XLAL_EINVAL, "Argument vector \"%s\" is not well-formed!", input);
  /* now allocate memory for results: */
  *strings  = (char**)  XLALMalloc(sizeof(char*) * (*n));
  for (i=0; i<(*n); ++i) (*strings)[i] = (char*) XLALMalloc(sizeof(char)*512);
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
        XLAL_PRINT_WARNING("Character \"%s\" argument too long!", (*strings)[k]);
      }
      else {
        (*strings)[k][l] = input[i];
        ++l;
      }
    }
    ++i;
  }
}

ProcessParamsTable *LALInferenceParseCommandLine(int argc, char *argv[])
/* parse command line and set up & fill in 'ProcessParamsTable' linked list.          */
/* If no command line arguments are supplied, the 'ProcessParamsTable' still contains */
/* one empty entry.                                                                   */
{
  int i, state=1;
  int dbldash;
  ProcessParamsTable *head, *ptr=NULL;
  /* always (even for argc==1, i.e. no arguments) put one element in list: */
  head = (ProcessParamsTable*) XLALCalloc(1, sizeof(ProcessParamsTable));
  XLALStringCopy(head->program, argv[0], sizeof(CHAR)*LIGOMETA_PROGRAM_MAX);
  ptr = head;
  i=1;
  while ((i<argc) & (state<=3)) {
    /* check for a double-dash at beginning of argument #i: */
    dbldash = ((argv[i][0]=='-') && (argv[i][1]=='-'));
    /* react depending on current state: */
    if (state==1){ /* ('state 1' means handling very 1st argument) */
      if (dbldash) {
        XLALStringCopy(head->param, argv[i], sizeof(CHAR)*LIGOMETA_PARAM_MAX);
        XLALStringCopy(ptr->type, "string", sizeof(CHAR)*LIGOMETA_TYPE_MAX);
        state = 2;
      }
      else { /* (very 1st argument needs to start with "--...") */
        // TODO: Perhaps this should be a warning?
        XLAL_ERROR_NULL(XLAL_EINVAL, "Orphaned first command line argument: \"%s\".", argv[i]);
        state = 4;
      }
    }
    else if (state==2) { /* ('state 2' means last entry was a parameter starting with "--") */
      if (dbldash) {
        ptr->next = (ProcessParamsTable*) XLALCalloc(1, sizeof(ProcessParamsTable));
        ptr = ptr->next;
        XLALStringCopy(ptr->program, argv[0],
sizeof(CHAR)*LIGOMETA_PROGRAM_MAX);
        XLALStringCopy(ptr->param, argv[i], sizeof(CHAR)*LIGOMETA_PARAM_MAX);
        XLALStringCopy(ptr->type, "string", sizeof(CHAR)*LIGOMETA_TYPE_MAX);
      }
      else {
        state = 3;
        XLALStringCopy(ptr->value, argv[i], sizeof(CHAR)*LIGOMETA_VALUE_MAX);
      }
    }
    else if (state==3) { /* ('state 3' means last entry was a value) */
      if (dbldash) {
        ptr->next = (ProcessParamsTable*) XLALCalloc(1, sizeof(ProcessParamsTable));
        ptr = ptr->next;
        XLALStringCopy(ptr->program, argv[0],
                       sizeof(CHAR)*LIGOMETA_PROGRAM_MAX);
        XLALStringCopy(ptr->param, argv[i], sizeof(CHAR)*LIGOMETA_PARAM_MAX);
        XLALStringCopy(ptr->type, "string", sizeof(CHAR)*LIGOMETA_TYPE_MAX);
        state = 2;
      }
      else {
        // TODO: Perhaps this should be a warning?
        XLAL_ERROR_NULL(XLAL_EINVAL, "Orphaned first command line argument: \"%s\".", argv[i]);
        state = 4;
      }
    }
    ++i;
  }

  return head;
}

ProcessParamsTable *LALInferenceParseCommandLineStringVector(LALStringVector* args){
  /* Version of the below function which can be easily exposed to Python via SWIG. */
  return LALInferenceParseCommandLine(args->length,args->data);
}

char* LALInferencePrintCommandLine(ProcessParamsTable *procparams)
{
  ProcessParamsTable *this=procparams;
  INT8 len=14; //number of characters of the "Command line: " string.
  while (this!=NULL) {
    len+=strlen(this->param);
    len+=strlen(this->value);
    len+=2;
    this=this->next;
  }// Now we know how long the buffer has to be.
  char * str = (char*) XLALCalloc(len+1,sizeof(char));
  if (str==NULL) {
    XLALPrintError("Calloc error, str is NULL (in %s, line %d)\n",__FILE__, __LINE__);
		XLAL_ERROR_NULL(XLAL_ENOMEM);
  }
  
  this=procparams;
  strcpy (str,"Command line: ");
  //strcat (str,this->program);
  while (this!=NULL) {
    strcat (str," ");
    strcat (str,this->param);
    strcat (str," ");
    strcat (str,this->value);
    this=this->next;
  }
  return str;
}

void LALInferenceExecuteFT(LALInferenceModel *model)
/* Execute (forward, time-to-freq) Fourier transform.  Contents of
model->timeh... are windowed and FT'ed, results go into
model->freqh...  

NOTE: the windowing is performed *in-place*, so do not call more than
once on a given timeModel!
*/
{
  UINT4 i;
  double norm;
  int errnum; 
  
  if (model==NULL) {
    fprintf(stderr," ERROR: model is a null pointer at LALInferenceExecuteFT, exiting!.\n");
    XLAL_ERROR_VOID(XLAL_EFAULT);
  }

  if(!model->freqhPlus && !model->timehPlus){
    XLALPrintError("freqhPlus and timeqhPlus are NULL at LALInferenceExecuteFT, exiting!");
    XLAL_ERROR_VOID(XLAL_EFAULT);
  }
  else if(!model->freqhCross && !model->timehCross){
    XLALPrintError("freqhCross and timeqhCross are NULL at LALInferenceExecuteFT, exiting!");
    XLAL_ERROR_VOID(XLAL_EFAULT);
  }
 
  /* h+ */
  if(!model->freqhPlus){     
      XLAL_TRY(model->freqhPlus=(COMPLEX16FrequencySeries *)XLALCreateCOMPLEX16FrequencySeries("freqData",&(model->timehPlus->epoch),0.0,model->deltaF,&lalDimensionlessUnit,model->freqLength),errnum);

    if (errnum){
      XLALPrintError("Could not create COMPLEX16FrequencySeries in LALInferenceExecuteFT");
      XLAL_ERROR_VOID(errnum);
    }
  }

  if (!model->window || !model->window->data){
    XLALPrintError("model->window is NULL at LALInferenceExecuteFT: Exiting!");
    XLAL_ERROR_VOID(XLAL_EFAULT);
  }

  XLAL_TRY(XLALDDVectorMultiply(model->timehPlus->data, model->timehPlus->data, model->window->data), errnum);

  if (errnum){
    XLALPrintError("Could not window time-series in LALInferenceExecuteFT");
    XLAL_ERROR_VOID(errnum);
  }
  	
  if (!model->timeToFreqFFTPlan){
    XLALPrintError("model->timeToFreqFFTPlan is NULL at LALInferenceExecuteFT: Exiting!");
    XLAL_ERROR_VOID(XLAL_EFAULT);
  }

  XLAL_TRY(XLALREAL8TimeFreqFFT(model->freqhPlus, model->timehPlus, model->timeToFreqFFTPlan), errnum);
  
  if (errnum){
    XLALPrintError("Could not h_plus FFT time-series");
    XLAL_ERROR_VOID(errnum);
  }
  			    
  /* hx */
  if(!model->freqhCross){ 
    XLAL_TRY(model->freqhCross=(COMPLEX16FrequencySeries *)XLALCreateCOMPLEX16FrequencySeries("freqData",&(model->timehCross->epoch),0.0,model->deltaF,&lalDimensionlessUnit,model->freqLength),errnum);
  
    if (errnum){	
      XLALPrintError("Could not create COMPLEX16FrequencySeries in LALInferenceExecuteFT");
      XLAL_ERROR_VOID(errnum);		
    }
  }

  XLAL_TRY(XLALDDVectorMultiply(model->timehCross->data, model->timehCross->data, model->window->data), errnum);

  if (errnum){
    XLALPrintError("Could not window time-series in LALInferenceExecuteFT");
    XLAL_ERROR_VOID(errnum);
  }
  	 
  XLAL_TRY(XLALREAL8TimeFreqFFT(model->freqhCross, model->timehCross, model->timeToFreqFFTPlan), errnum);
  
  if (errnum){
    XLALPrintError("Could not FFT h_cross time-series");
    XLAL_ERROR_VOID(errnum);
  }   

  norm=sqrt(model->window->data->length/model->window->sumofsquares);
  
  for(i=0;i<model->freqhPlus->data->length;i++){
    model->freqhPlus->data->data[i] *= ((REAL8) norm);
    model->freqhCross->data->data[i] *= ((REAL8) norm);
  }
}

void LALInferenceExecuteInvFT(LALInferenceModel *model)
/* Execute inverse (freq-to-time) Fourier transform. */
/* Results go into 'model->timeh...'                 */
{
  if (model->freqToTimeFFTPlan==NULL) {
    XLAL_ERROR_VOID(XLAL_EFAULT, "Encountered unallocated \"freqToTimeFFTPlan\".");
  }

  /*  h+ :  */
  if (model->timehPlus==NULL) {
    XLAL_ERROR_VOID(XLAL_EFAULT, "Encountered unallocated \"timehPlus\".");
  }
  if (model->freqhPlus==NULL) {
    XLAL_ERROR_VOID(XLAL_EFAULT, "Encountered unallocated \"freqhPlus\".");
  }

  XLALREAL8FreqTimeFFT(model->timehPlus, model->freqhPlus, model->freqToTimeFFTPlan);

  /*  hx :  */
  if (model->timehCross==NULL) {
    XLAL_ERROR_VOID(XLAL_EFAULT, "Encountered unallocated \"timehCross\".");
  }
  if (model->freqhCross==NULL) {
    XLAL_ERROR_VOID(XLAL_EFAULT, "Encountered unallocated \"freqhCross\".");
  }

  XLALREAL8FreqTimeFFT(model->timehCross, model->freqhCross, model->freqToTimeFFTPlan);
}

void LALInferenceProcessParamLine(FILE *inp, char **headers, LALInferenceVariables *vars) {
  size_t i;

  for (i = 0; headers[i] != NULL; i++) {
    double param;
    int nread;

    nread = fscanf(inp, " %lg ", &param);
    if (nread != 1) {
      XLAL_ERROR_VOID(XLAL_EFAILED, "Could not read the value of the %zu parameter in the row.", i);
    }

    LALInferenceAddVariable(vars, headers[i], &param, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  }
}

/* This function has a Memory Leak!  You cannot free the allocated
   header buffer (of length MAXSIZE).  Don't call it too many times!
   (It's only expected to be called once to initialize the
   differential evolution array, so this should be OK. */
char **LALInferenceGetHeaderLine(FILE *inp) {
  const size_t MAXSIZE=1024;
  const char *delimiters = " \n\t";
  char *header = XLALMalloc(MAXSIZE*sizeof(char));
  char **colNames = NULL;  /* Will be filled in with the column names,
                              terminated by NULL. */
  size_t colNamesLen=0, colNamesMaxLen=0;
  char *colName = NULL;

  if (!fgets(header, MAXSIZE, inp)) {
    /* Some error.... */
    XLAL_ERROR_NULL(XLAL_EFAILED, "Error reading header line from file.");
  } else if (strlen(header) >= MAXSIZE-1) {
    /* Probably ran out of space before reading the entire line. */
    XLAL_ERROR_NULL(XLAL_EFAILED, "Header line too long (more than %zu chars).", MAXSIZE - 1);
  }

  /* Sure hope we read the whole line. */
  colNamesMaxLen=2;
  colNames=(char **)XLALMalloc(2*sizeof(char *));

  if (!colNames) {
    XLAL_ERROR_NULL(XLAL_ENOMEM, "Failed to allocate memory for colNames.");
  }

  colName=strtok(header, delimiters);
  strcpy(colNames[0],colNameToParamName(colName));
  //colNames[0] = colNameToParamName(colName); /* switched to strcpy() to avoid warning: assignment discards qualifiers from pointer target type */
  colNamesLen=1;
  do {
    colName=strtok(NULL, delimiters);

    strcpy(colNames[colNamesLen],colNameToParamName(colName));
    colNamesLen++;

    /* Expand if necessary. */
    if (colNamesLen >= colNamesMaxLen) {
      colNamesMaxLen *= 2;
      colNames=XLALRealloc(colNames, colNamesMaxLen*sizeof(char *));
      if (!colNames) {
        XLAL_ERROR_NULL(XLAL_ENOMEM, "Failed to XLALReallocate memory for colNames.");
      }
    }

  } while (colName != NULL);

  /* Trim down to size. */
  colNames=XLALRealloc(colNames, colNamesLen*sizeof(char *));

  return colNames;
}

char *colNameToParamName(const char *colName) {
  char *retstr=NULL;
  if (colName == NULL) {
    return NULL;
  }
  else if (!strcmp(colName, "dist")) {
    retstr=XLALStringDuplicate("distance");
  }

  else if (!strcmp(colName, "ra")) {
    retstr=XLALStringDuplicate("rightascension");
  }

  else if (!strcmp(colName, "iota")) {
    retstr=XLALStringDuplicate("inclination");
  }

  else if (!strcmp(colName, "psi")) {
    retstr=XLALStringDuplicate("polarisation");
  }

  else if (!strcmp(colName, "mc")) {
    retstr=XLALStringDuplicate("chirpmass");
  }

  else if (!strcmp(colName, "phi_orb")) {
    retstr=XLALStringDuplicate("phase");
  }

  else if (!strcmp(colName, "eta")) {
    retstr=XLALStringDuplicate("massratio");
  }

  else if (!strcmp(colName, "q")) {
    retstr=XLALStringDuplicate("asym_massratio");
  }

  else if (!strcmp(colName, "dec")) {
    retstr=XLALStringDuplicate("declination");
  }

  /* Note the 1 <--> 2 swap between the post-proc world and the LI world. */
  else if (!strcmp(colName, "phi1")) {
    retstr=XLALStringDuplicate("phi_spin2");
  }

  else if (!strcmp(colName, "phi2")) {
    retstr=XLALStringDuplicate("phi_spin1");
  }

  else if (!strcmp(colName, "theta1")) {
    retstr=XLALStringDuplicate("theta_spin2");
  }

  else if (!strcmp(colName, "theta2")) {
    retstr=XLALStringDuplicate("theta_spin1");
  }

  else if (!strcmp(colName, "a1")) {
    retstr=XLALStringDuplicate("a_spin2");
  }
    
  else if (!strcmp(colName, "a2")) {
    retstr=XLALStringDuplicate("a_spin1");
  }
  else retstr=XLALStringDuplicate(colName);
  return retstr;
}

void LALInferenceSortVariablesByName(LALInferenceVariables *vars)
{
  LALInferenceVariables tmp;
  tmp.head=NULL;
  tmp.dimension=0;
  LALInferenceVariableItem *thisitem,*ptr;
  LALInferenceVariables *new=XLALCalloc(1,sizeof(*new));
  if(!vars){
    XLAL_ERROR_VOID(XLAL_EFAULT, "Received null input pointer.");
  }
  while(vars->head)
  {
    thisitem=vars->head;
    for (ptr=thisitem->next;ptr;ptr=ptr->next){
      if(strcmp(ptr->name,thisitem->name)<0)
        thisitem=ptr;
    }
    LALInferenceAddVariable(&tmp, thisitem->name, thisitem->value, thisitem->type, thisitem->vary);
    LALInferenceRemoveVariable(vars,thisitem->name);
  }
  vars->head=tmp.head;
  vars->dimension=tmp.dimension;
  return;
}

/**
 * Append the sample to a file. file pointer is stored in state->algorithmParams as a
 * LALInferenceVariable called "outfile", as a void ptr.
 * Caller is responsible for opening and closing file.
 * Variables are alphabetically sorted before being written
 */
void LALInferenceLogSampleToFile(LALInferenceRunState *state, LALInferenceVariables *vars)
{
  FILE *outfile=NULL;
  if(LALInferenceCheckVariable(state->algorithmParams,"outfile"))
    outfile=*(FILE **)LALInferenceGetVariable(state->algorithmParams,"outfile");
  /* Write out old sample */
  if(outfile==NULL) return;
  LALInferenceSortVariablesByName(vars);
  LALInferencePrintSample(outfile,vars);
  fprintf(outfile,"\n");
}

/**
 * Append the sample to an array which can be later processed by the user.
 * Array is stored as a C array in a LALInferenceVariable in state->algorithmParams
 * called "outputarray". Number of items in the array is stored as "N_outputarray".
 * Will create the array and store it in this way if it does not exist.
 * DOES NOT FREE ARRAY, user must clean up after use.
 * Also outputs sample to disk if possible
 */
void LALInferenceLogSampleToArray(LALInferenceRunState *state, LALInferenceVariables *vars)
{
  LALInferenceVariables *output_array=NULL;
  UINT4 N_output_array=0;
  LALInferenceSortVariablesByName(vars);
  LALInferenceLogSampleToFile(state,vars);

  /* Set up the array if it is not already allocated */
  if(LALInferenceCheckVariable(state->algorithmParams,"outputarray"))
    output_array=*(LALInferenceVariables **)LALInferenceGetVariable(state->algorithmParams,"outputarray");
  else
    LALInferenceAddVariable(state->algorithmParams,"outputarray",&output_array,LALINFERENCE_void_ptr_t,LALINFERENCE_PARAM_OUTPUT);

  if(LALInferenceCheckVariable(state->algorithmParams,"N_outputarray"))
    N_output_array=*(INT4 *)LALInferenceGetVariable(state->algorithmParams,"N_outputarray");
  else
    LALInferenceAddVariable(state->algorithmParams,"N_outputarray",&N_output_array,LALINFERENCE_INT4_t,LALINFERENCE_PARAM_OUTPUT);

  /* Expand the array for new sample */
  output_array=XLALRealloc(output_array, (N_output_array+1) *sizeof(LALInferenceVariables));
  if(!output_array){
    XLAL_ERROR_VOID(XLAL_EFAULT, "Unable to allocate array for samples.");
  }
  else
  {
    /* Save sample and update */
    memset(&(output_array[N_output_array]),0,sizeof(LALInferenceVariables));
    LALInferenceCopyVariables(vars,&output_array[N_output_array]);
    N_output_array++;

    LALInferenceSetVariable(state->algorithmParams,"outputarray",&output_array);
    LALInferenceSetVariable(state->algorithmParams,"N_outputarray",&N_output_array);
  }
  return;
}

void LALInferenceMcEta2Masses(double mc, double eta, double *m1, double *m2)
/*  Compute individual companion masses (m1, m2)   */
/*  for given chirp mass (m_c) & mass ratio (eta)  */
/*  (note: m1 >= m2).                              */
{
  double root = sqrt(0.25-eta);
  double fraction = (0.5+root) / (0.5-root);
  *m2 = mc * (pow(1+fraction,0.2) / pow(fraction,0.6));
  *m1 = mc * (pow(1+1.0/fraction,0.2) / pow(1.0/fraction,0.6));
  return;
}

void LALInferenceMcQ2Masses(double mc, double q, double *m1, double *m2)
/*  Compute individual companion masses (m1, m2)   */
/*  for given chirp mass (m_c) & asymmetric mass   */
/*  ratio (q).  note: q = m2/m1, where m1 >= m2    */
{
  double factor = mc * pow(1 + q, 1.0/5.0);
  *m1 = factor * pow(q, -3.0/5.0);
  *m2 = factor * pow(q, +2.0/5.0);
  return;
}
void LALInferenceQ2Eta(double q, double *eta)
/*  Compute symmetric mass ratio eta from the     */
/*  asymmetric mass ratio q.                      */
{
  *eta = q/((1+q)*(1+q));
  return;
}

void LALInferenceLambdaTsEta2Lambdas(REAL8 lambdaT, REAL8 dLambdaT, REAL8 eta, REAL8 *lambda1, REAL8 *lambda2){
  REAL8 a=(8./13.)*(1.+7.*eta-31.*eta*eta);
  REAL8 b=(8./13.)*sqrt(1.-4.*eta)*(1.+9.*eta-11.*eta*eta);
  REAL8 c=(1./2.)*sqrt(1.-4.*eta)*(1.-13272.*eta/1319.+8944.*eta*eta/1319.);
  REAL8 d=(1./2.)*(1.-15910.*eta/1319.+32850.*eta*eta/1319.+3380.*eta*eta*eta/1319.);
  *lambda1=((c-d)*lambdaT-(a-b)*dLambdaT)/(2.*(b*c-a*d));
  *lambda2=((c+d)*lambdaT-(a+b)*dLambdaT)/(2.*(a*d-b*c));
  return;
}

static void deleteCell(LALInferenceKDTree *cell) {
  if (cell == NULL) {
    return; /* Our work here is done. */
  } else {
    size_t i;

    deleteCell(cell->left);
    deleteCell(cell->right);

    if (cell->lowerLeft != NULL) XLALFree(cell->lowerLeft);
    if (cell->upperRight != NULL) XLALFree(cell->upperRight);
    if (cell->ptsMean != NULL) XLALFree(cell->ptsMean);
    if (cell->eigenMin != NULL) XLALFree(cell->eigenMin);
    if (cell->eigenMax != NULL) XLALFree(cell->eigenMax);

    for (i = 0; i < cell->ptsSize; i++) {
      if (cell->pts[i] != NULL) XLALFree(cell->pts[i]);
    }
    if (cell->pts != NULL) XLALFree(cell->pts);

    if (cell->ptsCov != NULL) {
      for (i = 0; i < cell->dim; i++) {
        XLALFree(cell->ptsCov[i]);
      }
      XLALFree(cell->ptsCov);
    }

    if (cell->ptsCovEigenVects != NULL) {
      for (i = 0; i < cell->dim; i++) {
        XLALFree(cell->ptsCovEigenVects[i]);
      }
      XLALFree(cell->ptsCovEigenVects);
    }

    XLALFree(cell);

    return;
  }
}

void LALInferenceKDTreeDelete(LALInferenceKDTree *tree) {
  deleteCell(tree);
}

typedef enum {
  LEFT,
  RIGHT,
  TOP
} cellType;

static LALInferenceKDTree *newCell(size_t ndim, REAL8 *lowerLeft, REAL8 *upperRight, size_t level, cellType type) {
  LALInferenceKDTree *cell = XLALCalloc(1, sizeof(LALInferenceKDTree));
  size_t i;

  cell->lowerLeft = XLALCalloc(ndim, sizeof(REAL8));
  cell->upperRight = XLALCalloc(ndim, sizeof(REAL8));
  cell->ptsMean = XLALCalloc(ndim, sizeof(REAL8));
  cell->eigenMin = XLALCalloc(ndim, sizeof(REAL8));
  cell->eigenMax = XLALCalloc(ndim, sizeof(REAL8));

  cell->pts = XLALCalloc(1, sizeof(REAL8 *));
  cell->ptsSize = 1;
  cell->npts = 0;
  cell->dim = ndim;

  cell->ptsCov = XLALCalloc(ndim, sizeof(REAL8 *));
  cell->ptsCovEigenVects = XLALCalloc(ndim, sizeof(REAL8 *));
  for (i = 0; i < ndim; i++) {
    cell->ptsCov[i] = XLALCalloc(ndim, sizeof(REAL8));
    cell->ptsCovEigenVects[i] = XLALCalloc(ndim, sizeof(REAL8));
  }
  
  memcpy(cell->upperRight, upperRight, ndim*sizeof(REAL8));
  memcpy(cell->lowerLeft, lowerLeft, ndim*sizeof(REAL8));
  if (type == LEFT) {
    cell->upperRight[level] = 0.5*(lowerLeft[level] + upperRight[level]);
  } else if (type == RIGHT) {
    cell->lowerLeft[level] = 0.5*(lowerLeft[level] + upperRight[level]);
  } else {
    /* Do not change lowerLeft or upperRight, since this is the top-level cell. */
  }

  return cell;
}

LALInferenceKDTree *LALInferenceKDEmpty(REAL8 *lowerLeft, REAL8 *upperRight, size_t ndim) {
  LALInferenceKDTree *cell = newCell(ndim, lowerLeft, upperRight, 0, TOP);
  return cell;
}

static int equalPoints(REAL8 *a, REAL8 *b, size_t n) {
  size_t i;
  for (i = 0; i < n; i++) {
    if (a[i] != b[i]) return 0;
  }
  return 1;
}

static int cellAllEqualPoints(LALInferenceKDTree *cell) {
  if (cell->npts <= 1) {
    return 1;
  } else {
    size_t i;
    REAL8 *pt0 = cell->pts[0];
    
    for (i = 1; i < cell->npts; i++) {
      if (!equalPoints(pt0, cell->pts[i], cell->dim)) return 0;
    }

    return 1;
  }
}

static void addPtToCellPts(LALInferenceKDTree *cell, REAL8 *pt) {
  size_t ptsSize = cell->ptsSize;
  size_t npts = cell->npts;
  size_t dim = cell->dim;
  
  /* copy previous points */
  if ( npts == ptsSize ){
    REAL8 tmpArr[npts][dim];
    
    /* copy points from cell */
    for( UINT4 i=0; i < npts; i++ ){
      for( UINT4 j=0; j < dim; j++ ){
        tmpArr[i][j] = cell->pts[i][j];
      }
      XLALFree(cell->pts[i]); /* free column */
    }
    
    /* free array */
    XLALFree(cell->pts);
  
    /* expand array */
    cell->pts = XLALCalloc(2*ptsSize, sizeof(REAL8 *));
  
    /* copy vector into array */
    for( UINT4 i=0; i < 2*ptsSize; i++ ){
      cell->pts[i] = XLALCalloc(dim, sizeof(REAL8));
      
      if (i < npts){
        for( UINT4 j=0; j < dim; j++ )
          cell->pts[i][j] = tmpArr[i][j];
      } 
    }
    
    cell->ptsSize *= 2;
  }
  
  if ( npts == 0 ) cell->pts[npts] = XLALCalloc(dim, sizeof(REAL8));
  
  /* add new point */
  for( UINT4 i = 0; i < dim; i++ )
    cell->pts[npts][i] = pt[i];
  
  cell->npts += 1;
  cell->eigenFrameStale = 1;
}

/* The following routine handles inserting new points into the tree.
   It is organized recursively, inserting into cells based on the
   following three cases:

   1. There are no points in the cell.  Then insert the given point,
   and stop; we're at a leaf.

   2. Both of the cell's sub-cells are NULL.  Then the cell is a leaf
   (with some number of points), and the first order of business is to
   push the existing points down a level, then insert into this cell,
   and the appropriate sub-cell.

   3. This cell has non-null sub-cells.  Then it is not a leaf cell,
   and we should just insert into this cell, and add the point to the
   appropriate sub-cell.

*/
static int insertIntoCell(LALInferenceKDTree *cell, REAL8 *pt, size_t level) {
  size_t dim = cell->dim;
  size_t nextLevel = (level+1)%dim;
  level = level%dim;

  if (cell->npts == 0) {
    /* Insert this point into the cell, and quit. */
    addPtToCellPts(cell, pt);
    return XLAL_SUCCESS;
  } else if (cell->left == NULL && cell->right == NULL) {    
    /* This cell is a leaf node.  Insert the point, then (unless the
       cell stores many copies of the same point), push everything
       down a level. */
    addPtToCellPts(cell, pt);

    if (cellAllEqualPoints(cell)) {
      /* Done, since there are many copies of the same point---cell
         remains a leaf. */
      return XLAL_SUCCESS;
    } else {
      size_t i;
      REAL8 mid = 0.5*(cell->lowerLeft[level] + cell->upperRight[level]);
      
      cell->left = newCell(dim, cell->lowerLeft, cell->upperRight, level, LEFT);
      cell->right = newCell(dim, cell->lowerLeft, cell->upperRight, level, RIGHT);

      for (i = 0; i < cell->npts; i++) {
        REAL8 *lowerPt = cell->pts[i];
        
        if (lowerPt[level] <= mid) {
          insertIntoCell(cell->left, lowerPt, nextLevel);
        } else {
          insertIntoCell(cell->right, lowerPt, nextLevel);
        }
      }

      return XLAL_SUCCESS;
    }
  } else {    
    /* This is not a leaf cell, so insert, and then move down the tree. */
    REAL8 mid = 0.5*(cell->lowerLeft[level] + cell->upperRight[level]);

    addPtToCellPts(cell, pt);
    
    if (pt[level] <= mid) {
      return insertIntoCell(cell->left, pt, nextLevel);
    } else {
      return insertIntoCell(cell->right, pt, nextLevel);
    }
  }
}

static int inBounds(REAL8 *pt, REAL8 *low, REAL8 *high, size_t n) {
  size_t i;

  for (i = 0; i < n; i++) {
    if (pt[i] < low[i] || pt[i] > high[i]) return 0;
  }

  return 1;
}

static void computeMean(LALInferenceKDTree *cell) {
  REAL8 **pts = cell->pts;
  REAL8 *mean = cell->ptsMean;
  size_t dim = cell->dim;
  size_t npts = cell->npts;
  size_t i;

  for (i = 0; i < dim; i++) {
    mean[i] = 0.0;
  }

  for (i = 0; i < npts; i++) {
    size_t j;
    for (j = 0; j < dim; j++) {
      mean[j] += pts[i][j];
    }
  }

  for (i = 0; i < dim; i++) {
    mean[i] /= npts;
  }
}

static void computeCovariance(LALInferenceKDTree *cell) {
  REAL8 **cov = cell->ptsCov;
  REAL8 **pts = cell->pts;
  REAL8 *mu = cell->ptsMean;
  size_t npts = cell->npts;
  size_t dim = cell->dim;
  size_t i;

  for (i = 0; i < dim; i++) {
    size_t j;
    for (j = 0; j < dim; j++) {
      cov[i][j] = 0.0;
    }
  }

  for (i = 0; i < npts; i++) {
    size_t j;

    REAL8 *pt = pts[i];

    for (j = 0; j < dim; j++) {
      size_t k;
      REAL8 ptj = pt[j];
      REAL8 muj = mu[j];

      for (k = 0; k < dim; k++) {
        REAL8 ptk = pt[k];
        REAL8 muk = mu[k];

        cov[j][k] += (ptj - muj)*(ptk-muk);
      }
    }
  }

  for (i = 0; i < cell->dim; i++) {
    size_t j;
    for (j = 0; j < cell->dim; j++) {
      cov[i][j] /= (npts - 1);
    }
  }
}

static void computeEigenVectorsCleanup(gsl_matrix *A, gsl_matrix *evects, gsl_vector *evals,
                                       gsl_eigen_symmv_workspace *ws) {
  if (A != NULL) gsl_matrix_free(A);
  if (evects != NULL) gsl_matrix_free(evects);
  if (evals != NULL) gsl_vector_free(evals);
  if (ws != NULL) gsl_eigen_symmv_free(ws);
}

static void computeEigenVectors(LALInferenceKDTree *cell) {
  size_t dim = cell->dim;
  gsl_matrix *A = gsl_matrix_alloc(dim, dim);
  gsl_matrix *evects = gsl_matrix_alloc(dim, dim);
  gsl_vector *evals = gsl_vector_alloc(dim);
  gsl_eigen_symmv_workspace *ws = gsl_eigen_symmv_alloc(dim);
  
  REAL8 **covEVs = cell->ptsCovEigenVects;
  REAL8 **cov = cell->ptsCov;

  size_t i;
  int status;

  if (A == NULL || evects == NULL || evals == NULL || ws == NULL) {
    computeEigenVectorsCleanup(A, evects, evals, ws);
    XLAL_ERROR_VOID(XLAL_ENOMEM);
  }

  /* Copy covariance matrix into A. */
  for (i = 0; i < dim; i++) {
    size_t j;

    for (j = 0; j < dim; j++) {
      gsl_matrix_set(A, i, j, cov[i][j]);
    }
  }

  /* Compute evecs. */
  if ((status = gsl_eigen_symmv(A, evals, evects, ws)) != GSL_SUCCESS) {
    computeEigenVectorsCleanup(A, evects, evals, ws);
    XLAL_ERROR_VOID(status, "gsl error");
  }

  /* Copy eigenvector matrix into covEVs; [i][j] is the jth component
     of the ith eigenvector. */
  for (i = 0; i < dim; i++) {
    size_t j;
    
    for (j = 0; j < dim; j++) {
      covEVs[i][j] = gsl_matrix_get(evects, j, i);
    }
  }

  computeEigenVectorsCleanup(A, evects, evals, ws);
}

/* xe = eigenvs^T x */
static void toEigenFrame( size_t dim,  REAL8 **eigenvs,  REAL8 *x, REAL8 *xe) {
  size_t j;

  memset(xe, 0, dim*sizeof(REAL8));

  for (j = 0; j < dim; j++) {
    size_t i;
    REAL8 xj = x[j];
    REAL8 *evj = eigenvs[j];
    for (i = 0; i < dim; i++) {
      xe[i] += evj[i]*xj;
    }
  }
}

/* x = eigenvs xe */
static void fromEigenFrame( size_t dim,  REAL8 **eigenvs,  REAL8 *xe, REAL8 *x) {
  size_t i;
  
  memset(x, 0, dim*sizeof(REAL8));

  for (i = 0; i < dim; i++) {
    size_t j;
     REAL8 *evi = eigenvs[i];
    for (j = 0; j < dim; j++) {
      x[i] += evi[j]*xe[j];
    }
  }
}

static void computeEigenMinMax(LALInferenceKDTree *cell) {
  REAL8 **pts = cell->pts;
  size_t dim = cell->dim;
  size_t npts = cell->npts;
  REAL8 **eigenvs = cell->ptsCovEigenVects;
  REAL8 *mu = cell->ptsMean;

  REAL8 *min = cell->eigenMin;
  REAL8 *max = cell->eigenMax;

  REAL8 *xe = NULL, *x = NULL;

  size_t i;

  xe = XLALCalloc(dim, sizeof(REAL8));
  if (xe == NULL) XLAL_ERROR_VOID(XLAL_ENOMEM);

  x = XLALCalloc(dim, sizeof(REAL8));
  if (x == NULL) {
    XLALFree(xe);
    XLAL_ERROR_VOID(XLAL_ENOMEM);
  }

  for (i = 0; i < dim; i++) {
    min[i] = 1.0/0.0;
    max[i] = -1.0/0.0;
  }

  for (i = 0; i < npts; i++) {
    size_t j;

    for (j = 0; j < dim; j++) {
      x[j] = pts[i][j] - mu[j];
    }

    toEigenFrame(dim, eigenvs, x, xe);

    for (j = 0; j < dim; j++) {
      if (xe[j] < min[j]) min[j] = xe[j];
      if (xe[j] > max[j]) max[j] = xe[j];
    }
  }

  XLALFree(x);
  XLALFree(xe);
}

static void updateEigenSystem(LALInferenceKDTree *cell) {
  computeMean(cell);
  computeCovariance(cell);
  computeEigenVectors(cell);
  computeEigenMinMax(cell);

  cell->eigenFrameStale = 0;
}

int LALInferenceKDAddPoint(LALInferenceKDTree *tree, REAL8 *pt) {
  if (tree == NULL) XLAL_ERROR(XLAL_EINVAL, "given NULL tree");

  if (!inBounds(pt, tree->lowerLeft, tree->upperRight, tree->dim))
    XLAL_ERROR(XLAL_EINVAL, "given point that is not in global tree bounds");
  
  return insertIntoCell(tree, pt, 0);
}

static LALInferenceKDTree *doFindCell(LALInferenceKDTree *cell, REAL8 *pt, size_t dim, size_t Npts, size_t level) {
  if (cell == NULL) {
    /* If we encounter a NULL cell, then pass it up the chain. */
    return cell;
  } else if (cell->npts == 0) {
    return NULL;
  } else if (cell->npts == 1 || cell->npts < Npts) {
    return cell;
  } else {
    REAL8 mid = 0.5*(cell->lowerLeft[level] + cell->upperRight[level]);

    if (pt[level] <= mid) {
      LALInferenceKDTree *maybeCell = doFindCell(cell->left, pt, dim, Npts, (level+1)%dim);
      if (maybeCell == NULL) {
        return cell; /* If a NULL comes up from below, then this cell
                        is the one with the fewest points containing
                        pt. */
      } else {
        return maybeCell;
      }
    } else {
      LALInferenceKDTree *maybeCell = doFindCell(cell->right, pt, dim, Npts, (level+1)%dim);
      if (maybeCell == NULL) {
        return cell;
      } else {
        return maybeCell;
      }
    }
  }
}

LALInferenceKDTree *LALInferenceKDFindCell(LALInferenceKDTree *tree, REAL8 *pt, size_t Npts) {
  return doFindCell(tree, pt, tree->dim, Npts, 0);
}

double LALInferenceKDLogCellVolume(LALInferenceKDTree *cell) {
  size_t ndim = cell->dim;
  size_t i;
  REAL8 logVol = 0.0;
  for (i = 0; i < ndim; i++) {
    logVol += log(cell->upperRight[i] - cell->lowerLeft[i]);
  }

  return logVol;
}

double LALInferenceKDLogCellEigenVolume(LALInferenceKDTree *cell) {
  double logVol = 0.0;
  size_t ndim = cell->dim;
  size_t i;

  if (cell->eigenFrameStale) {
    updateEigenSystem(cell);
  }
  
  for (i = 0; i < ndim; i++) {
    logVol += log(cell->eigenMax[i] - cell->eigenMin[i]);
  }

  return logVol;
}

void LALInferenceKDVariablesToREAL8(LALInferenceVariables *params, REAL8 *pt, LALInferenceVariables *templt) {
  LALInferenceVariableItem *templateItem = templt->head;
  size_t i = 0;
  while (templateItem != NULL) {
    if (LALInferenceCheckVariableNonFixed(templt, templateItem->name)) {
      pt[i] = *(REAL8 *)LALInferenceGetVariable(params, templateItem->name);
      i++;
    }
    templateItem = templateItem->next;
  }
}

void LALInferenceKDREAL8ToVariables(LALInferenceVariables *params, REAL8 *pt, LALInferenceVariables *templt) {
  LALInferenceVariableItem *templateItem = templt->head;
  size_t i = 0;
  while (templateItem != NULL) {
    if (LALInferenceCheckVariableNonFixed(templt, templateItem->name)) {
      LALInferenceSetVariable(params, templateItem->name, &(pt[i]));
      i++;
    }
    templateItem = templateItem->next;
  }
}

void LALInferenceKDDrawEigenFrame(gsl_rng *rng, LALInferenceKDTree *tree, REAL8 *pt, size_t Npts) {
  LALInferenceKDTree *topCell = tree;

  if (topCell == NULL || topCell->npts == 0) XLAL_ERROR_VOID(XLAL_EINVAL, "cannot draw from empty cell");

  LALInferenceKDTree *cell = LALInferenceKDFindCell(tree, topCell->pts[gsl_rng_uniform_int(rng, topCell->npts)], Npts);

  if (cell->npts == 1) {
    /* If there is one point, then the covariance matrix is undefined,
       and we draw from the rectangular area. */
    size_t i;
    
    for (i = 0; i < cell->dim; i++) {
      pt[i] = cell->lowerLeft[i] + gsl_rng_uniform(rng)*(cell->upperRight[i] - cell->lowerLeft[i]);
    }
  } else {
    REAL8 *ept = XLALCalloc(cell->dim, sizeof(REAL8));
    size_t i;

    if (cell->eigenFrameStale) {
      updateEigenSystem(cell);
    }

    do {
      for (i = 0; i < cell->dim; i++) {
        ept[i] = cell->eigenMin[i] + gsl_rng_uniform(rng)*(cell->eigenMax[i] - cell->eigenMin[i]);
      }
      
      fromEigenFrame(cell->dim, cell->ptsCovEigenVects, ept, pt);
      
      for (i = 0; i < cell->dim; i++) {
        pt[i] += cell->ptsMean[i];
      }
    } while (!inBounds(pt, cell->lowerLeft, cell->upperRight, cell->dim));

    XLALFree(ept);
  }
}

static int inEigenBox(LALInferenceKDTree *cell,  REAL8 *pt) {
  REAL8 *shiftedPt, *ept;
  size_t i;

  if (cell->eigenFrameStale) {
    updateEigenSystem(cell);
  }

  shiftedPt = XLALCalloc(cell->dim, sizeof(REAL8));
  ept = XLALCalloc(cell->dim, sizeof(REAL8));

  for (i = 0; i < cell->dim; i++) {
    shiftedPt[i] = pt[i] - cell->ptsMean[i];
  }

  toEigenFrame(cell->dim, cell->ptsCovEigenVects, shiftedPt, ept);

  for (i = 0; i < cell->dim; i++) {
    if (ept[i] < cell->eigenMin[i] || ept[i] > cell->eigenMax[i]) {
      XLALFree(shiftedPt);
      XLALFree(ept);
      return 0;
    }
  }

  XLALFree(shiftedPt);
  XLALFree(ept);
  return 1;
}

REAL8 LALInferenceKDLogProposalRatio(LALInferenceKDTree *tree, REAL8 *current, 
                                     REAL8 *proposed, size_t Npts) {
  LALInferenceKDTree *currentCell = LALInferenceKDFindCell(tree, current, Npts);
  LALInferenceKDTree *proposedCell = LALInferenceKDFindCell(tree, proposed, Npts);

  LALInferenceKDTree *topCell = tree;

  size_t npts = topCell->npts;

  REAL8 logCurrentVolume, logProposedVolume;

  if (currentCell->npts > 1) {
    if (currentCell->eigenFrameStale) {
      updateEigenSystem(currentCell);
    }
    
    if (!inEigenBox(currentCell, current)) {
      return log(0.0); /* If the current point is not in the eigen box
                          of the current cell, then we cannot jump
                          back there.  */
    } else {
      logCurrentVolume = LALInferenceKDLogCellEigenVolume(currentCell);
    }
  } else {
    logCurrentVolume = LALInferenceKDLogCellVolume(currentCell);
  }

  if (proposedCell->npts > 1) {
    /* Since we just proposed out of this cell, its eigenframe should
       *not* be stale. */
    if (proposedCell->eigenFrameStale) {
      XLAL_ERROR_REAL8(XLAL_EINVAL, "proposed cell eigen-frame is stale");
    }

    logProposedVolume = LALInferenceKDLogCellEigenVolume(proposedCell);
  } else {
    logProposedVolume = LALInferenceKDLogCellVolume(proposedCell);
  }

  REAL8 logCurrentCellFactor = log((REAL8)currentCell->npts / npts);
  REAL8 logProposedCellFactor = log((REAL8)proposedCell->npts / npts);

  return logCurrentCellFactor + logCurrentVolume - logProposedCellFactor - logProposedVolume;
}

UINT4 LALInferenceCheckPositiveDefinite( 
                          gsl_matrix       *matrix,
                          UINT4            dim
                          )
{
    gsl_matrix  *m     = NULL;
    gsl_vector  *eigen = NULL;
    gsl_eigen_symm_workspace *workspace = NULL;
    UINT4 i;
    
    /* copy input matrix */
    m =  gsl_matrix_alloc( dim,dim ); 
    gsl_matrix_memcpy( m, matrix);  
    
    /* prepare variables */
    eigen = gsl_vector_alloc ( dim );
    workspace = gsl_eigen_symm_alloc ( dim );
    
    /* compute the eigen values */
    gsl_eigen_symm ( m,  eigen, workspace );
    
    /* test the result */
    for (i = 0; i < dim; i++)
    {
        /* printf("diag: %f | eigen[%d]= %f\n", gsl_matrix_get( matrix,i,i), i, eigen->data[i]);*/
        if (eigen->data[i]<0) 
        {
            printf("NEGATIVE EIGEN VALUE!!! PANIC\n");
            return 0;
        }
    }
    
    /* freeing unused stuff */
    gsl_eigen_symm_free( workspace);
    gsl_matrix_free(m);
    gsl_vector_free(eigen);
    
    return 1;
}

/* Reference: http://www.mail-archive.com/help-gsl@gnu.org/msg00631.html*/
void
XLALMultiNormalDeviates(
                        REAL4Vector *vector,
                        gsl_matrix *matrix,
                        UINT4 dim,
                        RandomParams *randParam
                        )
{
    UINT4 i=0;
    gsl_matrix *work=NULL;
    gsl_vector *result = NULL;

    /* check input arguments */
    if (!vector || !matrix || !randParam)
        XLAL_ERROR_VOID( XLAL_EFAULT );

    if (dim<1)
        XLAL_ERROR_VOID( XLAL_EINVAL );

    /* copy matrix into workspace */
    work =  gsl_matrix_alloc(dim,dim);
    gsl_matrix_memcpy( work, matrix );

    /* compute the cholesky decomposition */
    gsl_linalg_cholesky_decomp(work);

    /* retrieve the normal distributed random numbers (LAL procedure) */
    XLALNormalDeviates( vector, randParam );

    /* store this into a gsl vector */
    result = gsl_vector_alloc ( (int)dim );
    for (i = 0; i < dim; i++)
    {
        gsl_vector_set (result, i, vector->data[i]);
    }

    /* compute the matrix-vector multiplication */
    gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, result);

    /* recopy the results */
    for (i = 0; i < dim; i++)
    {
        vector->data[i]=gsl_vector_get (result, i);
    }

    /* free unused stuff */
    gsl_matrix_free(work);
    gsl_vector_free(result);

}

void
XLALMultiStudentDeviates(
                         REAL4Vector  *vector,
                         gsl_matrix   *matrix,
                         UINT4         dim,
                         UINT4         n,
                         RandomParams *randParam
                         )
{
    REAL4Vector *dummy=NULL;
    REAL4 chi=0.0, factor;
    UINT4 i;

    /* check input arguments */
    if (!vector || !matrix || !randParam)
        XLAL_ERROR_VOID( XLAL_EFAULT );

    if (dim<1)
        XLAL_ERROR_VOID( XLAL_EINVAL );

    if (n<1)
        XLAL_ERROR_VOID( XLAL_EINVAL );


    /* first draw from MVN */
        XLALMultiNormalDeviates( vector, matrix, dim, randParam);

    /* then draw from chi-square with n degrees of freedom;
     this is the sum d_i*d_i with d_i drawn from a normal 
     distribution. */
        dummy = XLALCreateREAL4Vector( n );
        XLALNormalDeviates( dummy, randParam );

    /* calculate the chisquare distributed value */
    for (i=0; i<n; i++)
    {
        chi+=dummy->data[i]*dummy->data[i];
    }

    /* destroy the helping vector */
        XLALDestroyREAL4Vector( dummy );

    /* now, finally, calculate the distribution value */
    factor=sqrt(n/chi);
    for (i=0; i<dim; i++)
    {
        vector->data[i]*=factor;
    }

}

/* Calculate shortest angular distance between a1 and a2 */
REAL8 LALInferenceAngularDistance(REAL8 a1, REAL8 a2){
    double raw = (a2>a1 ? a2-a1 : a1-a2);
    return(raw>LAL_PI ? 2.0*LAL_PI - raw : raw);
}

/* Calculate the variance of a modulo-2pi distribution */
REAL8 LALInferenceAngularVariance(LALInferenceVariables **list,const char *pname, int N){
        int i=0;
        REAL8 ang_mean=0.0;
        REAL8 var=0.0;
        REAL8 ms,mc;
        /* Calc mean */
        for(i=0,ms=0.0,mc=0.0;i<N;i++) {
                ms+=sin(*(REAL8 *)LALInferenceGetVariable(list[i],pname));
                mc+=cos(*(REAL8 *)LALInferenceGetVariable(list[i],pname));
        }
        ms/=N; mc/=N;
        ang_mean=atan2(ms,mc);
        ang_mean = ang_mean<0? 2.0*LAL_PI + ang_mean : ang_mean;
        /* calc variance */
        for(i=0;i<N;i++) var+=LALInferenceAngularDistance(*(REAL8 *)LALInferenceGetVariable(list[i],pname),ang_mean)*LALInferenceAngularDistance(*(REAL8 *)LALInferenceGetVariable(list[i],pname),ang_mean);
        return(var/(REAL8)N);
}

/* Sanity check the data structures and print any encountered errors */
INT4 LALInferenceSanityCheck(LALInferenceRunState *state)
{
  INT4 retcode=0;
  if(!state) {
	fprintf(stderr,"NULL state pointer!\n");
	return(1);
  }
  
  LALInferenceIFOData *data=state->data;
  if(!data) {
	fprintf(stderr,"NULL data pointer!\n");
        return(1);
  }
  while(data){
    retcode=0;
    fprintf(stderr,"Checking %s:\n",data->name);
    if(data->timeData) {
      fprintf(stderr,"Checking timeData: ");
      if(!(retcode|=checkREAL8TimeSeries(data->timeData))) fprintf(stderr," OK\n");
    }
    if(data->whiteTimeData) {
      fprintf(stderr,"Checking whiteTimeData: ");
      if(!(retcode|=checkREAL8TimeSeries(data->whiteTimeData))) fprintf(stderr," OK\n");
    }
    if(data->windowedTimeData) {
      fprintf(stderr,"Checking windowedTimeData: ");
      if(!(retcode|=checkREAL8TimeSeries(data->windowedTimeData))) fprintf(stderr," OK\n");
    }
    if(data->freqData) {
      fprintf(stderr,"Checking freqData: ");
      if(!(retcode|=checkCOMPLEX16FrequencySeries(data->freqData))) fprintf(stderr," OK\n");
    }
    if(data->whiteFreqData) {
      fprintf(stderr,"Checking whiteFreqData: ");
      if(!(retcode|=checkCOMPLEX16FrequencySeries(data->whiteFreqData))) fprintf(stderr," OK\n");
    }
    if(data->oneSidedNoisePowerSpectrum) {
      fprintf(stderr,"Checking oneSidedNoisePowerSpectrum: ");
      if(!(retcode|=checkREAL8FrequencySeries(data->oneSidedNoisePowerSpectrum))) fprintf(stderr," OK\n");
    }
    data=data->next;
  }

  LALInferenceModel *model = state->model;
  if(!model) {
	fprintf(stderr,"NULL model pointer!\n");
        return(1);
  }
  if(model->timehPlus) {
    fprintf(stderr,"Checking timehPlus: ");
    if(!(retcode|=checkREAL8TimeSeries(model->timehPlus))) fprintf(stderr," OK\n");
  }
  if(model->timehCross) {
    fprintf(stderr,"Checking timehCross: ");
    if(!(retcode|=checkREAL8TimeSeries(model->timehCross))) fprintf(stderr," OK\n");
  }
  if(model->freqhPlus) {
    fprintf(stderr,"Checking freqhPlus: ");
    if(!(retcode|=checkCOMPLEX16FrequencySeries(model->freqhPlus))) fprintf(stderr," OK\n");
  }
  if(model->freqhCross) {
    fprintf(stderr,"Checking freqhCross: ");
    if(!(retcode|=checkCOMPLEX16FrequencySeries(model->freqhCross))) fprintf(stderr," OK\n");
  }

  return(retcode);
}

static INT4 checkREAL8Value(REAL8 val);

static INT4 checkREAL8TimeSeries(REAL8TimeSeries *series)
{
  UINT4 i;
  INT4 retcode=0;
  if(!series) fprintf(stderr,"Null REAL8TimeSeries *\n");
  else {
    if(!series->data) fprintf(stderr,"NULL REAL8Sequence structure in REAL8TimeSeries\n");
    else {
      if(!series->data->data) fprintf(stderr,"NULL REAL8[] in REAL8Sequence\n");
      else {
       for(i=0;i<series->data->length;i++) {if(checkREAL8Value(series->data->data[i])) {if(!retcode) fprintf(stderr,"Found value %lf at index %i\n",series->data->data[i],i); retcode+=1;}}
      }
    }
  }
  if(retcode>1) fprintf(stderr,"and %i more times\n",retcode);
  return(retcode);
}

static INT4 checkREAL8FrequencySeries(REAL8FrequencySeries *series)
{
  UINT4 i;
  INT4 retcode=0;
  if(!series) fprintf(stderr,"Null REAL8FrequencySeries *\n");
  else {
    if(!series->data) fprintf(stderr,"NULL REAL8Sequence structure in REAL8FrequencySeries\n");
    else {
      if(!series->data->data) fprintf(stderr,"NULL REAL8[] in REAL8Sequence\n");
      else {
       for(i=0;i<series->data->length;i++) {if(checkREAL8Value(series->data->data[i])) {if(!retcode) fprintf(stderr,"Found value %lf at index %i\n",series->data->data[i],i); retcode+=1;}}
      }
    }
  }
  if(retcode>1) fprintf(stderr,"and %i more times\n",retcode);
  return(retcode);
}

static INT4 checkCOMPLEX16FrequencySeries(COMPLEX16FrequencySeries *series)
{
  UINT4 i;
  INT4 retcode=0;
  if(!series) fprintf(stderr,"Null COMPLEX16FrequencySeries *\n");
  else {
    if(!series->data) fprintf(stderr,"NULL COMPLEX16Sequence structure in COMPLEX16FrequencySeries\n");
    else {
      if(!series->data->data) fprintf(stderr,"NULL REAL8[] in COMPLEX16Sequence\n");
      else {
       for(i=0;i<series->data->length;i++) {if(checkREAL8Value(creal(series->data->data[i]))) {if(!retcode) fprintf(stderr,"Found real value %lf at index %i\n",creal(series->data->data[i]),i); retcode+=1;}
					if(checkREAL8Value(cimag(series->data->data[i]))) {if(!retcode) fprintf(stderr,"Found imag value %lf at index %i\n",cimag(series->data->data[i]),i); retcode+=1;} }
      }
    }
  }
  if(retcode>1) fprintf(stderr,"and %i more times\n",retcode);
  return(retcode);
}

static INT4 checkREAL8Value(REAL8 val)
{
  if(isinf(val)) return 1;
  if(isnan(val)) return 1;
  return 0;
}

void LALInferenceDumpWaveforms(LALInferenceModel *model, const char *basefilename)
{
    UINT4 i;
    FILE *dumpfile=NULL;
    char basename[1024]="template_dump";
    char filename[1024]="";
    if(basefilename!=NULL)
    {
        sprintf(basename,"%s",basefilename);
    }

    if(model->timehPlus && model->timehCross){
        sprintf(filename,"%s_time.txt",basename);
        dumpfile=fopen(filename,"w");
        REAL8 epoch = model->timehPlus->epoch.gpsSeconds + 1e-9*model->timehPlus->epoch.gpsNanoSeconds;
        REAL8 dt=model->timehPlus->deltaT;
        for(i=0;i<model->timehPlus->data->length;i++) fprintf(dumpfile,"%10.20e %10.20e %10.20e\n",epoch+i*dt,model->timehPlus->data->data[i],model->timehCross->data->data[i]);
        fclose(dumpfile);
        fprintf(stdout,"Dumped file %s\n",filename);
    }
    if(model->freqhPlus && model->freqhCross){
        sprintf(filename,"%s_freq.txt",basename);
        dumpfile=fopen(filename,"w");
        REAL8 fLow=model->freqhPlus->f0;
        REAL8 df=model->freqhPlus->deltaF;
        for(i=0;i<model->freqhPlus->data->length;i++) fprintf(dumpfile,"%10.20e %10.20e %10.20e %10.20e %10.20e\n",fLow+i*df,creal(model->freqhPlus->data->data[i]), cimag(model->freqhPlus->data->data[i]),creal(model->freqhCross->data->data[i]), cimag(model->freqhCross->data->data[i]));
        fclose(dumpfile);
        fprintf(stdout,"Dumped file %s\n",filename);
    }
}

static void REAL8Vector_fwrite(FILE *f, REAL8Vector *vec);
static void REAL8Vector_fwrite(FILE *f, REAL8Vector *vec)
{
  fwrite(&(vec->length),sizeof(UINT4),1,f);
  fwrite(vec->data,sizeof(REAL8),vec->length,f);
}

static void UINT4Vector_fwrite(FILE *f, UINT4Vector *vec);
static void UINT4Vector_fwrite(FILE *f, UINT4Vector *vec) 
{
  fwrite(&(vec->length),sizeof(vec->length),1,f);
  fwrite(vec->data,sizeof(vec->data[0]),vec->length,f);
}

static REAL8Vector * REAL8Vector_fread(FILE *f);
static REAL8Vector * REAL8Vector_fread(FILE *f)
{
  REAL8Vector *out=NULL;
  UINT4 size;
  fread(&size,sizeof(size),1,f);
  out=XLALCreateREAL8Vector(size);
  fread(out->data,sizeof(REAL8),size,f);
  return out;
}

static UINT4Vector * UINT4Vector_fread(FILE *f);
static UINT4Vector * UINT4Vector_fread(FILE *f)
{
  UINT4 size = 0;
  UINT4Vector *vec = NULL;

  fread(&size,sizeof(size),1,f);
  vec = XLALCreateUINT4Vector(size);
  fread(vec->data, sizeof(UINT4), size, f);
  return vec;
}

int LALInferenceWriteVariablesBinary(FILE *file, LALInferenceVariables *vars)
{
  int i=0;
  char termchar='\n';
  if(!vars) return -1;
  LALInferenceVariableItem *item=vars->head;
  /* Write initial info (number of dimensions) */
  fwrite(&(vars->dimension),sizeof(vars->dimension),1,file);
  /* Write each item */
  for(i=0;item;i++)
  {
    /* Name */
    fputs(item->name,file);
    fwrite(&termchar,sizeof(char),1,file);
    fwrite(&(item->type),sizeof(item->type),1,file);
    fwrite(&(item->vary),sizeof(item->vary),1,file);
    switch(item->type)
    {
      case LALINFERENCE_gslMatrix_t:
      {
	gsl_matrix *matrix=*(gsl_matrix **)item->value;
	fwrite(&(matrix->size1),sizeof(matrix->size1),1,file);
	fwrite(&(matrix->size2),sizeof(matrix->size2),1,file);
	gsl_matrix_fwrite(file,matrix);
	break;
      }
      case LALINFERENCE_REAL8Vector_t:
      {
	REAL8Vector *vec=*(REAL8Vector **)item->value;
	REAL8Vector_fwrite(file,vec);
	break;
      }
    case LALINFERENCE_UINT4Vector_t:
      {
	UINT4Vector *vec = *(UINT4Vector **)item->value;
	UINT4Vector_fwrite(file, vec);
	break;
      }
    case LALINFERENCE_string_t:
      {
	char *value = *((char **)item->value);
	size_t len = strlen(value);
	fwrite(&len, sizeof(size_t),1, file);
	fwrite(value, sizeof(char), len, file);
	break;
      }
    case LALINFERENCE_MCMCrunphase_ptr_t:
      {
	LALInferenceMCMCRunPhase *ph = *((LALInferenceMCMCRunPhase **)item->value);
	fwrite(ph, sizeof(LALInferenceMCMCRunPhase), 1, file);
	break;
      }
    case LALINFERENCE_void_ptr_t:
      {
	/* Write void_ptr as NULL, so fails if used without
	   initialization on restart. */
	void *out = NULL;
	fwrite(&out,sizeof(void*),1,file);
	break;
      }
      default:
      {
	fwrite(item->value,LALInferenceTypeSize[item->type],1,file);
	break;
      }
    }
    item=item->next;
  }
  return i;
}

LALInferenceVariables *LALInferenceReadVariablesBinary(FILE *stream)
{
  UINT4 j;
  UINT4 dim;
  LALInferenceVariables *vars=XLALCalloc(1,sizeof(LALInferenceVariables));

  /* Number of variables to read */
  fread(&dim, sizeof(vars->dimension), 1, stream);

  /* Now read them in */
  for(;dim>0;dim--)
  {
    char name[VARNAME_MAX];
    LALInferenceVariableType type;
    LALInferenceParamVaryType vary;
    fgets(name,sizeof(name),stream);
    
    for(j=0;j<sizeof(name);j++) if(name[j]=='\n') {name[j]='\0'; break;}
    if(j==sizeof(name))
    {
      fprintf(stderr,"ERROR reading saved variable!");
      return(NULL);
    }
    fread(&type,sizeof(type),1,stream);
    fread(&vary,sizeof(vary),1,stream);
    switch(type)
    {
      case  LALINFERENCE_gslMatrix_t:
      {
	size_t size1,size2;
	fread(&size1,sizeof(size1),1,stream);
	fread(&size2,sizeof(size2),1,stream);
	gsl_matrix *matrix=gsl_matrix_alloc(size1,size2);
	gsl_matrix_fread(stream,matrix);
	LALInferenceAddVariable(vars,name,&matrix,type,vary);

	break;
      }
      case LALINFERENCE_REAL8Vector_t:
      {
	REAL8Vector *v=REAL8Vector_fread(stream);
	LALInferenceAddVariable(vars,name,&v,type,vary);

	break;
      }
    case LALINFERENCE_UINT4Vector_t:
      {
	UINT4Vector *vec = UINT4Vector_fread(stream);
	LALInferenceAddVariable(vars,name,&vec,type,vary);
	break;
      }
    case LALINFERENCE_string_t:
      {
	size_t len = 0;
	char *string = NULL;

	fread(&len, sizeof(size_t), 1, stream);
	string = XLALCalloc(sizeof(char), len+1); /* One extra character: '\0' */
	fread(string, sizeof(char), len, stream);
	LALInferenceAddVariable(vars,name,&string,type,vary);
      }
    case LALINFERENCE_MCMCrunphase_ptr_t:
      {
	LALInferenceMCMCRunPhase *ph = XLALCalloc(sizeof(LALInferenceMCMCRunPhase),1);
	fread(ph, sizeof(LALInferenceMCMCRunPhase), 1, stream);
	LALInferenceAddVariable(vars,name,&ph,type,vary);
	break;
      }
    case LALINFERENCE_void_ptr_t:
      {
	void *ptr = NULL;
	fread(&ptr,sizeof(void *), 1, stream);
	LALInferenceAddVariable(vars,name,&ptr,type,vary);
	break;
      }
      default:
      {
	void *value=NULL;
	UINT4 storagesize=LALInferenceTypeSize[type];
	value=XLALCalloc(1,storagesize);
	fread(value, storagesize, 1, stream);
	LALInferenceAddVariable(vars,name,value,type,vary);
      }
    }
  }
  return vars;
}

int LALInferenceWriteVariablesArrayBinary(FILE *file, LALInferenceVariables **vars, UINT4 N)
{
  UINT4 i=0;
  for(i=0;i<N;i++) LALInferenceWriteVariablesBinary(file, vars[i]);
  return N;
}

int LALInferenceReadVariablesArrayBinary(FILE *file, LALInferenceVariables **vars, UINT4 N)
{
  UINT4 i=0;
  for(i=0;i<N;i++){
    vars[i]=LALInferenceReadVariablesBinary(file);
  }
  return N;
}

int LALInferenceWriteRunStateBinary(FILE *file, LALInferenceRunState *runState)
{
  int flag=0;
  fwrite(&(runState->differentialPointsLength),sizeof(runState->differentialPointsLength),1,file);
  fwrite(&(runState->differentialPointsSize),sizeof(runState->differentialPointsSize),1,file);
  fwrite(&(runState->currentLikelihood),sizeof(runState->currentLikelihood),1,file);
  fwrite(&(runState->currentPrior),sizeof(runState->currentPrior),1,file);
  flag|=gsl_rng_fwrite (file , runState->GSLrandom);
  flag|=LALInferenceWriteVariablesBinary(file, runState->currentParams);
  flag|=LALInferenceWriteVariablesBinary(file, runState->priorArgs);
  flag|=LALInferenceWriteVariablesBinary(file, runState->proposalArgs);
  flag|=LALInferenceWriteVariablesBinary(file, runState->proposalStats);
  flag|=LALInferenceWriteVariablesBinary(file, runState->algorithmParams);
  // Currently live points are the same as differential points buffer
  //flag|=LALInferenceWriteVariablesArrayBinary(file, state->livePoints, state->differentialPointsLength);
  flag|=LALInferenceWriteVariablesArrayBinary(file, runState->differentialPoints, runState->differentialPointsLength);
  return flag;
}

int LALInferenceReadRunStateBinary(FILE *file, LALInferenceRunState *runState)
{
  fread(&(runState->differentialPointsLength),sizeof(runState->differentialPointsLength),1,file);
  fread(&(runState->differentialPointsSize),sizeof(runState->differentialPointsSize),1,file);
  runState->differentialPoints=XLALCalloc(runState->differentialPointsSize,sizeof(LALInferenceVariables *));
  fread(&(runState->currentLikelihood),sizeof(runState->currentLikelihood),1,file);
  fread(&(runState->currentPrior),sizeof(runState->currentPrior),1,file);
  
  gsl_rng_fread(file, runState->GSLrandom);
  runState->currentParams=LALInferenceReadVariablesBinary(file);
  runState->priorArgs=LALInferenceReadVariablesBinary(file);
  runState->proposalArgs=LALInferenceReadVariablesBinary(file);
  runState->proposalStats=LALInferenceReadVariablesBinary(file);
  runState->algorithmParams=LALInferenceReadVariablesBinary(file);
  LALInferenceReadVariablesArrayBinary(file, runState->differentialPoints,runState->differentialPointsLength);
  runState->livePoints=runState->differentialPoints;
  
  return 0;
}

void LALInferenceAddINT4Variable(LALInferenceVariables * vars, const char * name, INT4 value, LALInferenceParamVaryType vary)
/* Typed version of LALInferenceAddVariable for INT4 values.*/
{
  LALInferenceAddVariable(vars,name,(void*)&value,LALINFERENCE_INT4_t,vary);
}

INT4 LALInferenceGetINT4Variable(LALInferenceVariables * vars, const char * name)
/* Typed version of LALInferenceGetVariable for INT4 values.*/
{

  if(LALInferenceGetVariableType(vars,name)!=LALINFERENCE_INT4_t){
    XLAL_ERROR(XLAL_ETYPE);
  }

  INT4* rvalue=(INT4*)LALInferenceGetVariable(vars,name);

  return *rvalue;
}

void LALInferenceSetINT4Variable(LALInferenceVariables* vars,const char* name,INT4 value){
  LALInferenceSetVariable(vars,name,(void*)&value);
}

void LALInferenceAddINT8Variable(LALInferenceVariables * vars, const char * name, INT8 value, LALInferenceParamVaryType vary)
/* Typed version of LALInferenceAddVariable for INT8 values.*/
{
  LALInferenceAddVariable(vars,name,(void*)&value,LALINFERENCE_INT8_t,vary);
}

INT8 LALInferenceGetINT8Variable(LALInferenceVariables * vars, const char * name)
/* Typed version of LALInferenceGetVariable for INT8 values.*/
{

  if(LALInferenceGetVariableType(vars,name)!=LALINFERENCE_INT8_t){
    XLAL_ERROR(XLAL_ETYPE);
  }

  INT8* rvalue=(INT8*)LALInferenceGetVariable(vars,name);

  return *rvalue;
}

void LALInferenceSetINT8Variable(LALInferenceVariables* vars,const char* name,INT8 value){
  LALInferenceSetVariable(vars,name,(void*)&value);
}

void LALInferenceAddUINT4Variable(LALInferenceVariables * vars, const char * name, UINT4 value, LALInferenceParamVaryType vary)
/* Typed version of LALInferenceAddVariable for UINT4 values.*/
{
  LALInferenceAddVariable(vars,name,(void*)&value,LALINFERENCE_UINT4_t,vary);
}

UINT4 LALInferenceGetUINT4Variable(LALInferenceVariables * vars, const char * name)
/* Typed version of LALInferenceGetVariable for UINT4 values.*/
{

  if(LALInferenceGetVariableType(vars,name)!=LALINFERENCE_UINT4_t){
    XLAL_ERROR(XLAL_ETYPE);
  }

  UINT4* rvalue=(UINT4*)LALInferenceGetVariable(vars,name);

  return *rvalue;
}

void LALInferenceSetUINT4Variable(LALInferenceVariables* vars,const char* name,UINT4 value){
  LALInferenceSetVariable(vars,name,(void*)&value);
}

void LALInferenceAddREAL4Variable(LALInferenceVariables * vars, const char * name, REAL4 value, LALInferenceParamVaryType vary)
/* Typed version of LALInferenceAddVariable for REAL4 values.*/
{
  LALInferenceAddVariable(vars,name,(void*)&value,LALINFERENCE_REAL4_t,vary);
}

REAL4 LALInferenceGetREAL4Variable(LALInferenceVariables * vars, const char * name)
/* Typed version of LALInferenceGetVariable for REAL4 values.*/
{

  if(LALInferenceGetVariableType(vars,name)!=LALINFERENCE_REAL4_t){
    XLAL_ERROR_REAL4(XLAL_ETYPE);
  }

  REAL4* rvalue=(REAL4*)LALInferenceGetVariable(vars,name);

  return *rvalue;
}

void LALInferenceSetREAL4Variable(LALInferenceVariables* vars,const char* name,REAL4 value){
  LALInferenceSetVariable(vars,name,(void*)&value);
}

void LALInferenceAddREAL8Variable(LALInferenceVariables * vars, const char * name, REAL8 value, LALInferenceParamVaryType vary)
/* Typed version of LALInferenceAddVariable for REAL8 values.*/
{
  LALInferenceAddVariable(vars,name,(void*)&value,LALINFERENCE_REAL8_t,vary);
}

REAL8 LALInferenceGetREAL8Variable(LALInferenceVariables * vars, const char * name)
/* Typed version of LALInferenceGetVariable for REAL8 values.*/
{

  if(LALInferenceGetVariableType(vars,name)!=LALINFERENCE_REAL8_t){
    XLAL_ERROR_REAL8(XLAL_ETYPE);
  }

  REAL8* rvalue=(REAL8*)LALInferenceGetVariable(vars,name);

  return *rvalue;
}

void LALInferenceSetREAL8Variable(LALInferenceVariables* vars,const char* name,REAL8 value){
  LALInferenceSetVariable(vars,name,(void*)&value);
}

void LALInferenceAddCOMPLEX8Variable(LALInferenceVariables * vars, const char * name, COMPLEX8 value, LALInferenceParamVaryType vary)
/* Typed version of LALInferenceAddVariable for COMPLEX8 values.*/
{
  LALInferenceAddVariable(vars,name,(void*)&value,LALINFERENCE_COMPLEX8_t,vary);
}

COMPLEX8 LALInferenceGetCOMPLEX8Variable(LALInferenceVariables * vars, const char * name)
/* Typed version of LALInferenceGetVariable for COMPLEX8 values.*/
{

  COMPLEX8* rvalue=(COMPLEX8*)LALInferenceGetVariable(vars,name);

  return *rvalue;
}

void LALInferenceSetCOMPLEX8Variable(LALInferenceVariables* vars,const char* name,COMPLEX8 value){
  LALInferenceSetVariable(vars,name,(void*)&value);
}

void LALInferenceAddCOMPLEX16Variable(LALInferenceVariables * vars, const char * name, COMPLEX16 value, LALInferenceParamVaryType vary)
/* Typed version of LALInferenceAddVariable for COMPLEX16 values.*/
{
  LALInferenceAddVariable(vars,name,(void*)&value,LALINFERENCE_COMPLEX16_t,vary);
}

COMPLEX16 LALInferenceGetCOMPLEX16Variable(LALInferenceVariables * vars, const char * name)
/* Typed version of LALInferenceGetVariable for COMPLEX16 values.*/
{

  COMPLEX16* rvalue=(COMPLEX16*)LALInferenceGetVariable(vars,name);

  return *rvalue;
}

void LALInferenceSetCOMPLEX16Variable(LALInferenceVariables* vars,const char* name,COMPLEX16 value){
  LALInferenceSetVariable(vars,name,(void*)&value);
}

void LALInferenceAddgslMatrixVariable(LALInferenceVariables * vars, const char * name, gsl_matrix* value, LALInferenceParamVaryType vary)
/* Typed version of LALInferenceAddVariable for gsl_matrix values.*/
{
  LALInferenceAddVariable(vars,name,(void*)value,LALINFERENCE_gslMatrix_t,vary);
}

gsl_matrix* LALInferenceGetgslMatrixVariable(LALInferenceVariables * vars, const char * name)
/* Typed version of LALInferenceGetVariable for gsl_matrix values.*/
{

  if(LALInferenceGetVariableType(vars,name)!=LALINFERENCE_gslMatrix_t){
    XLAL_ERROR_NULL(XLAL_ETYPE);
  }

  gsl_matrix* rvalue=(gsl_matrix*)LALInferenceGetVariable(vars,name);

  return rvalue;
}

void LALInferenceSetgslMatrixVariable(LALInferenceVariables* vars,const char* name,gsl_matrix* value){
  LALInferenceSetVariable(vars,name,(void*)&value);
}

void LALInferenceAddREAL8VectorVariable(LALInferenceVariables * vars, const char * name, REAL8Vector* value, LALInferenceParamVaryType vary)
/* Typed version of LALInferenceAddVariable for REAL8Vector values.*/
{
  LALInferenceAddVariable(vars,name,(void*)value,LALINFERENCE_REAL8Vector_t,vary);
}

REAL8Vector* LALInferenceGetREAL8VectorVariable(LALInferenceVariables * vars, const char * name)
/* Typed version of LALInferenceGetVariable for REAL8Vector values.*/
{

  if(LALInferenceGetVariableType(vars,name)!=LALINFERENCE_REAL8Vector_t){
    XLAL_ERROR_NULL(XLAL_ETYPE);
  }

  REAL8Vector* rvalue=(REAL8Vector*)LALInferenceGetVariable(vars,name);

  return rvalue;
}

void LALInferenceSetREAL8VectorVariable(LALInferenceVariables* vars,const char* name,REAL8Vector* value){
  LALInferenceSetVariable(vars,name,(void*)&value);
}

void LALInferenceAddUINT4VectorVariable(LALInferenceVariables * vars, const char * name, UINT4Vector* value, LALInferenceParamVaryType vary)
/* Typed version of LALInferenceAddVariable for UINT4Vector values.*/
{
  LALInferenceAddVariable(vars,name,(void*)value,LALINFERENCE_UINT4Vector_t,vary);
}

UINT4Vector* LALInferenceGetUINT4VectorVariable(LALInferenceVariables * vars, const char * name)
/* Typed version of LALInferenceGetVariable for UINT4Vector values.*/
{

  if(LALInferenceGetVariableType(vars,name)!=LALINFERENCE_UINT4Vector_t){
    XLAL_ERROR_NULL(XLAL_ETYPE);
  }

  UINT4Vector* rvalue=(UINT4Vector*)LALInferenceGetVariable(vars,name);

  return rvalue;
}

void LALInferenceSetUINT4VectorVariable(LALInferenceVariables* vars,const char* name,UINT4Vector* value){
  LALInferenceSetVariable(vars,name,(void*)&value);
}

void LALInferenceAddMCMCrunphase_ptrVariable(LALInferenceVariables * vars, const char * name, LALInferenceMCMCRunPhase* value, LALInferenceParamVaryType vary)
/* Typed version of LALInferenceAddVariable for LALInferenceMCMCRunPhase values.*/
{
  LALInferenceAddVariable(vars,name,(void*)value,LALINFERENCE_MCMCrunphase_ptr_t,vary);
}

LALInferenceMCMCRunPhase* LALInferenceGetMCMCrunphase_ptrVariable(LALInferenceVariables * vars, const char * name)
/* Typed version of LALInferenceGetVariable for LALInferenceMCMCRunPhase values.*/
{

  if(LALInferenceGetVariableType(vars,name)!=LALINFERENCE_MCMCrunphase_ptr_t){
    XLAL_ERROR_NULL(XLAL_ETYPE);
  }

  LALInferenceMCMCRunPhase* rvalue=(LALInferenceMCMCRunPhase*)LALInferenceGetVariable(vars,name);

  return rvalue;
}

void LALInferenceSetMCMCrunphase_ptrVariable(LALInferenceVariables* vars,const char* name,LALInferenceMCMCRunPhase* value){
  LALInferenceSetVariable(vars,name,(void*)&value);
}

void LALInferenceAddstringVariable(LALInferenceVariables * vars, const char * name, CHAR* value, LALInferenceParamVaryType vary)
/* Typed version of LALInferenceAddVariable for CHAR values.*/
{
  LALInferenceAddVariable(vars,name,(void*)value,LALINFERENCE_string_t,vary);
}

CHAR* LALInferenceGetstringVariable(LALInferenceVariables * vars, const char * name)
/* Typed version of LALInferenceGetVariable for CHAR values.*/
{

  if(LALInferenceGetVariableType(vars,name)!=LALINFERENCE_string_t){
    XLAL_ERROR_NULL(XLAL_ETYPE);
  }

  CHAR* rvalue=(CHAR*)LALInferenceGetVariable(vars,name);

  return rvalue;
}

void LALInferenceSetstringVariable(LALInferenceVariables* vars,const char* name,CHAR* value){
  LALInferenceSetVariable(vars,name,(void*)&value);
}

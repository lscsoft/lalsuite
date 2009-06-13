/*
*  Copyright (C) 2007 Bernd Machenschalk, Philip Charlton
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

/* -*- mode: c; c-basic-offset: 2; -*- */

/*
  This is the source code for the official implementation of the
  Fast Chirp transform.
*/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>

#include <math.h>

#include "fct_fft.h"
#include "fct.h"

#include <lal/LALRCSID.h>
NRCSID (FCTC,"$Id$");

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795029L
#endif

fct_malloc_type_function fct_malloc_hook = 0;
fct_calloc_type_function fct_calloc_hook = 0;
fct_free_type_function fct_free_hook = 0;

/*Function implementations*/
STORAGE_CLASS
void* fct_malloc(size_t size)
{
  if (fct_malloc_hook != 0)
  {
    return fct_malloc_hook(size);
  }

  return malloc(size);
}

STORAGE_CLASS
void* fct_calloc(size_t nmemb, size_t size)
{
  if (fct_calloc_hook != 0)
  {
    return fct_calloc_hook(nmemb, size);
  }

  return calloc(nmemb, size);
}

STORAGE_CLASS
void fct_free(void* p)
{
  if (fct_free_hook != 0)
  {
    fct_free_hook(p);
    return;
  }

  free(p);
}

static
const char** setup_errstr(void)
{
  static const char* errstr[FCT_EUNKNOWN + 1];

  errstr[FCT_ENONE]      = "Not an error";
  errstr[FCT_ENULL_PLAN] = "The fct_plan pointer is null";
  errstr[FCT_ENOT_PWR2]  = "The fct data length is not a power of 2";
  errstr[FCT_EMEM]       = "Out of memory";
  errstr[FCT_EFFTW_PLAN] = "Unable to create FFTW plan";
  errstr[FCT_EDIM]       = "Dimension must be > 0 and < number of dimensions";
  errstr[FCT_ENULL_LOC]  = "A pointer to an array of start or end locations "
    "was null";

  errstr[FCT_ENONNULL_MASK]   = "Mask is unimplemented";
  errstr[FCT_ENULL_DATACUBE]  = "A pointer to a data cube is null";
  errstr[FCT_EDATACUBE_MODE]  = "Invalid mode for adding a data cube";
  errstr[FCT_EDATACUBE_RANGE] = "Invalid range for data cube";

  errstr[FCT_EOFAC] = "Oversampling factor must be != 0";

  errstr[FCT_EUNKNOWN] = "Unknown error";

  return &errstr[0];
}

STORAGE_CLASS
const char* fct_strerror(int fct_errno)
{
  static const char** errstr = 0;

  if (errstr == 0)
  {
    errstr = setup_errstr();
  }

  if ((fct_errno < 0) || (fct_errno > FCT_EUNKNOWN))
  {
    fct_errno = FCT_EUNKNOWN;
  }

  return errstr[fct_errno];
}

STORAGE_CLASS
fct_plan *fct_init_plan(int data_length,int number_of_dimensions,
			int dimension_0_stride, fct_status* const status)
/*This function partially initializes the FCT. It must be called before
  any other fct_ function. It generates an fct_plan structure which tells
  the other fct_ functions how to operate on the data. After calling
  fct_init_plan, the user must specify both the phase functions using
  fct_set_phase_function and the FCT parameters over which to calculate
  the FCT using fct_add_indicies.*/
{
  int k = 0;
  fct_plan *plan = 0;

  /*Allocate memory for the plan structure*/
  plan = (fct_plan *)fct_malloc(sizeof(fct_plan));

  /*Check to see if allocation was successful. Exit otherwise.*/
  if (plan == NULL) {
    FCT_ERROR(status, FCT_EMEM);
    return 0;
  }

  /*Set plan structure values*/
  plan->data_length             = data_length;
  plan->number_of_dimensions    = number_of_dimensions;
  plan->dimension_0_stride      = dimension_0_stride;
  plan->offset                  = 0;
  plan->dt                      = 1;
  plan->max_number_of_segments  = 1;
  plan->parameter_list_start    = (data_cube *)NULL;
  plan->parameter_list_end      = (data_cube *)NULL;

  /*Set up the FFTs*/
  fct_setup_fft(plan, status);
  if (status->fct_errno != 0)
  {
    return 0;
  }

  /*Allocate memory for the phase_func function array*/
  plan->phase_functions
    = fct_calloc(number_of_dimensions, sizeof(function_pointer));

  /*Check to make sure the allocation worked.*/
  if (plan->phase_functions == 0) {
    FCT_ERROR(status, FCT_EMEM);
    return 0;
  }

  /* Allocate memory for the phase_func oversampling factor array */
  plan->ofac = fct_calloc(number_of_dimensions, sizeof(*(plan->ofac)));

  /*Check to make sure the allocation worked.*/
  if (plan->ofac == 0) {
    FCT_ERROR(status, FCT_EMEM);
    return 0;
  }

  /* Initialise */
  for (k = 0; k < number_of_dimensions; ++k)
  {
    plan->ofac[k] = 1;
  }

  return plan;
}

STORAGE_CLASS
void fct_set_units(fct_plan *plan, float offset, float delta,
		   fct_status* const status)
/*Sets the conversion from data sample to the argument of the
  phase function. For the jth data point, the phase functions are
  evaluated at offset + j*delta. The default values are offset = 0 and
  delta = 1.
*/
{
  /*Check to make sure the plan structure is valid.*/
  if (plan == 0){
    FCT_ERROR(status, FCT_ENULL_PLAN);
    return;
  }

  plan->offset = offset;
  plan->dt     = delta;
}


STORAGE_CLASS
void fct_set_max_segments(fct_plan *plan, int max, fct_status* const status)
/*Set the maximum number of FFTs to be sent to
  the FFT engine. For certain FFT packages, this
  can increase the efficiency of the FFT calculations.
*/
{
  /*Check to make sure the plan structure is valid.*/
  if (plan == 0) {
    FCT_ERROR(status, FCT_ENULL_PLAN);
    return;
  }

  plan->max_number_of_segments = max;
}


STORAGE_CLASS
void fct_set_phase_function(fct_plan *plan,
			    int dimension,
			    function_pointer func,
			    fct_status* const status)
/*This function sets the phase function pointed to by func to the array of
  phase functions needed by the FCT. The phase function specified will be used
  to calculate the phase along the dimension specified by the variable
  dimension. Before adding the
  function, this routine checks for a valid plan pointer and dimension.*/
{
  /*Check to make sure the plan structure is valid.*/
  if (plan == 0) {
    FCT_ERROR(status, FCT_ENULL_PLAN);
    return;
  }

  /*
    Check to make sure the dimension is valid. Note that for some reason
    the indexing starts at 1, the 0th element is wasted.
   */
  if (dimension <= 0 || (dimension >= plan->number_of_dimensions)) {
    FCT_ERROR(status, FCT_EDIM);
    return;
  }

  plan->phase_functions[dimension] = func;
}

STORAGE_CLASS
void fct_set_oversampling_factor(fct_plan *plan,
				 int dimension,
				 int ofac,
				 fct_status* const status)
/*This function sets the oversampling factor for the phase function for the
  given dimension. Note that oversampling factors are initialised to 1 by
  default.
  Before adding the function, this routine checks for a valid plan pointer
  and dimension.*/
{
  /*Check to make sure the plan structure is valid.*/
  if (plan == 0) {
    FCT_ERROR(status, FCT_ENULL_PLAN);
    return;
  }

  /*
    Check to make sure the dimension is valid. Note that for some reason
    the indexing starts at 1, the 0th element is wasted.
   */
  if (dimension <= 0 || (dimension >= plan->number_of_dimensions)) {
    FCT_ERROR(status, FCT_EDIM);
    return;
  }

  /*
    Check to make sure the oversampling factor is valid.
   */
  if (ofac == 0) {
    FCT_ERROR(status, FCT_EOFAC);
    return;
  }

  plan->ofac[dimension] = ofac;
}


STORAGE_CLASS
void fct_add_data_cube(fct_plan *const plan,
		       const int *const start_locations,
		       const int *const end_locations,
		       const int *stride,
		       const float *const mask,
		       const int mode,
		       fct_status* const status)
{
  int       k = 0;
  data_cube *parameters = 0;

  /*Check to make sure the plan structure is valid.*/
  if (plan == 0) {
    FCT_ERROR(status, FCT_ENULL_PLAN);
    return;
  }

  /* Mask is unimplemented and should always be zero */
  if (mask != 0) {
    FCT_ERROR(status, FCT_ENONNULL_MASK);
    return;
  }

  /*allocate a new data_cube structure*/

  parameters = (data_cube *)fct_malloc(sizeof(data_cube));

  /*Check to make sure the memory was allocated.*/
  if (parameters == 0) {
    FCT_ERROR(status, FCT_EMEM);
    return;
  }

  /*Allocate the memory needed within the data_cube structure*/
  parameters->start_locations      = (int *)fct_malloc(sizeof(int)*(plan->number_of_dimensions));
  parameters->end_locations        = (int *)fct_malloc(sizeof(int)*(plan->number_of_dimensions));
  parameters->stride               = (int *)fct_malloc(sizeof(int)*(plan->number_of_dimensions));
  parameters->points_per_dimension = (int *)fct_malloc(sizeof(int)*(plan->number_of_dimensions));

  /* Check the allocations */
  if (   (parameters->start_locations == 0)
      || (parameters->end_locations == 0)
      || (parameters->stride == 0)
      || (parameters->points_per_dimension == 0))
  {
    FCT_ERROR(status, FCT_EMEM);
    return;
  }

  /*Set the next and previous pointers to NULL*/
  parameters->next     = (data_cube *)NULL;
  parameters->previous = (data_cube *)NULL;

  /*Set the stride parameters if stride is not null, otherwise set them to 1.*/
  if (stride == (int *)NULL) {
    for(k=1;k<plan->number_of_dimensions;++k){
      parameters->stride[k] = 1;
    }
  } else {
    for(k=1;k<plan->number_of_dimensions;++k){
      parameters->stride[k] = stride[k];
    }
  }

  /*Set the dimension 0 parameters to include every point along that
    dimension. Currently, this is not user definable.*/
  parameters->stride[0]          = plan->dimension_0_stride;
  parameters->start_locations[0] = 0;
  parameters->end_locations[0]   = plan->data_length;

  /*Now use the mode parameter to determine how to initialize the other parameters*/
  switch(mode) {

  case FCT_CALCULATE_ALL:
    /*The start and end locations are set to
      their minimal and maximal values.*/
    for (k=1;k<plan->number_of_dimensions; ++k){
      parameters->start_locations[k] = 0;
      parameters->end_locations[k]   = plan->data_length;
    }

    break;

  case FCT_SPECIFY_RANGES:
    /*The start and end location are user specified.*/
    /*First, check to make sure the user put in valid arrays*/
    if (start_locations == (int *)NULL){
      FCT_ERROR(status, FCT_ENULL_LOC);
      return;
    }

    if (end_locations == (int *)NULL){
      FCT_ERROR(status, FCT_ENULL_LOC);
      return;
    }

    for (k = 1; k<plan->number_of_dimensions; ++k) {
      if (end_locations[k] <= start_locations[k])
      {
	FCT_ERROR(status, FCT_EDATACUBE_RANGE);
	return;
      }
      parameters->start_locations[k] = start_locations[k];
      parameters->end_locations[k]   = end_locations[k];
    }

    break;

  default:
    /*If neither of the above modes were specified, then give an error */
    FCT_ERROR(status, FCT_EDATACUBE_MODE);
    return;
  }


/*Calculate the number of points per dimension and the total number of bytes
  needed for the FCT of this data cube.*/

  parameters->total_length = 1;

  for (k=0; k < plan->number_of_dimensions; ++k){

    parameters->points_per_dimension[k]
      = ceil((float)(parameters->end_locations[k]
		   - parameters->start_locations[k])/(parameters->stride[k]));
    parameters->total_length *= parameters->points_per_dimension[k];
  }

  parameters->total_length *= 2;
  /*Calculate the kernal_length and the pre_kernal length
  parameters->kernal_length     = (parameters->total_length * (plan->data_length
                                     /parameters->points_per_dimension[0]));

  */
  /*Calculate the pre_kernal_length*/
  parameters->pre_kernal_length = 0;
  for(k=1; k < plan->number_of_dimensions; ++k){
    parameters->pre_kernal_length += parameters->points_per_dimension[k];
  }

  parameters->pre_kernal_length *= 2*plan->data_length;

/*Now, add this data_cube to the parameter_list.*/

  if (plan->parameter_list_start == NULL) {
    plan->parameter_list_start = parameters;
    plan->parameter_list_end   = parameters;
  } else {
    plan->parameter_list_end->next = parameters;
    parameters->previous           = plan->parameter_list_end;
    plan->parameter_list_end       = parameters;
  }
}


STORAGE_CLASS
void fct_remove_data_cube(fct_plan *plan, fct_status* const status)
/*Removes the most recently added member to the
  data cube list*/
{
  data_cube* end_predecessor = 0;

  /*Check to make sure the plan structure is valid.*/
  if (plan == 0) {
    FCT_ERROR(status, FCT_ENULL_PLAN);
    return;
  }

  if (plan->parameter_list_start == NULL)
  {
    /* If it's an empty list, just return */
    return;
  } else {
    /* Get the predecessor of the last node */
    end_predecessor = plan->parameter_list_end->previous;

    /* Destroy the old node */
    fct_destroy_data_cube(plan->parameter_list_end, status);
    if (status->fct_errno != 0)
    {
      return;
    }

    /* Make the "new" last node be the predecessor of the original one */
    plan->parameter_list_end = end_predecessor;

    /*
      Set the successor of the last node to zero.

      :NOTE: A special case is when there was only one node on the list.
      Then the "new" last node will be null and dereferencing it will be
      bad. We need to check if the last node is null before referring to
      it's successor. Also, if the last node is null we must also set
      the first node to null.
    */
    if (plan->parameter_list_end == 0)
    {
      plan->parameter_list_start = 0;
    }
    else
    {
      plan->parameter_list_end->next = 0;
    }
  }
}


static
unsigned long fct_data_cube_size(fct_plan *plan, data_cube *cube,
				 fct_status* const status)
/*Calculates the number of real output points generated by
  the specified data cube.*/
{
  unsigned long total = 0;
  int           k = 0;

  /*Check to make sure the plan structure is valid.*/
  if (plan == 0) {
    FCT_ERROR(status, FCT_ENULL_PLAN);
    return 0;
  }

  /*Check to make sure the data_cube structure is valid. If not, return 0.*/
  if (cube == (data_cube *)NULL) {
    return (0);
  }

  total = 1;
  for (k=0;k<plan->number_of_dimensions;++k) {
    total *= ceil((float)(cube->end_locations[k] - cube->start_locations[k])/
		  (cube->stride[k]));
  }

  /*Since total equals the total number of points, 2*total is the
    total number of real data points since each point is complex.*/

  return (2*total);
}

STORAGE_CLASS
unsigned long fct_output_data_size(fct_plan *plan, fct_status* const status)
/*Using the current plan struct, this function calculates how large the
  output data array will be in real data points. Note, the FCT output is
  a series of complex numbers. This function takes into account the factor
  of two conversion to obtain the number of real numbers.*/
{
  unsigned long total = 0;
  data_cube     *current = 0;


  /*Check to make sure the plan structure is valid.*/
  if (plan == 0) {
    FCT_ERROR(status, FCT_ENULL_PLAN);
    return 0;
  }

  /*Check to make sure that a least one data_cube has been added to the list*/
  if (plan->parameter_list_start == (data_cube *)NULL) {
    FCT_ERROR(status, FCT_ENULL_DATACUBE);
    return 0;
  }

  /*Calculate the number of output points that will be
    generated by the first data_cube.*/
  /*total = fct_data_cube_size(plan,plan->parameter_list_start);*/
  total = plan->parameter_list_start->total_length;

  /*Now, scan through the rest of the list and calculate the
    total number of data points.*/
  current = plan->parameter_list_start;

  while (current->next != NULL) {
    total  += current->next->total_length;
    current = current->next;
  }

  return total;
}


STORAGE_CLASS
void fct_calculate(fct_plan *plan,
		   const fct_real* const input_data,
		   fct_real* const output_data,
		   fct_status* const status)
/*Initializes the FCT using the information in the plan file.*/
{
  /* kp is used to count # parameters */

  int           pre_ker_loc,output_real_loc,output_imag_loc,j,k,kp,dim,seg=0;
  int           input_imag_loc;
  data_cube     *current;
  long          *k_loc;
  unsigned long *index_strides;
  unsigned long *pre_k_index_offsets;
  unsigned long output_location;
  double        partial_phase;
  long          rint_phase = 0;
  fct_real      *pre_kernal;
  fct_real      **pk_real;
  fct_real      **pk_imag;
  fct_real      *work_space,*current_workspace;
  double        wpr,wpi,wtemp,wtempi;  /* wpr, wpi, wtemp(i) used  the trigonometri recurrence */
  fct_real      hold;

  /*Check to make sure the plan structure is valid.*/
  if (plan == 0) {
    FCT_ERROR(status, FCT_ENULL_PLAN);
    return;
  }

  /*Allocate temporary data arrays*/
  k_loc                = (long  *)fct_malloc(sizeof(long)*plan->number_of_dimensions);
  index_strides        = (unsigned long  *)fct_malloc(sizeof(unsigned long)*plan->number_of_dimensions);
  pre_k_index_offsets  = (unsigned long  *)fct_malloc(sizeof(unsigned long)*plan->number_of_dimensions);
  pk_real              = (fct_real **)fct_malloc(sizeof(fct_real *)*plan->number_of_dimensions);
  pk_imag              = (fct_real **)fct_malloc(sizeof(fct_real *)*plan->number_of_dimensions);

  /*Make sure memory has been allocated*/
  if (k_loc              == NULL ||
     index_strides       == NULL ||
     pre_k_index_offsets == NULL ||
     pk_real             == NULL ||
     pk_imag             == NULL)
  {
    FCT_ERROR(status, FCT_EMEM);
    return;
  }

  work_space = (fct_real *)fct_malloc(sizeof(fct_real)*
	       (int)ceil((float)2* plan->max_number_of_segments *
	       plan->data_length / plan->dimension_0_stride));

  /*Check to make sure the workspace is allocated*/
  if (work_space == NULL) {
    FCT_ERROR(status, FCT_EMEM);
    return;
  }

    /*Calculate the kernal product (K=exp(k1*phi_1)*exp(k2*phi_2)*...*), dot product (h*K), and FFTW (h*k) at the same loop  for each data cube.*/

  current        = plan->parameter_list_start;
  output_location = 0;
  /*Loop over each entry in the parameter_list.*/


  while (current != (data_cube *)NULL)
  {
    /*Zero the workspace.*/
    for (j=0;j < 2*current->points_per_dimension[0]*
	   plan->max_number_of_segments; ++j) {
      work_space[j] = 0;
    }

    /*Allocate the memory for this data cube's pre_kernal*/
    pre_kernal = fct_malloc(sizeof(fct_real)*current->pre_kernal_length);

    if (pre_kernal == NULL) {
      FCT_ERROR(status, FCT_EMEM);
      return;
    }


    /*Calculate the index_strides and pre_k_index_offsets arrays*/
    index_strides[0]       = current->points_per_dimension[0];
    pre_k_index_offsets[0] = 0;
    for (dim = 1; dim<plan->number_of_dimensions-1; ++dim) {
      index_strides[dim] = index_strides[dim-1] * current->points_per_dimension[dim];
      pre_k_index_offsets[dim] = pre_k_index_offsets[dim-1] +
	plan->data_length * current->points_per_dimension[dim];
    }

      /*Calculate the pre_kernal array*/
      for (j=0; j < plan->data_length; ++j) {

	for (dim=0;dim < plan->number_of_dimensions-1;++dim) {
	  /* remove rint here */
	  /* original one is : */

	  /*
	    Can't use rint() because it's non-ANSI. For positive
	    numbers, adding 0.5 and copying to an int is equivalent.
	    Doesn't work properly for negative numbers though.
	  */

	  rint_phase = plan->data_length*
	    (*plan->phase_functions[dim+1])(plan->offset + plan->dt*j) + 0.5;

	  partial_phase = 2.0*M_PI*rint_phase/
	    (plan->ofac[dim+1]*plan->data_length);

	  /* Original function with call to rint()

	  partial_phase = (2.0*M_PI/plan->data_length)*
	    rint(plan->data_length *
		 (*plan->phase_functions[dim+1])(plan->offset + plan->dt*j));

	  */

	  /* it creates problems for phase_functions such as 1/f, the k value and number of data points are restricted if we use the approximations */

/*
          partial_phase   = (2.0*M_PI/plan->data_length)*
	     (plan->data_length * (*plan->phase_functions[dim+1])
	        (plan->offset + plan->dt*j));
*/
	  /* debug */
	  /*	  printf(" data_length= %i data_length *(phase_functions[dim+1])(plan->offset + dt*j=%f\n", plan->data_length, rint(plan->data_length *(*plan->phase_functions[dim+1])(plan->offset + plan->dt*j))); */

	  /*  printf("partial phase at j=%i is %g\n", j, partial_phase); */

	  /* partial_phase here is equivalent to the incrementing theta */

	  pk_real[0] = &pre_kernal[2*(j+pre_k_index_offsets[dim])];
	  pk_imag[0] = pk_real[0] + 1;

	  /* new change for a trigonometric recurrence,
	     wpr(i)=incrementing theta */
	  wtemp=sin(partial_phase/2.*current->stride[dim+1]);
	  wpr=-2*wtemp*wtemp;
	  wpi=sin(partial_phase*current->stride[dim+1]);

	  /* for k=0 , wtemp=theta_old*/

	  *pk_real[0] = cos(current->start_locations[dim+1]*partial_phase);
	  *pk_imag[0] = sin(current->start_locations[dim+1]*partial_phase);

	  wtemp  = *pk_real[0];
	  wtempi = *pk_imag[0];

	  pk_real[0] += 2*plan->data_length;
	  pk_imag[0] += 2*plan->data_length;

	  for (k = 1;k < current->points_per_dimension[dim+1]; ++k) {

	    *pk_real[0]    = wtemp*wpr-wtempi*wpi+wtemp;
	    *pk_imag[0]    = wtempi*wpr+wtemp*wpi+wtempi;

	    wtemp = *pk_real[0];
	    wtempi = *pk_imag[0];

	    pk_real[0] += 2*plan->data_length;
	    pk_imag[0] += 2*plan->data_length;

	  }
	}

      }

      /*for(j=0;j<8;++j){
	printf("(%lf,%lf)\n",pre_kernal*/

      /*Now going through the entire parameter space */

      /*Initialize the k_loc array*/
      for (dim = 0; dim < plan->number_of_dimensions; ++dim) {
	k_loc[dim]=0;
      }

      /* merge fct_calculate in */
      kp = 0;
      while (k_loc[plan->number_of_dimensions-1] == 0) {

	seg = kp % plan->max_number_of_segments;
	current_workspace = &(work_space)[2*seg*current-> points_per_dimension[0]];

	kp++; /* counting # of parameters going through */

	/*	  kernal_loc = 0; */
	for (dim = 0; dim < plan->number_of_dimensions - 1; ++dim) {
	  /*	    kernal_loc      += index_strides[dim]*k_loc[dim]; */
	  pre_ker_loc      = (pre_k_index_offsets[dim] + k_loc[dim]*plan->data_length);
	  pk_real[dim]     = &pre_kernal[2*pre_ker_loc];
	  pk_imag[dim]     = &pre_kernal[2*pre_ker_loc+1];
	}

	for (j=0; j < 2*plan->data_length; j += 2) {
	  output_real_loc = 2*((j/2) % current->points_per_dimension[0]);
	  output_imag_loc = output_real_loc + 1;

	  input_imag_loc = j + 1;

	  wtemp  = pk_real[0][j];
	  wtempi = pk_imag[0][j];

	  /* debug */
	  /*	    printf("pk[0][j]= %f %f\n",pk_real[0][j],pk_imag[0][j]); */

	  for (dim=1; dim < plan->number_of_dimensions-1; ++dim) {
	    hold  = wtemp;
	    wtemp = wtemp*pk_real[dim][j] - wtempi*pk_imag[dim][j];

	    wtempi = hold*pk_imag[dim][j] + wtempi*pk_real[dim][j];
	  }

	  current_workspace[output_real_loc]
	    += wtemp*input_data[j] - wtempi*input_data[input_imag_loc];

	  current_workspace[output_imag_loc]
	    += wtemp*input_data[input_imag_loc] + wtempi*input_data[j];

	  /*  debug */

	  /*	    	    		printf("output_real_loc= %i current_workspace[output_real_loc]= %f current_workspace[output_imag_loc]= %g\n" , output_real_loc, current_workspace[output_real_loc] , current_workspace[output_imag_loc]);
	     */
	  /*	    printf("%f +i%f\n",current_workspace[output_real_loc],current_workspace[output_imag_loc]);  */


	}
	if (seg == plan->max_number_of_segments-1) {
	  /*FFT the data*/
	  fct_fft(plan, work_space,
	    &output_data[2*output_location*(current->points_per_dimension[0])],
            plan->max_number_of_segments);

	  /*Zero the workspace again.*/
	  for (j = 0; j < 2*current->points_per_dimension[0]*plan->max_number_of_segments; ++j) {
	    work_space[j] = 0;
	  }
	  output_location += plan->max_number_of_segments;
	}

	/*Increment the k_loc array*/
	++k_loc[0];

	for (dim = 0; dim < plan->number_of_dimensions-1; ++dim) {
	  if (k_loc[dim] >= current->points_per_dimension[dim+1]) {
	    k_loc[dim] = 0;
	    ++k_loc[dim+1];
	  }
	}
      }

      /*FFT the remaining data segments*/
      if (seg != plan->max_number_of_segments - 1) {
	fct_fft(plan,work_space,
     &output_data[2*output_location*(current->points_per_dimension[0])],seg+1);
	output_location += seg + 1;
      }

      /*Free the previous pre_kernal array*/
      fct_free(pre_kernal);
      current = current->next;

    }

  /*Delete temporary data arrays.*/

  fct_free(k_loc);
  fct_free(pre_k_index_offsets);
  fct_free(index_strides);
  /*  fct_free(pre_kernal);*/
  fct_free(pk_real);
  fct_free(pk_imag);
  fct_free(work_space);
}




STORAGE_CLASS
void fct_destroy_plan(fct_plan *plan, fct_status* const status)
/*Deallocates the memory within a fct_plan structure including*/
/*the structure itself.*/
{
  data_cube *current = 0, *next = 0;

  /*Check to make sure the plan structure is valid.*/
  if (plan == 0) {
    FCT_ERROR(status, FCT_ENULL_PLAN);
    return;
  }

  fct_free(plan->phase_functions);
  fct_free(plan->ofac);
  fct_destroy_fft_plan(plan);

  current = plan->parameter_list_start;

  while(current != NULL){
    next = current->next;
    fct_destroy_data_cube(current, status);
    if (status->fct_errno != 0)
    {
      return;
    }
    current = next;
  }
  fct_free(plan);
}


STORAGE_CLASS
void fct_destroy_data_cube(data_cube *cube, fct_status* const status)
/*Deallocates the memory within a data_cube structure including*/
/*the structure itself.*/
{
  /*First, make sure the pointer is not NULL.*/
  if (cube == (data_cube *)NULL) {
    FCT_ERROR(status, FCT_ENULL_DATACUBE);
    return;
  }

  fct_free(cube->start_locations);
  fct_free(cube->end_locations);
  fct_free(cube->stride);
  fct_free(cube->points_per_dimension);
  /*fct_free(cube->mask);*/
  fct_free(cube);

}

STORAGE_CLASS
int writefct(const fct_real* const out, const int n, const int m,
	     const char* const filename)
{
    /* The size of each element in the out array */
    const float out_size = sizeof(*out);

    const fct_real M = m;
    const fct_real N = n;

    /* Open the file */
    const int fd_output = creat(filename, 0644);

    if (fd_output != 0)
    {
	/* Write info about endian-ness and local sizeof(*out) */
	write(fd_output, &out_size, sizeof(out_size));

	/* Write the dimensions of the fct */
	write(fd_output, &M, sizeof(M));
	write(fd_output, &N, sizeof(N));

	/* Write the body of the fct */
	write(fd_output, out, sizeof(*out)*2*n*m);

	close(fd_output);
    }
    else
    {
	return -1;
    }

    return 0;
}

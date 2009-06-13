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

/*This is the header file for the official implentation of the
  Fast Chirp Transform (FCT). This implementation was written by
  Fredrick A. Jenet at the Calfornia Institute of Technology.*/
/*Test Compile line: cc fct.c testcode.c fct_fft.c -lm libfftw.a -o testcode*/

#ifndef FCT_H
#define FCT_H

#include <lal/LALRCSID.h>
NRCSID (FCTH,"$Id$");

#define FCT_CALCULATE_ALL  0
#define FCT_SPECIFY_RANGES 1

#define FCT_INIT_DEPRECATED

/* For compatability with the LIGO Analysis Library (LAL) */
#ifdef LALFCT
#define STORAGE_CLASS static
#else
#define STORAGE_CLASS
#endif

typedef float    fct_real;
typedef fct_real (*function_pointer)(fct_real);

/*
  Typedefs for replacing malloc- and free-like functions
*/
typedef void* (*fct_malloc_type_function)(size_t size);
typedef void* (*fct_calloc_type_function)(size_t nmemb, size_t size);
typedef void  (*fct_free_type_function)(void* p);

/*
  Forward declarations
*/
struct fct_fft_plan_;
typedef struct fct_fft_plan_ fct_fft_plan;

/*
  Set the following external variables to the address of a user-defined
  replacement for malloc, calloc or free
  eg. to use normal system malloc()

      fct_malloc_hook = malloc;
*/
extern fct_malloc_type_function fct_malloc_hook;
extern fct_calloc_type_function fct_calloc_hook;
extern fct_free_type_function   fct_free_hook;

typedef struct data_cube_ {
/* This structure defines the region
   over which the FCT will be calculated. */

  int               *start_locations;     /*An array with a
					    dimensionality equal to
					    fct_plan.number_of_dimensions -1
					    which sets the lower bounds of
					    the parameter ranges. */

  int               *end_locations;       /*An array with a
					    dimensionality equal to
					    fct_plan.number_of_dimensions -1
					    which sets the upper bounds of
					    the parameter ranges. */

  int               *stride;              /*An array with a
					    dimensionality equal to
					    fct_plan.number_of_dimensions -1
					    which sets the stride of the
					    parameter ranges. To evalute all points
					    between start_locations[1] and end_locations[1],
					    set stride[1] = 1. To evaluate every other point, set
					    stride[1] = 2.*/

  int               *points_per_dimension;/*The number of data points per dimension taking
					    the stride parameter into account.*/

  unsigned long     total_length;         /*The total length in bytes of the output when the
					    FCT is calculated in the region defined
					    by this data_cube.*/

  /*   unsigned long     kernal_length;        The length of the kernal in bytes for this
					    data_cube.*/

  unsigned long     pre_kernal_length;    /*The length of the pre_kernal in
					    bytes for this data_cube*/

  /*   fct_real          *kernal;              A parameter dependent
					    set of internally calculated
					    complex numbers. It is calculated
					    using the phase functions
					    supplied by the user.*/

  fct_real          *mask;                /*A parameter dependent mask
					    to be applied to the data.
					    (Not implemented yet)*/

  struct data_cube_ *previous;            /*The previous data_cube structure in
					    the linked list.*/

  struct data_cube_ *next;                /*The next data_cube structure in
					    the linked list.*/

} data_cube;


typedef struct fct_plan_
{
/*This structure holds the necessary user supplied information
  needed in order to calculate the FCT*/


  int                data_length;              /*Number of complex samples in the
						 input data set. Currently, this must
						 be a power of 2.*/


  int                number_of_dimensions;     /*Total number of FCT dimensions or
						 conjugate parameters.*/


  int                dimension_0_stride;       /*Set the stride of the dimension 0
						 parameter range. Currently, this number
						 must be a power of 2. To calculate all
						 points, it should be set to 1.*/

  float              offset;                   /*Used in converting from sample number
						 to time. If i is the sample number,
						 the phase functions will be evaluated
						 at phase_function(offset + i*dt). The
						 default value is 0. Use fct_set_units() to
						 change this value.*/

  float              dt;                       /*Used in converting from sample number to
						 time. If i is the sample number,
						 the phase functions will be evaluated
						 at phase_function(offset + i*dt). The
						 default value is 1. Use fct_set_units() to
						 change this value.*/

  function_pointer  *phase_functions;          /*The array of phase functions
						 used to calculate the FCT.*/

  int               *ofac;                     /* Oversampling factors for
                                                  each phase function - set
                                                  all to 1 for "standard" FCT
					       */

  data_cube         *parameter_list_start;     /*The start of a linked list of
						 structures defining the region
						 in FCT parameter space to be
						 calculated by this FCT.
						 */

  data_cube         *parameter_list_end;       /*A pointer to the last element of
						 the parameter_list.*/

  fct_fft_plan      *fft_plan;                 /*This structure contains the
						 setup information for the
						 fft package being used
						 with the FCT.*/


  fct_real           *work_space;

  int                max_number_of_segments;   /*This is the maximum number
						 of FFT that will be sent
						 to the FFT engine at any given
						 time. Setting this number
						 greater than one can greatly
						 increase performance.*/

  int                error;                    /*If any errors occur while
						 setting up the FCT or while
						 processing the data, error
						 will be set to a non-zero
						 number.*/

} fct_plan;


/* An enumeration type for holding error codes */
enum FCTErrNo {
  FCT_ENONE,
  FCT_ENULL_PLAN,
  FCT_ENOT_PWR2,
  FCT_EMEM,
  FCT_EFFTW_PLAN,
  FCT_EDIM,
  FCT_ENULL_LOC,
  FCT_ENONNULL_MASK,
  FCT_ENULL_DATACUBE,
  FCT_EDATACUBE_MODE,
  FCT_EDATACUBE_RANGE,
  FCT_EOFAC,
  FCT_EUNKNOWN
};


/* A structure to hold information about internal FCT errors */
typedef struct fct_status_
{
  int fct_errno;  /* Non-zero value indicates an error */
  const char* file;
  int line;
} fct_status;

#define FCT_ERROR(s, errno) \
{ \
  s->fct_errno = errno; \
  s->file = __FILE__; \
  s->line = __LINE__; \
}


/* Function prototypes */

STORAGE_CLASS
const char* fct_strerror(int fct_errno);
/*
  Passing an error code from the FCTErrNo enumerated type returns a
  pointer to a string holding the corresponding error message
*/

/* Replacements for the standard malloc, calloc and free */
STORAGE_CLASS
void* fct_malloc(size_t size);

STORAGE_CLASS
void* fct_calloc(size_t nmemb, size_t size);

STORAGE_CLASS
void  fct_free(void* p);

STORAGE_CLASS
fct_plan            *fct_init_plan               (int data_length,
						  int number_of_dimensions,
						  int dimension_0_stride,
						  fct_status* const status);
/*This function partially initializes the FCT. It must be called before
  any other fct_ function. It generates an fct_plan structure which tells
  the other fct_ functions how to operate on the data. After calling
  fct_init_plan, the user must specify both the phase functions using
  fct_set_phase_function and the FCT parameters over which to calculate
  the FCT using fct_add_indicies.

  The returned value is a pointer to an fct_plan_ structure which is
  allocated on the heap. The pointer is freed by calling fct_destroy_plan */

STORAGE_CLASS
void                fct_set_units                (fct_plan *plan,
						  float offset,
						  float delta,
						  fct_status* const status);
/*Sets the conversion from data sample to the argument of the
  phase function. For the jth data point, the phase functions are
  evaluated at offset + j*delta. The default values are offset = 0 and
  delta = 1.
*/

STORAGE_CLASS
void                fct_set_max_segments         (fct_plan *plan, int max,
						  fct_status* const status);
/*Set the maximum number of FFTs to be sent to
  the FFT engine. For certain FFT packages, this
  can increase the efficiency of the FFT calculations.
*/


STORAGE_CLASS
void                fct_set_phase_function       (fct_plan *plan,
						  int dimension,
						  function_pointer func,
						  fct_status* const status);
/*This function adds the phase function pointed to by func to the array of
  phase functions needed by the FCT. The phase function specified will be used
  to calculate the phase along the dimension speficied by the variable
  dimension. Before adding the function, this routine checks for a valid
  plan pointer and dimension.*/

STORAGE_CLASS
void fct_set_oversampling_factor(fct_plan *plan,
				 int dimension,
				 int ofac,
				 fct_status* const status);
/*This function sets the oversampling factor for the phase function for the
  given dimension. Note that oversampling factors are initialised to 1 by
  default.
  Before adding the function, this routine checks for a valid plan pointer
  and dimension.*/


STORAGE_CLASS
void                fct_add_data_cube            (fct_plan *const plan,
						  const int * const start_locations,
						  const int * const end_locations,
						  const int *stride,
						  const float *const mask,
						  const int mode,
						  fct_status* const status);

STORAGE_CLASS
void                fct_remove_data_cube         (fct_plan *plan,
						  fct_status* const status);
/*Removes the most recently added member to the
  data cube list*/

STORAGE_CLASS
unsigned long       fct_output_data_size         (fct_plan *plan,
						  fct_status* const status);

STORAGE_CLASS
void                fct_calculate                (fct_plan *plan,
						  const fct_real* const input_data,
						  fct_real* const output_data,
						  fct_status* const status);

/*Initializes the FCT using the information in the plan file.

void                fct_calculate                (fct_plan *plan,
						  const fct_real *input_data,
						  fct_real *output_data,
						  fct_status* const status);
*/

/*This function calculates the FCT of the input_data array and places it into*/
/*the output_data array. Proper initialization must take place before this function*/
/*can be called.*/

STORAGE_CLASS
void                fct_destroy_plan             (fct_plan *plan,
						  fct_status* const status);
/*Deallocates the memory within a fct_plan structure includingg*/
/*the structure itself.*/

STORAGE_CLASS
void                fct_destroy_data_cube        (data_cube *cube,
						  fct_status* const status);
/*Deallocates the memory within a data_cube structure including*/
/*the structure itself.*/

STORAGE_CLASS
int writefct(const fct_real* const out, const int n, const int m,
	     const char* const filename);
/*
  Write a 2-d FCT to the given file.
  If successful, 0 is returned, otherwise -1
*/

#endif





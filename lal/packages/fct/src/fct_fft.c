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
  This file contains the fft routines used in the official FCT implementation.
*/

#include "fct_fft.h"

#include <lal/LALRCSID.h>
NRCSID (FCT_FFTC,"$Id$");

STORAGE_CLASS
void fct_setup_fft(fct_plan *plan, fct_status* const status)
{
  fct_fft_plan *fft_plan;
  int          fft_length;

  /*Check to make sure the plan structure is valid.*/
  if (plan == (fct_plan *)NULL) {
    status->fct_errno = FCT_ENULL_PLAN;
    return;
  }

  /*Calcuate the fft length needed and see if it is compatible with
    the FFT package.*/
  fft_length = plan->data_length/plan->dimension_0_stride;

  if (!powerof2(fft_length)) {
    status->fct_errno = FCT_ENOT_PWR2;
    return;
  }

  /*Allocate the memory for the fct_fft_plan structure.*/
  fft_plan = (fct_fft_plan *)fct_malloc(sizeof(fct_fft_plan));

  if (fft_plan == NULL) {
      status->fct_errno = FCT_EMEM;
      return;
  }

  /*Create the plan for fftw*/
  fft_plan->plan = fftw_create_plan(fft_length, FFTW_BACKWARD, 0);
  if (fft_plan->plan == 0)
  {
      status->fct_errno = FCT_EFFTW_PLAN;
      return;
  }

  /*Set the fft_length variable in the fft_plan structure*/
  fft_plan->fft_length = fft_length;

  /*Copy the fct_fft_plan pointer in the fct_plan structure*/
  plan->fft_plan = fft_plan;
}


STORAGE_CLASS
void fct_fft(fct_plan *plan, fct_real *input, fct_real *output, int number)
{
  fftw(plan->fft_plan->plan, number,
       (fftw_complex *)input, 1, plan->fft_plan->fft_length,
       (fftw_complex *)output, 1, plan->fft_plan->fft_length);

  /*fftw(fftw_plan plan, int howmany,
          fftw_complex *in, int istride, int idist,
          fftw_complex *out, int ostride, int odist);*/
}

STORAGE_CLASS
void fct_destroy_fft_plan(fct_plan *plan)
{
  /*Destroy FFTW plan*/
  fftw_destroy_plan(plan->fft_plan->plan);

  /*Free the fct_fft_plan memory*/
  fct_free(plan->fft_plan);
}

STORAGE_CLASS
int powerof2(int data)
/*Check to see if data is a power of 2.
1->yes, 0->no.
*/
{
  int sum = 0;

  while (data) {
    sum += data & 1;
    data = data >> 1;
  }

  if (sum <= 1) {
    return (1);
  } else {
    return (0);
  }

}


/* -*- mode: c; c-basic-offset: 2; -*- */

/* The header file for the fct_fft routines. */

#ifndef FCT_FFT_H
#define FCT_FFT_H

#ifdef HAVE_SFFTW_H
#include <sfftw.h>
#elif HAVE_FFTW_H
#include <fftw.h>
#else
#error "don't have either sfftw.h or fftw.h"
#endif

#include "fct.h"

struct fct_fft_plan_
{
  int                fft_length;
  fftw_plan          plan;
} /* fct_fft_plan */;

/* Prototype FFT functions */

STORAGE_CLASS
void fct_setup_fft(fct_plan *plan, fct_status* const status);

STORAGE_CLASS
void fct_fft(fct_plan *plan,
	     fct_real *input, fct_real *output, int number);

STORAGE_CLASS
void fct_destroy_fft_plan(fct_plan *plan);

STORAGE_CLASS
int powerof2(int data);

#endif









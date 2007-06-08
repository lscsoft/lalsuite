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

#include <lal/LALRCSID.h>
NRCSID (FCT_FFTH,"$Id$");

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









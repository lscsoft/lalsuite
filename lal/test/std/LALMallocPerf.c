/*
 *  Copyright (C) 2016 Karl Wette
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

/**
 * \file
 * \ingroup LALMalloc_h
 * \brief Tests the performance of the routines in \ref LALMalloc_h.
 */

/** \cond DONT_DOXYGEN */

#include <stdio.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <lal/LALStdlib.h>
#include <lal/LALMalloc.h>
#include <lal/LogPrintf.h>

int main(void) {

  setvbuf(stdout, NULL, _IONBF, 0);

  printf("LALMallocPerf: LAL memory tracking is %s\n", (lalDebugLevel & LALMEMTRKBIT) ? "ON" : "OFF");

  const int n = 1 << 17;
  printf("LALMallocPerf: Testing performance with %i allocates\n", n);

  {
    void *x[n];
    for (int i = 0; i < n; ++i) {
      x[i] = XLALMalloc(sizeof(int));
    }
    printf("LALMallocPerf: Deallocate in same order as allocate:\t");
    const REAL8 t0 = XLALGetCPUTime();
    for (int i = 0; i < n; ++i) {
      XLALFree(x[i]);
    }
    const REAL8 t = XLALGetCPUTime() - t0;
    printf("%g sec (%e sec/deallocate)\n", t, t/n);
  }

  {
    void *x[n];
    for (int i = 0; i < n; ++i) {
      x[i] = XLALMalloc(sizeof(int));
    }
    printf("LALMallocPerf: Deallocate in reverse order as allocate:\t");
    const REAL8 t0 = XLALGetCPUTime();
    for (int i = n - 1; i >= 0; --i) {
      XLALFree(x[i]);
    }
    const REAL8 t = XLALGetCPUTime() - t0;
    printf("%g sec (%e sec/deallocate)\n", t, t/n);
  }

  {
    void *x[n];
    int ii[n];
    for (int i = 0; i < n; ++i) {
      x[i] = XLALMalloc(sizeof(int));
      ii[i] = i;
    }
    gsl_rng_env_setup();
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
    gsl_ran_shuffle(r, ii, n, sizeof(int));
    gsl_rng_free(r);
    printf("LALMallocPerf: Deallocate in random order:\t\t");
    const REAL8 t0 = XLALGetCPUTime();
    for (int i = 0; i < n; ++i) {
      XLALFree(x[ii[i]]);
    }
    const REAL8 t = XLALGetCPUTime() - t0;
    printf("%g sec (%e sec/deallocate)\n", t, t/n);
  }

  LALCheckMemoryLeaks();

  return EXIT_SUCCESS;

}

/** \endcond */

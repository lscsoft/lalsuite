/*
 * Copyright (C) 2017  Leo Singer
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA  02111-1307  USA
 */


#ifndef OMP_INTERRUPTIBLE_H
#define OMP_INTERRUPTIBLE_H

#include <signal.h>
#include <stdlib.h>


static __thread int *omp_interruptible_flag_ptr = NULL;
static __thread void (*omp_interruptible_old_handler)(int) = NULL;


static void omp_interruptible_restore_handler(int sig)
{
    signal(sig, omp_interruptible_old_handler);
    omp_interruptible_old_handler = NULL;
    omp_interruptible_flag_ptr = NULL;
}


static void omp_interruptible_handler(int sig)
{
    *omp_interruptible_flag_ptr = 1;
    #pragma omp flush
    omp_interruptible_restore_handler(sig);
    raise(sig);
}


static void omp_interruptible_set_handler(int sig, int *flag_ptr)
{
    omp_interruptible_flag_ptr = flag_ptr;
    *omp_interruptible_flag_ptr = 0;
    omp_interruptible_old_handler = signal(sig, omp_interruptible_handler);
}


#define OMP_BEGIN_INTERRUPTIBLE { \
    int omp_was_interrupted; \
    omp_interruptible_set_handler(SIGINT, &omp_was_interrupted);


#define OMP_END_INTERRUPTIBLE \
    omp_interruptible_restore_handler(SIGINT); \
}


#define OMP_WAS_INTERRUPTED omp_was_interrupted


#if _OPENMP
#define OMP_EXIT_LOOP_EARLY continue;
#else
#define OMP_EXIT_LOOP_EARLY break;
#endif


#endif /* OMP_INTERRUPTIBLE_H */

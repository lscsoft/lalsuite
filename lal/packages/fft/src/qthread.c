/*
*  Copyright (C) 2007 Bernd Machenschalk, Duncan Brown
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

/*-----------------------------------------------------------------------
 *
 * File Name: qthread.c
 *
 * Author: Brown, D. A.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <pthread.h>
#include <bits/local_lim.h>

#include <lal/LALRCSID.h>
NRCSID (QTHREADC,"$Id$");

/* lal debug flag for verbosity level */
extern int lalDebugLevel;

/* dummy symbol to ensure that we link this to the fft routines */
int dummy_have_qthread;

/* the following functions are needed to resolve symbols in     */
/* executables linked to the Intel Math Kernel Library DFT      */
/* routines. The function have an implementation that should    */
/* work as long as they are not actually run in threadded       */
/* code. These routines are called by the DFT routines          */
/* during execution, even if the MKL environment variable       */
/* MKL_SERIAL=YES                                               */

/* once-only execution */
enum { NEVER = 0, IN_PROGRESS = 1, DONE = 2 };

int pthread_once(pthread_once_t * once_control, void (*init_routine)(void))
{
  /* call the function init_routine() if once_control is not set to done */
  if ( lalDebugLevel & 4 )
    fprintf( stdout, "calling pthread_once(%p,%p)\n",
        once_control, init_routine );

  if ( *once_control == DONE ) return 0;
  init_routine();
  *once_control = DONE;
  return 0;
}

/* thread specific keys */
typedef void (*destr_function)(void *);

struct pthread_key_struct
{
  int in_use;                   /* already allocated? */
  destr_function destr;         /* destruction routine */
};

/* table of keys */
static struct pthread_key_struct pthread_keys[PTHREAD_KEYS_MAX] =
{ { 0, NULL } };

/* value in keys */
static void * pthread_key_values[PTHREAD_KEYS_MAX];

/* create a new key */
int pthread_key_create (pthread_key_t * key, destr_function destr)
{
  int i;

  if ( lalDebugLevel & 4 )
    fprintf( stdout, "calling pthread_key_create(%p,%p)\n", key, destr );

  for ( i = 0; i < PTHREAD_KEYS_MAX; i++ )
  {
    if (! pthread_keys[i].in_use)
    {
      /* Mark key in use */
      pthread_keys[i].in_use = 1;
      pthread_keys[i].destr = destr;
      *key = i;

      /* set the value of the key to null */
      pthread_key_values[i] = NULL;

      return 0;
    }
  }
  return EAGAIN;
}

/* delete a key */
int pthread_key_delete(pthread_key_t key)
{
  if ( lalDebugLevel & 4 )
    fprintf( stdout, "calling pthread_key_delete(%d)\n", key );

  if ( key >= PTHREAD_KEYS_MAX || ! pthread_keys[key].in_use )
  {
    return EINVAL;
  }

  pthread_keys[key].in_use = 0;
  pthread_keys[key].destr = NULL;

  /* set the value of the key to null */
  pthread_key_values[key] = NULL;

  return 0;
}

/* set key value */
int pthread_setspecific(pthread_key_t key, const void * pointer)
{
  if ( lalDebugLevel & 4 )
    fprintf( stdout, "calling pthread_setspecific(%d,%p)\n", key, pointer );

  if ( key >= PTHREAD_KEYS_MAX || ! pthread_keys[key].in_use )
  {
    return EINVAL;
  }
  pthread_key_values[key] = (void *) pointer;

  return 0;
}

/* get key value */
void * pthread_getspecific(pthread_key_t key)
{
  if ( lalDebugLevel & 4 )
    fprintf( stdout, "calling pthread_getspecific(%d)\n", key );

  if ( key >= PTHREAD_KEYS_MAX || ! pthread_keys[key].in_use )
  {
    return NULL;
  }

  return pthread_key_values[key];
}


/* these functions are needed to resolve the symbols when       */
/* linking in the DFT routines from the Intel Math Kernel       */
/* Library but are not actually called during execution if      */
/* the environment variable MKL_SERIAL=YES. They throw an       */
/* abort() since there is no implementation defined here.       */

int pthread_cancel(pthread_t thread)
{
  fprintf( stderr, "attempt to call pthread_cancel(%ld)\n", thread );
  abort();
  return 0;
}

int pthread_create(pthread_t *thread, const pthread_attr_t *attr,
    void * (*start_routine)(void *), void *arg)
{
  fprintf( stderr, "attempt to call pthread_create(%p,%p,%p,%p)\n",
      thread, attr, start_routine, arg );
  abort();
  return 0;
}

int pthread_join(pthread_t thread_id, void ** thread_return)
{
  fprintf( stderr, "attempt to call pthread_join(%ld,%p)\n",
      thread_id, thread_return );
  abort();
  return 0;
}

int pthread_sigmask(int how, const sigset_t * newmask, sigset_t * oldmask)
{
  fprintf( stderr, "attempt to call pthread_sigmask(%d,%p,%p)\n",
      how, newmask, oldmask );
  abort();
  return 0;
}

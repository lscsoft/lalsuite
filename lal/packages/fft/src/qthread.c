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
  fprintf( stderr, "calling pthread_once(%p,%p)\n", 
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

  fprintf( stderr, "calling pthread_key_create(%p,%p)\n", key, destr );

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
  fprintf( stderr, "calling pthread_key_delete(%d)\n", key );

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
  fprintf( stderr, "calling pthread_setspecific(%d,%p)\n", key, pointer );

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
  fprintf( stderr, "calling pthread_getspecific(%d)\n", key );

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


/* these symbols are unresolved in libguide.a but do not need   */
/* to be resolved for execultables that only use the Intel      */
/* Math Kernel Library DFT routines. They are here for          */
/* completeness, but are commented out since they are not       */
/* needed to compile MKL DFT only code.                         */

#if 0
int pthread_attr_destroy(pthread_attr_t *attr)
{
  fprintf( stderr, "attempt to call pthread_attr_destroy(%p)\n", attr );
  abort();
  return 0;
}

int pthread_attr_getstacksize(const pthread_attr_t *attr, size_t *stacksize)
{
  fprintf( stderr, "attempt to call pthread_attr_getstacksize(%p,%p)\n",
      attr, stacksize );
  abort();
  return 0;
}

int pthread_attr_init(pthread_attr_t *attr)
{
  fprintf( stderr, "attempt to call pthread_attr_init(%p)\n", attr );
  abort();
  return 0;
}

int pthread_attr_setdetachstate(pthread_attr_t *attr, int detachstate)
{
  fprintf( stderr, "attempt to call pthread_attr_setdetachstate(%p,%d)\n",
      attr, detachstate );
  abort();
  return 0;
}

int pthread_attr_setstacksize(pthread_attr_t *attr, size_t stacksize)
{
  fprintf( stderr, "attempt to call pthread_attr_setstacksize(%p,%d)\n",
      attr, stacksize );
  abort();
  return 0;
}

int pthread_cond_destroy(pthread_cond_t *cond)
{
  fprintf( stderr, "attempt to call pthread_cond_destroy(%p)\n",
      cond );
  abort();
  return 0;
}

int pthread_cond_init(pthread_cond_t *cond,
    const pthread_condattr_t *cond_attr)
{
  fprintf( stderr, "attempt to call pthread_cond_init(%p,%p)\n",
      cond, cond_attr );
  abort();
  return 0;
}

int pthread_cond_signal(pthread_cond_t *cond)
{
  fprintf( stderr, "attempt to call pthread_cond_signal(%p)\n",
      cond );
  abort();
  return 0;
}

int pthread_cond_timedwait(pthread_cond_t *cond, pthread_mutex_t *mutex,
    const struct timespec * abstime)
{
  fprintf( stderr, "attempt to call pthread_cond_timedwait(%p,%p,%p)\n",
      cond, mutex, abstime );
  abort();
  return 0;
}

int pthread_cond_wait(pthread_cond_t *cond, pthread_mutex_t *mutex)
{
  fprintf( stderr, "attempt to call pthread_cond_timedwait(%p,%p)\n",
      cond, mutex );
  abort();
  return 0;
}

int pthread_condattr_init(pthread_condattr_t *attr)
{
  fprintf( stderr, "attempt to call pthread_condattr_init(%p)\n", attr );
  abort();
  return 0;
}

int pthread_condattr_destroy(pthread_condattr_t *attr)
{
  fprintf( stderr, "attempt to call pthread_condattr_destroy(%p)\n", attr );
  abort();
  return 0;
}

int pthread_mutex_destroy(pthread_mutex_t * mutex)
{
  fprintf( stderr, "attempt to call pthread_mutex_destroy(%p)\n", mutex );
  abort();
  return 0;
}

int pthread_mutex_init(pthread_mutex_t * mutex,
    const pthread_mutexattr_t * mutex_attr)
{
  fprintf( stderr, "attempt to call pthread_mutex_init(%p,%p)\n", 
      mutex, mutex_attr );
  abort();
  return 0;
}

int pthread_mutex_lock(pthread_mutex_t * mutex)
{
  fprintf( stderr, "attempt to call pthread_mutex_lock(%p)\n", mutex );
  abort();
  return 0;
}

int pthread_mutex_unlock(pthread_mutex_t * mutex)
{
  fprintf( stderr, "attempt to call pthread_mutex_unlock(%p)\n", mutex );
  abort();
  return 0;
}

int pthread_mutexattr_init(pthread_mutexattr_t *attr)
{
  fprintf( stderr, "attempt to call pthread_mutexattr(%p)\n", attr );
  abort();
  return 0;
}

int pthread_mutexattr_destroy(pthread_mutexattr_t *attr)
{
  fprintf( stderr, "attempt to call pthread_mutexattr_destroy(%p)\n", attr );
  abort();
  return 0;
}

pthread_t pthread_self(void)
{
  fprintf( stderr, "attempt to call pthread_self()\n" );
  abort();
  return 0;
}

int pthread_setcancelstate(int state, int * oldstate)
{
  fprintf( stderr, "attempt to call pthread_setcancelstate(%d,%p)\n",
      state, oldstate );
  abort();
  return 0;
}

int pthread_setcanceltype(int type, int * oldtype)
{
  fprintf( stderr, "attempt to call pthread_setcanceltype(%d,%p)\n",
      type, oldtype );
  abort();
  return 0;
}
#endif

/*----------------------------------------------------------------------- 
 * 
 * File Name: Window.c
 * 
 * Author: Allen, B.
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * Window
 * 
 * SYNOPSIS 
 * void Window (Status *, REAL4Vector *vector, WindowParams parameters);
 * 
 * DESCRIPTION 
 * Create a Window function of specified type in vector. 
 * 
 * DIAGNOSTICS 
 * Illegal length, vector == NULL, *vector != NULL, malloc failure
 *
 * CALLS
 * 
 * NOTES
 * 
 *-----------------------------------------------------------------------
 */

#include <math.h>
#include "LALStdlib.h"
#include "Window.h"

NRCSID (WINDOW, "$Id$");

static char *WindowTypeNames[] = WINDOWNAMELIST;

void
Window (Status *status, REAL4Vector *vector, WindowParams *parameters) 
{
  INT4 i;
  INT4 length;
  WindowType windowtype;
  REAL8 wss;    /* window summed and squared */
  REAL8 win;
  REAL8 x,y,z;

  /* Initialize status structure   */
  INITSTATUS(status,"Window",WINDOW);

  /* Check that parameter block is there. */ 
  ASSERT(parameters!=NULL,status,WINDOW_NULLPARAM,WINDOW_MSGNULLPARAM);

  /* check that the vector is not null */
  ASSERT(vector!=NULL,status,WINDOW_NULLHANDLE,WINDOW_MSGNULLHANDLE);

  /* Check that window length is reasonable. */ 
  length=parameters->length;
  ASSERT(length>0,status,WINDOW_ELENGTH,WINDOW_MSGELENGTH);

  /* Make sure that window is of a known type */
  windowtype=parameters->type;
  ASSERT(windowtype>=Rectangular && windowtype<NUMBERWINDOWTYPES,status,
         WINDOW_TYPEUNKNOWN,WINDOW_MSGTYPEUNKNOWN);

  /* vector is apparently already allocated.  Check length, data area */
  ASSERT(vector->length==length,status,
         WINDOW_WRONGLENGTH,WINDOW_MSGWRONGLENGTH);
  ASSERT(vector->data!=NULL,status,WINDOW_NULLDATA,WINDOW_MSGNULLDATA);

  wss=0.0;
  for (i=0;i<length;i++)
  {
    x=(2.0*M_PI*i)/length;
    y=fabs(2.0*i/length-1.0);

    switch (windowtype)
    {
    /* rectangular (no) window */
    case Rectangular:
      win=1.0;
      break;

    /* Hann window */
    case Hann:
      win=0.5*(1.0-cos(x));
      break;

    /* Welch window */
    case Welch:
      win=1.0-y*y;
      break;

    /* Bartlett window */
    case Bartlett:
      win=1.0-y;
      break;

    /* Parzen window */
    case Parzen:
      z=1.0-y;
      if (y<=0.5)
        win=1.0-6.0*y*y*z;
      else
        win=2.0*z*z*z;
      break;

    /* Papoulis window */
    case Papoulis:
      win=1.0/M_PI*sin(M_PI*y)+(1.0-y)*cos(M_PI*y);
      break;

    case Hamming:
      win=1.0-0.46*(1.0+cos(x));
      break;

    /* Default case -- this will NEVER happen -- it is trapped above! */
    default:
      ASSERT(0,status,WINDOW_TYPEUNKNOWN,WINDOW_MSGTYPEUNKNOWN);
      break;
    }
    wss+=win*win;
    vector->data[i]=(REAL4)win;
  }
  parameters->sumofsquares=wss;
  parameters->windowname=WindowTypeNames[parameters->type];

  RETURN(status);
}

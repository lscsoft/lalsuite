/*-----------------------------------------------------------------------
 *
 * File Name: Window.h
 *
 * Author: Allen, Bruce ballen@dirac.phys.uwm.edu
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 * Window.h
 *
 * SYNOPSIS
 * #include "Window.h"
 *
 * DESCRIPTION
 * Error codes, typedefs and prototypes required for use of the Window
 * family of LAL functions
 *
 * DIAGNOSTICS
 *
 *-----------------------------------------------------------------------
 */

#ifndef _WINDOW_H
#define _WINDOW_H

#include "LALDatatypes.h"

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID (WINDOWH, "$Id$");

#define WINDOW_NULLPARAM 1
#define WINDOW_NULLVECTOR 2
#define WINDOW_EALLOCATE 4
#define WINDOW_ELENGTH 8
#define WINDOW_TYPEUNKNOWN 16
#define WINDOW_NULLHANDLE 32
#define WINDOW_WRONGLENGTH 64
#define WINDOW_NULLDATA 128

#define WINDOW_MSGNULLPARAM "null input parameter structure pointer"
#define WINDOW_MSGNULLVECTOR "null output vector pointer"
#define WINDOW_MSGEALLOCATE "unable to allocate vector to store window"
#define WINDOW_MSGELENGTH "length of window is <=0, must be positive"
#define WINDOW_MSGTYPEUNKNOWN "window is of unknown type"
#define WINDOW_MSGNULLHANDLE "input vector is null"
#define WINDOW_MSGWRONGLENGTH "input vector is the wrong length"
#define WINDOW_MSGNULLDATA "data area of input vector is null"

/* Define the types of available windows */
/* WARNING: additional window types must be added just before */
/* NumberWindowTypes, and after the existing window types */
typedef enum
{
    Rectangular,Hann,Welch,Bartlett,Parzen,Papoulis,Hamming,
    /* add any new window types just before this comment */
    NumberWindowTypes
}
WindowType;

#define WINDOWNAMELIST {"Rectangular","Hann","Welch","Bartlett","Parzen","Papoulis","Hamming"}

typedef struct
tagLALWindowParams
{
  INT4		length;		/* length of window */
  WindowType 	type;		/* type of window */
  REAL8   	sumofsquares;	/* sum of window squared  (returned) */
  CHAR*		windowname;	/* pointer to a char string with window name */
}
LALWindowParams;

void LALWindow(Status * ,REAL4Vector *, LALWindowParams *);

#ifdef  __cplusplus
}
#endif

#endif /* _WINDOW_H */

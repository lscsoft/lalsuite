/************************************************* <lalVerbatim file=WindowCV>
Authors: Allen, B., Brown, D. A., and Creighton, T.
$Id$
****************************************************** </lalVerbatim>*/
/****************************************************** <lalLaTeX>

\subsection{Module \texttt{Window.c}}
\label{ss:Window.c}

Creates and destroys a window structure.

\subsubsection*{Prototypes}
\input{WindowCP}
\idx{LALWindow()}
Note that the the function \verb|LALWindow()| is depricated and will soon be
deleted. Windows should be created and destroyed by calles to the 
\verb|LALCreateREAL4Window()| and \verb|LALDestroyREAL4Window()| functions.

\subsubsection*{Description}
These functions create or destroy a time-domain window function in a vector of
specified length.  Note that this function was not written to be
particularly efficient.  If you need a window lots of times, calculate
it once then save it, please.  See the \verb|Window.h| header
discussion for a list of available windows and their definitions.

These window functions are shown in Fig.~\ref{f:window-t}.

****************************************************** </lalLaTeX> */

#include <math.h>
#include <lal/LALConstants.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/Window.h>

NRCSID (WINDOWC, "$Id$");

static const char *WindowTypeNames[] = WINDOWNAMELIST;

static REAL8 BesselI0( REAL8 x )
{
  /*
   * 
   * Calculates the value of the 0th order, modified Bessel function of the 
   * first kind using a power series expansion. (See "Handbook of Math.
   * Functions," eds. Abramowitz and Stegun, 9.6.12 for details.) NOTE: the
   * accuracy of the expansion is chosen to be 2e-9. Stolen from Philip.
   *
   */


  REAL8 ds = 1.0;
  REAL8 d  = 0.0;
  REAL8 s  = 1.0;

  do 
  { 
    d  += 2.0; 
    ds *= x*x/(d*d); 
    s  += ds; 
  } 
  while (ds > 2e-9); 

  return s; 
} 

/* <lalVerbatim file="WindowCP"> */
void 
LALWindow( LALStatus       *status, 
           REAL4Vector     *vector, 
           LALWindowParams *parameters ) 
     /* </lalVerbatim> */
{
  UINT4 i;
  UINT4 length;
  INT4 windowtype;
  REAL8 wss;    /* window summed and squared */
  REAL8 win;
  REAL8 x,y,z;
  REAL8 beta = 0, betaI0 = 0;

  /* Initialize status structure   */
  INITSTATUS(status,"LALWindow Function",WINDOWC);

  /* Check that parameter block is there. */ 
  ASSERT(parameters!=NULL,status,WINDOWH_ENULLPARAM,WINDOWH_MSGENULLPARAM);

  /* check that the vector is not null */
  ASSERT(vector!=NULL,status,WINDOWH_ENULLHANDLE,WINDOWH_MSGENULLHANDLE);

  /* Check that window length is reasonable. */ 
  length=parameters->length;
  ASSERT(length>0,status,WINDOWH_EELENGTH,WINDOWH_MSGEELENGTH);

  /* Make sure that window is of a known type */
  windowtype=parameters->type;
  ASSERT(windowtype>=Rectangular && windowtype<NumberWindowTypes,status,
         WINDOWH_ETYPEUNKNOWN,WINDOWH_MSGETYPEUNKNOWN);

  /* vector is apparently already allocated.  Check length, data area */
  ASSERT(vector->length==length,status,
         WINDOWH_EWRONGLENGTH,WINDOWH_MSGEWRONGLENGTH);
  ASSERT(vector->data!=NULL,status,WINDOWH_ENULLDATA,WINDOWH_MSGENULLDATA);

  /* check that if the case of a Kaiser or Creighton window, beta is
     positive */
  if ( windowtype == Kaiser || windowtype == Creighton )
  {
    ASSERT(parameters->beta >= 0,status, WINDOWH_EBETA,WINDOWH_MSGEBETA);
    beta = parameters->beta;
    if ( windowtype == Kaiser )
      betaI0 = BesselI0( beta );
  }

  wss=0.0;
  for (i=0;i<length;i++)
  {
    x=(2.0*LAL_PI*i)/length;
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
      win=1.0/LAL_PI*sin(LAL_PI*y)+(1.0-y)*cos(LAL_PI*y);
      break;

    case Hamming:
      win=1.0-0.46*(1.0+cos(x));
      break;

    case Kaiser:
      {
        REAL8 kai = (i - (length-1.0)/2.0)*2.0/((REAL8)length-1.0);
        win = BesselI0( beta * sqrt(1.0 - kai*kai) )/betaI0;
      }
      break;

    case Creighton:
      if ( y == 0.0 )
	win = 1.0;
      else if ( y >= 1.0 )
	win = 0.0;
      else
	win = exp( beta/( 1.0 - 1.0/(y*y) ) );
      break;

    /* Default case -- this will NEVER happen -- it is trapped above! */
    default:
      ABORT(status,WINDOWH_ETYPEUNKNOWN,WINDOWH_MSGETYPEUNKNOWN);
      break;
    }
    wss+=win*win;
    vector->data[i]=(REAL4)win;
  }
  parameters->sumofsquares=wss;
  parameters->windowname=WindowTypeNames[parameters->type];

  RETURN(status);
}


/*
 *
 * XLAL Routines.
 *
 */


REAL4Window * XLALCreateREAL4Window( UINT4 length, WindowType type, REAL4 beta )
{
  static const char *func = "XLALCreateREAL4Window";
  REAL8 betaI0 = 0;     /* I_0( beta ) */
  REAL8 dy, pidy, y, w; /* dy=2/N, pidy = pi*dy, y=|j*dy - 1|, w=w(y) */
  REAL8 z;              /* generic intermediate variable */
  REAL8 wss = 0.0;      /* sum of squares of window values */
  REAL4 *data1, *data2; /* pointers to window data */
  REAL4Window *window;
  INT4 i, j;

  if ( ! length )
    XLAL_ERROR_NULL( func, XLAL_EBADLEN );

  if ( type >= NumberWindowTypes )
    XLAL_ERROR_NULL( func, XLAL_ETYPE );
  
  /* Set and check value of beta, if it is used. */
  if ( type == Kaiser || type == Creighton)
  {
    if ( beta < 0 )
      XLAL_ERROR_NULL( func, XLAL_EDOM );
    if ( type == Kaiser )
      betaI0 = BesselI0( beta );
  }

  window = LALCalloc( 1, sizeof( *window ) );
  if ( ! window )
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );

  window->data = XLALCreateREAL4Vector( length );
  if ( ! window->data )
  {
    LALFree( window );
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );
  }

  window->type = type;
  window->beta = beta;
  LALSnprintf( window->windowname, sizeof( window->windowname ), "%s",
      WindowTypeNames[window->type] );

  dy    = 2.0/(REAL4)( length );
  pidy  = LAL_PI*dy;
  data1 = window->data->data;
  data2 = window->data->data + length - 1;

  /* Set initial datum y = -1. */
  switch ( type ) {
  case Rectangular:
    *(data1++) = w = 1.0;
    break;
  case Hamming:
    *(data1++) = w = 0.08;
    break;
  case Kaiser:
    *(data1++) = w = 1.0/betaI0;
    break;
  default:
    *(data1++) = w = 0.0;
    break;
  }
  wss = 0.5*w*w;

  /* Set (symmetric) data for 0 < |y| < 1. */
  j = 1;
  i = ( length - 1 )/2;
  switch ( type ) {

    case Rectangular:
      while ( i-- ) {
        *(data1++) = *(data2--) = 1.0;
        wss += 1.0;
      }
      break;

    case Hann:
      while ( i-- ) {
        *(data1++) = *(data2--) = w = 0.5*( 1.0 - cos( pidy*(j++) ) );
        wss += w*w;
      }
      break;

    case Welch:
      while ( i-- ) {
        y = dy*(j++) - 1.0;
        *(data1++) = *(data2--) = w = 1.0 - y*y;
        wss += w*w;
      }
      break;

    case Bartlett:
      while ( i-- ) {
        *(data1++) = *(data2--) = w = dy*(j++);
        wss += w*w;
      }
      break;

    case Parzen:
      i /= 2;
      while ( i-- ) {
        z = dy*(j++);
        *(data1++) = *(data2--) = w = 2.0*z*z*z;
        wss += w*w;
      }
      i = ( length + 1 )/4;
      while ( i-- ) {
        y = 1.0 - dy*(j++);
        *(data1++) = *(data2--) = w = 1.0 - 6.0*y*y*( 1.0 - y );
        wss += w*w;
      }
      break;

    case Papoulis:
      while ( i-- ) {
        y = 1.0 - dy*(j++);
        z = LAL_PI*y;
        *(data1++) = *(data2--) = w = LAL_1_PI*sin( z )
          + ( 1.0 - y )*cos( z );
        wss += w*w;
      }
      break;

    case Hamming:
      while ( i-- ) {
        *(data1++) = *(data2--) = w = 1.0
          - 0.46*( 1.0 + cos( pidy*(j++) ) );
        wss += w*w;
      }
      break;

    case Kaiser:
      while ( i-- ) {
        y = 1.0 - dy*(j++);
        *(data1++) = *(data2--) = w = BesselI0( beta*sqrt( 1.0 - y*y ) )/betaI0;
        wss += w*w;
      }
      break;

    case Creighton:
      while ( i-- ) {
        y = 1.0 - dy*(j++);
        *(data1++) = *(data2--) = w = exp( beta/( 1.0 - 1.0/(y*y) ) );
        wss += w*w;
      }
      break;

      /* Default case -- this will NEVER happen -- it is trapped above! */
    default:
      XLALDestroyREAL4Vector( window->data );
      LALFree( window );
      XLAL_ERROR_NULL( func, XLAL_ETYPE );
  }

  /* NOTE: At present, all windows are symmetric.  If asymmetric
     windows are ever added, this will need to be changed. */
  wss *= 2.0;

  /* NOTE: At present, all windows are normalized to 1 at y=0.  If
     non-normalized windows are ever added, this will need to be
     changed. */
  if ( data1 == data2 ) {
    *data1 = 1.0;
    wss += 1.0;
  }

  /* Set sum of squares and exit. */
  window->sumofsquares = wss;
  return window;
}


REAL8Window * XLALCreateREAL8Window( UINT4 length, WindowType type, REAL4 beta )
{
  static const char *func = "XLALCreateREAL8Window";
  REAL8 betaI0 = 0;     /* I_0( beta ) */
  REAL8 dy, pidy, y, w; /* dy=2/N, pidy = pi*dy, y=|j*dy - 1|, w=w(y) */
  REAL8 z;              /* generic intermediate variable */
  REAL8 wss = 0.0;      /* sum of squares of window values */
  REAL8 *data1, *data2; /* pointers to window data */
  REAL8Window *window;
  INT4 i, j;    /* length of window, and indecies */

  if ( ! length )
    XLAL_ERROR_NULL( func, XLAL_EBADLEN );

  if ( type >= NumberWindowTypes )
    XLAL_ERROR_NULL( func, XLAL_ETYPE );

  /* Set and check value of beta, if it is used. */
  if ( type == Kaiser || type == Creighton)
  {
    if ( beta < 0 )
      XLAL_ERROR_NULL( func, XLAL_EDOM );
    if ( type == Kaiser )
      betaI0 = BesselI0( beta );
  }
  
  window = LALCalloc( 1, sizeof( *window ) );
  if ( ! window )
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );

  window->data = XLALCreateREAL8Vector( length );
  if ( ! window->data )
  {
    LALFree( window );
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );
  }

  window->type = type;
  window->beta = beta;
  LALSnprintf( window->windowname, sizeof( window->windowname ), "%s",
      WindowTypeNames[window->type] );

  dy    = 2.0/(REAL8)( length );
  pidy  = LAL_PI*dy;
  data1 = window->data->data;
  data2 = window->data->data + length - 1;

  /* Set initial datum y = -1. */
  switch ( type ) {
  case Rectangular:
    *(data1++) = w = 1.0;
    break;
  case Hamming:
    *(data1++) = w = 0.08;
    break;
  case Kaiser:
    *(data1++) = w = 1.0/betaI0;
    break;
  default:
    *(data1++) = w = 0.0;
    break;
  }
  wss = 0.5*w*w;

  /* Set (symmetric) data for 0 < |y| < 1. */
  j = 1;
  i = ( length - 1 )/2;
  switch ( type ) {

    case Rectangular:
      while ( i-- ) {
        *(data1++) = *(data2--) = 1.0;
        wss += 1.0;
      }
      break;

    case Hann:
      while ( i-- ) {
        *(data1++) = *(data2--) = w = 0.5*( 1.0 - cos( pidy*(j++) ) );
        wss += w*w;
      }
      break;

    case Welch:
      while ( i-- ) {
        y = dy*(j++) - 1.0;
        *(data1++) = *(data2--) = w = 1.0 - y*y;
        wss += w*w;
      }
      break;

    case Bartlett:
      while ( i-- ) {
        *(data1++) = *(data2--) = w = dy*(j++);
        wss += w*w;
      }
      break;

    case Parzen:
      i /= 2;
      while ( i-- ) {
        z = dy*(j++);
        *(data1++) = *(data2--) = w = 2.0*z*z*z;
        wss += w*w;
      }
      i = ( length + 1 )/4;
      while ( i-- ) {
        y = 1.0 - dy*(j++);
        *(data1++) = *(data2--) = w = 1.0 - 6.0*y*y*( 1.0 - y );
        wss += w*w;
      }
      break;

    case Papoulis:
      while ( i-- ) {
        y = 1.0 - dy*(j++);
        z = LAL_PI*y;
        *(data1++) = *(data2--) = w = LAL_1_PI*sin( z )
          + ( 1.0 - y )*cos( z );
        wss += w*w;
      }
      break;

    case Hamming:
      while ( i-- ) {
        *(data1++) = *(data2--) = w = 1.0
          - 0.46*( 1.0 + cos( pidy*(j++) ) );
        wss += w*w;
      }
      break;

    case Kaiser:
      while ( i-- ) {
        y = 1.0 - dy*(j++);
        *(data1++) = *(data2--) = w = BesselI0( beta*sqrt( 1.0 - y*y ) )/betaI0;
        wss += w*w;
      }
      break;

    case Creighton:
      while ( i-- ) {
        y = 1.0 - dy*(j++);
        *(data1++) = *(data2--) = w = exp( beta/( 1.0 - 1.0/(y*y) ) );
        wss += w*w;
      }
      break;

      /* Default case -- this will NEVER happen -- it is trapped above! */
    default:
      XLALDestroyREAL8Vector( window->data );
      LALFree( window );
      XLAL_ERROR_NULL( func, XLAL_ETYPE );
  }

  /* NOTE: At present, all windows are symmetric.  If asymmetric
     windows are ever added, this will need to be changed. */
  wss *= 2.0;

  /* NOTE: At present, all windows are normalized to 1 at y=0.  If
     non-normalized windows are ever added, this will need to be
     changed. */
  if ( data1 == data2 ) {
    *data1 = 1.0;
    wss += 1.0;
  }

  /* Set sum of squares and exit. */
  window->sumofsquares = wss;
  return window;
}


REAL4Window * XLALCreateRectangularREAL4Window( UINT4 length )
{
  static const char *func = "XLALCreateRectangularREAL4Window";
  REAL4Window *window;
  window = XLALCreateREAL4Window( length, Rectangular, 0 );
  if ( ! window )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return window;
}

REAL4Window * XLALCreateHannREAL4Window( UINT4 length )
{
  static const char *func = "XLALCreateHannREAL4Window";
  REAL4Window *window;
  window = XLALCreateREAL4Window( length, Hann, 0 );
  if ( ! window )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return window;
}

REAL4Window * XLALCreateWelchREAL4Window( UINT4 length )
{
  static const char *func = "XLALCreateWelchREAL4Window";
  REAL4Window *window;
  window = XLALCreateREAL4Window( length, Welch, 0 );
  if ( ! window )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return window;
}

REAL4Window * XLALCreateBartlettREAL4Window( UINT4 length )
{
  static const char *func = "XLALCreateBartlettREAL4Window";
  REAL4Window *window;
  window = XLALCreateREAL4Window( length, Bartlett, 0 );
  if ( ! window )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return window;
}

REAL4Window * XLALCreateParzenREAL4Window( UINT4 length )
{
  static const char *func = "XLALCreateParzenREAL4Window";
  REAL4Window *window;
  window = XLALCreateREAL4Window( length, Parzen, 0 );
  if ( ! window )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return window;
}

REAL4Window * XLALCreatePapoulisREAL4Window( UINT4 length )
{
  static const char *func = "XLALCreatePapoulisREAL4Window";
  REAL4Window *window;
  window = XLALCreateREAL4Window( length, Papoulis, 0 );
  if ( ! window )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return window;
}

REAL4Window * XLALCreateHammingREAL4Window( UINT4 length )
{
  static const char *func = "XLALCreateHammingREAL4Window";
  REAL4Window *window;
  window = XLALCreateREAL4Window( length, Hamming, 0 );
  if ( ! window )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return window;
}

REAL4Window * XLALCreateKaiserREAL4Window( UINT4 length, REAL4 beta )
{
  static const char *func = "XLALCreateKaiserREAL4Window";
  REAL4Window *window;
  window = XLALCreateREAL4Window( length, Kaiser, beta );
  if ( ! window )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return window;
}

REAL4Window * XLALCreateCreightonREAL4Window( UINT4 length, REAL4 beta )
{
  static const char *func = "XLALCreateCreightonREAL4Window";
  REAL4Window *window;
  window = XLALCreateREAL4Window( length, Creighton, beta );
  if ( ! window )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return window;
}


REAL8Window * XLALCreateRectangularREAL8Window( UINT4 length )
{
  static const char *func = "XLALCreateRectangularREAL8Window";
  REAL8Window *window;
  window = XLALCreateREAL8Window( length, Rectangular, 0 );
  if ( ! window )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return window;
}

REAL8Window * XLALCreateHannREAL8Window( UINT4 length )
{
  static const char *func = "XLALCreateHannREAL8Window";
  REAL8Window *window;
  window = XLALCreateREAL8Window( length, Hann, 0 );
  if ( ! window )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return window;
}

REAL8Window * XLALCreateWelchREAL8Window( UINT4 length )
{
  static const char *func = "XLALCreateWelchREAL8Window";
  REAL8Window *window;
  window = XLALCreateREAL8Window( length, Welch, 0 );
  if ( ! window )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return window;
}

REAL8Window * XLALCreateBartlettREAL8Window( UINT4 length )
{
  static const char *func = "XLALCreateBartlettREAL8Window";
  REAL8Window *window;
  window = XLALCreateREAL8Window( length, Bartlett, 0 );
  if ( ! window )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return window;
}

REAL8Window * XLALCreateParzenREAL8Window( UINT4 length )
{
  static const char *func = "XLALCreateParzenREAL8Window";
  REAL8Window *window;
  window = XLALCreateREAL8Window( length, Parzen, 0 );
  if ( ! window )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return window;
}

REAL8Window * XLALCreatePapoulisREAL8Window( UINT4 length )
{
  static const char *func = "XLALCreatePapoulisREAL8Window";
  REAL8Window *window;
  window = XLALCreateREAL8Window( length, Papoulis, 0 );
  if ( ! window )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return window;
}

REAL8Window * XLALCreateHammingREAL8Window( UINT4 length )
{
  static const char *func = "XLALCreateHammingREAL8Window";
  REAL8Window *window;
  window = XLALCreateREAL8Window( length, Hamming, 0 );
  if ( ! window )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return window;
}

REAL8Window * XLALCreateKaiserREAL8Window( UINT4 length, REAL4 beta )
{
  static const char *func = "XLALCreateKaiserREAL8Window";
  REAL8Window *window;
  window = XLALCreateREAL8Window( length, Kaiser, beta );
  if ( ! window )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return window;
}

REAL8Window * XLALCreateCreightonREAL8Window( UINT4 length, REAL4 beta )
{
  static const char *func = "XLALCreateCreightonREAL8Window";
  REAL8Window *window;
  window = XLALCreateREAL8Window( length, Creighton, beta );
  if ( ! window )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return window;
}


void XLALDestroyREAL4Window( REAL4Window *window )
{
  if ( window )
  {
    if ( window->data )
      XLALDestroyREAL4Vector( window->data );
    LALFree( window );
  }
  return;
}


void XLALDestroyREAL8Window( REAL8Window *window )
{
  if ( window )
  {
    if ( window->data )
      XLALDestroyREAL8Vector( window->data );
    LALFree( window );
  }
  return;
}



/*
 *
 * LAL Routines.
 *
 */




/* <lalVerbatim file="WindowCP"> */
void LALCreateREAL4Window ( LALStatus       *status,
			    REAL4Window     **output,
			    LALWindowParams *params )
{ /* </lalVerbatim> */

  /*
   *
   * FIXME: should call the XLAL routine!!!
   *
   */

  INT4 length, i, j;    /* length of window, and indecies */
  REAL8 beta, betaI0;   /* beta parameter, and I_0( beta ) */
  REAL8 dy, pidy, y, w; /* dy=2/N, pidy = pi*dy, y=|j*dy - 1|, w=w(y) */
  REAL8 z;              /* generic intermediate variable */
  REAL8 wss = 0.0;      /* sum of squares of window values */
  REAL4 *data1, *data2; /* pointers to window data */
  WindowType type;      /* window type */

  INITSTATUS( status, "LALCreateREAL4Window", WINDOWC );
  ATTATCHSTATUSPTR( status );

  /* Dummy assignments so compiler won't complain. */
  beta = betaI0 = 0.0;

  /* Check input parameters. */
  ASSERT( output, status, WINDOWH_ENULLHANDLE, WINDOWH_MSGENULLHANDLE );
  ASSERT( !*output, status, WINDOWH_ENNUL, WINDOWH_MSGENNUL );
  ASSERT( params, status, WINDOWH_ENULLPARAM, WINDOWH_MSGENULLPARAM );
  type = params->type;
  ASSERT( type < NumberWindowTypes, status, WINDOWH_ETYPEUNKNOWN,
	  WINDOWH_MSGETYPEUNKNOWN );
  ASSERT( (int)type >= 0, status, WINDOWH_ETYPEUNKNOWN,
	  WINDOWH_MSGETYPEUNKNOWN );
  length = params->length;
  ASSERT( length > 0, status, WINDOWH_EELENGTH, WINDOWH_MSGEELENGTH );
  dy = 2.0/(REAL4)( length );
  pidy = LAL_PI*dy;

  /* Set and check value of beta, if it is used. */
  if ( type == Kaiser || type == Creighton) {
    beta = params->beta;
    ASSERT( beta >= 0, status, WINDOWH_EBETA, WINDOWH_MSGEBETA );
    if ( type == Kaiser )
      betaI0 = BesselI0( beta );
  }

  /* Allocate the storage for the window vector */
  *output = (REAL4Window *) LALCalloc( 1, sizeof(REAL4Window) );
  if ( ! *output ) {
    ABORT( status, WINDOWH_EEALLOCATE, WINDOWH_MSGEEALLOCATE );
  }
  LALSCreateVector( status->statusPtr, &((*output)->data), length );
  BEGINFAIL( status ) {
    LALFree( *output );
    *output = NULL;
  } ENDFAIL( status );

  /* Create window. */
  (*output)->type = type;
  (*output)->beta = beta;
  LALSnprintf( (*output)->windowname, LALNameLength, "%s",
	       WindowTypeNames[type] );
  params->windowname = WindowTypeNames[type];
  data1 = (*output)->data->data;
  data2 = (*output)->data->data + length - 1;

  /* Set initial datum y = -1. */
  switch ( type ) {
  case Rectangular:
    *(data1++) = w = 1.0;
    break;
  case Hamming:
    *(data1++) = w = 0.08;
    break;
  case Kaiser:
    *(data1++) = w = 1.0/betaI0;
    break;
  default:
    *(data1++) = w = 0.0;
    break;
  }
  wss = 0.5*w*w;

  /* Set (symmetric) data for 0 < |y| < 1. */
  j = 1;
  i = ( length - 1 )/2;
  switch ( type ) {

  case Rectangular:
    while ( i-- ) {
      *(data1++) = *(data2--) = 1.0;
      wss += 1.0;
    }
    break;

  case Hann:
    while ( i-- ) {
      *(data1++) = *(data2--) = w = 0.5*( 1.0 - cos( pidy*(j++) ) );
      wss += w*w;
    }
    break;

  case Welch:
    while ( i-- ) {
      y = dy*(j++) - 1.0;
      *(data1++) = *(data2--) = w = 1.0 - y*y;
      wss += w*w;
    }
    break;

  case Bartlett:
    while ( i-- ) {
      *(data1++) = *(data2--) = w = dy*(j++);
      wss += w*w;
    }
    break;

  case Parzen:
    i /= 2;
    while ( i-- ) {
      z = dy*(j++);
      *(data1++) = *(data2--) = w = 2.0*z*z*z;
      wss += w*w;
    }
    i = ( length + 1 )/4;
    while ( i-- ) {
      y = 1.0 - dy*(j++);
      *(data1++) = *(data2--) = w = 1.0 - 6.0*y*y*( 1.0 - y );
      wss += w*w;
    }
    break;

  case Papoulis:
    while ( i-- ) {
      y = 1.0 - dy*(j++);
      z = LAL_PI*y;
      *(data1++) = *(data2--) = w = LAL_1_PI*sin( z )
	+ ( 1.0 - y )*cos( z );
      wss += w*w;
    }
    break;

  case Hamming:
    while ( i-- ) {
      *(data1++) = *(data2--) = w = 1.0
	- 0.46*( 1.0 + cos( pidy*(j++) ) );
      wss += w*w;
    }
    break;

  case Kaiser:
    while ( i-- ) {
      y = 1.0 - dy*(j++);
      *(data1++) = *(data2--) = w = BesselI0( beta*sqrt( 1.0 - y*y ) )
	/betaI0;
      wss += w*w;
    }
    break;

  case Creighton:
    while ( i-- ) {
      y = 1.0 - dy*(j++);
      *(data1++) = *(data2--) = w = exp( beta/( 1.0 - 1.0/(y*y) ) );
      wss += w*w;
    }
    break;

    /* Default case -- this will NEVER happen -- it is trapped above! */
  default:
    TRY( LALSDestroyVector( status, &((*output)->data) ), status );
    LALFree( *output );
    *output = NULL;
    ABORT( status, WINDOWH_ETYPEUNKNOWN, WINDOWH_MSGETYPEUNKNOWN );
    break;
  }

  /* NOTE: At present, all windows are symmetric.  If asymmetric
     windows are ever added, this will need to be changed. */
  wss *= 2.0;

  /* NOTE: At present, all windows are normalized to 1 at y=0.  If
     non-normalized windows are ever added, this will need to be
     changed. */
  if ( data1 == data2 ) {
    *data1 = 1.0;
    wss += 1.0;
  }

  /* Set sum of squares and exit. */
  (*output)->sumofsquares = params->sumofsquares = wss;
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="WindowCP"> */
void LALCreateREAL8Window ( LALStatus       *status,
			    REAL8Window     **output,
			    LALWindowParams *params )
{ /* </lalVerbatim> */

  /*
   *
   * FIXME: should call the XLAL routine!!!
   *
   */

  INT4 length, i, j;    /* length of window, and indecies */
  REAL8 beta, betaI0;   /* beta parameter, and I_0( beta ) */
  REAL8 dy, pidy, y, w; /* dy=2/N, pidy = pi*dy, y=|j*dy - 1|, w=w(y) */
  REAL8 z;              /* generic intermediate variable */
  REAL8 wss = 0.0;      /* sum of squares of window values */
  REAL8 *data1, *data2; /* pointers to window data */
  WindowType type;      /* window type */

  INITSTATUS( status, "LALCreateREAL8Window", WINDOWC );
  ATTATCHSTATUSPTR( status );

  /* Dummy assignments so compiler won't complain. */
  beta = betaI0 = 0.0;

  /* Check input parameters. */
  ASSERT( output, status, WINDOWH_ENULLHANDLE, WINDOWH_MSGENULLHANDLE );
  ASSERT( !*output, status, WINDOWH_ENNUL, WINDOWH_MSGENNUL );
  ASSERT( params, status, WINDOWH_ENULLPARAM, WINDOWH_MSGENULLPARAM );
  type = params->type;
  ASSERT( type < NumberWindowTypes, status, WINDOWH_ETYPEUNKNOWN,
	  WINDOWH_MSGETYPEUNKNOWN );
  length = params->length;
  ASSERT( length > 0, status, WINDOWH_EELENGTH, WINDOWH_MSGEELENGTH );
  dy = 2.0/(REAL8)( length );
  pidy = LAL_PI*dy;

  /* Set and check value of beta, if it is used. */
  if ( type == Kaiser || type == Creighton) {
    beta = params->beta;
    ASSERT( beta >= 0, status, WINDOWH_EBETA, WINDOWH_MSGEBETA );
    if ( type == Kaiser )
      betaI0 = BesselI0( beta );
  }

  /* Allocate the storage for the window vector */
  *output = (REAL8Window *) LALCalloc( 1, sizeof(REAL8Window) );
  if ( ! *output ) {
    ABORT( status, WINDOWH_EEALLOCATE, WINDOWH_MSGEEALLOCATE );
  }
  LALDCreateVector( status->statusPtr, &((*output)->data), length );
  BEGINFAIL( status ) {
    LALFree( *output );
    *output = NULL;
  } ENDFAIL( status );

  /* Create window. */
  (*output)->type = type;
  (*output)->beta = beta;
  LALSnprintf( (*output)->windowname, LALNameLength, "%s",
	       WindowTypeNames[type] );
  params->windowname = WindowTypeNames[type];
  data1 = (*output)->data->data;
  data2 = (*output)->data->data + length - 1;

  /* Set initial datum y = -1. */
  switch ( type ) {
  case Rectangular:
    *(data1++) = w = 1.0;
    break;
  case Hamming:
    *(data1++) = w = 0.08;
    break;
  case Kaiser:
    *(data1++) = w = 1.0/betaI0;
    break;
  default:
    *(data1++) = w = 0.0;
    break;
  }
  wss = 0.5*w*w;

  /* Set (symmetric) data for 0 < |y| < 1. */
  j = 1;
  i = ( length - 1 )/2;
  switch ( type ) {

  case Rectangular:
    while ( i-- ) {
      *(data1++) = *(data2--) = 1.0;
      wss += 1.0;
    }
    break;

  case Hann:
    while ( i-- ) {
      *(data1++) = *(data2--) = w = 0.5*( 1.0 - cos( pidy*(j++) ) );
      wss += w*w;
    }
    break;

  case Welch:
    while ( i-- ) {
      y = dy*(j++) - 1.0;
      *(data1++) = *(data2--) = w = 1.0 - y*y;
      wss += w*w;
    }
    break;

  case Bartlett:
    while ( i-- ) {
      *(data1++) = *(data2--) = w = dy*(j++);
      wss += w*w;
    }
    break;

  case Parzen:
    i /= 2;
    while ( i-- ) {
      z = dy*(j++);
      *(data1++) = *(data2--) = w = 2.0*z*z*z;
      wss += w*w;
    }
    i = ( length + 1 )/4;
    while ( i-- ) {
      y = 1.0 - dy*(j++);
      *(data1++) = *(data2--) = w = 1.0 - 6.0*y*y*( 1.0 - y );
      wss += w*w;
    }
    break;

  case Papoulis:
    while ( i-- ) {
      y = 1.0 - dy*(j++);
      z = LAL_PI*y;
      *(data1++) = *(data2--) = w = LAL_1_PI*sin( z )
	+ ( 1.0 - y )*cos( z );
      wss += w*w;
    }
    break;

  case Hamming:
    while ( i-- ) {
      *(data1++) = *(data2--) = w = 1.0
	- 0.46*( 1.0 + cos( pidy*(j++) ) );
      wss += w*w;
    }
    break;

  case Kaiser:
    while ( i-- ) {
      y = 1.0 - dy*(j++);
      *(data1++) = *(data2--) = w = BesselI0( beta*sqrt( 1.0 - y*y ) )
	/betaI0;
      wss += w*w;
    }
    break;

  case Creighton:
    while ( i-- ) {
      y = 1.0 - dy*(j++);
      *(data1++) = *(data2--) = w = exp( beta/( 1.0 - 1.0/(y*y) ) );
      wss += w*w;
    }
    break;

    /* Default case -- this will NEVER happen -- it is trapped above! */
  default:
    TRY( LALDDestroyVector( status, &((*output)->data) ), status );
    LALFree( *output );
    *output = NULL;
    ABORT( status, WINDOWH_ETYPEUNKNOWN, WINDOWH_MSGETYPEUNKNOWN );
    break;
  }

  /* NOTE: At present, all windows are symmetric.  If asymmetric
     windows are ever added, this will need to be changed. */
  wss *= 2.0;

  /* NOTE: At present, all windows are normalized to 1 at y=0.  If
     non-normalized windows are ever added, this will need to be
     changed. */
  if ( data1 == data2 ) {
    *data1 = 1.0;
    wss += 1.0;
  }

  /* Set sum of squares and exit. */
  (*output)->sumofsquares = params->sumofsquares = wss;
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="WindowCP"> */
void LALDestroyREAL4Window ( LALStatus *status, REAL4Window **output )
/* </lalVerbatim> */
{

  /*
   *
   * FIXME: should call the XLAL routine!!!
   *
   */

  INITSTATUS( status, "LALDestroyREAL4Window", WINDOWC );
  ATTATCHSTATUSPTR( status );

  ASSERT( output, status, WINDOWH_ENULL, WINDOWH_MSGENULL );
  ASSERT( *output, status, WINDOWH_ENULL, WINDOWH_MSGENULL );

  /* destroy the window */
  LALSDestroyVector( status->statusPtr, &((*output)->data) );
  CHECKSTATUSPTR( status );
  LALFree( *output );
  *output = NULL;

  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="WindowCP"> */
void LALDestroyREAL8Window ( LALStatus *status, REAL8Window **output )
/* </lalVerbatim> */
{

  /*
   *
   * FIXME: should call the XLAL routine!!!
   *
   */

  INITSTATUS( status, "LALDestroyREAL8Window", WINDOWC );
  ATTATCHSTATUSPTR( status );

  ASSERT( output, status, WINDOWH_ENULL, WINDOWH_MSGENULL );
  ASSERT( *output, status, WINDOWH_ENULL, WINDOWH_MSGENULL );

  /* destroy the window */
  LALDDestroyVector( status->statusPtr, &((*output)->data) );
  CHECKSTATUSPTR( status );
  LALFree( *output );
  *output = NULL;

  DETATCHSTATUSPTR( status );
  RETURN( status );
}

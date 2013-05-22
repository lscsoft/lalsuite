/*
*  Copyright (C) 2007 Jolien Creighton, Teviet Creighton
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

#include <complex.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <math.h>
#include <lal/IIRFilter.h>

/**
   \addtogroup CreateIIRFilter_c
   \author Creighton, T. D.

   \brief Creates IIR filter objects.

\heading{Description}

These functions create an object <tt>**output</tt> of type
<tt>\<datatype\>IIRFilter</tt>, where <tt>\<datatype\></tt> is \c REAL4 or
\c REAL8.  The filter coefficients are computed from the zeroes,
poles, and gain of an input object <tt>*input</tt> of type
\c COMPLEX8ZPGFilter or \c COMPLEX16ZPGFilter, respectively;
the sampling time interval is taken directly from
<tt>input->deltaT</tt>.  The ZPG filter should express the factored
transfer function in the \f$z=\exp(2\pi if)\f$ plane.  Initially the
output handle must be a valid handle (\c output\f$\neq\f$\c NULL)
but should not point to an existing object
(<tt>*output</tt>=\c NULL)

\heading{Algorithm}

An IIR filter is a real time-domain filter, which imposes certain
constraints on the zeros, poles, and gain of the filter transfer
function.  The function <tt>Create\<datatype\>IIRFilter()</tt> deals with
the constraints either by aborting if they are not met, or by
adjusting the filter response so that they are met.  In the latter
case, warning messages will be issued if the external parameter
\c lalDebugLevel is set to allow such messages.  The specific
constraints, and how they are dealt with, are as follows:

First, the filter must be \e causal; that is, the output at any
time can be generated entirely from the input at previous times.  In
practice this means that the number of (finite) poles in the \f$z\f$ plane
must equal or exceed the number of (finite) zeros.  If this is not the
case, <tt>Create\<datatype\>IIRFilter()</tt> will add additional poles at
\f$z=0\f$.  Effectively this is just adding a delay to the filter response
in order to make it causal.

Second, the filter should be \e stable, which means that all poles
should be located on or within the circle \f$|z|=1\f$.  This is not
enforced by <tt>Create\<datatype\>IIRFilter()</tt>, which can be used to
make unstable filters; however, warnings will be issued.  (In some
sense the first condition is a special case of this one, since a
transfer function with more zeros than poles actually has
corresponding poles at infinity.)

Third, the filter must be <em>physically realizable</em>; that is, the
transfer function should expand to a rational function of \f$z\f$ with
real coefficients.  Necessary and sufficient conditions for this are
that the gain be real, and that all zeros and poles either be real or
come in complex conjugate pairs.  The routine
<tt>Create\<datatype\>IIRFilter()</tt> enforces this by using only the
real part of the gain, and only the real or positive-imaginary zeros
and poles; it assumes that the latter are paired with
negative-imaginary conjugates.  The routine will abort if this
assumption results in a change in the given number of zeros or poles.
If \c lalDebugLevel is set to allow warnings, the routine will
actually check to see that each pair of nonreal poles or zeros are in
fact complex conjugates, and will issue a warning if an unmatched pair
is detected; however, the algorithm will then simply proceed as if the
negative-imaginary points were relocated to the "correct" positions.

The code associated with the warning messages is potentially rather
cumbersome for production algorithms; therefore, the value of
\c lalDebugLevel is tested before performing any other tests
associated with warning messages.  Furthermore, this code block is
surrounded with compiler directives to exclude the code entirely if
the module is compiled with the \c NDEBUG flag set.

*/
/*@{*/

/** \see See \ref CreateIIRFilter_c for documentation */
REAL4IIRFilter *XLALCreateREAL4IIRFilter( COMPLEX8ZPGFilter *input )
{
  REAL4IIRFilter *output;
  INT4 i;          /* Index counter for zeros and poles. */
  INT4 numZeros;   /* The number of zeros. */
  INT4 numPoles;   /* The number of poles. */
  INT4 numDirect;  /* The number of direct filter coefficients. */
  INT4 numRecurs;  /* The number of recursive filter coefficients. */
  INT4 num;        /* An extra counter for error checking. */
  COMPLEX8 *zeros; /* The zeros of the transfer function. */
  COMPLEX8 *poles; /* The poles of the transfer function. */
  REAL4 *direct;   /* The direct filter coefficients. */
  REAL4 *recurs;   /* The recursive filter coefficients. */
  REAL4 *history;  /* The filter history. */

  /* Make sure all the input structures have been initialized. */
  if ( ! input )
    XLAL_ERROR_NULL( XLAL_EFAULT );
  if ( ! input->zeros || ! input->poles
      || ! input->zeros->data || ! input->poles->data )
    XLAL_ERROR_NULL( XLAL_EINVAL );

  numZeros=input->zeros->length;
  numPoles=input->poles->length;
  zeros=input->zeros->data;
  poles=input->poles->data;

  /* Check that zeros are appropriately paired.  Also, keep track of
     the number of zeros at z=0, since these will reduce the number of
     direct coefficients required. */
  numDirect=1;
  for(i=0,num=0;i<numZeros;i++)
    if(cimagf(zeros[i])==0.0){
      num+=1;
      if(crealf(zeros[i])==0.0)
	numDirect-=1;
    }
    else if(cimagf(zeros[i])>0.0){
      num+=2;

#ifndef NDEBUG
      if(lalDebugLevel&LALWARNING){
	/* Check to see that another zero is an actual conjugate.
           This is not a foolproof test, as it can be fooled by
           multiple zeros at the same location. */
	INT4 j=0;
	INT4 k=0;
	REAL4 x = crealf(zeros[i]) - crealf(zeros[0]);
	REAL4 y = cimagf(zeros[i]) + cimagf(zeros[0]);
	REAL4 sep = x*x + y*y;
	for(j=1;j<numZeros;j++){
	  x=crealf(zeros[i])-crealf(zeros[j]);
	  y=cimagf(zeros[i])+cimagf(zeros[j]);
	  if(sep>x*x+y*y){
	    sep=x*x+y*y;
	    k=j;
	  }
	}
	if(sep>LAL_REAL4_EPS){
          XLALPrintWarning( "XLAL Warning - %s: ", __func__ );
          XLALPrintWarning("Complex zero has no conjugate pair\n");
	  XLALPrintWarning("\tUnmatched zero z_%i = %.8e + i*%.8e\n",i,
              crealf(zeros[i]),cimagf(zeros[i]));
	  XLALPrintWarning("\tNearest pair   z_%i = %.8e + i*%.8e\n",k,
              crealf(zeros[k]),cimagf(zeros[k]));
	}
      }
#endif
    }

  if ( num != numZeros )
    XLAL_ERROR_NULL( XLAL_EINVAL );

  /* Check that poles are appropriately paired.  Also, keep track of
     the number of poles at z=0, since these will reduce the number of
     recursive coefficients required. */
  numRecurs=1;
  for(i=0,num=0;i<numPoles;i++)
    if(cimagf(poles[i])==0.0){
      num+=1;
      if(crealf(poles[i])==0.0)
	numRecurs-=1;
    }
    else if(cimagf(poles[i])>0.0){
      num+=2;

#ifndef NDEBUG
      if(lalDebugLevel&LALWARNING){
	/* Check to see that another pole is an actual conjugate.
           This is not a foolproof test, as it can be fooled by
           multiple poles at the same location. */
	INT4 j=0;
	INT4 k=0;
	REAL4 x = crealf(poles[i]) - crealf(poles[0]);
	REAL4 y = cimagf(poles[i]) + cimagf(poles[0]);
	REAL4 sep = x*x + y*y;
	for(j=1;j<numPoles;j++){
	  x=crealf(poles[i])-crealf(poles[j]);
	  y=cimagf(poles[i])+cimagf(poles[j]);
	  if(sep>x*x+y*y){
	    sep=x*x+y*y;
	    k=j;
	  }
	}
	if(sep>LAL_REAL4_EPS){
          XLALPrintWarning( "XLAL Warning - %s: ", __func__ );
          XLALPrintWarning("Complex pole has no conjugate pair\n");
	  XLALPrintWarning("\tUnmatched pole p_%i = %.8e + i*%.8e\n",i,
              crealf(poles[i]),cimagf(poles[i]));
	  XLALPrintWarning("\tNearest pair   p_%i = %.8e + i*%.8e\n",k,
              crealf(poles[k]),cimagf(poles[k]));
	}
      }
#endif
    }

  if ( num != numPoles )
    XLAL_ERROR_NULL( XLAL_EINVAL );

#ifndef NDEBUG
  if(lalDebugLevel&LALWARNING){
    /* Issue a warning if the gain is nonreal. */
    if(fabs(cimag(input->gain))>fabs(LAL_REAL4_EPS*creal(input->gain))){
      XLALPrintWarning( "XLAL Warning - %s: ", __func__ );
      XLALPrintWarning("Gain is non-real\n");
      XLALPrintWarning("\tg = %.8e + i*%.8e\n", crealf(input->gain), cimagf(input->gain));
    }
    /* Issue a warning if there are any ``removeable'' poles. */
    for(i=0;i<numPoles;i++){
      INT4 j=0;
      for(;j<numZeros;j++)
	if((crealf(poles[i])==crealf(zeros[j]))&&(cimagf(poles[i])==cimagf(zeros[j]))){
          XLALPrintWarning( "XLAL Warning - %s: ", __func__ );
	  XLALPrintWarning("Removeable pole\n");
	  XLALPrintWarning("\tp_%i = z_%i = %.8e + i*%.8e\n",i,j,
              crealf(poles[i]),cimagf(poles[i]));
	}
    }
    /* Issue a warning if extra factors of 1/z will be applied. */
    if(numPoles<numZeros){
      XLALPrintWarning( "XLAL Warning - %s: ", __func__ );
      XLALPrintWarning("Filter has more zeros than poles\n");
      XLALPrintWarning("\t%i poles added at complex origin\n",
          numZeros-numPoles);
    }
    /* Issue a warning if any poles are outside |z|=1. */
    for(i=0;i<numPoles;i++){
      REAL4 zAbs=crealf(poles[i])*crealf(poles[i])+cimagf(poles[i])*cimagf(poles[i]);
      if(zAbs>1.0){
        XLALPrintWarning( "XLAL Warning - %s: ", __func__ );
	XLALPrintWarning("Filter has pole outside of unit circle\n");
	XLALPrintWarning("\tp_%i = %.8e + i*%.8e, |p_%i| = %.8e\n",i,
		      crealf(poles[i]),cimagf(poles[i]),i,zAbs);
      }
    }
  }
#endif

  /* Everything seems okay, so initialize the filter. */
  output=LALCalloc(1,sizeof(*output));
  if ( ! output )
    XLAL_ERROR_NULL( XLAL_ENOMEM );
  num = (numPoles>=numZeros) ? numZeros : numPoles;
  numDirect+=num;
  numRecurs+=num;

  output->deltaT=input->deltaT;
  output->directCoef = XLALCreateREAL4Vector( numDirect );
  output->recursCoef = XLALCreateREAL4Vector( numRecurs );
  if ( ! output->directCoef || ! output->recursCoef )
  {
    XLALDestroyREAL4IIRFilter( output );
    XLAL_ERROR_NULL( XLAL_EFUNC );
  }
  direct=output->directCoef->data;
  recurs=output->recursCoef->data;

  /* Expand the denominator as a polynomial in z. */
  *recurs=-1.0;
  for(i=1;i<numRecurs;i++)
    recurs[i]=0.0;
  for(i=0;i<numPoles;i++){
    INT4 j=numRecurs-1;
    REAL4 x=crealf(poles[i]);
    REAL4 y=cimagf(poles[i]);
    if(y==0.0)
      for(;j>0;j--)
	recurs[j]-=x*recurs[j-1];
    else if(y>0.0){
      for(;j>1;j--)
	recurs[j]-=2.0*x*recurs[j-1]-(x*x+y*y)*recurs[j-2];
      recurs[j]-=2.0*x*recurs[j-1];
    }
  }

  /* Expand the numerator as a polynomial in z. */
  for(i=0;i<numDirect;i++)
    direct[i]=0.0;
  direct[num-numZeros]=crealf(input->gain);
  for(i=0;i<numZeros;i++){
    INT4 j=numDirect-1;
    REAL4 x=crealf(zeros[i]);
    REAL4 y=cimagf(zeros[i]);
    if(y==0.0)
      for(;j>num-numZeros;j--)
	direct[j]-=x*direct[j-1];
    else if(y>0.0){
      for(;j>num-numZeros+1;j--)
	direct[j]-=2.0*x*direct[j-1]-(x*x+y*y)*direct[j-2];
      direct[j]-=2.0*x*direct[j-1];
    }
  }

  /* Initialize filter history. */
  if(numDirect>=numRecurs)
    num=numDirect-1;
  else
    num=numRecurs-1;
  output->history = XLALCreateREAL4Vector( numRecurs );
  if ( ! output->history )
  {
    XLALDestroyREAL4IIRFilter( output );
    XLAL_ERROR_NULL( XLAL_EFUNC );
  }
  history=output->history->data;
  for(i=0;i<num;i++)
    history[i]=0.0;

  /* Normal exit */
  return output;
}

/** \see See \ref CreateIIRFilter_c for documentation */
REAL8IIRFilter *XLALCreateREAL8IIRFilter( COMPLEX16ZPGFilter *input )
{
  REAL8IIRFilter *output;
  INT4 i;          /* Index counter for zeros and poles. */
  INT4 numZeros;   /* The number of zeros. */
  INT4 numPoles;   /* The number of poles. */
  INT4 numDirect;  /* The number of direct filter coefficients. */
  INT4 numRecurs;  /* The number of recursive filter coefficients. */
  INT4 num;        /* An extra counter for error checking. */
  COMPLEX16 *zeros; /* The zeros of the transfer function. */
  COMPLEX16 *poles; /* The poles of the transfer function. */
  REAL8 *direct;   /* The direct filter coefficients. */
  REAL8 *recurs;   /* The recursive filter coefficients. */
  REAL8 *history;  /* The filter history. */

  /* Make sure all the input structures have been initialized. */
  if ( ! input )
    XLAL_ERROR_NULL( XLAL_EFAULT );
  if ( ! input->zeros || ! input->poles
      || ! input->zeros->data || ! input->poles->data )
    XLAL_ERROR_NULL( XLAL_EINVAL );

  numZeros=input->zeros->length;
  numPoles=input->poles->length;
  zeros=input->zeros->data;
  poles=input->poles->data;

  /* Check that zeros are appropriately paired.  Also, keep track of
     the number of zeros at z=0, since these will reduce the number of
     direct coefficients required. */
  numDirect=1;
  for(i=0,num=0;i<numZeros;i++)
    if(cimag(zeros[i])==0.0){
      num+=1;
      if(creal(zeros[i])==0.0)
	numDirect-=1;
    }
    else if(cimag(zeros[i])>0.0){
      num+=2;

#ifndef NDEBUG
      if(lalDebugLevel&LALWARNING){
	/* Check to see that another zero is an actual conjugate.
           This is not a foolproof test, as it can be fooled by
           multiple zeros at the same location. */
	INT4 j=0;
	INT4 k=0;
	REAL8 x = creal(zeros[i]) - creal(zeros[0]);
	REAL8 y = cimag(zeros[i]) + cimag(zeros[0]);
	REAL8 sep = x*x + y*y;
	for(j=1;j<numZeros;j++){
	  x=creal(zeros[i])-creal(zeros[j]);
	  y=cimag(zeros[i])+cimag(zeros[j]);
	  if(sep>x*x+y*y){
	    sep=x*x+y*y;
	    k=j;
	  }
	}
	if(sep>LAL_REAL8_EPS){
          XLALPrintWarning( "XLAL Warning - %s: ", __func__ );
          XLALPrintWarning("Complex zero has no conjugate pair\n");
	  XLALPrintWarning("\tUnmatched zero z_%i = %.8e + i*%.8e\n",i,
              creal(zeros[i]),cimag(zeros[i]));
	  XLALPrintWarning("\tNearest pair   z_%i = %.8e + i*%.8e\n",k,
              creal(zeros[k]),cimag(zeros[k]));
	}
      }
#endif
    }

  if ( num != numZeros )
    XLAL_ERROR_NULL( XLAL_EINVAL );

  /* Check that poles are appropriately paired.  Also, keep track of
     the number of poles at z=0, since these will reduce the number of
     recursive coefficients required. */
  numRecurs=1;
  for(i=0,num=0;i<numPoles;i++)
    if(cimag(poles[i])==0.0){
      num+=1;
      if(creal(poles[i])==0.0)
	numRecurs-=1;
    }
    else if(cimag(poles[i])>0.0){
      num+=2;

#ifndef NDEBUG
      if(lalDebugLevel&LALWARNING){
	/* Check to see that another pole is an actual conjugate.
           This is not a foolproof test, as it can be fooled by
           multiple poles at the same location. */
	INT4 j=0;
	INT4 k=0;
	REAL8 x = creal(poles[i]) - creal(poles[0]);
	REAL8 y = cimag(poles[i]) + cimag(poles[0]);
	REAL8 sep = x*x + y*y;
	for(j=1;j<numPoles;j++){
	  x=creal(poles[i])-creal(poles[j]);
	  y=cimag(poles[i])+cimag(poles[j]);
	  if(sep>x*x+y*y){
	    sep=x*x+y*y;
	    k=j;
	  }
	}
	if(sep>LAL_REAL8_EPS){
          XLALPrintWarning( "XLAL Warning - %s: ", __func__ );
          XLALPrintWarning("Complex pole has no conjugate pair\n");
	  XLALPrintWarning("\tUnmatched pole p_%i = %.8e + i*%.8e\n",i,
              creal(poles[i]),cimag(poles[i]));
	  XLALPrintWarning("\tNearest pair   p_%i = %.8e + i*%.8e\n",k,
              creal(poles[k]),cimag(poles[k]));
	}
      }
#endif
    }

  if ( num != numPoles )
    XLAL_ERROR_NULL( XLAL_EINVAL );

#ifndef NDEBUG
  if(lalDebugLevel&LALWARNING){
    /* Issue a warning if the gain is nonreal. */
    if(fabs(cimag(input->gain))>fabs(LAL_REAL8_EPS*creal(input->gain))){
      XLALPrintWarning( "XLAL Warning - %s: ", __func__ );
      XLALPrintWarning("Gain is non-real\n");
      XLALPrintWarning("\tg = %.8e + i*%.8e\n", creal(input->gain), cimag(input->gain));
    }
    /* Issue a warning if there are any ``removeable'' poles. */
    for(i=0;i<numPoles;i++){
      INT4 j=0;
      for(;j<numZeros;j++)
	if((creal(poles[i])==creal(zeros[j]))&&(cimag(poles[i])==cimag(zeros[j]))){
          XLALPrintWarning( "XLAL Warning - %s: ", __func__ );
	  XLALPrintWarning("Removeable pole\n");
	  XLALPrintWarning("\tp_%i = z_%i = %.8e + i*%.8e\n",i,j,
              creal(poles[i]),cimag(poles[i]));
	}
    }
    /* Issue a warning if extra factors of 1/z will be applied. */
    if(numPoles<numZeros){
      XLALPrintWarning( "XLAL Warning - %s: ", __func__ );
      XLALPrintWarning("Filter has more zeros than poles\n");
      XLALPrintWarning("\t%i poles added at complex origin\n",
          numZeros-numPoles);
    }
    /* Issue a warning if any poles are outside |z|=1. */
    for(i=0;i<numPoles;i++){
      REAL8 zAbs=creal(poles[i])*creal(poles[i])+cimag(poles[i])*cimag(poles[i]);
      if(zAbs>1.0){
        XLALPrintWarning( "XLAL Warning - %s: ", __func__ );
	XLALPrintWarning("Filter has pole outside of unit circle\n");
	XLALPrintWarning("\tp_%i = %.8e + i*%.8e, |p_%i| = %.8e\n",i,
		      creal(poles[i]),cimag(poles[i]),i,zAbs);
      }
    }
  }
#endif

  /* Everything seems okay, so initialize the filter. */
  output=LALCalloc(1,sizeof(*output));
  if ( ! output )
    XLAL_ERROR_NULL( XLAL_ENOMEM );
  num = (numPoles>=numZeros) ? numZeros : numPoles;
  numDirect+=num;
  numRecurs+=num;

  output->deltaT=input->deltaT;
  output->directCoef = XLALCreateREAL8Vector( numDirect );
  output->recursCoef = XLALCreateREAL8Vector( numRecurs );
  if ( ! output->directCoef || ! output->recursCoef )
  {
    XLALDestroyREAL8IIRFilter( output );
    XLAL_ERROR_NULL( XLAL_EFUNC );
  }
  direct=output->directCoef->data;
  recurs=output->recursCoef->data;

  /* Expand the denominator as a polynomial in z. */
  *recurs=-1.0;
  for(i=1;i<numRecurs;i++)
    recurs[i]=0.0;
  for(i=0;i<numPoles;i++){
    INT4 j=numRecurs-1;
    REAL8 x=creal(poles[i]);
    REAL8 y=cimag(poles[i]);
    if(y==0.0)
      for(;j>0;j--)
	recurs[j]-=x*recurs[j-1];
    else if(y>0.0){
      for(;j>1;j--)
	recurs[j]-=2.0*x*recurs[j-1]-(x*x+y*y)*recurs[j-2];
      recurs[j]-=2.0*x*recurs[j-1];
    }
  }

  /* Expand the numerator as a polynomial in z. */
  for(i=0;i<numDirect;i++)
    direct[i]=0.0;
  direct[num-numZeros]=creal(input->gain);
  for(i=0;i<numZeros;i++){
    INT4 j=numDirect-1;
    REAL8 x=creal(zeros[i]);
    REAL8 y=cimag(zeros[i]);
    if(y==0.0)
      for(;j>num-numZeros;j--)
	direct[j]-=x*direct[j-1];
    else if(y>0.0){
      for(;j>num-numZeros+1;j--)
	direct[j]-=2.0*x*direct[j-1]-(x*x+y*y)*direct[j-2];
      direct[j]-=2.0*x*direct[j-1];
    }
  }

  /* Initialize filter history. */
  if(numDirect>=numRecurs)
    num=numDirect-1;
  else
    num=numRecurs-1;
  output->history = XLALCreateREAL8Vector( numRecurs );
  if ( ! output->history )
  {
    XLALDestroyREAL8IIRFilter( output );
    XLAL_ERROR_NULL( XLAL_EFUNC );
  }
  history=output->history->data;
  for(i=0;i<num;i++)
    history[i]=0.0;

  /* Normal exit */
  return output;
}


/** Deprecated.
 * \deprecated Use XLALCreateREAL4IIRFilter() instead */
void LALCreateREAL4IIRFilter( LALStatus         *stat,
			      REAL4IIRFilter    **output,
			      COMPLEX8ZPGFilter *input )
{
  INITSTATUS(stat);

  /* Make sure all the input structures have been initialized. */
  ASSERT(input,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(input->zeros,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(input->zeros->data,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(input->poles,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(input->poles->data,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);

  /* Make sure that the output handle exists, but points to a null
     pointer. */
  ASSERT(output,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(!*output,stat,IIRFILTERH_EOUT,IIRFILTERH_MSGEOUT);

  *output = XLALCreateREAL4IIRFilter( input );
  if ( ! *output )
  {
    int code = xlalErrno & ~XLAL_EFUNC;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EINVAL:
        ABORT(stat,IIRFILTERH_EPAIR,IIRFILTERH_MSGEPAIR);
      case XLAL_ENOMEM:
        ABORT(stat,IIRFILTERH_EMEM,IIRFILTERH_MSGEMEM);
      default:
        ABORTXLAL(stat);
    }
  }

  RETURN(stat);
}


/** Deprecated.
 * \deprecated Use XLALCreateREAL8IIRFilter() instead */
void LALCreateREAL8IIRFilter( LALStatus          *stat,
			      REAL8IIRFilter     **output,
			      COMPLEX16ZPGFilter *input )
{
  INITSTATUS(stat);

  /* Make sure all the input structures have been initialized. */
  ASSERT(input,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(input->zeros,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(input->zeros->data,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(input->poles,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(input->poles->data,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);

  /* Make sure that the output handle exists, but points to a null
     pointer. */
  ASSERT(output,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(!*output,stat,IIRFILTERH_EOUT,IIRFILTERH_MSGEOUT);

  *output = XLALCreateREAL8IIRFilter( input );
  if ( ! *output )
  {
    int code = xlalErrno & ~XLAL_EFUNC;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EINVAL:
        ABORT(stat,IIRFILTERH_EPAIR,IIRFILTERH_MSGEPAIR);
      case XLAL_ENOMEM:
        ABORT(stat,IIRFILTERH_EMEM,IIRFILTERH_MSGEMEM);
      default:
        ABORTXLAL(stat);
    }
  }

  RETURN(stat);
}
/*@}*/

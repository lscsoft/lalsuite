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


#include<lal/LALStdlib.h>
#include<lal/PulsarTimes.h>

/**
 * \file
 * \author Creighton, T. D.
 * \ingroup PulsarTimes_h
 * \brief Computes the rotation-synchronized time coordinate for an object with smoothly-varying spin.
 *
 * ### Description ###
 *
 * These routines compute the value and derivatives of a time
 * transformation \f$t_s(t)\f$ from a physical (inertial) time coordinate to
 * one that is synchronized with the rotation of an object (such as a
 * pulsar).  That is, the rotation period is constant in the new time
 * coordinate, while the actual physical rotation rate may be changing
 * gradually.  The frequency drift \f$f(t)\f$ is characterized by
 * frequency-normalized Taylor coefficients \f$f_k=(k!f)^{-1}\partial^k f/\partial t^k\f$.
 *
 * The routines obey the calling convention presented in the module
 * \ref PulsarTimes_h, with \f$n\f$ variable parameters \f$\lambda^k=f_k\f$
 * (measured in \f$\mathrm{Hz}^k\f$), where \f$k=1,\ldots,n\f$.  The only
 * constant parameter field used by these routines is
 * <tt>constants->t0</tt>, which is the time when the Taylor coefficients
 * \f$f_k\f$ are evaluated.  \f$t_s\f$ is defined to be zero at that time.
 *
 * <tt>*variables</tt> can be a vector of arbitrary length; for
 * consistency, <tt>*dtSpin</tt> must be one element longer.
 *
 * ### Algorithm ###
 *
 * If the frequency of a rotating body is varying in some suitably smooth
 * manner, then it can be represented as a Taylor series:
 * \f[
 * f(t) = f_0 \left(1 + \sum_{k=1}^n f_k [t-t_0]^k \right) \; .
 * \f]
 * It is worth pointing out that the parameters \f$f_k\f$ are
 * frequency-normalized Taylor coefficients, so their magnitudes
 * represent, not the frequency itself, but the timescale on which the
 * frequency is varying.  That is, if the frequency is changing on some
 * timescale \f$\mathcal{T}\f$, then \f$f_k\f$ is typically of order
 * \f$\mathcal{T}^{-k}\f$.
 *
 * The spin-synchronized time \f$t_s(t)\f$ counts off intervals of constant
 * rotational phase, so it must be proportional to \f$\int f(t)\,dt\f$.  The
 * integration constant is set by the convention that \f$t_s(t_0)=0\f$.  The
 * constant of proportionality is set by demanding that, at \f$t=t_0\f$, both
 * \f$t_s\f$ and \f$t\f$ measure time at the same rate.  This gives us the
 * transformation:
 * \f[
 * t_s(t) = t-t_0 + \sum_{k=1}^n \frac{f_k}{k+1} (t-t_0)^{k+1} \; .
 * \f]
 * The time derivative of this transformation is of course just the the
 * factor in the equation for \f$f(t)\f$ that we started with:
 * \f[
 * \frac{\partial t_s(t)}{\partial t} = 1 + \sum_{k=1}^n f_k [t-t_0]^k \; .
 * \f]
 * The derivatives with respect to \f$f_k\f$ are similarly trivial:
 * \f[
 * \frac{\partial t_s(t)}{\partial f_k} = \frac{(t-t_0)^{k+1}}{k+1} \; .
 * \f]
 *
 * ### Uses ###
 *
 * \code
 * lalDebugLevel
 * \endcode
 *
 */
/*@{*/

void
LALTSpin( LALStatus             *stat,
	  REAL8                 *tSpin,
	  REAL8Vector           *variables,
	  PulsarTimesParamStruc *constants )
{
  INT4 i;      /* An index. */
  INT4 k;      /* Another index. */
  REAL8 t;     /* A time variable, t-t0. */
  REAL8 tk;    /* t raised to some power. */
  REAL8 ts;    /* Another time variable storing ts(t). */
  REAL8 *data; /* Pointer to variables->data. */

  INITSTATUS(stat);

  /* This function may be called a lot.  Do error checking only in
     debug mode. */
#ifndef NDEBUG
  if(lalDebugLevel){

    /* Make sure parameter structures and their fields exist. */
    ASSERT(tSpin,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
    ASSERT(variables,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
    ASSERT(variables->data,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
    ASSERT(constants,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
  }
#endif

  /* Set some temporary variables. */
  data=variables->data;
  ts=t=tk=*(data++)-constants->t0;

  /* We'll compute ts and all its derivatives at once, so as to avoid
     repeated calculation of the powers of t. */
  i=variables->length-1;
  k=1;
  while(--i)
    ts+=*(data++)*(tk*=t)/(++k);

  /* Assign the output, and exit. */
  *tSpin=ts;
  RETURN(stat);
}


void
LALDTSpin( LALStatus             *stat,
	   REAL8Vector           *dtSpin,
	   REAL8Vector           *variables,
	   PulsarTimesParamStruc *constants )
{
  INT4 i;       /* An index. */
  INT4 k;       /* Another index. */
  REAL8 t;      /* A time variable, t-t0. */
  REAL8 tk;     /* t raised to some power. */
  REAL8 ts;     /* Another time variable storing ts(t). */
  REAL8 dts;    /* A final variable storing dts/dt. */
  REAL8 *data1; /* Pointer to variables->data. */
  REAL8 *data2; /* Pointer to dtSpin->data. */

  INITSTATUS(stat);

  /* This function may be called a lot.  Do error checking only in
     debug mode. */
#ifndef NDEBUG
  if(lalDebugLevel){

    /* Make sure parameter structures and their fields exist. */
    ASSERT(dtSpin,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
    ASSERT(dtSpin->data,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
    ASSERT(variables,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
    ASSERT(variables->data,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
    ASSERT(constants,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);

    /* Make sure array sizes are consistent. */
    ASSERT(dtSpin->length==variables->length+1,stat,
	   PULSARTIMESH_EBAD,PULSARTIMESH_MSGEBAD);
  }
#endif

  /* Set some temporary variables. */
  data1=variables->data;
  data2=dtSpin->data+2;
  ts=t=tk=*(data1++)-constants->t0;
  dts=1.0;

  /* We'll compute ts and all its derivatives at once, so as to avoid
     repeated calculation of the powers of t. */
  i=variables->length;
  k=1;
  while(--i){
    REAL8 fk=*(data1++);
    dts+=fk*tk;
    ts+=fk*(*(data2++)=(tk*=t)/(++k));
    /* I just love this instruction!  Here it is in expanded form:

       tk = tk*t;
       k = k + 1;
       *data2 = tk/k;
       ts = ts + fk*tk/k;
       data2 = data2 + 1;

       The reason for the obfuscated C code is to eliminate all
       redundant computations or variable access. */
  }

  /* Assign the first two output fields, and exit. */
  data2=dtSpin->data;
  *(data2++)=ts;
  *data2=dts;
  RETURN(stat);
}
/*@}*/

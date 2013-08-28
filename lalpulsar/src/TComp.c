/*
*  Copyright (C) 2007 Jolien Creighton
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
 * \brief Computes the composition of two time transformations.
 *
 * ### Description ###
 *
 * These routines compute the value and derivatives of a time
 * transformation \f$t_c(t)\f$ that is the composition of two other
 * transformations \f$t_1\f$ and \f$t_2\f$; that is, \f$t_c(t)=t_2(t_1(t))\f$.  More
 * precisely, the transformation is
 * \f$t_c(t,\vec\lambda_{(1)}\oplus\vec\lambda_{(2)}) =
 * t_2[t_1(t,\vec\lambda_{(1)}),\vec\lambda_{(2)}]\f$.  Note that
 * \f$\vec\lambda_{(1)}\f$ and \f$\vec\lambda_{(2)}\f$ are assumed to represent
 * \e independent sets of parameters.  If there is any overlap
 * between the parameter sets, LALDTComp() will \e not correctly
 * compute the derivatives of \f$t_c(t)\f$ (although the \e value of
 * \f$t_c\f$ will still be correct).
 *
 * The routines obey the calling convention presented in
 * \ref PulsarTimes_h.  The contents of <tt>*variables</tt> are, firstly,
 * the time \f$t\f$ that will be sent to \f$t_1(t,\vec\lambda_{(1)})\f$; next,
 * the \f$n\f$ parameters \f$\lambda^1,\ldots,\lambda^n\f$ that will be sent to
 * \f$t_1(t,\vec\lambda_{(1)})\f$ as
 * \f$\lambda_{(1)}^1,\ldots,\lambda_{(1)}^n\f$; last, the \f$m\f$ parameters
 * \f$\lambda^{n+1},\ldots,\lambda^{n+m}\f$ that will be sent to
 * \f$t_2(t_1,\vec\lambda_{(2)})\f$ as
 * \f$\lambda_{(2)}^1,\ldots,\lambda_{(2)}^m\f$.  Here \f$n\f$ and \f$m\f$ are the
 * number of variable parameters expected by the transformations \f$t_1\f$
 * and \f$t_2\f$, so that <tt>variables->length</tt>\f$=n+m+1\f$.
 *
 * The constant parameter fields used by these routines are:
 * <dl>
 * <dt><tt>constants->t1</tt></dt><dd> A function pointer to the function that evaluates \f$t_1(t)\f$.</dd>
 * <dt><tt>constants->t2</tt></dt><dd> A function pointer to the function that evaluates \f$t_2(t)\f$.</dd>
 * <dt><tt>constants->dt1</tt></dt><dd> A function pointer to the function that evaluates \f$t_1(t)\f$ \e and its derivatives.</dd>
 * <dt><tt>constants->dt2</tt></dt><dd> A function pointer to the function that evaluates \f$t_2(t)\f$ \e and its derivatives.</dd>
 * <dt><tt>constants->constants1</tt></dt><dd> A pointer to the constant parameters used by <tt>constants->t1</tt> and <tt>constants->dt1</tt>.</dd>
 * <dt><tt>constants->constants2</tt></dt><dd> A pointer to the constant parameters used by <tt>constants->t2</tt> and <tt>constants->dt2</tt>.</dd>
 * <dt><tt>constants->nArgs</tt></dt><dd> The number \f$n\f$ of variable parameters to be sent to \f$t_1(t)\f$.</dd>
 * </dl>
 *
 * Note that the number of variable parameters to be sent to \f$t_2(t)\f$ is
 * not specified in \c constants; after sending the first
 * <tt>constants->nArgs</tt> of them to \f$t_1(t)\f$, the remaining parameters
 * (however many they are) are sent to \f$t_2(t)\f$.  This is particularly
 * useful for pulsar timing routines, where the last function in the
 * composition chain is often a transformation that corrects for the
 * pulsar spindown, using an arbitrary number of spindown parameters.
 * The number of spindown parameters desired is then specified
 * unambiguously by setting <tt>variables->length</tt>.  Note however that
 * <tt>*dtComp</tt> must always have a length exactly one greater than
 * <tt>*variables</tt>.
 *
 * ### Algorithm ###
 *
 * Computing the value of \f$t_c\f$ is trivial:
 * \f[
 * t_c(t) = t_2(t_1(t)) \; .
 * \f]
 * The only trickiness is in handling the parameters, which is done using
 * a local REAL8Vector \a variables.  This vector is not given its
 * own memory; instead, its \c data field is pointed at either the
 * first or the last block of parameters in <tt>variable->data</tt>.
 *
 * Computing the derivatives of \f$t_c\f$ is not much trickier.  For the time
 * variable and the first \f$n\f$ parameters, the chain rule gives us:
 * \f{eqnarray}{
 * \frac{\partial t_c(t)}{\partial t} & = &
 * \frac{\partial t_2(t_1)}{\partial t_1}
 * \frac{\partial t_1(t)}{\partial t} \; , \nonumber\\
 * \frac{\partial t_c(t)}{\partial\lambda^i} & = &
 * \frac{\partial t_2(t_1)}{\partial t_1}
 * \frac{\partial t_1(t)}{\partial\lambda_{(1)}^i} \; ,
 * \quad i=1,\ldots,n \; . \nonumber
 * \f}
 * For the remaining parameters,
 * \f[
 * \frac{\partial t_c(t)}{\partial\lambda^j} =
 * \frac{\partial t_2(t_1)}{\partial\lambda_{(2)}^{j-n}} \; ,
 * \quad j=n+1,\ldots \; .
 * \f]
 *
 * As noted in the module description, the derivatives will not be
 * evaluated correctly if there is overlap between the two parameter sets
 * \f$\vec\lambda_{(1)}\f$ and \f$\vec\lambda_{(2)}\f$.  In particular, if some
 * variable \f$\alpha\f$ is represented both by \f$\lambda^i=\lambda_{(1)}^i\f$
 * and \f$\lambda^{n+j}=\lambda_{(2)}^j\f$, the value of \f$\partial
 * t_c/\partial\alpha\f$ is neither given by \f$\partial
 * t_c/\partial\lambda^i\f$ nor by \f$\partial t_c/\partial\lambda^{n+j}\f$,
 * but by:
 * \f[
 * \frac{\partial t_c}{\partial\alpha} =
 * \frac{\partial t_c}{\partial\lambda^{n+j}} +
 * \frac{\partial t_c}{\partial\lambda^{i}}
 * \frac{\partial t_1}{\partial\lambda_{(1)}^i} \; .
 * \f]
 * While this is not especially difficult to evaluate, it is impossible
 * to code without giving LALDTComp() some way of knowing
 * \e which parameters in \f$\vec\lambda_{(1)}\f$ and \f$\vec\lambda_{(2)}\f$
 * represent the same physical quantity.  Such a scheme is not
 * implemented at present.
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
LALTComp( LALStatus             *stat,
	  REAL8                 *tComp,
	  REAL8Vector           *variables,
	  PulsarTimesParamStruc *constants )
{
  INT4 n;     /* Number of variables to be sent to t1(). */
  INT4 m;     /* Number of variables to be sent to t2(). */
  REAL8 temp; /* Temporary storage variable. */
  REAL8 t1;   /* Value of t_1(t). */
  REAL8Vector variablesIn; /* Variables to be sent to t1,2(). */

  INITSTATUS(stat);
  ATTATCHSTATUSPTR(stat);

  /* This function may be called a lot.  Do error checking only in
     debug mode. */
#ifndef NDEBUG
  if(lalDebugLevel){

    /* Make sure parameter structures and their fields exist. */
    ASSERT(tComp,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
    ASSERT(variables,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
    ASSERT(variables->data,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
    ASSERT(constants,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
    ASSERT(constants->t1,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
    ASSERT(constants->t2,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);

    /* Make sure array sizes are consistent. */
    ASSERT(variables->length>constants->nArgs,stat,
	   PULSARTIMESH_EBAD,PULSARTIMESH_MSGEBAD);
  }
#endif

  /* Set up the vectors to pass into constants->t1(). */
  n=constants->nArgs;
  variablesIn.length=n+1;
  variablesIn.data=variables->data;

  /* Compute t_1 and its derivatives. */
  (constants->t1)(stat->statusPtr,&t1,&variablesIn,
		  constants->constants1);
  CHECKSTATUSPTR(stat);

  /* Set up the vectors to pass into constants->t2().  We set
     variablesIn so that we will read from the correct block
     *variables.  Unfortunately, this will require us to put t1 into
     variablesIn->data[0] with the value t1, so we should save the
     data previously there. */
  m=variables->length-n-1;
  variablesIn.length=m+1;
  temp=*(variablesIn.data=variables->data+n);
  *(variablesIn.data)=t1;

  /* Compute t_2 and its derivatives.  Restore the *variables before
     checking the status pointer; otherwise, an input field of TComp()
     could be modified on return, which is a violation of LAL
     standards. */
  (constants->t2)(stat->statusPtr,tComp,&variablesIn,
		  constants->constants2);
  *(variablesIn.data)=temp;
  CHECKSTATUSPTR(stat);

  /* The value of t_2 should already be in its proper place, so now
     we're done. */
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}


void
LALDTComp( LALStatus             *stat,
	   REAL8Vector           *dtComp,
	   REAL8Vector           *variables,
	   PulsarTimesParamStruc *constants )
{
  INT4 n;       /* Number of variables to be sent to dt1(). */
  INT4 m;       /* Number of variables to be sent to dt2(). */
  REAL8 temp1;  /* Temporary storage variable. */
  REAL8 temp2;  /* Another temporary storage variable. */
  REAL8 temp3;  /* Yet another temporary storage variable. */
  REAL8 dt2dt1; /* Derivative of t_2(t_1) with respect to t_1. */
  REAL8 *data;  /* A multipurpose pointer to vector data. */
  REAL8Vector variablesIn; /* Variables to be sent to dt1,2(). */
  REAL8Vector dtOut;       /* Derivatives returned by dt1,2(). */

  INITSTATUS(stat);
  ATTATCHSTATUSPTR(stat);

  /* This function may be called a lot.  Do error checking only in
     debug mode. */
#ifndef NDEBUG
  if(lalDebugLevel){

    /* Make sure parameter structures and their fields exist. */
    ASSERT(dtComp,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
    ASSERT(dtComp->data,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
    ASSERT(variables,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
    ASSERT(variables->data,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
    ASSERT(constants,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
    ASSERT(constants->dt1,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
    ASSERT(constants->dt2,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);

    /* Make sure array sizes are consistent. */
    ASSERT(dtComp->length==variables->length+1,stat,
	   PULSARTIMESH_EBAD,PULSARTIMESH_MSGEBAD);
    ASSERT(variables->length>constants->nArgs,stat,
	   PULSARTIMESH_EBAD,PULSARTIMESH_MSGEBAD);
  }
#endif

  /* Set up the vectors to pass into constants->dt1().  The output
     will be stored temporarily in the space allocated to *dtComp. */
  n=constants->nArgs;
  dtOut.length=n+2;
  dtOut.data=dtComp->data;
  variablesIn.length=n+1;
  variablesIn.data=variables->data;

  /* Compute t_1 and its derivatives. */
  (constants->dt1)(stat->statusPtr,&dtOut,&variablesIn,
		   constants->constants1);
  CHECKSTATUSPTR(stat);

  /* Set up the vectors to pass into constants->dt2().  We set up
     dtOut and variablesIn so that we will write to and read from the
     correct blocks of *dtComp and *variables, respectively.
     Unfortunately, this will involve overwriting the data at
     dtOut->data[0] and variablesIn->data[0], so we should save
     these. */
  m=variables->length-n-1;
  dtOut.length=m+2;
  temp1=*(dtOut.data=dtComp->data+n);
  temp2=*(dtOut.data+1);
  variablesIn.length=m+1;
  temp3=*(variablesIn.data=variables->data+n);
  *(variablesIn.data)=*(dtComp->data);

  /* Compute t_2 and its derivatives.  Unmangle the arrays before
     checking the status pointer; otherwise, the input field
     *variables of DTComp() could have been modified on return, which
     is a violation of LAL standards. */
  (constants->dt2)(stat->statusPtr,&dtOut,&variablesIn,
		   constants->constants2);
  *(dtComp->data)=*(dtOut.data);
  dt2dt1=*(dtOut.data+1);
  *(dtOut.data)=temp1;
  *(dtOut.data+1)=temp2;
  *(variablesIn.data)=temp3;
  CHECKSTATUSPTR(stat);

  /* Apply the chain rule to the derivatives of t_1. */
  data=dtComp->data+1;
  n++;
  while(n--)
    *(data++)*=dt2dt1;

  /* Everything else should already be in its proper place, so now
     we're done. */
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}
/*@}*/

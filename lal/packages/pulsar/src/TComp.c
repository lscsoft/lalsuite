/**************************************** <lalVerbatim file="TCompCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{TComp.c}}
\label{ss:TComp.c}

Computes the composition of two time transformations.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{TCompCP}
\idx{LALTComp()}
\idx{LALDTComp()}

\subsubsection*{Description}

These routines compute the value and derivatives of a time
transformation $t_c(t)$ that is the composition of two other
transformations $t_1$ and $t_2$; that is, $t_c(t)=t_2(t_1(t))$.  More
precisely, the transformation is
$t_c(t,\vec\lambda_{(1)}\oplus\vec\lambda_{(2)}) =
t_2[t_1(t,\vec\lambda_{(1)}),\vec\lambda_{(2)}]$.  Note that
$\vec\lambda_{(1)}$ and $\vec\lambda_{(2)}$ are assumed to represent
\emph{independent} sets of parameters.  If there is any overlap
between the parameter sets, \verb@DTComp()@ will \emph{not} correctly
compute the derivatives of $t_c(t)$ (although the \emph{value} of
$t_c$ will still be correct).

The routines obey the calling convention presented in the header
\verb@PulsarTimes.h@.  The contents of \verb@*variables@ are, firstly,
the time $t$ that will be sent to $t_1(t,\vec\lambda_{(1)})$; next,
the $n$ parameters $\lambda^1,\ldots,\lambda^n$ that will be sent to
$t_1(t,\vec\lambda_{(1)})$ as
$\lambda_{(1)}^1,\ldots,\lambda_{(1)}^n$; last, the $m$ parameters
$\lambda^{n+1},\ldots,\lambda^{n+m}$ that will be sent to
$t_2(t_1,\vec\lambda_{(2)})$ as
$\lambda_{(2)}^1,\ldots,\lambda_{(2)}^m$.  Here $n$ and $m$ are the
number of variable parameters expected by the transformations $t_1$
and $t_2$, so that \verb@variables->length@$=n+m+1$.

The constant parameter fields used by these routines are:
\begin{description}
\item[\texttt{constants->t1}] A function pointer to the function that
evaluates $t_1(t)$.

\item[\texttt{constants->t2}] A function pointer to the function that
evaluates $t_2(t)$.

\item[\texttt{constants->dt1}] A function pointer to the function that
evaluates $t_1(t)$ \emph{and} its derivatives.

\item[\texttt{constants->dt2}] A function pointer to the function that
evaluates $t_2(t)$ \emph{and} its derivatives.

\item[\texttt{constants->constants1}] A pointer to the constant
parameters used by \verb@constants->t1@ and \verb@constants->dt1@.

\item[\texttt{constants->constants2}] A pointer to the constant
parameters used by \verb@constants->t2@ and \verb@constants->dt2@.

\item[\texttt{constants->nArgs}] The number $n$ of variable parameters to be sent to $t_1(t)$.
\end{description}

Note that the number of variable parameters to be sent to $t_2(t)$ is
not specified in \verb@constants@; after sending the first
\verb@constants->nArgs@ of them to $t_1(t)$, the remaining parameters
(however many they are) are sent to $t_2(t)$.  This is particularly
useful for pulsar timing routines, where the last function in the
composition chain is often a transformation that corrects for the
pulsar spindown, using an arbitrary number of spindown parameters.
The number of spindown parameters desired is then specified
unambiguously by setting \verb@variables->length@.  Note however that
\verb@*dtComp@ must always have a length exactly one greater than
\verb@*variables@.

\subsubsection*{Algorithm}

Computing the value of $t_c$ is trivial:
$$
t_c(t) = t_2(t_1(t)) \; .
$$
The only trickiness is in handling the parameters, which is done using
a local \verb@REAL8Vector variablesIn@.  This vector is not given its
own memory; instead, its \verb@data@ field is pointed at either the
first or the last block of parameters in \verb@variable->data@.

Computing the derivatives of $t_c$ is not much trickier.  For the time
variable and the first $n$ parameters, the chain rule gives us:
\begin{eqnarray}
\frac{\partial t_c(t)}{\partial t} & = &
	\frac{\partial t_2(t_1)}{\partial t_1}
	\frac{\partial t_1(t)}{\partial t} \; , \nonumber\\
\frac{\partial t_c(t)}{\partial\lambda^i} & = &
	\frac{\partial t_2(t_1)}{\partial t_1}
	\frac{\partial t_1(t)}{\partial\lambda_{(1)}^i} \; ,
		\quad i=1,\ldots,n \; . \nonumber
\end{eqnarray}
For the remaining parameters,
$$
\frac{\partial t_c(t)}{\partial\lambda^j} =
	\frac{\partial t_2(t_1)}{\partial\lambda_{(2)}^{j-n}} \; ,
		\quad j=n+1,\ldots \; .
$$

As noted in the module description, the derivatives will not be
evaluated correctly if there is overlap between the two parameter sets
$\vec\lambda_{(1)}$ and $\vec\lambda_{(2)}$.  In particular, if some
variable $\alpha$ is represented both by $\lambda^i=\lambda_{(1)}^i$
and $\lambda^{n+j}=\lambda_{(2)}^j$, the value of $\partial
t_c/\partial\alpha$ is neither given by $\partial
t_c/\partial\lambda^i$ nor by $\partial t_c/\partial\lambda^{n+j}$,
but by:
$$
\frac{\partial t_c}{\partial\alpha} =
	\frac{\partial t_c}{\partial\lambda^{n+j}} +
	\frac{\partial t_c}{\partial\lambda^{i}}
		\frac{\partial t_1}{\partial\lambda_{(1)}^i} \; .
$$
While this is not especially difficult to evaluate, it is impossible
to code without giving \verb@DTComp()@ some way of knowing
\emph{which} parameters in $\vec\lambda_{(1)}$ and $\vec\lambda_{(2)}$
represent the same physical quantity.  Such a scheme is not
implemented at present.

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{TCompCV}}

******************************************************* </lalLaTeX> */

#include<lal/LALStdlib.h>
#include<lal/PulsarTimes.h>

NRCSID(TCOMPC,"$Id$");

/* <lalVerbatim file="TCompCP"> */
void
LALTComp( LALStatus             *stat,
	  REAL8                 *tComp,
	  REAL8Vector           *variables,
	  PulsarTimesParamStruc *constants )
{ /* </lalVerbatim> */
  INT4 n;     /* Number of variables to be sent to t1(). */
  INT4 m;     /* Number of variables to be sent to t2(). */
  REAL8 temp; /* Temporary storage variable. */
  REAL8 t1;   /* Value of t_1(t). */
  REAL8Vector variablesIn; /* Variables to be sent to t1,2(). */

  INITSTATUS(stat,"TComp",TCOMPC);
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


/* <lalVerbatim file="TCompCP"> */
void
LALDTComp( LALStatus             *stat,
	   REAL8Vector           *dtComp,
	   REAL8Vector           *variables,
	   PulsarTimesParamStruc *constants )
{ /* </lalVerbatim> */
  INT4 n;       /* Number of variables to be sent to dt1(). */
  INT4 m;       /* Number of variables to be sent to dt2(). */
  REAL8 temp1;  /* Temporary storage variable. */
  REAL8 temp2;  /* Another temporary storage variable. */
  REAL8 temp3;  /* Yet another temporary storage variable. */
  REAL8 dt2dt1; /* Derivative of t_2(t_1) with respect to t_1. */
  REAL8 *data;  /* A multipurpose pointer to vector data. */
  REAL8Vector variablesIn; /* Variables to be sent to dt1,2(). */
  REAL8Vector dtOut;       /* Derivatives returned by dt1,2(). */

  INITSTATUS(stat,"DTComp",TCOMPC);
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

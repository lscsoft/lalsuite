/**************************************** <lalVerbatim file="TSpinCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{TSpin.c}}
\label{ss:TSpin.c}

Computes the rotation-synchronized time coordinate for an object with
smoothly-varying spin.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{TSpinCP}
\index{\texttt{LALTSpin()}}
\index{\texttt{LALDTSpin()}}

\subsubsection*{Description}

These routines compute the value and derivatives of a time
transformation $t_s(t)$ from a physical (inertial) time coordinate to
one that is synchronized with the rotation of an object (such as a
pulsar).  That is, the rotation period is constant in the new time
coordinate, while the actual physical rotation rate may be changing
gradually.  The frequency drift $f(t)$ is characterized by
frequency-normalized Taylor coefficients $f_k=(k!f)^{-1}\partial^k
f/\partial t^k$.

The routines obey the calling convention presented in the header
\verb@PulsarTimes.h@, with $n$ variable parameters $\lambda^k=f_k$
(measured in $\mathrm{Hz}^k$), where $k=1,\ldots,n$.  The only
constant parameter field used by these routines is
\verb@constants->t0@, which is the time when the Taylor coefficients
$f_k$ are evaluated.  $t_s$ is defined to be zero at that time.

\verb@*variables@ can be a vector of arbitrary length; for
consistency, \verb@*dtSpin@ must be one element longer.

\subsubsection*{Algorithm}

If the frequency of a rotating body is varying in some suitably smooth
manner, then it can be represented as a Taylor series:
$$
f(t) = f_0 \left(1 + \sum_{k=1}^n f_k [t-t_0]^k \right) \; .
$$
It is worth pointing out that the parameters $f_k$ are
frequency-normalized Taylor coefficients, so their magnitudes
represent, not the frequency itself, but the timescale on which the
frequency is varying.  That is, if the frequency is changing on some
timescale $\mathcal{T}$, then $f_k$ is typically of order
$\mathcal{T}^{-k}$.

The spin-synchronized time $t_s(t)$ counts off intervals of constant
rotational phase, so it must be proportional to $\int f(t)\,dt$.  The
integration constant is set by the convention that $t_s(t_0)=0$.  The
constant of proportionality is set by demanding that, at $t=t_0$, both
$t_s$ and $t$ measure time at the same rate.  This gives us the
transformation:
$$
t_s(t) = t-t_0 + \sum_{k=1}^n \frac{f_k}{k+1} (t-t_0)^{k+1} \; .
$$
The time derivative of this transformation is of course just the the
factor in the equation for $f(t)$ that we started with:
$$
\frac{\partial t_s(t)}{\partial t} = 1 +
	\sum_{k=1}^n f_k [t-t_0]^k \; .
$$
The derivatives with respect to $f_k$ are similarly trivial:
$$
\frac{\partial t_s(t)}{\partial f_k} = \frac{(t-t_0)^{k+1}}{k+1} \; .
$$

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{TSpinCV}}

******************************************************* </lalLaTeX> */

#include<lal/LALStdlib.h>
#include<lal/PulsarTimes.h>

NRCSID(TSPINC,"$Id$");

/* <lalVerbatim file="TSpinCP"> */
void
LALTSpin( LALStatus             *stat,
	  REAL8                 *tSpin,
	  REAL8Vector           *variables,
	  PulsarTimesParamStruc *constants )
{ /* </lalVerbatim> */
  INT4 i;      /* An index. */
  INT4 k;      /* Another index. */
  REAL8 t;     /* A time variable, t-t0. */
  REAL8 tk;    /* t raised to some power. */
  REAL8 ts;    /* Another time variable storing ts(t). */
  REAL8 *data; /* Pointer to variables->data. */

  INITSTATUS(stat,"TSpin",TSPINC);

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


/* <lalVerbatim file="TSpinCP"> */
void
LALDTSpin( LALStatus             *stat,
	   REAL8Vector           *dtSpin,
	   REAL8Vector           *variables,
	   PulsarTimesParamStruc *constants )
{ /* </lalVerbatim> */
  INT4 i;       /* An index. */
  INT4 k;       /* Another index. */
  REAL8 t;      /* A time variable, t-t0. */
  REAL8 tk;     /* t raised to some power. */
  REAL8 ts;     /* Another time variable storing ts(t). */
  REAL8 dts;    /* A final variable storing dts/dt. */
  REAL8 *data1; /* Pointer to variables->data. */
  REAL8 *data2; /* Pointer to dtSpin->data. */

  INITSTATUS(stat,"DTSpin",TSPINC);

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

/**** <lalVerbatim file="ODEHV">
 * Author: J. D. E. Creighton
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \section{Header \texttt{ODE.h}}
 *
 * Routines for solving ordinary differential equations (ODEs).
 *
 * \subsection*{Synopsis}
 * \begin{verbatim}
 * #include <lal/ODE.h>
 * \end{verbatim}
 *
 * These routines solve ordinary differential equations (ODEs) of the form:
 * \begin{displaymath}
 *   \dot{\mathbf{x}} = {\mathbf{f}}({\mathbf{x}},t,\ldots)
 * \end{displaymath}
 * where $\mathbf{f}$ is a specified vector-valued function.
 *
 **** </lalLaTeX> */

#ifndef _ODE_H
#define _ODE_H

#include <lal/LALDatatypes.h>

#ifdef __cplusplus
extern "C" {
#pragma }
#endif

NRCSID( ODEH, "$Id$" );

/**** <lalLaTeX>
 *
 * \subsection*{Error conditions}
 *
 **** </lalLaTeX> */
/**** <lalErrTable> */
#define ODEH_ENULL 001
#define ODEH_ESAME 002
#define ODEH_ESIZE 004
#define ODEH_ESZMM 010
#define ODEH_ENSTP 020
#define ODEH_MSGENULL "Null pointer."
#define ODEH_MSGESAME "Same data pointer."
#define ODEH_MSGESIZE "Invalid size."
#define ODEH_MSGESZMM "Size mismatch."
#define ODEH_MSGENSTP "Step number mismatch."
/**** </lalErrTable> */

/**** <lalLaTeX>
 *
 * \subsection*{Structures}
 * \idx[Type]{REAL4ODEIndep}
 * \idx[Type]{REAL4ODEParams}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct
tagREAL4ODEIndep
{
  REAL4  t;
  void  *aux;
}
REAL4ODEIndep;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 *
 * The independent variables of the ODE (parameters to the ODE function).
 * The fields are:
 * \begin{description}
 * \item[\texttt{t}] The independent parameter (e.g., time) that is evolved.
 * \item[\texttt{aux}] Storage for auxiliary variables used internally in the
 *   ODE routine.
 * \end{description}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct
tagREAL4ODEParams
{
  void ( *ode )( LALStatus *, REAL4Vector *, REAL4Vector *, REAL4ODEIndep * );
  REAL4ODEIndep       *indep;
  REAL4                tstep;
  REAL4Vector         *xdot;
  REAL4Vector         *xerr;
  REAL4                eps;
  REAL4VectorSequence *dx;
}
REAL4ODEParams;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 *
 * The parameters for the ODE step integrator.  The fields are:
 * \begin{description}
 * \item[\texttt{ode}] Pointer to the function that computes the RHS of the ODE.
 * \item[\texttt{indep}] The independent variables of used by this function.
 * \item[\texttt{tstep}] The suggested time step to use.
 * \item[\texttt{xdot}] The value of the LHS of the ODE at the current time.
 * \item[\texttt{xerr}] The estimated errors from the last ODE step.
 * \item[\texttt{eps}] The allowed fractional error per step.  If \verb+eps+ is
 *   less than \verb+LAL_REAL4_EPS+, then the latter is used instead.
 * \item[\texttt{dx}] Workspace storage for use in the the step integrator.
 *   The length of the sequence depends on which integrator is used.
 * \end{description}
 *
 * \vfill{\footnotesize\input{ODEHV}}
 * \newpage\input{ODEC}
 * \newpage\input{ODETestC}
 *
 **** </lalLaTeX> */

void LALSRungeKutta4(
    LALStatus      *status,
    REAL4Vector    *output,
    REAL4Vector    *input,
    REAL4ODEParams *params
    );

void LALSRungeKutta5(
    LALStatus      *status,
    REAL4Vector    *output,
    REAL4Vector    *input,
    REAL4ODEParams *params
    );

void LALSRungeKutta5Adapt(
    LALStatus      *status,
    REAL4Vector    *output,
    REAL4Vector    *input,
    REAL4ODEParams *params
    );

#ifdef __cplusplus
#pragma {
}
#endif

#endif /* _ODE_H */

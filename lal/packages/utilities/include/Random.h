/**** <lalVerbatim file="RandomHV">
 * Author: Creighton, J. D. E. and Tibbits, M. M.
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \section{Header \texttt{Random.h}}
 * \label{s:Random.h}
 * 
 * Generates random numbers.
 * 
 * \subsection*{Synopsis}
 * \begin{verbatim}
 * #include <lal/Random.h>
 * \end{verbatim}
 * 
 * \noindent This header covers the routines for generating random numbers.
 * 
 **** </lalLaTeX> */

#ifndef _RANDOM_H
#define _RANDOM_H

#include <lal/LALDatatypes.h>

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID (RANDOMH, "$Id$");


/**** <lalLaTeX>
 *
 * \subsection*{Error conditions}
 *
 **** </lalLaTeX> */

/* <lalErrTable> */
#define RANDOMH_ENULL 1
#define RANDOMH_ENNUL 2
#define RANDOMH_ESIZE 4
#define RANDOMH_ELNTH 8
#define RANDOMH_ESEGZ 16
#define RANDOMH_ENUMZ 32
#define RANDOMH_EALOC 64
#define RANDOMH_EINIT 128 
#define RANDOMH_EZERO 256
#define RANDOMH_ESEED 512
#define RANDOMH_MSGENULL "Null pointer"
#define RANDOMH_MSGENNUL "Non-null pointer"
#define RANDOMH_MSGESIZE "Invalid size"
#define RANDOMH_MSGELNTH "Must have more than one data point"
#define RANDOMH_MSGESEGZ "Invalid number of segments"
#define RANDOMH_MSGENUMZ "Invalid number of points in segment"
#define RANDOMH_MSGEALOC "Memory Allocation Error"
#define RANDOMH_MSGEINIT "Params must be initialized with CreateParams first"
#define RANDOMH_MSGEZERO "Output Vector length must be greater than zero"
#define RANDOMH_MSGESEED "Improper seed value"
/* </lalErrTable> */

/**** <lalLaTeX>
 * 
 * \subsection*{Structures}
 * \idx[Type]{RandomParams}
 * \idx[Type]{MTRandomParams}
 * 
 * \begin{verbatim}
 * typedef struct tagRandomParams RandomParams;
 * \end{verbatim}
 * 
 * This structure contains the parameters necessary for generating the current
 * sequence of random numbers (based on the initial seed).  The contents should
 * not be manually adjusted.
 *
 * \begin{verbatim}
 * typedef struct tagMTRandomParams MTRandomParams;
 * \end{verbatim}
 * 
 * This structure contains the parameters necessary for generating the current
 * sequence of Mersenne twiser random numbers (based on the initial seed).  The
 * contents should not be manually adjusted.
 * 
 **** </lalLaTeX> */

typedef struct
tagRandomParams
{
  INT4 i;
  INT4 y;
  INT4 v[32];
}
RandomParams;

typedef struct tagMTRandomParams MTRandomParams;

/**** <lalLaTeX>
 * 
 * \newpage\input{RandomC}
 * \newpage\input{RandomTestC}
 * \newpage\input{MersenneRandomC}
 * \newpage\input{MersenneRandomTestC}
 *
 **** </lalLaTeX> */

INT4 XLALBasicRandom( INT4 i );
RandomParams * XLALCreateRandomParams( INT4 seed );
void XLALDestroyRandomParams( RandomParams *params );
REAL4 XLALUniformDeviate( RandomParams *params );
int XLALNormalDeviates( REAL4Vector *deviates, RandomParams *params );

void
LALCreateRandomParams (
    LALStatus        *status,
    RandomParams **params,
    INT4           seed
    );

void
LALDestroyRandomParams (
    LALStatus        *status,
    RandomParams **params
    );

void
LALUniformDeviate (
    LALStatus    *status,
    REAL4        *deviate,
    RandomParams *params
    );

void
LALNormalDeviates (
    LALStatus    *status,
    REAL4Vector  *deviates,
    RandomParams *params
    );

void LALMersenneRandom(
    LALStatus      *status,
    REAL8          *output,
    MTRandomParams *params
    );

void LALMersenneRandomVector(
    LALStatus      *status,
    REAL8Vector    *output,
    MTRandomParams *params
    );

void LALCreateMTRandomParams(
    LALStatus       *status,
    REAL8            seed,
    MTRandomParams **params
    );

void LALDestroyMTRandomParams(
    LALStatus       *status,
    MTRandomParams **params
    );

#ifdef  __cplusplus
}
#endif

#endif /* _RANDOM_H */

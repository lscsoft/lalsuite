/*<lalVerbatim file="AstroOmegaHV">
Author: Regimbau Tania
$Id$
</lalVerbatim> */

/*<lalLaTeX>

\section{Header \texttt{AstroOmega.h}}
\label{s:AstroOmega.h}

compute the energy density spectrum of stochastic backgrounds produced
by cosmological population of astrophysical sources.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/AstroOmega.h>
\end{verbatim}
</lalLaTeX> */

#ifndef _ASTROOMEGA_H
#define _ASTROOMEGA_H
#include <stdio.h>
#include <math.h>
#include <lal/LALConfig.h>
#include <lal/LALStdlib.h>
#include <lal/Integrate.h>

#ifdef __cplusplus
extern "C" {
#endif
NRCSID (ASTROOMEGAH, "$Id$");

/*<lalLaTeX>
\subsection*{Error conditions}
the errors that may occur in this module are integration errors already defined in Integrate.h  

\subsection*{Structures}
These are function pointers corresponding to the spectral energy density of a single source.
\begin{verbatim}
typedef void (REAL8LALSDensity) (REAL8 *output, REAL8 input);
\end{verbatim}
These are input structures corresponding to the model parameters (the cosmological model parameters and the source model parameters)

\begin{verbatim}
cosmological model parameters:
typedef struct
tagAstroOmegaCosmoParams
 {
   REAL8   ho; Hubble parameter
   REAL8   density_matter; density parameter of matter
   REAL8   density_vacuum; density parameter of vacuum
   REAL8   density_k; density parameter of curvature
 }
AstroOmegaCosmoParams;
\end{verbatim}

source parameters 
typedef struct
tagAstroOmegaSourceParams
 { 
   REAL8LALSDensity   *SDensitySource; single spectral energy density
   REAL8              numax; frequency cutoff in the source frame
   REAL8              lambda; mass fraction of source progenitors expressed in inverse solar masses.
 }
AstroOmegaSourceParams;


model parameters (cosmological + source)
typedef struct
tagAstroOmegaParams
 {
   AstroOmegaCosmoParams          cosmoparams;
   AstroOmegaSourceParams         sourceparams;
   void                           *extraparams;  
 }
AstroOmegaParams;
\end{verbatim}


\vfill{\footnotesize\input{AstroOmegaHV}}
\newpage\input{AstroOmegaC}
</lalLaTeX> */

/*type corresponding to the spectral energy density of a single source*/
typedef void (REAL8LALSDensity) (REAL8 *output, REAL8 input);

/*MODEL PARAMETERS*/

/*cosmological model*/
typedef struct
tagAstroOmegaCosmoParams
 {
   REAL8   ho;
   REAL8   density_matter;
   REAL8   density_vacuum;
   REAL8   density_k;
 }
AstroOmegaCosmoParams;
/*source model*/
/*in the general case, the user should define previously the single spectral energy density*/
typedef struct
tagAstroOmegaSourceParams
 { 
   REAL8LALSDensity   *SDensitySource;
   REAL8              numax;
   REAL8              lambda;
 }
AstroOmegaSourceParams;


typedef struct
tagAstroOmegaParams
 {
   AstroOmegaCosmoParams          cosmoparams;
   AstroOmegaSourceParams         sourceparams;
   void                           *extraparams;  
 }
AstroOmegaParams;



/*functions returning $\Omega _{gw}(\nu _{o})$*/

void
LALAstroOmega (
    LALStatus    *status,
    REAL8        *result,
    REAL8         nu,
    void         *params
    );



#ifdef __cplusplus
}
#endif

#endif /* _ASTROOMEGA_H */

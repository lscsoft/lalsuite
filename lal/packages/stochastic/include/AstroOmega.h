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
the errors that may occurs in this module are integration errors already definein Integrate.h  

\subsection*{Structures}
These are function pointer corresponding to the spectral energy density of a single source
\begin{verbatim}
typedef void (REAL8LALSDensity) (REAL8 *output, REAL8 input);
\end{verbatim}
These are input structures corresponding to the model parameters (the cosmological model parameters and the source model parameters)
cosmological model parameters:
\begin{verbatim}
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
source model parameters
In the general case, the user should define previously the single spectral energy density 
\begin{verbatim}
typedef struct
tagAstroOmegaGeneralSourceParams
 { 
   REAL8LALSDensity   *SDensitySource; single spectral energy density
   REAL8              numax; frequency cutoff in the source frame
   REAL8              fact; groupes various factors, see AstroOmegatemplates.c for details
 }
AstroOmegaGeneralSourceParams;

typedef struct
tagAstroOmegaTemplateSourceParams
 { 
   REAL8   numax;
   REAL8   fact;
 }
AstroOmegaTemplateSourceParams;
\end{verbatim}
model parameters (cosmological + source model)
\begin{verbatim}
typedef struct
tagAstroOmegaGeneralParams
 {
   AstroOmegaCosmoParams          cosmoparams;
   AstroOmegaGeneralSourceParams  gsourceparams;
   void                 *extraparams;  
 }
AstroOmegaGeneralParams;

typedef struct
tagAstroOmegaTemplatesParams
 {
   AstroOmegaCosmoParams             cosmoparams;
   TemplatesSourceParams   tsourceparams;   
   void                    *extraparams; 
}
AstroOmegaTemplatesParams;
\end{verbatim}

\vfill{\footnotesize\input{AstroOmegaHV}}
%\newpage\input{AstroOmegaC}
\newpage\input{AstroOmegaGeneralC}
\newpage\input{AstroOmegaTemplatesC}
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
tagAstroOmegaGeneralSourceParams
 { 
   REAL8LALSDensity   *SDensitySource;
   REAL8              numax;
   REAL8              fact;
 }
AstroOmegaGeneralSourceParams;

typedef struct
tagAstroOmegaTemplateSourceParams
 { 
   REAL8   numax;
   REAL8   fact;
 }
AstroOmegaTemplatesSourceParams;

typedef struct
tagAstroOmegaGeneralParams
 {
   AstroOmegaCosmoParams          cosmoparams;
   AstroOmegaGeneralSourceParams  gsourceparams;
   void                 *extraparams;  
 }
AstroOmegaGeneralParams;

typedef struct
tagAstroOmegaTemplatesParams
 {
   AstroOmegaCosmoParams             cosmoparams;
   AstroOmegaTemplatesSourceParams   tsourceparams;   
   void                    *extraparams; 
}
AstroOmegaTemplatesParams;

/*functions returning $\Omega _{gw}(\nu _{o})$*/

void
LALAstroOmegaSource (
    LALStatus    *status,
    REAL8        *result,
    REAL8         nu,
    void         *params
    );
/*for rotating pulsars, r and bar modes or binary system,the following functions for which the single spectral density has already been intruduced can be used directly*/
void
LALAstroOmegaPulsar (
    LALStatus    *status,
    REAL8        *result,
    REAL8         nu,
    void         *params
    );


void
LALAstroOmegaModes (
    LALStatus    *status,
    REAL8        *result,
    REAL8        nu,
    void         *params
    );

void
LALAstroOmegaBinary (
    LALStatus    *status,
    REAL8        *result,
    REAL8        nu,
    void         *params
    );


#ifdef __cplusplus
}
#endif

#endif /* _ASTROOMEGA_H */

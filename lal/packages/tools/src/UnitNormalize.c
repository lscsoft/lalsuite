/************************************ <lalVerbatim file="UnitNormalizeCV">
Author: Whelan, J. T.
$Id$  
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{UnitNormalize.c}}
\label{tools:ss:UnitNormalize.c}

Brings an \texttt{LALUnit} structure into standard form by reducing
all of the rational exponents into LCD form.

\subsubsection*{Prototypes}
\input{UnitNormalizeCP}
\idx{LALUnitsNormalize()}

\subsubsection*{Description}
 
Since the \texttt{LALUnit} structure stores the rational powers of the
fundamental units as numerator and denominator, it is possible to
represent the same units in different ways, \textit{e.g.}, m$^2$
versus m$^{4/2}$.  This function reduces all of those fractions to
convert the structure to its simplest form.

\subsubsection*{Algorithm}

The rational powers are reduced using Euclid's
algorithm\cite{tools:Geddes:1992}.

\subsubsection*{Uses}

None.

\subsubsection*{Notes}

Note that the functions \texttt{LALUnitRaise()},
\texttt{LALUnitMultiply()}, and \texttt{LALUnitCompare()} all call
\texttt{LALUnitNormalize()} themselves, so there is usually no need to
call it explicitly.

\vfill{\footnotesize\input{UnitNormalizeCV}}

******************************************************* </lalLaTeX> */ 
/**************************************** <lalLaTeX file="UnitNormalizeCB">
\bibitem{tools:Geddes:1992}
K.~O.~Geddes, S.~R.~Czapor, and G.~Labahn, \textit{Algorithms for
computer algebra}.  (Kluwer Academic, Boston, 1992)
******************************************************* </lalLaTeX> */ 

#define TRUE 1
#define FALSE 0

#include <lal/LALStdlib.h>
#include <lal/Units.h>

NRCSID( UNITNORMALIZEC, "$Id$" );

/* this is an implementation of the Euclidean Algorithm */
static UINT2 
gcd(INT2 numer, UINT2 denom) 
{
   UINT2 next_numer,next_denom,remainder;
   next_numer=abs(numer);
   next_denom=denom;
   while(next_denom != 0){
      remainder=next_numer%next_denom;
      next_numer=next_denom;
      next_denom=remainder;
   }
   return abs(next_numer);
}

/* <lalVerbatim file="UnitNormalizeCP"> */
void
LALUnitNormalize (LALStatus *status, LALUnit *output, const LALUnit *input)
/* </lalVerbatim> */
     /* Reduce all the rational powers in an LALUnit structure */
{
  UINT2 commonFactor;
  INT2  numer;
  UINT2 denom;
  UINT2 i;

  INITSTATUS( status, "LALUnitNormalize", UNITNORMALIZEC );

  ASSERT( input != NULL, status, UNITSH_ENULLPIN, UNITSH_MSGENULLPIN );

  ASSERT( output != NULL, status, UNITSH_ENULLPOUT, UNITSH_MSGENULLPOUT );

  output->powerOfTen = input->powerOfTen;
  for (i=0; i<LALNumUnits; ++i)
  {
    numer = input->unitNumerator[i];
    denom = input->unitDenominatorMinusOne[i] + 1;
    commonFactor = gcd (numer, denom);
    output->unitNumerator[i] = numer / commonFactor;
    output->unitDenominatorMinusOne[i] = denom / commonFactor - 1;
  } /* for i */
}

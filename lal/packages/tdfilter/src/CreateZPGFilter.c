/****************************** <lalVerbatim file="CreateZPGFilterCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{CreateZPGFilter.c}}
\label{ss:CreateZPGFilter.c}

Creates ZPG filter objects.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{CreateZPGFilterCP}
\idx{LALCreateCOMPLEX8ZPGFilter()}
\idx{LALCreateCOMPLEX16ZPGFilter()}

\subsubsection*{Description}

These functions create an object \verb@**output@, of type
\verb@COMPLEX8ZPGFilter@ or \verb@COMPLEX16ZPGFilter@, having
\verb@numZeros@ zeros and \verb@numPoles@ poles.  The values of those
zeros and poles are not set by these routines (in general they will
start out as garbage).  The handle passed into the functions must be a
valid handle (i.e.\ \verb@output@$\neq$\verb@NULL@), but must not
point to an existing object (\i.e.\ \verb@*output@=\verb@NULL@).

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
LALMalloc()                     LALFree()
LALCCreateVector()              LALCDestroyVector()
LALZCreateVector()              LALZDestroyVector()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{CreateZPGFilterCV}}

******************************************************* </lalLaTeX> */

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/ZPGFilter.h>

NRCSID(CREATEZPGFILTERC,"$Id$");


COMPLEX8ZPGFilter *XLALCreateCOMPLEX8ZPGFilter( INT4 numZeros, INT4 numPoles )
{
  static const char *func = "XLALCreateCOMPLEX8ZPGFilter";
  COMPLEX8ZPGFilter *output;
  if ( numZeros < 0 || numPoles < 0 )
    XLAL_ERROR_NULL( func, XLAL_EINVAL );
  output = LALCalloc( 1, sizeof(*output) );
  if ( ! output )
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );

  /* Allocate the data fields.  If the number of poles or zeros is 0,
     the corresponding field(s) should remain null. */
  if ( numZeros > 0 )
    if ( ! ( output->zeros = XLALCreateCOMPLEX8Vector( numZeros ) ) )
    {
      XLALDestroyCOMPLEX8ZPGFilter( output );
      XLAL_ERROR_NULL( func, XLAL_EFUNC );
    }
  if ( numPoles > 0 )
    if ( ! ( output->poles = XLALCreateCOMPLEX8Vector( numPoles ) ) )
    {
      XLALDestroyCOMPLEX8ZPGFilter( output );
      XLAL_ERROR_NULL( func, XLAL_EFUNC );
    }

  return output;
}

COMPLEX16ZPGFilter *XLALCreateCOMPLEX16ZPGFilter( INT4 numZeros, INT4 numPoles )
{
  static const char *func = "XLALCreateCOMPLEX16ZPGFilter";
  COMPLEX16ZPGFilter *output;
  if ( numZeros < 0 || numPoles < 0 )
    XLAL_ERROR_NULL( func, XLAL_EINVAL );
  output = LALCalloc( 1, sizeof(*output) );
  if ( ! output )
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );

  /* Allocate the data fields.  If the number of poles or zeros is 0,
     the corresponding field(s) should remain null. */
  if ( numZeros > 0 )
    if ( ! ( output->zeros = XLALCreateCOMPLEX16Vector( numZeros ) ) )
    {
      XLALDestroyCOMPLEX16ZPGFilter( output );
      XLAL_ERROR_NULL( func, XLAL_EFUNC );
    }
  if ( numPoles > 0 )
    if ( ! ( output->poles = XLALCreateCOMPLEX16Vector( numPoles ) ) )
    {
      XLALDestroyCOMPLEX16ZPGFilter( output );
      XLAL_ERROR_NULL( func, XLAL_EFUNC );
    }

  return output;
}



/* <lalVerbatim file="CreateZPGFilterCP"> */
void
LALCreateCOMPLEX8ZPGFilter( LALStatus         *stat,
			    COMPLEX8ZPGFilter **output,
			    INT4              numZeros,
			    INT4              numPoles )
{ /* </lalVerbatim> */
  INITSTATUS(stat,"LALCreateCOMPLEX8ZPGFilter",CREATEZPGFILTERC);

  /* Make sure that the output handle exists, but points to a null
     pointer. */
  ASSERT(output,stat,ZPGFILTERH_ENUL,ZPGFILTERH_MSGENUL);
  ASSERT(!*output,stat,ZPGFILTERH_EOUT,ZPGFILTERH_MSGEOUT);

  /* Make sure that numZeros and numPoles are non-negative. */
  ASSERT(numZeros>=0,stat,ZPGFILTERH_EBAD,ZPGFILTERH_MSGEBAD);
  ASSERT(numPoles>=0,stat,ZPGFILTERH_EBAD,ZPGFILTERH_MSGEBAD);

  /* Create the output structure. */
  *output = XLALCreateCOMPLEX8ZPGFilter( numZeros, numPoles );
  if ( ! *output )
  {
    ABORT(stat,ZPGFILTERH_EMEM,ZPGFILTERH_MSGEMEM);
  }

  /* Normal exit */
  RETURN(stat);
}


/* <lalVerbatim file="CreateZPGFilterCP"> */
void
LALCreateCOMPLEX16ZPGFilter( LALStatus          *stat,
			     COMPLEX16ZPGFilter **output,
			     INT4               numZeros,
			     INT4               numPoles )
{ /* </lalVerbatim> */
  INITSTATUS(stat,"LALCreateCOMPLEX16ZPGFilter",CREATEZPGFILTERC);

  /* Make sure that the output handle exists, but points to a null
     pointer. */
  ASSERT(output,stat,ZPGFILTERH_ENUL,ZPGFILTERH_MSGENUL);
  ASSERT(!*output,stat,ZPGFILTERH_EOUT,ZPGFILTERH_MSGEOUT);

  /* Make sure that numZeros and numPoles are non-negative. */
  ASSERT(numZeros>=0,stat,ZPGFILTERH_EBAD,ZPGFILTERH_MSGEBAD);
  ASSERT(numPoles>=0,stat,ZPGFILTERH_EBAD,ZPGFILTERH_MSGEBAD);

  /* Create the output structure. */
  *output = XLALCreateCOMPLEX16ZPGFilter( numZeros, numPoles );
  if ( ! *output )
  {
    ABORT(stat,ZPGFILTERH_EMEM,ZPGFILTERH_MSGEMEM);
  }

  /* Normal exit */
  RETURN(stat);
}

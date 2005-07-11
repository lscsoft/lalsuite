/** \file
 * \ingroup SkyCoordinates
 * \author Creighton, T. D.
 * \date $Date$
 * \brief Automatically converts among sky coordinate systems.
 * 
 * $Id$
 * 

\par Description

The function \verb+LALConvertSkyCoordinates()+ transforms the contents
of <tt>*input</tt> to the system  
specified in <tt>*params</tt>, storing the result in <tt>*output</tt>
(which may point to the same object as <tt>*input</tt> for an in-place
transformation).  The routine makes calls to the functions in
<tt>CelestialCoordinates.c</tt> and <tt>TerrestrialCoordinates.c</tt> as
required; the <tt>*params</tt> object must store any data fields
required by these functions, or an error will occur.


The function <tt>LALNormalizeSkyPosition()</tt> ``normalizes'' any given
(spherical) sky-position (in radians), which means it projects the
angles into \f$[0, 2\pi) \times [-\pi/2, \pi/2]\f$ if they lie outside.


\par Algorithm

<tt>LALConvertSkyCoordinates()</tt> is structured as a simple loop over
transformations, each of which moves the output sky position one step
closer to the desired final coordinates system.  The usual ``flow'' of
the algorithm is: 
\image html inject_ConvFlow.png
\latexonly
\begin{center}
horizon
\makebox[0pt][l]{\raisebox{0.4ex}{$\rightarrow$}}%
\makebox[0pt][l]{\raisebox{-0.2ex}{$\leftarrow$}}
\quad geographic
\makebox[0pt][l]{\raisebox{0.4ex}{$\rightarrow$}}%
\makebox[0pt][l]{\raisebox{-0.2ex}{$\leftarrow$}}
\quad equatorial
\makebox[0pt][l]{\raisebox{2.2ex}{$\nearrow$}}%
\makebox[0pt][l]{\raisebox{1.8ex}{$\,\swarrow$}}%
\makebox[0pt][l]{\raisebox{-1.8ex}{$\,\searrow$}}%
\makebox[0pt][l]{\raisebox{-2.2ex}{$\nwarrow$}}
\quad
\makebox[0pt][l]{\raisebox{4ex}{ecliptic}}%
\makebox[0pt][l]{\raisebox{-4ex}{Galactic}}
\end{center}
\endlatexonly
although one can also convert directly between equatorial and horizon
coordinate systems if <tt>params->zenith</tt> is given in equatorial
coordinates (i.e.\ if its longitudinal coordinate is the local mean
sidereal time rather than the geographic longitude of the observer).
This leads to the only error checking done within this function: when
transforming to horizon coordinates, it checks that
<tt>params->zenith</tt> is either in sky-fixed equatorial or Earth-fixed
geographic coordinates.  Other than this, error checking is left to
the secondary function call; if a parameter is absent or poorly
formatted, the called function will return an error.

\par Uses
\code 
LALHorizonToSystem()            LALSystemToHorizon()
LALGeographicToEquatorial()     LALEquatorialToGeographic()
LALEquatorialToEcliptic()       LALEclipticToEquatorial()
LALEquatorialToGalactic()       LALGalacticToEquatorial()
\endcode

*/

/*---------- laldoc-version of documentation follows ---------- */

/******************************* <lalVerbatim file="SkyCoordinatesCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{SkyCoordinates.c}}
\label{ss:SkyCoordinates.c}

Automatically converts among sky coordinate systems.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{SkyCoordinatesCP}
\idx{LALConvertSkyCoordinates()}
\idx{LALNormalizeSkyPosition()}

\subsubsection*{Description}

The function \verb+LALConvertSkyCoordinates()+ transforms the contents
of <tt>*input</tt> to the system  
specified in <tt>*params</tt>, storing the result in <tt>*output</tt>
(which may point to the same object as <tt>*input</tt> for an in-place
transformation).  The routine makes calls to the functions in
<tt>CelestialCoordinates.c</tt> and <tt>TerrestrialCoordinates.c</tt> as
required; the <tt>*params</tt> object must store any data fields
required by these functions, or an error will occur.


The function \verb+LALNormalizeSkyPosition()+ ``normalizes'' any given
(spherical) sky-position (in radians), which means it projects the
angles into $[0, 2\pi) \times [-\pi/2, \pi/2]$ if they lie outside.


\subsubsection*{Algorithm}

\verb+LALConvertSkyCoordinates()+ is structured as a simple loop over
transformations, each 
of which moves the output sky position one step closer to the desired
final coordinates system.  The usual ``flow'' of the algorithm is:
\begin{center}
horizon
\makebox[0pt][l]{\raisebox{0.4ex}{$\rightarrow$}}%
\makebox[0pt][l]{\raisebox{-0.2ex}{$\leftarrow$}}
\quad geographic
\makebox[0pt][l]{\raisebox{0.4ex}{$\rightarrow$}}%
\makebox[0pt][l]{\raisebox{-0.2ex}{$\leftarrow$}}
\quad equatorial
\makebox[0pt][l]{\raisebox{2.2ex}{$\nearrow$}}%
\makebox[0pt][l]{\raisebox{1.8ex}{$\,\swarrow$}}%
\makebox[0pt][l]{\raisebox{-1.8ex}{$\,\searrow$}}%
\makebox[0pt][l]{\raisebox{-2.2ex}{$\nwarrow$}}
\quad
\makebox[0pt][l]{\raisebox{4ex}{ecliptic}}%
\makebox[0pt][l]{\raisebox{-4ex}{Galactic}}
\end{center}
although one can also convert directly between equatorial and horizon
coordinate systems if <tt>params->zenith</tt> is given in equatorial
coordinates (i.e.\ if its longitudinal coordinate is the local mean
sidereal time rather than the geographic longitude of the observer).
This leads to the only error checking done within this function: when
transforming to horizon coordinates, it checks that
<tt>params->zenith</tt> is either in sky-fixed equatorial or Earth-fixed
geographic coordinates.  Other than this, error checking is left to
the secondary function call; if a parameter is absent or poorly
formatted, the called function will return an error.

\subsubsection*{Uses}
\begin{verbatim}
LALHorizonToSystem()            LALSystemToHorizon()
LALGeographicToEquatorial()     LALEquatorialToGeographic()
LALEquatorialToEcliptic()       LALEclipticToEquatorial()
LALEquatorialToGalactic()       LALGalacticToEquatorial()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{SkyCoordinatesCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/SkyCoordinates.h>

#define LAL_ALPHAGAL (3.366032942)
#define LAL_DELTAGAL (0.473477302)
#define LAL_LGAL     (0.576)

NRCSID( SKYCOORDINATESC, "$Id$" );

/* <lalVerbatim file="SkyCoordinatesCP"> */
void
LALConvertSkyCoordinates( LALStatus        *stat,
			  SkyPosition      *output,
			  SkyPosition      *input,
			  ConvertSkyParams *params )
{ /* </lalVerbatim> */
  SkyPosition temp; /* temporary sky position (duh) */

  INITSTATUS( stat, "LALConvertSkyCoordinates", SKYCOORDINATESC );
  ATTATCHSTATUSPTR( stat );

  /* Make sure parameter structures exist. */
  ASSERT( input, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );
  ASSERT( output, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );
  ASSERT( params, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );

  /* Start looping! */
  temp = *input;
  while ( temp.system != params->system ) {

    if ( temp.system == COORDINATESYSTEM_HORIZON )
      /* Only one possible direction. */
      TRY( LALHorizonToSystem( stat->statusPtr, &temp, &temp,
			       params->zenith ), stat );

    else if ( temp.system == COORDINATESYSTEM_GEOGRAPHIC ) {
      /* Two possible directions.  But, if headed towards the horizon
         coordinate system, make sure we can get there! */
      if ( params->system == COORDINATESYSTEM_HORIZON ) {
	if ( params->zenith ) {
	  if ( params->zenith->system == COORDINATESYSTEM_GEOGRAPHIC )
	    TRY( LALSystemToHorizon( stat->statusPtr, &temp, &temp,
				     params->zenith ), stat );
	  else if ( params->zenith->system == COORDINATESYSTEM_EQUATORIAL )
	    TRY( LALGeographicToEquatorial( stat->statusPtr, &temp, &temp,
					    params->gpsTime ), stat );
	  else
	    ABORT( stat, SKYCOORDINATESH_ESYS, SKYCOORDINATESH_MSGESYS );
	} else
	  ABORT( stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );
      } else
	TRY( LALGeographicToEquatorial( stat->statusPtr, &temp, &temp,
					params->gpsTime ), stat );
    }

    else if ( temp.system == COORDINATESYSTEM_EQUATORIAL ) {
      /* Up to four possible directions, depending on
         params->zenith->system. */
      if ( ( params->system == COORDINATESYSTEM_HORIZON ) &&
	   ( params->zenith ) &&
	   ( params->zenith->system == COORDINATESYSTEM_EQUATORIAL ) )
	TRY( LALSystemToHorizon( stat->statusPtr, &temp, &temp,
				 params->zenith ), stat );
      else if ( params->system == COORDINATESYSTEM_ECLIPTIC )
	TRY( LALEquatorialToEcliptic( stat->statusPtr, &temp, &temp ),
	     stat );
      else if ( params->system == COORDINATESYSTEM_GALACTIC )
	TRY( LALEquatorialToGalactic( stat->statusPtr, &temp, &temp ),
	     stat );
      else
	TRY( LALEquatorialToGeographic( stat->statusPtr, &temp, &temp,
					params->gpsTime ), stat );
    }

    else if ( temp.system == COORDINATESYSTEM_ECLIPTIC ) {
      /* Only one possible direction. */
      TRY( LALEclipticToEquatorial( stat->statusPtr, &temp, &temp ),
	   stat );
    }

    else if ( temp.system == COORDINATESYSTEM_GALACTIC ) {
      /* Only one possible direction. */
      TRY( LALGalacticToEquatorial( stat->statusPtr, &temp, &temp ),
	   stat );
    }

    else
      ABORT( stat, SKYCOORDINATESH_ESYS, SKYCOORDINATESH_MSGESYS );
  }

  /* We've gotten to the correct coordinate system.  Set the output
     and return. */
  *output = temp;
  DETATCHSTATUSPTR( stat );
  RETURN( stat );

} /* LALConvertSkyCoordinates */



/** If sky-position is not in the canonical range 
 * \f$(\alpha,\delta)\in [0,2\pi]\times[-\pi/2, \pi/2]\f$, normalize it
 * by mapping it into this coordinate-interval.
 * Based on Alicia's function with some additional "unwinding" added.  
 */
void
LALNormalizeSkyPosition (LALStatus *stat, 
			 SkyPosition *posOut, 		/**< [out] normalized sky-position */
			 const SkyPosition *posIn)	/**< [in] general sky-position */
{
  SkyPosition tmp;	/* allow posOut == posIn */

  INITSTATUS( stat, "NormalizeSkyPosition", SKYCOORDINATESC);
  
  ASSERT (posIn, stat, SKYCOORDINATESH_ENUL ,  SKYCOORDINATESH_MSGENUL );
  ASSERT (posOut, stat, SKYCOORDINATESH_ENUL ,  SKYCOORDINATESH_MSGENUL );

  tmp = *posIn;
  
  /* FIRST STEP: completely "unwind" positions, i.e. make sure that 
   * [0 <= alpha < 2pi] and [-pi < delta <= pi] */
  /* normalize longitude */
  while (tmp.longitude < 0)
    tmp.longitude += LAL_TWOPI;
  while (tmp.longitude >= LAL_TWOPI)
    tmp.longitude -= LAL_TWOPI;

  /* pre-normalize (unwind) latitude */
  while (tmp.latitude <= -LAL_PI)
    tmp.latitude += LAL_TWOPI;
  while (tmp.latitude > LAL_TWOPI)
    tmp.latitude -= LAL_TWOPI;

  /* SECOND STEP: get latitude into canonical interval [-pi/2 <= delta <= pi/2 ] */
  /* this requires also a change in longitude by adding/subtracting PI */
  if (tmp.latitude > LAL_PI_2)
    {
      tmp.latitude = LAL_PI - tmp.latitude;
      if (tmp.longitude < LAL_PI)
	{
	  tmp.longitude += LAL_PI;
	}
      else
	{
	  tmp.longitude -= LAL_PI;
	}
    }

  if (tmp.latitude < -LAL_PI_2)
    {
      tmp.latitude = -LAL_PI - tmp.latitude;
      if (tmp.longitude < LAL_PI)
	{
	  tmp.longitude += LAL_PI;
	}
      else
	{
	  tmp.longitude -= LAL_PI;
	}
    }

  *posOut = tmp;

  RETURN (stat);

} /* LALNormalizeSkyPosition() */



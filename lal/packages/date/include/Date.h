/*
<lalVerbatim file="DateHV">

Author: David Chin <dwchin@umich.edu> +1-734-730-1274
$Id$

</lalVerbatim>
*/

/*
<lalLaTeX>

\section{Header \texttt{Date.h}}
\label{s:Date.h}

Provides routines for manipulating date and time information.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/Date.h>
\end{verbatim}

This header covers routines for manipulating date and time
information.  The various time systems are discussed in~cite{esaa:1992}.


</lalLaTeX>
*/

#ifndef _DATE_H
#define _DATE_H

#ifndef _REENTRANT
#   define _REENTRANT
#endif

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <lal/LALRCSID.h>

#include <lal/LALConstants.h>
#include <lal/LALAtomicDatatypes.h>
#include <lal/LALDatatypes.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALStdlib.h>

#include <lal/DetectorSite.h>

#ifdef  __cplusplus
extern "C"
{
#endif

NRCSID (DATEH, "$Id$");

/*
<lalLaTeX>
\subsection*{Error conditions}
</lalLaTeX>
*/

/*
<lalErrTable> 
*/
#define JULIAN_ENULLINPUT    1
#define JULIAN_ENULLOUTPUT   2
#define JULIAN_EDATETOOEARLY 3

#define JULIAN_MSGENULLINPUT "Input is NULL"
#define JULIAN_MSGENULLOUTPUT "Output is NULL"
#define JULIAN_MSGEDATETOOEARLY "Date too early: Julian Day can only be computed for dates >= 1900-Mar"

#define UTOGPS_ENULLINPUT   1
#define UTOGPS_ENULLOUTPUT  2
    
#define UTOGPS_MSGENULLINPUT "Input is NULL"
#define UTOGPS_MSGENULLOUTPUT "Output is NULL"
    

#define UTIME_ENULLINPUT  1
#define UTIME_ENULLOUTPUT 2
#define UTIME_ERANGE      3

#define UTIME_MSGENULLINPUT "Input is NULL"
#define UTIME_MSGENULLOUTPUT "Output is NULL"
#define UTIME_MSGERANGE "Input time out of range: 0 <= utcSeconds <= 946684823"


#define DATESTRING_ENULLINPUT    1
#define DATESTRING_ENULLOUTPUT   2
#define DATESTRING_EBUFFTOOSMALL 3

#define DATESTRING_MSGENULLINPUT "Input is NULL"
#define DATESTRING_MSGENULLOUTPUT "Output is NULL"
#define DATESTRING_MSGEBUFFTOOSMALL "Output timestamp string too small: min. size = 26"
    
#define SECSTOLALDATE_ENULLOUTPUT 1
    
#define SECSTOLALDATE_MSGENULLOUTPUT "Output is NULL"


#define LMST1_ENULLINPUT  1
#define LMST1_ENULLOUTPUT 2
    
#define GMST1_ENULLINPUT  1
#define GMST1_ENULLOUTPUT 2

#define GPSTOLMST1_ENULLINPUT 1
#define GPSTOLMST1_ENULLOUTPUT 2

#define GPSTOGMST1_ENULLINPUT 1
#define GPSTOGMST1_ENULLOUTPUT 2    


#define LMST1_MSGENULLINPUT "Input is NULL"
#define LMST1_MSGENULLOUTPUT "Output is NULL"
    
#define GMST1_MSGENULLINPUT "Input is NULL"
#define GMST1_MSGENULLOUTPUT "Output is NULL"

#define GPSTOLMST1_MSGENULLINPUT "Input is NULL"
#define GPSTOLMST1_MSGENULLOUTPUT "Output is NULL"
    
#define GPSTOGMST1_MSGENULLINPUT "Input is NULL"
#define GPSTOGMST1_MSGENULLOUTPUT "Output is NULL"    

/*
</lalErrTable>

<lalLaTeX>


\subsection*{Structures}


</lalLaTeX>
*/

/*
<lalLaTeX>

\vfill{\footnotesize\input{DateHV}}

</lalLaTeX>
*/

/*
<lalLaTeX>

\subsection*{Types}

\subsubsection*{Enumeration \texttt{LALMSTUnits}}
\index{\texttt{LALMSTUnits}}

This enumerated type is used as a parameter for Mean Sidereal Time
routines to specify the units in which to return the Mean Sidereal
Time. The allowed values are:

\medskip\noindent
\begin{tabular}{ll}
  \verb@MST_SEC@ & arc-seconds \\
  \verb@MST_HRS@ & arc-hours (\textit{i.e.} units of Right Ascension)\\
  \verb@MST_DEG@ & degrees \\
  \verb@MST_RAD@ & radians
\end{tabular}
\bigskip

</lalLaTeX>
*/

typedef enum
{
  MST_SEC,       /* arc seconds */
  MST_HRS,       /* arc hours (i.e. units of Right Ascension) */
  MST_DEG,       /* degrees */
  MST_RAD,       /* radians */
} LALMSTUnits;


/*
<lalLaTeX>

\subsubsection*{Structure \texttt{LALUnixDate}}
\index{\texttt{LALUnixDate}}

This structure is just the standard Unix \texttt{tm} structure.

</lalLaTeX>
*/

/*
 * The standard Unix tm structure
 */
typedef struct
tm
LALUnixDate;

/*
<lalLaTeX>

\subsubsection*{Structure \texttt{LIGOTimeUnix}}
\index{\texttt{LIGOTimeUnix}}

This structure is the Unix-epoch analog of \texttt{LIGOTimeGPS}.  It
store the number of seconds and nanoseconds elapsed since the Unix
epoch (1970-Jan-01 00:00:00). The fileds are:

\begin{description}
\item[\texttt{INT4 unixSeconds}] The integral number of seconds
  elapsed since the Unix epoch
\item[\texttt{INT4 unixNanoSeconds}] The residual number of
  nanoseconds that have to be added to \texttt{unixSeconds} to bring us
  up to the time in question
\end{description}

</lalLaTeX>
*/

/*
 * This time object is exactly like LIGOTimeGPS, except for the name.
 * This measures the amount of time elapsed since the Unix time epoch
 */
typedef struct
tagLIGOTimeUnix
{
    INT4 unixSeconds;
    INT4 unixNanoSeconds;
}
LIGOTimeUnix;

/*
<lalLaTeX>


\subsubsection{Structure \texttt{LALTimeInterval}}

This structure is used for storing intervals of \texttt{LIGOTimeGPS}
and \texttt{LIGOTimeUnix} times.  The fields are:

\begin{description}
\item[\texttt{REAL4 seconds}] Integral part of the time interval
\item[\texttt{REAL4 nanoSeconds}] Residual nanoseconds (\textit{i.e.}
  fractional part, in nanoseconds)
\end{description}

</lalLaTeX>
*/

/*
 * This time object is for time intervals, i.e. no reference epoch implied
 */
typedef struct
tagLALTimeInterval
{
    INT4 seconds;
    INT4 nanoSeconds;
}
LALTimeInterval;
    

/*
 * Encode timezone information
 */
typedef struct
tagLALTimezone
{
    INT4 secondsWest; /* seconds West of UTC */
    INT4 dst;         /* Daylight Savings Time correction to apply */
}
LALTimezone;    

/*
 * Date and time structure
 */
typedef struct
tagLALDate
{
    LALUnixDate unixDate;
    INT4        residualNanoSeconds; /* residual nanoseconds */
    LALTimezone timezone;            /* timezone information */
}
LALDate;

/*
 * Place and time structures
 */
/* First, with GPS */
typedef struct
tagLALPlaceAndGPS
{
    LALDetector *detector;   /* pointer to a detector */
    LIGOTimeGPS *gps;        /* pointer to GPS time */
}
LALPlaceAndGPS;

/* Second, with Date-Time */
typedef struct
tagLALPlaceAndDate
{
    LALDetector *detector;   /* pointer to a detector */
    LALDate     *date;       /* pointer to a date */
}
LALPlaceAndDate;

/*
 * 9. Functions Declarations (i.e., prototypes).
 */
void LALJulianDay(LALStatus*,
                  INT4*,
                  const LALDate*);

void LALModJulianDay(LALStatus*,
                     REAL8*,
                     const LALDate*);

void LALJulianDate(LALStatus*,
                   REAL8*,
                   const LALDate*);

void LALModJulianDate(LALStatus*,
                      REAL8*,
                      const LALDate*);

void LALUtoGPS(LALStatus*,
               LIGOTimeGPS*,
               const LIGOTimeUnix*);

void LALGPStoU(LALStatus*,
               LIGOTimeUnix*,
               const LIGOTimeGPS*);

void LALUtime(LALStatus*,
              LALDate*,
              const LIGOTimeUnix*);

void LALDateString(LALStatus*,
                   CHARVector*,
                   const LALDate*);

void LALGMST1(LALStatus*,
              REAL8*,
              const LALDate*,
              LALMSTUnits);

void LALGPStoGMST1( LALStatus*,
                    REAL8*,
                    const LIGOTimeGPS*,
                    LALMSTUnits);

void LALLMST1(LALStatus*,
              REAL8*,
              const LALPlaceAndDate*,
              LALMSTUnits);

void LALGPStoLMST1( LALStatus*,
                    REAL8*,
                    const LALPlaceAndGPS*,
                    LALMSTUnits);

void LALSecsToLALDate(LALStatus*,
                      LALDate*,
                      REAL8);

#ifdef  __cplusplus
}
#endif

#endif /* _DATE_H */

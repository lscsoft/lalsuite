/* <lalVerbatim file="DateHV">

Author: David Chin <dwchin@umich.edu> +1-734-730-1274
$Id$

</lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{Date.h}}
\label{s:Date.h}

Provides routines for manipulating date and time information.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/Date.h>
\end{verbatim}

This header covers routines for manipulating date and time
information.  The various time systems are discussed in~\cite{esaa:1992}.


</lalLaTeX> */

#ifndef _DATE_H
#define _DATE_H

/* the following two preprocessor defines are to include the prototypes for
 * gmtime_r() and asctime_r() from /usr/include/time.h
 * HOWEVER, they do no good if -ansi is used in gcc: warnings are generated
 * that the prototypes have not been seen */

/* HP-UX and Solaris */
#ifndef _REENTRANT
#   define _REENTRANT
#endif

#ifndef _POSIX_PTHREAD_SEMANTICS
#   define _POSIX_PTHREAD_SEMANTICS
#endif

/* Linux */
#ifndef __USE_POSIX
#   define __USE_POSIX
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

/* <lalLaTeX>
\subsection*{Error conditions}
</lalLaTeX> */

/* <lalErrTable> */
#define DATEH_ENULLINPUT       1
#define DATEH_ENULLOUTPUT      2
#define DATEH_EDATETOOEARLY    3
#define DATEH_ERANGEGPSABS     4
#define DATEH_EBUFFTOOSMALL    5
#define DATEH_EASCTIMEFAIL     6
#define DATEH_EGPSDATETOOEARLY 7
#define DATEH_EFORMATPARAMOUTOFRANGE 8
#define DATEH_EACCPARAMOUTOFRANGE    9

#define DATEH_MSGENULLINPUT "Input is NULL"
#define DATEH_MSGENULLOUTPUT "Output is NULL"
#define DATEH_MSGEDATETOOEARLY "Date too early: Julian Day can only be computed for dates >= 1900-03-01"
#define DATEH_MSGERANGEGPSABS "Input time out of range: only able to accurately convert times between 1980-Jan-06 00:00:00 UTC (GPS 0) and 2002-Mar-31 23:59:00 UTC (GPS 701654353)"
#define DATEH_MSGEBUFFTOOSMALL "Output timestamp string too small: min. size = 26"
#define DATEH_MSGEASCTIMEFAIL "asctimeUNDERSCOREr() failed"
#define DATEH_MSGEGPSDATETOOEARLY "Date too early: GPS time only defined for times on or after 1980-Jan-06 00:00:00 UTC"
#define DATEH_MSGEFORMATPARAMOUTOFRANGE "Format parameter out of range: must be one of LALLEAPSECunderscoreTAIUTC or LALLEAPSECunderscoreGPSUTC"
#define DATEH_MSGEACCPARAMOUTOFRANGE "Accuracy parameter out of range: must be one of LALLEAPSECunderscoreSTRICT or LALLEAPSECunderscoreLOOSE"
  

/* </lalErrTable>

<lalLaTeX>


\subsection*{Structures}


</lalLaTeX> */

/* <lalLaTeX>

\vfill{\footnotesize\input{DateHV}}

</lalLaTeX> */

/* <lalLaTeX>

\subsection*{Types}

\subsubsection*{Enumeration \texttt{LALMSTUnits}}
\idx[Type]{LALMSTUnits}

This enumerated type is used as a parameter for Mean Sidereal Time
routines to specify the units in which to return the Mean Sidereal
Time. The allowed values are:

\medskip\noindent
\begin{tabular}{ll}
  \verb+MST_SEC+ & arc-seconds \\
  \verb+MST_HRS+ & arc-hours (\textit{i.e.} units of Right Ascension)\\
  \verb+MST_DEG+ & degrees \\
  \verb+MST_RAD+ & radians
\end{tabular}
\bigskip

</lalLaTeX> */

typedef enum
{
  MST_SEC,       /* arc seconds */
  MST_HRS,       /* arc hours (i.e. units of Right Ascension) */
  MST_DEG,       /* degrees */
  MST_RAD        /* radians */
} LALMSTUnits;

/* <lalLaTeX>
\subsubsection*{Enumeration \texttt{LALLeapSecAccuracy}}
\idx[Type]{LALLeapSecAccuracy}

This enumerated type is used as a parameter for \texttt{LALGPStoUTC()},
\texttt{LALUTCtoGPS()}, and \texttt{LALLeapSecs()} to specify if complete
accuracy is required in use of leap seconds.  The allowed values are:

\medskip\noindent
\begin{tabular}{ll}
  \verb+LALLEAPSEC_LOOSE+ & may miss leap seconds \\
  \verb+LALLEAPSEC_STRICT+ & require all leap seconds
\end{tabular}
\bigskip

If strict accuracy is selected, the code will \texttt{ABORT} if leap second
data is not current.  Otherwise, a warning will be printed, and the code
will continue execution.

</lalLaTeX> */
typedef enum
{
  LALLEAPSEC_LOOSE,
  LALLEAPSEC_STRICT
} LALLeapSecAccuracy;

/* <lalLaTeX>

\subsubsection*{Enumeration \texttt{LALLeapSecFormat}}
\idx[Type]{LALLeapSecFormat}

This enumerated type is used as a parameter for \texttt{LALLeapSecs()} to
specify whether TAI-UTC or GPS-UTC should be returned.  TAI-UTC is the
total number of leap seconds added to UTC since the TAI epoch.  GPS-UTC is
the total number of leap seconds added since the GPS epoch.  These two
quantities are related by:  TAI-UTC = GPS-UTC + 19.

\medskip\noindent
\begin{tabular}{ll}
  \verb+LALLEAPSEC_TAIUTC+ & return TAI-UTC \\
  \verb+LALLEAPSEC_GPSUTC+ & return GPS-UTC
\end{tabular}
\bigskip

</lalLaTeX> */

typedef enum
{
  LALLEAPSEC_TAIUTC,
  LALLEAPSEC_GPSUTC
} LALLeapSecFormat;


/* <lalLaTeX>

\subsubsection*{Structure \texttt{LALUnixDate}}
\idx[Type]{LALUnixDate}

This structure is just the standard Unix \texttt{tm} structure.  We shall
{\em always} ignore the daylight savings time field, \verb+tm_isdst+.

</lalLaTeX> */

/*
 * The standard Unix tm structure
 */
typedef struct
tm
LALUnixDate;


/* <lalLaTeX>


\subsubsection{Structure \texttt{LALTimeInterval}}
\idx[Type]{LALTimeInterval}

This structure is used for storing intervals of \texttt{LIGOTimeGPS}
and \texttt{LIGOTimeUnix} times.  The fields are:

\begin{description}
\item{\texttt{INT4 seconds}} Integral part of the time interval
\item{\texttt{INT8 nanoSeconds}} Residual nanoseconds (\textit{i.e.}
  fractional part, in nanoseconds)
\end{description}

</lalLaTeX> */

/*
 * This time object is for time intervals, i.e. no reference epoch implied
 */
typedef struct
tagLALTimeInterval
{
    INT4 seconds;
    INT8 nanoSeconds;
}
LALTimeInterval;

/* <lalLaTeX>


\subsubsection{Structure \texttt{LALDate}}
\idx[Type]{LALDate}

This structure is an extension of \texttt{LALUnixDate} to include residual
nanosecond information.  The fields are:

\begin{description}
\item{\texttt{LALUnixDate unixDate}} Unix date in \texttt{struct tm}
  format 
\item{\texttt{INT4 residualNanoSeconds}} Residual nanoseconds
\end{description}
</lalLaTeX> */

/*
 * Date and time structure
 */
typedef struct
tagLALDate
{
    LALUnixDate unixDate;
    INT4        residualNanoSeconds; /* residual nanoseconds */
}
LALDate;


/* <lalLaTeX>


\subsubsection{Structure \texttt{LALPlaceAndGPS}}
\idx[Type]{LALPlaceAndGPS}

This structure stores pointers to a \texttt{LALDetector} and a
\texttt{LIGOTimeGPS}. Its sole purpose is to aggregate these
structures for passing to functions.  The fields are:

\begin{description}
\item{\verb+LALDetector *p_detector+} Pointer to a detector
\item{\verb+LIGOTimeGPS *p_gps+} Pointer to a GPS time structure
\end{description}

</lalLaTeX> */

/*
 * Place and time structures
 */
/* First, with GPS */
typedef struct
tagLALPlaceAndGPS
{
    LALDetector *p_detector;   /* pointer to a detector */
    LIGOTimeGPS *p_gps;        /* pointer to GPS time */
}
LALPlaceAndGPS;

/* <lalLaTeX>


\subsubsection{Structure \texttt{LALPlaceAndDate}}
\idx[Type]{LALPlaceAndDate}

Like \texttt{LALPlaceAndGPS}, this structure aggregates a pointer to a
detector and a pointer to a date.  This is another convenience
structure, used in calling \texttt{LALLMST1()}.  The fields are:

\begin{description}
\item{\verb+LALDetector *p_detector+} Pointer to a detector
\item{\verb+LALDate *p_date+} Pointer to a date
\end{description}

</lalLaTeX> */

/* Second, with Date-Time */
typedef struct
tagLALPlaceAndDate
{
    LALDetector *p_detector;   /* pointer to a detector */
    LALDate     *p_date;       /* pointer to a date */
}
LALPlaceAndDate;

/* <lalLaTeX>

\subsubsection{Structure \texttt{LALLeapSecFormatAndAcc}}
\idx[Type]{LALLeapSecFormatAndAcc}

This structure aggregates the \texttt{LALLeapSecFormat} and
\texttt{LALLeapSecAccuracy} parameters for passing to
\texttt{LALLeapSecs()}.

The \texttt{format} field specifies whether \texttt{LALLeapSecs()} returns
TAI-UTC or GPS-UTC.  The \texttt{accuracy} field specifies whether a
warning/error should be produced if the function is given an input GPS time
that may result in a leap second not being accounted for.

</lalLaTeX> */

typedef struct
tagLALLeapSecFormatAndAcc
{
  LALLeapSecFormat   format;
  LALLeapSecAccuracy accuracy;
}
LALLeapSecFormatAndAcc;


/* 
 * Function prototypes
 */

/* <lalLaTeX>
\newpage\input{JulianC}
</lalLaTeX> */

void LALJulianDay(LALStatus     *status,
                  INT4          *jDay,
                  const LALDate *date);

void LALJulianDate (LALStatus     *status,
                    REAL8         *jDateOut,
                    const LALDate *date);

void LALModJulianDate (LALStatus     *status,
                       REAL8         *modJDate,
                       const LALDate *date);


/* <lalLaTeX>
\newpage\input{DateStringC}
</lalLaTeX> */

void LALDateString (LALStatus     *status,
                    CHARVector    *timestamp,
                    const LALDate *date);


/* <lalLaTeX>
\newpage\input{LMST1C}
</lalLaTeX> */

void LALGMST1 (LALStatus     *status,
               REAL8         *gmst,        /* output - GMST1 */
               const LALDate *date,        /* input  - date and time */
               LALMSTUnits    outunits);   /* GMST1 units */

void LALGPStoGMST1( LALStatus         *status,
                    REAL8             *gmst,      /* output - GMST1 */
                    const LIGOTimeGPS *gps,       /* input - GPS time */
                    LALMSTUnits        outunits); /* GMST1 units */

void LALLMST1 (LALStatus             *status,
               REAL8                 *lmst,          /* output - LMST1 */
               const LALPlaceAndDate *placeAndDate,  /* input -  location
                                                        and date */ 
               LALMSTUnits            outunits);     /* LMST1 units */

void LALGPStoLMST1( LALStatus             *status,
                    REAL8                 *lmst,        /* output - LMST1 */
                    const LALPlaceAndGPS  *placeAndGps, /* input - location and
                                                           GPS */  
                    LALMSTUnits            outunits);   /* LMST1 units */

/* <lalLaTeX>
\newpage\input{SecsToLALDateC}
</lalLaTeX> */

void LALSecsToLALDate(LALStatus *status,
                      LALDate   *date,
                      REAL8     );

/* <lalLaTeX>
\newpage\input{GPStoUTCC}
</lalLaTeX> */
void
LALGPStoUTC (LALStatus                *status,
             LALDate                  *pUtcDate,
             const LIGOTimeGPS        *pGpsTime,
             const LALLeapSecAccuracy *pAccuracy);

void
LALUTCtoGPS (LALStatus                *status,
             LIGOTimeGPS              *pGpsTime,
             const LALDate            *pUtcDate,
             const LALLeapSecAccuracy *pAccuracy);


void
LALLeapSecs (LALStatus                    *status,
             INT4                         *p_leapSecs,
             const LIGOTimeGPS            *p_gpsTime,
             const LALLeapSecFormatAndAcc *p_formatAndAcc);

/* The following 2 functions are from S.J. Berukoff, included at his request */
/* <lalLaTeX>
\newpage\input{GPStoFloatC}
</lalLaTeX> */
void LALGPStoFloat (LALStatus *status,
                    REAL8     *output, 
                    const LIGOTimeGPS *input);

void LALFloatToGPS(LALStatus *status, 
                   LIGOTimeGPS *output, 
                   const REAL8 *input);

/* <lalLaTeX>
\newpage\input{GPStoINT8C}
</lalLaTeX> */
void
LALINT8toGPS ( 
    LALStatus          *status,
    LIGOTimeGPS        *output, 
    const INT8         *input 
    );

void
LALGPStoINT8 ( 
    LALStatus          *status,
    INT8               *output, 
    const LIGOTimeGPS  *input 
    );





#ifdef  __cplusplus
}
#endif

#endif /* _DATE_H */

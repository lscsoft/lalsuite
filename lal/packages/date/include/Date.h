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
 * gmtime_r() and asctime_r() from /usr/include/time.h */

/* HP-UX and Solaris */
#ifndef _REENTRANT
#   define _REENTRANT
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
#define DATEH_ENULLINPUT     1
#define DATEH_ENULLOUTPUT    2
#define DATEH_EDATETOOEARLY  3
#define DATEH_ERANGEGPSTOUTC 4
#define DATEH_ERANGEGPSABS   5
#define DATEH_EBUFFTOOSMALL  6
#define DATEH_EASCTIMEFAIL   7

#define DATEH_MSGENULLINPUT "Input is NULL"
#define DATEH_MSGENULLOUTPUT "Output is NULL"
#define DATEH_MSGEDATETOOEARLY "Date too early: Julian Day can only be computed for dates >= 1900-03-01"
#define DATEH_MSGERANGEGPSTOUTC "Input time out of range: only able to accurately convert times between 1980-Jan-06 00:00:00 UTC (GPS 0) and 2002-Mar-31 23:59:00 UTC (GPS 701654353)"
#define DATEH_MSGERANGEGPSABS "Input time out of range: cannot convert times before 1972-Jan-01 00:00:00 UTC (GPS -252892800)"
#define DATEH_MSGEBUFFTOOSMALL "Output timestamp string too small: min. size = 26"
#define DATEH_MSGEASCTIMEFAIL "asctimeUNDERSCOREr() failed"
    

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
\index{\texttt{LALMSTUnits}}

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
  MST_RAD,       /* radians */
} LALMSTUnits;

/* <lalLaTeX>
\subsubsection*{Enumeration \texttt{LALLeapSecAccuracy}}
\index{\texttt{LALLeapSecAccuracy}}

This enumerated type is used as a parameter for \texttt{LALUtime()} to
specify if complete accuracy is required in use of leap seconds.  The
allowed values are:

\medskip\noindent
\begin{tabular}{ll}
  \verb+LALLEAPSEC_LOOSE+ & may miss leap seconds \\
  \verb+LALLEAPSEC_STRICT+ & require all leap seconds
\end{tabular}
\bigskip

</lalLaTeX> */
typedef enum
{
  LALLEAPSEC_LOOSE,
  LALLEAPSEC_STRICT
} LALLeapSecAccuracy;
    


/* <lalLaTeX>

\subsubsection*{Structure \texttt{LALUnixDate}}
\index{\texttt{LALUnixDate}}

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
\index{\texttt{LALTimeInterval}}

This structure is used for storing intervals of \texttt{LIGOTimeGPS}
and \texttt{LIGOTimeUnix} times.  The fields are:

\begin{description}
\item[\texttt{REAL4 seconds}] Integral part of the time interval
\item[\texttt{REAL4 nanoSeconds}] Residual nanoseconds (\textit{i.e.}
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
    INT4 nanoSeconds;
}
LALTimeInterval;

/* <lalLaTeX>


\subsubsection{Structure \texttt{LALDate}}
\index{\texttt{LALDate}}

This structure is an extension of \texttt{LALUnixDate} to include residual
nanosecond information.  The fields are:

\begin{description}
\item[\texttt{LALUnixDate unixDate}] Unix date in \texttt{struct tm}
  format 
\item[\texttt{INT4 residualNanoSeconds}] Residual nanoseconds
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
\index{\texttt{LALPlaceAndGPS}}

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
\index{\texttt{LALPlaceAndDate}}

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
               REAL8         *gmst,     /* output - GMST1 */
               const LALDate *date,     /* input  - date and time */
               LALMSTUnits    outunits);   /* GMST1 units */

void LALGPStoGMST1( LALStatus         *status,
                    REAL8             *gmst,   /* output - GMST1 */
                    const LIGOTimeGPS *gps,    /* input - GPS time */
                    LALMSTUnits        outunits); /* GMST1 units */

void LALLMST1 (LALStatus             *status,
               REAL8                 *lmst,            /* output - LMST1 */
               const LALPlaceAndDate *place_and_date,  /* input -
                                                            location
                                                            and date */ 
               LALMSTUnits            outunits);         /* LMST1 units */

void LALGPStoLMST1( LALStatus             *status,
                    REAL8                 *lmst,      /* output - LMST1 */
                    const LALPlaceAndGPS  *place_and_gps, /* input -
                                                               location and
                                                               GPS */  
                    LALMSTUnits            outunits);       /* LMST1 units */

/* <lalLaTeX>
\newpage\input{SecsToLALDateC}
</lalLaTeX> */

void LALSecsToLALDate(LALStatus*,
                      LALDate*,
                      REAL8);



/* FOOBAR! Put LALLATEX stuff here for GPStoUTC */
void
LALGPStoUTC (LALStatus                *status,
             LALDate                  *p_utcDate,
             const LIGOTimeGPS        *p_gpsTime,
             const LALLeapSecAccuracy *p_accuracy);

void
LALUTCtoGPS (LALStatus                *status,
             LIGOTimeGPS              *p_gpsTime,
             const LALDate            *p_utcDate,
             const LALLeapSecAccuracy *p_accuracy);


#ifdef  __cplusplus
}
#endif

#endif /* _DATE_H */

/* <lalVerbatim file="DateStringCV">
Author: David Chin <dwchin@umich.edu> +1-734-730-1274
$Id$
</lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{DateString.c}}
\label{ss:DateString.c}

Returns a formatted string for the date and time in ISO 8601 format,
given the date in an \texttt{LALDate} structure.

\subsection*{Prototypes}
\vspace{0.1in}
\input{DateStringCP}
\index{\texttt{LALDateString()}}

\subsubsection*{Description}

Returns a formatted string for the date and time in ISO 8601 format,
given the date in an \texttt{LALDate} structure.  A date and time in
ISO 8601 format looks like 2001-03-04 17:03:52 for March 4, 2001,
5:03:52 pm.  The advantage of this format, besides avoiding Y2K
issues, is that a numerical-order sort of dates will result in a
chronologically ordered sort.  This routine is an interface to the
standard Unix function \texttt{strftime(3C)}.

\subsubsection*{Algorithms}

\subsubsection*{Uses}

Suppose we would like to form a timestamp string for the current time.
The following program (taken from GRASP~\cite{grasp:194}) would
accomplish this:

\begin{verbatim}
void printone(Status *status, const LIGOTimeUnix *time1)
{
    LIGOTimeGPS  gpstime;
    LALDate      laldate;
    LIGOTimeUnix tmp;
    LALUnixDate  utimestruct;
    CHARVector *utc = NULL;

    INITSTATUS (status, "printone", TESTUTOGPSC);

    CHARCreateVector(status, &utc, (UINT4)64);

    UtoGPS(status, &gpstime, time1);

    Utime(status, &laldate, time1);
    DateString(status, utc, &laldate);

    printf("%s\n", utc->data);

    CHARDestroyVector(status, &utc);
    
    RETURN (status);
}
\end{verbatim}

\subsubsection*{Notes}

See the official ISO document at
\url{http://www.iso.ch/markete/8601.pdf}, and
an overview of the format at
\url{http://www.cl.cam.ac.uk/~mgk25/iso-time.html}. 

</lalLaTeX> */

#include <lal/LALRCSID.h>
#include <lal/Date.h>


NRCSID (DATESTRINGC, "$Id$");


/* <lalVerbatim file="DateStringCP"> */
void
LALDateString (LALStatus     *status,
               CHARVector    *timestamp,
               const LALDate *date)
{ /* </lalVerbatim> */
  
  INITSTATUS (status, "LALDateString", DATESTRINGC);

  /*
   * Check pointer to input variable
   */
  ASSERT (date != (LALDate *)NULL, status,
          DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  /*
   * Check pointer to output variable.
   */
  ASSERT (timestamp != (CHARVector *)NULL, status,
          DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT);

  /*
   * Check that timestamp buffer is large enough
   */
  ASSERT (timestamp->length >= 26, status,
          DATEH_EBUFFTOOSMALL, DATEH_MSGEBUFFTOOSMALL);

  /*
   * Use strftime (3) to form ISO8601-format time stamp, plus day name
   */
  strftime(timestamp->data, timestamp->length,
           "%Y-%m-%d %H:%M:%S UTC %a", &(date->unixDate));

  RETURN (status);
} /* END LALDateString() */

    

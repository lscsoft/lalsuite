/*
*  Copyright (C) 2007 Curt Cutler, Jolien Creighton, Reinhard Prix, Teviet Creighton, Bernd Machenschalk
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/**************************** <lalVerbatim file="LALInitBarycenterCV">
Author: Cutler, C.
$Id$
*********************************** </lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{LALInitBarycenter.c}}
\label{ss:LALInitBarycenter.c}

Reads Earth and Sun position information from data files.

\subsubsection*{Prototypes}
\input{LALInitBarycenterCP}
\idx{LALInitBarycenter()}

\subsubsection*{Description}

\verb@LALInitBarycenter()@ fills the contents of \verb@edat@ from data
read from data files.  See \verb@LALBarycenter.h@ in the \verb@pulsar@
package for the definition of the \verb@EphemerisData@ structure.

The function reads in two data files (specified in the
\verb@edat->ephiles@ structure) that contain the position, velocity,
and acceleration of the Earth and Sun, respectively, at regular
intervals througout the specified year. E.g., for 1998, the two files
are \verb@earth98.dat@ and \verb@sun98.dat@.  These files are derived
from the JPL DE405 ephemeris and are provided by Cutler.  The first
line of these files specifies the start time, sampling interval, and
number of datapoints stored in each file, which are used to allocate
data arrays \verb@edat->ephemE@ and \verb@edat->ephemS@ of appropriate
length.  \verb@LALInitBarycenter()@ should be called once near the
beginning of the larger analysis package, and the fields
\verb@edat->ephemE@ and \verb@edat->ephemS@ should be freed with
\verb@LALFree()@ near the end.  See the \verb@LALBarycenterTest@
program in the \verb@pulsar@ package for an illustration of how this
routine is used.


\subsubsection*{Uses}
\begin{verbatim}
LALOpenDataFile()
LALMalloc()
LALFree()
\end{verbatim}

\subsubsection*{Notes}
\vfill{\footnotesize\input{LALInitBarycenterCV}}

</lalLaTeX> */

#include <lal/FileIO.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>

NRCSID(LALINITBARYCENTERC,"$Id$");

#define ERRMSGLEN 512
CHAR errmsg[ERRMSGLEN];	/* string-buffer for more explicit error-messages */


/* <lalVerbatim file="LALInitBarycenterCP"> */
void
LALInitBarycenter(LALStatus *stat, EphemerisData *edat)
{ /* </lalVerbatim> */

    FILE *fp1, *fp2; /* fp1 is table of Earth location; fp2 is for Sun*/
    CHAR dummy;
    INT4 j; /*dummy index*/
    INT4 gpsYr; /*gpsYr + leap is the time on the GPS clock
                          at first instant of new year, UTC; equivalently
                          leap is # of leap secs added between Jan.6, 1980 and
                          Jan. 2 of given year */
    INT4 ret; /* return value for checking */

    INITSTATUS(stat,"LALInitBarycenter",LALINITBARYCENTERC);
    ATTATCHSTATUSPTR(stat);

    /* open earth file */
    fp1 = LALOpenDataFile(edat->ephiles.earthEphemeris);

    /* check that we could open the file */
    if ( fp1 == NULL ) {
      snprintf (errmsg, ERRMSGLEN, "%s '%s'\n", LALINITBARYCENTERH_MSGEOPEN, edat->ephiles.earthEphemeris);
      errmsg[ERRMSGLEN-1] = '\0';
      ABORT (stat, LALINITBARYCENTERH_EOPEN, errmsg);
    }

    /* read first line */
    ret = fscanf(fp1,"%d %le %d\n", &gpsYr, &edat->dtEtable, &edat->nentriesE);
    if (ret != 3) {
      fclose(fp1);
      LALPrintError("couldn't parse first line of %s: %d\n", edat->ephiles.earthEphemeris, ret);
      ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
    }

    /* allocate memory for ephemeris info */
    edat->ephemE  = (PosVelAcc *)LALMalloc(edat->nentriesE*sizeof(PosVelAcc));
    if (edat->ephemE == NULL) {
      fclose(fp1);
      ABORT(stat, LALINITBARYCENTERH_EMEM, LALINITBARYCENTERH_MSGEMEM);
    }

    /* first column in earth.dat or sun.dat is gps time--one long integer
       giving the number of secs that have ticked since start of GPS epoch
       +  on 1980 Jan. 6 00:00:00 UTC
    */

    /* read the remaining lines */
    for (j=0; j < edat->nentriesE; ++j) {
      ret = fscanf(fp1,"%le %le %le %le %le %le %le %le %le %le\n",
		   &edat->ephemE[j].gps,
		   &edat->ephemE[j].pos[0], &edat->ephemE[j].pos[1], &edat->ephemE[j].pos[2],
		   &edat->ephemE[j].vel[0], &edat->ephemE[j].vel[1], &edat->ephemE[j].vel[2],
		   &edat->ephemE[j].acc[0], &edat->ephemE[j].acc[1], &edat->ephemE[j].acc[2]);

      /* check number of scanned items */
      if (ret != 10) {
	fclose(fp1);
	LALFree(edat->ephemE);
	LALPrintError("Couldn't parse line %d of %s: %d\n", j+2, edat->ephiles.earthEphemeris, ret);
	ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
      }

      /* check timestamps */
      if(j == 0) {
	if (gpsYr - edat->ephemE[j].gps > 3600 * 24 * 365) {
	  LALPrintError("Wrong timestamp in line %d of %s: %d/%le\n",
			j+2, edat->ephiles.earthEphemeris, gpsYr, edat->ephemE[j].gps);
	  fclose(fp1);
	  LALFree(edat->ephemE);
	  ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
	}
      } else {
	if (edat->ephemE[j].gps != edat->ephemE[j-1].gps + edat->dtEtable) {
	  LALPrintError("Wrong timestamp in line %d of %s: %le/%le\n",
			j+2, edat->ephiles.earthEphemeris, edat->ephemE[j].gps, edat->ephemE[j-1].gps + edat->dtEtable);
	  fclose(fp1);
	  LALFree(edat->ephemE);
	  ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
	}
      }

      /* check position, velocity and acceleration */
#ifndef SQR
#define SQR(x) ((x)*(x))
#endif
      {
	REAL8 length;
	length = sqrt(SQR(edat->ephemE[j].pos[0]) + SQR(edat->ephemE[j].pos[1]) + SQR(edat->ephemE[j].pos[2]));
	if (abs(499.0 - length) > 25) /* 5% */ {
	  LALPrintError("earth position out of range in line %d of %s: %le %le %le: %le\n",
			j+2, edat->ephiles.earthEphemeris,
			edat->ephemE[j].pos[0], edat->ephemE[j].pos[1], edat->ephemE[j].pos[2], length);
	  fclose(fp1);
	  LALFree(edat->ephemE);
	  ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
	}
	length = sqrt(SQR(edat->ephemE[j].vel[0]) + SQR(edat->ephemE[j].vel[1]) + SQR(edat->ephemE[j].vel[2]));
	if (abs(1e-4 - length) > 1e-5) /* 10% */ {
	  LALPrintError("earth velocity out of range in line %d of %s: %le %le %le: %le\n",
			j+2, edat->ephiles.earthEphemeris,
			edat->ephemE[j].vel[0], edat->ephemE[j].vel[1], edat->ephemE[j].vel[2], length);
	  fclose(fp1);
	  LALFree(edat->ephemE);
	  ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
	}
	length = sqrt(SQR(edat->ephemE[j].acc[0]) + SQR(edat->ephemE[j].acc[1]) + SQR(edat->ephemE[j].acc[2]));
	if (abs(2e-11 - length) > 3e-12) /* 15% */ {
	  LALPrintError("earth acceleration out of range in line %d of %s: %le %le %le: %le\n",
			j+2, edat->ephiles.earthEphemeris,
			edat->ephemE[j].acc[0], edat->ephemE[j].acc[1], edat->ephemE[j].acc[2], length);
	  fclose(fp1);
	  LALFree(edat->ephemE);
	  ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
	}
      }

      /* debug
      {
	REAL8 length;
	fprintf(stderr,"earth line: %d:  ");
	length = SQR(edat->ephemE[j].pos[0]) + SQR(edat->ephemE[j].pos[1]) + SQR(edat->ephemE[j].pos[2]);
	fprintf(stderr,"pos: %le (%le), ", sqrt(length), length);
	length = SQR(edat->ephemE[j].vel[0]) + SQR(edat->ephemE[j].vel[1]) + SQR(edat->ephemE[j].vel[2]);
	fprintf(stderr,"vel: %le (%le), ", sqrt(length), length);
	length = SQR(edat->ephemE[j].acc[0]) + SQR(edat->ephemE[j].acc[1]) + SQR(edat->ephemE[j].acc[2]);
	fprintf(stderr,"acc: %le (%le)\n", sqrt(length), length);
      }
      */
    }

    if (fscanf(fp1,"%c",&dummy) != EOF) {
      LALPrintError("Garbage at end of ephemeris file %s\n", edat->ephiles.earthEphemeris);
      fclose(fp1);
      LALFree(edat->ephemE);
      ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
    }

    /* close earth file */
    fclose(fp1);


    /* open sun file */
    fp2 = LALOpenDataFile(edat->ephiles.sunEphemeris);

    /* check that we could open the file */
    if ( fp2 == NULL ) {
      LALFree(edat->ephemE);
      snprintf (errmsg, ERRMSGLEN, "%s '%s'\n", LALINITBARYCENTERH_MSGEOPEN, edat->ephiles.sunEphemeris);
      errmsg[ERRMSGLEN-1] = 0;
      ABORT (stat, LALINITBARYCENTERH_EOPEN, errmsg);
    }

    /* read first line */
    ret = fscanf(fp2,"%d %le %d\n", &gpsYr, &edat->dtStable, &edat->nentriesS);
    if (ret != 3) {
      LALFree(edat->ephemE);
      fclose(fp2);
      LALPrintError("Couldn't parse first line of %s: %d\n", edat->ephiles.sunEphemeris, ret);
      ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
    }

    /* allocate memory for ephemeris info */
    edat->ephemS  = (PosVelAcc *)LALMalloc(edat->nentriesS*sizeof(PosVelAcc));
    if (edat->ephemS == NULL) {
      fclose(fp2);
      LALFree(edat->ephemE);
      ABORT(stat, LALINITBARYCENTERH_EMEM, LALINITBARYCENTERH_MSGEMEM);
    }

    /* read the remaining lines */
    for (j=0; j < edat->nentriesS; ++j) {
      ret = fscanf(fp2,"%le %le %le %le %le %le %le %le %le %le\n",
		   &edat->ephemS[j].gps,
		   &edat->ephemS[j].pos[0], &edat->ephemS[j].pos[1], &edat->ephemS[j].pos[2],
		   &edat->ephemS[j].vel[0], &edat->ephemS[j].vel[1], &edat->ephemS[j].vel[2],
		   &edat->ephemS[j].acc[0], &edat->ephemS[j].acc[1], &edat->ephemS[j].acc[2]);

      /* check number of scanned items */
      if (ret != 10) {
	fclose(fp2);
	LALFree(edat->ephemE);
	LALFree(edat->ephemS);
	LALPrintError("Couldn't parse line %d of %s: %d\n", j+2, edat->ephiles.sunEphemeris, ret);
	ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
      }

      /* check timestamps */
      if(j == 0) {
	if (gpsYr - edat->ephemS[j].gps > 3600 * 24 * 365) {
	  LALPrintError("Wrong timestamp in line %d of %s: %d/%le\n",
			j+2, edat->ephiles.sunEphemeris, gpsYr, edat->ephemS[j].gps);
	  fclose(fp2);
	  LALFree(edat->ephemE);
	  LALFree(edat->ephemS);
	  ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
	}
      } else {
	if (edat->ephemS[j].gps != edat->ephemS[j-1].gps + edat->dtStable) {
	  LALPrintError("Wrong timestamp in line %d of %s: %le/%le\n",
			j+2, edat->ephiles.sunEphemeris, edat->ephemS[j].gps, edat->ephemS[j-1].gps + edat->dtStable);
	  fclose(fp2);
	  LALFree(edat->ephemE);
	  LALFree(edat->ephemS);
	  ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
	}
      }

      /* check position, velocity and acceleration */
      {
	REAL8 length;
	length = sqrt(SQR(edat->ephemS[j].pos[0]) + SQR(edat->ephemS[j].pos[1]) + SQR(edat->ephemS[j].pos[2]));
	if ((1 > length) || (length > 10)) {
	  LALPrintError("sun position out of range in line %d of %s: %f %f %f: %f\n",
			j+2, edat->ephiles.earthEphemeris,
			edat->ephemS[j].pos[0], edat->ephemS[j].pos[1], edat->ephemS[j].pos[2], length);
	  fclose(fp2);
	  LALFree(edat->ephemS);
	  ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
	}
	length = sqrt(SQR(edat->ephemS[j].vel[0]) + SQR(edat->ephemS[j].vel[1]) + SQR(edat->ephemS[j].vel[2]));
	if ((1e-8 > length) || (length > 1e-7)) {
	  LALPrintError("sun velocity out of range in line %d of %s: %f %f %f: %f\n",
			j+2, edat->ephiles.earthEphemeris,
			edat->ephemS[j].vel[0], edat->ephemS[j].vel[1], edat->ephemS[j].vel[2], length);
	  fclose(fp2);
	  LALFree(edat->ephemS);
	  ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
	}
	length = sqrt(SQR(edat->ephemS[j].acc[0]) + SQR(edat->ephemS[j].acc[1]) + SQR(edat->ephemS[j].acc[2]));
	if ((1e-16 > length) || (length > 1e-14)) {
	  LALPrintError("sun acceleration out of range in line %d of %s: %f %f %f: %f\n",
			j+2, edat->ephiles.earthEphemeris,
			edat->ephemS[j].acc[0], edat->ephemS[j].acc[1], edat->ephemS[j].acc[2], length);
	  fclose(fp2);
	  LALFree(edat->ephemS);
	  ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
	}
      }
      /* debug
      {
	REAL8 length;
	fprintf(stderr,"sun line: %d:  ");
	length = SQR(edat->ephemS[j].pos[0]) + SQR(edat->ephemS[j].pos[1]) + SQR(edat->ephemS[j].pos[2]);
	fprintf(stderr,"pos: %le (%le), ", sqrt(length), length);
	length = SQR(edat->ephemS[j].vel[0]) + SQR(edat->ephemS[j].vel[1]) + SQR(edat->ephemS[j].vel[2]);
	fprintf(stderr,"vel: %le (%le), ", sqrt(length), length);
	length = SQR(edat->ephemS[j].acc[0]) + SQR(edat->ephemS[j].acc[1]) + SQR(edat->ephemS[j].acc[2]);
	fprintf(stderr,"acc: %le (%le)\n", sqrt(length), length);
      }
      */
    }

    if (fscanf(fp2,"%c",&dummy) != EOF) {
      LALPrintError("Garbage at end of ephemeris file %s\n", edat->ephiles.sunEphemeris);
      fclose(fp2);
      LALFree(edat->ephemE);
      LALFree(edat->ephemS);
      ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
    }

    /* close the file */
    fclose(fp2);

    /* successful return */
    DETATCHSTATUSPTR(stat);
    RETURN(stat);
}

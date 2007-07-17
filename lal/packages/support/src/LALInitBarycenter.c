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

    INT4 j; /*dummy index*/
    INT4 gpsYr; /*gpsYr + leap is the time on the GPS clock
                          at first instant of new year, UTC; equivalently
                          leap is # of leap secs added between Jan.6, 1980 and
                          Jan. 2 of given year */
    INT4 ret; /* return value for checking */

    INITSTATUS(stat,"LALInitBarycenter",LALINITBARYCENTERC);
    ATTATCHSTATUSPTR(stat);

    fp1 = LALOpenDataFile(edat->ephiles.earthEphemeris);  
    fp2 = LALOpenDataFile(edat->ephiles.sunEphemeris);  

    /*
    fp1 = fopen(edat->ephiles.earthEphemeris,"r");  
    fp2 = fopen(edat->ephiles.sunEphemeris,"r");  
    */

    /* CHECK THAT fp1 and fp2 are not NULL: */
    if ( ( fp1 == NULL ) || ( fp2 == NULL ) ) 
      {
	LALSnprintf (errmsg, ERRMSGLEN, "%s '%s'\n", LALINITBARYCENTERH_MSGEOPEN, 
		     fp1 == NULL ? edat->ephiles.earthEphemeris : edat->ephiles.sunEphemeris);
	errmsg[ERRMSGLEN-1] = 0;
	if (fp1) close (fp1);
	if (fp2) close (fp2);
	ABORT (stat, LALINITBARYCENTERH_EOPEN, errmsg);
      }
    
    /* reading first line of each file */

    ret = fscanf(fp1,"%d %le %d\n", &gpsYr, &edat->dtEtable, &edat->nentriesE);
    if (ret != 3) {
      fclose(fp1);
      fclose(fp2);
      LALPrintError("couldn't parse first line of %s: %d\n", edat->ephiles.earthEphemeris, ret);
      ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
    }
    ret = fscanf(fp2,"%d %le %d\n", &gpsYr, &edat->dtStable, &edat->nentriesS);
    if (ret != 3) {
      fclose(fp1);
      fclose(fp2);
      LALPrintError("couldn't parse first line of %s: %d\n", edat->ephiles.sunEphemeris, ret);
      ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
    }
    
    edat->ephemE  = (PosVelAcc *)LALMalloc(edat->nentriesE*sizeof(PosVelAcc)); 
    if (edat->ephemE == NULL) {
      ABORT(stat, LALINITBARYCENTERH_EMEM, LALINITBARYCENTERH_MSGEMEM);
    }
    edat->ephemS  = (PosVelAcc *)LALMalloc(edat->nentriesS*sizeof(PosVelAcc)); 
    if (edat->ephemS == NULL) {
      LALFree(edat->ephemE);
      ABORT(stat, LALINITBARYCENTERH_EMEM, LALINITBARYCENTERH_MSGEMEM);
    }

    
    /* first column in earth.dat or sun.dat is gps time--one long integer
       giving the number of secs that have ticked since start of GPS epoch
       +  on 1980 Jan. 6 00:00:00 UTC 
    */

    for (j=0; j < edat->nentriesE; ++j) {
      ret = fscanf(fp1,"%le %le %le %le %le %le %le %le %le %le\n",
		   &edat->ephemE[j].gps,    &edat->ephemE[j].pos[0],
		   &edat->ephemE[j].pos[1], &edat->ephemE[j].pos[2],
		   &edat->ephemE[j].vel[0], &edat->ephemE[j].vel[1],
		   &edat->ephemE[j].vel[2], &edat->ephemE[j].acc[0],
		   &edat->ephemE[j].acc[1], &edat->ephemE[j].acc[2]);
      if (ret != 10) {
	fclose(fp1);
	fclose(fp2);
	LALFree(edat->ephemE);
	LALFree(edat->ephemS);
	LALPrintError("couldn't parse line %d of %s: %d\n", j+2, edat->ephiles.earthEphemeris, ret);
	ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
      }
    }

    for (j=0; j < edat->nentriesS; ++j) {
      ret = fscanf(fp2,"%le %le %le %le %le %le %le %le %le %le\n",
		   &edat->ephemS[j].gps,    &edat->ephemS[j].pos[0],
		   &edat->ephemS[j].pos[1], &edat->ephemS[j].pos[2],
		   &edat->ephemS[j].vel[0], &edat->ephemS[j].vel[1],
		   &edat->ephemS[j].vel[2], &edat->ephemS[j].acc[0],
		   &edat->ephemS[j].acc[1], &edat->ephemS[j].acc[2]);
      if (ret != 10) {
	fclose(fp1);
	fclose(fp2);
	LALFree(edat->ephemE);
	LALFree(edat->ephemS);
	LALPrintError("couldn't parse line %d of %s: %d\n", j+2, edat->ephiles.sunEphemeris, ret);
	ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
      }
    }
       
       
    /* Test to make sure last entries for gpsE,S are 
       reasonable numbers; specifically, checking that they are 
       of same order as t2000 
       --within factor e  
       Machine dependent: to be fixed!
      
       if ( (fabs(log(1.e0*(edat->ephemE[edat->nentriesE -1].gps)/t2000)) > 1.e0) 
       ||(fabs(log(1.e0*(edat->ephemS[edat->nentriesS -1].gps)/t2000)) > 1.e0) ){
       LALFree(edat->ephemE);
       LALFree(edat->ephemS);
       ABORT(stat,LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE); 
       }
    */

    fclose(fp1); 
    fclose(fp2);       
    
    /*curt: is below the right return???*/
    DETATCHSTATUSPTR(stat);
    RETURN(stat);
}

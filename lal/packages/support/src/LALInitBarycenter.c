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
		 
/* <lalVerbatim file="LALInitBarycenterCP"> */
void 
LALInitBarycenter(LALStatus *stat, EphemerisData *edat)
{ /* </lalVerbatim> */

    FILE *fp1, *fp2; /* fp1 is table of Earth location; fp2 is for Sun*/

    INT4 j; /*dummy index*/
    INT4 gpsYr;  /*gpsYr + leap is the time on the GPS clock
                          at first instant of new year, UTC; equivalently
                          leap is # of leap secs added between Jan.6, 1980 and
                          Jan. 2 of given year */

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
	CHAR msg[512];
	LALSnprintf (msg, 512, "%s '%s'\n", LALINITBARYCENTERH_MSGEOPEN, 
		     fp1 == NULL ? edat->ephiles.earthEphemeris : edat->ephiles.sunEphemeris);
	msg[511] = 0;
	ABORT (stat, LALINITBARYCENTERH_EOPEN, msg);
      }
    
/*reading first line of each file */

	 fscanf(fp1,"%d %le %d\n", &gpsYr, &edat->dtEtable, &edat->nentriesE);
	 fscanf(fp2,"%d %le %d\n", &gpsYr, &edat->dtStable, &edat->nentriesS);
    
     edat->ephemE  = (PosVelAcc *)LALMalloc(edat->nentriesE*sizeof(PosVelAcc)); 
     edat->ephemS  = (PosVelAcc *)LALMalloc(edat->nentriesS*sizeof(PosVelAcc)); 

/*first column in earth.dat or sun.dat is gps time--one long integer
  giving the number of secs that have ticked since start of GPS epoch
  +  on 1980 Jan. 6 00:00:00 UTC 
*/
       for (j=0; j < edat->nentriesE; ++j){

	 /* check return value of fscanf */
	 fscanf(fp1,"%le %le %le %le %le %le %le %le %le %le\n",
		&edat->ephemE[j].gps,&edat->ephemE[j].pos[0],
		&edat->ephemE[j].pos[1],&edat->ephemE[j].pos[2],
		&edat->ephemE[j].vel[0],&edat->ephemE[j].vel[1],
		&edat->ephemE[j].vel[2],&edat->ephemE[j].acc[0],
		&edat->ephemE[j].acc[1],
		&edat->ephemE[j].acc[2] );
       }

       for (j=0; j < edat->nentriesS; ++j){
	 fscanf(fp2,"%le %le %le %le %le %le %le %le %le %le\n",
		&edat->ephemS[j].gps,&edat->ephemS[j].pos[0],
		&edat->ephemS[j].pos[1],&edat->ephemS[j].pos[2],
		&edat->ephemS[j].vel[0],&edat->ephemS[j].vel[1],
		&edat->ephemS[j].vel[2],&edat->ephemS[j].acc[0],
		&edat->ephemS[j].acc[1],
		&edat->ephemS[j].acc[2] );
       }
       
       
/*Test to make sure last entries for gpsE,S are 
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





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
\index{\texttt{LALInitBarycenter()}}

\subsubsection*{Description}

\verb@LALInitBarycenter()@ reads in two data files (specified in the
\verb@EphemerisFilenames@ structure) that contain the position, velocity,
and acceleration of the Earth and Sun, resp., at regular intervals
througout the specified year. E.g., for 1998, the two files
are \verb@earth98.dat@ and \verb@sun98.dat@. These files are derived from the
JPL DE405 ephemeris and are provided by Cutler. 
The information gets stored in the \verb@EphemerisData@ structure.
This function should be called once
near the beginning of the larger analysis package, as
illustrated in \verb@LALBarycenterTest.c@.


\subsubsection*{Uses}
\begin{verbatim}
LALFopen()
LALFclose()
\end{verbatim}

\subsubsection*{Notes}
\vfill{\footnotesize\input{LALInitBarycenterCV}}

</lalLaTeX> */

#include "LALInitBarycenter.h"
#include <lal/FileIO.h>
		 
NRCSID(LALINITBARYCENTERC,"$Id$");
		 
/* <lalVerbatim file="LALInitBarycenterCP"> */
void 
LALInitBarycenter(LALStatus *stat, EphemerisData *edat)
{ /* </lalVerbatim> */

    FILE *fp1, *fp2; /* fp1 is table of Earth location; fp2 is for Sun*/

    INT4 j; /*dummy index*/
    INT4 gpsYr, leap;  /*gpsYr + leap is the time on the GPS clock
                          at first instant of new year, UTC; equivalently
                          leap is # of leap secs added between Jan.6, 1980 and
                          Jan. 2 of given year */
    INT4 nentriesE, nentriesS; /*number of rows in the earth, sun ephem */
    INT4 t2000 = 630720013;

    INITSTATUS(stat,"LALInitBarycenter",LALINITBARYCENTERC);
    ATTATCHSTATUSPTR(stat);

    fp1 = LALOpenDataFile(edat->ephiles.earthEphemeris);  
    fp2 = LALOpenDataFile(edat->ephiles.sunEphemeris);  

    /* CHECK THAT fp1 and fp2 are not NULL: */
    if ( ( fp1 == NULL ) || ( fp2 == NULL ) ) {
      ABORT (stat, LALINITBARYCENTERH_EOPEN, LALINITBARYCENTERH_MSGEOPEN);
    }
    
/*reading first line of each file */

	 fscanf(fp1,"%d %d %d\n", &gpsYr, &leap, &nentriesE);
         edat->leap = leap; 
	 fscanf(fp2,"%d %d %d\n", &gpsYr, &leap, &nentriesS);
    

/*first column in earth.dat or sun.dat is gps time--one long integer
  giving the number of secs that have ticked since start of GPS epoch
  +  on 1980 Jan. 6 00:00:00 UTC 
  right now it's hard-wired in that entries in earth.dat and sun.dat
  are exactly one hour apart
*/
       for (j=0; j < nentriesE; ++j){

	 /* check return value of fscanf */
	 fscanf(fp1,"%le %le %le %le %le %le %le %le %le %le\n",
		&edat->gpsE[j],&edat->position[j][0],
		&edat->position[j][1],&edat->position[j][2],
		&edat->velocity[j][0],&edat->velocity[j][1],
		&edat->velocity[j][2],&edat->acceleration[j][0],
		&edat->acceleration[j][1],
		&edat->acceleration[j][2] );
       }

       for (j=0; j < nentriesS; ++j){
	 fscanf(fp2,"%le %le %le %le %le %le %le %le %le %le\n",
		&edat->gpsS[j],&edat->sunPos[j][0],
		&edat->sunPos[j][1],&edat->sunPos[j][2],
		&edat->sunVel[j][0],&edat->sunVel[j][1],
		&edat->sunVel[j][2],
		&edat->sunAccel[j][0],&edat->sunAccel[j][1],
		&edat->sunAccel[j][2] );
       }
       
       /* below just for testing
	  printf("gpsE[nentriesE -1], gpsS[nentriesS -1] = %25.17e 
          %25.17e \n", edat->gpsE[nentriesE -1], 
          edat->gpsS[nentriesS -1]);
	  
	  printf("gpsE[nentriesE -1], gpsS[nentriesS -1] = %25.17e 
          %25.17e \n", (*edat).gpsE[nentriesE -1], 
          (*edat).gpsS[nentriesS -1]);
	  
	  for (j= 0; j<nentriesS;++j){
	  printf("j, gpsS[j] =%d %25.17e \n", j, edat->gpsS[j]); }
       */
       
       /*curt: as test, make sure last entries for gpsE,s are 
  reasonable numbers;
  specifically, checking that they are of same order as t2000 
  --within factor e


        printf("edat->gpsE[nentriesE -1]/t2000,
         1.e0*(edat->gpsE[nentriesE -1])/t2000,
           log(1.e0*(edat->gpsE[nentriesE -1])/t2000),
           fabs( log(1.e0*(edat->gpsE[nentriesE -1])/t2000)) = %25.17e 
          %25.17e %25.17e %25.17e \n", 
         edat->gpsE[nentriesE -1]/t2000,
         1.e0*(edat->gpsE[nentriesE -1])/t2000,
           log(1.e0*(edat->gpsE[nentriesE -1])/t2000),
           fabs( log(1.e0*(edat->gpsE[nentriesE -1])/t2000))
          );
*/
        ASSERT (fabs(log(1.e0*(edat->gpsE[nentriesE -1])/t2000)) < 1.e0, 
         stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);

        ASSERT (fabs(log(1.e0*(edat->gpsS[nentriesS -1])/t2000)) < 1.e0, 
         stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
	
	fclose(fp1); 
        fclose(fp2);       
	
	/*curt: is below the right return???*/
	DETATCHSTATUSPTR(stat);
	RETURN(stat);
}

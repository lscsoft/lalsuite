/************************************************************************************/
/* These function takes the GPS values of periapse passage time and the observation */
/* start time and the orbital period and returns a new value for the periapse       */
/* passage time.  The new time is that which is the closest occurance of the        */
/* periapse to the start time.  This is done so that the same values for the time   */
/* of periapse passage can be used for two different searches using different data  */
/* streches and/or detectors.  The seconds function does the reverse of this        */
/* procedure.  So given an original range of periapse passage times and the orbital */
/* period it converts a given periapse passage time to within this range.           */
/*                                                                                  */
/*			           C. Messenger                                     */
/*                                                                                  */
/*                         BIRMINGHAM UNIVERISTY -  2005                            */
/************************************************************************************/
#include <unistd.h>
#include <sys/types.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <unistd.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALDatatypes.h>
#include <lal/Date.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/LALStdlib.h>

int PeriapseShift(LIGOTimeGPS, LIGOTimeGPS *,LIGOTimeGPS, REAL8,INT4 *);
int PeriapseShiftBack(LIGOTimeGPS, LIGOTimeGPS, LIGOTimeGPS,LIGOTimeGPS *, REAL8,INT4 *);

int PeriapseShift(LIGOTimeGPS Tpold, LIGOTimeGPS *Tp,LIGOTimeGPS Tstart, REAL8 Period,INT4 *NORB)
{

  static LALStatus status;
  REAL8 tdiff;

  /*printf("recieved the following in PeriapseShift\n");
  printf("Tpold = %d %d\n",Tpold.gpsSeconds,Tpold.gpsNanoSeconds);
  printf("Tstart = %d %d\n",Tstart.gpsSeconds,Tstart.gpsNanoSeconds);
  printf("Period = %lf\n",Period);*/
  

  /* Do some error checking first */
  if ((*Tp).gpsSeconds<0) {
    fprintf(stderr,"ERROR : Negative GPS seconds component in Periapse passage time !!\n");
    exit(1);
  }
  if ((*Tp).gpsNanoSeconds<0) {
    fprintf(stderr,"ERROR : Negative GPS nanoseconds component in Periapse passage time !!\n");
    exit(1);
  }
  if ((*Tp).gpsNanoSeconds>999999999) {
    fprintf(stderr,"ERROR : GPS nanoseconds component > 10^9 in Periapse passage time !!\n");
    exit(1);
  }
  if (Tstart.gpsSeconds<0) {
    fprintf(stderr,"ERROR : Negative GPS seconds component in observation start time !!\n");
    exit(1);
  }
  if (Tstart.gpsNanoSeconds<0) {
    fprintf(stderr,"ERROR : Negative GPS nanoseconds component in observation start time !!\n");
    exit(1);
  }
  if (Tstart.gpsNanoSeconds>999999999) {
    fprintf(stderr,"ERROR : GPS nanoseconds component > 10^9 in observation start time !!\n");
    exit(1);
  }
  if (Period<=0.0) {
    fprintf(stderr,"ERROR : Negative value of orbital period !!\n");
    exit(1);
  }

  /* calculate the number of periods between periapse and start */
  LALDeltaFloatGPS(&status,&tdiff,&Tstart,&Tpold);
  
  /* deal with before periapse and after periapse differently */
  if (tdiff>=0.0) (*NORB)=(INT4)floor(tdiff/Period);
  else if (tdiff<0.0) (*NORB)=(INT4)ceil(tdiff/Period);
  
  /* reset tdiff to be a multiple of periods */
  tdiff=(REAL8)(*NORB)*Period;
  
  /* shift the periapse time so that 0 < Tstart-Tp < Period */
  LALAddFloatToGPS(&status,Tp,&Tpold,tdiff);  

  return 0;

}

/**********************************************************************************/

int PeriapseShiftBack(LIGOTimeGPS TpMIN, LIGOTimeGPS TpMAX, LIGOTimeGPS TpIN,LIGOTimeGPS *TpOUT, REAL8 Period,INT4 *NORB)
{

  static LALStatus status;
  REAL8 diff;
  REAL8 diffMIN;
  REAL8 diffMAX;
  REAL8 shift;
  LALGPSCompareResult comp;
  LALGPSCompareResult compMIN;
  LALGPSCompareResult compMAX;
  LALTimeInterval interval;
  INT4 n;

  /*printf("In PeriapseShiftBack now\n");
  printf("TpMIN is %d %d\n",TpMIN.gpsSeconds,TpMIN.gpsNanoSeconds);
  printf("TpMAX is %d %d\n",TpMAX.gpsSeconds,TpMAX.gpsNanoSeconds);
  printf("TpIN is %d %d\n",TpIN.gpsSeconds,TpIN.gpsNanoSeconds);
  printf("Period is %lf\n",Period);*/
  

  /* Do some error checking first */
  if (TpMIN.gpsSeconds<0) {
    fprintf(stderr,"ERROR : Negative GPS seconds component in minimum Periapse passage time !!\n");
    exit(1);
  }
  if (TpMIN.gpsNanoSeconds<0) {
    fprintf(stderr,"ERROR : Negative GPS nanoseconds component in minimum Periapse passage time !!\n");
    exit(1);
  }
  if (TpMIN.gpsNanoSeconds>999999999) {
    fprintf(stderr,"ERROR : GPS nanoseconds component > 10^9 in minimum Periapse passage time !!\n");
    exit(1);
  }
  if (TpMAX.gpsSeconds<0) {
    fprintf(stderr,"ERROR : Negative GPS seconds component in maximum Periapse passage time !!\n");
    exit(1);
  }
  if (TpMAX.gpsNanoSeconds<0) {
    fprintf(stderr,"ERROR : Negative GPS nanoseconds component in maximum Periapse passage time !!\n");
    exit(1);
  }
  if (TpMAX.gpsNanoSeconds>999999999) {
    fprintf(stderr,"ERROR : GPS nanoseconds component > 10^9 in maximum Periapse passage time !!\n");
    exit(1);
  }
  if (TpIN.gpsSeconds<0) {
    fprintf(stderr,"ERROR : Negative GPS seconds component in input Periapse passage time !!\n");
    exit(1);
  }
  if (TpIN.gpsNanoSeconds<0) {
    fprintf(stderr,"ERROR : Negative GPS nanoseconds component in input Periapse passage time !!\n");
    exit(1);
  }
  if (TpIN.gpsNanoSeconds>999999999) {
    fprintf(stderr,"ERROR : GPS nanoseconds component > 10^9 in input Periapse passage time !!\n");
    exit(1);
  }
  LALCompareGPS(&status,&comp,&TpMIN,&TpMAX);
  if (comp==LALGPS_LATER) {
    fprintf(stderr,"ERROR : Minimum GPS periapse passage time > Maximum GPS periapse passage time !!\n");
    exit(1);
  }
  LALDeltaFloatGPS(&status,&diff,&TpMIN,&TpMAX);
  if (fabs(diff)>=Period) {
    fprintf(stderr,"ERROR : Range in possible periapse passage times > orbital period !!\n");
    exit(1);
  }
  if (Period<=0.0) {
    fprintf(stderr,"ERROR : Negative value of orbital period !!\n");
    exit(1);
  }

  /* calculate the integer number of periods between periapse IN and a point in range */
  LALDeltaFloatGPS(&status,&diffMIN,&TpIN,&TpMIN);
  LALDeltaFloatGPS(&status,&diffMAX,&TpIN,&TpMAX);
  LALCompareGPS(&status,&compMIN,&TpMIN,&TpIN);
  LALCompareGPS(&status,&compMAX,&TpMAX,&TpIN);

  /* if we have set a non-null value of NORB */
  if (NORB!=NULL) {
    /* if NORB is positive */
    if ((*NORB)>=0) {
      shift=(REAL8)(*NORB)*Period;
      LALFloatToInterval(&status,&interval,&shift);
      LALDecrementGPS(&status,TpOUT,&TpIN,&interval);
      return 0;
    }
    /* if NORB is negative */
    else if ((*NORB)<0) {
      shift=(-1.0)*(REAL8)(*NORB)*Period;
      LALFloatToInterval(&status,&interval,&shift);
      LALIncrementGPS(&status,TpOUT,&TpIN,&interval);
      return 0;
      } 
    else {
      printf("error here with NORB variable\n");
      exit(1);
    }
  }

  /* if the input time lies between the min and max times then output with no change */
  if ((compMIN==LALGPS_EARLIER)&&(compMAX==LALGPS_LATER)) {
    (*TpOUT).gpsSeconds=TpIN.gpsSeconds;
    (*TpOUT).gpsNanoSeconds=TpIN.gpsNanoSeconds;
    return 0;
  }

  /* if both min and max times are earlier than the input time */
  else if ((compMIN==LALGPS_EARLIER)&&(compMAX==LALGPS_EARLIER)) {
    n=(INT4)floor(diffMIN/Period);
    shift=(REAL8)n*Period;
    LALFloatToInterval(&status,&interval,&shift);
    LALDecrementGPS(&status,TpOUT,&TpIN,&interval);
    return 0;
  }

  /* if both min and max times are later than the input time */
  else if ((compMIN==LALGPS_LATER)&&(compMAX==LALGPS_LATER)) {
    n=(INT4)floor(fabs(diffMAX)/Period);
    shift=(REAL8)n*Period;
    LALFloatToInterval(&status,&interval,&shift);
    LALIncrementGPS(&status,TpOUT,&TpIN,&interval);
    return 0;
  }

  else {
    fprintf(stderr,"ERROR : Something has gone wrong !!\n");
    exit(1);
  }

}

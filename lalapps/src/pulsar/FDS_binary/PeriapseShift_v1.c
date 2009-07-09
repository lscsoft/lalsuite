/*
*  Copyright (C) 2007 Chris Messenger
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
int PeriapseShiftBack(LIGOTimeGPS, LIGOTimeGPS,LIGOTimeGPS,LIGOTimeGPS *, REAL8,INT4);

int PeriapseShift(LIGOTimeGPS Tpold, LIGOTimeGPS *Tp,LIGOTimeGPS Tstart, REAL8 Period,INT4 *NORB)
{
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

  /*printf("Tpold is %d %d\n",Tpold.gpsSeconds,Tpold.gpsNanoSeconds);*/

  /* calculate the number of periods between periapse and start */
  tdiff = XLALGPSDiff(&Tstart,&Tpold);

  /*printf("tdiff is %f\n",tdiff);*/

  /* deal with before periapse and after periapse differently */
  if (tdiff>=0.0) (*NORB)=(INT4)floor(tdiff/Period);
  else if (tdiff<0.0) (*NORB)=(INT4)floor(tdiff/Period);
  
  /*printf("NORB is %d\n",(*NORB));*/

  /* reset tdiff to be a multiple of periods */
  tdiff=(REAL8)(*NORB)*Period;
  
  /*printf("tdiff now %f\n",tdiff);*/

  /* shift the periapse time so that 0 < Tstart-Tp < Period */
  *Tp = Tpold;
  XLALGPSAdd(Tp, tdiff);

  /*printf("returning Tpnew as %d %d\n",Tp->gpsSeconds,Tp->gpsNanoSeconds);*/

  return 0;

}

/**********************************************************************************/

int PeriapseShiftBack(LIGOTimeGPS Tpstart, LIGOTimeGPS Tp0, LIGOTimeGPS TpIN,LIGOTimeGPS *TpOUT, REAL8 Period,INT4 NORB)
{

  static LALStatus status;
  REAL8 shift;
  REAL8 twoPeriod;
  LALGPSCompareResult compSHIFT,compSHIFTone,compSHIFTtwo,compSHIFTthree;
  LALGPSCompareResult compone,comptwo,compthree,compfour;
  LALTimeInterval interval;
  LIGOTimeGPS tempone,temptwo,tempthree,tempfour;
  LALTimeInterval intervalone,intervaltwo;
  REAL8 halfperiod;
  LIGOTimeGPS tempOUT;

  /*printf("In PeriapseShiftBack now\n");
  printf("Tpstart is %d %d\n",Tpstart.gpsSeconds,Tpstart.gpsNanoSeconds);
  printf("TpIN is %d %d\n",TpIN.gpsSeconds,TpIN.gpsNanoSeconds);
  printf("Tp0 is %d %d\n",Tp0.gpsSeconds,Tp0.gpsNanoSeconds);
  printf("Period is %f\n",Period);*/
  

  /* Do some error checking first */
  if (Tpstart.gpsSeconds<0) {
    fprintf(stderr,"ERROR : Negative GPS seconds component in start time !!\n");
    exit(1);
  }
  if (Tpstart.gpsNanoSeconds<0) {
    fprintf(stderr,"ERROR : Negative GPS nanoseconds component in start time !!\n");
    exit(1);
  }
  if (Tp0.gpsSeconds<0) {
    fprintf(stderr,"ERROR : Negative GPS seconds component in Tp0 !!\n");
    exit(1);
  }
  if (Tp0.gpsNanoSeconds<0) {
    fprintf(stderr,"ERROR : Negative GPS nanoseconds component in Tp0 !!\n");
    exit(1);
  }
  if (Tp0.gpsNanoSeconds>999999999) {
    fprintf(stderr,"ERROR : GPS nanoseconds component > 10^9 in Tp0 !!\n");
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
  if (Period<=0.0) {
    fprintf(stderr,"ERROR : Negative value of orbital period !!\n");
    exit(1);
  }
   
 
  /* check whether the input Tp has been shifted out of range of 1 period before the start time */
  
  /* calculate iterval equal to plus one period */
  LALFloatToInterval(&status,&interval,&Period);

  /* calculate interval equal to plus two periods */
  twoPeriod=2.0*Period;
  LALFloatToInterval(&status,&intervaltwo,&twoPeriod);
  
  /* minus a single period from the start time */
  LALDecrementGPS(&status,&tempone,&Tpstart,&interval);
  
  /* add a single period to the start time */
  LALIncrementGPS(&status,&temptwo,&Tpstart,&interval);
  
  /* add two periods from to the start time */
  LALDecrementGPS(&status,&tempthree,&Tpstart,&intervaltwo);
  
  /* compare the input time to the start time */
  LALCompareGPS(&status,&compSHIFT,&Tpstart,&TpIN);
  /*printf("start : input -> %d\n",compSHIFT);*/

  /* compare the input time to the start-period */
  LALCompareGPS(&status,&compSHIFTone,&tempone,&TpIN);
  /*printf("start-period : input -> %d\n",compSHIFTone);*/

  /* compare the input time to the start+period */
  LALCompareGPS(&status,&compSHIFTtwo,&temptwo,&TpIN);
  /*printf("start+period : input -> %d\n",compSHIFTtwo);*/

  /* compare the input time to the start-2*period */
  LALCompareGPS(&status,&compSHIFTthree,&tempthree,&TpIN);
  /*printf("start-2period : input -> %d\n",compSHIFTthree);*/

  /* if input time is after the start time by less than one orbit */
  if ((compSHIFT==LALGPS_EARLIER)&&(compSHIFTtwo==LALGPS_LATER)) {
    printf("input time is after the start time and less than starttime+period -> NORB=NORB+1\n");
    NORB=NORB+1;
  }

  /* if input time is before the start time by more than one orbit */
  else if ((compSHIFTone==LALGPS_LATER)&&(compSHIFTthree==LALGPS_EARLIER)) {
    printf("input time is before the start time-period and after the starttime-2period -> NORB=NORB-1\n");
    NORB=NORB-1;
  }

  /* if input time is after the start time by more than one orbit */
  else if (compSHIFTtwo==LALGPS_EARLIER) {
    fprintf(stderr,"ERROR in PeripaseShiftBack : input time is after the start time by more than one orbital period !!!\n");
    exit(1);
  }

  /* if input time is before the start time by more than two orbits */
  else if (compSHIFTthree==LALGPS_LATER) {
    fprintf(stderr,"ERROR in PeripaseShiftBack : input time is before the start time by more than two orbital periods !!!\n");
    exit(1);
  }
    
  /* if NORB is positive */
  if (NORB>=0) {
    shift=(REAL8)(NORB)*Period;
    LALFloatToInterval(&status,&interval,&shift);
    LALDecrementGPS(&status,TpOUT,&TpIN,&interval);
    /*printf("in shift back : NORB = %d returning %d %d\n",NORB,TpOUT->gpsSeconds,TpOUT->gpsNanoSeconds);*/
  }
  /* if NORB is negative */
  else if (NORB<0) {
    shift=(-1.0)*(REAL8)(NORB)*Period;
    LALFloatToInterval(&status,&interval,&shift);
    LALIncrementGPS(&status,TpOUT,&TpIN,&interval);
    /*printf("in shift back : NORB = %d returning %d %d\n",NORB,TpOUT->gpsSeconds,TpOUT->gpsNanoSeconds);*/
  } 
  else {
    printf("error here with NORB variable\n");
    exit(1);
  }

  /* test if TpOUT is with +/- 1/2 period from Tp0 */
  
  /* calculate iterval equal to plus half period */
  halfperiod=Period/2.0;
  LALFloatToInterval(&status,&intervalone,&halfperiod);
  LALFloatToInterval(&status,&intervaltwo,&Period);

  /* make Tp0 +/- half period */
  LALDecrementGPS(&status,&tempone,&Tp0,&intervalone);
  LALIncrementGPS(&status,&temptwo,&Tp0,&intervalone);
  LALDecrementGPS(&status,&tempthree,&Tp0,&intervaltwo);
  LALIncrementGPS(&status,&tempfour,&Tp0,&intervaltwo);

  /* test if in the right range */
  LALCompareGPS(&status,&compone,TpOUT,&tempone);
  LALCompareGPS(&status,&comptwo,TpOUT,&temptwo);
  LALCompareGPS(&status,&compthree,TpOUT,&tempthree);
  LALCompareGPS(&status,&compfour,TpOUT,&tempfour);

  if ((compone==LALGPS_LATER)&&(comptwo==LALGPS_EARLIER)) return 0;
  else if ((compone==LALGPS_EARLIER)&&(compthree==LALGPS_LATER)) {
    tempOUT.gpsSeconds=TpOUT->gpsSeconds;
    tempOUT.gpsNanoSeconds=TpOUT->gpsNanoSeconds;
    LALIncrementGPS(&status,TpOUT,&tempOUT,&intervaltwo);
    return 0;
  }
  else if ((comptwo==LALGPS_LATER)&&(compfour==LALGPS_EARLIER)) {
    tempOUT.gpsSeconds=TpOUT->gpsSeconds;
    tempOUT.gpsNanoSeconds=TpOUT->gpsNanoSeconds;
    LALDecrementGPS(&status,TpOUT,&tempOUT,&intervaltwo);
    return 0;
  }
  else if (compthree==LALGPS_EARLIER) {
    printf("ERROR : output Tp is earlier than Tp0 by more than a period !!!\n");
    exit(1);
  }
  else if (compfour==LALGPS_LATER) {
    printf("ERROR : output Tp is later than Tp0 by more than a period !!!\n");
    exit(1);
  }

  return 0;

}

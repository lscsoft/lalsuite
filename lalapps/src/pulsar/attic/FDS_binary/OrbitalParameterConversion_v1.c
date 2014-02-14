/*
*  Copyright (C) 2007 Chris Messenger, Reinhard Prix
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
/* The functions below allow conversion between the different forms of the orbital  */
/* parameters required to describe the orbital motion of a source in a circular     */
/* orbit.  They allow conversion between the standard representation of (semi-major */
/* axis, periapse passage time) and the parameterised representation (X,Y).         */
/*                                                                                  */
/*			           C. Messenger                                     */
/*                                                                                  */
/*                         BIRMINGHAM UNIVERISTY -  2004                            */
/************************************************************************************/

#include "GenerateBinaryMesh_v1.h"

/* define the following so they are visible to the root finding routine */
static REAL8 sma_GLOBAL;
static REAL8 period_GLOBAL;
LALStatus empty_status;

static void OrbPhaseFunc(LALStatus *, REAL8 *, REAL8, void *);

int ConvertXYtoRTperi(XYLocation *XYloc, RTPLocation *RTPloc)
{

  /* this routine takes a location in XY space and converts it to a location in RT space */   
  REAL8 alpha;
  REAL8 Torb;
  REAL8 Tlight;
  REAL8 deltaT;
  REAL8 sma;
  REAL8 period;
  LIGOTimeGPS tstartSSB;
  REAL8 X;
  REAL8 Y;
  REAL8 ecc;
  REAL8 argp;
  LIGOTimeGPS tperi; 

  /* read in the XY location structure */
  if (XYloc->ecc!=0.0) {
    fprintf(stderr,"WARNING : Non-zero eccentricity being used in a function for circular orbits !!\n");
  }
  X=XYloc->X;
  Y=XYloc->Y;
  tstartSSB.gpsSeconds=XYloc->tstartSSB.gpsSeconds;
  tstartSSB.gpsNanoSeconds=XYloc->tstartSSB.gpsNanoSeconds;
  period=XYloc->period;
  ecc=XYloc->ecc;
  argp=XYloc->argp;

  /* convert to R and alpha easily */
  sma=sqrt(X*X+Y*Y);
  alpha=atan2(Y,X);
  /* keep the range of alpha as 0 -> 2pi */
  /* if (alpha<0.0) alpha=alpha+LAL_TWOPI; */ 
  /*printf("alpha is %f\n",alpha);*/

  /* Convert alpha and R to time since periapse as measured in SSB */
  Torb=alpha*(period)/LAL_TWOPI;
  Tlight=sma*sin(alpha);
  deltaT=Torb+Tlight;
  /*printf("deltaT is %f\n",deltaT);*/
  
  if (alpha>=0.0) {
    tperi = tstartSSB;
    XLALGPSAdd(&tperi, -deltaT);
  }
  else if (alpha<0.0) {
    deltaT=(-1.0)*deltaT;
    tperi = tstartSSB;
    XLALGPSAdd(&tperi, deltaT);
  }
  else {
    printf("error with alpha is parameter conversion\n");
    exit(1);
  }

  /* fill in RTPLocation structure */
  RTPloc->sma=sma;
  RTPloc->period=period;
  RTPloc->tperi.gpsSeconds=tperi.gpsSeconds;
  RTPloc->tperi.gpsNanoSeconds=tperi.gpsNanoSeconds;
  RTPloc->ecc=ecc;
  RTPloc->argp=argp;


  return 0;

}

/*******************************************************************************/

int ConvertRTperitoXY(RTPLocation *RTPloc, XYLocation *XYloc, REAL8 *alpha)
{

  /* this routine takes a point in RT space and converts it to XY space */
  LALStatus status = empty_status;
  DFindRootIn input;
  REAL8 deltaTorb;
  REAL8 alphatemp;
  REAL8 Xtemp;
  REAL8 Ytemp;
  LIGOTimeGPS tstartSSB;
  REAL8 ecc;
  REAL8 argp;
  LIGOTimeGPS tperi; 

  /* read in the XY location structure */
  if (RTPloc->ecc!=0.0) {
    fprintf(stderr,"WARNING : Non-zero eccentricity being used in a function for circular orbits !!\n");
  }
  sma_GLOBAL=RTPloc->sma;
  period_GLOBAL=RTPloc->period;
  tperi.gpsSeconds=RTPloc->tperi.gpsSeconds;
  tperi.gpsNanoSeconds=RTPloc->tperi.gpsNanoSeconds;
  tstartSSB.gpsSeconds=RTPloc->tstartSSB.gpsSeconds;
  tstartSSB.gpsNanoSeconds=RTPloc->tstartSSB.gpsNanoSeconds;
  ecc=RTPloc->ecc;
  argp=RTPloc->argp;

  /* calculate difference between the times using LAL functions */
  deltaTorb = XLALGPSDiff(&tstartSSB,&tperi);

  /* begin root finding procedure */
  input.function = OrbPhaseFunc;
  input.xmin = 0.0;
  input.xmax = LAL_TWOPI;
  input.xacc = 1E-12;
  /* expand domain until a root is bracketed */
  LALDBracketRoot(&status,&input,&deltaTorb); 
  /* bisect domain to find initial phase at observation start */
  LALDBisectionFindRoot(&status,&alphatemp,&input,&deltaTorb); 

  Xtemp=(sma_GLOBAL)*cos(alphatemp);
  Ytemp=(sma_GLOBAL)*sin(alphatemp);

  /* fill up output XY location structure */
  XYloc->X=Xtemp;
  XYloc->Y=Ytemp;
  XYloc->period=period_GLOBAL;
  XYloc->tstartSSB.gpsSeconds=tstartSSB.gpsSeconds;
  XYloc->tstartSSB.gpsNanoSeconds=tstartSSB.gpsNanoSeconds;
  XYloc->ecc=ecc;
  XYloc->argp=argp;

  *alpha=alphatemp;

  return 0;

}

/*******************************************************************************/

static void OrbPhaseFunc(LALStatus *status, REAL8 *y, REAL8 alpha, void *y_0)
{
  INITSTATUS(status);
  ASSERT(y_0,status, 1, "Null pointer");
  /* this is the transendental function we need to solve to find the true initial phase */
  /* note that it also includes a retarded time delay */
  *y = *(REAL8 *)y_0*(-1.0) + alpha*(period_GLOBAL/LAL_TWOPI)+sma_GLOBAL*sin(alpha);
  RETURN(status);
}

/*******************************************************************************/


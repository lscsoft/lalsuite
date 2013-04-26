/*
 *  xmlToLALInference.c:  translation tool between SimInspiral table 
 *  in a LIGOLW xml file and lalinference_mcmc command line.
 *
 *  Copyright (C) 2012 Vivien Raymond
 *
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
 *
 **********************************************************************
 *
 * Usage:
 *     ./xmlToLALInference <xml file> <injection number (0 is first)>
 *
 */


#include <ctype.h>
#include <getopt.h>
//#include <lalapps.h>
#include <lal/Date.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LALConstants.h>
#include <lal/LIGOMetadataInspiralUtils.h>
//#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/LIGOLwXML.h>
#include <lal/Random.h>
#include <lal/AVFactories.h>
#include <lal/InspiralInjectionParams.h>
#include <lal/LALDetectors.h>
#include <lal/LALSimulation.h>
//#include <processtable.h>
#include <lal/RingUtils.h>
//#include <LALAppsVCSInfo.h>

//#include "inspiral.h"


int main( int argc, char *argv[] )
{
  UINT4 Ninj=0, i=0;
	SimInspiralTable *injTable=NULL;
  REAL8 a_spin1, a_spin2, theta_spin1, theta_spin2, phi_spin1, phi_spin2;
  double GPSdouble, gmst, ra, psi, q;
  LIGOTimeGPS GPSlal;

  
  Ninj=SimInspiralTableFromLIGOLw(&injTable,argv[1],0,0);
  
  if(Ninj<atoi(argv[2])) fprintf(stderr,"Error reading event %d from %s\n",atoi(argv[2]),argv[1]);
  
  while(i<atoi(argv[2])) {i++; injTable = injTable->next;} /* Select event */
  
  q = injTable->mass2 / injTable->mass1;
  
  if(q>1.0) q=1.0/q;
  
  a_spin1 = sqrt(injTable->spin1x*injTable->spin1x+injTable->spin1y*injTable->spin1y+injTable->spin1z*injTable->spin1z);
  a_spin2 = sqrt(injTable->spin2x*injTable->spin2x+injTable->spin2y*injTable->spin2y+injTable->spin2z*injTable->spin2z);
  
  if(a_spin1 == 0.0){
    theta_spin1 = 0.000000000001;
    phi_spin1 = 0.0;
  }else{
    theta_spin1 = acos(injTable->spin1z / a_spin1);
    phi_spin1 = atan2(injTable->spin1y,injTable->spin1x);
    if (phi_spin1 < 0.0) phi_spin1+=2*LAL_PI;
  }
  if(a_spin2 == 0.0){
    theta_spin2 = 0.000000000001;
    phi_spin2 = 0.0;
  }else{
    theta_spin2 = acos(injTable->spin2z / a_spin2);
    phi_spin2 = atan2(injTable->spin2y,injTable->spin2x);
    if (phi_spin2 < 0.0) phi_spin2+=2*LAL_PI;
  }
  
  if(injTable->polarization>=LAL_PI){
    psi = injTable->polarization - LAL_PI;
  }else{
    psi = injTable->polarization;
  }
  printf("--phi2 %.12lf --theta2 %.12lf --a2 %.12lf --phi1 %.12lf --theta1 %.12lf --a1 %.12lf --iota %.12lf --psi %.12lf --dec %.12lf --ra %.12lf --dist %.12lf --phi %.12lf --time %.12lf --q %.12lf --mc %.12lf\n", phi_spin2, theta_spin2, a_spin2, phi_spin1, theta_spin1, a_spin1, injTable->inclination, psi, injTable->latitude, injTable->longitude, injTable->distance, injTable->coa_phase, XLALGPSGetREAL8(&(injTable->geocent_end_time)) , q, injTable->mchirp);

}

/*
 * Copyright (C) 2006 C. Messenger
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

/*********************************************************************************/
/**
 * \author C. Messenger
 * \file
 * \brief
 * A collection of commonly used utilities in the SideBand directory
 *
 */

/* System includes */
#include <stdio.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#include <time.h>
#endif

/* LAL-includes */
#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/AVFactories.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lalapps.h>
#include "SideBand.h"

/*---------- DEFINES ----------*/
#define TRUE (1==1)
#define FALSE (1==0)
#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10
#define MINX 1e-3

/*----- Error-codes -----*/
#define SIDEBANDUTILSC_ENULL 	        	1
#define SIDEBANDUTILSC_ESYS     		2
#define SIDEBANDUTILSC_EINPUT   		3
#define SIDEBANDUTILSC_EXLAL		        4

#define SIDEBANDUTILSC_MSGENULL 		"Arguments contained an unexpected null pointer"
#define SIDEBANDUTILSC_MSGESYS	        	"System call failed (probably file IO)"
#define SIDEBANDUTILSC_MSGEINPUT   		"Invalid input"
#define SIDEBANDUTILSC_MSGEXLAL	        	"XLALFunction-call failed"

/*---------- Global variables ----------*/
extern int vrbflg;		/**< defined in lalapps.c */

/****************************************************************************************/
/****************************************************************************************/
/* Read in the time stamps file */
/****************************************************************************************/
/****************************************************************************************/
void ReadTimeStamps(LALStatus *status,
		    CHAR *timestampsfile,
		    INT4 tsft,
		    SideBandTemplateParams **TParams)
{
  
  FILE *fp = NULL;
  INT4 i,k;
  INT4 N;
  fpos_t pos;
  INT4 dum;
  INT4 Ngap;
 

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  
  ASSERT ( *TParams, status, SIDEBANDUTILSC_ENULL,SIDEBANDUTILSC_MSGENULL );
  
  /* open and read parameter range file */
  if ((fp = fopen(timestampsfile,"r"))==NULL) {
    XLALPrintError("\nError opening file '%s' for reading..\n\n",timestampsfile);
    /* return (READTIMESTAMPSC_ESYS); */
  }
  
  fgetpos(fp,&pos);  /* record current position in file */
  /* count the number of lines in the file */
  i = 0;
  while (fscanf(fp,"%d %d",&dum,&dum)!=EOF) i++;
  N = i;

  /* allocate memory */
  (*TParams)->timestamps = XLALCreateINT4Vector(N);
 
  /* read again into memory */
  fsetpos(fp,&pos);
  i = 0;
  while (fscanf(fp,"%d %d",&((*TParams)->timestamps->data[i]),&dum)!=EOF) i++;

  fclose(fp);

  /* determine number of gaps */
  Ngap = 0;
  for (i=0;i<N-2;i++) {
    if ((*TParams)->timestamps->data[i]+tsft!=(*TParams)->timestamps->data[i+1]) Ngap++;
  }
  printf("Ngap = %d\n",Ngap);

  /* allocate memory to gap vector */
  (*TParams)->gapvectorstart = XLALCreateINT4Vector(Ngap+1);
  (*TParams)->gapvectorend = XLALCreateINT4Vector(Ngap+1);
  
  /* fill in gap vector */
  (*TParams)->gapvectorstart->data[0] = (*TParams)->timestamps->data[0];
  k = 0;
  for (i=1;i<N-2;i++) {
    if ((*TParams)->timestamps->data[i]+tsft!=(*TParams)->timestamps->data[i+1]) {
      (*TParams)->gapvectorend->data[k] = (*TParams)->timestamps->data[i]+tsft;
      (*TParams)->gapvectorstart->data[k+1] = (*TParams)->timestamps->data[i+1];
      k++;
    }
  }
  (*TParams)->gapvectorend->data[Ngap] = (*TParams)->timestamps->data[N-1]+tsft;

  /* compute tobs */
  (*TParams)->Tobs = 0.0;
  for (i=0;i<Ngap+1;i++) {
    (*TParams)->Tobs += (*TParams)->gapvectorend->data[i] - (*TParams)->gapvectorstart->data[i];
  }
  printf("tobs = %f\n",(*TParams)->Tobs);

  /* fill in other Tparams */
  (*TParams)->T = (*TParams)->timestamps->data[N-1]+tsft - (*TParams)->timestamps->data[0];
  printf("T = %f\n",(*TParams)->T);
  (*TParams)->tstart.gpsSeconds = (*TParams)->timestamps->data[0];
  (*TParams)->tstart.gpsNanoSeconds = 0;
  (*TParams)->nsft = N;

  DETATCHSTATUSPTR ( status );
  RETURN ( status );

} 

/****************************************************************************************/
/****************************************************************************************/
/* Compute the window function due to gappy finite data set */
/****************************************************************************************/
/****************************************************************************************/
void ComputeSideBandWindow(LALStatus *status,
			   ABCcoefficients *ABCco,
			   CHAR *outfile,
			   SideBandTemplateParams **TParams)
{

  INT4 i,j;
  FILE *fp = NULL;
  REAL8 T0;
  REAL8 W = LAL_TWOPI/LAL_DAYSID_SI;
  REAL8 w0 = ABCco->omega0;
  REAL8 x;
  UINT4 k;
  REAL8 et,st;
  COMPLEX16 ae,as,be,bs;
  REAL8 ddf = (*TParams)->dfwindow/1000;
  
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  
  ASSERT (*TParams, status, SIDEBANDUTILSC_ENULL,SIDEBANDUTILSC_MSGENULL );
 
  if ((fp = fopen(outfile,"w"))==NULL) {
    XLALPrintError("\nError opening file '%s' for writing..\n\n",outfile);
    /* return (SIDEBANDUTILSC_ESYS); */
  }
  
  printf("number of data stretches = %d\n",(*TParams)->gapvectorstart->length);
  T0 = (*TParams)->gapvectorstart->data[0];
  printf("T0 = %6.12f\n",T0);

  /* allocate memory */
  (*TParams)->wa = XLALCreateCOMPLEX16Vector(2*(*TParams)->windowrange);
  (*TParams)->wb = XLALCreateCOMPLEX16Vector(2*(*TParams)->windowrange);
  
  /* loop over highly resolved frequency bins */
  for (j=0;j<2*(*TParams)->windowrange;j++) {
      
    /* initialise */
    (*TParams)->wa->data[j].real_FIXME = 0.0;
    (*TParams)->wa->data[j].imag_FIXME = 0.0;
    (*TParams)->wb->data[j].real_FIXME = 0.0;
    (*TParams)->wb->data[j].imag_FIXME = 0.0;
    
    /* define frequency */
    x = (-1.0)*LAL_TWOPI*(REAL8)(*TParams)->windowrange*(*TParams)->dfwindow + LAL_TWOPI*(REAL8)j*(*TParams)->dfwindow;

    /* loop over each segment of data */
    for (i=0;i<(INT4)(*TParams)->gapvectorstart->length;i++) {
     
      et = (*TParams)->gapvectorend->data[i] - T0;
      st = (*TParams)->gapvectorstart->data[i] - T0;


      for (k=0;k<=2;k++) {   

	if ((fabs(k*W-x)>ddf)&&(fabs(k*W+x)>ddf)) {
	  
	  ae.real_FIXME = -0.5*(1.0/(k*W+x))*(-ABCco->apco[k]*cos(k*(w0-W*et)-x*et)+ABCco->aco[k]*sin(k*(w0-W*et)-x*et))
	    -0.5*(1.0/(k*W-x))*(-ABCco->apco[k]*cos(k*(w0-W*et)+x*et)+ABCco->aco[k]*sin(k*(w0-W*et)+x*et));
	  ae.imag_FIXME = -0.5*(1.0/(k*W+x))*(-ABCco->aco[k]*cos(k*(w0-W*et)-x*et)-ABCco->apco[k]*sin(k*(w0-W*et)-x*et))
	    -0.5*(1.0/(k*W-x))*(ABCco->aco[k]*cos(k*(w0-W*et)+x*et)+ABCco->apco[k]*sin(k*(w0-W*et)+x*et));	    
	  as.real_FIXME = -0.5*(1.0/(k*W+x))*(-ABCco->apco[k]*cos(k*(w0-W*st)-x*st)+ABCco->aco[k]*sin(k*(w0-W*st)-x*st))
	    -0.5*(1.0/(k*W-x))*(-ABCco->apco[k]*cos(k*(w0-W*st)+x*st)+ABCco->aco[k]*sin(k*(w0-W*st)+x*st));
	  as.imag_FIXME =-0.5*(1.0/(k*W+x))*(-ABCco->aco[k]*cos(k*(w0-W*st)-x*st)-ABCco->apco[k]*sin(k*(w0-W*st)-x*st))
	    -0.5*(1.0/(k*W-x))*(ABCco->aco[k]*cos(k*(w0-W*st)+x*st)+ABCco->apco[k]*sin(k*(w0-W*st)+x*st));
	  
	  be.real_FIXME = -0.5*(1.0/(k*W+x))*(-ABCco->bpco[k]*cos(k*(w0-W*et)-x*et)+ABCco->bco[k]*sin(k*(w0-W*et)-x*et))
	    -0.5*(1.0/(k*W-x))*(-ABCco->bpco[k]*cos(k*(w0-W*et)+x*et)+ABCco->bco[k]*sin(k*(w0-W*et)+x*et));
	  be.imag_FIXME = -0.5*(1.0/(k*W+x))*(-ABCco->bco[k]*cos(k*(w0-W*et)-x*et)-ABCco->bpco[k]*sin(k*(w0-W*et)-x*et))
	    -0.5*(1.0/(k*W-x))*(ABCco->bco[k]*cos(k*(w0-W*et)+x*et)+ABCco->bpco[k]*sin(k*(w0-W*et)+x*et));	    
	  bs.real_FIXME = -0.5*(1.0/(k*W+x))*(-ABCco->bpco[k]*cos(k*(w0-W*st)-x*st)+ABCco->bco[k]*sin(k*(w0-W*st)-x*st))
	    -0.5*(1.0/(k*W-x))*(-ABCco->bpco[k]*cos(k*(w0-W*st)+x*st)+ABCco->bco[k]*sin(k*(w0-W*st)+x*st));
	  bs.imag_FIXME = -0.5*(1.0/(k*W+x))*(-ABCco->bco[k]*cos(k*(w0-W*st)-x*st)-ABCco->bpco[k]*sin(k*(w0-W*st)-x*st))
	    -0.5*(1.0/(k*W-x))*(ABCco->bco[k]*cos(k*(w0-W*st)+x*st)+ABCco->bpco[k]*sin(k*(w0-W*st)+x*st));
	  
	}
	else if ((fabs(k*W-x)<ddf)&&(fabs(k*W+x)>ddf)) {
	 
	  
	  ae.real_FIXME = -0.5*(1.0/(k*W+x))*(-ABCco->apco[k]*cos(k*(w0-W*et)-x*et)+ABCco->aco[k]*sin(k*(w0-W*et)-x*et))
	    -0.5*et*(-ABCco->apco[k]*sin(k*w0)-ABCco->aco[k]*cos(k*w0));
	  ae.imag_FIXME = -0.5*(1.0/(k*W+x))*(-ABCco->aco[k]*cos(k*(w0-W*et)-x*et)-ABCco->apco[k]*sin(k*(w0-W*et)-x*et))
	    -0.5*et*(ABCco->aco[k]*sin(k*w0)-ABCco->apco[k]*cos(k*w0));
	  as.real_FIXME =-0.5*(1.0/(k*W+x))*(-ABCco->apco[k]*cos(k*(w0-W*st)-x*st)+ABCco->aco[k]*sin(k*(w0-W*st)-x*st))
	    -0.5*st*(-ABCco->apco[k]*sin(k*w0)-ABCco->aco[k]*cos(k*w0));
	  as.imag_FIXME = -0.5*(1.0/(k*W+x))*(-ABCco->aco[k]*cos(k*(w0-W*st)-x*st)-ABCco->apco[k]*sin(k*(w0-W*st)-x*st))
	    -0.5*st*(ABCco->aco[k]*sin(k*w0)-ABCco->apco[k]*cos(k*w0));
	  
	  be.real_FIXME = -0.5*(1.0/(k*W+x))*(-ABCco->bpco[k]*cos(k*(w0-W*et)-x*et)+ABCco->bco[k]*sin(k*(w0-W*et)-x*et))
	    -0.5*et*(-ABCco->bpco[k]*sin(k*w0)-ABCco->bco[k]*cos(k*w0));
	  be.imag_FIXME = -0.5*(1.0/(k*W+x))*(-ABCco->bco[k]*cos(k*(w0-W*et)-x*et)-ABCco->bpco[k]*sin(k*(w0-W*et)-x*et))
	    -0.5*et*(ABCco->bco[k]*sin(k*w0)-ABCco->bpco[k]*cos(k*w0));
	  bs.real_FIXME =-0.5*(1.0/(k*W+x))*(-ABCco->bpco[k]*cos(k*(w0-W*st)-x*st)+ABCco->bco[k]*sin(k*(w0-W*st)-x*st))
	    -0.5*st*(-ABCco->bpco[k]*sin(k*w0)-ABCco->bco[k]*cos(k*w0));
	  bs.imag_FIXME = -0.5*(1.0/(k*W+x))*(-ABCco->bco[k]*cos(k*(w0-W*st)-x*st)-ABCco->bpco[k]*sin(k*(w0-W*st)-x*st))
	    -0.5*st*(ABCco->bco[k]*sin(k*w0)-ABCco->bpco[k]*cos(k*w0));
	  
	}
	else if ((fabs(k*W-x)>ddf)&&(fabs(k*W+x)<ddf)) {
	 
	  
	  ae.real_FIXME = -0.5*et*(-ABCco->apco[k]*sin(k*w0)-ABCco->aco[k]*cos(k*w0))
	    -0.5*(1.0/(k*W-x))*(-ABCco->apco[k]*cos(k*(w0-W*et)+x*et)+ABCco->aco[k]*sin(k*(w0-W*et)+x*et));
	  ae.imag_FIXME = -0.5*et*(-ABCco->aco[k]*sin(k*w0)+ABCco->apco[k]*cos(k*w0))
	    -0.5*(1.0/(k*W-x))*(ABCco->aco[k]*cos(k*(w0-W*et)+x*et)+ABCco->apco[k]*sin(k*(w0-W*et)+x*et));	    
	  as.real_FIXME = -0.5*st*(-ABCco->apco[k]*sin(k*w0)-ABCco->aco[k]*cos(k*w0))
	    -0.5*(1.0/(k*W-x))*(-ABCco->apco[k]*cos(k*(w0-W*st)+x*st)+ABCco->aco[k]*sin(k*(w0-W*st)+x*st));
	  as.imag_FIXME = -0.5*st*(-ABCco->aco[k]*sin(k*w0)+ABCco->apco[k]*cos(k*w0))
	    -0.5*(1.0/(k*W-x))*(ABCco->aco[k]*cos(k*(w0-W*st)+x*st)+ABCco->apco[k]*sin(k*(w0-W*st)+x*st));
	  
	  be.real_FIXME = -0.5*et*(-ABCco->bpco[k]*sin(k*w0)-ABCco->bco[k]*cos(k*w0))
	    -0.5*(1.0/(k*W-x))*(-ABCco->bpco[k]*cos(k*(w0-W*et)+x*et)+ABCco->bco[k]*sin(k*(w0-W*et)+x*et));
	  be.imag_FIXME = -0.5*et*(-ABCco->bco[k]*sin(k*w0)+ABCco->bpco[k]*cos(k*w0))
	    -0.5*(1.0/(k*W-x))*(ABCco->bco[k]*cos(k*(w0-W*et)+x*et)+ABCco->bpco[k]*sin(k*(w0-W*et)+x*et));	    
	  bs.real_FIXME = -0.5*st*(-ABCco->bpco[k]*sin(k*w0)-ABCco->bco[k]*cos(k*w0))
	    -0.5*(1.0/(k*W-x))*(-ABCco->bpco[k]*cos(k*(w0-W*st)+x*st)+ABCco->bco[k]*sin(k*(w0-W*st)+x*st));
	  bs.imag_FIXME = -0.5*st*(-ABCco->bco[k]*sin(k*w0)+ABCco->bpco[k]*cos(k*w0))
	    -0.5*(1.0/(k*W-x))*(ABCco->bco[k]*cos(k*(w0-W*st)+x*st)+ABCco->bpco[k]*sin(k*(w0-W*st)+x*st));
	  
	}
	else {
		  
	  ae.real_FIXME = -et*(-ABCco->apco[k]*sin(k*w0)-ABCco->aco[k]*cos(k*w0));
	  ae.imag_FIXME = 0.0;
	  as.real_FIXME = -st*(-ABCco->apco[k]*sin(k*w0)-ABCco->aco[k]*cos(k*w0));
	  as.imag_FIXME = 0.0;
	  
	  be.real_FIXME = -et*(-ABCco->bpco[k]*sin(k*w0)-ABCco->bco[k]*cos(k*w0));
	  be.imag_FIXME = 0.0;
	  bs.real_FIXME = -st*(-ABCco->bpco[k]*sin(k*w0)-ABCco->bco[k]*cos(k*w0));
	  bs.imag_FIXME = 0.0;
	  
	}	  
	
	(*TParams)->wa->data[j].real_FIXME += creal(ae) - creal(as);
	(*TParams)->wa->data[j].imag_FIXME += cimag(ae) - cimag(as);
	(*TParams)->wb->data[j].real_FIXME += creal(be) - creal(bs);
	(*TParams)->wb->data[j].imag_FIXME += cimag(be) - cimag(bs);
	
      }
      
    }
    
    fprintf(fp,"%6.12f %6.12f %6.12f %6.12f %6.12f\n",x,creal((*TParams)->wa->data[j]),cimag((*TParams)->wa->data[j]),creal((*TParams)->wb->data[j]),cimag((*TParams)->wb->data[j]));
    
    
  }
  printf("W = %6.12f\n",W);
  fclose(fp);
  
  
  DETATCHSTATUSPTR ( status );
  RETURN ( status );

}


/***********************************************************************************/
/***********************************************************************************/
/* Computes the likelihood associated with the model and the data given a set of   */
/* model parameters */
/***********************************************************************************/
/***********************************************************************************/
void ComputeSideBandLikelihood(LALStatus *status,
			       SideBandMCMCVector *lambda,
			       SideBandDataSet *Data,
			       SideBandTemplate **Template,
			       SideBandTemplateParams *TParams)
{
  
  REAL8 Li = 0.0;                       /* the likelihood for the imaginary part of h(f) */
  REAL8 Lr = 0.0;                       /* the likelihood for the real part of h(f) */
  BinarySourceParams BSParams;          /* structure for binary source parameters */
  INT4 i;
 
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT (lambda,status,SIDEBANDUTILSC_ENULL,SIDEBANDUTILSC_MSGENULL );
  ASSERT (*Template,status,SIDEBANDUTILSC_ENULL,SIDEBANDUTILSC_MSGENULL );
  ASSERT (TParams,status,SIDEBANDUTILSC_ENULL,SIDEBANDUTILSC_MSGENULL );
  ASSERT (Data,status,SIDEBANDUTILSC_ENULL,SIDEBANDUTILSC_MSGENULL );

  /* fill in binary source parameters */
  BSParams.OrbitalSemiMajorAxis = lambda->a;
  BSParams.OrbitalPeriod = lambda->period;
  BSParams.OrbitalEccentricity = lambda->e; 
  BSParams.ArgumentofPeriapse = lambda->argp;  
  BSParams.TimeofSSBPeriapsePassage.gpsSeconds = lambda->tp.gpsSeconds;
  BSParams.TimeofSSBPeriapsePassage.gpsNanoSeconds = lambda->tp.gpsNanoSeconds;
  BSParams.alpha = TParams->alpha;
  BSParams.delta = TParams->delta;
  BSParams.f0 = lambda->f0;
  BSParams.phi0 = lambda->phi0;
  BSParams.psi = lambda->psi;
  BSParams.cosi = lambda->cosi;
  BSParams.h0 = lambda->h0;

  if (lalDebugLevel) printf ("\nFilled in the template parameter structures.\n");

  /* call the template generating function */
  GenerateSideBandTemplate(status->statusPtr,&BSParams,TParams,Template);
  
  /* for (i=0;i<(INT4)Data->freq->length;i++) printf("%6.12f %6.12f %6.12f %6.12f %6.12f\n",
						  Data->freq->data[i],Data->fourier->data[i].re,Data->fourier->data[i].im,
						  (*Template)->fourier->data[i].re,(*Template)->fourier->data[i].im);  
  
  */

  /* cycle over the data points */
  for (i=0;i<(INT4)Data->freq->length;i++) {

    /* compute log likelihood for real and inaginary parts of h(f) */
    Lr = Lr + (-0.5)*((creal(Data->fourier->data[i])-creal((*Template)->fourier->data[i]))*(creal(Data->fourier->data[i])-creal((*Template)->fourier->data[i])));
    Li = Li + (-0.5)*((cimag(Data->fourier->data[i])-cimag((*Template)->fourier->data[i]))*(cimag(Data->fourier->data[i])-cimag((*Template)->fourier->data[i])));
   
  } 
 
  /* combine real and imaginary log likelihoods */
  lambda->logL = 2.0*(Lr + Li)/(TParams->sqrtSh*TParams->sqrtSh);
  /* printf("Lr = %6.12f Li = %6.12f reduced logL = %6.12f\n",Lr,Li,lambda->logL); */
  
  DETATCHSTATUSPTR (status);
  RETURN(status);

}

/***********************************************************************************/
/***********************************************************************************/
/** Load Ephemeris from ephemeris data-files  */
/***********************************************************************************/
/***********************************************************************************/
void InitEphemeris (LALStatus * status,   	/**< pointer to LALStatus structure */
		    EphemerisData *edat,	/**< [out] the ephemeris-data */
		    const CHAR *ephemDir,	/**< directory containing ephems */
		    const CHAR *ephemYear	/**< which years do we need? */
		    )
{

#define FNAME_LENGTH 1024
  
  CHAR EphemEarth[FNAME_LENGTH];	/* filename of earth-ephemeris data */
  CHAR EphemSun[FNAME_LENGTH];	/* filename of sun-ephemeris data */

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT ( edat, status, SIDEBANDUTILSC_ENULL,SIDEBANDUTILSC_MSGENULL );
  ASSERT ( ephemYear, status,SIDEBANDUTILSC_ENULL, SIDEBANDUTILSC_MSGENULL );

  if ( ephemDir )
    {
      snprintf(EphemEarth, FNAME_LENGTH, "%s/earth%s.dat", ephemDir, ephemYear);
      snprintf(EphemSun, FNAME_LENGTH, "%s/sun%s.dat", ephemDir, ephemYear);
    }
  else
    {
      snprintf(EphemEarth, FNAME_LENGTH, "earth%s.dat", ephemYear);
      snprintf(EphemSun, FNAME_LENGTH, "sun%s.dat",  ephemYear);
    }
  EphemEarth[FNAME_LENGTH-1]=0;
  EphemSun[FNAME_LENGTH-1]=0;
  
  /* NOTE: the 'ephiles' are ONLY ever used in LALInitBarycenter, which is
   * why we can use local variables (EphemEarth, EphemSun) to initialize them.
   */
  edat->ephiles.earthEphemeris = EphemEarth;
  edat->ephiles.sunEphemeris = EphemSun;

  TRY ( LALInitBarycenter(status->statusPtr, edat), status);

  DETATCHSTATUSPTR ( status );
  RETURN ( status );

} /* InitEphemeris() */


/***********************************************************************************/
/***********************************************************************************/
/* Read in the prior parameter ranges and proposal sizes  */
/***********************************************************************************/
/***********************************************************************************/
void ReadSideBandPriors(LALStatus *status,
			CHAR *rangefile,
			SideBandMCMCRanges *ranges,
			SideBandMCMCJumpProbs *jumpsizes
                        )
{
  
  FILE *fprange = NULL;
  CHAR ymin[32],ymin2[32];
  CHAR ymax[32],ymax2[32];
  CHAR yjump1[32],yjump2[32],yjump3[32];
  CHAR yprob1[32],yprob2[32],yprob3[32];
  CHAR y[32];
  CHAR line[352];

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT (ranges,status,SIDEBANDUTILSC_ENULL,SIDEBANDUTILSC_MSGENULL );
  ASSERT (jumpsizes,status,SIDEBANDUTILSC_ENULL,SIDEBANDUTILSC_MSGENULL );
  
  /* open and read parameter range file */
  if ((fprange = fopen(rangefile,"r"))==NULL) {
    XLALPrintError("\nError opening file '%s' for reading..\n\n",rangefile);
  }

  while (fgets(line,sizeof(line),fprange)!=NULL) {
    
    /* read in line from file */
    sscanf(line,"%s %s %s %s %s %s %s %s %s %s %s",y,ymin,ymin2,ymax,ymax2,yjump1,yprob1,yjump2,yprob2,yjump3,yprob3);
    
    /* deal with frequency */
    if (!strcmp(y,"f0")) {
      ranges->f0min = (REAL8)atof(ymin);
      ranges->f0max = (REAL8)atof(ymax);
      jumpsizes->f0.jump[0] = (REAL8)atof(yjump1);
      jumpsizes->f0.prob[0] = (REAL8)atof(yprob1);
      jumpsizes->f0.jump[1] = (REAL8)atof(yjump2);
      jumpsizes->f0.prob[1] = (REAL8)atof(yprob2);
      jumpsizes->f0.jump[2] = (REAL8)atof(yjump3);
      jumpsizes->f0.prob[2] = (REAL8)atof(yprob3);
    }
    /* deal with orbital period */
    if (!strcmp(y,"P")) {
      ranges->periodmin = (REAL8)atof(ymin);
      ranges->periodmax = (REAL8)atof(ymax);
      jumpsizes->period.jump[0] = (REAL8)atof(yjump1);
      jumpsizes->period.prob[0] = (REAL8)atof(yprob1);
      jumpsizes->period.jump[1] = (REAL8)atof(yjump2);
      jumpsizes->period.prob[1] = (REAL8)atof(yprob2);
      jumpsizes->period.jump[2] = (REAL8)atof(yjump3);
      jumpsizes->period.prob[2] = (REAL8)atof(yprob3);
    }
    /* deal with semimajor axis */
    if (!strcmp(y,"a")) {
      ranges->amin = (REAL8)atof(ymin);
      ranges->amax = (REAL8)atof(ymax);
      jumpsizes->a.jump[0] = (REAL8)atof(yjump1);
      jumpsizes->a.prob[0] = (REAL8)atof(yprob1);
      jumpsizes->a.jump[1] = (REAL8)atof(yjump2);
      jumpsizes->a.prob[1] = (REAL8)atof(yprob2);
      jumpsizes->a.jump[2] = (REAL8)atof(yjump3);
      jumpsizes->a.prob[2] = (REAL8)atof(yprob3);
    }
    /* deal with tp */
    else if (!strcmp(y,"tp")) {
      ranges->tpmin.gpsSeconds = (INT4)atoi(ymin);
      ranges->tpmin.gpsNanoSeconds = (INT4)atoi(ymin2);
      ranges->tpmax.gpsSeconds = (INT4)atoi(ymax);
      ranges->tpmax.gpsNanoSeconds = (INT4)atoi(ymax2);
      jumpsizes->x.jump[0] = (REAL8)atof(yjump1);
      jumpsizes->x.prob[0] = (REAL8)atof(yprob1);
      jumpsizes->x.jump[1] = (REAL8)atof(yjump2);
      jumpsizes->x.prob[1] = (REAL8)atof(yprob2);
      jumpsizes->x.jump[2] = (REAL8)atof(yjump3);
      jumpsizes->x.prob[2] = (REAL8)atof(yprob3);
    }
    /* deal with argp */
    else if (!strcmp(y,"argp")) {
      ranges->argpmin = (REAL8)atof(ymin);
      ranges->argpmax = (REAL8)atof(ymax);
      jumpsizes->y.jump[0] = (REAL8)atof(yjump1);
      jumpsizes->y.prob[0] = (REAL8)atof(yprob1);
      jumpsizes->y.jump[1] = (REAL8)atof(yjump2);
      jumpsizes->y.prob[1] = (REAL8)atof(yprob2);
      jumpsizes->y.jump[2] = (REAL8)atof(yjump3);
      jumpsizes->y.prob[2] = (REAL8)atof(yprob3);
    }
    /* deal with eccentricity */
    else if (!strcmp(y,"e")) {
      ranges->emin = (REAL8)atof(ymin);
      ranges->emax = (REAL8)atof(ymax);
      jumpsizes->e.jump[0] = (REAL8)atof(yjump1);
      jumpsizes->e.prob[0] = (REAL8)atof(yprob1);
      jumpsizes->e.jump[1] = (REAL8)atof(yjump2);
      jumpsizes->e.prob[1] = (REAL8)atof(yprob2);
      jumpsizes->e.jump[2] = (REAL8)atof(yjump3);
      jumpsizes->e.prob[2] = (REAL8)atof(yprob3);
    }
    /* deal with h0 */
    else if (!strcmp(y,"h0")) {
      ranges->h0min = (REAL8)atof(ymin);
      ranges->h0max = (REAL8)atof(ymax);
      jumpsizes->h0.jump[0] = (REAL8)atof(yjump1);
      jumpsizes->h0.prob[0] = (REAL8)atof(yprob1);
      jumpsizes->h0.jump[1] = (REAL8)atof(yjump2);
      jumpsizes->h0.prob[1] = (REAL8)atof(yprob2);
      jumpsizes->h0.jump[2] = (REAL8)atof(yjump3);
      jumpsizes->h0.prob[2] = (REAL8)atof(yprob3);
    }
    /* deal with cosi */
    else if (!strcmp(y,"cosi")) {
      ranges->cosimin = (REAL8)atof(ymin);
      ranges->cosimax = (REAL8)atof(ymax);
      jumpsizes->cosi.jump[0] = (REAL8)atof(yjump1);
      jumpsizes->cosi.prob[0] = (REAL8)atof(yprob1);
      jumpsizes->cosi.jump[1] = (REAL8)atof(yjump2);
      jumpsizes->cosi.prob[1] = (REAL8)atof(yprob2);
      jumpsizes->cosi.jump[2] = (REAL8)atof(yjump3);
      jumpsizes->cosi.prob[2] = (REAL8)atof(yprob3);
    }
    /* deal with psi */
    else if (!strcmp(y,"psi")) {
      ranges->psimin = (REAL8)atof(ymin);
      ranges->psimax = (REAL8)atof(ymax);
      jumpsizes->psi.jump[0] = (REAL8)atof(yjump1);
      jumpsizes->psi.prob[0] = (REAL8)atof(yprob1);
      jumpsizes->psi.jump[1] = (REAL8)atof(yjump2);
      jumpsizes->psi.prob[1] = (REAL8)atof(yprob2);
      jumpsizes->psi.jump[2] = (REAL8)atof(yjump3);
      jumpsizes->psi.prob[2] = (REAL8)atof(yprob3);
    } 
    /* deal with phi0 */
    else if (!strcmp(y,"phi0")) {
      ranges->phi0min = (REAL8)atof(ymin);
      ranges->phi0max = (REAL8)atof(ymax);
      jumpsizes->phi0.jump[0] = (REAL8)atof(yjump1);
      jumpsizes->phi0.prob[0] = (REAL8)atof(yprob1);
      jumpsizes->phi0.jump[1] = (REAL8)atof(yjump2);
      jumpsizes->phi0.prob[1] = (REAL8)atof(yprob2);
      jumpsizes->phi0.jump[2] = (REAL8)atof(yjump3);
      jumpsizes->phi0.prob[2] = (REAL8)atof(yprob3);
    }
    
  }
  fclose(fprange);
  
  DETATCHSTATUSPTR ( status );
  RETURN ( status );

}

/***********************************************************************************/
/***********************************************************************************/
/** select subset of frequencies from input frequencies based on parameter ranges  */
/***********************************************************************************/
/***********************************************************************************/
void SelectSideBandFrequencies (LALStatus * status,   
				SideBandDataSet **fulldata,
				SelectSideBandFrequencyParams *sfparams,
				SideBandDataSet **reddata)
     
{

  REAL8 dPmin,dPmax;
  REAL8 dSide = 1.0/LAL_DAYSID_SI;
  REAL8 minrange, maxrange;
  UINT4 nam = 1;
  UINT4 M;
  REAL8 x,y;
  INT4 i;
  INT4 j;
  UINT4 k = 0;
  INT4 mmin,mmax;
  REAL8 f;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  
  ASSERT ( fulldata, status, SIDEBANDUTILSC_ENULL,SIDEBANDUTILSC_MSGENULL );
  ASSERT ( reddata, status,SIDEBANDUTILSC_ENULL, SIDEBANDUTILSC_MSGENULL );
  
  mmin = sfparams->mmin;
  mmax = sfparams->mmax;
  dPmin = (1.0/sfparams->ranges.periodmax);
  dPmax = (1.0/sfparams->ranges.periodmin);
  M = mmax - mmin + 1;
  if (sfparams->am) nam = 5;

  /* allocate maximum required memory to reduced data */
  (*reddata)->freq = NULL;
  (*reddata)->fourier = NULL;
  (*reddata)->freq = XLALCreateREAL8Vector((*fulldata)->freq->length);
  (*reddata)->fourier = XLALCreateCOMPLEX16Vector((*fulldata)->freq->length);

  /* allocate memory to range vectors */
  sfparams->minf = NULL;
  sfparams->maxf = NULL;
  sfparams->minf = XLALCreateREAL8Vector(nam*M);
  sfparams->maxf = XLALCreateREAL8Vector(nam*M);

  /* generate list of frequency ranges in which the signal may fall */
  /* based on f0 error and Period error including AM sidebands*/
  for (i=mmin;i<0;i++) {
    for (j=0;j<(INT4)nam;j++) {
      sfparams->minf->data[5*(i-mmin)+j] = sfparams->ranges.f0min + i*dPmax + (j-2)*dSide - sfparams->df/2.0;
      sfparams->maxf->data[5*(i-mmin)+j] = sfparams->ranges.f0max + i*dPmin + (j-2)*dSide + sfparams->df/2.0; 
    }
  }
  for (j=0;j<(INT4)nam;j++) {
    sfparams->minf->data[5*(-mmin)+j] = sfparams->ranges.f0min + (j-2)*dSide - sfparams->df/2.0;
    sfparams->maxf->data[5*(-mmin)+j] = sfparams->ranges.f0max + (j-2)*dSide + sfparams->df/2.0;
  }
  for (i=1;i<=mmax;i++) {
    for (j=0;j<(INT4)nam;j++) {
      sfparams->minf->data[5*(-mmin+i)+j] = sfparams->ranges.f0min + i*dPmin + (j-2)*dSide - sfparams->df/2.0;
      sfparams->maxf->data[5*(-mmin+i)+j] = sfparams->ranges.f0max + i*dPmax + (j-2)*dSide + sfparams->df/2.0; 
    }
  }
  minrange = sfparams->minf->data[0];
  maxrange = sfparams->maxf->data[nam*M-1];

  for (i=0;i<(INT4)(nam*M-1);i++) printf("range = %6.12f -> %6.12f\n",sfparams->minf->data[i],sfparams->maxf->data[i]);

  /* loop over the input frequencies */
  for (j=0;j<(INT4)(*fulldata)->freq->length;j++) {

    f = (*fulldata)->freq->data[j];

    if ((f>minrange)&&(f<=maxrange)) { 
      
      x = creal((*fulldata)->fourier->data[j]);
      y = cimag((*fulldata)->fourier->data[j]);
      
      /* loop over the ranges */
      for (i=0;i<(INT4)(nam*M);i++) {

	/* if within a range then record it */
	if ((f>sfparams->minf->data[i])&&(f<=sfparams->maxf->data[i])) {
	  (*reddata)->freq->data[k] = (*fulldata)->freq->data[j];
	  (*reddata)->fourier->data[k].real_FIXME = x;
	  (*reddata)->fourier->data[k].imag_FIXME = y;

	  printf("found %6.12f between %6.12f -> %6.12f\n",(*reddata)->freq->data[k],sfparams->minf->data[i],sfparams->maxf->data[i]);

	  k++;
	}
	
      }
    
    }

  }

  /* resize the reduced vector */
  (*reddata)->freq = XLALResizeREAL8Vector((*reddata)->freq,k);
  (*reddata)->fourier = XLALResizeCOMPLEX16Vector((*reddata)->fourier,k);
  printf("Reduced data length = %d\n",k);

  DETATCHSTATUSPTR ( status );
  RETURN ( status );

} /* SelectFrequencies() */

/***********************************************************************************/
/***********************************************************************************/
/** Read input sky position demodulated data from file */
/***********************************************************************************/
/***********************************************************************************/
void ReadSideBandData (LALStatus * status,   
		       ReadSideBandDataParams *params,
		       SideBandDataSet **fulldata
		       )
{

  FILE *fp;
  BOOLEAN header = 1;                  /* flag for defining whether we are still reading fstat file header */
  CHAR line[1024];                      /* buffer for reading lines of fourier data into */
  INT4 headercount = 0;                /* counts how many lines in the header */
  INT4 datacount = 0;                  /* counts how many data points we have */
  fpos_t pos;                          /* records the position in the fstat file */
  CHAR x;                              /* hold first character of each Fstat file line to test if header data */
  REAL8 ftemp1, ftemp2;                /* temporary frequency variables */
  REAL8 dumf;                          /* REAL8 dummy variable */
  UINT4 i;
  UINT4 Ndata;
  REAL8 fre,fim;
  REAL8 norm;
  REAL8 sum = 0.0;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* open and read file */
  if ((fp = fopen(params->file,"r"))==NULL) {
    XLALPrintError("\nError opening file '%s' for reading..\n\n",params->file);
    /* return (SIDEBANDUTILSC_ESYS); */
  }

  /* open file to read header */
  while (header) {  
    
    fgetpos(fp,&pos);  /* record current position in file */
    if (fgets(line,1024,fp)==NULL) header = 0;  /* read current file line */
    printf("line is %s\n",line);
    x=line[0];
    /* if first character is a % then add to header count */
    printf("%s\n",&x);
    if (strncmp(&x,"%",1)==0) {
      headercount++;
      printf("added to headercount %d\n",headercount);
    }
    else {
      header = 0;
      printf("past header\n");
    }
  }
  
  /* return to end of header and read first 2 lines to assess frequency resolution */
  fsetpos(fp,&pos);
  int count;
  count = fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                 &ftemp1,&dumf,&dumf,&dumf,&dumf,&dumf,&dumf,&dumf,&dumf,&dumf,&dumf);
  if ( count != 11 ) {
    XLALPrintError ("\n fscanf() failed to read 11 items from 'fp'\n" );
    ABORT ( status, SIDEBANDUTILSC_EINPUT, SIDEBANDUTILSC_MSGEINPUT );
  }
  
  count = fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                 &ftemp2,&dumf,&dumf,&dumf,&dumf,&dumf,&dumf,&dumf,&dumf,&dumf,&dumf);
  if ( count != 11 ) {
    XLALPrintError ("\n fscanf() failed to read 11 items from 'fp'\n" );
    ABORT ( status, SIDEBANDUTILSC_EINPUT, SIDEBANDUTILSC_MSGEINPUT );
  }

  params->df = ftemp2 - ftemp1;
  
  printf("ftemp1 = %6.12f ftemp2 = %6.12f\n",ftemp1,ftemp2);
  printf("df = %6.20f\n",params->df);
    
  /* return to end of header and read data to determine frequency range */
  fsetpos(fp,&pos);
  while (fgets(line,512,fp)!=NULL) {
    x=line[0];
    if (strncmp(&x,"%",1)!=0) {
      sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	     &ftemp1,&dumf,&dumf,&dumf,&dumf,&dumf,&dumf,&dumf,&dumf,&dumf,&dumf);
      
      /* if the frequency is in the range then add to data count */
      if ((ftemp1>=params->minf)&&(ftemp1<=params->maxf)) datacount++;
    
    }
  }
  Ndata = datacount;
  printf("Ndata = %d\n",Ndata);
  
  /* allocate memory */
  (*fulldata)->fourier = NULL;
  (*fulldata)->freq = NULL;
  (*fulldata)->fourier = XLALCreateCOMPLEX16Vector(Ndata);
  (*fulldata)->freq = XLALCreateREAL8Vector(Ndata);

  /* using this normalisation we convert from the Power spectrum to one sided Power spectral density */ 
  norm = 2.0/sqrt(params->Tobs);

  /* return to start of data and read in full fa and fb data */
  fsetpos(fp,&pos);
  i = 0;
  while(fgets(line,512,fp)!=NULL) {
    x=line[0];
    if (strncmp(&x,"%",1)!=0) {
      sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	     &ftemp1,&dumf,&dumf,&dumf,&dumf,&dumf,&fre,&fim,&dumf,&dumf,
	     &dumf);
      if ((ftemp1>=params->minf)&&(ftemp1<=params->maxf)) {
	(*fulldata)->freq->data[i] = ftemp1;

	/* normalise data by the input sqrt noise spectral density */
	(*fulldata)->fourier->data[i].real_FIXME = fre*norm;
	(*fulldata)->fourier->data[i].imag_FIXME = fim*norm;
	/* printf("Reading data : %6.12e %6.12e\n",(*fulldata)->fourier->data[i].re,(*fulldata)->fourier->data[i].im); */
	sum += (fre*fre+fim*fim)*norm*norm;
	i++;
      }
    }
  }
  printf("av = %e\n",sqrt(sum/(REAL8)i));
  fclose(fp);

  DETATCHSTATUSPTR ( status );
  RETURN ( status );

} /* ReadFourierData() */

/***********************************************************************************/
/***********************************************************************************/
/* Estimate noise from input data (avoiding frequency bins containing the signal) */
/***********************************************************************************/
/***********************************************************************************/
void EstimateSideBandNoise (LALStatus * status,   
			    SideBandDataSet **fulldata,
			    EstimateSideBandNoiseParams *ENparams
			    )
{

  REAL8 sum = 0.0;
  REAL8 f;
  INT4 i,j;
  INT4 k = 0;
  BOOLEAN flag;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT ( *fulldata, status, SIDEBANDUTILSC_ENULL,SIDEBANDUTILSC_MSGENULL );

  /* loop over full data set frequencies */
  for (i=0;i<(INT4)(*fulldata)->freq->length;i++) {
    
    f = (*fulldata)->freq->data[i];
    flag = 0;

    /* if data is within desired band */
    if ((f>ENparams->minfreq)&&(f<=ENparams->maxfreq)) {
      
      /* loop over signal data ranges */
      j = 0;
      while ((j<(INT4)ENparams->minf->length)&&(flag==0)) {
	
	/* if data is within signal data range */
	if ((f>ENparams->minf->data[j]-ENparams->safedf)&&(f<=ENparams->maxf->data[j]+ENparams->safedf)) flag = 1;

	j++;

      }

      if (flag==0) {
	sum += creal((*fulldata)->fourier->data[i])*creal((*fulldata)->fourier->data[i]) + cimag((*fulldata)->fourier->data[i])*cimag((*fulldata)->fourier->data[i]);
	k++;
      }
	
    }

  }

  if (k>=ENparams->Nthresh) {
    printf("average Sh = %e\n",sum/(REAL8)k);
    printf("sum = %e k = %d\n",sum,k);
    ENparams->sqrtSh = sqrt(sum/(REAL8)k);
    printf("Computed sqrtSh = %6.12e from %d points\n",ENparams->sqrtSh,k);
  }
  else {
    printf("ERROR : Not enough data points to compute noise floor, exiting.\n");
    exit(0);
  }


  /* normalise data */
  /* sum = 0.0; */
  /* loop over full data set frequencies */
  /* for (i=0;i<(*fulldata)->freq->length;i++) {
    (*fulldata)->fourier->data[i].re *= sqrt(2.0)/ENparams->sqrtSh;
    (*fulldata)->fourier->data[i].im *= sqrt(2.0)/ENparams->sqrtSh;
    printf("Normalised data : %6.12f %6.12f\n",(*fulldata)->fourier->data[i].re,(*fulldata)->fourier->data[i].im); 
    sum += (*fulldata)->fourier->data[i].re*(*fulldata)->fourier->data[i].re;
  } */


  DETATCHSTATUSPTR ( status );
  RETURN ( status );

} /* EstimateSideBandNoise() */


/***********************************************************************************/
/***********************************************************************************/
/* Compute the signal waveform at requested frequencies */
/***********************************************************************************/
/***********************************************************************************/
void GenerateSideBandTemplate (LALStatus *status,   			/**< pointer to LALStatus structure */
			       BinarySourceParams *BSParams,	/**< Structure containing Binary source parameters */
			       SideBandTemplateParams *TParams,	/**< structure containing template parameters */
			       SideBandTemplate **Template	        /**< Output frequency domain template */
			       )
{
  
  REAL8 h0;
  REAL8 psi;
  REAL8 phi0;
  REAL8 cosi;
  REAL8 A1,A2,A3,A4;
  REAL8 b;
  REAL8 T;
  REAL8 tp;
  REAL8 f0;
  INT4 i;
  REAL8 freq;
  INT4 n;
  REAL8 ftemp;
  REAL8 y;
  REAL8 x;
  REAL8 P;
  REAL8 argp;
  INT4 nmax;
  INT4 nmin;
  REAL8 kappa;
  REAL8 gam;
  REAL8 R1,R2;
  REAL8 e;
  REAL8 a;
  INT4 p;
  REAL8 bn;
  REAL8 bp;
  INT4 q;
  INT4 winindex;
  REAL8 cosy,siny;
  REAL8 tdiff;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* check input results vectors */ 
  ASSERT((*Template)->fourier,status,SIDEBANDUTILSC_ENULL,SIDEBANDUTILSC_MSGENULL);
  ASSERT(TParams->wa,status,SIDEBANDUTILSC_ENULL,SIDEBANDUTILSC_MSGENULL);
  ASSERT(TParams->wb,status,SIDEBANDUTILSC_ENULL,SIDEBANDUTILSC_MSGENULL);

  /* check input frequency vector */
  ASSERT(TParams->freqsamples,status,SIDEBANDUTILSC_ENULL,SIDEBANDUTILSC_MSGENULL); 

  if (lalDebugLevel) printf ("\nChecked input structures in SideBandTemplate.\n");
  
  /* compute sideband spacing */
  /* currently unused: REAL8 dfp = 1.0/BSParams->OrbitalPeriod; */
     
  if (lalDebugLevel) printf ("\nInside fourier option in SideBandTemplate.\n");

  /* redefine nuisance parameters */
  h0 = BSParams->h0;
  psi = BSParams->psi;
  phi0 = BSParams->phi0;
  cosi = BSParams->cosi;
  /* printf("reftime is %d %d\n",TParams->reftime.gpsSeconds,TParams->reftime.gpsNanoSeconds); 
     printf("starttime is %d %d\n",TParams->tstart.gpsSeconds,TParams->tstart.gpsNanoSeconds); 
     printf("phi0 = %f\n",phi0); 
     printf("temp = %f\n",temp); */
  
  /* redefine orbital parameters */
  e = BSParams->OrbitalEccentricity;
  P = BSParams->OrbitalPeriod;
  argp = BSParams->ArgumentofPeriapse;
  a = BSParams->OrbitalSemiMajorAxis;
  f0 =  BSParams->f0;

  /* modify initial phase to account for phase reference time */
  tdiff = TParams->tstart.gpsSeconds+1e-9*TParams->tstart.gpsNanoSeconds-TParams->reftime.gpsSeconds-1e-9*TParams->reftime.gpsNanoSeconds;
  phi0 = fmod((phi0 + LAL_TWOPI*tdiff*f0),LAL_TWOPI);

  /* redefine some detector and dataset parameters */
  T = TParams->Tobs;
  /* currently unused: Sh = TParams->sqrtSh*TParams->sqrtSh; */
  /* currently unused: REAL8 sPe = 86400.0*(1.0/((1.0/365.2524)+1.0)); // define the sidereal day */
  /* currently unused: REAL8 w0 = TParams->ABCco->omega0; */
  
  /* define the eccentricity expansion coefficients */
  R1 = LAL_TWOPI*BSParams->f0*a*sqrt(1.0-e*e*cos(argp)*cos(argp));
  R2 = R1*e/2.0;
  
  /* define some intermediate variables */
  kappa = LAL_TWOPI*BSParams->f0*a*sin(argp);
  gam = atan2(tan(argp),sqrt(1-e*e));
  
  /* make sure gam is correctly defined */
  if ((argp>=LAL_PI/2.0)&&(argp<3.0*LAL_PI/2.0)) gam += LAL_PI;
  else if (argp>=3.0*LAL_PI/2.0) gam += LAL_TWOPI;
  /* printf("gam = %f\n",gam);*/ 
  
  /* printf("BinaryFDTemplate : A = %f B = %f C = %f D = %f\n",A,B,C,D);
     printf("P = %6.12f argp = %f sPe = %6.12f w0 = %f\n",P,argp,sPe,w0);
     printf("h0 = %f psi = %f phi0 = %f cosi = %f\n",h0,psi,phi0,cosi); */
  
  /* test phi0 shift possibility */
  /* phi0 = fmod((phi0 - (TParams->tstart.gpsSeconds+1e-9*TParams->tstart.gpsNanoSeconds - TParams->tstartSSB.gpsSeconds+1e-9*TParams->tstartSSB.gpsNanoSeconds)*f0),LAL_TWOPI);*/
  
  /* compute JKS signal amplitudes */
  A1 = h0*(0.5*(1.0+cosi*cosi)*cos(2.0*psi)*cos(phi0)-cosi*sin(2.0*psi)*sin(phi0));
  A2 = h0*(0.5*(1.0+cosi*cosi)*sin(2.0*psi)*cos(phi0)+cosi*cos(2.0*psi)*sin(phi0));
  A3 = h0*(-0.5*(1+cosi*cosi)*cos(2.0*psi)*sin(phi0)-cosi*sin(2.0*psi)*cos(phi0));
  A4 = h0*(-0.5*(1+cosi*cosi)*sin(2.0*psi)*sin(phi0)+cosi*cos(2.0*psi)*cos(phi0));
    
  /* printf("A1 = %6.12f A2 = %6.12f A3 = %6.12f A4 = %6.12f\n",A1,A2,A3,A4); */
  
  if (lalDebugLevel) printf ("\nComputed signal amplitudes in BinaryFDTemplate.\n");
  
  /* compute orbital phase (for circular orbit only) */
  /* tp = XLALGPSDiff(&(TParams->tstart),&(BSParams->TimeofSSBPeriapsePassage)); */            /* compute the time since periapse passage */
  /* tp = XLALGPSDiff(&(BSParams->TimeofSSBPeriapsePassage),&(TParams->tstart));   */          /* compute the time since periapse passage */
  /* printf("tpSSB = %d %d\n",BSParams->TimeofSSBPeriapsePassage.gpsSeconds,BSParams->TimeofSSBPeriapsePassage.gpsNanoSeconds);
     printf("tstartSSB = %d %d\n",TParams->tstartSSB.gpsSeconds,TParams->tstartSSB.gpsNanoSeconds); */
  
  /* ------------------------------------------------------------------------------------------------------------------------------------------- */
  /* changed the following line to include light time travel from SSB to detector because when we compute Fstat we do not barycenter the binary */
  /* phase component, only the intrinsic phase component */
  tp = BSParams->TimeofSSBPeriapsePassage.gpsSeconds + 1e-9*BSParams->TimeofSSBPeriapsePassage.gpsNanoSeconds 
    - TParams->tstart.gpsSeconds - 1e-9*TParams->tstart.gpsNanoSeconds; 
  /* ------------------------------------------------------------------------------------------------------------------------------------------- */
  /* omega0 = LAL_TWOPI*fmod(tdiff,BSParams->OrbitalPeriod)/BSParams->OrbitalPeriod; */               /* compute the corresponding initial orbital phase */ 
  
  /* printf("new tp = %6.12f\n",tp); */
  
  /* initialise the results vectors */
  for (i=0;i<(INT4)(*Template)->fourier->length;i++) {
    (*Template)->fourier->data[i].real_FIXME = 0.0;
    (*Template)->fourier->data[i].imag_FIXME = 0.0;
  }
  
  if (lalDebugLevel) printf ("\nInitialised results vectors in BinaryFDTemplate.\n");
  
  /* printf("in BinaryFDTemplate e = %6.12f\n",e); */
  
  /* printf("freqsamples->length = %d\n",TParams->freqsamples->length); */
  
  /* loop over each input frequency sample */
  for (i=0;i<(INT4)TParams->freqsamples->length;i++) {
    
    /* define closest index to desired frequency and record new frequency */
    freq = TParams->freqsamples->data[i];
    /* printf("freq = %f\n",freq); */
    
    /* determine whether we are doing this locally or not */
    if (TParams->local) {
      nmin = floor((freq - f0)*P + 0.5);
      nmax = nmin;
    }
    else {
      /* the extent to which the window has been calculated out to */
      nmax = floor((freq - f0)*P + 0.5) + floor(TParams->windowrange*TParams->dfwindow*P) - 1;
      nmin = floor((freq - f0)*P + 0.5) - floor(TParams->windowrange*TParams->dfwindow*P) + 1;
    }
    /* printf("freq = %9.12f f0 = %9.12f P = %f windowrange = %d dfhighres = %e\n",freq,f0,P,TParams->windowrange,TParams->dfhighres);
       printf("loop %d of %d nmin = %d nmax = %d\n",i,TParams->freqsamples->length,nmin,nmax); */
    
    /* now loop over the sidebands we wish to take into account for this particular frequency */
    /* setting local means that we only use the nearest sideband (zeroth order e) */
    for (n=nmin;n<=nmax;n++) {
      
      /* define central frequency of this binary sideband */
      ftemp = f0 + ((REAL8)n/P);
      /* printf("ftemp = %9.12f freq = %9.12f\n",ftemp,freq); */
      
      /* now loop over the extent that the eccentricity spreads power amongst nearby sidebands (first order e)*/
      for (p=-TParams->pmax;p<=TParams->pmax;p++) {
	
	/* compute bessel functions */
	bn = bessj(n-2*p,R1);
	bp = bessj(p,R2);
	b = 2.0*bn*bp/sqrt(T);
	
	/* define phase summation index */
	q = n - p;
	
	/* compute phase of this sideband */
	y = n*(-LAL_TWOPI*tp/P) + q*(gam + LAL_PI) + 3.0*kappa*e/2.0;
	
	/* compute y phase factor */
	cosy = cos(y);
	siny = sin(y);
	
	/* compute signal profile at this frequency */
	x = freq - ftemp;
	/* printf("x = %9.12f\n",x);
	   printf("highresdf = %9.12f\n",TParams->dfhighres); */
	
	/* compute window index */
	winindex = floor(0.5 + x/TParams->dfwindow) + TParams->windowrange;
	/* printf("winindex = %d TParams->windowrange = %d\n",winindex,TParams->windowrange); */ 
	/* compute Rn and Sn */  
	if ((winindex<2*TParams->windowrange)&&(winindex>=0)) {
	  (*Template)->fourier->data[i].real_FIXME = creal((*Template)->fourier->data[i]) 
	    + 0.5*b*(creal(TParams->wa->data[winindex])*(A1*cosy + A3*siny) 
		     + cimag(TParams->wa->data[winindex])*(A3*cosy - A1*siny)
		     + creal(TParams->wb->data[winindex])*(A2*cosy + A4*siny) 
		     + cimag(TParams->wb->data[winindex])*(A4*cosy - A2*siny));
	  (*Template)->fourier->data[i].imag_FIXME = cimag((*Template)->fourier->data[i]) 
	    + 0.5*b*(creal(TParams->wa->data[winindex])*(A1*siny - A3*cosy) 
		     + cimag(TParams->wa->data[winindex])*(A1*cosy + A3*siny)
		     + creal(TParams->wb->data[winindex])*(A2*siny - A4*cosy) 
		     + cimag(TParams->wb->data[winindex])*(A2*cosy + A4*siny));
	}	
	
      } /* end loop over first order eccentricity terms */
      
    } /* end loop over binary sidebands */
    
  } /* end loop over frequency samples */
   
  DETATCHSTATUSPTR (status);
  RETURN (status);
  
}

/***********************************************************************************/
/***********************************************************************************/
/* The following three functions are borrowed from numerical recipes in */
/* order to compute the bessel function of the first kind */
/***********************************************************************************/
/***********************************************************************************/
REAL4 bessj(INT4 n, REAL4 x)
{

	/* void nrerror(char error_text[]); */
	INT4 j,jsum,m,sign;
	REAL4 ax,bj,bjm,bjp,sum,tox,ans;

	sign = 1;
	if (n<0) {
	  sign = pow((-1.0),n);
	  n = fabs(n);
	}
	
	if (n==0) return sign*bessj0(x);
	if (n==1) return sign*bessj1(x);

	ax=fabs(x);
	if (ax == 0.0)
		return 0.0;
	else if (ax > (REAL4) n) {
		tox=2.0/ax;
		bjm=bessj0(ax);
		bj=bessj1(ax);
		for (j=1;j<n;j++) {
			bjp=j*tox*bj-bjm;
			bjm=bj;
			bj=bjp;
		}
		ans=bj;
	} else {
		tox=2.0/ax;
		m=2*((n+(INT4) sqrt(ACC*n))/2);
		jsum=0;
		bjp=ans=sum=0.0;
		bj=1.0;
		for (j=m;j>0;j--) {
			bjm=j*tox*bj-bjp;
			bjp=bj;
			bj=bjm;
			if (fabs(bj) > BIGNO) {
				bj *= BIGNI;
				bjp *= BIGNI;
				ans *= BIGNI;
				sum *= BIGNI;
			}
			if (jsum) sum += bj;
			jsum=!jsum;
			if (j == n) ans=bjp;
		}
		sum=2.0*sum-bj;
		ans /= sum;
	}
	return x < 0.0 && (n & 1) ? -ans*sign : ans*sign;
}

REAL4 bessj0(REAL4 x)
{
	REAL4 ax,z;
	REAL8 xx,y,ans,ans1,ans2;

	if ((ax=fabs(x)) < 8.0) {
		y=x*x;
		ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
			+y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
		ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
			+y*(59272.64853+y*(267.8532712+y*1.0))));
		ans=ans1/ans2;
	} else {
		z=8.0/ax;
		y=z*z;
		xx=ax-0.785398164;
		ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
			+y*(-0.2073370639e-5+y*0.2093887211e-6)));
		ans2 = -0.1562499995e-1+y*(0.1430488765e-3
			+y*(-0.6911147651e-5+y*(0.7621095161e-6
			-y*0.934935152e-7)));
		ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
	}
	return ans;
}

REAL4 bessj1(REAL4 x)
{
	REAL4 ax,z;
	REAL8 xx,y,ans,ans1,ans2;

	if ((ax=fabs(x)) < 8.0) {
		y=x*x;
		ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
			+y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
		ans2=144725228442.0+y*(2300535178.0+y*(18583304.74
			+y*(99447.43394+y*(376.9991397+y*1.0))));
		ans=ans1/ans2;
	} else {
		z=8.0/ax;
		y=z*z;
		xx=ax-2.356194491;
		ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
			+y*(0.2457520174e-5+y*(-0.240337019e-6))));
		ans2=0.04687499995+y*(-0.2002690873e-3
			+y*(0.8449199096e-5+y*(-0.88228987e-6
			+y*0.105787412e-6)));
		ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
		if (x < 0.0) ans = -ans;
	}
	return ans;
}

/***********************************************************************************/
/***********************************************************************************/
/* ComputeABCcoefficients computes the coefficents of the expansion of a(t), */
/* b(t), a(t)*a(t), b(t)*b(t) and a(t)*b(t) in terms of harmonics of the */
/* earth spin frequency.*/
/***********************************************************************************/
/***********************************************************************************/
void
ComputeABCcoefficients (LALStatus *status,   			/**< pointer to LALStatus structure */
			ABCcoParams *abcparams,     	/**< structure containing the GPS observation start time and the sky position */
			LALDetector *detector,	        /**< structure containing the detector location and orientation parameters */
			ABCcoefficients **abcco	        /**< ABC expansion coefficients */
			)
{
  
  REAL8 alpha;
  REAL8 delta;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* check input results vectors */
  ASSERT(detector,status,SIDEBANDUTILSC_ENULL,SIDEBANDUTILSC_MSGENULL);
  ASSERT((*abcco),status,SIDEBANDUTILSC_ENULL,SIDEBANDUTILSC_MSGENULL);
  ASSERT(abcparams,status,SIDEBANDUTILSC_ENULL,SIDEBANDUTILSC_MSGENULL);

  printf("%6.12f %6.12f %6.12f\n",detector->response[0][0],detector->response[0][1],detector->response[0][2]);
  printf("%6.12f %6.12f %6.12f\n",detector->response[1][0],detector->response[1][1],detector->response[1][2]);
  printf("%6.12f %6.12f %6.12f\n",detector->response[2][0],detector->response[2][1],detector->response[2][2]);

  /* redefine some parameters for conciseness */
  alpha = abcparams->alpha;
  delta = abcparams->delta;
  printf("alpha = %f delta = %f\n",alpha,delta);

  /* compute sidereal time at detector at observation start (equal to alpha - GMST(t0)) */
  (*abcco)->omega0 = alpha - XLALGreenwichMeanSiderealTime(&(abcparams->tstart));

  if (lalDebugLevel) printf ("\nFinished computing omega0 = %f\n",(*abcco)->omega0);

  /* here we will compute the coefficients following Anderson et al 2001 */
  (*abcco)->aco[0] = 0.5*cos(delta)*cos(delta)*(detector->response[0][0] + detector->response[1][1] - 2.0*detector->response[2][2]);
  (*abcco)->aco[1] = 2.0*detector->response[0][2]*sin(delta)*cos(delta);
  (*abcco)->aco[2] = 0.5*(detector->response[1][1]-detector->response[0][0])*(1.0 + sin(delta)*sin(delta));
  
  (*abcco)->apco[0] = 0.0;
  (*abcco)->apco[1] = 2.0*detector->response[1][2]*sin(delta)*cos(delta);
  (*abcco)->apco[2] = (-1.0)*detector->response[0][1]*(1.0 + sin(delta)*sin(delta));
  
  (*abcco)->bco[0] = 0.0;
  (*abcco)->bco[1] = -2.0*detector->response[1][2]*cos(delta);
  (*abcco)->bco[2] = 2.0*detector->response[1][0]*sin(delta);
  (*abcco)->bpco[0] = 0.0;
  (*abcco)->bpco[1] = 2.0*detector->response[2][0]*cos(delta);
  (*abcco)->bpco[2] = -1.0*(detector->response[0][0]-detector->response[1][1])*sin(delta); 

  if (lalDebugLevel) printf ("\nFinished computing a(t) and b(t) coefficients.\n");
 
  /* compute expansion coefficients for a(t)*a(t) */
  (*abcco)->uco[0] = (*abcco)->aco[0]*(*abcco)->aco[0] 
    + 0.5*((*abcco)->aco[1]*(*abcco)->aco[1] + (*abcco)->aco[2]*(*abcco)->aco[2] 
	   + (*abcco)->apco[1]*(*abcco)->apco[1] + (*abcco)->apco[2]*(*abcco)->apco[2]);
  (*abcco)->uco[1] = (*abcco)->aco[1]*(*abcco)->aco[2] 
    + (*abcco)->apco[1]*(*abcco)->apco[2] + 2.0*(*abcco)->aco[0]*(*abcco)->aco[1];
  (*abcco)->uco[2]  = 0.5*((*abcco)->aco[1]*(*abcco)->aco[1] 
				 - (*abcco)->apco[1]*(*abcco)->apco[1]) + 2.0*(*abcco)->aco[0]*(*abcco)->aco[2];
  (*abcco)->uco[3] = (*abcco)->aco[1]*(*abcco)->aco[2] - (*abcco)->apco[1]*(*abcco)->apco[2];
  (*abcco)->uco[4] = 0.5*((*abcco)->aco[2]*(*abcco)->aco[2] - (*abcco)->apco[2]*(*abcco)->apco[2]);
  (*abcco)->upco[0] = 0.0;
  (*abcco)->upco[1] = ((*abcco)->aco[1]*(*abcco)->apco[2] 
			     - (*abcco)->aco[2]*(*abcco)->apco[1]) + 2.0*(*abcco)->aco[0]*(*abcco)->apco[1];
  (*abcco)->upco[2] = (*abcco)->aco[1]*(*abcco)->apco[1] + 2.0*(*abcco)->aco[0]*(*abcco)->apco[2];
  (*abcco)->upco[3] = ((*abcco)->aco[1]*(*abcco)->apco[2] + (*abcco)->aco[2]*(*abcco)->apco[1]);
  (*abcco)->upco[4] = (*abcco)->aco[2]*(*abcco)->apco[2];

  if (lalDebugLevel) printf ("\nFinished computing a(t)*a(t) coefficients.\n");
 
  /* compute expansion coefficients for b(t)*b(t) */
  (*abcco)->wco[0] = (*abcco)->bco[0]*(*abcco)->bco[0] 
    + 0.5*((*abcco)->bco[1]*(*abcco)->bco[1] 
	   + (*abcco)->bco[2]*(*abcco)->bco[2] + (*abcco)->bpco[1]*(*abcco)->bpco[1] 
	   + (*abcco)->bpco[2]*(*abcco)->bpco[2]);
  printf("%f\n",(*abcco)->wco[0]);
  (*abcco)->wco[1] = (*abcco)->bco[1]*(*abcco)->bco[2] 
    + (*abcco)->bpco[1]*(*abcco)->apco[2] + 2.0*(*abcco)->bco[0]*(*abcco)->bco[1];
  printf("%f\n",(*abcco)->wco[1]);
  (*abcco)->wco[2]  = 0.5*((*abcco)->bco[1]*(*abcco)->bco[1] 
				 - (*abcco)->bpco[1]*(*abcco)->bpco[1]) + 2.0*(*abcco)->bco[0]*(*abcco)->bco[2];
  printf("%f\n",(*abcco)->wco[2]);
  (*abcco)->wco[3] = (*abcco)->bco[1]*(*abcco)->bco[2] - (*abcco)->bpco[1]*(*abcco)->bpco[2];
  printf("%f\n",(*abcco)->wco[3]);
  (*abcco)->wco[4] = 0.5*((*abcco)->bco[2]*(*abcco)->bco[2] - (*abcco)->bpco[2]*(*abcco)->bpco[2]);
  printf("%f\n",(*abcco)->wco[4]);

  (*abcco)->wpco[0] = 0.0;
  printf("%f\n",(*abcco)->wpco[0]);
  (*abcco)->wpco[1] = ((*abcco)->bco[1]*(*abcco)->bpco[2] 
			     - (*abcco)->bco[2]*(*abcco)->bpco[1]) + 2.0*(*abcco)->bco[0]*(*abcco)->bpco[1];
  printf("%f\n",(*abcco)->wpco[1]);
  (*abcco)->wpco[2] = (*abcco)->bco[1]*(*abcco)->bpco[1] + 2.0*(*abcco)->bco[0]*(*abcco)->bpco[2];
  printf("%f\n",(*abcco)->wpco[2]);
  (*abcco)->wpco[3] = ((*abcco)->bco[1]*(*abcco)->bpco[2] + (*abcco)->bco[2]*(*abcco)->bpco[1]);
  printf("%f\n",(*abcco)->wpco[3]);
  (*abcco)->wpco[4] = (*abcco)->bco[2]*(*abcco)->bpco[2];
  printf("%f\n",(*abcco)->wpco[4]);
  
  if (lalDebugLevel) printf ("\nFinished computing b(t)*b(t) coefficients.\n");

  /* compute expansion coefficients for a(t)*b(t) */
  (*abcco)->vco[0] = 0.5*((*abcco)->aco[1]*(*abcco)->bco[1] 
				+ (*abcco)->aco[2]*(*abcco)->bco[2] + (*abcco)->apco[1]*(*abcco)->bpco[1] 
				+ (*abcco)->apco[2]*(*abcco)->bpco[2]);
  (*abcco)->vco[1] = (*abcco)->aco[0]*(*abcco)->bco[1] 
    + 0.5*((*abcco)->aco[1]*(*abcco)->bco[2] + (*abcco)->aco[2]*(*abcco)->bco[1] 
	   + (*abcco)->apco[1]*(*abcco)->bpco[2] + (*abcco)->apco[2]*(*abcco)->bpco[1]);
  (*abcco)->vco[2]  = (*abcco)->aco[0]*(*abcco)->bco[2] 
    + 0.5*((*abcco)->aco[1]*(*abcco)->bco[1] - (*abcco)->apco[1]*(*abcco)->bpco[1]);
  (*abcco)->vco[3] = 0.5*((*abcco)->aco[1]*(*abcco)->bco[2] 
				+ (*abcco)->aco[2]*(*abcco)->bco[1] - (*abcco)->apco[1]*(*abcco)->bpco[2] 
				- (*abcco)->apco[2]*(*abcco)->bpco[1]);
  (*abcco)->vco[4] = 0.5*((*abcco)->aco[2]*(*abcco)->bco[2] - (*abcco)->apco[2]*(*abcco)->bpco[2]);
  (*abcco)->vpco[0] = 0.0;
  (*abcco)->vpco[1] = (*abcco)->aco[0]*(*abcco)->bpco[1] 
    + 0.5*((*abcco)->aco[1]*(*abcco)->bpco[2] - (*abcco)->aco[2]*(*abcco)->bpco[1] 
	   - (*abcco)->apco[1]*(*abcco)->bco[2] + (*abcco)->apco[2]*(*abcco)->bco[1]);
  (*abcco)->vpco[2] = (*abcco)->aco[0]*(*abcco)->bpco[2] 
    + 0.5*((*abcco)->aco[1]*(*abcco)->bpco[1] + (*abcco)->apco[1]*(*abcco)->bco[1]);
  (*abcco)->vpco[3] = 0.5*((*abcco)->aco[1]*(*abcco)->bpco[2] 
				 + (*abcco)->aco[2]*(*abcco)->bpco[1] + (*abcco)->apco[1]*(*abcco)->bco[1] 
				 + (*abcco)->apco[2]*(*abcco)->bco[1]);
  (*abcco)->vpco[4] = 0.5*((*abcco)->aco[2]*(*abcco)->bpco[2] + (*abcco)->apco[2]*(*abcco)->bco[2]);

  if (lalDebugLevel) printf ("\nFinished computing a(t)*b(t) coefficients.\n");


  DETATCHSTATUSPTR (status);
  RETURN (status);

}

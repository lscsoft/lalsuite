/*
*  Copyright (C) 2007 Alexander Dietz, Duncan Brown, Thomas Cokelaer
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

/*-----------------------------------------------------------------------
 *
 * File Name: FindChirpBCVCFilter.c
 *
 * Author: Brown, D. A.,  Messaritaki E. and Cokelaer, T.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpBCVCFilterCV">
Author: Thomas Cokelaer
$Id$
</lalVerbatim>

<lalLaTeX>
\input{FindChirpBCVCFilterCDoc}

\vfill{\footnotesize\input{FindChirpBCVCFilterCV}}
</lalLaTeX>
#endif

#include <math.h>
#include <lal/LALErrno.h>
#include <lal/XLALError.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpBCV.h>

double rint(double x);

NRCSID (FINDCHIRPBCVCFILTERC, "$Id$");


/* <lalVerbatim file="FindChirpBCVCFilterCP"> */
void
LALFindChirpBCVCFilterSegment (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params
    )
/* </lalVerbatim> */
{
  UINT4                 j, k, kFinal/*, kOpt*/;
  UINT4                 numPoints;
  UINT4                 deltaEventIndex = 0;
  UINT4                 ignoreIndex;
  REAL4                 myfmin;
  REAL4                 deltaT, deltaF;
  REAL4                 norm;
  REAL4                 modqsqThresh;
  REAL4                 rhosqThresh;
  REAL4                 mismatch;
  REAL4                 chisqThreshFac = 0;
  /* REAL4                 modChisqThresh; */
  UINT4                 numChisqBins = 0;
  UINT4                 eventStartIdx = 0;
  UINT4                *chisqBin      = NULL;
  UINT4                *chisqBinBCV   = NULL;
  REAL4                 chirpTime     = 0;
  REAL4                *tmpltPower    = NULL;
  REAL4                *tmpltPowerBCV = NULL;
  COMPLEX8             *qtilde        = NULL;
  COMPLEX8             *qtildeBCV     = NULL;
  COMPLEX8             *q             = NULL;
  COMPLEX8             *qBCV          = NULL;
  COMPLEX8             *inputData     = NULL;
  COMPLEX8             *inputDataBCV = NULL;
  COMPLEX8             *tmpltSignal   = NULL;
  SnglInspiralTable    *thisEvent     = NULL;
  REAL4                 a1 = 0.0;
  REAL4                 b1 = 0.0;
  REAL4                 b2 = 0.0;
  REAL4                 templateNorm;
  REAL4                 m, psi0, psi3, fFinal;



  /* some declarations for the constrained  BCV codet (Thomas)*/
  REAL4         alphaUnity;     /* intermediate variable to get thetab  */
  REAL4         thetab;         /* critical angle for constraint        */
  REAL4         thetav;         /* angle between V1 and V2              */
  REAL4         V0;             /* intermediate quantity to compute rho */
  REAL4         V1;             /* intermediate quantity to compute rho */
  REAL4         V2;             /* intermediate quantity to compute rho */
  REAL4         alphaC;         /* constraint alpha value               */
  REAL4         alphaU;         /* unconstraint alpha value for fun     */
  REAL4         rhosqConstraint;
  REAL4         rhosqUnconstraint;      /* before constraint            */
  REAL4         deltaTPower2By3;        /* alias variable               */
  REAL4         deltaTPower1By6;        /* alias variable               */
  REAL4         ThreeByFive = 3.0/5.0;
  REAL4         FiveByThree = 5.0/3.0;
  REAL4         mchirp;
  REAL4         eta;
  REAL4         SQRTNormA1;             /* an alias                     */

  /* Most of the code below is identical and should be identical to
   * FindChirpBCVFilter.c execpt when entering the constraint code.
   * */


  INITSTATUS( status, "LALFindChirpBCVCFilter", FINDCHIRPBCVCFILTERC );
  ATTATCHSTATUSPTR( status );

  /*
   *
   * check that the arguments are reasonable
   *
   */

  /* make sure the output handle exists, but points to a null pointer */
  ASSERT( eventList, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( !*eventList, status, FINDCHIRPH_ENNUL, FINDCHIRPH_MSGENNUL );

  /* make sure that the parameter structure exists */
  ASSERT( params, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the filter parameters are reasonable */
  ASSERT( params->deltaT > 0, status,
      FINDCHIRPH_EDTZO, FINDCHIRPH_MSGEDTZO );
  ASSERT( params->rhosqThresh >= 0, status,
      FINDCHIRPH_ERHOT, FINDCHIRPH_MSGERHOT );
  ASSERT( params->chisqThresh >= 0, status,
      FINDCHIRPH_ECHIT, FINDCHIRPH_MSGECHIT );

  /* check that the fft plan exists */
  ASSERT( params->invPlan, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the workspace vectors exist */
  ASSERT(params->qVec, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qVec->data, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qtildeVec, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qtildeVec->data,status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL);
  ASSERT(params->qVecBCV, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qVecBCV->data, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qtildeVecBCV, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qtildeVecBCV->data,status, FINDCHIRPH_ENULL,
      FINDCHIRPH_MSGENULL);

  /* check that the chisq parameter and input structures exist */
  ASSERT( params->chisqParams, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( params->chisqInput,   status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( params->chisqInputBCV,status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* if a rhosqVec vector has been created, check we can store data in it */
  if ( params->rhosqVec )
  {
    ASSERT( params->rhosqVec->data->data, status,
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
    ASSERT( params->rhosqVec->data, status,
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  }

  /* if a chisqVec vector has been created, check we can store data in it */
  if ( params->chisqVec )
    {
      ASSERT( params->chisqVec->data, status,
          FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  }

  /* make sure that the input structure exists */
  ASSERT( input, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* make sure that the input structure contains some input */
  ASSERT( input->fcTmplt, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( input->segment, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* make sure the filter has been initialized for the correct approximant */
  if ( params->approximant != BCV )
  {
    ABORT( status, FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
  }

  /* make sure that the template and the segment are both BCVC */
  ASSERT( input->fcTmplt->tmplt.approximant == BCV, status,
      FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
  ASSERT( input->segment->approximant == BCV, status,
      FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );


  /*
   *
   * point local pointers to input and output pointers
   *
   */


  /* workspace vectors */
  q    = params->qVec->data;
  qBCV = params->qVecBCV->data;
  qtilde    = params->qtildeVec->data;
  qtildeBCV = params->qtildeVecBCV->data;


  /* template and data */
  inputData     = input->segment->data->data->data;
  inputDataBCV  = input->segment->dataBCV->data->data;
  tmpltSignal   = input->fcTmplt->data->data;
  templateNorm  = input->fcTmplt->tmpltNorm;
  deltaT        = params->deltaT;

  /* some aliases for later use (Thomas) */
  deltaTPower2By3 = pow(params->deltaT, 2./3.);
  deltaTPower1By6 = pow(params->deltaT, 1./6.);

  /* the length of the chisq bin vec is the number of bin   */
  /* _boundaries_ so the number of chisq bins is length - 1 */

  if ( input->segment->chisqBinVec->length )
    {
    /*
     * at this point, numChisqBins is only used as a parameter
     * on the basis of which we decide whether we will do a chisq test or not.
     * the actual number of chisq bins is:
     */
    numChisqBins = input->segment->chisqBinVec->length - 1;
    chisqBin    = input->segment->chisqBinVec->data;
    chisqBinBCV = input->segment->chisqBinVecBCV->data;
    tmpltPower    = input->segment->tmpltPowerVec->data;
    tmpltPowerBCV = input->segment->tmpltPowerVecBCV->data;
  }


  /* number of points in a segment */
  if (params->qVec->length != params->qVecBCV->length)
  {
    ABORT(status, FINDCHIRPBCVH_EQLEN, FINDCHIRPBCVH_MSGEQLEN);
  }
  numPoints = params->qVec->length;

  /*
   * template parameters, since FindChirpBCVCFilterSegment is run
   * for every template
   */
  psi0   = input->fcTmplt->tmplt.psi0;
  psi3   = input->fcTmplt->tmplt.psi3;
  fFinal = input->fcTmplt->tmplt.fFinal;

  {
    /* Calculate deltaEventIndex : the only acceptable clustering */
    /* is "window" method, for BCV                                */
    if ( params->clusterMethod == FindChirpClustering_window )
    {
      deltaEventIndex=(UINT4) rint((params->clusterWindow/params->deltaT)+1.0);
    }
    else if ( params->clusterMethod == FindChirpClustering_tmplt )
    {
      ABORT( status, FINDCHIRPBCVH_ECLUW, FINDCHIRPBCVH_MSGECLUW );
    }


    /* ignore corrupted data at start and end */
    ignoreIndex = ( input->segment->invSpecTrunc / 2 ) + deltaEventIndex;

    if ( lalDebugLevel & LALINFO )
    {
      CHAR newinfomsg[256];


      snprintf( newinfomsg, sizeof(newinfomsg) / sizeof(*newinfomsg),
          "chirp time = %e seconds => %d points\n"
          "invSpecTrunc = %d => ignoreIndex = %d\n",
          chirpTime, deltaEventIndex,
          input->segment->invSpecTrunc, ignoreIndex );
      LALInfo( status, newinfomsg );
    }

    /* XXX check that we are not filtering corrupted data XXX */
    /* XXX this is hardwired to 1/4 segment length        XXX */
    if ( ignoreIndex > numPoints / 4 )
    {
      ABORT( status, FINDCHIRPH_ECRUP, FINDCHIRPH_MSGECRUP );
    }
    /* XXX reset ignoreIndex to one quarter of a segment XXX */
    ignoreIndex = numPoints / 4;
 }

  if ( lalDebugLevel & LALINFO )
  {
    CHAR newinfomsg[256];

    snprintf( newinfomsg, sizeof(newinfomsg) / sizeof(*newinfomsg),
        "filtering from %d to %d\n",
        ignoreIndex, numPoints - ignoreIndex );
    LALInfo( status, newinfomsg );
  }

  /* k that corresponds to fFinal */
  deltaF = 1.0 / ( (REAL4) params->deltaT * (REAL4) numPoints );
  kFinal = fFinal / deltaF < numPoints/2 ? floor(fFinal / deltaF)
    : floor(numPoints/2);
  myfmin = input->segment->fLow;

  /* assign the values to a1, b1 and b2 */
  a1 = input->segment->a1->data[kFinal];
  b1 = input->segment->b1->data[kFinal];
  b2 = input->segment->b2->data[kFinal];


  /*
   * if one or more of a1, b1 and b2 are 0, output a message
   */

  if ( !fabs(a1) || !fabs(b1) || !fabs(b2) )
  {
    if ( lalDebugLevel & LALINFO )
    {
       CHAR newinfomsg[256];
       snprintf( newinfomsg, sizeof(newinfomsg) / sizeof(*newinfomsg),
              "a1 = %e b1 = %e b2 = %e\n"
              "fFinal = %e deltaF = %e numPoints = %d => kFinal = %d\n",
               a1, b1, b2, fFinal, deltaF, numPoints, kFinal );
       LALInfo( status, newinfomsg );
    }
  }

  /*
   *
   * compute qtilde, qtildeBCVC, and q, qBCVC
   * using the correct combination of inputData and inputDataBCVC
   *
   */


  memset( qtilde,    0, numPoints * sizeof(COMPLEX8) );
  memset( qtildeBCV, 0, numPoints * sizeof(COMPLEX8) );

  /* qtilde positive frequency, not DC or nyquist */
  for ( k = 1; k < numPoints/2; ++k )
  {
    REAL4 r    = a1 * inputData[k].re;
    REAL4 s    = a1 * inputData[k].im;
    REAL4 rBCV = b1 * inputData[k].re + b2 * inputDataBCV[k].re;
    REAL4 sBCV = b1 * inputData[k].im + b2 * inputDataBCV[k].im;
    REAL4 x = tmpltSignal[k].re;
    REAL4 y = 0.0 - tmpltSignal[k].im; /* note complex conjugate */

    qtilde[k].re = r * x - s * y ;
    qtilde[k].im = r * y + s * x ;
    qtildeBCV[k].re = rBCV * x - sBCV * y ;
    qtildeBCV[k].im = rBCV * y + sBCV * x ;
  }

  /* inverse fft to get q, and qBCV */
  LALCOMPLEX8VectorFFT( status->statusPtr, params->qVec,
      params->qtildeVec, params->invPlan );
  CHECKSTATUSPTR( status );
  LALCOMPLEX8VectorFFT( status->statusPtr, params->qVecBCV,
      params->qtildeVecBCV, params->invPlan );
  CHECKSTATUSPTR( status );



  /*
   *
   * calculate signal to noise squared
   *
   */



  /* if full snrsq vector is required, set it to zero */
  if ( params->rhosqVec )
    memset( params->rhosqVec->data->data, 0, numPoints * sizeof( REAL4 ) );

  rhosqThresh = params->rhosqThresh;

  norm = deltaT / ((REAL4) numPoints) ;
  /* notice difference from corresponding factor in the sp templates: */
  /* no factor of 4 (taken care of in inputData and inputDataBCV      */
  /* and no segnorm, since we already multiplied by a1, b1 and b2.    */


  /* normalized snr threhold */
  modqsqThresh = rhosqThresh / norm ;

  /* we threshold on the "modified" chisq threshold computed from       */
  /*   chisqThreshFac = delta^2 * norm / p                              */
  /*   rho^2 = norm * modqsq                                            */
  /*                                                                    */
  /* So we actually threshold on                                        */
  /*                                                                    */
  /*    r^2 < chisqThresh * ( 1 + modqsq * chisqThreshFac )             */
  /*                                                                    */
  /* which is the same as thresholding on                               */
  /*    r^2 < chisqThresh * ( 1 + rho^2 * delta^2 / p )                 */
  /* and since                                                          */
  /*    chisq = p r^2                                                   */
  /* this is equivalent to thresholding on                              */
  /*    chisq < chisqThresh * ( p + rho^2 delta^2 )                     */
  /*                                                                    */
  /* The raw chisq is stored in the database. this quantity is chisq    */
  /* distributed with 2p-2 degrees of freedom.                          */

  mismatch = 1.0 - input->fcTmplt->tmplt.minMatch;
  if ( !numChisqBins )
  {
    chisqThreshFac = norm * mismatch * mismatch / (REAL4) numChisqBins;
  }

  /* if full snrsq vector is required, store the snrsq */
  if ( params->rhosqVec )
  {
    memcpy( params->rhosqVec->name, input->segment->data->name,
        LALNameLength * sizeof(CHAR) );
    memcpy( &(params->rhosqVec->epoch), &(input->segment->data->epoch),
        sizeof(LIGOTimeGPS) );
    params->rhosqVec->deltaT = input->segment->deltaT;

    for ( j = 0; j < numPoints; ++j )
    {
      REAL4 modqsqSP  = q[j].re * q[j].re + q[j].im * q[j].im ;
      REAL4 modqsqBCV = qBCV[j].re * qBCV[j].re + qBCV[j].im * qBCV[j].im ;
      REAL4 ImProd = 2.0 * ( - q[j].re * qBCV[j].im + qBCV[j].re * q[j].im ) ;

      REAL4 newmodqsq = ( 0.5 * sqrt( modqsqSP + modqsqBCV + ImProd ) +
          0.5 * sqrt( modqsqSP + modqsqBCV - ImProd ) ) *
        ( 0.5 * sqrt( modqsqSP + modqsqBCV + ImProd ) +
          0.5 * sqrt( modqsqSP + modqsqBCV - ImProd ) ) ;

      params->rhosqVec->data->data[j] = norm * newmodqsq;
    }
  }

  /*
   * Here starts the actual constraint code.
   * */

  /* First for debugging purpose. dont erase  (Thomas)*/

/*
   {
    FILE  *V0File, *V1File,*V2File,
    *alphaFile,  *rho1File, *rho2File, *rho3File,
    *phaseFile, *phaseFileE, *alphaFileE,  *rhoFile, *thetavFile;

    fprintf(stderr, "a1=%e b1=%e b2 =%e fFinal=%e deltTa=%e\n", a1, b1, b2, fFinal, params->deltaT);

   alphaUnity = pow(fFinal, -2./3.)*pow(params->deltaT,-2./3.);


   thetab   = -(a1 * alphaUnity)/(b2+b1*alphaUnity);
   thetab   = atan(thetab);
   fprintf(stderr, "alphaMax -->thetab = %e\n", thetab);

   V0File=fopen("V0File.dat","w");
   V1File=fopen("V1File.dat","w");
   V2File=fopen("V2File.dat","w");
   rho1File=fopen("rho1File.dat","w");
   rho2File=fopen("rho2File.dat","w");
   rho3File=fopen("rho3File.dat","w");
   alphaFile=fopen("alphaFile.dat","w");
   thetavFile=fopen("thetavFile.dat","w");
   phaseFile=fopen("phaseFile.dat","w");

   phaseFileE=fopen("phaseFileE.dat","w");
   alphaFileE=fopen("alphaFileE.dat","w");

   rhoFile=fopen("rhoFile.dat","w");

   for ( j = 0; j < numPoints; ++j )
     {
       REAL4 K1 = q[j].re;
       REAL4 K2 = qBCV[j].re;
       REAL4 K3 = q[j].im;
       REAL4 K4 = qBCV[j].im;

       REAL4 V0 = q[j].re * q[j].re
	 + qBCV[j].re * qBCV[j].re
	 + q[j].im * q[j].im
	 + qBCV[j].im * qBCV[j].im ;

       REAL4 V1 = q[j].re * q[j].re
	 + q[j].im * q[j].im
	 - qBCV[j].im * qBCV[j].im
	 -  qBCV[j].re * qBCV[j].re;

       REAL4 V2 =  2 * ( q[j].re * qBCV[j].re + qBCV [j].im * q[j].im);

       REAL4 rhosqUnconstraint = 0.5 * (V0+sqrt(V1*V1+V2*V2));

       REAL4 thetav = atan2(V2,V1);

       REAL4 Num1 = K2 + K3 ;
       REAL4 Num2 = K2 - K3 ;
       REAL4 Den1 = K1 - K4;
       REAL4 Den2 = K1 + K4;

       REAL4 InvTan1 = (REAL4) atan2(Num1, Den1);
       REAL4 InvTan2 = (REAL4) atan2(Num2, Den2);

       REAL4 omega = 0.5 * InvTan1 + 0.5 * InvTan2 ;

       fprintf(V0File,"%e\n",V0);
       fprintf(V1File,"%e\n",V1);
       fprintf(V2File,"%e\n",V2);


       fprintf(phaseFileE,"%e\n",0.5 * InvTan1 + 0.5 * InvTan2);


       fprintf(alphaFileE,"%e\n", (- b2 * tan(omega)
				   / ( a1 + b1 * tan(omega)))
	       *pow(params->deltaT*fFinal, 2.0/3.0) );


       fprintf(alphaFile,"%e\n", -(b2 * tan(.5*thetav))
	       / (a1 + b1* tan(.5*thetav))*pow(params->deltaT*fFinal, 2.0/3.0));
       fprintf(phaseFile,"%e\n", thetav);


       fprintf(rho1File,"%e\n",sqrt(rhosqUnconstraint*norm));
       fprintf(rho2File,"%e\n",sqrt(.5*(V0 + V1)*norm));
       fprintf(rho3File,"%e\n",sqrt((V0+V1*cos(2*thetab)+V2*sin(2*thetab))/2.*norm));

       if (thetab >= 0){
	 if ( 0  <= thetav && thetav <= 2 * thetab){
	   rhosqConstraint  = sqrt(0.5 * (V0+sqrt(V1*V1+V2*V2)));
	   fprintf(rhoFile,"%e\n",rhosqConstraint);

	 }
	 else if (thetab-LAL_PI <= thetav && thetav < 0) {
	   rhosqConstraint = sqrt((V0 + V1)/2.);
	   fprintf(rhoFile,"%e\n",rhosqConstraint);
	 }
	 else if( (2*thetab  < thetav && thetav<LAL_PI )
		  || thetab-LAL_PI>thetav){
	   rhosqConstraint =sqrt((V0+V1*cos(2*thetab)+V2*sin(2*thetab))/2.);;
	   fprintf(rhoFile,"%e\n",rhosqConstraint);

	 }
	 else
	   {
	     fprintf(stderr,"must not enter here  thetav = %e thetab=%e\n ", thetav , thetab);
	     exit(0);
	   }
       }
       else{
	 if ( 2*thetab  <= thetav && thetav <= 0){
	   rhosqConstraint  = 0.5 * (V0+sqrt(V1*V1+V2*V2));
	   fprintf(rhoFile,"%e\n",sqrt(norm*rhosqConstraint));
	 }
	 else if (0 < thetav &&  thetav  <= LAL_PI +thetab) {
	   rhosqConstraint = (V0 + V1)/2.;
	   fprintf(rhoFile,"%e\n",sqrt(norm*rhosqConstraint));
	 }
        else if( (-LAL_PI-1e-4  <= thetav && thetav < 2*thetab ) ||
	       (LAL_PI +thetab <= thetav && thetav <= LAL_PI+1e-4))
        {
	   rhosqConstraint =(V0+V1*cos(2*thetab)+V2*sin(2*thetab))/2.;
	   fprintf(rhoFile,"%e\n",sqrt(norm*rhosqConstraint));
	 }
	 else
	   {
	     fprintf(stderr,"must not enter herethetav = %e thetab=%e %e %e %d\n ",thetav , thetab, V1, V2);
	     fprintf(rhoFile,"%e\n",-1);
	   }
       }

     }
   fclose(rhoFile);
   fclose(alphaFileE);
   fclose(thetavFile);

   fclose(V0File);
   fclose(V1File);
   fclose(V2File);

   fclose(alphaFile);
   fclose(phaseFile);
   fclose(phaseFileE);
   fclose(rho1File);
   fclose(rho2File);
   fclose(rho3File);
 }

*/


  /* let us create some aliases which are going to be constantly used in the
   * storage of the results. */
  mchirp = (1.0 / LAL_MTSUN_SI) * LAL_1_PI *
        pow( 3.0 / 128.0 / psi0 , ThreeByFive );
  m =  fabs(psi3) /
    (16.0 * LAL_MTSUN_SI * LAL_PI * LAL_PI * psi0) ;
  eta = 3.0 / (128.0*psi0 *
       pow( (m*LAL_MTSUN_SI*LAL_PI), FiveByThree) );
  SQRTNormA1 = sqrt(norm / a1);


  /* BCV Constraint code (Thomas)*/
  /* first we compute what alphaF=1 corresponds to */
  alphaUnity = pow(fFinal, -2./3.) * pow(params->deltaT,-2./3.);
  /* and what is its corresponding angle thetab */
  thetab   = -(a1 * alphaUnity) / (b2 + b1*alphaUnity);
  thetab   = atan(thetab);

  if ( lalDebugLevel & LALINFO )
    {
      CHAR newinfomsg[256];
      snprintf( newinfomsg, sizeof(newinfomsg) / sizeof(*newinfomsg),
		   "thetab = %e and alphaUnity = %e\n",
		   thetab, alphaUnity);
      LALInfo( status, newinfomsg );
    }

  /* look for an event in the filter output */
  for ( j = ignoreIndex; j < numPoints - ignoreIndex; ++j )
  {
    /* First, we compute the quantities V0, V1 and V2 */
    V0 = q[j].re * q[j].re
      + qBCV[j].re * qBCV[j].re
      + q[j].im * q[j].im
      + qBCV[j].im * qBCV[j].im ;

    V1 = q[j].re * q[j].re
      + q[j].im * q[j].im
      - qBCV[j].im * qBCV[j].im
      -  qBCV[j].re * qBCV[j].re;

    V2 =  2 * ( q[j].re * qBCV[j].re + qBCV [j].im * q[j].im);

    /* and finally the unconstraint SNR which is equivalent to
     * the unconstrained BCV situation */
    rhosqUnconstraint = 0.5*( V0 + sqrt(V1*V1 + V2*V2));

    /* We also get the angle between V1 and V2 vectors used later in the
     * constraint part of BCV filtering. */
    thetav = atan2(V2,V1);

    /* Now, we can compare the angle with thetab. Taking care of the angle is
     * quite important.  */
    if (thetab >= 0){
      if ( 0  <= thetav && thetav <= 2 * thetab){
        rhosqConstraint = rhosqUnconstraint;
      }
      else if (thetab-LAL_PI <= thetav && thetav < 0) {
        rhosqConstraint = (V0 + V1)/2.;
      }
     else if( (2*thetab  < thetav && thetav<=LAL_PI+1e-4 )
		 || (-LAL_PI-1e-4<=thetav && thetav < -LAL_PI+thetab)){
        rhosqConstraint =(V0+V1*cos(2*thetab)+V2*sin(2*thetab))/2.;;
      }
      else
      {
	CHAR newinfomsg[256];
	snprintf( newinfomsg, sizeof(newinfomsg) / sizeof(*newinfomsg),
		     "thetab = %e and thetav = %e\n"
		     "thetav not in the range allowed...V1= %e and V2 = %e\n",
		     thetab, thetav, V1, V2 );
	LALInfo( status, newinfomsg );
        ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
      }
    }
    else{
      if ( 2*thetab  <= thetav && thetav <= 0){
        rhosqConstraint = rhosqUnconstraint;
      }
      else if (0 < thetav &&  thetav  <= LAL_PI+thetab ) {
        rhosqConstraint = (V0 + V1)/2.;
      }
      else if( (-LAL_PI-1e-4  <= thetav && thetav < 2*thetab ) ||
	       (LAL_PI +thetab <= thetav && thetav <= LAL_PI+1e-4))
	{
	  rhosqConstraint =(V0+V1*cos(2*thetab)+V2*sin(2*thetab))/2.;
	}
      else
	{
	  CHAR newinfomsg[256];
	  snprintf( newinfomsg, sizeof(newinfomsg) / sizeof(*newinfomsg),
		       "thetab = %e and thetav = %e\n"
		       "thetav not in the range allowed...V1= %e and V2 = %e\n",
		       thetab, thetav, V1, V2 );
	  LALInfo( status, newinfomsg );
	  ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
	}
    }

    /* If one want to check that the code is equivalent to BCVFilter.c, just
     * uncomment the following line.
     */
    /*rhosqConstraint = rhosqUnconstraint; */

    if ( rhosqConstraint > modqsqThresh )
    {
      /* alpha computation needed*/
      alphaU  =   -(b2 * tan(.5*thetav))
        / (a1 + b1* tan(.5*thetav));

      /* I decided to store both constraint and unconstraint alpha */
      if (alphaU > alphaUnity) {
        alphaC = alphaUnity * deltaTPower2By3 ;
         alphaU*= deltaTPower2By3;
      }
      else if (alphaU < 0  ){
        alphaC = 0;
         alphaU*= deltaTPower2By3;
      }
      else {
         alphaU*= deltaTPower2By3;
         alphaC = alphaU;
      }


      if ( ! *eventList )
      {
        /* store the start of the crossing */
        eventStartIdx = j;

        /* if this is the first event, start the list */
        thisEvent = *eventList = (SnglInspiralTable *)
          LALCalloc( 1, sizeof(SnglInspiralTable) );
        if ( ! thisEvent )
        {
          ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
        }

        /* record the data that we need for the clustering algorithm */
        thisEvent->end_time.gpsSeconds = j;
        thisEvent->snr   = rhosqConstraint;
        thisEvent->alpha = alphaC;
        /* use some variable which are not filled in BCV filterinf to keep
         * track of different intermediate values. */
        thisEvent->tau0  = V0;
        thisEvent->tau2  = V1;
        thisEvent->tau3  = V2;
        thisEvent->tau4  = rhosqUnconstraint;
        thisEvent->tau5  = alphaU;
      }
      else if ( ! params->clusterMethod == FindChirpClustering_none &&
          j <= thisEvent->end_time.gpsSeconds + deltaEventIndex &&
          rhosqConstraint > thisEvent->snr )
      {
        /* if this is the same event, update the maximum */
        thisEvent->end_time.gpsSeconds = j;
        thisEvent->snr = rhosqConstraint;
        thisEvent->alpha = alphaC;
        thisEvent->tau0  = V0;
        thisEvent->tau2  = V1;
        thisEvent->tau3  = V2;
        thisEvent->tau4  = rhosqUnconstraint;
        thisEvent->tau5  = alphaU;
      }
      else if (j > thisEvent->end_time.gpsSeconds + deltaEventIndex ||
          params->clusterMethod == FindChirpClustering_none )
      {
        /* clean up this event */
        SnglInspiralTable *lastEvent;
        INT8               timeNS;
        INT4               timeIndex = thisEvent->end_time.gpsSeconds;

        /* set the event LIGO GPS time of the event */
        timeNS = 1000000000L *
          (INT8) (input->segment->data->epoch.gpsSeconds);
        timeNS += (INT8) (input->segment->data->epoch.gpsNanoSeconds);
        timeNS += (INT8) (1e9 * timeIndex * deltaT);
        thisEvent->end_time.gpsSeconds = (INT4) (timeNS/1000000000L);
        thisEvent->end_time.gpsNanoSeconds = (INT4) (timeNS%1000000000L);
        thisEvent->end_time_gmst = fmod(XLALGreenwichMeanSiderealTime(
            &(thisEvent->end_time)), LAL_TWOPI) * 24.0 / LAL_TWOPI;	/* hours */
        ASSERT( !XLAL_IS_REAL8_FAIL_NAN(thisEvent->end_time_gmst), status, LAL_FAIL_ERR, LAL_FAIL_MSG );

        /* set the impuse time for the event */
        thisEvent->template_duration = (REAL8) chirpTime;

        /* record the ifo and channel name for the event */
        strncpy( thisEvent->ifo, input->segment->data->name,
            2 * sizeof(CHAR) );
        strncpy( thisEvent->channel, input->segment->data->name + 3,
            (LALNameLength - 3) * sizeof(CHAR) );
        thisEvent->impulse_time = thisEvent->end_time;

        /* copy the template into the event */
        thisEvent->psi0   = (REAL4) input->fcTmplt->tmplt.psi0;
        thisEvent->psi3   = (REAL4) input->fcTmplt->tmplt.psi3;

        /* chirp mass in units of M_sun */
        thisEvent->mchirp   = mchirp;
        thisEvent->eta      = eta;
        thisEvent->f_final  = (REAL4) input->fcTmplt->tmplt.fFinal ;

        /* set the type of the template used in the analysis */
        snprintf( thisEvent->search, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
            "FindChirpBCVC" );

        thisEvent->chisq     = 0;
        thisEvent->chisq_dof = 0;

        thisEvent->sigmasq = SQRTNormA1;
        thisEvent->eff_distance =
          input->fcTmplt->tmpltNorm / norm / thisEvent->snr;
        thisEvent->eff_distance = sqrt( thisEvent->eff_distance ) / deltaTPower1By6;

        thisEvent->snr *= norm;
        thisEvent->snr = sqrt( thisEvent->snr );

        thisEvent->tau4 = sqrt( thisEvent->tau4 *  norm );

        /* compute the time since the snr crossing */
        thisEvent->event_duration =  (REAL8) timeIndex - (REAL8) eventStartIdx;
        thisEvent->event_duration *= (REAL8) deltaT;

        /* store the start of the crossing */
        eventStartIdx = j;

        /* allocate memory for the newEvent */
        lastEvent = thisEvent;

        lastEvent->next = thisEvent = (SnglInspiralTable *)
          LALCalloc( 1, sizeof(SnglInspiralTable) );
        if ( ! lastEvent->next )
        {
          ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
        }

        /* stick minimal data into the event */
        thisEvent->end_time.gpsSeconds = j;
        thisEvent->snr = rhosqConstraint;
        thisEvent->alpha = alphaC;
        thisEvent->tau0  = V0;
        thisEvent->tau2  = V1;
        thisEvent->tau3  = V2;
        thisEvent->tau4  = rhosqUnconstraint;
        thisEvent->tau5  = alphaU;
      }
    }
  }



  /*
   *
   * clean up the last event if there is one
   *
   */

 if ( thisEvent )
  {
    INT8           timeNS;
    INT4           timeIndex = thisEvent->end_time.gpsSeconds;

    /* set the event LIGO GPS time of the event */
    timeNS = 1000000000L *
      (INT8) (input->segment->data->epoch.gpsSeconds);
    timeNS += (INT8) (input->segment->data->epoch.gpsNanoSeconds);
    timeNS += (INT8) (1e9 * timeIndex * deltaT);
    thisEvent->end_time.gpsSeconds = (INT4) (timeNS/1000000000L);
    thisEvent->end_time.gpsNanoSeconds = (INT4) (timeNS%1000000000L);
    thisEvent->end_time_gmst = fmod(XLALGreenwichMeanSiderealTime(
        &(thisEvent->end_time)), LAL_TWOPI) * 24.0 / LAL_TWOPI;	/* hours */
    ASSERT( !XLAL_IS_REAL8_FAIL_NAN(thisEvent->end_time_gmst), status, LAL_FAIL_ERR, LAL_FAIL_MSG );

    /* set the impuse time for the event */
    thisEvent->template_duration = (REAL8) chirpTime;

    /* record the ifo name for the event */
    strncpy( thisEvent->ifo, input->segment->data->name,
        2 * sizeof(CHAR) );
    strncpy( thisEvent->channel, input->segment->data->name + 3,
        (LALNameLength - 3) * sizeof(CHAR) );
    thisEvent->impulse_time = thisEvent->end_time;


    /* copy the template into the event */
    thisEvent->psi0   = (REAL4) input->fcTmplt->tmplt.psi0;
    thisEvent->psi3   = (REAL4) input->fcTmplt->tmplt.psi3;

    /* chirp mass in units of M_sun */
    thisEvent->mchirp = mchirp;
    thisEvent->f_final  = (REAL4) input->fcTmplt->tmplt.fFinal;
    thisEvent->eta = eta;



    /* set the type of the template used in the analysis */
    snprintf( thisEvent->search, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
        "FindChirpBCVC" );

    /* set snrsq, chisq, sigma and effDist for this event */
    if ( input->segment->chisqBinVec->length )
    {
      /* we store chisq distributed with 2p - 2 degrees of freedom */
      /* in the database. params->chisqVec->data = r^2 = chisq / p */
      /* so we multiply r^2 by p here to get chisq                 */
      thisEvent->chisq =
        params->chisqVec->data[timeIndex] * (REAL4) numChisqBins;
      thisEvent->chisq_dof =  2 * numChisqBins; /* double for BCV */
    }
    else
    {
      thisEvent->chisq     = 0;
      thisEvent->chisq_dof = 0;
    }
    thisEvent->sigmasq = SQRTNormA1;
    thisEvent->eff_distance =
      input->fcTmplt->tmpltNorm / norm / thisEvent->snr;
    thisEvent->eff_distance = sqrt( thisEvent->eff_distance ) / deltaTPower1By6;

    thisEvent->snr *=  norm ;
    thisEvent->snr = sqrt( thisEvent->snr );

    thisEvent->tau4 = sqrt( thisEvent->tau4 * norm );

    /* compute the time since the snr crossing */
    thisEvent->event_duration = (REAL8) timeIndex - (REAL8) eventStartIdx;
    thisEvent->event_duration *= (REAL8) deltaT;
  }


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

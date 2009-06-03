/*
*  Copyright (C) 2007 Jolien Creighton, David McKechan, B.S. Sathyaprakash, Thomas Cokelaer
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

/*<lalVerbatim file="LALInspiralParseParametersCV">
  Author:  Cokelaer, T.
  $Id: LALInspiralParseParameters.c,v 1.24 2009/04/01 00:45:29 ajith Exp $
  </lalVerbatim>  */



/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralParseParameters.c}}
DOC IN PROGRESS


Module to work with the inspiralTemplate Structure.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralITStructureParseParametersCP}
\input{LALInspiralITStructurePrintCP}
\input{LALInspiralITStructureSetDefaultCP}
\input{LALInspiralITStructureHelpCP}


\idx{LALInspiralITStructureParseParameters}
\idx{LALInspiralITStructurePrint}
\idx{LALInspiralITStructureSetDefault}
\idx{LALInspiralITStructureHelp}


\idx{LALInspiralParseParametersForInjection}




\subsubsection*{Description}
This module is a set of functions to play with the inspiralTemplate structure of the
inspiral package. It allows to set default values to the inspiral strcuture, to parse
parameters from the inspiral strcuture and to print the inspiral structure.

Has to check and finalized...

\begin{itemize}
\item  The \texttt{LALInspiralITSructureParseParameters} function allows
the user to parse string with respect to that structure. Each variable in the
inspiralTemplate structure might be parse with a string like "--(name of the variable)+(value)"
 i.e. \textit{--approximant TaylorT1}. Each argument starts with a double dash character
followed by a key word which is exactly as written in the InspiralTemplate Structure such as
 --order, --mass1, --mass2, --fCutoff ...

Once the string is parsed, the checking function is called.

\item The \texttt{LALInspiralITStructurePrint} function will print on the stdout the value
of the InspiralTemplate structure.
\item The \texttt{LALInspiralITStructureSetDefault} set default values to the variables.
Those values are written in the C-code.
\item \texttt{LALInspiralITStructureHelp}
\end{itemize}


\subsubsection*{Algorithm}
None

\vfill{\footnotesize\input{LALInspiralParseParametersCV}}

</lalLaTeX>  */

#include <lal/LALInspiral.h>
#include <lal/StringInput.h>



#define  INSPIRALTEMPLATE_APPROXIMANT 	TaylorT3
#define  INSPIRALTEMPLATE_ORDER 	twoPN
#define  INSPIRALTEMPLATE_AMPORDER 	twoPN
#define  INSPIRALTEMPLATE_MASS1  	10.
#define  INSPIRALTEMPLATE_MASS2 	10.
#define  INSPIRALTEMPLATE_FCUTOFF 	1000.
#define  INSPIRALTEMPLATE_FLOWER 	40.
#define  INSPIRALTEMPLATE_TSAMPLING 	2048.
#define  INSPIRALTEMPLATE_DISTANCE 	1.   /*MPC*/
#define  INSPIRALTEMPLATE_SIGNALAMPLITUDE 1.
#define  INSPIRALTEMPLATE_STARTPHASE    0.
#define  INSPIRALTEMPLATE_STARTTIME     0.

#define  INSPIRALTEMPLATE_THETA 	0.
#define  INSPIRALTEMPLATE_ZETA2 	0.
#define  INSPIRALTEMPLATE_OMEGAS	0.

#define  INSPIRALTEMPLATE_ALPHA 	0.
#define  INSPIRALTEMPLATE_PSI0 		100000.
#define  INSPIRALTEMPLATE_PSI3 		-1000.

#define  INSPIRALTEMPLATE_ALPHA1 	0.
#define  INSPIRALTEMPLATE_ALPHA2 	0.
#define  INSPIRALTEMPLATE_ALPHA3 	0.
#define  INSPIRALTEMPLATE_ALPHA4 	0.
#define  INSPIRALTEMPLATE_ALPHA5 	0.
#define  INSPIRALTEMPLATE_ALPHA6 	0.
#define  INSPIRALTEMPLATE_BETA 		0.

#define  INSPIRALTEMPLATE_INCLINATION 	0.
#define  INSPIRALTEMPLATE_ECCENTRICITY 	0.
#define  INSPIRALTEMPLATE_ORBITTHETA0 	0.7
#define  INSPIRALTEMPLATE_ORBITPHI0 	0.1
#define  INSPIRALTEMPLATE_SPIN1X        0.1
#define  INSPIRALTEMPLATE_SPIN1Y        0.2
#define  INSPIRALTEMPLATE_SPIN1Z        0.3
#define  INSPIRALTEMPLATE_SPIN2X        0.4
#define  INSPIRALTEMPLATE_SPIN2Y        0.5
#define  INSPIRALTEMPLATE_SPIN2Z        0.6


#define  INSPIRALTEMPLATE_SOURCETHETA 	1.
#define  INSPIRALTEMPLATE_SOURCEPHI 	2.

NRCSID (LALINSPIRALPARSEPARAMETERSC, "");


/* <lalVerbatim file="LALInspiralITStructureParseParametersCP"> */
void LALInspiralITStructureParseParameters(LALStatus *status,
					   UINT4 argc,
					   CHAR **argv,
					   InspiralTemplate *params)
/* </lalVerbatim> */
{
  UINT4 i	= 1;

  INITSTATUS( status, "LALInspiralParseParameters", LALINSPIRALPARSEPARAMETERSC);
  ATTATCHSTATUSPTR( status );

  while(i <argc)
    {
      if (strcmp(argv[i], "--approximant")==0)
	{
	  params->massChoice  = m1Andm2;
	  if (strcmp(argv[++i],"AmpCorPPN")==0){
	    params->approximant = AmpCorPPN; }
          else if (strcmp(argv[i],"TaylorN")==0){
	    params->approximant = TaylorN; }
          else if (strcmp(argv[i],"TaylorEt")==0){
	    params->approximant = TaylorEt; }
          else if (strcmp(argv[i],"TaylorT4")==0){
	    params->approximant = TaylorT4; }
          else if (strcmp(argv[i],"TaylorT1")==0){
	    params->approximant = TaylorT1; }
	  else if (strcmp(argv[i],"TaylorT2")==0){
	    params->approximant = TaylorT2;}
	  else if (strcmp(argv[i],"TaylorT3")==0){
	    params->approximant = TaylorT3;}
	  else if (strcmp(argv[i],"TaylorF1")==0){
	    params->approximant = TaylorF1;}
	  else if (strcmp(argv[i],"TaylorF2")==0){
	    params->approximant = TaylorF2;}
	  else if (strcmp(argv[i],"PadeT1")==0){
	    params->approximant = PadeT1;}
	  else if (strcmp(argv[i],"PadeF1")==0){
	    params->approximant = PadeF1;}
	  else if (strcmp(argv[i],"EOB")==0){
	    params->approximant = EOB;}
	  else if (strcmp(argv[i],"EOBNR")==0){
	    params->approximant = EOBNR;}
	  else if (strcmp(argv[i],"IMRPhenomA")==0){
	    params->approximant = IMRPhenomA;}
	  else if (strcmp(argv[i],"SpinTaylor")==0){
	    params->approximant = SpinTaylor;}
	  else if (strcmp(argv[i],"BCV")==0){
	    params->approximant = BCV;
	    params->massChoice  = psi0Andpsi3;
	  }
	  else if (strcmp(argv[i],"BCVSpin")==0){
	    params->massChoice  = psi0Andpsi3;
	    params->approximant = BCVSpin;
	  }
	  else if (strcmp(argv[i],"SpinTaylorT3")==0){
	    params->approximant = SpinTaylorT3; }
          else if (strcmp(argv[i],"Eccentricity")==0){
	    params->approximant = Eccentricity; }
	  else {fprintf(stderr,"Approximant not found in ParseParameter function\n");} /*is it correct ? */
	}/* SpinTaylor is not available here only for inject package*/
      else if (strcmp(argv[i], "--fcutoff")==0){
	params->fCutoff = atof(argv[++i]);}
      else if (strcmp(argv[i], "--order")==0){
	params->order = atoi(argv[++i]);}
      else if (strcmp(argv[i], "--amp-order")==0){
	params->ampOrder = atoi(argv[++i]);}
      else if (strcmp(argv[i], "--nStartPad")==0){
	params->nStartPad = atoi(argv[++i]);}
      else if (strcmp(argv[i], "--nEndPad")==0){
	params->nEndPad = atoi(argv[++i]);}
      else if (strcmp(argv[i],"--fFinal")==0){
        params->fFinal = atof(argv[++i]);}
      else if (strcmp(argv[i],"--mass1")==0){
	params->mass1 = atof(argv[++i]);}
      else if (strcmp(argv[i],"--mass2")==0){
	params->mass2 = atof(argv[++i]);}
      else if (strcmp(argv[i],"--massChoice")==0){
	if (strcmp(argv[++i],"masses")==0){
	  params->massChoice = m1Andm2;}
	else if (strcmp(argv[i],"psi")==0){
	  params->massChoice = psi0Andpsi3;}
      }
      else if (strcmp(argv[i],"--fLower")==0){
	params->fLower = atof(argv[++i]);}
      else if (strcmp(argv[i],"--tSampling")==0){
	params->tSampling = atof(argv[++i]); }
      else if (strcmp(argv[i],"--distance")==0){
	params->distance = atof(argv[++i]); }
      else if (strcmp(argv[i],"--startPhase")==0){
	params->startPhase = atof(argv[++i]); }
      else if (strcmp(argv[i],"--startTime")==0){
	params->startTime = atof(argv[++i]); }
      else if (strcmp(argv[i],"--eccentricity")==0){
	params->eccentricity = atof(argv[++i]); }
      else if (strcmp(argv[i],"--theta")==0){
	params->Theta = atof(argv[++i]); }
      else if (strcmp(argv[i],"--zeta2")==0){
	params->Zeta2 = atof(argv[++i]); }
      else if (strcmp(argv[i],"--alpha")==0){
	params->alpha = atof(argv[++i]); }
      else if (strcmp(argv[i],"--alpha1")==0){
	params->alpha1 = atof(argv[++i]); }
      else if (strcmp(argv[i],"--alpha2")==0){
	params->alpha2 = atof(argv[++i]); }
      else if (strcmp(argv[i],"--alpha3")==0){
	params->alpha3 = atof(argv[++i]); }
      else if (strcmp(argv[i],"--alpha4")==0){
	params->alpha4 = atof(argv[++i]); }
      else if (strcmp(argv[i],"--alpha5")==0){
	params->alpha5 = atof(argv[++i]); }
      else if (strcmp(argv[i],"--alpha6")==0){
	params->alpha6 = atof(argv[++i]); }
      else if (strcmp(argv[i],"--beta")==0){
	params->beta = atof(argv[++i]); }
      else if (strcmp(argv[i],"--psi0")==0){
	params->psi0 = atof(argv[++i]); }
      else if (strcmp(argv[i],"--psi3")==0){
	params->psi3 = atof(argv[++i]); }
      else if (strcmp(argv[i],"--inclination")==0){
	params->inclination = atof(argv[++i]); }
      else if (strcmp(argv[i],"--sourceTheta")==0){
	params->sourceTheta = atof(argv[++i]); }
      else if (strcmp(argv[i],"--sourcePhi")==0){
	params->sourcePhi = atof(argv[++i]); }
      /* we need also the spin .Has to be checked*/
      else if (strcmp(argv[i],"--spin1x")==0){
	params->spin1[0]= atof(argv[++i]);}
      else if (strcmp(argv[i],"--spin1y")==0){
	params->spin1[1]= atof(argv[++i]); }
      else if (strcmp(argv[i],"--spin1z")==0){
	params->spin1[2]= atof(argv[++i]); }
      else if (strcmp(argv[i],"--spin2x")==0){
	params->spin2[0]= atof(argv[++i]);}
      else if (strcmp(argv[i],"--spin2y")==0){
	  params->spin2[1]= atof(argv[++i]);}
      else if  (strcmp(argv[i],"--spin2z")==0){
	  params->spin2[2]= atof(argv[++i]);}
      i++;
   }


  DETATCHSTATUSPTR(status);
  RETURN(status);
}


/* <lalVerbatim file="LALInspiralITStructurePrintCP"> */
void LALInspiralITStructurePrint(LALStatus *status,
				 InspiralTemplate params)
/* </lalVerbatim> */
{
  /* later */
  printf("# approximant = %-15.12d\n", params.approximant);
  printf("# order       = %-15.12d\n", params.order);
  printf("# ampOrder    = %-15.12d\n", params.ampOrder);
  printf("# mass1       = %-15.12f\n", params.mass1);
  printf("# mass2       = %-15.12f\n", params.mass2);
  printf("# fFinal      = %-15.12f\n", params.fFinal);
  printf("# fCutoff     = %-15.12f\n", params.fCutoff);
  printf("# fLower      = %-15.12f\n", params.fLower);
  printf("# tSampling   = %-15.12f\n", params.tSampling);
  printf("# distance    = %-15.12f\n", params.distance);
  printf("# startPhase  = %-15.12f\n", params.startPhase);
  printf("# startTime   = %-15.12f\n", params.startTime);
  printf("# nStartPad   = %-15.12d\n", params.nStartPad);

  printf("# zeta2       = %-15.12f\n", params.Zeta2);
  printf("# omegaS      = %-15.12f\n", params.OmegaS);

  printf("# alpha       = %-15.12f\n", params.alpha);
  printf("# psi0        = %-15.12f\n", params.psi0);
  printf("# psi3        = %-15.12f\n", params.psi3);
  printf("# alpha1      = %-15.12f\n", params.alpha1);
  printf("# alpha2      = %-15.12f\n", params.alpha2);
  printf("# alpha3      = %-15.12f\n", params.alpha3);
  printf("# alpha4      = %-15.12f\n", params.alpha4);
  printf("# alpha5      = %-15.12f\n", params.alpha5);
  printf("# alpha6      = %-15.12f\n", params.alpha6);
  printf("# beta        = %-15.12f\n", params.beta);

  printf("# inclination = %-15.12f\n", params.inclination);
  printf("# sourceTheta = %-15.12f\n", params.sourceTheta);
  printf("# sourcePhi   = %-15.12f\n", params.sourcePhi);
  printf("# spin1 x     = %-15.12f\n", params.spin1[0]);
  printf("# spin1 y     = %-15.12f\n", params.spin1[1]);
  printf("# spin1 z     = %-15.12f\n", params.spin1[2]);
  printf("# spin2 x     = %-15.12f\n", params.spin2[0]);
  printf("# spin2 y     = %-15.12f\n", params.spin2[1]);
  printf("# spin2 z     = %-15.12f\n", params.spin2[2]);


  printf("# eccentricity= %-15.12f\n", params.eccentricity);



/* Paramters which are computed using LALInspiralParameterCalc */

  printf("# chirpMass   = %-15.12f\n", params.chirpMass);
  printf("# eta         = %-15.12f\n", params.eta);
  printf("# totalMass   = %-15.12f\n", params.totalMass);
  printf("# fFinal      = %-15.12f\n", params.fFinal);
  printf("# t0          = %-15.12f\n", params.t0);
  printf("# t2          = %-15.12f\n", params.t2);
  printf("# t3          = %-15.12f\n", params.t3);
  printf("# t4          = %-15.12f\n", params.t4);
  printf("# t5          = %-15.12f\n", params.t5);
  printf("# tC          = %-15.12f\n", params.tC);

  printf("# massChoice  = %-15.12d\n", params.massChoice);

  RETURN(status);
}


/**/
/* <lalVerbatim file="LALInspiralITStructureSetDefaultCP"> */
void LALInspiralITStructureSetDefault(LALStatus *status,
				      InspiralTemplate *params)
/* </lalVerbatim> */
{

  params->approximant  = INSPIRALTEMPLATE_APPROXIMANT;
  params->order        = INSPIRALTEMPLATE_ORDER ;
  params->ampOrder     = INSPIRALTEMPLATE_AMPORDER ;
  params->mass1        = INSPIRALTEMPLATE_MASS1 ;
  params->mass2        = INSPIRALTEMPLATE_MASS2 ;
  params->fCutoff      = INSPIRALTEMPLATE_FCUTOFF;
  params->fLower       = INSPIRALTEMPLATE_FLOWER;
  params->tSampling    = INSPIRALTEMPLATE_TSAMPLING;
  params->distance     = INSPIRALTEMPLATE_DISTANCE;

  params->signalAmplitude = INSPIRALTEMPLATE_SIGNALAMPLITUDE;
  params->startPhase   = INSPIRALTEMPLATE_STARTPHASE;
  params->startTime    = INSPIRALTEMPLATE_STARTTIME;

  params->Theta        = INSPIRALTEMPLATE_THETA;
  params->Zeta2        = INSPIRALTEMPLATE_ZETA2;


  params->alpha        = INSPIRALTEMPLATE_ALPHA;
  params->psi0         = INSPIRALTEMPLATE_PSI0;
  params->psi3         = INSPIRALTEMPLATE_PSI3;

  params->alpha1       = INSPIRALTEMPLATE_ALPHA1;
  params->alpha2       = INSPIRALTEMPLATE_ALPHA2;
  params->alpha3       = INSPIRALTEMPLATE_ALPHA3;
  params->alpha4       = INSPIRALTEMPLATE_ALPHA4;
  params->alpha5       = INSPIRALTEMPLATE_ALPHA5;
  params->alpha6       = INSPIRALTEMPLATE_ALPHA6;
  params->beta         = INSPIRALTEMPLATE_BETA;

  params->inclination  = INSPIRALTEMPLATE_INCLINATION;
  params->sourceTheta  = INSPIRALTEMPLATE_SOURCETHETA;
  params->sourcePhi    = INSPIRALTEMPLATE_SOURCEPHI;
  params->spin1[0]     = INSPIRALTEMPLATE_SPIN1X;
  params->spin1[1]     = INSPIRALTEMPLATE_SPIN1Y;
  params->spin1[2]     = INSPIRALTEMPLATE_SPIN1Z;
  params->spin2[0]     = INSPIRALTEMPLATE_SPIN2X;
  params->spin2[1]     = INSPIRALTEMPLATE_SPIN2Y;
  params->spin2[2]     = INSPIRALTEMPLATE_SPIN2Z;

  params->eccentricity = INSPIRALTEMPLATE_ECCENTRICITY;

  params->ieta         = 1.;
  params->OmegaS       = INSPIRALTEMPLATE_OMEGAS;
  params->nStartPad    = 0;
  params->nEndPad      = 0;
  params->massChoice   = m1Andm2;

  RETURN(status);
}


/* <lalVerbatim file="LALInspiralITStructureHelpCP"> */
void LALInspiralITStructureHelp()
/* </lalVerbatim> */
{

  fprintf(stderr,"InspiralTemplate Structure; parsing arguments\n");
  fprintf(stderr,"--approximant (TaylorT1, TaylorT2, TaylorT3, EOB, BCV, BCVSpin,\nPadeT1, TaylorEt, TaylorT4, TaylorN, AmpCorPPN)\n");
  fprintf(stderr,"--order       (0, 1, 2, 3, 4, 5, 6 (i.e. 4==twoPN)\n");
  fprintf(stderr,"--ampOrder    (0, 1, 2, 3, 4, 5 (i.e. 4==twoPN)\n");
  fprintf(stderr,"--mass1       (in solar mass)\n");
  fprintf(stderr,"--mass2       (in solar mass)\n");
  fprintf(stderr,"--fLower      \n");
  fprintf(stderr,"--tSampling   \n");
  fprintf(stderr,"--distance    \n");
  fprintf(stderr,"--startPhase  \n");
  fprintf(stderr,"--startTime   \n");
  fprintf(stderr,"--fCutoff     \n");

  fprintf(stderr,"--zeta2       (EOB 3PN)\n");
  fprintf(stderr,"--omegaS      (EOB 3PN)\n");

  fprintf(stderr,"--alpha       (BCV must be > 0)\n");
  fprintf(stderr,"--psi0        (BCV must be > 0)\n");
  fprintf(stderr,"--psi3        (BCV must be < 0)\n");
  fprintf(stderr,"--alpha1      (BCVSpin)\n");
  fprintf(stderr,"--alpha2      (BCVSpin)\n");
  fprintf(stderr,"--alpha3      (BCVSpin)\n");
  fprintf(stderr,"--alpha4      (BCVSpin)\n");
  fprintf(stderr,"--alpha5      (BCVSpin)\n");
  fprintf(stderr,"--alpha6      (BCVSpin)\n");
  fprintf(stderr,"--beta        (BCVSpin)\n");

  fprintf(stderr,"--eccentricity \n");

  fprintf(stderr,"--inclination \n");
  fprintf(stderr,"--sourceTheta \n");
  fprintf(stderr,"--sourcePhi   \n");

  fprintf(stderr,"--eccentricity\n");

}

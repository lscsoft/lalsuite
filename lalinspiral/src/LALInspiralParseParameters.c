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

/**
  \author  Cokelaer, T.
  \file
  \ingroup LALInspiral_h

  \brief Module to work with the ::InspiralTemplate Structure.

\heading{Description}
This module is a set of functions to play with the ::InspiralTemplate structure of the
inspiral package. It allows to set default values to the inspiral strcuture, to parse
parameters from the inspiral strcuture and to print the inspiral structure.

Has to check and finalized...

<ul>
<li>  The \c LALInspiralITStructureParseParameters() function allows
the user to parse string with respect to that structure. Each variable in the
inspiralTemplate structure might be parse with a string like "--(name of the variable)+(value)"
 i.e. <em>--approximant TaylorT1</em>. Each argument starts with a double dash character
followed by a key word which is exactly as written in the InspiralTemplate Structure such as
 --order, --mass1, --mass2, --fCutoff ...

Once the string is parsed, the checking function is called.</li>

<li> The \c LALInspiralITStructurePrint() function will print on the stdout the value
of the InspiralTemplate structure.</li>
<li> The \c LALInspiralITStructureSetDefault() set default values to the variables.
Those values are written in the C-code.</li>
<li> \c LALInspiralITStructureHelp()</li>
</ul>


\heading{Algorithm}
None

*/

#include <lal/LALInspiral.h>
#include <lal/StringInput.h>



#define  INSPIRALTEMPLATE_APPROXIMANT           TaylorT4
#define  INSPIRALTEMPLATE_ORDER                 LAL_PNORDER_THREE_POINT_FIVE
#define  INSPIRALTEMPLATE_AMPORDER              LAL_PNORDER_NEWTONIAN
#define  INSPIRALTEMPLATE_MASS1                 10.
#define  INSPIRALTEMPLATE_MASS2                 10.
#define  INSPIRALTEMPLATE_FCUTOFF               1000.
#define  INSPIRALTEMPLATE_FLOWER                40.
#define  INSPIRALTEMPLATE_TSAMPLING             2048.
#define  INSPIRALTEMPLATE_DISTANCE              1.   /*MPC*/
#define  INSPIRALTEMPLATE_SIGNALAMPLITUDE       1.
#define  INSPIRALTEMPLATE_STARTPHASE            0.
#define  INSPIRALTEMPLATE_STARTTIME             0.

#define  INSPIRALTEMPLATE_THETA                 0.
#define  INSPIRALTEMPLATE_ZETA2                 0.
#define  INSPIRALTEMPLATE_OMEGAS                0.

#define  INSPIRALTEMPLATE_ALPHA                 0.
#define  INSPIRALTEMPLATE_PSI0                  100000.
#define  INSPIRALTEMPLATE_PSI3                  -1000.

#define  INSPIRALTEMPLATE_ALPHA1                0.
#define  INSPIRALTEMPLATE_ALPHA2                0.
#define  INSPIRALTEMPLATE_ALPHA3                0.
#define  INSPIRALTEMPLATE_ALPHA4                0.
#define  INSPIRALTEMPLATE_ALPHA5                0.
#define  INSPIRALTEMPLATE_ALPHA6                0.
#define  INSPIRALTEMPLATE_BETA                  0.

#define  INSPIRALTEMPLATE_INCLINATION           0.1
#define  INSPIRALTEMPLATE_ECCENTRICITY          0.
#define  INSPIRALTEMPLATE_ORBITTHETA0           0.0
#define  INSPIRALTEMPLATE_ORBITPHI0             0.0
#define  INSPIRALTEMPLATE_SPIN1X                0.0
#define  INSPIRALTEMPLATE_SPIN1Y                0.0
#define  INSPIRALTEMPLATE_SPIN1Z                0.0
#define  INSPIRALTEMPLATE_SPIN2X                0.0
#define  INSPIRALTEMPLATE_SPIN2Y                0.0
#define  INSPIRALTEMPLATE_SPIN2Z                0.0

#define  INSPIRALTEMPLATE_CHI                   0.
#define  INSPIRALTEMPLATE_KAPPA                 0.

#define  INSPIRALTEMPLATE_SOURCETHETA           0.
#define  INSPIRALTEMPLATE_SOURCEPHI             0.
#define  INSPIRALTEMPLATE_POLARISATIONANGLE     0.

#define INSPIRALTEMPLATE_INTERACTION	   LAL_SIM_INSPIRAL_INTERACTION_ALL_SPIN
#define INSPIRALTEMPLATE_AXISCHOICE             LAL_SIM_INSPIRAL_FRAME_AXIS_VIEW
#define INSPIRALTEMPLATE_FIXEDSTEP              0
#define INSPIRALTEMPLATE_INSPIRALONLY           0

void LALInspiralITStructureParseParameters(LALStatus *status,
					   UINT4 argc,
					   CHAR **argv,
					   InspiralTemplate *params)

{
  UINT4 i	= 1;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  while(i <argc)
    {
      if (strcmp(argv[i],"--approximant")==0)
	  {
	  params->massChoice  = m1Andm2;
	  if (strcmp(argv[++i],"AmpCorPPN")==0){
	    params->approximant = AmpCorPPN; }
	  else if (strcmp(argv[i],"GeneratePPN")==0){
	    params->approximant = GeneratePPN; }
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
	  else if (strcmp(argv[i],"TaylorF2RedSpin")==0){
	    params->approximant = TaylorF2RedSpin;}
	  else if (strcmp(argv[i],"PadeT1")==0){
	    params->approximant = PadeT1;}
	  else if (strcmp(argv[i],"PadeF1")==0){
	    params->approximant = PadeF1;}
	  else if (strcmp(argv[i],"EOB")==0){
	    params->approximant = EOB;}
	  else if (strcmp(argv[i],"EOBNR")==0){
	    params->approximant = EOBNR;}
          else if (strcmp(argv[i],"EOBNRv2")==0){
            params->approximant = EOBNRv2;}
          else if (strcmp(argv[i],"EOBNRv2HM")==0){
            params->approximant = EOBNRv2HM;}
	  else if (strcmp(argv[i],"IMRPhenomA")==0){
	    params->approximant = IMRPhenomA;}
	  else if (strcmp(argv[i],"IMRPhenomB")==0){
	    params->approximant = IMRPhenomB;}
	  else if (strcmp(argv[i],"IMRPhenomFA")==0){
	    params->approximant = IMRPhenomFA;}
	  else if (strcmp(argv[i],"IMRPhenomFB")==0){
	    params->approximant = IMRPhenomFB;}
	  else if (strcmp(argv[i],"SpinTaylor")==0){
	    params->approximant = SpinTaylor;}
	  else if (strcmp(argv[i],"SpinTaylorFrameless")==0){
	    params->approximant = SpinTaylorFrameless;}
	  else if (strcmp(argv[i],"SpinTaylorT3")==0){
	    params->approximant = SpinTaylorT3; }
	  else if (strcmp(argv[i],"SpinTaylorT4")==0){
	    params->approximant = SpinTaylorT4; }
	  else if (strcmp(argv[i],"SpinQuadTaylor")==0){
	    params->approximant = SpinQuadTaylor;}
	  else if (strcmp(argv[i],"PhenSpinTaylorRD")==0){
	    params->approximant = PhenSpinTaylorRD;}
	  else if (strcmp(argv[i],"FindChirpSP")==0){
	    params->approximant = FindChirpSP;}
	  else if (strcmp(argv[i],"FindChirpPTF")==0){
	    params->approximant = FindChirpPTF;}
	  else if (strcmp(argv[i],"BCV")==0){
	    params->approximant = BCV;
	    params->massChoice  = psi0Andpsi3;
	  }
	  else if (strcmp(argv[i],"BCVSpin")==0){
	    params->massChoice  = psi0Andpsi3;
	    params->approximant = BCVSpin;
	  }
	  else if (strcmp(argv[i],"BCVC")==0){
	    params->approximant = BCVC;
	    params->massChoice  = psi0Andpsi3;
	  }
	  else if (strcmp(argv[i],"Eccentricity")==0){
	    params->approximant = Eccentricity; }
	  else if (strcmp(argv[i],"NumRel")==0){
	    params->approximant = NumRel;}
	  else if (strcmp(argv[i],"FrameFile")==0){
	    params->approximant = FrameFile;}
	  else {fprintf(stderr,"Approximant not found in ParseParameter function\n");}
	  }
      else if (strcmp(argv[i],"--order")==0){
	params->order = (LALPNOrder) atoi(argv[++i]);}
      else if (strcmp(argv[i],"--amp-order")==0){
	params->ampOrder = (LALPNOrder) atoi(argv[++i]);}
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
      else if (strcmp(argv[i],"--fcutoff")==0){
	params->fCutoff = atof(argv[++i]);}
      else if (strcmp(argv[i],"--fFinal")==0){
	params->fFinal = atof(argv[++i]);}
      else if (strcmp(argv[i],"--distance")==0){
	params->distance = atof(argv[++i]); }
      else if (strcmp(argv[i],"--startPhase")==0){
	params->startPhase = atof(argv[++i]); }
      else if (strcmp(argv[i],"--startTime")==0){
	params->startTime = atof(argv[++i]); }
      else if (strcmp(argv[i],"--inclination")==0){
	params->inclination = atof(argv[++i]); }
      else if (strcmp(argv[i],"--sourceTheta")==0){
	params->sourceTheta = atof(argv[++i]); }
      else if (strcmp(argv[i],"--sourcePhi")==0){
	params->sourcePhi = atof(argv[++i]); }
      else if (strcmp(argv[i],"--polarisationAngle")==0){
	params->polarisationAngle = atof(argv[++i]); }
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
      else if (strcmp(argv[i],"--tSampling")==0){
	params->tSampling = atof(argv[++i]); }
      else if (strcmp(argv[i],"--nStartPad")==0){
	params->nStartPad = atoi(argv[++i]);}
      else if (strcmp(argv[i],"--nEndPad")==0){
	params->nEndPad = atoi(argv[++i]);}
      else if (strcmp(argv[i],"--eccentricity")==0){
	params->eccentricity = atof(argv[++i]); }
      else if (strcmp(argv[i],"--ieta")==0){
	params->ieta = atof(argv[++i]); }
      else if (strcmp(argv[i],"--Theta")==0){
	params->Theta = atof(argv[++i]); }
      else if (strcmp(argv[i],"--Zeta2")==0){
	params->Zeta2 = atof(argv[++i]); }
      else if (strcmp(argv[i],"--OmegaS")==0){
	params->OmegaS = atof(argv[++i]); }
      else if (strcmp(argv[i],"--alpha")==0){
	params->alpha = atof(argv[++i]); }
      else if (strcmp(argv[i],"--psi0")==0){
	params->psi0 = atof(argv[++i]); }
      else if (strcmp(argv[i],"--psi3")==0){
	params->psi3 = atof(argv[++i]); }
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
      else if (strcmp(argv[i],"--orbitTheta0")==0){
	params->orbitTheta0 = atof(argv[++i]); }
      else if (strcmp(argv[i],"--orbitPhi0")==0){
	params->orbitPhi0 = atof(argv[++i]); }
      else if (strcmp(argv[i],"--chi")==0){
	params->chi = atof(argv[++i]); }
      else if (strcmp(argv[i],"--kappa")==0){
	params->kappa = atof(argv[++i]); }
      else if (strcmp(argv[i],"--interaction")==0)
      {
        if (strcmp(argv[++i],"NO")==0){
          params->interaction = LAL_SIM_INSPIRAL_INTERACTION_NONE; }
        else if (strcmp(argv[i],"SO15")==0){
          params->interaction = LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_15PN; }
        else if (strcmp(argv[i],"SS")==0){
          params->interaction = LAL_SIM_INSPIRAL_INTERACTION_SPIN_SPIN_2PN; }
        else if (strcmp(argv[i],"SSself")==0){
          params->interaction = LAL_SIM_INSPIRAL_INTERACTION_SPIN_SPIN_SELF_2PN; }
        else if (strcmp(argv[i],"QM")==0){
          params->interaction = LAL_SIM_INSPIRAL_INTERACTION_QUAD_MONO_2PN; }
        else if (strcmp(argv[i],"SO25")==0){
          params->interaction = LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_25PN; }
        else if (strcmp(argv[i],"SO30")==0){
          params->interaction = LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_3PN; }
        else if (strcmp(argv[i],"AllSpin")==0){
          params->interaction = LAL_SIM_INSPIRAL_INTERACTION_ALL_SPIN; }
        else if (strcmp(argv[i],"Tidal")==0){
          params->interaction = LAL_SIM_INSPIRAL_INTERACTION_TIDAL_5PN; }
        else if (strcmp(argv[i],"Tidal6")==0){
          params->interaction = LAL_SIM_INSPIRAL_INTERACTION_TIDAL_6PN; }
        else if (strcmp(argv[i],"All")==0){
          params->interaction = LAL_SIM_INSPIRAL_INTERACTION_ALL; }
        else
          fprintf(stderr,"Invalid choice of --interaction\n");
      }
      else if (strcmp(argv[i],"--axisChoice")==0)
      {
        if (strcmp(argv[++i],"TotalJ")==0){
          params->axisChoice = LAL_SIM_INSPIRAL_FRAME_AXIS_TOTAL_J; }
        else if (strcmp(argv[i],"View")==0){
          params->axisChoice = LAL_SIM_INSPIRAL_FRAME_AXIS_VIEW; }
        else if (strcmp(argv[i],"OrbitalL")==0){
          params->axisChoice = LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L; }
        else
          fprintf(stderr,"Invalid choice of --axisChoice\n");
      }
      else if (strcmp(argv[i],"--fixedStep")==0){
        params->fixedStep = atoi(argv[++i]); }
      else if (strcmp(argv[i],"--inspiralOnly")==0){
        params->inspiralOnly = atoi(argv[++i]); }
      i++;
   }


  DETATCHSTATUSPTR(status);
  RETURN(status);
}



void LALInspiralITStructurePrint(LALStatus *status,
				 InspiralTemplate params)

{
  /* Print values of all parameters to screen */
  printf("# approximant         = %-15.12d\n", params.approximant);
  printf("# order               = %-15.12d\n", params.order);
  printf("# ampOrder            = %-15.12d\n", params.ampOrder);
  printf("# mass1               = %-15.12f\n", params.mass1);
  printf("# mass2               = %-15.12f\n", params.mass2);
  printf("# fLower              = %-15.12f\n", params.fLower);
  printf("# fCutoff             = %-15.12f\n", params.fCutoff);
  printf("# distance            = %-15.12f\n", params.distance);
  printf("# startPhase          = %-15.12f\n", params.startPhase);
  printf("# startTime           = %-15.12f\n", params.startTime);
  printf("# inclination         = %-15.12f\n", params.inclination);
  printf("# sourceTheta         = %-15.12f\n", params.sourceTheta);
  printf("# sourcePhi           = %-15.12f\n", params.sourcePhi);
  printf("# polarisationAngle   = %-15.12f\n", params.polarisationAngle);
  printf("# spin1 x             = %-15.12f\n", params.spin1[0]);
  printf("# spin1 y             = %-15.12f\n", params.spin1[1]);
  printf("# spin1 z             = %-15.12f\n", params.spin1[2]);
  printf("# spin2 x             = %-15.12f\n", params.spin2[0]);
  printf("# spin2 y             = %-15.12f\n", params.spin2[1]);
  printf("# spin2 z             = %-15.12f\n", params.spin2[2]);
  printf("# tSampling           = %-15.12f\n", params.tSampling);
  printf("# nStartPad           = %-15.12d\n", params.nStartPad);
  printf("# nEndPad             = %-15.12d\n", params.nStartPad);
  printf("# eccentricity        = %-15.12f\n", params.eccentricity);
  printf("# ieta                = %-15.12d\n", params.ieta);
  printf("# Theta               = %-15.12f\n", params.Theta);
  printf("# Zeta2               = %-15.12f\n", params.Zeta2);
  printf("# OmegaS              = %-15.12f\n", params.OmegaS);
  printf("# alpha               = %-15.12f\n", params.alpha);
  printf("# psi0                = %-15.12f\n", params.psi0);
  printf("# psi3                = %-15.12f\n", params.psi3);
  printf("# alpha1              = %-15.12f\n", params.alpha1);
  printf("# alpha2              = %-15.12f\n", params.alpha2);
  printf("# alpha3              = %-15.12f\n", params.alpha3);
  printf("# alpha4              = %-15.12f\n", params.alpha4);
  printf("# alpha5              = %-15.12f\n", params.alpha5);
  printf("# alpha6              = %-15.12f\n", params.alpha6);
  printf("# beta                = %-15.12f\n", params.beta);
  printf("# chi                 = %-15.12f\n", params.chi);
  printf("# kappa               = %-15.12f\n", params.kappa);
  printf("# orbitTheta0         = %-15.12f\n", params.orbitTheta0);
  printf("# orbitPhi0           = %-15.12f\n", params.orbitPhi0);
  printf("# interaction         = %-15.12d\n", params.interaction);
  printf("# axisChoice          = %-15.12d\n", params.axisChoice);
  printf("# fixedStep           = %-15.12d\n", params.fixedStep);
  printf("# inspiralOnly        = %-15.12d\n", params.inspiralOnly);

/* Paramters which are computed using LALInspiralParameterCalc */

  printf("# chirpMass           = %-15.12f\n", params.chirpMass);
  printf("# eta                 = %-15.12f\n", params.eta);
  printf("# totalMass           = %-15.12f\n", params.totalMass);
  printf("# mu                  = %-15.12f\n", params.mu);
  printf("# fFinal              = %-15.12f\n", params.fFinal);
  printf("# t0                  = %-15.12f\n", params.t0);
  printf("# t2                  = %-15.12f\n", params.t2);
  printf("# t3                  = %-15.12f\n", params.t3);
  printf("# t4                  = %-15.12f\n", params.t4);
  printf("# t5                  = %-15.12f\n", params.t5);
  printf("# t6                  = %-15.12f\n", params.t6);
  printf("# t7                  = %-15.12f\n", params.t7);
  printf("# tC                  = %-15.12f\n", params.tC);
  printf("# signalAmplitude     = %-15.12f\n", params.signalAmplitude);
  printf("# vFinal              = %-15.12f\n", params.vFinal);
  printf("# end_time (s)        = %-15.12d\n", params.end_time.gpsSeconds);
  printf("# end_time (ns)       = %-15.12d\n", params.end_time.gpsNanoSeconds);
  printf("# massChoice          = %-15.12d\n", params.massChoice);

  RETURN(status);
}


/**/

void LALInspiralITStructureSetDefault(LALStatus *status,
				      InspiralTemplate *params)

{

  params->approximant                   = INSPIRALTEMPLATE_APPROXIMANT;
  params->order                         = INSPIRALTEMPLATE_ORDER ;
  params->ampOrder                      = INSPIRALTEMPLATE_AMPORDER ;
  params->mass1                         = INSPIRALTEMPLATE_MASS1 ;
  params->mass2                         = INSPIRALTEMPLATE_MASS2 ;
  params->massChoice                    = m1Andm2;
  params->fLower                        = INSPIRALTEMPLATE_FLOWER;
  params->fCutoff                       = INSPIRALTEMPLATE_FCUTOFF;
  params->fFinal                        = INSPIRALTEMPLATE_FCUTOFF;
  params->distance                      = INSPIRALTEMPLATE_DISTANCE;
  params->startPhase                    = INSPIRALTEMPLATE_STARTPHASE;
  params->startTime                     = INSPIRALTEMPLATE_STARTTIME;
  params->inclination                   = INSPIRALTEMPLATE_INCLINATION;
  params->sourceTheta                   = INSPIRALTEMPLATE_SOURCETHETA;
  params->sourcePhi                     = INSPIRALTEMPLATE_SOURCEPHI;
  params->polarisationAngle             = INSPIRALTEMPLATE_POLARISATIONANGLE;
  params->spin1[0]                      = INSPIRALTEMPLATE_SPIN1X;
  params->spin1[1]                      = INSPIRALTEMPLATE_SPIN1Y;
  params->spin1[2]                      = INSPIRALTEMPLATE_SPIN1Z;
  params->spin2[0]                      = INSPIRALTEMPLATE_SPIN2X;
  params->spin2[1]                      = INSPIRALTEMPLATE_SPIN2Y;
  params->spin2[2]                      = INSPIRALTEMPLATE_SPIN2Z;
  params->tSampling                     = INSPIRALTEMPLATE_TSAMPLING;
  params->nStartPad                     = 0;
  params->nEndPad                       = 0;
  params->eccentricity                  = INSPIRALTEMPLATE_ECCENTRICITY;
  params->Theta                         = INSPIRALTEMPLATE_THETA;
  params->Zeta2                         = INSPIRALTEMPLATE_ZETA2;
  params->OmegaS                        = INSPIRALTEMPLATE_OMEGAS;
  params->alpha                         = INSPIRALTEMPLATE_ALPHA;
  params->psi0                          = INSPIRALTEMPLATE_PSI0;
  params->psi3                          = INSPIRALTEMPLATE_PSI3;
  params->alpha1                        = INSPIRALTEMPLATE_ALPHA1;
  params->alpha2                        = INSPIRALTEMPLATE_ALPHA2;
  params->alpha3                        = INSPIRALTEMPLATE_ALPHA3;
  params->alpha4                        = INSPIRALTEMPLATE_ALPHA4;
  params->alpha5                        = INSPIRALTEMPLATE_ALPHA5;
  params->alpha6                        = INSPIRALTEMPLATE_ALPHA6;
  params->beta                          = INSPIRALTEMPLATE_BETA;
  params->orbitTheta0                   = INSPIRALTEMPLATE_ORBITTHETA0;
  params->orbitPhi0                     = INSPIRALTEMPLATE_ORBITPHI0;
  params->chi                           = INSPIRALTEMPLATE_CHI;
  params->kappa                         = INSPIRALTEMPLATE_KAPPA;
  params->interaction               = INSPIRALTEMPLATE_INTERACTION;
  params->axisChoice                    = INSPIRALTEMPLATE_AXISCHOICE;
  params->fixedStep                     = INSPIRALTEMPLATE_FIXEDSTEP;
  params->inspiralOnly                  = INSPIRALTEMPLATE_INSPIRALONLY;
  params->signalAmplitude               = INSPIRALTEMPLATE_SIGNALAMPLITUDE;
  params->ieta                          = 1.;

  RETURN(status);
}



void LALInspiralITStructureHelp()

{

  fprintf(stderr,"InspiralTemplate Structure; parsing arguments\n");
  fprintf(stderr,"--approximant (TaylorT1, TaylorT2, TaylorT3, TaylorT4,\n");
  fprintf(stderr,"\t\tTaylorF1, TaylorF2, TaylorF2RedSpin, PadeT1, PadeF1, BCV, BCVSpin\n");
  fprintf(stderr,"\t\tBCVC, SpinTaylorT3, SpinTaylorFrameless, SpinTaylor,\n");
  fprintf(stderr,"\t\tSpinQuadTaylor, PhenSpinTaylorRD,\n");
  fprintf(stderr,"\t\tFindChirpSP, FindChirpPTF, GeneratePPN, AmpCorPPN,\n");
  fprintf(stderr,"\t\tFrameFile, NumRel, Eccentricity, EOB, EOBNR,\n");
  fprintf(stderr,"\t\tIMRPhenomA, IMRPhenomB, IMRPhenomFA, IMRPhenomFB,\n");
  fprintf(stderr,"\t\tTaylorEt, TaylorN)\n");
  fprintf(stderr,"--order       (0, 1, 2, 3, 4, 5, 6, 7, 8 (i.e. 4==twoPN)\n");
  fprintf(stderr,"--ampOrder    (0, 1, 2, 3, 4, 5 (i.e. 4==twoPN)\n");
  fprintf(stderr,"--mass1       (in solar mass)\n");
  fprintf(stderr,"--mass2       (in solar mass)\n");
  fprintf(stderr,"--massChoice  ('psi' for BCV(Spin) or 'masses' otherwise)\n");
  fprintf(stderr,"--fLower      (in Hz)\n");
  fprintf(stderr,"--fCutoff     (in Hz)\n");
  fprintf(stderr,"--fFinal      (in Hz - not usually needed)\n");
  fprintf(stderr,"--distance    (in Mpc)\n");
  fprintf(stderr,"--startPhase  \n");
  fprintf(stderr,"--startTime   \n");
  fprintf(stderr,"--inclination (angle between L and line of sight in rad.)\n");
  fprintf(stderr,"--sourceTheta (source sky position zenith angle)\n");
  fprintf(stderr,"--sourcePhi   (source sky position azimuth angle)\n");
  fprintf(stderr,"--polarisationAngle   (i.e. 'psi' in antenna patterns)\n");
  fprintf(stderr,"--spin1x      (Vector components for spin for mass1)\n");
  fprintf(stderr,"--spin1y      (Only used by spinning waveforms)\n");
  fprintf(stderr,"--spin1z      (Kerr limit: s1x^2 + s1y^2 + s1z^2 <= 1)\n");
  fprintf(stderr,"--spin2x      (Vector components for spin for mass2)\n");
  fprintf(stderr,"--spin2y      (Only used by spinning waveforms)\n");
  fprintf(stderr,"--spin2z      (Kerr limit: s1x^2 + s1y^2 + s1z^2 <= 1)\n");
  fprintf(stderr,"--tSampling   (sample rate in Hz)\n");
  fprintf(stderr,"--nStartPad   (# samples of zero padding at start)\n");
  fprintf(stderr,"--nEndPad     (# samples of zero padding at end)\n");
  fprintf(stderr,"--eccentricity \n");
  fprintf(stderr,"--ieta        (set to 0 to turn off eta dependence)\n");
  fprintf(stderr,"--Theta       (EOB 3PN - not used anymore)\n");
  fprintf(stderr,"--Zeta2       (EOB 3PN - not used anymore)\n");
  fprintf(stderr,"--OmegaS      (EOB 3PN - not used anymore)\n");
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
  fprintf(stderr,"--chi         (PTF spin magnitude)\n");
  fprintf(stderr,"--kappa       (PTF spin angle cosine)\n");
  fprintf(stderr,"--orbitTheta0 (initial orientation of L - not used)\n");
  fprintf(stderr,"--orbitPhi0   (initial orientation of L - not used)\n");
  fprintf(stderr,"--interaction (used by PhenSpinTaylorRD to control spin effects included)\n");
  fprintf(stderr,"              NO - no spin effects\n");
  fprintf(stderr,"              SO - spin-orbit effects (default)\n");
  fprintf(stderr,"              SS - spin-spin effects\n");
  fprintf(stderr,"              SSself - self spin-spin effects\n");
  fprintf(stderr,"              QM - quadrupole-monopole effects\n");
  fprintf(stderr,"              Allspin - all of the above effects\n");
  fprintf(stderr,"--inputAxis   (used by PhenSpinTaylorRD to set frame z-axis)\n");
  fprintf(stderr,"              TotalJ - z-axis along initial total angular momentum\n");
  fprintf(stderr,"              View - z-axis along line of sight\n");
  fprintf(stderr,"              OrbitalL - z-axis along initial orbital angular momentum\n");
  fprintf(stderr,"--fixedStep   (used by PhenSpinTaylorRD to set integrator)\n");
  fprintf(stderr,"              1     - use fixed step integrator\n");
  fprintf(stderr,"              other - use adaptive step integrator (default)\n");
  fprintf(stderr,"--inspiralOnly     (used by PhenSpinTaylorRD to attach RD or not)\n");
  fprintf(stderr,"              1     - inspiral-only waveform\n");
  fprintf(stderr,"              other - attach ringdown for IMR waveform (default)\n");


}

/**
 * XLAL function to determine adaptive integration flag from a string.  Returns
 * 1 if string contains 'fixedStep', otherwise returns 0 to signal
 * adaptive integration should be used.
 */
int XLALGetAdaptiveIntFromString(const CHAR *inString)
{
  if (strstr(inString, "fixedStep"))
    return 1;
  else
    return 0;
}

/**
 * XLAL function to determine inspiral-only flag from a string.  Returns
 * 1 if string contains 'inspiralOnly', otherwise returns 0 to signal
 * full inspiral-merger-ringdown waveform should be generated.
 */
int XLALGetInspiralOnlyFromString(const CHAR *inString)
{
  if (strstr(inString, "inspiralOnly"))
    return 1;
  else
    return 0;
}

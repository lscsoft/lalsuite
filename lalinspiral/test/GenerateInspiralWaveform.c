/*
*  Copyright (C) 2007 Anand Sengupta, Thomas Cokelaer, Evan Ochsner
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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

/**
 * \author Sathyaprakash, B. S., Cokelaer T.
 * \file
 * \ingroup LALInspiral_h
 *
 * \brief Test routine for wave generation codes.
 *
 * To get some help just type the name of the executable and the option --h
 *
 * Basically, you can provide all the arguments from the InspiralTemplate structure such as
 * --approximant, --order ....
 *
 */


#define LALGENERATEINSPIRALWAVEFORMC_ENORM 0
#define LALGENERATEINSPIRALWAVEFORMC_ESUB  1
#define LALGENERATEINSPIRALWAVEFORMC_EARG  2
#define LALGENERATEINSPIRALWAVEFORMC_EVAL  3
#define LALGENERATEINSPIRALWAVEFORMC_EFILE 4
#define LALGENERATEINSPIRALWAVEFORMC_EMEM  5

#define LALGENERATEINSPIRALWAVEFORMC_MSGENORM "Normal exit"
#define LALGENERATEINSPIRALWAVEFORMC_MSGESUB  "Subroutine failed"
#define LALGENERATEINSPIRALWAVEFORMC_MSGEARG  "Error parsing arguments"
#define LALGENERATEINSPIRALWAVEFORMC_MSGEVAL "Input argument out of valid range"
#define LALGENERATEINSPIRALWAVEFORMC_MSGEFILE "Could not open file"
#define LALGENERATEINSPIRALWAVEFORMC_MSGEMEM  "Out of memory"

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

#define INSPIRALTEMPLATE_INTERACTION	   LAL_INSPIRAL_INTERACTION_ALL_SPIN
#define INSPIRALTEMPLATE_AXISCHOICE             LAL_SIM_INSPIRAL_FRAME_AXIS_VIEW
#define INSPIRALTEMPLATE_FIXEDSTEP              0
#define INSPIRALTEMPLATE_INSPIRALONLY           0

#include <stdio.h>
#include <lal/LALInspiral.h>
#include <lal/LALNoiseModels.h>
#include <lal/RealFFT.h>
#include <lal/AVFactories.h>
#include <lal/Random.h>
#include <lal/GenerateInspiral.h>
#include <lal/StringInput.h>

#define ERROR( code, msg, statement )                                \
do                                                                   \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
         "        %s %s\n", (code), program, __FILE__,       \
         __LINE__, "$Id$", statement ? statement :  \
                 "", (msg) );                                        \
}                                                                    \
while (0)

#define WARNING( statement )                                         \
do                                                                   \
if ( lalDebugLevel & LALWARNING )                                    \
{                                                                    \
  LALPrintError( "Warning[0]: program %s, file %s, line %d, %s\n"    \
         "        %s\n", program, __FILE__, __LINE__,        \
         "$Id$", (statement) );                         \
}                                                                    \
while (0)

#define INFO( statement )                                            \
do                                                                   \
if ( lalDebugLevel & LALINFO )                                       \
{                                                                    \
  LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"       \
         "        %s\n", program, __FILE__, __LINE__,        \
         "$Id$", (statement) );                         \
}                                                                    \
while (0)

#define SUB( func, statusptr )                                       \
do                                                                   \
if ( (func), (statusptr)->statusCode )                               \
{                                                                    \
  ERROR( LALGENERATEINSPIRALWAVEFORMC_ESUB, LALGENERATEINSPIRALWAVEFORMC_MSGESUB,                      \
         "Function call \"" #func "\" failed:" );                    \
  exit( LALGENERATEINSPIRALWAVEFORMC_ESUB );                                          \
}                                                                    \
while (0)

typedef struct{
  INT4 PrintParameters;
  INT4 output;
  INT4 psd_data_file;
  INT4 taper;
  INT4 realImag;
  char tag[200];
  char waveformString[LIGOMETA_WAVEFORM_MAX];
} OtherParamIn;

char *program;

void printf_timeseries (FILE *f1, UINT4 n, REAL4 *signal1, REAL8 delta);
void printf_timeseries2 (UINT4 n, REAL4 *signal1, REAL4 *signal2, REAL8 delta);
void ParseParameters(UINT4 argc, CHAR **argv, OtherParamIn *otherIn);
void LALGenerateInspiralWaveformHelp(void);
void readPSD(REAL8 *psd, REAL4 Df, UINT4 length);
void buildhoft(LALStatus *status, REAL4Vector *wave, 
        InspiralTemplate *params, OtherParamIn *otherIn);


static void LALInspiralITStructureParseParameters(LALStatus *status,
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
          params->interaction = LAL_INSPIRAL_INTERACTION_NONE; }
        else if (strcmp(argv[i],"SO15")==0){
          params->interaction = LAL_INSPIRAL_INTERACTION_SPIN_ORBIT_15PN; }
        else if (strcmp(argv[i],"SS")==0){
          params->interaction = LAL_INSPIRAL_INTERACTION_SPIN_SPIN_2PN; }
        else if (strcmp(argv[i],"SSself")==0){
          params->interaction = LAL_INSPIRAL_INTERACTION_SPIN_SPIN_SELF_2PN; }
        else if (strcmp(argv[i],"QM")==0){
          params->interaction = LAL_INSPIRAL_INTERACTION_QUAD_MONO_2PN; }
        else if (strcmp(argv[i],"SO25")==0){
          params->interaction = LAL_INSPIRAL_INTERACTION_SPIN_ORBIT_25PN; }
        else if (strcmp(argv[i],"SO30")==0){
          params->interaction = LAL_INSPIRAL_INTERACTION_SPIN_ORBIT_3PN; }
        else if (strcmp(argv[i],"AllSpin")==0){
          params->interaction = LAL_INSPIRAL_INTERACTION_ALL_SPIN; }
        else if (strcmp(argv[i],"Tidal")==0){
          params->interaction = LAL_INSPIRAL_INTERACTION_TIDAL_5PN; }
        else if (strcmp(argv[i],"Tidal6")==0){
          params->interaction = LAL_INSPIRAL_INTERACTION_TIDAL_6PN; }
        else if (strcmp(argv[i],"All")==0){
          params->interaction = LAL_INSPIRAL_INTERACTION_ALL; }
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



static void LALInspiralITStructurePrint(LALStatus *status,
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

static void LALInspiralITStructureSetDefault(LALStatus *status,
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



static void LALInspiralITStructureHelp(void)

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


/* --- Main part --- */
int main (int argc , char **argv) {
  REAL4Vector *signal1 = NULL;   /* storing waveforms */
  REAL4Vector *signal2 = NULL;   /* storing waveforms */
  REAL4Vector *dummy = NULL;     /* stores h+ when hx desired */
  REAL8Vector *psd=NULL;
  static LALStatus status;
  InspiralTemplate params;       /* the parameters */
  REAL8 dt, cosI, ampFac, plusFac, crossFac, dist;
  UINT4 n,i;
  UINT4 nby2;
  InspiralInit paramsInit;
  expnCoeffs ak;

  RealFFTPlan *revp =NULL;
  RealFFTPlan *frwd =NULL;
  char name1[256];

  static OtherParamIn otherIn; /* some extra parameters to parse */
  FILE *f1=NULL, *f2=NULL, *f3=NULL, *f4=NULL; /* output file pointers */

  program = *argv;

  /* ---  we start real computation here --- */
  otherIn.PrintParameters = 0; /* by default we don't print the parameters */
  otherIn.output = 0;          /* by default we don't output the waveforms */
  otherIn.realImag = 0;        /* by default output FD waveforms as |h(f)| */
  strncpy(otherIn.tag, "1", sizeof(otherIn.tag));/*default tag for file names*/
  ParseParameters(argc, argv, &otherIn);/* let's parse user parameters */
  SUB( LALInspiralITStructureSetDefault(&status, &params),
       &status);

  SUB( LALInspiralITStructureParseParameters(&status, argc, argv, &params),
       &status); /* parse user inspiral template parameters */


  SUB( LALInspiralInit(&status, &params, &paramsInit), &status);

  /**************************************************************
   * The following are used only when the psd is read from a
   * file as, for example, the s5 psd
   *************************************************************/

  if (otherIn.psd_data_file)
  {
    n = params.tSampling * 8;
    if (paramsInit.nbins < n)
    {
      paramsInit.nbins = n;
    }
    else
    {
      fprintf(stderr, "Length is not enough\n");
      exit(0);
    }
  }


  if(otherIn.output)
  {
    sprintf(name1, "wave-TD-%s.dat", otherIn.tag); f1 = fopen(name1, "w");
    sprintf(name1, "wave-OT-%s.dat", otherIn.tag); f2 = fopen(name1, "w");
    sprintf(name1, "wave-NW-%s.dat", otherIn.tag); f3 = fopen(name1, "w");
    sprintf(name1, "wave-FD-%s.dat", otherIn.tag); f4 = fopen(name1, "w");
  }



  if (otherIn.PrintParameters)
  {
    fprintf(stderr, "the inspiral structure (your parameters) before the call to the waveform generation:\n");
    SUB( LALInspiralITStructurePrint(&status, params),  &status);
  }

  /* force those parameters */
  dt     = 1./params.tSampling;
  n      = paramsInit.nbins;
  nby2   = n/2;

  if( params.approximant  == AmpCorPPN )
  {
    n *= pow(1./(1.+params.ampOrder/2.), -8./3.);
  }

  if (otherIn.PrintParameters)
  {
    fprintf(stderr, "#Testing Inspiral Signal Generation Codes:\n");
    fprintf(stderr, "#Signal n=%d, t0=%e, t2=%e, t3=%e\n", 
        n, params.t0, params.t2, params.t3);
    fprintf(stderr,"#size in bins %d\n",n);
    fprintf(stderr,"#size in seconds %f\n",params.tC);
  }

  SUB( LALSCreateVector(&status, &(signal1), n), &status);
  SUB( LALSCreateVector(&status, &(signal2), n), &status);
  SUB( LALSCreateVector(&status, &(dummy), n), &status);
  frwd = XLALCreateForwardREAL4FFTPlan(n, 0);
  revp = XLALCreateReverseREAL4FFTPlan(n, 0);
  LALDCreateVector(&status, &psd, (nby2+1));

  params.ieta = 1;

  LALInspiralSetup (&status, &ak, &params);


  switch (params.approximant)
  {
    /* These FD approximants generate h(f) and IFFT to get h(t) */
    case PadeF1:
    case TaylorF1:
    case TaylorF2:
    case TaylorF2RedSpin:
    case FindChirpSP:
    case BCV:
    case BCVC:
    case BCVSpin:
    case IMRPhenomFA:
    case IMRPhenomFB:
      switch( otherIn.output )
      {
        default:
        case 0:
        case 1:
          cosI = cos(params.inclination);
          params.distance *= LAL_PC_SI*1e6;
          ampFac  = -2.0 * params.mu * LAL_MRSUN_SI/params.distance;
          plusFac = ampFac * (1.0 + cosI*cosI);
          params.signalAmplitude = plusFac;
          SUB(LALInspiralWave(&status, signal2, &params), &status);
          XLALREAL4VectorFFT(signal1, signal2, revp);
          break;
        case 2:
          cosI = cos(params.inclination);
          params.distance *= LAL_PC_SI*1e6;
          ampFac  = -2.0 * params.mu * LAL_MRSUN_SI/params.distance;
          crossFac = ampFac * 2.0 * cosI;
          params.signalAmplitude = crossFac;
          SUB(LALInspiralWaveTemplates(&status, dummy, signal2, &params), 
              &status);
          XLALREAL4VectorFFT(signal1, signal2, revp);
          break;
        case 3:
          buildhoft(&status, signal2, &params, &otherIn);
          XLALREAL4VectorFFT(signal1, signal2, revp);
          break;
      }
      break;
    /* For default case, print a warning, but continue */
    /* If a new waveform was implemented, it may succeed */
    /* Otherwise, the call to LALInspiral, etc. will still fail */
    default:
      fprintf(stderr, "Warning! approximant appears to be unknown\n");
      /* fallthrough */
    /* These TD approximants generate h(t), h+(t), hx(t) directly */
    case TaylorT1:
    case TaylorT2:
    case TaylorT3:
    case TaylorT4:
    case TaylorEt:
    case TaylorN:
    case EOB:
    case EOBNR:
    case EOBNRv2:
    case EOBNRv2HM:
    case PadeT1:
    case GeneratePPN:
    case AmpCorPPN:
    case SpinTaylorT3:
    case SpinTaylor:
    case SpinTaylorFrameless:
    case SpinQuadTaylor:
    case PhenSpinTaylorRD:
    case IMRPhenomA:
    case IMRPhenomB:
    case Eccentricity:
      switch( otherIn.output )
      {
        default:
        case 0:
        case 1:
          cosI = cos(params.inclination);
          params.distance *= LAL_PC_SI*1e6;
          ampFac  = -2.0 * params.mu * LAL_MRSUN_SI/params.distance;
          plusFac = ampFac * (1.0 + cosI*cosI);
          params.signalAmplitude = plusFac;
          SUB(LALInspiralWave(&status, signal1, &params), &status);
          break;
        case 2:
          cosI = cos(params.inclination);
          params.distance *= LAL_PC_SI*1e6;
          ampFac  = -2.0 * params.mu * LAL_MRSUN_SI/params.distance;
          crossFac = ampFac * 2.0 * cosI;
          params.signalAmplitude = crossFac;
          SUB(LALInspiralWaveTemplates(&status, dummy, signal1, &params), 
              &status);
          break;
        case 3:
          buildhoft(&status, signal1, &params, &otherIn);
          break;
      }
      if(otherIn.taper) /* Taper if requested */
      {
        LALSimInspiralApplyTaper bookends;
        bookends = 0;
        if (otherIn.taper==1) bookends = LAL_SIM_INSPIRAL_TAPER_START;
        if (otherIn.taper==2) bookends = LAL_SIM_INSPIRAL_TAPER_END;
        if (otherIn.taper==3) bookends = LAL_SIM_INSPIRAL_TAPER_STARTEND;
        XLALSimInspiralREAL4WaveTaper(signal1, bookends);
      }
      break;
  }
  if(otherIn.output) printf_timeseries (f1, n, signal1->data, dt);
  {
    REAL8 df, f, hSq, sSq, hMag, sMag, sRe, sIm, rho, rhosq=0, rhoDet=8.;

    XLALREAL4VectorFFT(signal2, signal1, frwd);
    df = 1./(n * dt);
    if (otherIn.psd_data_file)
    {
      readPSD(psd->data, df, nby2);
    }
    else
    {
      LALNoiseSpectralDensity(&status, psd, &LALLIGOIPsd, df);
    }
    for (i=1; i<nby2; i++)
    {
      UINT4 j = n-i;
      f = (double)i*df;
      if (psd->data[i])
      {
        sRe  = signal2->data[i] * dt;
        sIm  = signal2->data[j] * dt;
        hSq  = sRe*sRe + sIm*sIm;
        hMag = sqrt(hSq);
        sSq  = hSq / (psd->data[i]) ;
        sMag = hMag / (psd->data[i]) ;
        if (f>params.fLower)
        {
          if (params.approximant != EOBNR && params.approximant != EOBNRv2
               && params.approximant != EOBNRv2HM && params.fFinal > f)
            rhosq += sSq;
          else
            rhosq += sSq;
        }
        signal2->data[i] = sRe / sqrt(psd->data[i]);
        signal2->data[j] = sIm / sqrt(psd->data[i]);
        if( otherIn.realImag == 1 )
        {
          if(otherIn.output) fprintf(f3, "%e %e %e %e\n", f, 
            sRe/(psd->data[i]), sIm/(psd->data[i]), psd->data[i]);
          if(otherIn.output) fprintf(f4, "%e %e %e\n", f, sRe, sIm);
        }
        else
        {
          if(otherIn.output) fprintf(f3, "%e %e %e\n", 
            f, sMag, psd->data[i]);
          if(otherIn.output) fprintf(f4, "%e %e\n", f, hMag);
        }
      }
      else
      {
        signal2->data[i] = 0.;
        signal2->data[j] = 0.;
        sSq = 0.;
      }
    }
    /* Above, 'rhosq' = \Sum_i | h(f_i) |^2 / Sn(f_i)
     * Therefore, SNR^2 = 4 \int | h(f) |^2 / Sn(f) df ~= 4 * rhosq * df
     * so SNR = rho = sqrt(4 * rhosq * df)
     */
    rho = sqrt(4. * rhosq * df);
    /* Distance in Mpc at which the SNR is 8 */
    dist = (params.distance/(LAL_PC_SI*1e6))*(rho/rhoDet);
    if( otherIn.PrintParameters )
    {
      fprintf(stderr, "mass1: %e, mass2: %e, # samples: %d,\nSNR: %e, " 
          "distance at which SNR=8: %e\n", params.mass1, params.mass2, 
          signal2->length, rho, dist);
    }
    signal2->data[0] = 0.;
    signal2->data[nby2] = 0.;
    XLALREAL4VectorFFT(signal1, signal2, revp);
    if(otherIn.output) printf_timeseries (f2, n, signal1->data, dt);
  }

  XLALDestroyREAL4FFTPlan (revp);
  XLALDestroyREAL4FFTPlan (frwd);
  SUB( LALSDestroyVector(&status, &signal1), &status);
  SUB( LALSDestroyVector(&status, &signal2), &status);
  SUB( LALDDestroyVector(&status, &psd), &status);

  if (otherIn.PrintParameters)
  {
    fprintf(stderr, "fFinal = %f Hz tFinal = %f seconds\n" , 
        params.fFinal, params.tC);
    fprintf(stderr, "the inspiral structure after the call "
        "to the waveform generation:\n");
    SUB( LALInspiralITStructurePrint(&status, params),  &status);
  }
  /*
  LALCheckMemoryLeaks();
  */
  if(otherIn.output)
  {
    fclose(f1);
    fclose(f2);
    fclose(f3);
    fclose(f4);
  }

  return 0;
}

void
ParseParameters( UINT4              argc,
                 CHAR             **argv,
                 OtherParamIn     *otherIn)
{
  UINT4 i = 1, i1 = 0;
  INT4 order = 0;

  while(i < argc)
    {
      if ( strcmp(argv[i],    "--verbose")     == 0 )
      {
        otherIn->PrintParameters = 1;
      }
      else if ( strcmp(argv[i],    "--output")     == 0 )
      {
        otherIn->output = atoi(argv[++i]);
      }
      else if ( strcmp(argv[i],    "--taper")     == 0 )
      {
        otherIn->taper = atoi(argv[++i]);
      }
      else if ( strcmp(argv[i],    "--psd-data-file") == 0 )
      {
        otherIn->psd_data_file = 1;
      }
      else if( strcmp(argv[i],    "--h")     == 0 )
      {
        LALGenerateInspiralWaveformHelp();
      }
      else if( strcmp(argv[i],"-h")     == 0 )
      {
        LALGenerateInspiralWaveformHelp();
      }
      else if( strcmp(argv[i],"-help")     == 0 )
      {
        LALGenerateInspiralWaveformHelp();
      }
      else if( strcmp(argv[i],"--help")    == 0 )
      {
        LALGenerateInspiralWaveformHelp();
      }
      else if( strcmp(argv[i],"--tag") == 0 )
      {
        strncpy(otherIn->tag, argv[++i], sizeof(otherIn->tag) - 1);
      }
      else if( strcmp(argv[i],"--approximant") == 0 )
      {
        i1 = i + 1;
      }
      else if( strcmp(argv[i],"--order") == 0 )
      {
        order = atoi(argv[++i]);
      }
      else if( strcmp(argv[i],"--real-imag") == 0 )
      {
        otherIn->realImag = atoi(argv[++i]);
      }
      i++;
    }
    if( otherIn->output == 3 )
    { /* concatenate approximant and order in a string for SimInspiralTable */
      strcpy( otherIn->waveformString, argv[i1] );
      switch( order )
      {
        case 8:
          strcat( otherIn->waveformString, "pseudoFourPN" );
          break;
        case 7:
          strcat( otherIn->waveformString, "threePointFivePN" );
          break;
        case 6:
          strcat( otherIn->waveformString, "threePN" );
          break;
        case 5:
          strcat( otherIn->waveformString, "twoPointFivePN" );
          break;
        case 4:
          strcat( otherIn->waveformString, "twoPN" );
          break;
        case 3:
          strcat( otherIn->waveformString, "onePointFivePN" );
          break;
        case 2:
          strcat( otherIn->waveformString, "onePN" );
          break;
        case 1:
          strcat( otherIn->waveformString, "oneHalfPN" );
          break;
        case 0:
	  printf(" WARNING: you have chose Newtonian order\n");
          strcat( otherIn->waveformString, "newtonian" );
          break;
        default:
          printf("Invalid PN order requested, using 3.5PN\n");
          strcat( otherIn->waveformString, "threePointFivePN" );
          break;
      }
    }
}


void LALGenerateInspiralWaveformHelp(void)
{

  fprintf(stderr,"GenerateInspiralWaveform Help\n");
  fprintf(stderr,"---------------------------------------------------------\n");
  fprintf(stderr,"All unspecified parameters use default values.\n");
  fprintf(stderr,"If you don't know what a parameter does it's");
  fprintf(stderr," (probably) safe to omit it.\n");
  fprintf(stderr,"---------------------------------------------------------\n");
  fprintf(stderr,"--h for help\n");
  fprintf(stderr,"--verbose to print Inspiral Template parameters\n");
  fprintf(stderr,"--tag TAG adds 'TAG' identifier to file names\n");
  fprintf(stderr,"--psd-data-file file.dat reads noise curve from file.dat (NOT WORKING YET!)\n");
  fprintf(stderr,"\t Initial LIGO design PSD used if no PSD file is given\n");
  fprintf(stderr,"--output=N to generate time and freq. domain waveforms:\n");
  fprintf(stderr,"         N = 0 - do not output any files (default)\n");
  fprintf(stderr,"         N = 1 - output plus polarization h+\n");
  fprintf(stderr,"         N = 2 - output cross polarization hx\n");
  fprintf(stderr,"         N = 3 - output measured strain h = F+ h+ + Fx hx\n");
  fprintf(stderr,"For N != 0, the waveform will be written to these files:\n");
  fprintf(stderr,"         wave-TD-TAG.dat (time domain)\n");
  fprintf(stderr,"         wave-OT-TAG.dat (TD (Wiener) optimal template)\n");
  fprintf(stderr,"         wave-FD-TAG.dat (frequency domain)\n");
  fprintf(stderr,"         wave-NW-TAG.dat (noise-weighted FD)\n");
  fprintf(stderr,"Note: Not all approximants support all types of output\n");
  fprintf(stderr,"         N = 1 calls LALInspiralWave()\n");
  fprintf(stderr,"         N = 2 calls LALInspiralWaveTemplates()\n");
  fprintf(stderr,"         N = 3 calls LALInspiralWaveForInjection() via LALGenerateInspiral()\n");
  fprintf(stderr,"If the approximant you want fails, make it callable from these functions\n");
  fprintf(stderr,"--real-imag=N controls output of complex FD and NW waveforms:\n");
  fprintf(stderr,"         N = 0 - output | h(f) | (default)\n");
  fprintf(stderr,"         N = 1 - output Re(h(f)) and Im(h(f))\n");
  fprintf(stderr,"---------------------------------------------------------\n");
  LALInspiralITStructureHelp();

}


void printf_timeseries (FILE *f1, UINT4 n, REAL4 *signal1, REAL8 delta)
{
  UINT4 i=0;

  /*
  FILE *outfile1;
  outfile1=fopen("wave1.dat","w");
  */
  do
     /* fprintf (outfile1,"%e %e\n", i*delta+t0, *(signal1+i)  );
      */
     fprintf (f1,"%e %e\n", i*delta, signal1[i]  );
  while (n-++i);

  /*
  fprintf(outfile1,"&\n");
  fclose(outfile1);
  */
}



void printf_timeseries2 (UINT4 n, REAL4 *signal1, REAL4 *signal2, REAL8 delta)
{
  UINT4 i=0;
  FILE *outfile1;

  outfile1=fopen("wave2.dat","w");

  do
     /* fprintf (outfile1,"%e %e\n", i*delta+t0, *(signal+i)  );
      */
     fprintf (stdout,"%e %e %e\n", i*delta, *(signal1+i), *(signal2+i)  );
  while (n-++i);

  fclose(outfile1);
}


void readPSD(REAL8 *psd, REAL4 Df, UINT4 length)
{
  FILE *f1;
  UINT4 i=1, j=0;
  REAL4 f;
  REAL4 x;

  psd[0] = 0;
  f1 = fopen("lho4k_070318_strain.txt", "r");

  if (fscanf(f1, "%e %e", &f, &x) == EOF)
  {
    fprintf(stderr, "Problem reading file lho4k_070318_strain.txt stopping\n");
    exit(0);
  }
  if (f!=Df)
  {
    fprintf(stderr, "Incompatible frequency resolution, exiting\n");
    exit(0);
  }

  psd[i] = x*x;

  while(fscanf(f1, "%e %e", &f, &x) != EOF && i<length)
  {
    ++i;
    psd[i] = x*x;
  }

  for (j=i; j<length; j++)
  {
    psd[j] = 0;
  }
  fclose(f1);

  return;
}

void buildhoft(LALStatus *status, REAL4Vector *wave, 
        InspiralTemplate *params, OtherParamIn *otherIn)
{
  REAL4 Fp, Fc, hp, hc, A1, A2, cosshift, sinshift;
  REAL4 cosphi, sinphi, theta, phi, psi;
  UINT4 i, len;
  CoherentGW waveform;
  memset( &waveform, 0, sizeof(CoherentGW) );
  PPNParamStruc ppnParams;
  memset( &ppnParams, 0, sizeof(PPNParamStruc) );
  ppnParams.deltaT   = 1./params->tSampling;
  ppnParams.lengthIn = 0;
  ppnParams.ppn      = NULL;

  /* Construct a SimInspiralTable... */
  SimInspiralTable simTable;
  memset( &simTable, 0, sizeof(SimInspiralTable) );
  /* ... and populate it with desired parameter values */
  /*simTable.waveform = otherIn->waveformString;*/
  memcpy( simTable.waveform, otherIn->waveformString, 
          sizeof(otherIn->waveformString) );

  if (strstr(simTable.waveform,"PhenSpinTaylorRD")) {
    if (params->axisChoice==LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L) {
      strcat(simTable.waveform,"OrbitalL");}
    else if (params->axisChoice==LAL_SIM_INSPIRAL_FRAME_AXIS_TOTAL_J) {
      strcat(simTable.waveform,"TotalJ");
    }
    if (params->inspiralOnly==1) {
      strcat(simTable.waveform,"inspiralOnly");
    }
    if (params->fixedStep==1) {
      strcat(simTable.waveform,"fixedStep");
    }
  }
  switch (params->interaction) {
  case LAL_INSPIRAL_INTERACTION_NONE:
    strcat(simTable.waveform,"NO");
    break;
  case LAL_INSPIRAL_INTERACTION_SPIN_ORBIT_15PN:
    strcat(simTable.waveform,"SO15PN");
    break;
  case LAL_INSPIRAL_INTERACTION_SPIN_ORBIT_25PN:
    strcat(simTable.waveform,"SO25PN");
    break;
  case LAL_INSPIRAL_INTERACTION_SPIN_ORBIT_3PN:
    strcat(simTable.waveform,"SO");
    break;
  case LAL_INSPIRAL_INTERACTION_SPIN_SPIN_2PN:
    strcat(simTable.waveform,"SS");
    break;
  case LAL_INSPIRAL_INTERACTION_SPIN_SPIN_SELF_2PN:
    strcat(simTable.waveform,"SELF");
    break;
  case LAL_INSPIRAL_INTERACTION_QUAD_MONO_2PN:
    strcat(simTable.waveform,"QM");
    break;
  case LAL_INSPIRAL_INTERACTION_ALL_SPIN:
    strcat(simTable.waveform,"ALL_SPIN");
    break;
  case LAL_INSPIRAL_INTERACTION_TIDAL_5PN:
	strcat(simTable.waveform,"TIDAL5PN");
	break;
  case LAL_INSPIRAL_INTERACTION_TIDAL_6PN:
	strcat(simTable.waveform,"TIDAL");
	break;
  case LAL_INSPIRAL_INTERACTION_ALL:
	strcat(simTable.waveform,"ALL");
	break;
		  
  }

  simTable.mass1 = params->mass1;
  simTable.mass2 = params->mass2;
  simTable.eta = params->eta;
  simTable.mchirp = params->chirpMass;
  simTable.distance = params->distance;
  simTable.longitude = params->sourceTheta;
  simTable.latitude = params->sourcePhi;
  simTable.inclination = params->inclination;
  simTable.coa_phase = params->startPhase; /* Is this what I want?? */
  simTable.polarization = params->polarisationAngle;
  simTable.psi0 = params->psi0;
  simTable.psi3 = params->psi3;
  simTable.alpha = params->alpha;
  simTable.alpha1 = params->alpha1;
  simTable.alpha2 = params->alpha2;
  simTable.alpha3 = params->alpha3;
  simTable.alpha4 = params->alpha4;
  simTable.alpha5 = params->alpha5;
  simTable.alpha6 = params->alpha6;
  simTable.beta = params->beta;
  simTable.spin1x = params->spin1[0];
  simTable.spin1y = params->spin1[1];
  simTable.spin1z = params->spin1[2];
  simTable.spin2x = params->spin2[0];
  simTable.spin2y = params->spin2[1];
  simTable.spin2z = params->spin2[2];
  simTable.theta0 = params->orbitTheta0;
  simTable.phi0 = params->orbitPhi0;
  simTable.f_lower = params->fLower;
  simTable.f_final = params->fFinal;
  simTable.amp_order = params->ampOrder;
  sprintf(simTable.taper, "TAPER_NONE");
  /* We'll taper (or not) later in code, so just set TAPER_NONE here */
    
  theta = params->sourceTheta;
  phi   = params->sourcePhi;
  psi   = params->polarisationAngle;
  Fp    = 0.5 * (1. + cos(theta)*cos(theta) ) * cos(2.*phi) * cos(2.*psi)
              - cos(theta) * sin(2.*phi) * sin(2.*psi);
  Fc    = 0.5 * (1. + cos(theta)*cos(theta) ) * cos(2.*phi) * sin(2.*psi)
              + cos(theta) * sin(2.*phi) * cos(2.*psi);

  /* This function fills a CoherentGW with info to construct h(t) */
  LALGenerateInspiral(status, &waveform, &simTable, &ppnParams);

    
  if( waveform.h == NULL ) /* build h(t) from a, phi, shift */
  {
    len = waveform.phi->data->length;
    /* Some approximants do not set waveform.shift. We use shift(t) if set. */
    /* Otherwise, we assume the shift is a constant 0 */
    if( waveform.shift == NULL )
    {
      for(i = 0; i < len; i++)
      {
        A1 = waveform.a->data->data[2*i];
        A2 = waveform.a->data->data[2*i+1];
        cosshift = 1.;
        sinshift = 0.;
        cosphi = cos( waveform.phi->data->data[i] );
        sinphi = sin( waveform.phi->data->data[i] );
        hp = A1 * cosshift * cosphi - A2 * sinshift * sinphi;
        hc = A1 * sinshift * cosphi + A2 * cosshift * sinphi;
        wave->data[i] = Fp * hp + Fc * hc;
      }
    }
    else
    {
      for(i = 0; i < len; i++)
      {
        A1 = waveform.a->data->data[2*i];
        A2 = waveform.a->data->data[2*i+1];
        cosshift = cos( waveform.shift->data->data[i] );
        sinshift = sin( waveform.shift->data->data[i] );
        cosphi = cos( waveform.phi->data->data[i] );
        sinphi = sin( waveform.phi->data->data[i] );
        hp = A1 * cosshift * cosphi - A2 * sinshift * sinphi;
        hc = A1 * sinshift * cosphi + A2 * cosshift * sinphi;
        wave->data[i] = Fp * hp + Fc * hc;
      }
    }
  }
  else /* build h(t) from h+ and hx in waveform->h */
  {
    len=waveform.h->data->length < wave->length ? waveform.h->data->length : wave->length;
    for(i = 0; i < len; i++)
      {
	wave->data[i] = Fp * waveform.h->data->data[2*i] 
	  + Fc * waveform.h->data->data[2*i+1];
      }
    for (i=len;i<wave->length;i++) wave->data[i]=0.;
  }

  return;
}

/*
*  Copyright (C) 2007 Stas Babak, Drew Keppel, Duncan Brown, Eirini Messaritaki, Gareth Jones, Thomas Cokelaer, Laszlo Vereb
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

#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>
#include <lal/GenerateInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/SeqFactories.h>
#include <lal/Units.h>

/** \see See \ref GenerateInspiral_h for documentation */
void
LALGenerateInspiral(
    LALStatus		*status,	/**< UNDOCUMENTED */
    CoherentGW		*waveform,	/**< UNDOCUMENTED */
    SimInspiralTable	*thisEvent,	/**< UNDOCUMENTED */
    PPNParamStruc	*ppnParams	/**< UNDOCUMENTED */
    )

{
  LALPNOrder        order;              /* Order of the model             */
  Approximant       approximant;        /* And its approximant value      */
  InspiralTemplate  inspiralParams;     /* structure for inspiral package */
  CHAR              warnMsg[1024];
  int               oldxlalErrno;       /* store old xlal error number    */

  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  ASSERT(thisEvent, status,
      GENERATEINSPIRALH_ENULL, GENERATEINSPIRALH_MSGENULL);

  /* read the event waveform approximant and order */
  oldxlalErrno = xlalErrno;
  xlalErrno = 0;
  if (XLALGetApproximantFromString(thisEvent->waveform, &approximant) == XLAL_FAILURE)
    ABORTXLAL(status);

  if (XLALGetOrderFromString(thisEvent->waveform, &order) == XLAL_FAILURE)
    ABORTXLAL(status);
  xlalErrno = oldxlalErrno;

  /* when entering here, approximant is in principle well defined.  */
  /* We dont need any else if or ABORT in the if statement.         */
  /* Depending on the apporixmant we use inject or inspiral package */
  if ( approximant == GeneratePPN )
  {
    /* fill structure with input parameters */
    oldxlalErrno = xlalErrno;
    xlalErrno = 0;
    if (XLALGenerateInspiralPopulatePPN(ppnParams, thisEvent) == XLAL_FAILURE)
      ABORTXLAL(status);
    xlalErrno = oldxlalErrno;

    /* generate PPN waveform */
    LALGeneratePPNInspiral(status->statusPtr, waveform, ppnParams);
    CHECKSTATUSPTR(status);
  }
  else if ( approximant == AmpCorPPN )
  {
    int i;

    /* fill structure with input parameters */
    oldxlalErrno = xlalErrno;
    xlalErrno = 0;
    if (XLALGenerateInspiralPopulatePPN(ppnParams, thisEvent) == XLAL_FAILURE)
      ABORTXLAL(status);
    xlalErrno = oldxlalErrno;

    /* PPN parameter. */
    ppnParams->ppn = NULL;
    LALSCreateVector( status->statusPtr, &(ppnParams->ppn), order + 1 );
    ppnParams->ppn->length = order + 1;

    ppnParams->ppn->data[0] = 1.0;
    if ( order > 0 )
    {
      ppnParams->ppn->data[1] = 0.0;
      for ( i = 2; i <= (INT4)( order ); i++ )
      {
        ppnParams->ppn->data[i] = 1.0;
      }
    }

    /* generate PPN waveform */
    LALGeneratePPNAmpCorInspiral(status->statusPtr, waveform, ppnParams);
    CHECKSTATUSPTR(status);

    LALSDestroyVector(status->statusPtr, &(ppnParams->ppn) );
    CHECKSTATUSPTR(status);
  }
  else
  {
    inspiralParams.approximant = approximant;
    inspiralParams.order       = order;
    if ((approximant == SpinQuadTaylor)||(approximant == PhenSpinTaylorRD)) {
		xlalErrno = 0;
		if (XLALGetInteractionFromString(&inspiralParams.interaction, thisEvent->waveform) == XLAL_FAILURE) {
			ABORTXLAL(status);
		}
	}

	if (approximant == PhenSpinTaylorRD) {
	  xlalErrno = 0;
	  if ( (XLALGetAxisChoiceFromString(&inspiralParams.axisChoice, thisEvent->waveform) == XLAL_FAILURE) || 
	       (XLALGetAdaptiveIntFromString(&inspiralParams.fixedStep, thisEvent->waveform) == XLAL_FAILURE) || 
	       (XLALGetInspiralOnlyFromString(&inspiralParams.inspiralOnly, thisEvent->waveform) == XLAL_FAILURE ) ) {
	    ABORTXLAL(status);
	  }
	}

    /* We fill ppnParams */
    oldxlalErrno = xlalErrno;
    xlalErrno = 0;
    if (XLALGenerateInspiralPopulatePPN(ppnParams, thisEvent) == XLAL_FAILURE)
      ABORTXLAL(status);
    xlalErrno = oldxlalErrno;

    /* we fill inspiralParams structure as well.*/
    oldxlalErrno = xlalErrno;
    xlalErrno = 0;
    if (XLALGenerateInspiralPopulateInspiral(&inspiralParams, thisEvent, ppnParams) == XLAL_FAILURE)
      ABORTXLAL(status);
    xlalErrno = oldxlalErrno;

    /* the waveform generation itself */
    LALInspiralWaveForInjection(status->statusPtr, waveform, &inspiralParams,
        ppnParams);
    /* we populate the simInspiral table with the fFinal needed for
       template normalisation. */
    thisEvent->f_final = inspiralParams.fFinal;
    CHECKSTATUSPTR(status);
  }

  /* If no waveform has been generated. (AmpCorPPN and PhenSpinTaylorRD and SpinTaylorFrameless fill waveform.h) */
  if ( waveform->a == NULL && approximant != AmpCorPPN && approximant != PhenSpinTaylorRD && approximant != SpinTaylorFrameless 
       && approximant != EOBNRv2 && approximant != EOBNRv2HM )
  {
    snprintf( warnMsg, sizeof(warnMsg)/sizeof(*warnMsg),
        "No waveform generated (check lower frequency)\n");
    LALInfo( status, warnMsg );
    ABORT( status, LALINSPIRALH_ENOWAVEFORM, LALINSPIRALH_MSGENOWAVEFORM );
  }
  if ( waveform->h == NULL && ( approximant == AmpCorPPN || approximant == PhenSpinTaylorRD || approximant == SpinTaylorFrameless
       || approximant == EOBNRv2 || approximant == EOBNRv2HM ) )
  {
    snprintf( warnMsg, sizeof(warnMsg)/sizeof(*warnMsg),
             "No waveform generated (check lower frequency)\n");
    LALInfo( status, warnMsg );
    ABORT( status, LALINSPIRALH_ENOWAVEFORM, LALINSPIRALH_MSGENOWAVEFORM );
  }


  /* If sampling problem. (AmpCorPPN may not be compatible) */
  if ( ppnParams->dfdt > 2.0 && approximant != AmpCorPPN )
  {
    snprintf( warnMsg, sizeof(warnMsg)/sizeof(*warnMsg),
        "Waveform sampling interval is too large:\n"
        "\tmaximum df*dt = %f", ppnParams->dfdt );
    LALInfo( status, warnMsg );
    ABORT( status, GENERATEINSPIRALH_EDFDT, GENERATEINSPIRALH_MSGEDFDT );
  }

  /* Some info should add everything (spin and so on) */
  snprintf( warnMsg, sizeof(warnMsg)/sizeof(*warnMsg),
      "Injected waveform parameters:\n"
      "ppnParams->mTot\t= %"LAL_REAL4_FORMAT"\n"
      "ppnParams->eta\t= %"LAL_REAL4_FORMAT"\n"
      "ppnParams->d\t= %"LAL_REAL4_FORMAT"\n"
      "ppnParams->inc\t= %"LAL_REAL4_FORMAT"\n"
      "ppnParams->phi\t= %"LAL_REAL4_FORMAT"\n"
      "ppnParams->psi\t= %"LAL_REAL4_FORMAT"\n"
      "ppnParams->fStartIn\t= %"LAL_REAL4_FORMAT"\n"
      "ppnParams->fStopIn\t= %"LAL_REAL4_FORMAT"\n"
      "ppnParams->position.longitude\t= %"LAL_REAL8_FORMAT"\n"
      "ppnParams->position.latitude\t= %"LAL_REAL8_FORMAT"\n"
      "ppnParams->position.system\t= %d\n"
      "ppnParams->epoch.gpsSeconds\t= %"LAL_INT4_FORMAT"\n"
      "ppnParams->epoch.gpsNanoSeconds\t= %"LAL_INT4_FORMAT"\n"
      "ppnParams->tC\t= %"LAL_REAL8_FORMAT"\n"
      "ppnParams->dfdt\t =%"LAL_REAL4_FORMAT"\n",
      ppnParams->mTot,
      ppnParams->eta,
      ppnParams->d,
      ppnParams->inc,
      ppnParams->phi,
      ppnParams->psi,
      ppnParams->fStartIn,
      ppnParams->fStopIn,
      ppnParams->position.longitude,
      ppnParams->position.latitude,
      ppnParams->position.system,
      ppnParams->epoch.gpsSeconds,
      ppnParams->epoch.gpsNanoSeconds,
      ppnParams->tc,
      ppnParams->dfdt );
  LALInfo( status, warnMsg );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/** \see See \ref GenerateInspiral_h for documentation */
int
XLALGetOrderFromString(
    CHAR       * restrict thisEvent,
    LALPNOrder * restrict order
    )
{

#ifndef LAL_NDEBUG
  if ( !thisEvent )
    XLAL_ERROR( XLAL_EFAULT );

  if ( !order )
    XLAL_ERROR( XLAL_EFAULT );
#endif

  if ( strstr(thisEvent, "newtonian") )
  {
    *order = LAL_PNORDER_NEWTONIAN;
  }
  else if ( strstr(thisEvent, "oneHalfPN") )
  {
    *order = LAL_PNORDER_HALF;
  }
  else if ( strstr(thisEvent, "onePN") )
  {
    *order = LAL_PNORDER_ONE;
  }
  else if ( strstr(thisEvent, "onePointFivePN") )
  {
    *order = LAL_PNORDER_ONE_POINT_FIVE;
  }
  else if ( strstr(thisEvent, "twoPN") )
  {
    *order = LAL_PNORDER_TWO;
  }
  else if ( strstr(thisEvent, "twoPointFivePN") )
  {
    *order = LAL_PNORDER_TWO_POINT_FIVE;
  }
  else if (strstr(thisEvent, "threePN") )
  {
    *order = LAL_PNORDER_THREE;
  }
  else if ( strstr(thisEvent, "threePointFivePN") )
  {
    *order = LAL_PNORDER_THREE_POINT_FIVE;
  }
  else if ( strstr(thisEvent, "pseudoFourPN") )
  {
    *order = LAL_PNORDER_PSEUDO_FOUR;
  }
  else
  {
    XLALPrintError( "Cannot parse order from string: %s\n", thisEvent );
    XLAL_ERROR( XLAL_EINVAL );
  }

  return XLAL_SUCCESS;
}


/**	Convert a string provided by the #CoherentGW structure in order to retrieve
 *	the approximant of the waveform to generate.
 *	@param[out]	inter	: the level of the spin interaction
 *	@param[in]	thisEvent	: string containing the spin interaction
 *	@return error code
 */
int XLALGetInteractionFromString(LALSimInspiralInteraction *inter, CHAR *thisEvent) {
  if (strstr(thisEvent, "NO")) {
    *inter = LAL_SIM_INSPIRAL_INTERACTION_NONE;
  } else if (strstr(thisEvent, "SO15")) {
    *inter = LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_15PN;
  } else if (strstr(thisEvent,"SS")) {
    *inter = LAL_SIM_INSPIRAL_INTERACTION_SPIN_SPIN_2PN;
  } else if (strstr(thisEvent,"SELF")) {
    *inter = LAL_SIM_INSPIRAL_INTERACTION_SPIN_SPIN_SELF_2PN;
  } else if (strstr(thisEvent, "QM")) {
    *inter = LAL_SIM_INSPIRAL_INTERACTION_QUAD_MONO_2PN;
  } else if (strstr(thisEvent, "SO25")) {
    *inter = LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_25PN;
  } else if (strstr(thisEvent, "SO")) {
    *inter = LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_3PN;
  } else if (strstr(thisEvent, "ALL_SPIN")) {
    *inter = LAL_SIM_INSPIRAL_INTERACTION_ALL_SPIN;
  } else if (strstr(thisEvent, "TIDAL5PN")) {
    *inter = LAL_SIM_INSPIRAL_INTERACTION_TIDAL_5PN;
  } else if (strstr(thisEvent, "TIDAL")) {
    *inter = LAL_SIM_INSPIRAL_INTERACTION_TIDAL_6PN;
  } else if (strstr(thisEvent, "ALL")){
    *inter = LAL_SIM_INSPIRAL_INTERACTION_ALL;
  } else {
    XLALPrintError( "Cannot parse LALSimInspiralInteraction from string: %s\n Please add 'ALL' to the above string for including all spin interactions\n", thisEvent );
    XLAL_ERROR( XLAL_EINVAL );
  }

  return XLAL_SUCCESS;
}

/** \see See \ref GenerateInspiral_h for documentation */
int XLALGetAxisChoiceFromString(InputAxis *axisChoice, CHAR *thisEvent) {
  if (strstr(thisEvent, "TotalJ")) {
    *axisChoice = TotalJ;
  } else if  (strstr(thisEvent, "OrbitalL")) {
    *axisChoice = OrbitalL;
  }
  else  
    *axisChoice = View;
  return XLAL_SUCCESS;
}

/** \see See \ref GenerateInspiral_h for documentation */
int XLALGetAdaptiveIntFromString(UINT4 *fixedStep, CHAR *thisEvent) {
  if (strstr(thisEvent, "fixedStep")) {
    *fixedStep = 1;
  } else 
    *fixedStep = 0;
  return XLAL_SUCCESS;
}

/** \see See \ref GenerateInspiral_h for documentation */
int XLALGetInspiralOnlyFromString(UINT4 *inspiralOnly, CHAR *thisEvent) {
  if (strstr(thisEvent, "inspiralOnly")) {
    *inspiralOnly = 1;
  }
  else
    *inspiralOnly = 0;
  return XLAL_SUCCESS;
}

/** \see See \ref GenerateInspiral_h for documentation */
int
XLALGetApproximantFromString(
    CHAR        * restrict thisEvent,
    Approximant * restrict approximant
    )
{
  /* Function to search for the approximant into a string */


#ifndef LAL_NDEBUG
  if ( !thisEvent )
    XLAL_ERROR( XLAL_EFAULT );

  if ( !approximant )
    XLAL_ERROR( XLAL_EFAULT );
#endif

  if ( strstr(thisEvent, "TaylorT1" ) )
  {
    *approximant = TaylorT1;
  }
  else if ( strstr(thisEvent, "TaylorT2" ) )
  {
    *approximant = TaylorT2;
  }
  else if ( strstr(thisEvent, "TaylorF2" ) )
  {
    *approximant = TaylorF2;
  }
  else if ( strstr(thisEvent, "TaylorT3" ) )
  {
    *approximant = TaylorT3;
  }
  else if ( strstr(thisEvent, "EOBNRv2HM" ) )
  {
    *approximant = EOBNRv2HM;
  }
  else if ( strstr(thisEvent, "EOBNRv2" ) )
  {
    *approximant = EOBNRv2;
  }
  else if ( strstr(thisEvent, "EOBNR" ) )
  {
    *approximant = EOBNR;
  }
  else if ( strstr(thisEvent, "EOB" ) )
  {
    *approximant = EOB;
  }
  else if ( strstr(thisEvent, "PhenSpinTaylorRD" ) )
  {
    *approximant = PhenSpinTaylorRD;
  }
  else if ( strstr(thisEvent, "SpinTaylorT4" ) )
  {
    *approximant = SpinTaylorT4;
  }
  else if ( strstr(thisEvent, "SpinTaylorFrameless" ) )
  {
	  *approximant = SpinTaylorFrameless;
  }
  else if ( strstr(thisEvent, "SpinTaylorT3" ) )
  {
    *approximant = SpinTaylorT3;
  }
  else if ( strstr(thisEvent, "SpinTaylor" ) )
  {
    *approximant = SpinTaylor;
  }
  else if ( strstr(thisEvent, "SpinQuadTaylor" ) )
  {
	*approximant = SpinQuadTaylor;
  }
  else if ( strstr(thisEvent, "PadeT1" ) )
  {
    *approximant = PadeT1;
  }
  else if ( strstr(thisEvent, "AmpCorPPN" ) )
  {
    *approximant = AmpCorPPN;
  }
  else if ( strstr(thisEvent, "GeneratePPN" ) )
  {
    *approximant = GeneratePPN;
  }
  else if ( strstr(thisEvent, "TaylorT4" ) )
  {
    *approximant = TaylorT4;
  }
  else if ( strstr(thisEvent, "NumRel" ) )
  {
    *approximant = NumRel;
  }
  else if ( strstr(thisEvent, "Ninja2" ) )
  {
    *approximant = NumRelNinja2;
  }
  else if ( strstr(thisEvent, "IMRPhenomA" ) )
  {
    *approximant = IMRPhenomA;
  }
  else if ( strstr(thisEvent, "IMRPhenomB" ) )
  {
    *approximant = IMRPhenomB;
  }
  else
  {
    XLALPrintError( "Cannot parse approximant from string: %s \n", thisEvent );
    XLAL_ERROR( XLAL_EINVAL );
  }

  return XLAL_SUCCESS;
}

/** \see See \ref GenerateInspiral_h for documentation */
int
XLALGetTaperFromString(
    LALSimInspiralApplyTaper * restrict taper,
    CHAR                     * restrict thisEvent
    )
{

  if ( ! strcmp( "TAPER_NONE", thisEvent ) )
  {
    *taper = LAL_SIM_INSPIRAL_TAPER_NONE;
  }
  else if ( ! strcmp( "TAPER_START", thisEvent ) )
  {
    *taper = LAL_SIM_INSPIRAL_TAPER_START;
  }
  else if ( ! strcmp( "TAPER_END", thisEvent ) )
  {
    *taper = LAL_SIM_INSPIRAL_TAPER_END;
  }
  else if ( ! strcmp( "TAPER_STARTEND", thisEvent ) )
  {
    *taper = LAL_SIM_INSPIRAL_TAPER_STARTEND;
  }
  else
  {
    XLALPrintError( "Invalid injection tapering option specified: %s\n", thisEvent );
    XLAL_ERROR( XLAL_EINVAL );
  }

  return XLAL_SUCCESS;
}

/** \see See \ref GenerateInspiral_h for documentation */
int
XLALGenerateInspiralPopulatePPN(
    PPNParamStruc    * restrict ppnParams,
    SimInspiralTable * restrict thisEvent
    )
{
#ifndef LAL_NDEBUG
  if ( !ppnParams || !thisEvent )
    XLAL_ERROR( XLAL_EFAULT );
#endif

  /* input fields */
  ppnParams->mTot     = thisEvent->mass1 + thisEvent->mass2;
  ppnParams->eta      = thisEvent->eta;
  ppnParams->d        = thisEvent->distance* 1.0e6 * LAL_PC_SI; /*in Mpc*/
  ppnParams->inc      = thisEvent->inclination;
  ppnParams->phi      = thisEvent->coa_phase;
  ppnParams->ampOrder = thisEvent->amp_order;

  /* frequency cutoffs */
  if ( thisEvent->f_lower > 0 )
  {
    ppnParams->fStartIn = thisEvent->f_lower;
  }
  else
  {
    XLALPrintError( 
        "f_lower must be specified in the injection file generation.\n" );
    XLAL_ERROR( XLAL_EINVAL );
  }
  ppnParams->fStopIn  = -1.0 /
    ( sqrt(216.0) * LAL_PI * ppnParams->mTot * LAL_MTSUN_SI);

  /* passed fields */
  ppnParams->position.longitude   = thisEvent->longitude;
  ppnParams->position.latitude    = thisEvent->latitude;
  ppnParams->position.system      = COORDINATESYSTEM_EQUATORIAL;
  ppnParams->psi                  = thisEvent->polarization;
  ppnParams->epoch.gpsSeconds     = 0;
  ppnParams->epoch.gpsNanoSeconds = 0;

  return XLAL_SUCCESS;
}


/** \see See \ref GenerateInspiral_h for documentation */
int
XLALGenerateInspiralPopulateInspiral(
    InspiralTemplate * restrict inspiralParams,
    SimInspiralTable * restrict thisEvent,
    PPNParamStruc    * restrict ppnParams
    )
{

#ifndef LAL_NDEBUG
  if ( !inspiralParams || !thisEvent || !ppnParams )
    XLAL_ERROR( XLAL_EFAULT );
#endif

  /* --- Let's fill the inspiral structure now --- */
  inspiralParams->mass1	  =  thisEvent->mass1;  	/* masses 1 */
  inspiralParams->mass2	  =  thisEvent->mass2;  	/* masses 2 */
  inspiralParams->fLower  =  ppnParams->fStartIn; /* lower cutoff frequency */
  inspiralParams->fFinal  =  thisEvent->f_final;
  inspiralParams->fCutoff = 1./ (ppnParams->deltaT)/2.-1;

  /* -1 to be  in agreement with the inspiral assert. */
  inspiralParams->tSampling	  = 1./ (ppnParams->deltaT); /* sampling*/
  inspiralParams->signalAmplitude = 1.;
  inspiralParams->distance	  =  thisEvent->distance * LAL_PC_SI * 1e6;

  /* distance in Mpc */
  inspiralParams->startTime	  =  0.0;
  inspiralParams->startPhase	  =  thisEvent->coa_phase;

  inspiralParams->OmegaS = GENERATEINSPIRAL_OMEGAS;/* EOB 3PN contribution */
  inspiralParams->Theta	 = GENERATEINSPIRAL_THETA; /* EOB 3PN contribution */
  inspiralParams->Zeta2	 = GENERATEINSPIRAL_ZETA2; /* EOB 3PN contribution */

  inspiralParams->alpha	 = -1.;      /* bcv useless for the time being */
  inspiralParams->psi0	 = -1.;      /* bcv useless for the time being */
  inspiralParams->psi3	 = -1.;      /* bcv useless for the time being */
  inspiralParams->alpha1 = -1.;      /* bcv useless for the time being */
  inspiralParams->alpha2 = -1.;      /* bcv useless for the time being */
  inspiralParams->beta   = -1.;      /* bcv useless for the time being */

  /* inclination of the binary */
  /* inclination cannot be equal to zero for SpinTaylor injections */
  if ( inspiralParams->approximant == SpinTaylor && thisEvent->inclination == 0 )
  {
    XLALPrintError( "Inclination cannot be exactly zero for SpinTaylor approximant.\n");
    XLAL_ERROR( XLAL_EINVAL );
  }
  inspiralParams->inclination =  thisEvent->inclination;

  inspiralParams->ieta	    =  1;
  inspiralParams->nStartPad =  0;
  /* increased end padding from zero so that longer waveforms do not
  have errors due to underestimation of number of bins requred */
  inspiralParams->nEndPad   =  16384;

  inspiralParams->massChoice  = m1Andm2;

  /* spin parameters */
  inspiralParams->sourceTheta = GENERATEINSPIRAL_SOURCETHETA;
  inspiralParams->sourcePhi   = GENERATEINSPIRAL_SOURCEPHI;
  inspiralParams->spin1[0]    = thisEvent->spin1x;
  inspiralParams->spin2[0]    = thisEvent->spin2x;
  inspiralParams->spin1[1]    = thisEvent->spin1y;
  inspiralParams->spin2[1]    = thisEvent->spin2y;
  inspiralParams->spin1[2]    = thisEvent->spin1z;
  inspiralParams->spin2[2]    = thisEvent->spin2z;

  inspiralParams->orbitTheta0 = thisEvent->theta0;
  inspiralParams->orbitPhi0   = thisEvent->phi0;
  inspiralParams->qmParameter[0] = thisEvent->qmParameter1;
  inspiralParams->qmParameter[1] = thisEvent->qmParameter2;

  return XLAL_SUCCESS;
}

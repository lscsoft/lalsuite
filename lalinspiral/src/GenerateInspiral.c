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
  approximant = XLALGetApproximantFromString(thisEvent->waveform);
  if ( (int) approximant == XLAL_FAILURE)
    ABORTXLAL(status);

  order = XLALGetOrderFromString(thisEvent->waveform);
  if ( (int) order == XLAL_FAILURE)
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
      inspiralParams.interaction = XLALGetInteractionFromString(
          thisEvent->waveform);
      if ( (int) inspiralParams.interaction == XLAL_FAILURE)
        ABORTXLAL(status);
    }

    if (approximant == PhenSpinTaylorRD) {
      xlalErrno = 0;
      /* These next three functions cannot fail - they return a default value
         if no target string is present - so we don't check for failure. */
      inspiralParams.axisChoice = XLALGetFrameAxisFromString(
          thisEvent->waveform);
      inspiralParams.fixedStep = XLALGetAdaptiveIntFromString(
          thisEvent->waveform);
      inspiralParams.inspiralOnly = XLALGetInspiralOnlyFromString(
          thisEvent->waveform);
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

  /* If no waveform has been generated.*/
  if ( waveform->a == NULL && waveform->h == NULL)
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
  inspiralParams->ampOrder = ppnParams->ampOrder;

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

/*
*  Copyright (C) 2007 Drew Keppel, Duncan Brown, Gareth Jones, Peter Shawhan, Thomas Cokelaer, Laszlo Vereb
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

#ifndef _GENERATEINSPIRAL_H
#define _GENERATEINSPIRAL_H

#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/GeneratePPNInspiral.h>

#include <lal/LIGOMetadataTables.h>
#include <lal/SeqFactories.h>

#include <lal/Units.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * \addtogroup GenerateInspiral_h
 * \author Cokelaer, T.
 *
 * \brief %Header file for the inspiral injection interface code.
 *
 * The code contained in GenerateInspiral.c is an interface between the
 * injection package and the inspiral package. More precisely, the
 * function GenerateInspiral.c is used within the FindChirpSimulation.c
 * file of the FindChirp package in order to inject waveforms into real
 * data. The injection is done through the inject package in order to
 * take into account the interferometer position, binary orientation ...
 *
 * GenerateInspiral has the capability of injecting both waveform designed
 * within the inspiral package (TaylorT1, T2, T3, PadeT1, EOB, and spinning
 * waveform) and the inject package (so-called PPN waveform).
 * also a test code as well which allows to check the output of
 * code. It is called InjectionInterfaceTest.c
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/GenerateInspiral.h>
 * \endcode
 *
 * <dl>
 * <dt><tt>LALGenerateInspiral()</tt></dt><dd> create an inspiral binary
 * waveform generated either by the \ref pkg_inspiral (#EOB,
 * #EOBNR, #PadeT1, #TaylorT1, #TaylorT2, #TaylorT3, #SpinTaylor, #PhenSpinTaylorRD, #SpinQuadTaylor)
 * or the \c inject package (#GeneratePPN).  It is used in the module
 * \c FindChirpSimulation in \ref pkg_findchirp.
 *
 * There are three  parsed arguments
 * <ul>
 * <li> a ::CoherentGW  structure which stores amplitude,
 * frequency and phase of the  waveform (output)</li>
 * <li> a \c thisEvent  structure which provides some
 * waveform parameters (input)</li>
 * <li> a \c PPNParamStruc which gives some input
 * parameters needed by the GeneratePPN waveform  generation. That
 * arguments is also used as an output by all the different
 * approximant  (output/input).</li>
 * </ul>
 *
 * The input must be composed of a valid thisEvent structure as well as
 * the  variable deltaT of the PPNParamStruc. All others variables
 * of the PPNParamStruc are populated within that function.</dd>
 *
 * <dt><tt>XLALGenerateInspiralPopulatePPN()</tt></dt><dd> Populate the
 * PPNParamsStruc with the input argument \c thisEvent. That
 * structure is used by both inspiral waveforms inject waveforms.</dd>
 *
 * <dt><tt>XLALGenerateInspiralPopulateInspiral()</tt></dt><dd>  Populate the
 * InspiralTemplate structure if the model chosen belongs to the
 * inspiral package.
 * </dd>
 * </dl>
 *
 * \heading{Notes}
 * Inject only time-domain waveforms for the time being such as GeneratePPN,
 * TaylorT1, TaylorT2, TaylorT3, PadeT1 and EOB , SpinTaylor, PhenSpinTaylorRD.
 *
 */
/*@{*/

/**\name Error Codes */
/*@{*/
#define GENERATEINSPIRALH_ENORM 0	/**< Normal exit */
#define GENERATEINSPIRALH_ENULL 1	/**< Null pointer */
#define GENERATEINSPIRALH_EDFDT 2	/**< Waveform sampling interval is too large */
#define GENERATEINSPIRALH_EZERO 3	/**< inclination zero for SpinTaylor waveform */
/*@}*/

/** \cond DONT_DOXYGEN */
#define GENERATEINSPIRALH_MSGENORM "Normal exit"
#define GENERATEINSPIRALH_MSGENULL "Null pointer"
#define GENERATEINSPIRALH_MSGEDFDT "Waveform sampling interval is too large"
#define GENERATEINSPIRALH_MSGEZERO "inclination zero for SpinTaylor waveform"
/** \endcond */


/**
 * \name Parameter for the EOB at 3PN.
 * In principle, the three following parameter should be set to zero.
 */
/*@{*/
#define GENERATEINSPIRAL_ZETA2       0.
#define GENERATEINSPIRAL_OMEGAS      0.
#define GENERATEINSPIRAL_THETA       0.
/*@}*/

/** \name For the spinning case, might be changed later or include in the injection itself */
/*@{*/
#define GENERATEINSPIRAL_SOURCETHETA 1.
#define GENERATEINSPIRAL_SOURCEPHI   2.
/*@}*/

/** Default low freqnecy cutoff for injections */
#define GENERATEINSPIRAL_DEFAULT_FLOWER 40


void
LALGenerateInspiral(
    LALStatus        *status,
    CoherentGW       *waveform,
    SimInspiralTable *params,
    PPNParamStruc    *ppnParamsInputOutput
    );

int
XLALGenerateInspiralPopulatePPN(
    PPNParamStruc    * restrict ppnParams,
    SimInspiralTable * restrict thisEvent
    );

int
XLALGenerateInspiralPopulateInspiral(
    InspiralTemplate * restrict inspiralParams,
    SimInspiralTable * restrict thisEvent,
    PPNParamStruc    * restrict ppnParams
    );

/*@}*/

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _GENERATEINSPIRAL_H */

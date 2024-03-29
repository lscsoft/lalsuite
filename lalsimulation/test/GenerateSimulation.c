/*
*  Copyright (C) 2011 Nickolas Fotopoulos, Evan Ochsner, 2016 Riccardo Sturani
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

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <string.h>

#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>
#include <lal/XLALError.h>
#include <lal/LALAdaptiveRungeKuttaIntegrator.h>
#include <lal/LALSimInspiralWaveformParams.h>

/* internal storage is in SI units! */
typedef struct GSParams {
    Approximant approximant;  /**< waveform family or "approximant" */
    LALSimulationDomain domain; /**< flag for time or frequency domain waveform */
    REAL8 phiRef;             /**< phase at fRef */
    REAL8 fRef;               /**< reference frequency */
    REAL8 deltaT;             /**< sampling interval */
    REAL8 deltaF;             /**< frequency resolution */
    REAL8 m1;                 /**< mass of companion 1 */
    REAL8 m2;                 /**< mass of companion 2 */
    REAL8 f_min;              /**< start frequency */
    REAL8 f_max;              /**< end frequency */
    REAL8 distance;           /**< distance of source */
    REAL8 inclination;        /**< inclination of L relative to line of sight */
    REAL8 s1x;                /**< (x,y,z) components of spin of m1 body */
    REAL8 s1y;                /**< z-axis along line of sight, L in x-z plane */
    REAL8 s1z;                /**< dimensionless spin, Kerr bound: |s1| <= 1 */
    REAL8 s2x;                /**< (x,y,z) component ofs spin of m2 body */
    REAL8 s2y;                /**< z-axis along line of sight, L in x-z plane */
    REAL8 s2z;                /**< dimensionless spin, Kerr bound: |s2| <= 1 */
    REAL8 longAscNodes;       /**< longitude of ascending nodes 0<= omega < 2 pi */
    REAL8 ecc;                /**< eccentricity 0<= ecc < 1 */
    REAL8 meanPerAno;         /**< mean periastron anomaly 0<= psi < 2 Pi */
    char outname[256];        /**< file to which output should be written */
    LALDict *params;          /**<Container for all accessory parameters */
    int ampPhase;
    int verbose;
    // New interface specific
    int newinterface;         /**< 0 or 1: Old or New interface respectively */
    int condition;            /**< Conditioning waveform  for (Inverse)Fourier Transform */
    char module[256];         /**< Name of the External python module */
    char object[256];         /**< Name of the python class defined in python-module which includes the methods generate_fd/td_waveform **/
    int m1_newinterface;      /**< Check if mass1 is an input argument */
    int m2_newinterface;      /**< Check if mass2 is an input argument */
} GSParams;


const char * usage =
"Generate a simulation using the lalsimulation library\n\n"
"The following options can be given (will assume a default value if omitted):\n"
"--domain DOM               'TD' for time domain (default) or 'FD' for frequency\n"
"                           domain; not all approximants support all domains\n"
"--amp-phase                If given, will output:\n"
"                           |h+ - i hx|, Arg(h+ - i hx) (TD) or\n"
"                           |h+(f)|, Arg(h+(f)), |hx(f)|, Arg(hx(f)) (FD)\n"
"                           If not given, will output h+ and hx (TD and FD)\n"
"--approximant APPROX       Supported TD approximants:\n"
"                             TaylorT1 (default)\n"
"                             TaylorT2\n"
"                             TaylorT3\n"
"                             TaylorT4\n"
"                             TaylorEt\n"
"                             EccentricTD\n"
"                             IMRPhenomA\n"
"                             IMRPhenomB\n"
"                             IMRPhenomC\n"
"                             IMRPhenomD\n"
"                             IMRPhenomHM\n"
"                             IMRPhenomPv2\n"
"                             IMRPhenomPv3\n"
"                             IMRPhenomPv3HM\n"
"                             EOBNRv2\n"
"                             EOBNRv2HM\n"
"                             SEOBNRv1\n"
"                             SEOBNRv2\n"
"                             SEOBNRv3\n"
"                             SEOBNRv4\n"
"                             SEOBNRv4HM\n"
"                             SEOBNRv4P\n"
"                             SEOBNRv2T\n"
"                             SEOBNRv4T\n"
"                             TEOBResumS\n"
"                             SpinTaylorT4\n"
"                             SpinTaylorT5\n"
"                             PhenSpinTaylor\n"
"                             PhenSpinTaylorRD\n"
"                             SpinDominatedWf\n"
"                             HGimri\n"
"                             NRHybSur3dq8\n"
"                             NR_hdf5\n"
"                           Supported FD approximants:\n"
"                             IMRPhenomA\n"
"                             IMRPhenomB\n"
"                             IMRPhenomC\n"
"                             IMRPhenomD\n"
"                             IMRPhenomP\n"
"                             IMRPhenomPv2\n"
"                             IMRPhenomPv3\n"
"                             IMRPhenomPv3HM\n"
"                             EOBNRv2_ROM\n"
"                             EOBNRv2HM_ROM\n"
"                             SEOBNRv1_ROM_EffectiveSpin\n"
"                             SEOBNRv1_ROM_DoubleSpin\n"
"                             SEOBNRv2_ROM_EffectiveSpin\n"
"                             SEOBNRv2_ROM_DoubleSpin\n"
"                             SEOBNRv4_ROM\n"
"                             SEOBNRv4HM_ROM\n"
"                             SEOBNRv5_ROM\n"
"                             SEOBNRv5HM_ROM\n"
"                             TaylorF2\n"
"                             TaylorF2Ecc\n"
"                             SpinTaylorF2\n"
"                             TaylorR2F4\n"
"                             SpinTaylorT4Fourier\n"
"                             SpinTaylorT5Fourier\n"
"                             TaylorF2RedSpin\n"
"                             TaylorF2RedSpinTidal\n"
"                           NOTE: Other approximants may be available if the\n"
"                           developer forgot to edit this help message\n"
"--phase-order ORD          Twice PN order of phase (default ORD=7 <==> 3.5PN)\n"
"--amp-order ORD            Twice PN order of amplitude (default 0 <==> Newt.)\n"
"--phiRef PHIREF            Phase at the reference frequency (default 0)\n"
"--fRef FREF                Reference frequency in Hz\n"
"                           (default: 0)\n"
"--sample-rate SRATE        Sampling rate of TD approximants in Hz (default 4096)\n"
"--deltaF DF                Frequency bin size for FD approximants in Hz (default 1/8)\n"
"--m1 M1                    Mass of the 1st object in solar masses (default 10)\n"
"--m2 M2                    Mass of the 2nd object in solar masses (default 1.4)\n"
"--inclination IOTA         Angle in radians between line of sight (N) and \n"
"                           orbital angular momentum (L) at the reference\n"
"                           (default 0, face on)\n"
"--spin1x S1X               Vector components for spin of mass1 (default all 0)\n"
"--spin1y S1Y               z-axis=line of sight, L in x-z plane at reference\n"
"--spin1z S1Z               Kerr limit: s1x^2 + s1y^2 + s1z^2 <= 1\n"
"--spin2x S2X               Vector components for spin of mass2 (default all 0)\n"
"--spin2y S2Y               z-axis=line of sight, L in x-z plane at reference\n"
"--spin2z S2Z               Kerr limit: s2x^2 + s2y^2 + s2z^2 <= 1\n"
"--tidal-lambda1 L1         (tidal deformability of mass 1) / (mass of body 1)^5\n"
"                           (~128-2560 for NS, 0 for BH) (default 0)\n"
"--tidal-lambda2 L2         (tidal deformability of mass 2) / (mass of body 2)^5\n"
"                           (~128-2560 for NS, 0 for BH) (default 0)\n"
"--ecc_order eccO           PN order for eccentric correction term\n"
"                           (-1 for maximum order) (default -1)\n"
"--f_ecc f                  reference frequency for initial eccentricity\n"
"                           (10Hz) (default 10)\n"
"--tidal-lambda-octu1 L31         (octupolar tidal deformability of mass 1) / (mass of body 1)^7\n"
"                           (0 for BH) (default 0)\n"
"--tidal-lambda-octu2 L32         (octupolar tidal deformability of mass 2) / (mass of body 2)^7\n"
"                           (0 for BH) (default 0)\n"
"--tidal-quadfmode1 W21      dimensionless quadrupolar f-mode angular frequency of mass 1, normalized to mass 1\n"
"--tidal-quadfmode2 W22      dimensionless quadrupolar f-mode angular frequency of mass 2, normalized to mass 2\n"
"--tidal-octufmode1 W31      dimensionless octupolar f-mode angular frequency of mass 1, normalized to mass 1\n"
"--tidal-octufmode2 W32      dimensionless octupolar f-mode angular frequency of mass 2, normalized to mass 2\n"
"--quad-mon1 kappa1          dimensionless quadrupole-monopole parameter of mass 1, value 1 for a black hole\n"
"--quad-mon2 kappa2          dimensionless quadrupole-monopole parameter of mass 2, value 1 for a black hole\n"
"--spin-order ORD           Twice PN order of spin effects\n"
"                           (default ORD=-1 <==> All spin effects)\n"
"--tidal-order ORD          Twice PN order of tidal effects\n"
"                           (default ORD=-1 <==> All tidal effects)\n"
"--f-min FMIN               Lower frequency to start waveform in Hz (default 40)\n"
"--f-max FMAX               Frequency at which to stop waveform in Hz\n"
"                           (default: generate as much as possible)\n"
"--distance D               Distance in Mpc (default 100)\n"
"--long-asc-nodes omega     Longitude of ascending nodes in radians (default 0)\n"
"--eccentricity ecc         Eccentricity (default 0)\n"
"--mean-per-ano psi         Mean periastron anomaly in radians (default 0)\n"
"--axis AXIS                for PhenSpin: 'View' (default), 'TotalJ', 'OrbitalL'\n"
"--nr-file                  for NR_hdf5\n"
"--nonGRpar NAME VALUE      add the nonGRparam with name 'NAME' and value 'VALUE'\n"
"                           Supported names:\n"
"                             NonGRPhi1\n"
"                             NonGRPhi2\n"
"                             NonGRPhi3\n"
"                             NonGRPhi4\n"
"                             NonGRDChi0\n"
"                             NonGRDChi1\n"
"                             NonGRDChi2\n"
"                             NonGRDChi3\n"
"                             NonGRDChi4\n"
"                             NonGRDChi5\n"
"                             NonGRDChi5L\n"
"                             NonGRDChi6\n"
"                             NonGRDChi6L\n"
"                             NonGRDChi7\n"
"--higher-modes VALUE       specify l modes with value 'VALUE' (L2 or RESTRICTED is default)\n"
"--outname FNAME            Output to file FNAME (default 'simulation.dat')\n"
"--verbose                  If included, add verbose output\n"
"--phenomXHMMband float     Threshold parameter for the Multibanding of IMRPhenomXHM. By default set to 10^-3. If set to 0 then do not use multibanding.\n"
"--modesList string         List of modes to be used by the model, e.g. --modesList '2,2, 2,-2, 2,1, 2,-1'. \n"
"                           To use all the modes available in the mode do not add --modesList. Do not use if you do not know which modes the model returns!\n"
"--phenomXPHMMband float    Threshold parameter for the Multibanding of the Euler angles in IMRPhenomXPHM. \n"
"                           Options and default values can be found in https://lscsoft.docs.ligo.org/lalsuite/lalsimulation/group___l_a_l_sim_i_m_r_phenom_x__c.html.\n"
"--phenomXPHMUseModes int   Switch between the two polarization functions IMRPhenomXPHM and IMRPhenomXPHMFromModes.\n"
"                           Options and default values can be found in https://lscsoft.docs.ligo.org/lalsuite/lalsimulation/group___l_a_l_sim_i_m_r_phenom_x__c.html.\n"
"--phenomXPrecVersion int   Choose precessing version for the Euler angles.\n"
"                           Options and default values can be found in https://lscsoft.docs.ligo.org/lalsuite/lalsimulation/group___l_a_l_sim_i_m_r_phenom_x__c.html.\n"
"--phenomXPFinalSpinMod int Choose final spin prescription.\n"
"                           Options and default values can be found in https://lscsoft.docs.ligo.org/lalsuite/lalsimulation/group___l_a_l_sim_i_m_r_phenom_x__c.html.\n"
"--newinterface             If included, it will use XLALSimInspiralGenerateT(F)DWaveform instead of XLALSimInspiralChooseT(F)DWaveform\n"
"                           Can use flexible parameters: --total_mass --chirp_mass --mass_ratio ... --spin1_norm --spin1_tilt --spin1_phi ...\n"
"--python-module char       Name of the External python module\n"
"--python-object char       Name of the python class defined in python-module which includes the methods generate_fd/td_waveform\n"
;

/* Parse command line, sanity check arguments, and return a newly
 * allocated GSParams object */
static GSParams *parse_args(ssize_t argc, char **argv) {
    ssize_t i;
    GSParams *params;
    params = (GSParams *) XLALMalloc(sizeof(GSParams));
    memset(params, 0, sizeof(GSParams));
    params->params=XLALCreateDict();
    /* Set default values to the arguments */
    params->approximant = TaylorT1;
    params->domain = LAL_SIM_DOMAIN_TIME;
    XLALSimInspiralWaveformParamsInsertPNPhaseOrder(params->params, -1);
    XLALSimInspiralWaveformParamsInsertPNAmplitudeOrder(params->params, -1);
    params->phiRef = 0.;
    params->deltaT = 1./4096.;
    params->deltaF = 0.125;
    params->m1 = 10. * LAL_MSUN_SI;
    params->m2 = 1.4 * LAL_MSUN_SI;
    params->f_min = 40.;
    params->fRef = 0.;
    params->f_max = 0.; /* Generate as much as possible */
    params->distance = 100. * 1e6 * LAL_PC_SI;
    params->inclination = 0.;
    params->s1x = 0.;
    params->s1y = 0.;
    params->s1z = 0.;
    params->s2x = 0.;
    params->s2y = 0.;
    params->s2z = 0.;
    params->longAscNodes = 0.;
    params->ecc = 0.;
    params->meanPerAno = 0.;
    params->m1_newinterface = 0;
    params->m2_newinterface = 0;
    strcpy(params->module, "");
    strcpy(params->object, "");
    XLALSimInspiralWaveformParamsInsertTidalLambda1(params->params, 0.);
    XLALSimInspiralWaveformParamsInsertTidalLambda2(params->params, 0.);
    /* Note: given a LALDict, SEOBNRv2T/v4T will check wether other parameters ( TidalOctupolarLambda, TidalQuadrupolarFMode, TidalOctupolarFMode, dQuadMon) are present in the LALDict */
    /* If present, this value is looked up and used, if absent, the parameter is computed from quadrupolar lambda through universal relations */
    /* We therefore do not Insert default values for these parameters, they will be inserted in LALDict only if specified on the commandline */
    XLALSimInspiralWaveformParamsInsertPNEccentricityOrder(params->params, -1);
    XLALSimInspiralWaveformParamsInsertEccentricityFreq(params->params, 10.0);
    XLALSimInspiralWaveformParamsInsertTidalOctupolarLambda1(params->params, 0.);
    XLALSimInspiralWaveformParamsInsertTidalOctupolarLambda2(params->params, 0.);
    XLALSimInspiralWaveformParamsInsertTidalQuadrupolarFMode1(params->params, 0.);
    XLALSimInspiralWaveformParamsInsertTidalQuadrupolarFMode2(params->params, 0.);
    XLALSimInspiralWaveformParamsInsertTidalOctupolarFMode1(params->params, 0.);
    XLALSimInspiralWaveformParamsInsertTidalOctupolarFMode2(params->params, 0.);
    snprintf(params->outname, sizeof(params->outname), "simulation.dat"); /* output to this file */
    params->verbose = 0; /* No verbosity */
    /* consume command line */
    for (i = 1; i < argc; ++i) {
        if ((strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "--help") == 0)) {
            printf("%s", usage);
            XLALFree(params);
            exit(0);
        } else if (strcmp(argv[i], "--verbose") == 0) {
            params->verbose = 1;
        } else if (strcmp(argv[i], "--newinterface") == 0) {
            params->newinterface = 1;
        } else if (strcmp(argv[i], "--amp-phase") == 0) {
	    params->ampPhase = 1;
        } else if ( ( i == argc ) || ( !argv[i+1] ) ) {
          XLALPrintError("Error: value required for option %s\n", argv[i]);
        } else if (strcmp(argv[i], "--approximant") == 0) {
            params->approximant = XLALSimInspiralGetApproximantFromString(argv[++i]);
            if ( (int) params->approximant == XLAL_FAILURE) {
                XLALPrintError("Error: invalid value %s for --interaction-flag\n", argv[i]);
                goto fail;
            }
        } else if (strcmp(argv[i], "--domain") == 0) {
            i++;
            if (strcmp(argv[i], "TD") == 0)
                params->domain = LAL_SIM_DOMAIN_TIME;
            else if (strcmp(argv[i], "FD") == 0)
                params->domain = LAL_SIM_DOMAIN_FREQUENCY;
            else {
                XLALPrintError("Error: Unknown domain\n");
                goto fail;
            }
        } else if (strcmp(argv[i], "--phase-order") == 0) {
	    XLALSimInspiralWaveformParamsInsertPNPhaseOrder(params->params, atoi(argv[++i]));
        } else if (strcmp(argv[i], "--amp-order") == 0) {
	    XLALSimInspiralWaveformParamsInsertPNAmplitudeOrder(params->params, atoi(argv[++i]));
        } else if (strcmp(argv[i], "--phiRef") == 0) {
            params->phiRef = atof(argv[++i]);
        } else if (strcmp(argv[i], "--fRef") == 0) {
            params->fRef = atof(argv[++i]);
        } else if (strcmp(argv[i], "--sample-rate") == 0) {
            params->deltaT = 1./atof(argv[++i]);
        } else if (strcmp(argv[i], "--deltaF") == 0) {
            params->deltaF = atof(argv[++i]);
        } else if (strcmp(argv[i], "--m1") == 0) {
            params->m1 = atof(argv[++i]) * LAL_MSUN_SI;
            params->m1_newinterface = 1;
        } else if (strcmp(argv[i], "--m2") == 0) {
            params->m2 = atof(argv[++i]) * LAL_MSUN_SI;
            params->m2_newinterface = 1;
        } else if (strcmp(argv[i], "--spin1x") == 0) {
            params->s1x = atof(argv[++i]);
        } else if (strcmp(argv[i], "--spin1y") == 0) {
            params->s1y = atof(argv[++i]);
        } else if (strcmp(argv[i], "--spin1z") == 0) {
            params->s1z = atof(argv[++i]);
        } else if (strcmp(argv[i], "--spin2x") == 0) {
            params->s2x = atof(argv[++i]);
        } else if (strcmp(argv[i], "--spin2y") == 0) {
            params->s2y = atof(argv[++i]);
        } else if (strcmp(argv[i], "--spin2z") == 0) {
            params->s2z = atof(argv[++i]);
        } else if (strcmp(argv[i], "--tidal-lambda1") == 0) {
	    XLALSimInspiralWaveformParamsInsertTidalLambda1(params->params, atof(argv[++i]));
        } else if (strcmp(argv[i], "--tidal-lambda2") == 0) {
	    XLALSimInspiralWaveformParamsInsertTidalLambda2(params->params, atof(argv[++i]));
        } else if (strcmp(argv[i], "--ecc_order") == 0) {
	    XLALSimInspiralWaveformParamsInsertPNEccentricityOrder(params->params, atoi(argv[++i]));
        } else if (strcmp(argv[i], "--f_ecc") == 0) {
	    XLALSimInspiralWaveformParamsInsertEccentricityFreq(params->params, atof(argv[++i]));
        } else if (strcmp(argv[i], "--tidal-lambda-octu1") == 0) {
        XLALSimInspiralWaveformParamsInsertTidalOctupolarLambda1(params->params, atof(argv[++i]));
        } else if (strcmp(argv[i], "--tidal-lambda-octu2") == 0) {
        XLALSimInspiralWaveformParamsInsertTidalOctupolarLambda2(params->params, atof(argv[++i]));
        } else if (strcmp(argv[i], "--tidal-quadfmode1") == 0) {
        XLALSimInspiralWaveformParamsInsertTidalQuadrupolarFMode1(params->params, atof(argv[++i]));
        } else if (strcmp(argv[i], "--tidal-quadfmode2") == 0) {
        XLALSimInspiralWaveformParamsInsertTidalQuadrupolarFMode2(params->params, atof(argv[++i]));
        } else if (strcmp(argv[i], "--tidal-octufmode1") == 0) {
        XLALSimInspiralWaveformParamsInsertTidalOctupolarFMode1(params->params, atof(argv[++i]));
        } else if (strcmp(argv[i], "--tidal-octufmode2") == 0) {
        XLALSimInspiralWaveformParamsInsertTidalOctupolarFMode2(params->params, atof(argv[++i]));
        } else if (strcmp(argv[i], "--quad-mon1") == 0) {
        XLALSimInspiralWaveformParamsInsertdQuadMon1(params->params, atof(argv[++i]) - 1.); /* dQuadMon is kappa-1 */
        } else if (strcmp(argv[i], "--quad-mon2") == 0) {
        XLALSimInspiralWaveformParamsInsertdQuadMon2(params->params, atof(argv[++i]) - 1.); /* dQuadMon is kappa-1 */
        } else if (strcmp(argv[i], "--spin-order") == 0) {
	    XLALSimInspiralWaveformParamsInsertPNSpinOrder(params->params, atoi(argv[++i]));
        } else if (strcmp(argv[i], "--tidal-order") == 0) {
	    XLALSimInspiralWaveformParamsInsertPNTidalOrder(params->params, atoi(argv[++i]));
        } else if (strcmp(argv[i], "--f-min") == 0) {
            params->f_min = atof(argv[++i]);
        } else if (strcmp(argv[i], "--f-max") == 0) {
            params->f_max = atof(argv[++i]);
        } else if (strcmp(argv[i], "--distance") == 0) {
            params->distance = atof(argv[++i]) * 1e6 * LAL_PC_SI;
        } else if (strcmp(argv[i], "--inclination") == 0) {
            params->inclination = atof(argv[++i]);
        } else if (strcmp(argv[i], "--long-asc-nodes") == 0) {
            params->longAscNodes = atof(argv[++i]);
        } else if (strcmp(argv[i], "--eccentricity") == 0) {
            params->ecc = atof(argv[++i]);
        } else if (strcmp(argv[i], "--mean-per-ano") == 0) {
            params->meanPerAno = atof(argv[++i]);
        } else if (strcmp(argv[i], "--axis") == 0) {
	    XLALSimInspiralWaveformParamsInsertFrameAxis(params->params, XLALGetFrameAxisFromString(argv[++i]) );
            if ( (int) XLALSimInspiralWaveformParamsLookupFrameAxis(params->params)
                    == (int) XLAL_FAILURE) {
                XLALPrintError("Error: invalid value %s for --axis\n", argv[i]);
                goto fail;
            }
        } else if (strcmp(argv[i], "--modes") == 0) {
	    XLALSimInspiralWaveformParamsInsertModesChoice(params->params, XLALSimInspiralGetHigherModesFromString(argv[++i]));
            if ( (int) XLALSimInspiralWaveformParamsLookupModesChoice(params->params) == (int) XLAL_FAILURE) {
                XLALPrintError("Error: invalid value %s for --modes\n", argv[i]);
                goto fail;
            }
        } else if (strcmp(argv[i], "--nonGRpar") == 0) {
	    char name[100];
	    strcpy(name,argv[++i]);
            if ( ( i == argc ) || ( !argv[i+1] ) ) {
              XLALPrintError("Error: 'name value' pair required for option %s\n", argv[i-1]);
            } else if(strcmp(name,"Phi1")==0) {
	      XLALSimInspiralWaveformParamsInsertNonGRPhi1(params->params,atof(argv[++i]));
	    }
        } else if (strcmp(argv[i], "--outname") == 0) {
            snprintf(params->outname, sizeof(params->outname), "%s", argv[++i]);
        } else if (strcmp(argv[i], "--nr-file") == 0) {
	  XLALSimInspiralWaveformParamsInsertNumRelData(params->params, argv[++i]);
        } else if (strcmp(argv[i], "--modesList") == 0) {
            char *current_mode;
            int modes[50], iii = 0;
            while ((current_mode = strtok(iii ? NULL : argv[++i], ",")) != NULL){
              modes[iii++] = atoi(current_mode);
            }
            LALValue *ModeArray = XLALSimInspiralCreateModeArray();
            for (int ii=0; ii<iii; ii+=2){
               XLALSimInspiralModeArrayActivateMode(ModeArray, modes[ii], modes[ii+1]);
               XLALSimInspiralWaveformParamsInsertModeArray(params->params, ModeArray);
             }
             XLALDestroyValue(ModeArray);

        }else if(strcmp(argv[i], "--phenomXHMMband") == 0){
            XLALSimInspiralWaveformParamsInsertPhenomXHMThresholdMband(params->params, atof(argv[++i]));
        }else if(strcmp(argv[i], "--phenomXPHMMband") == 0){
            XLALSimInspiralWaveformParamsInsertPhenomXPHMThresholdMband(params->params, atof(argv[++i]));
        }else if(strcmp(argv[i], "--phenomXPHMUseModes") == 0){
            XLALSimInspiralWaveformParamsInsertPhenomXPHMThresholdMband(params->params, atoi(argv[++i]));
        }else if(strcmp(argv[i], "--phenomXPrecVersion") == 0){
            XLALSimInspiralWaveformParamsInsertPhenomXPrecVersion(params->params, atoi(argv[++i]));
        }else if(strcmp(argv[i], "--phenomXPFinalSpinMod") == 0){
            XLALSimInspiralWaveformParamsInsertPhenomXPFinalSpinMod(params->params, atoi(argv[++i]));
        }else if(strcmp(argv[i], "--phenomXHMRelease") == 0){
            XLALSimInspiralWaveformParamsInsertPhenomXHMReleaseVersion(params->params, atoi(argv[++i]));
        /* New interface specific */
        } else if (strcmp(argv[i], "--total_mass") == 0) {
            XLALSimInspiralWaveformParamsInsertTotalMass(params->params, atof(argv[++i]) * LAL_MSUN_SI);
        } else if (strcmp(argv[i], "--chirp_mass") == 0) {
            XLALSimInspiralWaveformParamsInsertChirpMass(params->params, atof(argv[++i]) * LAL_MSUN_SI);
        } else if (strcmp(argv[i], "--mass_difference") == 0) {
            XLALSimInspiralWaveformParamsInsertMassDifference(params->params, atof(argv[++i]) * LAL_MSUN_SI);
        } else if (strcmp(argv[i], "--reduced_mass") == 0) {
            XLALSimInspiralWaveformParamsInsertReducedMass(params->params, atof(argv[++i]) * LAL_MSUN_SI);
        } else if (strcmp(argv[i], "--spin1_norm") == 0) {
            XLALSimInspiralWaveformParamsInsertSpin1norm(params->params, atof(argv[++i]));
        } else if (strcmp(argv[i], "--spin1_tilt") == 0) {
            XLALSimInspiralWaveformParamsInsertSpin1tilt(params->params, atof(argv[++i]));
        } else if (strcmp(argv[i], "--spin1_phi") == 0) {
            XLALSimInspiralWaveformParamsInsertSpin1phi(params->params, atof(argv[++i]));
        } else if (strcmp(argv[i], "--spin2_norm") == 0) {
            XLALSimInspiralWaveformParamsInsertSpin2norm(params->params, atof(argv[++i]));
        } else if (strcmp(argv[i], "--spin2_tilt") == 0) {
            XLALSimInspiralWaveformParamsInsertSpin2tilt(params->params, atof(argv[++i]));
        } else if (strcmp(argv[i], "--spin2_phi") == 0) {
            XLALSimInspiralWaveformParamsInsertSpin2phi(params->params, atof(argv[++i]));
        } else if (strcmp(argv[i], "--condition") == 0) {
            params->condition = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--python-module") == 0) {
            snprintf(params->module, sizeof(params->module), "%s", argv[++i]);
        } else if (strcmp(argv[i], "--python-object") == 0) {
            snprintf(params->object, sizeof(params->object), "%s", argv[++i]);
        }else {
            XLALPrintError("Error: invalid option: %s\n", argv[i]);
            goto fail;
        }
    }

    return params;

    fail:
    printf("%s", usage);
    XLALFree(params);
    exit(1);
}

/* Function to "unwind" a phase variable with a branch cut */
static int unwind_phase(REAL8 phiUW[], REAL8 phi[],
        size_t len, REAL8 thresh) {
    int cnt = 0; // # of times wrapped around branch cut
    size_t i;
    phiUW[0] = phi[0];
    for(i=1; i<len; i++) {
        if(phi[i-1] - phi[i] > thresh) // phase wrapped forward
            cnt += 1;
        else if(phi[i] - phi[i-1] > thresh) // phase wrapped backwards
            cnt -= 1;
        phiUW[i] = phi[i] + cnt * LAL_TWOPI;
    }
    return 0;
}

static int dump_FD(FILE *f, COMPLEX16FrequencySeries *hptilde,
        COMPLEX16FrequencySeries *hctilde) {
    size_t i;
    COMPLEX16 *dataPtr1 = hptilde->data->data;
    COMPLEX16 *dataPtr2 = hctilde->data->data;
    if (hptilde->data->length != hctilde->data->length) {
        XLALPrintError("Error: hptilde and hctilde are not the same length\n");
        return 1;
    } else if (hptilde->deltaF != hctilde->deltaF) {
        XLALPrintError("Error: hptilde and hctilde do not have the same freq. bin size\n");
        return 1;
    }

    fprintf(f, "# f hptilde.re hptilde.im hctilde.re hctilde.im\n");
    for (i=0; i < hptilde->data->length; i++)
        fprintf(f, "%.16e %.16e %.16e %.16e %.16e\n",
                hptilde->f0 + i * hptilde->deltaF,
                creal(dataPtr1[i]), cimag(dataPtr1[i]), creal(dataPtr2[i]), cimag(dataPtr2[i]));
    return 0;
}

static int dump_FD2(FILE *f, COMPLEX16FrequencySeries *hptilde,
        COMPLEX16FrequencySeries *hctilde) {
    size_t i;
    REAL8 threshold=5.; // Threshold to determine phase wrap-around
    COMPLEX16 *dataPtr1 = hptilde->data->data;
    COMPLEX16 *dataPtr2 = hctilde->data->data;
    if (hptilde->data->length != hctilde->data->length) {
        XLALPrintError("Error: hptilde and hctilde are not the same length\n");
        return 1;
    } else if (hptilde->deltaF != hctilde->deltaF) {
        XLALPrintError("Error: hptilde and hctilde do not have the same freq. bin size\n");
        return 1;
    }
    REAL8 amp1[hptilde->data->length], amp2[hptilde->data->length];
    REAL8 phase1[hptilde->data->length], phase2[hptilde->data->length];
    REAL8 phaseUW1[hptilde->data->length], phaseUW2[hptilde->data->length];
    for (i=0; i < hptilde->data->length; i++)
    {
        amp1[i] = sqrt(creal(dataPtr1[i])*creal(dataPtr1[i])
                + cimag(dataPtr1[i])*cimag(dataPtr1[i]));
        phase1[i] = atan2(cimag(dataPtr1[i]), creal(dataPtr1[i]));
        amp2[i] = sqrt(creal(dataPtr2[i])*creal(dataPtr2[i])
                + cimag(dataPtr2[i])*cimag(dataPtr2[i]));
        phase2[i] = atan2(cimag(dataPtr2[i]), creal(dataPtr2[i]));
    }
    unwind_phase(phaseUW1, phase1, hptilde->data->length, threshold);
    unwind_phase(phaseUW2, phase2, hptilde->data->length, threshold);

    fprintf(f, "# f amp_+ phase_+ amp_x phase_x\n");
    for (i=0; i < hptilde->data->length; i++)
        fprintf(f, "%25.16e %25.16e %25.16e %25.16e %25.16e\n",
                hptilde->f0 + i * hptilde->deltaF,
                amp1[i], phaseUW1[i], amp2[i], phaseUW2[i]);
    return 0;
}

static int dump_TD(FILE *f, REAL8TimeSeries *hplus, REAL8TimeSeries *hcross) {
    size_t i;
    REAL8 t0 = XLALGPSGetREAL8(&(hplus->epoch));
    if (hplus->data->length != hcross->data->length) {
        XLALPrintError("Error: hplus and hcross are not the same length\n");
        return 1;
    } else if (hplus->deltaT != hcross->deltaT) {
        XLALPrintError("Error: hplus and hcross do not have the same sample rate\n");
        return 1;
    }

    fprintf(f, "# t hplus hcross\n");
    for (i=0; i < hplus->data->length; i++)
        fprintf(f, "%25.16e %25.16e %25.16e\n", t0 + i * hplus->deltaT,
                hplus->data->data[i], hcross->data->data[i]);
    return 0;
}

static int dump_TD2(FILE *f, REAL8TimeSeries *hplus, REAL8TimeSeries *hcross) {
    size_t i;
    REAL8 t0 = XLALGPSGetREAL8(&(hplus->epoch));
    REAL8 threshold=5.; // Threshold to determine phase wrap-around
    REAL8 *dataPtr1 = hplus->data->data;
    REAL8 *dataPtr2 = hcross->data->data;
    if (hplus->data->length != hcross->data->length) {
        XLALPrintError("Error: hplus and hcross are not the same length\n");
        return 1;
    } else if (hplus->deltaT != hcross->deltaT) {
        XLALPrintError("Error: hplus and hcross do not have the same sample rate\n");
        return 1;
    }
    REAL8 amp[hplus->data->length];
    REAL8 phase[hplus->data->length];
    REAL8 phaseUW[hplus->data->length];
    for (i=0; i < hplus->data->length; i++)
    {
        amp[i] = sqrt( dataPtr1[i]*dataPtr1[i] + dataPtr2[i]*dataPtr2[i]);
        phase[i] = atan2(-dataPtr2[i], dataPtr1[i]);
    }
    unwind_phase(phaseUW, phase, hplus->data->length, threshold);

    fprintf(f, "# t amp phase\n");
    for (i=0; i < hplus->data->length; i++)
        fprintf(f, "%25.16e %25.16e %25.16e\n", t0 + i * hplus->deltaT,
                amp[i], phaseUW[i]);
    return 0;
}

/*
 * main
 */
int main (int argc , char **argv) {
    FILE *f;
    int status;
    int start_time;
    COMPLEX16FrequencySeries *hptilde = NULL, *hctilde = NULL;
    REAL8TimeSeries *hplus = NULL;
    REAL8TimeSeries *hcross = NULL;
    GSParams *params;

    /* set us up to fail hard */
    XLALSetErrorHandler(XLALAbortErrorHandler);

    /* parse commandline */
    params = parse_args(argc, argv);

    /* generate waveform */
    start_time = time(NULL);
    if (params->newinterface == 0) // Use old interface
    {
        switch (params->domain) {
            case LAL_SIM_DOMAIN_FREQUENCY:
                XLALSimInspiralChooseFDWaveform(&hptilde, &hctilde,
                        params->m1, params->m2, params->s1x,
                        params->s1y, params->s1z, params->s2x, params->s2y,
                        params->s2z, params->distance, params->inclination,
                        params->phiRef, params->longAscNodes, params->ecc, params->meanPerAno,
                        params->deltaF, params->f_min, params->f_max, params->fRef,
                        params->params, params->approximant);
                break;
            case LAL_SIM_DOMAIN_TIME:
                XLALSimInspiralChooseTDWaveform(&hplus, &hcross,
                        params->m1, params->m2, params->s1x,
                        params->s1y, params->s1z, params->s2x, params->s2y,
                        params->s2z, params->distance, params->inclination,
                        params->phiRef, params->longAscNodes, params->ecc, params->meanPerAno,
                        params->deltaT, params->f_min, params->fRef,
                        params->params,
                        params->approximant);
                break;
            default:
                XLALPrintError("Error: domain must be either TD or FD\n");
        }
    }
    else // Use new interface
    {
        /* Generator */
        LALDict *generator_params = XLALCreateDict();
        if (params->condition){
          XLALDictInsertINT4Value(generator_params, "condition", params->condition);
        }
        if (strcmp(params->module, "") != 0){
            XLALDictInsertStringValue(generator_params, "module", params->module);
        }
        if (strcmp(params->object, "") != 0){
          XLALDictInsertStringValue(generator_params, "object", params->object);
          XLALDictInsertINT4Value(generator_params, "approximant", IMRPhenomXPHM);
        }
        LALSimInspiralGenerator *generator = XLALSimInspiralChooseGenerator(params->approximant, generator_params);
        XLAL_CHECK(generator, XLAL_EFUNC);

        /* Parse waveform params into laldictionary */
        /* If none of the optional mass arguments are used, we use the component masses which are stored in the params object */
        const char *mass_keys[6] = {"total_mass", "chirp_mass", "mass_difference", "reduced_mass", "mass_ratio", "sym_mass_ratio"};
        UINT2 key_exists = 0;
        for (size_t j = 0; j < sizeof(mass_keys)/sizeof(*mass_keys); ++j){
          if (XLALDictContains(params->params, mass_keys[j]) == 1) {
            key_exists = 1;
            break;
          }
        }
        if (key_exists == 0){
          XLALSimInspiralWaveformParamsInsertMass1(params->params, params->m1);
          XLALSimInspiralWaveformParamsInsertMass2(params->params, params->m2);
        }
        else{
          /* Insert component masses if they are in the command line */
          if (params->m1_newinterface == 1) XLALSimInspiralWaveformParamsInsertMass1(params->params, params->m1);
          if (params->m2_newinterface == 1) XLALSimInspiralWaveformParamsInsertMass2(params->params, params->m2);
        }
        /* If any of the optional spin arguments are used, we use the cartesian spins which are stored in the params object */
        if (XLALDictContains(params->params, "spin1_norm") == 0 && XLALDictContains(params->params, "spin1_tilt") == 0 && XLALDictContains(params->params, "spin1_phi") == 0){
            XLALSimInspiralWaveformParamsInsertSpin1x(params->params, params->s1x);
            XLALSimInspiralWaveformParamsInsertSpin1y(params->params, params->s1y);
            XLALSimInspiralWaveformParamsInsertSpin1z(params->params, params->s1z);
        }
        if (XLALDictContains(params->params, "spin2_norm") == 0 && XLALDictContains(params->params, "spin2_tilt") == 0 && XLALDictContains(params->params, "spin2_phi") == 0){
            XLALSimInspiralWaveformParamsInsertSpin2x(params->params, params->s2x);
            XLALSimInspiralWaveformParamsInsertSpin2y(params->params, params->s2y);
            XLALSimInspiralWaveformParamsInsertSpin2z(params->params, params->s2z);
        }
        /* Insert the rest of waveform params */
        XLALSimInspiralWaveformParamsInsertDeltaF(params->params, params->deltaF);
        XLALSimInspiralWaveformParamsInsertDeltaT(params->params, params->deltaT);
        XLALSimInspiralWaveformParamsInsertF22Ref(params->params, params->fRef);
        XLALSimInspiralWaveformParamsInsertF22Start(params->params, params->f_min);
        XLALSimInspiralWaveformParamsInsertFMax(params->params, params->f_max);
        XLALSimInspiralWaveformParamsInsertRefPhase(params->params, params->phiRef);
        XLALSimInspiralWaveformParamsInsertDistance(params->params, params->distance);
        XLALSimInspiralWaveformParamsInsertInclination(params->params, params->inclination);
        XLALSimInspiralWaveformParamsInsertLongAscNodes(params->params, params->longAscNodes);
        XLALSimInspiralWaveformParamsInsertEccentricity(params->params, params->ecc);
        XLALSimInspiralWaveformParamsInsertMeanPerAno(params->params, params->meanPerAno);

        /* generate waveform */
        switch (params->domain) {
          case LAL_SIM_DOMAIN_FREQUENCY:
              XLALSimInspiralGenerateFDWaveform(&hptilde, &hctilde, params->params, generator);
              break;
          case LAL_SIM_DOMAIN_TIME:
              XLALSimInspiralGenerateTDWaveform(&hplus, &hcross, params->params, generator);
              break;
          default:
              XLALPrintError("Error: domain must be either TD or FD\n");
        }
        XLALDestroyDict(generator_params);
        XLALDestroySimInspiralGenerator(generator);
    }
    if (params->verbose)
        XLALPrintInfo("Generation took %.0f seconds\n",
                difftime(time(NULL), start_time));
    if (((params->domain == LAL_SIM_DOMAIN_FREQUENCY) && (!hptilde || !hctilde)) ||
        ((params->domain == LAL_SIM_DOMAIN_TIME) && (!hplus || !hcross))) {
        XLALPrintError("Error: waveform generation failed\n");
        goto fail;
    }

    /* dump file */
    if ( strlen(params->outname) > 0 ) {
      f = fopen(params->outname, "w");
      if (f==NULL) {
        printf("**ERROR** Impossible to write file %s\n",params->outname);
        exit(1);
      }
      else {
        if (params->domain == LAL_SIM_DOMAIN_FREQUENCY)
          if (params->ampPhase == 1)
            status = dump_FD2(f, hptilde, hctilde);
          else
            status = dump_FD(f, hptilde, hctilde);
        else
          if (params->ampPhase == 1)
            status = dump_TD2(f, hplus, hcross);
          else
            status = dump_TD(f, hplus, hcross);
        fclose(f);
      }
      if (status) goto fail;
    }

    /* clean up */
    XLALDestroyDict(params->params);
    XLALDestroyREAL8TimeSeries(hplus);
    XLALDestroyREAL8TimeSeries(hcross);
    XLALFree(params);

    XLALDestroyCOMPLEX16FrequencySeries(hptilde);
    XLALDestroyCOMPLEX16FrequencySeries(hctilde);
    return 0;

    fail:
    XLALDestroyDict(params->params);
    XLALFree(params);
    XLALDestroyCOMPLEX16FrequencySeries(hptilde);
    XLALDestroyCOMPLEX16FrequencySeries(hctilde);
    return 1;
}

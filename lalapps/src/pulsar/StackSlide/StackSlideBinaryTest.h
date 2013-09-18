/*
*  Copyright (C) 2007 Virginia Re
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

/* REVISIONS: */
/* 04/26/04 gam; Change LALStackSlide to StackSlide and LALSTACKSLIDE to STACKSLIDE for initial entry to LALapps. */
/* 06/05/04 gam; Add gpsStartTimeSec and gpsStartTimeNan to StackSlideSkyParams; set these to epoch that gives T0 at SSB. */

/**
 * \author  Landry, M., and Mendell, G.
 *
 * ### Header \ref ComputeSky.h ###
 *
 * Computes phase coefficients necessary for a correct demodulation.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/ComputeSky.h>
 * \endcode
 *
 * This is a short summary of the analytical calculations which form the basis for the code in this routine.
 *
 * Recall that a demodulated Fourier Transform (DeFT) is given by
 * \f{equation}{
 * \label{e4}
 * \hat{x}_b({\vec{\lambda}})=
 * \sum_{\alpha =0}^{M-1}\sum_{k=0}^{N-1}\tilde{x}_{\alpha k}\left[\frac{1}{N}\sum_{j=0}^{N-1}e^{-2\pi i(\Phi_{\alpha jb}(\vec{\lambda})-\frac{jk}{N})}\right]
 * \f}
 * The index \f$b\f$ defines the DeFT frequency bin, the index \f$\alpha\f$ loops through
 * the SFTs that build the DeFT, \f$k\f$ runs on all the SFT frequency bins, and \f$j\f$
 * is a time index that runs on each SFT.  As shown in section
 * \TODOref{s_LALDemod_h}, the next step in the development of the demodulation
 * technique involves Taylor expanding the phase model about the temporal
 * midpoint of each short segment of data, while retaining only first order
 * terms.  The Taylor expansion of \f$\Phi (t)\f$ about the temporal midpoint
 * \f$t_{\alpha,1/2}\f$ is
 * \f{equation}{
 * \label{taylor2}
 * \Phi_{\alpha}(t) = \Phi(t_{\alpha,1/2})+\left[t-t_{\alpha,1/2}\right]\frac{d\Phi}{dt}(t_{\alpha,1/2}) \\
 * \f}
 * For each value of \f$\alpha\f$, this expression consist of either constant or linear terms in time.  With the particular time discretization chosen in this code, \f$t=t_{0}+(N\alpha+j)\ T_{obs}/NM\f$, we have
 * \f{equation}{
 * \label{time}
 * \left[t-t_{\alpha,1/2}\right]=\frac{\ T_{obs}}{M}\left(\frac{j}{N}-\frac{1}{2}\right)=\mathcal{T}_{s}\left(\frac{j}{N}-\frac{1}{2}\right),
 * \f}
 * where \f$\mathcal{T}_{s}\f$ is the short time baseline of the \f$M\f$ short FTs.  On
 * the other hand, the phase can also be expressed as a function of SSB time \f$T\f$
 * (i.e. the time at the solar system barycenter).  We will assume the source to
 * be at rest in this reference frame.  Now, if one adopts the notation \f$\Delta
 * T_{\alpha}\equiv\left[T(t_{\alpha,1/2})-
 * T(t_{0})\right]\f$ and \f$\dot{T}_{\alpha}\equiv
 * dT/dt(t_{\alpha,1/2})\f$
 * the phase terms in the above equation are (neglecting constants)
 * \f{eqnarray}{
 * \label{phi}
 * \Phi(t_{\alpha,1/2})                     & = & f_{0}\Delta T_{\alpha}+\frac{1}{2}f_{1}\Delta T_{\alpha}^{2}
 * +\frac{1}{3}f_{2}\Delta T_{\alpha}^{3}+\frac{1}{4}f_{3}\Delta T_{\alpha}^{4}+\frac{1}{5}f_{4}\Delta T_{\alpha}^{5}
 * +\frac{1}{6}f_{5}\Delta T_{\alpha}^{6} \\
 * &   & \\
 * \label{dphi}
 * \frac{d\Phi}{dt}(t_{\alpha,1/2})         & = & \dot{T}_{\alpha}\left(f_{0}+ f_{1}\Delta T_{\alpha}
 * +f_{2}\Delta T_{\alpha}^{2}+f_{3}\Delta T_{\alpha}^{3}
 * +f_{4}\Delta T_{\alpha}^{4}+f_{5}\Delta T_{\alpha}^{5}\right).
 * \f}
 * These constants, for each value of \f$\alpha\f$, require \f$\dot{T}_{\alpha}\f$ and
 * \f$\Delta T_{\alpha}\f$, which are calculated by a suitable timing routine.  For
 * this demodulation package, this timing routine is provided by <tt>tdb()</tt>.
 * Thus, for a given sky position, the timing routine will be called once for
 * each short time chunk, each call returning a specific  \f$\dot{T}_{\alpha}\f$ and
 * \f$\Delta T_{\alpha}\f$.  By substituting Eq.s\TODOref{time},\TODOref{phi} and
 * \TODOref{dphi} in \eqref{taylor2} and grouping together the terms in \f$j\f$ (linear
 * in \f$t\f$) in order to save computations, we have
 * \f{equation}{
 * \label{phasecalc}
 * \Phi_{\alpha}(t)=\sum_{s=0}^{n_{spin}}f_{s}A_{s\alpha}+\frac{j}{N}\sum_{s=0}^{n_{spin}}f_{s}B_{s\alpha},
 * \f}
 * where \f$n_{spin}\f$ is the maximum order of spindown parameter.  Rather than
 * store the values of \f$\dot{T}_{\alpha}\f$ and \f$\Delta T_{\alpha}\f$ for each value
 * of \f$\alpha\f$, it is more efficient to calculate the constants \f$A_{s\alpha}\f$ and
 * \f$B_{s\alpha}\f$ only once, and then use these values for every spindown
 * parameter set used when searching in a given sky position.  Analytical
 * formulae for these constants are easily derived:
 * \f{equation}{
 * A_{s \alpha}=\frac{1}{s+1}\Delta T_{\alpha}^{s+1}-\frac{1}{2}\mathcal{T}_{SFT}\dot{T}_{\alpha}\Delta T_{\alpha}^{s}
 * \f}
 * \f{equation}{
 * B_{s \alpha}=\mathcal{T}_{SFT}\dot{T}_{\alpha}\Delta T_{\alpha}^{s}
 * \f}
 *
 */

#ifndef _STACKSLIDE_H
#define _STACKSLIDE_H

#include <stdio.h>
#include <stdlib.h>
#include <lal/FileIO.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <math.h>
#include <string.h> 
#include <lal/LALConstants.h>
#include <lal/StreamInput.h>
#include <lal/SeqFactories.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/Random.h>
#include <getopt.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/LALDatatypes.h>
#include <lal/FindRoot.h>
#include <lal/DetResponse.h>
#include <lal/DetectorSite.h>
#include <lal/VectorOps.h>


#ifdef __cplusplus
extern "C" {
#endif

  
/* Author-defined error codes */

/**\name Error Codes */ /*@{*/
#define STACKSLIDEH_ENULL 1
#define STACKSLIDEH_ENNUL 2
#define STACKSLIDEH_ENEGA 4
#define STACKSLIDEH_MSGENULL "Null Pointer"
#define STACKSLIDEH_MSGENNUL "Non-Null Pointer"
#define STACKSLIDEH_MSGENEGA "Bad Negative Value"
#define STACKSLIDECOMPUTESKYH_ENULL 6
#define STACKSLIDECOMPUTESKYH_ENNUL 8
#define STACKSLIDECOMPUTESKYH_ENEGA 10
#define STACKSLIDECOMPUTESKYH_MSGENULL "Null Pointer in StackSlideComputeSky"
#define STACKSLIDECOMPUTESKYH_MSGENNUL "Non-Null Pointer in StackSlideComputeSky"
#define STACKSLIDECOMPUTESKYH_MSGENEGA "Bad Negative Value in StackSlideComputeSky"
/*@}*/

/**
 *
 * ### Structures ###
 *
 * \code
 * struct CSParams
 * \endcode
 * \c CSParams
 *
 * This structure contains the parameters for the <tt>ComputeSky()</tt> routine.  The parameters are:
 *
 * <dl>
 * <dt><tt>INT8 spinDwnOrder</tt></dt><dd> The maximal number of spindown parameters per spindown parameter set.</dd>
 * <dt><tt>INT8 mObsSFT</tt></dt><dd> The number of SFTs in the observation time.</dd>
 * <dt><tt>REAL8 tSFT</tt></dt><dd> The timescale of one SFT.</dd>
 * <dt><tt>LIGOTimeGPS *tGPS</tt></dt><dd> An array containing the GPS times of the first datum from each SFT.</dd>
 * <dt><tt>REAL8 *skyPos</tt></dt><dd> The array containing the sky patch coordinates.</dd>
 * <dt><tt>CHAR *sw</tt></dt><dd> A switch which turns modulation on/off. </dd>
 * <dt><tt>void (*funcName)(REAL8 , REAL8 , REAL8 , REAL8 *, REAL8 *, const CHAR *sw)</tt></dt><dd> A function pointer, to make the use of different timing routines easy.</dd>
 * </dl>
 *
 */

/**\name Error Codes */ /*@{*/
#define COMPUTESKYBINARYH_ENULL 1
#define COMPUTESKYBINARYH_ENNUL 2
#define COMPUTESKYBINARYH_ERANG 3
#define COMPUTESKYBINARYH_ENEGA 4
#define COMPUTESKYBINARYH_MSGENULL "Null Pointer"
#define COMPUTESKYBINARYH_MSGENNUL "Non-Null Pointer"
#define COMPUTESKYBINARYH_MSGERANG "Input parameter out of range"
#define COMPUTESKYBINARYH_MSGENEGA "Bad Negative Value"
/*@}*/

  
#define ACC 1e-9

typedef struct
tagStackSlideParams /* substituted tagStackSlideBinaryParams*/
{
	REAL8 **skyPosData;  
	REAL8 **freqDerivData;  
	INT4 numSkyPosTotal;
	INT4 numFreqDerivTotal;
	REAL8 f0STK;
	REAL8 f0SUM;
	REAL8 tSTK;
	REAL8 tSUM;
	INT4  nBinsPerSUM;
	INT4  numSTKs;
	REAL8 dfSUM;
	UINT4 gpsStartTimeSec;
	UINT4 gpsStartTimeNan;
	INT4 numSpinDown;
	EphemerisData *edat;
	/* INT4 patchNumber; */
	BarycenterInput *baryinput;
	}
StackSlideParams; /*substituted StackSlideBinaryParams*/



typedef struct
tagStackSlideBinarySkyParams
{
	INT8		spinDwnOrder;	/* max spindown parameter order */
	INT8		mObsSFT;	/* number of coherent timescales */
	REAL8		tSFT;		/* timescale of SFT */
	/*LIGOTimeGPS	*tGPS;	*/	/* GPS time of 1st data sample of each SFT */
	UINT4 gpsStartTimeSec;          /* 06/05/04 gam; set these to epoch that gives T0 at SSB. */
	UINT4 gpsStartTimeNan;          /* 06/05/04 gam; set these to epoch that gives T0 at SSB. */
	REAL8 ap0; /*central value for sma from command line*/
	INT4 Tperi0; /*central value for periapse time passage*/
	REAL8 		*skyPos; 	/* array of sky positions */
        REAL8 		SemiMajorAxis;  /* orbital radius of binary (in sec) */
        REAL8           OrbitalPeriod;         /* Period of binary (in sec) */
        REAL8           OrbitalEccentricity; /* Orbital eccentricy */
        REAL8           ArgPeriapse;    /* Argument of Periapse */
        LIGOTimeGPS     TperiapseSSB;   /* Instance of periapse passage measured in the SSB frame */
	UINT4 gpsTperiapseSSBSec;
	BarycenterInput *baryinput;	
	EmissionTime *emit;
	EarthState *earth;
	EphemerisData *edat;
	REAL8 alpha;
	REAL8 delta;
	REAL8 dInv;
       /* CHAR *sunEdatFile;
	CHAR *earthEdatFile;*/
}
StackSlideBinarySkyParams;

typedef struct
tagTdotsAndDeltaTs
{
	REAL8		*vecTDots;	/* 1-d array of (dT/dt)'s for frequency calculation */
	REAL8		**vecDeltaTs;	/* 2-d array of (T-T_0)'s for frequency calculation */
}
TdotsAndDeltaTs;

int ReadCommandLineArgs(int argc, char *argv[]);

void StackSlideBinary(	LALStatus *status, 
			/*REAL4FrequencySeries **SUMData, */
			/*REAL4FrequencySeries **STKData,*/
		        REAL8 **SUMData,
			REAL8 **fakeSFT,
			TdotsAndDeltaTs *pTdotsAndDeltaTs,
			StackSlideParams *SSparams);

int SetupBaryInput();
	
int GetSSBTime(LIGOTimeGPS *tdet, LIGOTimeGPS *tssb);


void StackSlideComputeSkyBinary (LALStatus 	*status, 
			TdotsAndDeltaTs 	*pTdotsAndDeltaTs, 
			INT8 		iSkyCoh, 
			StackSlideBinarySkyParams 	*params);
		

int FreeMem(void);
#ifdef __cplusplus
}
#endif

#endif /* _STACKSLIDE_H */

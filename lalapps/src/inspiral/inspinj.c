/*
 *  Copyright (C) 2007 Chad Hanna, Alexander Dietz, Duncan Brown, Gareth Jones, Jolien Creighton, Nickolas Fotopoulos, Patrick Brady, Stephen Fairhurst, Tania Regimbau, Salvatore Vitale
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



/*-----------------------------------------------------------------------
 *
 * File Name: inspinj.c
 *
 * Author: Brown, D. A., Creighton, J. D. E. and Dietz A. IPN contributions from Predoi, V.
 *
 *
 *-----------------------------------------------------------------------
 */

/**
 * \file
 * \ingroup lalapps_inspiral
 *
 * <dl>
 * <dt>Name</dt><dd>
 * \c lalapps_inspinj --- produces inspiral injection data files.</dd>
 *
 * <dt>Synopsis</dt><dd>
 * \c lalapps_inspinj
 *
 * [<tt>--help</tt>]
 * <tt>--source-file</tt> \c sfile
 * <tt>--mass-file</tt> \c mfile
 *
 * [<tt>--gps-start-time</tt> \c tstart]
 * [<tt>--gps-end-time</tt> \c tend]
 *
 * [<tt>--time-step</tt> \c tstep]
 * [<tt>--time-interval</tt> \c tinterval]
 *
 * [<tt>--seed</tt> \c seed]
 * [<tt>--waveform</tt> \c wave]
 * [<tt>--lal-eff-dist</tt>]
 * [<tt>--usertag</tt> \c tag]
 *
 * [<tt>--tama-output</tt>]
 * [<tt>--write-eff-dist</tt>]
 *
 * [<tt>--ilwd</tt>]</dd>
 *
 * <dt>Description</dt><dd>
 * \c lalapps_inspinj
 * generates a number of inspiral  parameters suitable  for using in a Monte
 * Carlo injection to test the efficiency of a inspiral search.  The  various
 * parameters (detailed  below)  are randomly chosen and are appropriate for a
 * particular population of binary neutron stars  whose spatial  distribution
 * includes the Milky Way and a number of extragalactic objects that are  input
 * in  a  datafile.  The  possible  mass pairs for the binary neutron star com-
 * panions are also specified in a (different) datafile.
 *
 * The output of this program  is  a  list  of  the  injected events,  starting
 * at  the specified start time and ending at the specified end time.  One
 * injection with random inspiral parameters will be made every specified time
 * step, and will be randomly placed within the specified time interval.
 * The output is written to a file name in the standard inspiral pipeline format:
 *
 * \code
 * HL-INJECTIONS_USERTAG_SEED-GPSSTART-DURATION.xml
 * \endcode
 *
 * where \c USERTAG is \c tag as specfied on the command line,
 * \c SEED is the  value  of  the random number seed chosen and
 * \c GPSSTART and \c DURATION describes the GPS time interval that
 * the file covers. The file is in the standard LIGO lightweight XML format
 * containing a \c sim_inspiral table that describes the injections.
 * In addition, an ASCII log file called <tt>injlog.txt</tt> is also written.
 * If a <tt>--user-tag</tt> is not specified on the command line, the
 * \c _USERTAG part of the filename will be omitted.</dd>
 *
 * <dt>Options</dt><dd>
 * <ul>
 * <li><tt>--help</tt>: Print a help message.</li>
 *
 * <li><tt>--source-file</tt> \c sfile:
 * Optional. Data file containing spatial distribution of  extragalactic  objects.
 * Default  is  the file <tt>inspsrcs.dat</tt> provided by LALApps. If that file is
 * empty, all signals are in the Milky Way.</li>
 *
 * <li><tt>--mass-file</tt> \c mfile:
 * Optional. Data file containing mass pairs  for  the binary  neutron  star
 * companions.   Default is the file <tt>BNSMasses.dat</tt> provided by LALApps.</li>
 *
 * <li><tt>--gps-start-time</tt> \c tstart:
 * Optional.  Start time of the injection data to be created. Defaults to the
 * start of S2, Feb 14 2003 16:00:00 UTC (GPS time 729273613)</li>
 *
 * <li><tt>--gps-end-time</tt> \c tend:
 * Optional. End time of the injection data to be created. Defaults to the end of
 * S2, Apr 14 2003 15:00:00 UTC (GPS time 734367613).</li>
 *
 * <li><tt>--time-step</tt> \c tstep:
 * Optional. Sets the time step interval between injections. The injections will
 * occur with an average spacing of \c tstep seconds. Defaults to
 * \f$2630/\pi\f$.</li>
 *
 * <li><tt>--time-interval</tt> \c tinterval:
 * Optional. Sets the time interval during which an injection can occur.
 * Injections are uniformly distributed over the interval.  Setting \c tstep
 * to \f$6370\f$ and \c tinterval to 600 guarantees there will be one injection
 * into each playground segment and they will be randomly distributed within the
 * playground times - taken the fact that your gps start time coincides with start of a playground segment.</li>
 *
 * <li><tt>--seed</tt> \c seed:
 * Optional. Seed the random number generator with the integer \c seed.
 * Defaults to \f$1\f$.</li>
 *
 * <li><tt>--waveform</tt> \c wave:
 * Optional. The string \c wave will be written into the \c waveform
 * column of the \c sim_inspiral table output. This is used by the
 * inspiral code to determine which type of waveforms it should inject into the
 * data. Defaults is \c GeneratePPNtwoPN.</li>
 *
 * <li><tt>--lal-eff-dist</tt>:
 * Optional.  If this option is specified, the effective distance will be
 * calculated using routines from LAL.  Otherwise, the default behaviour is to
 * use an independent method contained in inspinj.c.  There is good agreement
 * between these two methods, see below for more details.</li>
 *
 * <li><tt>--user-tag</tt> \c string: Optional. Set the user tag for this
 * job to be \c string. May also be specified on the command line as
 * <tt>-userTag</tt> for LIGO database compatibility.</li>
 *
 * <li><tt>--tama-output</tt>:
 * Optional.  If this option is given, \c lalapps_inspinj also produces a
 * text output file:
 *
 * \code
 * HLT-INJECTIONS_USERTAG_SEED-GPSSTART-DURATION.txt
 * \endcode
 *
 * which contains the following fields:
 *
 * <ul>
 * <li> geocentric end time</li>
 * <li> Hanford end time</li>
 * <li> Livingston end time</li>
 * <li> TAMA end time</li>
 * <li> total mass, \f$M_{\mathrm{TOT}}\f$</li>
 * <li> mass ratio, \f$\eta\f$</li>
 * <li> distance to source (in kpc)</li>
 * <li> longitude</li>
 * <li> latitude</li>
 * <li> inclination</li>
 * <li> coalescence phase</li>
 * <li> polarization</li>
 * <li> TAMA polarization</li>
 * <li> end time GMST</li>
 * </ul>
 *
 * In the above, all times are recorded as double precision real numbers and all
 * angles are in radians.  The TAMA polarization is calculated using
 *
 * \f{equation}{
 *   \tan( \psi_{T} ) = \frac{ x \cdot T_{z} }{ y \cdot T_{z} } \, .
 * \f}
 *
 * Here x and y are the x,y axes of the radiation frame expressed in earth fixed
 * coordinates \eqref{xrad}, \eqref{yrad}.  \f$T_{z}\f$ is a unit vector in earth fixed
 * coordinates which is orthogonal to the two arms of the TAMA detector
 * \eqref{tarm}.  It is given by
 *
 * \f{equation}{
 *   T_{z} = ( -0.6180, +0.5272, +0.5832 )
 * \f}</li>
 *
 * <li><tt>--write-eff-dist</tt>: Optional.  If this option is given, three extra
 * columns are added to the TAMA output file described above.  They are
 * <ul>
 * <li> Hanford effective distance (kpc)</li>
 * <li> Livingston effective distance (kpc)</li>
 * <li> TAMA effective distance (kpc)</li>
 * </ul>
 *
 * These entries are added to the list immediately after TAMA end time and before
 * total mass.</li>
 *
 * <li><tt>--ilwd</tt>: Optional. If this option is given,
 * \c lalapps_inspinj also produces two ILWD-format files, injepochs.ilwd and
 * injparams.ilwd, that contain, respectively, the  GPS  times  suitable for
 * inspiral injections, and the intrinsic inspiral signal parameters to be used
 * for  those injections.
 *
 * The  file  injepochs.ilwd  contains  a sequence of integer pairs representing
 * the injection GPS time in  seconds  and residual  nano-seconds.   The file
 * injparams.ilwd contains the intrinsic binary parameters for each injection,
 * which is  a  sequence  of  eight  real  numbers representing (in order) (1) the
 * total mass of the binary system  (in  solar masses),  (2)  the  dimensionless
 * reduced mass --- reduced mass per unit total mass --- in the range from  0
 * (extreme mass  ratio)  to  0.25 (equal masses), (3) the distance to the system
 * in meters, (4) the inclination  of  the  binary system  orbit  to the plane of
 * the sky in radians, (5) the coalescence phase in radians, (6)  the  longitude
 * to  the direction  of  the  source in radians, (7) the latitude to the
 * direction of the source in radians, (8) and the polar- ization angle of the
 * source in radians.</li>
 * </ul></dd>
 *
 * <dt>Example</dt><dd>
 * \code
 * lalapps_inspinj --seed 45\
 * --source-file inspsrcs.dat --mass-file BNSMasses.dat
 * \endcode</dd>
 *
 * <dt>Algorithm</dt><dd>
 *
 * The algorithm for computing the effective distance will be described in some
 * detail below.  The method is to compute both the strain due to the inspiral
 * and the detector response in the earth fixed frame.  This frame is such that
 * the z-axis points from the earth's centre to the North Pole, the x-axis points
 * from the centre to the intersection of the equator and the prime meridian and
 * the y-axis is chosen to complete the orthonormal basis.  The coordinates of
 * the injection are specified by longitude (or right ascension) \f$\alpha\f$ and
 * latitude (or declination) \f$\delta\f$.  The polarization is appropriate for
 * transferring from the radiation to earth fixed frame.  These are then
 * converted to the earth fixed frame by
 *
 * \f{eqnarray}{
 *   \theta &=& \frac{\pi}{2} - \delta \\
 *   \phi &=& \alpha - \textrm{gmst} \, .
 * \f}
 *
 * Here, gmst is the Greenwich Mean sidereal time of the injection.  The axes of
 * the radiation frame (x,y,z) can be expressed in terms of the earth fixed
 * coordinates as:
 *
 * \f{eqnarray}{
 *   x(1) &=& +( \sin( \phi ) \cos( \psi ) - \sin( \psi ) \cos( \phi )
 *       \cos( \theta ) ) \nonumber \\
 *   x(2) &=& -( \cos( \phi ) \cos( \psi ) + \sin( \psi ) \sin( \phi )
 *       \cos( \theta ) ) \nonumber \\
 *   x(3) &=& \sin( \psi ) \sin( \theta ) \label{xrad}\\
 *   y(1) &=& -( \sin( \phi ) \sin( \psi ) + \cos( \psi ) \cos( \phi )
 *       \cos( \theta ) ) \nonumber\\
 *   y(2) &=& +( \cos( \phi ) \sin( \psi ) - \cos( \psi ) \sin( \phi )
 *       \cos( \theta ) ) \nonumber \\
 *   y(3) &=& \cos( \psi ) \sin( \theta ) \label{yrad}
 * \f}
 *
 * Making use of these expressions, we can express the gravitational wave strain in
 * earth fixed coordinates as
 *
 * \f{equation}{\label{hij}
 *   h_{ij} = ( h^{+}(t) e^{+}_{ij} ) + (h^{\times}(t) e^{\times}_{ij})
 * \f}
 *
 * where
 *
 * \f{equation}{
 *   e^{+}_{ij} = x_{i} * x_{j} - y_{i} * y_{j} \qquad \mathrm{and} \qquad
 *   e^{\times}_{ij} = x_{i} * y_{j} + y_{i} * x_{j}.
 * \f}
 *
 * For the case of a binary inspiral signal, the two polarizations \f$h^{+}\f$
 * and \f$h^{\times}\f$ of the gravitational wave are given by
 *
 * \f{eqnarray}{
 *   h^{+}(t)  &=& \frac{A}{r}  ( 1 + \cos^2 ( \iota ) ) * \cos( \Phi(t) ) \\
 *   h^{\times}(t) &=& \frac{A}{r} * ( 2 \cos( \iota )   ) * \sin( \Phi(t) )
 * \f}
 *
 * where \f$A\f$ is a mass and frequency dependent amplitude factor, \f$r\f$ is the
 * physical distance at which the injection is located and \f$\iota\f$ is the
 * inclination angle.
 *
 * Next, we can write the detector response function as
 *
 * \f{equation}{
 *   d^{ij} = \left(\frac{1}{2} \right) \left( n_{x}^{i} n_{x}^{j}
 *       - n_{y}^{i} n_{y}^{j} \right) \, .
 * \f}
 *
 * Here, \f$n_{x}\f$ and \f$n_{y}\f$ are unit vectors directed along the arms of the
 * detector.  Specifically, for the Hanford, Livingston, GEO, TAMA and Virgo
 * detectors we use:
 *
 * \f{eqnarray}{
 *   H_{x} &=& ( -0.2239, +0.7998, +0.5569 ) \nonumber \\
 *   H_{y} &=& ( -0.9140, +0.0261, -0.4049 ) \\
 *   L_{x} &=& ( -0.9546, -0.1416, -0.2622 ) \nonumber \\
 *   L_{y} &=& ( +0.2977, -0.4879, -0.8205 ) \\
 *   G_{x} &=& ( -0.6261, -0.5522, +0.5506 ) \nonumber \\
 *   G_{y} &=& ( -0.4453, +0.8665, +0.2255 ) \\
 *   T_{x} &=& ( +0.6490, +0.7608, +0.0000 ) \nonumber \\
 *   T_{y} &=& ( -0.4437, +0.3785, -0.8123 ) \label{tarm} \\
 *   V_{x} &=& ( -0.7005, +0.2085, +0.6826 ) \nonumber \\
 *   V_{y} &=& ( -0.0538, -0.9691, +0.2408 )
 * \f}
 *
 * The response of an interferometric detector with arm locations given by \f$n_{x}\f$
 * and \f$n_{y}\f$ to an inspiralling binary system described by \eqref{hij} is
 *
 * \f{eqnarray}{
 *   h(t) &=& h^{+}(t) ( d^{ij} e^{+}_{ij} )
 *     + h^{\times}(t) ( d^{ij} e^{\times}_{ij} ) \nonumber \\
 *       &=&
 *     \left(\frac{A}{r}\right) \left[
 * 	( 1 + \cos^2 ( \iota ) ) F_{+} \cos( \Phi(t)) +
 *         2 \cos( \iota ) F_{\times} \sin( \Phi(t) ) \right] \, ,
 * \f}
 *
 * where we have introduced
 *
 * \f{equation}{
 *   F_{+} = d^{ij} e^{+}_{ij} \qquad \mathrm{and} \qquad
 *   F_{\times} = d^{ij} e^{\times}_{ij}
 * \f}
 *
 * Finally, to calculate the effective distance, we note that the two contributions
 * to \f$h(t)\f$ are \f$\pi/2\f$ radians out of phase, and hence orthogonal.  Thus, we can
 * compute the effective distance to be:
 *
 * \f{equation}{
 *   D_{\mathrm{eff}} = r / \left( \frac{ (1 + \cos^2(\iota))^2 }{4} F_{+}^{2} +
 *       cos^{2}(\iota) F_{\times}^{2} \right)
 * \f}
 *
 * \anchor eff_dist_comparison
 * \image html effective_distance_comparison.png "Comparison of effective distance computed by inspinj.c and LAL routines"
 *
 * The algorithm to calculate effective distances described above is completely
 * contained within inspinj.c.  There is an independent method of computing
 * effective distances can also be called by inspinj.  It is contained in the LAL
 * function <tt>LALPopulateSimInspiralSiteInfo()</tt>.  This function populates
 * the site end time and effective distance for all the interferomter sites.  It
 * makes use of LAL functionality in the tools and date packages.  These same
 * functions are used when generating the injection waveform which is added to
 * the data stream (in lalapps_inspiral).  As a check that these two
 * calculations produce the same effective distance, lalapps_inspinj was run
 * twice, once with the <tt>--lal-eff-dist</tt> option and once without.
 * \ref eff_dist_comparison "This figure" shows the fractional difference in effective
 * distance between the two methods for a set of injections.  We see that the
 * distances agree within 1\
 * occuring for the largest effective disances, i.e.  close to the dead spot of
 * the instrument.  For injections which initial LIGO is sensitive to, the
 * accuracy is few \f$\times 10^{-4}\f$.  </dd>
 *
 * <dt>Environment</dt><dd>
 * <ul>
 * <li>LALAPPS_DATA_PATH: Directory to look for the default mass
 * file <tt>BNSMasses.dat</tt> and the default source file <tt>inspsrcs.dat</tt>.</li>
 * </ul></dd>
 *
 * <dt>Author</dt><dd>
 * Jolien Creighton, Patrick Brady, Duncan Brown</dd>
 * </dl>
 */

#include "config.h"

#include <ctype.h>
#include <lal/Date.h>
#include <lal/LALgetopt.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOMetadataInspiralUtils.h>
#include <lal/Random.h>
#include <lal/AVFactories.h>
#include <lal/InspiralInjectionParams.h>
#include <lal/LALDetectors.h>
#include <lal/LALSimulation.h>
#include <lal/RingUtils.h>
#include <LALAppsVCSInfo.h>
#include <lal/LALDatatypes.h>
#include <lal/FrequencySeries.h>
#include "inspiral.h"
#include "LALExtTriggerTableFromLIGOLw.h"

#define CVS_REVISION "$Revision$"
#define CVS_ID_STRING "$Id$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define CVS_NAME_STRING "$Name$"
#define PROGRAM_NAME "inspinj"

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
  this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
calloc( 1, sizeof(ProcessParamsTable) ); \
snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
    PROGRAM_NAME ); \
snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
    long_options[option_index].name ); \
snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#define AXIS_MAX 12

/*
 *  *********************************
 *  Definition of the prototypes
 *  *********************************
 */
extern int vrbflg;
ProcessParamsTable *next_process_param( const char *name, const char *type,
    const char *fmt, ... );
void read_mass_data( char *filename );
void read_time_data( char *filename );
void read_nr_data( char* filename );
void read_source_data( char* filename );
void sourceComplete(void);
void drawFromSource( REAL8 *rightAscension,
    REAL8 *declination,
    REAL8 *distance,
    CHAR  name[LIGOMETA_SOURCE_MAX] );
void read_IPN_grid_from_file( char *fname );
void drawFromIPNsim( REAL8 *rightAscension,
    REAL8 *declination  );
void drawLocationFromExttrig( SimInspiralTable* table );
void drawMassFromSource( SimInspiralTable* table );
void drawMassSpinFromNR( SimInspiralTable* table );
void drawMassSpinFromNRNinja2( SimInspiralTable* table );

void adjust_snr(SimInspiralTable *inj, REAL8 target_snr, const char *ifos);
void adjust_snr_real8(SimInspiralTable *inj, REAL8 target_snr, const char *ifos);
REAL8 network_snr(const char *ifos, SimInspiralTable *inj);
REAL8 snr_in_ifo(const char *ifo, SimInspiralTable *inj);
REAL8 network_snr_real8(const char *ifos, SimInspiralTable *inj);
REAL8 snr_in_ifo_real8(const char *ifo, SimInspiralTable *inj);

REAL8 snr_in_psd_real8(const char *ifo, REAL8FrequencySeries *psd, REAL8 start_freq, SimInspiralTable *inj);
REAL8 network_snr_with_psds_real8(int num_ifos, const char **ifo_list, REAL8FrequencySeries **psds, REAL8 *start_freqs, SimInspiralTable *inj);
void adjust_snr_with_psds_real8(SimInspiralTable *inj, REAL8 target_snr, int num_ifos, const char **ifo_list, REAL8FrequencySeries **psds, REAL8 *start_freqs);

REAL8 probability_redshift(REAL8 rshift);
REAL8 luminosity_distance(REAL8 rshift);
REAL8 mean_time_step_sfr(REAL8 zmax, REAL8 rate_local);
REAL8 drawRedshift(REAL8 zmin, REAL8 zmax, REAL8 pzmax);
REAL8 redshift_mass(REAL8 mass, REAL8 z);
static void scale_lalsim_distance(SimInspiralTable *inj,char ** IFOnames, REAL8FrequencySeries **psds,REAL8 *start_freqs,LoudnessDistribution dDistr);
static REAL8 draw_uniform_snr(REAL8 snrmin,REAL8 snrmax);
static REAL8 draw_log10_snr(REAL8 snrmin,REAL8 snrmax);
static REAL8 draw_volume_snr(REAL8 snrmin,REAL8 snrmax);
/*
 *  *************************************
 *  Defining of the used global variables
 *  *************************************
 */

lalinspiral_time_distribution tDistr;
LoudnessDistribution          dDistr;
SkyLocationDistribution       lDistr;
MassDistribution              mDistr;
InclDistribution              iDistr;
SpinDistribution              spinDistr = uniformSpinDist;

SimInspiralTable *simTable;
SimRingdownTable *simRingTable;

char *massFileName = NULL;
char *nrFileName = NULL;
char *sourceFileName = NULL;
char *outputFileName = NULL;
char *injtimesFileName = NULL;
char *exttrigFileName = NULL;
char *IPNSkyPositionsFile = NULL;

INT4 outCompress = 0;
INT4 ninjaMass   = 0;
INT4 real8Ninja2 = 0;
INT4 ninjaSNR    = 0;
INT4 haveLoudness= 0;

REAL4 mwLuminosity = -1;
REAL4 minD = -1;
REAL4 maxD = -1;
REAL4 minZ = -1;
REAL4 maxZ = -1;
REAL4 localRate = -1.0;
REAL4 minSNR    = -1;
REAL4 maxSNR    = -1;
char *ifos      = NULL;

REAL4 minMass1  = -1;
REAL4 maxMass1  = -1;
REAL4 minMass2  = -1;
REAL4 maxMass2  = -1;
REAL4 minMtotal = -1;
REAL4 maxMtotal = -1;
REAL4 meanMass1 = -1.0;
REAL4 meanMass2 = -1.0;
REAL4 massStdev1= -1.0;
REAL4 massStdev2= -1.0;
REAL4 minMassRatio=-1.0;
REAL4 maxMassRatio=-1.0;

REAL4 inclStd=-1.0;
REAL4 fixed_inc=-1.0;
REAL4 max_inc=LAL_PI/2.0;
int coaPhaseFixed = 0;
REAL4 fixedCoaPhase = 0;
REAL4 psi=-1.0;
REAL4 longitude=181.0;
REAL4 latitude=91.0;
REAL4 epsAngle=1e-7;
int spinInjections=-1;
int spinAligned=-1;
REAL4 minSpin1=-1.0;
REAL4 maxSpin1=-1.0;
REAL4 meanSpin1=0.0;
REAL4 Spin1Std=0.0;
REAL4 minSpin2=-1.0;
REAL4 maxSpin2=-1.0;
REAL4 meanSpin2=0.0;
REAL4 Spin2Std=0.0;
REAL4 minKappa1=-1.0;
REAL4 maxKappa1=1.0;
REAL4 minabsKappa1=0.0;
REAL4 maxabsKappa1=1.0;
REAL4 fixedMass1=-1.0;
REAL4 fixedMass2=-1.0;
INT4  pntMass1=1;
INT4  pntMass2=1;
REAL4 deltaMass1=-1;
REAL4 deltaMass2=-1;
INT4 bandPassInj = 0;
LALSimInspiralApplyTaper taperInj = LAL_SIM_INSPIRAL_TAPER_NONE;
AlignmentType alignInj = notAligned;
REAL8 redshift;

REAL8 single_IFO_SNR_threshold=0.0;
char ** ifonames=NULL;
int numifos=0;

static LALStatus status;
static RandomParams* randParams=NULL;
INT4 numExtTriggers = 0;
ExtTriggerTable   *exttrigHead = NULL;

int num_source;
int numSkyPoints;
int galaxynum;
struct {
  char   name[LIGOMETA_SOURCE_MAX];
  REAL8 ra;
  REAL8 dec;
  REAL8 dist;
  REAL8 lum;
  REAL8 fudge;
} *source_data, *old_source_data,*temparray, *skyPoints;

char MW_name[LIGOMETA_SOURCE_MAX] = "MW";
REAL8* fracVec  =NULL;
REAL8* ratioVec = NULL;
REAL8 norm=0;

int num_mass;
struct {
  REAL8 mass1;
  REAL8 mass2;
} *mass_data;

int n_times;
LIGOTimeGPS* inj_times;

struct FakeGalaxy{
char name[LIGOMETA_SOURCE_MAX];
REAL8 ra;
REAL8 dec;
REAL8 lum;
REAL8 dist;
REAL8 fudge;
struct FakeGalaxy *next; };
int srcComplete = 0;
int makeCatalog = 0;
REAL8 srcCompleteDist;

int num_nr = 0;
int i = 0;
SimInspiralTable **nrSimArray = NULL;

/*
 *  *****************************************************
 *  Functions implementing SFR distribution over redshift
 *  *****************************************************
 */

REAL8 probability_redshift(REAL8 rshift)
{
  REAL8 pz;

  pz = -0.000429072589677+(rshift*(-0.036349728568888+(rshift*(0.860892111762314
     +(rshift*(-0.740935488674010+rshift*(0.265848831356864+rshift*(-0.050041573542298
     +rshift*(0.005184554232421+rshift*(-0.000281450045300+rshift*0.000006400690921))))))))));

  return pz;
}

REAL8 luminosity_distance(REAL8 rshift)
{
  REAL8 dL;

        dL = -2.89287707063171+(rshift*(4324.33492012756+(rshift*(3249.74193862773
           +(rshift*(-1246.66339928289+rshift*(335.354613407693+rshift*(-56.1194965448065
       +rshift*(5.20261234121263+rshift*(-0.203151569744028))))))))));

  return dL;
}

REAL8 mean_time_step_sfr(REAL8 zmax, REAL8 rate_local)
{
  REAL8 logzmax,loglambda,step;

  logzmax=log10(zmax);
  loglambda = -0.039563*pow(logzmax,6.)-0.15282*pow(logzmax,5.)-0.017596*pow(logzmax,4.)
            + 0.67193*pow(logzmax,3.)+1.1347*pow(logzmax,2.)-2.3543*logzmax+ 2.0228;
  step=pow(10.,loglambda)/rate_local;

  return step;
}

REAL8 drawRedshift(REAL8 zmin, REAL8 zmax, REAL8 pzmax)
{
  REAL8 test,z,p;
  do
  {
    test = pzmax * XLALUniformDeviate(randParams);
    z = (zmax-zmin) * XLALUniformDeviate(randParams)+zmin;
    p = probability_redshift(z);
  }
  while (test>p);

  return z;
}

REAL8 redshift_mass(REAL8 mass, REAL8 z)
{
  REAL8 mz;
  mz = mass * (1.+z);

  return mz;
}


/*************************************************************
 * Routines that calculate/adjust SNRs for REAL4 NINJA-1
 * injections.  In principle these are obsolete and could be
 * deleted.
 *************************************************************/
REAL8 snr_in_ifo(const char *ifo, SimInspiralTable *inj)
{
  REAL8 this_snr;
  REAL4TimeVectorSeries *tempStrain=NULL;

  AddNumRelStrainModes( &status, &tempStrain, inj);

  this_snr = calculate_ligo_snr_from_strain( tempStrain, inj, ifo);

  XLALDestroyREAL4VectorSequence (tempStrain->data);
  tempStrain->data = NULL;
  LALFree(tempStrain);
  tempStrain = NULL;

  return this_snr;
}

REAL8 network_snr(const char *ifo_list, SimInspiralTable *inj)
{
  char *tmp;
  char *ifo;
  REAL8 snr_total = 0.0;
  REAL8 this_snr;

  tmp = LALCalloc(1, strlen(ifos) + 1);
  strcpy(tmp, ifo_list);

  ifo = strtok (tmp,",");
  while (ifo != NULL)
  {
    this_snr   = snr_in_ifo(ifo, inj);
    snr_total += this_snr * this_snr;
    ifo        = strtok (NULL, ",");
  }

  LALFree(tmp);

  return sqrt(snr_total);
}

void adjust_snr(SimInspiralTable *inj, REAL8 target_snr, const char *ifo_list)
{
  /* Vars for calculating SNRs */
  REAL8 this_snr;
  REAL8 UNUSED low_snr, UNUSED high_snr;
  REAL8 low_dist,high_dist;

  this_snr = network_snr(ifo_list, inj);

  if (this_snr > target_snr)
  {
    high_snr  = this_snr;
    high_dist = inj->distance;

    while (this_snr > target_snr)
    {
      inj-> distance = inj->distance * 3.0;
      this_snr       = network_snr(ifo_list, inj);
    }
    low_snr  = this_snr;
    low_dist = inj->distance;
  } else {
    low_snr  = this_snr;
    low_dist = inj->distance;

    while (this_snr < target_snr)
    {
      inj->distance = (inj->distance) / 3.0;
      this_snr      = network_snr(ifo_list, inj);
    }
    high_snr  = this_snr;
    high_dist = inj->distance;
  }

  while ( fabs(target_snr - this_snr) > 1.0 )
  {
    inj->distance = (high_dist + low_dist) / 2.0;
    this_snr = network_snr(ifo_list, inj);

    if (this_snr > target_snr)
    {
      high_snr  = this_snr;
      high_dist = inj->distance;
    } else {
      low_snr  = this_snr;
      low_dist = inj->distance;
    }
  }
}


/*************************************************************
 * Routines that calculate/adjust SNRs for REAL8 NINJA-2
 * injections, using the initial LIGO/Virgo noise curves from
 * the noisemodels package.  In principle these could be replaced
 * with the next group, which will default to these noise curves
 * when alternatives are not provided.
 *************************************************************/
REAL8 snr_in_ifo_real8(const char *ifo, SimInspiralTable *inj)
{
  REAL8       this_snr;
  REAL8TimeSeries *strain = NULL;

  strain   = XLALNRInjectionStrain(ifo, inj);
  this_snr = calculate_ligo_snr_from_strain_real8(strain, ifo);

  XLALDestroyREAL8TimeSeries (strain);

  return this_snr;
}


REAL8 network_snr_real8(const char *ifo_list, SimInspiralTable *inj)
{
  char *tmp;
  char *ifo;

  REAL8 snr_total = 0.0;
  REAL8 this_snr;

  tmp = LALCalloc(1, strlen(ifos) + 1);
  strcpy(tmp, ifo_list);
  ifo = strtok (tmp,",");

  while (ifo != NULL)
  {
    this_snr   = snr_in_ifo_real8(ifo, inj);
    snr_total += this_snr * this_snr;
    ifo        = strtok (NULL, ",");
  }

  LALFree(tmp);

  return sqrt(snr_total);
}


void adjust_snr_real8(
    SimInspiralTable *inj,
    REAL8             target_snr,
    const char       *ifo_list)
{
  /* Vars for calculating SNRs */
  REAL8 this_snr;
  REAL8 UNUSED low_snr, UNUSED high_snr;
  REAL8 low_dist,high_dist;

  this_snr = network_snr_real8(ifo_list, inj);

  if (this_snr > target_snr)
  {
    high_snr  = this_snr;
    high_dist = inj->distance;

    while (this_snr > target_snr)
    {
      inj-> distance = inj->distance * 3.0;
      this_snr       = network_snr_real8(ifo_list, inj);
    }
    low_snr  = this_snr;
    low_dist = inj->distance;
  }
  else
  {
    low_snr  = this_snr;
    low_dist = inj->distance;

    while (this_snr < target_snr)
    {
      inj->distance = (inj->distance) / 3.0;
      this_snr      = network_snr_real8(ifo_list, inj);
    }
    high_snr  = this_snr;
    high_dist = inj->distance;
  }

  while ( fabs(target_snr - this_snr) > 1.0 )
  {
    inj->distance = (high_dist + low_dist) / 2.0;
    this_snr = network_snr_real8(ifo_list, inj);

    if (this_snr > target_snr)
    {
      high_snr  = this_snr;
      high_dist = inj->distance;
    }
    else
    {
      low_snr  = this_snr;
      low_dist = inj->distance;
    }
  }
}


/*************************************************************
 * Routines that calculate/adjust SNRs for REAL8 NINJA-2
 * injections, using arbitrary LIGO/Virgo noise curves given
 * in files.
 *************************************************************/
REAL8 snr_in_psd_real8(
    const char           *ifo,
    REAL8FrequencySeries *psd,
    REAL8                 start_freq,
    SimInspiralTable     *inj)
{
  REAL8            this_snr;
  REAL8TimeSeries *strain = NULL;

  strain   = XLALNRInjectionStrain(ifo, inj);
  this_snr = calculate_snr_from_strain_and_psd_real8(strain, psd, start_freq, ifo);

  XLALDestroyREAL8TimeSeries (strain);

  return this_snr;
}

REAL8 network_snr_with_psds_real8(
    int                    num_ifos,
    const char           **ifo_list,
    REAL8FrequencySeries **psds,
    REAL8                 *start_freqs,
    SimInspiralTable      *inj)
{
  REAL8 snr_total = 0.0;
  REAL8 this_snr;

  for (i=0; i< num_ifos; i++)
  {
    this_snr   = snr_in_psd_real8(ifo_list[i], psds[i], start_freqs[i], inj);
    snr_total += this_snr * this_snr;
  }

  return sqrt(snr_total);
}

void adjust_snr_with_psds_real8(
    SimInspiralTable      *inj,
    REAL8                  target_snr,
    int                    num_ifos,
    const char           **ifo_list,
    REAL8FrequencySeries **psds,
    REAL8                 *start_freqs)
{
  /* Vars for calculating SNRs */
  REAL8 this_snr;

  this_snr = network_snr_with_psds_real8(num_ifos, ifo_list, psds, start_freqs, inj);
  inj->distance = inj->distance * (this_snr/target_snr);
}


/*
 *
 * code to step forward in the process table
 *
 */
ProcessParamsTable *next_process_param( const char *name, const char *type,
    const char *fmt, ... )
{
  ProcessParamsTable *pp;
  va_list ap;
  pp = calloc( 1, sizeof( *pp ) );
  if ( ! pp )
  {
    perror( "next_process_param" );
    exit( 1 );
  }
  strncpy( pp->program, PROGRAM_NAME, LIGOMETA_PROGRAM_MAX );
  snprintf( pp->param, LIGOMETA_PARAM_MAX, "--%s", name );
  strncpy( pp->type, type, LIGOMETA_TYPE_MAX - 1 );
  va_start( ap, fmt );
  vsnprintf( pp->value, LIGOMETA_VALUE_MAX, fmt, ap );
  va_end( ap );
  return pp;
}

/*
 *
 * print-out of the usage
 *
 */
static void print_usage(char *program)
{
  fprintf(stderr,
      "%s [options]\n"\
      "The following options are recognized.  Options not surrounded in []\n"\
      "are required. Defaults are shown in brackets\n", program );
  fprintf(stderr,
      " [--help ]                 display this message\n"\
      " [--verbose]               print progress information\n"\
      " [--user-tag] usertag      set the usertag \n"\
      " [--output ] name          overwrite the standard file naming convention\n"\
      " [--write-compress]        write a compressed xml file\n\n");\
  fprintf(stderr,
      "Waveform details:\n"\
      " [--seed] randomSeed       seed for random number generator (default : 1)\n"\
      "  --f-lower freq           lower cut-off frequency.\n"\
      "  --waveform wfm           set waveform type to wfm\n"\
      "  --amp-order              set PN order in amplitude\n\n");
  fprintf(stderr,
      "Time distribution information:\n"\
      "  --gps-start-time start   GPS start time for injections\n"\
      "  --gps-end-time end       GPS end time for injections\n"\
      "  --ipn-gps-time IPNtime   GPS end time for IPN trigger\n"\
      "  --t-distr timeDist       set the time step distribution of injections\n"\
      "                           fixed: fixed time step\n"\
      "                           uniform: uniform distribution\n"\
      "                           exponential: exponential distribution for Poisson process\n"\
      "  [--time-step] step       space injections by average of step seconds\n"\
      "                           (suggestion : 2630 / pi seconds)\n"\
      "  [--time-interval] int    distribute injections in an interval, int s\n"\
      "                           (default : 0 seconds)\n\n");
  fprintf(stderr,
      "Source distribution information:\n"\
      "  --l-distr  locDist       set the source location distribution,\n"\
      "                           locDist must be one of:\n"\
      "                           source: use locations from source-file\n"\
      "                           exttrig: use external trigger file\n"\
      "                           random: uses random locations\n"\
      "                           fixed: set fixed location\n"\
      "                           ipn: random locations from IPN skypoints\n"\
      " [--longitude] longitude   read longitude if fixed value (degrees)\n"
      " [--latitude] latitude     read latitude if fixed value (degrees)\n"
      " [--d-distr] distDist      use a distribution over physical distance\n"\
      "                           source: take distance from galaxy source file\n"\
      "                           uniform: uniform distribution in distance\n"\
      "                           distancesquared: uniform distribution in distance^2\n"\
      "                           log10: uniform distribution in log10(d) \n"\
      "                           volume: uniform distribution in volume\n"\
      "                           sfr: distribution derived from the SFR\n"\
      " [--min-distance] DMIN     set the minimum (chirp) distance to DMIN kpc\n"\
      " [--max-distance] DMAX     set the maximum (chirp) distance to DMAX kpc\n"\
      "                           min/max distance required if d-distr not 'source'\n"\
      " [--source-file] sources   read source parameters from sources\n"\
      "                           requires enable/disable milkyway\n"\
      " [--sourcecomplete] distance \n"
      "                           complete galaxy catalog out to distance (kPc)\n"\
      " [--make-catalog]          create a text file of the completed galaxy catalog\n"\
      " [--enable-milkyway] lum   enables MW injections, set MW luminosity\n"\
      " [--disable-milkyway]      disables Milky Way injections\n"\
      " [--dchirp-distr]          use a distribution over chirp distance\n"\
      "                           (normalized to a 1.4,1.4 Msun binary)\n"\
      " [--z-distr]               use a distribution over redshift\n"\
      "                           currently only 'sfr' is supported\n"\
      " [--local-rate] rho        set the local coalescence rate for --z-distr sfr\n"\
      "                           (suggestion: 1 per Mpc^3 per Myr)\n"\
      " [--min-z]                 set the minimum redshift: at least 0.2 for sfr\n"\
      " [--max-z]                 set the maximum redshift: at most 1.0 for sfr\n"\
      " [--snr-distr]             use a distribution over expected (optimal) network SNR\n"\
      "                           uniform: uniform in SNR, log10: uniform in log10(SNR)\n"\
      "                           volume: uniform in 1/SNR^3\n"\
      "                           ( Setting max-snr == min-snr will allow you to choose a fixed SNR )\n"\
      " [--ninja-snr]             use a NINJA waveform SNR calculation (if not set, use LALSimulation)\n"\
      " [--min-snr] SMIN          set the minimum network snr\n"\
      " [--max-snr] SMAX          set the maximum network snr\n"\
      " [--min-coinc-snr] sm      Set the minimum SNR in two IFOs. Neglected if a single IFO is used\n"\
      " [--ligo-psd] filename     Ascii, tab-separated file of frequency, value pairs to use for LIGO PSD in snr computation\n"\
      " [--ligo-fake-psd] PSD     LALsimulation PSD fit to use instead of a file. Allowed values: LALLIGO, LALAdLIGO\n"\
      " [--ligo-start-freq] freq  Frequency in Hz to use for LIGO snr computation\n"\
      " [--virgo-psd] filename    Ascii, tab-separated file of frequency, value pairs to use for Virgo PSD in snr computation\n"\
      " [--virgo-fake-psd] PSD    LALsimulation PSD fit to use instead of a file. Allowed values: LALVirgo, LALAdVirgo\n"\
      " [--virgo-start-freq] freq Frequency in Hz to use for Virgo snr computation\n"\
      " [--ifos] ifos             Comma-separated list of ifos to include in network SNR\n\n"\
      "  --i-distr INCDIST        set the inclination distribution, must be either\n"\
      "                           uniform: distribute uniformly over arccos(i)\n"\
      "                           gaussian: gaussian distributed in (i)\n"\
      "                           fixed: no distribution, fixed values of (i)\n"\
      " [--polarization] psi      set the polarization angle for all injections (degrees)\n"\
      " [--incl-std]  inclStd     std dev for gaussian inclination dist\n"\
      " [--fixed-inc]  fixed_inc  value for the fixed inclination angle (in degrees) if '--i-distr fixed' is chosen.\n"\
      " [--max-inc]  max_inc      value for the maximum inclination angle (in degrees) if '--i-distr uniform' is chosen. \n"\
      " [--coa-phase-distr] cDist set the coalescence phase distribution,\n"\
      "                           cDist must be one of:\n"\
      "                           uniform: use random, uniformly distributed coalescence phase [default]\n"\
      "                           fixed: set fixed coalescence phase\n"\
      " [--fixed-coa-phase] phase set the coalescence phase (in degrees) for all injections if --coa-phase-distr=fixed\n"\
      " [--ipn-file] ipnskypoints read IPN sky points from file\n"\
      " [--exttrig-file] exttrig  XML file containing external trigger\n\n");
  fprintf(stderr,
      "Mass distribution information:\n"\
      "  --m-distr massDist       set the mass distribution of injections\n"\
      "                           must be one of:\n"\
      "                           source: using file containing list of mass pairs\n"\
      "                           nrwaves: using xml file with list of NR waveforms\n"\
      "                           (requires setting max/min total masses)\n"\
      "                           totalMass: uniform distribution in total mass\n"\
      "                           componentMass: uniform in m1 and m2\n"\
      "                           gaussian: gaussian mass distribution\n"\
      "                           log: log distribution in component mass\n"\
      "                           totalMassRatio: uniform distribution in total mass and\n"\
      "                           mass ratio m1 /m2\n"\
      "                           logTotalMassUniformMassRatio: log distribution in total mass\n"\
      "                           and uniform in mass ratio\n"\
      "                           totalMassFraction: uniform distribution in total mass and\n"\
      "                           in m1 /(m1+m2)\n"\
      "                           m1m2SquareGrid: component masses on a square grid\n"\
      "                           fixMasses: fix m1 and m2 to specific values\n"\
      " [--ninja2-mass]           use the NINJA 2 mass-selection algorithm\n"\
      " [--real8-ninja2]          when distributing by SNR for NINJA2, assume frames are REAL8\n"\
      " [--mass-file] mFile       read population mass parameters from mFile\n"\
      " [--nr-file] nrFile        read mass/spin parameters from xml nrFile\n"\
      " [--min-mass1] m1min       set the minimum component mass to m1min\n"\
      " [--max-mass1] m1max       set the maximum component mass to m1max\n"\
      " [--min-mass2] m2min       set the min component mass2 to m2min\n"\
      " [--max-mass2] m2max       set the max component mass2 to m2max\n"\
      " [--min-mtotal] minTotal   sets the minimum total mass to minTotal\n"\
      " [--max-mtotal] maxTotal   sets the maximum total mass to maxTotal\n"\
      " [--fixed-mass1] fixMass1  set mass1 to fixMass1\n"\
      " [--fixed-mass2] fixMass2  set mass2 to fixMass2\n"\
      " [--mean-mass1] m1mean     set the mean value for mass1\n"\
      " [--stdev-mass1] m1std     set the standard deviation for mass1\n"\
      " [--mean-mass2] m2mean     set the mean value for mass2\n"\
      " [--stdev-mass2] m2std     set the standard deviation for mass2\n"\
      " [--min-mratio] minr       set the minimum mass ratio\n"\
      " [--max-mratio] maxr       set the maximum mass ratio\n"\
      " [--mass1-points] m1pnt    set the number of grid points in the m1 direction if '--m-distr=m1m2SquareGrid'\n"\
      " [--mass2-points] m2pnt    set the number of grid points in the m2 direction if '--m-distr=m1m2SquareGrid'\n\n");
  fprintf(stderr,
      "Spin distribution information:\n"\
      "  --disable-spin           disables spinning injections\n"\
      "  --enable-spin            enables spinning injections\n"\
      "                           One of these is required.\n"\
      "  [--spin-gaussian]        enable gaussian spin distribution\n"\
      "  --aligned                enforces the spins to be along the direction\n"\
      "                           of orbital angular momentum. Spin z-components are the only non-vanishing (unless '--axis-choice view' convention is chosen)\n"\
      "  [--axis-choice] choice   frame axis choice: 'angmomentum' (default) or 'view' to define convention for spin aligned case\n"\
      "  [--min-spin1] spin1min   Set the minimum spin1 to spin1min (0.0)\n"\
      "  [--max-spin1] spin1max   Set the maximum spin1 to spin1max (0.0)\n"\
      "  [--mean-spin1] spin1mean Set the mean for |spin1| distribution\n"\
      "  [--stdev-spin1] spin1std Set the standard deviation for |spin1|\n"\
      "  [--min-spin2] spin2min   Set the minimum spin2 to spin2min (0.0)\n"\
      "  [--max-spin2] spin2max   Set the maximum spin2 to spin2max (0.0)\n"\
      "  [--mean-spin2] spin2mean Set the mean for |spin2| distribution\n"\
      "  [--stdev-spin2] spin2std Set the standard deviation for |spin2|\n"\
      "  [--min-kappa1] kappa1min Set the minimum cos(S1.L_N) to kappa1min (-1.0)\n"\
      "  [--max-kappa1] kappa1max Set the maximum cos(S1.L_N) to kappa1max (1.0)\n"\
      "  [--min-abskappa1] abskappa1min \n"\
      "                           Set the minimum absolute value of cos(S1.L_N)\n"\
      "                           to abskappa1min (0.0)\n"\
      "  [--max-abskappa1] abskappa1max \n"\
      "                           Set the maximum absolute value of cos(S1.L_N) \n"\
      "                           to abskappa1max (1.0)\n\n");
  fprintf(stderr,
      "Tapering the injection waveform:\n"\
      "  [--taper-injection] OPT  Taper the inspiral template using option OPT\n"\
      "                            (start|end|startend) \n"\
      "  [--band-pass-injection]  sets the tapering method of the injected waveform\n\n");
}


/*
 *
 * functions to read source masses
 *
 */

  void
read_mass_data( char* filename )
{
  char line[256];
  FILE   *fp;
  int n = 0;

  fp=fopen( filename, "r" );
  if ( ! fp )
  {
    perror( "read_mass_data" );
    fprintf( stderr,
        "Error while trying to open file %s\n",
        filename );
    exit( 1 );
  }

  /* count the number of lines in the file */
  num_mass=0;
  while ( fgets( line, sizeof( line ), fp ) )
    ++num_mass;

  /* alloc space for the data */
  mass_data = LALCalloc( num_mass, sizeof(*mass_data) );
  if ( !mass_data )
  {
    fprintf( stderr, "Allocation error for mass_data\n" );
    exit( 1 );
  }

  /* 'rewind' the file */
  rewind( fp );

  /* read the file finally */
  while ( fgets( line, sizeof( line ), fp ) )
  {
    sscanf( line, "%le %le", &mass_data[n].mass1, &mass_data[n].mass2 );
    n++;
  }

  /* close the file */
  fclose( fp );
}

void read_time_data( char* filename)
{
  char line[256];
  FILE   *fp;
  int n = 0;
  INT4 this_time = 0;

  fp=fopen( filename, "r" );
  if ( ! fp )
    {
      perror( "read_time_data" );
      fprintf( stderr,
	       "Error while trying to open file %s\n",
	       filename );
      exit( 1 );
    }

  /* count the number of lines in the file */
  n_times=0;
  while ( fgets( line, sizeof( line ), fp ) )
    ++n_times;

  /* alloc space for the data */
  inj_times = LALCalloc( n_times, sizeof(*inj_times) );
  if ( !inj_times )
    {
      fprintf( stderr, "Allocation error for inj_times\n" );
      exit( 1 );
    }

  /* 'rewind' the file */
  rewind( fp );

  /* read the file finally */
  while ( fgets( line, sizeof( line ), fp ) )
    {
      sscanf( line, "%d", &this_time);
      if ( this_time < 441417609 )
	{
	  fprintf( stderr, "invalid injection time %d:\n"
		   "GPS start time is prior to "
		   "Jan 01, 1994  00:00:00 UTC:\n"
		   "(%d specified)\n",
		   n, this_time );
	  exit( 1 );
	}
      inj_times[n].gpsSeconds = this_time;
      inj_times[n].gpsNanoSeconds = 0;
      // printf("%d Time: %d\t%d\n", n, inj_times[n].gpsSeconds, inj_times[n].gpsNanoSeconds);
      n++;
    }

  /* close the file */
  fclose( fp );
}

  void
read_nr_data( char* filename )
{
  SimInspiralTable  *nrSimHead = NULL;
  SimInspiralTable  *thisEvent= NULL;
  INT4               j = 0;

  nrSimHead = XLALSimInspiralTableFromLIGOLw( filename );
  if ( !nrSimHead )
  {
    fprintf( stderr, "error: unable to read sim_inspiral table from %s\n",
        filename );
    exit( 1 );
  }
  for(num_nr = 0, thisEvent = nrSimHead; thisEvent; num_nr++, thisEvent = thisEvent->next);
  if ( !num_nr )
  {
    fprintf( stderr, "error: zero events in sim_inspiral table from %s\n",
        filename );
  }

  /* allocate an array of pointers */
  nrSimArray = (SimInspiralTable ** )
    LALCalloc( num_nr, sizeof(SimInspiralTable *) );

  if ( !nrSimArray )
  {
    fprintf( stderr, "Allocation error for nr simulations\n" );
    exit( 1 );
  }

  for( j = 0, thisEvent=nrSimHead; j < num_nr;
      ++j, thisEvent = thisEvent->next )
  {
    nrSimArray[j] = thisEvent;
    if (j > 0)
    {
      nrSimArray[j-1]->next = NULL;
    }
  }
}


/*
 *
 * functions to read source distribution
 *
 */

  void
read_source_data( char* filename )
{
  char line[256];
  FILE *fp;
  int j, k;

  fp = fopen (filename, "r" );
  if ( ! fp )
  {
    perror( "read_source_data" );
    fprintf( stderr, "Could not find file %s\n", filename );
    exit( 1 );
  }

  /* count the number of entries in this file */
  num_source = 0;
  while ( fgets( line, sizeof( line ), fp ) )
    if ( line[0] == '#' )
      continue;
    else
      ++num_source;

  /* rewind the file */
  rewind( fp );

  /* allocate space */
  source_data = LALCalloc( num_source, sizeof( *source_data ) );
  if ( ! source_data )
  {
    fprintf( stderr, "Allocation error for source_data\n" );
    exit( 1 );
  }

  j = 0;
  while ( fgets( line, sizeof( line ), fp ) )
    if ( line[0] == '#' )
      continue;
    else
    {
      char ra_sgn, dec_sgn;
      REAL8 ra_h, ra_m, dec_d, dec_m;
      int c;

      c = sscanf( line, "%s %c%le:%le %c%le:%le %le %le %le",
          source_data[j].name, &ra_sgn, &ra_h, &ra_m, &dec_sgn, &dec_d, &dec_m,
          &source_data[j].dist, &source_data[j].lum, &source_data[j].fudge );
      if ( c != 10 )
      {
        fprintf( stderr, "error parsing source datafile %s\n", sourceFileName );
        exit( 1 );
      }

      /* by convention, overall sign is carried only on hours/degrees entry */
      source_data[j].ra  = ( ra_h + ra_m / 60.0 ) * LAL_PI / 12.0;
      source_data[j].dec = ( dec_d + dec_m / 60.0 ) * LAL_PI / 180.0;

      if ( ra_sgn == '-' )
        source_data[j].ra *= -1;
      if ( dec_sgn == '-' )
        source_data[j].dec *= -1;
      ++j;
    }

  /* close file */
  fclose( fp );

  /* generate ratio and fraction vectors */

  ratioVec = calloc( num_source, sizeof( REAL8 ) );
  fracVec  = calloc( num_source, sizeof( REAL8  ) );
  if ( !ratioVec || !fracVec )
  {
    fprintf( stderr, "Allocation error for ratioVec/fracVec\n" );
    exit( 1 );
  }

  /* MW luminosity might be zero */
  norm = mwLuminosity;

  /* calculate the fractions of the different sources */
  for ( k = 0; k < num_source; ++k )
    norm += ratioVec[k] = source_data[k].lum * source_data[k].fudge;
  fracVec[0] = ratioVec[0] / norm;
  for ( k = 1; k < num_source; ++k )
    fracVec[k] = fracVec[k-1] + ratioVec[k] / norm;
}

/*
 *
 * Function to read IPN sky simulations from text file given file - read(file,ra,dec)
 *
 */

  void
read_IPN_grid_from_file( char *fname )
{
  UINT4    j;                      /* counter */
  char     line[256];              /* string holders */
  FILE     *data;                  /* file object */

  /* read file */
  data = fopen(fname, "r");
  if ( ! data )
  {
    fprintf( stderr, "Could not find file %s\n", fname );
    exit( 1 );
  }

  /* find number of lines */
  numSkyPoints = 0;
  while ( fgets( line, sizeof( line ), data ) )
    ++numSkyPoints;

  /* seek to start of file again */
  fseek(data, 0, SEEK_SET);

  /* assign memory for sky points */
  skyPoints = LALCalloc(numSkyPoints, sizeof(*skyPoints));
  if ( ! skyPoints )
  {
    fprintf( stderr, "Allocation error for skyPoints\n" );
    exit( 1 );
  }

  j = 0;
  while ( fgets( line, sizeof( line ), data ) )
  {
    REAL8 ra, dec;
    int c;

    c = sscanf( line, "%le %le", &ra, &dec );
    if ( c != 2 )
    {
      fprintf( stderr, "error parsing IPN sky points datafile %s\n", IPNSkyPositionsFile );
      exit( 1 );
    }

    /* convert to radians */
    skyPoints[j].ra  = ra * ( LAL_PI / 180.0 );  /* from degrees (IPN file) to radians */
    skyPoints[j].dec = dec * ( LAL_PI / 180.0 );
    ++j;
  }

  /* close file */
  fclose( data );
}

/*
 *
 * Function to complete galaxy catalog
 *
 */

  void
sourceComplete(void)
{
  /*  Catalog Completion Constants */
  REAL8 Mbstar       = -20.45;
  /* Mbstar = magnitude at which the number of galaxies begins to fall off exponentially, */
  /* corrected for reddening (to agree with the lum density of 0.0198) */
  REAL8 phistar      = 0.0081/0.92;     /* normalization constant */
  REAL8 alpha        = -0.9;            /* determines slope at faint end of luminosity function */
  REAL8 initDistance = 0.0;             /* minimum Distance for galaxy catalog*/
  REAL8 DeltaD       = 100.0;           /* Distance step for stepping through galaxy catalog (kPc) */
  REAL8 maxDistance  = srcCompleteDist; /* Distance to which you want to correct the catalog (kPc)*/
  REAL8 M_min        = -12.0;           /* minimum blue light magnitude */
  REAL8 M_max        = -25;             /* maximum blue light magnitude */
  REAL8 edgestep     = 0.1;             /* magnitude bin size */

  /*  Vectors  */
  REAL8Vector *phibins     = NULL; /* Magnitude bins for calculating Schechter function */
  REAL8Vector *Distance    = NULL; /* Distances from initDistance to maxDistance in steps of DeltaD */
  REAL8Vector *phi         = NULL; /* Schechter magnitude function */
  REAL8Vector *phiN        = NULL; /* Number of expected galaxies in each magnitude bin */
  REAL8Vector *N           = NULL; /* Actual number of galaxies in each magnitude bin */
  REAL8Vector *pN          = NULL; /* Running tally of the fake galaxies added to the catalog */
  REAL8Vector *Corrections = NULL; /* Number of galaxies to be added in each magnitude bin */

  /* Other Variables */
  int edgenum    = (int) ceil((M_min-M_max)/edgestep); /* Number of magnitude bins */
  const char *galaxyname = "Fake";       /* Beginning of name for all added (non-real) galaxies */
  int distnum    = (maxDistance-initDistance)/DeltaD;  /* Number of elements in Distance vector */
  int k_at_25Mpc = floor((25000-initDistance)/DeltaD); /* Initial index for Distance vector - no galaxies added before 25Mpc */
  int j,k,q;            /* Indices for loops */
  REAL8 mag;            /* Converted blue light luminosity of each galaxy */
  int mag_index;        /* Index of each galaxy when binning by magnitude */
  FILE *fp;             /* File for output of corrected galaxy catalog */
  REAL8 pow1     = 0.0; /* Used to calculate Schechter function */
  REAL8 pow2     = 0.0; /* Used to calculate Schechter function */

  REAL8 UNUSED shellLum = 0.0;

  /* Parameters for generating random sky positions */
  SimInspiralTable     *randPositionTable;
  static RandomParams*  randPositions = NULL;
  int rand_skylocation_seed = 3456;

  /* Set up linked list for added galaxies*/
  struct FakeGalaxy *myFakeGalaxy;
  struct FakeGalaxy *head; /*=myFakeGalaxy;*/
  struct FakeGalaxy *saved_next;

  /* Create the Vectors */
  phibins     = XLALCreateREAL8Vector(edgenum);
  Distance    = XLALCreateREAL8Vector(distnum+1);
  phi         = XLALCreateREAL8Vector(edgenum);
  phiN        = XLALCreateREAL8Vector(edgenum);
  N           = XLALCreateREAL8Vector(edgenum);
  pN          = XLALCreateREAL8Vector(edgenum);
  Corrections = XLALCreateREAL8Vector(edgenum);

  /* Initialize sky location parameters and FakeGalaxy linked list */
  randPositionTable = calloc(1, sizeof(SimInspiralTable));
  LALCreateRandomParams( &status, &randPositions, rand_skylocation_seed);
  galaxynum = 0;
  myFakeGalaxy = (struct FakeGalaxy*) calloc(1, sizeof(struct FakeGalaxy));
  head = myFakeGalaxy;

  /* Initialize the vectors */
  for (j=0; j<edgenum; j++)
  {
    phibins->data[j] = M_max+j*edgestep;
    phiN->data[j] = 0;
    N->data[j] = 0;
    pN->data[j] = 0;
    Corrections->data[j] = 0;

    /* Calculate the theoretical blue light magnitude in each magnitude bin */
    pow1 = -1*pow(10, (-0.4*(phibins->data[j]-Mbstar)));
    pow2 = pow(10, (-0.4*(phibins->data[j]-Mbstar)));
    phi->data[j] = 0.92*phistar*exp(pow1)*pow(pow2, alpha+1);
  }

  /* Initialize the Distance array */
  for (j=0; j<=distnum; j++)
  {
    Distance->data[j] = initDistance+j*DeltaD;
  }

  /* Iterate through Distance vector and bin galaxies according to magnitude at each distance */
  for (k = k_at_25Mpc; k<distnum; k++)
  {

    /* Reset N to zero before you count the galaxies with distances less than the current Distance */
    for (q = 0; q<edgenum;q++)
    {
      N->data[q]=0;
    }

    /* Count the number of galaxies in the spherical volume with radius Distance->data[k+1] and bin them in magnitude */
    for( q = 0; q<num_source; q++)
    {
      if ( (source_data[q].dist<=Distance->data[k+1]) )
      {
        /* Convert galaxy luminosity to blue light magnitude */
        mag = -2.5*(log10(source_data[q].lum)+7.808);
        /* Calculate which magnitude bin it falls in */
        mag_index = (int) floor((mag-M_max)/edgestep);
        /* Create a histogram array of the number of galaxies in each magnitude bin */
        if (mag_index >= 0 && mag_index<edgenum)
        {
          N->data[mag_index] += 1.0;
        }
        else printf("WARNING GALAXY DOESNT FIT IN BIN\n");
      }
    }

    /* Add galaxies to the catalog based on the difference between the expected number of galaxies and the number of galaxies in the catalog */
    for (j = 0; j<edgenum; j++)
    {
      /* Number of galaxies expected in the spherical volume with radius Distance->data[k+1] */
      phiN->data[j] =edgestep*phi->data[j]*(4.0/3.0)*LAL_PI*(pow(Distance->data[k+1]/1000.0,3));
      /* Difference between the counted number of galaxies and the expected number of galaxies */
      Corrections->data[j] = phiN->data[j] - N->data[j] - pN->data[j];
      /* If there are galaxies missing, add them */
      if (Corrections->data[j]>0.0)
      {
        for (q=0;q<floor(Corrections->data[j]);q++)
        {
          randPositionTable = XLALRandomInspiralSkyLocation( randPositionTable, randPositions);
          myFakeGalaxy->dist = Distance->data[k+1];
          myFakeGalaxy->ra = randPositionTable->longitude;
          myFakeGalaxy->dec = randPositionTable->latitude;
          myFakeGalaxy->fudge = 1;
          sprintf(myFakeGalaxy->name, "%s%d", galaxyname, galaxynum);
          myFakeGalaxy->lum = pow(10.0, (phibins->data[j]/(-2.5)-7.808));
          myFakeGalaxy->next = (struct FakeGalaxy*) calloc(1,sizeof(struct FakeGalaxy));
          myFakeGalaxy = myFakeGalaxy->next;
          galaxynum++;
          pN->data[j] += 1.0;
        }
      }
    }
  }

  /* Combine source_data (original catalog) and FakeGalaxies into one array */
  temparray = calloc((num_source+galaxynum), sizeof(*source_data));
  if ( !temparray )
  {     fprintf( stderr, "Allocation error for temparray\n" );
    exit( 1 );
  }

  for (j=0;j<num_source;j++)
  {
    temparray[j].dist = source_data[j].dist;
    temparray[j].lum = source_data[j].lum;
    sprintf(temparray[j].name, "%s", source_data[j].name);
    temparray[j].ra = source_data[j].ra;
    temparray[j].dec = source_data[j].dec;
    temparray[j].fudge = source_data[j].fudge;
  }
  myFakeGalaxy = head;
  for (j=num_source;j<(num_source+galaxynum);j++)
  {
    temparray[j].dist = myFakeGalaxy->dist;
    temparray[j].lum = myFakeGalaxy->lum;
    sprintf(temparray[j].name, "%s", myFakeGalaxy->name);
    temparray[j].ra = myFakeGalaxy->ra;
    temparray[j].dec = myFakeGalaxy->dec;
    temparray[j].fudge = myFakeGalaxy->fudge;
    myFakeGalaxy = myFakeGalaxy->next;
  }
  myFakeGalaxy->next = NULL;

  /* Point old_source_data at source_data */
  old_source_data = source_data;

  /* Point source_data at the new array*/
  source_data = temparray;
  shellLum = 0;

  if (makeCatalog == 1)
  {
    /* Write the corrected catalog to a file */
    fp = fopen("correctedcatalog.txt", "w+");
    for (j=0; j<(num_source+galaxynum);j++) {
      fprintf(fp, "%s %g %g %g %g %g \n", source_data[j].name, source_data[j].ra,
        source_data[j].dec, source_data[j].dist, source_data[j].lum, source_data[j].fudge );
    }
    fclose(fp);
  }

  /* Recalculate some variables from read_source_data that will have changed due to the addition of fake galaxies */
  ratioVec = (REAL8*) calloc( (num_source+galaxynum), sizeof( REAL8 ) );
  fracVec  = (REAL8*) calloc( (num_source+galaxynum), sizeof( REAL8 ) );
  if ( !ratioVec || !fracVec )
  {
    fprintf( stderr, "Allocation error for ratioVec/fracVec\n" );
    exit( 1 );
  }

  /* MW luminosity might be zero */
  norm = mwLuminosity;

  /* calculate the fractions of the different sources */
  for ( i = 0; i <(num_source+galaxynum); ++i )
    norm += ratioVec[i] = source_data[i].lum * source_data[i].fudge;
  fracVec[0] = ratioVec[0] / norm;
  for ( i = 1; i <(num_source+galaxynum); ++i )
    fracVec[i] = fracVec[i-1] + ratioVec[i] / norm;

  /* Free some stuff */
  myFakeGalaxy = head;
  for (j=0; j<galaxynum; j++) {
    saved_next = myFakeGalaxy->next;
    free(myFakeGalaxy);
    myFakeGalaxy = saved_next;
  }
  LALFree(old_source_data);
  LALFree( skyPoints );
  LALDestroyRandomParams(&status, &randPositions);

  XLALDestroyREAL8Vector(phibins);
  XLALDestroyREAL8Vector(Corrections);
  XLALDestroyREAL8Vector(Distance);
  XLALDestroyREAL8Vector(phi);
  XLALDestroyREAL8Vector(phiN);
  XLALDestroyREAL8Vector(N);
  XLALDestroyREAL8Vector(pN);
}

/*
 *
 * functions to draw masses from source distribution
 *
 */

  void
drawMassFromSource( SimInspiralTable* table )
{
  REAL4 m1, m2, eta;
  int mass_index=0;

  /* choose masses from the mass-list */
  mass_index = (int)( num_mass * XLALUniformDeviate( randParams ) );
  m1 = redshift_mass(mass_data[mass_index].mass1,redshift);
  m2 = redshift_mass(mass_data[mass_index].mass2,redshift);

  eta=m1 * m2 / ( ( m1 + m2 ) * ( m1 + m2 ) );
  table->mass1 = m1;
  table->mass2 = m2;
  table->eta = eta;
  table->mchirp = pow( eta, 0.6) * (m1 + m2);
}


/*
 *
 * functions to draw masses and spins from NR distribution
 *
 */

  void
drawMassSpinFromNR( SimInspiralTable* table )
{
  int mass_index=0;

  /* choose masses from the mass-list */
  mass_index = (int)( num_nr * XLALUniformDeviate( randParams ) );
  XLALRandomNRInjectTotalMass( table, randParams, minMtotal, maxMtotal,
      nrSimArray[mass_index]);
}


  void
drawMassSpinFromNRNinja2( SimInspiralTable* inj )
{
  /* For ninja2 we first select a mass, then find */
  /* a waveform that can be injected at that mass */

  int j,k;
  REAL8 startFreq, startFreqHz, massTotal;
  int indx,tmp,*indicies;

  /* Permute the indicies in a random order      */
  /* This lets us check each available waveform  */
  /* once and lets us know when no option works  */
  indicies = (int *) LALCalloc( num_nr, sizeof(int) );

  for ( j = 0; j < num_nr; j++ )
    indicies[j] = j;

  for ( j = 0; j < num_nr; j++ )
  {
    indx           = (int) ( (num_nr-j) * XLALUniformDeviate( randParams ) ) + j;
    tmp            = indicies[j];
    indicies[j]    = indicies[indx];
    indicies[indx] = tmp;
  }

  massTotal = (maxMtotal - minMtotal) * XLALUniformDeviate( randParams ) + minMtotal;

  for ( j = 0; j < num_nr; j++ )
  {
    k           = indicies[j];
    if (nrSimArray[k]->f_lower > 0.0000001)
      startFreq = nrSimArray[k]->f_lower;
    else
      startFreq   = start_freq_from_frame_url(nrSimArray[k]->numrel_data);
    /*startFreqHz = startFreq / (LAL_TWOPI * massTotal * LAL_MTSUN_SI);*/
		startFreqHz = startFreq / (massTotal);
    /* if this startFreqHz makes us happy, inject it */
    if (startFreqHz <= inj->f_lower)
    {
      /* This is a copy of XLALRandomNRInjectTotalMass without  */
      /* the random mass selection.  TODO: refactor that method */
      inj->eta    = nrSimArray[k]->eta;
      inj->mchirp = massTotal * pow(inj->eta, 3.0/5.0);

      /* set mass1 and mass2 */
      inj->mass1 = (massTotal / 2.0) * (1 + pow( (1 - 4 * inj->eta), 0.5) );
      inj->mass2 = (massTotal / 2.0) * (1 - pow( (1 - 4 * inj->eta), 0.5) );

      /* copy over the spin parameters */
      inj->spin1x = nrSimArray[k]->spin1x;
      inj->spin1y = nrSimArray[k]->spin1y;
      inj->spin1z = nrSimArray[k]->spin1z;
      inj->spin2x = nrSimArray[k]->spin2x;
      inj->spin2y = nrSimArray[k]->spin2y;
      inj->spin2z = nrSimArray[k]->spin2z;

      /* copy over the numrel information */
      inj->numrel_mode_min = nrSimArray[k]->numrel_mode_min;
      inj->numrel_mode_max = nrSimArray[k]->numrel_mode_max;
      snprintf(inj->numrel_data, LIGOMETA_STRING_MAX, "%s",
               nrSimArray[k]->numrel_data);

      XLALFree(indicies);
      return;
    }
  }

  /* If we hit the end of the list, oops */
  XLALFree(indicies);
  /* should throw an error here... */
  fprintf(stderr,"No waveform could be injected at MTotal=%f Msun\n", massTotal);
}

/*
 *
 * functions to draw sky location from source distribution
 *
 */

  void
drawFromSource(
    REAL8 *rightAscension,
    REAL8 *declination,
    REAL8 *distance,
    CHAR   name[LIGOMETA_SOURCE_MAX] )
{
  REAL4 u;
  int j;

  u=XLALUniformDeviate( randParams );

  /* draw from the source table */
  for ( j = 0; j < num_source; ++j )
  {
    if ( u < fracVec[j] )
    {
      /* put the parameters */
      *rightAscension = source_data[j].ra;
      *declination    = source_data[j].dec;
      *distance = source_data[j].dist/1000.0;
      memcpy( name, source_data[j].name,
          sizeof(CHAR) * LIGOMETA_SOURCE_MAX );
      return;
    }
  }

  /* now then, draw from MilkyWay
   * WARNING: This sets location AND distance */
  XLALRandomInspiralMilkywayLocation( rightAscension, declination, distance,
      randParams );
  memcpy( name, MW_name, sizeof(CHAR) * 30 );
}

/*
 *
 * function to draw IPN sky location from IPN simulation points
 *
 */

  void
drawFromIPNsim(
    REAL8 *rightAscension,
    REAL8 *declination )
{
  REAL4 u;
  INT4 j;

  u=XLALUniformDeviate( randParams );
  j=( int ) (u*numSkyPoints);

  /* draw from the IPN source table */
    if ( j < numSkyPoints )
    {
      /* put the parameters */
      *rightAscension = skyPoints[j].ra;
      *declination    = skyPoints[j].dec;
      return;
    }
}

/*
 *
 * function to draw sky location from exttrig source file
 *
 */

  void
drawLocationFromExttrig( SimInspiralTable* table )
{
  LIGOTimeGPS timeGRB;  /* real time of the GRB */
  REAL4 ra_rad, de_rad;
  REAL8 gmst1, gmst2;

  /* convert the position (stored as degree) to radians first */
  ra_rad = exttrigHead->event_ra  * LAL_PI_180;
  de_rad = exttrigHead->event_dec * LAL_PI_180;

  /* populate the time structures */
  timeGRB.gpsSeconds     = exttrigHead->start_time;
  timeGRB.gpsNanoSeconds = exttrigHead->start_time_ns;

  gmst1 = XLALGreenwichMeanSiderealTime(&timeGRB);
  gmst2 = XLALGreenwichMeanSiderealTime(&table->geocent_end_time);

  /* populate the table */
  table->longitude = ra_rad- gmst1 + gmst2;
  table->latitude  = de_rad;
}


/*
 *
 * generate all parameters (sky position and angles) for a random inspiral
 *
 */

int main( int argc, char *argv[] )
{
  LIGOTimeGPS gpsStartTime = {-1,0};
  LIGOTimeGPS gpsEndTime = {-1,0};
  LIGOTimeGPS IPNgpsTime = {-1,0};
  LIGOTimeGPS currentGpsTime;
  long gpsDuration;

  REAL8 meanTimeStep = -1;
  REAL8 timeInterval = 0;
  REAL4 fLower = -1;
  UINT4 useChirpDist = 0;
  REAL4 minMass10, maxMass10, minMass20, maxMass20, minMtotal0, maxMtotal0,
      meanMass10, meanMass20, massStdev10, massStdev20; /* masses at z=0 */
  REAL8 pzmax=0; /* maximal value of the probability distribution of the redshift */
  INT4 ncount;
  size_t ninj;
  int rand_seed = 1;

  /* waveform */
  CHAR waveform[LIGOMETA_WAVEFORM_MAX];
  CHAR axisChoiceString[]="angmomentum";
  CHAR dummy[256];
  INT4 amp_order = -1;
  /* xml output data */
  CHAR                  fname[256];
  CHAR                 *userTag = NULL;
  ProcessTable         *proctable;
  ProcessParamsTable   *procparams;
  SimInspiralTable     *injections;
  SimRingdownTable     *ringparams;
  ProcessParamsTable   *this_proc_param;
  LIGOLwXMLStream       *xmlfp;

  REAL8 drawnDistance = 0.0;
  REAL8 drawnRightAscension = 0.0;
  REAL8 drawnDeclination = 0.0;
  CHAR  drawnSourceName[LIGOMETA_SOURCE_MAX];
  REAL8 IPNgmst1 = 0.0;
  REAL8 IPNgmst2 = 0.0;

  REAL8 targetSNR;

  CHAR *ligoPsdFileName   = NULL;
  REAL8 ligoStartFreq     = -1;
  CHAR *virgoPsdFileName  = NULL;
  REAL8 virgoStartFreq    = -1;
  CHAR *ligoFakePsd=NULL;
  CHAR *virgoFakePsd=NULL;
  REAL8FrequencySeries *ligoPsd  = NULL;
  REAL8FrequencySeries *virgoPsd = NULL;
  XLAL_INIT_MEM(status);

  /* LALgetopt arguments */
  struct LALoption long_options[] =
  {
    {"help",                    no_argument,       0,                'h'},
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"source-file",             required_argument, 0,                'f'},
    {"mass-file",               required_argument, 0,                'm'},
    {"nr-file",                 required_argument, 0,                'c'},
    {"exttrig-file",            required_argument, 0,                'E'},
    {"f-lower",                 required_argument, 0,                'F'},
    {"gps-start-time",          required_argument, 0,                'a'},
    {"gps-end-time",            required_argument, 0,                'b'},
    {"ipn-gps-time",            required_argument, 0,                '"'},
    {"t-distr",                 required_argument, 0,                '('},
    {"time-step",               required_argument, 0,                't'},
    {"time-interval",           required_argument, 0,                'i'},
    {"time-file",               required_argument, 0,               1035},
    {"seed",                    required_argument, 0,                's'},
    {"waveform",                required_argument, 0,                'w'},
    {"amp-order",               required_argument, 0,                'q'},
    {"user-tag",                required_argument, 0,                'Z'},
    {"userTag",                 required_argument, 0,                'Z'},
    {"m-distr",                 required_argument, 0,                'd'},
    {"min-mass1",               required_argument, 0,                'j'},
    {"max-mass1",               required_argument, 0,                'k'},
    {"min-mass2",               required_argument, 0,                'J'},
    {"max-mass2",               required_argument, 0,                'K'},
    {"min-mtotal",              required_argument, 0,                'A'},
    {"max-mtotal",              required_argument, 0,                'L'},
    {"fixed-mass1",             required_argument, 0,                ']'},
    {"fixed-mass2",             required_argument, 0,                '['},
    {"mean-mass1",              required_argument, 0,                'n'},
    {"mean-mass2",              required_argument, 0,                'N'},
    {"ninja2-mass",             no_argument,       &ninjaMass,         1},
    {"real8-ninja2",            no_argument,       &real8Ninja2,       1},
    {"mass1-points",            required_argument, 0,                ':'},
    {"mass2-points",            required_argument, 0,                ';'},
    {"stdev-mass1",             required_argument, 0,                'o'},
    {"stdev-mass2",             required_argument, 0,                'O'},
    {"min-mratio",              required_argument, 0,                'x'},
    {"max-mratio",              required_argument, 0,                'y'},
    {"d-distr",                 required_argument, 0,                'e'},
    {"min-distance",            required_argument, 0,                'p'},
    {"max-distance",            required_argument, 0,                'r'},
    {"dchirp-distr",            required_argument, 0,                ','},
    {"z-distr",                 required_argument, 0,                '5'},
    {"min-z",                   required_argument, 0,                '6'},
    {"max-z",                   required_argument, 0,                '7'},
    {"local-rate",              required_argument, 0,                ')'},
    {"snr-distr",               required_argument, 0,                '1'},
    {"min-snr",                 required_argument, 0,                '2'},
    {"max-snr",                 required_argument, 0,                '3'},
    {"min-coinc-snr",           required_argument, 0,                1707},
    {"ifos",                    required_argument, 0,                '4'},
    {"ninja-snr",               no_argument,       &ninjaSNR,          1},
    {"ligo-psd",                required_argument, 0,                500},
    {"ligo-fake-psd",           required_argument, 0,                501},
    {"ligo-start-freq",         required_argument, 0,                502},
    {"virgo-psd",               required_argument, 0,                600},
    {"virgo-fake-psd",          required_argument, 0,                601},
    {"virgo-start-freq",        required_argument, 0,                602},
    {"l-distr",                 required_argument, 0,                'l'},
    {"longitude",               required_argument, 0,                'v'},
    {"latitude",                required_argument, 0,                'z'},
    {"i-distr",                 required_argument, 0,                'I'},
    {"incl-std",                required_argument, 0,                'B'},
    {"fixed-inc",               required_argument, 0,                'C'},
    {"max-inc",                 required_argument, 0,               1001},
    {"polarization",            required_argument, 0,                'S'},
    {"coa-phase-distr",         required_argument, 0,               1007},
    {"fixed-coa-phase",         required_argument, 0,               1008},
    {"sourcecomplete",          required_argument, 0,                'H'},
    {"make-catalog",            no_argument,       0,                '.'},
    {"enable-milkyway",         required_argument, 0,                'M'},
    {"disable-milkyway",        no_argument,       0,                'D'},
    {"min-spin1",               required_argument, 0,                'g'},
    {"min-kappa1",              required_argument, 0,                'Q'},
    {"max-kappa1",              required_argument, 0,                'R'},
    {"min-abskappa1",           required_argument, 0,                'X'},
    {"max-abskappa1",           required_argument, 0,                'Y'},
    {"max-spin1",               required_argument, 0,                'G'},
    {"min-spin2",               required_argument, 0,                'u'},
    {"max-spin2",               required_argument, 0,                'U'},
    {"output",                  required_argument, 0,                'P'},
    {"version",                 no_argument,       0,                'V'},
    {"enable-spin",             no_argument,       0,                'T'},
    {"disable-spin",            no_argument,       0,                'W'},
    {"aligned",                 no_argument,       0,                '@'},
    {"axis-choice",             required_argument, 0,               1009},
    {"write-compress",          no_argument,       &outCompress,       1},
    {"taper-injection",         required_argument, 0,                '*'},
    {"band-pass-injection",     no_argument,       0,                '}'},
    {"ipn-file",                required_argument, 0,                '^'},
    {"spin-gaussian",           no_argument,       0,                 1002},
    {"stdev-spin1",             required_argument, 0,                 1003},
    {"stdev-spin2",             required_argument, 0,                 1004},
    {"mean-spin1",              required_argument, 0,                 1005},
    {"mean-spin2",              required_argument, 0,                 1006},
    {0, 0, 0, 0}
  };
  int c;

  /* set up initial debugging values */
  lal_errhandler = LAL_ERR_EXIT;

  /* create the process and process params tables */
  proctable = (ProcessTable *) calloc( 1, sizeof(ProcessTable) );
  XLALGPSTimeNow(&(proctable->start_time));
  XLALPopulateProcessTable(proctable, PROGRAM_NAME, lalAppsVCSIdentInfo.vcsId,
      lalAppsVCSIdentInfo.vcsStatus, lalAppsVCSIdentInfo.vcsDate, 0);
  snprintf( proctable->comment, LIGOMETA_COMMENT_MAX, " " );
  this_proc_param = procparams = (ProcessParamsTable *)
    calloc( 1, sizeof(ProcessParamsTable) );

  /* clear the waveform field */
  memset( waveform, 0, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR) );

  /* parse the arguments */
  while ( 1 )
  {
    /* LALgetopt_long stores long option here */
    int option_index = 0;
    long int gpsinput;
    size_t LALoptarg_len;

    c = LALgetopt_long_only( argc, argv,
        "hf:m:a:b:t:s:w:i:M:*", long_options, &option_index );

    /* detect the end of the options */
    if ( c == - 1 )
    {
      break;
    }

    switch ( c )
    {
      case 0:
        /* if this option set a flag, do nothing else now */
        if ( long_options[option_index].flag != 0 )
        {
          break;
        }
        else
        {
          fprintf( stderr, "error parsing option %s with argument %s\n",
              long_options[option_index].name, LALoptarg );
          exit( 1 );
        }
        break;

      case 'f':
        LALoptarg_len = strlen( LALoptarg ) + 1;
        sourceFileName = calloc( 1, LALoptarg_len * sizeof(char) );
        memcpy( sourceFileName, LALoptarg, LALoptarg_len * sizeof(char) );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "string",
              "%s", LALoptarg );
        break;

      case 'm':
        LALoptarg_len = strlen( LALoptarg ) + 1;
        massFileName = calloc( 1, LALoptarg_len * sizeof(char) );
        memcpy( massFileName, LALoptarg, LALoptarg_len * sizeof(char) );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "string",
              "%s", LALoptarg );
        break;

      case 'c':
        LALoptarg_len = strlen( LALoptarg ) + 1;
        nrFileName = calloc( 1, LALoptarg_len * sizeof(char) );
        memcpy( nrFileName, LALoptarg, LALoptarg_len * sizeof(char) );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "string",
              "%s", LALoptarg );
        break;

      case 1035:
        LALoptarg_len = strlen( LALoptarg ) + 1;
        injtimesFileName = calloc( 1, LALoptarg_len * sizeof(char) );
        memcpy( injtimesFileName, LALoptarg, LALoptarg_len * sizeof(char) );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "string",
              "%s", LALoptarg );
        break;

      case 'E':
        LALoptarg_len = strlen( LALoptarg ) + 1;
        exttrigFileName = calloc( 1, LALoptarg_len * sizeof(char) );
        memcpy( exttrigFileName, LALoptarg, LALoptarg_len * sizeof(char) );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "string",
              "%s", LALoptarg );
        break;

      case 'F':
        fLower = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "float",
              "%f", fLower );
        break;

      case 'a':
        gpsinput = atol( LALoptarg );
        if ( gpsinput < 441417609 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GPS start time is prior to "
              "Jan 01, 1994  00:00:00 UTC:\n"
              "(%ld specified)\n",
              long_options[option_index].name, gpsinput );
          exit( 1 );
        }
        gpsStartTime.gpsSeconds = gpsinput;
        gpsStartTime.gpsNanoSeconds = 0;
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "int",
              "%ld", gpsinput );
        break;

      case 'b':
        gpsinput = atol( LALoptarg );
        if ( gpsinput < 441417609 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GPS start time is prior to "
              "Jan 01, 1994  00:00:00 UTC:\n"
              "(%ld specified)\n",
              long_options[option_index].name, gpsinput );
          exit( 1 );
        }
        gpsEndTime.gpsSeconds = gpsinput;
        gpsEndTime.gpsNanoSeconds = 0;
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "int",
              "%ld", gpsinput );
        break;

      case '"':
        gpsinput = atol( LALoptarg );
        if ( gpsinput < 441417609 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GPS start time is prior to "
              "Jan 01, 1994  00:00:00 UTC:\n"
              "(%ld specified)\n",
              long_options[option_index].name, gpsinput );
          exit( 1 );
        }
        IPNgpsTime.gpsSeconds = gpsinput;
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "int",
              "%ld", gpsinput );
        break;

      case 's':
        rand_seed = atoi( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "int",
              "%d", rand_seed );
        break;

      case '(':
        LALoptarg_len = strlen( LALoptarg ) + 1;
        memcpy( dummy, LALoptarg, LALoptarg_len );

        if (!strcmp(dummy, "fixed"))
        {
          tDistr=LALINSPIRAL_FIXED_TIME_DIST;
        }
        else if (!strcmp(dummy, "uniform"))
        {
          tDistr=LALINSPIRAL_UNIFORM_TIME_DIST;
        }
        else if (!strcmp(dummy, "exponential"))
        {
          tDistr=LALINSPIRAL_EXPONENTIAL_TIME_DIST;
        }
        else if (!strcmp(dummy, "file"))
        {
          tDistr=LALINSPIRAL_FILE_TIME_DIST;
        }
        else
        {
          tDistr=LALINSPIRAL_UNKNOWN_TIME_DIST;
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown time distribution: %s must be one of\n"
              "fixed, uniform or exponential\n",
              long_options[option_index].name, LALoptarg );
          exit( 1 );
        }
        break;

      case ')':
        localRate = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "float",
              "%le", localRate );
        if ( !(localRate > 0.) )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "local coalescence rate must be positive"
              "(%f specified)\n",
              long_options[option_index].name, localRate );
          exit( 1 );
        }
        break;

      case 't':
        meanTimeStep = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "float",
              "%le", meanTimeStep );
        break;

      case 'i':
        timeInterval = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "float",
              "%le", timeInterval );
        break;

      case 'w':
        snprintf( waveform, LIGOMETA_WAVEFORM_MAX, "%s", LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "string",
              "%s", LALoptarg );
        break;

      case 'q':
        amp_order = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "int",
              "%ld", amp_order );
      break;

      case 'M':
        /* set the luminosity of the Milky Way */
        mwLuminosity = atof( LALoptarg );
        if ( mwLuminosity < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "Milky Way luminosity must be positive"
              "(%f specified)\n",
              long_options[option_index].name, mwLuminosity );
          exit( 1 );
        }

        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "float",
              "%le", mwLuminosity );
        break;

      case 'D':
        /* set the luminosity of the Milky Way */
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "string",
              "" );
        mwLuminosity = 0;
        break;

      case 'Z':
        /* create storage for the usertag */
        LALoptarg_len = strlen( LALoptarg ) + 1;
        userTag = (CHAR *) calloc( LALoptarg_len, sizeof(CHAR) );
        memcpy( userTag, LALoptarg, LALoptarg_len );

        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s",
            PROGRAM_NAME );
        snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--userTag" );
        snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%s",
            LALoptarg );
        break;

      case 'd':
        LALoptarg_len = strlen( LALoptarg ) + 1;
        memcpy( dummy, LALoptarg, LALoptarg_len );
        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s",
            PROGRAM_NAME );
        snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--m-distr" );
        snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%s",
            LALoptarg );

        if (!strcmp(dummy, "source"))
        {
          mDistr=massFromSourceFile;
        }
        else if (!strcmp(dummy, "nrwaves"))
        {
          mDistr=massFromNRFile;
        }
        else if (!strcmp(dummy, "totalMass"))
        {
          mDistr=uniformTotalMass;
        }
        else if (!strcmp(dummy, "componentMass"))
        {
          mDistr=uniformComponentMass;
        }
        else if (!strcmp(dummy, "gaussian"))
        {
          mDistr=gaussianMassDist;
        }
        else if (!strcmp(dummy, "log"))
        {
          mDistr=logComponentMass;
        }
        else if (!strcmp(dummy, "totalMassRatio"))
        {
          mDistr=uniformTotalMassRatio;
        }
        else if (!strcmp(dummy, "logTotalMassUniformMassRatio"))
        {
          mDistr=logMassUniformTotalMassRatio;
        }
        else if (!strcmp(dummy, "m1m2SquareGrid"))
        {
          mDistr=m1m2SquareGrid;
        }
        else if (!strcmp(dummy, "fixMasses"))
        {
          mDistr=fixMasses;
        }
        else if (!strcmp(dummy, "totalMassFraction"))
        {
          mDistr=uniformTotalMassFraction;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown mass distribution: %s must be one of\n"
              "(source, nrwaves, totalMass, componentMass, gaussian, log,\n"
              "totalMassRatio, totalMassFraction, logTotalMassUniformMassRatio,\n"
              "m1m2SquareGrid, fixMasses)\n",
              long_options[option_index].name, LALoptarg );
          exit( 1 );
        }
        break;

      case 'j':
        minMass1 = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%le", minMass1 );
        break;

      case 'k':
        maxMass1 = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%le", maxMass1 );
        break;

      case 'J':
        minMass2 = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%le", minMass2 );
        break;

      case 'K':
        maxMass2 = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%le", maxMass2 );
        break;

      case 'A':
        minMtotal = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%le", minMtotal );
        break;

      case 'L':
        maxMtotal = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%le", maxMtotal );
        break;

      case 'n':
        meanMass1 = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%le", meanMass1 );
        break;

      case 'N':
        meanMass2 = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%le", meanMass2 );
        break;

      case 'o':
        massStdev1 = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%le", massStdev1 );
        break;

      case 'O':
        massStdev2 = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%le", massStdev2 );
        break;

      case 'x':
        minMassRatio = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%le", minMassRatio );
        break;

      case 'y':
        maxMassRatio = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%le", maxMassRatio );
        break;

      case ':':
        pntMass1 = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "int", "%d", pntMass1 );
        break;

      case ';':
        pntMass2 = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "int", "%d", pntMass2 );
        break;

      case ']':
        fixedMass1 = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%f", fixedMass1 );
        break;

      case '[':
        fixedMass2 = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%f", fixedMass2 );
        break;

      case 'e':
        LALoptarg_len = strlen( LALoptarg ) + 1;
        memcpy( dummy, LALoptarg, LALoptarg_len );
        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s",
            PROGRAM_NAME );
        snprintf( this_proc_param->param,LIGOMETA_PARAM_MAX,"--d-distr" );
        snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        snprintf( this_proc_param->value,LIGOMETA_VALUE_MAX,"%s", LALoptarg );
        haveLoudness += 1;  /* counter to check for clashing options */

        if (!strcmp(dummy, "source"))
        {
          dDistr=distFromSourceFile;
        }
        else if (!strcmp(dummy, "uniform"))
        {
          dDistr=uniformDistance;
        }
        else if (!strcmp(dummy, "distancesquared"))
        {
          dDistr=uniformDistanceSquared;
        }
        else if (!strcmp(dummy, "log10"))
        {
          dDistr=uniformLogDistance;
        }
        else if (!strcmp(dummy, "volume"))
        {
          dDistr=uniformVolume;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown distance distribution: "
              "%s, must be one of (uniform, distancesquared, volume, log10, source)\n",
              long_options[option_index].name, LALoptarg );
          exit( 1 );
        }
        break;

      case ',':
        LALoptarg_len = strlen( LALoptarg ) + 1;
        memcpy( dummy, LALoptarg, LALoptarg_len );
        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s",
            PROGRAM_NAME );
        snprintf( this_proc_param->param,LIGOMETA_PARAM_MAX,"--dchirp-distr" );
        snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        snprintf( this_proc_param->value,LIGOMETA_VALUE_MAX,"%s", LALoptarg );
        haveLoudness += 1; /* counter to check for clashing options */
        useChirpDist = 1;

        if (!strcmp(dummy, "uniform"))
        {
          dDistr=uniformDistance;
        }
        else if (!strcmp(dummy, "distancesquared"))
        {
          dDistr=uniformDistanceSquared;
        }
        else if (!strcmp(dummy, "log10"))
        {
          dDistr=uniformLogDistance;
        }
        else if (!strcmp(dummy, "volume"))
        {
          dDistr=uniformVolume;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown distribution: "
              "%s, must be one of (uniform, distancesquared, volume, log10)\n",
              long_options[option_index].name, LALoptarg );
          exit( 1 );
        }
        break;

      case 'p':
        /* minimum distance from earth */
        minD = (REAL4) atof( LALoptarg );
        if ( minD <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "minimum distance must be > 0: "
              "(%f kpc specified)\n",
              long_options[option_index].name, minD );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%e", minD );
        break;

      case 'r':
        /* max distance from earth */
        maxD = (REAL4) atof( LALoptarg );
        if ( maxD <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "maximum distance must be greater than 0: "
              "(%f kpc specified)\n",
              long_options[option_index].name, maxD );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%e", maxD );
        break;

      case '5':
        LALoptarg_len = strlen( LALoptarg ) + 1;
        memcpy( dummy, LALoptarg, LALoptarg_len );
        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s",
            PROGRAM_NAME );
        snprintf( this_proc_param->param,LIGOMETA_PARAM_MAX,"--z-distr" );
        snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        snprintf( this_proc_param->value,LIGOMETA_VALUE_MAX,"%s", LALoptarg );
        haveLoudness += 1; /* counter to check for clashing options */

        if (!strcmp(dummy, "sfr"))
        {
          dDistr = starFormationRate;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown redshift distribution: "
              "%s, must be sfr (other distributions may be implemented in future)\n",
              long_options[option_index].name, LALoptarg );
          exit( 1 );
        }
        break;

      case '6':
        minZ = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
            next_process_param( long_options[option_index].name,
            "float", "%le", minZ );
        if ( minZ < 0 )
        {
          fprintf(stderr,"invalid argument to --%s:\n"
                  "%s must not be less than 0.\n",
                  long_options[option_index].name, LALoptarg );
          exit( 1 );
        }
        break;

      case '7':
        maxZ = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
            next_process_param( long_options[option_index].name,
            "float", "%le", maxZ );
        if ( maxZ < 0 )
        {
          fprintf(stderr,"invalid argument to --%s:\n"
                  "%s must not be less than 0.\n",
                  long_options[option_index].name, LALoptarg );
          exit( 1 );
        }
        break;

      case '1':
        LALoptarg_len = strlen( LALoptarg ) + 1;
        memcpy( dummy, LALoptarg, LALoptarg_len );
        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s",
            PROGRAM_NAME );
        snprintf( this_proc_param->param,LIGOMETA_PARAM_MAX,"--snr-distr" );
        snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        snprintf( this_proc_param->value,LIGOMETA_VALUE_MAX,"%s", LALoptarg );
        haveLoudness += 1; /* counter to check for clashing options */

        if (!strcmp(dummy, "uniform"))
        {
          dDistr=uniformSnr;
        }
        else if (!strcmp(dummy, "log10"))
        {
          dDistr=uniformLogSnr;
        }
        else if (!strcmp(dummy, "volume"))
        {
          dDistr=uniformVolumeSnr;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown SNR distribution: "
              "%s, must be uniform, log10, or volume \n",
              long_options[option_index].name, LALoptarg );
          exit( 1 );
        }
        break;

      case '2':
        minSNR = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
            next_process_param( long_options[option_index].name,
            "float", "%le", minSNR );
        if ( minSNR < 2 )
        {
          fprintf(stderr,"invalid argument to --%s:\n"
                  "%s must be greater than 2\n",
                  long_options[option_index].name, LALoptarg );
          exit( 1 );
        }
        break;

      case '3':
        maxSNR = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
            next_process_param( long_options[option_index].name,
            "float", "%le", maxSNR );
        if ( maxSNR < 2 )
        {
          fprintf(stderr,"invalid argument to --%s:\n"
                  "%s must be greater than 2\n",
                  long_options[option_index].name, LALoptarg );
          exit( 1 );
        }
        break;

      case '4':
        LALoptarg_len = strlen( LALoptarg ) + 1;
        ifos       = calloc( 1, LALoptarg_len * sizeof(char) );
        memcpy( ifos, LALoptarg, LALoptarg_len * sizeof(char) );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "string",
              "%s", LALoptarg );
        break;

      case 'l':
        LALoptarg_len = strlen( LALoptarg ) + 1;
        memcpy( dummy, LALoptarg, LALoptarg_len );
        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s",
            PROGRAM_NAME );
        snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--l-distr" );
        snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%s",
            LALoptarg );

        if (!strcmp(dummy, "source"))
        {
          lDistr=locationFromSourceFile;
        }
        else if (!strcmp(dummy, "exttrig"))
        {
          lDistr=locationFromExttrigFile;
        }
        else if (!strcmp(dummy, "random"))
        {
          lDistr=uniformSkyLocation;
        }
        else if (!strcmp(dummy, "fixed"))
        {
          lDistr=fixedSkyLocation;
        }
        else if (!strcmp(dummy, "ipn"))
        {
          lDistr=locationFromIPNFile;
        }
	else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown location distribution: "
              "%s must be one of (source, random)\n",
              long_options[option_index].name, LALoptarg );
          exit( 1 );
        }

        break;

      case 'H':
        /* Turn on galaxy catalog completion function */
        srcComplete = 1;
        srcCompleteDist = (REAL8) atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "string", "%s", LALoptarg );
        break;

      case '.':
        /* Create a text file of completed catalog */
        makeCatalog = 1;
        break;

      case 'v':
        /* fixed location (longitude) */
        longitude =  atof( LALoptarg )*LAL_PI_180 ;
        if (longitude <= (  LAL_PI + epsAngle ) && \
            longitude >= ( -LAL_PI - epsAngle ))
        {
          this_proc_param = this_proc_param->next =
            next_process_param( long_options[option_index].name,
                "float", "%e", longitude );
        }
        else
        {
          fprintf(stderr,"invalid argument to --%s:\n"
                  "%s must be between -180. and 180. degrees\n",
                  long_options[option_index].name, LALoptarg );
          exit( 1 );
        }
        break;

      case 'z':
        /* fixed location (latitude) */
        latitude = (REAL4) atof( LALoptarg )*LAL_PI_180;
        if (latitude <= (  LAL_PI/2. + epsAngle ) && \
            latitude >= ( -LAL_PI/2. - epsAngle ))
        {
          this_proc_param = this_proc_param->next =
            next_process_param( long_options[option_index].name,
                "float", "%e", latitude );
        }
        else
        {
          fprintf(stderr,"invalid argument to --%s:\n"
                  "%s must be between -90. and 90. degrees\n",
                  long_options[option_index].name, LALoptarg );
          exit( 1 );
        }
        break;

      case 'I':
        LALoptarg_len = strlen( LALoptarg ) + 1;
        memcpy( dummy, LALoptarg, LALoptarg_len );
        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s",
            PROGRAM_NAME );
        snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--i-distr" );
        snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%s",
            LALoptarg );

        if (!strcmp(dummy, "uniform"))
        {
          iDistr=uniformInclDist;
        }
        else if (!strcmp(dummy, "gaussian"))
        {
          iDistr=gaussianInclDist;
        }
        else if (!strcmp(dummy, "fixed"))
        {
          iDistr=fixedInclDist;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown inclination distribution: "
              "%s must be one of (uniform, gaussian, fixed)\n",
              long_options[option_index].name, LALoptarg );
          exit( 1 );
        }
        break;

      case 'B':
        /* gaussian width for inclination */
        inclStd = (REAL4) atof( LALoptarg );
        if ( inclStd <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "inclination gaussian width must be greater than 0: "
              "(%f specified)\n",
              long_options[option_index].name, inclStd );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%e", inclStd );
        break;

      case 'C':
        /* fixed angle of inclination */
        fixed_inc = (REAL4) atof( LALoptarg )/180.*LAL_PI;
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%e", fixed_inc );
        break;

      case 1001:
        /* maximum angle of inclination */
        max_inc = (REAL4) atof( LALoptarg )/180.*LAL_PI;
        if ( (atof(LALoptarg) < 0.) || (atof(LALoptarg) >= 180.) ) {
          fprintf( stderr, "invalid argument to --%s:\n"
              "maximum inclination angle must be between 0 and 180 degrees:"
              "(%s specified)\n",
              long_options[option_index].name, LALoptarg );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%e", max_inc );
        break;

      case 1007:
        /* coalescence phase distribution */
        if ( strcmp( LALoptarg, "uniform" ) == 0)
          coaPhaseFixed = 0;
        else if ( strcmp( LALoptarg, "fixed" ) == 0)
          coaPhaseFixed = 1;
        else {
          fprintf( stderr, "invalid argument to --%s:\n"
              "must either uniform or fixed (%s specified)\n",
              long_options[option_index].name, LALoptarg );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s",
            PROGRAM_NAME );
        snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--coa-phase-distr" );
        snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%s",
            LALoptarg );
        break;

     case 1008:
        /* fixed coalescence phase */
        fixedCoaPhase = (REAL4) atof( LALoptarg );
        if ( (fixedCoaPhase < 0.) || (fixedCoaPhase >= 360.) ) {
          fprintf( stderr, "invalid argument to --%s:\n"
              "fixed coalescence phase must be between 0 and 360 degrees:"
              "(%s specified)\n",
              long_options[option_index].name, LALoptarg );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%e", fixedCoaPhase );
        fixedCoaPhase *= LAL_PI / 180.;
        break;

      case 'S':
        /* set the polarization angle */
        psi = (REAL4) atof( LALoptarg )/180.*LAL_PI;
        if ( (atof(LALoptarg) < 0.) || (atof(LALoptarg) >= 360.) ) {
          fprintf( stderr, "invalid argument to --%s:\n"
              "polarization angle must be between 0 and 360 degrees: "
              "(%s specified)\n",
              long_options[option_index].name, LALoptarg );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%e", psi );
        break;

      case 'P':
        LALoptarg_len = strlen( LALoptarg ) + 1;
        outputFileName = calloc( 1, LALoptarg_len * sizeof(char) );
        memcpy( outputFileName, LALoptarg, LALoptarg_len * sizeof(char) );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "string", "%s", LALoptarg );
        break;

      case 500:  /* LIGO psd file */
        LALoptarg_len      = strlen( LALoptarg ) + 1;
        ligoPsdFileName = calloc( 1, LALoptarg_len * sizeof(char) );
        memcpy( ligoPsdFileName, LALoptarg, LALoptarg_len * sizeof(char) );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "string", "%s", LALoptarg );
        break;

      case 501:  /* LIGO fake LALSim PSD */
        LALoptarg_len      = strlen( LALoptarg ) + 1;
        ligoFakePsd = calloc( 1, LALoptarg_len * sizeof(char) );
        memcpy( ligoFakePsd, LALoptarg, LALoptarg_len * sizeof(char) );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "string", "%s", LALoptarg );
        break;

      case 502:  /* LIGO start frequency */
        ligoStartFreq = (REAL8) atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%f", ligoStartFreq );
        break;

      case 600:  /* Virgo psd file */
        LALoptarg_len       = strlen( LALoptarg ) + 1;
        virgoPsdFileName = calloc( 1, LALoptarg_len * sizeof(char) );
        memcpy( virgoPsdFileName, LALoptarg, LALoptarg_len * sizeof(char) );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "string", "%s", LALoptarg );
        break;

      case 601:  /* Virgo fake LALSim PSD */
        LALoptarg_len      = strlen( LALoptarg ) + 1;
        virgoFakePsd = calloc( 1, LALoptarg_len * sizeof(char) );
        memcpy( virgoFakePsd, LALoptarg, LALoptarg_len * sizeof(char) );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "string", "%s", LALoptarg );
        break;

      case 602:  /* Virgo start frequency */
        virgoStartFreq = (REAL8) atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%f", virgoStartFreq );
        break;

      case 1707: /* Set min coincident SNR in two IFOs */
        single_IFO_SNR_threshold=(REAL8) atof(LALoptarg);
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%e", single_IFO_SNR_threshold );
        break;

      case 'g':
        minSpin1 = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%le", minSpin1 );
        break;

      case 'G':
        maxSpin1 = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%le", maxSpin1 );
        break;

      case 'Q':
        minKappa1 = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%le", minKappa1 );
        break;

      case 'R':
        maxKappa1 = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%le", maxKappa1 );
        break;

      case 'X':
        minabsKappa1 = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%le", minabsKappa1 );
        break;

      case 'Y':
        maxabsKappa1 = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%le", maxabsKappa1 );
        break;

      case 'u':
        minSpin2 = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%le", minSpin2 );
        break;

      case 'U':
        maxSpin2 = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%le", maxSpin2 );
        break;

      case 'V':
        /* print version information and exit */
        fprintf( stdout, "LIGO/LSC inspiral injection engine\n");
        XLALOutputVCSInfo(stderr, lalAppsVCSInfoList, 0, "%% ");
        exit( 0 );
        break;

      case 'T':
        /* enable spinning injections */
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "string",
              "" );
        spinInjections = 1;
        break;

      case 'W':
        /* disable spinning injections */
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "string",
              "" );
        spinInjections = 0;
        break;

      case '@':
        /* enforce aligned spins */
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "string",
              "" );
        spinAligned = 1;
        break;

      case 1009:
        /* frame axis choice */
        snprintf( axisChoiceString, AXIS_MAX, "%s", LALoptarg );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "string",
              "%s", LALoptarg );
        break;

      case '}':
        /* enable band-passing */
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "string",
              "" );
        bandPassInj = 1;
        break;

      case '*':
        /* Set injection tapering */
        if ( ! strcmp( "start", LALoptarg ) )
        {
            taperInj = LAL_SIM_INSPIRAL_TAPER_START;
        }
        else if ( ! strcmp( "end", LALoptarg ) )
        {
            taperInj = LAL_SIM_INSPIRAL_TAPER_END;
        }
        else if ( ! strcmp( "startend", LALoptarg ) )
        {
            taperInj = LAL_SIM_INSPIRAL_TAPER_STARTEND;
        }
        else
        {
            fprintf( stderr, "invalid argument to --%s:\n"
                    "unknown option specified: %s\n"
                    "(Must be one of start|end|startend)\n",
                    long_options[option_index].name, LALoptarg );
        }
        this_proc_param = this_proc_param->next =
                next_process_param( long_options[option_index].name,
                        "string", LALoptarg );
        break;

      case 'h':
        print_usage(argv[0]);
        exit( 0 );
        break;

      case '?':
        print_usage(argv[0]);
        exit( 1 );
        break;

      case '^':
        LALoptarg_len = strlen( LALoptarg ) + 1;
        IPNSkyPositionsFile = calloc( 1, LALoptarg_len * sizeof(char) );
        memcpy( IPNSkyPositionsFile, LALoptarg, LALoptarg_len * sizeof(char) );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "string",
              "%s", LALoptarg );
        break;

      case 1002:
        spinDistr = gaussianSpinDist;
        this_proc_param = this_proc_param->next =
        next_process_param( long_options[option_index].name, "string", "" );
        break;

      case 1003:
        Spin1Std = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
        next_process_param( long_options[option_index].name,
          "float", "%le", Spin1Std );
        break;

      case 1004:
        Spin2Std = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
        next_process_param( long_options[option_index].name,
          "float", "%le", Spin2Std );
        break;
      case 1005:
        meanSpin1 = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
        next_process_param( long_options[option_index].name,
          "float", "%le", meanSpin1 );
        break;
      case 1006:
        meanSpin2 = atof( LALoptarg );
        this_proc_param = this_proc_param->next =
        next_process_param( long_options[option_index].name,
          "float", "%le", meanSpin2 );
        break;


      default:
        fprintf( stderr, "unknown error while parsing options\n" );
        print_usage(argv[0]);
        exit( 1 );
    }
  }

  /* must set MW flag */
  if ( mwLuminosity < 0  && dDistr == distFromSourceFile )
  {
    fprintf( stderr,
        "Must specify either --enable-milkyway LUM or --disable-milkyway\n"\
        " when using --d-distr=source\n" );
    exit( 1 );
  }

  if ( (gpsStartTime.gpsSeconds==-1 || gpsEndTime.gpsSeconds==-1) && tDistr != LALINSPIRAL_FILE_TIME_DIST)
  {
    fprintf( stderr,
        "Must specify both --gps-start-time and --gps-end-time.\n");
    exit( 1 );
  }

  gpsDuration=gpsEndTime.gpsSeconds-gpsStartTime.gpsSeconds;

  if ( (dDistr == unknownLoudnessDist) || (haveLoudness != 1) )
  {
    fprintf(stderr,"Must specify exactly one distribution out of\n"\
        "--d-distr, --dchirp-distr, --z-distr or --snr-distr.\n");
    exit( 1 );
  }

  if ( lDistr == unknownLocationDist )
  {
    fprintf(stderr,"Must specify a location distribution (--l-distr).\n");
    exit( 1 );
  }

  if ( lDistr == fixedSkyLocation && longitude == 181. )
  {
    fprintf(stderr,
        "Must specify both --longitude and --latitude when using \n"\
        "--l-distr=fixed\n");
    exit( 1 );
  }

  if ( lDistr == fixedSkyLocation && latitude == 91. )
  {
    fprintf(stderr,
        "Must specify both --longitude and --latitude when using \n"\
        "--l-distr=fixed\n");
    exit( 1 );
  }

  if ( mDistr == unknownMassDist )
  {
    fprintf(stderr,"Must specify a mass distribution (--m-distr).\n");
    exit( 1 );
  }

  if ( iDistr == unknownInclDist )
  {
    fprintf(stderr,"Must specify an inclination distribution (--i-distr).\n");
    exit( 1 );
  }

  /* if using source file, check that file and MW choice selected */
  if ( dDistr==distFromSourceFile || lDistr==locationFromSourceFile )
  {
    if ( ! sourceFileName )
    {
      fprintf( stderr,
          "Must specify --source-file when using --d-distr or --l-distr source \n" );
      exit( 1 );
    }

    if ( ( dDistr == distFromSourceFile ) && ( minD>0.0 || maxD>0.0 ) )
    {
      fprintf( stderr,
        "Cannot specify --min-distance or --max-distance\n"\
            "if --d-distr=source\n");
      exit( 1 );
    }

    /* read the source distribution here */
    read_source_data( sourceFileName );

    /* complete the galaxy catalog */
    if (srcComplete == 1)
    {
    sourceComplete();
    }
  }

  /* if using IPN sky points file, check that file exists and read it */
  if ( lDistr==locationFromIPNFile )
  {
    if ( ! IPNSkyPositionsFile )
    {
      fprintf( stderr,
          "Must specify --ipn-file when using IPN sky points distribution\n" );
      exit( 1 );
    }

    /* read the source distribution here */
   read_IPN_grid_from_file( IPNSkyPositionsFile );
  }

  /* check compatibility of distance/loudness options */
  if ( ( dDistr == uniformDistance || dDistr == uniformDistanceSquared ||
      dDistr == uniformLogDistance || dDistr == uniformVolume ) &&
      ( minD<=0.0 || maxD<=0.0 ) )
  {
    fprintf( stderr,
      "Positive minimum and maximum distances must be specified\n");
    exit( 1 );
  }
  if ( dDistr == uniformDistance || dDistr == uniformDistanceSquared ||
      dDistr == uniformLogDistance || dDistr == uniformVolume ||
      dDistr == distFromSourceFile )
  {
    if ( minZ>0.0 || maxZ>0.0 || localRate>0.0 || ninjaSNR ||
        minSNR>0.0 || maxSNR>0.0 || ifos!=NULL )
    {
      fprintf( stderr,
        "One or more options on redshift or SNR are incompatible\n"\
           "with --d-distr or --dchirp-distr !\n");
      exit( 1 );
    }
  }
  if ( dDistr == starFormationRate )
  {
    if ( minD>0.0 || maxD>0.0 || ninjaSNR || minSNR>0.0 || maxSNR>0.0 || ifos!=NULL )
    {
      fprintf( stderr,
          "One or more options on distance or SNR are incompatible\n"\
              "with --z-distr !\n");
      exit( 1 );
    }
    if ( minZ<0.2 || maxZ>1.0 )
    {
      fprintf( stderr,
          "Redshift can only take values between 0.2 and 1 for --z-distr=sfr\n");
      exit( 1 );
    }
    if ( localRate<=0. )
    {
      fprintf( stderr,
          "Local coalescence rate must be positive for --z-distr=sfr\n");
      exit( 1 );
    }
    if ( meanTimeStep>=0. )
    {
      fprintf( stderr, "Time step cannot be specified for --z-distr=sfr\n"\
              "(it is calculated from local coalescence rate)\n");
      exit( 1 );
    }
  }
  if ( ( dDistr == uniformSnr || dDistr == uniformLogSnr ||
      dDistr == uniformVolumeSnr ) &&
      ( minD>0.0 || maxD>0.0 || minZ>0.0 || maxZ>0.0 || localRate>0.0 ) )
  {
    fprintf( stderr,
        "One or more options on distance or redshift are incompatible\n"\
            "with --snr-distr !\n");
    exit( 1 );
  }

  /* SNR distribution options */

  if ( dDistr == uniformSnr || dDistr == uniformLogSnr ||
      dDistr == uniformVolumeSnr )
  {
    /* basic input checks */
    if ( minSNR == -1 || maxSNR == -1 || ifos == NULL )
    {
      fprintf( stderr,
        "Must provide all of --min-snr, --max-snr and --ifos to distribute by SNR\n" );
      exit( 1 );
    }
    if (single_IFO_SNR_threshold > maxSNR)
    {
      fprintf(stderr,
        "Minimum coincident SNR should not be larger than maximum network SNR. Exiting...\n");
      exit(1);
    }

    /* Check that each ifo has its PSD file or fakePSD */
    char *tmp, *ifo;

    /* Get the number and names of IFOs */
    tmp = LALCalloc(1, strlen(ifos) + 1);
    strcpy(tmp, ifos);
    ifo = strtok(tmp,",");

    while (ifo != NULL)
    {
      numifos += 1;
      ifo = strtok(NULL, ",");
    }

    ifonames=realloc(ifonames,(numifos+2)*sizeof(CHAR **));
    strcpy(tmp, ifos);
    ifo = strtok(tmp,",");
    ifonames[0]=malloc(strlen(ifo)+1);
    sprintf(ifonames[0],"%s",ifo);

    i=1;
    while (ifo != NULL)
    {
      ifo = strtok(NULL, ",");
      if (ifo!=NULL)
      {
        ifonames[i]=malloc(strlen(ifo)+1);
        sprintf(ifonames[i],"%s",ifo);
      }
      i++;
    }

    ifonames[numifos]=NULL;
    i=0;
    while (ifonames[i]!= NULL)
    {
      if (!strcmp(ifonames[i],"H1") || !strcmp(ifonames[i],"L1"))
      {
        /* Check either PSD file or fakePSD are given */
        if(!( ligoFakePsd || ligoPsdFileName ))
        {
          fprintf( stderr,
           "Must provide PSD file or the name of analytic PSD for LIGO if --snr-distr is given and H1 or L1 are in --ifos. \n" );
          exit( 1 );
        }
        /* Check we didn't give both fake PSD and filename*/
        if ( ligoFakePsd && ligoPsdFileName )
        {
          fprintf( stderr,"Must provide only one between --ligo-psd and --ligo-fake-psd \n" );
          exit( 1 );
        }
          /* Check flow for SNR calculation was given */
        if (ligoStartFreq < 0)
        {
          fprintf( stderr, "Must specify --ligo-start-freq together with --ligo-psd or --ligo-fake-psd.\n");
          exit( 1 );
        }
            /* Ninja only work with PSD file, check user provided it */
        if (ninjaSNR && !(ligoPsdFileName))
        {
          fprintf( stderr, "Ninja injections do not support SNR calculation with simulated PSD.\n"
            "Please provide PSD file for LIGO (with --ligo-psd filename.dat). Exiting...\n");
          exit( 1 );
        }
      }
      else if (!strcmp(ifonames[i],"V1"))
      {
        /* Check either PSD file or fakePSD are given */
        if (!(virgoFakePsd || virgoPsdFileName ))
        {
          fprintf( stderr,
            "Must provide PSD file or the name of analytic PSD for Virgo if --snr-distr is given and V1 is in --ifos. \n" );
          exit( 1 );
        }
        /* Check we didn't give both fake PSD and filename*/
        if ( virgoFakePsd && virgoPsdFileName )
        {
          fprintf( stderr,"Must provide only one between --virgo-psd and --virgo-fake-psd \n" );
          exit( 1 );
        }
        /* Check flow for SNR calculation was given */
        if (virgoStartFreq < 0)
        {
          fprintf( stderr,"Must specify --virgo-start-freq with --virgo-psd or --virgo-fake-psd.\n");
          exit( 1 );
        }
        /* Ninja only work with PSD file, check user provided it */
        if (ninjaSNR && !(virgoPsdFileName))
        {
          fprintf( stderr, "Ninja injections do not support SNR calculation with simulated PSD.\n"
            "Please provide PSD file for Virgo (with --virgo-psd filename.dat). Exiting...\n");
          exit( 1 );
        }
      }
      i++;
    }

    if (tmp) LALFree(tmp);
    if (ifo) LALFree(ifo);

    if ( maxSNR < minSNR )
    {
      fprintf( stderr, "max SNR must be greater than min SNR\n");
      exit( 1 );
    }
    if (single_IFO_SNR_threshold<0.0)
    {
      fprintf( stderr,
        "The single IFO SNR threshold must be positive. Exiting...\n" );
      exit( 1 );
    }

    /* Check custom PSDs */
    if (ligoPsdFileName) {
      if (XLALPsdFromFile(&ligoPsd, ligoPsdFileName) != XLAL_SUCCESS)
      {
        fprintf(stderr, "Unable to load PSD file %s.\n", ligoPsdFileName);
        exit( 1 );
      }
      /* We're done with the filename */
      free(ligoPsdFileName);
    }

    if (virgoPsdFileName) {
      if (XLALPsdFromFile(&virgoPsd, virgoPsdFileName) != XLAL_SUCCESS)
      {
        fprintf(stderr, "Unable to load PSD file %s.\n", virgoPsdFileName);
        exit( 1 );
      }
      /* We're done with the filename */
      free(virgoPsdFileName);
    }
  }

  /* check if the source file is specified for distance but NOT for
     location */
  if ( dDistr==distFromSourceFile && lDistr!=locationFromSourceFile )
  {
    fprintf( stderr,
        "WARNING: source file specified for distance "
        "but NOT for location. This might give strange distributions\n" );
  }

  /* check if the location file is specified for location but NOT for
   * distances: GRB case */
  if ( dDistr!=distFromSourceFile && lDistr==locationFromSourceFile &&
      mwLuminosity>0.0 )
  {
    fprintf( stderr,
        "WARNING: source file specified for locations "
        "but NOT for distances, while Milky Way injections "
        "are allowed. This might give strange distributions\n" );
  }

  /* read in the data from the external trigger file */
  if ( lDistr == locationFromExttrigFile && !exttrigFileName )
  {
    fprintf( stderr,
        "If --l-distr exttrig is specified, must specify "
        "external trigger XML file using --exttrig-file.\n");
    exit( 1 );
  }
  if ( lDistr == locationFromExttrigFile && exttrigFileName )
  {
    numExtTriggers=LALExtTriggerTableFromLIGOLw( &exttrigHead, exttrigFileName,
        0, 1);
    fprintf(stderr,
              "Number of triggers read from the external trigger file: %d\n",
               numExtTriggers);

    if (numExtTriggers>1)
    {
      fprintf(stderr,
                "WARNING: Only 1 external trigger expected in the file '%s'",
                 exttrigFileName );
    }
    if (numExtTriggers==0)
    {
      fprintf(stderr,
                "ERROR: No external trigger found in file '%s'",
                 exttrigFileName );
      exit(1);
    }
  }

  /* check inclination distribution */
  if ( ( iDistr == gaussianInclDist ) && ( inclStd < 0.0 ) )
  {
    fprintf( stderr,
        "Must specify width for gaussian inclination distribution; \n"
        "use --incl-std.\n" );
    exit( 1 );
  }
  if ( ( iDistr == fixedInclDist ) && ( fixed_inc < 0. ) )
  {
    fprintf( stderr,
        "Must specify an inclination if you want it fixed; \n"
        "use --fixed-inc.\n" );
    exit( 1 );
  }

  /* require --f-lower be explicit */
  if ( fLower <= 0.0 )
  {
    fprintf( stderr, "--f-lower must be specified and non-zero\n" );
    exit( 1 );
  }

  /* check files have been specified correctly for mass distributions */
  if ( !massFileName && mDistr==massFromSourceFile )
  {
    fprintf( stderr,
        "Must specify either a file contining the masses (--mass-file) \n"
        "or choose another mass-distribution (--m-distr).\n" );
    exit( 1 );
  }
  if ( !nrFileName && mDistr==massFromNRFile )
  {
    fprintf( stderr,
        "Must specify either a file contining the masses (--nr-file) \n"
        "or choose another mass-distribution (--m-distr).\n" );
    exit( 1 );
  }

  /* read the masses from the mass file here, check for junk options */
  if ( massFileName && mDistr==massFromSourceFile )
  {
    if ( minMass1>0.0 || minMass2>0.0 || maxMass1>0.0 || maxMass2>0.0 ||
         minMtotal>0.0 || maxMtotal>0.0 || meanMass1>0.0 || meanMass2>0.0 ||
         massStdev1>0.0 || massStdev2>0.0 || minMassRatio>0.0 || maxMassRatio>0.0 ||
         pntMass1>1 || pntMass2>1 || fixedMass1>0.0 || fixedMass2>0.0 )
    {
      fprintf( stderr,
        "One or more mass distribution options are incompatible with \n"
        "using a file containing masses (--mass-file).\n" );
      exit( 1 );
    }
    else
    {
      read_mass_data( massFileName );
    }
  }

  /* NR option requires min- & max-mtotal, all other options are junk */
  if ( nrFileName && mDistr==massFromNRFile )
  {
    if ( minMtotal<=0.0 || maxMtotal<=0.0 )
    {
      fprintf( stderr,
        "Must specify positive min-mtotal and max-mtotal when using \n"
        "a file containing NR injection masses (--nr-file).\n" );
      exit( 1 );
    }
    else if ( minMass1>0.0 || minMass2>0.0 || maxMass1>0.0 || maxMass2>0.0 ||
         meanMass1>0.0 || meanMass2>0.0 ||
         massStdev1>0.0 || massStdev2>0.0 || minMassRatio>0.0 || maxMassRatio>0.0 ||
         pntMass1>1 || pntMass2>1 || fixedMass1>0.0 || fixedMass2>0.0 )
    {
      fprintf( stderr,
        "One or more mass distribution options are incompatible with \n"
        "using a file containing masses (--nr-file).\n" );
      exit( 1 );
    }
    else
    {
      read_nr_data ( nrFileName );
    }
  }

  /* inverse logic check */
  if ( ( massFileName && mDistr!=massFromSourceFile ) ||
    ( nrFileName && mDistr!=massFromNRFile ) )
  {
    fprintf( stderr,
      "Cannot specify a source or NR injection mass file for your choice \n"
      "of --m-distr.\n" );
    exit( 1 );
  }

  /* check options for component-based distributions */
  if ( mDistr==uniformTotalMass || mDistr==uniformComponentMass ||
       mDistr==logComponentMass || mDistr==gaussianMassDist ||
       mDistr==m1m2SquareGrid )
  {
    /* first check: require component mass ranges */
    if ( minMass1<=0.0 || minMass2<=0.0 || maxMass1<=0.0 || maxMass2<=0.0 )
    {
      fprintf( stderr,
        "Must specify positive minimum and maximum component masses for \n"
        "your choice of --m-distr.\n" );
      exit( 1 );
    }
    /* second check: exclude junk options */
    if ( minMassRatio>=0.0 || maxMassRatio>=0.0 ||
         fixedMass1>=0.0 || fixedMass2>=0.0 )
    {
      fprintf( stderr,
        "Cannot specify --min-mratio, --max-mratio, or fixed m1,m2 for \n"
        "your choice of --m-distr.\n" );
      exit( 1 );
    }
    /* if max-mtotal and min-mtotal values were not specified, assign them */
    if ( minMtotal<0.0 )
      { minMtotal = minMass1 + minMass2; }
    if ( maxMtotal<0.0 )
      { maxMtotal = maxMass1 + maxMass2; }
    /* third check: proper values of min-mtotal and max-mtotal */
    if ( maxMtotal<(minMass1 + minMass2) || minMtotal>(maxMass1 + maxMass2) )
    {
      fprintf( stderr,
        "Maximum (minimum) total mass must be larger (smaller) than \n"
        "minMass1+minMass2 (maxMass1+maxMass2). Check your arithmetic\n");
      exit( 1 );
    }
  }

  /* check options for mtotal/q-based distributions */
  if ( mDistr==uniformTotalMassRatio || mDistr==logMassUniformTotalMassRatio ||
       mDistr==uniformTotalMassFraction )
  {
    /* first check: required options */
    if ( minMassRatio<=0.0 || maxMassRatio<=0.0 ||
         minMtotal<=0.0 || maxMtotal<=0.0 )
    {
      fprintf( stderr,
        "Must specify positive min-mratio and max-mratio and min/max mtotal for \n"
        "your choice of --m-distr.\n" );
      exit( 1 );
    }
    /* second check: exclude junk options */
    if ( minMass1>=0.0 || minMass2>=0.0 || maxMass1>=0.0 || maxMass2>=0.0 ||
         meanMass1>=0.0 || meanMass2>=0.0 || massStdev1>=0.0 || massStdev2>=0.0 ||
         pntMass1>1 || pntMass2>1 || fixedMass1>=0.0 || fixedMass2>=0.0 )
    {
      fprintf( stderr,
        "Cannot specify options related to component masses for your choice \n"
        "of --m-distr.\n" );
      exit( 1 );
    }
  }

  /* check for gaussian mass distribution parameters */
  if ( mDistr==gaussianMassDist )
  {
    if ( meanMass1 <= 0.0 || massStdev1 <= 0.0 ||
         meanMass2 <= 0.0 || massStdev2 <= 0.0 )
    {
      fprintf( stderr,
        "Must specify positive --mean-mass1/2 and --stdev-mass1/2 if choosing \n"
        " --m-distr gaussian\n" );
      exit( 1 );
    }
    if ( minMass1==maxMass1 || minMass2==maxMass2 )
    {
      fprintf( stderr,
        "Must specify a nonzero range of mass1 and mass2 if choosing \n"
        " --m-distr gaussian\n" );
      exit( 1 );
    }
  }
  /* inverse logic check for junk options */
  if ( ( meanMass1>=0.0 || meanMass2>=0.0 || massStdev1>=0.0 || massStdev2>=0.0 )
        && ( mDistr!=gaussianMassDist ) )
  {
    fprintf( stderr,
      "Cannot specify --mean-mass1/2 or --stdev-mass1/2 unless choosing \n"
      " --m-distr gaussian\n" );
    exit( 1 );
  }

  /* checks for m1m2 mass grid options */
  if ( mDistr==m1m2SquareGrid )
  {
    if ( pntMass1<2 || pntMass2<2 )
    {
      fprintf( stderr,
        "Values of --mass1-points and --mass2-points must be specified \n"
        "and >= 2 if choosing --m-distr m1m2SquareGrid\n" );
      exit( 1 );
    }
    else
    {
      deltaMass1 = ( maxMass1 - minMass1 ) / (REAL4) ( pntMass1 -1 );
      deltaMass2 = ( maxMass2 - minMass2 ) / (REAL4) ( pntMass2 -1 );
    }
  }
  /* inverse logic check */
  if ( ( pntMass1>1 || pntMass2>1 ) && ( mDistr!=m1m2SquareGrid ) )
  {
    fprintf( stderr,
      "Cannot specify --mass1-points or mass2-points unless choosing \n"
      " --m-distr m1m2SquareGrid\n" );
    exit( 1 );
  }

  /* checks for fixed-mass injections */
  if ( mDistr==fixMasses )
  {
    if ( fixedMass1<=0.0 || fixedMass2<=0.0 )
    {
      fprintf( stderr, "--fixed-mass1 and --fixed-mass2 must be specified "
        "and positive if choosing --m-distr fixMasses\n" );
      exit( 1 );
    }
    /* exclude junk options */
    if ( minMass1>=0.0 || minMass2>=0.0 || maxMass1>=0.0 || maxMass2>=0.0 ||
         minMtotal>=0.0 || maxMtotal>=0.0 || meanMass1>=0.0 || meanMass2>=0.0 ||
         massStdev1>=0.0 || massStdev2>=0.0 || minMassRatio>=0.0 ||
         maxMassRatio>=0.0 || pntMass1>1 || pntMass2>1 )
    {
      fprintf( stderr,
        "One or more mass options are incompatible with --m-distr fixMasses, \n"
        "only --fixed-mass1 and --fixed-mass2 may be specified\n" );
      exit( 1 );
    }
  }

  /* check if waveform is specified */
  if ( !*waveform )
  {
    fprintf( stderr, "No waveform specified (--waveform).\n" );
    exit( 1 );
  }

  if ( spinInjections==-1 && mDistr != massFromNRFile )
  {
    fprintf( stderr,
        "Must specify --disable-spin or --enable-spin\n"
        "unless doing NR injections\n" );
    exit( 1 );
  }

  if ( spinInjections==0 && spinAligned==1 )
  {
    fprintf( stderr,
        "Must enable spin to obtain aligned spin injections.\n" );
    exit( 1 );
  }

  if ( spinInjections==1 && spinAligned==-1 &&
      ( !strncmp(waveform, "IMRPhenomB", 10) || !strncmp(waveform, "IMRPhenomC", 10) ) )
  {
    fprintf( stderr,
        "Spinning IMRPhenomB or -C injections must have the --aligned option.\n" );
    exit( 1 );
  }

  if ( spinInjections==1 )
  {
    if ( spinDistr == unknownSpinDist ) /* Not currently used */
    {
      fprintf(stderr,"Must specify a spin magnitude distribution (--spin-distr).\n");
      exit( 1 );
    }

    /* check that spin stddev is positive */
    if ( spinDistr==gaussianSpinDist && (Spin1Std <= 0.0 || Spin2Std <= 0.0))
    {
        fprintf( stderr,
            "Must specify positive |spin| standard deviations when using"
            " --spin-gaussian\n" );
        exit( 1 );
    }

    /* check that spins are in range 0 - 1 */
    if (minSpin1 < 0. || minSpin2 < 0. || maxSpin1 > 1. || maxSpin2 >1.)
    {
      fprintf( stderr,
          "Spins can only take values between 0 and 1.\n" );
      exit( 1 );
    }

    /* check max and mins are the correct way around */
    if (minSpin1 > maxSpin1 || minSpin2 > maxSpin2 )
    {
      fprintf( stderr,
          "Minimal spins must be less than maximal spins.\n" );
      exit( 1 );
    }

    /* check that spin means are within a reasonable range */
    if (spinDistr==gaussianSpinDist && ( minSpin1 - meanSpin1 > 2.0*Spin1Std || meanSpin1 - maxSpin1 > 2.0*Spin2Std ))
    {
      fprintf(stderr,"Mean of |spin1| distribution is way out of range.\n");
      exit( 1 );
	}
    if (spinDistr==gaussianSpinDist && ( minSpin2 - meanSpin2 > 2.0*Spin2Std || meanSpin2 - maxSpin2 > 2.0*Spin2Std ))
    {
      fprintf(stderr,"Mean of |spin2| distribution is way out of range.\n");
      exit( 1 );
	}

    /* check that selection criteria for kappa are unique */
    if ( (minKappa1 > -1.0 || maxKappa1 < 1.0) &&
        (minabsKappa1 > 0.0 || maxabsKappa1 < 1.0) )
    {
      fprintf( stderr,
          "Either the options [--min-kappa1,--max-kappa1] or\n"
          "[--min-abskappa1,--max-abskappa1] can be specified\n" );
      exit( 1 );
    }

    /* check that kappa is in range */
    if (minKappa1 < -1.0 || maxKappa1 > 1.0)
    {
      fprintf( stderr,
          "Kappa can only take values between -1 and +1\n" );
      exit( 1 );
    }

    /* check that kappa min-max are set correctly */
    if (minKappa1 > maxKappa1)
    {
      fprintf( stderr,
          "Minimal kappa must be less than maximal kappa\n" );
      exit( 1 );
    }

    /* check that abskappa is in range */
    if (minabsKappa1 < 0.0 || maxabsKappa1 > 1.0)
    {
      fprintf( stderr,
          "The absolute value of kappa can only take values between 0 and +1\n" );
      exit( 1 );
    }

    /* check that kappa min-max are set correctly */
    if (minabsKappa1 > maxabsKappa1)
    {
      fprintf( stderr,
          "Minimal kappa must be less than maximal kappa\n" );
      exit( 1 );
    }
  }

  if ( dDistr == starFormationRate )
  {
    /* recalculate mean time step from the SFR  */
    meanTimeStep = mean_time_step_sfr(maxZ, localRate);
  }

  if (meanTimeStep<=0 && tDistr != LALINSPIRAL_FILE_TIME_DIST)
  {
    fprintf( stderr,
             "Minimum time step value must be larger than zero\n" );
    exit( 1 );
  }

  if (!injtimesFileName && tDistr == LALINSPIRAL_FILE_TIME_DIST)
    {
      fprintf(stderr, "No filename for injection GPStimes is given. Use --time-file.\n");
    }

  if ( injtimesFileName && tDistr != LALINSPIRAL_FILE_TIME_DIST )
    {
      fprintf( stderr,
	       "Cannot specify an injection times file for your choice of --t-distr.\n" );
      exit( 1 );
    }

  if (timeInterval > 0. && (tDistr == LALINSPIRAL_EXPONENTIAL_TIME_DIST || tDistr == LALINSPIRAL_FILE_TIME_DIST) )
  {
    fprintf( stderr,
         "time interval must be zero\n" );
    exit( 1 );
  }

  if ( injtimesFileName && tDistr == LALINSPIRAL_FILE_TIME_DIST)
    {
      if (meanTimeStep > 0.)
	{
	  fprintf(stderr, "Minimum time step value must be larger than zero\n" );
	  exit(1);
	}
      // printf("Reading injection times from file %s\n", injtimesFileName);
      read_time_data(injtimesFileName);
    }

  if ( userTag && outCompress )
  {
    snprintf( fname, sizeof(fname), "HL-INJECTIONS_%d_%s-%d-%ld.xml.gz",
        rand_seed, userTag, gpsStartTime.gpsSeconds, gpsDuration );
  }
  else if ( userTag && !outCompress )
  {
    snprintf( fname, sizeof(fname), "HL-INJECTIONS_%d_%s-%d-%ld.xml",
        rand_seed, userTag, gpsStartTime.gpsSeconds, gpsDuration );
  }
  else if ( !userTag && outCompress )
  {
    snprintf( fname, sizeof(fname), "HL-INJECTIONS_%d-%d-%ld.xml.gz",
        rand_seed, gpsStartTime.gpsSeconds, gpsDuration );
  }
  else
  {
    snprintf( fname, sizeof(fname), "HL-INJECTIONS_%d-%d-%ld.xml",
        rand_seed, gpsStartTime.gpsSeconds, gpsDuration );
  }
  if ( outputFileName )
  {
    snprintf( fname, sizeof(fname), "%s",
        outputFileName);
  }

  /* increment the random seed by the GPS start time:*/
  rand_seed += gpsStartTime.gpsSeconds;

  /* set up the LAL random number generator */
  LALCreateRandomParams( &status, &randParams, rand_seed );

  this_proc_param = procparams;
  procparams = procparams->next;
  free( this_proc_param );

  /* create the first injection */
  simTable = injections = (SimInspiralTable *)
    calloc( 1, sizeof(SimInspiralTable) );

  simRingTable = ringparams = (SimRingdownTable *)
    calloc( 1, sizeof(SimRingdownTable) );

  /* set redshift to zero */
  redshift=0.;

  /* set mass distribution parameters to their value at z = 0 */
  minMass10 = minMass1;
  maxMass10 = maxMass1;
  minMass20 = minMass2;
  maxMass20 = maxMass2;
  minMtotal0 = minMtotal;
  maxMtotal0 = maxMtotal;
  meanMass10 = meanMass1;
  meanMass20 = meanMass2;
  massStdev10 = massStdev1;
  massStdev20 = massStdev2;

  /* calculate the maximal value of the probability distribution of the redshift */
  if (dDistr == starFormationRate)
  {
    pzmax = probability_redshift(maxZ);
  }

  /* loop over parameter generation until end time is reached */
  ninj = 0;
  ncount = 0;
  currentGpsTime = gpsStartTime;
  if (tDistr == LALINSPIRAL_FILE_TIME_DIST){
    currentGpsTime.gpsSeconds = inj_times[0].gpsSeconds;
    currentGpsTime.gpsNanoSeconds = inj_times[0].gpsNanoSeconds;
  }

  while ( 1 )
  {
    /* increase counter */
    ninj++;

    /* store time in table */
    simTable=XLALRandomInspiralTime( simTable, randParams,
        currentGpsTime, timeInterval );

    /* populate waveform and other parameters */
    memcpy( simTable->waveform, waveform,
        sizeof(CHAR) * LIGOMETA_WAVEFORM_MAX );
    simTable->f_lower = fLower;
    simTable->amp_order = amp_order;

    /* draw redshift and apply to mass parameters */
    if (dDistr==starFormationRate)
    {
      redshift = drawRedshift(minZ,maxZ,pzmax);

      minMass1 = redshift_mass(minMass10, redshift);
      maxMass1 = redshift_mass(maxMass10, redshift);
      meanMass1 = redshift_mass(meanMass10, redshift);
      massStdev1 = redshift_mass(massStdev10, redshift);
      minMass2 = redshift_mass(minMass20, redshift);
      maxMass2 = redshift_mass(maxMass20, redshift);
      meanMass2 = redshift_mass(meanMass20, redshift);
      massStdev2 = redshift_mass(massStdev20, redshift);
      minMtotal = redshift_mass(minMtotal0, redshift);
      maxMtotal = redshift_mass(maxMtotal0, redshift);
    }

    /* populate masses */
    if ( mDistr==massFromSourceFile )
    {
      drawMassFromSource( simTable );
    }
    else if ( mDistr==massFromNRFile )
    {
      if (ninjaMass)
        drawMassSpinFromNRNinja2( simTable );
      else
        drawMassSpinFromNR( simTable );
    }
    else if ( mDistr==gaussianMassDist )
    {
      simTable=XLALGaussianInspiralMasses( simTable, randParams,
          minMass1, maxMass1,
          meanMass1, massStdev1,
          minMass2, maxMass2,
          meanMass2, massStdev2);
    }
    else if ( mDistr==uniformTotalMassRatio )
    {
      simTable=XLALRandomInspiralTotalMassRatio(simTable, randParams,
          mDistr, minMtotal, maxMtotal, minMassRatio, maxMassRatio );
    }
    else if ( mDistr==logMassUniformTotalMassRatio )
    {
      simTable=XLALRandomInspiralTotalMassRatio(simTable, randParams,
          mDistr, minMtotal, maxMtotal, minMassRatio, maxMassRatio );
    }
    else if ( mDistr==m1m2SquareGrid )
    {
      simTable=XLALm1m2SquareGridInspiralMasses( simTable, minMass1, minMass2,
          minMtotal, maxMtotal, deltaMass1, deltaMass2, pntMass1, pntMass2,
          ninj, &ncount);
    }
    else if ( mDistr==fixMasses )
    {
      simTable=XLALFixedInspiralMasses( simTable, fixedMass1, fixedMass2);
    }
    else if ( mDistr==uniformTotalMassFraction )
    {
      simTable=XLALRandomInspiralTotalMassFraction(simTable, randParams,
          mDistr, minMtotal, maxMtotal, minMassRatio, maxMassRatio );
    }
    else {
      simTable=XLALRandomInspiralMasses( simTable, randParams, mDistr,
          minMass1, maxMass1,
          minMass2, maxMass2,
          minMtotal, maxMtotal);
    }

    /* draw location and distances */
    drawFromSource( &drawnRightAscension, &drawnDeclination, &drawnDistance,
        drawnSourceName );
    drawFromIPNsim( &drawnRightAscension, &drawnDeclination );

    /* populate distances */
    if ( dDistr == distFromSourceFile )
    {
      if ( maxD > 0 )
      {
        while ( drawnDistance > maxD/1000.0 )
        {
          drawFromSource( &drawnRightAscension, &drawnDeclination,
                          &drawnDistance, drawnSourceName );
        }
      }
      simTable->distance = drawnDistance;
    }
    else if ( dDistr == starFormationRate )
    {
       /* fit of luminosity distance  between z=0-1, in Mpc for h0=0.7, omega_m=0.3, omega_v=0.7*/
       simTable->distance = luminosity_distance(redshift);
    }
    else if (dDistr== uniformDistance || dDistr== uniformDistanceSquared || dDistr== uniformLogDistance || dDistr==uniformVolume)
    {
      simTable=XLALRandomInspiralDistance(simTable, randParams,
          dDistr, minD/1000.0, maxD/1000.0);
    }
    else if (dDistr==uniformSnr || dDistr==uniformLogSnr || dDistr==uniformVolumeSnr)
    {
      /* Set distance to just any value, e.g. 100, which will be used to scale the SNR */
      simTable->distance=100.0;
    }
    /* Possible errors (i.e. no dDistr given) should have already been caught */
    /* check just to be sure in case someone adds new LoudnessDistribution    */
    else
    {
      fprintf(stderr,"Error while generating the distance of the event.\n");
      exit(1);
    }
    /* Scale by chirp mass if desired, relative to a 1.4,1.4 object */
    if ( useChirpDist )
    {
      REAL4 scaleFac;
      scaleFac = simTable->mchirp/(2.8*pow(0.25,0.6));
      simTable->distance = simTable->distance*pow(scaleFac,5./6.);
    }

    /* populate location */
    if ( lDistr == locationFromSourceFile )
    {
      simTable->longitude = drawnRightAscension;
      simTable->latitude  = drawnDeclination;
      memcpy( simTable->source, drawnSourceName,
          sizeof(CHAR) * LIGOMETA_SOURCE_MAX );
    }
    else if ( lDistr == locationFromExttrigFile )
    {
      drawLocationFromExttrig( simTable );
    }
    else if ( lDistr == fixedSkyLocation)
    {
      simTable->longitude = longitude;
      simTable->latitude = latitude;
    }
    else if ( lDistr == uniformSkyLocation )
    {
      simTable=XLALRandomInspiralSkyLocation(simTable, randParams);
    }
    else if ( lDistr == locationFromIPNFile )
    {
      IPNgmst1 = XLALGreenwichMeanSiderealTime(&IPNgpsTime);
      IPNgmst2 = XLALGreenwichMeanSiderealTime(&simTable->geocent_end_time);
      simTable->longitude = drawnRightAscension - IPNgmst1 + IPNgmst2;
      simTable->latitude  = drawnDeclination;
    }
    else
    {
      fprintf( stderr,
               "Unknown location distribution specified. Possible choices: "
               "source, exttrig, random or fixed\n" );
      exit( 1 );
    }

    /* populate polarization, inclination, and coa_phase */
    do
    {
      simTable=XLALRandomInspiralOrientation(simTable, randParams,
                                             iDistr, inclStd);
    } while ( (fabs(cos(simTable->inclination))<cos(max_inc)) );

    /* override inclination */
    if ( iDistr == fixedInclDist )
    {
      simTable->inclination = fixed_inc;
    }

    /* override polarization angle */
    if ( psi != -1.0 )
    {
      simTable->polarization = psi;
    }

    /* override coalescence phase */
    if ( coaPhaseFixed )
    {
      simTable->coa_phase = fixedCoaPhase;
    }

    /* populate spins, if required */
    if (spinInjections==1) {
      if ( spinAligned==1 ) {
	if ( !strncmp(axisChoiceString,"angmomentum", 11) )
	  alignInj = alongzAxis;
	else {
	  if ( !strncmp(axisChoiceString,"view", 4) )
	    alignInj = inxzPlane;
	  else {
	    fprintf( stderr, "Unknown axis-choice specification: 'angmomentum' and 'view' allowed, %s given.\n",axisChoiceString);
	    exit( 1 );
	  }
	}
      }
      else
	alignInj = notAligned;

      simTable = XLALRandomInspiralSpins( simTable, randParams,
					  minSpin1, maxSpin1,
					  minSpin2, maxSpin2,
					  minKappa1, maxKappa1,
					  minabsKappa1, maxabsKappa1,
					  alignInj, spinDistr,
					  meanSpin1, Spin1Std,
					  meanSpin2, Spin2Std );
    }

    /* adjust SNR to desired distribution using NINJA calculation */
    if ( ifos != NULL && ninjaSNR )
    {
      if (dDistr==uniformSnr){
        targetSNR=draw_uniform_snr(minSNR,maxSNR);
      }
      else if(dDistr==uniformLogSnr){
        targetSNR=draw_log10_snr(minSNR,maxSNR);
      }
      else if (dDistr==uniformVolumeSnr){
        targetSNR=draw_volume_snr(minSNR,maxSNR);
      }
      else{
        fprintf(stderr,"Allowed values for --snr-distr are uniform, log10 and volume. Exiting...\n");
        exit(1);
      }

      if (! real8Ninja2)
      {
        adjust_snr(simTable, targetSNR, ifos);
      }
      else
      {
        REAL8 *start_freqs;
        const char  **ifo_list;
        REAL8FrequencySeries **psds;
        int count, num_ifos = 0;
        char *tmp, *ifo;

        tmp = LALCalloc(1, strlen(ifos) + 1);
        strcpy(tmp, ifos);
        ifo = strtok (tmp,",");

        while (ifo != NULL)
        {
          num_ifos += 1;
          ifo       = strtok (NULL, ",");
        }

        start_freqs = (REAL8 *) LALCalloc(num_ifos, sizeof(REAL8));
        ifo_list    = (const char **) LALCalloc(num_ifos, sizeof(char *));
        psds        = (REAL8FrequencySeries **) LALCalloc(num_ifos, sizeof(REAL8FrequencySeries *));

        strcpy(tmp, ifos);
        ifo   = strtok (tmp,",");
        count = 0;

        while (ifo != NULL)
        {
          ifo_list[count] = ifo;

          if (ifo_list[count][0] == 'V')
          {
            start_freqs[count] = virgoStartFreq;
            psds[count]        = virgoPsd;
          }
          else
          {
            start_freqs[count] = ligoStartFreq;
            psds[count]        = ligoPsd;
          }
          count++;
          ifo = strtok (NULL, ",");
        }

        adjust_snr_with_psds_real8(simTable, targetSNR, num_ifos, ifo_list, psds, start_freqs);

        LALFree(start_freqs);
        LALFree(ifo_list);
        LALFree(psds);
        LALFree(tmp);
      }
    }

    /* adjust SNR to desired distribution using LALSimulation WF Generator */
    if (ifos!=NULL && !ninjaSNR)
    {
      char *ifo;
      REAL8 *start_freqs;
      REAL8FrequencySeries **psds;
      i=1;

      /*reset counter */
      ifo = ifonames[0];
      i = 0;
      /* Create variables for PSDs and starting frequencies */
      start_freqs = (REAL8 *) LALCalloc(numifos+1, sizeof(REAL8));
      psds        = (REAL8FrequencySeries **) LALCalloc(numifos+1, sizeof(REAL8FrequencySeries *));

      /* Hardcoded values of srate and segment length. If changed here they must also be changed in inspiralutils.c/calculate_lalsim_snr */
      REAL8 srate = 4096.0;
      /* Increase srate for EOB WFs */
      const char* WF=simTable->waveform;
      if (strstr(WF,"EOB"))
        srate = 8192.0;
      /* We may want to increase the segment length when starting at low frequencies */
      REAL8 segment = 64.0;
      size_t seglen = (size_t) segment*srate;

      /* Fill psds and start_freqs */
      /* If the user did not provide files for the PSDs, use XLALSimNoisePSD to fill in ligoPsd and virgoPsd */
      while(ifo !=NULL)
      {
        if(!strcmp("V1",ifo))
        {
          start_freqs[i]=virgoStartFreq;
          if (!virgoPsd)
          {
            virgoPsd=XLALCreateREAL8FrequencySeries("VPSD", &(simTable->geocent_end_time), 0, 1.0/segment, &lalHertzUnit, seglen/2+1);
            get_FakePsdFromString(virgoPsd,virgoFakePsd, virgoStartFreq);
          }
          if (!virgoPsd) fprintf(stderr,"Failed to produce Virgo PSD series. Exiting...\n");
          psds[i]=virgoPsd;
        }
        else if (!strcmp("L1",ifo) || !strcmp("H1",ifo))
        {
          start_freqs[i]=ligoStartFreq;
          if (!ligoPsd)
          {
            ligoPsd=XLALCreateREAL8FrequencySeries("LPSD", &(simTable->geocent_end_time), 0, 1.0/segment, &lalHertzUnit, seglen/2+1);
            get_FakePsdFromString(ligoPsd,ligoFakePsd,ligoStartFreq);
          }
          if (!ligoPsd) fprintf(stderr,"Failed to produce LIGO PSD series. Exiting...\n");
          psds[i]=ligoPsd;
        }
        else
        {
          fprintf(stderr,"Unknown IFO. Allowed IFOs are H1,L1 and V1. Exiting...\n");
          exit(-1);
        }
        i++;
        ifo=ifonames[i];
      }

      /* If exactly one detector is specified, turn the single IFO snr check off */
      if (numifos<2)
      {
        fprintf(stdout,"Warning: You are using less than 2 IFOs. Disabling the single IFO SNR threshold check...\n");
        single_IFO_SNR_threshold=0.0;
      }

      /* This function draws a proposed SNR and rescales the distance accordingly */
      scale_lalsim_distance(simTable, ifonames, psds, start_freqs, dDistr);

      /* Clean  */
      if (psds) LALFree(psds);
      if (start_freqs) LALFree(start_freqs);
    }

    /* populate the site specific information: end times and effective distances */
    LALPopulateSimInspiralSiteInfo( &status, simTable );

    /* populate the taper options */
    {
        switch (taperInj)
        {
            case LAL_SIM_INSPIRAL_TAPER_NONE:
                 snprintf( simTable->taper, LIGOMETA_INSPIRALTAPER_MAX,
                         "%s", "TAPER_NONE");
                 break;
            case LAL_SIM_INSPIRAL_TAPER_START:
                 snprintf( simTable->taper, LIGOMETA_INSPIRALTAPER_MAX,
                         "%s", "TAPER_START");
                 break;
            case LAL_SIM_INSPIRAL_TAPER_END:
                 snprintf( simTable->taper, LIGOMETA_INSPIRALTAPER_MAX,
                         "%s", "TAPER_END");
                 break;
            case LAL_SIM_INSPIRAL_TAPER_STARTEND:
                 snprintf( simTable->taper, LIGOMETA_INSPIRALTAPER_MAX,
                         "%s", "TAPER_STARTEND");
                 break;
            default: /* Never reach here */
                 fprintf( stderr, "unknown error while populating sim_inspiral taper options\n" );
                 exit(1);
        }

    }

    /* populate the bandpass options */
    simTable->bandpass = bandPassInj;


    /* increment current time, avoiding roundoff error;
       check if end of loop is reached */
    if (tDistr == LALINSPIRAL_EXPONENTIAL_TIME_DIST)
    {
      XLALGPSAdd( &currentGpsTime, -(REAL8)meanTimeStep * log( XLALUniformDeviate(randParams) ) );
    }
    else if (tDistr == LALINSPIRAL_FILE_TIME_DIST)
    {
      if (ninj >= (size_t) n_times)
	break;
      currentGpsTime.gpsSeconds = inj_times[ninj].gpsSeconds;
      currentGpsTime.gpsNanoSeconds = inj_times[ninj].gpsNanoSeconds;
    }
    else
    {
      currentGpsTime = gpsStartTime;
      XLALGPSAdd( &currentGpsTime, ninj * meanTimeStep );
    }
    if ( XLALGPSCmp( &currentGpsTime, &gpsEndTime ) >= 0 && tDistr!=LALINSPIRAL_FILE_TIME_DIST )
      break;

  /* allocate and go to next SimInspiralTable */
    simTable = simTable->next = (SimInspiralTable *)
      calloc( 1, sizeof(SimInspiralTable) );
    simRingTable = simRingTable->next = (SimRingdownTable *)
      calloc( 1, sizeof(SimRingdownTable) );

  }

  /* destroy the structure containing the random params */
  LAL_CALL(  LALDestroyRandomParams( &status, &randParams ), &status);

  /* If we read from an external trigger file, free our external trigger.
     exttrigHead is guaranteed to have no children to free. */
  if ( exttrigHead != NULL ) {
    LALFree(exttrigHead);
  }

  /* destroy the NR data */
  if ( num_nr )
  {
    for( i = 0; i < num_nr; i++ )
    {
      LALFree( nrSimArray[i] );
    }
    LALFree( nrSimArray );
  }

  xmlfp = XLALOpenLIGOLwXMLFile( fname );
  if (!xmlfp) XLAL_ERROR(XLAL_EIO);

  XLALGPSTimeNow(&(proctable->end_time));

  int retcode = XLALWriteLIGOLwXMLProcessTable(xmlfp, proctable);
  if (retcode != XLAL_SUCCESS)
  {
    XLAL_ERROR(retcode);
  }

  if ( procparams )
  {
    retcode = XLALWriteLIGOLwXMLProcessParamsTable(xmlfp, procparams);
    if (retcode != XLAL_SUCCESS)
    {
      XLAL_ERROR(retcode);
    }
  }

  XLALSimInspiralAssignIDs ( injections, 0, 0 );
  if ( injections )
  {
    retcode = XLALWriteLIGOLwXMLSimInspiralTable(xmlfp, injections);
    if ( retcode != XLAL_SUCCESS )
    {
        XLAL_ERROR(retcode);
    }
  }

  retcode = XLALCloseLIGOLwXMLFile ( xmlfp );
  if (retcode != XLAL_SUCCESS)
  {
    XLAL_ERROR(retcode);
  }


  if (source_data)
    LALFree(source_data);
  if (mass_data)
    LALFree(mass_data);
  if (skyPoints)
    LALFree(skyPoints);

  if ( ligoPsd )
    XLALDestroyREAL8FrequencySeries( ligoPsd );
  if ( virgoPsd )
    XLALDestroyREAL8FrequencySeries( virgoPsd );
  if (ifonames) LALFree(ifonames);

  LALCheckMemoryLeaks();
  return 0;
}

static void scale_lalsim_distance(SimInspiralTable *inj,
                                  char **IFOnames,
                                  REAL8FrequencySeries **psds,
                                  REAL8 *start_freqs,
                                  LoudnessDistribution snrDistr)
{
  REAL8 proposedSNR=0.0;
  REAL8 local_min=0.0;
  REAL8 net_snr=0.0;
  REAL8 *SNRs=NULL;
  UINT4 above_threshold=0;
  REAL8 ratio=1.0;
  UINT4 j=0;
  UINT4 num_ifos=0;
  /* If not already done, set distance to 100Mpc, just to have something while calculating the actual SNR */
  if (inj->distance<=0)
    inj->distance=100.0;

  if (IFOnames ==NULL)
  {
    fprintf(stderr,"scale_lalsim_distance() called with IFOnames=NULL. Exiting...\n");
    exit(1);
  }
  char * ifo=IFOnames[0];

  /* Get the number of IFOs from IFOnames*/
  while(ifo !=NULL)
  {
    num_ifos++;
    ifo=IFOnames[num_ifos];
  }

  SNRs=calloc(num_ifos+1 ,sizeof(REAL8));
  /* Calculate the single IFO and network SNR for the dummy distance of 100Mpc */
  for (j=0;j<num_ifos;j++)
  {
    SNRs[j]=calculate_lalsim_snr(inj,IFOnames[j],psds[j],start_freqs[j]);
    net_snr+=SNRs[j]*SNRs[j];
  }
  net_snr=sqrt(net_snr);

  local_min=minSNR;
  /* Draw a proposed new SNR. Check that two or more IFOs are */
  /* above threshold (if given and if num_ifos>=2)            */
  do
  {
    above_threshold=num_ifos;
    /* Generate a new SNR from given distribution */
    if (snrDistr==uniformSnr)
    {
      proposedSNR=draw_uniform_snr(local_min,maxSNR);
    }
    else if(snrDistr==uniformLogSnr)
    {
      proposedSNR=draw_log10_snr(local_min,maxSNR);
    }
    else if (snrDistr==uniformVolumeSnr)
    {
      proposedSNR=draw_volume_snr(local_min,maxSNR);
    }
    else
    {
      fprintf(stderr,"Allowed values for snr-distr are uniform, uniformLogSnr and uniformVolumeSnr. Exiting...\n");
      exit(1);
    }

    if (vrbflg) { printf("proposed SNR %lf. Proposed new dist %lf \n",
          proposedSNR,inj->distance*net_snr/proposedSNR); }

    ratio=net_snr/proposedSNR;

    /* Check that single ifo SNRs above threshold in at least two IFOs */
    for (j=0;j<num_ifos;j++)
    {
      if (SNRs[j]<single_IFO_SNR_threshold*ratio)
      above_threshold--;
    }
    /* Set the min to the proposed SNR, so that next drawing for */
    /* this event (if necessary) will give higher SNR            */
    local_min=proposedSNR;

    /* We hit the upper bound of the Network SNR. It is simply not possible to  */
    /* have >2 IFOs with single IFO above threshold without getting the network */
    /* SNR above maxSNR.                                                        */
    /* Use the last proposed value (~maxSNR) and continue                       */
    if (maxSNR-proposedSNR<0.1 && single_IFO_SNR_threshold>0.0)
    {
      fprintf(stdout,"WARNING: Could not get two or more IFOs having SNR>%.1f without\n"
        "making the network SNR larger that its maximum value %.1f. Setting SNR to %lf.\n",
        single_IFO_SNR_threshold,maxSNR,proposedSNR);

      /* set above_threshold to 3 to go out */
      above_threshold=3;
    }
    else if (vrbflg)
    {
      if (above_threshold<2 && num_ifos>=2)
        fprintf(stdout,"WARNING: Proposed SNR does not get two or more IFOs having SNR>%.1f. Re-drawing... \n",single_IFO_SNR_threshold);
    }

  } while(!(above_threshold>=2) && num_ifos>=2);

  inj->distance=inj->distance*ratio;

  /* We may want the resulting SNR printed to a file. If so, uncomment those lines
  char SnrName[200];
  sprintf(SnrName,"SNRfile.txt");
  FILE * snrout = fopen(SnrName,"a");
  fprintf(snrout,"%lf\n",proposedSNR);
  fclose(snrout);
  */

  if (SNRs) free(SNRs);

}

static REAL8 draw_volume_snr(REAL8 minsnr,REAL8 maxsnr)
{
  REAL8 proposedSNR=0.0;
  proposedSNR=1.0/(maxsnr*maxsnr*maxsnr) +
    (1.0/(minsnr*minsnr*minsnr)-1.0/(maxsnr*maxsnr*maxsnr))*XLALUniformDeviate(randParams);
  proposedSNR=1.0/cbrt(proposedSNR);
  return proposedSNR;
}

static REAL8 draw_uniform_snr(REAL8 minsnr,REAL8 maxsnr)
{
  REAL8 proposedSNR=0.0;
  proposedSNR=minsnr+ (maxsnr-minsnr)*XLALUniformDeviate(randParams);
  return proposedSNR;
}

static REAL8 draw_log10_snr(REAL8 minsnr,REAL8 maxsnr)
{
  REAL8 proposedlogSNR=0.0;
  REAL8 logminsnr=log10(minsnr);
  REAL8 logmaxsnr=log10(maxsnr);
  proposedlogSNR=logminsnr+ (logmaxsnr-logminsnr)*XLALUniformDeviate(randParams);
  return pow(10.0,proposedlogSNR);
}

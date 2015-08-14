/*
 *  Copyright (C) 2014 Michael Puerrer, John Veitch
 *  Reduced Order Model for SEOBNR
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

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <stdbool.h>
#include <alloca.h>
#include <string.h>
#include <libgen.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_spline.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/FrequencySeries.h>
#include <lal/Date.h>
#include <lal/StringInput.h>
#include <lal/Sequence.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>

#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>
#include "LALSimIMRSEOBNRROMUtilities.c"

#include <lal/LALConfig.h>
#ifdef LAL_PTHREAD_LOCK
#include <pthread.h>
#endif


/********* Input data for spline basis points **************/
#define nk_amp 200  // number of SVD-modes == number of basis functions for amplitude
#define nk_phi 250 // number of SVD-modes == number of basis functions for phase

// Check whether we actually need the full precison
static const double gA[] = {0.0001,0.00010412,0.000108409744,0.00011287622545279999,0.00011752672594145534, 0.00012236882705024332,
   0.00012741042272471334,0.00013265973214097152,0.00013812531310517954,0.00014381607600511294,0.0001497412983365236,
   0.00015591063982798835,0.00016233415818890147,0.0001690223255062842,0.00017598604531714312,0.0001832366703842094,
   0.00019078602120403885,0.00019864640527764526,0.00020683063717508426,0.00021535205942669773,0.00022422456427507767,
   0.00023346261632321087,0.00024308127611572715,0.0002530962246916951,0.00026352378914899294,0.00027438096926193147,
   0.00028568546519552304,0.00029745570636157857,0.0003097108814636756,0.000322470969779979,0.0003357567737349142,
   0.00034958995281279267,0.00036399305886867974,0.00037898957289406936,0.000394603943297305,0.000410861625761154,
   0.00042778912474251355,0.0004454140366819051,0.0004637650949931996,0.0004828722169069194,0.0005027665522434845,
   0.0005234805341959161,0.0005450479322047879,0.0005675039070116251,0.000590885067980504,0.0006152295327813008,
   0.0006405769895318904,0.0006669687615006042,0.0006944478744744291,0.0007230591269027756,0.00075284916293117,
   0.0007838665484439342,0.0008161618502398242,0.000849787718469705,0.0008847989724706569,0.0009212526901364479,
   0.0009592083009700695,0.0009987276829700365,0.001039875263508402,0.0010827181243649481,0.001127326111088784,
   0.0011737719468656418,0.0012221313510765062,0.0012724831627408582,0.0013249094690457816,0.0013794957391704678,
   0.001436330963624291,0.0014955077993256119,0.001557122720657827,0.0016212761767489296,0.0016880727552309855,
   0.0017576213527465022,0.0018300353524796581,0.00190543280900182,0.001983936640732695,0.002065674830330882,
   0.0021507806333405143,0.0022393927954341437,0.0023316557786060305,0.002427719996684599,0.0025277420605480045,
   0.002631885033442582,0.0027403186968204163,0.0028532198271294176,0.0029707724840071495,0.003093168310348244,
   0.0032206068447345917,0.003353295846737657,0.0034914516356232485,0.0036352994430109264,0.0037850737800629764,
   0.003941018819801571,0.0041033887951773965,0.004272448413538705,0.0044484732881765,0.004631750387649371,0.004822578503620526,
   0.005021268737969691,0.0052281450099740424,0.005443544584384973,0.005667818621261634,0.005901332748457613,
   0.0061444676576940666,0.006397619725191062,0.006661201657868934,0.006935643166173134,0.007221391664619468,0.00751891300120179,
   0.007828692216851304,0.008151234336185578,0.008487065190836423,0.008836732276698884,0.009200805646498878,0.009579878839134632,
   0.009974569847306979,0.010385522125016027,0.010813405636566686,0.011258917948793233,0.011722785368283514,0.012205764125456795,
   0.012708641607425615,0.01323223764165155,0.013777405832487594,0.014345034952786082,0.014936050392840869,0.015551415669025912,
   0.01619213399458978,0.01685924991516688,0.017553851011671756,0.018277069673352634,0.019030084943894764,0.019814124443583228,
   0.020630466370658858,0.021480441585130003,0.02236543577843736,0.023286891732508978,0.024246311671888347,0.025245259712770148,
   0.026285364412936277,0.02736832142674925,0.02849589626953132,0.029669927195836013,0.030892328196304455,0.0321650921179922,
   0.03349029391325348,0.03487009402247952,0.036306741896205676,0.03780257966232935,0.03936004594441732,0.04098167983732731,
   0.0426701250466252,0.04442813419854616,0.04625857332752626,0.04816442654862034,0.0501488009224235,0.052214931520427346,
   0.05436618669906895,0.056606073591070595,0.05893824382302271,0.06136649946853124,0.06389479924663473,0.06652726497559608,
   0.06926818829259064,0.07212203765024537,0.07509346560143548,0.07818731638421463,0.08140863381924426,0.08476266953259713,
   0.08825489151734013,0.09189099304785454,0.09567690196142616,0.09961879032223692,0.10372308448351308,0.10799647556423382,
   0.11244593035748025,0.11707870268820844,0.12190234523896262,0.12692472186280787,0.13215402040355556,0.13759876604418206,
   0.14326783520520237,0.1491704700156567,0.15531629338030176,0.1617153246675702,0.16837799604387407,0.1753151694808817,
   0.182538154463494,0.19005872642738997,0.19788914595619844,0.20604217876959383,0.2145311165349011,0.223369798536139,
   0.23257263423582794,0.24215462676634406,0.25213139738911744,0.26251921096154907,0.2733350024531649,0.2845964045542353,
   0.2963217764218698,0.3};

static const double gPhi[] = {0.0001,0.0001011603972084032,0.00010233878267233898,0.00010353550576343408,0.00010475092401703343,0.00010598540335537165,
   0.00010723931831773738,0.0001085130522978776,0.00010980699778889891,0.00011112155663593303,0.00011245714029684425,
   0.00011381417011126748,0.00011519307757827774,0.00011659430464300309,0.00011801830399250643,0.00011946553936127435,
   0.00012093648584666495,0.00012243163023468122,0.0001239514713364512,0.00012549652033581237,0.00012706730114841344,
   0.00012866435079276461,0.0001302882197736849,0.00013193947247861382,0.0001336186875872749,0.00013532645849519842,
   0.00013706339375163268,0.00013883011751239573,0.00014062727000824248,0.00014245550802934775,0.00014431550542653068,
   0.00014620795362987388,0.00014813356218541844,0.00015009305931064617,0.0001520871924694914,0.00015411672896765733,
   0.00015618245656904637,0.00015828518413414986,0.00016042574228128018,0.00016260498407156794,0.00016482378571868855,
   0.00016708304732432604,0.00016938369364042784,0.0001717266748593521,0.00017411296743306026,0.00017654357492255984,
   0.000179019528878859,0.0001815418897567524,0.00018411174786282027,0.0001867302243390864,0.0001893984721838501,
   0.0001921176773112781,0.00019488905965141826,0.00019771387429237596,0.00020059341266647862,0.00020352900378234071,
   0.0002065220155048353,0.0002095738558850754,0.0002126859745426121,0.00021585986410216392,0.00021909706168730754,
   0.0002223991504736801,0.00022576776130437057,0.0002292045743703129,0.00023271132095863427,0.00023628978527206272,
   0.00023994180632265591,0.00024366927990327993,0.0002474741606404437,0.00025135846413228087,0.00025532426917566865,
   0.000259373720086681,0.00026350902911879455,0.00026773247898349915,0.0002720464254782114,0.00027645330022665123,
   0.0002809556135371191,0.0002855559573844033,0.00029025700852135906,0.00029506153172652994,0.00029997238319453087,
   0.00030499251407628374,0.00031012497417658787,0.0003153729158169264,0.0003207395978718508,0.0003262283899877574,
   0.00033184277699336776,0.0003375863635117564,0.0003434628787843332,0.000349476181717788,0.0003556302661656416,
   0.00036192926645672906,0.00036837746318366013,0.0003749792892650752,0.00038173933629633125,0.00038866236120412836,
   0.0003957532932215169,0.00040301724120071894,0.00041045950128225665,0.0004180855649400099,0.0004259011274230334,
   0.00043391209661625043,0.0004421246023435188,0.0004505450061380353,0.0004591799115066175,0.00046803617471608764,
   0.0004771209161317804,0.0004864415321401282,0.0004960057076893408,0.0005058214294844088,0.0005158969998750328,
   0.0005262410514776242,0.0005368625625752512,0.0005477708733423325,0.0005589757029440297,0.000570487167563663,
   0.0005823157994151132,0.0005944725668010725,0.0006069688952822128,0.0006198166900268609,0.0006330283594156419,
   0.0006466168399807966,0.0006605956227655398,0.0006749787811949231,0.0006897810005562521,0.0007050176091942138,
   0.0007207046115335476,0.0007368587230503925,0.0007534974073224147,0.0007706389152975295,0.0007883023269315418,
   0.0008065075953564147,0.0008252755937532143,0.0008446281651171619,0.0008645881751167426,0.0008851795682645893,
   0.0009064274276349821,0.0009283580383814316,0.0009509989553280604,0.0009743790749305488,0.0009985287119264209,
   0.001023479681020617,0.0010492653839808423,0.0010759209025483307,0.0011034830976036844,0.0011319907150646252,
   0.001161484499033161,0.0011920073127541697,0.001223604267996143,0.0012563228635182632,0.001290213133346581,
   0.0013253278056463969,0.0013617224730486212,0.0013994557753655948,0.0014385895957173568,0.0014791892711835141,
   0.0015213238191996614,0.0015650661810317963,0.0016104934837885981,0.001657687322571147,0.0017067340645141838,
   0.0017577251766440938,0.0018107575796683735,0.0018659340300216179,0.0019233635327265298,0.0019831617878879122,
   0.0020454516739262497,0.002110363770978921,0.002178036928255416,0.0022486188795328055,0.0023222669114244575,
   0.0023991485895546133,0.002479442548330867,0.0025633393506336675,0.0026510424244457313,0.002742769084235,
   0.0028387516457943815,0.0029392386442435347,0.0030444961660280956,0.003154809307027991,0.0032704837703296783,
   0.0033918476188513494,0.0035192531998631655,0.00365307926054877,0.0037937332761471616,0.003941654014939118,
   0.00409731436745065,0.0042612244707968,0.004433935163151941,0.0046160418079890615,0.004808188533075856,0.0050110729353623155,
   0.005225451309975515,0.0054521444697090025,0.0056920442308420974,0.005946120652068382,0.006215430126014334,
   0.006501124437600008,0.00680446092070622,0.007126813864712588,0.007469687345993064,0.007834729687044038,0.00822374977835156,
   0.00863873553631827,0.00908187481570728,0.009555579148505829,0.01006251074455828,0.0106056132648355,0.011188146968337661,
   0.011813728941498652,0.012486379248439126,0.013210573996293905,0.013991306498080315,0.014834157943620526,0.015745379266211335,
   0.016731986230784316,0.017801870183039872,0.018963927407269193,0.02022821066724284,0.021606107280277058,0.023110549038761655,
   0.024756260496891908,0.02656005364907936,0.028541178926603447,0.030721744843400595,0.03312722167935155,0.03578704849750366,
   0.03873536781388351,0.04201191872880512,0.045663127765211343,0.049743447693639815,0.05431700914731336,0.05945966907716407,
   0.06526156578101447,0.07183032477168075,0.07929510653451066,0.08781175113621649,0.0975693627088662,0.10879879928126432,0.115,
   0.12178370533810129,0.13,0.1368749683040054,0.145,0.15450982970874225,0.17523738873197545,0.19975298003522046,
   0.22894501461505934,0.26395954156353574,0.3};

#ifdef LAL_PTHREAD_LOCK
static pthread_once_t SEOBNRv2ROMEffectiveSpin_is_initialized = PTHREAD_ONCE_INIT;
#endif

/*************** type definitions ******************/

typedef struct tagSEOBNRROMdata_coeff
{
  gsl_vector* c_amp;
  gsl_vector* c_phi;
} SEOBNRROMdata_coeff;

struct tagSEOBNRROMdata
{
  UINT4 setup;
  gsl_vector* cvec_amp;
  gsl_vector* cvec_phi;
  gsl_matrix *Bamp;
  gsl_matrix *Bphi;
  gsl_vector* cvec_amp_pre;
};
typedef struct tagSEOBNRROMdata SEOBNRROMdata;

static SEOBNRROMdata __lalsim_SEOBNRv2ROMSS_data;

typedef struct tagSplineData
{
  gsl_bspline_workspace *bwx;
  gsl_bspline_workspace *bwy;
  int ncx, ncy;
} SplineData;

/**************** Internal functions **********************/

static void SEOBNRv2ROMEffectiveSpin_Init_LALDATA(void);
static int SEOBNRv2ROMEffectiveSpin_Init(const char dir[]);
static bool SEOBNRv2ROMEffectiveSpin_IsSetup(void);

static int SEOBNRROMdata_Init(SEOBNRROMdata *romdata, const char dir[]);
static void SEOBNRROMdata_Cleanup(SEOBNRROMdata *romdata);

static void SEOBNRROMdata_coeff_Init(SEOBNRROMdata_coeff **romdatacoeff);
static void SEOBNRROMdata_coeff_Cleanup(SEOBNRROMdata_coeff *romdatacoeff);

static size_t NextPow2(const size_t n);
static void SplineData_Destroy(SplineData *splinedata);
static void SplineData_Init(SplineData **splinedata);

static int read_vector(const char dir[], const char fname[], gsl_vector *v);
static int read_matrix(const char dir[], const char fname[], gsl_matrix *m);

static int load_data(const char dir[], gsl_vector *cvec_amp, gsl_vector *cvec_phi, gsl_matrix *Bamp, gsl_matrix *Bphi, gsl_vector *cvec_amp_pre);

static int TP_Spline_interpolation_2d(
  REAL8 q,                  // Input: q-value for which projection coefficients should be evaluated
  REAL8 chi,                // Input: chi-value for which projection coefficients should be evaluated
  gsl_vector *cvec_amp,     // Input: data for spline coefficients for amplitude
  gsl_vector *cvec_phi,     // Input: data for spline coefficients for phase
  gsl_vector *cvec_amp_pre, // Input: data for spline coefficients for amplitude prefactor
  gsl_vector *c_amp,        // Output: interpolated projection coefficients for amplitude
  gsl_vector *c_phi,        // Output: interpolated projection coefficients for phase
  REAL8 *amp_pre            // Output: interpolated amplitude prefactor
);


/**
 * Core function for computing the ROM waveform.
 * Interpolate projection coefficient data and evaluate coefficients at desired (q, chi).
 * Construct 1D splines for amplitude and phase.
 * Compute strain waveform from amplitude and phase.
*/
static int SEOBNRv2ROMEffectiveSpinCore(
  COMPLEX16FrequencySeries **hptilde,
  COMPLEX16FrequencySeries **hctilde,
  double phiRef,
  double fRef,
  double distance,
  double inclination,
  double Mtot_sec,
  double eta,
  double chi,
  const REAL8Sequence *freqs, /* Frequency points at which to evaluate the waveform (Hz) */
  double deltaF
  /* If deltaF > 0, the frequency points given in freqs are uniformly spaced with
   * spacing deltaF. Otherwise, the frequency points are spaced non-uniformly.
   * Then we will use deltaF = 0 to create the frequency series we return. */
);

// Auxiliary function to perform setup of phase spline for t(f) and f(t) functions
static int SEOBNRv2ROMEffectiveSpinTimeFrequencySetup(
  gsl_spline **spline_phi,                      // phase spline
  gsl_interp_accel **acc_phi,                   // phase spline accelerator
  REAL8 *Mf_final,                              // ringdown frequency in Mf
  REAL8 *Mtot_sec,                              // total mass in seconds
  REAL8 m1SI,                                   // Mass of companion 1 (kg)
  REAL8 m2SI,                                   // Mass of companion 2 (kg)
  REAL8 chi                                     // Effective aligned spin
);

UNUSED static REAL8 Interpolate_Coefficent_Tensor(
  gsl_vector *v,
  REAL8 eta,
  REAL8 chi1,
  REAL8 chi2,
  int ncy,
  int ncz,
  gsl_bspline_workspace *bwx,
  gsl_bspline_workspace *bwy,
  gsl_bspline_workspace *bwz
);

/********************* Definitions begin here ********************/

/** Setup SEOBNRv2ROMEffectiveSpin model using data files installed in dir
 */
static int SEOBNRv2ROMEffectiveSpin_Init(const char dir[]) {
  if(__lalsim_SEOBNRv2ROMSS_data.setup)
    XLAL_ERROR(XLAL_EFAILED, "Error: SEOBNRROMdata was already set up!");

  SEOBNRROMdata_Init(&__lalsim_SEOBNRv2ROMSS_data, dir);

  if(__lalsim_SEOBNRv2ROMSS_data.setup) {
    return(XLAL_SUCCESS);
  }
  else {
    return(XLAL_EFAILED);
  }
}

/** Helper function to check if the SEOBNRv2ROMEffectiveSpin model has been initialised */
static bool SEOBNRv2ROMEffectiveSpin_IsSetup(void) {
  if(__lalsim_SEOBNRv2ROMSS_data.setup)
    return true;
  else
    return false;
}

// Read binary ROM data for basis functions and coefficients
static int load_data(const char dir[], gsl_vector *cvec_amp, gsl_vector *cvec_phi, gsl_matrix *Bamp, gsl_matrix *Bphi, gsl_vector *cvec_amp_pre) {
  // Load binary data for amplitude and phase spline coefficients as computed in Mathematica
  int ret = XLAL_SUCCESS;
  ret |= read_vector(dir, "SEOBNRv2ROM_SS_Amp_ciall.dat", cvec_amp);
  ret |= read_vector(dir, "SEOBNRv2ROM_SS_Phase_ciall.dat", cvec_phi);
  ret |= read_matrix(dir, "SEOBNRv2ROM_SS_Bamp_bin.dat", Bamp);
  ret |= read_matrix(dir, "SEOBNRv2ROM_SS_Bphase_bin.dat", Bphi);
  ret |= read_vector(dir, "SEOBNRv2ROM_SS_AmpPrefac_ci.dat", cvec_amp_pre);
  return(ret);
}

static void SplineData_Init( SplineData **splinedata )
{
  if(!splinedata) exit(1);
  if(*splinedata) SplineData_Destroy(*splinedata);

  (*splinedata)=XLALCalloc(1,sizeof(SplineData));

  const int ncx = 50;    // points in eta + 2
  const int ncy = 50;    // points in chi + 2
  (*splinedata)->ncx = ncx;
  (*splinedata)->ncy = ncy;

  // Set up B-spline basis for desired knots
  double etavec[] = {0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, \
    0.02, 0.0225, 0.025, 0.0275, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, \
    0.06, 0.065, 0.07, 0.075, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, \
    0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.245, \
    0.247, 0.248, 0.2495, 0.2498, 0.2499, 0.25};

 double chivec[] = {-1., -0.9, -0.8, -0.67115, -0.5423, -0.430639, -0.318978, -0.234489, \
    -0.15, -0.075, 0., 0.130167, 0.260333, 0.34325, 0.426167, 0.514612, \
    0.603056, 0.647278, 0.6915, 0.723561, 0.755622, 0.782156, 0.808689, \
    0.829695, 0.8507, 0.861756, 0.872811, 0.888289, 0.903767, 0.904873, \
    0.905978, 0.911506, 0.917033, 0.922561, 0.928089, 0.931406, 0.934722, \
    0.935828, 0.936933, 0.942461, 0.947989, 0.955728, 0.963467, 0.9701, \
    0.976733, 0.983367, 0.987, 0.99};

  const size_t nbreak_x = ncx-2;  // must have nbreak = n -2 for cubic splines
  const size_t nbreak_y = ncy-2;  // must have nbreak = n -2 for cubic splines

  // allocate a cubic bspline workspace (k = 4)
  gsl_bspline_workspace *bwx = gsl_bspline_alloc(4, nbreak_x);
  gsl_bspline_workspace *bwy = gsl_bspline_alloc(4, nbreak_y);

  // set breakpoints (and thus knots by hand)
  gsl_vector *breakpts_x = gsl_vector_alloc(nbreak_x);
  gsl_vector *breakpts_y = gsl_vector_alloc(nbreak_y);
  for (UINT4 i=0; i<nbreak_x; i++)
    gsl_vector_set(breakpts_x, i, etavec[i]);
  for (UINT4 j=0; j<nbreak_y; j++)
    gsl_vector_set(breakpts_y, j, chivec[j]);

  gsl_bspline_knots(breakpts_x, bwx);
  gsl_bspline_knots(breakpts_y, bwy);

  gsl_vector_free(breakpts_x);
  gsl_vector_free(breakpts_y);

  (*splinedata)->bwx=bwx;
  (*splinedata)->bwy=bwy;
}

static void SplineData_Destroy(SplineData *splinedata)
{
  if(!splinedata) return;
  if(splinedata->bwx) gsl_bspline_free(splinedata->bwx);
  if(splinedata->bwy) gsl_bspline_free(splinedata->bwy);
  XLALFree(splinedata);
}

// Interpolate projection coefficients for amplitude and phase over the parameter space (q, chi).
// The multi-dimensional interpolation is carried out via a tensor product decomposition.
static int TP_Spline_interpolation_2d(
  REAL8 eta,                // Input: eta-value for which projection coefficients should be evaluated
  REAL8 chi,                // Input: chi-value for which projection coefficients should be evaluated
  gsl_vector *cvec_amp,     // Input: data for spline coefficients for amplitude
  gsl_vector *cvec_phi,     // Input: data for spline coefficients for phase
  gsl_vector *cvec_amp_pre, // Input: data for spline coefficients for amplitude prefactor
  gsl_vector *c_amp,        // Output: interpolated projection coefficients for amplitude
  gsl_vector *c_phi,        // Output: interpolated projection coefficients for phase
  REAL8 *amp_pre            // Output: interpolated amplitude prefactor
) {

  SplineData *splinedata=NULL;
  SplineData_Init(&splinedata);
  gsl_bspline_workspace *bwx=splinedata->bwx;
  gsl_bspline_workspace *bwy=splinedata->bwy;

  int ncx = splinedata->ncx; // points in eta
  int ncy = splinedata->ncy; // points in chi
  int N = ncx*ncy;  // size of the data matrix for one SVD-mode

  // Evaluate the TP spline for all SVD modes - amplitude
  for (int k=0; k<nk_amp; k++) { // For each SVD mode
    gsl_vector v = gsl_vector_subvector(cvec_amp, k*N, N).vector; // Pick out the coefficient matrix corresponding to the k-th SVD mode.
    REAL8 csum = Interpolate_Coefficent_Matrix(&v, eta, chi, ncx, ncy, bwx, bwy);
    gsl_vector_set(c_amp, k, csum);
  }

  // Evaluate the TP spline for all SVD modes - phase
  for (int k=0; k<nk_phi; k++) {  // For each SVD mode
    gsl_vector v = gsl_vector_subvector(cvec_phi, k*N, N).vector; // Pick out the coefficient matrix corresponding to the k-th SVD mode.
    REAL8 csum = Interpolate_Coefficent_Matrix(&v, eta, chi, ncx, ncy, bwx, bwy);
    gsl_vector_set(c_phi, k, csum);
  }

  // Evaluate the TP spline for the amplitude prefactor
  *amp_pre = Interpolate_Coefficent_Matrix(cvec_amp_pre, eta, chi, ncx, ncy, bwx, bwy);

  SplineData_Destroy(splinedata);

  return(0);
}

/* Set up a new ROM model, using data contained in dir */
static int SEOBNRROMdata_Init(SEOBNRROMdata *romdata, const char dir[]) {
  // set up ROM
  int ncx = 50;     // points in eta + 2
  int ncy = 50;     // points in chi + 2
  int N = ncx*ncy;  // size of the data matrix for one SVD-mode

  int ret = XLAL_FAILURE;

  /* Create storage for structures */
  if(romdata->setup)
  {
    XLAL_PRINT_WARNING("WARNING: You tried to setup the SEOBNRv2ROMEffectiveSpin model that was already initialised. Ignoring");
    return (XLAL_FAILURE);
  }

  gsl_set_error_handler(&err_handler);
  (romdata)->cvec_amp = gsl_vector_alloc(N*nk_amp);
  (romdata)->cvec_phi = gsl_vector_alloc(N*nk_phi);
  (romdata)->Bamp = gsl_matrix_alloc(nk_amp, nk_amp);
  (romdata)->Bphi = gsl_matrix_alloc(nk_phi, nk_phi);
  (romdata)->cvec_amp_pre = gsl_vector_alloc(N);
  ret=load_data(dir, (romdata)->cvec_amp, (romdata)->cvec_phi, (romdata)->Bamp, (romdata)->Bphi, (romdata)->cvec_amp_pre);

  if(XLAL_SUCCESS==ret) romdata->setup=1;
  else SEOBNRROMdata_Cleanup(romdata);

  return (ret);
}


/* Deallocate contents of the given SEOBNRROMdata structure */
static void SEOBNRROMdata_Cleanup(SEOBNRROMdata *romdata) {
  if(romdata->cvec_amp) gsl_vector_free(romdata->cvec_amp);
  if(romdata->cvec_phi) gsl_vector_free(romdata->cvec_phi);
  if(romdata->Bamp) gsl_matrix_free(romdata->Bamp);
  if(romdata->Bphi) gsl_matrix_free(romdata->Bphi);
  if(romdata->cvec_amp_pre) gsl_vector_free(romdata->cvec_amp_pre);
  romdata->setup=0;
}

/* Structure for internal use */
static void SEOBNRROMdata_coeff_Init(SEOBNRROMdata_coeff **romdatacoeff) {

  if(!romdatacoeff) exit(1);
  /* Create storage for structures */
  if(!*romdatacoeff) *romdatacoeff=XLALCalloc(1,sizeof(SEOBNRROMdata_coeff));
  else
  {
    SEOBNRROMdata_coeff_Cleanup(*romdatacoeff);
  }

  (*romdatacoeff)->c_amp = gsl_vector_alloc(nk_amp);
  (*romdatacoeff)->c_phi = gsl_vector_alloc(nk_phi);
}

/* Deallocate contents of the given SEOBNRROMdata_coeff structure */
static void SEOBNRROMdata_coeff_Cleanup(SEOBNRROMdata_coeff *romdatacoeff) {
  if(romdatacoeff->c_amp) gsl_vector_free(romdatacoeff->c_amp);
  if(romdatacoeff->c_phi) gsl_vector_free(romdatacoeff->c_phi);
  XLALFree(romdatacoeff);
}

/* Return the closest higher power of 2  */
static size_t NextPow2(const size_t n) {
  return 1 << (size_t) ceil(log2(n));
}

/**
 * Core function for computing the ROM waveform.
 * Interpolate projection coefficient data and evaluate coefficients at desired (q, chi).
 * Construct 1D splines for amplitude and phase.
 * Compute strain waveform from amplitude and phase.
*/
static int SEOBNRv2ROMEffectiveSpinCore(
  COMPLEX16FrequencySeries **hptilde,
  COMPLEX16FrequencySeries **hctilde,
  double phiRef,
  double fRef,
  double distance,
  double inclination,
  double Mtot_sec,
  double eta,
  double chi,
  const REAL8Sequence *freqs_in, /* Frequency points at which to evaluate the waveform (Hz) */
  double deltaF
  /* If deltaF > 0, the frequency points given in freqs are uniformly spaced with
   * spacing deltaF. Otherwise, the frequency points are spaced non-uniformly.
   * Then we will use deltaF = 0 to create the frequency series we return. */
)
{
  /* Check output arrays */
  if(!hptilde || !hctilde)
    XLAL_ERROR(XLAL_EFAULT);
  SEOBNRROMdata *romdata=&__lalsim_SEOBNRv2ROMSS_data;
  if(*hptilde || *hctilde)
    XLAL_ERROR(XLAL_EFAULT, "(*hptilde) and (*hctilde) are supposed to be NULL, but got %p and %p",(*hptilde),(*hctilde));
  int retcode=0;

  // 'Nudge' parameter values to allowed boundary values if close by
  if (eta > 0.25) nudge(&eta, 0.25, 1e-6);
  if (eta < 0.01) nudge(&eta, 0.01, 1e-6);
  if (chi < -1.0) nudge(&chi, -1.0, 1e-6);
  if (chi > 0.99) nudge(&chi, 0.99, 1e-6);

  if (chi < -1.0 || chi > 0.99)
    XLAL_ERROR(XLAL_EDOM, "XLAL Error - %s: chi (%f) smaller than -1 or larger than 0.99!\nSEOBNRv2ROMEffectiveSpin is only available for spins in the range -1 <= a/M <= 0.99.\n", __func__,chi);

  if (eta < 0.01 || eta > 0.25)
    XLAL_ERROR(XLAL_EDOM, "XLAL Error - %s: eta (%f) smaller than 0.01 or unphysical!\nSEOBNRv2ROMEffectiveSpin is only available for eta in the range 0.01 <= eta <= 0.25.\n", __func__,eta);

  /* Find frequency bounds */
  if (!freqs_in) XLAL_ERROR(XLAL_EFAULT);
  double fLow  = freqs_in->data[0];
  double fHigh = freqs_in->data[freqs_in->length - 1];

  if(fRef==0.0)
    fRef=fLow;

  /* Convert to geometric units for frequency */
  double Mf_ROM_min = fmax(gA[0], gPhi[0]);               // lowest allowed geometric frequency for ROM
  double Mf_ROM_max = fmin(gA[nk_amp-1], gPhi[nk_phi-1]); // highest allowed geometric frequency for ROM
  double fLow_geom = fLow * Mtot_sec;
  double fHigh_geom = fHigh * Mtot_sec;
  double fRef_geom = fRef * Mtot_sec;
  double deltaF_geom = deltaF * Mtot_sec;

  // Enforce allowed geometric frequency range
  if (fLow_geom < Mf_ROM_min)
    XLAL_ERROR(XLAL_EDOM, "Starting frequency Mflow=%g is smaller than lowest frequency in ROM Mf=%g. Starting at lowest frequency in ROM.\n", fLow_geom, Mf_ROM_min);
  if (fHigh_geom == 0)
    fHigh_geom = Mf_ROM_max;
  else if (fHigh_geom > Mf_ROM_max) {
	  XLALPrintWarning("Maximal frequency Mf_high=%g is greater than highest ROM frequency Mf_ROM_Max=%g. Using Mf_high=Mf_ROM_Max.", fHigh_geom, Mf_ROM_max);
	  fHigh_geom = Mf_ROM_max;
  }
  else if (fHigh_geom < Mf_ROM_min)
    XLAL_ERROR(XLAL_EDOM, "End frequency %g is smaller than starting frequency %g!\n", fHigh_geom, fLow_geom);
  if (fRef_geom > Mf_ROM_max) {
	  XLALPrintWarning("Reference frequency Mf_ref=%g is greater than maximal frequency in ROM Mf=%g. Starting at maximal frequency in ROM.\n", fRef_geom, Mf_ROM_max);
    fRef_geom = Mf_ROM_max; // If fref > fhigh we reset fref to default value of cutoff frequency.
  }
  if (fRef_geom < Mf_ROM_min) {
    XLAL_PRINT_WARNING("Reference frequency Mf_ref=%g is smaller than lowest frequency in ROM Mf=%g. Starting at lowest frequency in ROM.\n", fLow_geom, Mf_ROM_min);
    fRef_geom = Mf_ROM_min;
  }

  /* Internal storage for w.f. coefficiencts */
  SEOBNRROMdata_coeff *romdata_coeff=NULL;
  SEOBNRROMdata_coeff_Init(&romdata_coeff);
  REAL8 amp_pre;

  /* Interpolate projection coefficients and evaluate them at (eta,chi) */
  retcode=TP_Spline_interpolation_2d(
    eta,                       // Input: eta-value for which projection coefficients should be evaluated
    chi,                       // Input: chi-value for which projection coefficients should be evaluated
    romdata->cvec_amp,         // Input: data for spline coefficients for amplitude
    romdata->cvec_phi,         // Input: data for spline coefficients for phase
    romdata->cvec_amp_pre,     // Input: data for spline coefficients for amplitude prefactor
    romdata_coeff->c_amp,      // Output: interpolated projection coefficients for amplitude
    romdata_coeff->c_phi,      // Output: interpolated projection coefficients for phase
    &amp_pre                   // Output: interpolated amplitude prefactor
  );

  if(retcode!=0) {
    SEOBNRROMdata_coeff_Cleanup(romdata_coeff);
    XLAL_ERROR(retcode, "Parameter-space interpolation failed.");
  }

  // Compute function values of amplitude an phase on sparse frequency points by evaluating matrix vector products
  // amp_pts = B_A^T . c_A
  // phi_pts = B_phi^T . c_phi
  gsl_vector* amp_f = gsl_vector_alloc(nk_amp);
  gsl_vector* phi_f = gsl_vector_alloc(nk_phi);
  gsl_blas_dgemv(CblasTrans, 1.0, romdata->Bamp, romdata_coeff->c_amp, 0.0, amp_f);
  gsl_blas_dgemv(CblasTrans, 1.0, romdata->Bphi, romdata_coeff->c_phi, 0.0, phi_f);

  // Setup 1d splines in frequency
  gsl_interp_accel *acc_amp = gsl_interp_accel_alloc();
  gsl_spline *spline_amp = gsl_spline_alloc(gsl_interp_cspline, nk_amp);
  gsl_spline_init(spline_amp, gA, gsl_vector_const_ptr(amp_f,0), nk_amp);

  gsl_interp_accel *acc_phi = gsl_interp_accel_alloc();
  gsl_spline *spline_phi = gsl_spline_alloc(gsl_interp_cspline, nk_phi);
  gsl_spline_init(spline_phi, gPhi, gsl_vector_const_ptr(phi_f,0), nk_phi);

  size_t npts = 0;
  LIGOTimeGPS tC = {0, 0};
  UINT4 offset = 0; // Index shift between freqs and the frequency series
  REAL8Sequence *freqs = NULL;
  if (deltaF > 0)  { // freqs contains uniform frequency grid with spacing deltaF; we start at frequency 0
    /* Set up output array with size closest power of 2 */
    npts = NextPow2(fHigh_geom / deltaF_geom) + 1;
    if (fHigh_geom < fHigh * Mtot_sec) /* Resize waveform if user wants f_max larger than cutoff frequency */
      npts = NextPow2(fHigh * Mtot_sec / deltaF_geom) + 1;

    XLALGPSAdd(&tC, -1. / deltaF);  /* coalesce at t=0 */
    *hptilde = XLALCreateCOMPLEX16FrequencySeries("hptilde: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, npts);
    *hctilde = XLALCreateCOMPLEX16FrequencySeries("hctilde: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, npts);

    // Recreate freqs using only the lower and upper bounds
    UINT4 iStart = (UINT4) ceil(fLow_geom / deltaF_geom);
    UINT4 iStop = (UINT4) ceil(fHigh_geom / deltaF_geom);
    freqs = XLALCreateREAL8Sequence(iStop - iStart);
    for (UINT4 i=iStart; i<iStop; i++)
      freqs->data[i-iStart] = i*deltaF_geom;

    offset = iStart;
  } else { // freqs contains frequencies with non-uniform spacing; we start at lowest given frequency
    npts = freqs_in->length;
    *hptilde = XLALCreateCOMPLEX16FrequencySeries("hptilde: FD waveform", &tC, fLow, 0, &lalStrainUnit, npts);
    *hctilde = XLALCreateCOMPLEX16FrequencySeries("hctilde: FD waveform", &tC, fLow, 0, &lalStrainUnit, npts);
    offset = 0;

    freqs = XLALCreateREAL8Sequence(freqs_in->length);
    for (UINT4 i=0; i<freqs_in->length; i++)
      freqs->data[i] = freqs_in->data[i] * Mtot_sec;
  }

  if (!(*hptilde) || !(*hctilde)) {
      XLALDestroyREAL8Sequence(freqs);
      gsl_spline_free(spline_amp);
      gsl_spline_free(spline_phi);
      gsl_interp_accel_free(acc_amp);
      gsl_interp_accel_free(acc_phi);
      gsl_vector_free(amp_f);
      gsl_vector_free(phi_f);
      SEOBNRROMdata_coeff_Cleanup(romdata_coeff);
      XLAL_ERROR(XLAL_EFUNC, "Waveform allocation failed.");
  }
  memset((*hptilde)->data->data, 0, npts * sizeof(COMPLEX16));
  memset((*hctilde)->data->data, 0, npts * sizeof(COMPLEX16));

  XLALUnitMultiply(&(*hptilde)->sampleUnits, &(*hptilde)->sampleUnits, &lalSecondUnit);
  XLALUnitMultiply(&(*hctilde)->sampleUnits, &(*hctilde)->sampleUnits, &lalSecondUnit);

  COMPLEX16 *pdata=(*hptilde)->data->data;
  COMPLEX16 *cdata=(*hctilde)->data->data;

  REAL8 cosi = cos(inclination);
  REAL8 pcoef = 0.5*(1.0 + cosi*cosi);
  REAL8 ccoef = cosi;

  REAL8 s = 1.0/2.0; // Scale polarization amplitude so that strain agrees with FFT of SEOBNRv2
  double Mtot = Mtot_sec / LAL_MTSUN_SI;
  double amp0 = Mtot * amp_pre * Mtot_sec * LAL_MRSUN_SI / (distance); // Correct overall amplitude to undo mass-dependent scaling used in single-spin ROM

  // Evaluate reference phase for setting phiRef correctly
  double phase_change = gsl_spline_eval(spline_phi, fRef_geom, acc_phi) - 2*phiRef;

  // Assemble waveform from aplitude and phase
  for (UINT4 i=0; i<freqs->length; i++) { // loop over frequency points in sequence
    double f = freqs->data[i];
    if (f > Mf_ROM_max) continue; // We're beyond the highest allowed frequency; since freqs may not be ordered, we'll just skip the current frequency and leave zero in the buffer
    int j = i + offset; // shift index for frequency series if needed
    double A = gsl_spline_eval(spline_amp, f, acc_amp);
    double phase = gsl_spline_eval(spline_phi, f, acc_phi) - phase_change;
    COMPLEX16 htilde = s*amp0*A * cexp(I*phase);
    pdata[j] =      pcoef * htilde;
    cdata[j] = -I * ccoef * htilde;
  }

  /* Correct phasing so we coalesce at t=0 (with the definition of the epoch=-1/deltaF above) */

  // Get SEOBNRv2 ringdown frequency for 22 mode
  // XLALSimInspiralGetFinalFreq wants masses in SI units, so unfortunately we need to convert back
  double q = (1.0 + sqrt(1.0 - 4.0*eta) - 2.0*eta) / (2.0*eta);
  double Mtot_SI = Mtot_sec / LAL_MTSUN_SI * LAL_MSUN_SI;
  double m1_SI = Mtot_SI * 1.0/(1.0+q);
  double m2_SI = Mtot_SI * q/(1.0+q);
  double Mf_final = XLALSimInspiralGetFinalFreq(m1_SI, m2_SI, 0,0,chi, 0,0,chi, SEOBNRv2) * Mtot_sec;

  UINT4 L = freqs->length;
  // prevent gsl interpolation errors
  if (Mf_final > freqs->data[L-1])
    Mf_final = freqs->data[L-1];
  if (Mf_final < freqs->data[0])
  {
      XLALDestroyREAL8Sequence(freqs);
      gsl_spline_free(spline_amp);
      gsl_spline_free(spline_phi);
      gsl_interp_accel_free(acc_amp);
      gsl_interp_accel_free(acc_phi);
      gsl_vector_free(amp_f);
      gsl_vector_free(phi_f);
      SEOBNRROMdata_coeff_Cleanup(romdata_coeff);
      XLAL_ERROR(XLAL_EDOM, "f_ringdown < f_min");
  }

  // Time correction is t(f_final) = 1/(2pi) dphi/df (f_final)
  // We compute the dimensionless time correction t/M since we use geometric units.
  REAL8 t_corr = gsl_spline_eval_deriv(spline_phi, Mf_final, acc_phi) / (2*LAL_PI);

  // Now correct phase
  for (UINT4 i=0; i<freqs->length; i++) { // loop over frequency points in sequence
    double f = freqs->data[i] - fRef_geom;
    int j = i + offset; // shift index for frequency series if needed
    pdata[j] *= cexp(-2*LAL_PI * I * f * t_corr);
    cdata[j] *= cexp(-2*LAL_PI * I * f * t_corr);
  }

  XLALDestroyREAL8Sequence(freqs);

  gsl_spline_free(spline_amp);
  gsl_spline_free(spline_phi);
  gsl_interp_accel_free(acc_amp);
  gsl_interp_accel_free(acc_phi);
  gsl_vector_free(amp_f);
  gsl_vector_free(phi_f);
  SEOBNRROMdata_coeff_Cleanup(romdata_coeff);

  return(XLAL_SUCCESS);
}

/**
 * @addtogroup LALSimIMRSEOBNRROM_c
 *
 * @{
 *
 * @name SEOBNRv2 Reduced Order Model (Effective Spin)
 *
 * @author Michael Puerrer, John Veitch
 *
 * @brief C code for SEOBNRv2 reduced order model (equal spin version).
 * See CQG 31 195010, 2014, arXiv:1402.4146 for details.
 *
 * This is a frequency domain model that approximates the time domain SEOBNRv2 model with equal spins.
 *
 * The binary data files are available at https://dcc.ligo.org/T1400701-v1.
 * Put the untared data into a location in your LAL_DATA_PATH.
 *
 * @note Note that due to its construction the iFFT of the ROM has a small (~ 20 M) offset
 * in the peak time that scales with total mass as compared to the time-domain SEOBNRv2 model.
 *
 * @note Parameter ranges:
 *   * 0.01 <= eta <= 0.25
 *   * -1 <= chi <= 0.99
 *   * Mtot >= 1.4Msun
 *
 *  Equal spin chi = chi1 = chi2 in terms of aligned component spins chi1, chi2.
 *  Symmetric mass-ratio eta = m1*m2/(m1+m2)^2.
 *  Total mass Mtot.
 *
 * @{
 */


/**
 * Compute waveform in LAL format at specified frequencies for the SEOBNRv2_ROM_EffectiveSpin model.
 *
 * XLALSimIMRSEOBNRv2ROMEffectiveSpin() returns the plus and cross polarizations as a complex
 * frequency series with equal spacing deltaF and contains zeros from zero frequency
 * to the starting frequency and zeros beyond the cutoff frequency in the ringdown.
 *
 * In contrast, XLALSimIMRSEOBNRv2ROMEffectiveSpinFrequencySequence() returns a
 * complex frequency series with entries exactly at the frequencies specified in
 * the sequence freqs (which can be unequally spaced). No zeros are added.
 *
 * If XLALSimIMRSEOBNRv2ROMEffectiveSpinFrequencySequence() is called with frequencies that
 * are beyond the maxium allowed geometric frequency for the ROM, zero strain is returned.
 * It is not assumed that the frequency sequence is ordered.
 *
 * This function is designed as an entry point for reduced order quadratures.
 */
int XLALSimIMRSEOBNRv2ROMEffectiveSpinFrequencySequence(
  struct tagCOMPLEX16FrequencySeries **hptilde, /**< Output: Frequency-domain waveform h+ */
  struct tagCOMPLEX16FrequencySeries **hctilde, /**< Output: Frequency-domain waveform hx */
  const REAL8Sequence *freqs,                   /**< Frequency points at which to evaluate the waveform (Hz), need to be strictly monotonically increasing */
  REAL8 phiRef,                                 /**< Orbital phase at reference time */
  REAL8 fRef,                                   /**< Reference frequency (Hz); 0 defaults to fLow */
  REAL8 distance,                               /**< Distance of source (m) */
  REAL8 inclination,                            /**< Inclination of source (rad) */
  REAL8 m1SI,                                   /**< Mass of companion 1 (kg) */
  REAL8 m2SI,                                   /**< Mass of companion 2 (kg) */
  REAL8 chi)                                    /**< Effective aligned spin */
{

  /* Get masses in terms of solar mass */
  double mass1 = m1SI / LAL_MSUN_SI;
  double mass2 = m2SI / LAL_MSUN_SI;
  double Mtot = mass1+mass2;
  double eta = mass1 * mass2 / (Mtot*Mtot);    /* Symmetric mass-ratio */
  double Mtot_sec = Mtot * LAL_MTSUN_SI; /* Total mass in seconds */

  if (!freqs) XLAL_ERROR(XLAL_EFAULT);

  // Load ROM data if not loaded already
#ifdef LAL_PTHREAD_LOCK
  (void) pthread_once(&SEOBNRv2ROMEffectiveSpin_is_initialized, SEOBNRv2ROMEffectiveSpin_Init_LALDATA);
#else
  SEOBNRv2ROMEffectiveSpin_Init_LALDATA();
#endif

  if(!SEOBNRv2ROMEffectiveSpin_IsSetup()) XLAL_ERROR(XLAL_EFAILED,"Error setting up SEOBNRv2ROMEffectiveSpin - check your $LAL_DATA_PATH\n");

  // Call the internal core function with deltaF = 0 to indicate that freqs is non-uniformly
  // spaced and we want the strain only at these frequencies
  int retcode = SEOBNRv2ROMEffectiveSpinCore(hptilde,hctilde,
            phiRef, fRef, distance, inclination, Mtot_sec, eta, chi, freqs, 0);

  return(retcode);
}

/**
 * Compute waveform in LAL format for the SEOBNRv2_ROM_EffectiveSpin model.
 *
 * Returns the plus and cross polarizations as a complex frequency series with
 * equal spacing deltaF and contains zeros from zero frequency to the starting
 * frequency fLow and zeros beyond the cutoff frequency fHigh to the next power of 2 in
 * the size of the frequency series.
 */
int XLALSimIMRSEOBNRv2ROMEffectiveSpin(
  struct tagCOMPLEX16FrequencySeries **hptilde, /**< Output: Frequency-domain waveform h+ */
  struct tagCOMPLEX16FrequencySeries **hctilde, /**< Output: Frequency-domain waveform hx */
  REAL8 phiRef,                                 /**< Orbital phase at reference time */
  REAL8 deltaF,                                 /**< Sampling frequency (Hz) */
  REAL8 fLow,                                   /**< Starting GW frequency (Hz) */
  REAL8 fHigh,                                  /**< End frequency; 0 defaults to Mf=0.14 */
  REAL8 fRef,                                   /**< Reference frequency (Hz); 0 defaults to fLow */
  REAL8 distance,                               /**< Distance of source (m) */
  REAL8 inclination,                            /**< Inclination of source (rad) */
  REAL8 m1SI,                                   /**< Mass of companion 1 (kg) */
  REAL8 m2SI,                                   /**< Mass of companion 2 (kg) */
  REAL8 chi)                                    /**< Effective aligned spin */
{

  /* Get masses in terms of solar mass */
  double mass1 = m1SI / LAL_MSUN_SI;
  double mass2 = m2SI / LAL_MSUN_SI;
  double Mtot = mass1+mass2;
  double eta = mass1 * mass2 / (Mtot*Mtot);    /* Symmetric mass-ratio */
  double Mtot_sec = Mtot * LAL_MTSUN_SI; /* Total mass in seconds */

  if(fRef==0.0)
    fRef=fLow;

  // Load ROM data if not loaded already
#ifdef LAL_PTHREAD_LOCK
  (void) pthread_once(&SEOBNRv2ROMEffectiveSpin_is_initialized, SEOBNRv2ROMEffectiveSpin_Init_LALDATA);
#else
  SEOBNRv2ROMEffectiveSpin_Init_LALDATA();
#endif

  if(!SEOBNRv2ROMEffectiveSpin_IsSetup()) XLAL_ERROR(XLAL_EFAILED,"Error setting up SEOBNRv2ROMEffectiveSpin data - check your $LAL_DATA_PATH\n");

  // Use fLow, fHigh, deltaF to compute freqs sequence
  // Instead of building a full sequency we only transfer the boundaries and let
  // the internal core function do the rest (and properly take care of corner cases).
  REAL8Sequence *freqs = XLALCreateREAL8Sequence(2);
  freqs->data[0] = fLow;
  freqs->data[1] = fHigh;

  int retcode = SEOBNRv2ROMEffectiveSpinCore(hptilde, hctilde,
            phiRef, fRef, distance, inclination, Mtot_sec, eta, chi, freqs, deltaF);

  XLALDestroyREAL8Sequence(freqs);

  return(retcode);
}

/**
 * Compute the 'time' elapsed in the ROM waveform from a given starting frequency until the ringdown.
 * UNREVIEWED!
 *
 * The notion of elapsed 'time' (in seconds) is defined here as the difference of the
 * frequency derivative of the frequency domain phase between the ringdown frequency
 * and the starting frequency ('frequency' argument). This notion of time is similar to the
 * chirp time, but it includes both the inspiral and the merger ringdown part of SEOBNRv2.
 *
 * The allowed frequency range for the starting frequency in geometric frequency is [0.0001, 0.3].
 * The SEOBNRv2 ringdown frequency can be obtained by calling XLALSimInspiralGetFinalFreq().
 *
 * See XLALSimIMRSEOBNRv2ROMEffectiveSpinFrequencyOfTime() for the inverse function.
 */
int XLALSimIMRSEOBNRv2ROMEffectiveSpinTimeOfFrequency(
  REAL8 *t,         /**< Output: time (s) at frequency */
  REAL8 frequency,  /**< Frequency (Hz) */
  REAL8 m1SI,       /**< Mass of companion 1 (kg) */
  REAL8 m2SI,       /**< Mass of companion 2 (kg) */
  REAL8 chi         /**< Equal aligned spin (chi = chi1 = chi2)*/
)
{
  // Set up phase spline
  gsl_spline *spline_phi;
  gsl_interp_accel *acc_phi;
  double Mf_final, Mtot_sec;
  int ret = SEOBNRv2ROMEffectiveSpinTimeFrequencySetup(&spline_phi, &acc_phi, &Mf_final, &Mtot_sec, m1SI, m2SI, chi);
  if(ret != 0)
    XLAL_ERROR(ret);

  // ROM frequency bounds in Mf
  double Mf_ROM_min = gPhi[0];
  double Mf_ROM_max = gPhi[nk_phi-1];

  // Time correction is t(f_final) = 1/(2pi) dphi/df (f_final)
  double t_corr = gsl_spline_eval_deriv(spline_phi, Mf_final, acc_phi) / (2*LAL_PI); // t_corr / M
  XLAL_PRINT_INFO("t_corr[s] = %g\n", t_corr * Mtot_sec);

  double Mf = frequency * Mtot_sec;
  if (Mf < Mf_ROM_min || Mf > Mf_ROM_max)
  {
      gsl_spline_free(spline_phi);
      gsl_interp_accel_free(acc_phi);
      XLAL_ERROR(XLAL_EDOM, "Frequency %g is outside allowed frequency range.\n", frequency);
  }
  // Compute time relative to origin at merger
  double time_M = gsl_spline_eval_deriv(spline_phi, frequency * Mtot_sec, acc_phi) / (2*LAL_PI) - t_corr;
  *t = time_M * Mtot_sec;

  gsl_spline_free(spline_phi);
  gsl_interp_accel_free(acc_phi);

  return(XLAL_SUCCESS);
}

/**
 * Compute the starting frequency so that the given amount of 'time' elapses in the ROM waveform
 * from the starting frequency until the ringdown.
 * UNREVIEWED!
 *
 * The notion of elapsed 'time' (in seconds) is defined here as the difference of the
 * frequency derivative of the frequency domain phase between the ringdown frequency
 * and the starting frequency ('frequency' argument). This notion of time is similar to the
 * chirp time, but it includes both the inspiral and the merger ringdown part of SEOBNRv2.
 *
 * If the frequency that corresponds to the specified elapsed time is lower than the
 * geometric frequency Mf=0.0001 (ROM starting frequency) or above half of the SEOBNRv2
 * ringdown frequency an error is thrown.
 * The SEOBNRv2 ringdown frequency can be obtained by calling XLALSimInspiralGetFinalFreq().
 *
 * See XLALSimIMRSEOBNRv2ROMEffectiveSpinTimeOfFrequency() for the inverse function.
 */
int XLALSimIMRSEOBNRv2ROMEffectiveSpinFrequencyOfTime(
  REAL8 *frequency,   /**< Output: Frequency (Hz) */
  REAL8 t,            /**< Time (s) at frequency */
  REAL8 m1SI,         /**< Mass of companion 1 (kg) */
  REAL8 m2SI,         /**< Mass of companion 2 (kg) */
  REAL8 chi           /**< Effective aligned spin */
)
{
  // Set up phase spline
  gsl_spline *spline_phi;
  gsl_interp_accel *acc_phi;
  double Mf_final, Mtot_sec;
  int ret = SEOBNRv2ROMEffectiveSpinTimeFrequencySetup(&spline_phi, &acc_phi, &Mf_final, &Mtot_sec, m1SI, m2SI, chi);
  if(ret != 0)
    XLAL_ERROR(ret);

  // ROM frequency bounds in Mf
  double Mf_ROM_min = gPhi[0];

  // Time correction is t(f_final) = 1/(2pi) dphi/df (f_final)
  double t_corr = gsl_spline_eval_deriv(spline_phi, Mf_final, acc_phi) / (2*LAL_PI); // t_corr / M
  XLAL_PRINT_INFO("t_corr[s] = %g\n", t_corr * Mtot_sec);

  // Assume for now that we only care about f(t) *before* merger so that f(t) - f_ringdown >= 0.
  // Assume that we only need to cover the frequency range [f_min, f_ringdown/2].
  int N = 20;
  double log_f_pts[N];
  double log_t_pts[N];
  double log_f_min   = log(Mf_ROM_min);
  double log_f_rng_2 = log(Mf_final/2.0);
  double dlog_f = (log_f_rng_2 - log_f_min) / (N-1);

  // Set up data in log-log space
  for (int i=0; i<N; i++) {
    log_f_pts[i] = log_f_rng_2 - i*dlog_f; // gsl likes the x-values to be monotonically increasing
    // Compute time relative to origin at merger
    double time_M = gsl_spline_eval_deriv(spline_phi, exp(log_f_pts[i]), acc_phi) / (2*LAL_PI) - t_corr;
    log_t_pts[i] = log(time_M * Mtot_sec);
  }

  // Check whether time is in bounds
  double t_rng_2 = exp(log_t_pts[0]);   // time of f_ringdown/2
  double t_min   = exp(log_t_pts[N-1]); // time of f_min
  if (t < t_rng_2 || t > t_min)
  {
      gsl_spline_free(spline_phi);
      gsl_interp_accel_free(acc_phi);
      XLAL_ERROR(XLAL_EDOM, "The frequency of time %g is outside allowed frequency range.\n", t);
  }

  // create new spline for data
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, N);
  gsl_spline_init(spline, log_t_pts, log_f_pts, N);

  *frequency = exp(gsl_spline_eval(spline, log(t), acc)) / Mtot_sec;

  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  gsl_spline_free(spline_phi);
  gsl_interp_accel_free(acc_phi);

  return(XLAL_SUCCESS);
}

/** @} */
/** @} */

// Auxiliary function to perform setup of phase spline for t(f) and f(t) functions
static int SEOBNRv2ROMEffectiveSpinTimeFrequencySetup(
  gsl_spline **spline_phi,                      // phase spline
  gsl_interp_accel **acc_phi,                   // phase spline accelerator
  REAL8 *Mf_final,                              // ringdown frequency in Mf
  REAL8 *Mtot_sec,                              // total mass in seconds
  REAL8 m1SI,                                   // Mass of companion 1 (kg)
  REAL8 m2SI,                                   // Mass of companion 2 (kg)
  REAL8 chi                                     // Effective aligned spin
)
{
  /* Get masses in terms of solar mass */
  double mass1 = m1SI / LAL_MSUN_SI;
  double mass2 = m2SI / LAL_MSUN_SI;
  double Mtot = mass1 + mass2;
  double eta = mass1 * mass2 / (Mtot*Mtot);    /* Symmetric mass-ratio */
  *Mtot_sec = Mtot * LAL_MTSUN_SI; /* Total mass in seconds */

  // 'Nudge' parameter values to allowed boundary values if close by
  nudge(&eta, 0.25, 1e-6);
  nudge(&eta, 0.01, 1e-6);
  nudge(&chi, -1.0, 1e-6);
  nudge(&chi, 0.99, 1e-6);

  if (chi < -1.0 || chi > 0.99)
    XLAL_ERROR(XLAL_EDOM, "XLAL Error - %s: chi (%f) smaller than -1 or larger than 0.99!\nSEOBNRv2ROMEffectiveSpin is only available for spins in the range -1 <= a/M <= 0.99.\n", __func__,chi);

  if (eta < 0.01 || eta > 0.25)
    XLAL_ERROR(XLAL_EDOM, "XLAL Error - %s: eta (%f) smaller than 0.01 or unphysical!\nSEOBNRv2ROMEffectiveSpin is only available for spins in the range 0.01 <= eta <= 0.25.\n", __func__,eta);

  // Load ROM data if not loaded already
#ifdef LAL_PTHREAD_LOCK
  (void) pthread_once(&SEOBNRv2ROMEffectiveSpin_is_initialized, SEOBNRv2ROMEffectiveSpin_Init_LALDATA);
#else
  SEOBNRv2ROMEffectiveSpin_Init_LALDATA();
#endif

  SEOBNRROMdata *romdata=&__lalsim_SEOBNRv2ROMSS_data;

  /* Internal storage for w.f. coefficiencts */
  SEOBNRROMdata_coeff *romdata_coeff=NULL;
  SEOBNRROMdata_coeff_Init(&romdata_coeff);
  REAL8 amp_pre;

  /* Interpolate projection coefficients and evaluate them at (eta,chi) */
  int retcode=TP_Spline_interpolation_2d(
    eta,                       // Input: eta-value for which projection coefficients should be evaluated
    chi,                       // Input: chi-value for which projection coefficients should be evaluated
    romdata->cvec_amp,         // Input: data for spline coefficients for amplitude
    romdata->cvec_phi,         // Input: data for spline coefficients for phase
    romdata->cvec_amp_pre,     // Input: data for spline coefficients for amplitude prefactor
    romdata_coeff->c_amp,      // Output: interpolated projection coefficients for amplitude
    romdata_coeff->c_phi,      // Output: interpolated projection coefficients for phase
    &amp_pre                   // Output: interpolated amplitude prefactor
  );

  if(retcode!=0) {
    SEOBNRROMdata_coeff_Cleanup(romdata_coeff);
    XLAL_ERROR(retcode);
  }

  // Compute function values of phase on sparse frequency points by evaluating matrix vector products
  // phi_pts = B_phi^T . c_phi
  gsl_vector* phi_f = gsl_vector_alloc(nk_phi);
  gsl_blas_dgemv(CblasTrans, 1.0, romdata->Bphi, romdata_coeff->c_phi, 0.0, phi_f);

  // Setup 1d phase spline in frequency
  *acc_phi = gsl_interp_accel_alloc();
  *spline_phi = gsl_spline_alloc(gsl_interp_cspline, nk_phi);
  gsl_spline_init(*spline_phi, gPhi, gsl_vector_const_ptr(phi_f,0), nk_phi);

  // Get SEOBNRv2 ringdown frequency for 22 mode
  double q = (1.0 + sqrt(1.0 - 4.0*eta) - 2.0*eta) / (2.0*eta);
  double Mtot_SI = *Mtot_sec / LAL_MTSUN_SI * LAL_MSUN_SI;
  double m1_SI = Mtot_SI * 1.0/(1.0+q);
  double m2_SI = Mtot_SI * q/(1.0+q);
  *Mf_final = XLALSimInspiralGetFinalFreq(m1_SI, m2_SI, 0,0,chi, 0,0,chi, SEOBNRv2) * (*Mtot_sec);

  gsl_vector_free(phi_f);
  SEOBNRROMdata_coeff_Cleanup(romdata_coeff);

  return(XLAL_SUCCESS);
}


/** Setup SEOBNRv2ROMEffectiveSpin model using data files installed in $LAL_DATA_PATH
 */
static void SEOBNRv2ROMEffectiveSpin_Init_LALDATA(void)
{
  if (SEOBNRv2ROMEffectiveSpin_IsSetup()) return;

  // If we find one ROM datafile in a directory listed in LAL_DATA_PATH,
  // then we expect the remaining datafiles to also be there.
  char datafile[] = "SEOBNRv2ROM_SS_Phase_ciall.dat";

  char *path = XLALFileResolvePathLong(datafile, PKG_DATA_DIR);
  if (path==NULL)
    XLAL_ERROR_VOID(XLAL_EIO, "Unable to resolve data file %s in $LAL_DATA_PATH\n", datafile);
  char *dir = dirname(path);
  int ret = SEOBNRv2ROMEffectiveSpin_Init(dir);
  XLALFree(path);

  if(ret!=XLAL_SUCCESS)
    XLAL_ERROR_VOID(XLAL_FAILURE, "Unable to find SEOBNRv2ROMEffectiveSpin data files in $LAL_DATA_PATH\n");
}

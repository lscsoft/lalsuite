/*
 *  LALInferenceLikelihood.c:  Bayesian Followup likelihood functions
 *
 *  Copyright (C) 2009 Ilya Mandel, Vivien Raymond, Christian Roever,
 *  Marc van der Sluys and John Veitch, Will M. Farr, Salvatore Vitale
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

#include <complex.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferencePrior.h>
#include <lal/LALInference.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/Sequence.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeFreqFFT.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_dawson.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_complex_math.h>
#include <lal/LALInferenceTemplate.h>

/* Scaling used for the burst analytic likelihood parameters */
// sqrt(bCM)/(scaling) is the sigma
  static const REAL8 burst_scaling[9] = {
    0.1,
    1.0,
    1.0,
    100.0,
    10.0/M_PI,
    10.0/M_PI,
    10.0/M_PI,
    10.0/M_PI,
    10.0};

/* burst covariance matrix used in analyic likelihoods */
/* Unimodal should give evidence -21.3,bimodal  -25.9 */
static const REAL8 bCM[9][9] = {{0.0118334553095770112635110,0.0101244893662672322265372,0.0164254129963652163726184,0.0112088766770308614906249,0.0148945067729633826014712,0.0176503111361420127189970,0.0122199881064022214394171,0.0139749282535349579614792,0.0},
{0.0101244893662672322265372,0.0223414853900214364912369,0.0223978018917797665199299,0.0155692879011910430275822,0.0201008529429579502201264,0.0229538138342478409414937,0.0180642089428136622120125,0.0155199853453866446623133,0.0},
{0.0164254129963652163726184,0.0223978018917797665199299,0.0372071319249132337336761,0.0269434973389567240797948,0.0277259617613290071380661,0.0321217958436038619751685,0.0325568864053485881870920,0.0283379829772117813879717,0.0},
{0.0112088766770308614906249,0.0155692879011910430275822,0.0269434973389567240797948,0.0219637138989302732605680,0.0202889296069694510804560,0.0225846698886396080041550,0.0238719578697630489816373,0.0197015389277761104880327,0.0},
{0.0148945067729633826014712,0.0201008529429579502201264,0.0277259617613290071380661,0.0202889296069694510804560,0.0310842562371710651181189,0.0334804433939082934923448,0.0233839662218643419555608,0.0266660912253427265228289,0.0},
{0.0176503111361420127189970,0.0229538138342478409414937,0.0321217958436038619751685,0.0225846698886396080041550,0.0334804433939082934923448,0.0380576501276081793911921,0.0264526579842353157245860,0.0312117893830374908137326,0.0},
{0.0122199881064022214394171,0.0180642089428136622120125,0.0325568864053485881870920,0.0238719578697630489816373,0.0233839662218643419555608,0.0264526579842353157245860,0.0305112248863869672810267,0.0256577816309597056543268,0.0},
{0.0139749282535349579614792,0.0155199853453866446623133,0.0283379829772117813879717,0.0197015389277761104880327,0.0266660912253427265228289,0.0312117893830374908137326,0.0256577816309597056543268,0.0308105179275282164974570,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.121}
};

const char *LALInferenceAnalyticNamesCBC[] = {"mass1","mass2","costheta_jn","phase","polarisation","rightascension","declination","logdistance","time","a_spin1","a_spin2","tilt_spin1", "tilt_spin2", "phi12", "phi_jl"};

const REAL8 LALInferenceAnalyticMeansCBC[] =
{
  16.0, 7.0, 0.0, LAL_PI, LAL_PI, LAL_PI, 0.0, 4.0, 0.0, 0.5, 0.5, LAL_PI/2.0, LAL_PI/2.0, LAL_PI, LAL_PI
};

/* Scaling used for the CBC analytic likelihood parameters */
const REAL8 scaling[15] = {
  1.0,
  1.0,
  /*20.0/M_PI, */ 20.0,
  10.0/M_PI,
  20.0/M_PI,
  10.0/M_PI,
  10.0/M_PI,
  /*0.1,*/ 10.0,
  10.0,
  10.0,
  10.0,
  20.0/M_PI,
  20.0/M_PI,
  10.0/M_PI,
  10.0/M_PI};

/* Covariance matrix for use in analytic likelihoods */
const REAL8 CM[15][15] = {{0.045991865933182365, -0.005489748382557155, -0.01025067223674548, 0.0020087713726603213, -0.0032648855847982987, -0.0034218261781145264, -0.0037173401838545774, -0.007694897715679858, 0.005260905282822458, 0.0013607957548231718, 0.001970785895702776, 0.006708452591621081, -0.005107684668720825, 0.004402554308030673, -0.00334987648531921},
  {-0.005489748382557152, 0.05478640427684032, -0.004786202916836846, -0.007930397407501268, -0.0005945107515129139, 0.004858466255616657, -0.011667819871670204, 0.003169780190169035, 0.006761345004654851, -0.0037599761532668475, 0.005571796842520162, -0.0071098291510566895, -0.004477773540640284, -0.011250694688474672, 0.007465228985669282},
  {-0.01025067223674548, -0.004786202916836844, 0.044324704403674524, -0.0010572820723801645, -0.009885693540838514, -0.0048321205972943464, -0.004576186966267275, 0.0025107211483955676, -0.010126911913571181, 0.01153595152487264, 0.005773054728678472, 0.005286558230422045, -0.0055438798694137734, 0.0044772210361854765, -0.00620416958073918},
  {0.0020087713726603213, -0.007930397407501268, -0.0010572820723801636, 0.029861342087731065, -0.007803477134405363, -0.0011466944120756021, 0.009925736654597632, -0.0007664415942207051, -0.0057593957402320385, -0.00027248233573270216, 0.003885350572544307, 0.00022362281488693097, 0.006609741178005571, -0.003292722856107429, -0.005873218251875897},
  {-0.0032648855847982987, -0.0005945107515129156, -0.009885693540838514, -0.007803477134405362, 0.0538403407841302, -0.007446654755103316, -0.0025216534232170153, 0.004499568241334517, 0.009591034277125177, 0.00008612746932654223, 0.003386296829505543, -0.002600737873367083, 0.000621148057330571, -0.006603857049454523, -0.009241221870695395},
  {-0.0034218261781145264, 0.004858466255616657, -0.004832120597294347, -0.0011466944120756015, -0.007446654755103318, 0.043746559133865104, 0.008962713024625965, -0.011099652042761613, -0.0006620240117921668, -0.0012591530037708058, -0.006899982952117269, 0.0019732354732442878, -0.002445676747004324, -0.006454778807421816, 0.0033303577606412765},
  {-0.00371734018385458, -0.011667819871670206, -0.004576186966267273, 0.009925736654597632, -0.0025216534232170153, 0.008962713024625965, 0.03664582756831382, -0.009470328827284009, -0.006213741694945105, 0.007118775954484294, -0.0006741237990418526, -0.006003374957986355, 0.005718636997353189, -0.0005191095254772077, -0.008466566781233205},
  {-0.007694897715679857, 0.0031697801901690347, 0.002510721148395566, -0.0007664415942207059, 0.004499568241334515, -0.011099652042761617, -0.009470328827284016, 0.057734267068088, 0.005521731225009532, -0.017008048805405164, 0.006749693090695894, -0.006348460110898, -0.007879244727681924, -0.005321753837620446, 0.011126783289057604},
  {0.005260905282822458, 0.0067613450046548505, -0.010126911913571181, -0.00575939574023204, 0.009591034277125177, -0.0006620240117921668, -0.006213741694945106, 0.005521731225009532, 0.04610670018969681, -0.010427010812879566, -0.0009861561285861987, -0.008896020395949732, -0.0037627528719902485, 0.00033704453138913093, -0.003173552163182467},
  {0.0013607957548231744, -0.0037599761532668475, 0.01153595152487264, -0.0002724823357326985, 0.0000861274693265406, -0.0012591530037708062, 0.007118775954484294, -0.01700804880540517, -0.010427010812879568, 0.05909125052583998, 0.002192545816395299, -0.002057672237277737, -0.004801518314458135, -0.014065326026662672, -0.005619012077114913},
  {0.0019707858957027763, 0.005571796842520162, 0.005773054728678472, 0.003885350572544309, 0.003386296829505542, -0.006899982952117272, -0.0006741237990418522, 0.006749693090695893, -0.0009861561285862005, 0.0021925458163952988, 0.024417715762416557, -0.003037163447600162, -0.011173674374382736, -0.0008193127407211239, -0.007137012700864866},
  {0.006708452591621083, -0.0071098291510566895, 0.005286558230422046, 0.00022362281488693216, -0.0026007378733670806, 0.0019732354732442886, -0.006003374957986352, -0.006348460110897999, -0.008896020395949732, -0.002057672237277737, -0.003037163447600163, 0.04762367868805726, 0.0008818947598625008, -0.0007262691465810616, -0.006482422704208912},
  {-0.005107684668720825, -0.0044777735406402895, -0.005543879869413772, 0.006609741178005571, 0.0006211480573305693, -0.002445676747004324, 0.0057186369973531905, -0.00787924472768192, -0.003762752871990247, -0.004801518314458137, -0.011173674374382736, 0.0008818947598624995, 0.042639958466440225, 0.0010194948614718209, 0.0033872675386130637},
  {0.004402554308030674, -0.011250694688474675, 0.004477221036185477, -0.003292722856107429, -0.006603857049454523, -0.006454778807421815, -0.0005191095254772072, -0.005321753837620446, 0.0003370445313891318, -0.014065326026662679, -0.0008193127407211239, -0.0007262691465810616, 0.0010194948614718226, 0.05244900188599414, -0.000256550861960499},
  {-0.00334987648531921, 0.007465228985669282, -0.006204169580739178, -0.005873218251875899, -0.009241221870695395, 0.003330357760641278, -0.008466566781233205, 0.011126783289057604, -0.0031735521631824654, -0.005619012077114915, -0.007137012700864866, -0.006482422704208912, 0.0033872675386130632, -0.000256550861960499, 0.05380987317762257}};



/* ============ Likelihood computations: ========== */


static void extractDimensionlessVariableVector(LALInferenceVariables *currentParams, REAL8 *x, INT4 mode) {

  REAL8 mean[15];
  UINT4 i=0;
  memset(x, 0, 15*sizeof(REAL8));
  memset(mean, 0, 15*sizeof(REAL8));

  if (mode==0) { /* Univariate */
    for(i=0;i<15;i++) mean[i]=LALInferenceAnalyticMeansCBC[i];
  } else if (mode==1) {
    for(i=0;i<15;i++) mean[i]=LALInferenceAnalyticMeansCBC[i] - 4.0/scaling[i]*sqrt(CM[i][i]);
  } else if (mode==2) {
    for(i=0;i<15;i++) mean[i]=LALInferenceAnalyticMeansCBC[i] + 4.0/scaling[i]*sqrt(CM[i][i]);
  } else {
    printf("Error!  Unrecognized mode in analytic likelihood!\n");
    exit(1);
  }
  
  for(i=0;i<15;i++)
  {
    x[i]=scaling[i]* (LALInferenceGetREAL8Variable(currentParams,LALInferenceAnalyticNamesCBC[i]) - mean[i]);
  }
  
}

static void extractBurstDimensionlessVariableVector(LALInferenceVariables *currentParams, REAL8 *x, INT4 mode) {
  REAL8 frequency=100.,q=1.0,loghrss, psi, ra, dec, t,alpha=0.0,polar_eccentricity=0.0;
  const UINT4 nparams=9;
  REAL8 mean[nparams];
  
  memset(x, 0, nparams*sizeof(REAL8));
  memset(mean, 0, nparams*sizeof(REAL8));

  if (mode==0) {
    mean[0] = 211.0; // freq
    mean[1] = 6.0;  //quality
    mean[2] = -46.; //loghrss
    mean[3] = 0.001; //time
    mean[4] = M_PI; //ra
    mean[5] = 0.0;   //dec
    mean[6] = 0.7; //psi
    mean[7] =0.5; // alpha (polar_angle)
    mean[8] =0.25;    // polar_eccentricity 
  } else if (mode==1) {
    mean[0] = 211.0;
    mean[1] = 6.0;
    mean[2] = -46.;
    mean[3] = 0.001;
    mean[4] = M_PI;
    mean[5] = 0.0;
    mean[6] = 0.7;
    mean[7] =0.5;
    mean[8]= 0.25;
  } else if (mode==2) {
    /* set means of second mode to be 8 sigma from first mode */
    mean[0] = 211 + 8./burst_scaling[0]*sqrt(bCM[0][0]);
    mean[1] = 6.0 + 8./burst_scaling[1]*sqrt(bCM[1][1]);
    mean[2] = -46. + 8./burst_scaling[2]*sqrt(bCM[2][2]);
    mean[3] = 0.001 + 8./burst_scaling[3]*sqrt(bCM[3][3]);
    mean[4] = M_PI + 8./burst_scaling[4]*sqrt(bCM[4][4]);
    mean[5] = 0.0 + 8./burst_scaling[5]*sqrt(bCM[5][5]);
    mean[6] = 0.7 + 8./burst_scaling[7]*sqrt(bCM[6][6]);
    mean[7] = 0.5 + 8./burst_scaling[7]*sqrt(bCM[7][7]);
    mean[8] = 0.25 + 8./burst_scaling[8]*sqrt(bCM[8][8]);
  } else {
    printf("Error!  Unrecognized mode in analytic likelihood!\n");
    exit(1);
  }
 
  loghrss = LALInferenceGetREAL8Variable(currentParams, "loghrss");
  frequency = LALInferenceGetREAL8Variable(currentParams, "frequency");
  q = LALInferenceGetREAL8Variable(currentParams, "quality");
  psi = LALInferenceGetREAL8Variable(currentParams, "polarisation");
  alpha = LALInferenceGetREAL8Variable(currentParams, "alpha");
  ra = LALInferenceGetREAL8Variable(currentParams, "rightascension");
  dec = LALInferenceGetREAL8Variable(currentParams, "declination");
  t = LALInferenceGetREAL8Variable(currentParams, "time");
  polar_eccentricity= LALInferenceGetREAL8Variable(currentParams, "polar_eccentricity");
  x[0] = burst_scaling[0] * (frequency    - mean[0]);
  x[1] = burst_scaling[1] * (q   - mean[1]);
  x[2] = burst_scaling[2] * (loghrss  - mean[2]);
  x[3] = burst_scaling[3] * (t   - mean[3]);
  x[4] = burst_scaling[4] * (ra   - mean[4]);
  x[5] = burst_scaling[5] * (dec    - mean[5]);
  x[6] = burst_scaling[6] * (psi     - mean[6]);
  x[7] = burst_scaling[7] * (alpha     - mean[7]);
  x[8] = burst_scaling[8] * (polar_eccentricity - mean[8]);  
}

REAL8 LALInferenceCorrelatedAnalyticLogLikelihood(LALInferenceVariables *currentParams, 
                                                  LALInferenceIFOData *data, 
                                                  LALInferenceModel *model) {
  INT4 tmpdim=0;
  int cbc=1;
  const REAL8 *cm=NULL;
  int ifo=0;
  while(data){
    model->ifo_loglikelihoods[ifo]=0.0;
    ifo++;
    data=data->next;
  }

  if (LALInferenceCheckVariable(currentParams, "logdistance")){
    /* We are dealing with spinning CBC. Set dimensions and CVM accordingly*/
    tmpdim = 15;
    cm=&(CM[0][0]);
  }
  else if ( LALInferenceCheckVariable(currentParams, "hrss") ||LALInferenceCheckVariable(currentParams, "loghrss")){
    /* We are dealing with a  burst. Set dimensions and CVM accordinly*/
    tmpdim = 9;
    cm=&(bCM[0][0]);
    cbc=0;
  }
  const INT4 DIM = tmpdim;
  gsl_matrix *LUCM = NULL;
  gsl_permutation *LUCMPerm = NULL;
  INT4 mode = 0;

  REAL8 x[DIM];
  REAL8 xOrig[DIM];

  gsl_vector_view xView = gsl_vector_view_array(x, DIM);

  if (cbc==1)
    extractDimensionlessVariableVector(currentParams, x, mode);
  else
    extractBurstDimensionlessVariableVector(currentParams, x, mode);

  memcpy(xOrig, x, DIM*sizeof(REAL8));

  if (LUCM==NULL) {
    gsl_matrix_const_view CMView = gsl_matrix_const_view_array(cm, DIM, DIM);
    int signum;

    LUCM = gsl_matrix_alloc(DIM, DIM);
    LUCMPerm = gsl_permutation_alloc(DIM);

    gsl_matrix_memcpy(LUCM, &(CMView.matrix));

    gsl_linalg_LU_decomp(LUCM, LUCMPerm, &signum);
  }

  gsl_linalg_LU_svx(LUCM, LUCMPerm, &(xView.vector));

  INT4 i;
  REAL8 sum = 0.0;
  for (i = 0; i < DIM; i++) {
    sum += xOrig[i]*x[i];
  }
  if(LUCM) gsl_matrix_free(LUCM);
  if(LUCMPerm) gsl_permutation_free(LUCMPerm);
  return -sum/2.0;
}

REAL8 LALInferenceBimodalCorrelatedAnalyticLogLikelihood(LALInferenceVariables *currentParams,
                                                  LALInferenceIFOData *data,
                                                LALInferenceModel *model) {
  INT4 tmpdim=0;
  int cbc=1;
  const REAL8 *cm=NULL;
  int ifo=0;
  while(data){
    model->ifo_loglikelihoods[ifo]=0.0;
    ifo++;
    data=data->next;
  }

  if (LALInferenceCheckVariable(currentParams, "logdistance")){
    /* We are dealing with spinning CBC. Set dimensions and CVM accordingly*/
    tmpdim = 15;
    cm=&(CM[0][0]);
  }
  else if ( LALInferenceCheckVariable(currentParams, "hrss") ||LALInferenceCheckVariable(currentParams, "loghrss")){
     /* We are dealing with a burst. Set dimensions and CVM accordingly*/
    tmpdim = 9;
    cm=&(bCM[0][0]);
    cbc=0;
  }
  const INT4 MODES = 2;
  const INT4 DIM = tmpdim;
  INT4 i, mode;
  REAL8 sum = 0.0;
  REAL8 a, b;
  gsl_matrix *LUCM = NULL;
  gsl_permutation *LUCMPerm = NULL;
  gsl_vector_view xView;

  REAL8 x[DIM];
  REAL8 xOrig[DIM];
  REAL8 exps[MODES];

  if (LUCM==NULL) {
    gsl_matrix_const_view CMView = gsl_matrix_const_view_array(cm, DIM, DIM);
    int signum;

    LUCM = gsl_matrix_alloc(DIM, DIM);
    LUCMPerm = gsl_permutation_alloc(DIM);

    gsl_matrix_memcpy(LUCM, &(CMView.matrix));

    gsl_linalg_LU_decomp(LUCM, LUCMPerm, &signum);
  }

  for(mode = 1; mode < 3; mode++) {
    xView = gsl_vector_view_array(x, DIM);
    if (cbc==1)
      extractDimensionlessVariableVector(currentParams, x, mode);
    else
      extractBurstDimensionlessVariableVector(currentParams, x, mode);

    memcpy(xOrig, x, DIM*sizeof(REAL8));

    gsl_linalg_LU_svx(LUCM, LUCMPerm, &(xView.vector));

    sum = 0.0;
    for (i = 0; i < DIM; i++) {
      sum += xOrig[i]*x[i];
    }
    exps[mode-1] = -sum/2.0;
  }

  /* Assumes only two modes used from here on out */
  if (exps[0] > exps[1]) {
    a = exps[0];
    b = exps[1];
  } else {
    a = exps[1];
    b = exps[0];
  }
  /* attempt to keep returned values finite */
  if(LUCM) gsl_matrix_free(LUCM);
  if(LUCMPerm) gsl_permutation_free(LUCMPerm);
  return a + log1p(exp(b-a));
}

REAL8 LALInferenceRosenbrockLogLikelihood(LALInferenceVariables *currentParams,
                                          LALInferenceIFOData *data,
                                          LALInferenceModel *model) {
  const INT4 DIM = 15;
  REAL8 x[DIM];

  REAL8 sum = 0;
  INT4 mode = 0;
  INT4 i;
  int ifo=0;
  while(data){
    model->ifo_loglikelihoods[ifo]=0.0;
    ifo++;
    data=data->next;
  }
  extractDimensionlessVariableVector(currentParams, x, mode);

  for (i = 0; i < DIM; i++) x[i] += 1.0;

  for (i = 0; i < DIM-1; i++) {
    REAL8 oneMX = 1.0 - x[i];
    REAL8 dx = x[i+1] - x[i]*x[i];

    sum += oneMX*oneMX + 100.0*dx*dx;
  }

  return -sum;
}

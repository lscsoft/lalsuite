/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpPTFWaveform.c
 *
 * Author: Brown, D. A., and Fazi, D.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0 
<lalVerbatim file="FindChirpPTFWaveformCV">
Author: Brown, D. A., and Fazi, D.
$Id$
</lalVerbatim> 

<lalLaTeX>
\subsection{Module \texttt{FindChirpPTFWaveform.c}}
\label{ss:FindChirpPTFDerviatives.c}

Provides functions to create physical template family templates in a
form that can be used by the \texttt{FindChirpPTFFilter()} function.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{FindChirpPTFTemplateCP}
\idx{LALFindChirpPTFWaveform()}

The function \texttt{LALFindChirpPTFWaveform()} creates vectors containing the
evolution of the dynamical variables needed by the physical template family
template as described by the algorithm below.

\subsubsection*{Algorithm}

Blah.

\subsubsection*{Uses}
\begin{verbatim}
LALCalloc()
LALFree()
LALCreateVector()
LALDestroyVector()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{FindChirpPTFWaveformCV}}
</lalLaTeX> 
#endif

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>

NRCSID(FINDCHIRPPTFWAVEFORMC, "$Id$");

/* define a structure so that the ptf waveform parameters */ 
/* can be used by the GSL ODE integrator                  */
typedef struct
{
  /* input parameters which control evolution */
  REAL4Vector* orbital;       /* pn param of orbital evolution                */
  REAL4Vector* spin;          /* pn params of spin compt of orbital evolution */
  float        S1_spin_orbit; /* coeff of S-O term in first spin evolution    */
  float        LNhat;         /* coeff which multiples dS_dt to get dLNhat_dt */
  float        mag_S1;        /* magnitude of the spin = chi1 * m1 *m1        */
  /* output parameters used to monitor evolution */
  double       LNhat_dot_S1;  /* dot product of LNhat and S1                  */
}
ptf_evolution_params_t;


/* function that computes the derivatives of the dynamical */
/* variables for the GSL ODE integrator                    */
int ptf_waveform_derivatives( 
    double t, const double y[], double dydt[], void* params )
{
  /* equation numbers in description of variables and algorithms refer to
   * Pan, Buonanno, Chan and Vallisneri, Phys. Rev. D 69, 104017 (2004)
   */

  /* In this code the first body is the larger, spinning body with spin
   * magnitude given by |S1| = chi1 * m1^2 with 0 \le chi1 \le 1          
   */

  /* post newtonian coeffients which are independent of time */
  ptf_evolution_params_t* pn_params = (ptf_evolution_params_t*) params;

  /* evolution variables */
  /* y[0] stores Phi, the gravitational wave phase as in Eq. (15), but it   */
  /* not needed in the evolution equations for the precession convention    */
  const double omega  = y[1]; /* omega as in Eq. (5)                        */
  const double S1x    = y[2]; /* x-cmpt of first bodies spin in Eq. (6)     */
  const double S1y    = y[3]; /* y-cmpt of first bodies spin in Eq. (6)     */
  const double S1z    = y[4]; /* z-cmpt of first bodies spin in Eq. (6)     */
  const double LNhatx = y[5]; /* x-cmpt of orb plane normal in Eq. (7)      */
  const double LNhaty = y[6]; /* y-cmpt of orb plane normal in Eq. (7)      */
  const double LNhatz = y[7]; /* z-cmpt of orb plane normal in Eq. (7)      */
  const double e1x    = y[8]; /* x-cmpt of orb plane basis vector 1 Eq.(13) */
  const double e1y    = y[9]; /* x-cmpt of orb plane basis vector 1 Eq.(13) */
  const double e1z   = y[10]; /* x-cmpt of orb plane basis vector 1 Eq.(13) */
  const double e2x   = y[11]; /* x-cmpt of orb plane basis vector 1 Eq.(13) */
  const double e2y   = y[12]; /* x-cmpt of orb plane basis vector 1 Eq.(13) */
  const double e2z   = y[13]; /* x-cmpt of orb plane basis vector 1 Eq.(13) */

  /* powers of omega used in post-Newtonian expansion */
  const double omega_1_3 = pow( omega, 1.0/3.0 );
  const double omega_2_3  = omega_1_3 * omega_1_3;
  const double omega_4_3  = omega * omega_1_3;
  const double omega_5_3  = omega_4_3 * omega_1_3;
  const double omega_7_3  = omega_5_3 * omega_2_3;
  const double omega_6_3  = omega * omega;
  const double omega_11_3 = omega_7_3 * omega_4_3;

  /* coefficients of the cross products in Eqs. (6) and (7) */
  const double S1_spin_orbit_coeff = omega_5_3 * pn_params->S1_spin_orbit;
  const double LNhat_coeff         = omega_1_3 * pn_params->LNhat;

  /* compute the cross product of LNhat and S1 */
  const double LNhat_cross_S1_x = LNhaty * S1z - LNhatz * S1y;
  const double LNhat_cross_S1_y = LNhatz * S1x - LNhatx * S1z;
  const double LNhat_cross_S1_z = LNhatx * S1y - LNhaty * S1x;

  /* dot product of LNhat and S1: pass back so we can check it is constant */
  const double LNhat_dot_S1 = LNhatx * S1x + LNhaty * S1y + LNhatz * S1z;
  pn_params->LNhat_dot_S1 = LNhat_dot_S1;

  /* OmegaL as defined in Eq. (7) */
  const double OmegaLx = - LNhat_coeff * S1_spin_orbit_coeff * S1x;
  const double OmegaLy = - LNhat_coeff * S1_spin_orbit_coeff * S1y;
  const double OmegaLz = - LNhat_coeff * S1_spin_orbit_coeff * S1z;

  /* dot product of OmegaL and LNhat needed in Eq. (14) */
  const double OmegaL_dot_LNhat = 
    OmegaLx * LNhatx + OmegaLy * LNhaty + OmegaLz * LNhatz;

  /* Omegae as defined by Eq. (14) */
  const double Omegaex = OmegaLx - OmegaL_dot_LNhat * LNhatx;
  const double Omegaey = OmegaLy - OmegaL_dot_LNhat * LNhaty;
  const double Omegaez = OmegaLz - OmegaL_dot_LNhat * LNhatz;

  /* compute the derivatives of the spin precession given by Eq. (6) */
  const double dS1x_dt = S1_spin_orbit_coeff * LNhat_cross_S1_x;
  const double dS1y_dt = S1_spin_orbit_coeff * LNhat_cross_S1_y;
  const double dS1z_dt = S1_spin_orbit_coeff * LNhat_cross_S1_z;

  /* compute the derivatives of the orbital precession given by Eq. (7) */
  const double dLNhatx_dt = LNhat_coeff * dS1x_dt;
  const double dLNhaty_dt = LNhat_coeff * dS1y_dt;
  const double dLNhatz_dt = LNhat_coeff * dS1z_dt;

  /* compute the of derivatives of the orbital plane basis given by Eq.(13) */
  const double de1x_dt = Omegaey * e1z - Omegaez * e1y;
  const double de1y_dt = Omegaez * e1x - Omegaex * e1z;
  const double de1z_dt = Omegaex * e1y - Omegaey * e1x;
  const double de2x_dt = Omegaey * e2z - Omegaez * e2y;
  const double de2y_dt = Omegaez * e2x - Omegaex * e2z;
  const double de2z_dt = Omegaex * e2y - Omegaey * e2x;

   /* compute the derivative of the orbital phase given by Eq. (5) */
  const double domega_dt = omega_11_3 * (
      /* contribution due to purely orbital evolution */
      + pn_params->orbital->data[0]                            /*   0 */
      + pn_params->orbital->data[1] * omega_1_3                /* 0.5 */
      + pn_params->orbital->data[2] * omega_2_3                /* 1.0 */
      + pn_params->orbital->data[3] * omega                    /* 1.5 */
      + pn_params->orbital->data[4] * omega_4_3                /* 2.0 */
      + pn_params->orbital->data[5] * omega_5_3                /* 2.5 */
      + pn_params->orbital->data[6] * omega_6_3                /* 3.0 */
      + pn_params->orbital->data[7] * omega_6_3 * log( omega ) /* 3.0 */
      + pn_params->orbital->data[8] * omega_7_3                /* 3.5 */
      /* contribution at 1.5 pN due spin-orbit interaction */
      + (pn_params->orbital->data[0]
        * pn_params->spin->data[0] * LNhat_dot_S1) * omega
      );

  /* compute the derivative of the gravitational wave phase the gw */
  /* phase evolution is purely orbital as we are working in the    */
  /* precessing converion: see the discussion after Eq. (14)       */
  const double dPhi_dt = omega;
     
  /* copy derivatives into output array */
  dydt[0]  = dPhi_dt;
  dydt[1]  = domega_dt;
  dydt[2]  = dS1x_dt;
  dydt[3]  = dS1y_dt;
  dydt[4]  = dS1z_dt;
  dydt[5]  = dLNhatx_dt;
  dydt[6]  = dLNhaty_dt;
  dydt[7]  = dLNhatz_dt;
  dydt[8]  = de1x_dt;
  dydt[9]  = de1y_dt;
  dydt[10] = de1z_dt;
  dydt[11] = de2x_dt;
  dydt[12] = de2y_dt;
  dydt[13] = de2z_dt;

  return GSL_SUCCESS;
}


REAL4Vector*
XLALPTFOmegaPNCoeffsOrbital( REAL4 m1, REAL4 m2 )
{
  static const char* func = "XLALPTFOmegaPNCoeffsOrbital";
  float m_total = m1 + m2;
  float eta = (m1 * m2) / (m_total * m_total);
  const UINT4 max_pn_order = 9;
  REAL4Vector* c_vec = XLALCreateREAL4Vector( max_pn_order );
  REAL4* c;

  if ( ! c_vec )
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );
  else
    c = c_vec->data;

  /* 0 pN coefficient which multiplies all the other terms */
  c[0] = (96.0/5.0) * eta;

  /* 0.5 pN correction */
  c[1] = 0;

  /* 1 pN correction */
  c[2] = c[0] * (-1.0/336.0) * (743.0 + 924.0 * eta);

  /* 1.5 pN correction */
  c[3] = c[0] * 4.0 * LAL_PI;

  /* 2 pN correction */
  c[4] = c[0] * (34103.0 + 122949.0 * eta + 59472.0 * eta * eta) 
    / 18144.0;

  /* 2.5 pN correction */
  c[5] = c[0] * (-1.0/672.0) * (4159.0 + 15876.0 * eta) * LAL_PI;

  /* 3 pN correction excluding log((M omega)^2/3) term */
  c[6] = c[0] * ( 16447322263.0/139708800.0 
      - 1712.0/105.0 * LAL_GAMMA
      + 16.0/3.0 * LAL_PI * LAL_PI
      + ( -273811877.0/1088640.0 + 451.0/48.0 * LAL_PI * LAL_PI 
        - (88.0 * 1039.0)/(3.0 * 4620.0) ) * eta
      + 541.0/896.0 * eta * eta
      - 5605.0/2592.0 * eta * eta * eta 
      - 856.0/105.0 * log(16.0) );

  /* 3 pN correction of log((M omega)^2/3) term */
  c[7] = c[0] * (-1712.0/315.0);

  /* 3.5 pN correction */
  c[8] = c[0] * (-13245.0 + 717350.0 * eta + 731960.0 * eta * eta) 
    * (LAL_PI/12096.0);

  return c_vec;
}


REAL4Vector*
XLALPTFOmegaPNCoeffsSpin( REAL4 m1, REAL4 m2, 
    REAL4 chi1, REAL4 chi2, 
    REAL4 Q1, REAL4 Q2 )
{
  static const char* func = "XLALPTFOmegaPNCoeffsSpin";
  const float m_total = m1 + m2;
  const float m1_5 = m1 * m1 * m1 * m1 * m1;
  const float m2_5 = m2 * m2 * m2 * m2 * m2;
  const UINT4 max_spin_order = 6;
  REAL4Vector* c_vec = XLALCreateREAL4Vector( max_spin_order );
  REAL4* c;

  if ( ! c_vec )
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );
  else
    c = c_vec->data;

  /* 1.5 pN spin-orbit interaction from body 1 */
  c[0] = (-1.0 / (12.0 * m_total * m_total)) * ( 113.0 + 75.0 * (m2/m1) );

  /* 1.5 pN spin-orbit interaction from body 2 */
  c[1] = (-1.0 / (12.0 * m_total * m_total)) * ( 113.0 + 75.0 * (m1/m2) );

  /* 2 pN spin-spin interaction: coefficent of S_1 . S_2 */
  c[2] = -247.0 / (48.0 * m_total * m_total * m1 * m2);

  /* 2 pN spin-spin interaction: coefficent of (L . S_1)(L . S_2) */
  c[3] = 721.0 / (48.0 * m_total * m_total * m1 * m2);

  /* 2 pN Kerr quadrupole-monopole coupling for first body */
  if ( chi1 )
    c[4] = ( -5.0 * Q1 ) 
      / ( 2.0 * m1_5 * chi1 * chi1 * m_total * m_total );
  else
    c[4] = 0;

  /* 2 pN Kerr quadrupole-monopole coupling for second body */
  if ( chi2 )
    c[5] = ( -5.0 * Q2 ) 
      / ( 2.0 * m2_5 * chi2 * chi2 * m_total * m_total );
  else
    c[5] = 0;

  return c_vec;
}


REAL4Vector*
XLALPTFOmegaPNCoeffsEnergy( REAL4 m1, REAL4 m2,
    REAL4 chi1, REAL4 chi2,
    REAL4 Q1, REAL4 Q2 )
{
  /* These coefficients are derived from Eqs. (11) and (12) and (13) of */
  /* Buonanno, Chen and Vallisneri, Phys. Rev. D 67, 104025 (BCV2)      */
  static const char* func = "XLALPTFOmegaPNCoeffsEnergy";
  const float m_total = m1 + m2;
  const float mu = m1 * m2 / m_total;
  const float eta = mu / m_total;
  const float m1_5 = m1 * m1 * m1 * m1 * m1;
  const float m2_5 = m2 * m2 * m2 * m2 * m2;
  const UINT4 max_pn_order = 9;
  REAL4Vector* c_vec = XLALCreateREAL4Vector( max_pn_order );
  REAL4* c;

  if ( ! c_vec )
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );
  else
    c = c_vec->data;

  /* 0 pN coefficient which multiplies all the other terms */
  c[0] = -eta/3.0;

  /* 0.5 pN correction */
  c[1] = 0;

  /* 1 pN correction */
  c[2] = c[0] * (-1.0/6.0) * (9.0 + eta);

  /* 1.5 pN correction */
  c[3] = 0;

  /* 2 pN correction */
  c[4] = c[0] * (3.0/24.0) * (-81.0 + 57.0 * eta - eta * eta);

  /* 2.5 pN correction */
  c[5] = 0;

  /* 3 pN correction */
  c[6] = c[0] * 4.0 
    * ((-675.0/64.0) 
    + ( 34445.0/576.0 - 205*LAL_PI*LAL_PI/96.0 ) * eta
    - 155.0/96.0 * eta * eta
    - 35.0/5184.0 * eta * eta * eta);

  /* 2.0 pN quadrupole-monopole coupling term from first body */
  if ( chi1 )
    c[7] = c[0] * (3.0 * Q1) 
      / (2.0 * m1_5 * chi1 * chi1 * m_total * m_total );
  else
    c[7] = 0;

  /* 2.0 pN quadrupole-monopole coupling term from second body */
  if ( chi2 )
    c[8] = c[0] * (3.0 * Q2) 
      / (2.0 * m2_5 * chi2 * chi2 * m_total * m_total );
  else
    c[8] = 0;

  return c_vec;
}


static float 
spin_so_coeff( float ma, float mb )
{
  float m_total = ma + mb;
  float eta = (ma * mb) / (m_total * m_total);
  return (eta/2.0) * ( 4.0 + 3.0 * ma/mb );
}


static float 
orbital_coeff( float m1, float m2 )
{
  float m_total = m1 + m2;
  float eta = (m1 * m2) / (m_total * m_total);
  return -1.0 / (eta * m_total * m_total);
}

static float stpn_orbital_energy( double omega,
    double LNhat_dot_S1, double LNhat_dot_S2, double S1_dot_S2, 
    float m1, float m2, float chi1, float chi2, REAL4Vector* pn_params )
{
  /* Function which computes the derivative of the orbital energy with     */
  /* respect to omega. This is algorithm here implements the (analytically */
  /* computed) derivative of BCV2 Eq. (12) with respect to omega (i.e. it  */
  /* is the quantity on the left hand side of BCV2 Eq. (18).               */

  /* mass parameters */
  const double m_total = m1 + m2;

  /* magnitude of spin vectors needed for qm coupling term in energy */
  const double mag_S1 = m1 * m1 * chi1;
  const double mag_S2 = m2 * m2 * chi2;

  /* LNhat dot S_effective defined in BCV2 Eq. (7) */
  const double LNhat_dot_Seff
    = ( 1.0 + (3.0 * m2) / (4.0 * m1) ) * LNhat_dot_S1
    + ( 1.0 + (3.0 * m1) / (4.0 * m2) ) * LNhat_dot_S2;

  /* powers of omega used in post-Newtonian expansion */
  const double omega_1_3 = pow( omega, 1.0/3.0 );
  const double omega_minus1_3 = 1.0 / omega_1_3;
  const double omega_2_3  = omega_1_3 * omega_1_3;
  const double omega_4_3  = omega * omega_1_3;
  const double omega_5_3  = omega_4_3 * omega_1_3;

  /* derivative of orbital energy with respect to time */
  const double energy
    /* contribution from non-spinning pN terms */
    = pn_params->data[0] * omega_minus1_3            /*   0 */
    + pn_params->data[1]                             /* 0.5 */
    + pn_params->data[2] * omega_1_3                 /* 1.0 */
    + pn_params->data[3] * omega_2_3                 /* 1.5 */
    + pn_params->data[4] * omega                     /* 2.0 */
    + pn_params->data[5] * omega_4_3                 /* 2.5 */
    + pn_params->data[6] * omega_5_3                 /* 3.0 */
    /* contribution from spin orbit terms at 1.5 pN order */
    + pn_params->data[0] 
    * 20.0 * LNhat_dot_Seff * omega_2_3 / (3.0 * m_total * m_total)
    /* XXX contribution from spin spin terms at 2.0 pN order */
    + (-1.0/2.0) * ( S1_dot_S2 - 3.0 * LNhat_dot_S1 * LNhat_dot_S2 ) 
    * omega_5_3 / (m_total * m_total)
    /* contribution from first quadrupole-monopole coupling at 2.0 pN */
    + pn_params->data[7] * omega 
    * ( 3.0 * LNhat_dot_S1 * LNhat_dot_S1 - mag_S1 * mag_S1 )
    /* contribution from second quadrupole-monopole coupling at 2.0 pN */
    + pn_params->data[8] * omega 
    * ( 3.0 * LNhat_dot_S2 * LNhat_dot_S2 - mag_S2 * mag_S2 );

  return (float) energy;
}


/* Function for the evaluation of the time evolution of the phase, the */
/* frequency and the e_i basis vectors in the precessing convention    */
/* <lalVerbatim file="FindChirpPTFWaveformCP"> */
INT4
XLALFindChirpPTFWaveform( 
    REAL4Vector         *PTFphi,
    REAL4Vector         *PTFomega_2_3,
    REAL4VectorSequence *PTFe1,
    REAL4VectorSequence *PTFe2,
    InspiralTemplate    *tmplt,
    REAL8                deltaT
    )
/* </lalVerbatim> */
{
  static const char* func = "XLALFindChirpPTFWaveform";
  UINT4  i, len;
  UINT4  N = PTFphi->length;
  INT4   errcode = 0;
  REAL8  deltaF = 1.0 / ( (REAL8) N * deltaT ); 
  double f_min = tmplt->fLower;
  double m1 = tmplt->mass1;
  double m2 = tmplt->mass2;
  double chi1 = tmplt->chi;
  double kappa = tmplt->kappa;
  double t, t_next;                                                       
  double step_size;                                                       
  double dE_dt, dE_dt_n_1, dE_dt_n_2;
  double N_steps;

  ptf_evolution_params_t pn_params;
  const ptf_evolution_params_t* pn_params_ptr = &pn_params;

  const double m_total = m1 + m2;  
  const double geometrized_m_total = m_total * LAL_MTSUN_SI;  
  const double freq_step = geometrized_m_total * LAL_PI;
  const double step      = deltaT / geometrized_m_total;                      
  const double omegam_to_hz = 1.0 / freq_step;                            
  const int    num_evolution_variables = 14;
  double       y[num_evolution_variables], dydt[num_evolution_variables];             

  /* Dynamical evolution variables and their derivatives */
  double Phi ;   /* gravitational wave phase in BCV2 Eq. (18)     */
  double omega;  /* omega as in Eq. (5)                           */
  double S1x;    /* x-cmpt of first bodies spin in Eq. (6)        */
  double S1y;    /* y-cmpt of first bodies spin in Eq. (6)        */
  double S1z;    /* z-cmpt of first bodies spin in Eq. (6)        */
  double LNhatx; /* x-cmpt of orb plane normal in Eq. (7)         */
  double LNhaty; /* y-cmpt of orb plane normal in Eq. (7)         */
  double LNhatz; /* z-cmpt of orb plane normal in Eq. (7)         */
  double e1x;    /* x-component of the basis vector e1 in Eq.(13) */
  double e1y;    /* y-component of the basis vector e1 in Eq.(13) */
  double e1z;    /* z-component of the basis vector e2 in Eq.(13) */
  double e2x;    /* x-component of the basis vector e2 in Eq.(13) */
  double e2y;    /* y-component of the basis vector e2 in Eq.(13) */
  double e2z;    /* z-component of the basis vector e2 in Eq.(13) */

  /* create the differential equation solver */
  const gsl_odeiv_step_type* solver_type
    = gsl_odeiv_step_rkf45;
  gsl_odeiv_step* solver_step 
    = gsl_odeiv_step_alloc( solver_type, num_evolution_variables );
  gsl_odeiv_control* solver_control
    = gsl_odeiv_control_standard_new( 1.0e-5, 1.0e-5, 1.0, 1.0 );
  gsl_odeiv_evolve* solver_evolve 
    = gsl_odeiv_evolve_alloc( num_evolution_variables );
  gsl_odeiv_system solver_system = { ptf_waveform_derivatives,
    NULL, num_evolution_variables, (void*) pn_params_ptr };

  /* orbital energy pn parameters and pn constants vector */
  REAL4Vector* orbital_energy_coeffs;
  orbital_energy_coeffs = XLALPTFOmegaPNCoeffsEnergy( m1, m2, chi1, 0, 0, 0 );

  /* set up post-Newtonian coefficents needed in evolution */
  pn_params.orbital       = XLALPTFOmegaPNCoeffsOrbital( m1, m2 );
  pn_params.spin          = XLALPTFOmegaPNCoeffsSpin( m1, m2, chi1, 0, 0.0, 0 );
  pn_params.S1_spin_orbit = spin_so_coeff( m2, m1 ); /* need m2/m1 */
  pn_params.LNhat         = orbital_coeff( m1, m2 );
  pn_params.mag_S1        = m1 * m1 * chi1;

  /* set up initial values for dynamical variables */
  /* in the precessing convention                  */

  Phi    = y[0]  = 0.0;
  omega  = y[1]  = f_min / omegam_to_hz;
  S1x    = y[2]  = sqrt( 1 - kappa * kappa ) * pn_params.mag_S1 ;
  S1y    = y[3]  = 0;
  S1z    = y[4]  = kappa * pn_params.mag_S1 ;  
  LNhatx = y[5]  = 0.0;
  LNhaty = y[6]  = 0.0;
  LNhatz = y[7]  = 1.0;
  e1x    = y[8]  = 1.0;
  e1y    = y[9]  = 0.0;
  e1z    = y[10] = 0.0;
  e2x    = y[11] = 0.0;
  e2y    = y[12] = 1.0;
  e2z    = y[13] = 0.0;

  /* start computing the waveform we use a while() loop which we   */
  /* break out of when one of three possible temination conditions */
  /* or two error conditions is reached                            */
  
  /* Zero out the dynamical variables so they don't contain garbage */
  memset ( PTFomega_2_3->data, 0, N * sizeof(REAL4));
  memset ( PTFphi->data, 0, N * sizeof(REAL4));
  memset ( PTFe1->data, 0, 3 * N * sizeof(REAL4));
  memset ( PTFe2->data, 0, 3 * N * sizeof(REAL4));
  
  i = 0; 
  t = 0;
  while ( 1 )
  {
    /* if we have run out of memory for the waveform, break out of the loop */
    if ( i >= N  )
    {
      XLALPrintError( "XLAL Error: output too short for PTF waveform" );
      errcode = XLAL_ENOMEM;
      break;
    }

    /* compute the gravitational waveform from the dynamical */
    /* variables and store it in the output structures       */
    PTFomega_2_3->data[i]    = (float) (pow( omega, 2.0/3.0 ));
    PTFphi->data[i]          = (float) (Phi);
    PTFe1->data[i]           = (float) (e1x);
    PTFe1->data[N + i]       = (float) (e1y);
    PTFe1->data[2 * N + i]   = (float) (e1z);
    PTFe2->data[i]           = (float) (e2x);
    PTFe2->data[N + i]       = (float) (e2y);
    PTFe2->data[2 * N + i]   = (float) (e2z);     

    /* advance the time (which is in units of total mass) */
    t = i * step;
    t_next = ++i * step;
    step_size = step;

    /* call the solver to get the next timestep */
    errcode = gsl_odeiv_evolve_apply( 
        solver_evolve, solver_control, solver_step, &solver_system,
        &t, t_next, &step_size, y );

    /* check for that the solver exited successfully */
    if ( errcode != GSL_SUCCESS )
    {
      XLALPrintError( "XLAL Error: GSL ODE integrator failure" );
      errcode = XLAL_EFAILED;
      break;
    }

    /* copy the output variables necessary to compute the gw */
    Phi    = y[0];
    omega  = y[1];
    LNhatx = y[5];
    LNhaty = y[6];
    LNhatz = y[7];
    e1x    = y[8];
    e1y    = y[9];
    e1z    = y[10];
    e2x    = y[11];
    e2y    = y[12];
    e2z    = y[13];

    /* exit with an error if any of the dynamical variables contain NaN */
    if ( isnan( Phi ) || isnan( omega ) || isnan( LNhatx ) || 
        isnan( LNhaty ) || isnan( LNhatz ) || isnan( e1x ) || isnan( e1y ) ||
        isnan( e1z ) || isnan( e2x ) || isnan( e2y ) || isnan( e2z ) )
    {
      /* check if we are close to the MECO */
      N_steps = ((i-2) * dE_dt_n_1 - (i-1) * dE_dt_n_2) / 
        ( dE_dt_n_1 -  dE_dt_n_2) - i + 1;

      if ( m2 > ( 0.5 * m1 - 1.4) && chi1 > (0.7 - kappa) )
      {
        if ( N_steps <= 5.0 && N_steps > 0.0) 
        {
          /* fprintf(stderr,"Secondary condition on MECO reached\n"); */
          break;
        }  
        else
        {  
          fprintf(stderr,"cycle %d\n",i);
          XLALPrintError( "XLAL Error: NaN in PTF dynamical variables\n" );
          errcode = XLAL_EFAILED;    
          break;
        }   

      }
      else  
      {  
        if ( N_steps <= 3.0 && N_steps > 0.0) 
        {
          /* secondary condition on MECO reached */
          break;
        }  
        else
        {   
          fprintf(stderr,"cycle %d\n",i);
          XLALPrintError( "XLAL Error: NaN in PTF dynamical variables\n" );
          errcode = XLAL_EFAILED;    
          break;
        }   
      }
    }

    /*  Store the last two values of dE/dt so as to be able to estimate   */
    /* how far from the MECO condition we are in case the code is failing */
    if ( i <= 1 )
    {  
      dE_dt_n_1 = stpn_orbital_energy( omega, pn_params.LNhat_dot_S1,
          0, 0, m1, m2, chi1, 0, orbital_energy_coeffs);
      dE_dt_n_2 = dE_dt_n_1 * 1.01;
    }  
    else if ( i > 1 )
    {  
      dE_dt_n_2 = dE_dt_n_1;
      dE_dt_n_1 = stpn_orbital_energy( omega, pn_params.LNhat_dot_S1,
          0, 0, m1, m2, chi1, 0, orbital_energy_coeffs);
    } 

    /* terminate if domega_dt is no longer positive as this means that */
    /* the adiabatic approximation has probably broken down            */
    ptf_waveform_derivatives( t, y, dydt, (void*) pn_params_ptr );
    if ( dydt[1] <= 0 ) break;

    /* terminate if the derivative of the orbital energy is zero or positive */
    /* this is the MECO termination condition discessed in BCV2 Eq. (13)     */
    if ( (dE_dt = stpn_orbital_energy( omega, pn_params.LNhat_dot_S1,
            0, 0, m1, m2, chi1, 0, orbital_energy_coeffs )) >= 0 ) break;

    /* If all check are ok set the final frequency */
    tmplt->fFinal = omega * omegam_to_hz;
    
  
  } /* end of evolution while ( 1 ) loop */
    
  /* free the memory used by the ode solver */
  gsl_odeiv_evolve_free( solver_evolve );
  gsl_odeiv_control_free( solver_control );
  gsl_odeiv_step_free( solver_step );

  /* free the memory used for the pn coefficients */
  XLALDestroyREAL4Vector( orbital_energy_coeffs );
  XLALDestroyREAL4Vector( pn_params.orbital );
  XLALDestroyREAL4Vector( pn_params.spin);

  /* Set the length of the template */
  tmplt->tC     = deltaT * (REAL8) i;
  
  /* shift the waveform so that the coalescence time */
  /* corresponds to the end of the segment           */
  len = N - i;
  
  /* Move the waveform at the end of the segment  */
  memmove( PTFomega_2_3->data + len, PTFomega_2_3->data, i * sizeof(REAL4) );
  memmove( PTFphi->data + len, PTFphi->data, i * sizeof(REAL4) );
  memmove( PTFe1->data + len, PTFe1->data, (2 * N + i) * sizeof(REAL4) );
  memmove( PTFe2->data + len, PTFe2->data, (2 * N + i) * sizeof(REAL4) );

  /* Set the waveform to zero at the beginning of the segment */
  memset ( PTFomega_2_3->data, 0, len * sizeof(REAL4));
  memset ( PTFphi->data, 0, len * sizeof(REAL4));
  memset ( PTFe1->data, 0, len * sizeof(REAL4));
  memset ( PTFe2->data, 0, len * sizeof(REAL4));

  /* the GSL success code is probably the same as LAL's but just in case... */
  if ( errcode == GSL_SUCCESS ) 
    errcode = XLAL_SUCCESS;
  else
    XLAL_ERROR( func, errcode );

  return errcode;
}

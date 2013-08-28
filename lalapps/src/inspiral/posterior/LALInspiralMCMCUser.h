/**
 * \author Dietz, A. & Veitch, J.
 * \file
 * \ingroup inspiral
 *
 * \brief Header file for the MCMC user code.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/LALInspiralMCMCUser.h>
 * \endcode
 *
 * This header file covers user routines that are used for the Markov Chain Monte Carlo algorithm tools.
 *
 */

#ifndef _LALINSPIRALMCMCUSER_H
#define _LALINSPIRALMCMCUSER_H

# include <math.h>
# include <stdio.h>
# include <stdlib.h>


# include <lal/LALStdlib.h>
# include <lal/LALConstants.h>
# include <lal/SimulateCoherentGW.h>
# include <lal/GeneratePPNInspiral.h>
# include <lal/LIGOMetadataTables.h>
# include <lal/LALDatatypes.h>
# include <lal/FindChirp.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_rng.h>


/*#ifdef  __cplusplus
extern "C" {
#endif*/

/**\name Error Codes */ /*@{*/
#define LALINSPIRALH_ENULL           1
#define LALINSPIRALH_EMEM            2
#define LALINSPIRALH_EDIV0           3
#define LALINSPIRALH_ESIZE           4
#define LALINSPIRALH_ECHOICE         5
#define LALINSPIRALH_EORDER          6
#define LALINSPIRALH_EAPPROXIMANT    7

#define LALINSPIRALH_MSGENULL         "Arguments contained an unexpected null pointer"
#define LALINSPIRALH_MSGEMEM          "Memory allocation error"
#define LALINSPIRALH_MSGEDIV0         "Division by zero"
#define LALINSPIRALH_MSGESIZE         "Invalid input range"
#define LALINSPIRALH_MSGECHOICE       "Invalid choice for an input parameter"
#define LALINSPIRALH_MSGEORDER        "unknown order specified"
#define LALINSPIRALH_MSGEAPPROXIMANT  "Invalid model"
/*@}*/

extern gsl_rng *RNG;
extern double timewindow;

double mc2mass1(double mc, double eta);

double mc2mass2(double mc, double eta);

double mc2mt(double mc, double eta);

double m2eta(double m1, double m2);

double m2mc(double m1, double m2);

double logJacobianMcEta(double mc, double eta);


  /* Example functions */
void MCMCInitTest(
  LALMCMCParameter *parameter,
  SnglInspiralTable *input);

REAL8 MCMCLikelihoodTest(
    LALMCMCInput      *inputMCMC,
    LALMCMCParameter  *parameter);

INT4 MCMCPriorTest(
   REAL4             *logPrior,
   LALMCMCInput      *inputMCMC,
   LALMCMCParameter  *parameter);

REAL8 NestPrior(LALMCMCInput *inputMCMC,LALMCMCParameter *parameter);

REAL8 NestPriorHighMass(LALMCMCInput *inputMCMC,LALMCMCParameter *parameter);

REAL8 GRBPrior(LALMCMCInput *inputMCMC,LALMCMCParameter *parameter);

REAL8 NestPriorSkyLoc(LALMCMCInput *inputMCMC, LALMCMCParameter *parameter);

void NestInitInjNINJA(LALMCMCParameter *parameter, void *iT);

void NestInitInjNINJAHighMass(LALMCMCParameter *parameter, void *iT);


int ParamInRange(LALMCMCParameter *parameter);

extern REAL4Vector *model;
extern REAL4Vector *Tmodel;
extern REAL8Vector **topdown_sum;
extern REAL8 *normalisations;

REAL8 MCMCLikelihoodMultiCoherentAmpCor(LALMCMCInput *inputMCMC, LALMCMCParameter *parameter);

REAL8 MCMCSTLikelihoodMultiCoherentF(LALMCMCInput *inputMCMC,LALMCMCParameter *parameter);

REAL8 MCMCLikelihoodMultiCoherentF(LALMCMCInput *inputMCMC, LALMCMCParameter *parameter);

REAL8 MCMCLikelihoodMultiCoherentF_PhenSpin(LALMCMCInput *inputMCMC,LALMCMCParameter *parameter);

REAL8  MCMCLikelihoodMultiCoherentF_test(LALMCMCInput *inputMCMC,LALMCMCParameter *parameter,const char *outfile);

void TaylorT_template(LALStatus *status,InspiralTemplate *template, LALMCMCParameter *parameter,LALMCMCInput *inputMCMC) ;

void TaylorF2_template(LALStatus *status,InspiralTemplate *template, LALMCMCParameter *parameter,LALMCMCInput *inputMCMC) ;

void PhenSpinTaylorRD_template(LALStatus *status,InspiralTemplate *template, LALMCMCParameter *parameter,LALMCMCInput *inputMCMC) ;

void SpinTaylor_template(LALStatus *status,InspiralTemplate *template, LALMCMCParameter *parameter,LALMCMCInput *inputMCMC) ;

void IMRPhenomFA_template(LALStatus *status,InspiralTemplate *template, LALMCMCParameter *parameter,LALMCMCInput *inputMCMC) ;

void IMRPhenomFB_template(LALStatus *status,InspiralTemplate *template, LALMCMCParameter *parameter,LALMCMCInput *inputMCMC) ;

void IMRPhenomB_template(LALStatus *status,InspiralTemplate *template, LALMCMCParameter *parameter,LALMCMCInput *inputMCMC) ;

void EOBNR_template(LALStatus *status,InspiralTemplate *template, LALMCMCParameter *parameter,LALMCMCInput *inputMCMC) ;
#endif /* _LALINSPIRALMCMCUSER_H */

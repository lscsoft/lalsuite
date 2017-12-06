/*
*  Copyright (C) 2012 Matthew Pitkin, Colin Gill, John Veitch
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

/* functions to create the likelihood for a pulsar search to be used with the
LALInference tools */

/**
 * \file
 * \ingroup lalapps_pulsar_HeterodyneSearch
 * \author Matthew Pitkin, John Veitch, Max Isi, Colin Gill
 *
 * \brief Parameter estimation code for known pulsar searches using the nested sampling algorithm.
 *
 * ### Description ###
 *
 * This code is used to perform parameter estimation and evidence calculation in targeted/semi-targeted searches for
 * gravitational waves from known pulsars. It may also be used to follow-up on signal candidates from semi-coherent all-sky
 * searches for unknown sources.
 *
 * It uses the Bayesian technique of 'Nested Sampling' to sample over a defined prior parameter space (unknown signal
 * parameters such as the gravitational wave amplitude). These samples can then be used to create posterior probability
 * density functions of the unknown parameters. The samples are also used to calculate the Bayesian evidence, also known as
 * the marginal likelihood, for a given signal model. This can be compared with other models, in particular the model that
 * the data is described by Gaussian noise alone.
 *
 * As input the code requires time domain data that has been heterodyned using the known (or close to) phase evolution of
 * the pulsar. The time domain input should consist of a three column text file containing the GPS time stamp of the data
 * point, the real part of the heterodyned data and the imaginary part of the heterodyned data, e.g.
 * \code
 * 900000000.000000  1.867532e-24  -7.675329e-25
 * 900000060.000000  2.783651e-24  3.654386e-25
 * ...
 * \endcode
 *
 * Most commonly such data will have a sample rate of 1/60 Hz, giving a bandwidth of the same amount, but the code can
 * accept any rate.
 *
 * The code also requires that you specify which parameters are to be searched over, and the prior ranges over these. Any
 * of the signal parameters can be searched over, including frequency, sky position and binary system parameters, although
 * the bandwidth of the data and search efficiency need to be taken into account.
 *
 * The 'Nested Sampling' algorithm (developed by \cite Skilling2006) used is that defined in \c LALinferenceNestedSampler.c
 * (see \cite VeitchVecchio2010). It is essentially an efficient way to perform the integral
 * \f[
 * Z = \int^{\mathbf{\theta}} p(d|\mathbf{\theta}) p(\mathbf{\theta}) \mathrm{d}\mathbf{\theta},
 * \f]
 * where \f$ \mathbf{\theta} \f$ is a vector of parameters, \f$ p(d|\mathbf{\theta}) \f$ is the likelihood of the data
 * given the parameters, and \f$ p(\mathbf{\theta}) \f$ is the prior on the parameters. It does this by changing the
 * multi-dimensional integral over N parameters into a one-dimensional integral
 * \f[
 * Z = \int^X L(X) \mathrm{d}X \approx \sum_i L(X_i) \Delta{}X_i,
 * \f]
 * where \f$ L(X) \f$ is the likelihood, and \f$ X \f$ is the prior mass. The algorithm will draw a number (\f$ N \f$) of
 * samples (live points) from the parameter priors, calculate the likelihood for each point and find the lowest likelihood
 * value. The lowest likelihood value will be added to the summation in the above equation, with \f$ \log{\Delta{}X_i}
 * \approx 1/N \f$ coming from the fact that the prior would be normalised to unity and therefore each point should occupy
 * an equal fraction and at each iteration the prior volume will decrease geometrically (for \f$\log{\Delta{}X_0} = 0\f$).
 * A new point is then drawn from the prior with the criterion that it has a higher likelihood than the previous lowest
 * point and substitutes that point. To draw the new point a Markov Chain Monte Carlo (MCMC) procedure is used. The procedure
 * is continued until a stopping criterion is reached, which in this case is that the remaining prior volume is less than the
 * \c tolerance value set (see below). The implementation of this can be seen in \cite VeitchVecchio2010 .
 *
 * ### Usage ###
 *
 * The usage format is given below and can also be found by running the code with
 * \code
 * lalapps_pulsar_parameter_estimation_nested --help
 * \endcode
 *
 * An example of running the code to search over the four unknown parameters \f$ h_0 \f$, \f$ \phi_0 \f$, \f$ \psi \f$
 * and \f$ \cos{\iota} \f$, for pulsar J0534-2200, given heterodyned time domain data from the H1 detector in the file
 * \c finehet_J0534-2200_H1, is:
 * \code
 * lalapps_pulsar_parameter_estimation_nested --detectors H1 --par-file J0534-2200.par --input-files finehet_J0534-2200_H1 --outfile ns_J0534-2200.hdf --prior-file prior_J0534-2200.txt --ephem-earth lscsoft/share/lalpulsar/earth05-09.dat --ephem-sun lscsoft/share/lalpulsar/sun05-09.dat --Nlive 1000 --Nmcmcinitial 0 --tolerance 0.25
 * \endcode
 * The \c par-file is a TEMPO(2)-style file containing the parameters of the pulsar used to perform the heterodyne (the
 * frequency parameters are the rotation frequency and therefore not necessarily the gravitational wave frequency) e.g.
 * \code
 * RA      12:54:31.87523895
 * DEC     -54:12:43.6572033
 * PMRA    1.7
 * PMDEC   2.8
 * POSEPOCH 54320.8531
 * F0  123.7896438753
 * F1  4.592e-15
 * PEPOCH 54324.8753
 * \endcode
 * The \c prior-file is a text file containing a list of the parameters to be searched over, the prior type ("uniform" or
 * "gaussian") and their given lower/mean and upper/standard deviation ranges e.g.
 * \code
 * h0 uniform 0 1e-21
 * phi0 uniform 0 6.283185307179586
 * cosiota uniform -1 1
 * psi uniform -0.785398163397448 0.785398163397448
 * \endcode
 * Note that if searching over frequency parameters the ranges specified in the \c prior-file should be given in terms of
 * the pulsars rotation frequency and not necessarily the gravitational wave frequency e.g. for a triaxial star emitting
 * gravitational waves at 100 Hz (which will be at twice the rotation frequency) if you wanted to search over 99.999 to
 * 100.001 Hz then you should used
 * \code
 * f0 uniform 49.9995 50.0005
 * \endcode
 *
 * An example of running the code as above, but this time on fake data created using the Advanced LIGO design noise curves
 * and with a signal injected into the data is:
 * \code
 * lalapps_pulsar_parameter_estimation_nested --fake-data AH1 --inject-file fake.par --par-file fake.par --outfile ns_fake.hdf --prior-file prior_fake.txt --ephem-earth lscsoft/share/lalpulsar/earth05-09.dat --ephem-sun lscsoft/share/lalpulsar/sun05-09.dat --Nlive 1000 --Nmcmcinitial 0 --tolerance 0.25
 * \endcode
 * In this case the \c inject-file parameter file must contain the values of \c h0, \c phi0, \c psi and \c cosiota,
 * otherwise these will be set to zero by default. The parameter files given for \c inject-file and \c par-file do not
 * have to be the same - the injection can be offset from the 'heterodyned' values, which will be reflected in the data. If
 * an \c inject-output file is also specified then the fake data containing the signal, and a fake signal-only data set,
 * will be output.
 */

#include "config.h"
#include "pulsar_parameter_estimation_nested.h"
#include "ppe_utils.h"
#include "ppe_init.h"
#include "ppe_readdata.h"
#include "ppe_inject.h"
#include "ppe_models.h"
#include "ppe_likelihood.h"
#include "ppe_testing.h"
#include "ppe_roq.h"

/* global variables */
LALStringVector *corlist = NULL;

INT4 main( INT4 argc, CHAR *argv[] ){
  ProcessParamsTable *param_table;
  LALInferenceRunState runState;
  REAL8 logZnoise = 0.;

  /* set error handler to abort in main function */
  XLALSetErrorHandler( XLALExitErrorHandler );

  /* Get ProcParamsTable from input arguments */
  param_table = LALInferenceParseCommandLine( argc, argv );
  runState.commandLine = param_table;

  /* Initialise the algorithm structures from the command line arguments */
  /* Include setting up random number generator etc */
  initialise_algorithm( &runState );

  /* read in data */
  read_pulsar_data( &runState );

  /* set algorithm to use Nested Sampling */
  runState.algorithm = &nested_sampling_algorithm_wrapper;
  runState.evolve = &LALInferenceNestedSamplingOneStep;

  /* set likelihood function */
  runState.likelihood = &pulsar_log_likelihood;

  /* set prior function */
  runState.prior = &priorFunction;

  /* Generate the lookup tables and read parameters from par file */
  setup_from_par_file( &runState );

  /* set signal model/template */
  runState.threads[0]->model->templt = &get_pulsar_model;

  /* add injections if requested */
  inject_signal( &runState );

  /* Initialise the prior distribution given the command line arguments */
  initialise_prior( &runState );

  /* create sum square of the data to speed up the likelihood calculation */
  sum_data( &runState );

  /* check whether using reduced order quadrature */
  generate_interpolant( &runState );

  gridOutput( &runState );

  /* get noise likelihood and add as variable to runState */
  logZnoise = noise_only_likelihood( &runState );
  LALInferenceAddVariable( runState.algorithmParams, "logZnoise", &logZnoise, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );

  /* Create live points array and fill initial parameters */
  if( !LALInferenceGetProcParamVal(param_table, "--compare-likelihoods") ){
    setup_live_points_array_wrapper( &runState );
  }

  /* output the live points sampled from the prior */
  outputPriorSamples( &runState );

  /* Initialise the MCMC proposal distribution */
  initialise_proposal( &runState );

  /* Set up threads */
  initialise_threads( &runState, 1 );

  if( !LALInferenceGetProcParamVal(param_table, "--compare-likelihoods") ){
    /* Call the nested sampling algorithm */
    runState.algorithm( &runState );
  }
  else{
    /* compare likelihoods from previous run */
    compare_likelihoods( &runState );
    return 0;
  }

  /* get SNR of highest likelihood point */
  get_loudest_snr( &runState );

  /* close timing file */
  if ( LALInferenceCheckVariable( runState.algorithmParams, "timefile" ) ){
    fclose(*(FILE**)LALInferenceGetVariable( runState.algorithmParams, "timefile" ));
  }

  return 0;
}

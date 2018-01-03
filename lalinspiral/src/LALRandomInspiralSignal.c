/*
*  Copyright (C) 2007 Stas Babak, David Churches, Duncan Brown, Jolien Creighton, B.S. Sathyaprakash, Anand Sengupta, Thomas Cokelaer
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

/**
 * \author Sathyaprakash, B. S., Thomas Cokelaer, Anand S. Sengupta
 * \file
 *
 * \brief Module to generate inspiral signals and simulated Gaussian noise
 *
 * Generate:<br>
 * (a) inspiral signals with random masses or chirp times
 * that have values within the parameter space specified by an input struct,
 *
 * (b) simulated Gaussian noise of PSD expected in a given interferometer
 *
 * (c) inspiral signal as in (a) but of a specified amplitude added
 * to simulated Gaussian noise as in (b).
 *
 * In all cases the returned vector is the Fourier transform of the relevant signal.
 *
 * ### Prototypes ###
 *
 * <tt>LALRandomInspiralSignal()</tt>
 *
 * ### Description ###
 *
 * The function receives input struct of type ::RandomInspiralSignalIn
 * whose members are
 * \code
 * typedef struct
 * tagRandomInspiralSignalIn
 * {
 * INT4 useed;
 * INT4 type;
 *
 * REAL8 mMin;
 * REAL8 mMax;
 * REAL8 MMax;
 * REAL8 SignalAmp;
 * REAL8 NoiseAmp;
 * REAL8 etaMin;
 * REAL8 t0Min;
 * REAL8 t0Max;
 * REAL8 tnMin;
 * REAL8 tnMax;
 *
 * InspiralTemplate param;
 * REAL8Vector psd;
 * RealFFTPlan *fwdp;
 * } RandomInspiralSignalIn;
 * \endcode
 *
 * Depending on the value of the parameter (<tt>randIn.type</tt>) this
 * code returns the Fourier transform of
 *
 * (a) a pure inspiral signal of a given type
 * (<tt>randIn.type=0</tt>),
 *
 * (b) simulated noise expected
 * in a chosen interferometer <tt>randIn.type=1</tt> or
 *
 * (c)
 * \f$\mathtt{SignalAmp}\times s+\mathtt{NoiseAmp}\times n\f$ (<tt>randIn.type=2</tt>),
 * where \f$s\f$ is normalised signal and \f$n\f$ random Gaussian noise whose PSD is
 * that expected in a given interferometer with zero mean and unit rms.
 *
 * User must specify the following quantities in the input structure
 *
 * <table class="doxtable" align="center">
 * <caption align="top" style="text-align: left; font-weight: normal;">Table: Input structure needed for the function LALRandomInspiralSignal().</caption>
 * <tr><th>Parameter</th><th>i/o</th><th>Comment</th></tr>
 * <tr><td><tt>INT4 useed</tt></td><td>input</td><td>Seed for the random number generator</td></tr>
 * <tr><td><tt>INT4 type</tt></td><td>input</td><td>Type of signal required to be generated</td></tr>
 * <tr><td><tt>InspiralTemplate p</tt></td><td>i/o</td><td>user must input certain params; others will be output</td></tr>
 * <tr><td><tt>p.startTime</tt></td><td></td><td>usually 0.</td></tr>
 * <tr><td><tt>p.startPhase</tt></td><td></td><td>\f$[0,\pi/2]\f$</td></tr>
 * <tr><td><tt>p.nStartPad</tt></td><td></td><td>number of zeros in the vector before the signal begins</td></tr>
 * <tr><td><tt>p.nEndPad</tt></td><td></td><td>number of zeros in the vector after the signal ends</td></tr>
 * <tr><td><tt>p.signalAmplitude</tt></td><td></td><td>usually 1</td></tr>
 * <tr><td><tt>p.ieta</tt></td><td></td><td>1 for comparable mass systems 0 for test mass model</td></tr>
 * <tr><td><tt>p.fLower</tt></td><td></td><td>lower frequency cutoff in Hz</td></tr>
 * <tr><td><tt>p.fCutoff</tt></td><td></td><td>upper frequency cutoff in Hz</td></tr>
 * <tr><td><tt>p.tSampling</tt></td><td></td><td>sampling rate in Hz</td></tr>
 * <tr><td><tt>p.order</tt></td><td></td><td>order of the PN approximant of the signal</td></tr>
 * <tr><td><tt>p.approximant</tt></td><td></td><td>PN approximation to be used for inspiral signal generation</td></tr>
 * <tr><td><tt>InputMasses massChoice</tt></td><td>input</td><td>space in which parameters are chosen; #m1Andm2, #totalMassAndEta, #totalMassUAndEta, #t02, #t03, #bhns</td></tr>
 * <tr><td><tt>REAL8Vector psd</tt></td><td>input</td><td>pre-computed power spectral density used for coloring the noise</td></tr>
 * <tr><td><tt>RealFFTPlan *fwdp</tt></td><td>input</td><td>pre-computed fftw plan to compute forward Fourier transform</td></tr>
 * <tr><td><tt>REAL8 mMin</tt></td><td>input</td><td>smallest component mass allowed</td></tr>
 * <tr><td><tt>REAL8 mMax</tt></td><td>input</td><td>largest component mass allowed   \c OR</td></tr>
 * <tr><td><tt>REAL8 MMax</tt></td><td>input</td><td>largest total mass allowed</td></tr>
 * <tr><td><tt>REAL8 SignalAmp</tt></td><td>input</td><td>amplitude of the signal (relevant only when <tt>type=2</tt>)</td></tr>
 * <tr><td><tt>REAL8 NoiseAmp</tt></td><td>input</td><td>amplitude of noise (relevant only when <tt>type=2</tt>)</td></tr>
 * <tr><td><tt>REAL8 etaMin</tt></td><td>input</td><td>smallest value of the symmetric mass ratio</td></tr>
 * <tr><td colspan="3"><center>Following chirp times are needed only if <tt>param.massChoice</tt> is #t02 \c or #t03 </center></td></tr>
 * <tr><td><tt>REAL8 t0Min</tt></td><td>input</td><td>smallest Newtonian chirp time</td></tr>
 * <tr><td><tt>REAL8 t0Max</tt></td><td>input</td><td>largest Newtonian chirp time</td></tr>
 * <tr><td><tt>REAL8 tnMin</tt></td><td>input</td><td>smallest 1 chirp time if <tt>param.massChoice=t02</tt></td></tr>
 * <tr><td></td><td></td><td>smallest 1.5 chirp time if <tt>param.massChoice=t03</tt></td></tr>
 * <tr><td><tt>REAL8 tnMax</tt></td><td>input</td><td>largest 1 chirp time  if <tt>param.massChoice=t02</tt></td></tr>
 * <tr><td></td><td></td><td>largest 1.5 chirp time  if <tt>param.massChoice=t03</tt></td></tr>
 * </table>
 *
 * When repeatedly called, the parameters of the signal will be
 * uniformly distributed in the space of
 *
 * (a) component masses in the range <tt>[randIn.mMin, randIn.mMax]</tt> if
 * <tt>param.massChoice=m1Andm2</tt>,
 *
 * (b) component masses greater than <tt>randIn.mMin</tt> and total mass
 * less than <tt>randIn.MMax</tt> if  <tt>param.massChoice=totalMassAndEta</tt>,
 *
 * (c) component masses greater than <tt>randIn.mMin</tt> and \c uniform total mass
 * less than <tt>randIn.MMax</tt> if  <tt>param.massChoice=totalMassUAndEta</tt>,
 *
 * (d) Newtonian and first post-Newtonian chirp times if
 * <tt>param.massChoice=t02</tt>,
 *
 * (e) Newtonian and 1.5 post-Newtonian chirp times if
 * <tt>param.massChoice=t03</tt> and.
 *
 * (f) component masses in the range <tt>[randIn.mMin, randIn.mMax]</tt> one of them being a neutron
 * start and the other a black hole (one above 3 solar mass and one below) if
 * <tt>param.massChoice=bhns</tt>. The function therefore checks the mass range validity
 * i.e. randIn.mMin must be less than 3 and randIn.mMax greater than 3.
 *
 * ### Algorithm ###
 *
 * No special algorithm, only a series of calls to pre-existing functions.
 *
 * ### Uses ###
 *
 * \code
 * random()
 * LALInspiralParameterCalc()
 * LALInspiralWave()
 * LALREAL4VectorFFT()
 * LALInspiralWaveNormaliseLSO()
 * LALCreateRandomParams()
 * LALNormalDeviates()
 * LALDestroyRandomParams()
 * LALREAL4VectorFFT()
 * LALColoredNoise()
 * LALAddVectors()
 * \endcode
 *
 * ### Notes ###
 *
 */
#include <lal/LALStdlib.h>
#include <lal/LALNoiseModelsInspiral.h>
#include <lal/Random.h>
#include <lal/GenerateInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/SkyCoordinates.h>

#define random() rand()
#define srandom( seed ) srand( seed )

static void GenerateRandomSpinTaylorParameters (
        LALStatus *status,
        RandomInspiralSignalIn *randIn
        );
static void GenerateRandomSkyPositionAndPolarisation (
        LALStatus *status,
        RandomInspiralSignalIn *randIn
        );

void LALRandomInspiralSignal
(
 LALStatus              *status,
 REAL4Vector            *signalvec,
 RandomInspiralSignalIn *randIn
 )
{

    REAL8                   maxTemp; /* temporary variable */
    INT4                    iMax;    /* temporary index    */
    UINT4                   indice;
    REAL8                   epsilon1, epsilon2, norm;
    REAL4Vector             noisy, buff;
    AddVectorsIn            addIn;
    INT4                    valid;
    static RandomParams     *randomparams;
    InspiralWaveNormaliseIn normin;

    INITSTATUS(status);
    ATTATCHSTATUSPTR(status);

    ASSERT (signalvec->data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
    ASSERT (randIn->psd.data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
    ASSERT (randIn->mMin > 0, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
    ASSERT (randIn->MMax > 2*randIn->mMin, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
    ASSERT (randIn->type >= 0, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
    ASSERT (randIn->type <= 2, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);

    buff.length = signalvec->length;
    if (!(buff.data = (REAL4*) LALCalloc(buff.length, sizeof(REAL4)) )) {
        ABORT (status, LALNOISEMODELSH_EMEM, LALNOISEMODELSH_MSGEMEM);
    }

    /* Use the seed to initialize random(). */
    srandom(randIn->useed);
    /* use the random number so generated as the next seed */
    randIn->useed = random();

    /* we need random parameters only if we need to generate a signal (i.e. type 0/2) */
    if (randIn->type==0 || randIn->type==2)
    {
        valid = 0;
        /* Keep generating random parameters until they
         * are located within the specified region */
        while (!valid)
        {
            epsilon1 = (float) random()/(float)RAND_MAX;
            epsilon2 = (float) random()/(float)RAND_MAX;
            switch (randIn->param.massChoice)
            {

                case bhns:
                    ASSERT(randIn->mMin<=3 && randIn->mMax>=3, status, 10,
                            "if massChoice is set to bhns, mass1 must be <= 3 "
                            "solar mass and mass2 >= 3 solar mass\n");
                    randIn->param.mass1 = randIn->mMin +
                            (3 - randIn->mMin) * epsilon1;
                    randIn->param.mass2 = 3  +
                            (randIn->mMax - 3) * epsilon2;
                    randIn->param.massChoice=m1Andm2;
                    LALInspiralParameterCalc(status->statusPtr, &(randIn->param));
                    CHECKSTATUSPTR(status);
                    randIn->param.massChoice=bhns;
                    break;

                case m1Andm2:
                    /*
                     * restriction is on the minimum and maximum individual
                     * masses of the two component stars.
                     */
                    randIn->param.mass1 = randIn->mMin +
                            (randIn->mMax - randIn->mMin) * epsilon1;
                    randIn->param.mass2 = randIn->mMin  +
                            (randIn->mMax - randIn->mMin) * epsilon2;
                    LALInspiralParameterCalc(status->statusPtr, &(randIn->param));
                    CHECKSTATUSPTR(status);
                    break;

                case minmaxTotalMass:
                    /*
                     * restriction is on the min and max Total mass. Should be
                     * check carefully. Right now, I think that it is almost
                     * uniformly distributed in total mass but seem to drop at
                     * very high mass. I'm not sure about the etaMin either. I
                     * might remove this code anyway soon. This is for a quick
                     * test for Craig Robinson and the high mass CBC search
                     * */
                    {
                     REAL8 etaMin ;
                     randIn->param.totalMass = randIn->MMin + (randIn->MMax - randIn->MMin)*epsilon1 ;

                    if (randIn->param.totalMass < (randIn->mMin+ randIn->mMax)) {
                      etaMin = (randIn->mMin / randIn->param.totalMass);
                      etaMin = etaMin - etaMin *etaMin;
                    }
                    else {
                      etaMin = (randIn->mMax / randIn->param.totalMass);
                      etaMin = etaMin - etaMin *etaMin;
                    }
                    randIn->param.eta = etaMin + epsilon2 * (.25 - etaMin);

                    }
                    randIn->param.massChoice = totalMassAndEta;
                    LALInspiralParameterCalc(status->statusPtr, &(randIn->param));
                    CHECKSTATUSPTR(status);
                    randIn->param.massChoice = minmaxTotalMass;
                    break;


                case totalMassAndEta:
                    /*
                     * restriction is on the total mass of the binary
                     * and the minimum mass of the component stars
                     */


                    randIn->param.mass1 = randIn->mMin
                            + (randIn->MMax - 2.*randIn->mMin) * epsilon1;
                    randIn->param.mass2 = randIn->mMin
                            + (randIn->MMax - randIn->param.mass1 - randIn->mMin) * epsilon2;
                    randIn->param.totalMass = randIn->param.mass1 + randIn->param.mass2 ;
                    randIn->param.eta = (randIn->param.mass1*randIn->param.mass2) / pow(randIn->param.totalMass,2.L);
                    LALInspiralParameterCalc(status->statusPtr, &(randIn->param));
                    CHECKSTATUSPTR(status);
                    break;

                case totalMassUAndEta:
                    /*
                     * restriction is on the total mass of the binary
                     * and the etaMin which depends on max total mass.
                     */
                    {
                        REAL4 etaMin;

                        randIn->param.totalMass =  2*randIn->mMin  +  epsilon1 * (randIn->MMax - 2 * randIn->mMin) ;

                        if (randIn->param.totalMass < (randIn->mMin+ randIn->mMax)) {
                            etaMin = (randIn->mMin / randIn->param.totalMass);
                            etaMin = etaMin - etaMin *etaMin;
                        }
                        else {
                            etaMin = (randIn->mMax / randIn->param.totalMass);
                            etaMin = etaMin - etaMin *etaMin;
                        }
                        randIn->param.eta = etaMin + epsilon2 * (.25 - etaMin);

                        LALInspiralParameterCalc(status->statusPtr, &(randIn->param));
                        CHECKSTATUSPTR(status);
                    }
                    break;

                case fixedMasses: /* the user has already given individual masses*/
                    randIn->param.massChoice = m1Andm2;
                    LALInspiralParameterCalc(status->statusPtr, &(randIn->param));
                    CHECKSTATUSPTR(status);
                    randIn->param.massChoice = fixedMasses;
                    break;

                case fixedPsi: /* the user has already given psi0/psi3*/
                    randIn->param.massChoice = psi0Andpsi3;
                    LALInspiralParameterCalc(status->statusPtr, &(randIn->param));
                    CHECKSTATUSPTR(status);
                    randIn->param.massChoice = fixedPsi;
                    break;

                case fixedTau: /* the user has already given tau0/tau3*/
                    randIn->param.massChoice = t03;
                    LALInspiralParameterCalc(status->statusPtr, &(randIn->param));
                    CHECKSTATUSPTR(status);
                    randIn->param.massChoice = fixedTau;
                    break;

                case t02:
                    /* chirptimes t0 and t2 are required in a specified range */
                    randIn->param.t0 = randIn->t0Min +
                            (randIn->t0Max - randIn->t0Min)*epsilon1;
                    randIn->param.t2 = randIn->tnMin +
                            (randIn->tnMax - randIn->tnMin)*epsilon2;
                    LALInspiralParameterCalc(status->statusPtr, &(randIn->param));
                    CHECKSTATUSPTR(status);
                    break;

                case t03:
                    /* chirptimes t0 and t3 are required in a specified range */
                    randIn->param.t0 = randIn->t0Min +
                            (randIn->t0Max - randIn->t0Min)*epsilon1;
                    randIn->param.t3 = randIn->tnMin +
                            (randIn->tnMax - randIn->tnMin)*epsilon2;
                    LALInspiralParameterCalc(status->statusPtr, &(randIn->param));
                    CHECKSTATUSPTR(status);
                    break;

                case psi0Andpsi3:
                    /* BCV parameters are required in a specified range */
                    randIn->param.psi0 = randIn->psi0Min +
                            (randIn->psi0Max - randIn->psi0Min)*epsilon1;
                    randIn->param.psi3 = randIn->psi3Min +
                            (randIn->psi3Max - randIn->psi3Min)*epsilon2;
                    break;

                case massesAndSpin:
                    /* masses, spin parameters and sky position needs to be set */

                    /* Set the random masses first */
                    randIn->param.mass1 = randIn->mMin +
                            (randIn->mMax - randIn->mMin) * epsilon1;
                    randIn->param.mass2 = randIn->mMin  +
                            (randIn->mMax - randIn->mMin) * epsilon2;

                    /* Set the random spin parameters */
                    GenerateRandomSpinTaylorParameters ( status->statusPtr, randIn );
                    CHECKSTATUSPTR(status);

                    /* Set the random sky position and polarisation angle */
                    GenerateRandomSkyPositionAndPolarisation ( status->statusPtr,randIn );
                    CHECKSTATUSPTR(status);

                    LALInspiralParameterCalc(status->statusPtr, &(randIn->param));
                    CHECKSTATUSPTR(status);
                    break;

                case t04:
                default:
                    /* if the choice of parameters is wrong abort the run */
                    ABORT (status, LALNOISEMODELSH_ECHOICE, LALNOISEMODELSH_MSGECHOICE);
                    break;
            }


            /* Validate the random parameters generated above */
            switch (randIn->param.massChoice)
            {
                case bhns:
                    valid=1;
                    break;
                case minmaxTotalMass:
                    if (
                            randIn->param.mass1 >= randIn->mMin &&
                            randIn->param.mass2 >= randIn->mMin &&
                            randIn->param.mass1 <= randIn->mMax &&
                            randIn->param.mass2 <= randIn->mMax &&
                            (randIn->param.eta > randIn->etaMin) &&
                            (randIn->param.mass1+randIn->param.mass2) < randIn->MMax  &&
                            (randIn->param.mass1+randIn->param.mass2) > randIn->MMin
                       )
                    {
                        valid = 1;
                    }
                     break;

                case fixedMasses:
                case fixedTau:
                  valid = 1;
                  break;
                case m1Andm2:
                case t03:
                case t02:
                    /*
                     * The following imposes a range in which min and
                     * max of component masses are restricted.
                     */
                    if (
                            randIn->param.mass1 >= randIn->mMin &&
                            randIn->param.mass2 >= randIn->mMin &&
                            randIn->param.mass1 <= randIn->mMax &&
                            randIn->param.mass2 <= randIn->mMax &&
                            randIn->param.eta <= 0.25 &&
                            randIn->param.eta >= randIn->etaMin
                       )
                    {
                        valid = 1;
                    }
                    break;

                case totalMassAndEta:
                case totalMassUAndEta:

                    /*
                     * The following imposes a range in which min of
                     * component masses and max total mass are restricted.
                     */
                    if (
                            randIn->param.mass1 >= randIn->mMin &&
                            randIn->param.mass2 >= randIn->mMin &&
                            randIn->param.totalMass <= randIn->MMax &&
                            randIn->param.eta <= 0.25 &&
                            randIn->param.eta >= randIn->etaMin &&
                            randIn->param.mass1 <= randIn->mMax &&
                            randIn->param.mass2 <= randIn->mMax
                       )

                    {
                        valid = 1;
                    }
                    break;
                case fixedPsi:
                    randIn->param.massChoice = psi0Andpsi3;
                    LALInspiralParameterCalc(status->statusPtr, &(randIn->param));
                    CHECKSTATUSPTR(status);
                    randIn->param.massChoice = fixedPsi;
                    valid = 1;
                case psi0Andpsi3:
                    /*
                     * the following makes sure that the BCV has
                     * a well defined end-frequency
                     */
                    randIn->param.massChoice = psi0Andpsi3;
                    LALInspiralParameterCalc(status->statusPtr, &(randIn->param));
                    CHECKSTATUSPTR(status);
                    valid = 1;
                    /*	       if (randIn->param.totalMass > 0.)
                               {
                               REAL8 fLR, fLSO, fend;
                               epsilon1 = (float) random()/(float)RAND_MAX;
                               fLR = 1.L/(LAL_PI * pow (3.L,1.5) * randIn->param.totalMass * LAL_MTSUN_SI);
                               fLSO = 1.L/(LAL_PI * pow (6.L,1.5) * randIn->param.totalMass * LAL_MTSUN_SI);
                               fend = fLSO + (fLR - fLSO) * epsilon1;
                               if (fend > randIn->param.tSampling/2. || fend < randIn->param.fLower) break;
                               randIn->param.fFinal = fend;
                               valid = 1;
                               }*/
                    break;

                case massesAndSpin:
                    if (
                            randIn->param.mass1 >= randIn->mMin &&
                            randIn->param.mass2 >= randIn->mMin &&
                            randIn->param.mass1 <= randIn->mMax &&
                            randIn->param.mass2 <= randIn->mMax &&
                            randIn->param.eta <= 0.25
                       )
                    {
                        valid = 1;
                    }
                    break;

                case t04:
                default:
                    ABORT (status, LALNOISEMODELSH_ECHOICE, LALNOISEMODELSH_MSGECHOICE);
                    break;
            }
        }
    }


    /* set up the structure for normalising the signal */
    normin.psd          = &(randIn->psd);
    normin.df           = randIn->param.tSampling / (REAL8) signalvec->length;
    normin.fCutoff      = randIn->param.fCutoff;
    normin.samplingRate = randIn->param.tSampling;

    switch (randIn->type)
    {
        case 0:

            /* First deal with the signal only case:
             * if the signal is generated in the Fourier domain no
             * need for Fourier transform
             */
            if (randIn->param.approximant == BCV ||
                    randIn->param.approximant == BCVSpin  ||
                    randIn->param.approximant == TaylorF1 ||
                    randIn->param.approximant == TaylorF2 ||
                    randIn->param.approximant == PadeF1)
            {
                LALInspiralWave(status->statusPtr, signalvec, &randIn->param);
                CHECKSTATUSPTR(status);
            }
            else /* Else - generate a time domain waveform */
            {
                /* Note that LALInspiralWave generates only the plus
                 * polarisation of the GW in the time/frequency domain.
                 * For SpinTaylor we need both plus and
                 * the cross polarisations - therefore for this case, we
                 * treat the waveform generation differently. For all other
                 * timedomain approximants, we recourse to calling
                 * LALInspiralWave () function.
                 */
                if (randIn->param.approximant == SpinTaylor)
                {
                    randIn->param.fFinal=0;
                    GenerateTimeDomainWaveformForInjection (status->statusPtr, &buff, &randIn->param);
                    CHECKSTATUSPTR(status);
                }
                else
                {
                    /* force to compute fFinal is it really necessary  ? */
                    randIn->param.fFinal=0;
                    LALInspiralWave(status->statusPtr, &buff, &randIn->param);
                    CHECKSTATUSPTR(status);
                }

                /* Once the time domain waveform has been generated, take its
                 * Fourier transform [i.e buff (t) ---> signal (f)].
                 */
                if (XLALREAL4VectorFFT(signalvec, &buff, randIn->fwdp) != 0)
                  ABORTXLAL(status);

            } /* End of else if time domain waveform */

            /* we might want to know where is the signalvec injected*/
            maxTemp  = 0;
            iMax     = 0;
            for ( indice = 0 ; indice< signalvec->length; indice++)
            {
                if (fabs(signalvec->data[indice]) > maxTemp){
                    iMax = indice;
                    maxTemp = fabs(signalvec->data[indice]);
                }
            }
            randIn->coalescenceTime = iMax;

            normin.fCutoff = randIn->param.fFinal;
            LALInspiralWaveNormaliseLSO(status->statusPtr, signalvec, &norm, &normin);
            CHECKSTATUSPTR(status);
            break;

        case 1:
            /*
             * next deal with the noise only case:
             */
            /*
                Old method of generating Gaussian noise
                LALGaussianNoise(status->statusPtr, &buff, &randIn->useed);
                */
            /*LALCreateRandomParams(status->statusPtr, &randomparams, randIn->useed);*/
            LALCreateRandomParams(status->statusPtr, &randomparams, randIn->useed);
            CHECKSTATUSPTR(status);
            LALNormalDeviates(status->statusPtr, &buff, randomparams);
            CHECKSTATUSPTR(status);
            LALDestroyRandomParams(status->statusPtr, &randomparams);
            CHECKSTATUSPTR(status);
            if (XLALREAL4VectorFFT(signalvec, &buff, randIn->fwdp) != 0)
              ABORTXLAL(status);
            LALColoredNoise(status->statusPtr, signalvec, randIn->psd);
            CHECKSTATUSPTR(status);

            /* multiply the noise vector by the correct normalisation factor */
            {
                double a2 = randIn->NoiseAmp * sqrt (randIn->param.tSampling)/2.L;
                UINT4 i;
                for (i=0; i<signalvec->length; i++) signalvec->data[i] *= a2;
            }
            break;

        default:
            /*
             * finally deal with the noise+signal only case:
             */
            noisy.length = signalvec->length;
            if (!(noisy.data = (REAL4*) LALMalloc(sizeof(REAL4)*noisy.length)))
            {
                if (buff.data != NULL) LALFree(buff.data);
                buff.data = NULL;
                ABORT (status, LALNOISEMODELSH_EMEM, LALNOISEMODELSH_MSGEMEM);
            }
            /*LALCreateRandomParams(status->statusPtr, &randomparams, randIn->useed);*/
            LALCreateRandomParams(status->statusPtr, &randomparams, randIn->useed);
            CHECKSTATUSPTR(status);
            LALNormalDeviates(status->statusPtr, &buff, randomparams);
            CHECKSTATUSPTR(status);
            LALDestroyRandomParams(status->statusPtr, &randomparams);
            CHECKSTATUSPTR(status);
            if (XLALREAL4VectorFFT(&noisy, &buff, randIn->fwdp) != 0)
              ABORTXLAL(status);
            LALColoredNoise(status->statusPtr, &noisy, randIn->psd);
            CHECKSTATUSPTR(status);

            if (randIn->param.approximant == BCV ||
                    randIn->param.approximant == BCVSpin  ||
                    randIn->param.approximant == TaylorF1 ||
                    randIn->param.approximant == TaylorF2 ||
                    randIn->param.approximant == PadeF1)
            {
                LALInspiralWave(status->statusPtr, &buff, &randIn->param);
                CHECKSTATUSPTR(status);
            }
            else
            {
                LALInspiralWave(status->statusPtr, signalvec, &randIn->param);
                CHECKSTATUSPTR(status);


                /* Now convert from time domain signal(t) ---> frequency
                 * domain waveform buff(f) i.e signal(t) -> buff(f)*/
                if (XLALREAL4VectorFFT(&buff, signalvec, randIn->fwdp) != 0)
                  ABORTXLAL(status);

            }

            /* we might want to know where is the signal injected*/
            maxTemp  = 0;
            iMax     = 0;
            for ( indice = 0 ; indice< signalvec->length; indice++)
            {
                if (fabs(signalvec->data[indice]) > maxTemp){
                    iMax = indice;
                    maxTemp = fabs(signalvec->data[indice]);
                }
            }
            randIn->coalescenceTime = iMax;

            normin.fCutoff = randIn->param.fFinal;
            LALInspiralWaveNormaliseLSO(status->statusPtr, &buff, &norm, &normin);
            CHECKSTATUSPTR(status);

            addIn.v1 = &buff;
            addIn.a1 = randIn->SignalAmp;
            addIn.v2 = &noisy;
            /* After experimenting we found that the following factor SamplingRate*sqrt(2)
             * is needed for noise amplitude to get the output of the signalless correlation
             * equal to unity; the proof of this is still needed
             */
            addIn.a2 = randIn->NoiseAmp * sqrt (randIn->param.tSampling)/2.L;
            LALAddVectors(status->statusPtr, signalvec, addIn);
            CHECKSTATUSPTR(status);
            if (noisy.data != NULL) LALFree(noisy.data);
            break;
    }

    /* Hash out signal for all frequencies less than or equal to the
     * lower cut-off frequency.
     */


    if (buff.data != NULL) LALFree(buff.data);

    DETATCHSTATUSPTR(status);
    RETURN(status);
}


/*----------- Functions static within this file -------------------*/

/* Function to generate random sky position and polarisation angle */
static void GenerateRandomSkyPositionAndPolarisation (
        LALStatus *status,
        RandomInspiralSignalIn *randIn )
{
    REAL8 u;
    REAL8 cosThetaMin, cosThetaMax;

    INITSTATUS(status);
    ATTATCHSTATUSPTR(status);

    ASSERT (randIn->sourceThetaMin <= randIn->sourceThetaMax ,
            status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
    ASSERT (randIn->sourceThetaMin >= 0.0 ,
            status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
    ASSERT (randIn->sourceThetaMax <= LAL_PI/2.0 ,
            status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
    ASSERT (randIn->sourcePhiMin <= randIn->sourcePhiMax ,
            status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
    ASSERT (randIn->sourcePhiMin >= 0.0 ,
            status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
    ASSERT (randIn->sourcePhiMax <= LAL_PI*2.0 ,
            status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
    ASSERT (randIn->polarisationAngleMin <= randIn->polarisationAngleMax ,
            status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
    ASSERT (randIn->polarisationAngleMin >= 0.0 ,
            status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
    ASSERT (randIn->polarisationAngleMax <= LAL_PI/2.0 ,
            status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);

    u = (float)(random())/(float)(RAND_MAX);
    cosThetaMax = cos(randIn->sourceThetaMax);
    cosThetaMin = cos(randIn->sourceThetaMin);
    randIn->param.sourceTheta = acos(cosThetaMin +
            u * (cosThetaMax - cosThetaMin) );
    u = (float)(random())/(float)(RAND_MAX);
    randIn->param.sourcePhi = randIn->sourcePhiMin +
            u * (randIn->sourcePhiMax - randIn->sourcePhiMin);
    u = (float)(random())/(float)(RAND_MAX);
    randIn->param.polarisationAngle = randIn->polarisationAngleMin +
            u * (randIn->polarisationAngleMax -
                    randIn->polarisationAngleMin);

    DETATCHSTATUSPTR(status);
    RETURN(status);
}


/* Function to generate random parameters for SpinTaylor injections */
static void GenerateRandomSpinTaylorParameters (
        LALStatus              *status,
        RandomInspiralSignalIn *randIn
        )
{
    REAL8 u;

    INITSTATUS(status);
    ATTATCHSTATUSPTR(status);

    /* Setting spin 1 */
    {
        REAL8  spin1Mag, r1, phi1;

        ASSERT (randIn->spin1min <= randIn->spin1max, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
        ASSERT (randIn->spin1min >= 0.0 , status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
        ASSERT (randIn->spin1max <= 1.0 , status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);

        u = (REAL8)(random())/(REAL8)(RAND_MAX);
        spin1Mag = randIn->spin1min +
                u * (randIn->spin1max - randIn->spin1min);
        u = (float)(random())/(float)(RAND_MAX);
        randIn->param.spin1[2] = (u - 0.5) * 2 * (spin1Mag);
        r1 = pow( ((spin1Mag * spin1Mag) - (randIn->param.spin1[2]
                        * randIn->param.spin1[2])) , 0.5);
        u = (float)(random())/(float)(RAND_MAX);
        phi1 = u * LAL_TWOPI;
        randIn->param.spin1[0]  = r1 * cos(phi1);
        randIn->param.spin1[1]  = r1 * sin(phi1);
    }

    /* Setting spin 2 */
    {
        REAL8 spin2Mag, r2, phi2;

        ASSERT (randIn->spin2min <= randIn->spin2max, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
        ASSERT (randIn->spin2min >= 0.0 , status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
        ASSERT (randIn->spin2max <= 1.0 , status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);

        u = (float)(random())/(float)(RAND_MAX);
        spin2Mag = randIn->spin2min +
                u * (randIn->spin2max - randIn->spin2min);
        u = (float)(random())/(float)(RAND_MAX);
        randIn->param.spin2[2] = (u - 0.5) * 2 * (spin2Mag);
        r2 = pow( ((spin2Mag * spin2Mag) - (randIn->param.spin2[2]
                        * randIn->param.spin2[2])) , 0.5);
        u = (float)(random())/(float)(RAND_MAX);
        phi2 = u * LAL_TWOPI;
        randIn->param.spin2[0]  = r2 * cos(phi2);
        randIn->param.spin2[1]  = r2 * sin(phi2);
    }

    /* Setting orbital parameters */
    {
        REAL8 cosTheta0Min, cosTheta0Max;

        ASSERT (randIn->theta0min <= randIn->theta0max , status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
        ASSERT (randIn->theta0min >= 0.0 , status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
        ASSERT (randIn->theta0max <= LAL_PI/2.0 , status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
        ASSERT (randIn->phi0min <= randIn->phi0max , status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
        ASSERT (randIn->phi0min >= 0.0 , status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
        ASSERT (randIn->phi0max <= LAL_PI*2.0 , status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);

        u = (float)(random())/(float)(RAND_MAX);
        cosTheta0Max = cos(randIn->theta0max);
        cosTheta0Min = cos(randIn->theta0min);
        randIn->param.orbitTheta0 = acos(cosTheta0Min +
                u * (cosTheta0Max - cosTheta0Min) );
        u = (float)(random())/(float)(RAND_MAX);
        randIn->param.orbitPhi0 = randIn->phi0min +
                u * (randIn->phi0max - randIn->phi0min);
    }

    /* Setting inclination angle */
    {
        ASSERT (randIn->inclinationMin <= randIn->inclinationMax , status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
        ASSERT (randIn->inclinationMin >= 0.0 , status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
        ASSERT (randIn->inclinationMax <= LAL_PI , status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);

        u = (float)(random())/(float)(RAND_MAX);
        randIn->param.inclination = randIn->inclinationMin +
                u * (randIn->inclinationMax - randIn->inclinationMin);
    }

    DETATCHSTATUSPTR(status);
    RETURN(status);
}

/* Function to generate time domain waveform for injection */
void GenerateTimeDomainWaveformForInjection (
        LALStatus              *status,
        REAL4Vector            *buff,
        InspiralTemplate       *param
        )
{

    CoherentGW          waveform;
    PPNParamStruc       ppnParams;
    INT4                nStartPad = 0;

    INITSTATUS(status);
    ATTATCHSTATUSPTR(status);


    memset( &waveform, 0, sizeof(CoherentGW) );

    ppnParams.deltaT               = 1.0 / param->tSampling;
    ppnParams.mTot                 = param->mass1 + param->mass2;
    ppnParams.eta                  = param->eta;
    ppnParams.d                    = param->distance;
    ppnParams.inc                  = param->inclination;
    ppnParams.phi                  = param->startPhase;
    ppnParams.fStartIn             = param->fLower;
    ppnParams.fStopIn              = -1.0 / (6.0 * sqrt(6.0) * LAL_PI * ppnParams.mTot * LAL_MTSUN_SI);
    ppnParams.position.longitude   = param->sourcePhi;
    ppnParams.position.latitude    = param->sourceTheta;
    ppnParams.position.system      = COORDINATESYSTEM_EQUATORIAL;
    ppnParams.psi                  = param->polarisationAngle;
    ppnParams.epoch.gpsSeconds     = 0;
    ppnParams.epoch.gpsNanoSeconds = 0;


    /* the waveform generation itself */

    /* Note that in the call to LALInspiralWaveForInjection,
     * the param.nStartPad will be reset to zero. We do
     * not want to lose the information. So we should save it
     * somewhere (in a temporary variable) before the function
     * call.
     */
    nStartPad = param->nStartPad;

    LALInspiralWaveForInjection(status->statusPtr, &waveform, param, &ppnParams);
    CHECKSTATUSPTR(status);

    /* Now reinstate nStartPad from saved value */
    param->nStartPad = nStartPad;

    /* Generate F+ and Fx and combine it with h+ and hx to get
     * the correct waveform. See equation from BCV2 (Eq 29,
     * 30).  */
    {
        REAL8 t, p, s, a1, a2, phi, phi0, shift;
        REAL8 fp, fc, hp, hc;
        UINT4 kk;

        t = cos(param->sourceTheta);
        p = 2. * param->sourcePhi;
        s = 2. * param->polarisationAngle;

        fp = 0.5*(1 + t*t)*cos(p)*cos(s) - t*sin(p)*sin(s);
        fc = 0.5*(1 + t*t)*cos(p)*sin(s) + t*sin(p)*cos(s);

        phi0 = waveform.phi->data->data[0];
        for (kk=0; kk < waveform.phi->data->length; kk++)
        {
            a1    = waveform.a->data->data[2*kk];
            a2    = waveform.a->data->data[2*kk+1];
/*            phi   = waveform.phi->data->data[kk] - phi0 -
 *            param->startPhase;*/
            phi   = waveform.phi->data->data[kk] - phi0;
            shift = waveform.shift->data->data[kk];
            hp    = a1*cos(shift)*cos(phi) - a2*sin(shift)*sin(phi);
            hc    = a1*sin(shift)*cos(phi) + a2*cos(shift)*sin(phi);
            buff->data[kk + param->nStartPad] = fp*hp + fc*hc;
        }

    LALSDestroyVectorSequence( status->statusPtr, &(waveform.a->data) );
    CHECKSTATUSPTR( status );
    LALSDestroyVector( status->statusPtr, &(waveform.f->data) );
    CHECKSTATUSPTR( status );
    LALDDestroyVector( status->statusPtr, &(waveform.phi->data) );
    CHECKSTATUSPTR( status );
    LALSDestroyVector( status->statusPtr, &(waveform.shift->data) );
    CHECKSTATUSPTR( status );
    LALFree( waveform.a );
    LALFree( waveform.f );
    LALFree( waveform.phi );
    LALFree( waveform.shift );
    }

    param->fFinal = ppnParams.fStop;
    DETATCHSTATUSPTR(status);
    RETURN(status);
}


INT4 XLALSimIMRSpinEOBWaveformTidal (COMPLEX16 * restrict hlm,
                                            /**< OUTPUT, hlm waveforms */
                                            REAL8Vector * restrict values,
                                            /**< dyanmical variables */
                                            const REAL8 v,
                                            /**< velocity */
                                            const INT4 l,
                                            /**< l mode index */
                                            const INT4 m,
                                            /**< m mode index */
                                            SpinEOBParams * restrict params
                                            /**< Spin EOB parameters */
);

INT4 XLALSimIMRSpinEOBFluxGetSpinFactorizedWaveform (COMPLEX16 * restrict hlm, REAL8Vector * restrict values, const REAL8 v, const REAL8 Hreal, const INT4 l, const INT4 m, SpinEOBParams * restrict params, INT4 use_optimized_v2,
                                       /**< Use optimized v2? */
                                REAL8 * vPhil2m2);

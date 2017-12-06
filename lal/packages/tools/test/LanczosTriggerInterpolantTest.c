/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA  02111-1307  USA
 *
 * Copyright (C) 2012 Leo Singer
 */


#include <complex.h>
#include <math.h>
#include <stdlib.h>

#include <lal/TriggerInterpolation.h>

int main(__attribute__((unused)) int argc, __attribute__((unused)) char **argv)
{
    int result;
    double tmax;
    COMPLEX16 ymax;

    LanczosTriggerInterpolant *interp = XLALCreateLanczosTriggerInterpolant(16);
    if (!interp)
        exit(EXIT_FAILURE);

    {
        /* Test some random sample data from original Python implementation. */
        const double expected_tmax = -0.35923390727111837;
        const double expected_re_ymax = 22.592995414852155;
        const double expected_im_ymax = -24.25569756344451;
        const COMPLEX16 y[] = {
            -14.04298458-12.42025719*I, -13.54512279-13.9975707*I ,
            -12.81987146-15.65632332*I, -11.83145220-17.25894236*I,
            -10.67978351-18.78930327*I,  -9.38321345-20.26279014*I,
            -7.93078349-21.75487168*I,  -6.20864157-23.23521056*I,
            -4.18656412-24.61039025*I,  -1.91352079-25.76366939*I,
            0.50485300-26.81349584*I,   3.24674313-27.88200749*I,
            6.61982798-28.75234467*I,  10.67117579-29.02191297*I,
            15.20752937-28.32539646*I,  19.86163875-26.31802974*I,
            23.95519327-22.82715585*I,  26.67973660-18.22160184*I,
            27.69204702-13.56075191*I,  27.61097719 -9.60342867*I,
            27.15105065 -6.36397975*I,  26.54148343 -3.43587179*I,
            25.59970437 -0.77676502*I,  24.42081565 +1.48165824*I,
            23.21572597 +3.26870707*I,  22.16738214 +4.89881504*I,
            21.06343508 +6.52136893*I,  19.79261510 +7.95655696*I,
            18.53300183 +9.1869077*I ,  17.27943443+10.35891335*I,
            15.85961229+11.44834434*I,  14.29135572+12.25702765*I,
            12.80306731+12.77700535*I};

        result = XLALCOMPLEX16ApplyLanczosTriggerInterpolant(interp, &tmax, &ymax, &y[16]);
        if (result)
            exit(EXIT_FAILURE);

        if (fabs(expected_tmax - tmax) > 1e-5)
            exit(EXIT_FAILURE);
        if (fabs((expected_re_ymax - creal(ymax)) / expected_re_ymax) > 1e-5)
            exit(EXIT_FAILURE);
        if (fabs((expected_im_ymax - cimag(ymax)) / expected_im_ymax) > 1e-5)
            exit(EXIT_FAILURE);
    }

    {
        /* Test some more random sample data from original Python implementation. */
        const double expected_tmax = 0.0080079263861989602;
        const double expected_re_ymax = 8.346533081471625;
        const double expected_im_ymax = 46.70076293771792;
        const COMPLEX16 y[] = {
            24.33546507 -4.37237844*I,  25.21430734 -2.94850093*I,
            26.24653016 -1.46698428*I,  27.43392439 +0.31641935*I,
            28.55556797 +2.4344684*I ,  29.53970968 +4.9329713*I ,
            30.24829181 +7.68385687*I,  30.74244848+10.63488387*I,
            31.05555847+13.79768684*I,  31.16910840+17.35620378*I,
            30.81906320+21.42675033*I,  29.71256983+25.80061299*I,
            27.85940204+30.3305615*I ,  25.15948657+34.9632285*I ,
            21.29005445+39.71231569*I,  15.69023426+44.03212162*I,
            8.40941881+46.68943372*I,   0.65145839+46.75721585*I,
            -6.03343803+44.61857557*I, -11.05937409+41.54119126*I,
            -14.90679982+38.44188372*I, -18.17229175+35.31228839*I,
            -20.85054167+32.08798444*I, -22.97315194+28.83554109*I,
            -24.64066997+25.66341162*I, -25.94581892+22.51742823*I,
            -26.82315799+19.46526862*I, -27.37637665+16.55731713*I,
            -27.69785243+13.83218299*I, -27.85474897+11.17681254*I,
            -27.73414707 +8.52204497*I, -27.23875261 +6.08054244*I,
            -26.58717550 +3.98060113*I};

        result = XLALCOMPLEX16ApplyLanczosTriggerInterpolant(interp, &tmax, &ymax, &y[16]);
        if (result)
            exit(EXIT_FAILURE);

        if (fabs(expected_tmax - tmax) > 1e-5)
            exit(EXIT_FAILURE);
        if (fabs((expected_re_ymax - creal(ymax)) / expected_re_ymax) > 1e-5)
            exit(EXIT_FAILURE);
        if (fabs((expected_im_ymax - cimag(ymax)) / expected_im_ymax) > 1e-5)
            exit(EXIT_FAILURE);
    }

    {
        /* Test yet more random sample data from original Python implementation. */
        const double expected_tmax = 0.35605452045274433;
        const double expected_re_ymax = 14.731795690833726;
        const double expected_im_ymax = 43.90822049035633;
        const COMPLEX16 y[] = {
            25.38672374-10.46990646*I,  26.98036717 -8.74105323*I,
            28.42577358 -6.75887934*I,  29.67590345 -4.50925463*I,
            30.62545698 -2.16923242*I,  31.44275377 +0.18575934*I,
            32.27276923 +2.57731966*I,  33.19687126 +5.33734012*I,
            33.89400960 +8.58856349*I,  34.13857476+12.22779037*I,
            33.87819938+16.02691928*I,  33.23756141+20.00079092*I,
            32.15025759+24.27606864*I,  30.38635264+28.89883078*I,
            27.61488144+33.82233139*I,  23.36306372+38.76765413*I,
            17.29380596+42.88182178*I,   9.89246651+44.9882819*I ,
            2.54059732+44.71650049*I,  -3.67718110+42.65873447*I,
            -8.48361774+40.01375385*I, -12.50375338+37.27989278*I,
            -16.04996471+34.31356422*I, -18.98265862+31.01719682*I,
            -21.19277639+27.6541388*I , -22.88668959+24.44207892*I,
            -24.24895445+21.30383581*I, -25.21392241+18.14916359*I,
            -25.71011223+15.13051498*I, -25.90832412+12.43072643*I,
            -26.04673459 +9.92719558*I, -26.08115387 +7.46027909*I,
            -25.90777937 +5.03328563*I};

        result = XLALCOMPLEX16ApplyLanczosTriggerInterpolant(interp, &tmax, &ymax, &y[16]);
        if (result)
            exit(EXIT_FAILURE);

        if (fabs(expected_tmax - tmax) > 1e-5)
            exit(EXIT_FAILURE);
        if (fabs((expected_re_ymax - creal(ymax)) / expected_re_ymax) > 1e-5)
            exit(EXIT_FAILURE);
        if (fabs((expected_im_ymax - cimag(ymax)) / expected_im_ymax) > 1e-5)
            exit(EXIT_FAILURE);
    }

    XLALDestroyLanczosTriggerInterpolant(interp);
    exit(EXIT_SUCCESS);
}

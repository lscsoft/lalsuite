/*                                           >y#
                                            ~'#o+
                                           '~~~md~
                '|+>#!~'':::::....        .~'~'cY#
            .+oy+>|##!~~~''':::......     ~:'':md! .
          #rcmory+>|#~''':::'::...::.::. :..'''Yr:...
        'coRRaamuyb>|!~'''::::':...........  .+n|.::..
       !maMMNMRYmuybb|!~'''':.........::::::: ro'..::..
      .cODDMYouuurub!':::...........:::~'.. |o>::...:..
      >BDNCYYmroyb>|#~:::::::::::::~':.:: :ob::::::::..
      uOCCNAa#'''''||':::.                :oy':::::::::.
    :rRDn!  :~::'y+::':  ... ...:::.     :ob':::::::::::.
   yMYy:   :>yooCY'.':.   .:'':......    ~u+~::::::::::::.
  >>:'. .~>yBDMo!.'': . .:'':.   .      >u|!:::::::::::::.
    ':'~|mYu#:'~'''. :.~:':...         yy>|~:::::::::::::..
    :!ydu>|!rDu::'. +'#~::!#'.~:     |r++>#':::::::::::::..
    mn>>>>>YNo:'': !# >'::::...  ..:cyb++>!:::::::::..:::...
    :ouooyodu:'': .!:.!:::.       yobbbb+>~::::::::....:....
     'cacumo~''' .'~ :~'.::.    :aybbbbbb>':::'~''::::....
      .mamd>'''. :~' :':'.:.   om>bbbyyyb>'.#b>|#~~~'':..
      .yYYo''': .:~' .'::'   .ny>+++byyoao!b+|||#!~~~''''''::.
      .#RUb:''. .:'' .:':   |a#|>>>>yBMdb #yb++b|':::::''':'::::::.
      .'CO!'''  .:'' .'    uu~##|+mMYy>+:|yyo+:::'::.         .::::::
      .:RB~''' ..::'.':   o>~!#uOOu>bby'|yB>.'::  '~!!!!!~':. ..  .::::
       :Rm''': ..:~:!:  'c~~+YNnbyyybb~'mr.':  !+yoy+>||!~'::.       :::.
      ..Oo''': .'' ~:  !+|BDCryuuuuub|#B!::  !rnYaocob|#!~'':.  ..    .::.
      . nB''': :  .'  |dNNduroomnddnuun::.  ydNAMMOary+>#~:.:::...      .:
       .uC~'''    :. yNRmmmadYUROMMBmm.:   bnNDDDMRBoy>|#~':....:.      .:
                 :' ymrmnYUROMAAAAMYn::. .!oYNDDMYmub|!~'::....:..     :
                 !'#booBRMMANDDDNNMO!:. !~#ooRNNAMMOOmuy+#!':::.......    :.
                .!'!#>ynCMNDDDDDNMRu.. '|:!raRMNAMOOdooy+|!~:::........   .:
                 : .'rdbcRMNNNNAMRB!:  |!:~bycmdYYBaoryy+|!~':::.::::.  ..
                 ..~|RMADnnONAMMRdy:. .>#::yyoroccruuybb>#!~'::':...::.
                  :'oMOMOYNMnybyuo!.  :>#::b+youuoyyy+>>|!~':.    :::::
                  ''YMCOYYNMOCCCRdoy##~~~: !b>bb+>>>||#~:..:::     ::::.
                  .:OMRCoRNAMOCROYYUdoy|>~:.~!!~!~~':...:'::::.   :::::.
                  ''oNOYyMNAMMMRYnory+|!!!:.....     ::.  :'::::::::::::
                 .:..uNabOAMMCOdcyb+|!~':::.          !!'.. :~:::::'''':.
                  .   +Y>nOORYauyy>!!'':....           !#~..  .~:''''''':.

****************  ____  _____  ______________________    ____     **************
***************  / __ )/   \ \/ / ____/ ___/_  __/   |  / __ \   ***************
**************  / __  / /| |\  / __/  \__ \ / / / /| | / /_/ /  ****************
*************  / /_/ / ___ |/ / /___ ___/ // / / ___ |/ _, _/  *****************
************  /_____/_/  |_/_/_____//____//_/ /_/  |_/_/ |_|  ******************
*/


/*
 * Copyright (C) 2013  Leo Singer
 *
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
 */


#ifndef BAYESTAR_SKY_MAP_H
#define BAYESTAR_SKY_MAP_H

#ifndef SWIG /* exclude from SWIG interface */

#include <complex.h>


/* Perform sky localization based on TDOAs, PHOAs, and amplitude. */
double *bayestar_sky_map_toa_phoa_snr(
    long *npix,
    /* Prior */
    double min_distance,            /* Minimum distance */
    double max_distance,            /* Maximum distance */
    int prior_distance_power,       /* Power of distance in prior */
    /* Detector network */
    double gmst,                    /* GMST (rad) */
    unsigned int nifos,             /* Number of detectors */
    unsigned long nsamples,         /* Length of autocorrelation sequence */
    double sample_rate,             /* Sample rate in seconds */
    const double complex **acors,   /* Autocorrelation sequences */
    const float (**responses)[3],   /* Detector responses */
    const double **locations,       /* Barycentered Cartesian geographic detector positions (m) */
    const double *horizons,         /* SNR=1 horizon distances for each detector */
    /* Observations */
    const double *toas,             /* Arrival time differences relative to network barycenter (s) */
    const double *phoas,            /* Phases on arrival */
    const double *snrs              /* SNRs */
);

double bayestar_log_likelihood_toa_snr(
    /* Parameters */
    double ra,                      /* Right ascension (rad) */
    double sin_dec,                 /* Sin(declination) */
    double distance,                /* Distance */
    double u,                       /* Cos(inclination) */
    double twopsi,                  /* Twice polarization angle (rad) */
    double t,                       /* Barycentered arrival time (s) */
    /* Detector network */
    double gmst,                    /* GMST (rad) */
    unsigned int nifos,             /* Number of detectors */
    unsigned long nsamples,         /* Length of autocorrelation sequence */
    double sample_rate,             /* Sample rate in seconds */
    const double complex **acors,   /* Autocorrelation sequences */
    const float (**responses)[3],   /* Detector responses */
    const double **locations,       /* Barycentered Cartesian geographic detector positions (m) */
    const double *horizons,         /* SNR=1 horizon distances for each detector */
    /* Observations */
    const double *toas,             /* Arrival time differences relative to network barycenter (s) */
    const double *snrs              /* SNRs */
);

double bayestar_log_likelihood_toa_phoa_snr(
    /* Parameters */
    double ra,                      /* Right ascension (rad) */
    double sin_dec,                 /* Sin(declination) */
    double distance,                /* Distance */
    double u,                       /* Cos(inclination) */
    double twopsi,                  /* Twice polarization angle (rad) */
    double t,                       /* Barycentered arrival time (s) */
    /* Detector network */
    double gmst,                    /* GMST (rad) */
    unsigned int nifos,             /* Number of detectors */
    unsigned long nsamples,         /* Length of autocorrelation sequence */
    double sample_rate,             /* Sample rate in seconds */
    const double complex **acors,   /* Autocorrelation sequences */
    const float (**responses)[3],   /* Detector responses */
    const double **locations,       /* Barycentered Cartesian geographic detector positions (m) */
    const double *horizons,         /* SNR=1 horizon distances for each detector */
    /* Observations */
    const double *toas,             /* Arrival time differences relative to network barycenter (s) */
    const double *phoas,            /* Phases on arrival */
    const double *snrs              /* SNRs */
);

/* Unit test suite. Return EXIT_SUCCESS if tests passed,
 * or otherwise EXIT_FAILURE. */
int bayestar_test(void);

#endif /* SWIG */

#endif /* BAYESTAR_SKY_MAP_H */

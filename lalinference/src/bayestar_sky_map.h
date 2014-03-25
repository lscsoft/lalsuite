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


/* Perform sky localization based on TDOAs alone. */
double *bayestar_sky_map_toa(
    long *npix, /* In/out: number of HEALPix pixels. */
    double gmst, /* Greenwich mean sidereal time in radians. */
    int nifos, /* Input: number of detectors. */
    const double **locs, /* Input: array of detector positions. */
    const double *toas, /* Input: array of times of arrival. */
    const double *w_toas /* Input: sum-of-squares weights, (1/TOA variance)^2. */
);

double bayestar_log_posterior_toa(
    double ra,
    double sin_dec,
    double gmst, /* Greenwich mean sidereal time in radians. */
    int nifos, /* Input: number of detectors. */
    const double **locations, /* Pointers to locations of detectors in Cartesian geographic coordinates. */
    const double *toas, /* Input: array of times of arrival with arbitrary relative offset. (Make toas[0] == 0.) */
    const double *w_toas /* Input: sum-of-squares weights, (1/TOA variance)^2. */
);

double bayestar_log_posterior_toa_snr(
    double ra,
    double sin_dec,
    double distance,
    double u,
    double twopsi,
    double gmst, /* Greenwich mean sidereal time in radians. */
    int nifos, /* Input: number of detectors. */
    const float (**responses)[3], /* Pointers to detector responses. */
    const double **locations, /* Pointers to locations of detectors in Cartesian geographic coordinates. */
    const double *toas, /* Input: array of times of arrival with arbitrary relative offset. (Make toas[0] == 0.) */
    const double *snrs, /* Input: array of SNRs. */
    const double *w_toas, /* Input: sum-of-squares weights, (1/TOA variance)^2. */
    const double *horizons, /* Distances at which a source would produce an SNR of 1 in each detector. */
    int prior_distance_power /* Use a prior of (distance)^(prior_distance_power) */
);

/* Perform sky localization based on TDOAs and amplitude. */
double *bayestar_sky_map_toa_snr(
    long *npix, /* In/out: number of HEALPix pixels. */
    double gmst, /* Greenwich mean sidereal time in radians. */
    int nifos, /* Input: number of detectors. */
    const float (**responses)[3], /* Pointers to detector responses. */
    const double **locations, /* Pointers to locations of detectors in Cartesian geographic coordinates. */
    const double *toas, /* Input: array of times of arrival with arbitrary relative offset. (Make toas[0] == 0.) */
    const double *snrs, /* Input: array of SNRs. */
    const double *w_toas, /* Input: sum-of-squares weights, (1/TOA variance)^2. */
    const double *horizons, /* Distances at which a source would produce an SNR of 1 in each detector. */
    double min_distance,
    double max_distance,
    int prior_distance_power /* Use a prior of (distance)^(prior_distance_power) */
);

/* Perform sky localization based on TDOAs, PHOAs, and amplitude. */
double *bayestar_sky_map_toa_phoa_snr(
    long *npix, /* In/out: number of HEALPix pixels. */
    double gmst, /* Greenwich mean sidereal time in radians. */
    int nifos, /* Input: number of detectors. */
    const float (**responses)[3], /* Pointers to detector responses. */
    const double **locations, /* Pointers to locations of detectors in Cartesian geographic coordinates. */
    const double *toas, /* Input: array of times of arrival with arbitrary relative offset. (Make toas[0] == 0.) */
    const double *phoas, /* Input: array of phases of arrival with arbitrary relative offset. (Make phoas[0] == 0.) */
    const double *snrs, /* Input: array of SNRs. */
    const double *w_toas, /* Input: sum-of-squares weights, (1/TOA variance)^2. */
    const double *w1s, /* Input: first moments of angular frequency. */
    const double *w2s, /* Input: second moments of angular frequency. */
    const double *horizons, /* Distances at which a source would produce an SNR of 1 in each detector. */
    double min_distance,
    double max_distance,
    int prior_distance_power /* Use a prior of (distance)^(prior_distance_power) */
);

#endif /* BAYESTAR_SKY_MAP_H */

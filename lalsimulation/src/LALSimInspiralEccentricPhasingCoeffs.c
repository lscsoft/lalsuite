/*
*  Copyright (C) 2017 Maria Haney
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

#include <lal/LALConstants.h>
#include <lal/LALAtomicDatatypes.h>

#include <math.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/*
 * The leading-order eccentric corrections to the TaylorF2 phasing function.
 * Provides eccentric contributions up to 3PN phase order in modified-harmonic gauge.
 * (This is a 3PN extension of Eqs. A3 and A4 in arXiv:1602.03081)
 */
static void UNUSED
XLALSimInspiralPNPhasing_eccF2(
        PNPhasingSeries *pfv19by3, /**< \todo UNDOCUMENTED */
        PNPhasingSeries *pfv25by3, /**< \todo UNDOCUMENTED */
        PNPhasingSeries *pfv28by3, /**< \todo UNDOCUMENTED */
        PNPhasingSeries *pfv31by3, /**< \todo UNDOCUMENTED */
        PNPhasingSeries *pfv34by3, /**< \todo UNDOCUMENTED */
        PNPhasingSeries *pfv37by3, /**< \todo UNDOCUMENTED */
        const REAL8 m1, /**< Mass of body 1, in Msol */
        const REAL8 m2, /**< Mass of body 2, in Msol */
        const REAL8 eccentricity, /**< eccentricity at reference epoch (here: initial eccentricity at f_min) */
        const REAL8 f_min /**< starting frequency (Hz) */
        )
{
    const REAL8 e0 = eccentricity;
    const REAL8 f0 = f_min;

    const REAL8 mtot = m1 + m2;
    const REAL8 eta = m1*m2/mtot/mtot;
    const REAL8 m_sec = mtot * LAL_MTSUN_SI;  /* total mass in seconds */
    const REAL8 piM = LAL_PI * m_sec;
    const REAL8 v0 = cbrt(piM*f0);

    const REAL8 eta2 = eta * eta;
    const REAL8 eta3 = eta2 * eta;
    const REAL8 v02 = v0 * v0;
    const REAL8 v03 = v02 * v0;
    const REAL8 v04 = v02 * v02;
    const REAL8 v05 = v04 * v0;
    const REAL8 v06 = v03 * v03;

    const REAL8 pfeN = 3.L/(128.L * eta) * e0 * e0 * pow(v0, 19.L/3.L);

    memset(pfv19by3, 0, sizeof(PNPhasingSeries));
    memset(pfv25by3, 0, sizeof(PNPhasingSeries));
    memset(pfv28by3, 0, sizeof(PNPhasingSeries));
    memset(pfv31by3, 0, sizeof(PNPhasingSeries));
    memset(pfv34by3, 0, sizeof(PNPhasingSeries));
    memset(pfv37by3, 0, sizeof(PNPhasingSeries));

    /* Leading-order eccentric corrections to non-spin phasing terms */

    /* v^(-19/3) coefficients at every phase order */

    pfv19by3->v[0] = -2355.L/1462.L;
    pfv19by3->v[2] = -128365.L/12432.L * eta - 2045665.L/348096.L; 
    pfv19by3->v[3] = 65561.L/4080.L * LAL_PI; 
    pfv19by3->v[4] = -10688155.L/294624.L * eta2 - 111064865.L/14141952.L - 165068815.L/4124736.L * eta; 
    pfv19by3->v[5] = (3873451.L/100548.L + 15803101.L/229824.L * eta) * LAL_PI; 
    pfv19by3->vlogv[6] = -734341.L/16800.L;
    pfv19by3->v[6] = (-21508213.L/276480.L + 103115.L/6144.L * eta) * LAL_PI * LAL_PI - 734341.L/16800.L * LAL_GAMMA 
			- 409265200567.L/585252864.L * eta - 4726688461.L/34836480.L * eta2 - 69237581.L/746496.L * eta3 
			+ 60634674069696877.L/112661176320000.L - 9663919.L/50400.L * log(2.) + 4602177.L/44800.L * log(3.);

    /* v^(-25/3) coefficients at every phase order */

    pfv25by3->v[2] = (154645.L/17544.L * eta - 2223905.L/491232.L) * v02;
    pfv25by3->v[4] = (4917245.L/1566432.L * eta - 5795368945.L/350880768.L + 25287905.L/447552.L * eta2) * v02;
    pfv25by3->v[5] = (185734313.L/4112640.L - 12915517.L/146880.L * eta) * LAL_PI * v02;
    pfv25by3->v[6] = (-314646762545.L/14255087616.L + 11585856665.L/98993664.L * eta2 
			- 1733730575525.L/24946403328.L * eta + 2105566535.L/10606464.L * eta3) * v02; 

    /* v^(-28/3) coefficients at every phase order */

    pfv28by3->v[3] = -295945.L/35088.L * LAL_PI * v03;
    pfv28by3->v[5] = (-771215705.L/25062912.L - 48393605.L/895104.L * eta) * LAL_PI * v03; 
    pfv28by3->v[6] = 24716497.L/293760.L * v03 * LAL_PI * LAL_PI; 

    /* v^(-31/3) coefficients at every phase order */

    pfv31by3->v[4] = (936702035.L/1485485568.L - 14251675.L/631584.L * eta2 + 3062285.L/260064.L * eta) * v04; 
    pfv31by3->v[6] = (2440991806915.L/1061063442432.L - 2330466575.L/16111872.L * eta3 
			- 1029307085.L/150377472.L * eta2 + 1781120054275.L/37895122944.L * eta) * v04; 

    /* v^(-34/3) coefficients at every phase order */

    pfv34by3->v[5] = (149064749.L/2210544.L * eta - 7063901.L/520128.L) * LAL_PI * v05;

    /* v^(-37/3) coefficients at every phase order */

    pfv37by3->v[6] = ((-96423905.L/5052672.L + 3121945.L/561408.L * eta) * LAL_PI * LAL_PI + 2425890995.L/68211072.L * eta3 
			+ 1898287.L/184212.L * log(2.) + 2603845.L/61404.L * log(v0) - 4165508390854487.L/16471063977984.L 
			+ 12246471.L/163744.L * log(3.) - 1437364085977.L/53477480448.L * eta 
			+ 2603845.L/61404.L * LAL_GAMMA + 4499991305.L/636636672.L * eta2) * v06;

    /* At the very end, multiply everything in the series by pfeN */

     for(int ii = 0; ii <= PN_PHASING_SERIES_MAX_ORDER; ii++)
     {
         pfv19by3->v[ii] *= pfeN;
         pfv19by3->vlogv[ii] *= pfeN;
         pfv19by3->vlogvsq[ii] *= pfeN;
     }

     for(int ii = 0; ii <= PN_PHASING_SERIES_MAX_ORDER; ii++)
     {
         pfv25by3->v[ii] *= pfeN;
         pfv25by3->vlogv[ii] *= pfeN;
         pfv25by3->vlogvsq[ii] *= pfeN;
     }

     for(int ii = 0; ii <= PN_PHASING_SERIES_MAX_ORDER; ii++)
     {
         pfv28by3->v[ii] *= pfeN;
         pfv28by3->vlogv[ii] *= pfeN;
         pfv28by3->vlogvsq[ii] *= pfeN;
     }

     for(int ii = 0; ii <= PN_PHASING_SERIES_MAX_ORDER; ii++)
     {
         pfv31by3->v[ii] *= pfeN;
         pfv31by3->vlogv[ii] *= pfeN;
         pfv31by3->vlogvsq[ii] *= pfeN;
     }

     for(int ii = 0; ii <= PN_PHASING_SERIES_MAX_ORDER; ii++)
     {
         pfv34by3->v[ii] *= pfeN;
         pfv34by3->vlogv[ii] *= pfeN;
         pfv34by3->vlogvsq[ii] *= pfeN;
     }

     for(int ii = 0; ii <= PN_PHASING_SERIES_MAX_ORDER; ii++)
     {
         pfv37by3->v[ii] *= pfeN;
         pfv37by3->vlogv[ii] *= pfeN;
         pfv37by3->vlogvsq[ii] *= pfeN;
     }
}

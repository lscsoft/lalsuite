/*
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

#ifndef _LIGOMETADATARINGDOWNUTILS_H
#define _LIGOMETADATARINGDOWNUTILS_H

#ifdef  __cplusplus
extern "C" {
#endif

/*
 *
 *  ringdown specific structures
 *
 */

typedef enum
taglalringdown_inject_type
{
  LALRINGDOWN_RING_INJECT,
  LALRINGDOWN_IMR_INJECT,
  LALRINGDOWN_IMR_RING_INJECT,
  LALRINGDOWN_EOBNR_INJECT,
  LALRINGDOWN_PHENOM_INJECT
}
lalringdown_inject_type;

typedef enum
taglalringdown_spectrum_type
{
  LALRINGDOWN_SPECTRUM_MEDIAN,
  LALRINGDOWN_SPECTRUM_MEDIAN_MEAN
}
lalringdown_spectrum_type;

typedef enum
taglalringdown_data_type
{
  LALRINGDOWN_DATATYPE_SIM,
  LALRINGDOWN_DATATYPE_ZERO,
  LALRINGDOWN_DATATYPE_UNCAL,
  LALRINGDOWN_DATATYPE_HT_REAL4,
  LALRINGDOWN_DATATYPE_HT_REAL8
}
lalringdown_data_type;

#ifdef  __cplusplus
}
#endif

#endif /* _LIGOMETADATARINGDOWNUTILS_H */

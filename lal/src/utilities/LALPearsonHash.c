/*
 *  Copyright (C) 2016 Karl Wette
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

#include <lal/LALHashFunc.h>

int XLALPearsonHash(
  void *hval,
  const size_t hval_len,
  const void *data,
  const size_t data_len
)
{
  XLAL_CHECK(hval != NULL, XLAL_EFAULT);
  XLAL_CHECK(hval_len > 0, XLAL_EINVAL);
  XLAL_CHECK(data != NULL, XLAL_EFAULT);
  XLAL_CHECK(data_len > 0, XLAL_EINVAL);

  /* Compute the hash value from the given data using a variant of
     Pearson hashing: https://en.wikipedia.org/wiki/Pearson_hashing */
  const unsigned char T[256] = {
    171,38,41,159,93,52,71,169,246,73,191,232,196,33,13,34,
    31,143,228,59,96,63,91,23,250,36,221,126,214,57,90,94,
    157,58,88,213,29,4,83,179,176,121,77,0,226,25,72,15,
    102,76,153,10,137,134,254,20,128,65,198,89,229,47,160,12,
    42,189,203,141,100,40,53,78,182,35,130,68,197,212,55,111,
    18,194,131,252,22,125,147,124,39,11,21,174,249,97,209,144,
    218,27,82,5,220,67,129,150,92,193,48,45,139,216,255,110,
    140,156,192,87,215,69,185,104,175,181,109,75,231,138,180,105,
    6,80,54,190,135,227,50,164,146,8,148,115,123,32,206,7,
    238,74,204,223,177,248,149,188,230,239,170,51,61,30,24,98,
    28,1,37,200,46,184,244,113,116,154,84,86,199,112,133,107,
    145,207,241,166,162,101,172,208,205,236,211,106,19,132,225,56,
    108,99,186,9,251,183,210,114,245,242,187,62,165,60,127,17,
    161,64,3,234,237,16,168,152,120,95,167,155,49,219,81,202,
    122,26,136,44,70,224,79,117,201,253,222,103,233,2,243,178,
    235,240,66,119,173,195,163,85,247,151,158,142,43,14,217,118,
  };
  unsigned char *h = (unsigned char *) hval;
  const unsigned char *x = (const unsigned char *) data;
  for (size_t j = 0; j < hval_len; ++j) {
    for (size_t i = 0; i < data_len; ++i) {
      h[j] = h[j] ^ (x[i] + j);
      h[j] = T[ h[j] ];
    }
  }

  return XLAL_SUCCESS;

}

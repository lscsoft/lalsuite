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

#include <stdlib.h>
#include <stdint.h>
#include <lal/LALHashFunc.h>

static int test_XLALPearsonHash(size_t hval_len, const int N, const int nc_ref)
{
  int n = 0, nc = 0;
  for (int i = 0; i < N; ++i) {
    UINT8 hi = 0;
    XLAL_CHECK(XLALPearsonHash(&hi, hval_len, &i, sizeof(i)) == XLAL_SUCCESS, XLAL_EFUNC);
    for (int j = 0; j <= i; ++j) {
      UINT8 hj = 0;
      XLAL_CHECK(XLALPearsonHash(&hj, hval_len, &j, sizeof(j)) == XLAL_SUCCESS, XLAL_EFUNC);
      XLAL_CHECK(i != j || hi == hj, XLAL_EFAILED);
      ++n;
      if (i != j && hi == hj) {
        ++nc;
      }
    }
  }
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j <= i; ++j) {
      UINT8 hij = 0, hji = 0;
      XLAL_CHECK(XLALPearsonHash(&hij, hval_len, &i, sizeof(i)) == XLAL_SUCCESS, XLAL_EFUNC);
      XLAL_CHECK(XLALPearsonHash(&hij, hval_len, &j, sizeof(j)) == XLAL_SUCCESS, XLAL_EFUNC);
      XLAL_CHECK(XLALPearsonHash(&hji, hval_len, &j, sizeof(j)) == XLAL_SUCCESS, XLAL_EFUNC);
      XLAL_CHECK(XLALPearsonHash(&hji, hval_len, &i, sizeof(i)) == XLAL_SUCCESS, XLAL_EFUNC);
      XLAL_CHECK(i != j || hij == hji, XLAL_EFAILED);
      ++n;
      if (i != j && hij == hji) {
        ++nc;
      }
    }
  }
  printf("XLALPearsonHash(hval=%zu) collisions in %i integer pairs: %i\n", hval_len, n, nc);
  XLAL_CHECK(nc == nc_ref, XLAL_EFAILED, "nc = %i != %i = nc_ref", nc, nc_ref);
  return XLAL_SUCCESS;
}

static int test_XLALCityHash(const char N, const int nc_ref)
{
  int n = 0, nc = 0;
  char buf1[5] = {0}, buf2[5] = {0};
  for (buf1[0] = 'a'; buf1[0] < 'a' + N; ++buf1[0]) {
    for (buf1[1] = 'a'; buf1[1] < 'a' + N; ++buf1[1]) {
      for (buf1[2] = 'a'; buf1[2] < 'a' + N; ++buf1[2]) {
        for (buf1[3] = 'a'; buf1[3] < 'a' + N; ++buf1[3]) {
          const UINT8 h1 = XLALCityHash64(buf1, strlen(buf1));
          for (buf2[0] = 'a'; buf2[0] < 'a' + N; ++buf2[0]) {
            for (buf2[1] = 'a'; buf2[1] < 'a' + N; ++buf2[1]) {
              for (buf2[2] = 'a'; buf2[2] < 'a' + N; ++buf2[2]) {
                for (buf2[3] = 'a'; buf2[3] < 'a' + N; ++buf2[3]) {
                  const UINT8 h2 = XLALCityHash64(buf2, strlen(buf2));
                  ++n;
                  if (strcmp(buf1, buf2) != 0 && h1 == h2) {
                    ++nc;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  printf("XLALCityHash() collisions in %i strings with %i letters: %i\n", n, N, nc);
  XLAL_CHECK(nc == nc_ref, XLAL_EFAILED, "nc = %i != %i = nc_ref", nc, nc_ref);
  return XLAL_SUCCESS;
}

int main(void)
{

  /* Test XLALPearsonHash() */
  XLAL_CHECK_MAIN(test_XLALPearsonHash(sizeof(UINT2), 1024, 8) == XLAL_SUCCESS, XLAL_EFAILED);
  XLAL_CHECK_MAIN(test_XLALPearsonHash(sizeof(UINT4), 1024, 0) == XLAL_SUCCESS, XLAL_EFAILED);
  XLAL_CHECK_MAIN(test_XLALPearsonHash(sizeof(UINT8), 1024, 0) == XLAL_SUCCESS, XLAL_EFAILED);

  /* Test XLALCityHash() */
  XLAL_CHECK_MAIN(test_XLALCityHash(5, 0) == XLAL_SUCCESS, XLAL_EFAILED);
  XLAL_CHECK_MAIN(test_XLALCityHash(7, 0) == XLAL_SUCCESS, XLAL_EFAILED);

  /* Check for memory leaks */
  LALCheckMemoryLeaks();

  return EXIT_SUCCESS;

}

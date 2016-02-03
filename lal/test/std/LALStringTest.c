/*
 *  Copyright (C) 2015 Karl Wette
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

#include <string.h>

#include <lal/LALString.h>
#include <lal/XLALError.h>
#include <lal/LALMalloc.h>

static inline int sign(int s) {
  if (s < 0) return -1;
  if (s > 0) return 1;
  return 0;
}

#define STRCMPCHK(XLALCMP, XLALS1, XLALS2, CCMP, CS1, CS2) do { \
    int xlalsign = sign(XLALCMP(XLALS1, XLALS2)); \
    int csign = sign(CCMP(CS1, CS2)); \
    XLAL_CHECK_MAIN( xlalsign == csign, XLAL_EFAILED, \
                     "Failed string comparison test: sign(%s('%s','%s')) = %i, sign(%s('%s','%s')) = %i\n", \
                     #XLALCMP, XLALS1, XLALS2, xlalsign, #CCMP, CS1, CS2, csign ); \
  } while(0)

#define STRNCMPCHK(N, XLALCMP, XLALS1, XLALS2, CCMP, CS1, CS2) do { \
    int xlalsign = sign(XLALCMP(XLALS1, XLALS2, N)); \
    int csign = sign(CCMP(CS1, CS2, N)); \
    XLAL_CHECK_MAIN( xlalsign == csign, XLAL_EFAILED, \
                     "Failed string comparison test: sign(%s('%s','%s',%zu)) = %i, sign(%s('%s','%s',%zu)) = %i\n", \
                     #XLALCMP, XLALS1, XLALS2, N, xlalsign, #CCMP, CS1, CS2, N, csign ); \
  } while(0)

int main( void )
{

  {
    char *s = NULL;
    XLAL_CHECK( strcmp( (s = XLALStringAppendFmt(s, "Hi %s!", "there")), "Hi there!" ) == 0, XLAL_EFAILED );
    XLAL_CHECK( strcmp( (s = XLALStringAppendFmt(s, " %i + %i = %i?", 2, 2, 4)), "Hi there! 2 + 2 = 4?" ) == 0, XLAL_EFAILED );
    XLALFree( s );
  }

  STRCMPCHK( XLALStringCaseCompare, "0gS1BGPuCG",   "aMYKHogF0r",   strcmp, "0gs1bgpucg",   "amykhogf0r"   );
  STRCMPCHK( XLALStringCaseCompare, "1XJ7G7dIjZ",   "ERB8R8FljJ",   strcmp, "1xj7g7dijz",   "erb8r8fljj"   );
  STRCMPCHK( XLALStringCaseCompare, "8n4EwTImL",    "kmuxDObvfI29", strcmp, "8n4ewtiml",    "kmuxdobvfi29" );
  STRCMPCHK( XLALStringCaseCompare, "bekJxzChXC",   "O3doTwMI7C",   strcmp, "bekjxzchxc",   "o3dotwmi7c"   );
  STRCMPCHK( XLALStringCaseCompare, "f42TpVwldV",   "F5qNN",        strcmp, "f42tpvwldv",   "f5qnn"        );
  STRCMPCHK( XLALStringCaseCompare, "FPpc9",        "w0uwzFMYnd",   strcmp, "fppc9",        "w0uwzfmynd"   );
  STRCMPCHK( XLALStringCaseCompare, "J4io4HKa62",   "s3erTmariX",   strcmp, "j4io4hka62",   "s3ertmarix"   );
  STRCMPCHK( XLALStringCaseCompare, "kdi7f2",       "aSiM7gvIm",    strcmp, "kdi7f2",       "asim7gvim"    );
  STRCMPCHK( XLALStringCaseCompare, "nwWd7qaipf",   "ZnDeR",        strcmp, "nwwd7qaipf",   "znder"        );
  STRCMPCHK( XLALStringCaseCompare, "oVTTSi",       "tG1",          strcmp, "ovttsi",       "tg1"          );
  STRCMPCHK( XLALStringCaseCompare, "t6fv04Bjc5",   "6tnnoIDZHY",   strcmp, "t6fv04bjc5",   "6tnnoidzhy"   );
  STRCMPCHK( XLALStringCaseCompare, "tW8Ng63utwwS", "x47v59r9e",    strcmp, "tw8ng63utwws", "x47v59r9e"    );
  STRCMPCHK( XLALStringCaseCompare, "yiMB0wmdkq",   "X1BpBeOFwy",   strcmp, "yimb0wmdkq",   "x1bpbeofwy"   );
  STRCMPCHK( XLALStringCaseCompare, "zaNDp",        "Kd0H0AVtaA",   strcmp, "zandp",        "kd0h0avtaa"   );
  STRCMPCHK( XLALStringCaseCompare, "ztZDUBu4rviQ", "gzM8YDLnzedn", strcmp, "ztzdubu4rviq", "gzm8ydlnzedn" );

  for (size_t n = 0; n < 12; ++n) {

    STRNCMPCHK( n, XLALStringNCaseCompare, "0gS1BGPuCG",   "aMYKHogF0r",   strncmp, "0gs1bgpucg",   "amykhogf0r"   );
    STRNCMPCHK( n, XLALStringNCaseCompare, "1XJ7G7dIjZ",   "ERB8R8FljJ",   strncmp, "1xj7g7dijz",   "erb8r8fljj"   );
    STRNCMPCHK( n, XLALStringNCaseCompare, "8n4EwTImL",    "kmuxDObvfI29", strncmp, "8n4ewtiml",    "kmuxdobvfi29" );
    STRNCMPCHK( n, XLALStringNCaseCompare, "bekJxzChXC",   "O3doTwMI7C",   strncmp, "bekjxzchxc",   "o3dotwmi7c"   );
    STRNCMPCHK( n, XLALStringNCaseCompare, "f42TpVwldV",   "F5qNN",        strncmp, "f42tpvwldv",   "f5qnn"        );
    STRNCMPCHK( n, XLALStringNCaseCompare, "FPpc9",        "w0uwzFMYnd",   strncmp, "fppc9",        "w0uwzfmynd"   );
    STRNCMPCHK( n, XLALStringNCaseCompare, "J4io4HKa62",   "s3erTmariX",   strncmp, "j4io4hka62",   "s3ertmarix"   );
    STRNCMPCHK( n, XLALStringNCaseCompare, "kdi7f2",       "aSiM7gvIm",    strncmp, "kdi7f2",       "asim7gvim"    );
    STRNCMPCHK( n, XLALStringNCaseCompare, "nwWd7qaipf",   "ZnDeR",        strncmp, "nwwd7qaipf",   "znder"        );
    STRNCMPCHK( n, XLALStringNCaseCompare, "oVTTSi",       "tG1",          strncmp, "ovttsi",       "tg1"          );
    STRNCMPCHK( n, XLALStringNCaseCompare, "t6fv04Bjc5",   "6tnnoIDZHY",   strncmp, "t6fv04bjc5",   "6tnnoidzhy"   );
    STRNCMPCHK( n, XLALStringNCaseCompare, "tW8Ng63utwwS", "x47v59r9e",    strncmp, "tw8ng63utwws", "x47v59r9e"    );
    STRNCMPCHK( n, XLALStringNCaseCompare, "yiMB0wmdkq",   "X1BpBeOFwy",   strncmp, "yimb0wmdkq",   "x1bpbeofwy"   );
    STRNCMPCHK( n, XLALStringNCaseCompare, "zaNDp",        "Kd0H0AVtaA",   strncmp, "zandp",        "kd0h0avtaa"   );
    STRNCMPCHK( n, XLALStringNCaseCompare, "ztZDUBu4rviQ", "gzM8YDLnzedn", strncmp, "ztzdubu4rviq", "gzm8ydlnzedn" );

  }

  {
    char s[] = "abc,def-,ghij k-lmn, ,-";
    char *p = s;
    XLAL_CHECK_MAIN( strcmp( XLALStringToken(&p, ",-", 0), "abc" ) == 0, XLAL_EFAILED );
    XLAL_CHECK_MAIN( strcmp( XLALStringToken(&p, ",-", 0), "def" ) == 0, XLAL_EFAILED );
    XLAL_CHECK_MAIN( strcmp( XLALStringToken(&p, ",-", 0), "ghij k" ) == 0, XLAL_EFAILED );
    XLAL_CHECK_MAIN( strcmp( XLALStringToken(&p, ",-", 0), "lmn" ) == 0, XLAL_EFAILED );
    XLAL_CHECK_MAIN( strcmp( XLALStringToken(&p, ",-", 0), " " ) == 0, XLAL_EFAILED );
    XLAL_CHECK_MAIN( XLALStringToken(&p, ",-", 0) == NULL, XLAL_EFAILED );
  }

  {
    char s[] = "abc,def-,ghij k-lmn, ,-";
    char *p = s;
    XLAL_CHECK_MAIN( strcmp( XLALStringToken(&p, ",-", 1), "abc" ) == 0, XLAL_EFAILED );
    XLAL_CHECK_MAIN( strcmp( XLALStringToken(&p, ",-", 1), "def" ) == 0, XLAL_EFAILED );
    XLAL_CHECK_MAIN( strcmp( XLALStringToken(&p, ",-", 1), "" ) == 0, XLAL_EFAILED );
    XLAL_CHECK_MAIN( strcmp( XLALStringToken(&p, ",-", 1), "ghij k" ) == 0, XLAL_EFAILED );
    XLAL_CHECK_MAIN( strcmp( XLALStringToken(&p, ",-", 1), "lmn" ) == 0, XLAL_EFAILED );
    XLAL_CHECK_MAIN( strcmp( XLALStringToken(&p, ",-", 1), " " ) == 0, XLAL_EFAILED );
    XLAL_CHECK_MAIN( strcmp( XLALStringToken(&p, ",-", 1), "" ) == 0, XLAL_EFAILED );
    XLAL_CHECK_MAIN( strcmp( XLALStringToken(&p, ",-", 1), "" ) == 0, XLAL_EFAILED );
    XLAL_CHECK_MAIN( XLALStringToken(&p, ",-", 1) == NULL, XLAL_EFAILED );
  }

  {
    char s[] = "abcbdeacd";
    XLAL_CHECK_MAIN( strcmp( XLALStringReplaceChar(s, 'a', 'b'), "bbcbdebcd" ) == 0, XLAL_EFAILED );
    XLAL_CHECK_MAIN( strcmp( XLALStringReplaceChar(s, 'b', 'c'), "ccccdeccd" ) == 0, XLAL_EFAILED );
    XLAL_CHECK_MAIN( strcmp( XLALStringReplaceChar(s, 'c', 'd'), "dddddeddd" ) == 0, XLAL_EFAILED );
    XLAL_CHECK_MAIN( strcmp( XLALStringReplaceChar(s, 'd', 'e'), "eeeeeeeee" ) == 0, XLAL_EFAILED );
  }

  LALCheckMemoryLeaks();

  return 0;

}

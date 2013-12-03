/*
*  Copyright (C) 2013 Jolien Creighton, Kipp Cannon
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
#include <ctype.h>
#include <string.h>
#include <config.h>
#include <lal/LALDebugLevel.h>
#include <lal/LALError.h>
#undef lalDebugLevel

#ifdef LAL_PTHREAD_LOCK
#include <pthread.h>
static pthread_once_t lalOnce = PTHREAD_ONCE_INIT;
#define LAL_ONCE(init) pthread_once(&lalOnce, (init))
#else
static int lalOnce = 1;
#define LAL_ONCE(init) (lalOnce ? (init)(), lalOnce = 0 : 0)
#endif

#if defined(NDEBUG) || defined(LAL_NDEBUG)
#define LAL_DEFAULT_DEBUG_LEVEL 0
#else
#define LAL_DEFAULT_DEBUG_LEVEL 1
#endif

static int lalDebugLevel = LAL_DEFAULT_DEBUG_LEVEL;

/*
 * Bad function, put here only for people who know what they are doing,
 * (and understand that what they are doing is wrong).
 */
void XLALClobberDebugLevel(int level);
void XLALClobberDebugLevel(int level)
{
    lalDebugLevel = level;
}

static int compare_token_to_name(const char* token, const size_t toklen, const char* name) {
  if (toklen != strlen(name)) {
    return 0;
  }
  for (size_t i = 0; i < toklen; ++i) {
    if (toupper(token[i]) != name[i]) {
      return 0;
    }
  }
  return 1;
}

static void XLALSetDebugLevel(void)
{
    int level;
    char *end;
    const char *env = getenv("LAL_DEBUG_LEVEL");
    if (env == NULL || *env == '\0')
        return;

    /* parse the environment string */

    /* see if it is an integer */
    level = (int) strtol(env, &end, 0);
    if (*end == '\0') { /* conversion successful */
        lalDebugLevel = level;
        return;
    }

    /* try to parse as a comma-separated list of debug level names */
    const char *const seps = ",";
    const char *token = env;
    do {
        size_t toklen = strcspn(token, seps);
        if (toklen > 0) {
            if (compare_token_to_name(token, toklen, "NDEBUG")) {
                level |= 0; /* no debugging */
            } else if (compare_token_to_name(token, toklen, "ERROR")) {
                level |= LALERRORBIT; /* enable error messages */
            } else if (compare_token_to_name(token, toklen, "WARNING")) {
                level |= LALWARNINGBIT; /* enable warning messages */
            } else if (compare_token_to_name(token, toklen, "INFO")) {
                level |= LALINFOBIT; /* enable info messages */
            } else if (compare_token_to_name(token, toklen, "TRACE")) {
                level |= LALTRACEBIT; /* enable tracing messages */
            } else if (compare_token_to_name(token, toklen, "MSGLVL1")) {
                level |= LALERRORBIT; /* enable error messages */
            } else if (compare_token_to_name(token, toklen, "MSGLVL2")) {
                level |= LALERRORBIT | LALWARNINGBIT; /* enable error and warning messages */
            } else if (compare_token_to_name(token, toklen, "MSGLVL3")) {
                level |= LALERRORBIT | LALWARNINGBIT | LALINFOBIT; /* enable error, warning, and info messages */
            } else if (compare_token_to_name(token, toklen, "MEMDBG")) {
                level |= LALMEMDBGBIT | LALMEMPADBIT | LALMEMTRKBIT; /* enable memory debugging tools */
            } else if (compare_token_to_name(token, toklen, "MEMTRACE")) {
                level |= LALTRACEBIT | LALMEMDBG | LALMEMINFOBIT; /* enable memory tracing tools */
            } else if (compare_token_to_name(token, toklen, "ALLDBG")) {
                level |= ~LALNDEBUG; /* enable all debugging */
            } else {
              /* not a valid debug level name */
                lalAbortHook("%s: could not parse LAL_DEBUG_LEVEL='%s'\n", __func__, env);
                return;
            }
        }
        token += toklen;
    } while (*(token++) != '\0');
    lalDebugLevel = level;
    return;

}

int XLALGetDebugLevel(void)
{
    LAL_ONCE(XLALSetDebugLevel);
    return lalDebugLevel;
}

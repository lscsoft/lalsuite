/*
*  Copyright (C) 2013, 2015 Karl Wette
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
#include <lal/LALString.h>
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

/* LAL_DEFAULT_DEBUG_LEVEL is defined in <config.h> by LAL_WITH_DEFAULT_DEBUG_LEVEL() in gnuscripts/lal.m4 */
static int lalDebugLevel = LAL_DEFAULT_DEBUG_LEVEL;

/*
 * Bad function, put here only for people who know what they are doing,
 * (and understand that what they are doing is wrong).
 */
void XLALClobberDebugLevel(int level)
{
    lalDebugLevel = level;
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
            if (XLALStringNCaseCompare("NDEBUG", token, toklen) == 0) {
                level |= 0; /* no debugging */
            } else if (XLALStringNCaseCompare("ERROR", token, toklen) == 0) {
                level |= LALERRORBIT; /* enable error messages */
            } else if (XLALStringNCaseCompare("WARNING", token, toklen) == 0) {
                level |= LALWARNINGBIT; /* enable warning messages */
            } else if (XLALStringNCaseCompare("INFO", token, toklen) == 0) {
                level |= LALINFOBIT; /* enable info messages */
            } else if (XLALStringNCaseCompare("TRACE", token, toklen) == 0) {
                level |= LALTRACEBIT; /* enable tracing messages */
            } else if (XLALStringNCaseCompare("MSGLVL1", token, toklen) == 0) {
                level |= LALERRORBIT; /* enable error messages */
            } else if (XLALStringNCaseCompare("MSGLVL2", token, toklen) == 0) {
                level |= LALERRORBIT | LALWARNINGBIT; /* enable error and warning messages */
            } else if (XLALStringNCaseCompare("MSGLVL3", token, toklen) == 0) {
                level |= LALERRORBIT | LALWARNINGBIT | LALINFOBIT; /* enable error, warning, and info messages */
            } else if (XLALStringNCaseCompare("MEMDBG", token, toklen) == 0) {
                level |= LALMEMDBGBIT | LALMEMPADBIT | LALMEMTRKBIT; /* enable memory debugging tools */
            } else if (XLALStringNCaseCompare("MEMTRACE", token, toklen) == 0) {
                level |= LALTRACEBIT | LALMEMDBG | LALMEMINFOBIT; /* enable memory tracing tools */
            } else if (XLALStringNCaseCompare("ALLDBG", token, toklen) == 0) {
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

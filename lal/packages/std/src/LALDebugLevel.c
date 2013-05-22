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
#include <config.h>
#include <lal/LALDebugLevel.h>
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

    /* could handle more complicated strings in the future? */
    /* for now, ignore environment strings that cannot be parsed */
    return;
}

int XLALGetDebugLevel(void)
{
    LAL_ONCE(XLALSetDebugLevel);
    return lalDebugLevel;
}

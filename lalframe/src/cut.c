/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Robert Adam Mercer
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALFrameU.h>
#include "utils.h"

#define FAILURE(...) do { fprintf(stderr, __VA_ARGS__); exit(99); } while (0)

char **chanvalloc(const char *str);
void chanvfree(char **chanv);
char *chanvstr(char **chanv, const char *str);

int main(int argc, char *argv[])
{
    char stdio[] = "-";
    char *defaultfilev[1] = { stdio };
    int filec = (argc == 2 ? 1 : argc - 2);
    char **filev = (argc == 2 ? defaultfilev : &argv[2]);
    char **chanv;
    int f;
    LALFrameUFrFile *output;

    if (argc == 1) {
        fprintf(stderr, "usage: %s list [file ...]\n", argv[0]);
        return 1;
    }

    if (argc == 2) {
        filec = 1;
        filev = defaultfilev;
    } else {
        filec = argc - 2;
        filev = &argv[2];
    }

    /* convert list of channel names to a vector of channel names */
    chanv = chanvalloc(argv[1]);

    output = XLALFrameUFrFileOpen(NULL, "w");
    if (!output)
        FAILURE("could not create output file\n");

    for (f = 0; f < filec; ++f) {
        char *fname = filev[f];
        LALFrameUFrFile *input;
        LALFrameUFrTOC *toc;
        size_t nframe;
        size_t pos;

        input = XLALFrameUFrFileOpen(fname, "r");
        if (!input)
            FAILURE("file %s not found\n", fname);

        toc = XLALFrameUFrTOCRead(input);
        if (!toc)
            FAILURE("no TOC found in file %s\n", fname);

        nframe = XLALFrameUFrTOCQueryNFrame(toc);
        if ((int)(nframe) <= 0)
            FAILURE("no frames found in file %s\n", fname);

        /* loop over frames in input file */
        for (pos = 0; pos < nframe; ++pos) {
            LALFrameUFrameH *frame;
            char **chanp;

            frame = framecpy(input, pos);
            copydetectors(frame, input);
            for (chanp = chanv; *chanp; ++chanp)
                copychannels(frame, input, pos, *chanp);
            XLALFrameUFrameHWrite(output, frame);
            XLALFrameUFrameHFree(frame);
        }

        XLALFrameUFrFileClose(input);
    }

    XLALFrameUFrFileClose(output);
    chanvfree(chanv);
    return 0;
}

/* creates a vector of channel names from a comma-delimited channel list */
char **chanvalloc(const char *str)
{
    char *s = strdup(str);
    char **chanv;
    int chanc = 1;
    int i;
    for (; (str = strchr(str, ',')); ++str)
        ++chanc;
    chanv = calloc(chanc + 1, sizeof(*chanv));
    if (chanc > 1) {
        chanv[0] = strdup(strtok(s, ","));
        for (i = 1; i < chanc; ++i)
            chanv[i] = strdup(strtok(NULL, ","));
        free(s);
    } else
        chanv[0] = s;
    return chanv;
}

/* frees a vector of channel names */
void chanvfree(char **chanv)
{
    if (chanv) {
        char **p;
        for (p = chanv; *p; ++p)
            free(*p);
        free(chanv);
    }
    return;
}

/* seeks the string str in a vector of channel names;
 * returns NULL if not found */
char *chanvstr(char **chanv, const char *str)
{
    while (*chanv && strcmp(*chanv, str))
        ++chanv;
    return *chanv;
}

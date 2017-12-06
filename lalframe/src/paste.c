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

int main(int argc, char *argv[])
{
    char stdio[] = "-";
    char *defaultfilev[1] = { stdio };
    int filec = (argc == 1 ? 1 : argc - 1);
    char **filev = (argc == 1 ? defaultfilev : &argv[1]);
    LALFrameUFrFile *output;
    LALFrameUFrFile **input;
    size_t *nframe;
    size_t nframemax = 0;
    size_t pos;
    int f;

    output = XLALFrameUFrFileOpen(stdio, "w");  /* output to stdout */
    if (!output)
        FAILURE("could not create output file\n");

    /* open all the input files and get number of frames in each */
    input = calloc(filec, sizeof(*input));
    nframe = calloc(filec, sizeof(*nframe));
    for (f = 0; f < filec; ++f) {
        LALFrameUFrTOC *toc;
        char *fname = filev[f];

        input[f] = XLALFrameUFrFileOpen(fname, "r");
        if (!input[f])
            FAILURE("file %s not found\n", fname);

        toc = XLALFrameUFrTOCRead(input[f]);
        if (!toc)
            FAILURE("no TOC found in file %s\n", fname);

        nframe[f] = XLALFrameUFrTOCQueryNFrame(toc);
        if ((int)(nframe[f]) <= 0)
            FAILURE("no frames found in file %s\n", fname);

        if (nframe[f] > nframemax)
            nframemax = nframe[f];
    }

    /* loop over frames */
    for (pos = 0; pos < nframemax; ++pos) {
        LALFrameUFrameH *frame = NULL;

        /* loop over files */
        for (f = 0; f < filec; ++f) {
            /* check to see that file has this frame */
            if (pos >= nframe[f])
                continue;       /* skip this file */

            /* TODO: check consistency of frame times? */
            if (frame == NULL)  /* first time: duplicate frame */
                frame = framecpy(input[f], pos);

            copydetectors(frame, input[f]);
            copychannels(frame, input[f], pos, NULL);
        }
        XLALFrameUFrameHWrite(output, frame);
        XLALFrameUFrameHFree(frame);
    }

    /* close files */
    for (f = 0; f < filec; ++f)
        XLALFrameUFrFileClose(input[f]);
    free(input);
    free(nframe);

    XLALFrameUFrFileClose(output);
    return 0;
}

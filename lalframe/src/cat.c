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

/**
 * @defgroup lalfr_cat lalfr-cat
 * @ingroup lalframe_programs
 *
 * @brief Concatenate frame files
 *
 * ### Synopsis
 *
 *     lalfr-cat [file ...]
 * 
 * ### Description
 * 
 * The `lalfr-cat` utility cuts reads frame files sequentially, writing them to
 * the standard output.  The file operands are  processed  in  command-line
 * order.  If file is a single dash (`-`) or absent, `lalfr-cat` reads from the
 * standard input.
 *
 * ### Example
 *
 * The command:
 *
 *     lalfr-cat file1.gwf file2.gwf > file3.gwf
 *
 * will concatenate `file1.gwf` and `file2.gwf` to the file `file3.gwf`.
 *
 * The command:
 *
 *     lalfr-cat file1.gwf - file2.gwf
 *
 * will output to standard output the contents of frame file `file1.gwf`, the
 * frame file received from standard input, and the contents of `file2.gwf`.
 * 
 * @sa @ref lalfr_split
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

    output = XLALFrameUFrFileOpen(stdio, "w");
    if (!output)
        FAILURE("could not create output file\n");

    while (filec-- > 0) {
        char *fname = *filev++;
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

            frame = framecpy(input, pos);
            copydetectors(frame, input);
            copychannels(frame, input, pos, NULL);

            XLALFrameUFrameHWrite(output, frame);
            XLALFrameUFrameHFree(frame);
        }

        XLALFrameUFrFileClose(input);
    }

    XLALFrameUFrFileClose(output);
    return 0;
}

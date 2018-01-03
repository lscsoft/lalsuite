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
 * @defgroup lalfr_split lalfr-split
 * @ingroup lalframe_programs
 *
 * @brief Split a frame file into component frames
 *
 * ### Synopsis
 *
 *     lalfr-split [file [name]]
 *     
 * ### Description
 *
 * The `lalfr-split` utility reads the given `file` containing multiple frames
 * and breaks it up into files each containing one frame.  If `file` is a
 * single dash (`-`) or absent, `lalfr-split` reads from the standard input.
 * If additional arguments are specified, the first is used as the name of the
 * input  file which is to be split.  If a second additional argument is
 * specified, it is used as a prefix for the names of  the  files  into which
 * the  file is split.  In this case, each file into which the file is split is
 * named by the prefix followed by a lexically ordered  suffix using  two
 * characters in the range `a-z` followed by the extension `.gwf`.  If the name
 * argument is not specified, the file is split into lexically ordered files
 * named with the prefix `x` and with suffixes as above.
 *
 * ### Examples
 *
 * If `file.gwf` contains three frames then the command:
 *
 *     lalfr-split file.gwf
 *
 * produces three files, `xaa.gwf`, `xab.gwf`, and `xac.gwf`, each containing
 * one of the three frames.
 *
 * @sa @ref lalfr_cat
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALFrameU.h>
#include "utils.h"

#define FAILURE(...) do { fprintf(stderr, __VA_ARGS__); exit(99); } while (0)

char *mkfname(const char *name, size_t pos);

int main(int argc, char *argv[])
{
    char stdio[] = "-";
    char *fname = NULL;
    char *basename = NULL;
    LALFrameUFrFile *input;
    LALFrameUFrTOC *toc;
    size_t nframe;
    size_t pos;

    if (argc == 1)      /* read from stdin */
        fname = stdio;
    if (argc == 2)      /* read from argv[1] */
        fname = argv[1];
    if (argc == 3) {    /* read from argv[1] with output prefix argv[2] */
        fname = argv[1];
        basename = argv[2];
    }
    if (argc > 3) {
        fprintf(stderr, "usage: %s [file [name]]\n", argv[0]);
        return 1;
    }

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
        LALFrameUFrFile *output;
        LALFrameUFrameH *frame;

        fname = mkfname(basename, pos);
        output = XLALFrameUFrFileOpen(fname, "w");
        frame = framecpy(input, pos);

        copydetectors(frame, input);
        copychannels(frame, input, pos, NULL);

        XLALFrameUFrameHWrite(output, frame);
        XLALFrameUFrameHFree(frame);
        XLALFrameUFrFileClose(output);
    }

    XLALFrameUFrFileClose(input);
    return 0;
}

char *mkfname(const char *name, size_t pos)
{
    static char fname[FILENAME_MAX];
    const char *ext = "gwf";
    char hi = 'a' + pos / 26;
    char lo = 'a' + pos % 26;
    snprintf(fname, sizeof(fname), "%s%c%c.%s", name ? name : "x", hi, lo,
        ext);
    return fname;
}

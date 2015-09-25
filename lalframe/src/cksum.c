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
 * @defgroup lalfr_cksum lalfr-cksum
 * @ingroup lalframe_programs
 *
 * @brief Validates the checksum of a frame file
 *
 * ### Synopsis
 *
 *     lalfr-cksum file [files ...]
 *
 * ### Description
 *
 * The `lalfr-cksum` utility validates the checksum on the specified files to
 * the standard output.	The file operands are processed in command-line order.
 * If file is a single dash (`-`) or absent, `lalfr-cksum` reads from the
 * standard input.
 *
 * ### Exit Status
 *
 * The `lalfr-cksum` utility exits 0 on success, and >0 if one or more of the
 * frame files have invalid checksums.
 *
 * ### Example
 *
 * The command:
 *
 *     lalfr-cksum file1.gwf file2.gwf
 *
 * will validate frame files `file1.gwf` and `file2.gwf`.
 */

#include <stdio.h>
#include <stdlib.h>
#include <lal/LALFrameU.h>

#define FAILURE(...) do { fprintf(stderr, __VA_ARGS__); exit(99); } while (0)

int main(int argc, char *argv[])
{
    int retval = 0;

    if (argc == 1) {
        fprintf(stderr, "usage: %s framefiles\n", argv[0]);
        return 1;
    }

    while (--argc > 0) {
        char *fname = *++argv;
        LALFrameUFrFile *frfile;
        int valid;

        frfile = XLALFrameUFrFileOpen(fname, "r");
        if (!frfile)
            FAILURE("file %s not found\n", fname);

        valid = XLALFrameUFileCksumValid(frfile);
        retval += !valid;

        fprintf(stdout, "%svalid checksum for %s\n", valid ? "" : "in",
            fname);

        XLALFrameUFrFileClose(frfile);
    }

    return retval;
}

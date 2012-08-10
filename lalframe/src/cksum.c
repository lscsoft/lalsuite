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

        valid = (XLALFrameUFileCksumValid(frfile) == 0);
        retval += !valid;

        fprintf(stdout, "%svalid checksum for %s\n", valid ? "" : "in", fname);

        XLALFrameUFrFileClose(frfile);
    }

    return retval;
}

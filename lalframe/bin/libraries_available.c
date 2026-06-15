/*
 * Copyright (C) 2026 Karl Wette
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA
 */

#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*
 * Print which frame libraries are available for use with LALFrame
 */

int main(int argc, char **argv) {

  /* Parse command line */
  if ( argc > 1 ) {
    fprintf(stderr,
            "Usage: %s [-h|--help]\n\n"
            "Print which frame libraries are available for use with LALFrame\n\n"
            "Set the environment variable LAL_FRAME_LIBRARY to one of the printed options to make LALFrame use that frame library\n\n"
            "Options:\n"
            "  --help       display this message and exit\n",
            argv[0]
      );
    if ( strcmp( argv[1], "-h" ) == 0 || strcmp( argv[1], "-help" ) == 0 || strcmp( argv[1], "--help" ) == 0 ) {
      return EXIT_SUCCESS;
    } else {
      return EXIT_FAILURE;
    }
  }

  /* Print if FrameL support if available */
#if defined HAVE_FRAMEL_H && defined HAVE_LIBFRAMEL
  printf("FrameL\n");
#endif

/* Print if FrameC support if available */
#if defined HAVE_FRAMECPPC_FRAMEC_H && defined HAVE_LIBFRAMECPPC
  printf("FrameC\n");
#endif

  return EXIT_SUCCESS;

}

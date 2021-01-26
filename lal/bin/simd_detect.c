/*
 * Copyright (C) 2015 Reinhard Prix, Karl Wette
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

/*
 * Utility for printing out detected SIMD extensions
 */

#include <stdio.h>
#include <config.h>

#include <lal/LALStdlib.h>
#include <lal/LALSIMD.h>

int main(int argc, char **argv) {

  /* Parse command line */
  if ( argc > 1 ) {
    fprintf(stderr,
            "Usage: %s [-h|--help]\n\n"
            "Print compiled and detected SIMD extensions\n\n"
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

  printf("%s was compiled with support for the following instruction sets:\n   %s %s\n",
         PACKAGE_STRING, XLALSIMDInstructionSetName(0), HAVE_SIMD_COMPILER);

  printf("This machine supports executing the following instruction sets:\n  ");
  for (LAL_SIMD_ISET iset = 0; XLALHaveSIMDInstructionSet(iset); ++iset) {
    printf(" %s", XLALSIMDInstructionSetName(iset));
  }
  printf("\n");

  return 0;

}

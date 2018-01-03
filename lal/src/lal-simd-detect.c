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
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA  02111-1307  USA
 */

/*
 * Utility for printing out detected SIMD extensions
 */

#include <stdio.h>
#include <config.h>

#include <lal/LALSIMD.h>

int main(void) {

  printf("%s was compiled with support for the following instruction sets:\n   %s %s\n",
         PACKAGE_STRING, XLALSIMDInstructionSetName(0), HAVE_SIMD_COMPILER);

  printf("This machine supports executing the following instruction sets:\n  ");
  for (LAL_SIMD_ISET iset = 0; XLALHaveSIMDInstructionSet(iset); ++iset) {
    printf(" %s", XLALSIMDInstructionSetName(iset));
  }
  printf("\n");

  return 0;

}

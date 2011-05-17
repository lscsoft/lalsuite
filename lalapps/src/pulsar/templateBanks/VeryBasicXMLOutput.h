/*
 *  Copyright (C) 2008 Karl Wette
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
 * \author Karl Wette
 * \file
 * \brief Some very basic XML output routines, might
          be useful for debugging or test programs
 */

#ifndef _VERYBASICXMLOUTPUT_H
#define _VERYBASICXMLOUTPUT_H

#include <stdio.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <lal/LALRCSID.h>
#include <lal/LALDatatypes.h>

#ifdef __cplusplus
extern "C" {
#endif

  NRCSID(VERYBASICXMLOUTPUTH, "$Id$");

  /**
   * Information for XML functions
   */
  typedef struct tagVeryBasicXMLOutput {

    /* XML file */
    FILE *file;

    /* Indent level */
    INT4 indent;

    /* Separator */
    CHAR separator;

  } VeryBasicXMLOutput;
  extern const VeryBasicXMLOutput empty_VeryBasicXMLOutput;

  /**
   * Functions
   */
  void XLAL_VBXMLO_Header(VeryBasicXMLOutput*, INT4, INT4);
  void XLAL_VBXMLO_Indent(VeryBasicXMLOutput*);
  void XLAL_VBXMLO_BeginTag(VeryBasicXMLOutput*, const char*);
  void XLAL_VBXMLO_EndTag(VeryBasicXMLOutput*, const char*);
  void XLAL_VBXMLO_Tag(VeryBasicXMLOutput*, const char*, const char*, ...);
  void XLAL_VBXMLO_Printf(VeryBasicXMLOutput*, const char*, ...);
  void XLAL_VBXMLO_gsl_vector(VeryBasicXMLOutput*, const char*, const char*, gsl_vector*);
  void XLAL_VBXMLO_gsl_vector_int(VeryBasicXMLOutput*, const char*, const char*, gsl_vector_int*);
  void XLAL_VBXMLO_gsl_matrix(VeryBasicXMLOutput*, const char*, const char*, gsl_matrix*);
  void XLAL_VBXMLO_gsl_matrix_int(VeryBasicXMLOutput*, const char*, const char*, gsl_matrix_int*);

#ifdef __cplusplus
}
#endif

#endif

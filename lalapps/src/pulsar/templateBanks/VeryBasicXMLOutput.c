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

#include "VeryBasicXMLOutput.h"

#include <stdarg.h>

NRCSID(VERYBASICXMLOUTPUTC, "$Id");

const VeryBasicXMLOutput empty_VeryBasicXMLOutput = {NULL, 0, ' '};

/**
 * Write XML header
 */
void XLAL_VBXMLO_Header(
			VeryBasicXMLOutput *xml, /**< [in] XML structure */
			INT4 major,              /**< [in] Major revision */
			INT4 minor               /**< [in] Minor revision */
			)
{

  if (!xml->file)
    return;

  /* Print header */
  fprintf(xml->file, "<?xml version=\"%i.%i\"?>\n", major, minor);

  /* Set indent */
  xml->indent = 0;

}

/**
 * Print the appropriate indent
 */
void XLAL_VBXMLO_Indent(
			VeryBasicXMLOutput *xml /**< [in] XML structure */
			)
{

  int i;

  if (!xml->file)
    return;

  for (i = 0; i < xml->indent; ++i)
    fprintf(xml->file, "  ");

}

/**
 * Begin an XML tag
 */
void XLAL_VBXMLO_BeginTag(
			  VeryBasicXMLOutput *xml, /**< [in] XML structure */
			  const char *name         /**< [in] Tag name */
			  )
{

  va_list va;

  if (!xml->file)
    return;

  XLAL_VBXMLO_Indent(xml);
  fprintf(xml->file, "<%s>\n", name);
  ++xml->indent;
  va_end(va);

}

/**
 * End an XML tag
 */
void XLAL_VBXMLO_EndTag(
			VeryBasicXMLOutput *xml, /**< [in] XML structure */
			const char *name         /**< [in] Tag name */
			)
{

  if (!xml->file)
    return;

  --xml->indent;
  XLAL_VBXMLO_Indent(xml);
  fprintf(xml->file, "</%s>\n", name);

}

/**
 * Print a complete XML tag
 */
void XLAL_VBXMLO_Tag(
		     VeryBasicXMLOutput *xml, /**< [in] XML structure */
		     const char *name,        /**< [in] Tag name */
		     const char *format,      /**< [in] Contents format */
		     ...                      /**< [in] Attribute and contents arguments */
		     )
{

  va_list va;

  if (!xml->file)
    return;

  XLAL_VBXMLO_BeginTag(xml, name);
  XLAL_VBXMLO_Indent(xml);
  va_start(va, format);
  vfprintf(xml->file, format, va);
  va_end(va);
  fprintf(xml->file, "\n");
  XLAL_VBXMLO_EndTag(xml, name);

}

/**
 * Print a formatted string
 */
void XLAL_VBXMLO_Printf(
			VeryBasicXMLOutput *xml, /**< [in] XML structure */
			const char *format,      /**< [in] Format */
			...                      /**< [in] Arguments */
			)
{

  va_list va;

  if (!xml->file)
    return;

  va_start(va, format);
  vfprintf(xml->file, format, va);
  va_end(va);

}

/**
 * Print a gsl_vector
 */
void XLAL_VBXMLO_gsl_vector(
			    VeryBasicXMLOutput *xml, /**< [in] XML structure */
			    const char *name,        /**< [in] Tag name */
			    const char *format,      /**< [in] Value format */
			    gsl_vector *v            /**< [in] Vector */
			    )
{

  size_t i;

  if (!xml->file)
    return;

  XLAL_VBXMLO_BeginTag(xml, name);
  XLAL_VBXMLO_Indent(xml);
  for (i = 0; i < v->size; ++i) {
    if (i > 0)
      fprintf(xml->file, "%c", xml->separator);
    fprintf(xml->file, format, gsl_vector_get(v, i));
  }
  fprintf(xml->file, "\n");
  XLAL_VBXMLO_EndTag(xml, name);

}

/**
 * Print a gsl_vector_int
 */
void XLAL_VBXMLO_gsl_vector_int(
				VeryBasicXMLOutput *xml, /**< [in] XML structure */
				const char *name,        /**< [in] Tag name */
				const char *format,      /**< [in] Value format */
				gsl_vector_int *v        /**< [in] Vector */
				)
{

  size_t i;

  if (!xml->file)
    return;

  XLAL_VBXMLO_BeginTag(xml, name);
  XLAL_VBXMLO_Indent(xml);
  for (i = 0; i < v->size; ++i) {
    if (i > 0)
      fprintf(xml->file, "%c", xml->separator);
    fprintf(xml->file, format, gsl_vector_int_get(v, i));
  }
  fprintf(xml->file, "\n");
  XLAL_VBXMLO_EndTag(xml, name);

}

/**
 * Print a gsl_matrix
 */
void XLAL_VBXMLO_gsl_matrix(
			    VeryBasicXMLOutput *xml, /**< [in] XML structure */
			    const char *name,        /**< [in] Tag name */
			    const char *format,      /**< [in] Value format */
			    gsl_matrix *m            /**< [in] Matrix */
			    )
{

  size_t i;

  if (!xml->file)
    return;

  XLAL_VBXMLO_BeginTag(xml, name);
  for (i = 0; i < m->size1; ++i) {
    gsl_vector_view row = gsl_matrix_row(m, i);
    XLAL_VBXMLO_gsl_vector(xml, "row", format, &row.vector);
  }
  XLAL_VBXMLO_EndTag(xml, name);

}

/**
 * Print a gsl_matrix_int
 */
void XLAL_VBXMLO_gsl_matrix_int(
				VeryBasicXMLOutput *xml, /**< [in] XML structure */
				const char *name,        /**< [in] Tag name */
				const char *format,      /**< [in] Value format */
				gsl_matrix_int *m        /**< [in] Matrix */
				)
{

  size_t i;

  if (!xml->file)
    return;

  XLAL_VBXMLO_BeginTag(xml, name);
  for (i = 0; i < m->size1; ++i) {
    gsl_vector_int_view row = gsl_matrix_int_row(m, i);
    XLAL_VBXMLO_gsl_vector_int(xml, "row", format, &row.vector);
  }
  XLAL_VBXMLO_EndTag(xml, name);

}

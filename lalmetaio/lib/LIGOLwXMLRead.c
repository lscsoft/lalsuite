/*
 * Copyright (C) 2007 Andres C. Rodriguez, Alexander Dietz, Duncan Brown,
 * Jolien Creighton, Kipp Cannon, Lisa M. Goggin, Patrick Brady, Robert
 * Adam Mercer, Saikat Ray-Majumder, Anand Sengupta, Stephen Fairhurst,
 * Xavier Siemens, Craig Robinson , Sean Seader, Thomas Cokelaer
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with with program; see the file COPYING. If not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301  USA
 */

/*-----------------------------------------------------------------------
 *
 * File Name: LIGOLwXMLRead.c
 *
 * Author: Brown, D. A.
 *
 *-----------------------------------------------------------------------
 */

/**
 * \author Cannon, K. C. and Brown, D. A. and Fairhurst, S.
 * \file
 * \ingroup lalmetaio_general
 *
 * \brief Routines to read tabular data from LIGO lightweight XML files.
 *
 * ### Description ###
 *
 * ### Algorithm ###
 *
 * None.
 *
 * ### Uses ###
 *
 * ### Notes ###
 *
 * %% Any relevant notes.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <metaio.h>

#include <lal/Date.h>
#include <lal/LALConstants.h>
#include <lal/LALStdio.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/XLALError.h>

/**
 * Test a LIGO Light Weight XML file for the presence of a specific table.
 * Returns > 0 if the document contains the table, 0 if the document does
 * not contain the table, and < 0 on error.
 *
 * BUGS:
 *
 * - This function can't tell the difference between a missing table and an
 * unparseable document.  This is a limitation in libmetaio.
 *
 * - This function parses the entire file to determine if the table is
 * present, which is slow.
 *
 * - This entire approach to XML I/O is the wrong way to go.  What's needed
 * is a "load document" function, and a "save document" function.  DO NOT
 * attempt to write such functions by using this function to test for
 * every possible table one-by-one and loading the ones that are found.
 * Put the time into writing a proper XML I/O layer!!
 */
int XLALLIGOLwHasTable(const char *filename, const char *table_name)
{
	struct MetaioParseEnvironment env;
	int has_table;

	/*
	 * open the file and find table
	 */

	if(MetaioOpenFile(&env, filename)) {
		XLALPrintError("%s(): error opening \"%s\": %s\n", __func__, filename, env.mierrmsg.data ? env.mierrmsg.data : "unknown reason");
		XLAL_ERROR(XLAL_EIO);
	}

	/*
	 * find table.  note:  parse errors are interpreted as "table is
	 * missing".  metaio provides no other mechanism for testing for
	 * the presence of a table.
	 */

	has_table = !MetaioOpenTableOnly(&env, table_name);
	MetaioClearErrno(&env);

	/*
	 * close
	 */

	if(MetaioClose(&env)) {
		XLALPrintError("%s(): error parsing document after %s table: %s\n", __func__, table_name, env.mierrmsg.data ? env.mierrmsg.data : "unknown reason");
		XLAL_ERROR(XLAL_EIO);
	}

	/*
	 * done
	 */

	return has_table;
}


/**
 * Convenience wrapper for MetaioFindColumn(), translating to XLAL-style
 * error reporting and printing useful error messages on failure.  Returns
 * the integer index of the column, or a negative integer if the column is
 * not found or has the wrong type.  If required is non-zero, then an XLAL
 * error is reported if the column is missing, but if required is zero then
 * no error is generated for missing columns.  When a column is found, it's
 * type is checked and an XLAL error is reported if it does not match the
 * requested type.  Passing METAIO_TYPE_UNKNOWN disables the column type
 * test.
 */
int XLALLIGOLwFindColumn(
	struct MetaioParseEnvironment *env,
	const char *name,
	unsigned int type,
	int required
)
{
	int pos = MetaioFindColumn(env, name);
	if(pos >= 0) {
		/* column was found, check type */
		if(type != METAIO_TYPE_UNKNOWN && env->ligo_lw.table.col[pos].data_type != type) {
			XLALPrintError("%s(): column \"%s\" has wrong type\n", __func__, name);
			XLAL_ERROR(XLAL_EDATA);
		}
	} else if(required) {
		/* required column is missing */
		XLALPrintError("%s(): missing required column \"%s\"\n", __func__, name);
		XLAL_ERROR(XLAL_EDATA);
	}
	return pos;
}


/**
 * Convenience function to extract the integer part of an ilwd:char ID
 * string with some error checking.  If either of ilwd_char_table_name or
 * ilwd_char_column_name is not NULL, then the corresponding portion of the
 * ilwd:char string must match it exactly.  The return value is the
 * recovered integer suffix or < 0 on failure.
 */


long long XLALLIGOLwParseIlwdChar(
	const struct MetaioParseEnvironment *env,
	int column_number,
	const char *ilwd_char_table_name,
	const char *ilwd_char_column_name
)
{
	char *fmt;
	const char *ilwd_char = env->ligo_lw.table.elt[column_number].data.lstring.data;
	long long id;

	/*
	 * 8 = 1 for the '\0', 2 for the ':' characters, and 5 for the
	 * "%%lld" string
	 */

	fmt = malloc(strlen(ilwd_char_table_name ? ilwd_char_table_name : "%*[^:]") + strlen(ilwd_char_column_name ? ilwd_char_column_name : "%*[^:]") + 8);
	if(!fmt)
		XLAL_ERROR(XLAL_ENOMEM);

	sprintf(fmt, "%s:%s:%%lld", ilwd_char_table_name ? ilwd_char_table_name : "%*[^:]", ilwd_char_column_name ? ilwd_char_column_name : "%*[^:]");

	if(sscanf(ilwd_char, fmt, &id) < 1) {
		free(fmt);
		XLALPrintError("%s(): invalid %s \"%s\" for %s\n", __func__, ilwd_char_column_name ? ilwd_char_column_name : "ID", ilwd_char, ilwd_char_table_name ? ilwd_char_table_name : "table");
		XLAL_ERROR(XLAL_EDATA);
	}

	free(fmt);

	return id;
}

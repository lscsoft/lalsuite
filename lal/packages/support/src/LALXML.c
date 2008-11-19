/*
*  Copyright (C) 2007 Jolien Creighton
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
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <lal/LALXML.h>
#include <lal/XLALError.h>

static void print_element_names(xmlNode *node)
{
	xmlNode *cur;
	for (cur = node; cur; cur = cur->next) {
		if (cur->type == XML_ELEMENT_NODE)
			printf("node type: Element, name: %s\n", cur->name);
		print_element_names(cur->children);
	}
	return;
}

int XLALXMLFilePrintElements(const char *fname)
{
	static const char file[] = "XLALXMLFilePrintElements";
	xmlDoc  *doc;

	/* make sure that the shared library is the same as the
	 * library version the code was compiled against */
	LIBXML_TEST_VERSION

	if (!(doc = xmlReadFile(fname, NULL, 0)))
		XLAL_ERROR(file, XLAL_EIO);
	print_element_names(xmlDocGetRootElement(doc));
	xmlFreeDoc(doc);
	xmlCleanupParser(); /* free global variables in parser */
	return 0;
}

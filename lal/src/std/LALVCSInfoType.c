/*
 * Copyright (C) 2014, 2016 Karl Wette
 * Copyright (C) 2009-2013 Adam Mercer
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA 02111-1307 USA
 */

#include <string.h>

#include <lal/LALStdlib.h>
#include <lal/LALString.h>
#include <lal/LALVCSInfoType.h>

char *XLALVCSInfoString(const LALVCSInfoList vcs_list, const int verbose, const char *prefix)
{

  /* check input */
  XLAL_CHECK_NULL( vcs_list != NULL, XLAL_EFAULT );

  /* generate VCS and build information string */
  char *str = NULL;
  for ( size_t i = 0; vcs_list[i] != NULL; ++i ) {
    if ( verbose ) {
      str = XLALStringAppendFmt( str,
                                 "%s-Version: %s\n"
                                 "%s-Id: %s\n"
                                 "%s-Date: %s\n"
                                 "%s-Branch: %s\n"
                                 "%s-Tag: %s\n"
                                 "%s-Status: %s\n"
                                 "%s-Configure-Args: %s\n"
                                 "%s-Configure-Date: %s\n"
                                 "%s-Build-Date: %s\n",
                                 vcs_list[i]->name, vcs_list[i]->version,
                                 vcs_list[i]->name, vcs_list[i]->vcsId,
                                 vcs_list[i]->name, vcs_list[i]->vcsDate,
                                 vcs_list[i]->name, vcs_list[i]->vcsBranch,
                                 vcs_list[i]->name, vcs_list[i]->vcsTag,
                                 vcs_list[i]->name, vcs_list[i]->vcsStatus,
                                 vcs_list[i]->name, vcs_list[i]->configureArgs,
                                 vcs_list[i]->name, vcs_list[i]->configureDate,
                                 vcs_list[i]->name, vcs_list[i]->buildDate
        );
    } else {
      str = XLALStringAppendFmt( str,
                                 "%s: %s (%s %s)\n",
                                 vcs_list[i]->name,
                                 vcs_list[i]->version,
                                 vcs_list[i]->vcsClean,
                                 vcs_list[i]->vcsId
        );
    }
    XLAL_CHECK_NULL( str != NULL, XLAL_EFUNC );
  }

  /* add prefix */
  if ( prefix != NULL ) {
    char *new_str = NULL;
    char *ptr = str;
    char *line = XLALStringToken( &ptr, "\n", 0 );
    while ( line != NULL ) {
      new_str = XLALStringAppendFmt( new_str, "%s%s\n", prefix, line );
      XLAL_CHECK_NULL( new_str != NULL, XLAL_EFUNC );
      line = XLALStringToken( &ptr, "\n", 0 );
    }
    XLALFree( str );
    str = new_str;
  }

  return str;

}

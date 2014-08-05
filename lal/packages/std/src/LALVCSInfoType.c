/*
 * LALVCSInfoType.c - LAL VCS Information Type
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
 *
 * Copyright (C) 2009-2013 Adam Mercer
 */

#include <string.h>
#include <lal/LALVCSInfoType.h>

/* function to compare two LALVCSInfo structures */
int XLALVCSInfoCompare(const LALVCSInfo *header, const LALVCSInfo *library)
{
  /* check for header/library mismatch */
  if (strcmp(header->name, library->name) || \
      strcmp(header->version, library->version) || \
      strcmp(header->vcsId, library->vcsId) || \
      strcmp(header->vcsDate, library->vcsDate) || \
      strcmp(header->vcsBranch, library->vcsBranch) || \
      strcmp(header->vcsTag, library->vcsTag) || \
      strcmp(header->vcsAuthor, library->vcsAuthor) || \
      strcmp(header->vcsCommitter, library->vcsCommitter) || \
      strcmp(header->vcsStatus, library->vcsStatus))
  {
    /* version mismatch */
    return 1;
  }

  return 0;
}

/*
 * vim: tw=0 ts=2 et
 */

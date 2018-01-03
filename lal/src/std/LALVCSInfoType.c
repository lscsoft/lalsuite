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

int XLALVCSInfoCompare(const LALVCSInfo *vcs1, const LALVCSInfo *vcs2)
{
  if (strcmp(vcs1->name, vcs2->name) || \
      strcmp(vcs1->version, vcs2->version) || \
      strcmp(vcs1->vcsId, vcs2->vcsId) || \
      strcmp(vcs1->vcsDate, vcs2->vcsDate) || \
      strcmp(vcs1->vcsBranch, vcs2->vcsBranch) || \
      strcmp(vcs1->vcsTag, vcs2->vcsTag) || \
      strcmp(vcs1->vcsAuthor, vcs2->vcsAuthor) || \
      strcmp(vcs1->vcsCommitter, vcs2->vcsCommitter) || \
      strcmp(vcs1->vcsStatus, vcs2->vcsStatus))
  {
    /* vcs1 != vcs2 */
    return 1;
  }

  /* vcs1 == vcs2 */
  return 0;
}

/*
 * vim: tw=0 ts=2 et
 */

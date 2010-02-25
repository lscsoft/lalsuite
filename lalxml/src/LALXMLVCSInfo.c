/*
 * LALXMLVCSInfo.c - LALXML VCS Information
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
 * Copyright (C) 2009,2010 Adam Mercer
 */

#include <lal/LALVCSInfo.h>

#include <lal/LALXMLConfig.h>
#include <lal/LALXMLVCSInfo.h>

/* global variables for vcs information */
const char *lalXMLVCSVersion = LALXML_VERSION;
const char *lalXMLVCSId = LALXML_VCS_ID;
const char *lalXMLVCSDate = LALXML_VCS_DATE;
const char *lalXMLVCSBranch = LALXML_VCS_BRANCH;
const char *lalXMLVCSTag = LALXML_VCS_TAG;
const char *lalXMLVCSAuthor = LALXML_VCS_AUTHOR;
const char *lalXMLVCSCommitter = LALXML_VCS_COMMITTER;
const char *lalXMLVCSStatus = LALXML_VCS_STATUS;

/* global variables for vcs information - identable */
const char *lalXMLVCSIdentId = LALXML_VCS_IDENT_ID;
const char *lalXMLVCSIdentDate = LALXML_VCS_IDENT_DATE;
const char *lalXMLVCSIdentBranch = LALXML_VCS_IDENT_BRANCH;
const char *lalXMLVCSIdentTag = LALXML_VCS_IDENT_TAG;
const char *lalXMLVCSIdentAuthor = LALXML_VCS_IDENT_AUTHOR;
const char *lalXMLVCSIdentCommitter = LALXML_VCS_IDENT_COMMITTER;
const char *lalXMLVCSIdentStatus = LALXML_VCS_IDENT_STATUS;

/* vcs information structure */
struct tagLALVCSInfo lalXMLVCSInfo = { \
  LALXML_VERSION, \
  LALXML_VCS_ID, \
  LALXML_VCS_DATE, \
  LALXML_VCS_BRANCH, \
  LALXML_VCS_TAG, \
  LALXML_VCS_AUTHOR, \
  LALXML_VCS_COMMITTER, \
  LALXML_VCS_STATUS \
};

/*
 * vim: tw=0 ts=2 et
 */

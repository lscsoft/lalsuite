/*
 * LALSimulationVCSInfo.c - LALSimulation VCS Information
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

#include <lal/LALSimulationConfig.h>
#include <lal/LALSimulationVCSInfo.h>

/* global variables for vcs information */
const char *lalSimulationVCSVersion = LALSIMULATION_VERSION;
const char *lalSimulationVCSId = LALSIMULATION_VCS_ID;
const char *lalSimulationVCSDate = LALSIMULATION_VCS_DATE;
const char *lalSimulationVCSBranch = LALSIMULATION_VCS_BRANCH;
const char *lalSimulationVCSTag = LALSIMULATION_VCS_TAG;
const char *lalSimulationVCSAuthor = LALSIMULATION_VCS_AUTHOR;
const char *lalSimulationVCSCommitter = LALSIMULATION_VCS_COMMITTER;
const char *lalSimulationVCSStatus = LALSIMULATION_VCS_STATUS;

/* global variables for vcs information - identable */
const char *lalSimulationVCSIdentId = LALSIMULATION_VCS_IDENT_ID;
const char *lalSimulationVCSIdentDate = LALSIMULATION_VCS_IDENT_DATE;
const char *lalSimulationVCSIdentBranch = LALSIMULATION_VCS_IDENT_BRANCH;
const char *lalSimulationVCSIdentTag = LALSIMULATION_VCS_IDENT_TAG;
const char *lalSimulationVCSIdentAuthor = LALSIMULATION_VCS_IDENT_AUTHOR;
const char *lalSimulationVCSIdentCommitter = LALSIMULATION_VCS_IDENT_COMMITTER;
const char *lalSimulationVCSIdentStatus = LALSIMULATION_VCS_IDENT_STATUS;

/* vcs information structure */
struct tagLALVCSInfo lalSimulationVCSInfo = { \
  LALSIMULATION_VCS_NAME, \
  LALSIMULATION_VERSION, \
  LALSIMULATION_VCS_ID, \
  LALSIMULATION_VCS_DATE, \
  LALSIMULATION_VCS_BRANCH, \
  LALSIMULATION_VCS_TAG, \
  LALSIMULATION_VCS_AUTHOR, \
  LALSIMULATION_VCS_COMMITTER, \
  LALSIMULATION_VCS_STATUS \
};

/*
 * vim: tw=0 ts=2 et
 */

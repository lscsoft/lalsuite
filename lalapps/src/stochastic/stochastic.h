/*
 * stochastic.h - SGWB Standalone Analysis Pipeline
 *
 * Copyright (C) 2002-2006,2009 Adam Mercer
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 *
 * $Id$
 */

#ifndef STOCHASTIC_H
#define STOCHASTIC_H

#include <lal/lalGitID.h>
#include <lalappsGitID.h>

/* xml process param table helper */
#define ADD_PROCESS_PARAM(pptype, format, ppvalue) \
	  this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
    calloc(1, sizeof(ProcessParamsTable)); \
  snprintf(this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
			      PROGRAM_NAME); \
  snprintf(this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
			      long_options[option_index].name); \
  snprintf(this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype); \
  snprintf(this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue);

/* window duration for PSD estimation */
#define PSD_WINDOW_DURATION 4

#endif /* STOCHASTIC_H */

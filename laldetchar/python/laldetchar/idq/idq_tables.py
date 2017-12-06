# Copyright (C) 2014 Reed Essick
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

## addtogroup pkg_py_laldetchar_idq
## Synopsis
# ~~~
# from laldetchar.idq import idq_tables
# ~~~
# \author Reed Essick (<reed.essick@ligo.org>)

"""Module to keep definitions of various tables for iDQ pipeline."""

from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import types as ligolwtypes
from glue.ligolw import ilwd
from laldetchar import git_version

__author__ = 'Reed Essick <reed.essick@ligo.org>'
__version__ = git_version.id
__date__ = git_version.date


## addtogroup pkg_py_laldetchar_idq_idq_tables
# @{

# contains definitions for iDQ xml table classes and objects
#
# =============================================================================
#
#                              glitch:table
#
# =============================================================================
#

class GlitchTable(table.Table):

    tableName = 'glitch:table'
    validcolumns = {
        'event_id': 'ilwd:char',
        'ifo': 'lstring',
        'gps': 'int_4s',
        'gps_ns': 'int_4s',
        'rank': 'real_4',
        'fap': 'real_4',
        'likelihood': 'real_4',
        }


class Glitch(object):

    __slots__ = GlitchTable.validcolumns.keys()


GlitchTable.RowType = Glitch


#
# =============================================================================
#
#                              ovl_data:table
#
# =============================================================================
#

class OVLDataTable(table.Table):

    tableName = 'ovl_data:table'
    validcolumns = {
        'event_id': 'ilwd:char',
        'ifo': 'lstring',
        'aux_channel': 'lstring',
        'veto_thr': 'real_4',
        'veto_win': 'real_4',
        }


class OVLData(object):

    __slots__ = OVLDataTable.validcolumns.keys()


OVLDataTable.RowType = OVLData

##@}


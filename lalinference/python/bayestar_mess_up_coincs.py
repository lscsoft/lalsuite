#
# Copyright (C) 2013-2017  Leo Singer
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
"""Mess up search pipeline output by performing nonphysical transformations of
the trigger parameters and SNR time series."""


# Standard library imports
import argparse
import sys


# LALSuite improts
from glue.ligolw.ligolw import LIGO_LW
from glue.ligolw.param import get_pyvalue
from glue.ligolw.array import get_array
from glue.ligolw.utils import load_fileobj, write_fileobj
from glue.ligolw.table import get_table
from glue.ligolw.lsctables import SnglInspiralTable
from lalinference.bayestar import command
from lalinference.bayestar.decorator import as_dict
from lalinference.bayestar.ligolw import LSCTablesAndSeriesContentHandler
import lal


# Available detectors
available_ifos = sorted(det.frDetector.prefix
    for det in lal.CachedDetectors)


# Transformations
def swap(xmldoc, ifo1, ifo2):
    """Swap a pair of detectors"""
    mapping = {ifo1: ifo2, ifo2: ifo1}
    sngls = get_table(xmldoc, SnglInspiralTable.tableName)
    for row in sngls:
        old_ifo = row.ifo
        try:
            new_ifo = mapping[old_ifo]
        except KeyError:
            pass
        else:
            row.ifo = new_ifo


def is_COMPLEX8TimeSeries(elem):
    return elem.tagName == LIGO_LW.tagName and getattr(elem, 'Name', None) == 'COMPLEX8TimeSeries'


@as_dict
def get_snr_series(xmldoc):
    for elem in xmldoc.getElements(is_COMPLEX8TimeSeries):
        try:
            event_id = get_pyvalue(elem, 'event_id')
            array = get_array(elem, 'snr')
        except ValueError:
            pass
        else:
            yield event_id, array
    

# Transformations
def conj(xmldoc, ifos):
    """Negate the phases on arrival in one or more detectors"""
    sngls = get_table(xmldoc, SnglInspiralTable.tableName)
    snrs = get_snr_series(xmldoc)
    for row in sngls:
        if row.ifo in ifos:
            row.coa_phase *= -1
        try:
            series = snrs[row.event_id]
        except KeyError:
            pass
        else:
            series.array[2, :] *= -1


# Command line interface
parser = command.ArgumentParser()
subparsers = parser.add_subparsers()
def add_parser(func):
    subparser = subparsers.add_parser(func.__name__, help=func.__doc__)
    subparser.set_defaults(func=func.__name__)
    return subparser

subparser = add_parser(swap)
subparser.add_argument('ifo1', choices=available_ifos)
subparser.add_argument('ifo2', choices=available_ifos)

subparser = add_parser(conj)
subparser.add_argument('ifos', choices=available_ifos, nargs='+')

args = parser.parse_args()


# Input document, process it, and write it out.
xmldoc, _ = load_fileobj(
    sys.stdin, contenthandler=LSCTablesAndSeriesContentHandler)
command.register_to_xmldoc(xmldoc, parser, args)
kwargs = dict(args.__dict__)
func = locals()[kwargs.pop('func')]
func(xmldoc, **kwargs)
write_fileobj(xmldoc, sys.stdout)

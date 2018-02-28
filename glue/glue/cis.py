# Copyright (C) 2013 Duncan Macleod
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
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

"""Interface to the LIGO Channel Information System

Queries can be made using any part of a channel name, returning a list
of entries (key,value dicts) of found channels::

>>> from glue import cis
>>> out = cis.query('PSL-ISS_PSB_OUT_DQ')
>>> print([channel['name'] for channel in out])
[u'H1:PSL-ISS_PDB_OUT_DQ', u'H2:PSL-ISS_PDB_OUT_DQ', u'L1:PSL-ISS_PDB_OUT_DQ']

(example accessed August 26 2013).
"""

import json
import numpy
from urllib2 import HTTPError

from glue.auth.saml import HTTPNegotiateAuthHandler

from . import (git_version, auth)

__author__ = "Duncan Macleod <duncan.macleod@ligo.org>"
__version__ = git_version.id
__date__ = git_version.date

CIS_API_URL = 'https://cis.ligo.org/api/channel'
CIS_DATA_TYPE = {4: numpy.float32}


def query(name, debug=False):
    """Query the Channel Information System for details on the given
    channel name

    Example::

    >>> from glue import cis
    >>> out = cis.query('PSL-ISS_PSB_OUT_DQ')
    >>> print([channel['name'] for channel in out])
    [u'H1:PSL-ISS_PDB_OUT_DQ', u'H2:PSL-ISS_PDB_OUT_DQ',
     u'L1:PSL-ISS_PDB_OUT_DQ']

    Parameters
    ----------
    name : str
        Name of the channel of interest, or part of channel to query

    Returns
    -------
    list
        list of (key, value) dicts recovered from the CIS

    Raises
    ------
    ValueError
        if no channels are found matching the given name
    """
    url = '%s/?q=%s' % (CIS_API_URL, name)
    more = True
    out = []
    while more:
        reply = _get(url, debug=debug)
        if reply['count'] == 0:
            raise ValueError("No channels found with name '%s'" % name)
        try:
           out.extend(reply[u'results'])
        except KeyError:
           pass
        more = reply.has_key('next') and reply['next'] is not None
        if more:
            url = reply['next']
        else:
            break
    out.sort(key=lambda c: c['name'])
    return out


def _get(url, debug=False):
    """Perform a GET query against the CIS
    """
    try:
        response = auth.request_ligodotorg(url, debug=debug)
    except HTTPError:
        raise ValueError("Channel named '%s' not found in Channel "
                         "Information System. Please double check the "
                         "name and try again." % url.strip('=')[-1])
    return json.loads(response.read())

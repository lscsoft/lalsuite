# Copyright (C) 2013 Ruslan Vaulin
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

from optparse import *

from laldetchar import git_version

__author__ = 'Ruslan Vaulin <ruslan.vaulin@ligo.org>'
__version__ = git_version.id
__date__ = git_version.date

description = \
    """This program reads the list of channels from KleineWelle config file and writes them into  selected_channels.txt file. This file is used as an input for iDQ pipeline.
"""

parser = OptionParser(version='Name: %%prog\n%s'
                      % git_version.verbose_msg, usage='%prog [options]'
                      , description=description)
parser.add_option('', '--input-file', type='string',
                  help='KleineWelle configuration file')

(opts, args) = parser.parse_args()

kw_file = open(opts.input_file, 'r')
file = open('selected_channels.txt', 'w')

for (i, line) in enumerate(kw_file):
    if i >= 5:
        words = line.split()
        channel_line = words[1].replace(':', '_') + '_' + words[2] \
            + '_' + words[3] + '\n'
        file.write(channel_line)

file.close()


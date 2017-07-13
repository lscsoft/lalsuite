#
# Copyright (C) 2017  Leo Singer
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Create a P-P plot to compare a posterior sample chain with a sky map.
"""


# Command line interface.
from argparse import FileType
from lalinference.bayestar import command
parser = command.ArgumentParser()
parser = command.ArgumentParser(parents=[command.figure_parser])
parser.add_argument(
    'skymap', metavar='SKYMAP.fits[.gz]', type=FileType('rb'),
    help='FITS sky map file')
parser.add_argument(
    'samples', metavar='SAMPLES.hdf5', type=FileType('rb'),
    help='HDF5 posterior sample chain file')
args = parser.parse_args()


# Late imports.
from matplotlib import pyplot as plt
from lalinference import io
import lalinference.plot
from lalinference.bayestar.postprocess import find_injection_moc


# Read input.
skymap = io.read_sky_map(args.skymap.name, moc=True)
chain = io.read_samples(args.samples.name)

# Calculate P-P plot.
result = find_injection_moc(skymap, chain['ra'], chain['dec'], chain['dist'])


# Make Matplotlib figure.
fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111, projection='pp_plot')
ax.add_diagonal()
ax.add_series(result.searched_prob, label='R.A., Dec.')
ax.add_series(result.searched_prob_dist, label='Distance')
ax.add_series(result.searched_prob_vol, label='Volume')
ax.set_xlabel('searched probability')
ax.set_ylabel('cumulative fraction of posterior samples')
ax.set_title(args.skymap.name)
ax.legend()
ax.grid()


# Show or save output.
args.output()

# Copyright (C) 2021 Rodrigo Tenorio
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

## \defgroup lalpulsar_py_lineFileParser LineFileParser
## \ingroup lalpulsar_python
"""
Parse identified and unidentified line files following the O3a convention.

The identified lines file contains detector spectral artifacts (lines and combs)
verified to be non-astrophysical in origin.
Each line of the file identifies an artifact with certain fundamental frequency,
(possibly non-identical) frequency wings,
the number of visible harmonics (if identified as a comb),
the scaling of the affected band with harmonics (constant or linearly scaling),
and a systematic off-set with respect to the specified frequency.

The unidentified lines file contains spectral artifacts that cannot yet be
convincingly associated with non-astrophysical detector noise.
These files only list the frequency of the most prominent peak.

This module provides a class, `LineFileParser`, to seamlessly read those files and return
the list of lines as a list of [left_most_line_frequency, right_most_line_frequency]
intervals.
This format is easier to use for example when applying line vetoes on CW outliers.

Input files are fed into the class using `LineFileParser.parse_identified_lines_csv` or
`LineFileParser.parse_unidentified_lines_csv`.
Intervals are accesible through attributes `LineFileParser.lines_left_side` and
`LineFileParser.lines_right_side`.
If several files are read, their resulting intervals are concatenated
(but not sorted, and duplicate or overlapping entries are kept as they are).
"""
## @{

import numpy as np

from . import git_version

__author__ = "Rodrigo Tenorio <rodrigo.tenorio@ligo.org>"
__version__ = git_version.id
__date__ = git_version.date


class LineFileParser:
    """
    # LineFileParser
    This class is taylored to work using the O3 - era line files.
    There are essentially two types of these files,
    their main difference being the number of columns.

    ## Identified lines file:

    - frequency: Column containing the central frequency of the line.
    - left_wing: Column containing the frequency limit towards
    lower frequencies of the line as a width (central frequency - low frequency limit) Hz.
    - right_wing: Column containing the frequency limit towards
    higher frequencies of the line as a width (high frequency limit - central frequency) Hz.
    - first_harmonic: Column with the index of the first visible harmonic of the line.
    - last_harmonic: Column with the index of the last visible harmonic of the line.
    - comb_type: Column containing line type
    (0=single line, 1=comb, or 2=comb whose width scales with the harmonic number).
    - offset: Column with the frequency line off-set.

    ## Unidentified lines file:

    - frequency: Column containing the central frequency of the line.
    """

    default_genfromtxt_kwargs = {"delimiter": ",", "skip_header": 4}
    identified_lines_keys = {
        "frequency": 0,
        "left_wing": 5,
        "right_wing": 6,
        "first_harmonic": 3,
        "last_harmonic": 4,
        "comb_type": 1,
        "offset": 2,
    }

    unidentified_lines_keys = {"frequency": 0}

    def __init__(self):
        """
        Read and expand line files into frequency intervals.
        If multiple files are read, resulting interval are concatenated.

        Intervals can be retrieved using the attributes `lines_left_side` and `lines_right_side`.

        """
        self.lines_left_side = None
        self.lines_right_side = None

        self._lines_are_set = False

    def _check_columns(self, columns, default_keys):
        columns = columns or default_keys
        if not all(key in columns for key in default_keys):
            raise ValueError(
                f"`columns` dictionary {columns} does not contain "
                f"all of the required keys: {list(default_keys.keys())}"
            )
        return columns

    def parse_identified_lines_csv(
        self, lines_file, columns=None, extra_wing_Hz=0.0, genfromtxt_kwargs=None
    ):
        """
        Parse a csv file containing lines of known origin (Advanced LIGO format).

        @param lines_file: Path to csv format lines file.
        @param columns: Dictionary with header fields as key
        and the corresponding (0-based) column index as value.
        If None, default ordering specified in the class attribute will be used.
        @param extra_wing_Hz: Extra wings to add at both sides of the resulting intervals, in Hz.
        @param genfromtxt_kwargs: kwargs to be passed to numpy.genfromtxt.
        Default is `delimiter=",", skip_header=4`.
        """

        columns = self._check_columns(columns, self.identified_lines_keys)

        expanded_lines_and_wings = self._get_identified_lines_center_left_right(
            lines_file, columns, genfromtxt_kwargs
        )

        lines_left_side, lines_right_side = self._add_frequency_wings(
            *expanded_lines_and_wings, extra_wing_Hz
        )

        self._set_lines(lines_left_side, lines_right_side)

    def parse_unidentified_lines_csv(
        self, lines_file, columns=None, extra_wing_Hz=0.0, genfromtxt_kwargs=None
    ):
        """
        Parse a csv file containing unidentified lines (Advanced LIGO format).

        @param lines_file: Path to csv format lines file.
        @param columns: Dictionary with header fields as key
        and the corresponding (0-based) column index as value.
        If None, default ordering specified in the class attribute will be used.
        @param extra_wing_Hz: Extra wings to add at both sides of the resulting intervals, in Hz.
        @param genfromtxt_kwargs: kwargs to be passed to numpy.genfromtxt.
        Default is `delimiter=",", skip_header=4`.
        """
        columns = self._check_columns(columns, self.unidentified_lines_keys)

        unidentified_lines = np.genfromtxt(
            lines_file,
            usecols=[columns[key] for key in self.unidentified_lines_keys],
            **(genfromtxt_kwargs or self.default_genfromtxt_kwargs),
        )

        lines_left_side, lines_right_side = self._add_frequency_wings(
            unidentified_lines, 0.0, 0.0, extra_wing_Hz
        )
        self._set_lines(lines_left_side, lines_right_side)

    def _get_identified_lines_center_left_right(
        self, lines_file, columns, genfromtxt_kwargs=None
    ):
        lines_with_wings = np.genfromtxt(
            lines_file,
            usecols=[columns[key] for key in self.identified_lines_keys],
            **(genfromtxt_kwargs or self.default_genfromtxt_kwargs),
        )
        return self._expand_harmonics(*lines_with_wings.T)

    def _set_lines(self, lines_left_side, lines_right_side):
        """
        Properly add left and right boundaries to the class attributes.
        That means concatenating instead of overwriting if a list of lines
        has been already read. This may happen e.g. when reading identified
        and unidentified line file.
        """
        if not self._lines_are_set:
            self.lines_left_side = lines_left_side
            self.lines_right_side = lines_right_side
            self._lines_are_set = True
        else:
            self.lines_left_side = np.concatenate(
                (self.lines_left_side, lines_left_side)
            )
            self.lines_right_side = np.concatenate(
                (self.lines_right_side, lines_right_side)
            )

    def _expand_harmonics(
        self,
        central_frequency,
        left_wing,
        right_wing,
        first_harmonic,
        last_harmonic,
        comb_type,
        offset,
    ):
        """
        Known line files contain only one of the harmonics of the lines. This method
        expands those harmonics, explicitely adding frequencies and wings to the line
        list in order to apply the veto.

        Comb type are used to properly re-scale harmonic wings. As of now (O3a),
        0 means line, 1 means non-scaling and 2 means scaling.

        Offsets shift the whole left-center-right structure as an overall adding term.
        """
        harmonics_per_line = (last_harmonic - first_harmonic + 1).astype(np.int32)
        total_number_of_lines = np.sum(harmonics_per_line, dtype=np.int32)

        expanded_central_frequency = np.zeros(total_number_of_lines)
        expanded_left_wing = np.zeros(total_number_of_lines)
        expanded_right_wing = np.zeros(total_number_of_lines)

        # If scaling comb (i.e. comb type 2), scale wings
        dont_scale_wings = comb_type != 2

        line_pointer = 0
        for line in range(len(central_frequency)):
            harmonic_index = np.arange(first_harmonic[line], last_harmonic[line] + 1)
            wing_scaling = (
                np.ones(harmonics_per_line[line])
                if dont_scale_wings[line]
                else harmonic_index
            )

            expanded_central_frequency[
                line_pointer : line_pointer + harmonics_per_line[line]
            ] = (harmonic_index * central_frequency[line] + offset[line])
            expanded_left_wing[
                line_pointer : line_pointer + harmonics_per_line[line]
            ] = (wing_scaling * left_wing[line])
            expanded_right_wing[
                line_pointer : line_pointer + harmonics_per_line[line]
            ] = (wing_scaling * right_wing[line])

            line_pointer += harmonics_per_line[line]

        return expanded_central_frequency, expanded_left_wing, expanded_right_wing

    def _add_frequency_wings(self, central_frequency, left_wing, right_wing, extra_Hz):
        """
        Given a line frequency and its wings, convert to a range of frequencies
        occupied by the line. Extra bins are added according to the specified
        input.
        """
        lines_left_side = central_frequency - left_wing - extra_Hz
        lines_right_side = central_frequency + right_wing + extra_Hz
        return lines_left_side, lines_right_side


## @}

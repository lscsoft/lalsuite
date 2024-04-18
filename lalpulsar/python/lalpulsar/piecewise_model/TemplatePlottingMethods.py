# Copyright (C) 2019--2023 Benjamin Grace
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

import random

import matplotlib.pyplot as plt
import numpy as np

import lalpulsar as lp

from . import BasisFunctions as bf
from . import SemicoherentMetricMethods as scmm


def metric_ellipse(g11, mutilde, tol=0.01):
    zz = [0]
    x = np.cos(zz[-1])
    y = np.sin(zz[-1])
    rr = [
        np.sqrt(
            mutilde / (g11[0, 0] * x**2 + 2 * g11[0, 1] * x * y + g11[1, 1] * y**2)
        )
    ]

    dz = 2 * np.pi / 1000

    while zz[-1] < 2 * np.pi:
        z = zz[-1] + dz

        x = np.cos(z)
        y = np.sin(z)
        r = np.sqrt(
            mutilde / (g11[0, 0] * x**2 + 2 * g11[0, 1] * x * y + g11[1, 1] * y**2)
        )
        e = abs(r - rr[-1]) / rr[-1]

        if e < 0.1 * tol:
            dz *= 1.5
        elif e > tol:
            dz /= 3
        else:
            zz.append(z)
            rr.append(r)

    zz[-1] = zz[0]
    rr[-1] = rr[0]
    xx = rr * np.cos(zz)
    yy = rr * np.sin(zz)

    return [xx, yy]


# Plots all the ellipses in the given parameter space for 2 and 3 dimensions
def plottemplates(temps, x_cross_section, dim1, dim2, metric, max_mismatch, show=True):

    dim_truth_array = np.array([False] * len(temps[0]))
    dim_truth_array[dim1] = True
    dim_truth_array[dim2] = True

    cross_section_dims = np.logical_not(dim_truth_array)

    g00 = metric[cross_section_dims, :][:, cross_section_dims]
    g10 = metric[dim_truth_array, :][:, cross_section_dims]
    g11 = metric[dim_truth_array, :][:, dim_truth_array]

    g11invg10 = np.linalg.inv(g11) @ g10

    # x_cross_section = temps[0]#int(np.floor(len(temps) / 2))]

    x0 = x_cross_section[cross_section_dims]

    validellipses = 0
    counter = 0

    fig = plt.figure(figsize=(10, 10))
    ax = plt.axes()

    for i, temp in enumerate(temps):
        if i > 100000000000000:
            break

        if i % 10000 == 0:
            print("Template plot counter is at: " + str(i))

        x0bar = temp[cross_section_dims]

        x1bar = temp[dim_truth_array]
        dx0 = x0 - x0bar

        x1tilde = -g11invg10 @ dx0

        mutilde = max_mismatch - dx0 @ g00 @ dx0 + x1tilde @ g11 @ x1tilde

        if mutilde < 0:
            continue

        validellipses += 1
        centre = x1bar + x1tilde

        ellipse = metric_ellipse(g11, mutilde)

        ax.plot(ellipse[0] + centre[0], ellipse[1] + centre[1], "k")
        plt.scatter(centre[0], centre[1], s=15, color="orange")

    print("Valid ellipses: " + str(validellipses))
    if show:
        plt.show()


# Plots all the ellipses in the given parameter space for 2 and 3 dimensions
def plottemplatesspecifyaxis(
    temps, x_cross_section, dim1, dim2, metric, max_mismatch, axis, show=False
):

    print("Max Mismatch")
    print(max_mismatch)
    dim_truth_array = np.array([False] * len(temps[0]))
    dim_truth_array[dim1] = True
    dim_truth_array[dim2] = True

    cross_section_dims = np.logical_not(dim_truth_array)

    g00 = metric[cross_section_dims, :][:, cross_section_dims]
    g10 = metric[dim_truth_array, :][:, cross_section_dims]
    g11 = metric[dim_truth_array, :][:, dim_truth_array]

    g11invg10 = np.linalg.inv(g11) @ g10

    # x_cross_section = temps[0]#int(np.floor(len(temps) / 2))]

    x0 = x_cross_section[cross_section_dims]

    validellipses = 0
    counter = 0

    for i, temp in enumerate(temps):
        if i > 100000000000000:
            break

        if i % 10000 == 0:
            print("Template plot counter is at: " + str(i))

        x0bar = temp[cross_section_dims]

        x1bar = temp[dim_truth_array]
        dx0 = x0 - x0bar

        x1tilde = -g11invg10 @ dx0

        mutilde = max_mismatch - dx0 @ g00 @ dx0 + x1tilde @ g11 @ x1tilde

        if mutilde < 0:
            continue

        validellipses += 1
        centre = x1bar + x1tilde

        ellipse = metric_ellipse(g11, mutilde)

        axis.plot(ellipse[0] + centre[0], ellipse[1] + centre[1], "k")
        axis.scatter(centre[0], centre[1], s=15, color="orange")

    print("Valid ellipses: " + str(validellipses))
    if show:
        plt.show()


# Plots all the ellipses in the given parameter space for 2 and 3 dimensions
def plottemplatesalldims(
    temps, x_cross_section, segment, metric, max_mismatch, show=True
):

    ax00_01 = plt.subplot(5, 5, 1)
    ax00_02 = plt.subplot(5, 5, 6)  # , sharex=ax00_01)
    ax01_02 = plt.subplot(5, 5, 7)  # , sharey=ax00_02)
    ax00_10 = plt.subplot(5, 5, 11)  # , sharex=ax00_01)
    ax01_10 = plt.subplot(5, 5, 12)  # , sharex=ax01_02, sharey=ax00_01)
    ax02_10 = plt.subplot(5, 5, 13)  # , sharey=ax00_10)
    ax00_11 = plt.subplot(5, 5, 16)  # , sharex=ax00_01)
    ax01_11 = plt.subplot(5, 5, 17)  # , sharex=ax01_02, sharey=ax00_11)
    ax02_11 = plt.subplot(5, 5, 18)  # , sharex=ax02_10, sharey=ax00_11)
    ax10_11 = plt.subplot(5, 5, 19)  # , sharey=ax00_11)
    ax00_12 = plt.subplot(5, 5, 21)  # , sharex=ax00_01)
    ax01_12 = plt.subplot(5, 5, 22)  # , sharex=ax01_02, sharey=ax00_12)
    ax02_12 = plt.subplot(5, 5, 23)  # , sharex=ax02_10, sharey=ax00_12)
    ax10_12 = plt.subplot(5, 5, 24)  # , sharex=ax10_11, sharey=ax00_12)
    ax11_12 = plt.subplot(5, 5, 25)  # , sharey=ax00_12)

    axislist = [
        ax00_01,
        ax00_02,
        ax01_02,
        ax00_10,
        ax01_10,
        ax02_10,
        ax00_11,
        ax01_11,
        ax02_11,
        ax10_11,
        ax00_12,
        ax01_12,
        ax02_12,
        ax10_12,
        ax11_12,
    ]
    dimslist = [
        [0, 1],
        [0, 2],
        [1, 2],
        [0, 3],
        [1, 3],
        [2, 3],
        [0, 4],
        [1, 4],
        [2, 4],
        [3, 4],
        [0, 5],
        [1, 5],
        [2, 5],
        [3, 5],
        [4, 5],
    ]

    for i, axis in enumerate(axislist):

        # This case corresponds with the two frequency parameters (one on knot i, the other on knot i + 1). In this case, both parameters have similar
        # ranges and leads to a large number of ellipses having cross sections through the given plane. So we skip this case to save time, as it is
        # mostly uninteresting and computationally expensive
        if dimslist[i][0] == 0 and dimslist[i][1] == 3:
            continue

        plottemplatesspecifyaxis(
            temps,
            x_cross_section,
            segment * 3 + dimslist[i][0],
            segment * 3 + dimslist[i][1],
            metric,
            max_mismatch,
            axis,
            show=False,
        )
        print("Done dims: " + str(dimslist[i]))

    ax00_12.set_xlabel("f_" + str(segment) + "_0")
    ax01_12.set_xlabel("f_" + str(segment) + "_1")
    ax02_12.set_xlabel("f_" + str(segment) + "_2")
    ax10_12.set_xlabel("f_" + str(segment + 1) + "_0")
    ax11_12.set_xlabel("f_" + str(segment + 1) + "_1")

    ax00_01.set_ylabel("f_" + str(segment) + "_1")
    ax00_02.set_ylabel("f_" + str(segment) + "_2")
    ax00_10.set_ylabel("f_" + str(segment + 1) + "_0")
    ax00_11.set_ylabel("f_" + str(segment + 1) + "_1")
    ax00_12.set_ylabel("f_" + str(segment + 1) + "_2")

    print("All dim's done")
    if show:
        plt.show()


# Shows a 2D projection of all templates randomly chosen to be within the parameter space bounds for axis which are all combinations of parameters on a given segment.
# Each plot also marks the median value for each parameter (orange dot) with a cross hair representing the metric spacing at that point. The metric spacing is usually
# quite small, and zooming in may be required to see it.
def parameter_space_visualiser(tbank, points, segment, show=True, log=True):

    templates = []

    dimensions = len(bf.knotslist) * tbank.s

    kmin = tbank.kmin
    kmax = tbank.kmax

    # Randomly generating templates within bounds
    for i in range(points):

        this_template = []

        for dim in range(dimensions):
            knot_num = int(np.floor(dim / 3))

            if knot_num > 0:
                seglength = bf.knotslist[knot_num] - bf.knotslist[knot_num - 1]
            else:
                seglength = 0

            # FIXME
            lowerbound = lp.PiecewiseParameterBounds(
                dim,
                this_template,
                -1,
                tbank.fmin,
                tbank.fmax,
                tbank.nmin,
                tbank.nmax,
                kmin,
                kmax,
                seglength,
            )
            upperbound = lp.PiecewiseParameterBounds(
                dim,
                this_template,
                1,
                tbank.fmin,
                tbank.fmax,
                tbank.nmin,
                tbank.nmax,
                kmin,
                kmax,
                seglength,
            )

            parameter = random.uniform(lowerbound, upperbound)

            this_template.append(parameter)

        templates.append(this_template)

    # Finding the median value of all randomly generated templates
    median_list = []

    for parameter in range(dimensions):

        parameter_list = []

        for template in templates:
            parameter_list.append(template[parameter])

        median_list.append(np.median(parameter_list))

    # Calculating the metric, metric spacing and cross hair start and end points.
    metric = scmm.metric(tbank.s)
    metric_diagonals = np.diagonal(metric)

    metric_spacing_list = [
        np.sqrt(tbank.maxmismatch / elem) for elem in metric_diagonals
    ]

    metric_spacing_line_points = []

    for i, elem in enumerate(median_list):
        line_start = elem - 0.5 * metric_spacing_list[i]
        line_end = elem + 0.5 * metric_spacing_list[i]

        if i % 3 != 0 and log:
            line_start = np.log10(np.abs(line_start))
            line_end = np.log10(np.abs(line_end))

        metric_spacing_line_points.append([line_start, line_end])

    # Defining all subplots
    ax00_01 = plt.subplot(5, 5, 1)
    ax00_02 = plt.subplot(5, 5, 6)  # , sharex=ax00_01)
    ax01_02 = plt.subplot(5, 5, 7)  # , sharey=ax00_02)
    ax00_10 = plt.subplot(5, 5, 11)  # , sharex=ax00_01)
    ax01_10 = plt.subplot(5, 5, 12)  # , sharex=ax01_02, sharey=ax00_01)
    ax02_10 = plt.subplot(5, 5, 13)  # , sharey=ax00_10)
    ax00_11 = plt.subplot(5, 5, 16)  # , sharex=ax00_01)
    ax01_11 = plt.subplot(5, 5, 17)  # , sharex=ax01_02, sharey=ax00_11)
    ax02_11 = plt.subplot(5, 5, 18)  # , sharex=ax02_10, sharey=ax00_11)
    ax10_11 = plt.subplot(5, 5, 19)  # , sharey=ax00_11)
    ax00_12 = plt.subplot(5, 5, 21)  # , sharex=ax00_01)
    ax01_12 = plt.subplot(5, 5, 22)  # , sharex=ax01_02, sharey=ax00_12)
    ax02_12 = plt.subplot(5, 5, 23)  # , sharex=ax02_10, sharey=ax00_12)
    ax10_12 = plt.subplot(5, 5, 24)  # , sharex=ax10_11, sharey=ax00_12)
    ax11_12 = plt.subplot(5, 5, 25)  # , sharey=ax00_12)

    axislist = [
        ax00_01,
        ax00_02,
        ax01_02,
        ax00_10,
        ax01_10,
        ax02_10,
        ax00_11,
        ax01_11,
        ax02_11,
        ax10_11,
        ax00_12,
        ax01_12,
        ax02_12,
        ax10_12,
        ax11_12,
    ]
    dimslist = [
        [0, 1],
        [0, 2],
        [1, 2],
        [0, 3],
        [1, 3],
        [2, 3],
        [0, 4],
        [1, 4],
        [2, 4],
        [3, 4],
        [0, 5],
        [1, 5],
        [2, 5],
        [3, 5],
        [4, 5],
    ]

    for i, axis in enumerate(axislist):

        x_axis = []
        y_axis = []

        current_dim = 0

        x_index = tbank.s * segment + dimslist[i][0]
        y_index = tbank.s * segment + dimslist[i][1]

        # Getting all template parameters for the appropriate dimensions
        for temp in templates:

            x_val = temp[x_index]
            y_val = temp[y_index]

            if x_index % 3 and log:
                x_val = np.log10(np.abs(x_val))

            if y_index % 3 and log:
                y_val = np.log10(np.abs(y_val))

            x_axis.append(x_val)
            y_axis.append(y_val)

            current_dim += 1

        axis.scatter(x_axis, y_axis)

        # Plotting the median value and metric spacing cross hair
        median_x = median_list[x_index]
        median_y = median_list[y_index]

        if x_index % 3 != 0 and log:
            median_x = np.log10(np.abs(median_x))
        if y_index % 3 != 0 and log:
            median_y = np.log10(np.abs(median_y))

        line_start_x = metric_spacing_line_points[x_index][0]
        line_end_x = metric_spacing_line_points[x_index][1]

        line_start_y = metric_spacing_line_points[y_index][0]
        line_end_y = metric_spacing_line_points[y_index][1]

        axis.scatter([median_x], [median_y], color="orange")
        axis.plot([line_start_x, line_end_x], [median_y, median_y], color="black")
        axis.plot([median_x, median_x], [line_start_y, line_end_y], color="black")

        print("Done dims: " + str(dimslist[i]))

    ax00_12.set_xlabel("f_" + str(segment) + "_0")
    ax01_12.set_xlabel("f_" + str(segment) + "_1")
    ax02_12.set_xlabel("f_" + str(segment) + "_2")
    ax10_12.set_xlabel("f_" + str(segment + 1) + "_0")
    ax11_12.set_xlabel("f_" + str(segment + 1) + "_1")

    ax00_01.set_ylabel("f_" + str(segment) + "_1")
    ax00_02.set_ylabel("f_" + str(segment) + "_2")
    ax00_10.set_ylabel("f_" + str(segment + 1) + "_0")
    ax00_11.set_ylabel("f_" + str(segment + 1) + "_1")
    ax00_12.set_ylabel("f_" + str(segment + 1) + "_2")

    print("All dim's done")
    if show:
        plt.show()

# Copyright (C) 2013 Reed Essick, Ruslan Vaulin
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your# option
# ) any later version.
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
# from laldetchar.idq import idq_summary_plots
# ~~~
# \author Reed Essick (<reed.essick@ligo.org>), Ruslan Vaulin (<ruslan.vaulin@ligo.org>)

"""Module for plotting routines for idq pipeline."""

import sys
from laldetchar.idq import pdf_estimation as pdf_e
import numpy as np
import matplotlib
matplotlib.use('Agg')

# import matplotlib.pyplot as plt

import pylab as plt
import glob
import matplotlib.cm
from matplotlib import mpl

import scipy.interpolate as scipy_interp
from laldetchar import git_version

__author__ = \
    'Reed Essick (<reed.essick@ligo.org>), Ruslan Vaulin <ruslan.vaulin@ligo.org>'
__version__ = git_version.id
__date__ = git_version.date

## addtogroup pkg_py_laldetchar_idq_idq_summary_plots
# @{

# Plotting routines for idq pipeline

###########################################################################
###########################################################################

chanlistD = {
    'channel': 0,
    'n_cln': 1,
    'n_gch': 2,
    'c_fap': 3,
    'c_eff': 4,
    'eff/fap': 5,
    'vsig': 6,
    }

#########################
# set standard figure parameters

fig_type = '.png'

# ROC

__roc_yscale = 'linear'
__roc_xscale = 'log'

__roc_ymin = 1e-6
__roc_ymax = 1.0
__roc_xmin = 1e-6
__roc_xmax = 1.0

__roc_xlabel = 'False Alarm Probability'
__roc_ylabel = 'Efficiency'

__roc_legend_loc = 'best'
__roc_grid = True
__roc_line_width = 2

# kde

__kde_num_samples = 100

# __kde_num_samples = 10
# __kde_scale = 0.05

__kde_s = 0.01
__kde_s = 0.1

__kde_c_color = 'green'
__kde_c_label = 'not-glitch'
__kde_g_color = 'blue'
__kde_g_label = 'glitch'

__kde_yscale = 'linear'
__kde_ymin = 0.0
__kde_xscale = 'linear'
__kde_xmin = 0.0

__kde_xlabel = 'rank'
__kde_ylabel = 'Probability Density'

__kde_legend_loc = 'best'
__kde_grid = True

# kde Likelihood

__kde_L_color = 'blue'
__kde_L_label = 'Likelihood ratio'

__kde_L_yscale = 'linear'
__kde_L_ymin = 0.0
__kde_L_xscale = 'linear'
__kde_L_xmin = 0.0

__kde_L_xlabel = 'rank'
__kde_L_ylabel = 'Likelihood Ratio'

__kde_L_legend_loc = 'best'
__kde_L_grid = True

# Livetime trending plot

__livetime_trend_xlabel = 'GPS time, days'
__livetime_trend_ylabel = 'Fractional livetime, %'

__livetime_trend_legend_loc = 'best'
__livetime_trend_grid = True

# Glitch trending plot

__glitchrate_trend_xlabel = 'GPS time, days'
__glitchrate_trend_ylabel = 'Trigger rate, Hz'

__glitchrate_trend_legend_loc = 'best'
__glitchrate_trend_grid = True
__glitchrate_trend_yscale = 'log'

# Clean trending plot

__cleanrate_trend_xlabel = 'GPS time, days'
__cleanrate_trend_ylabel = 'Trigger rate, Hz'

__cleanrate_trend_legend_loc = 'best'
__cleanrate_trend_grid = True
__cleanrate_trend_yscale = 'log'

# Efficiency trending

__eff_trend_xlabel = 'GPS time, days'
__eff_trend_ylabel = 'Efficiency'

__eff_trend_legend_loc = 'best'
__eff_trend_grid = True


def classifier_colors(classifier):
    """
....Mapping between classifiers and colors to be used in the plots
...."""

    colors = {}

    colors['ovl'] = 'c'
    colors['mvsc'] = 'g'
    colors['ann'] = 'b'
    colors['svm'] = 'r'

    return colors[classifier]


def classifier_labels(classifier):
    """
....Mapping between classifiers and labels to be used in the plots
...."""

    return classifier


def get_roc_from_file(ROCfile):
    """
....Reads roc file and return roc curve: (fap,eff)
...."""

    (ranks, c_cln, c_gch, tot_cln, tot_gch) = ROC_to_rcg(ROCfile)
    if tot_cln == 0 and tot_gch == 0:  # no data to plot
        return ([], [], 0, 0)
    tot_cln = float(tot_cln)
    tot_gch = float(tot_gch)
    fap = []
    eff = []
    for (ind, c) in enumerate(c_cln):
        g = c_gch[ind]
        if g == tot_gch and c == tot_cln:  # don't include the (1,1) point on the ROC curves
            continue
        if tot_cln > 0:
            if c > 0:
                fap.append(c / tot_cln)
            else:

                  # get rid of 0.0 fap estimates

                fap.append(0.1 / tot_cln)
        else:
            fap.append(0)
        if tot_gch > 0:
            eff.append(g / tot_gch)
        else:
            eff.append(0)

    # convert fap and eff into arrays

    fap = np.asarray(fap)
    eff = np.asarray(eff)
    return (fap, eff, tot_cln, tot_gch)


def get_eff_at_fap_threshold(fap_array, eff_array, fap_thresh):
    """
....Returns efficiency for a given value of fap.
...."""

    # get unique faps

    unique_faps = np.unique(fap_array)

    # keep highest eff (fap,eff) pair out of those with the same fap

    unique_effs = []
    for fap in unique_faps:
        unique_effs.append(np.max(eff_array[np.nonzero(fap_array
                           == fap)[0]]))
    unique_effs = np.asarray(unique_effs)

    # check if roc curve consist of a single point

    if len(unique_faps) == 1:
        if unique_faps <= fap_thresh:
            return unique_effs[0]
        else:
            return 0.0

    # generate linear interpolation curve of eff vs fap

    eff_vs_fap_curve = scipy_interp.interp1d(unique_faps, unique_effs,
            bounds_error=False, fill_value=-1.0)

    # compute efficiency at fap_thresh

    eff = eff_vs_fap_curve(np.array([fap_thresh]))[0]

    # check if fap_thresh is not outside of unique_faps

    if eff == -1.0:

        # set eff to either min or max of unique_effs depending on whether fap_thresh is lesser or greater of all unique_faps

        if fap_thresh > np.max(unique_faps):
            eff = unique_effs[-1]
        elif fap_thresh < np.min(unique_faps):
            eff = 0.0
        else:
            print 'Inconsistency Error: returned value for eff at fap_thresh = ' \
                + str(fap_thresh) \
                + ' was None, but fap_thresh appears to be within the range of fap_array'

    return eff


#########################
#
#
#########################

def residuals(
    x1,
    y1,
    x2,
    y2,
    num_smpls=100,
    log_space=False,
    ):
    """
  returns the residuals (y1-y2) at linearly spaced points in the intersection of x1 and x2
  returns num_smpls sample points, and if log_space == True, these are logarithmically spaced (default is linear spacing)
  """

    min_x = max([min(x1), min(x2)])  # define the intersection of x1 and x2
    max_x = min([max(x1), max(x2)])
    if min_x >= max_x:  # intersection is the empty set
        return ([], [])

  # ## generate interpolation objects

    c1 = scipy_interp.interp1d(x1, y1)  # linear interpolation for curve 1
    c2 = scipy_interp.interp1d(x2, y2)
    if log_space:
        smples = np.logspace(np.log10(min_x), np.log10(max_x),
                             num_smpls)
    else:
        smpls = np.linspace(min_x, max_x, num_smpls)
    residulas = []
    for smpl in smpls:
        residuals.append(c1(smpl) - c2(smple))
    return (smpls, np.array(residuals))


#########################

def ROC_to_rcg(ROCfile):
    """ 
  reads in a standard ROC file and returns lists of ranks, cumulative_cleans, cumulative_glitches 
  """

    ranks = []
    c_cln = []
    c_gch = []
    tot_cln = False
    tot_gch = False
    i = 0  # counter to find first lines not commented
    f = open(ROCfile, 'r')
    for line in f:
        if line[0] != '#':
            if i == 0:
                tot_gch = float(line)
                i += 1
            elif i == 1:
                tot_cln = float(line)
                i += 1
            else:
                fields = line.strip().split()
                ranks.append(float(fields[0]))
                c_cln.append(float(fields[1]))
                c_gch.append(float(fields[2]))

    return (ranks, c_cln, c_gch, tot_cln, tot_gch)


#########################

def ROC_to_uniformROC(ROCfile, num_samples=100):
    """
  reads in a standard ROC file and returns ranks, c_cln, c_gch, tot_cln, tot_gch with uniform sampling in rank space
  """

    print ROCfile
    (ranks, c_cln, c_gch, tot_cln, tot_gch) = ROC_to_rcg(ROCfile)
    N = len(ranks)

    uniform_ranks = list(np.linspace(1, 0, num_samples + 1))
    uniform_cln = []  # holder for c_cln corresponding to uniform_ranks
    uniform_gch = []

    ind = 0
    cc = 0  # current value of c_cln corresponding to 'ind'
    cg = 0
    for ur in uniform_ranks:
        while ind < N and ranks[ind] >= ur:
            cc = c_cln[ind]
            cg = c_gch[ind]
            ind += 1
        uniform_cln.append(cc)
        uniform_gch.append(cg)

    return (uniform_ranks, uniform_cln, uniform_gch, tot_cln, tot_gch)


#########################

def rcg_to_ROC(
    filename,
    ranks,
    c_cln,
    c_gch,
    tot_cln,
    tot_gch,
    ):
    """
  writes a standard ROC file into filename
  """

    f = open(filename, 'w')
    print >> f, tot_gch
    print >> f, tot_cln
    for (r, c, g) in zip(ranks, c_cln, c_gch):
        print >> f, r, c, g
    f.close()

    return filename


#########################

def ROC_to_pwg_kde(
    ROCfile,
    num_samples=100,
    scale=0.1,
    s=0.1,
    pilot_interp='NONE',
    pilot_x='NONE',
    pilot_y='NONE',
    filename=False,
    ):
    """ 
  generates a point-wise-gaussian pdf estimate from a ROC file. 
    essentially delegates to pdf_estimation.point_wise_gaussian_kde(), and includes all the same options
  returns estimates of pdf's: p(r|g), p(r|c)

  if filename: generates a pdf estimate summary file similar in structur to the standard ROC file
  """

    (ranks, c_cln, c_gch, tot_cln, tot_gch) = ROC_to_rcg(ROCfile)
    c_cln = [int(i) for i in c_cln]
    c_gch = [int(i) for i in c_gch]

  # find out locations of all observations in rank space

    c_observ = []
    g_observ = []
    n_cln = 0  # holders to determine the number of glitches and cleans at a given rank (from cumulative numbers)
    n_gch = 0
    for ind in range(len(ranks)):
        rank = ranks[ind]
        c_observ += [rank for i in range(c_cln[ind] - n_cln)]  # number of cleans associated with this rank
        g_observ += [rank for i in range(c_gch[ind] - n_gch)]
        n_cln = c_cln[ind]  # increment
        n_gch = c_gch[ind]

#  if n_cln != 0 and n_gch != 0:
    # generate point-wise-gaussian kde
#    max_ranks = max(ranks)
#    min_ranks = min(ranks)

    eval = np.linspace(0.0, 1.0, num_samples)

    if n_cln != 0:
        c_kde = np.array(pdf_e.point_wise_gaussian_kde(
            eval,
            c_observ,
            scale=scale,
            s=s,
            pilot_interp=pilot_interp,
            pilot_x=pilot_x,
            pilot_y=pilot_y,
            ))
    else:
        c_kde = np.ones((len(eval), ))
    if n_gch != 0:
        g_kde = np.array(pdf_e.point_wise_gaussian_kde(
            eval,
            g_observ,
            scale=scale,
            s=s,
            pilot_interp=pilot_interp,
            pilot_x=pilot_x,
            pilot_y=pilot_y,
            ))
    else:
        g_kde = np.ones((len(eval), ))

    if filename:  # write summary file
        f = open(filename, 'w')
        print >> f, \
            '''# line 1 : maximum rank
# line 2 : minimum rank
# line 3-: rank, clean_kde, glitch_kde'''
        print >> f, max_ranks
        print >> f, min_ranks
        for ind in range(eval):
            print >> f, eval[ind], c_kde[ind], g_kde[ind]
        f.close()

    return (eval, c_kde, g_kde, tot_cln, tot_gch)


#########################

def ROC_to_pwg_kde_likelihood(
    ROCfile,
    num_samples=100,
    scale=0.1,
    s=0.1,
    pilot_interp='NONE',
    pilot_x='NONE',
    pilot_y='NONE',
    filename=False,
    ):
    """ 
  generates an estimate of the likelihood ratio using point-wise-gaussian KDEs for 1D probility distributions.
  
  if filename: generates a summary file
  """

    (eval, c_kde, g_kde, tot_cln, tot_gch) = ROC_to_pwd_kde(
        ROCfile,
        num_samples=num_samples,
        scale=scale,
        s=0.1,
        pilot_interp=pilot_interp,
        pilot_x=pilot_x,
        pilot_y=pilot_y,
        filename=False,
        )

    L = []
    for (i, c) in enumerate(c_kde):
        if c == 0:
            L.append('inf')
        else:
            L.append(g_kde[i] / c)

    if filename:
        f = open(filename, 'w')
        print >> f, '# rank Likelihood_ratio'
        for (i, e) in enumerate(eval):
            print >> f, e, L[i]
        f.close()

    return (eval, L)


#########################

def ROC_to_ROC_plot(
    ROCfiles,
    labels=False,
    colors=False,
    figure_name=False,
    write=True,
    annotate=True,
    ):
    """ 
  generate an ROC plot from ROCfiles. adopts standard axis conventions defined at beginning of module. 
  expects ROCfiles to be either a list or string of filenames for standard ROC files.
  """

  # set input parameters

    if labels:
        if len(labels) != len(ROCfiles):
            print 'labels and ROCfiles must have the same length for labels to be plotted.'
            labels = False
    if colors:
        if len(colors) != len(ROCfiles):
            print 'colors and ROCfiles must have the same length for colors to be applied'
            colors = False

  # generate plot

    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)
    nums = []  # keeps track of all data

  # iterate through data

    for (i, ROCfile) in enumerate(ROCfiles):

    # load roc quantities from file

        (fap, eff, tot_cln, tot_gch) = get_roc_from_file(ROCfile)
        if not len(fap):

      # no data to plot

            continue
        nums.append([i, tot_cln, tot_gch])

    # plot

        if colors and labels:
            ax.step(
                fap,
                eff,
                where='post',
                color=colors[i],
                label=labels[i],
                linewidth=__roc_line_width,
                )
        elif colors:

#    ....ax.plot(fap, eff, color=colors[i], label=labels[i])

            ax.step(fap, eff, where='post', color=colors[i],
                    linewidth=__roc_line_width)
        elif labels:

#    ....ax.plot(fap, eff, color=colors[i])

            ax.step(fap, eff, where='post', label=labels[i],
                    linewidth=__roc_line_width)
        else:

#    ....ax.plot(fap, eff, label=labels[i])

            ax.step(fap, eff, where='post', linewidth=__roc_line_width)

#    ....ax.plot(fap, eff)

    ax.plot(np.logspace(np.log10(__roc_xmin), np.log10(__roc_xmax),
            1000), np.logspace(np.log10(__roc_xmin),
            np.log10(__roc_xmax), 1000), 'k--')

  # set plot to standard form

    if labels:
        ax.legend(loc=__roc_legend_loc)
    ax.set_yscale(__roc_yscale)
    ax.set_xscale(__roc_xscale)
    ax.set_ylim(ymin=__roc_ymin, ymax=__roc_ymax)
    ax.set_xlim(xmin=__roc_xmin, xmax=__roc_xmax)

    ax.grid(__roc_grid)

    ax.set_xlabel(__roc_xlabel)
    ax.set_ylabel(__roc_ylabel)

    if annotate and len(nums) > 0:
        s1 = ''
        s2 = ''
        good = True  # figure out if all ROC files have the same number of samples
        for line in nums[1:]:
            if line[1] != nums[0][1] or line[2] != nums[0][2]:
                good = False
                break
        if good:
            s1 += 'No. %s = %.0f' % (__kde_g_label, nums[0][2])
            s2 += 'No. %s = %.0f' % (__kde_c_label, nums[0][1])
        elif labels:
            s1 += 'No. %s = ' % __kde_g_label
            s2 += 'No. %s = ' % __kde_c_label
            for n in nums:
                s1 += '%s:%.0f ' % (labels[n[0]], n[2])
                s2 += '%s:%.0f ' % (labels[n[0]], n[1])

        ax.text(__roc_xmin, __roc_ymax, s1 + '\n' + s2, ha='left',
                va='bottom')

    if len(nums) == 0:
        ax.text(__roc_xmin, __roc_ymax, 'no data found', ha='left',
                va='bottom')

    if write:
        if not figure_name:
            figure_name = 'roc'
        plt.savefig(figure_name + fig_type)

    return fig


#########################

def ROC_to_pwg_kde_plot(
    ROCfile,
    figure_name=False,
    write=True,
    label=True,
    num_samples=__kde_num_samples,
    ):
    """ 
  generate an overlay of kde estimates from an ROC file. 
  ONLY PLOTS FROM A SINGLE ROCfile!
  """

    (eval, c_kde, g_kde, n_cln, n_gch) = ROC_to_pwg_kde(ROCfile,
            num_samples=num_samples, scale=__kde_scale, s=__kde_s,
            filename=False)

#  print ROCfile

    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)
    ax.plot(eval, c_kde, color=__kde_c_color, label=__kde_c_label)
    ax.plot(eval, g_kde, color=__kde_g_color, label=__kde_g_label)

#  ax.set_yscale(__kde_yscale)

    ax.set_ylim(ymin=__kde_ymin)

#  ax.set_xscale(__kde_xscale)

    ax.set_xlim(xmin=__kde_xmin)

    ax.set_xlabel(__kde_xlabel)
    ax.set_ylabel(__kde_ylabel)

    ax.legend(loc=__kde_legend_loc)
    ax.grid(__kde_grid)

    if label:
        lims = ax.axis()
        ax.text(lims[0], lims[-1], 'No. %s = %.0f\nNo. %s = %.0f'
                % (__kde_g_label, float(n_gch), __kde_c_label,
                float(n_cln)), ha='left', va='bottom')

    if write:
        if not figure_name:
            figure_name = 'kde'
        plt.savefig(figure_name + fig_type)

    return fig


#########################

def ROC_to_pwd_kde_likelihood_plot(
    ROCfiles,
    figure_name=False,
    write=True,
    num_samples=__kde_num_samples,
    ):
    """
  generate a mapping between rank and likelihood and make a plot using an ROC file.
  ONLY PLOTS FROM A SINGLE ROCfile!
  """

    (eval, L) = ROC_to_pwg_kde_likelihood(ROCfile,
            num_samples=num_samples, scale=__kde_scale, s=__kde_s,
            filename=False)

  # take care if there are any "inf" values in L, set them to be 2*max(L)

    max_L = max(L)
    for ind in range(len(L)):
        if L[ind] == 'inf':
            L[ind] = 2 * max_L

    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)

    ax.plot(eval, L, color=__kde_L_color, label=__kde_L_label)

#  ax.set_yscale(__kde_L_yscale)

    ax.set_ylim(ymin=__kde_L_ymin)

#  ax.set_xscale(__kde_L_xscale)

    ax.set_xlim(xmin=__kde_L_xmin)

    ax.set_xlabel(__kde_L_xlabel)
    ax.set_ylabel(__kde_L_ylabel)

    ax.legend(loc=__kde_L_legend_loc)
    ax.grid(__kde_L_grid)

    if write:
        if not figure_name:
            figure_name = 'kde_L'
        plt.savefig(figure_name + fig_type)

    return fig


###################################################
#
#         Trending plots for classifiers
#
####################################################

def stat_to_trends_plot(
    statfiles,
    classifiers,
    labels=False,
    colors=False,
    output_dir='.',
    figure_basename=False,
    write=True,
    annotate=True,
    ):
    """
....Reads in a list of stat.txt files containg live time and other statistics. Generates plots of livetime vs time, glitch/clean rate vs time. 
...."""

    statfiles.sort()

    # set input parameters

    if labels:
        if len(labels) != len(classifiers):
            print 'labels and classifiers must have the same length for labels to be plotted.'
            labels = False
    if colors:
        if len(colors) != len(classifiers):
            print 'colors and classifiers must have the same length for colors to be applied'
            colors = False

    # define dictionaries

    livetime_dict = {}
    r_glitches_dict = {}
    r_cleans_dict = {}
    for classifier in classifiers:
        livetime_dict[classifier] = []
        r_glitches_dict[classifier] = []
        r_cleans_dict[classifier] = []

    times = []
    durations = []
    for filename in statfiles:
        time = float(filename.split('-')[-2]) + float(filename.split('-'
                )[-1].split('.')[0]) / 2.
        duration = float(filename.split('-')[-1].split('.')[0])
        times.append(time)
        durations.append(duration)
        file = open(filename, 'r')

        for line in file:
            if line.strip() in classifiers:
                section = line.strip()
            if 'Livetime = ' in line:
                livetime_dict[section].append(float(line.split('='
                        )[-1].strip()))
            if 'Glitches = ' in line:
                r_glitches_dict[section].append(float(line.split('='
                        )[-1].strip()))
            if 'Cleans = ' in line:
                r_cleans_dict[section].append(float(line.split('='
                        )[-1].strip()))

    # last gps time

    last_gps_time = times[-1]

    # compute time differences from the last central time

    times = np.asarray(times) - times[-1]

    # convert into days

    times /= 86400.0

    # compute rates and fractional livetimes

    for classifier in classifiers:
        for i in range(len(times)):
            if not livetime_dict[classifier][i] == 0:
                r_glitches_dict[classifier][i] /= \
                    livetime_dict[classifier][i]
                r_cleans_dict[classifier][i] /= \
                    livetime_dict[classifier][i]
                livetime_dict[classifier][i] /= durations[i]
            else:
                if r_glitches_dict[classifier][i] != 0 \
                    or r_cleans_dict[classifier][i] != 0:
                    print 'ERROR: ' \
                        + str(r_glitches_dict[classifier][i]) \
                        + ' glitches  and ' \
                        + str(r_cleans_dict[classifier][i]) \
                        + ' cleans are incosistent with zero livetime.'
                    sys.exit(1)

    fig_paths = []

    # generate livetime trending plot

    fig = plt.figure()
    for (i, classifier) in enumerate(classifiers):
        if colors and labels:
            plt.plot(
                times,
                100.0 * np.asarray(livetime_dict[classifier]),
                color=colors[i],
                label=labels[i],
                marker='o',
                markersize=10.0,
                linestyle='-',
                )
        elif colors:
            plt.plot(
                times,
                100.0 * np.asarray(livetime_dict[classifier]),
                color=colors[i],
                marker='o',
                markersize=10.0,
                linestyle='-',
                )
        elif labels:
            plt.plot(
                times,
                100.0 * np.asarray(livetime_dict[classifier]),
                label=labels[i],
                marker='o',
                markersize=10.0,
                linestyle='-',
                )
        else:
            plt.plot(times, 100.0
                     * np.asarray(livetime_dict[classifier]), marker='o'
                     , markersize=10.0, linestyle='-')

    plt.xlabel(__livetime_trend_xlabel)
    plt.ylabel(__livetime_trend_ylabel)
    plt.title('Livetime at GPS ' + str(last_gps_time))
    l = plt.legend(loc=__livetime_trend_legend_loc, fancybox=True)
    if l:
        l.get_frame().set_alpha(0.5)

    plt.grid(__livetime_trend_grid)

    if write:
        if not figure_basename:
            figure_basename = 'stat_trending_'
        figure_name = output_dir + '/' + figure_basename \
            + 'livetime.png'
        plt.savefig(figure_name)
    plt.close()

    fig_paths.append((figure_name, 'livetime'))

    # generate glitch rate trending plot

    fig = plt.figure()
    for (i, classifier) in enumerate(classifiers):
        if not all(r == 0.0 for r in r_glitches_dict[classifier]):
            if colors and labels:
                plt.plot(
                    times,
                    r_glitches_dict[classifier],
                    color=colors[i],
                    label=labels[i],
                    marker='o',
                    markersize=10.0,
                    linestyle='-',
                    )
            elif colors:
                plt.plot(
                    times,
                    r_glitches_dict[classifier],
                    color=colors[i],
                    marker='o',
                    markersize=10.0,
                    linestyle='-',
                    )
            elif labels:
                plt.plot(
                    times,
                    r_glitches_dict[classifier],
                    label=labels[i],
                    marker='o',
                    markersize=10.0,
                    linestyle='-',
                    )
            else:
                plt.plot(times, r_glitches_dict[classifier], marker='o'
                         , markersize=10.0, linestyle='-')

    plt.xlabel(__glitchrate_trend_xlabel)
    plt.ylabel(__glitchrate_trend_ylabel)
    plt.yscale(__glitchrate_trend_yscale)
    plt.title('Glitch rate at GPS ' + str(last_gps_time))
    l = plt.legend(loc=__glitchrate_trend_legend_loc, fancybox=True)
    if l:
        l.get_frame().set_alpha(0.5)

    plt.grid(__glitchrate_trend_grid)

    if write:
        if not figure_basename:
            figure_basename = 'stat_trending_'
        figure_name = output_dir + '/' + figure_basename \
            + 'glitchrate.png'
        plt.savefig(figure_name)
    plt.close()

    fig_paths.append((figure_name, 'glitch rate'))

    # generate clean rate trending plot

    fig = plt.figure()
    for (i, classifier) in enumerate(classifiers):
        if not all(r == 0.0 for r in r_cleans_dict[classifier]):
            if colors and labels:
                plt.plot(
                    times,
                    r_cleans_dict[classifier],
                    color=colors[i],
                    label=labels[i],
                    marker='o',
                    markersize=10.0,
                    linestyle='-',
                    )
            elif colors:
                plt.plot(
                    times,
                    r_cleans_dict[classifier],
                    color=colors[i],
                    marker='o',
                    markersize=10.0,
                    linestyle='-',
                    )
            elif labels:
                plt.plot(
                    times,
                    r_cleans_dict[classifier],
                    label=labels[i],
                    marker='o',
                    markersize=10.0,
                    linestyle='-',
                    )
            else:
                plt.plot(times, r_cleans_dict[classifier], marker='o',
                         markersize=10.0, linestyle='-')

    plt.xlabel(__cleanrate_trend_xlabel)
    plt.ylabel(__cleanrate_trend_ylabel)
    plt.yscale(__cleanrate_trend_yscale)
    plt.title('Clean rate at GPS ' + str(last_gps_time))
    l = plt.legend(loc=__cleanrate_trend_legend_loc, fancybox=True)
    if l:
        l.get_frame().set_alpha(0.5)

    plt.grid(__cleanrate_trend_grid)

    if write:
        if not figure_basename:
            figure_basename = 'stat_trending_'
        figure_name = output_dir + '/' + figure_basename \
            + 'cleanrate.png'
        plt.savefig(figure_name)
    plt.close()

    fig_paths.append((figure_name, 'clean rate'))

    return fig_paths


def ROC_to_eff_trends_plot(
    ROCfiles_dict,
    FAP,
    classifiers,
    labels=False,
    colors=False,
    output_dir='.',
    figure_basename=False,
    write=True,
    annotate=True,
    ):
    """ 
....Reads roc files and generates plots of Eff at given value of FAP vs time. ROCfiles_dict is a dictinary. Its keys are classfiers and its values are lists of roc files.
...."""

    # get list of classfiers
    # classifiers = ROCfiles_dict.keys()
    # classifiers.sort()

    # sort roc files for each classifier to ensure correct time ordering

    for classifier in classifiers:
        ROCfiles_dict[classifier].sort()

    # set input parameters

    if labels:
        if len(labels) != len(classifiers):
            print 'labels and ROCfiles must have the same length for labels to be plotted.'
            labels = False
    if colors:
        if len(colors) != len(classifiers):
            print 'colors and ROCfiles must have the same length for colors to be applied'
            colors = False

    # check that lists of roc files are all of the same length

    for classifier in classifiers[1:]:
        if len(ROCfiles_dict[classifier]) \
            != len(ROCfiles_dict[classifiers[0]]):
            print 'Error: Lists of roc files are not the same length. Files are probably missing for one or more classifiers.'
            sys.exit(1)

    # define dictionary of efficiencies

    eff_dict = {}
    for classifier in ROCfiles_dict.keys():
        eff_dict[classifier] = []

    # compute central gps times corresponding to roc files

    times = []

    # time periods between classifiers should be synchronized
    # use the first (alphabetically) classifier to get central gps times

    for ROCfile in ROCfiles_dict[classifiers[0]]:
        time = float(ROCfile.split('-')[-2]) + float(ROCfile.split('-'
                )[-1].split('.')[0]) / 2.
        times.append(time)

    # get last gps time

    last_gps_time = times[-1]

    # compute time differences form the last central time

    times = np.asarray(times) - times[-1]

    # convert into days

    times /= 86400.0

    # get efficiencies from roc files

    for classifier in classifiers:
        for ROCfile in ROCfiles_dict[classifier]:

            # get roc curves from roc files

            (faps, effs, tot_cln, tot_gch) = get_roc_from_file(ROCfile)

            # compute efficiency corresponding to fap = FAP

            if len(faps):
                eff = get_eff_at_fap_threshold(faps, effs, FAP)
            else:
                eff = 0.0
            eff_dict[classifier].append(eff)

    fig_paths = []

    # make a effciency trending plot

    fig = plt.figure()
    for (i, classifier) in enumerate(classifiers):

        if colors and labels:
            plt.plot(
                times,
                eff_dict[classifier],
                color=colors[i],
                label=labels[i],
                marker='o',
                markersize=10.0,
                linestyle='-',
                )
        elif colors:
            plt.plot(
                times,
                eff_dict[classifier],
                color=colors[i],
                marker='o',
                markersize=10.0,
                linestyle='-',
                )
        elif labels:
            plt.plot(
                times,
                eff_dict[classifier],
                label=labels[i],
                marker='o',
                markersize=10.0,
                linestyle='-',
                )
        else:
            plt.plot(times, eff_dict[classifier], marker='o',
                     markersize=10.0, linestyle='-')

    plt.xlabel(__eff_trend_xlabel)
    plt.ylabel(__eff_trend_ylabel)
    plt.title('Efficiency for FAP = ' + str(FAP) + ' at GPS '
              + str(last_gps_time))
    plt.grid(__eff_trend_grid)

    l = plt.legend(loc=__eff_trend_legend_loc, fancybox=True)
    if l:
        l.get_frame().set_alpha(0.5)

    if write:
        if not figure_basename:
            figure_basename = '_trending_'
        figure_name = output_dir + '/' + figure_basename \
            + 'eff_vs_time.png'
        plt.savefig(figure_name)
    plt.close()

    fig_paths.append((figure_name, 'efficiency vs time'))

    return fig_paths


##################################################

def chanlist_trending(
    gps_start,
    gps_stop,
    glob_dir,
    classifier='ovl',
    figure_name=False,
    verbose=False,
    annotated=True,
    ):
    """
  builds a channel performance trending plot using chanlist files
  """

    ax_bounds = [0.025, 0.075, 0.5, 0.775]
    ax_cb_bounds = [0.025, 0.875, 0.5, 0.05]

  # figure out appropriat time variable:

    dur = gps_stop - gps_start
    for (tau, tau_name) in [(1, 'seconds'), (60, 'minutes'), (3600,
                            'hours'), (86400, 'days'), (7 * 86400,
                            'weeks')]:
        if dur / tau < 100:
            break

  # ## find only the channel lists within this range

    if verbose:
        print 'finding %s*channel_summary.txt files within [%f, %f]' \
            % (classifier, gps_start, gps_stop)
    chanlists = []
    for _dir in sorted(glob.glob(glob_dir + '/*_*')):
        try:
            (_start, _stop) = [int(l) for l in _dir.split('/'
                               )[-1].split('_')]
            if gps_start <= _start and gps_stop > _start or gps_start \
                < _stop and gps_stop >= _stop:
                chanlists.append(glob.glob(_dir + '/' + classifier
                                 + '*_channel_summary.txt')[0])  # expect only one channel_summary.txt per _dir per classifier
        except:
            pass

  # ## read in data from chanlists and build dictionary

    chans = {}
    for chanlist in chanlists:
        if verbose:
            print 'reading data from %s' % chanlist

    # get time range for this chan list

        (_start, _dur) = [int(l) for l in chanlist.split('/'
                          )[-1].split('_')[0].split('-')[-2:]]  # eg: ovl-1043366400-86400_channel_summary.txt

    # read in this chanlist

        c = [[line[0]] + [float(l) for l in line[1:]] for line in
             [_line.strip().split() for _line in open(chanlist, 'r'
             ).readlines()[1:]]]

        for chan_dat in c:
            channel = chan_dat[chanlistD['channel']]
            if chans.has_key(channel):
                chans[channel].append([_start, _start + _dur,
                        chan_dat[chanlistD['eff/fap']]])
            else:
                chans[channel] = [[_start, _start + _dur,
                                  chan_dat[chanlistD['eff/fap']]]]

  # ## build plot

    if verbose:
        print 'generating figure'
    color_map = matplotlib.cm.get_cmap()
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    fig = plt.figure()
    ax = fig.add_axes(ax_bounds)

    yticklabels = []
    time_min = 0
    ind = 0
    for (ind, channel) in enumerate([chan for chan in
                                    sorted(chans.keys()) if chan
                                    != 'none']):  # we don't care about the "none" channel
        yticklabels.append(channel)
        for (_start, _stop, eff_fap) in chans[channel]:
            if _start - gps_stop < time_min:
                time_min = _start - gps_stop
            rank = 1 - np.exp(-eff_fap / 100.0)
            ax.fill_between([(_start - gps_stop) / tau, (_stop
                            - gps_stop) / tau], [ind + 0.49, ind
                            + 0.49], [ind - 0.49, ind - 0.49],
                            facecolor=color_map(rank), edgecolor='none')
            if annotated:
                ax.text(((_start + _stop) / 2. - gps_stop) / tau, ind,
                        '%.2f' % rank, ha='center', va='center')

    ax.set_xlim(xmin=time_min / tau, xmax=0)
    ax.set_ylim(ymin=-0.5, ymax=ind + 0.5)

    ax.set_xlabel('%s relative to gps=%.0f' % (tau_name, gps_stop))
    ax.set_yticks(range(ind + 1))
    ax.set_yticklabels(yticklabels)
    ax.yaxis.tick_right()

#  plt.setp(plt.getp(ax, 'yticklabels'), fontsize=7.5)
#  plt.setp(plt.getp(ax, 'xticklabels'), fontsize=10)

    ax_colorbar = fig.add_axes(ax_cb_bounds)
    cb = mpl.colorbar.ColorbarBase(ax_colorbar, cmap=color_map,
                                   norm=norm, orientation='horizontal')
    cb.set_label(classifier + ' rank')
    ax_colorbar.xaxis.tick_top()
    ax_colorbar.xaxis.set_label_position('top')

    plt.setp(fig, figheight=2. + 0.54 * (ind + 1), figwidth=10)
    ax.grid(True)

    if figure_name:
        if verbose:
            print 'saving figure into %s' % figure_name
        plt.savefig(figure_name)
        return figure_name
    else:
        return fig


#########################
#
# given a range of data, a directory, and a channel name
#
#########################

def KWtrg_rates(filenames, channel):
    print 'WRITE ME'
    return False


#########################

def KWtrg_signif_hist(filenames, channel):
    print 'WRITE ME'
    return False


#########################
#
#
#########################

def build_html(filename, classifiers):
    """ writes a summary hmtl file """

    print "you should probably write this, `cause it don't do nothing."


##@}


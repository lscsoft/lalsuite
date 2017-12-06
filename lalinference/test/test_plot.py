from __future__ import division
from __future__ import print_function
import matplotlib
matplotlib.use('agg')
import numpy as np
import matplotlib.pyplot as plt
import sys

from pkg_resources import require, DistributionNotFound, VersionConflict
try:
    from astropy.tests.helper import pytest
    require('pytest_mpl')
except (ImportError, DistributionNotFound, VersionConflict):
    print('these tests require pytest and pytest-mpl', file=sys.stderr)
    raise SystemExit(77)

import lalinference.plot


def pp_plot():
    # Re-initialize the random seed to make the unit test repeatable
    np.random.seed(0)
    fig = plt.figure(figsize=(3, 3), dpi=72)
    ax = fig.add_subplot(111, projection='pp_plot')
    p_values = np.arange(1, 20) / 20
    return fig, ax, p_values


@pytest.mark.mpl_image_compare
def test_pp_plot_steps():
    """Test P--P plot with drawstyle='steps'."""
    fig, ax, p_values = pp_plot()
    ax.add_confidence_band(len(p_values))
    ax.add_diagonal()
    ax.add_lightning(len(p_values), 20, drawstyle='steps')
    ax.add_series(p_values, drawstyle='steps')
    return fig


@pytest.mark.mpl_image_compare
def test_pp_plot_lines():
    """Test P--P plot with drawstyle='steps'."""
    fig, ax, p_values = pp_plot()
    ax.add_confidence_band(len(p_values))
    ax.add_diagonal()
    ax.add_lightning(len(p_values), 20, drawstyle='lines')
    ax.add_series(p_values, drawstyle='lines')
    ax.add_diagonal()
    return fig


@pytest.mark.mpl_image_compare
def test_pp_plot_default():
    """Test P--P plot with drawstyle='steps'."""
    fig, ax, p_values = pp_plot()
    ax.add_confidence_band(len(p_values))
    ax.add_diagonal()
    ax.add_lightning(len(p_values), 20)
    ax.add_series(p_values)
    return fig


@pytest.mark.mpl_image_compare
def test_mollweide_axes():
    """Test HEALPix heat map on 'mollweide' axes"""
    fig = plt.figure(figsize=(6, 4), dpi=72)
    fig.add_subplot(111, projection='mollweide')
    lalinference.plot.healpix_heatmap(np.arange(12), nest=True)
    return fig


@pytest.mark.mpl_image_compare
def test_astro_mollweide_axes():
    """Test HEALPix heat map on 'astro hours mollweide' axes"""
    fig = plt.figure(figsize=(6, 4), dpi=72)
    fig.add_subplot(111, projection='astro hours mollweide')
    lalinference.plot.healpix_heatmap(np.arange(12), nest=True)
    return fig


if __name__ == '__main__':
    raise SystemExit(pytest.main(['-vv', '--mpl', __file__]))

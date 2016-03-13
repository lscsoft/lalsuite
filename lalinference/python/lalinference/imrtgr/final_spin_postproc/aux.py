# Auxiliary functions needed for Mfaf_posterior.py
# (c) Archisman Ghosh, 2013-10-23 and 2015-11-17

import numpy as np
from scipy import interpolate

# Simple functions to translate between bin edges, bin centers, bin widths
def width(edg_or_cen):
  return np.average(np.diff(edg_or_cen))

def cen_from_edg(edges):
  return (edges[1:]+edges[:-1])/2.

# Basic module for confidence calculations
class confidence(object):
  def __init__(self, counts):
    self.counts = counts
    # Sort in descending order in frequency
    self.counts_sorted = np.sort(counts.flatten())[::-1]
    # Get a normalized cumulative distribution from the mode
    self.norm_cumsum_counts_sorted = np.cumsum(self.counts_sorted) / np.sum(counts)
    # Set interpolations between heights, bins and levels
    self._set_interp()
  def _set_interp(self):
    self._length = len(self.counts_sorted)
    # height from index
    self._height_from_idx = interpolate.interp1d(np.arange(self._length), self.counts_sorted, bounds_error=False, fill_value=0.)
    # index from height
    self._idx_from_height = interpolate.interp1d(self.counts_sorted[::-1], np.arange(self._length)[::-1], bounds_error=False, fill_value=self._length)
    # level from index
    self._level_from_idx = interpolate.interp1d(np.arange(self._length), self.norm_cumsum_counts_sorted, bounds_error=False, fill_value=1.)
    # index from level
    self._idx_from_level = interpolate.interp1d(self.norm_cumsum_counts_sorted, np.arange(self._length), bounds_error=False, fill_value=self._length)
  def level_from_height(self, height):
    """
    Returns value for a given confidence level.
    """
    return self._level_from_idx(self._idx_from_height(height))
  def height_from_level(self, level):
    """
    Returns confidence level for a given value
    """
    return self._height_from_idx(self._idx_from_level(level))
  def _idx_above(self, level):
    """
    Returns indices of values above a given confidence level.
    """
    return np.where(self.counts>self.height_from_level(level))

# 1d confidence
class confidence_1d(confidence):
  def __init__(self, counts, edges):
    confidence.__init__(self, counts)
    self.centers = cen_from_edg(edges)
    self.width = width(edges)
  def edges(self, level=0.68):
    """
    Returns edges of the confidence interval (assumes connected).
    """
    valid_x = self.centers[self._idx_above(level)]
    return (min(valid_x)-self.width/2., max(valid_x)+self.width/2.)
  def range(self, level=0.68):
    """
    Returns size of 1d confidence region.
    """
    return np.size(self._idx_above(level))*self.width

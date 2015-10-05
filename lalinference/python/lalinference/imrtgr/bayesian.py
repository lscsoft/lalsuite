""" 
Module that does various computations in Bayesian framework with probability distributions, priors, posteriors. 

(c) Archisman Ghosh, 2013-10-23 
"""

import math, numpy as np
from scipy import integrate
from scipy import interpolate
from scipy.stats import gaussian_kde

# Simple functions involving bin edges, bin centers, bin widths.

def width(edg_or_cen, uniform=False):
  if uniform:
    return edg_or_cen[1]-edg_or_cen[0]
  return np.average(edg_or_cen[1:]-edg_or_cen[:-1])

def edg2cen(edges):
  return (edges[1:]+edges[:-1])/2.

def cen2edg(centers, uniform=False):
  bin_width = width(centers, uniform)
  return  np.array([centers[0]-bin_width/2.]+list((centers[1:]+centers[:-1])/2.)+[centers[-1]+bin_width/2.])

def nsigma_value(counts, confidence=1):
  """
  Function to get level of confidence value. Takes eithe n=1 or confidence=0.68 for 1-sigma confidence value. Default is set to 1-sigma. 
  
  Parameters
  ----------
  counts : a frequency distribution. can me n-dimensional.
  n : 1, 2, 3 for 1-sigma, 2-sigma, 3-sigma ..
  confidence: 0.682, 0.957, 0.997 for 1-sigma, 2-sigma, 3-sigma ..
  
  Returns
  -------
  Value of frequency above which the fraction specified by confidence lie.
  """
  if confidence >= 1:
    confidence = math.erf(confidence/np.sqrt(2)) # Get confidence corresponding to n-sigma value
  counts_sorted = np.sort(counts.flatten())[::-1] # Sort in descending order in frequency
  cumsum_counts_sorted = np.cumsum(counts_sorted) # Get a cumulative distribution from the mode
  #ns_val =  (counts_sorted[np.where(cumsum_counts_sorted<confidence*np.sum(counts))][-1]+counts_sorted[np.where(cumsum_counts_sorted>confidence*np.sum(counts))][0])/2. # Value within which confidence fraction of the points lie
  ns_val = counts_sorted[np.where(cumsum_counts_sorted>confidence*np.sum(counts))][0]
  return ns_val


class prob_distr(object):
  def __init__(self, data, centers=None, bin_width=None, binned=False, bins=50, baseline=0., uniform=True, kde=None, **kwargs):
    if binned:
      (n, centers) = (data, centers)
      #data = np.array(reduce(operator.add, [[centers[i]]*n[i] for i in range(len(centers))]))
    elif centers is not None:
      bins = cen2edg(centers, uniform=True)
      (n, bins) = np.histogram(data, bins=bins, **kwargs)
    else:
      (n, bins) = np.histogram(data, bins=bins, **kwargs)
      centers = edg2cen(bins)
    if bin_width is None:
      bin_width = width(centers, uniform=True)
    self.prob = n + baseline
    self.centers = centers
    self._nbins = np.size(centers)
    self._bin_width = bin_width
    self._normalize()
    if not binned:
      self.kde = gaussian_kde(data)
    elif kde is not None:
      if callable(kde):
        self.kde = kde
      else:
        self.kde = lambda x: gaussian_kde
    else:
      self.kde = lambda x: 1./(self.lower_edge()-self.upper_edge())
    #self._normalize_kde()
  def lower_edge(self):
    return self.centers[0]-self._bin_width/2.
  def upper_edge(self):
    return self.centers[-1]+self._bin_width/2.
  def get_prob(self):
    return (self.prob, self.centers)
  def _normalize(self):
    self.prob = self.prob/np.sum(self.prob*self._bin_width)
  #def self._normalize_kde():
    #self.kde = self.kde
  def __add__(pd_1, pd_2):
    sum_kde = lambda x: pd_1.kde(x) + pd_2.kde(x)
    norm_sum = integrate.quad(sum_kde, pd_1.lower_edge(), pd_1.upper_edge())[0]
    norm_sum_kde = lambda x: pd_1.kde(x) + pd_2.kde(x) / norm_sum
    return prob_distr(pd_1.prob + pd_2.prob, pd_1.centers, pd_1._bin_width, binned=True, kde=norm_sum_kde)
  def __mul__(pd_1, pd_2):
    prod_kde = lambda x: pd_1.kde(x) * pd_2.kde(x)
    norm_prod = integrate.quad(prod_kde, pd_1.lower_edge(), pd_1.upper_edge())[0]
    norm_prod_kde = lambda x: pd_1.kde(x) * pd_2.kde(x) / norm_prod
    return prob_distr(pd_1.prob * pd_2.prob, pd_1.centers, pd_1._bin_width, binned=True, kde=norm_prod_kde)
  def spline(self, padding=10):
    padded_centers = np.append(np.append(self.centers[0]+(np.arange(padding)-padding)*self._bin_width, self.centers), self.centers[-1]+(1+np.arange(padding))*self._bin_width)
    padded_prob = np.append(np.append(np.zeros(padding), self.prob), np.zeros(padding))
    return interpolate.InterpolatedUnivariateSpline(padded_centers, padded_prob)
  def setspline(self, padding=10.):
    self._spline = self.spline(padding)
  def interpval(self, x):
    if self.centers[0] <= x <= self.centers[-1]:
      return self._spline(x)
    return 0.
  def smooth(self, n_smooth=20, spline_x=None):
    if spline_x is None:
      spline_x = np.linspace(np.min(self.centers)-self._bin_width/2., np.max(self.centers)+self._bin_width/2., n_smooth*np.size(self.centers))
    #prob_ispline = interpolate.InterpolatedUnivariateSpline(self.centers, self.prob)
    #return prob_distr(prob_ispline(spline_x), spline_x, binned=True)
    self.setspline()
    return prob_distr(self._spline(spline_x), spline_x, binned=True)
  def nsigma_value(self, confidence=1):
    return nsigma_value(self.prob, confidence)
  def nsigma_range(self, confidence=1):
    ns_val = self.nsigma_value(confidence)
    ns_cen = self.centers[np.where(self.prob>=ns_val)]
    #return (min(ns_cen), max(ns_cen))
    return (min(ns_cen), max(ns_cen))
  def nsigma_error(self, confidence=1):
    return np.diff(self.nsigma_range(confidence))[0]

class prob_distr_2d(object):
  def __init__(self, data_x, data_y=None, xcenters=None, ycenters=None, binned=False, bins=100, kde=None, **kwargs):
    if data_y is None:
      binned = True
      (counts, xcenters, ycenters) = (data_x, xcenters, ycenters)
    elif xcenters is not None:
      #print 'At the right place'
      xedges=cen2edg(xcenters, uniform=True)
      yedges=cen2edg(ycenters, uniform=True)
      (counts, xedges, yedges) = np.histogram2d(data_x, data_y, bins=[xedges, yedges], normed=True, **kwargs)
    else:
      (counts, xedges, yedges) = np.histogram2d(data_x, data_y, bins=bins, normed=True, **kwargs)
      xcenters = edg2cen(xedges)
      ycenters = edg2cen(yedges)
    self.prob = counts
    self.xcenters = xcenters
    self.ycenters = ycenters
    self._xbin_width = width(xcenters, uniform=True)
    self._ybin_width = width(ycenters, uniform=True)
    self._normalize()
    if not binned:
      self.kde = gaussian_kde([data_x, data_y])
    elif kde is not None:
      if callable(kde):
        self.kde = kde
      else:
        self.kde = lambda x: gaussian_kde
    else:
      self.kde = lambda x: 1./(self.xcenters[-1]-self.xcenters[0]+self._xbin_width)/(self.ycenters[-1]-self.ycenters[0]+self._ybin_width)
  def lower_xedge(self):
    return self.xcenters[0]-self._xbin_width/2.
  def upper_xedge(self):
    return self.xcenters[-1]+self._xbin_width/2.
  def lower_yedge(self):
    return self.ycenters[0]-self._ybin_width/2.
  def upper_yedge(self):
    return self.ycenters[-1]+self._ybin_width/2.
  def get_prob(self):
    return (self.prob, self.xcenters, self.ycenters)
  def _normalize(self):
    self.prob = self.prob/np.sum(self.prob*self._xbin_width*self._ybin_width)
  def __add__(pd_1, pd_2):
    return prob_distr_2d(pd_1.prob + pd_2.prob, xcenters=pd_1.xcenters, ycenters=pd_1.ycenters, binned=True, kde=lambda x: pd_1.kde(x) + pd_2.kde(x))
  def __mul__(pd_1, pd_2):
    return prob_distr_2d(pd_1.prob * pd_2.prob, xcenters=pd_1.xcenters, ycenters=pd_1.ycenters, binned=True, kde=lambda x: pd_1.kde(x) * pd_2.kde(x))
  def spline(self):
    return interpolate.RectBivariateSpline(self.xcenters, self.ycenters, self.prob)
  def nsigma_value(self, confidence=1):
    return nsigma_value(self.prob, confidence)
  def nsigma_area(self, confidence=1, xweights=None, yweights=None):
    ns_val = self.nsigma_value(confidence)
    if xweights is None:
      xweights = self.xcenters*0.+width(self.xcenters)
    elif callable(xweights):
      xweights = np.vectorize(xweights)(self.xcenters)*width(self.xcenters)
    if yweights is None:
      yweights = self.ycenters*0.+width(self.ycenters)
    elif callable(yweights):
      yweights = np.vectorize(yweights)(self.ycenters)*width(self.ycenters)
    (XX, YY) = np.meshgrid(self.xcenters, self.ycenters, indexing='ij')
    (Xw, Yw) = np.meshgrid(xweights, yweights, indexing='ij')
    wh = np.where(self.prob>=ns_val)
    return sum(Xw[wh]*Yw[wh])
  def sky_area(self, confidence=1):
    """
    Arguments have to be in the order (ra, dec).
    """
    return self.nsigma_area(confidence=confidence, yweights=np.cos)

class prob_distr_3d(object):
  def __init__(self, data_x, data_y, data_z, bins=100, **kwargs):
    (counts, (xedges, yedges, zedges)) = np.histogramdd((data_x, data_y, data_z), bins=bins, normed=True, **kwargs)
    xcenters = edg2cen(xedges)
    ycenters = edg2cen(yedges)
    zcenters = edg2cen(zedges)
    self.prob = counts
    self.xcenters = xcenters
    self.ycenters = ycenters
    self.zcenters = zcenters
    self._xbin_width = width(xcenters, uniform=True)
    self._ybin_width = width(ycenters, uniform=True)
    self._zbin_width = width(zcenters, uniform=True)
    self._normalize()
  def get_prob(self):
    return (self.prob, self.xcenters, self.ycenters, self.zcenters)
  def _normalize(self):
    self.prob = self.prob/np.sum(self.prob*self._xbin_width*self._ybin_width*self._zbin_width)
  def nsigma_value(self, confidence=1):
    return nsigma_value(self.prob, confidence)
  def nsigma_volume(self, confidence=1, xweights=None, yweights=None, zweights=None):
    ns_val = self.nsigma_value(confidence)
    if xweights is None:
      xweights = self.xcenters*0.+width(self.xcenters)
    elif callable(xweights):
      xweights = np.vectorize(xweights)(self.xcenters)*width(self.xcenters)
    if yweights is None:
      yweights = self.ycenters*0.+width(self.ycenters)
    elif callable(yweights):
      yweights = np.vectorize(yweights)(self.ycenters)*width(self.ycenters)
    if zweights is None:
      zweights = self.zcenters*0.+width(self.zcenters)
    elif callable(zweights):
      zweights = np.vectorize(zweights)(self.zcenters)*width(self.zcenters)
    (XX, YY, ZZ) = np.meshgrid(self.xcenters, self.ycenters, self.zcenters, indexing='ij')
    (Xw, Yw, Zw) = np.meshgrid(xweights, yweights, zweights, indexing='ij')
    wh = np.where(self.prob>=ns_val)
    return sum(Xw[wh]*Yw[wh]*Zw[wh])
  def sky_volume(self, confidence=1):
    """
    Arguments have to be in the order (ra, dec, dist).
    """
    return self.nsigma_volume(confidence=confidence, yweights=np.cos, zweights=(lambda z: z**2))

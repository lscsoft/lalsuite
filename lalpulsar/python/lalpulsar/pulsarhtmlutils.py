# -*- coding: utf-8 -*-
#
#       pulsarhtmlutils.py
#
#       Copyright 2016
#       Matthew Pitkin <matthew.pitkin@ligo.org>
#
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

"""
Some helper classes and functions for outputing html and LaTeX pages
"""

from __future__ import print_function

import re
import numpy as np
import math
import six

from lalpulsar.pulsarpputils import rad_to_dms, rad_to_hms

# some parameter names for special LaTeX treatment in figures
paramlatexdict = {'H0': '$h_0$',
                  'COSIOTA': '$\cos{\iota}$',
                  'PSI': '$\psi$ (rad)',
                  'PHI0': '$\phi_0$ (rad)',
                  'RA': '$\\alpha$',
                  'DEC': '$\delta$',
                  'RAJ': '$\\alpha$',
                  'DECJ': '$\delta$',
                  'F0': '$f_0$ (Hz)',
                  'F1': '$\dot{f}$ (Hz/s)',
                  'F2': '$\\ddot{f}$ (Hz/s$^2$)',
                  'F3': '$f_3$ (Hz/s$^3$)',
                  'F4': '$f_4$ (Hz/s$^4$)',
                  'F5': '$f_5$ (Hz/s$^5$)',
                  'F6': '$f_6$ (Hz/s$^6$)',
                  'F7': '$f_7$ (Hz/s$^7$)',
                  'F8': '$f_8$ (Hz/s$^8$)',
                  'F9': '$f_9$ (Hz/s$^9$)',
                  'F10': '$f_{10}$ (Hz/s$^{10}$)',
                  'LOGL': '$\log{L}$',
                  'PMRA': 'p.m. $\\alpha$ (rad/s)',
                  'PMDEC': 'p.m. $\delta$ (rad/s)',
                  'PMDC': 'p.m. $\delta$ (rad/s)',
                  'PX': '$\pi$ (rad)',
                  'A1': '$a \sin{i}$ (lt s)',
                  'A1_2': '$(a \sin{i})_{2}$ (lt s)',
                  'A1_3': '$(a \sin{i})_{3}$ (lt s)',
                  'SINI': '$\sin{i}$',
                  'PB': '$P_b$ (s)',
                  'PB_2': '$(P_b)_2$ (s)',
                  'PB_3': '$(P_b)_3$ (s)',
                  'T0': '$T_0$ (s)',
                  'T0_2': '$(T_0)_2$ (s)',
                  'T0_3': '$(T_0)_3$ (s)',
                  'TASC': '$T_{\\textrm{asc}}$ (s)',
                  'OM': '$\omega_0$ (rad)',
                  'OM_2': '$(\omega_0)_2$ (rad)',
                  'OM_3': '$(\omega_0)_3$ (rad)',
                  'PBDT': '$\dot{P}$ (s/s)',
                  'PBDOT': '$\dot{P}$ (s/s)',
                  'GAMMA': '$\gamma$',
                  'E': '$e$',
                  'ECC': '$e$',
                  'ECC_2': '$e_2$',
                  'ECC_3': '$e_3$',
                  'FB0': '$(f_b)_0$ (Hz)',
                  'FB1': '$(f_b)_1$ (Hz)',
                  'FB2': '$(f_b)_2$ (Hz)',
                  'FB3': '$(f_b)_3$ (Hz)',
                  'FB4': '$(f_b)_4$ (Hz)',
                  'FB5': '$(f_b)_5$ (Hz)',
                  'FB6': '$(f_b)_6$ (Hz)',
                  'FB7': '$(f_b)_7$ (Hz)',
                  'FB8': '$(f_b)_8$ (Hz)',
                  'FB9': '$(f_b)_9$ (Hz)',
                  'M2': '$m_2$ (kg)',
                  'MTOT': '$M$ (kg)',
                  'ELL': '$\\varepsilon$',
                  'H95': '$h_0^{95\%}$',
                  'H0UL': '$h_0^{{{}\%}}$',
                  'Q22': '$Q_{22}$\,(kg\,m$^2$)',
                  'SDRAT': 'spin-down ratio',
                  'SDRAT95': '$h_0^{95\%}$/h_0^{\\rm sd}$',
                  'SDLIM': '$h_0^{\\rm sd}$',
                  'F0ROT': '$f_{\\rm rot}$ (Hz)',
                  'F0GW': '$f_{\\rm gw}$ (Hz)',
                  'F1ROT': '$\dot{f}_{\\rm rot}$ (Hz/s)',
                  'F1GW': '$\dot{f}_{\\rm gw}$ (Hz/s)',
                  'SDPOWRAT': 'power ratio (\%)',
                  'OMDOT': '$\dot{\omega}$',
                  'OMDT': '$\dot{\omega}$',
                  'EPS1': '$\\epsilon_1$',
                  'EPS2': '$\\epsilon_2$',
                  'C22': '$C_{22}$',
                  'C21': '$C_{21}$',
                  'C22UL': '$C_{{22}}^{{{}\%}}$',
                  'C21UL': '$C_{{21}}^{{{}\%}}$',
                  'PHI22': '$\phi_{22}$',
                  'PHI21': '$\phi_{21}$',
                  'I31': '$I_{31}$',
                  'I21': '$I_{21}$',
                  'I31UL': '$I_{{31}}^{{{}\%}}$',
                  'I21UL': '$I_{{21}}^{{{}\%}}$',
                  'LAMBDA': '$\lambda$ (rad)',
                  'COSTHETA': '$\cos{\\theta}$',
                  'DIST': 'distance (kpc)',
                  'SNR': '$\\rho$',
                  'BSN': '$\log{}_{10}\left(B_{\\textrm{SvN}}\\right)$',
                  'BCI': '$\log{}_{10}\left(B_{\\textrm{CvI}}\\right)$',
                  'BCIN': '$\log{}_{10}\left(B_{\\textrm{CvIN}}\\right)$'}


# html text to display for different parameter names
paramhtmldict = {'RAJ': '&alpha;',
                 'DECJ': '&delta;',
                 'RA': '&alpha;',
                 'DEC': '&delta;',
                 'F0': 'f<sub>0</sub> (Hz)',
                 'F1': 'f<sub>1</sub> (Hz/s)',
                 'F2': 'f<sub>2</sub> (Hz/s<sup>2</sup>)',
                 'F3': 'f<sub>3</sub> (Hz/s<sup>3</sup>)',
                 'F4': 'f<sub>4</sub> (Hz/s<sup>4</sup>)',
                 'F5': 'f<sub>5</sub> (Hz/s<sup>5</sup>)',
                 'F6': 'f<sub>6</sub> (Hz/s<sup>6</sup>)',
                 'F7': 'f<sub>7</sub> (Hz/s<sup>7</sup>)',
                 'F8': 'f<sub>8</sub> (Hz/s<sup>8</sup>)',
                 'F9': 'f<sub>9</sub> (Hz/s<sup>9</sup>)',
                 'F10': 'f<sub>10</sub> (Hz/s<sup>10</sup>)',
                 'F0ROT': 'f<sub>rotation</sub> (Hz)',
                 'F1ROT': 'Spin-down<sub>rotation</sub> (Hz/s)',
                 'F0GW': 'f<sub>GW</sub> (Hz)',
                 'F1GW': 'Spin-down<sub>GW</sub> (Hz/s)',
                 'PEPOCH': 'epoch (MJD)',
                 'A1': 'a sin<it>i</i> (lt s)',
                 'A1_2': '(a sin<it>i</i>)<sub>2</sub> (lt s)',
                 'A1_3': '(a sin<it>i</i>)<sub>3</sub> (lt s)',
                 'SINI': 'sin<it>i</it>$',
                 'E': '<it>e</it>',
                 'ECC': '<it>e</it>',
                 'ECC_2': '<it>e</it><sub>2</sub>',
                 'ECC_3': '<it>e</it><sub>3</sub>',
                 'EPS1': '&epsilon;<sub>1</sub>',
                 'EPS2': '&epsilon;<sub>2</sub>',
                 'T0': 'T<sub>0</sub> (MJD)',
                 'T0_2': '(T<sub>0</sub>)<sub>2</sub> (MJD)',
                 'T0_3': '(T<sub>0</sub>)<sub>3</sub> (MJD)',
                 'TASC': 'T<sub>asc</sub> (MJD)',
                 'OM': '&omega;<sub>0</sub> (deg)',
                 'OM_2': '(&omega;<sub>0</sub>)<sub>2</sub> (deg)',
                 'OM_3': '(&omega;<sub>0</sub>)<sub>3</sub> (deg)',
                 'M2': '<it>m</it><sub>2</sub> (kg)',
                 'MTOT': '<it>M</it> (kg)',
                 'PB': '<it>P</it><sub>b</sub> (days)',
                 'PB_2': '(<it>P</it><sub>b</sub>)<sub>2</sub> (days)',
                 'PB_3': '(<it>P</it><sub>b</sub>)<sub>3</sub> (days)',
                 'FB0': '(<it>f</it><sub>b</sub>)<sub>0</sub> (Hz)',
                 'FB1': '(<it>f</it><sub>b</sub>)<sub>1</sub> (Hz)',
                 'FB2': '(<it>f</it><sub>b</sub>)<sub>2</sub> (Hz)',
                 'FB3': '(<it>f</it><sub>b</sub>)<sub>3</sub> (Hz)',
                 'FB4': '(<it>f</it><sub>b</sub>)<sub>4</sub> (Hz)',
                 'FB5': '(<it>f</it><sub>b</sub>)<sub>5</sub> (Hz)',
                 'FB6': '(<it>f</it><sub>b</sub>)<sub>6</sub> (Hz)',
                 'FB7': '(<it>f</it><sub>b</sub>)<sub>7</sub> (Hz)',
                 'FB8': '(<it>f</it><sub>b</sub>)<sub>8</sub> (Hz)',
                 'FB9': '(<it>f</it><sub>b</sub>)<sub>9</sub> (Hz)',
                 'H0': 'h<sub>0</sub>',
                 'C21': 'C<sub>21</sub>',
                 'C21UL': 'C<sub>21</sub><sup>{}%</sup>',
                 'C22': 'C<sub>22</sub>',
                 'C22UL': 'C<sub>22</sub><sup>{}%</sup>',
                 'I21': 'I<sub>21</sub>',
                 'I21UL': 'I<sub>21</sub><sup>{}%</sup>',
                 'I31': 'I<sub>31</sub>',
                 'I31UL': 'I<sub>31</sub><sup>{}%</sup>',
                 'COSIOTA': 'cos&iota;',
                 'PSI': '&psi; (rad)',
                 'PHI0': '&phi;<sub>0</sub> (rad)',
                 'PHI21': '&phi;<sub>21</sub> (rad)',
                 'PHI22': '&phi;<sub>22</sub> (rad)',
                 'PMRA': 'p.m. &alpha; (rad/s)',
                 'PMDC': 'p.m. &delta; (rad/s)',
                 'PMDEC': 'p.m. &delta; (rad/s)',
                 'PX': '&pi; (rad)',
                 'DIST': 'Distance (kpc)',
                 'SDLIM': 'Spin-down limit',
                 'ELL': '&#949;',
                 'SDRAT': 'ratio',
                 'H95': 'h<sub>0</sub><sup>95%</sup>',
                 'H0UL': 'h<sub>0</sub><sup>{}%</sup>',
                 'H0PRIOR': 'h<sub>0</sub><sup>95%</sup> prior',
                 'SDPOWRAT': 'power ratio (%)',
                 'Q22': 'Q<sub>22</sub> (kg m<sup>2</sup>)',
                 'BSN': 'log<sub>10</sub>(B<sub>SvN</sub>)',
                 'BCI': 'log<sub>10</sub>(B<sub>CvI</sub>)',
                 'BCIN': 'log<sub>10</sub>(B<sub>CvIN</sub>)',
                 'MAXL': 'log<sub>10</sub>(max. L)',
                 'SNR': '&rho;'}


# function to return a float number to a given number of significant figures
def sigfig(x, sf):
  return round(x, -int(math.floor(math.log10(abs(x))) - (sf - 1)))


# function will return a string with the number (input as a string) to two decimal
# places if it is greater than |0.01|, or the number in exponent for to one decimal
# place it it is smaller
def dec_or_exp(f, dp=2, horl='html'):
  fv = float(f)
  if np.abs(fv) > 0.01 and np.abs(fv) < 1000.:
    if horl == 'html':
      return repr(round(fv, dp))
    else:
      return '$%s$' % repr(round(fv, dp))
  else:
    return exp_str(fv, dp, otype=horl)


# a class containing functions to output html parameter values in the appropriate format
class paramhtmldispfunc:
  def RAJ(f, stype='hms'): # default to output in hh:mm:ss.s format
    if stype == 'hms':
      return ra_str(f)
    else: # output in radians
      return dec_or_exp(f)
  def RA(f, stype='hms'): # default to output in hh:mm:ss.s format
    if stype == 'hms':
      return ra_str(f)
    else: # output in radians
      return dec_or_exp(f)
  def DECJ(f, stype='dms'): # default to output in dd:mm:ss.s format
    if stype == 'dms':
      return dec_str(f)
    else: # output in radians
      return dec_or_exp(f)
  def DEC(f, stype='dms'): # default to output in dd:mm:ss.s format
    if stype == 'dms':
      return dec_str(f)
    else: # output in radians
      return dec_or_exp(f)
  def PMRA(f): return dec_or_exp(f)
  def PMDEC(f): return dec_or_exp(f)
  def PMDC(f): return dec_or_exp(f)
  def F0(f, dp=5): # default to returning with 5
    return dec_or_exp(f, dp=dp)
  def F1(f): return exp_str(float(f), 2)
  def F2(f): return exp_str(float(f), 2)
  def F3(f): return exp_str(float(f), 2)
  def F4(f): return exp_str(float(f), 2)
  def F5(f): return exp_str(float(f), 2)
  def F6(f): return exp_str(float(f), 2)
  def F7(f): return exp_str(float(f), 2)
  def F8(f): return exp_str(float(f), 2)
  def F9(f): return exp_str(float(f), 2)
  def F10(f): return exp_str(float(f), 2)
  def PEPOCH(f): return '%.1f' % float(repr(44244. + (float(f)+51.184)/86400.)) # return epoch as an float (converted from GPS to MJD)
  def A1(f): return dec_or_exp(f)
  def E(f): return dec_or_exp(f)
  def EPS1(f): return dec_or_exp(f)
  def EPS2(f): return dec_or_exp(f)
  def M2(f): return dec_or_exp(f)
  def MTOT(f): return dec_or_exp(f)
  def SINI(f): return dec_or_exp(f)
  def T0(f, stype=None):
    if stype == 'diff':
      return dec_or_exp(repr((float(f)+51.184)/86400.), dp=2)
    else:
      return '%.2f' % float(repr(44244. + ((float(f)+51.184)/86400.)))   # convert from GPS to MJD for display
  def TASC(f, stype=None):
    if stype == 'diff':
      return dec_or_exp(repr((float(f)+51.184)/86400.), dp=2)
    else:
      return '%.2f' % float(repr(44244. + ((float(f)+51.184)/86400.))) # convert from GPS to MJD for display
  def OM(f): return dec_or_exp(repr(float(f)*180./math.pi), dp=1)               # convert from rads to deg
  def PB(f): return dec_or_exp(repr(float(f)/86400.))                           # convert from seconds to days
  def H0(f): return exp_str(float(f), 1)
  def H0UL(f): return exp_str(float(f), 1)
  def C21(f): return exp_str(float(f), 1)
  def C21UL(f): return exp_str(float(f), 1)
  def C22(f): return exp_str(float(f), 1)
  def C22UL(f): return exp_str(float(f), 1)
  def I21(f): return exp_str(float(f), 1)
  def I21UL(f): return exp_str(float(f), 1)
  def I31(f): return exp_str(float(f), 1)
  def I31UL(f): return exp_str(float(f), 1)
  def COSIOTA(f): return dec_or_exp(f)
  def PHI0(f): return dec_or_exp(f)
  def PHI22(f): return dec_or_exp(f)
  def PHI21(f): return dec_or_exp(f)
  def PSI(f): return dec_or_exp(f)
  def ELL(f): return exp_str(float(f), 1)
  def SDLIM(f): return exp_str(float(f), 1)
  def SDRAT(f):
    fsf = sigfig(float(f), 2) # get value rounded to 2 significant figure
    if fsf < 1.: # if spin-down ratio is less than 1
      return '%.2f' % fsf
    elif fsf < 10.: # if spin-down ratio is less than 10
      return '%.1f' % fsf
    else: # otherwise round to the nearest integer
      return '%d' % round(fsf)
  def DIST(f): return '%.1f' % float(f)
  def SDPOWRAT(f): return '%d' % round(float(f))
  def Q22(f): return dec_or_exp(f, dp=1) # quadrupole moment
  def F0ROT(f): return dec_or_exp(f)
  def F0GW(f): return dec_or_exp(f)
  def F1ROT(f): return exp_str(float(f), 1)
  def F1GW(f): return exp_str(float(f), 1)
  def PMRA(f): return dec_or_exp(f)
  def PMDC(f): return dec_or_exp(f)
  def BSN(f): return dec_or_exp(f)
  def BCI(f): return dec_or_exp(f)
  def BCIN(f): return dec_or_exp(f)
  def DEFAULTSTR(f): return f


# a class for outputting parameter values to the table
class paramlatexdispfunc:
  def ELL(f): return exp_str(float(f), 1, 'latex')
  def H95(f): return exp_str(float(f), 1, 'latex')
  def H0(f): return exp_str(float(f), 1, 'latex')
  def H0UL(f): return exp_str(float(f), 1, 'latex')
  def C21(f): return exp_str(float(f), 1, 'latex')
  def C21UL(f): return exp_str(float(f), 1, 'latex')
  def C22(f): return exp_str(float(f), 1, 'latex')
  def C22UL(f): return exp_str(float(f), 1, 'latex')
  def I21(f): return exp_str(float(f), 1, 'latex')
  def I21UL(f): return exp_str(float(f), 1, 'latex')
  def I31(f): return exp_str(float(f), 1, 'latex')
  def I31UL(f): return exp_str(float(f), 1, 'latex')
  def H0PRIOR(f): return exp_str(float(f), 1, 'latex')
  def SDLIM(f): return exp_str(float(f), 1, 'latex')
  def SDRAT(f):
    fsf = sigfig(float(f), 2) # get value rounded to 2 significant figure
    if fsf < 1.: # if spin-down ratio is less than 1
      return '%.2f' % fsf
    elif fsf < 10.: # if spin-down ratio is less than 10
      return '%.1f' % fsf
    else: # otherwise round to the nearest integer
      return '%d' % fsf
  def RAJ(f): return ra_str(f, 'latex') # RA in string format
  def DECJ(f): return dec_str(f, 'latex') # dec in string format
  def RA(f): return ra_str(f, 'latex') # RA in string format
  def DEC(f): return dec_str(f, 'latex') # dec in string format
  def PMRA(f): return dec_or_exp(f, horl='latex')
  def PMDEC(f): return dec_or_exp(f, horl='latex')
  def PMDC(f): return dec_or_exp(f, horl='latex')
  def M2(f): return dec_or_exp(f, horl='latex')
  def MTOT(f): return dec_or_exp(f, horl='latex')
  def DIST(f): return '%.1f' % float(f)
  def SDPOWRAT(f): return '%d' % int(f)
  def Q22(f): return exp_str(float(f), 1, 'latex') # quadrupole moment
  def F0(f): return '%.2f' % float(f)
  def F0ROT(f): return '%.2f' % float(f)
  def F0GW(f): return '%.2f' % float(f)
  def F1(f): return exp_str(float(f), 1, 'latex')
  def F2(f): return exp_str(float(f), 1, 'latex')
  def F3(f): return exp_str(float(f), 1, 'latex')
  def F4(f): return exp_str(float(f), 1, 'latex')
  def F5(f): return exp_str(float(f), 1, 'latex')
  def F6(f): return exp_str(float(f), 1, 'latex')
  def F7(f): return exp_str(float(f), 1, 'latex')
  def F8(f): return exp_str(float(f), 1, 'latex')
  def F9(f): return exp_str(float(f), 1, 'latex')
  def F10(f): return exp_str(float(f), 1, 'latex')
  def F1ROT(f): return exp_str(float(f), 1, 'latex')
  def F1GW(f): return exp_str(float(f), 1, 'latex')
  def BSN(f): return dec_or_exp(f, horl='latex')
  def BCI(f): return dec_or_exp(f, horl='latex')
  def BCIN(f): return dec_or_exp(f, horl='latex')
  def DEFAULTSTR(f): return f


class htmltag:
  """
  A class to create a html tag
  """
  def __init__(self, tag, tagtext="", tagclass="", tagid="", tagstyle="", newline=False):
    self._taginfo = {}
    self.set_tag(tag)           # the name of the tag
    self.set_tagtext(tagtext)    # the text to go into the tag
    self.set_tagclass(tagclass) # the tag class
    self.set_tagid(tagid)       # the tag id
    self.set_tagstyle(tagstyle) # the tag style
    self.set_tagextra("")       # any additional formatting

    if newline:
      self._newline = '\n'  # whether to add a new line after the tag
    else:
      self._newline = ''

    # set the tag format
    self._tagdata = '<{tag} class="{class}" id="{id}" style="{style}"{extra}>{text}</{tag}>'+self._newline

  def __iadd__(self, ttext):
    """
    Overload the += operator to append text to tagtext
    """
    if not isinstance(ttext, str) and not isinstance(ttext, unicode):
      raise ValueError("Error... appended text must be a string.")
    else:
      self.set_tagtext(self._taginfo['text']+ttext)
    return self

  @property
  def tag(self):
    return self._taginfo['tag']

  def set_tag(self, t):
    if not isinstance(t, str) and not isinstance(t, unicode):
      raise ValueError("Error... 'tag' must be a string.")
    else:
      self._taginfo['tag'] = t

  @property
  def tagclass(self):
    return self._taginfo['class']

  def set_tagclass(self, tclass):
    if not isinstance(tclass, str) and not isinstance(tclass, unicode):
      raise ValueError("Error... 'class' text must be a string.")
    else:
      self._taginfo['class'] = tclass

  @property
  def tagid(self):
    return self._taginfo['id']

  def set_tagid(self, tid):
    if not isinstance(tid, str) and not isinstance(tid, unicode):
      raise ValueError("Error... 'id' text must be a string.")
    else:
      self._taginfo['id'] = tid

  @property
  def tagstyle(self):
    return self._taginfo['style']

  def set_tagstyle(self, tstyle):
    if not isinstance(tstyle, str) and not isinstance(tstyle, unicode):
      raise ValueError("Error... 'style' text must be a string.")
    else:
      self._taginfo['style'] = tstyle

  @property
  def tagtext(self):
    return self._taginfo['text']

  def set_tagtext(self, ttext):
    if not isinstance(ttext, str) and not isinstance(ttext, unicode):
      raise ValueError("Error... tag text must be a string.")
    else:
      self._taginfo['text'] = ttext

  @property
  def tagextra(self):
    return self._taginfo['extra']

  def set_tagextra(self, textra):
    if not isinstance(textra, str) and not isinstance(textra, unicode):
      raise ValueError("Error... 'extra' tag text must be a string.")
    else:
      space = '' # add no space
      if len(textra) > 0:
        space = ' ' # add a space
      self._taginfo['extra'] = space + textra

  @property
  def taginfo(self):
    return self._taginfo

  @property
  def text(self):
    # return the full tag
    return self._tagdata.format(**self.taginfo)

  def __str__(self):
    return self.text


class atag(htmltag):
  """
  Class for a link tag
  """
  def __init__(self, link, linktext="", linkclass="", linkid="", linkstyle=""):
    """
    Input the link and the text that the link surrounds
    """
    if linktext == "": # if no link text is given then just use the link itself
      linktext = link

    htmltag.__init__(self, 'a', linktext, linkclass, linkid, linkstyle)
    self.set_tagextra('href="{}"'.format(link)) # add href


class htmltable(htmltag):
  """
  Class to make and return a html table
  """
  def __init__(self, tag='table', tableclass="", tableid="", tablestyle=""):
    htmltag.__init__(self, 'table', tagclass=tableclass, tagid=tableid, tagstyle=tablestyle)

    self._rows = []
    self._nrows = 0 # number of rows
    self._thisrow = -1 # the index of the current row given row additions

  @property
  def tabletext(self):
    innertable = ""
    for row in self._rows:
      rowtxt = ""
      for data in row['data']:
        td = 'td'
        if data['header']:
          td = 'th' # use header <th> tags for this value
        datatag = htmltag(td, data['text'], data['class'], data['id'], data['style'])
        datatag.set_tagextra('rowspan="{rowspan}" colspan="{colspan}"'.format(**data))
        rowtxt += datatag.text + " " # add space between <td> elements
      rowtag = htmltag('tr', rowtxt, row['class'], row['id'], row['style'], newline=True)
      innertable += rowtag.text
    self.set_tagtext(innertable)
    return self.text

  def addrow(self, rowclass="", rowid="", rowstyle=""):
    """
    Add a new empty row dictionary to the list and increment the current row index
    """
    row = {'data': [], 'class': rowclass, 'id': rowid, 'style': rowstyle}
    self._rows.append(row)
    self._nrows += 1   # number of rows
    self._thisrow += 1 # the index of the row that has just been added

  def deleterow(self, rowidx):
    """
    Delete a row
    """
    if rowidx > self._nrows-1:
      print("Warning... cannot delete row '%d'. Only %d row in table." % (rowdix, self._nrows))
    else:
      self._rows.pop(rowidx) # remove row
      self._nrows -= 1
      self._thisrow -= 1

  def adddata(self, datatext, dataclass="", dataid="", datastyle="", header=False, rowspan=0, colspan=0, rowidx=None):
    """
    Add table data <td> (or <th> is header is True) tags to a given row
    """
    if rowidx is None:
      rowidx = self._thisrow

    if rowidx > len(self._rows)-1:
      raise ValueError("Warning... row index is out of range.")
    else:
      td = {'text': datatext, 'class': dataclass, 'id': dataid, 'style': datastyle, 'header': header, 'rowspan': '', 'colspan': ''}
      if rowspan > 0:
        td['rowspan'] = str(int(rowspan))
      if colspan > 0:
        td['colspan'] = str(int(colspan))

      self._rows[rowidx]['data'].append(td)

  def __str__(self):
    return self.tabletxt


class latextable:
  """
  Class to make a return a LaTeX table
  """
  def __init__(self, ncolumns=1, columnalign='c', caption="", label="", floatval="h", preamble="", postamble=""):
    """
    Create a table environment with `ncolumns` columns positioned with `columnpos`
    """
    self._tableinfo = {} # dictionary containing the table data

    self.set_ncolumns(ncolumns)       # set number of columns
    self.set_columnalign(columnalign) # set column text alignment
    self.set_caption(caption)         # set table caption
    self.set_label(label)             # set table label
    self.set_floatval(floatval)       # set floating environment position
    self.set_preamble(preamble)       # set any preamble before the tabular environment
    self.set_postamble(postamble)     # set any text to go after the tabular environment

    self._rows = [] # a list of row data
    self._nrows = 0 # number of rows
    self._thisrow = -1 # the index of the current row given row additions

    # set the table format
    self._tableformat = '\\begin{{table}}{{{floatval}}}\n{preamble}\caption{{{caption}\label{{{label}}}}}\n\\begin{{tabular}}{{{columnalign}}}\n{table}\n\\end{{tabular}}\n{postamble}\n\\end{{table}}'
    self._tableinfo['data'] = ""

  def set_ncolumns(self, ncolumns):
    # set number of columns in table
    self._ncolumns = ncolumns

  def set_columnalign(self, columnalign):
    # set the alignment of data within each column (if a list then this should be equal in length to ncolumns)
    if isinstance(columnalign, list):
      if len(columnalign) != self._ncolumns:
        raise ValueError("Error... number of column alignments is not equal to the number of columns.")
      else:
        for ca in columnalign:
          if not isinstance(ca, str) and not isinstance(ca, unicode):
            raise TypeError("Error... columnalign must be a list of strings.")

        self._columnalign = columnalign
    else:
      if isinstance(columnalign, str) or isinstance(columnalign, unicode):
        if len(columnalign) == 1:
          self._columnalign = [columnalign for _ in range(self._ncolumns)]
        else:
          self._columnalign = [columnalign]
      else:
        raise TypeError("Error... columnalign must be a list or a string.")

    self._tableinfo['columnalign'] = ' '.join(self._columnalign)

  def set_label(self, label):
    # set the label
    self._tableinfo['label'] = label

  @property
  def label(self):
    return self._tableinfo['label']

  def set_caption(self, caption):
    # set the table caption
    self._tableinfo['caption'] = caption

  @property
  def caption(self):
    return self._tableinfo['caption']

  def set_floatval(self, floatval):
    # set the floating environment position
    self._tableinfo['floatval'] = floatval

  def set_preamble(self, preamble):
    # set any preamble before the tabular environment
    self._tableinfo['preamble'] = preamble

  @property
  def preamble(self):
    return self._tableinfo['preamble']

  def set_postamble(self, postamble):
    # set any text to go after the tabular environment
    self._tableinfo['postamble'] = postamble

  @property
  def postamble(self):
    return self._tableinfo['postamble']

  def addrow(self, underline=False):
    # add an empty row (underline will add a horizontal rule below the row)
    row = {'data': [], 'underline': underline}
    self._rows.append(row)
    self._nrows += 1   # number of rows
    self._thisrow += 1 # the index of the row that has just been added

  def addhrule(self, rowidx=None):
    # add horizontal rule
    if rowidx == None: # add horizontal rule as new row
      self.addrow(underline=True)
    else:
      self._rows[rowidx]['underline'] = True

  def adddata(self, datatxt, multicolumn=0, mcalign='c', rowidx=None):
    # add data to a row
    if rowidx is None:
      rowidx = self._thisrow

    if rowidx > len(self._rows)-1:
      raise ValueError("Warning... row index is out of range.")
    else:
      rowdata = {'text': datatxt}
      if multicolumn != 0:
        rowdata['multicolumn'] = multicolumn
        rowdata['mcalign'] = mcalign

      self._rows[rowidx]['data'].append(rowdata)

  @property
  def tabletext(self):
    # output the full table
    self._tableinfo['table'] = ""
    for i, row in enumerate(self._rows):
      rowtxt = []
      ncols = 0
      for data in row['data']:
        val = data['text']
        if isinstance(val, float) or isinstance(val, int): # if a number convert into a string
          val = repr(val)
        if 'multicolumn' in data:
          ncols += data['multicolumn']
          rowtxt.append('\multicolumn{%d}{%s}{%s} ' % (data['multicolumn'], data['mcalign'], val))
        else:
          ncols += 1
          rowtxt.append(val+' ')
      if ncols != self._ncolumns and (len(rowtxt) != 0 and row['underline'] == False):
        raise ValueError("Error... too many or too few inputs in row '%d'." % i)
      if len(rowtxt) != 0:
        self._tableinfo['table'] += '&'.join(rowtxt) + '\\\\\n'
      if row['underline']: # add horizontal rule below row
        self._tableinfo['table'] += '\hline\n'

    self._tabletext = self._tableformat.format(**self._tableinfo)
    return self._tabletext


# convert a floating point number into a string in X.X x 10^Z format
def exp_str(f, p=1, otype='html'):
  if p > 16:
    print("Precision must be less than 16 d.p.", file=sys.stderr)
    p = 16

  s = '%.16e' % f
  ssplit = s.split('e')
  if otype.lower() == 'html': # output html format
    return '%.*f&times;10<sup>%d</sup>' % (p, float(ssplit[0]), int(ssplit[1]))
  elif otype.lower() == 'latex': # output LaTeX format
    return '\\ensuremath{%.*f\!\\times\!10^{%d}}'  % (p, float(ssplit[0]), int(ssplit[1]))
  else:
    raise ValueError("Error... 'otype' must be 'html' or 'latex'.")


# convert a right ascension string in format 'hh:mm:ss.s' to a html/LaTeX string like H^h M^m S^s.ss
def ra_str(ra, otype='html'):
  if isinstance(ra, six.string_types):
    hms = ra.split(":")
  elif isinstance(ra, float):
    hms = [str(v) for v in rad_to_hms(ra)]
  else:
    raise ValueError("Error... ra must be a string or a float.")

  if len(hms) == 1:
    hms.append('0')
    hms.append('0')
  elif len(hms) == 2:
    hms.append('0')

  ss = ('%.2f' % float(hms[2])).split('.')

  if otype.lower() == 'html': # return html string
    return "%s<sup>h</sup>%s<sup>m</sup>%s<sup>s</sup>.%s" % (hms[0].zfill(2), hms[1].zfill(2), ss[0].zfill(2), ss[1].zfill(2) )
  elif otype.lower() == 'latex': # return LaTeX string
    return "$%s^{\\rm h}%s^{\\rm m}%s^{\\rm s}\!.%s$" % (hms[0].zfill(2), hms[1].zfill(2), ss[0].zfill(2), ss[1].zfill(2) )
  else:
    raise ValueError("Error... 'otype' input must be 'html' or 'latex'")


# convert a declination string in format 'dd:mm:ss.s' to a html/LaTeX string like dd^o mm' ss''.ss
def dec_str(dec, otype='html'):
  if isinstance(dec, six.string_types):
    dms = dec.split(":")
  elif isinstance(dec, float):
    dms = [str(v) for v in rad_to_dms(dec)]
  else:
    raise ValueError("Error... dec must be a string or a float.")

  if len(dms) == 1:
    dms.append('0')
    dms.append('0')
  elif len(dms) == 2:
    dms.append('0')

  ss = ('%.2f' % float(dms[2])).split('.')

  if otype.lower() == 'html': # html output
    return "%s&deg;%s'%s\".%s" % ((re.sub('\+', '', dms[0])).zfill(2), dms[1].zfill(2), ss[0].zfill(2), ss[1].zfill(2))
  elif otype.lower() == 'latex': # LaTeX output
    return "$%s^{\circ}%s'%s''\!.%s$" % ((re.sub('\+', '', dms[0])).zfill(2), dms[1].zfill(2), ss[0].zfill(2), ss[1].zfill(2))
  else:
    raise ValueError("Error... 'otype' must be 'html' or 'latex'.")


"""
CSS files for results pages
"""

# css file for the standard individual pulsar result page
result_page_css = """
/* create body style */
body {
  font-family: "Avant Garde", Avantegarde, Verdana, Geneva, "Trebuchet MS", sans-serif;
}

/* create header name style */
h1 {
  margin: 0px 0px 0px 0px;
  padding: 4px 4px 8px 4px;
  font-size: 20px;
  font-weight: bold;
  letter-spacing: 0px;
  font-family: "Avant Garde", Avantegarde, Verdana, Geneva, "Trebuchet MS", Sans-Serif;
  background-color: darkolivegreen;
}

h1 > a:link {
  color: white;
  text-shadow: 2px 2px 2px #1c2310;
  text-decoration: none;
}

h1 > a:visited {
  color: white;
  text-shadow: 2px 2px 2px #1c2310;
  text-decoration: none;
}

h1 > a:hover {
  color: white;
  text-shadow: 2px 2px 2px #7fa046;
  text-decoration: none;
}

h2 {
  margin: 0 0 8px 0;
  padding: 4px 4px 8px 4px;
  font-size: 16px;
  font-weight: bold;
  font-family: "Avant Garde", Avantegarde, Verdana, Geneva, "Trebuchet MS", Sans-Serif;
  text-shadow: 1px 1px 1px #333300;
  background-color: #666600;
  color: white
}

/* create footer style */
#footer {
  border-top: 1px solid #999;
  padding: 15px;
  font-family: monospace;
  text-align: left;
}

/* create a class for a posterior plot image */
.posplot{
  height: 400px;
  border: 0px solid #999;
}

/* create a class for a full joint posterior plot image */
.jointplot{
  height: 600px;
  border: 0px solid #999;
}

/* create a class for a background distribution plot image */
.backgroundplot{
  height: 400px;
  border: 0px solid #999;
}

/* create a class for a Bk data plot */
.dataplot{
  height: 275px;
}

/* create class for an amplitude spectral density plot */
.asdplot{
  height: 275px;
}

/* create class for an MCMC chain plot */
.chainplot{
  height: 275px;
}

/* style for links list */
.pagelinks {
  background-color: darkolivegreen;
  font-family: "Avant Garde", Avantegarde, Verdana, Geneva, "Trebuchet MS", Sans-Serif;
  overflow: hidden;
  color: white;
  font-size: 16px;
  padding: 0px 0px 3px 3px;
  margin: 0px 0px 8px 0px;
  text-shadow: 2px 2px 2px #1c2310;
}

div.pagelinks a:link {
  color: white;
  text-shadow: 2px 2px 2px #1c2310;
  text-decoration: none;
}

div.pagelinks a:visited {
  color: white;
  text-shadow: 2px 2px 2px #1c2310;
  text-decoration: none;
}

div.pagelinks a:hover {
  color: white;
  text-shadow: 2px 2px 2px #7fa046;
  text-decoration: none;
}

/* pulsar parameter table class */
.pulsartable {
  background-color: floralwhite;
  border: 0px solid;
  border-radius: 4px;
  box-shadow: 2px 2px 2px 2px #d8d8d8;
  -webkit-box-shadow: 2px 2px 2px 2px #d8d8d8;
  -moz-box-shadow: 2px 2px 2px 2px #d8d8d8;
  padding: 4px 4px 4px 4px;
}

/* upper limits table class */
.limitstable {
   background-color: floralwhite;
   border: 0px solid;
   border-radius: 4px;
   box-shadow: 2px 2px 2px 2px #d8d8d8;
   -webkit-box-shadow: 2px 2px 2px 2px #d8d8d8;
   -moz-box-shadow: 2px 2px 2px 2px #d8d8d8;
   padding: 4px 4px 4px 4px;
}

/* background evidence table class */
.evidencetable {
  background-color: floralwhite;
  border: 0px solid;
  border-radius: 4px;
  box-shadow: 2px 2px 2px 2px #d8d8d8;
  -webkit-box-shadow: 2px 2px 2px 2px #d8d8d8;
  -moz-box-shadow: 2px 2px 2px 2px #d8d8d8;
  padding: 4px 4px 4px 4px;
}

/* set defaults for table data and table headers */
table{
  border: 0px;
  border-collapse: collapse;
}

td{
  padding: 0px 8px 0px 8px;
  border: 0px;
}

th{
  padding: 0px 8px 0px 8px;
  border: 0px;
}

.leftborder{
  border-left: 1px solid #000;
}

.rightborder{
  border-right: 1px solid #000;
}

.topborder{
  border-top: 1px solid #000;
}

.bottomborder{
  border-bottom: 1px solid #000;
}

/* set text colour classes for detectors */
.H1{
  color: red;
}

.H2{
  color: cyan;
}

.L1{
  color: green;
}

.V1{
  color: blue;
}

.G1{
  color: magenta;
}

.Joint{
  color: black;
  font-weight: bold;
}
"""


# create css file text for the tables of all results
results_table_css = """
/* create body style */
body {
  font-family: Verdana, Geneva, "Trebuchet MS", sans-serif;
}

/* create footer style */
#footer {
  border-top: 1px solid #999;
  padding: 15px;
  font-family: monospace;
  text-align: left;
}

/* create link style */
a:link{
  color: #000000;
  text-decoration: none;
}

a:visited{
  color: #000000;
  text-decoration: none;
}

a:hover{
  color: #000000;
  text-decoration: none;
  text-shadow: 2px 2px 2px #ccc;
}

/* set defaults for table data and table headers */
table{
  border: 0px;
  border-collapse: collapse;
}

td{
  padding: 0px 8px 0px 8px;
  border: 0px;
}

th{
  padding: 0px 8px 0px 8px;
  border: 0px;
}

.leftborder{
  border-left: 1px solid #000;
}

.rightborder{
  border-right: 1px solid #000;
}

.topborder{
  border-top: 1px solid #000;
}

.bottomborder{
  border-bottom: 1px solid #000;
}

/* set text colour classes for detectors */
.H1{
  color: red;
}

.H2{
  color: cyan;
}

.L1{
  color: green;
}

.V1{
  color: blue;
}

.G1{
  color: magenta;
}

.Joint{
  color: black;
  font-weight: bold;
}
"""

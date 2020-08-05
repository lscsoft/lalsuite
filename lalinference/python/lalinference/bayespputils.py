# -*- coding: utf-8 -*-
#
#       bayespputils.py
#
#       Copyright 2010
#       Benjamin Aylott <benjamin.aylott@ligo.org>,
#       Benjamin Farr <bfarr@u.northwestern.edu>,
#       Will M. Farr <will.farr@ligo.org>,
#       John Veitch <john.veitch@ligo.org>,
#       Salvatore Vitale <salvatore.vitale@ligo.org>,
#       Vivien Raymond <vivien.raymond@ligo.org>
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

#===============================================================================
# Preamble
#===============================================================================

"""
This module contains classes and functions for post-processing the output
of the Bayesian parameter estimation codes.
"""

#standard library imports
import os
import sys
from math import cos,ceil,floor,sqrt,pi as pi_constant
from xml.dom import minidom
from operator import itemgetter

#related third party imports
import matplotlib
matplotlib.use('agg')
from .io import read_samples
import healpy as hp
import astropy.table
import numpy as np
np.random.seed(42)
from numpy import fmod
from matplotlib import pyplot as plt,lines as mpl_lines
from scipy import stats
from scipy import special
from scipy import signal
from scipy.optimize import newton
from scipy import interpolate
from scipy import integrate
from numpy import linspace
import random
import socket
from itertools import combinations
from .lalinference import LALInferenceHDF5PosteriorSamplesDatasetName as posterior_grp_name
import re
import six

try:
    import lalsimulation as lalsim
except ImportError:
    print('Cannot import lalsimulation SWIG bindings')
    raise

try:
    from .imrtgr.nrutils import bbh_average_fits_precessing
except ImportError:
    print('Cannot import lalinference.imrtgr.nrutils. Will suppress final parameter and peak luminosity calculations.')

from matplotlib.ticker import ScalarFormatter

try:
    hostname_short=socket.gethostbyaddr(socket.gethostname())[0].split('.',1)[1]
except:
    hostname_short='Unknown'
if hostname_short=='ligo.caltech.edu' or hostname_short=='cluster.ldas.cit': #The CIT cluster has troubles with the default 'cm' font. 'custom' has the least troubles, but does not include \odot
    matplotlib.rcParams.update(
                               {'mathtext.fontset' : "custom",
                               'mathtext.fallback_to_cm' : True
                               })

from xml.etree.cElementTree import Element, SubElement, tostring, XMLParser

#local application/library specific imports
import lal
from . import git_version

__author__="Ben Aylott <benjamin.aylott@ligo.org>, Ben Farr <bfarr@u.northwestern.edu>, Will M. Farr <will.farr@ligo.org>, John Veitch <john.veitch@ligo.org>, Vivien Raymond <vivien.raymond@ligo.org>"
__version__= "git id %s"%git_version.id
__date__= git_version.date

def replace_column(table, old, new):
    """Workaround for missing astropy.table.Table.replace_column method,
    which was added in Astropy 1.1.

    FIXME: remove this function when LALSuite depends on Astropy >= 1.1."""
    index = table.colnames.index(old)
    table.remove_column(old)
    table.add_column(astropy.table.Column(new, name=old), index=index)

def as_array(table):
    """Workaround for missing astropy.table.Table.as_array method,
    which was added in Astropy 1.0.

    FIXME: remove this function when LALSuite depends on Astropy >= 1.0."""
    try:
        return table.as_array()
    except:
        return table._data

#===============================================================================
# Constants
#===============================================================================
#Parameters which are not to be exponentiated when found
logParams=['logl','loglh1','loglh2','logll1','loglv1','deltalogl','deltaloglh1','deltalogll1','deltaloglv1','logw','logprior','logpost','nulllogl','chain_log_evidence','chain_delta_log_evidence','chain_log_noise_evidence','chain_log_bayes_factor']
#Parameters known to cbcBPP
relativePhaseParams=[ a+b+'_relative_phase' for a,b in combinations(['h1','l1','v1'],2)]
snrParams=['snr','optimal_snr','matched_filter_snr','coherence'] + ['%s_optimal_snr'%(i) for i in ['h1','l1','v1']] + ['%s_cplx_snr_amp'%(i) for i in ['h1','l1','v1']] + ['%s_cplx_snr_arg'%(i) for i in ['h1', 'l1', 'v1']] + relativePhaseParams
calAmpParams=['calamp_%s'%(ifo) for ifo in ['h1','l1','v1']]
calPhaseParams=['calpha_%s'%(ifo) for ifo in ['h1','l1','v1']]
calParams = calAmpParams + calPhaseParams
# Masses
massParams=['m1','m2','chirpmass','mchirp','mc','eta','q','massratio','asym_massratio','mtotal','mf','mf_evol','mf_nonevol']
#Spins
spinParamsPrec=['a1','a2','phi1','theta1','phi2','theta2','costilt1','costilt2','costheta_jn','cosbeta','tilt1','tilt1_isco','tilt2','tilt2_isco','phi_jl','theta_jn','phi12','phi12_isco','af','af_evol','af_nonevol','afz','afz_evol','afz_nonevol']
spinParamsAli=['spin1','spin2','a1z','a2z']
spinParamsEff=['chi','effectivespin','chi_eff','chi_tot','chi_p']
spinParams=spinParamsPrec+spinParamsEff+spinParamsAli
# Source frame params
cosmoParam=['m1_source','m2_source','mtotal_source','mc_source','redshift','mf_source','mf_source_evol','mf_source_nonevol','m1_source_maxldist','m2_source_maxldist','mtotal_source_maxldist','mc_source_maxldist','redshift_maxldist','mf_source_maxldist','mf_source_maxldist_evol','mf_source_maxldist_nonevol']
#Strong Field
ppEParams=['ppEalpha','ppElowera','ppEupperA','ppEbeta','ppElowerb','ppEupperB','alphaPPE','aPPE','betaPPE','bPPE']
tigerParams=['dchi%i'%(i) for i in range(8)] + ['dchi%il'%(i) for i in [5,6] ] + ['dxi%d'%(i+1) for i in range(6)] + ['dalpha%i'%(i+1) for i in range(5)] + ['dbeta%i'%(i+1) for i in range(3)] + ['dsigma%i'%(i+1) for i in range(4)]
bransDickeParams=['omegaBD','ScalarCharge1','ScalarCharge2']
massiveGravitonParams=['lambdaG']
lorentzInvarianceViolationParams=['log10lambda_a','lambda_a','log10lambda_eff','lambda_eff','log10livamp','liv_amp']
tidalParams=['lambda1','lambda2','lam_tilde','dlam_tilde','lambdat','dlambdat','lambdas','bluni']
fourPiecePolyParams=['logp1','gamma1','gamma2','gamma3']
spectralParams=['sdgamma0','sdgamma1','sdgamma2','sdgamma3']
energyParams=['e_rad', 'e_rad_evol', 'e_rad_nonevol', 'l_peak', 'l_peak_evol', 'l_peak_nonevol', 'e_rad_maxldist', 'e_rad_maxldist_evol', 'e_rad_maxldist_nonevol']
strongFieldParams=ppEParams+tigerParams+bransDickeParams+massiveGravitonParams+tidalParams+fourPiecePolyParams+spectralParams+energyParams+lorentzInvarianceViolationParams

#Extrinsic
distParams=['distance','distMPC','dist','distance_maxl']
incParams=['iota','inclination','cosiota']
polParams=['psi','polarisation','polarization']
skyParams=['ra','rightascension','declination','dec']
phaseParams=['phase', 'phi0','phase_maxl']
#Times
timeParams=['time','time_mean']
endTimeParams=['l1_end_time','h1_end_time','v1_end_time']
#others
statsParams=['logprior','logl','deltalogl','deltaloglh1','deltalogll1','deltaloglv1','deltaloglh2','deltaloglg1']
calibParams=['calpha_l1','calpha_h1','calpha_v1','calamp_l1','calamp_h1','calamp_v1']

## Greedy bin sizes for cbcBPP and confidence leves used for the greedy bin intervals
confidenceLevels=[0.67,0.9,0.95,0.99]

greedyBinSizes={'mc':0.025,'m1':0.1,'m2':0.1,'mass1':0.1,'mass2':0.1,'mtotal':0.1,'mc_source':0.025,'m1_source':0.1,'m2_source':0.1,'mtotal_source':0.1,'mc_source_maxldist':0.025,'m1_source_maxldist':0.1,'m2_source_maxldist':0.1,'mtotal_source_maxldist':0.1,'eta':0.001,'q':0.01,'asym_massratio':0.01,'iota':0.01,'cosiota':0.02,'time':1e-4,'time_mean':1e-4,'distance':1.0,'dist':1.0,'distance_maxl':1.0,'redshift':0.01,'redshift_maxldist':0.01,'mchirp':0.025,'chirpmass':0.025,'spin1':0.04,'spin2':0.04,'a1z':0.04,'a2z':0.04,'a1':0.02,'a2':0.02,'phi1':0.05,'phi2':0.05,'theta1':0.05,'theta2':0.05,'ra':0.05,'dec':0.05,'chi':0.05,'chi_eff':0.05,'chi_tot':0.05,'chi_p':0.05,'costilt1':0.02,'costilt2':0.02,'thatas':0.05,'costheta_jn':0.02,'beta':0.05,'omega':0.05,'cosbeta':0.02,'ppealpha':1.0,'ppebeta':1.0,'ppelowera':0.01,'ppelowerb':0.01,'ppeuppera':0.01,'ppeupperb':0.01,'polarisation':0.04,'rightascension':0.05,'declination':0.05,'massratio':0.001,'inclination':0.01,'phase':0.05,'tilt1':0.05,'tilt2':0.05,'phi_jl':0.05,'theta_jn':0.05,'phi12':0.05,'flow':1.0,'phase_maxl':0.05,'calamp_l1':0.01,'calamp_h1':0.01,'calamp_v1':0.01,'calpha_h1':0.01,'calpha_l1':0.01,'calpha_v1':0.01,'logdistance':0.1,'psi':0.1,'costheta_jn':0.1,'mf':0.1,'mf_evol':0.1,'mf_source':0.1,'mf_source_evol':0.1,'mf_source_nonevol':0.1,'mf_source_maxldist':0.1,'mf_source_maxldist_evol':0.1,'mf_source_maxldist_nonevol':0.1,'af':0.02,'af_evol':0.02,'af_nonevol':0.02,'afz':0.02,'afz_evol':0.01,'afz_nonevol':0.01,'e_rad':0.1,'e_rad_evol':0.1,'e_rad_nonevol':0.1,'e_rad_maxldist':0.1,'e_rad_maxldist_evol':0.1,'e_rad_maxldist_nonevol':0.1,'l_peak':0.1,'l_peak_evol':0.1,'l_peak_nonevol':0.1}
for s in snrParams:
    greedyBinSizes[s]=0.05
for derived_time in ['h1_end_time','l1_end_time','v1_end_time','h1l1_delay','l1v1_delay','h1v1_delay']:
    greedyBinSizes[derived_time]=greedyBinSizes['time']
for derived_phase in relativePhaseParams:
    greedyBinSizes[derived_phase]=0.05
for param in tigerParams + bransDickeParams + massiveGravitonParams + lorentzInvarianceViolationParams:
    greedyBinSizes[param]=0.01
for param in tidalParams:
    greedyBinSizes[param]=2.5
for param in fourPiecePolyParams:
    greedyBinSizes[param]=2.5
for param in spectralParams:
    greedyBinSizes[param]=2.5
    #Confidence levels
for loglname in statsParams:
    greedyBinSizes[loglname]=0.1

#Pre-defined ordered list of line styles for use in matplotlib contour plots.
__default_line_styles=['solid', 'dashed', 'dashdot', 'dotted']
#Pre-defined ordered list of matplotlib colours for use in plots.
__default_color_lst=['r','b','y','g','c','m']
#A default css string for use in html results pages.
__default_css_string="""
p,h1,h2,h3,h4,h5
{
font-family:"Trebuchet MS", Arial, Helvetica, sans-serif;
}

p
{
font-size:14px;
}

h1
{
font-size:20px;
}

h2
{
font-size:18px;
}

h3
{
font-size:16px;
}



table
{
font-family:"Trebuchet MS", Arial, Helvetica, sans-serif;
width:100%;
border-collapse:collapse;
}
td,th
{
font-size:12px;
border:1px solid #B5C1CF;
padding:3px 7px 2px 7px;
}
th
{
font-size:14px;
text-align:left;
padding-top:5px;
padding-bottom:4px;
background-color:#B3CEEF;
color:#ffffff;
}
#postable tr:hover
{
background: #DFF4FF;
}
#covtable tr:hover
{
background: #DFF4FF;
}
#statstable tr:hover
{
background: #DFF4FF;
}

img {
    max-width: 510px;
    max-height: 510px;
    width:100%;
    eight:100%;
}

.ppsection
{
border-bottom-style:double;
}

"""
__default_javascript_string='''
//<![CDATA[
function toggle_visibility(tbid,lnkid)
{

  if(document.all){document.getElementById(tbid).style.display = document.getElementById(tbid).style.display == 'block' ? 'none' : 'block';}

  else{document.getElementById(tbid).style.display = document.getElementById(tbid).style.display == 'table' ? 'none' : 'table';}

  document.getElementById(lnkid).value = document.getElementById(lnkid).value == '[-] Collapse' ? '[+] Expand' : '[-] Collapse';

 }
 //]]>

'''


#===============================================================================
# Function to return the correct prior distribution for selected parameters
#===============================================================================
def get_prior(name):
    distributions={
      'm1':'uniform',
      'm2':'uniform',
      'mc':None,
      'eta':None,
      'q':None,
      'mtotal':'uniform',
      'm1_source':None,
      'm2_source':None,
      'mtotal_source':None,
      'mc_source':None,
      'redshift':None,
      'm1_source_maxldist':None,
      'm2_source_maxldist':None,
      'mtotal_source_maxldist':None,
      'mc_source_maxldist':None,
      'redshift_maxldist':None,
      'mf':None,
      'mf_evol':None,
      'mf_nonevol':None,
      'mf_source':None,
      'mf_source_evol':None,
      'mf_source_nonevol':None,
      'mf_source_maxldist':None,
      'mf_source_maxldist_evol':None,
      'mf_source_maxldist_nonevol':None,
      'af':None,
      'af_evol':None,
      'af_nonevol':None,
      'afz':None,
      'afz_evol':None,
      'afz_nonevol':None,
      'e_rad':None,
      'e_rad_evol':None,
      'e_rad_nonevol':None,
      'e_rad_maxldist':None,
      'e_rad_maxldist_evol':None,
      'e_rad_maxldist_nonevol':None,
      'l_peak':None,
      'l_peak_evol':None,
      'l_peak_nonevol':None,
      'spin1':'uniform',
      'spin2':'uniform',
      'a1':'uniform',
      'a2':'uniform',
      'a1z':'uniform',
      'a2z':'uniform',
      'theta1':'uniform',
      'theta2':'uniform',
      'phi1':'uniform',
      'phi2':'uniform',
      'chi_eff':None,
      'chi_tot':None,
      'chi_p':None,
      'tilt1':None,
      'tilt2':None,
      'tilt1_isco':None,
      'tilt2_isco':None,
      'costilt1':'uniform',
      'costilt2':'uniform',
      'iota':'np.cos',
      'cosiota':'uniform',
      'time':'uniform',
      'time_mean':'uniform',
      'dist':'x**2',
      'distance_maxl':'x**2',
      'ra':'uniform',
      'dec':'np.cos',
      'phase':'uniform',
      'psi':'uniform',
      'theta_jn':'np.sin',
      'costheta_jn':'uniform',
      'beta':None,
      'cosbeta':None,
      'phi_jl':None,
      'phi12':None,
      'phi12_isco':None,
      'logl':None,
      'h1_end_time':None,
      'l1_end_time':None,
      'v1_end_time':None,
      'h1l1_delay':None,
      'h1v1_delay':None,
      'l1v1_delay':None,
      'lambdat' :None,
      'dlambdat':None,
      'lambda1' : 'uniform',
      'lambda2': 'uniform',
      'lam_tilde' : None,
      'dlam_tilde': None,
      'lambdas':'uniform',
      'bluni':'uniform',
      'logp1':None,
      'gamma1':None,
      'gamma2':None,
      'gamma3':None,
      'sdgamma0': None,
      'sdgamma1': None,
      'sdgamma2': None,
      'sdgamma3': None,
      'calamp_h1' : 'uniform',
      'calamp_l1' : 'uniform',
      'calpha_h1' : 'uniform',
      'calpha_l1' : 'uniform',
      'polar_eccentricity':'uniform',
      'polar_angle':'uniform',
      'alpha':'uniform'
    }
    try:
        return distributions(name)
    except:
        return None

#===============================================================================
# Function used to generate plot labels.
#===============================================================================
def plot_label(param):
    """
    A lookup table for plot labels.
    """
    m1_names = ['mass1', 'm1']
    m2_names = ['mass2', 'm2']
    mc_names = ['mc','mchirp','chirpmass']
    eta_names = ['eta','massratio','sym_massratio']
    q_names = ['q','asym_massratio']
    iota_names = ['iota','incl','inclination']
    dist_names = ['dist','distance']
    ra_names = ['rightascension','ra']
    dec_names = ['declination','dec']
    phase_names = ['phi_orb', 'phi', 'phase', 'phi0']
    gr_test_names = ['dchi%d'%i for i in range(8)]+['dchil%d'%i for i in [5,6]]+['dxi%d'%(i+1) for i in range(6)]+['dalpha%d'%(i+1) for i in range(5)]+['dbeta%d'%(i+1) for i in range(3)]+['dsigma%d'%(i+1) for i in range(4)]

    labels={
        'm1':r'$m_1\,(\mathrm{M}_\odot)$',
        'm2':r'$m_2\,(\mathrm{M}_\odot)$',
        'mc':r'$\mathcal{M}\,(\mathrm{M}_\odot)$',
        'eta':r'$\eta$',
        'q':r'$q$',
        'mtotal':r'$M_\mathrm{total}\,(\mathrm{M}_\odot)$',
        'm1_source':r'$m_{1}^\mathrm{source}\,(\mathrm{M}_\odot)$',
        'm2_source':r'$m_{2}^\mathrm{source}\,(\mathrm{M}_\odot)$',
        'mtotal_source':r'$M_\mathrm{total}^\mathrm{source}\,(\mathrm{M}_\odot)$',
        'mc_source':r'$\mathcal{M}^\mathrm{source}\,(\mathrm{M}_\odot)$',
        'redshift':r'$z$',
        'm1_source_maxldist':r'$m_{1}^\mathrm{source - maxLdist}\,(\mathrm{M}_\odot)$',
        'm2_source_maxldist':r'$m_{2}^\mathrm{source - maxLdist}\,(\mathrm{M}_\odot)$',
        'mtotal_source_maxldist':r'$M_\mathrm{total}^\mathrm{source - maxLdist}\,(\mathrm{M}_\odot)$',
        'mc_source_maxldist':r'$\mathcal{M}^\mathrm{source - maxLdist}\,(\mathrm{M}_\odot)$',
        'redshift_maxldist':r'$z^\mathrm{maxLdist}$',
        'mf':r'$M_\mathrm{final}\,(\mathrm{M}_\odot)$',
        'mf_evol':r'$M_\mathrm{final}^\mathrm{evol}\,(\mathrm{M}_\odot)$',
        'mf_nonevol':r'$M_\mathrm{final}^\mathrm{non-evol}\,(\mathrm{M}_\odot)$',
        'mf_source':r'$M_\mathrm{final}^\mathrm{source}\,(\mathrm{M}_\odot)$',
        'mf_source_evol':r'$M_\mathrm{final}^\mathrm{source, evol}\,(\mathrm{M}_\odot)$',
        'mf_source_nonevol':r'$M_\mathrm{final}^\mathrm{source, non-evol}\,(\mathrm{M}_\odot)$',
        'mf_source_maxldist':r'$M_\mathrm{final}^\mathrm{source - maxLdist}\,(\mathrm{M}_\odot)$',
        'mf_source_maxldist_evol':r'$M_\mathrm{final}^\mathrm{source, evol - maxLdist}\,(\mathrm{M}_\odot)$',
        'mf_source_maxldist_nonevol':r'$M_\mathrm{final}^\mathrm{source, non-evol - maxLdist}\,(\mathrm{M}_\odot)$',
        'af':r'$a_\mathrm{final}$',
        'af_evol':r'$a_\mathrm{final}^\mathrm{evol}$',
        'af_nonevol':r'$a_\mathrm{final}^\mathrm{non-evol}$',
        'afz':r'$a_{\mathrm{final}, z}$',
        'afz_evol':r'$a_{\mathrm{final}, z}^\mathrm{evol}$',
        'afz_nonevol':r'$a_{\mathrm{final}, z}^\mathrm{non-evol}$',
        'e_rad':r'$E_\mathrm{rad}\,(\mathrm{M}_\odot)$',
        'e_rad_evol':r'$E_\mathrm{rad}^\mathrm{evol}\,(\mathrm{M}_\odot)$',
        'e_rad_nonevol':r'$E_\mathrm{rad}^\mathrm{non-evol}\,(\mathrm{M}_\odot)$',
        'e_rad_maxldist':r'$E_\mathrm{rad}^\mathrm{maxLdist}\,(\mathrm{M}_\odot)$',
        'e_rad_maxldist_evol':r'$E_\mathrm{rad}^\mathrm{evol - maxLdist}\,(\mathrm{M}_\odot)$',
        'e_rad_maxldist_nonevol':r'$E_\mathrm{rad}^\mathrm{non-evol - maxLdist}\,(\mathrm{M}_\odot)$',
        'l_peak':r'$L_\mathrm{peak}\,(10^{56}\,\mathrm{ergs}\,\mathrm{s}^{-1})$',
        'l_peak_evol':r'$L_\mathrm{peak}^\mathrm{evol}\,(10^{56}\,\mathrm{ergs}\,\mathrm{s}^{-1})$',
        'l_peak_nonevol':r'$L_\mathrm{peak}^\mathrm{non-evol}\,(10^{56}\,\mathrm{ergs}\,\mathrm{s}^{-1})$',
        'spin1':r'$S_1$',
        'spin2':r'$S_2$',
        'a1':r'$a_1$',
        'a2':r'$a_2$',
        'a1z':r'$a_{1z}$',
        'a2z':r'$a_{2z}$',
        'theta1':r'$\theta_1\,(\mathrm{rad})$',
        'theta2':r'$\theta_2\,(\mathrm{rad})$',
        'phi1':r'$\phi_1\,(\mathrm{rad})$',
        'phi2':r'$\phi_2\,(\mathrm{rad})$',
        'chi_eff':r'$\chi_\mathrm{eff}$',
        'chi_tot':r'$\chi_\mathrm{total}$',
        'chi_p':r'$\chi_\mathrm{P}$',
        'tilt1':r'$t_1\,(\mathrm{rad})$',
        'tilt2':r'$t_2\,(\mathrm{rad})$',
        'tilt1_isco':r'$t_1^\mathrm{ISCO}\,(\mathrm{rad})$',
        'tilt2_isco':r'$t_2^\mathrm{ISCO}\,(\mathrm{rad})$',
        'costilt1':r'$\mathrm{cos}(t_1)$',
        'costilt2':r'$\mathrm{cos}(t_2)$',
        'iota':r'$\iota\,(\mathrm{rad})$',
        'cosiota':r'$\mathrm{cos}(\iota)$',
        'time':r'$t_\mathrm{c}\,(\mathrm{s})$',
        'time_mean':r'$<t>\,(\mathrm{s})$',
        'dist':r'$d_\mathrm{L}\,(\mathrm{Mpc})$',
        'distance_maxl':r'$d_\mathrm{L}^\mathrm{maxL}\,(\mathrm{Mpc})$',
        'ra':r'$\alpha$',
        'dec':r'$\delta$',
        'phase':r'$\phi\,(\mathrm{rad})$',
        'phase_maxl':r'$\phi^\mathrm{maxL}\,(\mathrm{rad})$',
        'psi':r'$\psi\,(\mathrm{rad})$',
        'theta_jn':r'$\theta_\mathrm{JN}\,(\mathrm{rad})$',
        'costheta_jn':r'$\mathrm{cos}(\theta_\mathrm{JN})$',
        'beta':r'$\beta\,(\mathrm{rad})$',
        'cosbeta':r'$\mathrm{cos}(\beta)$',
        'phi_jl':r'$\phi_\mathrm{JL}\,(\mathrm{rad})$',
        'phi12':r'$\phi_\mathrm{12}\,(\mathrm{rad})$',
        'phi12_isco':r'$\phi_\mathrm{12}^\mathrm{ISCO}\,(\mathrm{rad})$',
        'logl':r'$\mathrm{log}(\mathcal{L})$',
        'h1_end_time':r'$t_\mathrm{H}$',
        'l1_end_time':r'$t_\mathrm{L}$',
        'v1_end_time':r'$t_\mathrm{V}$',
        'h1l1_delay':r'$\Delta t_\mathrm{HL}$',
        'h1v1_delay':r'$\Delta t_\mathrm{HV}$',
        'l1v1_delay':r'$\Delta t_\mathrm{LV}$',
        'lambdat' : r'$\tilde{\Lambda}$',
        'dlambdat': r'$\delta \tilde{\Lambda}$',
        'lambda1' : r'$\lambda_1$',
        'lambda2': r'$\lambda_2$',
        'lam_tilde' : r'$\tilde{\Lambda}$',
        'dlam_tilde': r'$\delta \tilde{\Lambda}$',
        'logp1':r'$\log(p_1)$',
        'gamma1':r'$\Gamma_1$',
        'gamma2':r'$\Gamma_2$',
        'gamma3':r'$\Gamma_3$',
        'sdgamma0' : r'$\gamma_{0}$',
        'sdgamma1' : r'$\gamma_{1}$',
        'sdgamma2' : r'$\gamma_{2}$',
        'sdgamma3' : r'$\gamma_{3}$',
        'calamp_h1' : r'$\delta A_{H1}$',
        'calamp_l1' : r'$\delta A_{L1}$',
        'calpha_h1' : r'$\delta \phi_{H1}$',
        'calpha_l1' : r'$\delta \phi_{L1}$',
        'polar_eccentricity':r'$\epsilon_{polar}$',
        'polar_angle':r'$\alpha_{polar}$',
        'alpha':r'$\alpha_{polar}$',
        'dchi0':r'$d\chi_0$',
        'dchi1':r'$d\chi_1$',
        'dchi2':r'$d\chi_2$',
        'dchi3':r'$d\chi_3$',
        'dchi4':r'$d\chi_4$',
        'dchi5':r'$d\chi_5$',
        'dchi5l':r'$d\chi_{5}^{(l)}$',
        'dchi6':r'$d\chi_6$',
        'dchi6l':r'$d\chi_{6}^{(l)}$',
        'dchi7':r'$d\chi_7$',
        'dxi1':r'$d\xi_1$',
        'dxi2':r'$d\xi_2$',
        'dxi3':r'$d\xi_3$',
        'dxi4':r'$d\xi_4$',
        'dxi5':r'$d\xi_5$',
        'dxi6':r'$d\xi_6$',
        'dalpha1':r'$d\alpha_1$',
        'dalpha2':r'$d\alpha_2$',
        'dalpha3':r'$d\alpha_3$',
        'dalpha4':r'$d\alpha_4$',
        'dalpha5':r'$d\alpha_5$',
        'dbeta1':r'$d\beta_1$',
        'dbeta2':r'$d\beta_2$',
        'dbeta3':r'$d\beta_3$',
        'dsigma1':r'$d\sigma_1$',
        'dsigma2':r'$d\sigma_2$',
        'dsigma3':r'$d\sigma_3$',
        'dsigma4':r'$d\sigma_4$',
        'optimal_snr':r'$\rho^{opt}$',
        'h1_optimal_snr':r'$\rho^{opt}_{H1}$',
        'l1_optimal_snr':r'$\rho^{opt}_{L1}$',
        'v1_optimal_snr':r'$\rho^{opt}_{V1}$',
        'matched_filter_snr':r'$\rho^{MF}$',
        'lambdas':r'$\Lambda_S$',
        'bluni' : r'$BL_{uniform}$',
        'log10lambda_a':r'$\log\lambda_{\mathbb{A}} [\mathrm{m}]$',
        'log10lambda_eff':r'$\log\lambda_{eff} [\mathrm{m}]$',
        'lambda_eff':r'$\lambda_{eff} [\mathrm{m}]$',
        'lambda_a':r'$\lambda_{\mathbb{A}} [\mathrm{m}]$',
        'liv_amp':r'$\mathbb{A} [\mathrm{{eV}^{2-\alpha}}]$' ,
        'log10livamp':r'$\log \mathbb{A}[\mathrm{{eV}^{2-\alpha}}]$'
      }

    # Handle cases where multiple names have been used
    if param in m1_names:
        param = 'm1'
    elif param in m2_names:
        param = 'm2'
    elif param in mc_names:
        param = 'mc'
    elif param in eta_names:
        param = 'eta'
    elif param in q_names:
        param = 'q'
    elif param in iota_names:
        param = 'iota'
    elif param in dist_names:
        param = 'dist'
    elif param in ra_names:
        param = 'ra'
    elif param in dec_names:
        param = 'dec'
    elif param in phase_names:
        param = 'phase'

    try:
        label = labels[param]
    except KeyError:
        # Use simple string if no formated label is available for param
        label = param

    return label

#===============================================================================
# Class definitions
#===============================================================================

class PosteriorOneDPDF(object):
    """
    A data structure representing one parameter in a chain of posterior samples.
    The Posterior class generates instances of this class for pivoting onto a given
    parameter (the Posterior class is per-Sampler oriented whereas this class represents
    the same one parameter in successive samples in the chain).
    """
    def __init__(self,name,posterior_samples,injected_value=None,injFref=None,trigger_values=None,prior=None):
        """
        Create an instance of PosteriorOneDPDF based on a table of posterior_samples.

        @param name: A literal string name for the parameter.
        @param posterior_samples: A 1D array of the samples.
        @param injected_value: The injected or real value of the parameter.
        @param injFref: reference frequency for injection
        @param trigger_values: The trigger values of the parameter (dictionary w/ IFOs as keys).
        @param prior: The prior value corresponding to each sample.
        """
        self.__name=name
        self.__posterior_samples=np.array(posterior_samples)

        self.__injFref=injFref
        self.__injval=injected_value
        self.__trigvals=trigger_values
        self.__prior=prior

        return

    def __len__(self):
        """
        Container method. Defined as number of samples.
        """
        return len(self.__posterior_samples)

    def __getitem__(self,idx):
        """
        Container method . Returns posterior containing sample idx (allows slicing).
        """
        return PosteriorOneDPDF(self.__name, self.__posterior_samples[idx], injected_value=self.__injval, f_ref=self.__f_ref, trigger_values=self.__trigvals)

    @property
    def name(self):
        """
        Return the string literal name of the parameter.

        """
        return self.__name

    @property
    def mean(self):
        """
        Return the arithmetic mean for the marginal PDF on the parameter.

        """
        return np.mean(self.__posterior_samples)

    @property
    def median(self):
        """
        Return the median value for the marginal PDF on the parameter.

        """
        return np.median(self.__posterior_samples)

    @property
    def stdev(self):
        """
        Return the standard deviation of the marginal PDF on the parameter.

        """
        try:
            stdev = sqrt(np.var(self.__posterior_samples))
            if not np.isfinite(stdev):
                raise OverflowError
        except OverflowError:
            mean = np.mean(self.__posterior_samples)
            stdev = mean * sqrt(np.var(self.__posterior_samples/mean))
        return stdev

    @property
    def stacc(self):
        """
        Return the 'standard accuracy statistic' (stacc) of the marginal
        posterior of the parameter.

        stacc is a standard deviant incorporating information about the
        accuracy of the waveform recovery. Defined as the mean of the sum
        of the squared differences between the points in the PDF
        (x_i - sampled according to the posterior) and the true value
        (\f$x_{true}\f$).  So for a marginalized one-dimensional PDF:
        \f$stacc = \sqrt{\frac{1}{N}\sum_{i=1}^N (x_i-x_{\rm true})2}\f$

        """
        if self.__injval is None:
            return None
        else:
            return np.sqrt(np.mean((self.__posterior_samples - self.__injval)**2.0))

    @property
    def injval(self):
        """
        Return the injected value set at construction . If no value was set
        will return None .

        """
        return self.__injval

    @property
    def trigvals(self):
        """
        Return the trigger values set at construction. If no value was set
        will return None .

        """
        return self.__trigvals

    #@injval.setter #Python 2.6+
    def set_injval(self,new_injval):
        """
        Set the injected/real value of the parameter.

        @param new_injval: The injected/real value to set.
        """

        self.__injval=new_injval

    def set_trigvals(self,new_trigvals):
        """
        Set the trigger values of the parameter.

        @param new_trigvals: Dictionary containing trigger values with IFO keys.
        """

        self.__trigvals=new_trigvals

    @property
    def samples(self):
        """
        Return a 1D numpy.array of the samples.

        """
        return self.__posterior_samples

    def delete_samples_by_idx(self,samples):
        """
        Remove samples from posterior, analagous to numpy.delete but opperates in place.

        @param samples: A list containing the indexes of the samples to be removed.
        """
        self.__posterior_samples=np.delete(self.__posterior_samples,samples).reshape(-1,1)

    @property
    def gaussian_kde(self):
        """
        Return a SciPy gaussian_kde (representing a Gaussian KDE) of the samples.

        """
        from numpy import seterr as np_seterr
        from scipy import seterr as sp_seterr

        np_seterr(under='ignore')
        sp_seterr(under='ignore')
        try:
            return_value=stats.kde.gaussian_kde(np.transpose(self.__posterior_samples))
        except:
            exfile=open('exception.out','w')
            np.savetxt(exfile,self.__posterior_samples)
            exfile.close()
            raise

        return return_value

    @property
    def KL(self):
        """Returns the KL divergence between the prior and the posterior.
        It measures the relative information content in nats. The prior is evaluated
        at run time. It defaults to None. If None is passed, it just returns the information content
        of the posterior."
        """

        def uniform(x):
            return np.array([1./(np.max(x)-np.min(x)) for _ in x])

        posterior, dx = np.histogram(self.samples,bins=36,density=True)
        from scipy.stats import entropy
        # check the kind of prior and process the string accordingly
        prior = get_prior(self.name)
        if prior is None:
            raise ValueError
        elif prior=='uniform':
            prior+='(self.samples)'
        elif 'x' in prior:
            prior.replace('x','self.samples')
        elif not(prior.startswith('np.')):
            prior = 'np.'+prior
            prior+='(self.samples)'
        else:
            raise ValueError

        try:
            prior_dist = eval(prior)
        except:
            raise ValueError

        return entropy(posterior, qk=prior_dist)

    def prob_interval(self,intervals):
        """
        Evaluate probability intervals.

        @param intervals: A list of the probability intervals [0-1]
        """
        list_of_ci=[]
        samples_temp=np.sort(np.squeeze(self.samples))

        for interval in intervals:
            if interval<1.0:
                samples_temp
                N=np.size(samples_temp)
                #Find index of lower bound
                lower_idx=int(floor((N/2.0)*(1-interval)))
                if lower_idx<0:
                    lower_idx=0
                #Find index of upper bound
                upper_idx=N-int(floor((N/2.0)*(1-interval)))
                if upper_idx>N:
                    upper_idx=N-1

                list_of_ci.append((float(samples_temp[lower_idx]),float(samples_temp[upper_idx])))
            else:
                list_of_ci.append((None,None))

        return list_of_ci

class Posterior(object):
    """
    Data structure for a table of posterior samples .
    """
    def __init__(self,commonResultsFormatData,SimInspiralTableEntry=None,inj_spin_frame='OrbitalL', injFref=100,SnglInspiralList=None,name=None,description=None):
        """
        Constructor.

        @param commonResultsFormatData: A 2D array containing the posterior
            samples and related data. The samples chains form the columns.
        @param SimInspiralTableEntry: A SimInspiralTable row containing the injected values.
        @param SnglInspiralList: A list of SnglInspiral objects containing the triggers.
        @param inj_spin_frame: spin frame
        @param injFref: reference frequency
        @param name: optional name
        @param description: optional description

        """
        common_output_table_header,common_output_table_raw =commonResultsFormatData
        self._posterior={}
        self._injFref=injFref
        self._injection=SimInspiralTableEntry

        self._triggers=SnglInspiralList
        self._loglaliases=['deltalogl', 'posterior', 'logl','logL','likelihood']
        self._logpaliases=['logp', 'logP','prior','logprior','Prior','logPrior']

        common_output_table_header=[i.lower() for i in common_output_table_header]

        # Define XML mapping
        self._injXMLFuncMap={
                            'mchirp':lambda inj:inj.mchirp,
                            'chirpmass':lambda inj:inj.mchirp,
                            'mc':lambda inj:inj.mchirp,
                            'mass1':lambda inj:inj.mass1,
                            'm1':lambda inj:inj.mass1,
                            'mass2':lambda inj:inj.mass2,
                            'm2':lambda inj:inj.mass2,
                            'mtotal':lambda inj:float(inj.mass1)+float(inj.mass2),
                            'eta':lambda inj:inj.eta,
                            'q':self._inj_q,
                            'asym_massratio':self._inj_q,
                            'massratio':lambda inj:inj.eta,
                            'sym_massratio':lambda inj:inj.eta,
                            'time': lambda inj:float(inj.get_end()),
                            'time_mean': lambda inj:float(inj.get_end()),
                            'end_time': lambda inj:float(inj.get_end()),
                            'phi0':lambda inj:inj.phi0,
                            'phi_orb': lambda inj: inj.coa_phase,
                            'phase': lambda inj: inj.coa_phase,
                            'dist':lambda inj:inj.distance,
                            'distance':lambda inj:inj.distance,
                            'ra':self._inj_longitude,
                            'rightascension':self._inj_longitude,
                            'long':self._inj_longitude,
                            'longitude':self._inj_longitude,
                            'dec':lambda inj:inj.latitude,
                            'declination':lambda inj:inj.latitude,
                            'lat':lambda inj:inj.latitude,
                            'latitude':lambda inj:inj.latitude,
                            'psi': lambda inj: np.mod(inj.polarization, np.pi),
                            'f_ref': lambda inj: self._injFref,
                            'polarisation':lambda inj:inj.polarization,
                            'polarization':lambda inj:inj.polarization,
                            'h1_end_time':lambda inj:float(inj.get_end('H')),
                            'l1_end_time':lambda inj:float(inj.get_end('L')),
                            'v1_end_time':lambda inj:float(inj.get_end('V')),
                            'lal_amporder':lambda inj:inj.amp_order}

        # Add on all spin parameterizations
        for key, val in self._inj_spins(self._injection, frame=inj_spin_frame).items():
            self._injXMLFuncMap[key] = val

        for one_d_posterior_samples,param_name in zip(np.hsplit(common_output_table_raw,common_output_table_raw.shape[1]),common_output_table_header):

            self._posterior[param_name]=PosteriorOneDPDF(param_name.lower(),one_d_posterior_samples,injected_value=self._getinjpar(param_name),injFref=self._injFref,trigger_values=self._gettrigpar(param_name))

        if 'mchirp' in common_output_table_header and 'eta' in common_output_table_header \
        and (not 'm1' in common_output_table_header) and (not 'm2' in common_output_table_header):
            try:
                print('Inferring m1 and m2 from mchirp and eta')
                (m1,m2)=mc2ms(self._posterior['mchirp'].samples, self._posterior['eta'].samples)
                self._posterior['m1']=PosteriorOneDPDF('m1',m1,injected_value=self._getinjpar('m1'),trigger_values=self._gettrigpar('m1'))
                self._posterior['m2']=PosteriorOneDPDF('m2',m2,injected_value=self._getinjpar('m2'),trigger_values=self._gettrigpar('m2'))
            except KeyError:
                print('Unable to deduce m1 and m2 from input columns')


        logLFound=False

        for loglalias in self._loglaliases:

            if loglalias in common_output_table_header:
                try:
                    self._logL=self._posterior[loglalias].samples
                except KeyError:
                    print("No '%s' column in input table!"%loglalias)
                    continue
                logLFound=True

        if not logLFound:
            raise RuntimeError("No likelihood/posterior values found!")
        self._logP=None

        for logpalias in self._logpaliases:
            if logpalias in common_output_table_header:
                try:
                    self._logP=self._posterior[logpalias].samples
                except KeyError:
                    print("No '%s' column in input table!"%logpalias)
                    continue
                if not 'log' in logpalias:
                    self._logP=[np.log(i) for i in self._logP]

        if name is not None:
            self.__name=name

        if description is not None:
            self.__description=description

        return

    def extend_posterior(self):
        """
        Add some useful derived parameters (such as tilt angles, time delays, etc) in the Posterior object
        """
        injection=self._injection
        pos=self
        # Generate component mass posterior samples (if they didnt exist already)
        if 'mc' in pos.names:
            mchirp_name = 'mc'
        elif 'chirpmass' in pos.names:
            mchirp_name = 'chirpmass'
        else:
            mchirp_name = 'mchirp'

        if 'asym_massratio' in pos.names:
            q_name = 'asym_massratio'
        else:
            q_name = 'q'

        if 'sym_massratio' in pos.names:
            eta_name= 'sym_massratio'
        elif 'massratio' in pos.names:
            eta_name= 'massratio'
        else:
            eta_name='eta'

        if 'mass1' in pos.names and 'mass2' in pos.names:
            pos.append_mapping(('m1','m2'), lambda x,y:(x,y), ('mass1','mass2'))

        if (mchirp_name in pos.names and eta_name in pos.names) and \
        ('mass1' not in pos.names or 'm1' not in pos.names) and \
        ('mass2' not in pos.names or 'm2' not in pos.names):

            pos.append_mapping(('m1','m2'),mc2ms,(mchirp_name,eta_name))

        if (mchirp_name in pos.names and q_name in pos.names) and \
        ('mass1' not in pos.names or 'm1' not in pos.names) and \
        ('mass2' not in pos.names or 'm2' not in pos.names):

            pos.append_mapping(('m1','m2'),q2ms,(mchirp_name,q_name))
            pos.append_mapping('eta',q2eta,q_name)

        if ('m1' in pos.names and 'm2' in pos.names and not 'mtotal' in pos.names ):
            pos.append_mapping('mtotal', lambda m1,m2: m1+m2, ('m1','m2') )

        if('a_spin1' in pos.names): pos.append_mapping('a1',lambda a:a,'a_spin1')
        if('a_spin2' in pos.names): pos.append_mapping('a2',lambda a:a,'a_spin2')
        if('phi_spin1' in pos.names): pos.append_mapping('phi1',lambda a:a,'phi_spin1')
        if('phi_spin2' in pos.names): pos.append_mapping('phi2',lambda a:a,'phi_spin2')
        if('theta_spin1' in pos.names): pos.append_mapping('theta1',lambda a:a,'theta_spin1')
        if('theta_spin2' in pos.names): pos.append_mapping('theta2',lambda a:a,'theta_spin2')

        my_ifos=['h1','l1','v1']
        for ifo1,ifo2 in combinations(my_ifos,2):
            p1=ifo1+'_cplx_snr_arg'
            p2=ifo2+'_cplx_snr_arg'
            if p1 in pos.names and p2 in pos.names:
                delta=np.mod(pos[p1].samples - pos[p2].samples + np.pi ,2.0*np.pi)-np.pi
                pos.append(PosteriorOneDPDF(ifo1+ifo2+'_relative_phase',delta))

        # Ensure that both theta_jn and inclination are output for runs
        # with zero tilt (for runs with tilt, this will be taken care of
        # below when the old spin angles are computed as functions of the
        # new ones
        # Disabled this since the parameters are degenerate and causing problems
        #if ('theta_jn' in pos.names) and (not 'tilt1' in pos.names) and (not 'tilt2' in pos.names):
        #    pos.append_mapping('iota', lambda t:t, 'theta_jn')

        # Compute time delays from sky position
        try:
            if ('ra' in pos.names or 'rightascension' in pos.names) \
            and ('declination' in pos.names or 'dec' in pos.names) \
            and 'time' in pos.names:
                from lal import LIGOTimeGPS, TimeDelayFromEarthCenter
                from numpy import array
                detMap = {'H1': 'LHO_4k', 'H2': 'LHO_2k', 'L1': 'LLO_4k',
                        'G1': 'GEO_600', 'V1': 'VIRGO', 'T1': 'TAMA_300'}
                if 'ra' in pos.names:
                    ra_name='ra'
                else: ra_name='rightascension'
                if 'dec' in pos.names:
                    dec_name='dec'
                else: dec_name='declination'
                ifo_times={}
                my_ifos=['H1','L1','V1']
                for ifo in my_ifos:
                    inj_time=None
                    if injection:
                        inj_time=float(injection.get_end(ifo[0]))
                    location = lal.cached_detector_by_prefix[ifo].location
                    ifo_times[ifo]=array(list(map(lambda ra,dec,time: array([time[0]+TimeDelayFromEarthCenter(location,ra[0],dec[0],LIGOTimeGPS(float(time[0])))]), pos[ra_name].samples,pos[dec_name].samples,pos['time'].samples)))
                    loc_end_time=PosteriorOneDPDF(ifo.lower()+'_end_time',ifo_times[ifo],injected_value=inj_time)
                    pos.append(loc_end_time)
                for ifo1 in my_ifos:
                    for ifo2 in my_ifos:
                        if ifo1==ifo2: continue
                        delay_time=ifo_times[ifo2]-ifo_times[ifo1]
                        if injection:
                            inj_delay=float(injection.get_end(ifo2[0])-injection.get_end(ifo1[0]))
                        else:
                            inj_delay=None
                        time_delay=PosteriorOneDPDF(ifo1.lower()+ifo2.lower()+'_delay',delay_time,inj_delay)
                        pos.append(time_delay)
        except ImportError:
            print('Warning: Could not import lal python bindings, check you ./configured with --enable-swig-python')
            print('This means I cannot calculate time delays')

        #Calculate new spin angles
        new_spin_params = ['tilt1','tilt2','theta_jn','beta']
        if not set(new_spin_params).issubset(set(pos.names)):
            old_params = ['f_ref',mchirp_name,'eta','iota','a1','theta1','phi1']
            if 'a2' in pos.names: old_params += ['a2','theta2','phi2']
            try:
                pos.append_mapping(new_spin_params, spin_angles, old_params)
            except KeyError:
                print("Warning: Cannot find spin parameters.  Skipping spin angle calculations.")

        #Store signed spin magnitudes in separate parameters and make a1,a2 magnitudes
        if 'a1' in pos.names:
            if 'tilt1' in pos.names:
                pos.append_mapping('a1z', lambda a, tilt: a*np.cos(tilt), ('a1','tilt1'))
            else:
                pos.append_mapping('a1z', lambda x: x, 'a1')
                inj_az = None
                if injection is not None:
                    inj_az = injection.spin1z
                pos['a1z'].set_injval(inj_az)
                pos.pop('a1')
                pos.append_mapping('a1', lambda x: np.abs(x), 'a1z')

        if 'a2' in pos.names:
            if 'tilt2' in pos.names:
                pos.append_mapping('a2z', lambda a, tilt: a*np.cos(tilt), ('a2','tilt2'))
            else:
                pos.append_mapping('a2z', lambda x: x, 'a2')
                inj_az = None
                if injection is not None:
                    inj_az = injection.spin2z
                pos['a2z'].set_injval(inj_az)
                pos.pop('a2')
                pos.append_mapping('a2', lambda x: np.abs(x), 'a2z')

        #Calculate effective spin parallel to L
        if ('m1' in pos.names and 'a1z' in pos.names) and ('m2' in pos.names and 'a2z' in pos.names):
            pos.append_mapping('chi_eff', lambda m1,a1z,m2,a2z: (m1*a1z + m2*a2z) / (m1 + m2), ('m1','a1z','m2','a2z'))

        #If precessing spins calculate total effective spin
        if ('m1' in pos.names and 'a1' in pos.names and 'tilt1' in pos.names) and ('m2' in pos.names and 'a2' in pos.names and 'tilt2' in pos.names):
            pos.append_mapping('chi_tot', lambda m1,a1,m2,a2: (m1*a1 + m2*a2) / (m1 + m2), ('m1','a1','m2','a2'))

        #Calculate effective precessing spin magnitude
        if ('m1' in pos.names and 'a1' in pos.names and 'tilt1' in pos.names) and ('m2' in pos.names and 'a2' in pos.names and 'tilt2' in pos.names):
            pos.append_mapping('chi_p', chi_precessing, ['m1', 'a1', 'tilt1', 'm2', 'a2', 'tilt2'])

        # Calculate redshift from luminosity distance measurements
        if('distance' in pos.names):
            pos.append_mapping('redshift', calculate_redshift, 'distance')
        elif('dist' in pos.names):
            pos.append_mapping('redshift', calculate_redshift, 'dist')
        # If using the DistanceMarginalisation, compute the maxL redshift distribution from the maxL d_L
        elif('distance_maxl' in pos.names):
            pos.append_mapping('redshift_maxldist', calculate_redshift, 'distance_maxl')

        # Calculate source mass parameters
        if ('m1' in pos.names) and ('redshift' in pos.names):
            pos.append_mapping('m1_source', source_mass, ['m1', 'redshift'])

        if ('m2' in pos.names) and ('redshift' in pos.names):
            pos.append_mapping('m2_source', source_mass, ['m2', 'redshift'])

        if ('mtotal' in pos.names) and ('redshift' in pos.names):
            pos.append_mapping('mtotal_source', source_mass, ['mtotal', 'redshift'])

        if ('mc' in pos.names) and ('redshift' in pos.names):
            pos.append_mapping('mc_source', source_mass, ['mc', 'redshift'])

        # Calculate source mass parameters if DistanceMarginalisation was used, using the maxL distance and redshift
        if ('m1' in pos.names) and ('redshift_maxldist' in pos.names):
            pos.append_mapping('m1_source_maxldist', source_mass, ['m1', 'redshift_maxldist'])

        if ('m2' in pos.names) and ('redshift_maxldist' in pos.names):
            pos.append_mapping('m2_source_maxldist', source_mass, ['m2', 'redshift_maxldist'])

        if ('mtotal' in pos.names) and ('redshift_maxldist' in pos.names):
            pos.append_mapping('mtotal_source_maxldist', source_mass, ['mtotal', 'redshift_maxldist'])

        if ('mc' in pos.names) and ('redshift_maxldist' in pos.names):
            pos.append_mapping('mc_source_maxldist', source_mass, ['mc', 'redshift_maxldist'])

        # Calling functions testing Lorentz invariance violation
        if ('log10lambda_eff' in pos.names) and ('redshift' in pos.names):
            pos.append_mapping('log10lambda_a', lambda z,nonGR_alpha,wl,dist:np.log10(lambda_a(z, nonGR_alpha, 10**wl, dist)), ['redshift', 'nonGR_alpha', 'log10lambda_eff', 'dist'])
        if ('log10lambda_eff' in pos.names) and ('redshift' in pos.names):
            pos.append_mapping('log10livamp', lambda z,nonGR_alpha,wl,dist:np.log10(amplitudeMeasure(z, nonGR_alpha, 10**wl, dist)), ['redshift','nonGR_alpha','log10lambda_eff', 'dist'])
        if ('lambda_eff' in pos.names) and ('redshift' in pos.names):
            pos.append_mapping('lambda_a', lambda_a, ['redshift', 'nonGR_alpha', 'log10lambda_eff', 'dist'])
        if ('lambda_eff' in pos.names) and ('redshift' in pos.names):
            pos.append_mapping('liv_amp', amplitudeMeasure, ['redshift','nonGR_alpha','lambda_eff', 'dist'])

        #Calculate new tidal parameters
        new_tidal_params = ['lam_tilde','dlam_tilde']
        old_tidal_params = ['lambda1','lambda2','q']
        if 'lambda1' in pos.names or 'lambda2' in pos.names:
            try:
                pos.append_mapping(new_tidal_params, symm_tidal_params, old_tidal_params)
            except KeyError:
                print("Warning: Cannot find tidal parameters.  Skipping tidal calculations.")

        #If new spin params present, calculate old ones
        old_spin_params = ['iota', 'theta1', 'phi1', 'theta2', 'phi2', 'beta']
        new_spin_params = ['theta_jn', 'phi_jl', 'tilt1', 'tilt2', 'phi12', 'a1', 'a2', 'm1', 'm2', 'f_ref','phase']
        try:
            if pos['f_ref'].samples[0][0]==0.0:
                for name in ['flow','f_lower']:
                    if name in pos.names:
                        new_spin_params = ['theta_jn', 'phi_jl', 'tilt1', 'tilt2', 'phi12', 'a1', 'a2', 'm1', 'm2', name]
        except:
            print("No f_ref for SimInspiralTransformPrecessingNewInitialConditions().")
        if set(new_spin_params).issubset(set(pos.names)) and not set(old_spin_params).issubset(set(pos.names)):
            pos.append_mapping(old_spin_params, physical2radiationFrame, new_spin_params)

        #Calculate spin magnitudes for aligned runs
        if 'spin1' in pos.names:
            inj_a1 = inj_a2 = None
            if injection:
                inj_a1 = sqrt(injection.spin1x*injection.spin1x + injection.spin1y*injection.spin1y + injection.spin1z*injection.spin1z)
                inj_a2 = sqrt(injection.spin2x*injection.spin2x + injection.spin2y*injection.spin2y + injection.spin2z*injection.spin2z)

            try:
                a1_samps = abs(pos['spin1'].samples)
                a1_pos = PosteriorOneDPDF('a1',a1_samps,injected_value=inj_a1)
                pos.append(a1_pos)
            except KeyError:
                print("Warning: problem accessing spin1 values.")

            try:
                a2_samps = abs(pos['spin2'].samples)
                a2_pos = PosteriorOneDPDF('a2',a2_samps,injected_value=inj_a2)
                pos.append(a2_pos)
            except KeyError:
                print("Warning: no spin2 values found.")

        # For BBHs: Calculate mass and spin of final merged system, radiated energy, and peak luminosity in GWs

        # Only apply fits if this is a BBH run (with no tidal parameters)

        if len(np.intersect1d(pos.names,tidalParams)) == 0:

              # Set fits to consider (and average over)

              FinalSpinFits = ['HBR2016', 'UIB2016', 'HL2016']
              FinalMassFits = ['UIB2016', 'HL2016']
              LpeakFits     = ['UIB2016', 'HL2016']

              # If evolved spin angle samples are present, use those to compute the final mass and spin, peak luminosity, and radiated energy; also, use the _evol suffix in the aligned-spin case, since here the spin angles are trivially evolved

              spin_angle_suffix = ''
              evol_suffix = '_evol'

              if all([x in pos.names for x in ['tilt1_isco','tilt2_isco','phi12_isco']]):
                  spin_angle_suffix = '_isco'
              elif all([x in pos.names for x in ['tilt1','tilt2','phi12']]):
                  evol_suffix = '_nonevol'

              zero_vec = np.array([0.])

              tilt1_name = 'tilt1' + spin_angle_suffix
              tilt2_name = 'tilt2' + spin_angle_suffix
              phi12_name = 'phi12' + spin_angle_suffix
              mf_name = 'mf' + evol_suffix
              mf_source_name = 'mf_source' + evol_suffix
              mf_source_maxldist_name = 'mf_source_maxldist' + evol_suffix

              if ('m1' in pos.names) and ('m2' in pos.names):
                  if ('a1' in pos.names) and ('a2' in pos.names):
                      if (tilt1_name in pos.names) and (tilt2_name in pos.names) and (phi12_name in pos.names):
                          # Precessing case
                          print("Using averages of fit formulae for final mass, final spin, and peak luminosity (on masses and 3D spins).")
                          if evol_suffix == '_evol':
                              print("Applying these to *_isco evolved spin samples and outputting *_evol samples.")
                          else:
                              print("Applying these to unevolved spin samples and outputting *_nonevol samples.")
                          print("Final mass fits:", FinalMassFits, "; Final spin fits:", FinalSpinFits, "; Peak luminosity fits:", LpeakFits)
                          try:
                              pos.append_mapping('af' + evol_suffix, lambda m1, m2, chi1, chi2, tilt1, tilt2, phi12: bbh_average_fits_precessing(m1, m2, chi1, chi2, tilt1, tilt2, phi12, 'af', FinalSpinFits), ['m1', 'm2', 'a1', 'a2', tilt1_name, tilt2_name, phi12_name])
                          except Exception as e:
                              print("Could not calculate %s. The error was: %s"%('af' + evol_suffix, str(e)))
                          try:
                              pos.append_mapping('afz' + evol_suffix, lambda m1, m2, chi1, chi2, tilt1, tilt2: bbh_average_fits_precessing(m1, m2, chi1, chi2, tilt1, tilt2, zero_vec, 'afz', FinalSpinFits), ['m1', 'm2', 'a1', 'a2', tilt1_name, tilt2_name])
                          except Exception as e:
                              print("Could not calculate %s. The error was: %s"%('afz' + evol_suffix, str(e)))
                          try:
                              pos.append_mapping(mf_name, lambda m1, m2, chi1, chi2, tilt1, tilt2: bbh_average_fits_precessing(m1, m2, chi1, chi2, tilt1, tilt2, zero_vec, 'Mf', FinalMassFits), ['m1', 'm2', 'a1', 'a2', tilt1_name, tilt2_name])
                          except Exception as e:
                              print("Could not calculate %s. The error was: %s"%(mf_name, str(e)))
                          try:
                              pos.append_mapping('l_peak' + evol_suffix, lambda m1, m2, chi1, chi2, tilt1, tilt2: bbh_average_fits_precessing(m1, m2, chi1, chi2, tilt1, tilt2, zero_vec, 'Lpeak', LpeakFits), ['m1', 'm2', 'a1', 'a2', tilt1_name, tilt2_name])
                          except Exception as e:
                              print("Could not calculate %s. The error was: %s"%('l_peak' + evol_suffix, str(e)))
                      elif ('a1z' in pos.names) and ('a2z' in pos.names):
                          # Aligned-spin case
                          print("Using averages for final mass, final spin, and peak luminosity (on masses and projected spin components).")
                          print("Outputting *_evol samples because spin evolution is trivial in this nonprecessing case.")
                          print("Final mass fits:", FinalMassFits, "; Final spin fits:", FinalSpinFits, "; Peak luminosity fits:", LpeakFits)
                          try:
                              # Compute absolute values of spins and compute tilt angles to allow for negative spin values
                              pos.append_mapping('afz_evol', lambda m1, m2, chi1, chi2: bbh_average_fits_precessing(m1, m2, abs(chi1), abs(chi2), 0.5*np.pi*(1. - np.sign(chi1)), 0.5*np.pi*(1. - np.sign(chi2)), zero_vec, 'afz', FinalSpinFits), ['m1', 'm2', 'a1z', 'a2z'])
                          except Exception as e:
                              print("Could not calculate afz_evol. The error was: %s"%(str(e)))
                          try:
                              pos.append_mapping('af_evol', lambda a: abs(a), 'afz_evol')
                          except Exception as e:
                              print("Could not calculate af_evol. The error was: %s"%(str(e)))
                          try:
                              pos.append_mapping('mf_evol', lambda m1, m2, chi1, chi2: bbh_average_fits_precessing(m1, m2, abs(chi1), abs(chi2), 0.5*np.pi*(1. - np.sign(chi1)), 0.5*np.pi*(1. - np.sign(chi2)), zero_vec, 'Mf', FinalMassFits), ['m1', 'm2', 'a1z', 'a2z'])
                          except Exception as e:
                              print("Could not calculate mf_evol. The error was: %s"%(str(e)))
                          try:
                              pos.append_mapping('l_peak_evol', lambda m1, m2, chi1, chi2: bbh_average_fits_precessing(m1, m2, abs(chi1), abs(chi2), 0.5*np.pi*(1. - np.sign(chi1)), 0.5*np.pi*(1. - np.sign(chi2)), zero_vec, 'Lpeak', LpeakFits), ['m1', 'm2', 'a1z', 'a2z'])
                          except Exception as e:
                              print("Could not calculate l_peak_evol. The error was: %s"%(str(e)))
                      else:
                          print("Could not calculate final parameters or Lpeak. Found samples for a1 and a2 but not for tilt angles and phi12 or spin components (a1z and a2z).")
                  else:
                      # Nonspinning case
                      print("Using averages of fit formulae for final mass, final spin, and peak luminosity (on masses and zero spins).")
                      print("Outputting *_evol samples because spin evolution is trivial in this nonspinning case.")
                      print("Final mass fits:", FinalMassFits, "; Final spin fits:", FinalSpinFits, "; Peak luminosity fits:", LpeakFits)
                      try:
                          pos.append_mapping('afz_evol', lambda m1, m2: bbh_average_fits_precessing(m1, m2, zero_vec, zero_vec, zero_vec, zero_vec, zero_vec, 'afz', FinalSpinFits), ['m1', 'm2'])
                      except Exception as e:
                              print("Could not calculate afz_evol. The error was: %s"%(str(e)))
                      try:
                          pos.append_mapping('af_evol', lambda a: abs(a), 'afz_evol')
                      except Exception as e:
                              print("Could not calculate af_evol. The error was: %s"%(str(e)))
                      try:
                          pos.append_mapping('mf_evol', lambda m1, m2: bbh_average_fits_precessing(m1, m2, zero_vec, zero_vec, zero_vec, zero_vec, zero_vec, 'Mf', FinalMassFits), ['m1', 'm2'])
                      except Exception as e:
                              print("Could not calculate mf_evol. The error was: %s"%(str(e)))
                      try:
                          pos.append_mapping('l_peak_evol', lambda m1, m2: bbh_average_fits_precessing(m1, m2, zero_vec, zero_vec, zero_vec, zero_vec, zero_vec, 'Lpeak', LpeakFits), ['m1', 'm2'])
                      except Exception as e:
                          print("Could not calculate l_peak_evol. The error was: %s"%(str(e)))

              # Convert final mass to source frame
              if (mf_name in pos.names) and ('redshift' in pos.names):
                  try:
                      pos.append_mapping(mf_source_name, source_mass, [mf_name, 'redshift'])
                  except Exception as e:
                      print("Could not calculate final source frame mass. The error was: %s"%(str(e)))

              if (mf_name in pos.names) and ('redshift_maxldist' in pos.names):
                  try:
                      pos.append_mapping(mf_source_maxldist_name, source_mass, [mf_name, 'redshift_maxldist'])
                  except Exception as e:
                      print("Could not calculate final source frame mass using maxldist redshift. The error was: %s"%(str(e)))

              # Calculate radiated energy
              if ('mtotal_source' in pos.names) and (mf_source_name in pos.names):
                  try:
                      pos.append_mapping('e_rad' + evol_suffix, lambda mtot_s, mf_s: mtot_s-mf_s, ['mtotal_source', mf_source_name])
                  except Exception as e:
                      print("Could not calculate radiated energy. The error was: %s"%(str(e)))

              if ('mtotal_source_maxldist' in pos.names) and (mf_source_maxldist_name in pos.names):
                  try:
                      pos.append_mapping('e_rad_maxldist' + evol_suffix, lambda mtot_s, mf_s: mtot_s-mf_s, ['mtotal_source_maxldist', mf_source_maxldist_name])
                  except Exception as e:
                      print("Could not calculate radiated energy using maxldist redshift results. The error was: %s"%(str(e)))

    def bootstrap(self):
        """
        Returns a new Posterior object that contains a bootstrap
        sample of self.

        """
        names=[]
        samples=[]
        for name,oneDpos in self._posterior.items():
            names.append(name)
            samples.append(oneDpos.samples)

        samplesBlock=np.hstack(samples)

        bootstrapSamples=samplesBlock[:,:]
        Nsamp=bootstrapSamples.shape[0]

        rows=np.vsplit(samplesBlock,Nsamp)

        for i in range(Nsamp):
            bootstrapSamples[i,:]=random.choice(rows)

        return Posterior((names,bootstrapSamples),self._injection,self._triggers)

    def delete_samples_by_idx(self,samples):
        """
        Remove samples from all OneDPosteriors.

        @param samples: The indexes of the samples to be removed.
        """
        for name,pos in self:
            pos.delete_samples_by_idx(samples)
        return

    def delete_NaN_entries(self,param_list):
        """
        Remove samples containing NaN in request params.

        @param param_list: The parameters to be checked for NaNs.
        """
        nan_idxs = np.array(())
        nan_dict = {}
        for param in param_list:
            nan_bool_array = np.isnan(self[param].samples).any(1)
            idxs = np.where(nan_bool_array == True)[0]
            if len(idxs) > 0:
                nan_dict[param]=len(idxs)
                nan_idxs = np.append(nan_idxs, idxs)
        total_samps = len(self)
        nan_samps   = len(nan_idxs)
        if nan_samps is not 0:
            print("WARNING: removing %i of %i total samples due to NaNs:"% (nan_samps,total_samps))
            for param in nan_dict.keys():
                print("\t%i NaNs in %s."%(nan_dict[param],param))
            self.delete_samples_by_idx(nan_idxs)
        return

    @property
    def DIC(self):
        """Returns the Deviance Information Criterion estimated from the
        posterior samples.  The DIC is defined as -2*(<log(L)> -
        Var(log(L))); smaller values are "better."

        """

        return -2.0*(np.mean(self._logL) - np.var(self._logL))

    @property
    def injection(self):
        """
        Return the injected values.

        """

        return self._injection

    @property
    def triggers(self):
        """
        Return the trigger values .

        """

        return self._triggers

    def _total_incl_restarts(self, samples):
        total=0
        last=samples[0]
        for x in samples[1:]:
            if x < last:
                total += last
            last = x
        total += samples[-1]
        return total

    def longest_chain_cycles(self):
        """
        Returns the number of cycles in the longest chain

        """
        samps,header=self.samples()
        header=header.split()
        if not ('cycle' in header):
            raise RuntimeError("Cannot compute number of cycles in longest chain")

        cycle_col=header.index('cycle')
        if 'chain' in header:
            chain_col=header.index('chain')
            chain_indexes=np.unique(samps[:,chain_col])
            max_cycle=0
            for ind in chain_indexes:
                chain_cycle_samps=samps[ samps[:,chain_col] == ind, cycle_col ]
                max_cycle=max(max_cycle, self._total_incl_restarts(chain_cycle_samps))
            return int(max_cycle)
        else:
            return int(self._total_incl_restarts(samps[:,cycle_col]))

    #@injection.setter #Python 2.6+
    def set_injection(self,injection):
        """
        Set the injected values of the parameters.

        @param injection: A SimInspiralTable row object containing the injected parameters.
        """
        if injection is not None:
            self._injection=injection
            for name,onepos in self:
                new_injval=self._getinjpar(name)
                if new_injval is not None:
                    self[name].set_injval(new_injval)

    def set_triggers(self,triggers):
        """
        Set the trigger values of the parameters.

        @param triggers: A list of SnglInspiral objects.
        """
        if triggers is not None:
            self._triggers=triggers
            for name,onepos in self:
                new_trigvals=self._gettrigpar(name)
                if new_trigvals is not None:
                    self[name].set_trigvals(new_trigvals)


    def _getinjpar(self,paramname):
        """
        Map parameter names to parameters in a SimInspiralTable .
        """
        if self._injection is not None:
            for key,value in self._injXMLFuncMap.items():
                if paramname.lower().strip() == key.lower().strip():
                    try:
                        return self._injXMLFuncMap[key](self._injection)
                    except TypeError:
                        return self._injXMLFuncMap[key]
        return None

    def _gettrigpar(self,paramname):
        """
        Map parameter names to parameters in a SnglInspiral.
        """
        vals = None
        if self._triggers is not None:
            for key,value in self._injXMLFuncMap.items():
                if paramname.lower().strip() == key.lower().strip():
                    try:
                        vals = dict([(trig.ifo,self._injXMLFuncMap[key](trig)) for trig in self._triggers])
                    except TypeError:
                        return self._injXMLFuncMap[key]
                    except AttributeError:
                        return None
        return vals

    def __getitem__(self,key):
        """
        Container method . Returns posterior chain,one_d_pos, with name one_d_pos.name.
        """
        return self._posterior[key.lower()]

    def __len__(self):
        """
        Container method. Defined as number of samples.
        """
        return len(self._logL)

    def __iter__(self):
        """
        Container method. Returns iterator from self.forward for use in
        for (...) in (...) etc.
        """
        return self.forward()

    def forward(self):
        """
        Generate a forward iterator (in sense of list of names) over Posterior
        with name,one_d_pos.
        """
        current_item = 0
        while current_item < self.dim:
            name=list(self._posterior.keys())[current_item]
            pos=self._posterior[name]
            current_item += 1
            yield name,pos

    def bySample(self):
        """
        Generate a forward iterator over the list of samples corresponding to
        the data stored within the Posterior instance. These are returned as
        ParameterSamples instances.
        """
        current_item=0
        pos_array,header=self.samples
        while current_item < len(self):
            sample_array=(np.squeeze(pos_array[current_item,:]))
            yield PosteriorSample(sample_array, header, header)
            current_item += 1


    @property
    def dim(self):
        """
        Return number of parameters.
        """
        return len(self._posterior.keys())

    @property
    def names(self):
        """
        Return list of parameter names.
        """
        nameslist=[]
        for key,value in self:
            nameslist.append(key)
        return nameslist

    @property
    def means(self):
        """
        Return dict {paramName:paramMean} .
        """
        meansdict={}
        for name,pos in self:
            meansdict[name]=pos.mean
        return meansdict

    @property
    def medians(self):
        """
        Return dict {paramName:paramMedian} .
        """
        mediansdict={}
        for name,pos in self:
            mediansdict[name]=pos.median
        return mediansdict

    @property
    def stdevs(self):
        """
        Return dict {paramName:paramStandardDeviation} .
        """
        stdsdict={}
        for name,pos in self:
            stdsdict[name]=pos.stdev
        return stdsdict

    @property
    def name(self):
        """
        Return qualified string containing the 'name' of the Posterior instance.
        """
        return self.__name

    @property
    def description(self):
        """
        Return qualified string containing a 'description' of the Posterior instance.
        """
        return self.__description

    def append(self,one_d_posterior):
        """
        Container method. Add a new OneDParameter to the Posterior instance.
        """
        self._posterior[one_d_posterior.name]=one_d_posterior
        return

    def pop(self,param_name):
        """
        Container method.  Remove PosteriorOneDPDF from the Posterior instance.
        """
        return self._posterior.pop(param_name)

    def append_mapping(self, new_param_names, func, post_names):
        """
        Append posteriors pos1,pos2,...=func(post_names)
        """
        # deepcopy 1D posteriors to ensure mapping function doesn't modify the originals
        import copy
        #1D input
        if isinstance(post_names, str):
            old_post = copy.deepcopy(self[post_names])
            old_inj  = old_post.injval
            old_trigs  = old_post.trigvals
            if old_inj:
                new_inj = func(old_inj)
            else:
                new_inj = None
            if old_trigs:
                new_trigs = {}
                for IFO in old_trigs.keys():
                    new_trigs[IFO] = func(old_trigs[IFO])
            else:
                new_trigs = None

            samps = func(old_post.samples)
            new_post = PosteriorOneDPDF(new_param_names, samps, injected_value=new_inj, trigger_values=new_trigs)
            if new_post.samples.ndim is 0:
                print("WARNING: No posterior calculated for %s ..." % post.name)
            else:
                self.append(new_post)
        #MultiD input
        else:
            old_posts = [copy.deepcopy(self[post_name]) for post_name in post_names]
            old_injs = [post.injval for post in old_posts]
            old_trigs = [post.trigvals for post in old_posts]
            samps = func(*[post.samples for post in old_posts])
            #1D output
            if isinstance(new_param_names, str):
                if None not in old_injs:
                    inj = func(*old_injs)
                else:
                    inj = None
                if None not in old_trigs:
                    new_trigs = {}
                    for IFO in old_trigs[0].keys():
                        oldvals = [param[IFO] for param in old_trigs]
                        new_trigs[IFO] = func(*oldvals)
                else:
                    new_trigs = None
                new_post = PosteriorOneDPDF(new_param_names, samps, injected_value=inj, trigger_values=new_trigs)
                self.append(new_post)
            #MultiD output
            else:
                if None not in old_injs:
                    injs = func(*old_injs)
                else:
                    injs = [None for name in new_param_names]
                if None not in old_trigs:
                    new_trigs = [{} for param in range(len(new_param_names))]
                    for IFO in old_trigs[0].keys():
                        oldvals = [param[IFO] for param in old_trigs]
                        newvals = func(*oldvals)
                        for param,newval in enumerate(newvals):
                            new_trigs[param][IFO] = newval
                else:
                    new_trigs = [None for param in range(len(new_param_names))]
                if not samps: return() # Something went wrong
                new_posts = [PosteriorOneDPDF(new_param_name,samp,injected_value=inj,trigger_values=new_trigs) for (new_param_name,samp,inj,new_trigs) in zip(new_param_names,samps,injs,new_trigs)]
                for post in new_posts:
                    if post.samples.ndim is 0:
                        print("WARNING: No posterior calculated for %s ..." % post.name)
                    else:
                        self.append(post)
        return

    def _average_posterior(self, samples, post_name):
        """
        Returns the average value of the 'post_name' column of the
        given samples.
        """
        ap = 0.0
        for samp in samples:
            ap = ap + samp[post_name]
        return ap / len(samples)

    def _average_posterior_like_prior(self, samples, logl_name, prior_name, log_bias = 0):
        """
        Returns the average value of the posterior assuming that the
        'logl_name' column contains log(L) and the 'prior_name' column
        contains the prior (un-logged).
        """
        ap = 0.0
        for samp in samples:
            ap += np.exp(samp[logl_name]-log_bias)*samp[prior_name]
        return ap / len(samples)

    def _bias_factor(self):
        """
        Returns a sensible bias factor for the evidence so that
        integrals are representable as doubles.
        """
        return np.mean(self._logL)

    def di_evidence(self, boxing=64):
        """
        Returns the log of the direct-integration evidence for the
        posterior samples.
        """
        allowed_coord_names=["spin1", "spin2", "a1", "phi1", "theta1", "a2", "phi2", "theta2",
                             "iota", "psi", "ra", "dec",
                             "phi_orb", "phi0", "dist", "time", "mc", "mchirp", "chirpmass", "q"]
        samples,header=self.samples()
        header=header.split()
        coord_names=[name for name in allowed_coord_names if name in header]
        coordinatized_samples=[PosteriorSample(row, header, coord_names) for row in samples]
        tree=KDTree(coordinatized_samples)

        if "prior" in header and "logl" in header:
            bf = self._bias_factor()
            return bf + np.log(tree.integrate(lambda samps: self._average_posterior_like_prior(samps, "logl", "prior", bf), boxing))
        elif "prior" in header and "likelihood" in header:
            bf = self._bias_factor()
            return bf + np.log(tree.integrate(lambda samps: self._average_posterior_like_prior(samps, "likelihood", "prior", bf), boxing))
        elif "post" in header:
            return np.log(tree.integrate(lambda samps: self._average_posterior(samps, "post"), boxing))
        elif "posterior" in header:
            return np.log(tree.integrate(lambda samps: self._average_posterior(samps, "posterior"), boxing))
        else:
            raise RuntimeError("could not find 'post', 'posterior', 'logl' and 'prior', or 'likelihood' and 'prior' columns in output to compute direct integration evidence")

    def elliptical_subregion_evidence(self):
        """Returns an approximation to the log(evidence) obtained by
        fitting an ellipse around the highest-posterior samples and
        performing the harmonic mean approximation within the ellipse.
        Because the ellipse should be well-sampled, this provides a
        better approximation to the evidence than the full-domain HM."""
        allowed_coord_names=["spin1", "spin2", "a1", "phi1", "theta1", "a2", "phi2", "theta2",
                             "iota", "psi", "ra", "dec",
                             "phi_orb", "phi0", "dist", "time", "mc", "mchirp", "chirpmass", "q"]
        samples,header=self.samples()
        header=header.split()

        n=int(0.05*samples.shape[0])
        if not n > 1:
            raise IndexError

        coord_names=[name for name in allowed_coord_names if name in header]
        indexes=np.argsort(self._logL[:,0])

        my_samples=samples[indexes[-n:], :] # The highest posterior samples.
        my_samples=np.array([PosteriorSample(sample,header,coord_names).coord() for sample in my_samples])

        mu=np.mean(my_samples, axis=0)
        cov=np.cov(my_samples, rowvar=0)

        d0=None
        for mysample in my_samples:
            d=np.dot(mysample-mu, np.linalg.solve(cov, mysample-mu))
            if d0 is None:
                d0 = d
            else:
                d0=max(d0,d)

        ellipse_logl=[]
        ellipse_samples=[]
        for sample,logl in zip(samples, self._logL):
            coord=PosteriorSample(sample, header, coord_names).coord()
            d=np.dot(coord-mu, np.linalg.solve(cov, coord-mu))

            if d <= d0:
                ellipse_logl.append(logl)
                ellipse_samples.append(sample)

        if len(ellipse_samples) > 5*n:
            print('WARNING: ellpise evidence region encloses significantly more samples than %d'%n)

        ellipse_samples=np.array(ellipse_samples)
        ellipse_logl=np.array(ellipse_logl)

        ndim = len(coord_names)
        ellipse_volume=np.pi**(ndim/2.0)*d0**(ndim/2.0)/special.gamma(ndim/2.0+1)*np.sqrt(np.linalg.det(cov))

        try:
            prior_index=header.index('prior')
            pmu=np.mean(ellipse_samples[:,prior_index])
            pstd=np.std(ellipse_samples[:,prior_index])
            if pstd/pmu > 1.0:
                print('WARNING: prior variation greater than 100\% over elliptical volume.')
            approx_prior_integral=ellipse_volume*pmu
        except KeyError:
            # Maybe prior = 1?
            approx_prior_integral=ellipse_volume

        ll_bias=np.mean(ellipse_logl)
        ellipse_logl = ellipse_logl - ll_bias

        return np.log(approx_prior_integral) - np.log(np.mean(1.0/np.exp(ellipse_logl))) + ll_bias

    def harmonic_mean_evidence(self):
        """
        Returns the log of the harmonic mean evidence for the set of
        posterior samples.
        """
        bf = self._bias_factor()
        return bf + np.log(1/np.mean(1/np.exp(self._logL-bf)))

    def _posMaxL(self):
        """
        Find the sample with maximum likelihood probability. Returns value
        of likelihood and index of sample .
        """
        logl_vals=self._logL
        max_i=0
        max_logl=logl_vals[0]
        for i in range(len(logl_vals)):
            if logl_vals[i] > max_logl:
                max_logl=logl_vals[i]
                max_i=i
        return max_logl,max_i

    def _posMap(self):
        """
        Find the sample with maximum a posteriori probability. Returns value
        of posterior and index of sample .
        """
        logl_vals=self._logL
        if self._logP is not None:
            logp_vals=self._logP
        else:
            return None

        max_i=0
        max_pos=logl_vals[0]+logp_vals[0]
        for i in range(len(logl_vals)):
            if logl_vals[i]+logp_vals[i] > max_pos:
                max_pos=logl_vals[i]+logp_vals[i]
                max_i=i
        return max_pos,max_i

    def _print_table_row(self,name,entries):
        """
        Print a html table row representation of

        name:item1,item2,item3,...
        """

        row_str='<tr><td>%s</td>'%name
        for entry in entries:
            row_str+='<td>%s</td>'%entry
        row_str+='</tr>'
        return row_str

    @property
    def maxL(self):
        """
        Return the maximum likelihood probability and the corresponding
        set of parameters.
        """
        maxLvals={}
        max_logl,max_i=self._posMaxL()
        for param_name in self.names:
            maxLvals[param_name]=self._posterior[param_name].samples[max_i][0]

        return (max_logl,maxLvals)

    @property
    def maxP(self):
        """
        Return the maximum a posteriori probability and the corresponding
        set of parameters.
        """
        maxPvals={}
        max_pos,max_i=self._posMap()
        for param_name in self.names:
            maxPvals[param_name]=self._posterior[param_name].samples[max_i][0]

        return (max_pos,maxPvals)


    def samples(self):
        """
        Return an (M,N) numpy.array of posterior samples; M = len(self);
        N = dim(self) .
        """
        header_string=''
        posterior_table=[]
        for param_name,one_pos in self:
            column=np.array(one_pos.samples)
            header_string+=param_name+'\t'
            posterior_table.append(column)
        posterior_table=tuple(posterior_table)
        return np.column_stack(posterior_table),header_string

    def write_to_file(self,fname):
        """
        Dump the posterior table to a file in the 'common format'.
        """
        posterior_table, header_string = self.samples()
        np.savetxt(
            fname,
            posterior_table,
            comments='',
            delimiter='\t',
            header=header_string,
        )

    def gelman_rubin(self, pname):
        """
        Returns an approximation to the Gelman-Rubin statistic (see
        Gelman, A. and Rubin, D. B., Statistical Science, Vol 7,
        No. 4, pp. 457--511 (1992)) for the parameter given, accurate
        as the number of samples in each chain goes to infinity.  The
        posterior samples must have a column named 'chain' so that the
        different chains can be separated.
        """
        from numpy import seterr as np_seterr
        np_seterr(all='raise')

        if "chain" in self.names:
            chains=np.unique(self["chain"].samples)
            chain_index=self.names.index("chain")
            param_index=self.names.index(pname)
            data,header=self.samples()
            chainData=[data[ data[:,chain_index] == chain, param_index] for chain in chains]
            allData=np.concatenate(chainData)
            chainMeans=[np.mean(data) for data in chainData]
            chainVars=[np.var(data) for data in chainData]
            BoverN=np.var(chainMeans)
            W=np.mean(chainVars)
            sigmaHat2=W + BoverN
            m=len(chainData)
            VHat=sigmaHat2 + BoverN/m
            try:
                R = VHat/W
            except:
                print("Error when computer Gelman-Rubin R statistic for %s.  This may be a fixed parameter"%pname)
                R = np.nan
            return R
        else:
            raise RuntimeError('could not find necessary column header "chain" in posterior samples')

    def healpix_map(self, resol, nest=True):
        """Returns a healpix map in the pixel ordering that represents the
        posterior density (per square degree) on the sky.  The pixels
        will be chosen to have at least the given resolution (in
        degrees).

        """

        # Ensure that the resolution is twice the desired
        nside = 2
        while hp.nside2resol(nside, arcmin=True) > resol*60.0/2.0:
            nside *= 2

        ras = self['ra'].samples.squeeze()
        decs = self['dec'].samples.squeeze()

        phis = ras
        thetas = np.pi/2.0 - decs

        # Create the map in ring ordering
        inds = hp.ang2pix(nside, thetas, phis, nest=False)
        counts = np.bincount(inds)
        if counts.shape[0] < hp.nside2npix(nside):
            counts = np.concatenate((counts, np.zeros(hp.nside2npix(nside) - counts.shape[0])))

        # Smooth the map a bit (Gaussian sigma = resol)
        hpmap = hp.sphtfunc.smoothing(counts, sigma=resol*np.pi/180.0)

        hpmap = hpmap / (np.sum(hpmap)*hp.nside2pixarea(nside, degrees=True))

        if nest:
            hpmap = hp.reorder(hpmap, r2n=True)

        return hpmap

    def __str__(self):
        """
        Define a string representation of the Posterior class ; returns
        a html formatted table of various properties of posteriors.
        """
        return_val='<table border="1" id="statstable"><tr><th/>'
        column_names=['maP','maxL','stdev','mean','median','stacc','injection value']
        IFOs = []
        if self._triggers is not None:
            IFOs = [trig.ifo for trig in self._triggers]
            for IFO in IFOs:
                column_names.append(IFO+' trigger values')

        for column_name in column_names:
            return_val+='<th>%s</th>'%column_name

        return_val+='</tr>'

        for name,oned_pos in self:

            max_logl,max_i=self._posMaxL()
            maxL=oned_pos.samples[max_i][0]
            max_post,max_j=self._posMap()
            maP=oned_pos.samples[max_j][0]
            mean=str(oned_pos.mean)
            stdev=str(oned_pos.stdev)
            median=str(np.squeeze(oned_pos.median))
            stacc=str(oned_pos.stacc)
            injval=str(oned_pos.injval)
            trigvals=oned_pos.trigvals

            row = [maP,maxL,stdev,mean,median,stacc,injval]
            if self._triggers is not None:
                for IFO in IFOs:
                    try:
                        row.append(str(trigvals[IFO]))
                    except TypeError:
                        row.append(None)
            return_val+=self._print_table_row(name,row)

        return_val+='</table>'

        parser=XMLParser()
        parser.feed(return_val)
        Estr=parser.close()

        elem=Estr
        rough_string = tostring(elem, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        return_val=reparsed.toprettyxml(indent="  ")
        return return_val[len('<?xml version="1.0" ?>')+1:]


    #===============================================================================
    # Functions used to parse injection structure.
    #===============================================================================
    def _inj_m1(self,inj):
        """
        Return the mapping of (mchirp,eta)->m1; m1>m2 i.e. return the greater of the mass
        components (m1) calculated from the chirp mass and the symmetric mass ratio.

        @param inj: a custom type with the attributes 'mchirp' and 'eta'.
        """
        (mass1,mass2)=mc2ms(inj.mchirp,inj.eta)
        return mass1

    def _inj_m2(self,inj):
        """
        Return the mapping of (mchirp,eta)->m2; m1>m2 i.e. return the lesser of the mass
        components (m2) calculated from the chirp mass and the symmetric mass ratio.

        @param inj: a custom type with the attributes 'mchirp' and 'eta'.
        """
        (mass1,mass2)=mc2ms(inj.mchirp,inj.eta)
        return mass2

    def _inj_q(self,inj):
        """
        Return the mapping of (mchirp,eta)->q; m1>m2 i.e. return the mass ratio q=m2/m1.

        @param inj: a custom type with the attributes 'mchirp' and 'eta'.
        """
        (mass1,mass2)=mc2ms(inj.mchirp,inj.eta)
        return mass2/mass1

    def _inj_longitude(self,inj):
        """
        Return the mapping of longitude found in inj to the interval [0,2*pi).

        @param inj: a custom type with the attribute 'longitude'.
        """
        if inj.longitude>2*pi_constant or inj.longitude<0.0:
            maplong=2*pi_constant*(((float(inj.longitude))/(2*pi_constant)) - floor(((float(inj.longitude))/(2*pi_constant))))
            print("Warning: Injected longitude/ra (%s) is not within [0,2\pi)! Angles are assumed to be in radians so this will be mapped to [0,2\pi). Mapped value is: %s."%(str(inj.longitude),str(maplong)))
            return maplong
        else:
            return inj.longitude

    def _inj_spins(self, inj, frame='OrbitalL'):

        from lalsimulation import SimInspiralTransformPrecessingWvf2PE

        spins = {}
        f_ref = self._injFref

        if not inj:
            spins = {}

        else:
            axis = lalsim.SimInspiralGetFrameAxisFromString(frame)
            s1x=inj.spin1x
            s1y=inj.spin1y
            s1z=inj.spin1z
            s2x=inj.spin2x
            s2y=inj.spin2y
            s2z=inj.spin2z
            iota=inj.inclination
            m1, m2 = inj.mass1, inj.mass2
            mc, eta = inj.mchirp, inj.eta

            a1, theta1, phi1 = cart2sph(s1x, s1y, s1z)
            a2, theta2, phi2 = cart2sph(s2x, s2y, s2z)

            spins = {'a1':a1, 'theta1':theta1, 'phi1':phi1,
                     'a2':a2, 'theta2':theta2, 'phi2':phi2,
                     'iota':iota}
            # If spins are aligned, save the sign of the z-component
            if inj.spin1x == inj.spin1y == inj.spin2x == inj.spin2y == 0.:
                spins['a1z'] = inj.spin1z
                spins['a2z'] = inj.spin2z

            L  = orbital_momentum(f_ref, m1,m2, iota)
            S1 = np.hstack((s1x, s1y, s1z))
            S2 = np.hstack((s2x, s2y, s2z))

            zhat = np.array([0., 0., 1.])
            aligned_comp_spin1 = array_dot(S1, zhat)
            aligned_comp_spin2 = array_dot(S2, zhat)
            chi = aligned_comp_spin1 + aligned_comp_spin2 + \
                  np.sqrt(1. - 4.*eta) * (aligned_comp_spin1 - aligned_comp_spin2)
            S1 *= m1**2
            S2 *= m2**2
            J = L + S1 + S2

            beta  = array_ang_sep(J, L)
            spins['beta'] = beta
            spins['spinchi'] = chi
            # Huge caveat: SimInspiralTransformPrecessingWvf2PE assumes that the cartesian spins in the XML table  are given in the L frame, ie. in  a frame where L||z. While this is the default in inspinj these days, other possibilities exist.
            # Unfortunately, we don't have a function (AFIK), that transforms spins from an arbitrary  frame to an arbitrary frame, otherwise I'd have called it here to be sure we convert in the L frame.
            # FIXME: add that function here if it ever gets written. For the moment just check
            if not frame=='OrbitalL':
                print("I cannot calculate the injected values of the spin angles unless frame is OrbitalL. Skipping...")
                return spins
            # m1 and m2 here are NOT in SI, but in Msun, this is not a typo.
            theta_jn,phi_jl,tilt1,tilt2,phi12,chi1,chi2=SimInspiralTransformPrecessingWvf2PE(inj.inclination,inj.spin1x, inj.spin1y, inj.spin1z,inj.spin2x, inj.spin2y, inj.spin2z, m1, m2, f_ref, inj.coa_phase)
            spins['theta_jn']=theta_jn
            spins['phi12']=phi12
            spins['tilt1']=tilt1
            spins['tilt2']=tilt2
            spins['phi_jl']=phi_jl

            """
            #If everything is all right, this function should give back the cartesian spins. Uncomment to check
            print("Inverting ")
            iota_back,a1x_back,a1y_back,a1z_back,a2x_back,a2y_back,a2z_back = \
    lalsim.SimInspiralTransformPrecessingNewInitialConditions(theta_jn,phi_jl,tilt1,tilt2,phi12,chi1,chi2,m1*lal.MSUN_SI,m2*lal.MSUN_SI,f_ref,inj.coa_phase)
            print(a1x_back,a1y_back,a1z_back)
            print(a2x_back,a2y_back,a2z_back)
            print(iota_back)
            """

        return spins

class BurstPosterior(Posterior):
    """
    Data structure for a table of posterior samples .
    """
    def __init__(self,commonResultsFormatData,SimBurstTableEntry=None,injFref=None,SnglBurstList=None,name=None,description=None):
        """
        Constructor.

        @param commonResultsFormatData: A 2D array containing the posterior
            samples and related data. The samples chains form the columns.
        @param SimBurstTableEntry: A glue.ligolw.lscstables.SimBurst row containing the injected values.
        @param SnglBurstList: A list of SnglBurst objects containing the triggers.
        @param injFref: reference frequency in injection
        @param name: optional name for this Posterior
        @param description: optional description for this Posterior
        """
        common_output_table_header,common_output_table_raw =commonResultsFormatData
        self._posterior={}
        self._injFref=injFref
        self._injection=SimBurstTableEntry
        self._triggers=SnglBurstList
        self._loglaliases=['posterior', 'logl','logL','likelihood', 'deltalogl']
        self._logpaliases=['logp', 'logP','prior','logprior','Prior','logPrior']

        common_output_table_header=[i.lower() for i in common_output_table_header]

        # Define XML mapping
        self._injXMLFuncMap={
                            'f0':lambda inj:inj.frequency,
                            'frequency':lambda inj:inj.frequency,
                            'centre_frequency':lambda inj:inj.frequency,
                            'quality':lambda inj:inj.q,
                            'hrss':lambda inj:inj.hrss,
                            'loghrss':lambda inj:log(inj.hrss),
                            'polar_angle':lambda inj:inj.pol_ellipse_angle,
                            'pol_ellipse_angle':lambda inj:inj.pol_ellipse_angle,
                            'pol_ellipse_e':lambda inj:inj.pol_ellipse_e,
                            'alpha':lambda inj:inj.pol_ellipse_angle,
                            'polar_eccentricity':lambda inj:inj.pol_ellipse_e,
                            'eccentricity':lambda inj:inj.pol_ellipse_e,
                            'time': lambda inj:float(inj.get_end()),
                            'end_time': lambda inj:float(inj.get_end()),
                            'ra':self._inj_longitude,
                            'rightascension':self._inj_longitude,
                            'long':self._inj_longitude,
                            'longitude':self._inj_longitude,
                            'dec':lambda inj:inj.dec,
                            'declination':lambda inj:inj.dec,
                            'lat':lambda inj:inj.dec,
                            'latitude':lambda inj:inj.dec,
                            'psi': lambda inj: np.mod(inj.psi, np.pi),
                            'f_ref': lambda inj: self._injFref,
                            'polarisation':lambda inj:inj.psi,
                            'polarization':lambda inj:inj.psi,
                            'duration':lambda inj:inj.duration,
                            'h1_end_time':lambda inj:float(inj.get_end('H')),
                            'l1_end_time':lambda inj:float(inj.get_end('L')),
                            'v1_end_time':lambda inj:float(inj.get_end('V')),
                           }

        for one_d_posterior_samples,param_name in zip(np.hsplit(common_output_table_raw,common_output_table_raw.shape[1]),common_output_table_header):

            self._posterior[param_name]=PosteriorOneDPDF(param_name.lower(),one_d_posterior_samples,injected_value=self._getinjpar(param_name),injFref=self._injFref,trigger_values=self._gettrigpar(param_name))

        logLFound=False

        for loglalias in self._loglaliases:
            if loglalias in common_output_table_header:
                try:
                    self._logL=self._posterior[loglalias].samples
                except KeyError:
                    print("No '%s' column in input table!"%loglalias)
                    continue
                logLFound=True

        if not logLFound:
            raise RuntimeError("No likelihood/posterior values found!")

        self._logP=None
        for logpalias in self._logpaliases:
            if logpalias in common_output_table_header:
                try:
                    self._logP=self._posterior[logpalias].samples
                except KeyError:
                    print("No '%s' column in input table!"%logpalias)
                    continue
                if not 'log' in logpalias:
                    self._logP=[np.log(i) for i in self._logP]
        if name is not None:
            self.__name=name

        if description is not None:
            self.__description=description

        return
    #===============================================================================
    # Functions used to parse injection structure.
    #===============================================================================

    def _inj_longitude(self,inj):
        """
        Return the mapping of longitude found in inj to the interval [0,2*pi).

        @param inj: a custom type with the attribute 'longitude'.
        """
        if inj.ra>2*pi_constant or inj.ra<0.0:
            maplong=2*pi_constant*(((float(inj.ra)/(2*pi_constant)) - floor(((float(inj.ra))/(2*pi_constant)))))
            print("Warning: Injected longitude/ra (%s) is not within [0,2\pi)! Angles are assumed to be in radians so this will be mapped to [0,2\pi). Mapped value is: %s."%(str(inj.ra),str(maplong)))
            return maplong
        else:
            return inj.ra

class KDTree(object):
    """
    A kD-tree.
    """
    def __init__(self, objects):
        """
        Construct a kD-tree from a sequence of objects.  Each object
        should return its coordinates using obj.coord().
        """
        if len(objects) == 0:
            raise RuntimeError("cannot construct kD-tree out of zero objects---you may have a repeated sample in your list")
        elif len(objects) == 1:
            self._objects = objects[:]
            coord=self._objects[0].coord()
            self._bounds = coord,coord
        elif self._same_coords(objects):
            # All the same coordinates
            self._objects = [ objects[0] ]
            coord=self._objects[0].coord()
            self._bounds = coord,coord
        else:
            self._objects = objects[:]
            self._bounds = self._bounds_of_objects()
            low,high=self._bounds
            self._split_dim=self._longest_dimension()
            longest_dim = self._split_dim
            sorted_objects=sorted(self._objects, key=lambda obj: (obj.coord())[longest_dim])
            N = len(sorted_objects)
            bound=0.5*(sorted_objects[int(N/2)].coord()[longest_dim] + sorted_objects[int(N/2)-1].coord()[longest_dim])
            low = [obj for obj in self._objects if obj.coord()[longest_dim] < bound]
            high = [obj for obj in self._objects if obj.coord()[longest_dim] >= bound]
            if len(low)==0:
                # Then there must be multiple values with the same
                # coordinate as the minimum element of high
                low = [obj for obj in self._objects if obj.coord()[longest_dim]==bound]
                high = [obj for obj in self._objects if obj.coord()[longest_dim] > bound]
            self._left = KDTree(low)
            self._right = KDTree(high)

    def _same_coords(self, objects):
        """
        True if and only if all the given objects have the same
        coordinates.
        """
        if len(objects) <= 1:
            return True
        coords = [obj.coord() for obj in objects]
        c0 = coords[0]
        for ci in coords[1:]:
            if not np.all(ci == c0):
                return False
        return True

    def _bounds_of_objects(self):
        """
        Bounds of the objects contained in the tree.
        """
        low=self._objects[0].coord()
        high=self._objects[0].coord()
        for obj in self._objects[1:]:
            low=np.minimum(low,obj.coord())
            high=np.maximum(high,obj.coord())
        return low,high

    def _longest_dimension(self):
        """
        Longest dimension of the tree bounds.
        """
        low,high = self._bounds
        widths = high-low
        return np.argmax(widths)

    def objects(self):
        """
        Returns the objects in the tree.
        """
        return self._objects[:]

    def __iter__(self):
        """
        Iterator over all the objects contained in the tree.
        """
        return self._objects.__iter__()

    def left(self):
        """
        Returns the left tree.
        """
        return self._left

    def right(self):
        """
        Returns the right tree.
        """
        return self._right

    def split_dim(self):
        """
        Returns the dimension along which this level of the kD-tree
        splits.
        """
        return self._split_dim

    def bounds(self):
        """
        Returns the coordinates of the lower-left and upper-right
        corners of the bounding box for this tree: low_left, up_right
        """
        return self._bounds

    def volume(self):
        """
        Returns the volume of the bounding box of the tree.
        """
        v = 1.0
        low,high=self._bounds
        for l,h in zip(low,high):
            v = v*(h - l)
        return v

    def integrate(self,f,boxing=64):
        """
        Returns the integral of f(objects) over the tree.  The
        optional boxing parameter determines how deep to descend into
        the tree before computing f.
        """
        # if len(self._objects) <= boxing:
        #     return self.volume()*f(self._objects)
        # else:
        #     return self._left.integrate(f, boxing) + self._right.integrate(f, boxing)

        def x(tree):
            return tree.volume()*f(tree._objects)

        def y(a,b):
            return a+b

        return self.operate(x,y,boxing=boxing)

    def operate(self,f,g,boxing=64):
        """
        Operates on tree nodes exceeding boxing parameter depth.
        """
        if len(self._objects) <= boxing:
            return f(self)
        else:

            return g(self._left.operate(f,g,boxing),self._right.operate(f,g,boxing))


class KDTreeVolume(object):
    """
    A kD-tree suitable for splitting parameter spaces and counting hypervolumes.
    Is modified from the KDTree class so that bounding boxes are stored. This means that
    there are no longer gaps in the hypervolume once the samples have been split into groups.
    """
    def __init__(self, objects,boundingbox,dims=0):
        """
        Construct a kD-tree from a sequence of objects.  Each object
        should return its coordinates using obj.coord().
        the obj should also store the bounds of the hypervolume its found in.
        for non-leaf objects we need the name of the dimension split and value at split.
        """
        self._dimension = dims
        self._bounds = boundingbox
        self._weight = 1
        if len(objects) == 0: #for no objects - something is wrong, i think it can only happen in first call
            raise RuntimeError("cannot construct kD-tree out of zero objects---you may have a repeated sample in your list")
        elif len(objects) == 1: #1 object, have reached leaf of tree
            self._objects = objects[:]
        elif self._same_coords(objects): # When ALL samples have the same coordinates in all dimensions
            self._weight = len(objects)
            self._objects = [ objects[0] ] #need to modify kdtree_bin functions to use _weight to get correct number of samples
            coord=self._objects[0].coord()
        else: #construct next level of tree with multiple samples
            self._objects = objects[:]
            split_dim = self._dimension
            sorted_objects=sorted(self._objects, key=lambda obj: (obj.coord())[split_dim])
            N = len(sorted_objects)
            self._split_value = 0.5*(sorted_objects[int(N/2)].coord()[split_dim] + sorted_objects[int(N/2)-1].coord()[split_dim])
            bound = self._split_value
            low = [obj for obj in self._objects if obj.coord()[split_dim] < bound]
            high = [obj for obj in self._objects if obj.coord()[split_dim] >= bound]
            if len(low)==0:
                # Then there must be multiple values with the same
                # coordinate as the minimum element of 'high'
                low = [obj for obj in self._objects if obj.coord()[split_dim] == bound]
                high = [obj for obj in self._objects if obj.coord()[split_dim] > bound]
            leftBoundingbox = []
            rightBoundingbox = []
            for i in self._bounds:
                leftBoundingbox.append(list(i))
                rightBoundingbox.append(list(i))
            leftBoundingbox[1][split_dim] = bound
            rightBoundingbox[0][split_dim] = bound
            # designate the next dimension to use for split for sub-trees
            # if has got to the end of the list of dimensions then starts
            # again at dimension = 0
            if (split_dim < (len(self._objects[0].coord()) - 1)):
                child_dim = split_dim + 1
            else:
                child_dim = 0
            self._left = KDTreeVolume(low,leftBoundingbox,dims = child_dim)
            # added in a load of messing about incase there are repeated values in the currently checked dimension
            if (len(high) != 0):
                self._right = KDTreeVolume(high,rightBoundingbox,dims = child_dim)
            else:
                self._right = None

    def _same_coords(self, objects):
        """
        True if and only if all the given objects have the same
        coordinates.
        """
        if len(objects) <= 1:
            return True
        coords = [obj.coord() for obj in objects]
        c0 = coords[0]
        for ci in coords[1:]:
            if not np.all(ci == c0):
                return False
        return True

    def objects(self):
        """
        Returns the objects in the tree.
        """
        return self._objects[:]

    def __iter__(self):
        """
        Iterator over all the objects contained in the tree.
        """
        return self._objects.__iter__()

    def left(self):
        """
        Returns the left tree.
        """
        return self._left

    def right(self):
        """
        Returns the right tree.
        """
        return self._right

    def split_dim(self):
        """
        Returns the dimension along which this level of the kD-tree
        splits.
        """
        return self._split_dim

    def bounds(self):
        """
        Returns the coordinates of the lower-left and upper-right
        corners of the bounding box for this tree: low_left, up_right
        """
        return self._bounds

    def volume(self):
        """
        Returns the volume of the bounding box of the tree.
        """
        v = 1.0
        low,high=self._bounds
        for l,h in zip(low,high):
            v = v*(h - l)
        return v

    def integrate(self,f,boxing=64):
        """
        Returns the integral of f(objects) over the tree.  The
        optional boxing parameter determines how deep to descend into
        the tree before computing f.
        """
        def x(tree):
            return tree.volume()*f(tree._objects)

        def y(a,b):
            return a+b

        return self.operate(x,y,boxing=boxing)

    def operate(self,f,g,boxing=64):
        """
        Operates on tree nodes exceeding boxing parameter depth.
        """
        if len(self._objects) <= boxing:
            return f(self)
        else:
            return g(self._left.operate(f,g,boxing),self._right.operate(f,g,boxing))

    def search(self,coordinates,boxing = 64):
        """
        takes a set of coordinates and searches down through the tree untill it gets
        to a box with less than 'boxing' objects in it and returns the box bounds,
        number of objects in the box, and the weighting.
        """
        if len(self._objects) <= boxing:
            return self._bounds,len(self._objects),self._weight
        elif coordinates[self._dimension] < self._split_value:
            return self._left.search(coordinates,boxing)
        else:
            return self._right.search(coordinates,boxing)

    def fillNewTree(self,boxing = 64, isArea = False):
        """
        copies tree structure, but with KDSkeleton as the new nodes.
        """
        boxN = boxing
        if len(self._objects) <= boxN:
            newNode = KDSkeleton(self.bounds(), left_child = None , right_child = None)
            if isArea:
                newNode.setImportance(len(self._objects),skyArea(self.bounds()))
            else:
                newNode.setImportance(len(self._objects),self.volume())
            return newNode
        else:
            if isArea:
                newNode = KDSkeleton(self.bounds, left_child = self._left.fillNewTree(boxN,isArea=True), right_child = self._right.fillNewTree(boxN,isArea=True))
                newNode.setImportance(len(self._objects),skyArea(self.bounds()))
            else:
                newNode = KDSkeleton(self.bounds, left_child = self._left.fillNewTree(boxN), right_child = self._right.fillNewTree(boxN))
                newNode.setImportance(len(self._objects),self.volume())
            newNode.setSplit(self._dimension,self._split_value)
            return newNode

class KDSkeleton(object):
    """
    object to store the structure of a kd tree
    """

    def __init__(self, bounding_box, left_child = None, right_child = None):
        self._bounds = bounding_box
        #self._names = coordinate_names
        self._left = left_child
        self._right = right_child
        self._samples = 0
        self._splitValue = None
        self._splitDim = None
        self._importance = None
        self._volume = None

    def addSample(self):
        self._samples +=1

    def bounds(self):
        return self._bounds

    def search(self,coordinates):
        """
        takes a set of coordinates and searches down through the tree untill it gets
        to a box with less than 'boxing' objects in it and returns the box bounds,
        number of objects in the box, and the weighting.
        """
        if self._left is None:
            return self._bounds, self._samples, self._importance
        elif coordinates[self._splitDim] < self._splitValue:
            return self._left.search(coordinates)
        else:
            return self._right.search(coordinates)

    def setImportance(self, sampleNumber, volume):
        self._importance = sampleNumber/volume
        self._volume = volume

    def setSplit(self,dimension,value):
        self._splitDim = dimension
        self._splitValue = value


class PosteriorSample(object):
    """
    A single parameter sample object, suitable for inclusion in a
    kD-tree.
    """

    def __init__(self, sample_array, headers, coord_names):
        """
        Given the sample array, headers for the values, and the names
        of the desired coordinates, construct a parameter sample
        object.
        """
        self._samples=sample_array[:]
        self._headers=headers
        if not (len(sample_array) == len(self._headers)):
            print("Header length = ", len(self._headers))
            print("Sample length = ", len(sample_array))
            raise RuntimeError("parameter and sample lengths do not agree")
        self._coord_names=coord_names
        self._coord_indexes=[self._headers.index(name) for name in coord_names]

    def __getitem__(self, key):
        """
        Return the element with the corresponding name.
        """
        key=key.lower()
        if key in self._headers:
            idx=self._headers.index(key)
            return self._samples[idx]
        else:
            raise KeyError("key not found in posterior sample: %s"%key)

    def coord(self):
        """
        Return the coordinates for the parameter sample.
        """
        return self._samples[self._coord_indexes]





class AnalyticLikelihood(object):
    """
    Return analytic likelihood values.
    """

    def __init__(self, covariance_matrix_files, mean_vector_files):
        """
        Prepare analytic likelihood for the given parameters.
        """
        # Make sure files names are in a list
        if isinstance(covariance_matrix_files, str):
            covariance_matrix_files = [covariance_matrix_files]
        if isinstance(mean_vector_files, str):
            mean_vector_files = [mean_vector_files]

        covarianceMatrices = [np.loadtxt(csvFile, delimiter=',') for csvFile in covariance_matrix_files]
        num_matrices = len(covarianceMatrices)

        if num_matrices != len(mean_vector_files):
            raise RuntimeError('Must give a mean vector list for every covariance matrix')

        param_line = open(mean_vector_files[0]).readline()
        self._params = [param.strip() for param in param_line.split(',')]

        converter=lambda x: eval(x.replace('pi','%.32f'%pi_constant))  # converts fractions w/ pi (e.g. 3.0*pi/2.0)
        self._modes = []
        for i in range(num_matrices):
            CM = covarianceMatrices[i]
            vecFile = mean_vector_files[i]

            param_line = open(vecFile).readline()
            params = [param.strip() for param in param_line.split(',')]
            if set(params)!=set(self._params):
                raise RuntimeError('Parameters do not agree between mean vector files.')

            sigmas = dict(zip(params,np.sqrt(CM.diagonal())))
            colNums = range(len(params))
            converters = dict(zip(colNums,[converter for i in colNums]))
            meanVectors = np.loadtxt(vecFile, delimiter=',', skiprows=1, converters=converters)
            try:
                for vec in meanVectors:
                    means = dict(zip(params,vec))
                    mode = [(param, stats.norm(loc=means[param],scale=sigmas[param])) for param in params]
                    self._modes.append(dict(mode))
            except TypeError:
                means = dict(zip(params,meanVectors))
                mode = [(param, stats.norm(loc=means[param],scale=sigmas[param])) for param in params]
                self._modes.append(dict(mode))
        self._num_modes = len(self._modes)

    def pdf(self, param):
        """
        Return PDF function for parameter.
        """
        pdf = None
        if param in self._params:
            pdf = lambda x: (1.0/self._num_modes) * sum([mode[param].pdf(x) for mode in self._modes])
        return pdf

    def cdf(self, param):
        """
        Return PDF function for parameter.
        """
        cdf = None
        if param in self._params:
            cdf = lambda x: (1.0/self._num_modes) * sum([mode[param].cdf(x) for mode in self._modes])
        return cdf

    @property
    def names(self):
        """
        Return list of parameter names described by analytic likelihood function.
        """
        return self._params



#===============================================================================
# Web page creation classes (wrap ElementTrees)
#===============================================================================

class htmlChunk(object):
    """
    A base class for representing web content using ElementTree .
    """
    def __init__(self,tag,attrib=None,parent=None):

        self._html=Element(tag)#attrib={'xmlns':"http://www.w3.org/1999/xhtml"})
        if attrib:
            for attribname,attribvalue in attrib.items():
                self._html.attrib[attribname]=attribvalue
        if parent:
            parent.append(self._html)

    def toprettyxml(self):
        """
        Return a pretty-printed XML string of the htmlPage.
        """
        elem=self._html
        rough_string = tostring(elem)
        reparsed = minidom.parseString(rough_string)
        return reparsed.toprettyxml(indent="  ")

    def __str__(self):
        return self.toprettyxml()

    def write(self,string):
        parser=XMLParser()
        parser.feed(string)
        Estr=parser.close()
        self._html.append(Estr)

    def p(self,pstring):
        Ep=Element('p')
        Ep.text=pstring
        self._html.append(Ep)
        return Ep

    def h1(self,h1string):
        Ep=Element('h1')
        Ep.text=h1string
        self._html.append(Ep)
        return Ep
#
    def h5(self,h1string):
        Ep=Element('h5')
        Ep.text=h1string
        self._html.append(Ep)
        return Ep

    def h2(self,h2string):
        Ep=Element('h2')
        Ep.text=h2string
        self._html.append(Ep)
        return Ep

    def h3(self,h1string):
        Ep=Element('h3')
        Ep.text=h1string
        self._html.append(Ep)
        return Ep

    def br(self):
        Ebr=Element('br')
        self._html.append(Ebr)
        return Ebr

    def hr(self):
        Ehr=Element('hr')
        self._html.append(Ehr)
        return Ehr

    def a(self,url,linktext):
        Ea=Element('a')
        Ea.attrib['href']=url
        Ea.text=linktext
        self._html.append(Ea)
        return Ea

    def tab(self,idtable=None):
        args={}
        if idtable is not None:
            args={'id':idtable}

        Etab=Element('table',args)
        self._html.append(Etab)
        return Etab

    def insert_row(self,tab,label=None):

        """
        Insert row in table tab.
        If given, label used as id for the table tag
        """

        Etr=Element('tr')
        if label is not None:
            Etr.attrib['id']=label
        tab.append(Etr)
        return Etr

    def insert_td(self,row,td,label=None,legend=None):
        """
        Insert cell td into row row.
        Sets id to label, if given
        """

        Etd=Element('td')

        if type(td) is str:
            Etd.text=td
        else:
            td=tostring(td)
            td=minidom.parseString(td)
            td=td.toprettyxml(indent="  ")
            Etd.text=td
        if label is not None:
            Etd.attrib['id']=label
        if legend is not None:
            legend.a('#%s'%label,'%s'%label)
            legend.br()
        row.append(Etd)
        return Etd

    def append(self,element):
        self._html.append(element)


#
class htmlPage(htmlChunk):
    """
    A concrete class for generating an XHTML(1) document. Inherits from htmlChunk.
    """
    def __init__(self,title=None,css=None,javascript=None,toc=False):
        htmlChunk.__init__(self,'html',attrib={'xmlns':"http://www.w3.org/1999/xhtml"})
        self.doctype_str='<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">'

        self._head=SubElement(self._html,'head')
        Etitle=SubElement(self._head,'title')
        self._body=SubElement(self._html,'body')
        self._css=None
        self._jscript=None
        if title is not None:
            Etitle.text=str(title)
            self._title=SubElement(self._body,'h1')
            self._title.text=title
        if javascript is not None:
            self._jscript=SubElement(self._head,'script')
            self._jscript.attrib['type']="text/javascript"
            self._jscript.text=str(javascript)
        if css is not None:
            self._css=SubElement(self._head,'style')
            self._css.attrib['type']="text/css"
            self._css.text=str(css)

    def __str__(self):
        return self.doctype_str+'\n'+self.toprettyxml()

    def add_section(self,section_name,legend=None):
        newSection=htmlSection(section_name)
        self._body.append(newSection._html)
        if legend is not None:
            legend.a('#%s'%section_name,'%s'%section_name)
            legend.br()
        return newSection

    def add_collapse_section(self,section_name,legend=None,innertable_id=None,start_closed=True):
        """
        Create a section embedded into a table that can be collapsed with a button
        """
        newSection=htmlCollapseSection(section_name,table_id=innertable_id,start_closed=start_closed)
        self._body.append(newSection._html)
        if legend is not None:
            legend.a('#%s'%section_name,'%s'%section_name)
            legend.br()
        return newSection

    def add_section_to_element(self,section_name,parent):
        """
        Create a section which is not appended to the body of html, but to the parent Element
        """
        newSection=htmlSection(section_name,htmlElement=parent,blank=True)
        parent.append(newSection._html)
        return newSection


    @property
    def body():
        return self._body

    @property
    def head():
        return self._head


class htmlSection(htmlChunk):
    """
    Represents a block of html fitting within a htmlPage. Inherits from htmlChunk.
    """
    def __init__(self,section_name,htmlElement=None,blank=False):
        if not blank:
            htmlChunk.__init__(self,'div',attrib={'class':'ppsection','id':section_name},parent=htmlElement)
        else:
            htmlChunk.__init__(self,'div',attrib={'style':'"color:#000000"','id':section_name},parent=htmlElement)
        self.h2(section_name)


class htmlCollapseSection(htmlChunk):
    """
    Represents a block of html fitting within a htmlPage. Inherits from htmlChunk.
    """

    def __init__(self,section_name,htmlElement=None,table_id=None,start_closed=True):
        htmlChunk.__init__(self,'div',attrib={'class':'ppsection','id':section_name},parent=htmlElement)
        # if table id is none, generate a random id:
        if table_id is None:
            table_id=random.randint(1,10000000)
        self.table_id=table_id
        self._start_closed=start_closed

    def write(self,string):
        k=random.randint(1,10000000)
        if self._start_closed:
            st='<table border="0" align="center" cellpadding="5" cellspacing="0"><tr bgcolor="#4682B4" height="50"><td width="5%%"><font size="4" face="tahoma" color="white"><a href="#"> Top</a></font></td><td width="45%%"><font size="4" face="tahoma" color="white"><strong>%s</strong></font></td><td bgcolor="#4682B4" align="center" width="50%%"><input id="lnk%s" type="button" value="[+] Expand" onclick="toggle_visibility(\'%s\',\'lnk%s\');"></input></td></tr><tr><td colspan="7">'%(self._html.attrib['id'],k,self.table_id,k)
            string=string.replace('table ', 'table style="display: none" ')
        else:
            st='<table border="0" align="center" cellpadding="5" cellspacing="0"><tr bgcolor="#4682B4" height="50"><td width="5%%"><font size="4" face="tahoma" color="white"><a href="#"> Top</a></font></td><td width="45%%"><font size="4" face="tahoma" color="white"><strong>%s</strong></font></td><td bgcolor="#4682B4" align="center" width="50%%"><input id="lnk%s" type="button" value="[-] Collapse" onclick="toggle_visibility(\'%s\',\'lnk%s\');"></input></td></tr><tr><td colspan="7">'%(self._html.attrib['id'],k,self.table_id,k)
            string=string.replace('table ', 'table style="display: table" ')
        st+=string
        st+='</td></tr></table>'
        htmlChunk.write(self,st)

#===============================================================================
# Internal module functions
#===============================================================================


def _calculate_confidence_levels(hist, points, injBin, NSamples):
    """
    Returns (injectionconf, toppoints), where injectionconf is the
    confidence level of the injection, contained in the injBin and
    toppoints is a list of (pointx, pointy, ptindex, frac), with
    pointx and pointy the (x,y) coordinates of the corresponding
    element of the points array, ptindex the index of the point in the
    array, and frac the cumulative fraction of points with larger
    posterior probability.

    The hist argument should be a one-dimensional array that contains
    counts of sample points in each bin.

    The points argument should be a 2-D array storing the sky location
    associated with each bin; the first index runs from 0 to NBins -
    1, while the second index runs from 0 to 1.

    The injBin argument gives the bin index in which the injection is
    found.

    The NSamples argument is used to normalize the histogram counts
    into fractional probability.
    """

    histIndices=np.argsort(hist)[::-1]  # In decreasing order

    toppoints=[]
    frac=0.0
    injConf=None
    for i in histIndices:
        frac+=float(hist[i])/float(NSamples)
        toppoints.append((points[i,0], points[i,1], i, frac))
        if i == injBin:
            injConf=frac
            print('Injection found at confidence level %g'%injConf)

    return (injConf, toppoints)

def _greedy_bin(greedyHist,greedyPoints,injection_bin_index,bin_size,Nsamples,confidence_levels):
    """
    An interal function representing the common, dimensionally-independent part of the
    greedy binning algorithms.
    """

    #Now call confidence level C extension function to determine top-ranked pixels
    (injectionconfidence,toppoints)=_calculate_confidence_levels(greedyHist, greedyPoints, injection_bin_index, Nsamples)

    #Determine interval/area contained within given confidence intervals
    nBins=0
    confidence_levels.sort()
    reses={}
    toppoints=np.array(toppoints)
    for printcl in confidence_levels:
        nBins=np.searchsorted(toppoints[:,3], printcl) + 1

        if nBins >= len(toppoints):
            nBins=len(toppoints)-1

        accl=toppoints[nBins-1,3]

        reses[printcl]=nBins*bin_size

    #Find area
    injection_area=None
    if injection_bin_index and injectionconfidence:
        i=list(np.nonzero(np.asarray(toppoints)[:,2]==injection_bin_index))[0]
        injection_area=bin_size*(i+1)

    return toppoints,injectionconfidence,reses,injection_area
#
#### functions used in 2stage kdtree

def skyArea(bounds):
    return - (cos(pi_constant/2. - bounds[0][1])-cos(pi_constant/2. - bounds[1][1]))*(bounds[1][0] - bounds[0][0])

def random_split(items, fraction):
    size = int(len(items)*fraction)
    random.shuffle(items)
    return items[:size], items[size:]

def addSample(tree,coordinates):
    if tree._left is None:
        tree.addSample()
    elif coordinates[tree._splitDim] < tree._splitValue:
        addSample(tree._left,coordinates)
    else:
        addSample(tree._right,coordinates)

#
#===============================================================================
# Public module functions
#===============================================================================

def kdtree_bin_sky_volume(posterior,confidence_levels):

    confidence_levels.sort()

    class Harvester(list):

        def __init__(self):
            list.__init__(self)
            self.unrho=0.

        def __call__(self,tree):
            number_density=float(len(tree.objects()))/float(tree.volume())
            self.append([number_density,tree.volume(),tree.bounds()])
            self.unrho+=number_density

        def close_ranks(self):

            for i in range(len(self)):
                self[i][0]/=self.unrho

            return sorted(self,key=itemgetter(0))

    def h(a,b):
        pass

    samples,header=posterior.samples()
    header=header.split()
    coord_names=["ra","dec","dist"]
    coordinatized_samples=[PosteriorSample(row, header, coord_names) for row in samples]
    tree=KDTree(coordinatized_samples)

    a=Harvester()
    samples_per_bin=10
    tree.operate(a,h,boxing=samples_per_bin)

    b=a.close_ranks()
    b.reverse()

    acc_rho=0.
    acc_vol=0.
    cl_idx=0
    confidence_intervals={}
    for rho,vol,bounds in b:
        acc_rho+=rho
        acc_vol+=vol

        if acc_rho>confidence_levels[cl_idx]:
            confidence_intervals[acc_rho]=acc_vol
            cl_idx+=1
            if cl_idx==len(confidence_levels):
                break

    return confidence_intervals

def kdtree_bin_sky_area(posterior,confidence_levels,samples_per_bin=10):
    """
    takes samples and applies a KDTree to them to return confidence levels
    returns confidence_intervals - dictionary of user_provided_CL:calculated_area
            b - ordered list of KD leaves
            injInfo - if injection values provided then returns
                      [Bounds_of_inj_kd_leaf ,number_samples_in_box, weight_of_box,injection_CL ,injection_CL_area]
    Not quite sure that the repeated samples case is fixed, posibility of infinite loop.
    """
    confidence_levels.sort()
    from math import cos, pi
    class Harvester(list):
        """
        when called by kdtree.operate will be used to calculate the density of each bin (sky area)
        """
        def __init__(self):
            list.__init__(self)
            self.unrho=0.

        def __call__(self,tree):
            #area = (cos(tree.bounds()[0][1])-cos(tree.bounds()[1][1]))*(tree.bounds()[1][0] - tree.bounds()[0][0])
            area = - (cos(pi/2. - tree.bounds()[0][1])-cos(pi/2. - tree.bounds()[1][1]))*(tree.bounds()[1][0] - tree.bounds()[0][0])
            number_density=float(len(tree.objects()))/float(area)
            self.append([number_density,len(tree.objects()),area,tree.bounds()])
            self.unrho+=number_density

        def close_ranks(self):

            for i in range(len(self)):
                self[i][0]/=self.unrho

            return sorted(self,key=itemgetter(0))

    def h(a,b):
        pass

    peparser=PEOutputParser('common')

    samples,header=posterior.samples()
    header=header.split()
    coord_names=["ra","dec"]
    initial_dimensions = [[0.,-pi/2.],[2.*pi, pi/2.]]
    coordinatized_samples=[PosteriorSample(row, header, coord_names) for row in samples]
    tree=KDTreeVolume(coordinatized_samples,initial_dimensions)

    a=Harvester()
    tree.operate(a,h,boxing=samples_per_bin)
    totalSamples = len(tree.objects())
    b=a.close_ranks()
    b.reverse()
    samplecounter=0.0
    for entry in b:
        samplecounter += entry[1]
        entry[1] = float(samplecounter)/float(totalSamples)

    acc_rho=0.
    acc_vol=0.
    cl_idx=0

    #checks for injection and extract details of the node in the tree that the injection is found
    if posterior['ra'].injval is not None and posterior['dec'].injval is not None:
        injBound,injNum,injWeight = tree.search([posterior['ra'].injval,posterior['dec'].injval],boxing = samples_per_bin)
        injInfo = [injBound,injNum,injWeight]
        inj_area = - (cos(pi/2. - injBound[0][1])-cos(pi/2. - injBound[1][1]))*(injBound[1][0] - injBound[0][0])
        inj_number_density=float(injNum)/float(inj_area)
        inj_rho = inj_number_density / a.unrho
    else:
        injInfo = None
        inj_area = None
        inj_number_density=None
        inj_rho = None

    #finds the volume contained within the confidence levels requested by user
    confidence_intervals={}
    for rho,confidence_level,vol,bounds in b:
        acc_vol+=vol

        if confidence_level>confidence_levels[cl_idx]:
            print(str(confidence_level))
            print(acc_vol)
            confidence_intervals[confidence_levels[cl_idx]]=acc_vol
            cl_idx+=1
            if cl_idx==len(confidence_levels):
                break

    acc_vol = 0.
    for rho,sample_number,vol,bounds in b:
        acc_vol+=vol
    print('total area: ' + str(acc_vol))

    #finds the confidence level of the injection and the volume of the associated contained region
    inj_confidence = None
    inj_confidence_area = None
    if inj_rho is not None:
        acc_vol=0.
        for rho,confidence_level,vol,bounds in b:
            acc_vol+=vol
            if rho <= inj_rho:
                inj_confidence = confidence_level
                inj_confidence_area = acc_vol
                injInfo.append(inj_confidence)
                injInfo.append(inj_confidence_area)
                print('inj ' +str(vol))
                break
    return confidence_intervals, b, injInfo

def kdtree_bin(posterior,coord_names,confidence_levels,initial_boundingbox = None,samples_per_bin = 10):
    """
    takes samples and applies a KDTree to them to return confidence levels
    returns confidence_intervals - dictionary of user_provided_CL:calculated_volume
            b - ordered list of KD leaves
            initial_boundingbox - list of lists [upperleft_coords,lowerright_coords]
            injInfo - if injection values provided then returns
                      [Bounds_of_inj_kd_leaf ,number_samples_in_box, weight_of_box,injection_CL ,injection_CL_volume]
    Not quite sure that the repeated samples case is fixed, posibility of infinite loop.
    """
    confidence_levels.sort()
    print(confidence_levels)
    class Harvester(list):
        """
        when called by kdtree.operate will be used to calculate the density of each bin
        """
        def __init__(self):
            list.__init__(self)
            self.unrho=0.

        def __call__(self,tree):
            number_density=float(len(tree.objects()))/float(tree.volume())
            self.append([number_density,len(tree.objects()),tree.volume(),tree.bounds()])
            self.unrho+=number_density

        def close_ranks(self):

            for i in range(len(self)):
                self[i][0]/=self.unrho

            return sorted(self,key=itemgetter(0))

    def h(a,b):
        pass

    peparser=PEOutputParser('common')

    samples,header=posterior.samples()
    header=header.split()
    coordinatized_samples=[PosteriorSample(row, header, coord_names) for row in samples]

    #if initial bounding box is not provided, create it using max/min of sample coords.
    if initial_boundingbox is None:
        low=coordinatized_samples[0].coord()
        high=coordinatized_samples[0].coord()
        for obj in coordinatized_samples[1:]:
            low=np.minimum(low,obj.coord())
            high=np.maximum(high,obj.coord())
        initial_boundingbox = [low,high]

    tree=KDTreeVolume(coordinatized_samples,initial_boundingbox)

    a=Harvester()
    tree.operate(a,h,boxing=samples_per_bin)

    b=a.close_ranks()
    b.reverse()
    totalSamples = len(tree.objects())
    samplecounter=0.0
    for entry in b:
        samplecounter += entry[1]
        entry[1] = float(samplecounter)/float(totalSamples)

    acc_rho=0.
    acc_vol=0.
    cl_idx=0

    #check that there is an injection value for all dimension names
    def checkNone(listoParams):
        for param in listoParams:
            if posterior[param].injval is None:
                return False
        return True

    #checks for injection and extract details of the lnode in the tree that the injection is found
    if checkNone(coord_names):
        injBound,injNum,injWeight = tree.search([posterior[x].injval for x in coord_names],boxing = samples_per_bin)
        injInfo = [injBound,injNum,injWeight]
        #calculate volume of injections bin
        inj_volume = 1.
        low = injBound[0]
        high = injBound[1]
        for aCoord,bCoord in zip(low,high):
            inj_volume = inj_volume*(bCoord - aCoord)
        inj_number_density=float(injNum)/float(inj_volume)
        inj_rho = inj_number_density / a.unrho
        print(injNum,inj_volume,inj_number_density,a.unrho,injBound)
    else:
        injInfo = None
        inj_area = None
        inj_number_density=None
        inj_rho = None

    #finds the volume contained within the confidence levels requested by user
    confidence_intervals={}
    for rho,sample_number,vol,bounds in b:
        acc_vol+=vol

        if sample_number>confidence_levels[cl_idx]:
            confidence_intervals[confidence_levels[cl_idx]]=(acc_vol,sample_number)
            cl_idx+=1
            if cl_idx==len(confidence_levels):
                break

    acc_vol = 0.
    for rho,sample_number,vol,bounds in b:
        acc_vol+=vol

    #finds the confidence level of the injection and the volume of the associated contained region
    inj_confidence = None
    inj_confidence_area = None
    if inj_rho is not None:
        print('calculating cl')
        acc_vol=0.
        for rho,confidence_level,vol,bounds in b:
            acc_vol+=vol
            if rho <= inj_rho:
                inj_confidence = confidence_level
                inj_confidence_area = acc_vol
                injInfo.append(inj_confidence)
                injInfo.append(inj_confidence_area)
                break

    return confidence_intervals, b, initial_boundingbox,injInfo

def kdtree_bin2Step(posterior,coord_names,confidence_levels,initial_boundingbox = None,samples_per_bin = 10,injCoords = None,alternate = False, fraction = 0.5, skyCoords=False):
    """
    input: posterior class instance, list of confidence levels, optional choice of inital parameter space, samples per box in kdtree
    note initial_boundingbox is [[lowerbound of each param][upper bound of each param]], if not specified will just take limits of samples
    fraction is proportion of samples used for making the tree structure.
    returns: confidence_intervals, sorted list of kd objects, initial_boundingbox, injInfo
    where injInfo is [bounding box injection is found within, samples in said box, weighting of box (in case of repeated samples),inj_confidence, inj_confidence_area]
    """
    confidence_levels.sort()

    samples,header=posterior.samples()
    numberSamples = len(samples)
    if alternate == False:
        samplesStructure, samplesFill = random_split(samples,fraction)
    else:
        samplesStructure = samples[:int(numberSamples*fraction)]
        samplesFill = samples[int(numberSamples*fraction):]
    samplesFillLen = len(samplesFill)

    header=header.split()
    coordinatized_samples=[PosteriorSample(row, header, coord_names) for row in samplesStructure]
    #if initial bounding box is not provided, create it using max/min of sample coords.
    if skyCoords == True:
        initial_boundingbox = [[0,-pi_constant/2.],[2*pi_constant,pi_constant/2.]]
    if initial_boundingbox is None:
        low=coordinatized_samples[0].coord()
        high=coordinatized_samples[0].coord()
        for obj in coordinatized_samples[1:]:
            low=np.minimum(low,obj.coord())
            high=np.maximum(high,obj.coord())
        initial_boundingbox = [low,high]
    tree=KDTreeVolume(coordinatized_samples,initial_boundingbox)
    tree2fill = tree.fillNewTree(boxing=samples_per_bin, isArea = skyCoords)#set isArea True if looking at sky coords(modifies stored volume values

    columns = []
    for name in coord_names:
        columns.append(header.index(name))

    for sample in samplesFill:
        tempSample=[]
        for column in columns:
            tempSample.append(sample[column])
        addSample(tree2fill,tempSample)

    def getValues(tree,listing):
        if tree._left is None:
            listing.append([tree.bounds(),tree._importance,tree._samples,tree._volume])
        else:
            getValues(tree._left,listing)
            getValues(tree._right,listing)

    listLeaves = []
    getValues(tree2fill,listLeaves)

    clSamples = []
    for cl in confidence_levels:
        clSamples.append(samplesFillLen*cl)

    sortedLeavesList = sorted(listLeaves, key=lambda importance: importance[1])
    sortedLeavesList.reverse()
    runningTotalSamples = 0
    for i in range(len(sortedLeavesList)):
        runningTotalSamples += sortedLeavesList[i][2]
        sortedLeavesList[i].append(float(runningTotalSamples)/samplesFillLen,)


    level = 0
    countSamples = 0
    volume = 0
    lencl = len(clSamples)
    #finds confidence levels
    confidence_intervals={}
    interpConfAreas = {}
    countLeaves = 0
    for leaf in sortedLeavesList:
        countSamples += leaf[2]
        countLeaves += 1
        volume += leaf[3]
        if level < lencl and countSamples >= clSamples[level]:
            confidence_intervals[confidence_levels[level]]=(volume,float(countSamples)/samplesFillLen)
            interpConfAreas[confidence_levels[level]] = volume-leaf[3]*(countSamples-clSamples[level])/leaf[2]
            level +=1

    if injCoords is not None:
        injBound,injNum,injImportance = tree2fill.search(injCoords)
        injInfo = [injBound,injNum,injImportance]
    else:
        injInfo = None


    #finds the confidence level of the injection and the volume of the associated contained region
    inj_confidence = None
    inj_confidence_area = None
    if injInfo is not None:
        acc_vol=0.
        acc_cl=0.
        for leaf in sortedLeavesList:
            acc_vol+=leaf[3]
            acc_cl+=leaf[2]
            if leaf[1] <= injImportance:
                inj_confidence = float(acc_cl)/samplesFillLen
                inj_confidence_area = acc_vol
                injInfo.append(inj_confidence)
                injInfo.append(inj_confidence_area)
                break

    return sortedLeavesList, interpConfAreas, injInfo
    #check that there is an injection value for all dimension names
    def checkNone(listoParams):
        for param in listoParams:
            if posterior[param].injval is None:
                return False
        return True

def greedy_bin_two_param(posterior,greedy2Params,confidence_levels):
    """
    Determine the 2-parameter Bayesian Confidence Intervals using a greedy
    binning algorithm.

    @param posterior: an instance of the Posterior class.

    @param greedy2Params: a dict - {param1Name:param1binSize,param2Name:param2binSize} .

    @param confidence_levels: A list of floats of the required confidence intervals [(0-1)].
    """

    #Extract parameter names
    par1_name,par2_name=greedy2Params.keys()

    #Set posterior array columns
    par1pos=posterior[par1_name.lower()].samples
    par2pos=posterior[par2_name.lower()].samples

    #Extract bin sizes
    par1_bin=greedy2Params[par1_name]
    par2_bin=greedy2Params[par2_name]

    #Extract injection information
    par1_injvalue=posterior[par1_name.lower()].injval
    par2_injvalue=posterior[par2_name.lower()].injval

    #Create 2D bin array
    par1pos_min=min(par1pos)[0]
    par2pos_min=min(par2pos)[0]

    par1pos_max=max(par1pos)[0]
    par2pos_max=max(par2pos)[0]

    par1pos_Nbins= int(ceil((par1pos_max - par1pos_min)/par1_bin))+1

    par2pos_Nbins= int(ceil((par2pos_max - par2pos_min)/par2_bin))+1

    greedyHist = np.zeros(par1pos_Nbins*par2pos_Nbins,dtype='i8')
    greedyPoints = np.zeros((par1pos_Nbins*par2pos_Nbins,2))

    #Fill bin values
    par1_point=par1pos_min
    par2_point=par2pos_min
    for i in range(par2pos_Nbins):

        par1_point=par1pos_min
        for j in range(par1pos_Nbins):

            greedyPoints[j+par1pos_Nbins*i,0]=par1_point
            greedyPoints[j+par1pos_Nbins*i,1]=par2_point
            par1_point+=par1_bin
        par2_point+=par2_bin


    #If injection point given find which bin its in...
    injbin=None
    if par1_injvalue is not None and par2_injvalue is not None:

        par1_binNumber=int(floor((par1_injvalue-par1pos_min)/par1_bin))
        par2_binNumber=int(floor((par2_injvalue-par2pos_min)/par2_bin))

        injbin=int(par1_binNumber+par2_binNumber*par1pos_Nbins)
    elif par1_injvalue is None and par2_injvalue is not None:
        print("Injection value not found for %s!"%par1_name)

    elif par1_injvalue is not None and par2_injvalue is None:
        print("Injection value not found for %s!"%par2_name)

    #Bin posterior samples
    for par1_samp,par2_samp in zip(par1pos,par2pos):
        par1_samp=par1_samp[0]
        par2_samp=par2_samp[0]
        par1_binNumber=int(floor((par1_samp-par1pos_min)/par1_bin))
        par2_binNumber=int(floor((par2_samp-par2pos_min)/par2_bin))
        try:
            greedyHist[par1_binNumber+par2_binNumber*par1pos_Nbins]+=1
        except:
            raise RuntimeError("Problem binning samples: %i,%i,%i,%i,%i,%f,%f,%f,%f,%f,%f .")%(par1_binNumber,par2_binNumber,par1pos_Nbins,par2pos_Nbins,par1_binNumber+par2_binNumber*par1pos_Nbins,par1_samp,par1pos_min,par1_bin,par1_samp,par2pos_min,par2_bin)
    #Call greedy bins routine
    toppoints,injection_cl,reses,injection_area=\
                                _greedy_bin(
                                                greedyHist,
                                                greedyPoints,
                                                injbin,
                                                float(par1_bin*par2_bin),
                                                int(len(par1pos)),
                                                confidence_levels
                                            )

    return toppoints,injection_cl,reses,injection_area

def pol2cart(long,lat):
    """
    Utility function to convert longitude,latitude on a unit sphere to
    cartesian co-ordinates.
    """

    x=np.cos(lat)*np.cos(long)
    y=np.cos(lat)*np.sin(long)
    z=np.sin(lat)
    return np.array([x,y,z])
#

def sph2cart(r,theta,phi):
    """
    Utiltiy function to convert r,theta,phi to cartesian co-ordinates.
    """
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)
    return x,y,z


def cart2sph(x,y,z):
    """
    Utility function to convert cartesian coords to r,theta,phi.
    """
    r = np.sqrt(x*x + y*y + z*z)
    theta = np.arccos(z/r)
    phi = np.fmod(2*pi_constant + np.arctan2(y,x), 2*pi_constant)

    return r,theta,phi



def plot_sky_map(hpmap, outdir, inj=None, nest=True):
    """Plots a sky map from a healpix map, optionally including an
    injected position. This is a temporary map to display before
    ligo.skymap utility is used to generated a smoother one.

    :param hpmap: An array representing a healpix map (in nested
      ordering if ``nest = True``).

    :param outdir: The output directory.

    :param inj: If not ``None``, then ``[ra, dec]`` of the injection
      associated with the posterior map.

    :param nest: Flag indicating the pixel ordering in healpix.

    """

    fig = plt.figure(frameon=False, figsize=(8,6))
    hp.mollview(hpmap, nest=nest, min=0, max=np.max(hpmap), cmap='Greys', coord='E', fig=fig.number, title='Histogrammed skymap' )
    plt.grid(True,color='g',figure=fig)

    if inj is not None:
        theta = np.pi/2.0 - inj[1]
        hp.projplot(theta, inj[0], '*', markerfacecolor='white', markeredgecolor='black', markersize=10)

    plt.savefig(os.path.join(outdir, 'skymap.png'))

    return fig

def skymap_confidence_areas(hpmap, cls):
    """Returns the area (in square degrees) for each confidence level with
    a greedy binning algorithm for the given healpix map.

    """

    hpmap = hpmap / np.sum(hpmap) # Normalise to sum to one.

    hpmap = np.sort(hpmap)[::-1] # Sort from largest to smallest
    cum_hpmap = np.cumsum(hpmap)

    pixarea = hp.nside2pixarea(hp.npix2nside(hpmap.shape[0]))
    pixarea = pixarea*(180.0/np.pi)**2 # In square degrees

    areas = []
    for cl in cls:
        npix = np.sum(cum_hpmap < cl) # How many pixels to sum before cl?
        areas.append(npix*pixarea)

    return np.array(areas)

def skymap_inj_pvalue(hpmap, inj, nest=True):
    """Returns the greedy p-value estimate for the given injection.

    """

    nside = hp.npix2nside(hpmap.shape[0])
    hpmap = hpmap / np.sum(hpmap) # Normalise to sum to one

    injpix = hp.ang2pix(nside, np.pi/2.0-inj[1], inj[0], nest=nest)
    injvalue = hpmap[injpix]

    return np.sum(hpmap[hpmap >= injvalue])

#

def mc2ms(mc,eta):
    """
    Utility function for converting mchirp,eta to component masses. The
    masses are defined so that m1>m2. The rvalue is a tuple (m1,m2).
    """
    root = np.sqrt(0.25-eta)
    fraction = (0.5+root) / (0.5-root)
    invfraction = 1/fraction

    m2= mc * np.power((1+fraction),0.2) / np.power(fraction,0.6)

    m1= mc* np.power(1+invfraction,0.2) / np.power(invfraction,0.6)
    return (m1,m2)
#
#

def q2ms(mc,q):
    """
    Utility function for converting mchirp,q to component masses. The
    masses are defined so that m1>m2. The rvalue is a tuple (m1,m2).
    """
    factor = mc * np.power(1+q, 1.0/5.0);
    m1 = factor * np.power(q, -3.0/5.0);
    m2 = factor * np.power(q, 2.0/5.0);
    return (m1,m2)
#
#

def q2eta(q):
    """
    Utility function for converting q to eta. The
    rvalue is eta.
    """
    eta = q/((1+q)*(1+q))
    return np.clip(eta,0,0.25) # Explicitly cap eta at 0.25, in case it exceeds this slightly due to floating-point issues
#
#

def mc2q(mc,eta):
    """
    Utility function for converting mchirp,eta to new mass ratio q (m2/m1).
    """
    m1,m2 = mc2ms(mc,eta)
    q = m2/m1
    return q
#
#

def ang_dist(long1,lat1,long2,lat2):
    """
    Find the angular separation of (long1,lat1) and (long2,lat2), which are
        specified in radians.
    """

    x1=np.cos(lat1)*np.cos(long1)
    y1=np.cos(lat1)*np.sin(long1)
    z1=np.sin(lat1)
    x2=np.cos(lat2)*np.cos(long2)
    y2=np.cos(lat2)*np.sin(long2)
    z2=np.sin(lat2)
    sep=math.acos(x1*x2+y1*y2+z1*z2)
    return(sep)
#
#

def array_dot(vec1, vec2):
    """
    Calculate dot products between vectors in rows of numpy arrays.
    """
    if vec1.ndim==1:
        product = (vec1*vec2).sum()
    else:
        product = (vec1*vec2).sum(axis=1).reshape(-1,1)
    return product
#
#

def array_ang_sep(vec1, vec2):
    """
    Find angles between vectors in rows of numpy arrays.
    """
    vec1_mag = np.sqrt(array_dot(vec1, vec1))
    vec2_mag = np.sqrt(array_dot(vec2, vec2))
    return np.arccos(array_dot(vec1, vec2)/(vec1_mag*vec2_mag))
#
#

def array_polar_ang(vec):
    """
    Find polar angles of vectors in rows of a numpy array.
    """
    if vec.ndim==1:
        z = vec[2]
    else:
        z = vec[:,2].reshape(-1,1)
    norm = np.sqrt(array_dot(vec,vec))
    return np.arccos(z/norm)
#
#

def rotation_matrix(angle, direction):
    """
    Compute general rotation matrices for a given angles and direction vectors.
    """
    cosa = np.cos(angle)
    sina = np.sin(angle)
    direction /= np.sqrt(array_dot(direction,direction))
    #Assume calculating array of rotation matrices.
    try:
        nSamps = len(angle)
        R = np.array( [np.diag([i,i,i]) for i in cosa.flat] )
        R += np.array( [np.outer(direction[i],direction[i])*(1.0-cosa[i]) for i in range(nSamps)] )
        R += np.array( [np.array(   [[ 0.0,            -direction[i,2],    direction[i,1]],
                                     [ direction[i,2],  0.0,              -direction[i,0]],
                                     [-direction[i,1],  direction[i,0],    0.0          ]] ) * sina[i] for i in range(nSamps)] )
    #Only computing one rotation matrix.
    except TypeError:
        R = np.diag([cosa,cosa,cosa])
        R += np.outer(direction,direction) * (1.0 - cosa)
        R += np.array(   [[ 0.0,            -direction[2],    direction[1]],
                          [ direction[2],  0.0,              -direction[0]],
                          [-direction[1],  direction[0],    0.0          ]] ) * sina
    return R
#
#

def rotate_vector(R, vec):
    """
    Rotate vectors using the given rotation matrices.
    """
    if vec.ndim == 1:
        newVec = np.dot(R[i],vec[i])
    else:
        newVec = np.array( [np.dot(R[i],vec[i]) for i in range(len(vec))] )
    return newVec
#
#
def ROTATEZ(angle, vx, vy, vz):
    # This is the ROTATEZ in LALSimInspiral.c.
    tmp1 = vx*np.cos(angle) - vy*np.sin(angle);
    tmp2 = vx*np.sin(angle) + vy*np.cos(angle);
    return np.asarray([tmp1,tmp2,vz])

def ROTATEY(angle, vx, vy, vz):
    # This is the ROTATEY in LALSimInspiral.c
    tmp1 = vx*np.cos(angle) + vz*np.sin(angle);
    tmp2 = - vx*np.sin(angle) + vz*np.cos(angle);
    return np.asarray([tmp1,vy,tmp2])

def orbital_momentum(fref, m1,m2, inclination):
    """
    Calculate orbital angular momentum vector.
    Note: The units of Lmag are different than what used in lalsimulation.
    Mc must be called in units of Msun here.

    Note that if one wants to build J=L+S1+S2 with L returned by this function, S1 and S2
    must not get the Msun^2 factor.
    """
    eta = m1*m2/( (m1+m2)*(m1+m2) )
    Lmag = orbital_momentum_mag(fref, m1,m2,eta)
    Lx, Ly, Lz = sph2cart(Lmag, inclination, 0.0)
    return np.hstack((Lx,Ly,Lz))
#
#
def orbital_momentum_mag(fref, m1,m2,eta):
    v0 = np.power((m1+m2) *pi_constant * lal.MTSUN_SI * fref, 1.0/3.0)
    #1 PN Mtot*Mtot*eta/v
    PNFirst = (((m1+m2)**2)*eta)/v0
    PNSecond = 1+ (v0**2) * (3.0/2.0 +eta/6.0)
    Lmag= PNFirst*PNSecond
    return Lmag

def component_momentum(m, a, theta, phi):
    """
    Calculate BH angular momentum vector.
    """
    Sx, Sy, Sz = sph2cart(m**2 * a, theta, phi)
    return np.hstack((Sx,Sy,Sz))
#
#

def symm_tidal_params(lambda1,lambda2,q):
    """
    Calculate best tidal parameters [Eqs. (5) and (6) in Wade et al. PRD 89, 103012 (2014)]
    Requires q <= 1
    """
    lambdap = lambda1 + lambda2
    lambdam = lambda1 - lambda2

    # Check that q <= 1, as expected
    if np.any(q > 1):
        raise ValueError("q > 1, while this function requires q <= 1.")

    dmbym = (1. - q)/(1. + q) # Equivalent to sqrt(1 - 4*eta) for q <= 1

    eta = q2eta(q)

    lam_tilde = (8./13.)*((1.+7.*eta-31.*eta*eta)*lambdap + dmbym*(1.+9.*eta-11.*eta*eta)*lambdam)
    dlam_tilde = (1./2.)*(dmbym*(1.-13272.*eta/1319.+8944.*eta*eta/1319.)*lambdap + (1.-15910.*eta/1319.+32850.*eta*eta/1319.+3380.*eta*eta*eta/1319.)*lambdam)
    return lam_tilde, dlam_tilde

def spin_angles(fref,mc,eta,incl,a1,theta1,phi1,a2=None,theta2=None,phi2=None):
    """
    Calculate physical spin angles.
    """
    singleSpin = None in (a2,theta2,phi2)
    m1, m2 = mc2ms(mc,eta)
    L  = orbital_momentum(fref, m1,m2, incl)
    S1 = component_momentum(m1, a1, theta1, phi1)
    if not singleSpin:
        S2 = component_momentum(m2, a2, theta2, phi2)
    else:
        S2 = 0.0
    J = L + S1 + S2
    tilt1 = array_ang_sep(L,S1)
    if not singleSpin:
        tilt2 = array_ang_sep(L,S2)
    else:
        tilt2 = None
    theta_jn = array_polar_ang(J)
    beta  = array_ang_sep(J,L)
    return tilt1, tilt2, theta_jn, beta
#
def chi_precessing(m1, a1, tilt1, m2, a2, tilt2):
    """
    Calculate the magnitude of the effective precessing spin
    following convention from Phys. Rev. D 91, 024043   --   arXiv:1408.1810
    note: the paper uses naming convention where m1 < m2
    (and similar for associated spin parameters) and q > 1
    """
    q_inv = m1/m2
    A1 = 2. + (3.*q_inv/2.)
    A2 = 2. + 3./(2.*q_inv)
    S1_perp = a1*np.sin(tilt1)*m1*m1
    S2_perp = a2*np.sin(tilt2)*m2*m2
    Sp = np.maximum(A1*S2_perp, A2*S1_perp)
    chi_p = Sp/(A2*m1*m1)
    return chi_p

def calculate_redshift(distance,h=0.6790,om=0.3065,ol=0.6935,w0=-1.0):
    """
    Calculate the redshift from the luminosity distance measurement using the
    Cosmology Calculator provided in LAL.
    By default assuming cosmological parameters from arXiv:1502.01589 - 'Planck 2015 results. XIII. Cosmological parameters'
    Using parameters from table 4, column 'TT+lowP+lensing+ext'
    This corresponds to Omega_M = 0.3065, Omega_Lambda = 0.6935, H_0 = 67.90 km s^-1 Mpc^-1
    Returns an array of redshifts
    """
    def find_z_root(z,dl,omega):
        return dl - lal.LuminosityDistance(omega,z)

    omega = lal.CreateCosmologicalParameters(h,om,ol,w0,0.0,0.0)
    if isinstance(distance,float):
        z = np.array([newton(find_z_root,np.random.uniform(0.0,2.0),args = (distance,omega))])
    else:
        z = np.array([newton(find_z_root,np.random.uniform(0.0,2.0),args = (d,omega)) for d in distance[:,0]])
    return z.reshape(z.shape[0],1)

def source_mass(mass, redshift):
    """
    Calculate source mass parameter for mass m as:
    m_source = m / (1.0 + z)
    For a parameter m.
    """
    return mass / (1.0 + redshift)

## Following functions added for testing Lorentz violations
def integrand_distance(redshift,nonGR_alpha):
    """
    Calculate D_alpha integral; multiplicative factor put later
    D_alpha = integral{ ((1+z')^(alpha-2))/sqrt(Omega_m*(1+z')^3 +Omega_lambda) dz'} # eq.15 of arxiv 1110.2720
    """
    omega = lal.CreateCosmologicalParameters(0.6790,0.3065,0.6935,-1.0,0.0,0.0) ## Planck 2015 values
    omega_m = omega.om # matter density
    omega_l = omega.ol # dark energy density
    #lal.DestroyCosmologicalParameters(omega)
    return (1.0+redshift)**(nonGR_alpha-2.0)/(np.sqrt(omega_m*(1.0+redshift)**3.0 + omega_l))

def DistanceMeasure(redshift,nonGR_alpha):
    """
    D_alpha = ((1+z)^(1-alpha))/H_0 * D_alpha # from eq.15 of arxiv 1110.2720
    D_alpha calculated from integrand in above function
    """
    omega = lal.CreateCosmologicalParameters(0.6790,0.3065,0.6935,-1.0,0.0,0.0) ## Planck 2015 values
    H0 = omega.h*lal.H0FAC_SI ## Hubble constant in SI units
    dist = integrate.quad(integrand_distance, 0, redshift ,args=(nonGR_alpha))[0]
    dist *= (1.0 + redshift)**(1.0 - nonGR_alpha)
    dist /= H0
    #lal.DestroyCosmologicalParameters(omega)
    return dist*lal.C_SI ## returns D_alpha in metres

def lambda_a(redshift, nonGR_alpha, lambda_eff, distance):
    """
    Converting from the effective wavelength-like parameter to lambda_A:
    lambda_A = lambda_{eff}*(D_alpha/D_L)^(1/(2-alpha))*(1/(1+z)^((1-alpha)/(2-alpha)))
    """
    Dfunc = np.vectorize(DistanceMeasure)
    D_alpha = Dfunc(redshift, nonGR_alpha)
    dl = distance*lal.PC_SI*1e6  ## luminosity distane in metres
    return lambda_eff*(D_alpha/(dl*(1.0+redshift)**(1.0-nonGR_alpha)))**(1./(2.0-nonGR_alpha))

def amplitudeMeasure(redshift, nonGR_alpha, lambda_eff, distance):
    """
    Converting to Lorentz violating parameter "A" in dispersion  relation from lambda_A:
    A = (lambda_A/h)^(alpha-2) # eqn. 13 of arxiv 1110.2720
    """
    hPlanck = 4.13567e-15 # Planck's constant in eV.s
    ampFunc = np.vectorize(lambda_a)
    lambdaA = ampFunc(redshift, nonGR_alpha, lambda_eff, distance)/lal.C_SI # convert to seconds
    return (lambdaA/hPlanck)**(nonGR_alpha-2.0)
############################ changes for testing Lorentz violations made till here

def physical2radiationFrame(theta_jn, phi_jl, tilt1, tilt2, phi12, a1, a2, m1, m2, fref,phiref):
    """
    Wrapper function for SimInspiralTransformPrecessingNewInitialConditions().
    Vectorizes function for use in append_mapping() methods of the posterior class.
    """
    try:
        import lalsimulation as lalsim
    except ImportError:
        print('bayespputils.py: Cannot import lalsimulation SWIG bindings to calculate physical to radiation')
        print('frame angles, did you remember to use --enable-swig-python when ./configuring lalsimulation?')
        return None
    from numpy import shape
    transformFunc = lalsim.SimInspiralTransformPrecessingNewInitialConditions

    # Convert component masses to SI units
    m1_SI = m1*lal.MSUN_SI
    m2_SI = m2*lal.MSUN_SI

    # Flatten arrays
    ins = [theta_jn, phi_jl, tilt1, tilt2, phi12, a1, a2, m1_SI, m2_SI, fref,phiref]
    if len(shape(ins))>1:
        # ins is a list of lists (i.e. we are converting full posterior chains)
        try:
            for p,param in enumerate(ins):
                ins[p] = param.flatten()
        except:
            pass

        try:
            results = np.array([transformFunc(t_jn, p_jl, t1, t2, p12, a1, a2, m1_SI, m2_SI, f,phir) for (t_jn, p_jl, t1, t2, p12, a1, a2, m1_SI, m2_SI, f,phir) in zip(*ins)])
            iota = results[:,0].reshape(-1,1)
            spin1x = results[:,1].reshape(-1,1)
            spin1y = results[:,2].reshape(-1,1)
            spin1z = results[:,3].reshape(-1,1)
            spin2x = results[:,4].reshape(-1,1)
            spin2y = results[:,5].reshape(-1,1)
            spin2z = results[:,6].reshape(-1,1)
            a1,theta1,phi1 = cart2sph(spin1x,spin1y,spin1z)
            a2,theta2,phi2 = cart2sph(spin2x,spin2y,spin2z)

            mc = np.power(m1*m2,3./5.)*np.power(m1+m2,-1./5.)
            L  = orbital_momentum(fref, m1,m2, iota)
            S1 = np.hstack([m1*m1*spin1x,m1*m1*spin1y,m1*m1*spin1z])
            S2 = np.hstack([m2*m2*spin2x,m2*m2*spin2y,m2*m2*spin2z])
            J = L + S1 + S2
            beta = array_ang_sep(J,L)

            return iota, theta1, phi1, theta2, phi2, beta
        except: # Catch all exceptions, including failure for the transformFunc
            # Something went wrong, returning None
            return None

    elif len(shape(ins))<=1:
       # ins is a list of floats (i.e. we are converting the injected values) or empty
        try:
            for p,param in enumerate(ins):
                ins[p] = param
        except:
            pass

        try:
            results = np.array(transformFunc(theta_jn, phi_jl, tilt1, tilt2, phi12, a1, a2, m1_SI, m2_SI, fref,phiref))
            iota = results[0]
            spin1x = results[1]
            spin1y = results[2]
            spin1z = results[3]
            spin2x = results[4]
            spin2y = results[5]
            spin2z = results[6]
            a1,theta1,phi1 = cart2sph(spin1x,spin1y,spin1z)
            a2,theta2,phi2 = cart2sph(spin2x,spin2y,spin2z)

            mc = np.power(m1*m2,3./5.)*np.power(m1+m2,-1./5.)
            L  = orbital_momentum(fref, m1,m2, iota)
            S1 = m1*m1*np.hstack([spin1x,spin1y,spin1z])
            S2 = m2*m2*np.hstack([spin2x,spin2y,spin2z])
            J = L + S1 + S2
            beta = array_ang_sep(J,L)

            return iota, theta1, phi1, theta2, phi2, beta

        except: # Catch all exceptions, including failure for the transformFunc
            # Something went wrong, returning None
            return None
#
#

def plot_one_param_pdf_kde(fig,onedpos):

    from scipy import seterr as sp_seterr

    np.seterr(under='ignore')
    sp_seterr(under='ignore')
    pos_samps=onedpos.samples
    try:
        gkde=onedpos.gaussian_kde
    except np.linalg.linalg.LinAlgError:
        print('Singular matrix in KDE. Skipping')
    else:
        ind=np.linspace(np.min(pos_samps),np.max(pos_samps),101)
        kdepdf=gkde.evaluate(ind)
        plt.plot(ind,kdepdf,color='green')
    return

def plot_one_param_pdf_line_hist(fig,pos_samps):
    plt.hist(pos_samps,kdepdf)


def plot_one_param_pdf(posterior,plot1DParams,analyticPDF=None,analyticCDF=None,plotkde=False):
    """
    Plots a 1D histogram and (gaussian) kernel density estimate of the
    distribution of posterior samples for a given parameter.

    @param posterior: an instance of the Posterior class.

    @param plot1DParams: a dict; {paramName:Nbins}

    @param analyticPDF: an analytic probability distribution function describing the distribution.

    @param analyticCDF: an analytic cumulative distribution function describing the distribution.

    @param plotkde: Use KDE to smooth plot (default: False)
    """

    matplotlib.rcParams['text.usetex']=False

    param=list(plot1DParams.keys())[0].lower()
    histbins=list(plot1DParams.values())[0]

    pos_samps=posterior[param].samples
    injpar=posterior[param].injval
    trigvals=posterior[param].trigvals

    #myfig=plt.figure(figsize=(4,3.5),dpi=200)
    myfig=plt.figure(figsize=(6,4.5),dpi=150)
    #axes=plt.Axes(myfig,[0.2, 0.2, 0.7,0.7])
    axes=plt.Axes(myfig,[0.15,0.15,0.6,0.76])
    myfig.add_axes(axes)
    majorFormatterX=ScalarFormatter(useMathText=True)
    majorFormatterX.format_data=lambda data:'%.6g'%(data)
    majorFormatterY=ScalarFormatter(useMathText=True)
    majorFormatterY.format_data=lambda data:'%.6g'%(data)
    majorFormatterX.set_scientific(True)
    majorFormatterY.set_scientific(True)
    offset=0.0
    if param.find('time')!=-1:
        offset=floor(min(pos_samps))
        pos_samps=pos_samps-offset
        if injpar is not None:
            injpar=injpar-offset
        ax1_name=param+' + %i'%(int(offset))
    else: ax1_name=param

    (n, bins, patches)=plt.hist(pos_samps,histbins,density=True,facecolor='grey')
    Nchars=max(map(lambda d:len(majorFormatterX.format_data(d)),axes.get_xticks()))
    if Nchars>8:
        Nticks=3
    elif Nchars>5:
        Nticks=4
    elif Nchars>4:
        Nticks=6
    else:
        Nticks=6
    locatorX=matplotlib.ticker.MaxNLocator(nbins=Nticks)
    xmin,xmax=plt.xlim()
    if param=='rightascension' or param=='ra':
        locatorX=RALocator(min=xmin,max=xmax)
        majorFormatterX=RAFormatter()
    if param=='declination' or param=='dec':
        locatorX=DecLocator(min=xmin,max=xmax)
        majorFormatterX=DecFormatter()
    axes.xaxis.set_major_formatter(majorFormatterX)
    axes.yaxis.set_major_formatter(majorFormatterY)

    locatorX.view_limits(bins[0],bins[-1])
    axes.xaxis.set_major_locator(locatorX)
    if plotkde:  plot_one_param_pdf_kde(myfig,posterior[param])
    histbinSize=bins[1]-bins[0]
    if analyticPDF:
        (xmin,xmax)=plt.xlim()
        x = np.linspace(xmin,xmax,2*len(bins))
        plt.plot(x, analyticPDF(x+offset), color='r', linewidth=2, linestyle='dashed')
        if analyticCDF:
                # KS underflows with too many samples
            max_samps=1000
            from numpy.random import permutation
            samps = permutation(pos_samps)[:max_samps,:].flatten()
            D,p = stats.kstest(samps+offset, analyticCDF, mode='asymp')
            plt.title("%s: ks p-value %.3f"%(param,p))

    rbins=None

    if injpar is not None:
        # We will plot the injection if it is <5% outside the posterior
        delta_samps=max(pos_samps)-min(pos_samps)
        minrange=min(pos_samps)-0.05*delta_samps
        maxrange=max(pos_samps)+0.05*delta_samps
        if minrange<injpar and maxrange>injpar:

            plt.axvline(injpar, color='r', linestyle='-.', linewidth=4)

            #rkde=gkde.integrate_box_1d(min(pos[:,i]),getinjpar(injection,i))
            #print "r of injected value of %s (kde) = %f"%(param,rkde)

            #Find which bin the true value is in
        if min(pos_samps)<injpar and max(pos_samps)>injpar:
            bins_to_inj=(injpar-bins[0])/histbinSize
            injbinh=int(floor(bins_to_inj))
            injbin_frac=bins_to_inj-float(injbinh)

            #Integrate over the bins
            rbins=(sum(n[0:injbinh-1])+injbin_frac*n[injbinh])*histbinSize

    if trigvals is not None:
        for IFO in [IFO for IFO in trigvals.keys()]:
            trigval = trigvals[IFO]
            if min(pos_samps)<trigval and max(pos_samps)>trigval:
                if IFO=='H1': color = 'r'
                elif IFO=='L1': color = 'g'
                elif IFO=='V1': color = 'm'
                else: color = 'c'
                plt.axvline(trigval, color=color, linestyle='-.')
    #
    plt.grid()
    plt.xlabel(plot_label(ax1_name))
    plt.ylabel('Probability Density')

    # For RA and dec set custom labels and for RA reverse
    if(param.lower()=='ra' or param.lower()=='rightascension'):
        xmin,xmax=plt.xlim()
        plt.xlim(xmax,xmin)
    #if(param.lower()=='ra' or param.lower()=='rightascension'):
    #    locs, ticks = plt.xticks()
    #    newlocs, newticks = formatRATicks(locs)
    #    plt.xticks(newlocs,newticks,rotation=45)
    #if(param.lower()=='dec' or param.lower()=='declination'):
    #    locs, ticks = plt.xticks()
    #    newlocs, newticks = formatDecTicks(locs)
    #    plt.xticks(newlocs,newticks,rotation=45)

    return rbins,myfig#,rkde
#

class RALocator(matplotlib.ticker.MultipleLocator):
    """
    RA tick locations with some intelligence
    """
    def __init__(self,min=0.0,max=2.0*pi_constant):
        hour=pi_constant/12.0
        if(max-min)>12.0*hour:
            base=3.0*hour
        elif(max-min)>6.0*hour:
            base=2.0*hour
        # Put hour ticks if there are more than 3 hours displayed
        elif (max-min)>3.0*pi_constant/12.0:
            base=hour
        elif (max-min)>hour:
            base=hour/2.0
        else:
            base=hour/4.0

        matplotlib.ticker.MultipleLocator.__init__(self,base=base)

class DecLocator(matplotlib.ticker.MultipleLocator):
    """
    Dec tick locations with some intelligence
    """
    def __init__(self, min=-pi_constant/2.0,max=pi_constant/2.0):
        deg=pi_constant/180.0
        if (max-min)>60*deg:
            base=30.0*deg
        elif (max-min)>20*deg:
            base=10*deg
        elif (max-min)>10*deg:
            base=5*deg
        else:
            base=deg
        matplotlib.ticker.MultipleLocator.__init__(self,base=base)

class RAFormatter(matplotlib.ticker.FuncFormatter):
    def __init__(self,accuracy='auto'):
        matplotlib.ticker.FuncFormatter.__init__(self,getRAString)

class DecFormatter(matplotlib.ticker.FuncFormatter):
    def __init__(self,accuracy='auto'):
        matplotlib.ticker.FuncFormatter.__init__(self,getDecString)

def formatRATicks(locs, accuracy='auto'):
    """
    Format locs, ticks to RA angle with given accuracy
    accuracy can be 'hour', 'min', 'sec', 'all'
    returns (locs, ticks)
    'all' does no rounding, just formats the tick strings
    'auto' will use smallest appropriate units
    """
    newmax=max(locs)
    newmin=min(locs)
    if(accuracy=='auto'):
        acc='hour'
        if abs(newmax-newmin)<pi_constant/12.:
            acc='min'
        if abs(newmax-newmin)<pi_constant/(12.*60.):
            acc='sec'
    else:
        acc=accuracy

    if max(locs)>2*pi_constant: newmax=2.0*pi_constant
    if min(locs)<0.0: newmin=0.0
    locs=linspace(newmin,newmax,len(locs))

    roundlocs=map(lambda a: roundRadAngle(a, accuracy=acc), locs)

    newlocs=filter(lambda a:a>=0 and a<=2.0*pi_constant, roundlocs)
    return (list(newlocs), map(getRAString, list(newlocs) ) )

def formatDecTicks(locs, accuracy='auto'):
    """
    Format locs to Dec angle with given accuracy
    accuracy can be 'deg', 'arcmin', 'arcsec', 'all'
    'all' does no rounding, just formats the tick strings
    """
    newmax=max(locs)
    newmin=min(locs)
    if(accuracy=='auto'):
        acc='deg'
        if abs(newmax-newmin)<pi_constant/360.:
            acc='arcmin'
        if abs(newmax-newmin)<pi_constant/(360.*60.):
            acc='arcsec'
    else:
        acc=accuracy
    if newmax>0.5*pi_constant: newmax=0.5*pi_constant
    if newmin<-0.5*pi_constant: newmin=-0.5*pi_constant
    locs=linspace(newmin,newmax,len(locs))

    roundlocs=map(lambda a: roundRadAngle(a, accuracy=acc), locs)
    newlocs=filter(lambda a:a>=-pi_constant/2.0 and a<=pi_constant/2.0, roundlocs)
    return (list(newlocs), map(getDecString, list(newlocs) ) )

def roundRadAngle(rads,accuracy='all'):
    """
    round given angle in radians to integer hours, degrees, mins or secs
    accuracy can be 'hour'. 'deg', 'min', 'sec', 'all', all does nothing
    'arcmin', 'arcsec'
    """
    if accuracy=='all': return locs
    if accuracy=='hour': mult=24
    if accuracy=='deg': mult=360
    if accuracy=='min': mult=24*60
    if accuracy=='sec': mult=24*60*60
    if accuracy=='arcmin': mult=360*60
    if accuracy=='arcsec': mult=360*60*60
    mult=mult/(2.0*pi_constant)
    return round(rads*mult)/mult

def getRAString(radians,accuracy='auto'):
    secs=radians*12.0*3600/pi_constant
    hours, rem = divmod(secs, 3600 )
    mins,rem = divmod(rem, 60 )
    secs = rem
    if secs>=59.5:
        secs=secs-60
        mins=mins+1
    if mins>=59.5:
        mins=mins-60
        hours=hours+1
    if accuracy=='hour': return six.u(r'%ih'%(hours))
    if accuracy=='min': return six.u(r'%ih%im'%(hours,mins))
    if accuracy=='sec': return six.u(r'%ih%im%2.0fs'%(hours,mins,secs))
    else:
        if abs(fmod(secs,60.0))>=0.5: return(getRAString(radians,accuracy='sec'))
        if abs(fmod(mins,60.0))>=0.5: return(getRAString(radians,accuracy='min'))
        else: return(getRAString(radians,accuracy='hour'))

def getDecString(radians,accuracy='auto'):
    # LaTeX doesn't like unicode degree symbols etc
    if matplotlib.rcParams['text.usetex']:
        degsymb='$^\circ$'
        minsymb="'"
        secsymb="''"
    else:
        degsymb=six.unichr(0x0B0)
        minsymb=six.unichr(0x027)
        secsymb=six.unichr(0x2033)
    if(radians<0):
        radians=-radians
        sign=-1
    else: sign=+1
    deg,rem=divmod(radians,pi_constant/180.0)
    mins, rem = divmod(rem, pi_constant/(180.0*60.0))
    secs = rem * (180.0*3600.0)/pi_constant
    #if secs>=59.5:
    #    secs=secs-60.0
    #    mins=mins+1
    #if mins>=59.5:
    #    mins=mins-60.0
    #    deg=deg+1
    if (accuracy=='arcmin' or accuracy=='deg') and secs>30: mins=mins+1
    if accuracy=='deg' and mins>30: deg=deg+1
    if accuracy=='deg': return '%i'%(sign*deg)+degsymb
    if accuracy=='arcmin': return '%i%s%i%s'%(sign*deg,degsymb,mins,minsymb)
    if accuracy=='arcsec': return '%i%s%i%s%2.0f%s'%(sign*deg,degsymb,mins,minsymb,secs,secsymb)
    else:
    #    if abs(fmod(secs,60.0))>=0.5 and abs(fmod(secs,60)-60)>=0.5 : return(getDecString(sign*radians,accuracy='arcsec'))
    #    if abs(fmod(mins,60.0))>=0.5 and abs(fmod(mins,60)-60)>=0.5: return(getDecString(sign*radians,accuracy='arcmin'))
    #    else: return(getDecString(sign*radians,accuracy='deg'))
        return(getDecString(sign*radians,accuracy='deg'))

def plot_corner(posterior,levels,parnames=None):
    """
    Make a corner plot using the triangle module
    (See http://github.com/dfm/corner.py)
    @param posterior: The Posterior object
    @param levels: a list of confidence levels
    @param parnames: list of parameters to include
    """
    try:
        import corner
    except ImportError:
        try:
            import triangle as corner
        except ImportError:
            print('Cannot load corner module. Try running\n\t$ pip install corner')
            return None
    parnames=list(filter(lambda x: x in posterior.names, parnames))
    labels = [plot_label(parname) for parname in parnames]
    data = np.hstack([posterior[p].samples for p in parnames])
    if posterior.injection:
        injvals=[posterior[p].injval for p in parnames]
        myfig=corner.corner(data,labels=labels,truths=injvals,quantiles=levels,plot_datapoints=False,bins=20)
    else:
        myfig=corner.corner(data,labels=labels,quantiles=levels,plot_datapoints=False,bins=20)
    return(myfig)


def plot_two_param_kde_greedy_levels(posteriors_by_name,plot2DkdeParams,levels,colors_by_name,line_styles=__default_line_styles,figsize=(4,3),dpi=250,figposition=[0.2,0.2,0.48,0.75],legend='right',hatches_by_name=None,Npixels=50):
    """
    Plots a 2D kernel density estimate of the 2-parameter marginal posterior.

    @param posteriors_by_name: dictionary of Posterior objects, indexed by name

    @param plot2DkdeParams: a dict {param1Name:Nparam1Bins,param2Name:Nparam2Bins}

    @param levels: list of credible levels

    @param colors_by_name: dict of colors, indexed by name

    @param line_styles: list of line styles for the credible levels

    @param figsize: figsize to pass to matplotlib

    @param dpi: dpi to pass to matplotlib

    @param figposition: figposition to pass to matplotlib

    @param legend: position of legend

    @param hatches_by_name: dict of hatch styles indexed by name

    @param Npixels: Number of pixels in each axis of the 2D grid
    """

    from scipy import seterr as sp_seterr
    confidence_levels=levels

    # Reversed parameter order here for consistency with the other
    # plotting functions
    par2_name,par1_name=plot2DkdeParams.keys()
    xbin=plot2DkdeParams[par1_name]
    ybin=plot2DkdeParams[par2_name]
    levels= levels
    np.seterr(under='ignore')
    sp_seterr(under='ignore')

    fig=plt.figure(1,figsize=figsize,dpi=dpi)
    plt.clf()
    axes=fig.add_axes(figposition)
    name_list=[]

    #This fixes the precedence of line styles in the plot
    if len(line_styles)<len(levels):
        raise RuntimeError("Error: Need as many or more line styles to choose from as confidence levels to plot!")

    CSlst=[]
    for name,posterior in posteriors_by_name.items():
        print('Plotting '+name)
        name_list.append(name)
        par1_injvalue=posterior[par1_name].injval
        par2_injvalue=posterior[par2_name].injval

        par_trigvalues1=posterior[par1_name].trigvals
        par_trigvalues2=posterior[par2_name].trigvals
        xdat=posterior[par1_name].samples
        ydat=posterior[par2_name].samples
        a=np.squeeze(posterior[par1_name].samples)
        b=np.squeeze(posterior[par2_name].samples)
        offset=0.0
        if par1_name.find('time')!=-1:
            offset=floor(min(a))
            a=a-offset
        if par1_injvalue:
            par1_injvalue=par1_injvalue-offset
            ax1_name=par1_name+' + %i'%(int(offset))
        else: ax1_name=par1_name

        if par2_name.find('time')!=-1:
            offset=floor(min(b))
            b=b-offset
        if par2_injvalue:
            par2_injvalue=par2_injvalue-offset
            ax2_name=par2_name+' + %i'%(int(offset))
        else: ax2_name=par2_name

        samp=np.transpose(np.column_stack((xdat,ydat)))

        try:
            kde=stats.kde.gaussian_kde(samp)
            den=kde(samp)
        except:
            return None

        #grid_coords = np.append(x.reshape(-1,1),y.reshape(-1,1),axis=1)
        Nx=Npixels
        Ny=Npixels
        # Ugly hack to produce plots centred on the main mode:
        # Choose 1%-99% percentile range of the samples
        # Then zoom out by 0.95 to encompass the contour edges
        xsorted=np.sort(xdat,axis=None)
        ysorted=np.sort(ydat,axis=None)
        imin=int(0.01*len(xsorted))
        imax=int(0.99*len(xsorted))
        xmax=xsorted[imax]
        xmin=xsorted[imin]
        ymax=ysorted[imax]
        ymin=ysorted[imin]
        xmid=0.5*(xmin+xmax)
        ymid=0.5*(ymin+ymax)
        xmin=xmid+(xmin-xmid)/0.95
        xmax=xmid+(xmax-xmid)/0.95
        ymin=ymid+(ymin-ymid)/0.95
        ymax=ymid+(ymax-ymid)/0.95
        xax=np.linspace(xmin,xmax,Nx)
        yax=np.linspace(ymin,ymax,Ny)
        #print 'Plot limits %f-%f x %f-%f'%(xmin,xmax,ymin,ymax)

        # Original way which includes all samples
        #xax = np.linspace(np.min(xdat),np.max(xdat),Nx)
        #yax = np.linspace(np.min(ydat),np.max(ydat),Ny)

        # Do the actual plotting
        x,y = np.meshgrid(xax,yax)
        grid_coords = np.row_stack( (x.flatten(),y.flatten()) )
        z = kde(grid_coords)
        z = z.reshape(Nx,Ny)
        densort=np.sort(den)[::-1]
        values=[]
        Npts=xdat.shape[0]
        zvalues=[]
        for level in levels:
            ilevel = int(Npts*level + 0.5)
            if ilevel >= Npts:
                ilevel = Npts-1
            zvalues.append(densort[ilevel])
        CS=plt.contour(x, y, z, zvalues,colors=[colors_by_name[name]],linestyles=line_styles )
        CSlst.append(CS)

        if par1_injvalue is not None and par2_injvalue is not None:
            plt.plot([par1_injvalue],[par2_injvalue],'b*',scalex=False,scaley=False,markersize=12)

        if par_trigvalues1 is not None and par_trigvalues2 is not None:
            par1IFOs = set([IFO for IFO in par_trigvalues1.keys()])
            par2IFOs = set([IFO for IFO in par_trigvalues2.keys()])
            IFOs = par1IFOs.intersection(par2IFOs)
            for IFO in IFOs:
                if IFO=='H1': color = 'r'
                elif IFO=='L1': color = 'g'
                elif IFO=='V1': color = 'm'
                else: color = 'c'
                plt.plot([par_trigvalues1[IFO]],[par_trigvalues2[IFO]],color=color,marker='o',scalex=False,scaley=False)

    plt.xlabel(plot_label(par1_name))
    plt.ylabel(plot_label(par2_name))
    plt.grid()

    if len(name_list)!=len(CSlst):
        raise RuntimeError("Error number of contour objects does not equal number of names! Use only *one* contour from each set to associate a name.")

    full_name_list=[]
    dummy_lines=[]
    for plot_name in name_list:
        full_name_list.append(plot_name)
    if len(confidence_levels)>1:
        for ls_,cl in zip(line_styles[0:len(confidence_levels)],confidence_levels):
            dummy_lines.append(mpl_lines.Line2D(np.array([0.,1.]),np.array([0.,1.]),ls=ls_,color='k'))
            full_name_list.append('%s%%'%str(int(cl*100)))

    fig_actor_lst = [cs.collections[0] for cs in CSlst]
    fig_actor_lst.extend(dummy_lines)
    if legend is not None:
        try:
            twodcontour_legend=plt.figlegend(tuple(fig_actor_lst), tuple(full_name_list), loc='right',framealpha=0.1)
        except:
            twodcontour_legend=plt.figlegend(tuple(fig_actor_lst), tuple(full_name_list), loc='right')
        for text in twodcontour_legend.get_texts():
            text.set_fontsize('small')

    majorFormatterX=ScalarFormatter(useMathText=True)
    majorFormatterX.format_data=lambda data:'%.4g'%(data)
    majorFormatterY=ScalarFormatter(useMathText=True)
    majorFormatterY.format_data=lambda data:'%.4g'%(data)
    majorFormatterX.set_scientific(True)
    majorFormatterY.set_scientific(True)
    axes.xaxis.set_major_formatter(majorFormatterX)
    axes.yaxis.set_major_formatter(majorFormatterY)
    Nchars=max(map(lambda d:len(majorFormatterX.format_data(d)),axes.get_xticks()))
    if Nchars>8:
        Nticks=3
    elif Nchars>5:
        Nticks=4
    elif Nchars>4:
        Nticks=5
    else:
        Nticks=6
    locatorX=matplotlib.ticker.MaxNLocator(nbins=Nticks-1)
    if par1_name=='rightascension' or par1_name=='ra':
        (ramin,ramax)=plt.xlim()
        locatorX=RALocator(min=ramin,max=ramax)
        majorFormatterX=RAFormatter()
    if par1_name=='declination' or par1_name=='dec':
        (decmin,decmax)=plt.xlim()
        locatorX=DecLocator(min=decmin,max=decmax)
        majorFormatterX=DecFormatter()
    axes.xaxis.set_major_formatter(majorFormatterX)
    if par2_name=='rightascension' or par2_name=='ra':
        (ramin,ramax)=plt.ylim()
        locatorY=RALocator(ramin,ramax)
        axes.yaxis.set_major_locator(locatorY)
        majorFormatterY=RAFormatter()
    if par2_name=='declination' or par2_name=='dec':
        (decmin,decmax)=plt.ylim()
        locatorY=DecLocator(min=decmin,max=decmax)
        majorFormatterY=DecFormatter()
        axes.yaxis.set_major_locator(locatorY)

    axes.yaxis.set_major_formatter(majorFormatterY)
    #locatorX.view_limits(bins[0],bins[-1])
    axes.xaxis.set_major_locator(locatorX)

    fix_axis_names(plt,par1_name,par2_name)

    if(par1_name.lower()=='ra' or par1_name.lower()=='rightascension'):
        xmin,xmax=plt.xlim()
        if(xmin<0.0): xmin=0.0
        if(xmax>2.0*pi_constant): xmax=2.0*pi_constant
        plt.xlim(xmax,xmin)

    return fig


def plot_two_param_kde(posterior,plot2DkdeParams):
    """
    Plots a 2D kernel density estimate of the 2-parameter marginal posterior.

    @param posterior: an instance of the Posterior class.

    @param plot2DkdeParams: a dict {param1Name:Nparam1Bins,param2Name:Nparam2Bins}
    """

    from scipy import seterr as sp_seterr

    par1_name,par2_name=plot2DkdeParams.keys()
    Nx=plot2DkdeParams[par1_name]
    Ny=plot2DkdeParams[par2_name]

    xdat=posterior[par1_name].samples
    ydat=posterior[par2_name].samples

    par_injvalue1=posterior[par1_name].injval
    par_injvalue2=posterior[par2_name].injval

    par_trigvalues1=posterior[par1_name].trigvals
    par_trigvalues2=posterior[par2_name].trigvals

    np.seterr(under='ignore')
    sp_seterr(under='ignore')

    myfig=plt.figure(1,figsize=(6,4),dpi=200)
    myfig.add_axes(plt.Axes(myfig,[0.20,0.20,0.75,0.7]))
    plt.clf()

    xax=np.linspace(min(xdat),max(xdat),Nx)
    yax=np.linspace(min(ydat),max(ydat),Ny)
    x,y=np.meshgrid(xax,yax)

    samp=np.transpose(np.column_stack((xdat,ydat)))

    kde=stats.kde.gaussian_kde(samp)

    grid_coords = np.append(x.reshape(-1,1),y.reshape(-1,1),axis=1)

    z = kde(grid_coords.T)
    z = z.reshape(Nx,Ny)

    values=[]
    for level in levels:
        ilevel = int(Npts*level + 0.5)
        if ilevel >= Npts:
            ilevel = Npts-1
            zvalues.append(densort[ilevel])

            pp.contour(XS, YS, ZS, zvalues)

    asp=xax.ptp()/yax.ptp()
#    if(asp<0.8 or asp > 1.6): asp=1.4
    plt.imshow(z,extent=(xax[0],xax[-1],yax[0],yax[-1]),aspect=asp,origin='lower')
    plt.colorbar()

    if par_injvalue1 is not None and par_injvalue2 is not None:
        plt.plot([par_injvalue1],[par_injvalue2],'bo',scalex=False,scaley=False)

    if par_trigvalues1 is not None and par_trigvalues2 is not None:
        par1IFOs = set([IFO for IFO in par_trigvalues1.keys()])
        par2IFOs = set([IFO for IFO in par_trigvalues2.keys()])
        IFOs = par1IFOs.intersection(par2IFOs)
        for IFO in IFOs:
            if IFO=='H1': color = 'r'
            elif IFO=='L1': color = 'g'
            elif IFO=='V1': color = 'm'
            else: color = 'c'
            plt.plot([par_trigvalues1[IFO]],[par_trigvalues2[IFO]],color=color,marker='o',scalex=False,scaley=False)

    plt.xlabel(plot_label(par1_name))
    plt.ylabel(plot_label(par2_name))
    plt.grid()

    # For RA and dec set custom labels and for RA reverse
    if(par1_name.lower()=='ra' or par1_name.lower()=='rightascension'):
        xmin,xmax=plt.xlim()
        plt.xlim(xmax,xmin)
    #if(par1_name.lower()=='ra' or par1_name.lower()=='rightascension'):
    #        locs, ticks = plt.xticks()
    #        (newlocs, newticks)=formatRATicks(locs, ticks)
    #        plt.xticks(newlocs,newticks,rotation=45)
    #if(par1_name.lower()=='dec' or par1_name.lower()=='declination'):
    #        locs, ticks = plt.xticks()
    #        newlocs, newticks = formatDecTicks(locs)
    #        plt.xticks(newlocs,newticks,rotation=45)

    #if(par2_name.lower()=='ra' or par2_name.lower()=='rightascension'):
    #    ymin,ymax=plt.ylim()
    #    plt.ylim(ymax,ymin)
    #if(par2_name.lower()=='ra' or par2_name.lower()=='rightascension'):
    #    locs, ticks = plt.yticks()
    #    newlocs,newticks=formatRATicks(locs)
    #    plt.yticks(newlocs,newticks)
    #if(par2_name.lower()=='dec' or par2_name.lower()=='declination'):
    #    locs, ticks = plt.yticks()
    #    newlocs,newticks=formatDecTicks(locs)
    #    plt.yticks(newlocs,newticks)

    return myfig
#

def get_inj_by_time(injections,time):
    """
    Filter injections to find the injection with end time given by time +/- 0.1s
    """
    import itertools
    injection = itertools.ifilter(lambda a: abs(float(a.get_end()) - time) < 0.1, injections).next()
    return injection

def histogram2D(posterior,greedy2Params,confidence_levels):
    """
    Returns a 2D histogram and edges for the two parameters passed in greedy2Params, plus the actual discrete confidence levels
    imposed by the finite number of samples.
       H,xedges,yedges,Hlasts = histogram2D(posterior,greedy2Params,confidence_levels)
    @param posterior: Posterior instance
    @param greedy2Params: a dict ;{param1Name:param1binSize,param2Name:param2binSize}
    @param confidence_levels: a list of the required confidence levels to plot on the contour map.
    """

    par1_name,par2_name=greedy2Params.keys()
    par1_bin=greedy2Params[par1_name]
    par2_bin=greedy2Params[par2_name]
    par1_injvalue=posterior[par1_name.lower()].injval
    par2_injvalue=posterior[par2_name.lower()].injval
    a=np.squeeze(posterior[par1_name].samples)
    b=np.squeeze(posterior[par2_name].samples)
    par1pos_min=a.min()
    par2pos_min=b.min()
    par1pos_max=a.max()
    par2pos_max=b.max()
    par1pos_Nbins= int(ceil((par1pos_max - par1pos_min)/par1_bin))+1
    par2pos_Nbins= int(ceil((par2pos_max - par2pos_min)/par2_bin))+1
    H, xedges, yedges = np.histogram2d(a,b, bins=(par1pos_Nbins, par2pos_Nbins),density=True)
    temp=np.copy(H)
    temp=temp.ravel()
    confidence_levels.sort()
    Hsum=0
    Hlasts=[]
    idxes=np.argsort(temp)
    j=len(idxes)-1
    for cl in confidence_levels:
        while float(Hsum/np.sum(H))<cl:
        #ind = np.argsort(temp)
            max_i=idxes[j]
            j-=1
            val = temp[max_i]
            Hlast=val
            Hsum+=val
            temp[max_i]=0
        Hlasts.append(Hlast)
    return (H,xedges,yedges,Hlasts)

def plot_two_param_greedy_bins_contourf(posteriors_by_name,greedy2Params,confidence_levels,colors_by_name,figsize=(7,6),dpi=120,figposition=[0.3,0.3,0.5,0.5],legend='right',hatches_by_name=None):
    """
    @param posteriors_by_name A dictionary of posterior objects indexed by name
    @param greedy2Params: a dict ;{param1Name:param1binSize,param2Name:param2binSize}
    @param confidence_levels: a list of the required confidence levels to plot on the contour map.
    @param colors_by_name: dict of colors, indexed by name
    @param figsize: figsize to pass to matplotlib
    @param dpi: dpi to pass to matplotlib
    @param figposition: figposition to pass to matplotlib
    @param legend: Legend position to pass to matplotlib
    @param hatches_by_name: dict of hatch styles, indexed by name
    """
    fig=plt.figure(1,figsize=figsize,dpi=dpi)
    plt.clf()
    fig.add_axes(figposition)
    CSlst=[]
    name_list=[]
    par1_name,par2_name=greedy2Params.keys()
    for name,posterior in posteriors_by_name.items():
        name_list.append(name)
        H,xedges,yedges,Hlasts=histogram2D(posterior,greedy2Params,confidence_levels+[0.99999999999999])
        extent= [xedges[0], yedges[-1], xedges[-1], xedges[0]]
        CS2=plt.contourf(yedges[:-1],xedges[:-1],H,Hlasts,extend='max',colors=[colors_by_name[name]] ,alpha=0.3 )
        CS=plt.contour(yedges[:-1],xedges[:-1],H,Hlasts,extend='max',colors=[colors_by_name[name]] )
        CSlst.append(CS)

    plt.title("%s-%s confidence contours (greedy binning)"%(par1_name,par2_name)) # add a title
    plt.xlabel(plot_label(par2_name))
    plt.ylabel(plot_label(par1_name))
    if len(name_list)!=len(CSlst):
        raise RuntimeError("Error number of contour objects does not equal number of names! Use only *one* contour from each set to associate a name.")
    full_name_list=[]
    dummy_lines=[]
    for plot_name in name_list:
        full_name_list.append(plot_name)
        if len(confidence_levels)>1:
            for cl in confidence_levels+[1]:
                dummy_lines.append(mpl_lines.Line2D(np.array([0.,1.]),np.array([0.,1.]),color='k'))
                full_name_list.append('%s%%'%str(int(cl*100)))
        fig_actor_lst = [cs.collections[0] for cs in CSlst]
        fig_actor_lst.extend(dummy_lines)
    if legend is not None:
        try:
            twodcontour_legend=plt.figlegend(tuple(fig_actor_lst), tuple(full_name_list), loc='right',framealpha=0.1)
        except:
            twodcontour_legend=plt.figlegend(tuple(fig_actor_lst), tuple(full_name_list), loc='right')
        for text in twodcontour_legend.get_texts():
            text.set_fontsize('small')
    fix_axis_names(plt,par1_name,par2_name)
    return fig

def fix_axis_names(plt,par1_name,par2_name):
    """
    Fixes names of axes
    """
    return

    # For ra and dec set custom labels and for RA reverse
    if(par1_name.lower()=='ra' or par1_name.lower()=='rightascension'):
        ymin,ymax=plt.ylim()
        if(ymin<0.0): ylim=0.0
        if(ymax>2.0*pi_constant): ymax=2.0*pi_constant
        plt.ylim(ymax,ymin)
    if(par1_name.lower()=='ra' or par1_name.lower()=='rightascension'):
        locs, ticks = plt.yticks()
        newlocs, newticks = formatRATicks(locs)
        plt.yticks(newlocs,newticks)
    if(par1_name.lower()=='dec' or par1_name.lower()=='declination'):
        locs, ticks = plt.yticks()
        newlocs,newticks=formatDecTicks(locs)
        plt.yticks(newlocs,newticks)

    if(par2_name.lower()=='ra' or par2_name.lower()=='rightascension'):
        xmin,xmax=plt.xlim()
        if(xmin<0.0): xmin=0.0
        if(xmax>2.0*pi_constant): xmax=2.0*pi_constant
        plt.xlim(xmax,xmin)
    if(par2_name.lower()=='ra' or par2_name.lower()=='rightascension'):
        locs, ticks = plt.xticks()
        newlocs, newticks = formatRATicks(locs)
        plt.xticks(newlocs,newticks,rotation=45)
    if(par2_name.lower()=='dec' or par2_name.lower()=='declination'):
        locs, ticks = plt.xticks()
        newlocs, newticks = formatDecTicks(locs)
        plt.xticks(newlocs,newticks,rotation=45)
    return plt

def plot_two_param_greedy_bins_contour(posteriors_by_name,greedy2Params,confidence_levels,colors_by_name,line_styles=__default_line_styles,figsize=(4,3),dpi=250,figposition=[0.2,0.2,0.48,0.75],legend='right'):
    """
    Plots the confidence level contours as determined by the 2-parameter
    greedy binning algorithm.

    @param posteriors_by_name: A dict containing Posterior instances referenced by some id.

    @param greedy2Params: a dict ;{param1Name:param1binSize,param2Name:param2binSize}

    @param confidence_levels: a list of the required confidence levels to plot on the contour map.

    @param colors_by_name: A dict of colors cross-referenced to the above Posterior ids.

    @param legend: Argument for legend placement or None for no legend ('right', 'upper left', 'center' etc)

    @param line_styles: list of line styles for credible regions

    @param figsize: figsize to pass to matplotlib

    @param dpi: dpi to pass to matplotlib

    @param figposition: figposition to pass to matplotlib

    """

    fig=plt.figure(1,figsize=figsize,dpi=dpi)
    plt.clf()

    axes=fig.add_axes(figposition)

    #This fixes the precedence of line styles in the plot
    if len(line_styles)<len(confidence_levels):
        raise RuntimeError("Error: Need as many or more line styles to choose from as confidence levels to plot!")

    CSlst=[]
    name_list=[]
    for name,posterior in posteriors_by_name.items():

        name_list.append(name)
        #Extract parameter names
        par1_name,par2_name=greedy2Params.keys()
        #Extract bin sizes
        par1_bin=greedy2Params[par1_name]
        par2_bin=greedy2Params[par2_name]

        #Extract injection information
        par1_injvalue=posterior[par1_name.lower()].injval
        par2_injvalue=posterior[par2_name.lower()].injval

        #Extract trigger information
        par1_trigvalues=posterior[par1_name.lower()].trigvals
        par2_trigvalues=posterior[par2_name.lower()].trigvals

        a=np.squeeze(posterior[par1_name].samples)
        b=np.squeeze(posterior[par2_name].samples)

        #Create 2D bin array
        par1pos_min=a.min()
        par2pos_min=b.min()

        par1pos_max=a.max()
        par2pos_max=b.max()

        par1pos_Nbins= int(ceil((par1pos_max - par1pos_min)/par1_bin))+1
        par2pos_Nbins= int(ceil((par2pos_max - par2pos_min)/par2_bin))+1

        if par1_name.find('time')!=-1:
            offset=floor(min(a))
            a=a-offset
            if par1_injvalue:
                par1_injvalue=par1_injvalue-offset
            ax1_name=par1_name+' + %i'%(int(offset))
        else: ax1_name=par1_name

        if par2_name.find('time')!=-1:
            offset=floor(min(b))
            b=b-offset
            if par2_injvalue:
                par2_injvalue=par2_injvalue-offset
            ax2_name=par2_name+' + %i'%(int(offset))
        else: ax2_name=par2_name


        majorFormatterX=ScalarFormatter(useMathText=True)
        majorFormatterX.format_data=lambda data:'%.4g'%(data)
        majorFormatterY=ScalarFormatter(useMathText=True)
        majorFormatterY.format_data=lambda data:'%.4g'%(data)
        majorFormatterX.set_scientific(True)
        majorFormatterY.set_scientific(True)
        axes.xaxis.set_major_formatter(majorFormatterX)
        axes.yaxis.set_major_formatter(majorFormatterY)

        H, xedges, yedges = np.histogram2d(a,b, bins=(par1pos_Nbins, par2pos_Nbins),density=True)

        extent = [xedges[0], yedges[-1], xedges[-1], xedges[0]]

        temp=np.copy(H)
        temp=temp.ravel()
        confidence_levels.sort()
        Hsum=0
        Hlasts=[]
        idxes=np.argsort(temp)
        j=len(idxes)-1
        for cl in confidence_levels:
            while float(Hsum/np.sum(H))<cl:
                #ind = np.argsort(temp)
                max_i=idxes[j]
                j-=1
                val = temp[max_i]
                Hlast=val
                Hsum+=val
                temp[max_i]=0
            Hlasts.append(Hlast)
        CS=plt.contour(yedges[:-1],xedges[:-1],H,Hlasts,colors=[colors_by_name[name]],linestyles=line_styles)
        plt.grid()
        if(par1_injvalue is not None and par2_injvalue is not None):
            plt.plot([par2_injvalue],[par1_injvalue],'b*',scalex=False,scaley=False,markersize=12)
        if(par1_trigvalues is not None and par2_trigvalues is not None):
            par1IFOs = set([IFO for IFO in par1_trigvalues.keys()])
            par2IFOs = set([IFO for IFO in par2_trigvalues.keys()])
            IFOs = par1IFOs.intersection(par2IFOs)
            if IFO=='H1': color = 'r'
            elif IFO=='L1': color = 'g'
            elif IFO=='V1': color = 'm'
            else: color = 'c'
            plt.plot([par2_trigvalues[IFO]],[par1_trigvalues[IFO]],color=color,marker='*',scalex=False,scaley=False)
        CSlst.append(CS)

        Nchars=max(map(lambda d:len(majorFormatterX.format_data(d)),axes.get_xticks()))
        if Nchars>8:
            Nticks=3
        elif Nchars>5:
            Nticks=4
        elif Nchars>4:
            Nticks=5
        else:
            Nticks=6
        locatorX=matplotlib.ticker.MaxNLocator(nbins=Nticks-1)
        if par2_name=='rightascension' or par2_name=='ra':
            (ramin,ramax)=plt.xlim()
            locatorX=RALocator(min=ramin,max=ramax)
            majorFormatterX=RAFormatter()
        if par2_name=='declination' or par2_name=='dec':
            (decmin,decmax)=plt.xlim()
            locatorX=DecLocator(min=decmin,max=decmax)
            majorFormatterX=DecFormatter()
        axes.xaxis.set_major_formatter(majorFormatterX)
        if par1_name=='rightascension' or par1_name=='ra':
            (ramin,ramax)=plt.ylim()
            locatorY=RALocator(ramin,ramax)
            axes.yaxis.set_major_locator(locatorY)
            majorFormatterY=RAFormatter()
        if par1_name=='declination' or par1_name=='dec':
            (decmin,decmax)=plt.ylim()
            locatorY=DecLocator(min=decmin,max=decmax)
            majorFormatterY=DecFormatter()
            axes.yaxis.set_major_locator(locatorY)

        axes.yaxis.set_major_formatter(majorFormatterY)
        #locatorX.view_limits(bins[0],bins[-1])
        axes.xaxis.set_major_locator(locatorX)

    #plt.title("%s-%s confidence contours (greedy binning)"%(par1_name,par2_name)) # add a title
    plt.xlabel(plot_label(ax2_name))
    plt.ylabel(plot_label(ax1_name))

    if len(name_list)!=len(CSlst):
        raise RuntimeError("Error number of contour objects does not equal number of names! Use only *one* contour from each set to associate a name.")
    full_name_list=[]
    dummy_lines=[]

    for plot_name in name_list:
        full_name_list.append(plot_name)
    if len(confidence_levels)>1:
        for ls_,cl in zip(line_styles[0:len(confidence_levels)],confidence_levels):
            dummy_lines.append(mpl_lines.Line2D(np.array([0.,1.]),np.array([0.,1.]),ls=ls_,color='k'))
            full_name_list.append('%s%%'%str(int(cl*100)))

    fig_actor_lst = [cs.collections[0] for cs in CSlst]

    fig_actor_lst.extend(dummy_lines)

    if legend is not None:
        try:
            twodcontour_legend=plt.figlegend(tuple(fig_actor_lst), tuple(full_name_list), loc='right',framealpha=0.1)
        except:
            twodcontour_legend=plt.figlegend(tuple(fig_actor_lst), tuple(full_name_list), loc='right')
        for text in twodcontour_legend.get_texts():
            text.set_fontsize('small')


    # For ra and dec set custom labels and for RA reverse
    #if(par1_name.lower()=='ra' or par1_name.lower()=='rightascension'):
    #        ymin,ymax=plt.ylim()
    #        if(ymin<0.0): ylim=0.0
    #        if(ymax>2.0*pi_constant): ymax=2.0*pi_constant
    #        plt.ylim(ymax,ymin)
    #if(par1_name.lower()=='ra' or par1_name.lower()=='rightascension'):
    #        locs, ticks = plt.yticks()
    #        newlocs, newticks = formatRATicks(locs)
    #        plt.yticks(newlocs,newticks)
    #if(par1_name.lower()=='dec' or par1_name.lower()=='declination'):
    #        locs, ticks = plt.yticks()
    #        newlocs,newticks=formatDecTicks(locs)
    #        plt.yticks(newlocs,newticks)

    if(par2_name.lower()=='ra' or par2_name.lower()=='rightascension'):
        xmin,xmax=plt.xlim()
        if(xmin<0.0): xmin=0.0
        if(xmax>2.0*pi_constant): xmax=2.0*pi_constant
        plt.xlim(xmax,xmin)
    #if(par2_name.lower()=='ra' or par2_name.lower()=='rightascension'):
    #    locs, ticks = plt.xticks()
    #    newlocs, newticks = formatRATicks(locs)
    #    plt.xticks(newlocs,newticks,rotation=45)
    #if(par2_name.lower()=='dec' or par2_name.lower()=='declination'):
    #    locs, ticks = plt.xticks()
    #    newlocs, newticks = formatDecTicks(locs)
    #    plt.xticks(newlocs,newticks,rotation=45)

    return fig
#

def plot_two_param_greedy_bins_hist(posterior,greedy2Params,confidence_levels):
    """
    Histograms of the ranked pixels produced by the 2-parameter greedy
    binning algorithm colured by their confidence level.

    @param posterior: an instance of the Posterior class.

    @param greedy2Params: a dict ;{param1Name:param1binSize,param2Name:param2binSize}

    @param confidence_levels: list of confidence levels to show
    """

    from scipy import seterr as sp_seterr

    np.seterr(under='ignore')
    sp_seterr(under='ignore')

    #Extract parameter names
    par1_name,par2_name=greedy2Params.keys()
    #Extract bin sizes
    par1_bin=greedy2Params[par1_name]
    par2_bin=greedy2Params[par2_name]

    a=np.squeeze(posterior[par1_name].samples)
    b=np.squeeze(posterior[par2_name].samples)

    #Extract injection information
    par1_injvalue=posterior[par1_name.lower()].injval
    par2_injvalue=posterior[par2_name.lower()].injval

    #Create 2D bin array
    par1pos_min=a.min()
    par2pos_min=b.min()

    par1pos_max=a.max()
    par2pos_max=b.max()

    par1pos_Nbins= int(ceil((par1pos_max - par1pos_min)/par1_bin))+1
    par2pos_Nbins= int(ceil((par2pos_max - par2pos_min)/par2_bin))+1

    # Adjust for time parameter
    if par1_name.find('time')!=-1:
        offset=floor(min(a))
        a=a-offset
        if par1_injvalue:
            par1_injvalue=par1_injvalue-offset
        ax1_name=par1_name+' + %i'%(int(offset))
    else: ax1_name=par1_name

    if par2_name.find('time')!=-1:
        offset=floor(min(b))
        b=b-offset
        if par2_injvalue:
            par2_injvalue=par2_injvalue-offset
        ax2_name=par2_name+' + %i'%(int(offset))
    else: ax2_name=par2_name


    #Extract trigger information
    par1_trigvalues=posterior[par1_name.lower()].trigvals
    par2_trigvalues=posterior[par2_name.lower()].trigvals

    myfig=plt.figure()
    axes=plt.Axes(myfig,[0.3,0.3,0.95-0.3,0.90-0.3])
    myfig.add_axes(axes)

    #plt.clf()
    plt.xlabel(plot_label(ax2_name))
    plt.ylabel(plot_label(ax1_name))

    #bins=(par1pos_Nbins,par2pos_Nbins)
    bins=(50,50) # Matches plot_one_param_pdf

    majorFormatterX=ScalarFormatter(useMathText=True)
    majorFormatterX.format_data=lambda data:'%.4g'%(data)
    majorFormatterY=ScalarFormatter(useMathText=True)
    majorFormatterY.format_data=lambda data:'%.4g'%(data)
    majorFormatterX.set_scientific(True)
    majorFormatterY.set_scientific(True)
    axes.xaxis.set_major_formatter(majorFormatterX)
    axes.yaxis.set_major_formatter(majorFormatterY)
    H, xedges, yedges = np.histogram2d(a,b, bins,density=False)


    #Replace H with greedy bin confidence levels at each pixel...
    temp=np.copy(H)
    temp=temp.flatten()

    Hsum=0
    Hsum_actual=np.sum(H)

    idxes=np.argsort(temp)
    j=len(idxes)-1
    while Hsum<Hsum_actual:
            #ind = np.argsort(temp)
        max_i=idxes[j]
        j-=1
        val = temp[max_i]
        Hsum+=int(val)
        temp[max_i]=0

        #print Hsum,Hsum_actual
        H.flat[max_i]=1-float(Hsum)/float(Hsum_actual)

    extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
    plt.imshow(np.flipud(H), axes=axes, aspect='auto', extent=extent, interpolation='nearest',cmap='gray_r')
    plt.gca().autoscale_view()
    plt.colorbar()

    #plt.hexbin(a,b,cmap='gray_r',axes=axes )

    Nchars=max(map(lambda d:len(majorFormatterX.format_data(d)),axes.get_xticks()))
    if Nchars>8:
        Nticks=3
    elif Nchars>5:
        Nticks=4
    elif Nchars>4:
        Nticks=5
    else:
        Nticks=6
    locatorX=matplotlib.ticker.MaxNLocator(nbins=Nticks-1)
    (xmin,xmax)=plt.xlim()
    (ymin,ymax)=plt.ylim()
    if par2_name=='rightascension' or par2_name=='ra':
        locatorX=RALocator(min=xmin,max=xmax)
        majorFormatterX=RAFormatter()
    if par2_name=='declination' or par2_name=='dec':
        locatorX=DecLocator(min=xmin,max=xmax)
        majorFormatterX=DecFormatter()
    if par1_name=='rightascension' or par1_name=='ra':
        locatorY=RALocator(min=ymin,max=ymax)
        axes.yaxis.set_major_locator(locatorY)
        majorFormatterY=RAFormatter()
    if par1_name=='declination' or par1_name=='dec':
        locatorY=DecLocator(min=ymin,max=ymax)
        axes.yaxis.set_major_locator(locatorY)
        majorFormatterY=DecFormatter()

    axes.xaxis.set_major_formatter(majorFormatterX)
    axes.yaxis.set_major_formatter(majorFormatterY)
    #locatorX.view_limits(bins[0],bins[-1])
    axes.xaxis.set_major_locator(locatorX)

    if par1_injvalue is not None and par2_injvalue is not None:
        plt.plot([par2_injvalue],[par1_injvalue],'bo',scalex=False,scaley=False)

    if par1_trigvalues is not None and par2_trigvalues is not None:
        par1IFOs = set([IFO for IFO in par1_trigvalues.keys()])
        par2IFOs = set([IFO for IFO in par2_trigvalues.keys()])
        IFOs = par1IFOs.intersection(par2IFOs)
        if IFO=='H1': color = 'r'
        elif IFO=='L1': color = 'g'
        elif IFO=='V1': color = 'm'
        else: color = 'c'
        plt.plot([par2_trigvalues[IFO]],[par1_trigvalues[IFO]],color=color,marker='o',scalex=False,scaley=False)

    # For RA and dec set custom labels and for RA reverse
    #if(par1_name.lower()=='ra' or par1_name.lower()=='rightascension'):
    #        ymin,ymax=plt.ylim()
    #        plt.ylim(ymax,ymin)
    #if(par1_name.lower()=='ra' or par1_name.lower()=='rightascension'):
    #        locs, ticks = plt.yticks()
    #        newlocs, newticks = formatRATicks(locs)
    #        plt.yticks(newlocs,newticks)
    #if(par1_name.lower()=='dec' or par1_name.lower()=='declination'):
    #        locs, ticks = plt.yticks()
    #        newlocs, newticks = formatDecTicks(locs)
    #        plt.yticks(newlocs,newticks)

    if(par2_name.lower()=='ra' or par2_name.lower()=='rightascension'):
        xmin,xmax=plt.xlim()
        if(xmin)<0.0: xmin=0.0
        if(xmax>2.0*pi_constant): xmax=2.0*pi_constant
        plt.xlim(xmax,xmin)
    #if(par2_name.lower()=='ra' or par2_name.lower()=='rightascension'):
    #    locs, ticks = plt.xticks()
    #    newlocs, newticks = formatRATicks(locs)
    #    plt.xticks(newlocs,newticks,rotation=45)
    #if(par2_name.lower()=='dec' or par2_name.lower()=='declination'):
    #    locs, ticks = plt.xticks()
    #    newlocs, newticks = formatDecTicks(locs)
    #    plt.xticks(newlocs,newticks,rotation=45)

    return myfig

def greedy_bin_one_param(posterior,greedy1Param,confidence_levels):
    """
    Determine the 1-parameter Bayesian Confidence Interval using a greedy
    binning algorithm.

    @param posterior: an instance of the posterior class.

    @param greedy1Param: a dict; {paramName:paramBinSize}.

    @param confidence_levels: A list of floats of the required confidence intervals [(0-1)].
    """

    paramName=list(greedy1Param.keys())[0]
    par_bin=list(greedy1Param.values())[0]
    par_samps=posterior[paramName.lower()].samples

    parpos_min=min(par_samps)[0]
    parpos_max=max(par_samps)[0]

    par_point=parpos_min

    parpos_Nbins= int(ceil((parpos_max - parpos_min)/par_bin))+1

    greedyPoints=np.zeros((parpos_Nbins,2))
    # ...NB 2D so it can be put through same confidence level function
    greedyHist=np.zeros(parpos_Nbins,dtype='i8')

    #Bin up
    for i in range(parpos_Nbins):
        greedyPoints[i,0]=par_point
        greedyPoints[i,1]=par_point
        par_point+=par_bin

    for par_samp in par_samps:
        par_samp=par_samp[0]
        par_binNumber=int(floor((par_samp-parpos_min)/par_bin))
        try:
            greedyHist[par_binNumber]+=1
        except IndexError:
            print("IndexError: bin number: %i total bins: %i parsamp: %f "\
                %(par_binNumber,parpos_Nbins,par_samp))

    #Find injection bin
    injbin=None
    par_injvalue=posterior[paramName].injval
    if par_injvalue:
        par_binNumber=floor((par_injvalue-parpos_min)/par_bin)
        injbin=par_binNumber

    toppoints,injectionconfidence,reses,injection_area=_greedy_bin(greedyHist,greedyPoints,injbin,float(par_bin),int(len(par_samps)),confidence_levels)
    cl_intervals=[]
    confidence_levels.sort()
    for cl in confidence_levels:
        ind=np.nonzero(toppoints[:,-1]<cl)

        if len(ind[0]) > 1:
            cl_intervals.append((np.min(toppoints[ind,0]),np.max(toppoints[ind,0])))

        else:

            cl_intervals.append((toppoints[ind[0],0],toppoints[ind[0],0]))

    return toppoints,injectionconfidence,reses,injection_area,cl_intervals

#
def contigious_interval_one_param(posterior,contInt1Params,confidence_levels):
    """
    Calculates the smallest contigious 1-parameter confidence interval for a
    set of given confidence levels.

    @param posterior: an instance of the Posterior class.

    @param contInt1Params: a dict {paramName:paramBinSize}.

    @param confidence_levels: Required confidence intervals.

    """
    oneDContCL={}
    oneDContInj={}

    paramName=list(contInt1Params.keys())[0]
    par_bin=list(contInt1Params.values())[0]

    par_injvalue=posterior[paramName].injval

    par_samps=posterior[paramName].samples

    parpos_min=min(par_samps)
    parpos_max=max(par_samps)

    par_point=parpos_min
    parpos_Nbins= int(ceil((parpos_max - parpos_min)/par_bin))+1

    greedyHist=np.zeros(parpos_Nbins,dtype='i8')

    for par_samp in par_samps:
        par_binNumber=int(floor((par_samp-parpos_min)/par_bin))
        try:
            greedyHist[par_binNumber]+=1
        except IndexError:
            print("IndexError: bin number: %i total bins: %i parsamp: %f bin: %f - %f"\
                %(
                  par_binNumber,
                  parpos_Nbins,
                  par_samp,
                  greedyPoints[par_binNumber-1,0],
                  greedyPoints[par_binNumber-1,0]+par_bin
                  ))

    injbin=None
    #Find injection bin
    if par_injvalue:
        par_binNumber=floor((par_injvalue-parpos_min)/par_bin)
        injbin=par_binNumber

    j=0
    #print "Calculating contigious confidence intervals for %s..."%par_name
    len_par_samps=len(par_samps)

    injinterval=None

    #Determine smallest contigious interval for given confidence levels (brute force)
    while j < len(confidence_levels):
        confidence_level=confidence_levels[j]
        #Loop over size of interval
        max_left=0
        max_right=0

        for i in range(len(greedyHist)):

            max_frac=None
            left=0
            right=i

            #Slide interval
            while right<len(greedyHist):
                Npoints=sum(greedyHist[left:right])
                frac=float(Npoints)/float(len_par_samps)
                #print "left %i , right %i , frac %f"%(left,right,frac)

                if max_frac is None:
                    max_frac=frac
                    max_left=left
                    max_right=right
                else:
                    if frac>max_frac:
                        max_frac=frac
                        max_left=left
                        max_right=right

                left+=1
                right+=1

            if injbin is not None and injinterval is None:
                if injbin in range(max_left,max_right):
                    injinterval=(max_right-max_left)*par_bin
                    oneDContInj['interval']=injinterval
                    oneDContInj['confidence']=1-frac
            if max_frac > confidence_level:
                break

            max_frac=None

        if max_frac is None:
            print("Cant determine intervals at %f confidence!"%confidence_level)
        else:

            oneDContCL['left']=max_left*par_bin
            oneDContCL['right']=max_right*par_bin
            oneDContCL['width']=(max_right-max_left)*par_bin
            k=j
            while k+1<len(confidence_levels) :
                if confidence_levels[k+1]<max_frac:
                    j+=1
                k+=1
        j+=1

    return oneDContCL,oneDContInj
#
def burnin(data,spin_flag,deltaLogP,outputfile):

    pos,bayesfactor=_burnin(data,spin_flag,deltaLogP,outputfile)

    return pos,bayesfactor


class ACLError(Exception):
    def __init__(self, *args):
        super(ACLError, self).__init__(*args)


def autocorrelation(series):
    """Returns an estimate of the autocorrelation function of a given
    series.  Returns only the positive-lag portion of the ACF,
    normalized so that the zero-th element is 1."""
    x=series-np.mean(series)
    y=np.conj(x[::-1])

    acf=np.fft.ifftshift(signal.fftconvolve(y,x,mode='full'))

    N=series.shape[0]

    acf = acf[0:N]

    return acf/acf[0]


def autocorrelation_length_estimate(series, acf=None, M=5, K=2):
    """Attempts to find a self-consistent estimate of the
    autocorrelation length of a given series.

    If C(tau) is the autocorrelation function (normalized so C(0) = 1,
    for example from the autocorrelation procedure in this module),
    then the autocorrelation length is the smallest s such that

    1 + 2*C(1) + 2*C(2) + ... + 2*C(M*s) < s

    In words: the autocorrelation length is the shortest length so
    that the sum of the autocorrelation function is smaller than that
    length over a window of M times that length.

    The maximum window length is restricted to be len(series)/K as a
    safety precaution against relying on data near the extreme of the
    lags in the ACF, where there is a lot of noise.  Note that this
    implies that the series must be at least M*K*s samples long in
    order to get a reliable estimate of the ACL.

    If no such s can be found, raises ACLError; in this case it is
    likely that the series is too short relative to its true
    autocorrelation length to obtain a consistent ACL estimate."""

    if acf is None:
        acf=autocorrelation(series)
    acf[1:] *= 2.0

    imax=int(acf.shape[0]/K)

    # Cumulative sum and ACL length associated with each window
    cacf=np.cumsum(acf)
    s=np.arange(1, cacf.shape[0]+1)/float(M)

    # Find all places where cumulative sum over window is smaller than
    # associated ACL.
    estimates=np.flatnonzero(cacf[:imax] < s[:imax])

    if estimates.shape[0] > 0:
        # Return the first index where cumulative sum is smaller than
        # ACL associated with that index's window
        return s[estimates[0]]
    else:
        # Cannot find self-consistent ACL estimate.
        raise ACLError('autocorrelation length too short for consistent estimate')


def effectiveSampleSize(samples, Nskip=1):
    """
    Compute the effective sample size, calculating the ACL using only
    the second half of the samples to avoid ACL overestimation due to
    chains equilibrating after adaptation.
    """
    N = len(samples)
    acf = autocorrelation(samples[int(N/2):])
    try:
        acl = autocorrelation_length_estimate(samples[int(N/2):], acf=acf)
    except ACLError:
        acl = N
    Neffective = floor(N/acl)
    acl *= Nskip
    return (Neffective, acl, acf)


def readCoincXML(xml_file, trignum):
    triggers=None

    from ligo.lw import ligolw
    from ligo.lw import lsctables
    from ligo.lw import utils
    coincXML = utils.load_filename(xml_file, contenthandler = lsctables.use_in(ligolw.LIGOLWContentHandler))
    coinc = lsctables.CoincTable.get_table(coincXML)
    coincMap = lsctables.CoincMapTable.get_table(coincXML)
    snglInsps = lsctables.SnglInspiralTable.get_table(coincXML)

    if (trignum>len(coinc)):
        raise RuntimeError("Error: You asked for trigger %d, but %s contains only %d triggers" %(trignum,coincfile,len(tiggers)))
    else:
        coincEventID = coinc.getColumnByName('coinc_event_id')[trignum]
        eventIDs = [row.event_id for row in coincMap if row.coinc_event_id == coincEventID]
        triggers = [row for row in snglInsps if row.event_id in eventIDs]
    return triggers

#===============================================================================
# Parameter estimation codes results parser
#===============================================================================

def find_ndownsample(samples, nDownsample):
    """
    Given a list of files, threshold value, and a desired
    number of outputs posterior samples, return the skip number to
    achieve the desired number of posterior samples.
    """
    if nDownsample is None:
        print("Max ACL(s):")
        splineParams=["spcal_npts", "spcal_active","constantcal_active"]
        for i in np.arange(25):
            for k in lal.cached_detector_by_prefix:
                splineParams.append(k.lower()+'_spcal_freq_'+str(i))
                splineParams.append(k.lower()+'_spcal_logfreq_'+str(i))

        nonParams = ["logpost", "post", "cycle", "timestamp", "snrh1", "snrl1", "snrv1",
                     "margtime","margtimephi","margtime","time_max","time_min",
                     "time_mean", "time_maxl","sky_frame","psdscaleflag","logdeltaf","flow","f_ref",
                     "lal_amporder","lal_pnorder","lal_approximant","tideo","spino","signalmodelflag",
                     "temperature","nifo","nlocaltemps","ntemps","randomseed","samplerate","segmentlength","segmentstart",
                     "t0", "phase_maxl", "azimuth", "cosalpha", "lal_amporder", "bluni"] + logParams + snrParams + splineParams
        fixedParams = [p for p in samples.colnames if all(x==samples[p][0] for x in samples[p])]
        print("Fixed parameters: "+str(fixedParams))
        nonParams.extend(fixedParams)
        params = [p for p in samples.colnames if p.lower() not in nonParams]
        stride=np.diff(samples['cycle'])[0]
        results = np.array([np.array(effectiveSampleSize(samples[param])[:2]) for param in params])
        nEffs = results[:,0]
        nEffective = min(nEffs)
        ACLs  = results[:,1]
        maxACLind = np.argmax(ACLs)
        maxACL = ACLs[maxACLind]
        # Get index in header, which includes "non-params"
        print("%i (%s)." %(stride*maxACL,params[maxACLind]))

    nskip = 1
    if nDownsample is not None:
        if len(samples) > nDownsample:
            nskip *= floor(len(samples)/nDownsample)
            nskip = int(nskip)
    else:
        nEff = nEffective
        if nEff > 1:
            if len(samples) > nEff:
                nskip = int(ceil(len(samples)/nEff))
        else:
            nskip = np.nan
    return nskip

class PEOutputParser(object):
    """
    A parser for the output of Bayesian parameter estimation codes.

    TODO: Will be abstract class when LDG moves over to Python >2.6,
    inherited by each method .
    """
    def __init__(self,inputtype):
        if inputtype is 'mcmc_burnin':
            self._parser=self._mcmc_burnin_to_pos
        elif inputtype is 'ns':
            self._parser=self._ns_to_pos
        elif inputtype is 'common':
            self._parser=self._common_to_pos
        elif inputtype is 'fm':
            self._parser=self._followupmcmc_to_pos
        elif inputtype is "inf_mcmc":
            self._parser=self._infmcmc_to_pos
        elif inputtype is "xml":
            self._parser=self._xml_to_pos
        elif inputtype == 'hdf5':
            self._parser = self._hdf5_to_pos
        elif inputtype == 'hdf5s':
            self._parser = self._hdf5s_to_pos
        else:
            raise ValueError('Invalid value for "inputtype": %r' % inputtype)

    def parse(self,files,**kwargs):
        """
        Parse files.
        """
        return self._parser(files,**kwargs)

    def _infmcmc_to_pos(self,files,outdir=None,deltaLogP=None,fixedBurnins=None,nDownsample=None,oldMassConvention=False,**kwargs):
        """
        Parser for lalinference_mcmcmpi output.
        """
        if not (fixedBurnins is None):
            if not (deltaLogP is None):
                print("Warning: using deltaLogP criteria in addition to fixed burnin")
            if len(fixedBurnins) == 1 and len(files) > 1:
                print("Only one fixedBurnin criteria given for more than one output.  Applying this to all outputs.")
                fixedBurnins = np.ones(len(files),'int')*fixedBurnins[0]
            elif len(fixedBurnins) != len(files):
                raise RuntimeError("Inconsistent number of fixed burnin criteria and output files specified.")
            print("Fixed burning criteria: ",fixedBurnins)
        else:
            fixedBurnins = np.zeros(len(files))
        logPThreshold=-np.inf
        if not (deltaLogP is None):
            logPThreshold= - deltaLogP
            print("Eliminating any samples before log(Post) = ", logPThreshold)
        nskips=self._find_ndownsample(files, logPThreshold, fixedBurnins, nDownsample)
        if nDownsample is None:
            print("Downsampling to take only uncorrelated posterior samples from each file.")
            if len(nskips) == 1 and np.isnan(nskips[0]):
                print("WARNING: All samples in chain are correlated.  Downsampling to 10000 samples for inspection!!!")
                nskips=self._find_ndownsample(files, logPThreshold, fixedBurnins, 10000)
            else:
                for i in range(len(nskips)):
                    if np.isnan(nskips[i]):
                        print("%s eliminated since all samples are correlated.")
                    else:
                        print("Downsampling by a factor of ", nskips[0], " to achieve approximately ", nDownsample, " posterior samples")
        if outdir is None:
            outdir=''
        runfileName=os.path.join(outdir,"lalinfmcmc_headers.dat")
        postName="posterior_samples.dat"
        runfile=open(runfileName, 'w')
        outfile=open(postName, 'w')
        try:
            self._infmcmc_output_posterior_samples(files, runfile, outfile, logPThreshold, fixedBurnins, nskips, oldMassConvention)
        finally:
            runfile.close()
            outfile.close()
        return self._common_to_pos(open(postName,'r'))


    def _infmcmc_output_posterior_samples(self, files, runfile, outfile, logPThreshold, fixedBurnins, nskips=None, oldMassConvention=False):
        """
        Concatenate all the samples from the given files into outfile.
        For each file, only those samples past the point where the
        log(post) > logPThreshold are concatenated after eliminating
        fixedBurnin.
        """
        nRead=0
        outputHeader=False
        acceptedChains=0
        if nskips is None:
            nskips = np.ones(len(files),'int')
        for infilename,i,nskip,fixedBurnin in zip(files,range(1,len(files)+1),nskips,fixedBurnins):
            infile=open(infilename,'r')
            try:
                print("Writing header of %s to %s"%(infilename,runfile.name))
                runInfo,header=self._clear_infmcmc_header(infile)
                runfile.write('Chain '+str(i)+':\n')
                runfile.writelines(runInfo)
                print("Processing file %s to %s"%(infilename,outfile.name))
                write_fref = False
                if 'f_ref' not in header:
                    write_fref = True
                    f_ref=self._find_infmcmc_f_ref(runInfo)
                if oldMassConvention:
                    # Swap #1 for #2 because our old mass convention
                    # has m2 > m1, while the common convention has m1
                    # > m2
                    header=[self._swaplabel12(label) for label in header]
                if not outputHeader:
                    for label in header:
                        outfile.write(label)
                        outfile.write(" ")
                    if write_fref:
                        outfile.write("f_ref")
                    outfile.write(" ")
                    outfile.write("chain")
                    outfile.write("\n")
                    outputHeader=header
                iterindex=header.index("cycle")
                logpindex=header.index("logpost")
                output=False
                for line in infile:
                    line=line.lstrip()
                    lineParams=line.split()
                    iter=int(lineParams[iterindex])
                    logP=float(lineParams[logpindex])
                    if (iter > fixedBurnin) and (logP >= logPThreshold):
                        output=True
                    if output:
                        if nRead % nskip == 0:
                            for label in outputHeader:
                                # Note that the element "a1" in the
                                # *header* actually already
                                # corresponds to the "a2" *column* of
                                # the input because we switched the
                                # names above
                                outfile.write(lineParams[header.index(label)])
                                outfile.write("\t")
                            if write_fref:
                                outfile.write(f_ref)
                                outfile.write("\t")
                            outfile.write(str(i))
                            outfile.write("\n")
                        nRead=nRead+1
                if output: acceptedChains += 1
            finally:
                infile.close()
        print("%i of %i chains accepted."%(acceptedChains,len(files)))

    def _swaplabel12(self, label):
        if label[-1] == '1':
            return label[0:-1] + '2'
        elif label[-1] == '2':
            return label[0:-1] + '1'
        else:
            return label[:]

    def _find_max_logP(self, files):
        """
        Given a list of files, reads them, finding the maximum log(post)
        """
        maxLogP = -np.inf
        for inpname in files:
            infile=open(inpname, 'r')
            try:
                runInfo,header=self._clear_infmcmc_header(infile)
                logpindex=header.index("logpost")
                for line in infile:
                    line=line.lstrip().split()
                    logP=float(line[logpindex])
                    if logP > maxLogP:
                        maxLogP=logP
            finally:
                infile.close()
        print("Found max log(post) = ", maxLogP)
        return maxLogP

    def _find_ndownsample(self, files, logPthreshold, fixedBurnins, nDownsample):
        """
        Given a list of files, threshold value, and a desired
        number of outputs posterior samples, return the skip number to
        achieve the desired number of posterior samples.
        """
        nfiles = len(files)
        ntots=[]
        nEffectives = []
        if nDownsample is None: print("Max ACL(s):")
        for inpname,fixedBurnin in zip(files,fixedBurnins):
            infile = open(inpname, 'r')
            try:
                runInfo,header = self._clear_infmcmc_header(infile)
                header = [name.lower() for name in header]
                logpindex = header.index("logpost")
                iterindex = header.index("cycle")
                deltaLburnedIn = False
                fixedBurnedIn  = False
                adapting = True
                lines=[]
                ntot=0
                for line in infile:
                    line = line.lstrip().split()
                    iter = int(line[iterindex])
                    logP = float(line[logpindex])
                    if iter > fixedBurnin:
                        fixedBurnedIn = True
                    # If adaptation reset, throw out what was collected so far
                    elif fixedBurnedIn:
                        fixedBurnedIn = False
                        ntot = 0
                        lines = []
                    if logP > logPthreshold:
                        deltaLburnedIn = True
                    if iter > 0:
                        adapting = False
                    if fixedBurnedIn and deltaLburnedIn and not adapting:
                        ntot += 1
                        lines.append(line)
                ntots.append(ntot)
                if nDownsample is None:
                    try:
                        splineParams=["spcal_npts", "spcal_active","constantcal_active"]
                        for i in np.arange(5):
                            for k in ['h1','l1']:
                                splineParams.append(k+'_spcal_freq'+str(i))
                                splineParams.append(k+'_spcal_logfreq'+str(i))

                        nonParams = ["logpost", "cycle", "timestamp", "snrh1", "snrl1", "snrv1",
                                     "margtime","margtimephi","margtime","time_max","time_min",
                                     "time_mean", "time_maxl","sky_frame","psdscaleflag","logdeltaf","flow","f_ref",
                                     "lal_amporder","lal_pnorder","lal_approximant","tideo","spino","signalmodelflag",
                                     "temperature","nifo","nlocaltemps","ntemps","randomseed","samplerate","segmentlength","segmentstart",
                                     "t0", "phase_maxl", "azimuth", "cosalpha"] + logParams + snrParams + splineParams
                        nonParamsIdxs = [header.index(name) for name in nonParams if name in header]
                        samps = np.array(lines).astype(float)
                        fixedIdxs = np.where(np.amin(samps,axis=0)-np.amax(samps,axis=0) == 0.0)[0]
                        nonParamsIdxs.extend(fixedIdxs)
                        paramIdxs = [i for i in range(len(header)) if i not in nonParamsIdxs]
                        stride=samps[1,iterindex] - samps[0,iterindex]
                        results = np.array([np.array(effectiveSampleSize(samps[:,i])[:2]) for i in paramIdxs])
                        nEffs = results[:,0]
                        nEffectives.append(min(nEffs))
                        ACLs  = results[:,1]
                        maxACLind = np.argmax(ACLs)
                        maxACL = ACLs[maxACLind]
                        # Get index in header, which includes "non-params"
                        maxACLind = paramIdxs[maxACLind]
                        print("%i (%s) for chain %s." %(stride*maxACL,header[maxACLind],inpname))
                    except:
                        nEffectives.append(None)
                        print("Error computing effective sample size of %s!"%inpname)

            finally:
                infile.close()
        nskips = np.ones(nfiles)
        ntot = sum(ntots)
        if nDownsample is not None:
            if ntot > nDownsample:
                nskips *= int(floor(ntot/nDownsample))
        else:
            for i in range(nfiles):
                nEff = nEffectives[i]
                ntot = ntots[i]
                if nEff > 1:
                    if ntot > nEff:
                        nskips[i] = int(ceil(ntot/nEff))
                else:
                    nskips[i] = None
        return nskips

    def _find_infmcmc_f_ref(self, runInfo):
        """
        Searches through header to determine reference frequency of waveforms.
        If no fRef given, calls _find_infmcmc_f_lower to get the lower frequency
        bound, which is the default reference frequency for LALInference.
        """
        fRef = None
        runInfoIter = iter(runInfo)
        for line in runInfoIter:
            headers=line.lstrip().lower().split()
            try:
                fRefColNum = headers.index('f_ref') # strings get converted to all lower case
                info = runInfoIter.next().lstrip().lower().split()
                fRef = info[-1]#fRefColNum] # too many column names with spaces for this way to work. I just grab the last value. Hopefully we will update to xml output files and those messy headers will be gone.
                break
            except ValueError:
                continue

        # ***TEMPORARY*** If not in table, check command line.
        #   ...This is messy, but the only option for dealing with old headers
        if not fRef:
            runInfoIter = iter(runInfo)
            for line in runInfoIter:
                headers=line.lstrip().lower().split()
                try:
                    if headers[0]=="command":
                        try:
                            fRefInd = headers.index('--fref')+1
                            fRef = headers[fRefInd]
                        except ValueError:
                            pass
                        break
                except IndexError:
                    continue

        # If no fRef is found, use lower frequency bound
        if not fRef:
            fRef = self._find_infmcmc_f_lower(runInfo)

        return fRef

    def _find_infmcmc_f_lower(self, runInfo):
        """
        Searches through header to determine starting frequency of waveforms.
        Assumes same for all IFOs.
        """
        runInfo = iter(runInfo)
        for line in runInfo:
            headers=line.lstrip().lower().split()
            try:
                flowColNum = headers.index('f_low')
                IFOinfo = runInfo.next().lstrip().lower().split()
                f_lower = IFOinfo[flowColNum]
                break
            except ValueError:
                continue
        return f_lower

    def _clear_infmcmc_header(self, infile):
        """
        Reads lalinference_mcmcmpi file given, returning the run info and
        common output header information.
        """
        runInfo = []
        for line in infile:
            runInfo.append(line)
            headers=line.lstrip().lower().split()
            try:
                headers.index('cycle')
                break
            except ValueError:
                continue
        else:
            raise RuntimeError("couldn't find line with 'cycle' in LALInferenceMCMC input")
        return runInfo[:-1],headers


    def _mcmc_burnin_to_pos(self,files,spin=False,deltaLogP=None):
        """
        Parser for SPINspiral output .
        """
        raise NotImplementedError
        if deltaLogP is not None:
            pos,bayesfactor=burnin(data,spin,deltaLogP,"posterior_samples.dat")
            return self._common_to_pos(open("posterior_samples.dat",'r'))

    def _ns_to_pos(self,files,Nlive=None,Npost=None,posfilename='posterior_samples.dat'):
        """
        Parser for nested sampling output.
        files : list of input NS files
        Nlive : Number of live points
        Npost : Desired number of posterior samples
        posfilename : Posterior output file name (default: 'posterior_samples.dat')
        """
        try:
            from lalinference.nest2pos import draw_N_posterior_many,draw_posterior_many
        except ImportError:
            print("Need lalinference.nest2pos to convert nested sampling output!")
            raise

        if Nlive is None:
            raise RuntimeError("Need to specify number of live points in positional arguments of parse!")

        #posfile.write('mchirp \t eta \t time \t phi0 \t dist \t RA \t dec \t
        #psi \t iota \t likelihood \n')
        # get parameter list
        it = iter(files)

        # check if there's a file containing the parameter names
        parsfilename = (it.next()).strip('.gz')+'_params.txt'

        if os.path.isfile(parsfilename):
            print('Looking for '+parsfilename)

            if os.access(parsfilename,os.R_OK):

                with open(parsfilename,'r') as parsfile:
                    outpars=parsfile.readline()+'\n'
            else:
                raise RuntimeError('Cannot open parameters file %s!'%(parsfilename))

        else: # Use hardcoded CBC parameter names
            outpars='mchirp \t eta \t time \t phi0 \t dist \t RA \t \
            dec \t psi \t iota \t logl \n'

        # Find the logL column
        parsvec=outpars.split()
        logLcol=-1
        for i in range(len(parsvec)):
            if parsvec[i].lower()=='logl':
                logLcol=i
        if logLcol==-1:
            print('Error! Could not find logL column in parameter list: %s'%(outpars))
            raise RuntimeError

        inarrays=list(map(np.loadtxt,files))
        if Npost is None:
            pos=draw_posterior_many(inarrays,[Nlive for f in files],logLcols=[logLcol for f in files])
        else:
            pos=draw_N_posterior_many(inarrays,[Nlive for f in files],Npost,logLcols=[logLcol for f in files])

        with open(posfilename,'w') as posfile:

            posfile.write(outpars)

            for row in pos:
                for i in row:
                    posfile.write('%10.12e\t' %(i))
                posfile.write('\n')

        with open(posfilename,'r') as posfile:
            return_val=self._common_to_pos(posfile)

        return return_val

    def _followupmcmc_to_pos(self,files):
        """
        Parser for followupMCMC output.
        """
        return self._common_to_pos(open(files[0],'r'),delimiter=',')


    def _multinest_to_pos(self,files):
        """
        Parser for MultiNest output.
        """
        return self._common_to_pos(open(files[0],'r'))

    def _xml_to_pos(self,infile):
        """
        Parser for VOTable XML Using
        """
        from xml.etree import ElementTree as ET
        xmlns='http://www.ivoa.net/xml/VOTable/v1.1'
        try:
            register_namespace=ET.register_namespace
        except AttributeError:
            def register_namespace(prefix,uri):
                ET._namespace_map[uri]=prefix
        register_namespace('vot',xmlns)
        tree = ET.ElementTree()

        tree.parse(infile)
        # Find the posterior table
        tables = tree.findall('.//{%s}TABLE'%(xmlns))
        for table in tables:
            if table.get('utype')=='lalinference:results:posteriorsamples':
                return(self._VOTTABLE2pos(table))
        for table in tables:
            if table.get('utype')=='lalinference:results:nestedsamples':
                nsresource=[node for node in tree.findall('{%s}RESOURCE'%(xmlns)) if node.get('utype')=='lalinference:results'][0]
                return(self._VOTTABLE2pos(vo_nest2pos(nsresource)))
        raise RuntimeError('Cannot find "Posterior Samples" TABLE element in XML input file %s'%(infile))

    def _VOTTABLE2pos(self,table):
        """
        Parser for a VOT TABLE element with FIELDs and TABLEDATA elements
        """
        from xml.etree import ElementTree as ET
        xmlns='http://www.ivoa.net/xml/VOTable/v1.1'
        try:
            register_namespace=ET.register_namespace
        except AttributeError:
            def register_namespace(prefix,uri):
                ET._namespace_map[uri]=prefix
                register_namespace('vot',xmlns)

        header=[]
        for field in table.findall('./{%s}FIELD'%(xmlns)):
            header.append(field.attrib['name'])
        if(len(header)==0):
            raise RuntimeError('Unable to find FIELD nodes for table headers in XML table')
        data=table.findall('./{%s}DATA'%(xmlns))
        tabledata=data[0].find('./{%s}TABLEDATA'%(xmlns))
        llines=[]
        for row in tabledata:
            llines.append(np.array(list(map(lambda a:float(a.text),row))))
        flines=np.array(llines)
        for i in range(0,len(header)):
            if header[i].lower().find('log')!=-1 and header[i].lower() not in logParams and re.sub('log', '', header[i].lower()) not in [h.lower() for h in header] and header[i].lower() not in lorentzInvarianceViolationParams:
                print('exponentiating %s'%(header[i]))

                flines[:,i]=np.exp(flines[:,i])

                header[i]=re.sub('log', '', header[i], flags=re.IGNORECASE)
            if header[i].lower().find('sin')!=-1 and re.sub('sin', '', header[i].lower()) not in [h.lower() for h in header]:
                print('asining %s'%(header[i]))
                flines[:,i]=np.arcsin(flines[:,i])
                header[i]=re.sub('sin', '', header[i], flags=re.IGNORECASE)
            if header[i].lower().find('cos')!=-1 and re.sub('cos', '', header[i].lower()) not in [h.lower() for h in header]:
                print('acosing %s'%(header[i]))
                flines[:,i]=np.arccos(flines[:,i])
                header[i]=re.sub('cos', '', header[i], flags=re.IGNORECASE)
            header[i]=header[i].replace('(','')
            header[i]=header[i].replace(')','')
        print('Read columns %s'%(str(header)))
        return header,flines

    def _hdf5s_to_pos(self, infiles, fixedBurnins=None, deltaLogP=None, nDownsample=None, tablename=None, **kwargs):
        from astropy.table import vstack

        if fixedBurnins is None:
            fixedBurnins = np.zeros(len(infiles))

        if len(infiles) > 1:
            multiple_chains = True

        chains = []
        for i, [infile, fixedBurnin] in enumerate(zip(infiles, fixedBurnins)):
            chain = self._hdf5_to_table(infile, fixedBurnin=fixedBurnin, deltaLogP=deltaLogP, nDownsample=nDownsample, multiple_chains=multiple_chains, tablename=tablename, **kwargs)
            chain.add_column(astropy.table.Column(i*np.ones(len(chain)), name='chain'))
            chains.append(chain)

        # Apply deltaLogP criteria across chains
        if deltaLogP is not None:
            logPThreshold = -np.inf
            for chain in chains:
                if len(chain) > 0:
                    logPThreshold = max([logPThreshold, max(chain['logpost'])- deltaLogP])
            print("Eliminating any samples before log(L) = {}".format(logPThreshold))

        for i, chain in enumerate(chains):
            if deltaLogP is not None:
                above_threshold = np.arange(len(chain))[chain['logpost'] > logPThreshold]
                burnin_idx = above_threshold[0] if len(above_threshold) > 0 else len(chain)
            else:
                burnin_idx = 0
            chains[i] = chain[burnin_idx:]

        samples = vstack(chains)

        # Downsample one more time
        if nDownsample is not None:
            nskip = find_ndownsample(samples, nDownsample)
            samples = samples[::nskip]

        return samples.colnames, as_array(samples).view(float).reshape(-1, len(samples.columns))

    def _hdf5_to_table(self, infile, deltaLogP=None, fixedBurnin=None, nDownsample=None, multiple_chains=False, tablename=None, **kwargs):
        """
        Parse a HDF5 file and return an array of posterior samples ad list of
        parameter names. Equivalent to '_common_to_pos' and work in progress.
        """
        if not tablename:
            samples = read_samples(infile, tablename=posterior_grp_name)
        else:
            samples = read_samples(infile, tablename=tablename)
        params = samples.colnames

        for param in params:
            param_low = param.lower()
            if param_low.find('log') != -1 and param_low not in logParams and re.sub('log', '', param_low) not in [p.lower() for p in params] and param_low not in lorentzInvarianceViolationParams:
                print('exponentiating %s' % param)
                new_param = re.sub('log', '', param, flags=re.IGNORECASE)
                samples[new_param] = np.exp(samples[param])
                del samples[param]
                param = new_param
            if param_low.find('sin') != -1 and re.sub('sin', '', param_low) not in [p.lower() for p in params]:
                print('asining %s' % param)
                new_param = re.sub('sin', '', param, flags=re.IGNORECASE)
                samples[new_param] = np.arcsin(samples[param])
                del samples[param]
                param = new_param
            if param_low.find('cos') != -1 and re.sub('cos', '', param_low) not in [p.lower() for p in params]:
                print('acosing %s' % param)
                new_param = re.sub('cos', '', param, flags=re.IGNORECASE)
                samples[new_param] = np.arccos(samples[param])
                del samples[param]
                param = new_param

            if param != param.replace('(', ''):
                samples.rename_column(param, param.replace('(', ''))
            if param != param.replace(')', ''):
                samples.rename_column(param, param.replace(')', ''))

            #Make everything a float, since that's what's excected of a CommonResultsObj
            replace_column(samples, param, samples[param].astype(float))

        params = samples.colnames
        print('Read columns %s' % str(params))

        # MCMC burnin and downsampling
        if 'cycle' in params:
            if not (fixedBurnin is None):
                if not (deltaLogP is None):
                    print("Warning: using deltaLogP criteria in addition to fixed burnin")
                print("Fixed burning criteria: ",fixedBurnin)
            else:
                fixedBurnin = 0

            burned_in_cycles = np.arange(len(samples))[samples['cycle'] > fixedBurnin]
            burnin_idx = burned_in_cycles[0] if len(burned_in_cycles) > 0 else len(samples)
            samples = samples[burnin_idx:]

            logPThreshold=-np.inf
            if len(samples) > 0 and not (deltaLogP is None):
                logPThreshold = max(samples['logpost'])- deltaLogP
                print("Eliminating any samples before log(post) = ", logPThreshold)
                burnin_idx = np.arange(len(samples))[samples['logpost'] > logPThreshold][0]
                samples = samples[burnin_idx:]

            if len(samples) > 0:
                nskip = find_ndownsample(samples, nDownsample)
                if nDownsample is None:
                    print("Downsampling to take only uncorrelated posterior samples from each file.")
                    if np.isnan(nskip) and not multiple_chains:
                        print("WARNING: All samples in chain are correlated.  Downsampling to 10000 samples for inspection!!!")
                        nskip = find_ndownsample(samples, 10000)
                        samples = samples[::nskip]
                else:
                    if np.isnan(nskip):
                        print("WARNING: All samples in {} are correlated.".format(infile))
                        samples = samples[-1:]
                    else:
                        print("Downsampling by a factor of ", nskip, " to achieve approximately ", nDownsample, " posterior samples")
                        samples = samples[::nskip]

        return samples

    def _hdf5_to_pos(self, infile, fixedBurnins=None, deltaLogP=None, nDownsample=None, tablename=None, **kwargs):
        samples = self._hdf5_to_table(infile, fixedBurnin=fixedBurnins, deltaLogP=deltaLogP, nDownsample=nDownsample, tablename=tablename, **kwargs)

        return samples.colnames, as_array(samples).view(float).reshape(-1, len(samples.columns))

    def _common_to_pos(self,infile,info=[None,None]):
        """
        Parse a file in the 'common format' and return an array of posterior
        samples and list of parameter names. Will apply inverse functions to
        columns with names containing sin,cos,log.
        """

        [headerfile,delimiter]=info

        if headerfile==None:
            formatstr=infile.readline().lstrip()
        else:
            hf=open(headerfile,'r')
            formatstr=hf.readline().lstrip()
            hf.close()

        formatstr=formatstr.replace('#','')
        formatstr=formatstr.replace('"','')

        header=formatstr.split(delimiter)
        header[-1]=header[-1].rstrip('\n')
        nparams=len(header)
        llines=[]
        dec=re.compile(r'^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$|^inf$')

        for line_number,line in enumerate(infile):
            sline=line.split(delimiter)
            if sline[-1] == '\n':
                del(sline[-1])
            proceed=True
            if len(sline)<1:
                print('Ignoring empty line in input file: %s'%(sline))
                proceed=False
            elif len(sline)!=nparams:
                sys.stderr.write('WARNING: Malformed row %i, read %i elements but there is meant to be %i\n'%(line_number,len(sline),nparams))
                proceed=False

            for elemn,st in enumerate(sline):
                s=st.replace('\n','')
                if dec.search(s) is None:
                    print('Warning! Ignoring non-numeric data after the header: %s. Row = %i,Element=%i'%(s,line_number,elemn))
                    proceed=False
                elif s is '\n':
                    proceed=False

            if proceed:
                llines.append(list(map(float,sline)))

        flines=np.array(llines)

        if not flines.any():
            raise RuntimeError("ERROR: no lines read in!")


        for i in range(0,len(header)):
            if header[i].lower().find('log')!=-1 and header[i].lower() not in logParams and re.sub('log', '', header[i].lower()) not in [h.lower() for h in header] and header[i].lower() not in lorentzInvarianceViolationParams:
                print('exponentiating %s'%(header[i]))

                flines[:,i]=np.exp(flines[:,i])

                header[i]=re.sub('log', '', header[i], flags=re.IGNORECASE)
            if header[i].lower().find('sin')!=-1 and re.sub('sin', '', header[i].lower()) not in [h.lower() for h in header]:
                print('asining %s'%(header[i]))
                flines[:,i]=np.arcsin(flines[:,i])
                header[i]=re.sub('sin', '', header[i], flags=re.IGNORECASE)
            if header[i].lower().find('cos')!=-1 and re.sub('cos', '', header[i].lower()) not in [h.lower() for h in header]:
                print('acosing %s'%(header[i]))
                flines[:,i]=np.arccos(flines[:,i])
                header[i]=re.sub('cos', '', header[i], flags=re.IGNORECASE)
            header[i]=header[i].replace('(','')
            header[i]=header[i].replace(')','')
        print('Read columns %s'%(str(header)))
        return header,flines


def parse_converge_output_section(fo):
    result={}
    lines=fo.split('\n')
    chain_line=False
    for line in lines:

        if '[1]' in line:
            key=line.replace('[1]','').strip(' ').strip('"')
            result[key]={}
            out=result[key]
            continue
        if result is not {}:
            if 'chain' in line:
                chain_line=True
                continue
            if chain_line:
                chain_line=False
                key=line.strip('"').split()[1]
                out[key]=[]
                out2=out[key]
            else:
                try:
                    newline=line.strip('"').split()
                    if newline is not []:
                        out2.append(line.strip('"').split())
                except:
                    pass

    return result
#

def vo_nest2pos(nsresource,Nlive=None):
    """
    Parse a VO Table RESOURCE containing nested sampling output and
    return a VOTable TABLE element with posterior samples in it.
    This can be added to an existing tree by the user.
    Nlive will be read from the nsresource, unless specified
    """
    from xml.etree import ElementTree as ET
    import copy
    from math import log, exp
    xmlns='http://www.ivoa.net/xml/VOTable/v1.1'
    try:
        register_namespace=ET.register_namespace
    except AttributeError:
        def register_namespace(prefix,uri):
            ET._namespace_map[uri]=prefix
    register_namespace('vot',xmlns)

    postable=ET.Element("{%s}TABLE"%(xmlns),attrib={'name':'Posterior Samples','utype':'lalinference:results:posteriorsamples'})
    i=0
    nstables=[resource for resource in nsresource.findall("./{%s}TABLE"%(xmlns)) if resource.get("utype")=="lalinference:results:nestedsamples"]

    nstable=nstables[0]
    if Nlive is None:
        runstateResource = [resource for resource in nsresource.findall("./{%s}RESOURCE"%(xmlns)) if resource.get("utype")=="lalinference:state"][0]
        algTable = [table for table in runstateResource.findall("./{%s}TABLE"%(xmlns)) if table.get("utype")=="lalinference:state:algorithmparams"][0]
        Nlive = int ([param for param in algTable.findall("./{%s}PARAM"%(xmlns)) if param.get("name")=='Nlive'][0].get('value'))
        print('Found Nlive %i'%(Nlive))
    if Nlive is None:
        raise RuntimeError("Cannot find number of live points in XML table, please specify")
    logLcol = None
    for fieldnode in nstable.findall('./{%s}FIELD'%xmlns):
        if fieldnode.get('name') == 'logL':
            logLcol=i
        i=i+1
        postable.append(copy.deepcopy(fieldnode))
    for paramnode in nstable.findall('./{%s}PARAM'%(xmlns)):
        postable.append(copy.deepcopy(paramnode))
    if logLcol is None:
        RuntimeError("Unable to find logL column")
    posdataNode=ET.Element("{%s}DATA"%(xmlns))
    postabledataNode=ET.Element("{%s}TABLEDATA"%(xmlns))
    postable.append(posdataNode)
    posdataNode.append(postabledataNode)
    nstabledata=nstable.find('./{%s}DATA/{%s}TABLEDATA'%(xmlns,xmlns))
    logw=log(1.0 - exp(-1.0/float(Nlive)))
    weights=[]
    for row in nstabledata:
        logL=float(row[logLcol].text)
        weights.append(logL-logw)
        logw=logw-1.0/float(Nlive)
    mw=max(weights)
    weights = [w - mw for w in weights]
    for (row,weight) in zip(nstabledata,weights):
        if weight > log(random.random()):
            postabledataNode.append(copy.deepcopy(row))
    return postable

xmlns='http://www.ivoa.net/xml/VOTable/v1.1'

class VOT2HTML:
    def __init__(self):
        self.html=htmlSection("VOTable information")
        self.skiptable=0

    def start(self,tag,attrib):
        if tag=='{%s}TABLE'%(xmlns):
            if attrib['utype']=='lalinference:results:nestedsamples'\
            or attrib['utype']=='lalinference:results:posteriorsamples':
                self.skiptable=1
            else:
                self.skiptable=0
            self.tableouter=htmlChunk('div')
            self.tableouter.h2(attrib['name'])
            try:
                self.tableouter.p(attrib['utype'])
            except KeyError:
                pass
            self.fixedparams=htmlChunk('table',attrib={'class':'statstable'},parent=self.tableouter)
            self.table=htmlChunk('table',attrib={'class':'statstable'},parent=self.tableouter)
            self.tabheader=htmlChunk('tr',parent=self.table)
        if tag=='{%s}FIELD'%(xmlns):
            self.field=htmlChunk('th',{'name':attrib['name']},parent=self.tabheader)
        if tag=='{%s}TR'%(xmlns):
            self.tabrow=htmlChunk('tr',parent=self.table)
        if tag=='{%s}TD'%(xmlns):
            self.td=htmlChunk('td',parent=self.tabrow)
        if tag=='{%s}PARAM'%(xmlns):
            pnode=htmlChunk('tr',parent=self.fixedparams)
            namenode=htmlChunk('td',parent=pnode)
            namenode.p(attrib['name'])
            valnode=htmlChunk('td',parent=pnode)
            valnode.p(attrib['value'])

    def end(self,tag):
        if tag=='{%s}TABLE'%(xmlns):
            if not self.skiptable:
                self.html.append(self.tableouter._html)
        if tag=='{%s}FIELD'%(xmlns):
            self.field.p(self.data)
        if tag=='{%s}TD'%(xmlns):
            self.td.p(self.data)

    def data(self,data):
        self.data=data

    def close(self):
        return self.html.toprettyxml()

def _cl_width(cl_bound):
    """Returns (high - low), the width of the given confidence
    bounds."""

    return cl_bound[1] - cl_bound[0]

def _cl_count(cl_bound, samples):
    """Returns the number of samples within the given confidence
    bounds."""

    return np.sum((samples >= cl_bound[0]) & (samples <= cl_bound[1]))

def confidence_interval_uncertainty(cl, cl_bounds, posteriors):
    """Returns a tuple (relative_change, fractional_uncertainty,
    percentile_uncertainty) giving the uncertainty in confidence
    intervals from multiple posteriors.

    The uncertainty in the confidence intervals is the difference in
    length between the widest interval, formed from the smallest to
    largest values among all the cl_bounds, and the narrowest
    interval, formed from the largest-small and smallest-large values
    among all the cl_bounds.  Note that neither the smallest nor the
    largest confidence intervals necessarily correspond to one of the
    cl_bounds.

    The relative change relates the confidence interval uncertainty to
    the expected value of the parameter, the fractional uncertainty
    relates this length to the length of the confidence level from the
    combined posteriors, and the percentile uncertainty gives the
    change in percentile over the combined posterior between the
    smallest and largest confidence intervals.

    @param cl The confidence level (between 0 and 1).

    @param cl_bounds A list of (low, high) pairs giving the confidence
    interval associated with each posterior.

    @param posteriors A list of PosteriorOneDPDF objects giving the
    posteriors."""

    Ns=[p.samples.shape[0] for p in posteriors]
    Nsamplers=len(Ns)

    # Weight each sample within a run equally, and each run equally
    # with respect to the others
    all_samples = np.squeeze(np.concatenate([p.samples for p in posteriors], axis=0))
    weights = np.squeeze(np.concatenate([p.samples*0.0+1.0/(Nsamplers*N) for (N,p) in zip(Ns,posteriors)], axis=0))

    isort=np.argsort(all_samples)

    all_samples = all_samples[isort]
    weights = weights[isort]

    param_mean = np.average(all_samples, weights=weights)

    N=all_samples.shape[0]

    alpha = (1.0 - cl)/2.0

    wttotal = np.cumsum(weights)
    ilow = np.nonzero(wttotal >= alpha)[0][0]
    ihigh = np.nonzero(wttotal >= 1.0-alpha)[0][0]

    all_cl_bound = (all_samples[ilow], all_samples[ihigh])

    low_bounds = np.array([l for (l,h) in cl_bounds])
    high_bounds = np.array([h for (l,h) in cl_bounds])

    largest_cl_bound = (np.min(low_bounds), np.max(high_bounds))
    smallest_cl_bound = (np.max(low_bounds), np.min(high_bounds))

    if smallest_cl_bound[1] < smallest_cl_bound[0]:
        # Then the smallest CL is NULL
        smallest_cl_bound = (0.0, 0.0)

    ci_uncertainty = _cl_width(largest_cl_bound) - _cl_width(smallest_cl_bound)

    relative_change = ci_uncertainty/param_mean

    frac_uncertainty = ci_uncertainty/_cl_width(all_cl_bound)

    quant_uncertainty = float(_cl_count(largest_cl_bound, all_samples) - _cl_count(smallest_cl_bound, all_samples))/float(N)

    return (relative_change, frac_uncertainty, quant_uncertainty)


def plot_waveform(pos=None,siminspiral=None,event=0,path=None,ifos=['H1','L1','V1']):
    #import sim inspiral table content handler
    from ligo.lw import lsctables,ligolw
    from lalsimulation.lalsimulation import SimInspiralChooseTDWaveform,SimInspiralChooseFDWaveform
    from lalsimulation.lalsimulation import SimInspiralImplementedTDApproximants,SimInspiralImplementedFDApproximants
    from lal.lal import CreateREAL8TimeSeries,CreateForwardREAL8FFTPlan,CreateTukeyREAL8Window,CreateCOMPLEX16FrequencySeries,DimensionlessUnit,REAL8TimeFreqFFT
    from lal.lal import ComputeDetAMResponse, GreenwichMeanSiderealTime
    from lal.lal import LIGOTimeGPS
    from lal.lal import MSUN_SI as LAL_MSUN_SI
    from lal.lal import PC_SI as LAL_PC_SI
    import lalsimulation as lalsim
    from math import cos,sin,sqrt
    from ligo.lw import utils
    import os
    import numpy as np
    from numpy import arange
    if path is None:
        path=os.getcwd()
    if event is None:
        event=0
    colors_inj={'H1':'r','L1':'g','V1':'m','I1':'b','J1':'y'}
    colors_rec={'H1':'k','L1':'k','V1':'k','I1':'k','J1':'k'}
    # time and freq data handling variables
    srate=4096.0
    seglen=60.
    length=srate*seglen # lenght of 60 secs, hardcoded. May call a LALSimRoutine to get an idea
    deltaT=1/srate
    deltaF = 1.0 / (length* deltaT);

    # build window for FFT
    pad=0.4
    timeToFreqFFTPlan = CreateForwardREAL8FFTPlan(int(length), 1 );
    window=CreateTukeyREAL8Window(int(length),2.0*pad*srate/length);
    WinNorm = sqrt(window.sumofsquares/window.data.length);
    # time and freq domain strain:
    segStart=100000000
    strainT=CreateREAL8TimeSeries("strainT",segStart,0.0,1.0/srate,DimensionlessUnit,int(length));
    strainF= CreateCOMPLEX16FrequencySeries("strainF",segStart,   0.0,    deltaF, DimensionlessUnit,int(length/2. +1));

    f_min=25 # hardcoded default (may be changed below)
    f_ref=100 # hardcoded default (may be changed below)
    f_max=srate/2.0
    plot_fmax=f_max

    inj_strains=dict((i,{"T":{'x':None,'strain':None},"F":{'x':None,'strain':None}}) for i in ifos)
    rec_strains=dict((i,{"T":{'x':None,'strain':None},"F":{'x':None,'strain':None}}) for i in ifos)

    inj_domain=None
    rec_domain=None
    font_size=26
    if siminspiral is not None:
        skip=0
        try:
            xmldoc = utils.load_filename(siminspiral,contenthandler=lsctables.use_in(ligolw.LIGOLWContentHandler))
            tbl = lsctables.table.get_table(xmldoc, "sim_inspiral")
            if event>0:
                tbl=tbl[event]
            else:
                tbl=tbl[0]
        except:
            e = sys.exc_info()[0]
            print(e)
            print("Cannot read event %s from table %s. Won't plot injected waveform \n"%(event,siminspiral))
            skip=1
        if not skip:
            REAL8time=tbl.geocent_end_time+1e-9*tbl.geocent_end_time_ns
            GPStime=LIGOTimeGPS(REAL8time)
            M1=tbl.mass1
            M2=tbl.mass2
            D=tbl.distance
            m1=M1*LAL_MSUN_SI
            m2=M2*LAL_MSUN_SI
            phiRef=tbl.coa_phase

            f_min = tbl.f_lower
            s1x = tbl.spin1x
            s1y = tbl.spin1y
            s1z = tbl.spin1z
            s2x = tbl.spin2x
            s2y = tbl.spin2y
            s2z = tbl.spin2z

            r=D*LAL_PC_SI*1.0e6
            iota=tbl.inclination
            print("WARNING: Defaulting to inj_fref =100Hz to plot the injected WF. This is hardcoded since xml table does not carry this information\n")

            lambda1=0
            lambda2=0
            wf=str(tbl.waveform)

            injapproximant=lalsim.GetApproximantFromString(wf)
            amplitudeO=int(tbl.amp_order )
            phaseO=lalsim.GetOrderFromString(wf)

            waveFlags=lal.CreateDict()
            lalsim.SimInspiralWaveformParamsInsertPNAmplitudeOrder(waveFlags, amplitudeO)
            lalsim.SimInspiralWaveformParamsInsertPNPhaseOrder(waveFlags, phaseO)

            ra=tbl.longitude
            dec=tbl.latitude
            psi=tbl.polarization

            if SimInspiralImplementedFDApproximants(injapproximant):
                inj_domain='F'
                [plus,cross]=SimInspiralChooseFDWaveform(m1, m2, s1x, s1y, s1z,s2x,s2y,s2z,r, iota, phiRef,
                        0, 0, 0, # Non-circular binary parameters
                        deltaF, f_min, f_max, f_ref,
                        waveFlags, injapproximant)
            elif SimInspiralImplementedTDApproximants(injapproximant):
                inj_domain='T'
                [plus,cross]=SimInspiralChooseTDWaveform(m1, m2, s1x, s1y, s1z,s2x,s2y,s2z, r, iota, phiRef,
                        0, 0, 0, # Non-circular binary parameters
                        deltaT, f_min, f_ref,
                        waveFlags, injapproximant)
            else:
                print("\nThe approximant %s doesn't seem to be recognized by lalsimulation!\n Skipping WF plots\n"%injapproximant)
                return None

            for ifo in ifos:
                fp, fc = ComputeDetAMResponse(lal.cached_detector_by_prefix[ifo].response, ra, dec, psi, GreenwichMeanSiderealTime(REAL8time))
                if inj_domain=='T':
                    # strain is a temporary container for this IFO strain.
                    # Take antenna pattern into accout and window the data
                    for k in np.arange(strainT.data.length):
                        if k<plus.data.length:
                            strainT.data.data[k]=((fp*plus.data.data[k]+fc*cross.data.data[k]))
                        else:
                            strainT.data.data[k]=0.0
                        strainT.data.data[k]*=window.data.data[k]
                    # now copy in the dictionary only the part of strain which is not null (that is achieved using plus.data.length as length)
                    inj_strains[ifo]["T"]['strain']=np.array([strainT.data.data[k] for k in arange(plus.data.length)])
                    inj_strains[ifo]["T"]['x']=np.array([REAL8time - deltaT*(plus.data.length-1-k) for k in np.arange(plus.data.length)])

                    # Take the FFT
                    for j in arange(strainF.data.length):
                        strainF.data.data[j]=0.0
                    REAL8TimeFreqFFT(strainF,strainT,timeToFreqFFTPlan);
                    for j in arange(strainF.data.length):
                        strainF.data.data[j]/=WinNorm
                    # copy in the dictionary
                    inj_strains[ifo]["F"]['strain']=np.array([strainF.data.data[k] for k in arange(int(strainF.data.length))])
                    inj_strains[ifo]["F"]['x']=np.array([strainF.f0+ k*strainF.deltaF for k in arange(int(strainF.data.length))])
                elif inj_domain=='F':
                    for k in np.arange(strainF.data.length):
                        if k<plus.data.length:
                            strainF.data.data[k]=((fp*plus.data.data[k]+fc*cross.data.data[k]))
                        else:
                            strainF.data.data[k]=0.0
                    # copy in the dictionary
                    inj_strains[ifo]["F"]['strain']=np.array([strainF.data.data[k] for k in arange(int(strainF.data.length))])
                    inj_strains[ifo]["F"]['x']=np.array([strainF.f0+ k*strainF.deltaF for k in arange(int(strainF.data.length))])
    if pos is not None:

        # Select the maxP sample
        _,which=pos._posMap()

        if 'time' in pos.names:
            REAL8time=pos['time'].samples[which][0]
        elif 'time_maxl' in pos.names:
            REAL8time=pos['time_maxl'].samples[which][0]
        elif 'time_min' in pos.names and 'time_max' in pos.names:
            REAL8time=pos['time_min'].samples[which][0]+0.5*(pos['time_max'].samples[which][0]-pos['time_min'].samples[which][0])
        else:
            print("ERROR: could not find any time parameter in the posterior file. Not plotting the WF...\n")
            return None

        # first check we have approx in posterior samples, otherwise skip
        skip=0
        try:
            approximant=int(pos['LAL_APPROXIMANT'].samples[which][0])
            amplitudeO=int(pos['LAL_AMPORDER'].samples[which][0])
            phaseO=int(pos['LAL_PNORDER'].samples[which][0])
        except:
            skip=1
        if skip==0:
            GPStime=LIGOTimeGPS(REAL8time)

            q=pos['q'].samples[which][0]
            mc=pos['mc'].samples[which][0]
            M1,M2=q2ms(mc,q)
            if 'dist' in pos.names:
                D=pos['dist'].samples[which][0]
            elif 'distance' in pos.names:
                D=pos['distance'].samples[which][0]
            elif 'logdistance' in pos.names:
                D=exp(pos['distance'].samples[which][0])

            m1=M1*LAL_MSUN_SI
            m2=M2*LAL_MSUN_SI
            if 'phi_orb' in pos.names:
                phiRef=pos['phi_orb'].samples[which][0]
            elif 'phase' in pos.names:
                phiRef=pos['phase'].samples[which][0]
            elif 'phase_maxl' in pos.names:
                phiRef=pos['phase_maxl'].samples[which][0]
                print('INFO: phi_orb not estimated, using maximum likelihood value')
            else:
                print('WARNING: phi_orb not found in posterior files. Defaulting to 0.0 which is probably *not* what you want\n')
                phiRef=0.0

            try:
                for name in ['flow','f_lower']:
                    if name in pos.names:
                        f_min=pos[name].samples[which][0]
            except:
                pass

            try:
                for name in ['fref','f_ref','f_Ref','fRef']:
                    if name in pos.names:
                        fname=name

                Fref = np.unique(pos[fname].samples)
                if len(Fref) > 1:
                    print("ERROR: Expected f_ref to be constant for all samples.  Can't tell which value was injected! Defaulting to 100 Hz\n")
                    print(Fref)
                else:
                    f_ref = Fref[0]
            except ValueError:
                print("WARNING: Could not read fref from posterior file! Defaulting to 100 Hz\n")

            try:
                a = pos['a1'].samples[which][0]
                the = pos['theta_spin1'].samples[which][0]
                phi = pos['phi_spin1'].samples[which][0]
                s1x = (a * sin(the) * cos(phi));
                s1y = (a * sin(the) * sin(phi));
                s1z = (a * cos(the));
                a = pos['a2'].samples[which][0]
                the = pos['theta_spin2'].samples[which][0]
                phi = pos['phi_spin2'].samples[which][0]
                s2x = (a * sin(the) * cos(phi));
                s2y = (a * sin(the) * sin(phi));
                s2z = (a * cos(the));
                iota=pos['inclination'].samples[which][0]
            except:
                try:
                    iota, s1x, s1y, s1z, s2x, s2y, s2z=lalsim.SimInspiralTransformPrecessingNewInitialConditions(pos['theta_jn'].samples[which][0], pos['phi_JL'].samples[which][0], pos['tilt1'].samples[which][0], pos['tilt2'].samples[which][0], pos['phi12'].samples[which][0], pos['a1'].samples[which][0], pos['a2'].samples[which][0], m1, m2, f_ref, phiRef)
                except:
                    if 'a1z' in pos.names:
                        s1z=pos['a1z'].samples[which][0]
                    elif 'a1' in pos.names:
                        s1z=pos['a1'].samples[which][0]
                    else:
                        s1z=0
                    if 'a2z' in pos.names:
                        s2z=pos['a2z'].samples[which][0]
                    elif 'a2' in pos.names:
                        s2z=pos['a2'].samples[which][0]
                    else:
                        s2z=0
                    s1x=s1y=s2x=s2y=0.0
                    if 'inclination' in pos.names:
                        iota=pos['inclination'].samples[which][0]
                    else:
                        iota=pos['theta_jn'].samples[which][0]

            r=D*LAL_PC_SI*1.0e6

            approximant=int(pos['LAL_APPROXIMANT'].samples[which][0])
            amplitudeO=int(pos['LAL_AMPORDER'].samples[which][0])
            phaseO=int(pos['LAL_PNORDER'].samples[which][0])

            waveFlags=lal.CreateDict()
            lalsim.SimInspiralWaveformParamsInsertPNAmplitudeOrder(waveFlags, amplitudeO)
            lalsim.SimInspiralWaveformParamsInsertPNPhaseOrder(waveFlags, phaseO)
            if 'tideO' in pos.names:
                tidalO=int(pos['tideO'].samples[which][0])
                lalsim.SimInspiralWaveformParamsInsertPNTidalOrder(waveFlags, tidalO)
            if 'spinO' in pos.names:
                spinO=int(pos['spinO'].samples[which][0])
                lalsim.SimInspiralWaveformParamsInsertPNSpinOrder(waveFlags, spinO)
            if 'lambda1' in pos.names:
                lalsim.SimInspiralWaveformParamsInsertTidalLambda1(waveFlags, pos['lambda1'].samples[which][0])
            if 'lambda2' in pos.names:
                lalsim.SimInspiralWaveformParamsInsertTidalLambda2(waveFlags, pos['lambda2'].samples[which][0])

            if SimInspiralImplementedFDApproximants(approximant):
                rec_domain='F'
                [plus,cross]=SimInspiralChooseFDWaveform(m1, m2, s1x, s1y, s1z,s2x,s2y,s2z,r, iota, phiRef,
                        0, 0, 0, # Non-circular binary parameters
                        deltaF, f_min, f_max, f_ref,
                        waveFlags, approximant)
            elif SimInspiralImplementedTDApproximants(approximant):
                rec_domain='T'
                [plus,cross]=SimInspiralChooseTDWaveform(m1, m2, s1x, s1y, s1z,s2x,s2y,s2z, r, iota, phiRef,
                        0, 0, 0, # Non-circular binary parameters
                        deltaT, f_min, f_ref,
                        waveFlags, approximant)
            else:
                print("The approximant %s doesn't seem to be recognized by lalsimulation!\n Skipping WF plots\n"%approximant)
                return None

            ra=pos['ra'].samples[which][0]
            dec=pos['dec'].samples[which][0]
            psi=pos['psi'].samples[which][0]
            fs={}
            for ifo in ifos:
                fp, fc = ComputeDetAMResponse(lal.cached_detector_by_prefix[ifo].response, ra, dec, psi, GreenwichMeanSiderealTime(REAL8time))
                if rec_domain=='T':
                # strain is a temporary container for this IFO strain.
                # Take antenna pattern into accout and window the data
                    for k in np.arange(strainT.data.length):
                        if k<plus.data.length:
                            strainT.data.data[k]=((fp*plus.data.data[k]+fc*cross.data.data[k]))
                        else:
                            strainT.data.data[k]=0.0
                        strainT.data.data[k]*=window.data.data[k]
                    # now copy in the dictionary only the part of strain which is not null (that is achieved using plus.data.length as length)
                    rec_strains[ifo]["T"]['strain']=np.array([strainT.data.data[k] for k in arange(plus.data.length)])
                    rec_strains[ifo]["T"]['x']=np.array([REAL8time - deltaT*(plus.data.length-1-k) for k in np.arange(plus.data.length)])

                    # Take the FFT
                    for j in arange(strainF.data.length):
                        strainF.data.data[j]=0.0
                    REAL8TimeFreqFFT(strainF,strainT,timeToFreqFFTPlan);
                    for j in arange(strainF.data.length):
                        strainF.data.data[j]/=WinNorm
                    # copy in the dictionary
                    rec_strains[ifo]["F"]['strain']=np.array([strainF.data.data[k] for k in arange(int(strainF.data.length))])
                    rec_strains[ifo]["F"]['x']=np.array([strainF.f0+ k*strainF.deltaF for k in arange(int(strainF.data.length))])
                elif rec_domain=='F':
                    for k in np.arange(strainF.data.length):
                        if k<plus.data.length:
                            strainF.data.data[k]=((fp*plus.data.data[k]+fc*cross.data.data[k]))
                        else:
                            strainF.data.data[k]=0.0
                    # copy in the dictionary
                    rec_strains[ifo]["F"]['strain']=np.array([strainF.data.data[k] for k in arange(int(strainF.data.length))])
                    rec_strains[ifo]["F"]['x']=np.array([strainF.f0+ k*strainF.deltaF for k in arange(int(strainF.data.length))])

    myfig=plt.figure(1,figsize=(23,15))

    rows=len(ifos)
    cols=2

    #this variables decide which domain will be plotted on the left column of the plot.
    # only plot Time domain if both injections and recovery are TD
    global_domain="F"
    if rec_domain is not None and inj_domain is not None:
        if rec_domain=="T" and inj_domain=="T":
            global_domain="T"
    elif rec_domain is not None:
        if rec_domain=="T":
            global_domain="T"
    elif inj_domain is not None:
        if inj_domain=="T":
            global_domain="T"

    A,axes=plt.subplots(nrows=rows,ncols=cols,sharex=False,sharey=False)
    plt.setp(A,figwidth=23,figheight=15)
    for (r,i) in zip(np.arange(rows),ifos):
        for c in np.arange(cols):
            ax=axes[r]
            if type(ax)==np.ndarray:
                ax=ax[c]
            else:
                ax=axes[c]
            if rec_strains[i]["T"]['strain'] is not None or rec_strains[i]["F"]['strain'] is not None:
                if c==0:
                    if global_domain=="T":
                        ax.plot(rec_strains[i]["T"]['x'],rec_strains[i]["T"]['strain'],colors_rec[i],alpha=0.5,label='%s maP'%i)
                    else:
                        data=rec_strains[i]["F"]['strain']
                        f=rec_strains[i]["F"]['x']
                        mask=np.logical_and(f>=f_min,f<=plot_fmax)
                        ys=data
                        ax.semilogx(f[mask],ys[mask].real,'.-',color=colors_rec[i],alpha=0.5,label='%s maP'%i)
                else:
                    data=rec_strains[i]["F"]['strain']
                    f=rec_strains[i]["F"]['x']
                    mask=np.logical_and(f>=f_min,f<=plot_fmax)
                    ys=data
                    ax.loglog(f[mask],abs(ys[mask]),'--',color=colors_rec[i],alpha=0.5,linewidth=4)
                    ax.set_xlim([min(f[mask]),max(f[mask])])
                    ax.grid(True,which='both')
            if inj_strains[i]["T"]['strain'] is not None or inj_strains[i]["F"]['strain'] is not None:
                if c==0:
                    if global_domain=="T":
                        ax.plot(inj_strains[i]["T"]['x'],inj_strains[i]["T"]['strain'],colors_inj[i],alpha=0.5,label='%s inj'%i)
                    else:
                        data=inj_strains[i]["F"]['strain']
                        f=inj_strains[i]["F"]['x']
                        mask=np.logical_and(f>=f_min,f<=plot_fmax)
                        ys=data
                        ax.plot(f[mask],ys[mask].real,'.-',color=colors_inj[i],alpha=0.5,label='%s inj'%i)
                else:
                    data=inj_strains[i]["F"]['strain']
                    f=inj_strains[i]["F"]['x']
                    mask=np.logical_and(f>=f_min,f<=plot_fmax)
                    ys=data
                    ax.loglog(f[mask],abs(ys[mask]),'--',color=colors_inj[i],alpha=0.5,linewidth=4)
                    ax.set_xlim([min(f[mask]),max(f[mask])])
                    ax.grid(True,which='both')

            if r==0:
                if c==0:
                    if global_domain=="T":
                        ax.set_title(r"$h(t)$",fontsize=font_size)
                    else:
                        ax.set_title(r"$\Re[h(f)]$",fontsize=font_size)
                else:
                    ax.set_title(r"$|h(f)|$",fontsize=font_size)
            elif r==rows-1:
                if c==0:
                    if global_domain=="T":
                        ax.set_xlabel("time [s]",fontsize=font_size)
                    else:
                        ax.set_xlabel("frequency [Hz]",fontsize=font_size)
                else:
                    ax.set_xlabel("frequency [Hz]",fontsize=font_size)

            ax.legend(loc='best')
            ax.grid(True)

            #ax.tight_layout()
    A.savefig(os.path.join(path,'WF_DetFrame.png'),bbox_inches='tight')
    return inj_strains,rec_strains


def plot_psd(psd_files,outpath=None,f_min=30.):
    myfig2=plt.figure(figsize=(15,15),dpi=500)
    ax=plt.subplot(1,1,1)
    colors={'H1':'r','L1':'g','V1':'m','I1':'k','J1':'y'}

    if outpath is None:
        outpath=os.getcwd()
    tmp=[]
    for f in psd_files:
        if not os.path.isfile(f):
            print("PSD file %s has not been found and won't be plotted\n"%f)
        else:
            tmp.append(f)
    if tmp==[]:
        return None
    else:
        psd_files=tmp

    freqs = {}
    for f in psd_files:
        data=np.loadtxt(f)
        freq=data[:,0]
        data=data[:,1]
        idx=f.find('-PSD.dat')
        ifo=f[idx-2:idx]
        freqs[ifo.lower()] = freq
        fr=[]
        da=[]
        for (f,d) in zip(freq,data):
            if f>f_min and d!=0.0 and np.isfinite(d):
                fr.append(f)
                da.append(d)
        plt.loglog(fr,da,colors[ifo],label=ifo,alpha=0.5,linewidth=3)
    plt.xlim([min(fr),max(fr)])
    plt.xlabel("Frequency [Hz]",fontsize=26)
    plt.ylabel("PSD",fontsize=26)
    plt.legend(loc='best')
    plt.grid(which='both')
    try:
        plt.tight_layout()
        myfig2.savefig(os.path.join(outpath,'PSD.png'),bbox_inches='tight')
    except:
        myfig2.savefig(os.path.join(outpath,'PSD.png'))
    myfig2.clf()

    return freqs

cred_level = lambda cl, x: np.sort(x, axis=0)[int(cl*len(x))]

def cred_interval(x, cl=.9, lower=True):
    """Return location of lower or upper confidence levels
    Args:
        x: List of samples.
        cl: Confidence level to return the bound of.
        lower: If ``True``, return the lower bound, otherwise return the upper bound.
    """
    if lower:
        return cred_level((1.-cl)/2, x)
    else:
        return cred_level((1.+cl)/2, x)

def spline_angle_xform(delta_psi):
    """Returns the angle in degrees corresponding to the spline
    calibration parameters delta_psi.

    """
    rot = (2.0 + 1.0j*delta_psi)/(2.0 - 1.0j*delta_psi)

    return 180.0/np.pi*np.arctan2(np.imag(rot), np.real(rot))

def plot_spline_pos(logf, ys, nf=100, level=0.9, color='k', label=None, xform=None):
    """Plot calibration posterior estimates for a spline model in log space.
    Args:
        logf: The (log) location of spline control points.
        ys: List of posterior draws of function at control points ``logf``
        nx: Number of points to evaluate spline at for plotting.
        level: Credible level to fill in.
        color: Color to plot with.
        label: Label for plot.
        xform: Function to transform the spline into plotted values.
    """
    f = np.exp(logf)
    fs = np.linspace(f.min(), f.max(), nf)

    data = np.zeros((ys.shape[0], nf))

    if xform is None:
        zs = ys
    else:
        zs = xform(ys)

    mu = np.mean(zs, axis=0)
    lower_cl = mu - cred_interval(zs, level, lower=True)
    upper_cl = cred_interval(zs, level, lower=False) - mu
    plt.errorbar(np.exp(logf), mu, yerr=[lower_cl, upper_cl], fmt='.', color=color, lw=4, alpha=0.5, capsize=0)

    for i, samp in enumerate(ys):
        try:
            temp = interpolate.spline(logf, samp, np.log(fs))
        except AttributeError:   # scipy < 0.19.0
            calSpline = interpolate.InterpolatedUnivariateSpline(logf, samp, k=3, ext=2) #cubic spline (k=3), raises ValueError in extrapolation
            temp = calSpline(np.log(fs))
        if xform is None:
            data[i] = temp
        else:
            data[i] = xform(temp)

    line, = plt.plot(fs, np.mean(data, axis=0), color=color, label=label)
    color = line.get_color()
    plt.fill_between(fs, cred_interval(data, level), cred_interval(data, level, lower=False), color=color, alpha=.1, linewidth=0.1)
    plt.xlim(f.min()-.5, f.max()+50)

def plot_calibration_pos(pos, level=.9, outpath=None):
    fig, [ax1, ax2] = plt.subplots(2, 1, figsize=(15, 15), dpi=500)

    font_size = 32
    if outpath is None:
        outpath=os.getcwd()

    params = pos.names
    ifos = np.unique([param.split('_')[0] for param in params if 'spcal_freq' in param])
    for ifo in ifos:
        if ifo=='h1': color = 'r'
        elif ifo=='l1': color = 'g'
        elif ifo=='v1': color = 'm'
        else: color = 'c'

        # Assume spline control frequencies are constant
        freq_params = np.sort([param for param in params if
                               '{0}_spcal_freq'.format(ifo) in param])

        logfreqs = np.log([pos[param].median for param in freq_params])

        # Amplitude calibration model
        plt.sca(ax1)
        amp_params = np.sort([param for param in params if
                              '{0}_spcal_amp'.format(ifo) in param])
        if len(amp_params) > 0:
            amp = 100*np.column_stack([pos[param].samples for param in amp_params])
            plot_spline_pos(logfreqs, amp, color=color, level=level, label="{0} (mean, {1}%)".format(ifo.upper(), int(level*100)))

        # Phase calibration model
        plt.sca(ax2)
        phase_params = np.sort([param for param in params if
                                '{0}_spcal_phase'.format(ifo) in param])
        if len(phase_params) > 0:
            phase = np.column_stack([pos[param].samples for param in phase_params])

            plot_spline_pos(logfreqs, phase, color=color, level=level, label="{0} (mean, {1}%)".format(ifo.upper(), int(level*100)), xform=spline_angle_xform)

    ax1.tick_params(labelsize=.75*font_size)
    ax2.tick_params(labelsize=.75*font_size)
    try:
        plt.legend(loc='upper right', prop={'size':.75*font_size}, framealpha=0.1)
    except:
        plt.legend(loc='upper right', prop={'size':.75*font_size})
    ax1.set_xscale('log')
    ax2.set_xscale('log')

    ax2.set_xlabel('Frequency (Hz)', fontsize=font_size)
    ax1.set_ylabel('Amplitude (%)', fontsize=font_size)
    ax2.set_ylabel('Phase (deg)', fontsize=font_size)

    outp = os.path.join(outpath, 'calibration.png')
    try:
        fig.tight_layout()
        fig.savefig(outp, bbox_inches='tight')
    except:
        fig.savefig(outp)
    plt.close(fig)


def plot_burst_waveform(pos=None,simburst=None,event=0,path=None,ifos=['H1','L1','V1']):
    from lalinference.lalinference import SimBurstChooseFDWaveform,SimBurstChooseTDWaveform
    from lalinference.lalinference import SimBurstImplementedFDApproximants,SimBurstImplementedTDApproximants
    from lal.lal import CreateREAL8TimeSeries,CreateForwardREAL8FFTPlan,CreateTukeyREAL8Window,CreateCOMPLEX16FrequencySeries,DimensionlessUnit,REAL8TimeFreqFFT,CreateReverseREAL8FFTPlan
    from lal.lal import LIGOTimeGPS
    import lalinference as lalinf
    from lal import ComputeDetAMResponse, GreenwichMeanSiderealTime, LIGOTimeGPS

    from math import cos,sin,sqrt
    from ligo.lw import lsctables
    from ligo.lw import utils
    import os
    import numpy as np
    from numpy import arange,real,absolute,fabs,pi
    from matplotlib import pyplot as plt
    if path is None:
        path=os.getcwd()
    if event is None:
        event=0
    colors_inj={'H1':'r','L1':'g','V1':'m','I1':'b','J1':'y'}
    colors_rec={'H1':'k','L1':'k','V1':'k','I1':'k','J1':'k'}
    #import sim inspiral table content handler
    from ligo.lw import ligolw
    from ligo.lw import table
    class LIGOLWContentHandlerExtractSimBurstTable(ligolw.LIGOLWContentHandler):
        def __init__(self,document):
            ligolw.LIGOLWContentHandler.__init__(self,document)
            self.tabname=lsctables.SimBurstTable.tableName
            self.intable=False
            self.tableElementName=''
        def startElement(self,name,attrs):
            if attrs.has_key('Name') and table.Table.TableName(attrs['Name'])==self.tabname:
                self.tableElementName=name
                # Got the right table, let's see if it's the right event
                ligolw.LIGOLWContentHandler.startElement(self,name,attrs)
                self.intable=True
            elif self.intable: # We are in the correct table
                ligolw.LIGOLWContentHandler.startElement(self,name,attrs)
        def endElement(self,name):
            if self.intable: ligolw.LIGOLWContentHandler.endElement(self,name)
            if self.intable and name==self.tableElementName: self.intable=False

    lsctables.use_in(LIGOLWContentHandlerExtractSimBurstTable)

    # time and freq data handling variables
    srate=4096.0
    seglen=10.
    length=srate*seglen # lenght of 10 secs, hardcoded.
    deltaT=1/srate
    deltaF = 1.0 / (length* deltaT);

    # build window for FFT
    pad=0.4
    timeToFreqFFTPlan = CreateForwardREAL8FFTPlan(int(length), 1 );
    freqToTimeFFTPlan = CreateReverseREAL8FFTPlan(int(length), 1 );
    window=CreateTukeyREAL8Window(int(length),2.0*pad*srate/length);
    # A random GPS time to initialize arrays. Epoch will be overwritten with sensible times further down
    segStart=939936910.000
    strainFinj= CreateCOMPLEX16FrequencySeries("strainF",segStart,0.0,deltaF,DimensionlessUnit,int(length/2. +1));
    strainTinj=CreateREAL8TimeSeries("strainT",segStart,0.0,1.0/srate,DimensionlessUnit,int(length));
    strainFrec= CreateCOMPLEX16FrequencySeries("strainF",segStart,0.0,deltaF,DimensionlessUnit,int(length/2. +1));
    strainTrec=CreateREAL8TimeSeries("strainT",segStart,0.0,1.0/srate,DimensionlessUnit,int(length));
    GlobREAL8time=None
    f_min=25 # hardcoded default (may be changed below)
    f_ref=100 # hardcoded default (may be changed below)
    f_max=srate/2.0

    plot_fmax=2048
    plot_fmin=0.01
    plot_tmin=1e11
    plot_tmax=-1e11

    inj_strains=dict((i,{"T":{'x':None,'strain':None},"F":{'x':None,'strain':None}}) for i in ifos)
    rec_strains=dict((i,{"T":{'x':None,'strain':None},"F":{'x':None,'strain':None}}) for i in ifos)

    inj_domain=None
    rec_domain=None
    font_size=26
    if simburst is not None:
        skip=0
        try:
            xmldoc = utils.load_filename(simburst,contenthandler=LIGOLWContentHandlerExtractSimBurstTable)
            tbl = lsctables.SimBurstTable.get_table(xmldoc)
            if event>0:
                tbl=tbl[event]
            else:
                tbl=tbl[0]
        except:
            print("Cannot read event %s from table %s. Won't plot injected waveform \n"%(event,simburst))
            skip=1
        if not skip:
            REAL8time=tbl.time_geocent_gps+1e-9*tbl.time_geocent_gps_ns
            GPStime=LIGOTimeGPS(REAL8time)
            GlobREAL8time=(REAL8time)
            strainTinj.epoch=LIGOTimeGPS(round(GlobREAL8time,0)-seglen/2.)
            strainFinj.epoch=LIGOTimeGPS(round(GlobREAL8time,0)-seglen/2.)
            f0=tbl.frequency
            q=tbl.q
            dur=tbl.duration
            hrss=tbl.hrss
            polar_e_angle=tbl.pol_ellipse_angle
            polar_e_ecc=tbl.pol_ellipse_e

            BurstExtraParams=None
            wf=str(tbl.waveform)

            injapproximant=lalinf.GetBurstApproximantFromString(wf)
            ra=tbl.ra
            dec=tbl.dec
            psi=tbl.psi

            if SimBurstImplementedFDApproximants(injapproximant):
                inj_domain='F'
                [plus,cross]=SimBurstChooseFDWaveform(deltaF, deltaT, f0, q,dur, f_min, f_max,hrss,polar_e_angle ,polar_e_ecc,BurstExtraParams, injapproximant)
            elif SimBurstImplementedTDApproximants(injapproximant):
                inj_domain='T'
                [plus,cross]=SimBurstChooseTDWaveform(deltaT, f0, q,dur, f_min, f_max,hrss,polar_e_angle ,polar_e_ecc,BurstExtraParams, injapproximant)
            else:
                print("\nThe approximant %s doesn't seem to be recognized by lalinference!\n Skipping WF plots\n"%injapproximant)
                return None

            for ifo in ifos:
                fp, fc = ComputeDetAMResponse(lal.cached_detector_by_prefix[ifo].response, ra, dec, psi, GreenwichMeanSiderealTime(REAL8time))
                if inj_domain=='T':
                    # bin of ref time as seen in strainT
                    tCinstrain=np.floor(REAL8time-float(strainTinj.epoch))/deltaT
                    # bin of strainT where we need to start copying the WF over
                    #tSinstrain=floor(tCinstrain-float(plus.data.length)/2.)+1
                    tSinstrain=int(  (REAL8time-fabs(float(plus.epoch)) - fabs(float(strainTinj.epoch)))/deltaT)
                    rem=(REAL8time-fabs(float(plus.epoch)) - fabs(float(strainTinj.epoch)))/deltaT-tSinstrain
                    # strain is a temporary container for this IFO strain.
                    # Zero until tSinstrain
                    for k in np.arange(tSinstrain):
                        strainTinj.data.data[k]=0.0
                    # then copy plus/cross over
                    for k in np.arange(plus.data.length):
                        strainTinj.data.data[k+tSinstrain]=((fp*plus.data.data[k]+fc*cross.data.data[k]))
                    # Then zeros till the end (superfluous)
                    for k in np.arange(strainTinj.data.length- (tSinstrain +plus.data.length)):
                        strainTinj.data.data[k+tSinstrain+plus.data.length]=0.0
                    for k in np.arange(strainTinj.data.length):
                        strainTinj.data.data[k]*=window.data.data[k]
                    np.savetxt('file.out',zip(np.array([strainTinj.epoch + j*deltaT for j in arange(strainTinj.data.length)]),np.array([strainTinj.data.data[j] for j in arange(strainTinj.data.length)])))
                    # now copy in the dictionary
                    inj_strains[ifo]["T"]['strain']=np.array([strainTinj.data.data[j] for j in arange(strainTinj.data.length)])
                    inj_strains[ifo]["T"]['x']=np.array([strainTinj.epoch + j*deltaT for j in arange(strainTinj.data.length)])
                    # Take the FFT
                    for j in arange(strainFinj.data.length):
                        strainFinj.data.data[j]=0.0
                    REAL8TimeFreqFFT(strainFinj,strainTinj,timeToFreqFFTPlan);
                    twopit=2.*np.pi*(rem*deltaT)
                    for k in arange(strainFinj.data.length):
                        re = cos(twopit*deltaF*k)
                        im = -sin(twopit*deltaF*k)
                        strainFinj.data.data[k]*= (re + 1j*im);
                    # copy in the dictionary
                    inj_strains[ifo]["F"]['strain']=np.array([strainFinj.data.data[k] for k in arange(int(strainFinj.data.length))])
                    inj_strains[ifo]["F"]['x']=np.array([strainFinj.f0+ k*strainFinj.deltaF for k in arange(int(strainFinj.data.length))])
                elif inj_domain=='F':
                    for k in np.arange(strainFinj.data.length):
                        if k<plus.data.length:
                            strainFinj.data.data[k]=((fp*plus.data.data[k]+fc*cross.data.data[k]))
                        else:
                            strainFinj.data.data[k]=0.0
                    twopit=2.*np.pi*(REAL8time-float(strainFinj.epoch))
                    for k in arange(strainFinj.data.length):
                        re = cos(twopit*deltaF*k)
                        im = -sin(twopit*deltaF*k)
                        strainFinj.data.data[k]*= (re + 1j*im);
                    # copy in the dictionary
                    inj_strains[ifo]["F"]['strain']=np.array([strainFinj.data.data[k] for k in arange(int(strainFinj.data.length))])
                    inj_strains[ifo]["F"]['x']=np.array([strainFinj.f0+ k*strainFinj.deltaF for k in arange(int(strainFinj.data.length))])
                #update xlimits for plot, go 6 sigmas left and right of f0
                # This should work for SineGaussians
                if f0 is not None and f0 is not np.nan:
                    if q is not None and q is not np.nan:
                        sigmaF=f0/q
                        if f0-6.*sigmaF>plot_fmin:
                            plot_fmin=f0-6.*sigmaF
                        if f0+6.*sigmaF<plot_fmax:
                            plot_fmax=f0+6.*sigmaF
                        sigmaT=q/(2.*pi*f0)
                        if REAL8time-6.*sigmaT<plot_tmin:
                            plot_tmin=REAL8time-6.*sigmaT
                        if REAL8time+6.*sigmaT>plot_tmax:
                            plot_tmax=REAL8time+6.*sigmaT
                # Now handle gaussians. For gaussians f0 is nan (FD centered at f=0)
                if dur is not None and dur is not np.nan:
                    sigmaF=1./sqrt(2.)/pi/dur
                    if 0+6.*sigmaF<plot_fmax:
                        plot_fmax=0+6.*sigmaF
                    plot_fmin=0.0
                    sigmaT=dur/sqrt(2.)
                    if REAL8time-6.*sigmaT<plot_tmin:
                        plot_tmin=REAL8time-6.*sigmaT
                    if REAL8time+6.*sigmaT>plot_tmax:
                        plot_tmax=REAL8time+6.*sigmaT


    if pos is not None:

        # Select the maxP sample
        _,which=pos._posMap()

        if 'time' in pos.names:
            REAL8time=pos['time'].samples[which][0]
        elif 'time_maxl' in pos.names:
            REAL8time=pos['time_maxl'].samples[which][0]
        elif 'time_mean' in pos.names:
            REAL8time=pos['time_mean'].samples[which][0]
        elif 'time_min' in pos.names and 'time_max' in pos.names:
            REAL8time=pos['time_min'].samples[which][0]+0.5*(pos['time_max'].samples[which][0]-pos['time_min'].samples[which][0])
        else:
            print("ERROR: could not find any time parameter in the posterior file. Not plotting the WF...\n")
            return None

        # first check we have approx in posterior samples, otherwise skip
        skip=0

        try:
            approximant=int(pos['LAL_APPROXIMANT'].samples[which][0])
        except:
            skip=1
        if skip==0:
            GPStime=LIGOTimeGPS(REAL8time)
            if GlobREAL8time is None:
                GlobREAL8time=REAL8time
            strainTrec.epoch=LIGOTimeGPS(round(GlobREAL8time,0)-seglen/2.)
            strainFrec.epoch=LIGOTimeGPS(round(GlobREAL8time,0)-seglen/2.)
            if "duration" in pos.names:
                dur=pos["duration"].samples[which][0]
            else:
                dur=np.nan
            if "quality" in pos.names:
                q=pos['quality'].samples[which][0]
            else:
                q=np.nan
            if 'frequency' in pos.names:
                f0=pos['frequency'].samples[which][0]
            else:
                f0=np.nan
            try:
                hrss=pos['hrss'].samples[which][0]
            except:
                hrss=exp(pos['loghrss'].samples[which][0])
            if np.isnan(q) and not np.isnan(dur):
                q=sqrt(2)*pi*dur
            alpha=None
            if 'alpha' in pos.names:
                alpha=pos['alpha'].samples[which][0]
                polar_e_angle=alpha
                polar_e_ecc=pos['polar_eccentricity'].samples[which][0]
            elif 'polar_ellipse_angle' in pos.names:
                polar_e_angle=pos['polar_ellipse_angle'].samples[which][0]
                polar_e_ecc=pos['polar_eccentricity'].samples[which][0]

            BurstExtraParams=None
            #if alpha:
            #  BurstExtraParams=lalsim.SimBurstCreateExtraParam("alpha",alpha)

            if SimBurstImplementedFDApproximants(approximant):
                rec_domain='F'
                [plus,cross]=SimBurstChooseFDWaveform(deltaF, deltaT, f0, q,dur, f_min, f_max,hrss,polar_e_angle ,polar_e_ecc,BurstExtraParams, approximant)
            elif SimBurstImplementedTDApproximants(approximant):
                rec_domain='T'
                [plus,cross]=SimBurstChooseTDWaveform(deltaT, f0, q,dur, f_min, f_max,hrss,polar_e_angle ,polar_e_ecc,BurstExtraParams, approximant)
            else:
                print("The approximant %s doesn't seem to be recognized by lalinference!\n Skipping WF plots\n"%approximant)
                return None
            ra=pos['ra'].samples[which][0]
            dec=pos['dec'].samples[which][0]
            psi=pos['psi'].samples[which][0]
            fs={}
            for ifo in ifos:
                fp, fc = ComputeDetAMResponse(lal.cached_detector_by_prefix[ifo].response, ra, dec, psi, GreenwichMeanSiderealTime(REAL8time))
                if rec_domain=='T':
                    # bin of ref time as seen in strainT
                    tCinstrain=np.floor(REAL8time-float(strainTrec.epoch))/deltaT
                    # bin of strainT where we need to start copying the WF over
                    tSinstrain=int(  (REAL8time-fabs(float(plus.epoch)) - fabs(float(strainTrec.epoch)))/deltaT)
                    #tSinstrain=floor(tCinstrain-float(plus.data.length)/2.)+1
                    #reminder for fractions of bin, will be added back in the FD WF
                    rem=(REAL8time-fabs(float(plus.epoch)) - fabs(float(strainTrec.epoch)))/deltaT-tSinstrain

                    # strain is a temporary container for this IFO strain.
                    # Zero until tSinstrain
                    for k in np.arange(tSinstrain):
                        strainTrec.data.data[k]=0.0
                    # then copy plus/cross over
                    for k in np.arange(plus.data.length):
                        strainTrec.data.data[k+tSinstrain]=((fp*plus.data.data[k]+fc*cross.data.data[k]))
                    # Then zeros till the end (superfluous)
                    for k in np.arange(strainTrec.data.length- (tSinstrain +plus.data.length)):
                        strainTrec.data.data[k+tSinstrain+plus.data.length]=0.0
                    for k in np.arange(strainTrec.data.length):
                        strainTrec.data.data[k]*=window.data.data[k]
                    # now copy in the dictionary
                    rec_strains[ifo]["T"]['strain']=np.array([strainTrec.data.data[j] for j in arange(strainTrec.data.length)])
                    rec_strains[ifo]["T"]['x']=np.array([strainTrec.epoch + j*deltaT for j in arange(strainTrec.data.length)])
                    # Take the FFT
                    for j in arange(strainFrec.data.length):
                        strainFrec.data.data[j]=0.0
                    REAL8TimeFreqFFT(strainFrec,strainTrec,timeToFreqFFTPlan);
                    twopit=2.*np.pi*(rem*deltaT)
                    for k in arange(strainFrec.data.length):
                        re = cos(twopit*deltaF*k)
                        im = -sin(twopit*deltaF*k)
                        strainFrec.data.data[k]*= (re + 1j*im);
                    # copy in the dictionary
                    rec_strains[ifo]["F"]['strain']=np.array([strainFrec.data.data[k] for k in arange(int(strainFrec.data.length))])
                    rec_strains[ifo]["F"]['x']=np.array([strainFrec.f0+ k*strainFrec.deltaF for k in arange(int(strainFrec.data.length))])
                elif rec_domain=='F':
                    for k in np.arange(strainFrec.data.length):
                        if k<plus.data.length:
                            strainFrec.data.data[k]=((fp*plus.data.data[k]+fc*cross.data.data[k]))
                        else:
                            strainFrec.data.data[k]=0.0
                    twopit=2.*np.pi*(REAL8time-float(strainFrec.epoch))
                    for k in arange(strainFrec.data.length):
                        re = cos(twopit*deltaF*k)
                        im = -sin(twopit*deltaF*k)
                        strainFrec.data.data[k]*= (re + 1j*im);
                    # copy in the dictionary
                    rec_strains[ifo]["F"]['strain']=np.array([strainFrec.data.data[k] for k in arange(int(strainFrec.data.length))])
                    rec_strains[ifo]["F"]['x']=np.array([strainFrec.f0+ k*strainFrec.deltaF for k in arange(int(strainFrec.data.length))])
                #update xlimits for plot, go 6 sigmas left and right of f0
                # This should work for SineGaussians
                if f0 is not None and f0 is not np.nan:
                    if q is not None and q is not np.nan:
                        sigmaF=f0/q
                        if f0-6.*sigmaF>plot_fmin:
                            plot_fmin=f0-6.*sigmaF
                        if f0+6.*sigmaF<plot_fmax:
                            plot_fmax=f0+6.*sigmaF
                        sigmaT=q/(2.*pi*f0)
                        if REAL8time-6.*sigmaT<plot_tmin:
                            plot_tmin=REAL8time-6.*sigmaT
                        if REAL8time+6.*sigmaT>plot_tmax:
                            plot_tmax=REAL8time+6.*sigmaT
                # Now handle gaussians. For gaussians f0 is nan (FD centered at f=0)
                if dur is not None and dur is not np.nan:
                    sigmaF=1./sqrt(2.)/pi/dur
                    if 0+6.*sigmaF<plot_fmax:
                        plot_fmax=0+6.*sigmaF
                    plot_fmin=0.0
                    sigmaT=dur/sqrt(2.)
                    if REAL8time-6.*sigmaT<plot_tmin:
                        plot_tmin=REAL8time-6.*sigmaT
                    if REAL8time+6.*sigmaT>plot_tmax:
                        plot_tmax=REAL8time+6.*sigmaT

    myfig=plt.figure(1,figsize=(10,7))

    rows=len(ifos)
    cols=2

    #this variables decide which domain will be plotted on the left column of the plot.
    # only plot Time domain if both injections and recovery are TD
    global_domain="F"
    if rec_domain is not None and inj_domain is not None:
        if rec_domain=="T" and inj_domain=="T":
            global_domain="T"
    elif rec_domain is not None:
        if rec_domain=="T":
            global_domain="T"
    elif inj_domain is not None:
        if inj_domain=="T":
            global_domain="T"

    A,axes=plt.subplots(nrows=rows,ncols=cols,sharex=False,sharey=False)
    plt.setp(A,figwidth=10,figheight=7)
    for (r,i) in zip(np.arange(rows),ifos):
        for c in np.arange(cols):
            ax=axes[r]
            if type(ax)==np.ndarray:
                ax=ax[c]
            else:
                ax=axes[c]
            if rec_strains[i]["T"]['strain'] is not None or rec_strains[i]["F"]['strain'] is not None:
                if c==0:
                    if global_domain=="T":
                        ax.plot(rec_strains[i]["T"]['x'],rec_strains[i]["T"]['strain'],colors_rec[i],label='%s maP'%i,linewidth=5)
                        ax.set_xlim([plot_tmin,plot_tmax])
                        #ax.vlines(GlobREAL8time,0.9*min(rec_strains[i]["T"]['strain']),1.1*max(rec_strains[i]["T"]['strain']),'k')
                    else:
                        data=rec_strains[i]["F"]['strain']
                        f=rec_strains[i]["F"]['x']
                        mask=np.logical_and(f>=plot_fmin,f<=plot_fmax)
                        ys=data
                        ax.plot(f[mask],real(ys[mask]),'-',color=colors_rec[i],label='%s maP'%i,linewidth=5)
                        ax.set_xlim([plot_fmin,plot_fmax])
                else:
                    data=rec_strains[i]["F"]['strain']
                    f=rec_strains[i]["F"]['x']
                    mask=np.logical_and(f>=plot_fmin,f<=plot_fmax)
                    ys=data
                    ax.loglog(f[mask],absolute(ys[mask]),'--',color=colors_rec[i],linewidth=5)
                    ax.grid(True,which='both')
                    ax.set_xlim([plot_fmin,plot_fmax])
            if inj_strains[i]["T"]['strain'] is not None or inj_strains[i]["F"]['strain'] is not None:
                if c==0:
                    if global_domain=="T":
                        ax.plot(inj_strains[i]["T"]['x'],inj_strains[i]["T"]['strain'],colors_inj[i],label='%s inj'%i,linewidth=2)
                        ax.set_xlim([plot_tmin,plot_tmax])
                    else:
                        data=inj_strains[i]["F"]['strain']
                        f=inj_strains[i]["F"]['x']
                        mask=np.logical_and(f>=plot_fmin,f<=plot_fmax)
                        ys=data
                        ax.plot(f[mask],real(ys[mask]),'-',color=colors_inj[i],label='%s inj'%i,linewidth=2)
                        ax.set_xlim([plot_fmin,plot_fmax])
                else:
                    data=inj_strains[i]["F"]['strain']
                    f=inj_strains[i]["F"]['x']
                    mask=np.logical_and(f>=plot_fmin,f<=plot_fmax)
                    ys=data
                    ax.loglog(f[mask],absolute(ys[mask]),'--',color=colors_inj[i],linewidth=2)
                    ax.grid(True,which='both')
                    ax.set_xlim([plot_fmin,plot_fmax])
            if r==0:
                if c==0:
                    if global_domain=="T":
                        ax.set_title(r"$h(t)$",fontsize=font_size)
                    else:
                        ax.set_title(r"$\Re[h(f)]$",fontsize=font_size)
                else:
                    ax.set_title(r"$|h(f)|$",fontsize=font_size)
            elif r==rows-1:
                if c==0:
                    if global_domain=="T":
                        ax.set_xlabel("time [s]",fontsize=font_size)
                    else:
                        ax.set_xlabel("frequency [Hz]",fontsize=font_size)
                else:
                    ax.set_xlabel("frequency [Hz]",fontsize=font_size)

            ax.legend(loc='best')
            ax.grid(True)

            #ax.tight_layout()
    A.savefig(os.path.join(path,'WF_DetFrame.png'),bbox_inches='tight')
    return inj_strains, rec_strains

def make_1d_table(html,legend,label,pos,pars,noacf,GreedyRes,onepdfdir,sampsdir,savepdfs,greedy,analyticLikelihood,nDownsample):

    from numpy import unique, sort
    global confidenceLevels
    confidence_levels=confidenceLevels

    out={}
    if pars==[]:
        return out
    if set(pos.names)-set(pars)==set(pos.names):
        return out

    #2D plots list
    tabid='onedmargtable_'+label.lower()
    html_ompdf=html.add_collapse_section('1D marginal posterior PDFs (%s)'%label,legend=legend,innertable_id=tabid)
    #Table matter
    if not noacf:
        html_ompdf_write= '<table id="%s"><tr><th>Histogram and Kernel Density Estimate</th><th>Samples used</th><th>Autocorrelation</th></tr>'%tabid
    else:
        html_ompdf_write= '<table id="%s"><tr><th>Histogram and Kernel Density Estimate</th><th>Samples used</th></tr>'%tabid

    Nskip=0
    if 'chain' in pos.names:
        data,header=pos.samples()
        par_index=pos.names.index('cycle')
        chain_index=pos.names.index("chain")
        chains=unique(pos["chain"].samples)
        chainCycles = [sort(data[ data[:,chain_index] == chain, par_index ]) for chain in chains]
        chainNcycles = []
        chainNskips = []
        for cycles in chainCycles:
            if len(cycles) > 1:
                chainNcycles.append(cycles[-1] - cycles[0])
                chainNskips.append(cycles[1] - cycles[0])
            else:
                chainNcycles.append(1)
                chainNskips.append(1)
    elif 'cycle' in pos.names:
        cycles = sort(pos['cycle'].samples)
        if len(cycles) > 1:
            Ncycles = cycles[-1]-cycles[0]
            Nskip = cycles[1]-cycles[0]
        else:
            Ncycles = 1
            Nskip = 1

    printed=0
    for par_name in pars:
        par_name=par_name.lower()
        try:
            pos[par_name.lower()]
        except KeyError:
            #print "No input chain for %s, skipping binning."%par_name
            continue
        try:
            par_bin=GreedyRes[par_name]
        except KeyError:
            print("Bin size is not set for %s, skipping binning."%par_name)
            continue

        #print "Binning %s to determine confidence levels ..."%par_name
        binParams={par_name:par_bin}
        injection_area=None
        injection_area=None
        if greedy:
            if printed==0:
                print("Using greedy 1-d binning credible regions\n")
                printed=1
            toppoints,injectionconfidence,reses,injection_area,cl_intervals=greedy_bin_one_param(pos,binParams,confidence_levels)
        else:
            if printed==0:
                print("Using 2-step KDE 1-d credible regions\n")
                printed=1
            if pos[par_name].injval is None:
                injCoords=None
            else:
                injCoords=[pos[par_name].injval]
            _,reses,injstats=kdtree_bin2Step(pos,[par_name],confidence_levels,injCoords=injCoords)
            if injstats is not None:
                injectionconfidence=injstats[3]
                injection_area=injstats[4]
        #Generate 1D histogram/kde plots
        print("Generating 1D plot for %s."%par_name)
        out[par_name]=reses
        #Get analytic description if given
        pdf=cdf=None
        if analyticLikelihood:
            pdf = analyticLikelihood.pdf(par_name)
            cdf = analyticLikelihood.cdf(par_name)

        oneDPDFParams={par_name:50}
        try:
            rbins,plotFig=plot_one_param_pdf(pos,oneDPDFParams,pdf,cdf,plotkde=False)
        except:
            print("Failed to produce plot for %s."%par_name)
            continue

        figname=par_name+'.png'
        oneDplotPath=os.path.join(onepdfdir,figname)
        plotFig.savefig(oneDplotPath)
        if(savepdfs): plotFig.savefig(os.path.join(onepdfdir,par_name+'.pdf'))
        plt.close(plotFig)

        if rbins:
            print("r of injected value of %s (bins) = %f"%(par_name, rbins))

        ##Produce plot of raw samples
        myfig=plt.figure(figsize=(4,3.5),dpi=200)
        pos_samps=pos[par_name].samples
        if not ("chain" in pos.names):
            # If there is not a parameter named "chain" in the
            # posterior, then just produce a plot of the samples.
            plt.plot(pos_samps,'k.', markersize=5, alpha=0.5, linewidth=0.0, figure=myfig)
            maxLen=len(pos_samps)
        else:
            # If there is a parameter named "chain", then produce a
            # plot of the various chains in different colors, with
            # smaller dots.
            data,header=pos.samples()
            par_index=pos.names.index(par_name)
            chain_index=pos.names.index("chain")
            chains=unique(pos["chain"].samples)
            chainData=[data[ data[:,chain_index] == chain, par_index ] for chain in chains]
            chainDataRanges=[range(len(cd)) for cd in chainData]
            maxLen=max([len(cd) for cd in chainData])
            for rng, data in zip(chainDataRanges, chainData):
                plt.plot(rng, data, marker='.', markersize=1, alpha=0.5, linewidth=0.0,figure=myfig)
            plt.title("Gelman-Rubin R = %g"%(pos.gelman_rubin(par_name)))

            #dataPairs=[ [rng, data] for (rng,data) in zip(chainDataRanges, chainData)]
            #flattenedData=[ item for pair in dataPairs for item in pair ]
            #maxLen=max([len(data) for data in flattenedData])
            #plt.plot(array(flattenedData),marker=',',linewidth=0.0,figure=myfig)


        injpar=pos[par_name].injval

        if injpar is not None:
            # Allow injection to be 5% outside the posterior plot
            minrange=min(pos_samps)-0.05*(max(pos_samps)-min(pos_samps))
            maxrange=max(pos_samps)+0.05*(max(pos_samps)-min(pos_samps))
            if minrange<injpar and maxrange>injpar:
                plt.axhline(injpar, color='r', linestyle='-.',linewidth=4)
        myfig.savefig(os.path.join(sampsdir,figname.replace('.png','_samps.png')))
        if(savepdfs): myfig.savefig(os.path.join(sampsdir,figname.replace('.png','_samps.pdf')))
        plt.close(myfig)
        acfail=0
        if not (noacf):
            acffig=plt.figure(figsize=(4,3.5),dpi=200)
            if not ("chain" in pos.names):
                data=pos_samps[:,0]
                try:
                    (Neff, acl, acf) = effectiveSampleSize(data, Nskip)
                    lines=plt.plot(acf, 'k.', marker='.', markersize=1, alpha=0.5, linewidth=0.0, figure=acffig)
                    # Give ACL info if not already downsampled according to it
                    if nDownsample is None:
                        plt.title('Autocorrelation Function')
                    elif 'cycle' in pos.names:
                        last_color = lines[-1].get_color()
                        plt.axvline(acl/Nskip, linestyle='-.', color=last_color)
                        plt.title('ACL = %i   N = %i'%(acl,Neff))
                    acffig.savefig(os.path.join(sampsdir,figname.replace('.png','_acf.png')))
                    if(savepdfs): acffig.savefig(os.path.join(sampsdir,figname.replace('.png','_acf.pdf')))
                    plt.close(acffig)
                except:
                    # Ignore
                    acfail=1
                    pass
            else:
                try:
                    acls = []
                    Nsamps = 0.0;
                    for rng, data, Nskip, Ncycles in zip(chainDataRanges, chainData, chainNskips, chainNcycles):
                        (Neff, acl, acf) = effectiveSampleSize(data, Nskip)
                        acls.append(acl)
                        Nsamps += Neff
                        lines=plt.plot(acf,'k.', marker='.', markersize=1, alpha=0.5, linewidth=0.0, figure=acffig)
                        # Give ACL info if not already downsampled according to it
                        if nDownsample is not None:
                            last_color = lines[-1].get_color()
                            plt.axvline(acl/Nskip, linestyle='-.', color=last_color)
                    if nDownsample is None:
                        plt.title('Autocorrelation Function')
                    else:
                        plt.title('ACL = %i  N = %i'%(max(acls),Nsamps))
                    acffig.savefig(os.path.join(sampsdir,figname.replace('.png','_acf.png')))
                    if(savepdfs): acffig.savefig(os.path.join(sampsdir,figname.replace('.png','_acf.pdf')))
                    plt.close(acffig)
                except:
                    # Ignore
                    acfail=1
                    pass

        if not noacf:
            if not acfail:
                acfhtml='<td width="30%"><img width="100%" src="1Dsamps/'+figname.replace('.png', '_acf.png')+'"/></td>'
            else:
                acfhtml='<td>ACF generation failed!</td>'
            html_ompdf_write+='<tr><td width="30%"><img width="100%" src="1Dpdf/'+figname+'"/></td><td width="30%"><img width="100%" src="1Dsamps/'+figname.replace('.png','_samps.png')+'"/></td>'+acfhtml+'</tr>'
        else:
            html_ompdf_write+='<tr><td width="30%"><img width="100%" src="1Dpdf/'+figname+'"/></td><td width="30%"><img width="100%" src="1Dsamps/'+figname.replace('.png','_samps.png')+'"/></td></tr>'

    html_ompdf_write+='</table>'
    html_ompdf.write(html_ompdf_write)

    return out

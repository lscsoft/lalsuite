""" 
Compute the prior distribution in the mass and spin (Mf, af) of the final black hole 
corresponding to a uniform prior in component masses and spins. 

P. Ajith, 2015-09-20 

$Id:$
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np 
import lalinference.imrtgr.imrtgrutils as tgr
from scipy import interpolate as interp
import pickle, gzip
from optparse import OptionParser
from matplotlib import rc
import matplotlib
matplotlib.rc('text.latex', preamble = '\usepackage{txfonts}')

rc('text', usetex=True)
rc('font', family='serif')
rc('font', serif='times')
rc('mathtext', default='sf')
rc("lines", markeredgewidth=1)
rc("lines", linewidth=2)
rc('axes', labelsize=10)
rc("axes", linewidth=0.5)
rc('xtick', labelsize=8)
rc('ytick', labelsize=8)
rc('legend', fontsize=10)
rc('xtick.major', pad=6)
rc('ytick.major', pad=6)
rc('xtick.minor', size=5)
rc('ytick.minor', size=5)

def set_tick_sizes(ax, major, minor):
  for l in ax.get_xticklines() + ax.get_yticklines():
    l.set_markersize(major)
  for tick in ax.xaxis.get_minor_ticks() + ax.yaxis.get_minor_ticks():
    tick.tick1line.set_markersize(minor)
    tick.tick2line.set_markersize(minor)
  ax.xaxis.LABELPAD=10.
  ax.xaxis.OFFSETTEXTPAD=10.

# read inputs from command line 
parser = OptionParser()
parser.add_option("--fit-formula", dest="fit_formula",
                  help="fitting formula for the mass/spin of the final BH [options: 'nospin_Pan2011', 'nonprecspin_Healy2014']")
parser.add_option("--comp-mass-min", dest="comp_mass_min",
                  help="minimum value of the component mass prior [M_sun]")
parser.add_option("--comp-mass-max", dest="comp_mass_max",
                  help="maximum value of the component mass prior [M_sun]")
parser.add_option("--comp-spin-min", dest="comp_spin_min",
                  help="minimum value of the dimensionless spin prior [along the orb. ang. momentum]")
parser.add_option("--comp-spin-max", dest="comp_spin_max",
                  help="maximum value of the dimensionless spin prior [along the orb. ang. momentum]")
parser.add_option("--spin-angle-dist", dest="spin_angle_dist",
                  help="distribution of spin angles options: 'aligned', 'isotropic']")
parser.add_option("--Mf-lim", dest="Mf_lim",
                  help="the prior distribution will be calculated from -Mf_lim to Mf_lim [M_sun]")
parser.add_option("--af-lim", dest="af_lim",
                  help="the prior distribution will be calculated from -af_lim to af_lim")
parser.add_option("--num-threads", dest="num_threads",
                  help="number of threads to be used")
parser.add_option("--num-samples", dest="N_sampl",
                  help="number of prior samples to be used for computing the histogram")
parser.add_option("--num-bins", dest="N_bins",
                  help="number of bins to be used for computing the histogram")
(options, args) = parser.parse_args()

fit_formula = options.fit_formula
comp_mass_min = float(options.comp_mass_min)
comp_mass_max = float(options.comp_mass_max)
comp_spin_min = float(options.comp_spin_min)
comp_spin_max = float(options.comp_spin_max)
N_sampl = float(options.N_sampl)
N_bins = int(options.N_bins)
Mf_lim = float(options.Mf_lim)
af_lim = float(options.af_lim)
num_threads = int(options.num_threads)
if options.spin_angle_dist == None: 
  spin_angle_dist = 'aligned'
else: 
  spin_angle_dist = options.spin_angle_dist

# create the bins over which the histograms will be comptued 
Mf_bins = np.linspace(-Mf_lim, Mf_lim, N_bins)
af_bins = np.linspace(-af_lim, af_lim, N_bins)
P_Mfaf_pr = tgr.calc_Mfchif_prior(comp_mass_min, comp_mass_max, comp_spin_min, comp_spin_max, Mf_bins, af_bins, fit_formula, spin_angle_dist, N_sampl, num_threads)
print '... calculated the prior' 

# create an interpolation object and save it 
Mf_bins = (Mf_bins[:-1] + Mf_bins[1:])/2.
af_bins = (af_bins[:-1] + af_bins[1:])/2.
outfile = 'Prior_Mfaf_%s_comp_mass_min%2.1f_comp_mass_max%2.1f_comp_spin_min%2.1f_comp_spin_max%2.1f_%s'%(fit_formula, comp_mass_min, comp_mass_max, comp_spin_min, comp_spin_max, spin_angle_dist)
P_Mfaf_pr_interp_obj = interp.interp2d(Mf_bins, af_bins, P_Mfaf_pr, fill_value=0., bounds_error=False)
f = gzip.open(outfile+".pklz",'wb')
pickle.dump(P_Mfaf_pr_interp_obj, f)
print '... saved the interpolation object.' 

# read the interpolation object, reconstruct the data from the interpolation object. This is only used for estimating the error due to the interpolation 
f = gzip.open(outfile+".pklz",'rb')
P_Mfaf_pr_interp_obj = pickle.load(f)
P_Mfaf_pr_interp = P_Mfaf_pr_interp_obj(Mf_bins, af_bins)

# difference between the original and interpolated data 
interp_err = abs(P_Mfaf_pr-P_Mfaf_pr_interp)
print '... maximum difference between the original and interpolated data is %e' %np.max(interp_err)

plt.figure(figsize=(15,4))
plt.subplot(131)
plt.pcolormesh(Mf_bins, af_bins, P_Mfaf_pr, cmap='YlOrBr')
plt.xlabel('$M_f [M_{\odot}]$')
plt.ylabel('$a_f/M_f$')
plt.xlim(-Mf_lim, Mf_lim)
plt.ylim(-af_lim, af_lim)
plt.grid()
plt.colorbar()
plt.title('Original')

plt.subplot(132)
plt.pcolormesh(Mf_bins, af_bins, P_Mfaf_pr_interp, cmap='YlOrBr')
plt.xlabel('$M_f [M_{\odot}]$')
plt.ylabel('$a_f/M_f$')
plt.xlim(-Mf_lim, Mf_lim)
plt.ylim(-af_lim, af_lim)
plt.grid()
plt.colorbar()
plt.title('Interpolation')

plt.subplot(133)
plt.pcolormesh(Mf_bins, af_bins, interp_err, cmap='YlOrBr')
plt.xlabel('$M_f [M_{\odot}]$')
plt.ylabel('$a_f/M_f$')
plt.xlim(-Mf_lim, Mf_lim)
plt.ylim(-af_lim, af_lim)
plt.grid()
plt.clim(np.min(interp_err)-1e-16,np.max(interp_err)+1e-16)
plt.colorbar()
plt.title('Difference')
plt.savefig(outfile+'.png', dpi=300)

# save a zoomed version of the figure 
plt.subplot(131)
plt.xlim(50,150)
plt.ylim(0,1)
plt.subplot(132)
plt.xlim(50,150)
plt.ylim(0,1)
plt.subplot(133)
plt.xlim(50,150)
plt.ylim(0,1)
plt.savefig(outfile+'_zoom.png', dpi=300)

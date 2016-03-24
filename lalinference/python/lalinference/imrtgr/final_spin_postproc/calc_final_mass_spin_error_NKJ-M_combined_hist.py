"""
Estimate the errors in the fitting formulae used to compute the final mass/spin by comparing the final mass/spin estimated using the formula with NR results from the SXS catalog. 

P. Ajith, 2015-11-15, extended by NKJ-M, 12.2015-03.2016
"""

import matplotlib.pyplot as plt 
import plotsettings 
import numpy as np 
import nrutils_new as nr
from mpl_toolkits.mplot3d import Axes3D
import pneqns
#import pneqns_orig as pneqns
import os 

nr_data = 'SXS_catalog_params.txt'
#fit_tag = 'HLZ'
#fit_tag = 'HLZ_prec_extension_Mf'
fit_tag = 'HLZ_prec_extension_Minit'
#fit_tag = 'Husa_etal_prec_extension_Minit'
#fit_tag = 'Barausse_and_Rezzolla'

v_final = 6.**-0.5 # Final v for the evolution

evolve_spins = 1

v0L = 1

if evolve_spins:
	if v0L:
		evol_tag = 'evol_v0L'
	else:
		evol_tag = 'evol'
else:
	evol_tag = 'no_evol'

# Define final mass and spin functions

def final_mass(m1,m2,chi1z,chi2z,fit_tag):
	if fit_tag == 'HLZ':
		return nr.bbh_final_mass_non_precessing_Healyetal(m1,m2,chi1z,chi2z)
	if fit_tag == 'HLZ_prec_extension_Mf':
                return nr.bbh_final_mass_non_precessing_Healyetal(m1,m2,chi1z,chi2z)
	if fit_tag == 'HLZ_prec_extension_Minit':
		return nr.bbh_final_mass_non_precessing_Healyetal(m1,m2,chi1z,chi2z)
	if fit_tag == 'Husa_etal_prec_extension_Minit':
		return nr.bbh_final_mass_non_precessing_Husaetal(m1,m2,chi1z,chi2z)
	if fit_tag == 'Barausse_and_Rezzolla':
		return nr.bbh_final_mass_non_precessing_Healyetal(m1,m2,chi1z,chi2z)

def final_spin(m1,m2,chi1,chi2,tilt1,tilt2,phi12,fit_tag):
	if fit_tag == 'HLZ':
		return nr.bbh_final_spin_non_precessing_Healyetal(m1,m2,chi1*np.cos(tilt1),chi2*np.cos(tilt2))
	if fit_tag == 'HLZ_prec_extension_Mf':
		return nr.bbh_final_spin_precessing_Healyetal_extension_Mf(m1,m2,chi1,chi2,tilt1,tilt2,phi12)
	if fit_tag == 'HLZ_prec_extension_Minit':
		return nr.bbh_final_spin_precessing_Healyetal_extension_Minit(m1,m2,chi1,chi2,tilt1,tilt2,phi12)
	if fit_tag == 'Husa_etal_prec_extension_Minit':
		return nr.bbh_final_spin_precessing_Husaetal_extension_Minit(m1,m2,chi1,chi2,tilt1,tilt2,phi12)
	if fit_tag == 'Barausse_and_Rezzolla':
		return nr.bbh_final_spin_precessing_Barausse_and_Rezzolla(m1,m2,chi1,chi2,tilt1,tilt2,phi12)

# parameters for selecting all quasi-circular simulations and selecting out the precessing ones, as well
out_dir = 'SXS_allspins'
if fit_tag == 'HLZ':
	tag = 'allspins_%s_L0ADM_corr_test2_final_%s' %(evol_tag,fit_tag)
else:
	tag = 'allspins_incl_in-plane_spin_%s_L0ADM_corr_test2_final_%s' %(evol_tag,fit_tag)
MIN_IN_PLANE_SPINS = 1e-3
MAX_IN_PLANE_SPINS = 1.
MAX_Q = 18.
MAX_ECCENRICITY = 1e-3 

# make output directory and copy the run script there 
os.system('mkdir -p %s' %out_dir)
os.system('cp %s %s' %(__file__, out_dir))

# read the data 
SXS_run_tag = np.loadtxt(nr_data, dtype=str, skiprows=1, usecols=(0,), unpack=True)
q, chi1, chi2, chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z, eccentricity, omega0, num_orbs, af, af_x, af_y, af_z, Mf, m1, m2, Jx,	Jy,	Jz,	J = np.loadtxt(nr_data, usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23), unpack=True)

# select a subset of simulations 
idx = np.logical_and(np.logical_and(np.logical_and(abs(chi1_x) < MAX_IN_PLANE_SPINS, abs(chi2_x) < MAX_IN_PLANE_SPINS), q < MAX_Q), eccentricity < MAX_ECCENRICITY)

q, chi1, chi2, chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z, eccentricity, omega0, num_orbs, af, af_x, af_y, af_z, Mf, m1, m2, Jx,	Jy,	Jz,	J, SXS_run_tag = q[idx], chi1[idx], chi2[idx], chi1_x[idx], chi1_y[idx], chi1_z[idx], chi2_x[idx], chi2_y[idx], chi2_z[idx], eccentricity[idx], omega0[idx], num_orbs[idx], af[idx], af_x[idx], af_y[idx], af_z[idx], Mf[idx], m1[idx], m2[idx], Jx[idx],	Jy[idx],	Jz[idx],	J[idx], SXS_run_tag[idx]

idx_prec = np.logical_and(np.logical_or(abs(chi1_x) > MIN_IN_PLANE_SPINS, abs(chi2_x) > MIN_IN_PLANE_SPINS), eccentricity < MAX_ECCENRICITY)

# compute the orbital ang momentum 
Lx = Jx-(chi1_x*m1**2 + chi2_x*m2**2)
Ly = Jy-(chi1_y*m1**2 + chi2_y*m2**2)
Lz = Jz-(chi1_z*m1**2 + chi2_z*m2**2)
L = (Lx**2+Ly**2+Lz**2)**0.5

# unit vector along the initial ADM ang momentum 
Lx /= L
Ly /= L
Lz /= L

# rescale all results to unit M 
M = m1+m2 
eta = m1*m2/(M*M)
m1 /= M 
m2 /= M 
af /= Mf**2. 
Mf /= M 

if evolve_spins:
	# compute chi_i . L at Schwarzschild ISCO 
	chi1dL = np.zeros_like(q)
	chi2dL = np.zeros_like(q)
	phi12 = np.zeros_like(q)
	Sperpmag = np.zeros_like(q) 
	v_f = np.zeros_like(q)
	for i in range(len(q)):
		if v0L:
			aL = M[i]*eta[i]/L[i]
			v0 = aL + (1.5 + eta[i]/6)*aL*aL*aL # Consistent 1PN expression for v0 in terms of L
		else:
			v0 = omega0[i]**(1./3)
		print 'i = %d v0 = %f chi1z = %f chi2z = %f. #'  %(i,  v0, chi1_z[i], chi2_z[i]),  

		# in the case of non-precessing spins chi_i . L = chi_z 
		if abs(chi1_x[i]) < 1e-3 and abs(chi2_x[i]) < 1e-3:
			chi1dL[i], chi2dL[i], Sperpmag[i], v_f[i] = chi1_z[i], chi2_z[i], 0., v_final
		# in the case of precessing simulations estimate chi_i . L at v_isco by evolving PN precession eqns 
		else: 
			chi1x_v, chi1y_v, chi1z_v, chi2x_v, chi2y_v, chi2z_v, Lnx_v, Lny_v, Lnz_v, v_f[i] = pneqns.find_S_and_L_at_freq_dt(v0, m1[i]*10, m2[i]*10, chi1_x[i], chi1_y[i], chi1_z[i], chi2_x[i], chi2_y[i], chi2_z[i], Lx[i], Ly[i], Lz[i], v_final,25.)
			m1sq, m2sq = m1[i]*m1[i], m2[i]*m2[i]
			chi1dL[i] = chi1x_v*Lnx_v + chi1y_v*Lny_v + chi1z_v*Lnz_v
			chi2dL[i] = chi2x_v*Lnx_v + chi2y_v*Lny_v + chi2z_v*Lnz_v
			chi1inplanex = chi1x_v - chi1dL[i]*Lnx_v
			chi1inplaney = chi1y_v - chi1dL[i]*Lny_v
			chi1inplanez = chi1z_v - chi1dL[i]*Lnz_v
			chi2inplanex = chi2x_v - chi2dL[i]*Lnx_v
			chi2inplaney = chi2y_v - chi2dL[i]*Lny_v
			chi2inplanez = chi2z_v - chi2dL[i]*Lnz_v
			Sperpx = m1sq*chi1inplanex + m2sq*chi2inplanex
			Sperpy = m1sq*chi1inplaney + m2sq*chi2inplaney
			Sperpz = m1sq*chi1inplanez + m2sq*chi2inplanez
			Sperpmag[i] = (Sperpx*Sperpx + Sperpy*Sperpy + Sperpz*Sperpz)**0.5
			phi12[i] = np.arccos((chi1inplanex*chi2inplanex + chi1inplaney*chi2inplaney + chi1inplanez*chi2inplanez)/((chi1inplanex*chi1inplanex + chi1inplaney*chi1inplaney + chi1inplanez*chi1inplanez)*(chi2inplanex*chi2inplanex + chi2inplaney*chi2inplaney + chi2inplanez*chi2inplanez))**0.5)
		print 'v_f = %f: chi1dL = %f chi2dL = %f Sperpmag = %f' %(v_f[i], chi1dL[i], chi2dL[i], Sperpmag[i])

if not evolve_spins:

	chi1dL = chi1_x*Lx + chi1_y*Ly + chi1_z*Lz 
	chi2dL = chi2_x*Lx + chi2_y*Ly + chi2_z*Lz

	chi1inplanex = chi1_x - chi1dL*Lx
	chi1inplaney = chi1_y - chi1dL*Ly
	chi1inplanez = chi1_z - chi1dL*Lz

	chi2inplanex = chi2_x - chi2dL*Lx
	chi2inplaney = chi2_y - chi2dL*Ly
	chi2inplanez = chi2_z - chi2dL*Lz

	Sperpmag = ((m1*m1*chi1inplanex + m2*m2*chi2inplanex)**2. + (m1*m1*chi1inplaney + m2*m2*chi2inplaney)**2. + (m1*m1*chi1inplanez + m2*m2*chi2inplanez)**2.)**0.5

	phi12 = np.arccos((chi1inplanex*chi2inplanex + chi1inplaney*chi2inplaney + chi1inplanez*chi2inplanez)/((chi1inplanex*chi1inplanex + chi1inplaney*chi1inplaney + chi1inplanez*chi1inplanez)*(chi2inplanex*chi2inplanex + chi2inplaney*chi2inplaney + chi2inplanez*chi2inplanez))**0.5)

# tilts

tilt1 = np.arccos(chi1dL/chi1)
tilt2 = np.arccos(chi2dL/chi2)

# compute the final mass and spin 

Mf_fit = final_mass(m1, m2, chi1dL, chi2dL,fit_tag)
af_fit = final_spin(m1, m2, chi1, chi2, tilt1, tilt2,phi12,fit_tag)

# compute the errors in the fits 
Mf_err = Mf_fit-Mf 
Mf_err_norm = Mf_fit/Mf - 1.
af_err = abs(af_fit)-af 		# the af in the SXS catalogue is the magnitude of the final spin. Hence we take the abs of the final spin given by the fit formula 
af_err_norm = abs(af_fit)/af - 1.
print '# RMS normalized error in Mf = %3.2e. RMS normalized error in af = %3.2e' %(np.std(Mf_err_norm), np.std(af_err_norm))

# print the simulations with large fractional errors 
err_idx_vec = np.where(abs(af_err_norm) > np.std(af_err_norm))[0]
for err_idx in err_idx_vec: 
        print 'tag = %s q = %2.1f chi1 = [%2.1f, %2.1f, %2.1f] chi2 = [%2.1f, %2.1f, %2.1f] af = %f af_fit = %f af_err_norm = %f' %(SXS_run_tag[err_idx], q[err_idx], chi1_x[err_idx], chi1_y[err_idx], chi1_z[err_idx], chi2_x[err_idx], chi2_y[err_idx], chi2_z[err_idx], af[err_idx], af_fit[err_idx], af_err_norm[err_idx])

# print the simulations with large absolute errors 
err_idx_vec = np.where(abs(af_err) > np.std(af_err))[0]
for err_idx in err_idx_vec: 
	print 'tag = %s q = %2.1f chi1 = [%2.1f, %2.1f, %2.1f] chi2 = [%2.1f, %2.1f, %2.1f] af = %f af_fit = %f af_err = %f' %(SXS_run_tag[err_idx], q[err_idx], chi1_x[err_idx], chi1_y[err_idx], chi1_z[err_idx], chi2_x[err_idx], chi2_y[err_idx], chi2_z[err_idx], af[err_idx], af_fit[err_idx], af_err[err_idx])

# Compute histograms (first multiplying the errors by 100)

Mf_err *= 100
Mf_err_norm *= 100

af_err *= 100
af_err_norm *= 100

Mf_err_hist, Mf_err_hist_bins = np.histogram(Mf_err, bins=20)
Mf_err_hist_prec, Mf_err_hist_prec_bins = np.histogram(Mf_err[idx_prec], bins=Mf_err_hist_bins)

af_err_hist, af_err_hist_bins = np.histogram(af_err, bins=20)
af_err_hist_prec, af_err_hist_prec_bins = np.histogram(af_err[idx_prec], bins=af_err_hist_bins)

Mf_err_norm_hist, Mf_err_norm_hist_bins = np.histogram(Mf_err_norm, bins=20)
Mf_err_norm_hist_prec, Mf_err_norm_hist_prec_bins = np.histogram(Mf_err_norm[idx_prec], bins=Mf_err_norm_hist_bins)

af_err_norm_hist, af_err_norm_hist_bins = np.histogram(af_err_norm, bins=20)
af_err_norm_hist_prec, af_err_norm_hist_prec_bins = np.histogram(af_err_norm[idx_prec], bins=af_err_norm_hist_bins)

# Set up y axis limits

Mfymin, Mfymax = (1 - 0.1*np.sign(min(Mf_err)))*min(Mf_err), (1 + 0.1*np.sign(max(Mf_err)))*max(Mf_err)
Mfymin_norm, Mfymax_norm = (1 - 0.1*np.sign(min(Mf_err_norm)))*min(Mf_err_norm), (1 + 0.1*np.sign(max(Mf_err_norm)))*max(Mf_err_norm)

afymin, afymax = (1 - 0.3*np.sign(min(af_err)))*min(af_err), (1 + 0.1*np.sign(max(af_err)))*max(af_err)
afymin_norm, afymax_norm = (1 - 0.3*np.sign(min(af_err_norm)))*min(af_err_norm), (1 + 0.1*np.sign(max(af_err_norm)))*max(af_err_norm)

# Plotting histograms

plt.figure(figsize=(8,3.5))
plt.subplot(121)
plt.bar(Mf_err_hist_bins[:-1], Mf_err_hist, width=Mf_err_hist_bins[1] - Mf_err_hist_bins[0])
plt.bar(Mf_err_hist_bins[:-1], Mf_err_hist_prec, width=Mf_err_hist_bins[1] - Mf_err_hist_bins[0],color='purple')
plt.grid()
plt.xlabel('$\Delta M_f \\times 100$')
plt.ylabel('$N$')
plt.title('RMS error = %2.1e' %np.std(Mf_err/100), fontsize=8)
plt.subplot(122)
plt.bar(af_err_hist_bins[:-1], af_err_hist, width=af_err_hist_bins[1] - af_err_hist_bins[0])
plt.bar(af_err_hist_bins[:-1], af_err_hist_prec, width=af_err_hist_bins[1] - af_err_hist_bins[0],color='purple')
plt.grid()
plt.xlabel('$\Delta \chi_f \\times 100 $')
plt.ylabel('$N$')
plt.title('RMS error = %2.1e' %np.std(af_err/100), fontsize=8)
plt.tight_layout()
plt.savefig('%s/Mf_af_err_hist_combined_%s.png' %(out_dir, tag), dpi=200)
plt.close()

plt.figure(figsize=(8,3.5))
plt.subplot(121)
plt.bar(Mf_err_norm_hist_bins[:-1], Mf_err_norm_hist, width=Mf_err_norm_hist_bins[1] - Mf_err_norm_hist_bins[0])
plt.bar(Mf_err_norm_hist_bins[:-1], Mf_err_norm_hist_prec, width=Mf_err_norm_hist_bins[1] - Mf_err_norm_hist_bins[0],color='purple')
plt.grid()
plt.xlabel('$\Delta M_f/M_f \\times 100$')
plt.ylabel('$N$')
plt.title('RMS error = %2.1e' %np.std(Mf_err_norm/100), fontsize=8)
plt.subplot(122)
plt.bar(af_err_norm_hist_bins[:-1], af_err_norm_hist, width=af_err_norm_hist_bins[1] - af_err_norm_hist_bins[0])
plt.bar(af_err_norm_hist_bins[:-1], af_err_norm_hist_prec, width=af_err_norm_hist_bins[1] - af_err_norm_hist_bins[0],color='purple')
plt.grid()
plt.xlabel('$\Delta \chi_f/\chi_f \\times 100 $')
plt.ylabel('$N$')
plt.title('RMS error = %2.1e' %np.std(af_err_norm/100), fontsize=8)
plt.tight_layout()
plt.savefig('%s/Mf_af_err_norm_hist_combined_%s.png' %(out_dir, tag), dpi=200)
plt.close()

# Plotting absolute errors vs. various parameters

plt.figure(figsize=(8,3.5))
plt.subplot(121)
plt.scatter(q, Mf_err)
plt.scatter(q[idx_prec], Mf_err[idx_prec], color='purple')
plt.ylim(Mfymin,Mfymax)
plt.grid()
plt.xlabel('$q$')
plt.ylabel('$\Delta M_f \\times 100$')
plt.subplot(122)
plt.scatter(q, af_err)
plt.scatter(q[idx_prec], af_err[idx_prec], color='purple')
plt.ylim(afymin,afymax)
plt.ylim(afymin,afymax)
plt.grid()
plt.xlabel('$q$')
plt.ylabel('$\Delta \chi_f \\times 100$')
plt.tight_layout()
plt.savefig('%s/Mf_af_err_vs_q_combined_%s.png' %(out_dir, tag), dpi=200)
plt.close()

plt.figure(figsize=(8,3.5))
plt.subplot(121)
plt.scatter(chi1, Mf_err)
plt.scatter(chi1[idx_prec], Mf_err[idx_prec], color='purple')
plt.ylim(Mfymin,Mfymax)
plt.grid()
plt.xlabel('$\chi_1$')
plt.ylabel('$\Delta M_f \\times 100$')
plt.subplot(122)
plt.scatter(chi1, af_err)
plt.scatter(chi1[idx_prec], af_err[idx_prec], color='purple')
plt.ylim(afymin,afymax)
plt.grid()
plt.xlabel('$\chi_1$')
plt.ylabel('$\Delta \chi_f \\times 100$')
plt.tight_layout()
plt.savefig('%s/Mf_af_err_vs_chi1_combined_%s.png' %(out_dir, tag), dpi=200)
plt.close()

plt.figure(figsize=(8,3.5))
plt.subplot(121)
plt.scatter(chi2, Mf_err)
plt.scatter(chi2[idx_prec], Mf_err[idx_prec], color='purple')
plt.ylim(Mfymin,Mfymax)
plt.grid()
plt.xlabel('$\chi_2$')
plt.ylabel('$\Delta M_f \\times 100$')
plt.subplot(122)
plt.scatter(chi2, af_err)
plt.scatter(chi2[idx_prec], af_err[idx_prec], color='purple')
plt.ylim(afymin,afymax)
plt.grid()
plt.xlabel('$\chi_2$')
plt.ylabel('$\Delta \chi_f \\times 100$')
plt.tight_layout()
plt.savefig('%s/Mf_af_err_vs_chi2_combined_%s.png' %(out_dir, tag), dpi=200)
plt.close()

plt.figure(figsize=(8,3.5))
plt.subplot(121)
plt.scatter(chi1*np.sin(tilt1), Mf_err)
plt.scatter(chi1[idx_prec]*np.sin(tilt1[idx_prec]), Mf_err[idx_prec], color='purple')
plt.ylim(Mfymin,Mfymax)
plt.grid()
plt.xlabel('$\chi_{1\perp}$')
plt.ylabel('$\Delta M_f \\times 100$')
plt.subplot(122)
plt.scatter(chi1*np.sin(tilt1), af_err)
plt.scatter(chi1[idx_prec]*np.sin(tilt1[idx_prec]), af_err[idx_prec], color='purple')
plt.ylim(afymin,afymax)
plt.grid()
plt.xlabel('$\chi_{1\perp}$')
plt.ylabel('$\Delta \chi_f \\times 100$')
plt.tight_layout()
plt.savefig('%s/Mf_af_err_vs_chi1perp_combined_%s.png' %(out_dir, tag), dpi=200)
plt.close()

plt.figure(figsize=(8,3.5))
plt.subplot(121)
plt.scatter(chi2*np.sin(tilt2), Mf_err)
plt.scatter(chi2[idx_prec]*np.sin(tilt2[idx_prec]), Mf_err[idx_prec], color='purple')
plt.ylim(Mfymin,Mfymax)
plt.grid()
plt.xlabel('$\chi_{2\perp}$')
plt.ylabel('$\Delta M_f \\times 100$')
plt.subplot(122)
plt.scatter(chi2*np.sin(tilt2), af_err)
plt.scatter(chi2[idx_prec]*np.sin(tilt2[idx_prec]), af_err[idx_prec], color='purple')
plt.ylim(afymin,afymax)
plt.grid()
plt.xlabel('$\chi_{2\perp}$')
plt.ylabel('$\Delta \chi_f \\times 100$')
plt.tight_layout()
plt.savefig('%s/Mf_af_err_vs_chi2perp_combined_%s.png' %(out_dir, tag), dpi=200)
plt.close()

plt.figure(figsize=(8,3.5))
plt.subplot(121)
plt.scatter(Sperpmag, Mf_err)
plt.scatter(Sperpmag[idx_prec], Mf_err[idx_prec], color='purple')
plt.ylim(Mfymin,Mfymax)
plt.grid()
plt.xlabel('$\|\mathbf{S}_\perp\|/M^2$')
plt.ylabel('$\Delta M_f \\times 100$')
plt.tick_params(axis='x', labelsize=10)
plt.subplot(122)
plt.scatter(Sperpmag, af_err)
plt.scatter(Sperpmag[idx_prec], af_err[idx_prec], color='purple')
plt.ylim(afymin,afymax)
plt.grid()
plt.xlabel('$\|\mathbf{S}_\perp\|/M^2$')
plt.ylabel('$\Delta \chi_f \\times 100$')
plt.tick_params(axis='x', labelsize=10)
plt.tight_layout()
plt.savefig('%s/Mf_af_err_vs_Sperpmag_combined_%s.png' %(out_dir, tag), dpi=200)
plt.close()

# Plotting fractional errors vs. various parameters

plt.figure(figsize=(8,3.5))
plt.subplot(121)
plt.scatter(q, Mf_err_norm)
plt.scatter(q[idx_prec], Mf_err_norm[idx_prec], color='purple')
plt.ylim(Mfymin_norm,Mfymax_norm)
plt.grid()
plt.xlabel('$q$')
plt.ylabel('$\Delta M_f/M_f \\times 100$')
plt.subplot(122)
plt.scatter(q, af_err_norm)
plt.scatter(q[idx_prec], af_err_norm[idx_prec], color='purple')
plt.ylim(afymin_norm,afymax_norm)
plt.grid()
plt.xlabel('$q$')
plt.ylabel('$\Delta \chi_f/\chi_f \\times 100$')
plt.tight_layout()
plt.savefig('%s/Mf_af_err_norm_vs_q_combined_%s.png' %(out_dir, tag), dpi=200)
plt.close()

plt.figure(figsize=(8,3.5))
plt.subplot(121)
plt.scatter(chi1, Mf_err_norm)
plt.scatter(chi1[idx_prec], Mf_err_norm[idx_prec], color='purple')
plt.ylim(Mfymin_norm,Mfymax_norm)
plt.grid()
plt.xlabel('$\chi_1$')
plt.ylabel('$\Delta M_f/M_f \\times 100$')
plt.subplot(122)
plt.scatter(chi1, af_err_norm)
plt.scatter(chi1[idx_prec], af_err_norm[idx_prec], color='purple')
plt.ylim(afymin_norm,afymax_norm)
plt.grid()
plt.xlabel('$\chi_1$')
plt.ylabel('$\Delta \chi_f/\chi_f \\times 100$')
plt.tight_layout()
plt.savefig('%s/Mf_af_err_norm_vs_chi1_combined_%s.png' %(out_dir, tag), dpi=200)
plt.close()

plt.figure(figsize=(8,3.5))
plt.subplot(121)
plt.scatter(chi2, Mf_err_norm)
plt.scatter(chi2[idx_prec], Mf_err_norm[idx_prec], color='purple')
plt.ylim(Mfymin_norm,Mfymax_norm)
plt.grid()
plt.xlabel('$\chi_2$')
plt.ylabel('$\Delta M_f/M_f \\times 100$')
plt.subplot(122)
plt.scatter(chi2, af_err_norm)
plt.scatter(chi2[idx_prec], af_err_norm[idx_prec], color='purple')
plt.ylim(afymin_norm,afymax_norm)
plt.grid()
plt.xlabel('$\chi_2$')
plt.ylabel('$\Delta \chi_f/\chi_f \\times 100$')
plt.tight_layout()
plt.savefig('%s/Mf_af_err_norm_vs_chi2_combined_%s.png' %(out_dir, tag), dpi=200)
plt.close()

plt.figure(figsize=(8,3.5))
plt.subplot(121)
plt.scatter(chi1*np.sin(tilt1), Mf_err_norm)
plt.scatter(chi1[idx_prec]*np.sin(tilt1[idx_prec]), Mf_err_norm[idx_prec], color='purple')
plt.ylim(Mfymin_norm,Mfymax_norm)
plt.grid()
plt.xlabel('$\chi_{1\perp}$')
plt.ylabel('$\Delta M_f/M_f \\times 100$')
plt.subplot(122)
plt.scatter(chi1*np.sin(tilt1), af_err_norm)
plt.scatter(chi1[idx_prec]*np.sin(tilt1[idx_prec]), af_err_norm[idx_prec], color='purple')
plt.ylim(afymin_norm,afymax_norm)
plt.grid()
plt.xlabel('$\chi_{1\perp}$')
plt.ylabel('$\Delta \chi_f/\chi_f \\times 100$')
plt.tight_layout()
plt.savefig('%s/Mf_af_err_norm_vs_chi1perp_combined_%s.png' %(out_dir, tag), dpi=200)
plt.close()

plt.figure(figsize=(8,3.5))
plt.subplot(121)
plt.scatter(chi2*np.sin(tilt2), Mf_err_norm)
plt.scatter(chi2[idx_prec]*np.sin(tilt2[idx_prec]), Mf_err_norm[idx_prec], color='purple')
plt.ylim(Mfymin_norm,Mfymax_norm)
plt.grid()
plt.xlabel('$\chi_{2\perp}$')
plt.ylabel('$\Delta M_f/M_f \\times 100$')
plt.subplot(122)
plt.scatter(chi2*np.sin(tilt2), af_err_norm)
plt.scatter(chi2[idx_prec]*np.sin(tilt2[idx_prec]), af_err_norm[idx_prec], color='purple')
plt.ylim(afymin_norm,afymax_norm)
plt.grid()
plt.xlabel('$\chi_{2\perp}$')
plt.ylabel('$\Delta \chi_f/\chi_f \\times 100$')
plt.tight_layout()
plt.savefig('%s/Mf_af_err_norm_vs_chi2perp_combined_%s.png' %(out_dir, tag), dpi=200)
plt.close()

plt.figure(figsize=(8,3.5))
plt.subplot(121)
plt.scatter(Sperpmag, Mf_err_norm)
plt.scatter(Sperpmag[idx_prec], Mf_err_norm[idx_prec], color='purple')
plt.ylim(Mfymin_norm,Mfymax_norm)
plt.grid()
plt.xlabel('$\|\mathbf{S}_\perp\|/M^2$')
plt.ylabel('$\Delta M_f/M_f \\times 100$')
plt.tick_params(axis='x', labelsize=10)
plt.subplot(122)
plt.scatter(Sperpmag, af_err_norm)
plt.scatter(Sperpmag[idx_prec], af_err_norm[idx_prec], color='purple')
plt.ylim(afymin_norm,afymax_norm)
plt.grid()
plt.xlabel('$\|\mathbf{S}_\perp\|/M^2$')
plt.ylabel('$\Delta \chi_f/\chi_f \\times 100$')
plt.tick_params(axis='x', labelsize=10)
plt.tight_layout()
plt.savefig('%s/Mf_af_err_norm_vs_Sperpmag_combined_%s.png' %(out_dir, tag), dpi=200)
plt.close()

# Plotting versus final mass and spin

"""
plt.figure(figsize=(8,3.5))
plt.subplot(121)
plt.scatter(Mf, Mf_err)
plt.scatter(Mf[idx_prec], Mf_err[idx_prec], color='purple')
plt.ylim(0.9*min(Mf_err),1.1*max(Mf_err))
plt.grid()
plt.xlabel('$M_f/M$')
plt.ylabel('$\Delta M_f \\times 100$')
plt.subplot(122)
plt.scatter(Mf, af_err)
plt.scatter(Mf[idx_prec], af_err[idx_prec], color='purple')
plt.ylim(0.9*min(af_err),1.1*max(af_err))
plt.grid()
plt.xlabel('$M_f$')
plt.ylabel('$\Delta \chi_f \\times 100$')
plt.tight_layout()
plt.savefig('%s/Mf_af_err_vs_Mf_combined_%s.png' %(out_dir, tag), dpi=200)
plt.close()

plt.figure(figsize=(8,3.5))
plt.subplot(121)
plt.scatter(af, Mf_err)
plt.scatter(af[idx_prec], Mf_err[idx_prec], color='purple')
plt.ylim(0.9*min(Mf_err),1.1*max(Mf_err))
plt.grid()
plt.xlabel('$\chi_f$')
plt.ylabel('$\Delta M_f \\times 100$')
plt.subplot(122)
plt.scatter(af, af_err)
plt.scatter(af[idx_prec], af_err[idx_prec], color='purple')
plt.ylim(0.9*min(af_err),1.1*max(af_err))
plt.grid()
plt.xlabel('$\chi_f$')
plt.ylabel('$\Delta \chi_f \\times 100$')
plt.tight_layout()
plt.savefig('%s/Mf_af_err_vs_af_combined_%s.png' %(out_dir, tag), dpi=200)
plt.close()
"""

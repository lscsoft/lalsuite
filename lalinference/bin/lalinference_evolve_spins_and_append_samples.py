'''
Code to read in LALInference HDF5 posterior samples, evolve the spins up to a specified value (default: Schwarzschild ISCO), and write the resulting evolved spin samples to the HDF5 file.

Usage:
python lalinference_evolve_spins_and_append_samples.py --sample_file <directory path to the posterior samples> [--vfinal <final orbital velocity>]
'''
# NKJ-M and Anuradha Gupta, 05.2019

import numpy as np
import argparse
import lalsimulation as lalsim
import h5py

from numpy.linalg import norm
from lal import MSUN_SI, MTSUN_SI

import matplotlib
matplotlib.use('Agg')

from lalinference.io.hdf5 import read_samples, write_samples, extract_metadata
from astropy.table import Column, Table, hstack

from lalinference.bayespputils import mc2ms

# Define a convenience function
def tilts_and_phi12_from_Cartesian_spins_and_L(chi1_v, chi2_v, Ln_v):

        # norms and normalizing
        chi1_v_norm = norm(chi1_v)
        chi2_v_norm = norm(chi2_v)

        Ln_v /= norm(Ln_v)

        # dot products
        chi1dL_v = np.dot(chi1_v, Ln_v)
        chi2dL_v = np.dot(chi2_v, Ln_v)

        # in-plane spins

        chi1inplane = chi1_v - chi1dL_v*Ln_v
        chi2inplane = chi2_v - chi2dL_v*Ln_v

        # computing cosine of tilts and phi12
        cos_tilt1 = chi1dL_v/chi1_v_norm
        cos_tilt2 = chi2dL_v/chi2_v_norm
        cos_phi12 = np.dot(chi1inplane,chi2inplane)/(norm(chi1inplane)*norm(chi2inplane))

        # set quadrant of phi12
        phi12_evol_i = np.arccos(cos_phi12)
        if np.sign(np.dot(Ln_v,np.cross(chi1_v, chi2_v))) < 0:
            phi12_evol_i = 2.*np.pi - phi12_evol_i

        return np.arccos(cos_tilt1), np.arccos(cos_tilt2), phi12_evol_i

# Settings
dt = 0.1 # steps in time for the integration, in terms of the total mass of the binary

# Choose which approximant to use for spin evolution
approx = lalsim.GetApproximantFromString("SpinTaylorT5")

# Set up the parsing
parser = argparse.ArgumentParser(description = 'Evolve the LALInference posterior samples of spins and append to HDF5 file of samples')
parser.add_argument("--sample_file", help = "path to the HDF5 posterior samples file", required=True)
parser.add_argument("--vfinal", help = "final orbital velocity for the evolution (default is the Schwarzschild ISCO velocity 6**-0.5 ~= 0.408)", type=str, default="ISCO")
args = parser.parse_args()

hdf_pos_file = args.sample_file

if args.vfinal=='ISCO':
    v_final = 6.**-0.5
    label = '_isco'
else:
    try:
        v_final = float(args.vfinal)
    except ValueError:
        raise ValueError("vfinal has to either be ISCO or a fraction of the speed of light")
    label = '_evol_vfinal'

# Functions for table formatting

def _identity(x):
    return x

_colname_map = (('tilt_spin1' + label, 'tilt1' + label, _identity),
                 ('tilt_spin2' + label, 'tilt2' + label, _identity))

def _remap_colnames(table):
    for old_name, new_name, func in _colname_map:
        if old_name in table.colnames:
            table[new_name] = func(table.columns.pop(old_name))

# read the posterior samples and metadata
data = read_samples(hdf_pos_file)
metadata = {}
run_identifier = extract_metadata(hdf_pos_file, metadata)

# extract mass and spin samples
if ('mc' in data.dtype.names) and ('q' in data.dtype.names):
  q = np.atleast_1d(data['q'])
  eta = q/(1. + q)/(1. + q)
  m1, m2 = mc2ms(np.atleast_1d(data['mc']), eta) # We use the redshifted masses here, since the starting frequency is in the detector frame
else:
  raise ValueError("Chirp mass and mass ratio are not found in %s. There is likely something wrong with this posterior file."%hdf_pos_file)
if ('a1' in data.dtype.names) and ('a2' in data.dtype.names):
  chi1, chi2 = np.atleast_1d(data['a1']), np.atleast_1d(data['a2'])
else:
  raise ValueError("No spin magnitude samples were found in %s. There are no spins to evolve."%hdf_pos_file)
if ('tilt1' in data.dtype.names) and ('tilt2' in data.dtype.names):
  tilt1, tilt2 = np.atleast_1d(data['tilt1']), np.atleast_1d(data['tilt2'])
else:
  raise ValueError("No tilt angle samples were found in %s; cannot evolve spins."%hdf_pos_file)
if 'phi12' in data.dtype.names:
  phi12 = np.atleast_1d(data['phi12'])
else:
  raise ValueError("No phi12 samples were found in %s; cannot evolve spins."%hdf_pos_file)

# Extract start frequency appropriate for the waveform used in the analysis

if 'LAL_APPROXIMANT' not in data.dtype.names:
    raise ValueError('LAL_APPROXIMANT is missing--unable to determine which waveform was used to obtain the samples.')

spinfreq_enum = np.array([lalsim.SimInspiralGetSpinFreqFromApproximant(int(lal_approx)) for lal_approx in data['LAL_APPROXIMANT']])

if len(np.where(spinfreq_enum == lalsim.SIM_INSPIRAL_SPINS_CASEBYCASE)[0]) > 0:
        raise ValueError('Unable to evolve spins since approximant does not have a set frequency at which the spins are defined.')

f_start = np.where(np.array(spinfreq_enum == lalsim.SIM_INSPIRAL_SPINS_FLOW), data['flow'], data['f_ref'])

# Evolve spins
tilt1_evol = np.zeros_like(m1)
tilt2_evol = np.zeros_like(m1)
phi12_evol = np.zeros_like(m1)

for i, _ in enumerate(m1):
    # Only evolve spins if at least one spin magnitude is above 1e-3
    if np.logical_or(chi1[i] > 1e-3, chi2[i] > 1e-3):
        mtot_s = (m1[i] + m2[i])*MTSUN_SI # total mass in seconds
        f_final = v_final*v_final*v_final/(mtot_s*np.pi)

        _, _, chi1x_v_data, chi1y_v_data, chi1z_v_data, chi2x_v_data, chi2y_v_data, chi2z_v_data, Lnx_v_data, Lny_v_data, Lnz_v_data, _, _, _ = lalsim.SimInspiralSpinTaylorPNEvolveOrbit(deltaT=dt*mtot_s, m1=m1[i]*MSUN_SI, m2=m2[i]*MSUN_SI, fStart=f_start[i], fEnd=f_final, s1x=chi1[i]*np.sin(tilt1[i]), s1y=0., s1z=chi1[i]*np.cos(tilt1[i]), s2x=chi2[i]*np.sin(tilt2[i])*np.cos(phi12[i]), s2y=chi2[i]*np.sin(tilt2[i])*np.sin(phi12[i]), s2z=chi2[i]*np.cos(tilt2[i]), lnhatx=0., lnhaty=0., lnhatz=1., e1x=1., e1y=0., e1z=0., lambda1=0., lambda2=0., quadparam1=1., quadparam2=1., spinO=7, tideO=0, phaseO=7, lscorr=0, approx=approx)

        # Set index to take from array output by lalsim.SimInspiralSpinTaylorPNEvolveOrbit: -1 for evolving forward in time and 0 for evolving backward in time
        if f_start[i] <= f_final:
            idx_use = -1
        else:
            idx_use = 0

        chi1_v = np.array([chi1x_v_data.data.data[idx_use], chi1y_v_data.data.data[idx_use], chi1z_v_data.data.data[idx_use]])
        chi2_v = np.array([chi2x_v_data.data.data[idx_use], chi2y_v_data.data.data[idx_use], chi2z_v_data.data.data[idx_use]])

        Ln_v = np.array([Lnx_v_data.data.data[idx_use], Lny_v_data.data.data[idx_use], Lnz_v_data.data.data[idx_use]])

        tilt1_evol[i], tilt2_evol[i], phi12_evol[i] = tilts_and_phi12_from_Cartesian_spins_and_L(chi1_v, chi2_v, Ln_v)
    else:
        tilt1_evol[i], tilt2_evol[i], phi12_evol[i] = tilt1[i], tilt2[i], phi12[i]

# Setup for output

# Get the names of the columns
colnames = np.array(data.colnames)

tilt1_id = np.where(colnames=='tilt1')[0][0]
tilt2_id = np.where(colnames=='tilt2')[0][0]
phi12_id = np.where(colnames=='phi12')[0][0]
f_ref_id = np.where(colnames=='f_ref')[0][0]

# Get the meta for tilt1, tilt2, phi12
meta_tilt1 = tuple(data.columns.items())[tilt1_id][1].meta
meta_tilt2 = tuple(data.columns.items())[tilt2_id][1].meta
meta_phi12 = tuple(data.columns.items())[phi12_id][1].meta
meta_f_ref = tuple(data.columns.items())[f_ref_id][1].meta

# Create an astropy table with the evolved spin samples
tilts_evol = Table([Column(tilt1_evol, name='tilt1' + label, meta=meta_tilt1), Column(tilt2_evol, name='tilt2' + label, meta=meta_tilt2), Column(phi12_evol, name='phi12' + label, meta=meta_phi12)])

# Append the columns to the existing astropy table of samples from the HDF5 file
data_joined = hstack([data, tilts_evol])

if v_final != "ISCO":
    vfinal_col = Table([Column(v_final*np.ones_like(tilt1_evol), name='vfinal', meta=meta_f_ref)])
    data_joined = hstack([data_joined, vfinal_col])

_remap_colnames(data_joined)

f = h5py.File(hdf_pos_file, 'r')
path = '/lalinference/'+run_identifier+'/posterior_samples'
level = f[path]
arrt = level.attrs
names = np.array([list(arrt.items())[i][0] for i in range(len(list(arrt.items())))])
num_names = 0
for name in names:
    if 'NAME' in name:
        num_names += 1

# Updating the metadata
metadata[path]['FIELD_{0}_NAME'.format(num_names+1)]='tilt_spin1' + label
metadata[path]['FIELD_{0}_NAME'.format(num_names+2)]='tilt_spin2' + label
metadata[path]['FIELD_{0}_NAME'.format(num_names+3)]='phi12' + label

# Write the joined table
write_samples(data_joined, hdf_pos_file, path=path, metadata=metadata, overwrite=True)

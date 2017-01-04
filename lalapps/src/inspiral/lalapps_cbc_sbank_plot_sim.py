# Copyright (C) 2011  Nickolas Fotopoulos, Stephen Privitera
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

from __future__ import division
from math import ceil, floor, log10
import sys
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")

from lalsimulation import SimInspiralTaylorF2ReducedSpinComputeChi, SimIMRPhenomBComputeChi, SimIMRSEOBNRv4ROMTimeOfFrequency

from lal import MSUN_SI
from lalinspiral.sbank.waveforms import compute_mchirp
from lalinspiral.sbank.tau0tau3 import m1m2_to_tau0tau3

from h5py import File as H5File

def find_process_id(process_table, prog_name):
    for program, pid in process_table["program", "process_id"]:
        if program == prog_name:
            return pid
    return None

def find_process_param(pp_table, process_id, param):
    for pid, row_param, value in pp_table["process_id", "param", "value"]:
        if pid == process_id and row_param == param:
            return value

# grab data
match = None
for fname in sys.argv[1:]:
    with H5File(fname, "r") as h5file:

        if match is None:
            match = np.array(h5file["/match_map"])
            inj_arr = np.array(h5file["/sim_inspiral"])
            tmplt_arr = np.array(h5file["/sngl_inspiral"])

        else:
            match = np.concatenate((np.array(h5file["/match_map"]), match))
            tmplt_arr = np.concatenate((np.array(h5file["/sngl_inspiral"]), tmplt_arr))
            inj_arr = np.concatenate((np.array(h5file["/sim_inspiral"]), inj_arr))

inj_sigmasq = match["inj_sigmasq"]
match = match["match"]

order = match.argsort()[::-1]
match = match[order]
inj_arr = inj_arr[order]
tmplt_arr = tmplt_arr[order]

inj_approx = "FIXME"
tmplt_approx = "FIXME"
flow = 30 # FIXME

# derived quantities
# NB: Named fields in dtypes (see ligolw_table_to_array in lalapps_cbc_sbank_sim)
# allow dict-like access to columns: e.g. arr["mass1"].
min_match = match[int(floor(len(match)*0.9))]
inj_M = inj_arr["mass1"] + inj_arr["mass2"] #total mass
inj_eta = inj_arr["mass1"] * inj_arr["mass2"] / (inj_arr["mass1"] + inj_arr["mass2"])**2
tmplt_M = tmplt_arr["mass1"] + tmplt_arr["mass2"]
tmplt_eta = tmplt_arr["mass1"] * tmplt_arr["mass2"] / (tmplt_arr["mass1"] + tmplt_arr["mass2"])**2
inj_tau0, inj_tau3 = m1m2_to_tau0tau3(inj_arr["mass1"], inj_arr["mass2"], flow)
tmplt_tau0, tmplt_tau3 = m1m2_to_tau0tau3(tmplt_arr["mass1"], tmplt_arr["mass2"], flow)
tmplt_mchirp = compute_mchirp(tmplt_arr["mass1"], tmplt_arr["mass2"])
inj_mchirp = compute_mchirp(inj_arr["mass1"], inj_arr["mass2"])

#
# compute effective/reduced spin parameter for templates
#
if tmplt_approx == "TaylorF2RedSpin":
    tmplt_chi_label = "\chi_\mathrm{red}"
    tmplt_chi = [SimInspiralTaylorF2ReducedSpinComputeChi(float(row["mass1"]), float(row["mass2"]), float(row["spin1z"]), float(row["spin2z"])) for row in tmplt_arr]
else: # default to effective spin
    tmplt_chi_label = "\chi_\mathrm{eff}"
    tmplt_chi = [SimIMRPhenomBComputeChi(float(row["mass1"]), float(row["mass2"]), float(row["spin1z"]), float(row["spin2z"])) for row in tmplt_arr]

# FIXME hack
inj_dur = [SimIMRSEOBNRv4ROMTimeOfFrequency(flow, float(row["mass1"]) * MSUN_SI, float(row["mass2"]) * MSUN_SI, float(row["spin1z"]), float(row["spin2z"])) for row in inj_arr]


#
# compute effective/reduced spin parameter for injections
#
if inj_approx == "TaylorF2RedSpin" or "SpinTaylorT5" in inj_approx: # reduced spin
    inj_chi_label = "\chi_\mathrm{red}"
    inj_chi = [SimInspiralTaylorF2ReducedSpinComputeChi(float(row["mass1"]), float(row["mass2"]), float(row["spin1z"]), float(row["spin2z"])) for row in inj_arr]
elif "SpinTaylorT4" in inj_approx: # reduced spin (requires coordinate transformation)
    inj_chi_label = "\chi_\mathrm{red}"
    chi1 = [row["spin1x"] *np.sin(row["inclination"]) + row["spin1z"] *np.cos(row["inclination"]) for row in inj_arr]
    chi2 = [row["spin2x"] *np.sin(row["inclination"]) + row["spin2z"] *np.cos(row["inclination"]) for row in inj_arr]
    inj_chi = [SimIMRPhenomBComputeChi(float(row["mass1"]), float(row["mass2"]), s1z, s2z) for row, s1z, s2z in zip(inj_arr, chi1, chi2)]
else: # default to effective spin
    inj_chi_label = "\chi_\mathrm{eff}"
    inj_chi = [SimIMRPhenomBComputeChi(float(row["mass1"]), float(row["mass2"]), float(row["spin1z"]), float(row["spin2z"])) for row in inj_arr]

    # IMRPhenomP uses yet another coordinate convention
    if inj_approx == "IMRPhenomP":
        inj_chi *= np.cos(inj_arr["inclination"])


smallest_match = match.min()

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

name, ext =  os.path.splitext(os.path.basename(fname))

#
# binned, weighted (effective) fitting factor plots
#
fig = Figure()
ax = fig.add_subplot(111)
vals, xedges, yedges = np.histogram2d(inj_arr["mass1"], inj_arr["mass2"], bins=(10,10), weights= match**3 * inj_sigmasq**(3./2))
norm, xedges, yedges = np.histogram2d(inj_arr["mass1"], inj_arr["mass2"], bins=(10,10), weights=inj_sigmasq**(3./2))
vals = ( vals / np.array(norm + 1e-10, dtype="float") )**(1./3)
x, y= np.meshgrid(xedges, yedges)
coll = ax.pcolor(x, y, vals.T, vmin=0.9*vals[vals>0].min(), vmax=1)
ax.set_xlabel("$m_1$ ($M_\odot$)")
ax.set_ylabel("$m_2$ ($M_\odot$)")
ax.set_xticks((xedges[:-1] + xedges[1:])/2)
ax.set_yticks((yedges[:-1] + yedges[1:])/2)
ax.set_xticklabels(["%.1f" % round(k,1) for k in ax.get_xticks()])
ax.set_yticklabels(["%.1f" % round(k,1) for k in ax.get_yticks()])
ax.set_xlim([xedges.min(), xedges.max()])
ax.set_ylim([yedges.min(), yedges.max()])
fig.colorbar(coll, ax=ax).set_label("Effective Fitting Factor")
ax.grid(True)
canvas = FigureCanvas(fig)
fig.savefig(name+"_match_vs_injm1_vs_injm2_binned_lum_weighted.png")


fig = Figure()
ax = fig.add_subplot(111)
vals, xedges, yedges = np.histogram2d(inj_M, inj_eta, bins=(6,6), weights= match**3 * inj_sigmasq**(3./2))
norm, xedges, yedges = np.histogram2d(inj_M, inj_eta, bins=(6,6), weights=inj_sigmasq**(3./2))
vals = ( vals / np.array(norm + 1e-10, dtype="float") )**(1./3)
x, y= np.meshgrid(xedges, yedges)
coll = ax.pcolor(x, y, vals.T, vmin=0.9*vals[vals>0].min(), vmax=1)
ax.set_xlabel("$M_\mathrm{total}$ ($M_\odot$)")
ax.set_ylabel("$\eta$")
ax.set_xticks((xedges[:-1] + xedges[1:])/2)
ax.set_yticks((yedges[:-1] + yedges[1:])/2)
ax.set_xticklabels(["%.1f" % round(k,1) for k in ax.get_xticks()])
ax.set_yticklabels(["%.3f" % round(k,3) for k in ax.get_yticks()])
ax.set_xlim([xedges.min(), xedges.max()])
ax.set_ylim([yedges.min(), yedges.max()])
fig.colorbar(coll, ax=ax).set_label("Effective Fitting Factor")
ax.grid(True)
canvas = FigureCanvas(fig)
fig.savefig(name+"_match_vs_injM_vs_injeta_binned_lum_weighted.png")


#
# binned (unweighted) fitting factor plots
#
fig = Figure()
ax = fig.add_subplot(111)
coll = ax.hexbin(inj_M, inj_chi, C=match, gridsize=50, vmin=smallest_match, vmax=1)
ax.set_xlabel("Total Mass ($M_\odot$)")
ax.set_ylabel("$%s$" % inj_chi_label)
ax.set_xlim([min(inj_M), max(inj_M)])
ax.set_ylim([min(inj_chi), max(inj_chi)])
fig.colorbar(coll, ax=ax).set_label("Mean Fitting Factor")
ax.grid(True)
canvas = FigureCanvas(fig)
fig.savefig(name+"_match_vs_injM_vs_injchi_hexbin.png")


fig = Figure()
ax = fig.add_subplot(111)
coll = ax.hexbin(inj_arr["mass1"], inj_arr["mass2"], C=match, gridsize=50, vmin=smallest_match, vmax=1)
ax.set_xlabel("$m_1$ ($M_\odot$)")
ax.set_ylabel("$m_2$ ($M_\odot$)")
ax.set_xlim([min(inj_arr["mass1"]), max(inj_arr["mass1"])])
ax.set_ylim([min(inj_arr["mass2"]), max(inj_arr["mass2"])])
fig.colorbar(coll, ax=ax).set_label("Mean Fitting Factor")
ax.grid(True)
canvas = FigureCanvas(fig)
fig.savefig(name+"_match_vs_injm1_vs_injm2_hexbin.png")

#
# luminosity vs FF
#
fig = Figure()
ax = fig.add_subplot(111)
coll = ax.hexbin(inj_sigmasq**(0.5)/(8*inj_mchirp**(5./6)), inj_mchirp, gridsize=25, C=match)
ax.set_xlabel("SNR=8 Horizon Chirp Distance (Mpc)")
ax.set_ylabel("$\mathcal{M}_\mathrm{chirp}$")
fig.colorbar(coll, ax=ax).set_label("Mean Fitting Factor")
ax.grid(True)
canvas = FigureCanvas(fig)
fig.savefig(name+"_match_vs_injsigsq_vs_mc.png")


fig = Figure()
ax = fig.add_subplot(111)
ax.plot(inj_arr["mass1"], match, marker = '.', linestyle="None")
ax.set_xlabel("Injection $m_1$ ($M_\odot$)")
ax.set_ylabel("Fitting Factor Between Injection and Template Bank")
ax.grid(True)
canvas = FigureCanvas(fig)
fig.savefig(name+"_match_vs_injm1.png")

fig = Figure()
ax = fig.add_subplot(111)
ax.plot(inj_arr["mass2"], match, marker = ".", linestyle="None")
ax.set_xlabel("Injection $m_2$ ($M_\odot$)")
ax.set_ylabel("Fitting Factor Between Injection and Template Bank")
ax.grid(True)
canvas = FigureCanvas(fig)
fig.savefig(name+"_match_vs_injm2.png")

fig = Figure()
ax = fig.add_subplot(111)
ax.plot(inj_M, match, marker = ".",linestyle="None")
ax.set_xlabel("Injection $M$ ($M_\odot$)")
ax.set_ylabel("Fitting Factor Between Injection and Template Bank")
ax.grid(True)
canvas = FigureCanvas(fig)
fig.savefig(name+"_match_vs_injM.png")

fig = Figure()
ax = fig.add_subplot(111)
ax.plot(inj_mchirp, match, marker = ".", linestyle="None")
ax.set_xlabel("Injection $M_c$ ($M_\odot$)")
ax.set_ylabel("Fitting Factor Between Injection and Template Bank")
ax.grid(True)
canvas = FigureCanvas(fig)
fig.savefig(name+"_match_vs_injMc.png")

fig = Figure()
ax = fig.add_subplot(111)
ax.plot(inj_dur, match, marker = ".", linestyle="None")
ax.axvline(x=0.2, linewidth=2, color='k')
ax.set_xlabel("Injection Duration (s)")
ax.set_ylabel("Fitting Factor Between Injection and Template Bank")
ax.grid(True)
canvas = FigureCanvas(fig)
fig.savefig(name+"_match_vs_injdur.png")

fig = Figure()
ax = fig.add_subplot(111)
ax.plot(inj_eta, match , marker = ".", linestyle="None")
ax.set_xlabel("Injection $\eta$")
ax.set_ylabel("Match Between Injection and Template Bank")
ax.set_xlim([inj_eta.min(), inj_eta.max()])
ax.grid(True)
canvas = FigureCanvas(fig)
fig.savefig(name+"_match_vs_injeta.png")

fig = Figure()
ax = fig.add_subplot(111)
collection = ax.scatter(inj_arr["mass1"], tmplt_arr["mass1"], c=match, s=20, vmin=smallest_match, linewidth=0, alpha=0.5, vmax=1)
ax.set_xlabel("Injection $m_1$")
ax.set_ylabel("Best Matching Template $m_1$ ($M_\odot$)")
ax.grid(True)
fig.colorbar(collection, ax=ax).set_label("Fitting Factor")
canvas = FigureCanvas(fig)
fig.savefig(name+"_match_vs_tmpltm1_vs_injm1.png")

fig = Figure()
ax = fig.add_subplot(111)
collection = ax.scatter(inj_arr["mass2"], tmplt_arr["mass2"], c=match, s=20, vmin=smallest_match, linewidth=0, alpha=0.5, vmax=1)
ax.set_xlabel("Injection $m_2$")
ax.set_ylabel("Best Matching Template $m_2$ ($M_\odot$)")
ax.grid(True)
fig.colorbar(collection, ax=ax).set_label("Fitting Factor")
canvas = FigureCanvas(fig)
fig.savefig(name+"_match_vs_tmpltm2_vs_injm2.png")

fig = Figure()
ax = fig.add_subplot(111)
collection = ax.scatter(inj_chi, tmplt_chi, c=match, s=20, vmin=smallest_match, linewidth=0, alpha=0.5, vmax=1)
ax.set_xlabel("Injection $%s$" % inj_chi_label)
ax.set_ylabel("Best Matching Template $%s$" % tmplt_chi_label)
ax.grid(True)
fig.colorbar(collection, ax=ax).set_label("Fitting Factor")
canvas = FigureCanvas(fig)
fig.savefig(name+"_match_vs_tmpltchi_vs_injchi.png")

fig = Figure()
ax = fig.add_subplot(111)
collection = ax.scatter(inj_eta, tmplt_eta, c=match, s=20, vmin=smallest_match, linewidth=0, alpha=0.5, vmax=1)
ax.set_xlabel("Injection $\eta$")
ax.set_ylabel("Best Matching Template $eta$")
ax.set_xlim([inj_eta.min(), inj_eta.max()])
ax.set_ylim([tmplt_eta.min(), tmplt_eta.max()])
ax.grid(True)
fig.colorbar(collection, ax=ax).set_label("Fitting Factor")
canvas = FigureCanvas(fig)
fig.savefig(name+"_match_vs_injeta_vs_tmplteta.png")

fig = Figure()
ax = fig.add_subplot(111)
collection = ax.scatter(inj_M, tmplt_M, c=match, s=20, vmin=smallest_match, linewidth=0, alpha=0.5, vmax=1)
ax.set_xlabel("Injection total mass, $M$ ($M_\odot$)")
ax.set_ylabel("Best Matching Template $M$ ($M_\odot$)")
ax.grid(True)
fig.colorbar(collection, ax=ax).set_label("Fitting Factor")
canvas = FigureCanvas(fig)
fig.savefig(name+"_match_vs_injM_vs_tmpltM.png")

fig = Figure()
ax = fig.add_subplot(111)
collection = ax.scatter(inj_mchirp, tmplt_mchirp, c=match, s=20, vmin=smallest_match, linewidth=0, alpha=0.5, vmax=1)
ax.set_xlabel("Injection $\mathcal{M}_\mathrm{chirp}$ ($M_\odot$)")
ax.set_ylabel("Best Matching Template $\mathcal{M}_\mathrm{chirp}$ ($M_\odot$)")
ax.grid(True)
fig.colorbar(collection, ax=ax).set_label("Fitting Factor")
canvas = FigureCanvas(fig)
fig.savefig(name+"_match_vs_injmc_vs_tmpltmc.png")

fig = Figure()
ax = fig.add_subplot(111)
collection = ax.scatter(tmplt_arr["mass1"], tmplt_arr["mass2"], c=match, s=20, vmin=smallest_match, linewidth=0, alpha=0.5, vmax=1)
ax.set_xlabel("Template $m_1$ ($M_\odot$)")
ax.set_ylabel("Template $m_2$ ($M_\odot$)")
ax.grid(True)
fig.colorbar(collection, ax=ax).set_label("Fitting Factor")
canvas = FigureCanvas(fig)
fig.savefig(name+"_match_vs_tmpltm2_vs_tmpltm1.png")

fig = Figure()
ax = fig.add_subplot(111)
collection = ax.scatter(inj_arr["mass1"], inj_arr["mass2"], c=match, s=20, vmin=smallest_match, linewidth=0, alpha=0.5, vmax=1)
ax.set_xlabel("Injected $m_1$ ($M_\odot$)")
ax.set_ylabel("Injected $m_2$ ($M_\odot$)")
ax.grid(True)
fig.colorbar(collection, ax=ax).set_label("Fitting Factor")
canvas = FigureCanvas(fig)
fig.savefig(name+"_match_vs_injm2_vs_injm1.png")

fig = Figure()
ax = fig.add_subplot(111)
collection = ax.scatter(inj_arr["spin1z"], inj_arr["spin2z"], c=match, s=20, vmin=smallest_match, linewidth=0, alpha=0.5, vmax=1)
ax.set_xlabel("Injected $s_{1z}$")
ax.set_ylabel("Injected $s_{2z}$")
ax.grid(True)
fig.colorbar(collection, ax=ax).set_label("Fitting Factor")
canvas = FigureCanvas(fig)
fig.savefig(name+"_match_vs_injs1z_vs_injs2z.png")

fig = Figure()
ax = fig.add_subplot(111)
collection = ax.scatter(inj_tau0, inj_tau3, c=match, s=20, vmin=smallest_match, linewidth=0, alpha=0.5, vmax=1)
ax.set_xlabel(r"Injected $\tau_0$ (s)")
ax.set_ylabel(r"Injected $\tau_3$ (s)")
ax.set_title(r"Colorbar is Fitting Factor; assuming $f_\mathrm{low}=%d\,\mathrm{Hz}$" % flow)
ax.grid(True)
fig.colorbar(collection, ax=ax).set_label("Fitting Factor")
canvas = FigureCanvas(fig)
fig.savefig(name+"_match_vs_injtau0_vs_injtau3.png")

fig = Figure()
ax = fig.add_subplot(111)
collection = ax.scatter(inj_tau0, inj_chi, c=match, s=20, vmin=smallest_match, linewidth=0, alpha=0.5, vmax=1)
ax.set_xlabel(r"Injected $\tau_0$ (s)")
ax.set_ylabel(r"Injected $%s$" % inj_chi_label)
ax.set_title(r"Colorbar is Fitting Factor; assuming $f_\mathrm{low}=%d\,\mathrm{Hz}$" % flow)
ax.grid(True)
fig.colorbar(collection, ax=ax).set_label("Fitting Factor")
canvas = FigureCanvas(fig)
fig.savefig(name+"_match_vs_injtau0_vs_injchi.png")

fig = Figure()
ax = fig.add_subplot(111)
collection = ax.scatter(tmplt_tau0, tmplt_tau3, c=match, s=20, vmin=smallest_match, linewidth=0, alpha=0.5, vmax=1)
ax.set_xlabel(r"Template $\tau_0$ (s)")
ax.set_ylabel(r"Template $\tau_3$ (s)")
ax.set_title(r"Colorbar is Fitting Factor; assuming $f_\mathrm{low}=%d\,\mathrm{Hz}$" % flow)
ax.grid(True)
fig.colorbar(collection, ax=ax).set_label("Fitting Factor")
canvas = FigureCanvas(fig)
fig.savefig(name+"_match_vs_tmplttau0_vs_tmplttau3.png")

fig = Figure()
ax = fig.add_subplot(111)
collection = ax.scatter(tmplt_tau0, tmplt_chi, c=match, s=20, vmin=smallest_match, linewidth=0, alpha=0.5, vmax=1)
ax.set_xlabel(r"Template $\tau_0$ (s)")
ax.set_ylabel(r"Template $\chi$ (s)")
ax.set_title(r"Colorbar is Fitting Factor; assuming $f_\mathrm{low}=%d\,\mathrm{Hz}$" % flow)
ax.grid(True)
fig.colorbar(collection, ax=ax).set_label("Fitting Factor")
canvas = FigureCanvas(fig)
fig.savefig(name+"_match_vs_tmplttau0_vs_tmpltchi.png")

fig = Figure()
ax = fig.add_subplot(111)
collection = ax.scatter(tmplt_M, tmplt_chi, c=match, s=20, vmin=smallest_match, linewidth=0, alpha=0.5, vmax=1)
ax.set_xlabel(r"Template $M_{total}$ ($M_\odot$)")
ax.set_ylabel(r"Template $%s$" % tmplt_chi_label)
ax.set_title(r"Colorbar is Fitting Factor; assuming $f_\mathrm{low}=%d\,\mathrm{Hz}$" % flow)
ax.grid(True)
fig.colorbar(collection, ax=ax).set_label("Fitting Factor")
canvas = FigureCanvas(fig)
fig.savefig(name+"_match_vs_tmpltM_vs_tmpltchi.png")

fig = Figure()
ax = fig.add_axes((0.1, 0.1, 0.85, 0.85), xscale="log", yscale="log")
count, bin_edges = np.histogram(np.log(1 - match + 1e-5), bins=40) # fudge factor to avoid log(0)
bin_centers = np.exp(-(bin_edges[:-1] * bin_edges[1:])**0.5)
ax.plot(bin_centers, count, linewidth=2.5)
ax.plot((1 - min_match, 1 - min_match), (1, 10**ceil(log10(count.max()))), "k--")
ax.set_xlabel("mismatch (1 - fitting factor)")
ax.set_ylabel("Number")
ax.set_title(r'$N_\mathrm{inj}=%d$, 10th percentile=%.2f' % (len(match), min_match,))
ax.grid(True)
canvas = FigureCanvas(fig)
fig.savefig(name+"_match_hist.png")

fig = Figure()
ax = fig.add_subplot(111)
collection = ax.scatter(inj_M, inj_chi, c=match, s=20, vmin=smallest_match, linewidth=0, alpha=0.5, vmax=1)
ax.set_xlabel("Injected Total Mass")
ax.set_ylabel("Injected $%s$" % inj_chi_label)
ax.grid(True)
fig.colorbar(collection, ax=ax).set_label("Fitting Factor")
canvas = FigureCanvas(fig)
fig.savefig(name+"_match_vs_injM_vs_injchi.png")

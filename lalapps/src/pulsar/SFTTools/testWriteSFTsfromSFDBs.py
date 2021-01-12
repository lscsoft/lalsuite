# Copyright (C) 2020 Pep Covas, David Keitel, Rodrigo Tenorio
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

import os
import lalpulsar
import glob
import subprocess
import numpy as np
from numpy.testing import assert_allclose

# run this test with "make check TESTS=testWriteSFTsfromSFDBs"
# everything will be automatically put into testWriteSFTsfromSFDBs.testdir
testdir = "."


def validate_and_load_SFTs ( globstr, fmin, fmax ):
    sftnames = list(glob.glob(os.path.join(testdir,globstr)))
    assert len(sftnames) == 1
    sftname = sftnames[0]
    app = "lalapps_SFTvalidate"
    cl = app + " " + sftname
    print("Executing: " + cl)
    out = subprocess.check_output(cl, shell=True, stderr=subprocess.STDOUT, universal_newlines=True)
    print("{:s}: {}".format(app,out))
    catalog = lalpulsar.SFTdataFind(sftname, None)
    sfts = lalpulsar.LoadSFTs(catalog, fmin, fmax)
    return sfts, sftname


def plot_SFTs_comparison ( sfts1, sfts2, plotnamebase="plot", label1="SFTs1", label2="SFTs2", title="" ):
    minlen = min(sfts1.length,sfts2.length)
    if sfts1.length != sfts2.length:
        print("Inconsistent SFT vector lengths {:d}!={:d} for {:s}, {:s}. Will only plot the first {:d}.".format(sfts1.length,sfts2.length,label1,label2,minlen))
        # note: zip will do this automagically
    plotnamebase += "-"+label1+"-"+label2
    for k, (sft1,sft2) in enumerate(zip(sfts1.data,sfts2.data)):
        if sft1.epoch != sft2.epoch:
            raise ValueError("Inconsistent epochs {:d}!={:d} in SFT {:d} of {:s}, {:s}".format(sft1.epoch.gpsSeconds,sft2.epoch.gpsSeconds,k,label1,label2))
        title_this_epoch = title + ", epoch={:d}".format(sft1.epoch.gpsSeconds)
        plotnamebase_this_epoch = plotnamebase + "_{:d}".format(sft1.epoch.gpsSeconds)
        freqs1 = sft1.f0 + sft1.deltaF * np.arange(0,sft1.data.length)
        freqs2 = sft2.f0 + sft2.deltaF * np.arange(0,sft2.data.length)
        sft1_abs = np.abs(sft1.data.data)
        sft2_abs = np.abs(sft2.data.data)
        plt.plot(freqs1, sft1_abs, label=label1, color="black")
        plt.plot(freqs2, sft2_abs, label=label2, color="red")
        plt.xlabel("frequency [Hz]")
        plt.ylabel("|complex SFT bin|")
        plt.title(title_this_epoch)
        plt.legend()
        plotfile = os.path.join(testdir,plotnamebase_this_epoch+"-comparison.png")
        print("Plotting to " + plotfile)
        plt.savefig(plotfile)
        plt.close()
        plt.plot(freqs1, np.abs(sft1.data.data-sft2.data.data), color="black")
        plt.xlabel("frequency [Hz]")
        plt.ylabel("|{:s}-{:s}|".format(label1,label2))
        plt.title(title_this_epoch)
        plotfile = os.path.join(testdir,plotnamebase_this_epoch+"-absdiff.png")
        print("Plotting to " + plotfile)
        plt.savefig(plotfile)
        plt.close()
        plt.plot(freqs1, np.abs(sft1.data.data-sft2.data.data)/sft2_abs, color="black")
        plt.xlabel("frequency [Hz]")
        plt.ylabel("|{:s}-{:s}|/|{:s}|".format(label1,label2,label2))
        plt.title(title_this_epoch)
        plotfile = os.path.join(testdir,plotnamebase_this_epoch+"-reldiff.png")
        print("Plotting to " + plotfile)
        plt.savefig(plotfile)
        plt.close()
        plt.plot(sft2_abs, np.abs(sft1.data.data-sft2.data.data), ".", color="black")
        plt.xlabel("|{:s}|".format(label2))
        plt.ylabel("|{:s}-{:s}|".format(label1,label2))
        plt.title(title_this_epoch)
        plotfile = os.path.join(testdir,plotnamebase_this_epoch+"-absdiff-vs-abs2.png")
        print("Plotting to " + plotfile)
        plt.savefig(plotfile)
        plt.close()
        plt.plot(sft2_abs, np.abs(sft1.data.data-sft2.data.data)/sft2_abs, ".", color="black")
        plt.xlabel("|{:s}|".format(label2))
        plt.ylabel("|{:s}-{:s}|/|{:s}|".format(label1,label2,label2))
        plt.title(title_this_epoch)
        plotfile = os.path.join(testdir,plotnamebase_this_epoch+"-reldiff-vs-abs2.png")
        print("Plotting to " + plotfile)
        plt.savefig(plotfile)
        plt.close()


# files from reference tarball
ref_gwf = "V-V1_mfdv5-1257741529-2048.gwf"
ref_sft = "V-3_V1_1024SFT_mfdv5-1257741529-2048.sft"
ref_sfdb = "V1-mfdv5_20191114_043831.SFDB09"

# data parameters
IFOs = ["V1"]
fmin = 49.0
fmax = 51.0
startTime = 1257741529
frameDuration = 2048
# the following are only for the SFDBs and the SFTs we'll make to compare against,
# ref_sft has Tsft=frameDuration and Toverlap=0 instead!
Tsft = int(ref_sft.split("SFT")[0].split("_")[-1])
Toverlap = int(0.5*Tsft)


print("\n=== Running WriteSFTsfromSFDBs on reference SFDB file ===\n")
cl_SFDB = " ".join(["lalapps_WriteSFTsfromSFDBs",
                    "--file-pattern", ref_sfdb,
                    "--fmin", str(fmin),
                    "--fmax", str(fmax),
                    "--outSFTdir", testdir,
                   ])
print("Executing: " + cl_SFDB)
subprocess.check_call(cl_SFDB, shell=True)
sfts_from_sfdb, _ = validate_and_load_SFTs ( IFOs[0][0]+"*_from_SFDBs-*.sft", fmin, fmax )

# TODO: also test and validate --outSFTnames, --outSingleSFT=False modes

print("\n=== Reading reference SFTs from MFDv5 ===\n")
sfts_ref, _ = validate_and_load_SFTs ( ref_sft, fmin, fmax )
#catalog_ref = lalpulsar.SFTdataFind(ref_sft, None)
#sfts_ref = lalpulsar.LoadSFTs(catalog_ref, fmin, fmax)


print("=== Comparing outputs of WriteSFTsfromSFDBs and MFDv5: overall properties ===\n")
# compare overall properties
atol_freq = 1e-6
print("nSFTs: {} vs {}".format(sfts_from_sfdb.length,sfts_ref.length))  # Number of SFTs
print("Epochs: {} vs {}".format([sft.epoch for sft in sfts_from_sfdb.data], [sft.epoch for sft in sfts_ref.data]))
# lengths will not actually be the same:
# SFDBs contain a 4th SFT starting at 1257743065 and hanging over the end of the .gwf by half
if sfts_from_sfdb.length == sfts_ref.length + 1:
    print("Inconsistent SFT vector lengths {:d}=={:d}+1 from SFDBs vs referenc SFT. This is expected.".format(sfts_from_sfdb.length,sfts_ref.length))
assert sfts_from_sfdb.length == sfts_ref.length + 1
mean_sfts = np.zeros(sfts_ref.length) # will reuse this for test at the end
for k, (sft1,sft2) in enumerate(zip(sfts_from_sfdb.data,sfts_ref.data)):
    print("Comparing SFT number {:d}...".format(k))
    f0_sfdb = sft1.f0
    deltaF_sfdb = sft1.deltaF
    nBins_sfdb = sft1.data.length
    abs_sfdb = np.abs(sft1.data.data)
    max_sfdb = np.max(abs_sfdb)
    maxbin_sfdb = np.argmax(abs_sfdb)
    maxfreq_sfdb = f0_sfdb + deltaF_sfdb * maxbin_sfdb
    without_maxbin = np.arange(nBins_sfdb)!=maxbin_sfdb
    abs_without = abs_sfdb[without_maxbin]
    mean_sfdb = np.mean(abs_without)
    stdev_sfdb = mean_sfdb*np.std(abs_without/mean_sfdb)
    f0_sfts = sft2.f0
    deltaF_sfts = sft2.deltaF
    nBins_sfts = sft2.data.length
    abs_sfts = np.abs(sft2.data.data)
    max_sfts = np.max(abs_sfts)
    maxbin_sfts = np.argmax(abs_sfts)
    maxfreq_sfts = f0_sfts + deltaF_sfts * maxbin_sfts
    abs_without = abs_sfts[np.arange(nBins_sfts)!=maxbin_sfts]
    mean_sfts[k] = np.mean(abs_without)
    stdev_sfts = mean_sfts[k]*np.std(abs_without/mean_sfts[k])
    print("IFO:   {} vs {}".format(sft1.name,sft2.name))  # Name of IFO
    print("f0:    {} vs {}".format(f0_sfdb,f0_sfts))  # Starting frequency
    print("df:    {} vs {}".format(deltaF_sfdb,deltaF_sfts))  # Frequency spacing between bins
    print("nBins: {} vs {}".format(nBins_sfdb,nBins_sfts))  # Number of frequency bins in one SFT
    print("max in bin:  {} vs {}".format(maxbin_sfdb,maxbin_sfts))  # bin index for max of |SFT|
    print("freq of max: {:.4f} vs {:.4f}".format(maxfreq_sfdb, maxfreq_sfts))
    print("max |value|: {:.4e} vs {:.4e}".format(max_sfdb, max_sfts))
    print("background mean:  {:.4e} vs {:.4e}".format(mean_sfdb, mean_sfts[k]))
    print("background stdev: {:.4e} vs {:.4e}".format(stdev_sfdb, stdev_sfts))
    assert sft1.name == sft2.name
    assert abs(f0_sfdb-f0_sfts) < atol_freq
    assert abs(deltaF_sfdb-deltaF_sfts) < atol_freq
    assert nBins_sfdb == nBins_sfts
    assert maxbin_sfdb == maxbin_sfts
    assert abs(max_sfdb-max_sfts)/max_sfts < 0.01
    assert abs(maxfreq_sfdb-maxfreq_sfts) < atol_freq
    assert abs(mean_sfdb-mean_sfts[k])/mean_sfts[k] < 0.001
    assert abs(stdev_sfdb-stdev_sfts)/stdev_sfts < 0.001

print("\n=== Optional: comparison plots ===\n")
try:
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    plot_SFTs_comparison ( sfts_ref,
                           sfts_from_sfdb,
                           plotnamebase=IFOs[0],
                           label1="SFTs",
                           label2="SFDBs",
                           title="Tsft={:d}, Toverlap={:d}".format(Tsft,Toverlap)
                         )
except ImportError:
    print("Cannot import matplotlib to make comparison plots, continuing with numerical comparisons.")
except Exception as e:
    print("Failed to create comparison plots, continuing with numerical comparisons. Exception was:")
    print(e)


print("\n=== Comparing outputs of WriteSFTsfromSFDBs and MFDv5: actual per-bin data ===\n")
# ignoring here the overhanging 4th SFT from the SFDBs
# we use an atol threshold at 2% of the mean noise level
# because rtol tends to blow up in bins with very low noise
for k, (sft1,sft2) in enumerate(zip(sfts_from_sfdb.data,sfts_ref.data)):
    atol = 0.02*mean_sfts[k]
    maxdiff = max(np.abs(sft1.data.data-sft2.data.data))
    print("Comparing SFT number {:d}: max absolute difference is {:g} (allowed: {:g})".format(k,maxdiff,atol))
    assert_allclose(
        sft1.data.data,
        sft2.data.data,
        rtol=0,
        atol=atol,
        verbose=True,
    )

print("All SFDB checks passed (or skipped).")

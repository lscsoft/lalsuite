import sys
import os
from argparse import ArgumentParser

import numpy
import h5py

from glue.ligolw import utils, ligolw, lsctables
lsctables.use_in(ligolw.LIGOLWContentHandler)
from lalinference.rapid_pe import xmlutils
#import xmlutils

def append_lalinf_samples_to_hdf5_group(grp, samples):
    samp_grp = grp.create_group("samples")
    cols = set(samples.dtype.fields.keys())
    for col in cols:
        try:
            samples[col]
        except ValueError:
            continue

        ds = samp_grp.create_dataset(col, data=samples[col])

def obtain_tmplt_id(tmplt_bank, **kwargs):
    for t in tmplt_bank:
        if all([numpy.isclose(getattr(t, k), v) for k, v in kwargs.iteritems()]):
            return t.simulation_id
    raise ValueError("Cannot find template in bank")

VALID_PARAMS = set(("mass1", "mass2", "spin1x", "spin1y", "spin1z", "spin2x", "spin2y", "spin2z"))
def get_intr_prms_from_pp_table(table):
    intr_prms = {}
    for pp in table:
        if pp.param[2:] in VALID_PARAMS:
            intr_prms[pp.param[2:]] = float(pp.value)
    return intr_prms

KEY_MAPPING = {
    "snr": "logevidence",
    "ttotal": "n_samples",
    "tau0" : "n_eff",
    "event_duration": "int_var"
}

argp = ArgumentParser()
argp.add_argument("-o", "--output-file", help="Name of output HDF file. Will append if the file already exists.")
argp.add_argument("-t", "--tmplt-bank-file", help="Template bank / grid file used. Required.")
argp.add_argument("infiles", nargs="+", help="Files to convert and concatenate.")
args = argp.parse_args()

#outfname = infname.replace("xml.gz", "hdf")
#outfname = outfname.replace("xml", "hdf")

def rapidpe_to_hdf(basegrp, bankfile, samplefiles):
    # import (or look up template bank)
    #bank_xmldoc = utils.load_filename(args.tmplt_bank_file, contenthandler=ligolw.LIGOLWContentHandler)
    bank_xmldoc = utils.load_filename(bankfile, contenthandler=ligolw.LIGOLWContentHandler)
    try:
        tmplt_bank = lsctables.SimInspiralTable.get_table(bank_xmldoc)
    except:
        tmplt_bank = lsctables.SnglInspiralTable.get_table(bank_xmldoc)

    for sample_file in samplesfiles:
        print "Processing %s" % sample_file

        # Get tables from xmldoc
        xmldoc = utils.load_filename(sample_file, contenthandler=ligolw.LIGOLWContentHandler)
        process_params = lsctables.ProcessParamsTable.get_table(xmldoc)
        samples = lsctables.SimInspiralTable.get_table(xmldoc)
        run_result = lsctables.SnglInspiralTable.get_table(xmldoc)[0]

        # Get intrinsic parameter used
        intr_prms = get_intr_prms_from_pp_table(process_params)

        # Identify intrinisic ID
        tmplt_id = obtain_tmplt_id(tmplt_bank, **intr_prms)

        # Append sample data and metadata
        subgrp = base_grp.create_group(str(tmplt_id))
        xmlutils.append_samples_to_hdf5_group(subgrp, samples)

        def pp_table_to_dict(pptable):
            return dict([(pp.param, pp.value) for pp in pptable])

        run_info = pp_table_to_dict(process_params)
        for key in ("snr", "tau0", "ttotal"):
            run_info[KEY_MAPPING[key]] = getattr(run_result, key)

        xmlutils.append_metadata_to_hdf5_group(subgrp, run_info)

        xmldoc.unlink()

def txt_to_hdf(basegrp, samplefile):

    samplefile = samplefile[0]
    print "Processing %s" % samplefile

    samples = numpy.genfromtxt(samplefile, names=True)

    # Append sample data and metadata
    #subgrp = base_grp.create_group(str(tmplt_id))
    subgrp = base_grp.create_group("run_identifier")
    append_lalinf_samples_to_hdf5_group(subgrp, samples)

    return
    def pp_table_to_dict(pptable):
        return dict([(pp.param, pp.value) for pp in pptable])

    run_info = pp_table_to_dict(process_params)
    for key in ("snr", "tau0", "ttotal"):
        run_info[KEY_MAPPING[key]] = getattr(run_result, key)

    xmlutils.append_metadata_to_hdf5_group(subgrp, run_info)

    xmldoc.unlink()

pipeline = "lalinference"
if pipeline == "rapidpe":
    if os.path.exists(args.output_file):
        hfile = h5py.File(args.output_file, "a")
        base_grp = hfile["rapidpe"]
    else:
        hfile = h5py.File(args.output_file, "w")
        base_grp = hfile.create_group("rapidpe")
    rapidpe_to_hdf(basegrp, args.tmplt_bank_file, args.infiles)

elif pipeline == "lalinference":
    if os.path.exists(args.output_file):
        hfile = h5py.File(args.output_file, "a")
        base_grp = hfile["lalinference"]
    else:
        hfile = h5py.File(args.output_file, "w")
        base_grp = hfile.create_group("lalinference")

    txt_to_hdf(base_grp, args.infiles)

hfile.close()

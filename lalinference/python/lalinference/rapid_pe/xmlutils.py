# Copyright (C) 2012 Chris Pankow
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

import types
import sqlite3
from collections import namedtuple, Iterable

import numpy

from glue.lal import LIGOTimeGPS
from glue.ligolw import ligolw, lsctables, table, ilwd
from glue.ligolw.utils import process

__author__ = "Chris Pankow <chris.pankow@ligo.org>"

def assign_id(row, i):
    row.simulation_id = ilwd.ilwdchar("sim_inspiral_table:sim_inspiral:%d" % i)

CMAP = { "right_ascension": "longitude",
    "longitude":"longitude",
    "latitude":"latitude",
    "declination": "latitude",
    "inclination": "inclination",
    "polarization": "polarization",
    "t_ref": lambda r, t: r.set_time_geocent(LIGOTimeGPS(float(t))),
    "coa_phase": "coa_phase",
    "distance": "distance",
    "mass1": "mass1",
    "mass2": "mass2",
    "lam_tilde": "psi0",
    "dlam_tilde": "psi3",
    "psi0": "psi0",
    "psi3": "psi3",
    "spin1z": "spin1z",
    "spin2z": "spin2z",
    # SHOEHORN ALERT
    "sample_n": assign_id,
    "alpha1":"alpha1",
    "alpha2":"alpha2",
    "alpha3":"alpha3",
    "loglikelihood": "alpha1",
    "joint_prior": "alpha2",
    "joint_s_prior": "alpha3"
}

# FIXME: Find way to intersect given cols with valid cols when making table.
# Otherwise, we'll have to add them manually and ensure they all exist
sim_valid_req = ["process_id", "simulation_id", "alpha1", "alpha2", "alpha3"]
sim_valid_ext = ["inclination", "longitude", "latitude", "polarization", "geocent_end_time", "geocent_end_time_ns", "coa_phase", "distance"]
sim_valid_int = ["mass1", "mass2", "spin1z", "spin2z", "psi0", "psi3"]
sngl_valid_cols = ["process_id", "event_id", "snr", "tau0", "tau3"]
multi_valid_cols = ["process_id", "event_id", "snr"]

def append_samples_to_xmldoc(xmldoc, sampdict):
    write_cols = set(sim_valid_ext + sim_valid_int) & set(sampdict.keys())
    write_cols = list(write_cols) + sim_valid_req
    try: 
        si_table = table.get_table(xmldoc, lsctables.SimInspiralTable.tableName)
        new_table = False
    # Warning: This will also get triggered if there is *more* than one table
    except ValueError:
        si_table = lsctables.New(lsctables.SimInspiralTable, write_cols)
        new_table = True
    
    keys = sampdict.keys()
    # Just in case the key/value pairs don't come out synchronized
    values = numpy.array([sampdict[k] for k in keys], object)
    
    # Flatten the keys
    keys = reduce(list.__add__, [list(i) if isinstance(i, tuple) else [i] for i in keys])

    # Get the process
    # FIXME: Assumed that only we have appended information
    procid = table.get_table(xmldoc, lsctables.ProcessTable.tableName)[-1].process_id
    
    # map the samples to sim inspiral rows
    # NOTE :The list comprehension is to preserve the grouping of multiple 
    # parameters across the transpose operation. It's probably not necessary,
    # so if speed dictates, it can be reworked by flattening before arriving 
    # here
    for vrow in numpy.array(zip(*[vrow_sub.T for vrow_sub in values]), dtype=numpy.object):
        #si_table.append(samples_to_siminsp_row(si_table, **dict(zip(keys, vrow.flatten()))))
        vrow = reduce(list.__add__, [list(i) if isinstance(i, Iterable) else [i] for i in vrow])
        si_table.append(samples_to_siminsp_row(si_table, **dict(zip(keys, vrow))))
        si_table[-1].process_id = procid

    if new_table:
        xmldoc.childNodes[0].appendChild(si_table)
    return xmldoc

def append_likelihood_result_to_xmldoc(xmldoc, loglikelihood, neff=0, converged=False,**cols):
    try: 
        si_table = table.get_table(xmldoc, lsctables.SnglInspiralTable.tableName)
        new_table = False
        # NOTE: MultiInspiralTable has no spin columns
        #si_table = table.get_table(xmldoc, lsctables.MultiInspiralTable.tableName)
    # Warning: This will also get triggered if there is *more* than one table
    except ValueError:
        si_table = lsctables.New(lsctables.SnglInspiralTable, sngl_valid_cols + cols.keys())
        new_table = True
        # NOTE: MultiInspiralTable has no spin columns
        #si_table = lsctables.New(lsctables.MultiInspiralTable, multi_valid_cols + cols.keys())

    # Get the process
    # FIXME: Assumed that only we have appended information
    procid = table.get_table(xmldoc, lsctables.ProcessTable.tableName)[-1].process_id
    
    # map the samples to sim inspiral rows
    si_table.append(likelihood_to_snglinsp_row(si_table, loglikelihood, neff, converged,**cols))
    si_table[-1].process_id = procid

    if new_table:
        xmldoc.childNodes[0].appendChild(si_table)

    return xmldoc

def samples_to_siminsp_row(table, colmap={}, **sampdict):
    row = table.RowType()
    row.simulation_id = table.get_next_id()
    for key, col in CMAP.iteritems():
        if key not in sampdict:
            continue
        if isinstance(col, types.FunctionType):
            col(row, sampdict[key])
        else:
            setattr(row, col, sampdict[key])

    return row

def likelihood_to_snglinsp_row(table, loglikelihood, neff=0, converged=False, **cols):
    row = table.RowType()
    row.event_id = table.get_next_id()
    for col in cols:
            setattr(row, col, cols[col])
    row.snr = loglikelihood
    row.tau0 = neff
    row.tau3 = int(converged)

    return row

def db_identify_param(db_fname, process_id, param):
    """
    Extract the event time for a given process ID. This may fail in the case that the event time was not given on the command line (rather in a pointer to a XML file)
    NOTE: This is definitely not the best way to do this.
    """

    cmd_prm = "--" + param.replace("_", "-")

    sql = """select value from process_params where process_id = "%s" and param = "%s" """ % (str(process_id), cmd_prm)

    try:
        connection = sqlite3.connect(db_fname)
        result = list(connection.execute(sql))[0][0]
    finally:
        connection.close()
    return result

def db_to_samples(db_fname, tbltype, cols):
    """
    Pull samples from db_fname and return object that resembles a row from an XML table.
    """
    if "geocent_end_time" in cols:
        cols.append("geocent_end_time_ns")

    # FIXME: Get columns from db
    #if cols is None:
        #colsspec = "*"
    #else:
    colsspec = ", ".join(cols)

    if tbltype == lsctables.SimInspiralTable:
        sql = """select %s from sim_inspiral""" % colsspec
    elif tbltype == lsctables.SnglInspiralTable:
        sql = """select %s from sngl_inspiral""" % colsspec
    else:
        raise ValueError("Don't know SQL for table %s" % tbltype.tableName)

    Sample = namedtuple("Sample", cols)
    samples = []

    try:
        connection = sqlite3.connect(db_fname)
        connection.row_factory = sqlite3.Row
        for row in connection.execute(sql):
            # FIXME: UGH!
            res = dict(zip(cols, row))
            if "geocent_end_time" in res.keys():
                res["geocent_end_time"] += res["geocent_end_time_ns"]*1e-9
            if "simulation_id" in res.keys():
                res["simulation_id"] = ilwd.ilwdchar(res["simulation_id"])
            if "process_id" in res.keys():
                res["process_id"] = ilwd.ilwdchar(res["process_id"])

            samples.append(Sample(**res))
    finally:
        connection.close()

    return samples

#
# HDF5 I/O
#
def append_samples_to_hdf5_group(grp, samples, compress='gzip'):
    if "samples" in grp:
        samp_grp = grp["samples"]
    else:
        samp_grp = grp.create_group("samples")

    cols = set(samples.validcolumns.keys()) & set(dir(samples[0]))
    for col in cols:
        try:
            getattr(samples[0], col)
        except AttributeError:
            continue

        if samples.validcolumns[col] in ("ilwd:char",):
            raw_dat = map(int, [getattr(samples[i], col) for i, row in enumerate(samples)])
        else:
            raw_dat = [getattr(samples[i], col) for i, row in enumerate(samples)]

        if col in samp_grp:
            ds = samp_grp[col]
            ds.resize((len(ds) + len(samples),))
            ds[len(ds) - len(raw_dat):len(ds)] = raw_dat
        else:
            ds = samp_grp.create_dataset(col, data=raw_dat, maxshape=(None,), compression=compress)


def append_metadata_to_hdf5_group(grp, metadata):
    for name, val in metadata.iteritems():
        grp.attrs[name] = val or ""

# TESTING
import sys
if __file__ == sys.argv[0]:
    import numpy

    # Not used yet
    del CMAP["int_var"]
    del CMAP["int_val"]
    del CMAP["sample_n"]

    # Reworked to resemble usage in pipeline
    del CMAP["mass1"]
    del CMAP["mass2"]
    CMAP[("mass1", "mass2")] = ("mass1", "mass2")
    ar = numpy.random.random((len(CMAP), 10))
    samp_dict = dict(zip(CMAP, ar))
    ar = samp_dict[("mass1", "mass2")]
    samp_dict[("mass1", "mass2")] = numpy.array([ar, ar])
    del CMAP[("mass1", "mass2")]
    CMAP["mass1"] = "mass1"
    CMAP["mass2"] = "mass2"

    samp_dict["samp_n"] = numpy.array(range(0,10))
    CMAP["sample_n"] = "sample_n"
    
    xmldoc = ligolw.Document()
    xmldoc.appendChild(ligolw.LIGO_LW())
    process.register_to_xmldoc(xmldoc, sys.argv[0], {})

    append_samples_to_xmldoc(xmldoc, samp_dict)

    def gaussian(x, mu=0, std=1):
        return 1/numpy.sqrt(numpy.pi*2)/std * numpy.exp(-(x-mu)**2/2/std**2)

    m1m, m2m = 1.4, 1.5
    m1, m2 = numpy.random.random(2000).reshape(2,1000)*1.0+1.0
    loglikes = [gaussian(m1i, m1m)*gaussian(m2i, m2m) for m1i, m2i in zip(m1, m2)]
    #loglikelihood = - 7.5**2/2
    #append_likelihood_result_to_xmldoc(xmldoc, loglikelihood, **{"mass1": 1.4, "mass2": 1.4, "ifos": "H1,L1,V1"})
    for m1i, m2i, loglikelihood in zip(m1, m2, loglikes):
        append_likelihood_result_to_xmldoc(xmldoc, loglikelihood, **{"mass1": m1i, "mass2": m2i})

    from glue.ligolw import utils
    utils.write_filename(xmldoc, "iotest.xml.gz", gz=True)

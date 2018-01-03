#!/usr/bin/python
import sys
from glue.lal import CacheEntry
from glue.ligolw import lsctables, utils
for filename in (CacheEntry(line).path for line in file(sys.argv[1])):
	xmldoc = utils.load_filename(filename, gz = (filename or "stdin").endswith(".gz"))
	try:
		lsctables.table.get_table(xmldoc, lsctables.SnglInspiralTable.tableName)
	except ValueError:
		xmldoc.childNodes[-1].appendChild(lsctables.New(lsctables.SnglInspiralTable, columns = ("process_id", "ifo", "search", "channel", "end_time", "end_time_ns", "end_time_gmst", "impulse_time", "impulse_time_ns", "template_duration", "event_duration", "amplitude", "eff_distance", "coa_phase", "mass1", "mass2", "mchirp", "mtotal", "eta", "kappa", "chi", "tau0", "tau2", "tau3", "tau4", "tau5", "ttotal", "psi0", "psi3", "alpha", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6", "beta", "f_final", "snr", "chisq", "chisq_dof", "bank_chisq", "bank_chisq_dof", "cont_chisq", "cont_chisq_dof", "sigmasq", "rsqveto_duration", "Gamma0", "Gamma1", "Gamma2", "Gamma3", "Gamma4", "Gamma5", "Gamma6", "Gamma7", "Gamma8", "Gamma9", "event_id")))
		utils.write_filename(filename, xmldoc, gz = (filename or "stdout").endswith(".gz"))

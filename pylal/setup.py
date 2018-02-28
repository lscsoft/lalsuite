# setup for pylal


import os
from misc import generate_vcs_info as gvcsi
from distutils.core import setup, Extension
from distutils.command import install
from distutils.command import build
from distutils.command import build_py
from distutils.command import sdist
from distutils import log
import subprocess
import sys
import time
from numpy.lib.utils import get_include as numpy_get_include


class PkgConfig(object):
	def __init__(self, names):
		def stripfirsttwo(string):
			return string[2:]
		self.libs = list(map(stripfirsttwo, os.popen("pkg-config --libs-only-l %s" % names).read().split()))
		self.libdirs = list(map(stripfirsttwo, os.popen("pkg-config --libs-only-L %s" % names).read().split()))
		self.incdirs = list(map(stripfirsttwo, os.popen("pkg-config --cflags-only-I %s" % names).read().split()))
		self.extra_cflags = os.popen("pkg-config --cflags-only-other %s" % names).read().split()

gsl_pkg_config = PkgConfig("gsl")
lal_pkg_config = PkgConfig("lal")
lalsupport_pkg_config = PkgConfig("lalsupport")
lalburst_pkg_config = PkgConfig("lalburst")
# FIXME:  works for GCC only!!!
lal_pkg_config.extra_cflags += ["-std=c99"]
lalframe_pkg_config = PkgConfig("lalframe")
libframe_pkg_config = PkgConfig("libframe")
lalmetaio_pkg_config = PkgConfig("lalmetaio")
lalsimulation_pkg_config = PkgConfig("lalsimulation")
lalinspiral_pkg_config = PkgConfig("lalinspiral")

def remove_root(path, root):
	if root:
		return os.path.normpath(path).replace(os.path.normpath(root), "")
	return os.path.normpath(path)

def write_build_info():
	"""
	Get VCS info from misc/generate_vcs_info.py and add build information.
	Substitute these into misc/git_version.py.in to produce
	pylal/git_version.py.
	"""
	vcs_info = gvcsi.generate_git_version_info()

	# determine current time and treat it as the build time
	build_date = time.strftime('%Y-%m-%d %H:%M:%S +0000', time.gmtime())

	# determine builder
	retcode, builder_name = gvcsi.call_out(('git', 'config', 'user.name'))
	if retcode:
		builder_name = "Unknown User"
	retcode, builder_email = gvcsi.call_out(('git', 'config', 'user.email'))
	if retcode:
		builder_email = ""
	builder = "%s <%s>" % (builder_name, builder_email)

	sed_cmd = ('sed',
		'-e', 's/@ID@/%s/' % vcs_info.id,
		'-e', 's/@DATE@/%s/' % vcs_info.date,
		'-e', 's/@BRANCH@/%s/' % vcs_info.branch,
		'-e', 's/@TAG@/%s/' % vcs_info.tag,
		'-e', 's/@AUTHOR@/%s/' % vcs_info.author,
		'-e', 's/@COMMITTER@/%s/' % vcs_info.committer,
		'-e', 's/@STATUS@/%s/' % vcs_info.status,
		'-e', 's/@BUILDER@/%s/' % builder,
		'-e', 's/@BUILD_DATE@/%s/' % build_date,
		'misc/git_version.py.in')

	# FIXME: subprocess.check_call becomes available in Python 2.5
	sed_retcode = subprocess.call(sed_cmd,
		stdout=open('pylal/git_version.py', 'w'))
	if sed_retcode:
		raise gvcsi.GitInvocationError


class pylal_build_py(build_py.build_py):
	def run(self):
		# If we are building from tarball, do not update git version.
		# PKG-INFO is inserted into the tarball by the sdist target.
		if not os.path.exists("PKG-INFO"):
			# create the git_version module
			try:
				write_build_info()
				log.info("Generated pylal/git_version.py")
			except gvcsi.GitInvocationError:
				if not os.path.exists("pylal/git_version.py"):
					log.error("Not in git checkout or cannot find git executable.")
					sys.exit(1)

		# resume normal build procedure
		build_py.build_py.run(self)

class pylal_install(install.install):
	def run(self):
		pylal_prefix = remove_root(self.prefix, self.root)

		# Hardcode a check for system-wide installation;
		# in this case, don't make the user-env scripts.
		if pylal_prefix == sys.prefix:
			self.distribution.data_files = []
			install.install.run(self)
			return

		# create the user env scripts
		if self.install_purelib == self.install_platlib:
			pylal_pythonpath = self.install_purelib
		else:
			pylal_pythonpath = self.install_platlib + ":" + self.install_purelib

		pylal_install_scripts = remove_root(self.install_scripts, self.root)
		pylal_pythonpath = remove_root(pylal_pythonpath, self.root)
		pylal_install_platlib = remove_root(self.install_platlib, self.root)

		if not os.path.isdir("etc"):
			os.mkdir("etc")
		log.info("creating pylal-user-env.sh script")
		env_file = open(os.path.join("etc", "pylal-user-env.sh"), "w")
		print >> env_file, "# Source this file to access PYLAL"
		print >> env_file, "PYLAL_PREFIX=" + pylal_prefix
		print >> env_file, "export PYLAL_PREFIX"
		if self.distribution.scripts:
			print >> env_file, "PATH=" + pylal_install_scripts + ":${PATH}"
			print >> env_file, "export PATH"
		print >> env_file, "PYTHONPATH=" + pylal_pythonpath + ":${PYTHONPATH}"
		print >> env_file, "LD_LIBRARY_PATH=" + pylal_install_platlib + ":${LD_LIBRARY_PATH}"
		print >> env_file, "DYLD_LIBRARY_PATH=" + pylal_install_platlib + ":${DYLD_LIBRARY_PATH}"
		print >> env_file, "export PYTHONPATH LD_LIBRARY_PATH DYLD_LIBRARY_PATH"
		env_file.close()

		log.info("creating pylal-user-env.csh script")
		env_file = open(os.path.join("etc", "pylal-user-env.csh"), "w")
		print >> env_file, "# Source this file to access PYLAL"
		print >> env_file, "setenv PYLAL_PREFIX " + pylal_prefix
		if self.distribution.scripts:
			print >> env_file, "setenv PATH " + pylal_install_scripts + ":${PATH}"
		print >> env_file, "if ( $?PYTHONPATH ) then"
		print >> env_file, "  setenv PYTHONPATH " + pylal_pythonpath + ":${PYTHONPATH}"
		print >> env_file, "else"
		print >> env_file, "  setenv PYTHONPATH " + pylal_pythonpath
		print >> env_file, "endif"
		print >> env_file, "if ( $?LD_LIBRARY_PATH ) then"
		print >> env_file, "  setenv LD_LIBRARY_PATH " + pylal_install_platlib + ":${LD_LIBRARY_PATH}"
		print >> env_file, "else"
		print >> env_file, "  setenv LD_LIBRARY_PATH " + pylal_install_platlib
		print >> env_file, "endif"
		print >> env_file, "if ( $?DYLD_LIBRARY_PATH ) then"
		print >> env_file, "  setenv DYLD_LIBRARY_PATH " + pylal_install_platlib + ":${DYLD_LIBRARY_PATH}"
		print >> env_file, "else"
		print >> env_file, "  setenv DYLD_LIBRARY_PATH " + pylal_install_platlib
		print >> env_file, "endif"
		env_file.close()

		# now run the installer
		install.install.run(self)


class pylal_sdist(sdist.sdist):
	def run(self):
		# customize tarball contents
		self.distribution.data_files = []
		for root,dirs,files in os.walk("debian"):
			for file in files:
				self.distribution.data_files += [os.path.join(root,file)]
		self.distribution.data_files += ["pylal.spec"]

		# create the git_version module
		try:
			write_build_info()
			log.info("generated pylal/git_version.py")
		except gvcsi.GitInvocationError:
			if not os.path.exists("pylal/git_version.py"):
				log.error("Not in git checkout or cannot find git executable. Exiting.")
				sys.exit(1)

		# now run sdist
		sdist.sdist.run(self)

# FIXME:  all occurances of -DPY_SSIZE_T_CLEAN are a temporary hack to
# induce PyArg_ParseTuple() to exhibit 3.x behaviour in 2.x.  Remove when
# no longer needed

setup(
	name = "pylal",
	version = "0.10.0",
	author = "Kipp Cannon and Nickolas Fotopoulos",
	author_email = "lal-discuss@ligo.org",
	description = "Python LIGO Algorithm Library",
	url = "http://www.lsc-group.phys.uwm.edu/daswg/",
	license = "See file LICENSE",
	packages = [
		"pylal",
		"pylal.xlal",
		"pylal.xlal.datatypes",
		"pylal.dq"
	],
	cmdclass = {
		"build_py": pylal_build_py,
		"install": pylal_install,
		"sdist": pylal_sdist
	},
	ext_modules = [
		Extension(
			"pylal.Fr",
			["src/Fr.c"],
			# Use lalframe headers to silence warnings but link against libframe
			include_dirs = lalframe_pkg_config.incdirs + [numpy_get_include()],
			libraries = libframe_pkg_config.libs,
			library_dirs = libframe_pkg_config.libdirs,
			runtime_library_dirs = libframe_pkg_config.libdirs,
			extra_compile_args = lalframe_pkg_config.extra_cflags
		),
		Extension(
			"pylal.tools",
			["src/tools.c"],
			include_dirs = lal_pkg_config.incdirs + lalmetaio_pkg_config.incdirs + lalinspiral_pkg_config.incdirs,
			libraries = lal_pkg_config.libs + lalinspiral_pkg_config.libs,
			library_dirs = lal_pkg_config.libdirs + lalinspiral_pkg_config.libdirs,
			runtime_library_dirs = lal_pkg_config.libdirs + lalinspiral_pkg_config.libdirs,
			extra_compile_args = lal_pkg_config.extra_cflags
		),
		Extension(
			"pylal.xlal.datatypes.complex16fftplan",
			["src/xlal/datatypes/complex16fftplan.c"],
			include_dirs = lal_pkg_config.incdirs + [numpy_get_include(), "src/xlal/datatypes"],
			libraries = lal_pkg_config.libs,
			library_dirs = lal_pkg_config.libdirs,
			runtime_library_dirs = lal_pkg_config.libdirs,
			extra_compile_args = lal_pkg_config.extra_cflags
		),
		Extension(
			"pylal.xlal.datatypes.complex16frequencyseries",
			["src/xlal/datatypes/complex16frequencyseries.c"],
			include_dirs = lal_pkg_config.incdirs + [numpy_get_include(), "src/xlal/datatypes"],
			libraries = lal_pkg_config.libs,
			library_dirs = lal_pkg_config.libdirs,
			runtime_library_dirs = lal_pkg_config.libdirs,
			extra_compile_args = lal_pkg_config.extra_cflags
		),
		Extension(
			"pylal.xlal.datatypes.complex16timeseries",
			["src/xlal/datatypes/complex16timeseries.c"],
			include_dirs = lal_pkg_config.incdirs + [numpy_get_include(), "src/xlal/datatypes"],
			libraries = lal_pkg_config.libs,
			library_dirs = lal_pkg_config.libdirs,
			runtime_library_dirs = lal_pkg_config.libdirs,
			extra_compile_args = lal_pkg_config.extra_cflags
		),
		Extension(
			"pylal.xlal.datatypes.lalunit",
			["src/xlal/datatypes/lalunit.c"],
			include_dirs = lal_pkg_config.incdirs + ["src/xlal/datatypes"],
			libraries = lal_pkg_config.libs,
			library_dirs = lal_pkg_config.libdirs,
			runtime_library_dirs = lal_pkg_config.libdirs,
			extra_compile_args = lal_pkg_config.extra_cflags
		),
		Extension(
			"pylal.xlal.datatypes.ligotimegps",
			["src/xlal/datatypes/ligotimegps.c"],
			include_dirs = lal_pkg_config.incdirs + ["src/xlal/datatypes"],
			libraries = lal_pkg_config.libs,
			library_dirs = lal_pkg_config.libdirs,
			runtime_library_dirs = lal_pkg_config.libdirs,
			extra_compile_args = lal_pkg_config.extra_cflags
		),
		Extension(
			"pylal.xlal.datatypes.real8fftplan",
			["src/xlal/datatypes/real8fftplan.c"],
			include_dirs = lal_pkg_config.incdirs + [numpy_get_include(), "src/xlal/datatypes"],
			libraries = lal_pkg_config.libs,
			library_dirs = lal_pkg_config.libdirs,
			runtime_library_dirs = lal_pkg_config.libdirs,
			extra_compile_args = lal_pkg_config.extra_cflags
		),
		Extension(
			"pylal.xlal.datatypes.real8frequencyseries",
			["src/xlal/datatypes/real8frequencyseries.c"],
			include_dirs = lal_pkg_config.incdirs + [numpy_get_include(), "src/xlal/datatypes"],
			libraries = lal_pkg_config.libs,
			library_dirs = lal_pkg_config.libdirs,
			runtime_library_dirs = lal_pkg_config.libdirs,
			extra_compile_args = lal_pkg_config.extra_cflags
		),
		Extension(
			"pylal.xlal.datatypes.real8timeseries",
			["src/xlal/datatypes/real8timeseries.c"],
			include_dirs = lal_pkg_config.incdirs + [numpy_get_include(), "src/xlal/datatypes"],
			libraries = lal_pkg_config.libs,
			library_dirs = lal_pkg_config.libdirs,
			runtime_library_dirs = lal_pkg_config.libdirs,
			extra_compile_args = lal_pkg_config.extra_cflags
		),
		Extension(
			"pylal.xlal.datatypes.simburst",
			["src/xlal/datatypes/simburst.c", "src/xlal/misc.c"],
			include_dirs = lal_pkg_config.incdirs + lalmetaio_pkg_config.incdirs + ["src/xlal", "src/xlal/datatypes"],
			libraries = lal_pkg_config.libs,
			library_dirs = lal_pkg_config.libdirs,
			runtime_library_dirs = lal_pkg_config.libdirs,
			extra_compile_args = lal_pkg_config.extra_cflags + ["-DPY_SSIZE_T_CLEAN"]
		),
		Extension(
			"pylal.xlal.datatypes.siminspiraltable",
			["src/xlal/datatypes/siminspiraltable.c", "src/xlal/misc.c"],
			include_dirs = lal_pkg_config.incdirs + lalmetaio_pkg_config.incdirs+ ["src/xlal", "src/xlal/datatypes"],
			libraries = lal_pkg_config.libs,
			library_dirs = lal_pkg_config.libdirs,
			runtime_library_dirs = lal_pkg_config.libdirs,
			extra_compile_args = lal_pkg_config.extra_cflags + ["-DPY_SSIZE_T_CLEAN"]
		),
		Extension(
			"pylal.xlal.datatypes.snglburst",
			["src/xlal/datatypes/snglburst.c", "src/xlal/misc.c"],
			include_dirs = lal_pkg_config.incdirs + lalmetaio_pkg_config.incdirs + ["src/xlal", "src/xlal/datatypes"],
			libraries = lal_pkg_config.libs,
			library_dirs = lal_pkg_config.libdirs,
			runtime_library_dirs = lal_pkg_config.libdirs,
			extra_compile_args = lal_pkg_config.extra_cflags + ["-DPY_SSIZE_T_CLEAN"]
		),
		Extension(
			"pylal.xlal.datatypes.snglinspiraltable",
			["src/xlal/datatypes/snglinspiraltable.c", "src/xlal/misc.c"],
			include_dirs = lal_pkg_config.incdirs + lalmetaio_pkg_config.incdirs + ["src/xlal", "src/xlal/datatypes"],
			libraries = lal_pkg_config.libs,
			library_dirs = lal_pkg_config.libdirs,
			runtime_library_dirs = lal_pkg_config.libdirs,
			extra_compile_args = lal_pkg_config.extra_cflags + ["-DPY_SSIZE_T_CLEAN"]
		),
		Extension(
			"pylal.xlal.datatypes.snglringdowntable",
			["src/xlal/datatypes/snglringdowntable.c", "src/xlal/misc.c"],
			include_dirs = lal_pkg_config.incdirs + lalmetaio_pkg_config.incdirs + ["src/xlal", "src/xlal/datatypes"],
			libraries = lal_pkg_config.libs,
			library_dirs = lal_pkg_config.libdirs,
			runtime_library_dirs = lal_pkg_config.libdirs,
			extra_compile_args = lal_pkg_config.extra_cflags
		),
		Extension(
			"pylal.xlal.date",
			["src/xlal/date.c", "src/xlal/misc.c"],
			include_dirs = lal_pkg_config.incdirs + [numpy_get_include(), "src/xlal"],
			libraries = lal_pkg_config.libs,
			library_dirs = lal_pkg_config.libdirs,
			runtime_library_dirs = lal_pkg_config.libdirs,
			extra_compile_args = lal_pkg_config.extra_cflags
		),
		Extension(
			"pylal.xlal.fft",
			["src/xlal/fft.c", "src/xlal/misc.c"],
			include_dirs = lal_pkg_config.incdirs + ["src/xlal"],
			libraries = lal_pkg_config.libs,
			library_dirs = lal_pkg_config.libdirs,
			runtime_library_dirs = lal_pkg_config.libdirs,
			extra_compile_args = lal_pkg_config.extra_cflags
		),
		Extension(
			"pylal.xlal.lalburst",
			["src/xlal/lalburst.c", "src/xlal/misc.c"],
			include_dirs = lal_pkg_config.incdirs + lalmetaio_pkg_config.incdirs + lalsimulation_pkg_config.incdirs + lalburst_pkg_config.incdirs + [numpy_get_include(),"src/xlal"],
			libraries = lal_pkg_config.libs + lalsimulation_pkg_config.libs + lalburst_pkg_config.libs,
			library_dirs = lal_pkg_config.libdirs + lalsimulation_pkg_config.libdirs + lalburst_pkg_config.libdirs,
			runtime_library_dirs = lal_pkg_config.libdirs + lalsimulation_pkg_config.libdirs + lalburst_pkg_config.libdirs,
			extra_compile_args = lal_pkg_config.extra_cflags
		),
		Extension(
			"pylal.xlal.tools",
			["src/xlal/tools.c", "src/xlal/misc.c"],
			include_dirs = lal_pkg_config.incdirs + lalmetaio_pkg_config.incdirs + lalinspiral_pkg_config.incdirs + [numpy_get_include(), "src/xlal"],
			libraries = lal_pkg_config.libs + lalinspiral_pkg_config.libs,
			library_dirs = lal_pkg_config.libdirs + lalinspiral_pkg_config.libdirs,
			runtime_library_dirs = lal_pkg_config.libdirs + lalinspiral_pkg_config.libdirs,
			extra_compile_args = lal_pkg_config.extra_cflags
		),
		Extension(
			"pylal._spawaveform",
			["src/_spawaveform.c"],
			include_dirs = lal_pkg_config.incdirs + lalinspiral_pkg_config.incdirs + [numpy_get_include()],
			libraries = lal_pkg_config.libs + lalinspiral_pkg_config.libs,
			library_dirs = lal_pkg_config.libdirs + lalinspiral_pkg_config.libdirs,
			runtime_library_dirs = lal_pkg_config.libdirs + lalinspiral_pkg_config.libdirs,
			extra_compile_args = lal_pkg_config.extra_cflags
		),
		Extension(
			"pylal._bayespputils",
			["src/bayespputils.c","src/burnin.c"],
			include_dirs = [numpy_get_include(),'src/']
		),
		Extension(
			"pylal._stats",
			["src/_stats.c"],
			include_dirs = [numpy_get_include()],
			extra_compile_args = ["-std=c99"]
		),
		Extension(
			"pylal.cbc_network_efficiency",
			["src/cbc_network_efficiency.c"],
			include_dirs = [numpy_get_include()] + gsl_pkg_config.incdirs + lalsimulation_pkg_config.incdirs + lal_pkg_config.incdirs,
			library_dirs = gsl_pkg_config.libdirs + lalsimulation_pkg_config.libdirs + lal_pkg_config.libdirs,
			runtime_library_dirs = gsl_pkg_config.libdirs + lalsimulation_pkg_config.libdirs + lal_pkg_config.libdirs,
			libraries = gsl_pkg_config.libs + lalsimulation_pkg_config.libs + lal_pkg_config.libs,
			extra_compile_args = ["-std=c99", "-ffast-math", "-O3", "-Wno-unknown-pragmas"] + gsl_pkg_config.extra_cflags + lalsimulation_pkg_config.extra_cflags + lal_pkg_config.extra_cflags,  # , "-fopenmp"
			#extra_link_args=['-fopenmp']
		),
		Extension(
			"pylal.inspiral_metric",
			["src/inspiral_metric.c", "src/xlal/misc.c"],
			include_dirs = lal_pkg_config.incdirs + lalinspiral_pkg_config.incdirs + ["src/xlal/", "src/xlal/datatypes/"],
			libraries = lal_pkg_config.libs + lalinspiral_pkg_config.libs,
			library_dirs = lal_pkg_config.libdirs + lalinspiral_pkg_config.libdirs,
			runtime_library_dirs = lal_pkg_config.libdirs + lalinspiral_pkg_config.libdirs,
			extra_compile_args = lal_pkg_config.extra_cflags
		),
	],
	scripts = [
		os.path.join("bin", "analyseQscan.py"),
		os.path.join("bin", "distrib_fu_qscan_results.py"),
		os.path.join("bin", "submit_remote_scan.py"),
		os.path.join("bin", "fup_triggers.py"),
		os.path.join("bin", "galaxies_in_polygon"),
		os.path.join("bin", "lal_query_cache"),
		os.path.join("bin", "makeCheckList.py"),
		os.path.join("bin", "makeCheckListWiki.py"),
		os.path.join("bin", "followupQueryDQ.py"),
		os.path.join("bin", "followupQueryVeto.py"),
		os.path.join("bin", "followupRatioTest.py"),
		os.path.join("bin", "followupGetDataAsAscii.py"),
		os.path.join("bin", "followupGenerateDQBackground.py"),
		os.path.join("bin", "followupCustomFOM.py"),
		os.path.join("bin", "followupPDSurface.py"),		
		os.path.join("bin", "plotbank"),
		os.path.join("bin", "plotchannel"),
		os.path.join("bin", "plotcohsnr"),
		os.path.join("bin", "plotchiatimeseries"),
		os.path.join("bin", "plotdetresponse"),
		os.path.join("bin", "plotextrapolation"),
		os.path.join("bin", "plotgrbl"),
		os.path.join("bin", "plotlalseries"),
		os.path.join("bin", "plotnumgalaxies"),
		os.path.join("bin", "lalapps_compute_posterior"),
		os.path.join("bin", "plotulvsmass"),
		os.path.join("bin", "plotifar"),
		os.path.join("bin", "plotinspiral"),
		os.path.join("bin", "plotinspinj"),
		os.path.join("bin", "plotinspmissed"),
		os.path.join("bin", "plotnumtemplates"),
		os.path.join("bin", "plotinspiralrange"),
		os.path.join("bin", "plotcoincseglength"),
		os.path.join("bin", "plotsegments"),
		os.path.join("bin", "plotthinca"),
		os.path.join("bin", "plot_medianmax_sngl_inspiral"),
		os.path.join("bin", "plot_num_sngl_inspiral"),
		os.path.join("bin", "plot_tmpltbank_range"),
		os.path.join("bin", "pylal_cache_to_mvsc"),
		os.path.join("bin", "pylal_mvsc_player"),
		os.path.join("bin", "pylal_cbc_dq_page"),
		os.path.join("bin", "mvsc_plot"),
		os.path.join("bin", "mvsc_plot_ringdown"),
		os.path.join("bin", "mvsc_ROC_and_histograms"),
		os.path.join("bin", "mvsc_plot_cuts"),
		os.path.join("bin", "pylal_combine_posteriors"),
		os.path.join("bin", "pylal_followup_missed"),
		os.path.join("bin", "pylal_exttrig_summary"),
		os.path.join("bin", "pylal_grblikelihood"),
		os.path.join("bin", "pylal_grbUL"),
		os.path.join("bin", "pylal_grbtimeslide_stats"),
		os.path.join("bin", "pylal_exttrig_llmonitor"),
		os.path.join("bin", "pylal_exttrig_llsummary"),
		os.path.join("bin", "pylal_exttrig_llbox"),
		os.path.join("bin", "pylal_exttrig_llpage"),
		os.path.join("bin", "pylal_make_imr_summary_page"),
		os.path.join("bin", "pylal_recolor"),
		os.path.join("bin", "pylal_relic"),
		os.path.join("bin", "pylal_version"),
		os.path.join("bin", "plotethinca"),
		os.path.join("bin", "ploteffdistcut"),
		os.path.join("bin", "plotefficiency"),
		os.path.join("bin", "plotsnrchi"),
		os.path.join("bin", "frame_check"),
		os.path.join("bin", "IFOstatus_check"),
		os.path.join("bin", "plotsnrchisq_pipe"),
		os.path.join("bin", "plotmcmc.py"),
		os.path.join("bin", "plotspinmcmc.py"),
		os.path.join("bin", "plotinsppop"),
		os.path.join("bin", "lalapps_binj_pic"),
		os.path.join("bin", "lalapps_burca_tailor"),
		os.path.join("bin", "lalapps_cbc_plotroc"),
		os.path.join("bin", "lalapps_cbc_plotroc_ring"),
		os.path.join("bin", "lalapps_cbc_plotsummary"),
		os.path.join("bin", "lalapps_cbc_plot_likelihood_arrays"),
		os.path.join("bin", "lalapps_excesspowerfinal"),
		os.path.join("bin", "lalapps_farburst"),
		os.path.join("bin", "lalapps_followup_pipe"),
		os.path.join("bin", "lalapps_followup_page"),
		os.path.join("bin", "WOD_Bologna.py"),
		os.path.join("bin", "lalapps_likeliness"),
		os.path.join("bin", "lalapps_newcorse"),
		os.path.join("bin", "lalapps_cbc_volmc"),
		os.path.join("bin", "lalapps_cbc_svim"),
		os.path.join("bin", "lalapps_cbc_sink"),
		os.path.join("bin", "lalapps_plot_tisi"),
		os.path.join("bin", "lalapps_power_calc_likelihood"),
		os.path.join("bin", "lalapps_power_plot_binj"),
		os.path.join("bin", "lalapps_power_plot_burca"),
		os.path.join("bin", "lalapps_power_plot_burca2"),
		os.path.join("bin", "lalapps_power_plot_burst"),
		os.path.join("bin", "lalapps_power_plot_burstrate"),
		os.path.join("bin", "lalapps_remote_cache"),
		os.path.join("bin", "lalapps_run_sqlite"),
		os.path.join("bin", "lalapps_stringfinal"),
		os.path.join("bin", "lalapps_string_calc_likelihood"),
		os.path.join("bin", "lalapps_string_contour_plotter"),
		os.path.join("bin", "lalapps_string_meas_likelihood"),
		os.path.join("bin", "lalapps_string_plot_binj"),
		os.path.join("bin", "lalapps_string_plot_likelihood"),
		os.path.join("bin", "wscan_background.py"),
		os.path.join("bin", "wscan_bg_setup_log.py"),
		os.path.join("bin", "ligolw_binjfind"),
		os.path.join("bin", "ligolw_summmime"),
		os.path.join("bin", "ligolw_bucluster"),
		os.path.join("bin", "ligolw_bucut"),
		os.path.join("bin", "ligolw_burca"),
		os.path.join("bin", "ligolw_cafe"),
		os.path.join("bin", "ligolw_inspinjfind"),
		os.path.join("bin", "lalapps_cbc_injfind"),
		os.path.join("bin", "ligolw_rinca"),
		os.path.join("bin", "ligolw_segments"),
		os.path.join("bin", "ligolw_sicluster"),
		os.path.join("bin", "ligolw_thinca"),
		os.path.join("bin", "ligolw_sstinca"),
		os.path.join("bin", "ligolw_veto"),
		os.path.join("bin", "ligolw_cbc_hardware_inj_page"),
		os.path.join("bin", "ligolw_fr_to_science"),
		os.path.join("bin", "KW_veto_setup"),
		os.path.join("bin", "KW_veto_calc"),
		os.path.join("bin", "KW_veto_plots"),
		os.path.join("bin", "KW_veto_channelPage"),
		os.path.join("bin", "KW_veto_reportPage"),
		os.path.join("bin", "KW_veto_insert"),
		os.path.join("bin", "pylal_plot_inspiral_skymap"),
		os.path.join("bin", "upper_limit_results"),
		os.path.join("bin", "pylal_expose"),
		os.path.join("bin", "ligolw_cbc_expected_snrs"),
		os.path.join("bin", "ligolw_cbc_align_total_spin"),
		os.path.join("bin", "ligolw_cbc_dbsimplify"),
		os.path.join("bin", "ligolw_cbc_dbaddinj"),
		os.path.join("bin", "ligolw_cbc_printlc"),
		os.path.join("bin", "ligolw_cbc_cluster_coincs"),
		os.path.join("bin", "ligolw_cbc_cfar"),
		os.path.join("bin", "ligolw_cbc_jitter_skyloc"),
		os.path.join("bin", "ligolw_cbc_plotslides"),
		os.path.join("bin", "ligolw_cbc_plotsngl"),
		os.path.join("bin", "ligolw_cbc_plotifar"),
		os.path.join("bin", "ligolw_cbc_plotfm"),
		os.path.join("bin", "ligolw_cbc_plotvt"),
		os.path.join("bin", "lalapps_cbc_plotrates"),
		os.path.join("bin", "pylal_cbc_extended_background"),
		os.path.join("bin", "ligolw_cbc_compute_durations"),
		os.path.join("bin", "ligolw_cbc_repop_coinc"),
		os.path.join("bin", "ligolw_segments_compat"),
		os.path.join("bin", "extractCommand"),
		os.path.join("bin", "make_inspiral_summary_page"),
		os.path.join("bin", "mvsc_update_sql"),
		os.path.join("bin", "mvsc_get_doubles"),
		os.path.join("bin", "mvsc_get_doubles_ringdown"),
		os.path.join("bin", "mvsc_get_doubles_snglstage"),
		os.path.join("bin", "mvsc_dag"),
		os.path.join("bin", "mvsc_dag_ringdown"),
		os.path.join("bin", "mvsc_dag_snglstage"),
		os.path.join("bin", "post_process_pipe"),
		os.path.join("bin", "prepare_sendback.py"),
		os.path.join("bin", "qsub_wscan.sh"),
		os.path.join("bin", "qsub_wscanlite.sh"),
		os.path.join("bin", "search_volume_by_m1_m2"),
		os.path.join("bin", "search_upper_limit_by_m1_m2"),
		os.path.join("bin", "search_upper_limit_by_s1z_s2z"),
		os.path.join("bin", "search_volume_by_s1z_s2z"),
		os.path.join("bin", "search_upper_limit_by_m_chi"),
		os.path.join("bin", "search_volume_by_m_chi"),
		os.path.join("bin", "search_upper_limit_by_M"),
		os.path.join("bin", "imr_compare"),		
		os.path.join("bin", "imr_roc"),				
		os.path.join("bin", "ligolw_omega_to_coinc"),
		os.path.join("bin", "ligolw_cwb_to_coinc"),
		os.path.join("bin", "virgo_qscan_in2p3.py"),
		os.path.join("bin", "wscan_in2p3.sh"),
		os.path.join("bin", "wscanlite_in2p3.sh"),
		os.path.join("bin", "minifollowups"),
		os.path.join("bin", "ligolw_cbc_plotcumhist"),
		os.path.join("bin", "pylal_imr_missed_found"),
		os.path.join("bin", "pylal_imr_livetime"),
		os.path.join("bin", "lalapps_cbc_compute_rs"),
		os.path.join("bin", "lalapps_cbc_print_rs"),
		os.path.join("bin", "ligolw_cbc_printsims"),
		os.path.join("bin", "ligolw_cbc_printmissed"),
		os.path.join("bin", "run_skypoints.py"),
		os.path.join("bin", "make_skypoints_grids.py"),
		os.path.join("bin", "make_skypoints_rankings.py"),
		os.path.join("bin", "plot_skypoints.py"),
		os.path.join("bin", "pylal_cbc_select_hardware_injections"),
		os.path.join("bin", "coh_PTF_efficiency"),
		os.path.join("bin", "coh_PTF_html_summary"),
		os.path.join("bin", "coh_PTF_injfinder"),
		os.path.join("bin", "coh_PTF_injcombiner"),
		os.path.join("bin", "coh_PTF_sbv_plotter"),
		os.path.join("bin", "coh_PTF_trig_cluster"),
		os.path.join("bin", "coh_PTF_trig_combiner"),
		os.path.join("bin", "coh_PTF_inspiral_horizon"),
		os.path.join("bin", "pylal_geom_aligned_bank"),
		os.path.join("bin", "pylal_aligned_bank_cat"),
		os.path.join("bin", "pylal_geom_aligned_2dstack"),
		os.path.join("bin", "pylal_aligned_stoch_bank"),
		os.path.join("bin", "pylal_add_missing_sngl_inspiral_columns"),
		os.path.join("bin", "ring_post"),
		os.path.join("bin", "ligolw_rinca_to_coinc"),
		os.path.join("bin", "cbcBayesPostProc.py"),
                os.path.join("bin", "cbcBayesPlotSpinDisk.py"),
                os.path.join("bin", "cbcBayesBurstPostProc.py"),
		os.path.join("bin", "cbcBayesCompPos.py"),
		os.path.join("bin", "cbcBayesDIEvidence.py"),
		os.path.join("bin", "cbcBayesInjProc.py"),
		os.path.join("bin", "cbcBayesThermoInt.py"),
                os.path.join("bin", "cbcBayesBurstPPAnalysis.py"),
		os.path.join("bin", "cbcBayesPPAnalysis.py"),
                os.path.join("bin", "cbcBayesPosToSimBurst.py"),
		os.path.join("bin", "cbcBayesPosToSimInspiral.py"),
		os.path.join("bin", "pulsarBayesPostProc.py"),
		os.path.join("bin", "processing_kdtreeSky.py"),
		os.path.join("bin", "ligo_channel_query"),
		os.path.join("bin", "projectedDetectorTensor"),
		os.path.join("bin", "pylal_exttrig_dataquery"),
		os.path.join("bin", "pylal_exttrig_allquery"),
		os.path.join("bin", "cbcBayesDIEvidence.py"),
		os.path.join("bin", "cbcBayesInjProc.py"),
		os.path.join("bin", "ligo_channel_query"),
		os.path.join("bin", "pylal_exttrig_dataquery"),
		os.path.join("bin", "pylal_exttrig_allquery"),
		os.path.join("bin", "auxmvc_mvsc_results_plots.py"),
		os.path.join("bin", "auxmvc_generate_spr_files.py"),
		os.path.join("bin", "auxmvc_create_mvsc_dag"),
		os.path.join("bin", "auxmvc_create_ann_dag"),
		os.path.join("bin", "auxmvc_ROC_combiner.py"),
		os.path.join("bin", "auxmvc_plot_mvsc_channels_significance.py"),
		os.path.join("bin", "auxmvc_comparison_plots.py"),
		os.path.join("bin", "pylal_imr_search_volume"),
		os.path.join("bin", "pylal_imr_plot_search_volume"),
		os.path.join("bin", "pylal_imr_plot_ifar"),
		os.path.join("bin", "pylal_cbc_rankprod"),
		os.path.join("bin", "pylal_cbc_ulmc"),
		os.path.join("bin", "ligolw_split_inspiral"),
		os.path.join("bin", "ligolw_cbc_stfar"),
		os.path.join("bin", "ligolw_far_compare"),
		os.path.join("bin", "cbcBayesConvergence.py"),
		os.path.join("bin", "pylal_noise_budget"),
		os.path.join("bin", "pylal_summary_page"),
		os.path.join("bin", "pylal_plot_triggers"),
		os.path.join("bin", "ligolw_miinjfind"),
		os.path.join("bin", "pylal_coh_PTF_cluster"),
		os.path.join("bin", "pylal_coh_PTF_ihope"),
		os.path.join("bin", "pylal_coh_PTF_dqv"),
		os.path.join("bin", "pylal_coh_PTF_sbv"),
		os.path.join("bin", "pylal_coh_PTF_plot_results"),
		os.path.join("bin", "pylal_write_coh_PTF_page"),
		os.path.join("bin", "ligolw_dbinjfind"),
		os.path.join("bin", "lalapps_cbc_plothzdist"),
		os.path.join("bin", "pylal_plot_segments"),
		os.path.join("bin", "ligolw_fix_ids")
		],
	data_files = [ ("etc", [
		os.path.join("etc", "pylal-user-env.sh"),
		os.path.join("etc", "pylal-user-env.csh"),
		] ) ],
    classifiers = [
        'Development Status :: 7 - Inactive',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 2 :: Only',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Physics'
    ]
)

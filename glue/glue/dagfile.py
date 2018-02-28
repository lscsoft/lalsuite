# Copyright (C) 2011--2015  Kipp Cannon
# Copyright (C) 2004--2006  Brian Moe
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
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

"""
Machinery for reading, editing, and writing Condor DAG files.

When running DAGs on Condor compute clusters, very often one will wish to
re-run a portion of a DAG.  This can be done by marking all jobs except the
ones to be re-run as "DONE".  Unfortunately the Condor software suite lacks
an I/O library for reading and writing Condor DAG files, so there is no
easy way to edit DAG files except by playing games sed, awk, or once-off
Python or Perl scripts.  That's where this module comes in.  This module
will read a DAG file into an in-ram representation that is easily edited,
and allow the file to be written to disk again.

Example:

>>> from glue import dagfile
>>> dag = dagfile.DAG.parse(open("pipeline.dag"))
>>> dag.write(open("pipeline.dag", "w"))

Although it is possible to machine-generate an original DAG file using this
module and write it to disk, this module does not provide the tools
required to do any of the other tasks associated with pipeline
construction.  For example there is no facility here to generate or manage
submit files, data files, or any other files that are associated with a
full pipeline.  Only the DAG file itself is considered here.  For general
pipeline construction see the pipeline module.  The focus of this module is
on editing existing DAG files.

Developers should also consider doing any new pipeline development using
DAX files as the fundamental workflow description, instead of DAGs.  See
http://pegasus.isi.edu for more information.

A DAG file is loaded using the .parse() class method of the DAG class.
This parses the file-like object passed to it and returns an instance of
the DAG class representing the file's contents.  Once loaded, the nodes in
the DAG can all be found in the .nodes dictionary, whose keys are the node
names and whose values are the corresponding node objects.  Among each node
object's attributes are sets .children and .parents containing references
to the child and parent nodes (not their names) for each node.  Note that
every node must appear listed as a parent of each of its children, and vice
versa.  The other attributes of a DAG instance contain information about
the DAG, for example the CONFIG file or the DOT file, and so on.  All of
the data for each node in the DAG, for example the node's VARS value, its
initial working directory, and so on, can be found in the attributes of the
nodes themselves.  A DAG is written to a file using the .write() method of
the DAG object.
"""


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#


import re


__all__ = ["DAG", "JOB", "DATA", "SPLICE", "SUBDAG_EXTERNAL"]


#
# =============================================================================
#
#                               Progress Wrapper
#
# =============================================================================
#


class progress_wrapper(object):
	"""
	Progress report wrapper.  For internal use only.
	"""
	def __init__(self, f, callback):
		self.n = 0
		self.f = f
		self.callback = callback

	def __iadd__(self, dn):
		self.n += dn
		if self.callback is not None and not self.n % 7411:
			self.callback(self.f, self.n, False)
		return self

	def __del__(self):
		if self.callback is not None:
			self.callback(self.f, self.n, True)


class nofile(object):
	"""
	Object providing a no-op .write() method to fake a file.  For
	internal use only.
	"""
	def write(self, *args):
		pass


#
# =============================================================================
#
#                      The Contents of a Condor DAG File
#
# =============================================================================
#


class JOB(object):
	"""
	Representation of a JOB node in a Condor DAG.  JOB objects have the
	following attributes corresponding to information in the DAG file:

	.name
		The name of the node in the DAG.

	.filename
		The name of the submit file for the JOB.

	.directory
		The initial working directory for the JOB.  Set to None to
		omit from DAG (job's working directory will be chosen by
		Condor).

	.done
		Boolean indicating if the JOB is done or not.  See
		DAG.load_rescue() for more information.

	.noop
		Boolean indicating if the JOB is a no-op or not.

	.vars
		A dictionary of the name-->value pairs in the VARS line for
		the JOB.  Leave empty to omit VARS from DAG.

	.retry
		The number of retries for the job.  Set to None to omit
		from DAG.

	.retry_unless_exit_value
		The value of the UNLESS-EXIT suffix of the RETRY line.
		Set to None to omit from DAG.

	.priority
	.category
		The PRIORITY value and CATEGORY name for the node in the
		DAG.  Set to None to omit from the DAG.

	.parents
	.children
		Sets of the parent and child nodes of JOB.  The sets
		contain references to the node objects, not their names.

	.prescript
	.prescriptargs
	.postscript
	.postscriptargs
		The names and lists of arguments of the PRE and POST
		scripts.  Set to None to omit from DAG.

	.abort_dag_on_abortexitvalue
	.abort_dag_on_dagreturnvalue
		The ABORT-DAG-ON abort exit value and DAG return value for
		the JOB.  Set to None to omit from DAG.

	For more information about the function of these parameters, refer
	to the Condor documentation.
	"""
	keyword = "JOB"

	def __init__(self, name, filename, directory = None, done = False, noop = False):
		# information from the JOB line in the DAG file
		self.name = name
		self.filename = filename
		self.directory = directory
		self.done = done
		self.noop = noop

		# the VARS line in the DAG file.  orderless name, value
		# pairs
		self.vars = {}

		# the RETRY line in the DAG file
		self.retry = None
		self.retry_unless_exit_value = None

		# the PRIORITY and CATEGORY lines in the DAG file
		self.priority = None
		self.category = None

		# the parents and children of this node.  the sets contain
		# references to the parent and child objects, not their
		# names
		self.parents = set()
		self.children = set()

		# the names and arguments of the PRE and POST scripts, if
		# any
		self.prescript = None
		self.prescriptargs = None
		self.postscript = None
		self.postscriptargs = None

		# the ABORT-DAG-ON abort exit value and dag return value
		# for this job if they are set, or None if not
		self.abort_dag_on_abortexitvalue = None
		self.abort_dag_on_dagreturnvalue = None

	def write(self, f, progress = None):
		"""
		Write the lines describing this node to the file-like
		object f.  The object must provide a .write() method.

		If progress is not None, it will be incremented by 1 for
		every line written.
		"""
		# JOB ...
		f.write("%s %s %s" % (self.keyword, self.name, self.filename))
		if self.directory is not None:
			f.write(" DIR \"%s\"" % self.directory)
		if self.noop:
			f.write(" NOOP")
		if self.done:
			f.write(" DONE")
		f.write("\n")
		if progress is not None:
			progress += 1

		# PRIORITY ...
		if self.priority:
			f.write("PRIORITY %s %d\n" % (self.name, self.priority))
			if progress is not None:
				progress += 1

		# CATEGORY ...
		if self.category is not None:
			f.write("CATEGORY %s %s\n" % (self.name, self.category))
			if progress is not None:
				progress += 1

		# RETRY ...
		if self.retry:
			f.write("RETRY %s %d" % (self.name, self.retry))
			if self.retry_unless_exit_value is not None:
				f.write(" UNLESS-EXIT %d" % self.retry_unless_exit_value)
			f.write("\n")
			if progress is not None:
				progress += 1

		# VARS ...
		if self.vars:
			f.write("VARS %s" % self.name)
			for name, value in sorted(self.vars.items()):
				# apply escape rules to the value
				f.write(" %s=\"%s\"" % (name, value.replace("\\", "\\\\").replace("\"", "\\\"")))
			f.write("\n")
			if progress is not None:
				progress += 1

		# SCRIPT PRE ...
		if self.prescript is not None:
			f.write("SCRIPT PRE %s %s" % (self.name, self.prescript))
			if self.prescriptargs:
				f.write(" %s" % " ".join(self.prescriptargs))
			f.write("\n")
			if progress is not None:
				progress += 1

		# SCRIPT POST ...
		if self.postscript is not None:
			f.write("SCRIPT POST %s %s" % (self.name, self.postscript))
			if self.postscriptargs:
				f.write(" %s" % " ".join(self.postscriptargs))
			f.write("\n")
			if progress is not None:
				progress += 1

		# ABORT-DAG-ON ...
		if self.abort_dag_on_abortexitvalue is not None:
			f.write("ABORT-DAG-ON %s %d" % (self.name, self.abort_dag_on_abortexitvalue))
			if self.abort_dag_on_dagreturnvalue is not None:
				f.write(" RETURN %d" % self.abort_dag_on_dagreturnvalue)
			f.write("\n")
			if progress is not None:
				progress += 1

	# state
	@property
	def state(self):
		"""
		Get the state of the node.  One of 'wait', 'idle', 'run',
		'abort', 'stop', 'success', 'fail'.

		NOTE:  this feature is not implemented at this time.
		"""
		raise NotImplemented


class DATA(JOB):
	"""
	Representation of a Stork DATA node in a Condor DAG.
	"""
	keyword = "DATA"


class SUBDAG_EXTERNAL(JOB):
	"""
	Representation of a SUBDAG EXTERNAL node in a Condor DAG.
	"""
	keyword = "SUBDAG EXTERNAL"


class SPLICE(JOB):
	"""
	Representation of a SPLICE node in a Condor DAG.
	"""
	# NOTE:  although this is a subclass of the JOB class, splices
	# don't support most of the things that can be associated with
	# jobs, like VARS and so on, so don't set attributes that shouldn't
	# be set or you'll get a nonsense DAG.  In the future, more error
	# checking might be added to prevent mis-use
	keyword = "SPLICE"


class DAG(object):
	"""
	Representation of the contents of a Condor DAG file.

	BUGS:  the semantics of the "+" special character in category names
	is not understood.  For now, it is an error for a node's category
	to not be found verbatim in a MAXJOBS line.  The "+" character is a
	wildcard-like character used in the assignment of MAXJOBS values to
	job categories in splices;  see the Condor documentation for more
	information.
	"""

	#
	# lines in DAG files
	#

	dotpat = re.compile(r'^DOT\s+(?P<filename>\S+)(\s+(?P<options>.+))?', re.IGNORECASE)
	jobpat = re.compile(r'^JOB\s+(?P<name>\S+)\s+(?P<filename>\S+)(\s+DIR\s+(?P<directory>\S+))?(\s+(?P<noop>NOOP))?(\s+(?P<done>DONE))?', re.IGNORECASE)
	datapat = re.compile(r'^DATA\s+(?P<name>\S+)\s+(?P<filename>\S+)(\s+DIR\s+(?P<directory>\S+))?(\s+(?P<noop>NOOP))?(\s+(?P<done>DONE))?', re.IGNORECASE)
	subdagpat = re.compile(r'^SUBDAG\s+EXTERNAL\s+(?P<name>\S+)\s+(?P<filename>\S+)(\s+DIR\s+(?P<directory>\S+))?(\s+(?P<noop>NOOP))?(\s+(?P<done>DONE))?', re.IGNORECASE)
	splicepat = re.compile(r'^SPLICE\s+(?P<name>\S+)\s+(?P<filename>\S+)(\s+DIR\s+(?P<directory>\S+))?', re.IGNORECASE)
	prioritypat = re.compile(r'^PRIORITY\s+(?P<name>\S+)\s+(?P<value>\S+)', re.IGNORECASE)
	categorypat = re.compile(r'^CATEGORY\s+(?P<name>\S+)\s+(?P<category>\S+)', re.IGNORECASE)
	retrypat = re.compile(r'^RETRY\s+(?P<name>\S+)\s+(?P<retries>\S+)(\s+UNLESS-EXIT\s+(?P<retry_unless_exit_value>\S+))?', re.IGNORECASE)
	varspat = re.compile(r'^VARS\s+(?P<name>\S+)\s+(?P<vars>.+)', re.IGNORECASE)
	varsvaluepat = re.compile(r'(?P<name>\S+)\s*=\s*"(?P<value>.*?)(?<!\\)"', re.IGNORECASE)
	scriptpat = re.compile(r'^SCRIPT\s+(?P<type>(PRE)|(POST))\s(?P<name>\S+)\s+(?P<executable>\S+)(\s+(?P<arguments>.+))?', re.IGNORECASE)
	abortdagonpat = re.compile(r'^ABORT-DAG-ON\s+(?P<name>\S+)\s+(?P<exitvalue>\S+)(\s+RETURN\s+(?P<returnvalue>\S+))?', re.IGNORECASE)
	arcpat = re.compile(r'^PARENT\s+(?P<parents>.+?)\s+CHILD\s+(?P<children>.+)', re.IGNORECASE)
	maxjobspat = re.compile(r'^MAXJOBS\s+(?P<category>\S+)\s+(?P<value>\S+)', re.IGNORECASE)
	configpat = re.compile(r'^CONFIG\s+(?P<filename>\S+)', re.IGNORECASE)
	nodestatuspat = re.compile(r'^NODE_STATUS_FILE\s+(?P<filename>\S+)(\s+(?P<updatetime>\S+))?', re.IGNORECASE)
	jobstatepat = re.compile(r'^JOBSTATE_LOG\s+(?P<filename>\S+)', re.IGNORECASE)

	#
	# lines in rescue DAG files
	#

	donepat = re.compile(r'^DONE\s+(?P<name>\S+)', re.IGNORECASE)

	#
	# methods
	#

	def __init__(self):
		# node name --> JOB object mapping
		self.nodes = {}
		# category name --> integer max jobs value mapping.  all
		# categories are listed, that is it is an error for a JOB
		# in the DAG to claim to be in a category that cannot be
		# found in this dictionary.  categories that don't have a
		# MAXJOBS set for them use None as their max jobs value in
		# this dictionary.
		self.maxjobs = {}
		# filename or None
		self.config = None
		# filename or None
		self.dot = None
		# booleans, defaults match Condor's
		self.dotupdate = False
		self.dotoverwrite = True
		# filename or None
		self.dotinclude = None
		# filename and update time or None
		self.node_status_file = None
		self.node_status_file_updatetime = None
		# filename or None
		self.jobstate_log = None

	def reindex(self):
		"""
		Rebuild the .nodes index.  This is required if the names of
		nodes are changed.
		"""
		# the .nodes object has its contents replaced instead of
		# building a new object so that if external code is holding
		# a reference to it that code sees the new index as well
		nodes = dict((node.name, node) for node in self.nodes.values())
		if len(nodes) != len(self.nodes):
			raise ValueError("node names are not unique")
		self.nodes.clear()
		self.nodes.update(nodes)

	@classmethod
	def parse(cls, f, progress = None):
		"""
		Parse the file-like object f as a Condor DAG file.  Return
		a DAG object.  The file object must be iterable, yielding
		one line of text of the DAG file in each iteration.

		If the progress argument is not None, it should be a
		callable object.  This object will be called periodically
		and passed the f argument, the current line number, and a
		boolean indicating if parsing is complete.  The boolean is
		always False until parsing is complete, then the callable
		will be invoked one last time with the final line count and
		the boolean set to True.

		Example:

		>>> def progress(f, n, done):
		...	print "reading %s: %d lines\\r" % (f.name, n),
		...	if done:
		...		print
		...
		>>> dag = DAG.parse(open("pipeline.dag"), progress = progress)
		"""
		progress = progress_wrapper(f, progress)
		self = cls()
		arcs = []
		for n, line in enumerate(f, start = 1):
			# progress
			progress += 1
			# skip comments and blank lines
			line = line.strip()
			if not line or line.startswith("#"):
				continue
			# JOB ...
			m = self.jobpat.search(line)
			if m is not None:
				if m.group("name") in self.nodes:
					raise ValueError("line %d: duplicate JOB %s" % (n, m.group("name")))
				self.nodes[m.group("name")] = JOB(m.group("name"), m.group("filename"), directory = m.group("directory") and m.group("directory").strip("\""), done = bool(m.group("done")), noop = bool(m.group("noop")))
				continue
			# DATA ...
			m = self.datapat.search(line)
			if m is not None:
				if m.group("name") in self.nodes:
					raise ValueError("line %d: duplicate DATA %s" % (n, m.group("name")))
				self.nodes[m.group("name")] = DATA(m.group("name"), m.group("filename"), directory = m.group("directory") and m.group("directory").strip("\""), done = bool(m.group("done")), noop = bool(m.group("noop")))
				continue
			# SUBDAG EXTERNAL ...
			m = self.subdagpat.search(line)
			if m is not None:
				if m.group("name") in self.nodes:
					raise ValueError("line %d: duplicate SUBDAG EXTERNAL %s" % (n, m.group("name")))
				self.nodes[m.group("name")] = SUBDAG_EXTERNAL(m.group("name"), m.group("filename"), directory = m.group("directory") and m.group("directory").strip("\""), done = bool(m.group("done")), noop = bool(m.group("noop")))
				continue
			# SPLICE ...
			m = self.splicepat.search(line)
			if m is not None:
				if m.group("name") in self.nodes:
					raise ValueError("line %d: duplicate SPLICE %s" % (n, m.group("name")))
				self.nodes[m.group("name")] = SPLICE(m.group("name"), m.group("filename"), directory = m.group("directory") and m.group("directory").strip("\""))
				continue
			# VARS ...
			m = self.varspat.search(line)
			if m is not None:
				node = self.nodes[m.group("name")]
				# FIXME:  find a way to detect malformed name=value pairs
				for name, value in self.varsvaluepat.findall(m.group("vars")):
					if name in node.vars:
						raise ValueError("line %d: multiple variable %s for %s %s" % (n, name, node.keyword, node.name))
					# apply unescape rules to the value
					node.vars[name] = value.replace("\\\\", "\\").replace("\\\"", "\"")
				continue
			# PARENT ... CHILD ...
			m = self.arcpat.search(line)
			if m is not None:
				parents = m.group("parents").strip().split()
				children = m.group("children").strip().split()
				arcs.extend((parent, child) for parent in parents for child in children)
				continue
			# RETRY ...
			m = self.retrypat.search(line)
			if m is not None:
				node = self.nodes[m.group("name")]
				node.retry = int(m.group("retries"))
				node.retry_unless_exit_value = m.group("retry_unless_exit_value")
				continue
			# SCRIPT ...
			m = self.scriptpat.search(line)
			if m is not None:
				node = self.nodes[m.group("name")]
				if m.group("type").upper() == "PRE":
					if node.prescript is not None:
						raise ValueError("line %d: multiple SCRIPT PRE for %s %s" % (n, node.keyword, node.name))
					node.prescript = m.group("executable")
					if m.group("arguments") is not None:
						node.prescriptargs = m.group("arguments").split()
				elif m.group("type").upper() == "POST":
					if node.postscript is not None:
						raise ValueError("line %d: multiple SCRIPT POST for %s %s" % (n, node.keyword, node.name))
					node.postscript = m.group("executable")
					if m.group("arguments") is not None:
						node.postscriptargs = m.group("arguments").split()
				else:
					assert False	# impossible to get here
				continue
			# PRIORITY ...
			m = self.prioritypat.search(line)
			if m is not None:
				node = self.nodes[m.group("name")]
				if node.priority is not None:
					raise ValueError("line %d: multiple PRIORITY for %s %s" % (n, node.keyword, node.name))
				node.priority = int(m.group("value"))
				continue
			# CATEGORY ...
			m = self.categorypat.search(line)
			if m is not None:
				self.nodes[m.group("name")].category = m.group("category")
				continue
			# ABORT-DAG-ON ...
			m = self.abortdagonpat.search(line)
			if m is not None:
				node = self.nodes[m.group("name")]
				if node.abort_dag_on_abortexitvalue is not None:
					raise ValueError("line %d: multiple ABORT-DAG-ON for %s %s" % (n, node.keyword, node.name))
				node.abort_dag_on_abortexitvalue = int(m.group("exitvalue"))
				if m.group("returnvalue") is not None:
					node.abort_dag_on_dagreturnvalue = int(m.group("returnvalue"))
				continue
			# MAXJOBS ...
			m = self.maxjobspat.search(line)
			if m is not None:
				if m.group("category") in self.maxjobs:
					raise ValueError("line %d: multiple MAXJOBS for category %s" % (n, m.group("category")))
				self.maxjobs[m.group("category")] = int(m.group("value"))
				continue
			# DOT ...
			m = self.dotpat.search(line)
			if m is not None:
				self.dot = m.group("filename")
				options = (m.group("options") or "").split()
				while options:
					option = options.pop(0).upper()
					if option == "UPDATE":
						self.dotupdate = True
					elif option == "DONT-UPDATE":
						self.dotupdate = False
					elif option == "OVERWRITE":
						self.dotoverwrite = True
					elif option == "DONT-OVERWRITE":
						self.dotoverwrite = False
					elif option == "INCLUDE":
						try:
							self.dotinclude = options.pop(0)
						except IndexError:
							raise ValueError("line %d: missing filename for INCLUDE option of DOT" % n)
					else:
						raise ValueError("unrecognized option %s for DOT" % option)
				continue
			# CONFIG ...
			m = self.dotpat.search(line)
			if m is not None:
				if self.config is not None:
					raise ValueError("line %d: multiple CONFIG lines in dag file" % n)
				self.config = m.group("filename")
				continue
			# NODE_STATUS_FILE ...
			m = self.nodestatuspat.search(line)
			if m is not None:
				if self.node_status_file is not None:
					raise ValueError("line %d: multiple NODE_STATUS_FILE lines in dag file" % n)
				self.node_status_file = m.group("filename")
				if m.group(updatetime) is not None:
					self.node_status_file_updatetime = int(m.group("updatetime"))
				continue
			# JOBSTATE_LOG ...
			m = self.jobstatepat.search(line)
			if m is not None:
				# dagman allows more than one of these
				# statements, ignoring all but the first
				if self.jobstate_log is None:
					self.jobstate_log = m.group("filename")
				continue
			# error
			raise ValueError("line %d: invalid line in dag file: %s" % (n, line))
		# progress
		del progress
		# populate parent and child sets
		for parent, child in arcs:
			self.nodes[parent].children.add(self.nodes[child])
			self.nodes[child].parents.add(self.nodes[parent])
		# make sure all categories are known
		for node in self.nodes.values():
			if node.category is not None and node.category not in self.maxjobs:
				self.maxjobs[node.category] = None
		# done
		return self

	@classmethod
	def select_nodes_by_name(cls, dag, nodenames):
		"""
		Construct a new DAG object containing only the nodes whose
		names are in nodenames.

		Example:

		>>> names_to_rerun = set(["triggergen"])
		>>> dag = DAG.select_nodes_by_name(dag, names_to_rerun | dag.get_all_parent_names(names_to_rerun))

		NOTE:  the new DAG object is given references to the node
		(JOB, DATA, etc.) objects in the original DAG, not copies
		of them.  Therefore, editing the node objects, for example
		modifying their parent or child sets, will affect both
		DAGs.  To obtain an independent DAG with its own node
		objects, make a deepcopy of the object that is returned
		(see the copy module in the Python standard library for
		more information).

		Example:

		>>> import copy
		>>> dag = copy.deepcopy(DAG.select_nodes_by_name(dag, names_to_rerun | dag.get_all_parent_names(names_to_rerun)))
		"""
		self = cls()
		self.nodes = dict((name, node) for name, node in dag.nodes.items() if name in nodenames)
		self.maxjobs = dict((category, dag.maxjobs[category]) for category in set(node.category for node in self.nodes.values() if node.category is not None))
		self.config = dag.config
		self.node_status_file = dag.node_status_file
		self.node_status_file_updatetime = dag.node_status_file_updatetime
		self.jobstate_log = dag.jobstate_log
		self.dot = dag.dot
		self.dotupdate = dag.dotupdate
		self.dotoverwrite = dag.dotoverwrite
		self.dotinclude = dag.dotinclude
		return self

	def get_all_parent_names(self, names):
		"""
		Trace the DAG backward from the parents of the nodes whose
		names are given to the head nodes, inclusively, and return
		the set of the names of all nodes visited.

		Example:

		>>> all_parents = dag.get_all_parent_names(["triggergen"])
		"""
		all_parent_names = set()
		nodes_to_scan = set(self.nodes[name] for name in names)
		while nodes_to_scan:
			node = nodes_to_scan.pop()
			nodes_to_scan |= node.parents
			all_parent_names |= set(parent.name for parent in node.parents)
		return all_parent_names

	def get_all_child_names(self, names):
		"""
		Trace the DAG forward from the children of the nodes whose
		names are given to the leaf nodes, inclusively, and return
		the set of the names of all nodes visited.

		Example:

		>>> all_children = dag.get_all_child_names(["triggergen"])
		"""
		all_child_names = set()
		nodes_to_scan = set(self.nodes[name] for name in names)
		while nodes_to_scan:
			node = nodes_to_scan.pop()
			nodes_to_scan |= node.children
			all_child_names |= set(child.name for child in node.children)
		return all_child_names

	def check_edges(self):
		"""
		Check all graph edges for validity.  Checks that each of
		every node's children lists that node as a parent, and vice
		versa, and that all nodes listed in the parent and child
		sets of all nodes are contained in this DAG.  Raises
		ValueError if a problem is found, otherwise returns None.

		Example:

		>>> try:
		...	dag.check_edges()
		... except ValueError as e:
		...	print "edges are broken: %s" % str(e)
		... else:
		...	print "all edges are OK"
		...
		"""
		nodes = set(self.nodes.values())
		for node in nodes:
			for child in node.children:
				if node not in child.parents:
					raise ValueError("node %s is not a parent of its child %s" % (node.name, child.name))
				if child not in nodes:
					raise ValueError("node %s has child %s that is not in DAG" % (node.name, child.name))
			for parent in node.parents:
				if node not in parent.children:
					raise ValueError("node %s is not a child of its parent %s" % (node.name, parent.name))
				if parent not in nodes:
					raise ValueError("node %s has parent %s that is not in DAG" % (node.name, parent.name))

	def load_rescue(self, f, progress = None):
		"""
		Parse the file-like object f as a rescue DAG, using the
		DONE lines therein to set the job states of this DAG.

		In the past, rescue DAGs were full copies of the original
		DAG with the word DONE added to the JOB lines of completed
		jobs.  In version 7.7.2 of Condor, the default format of
		rescue DAGs was changed to a condensed format consisting of
		only the names of completed jobs and the number of retries
		remaining for incomplete jobs.  Currently Condor still
		supports the original rescue DAG format, but the user must
		set the DAGMAN_WRITE_PARTIAL_RESCUE config variable to
		false to obtain one.  This module does not directly support
		the new format, however this method allows a new-style
		rescue DAG to be parsed to set the states of the jobs in a
		DAG.  This, in effect, converts a new-style rescue DAG to
		an old-style rescue DAG, allowing the result to be
		manipulated as before.

		If the progress argument is not None, it should be a
		callable object.  This object will be called periodically
		and passed the f argument, the current line number, and a
		boolean indicating if parsing is complete.  The boolean is
		always False until parsing is complete, then the callable
		will be invoked one last time with the final line count and
		the boolean set to True.
		"""
		# set all jobs to "not done"
		for job in self.nodes.values():
			job.done = False
		# now load rescue DAG, updating done and retries states
		progress = progress_wrapper(f, progress)
		for n, line in enumerate(f):
			# lines are counted from 1, enumerate counts from 0
			n += 1
			# progress
			progress += 1
			# skip comments and blank lines
			line = line.strip()
			if not line or line.startswith("#"):
				continue
			# DONE ...
			m = self.donepat.search(line)
			if m is not None:
				self.nodes[m.group("name")].done = True
				continue
			# RETRY ...
			m = self.retrypat.search(line)
			if m is not None:
				node = self.nodes[m.group("name")]
				node.retry = int(m.group("retries"))
				node.retry_unless_exit_value = m.group("retry_unless_exit_value")
				continue
			# error
			raise ValueError("line %d: invalid line in rescue file: %s" % (n, line))
		# progress
		del progress

	def write(self, f, progress = None, rescue = None):
		"""
		Write the DAG to the file-like object f.  The object must
		provide a .write() method.  In the special case that the
		optional rescue argument is not None (see below) then f can
		be set to None and no DAG file will be written (just the
		rescue DAG will be written).

		If the progress argument is not None, it should be a
		callable object.  This object will be called periodically
		and passed the f argument, the current line number, and a
		boolean indicating if writing is complete.  The boolean is
		always False until writing is complete, then the callable
		will be invoked one last time with the final line count and
		the boolean set to True.

		Example:

		>>> def progress(f, n, done):
		...	print "writing %s: %d lines\\r" % (f.name, n),
		...	if done:
		...		print
		...
		>>> dag.write(open("pipeline.dag", "w"), progress = progress)

		NOTE:  when writing PARENT/CHILD graph edges, this method
		will silently skip any node names that are not in this
		DAG's graph.  This is a convenience to simplify writing
		DAGs constructed by the .select_nodes_by_name() class
		method.  If one wishes to check for broken parent/child
		links before writing the DAG use the .check_edges() method.

		If the optional rescue argument is not None, it must be a
		file-like object providing a .write() method and the DONE
		state of jobs will be written to this file instead of the
		.dag (in the .dag all jobs will be marked not done).

		Example:

		>>> dag.write(open("pipeline.dag", "w"), rescue = open("pipeline.dag.rescue001", "w"))

		NOTE:  it is left as an exercise for the calling code to
		ensure the name chosen for the rescue file is consistent
		with the naming convention assumed by condor_dagman when it
		starts up.
		"""
		# initialize proegress report wrapper
		progress = progress_wrapper(f, progress)

		# if needed, create a dummy object to allow .write() method
		# calls
		if f is None and rescue is not None:
			f = nofile()

		# DOT ...
		if self.dot is not None:
			f.write("DOT %s" % self.dot)
			if self.dotupdate:
				f.write(" UPDATE")
			if not self.dotoverwrite:
				f.write(" DONT-OVERWRITE")
			if self.dotinclude is not None:
				f.write(" INCLUDE %s" % self.dotinclude)
			f.write("\n")
			progress += 1

		# CONFIG ...
		if self.config is not None:
			f.write("CONFIG %s\n" % self.config)
			progress += 1

		# NODE_STATUS_FILE ...
		if self.node_status_file is not None:
			f.write("NODE_STATUS_FILE %s" % self.node_status_file)
			if self.node_status_file_updatetime is not None:
				f.write(" %d" % self.node_status_file_updatetime)
			f.write("\n")
			progress += 1

		# JOBSTATE_LOG ...
		if self.jobstate_log is not None:
			f.write("JOBSTATE_LOG %s\n" % self.jobstate_log)
			progress += 1

		# MAXJOBS ...
		if set(node.category for node in self.nodes.values() if node.category is not None) - set(self.maxjobs):
			raise ValueError("no MAXJOBS statement(s) for node category(ies) %s" % ", ".join(sorted(set(node.category for node in self.nodes.values() if node.category is not None) - set(self.maxjobs))))
		for name, value in sorted(self.maxjobs.items()):
			if value is not None:
				f.write("MAXJOBS %s %d\n" % (name, value))
				progress += 1

		# JOB/DATA/SUBDAG ... (and things that go with them)
		for name, node in sorted(self.nodes.items()):
			if rescue is not None:
				if node.done:
					rescue.write("DONE %s\n" % node.name)
				# save done state, then clear
				done = node.done
				node.done = False
			node.write(f, progress = progress)
			if rescue is not None:
				# restore done state
				node.done = done

		# PARENT ... CHILD ...
		names = set(self.nodes)
		parents_of = {}
		for name, node in self.nodes.items():
			parents_of.setdefault(frozenset(child.name for child in node.children) & names, set()).add(node.name)
		for children, parents in parents_of.items():
			if children:
				f.write("PARENT %s CHILD %s\n" % (" ".join(sorted(parents)), " ".join(sorted(children))))
				progress += 1

		# progress
		del progress

	def dot_source(self, title = "DAG", rename = False, colour = "black", bgcolour = "#a3a3a3", statecolours = {'wait': 'yellow', 'idle': 'yellow', 'run': 'lightblue', 'abort': 'red', 'stop': 'red', 'success': 'green', 'fail': 'red'}):
		"""
		Return a string containing DOT code to generate a
		visualization of the DAG graph.  See
		http://www.graphviz.org for more information.

		title provides a title for the graph.  If rename is True,
		instead of using the names of the nodes for the node names
		in the graph, numbers will be used instead.  The numbers
		are assigned to the nodes in alphabetical order by node
		name.  This might be required if the nodes have names that
		are incompatible with the DOT syntax.

		colour and bgcolour set the outline colour of the graph
		nodes and the background colour for the graph respectively.
		statecolours is a dictionary mapping node state (see the
		.state attribute of the JOB class and its derivatives) to a
		colour.  Set statecolours to None to disable state-based
		colouring of graph nodes.

		Example:

		>>> print dag.dot_source(statecolours = None)

		BUGS:  the JOB class does not implement the ability to
		retrieve the job state at this time, therefore it is always
		necessary to set statecolours to None.  This might change
		in the future.
		"""
		# set up renaming map

		if rename:
			namemap = dict((name, str(n)) for n, name in enumerate(sorted(self.nodes), start = 1))
		else:
			namemap = dict((name, name) for name in self.nodes)

		# generate dot code

		code = 'digraph "%s" {\nnode [color="%s", href="\N"];\ngraph [bgcolor="%s"];\n' % (title, colour, bgcolour)
		for node in self.nodes.values():
			if statecolours is not None:
				code += '"%s"[color="%s"];\n' % (namemap[node.name], statecolours[node.state])
			for child in node.children:
				code += '"%s" -> "%s";\n' % (namemap[node.name], namemap[child.name])
		code += '}\n'

		# done

		return code

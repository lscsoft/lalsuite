# Copyright (C) 2006--2014  Kipp Cannon
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


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#


"""
Library of utility code for LIGO Light Weight XML applications.
"""


import codecs
import gzip
from hashlib import md5
import warnings
import os
import urllib2
import urlparse
import signal
import stat
import sys


from glue import git_version
from .. import ligolw


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


__all__ = ["sort_files_by_size", "local_path_from_url", "load_fileobj", "load_filename", "load_url", "write_fileobj", "write_filename", "write_url"]


#
# =============================================================================
#
#                                 Input/Output
#
# =============================================================================
#


def sort_files_by_size(filenames, verbose = False, reverse = False):
	"""
	Return a list of the filenames sorted in order from smallest file
	to largest file (or largest to smallest if reverse is set to True).
	If a filename in the list is None (used by many glue.ligolw based
	codes to indicate stdin), its size is treated as 0.  The filenames
	may be any sequence, including generator expressions.
	"""
	if verbose:
		if reverse:
			sys.stderr.write("sorting files from largest to smallest ...\n")
		else:
			sys.stderr.write("sorting files from smallest to largest ...\n")
	return sorted(filenames, key = (lambda filename: os.stat(filename)[stat.ST_SIZE] if filename is not None else 0), reverse = reverse)


def local_path_from_url(url):
	"""
	For URLs that point to locations in the local filesystem, extract
	and return the filesystem path of the object to which they point.
	As a special case pass-through, if the URL is None, the return
	value is None.  Raises ValueError if the URL is not None and does
	not point to a local file.

	Example:

	>>> print local_path_from_url(None)
	None
	>>> local_path_from_url("file:///home/me/somefile.xml.gz")
	'/home/me/somefile.xml.gz'
	"""
	if url is None:
		return None
	scheme, host, path = urlparse.urlparse(url)[:3]
	if scheme.lower() not in ("", "file") or host.lower() not in ("", "localhost"):
		raise ValueError("%s is not a local file" % repr(url))
	return path


class RewindableInputFile(object):
	# The GzipFile class in Python's standard library is, in my
	# opinion, somewhat weak.  Instead of relying on the return values
	# from the file object's .read() method, GzipFile checks for EOF
	# using calls to .seek().  Furthermore, it uses .seek() instead of
	# buffering data internally as required.  This makes GzipFile
	# gratuitously unable to work with pipes, urlfile objects, and
	# anything else that does not support seeking (including the
	# MD5File class in this module).  To hack around this, this class
	# provides the buffering needed by GzipFile.  It also does proper
	# EOF checking, and uses the results to emulate the results of
	# GzipFile's .seek() games.
	#
	# By wrapping your file object in this class before passing it to
	# GzipFile, you can use GzipFile to read from non-seekable files.
	#
	# How GzipFile checks for EOF == call .tell() to get current
	# position, seek to end of file with .seek(0, 2), call .tell()
	# again and check if the number has changed from before, if it has
	# then we weren't at EOF so call .seek() with original position and
	# keep going.  ?!

	def __init__(self, fileobj, buffer_size = 1024):
		# the real source of data
		self.fileobj = fileobj
		# where the application thinks it is in the file.  this is
		# used to fake .tell() because file objects that don't
		# support seeking, like stdin, report IOError, and the
		# things returned by urllib don't have a .tell() method at
		# all.
		self.pos = 0
		# how many octets of the internal buffer to return before
		# getting more data
		self.reuse = 0
		# the internal buffer
		self.buf = buffer(" " * buffer_size)
		# flag indicating a .seek()-based EOF test is in progress
		self.gzip_hack_pretend_to_be_at_eof = False
		# avoid attribute look-ups
		self._next = self.fileobj.next
		self._read = self.fileobj.read
		self.close = self.fileobj.close

	def __iter__(self):
		return self

	def next(self):
		if self.gzip_hack_pretend_to_be_at_eof:
			return buffer()
		if self.reuse:
			buf = self.buf[-self.reuse:]
			self.reuse = 0
		else:
			buf = self._next()
			self.buf = (self.buf + buf)[-len(self.buf):]
		self.pos += len(buf)
		return buf

	def read(self, size = None):
		if self.gzip_hack_pretend_to_be_at_eof:
			return buffer()
		if self.reuse:
			if self.reuse < 0:
				buf = self._read(size - self.reuse)
				self.buf = (self.buf + buf)[-len(self.buf):]
				buf = buf[-self.reuse:]
				self.reuse = 0
			# size is None --> condition is False
			elif 0 <= size < self.reuse:
				buf = self.buf[-self.reuse:-self.reuse + size]
				self.reuse -= size
			else:
				buf = self.buf[-self.reuse:]
				self.reuse = 0
				# size is None --> condition is False
				if len(buf) < size:
					buf += self.read(size - len(buf))
		else:
			buf = self._read(size)
			self.buf = (self.buf + buf)[-len(self.buf):]
		self.pos += len(buf)
		return buf

	def seek(self, offset, whence = os.SEEK_SET):
		self.gzip_hack_pretend_to_be_at_eof = False
		if whence == os.SEEK_SET:
			if offset >= 0 and 0 <= self.pos + self.reuse - offset < len(self.buf):
				self.reuse += self.pos - offset
				self.pos = offset
			else:
				raise IOError("seek out of range")
		elif whence == os.SEEK_CUR:
			if self.reuse - len(self.buf) <= offset:
				self.reuse -= offset
				self.pos += offset
			else:
				raise IOError("seek out of range")
		elif whence == os.SEEK_END:
			if offset == 0:
				self.gzip_hack_pretend_to_be_at_eof = True
			else:
				raise IOError("seek out of range")

	def tell(self):
		if self.gzip_hack_pretend_to_be_at_eof:
			# check to see if we are at EOF by seeing if we can
			# read 1 character.  save it in the internal buffer
			# to not loose it.
			c = self._read(1)
			self.buf = (self.buf + c)[-len(self.buf):]
			self.reuse += len(c)
			if c:
				# this will not return the same answer as
				# when GzipFile called it before seeking to
				# EOF
				return self.pos + 1
		return self.pos

	def __enter__(self):
		return self

	def __exit__(self, *args):
		self.close()
		return False


class MD5File(object):
	def __init__(self, fileobj, md5obj = None, closable = True):
		self.fileobj = fileobj
		if md5obj is None:
			self.md5obj = md5()
		else:
			self.md5obj = md5obj
		self.closable = closable
		# avoid attribute look-ups
		self._update = self.md5obj.update
		try:
			self._next = self.fileobj.next
		except AttributeError:
			# replace our .next() method with something that
			# will raise a more meaningful exception if
			# attempted
			self.next = lambda *args, **kwargs: fileobj.next(*args, **kwargs)
		try:
			self._read = self.fileobj.read
		except AttributeError:
			# replace our .read() method with something that
			# will raise a more meaningful exception if
			# attempted
			self.read = lambda *args, **kwargs: fileobj.read(*args, **kwargs)
		try:
			self._write = self.fileobj.write
		except AttributeError:
			# replace our .write() method with something that
			# will raise a more meaningful exception if
			# attempted
			self.write = lambda *args, **kwargs: fileobj.write(*args, **kwargs)
		try:
			self.tell = self.fileobj.tell
		except AttributeError:
			self.tell = lambda *args, **kwargs: fileobj.tell(*args, **kwargs)
		try:
			self.flush = self.fileobj.flush
		except AttributeError:
			self.flush = lambda *args, **kwargs: fileobj.flush(*args, **kwargs)

	def __iter__(self):
		return self

	def next(self):
		buf = self._next()
		self._update(buf)
		return buf

	def read(self, *args):
		buf = self._read(*args)
		self._update(buf)
		return buf

	def write(self, buf):
		self._update(buf)
		return self._write(buf)

	def close(self):
		if self.closable:
			return self.fileobj.close()
		else:
			# at least make sure we're flushed
			self.flush()

	def __enter__(self):
		return self

	def __exit__(self, *args):
		self.close()
		return False


class SignalsTrap(object):
	default_signals = (signal.SIGTERM, signal.SIGTSTP)

	def __init__(self, trap_signals = default_signals):
		self.trap_signals = trap_signals

	def handler(self, signum, frame):
		self.deferred_signals.append(signum)

	def __enter__(self):
		self.oldhandlers = {}
		self.deferred_signals = []
		if self.trap_signals is None:
			return self
		for sig in self.trap_signals:
			self.oldhandlers[sig] = signal.getsignal(sig)
			signal.signal(sig, self.handler)
		return self

	def __exit__(self, *args):
		# restore original handlers
		for sig, oldhandler in self.oldhandlers.iteritems():
			signal.signal(sig, oldhandler)
		# send ourselves the trapped signals in order
		while self.deferred_signals:
			os.kill(os.getpid(), self.deferred_signals.pop(0))
		return False


def load_fileobj(fileobj, gz = None, xmldoc = None, contenthandler = None):
	"""
	Parse the contents of the file object fileobj, and return the
	contents as a LIGO Light Weight document tree.  The file object
	does not need to be seekable.

	If the gz parameter is None (the default) then gzip compressed data
	will be automatically detected and decompressed, otherwise
	decompression can be forced on or off by setting gz to True or
	False respectively.

	If the optional xmldoc argument is provided and not None, the
	parsed XML tree will be appended to that document, otherwise a new
	document will be created.  The return value is a tuple, the first
	element of the tuple is the XML document and the second is a string
	containing the MD5 digest in hex digits of the bytestream that was
	parsed.

	Example:

	>>> from glue.ligolw import ligolw
	>>> import StringIO
	>>> f = StringIO.StringIO('<?xml version="1.0" encoding="utf-8" ?><!DOCTYPE LIGO_LW SYSTEM "http://ldas-sw.ligo.caltech.edu/doc/ligolwAPI/html/ligolw_dtd.txt"><LIGO_LW><Table Name="demo:table"><Column Name="name" Type="lstring"/><Column Name="value" Type="real8"/><Stream Name="demo:table" Type="Local" Delimiter=",">"mass",0.5,"velocity",34</Stream></Table></LIGO_LW>')
	>>> xmldoc, digest = load_fileobj(f, contenthandler = ligolw.LIGOLWContentHandler)
	>>> digest
	'6bdcc4726b892aad913531684024ed8e'

	The contenthandler argument specifies the SAX content handler to
	use when parsing the document.  The contenthandler is a required
	argument.  See the glue.ligolw package documentation for typical
	parsing scenario involving a custom content handler.  See
	glue.ligolw.ligolw.PartialLIGOLWContentHandler and
	glue.ligolw.ligolw.FilteringLIGOLWContentHandler for examples of
	custom content handlers used to load subsets of documents into
	memory.
	"""
	fileobj = MD5File(fileobj)
	md5obj = fileobj.md5obj
	if gz or gz is None:
		fileobj = RewindableInputFile(fileobj)
		magic = fileobj.read(2)
		fileobj.seek(0, os.SEEK_SET)
		if gz or magic == '\037\213':
			fileobj = gzip.GzipFile(mode = "rb", fileobj = fileobj)
	if xmldoc is None:
		xmldoc = ligolw.Document()
	ligolw.make_parser(contenthandler(xmldoc)).parse(fileobj)
	return xmldoc, md5obj.hexdigest()


def load_filename(filename, verbose = False, **kwargs):
	"""
	Parse the contents of the file identified by filename, and return
	the contents as a LIGO Light Weight document tree.  stdin is parsed
	if filename is None.  Helpful verbosity messages are printed to
	stderr if verbose is True.  All other keyword arguments are passed
	to load_fileobj(), see that function for more information.  In
	particular note that a content handler must be specified.

	Example:

	>>> from glue.ligolw import ligolw
	>>> xmldoc = load_filename("demo.xml", contenthandler = ligolw.LIGOLWContentHandler, verbose = True)
	"""
	if verbose:
		sys.stderr.write("reading %s ...\n" % (("'%s'" % filename) if filename is not None else "stdin"))
	if filename is not None:
		fileobj = open(filename, "rb")
	else:
		fileobj = sys.stdin
	xmldoc, hexdigest = load_fileobj(fileobj, **kwargs)
	if verbose:
		sys.stderr.write("md5sum: %s  %s\n" % (hexdigest, (filename if filename is not None else "")))
	return xmldoc


def load_url(url, verbose = False, **kwargs):
	"""
	Parse the contents of file at the given URL and return the contents
	as a LIGO Light Weight document tree.  Any source from which
	Python's urllib2 library can read data is acceptable.  stdin is
	parsed if url is None.  Helpful verbosity messages are printed to
	stderr if verbose is True.  All other keyword arguments are passed
	to load_fileobj(), see that function for more information.  In
	particular note that a content handler must be specified.

	Example:

	>>> from os import getcwd
	>>> from glue.ligolw import ligolw
	>>> xmldoc = load_url("file://localhost/%s/demo.xml" % getcwd(), contenthandler = ligolw.LIGOLWContentHandler, verbose = True)
	"""
	if verbose:
		sys.stderr.write("reading %s ...\n" % (("'%s'" % url) if url is not None else "stdin"))
	if url is not None:
		scheme, host, path = urlparse.urlparse(url)[:3]
		if scheme.lower() in ("", "file") and host.lower() in ("", "localhost"):
			fileobj = open(path)
		else:
			fileobj = urllib2.urlopen(url)
	else:
		fileobj = sys.stdin
	xmldoc, hexdigest = load_fileobj(fileobj, **kwargs)
	if verbose:
		sys.stderr.write("md5sum: %s  %s\n" % (hexdigest, (url if url is not None else "")))
	return xmldoc


def write_fileobj(xmldoc, fileobj, gz = False, compresslevel = 3, **kwargs):
	"""
	Writes the LIGO Light Weight document tree rooted at xmldoc to the
	given file object.  Internally, the .write() method of the xmldoc
	object is invoked and any additional keyword arguments are passed
	to that method.  The file object need not be seekable.  The output
	data is gzip compressed on the fly if gz is True, and in that case
	the compresslevel parameter sets the gzip compression level (the
	default is 3).  The return value is a string containing the hex
	digits of the MD5 digest of the output bytestream.

	Example:

	>>> import sys
	>>> from glue.ligolw import ligolw
	>>> xmldoc = load_filename("demo.xml", contenthandler = ligolw.LIGOLWContentHandler)
	>>> digest = write_fileobj(xmldoc, sys.stdout)	# doctest: +NORMALIZE_WHITESPACE
	<?xml version='1.0' encoding='utf-8'?>
	<!DOCTYPE LIGO_LW SYSTEM "http://ldas-sw.ligo.caltech.edu/doc/ligolwAPI/html/ligolw_dtd.txt">
	<LIGO_LW>
		<Table Name="demo:table">
			<Column Type="lstring" Name="name"/>
			<Column Type="real8" Name="value"/>
			<Stream Delimiter="," Type="Local" Name="demo:table">
	"mass",0.5,"velocity",34
			</Stream>
		</Table>
	</LIGO_LW>
	>>> digest
	'37044d979a79409b3d782da126636f53'
	"""
	with MD5File(fileobj, closable = False) as fileobj:
		md5obj = fileobj.md5obj
		with fileobj if not gz else gzip.GzipFile(mode = "wb", fileobj = fileobj, compresslevel = compresslevel) as fileobj:
			with codecs.getwriter("utf_8")(fileobj) as fileobj:
				xmldoc.write(fileobj, **kwargs)
		return md5obj.hexdigest()


class tildefile(object):
	def __init__(self, filename):
		if not filename:
			raise ValueError(filename)
		self.filename = filename

	def __enter__(self):
		try:
			self.tildefilename = self.filename + "~"
			self.fobj = open(self.tildefilename, "w")
		except IOError:
			self.tildefilename = None
			self.fobj = open(self.filename, "w")
		return self.fobj

	def __exit__(self, exc_type, exc_val, exc_tb):
		self.fobj.close()
		del self.fobj

		#
		# only rename the "~" version to the final destination if
		# no exception has occurred.
		#

		if exc_type is None and self.tildefilename is not None:
			os.rename(self.tildefilename, self.filename)

		return False


def write_filename(xmldoc, filename, verbose = False, gz = False, with_mv = True, trap_signals = SignalsTrap.default_signals, **kwargs):
	"""
	Writes the LIGO Light Weight document tree rooted at xmldoc to the
	file name filename.  If filename is None the file is written to
	stdout, otherwise it is written to the named file.  Friendly
	verbosity messages are printed while writing the file if verbose is
	True.  The output data is gzip compressed on the fly if gz is True.
	If with_mv is True and filename is not None the filename has a "~"
	appended to it and the file is written to that name then moved to
	the requested name once the write has completed successfully.

	Internally, write_fileobj() is used to perform the write.  All
	additional keyword arguments are passed to write_fileobj().

	This function traps the signals in the trap_signals iterable during
	the write process (see SignalsTrap for the default signals), and it
	does this by temporarily installing its own signal handlers in
	place of the current handlers.  This is done to prevent Condor
	eviction during the write process.  When the file write is
	concluded the original signal handlers are restored.  Then, if
	signals were trapped during the write process, the signals are then
	resent to the current process in the order in which they were
	received.  The signal.signal() system call cannot be invoked from
	threads, and trap_signals must be set to None or an empty sequence
	if this function is used from a thread.

	Example:

	>>> write_filename(xmldoc, "demo.xml")	# doctest: +SKIP
	>>> write_filename(xmldoc, "demo.xml.gz", gz = True)	# doctest: +SKIP
	"""
	if verbose:
		sys.stderr.write("writing %s ...\n" % (("'%s'" % filename) if filename is not None else "stdout"))
	with SignalsTrap(trap_signals):
		if filename is None:
			hexdigest = write_fileobj(xmldoc, sys.stdout, gz = gz, **kwargs)
		else:
			if not gz and filename.endswith(".gz"):
				warnings.warn("filename '%s' ends in '.gz' but file is not being gzip-compressed" % filename, UserWarning)
			with (open if not with_mv else tildefile)(filename) as fileobj:
				hexdigest = write_fileobj(xmldoc, fileobj, gz = gz, **kwargs)
	if verbose:
		sys.stderr.write("md5sum: %s  %s\n" % (hexdigest, (filename if filename is not None else "")))


def write_url(xmldoc, url, **kwargs):
	"""
	Writes the LIGO Light Weight document tree rooted at xmldoc to the
	URL name url.

	NOTE:  only URLs that point to local files can be written to at
	this time.  Internally, write_filename() is used to perform the
	write.  All additional keyword arguments are passed to that
	function.  The implementation might change in the future,
	especially if support for other types of URLs is ever added.

	Example:

	>>> write_url(xmldoc, "file:///data.xml")	# doctest: +SKIP
	>>> write_url(xmldoc, "file:///data.xml.gz", gz = True)	# doctest: +SKIP
	"""
	return write_filename(xmldoc, local_path_from_url(url), **kwargs)

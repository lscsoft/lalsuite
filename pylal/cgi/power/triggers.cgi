#!/usr/bin/python

import cgi
import cgitb ; cgitb.enable()
import os
import popen2

import eventdisplay
import webplot


#
# Hack PlotDescription class from webplot.py to describe the output document
#

class Description(webplot.PlotDescription):
	def __init__(self):
		webplot.PlotDescription.__init__(self)
		self.set_format("xml")
		return self

	def parse_form(self):
		webplot.PlotDescription.parse_form(self)
		self.set_format("xml")
		return self


#
# How to run lalapps_lladd
#

class LLAddCommand(object):
	def __init__(self, urls, output = None):
		self._exec = "/home/kipp/local/bin/ligolw_add"
		self.urls = urls
		self.output = output

	def __str__(self):
		s = self._exec
		for url in self.urls:
			s += " \"" + url + "\""
		if self.output:
			s += " --output=\"" + self.output + "\""
		return s

def runlladd(command):
	if type(command) != LLAddCommand:
		raise ValueError, "invalid argument to runlladd(command): command must type LLAddCommand"
	child = popen2.Popen3(str(command), True)
	for line in child.childerr:
		pass
	result = reduce(str.__add__, child.fromchild, "")
	status = child.wait()
	if not os.WIFEXITED(status) or os.WEXITSTATUS(status):
		raise Exception, "failure running \"" + str(command) + "\""
	return result


#
# Generate and send trigger file.
#

description = Description().parse_form()

command = LLAddCommand(webplot.CacheURLs(eventdisplay.cache[description.instrument], description.segment), output = description.filename)

runlladd(command)

webplot.SendImage(description)

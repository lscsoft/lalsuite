# Copyright (C) 2009 Chad Hanna
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

from glue import markup
from glue.markup import oneliner as e
from glue import git_version

import subprocess
import os, sys, time, socket, glob, math
import shutil,urllib

__author__ = "Chad Hanna <channa@caltech.edu>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


###############################################################################
##### UTILITY FUNCTIONS #######################################################
###############################################################################

# THIS SHOULD MOST CERTAINLY NEVER BE USED :)
def web_path_to_url(path):
	host = socket.getfqdn()
	pl = path.rstrip('/').split('/')

	#FIXME add more hosts as you need them
	if 'ligo.caltech.edu' in host: return "https://ldas-jobs.ligo.caltech.edu/~" +pl[pl.index('public_html')-1] + '/' + '/'.join(pl[pl.index('public_html')+1:])
	if 'ligo-la.caltech.edu' in host: return "https://ldas-jobs.ligo-la.caltech.edu/~" +pl[pl.index('public_html')-1] + '/' + '/'.join(pl[pl.index('public_html')+1:])
	if 'ligo-wa.caltech.edu' in host: return "https://ldas-jobs.ligo-wa.caltech.edu/~" +pl[pl.index('public_html')-1] + '/' + '/'.join(pl[pl.index('public_html')+1:])
	if 'phys.uwm.edu' in host: return "https://ldas-jobs.phys.uwm.edu/~" + pl[pl.index('public_html')-1] + '/' + '/'.join(pl[pl.index('public_html')+1:])
	if 'phy.syr.edu' in host: return "https://sugar-jobs.phy.syr.edu/~" + pl[pl.index('public_html')-1] + '/' + '/'.join(pl[pl.index('public_html')+1:])
	if 'aei.uni-hannover.de' in host: return "https://atlas.atlas.aei.uni-hannover.de/~" + pl[pl.index('WWW')-1] + '/' + '/'.join(pl[pl.index('WWW')+1:])
	sys.stderr.write("WARNING: could not find web server, returning empty string\n")
	return ''


def create_toggle(filename="toggle.js"):
	"""
This function is just an alias to create a javascript for the toggle on/off. 
  
@return: nothing
	"""
	fname = open(filename, "w")
	fname.write("""
function get_url_vars() {
    var st = window.location.href.split('?'),
        obj = {}, //an object to store properties and values in
        eq,
        i;
    if (st[1]) { //if a ( ? ) was found in the split, use the second part after the ?
        st = unescape(st[1]).split('&'); //split st into array of strings containing url variables=values
        for (i = 0; i < st.length; i++) {
            eq = st[i].split('='); //get values from both sides of ( = ) sign
            obj[eq[0]] = eq[1]; //insert properties and values into object
        }
        return obj;
    }
    return false;
}

$(document).ready(function() {
$(".fancybox").fancybox();
});

window.onload = function () {
        var vars = get_url_vars(), prop;
        for (url in vars) {
                loadURL(url)
        }

};

function toggle2(showHideDiv, switchTextDiv) {
	var ele = document.getElementById(showHideDiv);
	var text = document.getElementById(switchTextDiv);
	if(ele.style.display == "block") {
    		ele.style.display = "none";
  	}
	else {
		ele.style.display = "block";
	}
}

function afterLoadFrame() {
	$('#iframecontent a[rel="external"]').attr('target','_blank');
	$('#iframecontent input').hide();
	$('#iframecontent p:first').hide(); 
	}

function loadFrame(sourceURL) {
        $("#iframecontent").load(sourceURL,{},afterLoadFrame);
        str = window.location.href.split('?')
        window.location.href = str[0] + "?" + sourceURL
        }

function loadURL(sourceURL) {
        $("#iframecontent").load(sourceURL,{},afterLoadFrame);
        }

function toggleAllOpen() {
	var tags = document.getElementsByTagName('div');
	for (t in tags)
		tags[t].style.display = "block";
	}
	""")
	fname.close()
	return filename

def script_dict(fname):
	script = {}
	tog = os.path.split(create_toggle(fname))[1]
	script[tog] = 'javascript'
	script['https://ldas-jobs.ligo-wa.caltech.edu/~detchar/html/fancybox/source/jquery.fancybox.pack.js?v=2.1.5'] = 'javascript'
	script['https://ldas-jobs.ligo-wa.caltech.edu/~detchar/html/jquery-1.10.2.min.js'] = 'javascript'
	return (script, [tog])


def which(prog):
	which = subprocess.Popen(['which',prog], stdout=subprocess.PIPE)
	out = which.stdout.read().strip()
	if not out:
		sys.stderr.write("ERROR: could not find %s in your path, have you built the proper software and source the proper env. scripts?\n" % (prog,prog))
		raise ValueError
		sys.exit(1)
	return out

def user_and_date():
	tmstr = "/".join([str(i) for i in time.gmtime()[0:3]])
        tmstr += " " + ":".join([str(i) for i in time.gmtime()[3:5]])
	return "%s - %s" % (os.environ['USER'], tmstr)

def image_glob(pat,cols=3,ignore_thumb=True, width=240):
	image_list = []
	for image in glob.glob(pat):
		if 'thumb' in image and ignore_thumb: continue
		# add the absolute path
		else: image_list.append(os.path.abspath(image))
	image_list.sort()
	plot_list = [_imagelinkcpy(plot,width=width) for plot in image_list]
	cols = int(cols)
	return [plot_list[i*cols:i*cols+cols] for i in range(int(math.ceil(len(plot_list) / float(cols))))]

def image_table_from_url_table(table):
	out = []
	for col in table:
		row = [_imagelinkcpy(url) for url in col]
	 	out.append(row)
        return out

def image_table_from_cache(inputcache,cols=3,ignore_thumb=True):
	image_list = []
	for image in inputcache:
		image_list.append(image.url)
	image_list.sort()
	plot_list = [_imagelinkcpy(plot) for plot in image_list]
	cols = int(cols)
	return [plot_list[i*cols:i*cols+cols] for i in range(int(math.ceil(len(plot_list) / float(cols))))]

def wiki_table_parse(file):
	#FIXME assumes table files of the form
	# === title ===
	# ||data||data||
	# ||data||data||

	tabs = []
	titles = []
	tab = []
	for line in open(file).readlines():
		if '===' in line:
			titles.append(line.replace("=",""))
			if tab: tabs.append(tab)
			tab = []
		if '||' in line: tab.append(line.split('||')[1:])
	tabs.append(tab)
	return tabs, titles
	
# PROBABLY DOES NOT EVER NEED TO BE USED DIRECTLY, BUT IS USED IN cbcpage
class _subpage_id(object):
	def __init__(self, id, link_text, closed_flag=0):
		self.id = id
		self.link_text = link_text
		self.closed_flag = closed_flag

###############################################################################
##### CBC WEB PAGE CLASSES ####################################################
###############################################################################

class _link(markup.page):
	def __init__(self, href="", text=""):
		markup.page.__init__(self, mode="strict_html")
		self.a(href=href)
		self.add(text)
		self.a.close()
	def get_content(self):
		return self.content

class _text(markup.page):
	def __init__(self, txt="", bold=False, italic=False):
		markup.page.__init__(self, mode="strict_html")
		if bold: self.b()
		if italic: self.i()
		self.add(txt)
		if bold: self.b.close()
		if italic: self.i.close()
	def get_content(self):
		return self.content

class _imagelink(markup.page):
	def __init__(self, imageurl, thumburl, tag="img", width=240):
		markup.page.__init__(self, mode="strict_html")
		self.add('<a class="fancybox" href=%s target="_blank"><img src=%s width=%d></a>' % (imageurl, thumburl, width))

	def get_content(self):
		return self.content
		
class _imagelinkcpy(markup.page):
	def __init__(self, imagepath, thumbpath=None, tag="img", width=240):
		markup.page.__init__(self, mode="strict_html")
		try: os.mkdir('Images')
		except: pass
		#So that you can give it a url
		#imagepath.replace('file://localhost','').strip()
		imagepath, headers = urllib.urlretrieve(imagepath)
		imgname = os.path.split(imagepath.rstrip('/'))[1]
		shutil.copy(imagepath, 'Images/')
		if not thumbpath:
			# FIXME we cannot assume convert is installed everywhere
			# Is there a python library to do this?
			thumbname = 'Images/' + "thumb_" + imgname
			command = 'convert Images/%s -resize %dx%d -antialias -sharpen 2x2 %s' % (imgname, width, width, thumbname)
			popen = subprocess.Popen(command.split())
			popen.communicate()
			status = popen.returncode
			imgname = 'Images/' + imgname
		else:
			thumbpath.replace('file://localhost','').strip()
			thumbname = os.path.split(thumbpath.rstrip('/'))[1]
			shutil.copy(thumbpath, 'Images/')
			thumbname = 'Images/' + thumbname
			imgname = 'Images/' + imgname
		self.add('<a class="fancybox" href=%s target="_blank"><img src=%s width=%d></a>' % (imgname, thumbname, width))

	def get_content(self):
		return self.content
	

class _table(markup.page):
	def __init__(self, two_d_data, title="", caption="", tag="table", num="1"):
		markup.page.__init__(self, mode="strict_html")
		self.add("<br>")
		if title: 
			self.b("%s. %s" %(num, title.upper()) )
	
		self.table()
		for row in two_d_data:
			self.add('<tr>')
			tdstr = ""
			for col in row:
				tdstr += "<td>%s</td>" % (str(col),)
			self.add(tdstr)
			self.add('</tr>')
		self.table.close()
		if self.caption: self.i("%s. %s" %(num, caption))
		self.add("<br>")

	def get_content(self):
		return self.content			

# PROBABLY DOES NOT EVER NEED TO BE USED DIRECTLY, BUT IS USED IN cbcpage
class _section(markup.page):
	def __init__(self, tag, title="", secnum="1", pagenum="1", level=2, open_by_default=False):
		markup.page.__init__(self, mode="strict_html")
		self.pagenum = pagenum
		self.secnum = secnum
		self._title = title
		self.sections = {}
		self.section_ids = []
		self.level = level 
		self.tag = tag
		self.id = tag + self.secnum
		self.tables = 0
		self.add('<div class="contenu"><h%d id="toggle_%s" onclick="javascript:toggle2(\'div_%s\', \'toggle_%s\');"> %s.%s %s </h%d>' % (level, self.id, secnum, self.id, pagenum, secnum, title, level) )
		if open_by_default:
			style = 'display:block;'
		else:
			style = 'display:none;'
		self.div(id="div_"+secnum , style=style)

	def add_section(self, tag, title="", open_by_default=False):
		secnum = "%s.%d" % (self.secnum, len(self.sections.values())+1)
		self.sections[tag] = _section(tag, title=title, secnum=secnum, pagenum=self.pagenum, level=self.level+1, open_by_default=open_by_default)
		self.section_ids.append([len(self.sections.values()), tag])
		return self.sections[tag]

	def get_content(self):
		self.section_ids.sort()
		out = self.content
		self.div.close()
		self.div.close()
		for num, key in self.section_ids:
			out.extend(self.sections[key].get_content())
		return out
	
	def add_table(self, two_d_data, title="", caption="", tag="table", num=0):
		self.tables += 1
		tabnum = "%s %s.%s.%s" %  ("Table", self.pagenum, self.secnum, str(self.tables))
		table = _table(two_d_data, title=title, caption=caption, tag="table", num=tabnum)
		self.content.extend(table.get_content())
		return self

	def add_link(self, **kwargs):
		link = _link(**kwargs)
		self.content.extend(link.get_content())
		return self

	def add_text(self, **kwargs):
		text = _text(**kwargs)
		self.content.extend(text.get_content())
		return self

# MAIN CBC WEB PAGE CLASS 
class cbcpage(markup.page):

	def __init__(self, title="cbc web page", path='./', css=["//versions.ligo.org/cgit/lalsuite/plain/glue/etc/cbcwebpage.css","https://ldas-jobs.ligo-wa.caltech.edu/~detchar/html/fancybox/source/jquery.fancybox.css?v=2.1.5"], script=None, pagenum=1, verbose=False):
		"""
		"""
		scdict = script_dict(fname='%s/%s' % (path,"toggle.js"))
		if not script: script = scdict[0]
		self.front = ""
		scriptfiles = scdict[1]
		self.verbose = verbose
		self._style = css
		self._title = title
		self._script = script
		self.path = path
		self.pagenum = pagenum

		markup.page.__init__(self, mode="strict_html")	
		self._escape = False
		doctype="""<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">"""
		doctype+="""\n<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">"""
		self.init(title=self._title, css=self._style, script=self._script , doctype=doctype)
		
		self.subpages = {}
		self.subpage_ids = [] 
		self.external_frames = []	
		
		self.sections = {}
		self.section_ids = []

		self.tables = 1
		self.fnames = scriptfiles

	def add_subpage(self, tag, title, link_text=None):
		""" 
		"""

		subpage_num = len(self.subpages.values()) + 1
		if not link_text: link_text=str(subpage_num)

		# tuple including number so it can be sorted later
		self.subpages[tag] = cbcpage(title=title,path=self.path,css=self._style,script=self._script,pagenum=subpage_num)
		self.subpages[tag].add('<table align=right><tr><td align=right onclick="javascript:toggleAllOpen();"><b>Toggle Open</b></td></tr></table>')
		self.subpage_ids.append( [subpage_num, _subpage_id(tag, link_text)] )
		return self.subpages[tag]

	def close_subpage(self,id=None):
		#SECTIONS WILL AUTOMATICALLY BE CLOSED IN WRITE METHOD IF NOT DONE EXPLICITELY
		self.subpage_ids.sort()
		if not id: id = subpage_ids[-1][1].id

		self.subpages[id].div.close()
		self.subpages[id].add("<!-- close div contenu-->")
		self.subpages[id].div.close()
		self.subpages[id].add_footer()

	def add_footer(self):
		#FIXME DO SOMETHING
		pass

	def add_external_frame(self, linkurl, linktext):
		self.external_frames.append([linkurl, linktext])

	def write(self, file_name="index", image = "https://www.lsc-group.phys.uwm.edu/ligovirgo/cbc/public/segments/S5/thomasLegacy.jpg", tag = "CBC"):

		if self.subpage_ids:

			self.div(id_="wrapper")
			self.div(id_="menubar")
			self.div(id_="menu")
			self.subpage_ids.sort()
			# do the sub pages
			for num,secid in self.subpage_ids:
				id = secid.id
				if secid.closed_flag == 0: self.close_subpage(id)
				secfname = file_name + "_" + id
				self.fnames.append(self.subpages[id].write(secfname))
				self.div(class_="menuitem")
				self.add('\t<a class="menulink" href="javascript:loadFrame(\'%s.html\');"> %d: %s </a>\n' % (secfname, num, secid.link_text) )
				self.div.close()
			for i, ext_frame in enumerate(self.external_frames):
				self.div(class_="menuitem")
				self.add('\t<a class="menulink" href=%s# onclick="javascript:loadFrame(\'%s\');"> %d: %s </a>\n' % (ext_frame[0], ext_frame[0], num+i, ext_frame[1]) )
				self.div.close()
			self.div.close()
			self.div.close()
			self.div(id_="ihope")
			self.add('<h2> %s </h2>'%tag)
			self.add('<img width=90 src="%s">'%image)
			self.div.close()
			self.div(id_='header')
			self.add('<h1>' + self._title  +' </h1>')
			self.add('<h3> ' + user_and_date() + ' </h3>')
			self.div.close()
			self.div(id_='iframecontent')
			if self.front: self.add(self.front)
			self.add('<p id="placeholder">Please select a report section on the left.</p>')
			self.div.close()
			self.div.close()

		# do the sections
		self.section_ids.sort()
		for num, key in self.section_ids:
			self.content.extend(self.sections[key].get_content())
		self.fnames.append('%s/%s.html' % (self.path, file_name))
		pagefile = file('%s/%s.html' % (self.path, file_name), 'w')
		pagefile.write(str(self))
		pagefile.close()
		return '%s/%s.html' % (self.path, file_name)
			
	def add_section(self, tag, title="", level=2, open_by_default=False):
		"""
		"""
		secnum = len(self.sections.values()) + 1
		self.section_ids.append([secnum, tag])
		self.sections[tag] = _section(title=title, tag=tag, secnum=str(secnum), pagenum=str(self.pagenum), level=level, open_by_default=open_by_default)
		return self.sections[tag]

	def add_table(self, two_d_data, title="", caption="", tag="table"):
		self.tables += 1
		table = _table(two_d_data, title=title, caption=caption, tag="table", num=str(self.pagenum) + " Table "+str(self.tables))
		self.content.extend(table.get_content())
		return self

	def add_link(self, **kwargs):
		link = _link(**kwargs)
		self.content.extend(link.get_content())
		return self

	def add_text(self, **kwargs):
		text = _text(**kwargs)
		self.content.extend(text.get_content())
		return self

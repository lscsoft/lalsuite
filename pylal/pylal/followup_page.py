import sys, os, socket, re
import glob, math
from glue import cbcwebpage
from glue import lal
from glue import segments
from pylal import fu_utils

#
# useful functions
#

def parse_plot_cache_for_some_images(cache, basepath, tag='.*.png', brk=None):
	out = []
	#FIXME make faster
	for l in open(cache).readlines():
		if re.match(tag, l):
			img = basepath + '/' + l
			out.append((img.strip(), img.replace('.png', '_thumb.png').strip()))
			if brk: return out[0]
	return out


def parse_plot_cache_for_image(cache, basepath, tag):
	return parse_plot_cache_for_some_images(cache, basepath, tag='.*.png', brk=True)


def parse_plot_cache_for_all_images(cache, basepath):
	return parse_plot_cache_for_some_images(cache, basepath)


def cache_parser(cachefile):
	coinc = {}
	f = open(cachefile)
	out_cache = []
	for l in f.readlines():
		if "COINC_INFO" in l:
			c = l.split()
			coinc.setdefault(c[1].replace('COINC_INFO_',''),[]).append(c[4].replace('file://localhost',''))
		else: out_cache.append(lal.CacheEntry(l))
	return coinc, out_cache

class Coinc(object):
	def __init__(self, coinc, search, cache):
		self.cache = cache
		# set up some useful cache views
		self.htqscan_cache = self.parse_cache_by_desc("WPIPELINE_FG_HT_"+search.upper())
		self.seismicqscan_cache = self.parse_cache_by_desc("WPIPELINE_FG_SEIS_RDS_"+search.upper())
		#(CVT)
		self.rdsqscan_cache = self.parse_cache_by_desc("WPIPELINE_FG_RDS_"+search.upper())
		self.plotsnrchisq_cache = self.parse_cache_by_desc("PLOTSNRCHISQ_PIPE__"+search.upper())
		self.plotchia_cache = self.parse_cache_by_desc("PLOTCHIATIMESERIES__"+search.upper())
		self.skymap_cache = self.parse_cache_by_desc("PYLAL_PLOT_INSPIRAL_SKYMAP__"+search.upper())
		self.analyze_qscan_ht_cache = self.parse_cache_by_desc("ANALYSEQSCAN.PY_FG_HT_"+search.upper())
		self.analyze_qscan_rds_cache = self.parse_cache_by_desc("ANALYSEQSCAN.PY_FG_RDS_"+search.upper())
		self.analyze_qscan_seis_cache = self.parse_cache_by_desc("ANALYSEQSCAN.PY_FG_SEIS_RDS_"+search.upper())
		self.flag_cache = self.parse_cache_by_desc("FOLLOWUPQUERYDQ.PY__"+search.upper())
		self.veto_cache = self.parse_cache_by_desc("FOLLOWUPQUERYVETO.PY__"+search.upper())

		f = open(coinc)
		line = f.readlines()
		d = line[1].split()
		self.dir = d[0]
		self.rank = d[1]
		self.cfar = d[2]
		self.coincsnr = d[3]
		self.ifos = d[4]
		self.instruments = d[5]
		self.coinctime = d[6]
		self.coincmass = d[7]
		self.time = {}
		self.snr = {}
		self.chisq = {}
		self.mass1 = {}
		self.mass2 = {}
		for l in line[3:]:
			#don't look at ifos not found in coincidence since the parameters are stolen from another ifo
			d = l.split()
                        if d[1].strip() not in self.ifos: d[3:] = ["0" for i in d[3:]]
			self.time[d[1]] = d[2]
			self.snr[d[1]] = d[3]
			self.chisq[d[1]] = d[4]
			self.mass1[d[1]] = d[5]
			self.mass2[d[1]] = d[6]

	def parse_cache_by_desc(self, tag, cache=None):
		if cache is None: cache = self.cache
		return [l for l in cache if tag in l.description]

	def parse_cache_by_time_and_ifo(self, time, ifo, cache=None):
		if cache is None: cache = self.cache
		return [l for l in cache if float(time) == float(l.segment[0]) and str(ifo) == str(l.observatory)]

	def write_param_table(self, page):
		page.add_section("param", "Parameter table for %s" % (self.coinctime,))
		params = [["<b>RANK</b>", "<b>CFAR</b>",  "<b>TIME</b>", "<b>SNR</b>", "<b>MASS</b>","<b>IFOS</b>","<b>INSTRUMENTS</b>"],[self.rank, self.cfar, self.coinctime, self.coincsnr, self.coincmass, self.ifos, self.instruments]]
		page.sections["param"].add_table(params, title="Coinc Parameter Table", caption="Coinc parameters for the event", tag="coincparamtable")

		params = [["<b>IFO</b>","<b>TIME</b>", "<b>SNR</b>",  "<b>CHISQ</b>", "<b>MASS1</b>", "<b>MASS2</b>"]]
		for ifo, data in self.time.items():
			params.append([ifo, self.time[ifo], self.snr[ifo], self.chisq[ifo], self.mass1[ifo], self.mass2[ifo]])
		page.sections["param"].add_table(params, title="Sngl Parameter Table", caption="Sngl parameters for the event", tag="snglparamtable")

	def add_htqscan(self, page):
		self._add_qscan(page, self.htqscan_cache, self.analyze_qscan_ht_cache, "h(t)")

	def add_seismicqscan(self, page):
		self._add_qscan(page, self.seismicqscan_cache, self.analyze_qscan_seis_cache, "SEISMIC")

	def add_rdsqscan(self, page):
		self._add_qscan(page, self.rdsqscan_cache, self.analyze_qscan_rds_cache, "RDS")

	def _add_qscan(self, page, thiscache, ancache, name):
		job_list = []
		for ifo, time in self.time.items():
			c = self.parse_cache_by_time_and_ifo(time, ifo, thiscache)
			if not c: continue
			else: job_list.append(c)
		if not job_list: return # get out of here if these jobs were not found

		page.add_section(name, "%s Qscan for %s" % (name, self.coinctime,))
		page.sections[name].div("This section gives the %s omega scans and plots that summarize the significance." % (name,))
		img_col = {}
		# since qscans are already by default on web space, they are handled differently
		for ifo, time in self.time.items():
			c = self.parse_cache_by_time_and_ifo(time, ifo, thiscache)
			if not c: continue # get out of here if these jobs were not found

			ca = self.parse_cache_by_time_and_ifo(time, ifo, ancache)
			page.sections[name].add("<a href=%s>LINK TO %s QSCAN</a><br>" % (cbcwebpage.web_path_to_url(c[0].url.replace('file://localhost','')),ifo))
			confile = '%s/configuration.txt' % (c[0].path,)
			try: qconf = fu_utils.omega_config_parser(confile)
			except ValueError:
				print >>sys.stderr, "File %s could not be parsed" % ( confile, )
				continue
			plots = qconf.to_plot_tuple()
			self._add_qscan_plots(page, plots, c, ca, img_col, ifo)

		self._finish_qscan(page, img_col, name)

	def _add_qscan_plots(self, page, plots, c, ca, img_col, ifo):
		# get the analyze qscan cache of images
		cfile = [l.path for l in ca if l.path.endswith('.cache')][0]
		for i,plot in enumerate(plots):
			img_col.setdefault(plot[0],{})
			# first the qscans
			thumb = plot[1].strip().replace(".png",".thumb.png")
			pat = c[0].url.replace('file://localhost','')+'/' + plot[1]
			img_glob = glob.glob(pat)
			pat = c[0].url.replace('file://localhost','')+'/' + thumb
			thumb_glob = glob.glob(pat)
			# now for the analyze qscan stuff
			basename = '/' + os.path.split(cfile)[0].lstrip('/')
			analyze_images = parse_plot_cache_for_some_images(cfile, basename, '.*%s.*' % (plot[0].replace(':','_'),))
			if analyze_images:
				img_glob.extend([im[0] for im in analyze_images])
				thumb_glob.extend([im[1] for im in analyze_images])
			for img, thmb in zip(img_glob, thumb_glob):
				img = '/' + img.lstrip('/')
				thmb = '/' + thmb.lstrip('/')
				img_col[plot[0]].setdefault(ifo,[]).append(cbcwebpage._imagelinkcpy(img, thmb)())

	def _finish_qscan(self, page, img_col, type):
		for name, plot in sorted(img_col.items()):
			#FIXME terrible hack to make nice sections
			if plot:
				chan = name[:6]
				try: page.sections[type].sections[chan]
				except: page.sections[type].add_section(chan,chan)
				page.sections[type].sections[chan].div('<br><hr><b>%s</b>' % (name,))
				for ifo, row in plot.items():
					title = '%s %s' % (ifo, name)
					page.sections[type].sections[chan].add_table([row],title,title + " qscans and plots summarizing significance")


	def add_plotsnrchisq(self, page):
		for ifo, time in self.time.items():
			# Parse plotting codes nearly useless "cache" file
			c = self.parse_cache_by_time_and_ifo(time, ifo, self.plotsnrchisq_cache)
			if not c: return # get out of here if these jobs were not found

		page.add_section("plotsnrchisq", "SNR, Chisq and template time series for %s" % (self.coinctime,))
		img_row = []
		ifo_row = []
		table = []
		for ifo, time in self.time.items():
			# Parse plotting codes nearly useless "cache" file
			c = self.parse_cache_by_time_and_ifo(time, ifo, self.plotsnrchisq_cache)
			cfile = c[0].url.replace('file://localhost','')
			path = '/' + os.path.split(cfile.rstrip('/').lstrip('/'))[0]
			#clist = open(cfile).readlines()
			plots = ['snr-','snr_zoom-', 'rsq-', 'rsq_zoom-', 'chisq-', 'chisq_zoom-', 'PSD-', 'fft_of_template_and_asd-', 'template-', 'white_template-']
			plot_list = []
			for plot in plots:
				img, thumb = parse_plot_cache_for_image(cfile, path, plot)
				plot_list.append(cbcwebpage._imagelinkcpy(img,thumb,plot))
			img_row.append(plot_list)
			ifo_row.append(ifo)
		table.append(ifo_row)
		for row in zip(*img_row): table.append(row)
		page.sections["plotsnrchisq"].add_table(table, "Plots of inspiral stuff", "Plots of snr, snrzoom, rsq, rsqzoom, chisq, chisqzoom, PSD, fft of templates and PSD, template and whitened template by ifo", tag="plotsnrchisq")

	def add_plotchia(self, page):
		c = self.parse_cache_by_time_and_ifo(self.coinctime, self.instruments, self.plotchia_cache)
		if not c: return # if the job didn't finish return

		page.add_section("plotchia", "Coherent Code Plots for %s" % (self.coinctime,))
		img_row = []
		ifo_row = []
		table = []
		plot_list = []

		cfile = c[0].url.replace('file://localhost','')
		path = '/' + os.path.split(cfile.rstrip('/').lstrip('/'))[0]
		#try: clist = open(cfile).readlines()
		#except:
		#	print >>sys.stderr, "couldn't find cachefile %s" % (cfile,)
		#	page.sections["plotchia"].add("<br><b>plot chia job did not finish correctly</b><br>")
		#	return

		for num, plot in enumerate(parse_plot_cache_for_all_images(cfile, path)):
			plot_list.append(cbcwebpage._imagelinkcpy(plot[0],plot[1],"chia"+str(num)))
		# group by 3s
		plot_list = [plot_list[i*3:i*3+3] for i in range(int(math.ceil(len(plot_list) / 3.)))]
		page.sections["plotchia"].add_table(plot_list, "Plots of coherent inspiral stuff", "all of plotchiatimeseries output", tag="plotchia")

	def add_skymap(self,page):
		c = self.parse_cache_by_time_and_ifo(self.coinctime, self.instruments, self.skymap_cache)
		if not c: return # if the job didn't finish return

		page.add_section("skymap", "Sky Map for %s" % (self.coinctime,))
		img_row = []
		ifo_row = []
		table = []
		plot_list = []

		c = self.parse_cache_by_time_and_ifo(self.coinctime, self.instruments, self.skymap_cache)
		cfile = c[0].url.replace('file://localhost','')
		path = '/' + os.path.split(cfile.rstrip('/').lstrip('/'))[0]
		#try: clist = open(cfile).readlines()
		#except:
		#	print >>sys.stderr, "couldn't find cachefile %s" % (cfile,)
		#	page.sections["skymap"].add("<br><b>skymap job did not finish correctly</b><br>")
		#	return

		for num, plot in enumerate(parse_plot_cache_for_all_images(cfile, path)):
			plot_list.append(cbcwebpage._imagelinkcpy(plot[0],plot[1],"skymap"+str(num)))
		# group by 3s
		plot_list = [plot_list[i*3:i*3+3] for i in range(int(math.ceil(len(plot_list) / 3.)))]
		page.sections["skymap"].add_table(plot_list, "Sky map", "lalapps skymap plot", tag="plotchia")

	def add_checklist(self, page):
		page.add_section("checklist", "Detection Checklist for %s" % (self.coinctime,))
		page.sections["checklist"].add("<a href=https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/followup_%s>DETECTION CHECKLIST FOR %s</a><br>" % (self.coinctime,self.coinctime))
		page.sections["checklist"].add('<i>NOTE IF PAGE DOES NOT EXIST CHOOSE "FUCheckListTemplate" FROM THE TEMPLATE SECTION<br>')

	def add_dq(self, page):
		page.add_section("DQ", "Data Quality for %s" % (self.coinctime,))
		page.sections["DQ"].div("This section gives vetoes and flags that were on")

		ca = self.parse_cache_by_time_and_ifo(self.coinctime, self.instruments, self.flag_cache)
		if ca and os.path.isfile(ca[0].path):
			tab, title = cbcwebpage.wiki_table_parse(ca[0].path)
			#FIXME HACK, may stop working
			for t in tab[0]: t[0] = t[0].replace('<rowbgcolor','</td></tr><tr bgcolor') + '<td>' + t[0]
			page.sections["DQ"].add_table(tab[0], 'dq flags', 'dq flags: Yellow denotes before, red during and green after')
		else: page.sections["DQ"].div("Job did not finish")

		ca = self.parse_cache_by_time_and_ifo(self.coinctime, self.instruments, self.veto_cache)
		if ca and os.path.isfile(ca[0].path):
			tab, title = cbcwebpage.wiki_table_parse(ca[0].path)
			#FIXME HACK, may stop working
			for t in tab[0]: t[0] = t[0].replace('<rowbgcolor','</td></tr><tr bgcolor') + '<td>' + t[0]
			page.sections["DQ"].add_table(tab[0], 'vetoes', 'vetoes: Yellow denotes before, red during and green after')
		else: page.sections["DQ"].div("Job did not finish")

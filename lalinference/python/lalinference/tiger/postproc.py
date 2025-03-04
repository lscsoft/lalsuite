"""
Postprocessing for the tiger pipeline
"""

__author__ = "Tjonnie Li"
__credits__ = ["Tjonnie Li"]
__maintainer__ = "Tjonnie Li"
__email__ = "tjonnie.li@ligo.org"
__status__ = "Production"

###############################################################################
#
# USER DEFINED DATA
#
###############################################################################

clusters = {'atlas':'titan2.atlas.aei.uni-hannover.de', \
		'nemo':'pcdev2.phys.uwm.edu', \
		'cit':'ldas-pcdev2.ligo.caltech.edu', \
		'llo':'ldas-pcdev2.ligo-la.caltech.edu', \
		'lho':'ldas-pcdev1.ligo-wa.caltech.edu'}

color_default = ['blue','red','green']
hatch_default = ['..','//','||']

###############################################################################
#
# LOAD LIBRARIES
#
###############################################################################

from configparser import ConfigParser
from pickle import dump, load
from datetime import datetime
from itertools import combinations
from itertools import cycle
from matplotlib import use, rcParams, __version__ as mpl_version
use('Agg')
from matplotlib.pyplot import clf, figure
from os import access, path, R_OK, makedirs
from scipy.stats import ks_2samp
from scipy.stats.mstats import mquantiles
from subprocess import Popen, PIPE
from sys import exit, stdout, version_info, float_info
from time import time
from numpy import sqrt, array, empty, hstack, min, max, reshape, shape, loadtxt, vstack, append, arange, random, column_stack, concatenate, savetxt, log, exp, size, zeros, argmax, argsort, sort, sum, subtract, array_split
from re import findall

py_version = version_info[:2]
if py_version < (2, 7):
	from optparse import OptionParser
else:
	from argparse import ArgumentParser

###############################################################################
#
# PLOTTING OPTIONS
#
###############################################################################

fig_width_pt = 3*246.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]

fontsize = 'font.size' if mpl_version >= '1.5.0' else 'text.fontsize'
params = {'backend': 'PDF',
          'axes.labelsize': 24,
          fontsize: 24,
          'legend.fontsize': 20,
          'xtick.labelsize': 24,
          'ytick.labelsize': 24,
          'axes.grid' : True,
          'text.usetex': True,
          'savefig.dpi' : 200,
          'lines.markersize' : 14,
          'figure.figsize': fig_size}

rcParams.update(params)

scriptfilename = path.realpath(__file__)
###############################################################################
#
# PLOTTING OPTIONS
#
###############################################################################

def main():
	"""
	-------------------------------------------------------------------------------
	"""


	###############################################################################
	#
	# ARGUMENT PARSING
	#
	###############################################################################

	# SETUP ARGUMENT PARSER
	if py_version < (2, 7):
		parser = OptionParser()
		parser.add_option('-c', '--config', metavar='FILE', dest='configfile', \
				type=str, help='configuration file')
		parser.add_option('-g', '--generate', metavar='FILE', dest='examplecfg', \
				type=str, help='generate example configuration file')
		parser.add_option('-p', '--preprocess', metavar='FILE', dest='preprocessfile', \
				type=str, help='preprocess data before sending to main script (NB: NOT WORKING!)')
		(args, options) = parser.parse_args()
	else:
		parser = ArgumentParser(description='Post-processing for TIGER')
		parser.add_argument('-c', '--config', metavar='FILE', dest='configfile', \
				type=str, help='configuration file')
		parser.add_argument('-g', '--generate', metavar='FILE', dest='examplecfg', \
				type=str, help='generate example configuration file')
		parser.add_argument('-p', '--preprocess', metavar='FILE', dest='preprocessfile', \
				type=str, help='preprocess data before sending to main script (NB: NOT WORKING!)')
		args = parser.parse_args()

	if (args.configfile == None) and (args.examplecfg == None) and (args.preprocessfile == None):
		exit("Specify either -c/--config, -g/--generate or -p/--preprocess")
	elif (args.configfile != None):
		stdout.write("****************************************\nTIGER post-process\n%s\n****************************************\n" % args.configfile)
		TigerPostProcess(args.configfile)
	elif (args.examplecfg != None):
		TigerCreateExampleConfigFile(args.examplecfg)
	elif (args.preprocessfile != None):
		TigerPreProcess(args.preprocessfile)
	else:
		exit('Unknown options - check input')
		exit(0)

def TigerPostProcess(configfile):
	"""
	Main TIGER post processing function. Calls the relevant function for source
	collection, odds ratio calculations, and plotting
	"""

	# CONFIGURATION FILE: CHECKING FILE
	stdout.write("Configuration file: checking file\n")
	if access(configfile, R_OK):
		config = ConfigParser()
		config.read(configfile)
	else:
		exit("Configuration file: checking file - file does not exist. Abort script\n")

	# CONFIGURATION FILE: READING
	stdout.write("Configuration file: reading\n")
	testcoef = config.get('run','hypotheses').split(',')
	localdest = config.get('run','localdest')
	N_sources=[int(i) for i in config.get('hist','nsources').split(',')]
	types = config.get('plot','types').split(',')
	engine=config.get('run','engine').split(',')
	runid = config.get('run','runid')
	labels=config.get('fetch','labels').split(',')

	latexlabels_str = config.get('fetch','latexlabels')
	matches=findall(r'\$(.+?)\$',latexlabels_str)
	latexlabels=["$%s$"%s for s in matches]

	# CREATE DESTINATION DIRECTOR
	ensure(localdest)

	###############################################################################
	#
	# FETCH DATA
	#
	###############################################################################

	# FETCHING DATA: GETTING OPTIONS
	fetch_type = config.get('fetch','type').split(',')
	fetch_loc = config.get('fetch','locs').split(',')
	fetch_snrthres = config.getfloat('cut','snrthres')
	fetch_bayesgrthres = config.getfloat('cut','bayesgrthres')
	fetch_seed = config.getint('fetch','seed')
	if len(fetch_type) != len(fetch_loc):
		exit('Locations and types are not of the same length - exit\n')

	# FETCHING DATA
	tigerruns = []
	for i in range(len(fetch_type)):
		if fetch_type[i] == 'source':
			# READING LOCATION FILES
			locfp = open(fetch_loc[i],'r')
			locations = [k.split() for k in locfp.readlines()]
			stdout.write("Fetching data: %s\n" % fetch_loc[i])
			# CREATING CLASS INSTANCES
			tigerruns.append(TigerSet(locations, labels[i],latexlabels[i],testcoef,engine[i]))
			# FETCH DATA
			tigerruns[-1].searchsources() # SEARCH FOR SOURCES
			tigerruns[-1].pullbayes() # PULL BAYESFACTORS
			tigerruns[-1].pullsnr() # PULL SNR
			#tigerruns[-1].preprocess() # preprocess data on clusters (NOT READY YET)
		elif fetch_type[i] == 'pickle':
			tigerruns.append(LoadPickledTigerSet(fetch_loc[i]))
			tigerruns[-1].latexlabel=latexlabels[i]
		else:
			exit('Fetch has to be either source or pickle. Check the configuration file - Abort')

	# SAVING DATA BEFORE CUT
	html_files_sources = []
	stdout.write('Saving data\n')
	for i in range(len(tigerruns)):
		html_files_sources.append(path.join(localdest,tigerruns[i].label+'.pickle'))
		tigerruns[i].savetopickle(html_files_sources[-1])
		stdout.write('... Saving data: %s\n' % html_files_sources[-1])
		html_files_sources.append(path.join(localdest,tigerruns[i].label+'.dat'))
		tigerruns[i].savetoascii(html_files_sources[-1])
		stdout.write('... Saving data: %s\n' % html_files_sources[-1])

	# APPLY CUTS
	stdout.write("Applying cuts\n")
	for tg in tigerruns:
		tg.applycut(st=fetch_snrthres, bt=fetch_bayesgrthres)

	# SHUFFLE REALISATIONS
	#if fetch_seed != 0:
		#print 'Shuffling realisation'
		#for tg in tigerruns:
			#tg.shufflesources(seed=fetch_seed)

	# SAVING DATA AFTER CUT
	#print 'Saving data - after cut'
	#for i in range(len(tigerruns)):
		#html_files_sources.append(path.join(localdest,tigerruns[i].label+'_cut'+'.pickle'))
		#tigerruns[i].savetopickle(html_files_sources[-1])
		#html_files_sources.append(path.join(localdest,tigerruns[i].label+'_cut'+'.dat'))
		#tigerruns[i].savetoascii(html_files_sources[-1])

	###############################################################################
	#
	# PRODUCE PLOTS
	#
	###############################################################################

	stdout.write("Producing plots\n")

	# FETCH OPTIONS
	plot_hist = config.getboolean('plot','hist')
	plot_snrVSodds = config.getboolean('plot','snrvsodds')
	plot_cumbayes = config.getboolean('plot','cumbayes')
	plot_cumfreq = config.getboolean('plot','cumfreq')
	hist_odds_bins=config.getint('hist','bins')

	# STATISTICS
	html_data_efficiencies = []
	betas = [0.95,0.99]
	if len(tigerruns)>1:
		for ns in N_sources:
			tmp = TigerCalculateEfficiency(tigerruns,N=ns,beta=betas)
			html_data_efficiencies.append(tmp)
		html_data_efficiencies = array(html_data_efficiencies)

		html_data_ks =[[TigerKSandPvalue(tigerruns[0],j,N=ns) for ns in N_sources] for j in tigerruns[1:] if len(tigerruns)>1]
	#for ns in N_sources:
		#html_data_ks.append(TigerKSandPvalue(tigerruns[0],tigerruns[1],N=ns))

	# HISTOGRAMS
	html_files_hist = []
	if plot_hist==True:
		stdout.write("... Histograms\n")
		# GET XLIM FOR HISTOGRAMS (IF SPECIFICIED)
		hist_xlim = [float(item) for item in config.get('hist','xlim').split(',')]
		if hist_xlim == [0.0,0.0]:
			# SETTING HISTOGRAM LIMITS (ADD 10% TO EACH EDGE)
			tmplist = empty(0)
			for tg in tigerruns:
				for ns in N_sources:
					tmplist = hstack((tmplist,tg.odds(ns)))

			# CALCULATE XLIM FROM [MIN,MAX] + 10 PERCENT
			xlim = array([min(tmplist),max(tmplist)])
			xlim = array([xlim[0]-0.1*(xlim[1]-xlim[0]),xlim[1]+0.1*(xlim[1]-xlim[0])])
			del tmplist
		else:
			xlim = hist_xlim

		# CREATE HISTOGRAMS
		for ns in N_sources:
			clf()
			fig = figure(ns)
			ax = fig.add_subplot(111)
			TigerCreateHistogram(tigerruns,ax,N=ns, xlim=xlim, bins=hist_odds_bins)
			for typ in types:
				html_files_hist.append(path.join(localdest,'hist_'+runid+"_"+str(ns)+"sources."+typ))
				fig.savefig(html_files_hist[-1],bbox_inches='tight')

	# SNR VS ODDS
	html_files_snrvsodds = []
	if plot_snrVSodds == True:
		stdout.write("... SNR vs. Odds\n")
		clf()
		fig = figure()
		ax = fig.add_subplot(111)
		TigerCreateSNRvsOdds(tigerruns, ax)
		for ext in types:
			html_files_snrvsodds.append(path.join(localdest,"snrVSodds_"+runid+"."+ext))
			fig.savefig(html_files_snrvsodds[-1], bbox_inches='tight')

	# CUMULATIVE FREQUENCY PLOTS
	html_files_cumfreq = [[]]*len(tigerruns)
	if plot_cumfreq == True:
		stdout.write("... Cumulative frequency\n")
		for i in range(len(tigerruns)):
			html_files_cumfreq[i] = []
			clf()
			fig = figure()
			ax = fig.add_subplot(111)
			TigerCreateCumFreq(tigerruns[i], ax)
			for ext in types:
				html_files_cumfreq[i].append(path.join(localdest,"cumfreq_"+tigerruns[i].label+"."+ext))
				fig.savefig(html_files_cumfreq[i][-1],bbox_inches='tight')

	# CUMULATIVE BAYES PLOTS
	html_files_cumbayes = [[]]*len(tigerruns)
	N_sources = array(N_sources)
	stdout.write("... Cumulative Bayes\n")
	if plot_cumbayes == True:
		N_sources_cats = array(N_sources[N_sources!=1])
		for i in range(len(tigerruns)):
			html_files_cumbayes[i] = [[]]*len(N_sources_cats)
			for j in range(len(N_sources_cats)):
				html_files_cumbayes[i][j] = []
				for k in range(tigerruns[i].nsources/N_sources_cats[j]):
					clf()
					fig = figure()
					ax = fig.add_subplot(111)
					TigerCreateCumBayes(tigerruns[i],ax,k,N_sources_cats[j])
					for ext in types:
						html_files_cumbayes[i][j].append(path.join(localdest,"cumbayes_"+tigerruns[i].label+"_catalog"+str(k)+"_sourcespercat_"+str(N_sources_cats[j])+"."+ext))
						fig.savefig(html_files_cumbayes[i][j][-1], bbox_inches='tight')

	###############################################################################
	#
	# MAKE HTML OUTPUT PAGE
	#
	###############################################################################

	# TODO - CREATE SIMPLE OUTPUT PAGE

	stdout.write("Creating html pages\n")
	htmlfile=open(path.join(localdest,'index.html'), 'w')
	htmltxt=\
	"""
	<!DOCTYPE html>
	<html>

	<head>
	<title> TIGER post-processing page </title>
	</head>

	<body>

	<h1>TIGER post-processing page</h1>
	"""

	timesec = time()
	timedate = datetime.fromtimestamp(timesec).strftime('%Y-%m-%d %H:%M:%S')
	htmltxt+="Last updated on "+timedate

	# ADD SOURCE LOCATIONS
	#htmltxt+=\
	#"""
	#<h2>Sources</h2>
	#"""
	#for hfs in html_files_sources:
		#hfsbase = path.basename(hfs)
		#htmltxt+="<a href='"+hfsbase+"'>"+hfsbase+"</a>&nbsp;&nbsp;&nbsp;&nbsp;"

	if len(tigerruns)>1:
		htmltxt+=\
		"""
		<h2>Statistics</h2>
		"""
		htmltxt+=\
		"""
		<h3>Efficiencies</h3>
		"""
		htmltxt+=\
		"""
		<table border="1">
		<tr>
		<td>catsize</td>
		"""
		for ns in N_sources:
			htmltxt+='<td>'+str(ns)+'</td>'
		htmltxt+='</tr>'

		background = 0
		for i in range(len(tigerruns)):
			if i != background: # hardcoded for the time being
				htmltxt+='<tr>'
				htmltxt+='<td>'+tigerruns[i].label+'</td>'
				for j in range(len(N_sources)):
					if i < background:
						effstr="%.2f "%html_data_efficiencies[j,i,0]+' | '+ '%.2f'% html_data_efficiencies[j,i,1]
					else:
						effstr="%.2f"%html_data_efficiencies[j,i-1,0]+' | '+'%.2f'%html_data_efficiencies[j,i-1,1]
					htmltxt+='<td>'+effstr+'</td>'
				htmltxt+='</tr>'
		htmltxt+='</table>'

		htmltxt+=\
		"""
		<h3>KS and p value</h3>
		"""
		htmltxt+=\
		"""
		<table border="1">
		<tr>
		<td></td>
		<td>ks stat</td>
		<td>p value</td>
		</tr>
		"""
		for i in range(len(tigerruns)-1):
			for j in range(len(N_sources)):
				htmltxt+='<tr>'
				htmltxt+='<td>'+'N='+str(N_sources[j])+'</td>'
				htmltxt+='<td>'+'%.2f'%html_data_ks[i][j][0]+'</td>'
				htmltxt+='<td>'+'%.2f'%html_data_ks[i][j][1]+'</td>'
				htmltxt+='</tr>'
		htmltxt+='</table>'

	# HISTOGRAMS
	htmltxt+=\
	"""
	<h2>Histograms</h2>
	"""
	for item in html_files_hist:
		if ".png" in item:
			itembase = path.basename(item)
			htmltxt+="<a href='"+itembase+"'>"+"<img src='"+itembase+"' width=512></a>"

	# SNR VS ODDS
	htmltxt+=\
	"""
	<h2>SNR vs odds</h2>
	"""
	for item in html_files_snrvsodds:
		if ".png" in item:
			itembase = path.basename(item)
			htmltxt+="<a href='"+itembase+"'>"+"<img src='"+itembase+"' width=512></a>"

	# PREPARE INDIVIDUAL PAGES
	htmltxt+=\
	"""
	<h2>Individual runs</h2>
	"""
	for i in range(len(tigerruns)):
		htmltxt+="<a href='"+tigerruns[i].label+".html"+"'>"+tigerruns[i].label+"</a>&nbsp;&nbsp;&nbsp;&nbsp;"

	# CLOSE
	htmltxt+=\
	"""
	</body>
	</html>
	"""

	# CREATE PAGES FOR INIDVIDUAL RUNS
	for i in range(len(tigerruns)):
		indhtmlfile=open(path.join(localdest,tigerruns[i].label+'.html'), 'w')
		indhtmltxt=\
		"""
		<!DOCTYPE html>
		<html>

		<head>
		<title> TIGER post-processing page </title>
		</head>

		<body>
		"""

		indhtmltxt+="<h1>"+tigerruns[i].label+"</h1>"

	# CUMFREQ
		indhtmltxt+=\
		"""
		<h2>Cumulative frequency highest Bayes factor</h2>
		"""
		for item in html_files_cumfreq[i]:
			if ".png" in item:
				itembase = path.basename(item)
				indhtmltxt+="<a href='"+itembase+"'>"+"<img src='"+itembase+"' width=512></a>"

		# CUMBAYES
		indhtmltxt+=\
		"""
		<h2>Cumulative Bayes factors (per catalog)</h2>
		"""
		N_sources_cats = array(N_sources[N_sources!=1])
		for j in range(len(html_files_cumbayes[i])):
			indhtmltxt+="<h2>"+str(N_sources_cats[j])+" sources per catalog"+"</h2>"
			for item in html_files_cumbayes[i][j]:
				if ".png" in item:
					itembase = path.basename(item)
					indhtmltxt+="<a href='"+itembase+"'>"+"<img src='"+itembase+"' width=512></a>"

		indhtmlfile.write(indhtmltxt)


	htmlfile.write(htmltxt)

	#############################################################################
	#
	# COPY TO REMOTE CLUSTER
	#
	#############################################################################


	#############################################################################
	#
	# EXIT MAIN
	#
	#############################################################################

	print('Script finished')
	"""
	------------------------------------------------------------------------------
	"""

###############################################################################
#
# TIGERRUN CLASS DEFINITIONS
#
###############################################################################

class TigerRun:
	"""
	Class for a specific TIGER run
	"""

	def __init__(self, cluster, directory, engine, subhyp):
		"""
		Initialisation of class
		"""
		self.cluster = cluster
		if directory[-1] != '/':
			self.directory = directory+'/'
		else:
			self.directory = directory
		self.engine = engine
		self.subhyp = subhyp
		self.bayesfiles = []
		self.snrfiles = []
		self.detectors = []
		self.posteriorfiles = []
		self.gps = []
		self.nsources = 0
		self.bayes = []
		self.snr = []
		self.posteriorsamples =[]

	def searchsources(self):
		"""
		Search for sources completed sources in the specified cluster and location
		"""

		stdout.write("... searching sources: %s\t%s\r" % (self.cluster,self.directory))
		stdout.flush()
		# GET BAYESFILES
		command = ""
		if self.cluster in clusters.keys():
			command += "gsissh -C "+clusters[self.cluster]+" '"
		command+="for i in `ls -a "+self.directory+"/"+self.subhyp[0]
		if self.engine == 'lalnest':
			command+="/nest/"
		elif self.engine == 'lalinference':
			command+="/posterior_samples/"
		command+="*_B.txt | xargs -n1 basename`; do "
		for sh in self.subhyp:
			if self.engine == 'lalnest':
				command+="test -f "+self.directory+"/"+sh
				command+="/nest/"
			elif self.engine == 'lalinference':
				command+="test -f "+self.directory+"/"+sh.replace('phi','chi')
				command+="/posterior_samples/"
			command+="$i && "
		command=command+"echo $i; done"
		if self.cluster in clusters.keys():
			command += "'"
		p = Popen(command, stdout=PIPE, stderr=PIPE,shell=True)
		lines = p.stdout.readlines()
		self.bayesfiles = array([k.strip('\n') for k in lines])
		del lines

		# GET DETECTOR INFORMATION
		if self.engine == 'lalnest':
			self.detectors = array([k.split('_')[2].split('.')[0] for k in self.bayesfiles])
		elif self.engine == 'lalinference':
			self.detectors = array([k.split('_')[1] for k in self.bayesfiles])

		# GET POSTERIOR FILES
		if self.engine == 'lalnest':
			self.posteriorfiles = array([k.replace('_B.txt','') for k in self.bayesfiles])
		elif self.engine == 'lalinference':
			self.posteriorfiles = array([k.replace('_B.txt','') for k in self.bayesfiles])

		# EXTRACT GPS TIMES
		if self.engine == 'lalnest':
			self.gps=array([k.split('_')[1].split('.')[0] for k in self.bayesfiles])
		elif self.engine == 'lalinference':
			self.gps=array([k.split('_')[2].split('-')[0] for k in self.bayesfiles])
		#self.gps=hstack((self.gps,gpstmp))

		# GET SNR FILES
		#snrfilename = "snr_%s_%.1f.dat"%(self.bayesfiles[k].split('_')[1],int(self.gps[k]))
		#self.snrfiles = array(["snr_%s_%.1f.dat"%(self.bayesfiles[k].split('_')[1],int(self.gps[k])) for k in range(len(self.bayesfiles))])
		#self.snrfiles = array(["snr_%s_%.1f.dat"%(x,int(y)) for x,y in zip(self.detectors, self.gps)])
		self.snrfiles = array([k.replace('_','-').replace('-B.txt','_snr.txt').replace('posterior','lalinferencenest-%s'%(k.split('-')[1].split('.')[0])).replace('%s'%(k.split('-')[0].split('_')[-1]),'%.1f'%(float(k.split('-')[0].split('_')[-1]))) for k in self.bayesfiles])
		#self.snrfiles = array(["snr_H1L1V1_"+item+".0.dat" for item in self.gps])

		self.nsources = len(self.bayesfiles)
		stdout.write("... searching sources: %s\t%s - %d sources\n" % (self.cluster,self.directory, self.nsources))
		#print "... ", self.cluster,self.directory, " - ", self.nsources, " sources found"

	def pullbayes(self):
		"""
		Download Bayes factors from remote location. Only works after sources are
		found by the searchsources function.
		"""
		stdout.write("... pulling Bayes factors: %s\t%s\r" % (self.cluster,self.directory))
		stdout.flush()
		if self.nsources == 0:
			stdout.write("... pulling Bayes factors: %s\t%s - nothing to fetch\n" % (self.cluster,self.directory))
		else:
			# new command
			files=""
			for item in self.bayesfiles:
				files+=" "+self.directory+"{"+",".join(self.subhyp)+"}"
				if self.engine == 'lalnest':
					files += "/nest/"
				elif self.engine == 'lalinference':
					files += "/posterior_samples/"
				files += item
			nsplit = 3 if self.nsources > 1000 else 1
			for fsplit in array_split(files.split(),nsplit):
				command="%s %s '%s %s'"%("gsissh -C",clusters[self.cluster],"cat"," ".join(fsplit)) if self.cluster in clusters.keys() else "%s %s"%("cat"," ".join(fsplit))
				p = Popen(command, stdout=PIPE, stderr=PIPE,shell=True)
				self.bayes.extend(loadtxt(p.stdout, usecols=[0]))
			self.bayes = reshape(self.bayes, (self.nsources,len(self.subhyp)))
			savetxt('test.dat',self.bayes);
		stdout.write("... pulling Bayes factors: %s\t%s - %d x %d dimension\n" % (self.cluster,self.directory, shape(self.bayes)[0],shape(self.bayes)[1]))

	def pullsnr(self):
		"""
		Download SNRs from remote location. Only works after sources are found by
		the searchsources function.
		"""
		stdout.write("... pulling SNRs: %s\t%s\n" % (self.cluster,self.directory))
		# GET SNRS
		if self.nsources == 0:
			stdout.write("... pulling SNRs: %s\t%s - nothing to fetch\n" % (self.cluster,self.directory))
		else:
			command = ""
			if self.cluster != 'local':
				command+="gsissh -C "+clusters[self.cluster]+" '"
			command+="cat"
			snrfilescomma=",".join(self.snrfiles)
			command += " "+self.directory+self.subhyp[0]+"/engine/"+"{"+snrfilescomma+"}"
			if self.cluster != 'local':
				command+="'"
			p = Popen(command, stdout=PIPE, stderr=PIPE,shell=True)
			#self.snr = array([float(k.strip('Network:').strip()) for k in p.stdout.readlines() if k.find('Network')!=-1])
			snrrawdata = p.stdout.readlines()
			if snrrawdata == []:
				stdout.write("... Warning: Cannot find SNR file. Writing zeros\n")
				self.snr = zeros(self.nsources)
				return
			count = 0
			for k in range(len(self.gps)):
				ndet = len(self.detectors[k])/2
				if ndet == 1:
					self.snr.append(float(snrrawdata[count+ndet-1].strip('%s:'%(self.detectors[k])).strip()))
					count+=1
				else:
					self.snr.append(float(snrrawdata[count+ndet].strip('Network:').strip()))
					count+=ndet+1
			self.snr = array(self.snr)
		stdout.write("... pulling SNRs: %s\t%s - %d sources\n" % (self.cluster,self.directory, len(self.snr)))

	def pullposteriors(self):
		"""
		Download posterior pdfs from remote location. Only works after sources are
		found by the searchsources function.
		NB: NOT READY FOR USE YET!
		"""
		if self.nsources == 0:
			print("Nothing to fetch")
		else:
			# FIND THE NUMBER OF LINE FIRST SO THAT ONE CAN DIFFERENTIATE FILES
			command = ""
			if self.cluster != 'local':
				command+="gsissh -C "+clusters[self.cluster]+" '"
			command+="wc -l"
			posteriorfilescomma=",".join(self.posteriorfiles)
			command += " "+self.directory+self.subhyp[0]+"/nest/"+"{"+posteriorfilescomma+"}"
			if self.cluster != 'local':
				command+="'"
			p = Popen(command, stdout=PIPE, stderr=PIPE,shell=True)
			posteriorfilesizes = loadtxt(p.stdout, usecols=[0])[:-1]
			print(shape(posteriorfilesizes))

			# PULLING THE POSTERIOR SAMPLES
			command = ""
			if self.cluster != 'local':
				command+="gsissh -C "+clusters[self.cluster]+" '"
			command+="cat"
			posteriorfilescomma=",".join(self.posteriorfiles)
			command += " "+self.directory+self.subhyp[0]+"/nest/"+"{"+posteriorfilescomma+"}"
			if self.cluster != 'local':
				command+="'"
			p = Popen(command, stdout=PIPE, stderr=PIPE,shell=True)
			self.posteriorsamples = loadtxt(p.stdout)
			print(shape(self.posteriorsamples))


	def applycut(self, st=8.0, bt=32.0):
		"""
		Apply a Bayes factor and SNR cut on the data
		"""
		# SET CONDITION
		cond = ((self.snr>st) | (self.snr==0.0)) & (self.bayes[:,0]>bt)

		# UPDATE DATA
		self.bayes = self.bayes[cond,:]
		self.snr = self.snr[cond]
		self.gps = self.gps[cond]
		self.nsources=shape(self.bayes)[0]

	def savetopickle(self,dest):
		"""
		Pull data from remote or local locations
		"""

		# SAVING TO PICKLE
		f = open(dest,'wb')
		dump(self,f,2)

class TigerSet:
	"""
	CLASS TO CONTAIN A SET OF RUNS UNDER A COMMON CHARACTERISTIC (E.G.
	BACKGROUND, FOREGROUND).
	"""
	def __init__(self,locations,label,latexlabel,testcoef,engine):
		"""
		INITIALISATION OF CLASS
		"""

		# CALCULATING SUBHYPOTHESES
		self.testcoef = testcoef # testing coefficients
		self.subhyp=[]
		for i in range(len(self.testcoef)):
			s = combinations(testcoef, i+1)
			self.subhyp.extend(["".join(j) for j in s])
		self.subhyp.insert(0,'GR')
		self.nsubhyp = len(self.subhyp)

		# SET ENGINE
		if engine in ['lalinference','lalnest']: # engine: either lalnest or lalinference
			self.engine = engine
		else:
			exit('Engine not recognised')

		# CREATE LOCATIONS
		self.locations = [TigerRun(loc[0].lower(),loc[1],self.engine,self.subhyp) for loc in locations if loc[0].lower() in clusters.keys()+['local']]

		# SET LABELS
		self.label = label
		self.latexlabel = latexlabel

		# DATA STRUCTURES TO BE FILLED USING FUNCTIONS
		self.nsources=0
		self.bayes = empty((0,len(self.subhyp)))
		self.snr = empty(0)

	def searchsources(self):
		"""
		FIND COMPLETED JOBS AND GET THE FILESNAMES
		"""

		for loc in self.locations:
			loc.searchsources()
			self.nsources += loc.nsources

		stdout.write("--> searching sources: Total number of sources - %d\n" % self.nsources)

	def pullbayes(self):
		"""
		PULL DATA FROM REMOTE OR LOCAL LOCATIONS
		"""

		for loc in self.locations:
			if loc.nsources > 0:
				loc.pullbayes()
				self.bayes=vstack((self.bayes,loc.bayes))
		stdout.write("--> pulling Bayes factors: matrix size - %d x %d\n" % (shape(self.bayes)[0], shape(self.bayes)[1]))

	def pullsnr(self):
		"""
		PULL DATA FROM REMOTE OR LOCAL LOCATIONS
		"""
		for loc in self.locations:
			if loc.nsources > 0:
				loc.pullsnr()
				self.snr=hstack((self.snr,loc.snr))
		stdout.write("--> pulling SNRs: total number of sources %d\n" % (len(self.snr)))

	def pullposteriors(self):
		"""
		PULLING THE POSTERIOR FILES
		"""
		for loc in self.locations:
			if loc.nsources > 0:
				loc.pullposteriors()

	def odds(self,N=1):
		"""
		CALCULATE THE TOTAL ODDS FROM DATA STRUCTURE
		"""
		odds = empty(0)
		for i in range(self.nsources/N):
			odds = append(odds,OddsFromBayes(self.bayes[i*N:(i+1)*N,:], len(self.testcoef)))
		return odds

	def applycut(self, st=8.0,bt=32.0):
		"""
		APPLY DETECTION CUTS (SNR AND BAYESGR)
		"""
		stdout.write("... %s - before cut: %d sources\n" % (self.label,self.nsources))

		for loc in self.locations:
			if loc.nsources > 0:
				loc.applycut(st=st, bt=bt)

		# RELOAD BAYES, SNR AND GPS
		self.bayes = empty((0,len(self.subhyp)))
		self.snr = empty(0)
		for loc in self.locations:
			if loc.nsources > 0:
				self.bayes=vstack((self.bayes,loc.bayes))
				self.snr=hstack((self.snr,loc.snr))
		self.nsources=shape(self.bayes)[0]

		stdout.write("... %s - after cut: %d sources\n" % (self.label,self.nsources))

	def shufflesources(self,seed):
		"""
		APPLY DETECTION CUTS (SNR AND BAYESGR)
		"""

		print("... Shuufling with seed", str(seed))

		# CREATE NEW ORDER OF SOURCES
		if (self.nsources != 0) and (self.nsources == shape(self.bayes)[0]):
			order = arange(self.nsources)
			random.seed(seed)
			random.shuffle(order)

			# UPDATE DATA
			self.bayes = self.bayes[order,:]
			self.snr = self.snr[order]

	def savetopickle(self,dest):
		"""
		SAVE DATA TO PICKLE FOR QUICK FUTURE LOADING
		"""

		# SAVING TO PICKLE
		f = open(dest,'wb')
		dump(self,f,2)

	def savetoascii(self,filename):
		"""
		SAVE DATA TO ASCII FILE FOR QUICK FUTURE LOADING.
		NB: BROKEN - DO NOT USE!
		"""
		# PRINT DATA FOR FUTURE USE
		savedata = empty((0,4+len(self.subhyp))) # cluster, directory, gps, bayes, ...., netsnr
		for loc in self.locations:
			if loc.nsources > 0:
				locsdata = column_stack(([loc.cluster]*loc.nsources, [loc.directory]*loc.nsources))
				savedatatmp = column_stack((locsdata,loc.gps,loc.bayes,loc.snr))
				savedata = vstack((savedata,savedatatmp))
		header = array(concatenate((["#cluster",'folder','gps'],self.subhyp,['netsnr'])))
		savedata = vstack((header,savedata))
		savetxt(filename,savedata,delimiter='\t',fmt='%s')

	def preprocess(self):
		"""
		PREPROCESS DATA ON THE CLUSTERS BEFORE SENDING DATA TO MASTER.
		NB: UNDER DEVELOPMENT
		"""
		for loc in self.locations:
			# CREATE PREPROCESSFILE
			print('Writing preprocessfile')
			preprocesstmp = 'tig_preprocess_tmp.txt'
			preprocessfp = open(preprocesstmp,'w')
			preprocessfp.write('local')
			preprocessfp.write('\t')
			preprocessfp.write(loc.directory)
			preprocessfp.write('\t')
			preprocessfp.write(self.engine)
			preprocessfp.write('\t')
			preprocessfp.write(str(loc.subhyp))
			preprocessfp.close()
			print('Writing preprocessfile - done')

			# UPLOAD FILE AND RUN: TESTING PHASE
			command = ""
			command += "gsiscp -C"
			command += " "
			command += scriptfilename
			command += " "
			command += preprocesstmp
			command += " "
			command += clusters[loc.cluster]+":~/"
			print(command)
			p = Popen(command, stdout=PIPE, stderr=PIPE,shell=True)
			p.communicate()
			if p.returncode == 0:
				command = "gsissh "+clusters[loc.cluster]+" 'python "+path.basename(scriptfilename)+" -p "+preprocesstmp+"'"
				print(command)
				p = Popen(command, stdout=PIPE, stderr=PIPE,shell=True)
				print(p.stdout.read())
				print(p.stderr.read())
			exit(0)


###############################################################################
#
# FUNCTIONS
#
###############################################################################

def LoadPickledTigerSet(filename):
	"""
	LOAD FROM PICKLE
	"""
	print('... Loading file from pickle',filename)
	fp = open(filename,'rb')
	return load(fp)

def ensure(f):
	"""
	CREATE DIRECTORY IF IT DOES NOT EXIST
	"""
	if not path.exists(f):
		makedirs(f)

def lnPLUS(x,y):
	"""
	FUNCTION TO ADD IN LOGARITHMIC SPACE
	"""

	output = 0.0
	if x>y:
		output = x+log(1.0+exp(y-x))
	else:
		output = y+log(1.0+exp(x-y))

	return output

def OddsFromBayes(matrix, Nhyp):
	"""
	CALCULATE THE TOTAL ODDS FROM DATA STRUCTURE
	"""

	TotalOdds = -float_info.max
	# ONE SMALLER THAN THE WIDTH OF THE MATRIX FOR GR
	for j in range(shape(matrix)[1]-1):
		TotalOdds = lnPLUS(TotalOdds, sum(matrix[:, j+1]-matrix[:, 0]))
	TotalOdds -= log(pow(2,Nhyp)-1)

	return TotalOdds

def TigerCreateHistogram(list_tigerruns, axis, N=1, xlim=None,bins=50):
	"""
	CREATE TIGER HISTOGRAMS FROM A LIST OF TIGERRUNS
	"""
	for i in range(len(list_tigerruns)):
		if N == 1:
			# PLOT HISTOGRAMS
			axis.hist(list_tigerruns[i].odds(N), \
					facecolor="none", \
					label=list_tigerruns[i].latexlabel+r"\n$\mathrm{("+str(size(list_tigerruns[i].odds(N)))+r"\; sources)}$", \
					histtype="stepfilled", \
					linewidth=1, \
					range=xlim, \
					bins = bins, \
					hatch = hatch_default[i], \
					edgecolor = color_default[i], \
					density=True)
			# ADD X-AXIS LABEL
			axis.set_xlabel(r"$\ln{O_{GR}^{modGR}}$")
			axis.set_ylabel(r"$P(\ln{O_{GR}^{modGR}})$")
		elif N < list_tigerruns[i].nsources:
			# PLOT HISTOGRAMS
			axis.hist(list_tigerruns[i].odds(N), \
					facecolor="none", \
					label=list_tigerruns[i].latexlabel+r"\n$\mathrm{("+str(size(list_tigerruns[i].odds(N)))+r"\; catalogs)}$", \
					histtype="stepfilled", \
					linewidth=1, \
					range=xlim, \
					bins = bins, \
					hatch = hatch_default[i], \
					edgecolor = color_default[i], \
					density=True)
			# ADD X-AXIS LABEL
			axis.set_xlabel(r"$\ln{\mathcal{O}_{GR}^{modGR}}$")
			axis.set_ylabel(r"$P(\ln{\mathcal{O}_{GR}^{modGR}})$")

		# ADD LEGEND
		handles, labels = axis.get_legend_handles_labels()
		axis.legend(handles, labels, loc='upper left')

def TigerCalculateEfficiency(list_tigerruns, N=1, beta=[0.95], background=0):
	"""
	CALCULATE EFFICIENCY FROM A LIST OF TIGERRUNS
	"""
	efficiencies = []
	OddsBeta=[mquantiles(list_tigerruns[background].odds(N),prob=[b]) for b in beta]
	efficiencies = empty((len(list_tigerruns)-1,len(beta)))
	for i in range(len(list_tigerruns)):
		if N>list_tigerruns[i].nsources:
			stdout.write("... Warning: Not sufficient events (%s) to calculate the efficiency for %s sources. Writing zeros\n"%(list_tigerruns[i].nsources,N))
			if i < background:
				efficiencies[i,:] = 0.0
			else:
				efficiencies[i-1,:] = 0.0
			continue
		if i != background:
			tmp = list_tigerruns[i].odds(N)
			for j in range(len(OddsBeta)):
				msk = tmp>OddsBeta[j]
				nmsk = tmp<OddsBeta[j]
				nabovebeta=len(tmp[msk])
				ntotal=len(tmp)
				eff=float(nabovebeta)/float(ntotal)
				if i < background:
					efficiencies[i,j] = eff
				else:
					efficiencies[i-1,j] = eff
	return efficiencies

def TigerKSandPvalue(tigerrun1,tigerrun2, N=1):
	"""
	CALCULATE KS STATISTICS FROM A LIST OF TIGERRUNS
	"""
	if tigerrun1.nsources < N or tigerrun2.nsources < N:
		stdout.write("... Warning: Not sufficient events (%s,%s) to calculate KS-statistic for %s sources. Writing zeros\n"%(tigerrun1.nsources, tigerrun2.nsources,N))
		return [0.0,0.0]
	data1 = tigerrun1.odds(N)
	data2 = tigerrun2.odds(N)
	ksstat, pvalue = ks_2samp(data1,data2)
	return [ksstat, pvalue]


def TigerCreateSNRvsOdds(list_tigerruns, axis):
	"""
	CREATE SNR VS ODDS DIAGRAMS FROM A LIST OF TIGERRUNS
	"""
	markers = ['x','o']
	colors = ['blue','red']
	axis.set_xlabel(r"$\mathrm{SNR}$")
	axis.set_ylabel(r"$\ln{O_{\mathrm{GR}}^{\mathrm{modGR}}}$")

	for i in range(len(list_tigerruns)):
		axis.scatter(list_tigerruns[i].snr, list_tigerruns[i].odds(), \
				color=colors[i], \
				marker=markers[i], \
				label=list_tigerruns[i].latexlabel+r"\n$\mathrm{("+str(list_tigerruns[i].nsources)+r"\; sources)}$",\
				s=64)

	# SET AXIS LIMITS
	axis.set_xlim(5.0,50.0)

	# ADD LEGEND
	handles, labels = axis.get_legend_handles_labels()
	axis.legend(handles, labels, loc='best')

	return None

def TigerCreateCumFreq(tigerrun, axis):
	"""
	CREATE CUMULATIVE FREQUENCY DIAGRAMS FROM A LIST OF TIGERRUNS
	"""
	xlim = [5.0,100.0]
	snrrange = arange(xlim[0], xlim[1], 5.)
	nsnrsteps = len(snrrange)
	freqs = zeros((nsnrsteps,tigerrun.nsubhyp))

	for i in range(nsnrsteps):
		if snrrange[i] > min(tigerrun.snr):
			snrmask = tigerrun.snr<snrrange[i]
			bayesmasked = tigerrun.bayes[snrmask,:]
			maxarray = list(argmax(bayesmasked, axis=1))
			for j in range(tigerrun.nsubhyp):
					freqs[i,j] = maxarray.count(j)

	# CHECK IF TOTAL ADDS UP TO TOTAL NUMBER OF SOURCES
	if sum(freqs[-1,:]) != tigerrun.nsources:
		print(r"Warning, some sources maybe out of the SNR range and are not plotted: SNR \in [%.1f:%.1f], %i sources missing" %(xlim[0], xlim[1], tigerrun.nsources-sum(freqs[-1,:])))

	# PLOT CUMULATIVE HIGHEST BAYES SIGNAL VERSUS NOISE
	rcParams['legend.fontsize']=18 # MANUAL RESET FONTSIZE
	linestyles = cycle(['-', '--', '-.', ':'])
	markerstyles = cycle(['+','*','.','<','d', '^', 'x', 's'])
	colorstyles = cycle(["b", "g","r","c","m","y","k"])

	for i in range(tigerrun.nsubhyp):
			if tigerrun.subhyp[i].split("_")[-1] =="GR":
				axis.plot(snrrange, freqs[:,i], label=r'$\mathcal{H}_{\mathrm{GR}}$', color="black", linewidth=3.0)
			else:
				if tigerrun.engine == 'lalnest':
					curlabel = "$H_{"+''.join(tigerrun.subhyp[i].split("_")[-1].split('dphi'))+"}$"
				elif tigerrun.engine == 'lalinference':
					curlabel="$H_{"+''.join(tigerrun.subhyp[i].split("_")[-1].split('dchi'))+"}$"
				axis.plot(snrrange, freqs[:,i], \
						label=curlabel,\
						marker=markerstyles.next(), \
						linestyle='--', \
						linewidth=1.0, \
						color=colorstyles.next())

	axis.set_xlabel(r"$\mathrm{SNR}$")
	axis.set_ylabel(r"$\mathrm{Cumulative\;frequency\;(lines)}$")
	axis.set_xlim(xlim)

	# ADD LEGEND
	handles, labels = axis.get_legend_handles_labels()
	axis.legend(handles, labels, loc='best',ncol=3)

	# DRAW HISTOGRAM
	ax2 = axis.twinx()
	ax2.hist(tigerrun.snr, range=xlim, bins=nsnrsteps, alpha=0.25, color="k")
	ax2.grid(False)
	ax2.set_ylabel(r'$\mathrm{\#\;of\;sources\;per\;SNR\;bin\;(hist)}$')
	ax2.set_xlim(xlim)

def TigerCreateCumBayes(tigerrun,axis,nthcat,nsourcespercat):
	"""
	CREATE CUMULATIVE BAYES FACTOR/ODDS DIAGRAMS FROM A LIST OF TIGERRUNS (FOR
	EACH CATALOG)
	"""

	# COUNT CUMULATIVE FREQUENCY OF HIGHEST BAYES (SIGNAL VS. NOISE)
	snrcat = tigerrun.snr[nthcat*nsourcespercat:(nthcat+1)*nsourcespercat]
	bayescat = tigerrun.bayes[nthcat*nsourcespercat:(nthcat+1)*nsourcespercat,:]
	sortorder = argsort(snrcat)
	snrcatsorted = sort(snrcat)
	bayescatsorted = bayescat[sortorder,:]
	cumodds = array([OddsFromBayes(bayescatsorted[:i+1,:],len(tigerrun.testcoef)) for i in range(nsourcespercat)])
	cumbayes = array([sum(subtract(bayescatsorted[:i+1,1:].T,bayescatsorted[:i+1,0]).T,axis=0) for i in range(nsourcespercat)])

	markerstyles = cycle(['+','*','.','<','d', '^', 'x', 's'])

	axis.set_xlabel(r"$\mathrm{SNR}$")
	axis.set_ylabel(r"$\ln{{}^{(cat)}B^{i_1 i_2 \ldots i_k}_{\mathrm{GR}}}$")

	# PLOT BAYESFACTORS
	for i in range(tigerrun.nsubhyp-1):

		if tigerrun.engine == 'lalnest':
			curlabel = "$H_{"+''.join(tigerrun.subhyp[i+1].split("_")[-1].split('dphi'))+"}$"
		elif tigerrun.engine == 'lalinference':
			curlabel="$H_{"+''.join(tigerrun.subhyp[i+1].split("_")[-1].split('dchi'))+"}$"

		axis.plot(snrcatsorted, cumbayes[:,i],\
				label=curlabel, \
				marker=markerstyles.next(), \
				linestyle='--', \
				linewidth=1.0)

	# PLOT ODDS RATIO
	axis.plot(snrcatsorted, cumodds, label=r"$\ln{O_{GR}^{modGR}}$", linewidth=2.0, color="black")

	axis.set_xlim(5.0,50.0)
	# ADD LEGEND
	handles, labels = axis.get_legend_handles_labels()
	axis.legend(handles, labels, loc='best',ncol=3)

def TigerCreateExampleConfigFile(filename):
	"""
	CREATE EXAMPLE CONFIGURATION FILE (TO BE USED WITH -G OPTION)
	"""
	basename = path.basename(filename)
	dirname = path.dirname(filename)
	if dirname == '':
		dirname = './'

	# CREATE DIRECTORY IF IT DOES NOT EXIST
	ensure(dirname)

	# OPEN FILE
	outfile = open(path.join(dirname,basename),'w')

	examplecfgtxt="""\
###############################################################################
#
# Configuration file for tiger_postproc.py
#
###############################################################################

#------------------------------------------------------------------------------
[run]
#------------------------------------------------------------------------------

# Engine for TIGER runs: lalnest or lalinference
engine=engine1,engine2

# Subhypotheses: Notation depends on engine.
#		lalnest: dphi1,dphi2,...
#		lalinference: dchi1,dchi2,...
hypotheses=dphi1,dphi2,dphi3

# Run identifier (used in filenames for plots)
runid=name_for_run

# Local destination folder for postproc files: created if non-existing
localdest=/home/user/public_html/tiger/inj_tt4spin_tt4all_rec_tf2spin/

#------------------------------------------------------------------------------
[fetch]
#------------------------------------------------------------------------------

# Fetch data from (one per type of runs)
#		- source : directly from output of engine
#		- pickle : saved state
type=source,source

# Location files.(one per type of runs)
# 	source: ascii file two columns: cluster directory
# 		cluster: (atlas, nemo, cit, lho, llo, or local)
# 	pickle: location to python pickle file
locs=/path/to/the/location/file.txt,/path/to/the/location/file.txt

# Labels (one for each type of runs)
labels=tt4spin,tt4all

# Latex labels for different runs for plot legends (one for each run type).
# Currently not working with commas in the latex itself!
latexlabels=$\\mathrm{TaylorT4+spin}$,$\\mathrm{TaylorT4+all}$

# Seed for shuffling the sources. seed=0 means no shuffling
# NB: NOT WORKING!
seed=0

#------------------------------------------------------------------------------
[cut]
#------------------------------------------------------------------------------

# SNR threshold
snrthres=8.0

# B^GR_noise cutoff
bayesgrthres=32.0

#------------------------------------------------------------------------------
[plot]
#------------------------------------------------------------------------------

# Plot histograms (see [hist] for plotting options)
hist=true

# Plot SNR vs odds scatter plot (see [snrvsodds] for plotting options)
snrvsodds=true

# Plot cumulative frequency of a subhypothesis having the highest Bayes factor
# (see [cumfreq] for plotting options)
cumfreq=true

# Plot cumulative Bayes factor, per subhypothesis, and the odds ratio, as a
# function of SNR, per catalog. (see [cumbayes] for plotting options).
cumbayes=false

# File output extensions: png required for html page
# Other supported formats: pdf, ps, eps, svg
types=png,pdf

#------------------------------------------------------------------------------
[hist]
#------------------------------------------------------------------------------

# Number of bins in the histograms
bins=50

# Number of sources per catalog for the histograms
nsources=1,15

# Limits on the x-axis. Use xlim=0,0 for automatic determination (range will be
# the same for ALL sizes of the catalog).
xlim=0,0

#------------------------------------------------------------------------------
[snrvsodds]
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
[cumbayes]
#------------------------------------------------------------------------------

"""

	outfile.write(examplecfgtxt)

def TigerPreProcess(preprocessfile):
	"""
	TIGER PREPROCESSOR FUNCTION
	NB: UNDER DEVELOPMENT
	"""
	if access(preprocessfile, R_OK):
		preprocessfp = open(preprocessfile,'r')
		args = preprocessfp.read().split('\t')
		subhyp= args[3][1:-1].replace("'","").replace(" ","").split(',')
		tgloc = TigerRun(args[0].lower(),args[1],args[2],subhyp)
		tgloc.searchsources()
		tgloc.pullbayes()
		tgloc.pullposteriors()
		tgloc.savetopickle('tigerloc_preprocess.pickle')
		exit(0)
	else:
			exit('Preprocess file cannot be found - abort')

###############################################################################
#
# START THE MAIN FUNCTION
#
###############################################################################

if __name__ == "__main__":
	# START THE MAIN FUNCTION IF RUN AS A SCRIPT. OTHERWISE, JUST LOADS THE CLASS
	# AND FUNCTION DEFINITIONS
	exit(main())

#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       cbcBayesPostProc.py
#
#       Copyright 2010
#       Benjamin Farr <bfarr@u.northwestern.edu>,
#       Will M. Farr <will.farr@ligo.org>,
#       2014
#       Salvatore Vitale <salvatore.vitale@ligo.mit.edu>
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import table
from glue.ligolw import utils
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as pp
import numpy as np
import optparse
import os
import os.path
import scipy.stats as ss
import string

class LIGOLWContentHandlerExtractSimBurstTable(ligolw.LIGOLWContentHandler):
    def __init__(self,document):
      ligolw.LIGOLWContentHandler.__init__(self,document)
      self.tabname=lsctables.SimBurstTable.tableName
      self.intable=False
      self.tableElementName=''
    def startElement(self,name,attrs):
      if attrs.has_key('Name') and attrs['Name']==self.tabname:
        self.tableElementName=name
        # Got the right table, let's see if it's the right event
        ligolw.LIGOLWContentHandler.startElement(self,name,attrs)
        self.intable=True
      elif self.intable: # We are in the correct table
        ligolw.LIGOLWContentHandler.startElement(self,name,attrs)
    def endElement(self,name):
      if self.intable: ligolw.LIGOLWContentHandler.endElement(self,name)
      if self.intable and name==self.tableElementName: self.intable=False

lsctables.use_in(LIGOLWContentHandlerExtractSimBurstTable)

posterior_name_to_sim_burst_extractor = {
    'frequency' : lambda sb: sb.frequency,
    'duration' : lambda sb: sb.duration,
    'quality' : lambda sb: sb.q,
    'hrss' : lambda sb: sb.hrss,
    'psi' : lambda sb: sb.psi,
    'time' : lambda sb: sb.time_geocent_gps + 1e-9*sb.time_geocent_gps_ns,
    'ra' : lambda sb: sb.ra,
    'dec' : lambda sb: sb.dec,
    'polar_angle':lambda sb:sb.pol_ellipse_angle,
    'polar_eccentricity':lambda sb:sb.pol_ellipse_e,
    'alpha':lambda sb:sb.pol_ellipse_angle
}


posterior_name_to_latex_name = {
    'frequency' : r'$f$',
    'quality' : r'$q$',
    'duration' : r'$\tau$',
    'hrss' : r'$hrss$',
    'time' : r'$t$',
    'ra' : r'$\alpha$',
    'dec' : r'$\delta$',
    'psi' : r'$\psi$',
    'polar_angle' : r'$\xi_a$',
    'polar_eccentricity' : r'$\xi_e$',
    'alpha' : r'$\zeta_{orb}$',


}

def fractional_rank(x, xs):
    """Returns the fraction of samples, ``xs``, that fall below the value
    ``x``.

    """
    nbelow = np.sum(xs < x)

    return float(nbelow)/float(xs.shape[0])

def pp_plot(ps, title=None, outfile=None):
    """Generates a p-p plot for the given ps.  

    :param ps: The p-values whose cumulative distribution is to be
      plotted.

    :param title: An (optional) title for the plot.

    :param outfile: An (optional) basename for the plot file; both
      basename.png and basename.pdf will be created.

    """
    
    ps = np.atleast_1d(ps)
    ps = np.sort(ps)
    ys = np.zeros(ps.shape[0]+2)
    ys[:-1] = np.linspace(0, 1, ps.shape[0]+1)
    ys[-1] = 1.0
    xs = np.zeros(ps.shape[0]+2)
    xs[-1] = 1.0
    xs[1:-1] = ps

    pp.figure(figsize=(6,6), dpi=100)

    for i in range(10):
        syn_ps = np.random.uniform(size=ps.shape[0])
        syn_xs = np.zeros(ps.shape[0]+2)
        syn_xs[-1] = 1.0
        syn_xs[1:-1] = np.sort(syn_ps)
        pp.plot(syn_xs, ys, '-', color='0.9')

    pp.plot(xs, ys, '-k')
    pp.plot(ys, ys, '--k')

    pp.xlabel(r'$p$')
    pp.ylabel(r'$P(p)$')
    if title is not None:
        pp.title(title)

    if outfile is not None:
        pp.savefig(outfile + '.png')
        pp.savefig(outfile + '.pdf')

def pp_kstest_pvalue(ps):
    """Returns the K-S p-value for the test of the given ``ps`` against a
    uniform distribution on [0,1].

    """
    stat, p = ss.kstest(ps, lambda x: x)

    return p

def read_posterior_samples(f):
    """Returns a named numpy array of the posterior samples in the file
    ``f``.

    """
    with open(f, 'r') as inp:
        header = inp.readline().split()
        dtype = np.dtype([(n, np.float) for n in header])
        data = np.loadtxt(inp, dtype=dtype)

    return data

def output_html(outdir, ks_pvalues, injnum,skypp=False):
    """Outputs the HTML page summarizing the results.

    """
    table_row_template = string.Template("""<tr> <td> ${name} </td>
    <td> ${pvalue} </td>
    <td> <img src="${name}.png" alt="${name} p-p plot" width="300" height="300" /> </td> <td> <a href="${name}.png">PNG</a> <a href="${name}.pdf">PDF</a> <a href="${name}-ps.dat">p-values</a> </td> </tr>

    """)

    html_template = string.Template("""<!DOCTYPE html> 
    <html>
    <head>
    <title> LALInference P-P Plots </title>
    </head>

    <body>

	<p>This page was generated with the output of ${injnum} simulations.</p>
	${linkstr}
	<br>
    <table border="1"> 
    <tr>
    <th> Parameter </th> <th> K-S p-value </th> <th> p-p Plot </th> <th> Links </th>
    </tr>

    ${tablerows}
    ${skypp}
    </table>

    </body>
    </html>

    """)

    # If this script is run with lalinference_pp_pipe then the following directory structure should exist
    links="<ul>"
    if os.path.isdir(os.path.join(outdir,'prior')):
        links+='<li><a href="prior/">Prior Samples used in this test</a>'
    if os.path.isdir(os.path.join(outdir,'injections')):
        links+='<li><a href="injections/">Posteriors for each injection</a>'
    links+='</ul>'

    tablerows = []
    for par, pv in ks_pvalues.items():
        tablerows.append(table_row_template.substitute(name=par, pvalue=str(pv)))
    tablerows = '\n'.join(tablerows)

    if skypp is False:
      skysub=''
    else:
      skysub='<tr> <td> sky </td> <td> See Plot </td> <td> <img src="sky_p-p.png" alt="sky p-p plot" width="300" height="300" /> </td> <td> <a href="sky_p-p.png">PNG</a> <a href="sky_p-p.pdf">PDF</a></td>  </tr>'

    html = html_template.substitute(tablerows=tablerows, injnum=str(injnum), linkstr=links,skypp=skysub)

    with open(os.path.join(outdir, 'index.html'), 'w') as out:
        out.write(html)

if __name__ == '__main__':
    USAGE='''%prog [options] posfile1.dat posfile2.dat ...
	            Generate PP analysis for a set of injections. posfiles must be in same order as injections.'''
    parser = optparse.OptionParser(USAGE)
    parser.add_option('--injXML', action='store', type='string', dest='injxml', 
                      help='sim_burst XML file for injections')
    parser.add_option('--outdir', action='store', type='string',
                      help='output directory')

    parser.add_option('--postsamples', action='store', type='string', 
                      default='posterior_samples.dat', 
                      help='filename for posterior samples files')

    parser.add_option('--par', action='append', default=[], type='string', 
                      help='parameter names for the p-p plot')
    parser.add_option('--skyPPfolder', action='store',dest='skypp',type='string',default=None,help='Path to folder containing png/pdf with 2D skyarea PP plots')

    (options, args) = parser.parse_args()

    injs = table.get_table(utils.load_filename(options.injxml,contenthandler=LIGOLWContentHandlerExtractSimBurstTable),
                           lsctables.SimBurstTable.tableName)

    if options.par == []:
        parameters = ['frequency', 'quality', 'hrss', 'ra', 'dec', 'psi', 'time', 'alpha','polar_eccentricity']
    else:
        parameters = options.par

    try:
        os.mkdir(options.outdir)
    except:
        pass

    pvalues = { }
    posfiles=args
    Ninj=0
    for index,posfile in enumerate(posfiles):
	    try:
	      psamples = read_posterior_samples(posfile)
	      #index = int(element)
	      true_params = injs[index]
	      Ninj+=1
	    except:
	      # Couldn't read the posterior samples or the XML.
              print "not found %s \n"%posfile
	      continue

	    for par in parameters:
	      try:
		samples = psamples[par]
		true_value = posterior_name_to_sim_burst_extractor[par](true_params)
		p = fractional_rank(true_value, samples)

		try:
		  pvalues[par].append(p)
		except:
		  pvalues[par] = [p]
	      except:
		# Couldn't read samples for parameter or injection
		continue

    # Generate plots, K-S tests
    ks_pvalues = {}
    for par, ps in pvalues.items():
        pp_plot(ps, title=posterior_name_to_latex_name[par], outfile=os.path.join(options.outdir, par))
        pp.clf()
        ks_pvalues[par] = pp_kstest_pvalue(ps)
        np.savetxt(os.path.join(options.outdir, par + '-ps.dat'), np.reshape(ps, (-1, 1)))

    skypp=False
    if options.skypp is not None:
      found=0
      if os.path.isdir(options.skypp):
        for i in ['png','pdf']:
          if os.path.isfile(os.path.join(options.skypp,'p-p.%s'%i)):
            inf=os.path.join(os.path.realpath(options.skypp),'p-p.%s'%i)
            outf=os.path.join(options.outdir,'sky_p-p.%s'%i)
            os.system('cp %s %s'%(inf,outf))
            found+=1
          else:
            print "could not find %s\n"%os.path.join(options.skypp,'p-p.%s'%i)
      else:
        print "skyPPfolder %s doesn't seem to be a valid folder or cannot be read. Skipping skyPP plot\n"%os.path.realpath(options.skypp)

      if found>0:
        skypp=True
        
    output_html(options.outdir, ks_pvalues, Ninj,skypp= skypp )


#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       cbcBayesPPAnalysis.py
#
#       Copyright 2013
#       Benjamin Farr <bfarr@u.northwestern.edu>,
#       Will M. Farr <will.farr@ligo.org>,
#       SalvatoreVitale <salvatore.vitale@ligo.org>
#       Vivien Raymond <vivien.raymond@ligo.org>
#
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
from pylal.SimInspiralUtils import ExtractSimInspiralTableLIGOLWContentHandler
lsctables.use_in(ExtractSimInspiralTableLIGOLWContentHandler)
from pylal import bayespputils as bppu

posterior_name_to_latex_name = {
    'm1' : r'$m_1$',
    'm2' : r'$m_2$',
    'eta' : r'$\eta$',
    'q' : r'$q$',
    'mc' : r'$\mathcal{M}$',
    'distance' : r'$d$',
    'time' : r'$t$',
    'ra' : r'$\alpha$',
    'dec' : r'$\delta$',
    'phi_orb' : r'$\phi_\mathrm{orb}$',
    'psi' : r'$\psi$',
    'iota' : r'$\iota$',
    'a1' : r'$a_1$',
    'a2' : r'$a_2$',
    'theta1' : r'$\theta_1$',
    'theta2' : r'$\theta_2$',
    'phi1' : r'$\phi_1$',
    'phi2' : r'$\phi_2$',
    'phi12':r'$\phi_{12}$',
    'phi_jl':r'$\phi_{jl}$',
    'theta_jn':r'$\theta_{jn}$',
    'tilt1':r'$\tau_1$',
    'tilt2':r'$\tau_2$',
    'dchi0':r'$d\chi_0$',
    'dchi1':r'$d\chi_1$',
    'dchi2':r'$d\chi_2$',
    'dchi3':r'$d\chi_3$',
    'dchi4':r'$d\chi_4$',
    'dchi5':r'$d\chi_5$',
    'dchi5l':r'$d\chi_5^l$',
    'dchi6':r'$d\chi_6$',
    'dchi6l':r'$d\chi_6^l$',
    'dchi7':r'$d\chi_7$',
    'dsigma1':r'$d\sigma_1$',
    'dsigma2':r'$d\sigma_2$',
    'dsigma3':r'$d\sigma_3$',
    'dsigma4':r'$d\sigma_4$',
    'dalpha1':r'$d\alpha_1$',
    'dalpha2':r'$d\alpha_2$',
    'dalpha3':r'$d\alpha_3$',
    'dalpha4':r'$d\alpha_4$',
    'dalpha5':r'$d\alpha_5$',
    'dbeta1':r'$d\beta_1$',
    'dbeta2':r'$d\beta_2$',
    'dbeta3':r'$d\beta_3$',
    'dbeta4':r'$d\beta_4$'
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

def read_posterior_samples(f,injrow):
    """Returns a bppu posterior sample object
    """
    peparser=bppu.PEOutputParser('common')
    commonResultsObj=peparser.parse(open(f,'r'))
    data = bppu.Posterior(commonResultsObj,SimInspiralTableEntry=injrow,injFref=100.0)
    # add tilts, comp masses, tidal...
    try:
      data.extend_posterior()
    except Exception,e:
      pass 
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
                      help='sim_inspiral XML file for injections')
    parser.add_option('--outdir', action='store', type='string',
                      help='output directory')

    parser.add_option('--postsamples', action='store', type='string', 
                      default='posterior_samples.dat', 
                      help='filename for posterior samples files')

    parser.add_option('--par', action='append', default=[], type='string', 
                      help='parameter names for the p-p plot')
    parser.add_option('--skyPPfolder', action='store',dest='skypp',type='string',default=None,help='Path to folder containing png/pdf with 2D skyarea PP plots')

    (options, args) = parser.parse_args()

    injs = table.get_table(utils.load_filename(options.injxml,contenthandler=ExtractSimInspiralTableLIGOLWContentHandler),lsctables.SimInspiralTable.tableName)

    if options.par == []:
        parameters = ['m1', 'm2', 'mc', 'eta', 'q',  'theta_jn', 'a1', 'a2', 'tilt1', 'tilt2', 'phi12', 'phi_jl', 'ra', 'dec', 'distance', 'time', 'phi_orb', 'psi']
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
	      psamples = read_posterior_samples(posfile,injs[index])
	      Ninj+=1
	    except:
	      # Couldn't read the posterior samples or the XML.
	      continue

	    for par in parameters:
	      try:
		samples = psamples[par].samples
                true_value=psamples[par].injval
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
        print "Trying to create the plot for",par,"."
        try:
          pp_plot(ps, title=posterior_name_to_latex_name[par], outfile=os.path.join(options.outdir, par))
          pp.clf()
          ks_pvalues[par] = pp_kstest_pvalue(ps)
          np.savetxt(os.path.join(options.outdir, par + '-ps.dat'), np.reshape(ps, (-1, 1)))
        except:
          print "Could not create the plot for",par,"!!!"

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

    output_html(options.outdir, ks_pvalues, Ninj ,skypp=skypp )


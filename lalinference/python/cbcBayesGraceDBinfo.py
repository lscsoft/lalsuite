# -*- coding: utf-8 -*-
#
#       cbcBayesGraceDBinfo
#
#       Copyright 2015
#       Salvatore Vitale <salvatore.vitale@ligo.org>
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

#===============================================================================
# Preamble
#===============================================================================

#standard library imports
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")

from lalinference import bayespputils as bppu
from lalinference import git_version

__author__="Salvatore Vitale <salvatore.vitale@ligo.org>"
__version__= "git id %s"%git_version.id
__date__= git_version.date

USAGE='''%prog --gid graceDBid --samples posterior_samples.dat

Upload to the pe section of the graceDB page gid information about a lalinference postprocessing run.
For the moment this is maximum a posteriori and stdev for the parameters in the string pars.

'''

def cbcBayesGraceDBinfo(gid=None,samples=None,skymap=None,analysis='LALInference', bcifile=None,bsnfile=None,email=None,message=None,server="https://gracedb.ligo.org/api/"):

  if gid is None or (samples is None and skymap is None):
    print("Must provide both a graceDB id and a posterior samples file or skymap file\n")
    sys.exit(1)

  import ligo.gracedb.rest
  import os
  if server is not None:
    g=ligo.gracedb.rest.GraceDb(server)
  else:
    g=ligo.gracedb.rest.GraceDb()
  if samples is not None:
    samples=os.path.realpath(samples)
    if '.hdf' in samples  or '.h5' in samples:
      peparser = bppu.PEOutputParser('hdf5')
      commonResultsObj=peparser.parse(samples)
    else:
      peparser=bppu.PEOutputParser('common')
      commonResultsObj=peparser.parse(open(samples,'r'))

    try:
      pos = bppu.BurstPosterior(commonResultsObj)
      pars=['frequency','quality','hrss']
      units={'frequency':'[Hz]','quality':'','hrss':'','loghrss':''}
    except:
      pos = bppu.Posterior(commonResultsObj)
      pars=['mchirp','q','distance']
    strs=[]
    outstr='<table><tr><th colspan=2 align=center>%s PE summary</th></tr>'%analysis

    for i in pars:
      if i in pos.names:
        _,which=pos._posMap()
        if i=='hrss':
          outstr+='<tr><td align=left>%s %s</td>'%(i,units[i])
          outstr+='<td align=left>%.3e &plusmn; %.3e</td></tr>'%(pos[i].samples[which][0],pos[i].stdev)
        else:
          outstr+='<tr><td align=left>%s %s</td>'%(i,units[i])
          outstr+='<td align=left>%.3f &plusmn; %.3f</td></tr>'%(pos[i].samples[which][0],pos[i].stdev)
    if bcifile is not None and os.path.isfile(bcifile):
      bci=np.loadtxt(bcifile)
    else: bci=None
    if bci is not None:
      outstr+='<tr><td align=left>logBCI</td>'
      outstr+='<td align=center>%.2f</td></tr>'%(bci)

    bsn=None
    if bsnfile is not None and os.path.isfile(bsnfile):
      bsn=np.loadtxt(bsnfile)
      bsn=bsn[0]
    else:
      try:
        import h5py
        with h5py.File(samples,'r') as h5grp:
          tmp=h5grp['lalinference']['lalinference_nest'].attrs
          bsn=tmp['log_bayes_factor']
      except Exception as e:
        print("Could not obtain BNS\n")
        print(e)

    if bsn is not None:
      outstr+='<tr><td align=left>logBSN</td>'
      outstr+='<td align=center>%.2f</td></tr>'%(bsn)
    outstr+='</table>'

    if email is not None and bci is not None:
      import os
      import smtplib
      address=email.split(',')
      SERVER="localhost"
      USER=os.environ['USER']
      import socket
      HOST=socket.getfqdn()#socket.gethostbyaddr(socket.gethostname())[0]
      pref=""
      if bci>3 and bci<6:
        pref="A promising"
      elif bci>6 and bci<10:
        pref="A very interesting"
      elif bci>10:
        pref="A SPECTACULAR"
      FROM="salvatore.vitale@"+HOST
      SUBJECT="%s LIB result page is ready at "%pref+HOST+" for graceID %s!"%(gid)
      TEXT="LIB run for graceID %s is done on "%gid+HOST+".\nThe BCI is %lf\n"%bci
      if bci>10:
        TEXT+="RUN!!!!!!!!!!\n"
      message="From: %s\nTo: %s\nSubject: %s\n\n%s"%(FROM,', '.join(address),SUBJECT,TEXT)
      try:
        import os
        os.system('echo "%s" | mail -s "%s" "%s"'%(TEXT,SUBJECT,', '.join(address)))
        server=smtplib.SMTP(SERVER)
        server.sendmail(FROM,address,message)
        server.quit()
      except:
        print("Cound not send email\n")

    g.writeLog(gid,outstr,filename=None,tagname='pe')
  elif skymap is not None:
    if bcifile is not None and os.path.isfile(bcifile):
      bci=np.loadtxt(bcifile)
    else: bci=None
    if bsnfile is not None and os.path.isfile(bsnfile):
      bsn=np.loadtxt(bsnfile)
    else: bsn=None
    tag=['sky_loc']
    """
    if bci is not None and bsn is not None:
      if bsn>5. and bci>2.:
        tag.append('lvem')
    """
    g.writeLog(gid,message,filename=skymap,tagname=tag)


if __name__=='__main__':

    from optparse import OptionParser
    parser=OptionParser(USAGE)
    parser.add_option("-g","--gid", dest="gid",help="GraceDB id", metavar="G123456",default=None)
    parser.add_option("-s","--samples",dest="samples",help="posterior_samples.hdf5/dat",default=None)
    parser.add_option("--analysis",help="Prefix to use for the graceDB entries. Should be the name of the analysis (default LALInference)",default='LALInference')
    parser.add_option("--bci",dest="bci",help="coherence test file: bci.dat",default=None)
    parser.add_option("--bsn",dest="bsn",help="evidence file: bsn.dat [Deprecated, now is read from samples.hdf5]",default=None)
    parser.add_option("--skymap",dest="skymap",help="FITS file skymap",default=None)
    parser.add_option("--message",dest="message",type='str', help="Message to go with skymap uplaod",default=None)
    parser.add_option('--email',dest='email',help="Will email when run is done.",default=None)
    parser.add_option('--server',dest='server',help="GraceDB server to interact with.",default="https://gracedb.ligo.org/api/")

    (opts,args)=parser.parse_args()
    if opts.gid is None:
      print("Must provide a graceDB id with --gid/-g ")
      sys.exit(1)
    if opts.samples is None and opts.skymap is None:
      print("Must provide lalinference posterior samples with --samples/-s or FITS skymap with --skymap ")
      sys.exit(1)
    cbcBayesGraceDBinfo(opts.gid, opts.samples,analysis=opts.analysis,bcifile=opts.bci,bsnfile=opts.bsn,email=opts.email,skymap=opts.skymap,message=opts.message,server=opts.server)

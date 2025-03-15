# -*- coding: utf-8 -*-
"""
Michalis Agathos (2014)

This script provides a set of good gps times for injecting a GW signal
into real multiple detector data.
The input consists of a set of science segment and veto segment ASCII files
having the 4-column format:

segid   gpsstart   gpsend   duration

The output consists of an ASCII file containing a list of gps injection times.
If the timeslides functionality is activated, an additional ASCII file is
output, containing timeslide information (in seconds) for each IFO, with
respect to the nominal times provided in the first file.
"""

__author__ = "Michalis Agathos"
__maintainer__ = "Michalis Agathos"
__email__ = "michalis.agathos@ligo.org"
__version__ = "1.0.0"
__status__ = "Production"

import matplotlib as mpl
import argparse
import os
import sys
import random as rd
mpl.use('Agg')
from pylab import *
py_version = sys.version_info[:2]
np_version = np.__version__

colordict={'H1':'r', 'L1':'g', 'V1':'b', 'H1L1':'c', 'L1H1':'c', 'H1V1':'m', 'V1H1':'m', 'L1V1':'y', 'V1L1':'y', 'H1L1V1':'k'}

'''Exclude segments that contain the following gpstimes'''
excludetimes = [966951349, 967054240, 968021010]

###########################
#
#  CLASS DEFINITIONS
#
###########################


class IFO:
  '''An interferometer. Can also refer to multiple interferometers to host doubles, triples etc.'''
  def __init__(self, name, segments, vetoes, minlen=1):
    self._name = name
    self._minlen = minlen

    if type(segments) is str:
      print("Reading segments for " + self._name + "")
      self._allsegments = self.readSegments(segments)
    elif type(segments) is ndarray:
      self._allsegments = segmentData(list(segments))
    elif  isinstance(segments, segmentData):
      self._allsegments = segments
    else:
      print("Cannot recognize segments!")
      return -1
    self._segments = self._allsegments.fits(self._minlen)

    if type(vetoes) is str:
      print("Reading veto segments for " + self._name + "")
      self._vetoes = self.readSegments(vetoes)
    elif type(vetoes) is ndarray:
      self._vetoes = segmentData(list(vetoes))
    elif isinstance(vetoes, segmentData):
      self._vetoes = vetoes
    else:
      print("Cannot recognize veto segments!")
      return -1

    self.setUnvetoed()
    print("Number of unvetoed segments that fit " + str(self._minlen) + ": " + str(self._unvetoed.length()))

  def setUnvetoed(self):
    '''This is getting a list of unvetoed segments that fit the minimum length'''
    self._unvetoed = self._segments.getUnvetoed(self._vetoes).fits(self._minlen)
    #self._unvetoed = self._segments.fits(self._minlen).getUnvetoed(self._vetoes).fits(self._minlen)

  def readSegments(self, segfname):
    print("Reading " + segfname)
    segdata = genfromtxt(segfname, dtype=int)
    return segmentData(list(segdata))

  def printUnvetoedToFile(self, outname):
    self._unvetoed.printToFile(outname)

  def getTrigTimes(self, whereInj='middle', interval=None, rmargin=2, lmargin=2, n=None, outfile=None):
    '''Returns a list of gps times on which injections can be made'''
    trigtimes = []
    if interval is None:
      l = self._minlen
    else:
      l = interval
    for seg in self._unvetoed._seglist:
      t = int(seg[1])
      if whereInj == 'often':
        '''Inject every l seconds until you reach rmargin seconds from end of segment.'''
        t += lmargin
        while t + rmargin <= seg[2]:
          if t >= seg[1] + lmargin:
            trigtimes.append(t)
          t += l
      elif whereInj == 'middle':
        '''Inject in the middle of minlen-long intervals.'''
        while t + self._minlen <= seg[2]:
          if (t + self._minlen//2 >= seg[1] + lmargin) and (t + self._minlen//2 + rmargin <= seg[2]):
            trigtimes.append(t + self._minlen//2)
          t += self._minlen
      else:
        print("Invalid whereInj argument. Output times will be empty.")

    if n is not None:
      trigtimes = trigtimes[:n]
    if outfile is None:
      return array(trigtimes)
    else:
      savetxt(outfile, array(trigtimes), fmt='%i')

  def plotCumulativeDurations(self, outfile, maxdur=None):
    print("Plotting segment lengths distribution to file " + outfile)
    fig = figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel("segment length")
    ax.set_ylabel("Cumulative # segments")
    sdur = array(self._segments._seglist)[:,-1]
    udur = array(self._unvetoed._seglist)[:,-1]
    if maxdur is not None:
      sdur = sdur[where(sdur <= maxdur)[0]]
      udur = udur[where(udur <= maxdur)[0]]
    ax.set_ylim((0,max([len(udur),len(sdur)])))
    ax.set_title("Segment length distribution for " + self._name)
    ax.hist(sdur, bins=sort(sdur), cumulative=True, histtype='stepfilled', alpha=0.3, label="All segments (" + str(self._segments.length()) + ")")
    ax.hist(udur, bins=sort(udur), cumulative=True, histtype='stepfilled', alpha=0.3, label="Unvetoed segments (" + str(self._unvetoed.length()) + ")")
    ax.legend(loc=4)
    fig.savefig(outfile)

  def plotSegments(self, outfile=None, lenperline=200000, segcolorlist=['b','r','g'], title=None):
    '''Plot segments in rows of lenperline seconds'''
    start = min(self._allsegments._start, self._vetoes._start)
    end = max(self._allsegments._end, self._vetoes._end)
    n = ceil((end-start)/lenperline)
    plot_params = {'figure.figsize': [fig_width,fig_height*2]}
    rcParams.update(plot_params)
    fig = figure()
    for k, segs in zip([0,1,2],[self._segments, self._vetoes, self._unvetoed]):
      s = array(segs._seglist)
      t = sort(hstack((s[:,1],s[:,1],s[:,2],s[:,2])))
      y = array(len(s)*[0,1,1,0])
      for i in arange(n):
        ax=fig.add_subplot(n,1,i+1)
        ax.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
        ax.tick_params(axis='x', which='both', top='on', bottom='off', labelbottom='off')
        ax.set_ylim(0,2)
        ax.xaxis.set_major_formatter(FormatStrFormatter('%0.0f'))
        t0 = start + i*lenperline
        t1 = t0 + lenperline
        ax.set_ylabel(int(t0), verticalalignment='top', rotation=45)
        ax.set_xlim(t0, t1)
        ax.fill_between(t, 0*y + 2*(k//2), y + 2*(1-y)*(k//2), alpha=0.3, color=segcolorlist[k])
    if title is not None:
      fig.suptitle(title)
    if outfile is None:
      return fig
    else:
      print('Plotting segments to ', outfile)
      fig.savefig(outfile)
      plot_params = {'figure.figsize': [fig_width,fig_height]}
      rcParams.update(plot_params)


#  def distToVeto(self, time):
#    '''Returns the distance to the closest veto segment in seconds'''
#    dist = self._unvetoed

class segment:
  '''A time segment in an IFO'''

  def __init__(self, data, gpsstart=None, gpsend=None):
    '''Creates a segment'''
    if gpsstart is not None and gpsend is not None:
      self._id = data[0]
      self._start = gpsstart
      self._end = gpsend
      self._len = max(gpsend-gpsstart, 0)
    elif gpsstart is None and gpsend is None:
      '''Overloading to generate segment from an array (id, start, end, length)'''
      if len(data) != 4:
        print("Segment data doesn't have the correct length!")
        return -1
      self._id = data[0]
      self._start = data[1]
      self._end = data[2]
      self._len = data[3]
    else:
      print("Wrong arguments to create segment!")
      return -1
#    if self._len != (self._end - self._start):
#      print("Error in segment data: inconsistent length! " + str(self._len) + " " + str(self._end - self._start))

  def intersectWithList(self, id0, other, sortedlist=True, startid=None):
    '''Intersect segment with list of (non-overlapping) segments'''
    newlist = []
    id = id0
    minid = startid
    '''Cut list short from the left'''
    for oseg in other[startid:]:
      if segment(oseg)._start > self._end and sortedlist:
        minid = max((0,segment(oseg)._id - 1))
        break
      newseg = self.intersectWith(segment(oseg))
      if newseg._start < newseg._end:
        newseg._id = id
        id += 1
        newlist.append(newseg.toArray())
    if startid is None:
      return newlist
    else:
      return (newlist, minid)

  def hasTime(self, time):
    '''Checks if time is in (open) segment'''
    if (time < self._start or time > self._end):
      dist = min(self._start - time, time - self._end)
    else:
      dist = 0
    return dist

  def intersectWith(self, other):
    '''Intersects with another segment'''
    newstart = max(self._start, other._start)
    newend = min(self._end, other._end)
    newseg = segment([self._id], newstart, newend)
    return newseg

  def toArray(self):
    a = array([self._id, self._start, self._end, self._len])
    return a

class segmentData:
  '''Data that holds segments for one or more IFOs.'''

  def __init__(self, seglist, gpsstart=None, gpsend=None):
    if gpsstart is None:
      gpsstart = seglist[0][1]
    if gpsend is None:
      gpsend = seglist[-1][2]
    sortedlist = list(array(seglist)[argsort(array(seglist)[:,1])])
    self._seglist = sortedlist
    self._start = gpsstart
    self._end = gpsend

  def length(self):
    return len(self._seglist)

  def fits(self, minlen):
    '''Returns the list of segments that fit a minimum required length.'''
    data = array(self._seglist)
    fitted = data[where(data[:,3] >= minlen)[0]]
    fitted = vstack( (arange(len(fitted)),fitted[:, 1:].T) ).T
    return segmentData(list(fitted), self._start, self._end)

  def intersectSegments(self, other):
    '''Returns the intersection with another list of segments.'''
    newdata = []
    startfrom = 0
    id0 = self._seglist[0][0]
    for s in self._seglist:
      exclude = False
      seg = segment(s)
      for et in excludetimes:
        exclude = exclude or (seg.hasTime(et) == 0)
      if not exclude:
        id = id0 + len(newdata)
        nd, minid= seg.intersectWithList(id, other._seglist, startid=startfrom)
   #     nd= seg.intersectWithList(id, other._seglist, startid=None)
        newdata += nd
        startfrom = minid
    return segmentData(newdata)

  def notSegments(self, start, end):
    '''Get the complement of the (non-overlapping) segments in this data'''
    notlist = []
    c = 0
    times = array(self._seglist)[:,[1,2]]
    # print(shape(times))
    times = times[argsort(times[:,1])]
    t1 = start
    for i in (arange(len(times))):
      t2 = times[i,0]
      if t2 > t1 and t1 >= start and t2 <= end:
        newseg = array([c, t1, t2, t2-t1]) #FIXME: check c initial value
        notlist.append(newseg)
        c += 1
      t1 = times[i,1]
    if t1 < end:
      newseg = array([c, t1, end, end-t1])
      notlist.append(newseg)
    return segmentData(notlist, start, end)

  def getUnvetoed(self, vetodata):
    '''Returns a segmentData object with unvetoed segments'''
    noveto = vetodata.notSegments(self._start, self._end)
    unvetoed = self.intersectSegments(noveto)
    return unvetoed

  def unionSegments(self, other):
    '''Take the union of segment lists. This is used for combining veto segments of different detectors.'''
    '''Here I use AUB = not(notA^notB)'''
    start = min(self._start, other._start)
    end = max(self._end, other._end)
    united = self.notSegments(start, end).intersectSegments(other.notSegments(start, end)).notSegments(start, end)
    return united

  def hasTime(self, time):
    '''Checks whether time is in the segments. Gives 0 if true and the distance of time to the segment list otherwise.'''
    dist = max(abs(self._start - time), abs(self._end - time))
    for s in self._seglist:
      dist = min(dist, s.hasTime(time))
    return dist

  def timeIn(self, time):
    '''If t is in a segment returns the distance to the edges, otherwise returns False'''
    sl = array(self._seglist)
    ds = time - sl[:,1]
    de = time - sl[:,2]
    sgn = ds*de    #should be negative if segment contains t
    ids = where(sgn < 0)[0]
    if len(ids) > 1:
      print("Time " + time + " contained in more than one segment!")
    if len(ids) == 1:
      return [ds[ids[0]], -de[ids[0]]]
    else:
      return False


  def printToFile(self, outfile):
    print("Printing segment list to file " + outfile)
    savetxt(outfile, array(self._seglist), fmt='%i')

  def plotSegments(self, outfile=None, lenperline=200000, segcolor='b', title=None):
    '''Plot segments in rows of lenperline seconds'''
    s = array(self._seglist)
    t = sort(hstack((s[:,1],s[:,1],s[:,2],s[:,2])))
    y = len(s)*[0,1,1,0]
    n = ceil((t[-1]-t[0])/lenperline)
    fig = figure()
    for i in arange(n):
      ax=fig.add_subplot(n,1,i+1)
      ax.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
      ax.tick_params(axis='x', which='both', top='on', bottom='off', labelbottom='off')
      ax.set_ylim(0,2)
      ax.xaxis.set_major_formatter(FormatStrFormatter('%0.0f'))
      t0 = t[0] + i*lenperline
      t1 = t0 + lenperline
      ax.set_ylabel(int(t0), verticalalignment='top', rotation=45)
      ax.set_xlim(t0, t1)
      ax.fill_between(t, 0*array(y), y, alpha=0.3, color=segcolor)
    if title is not None:
      fig.suptitle(title)
    if outfile is None:
      return fig
    else:
      print('Plotting segments to ', outfile)
      fig.savefig(outfile)


  def plotCumulativeDurations(self, outfile, maxdur=None):
    print("Plotting segment lengths distribution to file " + outfile)
    fig = figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel("segment length")
    ax.set_ylabel("# segments")
    dur = array(self._seglist)[:,-1]
    if maxdur is not None:
      dur = dur[where(dur <= maxdur)[0]]
    ax.hist(dur, bins=sort(dur), cumulative=True, histtype='stepfilled', alpha=0.3)
    fig.savefig(outfile)



###########################
#
#  OTHER FUNCTIONS
#
###########################


def ensure_dir(f):
    """
    CREATE FOLDER IF IT DOES NOT EXIST
    """
    if not os.path.exists(f):
        os.makedirs(f)


def getDoubles(ifo1, ifo2, unvetoed=False):
  '''
  Combines 2 interferometer objects into a new one, taking intersection of segments and union of vetoes
  Setting unvetoed=True will speed up the calculations by ignoring
  small/vetoed segments. Use with caution! (not good for plotting seglen distr)
  '''
  name1 = ifo1._name
  name2 = ifo2._name
  seg1 = ifo1._allsegments
  seg2 = ifo2._allsegments
  vet1 = ifo1._vetoes
  vet2 = ifo2._vetoes
  unv1 = ifo1._unvetoed
  unv2 = ifo2._unvetoed
  minl1 = ifo1._minlen
  minl2 = ifo2._minlen

  name12 = name1 + name2
  print('Combining ' + name1 + ' and ' + name2 + ' into ' + name12)
  if unvetoed:
    seg12 = unv1.intersectSegments(unv2)  # only uses unvetoed long segment list for each IFO
  else:
    seg12 = seg1.intersectSegments(seg2)  # uses full segment lists including small&vetoed segments
  vet12 = vet1.unionSegments(vet2)
  minl12 = max(minl1, minl2)
  ifo12 = IFO(name12, seg12, vet12, minl12)
#  if ifo12
  return ifo12


def getTriples(ifo1, ifo2, ifo3, unvetoed=False):
  '''Combines 3 IFO objects into one'''
  ifo12 = getDoubles(ifo1, ifo2, unvetoed)
  ifo123 = getDoubles(ifo12, ifo3, unvetoed)
  return ifo123


def generateTimeslides(tdict, n, ref=None, outfolder=None, verbose=False):
  '''
     Generate a list of injection times and timeslides,
     given a dictionary of (single) injection times
  '''
  ifos = tdict.keys()
  # nifo = len(ifos)
  if ref is None:
    ref = ifos[0]
  injtimes = []
  injdict = {}
  slidedict = {}
  for k in tdict.keys():
    injdict[k] = []
    slidedict[k] = []
  for i in arange(n):
    for ifo in ifos:
      injdict[ifo].append(rd.sample(tdict[ifo], 1)[0])

  injtimes = injdict[ref]
  for ifo in ifos:
    slidedict[ifo] = array(injdict[ifo]) - array(injdict[ref])

  if outfolder is None:
    return injtimes, slidedict
  else:
    label=''.join(ifos) + '_' + str(n)
    injtimesfile = os.path.join(outfolder,'injtimes_'+label+'.dat')
    slidefilename = os.path.join(outfolder,'timeslides_'+label+'.dat')
    slidefile = open(slidefilename, 'w')
    header = " ".join(ifos)
    slidefile.write(header+'\n')
    slidedata = array(slidedict.values()).T
    if verbose:
      print(injtimes)
      print(slidedict)
    savetxt(injtimesfile, array(injtimes), fmt='%i')
    savetxt(slidefile, slidedata, delimiter=' ', fmt='%i')
    slidefile.close()
    return (injtimesfile, slidefilename)



###########################
#
#  MAIN FUNCTION
#
###########################

if __name__ == "__main__":

###########################################
#
#  Parse arguments (argparser)
#
###########################################

  #parser = OptionParser()
  parser = argparse.ArgumentParser(description="Reads segment files and outputs injection times in veto-free segments.")


  parser.add_argument('-i', '--ifos', nargs=3, type=str, metavar='IFO', dest='ifos', help='List of IFOs to be considered (currently takes 3)') # FIXME: Variable length???
  parser.add_argument('-s', '--segfiles', nargs=3, type=str, metavar='FILE', dest='segf', help='path to the segment files for IFOS in the order entered above')
  parser.add_argument('-v', '--vetofiles', nargs=3, type=str, metavar='FILE', dest='vetof', help='path to the veto files for IFOS in the order entered above')
  parser.add_argument('-o', '--outfolder', type=str, metavar='FILE', dest='outf', help='path to the output folder', default='.')
  parser.add_argument('-l', '--length', type=int, metavar='INT', dest='length', help='length of signal segments in seconds', default=45.0)
  parser.add_argument('-P', '--psdlength', type=int, metavar='INT', dest='psdlength', help='minimum length for calculating PSD in seconds', default=1024)
  parser.add_argument('-N', '--Ninj', type=int, dest="Ninj", help="# injections", default=1000)
  parser.add_argument('-I', '--interval', type=int, metavar='INT', dest='interval', help='minimum space between injections in seconds (defaults to injecting only in the middle of segments)', default=None)
  parser.add_argument('-S', '--singles', action="store_true", dest="singles", help="process singles (default for timeslides)", default=False)
  parser.add_argument('-D', '--doubles', action="store_true", dest="doubles", help="process doubles (N/A for timeslides)", default=False)
  parser.add_argument('-T', '--triples', action="store_true", dest="triples", help="process triples (N/A for timeslides)", default=False)
  parser.add_argument('-t', '--timeslides', action="store_true", dest="timeslides", help="Enable timeslides. This automatically deactivates --triples and --doubles.", default=False)
  parser.add_argument('-p', '--plot', action="store_true", dest="plotsegdist", help="plot cumulative segment length distribution", default=False)
  parser.add_argument('-c', '--checkslides', action="store_true", dest="check", help="verbose check that all timeslid times are in unvetoed segments", default=False)
  parser.add_argument('-u', '--dumpunvetoed', action="store_true", dest="dumpunvetoed", help="dump list of unvetoed segments to file and plot segments", default=False)

  #TODO:
  # * add option to dump to pickle (IFOlist, doubleIFOs, tripleIFO)
  # * add option to read IFOs from pickle

  args = parser.parse_args()

  outfolder = str(args.outf)
  segfiles = args.segf
  vetofiles = args.vetof
  ifos = args.ifos
  nifos = len(ifos)
  Ninj = args.Ninj
  seglen = args.length
  psdlen = args.psdlength
  plotsegdist = args.plotsegdist
  dumpunvetoed = args.dumpunvetoed
  interval = args.interval
  singles = args.singles
  doubles = args.doubles
  triples = args.triples
  timeslides = args.timeslides
  check = args.check
  maxplot = 40000

  # Ensure that directories exist
  ensure_dir(outfolder)
  if dumpunvetoed:
    ensure_dir(os.path.join(outfolder, 'segments'))

  # Choose where to inject (currently middle of segment or every interval)
  if interval is None:
    whereInj = 'middle'
  else:
    whereInj = 'often'

  # Read segments data in rawlist and only keep long enough segments in fitlist for each of the IFOs
  timesdict = {}
  minlen = max(psdlen, seglen)

  # Read in segment files and veto files and create IFO list
  IFOlist = []
  for i in arange(nifos):
    IFOlist.append(IFO(ifos[i], segfiles[i], vetofiles[i], minlen))
    # Generate triggers for each IFO (for timeslides only)
    if timeslides:
      if plotsegdist:
        IFOlist[-1].plotCumulativeDurations(os.path.join(outfolder,'singleseg_' + IFOlist[-1]._name +'.png'), maxplot)
      timesdict[ifos[i]] = IFOlist[-1].getTrigTimes(whereInj=whereInj, interval=interval, lmargin=seglen)

  # Generate timeslide output
  if timeslides:
    ensure_dir(os.path.join(outfolder, 'timeslides'))
    generateTimeslides(timesdict, Ninj, ref='H1', outfolder=os.path.join(outfolder,'timeslides'))
    singles = True
    doubles = False
    triples = False

  if singles:
    for single in IFOlist:
      if plotsegdist:
        single.plotCumulativeDurations(os.path.join(outfolder, 'singleseg_' + single._name + '.png'), maxplot)
      if not timeslides:
        single.getTrigTimes(whereInj=whereInj, n=Ninj, interval=interval, lmargin=seglen, outfile=os.path.join(outfolder, 'injtimes_' + single._name + '_' + str(single._minlen) + '.dat'))
        # Check unvetoedness and prints the distance of trigger times to edges of unvetoed segments
        if check:
          print('checking single ' + single._name)
          for time in single.getTrigTimes(whereInj=whereInj, n=Ninj, interval=interval, lmargin=seglen):
            dist = single._unvetoed.timeIn(time)
            print(dist, dist > 0)
      # Dump unvetoed segment list to file and plot segments to file
      if dumpunvetoed:
        single._unvetoed.printToFile(os.path.join(outfolder, 'segments', single._name + '_unvetoed_segs.dat'))
        figseg = single._unvetoed.plotSegments(os.path.join(outfolder, 'segments', single._name + '_unvetoed_segs.pdf'), title=single._name + ' unvetoed segments', segcolor=colordict[single._name])


  # Combine IFOs into doubles
  if doubles:
    doubleIFOs = []
    doubleIFOs.append( getDoubles(IFOlist[0], IFOlist[1]) ) #, unvetoed=not plotsegdist) )  #FIXME: test it with plots!
    doubleIFOs.append( getDoubles(IFOlist[1], IFOlist[2]) ) #, unvetoed=not plotsegdist) )  # speedup is ~3s for full triples set
    doubleIFOs.append( getDoubles(IFOlist[0], IFOlist[2]) ) #, unvetoed=not plotsegdist) )  # not worth overall?

    # Generate output for doubles
    for double in doubleIFOs:
      double.getTrigTimes(whereInj=whereInj, n=Ninj, interval=interval, lmargin=seglen, outfile=os.path.join(outfolder, 'injtimes_' + double._name + '_' + str(double._minlen) +'.dat'))
      if plotsegdist:
          double.plotCumulativeDurations(os.path.join(outfolder,'doubleseg_' + double._name +'.png'), maxplot)
      # Check unvetoedness and prints the distance of trigger times to edges of unvetoed segments
      if check:
        print('checking double ' + double._name)
        for time in double.getTrigTimes(whereInj=whereInj, n=Ninj, interval=interval, lmargin=seglen):
          dist = double._unvetoed.timeIn(time)
          print(dist, dist > 0)
      # Dump unvetoed segment list to file and plot segments to file
      if dumpunvetoed:
        double._unvetoed.printToFile(os.path.join(outfolder, 'segments', double._name + '_unvetoed_segs.dat'))
        figseg = double._unvetoed.plotSegments(os.path.join(outfolder, 'segments', double._name + '_unvetoed_segs.pdf'), title=double._name + ' unvetoed segments')


  # Combine IFOs into triples
  if triples:
    # If doubles are combined half of the work is already done
    if doubles:
      tripleIFO = getDoubles(doubleIFOs[0], IFOlist[2]) #, unvetoed=not plotsegdist)
    else:
      tripleIFO = getTriples(IFOlist[0], IFOlist[1], IFOlist[2]) #, unvetoed=not plotsegdist)

    # Generate output for doubles
    tripleIFO.getTrigTimes(whereInj=whereInj, n=Ninj, interval=interval, lmargin=seglen, outfile=os.path.join(outfolder,'injtimes_' + tripleIFO._name + '_' + str(tripleIFO._minlen) +'.dat'))
    if plotsegdist:
      tripleIFO.plotCumulativeDurations(os.path.join(outfolder,'tripleseg_' + tripleIFO._name +'.png'), maxplot)
      # Check unvetoedness and prints the distance of trigger times to edges of unvetoed segments
    if check:
      print('checking triple ' + tripleIFO._name)
      for time in tripleIFO.getTrigTimes(whereInj=whereInj, n=Ninj, interval=interval, lmargin=seglen):
        dist = tripleIFO._unvetoed.timeIn(time)
        print(dist, dist > 0)
    # Dump unvetoed segment list to file and plot segments to file
    if dumpunvetoed:
      tripleIFO._unvetoed.printToFile(os.path.join(outfolder, 'segments', tripleIFO._name + '_unvetoed_segs.dat'))
      figseg = tripleIFO._unvetoed.plotSegments(os.path.join(outfolder, 'segments', tripleIFO._name + '_unvetoed_segs.pdf'), title=tripleIFO._name + ' unvetoed segments')

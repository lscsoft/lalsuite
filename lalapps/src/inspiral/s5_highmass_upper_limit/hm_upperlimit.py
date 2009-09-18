#!/usr/bin/python
import scipy
from scipy import interpolate
import numpy
try:
        import sqlite3
except ImportError:
        # pre 2.5.x
        from pysqlite2 import dbapi2 as sqlite3
from math import *
import sys
import glob
import copy
from optparse import OptionParser

from glue import segments
from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import dbtables
from glue.ligolw import utils
from glue.ligolw import table
from glue import segmentsUtils

from pylal import db_thinca_rings
from pylal import llwapp
from pylal import rate
from pylal import SimInspiralUtils
from pylal.xlal.date import LIGOTimeGPS

lsctables.LIGOTimeGPS = LIGOTimeGPS

class upper_limit(object):

  def __init__(self, flist, opts):
    self.far = {}
    self.segments = segments.segmentlistdict()
    self.non_inj_fnames = []
    self.inj_fnames = []
    #self.der_fit = None
    self.twoDMassBins = None
    #self.dBin = {}
    self.gw = None
    self.found = {}
    self.missed = {}
    self.wnfunc = None
    self.opts = opts
    if opts.bootstrap_iterations: self.bootnum = int(opts.bootstrap_iterations)
    else: self.bootnum = 100
    self.veto_segments = segments.segmentlistdict()
    self.zero_lag_segments = {}
    self.instruments = []
    self.livetime = {}
    self.minmass = None
    self.maxmass = None
    self.mintotal = None
    self.maxtotal = None

    for f in flist: 
      if opts.verbose: print >> sys.stderr, "Gathering stats from: %s...." % (f,)
      working_filename = dbtables.get_connection_filename(f, verbose = opts.verbose)
      connection = sqlite3.connect(working_filename)
      dbtables.DBTable_set_connection(connection)
      xmldoc = dbtables.get_xml(connection)

      # look for a sim table
      try:
        sim_inspiral_table = table.get_table(xmldoc, dbtables.lsctables.SimInspiralTable.tableName)
        self.inj_fnames.append(f)
        sim = True
      except ValueError:
        self.non_inj_fnames.append(f)
        sim = False

      if not sim: 
        if opts.veto_segments_name is not None: self.veto_segments = db_thinca_rings.get_veto_segments(connection, opts.veto_segments_name)
        self.get_instruments(connection)
        self.segments += db_thinca_rings.get_thinca_zero_lag_segments(connection, program_name = opts.live_time_program)
        self.get_far_thresholds(connection)
      else: 
        self.get_mass_ranges(connection)
      

      #connection.close()
      dbtables.discard_connection_filename(f, working_filename, verbose = opts.verbose)
      dbtables.DBTable_set_connection(None)      

    # FIXME Do these have to be done by instruments?
    self.segments -= self.veto_segments

    # compute far, segments and livetime by instruments
    for i in self.instruments:
      self.far[i] = min(self.far[i])
      # FIXME this bombs if any of the FARS are zero. maybe it should continue
      # and just remove that instrument combo from the calculation
      if self.far[i] == 0: 
        print >> sys.stderr, "Encountered 0 FAR in %s, ABORTING" % (i,)
        sys.exit(1)
      self.zero_lag_segments[i] = self.segments.intersection(i) - self.segments.union(set(self.segments.keys()) - i)
      # Livetime must have playground removed
      self.livetime[i] = float(abs(self.zero_lag_segments[i] - segmentsUtils.S2playground(self.segments.extent_all())))
      if opts.verbose: print >> sys.stderr, "%s FAR %e, livetime %f" % (",".join(sorted(list(i))), self.far[i], self.livetime[i])

    # get a 2D mass binning
    self.twoDMassBins = self.get_2d_mass_bins(self.minmass, self.maxmass, opts.mass_bins)

  def get_distance_bins(self, instruments, found=None, missed=None):
    if not found and not missed: found, missed = self.get_injections(instruments)
    if not found:  
      print >>sys.stderr,"Found no injections cannot compute distance bins ABORTING"
      sys.exit(1)
    #Give the bins some padding based on the errors
    maxdist = max([s.distance for s in found])
    mindist = min([s.distance for s in found])
    if (maxdist < 0) or (mindist < 0) or (mindist > maxdist):
      print >>sys.stderr, "minimum and maximum distances are screwy, maybe the distance errors given in the options don't make sense? ABORTING"
      sys.exit(1)
    self.dBin[instruments] = rate.LogarithmicBins(mindist,maxdist,self.opts.dist_bins)


  def set_instruments_to_calculate(self):
    if not opts.instruments: return self.instruments
    if opts.instruments in self.instruments:
      return frozenset(lsctables.instrument_set_from_ifos(i[0]))
    else:
      print >> sys.stderr, "Instruments %s do not exist in DB, nothing will be calculated" % (str(frozenset(lsctables.instrument_set_from_ifos(i[0]))),)  
      return []

  def get_mass_ranges(self, connection):
    query = 'SELECT MIN(mass1), MIN(mass2), MAX(mass1), MAX(mass2), MIN(mass1+mass2), MAX(mass1+mass2) FROM sim_inspiral;'
    for v in connection.cursor().execute(query): 
      if self.minmass: self.minmass = min([v[0], v[1], self.minmass])
      else: self.minmass = min(v[0], v[1])
      if self.maxmass: self.maxmass = max([v[2], v[3], self.maxmass])
      else: self.maxmass = max(v[2], v[3])
      if self.mintotal: self.mintotal = min([v[4], self.mintotal])
      else: self.mintotal = v[4]
      if self.maxtotal: self.maxtotal = max([v[5], self.maxtotal])
      else: self.maxtotal = v[5]

  def get_instruments(self, connection):
    for i in connection.cursor().execute('SELECT DISTINCT(instruments) FROM coinc_event'):
      self.instruments.append(frozenset(lsctables.instrument_set_from_ifos(i[0])))
    return self.instruments    

  def get_far_thresholds(self,connection):
    """
    return the false alarm rate of the most rare zero-lag coinc, and a
    dictionary of the thinca segments indexed by instrument.
    """
    live_time_program = opts.live_time_program
    verbose = opts.verbose
    if self.opts.verbose: print >>sys.stderr, "getting FAR thresholds..."
    # extract false alarm rate threshold
    query = 'CREATE TEMPORARY TABLE distinct_instruments AS SELECT DISTINCT(instruments) as instruments FROM coinc_event;'
    connection.cursor().execute(query)
    
    def create_is_playground_func( connection, playground_segs = segmentsUtils.S2playground(self.segments.extent_all()) ):
      """
      Construct the is_playground() SQL function.
      """
      connection.create_function("is_playground", 2, lambda seconds, nanoseconds: LIGOTimeGPS(seconds, nanoseconds) in playground_segs)

    create_is_playground_func(connection)

    query = 'SELECT distinct_instruments.instruments, (SELECT MIN(coinc_inspiral.combined_far) AS combined_far FROM coinc_inspiral JOIN coinc_event ON (coinc_inspiral.coinc_event_id == coinc_event.coinc_event_id) WHERE coinc_event.instruments == distinct_instruments.instruments AND NOT EXISTS(SELECT * FROM time_slide WHERE time_slide.time_slide_id == coinc_event.time_slide_id AND time_slide.offset != 0) AND NOT is_playground(coinc_inspiral.end_time, coinc_inspiral.end_time_ns) ) FROM distinct_instruments;'

    for inst, far in connection.cursor().execute(query):
      inst = frozenset(lsctables.instrument_set_from_ifos(inst))
      self.far.setdefault(inst,[])
      self.far[inst].append(far)
    query = 'DROP TABLE distinct_instruments'
    connection.cursor().execute(query)

  def get_volume_derivative(self,instruments):
    injfnames = self.inj_fnames
    FAR = self.far[instruments]
    zero_lag_segments = self.zero_lag_segments[instruments]
    gw = self.gw
    twoDMassBins = self.twoDMassBins

    #determine binning up front for infinite FAR
    found, missed = self.get_injections(instruments, FAR=float("inf"))
    dbin = rate.LogarithmicBins(min([l.distance for l in found]),max([l.distance for l in found]), int(self.opts.dist_bins))
    
    livetime = float(abs(zero_lag_segments))
    FARh = FAR*100000
    FARl = FAR*0.001
    nbins = 5
    FARS = rate.LogarithmicBins(FARl, FARh, nbins)
    vA = []
    vA2 = []
    for far in FARS.centres():
      print >>sys.stderr, "computing volume at FAR " + str(far)
      vAt, vA2t = self.twoD_SearchVolume(instruments, dbin=dbin, FAR = far, bootnum=1)
      # we need to compute derivitive of log according to ul paper
      vAt.array = scipy.log10(vAt.array + 0.001)
      vA.append(vAt)
    # the derivitive is calcuated with respect to FAR * t
    FARTS = rate.LogarithmicBins(FARl * livetime, FARh * livetime, nbins)
    return self._derivitave_fit(FARTS, FAR * livetime, vA, twoDMassBins)
  
  
  def _derivitave_fit(self, farts, FARt, vAs, twodbin):
    '''
       Relies on scipy spline fits for each mass bin
       to find the derivitave of the volume at a given
       FAR.  See how this works for a simple case where
       I am clearly giving it a parabola.  To high precision it calculates
       the proper derivitave. 
       A = [1, 4, 9, 16, 25, 36, 49, 64, 81, 100]
       B = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
       C = interpolate.splrep(B,A,s=0, k=4)
       interpolate.splev(5,C,der=1) 
       10.000
    '''
    dA = rate.BinnedArray(twodbin)
    for m1 in range(dA.array.shape[0]):
      for m2 in range(dA.array.shape[1]):
        da = []
        for f in farts.centres():
          da.append(vAs[farts[f]].array[m1][m2])
        fit = interpolate.splrep(farts.centres(),da,k=4)
        val = interpolate.splev(FARt,fit,der=1)
        #print val
        # FIXME this prevents negative derivitives arising from bad fits
        if val < 0: val = 0
        dA.array[m1][m2] = val # minus the derivitave
    return dA

  def get_injections(self, instruments, FAR=None):
    injfnames = self.inj_fnames
    if not FAR: FAR = self.far[instruments]
    zero_lag_segments = self.zero_lag_segments[instruments]
    verbose = self.opts.verbose
    """
    """
    def injection_was_made(geocent_end_time, geocent_end_time_ns, zero_lag_segments = zero_lag_segments):
      """
      return True if injection was made in the given segmentlist
      """
      return lsctables.LIGOTimeGPS(geocent_end_time, geocent_end_time_ns) in zero_lag_segments

    found = []
    missed = []
    print >>sys.stderr, ""
    for cnt, f in enumerate(injfnames):
      print >>sys.stderr, "getting injections below FAR: " + str(FAR) + ":\t%.1f%%\r" % (100.0 * cnt / len(injfnames),),
      working_filename = dbtables.get_connection_filename(f, tmp_path = None, verbose = verbose)
      connection = sqlite3.connect(working_filename)
      connection.create_function("injection_was_made", 2, injection_was_made)

      make_sim_inspiral = lsctables.table.get_table(dbtables.get_xml(connection), lsctables.SimInspiralTable.tableName)._row_from_cols

      for values in connection.cursor().execute("""
SELECT
  sim_inspiral.*,
  -- true if injection matched a coinc below the false alarm rate threshold
  EXISTS (
    SELECT
      *
    FROM
      coinc_event_map AS mapa
      JOIN coinc_event_map AS mapb ON (
        mapa.coinc_event_id == mapb.coinc_event_id
      )
      JOIN coinc_inspiral ON (
        mapb.table_name == "coinc_event"
        AND mapb.event_id == coinc_inspiral.coinc_event_id
      )
    WHERE
      mapa.table_name == "sim_inspiral"
      AND mapa.event_id == sim_inspiral.simulation_id
      AND coinc_inspiral.combined_far < ?
  )
FROM
  sim_inspiral
WHERE
  -- only interested in injections that were injected
  injection_was_made(sim_inspiral.geocent_end_time, sim_inspiral.geocent_end_time_ns)
      """, (FAR,)):
        sim = make_sim_inspiral(values)
        if values[-1]:
          found.append(sim)
        else:
          missed.append(sim)

      # done
      connection.close()
      dbtables.discard_connection_filename(f, working_filename, verbose = verbose)
      dbtables.DBTable_set_connection(None)

    print >>sys.stderr, "\nFound = %d Missed = %d" % (len(found), len(missed))
    return found, missed


  def trim_mass_space(self, eff, instruments, minthresh=0.0, minM=25.0, maxM=100.0):
    """
    restricts array to only have data within the mass space and sets everything
    outside the mass space to some canonical value, minthresh
    """
    twodbin = self.twoDMassBins
    x = eff.array.shape[0]
    y = eff.array.shape[1]
    c1 = twodbin.centres()[0]
    c2 = twodbin.centres()[1]
    numbins = 0
    for i in range(x):
      for j in range(y):
        if c1[i] > c2[j] or (c1[i] + c2[j]) > maxM or (c1[i]+c2[j]) < minM: eff.array[i][j] = minthresh
        else: numbins+=1
    print "found " + str(numbins) + " bins within total mass"

  def fix_masses(self, sims):
    """
    Function to duplicate the mass pairs to remove edge effects 
    on the equal mass line, takes a list of sim rows
    """
    sims2 = []
    for l in sims:
      l2 = copy.deepcopy(l)
      l2.mass1 = l.mass2
      l2.mass2 = l.mass1
      sims2.append(l2)
    sims.extend(sims2)

  def get_2d_mass_bins(self, low, high, bins):
    """
    Given the component mass range low, high of the search it will
    return 2D bins with size bins in each direction
    """
    mass1Bin = rate.LinearBins(low,high,bins)
    mass2Bin = rate.LinearBins(low,high,bins)
    twoDMB=rate.NDBins( (mass1Bin,mass2Bin) )
    return twoDMB
    
  def _scramble_pop(self, m, f):
    """
    A function to draw a new injection sample in the "boot strap" method 
    http://en.wikipedia.org/wiki/Bootstrapping_(statistics) 
    and included refereneces.
    This was used in the stack-a-flare search to get MC errors etc. 
    """
    inj = m+f
    ix = scipy.random.randint(0,len(inj), (len(inj),))
    #return new missed, found
    missed = [inj[i] for i in ix if i < len(m) ]
    found = [inj[i] for i in ix if i >=len(m) ]
    return missed, found

  def _scramble_dist(self, inj, relerr, syserr):
    """
    function to handle random calibration error.  Individually srambles the distances
    of injection by an error.
    """
    dist_array = numpy.zeros(len(inj))
    for i,sim in enumerate(inj):
      dist_array[i] = sim.distance * (1.0-syserr) * float(scipy.exp( relerr * scipy.random.standard_normal(1))) 
    return dist_array

  def live_time_array(self, instruments):
    """
    return an array of live times, note every bin will be the same :) it is just a 
    convenience.
    """
    live_time = rate.BinnedArray(self.twoDMassBins)
    live_time.array += 1.0
    live_time.array *= self.livetime[instruments]
    return live_time
    

  def twoD_SearchVolume(self, instruments, dbin=None, FAR=None, bootnum=None, derr=0.197, dsys=0.074):
    """ 
    Compute the search volume in the mass/mass plane, bootstrap
    and measure the first and second moment (assumes the underlying 
    distribution can be characterized by those two parameters) 
    This is gonna be brutally slow
    derr = (0.134**2+.103**2+.102**2)**.5 = 0.197 which is the 3 detector 
    calibration uncertainty in quadrature.  This is conservative since some injections
     will be H1L1 and have a lower error of .17
    the dsys is the DC offset which is the max offset of .074. 
    """

    if not FAR: FAR = self.far[instruments]
    found, missed = self.get_injections(instruments, FAR)
    twodbin = self.twoDMassBins
    wnfunc = self.gw
    livetime = self.livetime[instruments]
    if not bootnum: bootnum = self.bootnum

    if wnfunc: wnfunc /= wnfunc[(wnfunc.shape[0]-1) / 2, (wnfunc.shape[1]-1) / 2]

    x = twodbin.shape[0]
    y = twodbin.shape[1]
    z = int(self.opts.dist_bins)

    rArrays = []
    volArray=rate.BinnedArray(twodbin)
    volArray2=rate.BinnedArray(twodbin)
    #set up ratio arrays for each distance bin
    for k in range(z):
      rArrays.append(rate.BinnedRatios(twodbin))

    # Bootstrap to account for errors
    for n in range(bootnum):
      #initialize by setting these to zero
      for k in range(z):
        rArrays[k].numerator.array = numpy.zeros(rArrays[k].numerator.bins.shape)
        rArrays[k].denominator.array = numpy.zeros(rArrays[k].numerator.bins.shape)
      #Scramble the inj population and distances
      if bootnum > 1: 
        sm, sf = self._scramble_pop(missed, found)
        # I make a separate array of distances to speed up this calculation
        f_dist = self._scramble_dist(sf, derr, dsys)
      else: 
        sm, sf = missed, found
        f_dist = numpy.array([l.distance for l in found])
     
      # compute the distance bins
      if not dbin: 
        dbin = rate.LogarithmicBins(min(f_dist),max(f_dist), z)
      #else: print dbin.centres()
      

      # get rid of all missed injections outside the distance bins
      # to prevent binning errors
      sm, m_dist = self.cut_distance(sm, dbin)
      sf, f_dist = self.cut_distance(sf, dbin)


      for i, l in enumerate(sf):#found:
        tbin = rArrays[dbin[f_dist[i]]]
        tbin.incnumerator( (l.mass1, l.mass2) )
      for i, l in enumerate(sm):#missed:
        tbin = rArrays[dbin[m_dist[i]]]
        tbin.incdenominator( (l.mass1, l.mass2) )
    
      tmpArray2=rate.BinnedArray(twodbin) #start with a zero array to compute the mean square
      for k in range(z): 
        tbins = rArrays[k]
        tbins.denominator.array += tbins.numerator.array
        if wnfunc: rate.filter_array(tbins.denominator.array,wnfunc)
        if wnfunc: rate.filter_array(tbins.numerator.array,wnfunc)
        tbins.regularize()
        # logarithmic(d)
        integrand = 4.0 * pi * tbins.ratio() * dbin.centres()[k]**3 * dbin.delta
        volArray.array += integrand
        tmpArray2.array += integrand #4.0 * pi * tbins.ratio() * dbin.centres()[k]**3 * dbin.delta
        print >>sys.stderr, "bootstrapping:\t%.1f%% and Calculating smoothed volume:\t%.1f%%\r" % ((100.0 * n / bootnum), (100.0 * k / z)),
      tmpArray2.array *= tmpArray2.array
      volArray2.array += tmpArray2.array
    
    print >>sys.stderr, "" 
    #Mean and variance
    volArray.array /= bootnum
    volArray2.array /= bootnum
    volArray2.array -= volArray.array**2 # Variance
    volArray.array *= livetime
    volArray2.array *= livetime*livetime # this gets two powers of live time
    return volArray, volArray2
 

  def cut_distance(self, sims, dbin):
    """
    Exclude sims outside some distance range to avoid errors when binning
    """
    mnd = min(dbin.lower())
    mxd = max(dbin.upper())
    #print len(sims), len(m_dist)
    #m_dist_low = m_dist[m_dist >= mnd]
    #m_dist_low = m_dist_low[m_dist_low <= mxd]
    out =  [sim for sim in sims if mnd <= sim.distance <= mxd]
    return out, numpy.array([l.distance for l in out])

#### STUFF FOR PLAYGROUND ####

#playground_segs = segmentsUtils.S2playground(self.seglists.extent_all())

#def create_is_playground_func(connection, playground_segs):
#        """
#        Construct the is_playground() SQL function.
#        """
#        connection.create_function("is_playground", 2, lambda seconds, nanoseconds: LIGOTimeGPS(seconds, nanoseconds) in playground_segs)
 

######################## ACTUAL PROGRAM #######################################
###############################################################################
###############################################################################


def parse_command_line():
  parser = OptionParser(version = "%prog CVS $Id$", usage = "%prog [options] [file ...]", description = "%prog computes mass/mass upperlimit")
  parser.add_option("--instruments", metavar = "name[,name,...]", help = "Set the list of instruments.  Required.  Example \"H1,H2,L1\"")
  parser.add_option("--live-time-program", default = "thinca", metavar = "name", help = "Set the name of the program whose rings will be extracted from the search_summary table.  Default = \"thinca\".")
  parser.add_option("--output-name-tag", default = "", metavar = "name", help = "Set the file output name tag, real name is 2Dsearchvolume-<tag>-<ifos>.xml")
  parser.add_option("--bootstrap-iterations", default = 1000, metavar = "integer", type = "int", help = "Number of iterations to compute mean and variance of volume MUST BE GREATER THAN 1 TO GET USABLE NUMBERS, a good number is 10000")
  parser.add_option("--dist-bins", default = 50, metavar = "integer", type = "int", help = "Number of mass bins")
  parser.add_option("--d-err", default = 0.197, metavar = "float", type = "int", help = "random calibration error on distance")
  parser.add_option("--d-sys-err", default = 0.074, metavar = "float", type = "int", help = "systematic calibration error on distance (should use worst)")
  parser.add_option("--mass-bins", default = 11, metavar = "integer", type = "int", help = "Number of mass bins along 1 dimension (Note the total number of mass bins will generally be less than --mass-bins * --mass-bins once the actual parameter space is carved out)")
  parser.add_option("--veto-segments-name", default = "vetoes", help = "Set the name of the veto segments to use from the XML document.")
  parser.add_option("--verbose", action = "store_true", help = "Be verbose.")

  opts, filenames = parser.parse_args()

  if opts.instruments: opts.instruments = lsctables.instrument_set_from_ifos(opts.instruments)
  if not filenames: 
    print >>sys.stderr, "must specify at least one database file"
    sys.exit(1)
  return opts, filenames

# FIXME These values should probably be command line arguments or derived from the database
secs_in_year = 31556926.0

opts, filenames = parse_command_line()

#initialize an upper limit class
UL = upper_limit(filenames, opts)

#loop over the requested instruments
for instruments in UL.set_instruments_to_calculate():
  if opts.verbose: print >>sys.stderr, "calculating upper limit for %s" % (",".join(sorted(list(instruments))),)

  #compute volume derivitive
  dvA = UL.get_volume_derivative(instruments)

  #compute volume first and second moments
  vA, vA2 = UL.twoD_SearchVolume(instruments)

  # get an array of livetimes for convenience
  ltA = UL.live_time_array(instruments)

  # FIXME convert to years (use some lal or pylal thing in the future)
  vA.array /= secs_in_year
  vA2.array /= secs_in_year * secs_in_year #two powers for this squared quantity

  #Trim the array to have sane values outside the total mass area of interest
  try: minvol = scipy.unique(vA.array)[1]/10.0
  except: minvol = 0
  UL.trim_mass_space(dvA, instruments, minthresh=0.0, minM=UL.mintotal, maxM=UL.maxtotal)
  UL.trim_mass_space(vA, instruments, minthresh=minvol, minM=UL.mintotal, maxM=UL.maxtotal)
  UL.trim_mass_space(vA2, instruments, minthresh=0.0, minM=UL.mintotal, maxM=UL.maxtotal)

  #output an XML file with the result
  xmldoc = ligolw.Document()
  xmldoc.appendChild(ligolw.LIGO_LW())
  xmldoc.childNodes[-1].appendChild(rate.binned_array_to_xml(vA, "2DsearchvolumeFirstMoment"))
  xmldoc.childNodes[-1].appendChild(rate.binned_array_to_xml(vA2, "2DsearchvolumeSecondMoment"))
  xmldoc.childNodes[-1].appendChild(rate.binned_array_to_xml(dvA, "2DsearchvolumeDerivative"))

  # DONE with vA, so it is okay to mess it up...
  # Compute range 
  vA.array = (vA.array * secs_in_year / UL.livetime[instruments] / (4.0/3.0 * pi)) **(1.0/3.0)
  UL.trim_mass_space(vA, instruments, minthresh=0.0, minM=UL.mintotal, maxM=UL.maxtotal)
  xmldoc.childNodes[-1].appendChild(rate.binned_array_to_xml(vA, "2DsearchvolumeDistance"))

  # make a live time 
  UL.trim_mass_space(ltA, instruments, minthresh=0.0, minM=UL.mintotal, maxM=UL.maxtotal)
  xmldoc.childNodes[-1].appendChild(rate.binned_array_to_xml(ltA, "2DsearchvolumeLiveTime"))

  utils.write_filename(xmldoc, "2Dsearchvolume-%s-%s.xml" % (opts.output_name_tag, "".join(sorted(list(instruments)))))

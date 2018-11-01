
#flow DAG Class definitions for LALInference Pipeline
# (C) 2012 John Veitch, Vivien Raymond, Kiersten Ruisard, Kan Wang

import itertools
import glue
from glue import pipeline,segmentsUtils,segments
from glue.ligolw import ligolw, lsctables, utils
import os
import socket
from lalapps import inspiralutils
import uuid
import ast
import pdb
import string
from math import floor,ceil,log,pow
import sys
import random
from itertools import permutations
import shutil
import numpy as np
import math
from six.moves import range
from six import next
from functools import reduce

# We use the GLUE pipeline utilities to construct classes for each
# type of job. Each class has inputs and outputs, which are used to
# join together types of jobs into a DAG.

def guess_url(fslocation):
    """
    Try to work out the web address of a given path
    """
    SERVER="localhost"
    USER=os.environ['USER']
    HOST=socket.getfqdn()
    if 'public_html' in fslocation:
        k='public_html/'
    elif 'WWW' in fslocation:
        k='WWW/'
    elif 'www_html' in fslocation:
        k='www_html/'
    else:
        k=None
    if k is not None:
        (a,b)=fslocation.split(k)
        webpath=os.path.join('~%s'%(USER),b)
        onweb=True
    else:
        (c,d)=fslocation.split(USER,1)
        for k in ['public_html','WWW','www_html']:
            trypath=c+os.environ['USER']+'/'+k+d
            #Follow symlinks
            if os.path.realpath(trypath)==os.path.normpath(fslocation):
                #(a,b)=trypath.split(k)
                webpath=os.path.join('~%s'%(USER),d)
                onweb=True
                break
            else:
                webpath=fslocation
                onweb=False
    if 'atlas' in HOST:
        url="https://atlas1.atlas.aei.uni-hannover.de/"
    elif 'cit' in HOST or 'caltech' in HOST:
        url="https://ldas-jobs.ligo.caltech.edu/"
    elif 'ligo-wa' in HOST:
        url="https://ldas-jobs.ligo-wa.caltech.edu/"
    elif 'ligo-la' in HOST:
        url="https://ldas-jobs.ligo-la.caltech.edu/"
    elif 'uwm' in HOST or 'nemo' in HOST:
        url="https://ldas-jobs.phys.uwm.edu/"
    elif 'phy.syr.edu' in HOST:
        url="https://sugar-jobs.phy.syr.edu/"
    elif 'arcca.cf.ac.uk' in HOST:
        url="https://geo2.arcca.cf.ac.uk/"
    elif 'vulcan' in HOST:
        url="https://galahad.aei.mpg.de/"
    else:
        if onweb:
            url="http://%s/"%(HOST)
        else:
            url=HOST+':'
    url=url+webpath
    return(url)

class Event():
    """
    Represents a unique event to run on
    """
    new_id=itertools.count()
    def __init__(self,trig_time=None,SimInspiral=None,SimBurst=None,SnglInspiral=None,CoincInspiral=None,event_id=None,timeslide_dict=None,GID=None,ifos=None, duration=None,srate=None,trigSNR=None,fhigh=None,horizon_distance=None):
        self.trig_time=trig_time
        self.injection=SimInspiral
        self.burstinjection=SimBurst
        self.sngltrigger=SnglInspiral
        if timeslide_dict is None:
            self.timeslides={}
        else:
            self.timeslides=timeslide_dict
        self.GID=GID
        self.coinctrigger=CoincInspiral
        if ifos is None:
            self.ifos = []
        else:
            self.ifos = ifos
        self.duration = duration
        self.srate = srate
        self.trigSNR = trigSNR
        self.fhigh = fhigh
        self.horizon_distance = horizon_distance
        if event_id is not None:
            self.event_id=event_id
        else:
            self.event_id=next(Event.new_id)
        if self.injection is not None:
            self.trig_time=self.injection.get_end()
            if event_id is None: self.event_id=int(str(self.injection.simulation_id).split(':')[2])
        if self.burstinjection is not None:
            self.trig_time=self.burstinjection.get_end()
            if event_id is None: self.event_id=int(str(self.burstinjection.simulation_id).split(':')[2])
        if self.sngltrigger is not None:
            self.trig_time=self.sngltrigger.get_end()
            self.event_id=int(str(self.sngltrigger.event_id).split(':')[2])
        if self.coinctrigger is not None:
            self.trig_time=self.coinctrigger.end_time + 1.0e-9 * self.coinctrigger.end_time_ns
        if self.GID is not None:
            self.event_id=int(''.join(i for i in self.GID if i.isdigit()))
        self.engine_opts={}
    def set_engine_option(self,opt,val):
        """
        Can set event-specific options for the engine nodes
        using this option, e.g. ev.set_engine_option('time-min','1083759273')
        """
        self.engine_opts[opt]=val

dummyCacheNames=['LALLIGO','LALVirgo','LALAdLIGO','LALAdVirgo']

def readLValert(threshold_snr=None,gid=None,flow=20.0,gracedb="gracedb",basepath="./",downloadpsd=True,roq=False,service_url=None):
    """
    Parse LV alert file, containing coinc, sngl, coinc_event_map.
    and create a list of Events as input for pipeline
    Based on Chris Pankow's script
    """
    output=[]
    from glue.ligolw import utils as ligolw_utils
    from glue.ligolw import lsctables
    from glue.ligolw import ligolw
    from lal import series as lalseries
    import lal
    from lalsimulation import SimInspiralChirpTimeBound, GetApproximantFromString, IMRPhenomDGetPeakFreq
    from ligo.gracedb.rest import GraceDb, HTTPError
    try:
        from gstlal import reference_psd
    except ImportError:
        reference_psd = None
    cwd=os.getcwd()
    os.chdir(basepath)
    print("Download %s coinc.xml" % gid)
    if service_url is None:
        client = GraceDb()
    else:
        client = GraceDb(service_url=service_url)
    xmldoc = ligolw_utils.load_fileobj(client.files(gid, "coinc.xml"), contenthandler = lsctables.use_in(ligolw.LIGOLWContentHandler))[0]
    ligolw_utils.write_filename(xmldoc, "coinc.xml")
    coinc_events = lsctables.CoincInspiralTable.get_table(xmldoc)
    sngl_event_idx = dict((row.event_id, row) for row in lsctables.SnglInspiralTable.get_table(xmldoc))
    ifos = sorted(coinc_events[0].instruments)
    trigSNR = coinc_events[0].snr
    # Parse PSD
    srate_psdfile=16384
    fhigh=None
    psdfileobj = None
    if downloadpsd:
        print("Download %s psd.xml.gz" % gid)
        try:
            psdfileobj = client.files(gid, "psd.xml.gz")
        except HTTPError:
            print("Failed to download %s psd.xml.gz. lalinference will estimate the psd itself." % gid)
        if psdfileobj is not None:
            if reference_psd is not None:
                xmlpsd = ligolw_utils.load_fileobj(psdfileobj, contenthandler = lalseries.PSDContentHandler)[0]
                psd = lalseries.read_psd_xmldoc(xmlpsd)
                ligolw_utils.write_filename(xmlpsd, "psd.xml.gz", gz = True)
            else:
                open("psd.xml.gz", "wb").write(psdfileobj.read())
            psdasciidic = get_xml_psds(os.path.realpath("./psd.xml.gz"),ifos,os.path.realpath('./PSDs'),end_time=None)
            combine = np.loadtxt(psdasciidic[list(psdasciidic.keys())[0]])
            srate_psdfile = pow(2.0, ceil( log(float(combine[-1][0]), 2) ) ) * 2
    coinc_map = lsctables.CoincMapTable.get_table(xmldoc)
    for coinc in coinc_events:
        these_sngls = [sngl_event_idx[c.event_id] for c in coinc_map if c.coinc_event_id == coinc.coinc_event_id]
        dur=[]
        srate=[]
        horizon_distance=[]
        for e in these_sngls:
            if roq==False:
                chirplen = SimInspiralChirpTimeBound(flow, e.mass1 * lal.MSUN_SI, e.mass2 * lal.MSUN_SI, 0.0, 0.0)
                fstop = IMRPhenomDGetPeakFreq(e.mass1, e.mass2, 0.0, 0.0)
                dur.append(pow(2.0, ceil( log(max(8.0, chirplen + 2.0), 2) ) ) )
                srate.append(pow(2.0, ceil( log(fstop, 2) ) ) * 2)
            # determine horizon distance
            if threshold_snr is not None:
                if e.eff_distance is not None and not math.isnan(e.eff_distance):
                    if e.snr > threshold_snr:
                        horizon_distance.append(e.eff_distance * e.snr / threshold_snr)
                    else:
                        horizon_distance.append(2 * e.eff_distance)
                else:
                    if reference_psd is not None and psdfileobj is not None:
                        if not roq==False:
                            fstop = IMRPhenomDGetPeakFreq(e.mass1, e.mass2, 0.0, 0.0)
                        HorizonDistanceObj = reference_psd.HorizonDistance(f_min = flow, f_max = fstop, delta_f = 1.0 / 32.0, m1 = e.mass1, m2 = e.mass2)
                        horizon_distance.append(HorizonDistanceObj(psd[e.ifo], snr = threshold_snr)[0])
        if srate:
            if max(srate)<srate_psdfile:
                srate = max(srate)
            else:
                srate = srate_psdfile
                if psdfileobj is not None:
                    fhigh = srate_psdfile/2.0 * 0.95 # Because of the drop-off near Nyquist of the PSD from gstlal
        else:
            srate = None
        if dur:
            duration = max(dur)
        else:
            duration = None
        horizon_distance = max(horizon_distance) if len(horizon_distance) > 0 else None
        ev=Event(CoincInspiral=coinc, GID=gid, ifos = ifos, duration = duration, srate = srate,
                 trigSNR = trigSNR, fhigh = fhigh, horizon_distance=horizon_distance)
        output.append(ev)

    print("Found %d coinc events in table." % len(coinc_events))
    os.chdir(cwd)
    return output

def open_pipedown_database(database_filename,tmp_space):
    """
    Open the connection to the pipedown database
    """
    if not os.access(database_filename,os.R_OK):
        raise Exception('Unable to open input file: %s'%(database_filename))
    from glue.ligolw import dbtables
    import sqlite3
    working_filename=dbtables.get_connection_filename(database_filename,tmp_path=tmp_space)
    connection = sqlite3.connect(working_filename)
    if tmp_space:
        dbtables.set_temp_store_directory(connection,tmp_space)
    #dbtables.DBTable_set_connection(connection)
    return (connection,working_filename)

def get_zerolag_lloid(database_connection, dumpfile=None, gpsstart=None, gpsend=None, max_cfar=-1, min_cfar=-1):
    """
    Returns a list of Event objects
    from pipedown data base. Can dump some stats to dumpfile if given,
    and filter by gpsstart and gpsend to reduce the nunmber or specify
    max_cfar to select by combined FAR
    """
    output={}
    if gpsstart is not None: gpsstart=float(gpsstart)
    if gpsend is not None: gpsend=float(gpsend)
    # Get coincs
    get_coincs = "SELECT sngl_inspiral.end_time+sngl_inspiral.end_time_ns*1e-9,sngl_inspiral.ifo,coinc_event.coinc_event_id,sngl_inspiral.snr,sngl_inspiral.chisq,coinc_inspiral.combined_far \
            FROM sngl_inspiral join coinc_event_map on (coinc_event_map.table_name=='sngl_inspiral' and coinc_event_map.event_id ==\
            sngl_inspiral.event_id) join coinc_event on (coinc_event.coinc_event_id==coinc_event_map.coinc_event_id) \
            join coinc_inspiral on (coinc_event.coinc_event_id==coinc_inspiral.coinc_event_id) \
    WHERE coinc_event.time_slide_id=='time_slide:time_slide_id:1'\
            "
    if gpsstart is not None:
        get_coincs=get_coincs+' and coinc_inspiral.end_time+coinc_inspiral.end_time_ns*1.0e-9 > %f'%(gpsstart)
    if gpsend is not None:
        get_coincs=get_coincs+' and coinc_inspiral.end_time+coinc_inspiral.end_time_ns*1.0e-9 < %f'%(gpsend)
    if max_cfar !=-1:
        get_coincs=get_coincs+' and coinc_inspiral.combined_far < %f'%(max_cfar)
    if min_cfar != -1:
        get_coincs=get_coincs+' and coinc_inspiral.combined_far > %f'%(min_cfar)
    db_out=database_connection.cursor().execute(get_coincs)
    extra={}
    for (sngl_time, ifo, coinc_id, snr, chisq, cfar) in db_out:
        coinc_id=int(coinc_id.split(":")[-1])
        if not coinc_id in output.keys():
            output[coinc_id]=Event(trig_time=sngl_time,timeslide_dict={},event_id=int(coinc_id))
            extra[coinc_id]={}
        output[coinc_id].timeslides[ifo]=0
        output[coinc_id].ifos.append(ifo)
        extra[coinc_id][ifo]={'snr':snr,'chisq':chisq,'cfar':cfar}
    if dumpfile is not None:
        fh=open(dumpfile,'w')
        for co in output.keys():
            for ifo in output[co].ifos:
                fh.write('%s %s %s %s %s %s %s\n'%(str(co),ifo,str(output[co].trig_time),str(output[co].timeslides[ifo]),str(extra[co][ifo]['snr']),str(extra[co][ifo]['chisq']),str(extra[co][ifo]['cfar'])))
        fh.close()
    return output.values()

def get_zerolag_pipedown(database_connection, dumpfile=None, gpsstart=None, gpsend=None, max_cfar=-1, min_cfar=-1):
    """
    Returns a list of Event objects
    from pipedown data base. Can dump some stats to dumpfile if given,
    and filter by gpsstart and gpsend to reduce the nunmber or specify
    max_cfar to select by combined FAR
    """
    output={}
    if gpsstart is not None: gpsstart=float(gpsstart)
    if gpsend is not None: gpsend=float(gpsend)
    # Get coincs
    get_coincs = "SELECT sngl_inspiral.end_time+sngl_inspiral.end_time_ns*1e-9,sngl_inspiral.ifo,coinc_event.coinc_event_id,sngl_inspiral.snr,sngl_inspiral.chisq,coinc_inspiral.combined_far \
            FROM sngl_inspiral join coinc_event_map on (coinc_event_map.table_name=='sngl_inspiral' and coinc_event_map.event_id ==\
            sngl_inspiral.event_id) join coinc_event on (coinc_event.coinc_event_id==coinc_event_map.coinc_event_id) \
            join coinc_inspiral on (coinc_event.coinc_event_id==coinc_inspiral.coinc_event_id) \
            WHERE coinc_event.time_slide_id=='time_slide:time_slide_id:10049'\
            "
    if gpsstart is not None:
        get_coincs=get_coincs+' and coinc_inspiral.end_time+coinc_inspiral.end_time_ns*1.0e-9 > %f'%(gpsstart)
    if gpsend is not None:
        get_coincs=get_coincs+' and coinc_inspiral.end_time+coinc_inspiral.end_time_ns*1.0e-9 < %f'%(gpsend)
    if max_cfar !=-1:
        get_coincs=get_coincs+' and coinc_inspiral.combined_far < %f'%(max_cfar)
    if min_cfar != -1:
        get_coincs=get_coincs+' and coinc_inspiral.combined_far > %f'%(min_cfar)
    db_out=database_connection.cursor().execute(get_coincs)
    extra={}
    for (sngl_time, ifo, coinc_id, snr, chisq, cfar) in db_out:
        coinc_id=int(coinc_id.split(":")[-1])
        if not coinc_id in output.keys():
            output[coinc_id]=Event(trig_time=sngl_time,timeslide_dict={},event_id=int(coinc_id))
            extra[coinc_id]={}
        output[coinc_id].timeslides[ifo]=0
        output[coinc_id].ifos.append(ifo)
        extra[coinc_id][ifo]={'snr':snr,'chisq':chisq,'cfar':cfar}
    if dumpfile is not None:
        fh=open(dumpfile,'w')
        for co in output.keys():
            for ifo in output[co].ifos:
                fh.write('%s %s %s %s %s %s %s\n'%(str(co),ifo,str(output[co].trig_time),str(output[co].timeslides[ifo]),str(extra[co][ifo]['snr']),str(extra[co][ifo]['chisq']),str(extra[co][ifo]['cfar'])))
        fh.close()
    return output.values()

def get_timeslides_pipedown(database_connection, dumpfile=None, gpsstart=None, gpsend=None, max_cfar=-1):
    """
    Returns a list of Event objects
    with times and timeslide offsets
    """
    output={}
    if gpsstart is not None: gpsstart=float(gpsstart)
    if gpsend is not None: gpsend=float(gpsend)
    db_segments=[]
    sql_seg_query="SELECT search_summary.out_start_time, search_summary.out_end_time from search_summary join process on process.process_id==search_summary.process_id where process.program=='thinca'"
    db_out = database_connection.cursor().execute(sql_seg_query)
    for d in db_out:
        if d not in db_segments:
            db_segments.append(d)
    seglist=segments.segmentlist([segments.segment(d[0],d[1]) for d in db_segments])
    db_out_saved=[]
    # Get coincidences
    get_coincs="SELECT sngl_inspiral.end_time+sngl_inspiral.end_time_ns*1e-9,time_slide.offset,sngl_inspiral.ifo,coinc_event.coinc_event_id,sngl_inspiral.snr,sngl_inspiral.chisq,coinc_inspiral.combined_far \
                FROM sngl_inspiral join coinc_event_map on (coinc_event_map.table_name == 'sngl_inspiral' and coinc_event_map.event_id \
                == sngl_inspiral.event_id) join coinc_event on (coinc_event.coinc_event_id==coinc_event_map.coinc_event_id) join time_slide\
                on (time_slide.time_slide_id == coinc_event.time_slide_id and time_slide.instrument==sngl_inspiral.ifo)\
                join coinc_inspiral on (coinc_inspiral.coinc_event_id==coinc_event.coinc_event_id) where coinc_event.time_slide_id!='time_slide:time_slide_id:10049'"
    joinstr = ' and '
    if gpsstart is not None:
        get_coincs=get_coincs+ joinstr + ' coinc_inspiral.end_time+coinc_inspiral.end_time_ns*1e-9 > %f'%(gpsstart)
    if gpsend is not None:
        get_coincs=get_coincs+ joinstr+' coinc_inspiral.end_time+coinc_inspiral.end_time_ns*1e-9 <%f'%(gpsend)
    if max_cfar!=-1:
        get_coincs=get_coincs+joinstr+' coinc_inspiral.combined_far < %f'%(max_cfar)
    db_out=database_connection.cursor().execute(get_coincs)
    # Timeslide functionality requires obsolete pylal - will be removed
    import pylal
    from pylal import SnglInspiralUtils
    extra={}
    for (sngl_time, slide, ifo, coinc_id, snr, chisq, cfar) in db_out:
        coinc_id=int(coinc_id.split(":")[-1])
        seg=filter(lambda seg:sngl_time in seg,seglist)[0]
        slid_time = SnglInspiralUtils.slideTimeOnRing(sngl_time,slide,seg)
        if not coinc_id in output.keys():
            output[coinc_id]=Event(trig_time=slid_time,timeslide_dict={},event_id=int(coinc_id))
            extra[coinc_id]={}
        output[coinc_id].timeslides[ifo]=slid_time-sngl_time
        output[coinc_id].ifos.append(ifo)
        extra[coinc_id][ifo]={'snr':snr,'chisq':chisq,'cfar':cfar}
    if dumpfile is not None:
        fh=open(dumpfile,'w')
        for co in output.keys():
            for ifo in output[co].ifos:
                fh.write('%s %s %s %s %s %s %s\n'%(str(co),ifo,str(output[co].trig_time),str(output[co].timeslides[ifo]),str(extra[co][ifo]['snr']),str(extra[co][ifo]['chisq']),str(extra[co][ifo]['cfar'])))
        fh.close()
    return output.values()

def mkdirs(path):
    """
    Helper function. Make the given directory, creating intermediate
    dirs if necessary, and don't complain about it already existing.
    """
    if os.access(path,os.W_OK) and os.path.isdir(path): return
    else: os.makedirs(path)

def chooseEngineNode(name):
    if name=='lalinferencenest':
        return LALInferenceNestNode
    if name=='lalinferenceburst':
        return LALInferenceBurstNode
    if name=='lalinferencemcmc':
        return LALInferenceMCMCNode
    if name=='lalinferencedatadump':
        return LALInferenceDataDumpNode
    if name=='bayeswavepsd':
        return BayesWavePSDNode
    return EngineNode

def get_engine_name(cp):
    name=cp.get('analysis','engine')
    if name=='random':
        engine_list=['lalinferencenest','lalinferencemcmc']
        if cp.has_option('input','gid'):
            gid=cp.get('input','gid')
            engine_number=int(''.join(i for i in gid if i.isdigit())) % 2
        else:
            engine_number=random.randint(0,1)
        return engine_list[engine_number]
    else:
        return name

def scan_timefile(timefile):
    import re
    p=re.compile('[\d.]+')
    times=[]
    timefilehandle=open(timefile,'r')
    for time in timefilehandle:
        if not p.match(time):
            continue
        if float(time) in times:
            print('Skipping duplicate time %s'%(time))
            continue
        print('Read time %s'%(time))
        times.append(float(time))
    timefilehandle.close()
    return times

def get_xml_psds(psdxml,ifos,outpath,end_time=None):
    """
    Get a psd.xml.gz file and:
    1) Reads it
    2) Checks the psd file contains all the IFO we want to analyze
    3) Writes down the PSDs into an ascii file for each IFO in psd.xml.gz. The name of the file contains the trigtime (if given) and the IFO name.
    Input:
      psdxml: psd.xml.gz file
      ifos: list of ifos used for the analysis
      outpath: path where the ascii PSD will be written to
      (end_time): trigtime for this event. Will be used a part of the PSD file name
    """
    from glue.ligolw import utils as ligolw_utils
    try:
        from lal import series as lalseries
    except ImportError:
        print("ERROR, cannot import lal.series in bppu/get_xml_psds()\n")
        raise

    out={}
    if not os.path.isdir(outpath):
        os.makedirs(outpath)
    if end_time is not None:
        time=repr(float(end_time))
    else:
        time=''
    #check we don't already have ALL the psd files #
    got_all=1
    for ifo in ifos:
        path_to_ascii_psd=os.path.join(outpath,ifo+'_psd_'+time+'.txt')
        # Check we don't already have that ascii (e.g. because we are running parallel runs of the save event
        if os.path.isfile(path_to_ascii_psd):
            got_all*=1
        else:
            got_all*=0
    if got_all==1:
        #print "Already have PSD files. Nothing to do...\n"
        for ifo in ifos:
            out[ifo]=os.path.join(outpath,ifo+'_psd_'+time+'.txt')
        return out

    # We need to convert the PSD for one or more IFOS. Open the file
    if not os.path.isfile(psdxml):
        print("ERROR: impossible to open the psd file %s. Exiting...\n"%psdxml)
        sys.exit(1)
    xmlpsd =  lalseries.read_psd_xmldoc(ligolw_utils.load_filename(psdxml,contenthandler = lalseries.PSDContentHandler))
    # Check the psd file contains all the IFOs we want to analize
    for ifo in ifos:
        if not ifo in xmlpsd:
            print("ERROR. The PSD for the ifo %s does not seem to be contained in %s\n"%(ifo,psdxml))
            sys.exit(1)
    #loop over ifos in psd xml file
    for instrument in xmlpsd.keys():
        #name of the ascii file we are going to write the PSD into
        path_to_ascii_psd=os.path.join(outpath,instrument+'_psd_'+time+'.txt')
        # Check we don't already have that ascii (e.g. because we are running parallel runs of the save event
        if os.path.isfile(path_to_ascii_psd):
            continue
        # get data for the IFO
        ifodata=xmlpsd[instrument]
        #check data is not empty
        if ifodata is None:
            continue
        # we have data. Get psd array
        data=ifodata.data
        # Fill a two columns array of (freq, psd) and save it in the ascii file
        f0=ifodata.f0
        deltaF=ifodata.deltaF

        combine=[]

        for i in np.arange(len(data.data.data)) :
            combine.append([f0+i*deltaF,data.data.data[i]])
        np.savetxt(path_to_ascii_psd,combine)
        ifo=instrument
        # set node.psds dictionary with the path to the ascii files
        out[ifo]=os.path.join(outpath,ifo+'_psd_'+time+'.txt')
    return out

def get_trigger_chirpmass(gid=None,gracedb="gracedb",service_url=None):
    from glue.ligolw import lsctables
    from glue.ligolw import ligolw
    from glue.ligolw import utils as ligolw_utils
    from ligo.gracedb.rest import GraceDb
    cwd=os.getcwd()
    if service_url is None:
        client = GraceDb()
    else:
        client = GraceDb(service_url=service_url)
    xmldoc = ligolw_utils.load_fileobj(client.files(gid, "coinc.xml"), contenthandler = lsctables.use_in(ligolw.LIGOLWContentHandler))[0]
    ligolw_utils.write_filename(xmldoc, "coinc.xml")
    coinc_events = lsctables.CoincInspiralTable.get_table(xmldoc)
    sngl_event_idx = dict((row.event_id, row) for row in lsctables.SnglInspiralTable.get_table(xmldoc))
    coinc_map = lsctables.CoincMapTable.get_table(xmldoc)
    mass1 = []
    mass2 = []
    for coinc in coinc_events:
        these_sngls = [sngl_event_idx[c.event_id] for c in coinc_map if c.coinc_event_id == coinc.coinc_event_id]
        for e in these_sngls:
            mass1.append(e.mass1)
            mass2.append(e.mass2)
    # check that trigger masses are identical in each IFO
    assert len(set(mass1)) == 1
    assert len(set(mass2)) == 1

    mchirp = (mass1[0]*mass2[0])**(3./5.) / ( (mass1[0] + mass2[0])**(1./5.) )
    os.remove("coinc.xml")

    return mchirp

def get_roq_mchirp_priors(path, roq_paths, roq_params, key, gid=None,sim_inspiral=None, service_url=None):

    ## XML and GID cannot be given at the same time
    ## sim_inspiral must already point at the right row
    mc_priors = {}

    if gid is not None and sim_inspiral is not None:
        print("Error in get_roq_mchirp_priors, cannot use both gid and sim_inspiral\n")
        sys.exit(1)

    for roq in roq_paths:
        params=os.path.join(path,roq,'params.dat')
        roq_params[roq]=np.genfromtxt(params,names=True)
        mc_priors[roq]=[float(roq_params[roq]['chirpmassmin']),float(roq_params[roq]['chirpmassmax'])]
    ordered_roq_paths=[item[0] for item in sorted(roq_params.items(), key=key)][::-1]
    # below is to construct non-overlapping mc priors for multiple roq mass-bin runs
    '''i=0
    for roq in ordered_roq_paths:
      if i>0:
        # change min, just set to the max of the previous one since we have already aligned it in the previous iteration of this loop
        #mc_priors[roq][0]+= (mc_priors[roq_lengths[i-1]][1]-mc_priors[roq][0])/2.
        mc_priors[roq][0]=mc_priors[ordered_roq_paths[i-1]][1]
      if i<len(roq_paths)-1:
        mc_priors[roq][1]-= (mc_priors[roq][1]- mc_priors[ordered_roq_paths[i+1]][0])/2.
      i+=1'''
    if gid is not None:
        trigger_mchirp = get_trigger_chirpmass(gid=gid,service_url=service_url)
    elif sim_inspiral is not None:
        trigger_mchirp = sim_inspiral.mchirp
    else:
        trigger_mchirp = None

    return mc_priors, trigger_mchirp

def get_roq_component_mass_priors(path, roq_paths, roq_params, key, gid=None,sim_inspiral=None,service_url=None):

    ## XML and GID cannot be given at the same time
    ## sim_inspiral must already point at the right row
    m1_priors = {}
    m2_priors = {}

    if gid is not None and sim_inspiral is not None:
        print("Error in get_roq_mchirp_priors, cannot use both gid and sim_inspiral\n")
        sys.exit(1)

    for roq in roq_paths:
        params=os.path.join(path,roq,'params.dat')
        roq_params[roq]=np.genfromtxt(params,names=True)
        m1_priors[roq]=[float(roq_params[roq]['mass1min']),float(roq_params[roq]['mass1max'])]
        m2_priors[roq]=[float(roq_params[roq]['mass2min']),float(roq_params[roq]['mass2max'])]

    if gid is not None:
        trigger_mchirp = get_trigger_chirpmass(gid,service_url=service_url)
    elif sim_inspiral is not None:
        trigger_mchirp = sim_inspiral.mchirp
    else:
        trigger_mchirp = None

    return m1_priors, m2_priors, trigger_mchirp

def get_roq_mass_freq_scale_factor(mc_priors, trigger_mchirp, force_flow=None):
    mc_priors_keys_list = list(mc_priors.keys())
    mc_priors_keys_int = [int(seglen[:-1]) for seglen in mc_priors_keys_list]
    roq_min = mc_priors_keys_list[np.argmin(mc_priors_keys_int)]
    roq_max = mc_priors_keys_list[np.argmax(mc_priors_keys_int)]
    mc_max = mc_priors[roq_min][1]
    mc_min = mc_priors[roq_max][0]
    scale_factor = 1.
    if force_flow == None and trigger_mchirp != None:
        if trigger_mchirp >= mc_max:
            scale_factor = 2.**(floor(trigger_mchirp/mc_max))
        if trigger_mchirp <= mc_min:
            scale_factor = (2./3.2)**(ceil(trigger_mchirp/mc_min))
    elif force_flow != None:
        scale_factor = 20./force_flow
    return scale_factor

def create_pfn_tuple(filename,protocol='file://',site='local'):
    return( (os.path.basename(filename),protocol+os.path.abspath(filename),site) )

def mchirp_from_components(m1, m2):
    return (m1*m2)**(3.0/5.0) / (m1+m2)**(1.0/5.0)

def Query_ROQ_Bounds_Type(path, roq_paths):
    # Assume that parametrization of ROQ bounds is independent of seglen; just look at first one
    import numpy as np
    roq = roq_paths[0]
    params = os.path.join(path,roq,'params.dat')
    roq_params0 = np.genfromtxt(params,names=True)
    roq_names_set = set(roq_params0.dtype.names)
    component_mass_bounds_set = set(['mass1min', 'mass1max', 'mass2min', 'mass2max'])
    chirp_mass_q_bounds_set = set(['chirpmassmin', 'chirpmassmax', 'qmin', 'qmax'])
    if roq_names_set.issuperset(component_mass_bounds_set):
        roq_bounds = 'component_mass'
    elif roq_names_set.issuperset(chirp_mass_q_bounds_set):
        roq_bounds = 'chirp_mass_q'
    else:
        print('Invalid bounds for ROQ. Ether (m1,m2) or (mc,q) bounds are supported.')
        sys.exit(1)
    return roq_bounds

class LALInferencePipelineDAG(pipeline.CondorDAG):
    def __init__(self,cp,dax=False,site='local'):
        self.subfiles=[]
        self.config=cp
        self.engine=get_engine_name(cp)
        self.EngineNode=chooseEngineNode(self.engine)
        self.site=site
        if cp.has_option('paths','basedir'):
            self.basepath=cp.get('paths','basedir')
        else:
            self.basepath=os.getcwd()
            print('No basepath specified, using current directory: %s'%(self.basepath))
        mkdirs(self.basepath)
        print("Generating LALInference DAG in {0}".format(self.basepath))
        if dax:
            os.chdir(self.basepath)
        self.posteriorpath=os.path.join(self.basepath,'posterior_samples')
        mkdirs(self.posteriorpath)
        daglogdir=cp.get('paths','daglogdir')
        mkdirs(daglogdir)
        self.daglogfile=os.path.join(daglogdir,'lalinference_pipeline-'+str(uuid.uuid1())+'.log')
        super(LALInferencePipelineDAG,self).__init__(self.daglogfile,dax=dax)
        if cp.has_option('paths','cachedir'):
            self.cachepath=cp.get('paths','cachedir')
        else:
            self.cachepath=os.path.join(self.basepath,'caches')
        mkdirs(self.cachepath)
        if cp.has_option('paths','logdir'):
            self.logpath=cp.get('paths','logdir')
        else:
            self.logpath=os.path.join(self.basepath,'log')
        mkdirs(self.logpath)
        if cp.has_option('analysis','ifos'):
            self.ifos=ast.literal_eval(cp.get('analysis','ifos'))
        else:
            self.ifos=['H1','L1','V1']
        self.segments={}
        if cp.has_option('datafind','veto-categories'):
            self.veto_categories=cp.get('datafind','veto-categories')
        else: self.veto_categories=[]
        for ifo in self.ifos:
            self.segments[ifo]=[]
        self.computeroqweightsnode={}
        self.bayeslinenode={}
        self.bayeswavepsdnode={}
        self.dq={}
        self.frtypes=ast.literal_eval(cp.get('datafind','types'))
        self.channels=ast.literal_eval(cp.get('data','channels'))
        self.use_available_data=False
        self.webdir=cp.get('paths','webdir')
        if cp.has_option('analysis','dataseed'):
            self.dataseed=cp.getint('analysis','dataseed')
        else:
            self.dataseed=None
        # Set up necessary job files.
        self.prenodes={}
        self.datafind_job = pipeline.LSCDataFindJob(self.cachepath,self.logpath,self.config,dax=self.is_dax())
        self.datafind_job.add_opt('url-type','file')
        # If running on OSG use its datafind server
        if cp.has_option('analysis','osg') and cp.getboolean('analysis','osg'):
            self.datafind_job.add_opt('server','datafind.ligo.org')
        if cp.has_option('condor','accounting_group'):
            self.datafind_job.add_condor_cmd('accounting_group',cp.get('condor','accounting_group'))
        if cp.has_option('condor','accounting_group_user'):
            self.datafind_job.add_condor_cmd('accounting_group_user',cp.get('condor','accounting_group_user'))
        self.datafind_job.set_sub_file(os.path.abspath(os.path.join(self.basepath,'datafind.sub')))
        self.preengine_job = EngineJob(self.config, os.path.join(self.basepath,'prelalinference.sub'),self.logpath,engine='lalinferencedatadump',ispreengine=True,dax=self.is_dax())
        self.preengine_job.set_grid_site('local')
        self.preengine_job.set_universe('vanilla')
        if self.config.getboolean('analysis','roq'):
            self.computeroqweights_job = ROMJob(self.config,os.path.join(self.basepath,'computeroqweights.sub'),self.logpath,dax=self.is_dax())
            self.computeroqweights_job.set_grid_site('local')
        if self.config.has_option('condor','bayesline'):
            self.bayesline_job = BayesLineJob(self.config,os.path.join(self.basepath,'bayesline.sub'),self.logpath,dax=self.is_dax())
            self.bayesline_job.set_grid_site('local')
        self.bayeswavepsd_job={}
        if self.config.has_option('condor','bayeswave'):
            for ifo in self.ifos:
                self.bayeswavepsd_job[ifo] = BayesWavePSDJob(self.config,os.path.join(self.basepath,'bayeswavepsd_%s.sub'%(ifo)),self.logpath,dax=self.is_dax())
                self.bayeswavepsd_job[ifo].set_grid_site('local')
        # Need to create a job file for each IFO combination
        self.engine_jobs={}
        ifocombos=[]
        for N in range(1,len(self.ifos)+1):
            for a in permutations(self.ifos,N):
                ifocombos.append(a)
        for ifos in ifocombos:
            self.engine_jobs[ifos] = EngineJob(self.config, os.path.join(self.basepath,'engine_%s.sub'%(reduce(lambda x,y:x+y, map(str,ifos)))),self.logpath,engine=self.engine,dax=self.is_dax(), site=site)
        self.results_page_job = ResultsPageJob(self.config,os.path.join(self.basepath,'resultspage.sub'),self.logpath,dax=self.is_dax())
        self.results_page_job.set_grid_site('local')
        self.cotest_results_page_job = ResultsPageJob(self.config,os.path.join(self.basepath,'resultspagecoherent.sub'),self.logpath,dax=self.is_dax())
        self.cotest_results_page_job.set_grid_site('local')
        if self.engine=='lalinferencemcmc':
            self.combine_job = CombineMCMCJob(self.config,os.path.join(self.basepath,'combine_files.sub'),self.logpath,dax=self.is_dax())
            self.combine_job.set_grid_site('local')
            self.merge_job = MergeJob(self.config,os.path.join(self.basepath,'merge_runs.sub'),self.logpath,dax=self.is_dax(),engine='mcmc')
            self.merge_job.set_grid_site('local')
        else:
            self.merge_job = MergeJob(self.config,os.path.join(self.basepath,'merge_runs.sub'),self.logpath,dax=self.is_dax(),engine='nest')
            self.merge_job.set_grid_site('local')
        self.coherence_test_job = CoherenceTestJob(self.config,os.path.join(self.basepath,'coherence_test.sub'),self.logpath,dax=self.is_dax())
        self.coherence_test_job.set_grid_site('local')
        self.gracedbjob = GraceDBJob(self.config,os.path.join(self.basepath,'gracedb.sub'),self.logpath,dax=self.is_dax())
        self.gracedbjob.set_grid_site('local')
        self.mapjob = SkyMapJob(cp, os.path.join(self.basepath,'skymap.sub'), self.logpath)
        self.plotmapjob = PlotSkyMapJob(cp, os.path.join(self.basepath,'plotskymap.sub'),self.logpath)
        # Process the input to build list of analyses to do
        self.events=self.setup_from_inputs()

        # Sanity checking
        if len(self.events)==0:
            print('No input events found, please check your config if you expect some events')
        self.times=[e.trig_time for e in self.events]

        # Set up the segments
        if not (self.config.has_option('input','gps-start-time') and self.config.has_option('input','gps-end-time')) and len(self.times)>0:
            (mintime,maxtime)=self.get_required_data(self.times)
            if not self.config.has_option('input','gps-start-time'):
                self.config.set('input','gps-start-time',str(int(floor(mintime))))
            if not self.config.has_option('input','gps-end-time'):
                self.config.set('input','gps-end-time',str(int(ceil(maxtime))))
        self.add_science_segments()

        # Save the final configuration that is being used
        # first to the run dir
        conffilename=os.path.join(self.basepath,'config.ini')
        with open(conffilename,'w') as conffile:
            self.config.write(conffile)
        if self.config.has_option('paths','webdir'):
            mkdirs(self.config.get('paths','webdir'))
            with open(os.path.join(self.config.get('paths','webdir'),'config.ini'),'w') as conffile:
                self.config.write(conffile)

        # Generate the DAG according to the config given
        for event in self.events: self.add_full_analysis(event)
        if self.config.has_option('analysis','upload-to-gracedb'):
            if self.config.getboolean('analysis','upload-to-gracedb'):
                self.add_gracedb_FITSskymap_upload(self.events[0],engine=self.engine)
        self.dagfilename="lalinference_%s-%s"%(self.config.get('input','gps-start-time'),self.config.get('input','gps-end-time'))
        self.set_dag_file(os.path.join(self.basepath,self.dagfilename))
        if self.is_dax():
            self.set_dax_file(self.dagfilename)

    def add_full_analysis(self,event):
        if self.engine=='lalinferencenest' or  self.engine=='lalinferenceburst':
            result=self.add_full_analysis_lalinferencenest(event)
        elif self.engine=='lalinferencemcmc':
            result=self.add_full_analysis_lalinferencemcmc(event)
        else:
            raise Exception('Unknown engine {0}'.format(self.engine))
        return result

    def create_frame_pfn_file(self):
        """
        Create a pegasus cache file name, uses inspiralutils
        """
        import inspiralutils as iu
        gpsstart=self.config.get('input','gps-start-time')
        gpsend=self.config.get('input','gps-end-time')
        pfnfile=iu.create_frame_pfn_file(self.frtypes,gpsstart,gpsend)
        return pfnfile

    def get_required_data(self,times):
        """
        Calculate the data that will be needed to process all events
        """
        #psdlength = self.config.getint('input','max-psd-length')
        padding=self.config.getint('input','padding')
        if self.config.has_option('engine','seglen') or self.config.has_option('lalinference','seglen'):
            if self.config.has_option('engine','seglen'):
                seglen = int(np.ceil(self.config.getfloat('engine','seglen')))
            if self.config.has_option('lalinference','seglen'):
                seglen = self.config.getint('lalinference','seglen')

            if os.path.isfile(os.path.join(self.basepath,'psd.xml.gz')) or self.config.has_option('condor','bayesline') or self.config.has_option('condor','bayeswave'):
                psdlength = 0
                padding = 0
                self.config.set('input','padding',str(padding))
                if self.config.has_option('condor','bayeswave'):
                    if (np.log2(seglen)%1):
                        seglen = np.power(2., np.ceil(np.log2(seglen)))

            else:
                psdlength = 32*seglen
        else:
            seglen = max(e.duration for e in self.events)
            if os.path.isfile(os.path.join(self.basepath,'psd.xml.gz')) or self.config.has_option('condor','bayesline') or self.config.has_option('condor','bayeswave'):
                psdlength = 0
                padding = 0
                self.config.set('input','padding',str(padding))
                if self.config.has_option('condor','bayeswave'):
                    if (np.log2(seglen)%1):
                        seglen = np.power(2., np.ceil(np.log2(seglen)))
            else:
                psdlength = 32*seglen
        # Assume that the data interval is (end_time - seglen -padding , end_time + psdlength +padding )
        # -> change to (trig_time - seglen - padding - psdlength + 2 , trig_time + padding + 2) to estimate the psd before the trigger for online follow-up.
        # Also require padding before start time
        return (min(times)-padding-seglen-psdlength+2,max(times)+padding+2)

    def setup_from_times(self,times):
        """
        Generate a DAG from a list of times
        """
        for time in self.times:
            self.add_full_analysis(Event(trig_time=time))

    def select_events(self):
        """
        Read events from the config parser. Understands both ranges and comma separated events, or combinations
        eg. events=[0,1,5:10,21] adds to the analysis the events: 0,1,5,6,7,8,9,10 and 21
        """
        events=[]
        times=[]
        raw_events=self.config.get('input','events').replace('[','').replace(']','').split(',')
        for raw_event in raw_events:
            if ':' in raw_event:
                limits=raw_event.split(':')
                if len(limits) != 2:
                    print("Error: in event config option; ':' must separate two numbers.")
                    exit(0)
                low=int(limits[0])
                high=int(limits[1])
                if low>high:
                    events.extend(range(int(high),int(low)+1))
                elif high>low:
                    events.extend(range(int(low),int(high)+1))
            else:
                events.append(int(raw_event))
        return events

    def setup_from_inputs(self):
        """
        Scan the list of inputs, i.e.
        gps-time-file, injection-file, sngl-inspiral-file, coinc-inspiral-file, pipedown-database
        in the [input] section of the ini file.
        And process the events found therein
        """
        events=[]
        gpsstart=None
        gpsend=None
        if self.config.has_option('input','gps-start-time'):
            gpsstart=self.config.getfloat('input','gps-start-time')
        if self.config.has_option('input','gps-end-time'):
            gpsend=self.config.getfloat('input','gps-end-time')
        inputnames=['gps-time-file','burst-injection-file','injection-file','sngl-inspiral-file','coinc-inspiral-file','pipedown-db','gid','gstlal-db']
        ReadInputFromList=sum([ 1 if self.config.has_option('input',name) else 0 for name in inputnames])
        # If no input events given, just return an empty list (e.g. for PP pipeline)
        if ReadInputFromList!=1 and (gpsstart is None or gpsend is None):
            return []
        # Review: Clean up this section
        if self.config.has_option('input','events'):
            selected_events=self.config.get('input','events')
            print('Selected events %s'%(str(selected_events)))

            if selected_events=='all':
                selected_events=None
            else:
                selected_events=self.select_events()
        else:
            selected_events=None

        if(self.config.has_option('engine','correlatedGaussianLikelihood') or
           self.config.has_option('engine','bimodalGaussianLikelihood') or
           self.config.has_option('engine','rosenbrockLikelihood')):
            analytic_test = True
        else:
            analytic_test = False

        # No input file given, analyse the entire time stretch between gpsstart and gpsend
        if self.config.has_option('input','analyse-all-time') and self.config.getboolean('input','analyse-all-time')==True:
            print('Setting up for analysis of continuous time stretch %f - %f'%(gpsstart,gpsend))
            if self.config.has_option('engine','seglen'):
                seglen=self.config.getfloat('engine','seglen')
            else:
                print('ERROR: seglen must be specified in [engine] section when running without input file')
                sys.exit(1)
            if(self.config.has_option('input','segment-overlap')):
                overlap=self.config.getfloat('input','segment-overlap')
            else:
                overlap=32.;
            if(overlap>seglen):
                print('ERROR: segment-overlap is greater than seglen')
                sys.exit(1)
            # Now divide gpsstart - gpsend into jobs of seglen - overlap length
            t=gpsstart
            events=[]
            while(t<gpsend):
                ev=Event(trig_time=t+seglen-2)
                ev.set_engine_option('segment-start',str(t-overlap))
                if not analytic_test:
                    ev.set_engine_option('time-min',str(t))
                tMax=t + seglen - overlap
                if tMax>=gpsend:
                    tMax=gpsend
                if not analytic_test:
                    ev.set_engine_option('time-max',str(tMax))
                events.append(ev)
                t=tMax
            return events

        # ASCII list of GPS times
        if self.config.has_option('input','gps-time-file'):
            times=scan_timefile(self.config.get('input','gps-time-file'))
            if self.config.has_option('input','timeslides-ascii'):
            # The timeslides-ascii files contains one row per trigtime, and a column per IFO
            # Note: the IFO order is the same given in the ifos section of the [analysis] tag
                print("Reading timeslides from ascii file. Columns order is understood as follow:")
                for this_ifo,ifo in enumerate(self.ifos):
                    print("Column %d"%this_ifo + "= %s "%(ifo))
                dest=self.config.get('input','timeslides-ascii')
                if not os.path.isfile(dest):
                    print("ERROR the ascii file %s containing the timeslides does not exist\n"%dest)
                    exit(1)
                else:
                    from numpy import loadtxt
                    data=loadtxt(dest).reshape(-1,len(self.ifos))
                    if len(self.ifos)!= len(data[0,:]):
                        print("ERROR: ascii timeslide file must contain a column for each IFO used in the analysis!\n")
                        exit(1)
                    if len(times)!=len(data[:,0]):
                        print('ERROR: ascii timeslide must contain a row for each trigtime. Exiting...\n')
                        exit(1)
                    timeslides={}
                    for this_time,time in enumerate(times):
                        timeslides[this_time]={}
                        for this_ifo,ifo in enumerate(self.ifos):
                            timeslides[this_time][ifo]=data[this_time,this_ifo]
                events=[Event(trig_time=time,timeslide_dict=timeslides[i_time]) for i_time,time in enumerate(times)]
            else:
                events=[Event(trig_time=time) for time in times]
        # Siminspiral Table
        if self.config.has_option('input','injection-file'):
            from glue.ligolw import ligolw, lsctables, utils
            injTable = lsctables.SimInspiralTable.get_table(
                              utils.load_filename(self.config.get('input','injection-file'),
                                                  contenthandler=lsctables.use_in(ligolw.LIGOLWContentHandler)) )
            events=[Event(SimInspiral=inj) for inj in injTable]
            self.add_pfn_cache([create_pfn_tuple(self.config.get('input','injection-file'))])
        # SimBurst Table
        if self.config.has_option('input','burst-injection-file'):
            injfile=self.config.get('input','burst-injection-file')
            injTable=lsctables.SimBurstTable.get_table(utils.load_filename(injfile,contenthandler = lsctables.use_in(LIGOLWContentHandler)))
            events=[Event(SimBurst=inj) for inj in injTable]
            self.add_pfn_cache([create_pfn_tuple(self.config.get('input','burst-injection-file'))])
        # SnglInspiral Table
        if self.config.has_option('input','sngl-inspiral-file'):
            trigTable = lsctables.SnglInspiralTable.get_table(utils.load_filename(injfile, contenthandler = lsctables.use_in(ligolw.LIGOLWContentHandler)))
            events=[Event(SnglInspiral=trig) for trig in trigTable]
            self.add_pfn_cache([create_pfn_tuple(self.config.get('input','sngl-inspiral-file'))])
        if self.config.has_option('input','coinc-inspiral-file'):
            coincTable = lsctables.CoincInspiralTable.get_table(utils.load_filename(injfile, contenthandler = lsctables.use_in(ligolw.LIGOLWContentHandler)))
            events = [Event(CoincInspiral=coinc) for coinc in coincTable]
            self.add_pfn_cache([create_pfn_tuple(self.config.get('input','coinc-inspiral-file'))])
        # LVAlert CoincInspiral Table
        if self.config.has_option('input','gid'):
            gid=self.config.get('input','gid')
            flow=20.0
            if self.config.has_option('lalinference','flow'):
                flow=min(ast.literal_eval(self.config.get('lalinference','flow')).values())
            downloadgracedbpsd=True
            if self.config.has_option('input','ignore-gracedb-psd'):
                if self.config.getboolean('input','ignore-gracedb-psd'):
                    downloadgracedbpsd=False
            threshold_snr = None
            if not self.config.has_option('engine','distance-max') and self.config.has_option('input','threshold-snr'):
                threshold_snr=self.config.getfloat('input','threshold-snr')
            service_url = None
            if self.config.has_option('analysis','service-url'):
                service_url = self.config.get('analysis','service-url')
            events = readLValert(gid=gid,flow=flow,gracedb=self.config.get('condor','gracedb'),basepath=self.basepath,downloadpsd=downloadgracedbpsd,
                                 threshold_snr=threshold_snr,roq=self.config.getboolean('analysis','roq'),service_url=service_url)
        else: gid=None
        # pipedown-database
        if self.config.has_option('input','gstlal-db'):
            queryfunc=get_zerolag_lloid
            dbname=self.config.get('input','gstlal-db')
        elif self.config.has_option('input','pipedown-db'):
            queryfunc=get_zerolag_pipedown
            dbname=self.config.get('input','pipedown-db')
        else: dbname=None
        if dbname:
            db_connection = open_pipedown_database(dbname,None)[0]
            # Timeslides
            if self.config.has_option('input','time-slide-dump'):
                timeslidedump=self.config.get('input','time-slide-dump')
            else:
                timeslidedump=None
            if self.config.has_option('input','min-cfar'):
                mincfar=self.config.getfloat('input','min-cfar')
            else:
                mincfar=-1
            if self.config.has_option('input','max-cfar'):
                maxcfar=self.config.getfloat('input','max-cfar')
            else:
                maxcfar=-1
            if self.config.get('input','timeslides').lower()=='true':
                events=get_timeslides_pipedown(db_connection, gpsstart=gpsstart, gpsend=gpsend,dumpfile=timeslidedump,max_cfar=maxcfar)
            else:
                events=queryfunc(db_connection, gpsstart=gpsstart, gpsend=gpsend, dumpfile=timeslidedump,max_cfar=maxcfar,min_cfar=mincfar)
        if(selected_events is not None):
            used_events=[]
            for i in selected_events:
                e=events[i]
                e.event_id=i
                used_events.append(e)
            events=used_events
        if gpsstart is not None:
            events = filter(lambda e: not e.trig_time<gpsstart, events)
        if gpsend is not None:
            events = filter(lambda e: not e.trig_time>gpsend, events)
        return events

    def add_full_analysis_lalinferencenest(self,event):
        """
        Generate an end-to-end analysis of a given event (Event class)
        For LALinferenceNest code. Uses parallel runs if specified
        """
        evstring=str(event.event_id)
        if event.trig_time is not None:
            evstring=str(event.trig_time)+'-'+str(event.event_id)
        if self.config.has_option('analysis','nparallel'):
            Npar=self.config.getint('analysis','nparallel')
        else:
            Npar=4
        # Set up the parallel engine nodes
        enginenodes=[]
        bwpsdnodes={}
        for i in range(Npar):
            n,bwpsdnodes=self.add_engine_node(event,bwpsdnodes)
            if n is not None:
                if i>0:
                    n.add_var_arg('--dont-dump-extras')
                enginenodes.append(n)
        if len(enginenodes)==0:
            return False
        myifos=enginenodes[0].get_ifos()
        # Merge the results together
        pagedir=os.path.join(self.webdir,evstring,myifos)
        #pagedir=os.path.join(self.basepath,evstring,myifos)
        mkdirs(pagedir)
        mergenode=MergeNode(self.merge_job,parents=enginenodes,engine='nest')
        mergenode.set_pos_output_file(os.path.join(self.posteriorpath,'posterior_%s_%s.hdf5'%(myifos,evstring)))
        self.add_node(mergenode)
        # Call finalize to build final list of available data
        enginenodes[0].finalize()
        enginenodes[0].set_psd_files()
        enginenodes[0].set_snr_file()
        if self.config.getboolean('analysis','coherence-test') and len(enginenodes[0].ifos)>1:
            if self.site!='local':
                zipfilename='postproc_'+evstring+'.tar.gz'
            else:
                zipfilename=None
            respagenode=self.add_results_page_node(resjob=self.cotest_results_page_job,outdir=pagedir,parent=mergenode,gzip_output=zipfilename,ifos=enginenodes[0].ifos)
            respagenode.set_psd_files(enginenodes[0].get_psd_files())
            respagenode.set_snr_file(enginenodes[0].get_snr_file())
            if os.path.exists(self.basepath+'/coinc.xml'):
                respagenode.set_coinc_file(self.basepath+'/coinc.xml',self.config.get('input','gid'))
            mkdirs(os.path.join(self.basepath,'coherence_test'))
            par_mergenodes=[]
            for ifo in enginenodes[0].ifos:
                co_merge_job = MergeJob(self.config,os.path.join(self.basepath,'merge_runs_%s.sub'%(ifo)),self.logpath,dax=self.is_dax(),engine='nest')
                co_merge_job.set_grid_site('local')
                cotest_nodes=[]
                for i in range(Npar):
                    cot_node,bwpsdnodes=self.add_engine_node(event,bwpsdnodes,ifos=[ifo],co_test=True)
                    if cot_node is not None:
                        if i>0:
                            cot_node.add_var_arg('--dont-dump-extras')
                        cotest_nodes.append(cot_node)
                if len(cotest_nodes)==0:
                    return False
                for co in cotest_nodes:
                    co.set_psdstart(enginenodes[0].GPSstart)
                    co.set_psdlength(enginenodes[0].psdlength)
                    if co!=cotest_nodes[0]:
                        co.add_var_arg('--dont-dump-extras')
                    else:
                        co.set_psd_files()
                        co.set_snr_file()
                pmergenode=MergeNode(co_merge_job,parents=cotest_nodes,engine='nest')
                pmergenode.set_pos_output_file(os.path.join(self.posteriorpath,'posterior_%s_%s.hdf5'%(ifo,evstring)))
                self.add_node(pmergenode)
                par_mergenodes.append(pmergenode)
                presultsdir=os.path.join(pagedir,ifo)
                mkdirs(presultsdir)
                if self.site!='local':
                    pzipfilename='postproc_'+evstring+'_'+ifo+'.tar.gz'
                else:
                    pzipfilename=None
                subresnode=self.add_results_page_node(outdir=presultsdir,parent=pmergenode, gzip_output=pzipfilename,ifos=ifo)
                subresnode.set_psd_files(cotest_nodes[0].get_psd_files())
                subresnode.set_snr_file(cotest_nodes[0].get_snr_file())
                if os.path.exists(self.basepath+'/coinc.xml'):
                    subresnode.set_coinc_file(self.basepath+'/coinc.xml',self.config.get('input','gid'))
                if self.config.has_option('input','injection-file') and event.event_id is not None:
                    subresnode.set_injection(self.config.get('input','injection-file'),event.event_id)
                elif self.config.has_option('input','burst-injection-file') and event.event_id is not None:
                    subresnode.set_injection(self.config.get('input','burst-injection-file'),event.event_id)
            coherence_node=CoherenceTestNode(self.coherence_test_job,outfile=os.path.join(self.basepath,'coherence_test','coherence_test_%s_%s.dat'%(myifos,evstring)))
            coherence_node.add_coherent_parent(mergenode)
            for parmergenode in par_mergenodes:
                coherence_node.add_incoherent_parent(parmergenode)
            self.add_node(coherence_node)
            respagenode.add_parent(coherence_node)
            respagenode.set_bayes_coherent_incoherent(coherence_node.get_output_files()[0])
        else:
            if self.site!='local':
                zipfilename='postproc_'+evstring+'.tar.gz'
            else:
                zipfilename=None
            # Note: Disabled gzip_output for now. Possibly need it for future Pegasus use
            respagenode=self.add_results_page_node(outdir=pagedir,parent=mergenode,gzip_output=None,ifos=enginenodes[0].ifos)
            respagenode.set_psd_files(enginenodes[0].get_psd_files())
            respagenode.set_snr_file(enginenodes[0].get_snr_file())
            if os.path.exists(self.basepath+'/coinc.xml'):
                respagenode.set_coinc_file(self.basepath+'/coinc.xml',self.config.get('input','gid'))
        if self.config.has_option('input','injection-file') and event.event_id is not None:
            respagenode.set_injection(self.config.get('input','injection-file'),event.event_id)
        elif self.config.has_option('input','burst-injection-file') and event.event_id is not None:
            respagenode.set_injection(self.config.get('input','burst-injection-file'),event.event_id)

        if self.config.has_option('analysis','upload-to-gracedb'):
            if self.config.getboolean('analysis','upload-to-gracedb') and event.GID is not None:
                self.add_gracedb_start_node(event.GID,'LALInference',[sciseg.get_df_node() for sciseg in enginenodes[0].scisegs.values()])
                self.add_gracedb_log_node(respagenode,event.GID)
            elif self.config.has_option('analysis','ugid'):
                # LIB will want to upload info to gracedb but if we pass the gid in the usual way the pipeline
                # will try to pull inspiral-only XML tables from the gdb page, failing.
                # To avoid that, LIB will read the gracedDB id to upload info to as an ugid=ID option
                # in the analysis section.
                ugid=self.config.get('analysis','ugid')
                self.add_gracedb_start_node(ugid,'LIB',[sciseg.get_df_node() for sciseg in enginenodes[0].scisegs.values()])
                self.add_gracedb_log_node(respagenode,ugid)
        if self.config.has_option('condor','ligo-skymap-plot') and self.config.has_option('condor','ligo-skymap-from-samples'):
            if self.engine=='lalinferenceburst': prefix='LIB'
            else: prefix='LALInference'
            mapnode = SkyMapNode(self.mapjob, posfile = mergenode.get_pos_file(), parent=mergenode,
                    prefix= prefix, outdir=pagedir)
            plotmapnode = PlotSkyMapNode(self.plotmapjob, parent=mapnode, inputfits = mapnode.outfits, output=os.path.join(pagedir,'skymap.png'))
            self.add_node(mapnode)
            self.add_node(plotmapnode)
        return True

    def add_full_analysis_lalinferencemcmc(self,event):
        """
        Generate an end-to-end analysis of a given event
        For LALInferenceMCMC.
        """
        evstring=str(event.event_id)
        if event.trig_time is not None:
            evstring=str(event.trig_time)+'-'+str(event.event_id)
        if self.config.has_option('analysis','nparallel'):
            Npar=self.config.getint('analysis','nparallel')
        else:
            Npar=2
        enginenodes=[]
        bwpsdnodes={}
        for i in range(Npar):
            n,bwpsdnodes=self.add_engine_node(event,bwpsdnodes)
            if n is not None:
                if i>0:
                    n.add_var_arg('--dont-dump-extras')
                enginenodes.append(n)
        if len(enginenodes)==0:
            return False
        myifos=enginenodes[0].get_ifos()
        enginenodes[0].set_psd_files()
        enginenodes[0].set_snr_file()
        pagedir=os.path.join(self.webdir,evstring,myifos)
        mkdirs(pagedir)
        combinenodes=[]
        for i in range(Npar):
            combinenodes.append(CombineMCMCNode(self.combine_job,parents=[enginenodes[i]]))
            input_file = combinenodes[i].get_parent_posfile(enginenodes[i])
            input_file_split_index = input_file.find('lalinferencemcmc-')
            combinenodes[i].set_pos_output_file(input_file[:input_file_split_index]+'combine_'+input_file[input_file_split_index:])
            combinenodes[i].add_var_arg(input_file)
            number_of_mpi_jobs = self.config.getint('mpi','mpi_task_count')
            for j in range(1,number_of_mpi_jobs):
                combinenodes[i].add_var_arg(input_file+".%02d" % j)
            self.add_node(combinenodes[i])
        mergenode=MergeNode(self.merge_job,parents=combinenodes,engine='mcmc')
        mergenode.set_pos_output_file(os.path.join(self.posteriorpath,'posterior_%s_%s.hdf5'%(myifos,evstring)))
        if self.config.has_option('resultspage','deltaLogP'):
            mergenode.add_var_arg('--deltaLogP '+str(self.config.getfloat('resultspage','deltaLogP')))
        if self.config.has_option('resultspage','downsample'):
            mergenode.add_var_arg('--downsample '+str(self.config.getint('resultspage','downsample')))
        if self.config.has_option('resultspage','fixedBurnin'):
            mergenode.add_var_arg('--fixedBurnin '+str(self.config.getint('resultspage','fixedBurnin')))
        self.add_node(mergenode)
        respagenode=self.add_results_page_node(outdir=pagedir,parent=mergenode,gzip_output=None,ifos=enginenodes[0].ifos)
        respagenode.set_psd_files(enginenodes[0].get_psd_files())
        respagenode.set_snr_file(enginenodes[0].get_snr_file())
        if os.path.exists(self.basepath+'/coinc.xml'):
            respagenode.set_coinc_file(self.basepath+'/coinc.xml',self.config.get('input','gid'))
        if self.config.has_option('input','injection-file') and event.event_id is not None:
            respagenode.set_injection(self.config.get('input','injection-file'),event.event_id)
        if self.config.has_option('input','burst-injection-file') and event.event_id is not None:
            respagenode.set_injection(self.config.get('input','burst-injection-file'),event.event_id)
        if event.GID is not None:
            if self.config.has_option('analysis','upload-to-gracedb'):
                if self.config.getboolean('analysis','upload-to-gracedb'):
                    self.add_gracedb_start_node(event.GID,'LALInference',[sciseg.get_df_node() for sciseg in enginenodes[0].scisegs.values()])
                    self.add_gracedb_log_node(respagenode,event.GID)
        if self.config.has_option('condor','ligo-skymap-plot') and self.config.has_option('condor','ligo-skymap-from-samples'):
            mapnode = SkyMapNode(self.mapjob, posfile = mergenode.get_pos_file(), parent=mergenode,
                    prefix= 'LALInference', outdir=pagedir)
            plotmapnode = PlotSkyMapNode(self.plotmapjob, parent=mapnode, inputfits = mapnode.outfits, output=os.path.join(pagedir,'skymap.png'))
            self.add_node(mapnode)
            self.add_node(plotmapnode)

    def add_science_segments(self):
        # Query the segment database for science segments and
        # add them to the pool of segments
        if self.config.has_option('input','ignore-science-segments'):
            if self.config.getboolean('input','ignore-science-segments'):
                start=self.config.getfloat('input','gps-start-time')
                end=self.config.getfloat('input','gps-end-time')
                i=0
                for ifo in self.ifos:
                    sciseg=pipeline.ScienceSegment((i,start,end,end-start))
                    df_node=self.get_datafind_node(ifo,self.frtypes[ifo],int(sciseg.start()),int(sciseg.end()))
                    sciseg.set_df_node(df_node)
                    self.segments[ifo].append(sciseg)
                    i+=1
                return
        # Look up science segments as required
        segmentdir=os.path.join(self.basepath,'segments')
        mkdirs(segmentdir)
        curdir=os.getcwd()
        os.chdir(segmentdir)
        for ifo in self.ifos:
            (segFileName,dqVetoes)=inspiralutils.findSegmentsToAnalyze(self.config, ifo, self.veto_categories, generate_segments=True,\
              use_available_data=self.use_available_data , data_quality_vetoes=False)
            self.dqVetoes=dqVetoes
            segfile=open(segFileName)
            segs=segmentsUtils.fromsegwizard(segfile)
            segs.coalesce()
            segfile.close()
            for seg in segs:
                sciseg=pipeline.ScienceSegment((segs.index(seg),seg[0],seg[1],seg[1]-seg[0]))
                df_node=self.get_datafind_node(ifo,self.frtypes[ifo],int(sciseg.start()),int(sciseg.end()))
                sciseg.set_df_node(df_node)
                self.segments[ifo].append(sciseg)
        os.chdir(curdir)

    def get_datafind_node(self,ifo,frtype,gpsstart,gpsend):
        node=pipeline.LSCDataFindNode(self.datafind_job)
        node.set_observatory(ifo[0])
        node.set_type(frtype)
        node.set_start(gpsstart)
        node.set_end(gpsend)
        #self.add_node(node)
        return node

    def add_engine_node(self,event,bwpsd={},ifos=None,co_test=False,extra_options=None):
        """
        Add an engine node to the dag. Will find the appropriate cache files automatically.
        Will determine the data to be read and the output file.
        Will use all IFOs known to the DAG, unless otherwise specified as a list of strings
        """
        if ifos is None and len(event.ifos)>0:
            ifos=event.ifos
        if ifos is None:
            ifos=self.ifos
        end_time=event.trig_time
        if self.config.has_option('lalinference','seglen'):
            seglen=self.config.getfloat('lalinference','seglen')
        elif self.config.has_option('engine','seglen'):
            seglen=self.config.getfloat('engine','seglen')
        else:
            seglen=event.duration
        segstart=end_time+2-seglen
        segend=segstart+seglen
        myifos=set([])
        for ifo in ifos:
            for seg in self.segments[ifo]:
                if segstart >= seg.start() and segend <= seg.end():
                    myifos.add(ifo)
        ifos=myifos
        if len(ifos)==0:
            print('No data found for time %f - %f, skipping'%(segstart,segend))
            return

        computeroqweightsnode={}
        bayeslinenode={}
        bayeswavepsdnode=bwpsd
        prenode=LALInferenceDataDumpNode(self.preengine_job)
        node=self.EngineNode(self.engine_jobs[tuple(ifos)])
        roqeventpath=os.path.join(self.preengine_job.roqpath,str(event.event_id)+'/')
        if self.config.has_option('condor','bayesline') or self.config.has_option('condor','bayeswave') or self.config.getboolean('analysis','roq'):
            mkdirs(roqeventpath)
        node.set_trig_time(end_time)
        prenode.set_trig_time(end_time)
        randomseed=random.randint(1,2**31)
        node.set_seed(randomseed)
        prenode.set_seed(randomseed)
        original_srate=0
        srate=0
        if event.srate:
            original_srate=event.srate
        if self.config.has_option('lalinference','srate'):
            original_srate=self.config.getfloat('lalinference','srate')
        elif self.config.has_option('engine','srate'):
            original_srate=self.config.getfloat('engine','srate')
        if (np.log2(original_srate)%1):
            print('The srate given,'+str(original_srate)+' Hz, is not a power of 2.')
            print('For data handling purposes, the srate used will be '+str(np.power(2., np.ceil(np.log2(original_srate))))+' Hz')
            print('The inner products will still however be integrated up to '+str(original_srate/2.)+' Hz')
            srate = np.power(2., np.ceil(np.log2(original_srate)))
        else:
            srate = original_srate
        if srate is not 0:
            node.set_srate(int(np.ceil(srate)))
            prenode.set_srate(int(np.ceil(srate)))
        if original_srate is not 0:
            for ifo in ifos:
                node.fhighs[ifo]=str(original_srate/2.-1./seglen)
                prenode.fhighs[ifo]=str(original_srate/2.-1./seglen)
            node.set_srate(srate)
            prenode.set_srate(srate)
        if event.trigSNR:
            node.set_trigSNR(event.trigSNR)
        if event.horizon_distance:
            node.set_horizon_distance(event.horizon_distance)
        if self.dataseed:
            node.set_dataseed(self.dataseed+event.event_id)
            prenode.set_dataseed(self.dataseed+event.event_id)
        gotdata=0
        for ifo in ifos:
            if ifo in event.timeslides:
                slide=event.timeslides[ifo]
            else:
                slide=0
            for seg in self.segments[ifo]:
                if segstart >= seg.start() and segend <= seg.end():
                    if not self.config.has_option('lalinference','fake-cache'):
                        if self.config.has_option('condor','bayesline') or self.config.getboolean('analysis','roq'):
                            prenode.add_ifo_data(ifo,seg,self.channels[ifo],timeslide=slide)
                        gotdata+=node.add_ifo_data(ifo,seg,self.channels[ifo],timeslide=slide)
                    else:
                        fakecachefiles=ast.literal_eval(self.config.get('lalinference','fake-cache'))
                        if self.config.has_option('condor','bayesline') or self.config.getboolean('analysis','roq'):
                            prenode.add_fake_ifo_data(ifo,seg,fakecachefiles[ifo],self.channels[ifo],timeslide=slide)
                        gotdata+=node.add_fake_ifo_data(ifo,seg,fakecachefiles[ifo],self.channels[ifo],timeslide=slide)
        if self.config.has_option('lalinference','psd-xmlfile'):
            psdpath=os.path.realpath(self.config.get('lalinference','psd-xmlfile'))
            node.psds=get_xml_psds(psdpath,ifos,os.path.join(self.basepath,'PSDs'),end_time=end_time)
            prenode.psds=get_xml_psds(psdpath,ifos,os.path.join(self.basepath,'PSDs'),end_time=end_time)
            if len(ifos)==0:
                node.ifos=node.cachefiles.keys()
                prenode.ifos=prenode.cachefiles.keys()
            else:
                node.ifos=ifos
                prenode.ifos=ifos
            gotdata=1
        if self.config.has_option('input','gid'):
            if os.path.isfile(os.path.join(self.basepath,'psd.xml.gz')):
                psdpath=os.path.join(self.basepath,'psd.xml.gz')
                node.psds=get_xml_psds(psdpath,ifos,os.path.join(self.basepath,'PSDs'),end_time=None)
                prenode.psds=get_xml_psds(psdpath,ifos,os.path.join(self.basepath,'PSDs'),end_time=None)
        for ifo in ifos:
            prenode.flows[ifo]=str(20.0)
        if self.config.has_option('lalinference','flow'):
            node.flows=ast.literal_eval(self.config.get('lalinference','flow'))
            prenode.flows=ast.literal_eval(self.config.get('lalinference','flow'))
        if event.fhigh:
            for ifo in ifos:
                node.fhighs[ifo]=str(event.fhigh)
                prenode.fhighs[ifo]=str(event.fhigh)
        if self.config.has_option('lalinference','fhigh'):
            node.fhighs=ast.literal_eval(self.config.get('lalinference','fhigh'))
            prenode.fhighs=ast.literal_eval(self.config.get('lalinference','fhigh'))
        prenode.set_max_psdlength(self.config.getint('input','max-psd-length'))
        prenode.set_padding(self.config.getint('input','padding'))
        #prenode[ifo].set_output_file('/dev/null')
        prenode.add_var_arg('--outfile '+roqeventpath+'data-dump')
        prenode.add_var_arg('--data-dump')
        if self.config.has_option('lalinference','seglen'):
            p_seglen=self.config.getfloat('lalinference','seglen')
        elif self.config.has_option('engine','seglen'):
            p_seglen=self.config.getfloat('engine','seglen')
        else:
            p_seglen=event.duration
        prenode.set_seglen(p_seglen)
        if self.config.has_option('condor','bayeswave'):
            prenode.set_psdlength(p_seglen)
            if self.config.has_option('lalinference','seglen'):
                bw_seglen = self.config.getfloat('lalinference','seglen')
            elif self.config.has_option('engine','seglen'):
                bw_seglen = self.config.getfloat('engine','seglen')
            else:
                bw_seglen = event.duration
            if (np.log2(bw_seglen)%1):
                print('BayesWave only supports seglengths which are powers of 2, you have specified a seglength of '+str(bw_seglen)+' seconds.')
                print('Instead, a seglenth of '+str(np.power(2., np.ceil(np.log2(bw_seglen))))+'s will be used for the BayesWave PSD estimation.')
                print('The main LALInference job will stil use seglength '+str(bw_seglen)+' seconds.')
                bw_seglen = np.power(2., np.ceil(np.log2(bw_seglen)))
        if self.config.has_option('condor','bayeswave'):
            if (np.log2(srate)%1):
                print('BayesWave only supports srates which are powers of 2, you have specified a srate of '+str(srate)+' Hertz.')
                print('Instead, a srate of '+str(np.power(2., np.ceil(np.log2(srate))))+'Hz will be used for the BayesWave PSD estimation.')
                print('The main LALInference job will stil use srate '+str(srate)+' Hertz.')
                bw_srate = np.power(2., np.ceil(np.log2(srate)))
            else:
                bw_srate = srate
        # Add the nodes it depends on
        for ifokey, seg in node.scisegs.items():
            dfnode=seg.get_df_node()

            if 1==1:
                if self.config.has_option('condor','bayeswave') and not co_test:
                    for ifo in ifos:
                        if ifo not in bayeswavepsdnode:
                            bayeswavepsdnode[ifo]=self.add_bayeswavepsd_node(ifo)
                            bayeswavepsdnode[ifo].add_var_arg('--bayesLine')
                            bayeswavepsdnode[ifo].add_var_arg('--cleanOnly')
                            bayeswavepsdnode[ifo].add_var_arg('--checkpoint')
                            bayeswavepsdnode[ifo].add_var_arg('--outputDir '+roqeventpath)
                            bayeswavepsdnode[ifo].add_var_arg('--runName BayesWave_PSD_'+ifo)
                            bayeswavepsdnode[ifo].add_output_file(os.path.join(roqeventpath,'BayesWave_PSD_'+ifo+'_IFO0_psd.dat'))
                            bayeswavepsdnode[ifo].set_trig_time(end_time)
                            bayeswavepsdnode[ifo].set_seglen(bw_seglen)
                            bayeswavepsdnode[ifo].set_psdlength(bw_seglen)
                            bayeswavepsdnode[ifo].set_srate(bw_srate)
                            if event.timeslides.has_key(ifo):
                                slide=event.timeslides[ifo]
                            else:
                                slide=0
                            for seg in self.segments[ifo]:
                                if segstart >= seg.start() and segend <= seg.end():
                                    if not self.config.has_option('lalinference','fake-cache'):
                                        bayeswavepsdnode[ifo].add_ifo_data(ifo,seg,self.channels[ifo],timeslide=slide)
                                    else:
                                        fakecachefiles=ast.literal_eval(self.config.get('lalinference','fake-cache'))
                                        bayeswavepsdnode[ifo].add_fake_ifo_data(ifo,seg,fakecachefiles[ifo],self.channels[ifo],timeslide=slide)
                            if self.config.has_option('lalinference','flow'):
                                bayeswavepsdnode[ifo].flows[ifo]=np.power(2,np.floor(np.log2(ast.literal_eval(self.config.get('lalinference','flow'))[ifo])))
                                print('BayesWave requires f_low being a power of 2, therefore f_low for '+ifo+' has been changed from '+str(ast.literal_eval(self.config.get('lalinference','flow'))[ifo])+' to '+str(np.power(2,np.floor(np.log2(ast.literal_eval(self.config.get('lalinference','flow'))[ifo]))))+' Hz (for the BayesWave job only, in the main LALInference jobs f_low will still be '+str(ast.literal_eval(self.config.get('lalinference','flow'))[ifo])+' Hz)')
                            bayeswavepsdnode[ifo].set_seed(randomseed)
                            if self.dataseed:
                                bayeswavepsdnode[ifo].set_dataseed(self.dataseed+event.event_id)
                if self.config.has_option('condor','bayesline') or self.config.getboolean('analysis','roq'):
                    if gotdata and event.event_id not in self.prenodes.keys():
                        if prenode not in self.get_nodes():
                            self.add_node(prenode)
                            for ifo in ifos:
                                if self.config.getboolean('analysis','roq'):
                                    computeroqweightsnode[ifo]=self.add_rom_weights_node(ifo,prenode)
                                #self.add_node(computeroqweightsnode[ifo])
                                if self.config.has_option('input','injection-file'):
                                    freqDataFile=os.path.join(roqeventpath,'data-dump'+ifo+'-freqDataWithInjection.dat')
                                else:
                                    freqDataFile=os.path.join(roqeventpath,'data-dump'+ifo+'-freqData.dat')
                                prenode.add_output_file(freqDataFile)
                                prenode.add_output_file(os.path.join(roqeventpath,'data-dump'+ifo+'-PSD.dat'))
                                if self.config.has_option('condor','bayesline'):
                                    bayeslinenode[ifo]=self.add_bayesline_node(ifo,prenode)
                                    bayeslinenode[ifo].add_var_arg('-i '+freqDataFile)
                                    bayeslinenode[ifo].add_input_file(freqDataFile)
                                    bayeslinenode[ifo].add_var_arg('-o '+os.path.join(roqeventpath,'BayesLine_PSD_'+ifo+'.dat'))
                                    bayeslinenode[ifo].add_output_file(os.path.join(roqeventpath,'BayesLine_PSD_'+ifo+'.dat'))
                                if self.config.getboolean('analysis','roq'):
                                    computeroqweightsnode[ifo].add_var_arg('--fHigh '+str(prenode.fhighs[ifo]))
                                    computeroqweightsnode[ifo].add_var_arg('--data '+freqDataFile)
                                    computeroqweightsnode[ifo].add_input_file(freqDataFile)
                                    computeroqweightsnode[ifo].add_var_arg('--psd '+os.path.join(roqeventpath,'data-dump'+ifo+'-PSD.dat'))
                                    computeroqweightsnode[ifo].add_input_file(os.path.join(roqeventpath,'data-dump'+ifo+'-PSD.dat'))
                                    computeroqweightsnode[ifo].add_var_arg('--out '+roqeventpath)
                                    computeroqweightsnode[ifo].add_output_file(os.path.join(roqeventpath,'weights_quadratic_'+ifo+'.dat'))
                            #self.prenodes[seg.id()]=(prenode,computeroqweightsnode)
                            if self.config.has_option('condor','bayesline'):
                                self.prenodes[event.event_id]=(prenode,bayeslinenode)
                            if self.config.getboolean('analysis','roq'):
                                self.prenodes[event.event_id]=(prenode,computeroqweightsnode)

            if self.config.has_option('condor','bayesline') or self.config.getboolean('analysis','roq'):
                #node.add_parent(self.prenodes[seg.id()][1][ifokey])
                node.add_parent(self.prenodes[event.event_id][1][ifokey])

            if dfnode is not None and dfnode not in self.get_nodes():
                if not self.config.has_option('lalinference','fake-cache'):
                    self.add_node(dfnode)

        if gotdata:
            self.add_node(node)
        else:
            print('no data found for time %f'%(end_time))
            return None, bayeswavepsdnode
        if extra_options is not None:
            for opt in extra_options.keys():
                node.add_var_arg('--'+opt+' '+extra_options[opt])
        # Add control options
        if self.config.has_option('input','injection-file'):
            node.set_injection(self.config.get('input','injection-file'),event.event_id)
            prenode.set_injection(self.config.get('input','injection-file'),event.event_id)
            if self.config.has_option('condor','bayeswave') and bayeswavepsdnode:
                for ifo in ifos:
                    bayeswavepsdnode[ifo].set_injection(self.config.get('input','injection-file'),event.event_id)
        if self.config.has_option('input','burst-injection-file'):
            node.set_injection(self.config.get('input','burst-injection-file'),event.event_id)
            prenode.set_injection(self.config.get('input','burst-injection-file'),event.event_id)
        if self.config.has_option('lalinference','seglen'):
            node.set_seglen(self.config.getfloat('lalinference','seglen'))
        elif self.config.has_option('engine','seglen'):
            node.set_seglen(self.config.getfloat('engine','seglen'))
        else:
            node.set_seglen(event.duration)
        if self.config.has_option('condor','bayeswave') and bayeswavepsdnode:
            node.set_psdlength(bw_seglen)
        if self.config.has_option('input','psd-length'):
            node.set_psdlength(self.config.getint('input','psd-length'))
            prenode.set_psdlength(self.config.getint('input','psd-length'))
            if self.config.has_option('condor','bayeswave') and bayeswavepsdnode:
                for ifo in ifos:
                    bayeswavepsdnode[ifo].set_psdlength(self.config.getint('input','psd-length'))
        if self.config.has_option('input','psd-start-time'):
            node.set_psdstart(self.config.getfloat('input','psd-start-time'))
            prenode.set_psdstart(self.config.getfloat('input','psd-start-time'))
            if self.config.has_option('condor','bayeswave') and bayeswavepsdnode:
                for ifo in ifos:
                    bayeswavepsdnode[ifo].set_psdstart(self.config.getfloat('input','psd-start-time'))
        node.set_max_psdlength(self.config.getint('input','max-psd-length'))
        prenode.set_max_psdlength(self.config.getint('input','max-psd-length'))
        node.set_padding(self.config.getint('input','padding'))
        prenode.set_padding(self.config.getint('input','padding'))
        if self.config.has_option('condor','bayeswave') and bayeswavepsdnode:
            for ifo in ifos:
                bayeswavepsdnode[ifo].set_max_psdlength(self.config.getint('input','max-psd-length'))
                bayeswavepsdnode[ifo].set_padding(self.config.getint('input','padding'))
        out_dir=os.path.join(self.basepath,'engine')
        mkdirs(out_dir)
        node.set_output_file(os.path.join(out_dir,node.engine+'-'+str(event.event_id)+'-'+node.get_ifos()+'-'+str(node.get_trig_time())+'-'+str(node.id)))
        if self.config.getboolean('analysis','roq'):
            for ifo in ifos:
                node.add_var_arg('--'+ifo+'-roqweightsLinear '+os.path.join(roqeventpath,'weights_linear_'+ifo+'.dat'))
                node.add_input_file(os.path.join(roqeventpath,'weights_linear_'+ifo+'.dat'))
                node.add_var_arg('--'+ifo+'-roqweightsQuadratic '+os.path.join(roqeventpath,'weights_quadratic_'+ifo+'.dat'))
                node.add_input_file(os.path.join(roqeventpath,'weights_quadratic_'+ifo+'.dat'))
            node.add_var_arg('--roqtime_steps '+os.path.join(roqeventpath,'roq_sizes.dat'))
            node.add_input_file(os.path.join(roqeventpath,'roq_sizes.dat'))
            node.add_var_arg('--roq-times '+os.path.join(roqeventpath,'tcs.dat'))
            node.add_input_file(os.path.join(roqeventpath,'tcs.dat'))
            node.add_var_arg('--roqnodesLinear '+os.path.join(roqeventpath,'fnodes_linear.dat'))
            node.add_input_file(os.path.join(roqeventpath,'fnodes_linear.dat'))
            node.add_var_arg('--roqnodesQuadratic '+os.path.join(roqeventpath,'fnodes_quadratic.dat'))
            node.add_input_file(os.path.join(roqeventpath,'fnodes_quadratic.dat'))
        if self.config.has_option('condor','bayesline'):
            for ifo in ifos:
                node.psds[ifo]=os.path.join(roqeventpath,'BayesLine_PSD_'+ifo+'.dat')
                node.add_input_file(os.path.join(roqeventpath,'BayesLine_PSD_'+ifo+'.dat'))
                prenode.psds[ifo]=os.path.join(roqeventpath,'BayesLine_PSD_'+ifo+'.dat')
                prenode.add_input_file(os.path.join(roqeventpath,'BayesLine_PSD_'+ifo+'.dat'))
        if self.config.has_option('condor','bayeswave') and bayeswavepsdnode:
            for ifo in ifos:
                node.psds[ifo]=os.path.join(roqeventpath,'BayesWave_PSD_'+ifo+'_IFO0_psd.dat')
                node.add_input_file(os.path.join(roqeventpath,'BayesWave_PSD_'+ifo+'_IFO0_psd.dat'))
                prenode.psds[ifo]=os.path.join(roqeventpath,'BayesWave_PSD_'+ifo+'_IFO0_psd.dat')
                prenode.add_input_file(os.path.join(roqeventpath,'BayesWave_PSD_'+ifo+'_IFO0_psd.dat'))
        for (opt,arg) in event.engine_opts.items():
            node.add_var_opt(opt,arg)
        if self.config.has_option('condor','bayeswave') and self.engine is not 'bayeswave':
            for ifo in ifos:
                node.add_parent(bayeswavepsdnode[ifo])
                prenode.add_parent(bayeswavepsdnode[ifo])
        return node,bayeswavepsdnode

    def add_results_page_node(self,resjob=None,outdir=None,parent=None,extra_options=None,gzip_output=None,ifos=None):
        if resjob is None:
            resjob=self.results_page_job
        node=ResultsPageNode(resjob)
        if parent is not None:
            node.add_parent(parent)
            infile=parent.get_pos_file()
            node.add_file_arg(infile)
        node.set_output_path(outdir)
        if gzip_output is not None:
            node.set_gzip_output(gzip_output)
        if ifos is not None:
            if isinstance(ifos,list):
                pass
            else:
                ifos=[ifos]
            node.set_ifos(ifos)

        self.add_node(node)
        return node

    def add_gracedb_start_node(self,gid,name='',parent=None):

        node=GraceDBNode(self.gracedbjob,parent=parent,gid=gid,command='log',tag='pe')
        node.set_message(name+' online parameter estimation started.')
        self.add_node(node)
        return node

    def add_gracedb_log_node(self,respagenode,gid):
        nodes=[]
        node=GraceDBNode(self.gracedbjob,parent=respagenode,gid=gid,command='log',tag='pe')
        resurl=respagenode.webpath.replace(self.gracedbjob.basepath,self.gracedbjob.baseurl)
        #node.set_message('online parameter estimation results:  '+resurl+'/posplots.html')
        node.set_message("LALInference online parameter estimation finished. <a href="+resurl+"/posplots.html>results</a>")
        self.add_node(node)
        nodes.append(node)

        tag='pe'
        if self.config.has_option('analysis','add-lvem-tag'):
            if self.config.getboolean('analysis','add-lvem-tag'):
                tag='pe,lvem'

        node=GraceDBNode(self.gracedbjob,parent=respagenode,gid=gid,command='upload',tag=tag)
        node.set_filename(respagenode.webpath+'/corner/extrinsic.png')
        self.add_node(node)
        nodes.append(node)

        node=GraceDBNode(self.gracedbjob,parent=respagenode,gid=gid,command='upload',tag='pe')
        node.set_filename(respagenode.webpath+'/corner/intrinsic.png')
        self.add_node(node)
        nodes.append(node)

        node=GraceDBNode(self.gracedbjob,parent=respagenode,gid=gid,command='upload',tag='pe')
        node.set_filename(respagenode.webpath+'/corner/sourceFrame.png')
        self.add_node(node)
        nodes.append(node)

        return nodes

    def add_gracedb_FITSskymap_upload(self,event,engine=None):
        gid=event.GID
        if gid is None:
            return
        if engine=='lalinferenceburst':
            prefix='LIB'
        elif engine is None:
            prefix="skymap"
        else:
            prefix='LALInference'
        nodes=None
        if self.config.has_option('condor','skyarea'):
            if self.config.has_option('analysis','upload-to-gracedb'):
                if self.config.getboolean('analysis','upload-to-gracedb'):
                    tag='sky_loc'
                    if self.config.has_option('analysis','add-lvem-tag'):
                        if self.config.getboolean('analysis','add-lvem-tag'):
                            tag='sky_loc,lvem'
                    skynodes=filter(lambda x: isinstance(x,SkyMapNode) ,self.get_nodes())
                    nodes=[]
                    for sk in skynodes:
                        if len(sk.ifos)>1:
                            node=GraceDBNode(self.gracedbjob,parent=sk,gid=gid,tag=tag)
                            #for p in sk.__parents:
                            #  if isinstance(p,ResultPageNode):
                            #    resultpagenode=p
                            node.set_filename(sk.outfits)
                            node.set_message('%s FITS sky map'%prefix)
                            self.add_node(node)
                            nodes.append(node)
        return nodes

    def add_rom_weights_node(self,ifo,parent=None):
        #try:
        #node=self.computeroqweightsnodes[ifo]
                #except KeyError:
        node=ROMNode(self.computeroqweights_job,ifo,parent.seglen,parent.flows[ifo])
        self.computeroqweightsnode[ifo]=node
        if parent is not None:
            node.add_parent(parent)
        self.add_node(node)
        return node

    def add_bayesline_node(self,ifo,parent=None):
        node=BayesLineNode(self.bayesline_job)
        self.bayeslinenode[ifo]=node
        if parent is not None:
            node.add_parent(parent)
        self.add_node(node)
        return node

    def add_bayeswavepsd_node(self,ifo,parent=None):
        node=BayesWavePSDNode(self.bayeswavepsd_job[ifo])
        self.bayeswavepsdnode[ifo]=node
        if parent is not None:
            node.add_parent(parent)
        self.add_node(node)
        return node

class SingularityJob(pipeline.CondorDAGJob):
    """
    Wrapper class for running jobs via a singularity image
    """
    # Hard-coded frame location in cvmfs
    CVMFS_FRAMES="/cvmfs/oasis.opensciencegrid.org/ligo/frames/"
    image=None
    def __init__(self, cp, *args, **kwargs):
        self.requirements=[]
        # Dir in which the DAG will run
        # Execute from the basedir so all paths can be resolved
        if cp.has_option('analysis','singularity'):
            self.singularity = cp.getboolean('analysis','singularity')
        else:
            self.singularity = False
        if cp.has_option('analysis','osg'):
            self.osg=cp.getboolean('analysis','osg')
        else:
            self.osg=False
        if self.osg and not self.singularity:
            raise Exception("Running on the OSG requires singularity=True in [analysis] section of ini file")
        self.basedir = cp.get('paths','basedir')
        self.add_condor_cmd('initialdir',self.basedir)
        if not self.singularity:
            return
        if cp.has_option('condor','singularity'):
            self.singularity_path = cp.get('condor','singularity')
        else:
            self.singularity_path = 'singularity' # must be in $PATH
        if cp.has_option('singularity','image'):
            self.image = cp.get('singularity','image')
        else:
            print("ERROR: You requested a singularity run but did not specify \
                an image file in the [singularity] section of the config file")
            sys.exit(-1)
        # If running on the OSG with real data, use frames from CVMFS
        if cp.has_option('lalinference','fake-cache') or not self.osg:
            extra_paths=""
        else:
            extra_paths="--bind {cvmfs_frames}".format(cvmfs_frames = self.CVMFS_FRAMES)
            self.add_condor_cmd('+SingularityBindCVMFS','True')
            self.add_condor_cmd('use_x509userproxy','true')
        if cp.has_option('analysis','roq') and cp.getboolean('analysis','roq'):
            extra_paths+=" --bind {roqpath}".format(roqpath=cp.get('paths','roq_b_matrix_directory'))

        self.wrapper_string="""
            echo "Workspace on execute node $(hostname -f)"
            echo "PWD" ${{PWD}}
            echo "contents"
            ls -l
            set -e
            {executable} \\
                    "$@"
            """.format(singularity = self.singularity_path,
                    basedir=self.basedir,
                    extra_paths = extra_paths,
                    executable = super(SingularityJob,self).get_executable(),
                    image = self.image
                    )

        if self.osg:
            self.add_condor_cmd('+OpenScienceGrid','True')
            self.requirements.append('IS_GLIDEIN==True')
            # Add requested sites if specified
            if cp.has_option('condor','desired-sites'):
                self.add_condor_cmd('+DESIRED_Sites',cp.get('condor','desired-sites'))
        if self.singularity:
            self.requirements.append('HasSingularity == TRUE')
            self.add_condor_cmd('+SingularityImage','"{0}"'.format(self.image))
        # Add data transfer options
        self.add_condor_cmd('should_transfer_files','YES')
        self.add_condor_cmd('when_to_transfer_output','ON_EXIT_OR_EVICT')

    def write_script(self,path):
        """
        Write the wrapper script
        """
        f=open(path,'w')
        f.writelines('#!/usr/bin/env bash')
        f.writelines(self.wrapper_string)
        f.close()
        os.chmod(path,0o755)

    def write_sub_file(self):
        """
        Over-load CondorDAGJob.write_sub_file to write the wrapper script and
        set the exe to call it
        """
        true_exec = self.get_executable()
        self.add_condor_cmd('requirements','&&'.join('({0})'.format(r) for r in self.requirements))
        if self.singularity:
            # Write the wrapper script
            wrapper=os.path.splitext(self.get_sub_file())[0] +'_wrapper.sh'
            self.write_script( wrapper )
            # Over-write the executable to set the wrapper script
            self.set_executable( wrapper )
        # Call the parent method to do the rest
        super(SingularityJob,self).write_sub_file()
        # Put the true exe back just in case
        self.set_executable(true_exec)



class SingularityNode(pipeline.CondorDAGNode):
    """
    Class representing a node run via singularity
    """
    def add_output_file(self,filename):
        filename=os.path.relpath(filename,start=self.job().get_config('paths','basedir'))
        self.add_output_macro(filename)
        super(SingularityNode,self).add_output_file(filename)
    def add_input_file(self,filename):
        filename=os.path.relpath(filename,start=self.job().get_config('paths','basedir'))
        self.add_input_macro(filename)
        super(SingularityNode,self).add_input_file(filename)
    def add_checkpoint_file(self,filename):
        filename=os.path.relpath(filename,start=self.job().get_config('paths','basedir'))
        self.add_checkpoint_macro(filename)
        super(SingularityNode,self).add_checkpoint_file(filename)
    def add_file_opt(self, opt, filename, file_is_output_file=False):
        # The code option needs the path as seen inside singularity, i.e. relative to basedir
        relfile=os.path.relpath(filename,start=self.job().get_config('paths','basedir'))
        self.add_var_opt(opt,relfile)
        if file_is_output_file:
            self.add_output_file(filename)
        else:
            self.add_input_file(filename)

class EngineJob(SingularityJob,pipeline.AnalysisJob):
    def __init__(self,cp,submitFile,logdir,engine,ispreengine=False,dax=False,site=None, *args, **kwargs):
        self.ispreengine=ispreengine
        self.engine=engine
        basepath=cp.get('paths','basedir')
        if ispreengine is True:
            roqpath=os.path.join(basepath,'ROQdata')
            self.roqpath=roqpath
            mkdirs(roqpath)
            exe=cp.get('condor',self.engine)
            universe='vanilla'
        else:
            if self.engine=='lalinferencemcmc':
                exe=cp.get('condor','mpiwrapper')
                universe="vanilla"
            elif self.engine=='lalinferencenest' or self.engine=='lalinferenceburst':
                exe=cp.get('condor',self.engine)
                if site is not None and site!='local':
                    universe='vanilla'
                else:
                    # Run in the vanilla universe when using resume
                    if cp.has_option('engine','resume'):
                        universe='vanilla'
                    else:
                        # This will go away, since condor_compile will not be supported in future
                        universe='standard'
            else:
                print('LALInferencePipe: Unknown engine node type %s!'%(self.engine))
                sys.exit(1)

        pipeline.CondorDAGJob.__init__(self,universe,exe)
        pipeline.AnalysisJob.__init__(self,cp,dax=dax)
        SingularityJob.__init__(self, cp)
        # locations for file IO
        self.add_condor_cmd('transfer_input_files','caches,engine,$(macroinput)')
        self.add_condor_cmd('transfer_output_files','engine')

        if cp.has_option('condor','accounting_group'):
            self.add_condor_cmd('accounting_group',cp.get('condor','accounting_group'))
        if cp.has_option('condor','accounting_group_user'):
            self.add_condor_cmd('accounting_group_user',cp.get('condor','accounting_group_user'))
        try:
            hostname=socket.gethostbyaddr(socket.gethostname())[0]
        except:
            hostname='Unknown'
        requirements=''
        if cp.has_option('condor','queue'):
            self.add_condor_cmd('+'+cp.get('condor','queue'),'True')
            requirements='(TARGET.'+cp.get('condor','queue')+' =?= True)'
        if cp.has_option('condor','Requirements'):
            if requirements!='':
                requirements=requirements+' && '
            requirements=requirements+cp.get('condor','Requirements')
        if requirements!='':
            self.add_condor_cmd('Requirements',requirements)
        # Set grid site if needed
        if cp.has_option('engine','resume'):
            self.resume=True
        else:
            self.resume=False
        if site:
            self.set_grid_site(site)
            if site!='local':
                self.set_executable_installed(False)
        # Set the options which are always used
        self.set_sub_file(os.path.abspath(submitFile))
        self.add_condor_cmd('getenv','true')
        if self.engine=='lalinferencemcmc':
            self.binary=cp.get('condor',self.engine.replace('mpi',''))
            self.mpirun=cp.get('condor','mpirun')
            if cp.has_section('mpi'):
                if ispreengine is False:
                    self.machine_count=cp.get('mpi','machine-count')
                    self.machine_memory=cp.get('mpi','machine-memory')
                    if cp.has_option('mpi','mpi_task_count'):
                        self.mpi_task_count=cp.get('mpi','mpi_task_count')
                    else:
                        self.mpi_task_count=self.machine_count
                else:
                    self.machine_count=str(1)
                    self.machine_memory=cp.get('mpi','machine-memory')
                    self.mpi_task_count=self.machine_count
            else:
                self.machine_count=str(1)
                self.machine_memory=str(1024) # default value if the user did not specify something
                self.mpi_task_count=self.machine_count
            #self.add_condor_cmd('machine_count',machine_count)
            #self.add_condor_cmd('environment','CONDOR_MPI_PATH=%s'%(openmpipath))
            if hostname=='pcdev1.phys.uwm.edu':
                self.add_condor_cmd('Requirements','CAN_RUN_MULTICORE')
                self.add_condor_cmd('+RequiresMultipleCores','True')
            self.add_condor_cmd('request_cpus',self.machine_count)
            self.add_condor_cmd('request_memory',str(float(self.machine_count)*float(self.machine_memory)))
        if self.engine=='lalinferencenest':
            self.add_condor_cmd('request_memory','4000') # 4GB RAM for high SNR BNS
        if cp.has_section(self.engine):
            if not ispreengine:
                self.add_ini_opts(cp,self.engine)
        if  cp.has_section('engine'):
            self.add_ini_opts(cp,'engine')
        self.set_stdout_file(os.path.join(logdir,'lalinference-$(cluster)-$(process)-$(node).out'))
        self.set_stderr_file(os.path.join(logdir,'lalinference-$(cluster)-$(process)-$(node).err'))
        # For LALInferenceNest demand only 1 thread (to be tuned later)
        if self.engine=='lalinferencenest':
            self.add_condor_cmd('environment','OMP_NUM_THREADS=1')
        if cp.has_option('condor','notification'):
            self.set_notification(cp.get('condor','notification'))
            if cp.has_option('resultspage','email'):
                self.add_condor_cmd('notify_user',cp.get('resultspage','email'))

    def set_grid_site(self,site=None):
        """
        Over-load base class method to choose condor universe properly
        """
        if self.engine=='lalinferencenest' or self.engine=='lalinferenceburst':
            if site is not None and site!='local':
                self.set_universe('vanilla')
            else:
                if self.resume:
                    self.set_universe('vanilla')
                else:
                    self.set_universe('standard')
        else:
            self.set_universe('vanilla')
        pipeline.CondorDAGJob.set_grid_site(self,site)

class EngineNode(SingularityNode):
    new_id = itertools.count()
    def __init__(self,li_job):
        super(EngineNode,self).__init__(li_job)
        self.ifos=[]
        self.scisegs={}
        self.channels={}
        self.psds={}
        self.flows={}
        self.fhighs={}
        self.timeslides={}
        self.seglen=None
        self.psdlength=None
        self.padding=None
        self.maxlength=None
        self.psdstart=None
        self.snrfile=None
        self.psdfiles=None
        self.cachefiles={}
        if li_job.ispreengine is False:
            self.id=next(EngineNode.new_id)
        self.__finaldata=False
        self.fakedata=False
        self.lfns=[] # Local file names (for frame files and pegasus)

    def set_seglen(self,seglen):
        self.seglen=seglen

    def set_psdlength(self,psdlength):
        self.psdlength=psdlength

    def set_max_psdlength(self,psdlength):
        self.maxlength=psdlength

    def set_padding(self,padding):
        self.padding=padding

    def set_psdstart(self,psdstart):
        self.psdstart=psdstart

    def set_seed(self,seed):
        self.add_var_opt('randomseed',str(seed))

    def set_srate(self,srate):
        self.add_var_opt('srate',str(srate))

    def set_trigSNR(self,trigSNR):
        self.add_var_opt('trigger-snr',str(trigSNR))

    def set_horizon_distance(self,horizon_distance):
        self.add_var_opt('distance-max',str(horizon_distance))

    def set_dataseed(self,seed):
        self.add_var_opt('dataseed',str(seed))

    def get_ifos(self):
        return ''.join(map(str,self.ifos))

    def set_psd_files(self):
        #if node.psdspath is not None:
        try:
            #bambi store the outfile in fileroot
            pathroot=self.fileroot
        except:
            #mcmc and nest in postfile
            pathroot=self.posfile
            if pathroot[-3:]=='.00':
                pathroot=pathroot[:-3]
        st=""
        for i in self.ifos:
            tmpst="%s%s-PSD.dat,"%(pathroot,i)
            st+=tmpst
            self.add_output_file(tmpst[:-1])
        st=st[:-1]
        self.psdfiles=st

    def get_psd_files(self):
        return self.psdfiles

    def set_snr_file(self):
        try:
            #bambi store the outfile in fileroot
            pathroot=self.fileroot
        except:
            #mcmc and nest in postfile
            pathroot=self.posfile
            if pathroot[-3:]=='.00':
                pathroot=pathroot[:-3]
        st="%s_snr.txt"%pathroot
        self.add_output_file(st)
        self.snrfile=st

    def get_snr_file(self):
        return self.snrfile

    def set_trig_time(self,time):
        """
        Set the end time of the signal for the centre of the prior in time
        """
        self.__trigtime=float(time)
        self.add_var_opt('trigtime',str(time))

    def set_event_number(self,event):
        """
        Set the event number in the injection XML.
        """
        if event is not None:
            self.__event=int(event)
            self.add_var_opt('event',str(event))

    def set_injection(self,injfile,event):
        """
        Set a software injection to be performed.
        """
        self.add_file_opt('inj',injfile)
        self.set_event_number(event)

    def get_trig_time(self): return self.__trigtime

    def add_fake_ifo_data(self,ifo,sciseg,fake_cache_name,fake_channel_name,timeslide=0):
        """
        Dummy method to set up fake data without needing to run datafind
        """
        self.ifos.append(ifo)
        self.scisegs[ifo]=sciseg
        self.cachefiles[ifo]=fake_cache_name
        self.timeslides[ifo]=timeslide
        self.channels[ifo]=fake_channel_name
        self.fakedata=True
        return 1

    def add_ifo_data(self,ifo,sciseg,channelname,timeslide=0):
        if self.ifos != ifo:
            self.ifos.append(ifo)
        self.scisegs[ifo]=sciseg
        parent=sciseg.get_df_node()
        if parent is not None:
            self.add_parent(parent)
            df_output=parent.get_output()
            self.set_cache(df_output,ifo)
            #self.cachefiles[ifo]=parent.get_output_files()[0]
            #self.add_input_file(self.cachefiles[ifo])
            self.timeslides[ifo]=timeslide
            self.channels[ifo]=channelname
            return 1
        else: return 0

    def set_cache(self,filename,ifo):
        """
        Add a cache file from LIGODataFind. Based on same method from pipeline.AnalysisNode
        """
        #print 'Adding cache files %s'%(str(filename))
        if isinstance(filename,str): # A normal lal cache file
            self.cachefiles[ifo]=filename
            self.add_input_file(filename)
        elif isinstance(filename,list): # A list of LFNs (for DAX mode)
            self.add_var_opt('glob-frame-data',' ')
            if len(filename) == 0:
                raise pipeline.CondorDAGNodeError(
                    "LDR did not return any LFNs for query: check ifo and frame type")
            for lfn in filename:
                self.lfns.append(lfn)

    def finalize(self):
        if not self.__finaldata:
            self._finalize_ifo_data()
        pipeline.CondorDAGNode.finalize(self)

    def _finalize_ifo_data(self):
        """
        Add final list of IFOs and data to analyse to command line arguments.
        """
        for ifo in self.ifos:
            self.add_var_arg('--ifo '+ifo)
            if self.fakedata:
                self.add_var_opt('%s-cache'%(ifo),self.cachefiles[ifo])
            elif not self.lfns:
                self.add_file_opt('%s-cache'%(ifo),self.cachefiles[ifo])
            self.add_var_opt('%s-channel'%(ifo),self.channels[ifo])
            if self.flows: self.add_var_opt('%s-flow'%(ifo),self.flows[ifo])
            if self.fhighs: self.add_var_opt('%s-fhigh'%(ifo),self.fhighs[ifo])
            if self.psds:
                self.add_var_opt('%s-psd'%(ifo),self.psds[ifo])
                #self.add_input_file(self.psds[ifo])
            if any(self.timeslides): self.add_var_opt('%s-timeslide'%(ifo),self.timeslides[ifo])

        """ The logic here is the following:
                The CBC code starts from the earliest commont time, but that means that if you run on *the same trigtime* the PSD start and PSDlength you'll get will be different, depending on wheather you are running on only one event or several, and the exact position of the event you are interested in in the list of times.
                Instead for each event (including single IFO runs) we do:
                a) get its trigtime
                b) set PSDlengh=maxPSD (requested by the user or equal to 32seglen)
                c) go define GPSstart= trigtime - PSDlength - seglen - padding -2

                By definition this means that GPSstart+ PSDlengh with never overlap with trigtime. Furthermore running on the same event will lead to the same PSDstart and lenght, no matter of whether that is a one-event or multi-event run.
                We should check that the PSDstart so obtained is in science mode. This is what the while loop 9 lines below is meant for. However that part is not active yet because I need to learn how to use scisegs. That is not a problem right now since we do run with disable-science (Since the searches will already have checked that the ~1-2 minutes of time prior to the event are in science. It might be a problem if one runs with hour-long slides).
        """
        trig_time=self.get_trig_time()
        maxLength=self.maxlength
        offset=(maxLength+self.seglen+2+self.padding)
        self.GPSstart=trig_time-offset
        self.__GPSend=0
        length=maxLength
        dt=self.seglen/4.

        while(self.GPSstart+length>=trig_time):
        ### or self.GPSstart not in  Science) --><-- here we should also have checked that we are in science mode, but I'm not sure how to do that yet.
            self.GPSstart+=dt
            length-=dt

        if self.psdstart is not None:
            self.GPSstart=self.psdstart
            #print 'Over-riding start time to user-specified value %f'%(self.GPSstart)
            #if self.GPSstart<starttime or self.GPSstart>endtime:
            #  print 'ERROR: Over-ridden time lies outside of science segment!'
            #  raise Exception('Bad psdstart specified')
        self.add_var_opt('psdstart',str(self.GPSstart))
        if self.psdlength is None:
            self.psdlength=length
            if(self.psdlength>self.maxlength):
                self.psdlength=self.maxlength
        self.add_var_opt('psdlength',self.psdlength)
        self.add_var_opt('seglen',self.seglen)
        for lfn in self.lfns:
            a, b, c, d = lfn.split('.')[0].split('-')
            t_start = int(c)
            t_end = int(c) + int(d)
            data_end=max(self.GPSstart+self.psdlength,trig_time+2)
            if( t_start <= data_end and t_end>self.GPSstart):
            #if (t_start <= (self.GPSstart+self.psdlength or t_start <=trig_time+2 or t_end >=) \
            #    and ( (t_end <= (self.GPSstart+self.psdlength )) or (t_end <= trig_time+2) ))  :
                self.add_input_file(lfn)
        self.__finaldata=True

class LALInferenceNestNode(EngineNode):
    def __init__(self,li_job):
        super(LALInferenceNestNode,self).__init__(li_job)
        self.engine='lalinferencenest'
        self.outfilearg='outfile'

    def set_output_file(self,filename):
        self.nsfile=filename+'.hdf5'
        self.posfile=self.nsfile
        self.add_file_opt(self.outfilearg,self.nsfile,file_is_output_file=True)
        self.Bfilename=self.nsfile+'_B.txt'
        self.add_output_file(self.Bfilename)
        self.headerfile=self.nsfile+'_params.txt'
        self.add_output_file(self.headerfile)
        if self.job().resume:
            self.add_checkpoint_file(self.nsfile+'_resume')

    def get_ns_file(self):
        return self.nsfile

class LALInferenceBurstNode(LALInferenceNestNode):
    def __init__(self,li_job):
        super(LALInferenceNestNode,self).__init__(li_job)
        self.engine='lalinferenceburst'
        self.outfilearg='outfile'
    def set_injection(self,injfile,event):
        """
        Set a software injection to be performed.
        """
        self.add_file_opt('binj',injfile)
        self.set_event_number(event)

class LALInferenceMCMCNode(EngineNode):
    def __init__(self,li_job):
        super(LALInferenceMCMCNode,self).__init__(li_job)
        self.engine='lalinferencemcmc'
        self.outfilearg='outfile'
        self.add_var_opt('mpirun',li_job.mpirun)
        self.add_var_opt('np',str(li_job.mpi_task_count))
        self.add_var_opt('executable',li_job.binary)

    def set_output_file(self,filename):
        self.posfile=filename+'.hdf5'
        # Should also take care of the higher temperature outpufiles with
        # self.add_output_file, getting the number of files from machine_count
        self.add_file_opt(self.outfilearg,self.posfile,file_is_output_file=True)

    def get_pos_file(self):
        return self.posfile

class LALInferenceDataDumpNode(EngineNode):
    def __init__(self,li_job):
        super(LALInferenceDataDumpNode,self).__init__(li_job)
        self.engine='lalinferencedatadump'
        self.outfilearg='outfile'
    def set_output_file(self,filename):
        pass

class BayesWavePSDJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):
    """
    Class for a BayesWave job

    Make sure all necessary commands are given for O2 BayesWave
    """
    def __init__(self,cp,submitFile,logdir,dax=False):
        exe=cp.get('condor','bayeswave')
        pipeline.CondorDAGJob.__init__(self,"vanilla",exe)
        pipeline.AnalysisJob.__init__(self,cp,dax=dax)
        if cp.has_section('bayeswave'):
            self.add_ini_opts(cp,'bayeswave')
        if cp.has_option('condor','accounting_group'):
            self.add_condor_cmd('accounting_group',cp.get('condor','accounting_group'))
        if cp.has_option('condor','accounting_group_user'):
            self.add_condor_cmd('accounting_group_user',cp.get('condor','accounting_group_user'))
        requirements=''
        if cp.has_option('condor','queue'):
            self.add_condor_cmd('+'+cp.get('condor','queue'),'True')
            requirements='(TARGET.'+cp.get('condor','queue')+' =?= True)'
        if cp.has_option('condor','Requirements'):
            if requirements!='':
                requirements=requirements+' && '
            requirements=requirements+cp.get('condor','Requirements')
        if requirements!='':
            self.add_condor_cmd('Requirements',requirements)
        self.set_sub_file(submitFile)
        self.set_stdout_file(os.path.join(logdir,'bayeswavepsd-$(cluster)-$(process).out'))
        self.set_stderr_file(os.path.join(logdir,'bayeswavepsd-$(cluster)-$(process).err'))
        self.add_condor_cmd('getenv','True')
        self.ispreengine = False

class BayesWavePSDNode(EngineNode):
    def __init__(self,bayeswavepsd_job):
        super(BayesWavePSDNode,self).__init__(bayeswavepsd_job)
        self.engine='bayeswave'
        self.outfilearg='outfile'

    def set_output_file(self,filename):
        pass

class ResultsPageJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):
    def __init__(self,cp,submitFile,logdir,dax=False):
        exe=cp.get('condor','resultspage')
        pipeline.CondorDAGJob.__init__(self,"vanilla",exe)
        pipeline.AnalysisJob.__init__(self,cp,dax=dax) # Job always runs locally
        if cp.has_option('condor','accounting_group'):
            self.add_condor_cmd('accounting_group',cp.get('condor','accounting_group'))
        if cp.has_option('condor','accounting_group_user'):
            self.add_condor_cmd('accounting_group_user',cp.get('condor','accounting_group_user'))
        requirements=''
        if cp.has_option('condor','queue'):
            self.add_condor_cmd('+'+cp.get('condor','queue'),'True')
            requirements='(TARGET.'+cp.get('condor','queue')+' =?= True)'
        if cp.has_option('condor','Requirements'):
            if requirements!='':
                requirements=requirements+' && '
            requirements=requirements+cp.get('condor','Requirements')
        if requirements!='':
            self.add_condor_cmd('Requirements',requirements)
        self.set_sub_file(os.path.abspath(submitFile))
        self.set_stdout_file(os.path.join(logdir,'resultspage-$(cluster)-$(process).out'))
        self.set_stderr_file(os.path.join(logdir,'resultspage-$(cluster)-$(process).err'))
        self.add_condor_cmd('getenv','True')
        self.add_condor_cmd('request_memory','2000')
        self.add_ini_opts(cp,'resultspage')

        if cp.has_option('results','skyres'):
            self.add_opt('skyres',cp.get('results','skyres'))

class ResultsPageNode(pipeline.CondorDAGNode):
    def __init__(self,results_page_job,outpath=None):
        super(ResultsPageNode,self).__init__(results_page_job)
        if outpath is not None:
            self.set_output_path(path)
        self.__event=0
        self.ifos=None
        self.injfile=None

    def set_gzip_output(self,path):
        self.add_file_opt('archive',path,file_is_output_file=True)

    def set_output_path(self,path):
        self.webpath=path
        #self.add_file_opt('outpath',path,file_is_output_file=True)
        self.add_var_opt('outpath',path)
        #self.add_file_opt('archive','results.tar.gz',file_is_output_file=True)
        mkdirs(path)
        self.posfile=os.path.join(path,'posterior_samples.dat')
        self.add_output_file(self.posfile)

    def get_output_path(self):
        return self.webpath

    def set_injection(self,injfile,eventnumber):
        self.injfile=injfile
        self.add_file_opt('inj',injfile)
        self.set_event_number(eventnumber)

    def get_injection(self):
        return self.injfile

    def set_event_number(self,event):
        """
        Set the event number in the injection XML.
        """
        if event is not None:
            self.__event=int(event)
            self.add_var_arg('--eventnum '+str(event))

    def get_event_number(self):
        return self.__event

    def set_psd_files(self,st):
        if st is None:
            return
        for i in st.split(','):
            self.add_input_file(i)
        self.add_var_arg('--psdfiles %s'%st)

    def set_snr_file(self,st):
        if st is None:
            return
        self.add_file_opt('snr',st)

    def set_coinc_file(self,coinc,gid=None):
        if gid:
            self.__event=gid
        if coinc is None:
            return
        self.add_var_arg('--trig '+coinc)

    def add_engine_parent(self,node):
        """
        Add a parent node which is one of the engine nodes
        And automatically set options accordingly
        """
        self.add_parent(node)
        self.add_file_arg(node.get_pos_file())
        self.infiles.append(node.get_pos_file())

    def get_pos_file(self): return self.posfile

    def set_bayes_coherent_incoherent(self,bcifile):
        self.add_file_opt('bci',bcifile)

    def set_bayes_coherent_noise(self,bsnfile):
        self.add_file_opt('bsn',bsnfile)

    def set_header_file(self,headerfile):
        self.add_file_opt('header',headerfile)

    def set_ifos(self,ifos):
        self.ifos=ifos

class CoherenceTestJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):
    """
    Class defining the coherence test job to be run as part of a pipeline.
    """
    def __init__(self,cp,submitFile,logdir,dax=False):
        exe=cp.get('condor','coherencetest')
        pipeline.CondorDAGJob.__init__(self,"vanilla",exe)
        pipeline.AnalysisJob.__init__(self,cp,dax=dax)
        if cp.has_option('condor','accounting_group'):
            self.add_condor_cmd('accounting_group',cp.get('condor','accounting_group'))
        if cp.has_option('condor','accounting_group_user'):
            self.add_condor_cmd('accounting_group_user',cp.get('condor','accounting_group_user'))
        requirements=''
        if cp.has_option('condor','queue'):
            self.add_condor_cmd('+'+cp.get('condor','queue'),'True')
            requirements='(TARGET.'+cp.get('condor','queue')+' =?= True)'
        if cp.has_option('condor','Requirements'):
            if requirements!='':
                requirements=requirements+' && '
            requirements=requirements+cp.get('condor','Requirements')
        if requirements!='':
            self.add_condor_cmd('Requirements',requirements)
        self.add_opt('new-coherent-incoherent-noise','')
        self.add_condor_cmd('getenv','True')
        self.set_stdout_file(os.path.join(logdir,'coherencetest-$(cluster)-$(process).out'))
        self.set_stderr_file(os.path.join(logdir,'coherencetest-$(cluster)-$(process).err'))
        self.set_sub_file(os.path.abspath(submitFile))

class CoherenceTestNode(pipeline.CondorDAGNode):
    """
    Class defining the node for the coherence test
    """
    def __init__(self,coherencetest_job,outfile=None):
        super(CoherenceTestNode,self).__init__(coherencetest_job)
        self.incoherent_parents=[]
        self.coherent_parent=None
        self.finalized=False
        if outfile is not None:
            self.add_file_opt('outfile',outfile,file_is_output_file=True)

    def add_coherent_parent(self,node):
        """
        Add a parent node which is an engine node, and process its outputfiles
        """
        self.coherent_parent=node
        self.add_parent(node)
    def add_incoherent_parent(self,node):
        """
        Add a parent node which provides one of the single-ifo evidence values
        """
        self.incoherent_parents.append(node)
        self.add_parent(node)
    def finalize(self):
        """
        Construct command line
        """
        if self.finalized==True: return
        self.finalized=True
        self.add_file_arg(self.coherent_parent.get_pos_file())
        for inco in self.incoherent_parents:
            self.add_file_arg(inco.get_pos_file())

class MergeJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):
    """
    Class defining a job which merges several parallel nested sampling or MCMC jobs into a single file
    Input arguments:
    cp        - A configparser object containing the setup of the analysis
    submitFile    - Path to store the submit file
    logdir        - A directory to hold the stderr, stdout files of the merge runs
    dax      -  Is the job to be configured as a pegasus job, if so dax=True
    engine   - Set to either 'nest' or 'mcmc' for the appropriate behaviour
    """
    def __init__(self,cp,submitFile,logdir,dax=False,engine='nest'):
        if engine == 'mcmc':
            exe=cp.get('condor','mergeMCMCscript')
        else:
            exe=cp.get('condor','mergeNSscript')
        pipeline.CondorDAGJob.__init__(self,"vanilla",exe)
        pipeline.AnalysisJob.__init__(self,cp,dax=dax)
        if cp.has_option('condor','accounting_group'):
            self.add_condor_cmd('accounting_group',cp.get('condor','accounting_group'))
        if cp.has_option('condor','accounting_group_user'):
            self.add_condor_cmd('accounting_group_user',cp.get('condor','accounting_group_user'))
        requirements=''
        if cp.has_option('condor','queue'):
            self.add_condor_cmd('+'+cp.get('condor','queue'),'True')
            requirements='(TARGET.'+cp.get('condor','queue')+' =?= True)'
        if cp.has_option('condor','Requirements'):
            if requirements!='':
                requirements=requirements+' && '
            requirements=requirements+cp.get('condor','Requirements')
        if requirements!='':
            self.add_condor_cmd('Requirements',requirements)
        self.set_sub_file(os.path.abspath(submitFile))
        self.set_stdout_file(os.path.join(logdir,'merge-$(cluster)-$(process).out'))
        self.set_stderr_file(os.path.join(logdir,'merge-$(cluster)-$(process).err'))
        self.add_condor_cmd('getenv','True')
        if cp.has_option('merge','npos') and engine == 'nest':
            self.add_opt('npos',cp.get('merge','npos'))


class MergeNode(pipeline.CondorDAGNode):
    """
    Class defining the DAG node for a NS merge job
    Input arguments:
    merge_job = A MergeJob object
    parents = iterable of parent LALInferenceNest nodes (must have get_ns_file() method)
    engine   - Set to either 'nest' or 'mcmc' for the appropriate behaviour
    """
    def __init__(self,merge_job,parents=None,engine='nest'):
        super(MergeNode,self).__init__(merge_job)
        if parents is not None:
            for parent in parents:
                if engine == 'nest':
                    self.add_engine_parent(parent)
                else:
                    self.add_combine_parent(parent)

    def add_engine_parent(self,parent):
        self.add_parent(parent)
        self.add_file_arg(parent.get_ns_file())

    def add_combine_parent(self,parent):
        self.add_parent(parent)
        self.add_file_arg(parent.get_pos_file())

    def set_pos_output_file(self,file):
        self.add_file_opt('pos',file,file_is_output_file=True)
        self.posfile=file

    def get_pos_file(self): return self.posfile

class CombineMCMCJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):
    """
    Class defining a job which combines several parallel MCMC chains into a single hdf5 file
    Input arguments:
    cp        - A configparser object containing the setup of the analysis
    submitFile    - Path to store the submit file
    logdir        - A directory to hold the stderr, stdout files of the merge runs
    """
    def __init__(self,cp,submitFile,logdir,dax=False):
        exe=cp.get('condor','combinePTMCMCh5script')
        pipeline.CondorDAGJob.__init__(self,"vanilla",exe)
        pipeline.AnalysisJob.__init__(self,cp,dax=dax)
        if cp.has_option('condor','accounting_group'):
            self.add_condor_cmd('accounting_group',cp.get('condor','accounting_group'))
        if cp.has_option('condor','accounting_group_user'):
            self.add_condor_cmd('accounting_group_user',cp.get('condor','accounting_group_user'))
        self.set_sub_file(os.path.abspath(submitFile))
        self.set_stdout_file(os.path.join(logdir,'combine-$(cluster)-$(process).out'))
        self.set_stderr_file(os.path.join(logdir,'combine-$(cluster)-$(process).err'))
        self.add_condor_cmd('getenv','True')

class CombineMCMCNode(pipeline.CondorDAGNode):
    """
    Class defining the DAG node for a MCMC combine job
    Input arguments:
    combine_job = A CombineMCMCJob object
    parents = iterable of parent LALInferenceMCMC nodes (must have get_ns_file() method)
    """
    def __init__(self,combine_job,parents=None):
        super(CombineMCMCNode,self).__init__(combine_job)
        if parents is not None:
            for parent in parents:
                self.add_engine_parent(parent)

    def add_engine_parent(self,parent):
        self.add_parent(parent)
        self.add_file_arg(parent.get_pos_file())

    def get_parent_posfile(self,parent):
        return parent.get_pos_file()

    def set_pos_output_file(self,file):
        self.add_file_opt('outfile',file,file_is_output_file=True)
        self.posfile=file

    def get_pos_file(self): return self.posfile

class GraceDBJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):
    """
    Class for a gracedb job
    """
    def __init__(self,cp,submitFile,logdir,dax=False):
        exe=cp.get('condor','gracedb')
        pipeline.CondorDAGJob.__init__(self,"scheduler",exe)
        pipeline.AnalysisJob.__init__(self,cp,dax=dax)
        if cp.has_option('condor','accounting_group'):
            self.add_condor_cmd('accounting_group',cp.get('condor','accounting_group'))
        if cp.has_option('condor','accounting_group_user'):
            self.add_condor_cmd('accounting_group_user',cp.get('condor','accounting_group_user'))
        self.set_sub_file(os.path.abspath(submitFile))
        self.set_stdout_file(os.path.join(logdir,'gracedb-$(cluster)-$(process).out'))
        self.set_stderr_file(os.path.join(logdir,'gracedb-$(cluster)-$(process).err'))
        self.add_condor_cmd('getenv','True')
        self.basepath=cp.get('paths','webdir')
        self.baseurl=guess_url(self.basepath)

class GraceDBNode(pipeline.CondorDAGNode):
    """
    Run the gracedb executable to report the results
    """
    def __init__(self,gracedb_job,gid=None,parent=None,message=None,upfile=None,command='upload',tag=None):
            # Message need to be a string
            # Upfile is the full path of the file to be uploaded
        super(GraceDBNode,self).__init__(gracedb_job)
        if gid: self.set_gid(gid)
        if parent:
            if isinstance(parent, list):
                for p in parent:
                    self.add_parent(p)
            else:
                self.add_parent(parent)
        self.message=message
        self.filename=upfile
        self.command=command
        self.tag=tag
        self.__finalized=False

    def set_gid(self,gid):
        """
        Set the GraceDB ID to log to
        """
        self.gid=gid

    def set_message(self,message):
        self.message=message

    def set_filename(self,filename):
        self.filename=filename

    def finalize(self):
        if self.__finalized:
            return
        self.add_var_arg(self.command)
        if self.tag:
            self.add_var_arg('--tag-name='+self.tag)
        self.add_var_arg(str(self.gid))
        if self.filename:
            self.add_var_arg(self.filename+' ')
        if self.message:
            self.add_var_arg(self.message)
        self.__finalized=True

class ROMJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):
    """
    Class for a ROM compute weights job
    """
    def __init__(self,cp,submitFile,logdir,dax=False):
        time_step=0.000172895418228
        #This ensures that, for SNR < 100, a signal with bandwidth 20-4096 Hz will
        #have a resolved time posterior assuming a chirp-like frequency evolution
        #and aLIGO_ZDHP PSD
        dt=0.1
        exe=cp.get('condor','computeroqweights')
        pipeline.CondorDAGJob.__init__(self,"vanilla",exe)
        pipeline.AnalysisJob.__init__(self,cp,dax=dax)
        if cp.has_option('condor','accounting_group'):
            self.add_condor_cmd('accounting_group',cp.get('condor','accounting_group'))
        if cp.has_option('condor','accounting_group_user'):
            self.add_condor_cmd('accounting_group_user',cp.get('condor','accounting_group_user'))
        requirements=''
        if cp.has_option('condor','queue'):
            self.add_condor_cmd('+'+cp.get('condor','queue'),'True')
            requirements='(TARGET.'+cp.get('condor','queue')+' =?= True)'
        if cp.has_option('condor','Requirements'):
            if requirements!='':
                requirements=requirements+' && '
            requirements=requirements+cp.get('condor','Requirements')
        if requirements!='':
            self.add_condor_cmd('Requirements',requirements)
        self.set_sub_file(submitFile)
        self.set_stdout_file(os.path.join(logdir,'computeroqweights-$(cluster)-$(process).out'))
        self.set_stderr_file(os.path.join(logdir,'computeroqweights-$(cluster)-$(process).err'))
        self.add_condor_cmd('getenv','True')
        self.add_arg('-B '+str(cp.get('paths','roq_b_matrix_directory')))
        if cp.has_option('engine','dt'):
            dt=cp.getfloat('engine','dt')
        self.add_arg('-t '+str(dt))
        if cp.has_option('engine','time_step'):
            time_step=cp.get('engine','time_step')
        self.add_arg('-T '+str(time_step))
        if cp.has_option('condor','computeroqweights_memory'):
            computeroqweights_memory=str(cp.get('condor','computeroqweights_memory'))
        else:
            params = np.genfromtxt(str(cp.get('paths','roq_b_matrix_directory')+'/params.dat'), names=True)
            computeroqweights_memory=str(int(
            os.path.getsize(str(cp.get('paths','roq_b_matrix_directory')+'/B_linear.npy'))/(1024*1024)
            + 3*((params['fhigh']-params['flow'])*params['seglen'])*(float(dt+0.05)*2/float(time_step))*2*8/(1024*1024)
            + os.path.getsize(str(cp.get('paths','roq_b_matrix_directory')+'/B_quadratic.npy'))/(1024*1024)
            ) + 4096) # add 4gb of memory due to how matrix-copying is handled in lalapps_compute_roq_weights.py/numpy
            print('Requesting '+computeroqweights_memory+' of memory for computeroqweights.')
        self.add_condor_cmd('request_memory',computeroqweights_memory)

class ROMNode(pipeline.CondorDAGNode):
    """
    Run the ROM compute weights script
    """
    def __init__(self,computeroqweights_job,ifo,seglen,flow):
        super(ROMNode,self).__init__(computeroqweights_job)
        self.__finalized=False
        self.add_var_arg('--seglen '+str(seglen))
        self.add_var_arg('--fLow '+str(flow))
        self.add_var_arg('--ifo '+ifo)

    def finalize(self):
        if self.__finalized:
            return
        self.__finalized=True

class BayesLineJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):
    """
    Class for a BayesLine job
    """
    def __init__(self,cp,submitFile,logdir,dax=False):
        exe=cp.get('condor','bayesline')
        pipeline.CondorDAGJob.__init__(self,"vanilla",exe)
        pipeline.AnalysisJob.__init__(self,cp,dax=dax)
        if cp.has_option('condor','accounting_group'):
            self.add_condor_cmd('accounting_group',cp.get('condor','accounting_group'))
        if cp.has_option('condor','accounting_group_user'):
            self.add_condor_cmd('accounting_group_user',cp.get('condor','accounting_group_user'))
        requirements=''
        if cp.has_option('condor','queue'):
            self.add_condor_cmd('+'+cp.get('condor','queue'),'True')
            requirements='(TARGET.'+cp.get('condor','queue')+' =?= True)'
        if cp.has_option('condor','Requirements'):
            if requirements!='':
                requirements=requirements+' && '
            requirements=requirements+cp.get('condor','Requirements')
        if requirements!='':
            self.add_condor_cmd('Requirements',requirements)
        self.set_sub_file(submitFile)
        self.set_stdout_file(os.path.join(logdir,'bayesline-$(cluster)-$(process).out'))
        self.set_stderr_file(os.path.join(logdir,'bayesline-$(cluster)-$(process).err'))
        self.add_condor_cmd('getenv','True')

class BayesLineNode(pipeline.CondorDAGNode):
    """
    Run the BayesLine code
    """
    def __init__(self,bayesline_job):
        super(BayesLineNode,self).__init__(bayesline_job)
        self.__finalized=False

    def finalize(self):
        if self.__finalized:
            return
        self.__finalized=True


class SkyMapNode(pipeline.CondorDAGNode):
    def __init__(self, skymap_job, posfile=None, parent=None, objid=None, prefix=None, outdir=None):
        self.prefix=prefix
        super(SkyMapNode, self).__init__(skymap_job)
        self.objid=None
        self.outdir=None
        self.finalized=False
        if parent:
            self.add_parent(parent)
        if posfile:
            self.set_posfile(posfile)
        if objid:
            self.set_objid(objid)
        if outdir:
            self.set_outdir(outdir)
    def set_outdir(self, outdir):
        self.outdir=outdir
        if self.prefix:
            name = self.prefix+'.fits.gz'
        else:
            name = 'skymap.fits.gz'
        self.outfits = os.path.join(outdir,name)
    def set_posfile(self, posfile):
        self.posfile=posfile
    def get_outdir(self):
        return self.outdir
    def set_objid(self,objid):
        """
        Object ID for the fits file
        """
        self.objid = objid
    def finalize(self):
        """
        Construct command line
        """
        if self.finalized==True: return
        self.finalized=True
        self.add_input_file(self.posfile)
        self.add_file_arg(self.posfile)
        self.add_var_opt('fitsoutname',self.outfits)
        self.add_file_opt('outdir',self.outdir, file_is_output_file=True)
        if self.objid:
            self.add_var_opt('objid',self.objid)
        super(SkyMapNode,self).finalize()

class SkyMapJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):
    """
    Node to run ligo-skymap-from-samples
    """
    def __init__(self, cp, submitFile, logdir):
        exe=cp.get('condor','ligo-skymap-from-samples')
        pipeline.CondorDAGJob.__init__(self,"vanilla",exe)
        pipeline.AnalysisJob.__init__(self,cp)
        requirements=[]
        if cp.has_option('condor','accounting_group'):
            self.add_condor_cmd('accounting_group',cp.get('condor','accounting_group'))
        if cp.has_option('condor','accounting_group_user'):
            self.add_condor_cmd('accounting_group_user',cp.get('condor','accounting_group_user'))
        if cp.has_option('condor','queue'):
            self.add_condor_cmd('+'+cp.get('condor','queue'),'True')
            requirements.append('(TARGET.'+cp.get('condor','queue')+' =?= True)')
        if cp.has_option('condor','Requirements'):
            requirements.append(cp.get('condor','Requirements'))
        if requirements:
            self.add_condor_cmd('Requirements','&&'.join(requirements))
        self.set_sub_file(submitFile)
        self.set_stdout_file(os.path.join(logdir,'samples2map-$(cluster)-$(process).out'))
        self.set_stderr_file(os.path.join(logdir,'samples2map-$(cluster)-$(process).err'))
        # The user environment PYTHONPATH may be set to python2.7 version of lalsuite, so this is disabled
        #self.add_condor_cmd('getenv','True')
        # Add user-specified options from ini file
        self.add_ini_opts(cp,'ligo-skymap-from-samples')


class PlotSkyMapJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
    """
    Job to run ligo-skymap-plot
    """
    def __init__(self, cp, submitFile, logdir):
        exe=cp.get('condor','ligo-skymap-plot')
        pipeline.CondorDAGJob.__init__(self, "vanilla", exe)
        pipeline.AnalysisJob.__init__(self, cp)
        requirements=[]
        if cp.has_option('condor','accounting_group'):
            self.add_condor_cmd('accounting_group',cp.get('condor','accounting_group'))
        if cp.has_option('condor','accounting_group_user'):
            self.add_condor_cmd('accounting_group_user',cp.get('condor','accounting_group_user'))
        if cp.has_option('condor','queue'):
            self.add_condor_cmd('+'+cp.get('condor','queue'),'True')
            requirements.append('(TARGET.'+cp.get('condor','queue')+' =?= True)')
        if cp.has_option('condor','Requirements'):
            requirements.append(cp.get('condor','Requirements'))
        if requirements:
            self.add_condor_cmd('Requirements','&&'.join(requirements))
        self.set_sub_file(submitFile)
        self.set_stdout_file(os.path.join(logdir,'plotskymap-$(cluster)-$(process).out'))
        self.set_stderr_file(os.path.join(logdir,'plotskymap-$(cluster)-$(process).err'))
        # The user environment PYTHONPATH may be set to python2.7 version of lalsuite, so this is disabled
        # self.add_condor_cmd('getenv','True')
        # Add user-specified options from ini file
        self.add_ini_opts(cp,'ligo-skymap-plot')

class PlotSkyMapNode(pipeline.CondorDAGNode):
    def __init__(self, plotskymap_job, parent=None, inputfits = None, output=None):
        super(PlotSkyMapNode, self).__init__(plotskymap_job)
        if parent:
            self.add_parent(parent)
        if inputfits:
            self.set_input_fits(inputfits)
        if output:
            self.set_output(output)
        self.finalized=False
    def set_output(self, outfile):
        self.outfile=outfile
    def set_input_fits(self, fitsfile):
        self.fitsfile=fitsfile
    def get_output(self):
        return self.output
    def finalize(self):
        """
        Construct command line
        """
        if self.finalized==True: return
        self.finalized=True
        self.add_input_file(self.fitsfile)
        self.add_file_arg(self.fitsfile)
        self.add_file_opt('output',self.outfile, file_is_output_file=True)
        super(PlotSkyMapNode,self).finalize()

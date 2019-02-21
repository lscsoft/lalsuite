"""
Classes needed for ihope.
"""

from __future__ import print_function

__author__ = 'Stephen Fairhurst <sfairhur@gravity.phys.uwm.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'

##############################################################################
# import standard modules and append the lalapps prefix to the python path
import os, sys, copy, shutil, glob
from six.moves.configparser import NoSectionError
import optparse
import tempfile
import urllib
import subprocess

##############################################################################
# import the modules we need to build the pipeline
from ligo import segments
from ligo.segments import utils as segmentsUtils
from glue import pipeline
from glue import lal
from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue.ligolw.utils import process as ligolw_process
from glue.ligolw.utils import segments as ligolw_segments
from lalapps import inspiral

##############################################################################
# Functions used in setting up the dag:
def make_external_call(command, show_stdout=False, show_command=False):
    """
    Run a program on the shell and print informative messages on failure.
    """
    if show_command: print(command)

    p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
      shell=isinstance(command, str))
    out, err = p.communicate()

    if p.returncode != 0:
        print("External call failed.", file=sys.stderr)
        print("  stdout: %s" % out, file=sys.stderr)
        print("  stderr: %s" % err, file=sys.stderr)
        raise subprocess.CalledProcessError(p.returncode, command)
    if show_stdout:
        print(out)

##############################################################################
def mkdir( newdir ):
    """
    Create a directory, don't complain if it exists

    newdir = name of directory to be created
    """
    if os.path.isdir(newdir): pass
    elif os.path.isfile(newdir):
        raise OSError("a file with the same name as the desired " \
                      "dir, '%s', already exists." % newdir)
    else: os.mkdir(newdir)

##############################################################################
def hipe_cache(ifos, usertag, gps_start_time, gps_end_time):
    """
    return the name of the hipe cache

    ifos    = list of ifos
    usertag = usertag
    gps_start_time = start time of analysis
    gps_end_time   = end time of analysis
    """
    hipeCache = "".join(ifos) + "-INSPIRAL_HIPE"
    if usertag:
        hipeCache += "_" + usertag
    hipeCache += "-" + str(gps_start_time) + "-" + \
        str(gps_end_time - gps_start_time)  + ".cache"

    return hipeCache

##############################################################################
def hipe_pfn_glob_cache(globpat):
    """
    create and return the name of a pfn cache containing files that match
    globpat. This is needed to manage the .input files that hipe creates.

    cachename = the name of the pfn cache file
    globpat = the pattern to search for
    """
    cache_list = []
    for file in glob.glob(globpat):
        lfn = os.path.basename(file)
        pfn = "file://" + os.path.join(os.getcwd(),file)
        cache_list.append( (lfn, pfn, "local") )
    return cache_list

##############################################################################
def hipe_pfn_list_cache(files):
    """
    create and return the name of a pfn cache containing files in files.
    This is needed to manage the .input files that hipe creates.

    cachename = the name of the pfn cache file
    files = a list of files
    """
    cache_list = []
    for file in files:
        lfn = os.path.basename(file)
        pfn = "file://" + os.path.abspath(file)
        cache_list.append( (lfn, pfn, "local") )
    return cache_list


##############################################################################
def tmpltbank_cache(datafind_filename):
    """
    open and read the datafind hipe cache, find only the tmpltbank files
    return a lal.Cache class with only the tmpltbanks

    datafind_filename = datafind cache filename relative to current path
    """
    cache_file = open(datafind_filename)
    result_cache = lal.Cache.fromfile(cache_file)
    cache_file.close()

    return result_cache.sieve(description="TMPLTBANK")

##############################################################################
def symlink_tmpltbank(tmpltbank_cache, user_tag):
    """
    Symlinks the tmpltbank files from the datafind directory into the current
    directory
    Returns

    tmpltbank_cache = lal.Cache of the tmpltbank files in the datafind dir
    usertag = playground, full_data, or the name of the injection directory
    """
    new_dir = user_tag.lower()
    file_tag = user_tag.upper()

    for entry in tmpltbank_cache:
        old_file = os.path.basename(entry.path)
        new_file = old_file.replace("DATAFIND",file_tag)
        try: # Remove file if it already exists
            os.remove(new_file)
        except: pass
        os.symlink("../datafind/" + old_file, new_file)

    new_pfnlist = [url.replace("datafind",new_dir).replace("DATAFIND",file_tag)
        for url in tmpltbank_cache.pfnlist()]

    return lal.Cache.from_urls(new_pfnlist)

##############################################################################
# Function to set up the segments for the analysis
def science_segments(ifo, config, generate_segments = True):
    """
    generate the segments for the specified ifo
    @param ifo:    the name of the ifo
    @param config: the configParser object with analysis details
    @param generate_segments: whether or not to actually make segments
    """
    start = config.get("input","gps-start-time")
    end = config.get("input","gps-end-time")

    segFindFile = ifo + "-SCIENCE_SEGMENTS-" + start + "-" + \
        str(int(end) - int(start))

    segFindXML = segFindFile + ".xml"
    segFindFile = segFindFile + ".txt"

    # if not generating segments, all we need is the name of the segment file
    if not generate_segments: return segFindFile

    whichtoanalyze = ifo.lower() + "-analyze"
    print("For " + ifo + ", analyze " +config.get("segments",whichtoanalyze))

    # run segFind to determine science segments
    segFindCall = ' '.join([ config.get("condor", "segfind"),
          "--query-segments",
          "--segment-url", config.get("segfind", "segment-url"),
          "--gps-start-time", start,
          "--gps-end-time", end,
          "--include-segments",
             ''.join([ '"', config.get("segments", whichtoanalyze), '"' ]),
          "--output-file", segFindXML ])
    make_external_call(segFindCall)

    segsToTxtCall = ' '.join([ config.get("condor", "ligolw_print"),
          "--table segment --column start_time --column end_time",
          "--delimiter ' '",
          segFindXML,
          ">", segFindFile ])
    make_external_call(segsToTxtCall)

    return segFindFile

##############################################################################
# the hardware injection script part one to get segments for it to use later
def get_hwinj_segments(config,ifos,hw_inj_dir):

    os.chdir(hw_inj_dir)

    hwinjDefurl = config.get("hardware-injections", "hwinj-def-server-url")
    hwinjDefFile = config.get("hardware-injections", "hwinj-def-file")

    print("Downloading hardware injection list " + hwinjDefFile + " from " \
          + hwinjDefurl)
    hwinjDefFile, info = urllib.urlretrieve(hwinjDefurl + '/' + hwinjDefFile,
          hwinjDefFile)
    ifostr = ''
    for ifo in ifos:
        ifostr = ifostr + '--' + ifo.lower() + '-injection '

    hwinjpageCall = ' '.join(['../'+config.get("condor","hardware_inj_page"),
          "--gps-start-time",config.get("input","gps-start-time"),
          "--gps-end-time",config.get("input","gps-end-time"),
          "--segment-db",config.get("segfind","segment-url"),
          "--segment-dir","./",ifostr,
          "--source-xml",hwinjDefFile,"--get-segment-list"])

    make_external_call(hwinjpageCall)

    os.chdir("..")

#####################################################################
# Function to set lp the veto-category xml files from the vetoDefFile
def generate_veto_cat_files(config, vetoDefFile, generateVetoes):
    """
    Generate veto category xml files for each ifo using the
    veto-definer file.
    @param config: the configParser object with analysis details
    @param vetoDefFile: veto-definer xml file
    @param generateVetoes: generate the veto files
    """
    executable = config.get("condor", "segs_from_cats")
    start = config.get("input", "gps-start-time")
    end = config.get("input", "gps-end-time")

    genVetoCall = executable
    genVetoCall = ' '.join([ genVetoCall, "--cumulative-categories",
          "--segment-url", config.get("segfind", "segment-url"),
          "--veto-file", vetoDefFile,
          "--gps-start-time", start,
          "--gps-end-time", end ])

    if generateVetoes:
        print("Generating veto category xml files... this may take some time...")
        make_external_call(genVetoCall)

##############################################################################
# Function to convert a veto-category xml file to a txt file
def convert_veto_cat_xml_to_txt(config, veto_cat_file, output_file):
    """
    Converts a veto-category xml file into a veto-segments txt file.
    Note: if the veto_cat_file doesn't exist, or if there are no
    vetoes in that veto_cat_file, the output file is still written;
    it just will contain the header with no segments listed.

    @param config: the configParser object with analysis details
    @param veto_cat_file: name (and loction) of veto_cat_file to convert
    @param output_file: name (and location) of output txt file to save to
    """

    if not (os.path.isfile(veto_cat_file) and os.access(veto_cat_file,os.R_OK) ):
        print('Veto file not found or unreadable, skipping %s'%(veto_cat_file))
        return False

    print("Converting veto-category xml file %s to txt file %s" \
          %(veto_cat_file, output_file))
    sys.stdout.flush()
    # get xml-to-txt converter
    xml_to_txt_converter = config.get("condor", "ligolw_print")

    convertCall = ' '.join([ xml_to_txt_converter,
          "-t segment -c start_time -c end_time", veto_cat_file ])
    #tempfile = os.popen(convertCall)
    import subprocess
    tempfile = subprocess.Popen(convertCall, shell=True, stdout=subprocess.PIPE).stdout

    # get times from tempfile
    times = tempfile.readlines()

    tempfile.close()

    if times:
        vetoes_available = True
    else:
        vetoes_available = False

    # open output file for writing
    file = open(output_file, "w")
    file.write("# seg\tstart    \tstop     \tduration\n")

    for n,line in enumerate(times):
        (start,end) = line.split(',')
        duration = str(int(end) - int(start))
        newline = "\t".join([ str(n), start, end.strip(), duration ])
        file.write( newline + "\n" )

    file.close()

    return vetoes_available

##############################################################################
# Function to set up the segments for the analysis
def veto_segments(ifo, config, categories, generateVetoes):
    """
    generate veto segments for the given ifo

    @param ifo         : name of the ifo
    @param config      : the configParser object with analysis details
    @param categories  : list of veto categories
    @param generateVetoes: generate the veto files
    """

    start = config.getint("input","gps-start-time")
    end = config.getint("input","gps-end-time")
    vetoFiles = {}

    for category in categories:
        veto_cat_file = ifo + "-VETOTIME_CAT" + str(category) + "-" + \
          str(start) + "-" + str(end - start) + ".xml"

        vetoFile = ifo + "-CATEGORY_" + str(category) + "_VETO_SEGS-" + \
            str(start) + "-" + str(end - start) + ".txt"

        if generateVetoes:
            return_val = convert_veto_cat_xml_to_txt(config, veto_cat_file, vetoFile)
            if not return_val:
                print("No vetoes found for %s cat %i. %s will contain no segments." \
                  %(ifo, category, vetoFile))


        # if there are previous vetoes, generate combined
        if (category-1) in categories:
            combinedFile = ifo + "-COMBINED_CAT_" + str(category) + "_VETO_SEGS-" + \
                str(start) + "-" + \
                str(end - start) + ".txt"
            combinedFileXML = ifo + "-COMBINED_CAT_" + str(category) + "_VETO_SEGS-" + \
                str(start) + "-" + \
                str(end - start) + ".xml"
            vetoFiles[category] = combinedFile

            if generateVetoes:
                previousSegs = segmentsUtils.fromsegwizard(open(vetoFiles[category-1]))
                vetoSegs = segmentsUtils.fromsegwizard(open(vetoFile)).coalesce()
                vetoSegs |= previousSegs
                segmentsUtils.tosegwizard(file(combinedFile,"w"), vetoSegs)

                # create segment xml file
                xmldoc = ligolw.Document()
                xmldoc.appendChild(ligolw.LIGO_LW())
                xmldoc.childNodes[-1].appendChild(lsctables.New(lsctables.ProcessTable))
                xmldoc.childNodes[-1].appendChild(lsctables.New(lsctables.ProcessParamsTable))

                # add process table
                process = ligolw_process.append_process(xmldoc, program=__name__, version=__version__)

                gpssegs = segments.segmentlist()
                for seg in vetoSegs:
                    gpssegs.append(segments.segment(lsctables.LIGOTimeGPS(seg[0]), lsctables.LIGOTimeGPS(seg[1])))

                # add segment table
                segments_tables = ligolw_segments.LigolwSegments(xmldoc)
                segments_tables.add(ligolw_segments.LigolwSegmentList(active=gpssegs, instruments=set([ifo])))
                segments_tables.coalesce()
                segments_tables.optimize()
                segments_tables.finalize(process)

                # write file
                ligolw_process.set_process_end_time(process)
                fp = open(combinedFileXML, 'w')
                with utils.SignalsTrap():
                    utils.write_fileobj(xmldoc, fp, gz=False)
                fp.close()

        else: vetoFiles[category] = vetoFile

    return vetoFiles

##############################################################################
# Function to set up the segments for the analysis
def datafind_segments(ifo, config):
    """
    generate datafind segments for the given ifo

    @param ifo         : name of the ifo
    @param config: the configParser object with analysis details
    """
    ligoIfos = ["H1","H2","L1"]

    if ifo in ligoIfos: type = config.get("input","ligo-type")
    elif ifo == "G1": type = config.get("input","geo-type")
    elif ifo == "V1": type = config.get("input","virgo-type")

    if ifo == "V1": ifo_type = type
    else: ifo_type = ifo + "_" + type

    executable = config.get("condor", "datafind")
    start = config.getint("input","gps-start-time")
    end = config.getint("input","gps-end-time")


    dataFindFile = ifo_type + "-" + str(start) + "-" + \
        str(end - start) + ".txt"

    print("Running ligo_data_find to determine available data from " + type + \
          " frames for " + ifo)
    dataFindCall = executable
    for opt,arg in config.items("datafind"):
        dataFindCall += " --" + opt + " " + arg
    dataFindCall += " --observatory=" + ifo[0] + \
        " --type=" + ifo_type + \
        " --gps-start-time=" + str(start) + \
        " --gps-end-time=" + str(end) + " --show-times > " + \
        dataFindFile
    make_external_call(dataFindCall)
    dfSegs = segmentsUtils.fromsegwizard(file(dataFindFile)).coalesce()

    return dfSegs

##############################################################################
# Function to copy the segments files
def copyCategoryFiles(config,vetoes,directory,\
                      infile,analysisDirectory):
    """
    Copy the category files to the directory specified
    """
    outfile = analysisDirectory + "/" + directory + "/" \
              + os.path.basename(infile)
    rel_outfile = "../" + directory + "/" + os.path.basename(infile)
    if "veto-file" in vetoes:
        if infile == "":
            print(sys.stderr, "warning: " + vetoes + " left blank; proceeding "\
              "without DQ vetoes", file=sys.stderr)
            outfile += vetoes + "_BLANK.txt"
            rel_outfile += vetoes + "_BLANK.txt"
            open(outfile, "w").write("")  # touch
        else:
            shutil.copy(infile, outfile)
        config.set("segments", vetoes, rel_outfile)

##############################################################################
# Function to download the VetoDefFile
def downloadVetoDefFile(config,generate_segments):
    """
    Download the vetoDefFile from a html location
    @param config      : the configParser object with analysis details
    @param generate_segments : If False do not download, just return filename
    """

    vetoDefurl = config.get("segments", "veto-def-server-url")
    vetoDefFile = config.get("segments", "veto-def-file")

    if generate_segments:
        print("Downloading veto-definer file " + vetoDefFile + " from " \
            + vetoDefurl)
        vetoDefFile, info = urllib.urlretrieve(vetoDefurl + '/' + vetoDefFile,
            vetoDefFile)

    return vetoDefFile

##############################################################################
# Function to download the dqSegFiles
def downloadDqSegFiles(config,ifo,generate_segments):
    """
    Download the dqSegFiles from a html location
    @param config      : the configParser object with analysis details
    @param ifo         : name of the ifo
    @param generate_segments : If False do not download, just return filename
    """
    start = config.getint("input","gps-start-time")
    end = config.getint("input","gps-end-time")
    dqSegFile = ifo + "-DQ_SEGMENTS-" + str(start) + "-" + \
        str(end - start) + ".txt"
    if generate_segments:
        dq_url = config.get("segments","dq-server-url")
        dq_segdb_file = config.get("segments", ifo.lower() + '-dq-file')
        if dq_segdb_file == "":
            print("warning: no file provided to %s-dq-file; " \
              "running without data quality" % ifo.lower(), file=sys.stderr)
        else:
            print("Downloading DQ segment file " + dq_segdb_file + " from " \
                  + dq_url + " to " + dqSegFile + " ..." + sys.stdout.flush())
            dqSegFile, info = urllib.urlretrieve(dq_url + '/' + dq_segdb_file,
                  dqSegFile)
            print("done")
    return dqSegFile


##############################################################################
# Function to determine the segments to analyze
#(science segments, data quality, missing segments)
def findSegmentsToAnalyze(config, ifo, veto_categories, generate_segments=True,\
    use_available_data=False, data_quality_vetoes=False):
    """
    generate segments for the given ifo

    @param config      : the configParser object with analysis details
    @param ifo         : name of the ifo
    @param veto_categories UNDOCUMENTED
    @param generate_segments: whether the segment files should be generated
    @param use_available_data: restrict segments to data which actually available
    @param data_quality_vetoes: generate the cat2,3,4 DQ veto segments
    """

    start = config.getint("input","gps-start-time")
    end = config.getint("input","gps-end-time")

    # file names
    segFile = ifo + "-SELECTED_SEGS-" + str(start) + "-" + \
        str(end - start) + ".txt"
    segFileXML = ifo + "-SELECTED_SEGS-" + str(start) + "-" + \
        str(end - start) + ".xml"
    missedFile = ifo + "-MISSED_SEGS-" + str(start) + "-" + \
        str(end - start) + ".txt"
    missedFileXML = ifo + "-MISSED_SEGS-" + str(start) + "-" + \
        str(end - start) + ".xml"

    if generate_segments:
        print("Generating science segments for " + ifo + " ...",
              sys.stdout.flush())
    sciSegFile = science_segments(ifo, config, generate_segments)

    # generate vetoFiles
    if generate_segments:
        sciSegs = segmentsUtils.fromsegwizard(file(sciSegFile)).coalesce()
        print(" done.")
        print("Generating cat 1 veto segments for " + ifo + " ...",
              sys.stdout.flush())
        vetoFiles = veto_segments(ifo, config, [1], generate_segments)
        print("done")

        # remove cat 1 veto times
        if os.path.exists(vetoFiles[1]):
            with open(vetoFiles[1]) as f:
                vetoSegs = segmentsUtils.fromsegwizard(f).coalesce()
                sciSegs = sciSegs.__and__(vetoSegs.__invert__())

        if use_available_data:
            dfSegs = datafind_segments(ifo, config)
            analyzedSegs = sciSegs.__and__(dfSegs)
            missedSegs = sciSegs.__and__(dfSegs.__invert__())
            segmentsUtils.tosegwizard(file(missedFile,"w"), missedSegs)
            print("Writing " + ifo + " segments which cannot be analyzed to file " \
                + missedFile)
            print("Writing " + ifo + " segments which cannot be analyzed to file " \
                + missedFileXML)
            print("Not analyzing %d s, representing %.2f percent of time" %  \
               (missedSegs.__abs__(),
               100. * missedSegs.__abs__() / max(analyzedSegs.__abs__(), 0.1) ))

        else: analyzedSegs = sciSegs

        segmentsUtils.tosegwizard(file(segFile,"w"), analyzedSegs)

        # create segment xml file
        xmldoc = ligolw.Document()
        xmldoc.appendChild(ligolw.LIGO_LW())
        xmldoc.childNodes[-1].appendChild(lsctables.New(lsctables.ProcessTable))
        xmldoc.childNodes[-1].appendChild(lsctables.New(lsctables.ProcessParamsTable))

        # add process table
        process = ligolw_process.append_process(xmldoc, program=__name__, version=__version__)

        gpssegs = segments.segmentlist()
        for seg in analyzedSegs:
            gpssegs.append(segments.segment(lsctables.LIGOTimeGPS(seg[0]), lsctables.LIGOTimeGPS(seg[1])))

        # add segment table
        segments_tables = ligolw_segments.LigolwSegments(xmldoc)
        segments_tables.add(ligolw_segments.LigolwSegmentList(active=gpssegs, instruments=set([ifo])))
        segments_tables.coalesce()
        segments_tables.optimize()
        segments_tables.finalize(process)

        # write file
        ligolw_process.set_process_end_time(process)
        fp = open(segFileXML, 'w')
        with utils.SignalsTrap():
            utils.write_fileobj(xmldoc, fp, gz=False)
        fp.close()

        print("Writing " + ifo + " segments of total time " + \
            str(analyzedSegs.__abs__()) + "s to file: " + segFile)
        print("Writing " + ifo + " segments of total time " + \
            str(analyzedSegs.__abs__()) + "s to file: " + segFileXML)
        print("done")

    if data_quality_vetoes:
        print("Generating cat " + str(veto_categories) + " veto segments for " + ifo + "...")
        sys.stdout.flush()
    dqVetoes = veto_segments(ifo, config, veto_categories, data_quality_vetoes )
    if data_quality_vetoes: print("done")

    return tuple([segFile, dqVetoes])


##############################################################################
# Function to make sure we don't slide by more than the length of our data
def slide_sanity(config, playOnly = False):
    """
    check that the requested number of slides makes sense
    config    = config file
    playOnly  = only doing playground
    """
    # extract the num-slides
    try: numSlides = config.get("pipeline","num-slides")
    except: numSlides = ""

    maxLength = None
    if playOnly: maxLength = 600
    elif config.has_option("input", "max-thinca-segment"):
        maxLength = config.getint("input", "max-thinca-segment")

    if maxLength and numSlides:
        ifos = [x[0:2].upper() for x in config.options("ifo-details") \
            if "data" in x]
        maxSlide = max([config.getint("thinca-slide", ifo.lower() + "-slide") \
            for ifo in ifos])
        if (maxLength/2/maxSlide - 1) < int(numSlides):
            numSlides = str(maxLength/2/maxSlide - 1)
            print("Setting number of slides to " + numSlides + \
                  " to avoid double wrapping")

    return numSlides

##############################################################################
# Function to set up lalapps_inspiral_hipe
def hipe_setup(hipeDir, config, ifos, logPath, injSeed=None, dataFind = False, \
    tmpltBank = False, playOnly = False, vetoCat = None, vetoFiles = None, \
    dax = False, tmpltbankCache = None, local_exec_dir = None, \
    data_checkpoint = False, static_pfn_cache=None, reuse_data=False):
    """
    run lalapps_inspiral_hipe and add job to dag
    hipeDir   = directory in which to run inspiral hipe
    config    = config file
    logPath   = location where log files will be written
    injSeed   = injection file to use when running
    dataFind  = run the datafind step of the pipeline
    tmpltBank = run the template bank step (part of the datafind dag)
    vetoCat   = run this category of veto
    vetoFiles = dictionary of veto files
    tmpltbankCache = lal.Cache of template bank files
    """
    # remember what directory we are in to construct the pegasus exec dir
    ihopeDir = os.path.split(os.getcwd())[-1]

    # make the directory for running hipe
    mkdir(hipeDir)

    # create the hipe config parser, keep only relevant info
    hipecp = copy.deepcopy(config)
    if dataFind or tmpltBank: # Template generation and datafind share a number of options
        hipeSections = ["condor", "pipeline", "input", "datafind","data", \
            "ligo-data","inspiral","virgo-data", "condor-max-jobs", \
            "ligolw_add", "calibration"]
        if tmpltBank: # Template generation needs some extra options that datafind doesn't
            hipeSections.extend(["geo-data", "tmpltbank", \
                "tmpltbank-1", "tmpltbank-2", "h1-tmpltbank", "h2-tmpltbank", \
                "l1-tmpltbank", "v1-tmpltbank", "g1-tmpltbank"])
    elif vetoCat > 1:
        hipeSections = ["condor", "pipeline", "input", "data", "datafind", \
            "ligo-data", "virgo-data", "geo-data", "calibration", "tmpltbank", \
            "inspiral", "veto-inspiral", "g1-inspiral", "h1-inspiral", \
            "h2-inspiral", "l1-inspiral", "v1-inspiral", "ligolw_add", \
            "ligolw_cafe", "thinca", "thinca-2", "thinca-slide", "coire", \
            "coire-1", "coire-2", "coire-inj", "sire", "sire-inj", "condor-max-jobs"]
        if vetoCat == 4:
            hipeSections.extend(["h1-inspiral-coherent","h2-inspiral-coherent",
              "l1-inspiral-coherent","v1-inspiral-coherent","cohbank","trigbank-coherent",
              "inspiral-coherent", "chia", "cohinspbank", "chia-inj", "cohire", "cohire-inj"])
    else:
        hipeSections = ["condor", "pipeline", "input", "data", "datafind",\
            "ligo-data", "virgo-data", "geo-data", "calibration", "tmpltbank", \
            "tmpltbank-1", "tmpltbank-2", "h1-tmpltbank", "h2-tmpltbank", \
            "l1-tmpltbank", "v1-tmpltbank", "g1-tmpltbank", "no-veto-inspiral", \
            "veto-inspiral", "inspiral", "h1-inspiral", "h2-inspiral", \
            "l1-inspiral", "g1-inspiral", "v1-inspiral", "ligolw_add", \
            "ligolw_cafe", "thinca", "thinca-1", "thinca-2", "thinca-slide", \
            "trigbank", "sire", "sire-inj", "coire", "coire-1", "coire-2", \
            "coire-inj", "condor-max-jobs"]

    for seg in hipecp.sections():
        if not seg in hipeSections: hipecp.remove_section(seg)

    hipecp.remove_option("condor","hipe")
    hipecp.remove_option("condor","plot")
    hipecp.remove_option("condor","follow")

    # add the inspinj section
    hipecp.add_section("inspinj")

    # set the data type
    if playOnly:
        hipecp.set("pipeline", "playground-data-mask", "playground_only")
    else:
        hipecp.set("pipeline", "playground-data-mask", "all_data")

    # set the user-tag to be the same as the directory name
    if hipecp.get("pipeline","user-tag"):
        usertag = hipecp.get("pipeline", "user-tag") + "_" + hipeDir.upper()
    else:
        usertag = hipeDir.upper()

    # template banks should not have a usertag so they can be
    # passed between sub-dags more easily
    if hipeDir == "datafind":
        hipecp.remove_option("pipeline","user-tag")
    else:
        hipecp.set("tmpltbank","user-tag","")

    if vetoCat:
        if vetoCat > 1:
            # set the old usertag in inspiral and inspinj,
            # so that we pick up the correct xml inputs
            sections = ["inspiral"]
            for section in ["inspiral","inspinj"]:
                hipecp.set(section, "user-tag",usertag)

            # set the correct pipeline usertag
            usertag += "_CAT_" + str(vetoCat) + "_VETO"

            for cat in range(2, vetoCat+1):
                section = "coire-2-cat-" + str(cat)
                if config.has_section(section):
                    for opt, arg in config.items(section):
                        hipecp.set("coire-2",opt,arg)

        # add the veto files in the thinca section
        if hipecp.has_section("thinca-2"):
            # add the veto files in the thinca section
            for ifo in ifos:
                hipecp.set("thinca-2", ifo.lower() + "-veto-file", vetoFiles[ifo][vetoCat])
            hipecp.set("thinca-2", "do-veto", "")
        elif config.has_option("hipe-arguments","ringdown"):
            # add the veto files in the thinca section
            for ifo in ifos:
                hipecp.set("thinca", ifo.lower() + "-veto-file", vetoFiles[ifo][vetoCat])
            hipecp.set("thinca", "do-veto", "")
        else:
            # add a vetoes section
            hipecp.add_section("vetoes")
            hipecp.set("vetoes", "vetoes-file", vetoFiles["combined-veto-file"][vetoCat])
            hipecp.set("vetoes", "vetoes-name", "VETO_CAT%i_CUMULATIVE"%(vetoCat))

    # set the usertag
    if hipeDir == "datafind":
        # we don't set a usertag for the datafind directory
        pass
    else:
        hipecp.set("pipeline", "user-tag",usertag)

    if injSeed:
        # fail if the seed is not identical to its integer form
        # to prevent later problems with inspinj
        if not ( str(int(injSeed)) == injSeed ):
            print("Injection seed: " + injSeed + "\n", file=sys.stderr)
            print("Error: the injection seed must be an integer without "
                  "leading zeros! Exiting...", file=sys.stderr)
            sys.exit(1)

        # copy over the arguments from the relevant injection section
        for (name,value) in config.items(hipeDir):
            hipecp.set("inspinj",name,value)
        if config.has_option("hipe-arguments","ringdown"):
            injType = config.get(hipeDir,"injection-type")
            hipecp.set("inspiral","injection-type",injType)
            hipecp.remove_option("inspinj","injection-type")
            if injType == "RINGDOWN":
                executable = "../executables/lalapps_rinj"
                hipecp.set("condor","inspinj",executable)
        hipecp.remove_section(hipeDir)
        hipecp.set("input","injection-seed",injSeed)
        hipecp.set("input", "num-slides", "0")
        # set any extra inspiral arguments for the injection
        try:
            for item in config.items('-'.join([hipeDir,"inspiral"])):
                hipecp.set("inspiral",item[0],item[1])
        except NoSectionError:
            pass
    else:
        # add the time slide to the ini file
        numSlides = slide_sanity(config, playOnly)
        hipecp.set("input","num-slides",numSlides)

    # return to the directory, write ini file and run hipe
    os.chdir(hipeDir)
    iniFile = "inspiral_hipe_"
    iniFile += hipeDir
    if vetoCat > 1:
        iniFile += "_cat" + str(vetoCat) + "_veto"
    iniFile += ".ini"

    hipecp.write(file(iniFile,"w"))

    print("Running hipe in directory " + hipeDir)
    if dataFind or tmpltBank:
        print("Running datafind / template bank generation")
    elif injSeed:
        print("Injection seed: " + injSeed)
    else:
        print("No injections, " + str(hipecp.get("input","num-slides")) + \
              " time slides")
    if vetoCat:
        print("Running the category " + str(vetoCat) + " vetoes")
    print()

    # work out the hipe call:
    hipeCommand = config.get("condor","hipe")
    hipeCommand += " --log-path " + logPath
    hipeCommand += " --config-file " + iniFile
    if dax:
        hipeCommand += " --dax "
    if playOnly: hipeCommand += " --priority 10"
    for item in config.items("ifo-details"):
        hipeCommand += " --" + item[0] + " " + item[1]
    # describes the ckpt if statement
    if data_checkpoint:
        hipeCommand += " --data-checkpoint "

    def test_and_add_hipe_arg(hipeCommand, hipe_arg):
        if config.has_option("hipe-arguments",hipe_arg):
            hipeCommand += "--" + hipe_arg + " " + \
              config.get("hipe-arguments",hipe_arg)
        return(hipeCommand)

    if dataFind or tmpltBank:
        if dataFind:
            for hipe_arg in ["datafind","ringdown","write-script"]:
                hipeCommand = test_and_add_hipe_arg(hipeCommand,hipe_arg)
        if tmpltBank:
            for hipe_arg in ["template-bank","ringdown","write-script"]:
                hipeCommand = test_and_add_hipe_arg(hipeCommand,hipe_arg)
    elif vetoCat > 1:
        if config.has_option("hipe-arguments","ringdown"):
            hipe_args = ["coincidence", "ringdown","coire-coincidence",
              "summary-first-coinc-triggers","write-script"]
        elif config.has_option("hipe-arguments","second-coinc"):
            hipe_args = ["second-coinc", "coire-second-coinc",
              "summary-coinc-triggers", "sire-second-coinc",
              "summary-single-ifo-triggers","write-script"]
            if vetoCat == 4 and config.has_option("hipe-arguments","coherent-bank"):
                hipe_args.extend(["coherent-bank","coherent-inspiral","cohire",
                  "summary-coherent-inspiral-triggers"])
        elif vetoCat == 4 and config.has_option("hipe-arguments","coherent-bank"):
            hipe_args = ["coherent-bank","coherent-inspiral","cohire",
              "summary-coherent-inspiral-triggers","write-script"]
        else:
            hipe_args = ["coincidence","write-script"]
        for hipe_arg in hipe_args:
            hipeCommand = test_and_add_hipe_arg(hipeCommand,hipe_arg)
    else:
        omit = ["datafind", "template-bank", "disable-dag-categories", "disable-dag-priorities"]
        if config.has_option("hipe-arguments","coherent-bank"):
            omit.extend(["coherent-bank","coherent-inspiral","cohire",
              "summary-coherent-inspiral-triggers"])
        for (opt, arg) in config.items("hipe-arguments"):
            if opt not in omit:
                hipeCommand += "--" + opt + " " + arg

    hipeCommand = test_and_add_hipe_arg(hipeCommand,"disable-dag-categories")
    hipeCommand = test_and_add_hipe_arg(hipeCommand,"disable-dag-priorities")

    # run lalapps_inspiral_hipe
    make_external_call(hipeCommand)

    # link datafind
    if not dax and not dataFind and not tmpltBank and not vetoCat > 1:
        try:
            os.rmdir("cache")
            os.symlink("../datafind/cache", "cache")
        except: pass

    # symlink in the template banks needed by the inspiral jobs
    if not dax and tmpltbankCache:
        symlinkedCache = symlink_tmpltbank(tmpltbankCache, hipeDir)

    iniBase = iniFile.rstrip("ini")
    if hipeDir == "datafind":
        hipeDag = iniBase + "dag"
        hipeDax = iniBase + "dax"
    else:
        hipeDag = iniBase + usertag + ".dag"
        hipeDax = iniBase + usertag + ".dax"

    hipeJob = pipeline.CondorDAGManJob(hipeDag, hipeDir, hipeDax)
    hipeNode = pipeline.CondorDAGManNode(hipeJob)

    # pass a static pfn cache if we are given one
    if static_pfn_cache:
        hipeNode.set_static_pfn_cache(static_pfn_cache)

    # set the reuse data option
    if reuse_data:
        hipeNode.set_reduce_dax = reuse_data

    # grab the tmpltbank from hipecp file if provided in input section
    if hipecp.has_section("input") and \
        hipecp.has_option("input", "fixed-bank"):
        tmpltbankfile = hipecp.get("input", "fixed-bank")
        hipeJob.add_pfn_cache(hipe_pfn_list_cache( [tmpltbankfile] ))

    # add the maxjob categories to the dagman node class
    # FIXME pegasus should handle this in the dax schema itself
    for cp_opt in config.options('condor-max-jobs'):
        hipeNode.add_maxjobs_category(cp_opt,config.getint('condor-max-jobs',cp_opt))

    # collapse the short running jobs in all of the dags
    hipeNode.set_cluster_jobs('horizontal')

    # tell pegasus where ihope wants us to run the jobs
    hipeJob.set_pegasus_exec_dir(os.path.join(
      local_exec_dir, '/'.join(os.getcwd().split('/')[-1:])))

    if hipeDir == "datafind":
        hipeNode.add_output_file( hipe_cache(ifos, None, \
            hipecp.getint("input", "gps-start-time"), \
            hipecp.getint("input", "gps-end-time")) )
    else:
        hipeNode.set_user_tag(usertag)
        hipeNode.add_output_file( hipe_cache(ifos, usertag, \
            hipecp.getint("input", "gps-start-time"), \
            hipecp.getint("input", "gps-end-time")) )

    hipeJob.add_pfn_cache(hipe_pfn_glob_cache('../segments/*txt'))
    hipeJob.add_pfn_cache(hipe_pfn_glob_cache('../segments/*VETOTIME*xml'))

    if hipecp.has_section("ligolw_cafe"):
        num_slide_files = int(hipecp.get("ligolw_cafe", "num-slides-files"))
        slide_files = [hipecp.get("ligolw_cafe","slides-file-%i"%idx)
            for idx in range(num_slide_files)]
        hipeJob.add_pfn_cache(hipe_pfn_list_cache(slide_files))
        hipeJob.add_pfn_cache(hipe_pfn_glob_cache('*CAFE*.cache'))

    # return to the original directory
    os.chdir("..")

    return hipeNode

##############################################################################
# Function to remove meta option from lalapps_plot_hipe ini file
def remove_plot_meta_option(plotcp,patternType,plottingCode):
    programTag = patternType + '-program-tag'
    plotMeta = plottingCode + '-meta'
    plotcp.remove_option(plotMeta,programTag)
    patterns = (plotcp.get(plotMeta,'cache-patterns')).split(',')
    patterns.remove(patternType)
    newpatterns=','.join(patterns)
    plotcp.set(plotMeta,'cache-patterns',newpatterns)

##############################################################################
# Function to set up lalapps_plot_hipe
def plot_setup(plotDir, config, logPath, stage, injectionSuffix,
    zerolagSuffix, slideSuffix, bankSuffix, cacheFile, injdirType, tag = None,
    ifos = None, cat = 3):
    """
    run lalapps_plot_hipe and add job to dag
    plotDir   = directory in which to run inspiral hipe
    config    = config file
    logPath   = location where log files will be written
    stage     = which stage to run (first, second or both)
    injectionSuffix = the string to restrict to for injections
    zerolagSuffix   = the string to restrict to for zero lag
    slideSuffix     = the string to restrict to for time slides
    bankSuffix      = the string to restrict to for bank plots
    cacheFile       = the input cache file for plotting
    tag             = extra tag for naming
    """
    # make the directory for running hipe
    mkdir(plotDir)

    plotcp = copy.deepcopy(config)

    # set details for the common section
    plotcp.add_section("common")
    plotcp.set("common","gps-start-time", plotcp.get("input","gps-start-time") )
    plotcp.set("common","gps-end-time", plotcp.get("input","gps-end-time") )
    plotcp.set("common","output-path", ".")
    plotcp.set("common","enable-output","")

    plotSections = ["common", "pipeline", "condor",\
        "plotinspiral", "plotinspiral-meta", \
        "plotthinca", "plotthinca-meta", \
        "plotnumtemplates", "plotnumtemplates-meta", \
        "plotethinca", "plotethinca-meta", \
        "plotinspmissed", "plotinspmissed-meta", \
        "plotinspinj", "plotinspinj-meta", \
        "plotsnrchi", "plotsnrchi-meta", \
        "plotinspiralrange", "plotinspiralrange-meta", \
        "ploteffdistcut", "ploteffdistcut-meta"]

    for seg in plotcp.sections():
        if not seg in plotSections: plotcp.remove_section(seg)

    plotcp.remove_option("condor","hipe")
    plotcp.remove_option("condor","plot")
    plotcp.remove_option("condor","follow")

    # XXX Can't yet run the plotting codes in standard universe
    if plotcp.get("condor","universe") == "standard":
        plotcp.set("condor","universe","vanilla")

    # set the various suffixes in pipeline
    plotcp.set("pipeline","injection-suffix",injectionSuffix)
    plotcp.set("pipeline","inj-suffix",injectionSuffix)
    plotcp.set("pipeline","found-suffix",injectionSuffix)
    plotcp.set("pipeline","missed-suffix",injectionSuffix)
    plotcp.set("pipeline","trigbank-suffix",bankSuffix)
    plotcp.set("pipeline","zerolag-suffix",zerolagSuffix)
    plotcp.set("pipeline","trig-suffix",zerolagSuffix)
    plotcp.set("pipeline","coinc-suffix",zerolagSuffix)
    plotcp.set("pipeline","slide-suffix",slideSuffix)

    numSlides = slide_sanity(config, ("PLAYGROUND" in slideSuffix ))
    plotcp.set("pipeline","num-slides", numSlides)

    # Adding followup options to plotinspmissed
    analysisstart = plotcp.get("common","gps-start-time")
    analysisend = plotcp.get("common","gps-end-time")
    analysisduration = int(analysisend) - int(analysisstart)
    inspmissedVetoDir = "../segments"
    for ifo in ifos:
        if cat == 2:
            plotcp.set("plotinspmissed","followup-vetofile-" + ifo.lower(),
              inspmissedVetoDir + "/" + ifo + "-CATEGORY_" + str(cat) +
              "_VETO_SEGS-" + analysisstart + "-" + str(analysisduration) + ".txt")
        else:
            plotcp.set("plotinspmissed","followup-vetofile-" + ifo.lower(),
              inspmissedVetoDir + "/" + ifo + "-COMBINED_CAT_" + str(cat) +
              "_VETO_SEGS-" + analysisstart + "-" + str(analysisduration) + ".txt")

    # Adding followup option to plotinspmissed
    plotcp.set("plotinspmissed","followup-tag",injdirType)

    # Remove options if no slide or zero lag files are available.
    if "NONE_AVAILABLE" in slideSuffix:
        if plotcp.has_option('plotsnrchi-meta','slide-program-tag'):
            remove_plot_meta_option(plotcp,'slide','plotsnrchi')
        if plotcp.has_option('ploteffdistcut-meta','slide-program-tag'):
            remove_plot_meta_option(plotcp,'slide','ploteffdistcut')
        if plotcp.has_option('plotethinca-meta','slide-program-tag'):
            remove_plot_meta_option(plotcp,'slide','plotethinca')
    if "NONE_AVAILABLE" in zerolagSuffix:
        if plotcp.has_option('plotsnrchi-meta','trig-program-tag'):
            remove_plot_meta_option(plotcp,'trig','plotsnrchi')

    # set the user-tag
    if plotcp.get("pipeline","user-tag"):
        usertag = plotcp.get("pipeline","user-tag")
        plotcp.set("pipeline","input-user-tag",usertag)
        usertag += plotDir.upper()
    else:
        usertag = plotDir.upper()
        plotcp.set("pipeline","input-user-tag","")

    if tag: usertag += "_" + tag
    plotcp.set("pipeline","user-tag",usertag)

    plotcp.set("common","cache-file",cacheFile)

    # return to the directory, write ini file and run hipe
    os.chdir(plotDir)
    iniFile = "plot_hipe_"
    iniFile += plotDir
    if tag: iniFile += "_" + tag.lower()
    iniFile += ".ini"

    plotcp.write(file(iniFile,"w"))

    print("Running plot hipe in directory " + plotDir)
    print("Using zero lag sieve: " + zerolagSuffix)
    print("Using time slide sieve: " + slideSuffix)
    print("Using injection sieve: " + injectionSuffix)
    print("Using bank sieve: " + bankSuffix)
    print()

    # work out the hipe call:
    plotCommand = config.get("condor","plot")
    plotCommand += " --log-path " + logPath
    plotCommand += " --config-file " + iniFile
    plotCommand += " --priority 10"
    for item in config.items("ifo-details"):
        plotCommand += " --" + item[0] + " " + item[1]

    for item in config.items("plot-arguments"):
        plotCommand += " --" + item[0] + " " + item[1]

    if stage == "first" or stage == "both":
        plotCommand += " --first-stage"
    if stage == "second" or stage == "both":
        plotCommand += " --second-stage"

    # run lalapps_inspiral_hipe
    make_external_call(plotCommand)

    # make hipe job/node
    plotDag = iniFile.rstrip("ini") + usertag + ".dag"
    plotJob = pipeline.CondorDAGManJob(plotDag, plotDir)
    plotNode = pipeline.CondorDAGManNode(plotJob)
    plotNode.set_user_tag(usertag)

    # return to the original directory
    os.chdir("..")

    return plotNode

##############################################################################
# Function to set up and run pipedown

def pipedownSetup(dag,config,logPath,pipedownDir,\
                  cacheFile,parentNodes,playgroundOnly, \
                  optiontorunmvsc):
    """
    Set up and run pipedown
    dag      = the dag
    config   = the config object
    logPath  = Location where log files will be written
    pipedownDir = Directory to run pipedown in
    cacheFile = The input ihope cache file
    parentNodes = Name of parent dags to add
    playgroundOnly = Are we opening the box?
    """
    # Get necessary information from the config object
    gpsStart = config.get("input","gps-start-time")
    gpsEnd = config.get("input","gps-end-time")

    # Make directory

    mkdir(pipedownDir)
    os.chdir(pipedownDir)

    # Create the necessary ini file

    pipeCp = copy.deepcopy(config)
    pipeCp.set("condor","universe","vanilla")

    # Create input section

    pipeCp.remove_section("input")
    pipeCp.add_section("input")
    pipeCp.set("input","ihope-segments-directory","../segments")

    # Write ini file to folder

    iniFile = pipedownDir
    iniFile += ".ini"

    pipeCp.write(file(iniFile,"w"))

    # Set up the command to run pipedown

    pipeCommand = config.get("condor","pipedown")
    pipeCommand += " --log-path " + logPath
    pipeCommand += " --config-file " + iniFile
    pipeCommand += " --gps-start-time " + gpsStart
    pipeCommand += " --gps-end-time " + gpsEnd
    pipeCommand += " --ihope-cache " + cacheFile
    if optiontorunmvsc:
        pipeCommand += " --run-mvsc"
    if not playgroundOnly:
        pipeCommand += " --generate-all-data-plots"
        # Need to add playground command

    # run lalapps_pipedown
    make_external_call(pipeCommand)

    # make pipedown job/node
    pipeDag = iniFile.rstrip("ini") + "dag"
    pipeJob = pipeline.CondorDAGManJob(pipeDag, pipedownDir)
    pipeNode = pipeline.CondorDAGManNode(pipeJob)
    dag.add_node(pipeNode)
    if parentNodes:
        for thisDag in parentNodes:
            pipeNode.add_parent(thisDag)

    # return to the original directory
    os.chdir("..")

    return dag



##############################################################################
# Function to set up zero/slide plots:
def zeroSlidePlots(dag, plotDir, config, logPath, zerolagSuffix, slideSuffix,
    cacheFile, doDagCategories, parentDags = None, vetoParentDags = None,
    vetoCat = [2,3,4], ifos = []):
    """
    set up plots for zero lag and time slides
    dag       = the dag
    plotDir   = directory in to set up plots
    config    = config file
    logPath   = location where log files will be written
    zerolagSuffix   = the string to restrict to for zero lag
    slideSuffix     = the string to restrict to for time slides
    cacheFile       = the input cache file for plotting
    doCategories    = dag categories turned on
    parentDags      = the name of the parent dag to add
    vetoParentDags  = the name of the veto parent dag to add
    vetoCat         = veto category
    """
    # first stage
    plotcp = copy.deepcopy(config)
    plotcp.add_section("plot-arguments")
    plotcp.set("plot-arguments","plotinspiralrange","")
    plotcp.set("plot-arguments","plotnumtemplates","")
    plotcp.set("plot-arguments","plotinspiral","")
    plotcp.set("plot-arguments","plotthinca","")
    plotcp.set("plot-arguments","write-script","")
    # Add plotthinca scaling option for zero lag play full data slide plots
    if zerolagSuffix == "PLAYGROUND" and slideSuffix == "FULL_DATA":
        plotcp.set("plotthinca","zero-lag-playground","")

    plotNode = plot_setup(plotDir, plotcp, logPath, "first", \
        "", zerolagSuffix, slideSuffix, slideSuffix, cacheFile, "", ifos=ifos)
    if doDagCategories:
        plotNode.set_category('plotting')
    dag.add_node(plotNode)
    if parentDags:
        for thisDag in parentDags:
            plotNode.add_parent(thisDag)

    return dag

##############################################################################
# Function to set up inj/zero/slide plots:
def injZeroSlidePlots(dag, plotDir, config, logPath, injectionSuffix,
    zerolagSuffix, slideSuffix, cacheFile, doDagCategories, parentDags = None,
    vetoParentDags = None, vetoCat = [2,3,4], ifos = []):
    """
    set up plots for injections, zero lag and time slides
    dag       = the dag
    plotDir   = directory in to set up plots
    config    = config file
    logPath   = location where log files will be written
    injectionSuffix = the string to restrict to for injections
    zerolagSuffix   = the string to restrict to for zero lag
    slideSuffix     = the string to restrict to for time slides
    cacheFile       = the input cache file for plotting
    doCategories    = dag categories turned on
    parentDags      = the name of the parent dag to add
    vetoParentDags  = the name of the veto parent dag to add
    vetoCat         = veto category
    """
    # first stage
    plotcp = copy.deepcopy(config)
    plotcp.add_section("plot-arguments")
    plotcp.set("plot-arguments","plotinspinj","")
    plotcp.set("plot-arguments","plotinspmissed","")
    plotcp.set("plot-arguments","ploteffdistcut","")
    plotcp.set("plot-arguments","write-script","")
    injPlotNode = plot_setup( plotDir, \
        plotcp, logPath, "first", injectionSuffix, zerolagSuffix, \
        slideSuffix, injectionSuffix, cacheFile, injectionSuffix, ifos=ifos )
    if doDagCategories:
        injPlotNode.set_category('plotting')
    dag.add_node(injPlotNode)

    if parentDags:
        for thisDag in parentDags:
            injPlotNode.add_parent(thisDag)

    return dag

##############################################################################
# Function to set up a HW inj page job
def hwinj_page_setup(cp,ifos,veto_categories,hw_inj_dir):
    """
    run ligolw_cbc_hardware injection page, soring the input and output in
    the subdirectory hardware_injection_summary
    """
    hwInjNodes = []
    hwinj_length = cp.getint("input","gps-end-time") - cp.getint("input","gps-start-time")

    hwInjJob = inspiral.HWinjPageJob(cp)

    veto_categories.append(None)
    for veto in veto_categories:
        if cp.get("pipeline","user-tag"):
            usertag = cp.get("pipeline", "user-tag") + "_" + "FULL_DATA"
        else:
            usertag = "FULL_DATA"

        if veto:
            usertag += "_CAT_" + str(veto) + "_VETO"

        cacheFile = hipe_cache( ifos, usertag, cp.getint("input", "gps-start-time"), cp.getint("input", "gps-end-time") )

        if not os.path.isfile(os.path.join("full_data", cacheFile)):
            print("WARNING: Cache file FULL_DATA/" + cacheFile, file=sys.stderr)
            print("does not exist! This might cause later failures.", file=sys.stderr)

        outfilename = os.path.join(hw_inj_dir, ''.join(ifos) + '-HWINJ_SUMMARY')
        if veto:
            outfilename += '_CAT_' + str(veto)
        outfilename += '-' + cp.get("input","gps-start-time") + '-' + str(hwinj_length) + '.html'

        hwInjNode = inspiral.HWinjPageNode(hwInjJob)
        hwInjNode.set_start(cp.get("input","gps-start-time"))
        hwInjNode.set_end(cp.get("input","gps-end-time"))

        hwInjNode.set_input_cache(os.path.join('full_data', cacheFile))
        hwInjNode.set_cache_string('*SIRE_FIRST*')

        hwInjNode.set_source_xml(os.path.join(hw_inj_dir,cp.get("hardware-injections", "hwinj-def-file")))
        hwInjNode.set_segment_dir(hw_inj_dir)
        hwInjNode.set_output_file(outfilename)

        hwInjNode.add_var_opt('analyze-injections','')
        for ifo in ifos:
            hwInjNode.add_var_opt(ifo.lower()+'-injections','')

        hwInjNodes.append(hwInjNode)

    return hwInjNodes

##############################################################################
# Function to set up lalapps_followup_pipe
def followup_setup(followupDir, config, opts, hipeDir):
    """
    run lalapps_followup_pipe and add job to dag
    followupDir = directory to output the followup
    config    = config file
    """

    # make the directory for followup pipe
    mkdir(followupDir)

    # create the followup config parser, keep only relevant info
    followupcp = copy.deepcopy(config)
    followupSections = ["condor", "hipe-cache", "triggers", "datafind", \
        "q-datafind", "qscan", "q-hoft-datafind", "qscan-hoft", \
        "plots", "output", "seg"]

    for seg in followupcp.sections():
        if not seg in followupSections: followupcp.remove_section(seg)

    followupcp.remove_option("condor","hipe")
    followupcp.remove_option("condor","follow")

    # XXX this should be replaced by getting the information from the hipe cache
    # set the cache paths
    followupcp.add_section("hipe-cache")
    followupcp.set("hipe-cache", "hipe-cache-path", "hipe_cache")
    followupcp.set("hipe-cache", "science-run", "S5")

    for path in ["tmpltbank-path", "trigbank-path", "first-inspiral-path", \
         "first-coinc-path"]:
        followupcp.set("hipe-cache", path, "../" + hipeDir)

    # set the xml-glob
    followupcp.set("triggers", "xml-glob", "../" + hipeDir + "/*COIRE*H*xml")
    # to here XXX

    # correct paths to qscan config files
    for section in ["qscan", "qscan-hoft"]:
        for (opt, arg) in followupcp.items(section):
            if "config-file" in opt and arg[0] != "/":
                arg = "../../" + arg
                followupcp.set(section, opt, arg)

    # return to the directory, write ini file and run hipe
    os.chdir(followupDir)
    iniFile = "followup_pipe_" + followupDir + ".ini"
    followupcp.write(file(iniFile,"w"))

    # link datafind output from original hipe
    try: os.symlink("../datafind/cache", "hipe_cache")
    except: pass
    print("Running followup pipe in directory " + followupDir)

    # work out the followup_pipe call:
    followupCommand = config.get("condor","follow")
    followupCommand += " --log-path " + opts.log_path
    followupCommand += " --config-file " + iniFile

    for item in config.items("followup-arguments"):
        followupCommand += " --" + item[0] + " " + item[1]

    # set up a fake followup dag -- the real one can't be generated until the
    # analysis is done
    followupDag = iniFile.rstrip("ini") + "dag"
    f = open(followupDag,"w")
    f.write("\n")
    f.close()

    # add job to dag
    followupJob = pipeline.CondorDAGManJob(followupDag, followupDir)
    followupNode = pipeline.CondorDAGManNode(followupJob)


    # write the pre-script to run lalapps_followup_pipe at the appropriate time
    f = open(followupDag + ".pre","w")
    f.write("#! /bin/bash\n")
    f.write("cd followup\n")
    f.write(followupCommand)
    f.write("cd ..\n")
    f.close()
    os.chmod(followupDag + ".pre", 0o744)
    followupNode.set_pre_script(followupDir + "/" + followupDag + ".pre")

    # return to the original directory
    os.chdir("..")

    return followupNode

##############################################################################
# plot_hipe helpers

def determine_sieve_patterns(cp, plot_name, ifotag, usertag=None):
    """
    From the given plot configuration file, determine the --*-pattern values
    for a given plot name.  Returns as a dictionary of command-line option-
    value pairs.

    Example ini file:
    [pipeline]
    trig-suffix = PLAYGROUND
    bank-suffix = PLAYGROUND
    missed-suffix = INJ*

    [plotinspiral-meta]
    cache-patterns = trig,bank,missed
    trig-program-tag = SIRE
    bank-program-tag = TRIGBANK
    missed-program-tag = SIRE

    Example invocation:
    >>> from lalapps.inspiralutils import determine_sieve_patterns
    >>> from six.moves.configparser import ConfigParser
    >>> cp = ConfigParser()
    >>> cp.read("plot_hipe.ini")
    >>> print(determine_sieve_patterns(cp, "plotinspiral", "H1"))
    {'bank-pattern': 'TRIGBANK_H1*_PLAYGROUND', 'missed-pattern':
     'SIRE_H1*_INJ*', 'trig-pattern': 'SIRE_H1*_PLAYGROUND'}

    """
    meta_name = plot_name + "-meta"
    usertag = cp.get("pipeline", "input-user-tag") or usertag
    suffixes = dict([(opt[:-len("-suffix")],val) for opt,val in cp.items("pipeline") if opt.endswith("-suffix")])

    patterns = {}
    for pattern_name in cp.get(meta_name, "cache-patterns").split(","):
        program_tag = cp.get(meta_name, pattern_name + "-program-tag")
        pattern = program_tag + "_" + ifotag + "*"
        if usertag is not None:
            pattern += "_" + usertag
        suffix = suffixes.get(pattern_name)
        if suffix is not None:
            pattern += "_" + suffix
        patterns[pattern_name + "-pattern"] = pattern
    return patterns

###############################################################################
# This function will set up the necessary omega scan inputs in the current
# directory

def omega_scan_setup(cp,ifos):
    cp.set('omega-scans','do-omega-scan','')
    cp.set('omega-scans','omega-executable',cp.get('condor','omegascan'))
    print("Beginning set up of omega scan directory.")
    start = cp.get("input","gps-start-time")
    end = cp.get("input","gps-end-time")
    # First we set up the configuration files
    sampleFrequency = cp.get('omega-setup','sample-frequency')
    searchTRange = cp.get('omega-setup','search-time-range')
    searchFRange = cp.get('omega-setup','search-frequency-range')
    searchFRange = searchFRange.split(',')
    searchFRange = '[' + searchFRange[0] + ' ' + searchFRange[1] + ']'
    searchQRange = cp.get('omega-setup','search-q-range')
    searchQRange = searchQRange.split(',')
    searchQRange = '[' + searchQRange[0] + ' ' + searchQRange[1] + ']'
    searchMaxELoss = cp.get('omega-setup','search-max-energy-loss')
    whiteNoiseFAR = cp.get('omega-setup','white-noise-far')
    searchWindowDuration = cp.get('omega-setup','search-window-duration')
    plotNormERange = cp.get('omega-setup','plot-normalized-energy-range')
    plotNormERange = plotNormERange.split(',')
    plotNormERange = '[' + plotNormERange[0] + ' ' + plotNormERange[1] + ']'
    plotTimeRanges = cp.get('omega-setup','plot-time-ranges')
    plotTimeRanges = plotTimeRanges.split(',')
    plotTimeRanges = ' '.join(plotTimeRanges)
    plotTimeRanges = '[' + plotTimeRanges + ']'
    configFileText = []
    configFileHead = []
    configFileHead.append('''[Context,Context]''')
    configFileHead.append('')
    configFileHead.append('''[Parameters,Parameter Estimation]''')
    configFileHead.append('')
    configFileHead.append('''[Notes,Notes]''')
    configFileHead.append('')
    configFileHead.append('''[Omega scans,Omega scans]''')
    configFileHead.append('')
    configFileText.append('{')
    #This will be the line for the channel this is configFileText[1]
    configFileText.append('')
    #This will be the line for the frame this is configFileText[2]
    configFileText.append('')
    configFileText.append('  sampleFrequency:  ' + sampleFrequency)
    configFileText.append('  searchTimeRange:  ' + searchTRange)
    configFileText.append('  searchFrequencyRange:  ' + searchFRange)
    configFileText.append('  searchQRange:  ' + searchQRange)
    configFileText.append('  searchMaximumEnergyLoss:  ' + searchMaxELoss)
    configFileText.append('  whiteNoiseFalseRate:  ' + whiteNoiseFAR)
    configFileText.append('  alwaysPlotFlag:  1')
    configFileText.append('  searchWindowDuration:  ' + searchWindowDuration)
    configFileText.append('  plotTimeRanges:  ' + plotTimeRanges)
    configFileText.append('  plotFrequencyRange:  []')
    configFileText.append('  plotNormalizedEnergyRange:  ' + plotNormERange)
    configFileText.append('}')
    configFileText.append('')
    typeNames = {}
    for ifo in ifos:
        ifo = ifo.upper()
        if ifo == 'H1' or ifo == 'H2' or ifo == 'L1':
            channel_name = 'ligo-channel'
            type_name = 'ligo-type'
        elif ifo == 'V1':
            channel_name = 'virgo-channel'
            type_name = 'virgo-type'
        else:
            print("IFO " + ifo + " is not yet supported for omega scans in ihope")
            continue
        if cp.has_option('omega-setup',channel_name):
            channels = (cp.get('omega-setup',channel_name)).split(',')
            type = cp.get('omega-setup',type_name)
        else:
            channels = (cp.get('input',channel_name)).split(',')
            type = cp.get('input',type_name)
        if ifo == 'H1' or ifo == 'H2' or ifo == 'L1':
            type = ifo + '_' + type
        typeNames[ifo] = type
        configFileOut = '\n'.join(configFileHead)
        configFileOut += '\n'
        for channel in channels:
            if channel[0:3] == 'PEM' or channel[0:3] == 'GDS':
                configFileText[1] = "  channelName:  '" + ifo[0]+'0:'+channel + "'"
            else:
                configFileText[1] = "  channelName:  '" + ifo+':'+channel + "'"
            configFileText[2] = "  frameType:  '" + type + "'"
            configFileOut += '\n'.join(configFileText)
            configFileOut +='\n'
        outFile = open(ifo + '_omega_config.txt','w')
        outFile.write(configFileOut)
        outFile.close()
        cp.set('omega-scans',ifo.lower() + '-omega-config-file','../omega_setup/'+ifo + '_omega_config.txt')

    print("Created omega scan configuration files")

    # And we need to create the necessary frame caches
    if not os.path.isdir('cache'):
        os.mkdir('cache')
    dataFindCall = []
    dataFindCall.append(cp.get('condor','datafind'))
    for opt,val in cp.items('datafind'):
        if opt == 'gaps':
            continue
        dataFindCall.append('--' + opt + ' ' + val)
    dataFindCall.append('--lal-cache')
    dataFindCall.append('--gps-start-time ' + start)
    dataFindCall.append('--gps-end-time ' + end)
    for ifo in ifos:
        ifo = ifo.upper()
        ifoSpecific = []
        ifoSpecific.append('--observatory ' + ifo[0])
        ifoSpecific.append('--type ' + typeNames[ifo])
        ifoSpecific.append('--output cache/' + ifo + '_' + start + '_' + end + '_frames.cache')
        command = ' '.join(dataFindCall) + ' ' + ' '.join(ifoSpecific)
        make_external_call(command)
        convertCall = cp.get('condor','convertlalcache')
        convertCall += ' '
        convertCall += 'cache/' + ifo + '_' + start + '_' + end + '_frames.cache'
        convertCall += ' '
        convertCall += 'cache/' + ifo + '_' + start + '_' + end + '_frames.wcache'
        make_external_call(convertCall)
        cp.set('omega-scans',ifo.lower() + '-omega-frame-file','../omega_setup/' + 'cache/' + ifo + '_' + start + '_' + end + '_frames.wcache')

    print("Created omega scan frame files \n")


###############################################################################
# This function will...

def create_frame_pfn_file(ifos, gpsstart, gpsend, server='internal_ldr_port'):
    namer = "pegasus-pfn-cache-"+str(gpsstart)+"-"+ \
            str(gpsend)
    gwfname = namer+".pfn" # physical file location file
    # Deletes the gwfname file if it exists prior to the execution of this
    # function.
    try:
        os.unlink(gwfname)
    except:
        pass
    # Calls ligo_data_find and passes the output to a the file gwfname.
    # The output from this will be returned and used to create the Pegasus
    # cache file in the function create_pegasus_cache_file.
    for v in ifos.keys():
        # Calls a system command to create the file.
        ldfcommand = "ligo_data_find --gps-start-time "+str(gpsstart)+ \
        " --gps-end-time "+str(gpsend)+" --observatory "+v[0]+" --type "+ ifos[v]
        if server != 'internal_ldr_port':
            ldfcommand += " --server "+ server
        ldfcommand +=" --url-type=file >> "+ gwfname
        make_external_call(ldfcommand)

    return gwfname


def create_pegasus_cache_file(framename):
    rep = open(framename,"r") # Reads the file from gwfname to be formatted properly.
    cachename = framename.split(".")[0]+".cache" # split gwfname
    # Deletes the cachename file if it already exists so a fresh file is used
    # each time.
    try:
        os.unlink(cachename)
    except:
        pass
    cac = open(cachename,"a") # Opens cachename in append mode.
    # Writes a new file from the original, unformatted create_frame_pfn_file function
    # that provides the information about the frame files in a way that the workflow
    # can read them.
    for line in open(framename):
        line = line.replace("localhost","") # remove localhost
        line2 = line.replace("\n","") # remove new line characters
        # removes the localhost/ from the beginning of the file name
        kpline = line2 # creates a better variable name for line2
        formline = kpline.split("/") # splits the line on slashes to pull out filename
        cac.write(formline[-1]+" "+kpline+ " "+ 'pool="local"' + "\n") # .cache formatter
    rep.close()
    cac.close()
    return cachename

def get_data_options(cp,ifo_name):
    if ifo_name == 'G1':
        data_opts = 'geo-data'
        try: type = cp.get('input','geo-type')
        except: type = None
        channel = cp.get('input','geo-channel')
    elif ifo_name == 'V1':
        data_opts = 'virgo-data'
        try:
            type = cp.get('input','virgo-type')
            if ('NINJA' in type) or ('T1300121_' in type):
                type = ifo_name + '_' + type
        except: type = None
        channel = cp.get('input','virgo-channel')
    else:
        data_opts = 'ligo-data'
        try:
            type = cp.get('input','ligo-type')
            if (type == 'RDS_R_L4') or ('RDS_C' in type) or ('DMT_C' in type) or ('LDAS_C' in type) or ('NINJA' in type) or ('ER_' in type) or ('T1200307_' in type) or ('T0900288_' in type) or ('ILIGO_GAUSSIAN' in type):
                type = ifo_name + '_' + type
            if ("DMT_ERHOFT" in type):
                if (ifo_name == 'L1'): type = ifo_name + '_' + type
                else: type = 'H1H2_' + type
        except: type = None
        channel = cp.get('input','ligo-channel')

    return data_opts, type, channel

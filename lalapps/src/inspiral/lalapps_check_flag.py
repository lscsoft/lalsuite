"""
check_flag.in - test a DQ flag against a set of daily ihope triggers 
                for eficiency and deadtime

$Id $

This small utility reads in a file with a csv list of veto times, and information
from a daily ihope run, and determines the efficiency and deadtime of the given
veto.  In typical usage the veto times will come from a DQ flag, prepared with
a call like

  ligolw_segment_query --database --query-segments --include-segments flag_of_interest  \
    --gps-start-time ... --gps-end-time ... | \
  ligolw_print -t segment -c start_time -c end_time > vetoed.txt

although there might also be interest in trying manually-created veto times.

Usage: lalapps_checK-flag daily_ihope_dir ifo veto_category clustering

where veto_category is 0 (science), 1, 2, or 4 and clustering is UNCLUSTERED
100MILLISEC_CLUSTERED or 16SEC_CLUSTERED

"""

import sys
import os
from glue.segments import segment, segmentlist
from optparse import OptionParser
from glue import git_version



def parse_command_line():
    """
    Parse the command line, return an options object
    """

    parser = OptionParser(
        version = "Name: %%prog\n%s" % git_version.verbose_msg,
        usage   = "%prog -v|--veto-file veto_cat [--unclustered|--thirty-ms|--sixteen-sec] [--min-snr val|--min-new-snr val] --gps-start-time start --gps-end-time end segfile1 segfile2...",
        description = "Extracts all daily ihope triggers in the specified time range, with specified clustering, above the specified (new) snr.  For each file, reports the deadtime and efficiency of the veto"
	)
    
    parser.add_option("-b", "--basedir",       metavar = "basedir",       help = "base directory for datafiles (defaults to ~cbc/ihope_daily.")
    parser.add_option("-i", "--ifo",           metavar = "ifo",           help = "ifo.")
    parser.add_option("-v", "--veto-category", metavar = "veto_category", help = "veto category to filter (0,1,2,3,4) (required).")

    parser.add_option("-u", "--unclustered",   metavar = "unclustered",           action = "store_true",    help = "Use unclustered triggers.")
    parser.add_option("-t", "--thirty-ms",     metavar = "thirty_ms_clustered",   action = "store_true",    help = "Use 30-millisec triggers.")
    parser.add_option("-x", "--sixteen-sec",   metavar = "sixteen_sec_clustered", action = "store_true",    help = "Use 16-second triggers.")

    parser.add_option("-m", "--min-snr",       metavar = "min_snr",     help = "Only examine triggers with snr above this value.")
    parser.add_option("-n", "--min-new-snr",   metavar = "min_new_snr", help = "Only examine triggers with new snr above this value.")


    # Time options
    parser.add_option("-s", "--gps-start-time", metavar = "gps_start_time", help = "Start of GPS time range")
    parser.add_option("-e", "--gps-end-time",   metavar = "gps_end_time", help = "End of GPS time range")
    
    options, others = parser.parse_args()

    if not options.ifo:
        raise ValueError, "missing required argument --ifo"

    if not options.gps_end_time:
        raise ValueError, "missing required argument --gps-end_time"
   
    if not options.gps_start_time:
        raise ValueError, "missing required argument --gps-start_time"
   
    if not options.veto_category:
        raise ValueError, "missing required argument --veto-category"
   
    if len( [x for x in (options.unclustered, options.thirty_ms, options.sixteen_sec) if x] ) != 1:
        raise ValueError, "must provide one of [--unclustered | --thirty-ms | --sixteen-sec]"

    if len( [x for x in (options.min_snr, options.min_new_snr) if x] ) != 1:
        raise ValueError, "must provide exactly one of [--min-snr | --min-new-snr]"

    if len(others) == 0:
        raise ValueError, "must provide at least one file of segments"

    return options, others



def get_summary(basedir, ifo, cluster, cat, start_time, end_time):
    all_sum     = segmentlist([])
    cur_time    = start_time

    while cur_time < end_time:
        tstring  = os.popen('tconvert -f %Y%m/%Y%m%d ' + str(cur_time)).readlines()[0].strip()
        infile   = open('%s/%s/%s-0-SUMMARY_%s.csv'  % (basedir, tstring, ifo, cluster))
        lines    = [l.strip().split(',') for l in infile.readlines()]
        summary  = segmentlist([segment(int(l[0]), int(l[1])) for l in lines]).coalesce()
        all_sum  = all_sum + summary

        cur_time += 60 * 60 * 24

    all_sum = all_sum & segmentlist([segment(start_time, end_time)])

    return all_sum



def get_triggers(basedir, ifo, cluster, cat, start_time, end_time, snr_func, filter_func):
    all_triggers = []
    cur_time     = start_time

    while cur_time < end_time:
        tstring  = os.popen('tconvert -f %Y%m/%Y%m%d ' + str(cur_time)).readlines()[0].strip()
        infile   = open('%s/%s/%s-%s-INSPIRAL_%s.csv'  % (basedir, tstring, ifo, cat, cluster))
        data     = [ t.split(',') for t in infile ]
        triggers = [ (int(trigger[0]), snr_func(trigger)) for trigger in data ]
        triggers = [ t for t in triggers if filter_func(t) ]

        all_triggers.extend(triggers)

        cur_time += 60 * 60 * 24

    return all_triggers


def get_snr(t):
    return float(t[3])

def get_new_snr(d, index=6.0):
    snr       = float(d[3])
    chisq     = float(d[11])
    chisq_dof = float(d[12])

    rchisq = chisq/(2 * chisq_dof - 2)
    nhigh  = 2.

    if rchisq > 1.:
        return snr / ((1. + rchisq**(index/nhigh))/2)**(1./index)
    else:
        return snr

if __name__ == '__main__':
    options, others = parse_command_line()
    ifo   = options.ifo
    start = int(options.gps_start_time)
    end   = int(options.gps_end_time)
    cat   = int(options.veto_category)

    basedir  = options.basedir or "/archive/home/cbc/ihope_daily"
    snrifier = (options.min_snr and get_snr) or get_new_snr
    snr_lim  = float(options.min_snr or options.min_new_snr)
    cluster  = (options.unclustered and "UNCLUSTERED") or (options.thirty_ms and "30MILLISEC_CLUSTERED") or (options.sixteen_sec and "16SEC_CLUSTERED")
    

    filter_func = lambda x: x[0] >= start and x[0] < end and x[1] >= snr_lim

    trigs = get_triggers(basedir, ifo, cluster, cat, start, end, snrifier, filter_func)
    summ  = get_summary(basedir, ifo, cluster, cat, start, end)

    incount = len(trigs)

    if incount == 0:
        print "No triggers found"
        sys.exit(0)

    for filename in others:
        infile  = open(filename)
        lines   = [l.strip().split(',') for l in infile.readlines()]
        vetoed  = segmentlist([segment(int(l[0]), int(l[1])) for l in lines]).coalesce()

        new_summary = summ - vetoed
        new_trigs   = [t for t in trigs if t[0] not in vetoed]
        new_trigs   = sorted(new_trigs, cmp=lambda x,y: cmp(y[1],x[1]))
       
        outcount    = len(new_trigs) 
        efficiency  = (float(incount) - float(outcount)) / float(incount) * 100.0
        deadtime    = float(abs(summ) - abs(new_summary))  / float(abs(summ)) * 100.0

        print "File: ", filename
        print "Efficiency: %.2f" % efficiency
        print "Deadtime: %.2f" % deadtime
        print "Ratio: %s" % (deadtime > 0 and "%.2f" % (efficiency / deadtime) or 'NA')
        print "Loudest remaining trigger at %d with snr %.2f" % (new_trigs[0][0], new_trigs[0][1])
        print


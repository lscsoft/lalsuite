# DAG generation code for running LALInference pipeline
# (C) 2012 John Veitch, Vivien Raymond

from lalinference import lalinference_pipe_utils as pipe_utils
from lalapps import inspiralutils
from six.moves import configparser
from optparse import OptionParser,OptionValueError
import sys
import ast
import os
import uuid
from glue import pipeline
from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import utils as ligolw_utils 
from math import ceil

usage=""" %prog [options] config.ini
Setup a Condor DAG file to run the LALInference pipeline based on
the config.ini file.

The user can specify either an injection file to analyse, with the --inj option,
a list of SnglInspiralTable or CoincInspiralTable triggers with the --<x>-triggers options,
a GraceDB ID with the --gid option,
or an ASCII list of GPS times with the --gps-time-file option.

If none of the above options are given, the pipeline will analyse the entire
stretch of time between gps-start-time and gps-end-time, specified in the config.ini file.

The user must also specify and ini file which will contain the main analysis config.

"""
parser=OptionParser(usage)
parser.add_option("-r","--run-path",default=None,action="store",type="string",help="Directory to run pipeline in (default: $PWD)",metavar="RUNDIR")
parser.add_option("-p","--daglog-path",default=None,action="store",type="string",help="Path to directory to contain DAG log file. SHOULD BE LOCAL TO SUBMIT NODE",metavar="LOGDIR")
parser.add_option("-g","--gps-time-file",action="store",type="string",default=None,help="Text file containing list of GPS times to analyse",metavar="TIMES.txt")
parser.add_option("--gid",action="store",type="string",default=None,help="GraceDB ID")
parser.add_option("--service-url",action="store",type="string",default=None,help="GraceDB url from which xml files are downloaded")
parser.add_option("-I","--injections",action="store",type="string",default=None,help="List of injections to perform and analyse",metavar="INJFILE.xml")
parser.add_option("-B","--burst_injections",action="store",type="string",default=None,help="SimBurst table for LIB injections",metavar="INJFILE.xml")
parser.add_option("-P","--pipedown-db",action="store",type="string",default=None,help="Pipedown database to read and analyse",metavar="pipedown.sqlite")
parser.add_option("--condor-submit",action="store_true",default=False,help="Automatically submit the condor dag")

(opts,args)=parser.parse_args()

# We will store variations of runs to use here as a dictionary of
# ((section, option), ...) : [variation1, ...]
# Where (section, option) indicate the configparser option to vary
# and multiple section,option pairs are allowed
variations={}

def mkdirs(path):
    """
    Helper function. Make the given directory, creating intermediate
    dirs if necessary, and don't complain about it already existing.
    """
    if os.access(path,os.W_OK) and os.path.isdir(path): return
    else: os.makedirs(path)

def add_variations(cp, section, option, values=None, allowed_values=None):
    """
    Push some possible variations onto the stack.
    If only one value is specified then just store it in cp as usual
    cp : ConfigParser object to look in
    section : [section] in cp
    option : option in section
    values: If given, use instead of cp's
    allowed_values : if given, anything else will trigger an error
    """
    if not cp.has_section(section) and not cp.has_option(section,option):
        return

    if values is not None:
        vals = values
    else:
        vals = cp.get(section,option).split(',')
    if allowed_values is not None and any(not vals in allowed_values):
        raise ValueError("Unknown value for {section} . {option} {value}"
                         .format(section, option, \
                             ' '.join([v for v in vals if v not in allowed_options])
                                )
                         )
    if len(vals) >1:
        return {(section,option): vals}
    elif len(vals)==1:
        cp.set(section, option, vals[0])
        return {}
    else:
        print(("Found no variations of [{section}] {option}".format(section=section,
                                                                   option=option)))
        return {}

def generate_variations(master_cp, variations):
    """
    Generate config parser objects for each of the variations
    """
    cur_basedir=master_cp.get('paths','basedir')
    mkdirs(cur_basedir)
    # save the master file
    masterpath=os.path.join(cur_basedir,'config.ini')
    with open(masterpath,'w') as cpfile:
        master_cp.write(cpfile)

    cur_webdir=master_cp.get('paths','webdir')
    cur_dagdir=master_cp.get('paths','daglogdir')
    # If no variations, done
    if not variations:
        yield master_cp
        raise StopIteration()

    # Otherwise, vary over the next option
    (section, opt), vals = variations.popitem()
    for val in vals:
        # Read file back in to get a new object
        cp = configparser.ConfigParser()
        cp.optionxform = str
        cp.read(masterpath)
        cp.set(section,opt,val)
        # Append to the paths
        cp.set('paths','basedir',os.path.join(cur_basedir, \
                                        '{val}'.format(opt=opt,val=val)))
        cp.set('paths','webdir',os.path.join(cur_webdir, \
                                        '{val}'.format(opt=opt,val=val)))
        cp.set('paths','daglogdir',os.path.join(cur_dagdir, \
                                        '{val}'.format(opt=opt,val=val)))
        # recurse into remaining options
        for sub_cp in generate_variations(cp,variations):
            yield sub_cp

if len(args)!=1:
    parser.print_help()
    print('Error: must specify one ini file')
    sys.exit(1)

inifile=args[0]

cp=configparser.SafeConfigParser()
fp=open(inifile)
cp.optionxform = str
cp.readfp(fp)
fp.close()

# Set the base directory for the run
if opts.run_path is not None:
    cp.set('paths','basedir',os.path.abspath(opts.run_path))
if cp.get('paths','basedir') is None:
    print("Warning: no run dir set, using current dir")
    cp.set('paths','basedir',os.path.getcwd())
mkdirs(cp.get('paths','basedir'))

if opts.daglog_path is not None:
    cp.set('paths','daglogdir',os.path.abspath(opts.daglog_path))
else:
    cp.set('paths','daglogdir',os.path.abspath(cp.get('paths','basedir')))
daglogdir=cp.get('paths','daglogdir')
mkdirs(daglogdir)


# Set up from all the various input options
if opts.gps_time_file is not None:
    cp.set('input','gps-time-file',os.path.abspath(opts.gps_time_file))

if opts.injections is not None:
    cp.set('input','injection-file',os.path.abspath(opts.injections))

if opts.burst_injections is not None:
    if opts.injections is not None:
        print("ERROR: cannot pass both inspiral and burst tables for injection\n")
        sys.exit(1)
    cp.set('input','burst-injection-file',os.path.abspath(opts.burst_injections))

if opts.gid is not None:
    cp.set('input','gid',opts.gid)

if opts.service_url is not None:
    cp.set('analysis','service-url',opts.service_url)

if opts.pipedown_db is not None:
    cp.set('input','pipedown-db',os.path.abspath(opts.pipedown_db))

# Some sanity checking
approx='approx'
if not (cp.has_option('engine','approx') or cp.has_option('engine','approximant') ):
    print("Error: was expecting an 'approx' filed in the [engine] section\n")
    sys.exit(1)

# Build a list of allowed variations
variations.update(add_variations(cp, 'engine','approx'))
variations.update(add_variations(cp, 'analysis', 'engine'))
variations.update(add_variations(cp, 'analysis', 'roq'))

roq_paths=[]
def setup_roq(cp):
    """
    Generates cp objects with the different ROQs applied
    """
    use_roq=False
    if cp.has_option('paths','roq_b_matrix_directory') or cp.has_option('paths','computeroqweights'):
        if not cp.has_option('analysis','roq'):
            print("Warning: If you are attempting to enable ROQ by specifying roq_b_matrix_directory or computeroqweights,\n\
            please use analysis.roq in your config file in future. Enabling ROQ.")
            cp.set('analysis','roq',True)
    if not cp.getboolean('analysis','roq'):
        yield cp
        raise StopIteration()
    from numpy import genfromtxt, array
    path=cp.get('paths','roq_b_matrix_directory')
    if not os.path.isdir(path):
        print("The ROQ directory %s does not seem to exist\n"%path)
        sys.exit(1)
    use_roq=True
    roq_paths=os.listdir(path)
    roq_params={}
    roq_force_flow = None

    if cp.has_option('lalinference','roq_force_flow'):
        roq_force_flow = cp.getfloat('lalinference','roq_force_flow')
        print("WARNING: Forcing the f_low to ", str(roq_force_flow), "Hz")
        print("WARNING: Overwriting user choice of flow, srate, seglen, and (mc_min, mc_max and q-min) or (mass1_min, mass1_max, mass2_min, mass2_max)")

    def key(item): # to order the ROQ bases
        return float(item[1]['seglen'])

    coinc_xml_obj = None
    row=None

    # Get file object of coinc.xml
    if opts.gid is not None:
        from ligo.gracedb.rest import GraceDb
        gid=opts.gid
        cwd=os.getcwd()
        if cp.has_option('analysis', 'service-url'):
            client = GraceDb(cp.get('analysis', 'service-url'))
        else:
            client = GraceDb()
        coinc_xml_obj = ligolw_utils.load_fileobj(client.files(gid, "coinc.xml"), contenthandler = lsctables.use_in(ligolw.LIGOLWContentHandler))[0]
    elif cp.has_option('input', 'coinc-xml'):
        coinc_xml_obj = ligolw_utils.load_fileobj(open(cp.get('input', 'coinc-xml'), "rb"), contenthandler = lsctables.use_in(ligolw.LIGOLWContentHandler))[0]

    # Get sim_inspiral from injection file
    if cp.has_option('input','injection-file'):
        print("Only 0-th event in the XML table will be considered while running with ROQ\n")
        row = lsctables.SimInspiralTable.get_table(
                  ligolw_utils.load_filename(cp.get('input','injection-file'),contenthandler=lsctables.use_in(ligolw.LIGOLWContentHandler))
              )[0]

    roq_bounds = pipe_utils.Query_ROQ_Bounds_Type(path, roq_paths)
    if roq_bounds == 'chirp_mass_q':
        print('ROQ has bounds in chirp mass and mass-ratio')
        mc_priors, trigger_mchirp = pipe_utils.get_roq_mchirp_priors(
            path, roq_paths, roq_params, key, coinc_xml_obj=coinc_xml_obj, sim_inspiral=row
        )
    elif roq_bounds == 'component_mass':
        print('ROQ has bounds in component masses')
        # get component mass bounds, then compute the chirp mass that can be safely covered
        # further below we pass along the component mass bounds to the sampler, not the tighter chirp-mass, q bounds
        m1_priors, m2_priors, trigger_mchirp = pipe_utils.get_roq_component_mass_priors(
            path, roq_paths, roq_params, key, coinc_xml_obj=coic_xml_obj, sim_inspiral=row
        )
        mc_priors = {}
        for (roq,m1_prior), (roq2,m2_prior) in zip(m1_priors.items(), m2_priors.items()):
            mc_priors[roq] = sorted([pipe_utils.mchirp_from_components(m1_prior[1], m2_prior[0]), pipe_utils.mchirp_from_components(m1_prior[0], m2_prior[1])])

    if cp.has_option('lalinference','trigger_mchirp'):
        trigger_mchirp=float(cp.get('lalinference','trigger_mchirp'))
    roq_mass_freq_scale_factor = pipe_utils.get_roq_mass_freq_scale_factor(mc_priors, trigger_mchirp, roq_force_flow)
    if roq_mass_freq_scale_factor != 1.:
        print('WARNING: Rescaling ROQ basis, please ensure it is allowed with the model used.')

    # If the true chirp mass is unknown, add variations over the mass bins
    if opts.gid is not None or (opts.injections is not None or cp.has_option('input','injection-file')) or cp.has_option('lalinference','trigger_mchirp') or cp.has_option('input', 'coinc-xml'):

        for mc_prior in mc_priors:
            mc_priors[mc_prior] = array(mc_priors[mc_prior])
        # find mass bin containing the trigger
        trigger_bin = None
        for roq in roq_paths:
            if mc_priors[roq][0]*roq_mass_freq_scale_factor <= trigger_mchirp <= mc_priors[roq][1]*roq_mass_freq_scale_factor:
                trigger_bin = roq
                print('Prior in Mchirp will be ['+str(mc_priors[roq][0]*roq_mass_freq_scale_factor)+','+str(mc_priors[roq][1]*roq_mass_freq_scale_factor)+'] to contain the trigger Mchirp '+str(trigger_mchirp))
                break
        roq_paths = [trigger_bin]
    else:
        for mc_prior in mc_priors:
            mc_priors[mc_prior] = array(mc_priors[mc_prior])*roq_mass_freq_scale_factor

    # write the master configparser
    cur_basedir = cp.get('paths','basedir')
    masterpath=os.path.join(cur_basedir,'config.ini')
    with open(masterpath,'w') as cpfile:
        cp.write(cpfile)

    for roq in roq_paths:
        this_cp = configparser.ConfigParser()
        this_cp.optionxform = str
        this_cp.read(masterpath)
        basedir = this_cp.get('paths','basedir')
        for dirs in 'basedir','daglogdir','webdir':
            val = this_cp.get('paths',dirs)
            newval = os.path.join(val,roq)
            mkdirs(newval)
            this_cp.set('paths',dirs,newval)
        this_cp.set('paths','roq_b_matrix_directory',os.path.join(cp.get('paths','roq_b_matrix_directory'),roq))
        flow=roq_params[roq]['flow'] / roq_mass_freq_scale_factor
        srate=2.*roq_params[roq]['fhigh'] / roq_mass_freq_scale_factor
        if srate > 8192:
            srate = 8192

        seglen=roq_params[roq]['seglen'] * roq_mass_freq_scale_factor
        # params.dat uses the convention q>1 so our q_min is the inverse of their qmax
        this_cp.set('engine','srate',str(srate))
        this_cp.set('engine','seglen',str(seglen))
        if this_cp.has_option('lalinference','flow'):
            tmp=this_cp.get('lalinference','flow')
            tmp=eval(tmp)
            ifos=tmp.keys()
        else:
            tmp={}
            ifos=eval(this_cp.get('analysis','ifos'))
            for i in ifos:
                tmp[i]=flow
                this_cp.set('lalinference','flow',str(tmp))
        if roq_bounds == 'chirp_mass_q':
            mc_min=mc_priors[roq][0]*roq_mass_freq_scale_factor
            mc_max=mc_priors[roq][1]*roq_mass_freq_scale_factor
            # params.dat uses the convention q>1 so our q_min is the inverse of their qmax
            q_min=1./float(roq_params[roq]['qmax'])
            this_cp.set('engine','chirpmass-min',str(mc_min))
            this_cp.set('engine','chirpmass-max',str(mc_max))
            this_cp.set('engine','q-min',str(q_min))
            this_cp.set('engine','comp-min', str(max(roq_params[roq]['compmin'] * roq_mass_freq_scale_factor, mc_min * pow(1+q_min, 1./5.) * pow(q_min, 2./5.))))
            this_cp.set('engine','comp-max', str(mc_max * pow(1+q_min, 1./5.) * pow(q_min, -3./5.)))
        elif roq_bounds == 'component_mass':
            m1_min = m1_priors[roq][0]
            m1_max = m1_priors[roq][1]
            m2_min = m2_priors[roq][0]
            m2_max = m2_priors[roq][1]
            this_cp.set('engine','mass1-min',str(m1_min))
            this_cp.set('engine','mass1-max',str(m1_max))
            this_cp.set('engine','mass2-min',str(m2_min))
            this_cp.set('engine','mass2-max',str(m2_max))
        yield this_cp
    raise StopIteration()

# Create an outer dag to wrap the sub-dags
outerdaglog=os.path.join(daglogdir,'lalinference_multi_'+str(uuid.uuid1())+'.log')
outerdag=pipeline.CondorDAG(outerdaglog,dax=False)
outerdag.set_dag_file(os.path.join(cp.get('paths','basedir'),'multidag'))


master_cp=cp
# Iterate over variations and generate sub-dags
for cp in generate_variations(master_cp,variations):
    basepath=cp.get('paths','basedir')
    # Link injection file into place as paths outside basedir are inaccessible to containerised jobs
    if cp.has_option('input','injection-file'):
        injpath=cp.get('input','injection-file')
        myinjpath=os.path.join(basepath,os.path.basename(injpath))
        if os.path.abspath(myinjpath) != os.path.abspath(injpath):
            # If the injection file does not exist in the run dir, link it into place
            # Useful for singularity jobs which see only rundir
            if os.path.lexists(myinjpath):
                # If the path exists, see if it is a link to the current file
                # and if so, just update the config
                if os.path.islink(myinjpath) and os.path.realpath(myinjpath)==os.path.realpath(injpath):
                    cp.set('input','injection-file',myinjpath)
                else:
                    # Do not over-write the injection file
                    print(("Error: File {0} exists in run directory, not over-writing with \
                            {1}. Remove the existing file or create a fresh run directory".format(myinjpath,injpath)
                            ))
                    sys.exit(1)
            else:
                # The link doens't exist, so create it and update config
                try:
                    os.link(os.path.abspath(injpath), myinjpath)
                except:
                    from shutil import copyfile
                    copyfile(injpath,myinjpath)
                cp.set('input','injection-file',myinjpath)

    for this_cp in setup_roq(cp):
        # Create the DAG from the configparser object
        dag=pipe_utils.LALInferencePipelineDAG(this_cp,dax=False)
        dagjob=pipeline.CondorDAGManJob(os.path.join(this_cp.get('paths','basedir'),dag.get_dag_file()),
                                    this_cp.get('paths','basedir'))
        dagnode=pipeline.CondorDAGManNode(dagjob)
        outerdag.add_node(dagnode)
        dag.write_sub_files()
        dag.write_dag()
        dag.write_script()

outerdag.write_sub_files()
outerdag.write_dag()
outerdag.write_script()

# End of program
print('Successfully created DAG file.')

if opts.condor_submit:
    import subprocess
    from subprocess import Popen
    if cp.has_option('condor','notification'):
        x = subprocess.Popen(['condor_submit_dag','-dont_suppress_notification',outerdag.get_dag_file()])
    else:
        x = subprocess.Popen(['condor_submit_dag',outerdag.get_dag_file()])
    x.wait()
    if x.returncode==0:
        print('Submitted DAG file')
    else:
        print('Unable to submit DAG file')
else:
    print('To submit, run:\n\tcondor_submit_dag {0}'.format(outerdag.get_dag_file()))

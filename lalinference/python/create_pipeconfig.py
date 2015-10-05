#!/usr/bin/env python

"""
This script creates an injection and a .ini file to run LALInference

(C) Archisman Ghosh, Abhirup Ghosh, 2015-09-20
"""

import sys, lal, os, commands, numpy as np
import imrtestgr as tgr
import nr_fits as nr


os.system('. ${HOME}/.profile')
home = os.getenv('HOME')
user = commands.getoutput('whoami')
psd_path = os.path.join(home, 'src/lalsuite/lalsimulation/src')
lal_prefix = os.getenv('LAL_PREFIX')
glue_location = os.getenv('GLUE_LOCATION')
pylal_location = os.getenv('PYLAL_LOCATION')

# ### INJECTION ###

# Waveform info
#approximant = 'SEOBNRv2_ROM_DoubleSpinthreePointFivePN' # injection waveform 
approximant = 'IMRPhenomBpseudoFourPN' # injection waveform 
amp_order = -1 # amplitude PN order of the injection waveform
f_low = 25. # low-frequency cutoff (Hz)
waveform_info = '--waveform %s --amp-order %d --f-lower %f --taper-injection startend'%(approximant, amp_order, f_low)

# PSD info ### FIXME ###
ligo_psd = 'LALAdLIGO' #os.path.join(psd_path, 'LIGO-P1200087-v18-aLIGO_EARLY_HIGH.txt')
#virgo_psd = os.path.join(psd_path, 'LIGO-P1200087-v18-AdV_EARLY_HIGH.txt')
f_start = 30.
psd_info = '--ifos H1,L1 --ligo-fake-psd "%s" --ligo-start-freq %f'%(ligo_psd, f_start)

# Time info
gps_time = 1126285216 # O1 start time 
time_step = 2630./np.pi # time interval between nearby injections 
time_info = '--t-distr fixed --gps-start-time %f --gps-end-time %f --time-step %f'%(gps_time, gps_time, time_step)

# Parameters
M_list = [50., 75., 100., 150., 200.] # total mass (M_sun)
q_list = [1., 2., 4.]
inj_snr = 50
iota = 0. # inclination (degrees)
phi_c = 0. # coalescent phase (degrees)
psi = 0. # polarization (degrees)
alpha = 0. # right ascension (degrees)
delta = 0. # declination (degrees)

for M in M_list:
  for q in q_list:
    m1 = M/(1.+q)
    m2 = M*q/(1.+q)
    print 'm1 = %.1f,\tm2 = %.1f'%(m1, m2)
    param_info = '--m-distr fixMasses --fixed-mass1 %f --fixed-mass2 %f --snr-distr uniform --min-snr %f --max-snr %f --i-distr fixed --fixed-inc %f --coa-phase-distr fixed --fixed-coa-phase %f --polarization %f --l-distr fixed --latitude %f --longitude %f --disable-spin'%(m1, m2, inj_snr, inj_snr, iota, phi_c, psi, alpha, delta)

    # Output info
    inj_file = '../injections/example_inj_%s_snr%d_%s_%s'%(approximant, inj_snr,str(M),str(q)) # output file name for xml and txt file
    output_info = '--output %s.xml'%(inj_file)

    # Generate the injection xml file by calling lalapps_inspinj 
    run_command = 'lalapps_inspinj %s %s %s %s %s'%(waveform_info, psd_info, time_info, param_info, output_info)
    print(run_command)
    os.system(run_command)

    # Print some relevant columns of the xml file to an ASCII table with header 
    os.system('ligolw_print -d " " -c mass1 -c mass2 -c longitude -c latitude -c inclination -c polarization -c distance -c geocent_end_time -c geocent_end_time_ns -c coa_phase %s.xml > %s.tmp' %(inj_file, inj_file))
    os.system('echo \# mass1 mass2 longitude latitude inclination polarization distance geocent_end_time geocent_end_time_ns coa_phase > header.txt')
    os.system('cat header.txt %s.tmp > %s.txt' %(inj_file, inj_file))
    os.system('rm %s.tmp header.txt' %(inj_file))

# ### RECOVERY ###

ligo_psd = os.path.join(psd_path, 'LIGO-P1200087-v18-aLIGO_EARLY_HIGH.txt')
virgo_psd = os.path.join(psd_path, 'LIGO-P1200087-v18-AdV_EARLY_HIGH.txt')


def write_pipeline(out_folder, M, q, flow=30., fhigh=1023., nlive=1024, create_dag=True):
    os.system('mkdir -p %s'%(out_folder))

    inj_file = '../injections/example_inj_%s_snr%d_%s_%s'%(approximant, inj_snr,str(M),str(q)) # output file name for xml and txt file
    os.system('cp %s.xml %s/injection.xml'%(inj_file, out_folder))

    m1 = M/(1.+q)
    m2 = M*q/(1.+q)    

    ofile = open(os.path.join(out_folder, 'pipeconfig.ini'), 'w')
    
    ofile.write("[analysis]\n")
    ofile.write("ifos=['H1','L1']\n")
    ofile.write("engine=lalinferencenest\n")
    ofile.write("nparallel=4\n")
    ofile.write("coherence-test=False\n")
    ofile.write("upload-to-gracedb=False\n")
    ofile.write("dataseed=1\n")
    ofile.write("\n")
    
    projectdir = out_folder.replace('%s/Documents/Work/'%home, '')

    ofile.write("[paths]\n")
    ofile.write("webdir=%s/public_html/%s\n"%(home, projectdir))
    ofile.write("basedir=%s/public_html/%s\n"%(home, projectdir))
    ofile.write("baseurl=https://dogmatix.icts.res.in/~%s/%s\n"%(user, projectdir))
    ofile.write("\n")
    
    ofile.write("[input]\n")
    ofile.write("max-psd-length=1024\n")
    ofile.write("padding=16\n")
    ofile.write("timeslides=false\n")
    ofile.write("ignore-science-segments=True\n")
    ofile.write("events=all\n")
    ofile.write("\n")
    
    ofile.write("[datafind]\n")
    ofile.write("types={'H1':'LALSimAdLIGO','L1':'LALSimAdLIGO','V1':'LALSimAdVirgo','I1':'LALSimAdLIGO','J1':'LALSimAdLIGO'}\n")
    ofile.write("\n")
    
    ofile.write("[data]\n")
    ofile.write("channels={'H1':'LALSimAdLIGO','L1':'LALSimAdLIGO','V1':'LALSimAdVirgo','I1':'LALSimAdLIGO','J1':'LALSimAdLIGO'}\n")
    ofile.write("\n")
    
    ofile.write("[condor]\n")
    ofile.write("mpirun=/bin/true\n")
    ofile.write("datafind=/bin/true\n")
    ofile.write("gracedb=/bin/true\n")
    ofile.write("ligolw_print=%s/bin/ligolw_print\n"%(glue_location))
    ofile.write("lalinferencenest=%s/bin/lalinference_nest\n"%(lal_prefix))
    ofile.write("lalinferencemcmc=%s/bin/lalinference_mcmc\n"%(lal_prefix))
    ofile.write("mergescript=%s/bin/lalapps_nest2pos\n"%(lal_prefix))
    ofile.write("resultspage=%s/bin/cbcBayesPostProc.py\n"%(pylal_location))
    ofile.write("coherencetest=%s/bin/lalapps_coherence_test\n"%(lal_prefix))
    ofile.write("\n")
    
    ofile.write("[engine]\n")
    ofile.write("approx=%s\n"%(approximant))
    ofile.write("seglen=8\n")
    ofile.write("nlive=%d\n"%(nlive))
    ofile.write("srate=2048\n")
    #ofile.write("dt=0.1\n")
    #ofile.write("distance-min=1\n")
    #ofile.write("distance-max=50000\n")
    ofile.write("comp-min=5.0\n") # ### FIXME ###
    ofile.write("comp-max=%.2f\n"%(8.*m2)) # ### FIXME ###
    ofile.write("disable-spin=\n") # ### FIXME ###
    ofile.write("amporder=-1\n")
    ofile.write("margphi=\n")
    #ofile.write("marginal-d=\n")
    #ofile.write("use-logdistance=\n")
    ofile.write("resume=\n")
    ofile.write("progress=\n")
    ofile.write("H1-psd=%s\n"%(ligo_psd))
    ofile.write("L1-psd=%s\n"%(ligo_psd))
    #ofile.write("0noise=\n") # ### FIXME ###
    ofile.write("\n")
    
    ofile.write("[lalinference]\n")
    #ofile.write("fake-cache={'H1':'LALSimAdLIGO','L1':'LALSimAdLIGO','V1':'LALSimAdVirgo','I1':'LALSimAdLIGO','J1':'LALSimAdLIGO'}\n")
    ofile.write("fake-cache={'H1':'interp:%s','L1':'interp:%s','V1':'interp:%s'}\n"%(ligo_psd, ligo_psd, virgo_psd))
    ofile.write("flow={'H1':%.2f,'L1':%.2f,'V1':%.2f}\n"%(flow, flow, flow))
    ofile.write("fhigh={'H1':%.2f,'L1':%.2f,'V1':%.2f}\n"%(fhigh, fhigh, fhigh))
    ofile.write("\n")
    
    ofile.write("[merge]\n")
    #ofile.write("npos=5000\n")
    ofile.write("\n")
    
    ofile.write("[resultspage]\n")
    ofile.write("skyres=0.5\n")
    ofile.write("\n")
    
    ofile.close()
    
    if create_dag:
      os.system('%s/bin/lalinference_pipe -I %s/injection.xml -r %s -p %s %s/pipeconfig.ini'%(lal_prefix, out_folder, out_folder, out_folder, out_folder))

# Main
for M in M_list:
  for q in q_list:
    # ### IMR ###
    noise_type = 'LIGO-P1200087-v18-aLIGO_EARLY_HIGH'
    date = '2015-09-28'
    run_tag = 'nospin_Pan_etal_2011_seglen8'

    imr_folder = '%s/Documents/Work/testGR_IR/runs/simulations/gaussian_noise/%s/%s/%s/%s/IMR/%s_%s'%(home, noise_type, approximant, date, run_tag, M, q)
    insp_folder = '%s/Documents/Work/testGR_IR/runs/simulations/gaussian_noise/%s/%s/%s/%s/inspiral/%s_%s'%(home, noise_type, approximant, date, run_tag, M, q)
    ring_folder = '%s/Documents/Work/testGR_IR/runs/simulations/gaussian_noise/%s/%s/%s/%s/ringdown/%s_%s'%(home, noise_type, approximant, date, run_tag, M, q)

    m1 = M/(1.+q)
    m2 = M*q/(1.+q)
    chi1 = -0.0
    chi2 = 0.    
    
    fit_formula = 'nospin_Pan2011'
    #fit_formula = 'nonprecspin_Healy2014'

    # calculate the mass and spin of the final BH 
    Mf, af = tgr.calc_final_mass_spin(m1, m2, chi1, chi2, fit_formula)

    # calculate the Kerr ISCO freq 
    f_isco_Kerr = nr.calc_isco_freq(af)/(Mf*lal.MTSUN_SI)

    # calculate the dominant QNM freq 
    f_qnm = nr.calc_fqnm_dominant_mode(af)/(Mf*lal.MTSUN_SI)
 

    insp_fhigh = f_isco_Kerr
    ring_flow = f_isco_Kerr
    
    write_pipeline(imr_folder, M, q)
    write_pipeline(insp_folder, M, q, fhigh=insp_fhigh)
    write_pipeline(ring_folder, M, q, flow=ring_flow)

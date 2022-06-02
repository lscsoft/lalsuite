import numpy as np
import subprocess as sp

# from LALConstants.h
LAL_C_SI = 299792458e0
LAL_REARTH_SI = 6378136.6
LAL_AU_SI = 149597870700e0

# compute times for a variety of telescopes and sky positions
telescopes = (
    'H1', 'L1', 'V1', 'G1', 'T1',
    'GBT', 'PKS', 'JBO', 'AO', 'EFF', 'NRT', 'HOB', 'HART', 'VLA', 'WSRT'
)
gps_ssb = 1100000000.0
gps_dets = []
for telescope in telescopes:
    for ra in np.arange(0, 24, 8):
        cmd = f'lalapps_ssbtodetector --gps {gps_ssb} --ra {ra}:21:34.76 --dec -1:53:12.36 -t {telescope}'
        out = sp.check_output(cmd, shell=True)
        gps_dets.append(float(out))
gps_dets = np.array(gps_dets)

# check difference between SSB and detector times are sensible

dt_dets_ssb = abs(gps_dets - gps_ssb)
max_dt_dets_ssb = (LAL_AU_SI + LAL_REARTH_SI) / LAL_C_SI
assert max(dt_dets_ssb) <= max_dt_dets_ssb, f'{max(dt_dets_ssb)} <= {max_dt_dets_ssb}'

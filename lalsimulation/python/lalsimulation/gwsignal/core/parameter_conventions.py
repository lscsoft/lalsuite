import numpy as np
from astropy import units as u

# Define lists of allowed mass / spin parameter names
mass_params_old_ = ["mass1_det", "mass2_det",
                "total_mass_det", "chirp_mass_det",
                "mass_difference_det", "reduced_mass_det"]


mass_params_ = ["mass1", "mass2","total_mass", "chirp_mass",
                "mass_difference", "reduced_mass"]

mass_params = mass_params_+["mass_ratio", "sym_mass_ratio"]

spin_params = ["spin1x", "spin1y", "spin1z", "spin2x", "spin2y", "spin2z",
            "spin1_norm","spin1_tilt","spin1_phi","spin2_norm","spin2_tilt","spin2_phi"]

# Waveform generation parameters
gen_params = ["deltaT", "deltaF", "f22_start", "f_max", "phi_ref", "f22_ref","lmax"]

#Extrinsic parameters
extrinsic_params = ["distance", "inclination",  "longAscNodes", "meanPerAno"]

# Condition Parameters
condition_params = ["condition"]

#Tidal parameters
tidal_params = ["lambda1","lambda2","TidalOctupolarLambda1","TidalOctupolarLambda2",
                "TidalHexadecapolarLambda1","TidalHexadecapolarLambda2",
                "TidalQuadrupolarFMode1","TidalQuadrupolarFMode2",
                "TidalOctupolarFMode1","TidalOctupolarFMode2"]

#Non GR parameters
nongr_params = ["phi1","phi2","phi3","phi4","dchi0","dchi1","dchi2","dchi3","dchi4",
                "dchi5","dchi5l","dchi6","dchi6l","dchi7","dxi1","dxi2","dxi3","dxi4",
                "dxi5","dxi6","dsigma1","dsigma2","dsigma3","dsigma4","dalpha1",
                "dalpha2","dalpha3","dalpha4","dalpha5","dbeta1","dbeta2","dbeta3",
                "alphaPPE","betaPPE","alphaPPE0","betaPPE0","alphaPPE1","betaPPE1",
                "alphaPPE2","betaPPE2","alphaPPE3","betaPPE3","alphaPPE4","betaPPE4",
                "alphaPPE5","betaPPE5","alphaPPE6","betaPPE6","alphaPPE7","betaPPE7",
                "liv","log10lambda_eff","LIV_A_sign","nonGR_alpha","dchikappaS",
                "dchikappaA","domega220","dtau220","domega210","dtau210","domega330",
                "dtau330","domega440","dtau440","domega550","dtau550",]

#Other parameters

other_params = ["eccentricity"]

_other_params = ["modes","axis","NumRelData", "approximant"]

#Array-like parameters

arr_params = ["ModeArray", "ModeArrayJframe"]


full_parameter_list = np.concatenate([mass_params, spin_params, gen_params,extrinsic_params,
                                      other_params,_other_params,arr_params,tidal_params,nongr_params, condition_params])

#Units shared between dictionaries
common_units_dictionary = {
                        "spin1x" : u.dimensionless_unscaled,
                        "spin1y" : u.dimensionless_unscaled,
                        "spin1z" : u.dimensionless_unscaled,
                        "spin2x" : u.dimensionless_unscaled,
                        "spin2y" : u.dimensionless_unscaled,
                        "spin2z" : u.dimensionless_unscaled,
                        "spin1_norm":u.dimensionless_unscaled,
                        "spin2_norm":u.dimensionless_unscaled,
                        "spin1_tilt":u.rad,
                        "spin2_tilt":u.rad,
                        "spin1_phi":u.rad,
                        "spin2_phi":u.rad,
                        "deltaT":u.s,
                        "deltaF":u.Hz,
                        "f_min":u.Hz,
                        "f22_start":u.Hz,
                        "f22_ref":u.Hz,
                        "f_max":u.Hz,
                        "f_ref":u.Hz,
                        "phi_ref":u.rad,
                        "inclination":u.rad,
                        "eccentricity":u.dimensionless_unscaled,
                        "longAscNodes":u.rad,
                        "meanPerAno":u.rad
                        }

common_units_dictionary = {**common_units_dictionary,**{param:u.dimensionless_unscaled for param in other_params}}
common_units_dictionary = {**common_units_dictionary,**{param:u.dimensionless_unscaled for param in nongr_params}}
common_units_dictionary = {**common_units_dictionary,**{param:u.dimensionless_unscaled for param in tidal_params}}

#SI Units Dictionary
SI_units_dictionary = {**{mass:u.kg for mass in mass_params_},**{mass_dimensionless:u.dimensionless_unscaled for mass_dimensionless in mass_params if mass_dimensionless not in mass_params_}}
SI_units_dictionary = {**SI_units_dictionary,**{"distance":u.m},**common_units_dictionary}

#Cosmo Units Dictionary
Cosmo_units_dictionary = {**{mass:u.solMass for mass in mass_params_},**{mass_dimensionless:u.dimensionless_unscaled for mass_dimensionless in mass_params if mass_dimensionless not in mass_params_}}
Cosmo_units_dictionary = {**Cosmo_units_dictionary,**{"distance":u.pc},**common_units_dictionary}

units_dict  = {"S.I.":SI_units_dictionary,"Cosmo":Cosmo_units_dictionary}

#Non-Dimensional quantities

Non_d_list = np.concatenate([_other_params,arr_params])

default_dict = {"mass1" : 1.*u.solMass,
              "mass2" : 1.*u.solMass,
              "spin1x" : 0.*u.dimensionless_unscaled,
              "spin1y" : 0.*u.dimensionless_unscaled,
              "spin1z" : 0.*u.dimensionless_unscaled,
              "spin2x" : 0.*u.dimensionless_unscaled,
              "spin2y" : 0.*u.dimensionless_unscaled,
              "spin2z" : 0.*u.dimensionless_unscaled,
              "deltaT" : 1./64.*u.s,
              "f22_start" : 20.*u.Hz,
              "f22_ref" : 20.*u.Hz,
              "phi_ref": 0.*u.rad,
              "distance" : 1.*u.Mpc,
              "inclination" : 0.*u.rad,
              "eccentricity" : 0.*u.dimensionless_unscaled,
              "longAscNodes" : 0.*u.rad,
              "meanPerAno" : 0.*u.rad
              }
